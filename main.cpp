//
// main.cpp
//
// Fluid dynamics solver for PH30110 Computational Astrophysics.
// Solves Euler equations for:
//   Problem A: Cartesian shock tube (standard test case)
//   Problem B: Spherical shock tube (supernova remnant)
//
// Uses Lax-Friedrichs scheme with operator splitting for spherical geometry.
// Boundary conditions: outflow (non-reflecting) for Problem A,
// mixed reflecting/outflow for Problem B.
//
// Outputs CSV data to the directory containing the executable,
// ensuring files are found regardless of CWD or build directory name.
//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <future>
#include <unistd.h>
#include "constants.h"
#include "solver.h"
#include "analytical.h"
#include "problems.h"

namespace fs = std::filesystem;

/**
 * Resolve the directory where the executable lives.
 * Works on Linux/macOS (via argv[0]) and Windows.
 * Falls back to "." (CWD) if the path cannot be resolved.
 */
std::string resolve_exec_dir(int argc, char* argv[]) {
    std::string exe_path;
    if (argc > 0 && argv[0] != nullptr) {
        exe_path = argv[0];
    }
    // If argv[0] is just a name (no /), the executable is in PATH or CWD
    if (exe_path.find('/') == std::string::npos && exe_path.find('\\') == std::string::npos) {
        // Try /proc/self/exe on Linux
#ifdef __linux__
        char buf[4096];
        ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
        if (len > 0) {
            exe_path = std::string(buf, len);
        }
#endif
    }
    fs::path p(exe_path);
    fs::path dir = p.parent_path();
    if (dir.empty() || dir == ".") {
        return ".";
    }
    return dir.string();
}

int main(int argc, char* argv[]) {
    // Resolve output directory relative to executable location
    std::string exec_dir = resolve_exec_dir(argc, argv);

    std::cout << "Fluid Solver - PH30110 Computational Astrophysics" << std::endl;
    std::cout << "Grid: N = " << N_ZONES << ", Gamma = " << GAMMA << std::endl;
    std::cout << "Domain: [0, 1], dx = " << DX << std::endl;
    std::cout << "CFL Number = " << CFL_NUMBER << std::endl;
    std::cout << std::endl;

    std::string outputA = exec_dir + "/12345_problemA_results.csv";
    std::string outputB = exec_dir + "/12345_problemB_results.csv";
    std::string consA   = exec_dir + "/conservation_12345_problemA.csv";
    std::string consB   = exec_dir + "/conservation_12345_problemB.csv";
    std::string exactA  = exec_dir + "/exact_solution.csv";

    // Solve Problems A and B in parallel using std::async.
    // Each problem is independent, so they can run concurrently on separate threads.
    std::future<void> futA = std::async(std::launch::async, solveProblemA, outputA, exactA);
    std::future<void> futB = std::async(std::launch::async, solveProblemB, outputB);

    // Wait for both problems to complete
    futA.get();
    std::cout << std::endl;
    futB.get();

    std::cout << std::endl << "All simulations complete." << std::endl;
    std::cout << "Output directory: " << exec_dir << std::endl;
    return 0;
}
