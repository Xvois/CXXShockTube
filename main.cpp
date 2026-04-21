/**
 * @file main.cpp
 * @brief Fluid dynamics solver for PH30110 Computational Astrophysics.
 *
 * Solves the 1D Euler equations for two shock tube test cases:
 *   Problem A: Cartesian shock tube (standard test case, Woodward-Colella)
 *   Problem B: Spherical shock tube (supernova remnant model)
 *
 * Numerical methods:
 *   - Lax-Friedrichs scheme with Rusanov flux for spatial discretization
 *   - Strang operator splitting for spherical geometric source terms
 *   - CFL-limited timestep for stability
 *
 * Boundary conditions:
 *   - Problem A: outflow (non-reflecting) at both ends
 *   - Problem B: mixed — reflecting wall at r=0, outflow at r=1.0
 *
 * Output: CSV files containing simulation snapshots and conservation
 * diagnostics. All output files are placed in the directory containing
 * the executable, ensuring files are found regardless of CWD or
 * build directory name.
 *
 * Parallel execution: Problems A and B are solved concurrently using
 * std::async (one std::thread per problem).
 *
 * Compile and run:
 *   g++ -std=c++20 -O2 -fopenmp -o fluid_solver main.cpp src/*.cpp
 *   ./fluid_solver
 *
 * @note This is a fixed-grid solver with N_ZONES = 100 zones.
 *       Grid refinement requires recompilation.
 */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <future>
#include "constants.h"
#include "solver.h"
#include "analytical.h"
#include "problems.h"

using fluidsolver::solveProblemA;
using fluidsolver::solveProblemB;
using fluidsolver::N_ZONES;
using fluidsolver::GAMMA;
using fluidsolver::DX;
using fluidsolver::CFL_NUMBER;

namespace fs = std::filesystem;

/**
 * @brief Resolve the directory where the executable lives.
 *
 * Falls back to "." (CWD) if the path cannot be resolved.
 *
 * @param argc Argument count (unused except for checking argc > 0).
 * @param argv Argument vector (argv[0] contains the executable path).
 * @return Directory path string. Returns "." if the executable directory
 *         cannot be determined.
 */
std::string resolve_exec_dir(int argc, char* argv[]) {
    if (argc == 0 || argv[0] == nullptr) {
        return ".";
    }

    std::string exe_path = argv[0];
    fs::path p(exe_path);
    fs::path dir = p.parent_path();

    if (dir.empty() || dir == ".") {
        return ".";
    }
    return dir.string();
}

/**
 * @brief Main entry point: run Problem A and Problem B in parallel.
 *
 * Steps:
 *   1. Resolve output directory relative to executable location.
 *   2. Set output file paths for both problems.
 *   3. Launch Problem A solver on std::async (creates a new thread).
 *   4. Launch Problem B solver on std::async (creates a new thread).
 *   5. Wait for both to complete and print summary.
 *
 * @param argc Argument count (expected to be >= 1).
 * @param argv Argument vector (argv[0] contains executable path).
 * @return 0 on success.
 */
int main(int argc, char* argv[]) {
    // Resolve output directory relative to executable location
    std::string exec_dir = resolve_exec_dir(argc, argv);

    std::cout << "Fluid Solver - PH30110 Computational Astrophysics" << std::endl;
    std::cout << "Grid: N = " << N_ZONES << ", Gamma = " << GAMMA << std::endl;
    std::cout << "Domain: [0, 1], dx = " << DX << std::endl;
    std::cout << "CFL Number = " << CFL_NUMBER << std::endl;
    std::cout << std::endl;

    // Output filenames (simplified)
    std::string outputA = exec_dir + "/shock_A.csv";
    std::string outputB = exec_dir + "/shock_B.csv";
    std::string exactA  = exec_dir + "/exact_A.csv";
    std::string consA   = exec_dir + "/conservation_A.csv";
    std::string consB   = exec_dir + "/conservation_B.csv";

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
