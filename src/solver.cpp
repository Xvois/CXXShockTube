//
// solver.cpp
//
// Defines the simulation driver, output handling, and problem-specific solvers.
// Orchestrates time-stepping, CSV output, and exact-solution comparison.
//
// All output files are written to the directory of the output filename,
// ensuring CSVs are found regardless of CWD or build directory.
//

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include "solver.h"
#include "physics.h"
#include "numerical.h"
#include "problems.h"
#include "constants.h"
#include "analytical.h"

namespace fs = std::filesystem;

// =========================================================================
// Helper: extract directory from a path
// =========================================================================
static std::string dir_of(const std::string& path) {
    fs::path p(path);
    fs::path d = p.parent_path();
    if (d.empty()) return ".";
    return d.string();
}

// =========================================================================
// Output
// =========================================================================

/**
 * Writes simulation data to CSV.
 * Creates new file with headers if currentTime == 0; appends otherwise.
 * Computes specific internal energy: e = p / (rho * (gamma - 1))
 */
void outputCSV(const std::vector<Conserved>& grid, double currentTime, const std::string& filename) {
    std::ofstream file;
    if (currentTime == 0.0) {
        file.open(filename, std::ios::out);
        file << "time,x,density,velocity,pressure,internal_energy\n";
    } else {
        file.open(filename, std::ios::app);
    }

    for (int i = 0; i < N_ZONES; ++i) {
        double x = X_MIN + (i + 0.5) * DX;
        Primitive w = consToPrim(grid[i]);
        double e_internal = w.p / (w.rho * (GAMMA - 1.0));

        file << std::fixed << std::setprecision(6)
             << currentTime << "," << x << "," << w.rho << ","
             << w.v << "," << w.p << "," << e_internal << "\n";
    }
    file.close();
}

// =========================================================================
// Conservation Diagnostics
// =========================================================================

/**
 * Writes conservation diagnostics to CSV.
 * Computes integrated mass, momentum, and energy at each output timestep.
 * Creates new file with headers if currentTime == 0; appends otherwise.
 *
 * Used for verifying numerical conservation properties:
 *   - Problem A (Cartesian): near-exact conservation (ghost cell outflow)
 *   - Problem B (Spherical): some mass drift due to approximate source terms
 */
void outputConservation(const std::vector<Conserved>& grid, double currentTime,
                        const std::string& filename) {
    std::ofstream file;
    if (currentTime == 0.0) {
        file.open(filename, std::ios::out);
        file << "time,mass,momentum,energy\n";
    } else {
        file.open(filename, std::ios::app);
    }

    double mass = 0.0, momentum = 0.0, energy = 0.0;
    for (int i = 0; i < N_ZONES; ++i) {
        mass   += grid[i].mass;
        momentum += grid[i].mom;
        energy += grid[i].energy;
    }

    file << std::fixed << std::setprecision(6)
         << currentTime << "," << mass << "," << momentum << "," << energy << "\n";
    file.close();
}

// =========================================================================
// Problem A: Cartesian Shock Tube
// =========================================================================

/**
 * Solves Problem A: Cartesian shock tube.
 *
 * Left state (x < 0.3): rho = 1.0, p = 1.0, v = 0.75
 * Right state (x > 0.3): rho = 0.125, p = 0.1, v = 0.0
 * Snapshot target: t = 0.2
 *
 * Boundary conditions: outflow (non-reflecting) at both ends.
 * Ghost cells copy the interior state so waves exit freely.
 *
 * @param outputFilename    Path for the simulation CSV output file
 * @param exactOutputFile   Path for the exact solution CSV output file
 */
void solveProblemA(const std::string& outputFilename, const std::string& exactOutputFile) {
    std::cout << "=== Problem A: Cartesian Shock Tube ===" << std::endl;
    std::cout << "Left state (x < 0.3): rho=1, p=1, v=0.75" << std::endl;
    std::cout << "Right state (x > 0.3): rho=0.125, p=0.1, v=0" << std::endl;
    std::cout << "Target time: t = 0.2" << std::endl;

    std::vector<Conserved> grid = initialiseProblemA();

    double t = 0.0;
    const double T_FINAL = 0.2;
    const double OUTPUT_INTERVAL = 0.05;
    double next_output_t = 0.0;
    int step = 0;
    int output_count = 0;

    // Open conservation CSV once (headers at t=0)
    std::string consDir = dir_of(outputFilename);
    std::string consPath = consDir + "/conservation_12345_problemA.csv";
    {
        std::ofstream conservationFile(consPath, std::ios::out);
        conservationFile << "time,mass,momentum,energy\n";
        conservationFile.close();
    }

    while (t < T_FINAL) {
        if (t >= next_output_t - 1e-10) {
            outputCSV(grid, t, outputFilename);
            outputConservation(grid, t, consPath);
            next_output_t += OUTPUT_INTERVAL;
            output_count++;
        }

        if (t >= T_FINAL) break;

        double dt = calculateTimeStep(grid);

        // Clip dt to hit output times and final time exactly
        if (t + dt > next_output_t && next_output_t <= T_FINAL) {
            dt = next_output_t - t;
        }
        if (t + dt > T_FINAL) {
            dt = T_FINAL - t;
        }

        grid = updateLaxFriedrichs(grid, dt);
        t += dt;
        step++;
    }

    // Final output at T_FINAL
    outputCSV(grid, t, outputFilename);
    outputConservation(grid, t, consPath);
    output_count++;

    std::cout << "  Completed: " << step << " timesteps" << std::endl;
    std::cout << "  Outputs written: " << output_count << std::endl;
    std::cout << "  Final t value: " << t << std::endl;
    std::cout << "  Data saved to " << outputFilename << std::endl;

    // Compute and save exact solution for comparison at t = 0.2
    std::cout << "  Computing exact solution..." << std::endl;
    std::vector<Conserved> exact_grid(N_ZONES);
    Primitive L = {1.0, 0.75, 1.0};   // Left state
    Primitive R = {0.125, 0.0, 0.1};  // Right state
    for (int i = 0; i < N_ZONES; ++i) {
        double x = X_MIN + (i + 0.5) * DX;
        Primitive w_exact = sampleExactSolution(L, R, 0.3, x, T_FINAL);
        exact_grid[i] = primToCons(w_exact);
    }
    {
        std::ofstream file(exactOutputFile, std::ios::out);
        file << "time,x,density,velocity,pressure,internal_energy\n";
        for (int i = 0; i < N_ZONES; ++i) {
            double x = X_MIN + (i + 0.5) * DX;
            Primitive w = consToPrim(exact_grid[i]);
            double e_internal = w.p / (w.rho * (GAMMA - 1.0));
            file << std::fixed << std::setprecision(6)
                 << T_FINAL << "," << x << "," << w.rho << ","
                 << w.v << "," << w.p << "," << e_internal << "\n";
        }
        file.close();
    }
    std::cout << "  Exact solution saved to " << exactOutputFile << std::endl;
}

// =========================================================================
// Problem B: Spherical Shock Tube
// =========================================================================

/**
 * Solves Problem B: Spherical shock tube.
 *
 * Left state (r < 0.4): rho = 1.0, p = 1.0, v = 0.0
 * Right state (r > 0.4): rho = 0.125, p = 0.1, v = 0.0
 * Target time: t = 0.25
 *
 * Boundary conditions:
 *   - Inner boundary (r = 0): reflecting wall (v = 0)
 *   - Outer boundary (r = 1.0): outflow (non-reflecting)
 */
void solveProblemB(const std::string& outputFilename) {
    std::cout << "=== Problem B: Spherical Shock Tube ===" << std::endl;
    std::cout << "Left state (r < 0.4): rho=1, p=1, v=0" << std::endl;
    std::cout << "Right state (r > 0.4): rho=0.125, p=0.1, v=0" << std::endl;
    std::cout << "Target time: t = 0.25 (spherical)" << std::endl;

    std::vector<Conserved> grid = initialiseProblemB();

    double t = 0.0;
    const double T_FINAL = 0.25;
    const double OUTPUT_INTERVAL = 0.05;
    double next_output_t = 0.0;
    int step = 0;
    int output_count = 0;

    // Open conservation CSV once (headers at t=0)
    std::string consDir = dir_of(outputFilename);
    std::string consPath = consDir + "/conservation_12345_problemB.csv";
    {
        std::ofstream conservationFile(consPath, std::ios::out);
        conservationFile << "time,mass,momentum,energy\n";
        conservationFile.close();
    }

    while (t < T_FINAL) {
        if (t >= next_output_t - 1e-10) {
            outputCSV(grid, t, outputFilename);
            outputConservation(grid, t, consPath);
            next_output_t += OUTPUT_INTERVAL;
            output_count++;
        }

        if (t >= T_FINAL) break;

        double dt = calculateTimeStep(grid);

        if (t + dt > next_output_t && next_output_t <= T_FINAL) {
            dt = next_output_t - t;
        }
        if (t + dt > T_FINAL) {
            dt = T_FINAL - t;
        }

        grid = updateSphericalLaxFriedrichs(grid, dt);
        t += dt;
        step++;
    }

    outputCSV(grid, t, outputFilename);
    outputConservation(grid, t, consPath);
    output_count++;

    std::cout << "  Completed: " << step << " timesteps" << std::endl;
    std::cout << "  Outputs written: " << output_count << std::endl;
    std::cout << "  Final t value: " << t << std::endl;
    std::cout << "  Data saved to " << outputFilename << std::endl;
}
