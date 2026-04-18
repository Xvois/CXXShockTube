//
// main.cpp
//
// Fluid dynamics solver for PH30110 Computational Astrophysics
// Solves Euler equations for:
//   Problem A: Cartesian shock tube (standard test case)
//   Problem B: Spherical shock tube (supernova remnant)
//
// Uses Lax-Friedrichs scheme with operator splitting for spherical geometry.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "include/constants.h"
#include "include/physics.h"
#include "include/numerical.h"
#include "include/analytical.h"
#include "include/problems.h"

// =========================================================================
// Output Functions
// =========================================================================

/**
 * Writes simulation data to CSV file.
 * Creates new file if currentTime == 0, appends otherwise.
 */
void outputCSV(const std::vector<Conserved>& grid, double currentTime, const std::string& filename) {
    std::ofstream file;
    if (currentTime == 0.0) {
        file.open(filename, std::ios::out); // Overwrite/Create new
        file << "time,x,density,velocity,pressure,internal_energy\n";
    } else {
        file.open(filename, std::ios::app); // Append
    }

    for (int i = 0; i < N_ZONES; ++i) {
        double x = X_MIN + (i + 0.5) * DX;
        Primitive w = consToPrim(grid[i]);
        // Compute specific internal energy
        double e_internal = w.p / (w.rho * (GAMMA - 1.0));
        
        file << std::fixed << std::setprecision(6)
             << currentTime << "," << x << "," << w.rho << "," 
             << w.v << "," << w.p << "," << e_internal << "\n";
    }
    file.close();
}

// =========================================================================
// Problem A: Cartesian Shock Tube
// =========================================================================

void solveProblemA(const std::string& outputFilename) {
    std::cout << "=== Problem A: Cartesian Shock Tube ===" << std::endl;
    std::cout << "Left state (x < 0.3): ρ=1, p=1, v=0.75" << std::endl;
    std::cout << "Right state (x > 0.3): ρ=0.125, p=0.1, v=0" << std::endl;
    std::cout << "Target time: t = 0.2" << std::endl;

    // Initialise grid with Problem A initial conditions
    std::vector<Conserved> grid = initialiseProblemA();

    double t = 0.0;
    const double T_FINAL = 0.2;
    const double OUTPUT_INTERVAL = 0.05;
    double next_output_t = 0.0;
    int step = 0;

    while (t <= T_FINAL) {
        // Output snapshot
        if (t >= next_output_t - 1e-10) {
            outputCSV(grid, t, outputFilename);
            next_output_t += OUTPUT_INTERVAL;
        }

        if (t >= T_FINAL) break;

        // Calculate CFL-limited timestep
        double dt = calculateTimeStep(grid);

        // Clip dt to hit output intervals exactly
        if (t + dt > next_output_t && next_output_t <= T_FINAL) dt = next_output_t - t;
        if (t + dt > T_FINAL) dt = T_FINAL - t;

        // Advance one timestep using Lax-Friedrichs
        grid = updateLaxFriedrichs(grid, dt);
        t += dt;
        step++;
    }

    std::cout << "  Completed: " << step << " timesteps" << std::endl;
    std::cout << "  Data saved to " << outputFilename << std::endl;

    // Compute and save exact solution for comparison
    std::cout << "  Computing exact solution..." << std::endl;
    std::vector<Conserved> exact_grid(N_ZONES);
    Primitive L = {1.0, 0.75, 1.0};  // Left state
    Primitive R = {0.125, 0.0, 0.1};  // Right state
    for (int i = 0; i < N_ZONES; ++i) {
        double x = X_MIN + (i + 0.5) * DX;
        Primitive w_exact = sampleExactSolution(L, R, 0.3, x, T_FINAL);
        exact_grid[i] = primToCons(w_exact);
    }
    {
        std::ofstream file("exact_solution.csv", std::ios::out);
        file << "time,x,density,velocity,pressure,internal_energy" << std::endl;
        for (int i = 0; i < N_ZONES; ++i) {
            double x = X_MIN + (i + 0.5) * DX;
            Primitive w = consToPrim(exact_grid[i]);
            double e_internal = w.p / (w.rho * (GAMMA - 1.0));
            file << std::fixed << std::setprecision(6)
                 << T_FINAL << "," << x << "," << w.rho << ","
                 << w.v << "," << w.p << "," << e_internal << std::endl;
        }
        file.close();
    }
    std::cout << "  Exact solution saved to exact_solution.csv" << std::endl;
}

// =========================================================================
// Problem B: Spherical Shock Tube
// =========================================================================

void solveProblemB(const std::string& outputFilename) {
    std::cout << "=== Problem B: Spherical Shock Tube ===" << std::endl;
    std::cout << "Left state (r < 0.4): ρ=1, p=1, v=0" << std::endl;
    std::cout << "Right state (r > 0.4): ρ=0.125, p=0.1, v=0" << std::endl;
    std::cout << "Target time: t = 0.25 (spherical)" << std::endl;

    // Initialise grid with Problem B initial conditions
    std::vector<Conserved> grid = initialiseProblemB();

    double t = 0.0;
    const double T_FINAL = 0.25;
    const double OUTPUT_INTERVAL = 0.05;
    double next_output_t = 0.0;
    int step = 0;

    while (t <= T_FINAL) {
        // Output snapshot
        if (t >= next_output_t - 1e-10) {
            outputCSV(grid, t, outputFilename);
            next_output_t += OUTPUT_INTERVAL;
        }

        if (t >= T_FINAL) break;

        // Calculate CFL-limited timestep (uses sound speed for sphericity)
        double dt = calculateTimeStep(grid);

        // Clip dt to hit output intervals exactly
        if (t + dt > next_output_t && next_output_t <= T_FINAL) dt = next_output_t - t;
        if (t + dt > T_FINAL) dt = T_FINAL - t;

        // Advance one timestep using spherical Lax-Friedrichs with source terms
        // The spherical solver applies the 2p/r geometric source term using operator splitting
        grid = updateSphericalLaxFriedrichs(grid, dt);
        t += dt;
        step++;
    }

    std::cout << "  Completed: " << step << " timesteps" << std::endl;
    std::cout << "  Data saved to " << outputFilename << std::endl;
}

// =========================================================================
// Main Entry Point
// =========================================================================

int main() {
    std::cout << "Fluid Solver - PH30110 Computational Astrophysics" << std::endl;
    std::cout << "Grid: N = " << N_ZONES << ", Gamma = " << GAMMA << std::endl;
    std::cout << "Domain: [0, 1], dx = " << DX << std::endl;
    std::cout << std::endl;

    // Solve Problem A (Cartesian)
    solveProblemA("12345_problemA_results.csv");
    
    // Solve Problem B (Spherical)
    solveProblemB("12345_problemB_results.csv");
    
    std::cout << std::endl << "All simulations complete." << std::endl;
    return 0;
}
