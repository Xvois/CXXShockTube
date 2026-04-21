/**
 * @file solver.h
 * @brief Simulation driver, CSV output handling, and problem-specific solvers.
 *
 * Orchestrates time-stepping loops, writes simulation snapshots and conservation
 * diagnostics to CSV files, and computes exact Riemann solutions for comparison.
 *
 * Output files:
 *   - Problem A results: {filename}_problemA_results.csv
 *   - Problem B results: {filename}_problemB_results.csv
 *   - Conservation A:    conservation_problemA.csv
 *   - Conservation B:    conservation_problemB.csv
 *   - Exact solution A:  exact_solution_problemA.csv
 *
 * All CSV files are written with 6-digit precision.
 *
 * @note Output files are placed in the same directory as the output filename,
 *       not the working directory. This ensures files are found regardless of
 *       CWD or build directory name.
 */
#ifndef FLUIDSOLVER_SOLVER_H
#define FLUIDSOLVER_SOLVER_H

namespace fluidsolver {}

#include <string>
#include <vector>
#include "physics.h"

/**
 * @brief Writes simulation data to a CSV file.
 *
 * Creates a new file with headers if currentTime == 0;
 * appends to an existing file otherwise (for time-series output).
 *
 * Output columns: time, x, density, velocity, pressure, internal_energy
 *
 * Specific internal energy is computed as:
 *   e = p / (rho * (gamma - 1))
 *
 * This is the thermal energy per unit mass, excluding kinetic energy.
 *
 * @param grid        Current grid state (conserved variables).
 * @param currentTime Simulation time for this snapshot.
 * @param filename    Output file path.
 */
void outputCSV(const std::vector<Conserved>& grid, double currentTime, const std::string& filename);

/**
 * @brief Writes conservation diagnostics to CSV.
 *
 * Computes the volume-integrated mass, momentum, and total energy at each
 * output timestep to verify numerical conservation properties.
 *
 * Integrated quantities:
 *   M = sum(rho_i * dx)  — total mass
 *   P = sum(rho*v_i * dx)  — total momentum
 *   E = sum(E_i * dx)  — total energy
 *
 * For the Euler equations with outflow boundaries, these quantities
 * should be constant for waves remaining within the domain, and
 * should decrease only when waves carry mass/momentum/energy out.
 *
 * @param grid        Current grid state (conserved variables).
 * @param currentTime Simulation time for this snapshot.
 * @param filename    Output file path (created with headers if nonexistent).
 */
void outputConservation(const std::vector<Conserved>& grid, double currentTime, const std::string& filename);

/**
 * @brief Solves Problem A: Cartesian shock tube.
 *
 * Runs the full simulation loop:
 *   1. Initialise Problem A initial conditions
 *   2. Time-step using Lax-Friedrichs scheme with CFL-limited dt
 *   3. Output snapshots at fixed intervals (every N_OUTPUT steps)
 *   4. Write conservation diagnostics
 *   5. Write exact Riemann solution for comparison
 *
 * @param outputFilename    Path for the simulation CSV output file.
 * @param exactOutputFile   Path for the exact solution CSV output file.
 */
void solveProblemA(const std::string& outputFilename, const std::string& exactOutputFile);

/**
 * @brief Solves Problem B: spherical shock tube.
 *
 * Runs the full simulation loop with operator splitting for the
 * geometric source term:
 *   1. Initialise Problem B initial conditions
 *   2. Time-step using spherical Lax-Friedrichs with CFL-limited dt
 *   3. Output snapshots at fixed intervals
 *   4. Write conservation diagnostics
 *
 * @param outputFilename    Path for the simulation CSV output file.
 */
void solveProblemB(const std::string& outputFilename);

} // namespace fluidsolver

#endif // FLUIDSOLVER_SOLVER_H
