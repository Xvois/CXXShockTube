//
// solver.h
//
// Defines the simulation driver, output handling, and problem-specific solvers.
// Orchestrates time-stepping, CSV output, and exact-solution comparison.
//
#ifndef FLUIDSOLVER_SOLVER_H
#define FLUIDSOLVER_SOLVER_H

#include <string>
#include <vector>
#include "physics.h"

/**
 * Writes simulation data to a CSV file.
 * Creates a new file with headers if currentTime == 0;
 * appends to an existing file otherwise.
 *
 * @param grid        Current grid state (conserved variables)
 * @param currentTime Simulation time for this snapshot
 * @param filename    Output file path
 */
void outputCSV(const std::vector<Conserved>& grid, double currentTime, const std::string& filename);

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
 * @param outputFilename Path for the CSV output file
 */
void solveProblemA(const std::string& outputFilename);

/**
 * Solves Problem B: Spherical shock tube.
 *
 * Left state (r < 0.4): rho = 1.0, p = 1.0, v = 0.0
 * Right state (r > 0.4): rho = 0.125, p = 0.1, v = 0.0
 * Snapshot target: t = 0.25
 *
 * Spherical Euler equations with geometric source terms:
 *   dU/dt + dF/dr = S,  where S = (0, 2*p/r, 0)
 *
 * Boundary conditions:
 *   - Inner boundary (r = 0): reflecting wall (v = 0, physical symmetry)
 *   - Outer boundary (r = 1.0): outflow (non-reflecting)
 *
 * @param outputFilename Path for the CSV output file
 */
void solveProblemB(const std::string& outputFilename);

#endif // FLUIDSOLVER_SOLVER_H
