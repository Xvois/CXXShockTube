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
 * Writes conservation diagnostics to CSV.
 * Computes integrated mass, momentum, and energy at each output timestep.
 * Used for verifying numerical conservation properties.
 */
void outputConservation(const std::vector<Conserved>& grid, double currentTime, const std::string& filename);

/**
 * Solves Problem A: Cartesian shock tube.
 *
 * @param outputFilename    Path for the simulation CSV output file
 * @param exactOutputFile   Path for the exact solution CSV output file
 */
void solveProblemA(const std::string& outputFilename, const std::string& exactOutputFile);

/**
 * Solves Problem B: Spherical shock tube.
 *
 * @param outputFilename    Path for the simulation CSV output file
 */
void solveProblemB(const std::string& outputFilename);

#endif // FLUIDSOLVER_SOLVER_H
