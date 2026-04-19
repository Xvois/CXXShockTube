//
// numerical.h
//
// Declares the CFL-limited timestep calculator and the Lax-Friedrichs
// time-stepping update for both Cartesian and spherical geometries.
//

#ifndef FLUIDSOLVER_NUMERICAL_H
#define FLUIDSOLVER_NUMERICAL_H

#include <vector>
#include "physics.h"
#include "constants.h"

/**
 * Computes the maximum allowable timestep from the CFL condition.
 * Returns the maximum stable dt. If max_signal_speed is zero,
 * returns a safe default (1.0).
 */
double calculateTimeStep(const std::vector<Conserved>& grid);

/**
 * Advances the grid by one timestep using the conservative Lax-Friedrichs
 * scheme with outflow (non-reflecting) boundaries.
 */
std::vector<Conserved> updateLaxFriedrichs(const std::vector<Conserved>& grid, double dt);

#endif // FLUIDSOLVER_NUMERICAL_H
