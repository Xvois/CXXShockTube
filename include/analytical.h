//
// analytical.h
//
// Declares analytical solution helpers for the shock tube problems:
// sound speed, pressure function, star-pressure iteration, and
// exact solution sampling.
//

#ifndef FLUIDSOLVER_ANALYTICAL_H
#define FLUIDSOLVER_ANALYTICAL_H

#include "constants.h"
#include "physics.h"

/** Compute local sound speed: c = sqrt(gamma * p / rho) */
double calculateSoundSpeed(const Primitive& w);

/**
 * Evaluate the pressure-difference function for left/right waves.
 * Positive result (p_star > p_side) → shock wave;
 * Negative result (p_star <= p_side) → rarefaction wave.
 */
double pressureFunction(double p_star, double p_side, double rho_side, double cs_side);

/**
 * Solve for the pressure in the star region (Region 3) using
 * Newton-Raphson iteration on the two-sided Riemann problem.
 */
double findStarPressure(const Primitive& L, const Primitive& R);

/**
 * Sample the exact solution at a given position and time.
 * Returns the primitive variables at the specified point.
 */
Primitive sampleExactSolution(const Primitive& L, const Primitive& R,
                              double x_initial, double x_current, double t);

#endif // FLUIDSOLVER_ANALYTICAL_H
