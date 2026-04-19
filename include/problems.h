//
// problems.h
//
// Declares initial condition generators and the spherical solver.
//
// Problem A: Cartesian shock tube with outflow boundaries.
// Problem B: Spherical shock tube with geometric source terms
//            and mixed boundary conditions.
//

#ifndef FLUIDSOLVER_PROBLEMS_H
#define FLUIDSOLVER_PROBLEMS_H

#include "physics.h"
#include "constants.h"
#include <vector>

// --- Problem A: Cartesian Shock Tube ---

/**
 * Initialise Problem A: Cartesian shock tube.
 * Discontinuity at x = 0.3 with outflow boundaries at both ends.
 */
std::vector<Conserved> initialiseProblemA();

// --- Problem B: Spherical Shock Tube ---

/**
 * Initialise Problem B: Spherical shock tube.
 * Discontinuity at r = 0.4 with mixed boundaries:
 *   reflecting (r = 0) and outflow (r = 1.0).
 */
std::vector<Conserved> initialiseProblemB();

/**
 * Solve the spherical Euler equations using operator-splitting:
 *   1. Half-step of geometric source term S = (0, 2p/r, 0)
 *   2. Full-step Lax-Friedrichs flux update
 *   3. Remaining half-step of source term
 *   4. Boundary conditions (v=0 at inner, outflow at outer)
 *
 * Second-order in time, first-order in space.
 */
std::vector<Conserved> updateSphericalLaxFriedrichs(
    const std::vector<Conserved>& grid, double dt);

#endif // FLUIDSOLVER_PROBLEMS_H
