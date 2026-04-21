/**
 * @file problems.h
 * @brief Initial condition generators and problem-specific solvers.
 *
 * Declares two standard shock tube test cases for the Euler equations:
 *
 *   Problem A (Cartesian): Standard 1D shock tube with discontinuity at x = 0.3,
 *     outflow boundaries at both ends. Represents the classic Woodward-Colella
 *     "blast wave" setup used to verify numerical schemes.
 *
 *   Problem B (Spherical): Spherical shock tube with discontinuity at r = 0.4,
 *     mixed boundary conditions (reflecting at r = 0, outflow at r = 1.0).
 *     Includes geometric source terms arising from spherical symmetry.
 *
 * Spherical geometry adds a geometric source term to the Euler equations:
 *   dU/dt + dF/dr = S,   where S = (2*rho*v/r, 2*p/r, 0)
 *
 * This source term accounts for the expanding cross-sectional area of spherical
 * shells. It is handled via Strang operator splitting (half-step, full-step,
 * half-step) for second-order accuracy.
 *
 * Boundary conditions:
 *   - Outflow (non-reflecting): ghost cells copy interior state, allowing
 *     waves to exit the domain without reflection.
 *   - Reflecting (v = 0): inner boundary enforces zero velocity (symmetry).
 *
 * @note Problem A is the baseline test case. Problem B extends the physics
 *       to spherical geometry and tests the source-term treatment.
 */
#ifndef FLUIDSOLVER_PROBLEMS_H
#define FLUIDSOLVER_PROBLEMS_H

#include <vector>
#include "physics.h"
#include "constants.h"

namespace fluidsolver {

/**
 * @brief Initialise Problem A: Cartesian shock tube.
 *
 * Sets up a Riemann problem with discontinuity at x = 0.3:
 *   Left state  (x < 0.3): rho = 1.0,  v = 0.75,  p = 1.0
 *   Right state (x > 0.3): rho = 0.125, v = 0.0,   p = 0.1
 *
 * Boundary conditions: outflow (non-reflecting) at both ends.
 * Ghost cells copy the interior state at x=0 and x=1.
 *
 * @return Grid state as a vector of conserved variables, one per zone.
 */
std::vector<Conserved> initialiseProblemA();

/**
 * @brief Initialise Problem B: Spherical shock tube.
 *
 * Sets up a spherical Riemann problem with discontinuity at r = 0.4:
 *   Left state  (r < 0.4): rho = 1.0,  v = 0.0,   p = 1.0
 *   Right state (r > 0.4): rho = 0.125, v = 0.0,   p = 0.1
 *
 * Boundary conditions:
 *   - Inner boundary (r = 0): reflecting wall, v = 0 (physical symmetry).
 *     Source terms are skipped at r = 0 to avoid the 0/0 singularity.
 *   - Outer boundary (r = 1.0): outflow (non-reflecting), ghost cells copy interior.
 *
 * @return Grid state as a vector of conserved variables, one per zone.
 */
std::vector<Conserved> initialiseProblemB();

/**
 * @brief Solve Problem B: spherical Euler equations with operator splitting.
 *
 * Implements Strang splitting for the geometric source term:
 *   1. Half-step of source term:   S * dt/2
 *   2. Full-step Lax-Friedrichs flux update (no source): U* = U^n - dt*dF/dr
 *   3. Remaining half-step of source term:  S * dt/2
 *   4. Apply boundary conditions after each substep
 *
 * The operator splitting is second-order accurate in time.
 * Spatial discretization uses the conservative Lax-Friedrichs scheme,
 * first-order accurate but robust (monotone, no spurious oscillations).
 *
 * @param grid Current grid state (conserved variables).
 * @param dt Timestep size (must satisfy CFL condition).
 * @return Updated grid state after one timestep.
 */
std::vector<Conserved> updateSphericalLaxFriedrichs(
    const std::vector<Conserved>& grid, double dt);

}  // namespace fluidsolver

#endif // FLUIDSOLVER_PROBLEMS_H
