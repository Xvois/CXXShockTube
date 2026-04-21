/**
 * @file numerical.h
 * @brief CFL-limited timestep calculation and Lax-Friedrichs time-stepping
 *        for both Cartesian and spherical geometries.
 *
 * Numerical methods:
 *
 *   CFL-limited timestep:
 *     dt <= CFL_NUMBER * dx / max(|v| + c)
 *   where c = sqrt(gamma * p / rho) is the local sound speed.
 *   This ensures the numerical domain of dependence contains the physical
 *   domain of dependence for stability.
 *
 *   Conservative Lax-Friedrichs scheme:
 *     U_i^{n+1} = 0.5*(U_{i-1}^n + U_{i+1}^n) - dt/(2*dx) * (F_{i+1/2} - F_{i-1/2})
 *
 *   This scheme has a built-in numerical diffusion term:
 *     U_i^{n+1} = U_i^n - dt/(2*dx)*(F_{i+1/2} - F_{i-1/2)) + (dt^2/(2*dx^2))*(F_{i+1} - 2*F_i + F_{i-1})
 *   The artificial viscosity (second term) damps oscillations near shocks.
 *
 *   Boundary conditions: outflow (non-reflecting) via ghost cell copy.
 *
 * @note OpenMP parallelisation is applied to the inner loop (zone-by-zone
 *       update) when available. This provides speedup on multi-core systems.
 */
#ifndef FLUIDSOLVER_NUMERICAL_H
#define FLUIDSOLVER_NUMERICAL_H

#include <vector>
#include "physics.h"
#include "constants.h"

/**
 * @brief Computes the maximum allowable timestep from the CFL condition.
 *
 * Uses OpenMP parallel reduction (reduction(max:max_speed)) to compute signal
 * speeds across all zones simultaneously, then takes the maximum.
 *
 * The CFL condition ensures numerical stability:
 *   dt <= CFL_NUMBER * dx / max_signal_speed
 * where max_signal_speed = max(|v| + c) over all zones,
 * and c = sqrt(gamma * p / rho) is the local sound speed.
 *
 * @param grid Current grid state (conserved variables).
 * @return The maximum stable timestep dt. Returns 1.0 if max_signal_speed is zero
 *         (trivially at rest).
 * @note Uses OpenMP parallel reduction when compiled with OpenMP support.
 */
double calculateTimeStep(const std::vector<Conserved>& grid);

/**
 * @brief Advances the grid by one timestep using the conservative Lax-Friedrichs
 *        scheme with outflow (non-reflecting) boundaries.
 *
 * Algorithm:
 *   1. Create ghost cells at both ends by copying interior state.
 *   2. Compute interface fluxes F_{i+1/2} = F(U_i, U_{i+1}) for all interfaces.
 *   3. Apply Lax-Friedrichs update:
 *        U_i^{n+1} = 0.5*(U_{i-1}^n + U_{i+1}^n) - dt/(2*dx) * (F_{i+1/2} - F_{i-1/2})
 *   4. Boundary conditions: ghost cells retain copied interior values.
 *
 * The scheme is:
 *   - First-order accurate in space
 *   - Unconditionally stable for linear advection
 *   - Monotone (does not produce new extrema) for CFL <= 1
 *   - Robust near shocks (artificial viscosity damps oscillations)
 *
 * @param grid Current grid state (conserved variables).
 * @param dt Timestep size (must satisfy CFL condition).
 * @return Updated grid state after one timestep.
 * @note Ghost cells at i=-1 and i=N_ZONES are created by copying i=0 and i=N_ZONES-1.
 */
std::vector<Conserved> updateLaxFriedrichs(const std::vector<Conserved>& grid, double dt);

#endif // FLUIDSOLVER_NUMERICAL_H
