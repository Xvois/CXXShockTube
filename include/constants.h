/**
 * @file constants.h
 * @brief Global simulation parameters for the ideal-gas Euler equations solver.
 *
 * Defines domain geometry, numerical constants, and the CFL stability number
 * for the Lax-Friedrichs finite volume scheme.
 *
 * Physical model:
 *   - Ideal gas with adiabatic exponent gamma = 1.4 (monatomic gas)
 *   - 1D domain: x ∈ [0, 1] with N_ZONES = 100 uniform zones
 *   - Zone width: dx = 0.01
 *   - CFL number: 0.5 (ensures numerical stability)
 *
 * @note This is a fixed-grid solver. The grid is not adaptive — zone count
 *       and domain size are compile-time constants.
 */
#ifndef FLUIDSOLVER_CONSTANTS_H
#define FLUIDSOLVER_CONSTANTS_H

/** Adiabatic exponent for a monatomic ideal gas (gamma = 1.4). */
const double GAMMA  = 1.4;

/** Number of computational zones across the domain. */
const int    N_ZONES = 100;

/** Domain lower bound (x_min). */
const double X_MIN   = 0.0;

/** Domain upper bound (x_max). */
const double X_MAX   = 1.0;

/** Zone width: dx = (x_max - x_min) / N_ZONES = 0.01. */
const double DX      = (X_MAX - X_MIN) / N_ZONES;

/** CFL stability number. Timestep dt = CFL * dx / max_signal_speed. */
const double CFL_NUMBER = 0.5;

#endif // FLUIDSOLVER_CONSTANTS_H
