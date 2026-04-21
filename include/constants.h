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

namespace fluidsolver {

/** Adiabatic exponent for a monatomic ideal gas (gamma = 1.4). */
constexpr double GAMMA  = 1.4;

/** Number of computational zones across the domain. */
constexpr int    N_ZONES = 100;

/** Domain lower bound (x_min). */
constexpr double X_MIN   = 0.0;

/** Domain upper bound (x_max). */
constexpr double X_MAX   = 1.0;

/** Zone width: dx = (x_max - x_min) / N_ZONES = 0.01. */
constexpr double DX      = (X_MAX - X_MIN) / N_ZONES;

/** CFL stability number. Timestep dt = CFL_NUMBER * dx / max_signal_speed. */
constexpr double CFL_NUMBER = 0.5;

/** Minimum density threshold to avoid division by zero. */
constexpr double MIN_DENSITY = 1e-12;

/** Minimum pressure threshold to avoid numerical issues. */
constexpr double MIN_PRESSURE = 1e-12;

/** Small epsilon for floating-point comparisons. */
constexpr double EPSILON = 1e-10;

/** Enable OpenMP parallelization. */
#if defined(_OPENMP)
    #define OPENMP_AVAILABLE 1
#else
    #define OPENMP_AVAILABLE 0
#endif

} // namespace fluidsolver

#endif // FLUIDSOLVER_CONSTANTS_H
