/**
 * @file analytical.h
 * @brief Analytical solution helpers for the shock tube problems.
 *
 * Provides functions for computing the exact Riemann solution at any point
 * in space and time. The exact solution is used for verifying the numerical
 * scheme by comparing simulation results against the analytical prediction.
 *
 * The Riemann problem (step discontinuity) for the Euler equations produces
 * a self-similar solution U(x,t) = U((x - x0)/t) consisting of up to 5 waves:
 *
 *   Left state → [left wave] → Star region 1 → [contact] → Star region 2 → [right wave] → Right state
 *
 * Wave types (depending on pressure ratio p_star/p_side):
 *   - Shock wave (p_star > p_side): discontinuous jump governed by Rankine-Hugoniot
 *   - Rarefaction wave (p_star < p_side): continuous expansion fan (self-similar)
 *
 * The star region (region 3) has uniform pressure p_star, velocity v_star,
 * and density/rho determined by the wave type.
 *
 * Newton-Raphson iteration is used to solve the coupled equation for p_star.
 *
 * @note The exact solution assumes an ideal gas with gamma = 1.4.
 *       It is valid for t > 0 and provides a snapshot at any point in time.
 */
#ifndef FLUIDSOLVER_ANALYTICAL_H
#define FLUIDSOLVER_ANALYTICAL_H

#include "constants.h"
#include "physics.h"

namespace fluidsolver {

using fluidsolver::GAMMA;

/**
 * @brief Compute the local sound speed.
 *
 * For an ideal gas:
 *   c = sqrt(gamma * p / rho)
 *
 * @param w Primitive state (density, velocity, pressure).
 * @return Sound speed c.
 */
double calculateSoundSpeed(const Primitive& w);

/**
 * @brief Evaluate the pressure-difference function for left/right waves.
 *
 * Computes the integral of the momentum equation across a wave connecting
 * the initial state to the star region pressure p_star:
 *   f(p_star) = integral(p_star/p_side) (p_side / (rho*c^2)) dp * dp
 *
 * For a shock wave (p_star > p_side):
 *   f = (p_star - p_side) * sqrt(A / (p_star + B))
 *   where A = 2/(gamma+1)/rho_side, B = (gamma-1)/(gamma+1)*p_side
 *
 * For a rarefaction wave (p_star <= p_side):
 *   f = (2*c_side/(gamma-1)) * ((p_star/p_side)^((gamma-1)/(2*gamma)) - 1)
 *
 * Positive result indicates p_star > p_side (shock), negative indicates
 * p_star <= p_side (rarefaction).
 *
 * @param p_star Target star region pressure.
 * @param p_side Initial pressure on the side being evaluated.
 * @param rho_side Initial density on the side being evaluated.
 * @param cs_side Sound speed on the side being evaluated.
 * @return Pressure-difference function value.
 */
double pressureFunction(double p_star, double p_side, double rho_side, double cs_side);

/**
 * @brief Solve for the pressure in the star region (Region 3) using
 *        Newton-Raphson iteration on the two-sided Riemann problem.
 *
 * The star region pressure p_star satisfies:
 *   p_star = f_left(p_star) + f_right(p_star) + v_right - v_left = 0
 *
 * where f_left and f_right are the pressure-difference functions for
 * the left and right waves, respectively.
 *
 * Newton-Raphson iteration:
 *   p_star^{(k+1)} = p_star^{(k)} - f(p_star^{(k)}) / f'(p_star^{(k)})
 *
 * The derivative f' is computed numerically using a perturbation eps = p_star * 0.01.
 * Convergence is typically achieved in 3-5 iterations.
 *
 * @param L Left initial state (rho, v, p).
 * @param R Right initial state (rho, v, p).
 * @return Star region pressure p_star.
 * @note A physical floor of 1e-6 is applied to p_star to avoid negative pressures.
 */
double findStarPressure(const Primitive& L, const Primitive& R);

/**
 * @brief Sample the exact solution at a given position and time.
 *
 * Evaluates the self-similar Riemann solution at position x at time t:
 *   s = (x - x0) / t  (self-similar coordinate)
 *
 * Determines which region of the wave structure the point lies in:
 *   1. Left of left wave: U = U_L (undisturbed left state)
 *   2. In left wave (shock or rarefaction): jump or fan solution
 *   3. In star region (left of contact): p = p_star, v = v_star, rho determined by left wave
 *   4. In star region (right of contact): p = p_star, v = v_star, rho determined by right wave
 *   5. In right wave (shock or rarefaction): jump or fan solution
 *   6. Right of right wave: U = U_R (undisturbed right state)
 *
 * @param L Left initial state.
 * @param R Right initial state.
 * @param x_initial Discontinuity position x0.
 * @param x_current Position to evaluate.
 * @param t Time since initialization (must be > 0).
 * @return Primitive variables at the specified point.
 */
Primitive sampleExactSolution(const Primitive& L, const Primitive& R,
                              double x_initial, double x_current, double t);

}  // namespace fluidsolver

#endif // FLUIDSOLVER_ANALYTICAL_H
