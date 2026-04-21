/**
 * @file analytical.cpp
 * @brief Implements analytical solutions and helper functions for the shock tube
 *        problems: sound speed, pressure function, star-pressure iteration via
 *        Newton-Raphson, and exact solution sampling.
 *
 * Physical model:
 *   The Euler equations with a step discontinuity (Riemann problem) produce
 *   a self-similar solution U(x,t) = U((x - x0)/t) with up to 5 wave regions:
 *
 *   U_L → [wave 1] → Region 2 → [contact discontinuity] → Region 3 → [wave 2] → U_R
 *
 *   Wave 1 (left): shock if p_star > p_L, rarefaction if p_star < p_L
 *   Wave 2 (right): shock if p_star > p_R, rarefaction if p_star < p_R
 *   Star region (2 and 3): uniform p = p_star, v = v_star
 *
 * The pressure function f(p_star) integrates the momentum equation across each wave,
 * and the Newton-Raphson iteration solves f_L(p_star) + f_R(p_star) + v_R - v_L = 0.
 *
 * @note The exact solution is valid for t > 0 and assumes gamma = 1.4.
 */

#include "analytical.h"
#include "constants.h"

#include <cmath>

/**
 * @brief Compute the local sound speed.
 *
 * For an ideal gas:
 *   c = sqrt(gamma * p / rho)
 *
 * @param w Primitive state (density, velocity, pressure).
 * @return Sound speed c.
 */
double calculateSoundSpeed(const Primitive& w) {
    return std::sqrt(GAMMA * w.p / w.rho);
}

/**
 * @brief Evaluate the pressure-difference function for left/right waves.
 *
 * Computes the integral of the momentum equation across a wave connecting
 * the initial state to the star region pressure p_star.
 *
 * For a shock wave (p_star > p_side):
 *   f = (p_star - p_side) * sqrt(2 / ((gamma+1)*rho) / (p_star + (gamma-1)/(gamma+1)*p_side))
 *   This comes from the Rankine-Hugoniot relations across a shock.
 *
 * For a rarefaction wave (p_star <= p_side):
 *   f = (2*c/(gamma-1)) * ((p_star/p_side)^((gamma-1)/(2*gamma)) - 1)
 *   This comes from integrating the Riemann invariant along the fan.
 *
 * @param p_star Target star region pressure.
 * @param p_side Initial pressure on the side being evaluated.
 * @param rho_side Initial density on the side being evaluated.
 * @param cs_side Sound speed on the side being evaluated.
 * @return Pressure-difference function value.
 */
double pressureFunction(double p_star, double p_side, double rho_side, double cs_side) {
    double beta = (GAMMA - 1) / (2*GAMMA);
    if (p_star > p_side) { // Shock wave
        double A = 2.0 / ((GAMMA + 1.0) * rho_side);
        double B = (GAMMA - 1.0) / (GAMMA + 1.0) * p_side;
        return (p_star - p_side) * std::sqrt(A / (p_star + B));
    } else { // Rarefaction wave
        return (2.0 * cs_side / (GAMMA - 1.0)) * (std::pow(p_star / p_side, beta) - 1.0);
    }
}

/**
 * @brief Solve for the pressure in Region 3 of the shock tube,
 *        using a Newton-Raphson iteration.
 *
 * The star region pressure p_star satisfies the coupled equation:
 *   f_L(p_star) + f_R(p_star) + v_R - v_L = 0
 *
 * Algorithm:
 *   1. Initial guess: p_star = 0.5 * (p_L + p_R)
 *   2. Iterate Newton-Raphson (up to 20 iterations):
 *      a. Evaluate f(p_star) = f_L + f_R + v_R - v_L
 *      b. Numerical derivative: df = (f(p_star + eps) - f(p_star)) / eps
 *      c. Update: p_star -= f / df
 *      d. Physical floor: if p_star < 0, set p_star = 1e-6
 *   3. Return converged p_star
 *
 * @param L Left initial state (rho, v, p).
 * @param R Right initial state (rho, v, p).
 * @return Star region pressure p_star.
 */
double findStarPressure(const Primitive& L, const Primitive& R) {
    double p_star = 0.5 * (L.p + R.p); // Initial guess
    double cs_l = std::sqrt(GAMMA * L.p / L.rho);
    double cs_r = std::sqrt(GAMMA * R.p / R.rho);

    // Newton-Raphson iteration
    for (int i = 0; i < 20; ++i) {
        double f = pressureFunction(p_star, L.p, L.rho, cs_l) +
                   pressureFunction(p_star, R.p, R.rho, cs_r) + (R.v - L.v);

        // Numerical derivative (central difference with perturbation eps)
        double eps = p_star * 0.01;
        double f_eps = pressureFunction(p_star + eps, L.p, L.rho, cs_l) +
                       pressureFunction(p_star + eps, R.p, R.rho, cs_r) + (R.v - L.v);
        double df = (f_eps - f) / eps;

        p_star -= f / df;
        if (p_star < 0) p_star = 1e-6; // Physical floor
    }
    return p_star;
}

/**
 * @brief Sample the exact Riemann solution at a given position and time.
 *
 * Determines which wave region the point (x, t) lies in by computing the
 * self-similar coordinate s = (x - x0) / t and comparing against wave speeds.
 *
 * Algorithm:
 *   1. If t <= 0, return left or right state based on x relative to x0.
 *   2. Compute star pressure p_star using Newton-Raphson.
 *   3. Compute star velocity v_star from the Riemann invariants.
 *   4. Classify wave type (shock vs rarefaction) for left and right.
 *   5. Determine wave speeds and compare s against each wave speed.
 *   6. Return the appropriate state for the identified region.
 *
 * Region identification:
 *   s < s_shock_L:     Region 1 (left of left shock) — U_L
 *   s_shock_L < s < s_contact: Region 2 — p=p_star, v=v_star, rho from left wave
 *   s_contact < s < s_contact: Region 3 — p=p_star, v=v_star, rho from right wave
 *   s > s_shock_R:     Region 5 (right of right shock) — U_R
 *   Rarefaction waves: continuous interpolation through fan
 *
 * @param L Left initial state.
 * @param R Right initial state.
 * @param x_initial Discontinuity position x0.
 * @param x_current Position to evaluate.
 * @param t Time since initialization (must be > 0).
 * @return Primitive variables at the specified point.
 */
Primitive sampleExactSolution(const Primitive& L, const Primitive& R, double x_initial, double x_current, double t) {
    if (t <= 0) return (x_current < x_initial) ? L : R;

    double s = (x_current - x_initial) / t;
    double p_star = findStarPressure(L, R);
    double cs_l = std::sqrt(GAMMA * L.p / L.rho);
    double cs_r = std::sqrt(GAMMA * R.p / R.rho);

    // Calculate velocity in the star region
    double v_star = 0.5 * (L.v + R.v) + 0.5 * (pressureFunction(p_star, R.p, R.rho, cs_r) -
                                              pressureFunction(p_star, L.p, L.rho, cs_l));

    // Logic to determine which of the 5 regions we are in
    if (s < v_star) { // Left of contact discontinuity
        if (p_star <= L.p) {
            // Left wave is a rarefaction fan
            double s_left = L.v - cs_l;  // Head of rarefaction fan
            double s_right = v_star;     // Tail of fan (contact discontinuity)

            if (s < s_left) {
                // Before the fan — undisturbed left state
                return L;
            } else if (s < s_right) {
                // Inside the rarefaction fan — interpolate
                Primitive w{};
                double ratio = (s - s_left) / (s_right - s_left);
                w.v = L.v + ratio * (v_star - L.v);
                w.p = L.p * std::pow(1.0 - (GAMMA - 1) * std::abs(w.v - L.v) / (2 * cs_l), GAMMA / (GAMMA - 1));
                w.rho = L.rho * std::pow(w.p / L.p, 1.0 / GAMMA);
                return w;
            } else {
                // After the fan — star region
                Primitive w{};
                w.v = v_star;
                w.p = p_star;
                // Compute density from isentropic relation
                w.rho = L.rho * std::pow(p_star / L.p, 1.0 / GAMMA);
                return w;
            }
        } else {
            // Left wave is a shock
            // Compute shock speed using Rankine-Hugoniot
            double A = 2.0 / ((GAMMA + 1.0) * L.rho);
            double B = (GAMMA - 1.0) / (GAMMA + 1.0) * L.p;
            double s_shock = L.v + std::sqrt((p_star - L.p) * A / (1.0 + (p_star + B) / (L.p * (GAMMA + 1.0))));

            if (s < s_shock) {
                // Before shock — undisturbed left state
                return L;
            } else {
                // After shock — region 2 (star region left side)
                Primitive w{};
                w.v = v_star;
                w.p = p_star;
                // Density behind shock from RH relations
                w.rho = L.rho * (p_star + (GAMMA - 1.0)/(GAMMA + 1.0) * L.p) / ((GAMMA + 1.0)/(GAMMA - 1.0) * p_star - (GAMMA + 1.0)/(GAMMA - 1.0) * L.p);
                if (w.rho < 1e-12) w.rho = 1e-12; // Avoid division by zero
                return w;
            }
        }
    } else { // Right of contact discontinuity
        if (p_star <= R.p) {
            // Right wave is a rarefaction fan
            double s_right = R.v + cs_r;  // Head of right rarefaction
            double s_left = v_star;       // Tail of fan (contact discontinuity)

            if (s > s_right) {
                // After the fan — undisturbed right state
                return R;
            } else if (s > s_left) {
                // Inside the rarefaction fan — interpolate
                Primitive w{};
                double ratio = (s - s_left) / (s_right - s_left);
                w.v = R.v - ratio * (R.v - v_star);
                w.p = R.p * std::pow(1.0 - (GAMMA - 1) * std::abs(w.v - R.v) / (2 * cs_r), GAMMA / (GAMMA - 1));
                w.rho = R.rho * std::pow(w.p / R.p, 1.0 / GAMMA);
                return w;
            } else {
                // Before the fan — star region
                Primitive w{};
                w.v = v_star;
                w.p = p_star;
                w.rho = R.rho * std::pow(p_star / R.p, 1.0 / GAMMA);
                return w;
            }
        } else {
            // Right wave is a shock
            // Compute shock speed using Rankine-Hugoniot
            double A = 2.0 / ((GAMMA + 1.0) * R.rho);
            double B = (GAMMA - 1.0) / (GAMMA + 1.0) * R.p;
            double s_shock = R.v + std::sqrt((p_star - R.p) * A / (1.0 + (p_star + B) / (R.p * (GAMMA + 1.0))));

            if (s > s_shock) {
                // After shock — undisturbed right state
                return R;
            } else {
                // Before shock — region 3 (star region right side)
                Primitive w{};
                w.v = v_star;
                w.p = p_star;
                w.rho = R.rho * (p_star + (GAMMA - 1.0)/(GAMMA + 1.0) * R.p) / ((GAMMA + 1.0)/(GAMMA - 1.0) * p_star - (GAMMA + 1.0)/(GAMMA - 1.0) * R.p);
                if (w.rho < 1e-12) w.rho = 1e-12;
                return w;
            }
        }
    }
}
