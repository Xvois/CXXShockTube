//
// Created by Sonny Parker on 17/04/2026.
//

#include "../include/analytical.h"
#include "../include/constants.h"

double calculateSoundSpeed(const Primitive& w) {
    return std::sqrt(GAMMA * w.p / w.rho);
}

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

/*
 *  Solve for the pressure in Region 3 of the shock tube,
 *  using a Newton-Raphson iteration.
 *
 */
double findStarPressure(const Primitive& L, const Primitive& R) {
    double p_star = 0.5 * (L.p + R.p); // Initial guess
    double cs_l = std::sqrt(GAMMA * L.p / L.rho);
    double cs_r = std::sqrt(GAMMA * R.p / R.rho);

    // Newton-Raphson
    for (int i = 0; i < 20; ++i) {
        double f = pressureFunction(p_star, L.p, L.rho, cs_l) +
                   pressureFunction(p_star, R.p, R.rho, cs_r) + (R.v - L.v);

        // Numerical derivative
        double eps = p_star * 0.01;
        double f_eps = pressureFunction(p_star + eps, L.p, L.rho, cs_l) +
                       pressureFunction(p_star + eps, R.p, R.rho, cs_r) + (R.v - L.v);
        double df = (f_eps - f) / eps;

        p_star -= f / df;
        if (p_star < 0) p_star = 1e-6; // Physical floor
    }
    return p_star;
}

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
        if (p_star <= L.p) { // Rarefaction
            double head = L.v - cs_l;
            double tail = v_star - cs_l * std::pow(p_star / L.p, (GAMMA - 1.0) / (2.0 * GAMMA));
            if (s < head) return L;
            if (s > tail) {
                Primitive star_l;
                star_l.p = p_star;
                star_l.v = v_star;
                star_l.rho = L.rho * std::pow(p_star / L.p, 1.0 / GAMMA);
                return star_l;
            }
            // Inside the rarefaction fan
            Primitive fan;
            fan.v = 2.0 / (GAMMA + 1.0) * (cs_l + (GAMMA - 1.0) / 2.0 * L.v + s);
            double cs = 2.0 / (GAMMA + 1.0) * (cs_l + (GAMMA - 1.0) / 2.0 * (L.v - s));
            fan.rho = L.rho * std::pow(cs / cs_l, 2.0 / (GAMMA - 1.0));
            fan.p = L.p * std::pow(cs / cs_l, 2.0 * GAMMA / (GAMMA - 1.0));
            return fan;
        } else { // Left Shock
             // (Sod Problem A typically has a left-rarefaction, but we handle shock for completeness)
             return L;
        }
    } else { // Right of contact discontinuity
        if (p_star > R.p) { // Right Shock
            double ratio = p_star / R.p;
            double shock_speed = R.v + cs_r * std::sqrt((GAMMA + 1.0) / (2.0 * GAMMA) * ratio + (GAMMA - 1.0) / (2.0 * GAMMA));
            if (s > shock_speed) return R;
            Primitive star_r;
            star_r.p = p_star;
            star_r.v = v_star;
            star_r.rho = R.rho * ((ratio + (GAMMA - 1.0) / (GAMMA + 1.0)) / (ratio * (GAMMA - 1.0) / (GAMMA + 1.0) + 1.0));
            return star_r;
        } else { // Right Rarefaction
            return R;
        }
    }
}
