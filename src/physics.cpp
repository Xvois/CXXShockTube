/**
 * @file physics.cpp
 * @brief Implements primitive/conserved variable conversions and flux computation
 *        for the ideal-gas Euler equations (gamma = 1.4).
 *
 * The Euler equations describe conservation of mass, momentum, and energy
 * for an inviscid (non-viscous) compressible fluid:
 *   dU/dt + dF/dx = 0
 *
 * where U = (rho, rho*v, E) is the vector of conserved variables and
 * F = (rho*v, rho*v^2 + p, v*(E + p)) is the flux vector.
 *
 * For an ideal gas with adiabatic exponent gamma:
 *   E = p/(gamma-1) + 0.5*rho*v^2   (thermal + kinetic energy)
 *   p = (gamma-1)*(E - 0.5*rho*v^2)  (pressure from energy equation)
 */

#include "constants.h"
#include "physics.h"


/**
 * @brief Convert primitive variables to conserved variables.
 *
 * Applies the ideal-gas relations:
 *   mass   = rho
 *   mom    = rho * v
 *   energy = p/(gamma-1) + 0.5*rho*v^2
 *
 * @param prim Input primitive variables (rho, v, p).
 * @return Conserved variables (mass, mom, energy).
 */
/**
 * @brief Convert primitive variables to conserved variables.
 *
 * Applies the ideal-gas relations:
 *   mass   = rho
 *   mom    = rho * v
 *   energy = p/(gamma-1) + 0.5*rho*v^2
 *
 * @param prim Input primitive variables (rho, v, p).
 * @return Conserved variables (mass, mom, energy).
 */
Conserved primToCons(const Primitive& prim) {
    Conserved cons{};
    cons.mass = prim.rho;
    cons.mom = prim.rho * prim.v;
    cons.energy = prim.p / (GAMMA - 1.0) + 0.5 * prim.rho * prim.v * prim.v;
    return cons;
}

/**
 * @brief Convert conserved variables to primitive variables.
 *
 * Inverts the primToCons relations:
 *   rho = mass
 *   v   = mom / rho
 *   p   = (gamma-1) * (energy - 0.5*rho*v^2)
 *
 * @param u Input conserved variables (mass, mom, energy).
 * @return Primitive variables (rho, v, p).
 * @note Returns zero state if mass <= MIN_DENSITY to avoid division by zero.
 */
Primitive consToPrim(const Conserved& u) {
    // Check for zero or near-zero mass to avoid division by zero
    if (u.mass <= 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    
    Primitive w;
    w.rho = u.mass;
    w.v = u.mom / u.mass;
    w.p = (GAMMA - 1.0) * (u.energy - 0.5 * w.rho * w.v * w.v);
    return w;
}

/**
 * @brief Convert conserved variables to primitive variables.
 *
 * Inverts the primToCons relations:
 *   rho = mass
 *   v   = mom / rho
 *   p   = (gamma-1) * (energy - 0.5*rho*v^2)
 *
 * @param u Input conserved variables (mass, mom, energy).
 * @return Primitive variables (rho, v, p).
 * @note Returns zero state if mass <= MIN_DENSITY to avoid division by zero.
 */
Primitive consToPrim(const Conserved& u) {
    // Check for zero or near-zero mass to avoid division by zero
    if (u.mass <= 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    
    Primitive w;
    w.rho = u.mass;
    w.v = u.mom / u.mass;
    w.p = (GAMMA - 1.0) * (u.energy - 0.5 * w.rho * w.v * w.v);
    return w;
}

/**
 * @brief Compute the flux vector for the Euler equations.
 *
 * The flux is evaluated as F(U) using the primitive state w and
 * conserved state u:
 *   F_mass   = u.mom              (momentum carries mass)
 *   F_mom    = u.mom * w.v + w.p  (convective momentum + pressure)
 *   F_energy = w.v * (E + p)       (advective enthalpy flux)
 *
 * @param w Primitive state (used for pressure).
 * @param u Conserved state (used for convective fluxes).
 * @return Flux vector F = (mass_flux, momentum_flux, energy_flux).
 */
Conserved computeFlux(const Primitive& w, const Conserved& u) {
    Conserved flux{};
    flux.mass = u.mom;
    flux.mom = u.mom * w.v + w.p;
    flux.energy = w.v * (u.energy + w.p);
    return flux;
}
