/**
 * @file physics.h
 * @brief Primitive and conserved variable structs, conversions, and flux computation
 *        for the ideal-gas Euler equations.
 *
 * The Euler equations describe conservation of mass, momentum, and energy
 * for an inviscid (non-viscous) compressible fluid:
 *
 *   dU/dt + dF/dx = 0
 *
 * where U is the vector of conserved variables and F(U) is the flux vector.
 *
 * Variable representations:
 *   Primitive:   (rho, v, p)  — density, velocity, pressure
 *   Conserved:   (rho*u, rho*v, E)  — mass, momentum, total energy
 *
 * For an ideal gas with adiabatic exponent gamma:
 *   E = p/(gamma-1) + 0.5*rho*v^2   (thermal + kinetic energy)
 *   p = (gamma-1)*(E - 0.5*rho*v^2) (pressure from energy equation)
 *
 * The flux vector is:
 *   F = (rho*v, rho*v^2 + p, v*(E + p))
 *   mass flux:     rho*v
 *   momentum flux: rho*v^2 + p   (convective + pressure)
 *   energy flux:   v*(E + p)     (advective enthalpy flux)
 *
 * @note The Euler equations are hyperbolic PDEs. Characteristic speeds
 *       are v ± c (acoustic waves) and v (contact/discontinuity wave),
 *       where c = sqrt(gamma*p/rho) is the sound speed.
 */
#ifndef FLUIDSOLVER_PHYSICS_H
#define FLUIDSOLVER_PHYSICS_H

/**
 * @brief Primitive variables: density (rho), velocity (v), pressure (p).
 *
 * These are the physically measurable quantities at each grid point.
 * Used for output, boundary conditions, and analytical solutions.
 */
struct Primitive {
    /** Density (kg/m³ or arbitrary units). */
    double rho;
    /** Velocity (m/s or arbitrary units). */
    double v;
    /** Pressure (Pa or arbitrary units). */
    double p;
};

/**
 * @brief Conserved variables: mass (rho), momentum (rho*v), total energy (E).
 *
 * These quantities are conserved by the Euler equations in conservative form.
 * The finite volume scheme updates these directly to guarantee exact conservation.
 */
struct Conserved {
    /** Mass density (same as rho in primitive form). */
    double mass;
    /** Momentum density = rho * velocity. */
    double mom;
    /** Total energy density = thermal + kinetic. */
    double energy;
};

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
Conserved primToCons(const Primitive& prim);

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
 * @note If mass <= 0, returns zero state to avoid division by zero.
 */
Primitive consToPrim(const Conserved& u);

/**
 * @brief Compute the flux vector for the Euler equations.
 *
 * The flux is evaluated as F(U) using the primitive state w and conserved state u:
 *   F_mass   = u.mom             (momentum carries mass)
 *   F_mom    = u.mom * w.v + w.p (convective momentum + pressure)
 *   F_energy = w.v * (u.energy + w.p) (advective enthalpy)
 *
 * @param w Primitive state (used for pressure).
 * @param u Conserved state (used for convective fluxes).
 * @return Flux vector F = (mass_flux, momentum_flux, energy_flux).
 */
Conserved computeFlux(const Primitive& w, const Conserved& u);

#endif // FLUIDSOLVER_PHYSICS_H
