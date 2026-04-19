//
// physics.h
//
// Defines primitive and conserved variable structs, and their conversions,
// along with the flux computation for the ideal-gas Euler equations.
//
#ifndef FLUIDSOLVER_PHYSICS_H
#define FLUIDSOLVER_PHYSICS_H

/**
 * Primitive variables: density (rho), velocity (v), pressure (p)
 */
struct Primitive {
    double rho, v, p;
};

/**
 * Conserved variables: mass (rho), momentum (rho*v), total energy (E)
 */
struct Conserved {
    double mass, mom, energy;
};

/** Convert primitive → conserved variables */
Conserved primToCons(const Primitive& prim);

/** Convert conserved → primitive variables */
Primitive consToPrim(const Conserved& u);

/** Compute the flux vector given primitive and conserved state */
Conserved computeFlux(const Primitive& w, const Conserved& u);

#endif // FLUIDSOLVER_PHYSICS_H
