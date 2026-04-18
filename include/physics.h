//
// Created by Sonny Parker on 17/04/2026.
//

#ifndef FLUIDSOLVER_PHYSICS_H
#define FLUIDSOLVER_PHYSICS_H

struct Primitive {
    double rho, v, p;
};

struct Conserved {
    double mass, mom, energy;
};

Conserved primToCons(const Primitive& prim);
Primitive consToPrim(const Conserved& u);
Conserved computeFlux(const Primitive& w, const Conserved& u);

#endif //FLUIDSOLVER_PHYSICS_H