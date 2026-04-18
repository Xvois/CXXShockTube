//
// Created by Sonny Parker on 17/04/2026.
//

#include "../include/constants.h"
#include "../include/physics.h"


// --- Physics Functions ---

Conserved primToCons(const Primitive& prim) {
    Conserved cons{};
    cons.mass = prim.rho;
    cons.mom = prim.rho * prim.v;
    cons.energy = prim.p / (GAMMA - 1.0) + 0.5 * prim.rho * prim.v * prim.v;
    return cons;
}

Primitive consToPrim(const Conserved& u) {
    Primitive w;
    w.rho = u.mass;
    w.v = u.mom / u.mass;
    w.p = (GAMMA - 1.0) * (u.energy - 0.5 * w.rho * w.v * w.v);
    return w;
}

Conserved computeFlux(const Primitive& w, const Conserved& u) {
    Conserved flux{};
    flux.mass = u.mom;
    flux.mom = u.mom * w.v + w.p;
    flux.energy = w.v * (u.energy + w.p);
    return flux;
}
