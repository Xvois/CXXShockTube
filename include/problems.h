//
// Problems.h
//
// Defines initial conditions and geometry for Problem A (Cartesian) and Problem B (Spherical).
//
#ifndef FLUIDSOLVER_PROBLEMS_H
#define FLUIDSOLVER_PROBLEMS_H

#include "physics.h"
#include "constants.h"
#include <vector>

// --- Problem A: Cartesian Shock Tube ---
// Left state (x < 0.3): rho=1, p=1, v=0.75
// Right state (x > 0.3): rho=0.125, p=0.1, v=0
// Snapshot at t = 0.2

std::vector<Conserved> initialiseProblemA();

// --- Problem B: Spherical Shock Tube ---
// Left state (x < 0.4): rho=1, p=1, v=0
// Right state (x > 0.4): rho=0.125, p=0.1, v=0
// Snapshot at t = 0.25
// Uses spherical symmetry with r-dependence

std::vector<Conserved> initialiseProblemB();

// --- Problem B: Spherical Update (with source terms) ---
// Solves the spherical Euler equations with 1/r source term
// Conserved form: U_t + F(U)_r = S(U, r)
// where S = (-2*p/r, 0, 0) for spherical symmetry

std::vector<Conserved> updateSphericalLaxFriedrichs(const std::vector<Conserved>& grid, double dt);

#endif //FLUIDSOLVER_PROBLEMS_H
