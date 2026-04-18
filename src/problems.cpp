//
// Problems.cpp
//
// Implementation of Problem A (Cartesian shock tube) and Problem B (Spherical shock tube)
// with spherical source terms for handling r-dependence.
//

#include "../include/problems.h"
#include "../include/analytical.h"
#include <iostream>
#include <cmath>

// =========================================================================
// Problem A: Cartesian Shock Tube
// =========================================================================

std::vector<Conserved> initialiseProblemA() {
    std::vector<Conserved> grid(N_ZONES);
    for (int z = 0; z < N_ZONES; ++z) {
        double x = X_MIN + (z + 0.5) * DX;
        Primitive w{};
        if (x < 0.3) {
            w.rho = 1.0;
            w.p = 1.0;
            w.v = 0.75;
        } else {
            w.rho = 0.125;
            w.p = 0.1;
            w.v = 0.0;
        }
        grid[z] = primToCons(w);
    }
    return grid;
}

// =========================================================================
// Problem B: Spherical Shock Tube
// =========================================================================

std::vector<Conserved> initialiseProblemB() {
    std::vector<Conserved> grid(N_ZONES);
    for (int z = 0; z < N_ZONES; ++z) {
        double r = X_MIN + (z + 0.5) * DX;
        Primitive w{};
        if (r < 0.4) {
            w.rho = 1.0;
            w.p = 1.0;
            w.v = 0.0;
        } else {
            w.rho = 0.125;
            w.p = 0.1;
            w.v = 0.0;
        }
        grid[z] = primToCons(w);
    }
    return grid;
}

// =========================================================================
// Spherical Solver: Lax-Friedrichs with Source Terms
// =========================================================================

/*
 * Solves the spherical Euler equations using a split Lax-Friedrichs scheme.
 * 
 * The spherical Euler equations in conservation form are:
 *   ∂U/∂t + ∂F/∂r = S,
 * where U = (ρ, ρv, E) and S = (-2p/r, 0, 0) represents the geometric source term.
 *
 * We use a simple operator splitting approach:
 *   1. Apply source term S at half-time step
 *   2. Apply flux update (Lax-Friedrichs) at full-time step  
 *   3. Apply source term S at remaining half-time step
 *
 * For the source term, we apply it in non-conservative form for simplicity,
 * adding/subtracting pressure force to momentum and internal energy.
 */

std::vector<Conserved> updateSphericalLaxFriedrichs(const std::vector<Conserved>& grid, double dt) {
    std::vector<Conserved> next_grid = grid;
    double dtdx = dt / DX;
    
    // Apply source term at half-time step (modify in-place)
    // S = (-2*p/r, 0, 0) affects momentum and internal energy
    for (int i = 1; i < N_ZONES - 1; ++i) {  // Skip cell 0 (inner boundary)
        Primitive w = consToPrim(grid[i]);
        double r = X_MIN + (i + 0.5) * DX;
        double dt_half = 0.5 * dt;
        
        // Apply source term: dp/dt = -2p/r, d(ρv)/dt = 0, dE/dt = -2pv/r
        if (r > 1e-10) {  // Avoid division by zero
            double source_mom = -2.0 * w.p / r * dt_half;
            double source_E = -2.0 * w.p * w.v / r * dt_half;
            
            next_grid[i].mom += source_mom;
            next_grid[i].energy += source_E;
        }
    }
    
    // Apply Lax-Friedrichs flux update
    for (int i = 1; i < N_ZONES - 1; ++i) {
        Primitive w_l = consToPrim(next_grid[i - 1]);
        Primitive w_r = consToPrim(next_grid[i + 1]);
        Conserved f_l = computeFlux(w_l, next_grid[i - 1]);
        Conserved f_r = computeFlux(w_r, next_grid[i + 1]);
        
        next_grid[i].mass = 0.5 * (next_grid[i-1].mass + next_grid[i+1].mass) 
                           - 0.5 * dtdx * (f_r.mass - f_l.mass);
        next_grid[i].mom = 0.5 * (next_grid[i-1].mom + next_grid[i+1].mom) 
                          - 0.5 * dtdx * (f_r.mom - f_l.mom);
        next_grid[i].energy = 0.5 * (next_grid[i-1].energy + next_grid[i+1].energy) 
                             - 0.5 * dtdx * (f_r.energy - f_l.energy);
    }
    
    // Apply source term at remaining half-time step
    for (int i = 1; i < N_ZONES - 1; ++i) {
        Primitive w = consToPrim(next_grid[i]);
        double r = X_MIN + (i + 0.5) * DX;
        double dt_half = 0.5 * dt;
        
        if (r > 1e-10) {
            double source_mom = -2.0 * w.p / r * dt_half;
            double source_E = -2.0 * w.p * w.v / r * dt_half;
            
            next_grid[i].mom += source_mom;
            next_grid[i].energy += source_E;
        }
    }
    
    // Boundary conditions (simple reflection)
    // Inner boundary (r = 0): mirror the first cell
    next_grid[0] = next_grid[1];
    // Outer boundary: mirror the last cell  
    next_grid[N_ZONES - 1] = next_grid[N_ZONES - 2];
    
    return next_grid;
}
