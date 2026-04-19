//
// problems.cpp
//
// Implementation of Problem A (Cartesian shock tube) and Problem B
// (Spherical shock tube) with spherical source terms.
//
// Problem A: Standard 1D Cartesian shock tube with outflow boundaries.
// Problem B: Spherical shock tube with geometric source terms and mixed
//            boundary conditions (reflecting inner, outflow outer).
//
// The spherical Euler equations in conservative form are:
//   dU/dt + dF/dr = S
//
// where U = (rho, rho*v, E) is the vector of conserved variables,
// F = (rho*v, rho*v^2 + p, v*(E + p)) is the flux vector, and
// S is the geometric source vector arising from spherical symmetry.
//
// The correct source terms for spherical symmetry are derived from
// the divergence operator in spherical coordinates:
//   S = (2*rho*v/r, 2*p/r, 0)
//
// Physical interpretation:
//   - S_mass = +2*rho*v/r: mass flows through expanding spherical surface
//   - S_mom  = +2*p/r: pressure force from expanding geometry
//   - S_energy = 0: energy is conserved (work done appears via momentum)
//
// Boundary conditions: reflecting wall at inner boundary (r=0, physical symmetry)
// and outflow (non-reflecting) at outer boundary (r=1.0).
// At r = 0 (inner boundary), the source terms are singular (0/0),
// so we skip source term application at the first cell.
// At r=0, v=0 is enforced as a physical symmetry condition (no flow through origin).
// At the outer boundary, ghost cells copy the interior state to allow
// the shock to exit freely without reflection.
//

#include "problems.h"
#include "analytical.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#if OPENMP_AVAILABLE
#include <omp.h>
#endif

// =========================================================================
// Problem A: Cartesian Shock Tube
// =========================================================================

/**
 * Initializes Problem A: Cartesian shock tube.
 *
 * Discontinuity at x = 0.3:
 *   Left state (x < 0.3): rho = 1.0, v = 0.75, p = 1.0
 *   Right state (x > 0.3): rho = 0.125, v = 0.0, p = 0.1
 *
 * Boundary conditions: outflow (non-reflecting) at both ends.
 * Ghost cells copy the interior state so waves exit freely.
 *
 * Each zone is initialized with the cell-centered value.
 */
std::vector<Conserved> initialiseProblemA() {
    std::vector<Conserved> grid(N_ZONES);
    for (int z = 0; z < N_ZONES; ++z) {
        double x = X_MIN + (z + 0.5) * DX;  // Cell center
        Primitive w{};

        if (x < 0.3) {
            w.rho = 1.0;   // Left state
            w.p   = 1.0;
            w.v   = 0.75;
        } else {
            w.rho = 0.125;  // Right state
            w.p   = 0.1;
            w.v   = 0.0;
        }
        grid[z] = primToCons(w);
    }
    return grid;
}

// =========================================================================
// Problem B: Spherical Shock Tube
// =========================================================================

/**
 * Initializes Problem B: Spherical shock tube.
 *
 * Discontinuity at r = 0.4:
 *   Left state (r < 0.4): rho = 1.0, v = 0.0, p = 1.0
 *   Right state (r > 0.4): rho = 0.125, v = 0.0, p = 0.1
 *
 * Boundary conditions: mixed — reflecting wall at r=0 (physical symmetry),
 * outflow at r=1.0 (shock exits freely).
 *
 * This models the initial conditions of a supernova explosion:
 * a high-pressure core surrounded by low-density gas.
 */
std::vector<Conserved> initialiseProblemB() {
    std::vector<Conserved> grid(N_ZONES);
    for (int z = 0; z < N_ZONES; ++z) {
        double r = X_MIN + (z + 0.5) * DX;  // Cell center
        Primitive w{};

        if (r < 0.4) {
            w.rho = 1.0;   // Left state (high pressure core)
            w.p   = 1.0;
            w.v   = 0.0;
        } else {
            w.rho = 0.125;  // Right state (low pressure ambient)
            w.p   = 0.1;
            w.v   = 0.0;
        }
        grid[z] = primToCons(w);
    }
    return grid;
}

// =========================================================================
// Spherical Solver: Lax-Friedrichs with Conservative Source Terms
// =========================================================================

/**
 * Solves the spherical Euler equations using an operator-splitting scheme.
 *
 * The spherical Euler equations in conservative form are:
 *   dU/dt + dF/dr = S
 *
 * where U = (rho, rho*v, E) is the vector of conserved variables,
 * F = (rho*v, rho*v^2 + p, v*(E + p)) is the flux vector, and
 * S = (0, 2*p/r, 0) is the geometric source vector arising from spherical symmetry.
 *
 * We use Strang splitting (symmetric operator splitting):
 *   1. Apply half-step of source term: U* = U^n + (dt/2) * S(U^n)
 *   2. Apply full-step of flux update: U** = U* - dt * dF/dr (conservative LF)
 *   3. Apply remaining half-step of source term: U^{n+1} = U** + (dt/2) * S(U**)
 *   4. Apply boundary conditions: v=0 at inner boundary (physical symmetry),
 *      outflow at outer boundary (ghost cell copy).
 *
 * This is second-order accurate in time (O(dt^2)) despite first-order
 * in space (O(dx)).
 *
 * Source terms (correct for spherical symmetry):
 *   S_mass   = 0               (mass is conserved)
 *   S_mom    = +2*p/r         (pressure force from geometry)
 *   S_energy = 0               (energy is conserved)
 *
 * Boundary conditions:
 *   - Inner boundary (r=0): reflecting wall (v=0), physical symmetry.
 *   - Outer boundary (r=1.0): outflow (non-reflecting).
 * At r = 0 (first cell), source terms are skipped to avoid 0/0 singularity.
 *
 * @param grid  Current conserved variable state
 * @param dt    Timestep (must satisfy CFL condition)
 * @return      Updated conserved variable state
 */
std::vector<Conserved> updateSphericalLaxFriedrichs(const std::vector<Conserved>& grid, double dt) {
    std::vector<Conserved> next_grid = grid;
    double dtdx = dt / DX;

     // --- Step 1: Apply source term at half-time step ---
    // S = (0, 2*p/r, 0) applied at dt/2
    // Skip first cell (r ~ 0) to avoid division by zero.
    // OpenMP parallel for applies source terms to all interior cells simultaneously.
#if OPENMP_AVAILABLE
    #pragma omp parallel for
#endif
    for (int i = 1; i < N_ZONES - 1; ++i) {
        Primitive w = consToPrim(next_grid[i]);
        double r = X_MIN + (i + 0.5) * DX;
        if (r > 1e-10) {
            double source_mom = 2.0 * w.p / r * 0.5 * dt;
            next_grid[i].mom += source_mom;
        }
    }

    // --- Step 2: Compute maximum wave speed (Rusanov) using parallel reduction ---
    double max_a = 0.0;
#if OPENMP_AVAILABLE
    #pragma omp parallel for reduction(max:max_a)
#endif
    for (int i = 0; i < (int)next_grid.size(); ++i) {
        Primitive w = consToPrim(next_grid[i]);
        double c = calculateSoundSpeed(w);
        double speed = std::abs(w.v) + c;
        if (speed > max_a) max_a = speed;
    }

    // --- Step 3: Create ghost cells for outflow boundaries ---
    // Ghost cell at left (i=-1): copy cell 0
    // Ghost cell at right (i=N): copy cell N-1
    // Both boundaries use outflow — waves exit freely without reflection.
    std::vector<Conserved> ghost_grid(N_ZONES + 2);
    for (int i = 0; i < N_ZONES; ++i) {
        ghost_grid[i + 1] = next_grid[i];  // ghost_grid[1..N_ZONES] = next_grid[0..N-1]
    }
    // Left ghost cell: copy from cell 0 (outflow; v=0 enforced in Step 7)
    ghost_grid[0] = next_grid[0];
    // Right ghost cell: copy from cell N-1 (outflow — no modification)
    ghost_grid[N_ZONES + 1] = next_grid[N_ZONES - 1];

   // --- Step 4: Compute Rusanov numerical fluxes at all N_ZONES+1 interfaces ---
    // OpenMP parallel for computes fluxes at all interfaces simultaneously.
    std::vector<Conserved> F(N_ZONES + 1);
#if OPENMP_AVAILABLE
    #pragma omp parallel for
#endif
    for (int i = 0; i < (int)F.size(); ++i) {
        Primitive wL = consToPrim(ghost_grid[i]);
        Primitive wR = consToPrim(ghost_grid[i + 1]);
        Conserved fL = computeFlux(wL, ghost_grid[i]);
        Conserved fR = computeFlux(wR, ghost_grid[i + 1]);
        // Rusanov numerical flux at interface i
        F[i].mass     = 0.5 * (fL.mass + fR.mass)     - 0.5 * max_a * (ghost_grid[i + 1].mass     - ghost_grid[i].mass);
        F[i].mom      = 0.5 * (fL.mom + fR.mom)      - 0.5 * max_a * (ghost_grid[i + 1].mom      - ghost_grid[i].mom);
        F[i].energy   = 0.5 * (fL.energy + fR.energy) - 0.5 * max_a * (ghost_grid[i + 1].energy - ghost_grid[i].energy);
    }

     // --- Step 5: Update interior cells using conservative flux difference ---
    // OpenMP parallel for updates all cells simultaneously.
#if OPENMP_AVAILABLE
    #pragma omp parallel for
#endif
    for (int i = 0; i < (int)next_grid.size(); ++i) {
        next_grid[i].mass     = grid[i].mass     - dtdx * (F[i + 1].mass     - F[i].mass);
        next_grid[i].mom      = grid[i].mom      - dtdx * (F[i + 1].mom      - F[i].mom);
        next_grid[i].energy   = grid[i].energy   - dtdx * (F[i + 1].energy - F[i].energy);
    }

// OpenMP parallel for applies remaining half-step source term.
#if OPENMP_AVAILABLE
    #pragma omp parallel for
#endif
    for (int i = 1; i < N_ZONES - 1; ++i) {
        Primitive w = consToPrim(next_grid[i]);
        double r = X_MIN + (i + 0.5) * DX;
        if (r > 1e-10) {
            double source_mom = 2.0 * w.p / r * 0.5 * dt;
            next_grid[i].mom += source_mom;
        }
    }

    // --- Step 7: Apply boundary conditions ---
    // Inner boundary (r ~ 0): v=0 enforced as physical symmetry (no flow through origin).
    // Outer boundary (r ~ 1.0): outflow — no modification needed (ghost cell copy).
    //
    // Process for inner boundary only:
    //   1. Convert conserved → primitive (get v, ρ, p)
    //   2. Set v = 0 (physical symmetry at origin)
    //   3. Recompute pressure from internal energy
    //   4. Floor pressure if negative
    Primitive w0 = consToPrim(next_grid[0]);
    double E_old = next_grid[0].energy;
    w0.v = 0.0;                          // Physical symmetry: no flow through origin
    w0.p = (GAMMA - 1.0) * w0.rho * (E_old - 0.5 * w0.rho * w0.v * w0.v);
    if (w0.p < 0.0) w0.p = 0.0;           // Physical floor
    next_grid[0] = primToCons(w0);
    // Outer boundary: no modification (outflow — shock exits freely)

    return next_grid;
}
