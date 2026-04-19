//
// numerical.cpp
//
// Implements the conservative Lax-Friedrichs time-stepping scheme for the
// Euler equations on a 1D grid with outflow (non-reflecting) boundaries.
//
// The scheme is fully conservative: mass, momentum, and energy are conserved
// for waves that remain within the computational domain. Waves reaching the
// boundaries exit freely via ghost cell copy.
//

#include "numerical.h"
#include "analytical.h"

#include <cmath>

// =========================================================================
// CFL-limited timestep calculation
// =========================================================================

/**
 * Computes the maximum allowable timestep from the CFL condition.
 *
 * The CFL condition ensures numerical stability:
 *   dt <= CFL_NUMBER * dx / max_signal_speed
 * where max_signal_speed = max(|v| + c) over all zones,
 * and c = sqrt(gamma * p / rho) is the local sound speed.
 *
 * Returns the maximum stable dt. If max_signal_speed is zero
 * (trivially at rest), returns a safe default (1.0).
 */
double calculateTimeStep(const std::vector<Conserved>& grid) {
    double max_speed = 0.0;
    for (int i = 0; i < N_ZONES; ++i) {
        Primitive w = consToPrim(grid[i]);
        double cs = calculateSoundSpeed(w);
        double signal_speed = std::abs(w.v) + cs;
        if (signal_speed > max_speed) {
            max_speed = signal_speed;
        }
    }
    // Guard against division by zero (should never happen for valid data)
    if (max_speed < 1e-12) return 1.0;
    return CFL_NUMBER * DX / max_speed;
}

// =========================================================================
// Conservative Lax-Friedrichs update with outflow boundary conditions
// =========================================================================

/**
 * Advances the grid by one timestep using the conservative Lax-Friedrichs scheme.
 *
 * Algorithm (conservative interface flux form with outflow boundaries):
 *   1. Create ghost cells at both ends by copying interior states.
 *      Left ghost cell: identical state to cell 0.
 *      Right ghost cell: identical state to cell N-1.
 *      This implements outflow (non-reflecting) boundaries — waves that
 *      reach a boundary simply leave the domain without reflection.
 *   2. Compute primitive variables and Rusanov numerical fluxes at all interfaces.
 *      F_{i+1/2} = 0.5*(F(u_i) + F(u_{i+1})) - 0.5*max_a*(u_{i+1} - u_i)
 *   3. Update all interior cells: u_i^{n+1} = u_i^n - (dt/dx)*(F_{i+1/2} - F_{i-1/2})
 *      The ghost cell copy ensures no artificial reflection at boundaries.
 *
 * Conservation properties (outflow boundaries):
 *   - Mass: conserved (waves leaving the domain are physical, not artifacts)
 *   - Momentum: conserved (no spurious wall forces)
 *   - Energy: conserved (no artificial conversion at walls)
 *
 * @param grid  Current conserved variable state
 * @param dt    Timestep (must satisfy CFL condition)
 * @return      Updated conserved variable state
 */
std::vector<Conserved> updateLaxFriedrichs(const std::vector<Conserved>& grid, double dt) {
    std::vector<Conserved> next_grid = grid;
    double dtdx = dt / DX;

    // --- Step 1: Create ghost cells with outflow boundary conditions ---
    // Ghost cells copy the interior state, allowing waves to exit the domain
    // without reflection. When a wave reaches a boundary, the ghost cell
    // contains the same state as the adjacent interior cell, so the Rusanov
    // flux naturally carries the wave outward.
    std::vector<Conserved> ghost_grid(N_ZONES + 2);
    for (int i = 0; i < N_ZONES; ++i) {
        ghost_grid[i + 1] = grid[i];  // ghost_grid[1..N_ZONES] = grid[0..N-1]
    }
    // Left ghost cell (i=0): copy state from cell 0 → outflow
    // No modification of primitive variables — identical to interior cell.
    ghost_grid[0] = grid[0];
    // Right ghost cell (i=N_ZONES+1): copy state from cell N-1 → outflow
    // No modification of primitive variables — identical to interior cell.
    ghost_grid[N_ZONES + 1] = grid[N_ZONES - 1];

    // --- Step 2: Compute maximum wave speed (Rusanov) ---
    double max_a = 0.0;
    for (int i = 0; i < N_ZONES; ++i) {
        Primitive w = consToPrim(grid[i]);
        double c = calculateSoundSpeed(w);
        double a = std::abs(w.v) + c;
        if (a > max_a) max_a = a;
    }

    // --- Step 3: Compute Rusanov numerical fluxes at all N_ZONES+1 interfaces ---
    // F_{i}  = flux at interface between ghost_grid[i] and ghost_grid[i+1]
    // Indices: i=0  (left ghost - cell 0), i=1  (cell 0 - cell 1), ..., i=N_ZONES (cell N-1 - right ghost)
    std::vector<Conserved> F(N_ZONES + 1);
    for (int i = 0; i <= N_ZONES; ++i) {
        Primitive wL = consToPrim(ghost_grid[i]);
        Primitive wR = consToPrim(ghost_grid[i + 1]);
        Conserved fL = computeFlux(wL, ghost_grid[i]);
        Conserved fR = computeFlux(wR, ghost_grid[i + 1]);
        // Rusanov numerical flux at interface i
        F[i].mass     = 0.5 * (fL.mass + fR.mass)     - 0.5 * max_a * (ghost_grid[i + 1].mass     - ghost_grid[i].mass);
        F[i].mom      = 0.5 * (fL.mom + fR.mom)      - 0.5 * max_a * (ghost_grid[i + 1].mom      - ghost_grid[i].mom);
        F[i].energy   = 0.5 * (fL.energy + fR.energy) - 0.5 * max_a * (ghost_grid[i + 1].energy - ghost_grid[i].energy);
    }

    // --- Step 4: Update interior cells using conservative flux difference ---
    // u_i^{n+1} = u_i^n - (dt/dx) * (F_{i+1/2} - F_{i-1/2})
    // where F_{i-1/2} = F[i],  F_{i+1/2} = F[i+1]
    // Ghost cell copy at boundaries means no artificial reflection occurs.
    for (int i = 0; i < N_ZONES; ++i) {
        next_grid[i].mass     = grid[i].mass     - dtdx * (F[i + 1].mass     - F[i].mass);
        next_grid[i].mom      = grid[i].mom      - dtdx * (F[i + 1].mom      - F[i].mom);
        next_grid[i].energy   = grid[i].energy   - dtdx * (F[i + 1].energy - F[i].energy);
    }

    return next_grid;
}
