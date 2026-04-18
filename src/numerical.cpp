//
// Created by Sonny Parker on 17/04/2026.
//

#include "../include/numerical.h"
#include "../include/analytical.h"


double calculateTimeStep(const std::vector<Conserved>& grid) {
    double max_speed = 0.0;
    for (int i = 0; i < N_ZONES; ++i) {
        Primitive w = consToPrim(grid[i]);
        double cs = calculateSoundSpeed(w);
        double signal_speed = std::abs(w.v) + cs;
        if (signal_speed > max_speed) max_speed = signal_speed;
    }
    return CFL_NUMBER * DX / max_speed;
}

std::vector<Conserved> updateLaxFriedrichs(const std::vector<Conserved>& grid, double dt) {
    std::vector<Conserved> next_grid = grid;
    double dtdx = dt / DX;

    for (int i = 1; i < N_ZONES - 1; ++i) {
        Primitive w_l = consToPrim(grid[i - 1]);
        Primitive w_r = consToPrim(grid[i + 1]);
        Conserved f_l = computeFlux(w_l, grid[i - 1]);
        Conserved f_r = computeFlux(w_r, grid[i + 1]);

        next_grid[i].mass = 0.5 * (grid[i-1].mass + grid[i+1].mass) - 0.5 * dtdx * (f_r.mass - f_l.mass);
        next_grid[i].mom = 0.5 * (grid[i-1].mom + grid[i+1].mom) - 0.5 * dtdx * (f_r.mom - f_l.mom);
        next_grid[i].energy = 0.5 * (grid[i-1].energy + grid[i+1].energy) - 0.5 * dtdx * (f_r.energy - f_l.energy);
    }
    next_grid[0] = next_grid[1];
    next_grid[N_ZONES - 1] = next_grid[N_ZONES - 2];
    return next_grid;
}