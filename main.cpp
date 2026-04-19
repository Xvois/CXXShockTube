//
// main.cpp
//
// Fluid dynamics solver for PH30110 Computational Astrophysics.
// Solves Euler equations for:
//   Problem A: Cartesian shock tube (standard test case)
//   Problem B: Spherical shock tube (supernova remnant)
//
// Uses Lax-Friedrichs scheme with operator splitting for spherical geometry.
// Boundary conditions: outflow (non-reflecting) for Problem A,
// mixed reflecting/outflow for Problem B.
//

#include <iostream>

#include "constants.h"
#include "solver.h"

int main() {
    std::cout << "Fluid Solver - PH30110 Computational Astrophysics" << std::endl;
    std::cout << "Grid: N = " << N_ZONES << ", Gamma = " << GAMMA << std::endl;
    std::cout << "Domain: [0, 1], dx = " << DX << std::endl;
    std::cout << "CFL Number = " << CFL_NUMBER << std::endl;
    std::cout << std::endl;

    // Solve Problem A (Cartesian shock tube)
    solveProblemA("12345_problemA_results.csv");
    std::cout << std::endl;

    // Solve Problem B (Spherical shock tube)
    solveProblemB("12345_problemB_results.csv");

    std::cout << std::endl << "All simulations complete." << std::endl;
    return 0;
}
