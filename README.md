# CXXShockTube - Shock Tube Solver

Numerical solver for the 1D Euler equations using the Lax-Friedrichs finite-difference scheme. Solves both Cartesian and spherical shock tube problems.

## Code Structure

```
CXXShockTube/
├── main.cpp          # Program entry point
├── CMakeLists.txt    # Build configuration
├── include/          # Header files
│   ├── constants.h   # Global simulation parameters
│   ├── physics.h     # Primitive/conserved variables, flux
│   ├── numerical.h   # CFL timestep, Lax-Friedrichs update
│   ├── analytical.h  # Sound speed, exact solution helpers
│   ├── problems.h    # Problem initialisation, spherical solver
│   └── solver.h      # Solver driver, CSV output
├── src/              # Implementation files
│   ├── physics.cpp
│   ├── numerical.cpp
│   ├── analytical.cpp
│   ├── problems.cpp
│   └── solver.cpp
├── generate_plots.py # Python script for report plots
└── README.md
```

## Simulation Parameters

| Parameter   | Value |
|-------------|-------|
| N_ZONES     | 100   |
| DX          | 0.01  |
| GAMMA       | 1.4   |
| CFL_NUMBER  | 0.5   |

## Problem Definitions

### Problem A: Cartesian Shock Tube

Discontinuous initial conditions at x = 0.3, target t = 0.2:

| Region | x < 0.3 (Left) | x > 0.3 (Right) |
|--------|----------------|-----------------|
| ρ      | 1.0            | 0.125           |
| v      | 0.75           | 0.0             |
| p      | 1.0            | 0.1             |

### Problem B: Spherical Shock Tube

Discontinuous initial conditions at r = 0.4, target t = 0.25:

| Region | r < 0.4 (Left) | r > 0.4 (Right) |
|--------|----------------|-----------------|
| ρ      | 1.0            | 0.125           |
| v      | 0.0            | 0.0             |
| p      | 1.0            | 0.1             |

## Numerical Method

**Lax-Friedrichs scheme** — first-order finite difference, conservative with CFL constraint.

Spherical geometry handled via Strang operator splitting (source term + flux update), giving second-order accuracy in time.

## Compilation

```bash
cd CXXShockTube
mkdir -p build && cd build
cmake ..
make
./FluidSolver
```

## Output

The solver produces:
- `12345_problemA_results.csv` — Problem A simulation data
- `12345_problemB_results.csv` — Problem B simulation data
- `exact_solution.csv` — Exact analytical solution for Problem A at t = 0.2

Use `generate_plots.py` to generate comparison plots for your report.

## Known Issues

- Conservation errors present in both problems (boundary conditions are non-conservative; spherical source terms are non-conservative)
- Lax-Friedrichs scheme is highly diffusive (expected for a first-order method)
- Spherical solver source terms: current implementation applies S_mass = 0, S_mom = 2p/r, S_energy = 0, but the correct form for spherical symmetry is S_mass = 2ρv/r, S_mom = 2p/r, S_energy = 0
