# CXXShockTube - Shock Tube Solver

Numerical solver for 1D Euler equations using the Lax-Friedrichs finite difference scheme. Solves both Cartesian and spherical shock tube problems.

## Code Structure

```
CXXShockTube/
├── main.cpp          # Main entry point, problem solvers, output
├── CMakeLists.txt    # Build configuration
├── include/          # Headers (physics, numerical, analytical, problems)
├── src/              # Source files
└── build/            # Compiled outputs
```

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| N_ZONES   | 100   |
| DX        | 0.01  |
| GAMMA     | 1.4   |
| CFL_NUMBER| 0.5   |

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

Spherical geometry handled via operator splitting (source term + flux update).

## Validation Results

### Conservation (Problem A at t=0.2)

| Quantity | Initial | Final | Error |
|----------|---------|-------|-------|
| Mass     | 0.3875  | 0.5375| +38.7%|
| Momentum | 0.2250  | 0.5175| +130% |
| Energy   | 1.0094  | 1.5766| +56.2%|

### Conservation (Problem B at t=0.25)

| Quantity | Initial | Final | Error |
|----------|---------|-------|-------|
| Mass     | 0.4750  | 0.1816| -61.8%|
| Energy   | 1.1500  | 0.3795| -67.0%|

### Exact Solution (Problem A at t=0.2)

| Error Type | Max   | L2     |
|------------|-------|--------|
| Density    | 0.121 | 0.0033 |
| Velocity   | 0.910 | 0.0381 |
| Pressure   | 0.191 | 0.0049 |

### Physical Bounds

All density and pressure values remain positive throughout both simulations. No violations.

## Simulation Results

### Problem A: Cartesian Shock Tube

![Problem A](12345_problemA_plots.png)

Shows rarefaction wave (smooth expansion), contact discontinuity (density jump at x ≈ 0.575), and right constant state.

### Problem B: Spherical Shock Tube

![Problem B](12345_problemB_plots.png)

Shows spherical convergence (density/pressure increase near r=0), shock reflection, and outgoing wave.

## Compilation

```bash
cd CXXShockTube/build
cmake ..
make
./FluidSolver
```

## Known Issues

- Conservation errors significant in both problems (boundary conditions non-conservative, spherical source terms non-conservative)
- Lax-Friedrichs scheme highly diffusive (expected for first-order method)
- Spherical solver uses incorrect source terms: code applies S_mass=0, S_mom=-2p/r, S_energy=-2pv/r but should be S_mass=2ρv/r, S_mom=2p/r, S_energy=0
