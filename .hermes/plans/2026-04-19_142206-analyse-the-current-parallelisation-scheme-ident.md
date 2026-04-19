# Parallelisation Scheme Analysis: CXXShockTube

## Goal

Determine whether the current baseline parallelisation (`std::execution::seq`) is correct and sufficient for the integration techniques, and whether those techniques can actually benefit from parallelisation.

## Current Context

- **Project**: CXXShockTube — 1D Euler equation solver (Cartesian + spherical shock tubes)
- **Language**: C++20 strict mode
- **Target compilers**: GCC (any), Clang++, MSVC — must compile everywhere
- **Baseline**: `std::execution::seq` only (no `par`, no `par_unseq`, no OpenMP, no TBB)
- **Problem size**: N=100 zones (small, but correctness matters more than speed)

## Baseline Parallelisation: What's Already There

| Layer | Mechanism | Where | Algorithm |
|---|---|---|---|
| Intra-timestep | `std::execution::seq` | `numerical.cpp`, `problems.cpp` | `transform_reduce`, `for_each` across cells |
| Inter-problem | `std::async` + `launch::async` | `main.cpp`, `bench.cpp` | Problem A and B in separate OS threads |

Everything else (ghost cells, BCs, source terms) is sequential loops — trivially parallelisable but small overhead.

## Can the Integration Techniques Be Parallelised?

### Lax-Friedrichs (Cartesian — Problem A)

```
Step 1: Ghost cells — copy interior to boundaries
Step 2: Max wave speed (a = |v| + c) — reduce across all cells
Step 3: Rusanov flux F[i] at each interface — depends on old timestep values
Step 4: Cell update u_i^{n+1} — depends on F[i] and F[i+1]
```

**Parallelism analysis:**

| Step | Data dependencies | Parallelisable? |
|---|---|---|
| Ghost cells | None (copy only) | Trivially — but tiny overhead for N=100 |
| Max wave speed | **Reduce** across all cells | YES — already `transform_reduce` |
| Flux computation | Each F[i] depends on u_i, u_{i+1} from OLD timestep | YES — all fluxes independent |
| Cell update | Each u_i^{n+1} depends on F[i], F[i+1] | YES — all updates independent |

**Verdict: Embarrassingly parallel.** Steps 3 and 4 have no cross-cell data dependencies within the timestep. All fluxes can be computed simultaneously, then all updates can be applied simultaneously. The code already uses `std::for_each` across the full range — correct structure for parallelism.

### Strang Splitting (Spherical — Problem B)

```
Step 1: Half-step source term S = (0, 2p/r, 0)
Step 2: Max wave speed — reduce
Step 3: Ghost cells
Step 4: Rusanov fluxes
Step 5: Cell update
Step 6: Remaining half-step source term
Step 7: Boundary conditions
```

**Parallelism analysis:**

| Step | Data dependencies | Parallelisable? |
|---|---|---|
| Source term half-step | Each cell depends on its own (r, p) — no cross-cell deps | YES |
| Max wave speed | Reduce | YES — already `transform_reduce` |
| Ghost cells | None | Trivially |
| Flux computation | All independent (old timestep values) | YES |
| Cell update | All independent | YES |
| Source term half-step | Each cell depends on its own (r, p) — no cross-cell deps | YES |
| Boundary conditions | Single cell (cell 0) | No — but negligible |

**Verdict: Embarrassingly parallel.** The Strang splitting is sequential **between** steps (you can't do source term and flux at the same time), but **within** each step the work is independent across cells.

### The Bottleneck: Max Wave Speed Reduction

The `transform_reduce` for max wave speed (`max_a`) is the only reduction that spans the entire grid:

```cpp
double max_a = std::transform_reduce(
    std::execution::seq,
    grid.begin(), grid.end(),
    0.0,
    [](double a, double b) { return std::max(a, b); },
    [](const Conserved& zone) { /* signal speed */ }
);
```

With `std::execution::seq` this is sequential. For N=100 the overhead is tiny (~0.1ms), but the code is correct for parallel execution if someone upgrades to `par` later.

## What `std::execution::seq` Actually Does

- **Guarantees**: Sequential execution — no two iterations overlap in time
- **Compiler freedom**: GCC/Clang can still **auto-vectorise** (SIMD) the loops. The compiler knows the algorithm is parallel-safe (reduction + independent updates).
- **Runtime cost**: No thread creation overhead — exactly what we want for small N=100
- **MSVC**: `seq` is a no-op (always sequential) — correct and expected

**Bottom line**: `std::execution::seq` is the correct baseline. It:
1. Compiles on GCC, Clang, MSVC (all with C++20)
2. Produces correct results (no race conditions)
3. Allows compiler auto-vectorisation (SIMD) at compile time
4. Has zero runtime overhead for thread management

## Why `par`/`par_unseq` Are Not Suitable

| Policy | GCC 11+ | Clang 15+ | MSVC 17.2+ |
|---|---|---|---|
| `seq` | Yes (sequential) | Yes | Yes (no-op) |
| `unseq` | Yes (SIMD only) | Yes | No |
| `par` | Yes (multithread) | Yes | **No** (stub/missing) |
| `par_unseq` | Yes (multithread+SIMD) | Yes | **No** (missing) |

MSVC has no `par` or `par_unseq` support. Using them would break MSVC builds.

**Conclusion**: `std::execution::seq` is the **only** portable execution policy across all three compilers.

## Integration Technique Parallelisation Summary

| Technique | Within-step parallelism | With `seq` | With `par` |
|---|---|---|---|
| Lax-Friedrichs flux update | Embarrassingly parallel | Sequential (compiler may SIMD) | Multithreaded |
| Strang splitting source terms | Embarrassingly parallel | Sequential (compiler may SIMD) | Multithreaded |
| Max wave speed reduction | Reduce operation | Sequential | Multithreaded |
| Ghost cells / BCs | Trivial (copy) | Sequential | Unnecessary overhead |

**The integration is already structured for parallelism.** The `std::for_each` calls span the full grid with no cross-cell data dependencies. The only thing preventing actual speedup is the `seq` policy — which is intentional for portability.

## Files to Keep (No Changes Needed)

These are already correct for the baseline:

| File | What it does | Verdict |
|---|---|---|
| `src/numerical.cpp` | `transform_reduce` + `for_each` with `seq` | **Correct** |
| `src/problems.cpp` | `transform_reduce` + `for_each` with `seq` | **Correct** |
| `main.cpp` | `std::async` + `launch::async` for A/B | **Correct** |

## Risks (None with `seq`)

1. **Thread safety**: N/A — `seq` is sequential, no races possible
2. **Reduction correctness**: `std::max` is associative — correct regardless of execution policy
3. **Memory ordering**: `seq` guarantees sequential consistency — no atomics needed
4. **MSVC compatibility**: `seq` works, `par` does not — this is why `seq` is the right choice
5. **Small N=100**: Thread overhead would exceed computation time for `par` — `seq` is actually faster here

## Summary

The integration techniques are **already correctly structured** for parallel execution. The flux computation and cell update steps are embarrassingly parallel — each cell's computation is independent. With `std::execution::seq`:

- The loops run sequentially (correct, no races)
- The compiler may auto-vectorise to SIMD (GCC/Clang with `-O2` will)
- No thread creation overhead (good for small N=100)
- Compiles and runs on GCC, Clang, and MSVC

**No changes needed to the integration techniques.** The baseline (`seq`) is correct, portable, and sufficient.
