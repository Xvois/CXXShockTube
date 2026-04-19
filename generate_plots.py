"""generate_plots.py — Orchestrator for the CXXShockTube plotting system.

Entry point that imports from the modular plotting packages and generates
all figures in order. Run from the project root:

    python generate_plots.py

Output: PDF + PNG plots in the project root directory.
"""

import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

BANNER = """
╔═══════════════════════════════════════════════════════════╗
║    CXXShockTube — Comprehensive Plotting System          ║
║    PH30110 Computational Astrophysics                    ║
║    Lax-Friedrichs Euler Solver, N=100, γ=1.4            ║
╚═══════════════════════════════════════════════════════════╝
"""


def print_banner():
    print(BANNER)
    print(f"Working directory: {PROJECT_ROOT}")
    print(f"{'─' * 60}")


def main():
    print_banner()

    # Import all plotting modules
    import plot_core, plot_dynamics, plot_comparison
    import plot_numerical, plot_wave, plot_phase

    # Step 1: Dynamics — physical state snapshots
    print("\n[1/5] Generating dynamics (state snapshots)...")
    plot_dynamics.generate()

    # Step 2: Wave propagation — time evolution
    print("[2/5] Generating wave propagation plots...")
    plot_wave.generate()

    # Step 3: Cross-problem comparisons
    print("[3/5] Generating comparison plots...")
    plot_comparison.generate()

    # Step 4: Numerical diagnostics
    print("[4/5] Generating numerical diagnostics...")
    plot_numerical.generate()

    # Step 5: Phase-space / thermodynamic diagrams
    print("[5/5] Generating phase-space diagrams...")
    plot_phase.generate()

    print(f"\n{'═' * 60}")
    print("All plots generated successfully.")
    print("Run 'ls *.pdf *.png' to see the output files.")
    print(f"{'═' * 60}\n")


if __name__ == "__main__":
    main()
