"""plot_dynamics.py — Physical state snapshots (ρ, v, p, e) at final time."""

import numpy as np
from plot_core import (
    load_numeric, load_exact, save_plot,
    plt, PROJECT_ROOT
)

COLOURS = {"A": "#2196F3", "B": "#FF9800", "exact": "black"}
MARKERS = {"numerical": "-", "exact": "--"}


def _plot_panel(ax, x, num, exact, xlabel, ylim, problem_label=""):
    """Draw one subplot: density, velocity, pressure, or internal energy."""
    num_col = {"density": "density", "velocity": "velocity",
               "pressure": "pressure", "internal_energy": "internal_energy"}[xlabel]
    label = {"density": r"Density $\rho$", "velocity": r"Velocity $v$",
             "pressure": r"Pressure $p$",
             "internal_energy": r"Specific Internal Energy $e$"}[xlabel]

    ax.plot(x, num[num_col].values, COLOURS[problem_label],
            MARKERS["numerical"], linewidth=1.8,
            label=f"{problem_label} numerical")

    if exact is not None and len(exact) > 0:
        ax.plot(exact["x"].values, exact[num_col].values,
                COLOURS["exact"], MARKERS["exact"], linewidth=1.5,
                label=f"{problem_label} exact")

    ax.set_xlabel(f"Position $x$")
    ax.set_ylabel(label)
    ax.set_ylim(ylim)
    ax.legend(loc="best", frameon=True, fontsize=9)


def generate():
    """Generate Problem A and Problem B state snapshot plots."""
    # --- Problem A ---
    final_a, t_a, _ = load_numeric("12345_problemA_results.csv", "A")
    exact_a = load_exact("exact_solution.csv")
    x_a = final_a["x"].values

    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)
    rows = [[0, 0], [0, 1], [1, 0], [1, 1]]
    xlabel_map = ["density", "velocity", "pressure", "internal_energy"]
    ylims = [(0.08, 1.2), (0, 0.9), (0.05, 1.2), (0, 4.5)]

    for ax_i, (row, col), xlabel, ylim in zip(range(4), rows, xlabel_map, ylims):
        ax = fig.add_subplot(gs[row, col])
        _plot_panel(ax, x_a, final_a, exact_a, xlabel, ylim, "A")

    fig.suptitle("Problem A: Cartesian Shock Tube — Final State ($t=0.2$)",
                 fontsize=13, fontweight="bold", y=0.98)
    save_plot(fig, "12345_problemA_dynamics")

    # --- Problem B ---
    final_b, t_b, _ = load_numeric("12345_problemB_results.csv", "B")
    x_b = final_b["x"].values

    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)
    ylims_b = [(0.08, 1.2), (0, 0.85), (0.05, 1.2), (0, 4.0)]

    for ax_i, (row, col), xlabel, ylim in zip(range(4), rows, xlabel_map, ylims_b):
        ax = fig.add_subplot(gs[row, col])
        _plot_panel(ax, x_b, final_b, None, xlabel, ylim, "B")

    fig.suptitle("Problem B: Spherical Shock Tube — Final State ($t=0.25$)",
                 fontsize=13, fontweight="bold", y=0.98)
    save_plot(fig, "12345_problemB_dynamics")

    print("  ✓ Dynamics plots: problemA_dynamics, problemB_dynamics")
