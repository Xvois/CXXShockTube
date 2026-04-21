"""plot_dynamics.py — Physical state snapshots (rho, v, p, e) at final time."""

import numpy as np
from plot_core import (
    load_numeric, load_exact, save_plot,
    plt, PROJECT_ROOT
)

COLOURS = {"A": "#2196F3", "B": "#FF9800", "exact": "black"}


def _plot_panel(ax, x, num, exact, ylabel, ylim, problem_label="",
                show_legend=True, show_exact=True):
    """Draw one subplot: density, velocity, pressure, or internal energy."""
    num_col = {"density": "density", "velocity": "velocity",
               "pressure": "pressure", "internal_energy": "internal_energy"}[ylabel]
    label = {"density": r"Density $\rho$", "velocity": r"Velocity $v$",
             "pressure": r"Pressure $p$",
             "internal_energy": r"Specific Internal Energy $e$"}[ylabel]

    # Fix: pass color/format as keyword args, NOT positional
    # (plot(x, y, color, linestyle) creates a phantom second series)
    ax.plot(x, num[num_col].values, color=COLOURS[problem_label],
            linestyle="-", linewidth=2, marker="o", markevery=15,
            label=f"{problem_label} numerical")

    if show_exact and exact is not None and len(exact) > 0:
        ax.plot(exact["x"].values, exact[num_col].values,
                color=COLOURS["exact"], linestyle="--", linewidth=2,
                marker="x", markevery=15,
                label=f"{problem_label} exact")

    ax.set_xlabel("Position $x$")
    ax.set_ylabel(label)
    ax.set_ylim(ylim)
    if show_legend:
        ax.legend(loc="best", frameon=False, fontsize=9)


def generate():
    """Generate Problem A and Problem B state snapshot plots."""
    # --- Problem A ---
    final_a, t_a, _ = load_numeric("12345_problemA_results.csv", "A")
    exact_a = load_exact("exact_solution.csv")
    x_a = final_a["x"].values

    ylabel_map = ["density", "velocity", "pressure", "internal_energy"]

    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)

    for i, ylabel in enumerate(ylabel_map):
        row, col = i // 2, i % 2
        ax = fig.add_subplot(gs[row, col])
        num_col = ylabel
        # Compute ylim from both numerical and exact if available
        if exact_a is not None:
            vmin = min(final_a[num_col].min(), exact_a[num_col].min())
            vmax = max(final_a[num_col].max(), exact_a[num_col].max())
        else:
            vmin = final_a[num_col].min()
            vmax = final_a[num_col].max()
        margin = (vmax - vmin) * 0.08 if vmax > vmin else 0.1
        ylim = (vmin - margin, vmax + margin)
        _plot_panel(ax, x_a, final_a, exact_a, ylabel, ylim, "A",
                    show_legend=(i == 0), show_exact=True)

    save_plot(fig, "12345_problemA_dynamics")

    # --- Problem B ---
    final_b, t_b, _ = load_numeric("12345_problemB_results.csv", "B")
    x_b = final_b["x"].values

    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)

    for i, ylabel in enumerate(ylabel_map):
        row, col = i // 2, i % 2
        ax = fig.add_subplot(gs[row, col])
        num_col = ylabel
        vmin = final_b[num_col].min()
        vmax = final_b[num_col].max()
        margin = (vmax - vmin) * 0.08 if vmax > vmin else 0.1
        ylim = (vmin - margin, vmax + margin)
        _plot_panel(ax, x_b, final_b, None, ylabel, ylim, "B",
                    show_legend=(i == 0), show_exact=False)

    save_plot(fig, "12345_problemB_dynamics")

    print("  ✓ Dynamics plots: problemA_dynamics, problemB_dynamics")
