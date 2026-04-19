"""plot_wave.py — Time evolution and wave propagation plots."""

import numpy as np
from plot_core import load_numeric, save_plot, plt, PROJECT_ROOT


def _density_evolution(df, x_label, filename, title, problem_label,
                       discontinuity_x=0.3):
    """Plot density vs x at multiple times."""
    all_times = df["time"].unique()
    n_t = len(all_times)
    cmap = plt.get_cmap("viridis")
    colours = [cmap(i / max(n_t - 1, 1)) for i in range(n_t)]

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)

    for i, t_val in enumerate(all_times):
        subset = df[df["time"] == t_val]
        if len(subset) == 0:
            continue
        ax.plot(subset["x"].values, subset["density"].values,
                color=colours[i], linewidth=1.3,
                label=f"t={t_val:.2f}" if i % 2 == 0 or n_t <= 3 else None)

    ax.plot([discontinuity_x, discontinuity_x],
            ax.get_ylim()[0:2], "r--", linewidth=1.2, alpha=0.6,
            label=f"Initial discontinuity")

    ax.set_xlabel(f"Position {x_label}")
    ax.set_ylabel(r"Density $\rho$")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="best", fontsize=8, frameon=True)

    save_plot(fig, filename)


def generate():
    """Generate wave propagation evolution plots."""
    # Load full data for both problems
    _, _, full_a = load_numeric("12345_problemA_results.csv", "A")
    _, _, full_b = load_numeric("12345_problemB_results.csv", "B")

    # --- Problem A density evolution ---
    _density_evolution(
        full_a, "x",
        "12345_problemA_wave_propagation",
        r"Problem A: Density Evolution ($\gamma=1.4$, Cartesian)",
        "A", discontinuity_x=0.3
    )

    # --- Problem B density evolution ---
    _density_evolution(
        full_b, "r",
        "12345_problemB_wave_propagation",
        r"Problem B: Density Evolution ($\gamma=1.4$, Spherical)",
        "B", discontinuity_x=0.4
    )

    # --- Combined side-by-side ---
    all_times_a = sorted(full_a["time"].unique())
    all_times_b = sorted(full_b["time"].unique())

    fig = plt.figure(figsize=(14, 6))
    cmap = plt.get_cmap("viridis")

    # Problem A
    ax_a = fig.add_subplot(1, 2, 1)
    n_a = len(all_times_a)
    for i, t_val in enumerate(all_times_a):
        subset = full_a[full_a["time"] == t_val]
        c = cmap(i / max(n_a - 1, 1))
        lw = 2.5 if t_val == all_times_a[-1] else 1.2
        ax_a.plot(subset["x"].values, subset["density"].values,
                  color=c, linewidth=lw,
                  label=f"t={t_val:.2f}" if i % 2 == 0 else None)
    ax_a.axvline(0.3, color="red", linestyle="--", alpha=0.5)
    ax_a.set_xlabel("Position x")
    ax_a.set_ylabel(r"Density $\rho$")
    ax_a.set_title("Problem A: Cartesian", fontweight="bold")
    ax_a.legend(loc="best", fontsize=8)

    # Problem B
    ax_b = fig.add_subplot(1, 2, 2)
    n_b = len(all_times_b)
    for i, t_val in enumerate(all_times_b):
        subset = full_b[full_b["time"] == t_val]
        c = cmap(i / max(n_b - 1, 1))
        lw = 2.5 if t_val == all_times_b[-1] else 1.2
        ax_b.plot(subset["x"].values, subset["density"].values,
                  color=c, linewidth=lw,
                  label=f"t={t_val:.2f}" if i % 2 == 0 else None)
    ax_b.axvline(0.4, color="red", linestyle="--", alpha=0.5)
    ax_b.set_xlabel("Radial position r")
    ax_b.set_ylabel(r"Density $\rho$")
    ax_b.set_title("Problem B: Spherical", fontweight="bold")
    ax_b.legend(loc="best", fontsize=8)

    fig.suptitle("Wave Propagation: Density Evolution", fontweight="bold", fontsize=13)
    save_plot(fig, "12345_combined_wave_propagation")

    print("  ✓ Wave propagation plots: problemA/B_wave_propagation, combined_wave_propagation")
