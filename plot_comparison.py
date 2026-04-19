"""plot_comparison.py — Cross-problem comparisons (A vs B, numerical vs exact)."""

import numpy as np
from plot_core import load_numeric, load_exact, save_plot, plt, PROJECT_ROOT

COLOURS = {"A": "#2196F3", "B": "#FF9800", "exact": "black"}


def generate():
    """Generate cross-problem comparison plots."""
    final_a, t_a, _ = load_numeric("12345_problemA_results.csv", "A")
    final_b, t_b, _ = load_numeric("12345_problemB_results.csv", "B")
    exact_a = load_exact("exact_solution.csv")

    x_common = np.sort(np.intersect1d(final_a["x"].values, final_b["x"].values))
    num_a = final_a[final_a["x"].isin(x_common)].sort_values("x")
    num_b = final_b[final_b["x"].isin(x_common)].sort_values("x")
    exact_df = exact_a[exact_a["x"].isin(x_common)].sort_values("x") if exact_a is not None else None

    # --- Density comparison ---
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(num_a["x"].values, num_a["density"].values,
            COLOURS["A"], linewidth=2, label="Problem A (Cartesian)")
    ax.plot(num_b["x"].values, num_b["density"].values,
            COLOURS["B"], linewidth=2, label="Problem B (Spherical)")
    if exact_df is not None and len(exact_df) > 0:
        ax.plot(exact_df["x"].values, exact_df["density"].values,
                COLOURS["exact"], linewidth=1.5, linestyle="--", label="Problem A exact")
    ax.set_xlabel("Position x")
    ax.set_ylabel(r"Density $\rho$")
    ax.set_title("Density Comparison at Final Time", fontweight="bold")
    ax.legend(loc="best")
    ax.set_ylim(0, 1.2)
    save_plot(fig, "12345_density_comparison")

    # --- Pressure comparison ---
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(num_a["x"].values, num_a["pressure"].values,
            COLOURS["A"], linewidth=2, label="Problem A (Cartesian)")
    ax.plot(num_b["x"].values, num_b["pressure"].values,
            COLOURS["B"], linewidth=2, label="Problem B (Spherical)")
    if exact_df is not None and len(exact_df) > 0:
        ax.plot(exact_df["x"].values, exact_df["pressure"].values,
                COLOURS["exact"], linewidth=1.5, linestyle="--", label="Problem A exact")
    ax.set_xlabel("Position x")
    ax.set_ylabel(r"Pressure $p$")
    ax.set_title("Pressure Comparison at Final Time", fontweight="bold")
    ax.legend(loc="best")
    ax.set_ylim(0, 1.2)
    save_plot(fig, "12345_pressure_comparison")

    # --- Velocity comparison ---
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(num_a["x"].values, num_a["velocity"].values,
            COLOURS["A"], linewidth=2, label="Problem A (Cartesian)")
    ax.plot(num_b["x"].values, num_b["velocity"].values,
            COLOURS["B"], linewidth=2, label="Problem B (Spherical)")
    if exact_df is not None and len(exact_df) > 0:
        ax.plot(exact_df["x"].values, exact_df["velocity"].values,
                COLOURS["exact"], linewidth=1.5, linestyle="--", label="Problem A exact")
    ax.set_xlabel("Position x")
    ax.set_ylabel(r"Velocity $v$")
    ax.set_title("Velocity Comparison at Final Time", fontweight="bold")
    ax.legend(loc="best")
    ax.set_ylim(0, 0.9)
    save_plot(fig, "12345_velocity_comparison")

    print("  ✓ Comparison plots: density/pressure/velocity_comparison")
