"""plot_comparison.py — Cross-problem comparisons (A vs B, numerical vs exact)."""

import numpy as np
from plot_core import load_numeric, load_exact, save_plot, plt, PROJECT_ROOT

COLOURS = {"A": "#2196F3", "B": "#FF9800", "exact": "black"}

COL_LABEL_MAP = {
    "density": r"Density $\rho$",
    "pressure": r"Pressure $p$",
    "velocity": r"Velocity $v$",
}


def _compute_ylim(values_a, values_b, exact_values, margin_frac=0.1):
    """Compute safe y-limits from data with a margin."""
    all_vals = np.concatenate([values_a, values_b])
    if exact_values is not None and len(exact_values) > 0:
        all_vals = np.concatenate([all_vals, exact_values])
    vmin, vmax = all_vals.min(), all_vals.max()
    margin = (vmax - vmin) * margin_frac if vmax > vmin else 0.05
    return (vmin - margin, vmax + margin)


def _plot_comparison(ax, num_a, num_b, exact_df, col, xlabel="Position $x$"):
    """Plot a side-by-side comparison of two problems with optional exact."""
    label = COL_LABEL_MAP.get(col, f"{col.title()}")

    # Compute ylim from all data
    y_vals = np.concatenate([num_a[col].values, num_b[col].values])
    if exact_df is not None and len(exact_df) > 0:
        y_vals = np.concatenate([y_vals, exact_df[col].values])
    vmin, vmax = y_vals.min(), y_vals.max()
    margin = (vmax - vmin) * 0.08 if vmax > vmin else 0.05
    ax.set_ylim(vmin - margin, vmax + margin)

    ax.plot(num_a["x"].values, num_a[col].values,
            COLOURS["A"], linewidth=2, linestyle="-", marker="o", markevery=20,
            label="Problem A (Cartesian)")
    ax.plot(num_b["x"].values, num_b[col].values,
            COLOURS["B"], linewidth=2, linestyle="-", marker="s", markevery=20,
            label="Problem B (Spherical)")
    if exact_df is not None and len(exact_df) > 0:
        ax.plot(exact_df["x"].values, exact_df[col].values,
                COLOURS["exact"], linewidth=2, linestyle="--", marker="x", markevery=20,
                label="Problem A exact")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(label)
    ax.legend(loc="best", frameon=False)


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
    _plot_comparison(ax, num_a, num_b, exact_df, "density")
    save_plot(fig, "12345_density_comparison")

    # --- Pressure comparison ---
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    _plot_comparison(ax, num_a, num_b, exact_df, "pressure")
    save_plot(fig, "12345_pressure_comparison")

    # --- Velocity comparison ---
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    _plot_comparison(ax, num_a, num_b, exact_df, "velocity")
    save_plot(fig, "12345_velocity_comparison")

    print("  ✓ Comparison plots: density/pressure/velocity_comparison")
