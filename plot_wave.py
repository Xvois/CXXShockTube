"""plot_wave.py — Time evolution and wave propagation plots."""

import numpy as np
from plot_core import load_numeric, save_plot, plt, PROJECT_ROOT

# Consistent color palette for time steps
TIME_COLOURS = [
    "#1f77b4",  # blue  - t=0
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # grey
]

TIME_STYLES = ["-", "--", "-.", ":", "-", "--", "-.", ":"]
MARKER_LIST = ["o", "s", "^", "D", "v", "<", ">", "p"]


def _density_evolution(df, x_label, filename, discontinuity_x, share_time_steps=None):
    """Plot density vs x at multiple times."""
    all_times = df["time"].unique()
    sorted_times = sorted(all_times)
    n_t = len(sorted_times)

    # Build time-to-index mapping for consistency
    if share_time_steps is not None:
        time_map = {}
        for i, t in enumerate(sorted_times):
            time_map[t] = i
    else:
        time_map = {t: i for i, t in enumerate(sorted_times)}

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1)

    for t_val in sorted_times:
        idx = time_map[t_val]
        subset = df[df["time"] == t_val]
        if len(subset) == 0:
            continue
        lw = 2.5 if t_val == sorted_times[-1] else 1.5
        ax.plot(subset["x"].values, subset["density"].values,
                linewidth=lw, linestyle=TIME_STYLES[idx % len(TIME_STYLES)],
                marker=MARKER_LIST[idx % len(MARKER_LIST)], markevery=max(1, len(subset) // 10),
                color=TIME_COLOURS[idx % len(TIME_COLOURS)],
                label=f"t={t_val:.2f}")

    ax.axvline(discontinuity_x, color="red", linestyle="--", alpha=0.6, linewidth=1.5,
               label="Initial discontinuity")

    ax.set_xlabel(f"Position {x_label}")
    ax.set_ylabel(r"Density $\rho$")
    ax.legend(loc="best", frameon=False, fontsize=8, ncol=2)
    ax.tick_params(axis='both', which='major', labelsize=9)

    save_plot(fig, filename)


def generate():
    """Generate wave propagation evolution plots."""
    _, _, full_a = load_numeric("12345_problemA_results.csv", "A")
    _, _, full_b = load_numeric("12345_problemB_results.csv", "B")

    # --- Problem A density evolution ---
    _density_evolution(
        full_a, "x",
        "12345_problemA_wave_propagation", 0.3, share_time_steps=None
    )

    # --- Problem B density evolution ---
    _density_evolution(
        full_b, "r",
        "12345_problemB_wave_propagation", 0.4, share_time_steps=None
    )

    # --- Combined side-by-side ---
    all_times_a = sorted(full_a["time"].unique())
    all_times_b = sorted(full_b["time"].unique())

    fig = plt.figure(figsize=(14, 6))

    # Problem A
    ax_a = fig.add_subplot(1, 2, 1)
    for i, t_val in enumerate(all_times_a):
        subset = full_a[full_a["time"] == t_val]
        lw = 2.5 if t_val == all_times_a[-1] else 1.5
        ax_a.plot(subset["x"].values, subset["density"].values,
                  linewidth=lw, linestyle=TIME_STYLES[i % len(TIME_STYLES)],
                  marker=MARKER_LIST[i % len(MARKER_LIST)], markevery=max(1, len(subset) // 10),
                  color=TIME_COLOURS[i % len(TIME_COLOURS)],
                  label=f"t={t_val:.2f}")
    ax_a.axvline(0.3, color="red", linestyle="--", alpha=0.5)
    ax_a.set_xlabel("Position $x$")
    ax_a.set_ylabel(r"Density $\rho$")
    ax_a.legend(loc="best", frameon=False, fontsize=8, ncol=2)
    ax_a.tick_params(axis='both', which='major', labelsize=9)

    # Problem B
    ax_b = fig.add_subplot(1, 2, 2)
    for i, t_val in enumerate(all_times_b):
        subset = full_b[full_b["time"] == t_val]
        lw = 2.5 if t_val == all_times_b[-1] else 1.5
        ax_b.plot(subset["x"].values, subset["density"].values,
                  linewidth=lw, linestyle=TIME_STYLES[i % len(TIME_STYLES)],
                  marker=MARKER_LIST[i % len(MARKER_LIST)], markevery=max(1, len(subset) // 10),
                  color=TIME_COLOURS[i % len(TIME_COLOURS)],
                  label=f"t={t_val:.2f}")
    ax_b.axvline(0.4, color="red", linestyle="--", alpha=0.5)
    ax_b.set_xlabel("Radial position $r$")
    ax_b.set_ylabel(r"Density $\rho$")
    ax_b.legend(loc="best", frameon=False, fontsize=8, ncol=2)
    ax_b.tick_params(axis='both', which='major', labelsize=9)

    save_plot(fig, "12345_combined_wave_propagation")

    print("  ✓ Wave propagation plots: problemA/B_wave_propagation, combined_wave_propagation")
