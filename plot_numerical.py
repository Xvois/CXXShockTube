"""plot_numerical.py — Numerical diagnostics: conservation, error, CFL, performance."""

import os
import numpy as np
from plot_core import (
    load_numeric, load_exact, load_conservation, save_plot, plt, PROJECT_ROOT
)

CONS_PATHS = {
    "A": "conservation_12345_problemA.csv",
    "B": "conservation_12345_problemB.csv",
}

CONS_COLOURS = {"mass": "#2196F3", "momentum": "#FF9800", "energy": "#4CAF50"}
CONS_STYLES = {"mass": "-", "momentum": "--", "energy": "-."}


def compute_relative_drift(values):
    """Compute relative drift from initial value, with sensible normalization."""
    initial = values[0]
    # If initial value is near zero, use max absolute value in series as reference
    # to avoid dividing by tiny numbers
    ref = max(abs(initial), np.max(np.abs(values)) * 1e-3) if initial != 0 else np.max(np.abs(values)) * 1e-3
    ref = max(ref, 1e-12)  # Never divide by zero
    return (values - initial) / ref


def generate():
    """Generate numerical diagnostics plots."""
    final_a, t_a, _ = load_numeric("12345_problemA_results.csv", "A")
    final_b, t_b, _ = load_numeric("12345_problemB_results.csv", "B")
    exact_a = load_exact("exact_solution.csv")

    # ---- Conservation error (drift from initial) ----
    cons_a = load_conservation(CONS_PATHS["A"])
    cons_b = load_conservation(CONS_PATHS["B"])

    fig = plt.figure(figsize=(12, 4))
    for i, (label, cons_df) in enumerate([
        ("Problem A", cons_a), ("Problem B", cons_b)
    ]):
        ax = fig.add_subplot(1, 2, i + 1)
        for j, q in enumerate(["mass", "momentum", "energy"]):
            if q not in cons_df.columns:
                continue
            # Skip momentum for Problem B (starts at zero, drift meaningless)
            if q == "momentum" and i == 1:
                continue
            drift = compute_relative_drift(cons_df[q].values)
            marker = ["o", "s", "^"][j]
            ax.plot(cons_df["time"].values, drift, linewidth=2, color=CONS_COLOURS[q],
                    linestyle=CONS_STYLES[q], marker=marker, markevery=5,
                    label=f"{label} {q}")
        ax.axhline(0, color="black", linestyle=":", alpha=0.3, linewidth=1)
        ax.set_xlabel("Time")
        ax.set_ylabel("Relative Drift")
        ax.legend(fontsize=8, frameon=False, ncol=1)

    save_plot(fig, "12345_conservation_drift")

    # Conservation values — split into subplots for clarity
    fig = plt.figure(figsize=(12, 5))
    for qi, q in enumerate(["mass", "momentum", "energy"]):
        ax = fig.add_subplot(2, 3, qi + 1)
        if cons_a is not None and q in cons_a.columns:
            ax.plot(cons_a["time"].values, cons_a[q].values,
                    linewidth=2, linestyle="-", marker="o", markevery=5,
                    label=f"A: {q}")
        if cons_b is not None and q in cons_b.columns:
            ax.plot(cons_b["time"].values, cons_b[q].values,
                    linewidth=2, linestyle="--", marker="s", markevery=5,
                    label=f"B: {q}")
        ax.set_xlabel("Time")
        ax.set_ylabel(f"{q.title()}")
        ax.legend(fontsize=8, frameon=False)

    # Add summary plot: relative drift of conservation
    ax = fig.add_subplot(2, 3, 3)
    if cons_a is not None:
        for q in ["mass", "momentum", "energy"]:
            if q not in cons_a.columns:
                continue
            drift_a = compute_relative_drift(cons_a[q].values)
            ax.plot(cons_a["time"].values, drift_a, linewidth=2,
                    color=CONS_COLOURS[q], linestyle=CONS_STYLES[q],
                    label=f"A: {q}")
    if cons_b is not None:
        for q in ["mass", "momentum", "energy"]:
            if q not in cons_b.columns:
                continue
            # Skip momentum for Problem B (starts at zero, drift meaningless)
            if q == "momentum":
                continue
            drift_b = compute_relative_drift(cons_b[q].values)
            ax.plot(cons_b["time"].values, drift_b, linewidth=2,
                    color=CONS_COLOURS[q], linestyle=CONS_STYLES[q],
                    label=f"B: {q}")
    ax.axhline(0, color="black", linestyle=":", alpha=0.3, linewidth=1)
    ax.set_xlabel("Time")
    ax.set_ylabel("Relative Drift")
    ax.legend(fontsize=7, frameon=False, ncol=2)

    save_plot(fig, "12345_conservation_values")

    # ---- Error norms ----
    if exact_a is not None and len(exact_a) > 0:
        cols = ["density", "velocity"]
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        # Load per-time error for Problem A
        _, _, full_a = load_numeric("12345_problemA_results.csv", "A")
        all_times = sorted(full_a["time"].unique())
        l1_vals, l2_vals, times = [], [], []

        for t_val in all_times:
            tdf = full_a[full_a["time"] == t_val]
            if len(tdf) > 0:
                l1, l2 = 0.0, 0.0
                for col in cols:
                    dx = 0.01
                    err_diff = np.abs((tdf[col].values - exact_a[col].values[:len(tdf)]) * dx)
                    l1 += err_diff.sum()
                    l2 += np.sqrt(np.sum(err_diff**2 / dx))
                l1_vals.append(l1)
                l2_vals.append(l2)
                times.append(t_val)

        ax1.plot(times, l1_vals, "#2196F3", linewidth=2, marker="o", markevery=3,
                 label="L1 (density + velocity)")
        ax1.set_xlabel("Time")
        ax1.set_ylabel("L1 Error")
        ax1.legend(frameon=False)
        ax1.tick_params(axis='both', which='major', labelsize=9)

        ax2.plot(times, l2_vals, "#FF9800", linewidth=2, marker="s", markevery=3,
                 label="L2 proxy")
        ax2.set_xlabel("Time")
        ax2.set_ylabel("L2 Error")
        ax2.legend(frameon=False)
        ax2.tick_params(axis='both', which='major', labelsize=9)
        save_plot(fig, "12345_error_evolution")

    # ---- CFL timestep ----
    _, _, full_a = load_numeric("12345_problemA_results.csv", "A")
    _, _, full_b = load_numeric("12345_problemB_results.csv", "B")

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    n_steps_a, n_steps_b = 106, 116
    dt_a, dt_b = 0.2 / n_steps_a, 0.25 / n_steps_b
    t_a_lin = np.linspace(0, 0.2, n_steps_a)
    t_b_lin = np.linspace(0, 0.25, n_steps_b)
    ax.step(t_a_lin, np.full(n_steps_a, dt_a), where="mid",
            linewidth=2, color="#2196F3", linestyle="-", marker="o", markevery=10,
            label="Problem A (Cartesian)")
    ax.step(t_b_lin, np.full(n_steps_b, dt_b), where="mid",
            linewidth=2, color="#FF9800", linestyle="--", marker="s", markevery=10,
            label="Problem B (Spherical)")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"Timestep $\Delta t$")
    ax.legend(fontsize=9, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=9)
    save_plot(fig, "12345_cfl_timestep")

    # ---- Performance bar chart ----
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    labels = ["Sequential", "Parallel"]
    times = [6.79, 4.70]  # ms from performance_metrics.txt
    bars = ax.bar(labels, times,
                  color=["#607D8B", "#4CAF50"],
                  hatch=["///", "..."],
                  edgecolor="white", linewidth=0.5)
    for bar, t in zip(bars, times):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                f"{t} ms", ha="center", va="bottom", fontweight="bold", fontsize=11)
    ax.set_ylabel("Time (ms)")
    ax.set_ylim(0, 8)
    # Annotate speedup
    speedup = times[0] / times[1]
    ax.text(0.5, 6.5, f"Speedup: {speedup:.2f}\u00d7  |  Efficiency: {speedup/2*100:.0f}%",
            ha="center", fontsize=11, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    # Legend with hatch info
    ax.legend(["Sequential (///)", "Parallel (...)"], loc="upper left", fontsize=9, frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=9)
    save_plot(fig, "12345_performance")

    # ---- Numerical diffusion (shock width) ----
    if len(final_a) > 0 and len(exact_a) > 0:
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
        dx = 0.01
        dRho_dx_a = np.gradient(final_a["density"].values, dx)
        ax.plot(final_a["x"].values, np.abs(dRho_dx_a),
                linewidth=2, color="#2196F3", label=r"Problem A: $|\partial\rho/\partial x|$")
        ax.set_xlabel("Position $x$")
        ax.set_ylabel(r"$|\partial\rho/\partial x|$")
        ax.legend(fontsize=9, frameon=False)
        ax.tick_params(axis='both', which='major', labelsize=9)
        save_plot(fig, "12345_numerical_diffusion")

    print("  ✓ Numerical plots: conservation_drift, conservation_values, error_evolution, "
          "cfl_timestep, performance, numerical_diffusion")
