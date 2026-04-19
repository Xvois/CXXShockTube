"""plot_phase.py — Thermodynamic and phase-space diagrams."""

import numpy as np
from plot_core import load_numeric, save_plot, plt, PROJECT_ROOT

COLOURS = {"A": "#2196F3", "B": "#FF9800"}


def generate():
    """Generate phase-space and thermodynamic plots."""
    final_a, t_a, _ = load_numeric("12345_problemA_results.csv", "A")
    final_b, t_b, _ = load_numeric("12345_problemB_results.csv", "B")

    # ---- Pressure vs Density (Hugoniot-style) ----
    fig = plt.figure(figsize=(10, 5))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    ax_a.scatter(final_a["density"].values, final_a["pressure"].values,
                 s=8, c=final_a["x"].values, cmap="viridis", alpha=0.6, edgecolor="none")
    # Annotate key states
    n_states = 5
    for i in range(0, len(final_a), n_states):
        if i // n_states == 0:
            ax_a.plot(final_a["density"].values[i], final_a["pressure"].values[i],
                      "kx", markersize=8, label="Left state")
        if i // n_states == n_states - 1:
            ax_a.plot(final_a["density"].values[i], final_a["pressure"].values[i],
                      "ks", markersize=8, label="Right state")
    ax_a.set_xlabel(r"Density $\rho$")
    ax_a.set_ylabel(r"Pressure $p$")
    ax_a.set_title("Problem A: Phase Space (p, ρ)", fontweight="bold")
    ax_a.legend(fontsize=8)
    ax_a.grid(True, alpha=0.3)

    # Problem B: P vs rho
    ax_b.scatter(final_b["density"].values, final_b["pressure"].values,
                 s=8, c=final_b["x"].values, cmap="viridis", alpha=0.6, edgecolor="none")
    for i in range(0, len(final_b), n_states):
        if i // n_states == 0:
            ax_b.plot(final_b["density"].values[i], final_b["pressure"].values[i],
                      "kx", markersize=8, label="Left state")
        if i // n_states == n_states - 1:
            ax_b.plot(final_b["density"].values[i], final_b["pressure"].values[i],
                      "ks", markersize=8, label="Right state")
    ax_b.set_xlabel(r"Density $\rho$")
    ax_b.set_ylabel(r"Pressure $p$")
    ax_b.set_title("Problem B: Phase Space (p, ρ)", fontweight="bold")
    ax_b.legend(fontsize=8)
    ax_b.grid(True, alpha=0.3)

    fig.suptitle("Thermodynamic Phase Space: Pressure vs Density",
                 fontweight="bold", fontsize=13)
    save_plot(fig, "12345_phase_space_pd")

    # ---- Velocity vs Pressure ----
    fig = plt.figure(figsize=(10, 5))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    ax_a.scatter(final_a["velocity"].values, final_a["pressure"].values,
                 s=8, c=final_a["x"].values, cmap="viridis", alpha=0.6, edgecolor="none")
    ax_a.set_xlabel(r"Velocity $v$")
    ax_a.set_ylabel(r"Pressure $p$")
    ax_a.set_title("Problem A: Phase Space (v, p)", fontweight="bold")
    ax_a.grid(True, alpha=0.3)

    ax_b.scatter(final_b["velocity"].values, final_b["pressure"].values,
                 s=8, c=final_b["x"].values, cmap="viridis", alpha=0.6, edgecolor="none")
    ax_b.set_xlabel(r"Velocity $v$")
    ax_b.set_ylabel(r"Pressure $p$")
    ax_b.set_title("Problem B: Phase Space (v, p)", fontweight="bold")
    ax_b.grid(True, alpha=0.3)

    fig.suptitle("Thermodynamic Phase Space: Velocity vs Pressure",
                 fontweight="bold", fontsize=13)
    save_plot(fig, "12345_phase_space_vp")

    # ---- Energy partition ----
    # Compute kinetic and internal energy per cell
    def energy_partition(df):
        rho = df["density"].values
        v = df["velocity"].values
        p = df["pressure"].values
        dx = 0.01
        ke = 0.5 * rho * v**2 * dx  # kinetic energy density
        ie = p / (GAMMA - 1) * dx  # internal energy density
        return ke, ie, dx

    from plot_core import GAMMA
    ke_a, ie_a, dx_a = energy_partition(final_a)
    ke_b, ie_b, dx_b = energy_partition(final_b)

    fig = plt.figure(figsize=(10, 5))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    ax_a.fill_between(final_a["x"].values, 0, ke_a,
                       alpha=0.6, color="#2196F3", label="Kinetic")
    ax_a.fill_between(final_a["x"].values, ke_a, ke_a + ie_a,
                       alpha=0.6, color="#FF9800", label="Internal")
    ax_a.set_xlabel("Position x")
    ax_a.set_ylabel("Energy Density")
    ax_a.set_title("Problem A: Energy Partition", fontweight="bold")
    ax_a.legend(fontsize=9)

    ax_b.fill_between(final_b["x"].values, 0, ke_b,
                       alpha=0.6, color="#2196F3", label="Kinetic")
    ax_b.fill_between(final_b["x"].values, ke_b, ke_b + ie_b,
                       alpha=0.6, color="#FF9800", label="Internal")
    ax_b.set_xlabel("Radial position r")
    ax_b.set_ylabel("Energy Density")
    ax_b.set_title("Problem B: Energy Partition", fontweight="bold")
    ax_b.legend(fontsize=9)

    fig.suptitle("Energy Partition: Kinetic vs Internal",
                 fontweight="bold", fontsize=13)
    save_plot(fig, "12345_energy_partition")

    print("  ✓ Phase plots: phase_space_pd, phase_space_vp, energy_partition")
