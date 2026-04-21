"""plot_phase.py — Thermodynamic and phase-space diagrams."""

import numpy as np
from plot_core import load_numeric, save_plot, plt, PROJECT_ROOT, GAMMA

COLOURS = {"A": "#2196F3", "B": "#FF9800"}


def generate():
    """Generate phase-space and thermodynamic plots."""
    final_a, t_a, _ = load_numeric("12345_problemA_results.csv", "A")
    final_b, t_b, _ = load_numeric("12345_problemB_results.csv", "B")

    # ---- Pressure vs Density (Hugoniot-style) ----
    fig = plt.figure(figsize=(10, 5))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    im_a = ax_a.scatter(final_a["density"].values, final_a["pressure"].values,
                 s=15, c=final_a["x"].values, cmap="viridis", alpha=0.8,
                 edgecolors="#FFFFFF", linewidth=0.5)
    # Annotate key states
    n_states = 5
    for i in range(0, len(final_a), n_states):
        if i // n_states == 0:
            ax_a.plot(final_a["density"].values[i], final_a["pressure"].values[i],
                      "kx", markersize=10, linewidth=2, label="Left state")
        if i // n_states == n_states - 1:
            ax_a.plot(final_a["density"].values[i], final_a["pressure"].values[i],
                      "ks", markersize=10, linewidth=2, label="Right state")
    ax_a.set_xlabel(r"Density $\rho$")
    ax_a.set_ylabel(r"Pressure $p$")
    ax_a.legend(fontsize=8, frameon=False)
    cb_a = plt.colorbar(im_a, ax=ax_a, pad=0.02)
    cb_a.set_label("Position $x$")
    fig.tight_layout()

    # Problem B: P vs rho
    im_b = ax_b.scatter(final_b["density"].values, final_b["pressure"].values,
                 s=15, c=final_b["x"].values, cmap="viridis", alpha=0.8,
                 edgecolors="#FFFFFF", linewidth=0.5)
    for i in range(0, len(final_b), n_states):
        if i // n_states == 0:
            ax_b.plot(final_b["density"].values[i], final_b["pressure"].values[i],
                      "kx", markersize=10, linewidth=2, label="Left state")
        if i // n_states == n_states - 1:
            ax_b.plot(final_b["density"].values[i], final_b["pressure"].values[i],
                      "ks", markersize=10, linewidth=2, label="Right state")
    ax_b.set_xlabel(r"Density $\rho$")
    ax_b.set_ylabel(r"Pressure $p$")
    ax_b.legend(fontsize=8, frameon=False)
    cb_b = plt.colorbar(im_b, ax=ax_b, pad=0.02)
    cb_b.set_label("Position $x$")

    save_plot(fig, "12345_phase_space_pd")

    # ---- Velocity vs Pressure ----
    fig = plt.figure(figsize=(10, 5))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    im_a = ax_a.scatter(final_a["velocity"].values, final_a["pressure"].values,
                 s=15, c=final_a["x"].values, cmap="viridis", alpha=0.8,
                 edgecolors="#FFFFFF", linewidth=0.5)
    ax_a.set_xlabel(r"Velocity $v$")
    ax_a.set_ylabel(r"Pressure $p$")
    cb_a = plt.colorbar(im_a, ax=ax_a)
    cb_a.set_label("Position $x$")

    im_b = ax_b.scatter(final_b["velocity"].values, final_b["pressure"].values,
                 s=15, c=final_b["x"].values, cmap="viridis", alpha=0.8,
                 edgecolors="#FFFFFF", linewidth=0.5)
    ax_b.set_xlabel(r"Velocity $v$")
    ax_b.set_ylabel(r"Pressure $p$")
    cb_b = plt.colorbar(im_b, ax=ax_b)
    cb_b.set_label("Position $x$")

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

    ke_a, ie_a, dx_a = energy_partition(final_a)
    ke_b, ie_b, dx_b = energy_partition(final_b)

    fig = plt.figure(figsize=(10, 5))
    ax_a = fig.add_subplot(1, 2, 1)
    ax_b = fig.add_subplot(1, 2, 2)

    # Use patterned fills for accessibility: colour + hatch patterns
    ax_a.fill_between(final_a["x"].values, 0, ke_a,
                       alpha=0.8, color="#2196F3", hatch="/",
                       label="Kinetic")
    ax_a.fill_between(final_a["x"].values, ke_a, ke_a + ie_a,
                       alpha=0.8, color="#FF9800", hatch="\\",
                       label="Internal")
    ax_a.set_xlabel("Position $x$")
    ax_a.set_ylabel("Energy Density")
    ax_a.legend(fontsize=9, frameon=False)

    ax_b.fill_between(final_b["x"].values, 0, ke_b,
                       alpha=0.8, color="#2196F3", hatch="/",
                       label="Kinetic")
    ax_b.fill_between(final_b["x"].values, ke_b, ke_b + ie_b,
                       alpha=0.8, color="#FF9800", hatch="\\",
                       label="Internal")
    ax_b.set_xlabel("Radial position $r$")
    ax_b.set_ylabel("Energy Density")
    ax_b.legend(fontsize=9, frameon=False)

    save_plot(fig, "12345_energy_partition")

    print("  ✓ Phase plots: phase_space_pd, phase_space_vp, energy_partition")
