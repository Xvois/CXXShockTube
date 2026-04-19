import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import sys
import numpy as np

# Setup for PDF output
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['axes.titlesize'] = 12
matplotlib.rcParams['axes.labelsize'] = 10

GAMMA = 1.4


def load_numeric(filename):
    """Load numeric simulation data and return final snapshot."""
    df = pd.read_csv(filename)
    final_time = df['time'].max()
    data = df[df['time'] == final_time].copy()
    if 'internal_energy' not in data.columns:
        data['internal_energy'] = data['pressure'] / (data['density'] * (GAMMA - 1.0))
    return data, final_time, df


def load_exact(filename):
    """Load exact solution data if available."""
    if not os.path.exists(filename):
        return None
    try:
        df = pd.read_csv(filename)
        if 'internal_energy' not in df.columns:
            df['internal_energy'] = df['pressure'] / (df['density'] * (GAMMA - 1.0))
        return df
    except (KeyError, AssertionError):
        return None


def plot_comparison(numeric_file, exact_file, label, solver_label, exact_label, candidate="12345"):
    """Plot 1-panel: 2x2 grid of density, velocity, pressure, internal energy."""
    print(f"Loading data for {label}...")
    num_data, final_time, _ = load_numeric(numeric_file)
    exact_data = load_exact(exact_file) if exact_file else None

    has_exact = exact_data is not None and len(exact_data) > 0
    print(f"Generating {label} (t = {final_time:.2f}, exact = {has_exact})...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"{label} — Numerical vs Exact (t = {final_time:.2f})",
                 fontsize=14, fontweight='bold')

    plot_configs = [
        ('density', 'Density', 'tab:blue', axes[0, 0]),
        ('velocity', 'Velocity', 'tab:orange', axes[0, 1]),
        ('pressure', 'Pressure', 'tab:green', axes[1, 0]),
        ('internal_energy', 'Specific Internal Energy', 'tab:red', axes[1, 1]),
    ]

    for col, title, color, ax in plot_configs:
        if has_exact:
            ax.plot(exact_data['x'], exact_data[col], color='black', linewidth=1.0,
                    linestyle='--', alpha=0.7, label=exact_label, zorder=2)
        ax.plot(num_data['x'], num_data[col], color=color, linewidth=2.0,
                label=solver_label, zorder=3)
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.legend(loc='upper right', fontsize=9)

    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_{label.lower().replace(' ', '_').replace('(', '').replace(')', '')}.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    output_png = output_file.replace('.pdf', '.png')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")
    plt.close()


def main():
    """Main entry point — generate identical plots for Problem A and Problem B."""

    A_NUMERIC_PATH = "/home/sonny/CXXShockTube/build/12345_problemA_results.csv"
    A_EXACT_PATH = "/home/sonny/CXXShockTube/build/exact_solution.csv"
    B_NUMERIC_PATH = "/home/sonny/CXXShockTube/build/12345_problemB_results.csv"
    B_EXACT_PATH = None
    CANDIDATE = "12345"

    print("\n" + "=" * 55)
    print("\tGENERATING COURSEWORK PLOTS")
    print("=" * 55)

    # Problem A: Cartesian shock tube
    print("\n--- Problem A: Cartesian ---")
    plot_comparison(A_NUMERIC_PATH, A_EXACT_PATH,
                    "Problem A", "Lax-Friedrichs (N=100)", "Exact Solution", CANDIDATE)

    # Problem B: Spherical shock tube
    print("\n--- Problem B: Spherical ---")
    plot_comparison(B_NUMERIC_PATH, B_EXACT_PATH,
                    "Problem B", "Spherical Solver (N=100)", "Exact Solution", CANDIDATE)

    print("\n" + "=" * 55)
    print("\tALL PLOTS COMPLETE")
    print("=" * 55 + "\n")


if __name__ == "__main__":
    main()
