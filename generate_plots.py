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
    return data, final_time, df  # return full dataframe for time-lapse

def load_exact(filename):
    """Load exact solution data (handles both headerless and header CSVs)."""
    try:
        df = pd.read_csv(filename)
        if 'internal_energy' not in df.columns:
            df['internal_energy'] = df['pressure'] / (df['density'] * (GAMMA - 1.0))
        return df
    except (KeyError, AssertionError):
        pass
    # Fallback: headerless CSV
    df = pd.read_csv(filename, header=None, names=['time', 'x', 'density', 'velocity', 'pressure'])
    df['internal_energy'] = df['pressure'] / (df['density'] * (GAMMA - 1.0))
    return df


def plot_problem_a(numeric_file, exact_file, candidate="12345"):
    """Plot 1: Problem A final state vs exact solution."""
    print("Loading numeric data for Problem A...")
    num_data, final_time, _ = load_numeric(numeric_file)
    exact_data = load_exact(exact_file)
    print(f"Generating Problem A exact comparison (t = {final_time:.2f})...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Problem A: Cartesian Shock Tube — Numerical vs Exact (t = {final_time:.2f})",
                 fontsize=14, fontweight='bold')

    plot_configs = [
        ('density', 'Density', 'tab:blue', axes[0, 0]),
        ('velocity', 'Velocity', 'tab:orange', axes[0, 1]),
        ('pressure', 'Pressure', 'tab:green', axes[1, 0]),
        ('internal_energy', 'Specific Internal Energy', 'tab:red', axes[1, 1]),
    ]

    for col, title, color, ax in plot_configs:
        # Exact solution in background
        ax.plot(exact_data['x'], exact_data[col], color='black', linewidth=1.0,
                linestyle='--', alpha=0.7, label='Exact Solution', zorder=2)
        # Numerical solution on top
        ax.plot(num_data['x'], num_data[col], color=color, linewidth=2.0,
                label='Lax-Friedrichs (N=100)', zorder=3)
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.legend(loc='upper right', fontsize=9)

    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_problemA_exact_comparison.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

    output_png = f"{candidate}_problemA_exact_comparison.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")


def plot_problem_b(numeric_file, candidate="12345"):
    """Plot 2: Problem B spherical structure with annotations."""
    print("Loading numeric data for Problem B...")
    num_data, final_time, _ = load_numeric(numeric_file)
    print(f"Generating Problem B spherical structure (t = {final_time:.2f})...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Problem B: Spherical Shock Tube — Supernova Model (t = {final_time:.2f})",
                 fontsize=14, fontweight='bold')

    plot_configs = [
        ('density', 'Density', 'tab:blue', axes[0, 0]),
        ('velocity', 'Velocity', 'tab:orange', axes[0, 1]),
        ('pressure', 'Pressure', 'tab:green', axes[1, 0]),
        ('internal_energy', 'Specific Internal Energy', 'tab:red', axes[1, 1]),
    ]

    for col, title, color, ax in plot_configs:
        ax.plot(num_data['x'], num_data[col], color=color, linewidth=2.0,
                label='Spherical Solver (N=100)', zorder=2)
        # Mark discontinuity location
        ax.axvline(x=0.4, color='gray', linestyle=':', linewidth=1.5, alpha=0.7, label='Discontinuity (r=0.4)')
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.legend(loc='upper right', fontsize=9)

    # Add annotations for key features
    # Annotate origin and shock front on density plot
    ax_den = axes[0, 0]
    ax_den.annotate('Origin (r ≈ 0)', xy=(0.02, num_data['density'].min() * 0.95),
                    fontsize=9, color='gray', ha='left')
    # Find shock position (rightmost pressure peak)
    shock_idx = num_data['pressure'].idxmax()
    shock_x = num_data.loc[shock_idx, 'x']
    ax_den.annotate('Shock front', xy=(shock_x + 0.03, num_data.loc[shock_idx, 'density'] * 1.02),
                    fontsize=9, color='red', ha='left',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.0, alpha=0.7))

    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_problemB_spherical_structure.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

    output_png = f"{candidate}_problemB_spherical_structure.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")


def plot_wave_propagation(numeric_file, candidate="12345"):
    """Plot 3: Time-lapse density showing wave propagation for Problem A."""
    print("Loading full simulation data for wave propagation...")
    df = pd.read_csv(numeric_file)
    times = sorted(df['time'].unique())
    print(f"Generating wave propagation plot ({len(times)} snapshots)...")

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    fig.suptitle("Problem A: Density Evolution — Rarefaction Fan, Contact, and Shock Formation",
                 fontsize=13, fontweight='bold')

    # Color map for time series
    colors = plt.cm.viridis(np.linspace(0.2, 1.0, len(times)))
    linestyles = ['solid'] * len(times)

    for i, t in enumerate(times):
        data = df[df['time'] == t].copy()
        lw = 2.5 if i == len(times) - 1 else 1.5  # bold final timestep
        alpha = 0.7 if i < len(times) - 1 else 1.0
        ax.plot(data['x'], data['density'], color=colors[i], linewidth=lw,
                alpha=alpha, linestyle=linestyles[i],
                label=f't = {t:.2f}' if i == len(times) - 1 or i % 2 == 0 else None,
                zorder=2 if i == len(times) - 1 else 1)

    # Mark initial discontinuity
    ax.axvline(x=0.3, color='black', linestyle=':', linewidth=1.5, alpha=0.5, label='Initial discontinuity')

    # Annotate key features on final timestep
    final = df[df['time'] == times[-1]]
    # Rarefaction head: where density starts dropping
    rh_idx = int(0.15 / 0.01)  # approx x=0.15
    ax.annotate('Rarefaction\nhead', xy=(0.05, final['density'].iloc[rh_idx] * 1.15),
                fontsize=9, color='darkblue', ha='left',
                arrowprops=dict(arrowstyle='->', color='darkblue', lw=1.0, alpha=0.7))

    # Contact discontinuity: sharp density jump
    # Find the steepest gradient
    dense = final['density'].values
    dx = final['x'].values
    grad = np.abs(np.gradient(dense, dx))
    contact_idx = np.argmax(grad)
    if dx[contact_idx] > 0.3:
        ax.annotate('Contact\n', xy=(dx[contact_idx] + 0.04, dense[contact_idx]),
                    fontsize=9, color='darkgreen', ha='left',
                    arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1.0, alpha=0.7))

    # Shock front: rightmost pressure/density feature
    shock_idx = np.argmax(np.abs(np.gradient(final['pressure'].values, final['x'].values)))
    ax.annotate('Shock\nfront', xy=(final['x'].iloc[shock_idx] + 0.04, final['density'].iloc[shock_idx] * 1.1),
                fontsize=9, color='darkred', ha='left',
                arrowprops=dict(arrowstyle='->', color='darkred', lw=1.0, alpha=0.7))

    ax.set_xlabel("x", fontsize=11)
    ax.set_ylabel("Density (ρ)", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    ax.legend(loc='upper right', fontsize=9)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    output_file = f"{candidate}_wave_propagation.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

    output_png = f"{candidate}_wave_propagation.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")


def plot_error_analysis(numeric_file, exact_file, candidate="12345"):
    """Plot 4: Absolute and relative error between numerical and exact solutions."""
    print("Loading data for error analysis...")
    num_data, final_time, _ = load_numeric(numeric_file)
    exact_data = load_exact(exact_file)
    print(f"Computing error norms (t = {final_time:.2f})...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"Problem A: Numerical Error at t = {final_time:.2f} — First-Order Diffusion at Discontinuities",
                 fontsize=13, fontweight='bold')

    # --- Density error ---
    ax1 = axes[0]
    error_density = np.abs(num_data['density'].values - exact_data['density'].values)
    ax1.plot(num_data['x'].values, error_density, color='tab:blue', linewidth=2.0, zorder=2)
    ax1.fill_between(num_data['x'].values, error_density, alpha=0.3, color='tab:blue')
    ax1.set_title('Absolute Density Error', fontweight='bold')
    ax1.set_xlabel('x')
    ax1.set_ylabel('|ρ_num − ρ_exact|')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    # Print norms
    l1_rho = np.trapezoid(error_density, num_data['x'].values)
    l2_rho = np.sqrt(np.trapezoid(error_density**2, num_data['x'].values))
    linf_rho = np.max(error_density)
    ax1.text(0.98, 0.95,
             f'L1  = {l1_rho:.3e}\nL2  = {l2_rho:.3e}\nL∞ = {linf_rho:.3e}',
             transform=ax1.transAxes, fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             family='monospace')

    # --- Pressure relative error ---
    ax2 = axes[1]
    error_pressure = np.abs(num_data['pressure'].values - exact_data['pressure'].values)
    # Safe relative error: mask near-zero denominators
    mask = exact_data['pressure'].values > 1e-10
    rel_pressure = np.where(mask, error_pressure / exact_data['pressure'].values, 0.0)
    ax2.plot(num_data['x'].values, rel_pressure, color='tab:green', linewidth=2.0, zorder=2)
    ax2.fill_between(num_data['x'].values, rel_pressure, alpha=0.3, color='tab:green')
    ax2.set_title('Relative Pressure Error', fontweight='bold')
    ax2.set_xlabel('x')
    ax2.set_ylabel('|p_num − p_exact| / p_exact')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1)
    # Print norms
    l1_p = np.trapezoid(rel_pressure, num_data['x'].values)
    l2_p = np.sqrt(np.trapezoid(rel_pressure**2, num_data['x'].values))
    linf_p = np.max(rel_pressure)
    ax2.text(0.98, 0.95,
             f'L1  = {l1_p:.3e}\nL2  = {l2_p:.3e}\nL∞ = {linf_p:.3e}',
             transform=ax2.transAxes, fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             family='monospace')

    plt.tight_layout(rect=[0, 0.05, 1, 0.93])
    output_file = f"{candidate}_error_analysis.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

    output_png = f"{candidate}_error_analysis.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")

    # Print norms to stdout for report
    print("\n=== Error Norms ===")
    print(f"Density  — L1: {l1_rho:.6e}, L2: {l2_rho:.6e}, Linf: {linf_rho:.6e}")
    print(f"Pressure — L1: {l1_p:.6e}, L2: {l2_p:.6e}, Linf: {linf_p:.6e}")


def main():
    """Main entry point — generate all 4 plots."""

    A_NUMERIC_PATH = "/home/sonny/CXXShockTube/build/12345_problemA_results.csv"
    A_EXACT_PATH = "/home/sonny/CXXShockTube/build/exact_solution.csv"
    B_NUMERIC_PATH = "/home/sonny/CXXShockTube/build/12345_problemB_results.csv"
    CANDIDATE = "12345"

    print("\n" + "=" * 55)
    print("\tGENERATING COURSEWORK PLOTS")
    print("=" * 55)

    # Plot 1: Problem A exact comparison
    print("\n--- Plot 1: Problem A vs Exact ---")
    plot_problem_a(A_NUMERIC_PATH, A_EXACT_PATH, CANDIDATE)

    # Plot 2: Problem B spherical structure
    print("\n--- Plot 2: Problem B Spherical ---")
    plot_problem_b(B_NUMERIC_PATH, CANDIDATE)

    # Plot 3: Wave propagation (time-lapse)
    print("\n--- Plot 3: Wave Propagation ---")
    plot_wave_propagation(A_NUMERIC_PATH, CANDIDATE)

    # Plot 4: Error analysis
    print("\n--- Plot 4: Error Analysis ---")
    plot_error_analysis(A_NUMERIC_PATH, A_EXACT_PATH, CANDIDATE)

    print("\n" + "=" * 55)
    print("\tALL PLOTS COMPLETE")
    print("=" * 55 + "\n")


if __name__ == "__main__":
    main()
