import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import sys

# Setup for PDF output
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
    # Handle both 5-column (legacy) and 6-column (new) formats
    if 'internal_energy' not in data.columns:
        data['internal_energy'] = data['pressure'] / (data['density'] * (GAMMA - 1.0))
    return data, final_time

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
    """Generate plots for Problem A."""
    print(f"Loading numeric data from {numeric_file}...")
    num_data, final_time = load_numeric(numeric_file)
    
    has_exact = os.path.exists(exact_file)
    if has_exact:
        print(f"Loading exact data from {exact_file}...")
        exact_data = load_exact(exact_file)
    else:
        print("Warning: exact solution file not found")
    
    print(f"Generating plots for Problem A (t = {final_time:.2f})...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Shock Tube Problem A (t = {final_time:.2f}) "
                 f"Candidate: {candidate}",
                 fontsize=14, fontweight='bold')
    
    plot_configs = [
        ('density', 'Density', 'tab:blue', axes[0, 0]),
        ('velocity', 'Velocity', 'tab:orange', axes[0, 1]),
        ('pressure', 'Pressure', 'tab:green', axes[1, 0]),
        ('internal_energy', 'Specific Internal Energy', 'tab:red', axes[1, 1]),
    ]
    
    for col, title, color, ax in plot_configs:
        # Plot Exact Solution in background
        if has_exact:
            ax.plot(exact_data['x'], exact_data[col], 
                   color='#D3D3D3', linewidth=5, 
                   label='Exact Solution', zorder=1)
        
        # Plot Numerical Solution on top
        ax.plot(num_data['x'], num_data[col], 
               color=color, linewidth=1.5, 
               label='Lax-Friedrichs (N=100)', zorder=2)
        
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.legend(loc='upper right', fontsize=9)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_problemA_plots.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    # Also save as PNG
    output_png = f"{candidate}_problemA_plots.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")
    plt.close()

def plot_problem_b(numeric_file, candidate="12345"):
    """Generate plots for Problem B (spherical)."""
    print(f"Loading Problem B data from {numeric_file}...")
    num_data, final_time = load_numeric(numeric_file)
    
    print(f"Generating Problem B plots (t = {final_time:.2f})...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Spherical Shock Tube Problem B (t = {final_time:.2f}) "
                 f"Candidate: {candidate}",
                 fontsize=14, fontweight='bold')
    
    plot_configs = [
        ('density', 'Density', 'tab:blue', axes[0, 0]),
        ('velocity', 'Velocity', 'tab:orange', axes[0, 1]),
        ('pressure', 'Pressure', 'tab:green', axes[1, 0]),
        ('internal_energy', 'Specific Internal Energy', 'tab:red', axes[1, 1]),
    ]
    
    for col, title, color, ax in plot_configs:
        ax.plot(num_data['x'], num_data[col], 
               color=color, linewidth=1.5, 
               label='Spherical Solver (N=100)', zorder=2)
        
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.legend(loc='upper right', fontsize=9)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_problemB_plots.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    output_png = f"{candidate}_problemB_plots.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")
    plt.close()

def plot_convergence(numeric_file_100, numeric_file_400, candidate="12345"):
    """Plot convergence comparison (N=100 vs N=400)."""
    print("Loading convergence data...")
    num_100, _ = load_numeric(numeric_file_100)
    num_400, _ = load_numeric(numeric_file_400)
    
    print("Generating convergence plots...")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"Convergence Study: N=100 vs N=400 (t = {num_100['time'].max():.2f}) "
                 f"Candidate: {candidate}",
                 fontsize=14, fontweight='bold')
    
    for idx, (col, title) in enumerate([('density', 'Density'), 
                                          ('velocity', 'Velocity'), 
                                          ('pressure', 'Pressure')]):
        ax = axes[idx]
        ax.plot(num_100['x'], num_100[col], 
               color='tab:blue', linewidth=2, 
               label='N=100', alpha=0.7, zorder=2)
        ax.plot(num_400['x'], num_400[col], 
               color='tab:red', linewidth=1.5, 
               label='N=400', alpha=0.9, zorder=3)
        
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.legend(fontsize=9)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_convergence_plots.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    output_png = f"{candidate}_convergence_plots.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_png}")
    plt.close()

def plot_all_snapshots(numeric_file, candidate="12345"):
    """Plot all time snapshots for visual inspection."""
    print(f"Loading full simulation data...")
    df = pd.read_csv(numeric_file)
    times = sorted(df['time'].unique())
    
    print(f"Generating snapshot plots for {len(times)} time steps...")
    
    n_plots = len(times)
    ncols = 4
    nrows = (n_plots + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))
    axes = axes.flatten() if n_plots > 1 else [axes]
    
    fig.suptitle(f"All Snapshots (Candidate: {candidate})",
                 fontsize=14, fontweight='bold')
    
    for idx, t in enumerate(times):
        ax = axes[idx]
        data = df[df['time'] == t].copy()
        ax.plot(data['x'], data['density'], color='blue', linewidth=1.5)
        ax.set_title(f"t = {t:.2f}", fontsize=9)
        ax.set_xlabel("x")
        ax.set_ylabel("ρ")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
    
    # Hide unused subplots
    for idx in range(len(times), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    output_file = f"{candidate}_all_snapshots.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_comparison():
    """Main entry point - generate all plots."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate plots for fluid dynamics report")
    parser.add_argument("--numeric", type=str, default="12345_problemA_results.csv",
                       help="Numeric simulation output")
    parser.add_argument("--exact", type=str, default="exact_solution.csv",
                       help="Exact solution output")
    parser.add_argument("--candidate", type=str, default="12345",
                       help="Candidate number for filenames")
    parser.add_argument("--problem-b", type=str, default=None,
                       help="Problem B numeric output")
    parser.add_argument("--convergence", type=str, default=None,
                       help="Convergence data (N=400)")
    args = parser.parse_args()
    
    # Problem A plots
    plot_problem_a(args.numeric, args.exact, args.candidate)
    
    # Problem B plots (if available)
    if args.problem_b:
        plot_problem_b(args.problem_b, args.candidate)
    
    # Convergence plots (if available)
    if args.convergence:
        plot_convergence(args.numeric, args.convergence, args.candidate)
    
    # All snapshots
    plot_all_snapshots(args.numeric, args.candidate)

if __name__ == "__main__":
    plot_comparison()
