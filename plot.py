import pandas as pd
import matplotlib.pyplot as plt
import os

# Configuration
GAMMA = 1.4
NUMERIC_FILE = "cmake-build-debug/12345_problemA_results.csv"
EXACT_FILE = "cmake-build-debug/exact_solution.csv"

def plot_comparison():
    # 1. Load Numeric Data
    if not os.path.exists(NUMERIC_FILE):
        print(f"Error: {NUMERIC_FILE} not found.")
        return
    df_num = pd.read_csv(NUMERIC_FILE)

    # Get the final snapshot from the numeric simulation
    final_time = df_num['time'].max()
    num_data = df_num[df_num['time'] == final_time].copy()
    num_data['internal_energy'] = num_data['pressure'] / (num_data['density'] * (GAMMA - 1.0))

    # 2. Load Exact Data (if available)
    has_exact = os.path.exists(EXACT_FILE)
    if has_exact:
        exact_data = pd.read_csv(EXACT_FILE)
        # Usually exact solutions are saved at a specific time, but let's check
        if 'time' in exact_data.columns:
            exact_data = exact_data[exact_data['time'] == exact_data['time'].max()].copy()
        exact_data['internal_energy'] = exact_data['pressure'] / (exact_data['density'] * (GAMMA - 1.0))
    else:
        print("Warning: exact_solution.csv not found. Plotting numeric only.")

    # 3. Create 2x2 Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Shock Tube Problem A (t = {final_time:.2f})\nNumeric vs. Exact Solution", fontsize=16)

    plot_configs = [
        ('density', 'Density', 'blue', axes[0, 0]),
        ('velocity', 'Velocity', 'red', axes[0, 1]),
        ('pressure', 'Pressure', 'green', axes[1, 0]),
        ('internal_energy', 'Specific Internal Energy', 'purple', axes[1, 1])
    ]

    for col, title, color, ax in plot_configs:
        # Plot Exact Solution in background
        if has_exact:
            ax.plot(exact_data['x'], exact_data[col], color='#D3D3D3',
                    linewidth=5, label='Exact Solution', zorder=1)

        # Plot Numerical Solution on top
        ax.plot(num_data['x'], num_data[col], color=color,
                linewidth=1.5, label='Lax-Friedrichs', zorder=2)

        ax.set_title(title)
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)

    # Add a legend to the first plot
    axes[0, 0].legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    plot_comparison()