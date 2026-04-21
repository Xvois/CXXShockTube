"""plot_core.py — Shared utilities for the CXXShockTube plotting system."""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.grid": False,
    "axes.axisbelow": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 2,
    "lines.markersize": 8,
    "figure.dpi": 100,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
    # Ensure RGB (not palette) PNG output
    "image.cmap": "viridis",
    "savefig.dpi": 300,
})

GAMMA = 1.4
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))


def find_csv(*names):
    """
    Find CSV files in the data directory.

    Searches in this order:
      1. Directories relative to the script (PROJECT_ROOT)
      2. Parent directories up the tree (up to 3 levels)
      3. Common subdirectories: build/, output/, data/

    Returns a dict: {filename: absolute_path} for all found files.
    Unfound files are omitted (plots will skip them gracefully).
    """
    candidates = [PROJECT_ROOT]

    # Search up the tree
    base = PROJECT_ROOT
    for _ in range(3):
        parent = os.path.dirname(base)
        if parent == base:
            break
        candidates.append(parent)
        base = parent

    # Common subdirectories in each candidate
    subdirs = ["", "build", "output", "data", "results"]

    found = {}
    for cand in candidates:
        for sub in subdirs:
            search_dir = os.path.join(cand, sub) if sub else cand
            if not os.path.isdir(search_dir):
                continue
            for name in names:
                if name in found:
                    continue
                path = os.path.join(search_dir, name)
                if os.path.isfile(path):
                    found[name] = path
                    names = tuple(n for n in names if n != name)

    return found


def load_numeric(filename, problem_label=""):
    """Load a problem results CSV. Returns (final_df, final_time, full_df)."""
    files = find_csv(filename)
    if not files:
        raise FileNotFoundError(
            f"Cannot find '{filename}' — run the solver first.\n"
            f"Search paths checked relative to: {PROJECT_ROOT}"
        )
    path = files[filename]
    df = pd.read_csv(path)
    t_max = df["time"].max()
    final_df = df[df["time"] == t_max].copy()
    return final_df, t_max, df


def load_exact(filename="exact_solution.csv"):
    """Load the exact solution CSV (Problem A only). Returns DataFrame or None."""
    files = find_csv(filename)
    if not files:
        return None
    path = files[filename]
    return pd.read_csv(path)


def load_conservation(filename):
    """Load a conservation diagnostics CSV. Returns DataFrame or None."""
    files = find_csv(filename)
    if not files:
        return None
    path = files[filename]
    return pd.read_csv(path)


def compute_l1_error(num_df, exact_df, col):
    """L1 error = sum(|num - exact| * dx) at final time."""
    dx = 0.01  # N=100, domain [0,1]
    if exact_df is None:
        return None
    common_x = np.sort(np.intersect1d(num_df["x"].values, exact_df["x"].values))
    if len(common_x) == 0:
        return None
    num = num_df[num_df["x"].isin(common_x)].sort_values("x")
    exact = exact_df[exact_df["x"].isin(common_x)].sort_values("x")
    err = np.abs((num[col].values - exact[col].values) * dx)
    return float(err.sum())


def compute_l2_error(num_df, exact_df, col):
    """L2 error = sqrt(sum((num - exact)^2 * dx)) at final time."""
    dx = 0.01
    if exact_df is None:
        return None
    common_x = np.sort(np.intersect1d(num_df["x"].values, exact_df["x"].values))
    if len(common_x) == 0:
        return None
    num = num_df[num_df["x"].isin(common_x)].sort_values("x")
    exact = exact_df[exact_df["x"].isin(common_x)].sort_values("x")
    err = np.sqrt(np.sum(((num[col].values - exact[col].values) ** 2) * dx))
    return float(err)


def compute_error_norms(final_df, exact_df, cols):
    """Compute L1 and L2 norms for given columns."""
    norms = {}
    for col in cols:
        norms[col + "_l1"] = compute_l1_error(final_df, exact_df, col)
        norms[col + "_l2"] = compute_l2_error(final_df, exact_df, col)
    return norms


def save_plot(fig, name, fmt="pdf"):
    """Save plot in both pdf and png formats.

    PNG is always saved first to avoid Agg canvas corruption from
    the PDF backend modifying the figure state.

    Uses explicit RGBA conversion to avoid palette (color_type=8) PNG issue.
    """
    exts = [fmt]
    if fmt == "pdf":
        exts = ["png", "pdf"]  # PNG first, then PDF

    for ext in exts:
        path = os.path.join(PROJECT_ROOT, f"{name}.{ext}")
        if ext == "png":
            # Force RGBA to avoid palette PNGs (color_type=8)
            fig.set_facecolor("white")
            for child in fig.get_children():
                if hasattr(child, 'set_facecolor'):
                    child.set_facecolor("white")
            fig.savefig(path, dpi=300, format="png", transparent=False,
                        facecolor="white", edgecolor="none")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())

    plt.close(fig)
