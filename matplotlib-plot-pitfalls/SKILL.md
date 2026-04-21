---
name: matplotlib-plot-pitfalls
description: Common matplotlib pitfalls — phantom lines from plot() positional args, palette PNG corruption, tick label rendering issues, legend duplication bugs
category: data-science
---

# Matplotlib Plot Pitfalls

Common subtle matplotlib bugs that produce phantom lines, duplicate legends, blank plots, or invisible tick labels.

## Critical Bug: Phantom Lines from `plot()` Positional Args

### Problem
Calling `ax.plot(x, y, hex_color, linestyle_str, ...)` with the linestyle as a **4th positional argument** causes matplotlib to interpret it as a **second data series**:
```python
ax.plot([0, 1, 2], [0, 1, 2], "#2196F3", "-")  # BUG
```

Matplotlib parses this as: `plot(xdata1, ydata1, color, xdata2)` where:
- Series 1: `x=[0,1,2], y=[0,1,2], color=#2196F3` ✓
- Series 2: `x=[0.], y=['-'], color=default` ✗ (phantom line!)

### Symptoms
- **Duplicate legend entries** (same label twice)
- **Off-chart blue marker at (0,0)** — the phantom line's single point
- **Two lines in `ax.get_lines()`** when only one was intended
- Phantom line data: `x=[0.]`, `y=['-']` (the linestyle string itself)

### Fix
Always use **keyword arguments** for `color` and `linestyle`:
```python
ax.plot(x, y, color="#2196F3", linestyle="-", linewidth=2, marker="o")  # ✓
ax.plot(x, y, color="#2196F3", ls="-", ...)  # ✓ (short form)
ax.plot(x, y, c="#2196F3", linestyle="-", ...)  # ✓
```

### Debugging Checklist
```python
print(f"Lines: {len(ax.get_lines())}")  # Should equal number of data series
for j, line in enumerate(ax.get_lines()):
    print(f"  Line {j}: x={line.get_xdata()}, y={line.get_ydata()}, "
          f"label='{line.get_label()}'")
# If y=['-'], x=[0.], color=default → you hit this bug
```

## Palette PNG Corruption (color_type=8)

### Problem
PNGs save with `color_type=8` (indexed/palette) instead of RGB/RGBA, appearing as blank white boxes.

### Root Cause
Matplotlib Agg backend sometimes chooses a palette mode when figure transparency or `rcParams` interfere with RGBA conversion.

### Fix in `save_plot()`
```python
def save_plot(fig, name, fmt="pdf"):
    exts = ["png", fmt] if fmt == "pdf" else [fmt]  # PNG first
    for ext in exts:
        path = os.path.join(PROJECT_ROOT, f"{name}.{ext}")
        if ext == "png":
            fig.set_facecolor("white")
            for child in fig.get_children():
                if hasattr(child, 'set_facecolor'):
                    child.set_facecolor("white")
            fig.savefig(path, dpi=300, format="png",
                        transparent=False, facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight",
                        facecolor=fig.get_facecolor())
    plt.close(fig)
```

Key rules:
- **PNG first** — PDF backend corrupts Agg canvas state
- **Explicit `format="png"`** — forces Agg renderer
- **`transparent=False, facecolor="white"`** — prevents palette mode
- **Always call `plt.close(fig)`** — prevents memory leaks

## Tick Label Rendering Issues

### Problem
Y-axis tick marks appear as dashes but no numerical values render.

### Causes
1. `ax.yaxis.set_ticklabels()` called without a fixed `FixedLocator` → `UserWarning` + invisible labels
2. Tick formatter returns empty strings for some ticks
3. Figure save (`bbox_inches="tight"`) crops labels outside the figure bounds

### Fix
```python
# ❌ Don't do this (causes warnings + invisible labels):
ax.yaxis.set_ticklabels([f'{v:.2f}' for v in ax.get_yticks()])

# ✓ Do this instead — let matplotlib handle formatting:
# If you need custom formatting, use FormatStrFormatter:
from matplotlib.ticker import FormatStrFormatter
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
```

For cropped labels during save:
```python
plt.tight_layout()  # Before savefig
fig.savefig(path, dpi=300, bbox_inches="tight")
```

## Legend Duplication Patterns

### Pattern 1: Phantom line from plot() positional args
See section above.

### Pattern 2: Multiple plot() calls with same label
```python
ax.plot(x1, y1, label="data")  # First call
ax.plot(x2, y2, label="data")  # Second call with same label → duplicate!
```
Fix: Only call `ax.plot()` once per data series. Remove duplicate plotting.

### Pattern 3: `ax.legend()` called multiple times
Each call creates new legend entries. Call once per axes.

### Pattern 4: Line objects with identical `get_label()`
Check:
```python
for line in ax.get_lines():
    print(line.get_label())  # Look for duplicates
```

## Unified Verification Script

After regenerating plots, run this to catch common issues:
```python
from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
import matplotlib

matplotlib.use('Agg')

def verify_plot(path):
    img = Image.open(path)
    arr = np.array(img.convert("RGB"))
    non_white = np.sum(~((arr[:,:,0] > 240) & (arr[:,:,1] > 240) & (arr[:,:,2] > 240)))
    ratio = non_white / (arr.shape[0] * arr.shape[1]) * 100
    print(f"{path}: {ratio:.1f}% non-white pixels")
    assert ratio > 1.0, f"Plot is mostly blank ({ratio:.1f}% non-white)"

# Check all PNGs
for png in ["problemA_dynamics.png", "problemB_dynamics.png"]:
    verify_plot(png)
```

## Quick Reference: `ax.plot()` Correct Usage

```python
# ✓ Correct — all format params as keywords:
ax.plot(x, y, color="tab:blue", linestyle="-", linewidth=2,
        marker="o", markersize=8, markevery=15, label="data")

# ✓ Correct — shorthand color:
ax.plot(x, y, c="tab:blue", ls="-", lw=2, m="o", me=15, label="data")

# ✓ Correct — hex color:
ax.plot(x, y, color="#FF9800", linestyle="--", label="exact")

# ✗ WRONG — positional format string:
ax.plot(x, y, "#FF9800", "--", ...)  # Creates phantom line!

# ✗ WRONG — mixing positional and keyword:
ax.plot(x, y, "#FF9800", linestyle="--", ...)  # Also creates phantom line!
```

## Where to Check
- Any `plot_dynamics.py` or state snapshot scripts
- Any script calling `ax.plot()` with color hex codes as positional args
- Plot generation with multiple output formats (PDF + PNG)
- Scripts using `set_ticklabels()` without a FixedLocator
- Any plot with unexpected blue markers at (0,0) or off-chart points

## Related Skills
- **matplotlib-pdf-fix** — Blank PNG generation fix
- **cfd-plot-generator** — CFD/shock tube plotting (references this skill)
- **unified-plot-generator** — Consolidating multiple plots into one
- **modular-plotting-architecture** — Refactoring monolithic plot scripts
