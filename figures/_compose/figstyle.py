"""Shared matplotlib style for publication-figure composition.

Typography prefers Helvetica Neue and falls back to Helvetica, Arial, then
DejaVu Sans, so composition is portable across machines. Nature-scale type sizes
keep every composed figure consistent. Import and call apply() before composing.
"""
from __future__ import annotations

import matplotlib as mpl

FONT = "Helvetica Neue"

# Panel-letter style (Nature: lowercase bold, top-left of each panel).
LETTER = dict(fontsize=10, fontweight="bold", family=FONT, va="top", ha="left",
              color="#000000")


def apply() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 600,
        "figure.dpi": 150,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# Nature column widths (mm)
WIDTH_SINGLE_MM = 88.0
WIDTH_DOUBLE_MM = 180.0
MM = 1.0 / 25.4  # mm -> inches
