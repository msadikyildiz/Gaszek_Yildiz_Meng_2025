"""Regenerate Figure 4 panels G (R²) and H (ρ) compact + readable for the
right-hand column of the composite, in Helvetica Neue.

Same R²/ρ definition as Figure 4E-F (Pearson r between measured and
order-K-predicted fitness; R² = r², ρ = r) at AMP 781 / AZT 36 µg/mL, but
rendered square-ish with a wrapped axis label and inside legend so the text
stays legible at ~34 mm panel width (no downscaling in the composite).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figstyle

def _repo_root(p: Path) -> Path:
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


OUT = Path(__file__).resolve().parent
PARQUET = _repo_root(OUT) / "data" / "processed" / "Epistasis_Combined.parquet"
AMP_CONC, AZT_CONC, NORDERS = 781.0, 36.0, 13
COLOR_AMP, COLOR_AZT = "gray", "hotpink"


def curve(df, drug, conc):
    sub = df[(df.Drug == drug) & (df.Concentration == conc)]
    fit = sub["Fitness"].to_numpy()
    r2, rho = [], []
    for o in range(1, NORDERS + 1):
        pred = sub[f"Fitness_predicted for order {o}"].to_numpy()
        m = np.isfinite(fit) & np.isfinite(pred)
        r = np.corrcoef(fit[m], pred[m])[0, 1]
        r2.append(r * r); rho.append(r)
    return r2, rho


df = pd.read_parquet(PARQUET)
amp_r2, amp_rho = curve(df, "AMP", AMP_CONC)
azt_r2, azt_rho = curve(df, "AZT", AZT_CONC)
orders = list(range(1, NORDERS + 1))

figstyle.apply()


def panel(amp_vals, azt_vals, ylab, out):
    # Full box, dashed grid, framed legend, y-ticks every 0.2, all 13 x-ticks,
    # bold single-line axis labels.
    fig, ax = plt.subplots(figsize=(4.05, 2.05))
    ax.grid(True, linestyle="--", linewidth=0.6, color="0.82", zorder=0)
    ax.plot(orders, amp_vals, "o-", color=COLOR_AMP, lw=2.2, markersize=5.5,
            label="Ampicillin", zorder=3)
    ax.plot(orders, azt_vals, "o-", color=COLOR_AZT, lw=2.2, markersize=5.5,
            label="Aztreonam", zorder=3)
    ax.set_xlim(0.4, NORDERS + 0.6)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xticks(range(1, NORDERS + 1))
    ax.set_xticklabels(range(1, NORDERS + 1), rotation=45)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(labelsize=11, length=3)
    ax.set_xlabel("Highest-order Epistatic Terms Included", fontsize=12,
                  fontweight="bold", labelpad=4)
    ax.set_ylabel(ylab, fontsize=15, fontweight="bold", rotation=0,
                  labelpad=20, va="center")
    ax.legend(loc="lower right", fontsize=11, frameon=True, edgecolor="0.35",
              fancybox=False, borderpad=0.5)
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", out.name)


panel(amp_r2, azt_r2, r"$\mathbf{R^2}$", OUT / "panel_G.png")
panel(amp_rho, azt_rho, r"$\boldsymbol{\rho}$", OUT / "panel_H.png")
