"""Emit Figure 3's panels as vector PDFs (selectable text, vector line-art;
dense density/scatter layers rasterized in place at 600 dpi).

Reproduces the exact plotting of ``src/figure_scripts/figure_3a.py`` (panel A,
the 2x3 concentration-grid fitness comparison) and
``src/05_epistasis_figures.ipynb`` cells 11/12/22/23 (C/D fitness-vs-order
scatters, B sequence logos), reading the same precomputed parquet, so the
vector panels reproduce the manuscript figure.

Outputs go to ``figures/_compose/vector/fig3_[A,Bamp,Bazt,C,D].pdf`` for
``compose_figure3_vector.py`` to assemble.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))

import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

WT_SEQ = "LQMERMAGERTRN"
SIGN_POS = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
DEAD_PROFILE, WT_PROFILE = "X" * 13, "." * 13


def _repo_root(p: Path) -> Path:
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


HERE = Path(__file__).resolve().parent
ROOT = _repo_root(HERE)
OUT = HERE / "vector"
OUT.mkdir(exist_ok=True)

df = pd.read_parquet(ROOT / "data/processed/Epistasis_Combined.parquet")
amp_stats = pd.read_parquet(ROOT / "data/processed/amp_auc_long_df.parquet")
azt_stats = pd.read_parquet(ROOT / "data/processed/azt_auc_long_df.parquet")
AMP = df[(df.Drug == "AMP") & (df.Concentration == 781.0)].reset_index(drop=True)
AZT = df[(df.Drug == "AZT") & (df.Concentration == 36.0)].reset_index(drop=True)


def panel_A(out: Path) -> None:
    """2x3 grid: AMP(195,781) x AZT(36,108,324) fitness density (figure_3a.py)."""
    amp_concs, azt_concs = [195.0, 781.0], [36.0, 108.0, 324.0]
    colors = {"dead": "#0072B2", "wt": "#000000", "tested": "#FF0000", "clinical": "#FFFF00"}
    handles = [
        Line2D([0], [0], marker="o", ls="none", color=colors["dead"], label="Dead variant",
               markerfacecolor="none", ms=14, mew=3.5),
        Line2D([0], [0], marker="o", ls="none", color=colors["wt"], label="Wild-type",
               markerfacecolor="none", ms=14, mew=3.5),
        Line2D([0], [0], marker="D", ls="none", color=colors["tested"], label="Tested variants",
               ms=12, mec="k", mew=1),
        Line2D([0], [0], marker="^", ls="none", color=colors["clinical"], label="Clinical isolates",
               ms=13, mec="k", mew=1),
    ]
    extra = {"c1.1": "LQMKNTAGKRTRN", "c1.2": "PQMKNMAGKRMRN", "c1.3": "LKMKSMAGKRMRN",
             "c2.1": "LQMKNMAGKRMRN", "c2.2": "LKMKNMAGKRTRN", "c2.3": "LKMKNTAGKRTRN",
             "c3.1": "PQMKNTAGKRTRN"}
    try:
        clin = pd.read_csv(ROOT / "data/known_variants/encoded_variants.csv")[
            "Encoded_Sequence"].dropna().unique()
    except FileNotFoundError:
        clin = []

    rc = {"font.family": "sans-serif", "font.size": 24, "axes.titlesize": 23,
          "axes.labelsize": 24, "xtick.labelsize": 24, "ytick.labelsize": 24,
          "legend.fontsize": 24, "axes.linewidth": 1.2, "lines.linewidth": 1.5,
          "grid.color": "0.5", "grid.linestyle": "--", "grid.linewidth": 0.5, "grid.alpha": 0.3}
    with plt.rc_context(rc):
        fig, axes = plt.subplots(2, 3, figsize=(20, 14), constrained_layout=True)
        hb = None
        for i, a_conc in enumerate(amp_concs):
            for j, z_conc in enumerate(azt_concs):
                ax = axes[i, j]
                amp_fit = df.query("Drug=='AMP' & Concentration==@a_conc").Fitness
                azt_fit = df.query("Drug=='AZT' & Concentration==@z_conc").Fitness
                hb = ax.hist2d(amp_fit, azt_fit,
                               bins=[np.linspace(0, 4.5, 50), np.linspace(.95, 6.1, 50)],
                               cmap="BuPu", norm=LogNorm(), alpha=0.75)
                hb[3].set_rasterized(True)
                wt_a = amp_stats.query("mutant_profile==@WT_PROFILE and concentration==@a_conc")["median"].iat[0]
                dead_a = amp_stats.query("mutant_profile==@DEAD_PROFILE and concentration==@a_conc")["median"].iat[0]
                wt_z = azt_stats.query("mutant_profile==@WT_PROFILE and concentration==@z_conc")["median"].iat[0]
                dead_z = azt_stats.query("mutant_profile==@DEAD_PROFILE and concentration==@z_conc")["median"].iat[0]
                ax.plot(dead_a, dead_z, "o", color=colors["dead"], ms=12, mew=3, mfc="none", zorder=200)
                ax.plot(wt_a, wt_z, "o", color=colors["wt"], ms=12, mew=3, mfc="none", zorder=200)
                ex_a, ex_z = [], []
                for seq in extra.values():
                    a = df.query("Drug=='AMP' & Concentration==@a_conc & Genotype==@seq").Fitness
                    z = df.query("Drug=='AZT' & Concentration==@z_conc & Genotype==@seq").Fitness
                    if not a.empty and not z.empty:
                        ex_a.append(a.iat[0]); ex_z.append(z.iat[0])
                ax.plot(ex_a, ex_z, "D", color=colors["tested"], mec="k", ms=10, alpha=0.8)
                cl_a, cl_z = [], []
                for seq in clin:
                    a = df.query("Drug=='AMP' & Concentration==@a_conc & Genotype==@seq").Fitness
                    z = df.query("Drug=='AZT' & Concentration==@z_conc & Genotype==@seq").Fitness
                    if not a.empty and not z.empty:
                        cl_a.append(a.iat[0]); cl_z.append(z.iat[0])
                ax.plot(cl_a, cl_z, "^", color=colors["clinical"], mec="k", ms=11, alpha=0.8)
                ax.set_xlim(0, 4.6); ax.set_ylim(.95, 6.1)
                ax.set_xlabel(rf"Fitness in AMP {a_conc:.0f}$\mathrm{{\mu g/ml}}$")
                ax.set_ylabel(rf"Fitness in AZT {z_conc:.0f}$\mathrm{{\mu g/ml}}$")
                ax.grid(True)
                if j != 0:
                    ax.set_yticklabels([])
                if i != 0:
                    ax.set_xticklabels([])
                if i == 0 and j == 0:
                    ax.legend(handles=handles, loc="upper left")
        cbar_ax = fig.add_axes([1.01, 0.4, 0.03, 0.3], zorder=100)
        cbar = fig.colorbar(hb[3], cax=cbar_ax)
        cbar.ax.set_ylabel("Mutant Count", labelpad=0)
        fig.savefig(out, dpi=600, bbox_inches="tight")
        plt.close(fig)
    print("wrote", out.name)


def _logo(mask_source: str, title: str, out: Path) -> None:
    """Sequence logo of the top-1% genotypes in one drug (src/05 cells 21-23)."""
    if mask_source == "AMP":
        thr = np.percentile(AMP["Fitness"], 99); mask = AMP["Fitness"] >= thr
    else:
        thr = np.percentile(AZT["Fitness"], 99); mask = AZT["Fitness"] >= thr
    src = AMP if mask_source == "AMP" else AZT
    sequences = src[mask.values]["Genotype"].tolist()
    counts_df = logomaker.alignment_to_matrix(sequences, to_type="counts")
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    logo = logomaker.Logo(counts_df, ax=ax)
    logo.style_spines(visible=True)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    for ix, row in counts_df.iterrows():
        freq = row / row.sum()
        n_nonzero = (row > 0).sum()
        threshold = 0.55 if n_nonzero == 2 else 0.37 if n_nonzero == 3 else 0.275
        wt_aa = WT_SEQ[ix]
        for aa in freq.index:
            if aa == wt_aa:
                logo.style_single_glyph(p=ix, c=aa, color="black")
            elif freq[aa] > threshold:
                logo.style_single_glyph(p=ix, c=aa, color="red")
    ax.set_xticks(range(len(SIGN_POS)))
    ax.set_xticklabels(SIGN_POS, rotation=45, fontsize=14)
    ax.set_title(title, fontsize=16)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out.name)


def panel_CD(drug: str, out: Path) -> None:
    """Fitness vs number of mutated residues, jittered scatter + median (cells 11/12)."""
    data = AMP if drug == "AMP" else AZT
    color = "gray" if drug == "AMP" else "hotpink"
    alpha = 0.85 if drug == "AMP" else 0.95
    label = (r"Fitness in Ampicillin$_{781 \mu g/ml}$" if drug == "AMP"
             else r"Fitness in Aztreonam$_{36 \mu g/ml}$")
    np.random.seed(0)
    order = data["Epistatic Order"]
    jitter = np.where(order > 0, np.random.normal(0, 0.11, size=len(order)), 0)
    fig, ax = plt.subplots(figsize=(9, 2.7) if drug == "AMP" else (10, 3))
    sc = ax.scatter(order + jitter, data["Fitness"], s=3, alpha=alpha, color=color, label="Fitness")
    sc.set_rasterized(True)
    med = data.groupby("Epistatic Order")["Fitness"].median()
    for o, m in med.items():
        ax.plot(o, m, "ko", markersize=4 if drug == "AMP" else 3)
    ax.plot(med.index, med.values, "k-", alpha=1, linewidth=1, label="Median")
    ax.legend(loc="upper right", fontsize=10, frameon=True, edgecolor="black")
    ax.set_ylabel(label, fontsize=10)
    ax.set_xlabel("Number of Mutated Residues", fontsize=10)
    ax.set_xticks(range(14))
    ax.set_xticklabels(range(14), fontsize=10, rotation=45)
    ax.tick_params(axis="y", labelsize=10)
    plt.tight_layout()
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out.name)


if __name__ == "__main__":
    panel_A(OUT / "fig3_A.pdf")
    _logo("AMP", "Top 1% in Ampicillin, 781 μg/mL", OUT / "fig3_Bamp.pdf")
    _logo("AZT", "Top 1% in Aztreonam, 36 μg/mL", OUT / "fig3_Bazt.pdf")
    panel_CD("AMP", OUT / "fig3_C.pdf")
    panel_CD("AZT", OUT / "fig3_D.pdf")
