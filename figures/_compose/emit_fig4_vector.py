"""Emit Figure 4's eight panels as vector PDFs (selectable text, vector marks).

Panels C-H reproduce the exact plotting in ``src/05_epistasis_figures.ipynb``
(double-mutant heat-maps, measured-vs-predicted density grids, R2/rho-vs-order
curves), reading the same precomputed parquet, so the vector output reproduces
the manuscript raster panels. Panels A/B are single-mutant fitness distributions
(per-replicate strips) built from the replicate long-format table.

Text is embedded as TrueType (pdf.fonttype 42) so it stays editable. The
hist2d / heat-map cell meshes are light enough (<7k quads) to keep as vector.
Outputs go to ``figures/_compose/vector/panel_[A-H].pdf`` for
``compose_figure4_vector.py`` to assemble.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm
from scipy import stats

import figstyle

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def _repo_root(p: Path) -> Path:
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


HERE = Path(__file__).resolve().parent
ROOT = _repo_root(HERE)
OUT = HERE / "vector"
OUT.mkdir(exist_ok=True)

WT_SEQ = "LQMERMAGERTRN"
SIGN_POS_AMBLER = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
AZT_CONC, AMP_CONC = 36.0, 781.0

# ---- load data ----
df = pd.read_parquet(ROOT / "data/processed/Epistasis_Combined.parquet")
CMB = {
    "AMP": df[(df.Drug == "AMP") & (df.Concentration == AMP_CONC)].reset_index(drop=True),
    "AZT": df[(df.Drug == "AZT") & (df.Concentration == AZT_CONC)].reset_index(drop=True),
}
amp_stats_df = pd.read_parquet(ROOT / "data/processed/amp_auc_long_df.parquet")
azt_stats_df = pd.read_parquet(ROOT / "data/processed/azt_auc_long_df.parquet")
STATS = {"AMP": amp_stats_df, "AZT": azt_stats_df}
CONC = {"AMP": AMP_CONC, "AZT": AZT_CONC}


def _ref(stats_df, profile, conc, col):
    return stats_df[(stats_df["mutant_profile"] == profile)
                    & (stats_df["concentration"] == conc)][col].values[0]


WT_FIT = {d: _ref(STATS[d], "." * 13, CONC[d], "median") for d in ("AMP", "AZT")}
DEAD_FIT = {d: _ref(STATS[d], "X" * 13, CONC[d], "median") for d in ("AMP", "AZT")}
FIT_RANGE = {"AMP": (0.0, 4.75), "AZT": (1.25, 5.25)}


# ---- single-mutant tables ----
def _one_mut(drug: str) -> pd.DataFrame:
    sub = CMB[drug]
    m = sub[(sub["Epistatic Order"] == 1) | (sub["Epistatic Order"] == 0)].copy()
    m["Mutated Residues"] = m["Genotype"].apply(
        lambda x: SIGN_POS_AMBLER[next((i for i in range(len(x)) if x[i] != WT_SEQ[i]), -1)]
        if x != WT_SEQ else "WT")
    m["Mutated Amino Acids"] = m["Genotype"].apply(
        lambda x: next((x[i] for i in range(len(x)) if x[i] != WT_SEQ[i]), "WT")
        if x != WT_SEQ else "WT")
    m["Mutation"] = m.apply(
        lambda row: "WT" if row["Mutated Amino Acids"] == "WT"
        else f"{WT_SEQ[SIGN_POS_AMBLER.index(row['Mutated Residues'])]}"
             f"{row['Mutated Residues']}{row['Mutated Amino Acids']}", axis=1)
    return m


one_mut = {"AMP": _one_mut("AMP"), "AZT": _one_mut("AZT")}

# Mutation display order: WT, then positions ascending, then the designed
# amino-acid order within each position (src/05 cell 1 `intended` dict). This is
# deterministic (independent of parquet row order) and matches the single-mutant
# distribution panels; it also makes A/B and C/D read identically.
_MUT_AA_ORDER = [["P"], ["K"], ["L", "V"], ["K"], ["S", "H", "N"], ["T"], ["T"],
                 ["S"], ["K"], ["S", "C"], ["M"], ["L", "Q"], ["D"]]
MUTATIONS = ["WT"] + [f"{WT_SEQ[k]}{SIGN_POS_AMBLER[k]}{aa}"
                      for k in range(13) for aa in _MUT_AA_ORDER[k]]


# ---- double-mutant tables ----
def _zero_one_two(drug: str) -> pd.DataFrame:
    sub = CMB[drug]
    m = sub[sub["Epistatic Order"] <= 2].copy()
    m["Mutated Residues"] = m["Genotype"].apply(
        lambda x: [SIGN_POS_AMBLER[i] for i in range(len(x)) if x[i] != WT_SEQ[i]]
        if x != WT_SEQ else ["WT"])
    m["Mutated Amino Acids"] = m["Genotype"].apply(
        lambda x: [x[i] for i in range(len(x)) if x[i] != WT_SEQ[i]]
        if x != WT_SEQ else ["WT"])

    def _mstr(row):
        if row["Mutated Residues"] == ["WT"]:
            return "WT"
        return "; ".join(
            f"{WT_SEQ[SIGN_POS_AMBLER.index(pos)]}{pos}{aa}"
            for pos, aa in zip(row["Mutated Residues"], row["Mutated Amino Acids"]))

    m["Mutations"] = m.apply(_mstr, axis=1)
    return m


zero_one_two = {"AMP": _zero_one_two("AMP"), "AZT": _zero_one_two("AZT")}


# ============================ panels ============================
def _profile(genotype: str) -> str:
    return "".join("." if c == w else c for c, w in zip(genotype, WT_SEQ))


def panel_AB(drug: str, out: Path) -> None:
    """Single-mutant fitness distributions: 3 replicate dots + median bar."""
    figstyle.apply()
    dot = "0.42" if drug == "AMP" else "#c9327a"
    med = "#111111" if drug == "AMP" else "#c9327a"
    stats_df, conc = STATS[drug], CONC[drug]
    geno = dict(zip(one_mut[drug]["Mutation"], one_mut[drug]["Genotype"]))
    off = np.array([-0.16, 0.0, 0.16])

    fig, ax = plt.subplots(figsize=(6.4, 2.55))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="0.9", lw=0.8, zorder=0)
    for i, mut in enumerate(MUTATIONS):
        rows = stats_df[(stats_df["mutant_profile"] == _profile(geno[mut]))
                        & (stats_df["concentration"] == conc)]
        reps = rows[["replicate1", "replicate2", "replicate3"]].to_numpy().ravel()
        is_wt = mut == "WT"
        ax.scatter(i + off, reps, s=42, zorder=3,
                   facecolor="#8fe388" if is_wt else dot,
                   edgecolor="#2f7d2a" if is_wt else "0.15", linewidth=0.7,
                   alpha=1.0 if is_wt else 0.78)
        ax.plot([i - 0.28, i + 0.28], [np.median(reps)] * 2,
                color="#3a9c33" if is_wt else med, lw=2.4, zorder=4,
                solid_capstyle="butt")
    ax.axhline(WT_FIT[drug], ls=":", lw=1.3,
               color="0.45" if drug == "AMP" else "#c9327a", zorder=1)
    ax.set_xlim(-0.7, len(MUTATIONS) - 0.3)
    ax.set_ylim(0, None)
    ax.set_xticks(range(len(MUTATIONS)))
    ax.set_xticklabels(["WT" if m == "WT" else m for m in MUTATIONS],
                       rotation=45, ha="right", fontsize=11)
    ax.set_ylabel("Fitness", fontsize=17, labelpad=6)
    ax.tick_params(axis="y", labelsize=12)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", out.name)


def panel_heatmap(drug: str, out: Path) -> None:
    """Double-mutant fitness heat-map (matches Figure 4C/D)."""
    plt.rcParams["font.size"] = 8
    data = zero_one_two[drug]
    mutations = MUTATIONS
    size = len(mutations)
    matrix = np.full((size, size), np.nan)
    for _, row in data.iterrows():
        if row["Mutations"] == "WT":
            matrix[0, 0] = row["Fitness"]
        elif "; " not in row["Mutations"]:
            try:
                i, j = mutations.index("WT"), mutations.index(row["Mutations"])
                matrix[j, i] = matrix[i, j] = row["Fitness"]
            except ValueError:
                continue
        elif len(row["Mutations"].split("; ")) == 2:
            try:
                mut1, mut2 = row["Mutations"].split("; ")
                i, j = mutations.index(mut1), mutations.index(mut2)
                matrix[j, i] = matrix[i, j] = row["Fitness"]
            except ValueError:
                continue

    cmap = sns.color_palette("vlag", as_cmap=True)
    wt = WT_FIT[drug]
    plt.figure(figsize=(4, 4))
    heatmap = sns.heatmap(
        matrix, annot=False, cmap=cmap, center=wt, vmin=wt - 3, vmax=wt + 3,
        cbar_kws={"label": "Fitness", "use_gridspec": False, "location": "right",
                  "pad": 0.01, "fraction": 0.1, "shrink": 0.5, "anchor": (2.8, 0.9)})
    cbar = heatmap.collections[0].colorbar
    cbar.ax.set_ylabel("ΔAUC-Fitness", fontsize=10, labelpad=10)
    _deltas = np.array([-3, -2, -1, 0, 1, 2, 3])
    cbar.set_ticks(wt + _deltas)
    cbar.set_ticklabels([f"{d:+g}" if d != 0 else "WT" for d in _deltas])
    labels = ["none" if x == "WT" else x for x in mutations]
    plt.xticks(np.arange(size) + 0.5, labels, rotation=90, ha="center", va="bottom")
    ax = plt.gca()
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    plt.yticks(np.arange(size) + 0.5, labels, rotation=0)
    plt.xlabel("Mutation 2", fontsize=10, fontweight="bold", labelpad=-267)
    plt.ylabel("Mutation 1", fontsize=10, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print("wrote", out.name)


def panel_density(drug: str, out: Path) -> None:
    """Measured-vs-predicted density grid, orders 1-6 (matches Figure 4E/F)."""
    sns.set_style("ticks")
    sub = CMB[drug]
    bin_range = FIT_RANGE[drug]
    gridsize, max_order = 34, 6
    fig, axes = plt.subplots(1, max_order, figsize=(max_order * 2.25, 2.7))
    x_min, x_max = bin_range
    bins = np.linspace(bin_range[0], bin_range[1], gridsize)
    cmap = "Greys" if drug == "AMP" else "RdPu"
    for order in range(1, max_order + 1):
        ax = axes[order - 1]
        meas = sub["Fitness"]
        pred = sub[f"Fitness_predicted for order {order}"]
        slope, intercept, r_value, _, _ = stats.linregress(meas, pred)
        ax.hist2d(meas, pred, bins=[bins, bins], cmap=cmap, norm=LogNorm())
        ax.plot([x_min, x_max], [slope * x_min + intercept, slope * x_max + intercept],
                "b--", lw=1)
        ax.text(0.4, 0.9, f"Epistatic Order ≤ {order}", fontsize=10,
                fontweight="bold", ha="center", va="bottom", transform=ax.transAxes)
        ax.text(0.06, 0.80, f"$R^2 = {r_value ** 2:.2f}$", fontsize=10, color="blue",
                ha="left", va="bottom", transform=ax.transAxes)
        ax.set_xlabel("Fitness$_{measured}$", fontsize=12, labelpad=0, fontweight="bold")
        ax.set_xticks(np.arange(np.floor(x_min), np.ceil(x_max) + 0.5, 1))
        ax.set_yticks(np.arange(np.floor(x_min), np.ceil(x_max) + 0.5, 1))
        if order == 1:
            ax.set_ylabel("Fitness$_{predicted}$", fontsize=12, labelpad=0, fontweight="bold")
        else:
            ax.set_yticklabels([])
        ax.plot([x_min, x_max], [x_min, x_max], "k--", lw=1, alpha=0.95)
        ax.set_xlim(bin_range)
        ax.set_ylim(bin_range)
        ax.set_aspect("equal")
        ax.grid(False)
    plt.tight_layout(rect=[0, 0, 1, 1], h_pad=0.3, w_pad=0.12)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", out.name)


def panel_order(metric: str, out: Path) -> None:
    """R2 / rho vs highest-order epistatic terms (matches regen_gh.py)."""
    figstyle.apply()
    norders = 13
    orders = list(range(1, norders + 1))

    def curve(drug, conc):
        sub = df[(df.Drug == drug) & (df.Concentration == conc)]
        fit = sub["Fitness"].to_numpy()
        vals = []
        for o in range(1, norders + 1):
            pred = sub[f"Fitness_predicted for order {o}"].to_numpy()
            mask = np.isfinite(fit) & np.isfinite(pred)
            r = np.corrcoef(fit[mask], pred[mask])[0, 1]
            vals.append(r * r if metric == "r2" else r)
        return vals

    amp = curve("AMP", AMP_CONC)
    azt = curve("AZT", AZT_CONC)
    ylab = r"$\mathbf{R^2}$" if metric == "r2" else r"$\boldsymbol{\rho}$"
    fig, ax = plt.subplots(figsize=(4.05, 2.05))
    ax.grid(True, linestyle="--", linewidth=0.6, color="0.82", zorder=0)
    ax.plot(orders, amp, "o-", color="gray", lw=2.2, markersize=5.5,
            label="Ampicillin", zorder=3)
    ax.plot(orders, azt, "o-", color="hotpink", lw=2.2, markersize=5.5,
            label="Aztreonam", zorder=3)
    ax.set_xlim(0.4, norders + 0.6)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xticks(range(1, norders + 1))
    ax.set_xticklabels(range(1, norders + 1), rotation=45)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(labelsize=11, length=3)
    ax.set_xlabel("Highest-order Epistatic Terms Included", fontsize=12,
                  fontweight="bold", labelpad=4)
    ax.set_ylabel(ylab, fontsize=15, fontweight="bold", rotation=0, labelpad=20, va="center")
    ax.legend(loc="lower right", fontsize=11, frameon=True, edgecolor="0.35",
              fancybox=False, borderpad=0.5)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", out.name)


if __name__ == "__main__":
    panel_AB("AMP", OUT / "panel_A.pdf")
    panel_AB("AZT", OUT / "panel_B.pdf")
    panel_heatmap("AMP", OUT / "panel_C.pdf")
    panel_heatmap("AZT", OUT / "panel_D.pdf")
    panel_density("AMP", OUT / "panel_E.pdf")
    panel_density("AZT", OUT / "panel_F.pdf")
    panel_order("r2", OUT / "panel_G.pdf")
    panel_order("rho", OUT / "panel_H.pdf")
