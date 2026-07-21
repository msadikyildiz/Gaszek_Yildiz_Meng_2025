"""Supplementary Figure S3 - dose-response fitness profiles for a
representative subset of the TEM-1 CML library.

(A)/(B): raw genotype AUC ("auc") vs. antibiotic concentration for
Ampicillin and Aztreonam. Each panel shows the full library as a faint
gray background, the top-18 highest-AUC genotypes in color, and the spiked
wild-type (TEM-1^WT) and dead-control (TEM-1^dead) reference genotypes
highlighted.

Inputs
------
Reads the deposited genotype-AUC tables
    data/raw/{Ampicillin,Aztreonam}_auc_per_genotype.csv
one row per genotype, one column per "drug concentration replicate", sorted
descending by total AUC. These hold the same per-genotype AUC values written
by source_notebook/04-the-most-abundant-genotypes.ipynb (cell 17), which
derives them from raw barcode-count parquet files that are not part of the
public data deposit. The deposited CSVs and that intermediate have equal byte
counts (21,968,691 / 28,665,545 bytes); byte identity was not verified.

Method
------
Selects the top-18 genotypes by total AUC plus the WT and dead rows, plots a
1000-genotype random background sample, and highlights the references:
  1. TEM-1^WT / TEM-1^dead are drawn last in fixed black / blue on top of the
     tab20b qualitative palette used for the other 18 genotypes.
  2. The background sample (`np.random.choice`) is drawn with a fixed seed for
     reproducibility.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# --- paths -------------------------------------------------------------------
HERE = Path(__file__).resolve().parent


def _repo_root(p: Path) -> Path:
    """Walk up from this file to the repo root (dir with the manuscript parquet)."""
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


REPO = _repo_root(Path(__file__).resolve())
RAW = REPO / "data" / "raw"
FIGDIR = HERE / "figures"
SUPP = REPO / "figures" / "supplementary"
FIGDIR.mkdir(exist_ok=True, parents=True)
SUPP.mkdir(exist_ok=True, parents=True)

# --- constants -----------------------------------------------------------
SEED = 20240723  # fixed seed for the background-sample draw
N_TOP = 18
N_BACKGROUND = 1000

WT = "............."
DEAD = "XXXXXXXXXXXXX"
WT_LABEL = "TEM-1$^{WT}$"
DEAD_LABEL = "TEM-1$^{dead}$"
WT_COLOR = "black"
DEAD_COLOR = "#0072B2"  # matches the blue used for the dead reference elsewhere in this lab's figures

DRUG_FILES = {
    "Ampicillin": RAW / "Ampicillin_auc_per_genotype.csv",
    "Aztreonam": RAW / "Aztreonam_auc_per_genotype.csv",
}


def load_drug_table(path: Path) -> pd.DataFrame:
    """Load a genotype-AUC-per-concentration table.

    Already sorted descending by total AUC (this is exactly
    `df_sorted[drug]` from the source notebook's cell 17, written directly
    via `.to_csv()`, index preserved) -- so `df.index[:N_TOP]` reproduces
    the notebook's top-N selection without needing to re-derive the sort.
    """
    df = pd.read_csv(path, index_col=0)
    return df


def select_genotypes(df: pd.DataFrame, rng: np.random.Generator) -> tuple[list, list, list]:
    """Return (top_n_genotype_row_indices, wt_row_indices, dead_row_indices,
    background_row_indices) mirroring the source notebook's cell 18 exactly:
    `df.index[:18].to_list() + wt_index.to_list() + dead_index.to_list()`.
    """
    wt_idx = df.index[df["mut_profile_masked"] == WT].tolist()
    dead_idx = df.index[df["mut_profile_masked"] == DEAD].tolist()
    top_idx = df.index[:N_TOP].tolist()
    background_idx = rng.choice(df.index.to_numpy(), size=min(N_BACKGROUND, len(df)), replace=False)
    return top_idx, wt_idx, dead_idx, background_idx


def melt_long(df: pd.DataFrame, drug: str, row_indices) -> pd.DataFrame:
    """Wide (genotype x 'Drug conc rep' columns) -> long (conc, genotype, auc).

    Mirrors the source notebook's `slice_select` / cell-19 melt loop.
    """
    value_cols = [c for c in df.columns if c.startswith(drug + " ")]
    rows = []
    for col in value_cols:
        _, conc, _rep = col.split(" ")
        conc = float(conc)
        sub = df.loc[row_indices, ["mut_profile_masked", col]]
        for genotype, auc in zip(sub["mut_profile_masked"], sub[col]):
            rows.append({"conc": conc, "genotype": genotype, "auc": auc})
    return pd.DataFrame(rows)


def plot_panel(ax, long_background: pd.DataFrame, long_top: pd.DataFrame, drug: str) -> None:
    concs = sorted(long_background["conc"].unique())
    conc_labels = [f"{c:g}" for c in concs]

    # Background: full-library context, one very faint line per genotype.
    sns.pointplot(
        data=long_background, x="conc", y="auc", hue="genotype",
        order=concs, palette=["0.15"] * long_background["genotype"].nunique(),
        errorbar=None, alpha=0.05, linewidth=0.6, markersize=2,
        legend=False, ax=ax,
    )

    # Representative subset (excluding WT/dead, styled separately below).
    other = long_top[~long_top["genotype"].isin([WT, DEAD])]
    genotypes_in_rank_order = list(dict.fromkeys(other["genotype"]))  # already rank-ordered by AUC
    palette = sns.color_palette("tab20b", n_colors=max(len(genotypes_in_rank_order), 1))
    sns.pointplot(
        data=other, x="conc", y="auc", hue="genotype",
        order=concs, hue_order=genotypes_in_rank_order, palette=palette,
        errorbar="sd", alpha=0.9, linewidth=1.3, markersize=4, dodge=0.3,
        ax=ax,
    )

    # WT / dead reference genotypes, drawn last so they sit on top, in a
    # fixed color that does not depend on how many "other" genotypes there are.
    for genotype, color, label in ((WT, WT_COLOR, WT_LABEL), (DEAD, DEAD_COLOR, DEAD_LABEL)):
        sub = long_top[long_top["genotype"] == genotype]
        if sub.empty:
            continue
        stats = sub.groupby("conc")["auc"].agg(["mean", "std"]).reindex(concs)
        ax.errorbar(
            range(len(concs)), stats["mean"], yerr=stats["std"].fillna(0),
            marker="o", ms=6, lw=2.2, capsize=3, color=color, label=label, zorder=5,
        )

    ax.set_xticks(range(len(concs)))
    ax.set_xticklabels(conc_labels)
    ax.set_xlabel(f"{drug} concentration")
    ax.set_ylabel("auc")

    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, frameon=True)
    for text in leg.get_texts():
        if text.get_text() not in (WT_LABEL, DEAD_LABEL):
            text.set_fontfamily("monospace")


def main() -> None:
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(2, 1, figsize=(9, 12))
    panel_letters = ("A", "B")

    for ax, letter, (drug, path) in zip(axes, panel_letters, DRUG_FILES.items()):
        df = load_drug_table(path)
        top_idx, wt_idx, dead_idx, background_idx = select_genotypes(df, rng)
        selected_idx = top_idx + wt_idx + dead_idx

        long_background = melt_long(df, drug, background_idx)
        long_top = melt_long(df, drug, selected_idx)

        plot_panel(ax, long_background, long_top, drug)
        ax.text(-0.08, 1.03, f"({letter})", transform=ax.transAxes, fontsize=14, fontweight="bold")

        # Per-drug diagnostic copy (not the combined publication panel).
        fig_single, ax_single = plt.subplots(figsize=(7, 5))
        plot_panel(ax_single, long_background, long_top, drug)
        fig_single.tight_layout()
        fig_single.savefig(FIGDIR / f"dose_response_{drug.lower()}.png", dpi=150, bbox_inches="tight")
        plt.close(fig_single)

    fig.tight_layout()
    fig.savefig(SUPP / "figure_s03.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "figure_s03_combined.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {SUPP / 'figure_s03.png'}")


if __name__ == "__main__":
    main()
