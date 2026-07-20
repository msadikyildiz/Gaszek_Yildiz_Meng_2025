"""
S7 - Concentration-grid extension of Figure 4 (and the pairwise-epistasis panel
of Figure 5) for Reviewer #1 and Reviewer #3.

Produces three extended-figure grids per drug, each panelled over every
measured non-zero concentration:

  1. Pairwise epistasis heatmap   - one 19 x 19 panel per concentration.
                                    Uses the "Biochemical Definition" column
                                    that the manuscript's Figure 5E/F uses.
                                    (The task prompt calls this "Fig 4 C/D"
                                    by analogy.)

  2. Measured-vs-predicted density - one hist2d panel per (concentration,
                                    included epistatic order <= K) for
                                    K = 1, 2, ..., 6. Mirrors Fig 4 A/C.

  3. R^2 vs epistatic order       - one line per concentration, K in 1..13.
                                    Mirrors Fig 4 B/D but plots R^2 (per the
                                    S7 task brief) rather than RMSD.

All predictions come pre-computed from
  Gaszek_Yildiz_Meng_2025/data/processed/Epistasis_Combined.parquet
which stores `Fitness` (log10 AUC) alongside `Fitness_predicted for order K`
for K in 1..13 (see src/utils/calculate_epistasis.py). We do NOT re-fit the
partial linear-regression epistasis model - we match the manuscript exactly.

Concentrations analysed (locked in the S7 brief; zero-drug panels excluded):

    AMP: 3.1, 12.2, 48.8, 195.0, 781.0   microgram/mL    (5 panels per figure)
    AZT: 0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0        (7 panels per figure)

Colour scales are per-drug shared across panels so a panel-to-panel
comparison is legitimate. AMP and AZT get independent scales.
"""

from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm, LinearSegmentedColormap

# --- paths -------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / "data" / "processed" / "Epistasis_Combined.parquet"
FIGDIR = HERE / "figures"
DATADIR = HERE / "data"
SUPP = REPO / "figures" / "supplementary"
FIGDIR.mkdir(exist_ok=True, parents=True)
DATADIR.mkdir(exist_ok=True, parents=True)
SUPP.mkdir(exist_ok=True, parents=True)

# --- constants ---------------------------------------------------------------
AMBLER_POS = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
WT_LETTER = "LQMERMAGERTRN"
ALT_LETTERS = [
    ["P"], ["K"], ["L", "V"], ["K"], ["H", "N", "S"],
    ["T"], ["T"], ["S"], ["K"], ["C", "S"],
    ["M"], ["L", "Q"], ["D"],
]

AMP_CONCS = [3.1, 12.2, 48.8, 195.0, 781.0]
AZT_CONCS = [0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]
DRUG_CONCS = {"AMP": AMP_CONCS, "AZT": AZT_CONCS}
DRUG_NAME = {"AMP": "Ampicillin", "AZT": "Aztreonam"}

# Orders shown in the pred-vs-measured grid (cumulative "<= K" features).
PRED_ORDERS = [1, 2, 3, 4, 5, 6]
# Orders evaluated in the R^2 vs order plot (full range).
R2_ORDERS = list(range(1, 14))

# Per-drug colormaps and scalar limits for pairwise heatmap (Biochemical Def).
# The manuscript chose vmin/vmax = -3/+3 for AMP and -4/+4 for AZT in Fig 5E/F
# at the single "featured" concentration. Across the full concentration grid,
# the strongest pairwise epistasis values appear at the highest concentrations
# (sign-epistasis regime); we keep the manuscript's symmetric limits so that
# each panel is directly comparable to the main-text figure.
PAIRWISE_CMAP = "RdBu_r"
PAIRWISE_VLIM = {"AMP": 3.0, "AZT": 4.0}

# Per-drug colormap for the 2D density panels (matches manuscript).
DENSITY_CMAP = {"AMP": "Greys", "AZT": "RdPu"}

# R^2 vs order line-plot colours per concentration - one warm ramp per drug.
AMP_RAMP = plt.cm.Greys(np.linspace(0.35, 0.95, len(AMP_CONCS)))
AZT_RAMP = plt.cm.RdPu(np.linspace(0.35, 0.95, len(AZT_CONCS)))
CONC_COLOR = {
    ("AMP", c): AMP_RAMP[i] for i, c in enumerate(AMP_CONCS)
}
CONC_COLOR.update(
    {("AZT", c): AZT_RAMP[i] for i, c in enumerate(AZT_CONCS)}
)

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"


# --- helpers -----------------------------------------------------------------
def mutation_labels() -> list[str]:
    """Single-mutation strings 'L21P', 'Q39K', ... in position order."""
    out: list[str] = []
    for pos, wt, alts in zip(AMBLER_POS, WT_LETTER, ALT_LETTERS):
        for a in alts:
            out.append(f"{wt}{pos}{a}")
    return out


def genotype_to_mut_list(gt: str) -> list[str]:
    """Turn a 13-char genotype into ['WT'] or the list of point mutations."""
    muts: list[str] = []
    for i, (wt, obs) in enumerate(zip(WT_LETTER, gt)):
        if obs != wt:
            muts.append(f"{wt}{AMBLER_POS[i]}{obs}")
    return muts if muts else ["WT"]


def calc_r2_rmsd(y: np.ndarray, yhat: np.ndarray) -> tuple[float, float]:
    """R^2 (Pearson squared) and RMSD, ignoring non-finite pairs."""
    ok = np.isfinite(y) & np.isfinite(yhat)
    if ok.sum() < 2:
        return float("nan"), float("nan")
    y_ok, yhat_ok = y[ok], yhat[ok]
    r = np.corrcoef(y_ok, yhat_ok)[0, 1]
    rmsd = float(np.sqrt(np.mean((y_ok - yhat_ok) ** 2)))
    return float(r * r), rmsd


# --- data loading ------------------------------------------------------------
def load_drug_conc(df_all: pl.DataFrame, drug: str, conc: float) -> pl.DataFrame:
    return df_all.filter(
        (pl.col("Drug") == drug) & (pl.col("Concentration") == conc)
    )


# --- pairwise heatmap --------------------------------------------------------
MUTATIONS = ["WT"] + mutation_labels()  # 19 entries: WT + 18 point mutations


def pairwise_matrix(df_conc: pl.DataFrame) -> np.ndarray:
    """
    Build a 19x19 symmetric matrix of pairwise epistasis (Biochemical
    Definition). Off-diagonal (i,j) with i != j holds the pairwise-epistasis
    value for the double mutant {muts[i], muts[j]}. Row/column 0 ("WT")
    stores the order-1 biochemical value on the corresponding mutation row
    (single-mutant effect in WT background).
    """
    size = len(MUTATIONS)
    M = np.full((size, size), np.nan, dtype=np.float64)

    subset = df_conc.filter(pl.col("Epistatic Order").is_in([1, 2]))
    geno = subset["Genotype"].to_list()
    bioch = subset["Biochemical Definition"].to_numpy()

    for g, v in zip(geno, bioch):
        muts = genotype_to_mut_list(g)
        if len(muts) == 1:
            j = MUTATIONS.index(muts[0])
            M[0, j] = v
            M[j, 0] = v
        elif len(muts) == 2:
            i = MUTATIONS.index(muts[0])
            j = MUTATIONS.index(muts[1])
            M[i, j] = v
            M[j, i] = v
    return M


def plot_pairwise_grid(df_all: pl.DataFrame, drug: str, out_path: Path) -> pl.DataFrame:
    concs = DRUG_CONCS[drug]
    n = len(concs)
    ncols = n
    fig_w = 3.1 * ncols + 1.0
    fig, axes = plt.subplots(1, ncols, figsize=(fig_w, 3.6), dpi=150)
    axes = np.atleast_1d(axes)

    vlim = PAIRWISE_VLIM[drug]
    long_rows: list[dict] = []

    for col_i, conc in enumerate(concs):
        ax = axes[col_i]
        sub = load_drug_conc(df_all, drug, conc)
        M = pairwise_matrix(sub)

        im = ax.imshow(
            M,
            cmap=PAIRWISE_CMAP,
            vmin=-vlim,
            vmax=vlim,
            origin="upper",
        )
        title = f"{drug} {conc:g} µg/mL"
        ax.set_title(title, fontsize=9.5, pad=4)
        ax.set_xticks(range(len(MUTATIONS)))
        ax.set_yticks(range(len(MUTATIONS)))
        xlabels = ["none" if m == "WT" else m for m in MUTATIONS]
        ax.set_xticklabels(xlabels, rotation=90, fontsize=5.4)
        if col_i == 0:
            ax.set_yticklabels(xlabels, fontsize=5.4)
            ax.set_ylabel("Mutation 1", fontsize=9)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("Mutation 2", fontsize=9)
        ax.tick_params(axis="both", length=1.5, pad=1)

        # Dump long-form data for CSV output
        for i in range(len(MUTATIONS)):
            for j in range(i, len(MUTATIONS)):
                v = M[i, j]
                if np.isfinite(v):
                    long_rows.append({
                        "drug": drug,
                        "concentration": conc,
                        "mut_i": MUTATIONS[i],
                        "mut_j": MUTATIONS[j],
                        "biochemical_epistasis": float(v),
                    })

    cbar = fig.colorbar(im, ax=axes, shrink=0.78, pad=0.01, aspect=24)
    cbar.set_label("Pairwise epistasis\n(biochemical)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    fig.suptitle(
        f"{DRUG_NAME[drug]} - pairwise epistasis across all measured concentrations",
        fontsize=11, y=1.02,
    )
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")

    return pl.DataFrame(long_rows)


# --- measured-vs-predicted density grid --------------------------------------
def plot_pred_vs_measured(df_all: pl.DataFrame, drug: str, out_path: Path) -> None:
    concs = DRUG_CONCS[drug]
    nrows = len(concs)
    ncols = len(PRED_ORDERS)

    # Determine per-drug shared axis bounds (all concs + all orders + measured)
    # Use small padding so all density points have space.
    sub_all = df_all.filter(
        (pl.col("Drug") == drug) & pl.col("Concentration").is_in(concs)
    )
    measured = sub_all["Fitness"].to_numpy()
    pred_cols = [f"Fitness_predicted for order {k}" for k in PRED_ORDERS]
    all_vals = [measured]
    for c in pred_cols:
        all_vals.append(sub_all[c].to_numpy())
    all_arr = np.concatenate(all_vals)
    all_arr = all_arr[np.isfinite(all_arr)]
    lo = float(np.quantile(all_arr, 0.0005))
    hi = float(np.quantile(all_arr, 0.9995))
    pad = 0.12 * (hi - lo)
    bin_range = (lo - pad, hi + pad)
    bins = np.linspace(*bin_range, 34)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(1.85 * ncols, 1.85 * nrows),
        dpi=150,
    )
    if nrows == 1:
        axes = axes[None, :]

    for r, conc in enumerate(concs):
        sub = load_drug_conc(df_all, drug, conc)
        y = sub["Fitness"].to_numpy()
        for c_i, k in enumerate(PRED_ORDERS):
            ax = axes[r, c_i]
            yhat = sub[f"Fitness_predicted for order {k}"].to_numpy()
            ok = np.isfinite(y) & np.isfinite(yhat)
            ax.hist2d(
                y[ok], yhat[ok],
                bins=[bins, bins],
                cmap=DENSITY_CMAP[drug],
                norm=LogNorm(),
            )
            r2, rmsd = calc_r2_rmsd(y, yhat)
            ax.plot(
                [bin_range[0], bin_range[1]],
                [bin_range[0], bin_range[1]],
                "k--", lw=0.6, alpha=0.8,
            )
            ax.text(
                0.04, 0.96,
                f"R\u00b2 = {r2:.2f}",
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=8, color="#1f3f7a", fontweight="bold",
            )
            ax.set_xlim(*bin_range)
            ax.set_ylim(*bin_range)
            ax.set_aspect("equal", adjustable="box")
            ax.tick_params(axis="both", labelsize=7, length=2)

            if r == 0:
                ax.set_title(f"order \u2264 {k}", fontsize=9)
            if r == nrows - 1:
                ax.set_xlabel("measured", fontsize=8)
            else:
                ax.set_xticklabels([])
            if c_i == 0:
                ax.set_ylabel(
                    f"{conc:g} µg/mL\npredicted",
                    fontsize=8.5,
                )
            else:
                ax.set_yticklabels([])

    fig.suptitle(
        f"{DRUG_NAME[drug]} - measured vs. predicted fitness "
        f"across all concentrations and epistatic orders",
        fontsize=11, y=1.005,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.995))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- R^2 vs order ------------------------------------------------------------
def compute_r2_table(df_all: pl.DataFrame) -> pl.DataFrame:
    rows: list[dict] = []
    for drug, concs in DRUG_CONCS.items():
        for conc in concs:
            sub = load_drug_conc(df_all, drug, conc)
            y = sub["Fitness"].to_numpy()
            for k in R2_ORDERS:
                yhat = sub[f"Fitness_predicted for order {k}"].to_numpy()
                r2, rmsd = calc_r2_rmsd(y, yhat)
                rows.append({
                    "drug": drug,
                    "concentration": conc,
                    "order": k,
                    "r2": r2,
                    "rmsd": rmsd,
                    "n": int(np.isfinite(y).sum()),
                })
    return pl.DataFrame(rows)


def plot_r2_vs_order(r2_df: pl.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.3), dpi=150, sharey=True)

    for ax_i, drug in enumerate(["AMP", "AZT"]):
        ax = axes[ax_i]
        concs = DRUG_CONCS[drug]
        for conc in concs:
            sub = (
                r2_df
                .filter((pl.col("drug") == drug) & (pl.col("concentration") == conc))
                .sort("order")
            )
            x = sub["order"].to_numpy()
            y = sub["r2"].to_numpy()
            ax.plot(
                x, y,
                marker="o", markersize=4, linewidth=1.7,
                color=CONC_COLOR[(drug, conc)],
                label=f"{conc:g} µg/mL",
            )
        ax.set_xlabel("highest-order epistatic terms included", fontsize=10)
        if ax_i == 0:
            ax.set_ylabel("R²  (measured vs. predicted fitness)", fontsize=10)
        ax.set_xticks(R2_ORDERS)
        ax.set_xlim(0.7, 13.3)
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, alpha=0.25, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title(DRUG_NAME[drug], fontsize=11, fontweight="bold")
        leg = ax.legend(
            loc="lower right", fontsize=7.2, frameon=False,
            title=f"{DRUG_NAME[drug]} concentration",
            title_fontsize=8, handlelength=1.6,
        )
        leg._legend_box.align = "left"

    fig.suptitle(
        "R² of partial-linear-regression fit versus maximum epistatic order, per concentration",
        fontsize=11, y=1.03,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- main --------------------------------------------------------------------
def build_results_table(r2_df: pl.DataFrame) -> pl.DataFrame:
    """Pivot R^2 table into the headline drug x conc x order wide table."""
    wide = (
        r2_df
        .pivot(on="order", index=["drug", "concentration"], values="r2")
        .sort(["drug", "concentration"])
    )
    wide.columns = [
        f"r2_order_{c}" if c.isdigit() else c
        for c in wide.columns
    ]
    return wide


def main() -> None:
    print(f"reading {PARQUET}")
    df_all = pl.read_parquet(PARQUET)
    print(f"  {df_all.height} rows")

    # (1) Pairwise heatmap grids
    pair_frames: list[pl.DataFrame] = []
    for drug in ["AMP", "AZT"]:
        out = SUPP / f"figure_s11_{drug.lower()}.png"
        pair_frames.append(plot_pairwise_grid(df_all, drug, out))
    pair_long = pl.concat(pair_frames)
    pair_long.write_csv(DATADIR / "pairwise_mean_fitness.csv")
    print(f"wrote {DATADIR / 'pairwise_mean_fitness.csv'} ({pair_long.height} rows)")

    # (2) Measured-vs-predicted density grids
    for drug in ["AMP", "AZT"]:
        out = SUPP / f"figure_s12_{drug.lower()}.png"
        plot_pred_vs_measured(df_all, drug, out)

    # (3) R^2 vs order
    r2_long = compute_r2_table(df_all)
    r2_long.write_csv(DATADIR / "regression_r2_by_order.csv")
    print(f"wrote {DATADIR / 'regression_r2_by_order.csv'} ({r2_long.height} rows)")

    out = FIGDIR / "r2_vs_order.png"
    plot_r2_vs_order(r2_long, out)

    # Results table: compact R^2 by drug x concentration x order.
    wide = build_results_table(r2_long)
    wide.write_csv(HERE / "results_table.csv")
    print(f"wrote {HERE / 'results_table.csv'} ({wide.height} rows)")


if __name__ == "__main__":
    main()
