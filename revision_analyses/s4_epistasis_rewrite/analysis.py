"""
S4 - Epistasis section rewrite and two new extended figures.

Produces the analytical backing for the rewritten Epistasis section
(Reviewer 3 requested a global-to-local narration) plus the two new
extended figures called out in the draft:

    Ext Fig A - Order-decomposition: cumulative R^2 and incremental
                 dR^2 per additional epistatic order, one curve per
                 concentration, stratified by drug.

    Ext Fig E - Drug asymmetry at matched selective stringency:
                 R^2-vs-order, RMSD-vs-order, measured-fitness
                 dynamic-range histogram, pairwise |Biochemical-
                 epistasis| distribution. Primary pair is the main-
                 text concentrations (AMP 781 ug/mL, AZT 36 ug/mL),
                 with a secondary "low-stringency" pair (AMP 195
                 vs AZT 12) for robustness.

All numbers come from the manuscript's pre-computed partial linear-
regression predictions in
  Gaszek_Yildiz_Meng_2025/data/processed/Epistasis_Combined.parquet
(columns: Fitness, Fitness_predicted for order K, Biochemical
 Definition). We do NOT refit any model here.

Reuses derivations already computed in
  revision_analyses/s7_concentration_grid/data/regression_r2_by_order.csv
as a cross-check for R^2.

Outputs:
  data/order_decomposition.csv           - (drug, conc, order, R^2, dR^2, cum_R^2)
  data/matched_stringency_summary.csv    - (drug, conc_label, order, R^2, RMSD, n)
  data/source_data.xlsx                  - Excel, one sheet per figure
  figures/ext_order_decomposition_amp.png
  figures/ext_order_decomposition_azt.png
  figures/ext_order_decomposition_combined.png
  figures/ext_drug_asymmetry.png
  results_table.csv                      - Compact summary
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

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
S7_R2 = HERE.parent / "s7_concentration_grid" / "data" / "regression_r2_by_order.csv"
S7_PAIRWISE = HERE.parent / "s7_concentration_grid" / "data" / "pairwise_mean_fitness.csv"
FIGDIR = HERE / "figures"
DATADIR = HERE / "data"
FIGDIR.mkdir(exist_ok=True, parents=True)
DATADIR.mkdir(exist_ok=True, parents=True)

SEED = 20260420

# --- constants ---------------------------------------------------------------
AMP_CONCS = [3.1, 12.2, 48.8, 195.0, 781.0]
AZT_CONCS = [0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]
DRUG_CONCS = {"AMP": AMP_CONCS, "AZT": AZT_CONCS}
DRUG_NAME = {"AMP": "Ampicillin", "AZT": "Aztreonam"}
R2_ORDERS = list(range(1, 14))

# Primary matched-stringency pair (main-text concentrations) and secondary pair.
MATCHED_PAIRS = {
    "primary": {"AMP": 781.0, "AZT": 36.0, "label": "main-text (just-above-WT-MIC)"},
    "secondary": {"AMP": 195.0, "AZT": 12.0, "label": "alternative (~3x WT-MIC)"},
}

# Drug-colour convention (matches manuscript).
AMP_PRIMARY = "#2f2f2f"   # dark grey
AMP_SECONDARY = "#7a7a7a"  # mid grey
AZT_PRIMARY = "#c2185b"   # hotpink / RdPu darkest
AZT_SECONDARY = "#e898b9"  # RdPu mid
DRUG_PRIMARY = {"AMP": AMP_PRIMARY, "AZT": AZT_PRIMARY}
DRUG_SECONDARY = {"AMP": AMP_SECONDARY, "AZT": AZT_SECONDARY}

# Concentration ramps.
AMP_RAMP = plt.cm.Greys(np.linspace(0.35, 0.92, len(AMP_CONCS)))
AZT_RAMP = plt.cm.RdPu(np.linspace(0.35, 0.92, len(AZT_CONCS)))
CONC_COLOR = {("AMP", c): AMP_RAMP[i] for i, c in enumerate(AMP_CONCS)}
CONC_COLOR.update({("AZT", c): AZT_RAMP[i] for i, c in enumerate(AZT_CONCS)})

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"


# --- helpers -----------------------------------------------------------------
def calc_r2_rmsd(y: np.ndarray, yhat: np.ndarray) -> tuple[float, float]:
    """R^2 (Pearson squared) and RMSD, ignoring non-finite pairs."""
    ok = np.isfinite(y) & np.isfinite(yhat)
    if ok.sum() < 2:
        return float("nan"), float("nan")
    y_ok, yhat_ok = y[ok], yhat[ok]
    r = np.corrcoef(y_ok, yhat_ok)[0, 1]
    rmsd = float(np.sqrt(np.mean((y_ok - yhat_ok) ** 2)))
    return float(r * r), rmsd


def load_drug_conc(df_all: pl.DataFrame, drug: str, conc: float) -> pl.DataFrame:
    return df_all.filter(
        (pl.col("Drug") == drug) & (pl.col("Concentration") == conc)
    )


# --- Step 1: order-decomposition (cumulative + incremental dR^2) --------------
def compute_order_decomposition(df_all: pl.DataFrame) -> pl.DataFrame:
    """Return (drug, conc, order, r2, dr2, cum_r2, rmsd, n).

    Convention: R^2(0) = 0, so dR^2(K=1) = R^2(1).
    cum_r2 = R^2(K) = Pearson correlation squared at that order.
    dr2    = R^2(K) - R^2(K-1) (with R^2(0)=0).
    """
    rows: list[dict] = []
    for drug, concs in DRUG_CONCS.items():
        for conc in concs:
            sub = load_drug_conc(df_all, drug, conc)
            y = sub["Fitness"].to_numpy()
            prev_r2 = 0.0
            for k in R2_ORDERS:
                yhat = sub[f"Fitness_predicted for order {k}"].to_numpy()
                r2, rmsd = calc_r2_rmsd(y, yhat)
                rows.append({
                    "drug": drug,
                    "concentration": conc,
                    "order": k,
                    "r2": r2,
                    "cum_r2": r2,
                    "dr2": r2 - prev_r2,
                    "rmsd": rmsd,
                    "n": int(np.isfinite(y).sum()),
                })
                prev_r2 = r2
    return pl.DataFrame(rows)


# --- Extended Fig A ----------------------------------------------------------
def _plot_order_decomposition_panel(
    ax_line,
    ax_bar,
    df_drug: pl.DataFrame,
    drug: str,
) -> None:
    """Plot cumulative R^2 as lines (per conc) and ΔR^2 grouped bars.

    ax_line: matplotlib axis for cumulative-R^2 lines.
    ax_bar:  matplotlib axis for ΔR^2 grouped bars (shares x-axis with ax_line).
    """
    concs = DRUG_CONCS[drug]
    x = np.array(R2_ORDERS, dtype=float)

    # --- line plot: cumulative R^2 per concentration ------------------------
    for conc in concs:
        sub = (
            df_drug
            .filter(pl.col("concentration") == conc)
            .sort("order")
        )
        cum = sub["cum_r2"].to_numpy()
        ax_line.plot(
            x, cum,
            marker="o", markersize=4, linewidth=1.6,
            color=CONC_COLOR[(drug, conc)],
            label=f"{conc:g} µg/mL",
        )

    # Annotate R^2 = 0.95 and 0.99 gridlines (inline labels, above the line,
    # spaced along x so they don't stack where curves asymptote near K=10-13).
    for thr, ls, x_pos in [(0.95, ":", 3.5), (0.99, "--", 7.0)]:
        ax_line.axhline(thr, color="k", linewidth=0.6, linestyle=ls, alpha=0.4)
        ax_line.text(
            x_pos, thr + 0.008, f"R² = {thr:g}",
            ha="center", va="bottom", fontsize=7, color="0.3",
        )

    # Vertical guides at K = 1, 2, 3, 6.
    for k in (1, 2, 3, 6):
        ax_line.axvline(k, color="k", linewidth=0.5, linestyle=":", alpha=0.25)

    ax_line.set_xlim(0.7, 13.3)
    ax_line.set_ylim(-0.02, 1.05)
    ax_line.set_xticks(R2_ORDERS)
    ax_line.set_xticklabels([])
    ax_line.set_ylabel("cumulative R²", fontsize=10)
    ax_line.set_title(
        f"{DRUG_NAME[drug]} - cumulative R² and incremental ΔR² per epistatic order",
        fontsize=11, fontweight="bold", loc="left", pad=6,
    )
    ax_line.grid(True, axis="y", alpha=0.2, linewidth=0.5)
    ax_line.spines["top"].set_visible(False)
    ax_line.spines["right"].set_visible(False)

    leg = ax_line.legend(
        loc="lower right", fontsize=7.4, frameon=False,
        title=f"{drug} concentration",
        title_fontsize=8.0, handlelength=1.6,
        ncol=1, borderpad=0.2,
    )
    leg._legend_box.align = "left"

    # --- bar plot: incremental ΔR^2 grouped per order -----------------------
    n_concs = len(concs)
    bar_w = 0.85 / n_concs
    for i, conc in enumerate(concs):
        sub = (
            df_drug
            .filter(pl.col("concentration") == conc)
            .sort("order")
        )
        dr2 = sub["dr2"].to_numpy()
        offset = (i - (n_concs - 1) / 2.0) * bar_w
        ax_bar.bar(
            x + offset, dr2,
            width=bar_w,
            color=CONC_COLOR[(drug, conc)],
            edgecolor="none",
        )

    ax_bar.set_xlim(0.7, 13.3)
    ax_bar.set_ylim(0, max(
        0.35,
        float(df_drug["dr2"].max()) * 1.1,
    ))
    ax_bar.set_xticks(R2_ORDERS)
    ax_bar.set_xlabel("epistatic order K (highest term included)", fontsize=10)
    ax_bar.set_ylabel("incremental ΔR²", fontsize=9.5)
    ax_bar.grid(True, axis="y", alpha=0.2, linewidth=0.5)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    # Guides at K = 1, 2, 3, 6 to match line plot.
    for k in (1, 2, 3, 6):
        ax_bar.axvline(k, color="k", linewidth=0.5, linestyle=":", alpha=0.25)


def plot_ext_fig_a_single(df_order: pl.DataFrame, drug: str, out_path: Path) -> None:
    """Single-drug version of Ext Fig A (stand-alone PNG)."""
    fig = plt.figure(figsize=(7.2, 5.6), dpi=150)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.8, 1.0], hspace=0.09)
    ax_line = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[1, 0])

    df_drug = df_order.filter(pl.col("drug") == drug)
    _plot_order_decomposition_panel(ax_line, ax_bar, df_drug, drug)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


def plot_ext_fig_a_combined(df_order: pl.DataFrame, out_path: Path) -> None:
    """Two-drug stacked version of Ext Fig A (AMP top, AZT bottom)."""
    fig = plt.figure(figsize=(8.6, 10.2), dpi=150)
    # Two drug blocks with a large gap between them.
    gs_outer = fig.add_gridspec(2, 1, hspace=0.42)
    gs_amp = gs_outer[0].subgridspec(2, 1, height_ratios=[1.8, 1.0], hspace=0.09)
    gs_azt = gs_outer[1].subgridspec(2, 1, height_ratios=[1.8, 1.0], hspace=0.09)
    axes = {
        "AMP_line": fig.add_subplot(gs_amp[0, 0]),
        "AMP_bar":  fig.add_subplot(gs_amp[1, 0]),
        "AZT_line": fig.add_subplot(gs_azt[0, 0]),
        "AZT_bar":  fig.add_subplot(gs_azt[1, 0]),
    }
    for drug in ("AMP", "AZT"):
        df_drug = df_order.filter(pl.col("drug") == drug)
        _plot_order_decomposition_panel(
            axes[f"{drug}_line"], axes[f"{drug}_bar"], df_drug, drug,
        )

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- Step 2: Matched-stringency comparison (Ext Fig E) -----------------------
def compute_matched_stringency(
    df_all: pl.DataFrame, df_order: pl.DataFrame,
) -> tuple[pl.DataFrame, dict]:
    """Return long-form R^2/RMSD table + metadata for the two matched pairs."""
    rows: list[dict] = []
    extras: dict = {}

    for tier, spec in MATCHED_PAIRS.items():
        for drug in ("AMP", "AZT"):
            conc = spec[drug]
            sub_ord = (
                df_order
                .filter((pl.col("drug") == drug) & (pl.col("concentration") == conc))
                .sort("order")
            )
            for row in sub_ord.iter_rows(named=True):
                rows.append({
                    "tier": tier,
                    "tier_label": spec["label"],
                    "drug": drug,
                    "concentration": conc,
                    "order": row["order"],
                    "r2": row["r2"],
                    "rmsd": row["rmsd"],
                    "n": row["n"],
                })

            # Also snapshot fitness distribution (primary tier only).
            if tier == "primary":
                sub = load_drug_conc(df_all, drug, conc)
                y = sub["Fitness"].to_numpy()
                y = y[np.isfinite(y)]
                extras[f"fitness_{drug}_{conc}"] = y

                # |pairwise epistasis| distribution from S7 CSV.
                pair_df = pl.read_csv(S7_PAIRWISE)
                pair_df_drug = pair_df.filter(
                    (pl.col("drug") == drug) & (pl.col("concentration") == conc)
                )
                vals = pair_df_drug["biochemical_epistasis"].to_numpy()
                vals = vals[np.isfinite(vals)]
                extras[f"pairabs_{drug}_{conc}"] = np.abs(vals)

    return pl.DataFrame(rows), extras


def plot_ext_fig_e(
    matched_long: pl.DataFrame,
    extras: dict,
    df_order: pl.DataFrame,
    out_path: Path,
) -> None:
    """2x2 drug-asymmetry panel figure for matched-stringency comparison."""
    fig = plt.figure(figsize=(9.0, 7.4), dpi=150)
    gs = fig.add_gridspec(2, 2, hspace=0.33, wspace=0.28)
    ax_r2 = fig.add_subplot(gs[0, 0])
    ax_rmsd = fig.add_subplot(gs[0, 1])
    ax_fit = fig.add_subplot(gs[1, 0])
    ax_pair = fig.add_subplot(gs[1, 1])

    primary_pair = MATCHED_PAIRS["primary"]
    secondary_pair = MATCHED_PAIRS["secondary"]

    # Panel (a): R^2 vs order - primary pair solid, secondary pair dashed.
    for spec, ls, alpha in [(primary_pair, "-", 1.0), (secondary_pair, "--", 0.75)]:
        for drug in ("AMP", "AZT"):
            conc = spec[drug]
            sub = (
                matched_long
                .filter(
                    (pl.col("drug") == drug) & (pl.col("concentration") == conc)
                )
                .sort("order")
            )
            x = sub["order"].to_numpy()
            y = sub["r2"].to_numpy()
            ax_r2.plot(
                x, y, linestyle=ls, marker="o", markersize=3.2,
                linewidth=1.7, alpha=alpha,
                color=DRUG_PRIMARY[drug] if spec is primary_pair else DRUG_SECONDARY[drug],
                label=f"{drug} {conc:g} µg/mL" if spec is primary_pair else None,
            )
    ax_r2.set_xticks(R2_ORDERS)
    ax_r2.set_xlim(0.7, 13.3)
    ax_r2.set_ylim(-0.02, 1.02)
    ax_r2.set_xlabel("epistatic order K", fontsize=9.5)
    ax_r2.set_ylabel("R²", fontsize=10)
    ax_r2.set_title("(a) R² vs. epistatic order", fontsize=10.5, fontweight="bold", loc="left")
    ax_r2.grid(True, alpha=0.2, linewidth=0.5)
    ax_r2.spines["top"].set_visible(False)
    ax_r2.spines["right"].set_visible(False)
    ax_r2.legend(loc="lower right", fontsize=8, frameon=False, title="matched stringency", title_fontsize=8.2)

    # Panel (b): RMSD vs order.
    for spec, ls, alpha in [(primary_pair, "-", 1.0), (secondary_pair, "--", 0.75)]:
        for drug in ("AMP", "AZT"):
            conc = spec[drug]
            sub = (
                matched_long
                .filter(
                    (pl.col("drug") == drug) & (pl.col("concentration") == conc)
                )
                .sort("order")
            )
            ax_rmsd.plot(
                sub["order"].to_numpy(), sub["rmsd"].to_numpy(),
                linestyle=ls, marker="o", markersize=3.2, linewidth=1.7, alpha=alpha,
                color=DRUG_PRIMARY[drug] if spec is primary_pair else DRUG_SECONDARY[drug],
                label=f"{drug} {conc:g} µg/mL" if spec is primary_pair else None,
            )
    ax_rmsd.set_xticks(R2_ORDERS)
    ax_rmsd.set_xlim(0.7, 13.3)
    ax_rmsd.set_xlabel("epistatic order K", fontsize=9.5)
    ax_rmsd.set_ylabel("RMSD (log₁₀ AUC)", fontsize=10)
    ax_rmsd.set_title("(b) RMSD vs. epistatic order", fontsize=10.5, fontweight="bold", loc="left")
    ax_rmsd.grid(True, alpha=0.2, linewidth=0.5)
    ax_rmsd.spines["top"].set_visible(False)
    ax_rmsd.spines["right"].set_visible(False)
    ax_rmsd.legend(loc="upper right", fontsize=8, frameon=False, title="matched stringency", title_fontsize=8.2)

    # Panel (c): measured-fitness histogram at primary matched pair.
    amp_c = primary_pair["AMP"]
    azt_c = primary_pair["AZT"]
    amp_y = extras[f"fitness_AMP_{amp_c}"]
    azt_y = extras[f"fitness_AZT_{azt_c}"]
    bins = np.linspace(
        min(float(amp_y.min()), float(azt_y.min())) - 0.1,
        max(float(amp_y.max()), float(azt_y.max())) + 0.1,
        56,
    )
    ax_fit.hist(
        amp_y, bins=bins, alpha=0.55, color=AMP_PRIMARY,
        label=f"AMP {amp_c:g} µg/mL  (σ = {amp_y.std():.2f})",
    )
    ax_fit.hist(
        azt_y, bins=bins, alpha=0.55, color=AZT_PRIMARY,
        label=f"AZT {azt_c:g} µg/mL  (σ = {azt_y.std():.2f})",
    )
    ax_fit.set_xlabel("measured fitness (log₁₀ AUC)", fontsize=9.5)
    ax_fit.set_ylabel("count", fontsize=10)
    ax_fit.set_title(
        "(c) measured-fitness dynamic range",
        fontsize=10.5, fontweight="bold", loc="left",
    )
    ax_fit.legend(loc="upper left", fontsize=8, frameon=False)
    ax_fit.spines["top"].set_visible(False)
    ax_fit.spines["right"].set_visible(False)

    # Panel (d): |biochemical pairwise-epistasis| distribution at primary pair.
    amp_pair = extras[f"pairabs_AMP_{amp_c}"]
    azt_pair = extras[f"pairabs_AZT_{azt_c}"]
    lo = 0.0
    hi = float(np.quantile(np.concatenate([amp_pair, azt_pair]), 0.995))
    hi = max(hi, 0.5)
    bins_p = np.linspace(lo, hi * 1.05, 42)
    ax_pair.hist(
        amp_pair, bins=bins_p, alpha=0.55, color=AMP_PRIMARY,
        label=(
            f"AMP {amp_c:g} µg/mL  "
            f"(median = {np.median(amp_pair):.2f}, max = {amp_pair.max():.2f})"
        ),
    )
    ax_pair.hist(
        azt_pair, bins=bins_p, alpha=0.55, color=AZT_PRIMARY,
        label=(
            f"AZT {azt_c:g} µg/mL  "
            f"(median = {np.median(azt_pair):.2f}, max = {azt_pair.max():.2f})"
        ),
    )
    ax_pair.set_xlabel("|pairwise epistasis|  (Biochemical definition)", fontsize=9.5)
    ax_pair.set_ylabel("count", fontsize=10)
    ax_pair.set_title(
        "(d) pairwise-epistasis magnitude",
        fontsize=10.5, fontweight="bold", loc="left",
    )
    ax_pair.legend(loc="upper right", fontsize=8, frameon=False)
    ax_pair.spines["top"].set_visible(False)
    ax_pair.spines["right"].set_visible(False)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- Step 3: Source Data workbook --------------------------------------------
def write_source_data(
    df_order: pl.DataFrame,
    matched_long: pl.DataFrame,
    extras: dict,
    out_path: Path,
) -> None:
    """Assemble source_data.xlsx (one sheet per figure/table)."""
    amp_c = MATCHED_PAIRS["primary"]["AMP"]
    azt_c = MATCHED_PAIRS["primary"]["AZT"]

    sheets: dict[str, pd.DataFrame] = {}

    # Ext Fig A panels (one per drug) - long form with cum + dR^2.
    for drug in ("AMP", "AZT"):
        sub = (
            df_order
            .filter(pl.col("drug") == drug)
            .sort(["concentration", "order"])
            .rename({
                "concentration": "concentration_ug_per_mL",
                "order": "epistatic_order_K",
                "cum_r2": "cumulative_R2",
                "dr2": "incremental_dR2",
                "r2": "R2",
                "rmsd": "RMSD_log10_AUC",
            })
            .select(
                "drug",
                "concentration_ug_per_mL",
                "epistatic_order_K",
                "cumulative_R2",
                "incremental_dR2",
                "RMSD_log10_AUC",
                "n",
            )
        )
        sheets[f"Ext_Fig_A_{drug}"] = sub.to_pandas()

    # Ext Fig E - R^2 / RMSD curves at primary + secondary matched pairs.
    sheets["Ext_Fig_E_R2_curves"] = (
        matched_long
        .sort(["tier", "drug", "order"])
        .rename({
            "concentration": "concentration_ug_per_mL",
            "order": "epistatic_order_K",
            "r2": "R2",
            "rmsd": "RMSD_log10_AUC",
        })
        .to_pandas()
    )

    # Ext Fig E - fitness distributions at primary matched pair.
    amp_y = extras[f"fitness_AMP_{amp_c}"]
    azt_y = extras[f"fitness_AZT_{azt_c}"]
    fit_df = pd.DataFrame({
        "drug": ["AMP"] * len(amp_y) + ["AZT"] * len(azt_y),
        "concentration_ug_per_mL": [amp_c] * len(amp_y) + [azt_c] * len(azt_y),
        "measured_fitness_log10_AUC": np.concatenate([amp_y, azt_y]),
    })
    sheets["Ext_Fig_E_fitness_dist"] = fit_df

    # Ext Fig E - |pairwise epistasis| at primary matched pair.
    amp_pair = extras[f"pairabs_AMP_{amp_c}"]
    azt_pair = extras[f"pairabs_AZT_{azt_c}"]
    pair_df = pd.DataFrame({
        "drug": ["AMP"] * len(amp_pair) + ["AZT"] * len(azt_pair),
        "concentration_ug_per_mL": [amp_c] * len(amp_pair) + [azt_c] * len(azt_pair),
        "abs_biochemical_pairwise_epistasis":
            np.concatenate([amp_pair, azt_pair]),
    })
    sheets["Ext_Fig_E_pairwise_abs"] = pair_df

    # Summary of matched-stringency comparison at K=1,2,3,6.
    snapshot = (
        matched_long
        .filter(pl.col("order").is_in([1, 2, 3, 6]))
        .pivot(
            on="order",
            index=["tier", "drug", "concentration"],
            values=["r2", "rmsd"],
        )
        .sort(["tier", "drug", "concentration"])
    )
    sheets["Matched_stringency_summary"] = snapshot.to_pandas()

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for name, pdf in sheets.items():
            pdf.to_excel(writer, sheet_name=name[:31], index=False)
    print(f"wrote {out_path} (sheets: {list(sheets.keys())})")


# --- Step 4: Compact results table -------------------------------------------
def build_results_table(
    df_order: pl.DataFrame,
    matched_long: pl.DataFrame,
) -> pl.DataFrame:
    """Compact summary: matched-stringency R^2 and ΔR^2 at K=1,2,3,6,13."""
    rows: list[dict] = []
    for tier, spec in MATCHED_PAIRS.items():
        for drug in ("AMP", "AZT"):
            conc = spec[drug]
            sub = (
                df_order
                .filter(
                    (pl.col("drug") == drug) & (pl.col("concentration") == conc)
                )
                .sort("order")
            )
            r2_vec = {int(row["order"]): row["r2"] for row in sub.iter_rows(named=True)}
            rmsd_vec = {int(row["order"]): row["rmsd"] for row in sub.iter_rows(named=True)}
            row = {
                "tier": tier,
                "tier_label": spec["label"],
                "drug": drug,
                "concentration": conc,
                "R2_K1": round(r2_vec[1], 3),
                "R2_K2": round(r2_vec[2], 3),
                "R2_K3": round(r2_vec[3], 3),
                "R2_K6": round(r2_vec[6], 3),
                "R2_K13": round(r2_vec[13], 4),
                "RMSD_K1": round(rmsd_vec[1], 3),
                "RMSD_K6": round(rmsd_vec[6], 3),
                "RMSD_K13": round(rmsd_vec[13], 4),
            }
            rows.append(row)
    return pl.DataFrame(rows)


# --- main --------------------------------------------------------------------
def main() -> None:
    np.random.seed(SEED)
    print(f"reading {PARQUET}")
    df_all = pl.read_parquet(PARQUET)
    print(f"  {df_all.height} rows, {len(df_all.columns)} cols")

    # Order-decomposition table.
    df_order = compute_order_decomposition(df_all)
    df_order.write_csv(DATADIR / "order_decomposition.csv")
    print(f"wrote {DATADIR / 'order_decomposition.csv'} ({df_order.height} rows)")

    # Cross-check: our R^2 should equal S7's regression_r2_by_order.csv up to
    # floating-point noise.
    if S7_R2.exists():
        s7 = pl.read_csv(S7_R2)
        merged = (
            df_order
            .select("drug", "concentration", "order", "r2")
            .join(
                s7.select("drug", "concentration", "order", pl.col("r2").alias("r2_s7")),
                on=["drug", "concentration", "order"],
                how="inner",
            )
        )
        diff = (merged["r2"] - merged["r2_s7"]).abs().max()
        print(f"  cross-check max |R^2 - R^2_S7| = {diff:.3e}")
        assert diff < 1e-9, "S4 R^2 disagrees with S7 R^2"

    # Ext Fig A (single-drug + combined).
    plot_ext_fig_a_single(df_order, "AMP", FIGDIR / "ext_order_decomposition_amp.png")
    plot_ext_fig_a_single(df_order, "AZT", FIGDIR / "ext_order_decomposition_azt.png")
    plot_ext_fig_a_combined(df_order, FIGDIR / "ext_order_decomposition_combined.png")

    # Ext Fig E (matched stringency).
    matched_long, extras = compute_matched_stringency(df_all, df_order)
    matched_long.write_csv(DATADIR / "matched_stringency_summary.csv")
    print(f"wrote {DATADIR / 'matched_stringency_summary.csv'} ({matched_long.height} rows)")
    plot_ext_fig_e(matched_long, extras, df_order, FIGDIR / "ext_drug_asymmetry.png")

    # Source Data workbook.
    write_source_data(df_order, matched_long, extras, DATADIR / "source_data.xlsx")

    # Compact results table.
    results = build_results_table(df_order, matched_long)
    results.write_csv(HERE / "results_table.csv")
    print(f"wrote {HERE / 'results_table.csv'} ({results.height} rows)")
    print(results)


if __name__ == "__main__":
    main()
