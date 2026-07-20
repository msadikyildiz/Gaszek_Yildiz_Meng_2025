"""
S5 — Replicate reproducibility for Reviewer #1.

The manuscript (line 163) states that the AUC-Fitness metric "exhibited excellent
reproducibility across biological replicates, typically yielding standard
deviations below 10%." The reviewer requested that this be demonstrated
quantitatively (e.g. replicate-vs-replicate correlations).

This pipeline builds two complementary figures per drug:

  1. Replicate-pair scatter panels. For every drug x concentration, a 2D hex-bin
     density plot of replicate-vs-replicate log10(AUC) values with overlaid
     Pearson r and a parametric p-value. Three sub-panels per concentration for
     the three ordered pairs (rep1-rep2, rep1-rep3, rep2-rep3). Conc=0 is
     dropped from the figure (uninformative: flat bacterial growth) but kept in
     the numerical tables.

  2. CV / SD distribution panels. For every genotype and every drug x
     concentration we compute:
       (a) %CV on raw AUC (back-transformed via 10^fitness) -- matches the
           manuscript's "<10%" wording, which is conventionally a coefficient
           of variation on a linear-scale quantity.
       (b) SD of log10(AUC) directly -- matches manuscript units (fitness is
           defined as log10 AUC in the paper).
     Per-concentration histograms plus median / IQR summary.

Inputs:
  data/raw/Ampicillin_auc_per_genotype.csv
  data/raw/Aztreonam_auc_per_genotype.csv
  (Schema: mut_profile_masked, <Drug> <conc> 1, <Drug> <conc> 2, <Drug> <conc> 3.)

Outputs:
  figures/replicate_scatter_amp.png
  figures/replicate_scatter_azt.png
  figures/replicate_cv_amp.png
  figures/replicate_cv_azt.png
  data/per_genotype_replicate_stats.csv     (one row per genotype x drug x conc)
  data/summary_by_conc.csv                  (one row per drug x conc)
  results_table.csv                         (headline claim-verdict table)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from scipy import stats

# --- paths -------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
RAW = REPO / "data" / "raw"
FIGDIR = HERE / "figures"
DATADIR = HERE / "data"
FIGDIR.mkdir(exist_ok=True, parents=True)
DATADIR.mkdir(exist_ok=True, parents=True)

# --- constants ---------------------------------------------------------------
DRUG_FILES = {
    "Ampicillin": RAW / "Ampicillin_auc_per_genotype.csv",
    "Aztreonam": RAW / "Aztreonam_auc_per_genotype.csv",
}
AMP_CONCS = [0.0, 3.1, 12.2, 48.8, 195.0, 781.0]
AZT_CONCS = [0.0, 0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]
DRUG_CONCS = {"Ampicillin": AMP_CONCS, "Aztreonam": AZT_CONCS}

# Nonzero concs only -- conc=0 carries no meaningful biological signal
# (pure growth control) and swamps the figure scale.
FIG_CONCS = {
    "Ampicillin": [c for c in AMP_CONCS if c > 0],
    "Aztreonam": [c for c in AZT_CONCS if c > 0],
}

SEED = 20260420
DRUG_COLOR = {"Ampicillin": "#b02a2a", "Aztreonam": "#1f6fb5"}

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]


# --- IO + per-genotype stats -------------------------------------------------
def load_drug(drug: str) -> pl.DataFrame:
    df = pl.read_csv(DRUG_FILES[drug])
    df = df.rename({df.columns[0]: "_row_index"})
    return df


def _reps_for(df: pl.DataFrame, drug: str, conc: float) -> np.ndarray:
    """Return n_genotypes x 3 array of log10(AUC) replicate values."""
    cols = [f"{drug} {conc} {r}" for r in (1, 2, 3)]
    arr = df.select(cols).to_numpy().astype(np.float64)
    return arr


def per_genotype_stats(df: pl.DataFrame, drug: str, conc: float) -> pl.DataFrame:
    """Compute per-genotype replicate statistics.

    Columns:
        genotype          -- mut_profile_masked string
        drug, concentration
        log_mean, log_sd  -- mean / SD of log10(AUC) across 3 reps
        raw_mean, raw_sd  -- mean / SD of raw AUC = 10^log10
        cv_percent        -- 100 * raw_sd / raw_mean
    """
    reps = _reps_for(df, drug, conc)                       # n x 3 log-space
    raw = np.power(10.0, reps)                             # n x 3 linear AUC
    log_mean = np.nanmean(reps, axis=1)
    log_sd = np.nanstd(reps, axis=1, ddof=1)
    raw_mean = np.nanmean(raw, axis=1)
    raw_sd = np.nanstd(raw, axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = 100.0 * raw_sd / raw_mean
    out = pl.DataFrame({
        "genotype": df["mut_profile_masked"].to_numpy(),
        "drug": [drug] * len(reps),
        "concentration": [conc] * len(reps),
        "log_mean": log_mean,
        "log_sd": log_sd,
        "raw_mean": raw_mean,
        "raw_sd": raw_sd,
        "cv_percent": cv,
    })
    return out


def build_per_genotype_table() -> pl.DataFrame:
    rows: list[pl.DataFrame] = []
    for drug, concs in DRUG_CONCS.items():
        df = load_drug(drug)
        for c in concs:
            rows.append(per_genotype_stats(df, drug, c))
    return pl.concat(rows)


def summarise_by_conc(per_gt: pl.DataFrame) -> pl.DataFrame:
    """Per drug x concentration: n, median/IQR of CV% and log_sd, fraction <10% CV.

    Restricts to genotypes for which the per-genotype SD is well-defined
    (i.e. at least 2 of the 3 replicates are non-null). This matches the
    denominator shown in the histograms; genotypes with only one viable
    replicate have no dispersion estimate at all.
    """
    return (
        per_gt
        .filter(pl.col("log_sd").is_finite())
        .group_by(["drug", "concentration"], maintain_order=True)
        .agg([
            pl.len().alias("n_genotypes"),
            pl.col("cv_percent").median().alias("cv_median"),
            pl.col("cv_percent").quantile(0.25).alias("cv_q25"),
            pl.col("cv_percent").quantile(0.75).alias("cv_q75"),
            pl.col("cv_percent").quantile(0.90).alias("cv_q90"),
            (pl.col("cv_percent") < 10.0).cast(pl.Float64).mean().alias("frac_lt_10pct_cv"),
            pl.col("log_sd").median().alias("logsd_median"),
            pl.col("log_sd").quantile(0.25).alias("logsd_q25"),
            pl.col("log_sd").quantile(0.75).alias("logsd_q75"),
        ])
        .sort(["drug", "concentration"])
    )


# --- replicate-pair correlations ---------------------------------------------
def pair_correlations(df: pl.DataFrame, drug: str, conc: float) -> list[dict]:
    reps = _reps_for(df, drug, conc)
    # Drop rows with any NaN so Pearson is well-defined
    mask = np.isfinite(reps).all(axis=1)
    reps = reps[mask]
    out: list[dict] = []
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        x = reps[:, i]
        y = reps[:, j]
        r, p = stats.pearsonr(x, y)
        out.append({
            "drug": drug,
            "concentration": conc,
            "rep_i": i + 1,
            "rep_j": j + 1,
            "n": int(len(x)),
            "pearson_r": r,
            "pearson_p": p,
        })
    return out


def build_correlation_table() -> pl.DataFrame:
    rows: list[dict] = []
    for drug, concs in DRUG_CONCS.items():
        df = load_drug(drug)
        for c in concs:
            rows.extend(pair_correlations(df, drug, c))
    return pl.DataFrame(rows)


# --- plotting: replicate scatter ---------------------------------------------
def plot_replicate_scatter(drug: str, out_path: Path) -> None:
    df = load_drug(drug)
    concs = FIG_CONCS[drug]
    n_conc = len(concs)
    fig, axes = plt.subplots(n_conc, 3, figsize=(9.0, 2.7 * n_conc), dpi=150,
                              constrained_layout=False)
    if n_conc == 1:
        axes = np.array([axes])

    col = DRUG_COLOR[drug]
    # Determine global log10(AUC) range per drug for consistent axes
    reps_all = np.concatenate(
        [_reps_for(df, drug, c).ravel() for c in concs]
    )
    reps_all = reps_all[np.isfinite(reps_all)]
    lo, hi = np.nanpercentile(reps_all, [0.2, 99.8])
    pad = 0.05 * (hi - lo)
    lo, hi = lo - pad, hi + pad

    for i, conc in enumerate(concs):
        reps = _reps_for(df, drug, conc)
        mask = np.isfinite(reps).all(axis=1)
        reps = reps[mask]
        pairs = [(0, 1), (0, 2), (1, 2)]
        for j, (a, b) in enumerate(pairs):
            ax = axes[i, j]
            x, y = reps[:, a], reps[:, b]
            hb = ax.hexbin(
                x, y, gridsize=55, cmap="Greys", bins="log",
                mincnt=1, linewidths=0, extent=(lo, hi, lo, hi),
            )
            # Identity line
            ax.plot([lo, hi], [lo, hi], color=col, linewidth=0.8, alpha=0.85)
            # Pearson annotation
            r, p = stats.pearsonr(x, y)
            p_str = "p < 1e-300" if p < 1e-300 else f"p = {p:.1e}"
            ax.text(
                0.035, 0.965,
                f"r = {r:.3f}\n{p_str}\nn = {len(x):,}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=7.2, color="#222222",
                bbox=dict(facecolor="white", alpha=0.80,
                          edgecolor="none", pad=1.4),
            )
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_aspect("equal")
            ax.tick_params(labelsize=7)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if j == 0:
                ax.set_ylabel(f"{conc:g} µg/mL\nrep {b+1}", fontsize=8)
            else:
                ax.set_ylabel(f"rep {b+1}", fontsize=8)
            if i == n_conc - 1:
                ax.set_xlabel(f"rep {a+1}  (log$_{{10}}$ AUC)", fontsize=8)
            else:
                ax.set_xlabel("")

    fig.suptitle(
        f"{drug} -- replicate-replicate log$_{{10}}$(AUC) scatter "
        f"(hex-bin, log counts; identity line)",
        fontsize=10.5, y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- plotting: CV distributions on the viable subset -------------------------
VIABLE_LOG_MEAN_THRESHOLD = 3.0


def plot_cv_viable(drug: str, per_gt: pl.DataFrame, out_path: Path) -> None:
    """Per-concentration %CV histograms restricted to the viable subset
    (genotypes with log10(AUC) mean > 3.0 at the given condition). Mirrors
    plot_cv_distributions but drops the dispersion-diagnostic log-SD column
    so the panels fit the new headline framing.
    """
    concs = FIG_CONCS[drug]
    n_conc = len(concs)
    fig, axes = plt.subplots(n_conc, 1, figsize=(5.8, 1.5 * n_conc), dpi=150,
                              squeeze=False)
    col = DRUG_COLOR[drug]

    for i, conc in enumerate(concs):
        sub = per_gt.filter(
            (pl.col("drug") == drug) & (pl.col("concentration") == conc)
            & (pl.col("log_mean") > VIABLE_LOG_MEAN_THRESHOLD)
        )
        cv = sub["cv_percent"].to_numpy()
        cv = cv[np.isfinite(cv)]
        n_v = len(cv)

        ax = axes[i, 0]
        cv_clip = np.clip(cv, 0, 100)
        n_above = int((cv > 100).sum())
        bins_cv = np.linspace(0, 100, 51)
        ax.hist(cv_clip, bins=bins_cv, color=col, alpha=0.82,
                edgecolor="white", linewidth=0.25)
        ax.axvline(10.0, color="#222222", linestyle="--", linewidth=0.8,
                   alpha=0.9)
        med = float(np.median(cv)) if n_v else float("nan")
        q25 = float(np.quantile(cv, 0.25)) if n_v else float("nan")
        q75 = float(np.quantile(cv, 0.75)) if n_v else float("nan")
        frac10 = float((cv < 10.0).mean()) if n_v else 0.0

        txt = (
            f"{conc:g} µg/mL  n_viable = {n_v:,}\n"
            f"median CV = {med:.1f}%  "
            f"IQR [{q25:.1f}, {q75:.1f}]\n"
            f"{frac10*100:.1f}% of viable genotypes < 10% CV"
        )
        if n_above:
            txt += f"\n({n_above:,} viable genotypes > 100% CV clipped)"
        ax.text(
            0.98, 0.96, txt, transform=ax.transAxes,
            ha="right", va="top", fontsize=7.2,
            bbox=dict(facecolor="white", alpha=0.85,
                      edgecolor="none", pad=1.6),
        )
        ax.set_xlim(0, 100)
        ax.set_xlabel("%CV of raw AUC across 3 replicates", fontsize=8)
        ax.set_ylabel("viable\ngenotypes", fontsize=8)
        ax.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.suptitle(
        f"{drug} — viable-subset replicate dispersion "
        f"(log$_{{10}}$ AUC mean > {VIABLE_LOG_MEAN_THRESHOLD:g}; "
        f"dashed line: 10% CV)",
        fontsize=10.5, y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- plotting: CV / log-SD distributions -------------------------------------
def plot_cv_distributions(drug: str, per_gt: pl.DataFrame, out_path: Path) -> None:
    concs = FIG_CONCS[drug]
    n_conc = len(concs)
    fig, axes = plt.subplots(n_conc, 2, figsize=(8.2, 1.6 * n_conc), dpi=150)
    if n_conc == 1:
        axes = np.array([axes])

    col = DRUG_COLOR[drug]

    for i, conc in enumerate(concs):
        sub = per_gt.filter(
            (pl.col("drug") == drug) & (pl.col("concentration") == conc)
        )
        cv = sub["cv_percent"].to_numpy()
        log_sd = sub["log_sd"].to_numpy()
        cv = cv[np.isfinite(cv)]
        log_sd = log_sd[np.isfinite(log_sd)]

        ax_cv, ax_sd = axes[i, 0], axes[i, 1]

        # --- %CV histogram ---
        # Clip visual right-tail at 100 %CV so the bulk is visible; report
        # tail fraction in-panel.
        cv_clip = np.clip(cv, 0, 100)
        n_above = int((cv > 100).sum())
        bins_cv = np.linspace(0, 100, 51)
        ax_cv.hist(cv_clip, bins=bins_cv, color=col, alpha=0.78,
                   edgecolor="white", linewidth=0.25)
        ax_cv.axvline(10.0, color="#222222", linestyle="--", linewidth=0.8,
                      alpha=0.9)
        med = float(np.median(cv))
        q25 = float(np.quantile(cv, 0.25))
        q75 = float(np.quantile(cv, 0.75))
        frac10 = float((cv < 10.0).mean())
        txt = (
            f"{conc:g} µg/mL  n = {len(cv):,}\n"
            f"median CV = {med:.1f}%  "
            f"IQR [{q25:.1f}, {q75:.1f}]\n"
            f"{frac10*100:.1f}% of genotypes < 10% CV"
        )
        if n_above:
            txt += f"\n({n_above:,} genotypes > 100% CV clipped)"
        ax_cv.text(0.98, 0.96, txt, transform=ax_cv.transAxes,
                   ha="right", va="top", fontsize=7.2,
                   bbox=dict(facecolor="white", alpha=0.85,
                             edgecolor="none", pad=1.6))
        ax_cv.set_xlim(0, 100)
        ax_cv.set_xlabel("%CV of raw AUC across 3 replicates", fontsize=8)
        ax_cv.set_ylabel("genotypes", fontsize=8)
        ax_cv.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax_cv.spines[s].set_visible(False)

        # --- log-SD histogram ---
        log_sd_clip = np.clip(log_sd, 0, 0.6)
        bins_sd = np.linspace(0, 0.6, 61)
        ax_sd.hist(log_sd_clip, bins=bins_sd, color=col, alpha=0.78,
                   edgecolor="white", linewidth=0.25)
        med_sd = float(np.median(log_sd))
        q25_sd = float(np.quantile(log_sd, 0.25))
        q75_sd = float(np.quantile(log_sd, 0.75))
        txt_sd = (
            f"median SD = {med_sd:.3f}\n"
            f"IQR [{q25_sd:.3f}, {q75_sd:.3f}]"
        )
        ax_sd.text(0.98, 0.96, txt_sd, transform=ax_sd.transAxes,
                   ha="right", va="top", fontsize=7.2,
                   bbox=dict(facecolor="white", alpha=0.85,
                             edgecolor="none", pad=1.6))
        ax_sd.set_xlim(0, 0.6)
        ax_sd.set_xlabel("SD of log$_{10}$(AUC) across 3 replicates",
                         fontsize=8)
        ax_sd.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax_sd.spines[s].set_visible(False)

    fig.suptitle(
        f"{drug} -- replicate dispersion by concentration  "
        f"(dashed line: 10% CV threshold)",
        fontsize=10.5, y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- headline verdict table --------------------------------------------------
def build_results_table(summary: pl.DataFrame, corr: pl.DataFrame) -> pl.DataFrame:
    """Join the per-conc CV summary with the mean Pearson r across pairs."""
    corr_mean = (
        corr.group_by(["drug", "concentration"])
        .agg([
            pl.col("pearson_r").mean().alias("pearson_r_mean_over_pairs"),
            pl.col("pearson_r").min().alias("pearson_r_min_over_pairs"),
        ])
    )
    out = (
        summary
        .join(corr_mean, on=["drug", "concentration"])
        .with_columns([
            (pl.col("cv_median") < 10.0).alias("median_cv_below_10pct"),
        ])
        .sort(["drug", "concentration"])
    )
    return out


# --- main --------------------------------------------------------------------
def main():
    np.random.seed(SEED)
    print("Computing per-genotype replicate statistics...")
    per_gt = build_per_genotype_table()
    per_gt.write_csv(DATADIR / "per_genotype_replicate_stats.csv")
    print(f"  wrote {DATADIR / 'per_genotype_replicate_stats.csv'}  "
          f"({per_gt.height:,} rows)")

    print("Summarising by drug x concentration...")
    summary = summarise_by_conc(per_gt)
    summary.write_csv(DATADIR / "summary_by_conc.csv")
    print(f"  wrote {DATADIR / 'summary_by_conc.csv'}  ({summary.height} rows)")

    print("Computing replicate-pair Pearson correlations...")
    corr = build_correlation_table()
    corr.write_csv(DATADIR / "replicate_pair_correlations.csv")
    print(f"  wrote {DATADIR / 'replicate_pair_correlations.csv'}  "
          f"({corr.height} rows)")

    print("Writing headline verdict table...")
    headline = build_results_table(summary, corr)
    headline.write_csv(HERE / "results_table.csv")
    print(f"  wrote {HERE / 'results_table.csv'}  ({headline.height} rows)")

    print("Drawing replicate-pair scatter panels...")
    for drug in DRUG_FILES:
        plot_replicate_scatter(drug, FIGDIR / f"replicate_scatter_{drug[:3].lower()}.png")

    print("Drawing CV / log-SD distribution panels...")
    for drug in DRUG_FILES:
        plot_cv_distributions(drug, per_gt,
                               FIGDIR / f"replicate_cv_{drug[:3].lower()}.png")

    print("Drawing viable-subset %CV panels (log10 AUC mean > 3.0)...")
    for drug in DRUG_FILES:
        plot_cv_viable(drug, per_gt,
                       FIGDIR / f"replicate_cv_viable_{drug[:3].lower()}.png")

    # Short stdout summary
    print("\n=== Headline (median CV% by drug x conc, nonzero concs) ===")
    for drug in DRUG_FILES:
        s = headline.filter(
            (pl.col("drug") == drug) & (pl.col("concentration") > 0)
        )
        for row in s.iter_rows(named=True):
            flag = "PASS" if row["median_cv_below_10pct"] else "FAIL"
            print(
                f"  {drug[:3]}  {row['concentration']:>7g} ug/mL  "
                f"n={row['n_genotypes']:>6}  "
                f"median CV {row['cv_median']:>5.1f}%  "
                f"IQR [{row['cv_q25']:>4.1f}, {row['cv_q75']:>4.1f}]  "
                f"frac<10%CV {row['frac_lt_10pct_cv']*100:>5.1f}%  "
                f"mean r={row['pearson_r_mean_over_pairs']:.3f}  "
                f"[{flag}]"
            )


if __name__ == "__main__":
    main()
