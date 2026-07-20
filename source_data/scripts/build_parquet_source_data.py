#!/usr/bin/env python3
"""Extract the parquet-derived Source Data tables (Fig 1E, Fig 4E-F, Fig S12) from
Epistasis_Combined.parquet, matching exactly what each panel plots in
src/05_epistasis_figures.ipynb (Fig 4E-F) and analysis/s11_s12_concentration_grid (S12).

Outputs land in source_data/derived/ as CSVs that build_source_data.py folds
into source_data.xlsx.

Panel definitions verified against the notebook / analysis code:
- Fig 1E: ridgeline of AUC-fitness (Epistasis_Combined 'Fitness') per drug x concentration
  (AMP 6 concs incl 0; AZT 8 concs incl 0). Black dashed = TEM-1(WT) median; blue dashed =
  TEM-1(dead) median, both from *_auc_long_df.parquet ('.'*13 and 'X'*13 profiles).
- Fig 4E-F: measured ('Fitness') vs cumulative predicted ('Fitness_predicted for order K',
  K=1..6) at the reference concentration (AMP 781, AZT 36). R2 = Pearson squared.
- Fig S12: same measured-vs-predicted, K=1..6, across the s11/s12 concentration set (0-drug
  panels excluded). R2 = Pearson squared (calc_r2_rmsd). Cross-checked vs the analysis CSV.
"""
import csv
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

HERE = Path(__file__).resolve().parent


def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))


REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / "data/processed/Epistasis_Combined.parquet"
AMP_LONG = REPO / "data/processed/amp_auc_long_df.parquet"
AZT_LONG = REPO / "data/processed/azt_auc_long_df.parquet"
S7_R2 = REPO / "analysis/s11_s12_concentration_grid/data/regression_r2_by_order.csv"

OUT = HERE.parent / "derived"
OUT.mkdir(exist_ok=True)

ORDERS = [1, 2, 3, 4, 5, 6]
REF_CONC = {"AMP": 781.0, "AZT": 36.0}                 # Fig 4E-F reference doses
# concentration grid for S12 (0-drug excluded), matching analysis/s11_s12_concentration_grid.
S12_CONCS = {"AMP": [3.1, 12.2, 48.8, 195.0, 781.0],
             "AZT": [0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]}
WT_PROFILE = "." * 13
DEAD_PROFILE = "X" * 13


def pearson_r2_rmsd(m, p):
    r = np.corrcoef(m, p)[0, 1]
    rmsd = float(np.sqrt(np.mean((m - p) ** 2)))
    return float(r), float(r * r), rmsd


def write_csv(name, header, rows):
    path = OUT / name
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {name}  ({len(rows)} rows)")
    return path


def ref_medians(long_path):
    """Per-concentration WT and dead medians (the dashed reference lines)."""
    d = pl.read_parquet(long_path, columns=["mutant_profile", "concentration", "median"])
    wt = {r["concentration"]: r["median"]
          for r in d.filter(pl.col("mutant_profile") == WT_PROFILE).to_dicts()}
    dead = {r["concentration"]: r["median"]
            for r in d.filter(pl.col("mutant_profile") == DEAD_PROFILE).to_dicts()}
    return wt, dead


def build_fig1e(df):
    print("Fig 1E (ridgeline distributions):")
    edges = np.linspace(0.0, 6.0, 121)                 # 120 bins, width 0.05, matches axis 0-6
    centers = (edges[:-1] + edges[1:]) / 2
    refs = {"AMP": ref_medians(AMP_LONG), "AZT": ref_medians(AZT_LONG)}

    summary = []
    hist_cols = {}                                     # "AMP_0" -> counts
    col_order = []
    for drug in ("AMP", "AZT"):
        wt_map, dead_map = refs[drug]
        concs = sorted(df.filter(pl.col("Drug") == drug)["Concentration"].unique().to_list())
        for c in concs:
            fit = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == c))["Fitness"].to_numpy()
            fit = fit[np.isfinite(fit)]
            counts, _ = np.histogram(fit, bins=edges)
            key = f"{drug}_{c:g}"
            hist_cols[key] = counts
            col_order.append(key)
            # WT cross-check: Epistasis_Combined WT genotype fitness vs auc_long_df WT median
            wt_ec = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == c)
                              & (pl.col("Epistatic Order") == 0))["Fitness"]
            wt_ec = float(wt_ec[0]) if wt_ec.len() else float("nan")
            if not np.isnan(wt_ec) and abs(wt_ec - wt_map[c]) > 0.05:
                print(f"    WARN WT mismatch {drug} {c}: EC={wt_ec:.3f} vs long={wt_map[c]:.3f}")
            summary.append([drug, f"{c:g}", len(fit),
                            round(float(fit.min()), 4), round(float(fit.max()), 4),
                            round(float(fit.mean()), 4), round(float(np.median(fit)), 4),
                            round(float(fit.std(ddof=1)), 4),
                            round(wt_map[c], 4), round(dead_map[c], 4)])

    write_csv("fig1E_distribution_summary.csv",
              ["drug", "concentration_ug_ml", "n_genotypes", "fitness_min", "fitness_max",
               "fitness_mean", "fitness_median", "fitness_sd",
               "WT_median_black_dashed", "dead_median_blue_dashed"], summary)

    hrows = [[round(float(centers[i]), 4)] + [int(hist_cols[k][i]) for k in col_order]
             for i in range(len(centers))]
    write_csv("fig1E_histogram_counts.csv",
              ["fitness_bin_center"] + col_order, hrows)


def build_fig4ef(df):
    print("Fig 4E-F (measured vs predicted, reference dose):")
    summary, raw = [], []
    for drug in ("AMP", "AZT"):
        c = REF_CONC[drug]
        sub = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == c))
        m = sub["Fitness"].to_numpy()
        preds = {k: sub[f"Fitness_predicted for order {k}"].to_numpy() for k in ORDERS}
        gts = sub["Genotype"].to_list()
        for k in ORDERS:
            r, r2, rmsd = pearson_r2_rmsd(m, preds[k])
            lr = stats.linregress(m, preds[k])
            summary.append([drug, f"{c:g}", k, len(m), round(r, 4), round(r2, 4),
                            round(rmsd, 4), round(float(lr.slope), 4), round(float(lr.intercept), 4)])
        for i, g in enumerate(gts):
            raw.append([g, drug, f"{c:g}", round(float(m[i]), 4)]
                       + [round(float(preds[k][i]), 4) for k in ORDERS])
    write_csv("fig4EF_summary.csv",
              ["drug", "concentration_ug_ml", "max_epistatic_order", "n", "pearson_r",
               "R2_pearson_sq", "rmsd", "fit_slope", "fit_intercept"], summary)
    write_csv("fig4EF_measured_vs_predicted.csv",
              ["genotype", "drug", "concentration_ug_ml", "fitness_measured"]
              + [f"predicted_order_le_{k}" for k in ORDERS], raw)


def build_figs12(df):
    print("Fig S12 (measured vs predicted across concentrations):")
    # cross-check reference
    ref = {}
    with open(S7_R2) as fh:
        for row in csv.DictReader(fh):
            ref[(row["drug"], float(row["concentration"]), int(row["order"]))] = (
                float(row["r2"]), float(row["rmsd"]))
    summary, maxdev = [], 0.0
    for drug in ("AMP", "AZT"):
        for c in S12_CONCS[drug]:
            sub = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == c))
            m = sub["Fitness"].to_numpy()
            for k in ORDERS:
                p = sub[f"Fitness_predicted for order {k}"].to_numpy()
                _, r2, rmsd = pearson_r2_rmsd(m, p)
                if (drug, c, k) in ref:
                    dev = abs(r2 - ref[(drug, c, k)][0]) + abs(rmsd - ref[(drug, c, k)][1])
                    maxdev = max(maxdev, dev)
                summary.append([drug, f"{c:g}", k, len(m), round(r2, 4), round(rmsd, 4)])
    print(f"  max |deviation| vs analysis/s11_s12_concentration_grid regression_r2_by_order.csv = {maxdev:.2e}")
    assert maxdev < 1e-3, "S12 recompute does not match the analysis/s11_s12_concentration_grid CSV"
    write_csv("figS12_r2_by_order_summary.csv",
              ["drug", "concentration_ug_ml", "max_epistatic_order", "n", "R2_pearson_sq", "rmsd"],
              summary)


def main():
    df = pl.read_parquet(PARQUET)
    print(f"loaded {PARQUET.name}: {df.height} rows\n")
    build_fig1e(df)
    build_fig4ef(df)
    build_figs12(df)
    print("\ndone.")


if __name__ == "__main__":
    main()
