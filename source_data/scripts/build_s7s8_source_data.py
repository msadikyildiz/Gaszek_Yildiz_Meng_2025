#!/usr/bin/env python3
"""Source Data for Supplementary Figs S7 (monoculture IC50 vs mean AUC-fitness) and
S8 (per-concentration AUC-fitness vs IC50), from the exact validation parquets the
manuscript panels were plotted from (validation/src/fig_c_mic_vs_fitness.py):
AMP batch 20260407, AZT batch 20260307, xref_expanded_df.parquet.

S7 (plot_ic50_vs_auc): y = log10(IC50), x = mean_fitness, per variant, per drug + Pearson r/p.
S8 (plot_ic50_per_conc): per variant, y = log10(IC50), x = AUC-fitness at each concentration.
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
PROC = REPO / "validation/src/processed"
BATCH = {"AMP": "20260407", "AZT": "20260307"}
OUT = HERE.parent / "derived"
OUT.mkdir(exist_ok=True)


def write_csv(name, header, rows):
    with open(OUT / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {name}  ({len(rows)} rows)")


def load(drug):
    return pl.read_parquet(PROC / BATCH[drug] / "xref_expanded_df.parquet")


def r(x):
    return None if x is None else round(float(x), 4)


def pear(xs, ys):
    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return int(ok.sum()), None, None
    pr, pp = stats.pearsonr(x[ok], y[ok])
    return int(ok.sum()), round(float(pr), 4), f"{pp:.2e}"


def build():
    s7, s8, corr7, corr8 = [], [], [], []
    for drug in ("AMP", "AZT"):
        d = load(drug)
        for row in d.to_dicts():
            s7.append([drug, row["variant"], row["label"], row["genotype_13"], row["n_mutations"],
                       r(row["mean_fitness"]), r(row["log10_ic50"]), r(row["ic50_mean"]), r(row["mic_ugml"])])
        # S7 Pearson r/p on the pairs actually plotted (both defined)
        n, pr, pp = pear(d["mean_fitness"].to_list(), d["log10_ic50"].to_list())
        corr7.append([drug, n, pr, pp])
        # S8: per-concentration long form + per-conc Pearson (0-drug excluded; no S8 panel)
        fcols = sorted((c for c in d.columns if c.startswith("fitness_")),
                       key=lambda c: float(c.split("_", 1)[1]))
        li = d["log10_ic50"].to_list()
        for c in fcols:
            conc = c.split("_", 1)[1]
            if float(conc) == 0.0:
                continue
            n, pr, pp = pear(d[c].to_list(), li)
            corr8.append([drug, f"{float(conc):g}", n, pr, pp])
            for row in d.to_dicts():
                s8.append([drug, row["variant"], row["genotype_13"], f"{float(conc):g}",
                           r(row[c]), r(row["log10_ic50"]), r(row["ic50_mean"])])

    write_csv("figS7_ic50_vs_meanfitness.csv",
              ["drug", "variant", "label", "genotype_13", "n_mutations", "mean_fitness",
               "log10_ic50", "ic50_mean_ugml", "mic_ugml"], s7)
    write_csv("figS7_pearson_correlation.csv",
              ["drug", "n_variants", "pearson_r", "p_value"], corr7)
    write_csv("figS8_ic50_vs_fitness_per_conc.csv",
              ["drug", "variant", "genotype_13", "concentration_ugml", "auc_fitness_at_conc",
               "log10_ic50", "ic50_mean_ugml"], s8)
    write_csv("figS8_pearson_by_conc.csv",
              ["drug", "concentration_ugml", "n_variants", "pearson_r", "p_value"], corr8)
    print("  S7 corr:", corr7)
    print("  S8 corr:", corr8)


if __name__ == "__main__":
    build()
