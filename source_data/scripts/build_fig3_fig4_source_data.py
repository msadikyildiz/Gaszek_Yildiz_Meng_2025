#!/usr/bin/env python3
"""Extract Source Data for Fig 3 (A/B/C/D) and Fig 4 (A/B/C/D) from the parquet
matrices, matching what src/05_epistasis_figures.ipynb plots.

Verified panel definitions:
- Fig 3A: AZT-vs-AMP density-scatter at matched doses (AMP 195/781 x AZT 36/108/324) with
  WT / dead / 7 engineered follow-ups / clinical-isolate overlays (cell 10 + legend).
- Fig 3B: sequence logos from the top 1% of genotypes at AMP 781 (left) / AZT 36 (right);
  letter height = per-position residue frequency among that top 1% (cells 19-21).
- Fig 3C/D: fitness vs number of mutations at AMP 781 / AZT 36, with median-per-order line
  (cells 11-12).
- Fig 4A/B: single-mutant fitness at AMP 781 / AZT 36 with n=3 selection replicates (legend;
  replicate values from *_auc_long_df.parquet).
- Fig 4C/D: heatmap of AUC-Fitness for every single/double combination of the 18 mutations
  at AMP 781 / AZT 36 (create_summary_heatmap, cell 32 -> row['Fitness'], NOT epistasis).
"""
import csv
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent


def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))


REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / "data/processed/Epistasis_Combined.parquet"
LONG = {"AMP": REPO / "data/processed/amp_auc_long_df.parquet",
        "AZT": REPO / "data/processed/azt_auc_long_df.parquet"}
CLINICAL = REPO / "data/known_variants/encoded_variants.csv"
OUT = HERE.parent / "derived"
OUT.mkdir(exist_ok=True)

WT_SEQ = "LQMERMAGERTRN"
AMBLER = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
REF = {"AMP": 781.0, "AZT": 36.0}
FIG3A_CONCS = [("AMP", 195.0), ("AMP", 781.0), ("AZT", 36.0), ("AZT", 108.0), ("AZT", 324.0)]
FIG3A_COLS = [f"{d}_{c:g}" for d, c in FIG3A_CONCS]
FOLLOWUPS = {  # 7 engineered follow-up variants (cell 10), letter-encoded
    "c1.1": "LQMKNTAGKRTRN", "c1.2": "PQMKNMAGKRMRN", "c1.3": "LKMKSMAGKRMRN",
    "c2.1": "LQMKNMAGKRMRN", "c2.2": "LKMKNMAGKRTRN", "c2.3": "LKMKNTAGKRTRN",
    "c3.1": "PQMKNTAGKRTRN"}


def mut_label(seq):
    """Genotype (letter-encoded) -> 'WT' or '; '-joined mutation names (e.g. E104K)."""
    muts = [f"{WT_SEQ[i]}{AMBLER[i]}{seq[i]}" for i in range(13) if seq[i] != WT_SEQ[i]]
    return "; ".join(muts) if muts else "WT"


def write_csv(name, header, rows):
    with open(OUT / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {name}  ({len(rows)} rows)")


def fitness_lookup(df):
    """(drug, conc) -> {genotype: fitness}."""
    out = {}
    for d, c in FIG3A_CONCS:
        sub = df.filter((pl.col("Drug") == d) & (pl.col("Concentration") == c))
        out[(d, c)] = dict(zip(sub["Genotype"].to_list(), sub["Fitness"].to_list()))
    return out


def dead_medians():
    """(drug, conc) -> dead-control ('X'*13) median."""
    out = {}
    for d, c in FIG3A_CONCS:
        lg = pl.read_parquet(LONG[d], columns=["mutant_profile", "concentration", "median"])
        row = lg.filter((pl.col("mutant_profile") == "X" * 13) & (pl.col("concentration") == c))
        out[(d, c)] = float(row["median"][0]) if row.height else None
    return out


def build_fig3a(df):
    print("Fig 3A (overlay points + all-genotype density source):")
    look = fitness_lookup(df)
    dead = dead_medians()

    # all-genotype companion (55,296 x 5 doses)
    genos = df.filter((pl.col("Drug") == "AMP") & (pl.col("Concentration") == 781.0))["Genotype"].to_list()
    rows = [[g] + [round(look[(d, c)].get(g, float("nan")), 4) for d, c in FIG3A_CONCS] for g in genos]
    write_csv("fig3A_all_genotypes.csv", ["genotype"] + FIG3A_COLS, rows)

    # overlay points
    overlay = [("TEM-1(WT)", WT_SEQ, "reference")]
    for lbl, seq in FOLLOWUPS.items():
        overlay.append((lbl, seq, "engineered follow-up"))
    if CLINICAL.exists():
        with open(CLINICAL) as fh:
            for r in csv.DictReader(fh):
                seq = (r.get("Encoded_Sequence") or "").strip()
                if len(seq) == 13:
                    overlay.append((r["Variant"], seq, "clinical isolate"))
    orows = []
    for lbl, seq, kind in overlay:
        orows.append([lbl, seq, kind] + [round(look[(d, c)].get(seq, float("nan")), 4) for d, c in FIG3A_CONCS])
    # dead control (not a library genotype; medians from long-df)
    orows.append(["TEM-1(dead)", "X" * 13, "dead control"]
                 + [round(dead[(d, c)], 4) if dead[(d, c)] is not None else "" for d, c in FIG3A_CONCS])
    write_csv("fig3A_overlay_points.csv", ["label", "genotype", "category"] + FIG3A_COLS, orows)


def build_fig3b(df):
    print("Fig 3B (top-1% sequence-logo frequencies):")
    for drug in ("AMP", "AZT"):
        sub = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == REF[drug]))
        fit = sub["Fitness"].to_numpy()
        thr = np.percentile(fit, 99.0)                    # top 1%: Fitness >= 99th pct
        top = sub.filter(pl.col("Fitness") >= thr)["Genotype"].to_list()
        n = len(top)
        rows = []
        for pos_i, amb in enumerate(AMBLER):
            col = [g[pos_i] for g in top]
            counts = {}
            for aa in col:
                counts[aa] = counts.get(aa, 0) + 1
            # cell-19 colour rule: threshold by number of residues present at the position;
            # non-WT above threshold = "enriched" (red); WT = black; else grey.
            npres = len(counts)
            ct = 0.55 if npres == 2 else 0.37 if npres == 3 else 0.275
            for aa in sorted(counts, key=lambda a: -counts[a]):
                freq = counts[aa] / n
                color = "black" if aa == WT_SEQ[pos_i] else ("red" if freq > ct else "grey")
                rows.append([amb, WT_SEQ[pos_i], aa, counts[aa], round(freq, 5),
                             "yes" if aa == WT_SEQ[pos_i] else "no", color])
        write_csv(f"fig3B_top1pct_logo_freq_{drug.lower()}.csv",
                  ["ambler_position", "WT_residue", "residue", "count", "frequency", "is_WT",
                   "logo_color"], rows)
        print(f"    {drug}: top-1% n={n} genotypes (threshold Fitness>={thr:.3f})")


def build_fig3cd(df):
    print("Fig 3C/D (fitness vs number of mutations):")
    per = {}
    for drug in ("AMP", "AZT"):
        sub = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == REF[drug]))
        med = (sub.group_by("Epistatic Order").agg(pl.col("Fitness").median().alias("median"),
                                                    pl.len().alias("n")).sort("Epistatic Order"))
        per[drug] = {r["Epistatic Order"]: (r["median"], r["n"]) for r in med.to_dicts()}
    orders = sorted(set(per["AMP"]) | set(per["AZT"]))
    mrows = [[o, per["AMP"].get(o, ("", ""))[1], round(per["AMP"].get(o, ("", ""))[0], 4) if o in per["AMP"] else "",
              per["AZT"].get(o, ("", ""))[1], round(per["AZT"].get(o, ("", ""))[0], 4) if o in per["AZT"] else ""]
             for o in orders]
    write_csv("fig3CD_median_by_nmut.csv",
              ["n_mutations", "AMP_n", "AMP_median_fitness", "AZT_n", "AZT_median_fitness"], mrows)
    # per-genotype companion
    amp = df.filter((pl.col("Drug") == "AMP") & (pl.col("Concentration") == 781.0))
    azt = df.filter((pl.col("Drug") == "AZT") & (pl.col("Concentration") == 36.0))
    afit = dict(zip(amp["Genotype"].to_list(), amp["Fitness"].to_list()))
    zfit = dict(zip(azt["Genotype"].to_list(), azt["Fitness"].to_list()))
    nmut = dict(zip(amp["Genotype"].to_list(), amp["Epistatic Order"].to_list()))
    grows = [[g, nmut[g], round(afit[g], 4), round(zfit.get(g, float("nan")), 4)] for g in amp["Genotype"].to_list()]
    write_csv("fig3CD_all_genotypes.csv",
              ["genotype", "n_mutations", "AMP781_fitness", "AZT36_fitness"], grows)


def build_fig4ab():
    print("Fig 4A/B (single-mutant fitness, n=3 replicates):")
    rows = []
    for drug in ("AMP", "AZT"):
        lg = pl.read_parquet(LONG[drug])
        sub = lg.filter(pl.col("concentration") == REF[drug])
        for r in sub.to_dicts():
            prof = r["mutant_profile"]
            nonwt = [i for i in range(13) if prof[i] != "."]
            if prof == "." * 13:
                label = "WT"
            elif len(nonwt) == 1 and prof[nonwt[0]] != "X":
                i = nonwt[0]
                label = f"{WT_SEQ[i]}{AMBLER[i]}{prof[i]}"
            else:
                continue
            rows.append([label, drug, round(r["replicate1"], 4), round(r["replicate2"], 4),
                         round(r["replicate3"], 4), round(r["mean"], 4), round(r["std"], 4)])
    write_csv("fig4AB_single_mutant_replicates.csv",
              ["mutation", "drug", "replicate1", "replicate2", "replicate3", "mean", "sd"], rows)


def build_fig4cd(df):
    print("Fig 4C/D (single+double-mutant AUC-Fitness matrix):")
    rows = []
    for drug in ("AMP", "AZT"):
        sub = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == REF[drug])
                        & (pl.col("Epistatic Order") <= 2))
        for g, f in zip(sub["Genotype"].to_list(), sub["Fitness"].to_list()):
            lab = mut_label(g)
            parts = lab.split("; ")
            if lab == "WT":
                m1, m2 = "none", "none"
            elif len(parts) == 1:
                m1, m2 = parts[0], "none"
            else:
                m1, m2 = parts[0], parts[1]
            rows.append([drug, m1, m2, round(f, 4)])
    write_csv("fig4CD_double_mutant_fitness.csv", ["drug", "mutation_1", "mutation_2", "AUC_fitness"], rows)


def main():
    df = pl.read_parquet(PARQUET)
    print(f"loaded {PARQUET.name}: {df.height} rows\n")
    build_fig3a(df)
    build_fig3b(df)
    build_fig3cd(df)
    build_fig4ab()
    build_fig4cd(df)
    print("\ndone.")


if __name__ == "__main__":
    main()
