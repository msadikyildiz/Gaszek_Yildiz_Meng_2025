"""Reproduce Supplementary Figure S18 (global-peak robustness) end-to-end from this
repository's AUC data, using Devin Meng's `fitness_landscape_graph` library unchanged.

Pipeline:
  1. Per-concentration pair tables (median_diff) from data/raw/ AUC. Uses a vectorized
     equivalent of fitness_landscape_graph.pair_table.get_mutant_pairs (Meng's version does
     ~4M per-pair-per-concentration polars filters); verified byte-equal to his function on
     genotype subsets.
  2. Neutral-threshold graph sweep (0.15-0.45 step 0.01 x 8 AZT concentrations) via
     `python -m fitness_landscape_graph.build_graphs_parallel` (~15-20 min, 248 graphs).
  3. S18C neutral-threshold robustness matrix: GraphAnalyzer.has_global_peak(min_group_size=12)
     over the sweep.
  4. S18A/B peak fitness advantage: GraphAnalyzer.get_peak_genotypes(rank=0) on the AZT 12 and
     AZT 108 main graphs (neutral_threshold=0.40) -> FitnessAdvantageAnalyzer(max_distance=2).

Outputs (source_data/): figS18C_neutral_threshold_matrix.csv,
  figS18{A,B}_azt{12,108}_peak_advantage_boxstats.csv, and the peak-group genotype lists.

Requires the repo env plus `plotly` (imported by graph_analyzer). Run from this directory:
    python reproduce_s18.py            # full run
    python reproduce_s18.py --skip-build   # reuse an existing graph sweep under --work
"""
import argparse
import subprocess
import sys
from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent


def _repo_root(p):
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


REPO = _repo_root(HERE)
sys.path.insert(0, str(HERE))
from fitness_landscape_graph.preprocess import clean_nulls, convert_long_format, preprocess_data
from fitness_landscape_graph.graph_analyzer import GraphAnalyzer
from fitness_landscape_graph.fitness_advantage import FitnessAdvantageAnalyzer

INTENDED = {19:[".","P"],37:[".","K"],67:[".","L","V"],102:[".","K"],162:[".","S","H","N"],
            180:[".","T"],235:[".","T"],236:[".","S"],237:[".","K"],241:[".","S","C"],
            261:[".","M"],271:[".","L","Q"],272:[".","D"]}
ALLOWED = list(np.append(np.array(["".join(m) for m in product(*[INTENDED[p] for p in sorted(INTENDED)])]),
                         ["."*13, "X"*13]))
AZT_CONCS = ["Aztreonam 0.0","Aztreonam 0.44","Aztreonam 1.33","Aztreonam 4.0","Aztreonam 12.0","Aztreonam 36.0","Aztreonam 108.0","Aztreonam 324.0"]
AMP_CONCS = ["Ampicillin 0.0","Ampicillin 3.1","Ampicillin 12.2","Ampicillin 48.8","Ampicillin 195.0","Ampicillin 781.0"]
SWEEP_CONCS = [0.0, 0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]
MIN_GROUP_SIZE = 12


def _make_long(path, concs):
    df = pl.read_csv(path).filter(pl.col("mut_profile_masked").is_in(ALLOWED))
    return convert_long_format(clean_nulls(df, num_rep_threshold=2), concs)


def _pairs_lean(long_df):
    """Per-concentration pairs (mutant_profile1/2, concentration, median_diff) — the columns
    GraphBuilder(fitness_diff_col='median_diff') and FitnessAdvantageAnalyzer consume.
    Vectorized equivalent of pair_table.get_mutant_pairs (verified byte-equal on subsets)."""
    mp = sorted(long_df["mutant_profile"].unique().to_list())
    cdd = defaultdict(set)
    for i, m in enumerate(mp):
        for j in range(len(m)):
            cdd[m[:j] + "?" + m[j+1:]].add(i)
    ps = set()
    for idxs in cdd.values():
        for i in idxs:
            for j in idxs:
                if i < j and sum(a != b for a, b in zip(mp[i], mp[j], strict=True)) == 1:
                    ps.add((mp[i], mp[j]))
    pairs = sorted(ps)
    m1s = [p[0] for p in pairs]; m2s = [p[1] for p in pairs]; n = len(pairs)
    frames = []
    for conc in sorted(long_df["concentration"].unique().to_list()):
        sub = long_df.filter(pl.col("concentration") == conc)
        med = {r["mutant_profile"]: np.nanmedian([r["replicate1"], r["replicate2"], r["replicate3"]])
               for r in sub.iter_rows(named=True)}
        md = np.array([med[m] for m in m1s]) - np.array([med[m] for m in m2s])
        frames.append(pl.DataFrame({"mutant_profile1": m1s, "mutant_profile2": m2s,
                                    "concentration": [float(conc)]*n, "median_diff": md}))
    return pl.concat(frames)


def _box_summary(res):
    return (res.group_by(["concentration", "distance"]).agg([
        pl.len().alias("n"),
        pl.col("fitness_diff").median().alias("median"),
        pl.col("fitness_diff").quantile(0.25).alias("q25"),
        pl.col("fitness_diff").quantile(0.75).alias("q75"),
        pl.col("fitness_diff").min().alias("min"),
        pl.col("fitness_diff").max().alias("max"),
        pl.col("fitness_diff").mean().alias("mean")]).sort(["concentration", "distance"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/s18_work", help="scratch dir for pairs + graphs")
    ap.add_argument("--skip-build", action="store_true", help="reuse an existing sweep under --work")
    args = ap.parse_args()
    work = Path(args.work)
    graphs = work / "outputs/global-peak-robustness"
    out = HERE / "source_data"; out.mkdir(exist_ok=True)

    if not args.skip_build:
        # 1. pairs + base-path layout Meng's build expects
        (work / "data/raw/combined-auc").mkdir(parents=True, exist_ok=True)
        (work / "data/processed").mkdir(parents=True, exist_ok=True)
        for drug, raw, concs in [("ampicillin", "Ampicillin", AMP_CONCS), ("aztreonam", "Aztreonam", AZT_CONCS)]:
            (work / f"data/raw/combined-auc/genotype_auc_sorted_{drug}.csv").unlink(missing_ok=True)
            (work / f"data/raw/combined-auc/genotype_auc_sorted_{drug}.csv").symlink_to(
                REPO / f"data/raw/{raw}_auc_per_genotype.csv")
            _pairs_lean(_make_long(REPO / f"data/raw/{raw}_auc_per_genotype.csv", concs)).write_csv(
                work / f"data/processed/{drug[:3]}_pairs.csv")
        print("pairs + layout ready; building 248-graph sweep (~15-20 min)...", flush=True)
        subprocess.run([sys.executable, "-m", "fitness_landscape_graph.build_graphs_parallel",
                        "--base-path", str(work), "--output-dir", str(graphs),
                        "--neutral-thresholds", "0.15", "0.45", "0.01",
                        "--concentrations", *[str(c) for c in SWEEP_CONCS]],
                       cwd=str(HERE), check=True)

    # 3. S18C matrix
    rows = []
    for gp in sorted(graphs.glob("azt_c*_t*.graphml")):
        cs, ts = gp.stem.replace("azt_c", "").split("_t")
        has_peak, _, _ = GraphAnalyzer(str(gp)).has_global_peak(min_group_size=MIN_GROUP_SIZE)
        rows.append({"concentration": float(cs.replace("_", ".")),
                     "neutral_threshold": float(ts.replace("_", ".")), "has_global_peak": bool(has_peak)})
    pl.DataFrame(rows).sort(["concentration", "neutral_threshold"]).write_csv(
        out / "figS18C_neutral_threshold_matrix.csv")
    print(f"S18C: {len(rows)} (conc x threshold) rows written", flush=True)

    # 4. S18A/B peak fitness advantage
    pdata = preprocess_data(
        amp_path=str(work / "data/raw/combined-auc/genotype_auc_sorted_ampicillin.csv"),
        azt_path=str(work / "data/raw/combined-auc/genotype_auc_sorted_aztreonam.csv"),
        amp_pairs_path=str(work / "data/processed/amp_pairs.csv"),
        azt_pairs_path=str(work / "data/processed/azt_pairs.csv"),
        clean_nulls_flag=True)
    fa = FitnessAdvantageAnalyzer(pdata["azt"]["pairs"], pdata["azt"]["long"])
    for panel, conc in [("A", 12.0), ("B", 108.0)]:
        ga = GraphAnalyzer(str(graphs / f"azt_c{str(conc).replace('.', '_')}_t0_40.graphml"))
        group = ga.get_peak_genotypes(rank=0)
        res = fa.compute_fitness_advantage(group_genotypes=group, max_distance=2)
        _box_summary(res).write_csv(out / f"figS18{panel}_azt{int(conc)}_peak_advantage_boxstats.csv")
        pl.DataFrame({"genotype": sorted(group)}).write_csv(
            out / f"figS18{panel}_azt{int(conc)}_peak_group_genotypes.csv")
        print(f"S18{panel}: AZT {int(conc)} peak group={len(group)}, {res.height} comparisons", flush=True)
    print("done -> source_data/", flush=True)


if __name__ == "__main__":
    main()
