"""Reproduce the Figure 2 fitness-landscape graphs and export their Source Data.

Runs the graph builder (author: Devin Meng; see ./src/) on the repository's
design-filtered per-genotype AUC data and writes node/edge tables for the
plotted concentrations. Deterministic: the builder uses no RNG, so the output
matches the published Figure 2. Devin's original ``get_pairs`` does an O(n)
DataFrame filter per pair (cluster-scale slow); the neighbour-diff computation
here is an identical dict-lookup equivalent.

Run:  python extract_fig2_source_data.py   (from this directory)
Out:  ./source_data/fig2_nodes.csv, ./source_data/fig2_edges.csv
"""
import os, sys
# Fix hash-seed so set/dict iteration (hence supernode representative tie-breaks) is
# reproducible across process invocations; re-exec once if not already pinned.
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)
import time
from itertools import product
from collections import defaultdict
import numpy as np
import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                                   # so `import src.*` resolves
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))  # graph_builder -> graph_analysis -> src -> repo
from src.preprocess import clean_nulls, convert_long_format   # noqa: E402
from src.graph_builder import GraphBuilder                     # noqa: E402

AMP_CSV = os.path.join(REPO, "data", "raw", "Ampicillin_auc_per_genotype.csv")
AZT_CSV = os.path.join(REPO, "data", "raw", "Aztreonam_auc_per_genotype.csv")
AMP_CONCS = ["Ampicillin 0.0", "Ampicillin 3.1", "Ampicillin 12.2", "Ampicillin 48.8", "Ampicillin 195.0", "Ampicillin 781.0"]
AZT_CONCS = ["Aztreonam 0.0", "Aztreonam 0.44", "Aztreonam 1.33", "Aztreonam 4.0", "Aztreonam 12.0", "Aztreonam 36.0", "Aztreonam 108.0", "Aztreonam 324.0"]
PARAMS = dict(tiny_initial_threshold=0.02, large_edge_threshold=5.5, num_forbidden_pairs=1)
NEUTRAL = 0.4
FIG2 = {"AMP": (AMP_CSV, AMP_CONCS, [12.2, 48.8, 195.0]),
        "AZT": (AZT_CSV, AZT_CONCS, [12.0, 36.0, 108.0])}

_intended = {19: ['.', 'P'], 37: ['.', 'K'], 67: ['.', 'L', 'V'], 102: ['.', 'K'], 162: ['.', 'S', 'H', 'N'],
             180: ['.', 'T'], 235: ['.', 'T'], 236: ['.', 'S'], 237: ['.', 'K'], 241: ['.', 'S', 'C'],
             261: ['.', 'M'], 271: ['.', 'L', 'Q'], 272: ['.', 'D']}
POSSIBLE = ["".join(m) for m in product(*[_intended[p] for p in sorted(_intended)])] + [".............", "XXXXXXXXXXXXX"]


def load_long(csv, concs):
    df = pl.read_csv(csv).filter(pl.col("mut_profile_masked").is_in(POSSIBLE))
    df = clean_nulls(df, num_rep_threshold=2)
    long = convert_long_format(df, concs)
    return long.with_columns(
        pl.struct(['replicate1', 'replicate2', 'replicate3'])
        .map_elements(lambda s: np.nanmedian([s["replicate1"], s["replicate2"], s["replicate3"]]),
                      return_dtype=pl.Float64).alias('median'))


def neighbor_pairs(profiles):
    d = defaultdict(list)
    for i, m in enumerate(profiles):
        for j in range(len(m)):
            d[m[:j] + '?' + m[j + 1:]].append(i)
    pairs = set()
    for idxs in d.values():
        for a in idxs:
            for b in idxs:
                if a < b and sum(c1 != c2 for c1, c2 in zip(profiles[a], profiles[b])) == 1:
                    pairs.add((profiles[a], profiles[b]))
    return pairs


def per_conc_pairs(neighbors, long):
    med = {(r['mutant_profile'], r['concentration']): r['median'] for r in long.iter_rows(named=True)}
    concs = long['concentration'].unique().to_list()
    return pl.DataFrame([{'mutant_profile1': m1, 'mutant_profile2': m2, 'concentration': c,
                          'median_diff': med.get((m1, c), np.nan) - med.get((m2, c), np.nan)}
                         for (m1, m2) in neighbors for c in concs])


def main():
    out = os.path.join(HERE, "source_data")
    os.makedirs(out, exist_ok=True)
    nodes, edges = [], []
    for drug, (csv, concs, plotted) in FIG2.items():
        t0 = time.time()
        long = load_long(csv, concs)
        pdf = per_conc_pairs(neighbor_pairs(long['mutant_profile'].unique().to_list()), long)
        print(f"[{drug}] {long.shape}, {time.time()-t0:.1f}s", flush=True)
        for c in plotted:
            gb = GraphBuilder(long_table=long, pairwise_table=pdf)
            # serial peak-merge => deterministic labels (parallel merge is tie-break-nondeterministic);
            # same partition either way.
            g = gb.build_graph(concentration=c, neutral_threshold=NEUTRAL, fitness_col='median',
                               fitness_diff_col='median_diff', use_parallel_peak_merge=False,
                               rename_by_avg_fitness=True, merge_peaks=False, **PARAMS)
            for n, d in g.nodes(data=True):
                gm = d.get('group_mutants') or {n: d.get('fitness')}
                vals = list(gm.values())
                # y-position per the Fig 2 caption: peak node = max member fitness; connection node = mean
                node_fit = max(vals) if d.get('is_peak') else sum(vals) / len(vals)
                nodes.append({'drug': drug, 'concentration': c, 'node_id': n, 'fitness': node_fit,
                              'n_genotypes': d.get('group_size'), 'is_peak': d.get('is_peak'),
                              'contains_wildtype': d.get('contain_wildtype')})
            for u, v, d in g.edges(data=True):
                edges.append({'drug': drug, 'concentration': c, 'source': u, 'target': v,
                              'weight': d.get('weight'), 'count': d.get('count')})
            print(f"  [{drug} {c}] {g.number_of_nodes()} nodes, {g.number_of_edges()} edges", flush=True)
    pl.DataFrame(nodes).sort("drug", "concentration", "node_id").write_csv(os.path.join(out, "fig2_nodes.csv"))
    pl.DataFrame(edges).sort("drug", "concentration", "source", "target").write_csv(os.path.join(out, "fig2_edges.csv"))
    print(f"WROTE {len(nodes)} nodes + {len(edges)} edges to {out}", flush=True)


if __name__ == "__main__":
    main()
