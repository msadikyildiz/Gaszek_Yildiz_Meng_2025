"""Reproduce the Supplementary Figure S4 fitness-landscape graphs and export their
Source Data. Same coarse-graining as Figure 2 (Devin Meng's builder) but for the
full concentration series shown in S4 -- ampicillin 0 -> 781 and aztreonam 0 -> 324
ug/mL. The Gephi layout (x, y positions) is a manual, non-reproducible aesthetic
choice; the plotted quantities are the node fitness / size / peak status and the
edge weight / count, which this script writes deterministically.

Run:  python extract_figS4_source_data.py   (from this directory)
Out:  ./source_data/figS4_nodes.csv, ./source_data/figS4_edges.csv
"""
import os
import sys

# Pin hash seed so supernode representative tie-breaks are reproducible (mirrors
# extract_fig2_source_data.py); re-exec once if not already pinned.
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract_fig2_source_data import (  # noqa: E402  (shares the exact Fig 2 machinery)
    FIG2, NEUTRAL, PARAMS, load_long, neighbor_pairs, per_conc_pairs,
)
from src.graph_builder import GraphBuilder  # noqa: E402


def main():
    out = os.path.join(HERE, "source_data")
    os.makedirs(out, exist_ok=True)
    nodes, edges = [], []
    for drug, (csv, concs, _plotted) in FIG2.items():
        long = load_long(csv, concs)
        pdf = per_conc_pairs(neighbor_pairs(long["mutant_profile"].unique().to_list()), long)
        all_concs = [float(c.split(" ")[1]) for c in concs]  # every measured concentration (S4 shows all)
        for c in all_concs:
            gb = GraphBuilder(long_table=long, pairwise_table=pdf)
            g = gb.build_graph(concentration=c, neutral_threshold=NEUTRAL, fitness_col="median",
                               fitness_diff_col="median_diff", use_parallel_peak_merge=False,
                               rename_by_avg_fitness=True, merge_peaks=False, **PARAMS)
            for n, d in g.nodes(data=True):
                gm = d.get("group_mutants") or {n: d.get("fitness")}
                vals = list(gm.values())
                node_fit = max(vals) if d.get("is_peak") else sum(vals) / len(vals)
                nodes.append({"drug": drug, "concentration": c, "node_id": n, "fitness": node_fit,
                              "n_genotypes": d.get("group_size"), "is_peak": d.get("is_peak"),
                              "contains_wildtype": d.get("contain_wildtype")})
            for u, v, d in g.edges(data=True):
                edges.append({"drug": drug, "concentration": c, "source": u, "target": v,
                              "weight": d.get("weight"), "count": d.get("count")})
            print(f"  [{drug} {c}] {g.number_of_nodes()} nodes, {g.number_of_edges()} edges", flush=True)
    pl.DataFrame(nodes).sort("drug", "concentration", "node_id").write_csv(os.path.join(out, "figS4_nodes.csv"))
    pl.DataFrame(edges).sort("drug", "concentration", "source", "target").write_csv(os.path.join(out, "figS4_edges.csv"))
    print(f"WROTE {len(nodes)} nodes + {len(edges)} edges to {out}", flush=True)


if __name__ == "__main__":
    main()
