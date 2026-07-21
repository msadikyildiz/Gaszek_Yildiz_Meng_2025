# Fitness-landscape graph builder (Figure 2 + Supplementary Figure S18)

**Author of the graph-construction code (`src/`): Devin Meng.**

The `src/` package builds the coarse-grained fitness-landscape graphs behind
Figure 2 (per-concentration landscapes) and the Supplementary Figure S18
neutrality-cutoff analysis. Genotypes are nodes; single-mutation transitions are
edges oriented toward higher fitness; edges below the neutral threshold are
merged (union-find), then peak/connection nodes are identified.

## Files (`src/`, by Devin Meng)

- `preprocess.py` — raw genotype-AUC CSV → intended-mutant-filtered long table (per-concentration median).
- `pair_table_global.py` — global (integrated-AUC) fitness + Hamming-1 neighbour pairs.
- `graph_builder.py` — `GraphBuilder`: raw graph → neutral merging → peak detection (`build_graph`).
- `build_graph.py` — batch driver (writes GraphML).
- `make_logo.py` — sequence-logo / mutation helpers used by the builder.

Build parameters (from `build_graph.sh`): `neutral_threshold=0.4`,
`tiny_initial_threshold=0.02`, `large_edge_threshold=5.5`, `num_forbidden_pairs=1`.

## Reproducing the Figure 2 Source Data

`extract_fig2_source_data.py` runs `GraphBuilder` on this repository's
`data/raw/{Ampicillin,Aztreonam}_auc_per_genotype.csv` and writes node/edge
tables for the plotted concentrations (AMP 12.2/48.8/195, AZT 12/36/108 µg/mL):

```bash
python extract_fig2_source_data.py     # -> source_data/fig2_{nodes,edges}.csv
```

The construction is deterministic (no RNG), so the output matches the published
Figure 2. Each concentration partitions all 55,293 design genotypes across the
merged nodes; `is_peak = 1` marks a peak node (out-degree 0). These two CSVs are
the Figure 2 sheet in the paper's Source Data file. The driver's only departure
from `src/pair_table_global.get_pairs` is an O(1) dict lookup for neighbour
fitness differences in place of the original O(n) per-pair DataFrame filter — the
values are identical.

## Not included here

The Supplementary Figure S18 *peak-absorption diagnostic* (per-supernode fitness
advantage over 1-/2-mutation neighbours; global-peak presence across the
0.15–0.45 cutoff sweep) is a separate analysis on top of these graphs and is not
part of this builder — its Source Data is provided by Devin Meng.

## Batch drivers

The batch drivers `src/build_graph.py`, `src/pair_table_global.py`, and
`src/make_logo.py` are not used by the reproduction entry points
(`extract_fig2_source_data.py` for Figure 2, `extract_figS4_source_data.py` for
Supplementary Figure S4), which read only this repository's `data/`. Their `__main__` blocks default their base path to the repository root
(override with the `TEM1CML_ROOT` environment variable).
