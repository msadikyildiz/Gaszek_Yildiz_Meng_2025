# Supplementary Figure S18 — global-peak robustness (graph analysis)

Code behind **Supplementary Figure S18**, contributed by **Devin Meng**. See
`README_original_meng.md` for the original module notes.

## Panels

| Panel | Produced by | What it is |
|---|---|---|
| **S18A, S18B** | `fitness_advantage_analysis.ipynb` | fitness advantage of a peak supernode over its external neighbours |
| **S18C** | `neutral_threshold_robustness.ipynb` | global-peak presence across the neutral-threshold sweep (0.15–0.45) |
| **S18D** | Gephi | landscape-graph visualisation (same Gephi protocol as the main figures) |

Pipeline: `generate_graphs.sh` builds the aztreonam graphs across concentrations × neutral
thresholds (`python -m fitness_landscape_graph.build_graphs_parallel`, ~20 min), then the two
notebooks read those graphs.

## The `fitness_landscape_graph/` package

This is a **newer, fuller version of the Fig 2 graph code** in
`../graph_builder/`. The shared modules (`graph_builder.py`, `build_graph.py`,
`make_logo.py`, `pair_table_global.py`, `preprocess.py`) carry the **same core algorithm
and the same Fig 2 parameters** (`neutral_threshold=0.4`, `tiny_initial_threshold=0.02`,
`large_edge_threshold=5.5`, `num_forbidden_pairs=1`); the differences are refactoring —
a `--base-path` argument that removes the hardcoded cluster paths, the package rename
`src` → `fitness_landscape_graph`, and formatting. It **adds** the S18-specific modules
`fitness_advantage.py`, `graph_analyzer.py`, `build_graphs_parallel.py`, `pair_table.py`.

## Data

The pipeline reads **our** AUC-fitness data — no external Meng data is required:

- `pair_table_global.py` **generates** the global-fitness table and the `*_pairs.csv`
  pair tables itself (`calculate_normalized_fitness` + `get_pairs`); they are not inputs.
- `preprocess.py` reads the raw per-genotype AUC CSVs. Meng's code expects them as
  `genotype_auc_sorted_ampicillin.csv` / `..._aztreonam.csv` under `combined-auc/` — a
  sorted/renamed form of this repo's `data/raw/Ampicillin_auc_per_genotype.csv` /
  `Aztreonam_auc_per_genotype.csv`.

Dependencies (`polars`, `networkx`, `matplotlib`, `seaborn`, `logomaker`, `numpy`) are the
same set already used by the Fig 2 graph code — no new packages.

## Reproduced in-repo

`reproduce_s18.py` runs the whole pipeline on this repo's AUC data and writes the S18
Source Data into `source_data/`:

```bash
uv pip install plotly           # graph_analyzer's only extra dependency
python reproduce_s18.py         # ~15-20 min (the 248-graph neutral-threshold sweep)
python reproduce_s18.py --skip-build --work <dir>   # re-run only the analysis on an existing sweep
```

It (1) builds per-concentration pair tables from `data/raw/` via a vectorized equivalent of
`pair_table.get_mutant_pairs` (verified to produce identical `median_diff` and the same
Hamming-1 pair set as Meng's function on genotype subsets — his does ~4M
per-pair-per-concentration polars filters), (2) runs `fitness_landscape_graph.build_graphs_parallel`
(with `PYTHONHASHSEED=0` so supernode representative labels are deterministic) for the AZT
sweep, and (3) runs the S18C / S18A / S18B / S18D extraction with Meng's `GraphAnalyzer` +
`FitnessAdvantageAnalyzer` unchanged.

**Verified against the published `figure_s18.png`:** the S18C matrix reproduces the green/grey
pattern (conc 0/0.44/1.33 → no peak; 4 → 19/31; 36 → 1/31; 108 → 16/31 with the 0.42→0.43 gap
in panel D; 324 → 21/31), and the AZT-12 and AZT-108 peak fitness advantages peak at ~1.5 and
~2.35 at their respective concentrations, matching panels A and B.

Source Data (`source_data/`): `figS18C_neutral_threshold_matrix.csv`;
`figS18{A,B}_azt{12,108}_peak_advantage_boxstats.csv` (box statistics — the **full plotted
observations**, every group-member vs external-neighbour comparison including fliers, are the
regenerable `source_data/full/figS18{A,B}_*_all_comparisons.csv`, ~1.9M rows, gitignored and
destined for the release Source Data zip); `figS18D_azt108_{nodes,edges}.csv` (the AZT-108
t=0.42-vs-0.43 networks behind panel D; Gephi layout not reproducible); and the peak-group
genotype lists behind the panel A/B logos.

The two notebooks (`fitness_advantage_analysis.ipynb`, `neutral_threshold_robustness.ipynb`)
contain the same analysis and the figure-plotting code; `reproduce_s18.py` is the
repo-relative entry point. Panel D is laid out in Gephi (same protocol as the main
figures).
