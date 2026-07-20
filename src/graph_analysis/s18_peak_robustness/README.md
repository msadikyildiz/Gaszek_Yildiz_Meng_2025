# Supplementary Figure S18 — global-peak robustness (graph analysis)

Code behind **Supplementary Figure S18**, contributed by **Devin Meng**. Delivered as-is
and committed verbatim for provenance; `README_original_meng.md` is Meng's own note.

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

> Fig 2 output is expected to be unchanged under this version (identical algorithm and
> parameters), but has **not yet been re-verified** against the committed
> `fig2_nodes.csv` / `fig2_edges.csv`. Unifying the two graph packages onto this newer
> version (and re-verifying Fig 2) is a deliberate follow-up, not done here.

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

## Status — folded in as received; not yet runnable in-repo

To run it here and extract the S18 Source Data (A/B fitness-advantage values; C
neutral-threshold matrix), two ports are needed:

1. **Data path/format** — point `preprocess_data(amp_path=…, azt_path=…)` at this repo's
   `data/raw/` AUC CSVs (or emit the `genotype_auc_sorted_*.csv` form the code expects).
2. **Output paths** — both notebooks hardcode `output_dir =
   Path("/work/greencenter/s439821/…")`; repoint to a repo-relative
   `outputs/global-peak-robustness/` (gitignored), and fix `generate_graphs.sh`'s
   `PROJECT_ROOT` for this directory depth.

Once ported: `generate_graphs.sh` → the two notebooks → the S18 Source Data sheet can be
populated (it is a stub in `source_data/`).
