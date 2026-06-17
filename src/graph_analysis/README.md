# Graph-theoretic landscape and distribution analysis

The fitness-landscape graph and the distribution/summary-statistics analyses
behind several main-text and supplementary panels (the Figure 1 fitness
ridgelines, the Figure 2 landscape graphs, the Supplementary Figure S4
concentration series, and the Supplementary Figure S18 high-fitness tails). All
notebooks run directly from this repository's canonical processed data via
`repro_helpers.py`.

## Files

- `summary-stats.ipynb` — Polars/Matplotlib notebook, runs top to bottom:
  - Plot1/Plot2 — pairwise condition comparison + per-concentration fitness histograms
  - **Plot3 — joyplot of fitness across concentrations (Figure 1 ridgelines)**
  - density plots for global fitness; Plot4 — fitness with/without a given mutation
  - Plot5/Plot6 + heatmaps — amino-acid-state fitness, global-fitness heatmaps
  - Plot7/Plot8/Plot9 — global fitness and epistasis vs mutation count and vs fitness
  - replicate QC (per-genotype coefficient of variation) and growth-trajectory traces
- `finalize-2025-05-04.ipynb` — companion notebook (replicate-CV QC, the
  `sd-vs-mean` figures).
- `repro_helpers.py` — repository-backed data access used by both notebooks (see below).
- `reference_outputs/` — reference PNGs (joyplots, high-fitness AZT tails,
  sd-vs-mean) kept for side-by-side comparison against regenerated figures.

## How to run

```bash
uv run jupyter lab     # then run src/graph_analysis/summary-stats.ipynb
```

Both notebooks are verified end-to-end headless with the repository's `uv`
environment (`uv run jupyter nbconvert --to notebook --execute ...`). Figures the
notebooks save land in `./_generated/` (gitignored). Two dependencies beyond the
base pipeline: `joypy` (in `pyproject.toml`) and, optionally, `plotly` (only the
commented-out interactive cells; the import is guarded).

## repro_helpers.py

- `preprocess_data(...)` returns the per-drug long/wide AUC tables from
  `data/processed/{amp,azt}_auc_long_df.parquet`. Path arguments are accepted for
  signature compatibility and ignored.
- `calculate_normalized_fitness(long_df)` returns per-genotype **global fitness**,
  defined as the local median fitness at the representative concentration used
  throughout the paper (AMP 781 µg/ml, AZT 36 µg/ml — `amp_select_conc` /
  `azt_select_conc` in `05_epistasis_figures.ipynb`). At that concentration the
  local median is identical to the `Fitness` column of `Epistasis_Combined`
  (verified to 0.00e+00 across all 55,293 shared genotypes), and unlike the
  epistasis slice it retains the dead-mutant baseline several plots reference.
- `load_global_epistasis()` slices the per-drug global epistasis tables from
  `Epistasis_Combined.parquet` at that same concentration, matching how
  `05_epistasis_figures.ipynb` builds `Epistasis_Combined_{AMP,AZT}_auc_10`.
- On import it patches a `joypy` 0.2.6 / modern-pandas incompatibility
  (`flatten_axes` now returns a generator) so the ridgeline joyplots render.

## Fidelity note

Regenerated distributions match the `reference_outputs/` joyplots. The WT (black)
and dead-mutant (blue) annotation lines use this repository's canonical
per-concentration medians (e.g. AZT WT median at 36 µg/ml = 2.211, equal to the
published `Fitness`); the underlying densities are unchanged.

## Figure 2 (fitness-landscape graphs)

Figure 2 has no plotting-library code: the graph data is exported as GraphML and
laid out in **Gephi**. The settings:

```
Illustration graph (Gephi preview):
- node opacity 90
- node labels: proportional font size on; Apple Braille 4 Plain; no box
- edges: thickness 1; rescale weight min 3.0 max 10.0
- edge opacity 80 (set to 50 for high-AZT concentration, where edges are dense)
- edge arrow size 3; edge labels Apple Braille 20 Plain

Landscape graph (Gephi):
- node color = is_peak
- node size 15–50, exponential spline
- edge color deterministic
- global virtual-node fitness = 7.5; concentration-specific virtual-node fitness = 6
- layout ForceAtlas 2: scaling 150.0, prevent overlap
- rotate so the higher peak aligns to the right
- Network Splitter 3D: Z-Maximum 30

z-level per node = round( z_V * z_levels / z_max ),
  where z_V is the node's z value, z_levels the user-defined max number of
  levels, and z_max the largest value in the z column.
```

Node-size scaling in panels B and C (the six concentrations) was set so that
super-nodes carrying thousands of mutants remain visible.
