# Source Data

Code that regenerates the Nature Communications "Source Data" workbook: one
sheet per display item (main and supplementary figures), each holding the
underlying numerical values for that panel. This mirrors the journal's Source
Data requirement — reviewers and readers get the underlying table for every
figure, not just the rendered image. Two sheets are deterministic
re-computations rather than the original plotted arrays — Fig 5 (SHAP, re-fit
at seed 42; per-mutation mean|SHAP| Spearman rho 0.95-0.99 vs the original
unseeded run) and Fig 2 (landscape graph, re-derived from the design-filtered
data) — as each sheet's header note states.

## Layout

```
source_data/
  README.md          this file
  scripts/            builder scripts (assembler + 4 parquet builders + packager)
  derived/            per-figure CSVs written by the parquet builders (gitignored; regenerable)
  source_data.xlsx    assembled workbook (gitignored; regenerable)
  Source_Data.zip     release archive: workbook + all bulk CSVs (gitignored; regenerable)
```

`derived/`, `source_data.xlsx`, and `Source_Data.zip` are not committed — all
are regenerated deterministically from data already in this repository (see
"Reproducing" below). The Fig 6 / S6 and S18 CSVs live under their own analysis
folders (`src/evolutionary_statistics/source_data/`,
`src/graph_analysis/s18_peak_robustness/source_data/`); the packager gathers
them together. The scripts and this README are tracked.

## Scripts

| Script | Produces | Reads |
|---|---|---|
| `build_parquet_source_data.py` | Fig 1E (ridgeline distribution summary + histograms), Fig 4E-F (measured vs. predicted fit summary), Fig S12 (R² by order across the concentration grid) | `data/processed/Epistasis_Combined.parquet`, `*_auc_long_df.parquet`, `analysis/s11_s12_concentration_grid/data/regression_r2_by_order.csv` (cross-check) |
| `build_fig3_fig4_source_data.py` | Fig 3A-D (landscape signatures: overlay points, sequence logos, fitness-vs-order), Fig 4A-B (single-mutant replicates), Fig 4C-D (single/double-mutant fitness matrix) | `Epistasis_Combined.parquet`, `*_auc_long_df.parquet`, `data/known_variants/encoded_variants.csv` |
| `build_fig5_source_data.py` | Fig 5C (mean\|SHAP\| per mutation) + seed-robustness table | `Epistasis_Combined.parquet` (refits the published LightGBM + SHAP pipeline deterministically, seed 42; also verifies the attributions are not an artifact of the original notebooks' unseeded split) |
| `build_s7s8_source_data.py` | Fig S7 (monoculture IC50 vs. mean AUC-fitness), Fig S8 (per-concentration AUC-fitness vs. IC50) | `validation/src/processed/{20260407,20260307}/xref_expanded_df.parquet` |
| `build_source_data.py` | Assembles `source_data.xlsx`: one sheet per display item, pulling the CSVs above plus the pre-computed tables in `analysis/s09_s10_epistatic_order/` through `analysis/s17_dose_response_low_count/`, Fig 6/S6, and S18 | Everything above, via `derived/`, `analysis/*/data/*.csv`, and the Fig 6/S18 module `source_data/` folders |
| `package_source_data.py` | `Source_Data.zip` = `source_data.xlsx` + `tables/` (all bulk CSVs gathered from `derived/`, Fig 6/S6, and S18) | `source_data.xlsx` + the three CSV folders |

Each script is self-contained: it walks up from its own location to find the
repository root (the directory containing `data/processed/Epistasis_Combined.parquet`)
and resolves every input relative to that root — no machine-specific paths.

## Reproducing

From the repository root, with the environment set up (`uv sync`):

```bash
# 1. Per-figure CSVs from the processed parquets -> source_data/derived/
uv run python source_data/scripts/build_parquet_source_data.py
uv run python source_data/scripts/build_fig3_fig4_source_data.py
uv run python source_data/scripts/build_fig5_source_data.py
uv run python source_data/scripts/build_s7s8_source_data.py

# 2. Fig 2 landscape tables (separate generator; copy into derived/)
PYTHONHASHSEED=0 uv run python src/graph_analysis/graph_builder/extract_fig2_source_data.py
cp src/graph_analysis/graph_builder/source_data/fig2_{nodes,edges}.csv source_data/derived/

# 3. Fig 6 / S6 (DCA) and S18 (graph) source data — see their module READMEs
uv run --extra evo-stats python src/evolutionary_statistics/reproduce_fig6.py
uv run python src/graph_analysis/s18_peak_robustness/reproduce_s18.py   # committed CSVs suffice; this refreshes full/

# 4. Assemble the workbook, then package the release archive
uv run python source_data/scripts/build_source_data.py
uv run python source_data/scripts/package_source_data.py
```

`build_fig5_source_data.py`'s printed seed-robustness check and
`build_s7s8_source_data.py`'s printed Pearson r/n match the values documented in
`build_source_data.py`'s per-sheet notes (S7: AMP r=0.88 n=13, AZT r=0.80 n=18).

## Status — complete

`source_data.xlsx` has **20 of 20 data sheets populated**; the packaged
`Source_Data.zip` (workbook + 40 bulk CSVs under `tables/`) is the release
archive. The only non-data item is **Fig 1F** (structural render + ChemDraw
artwork — no plotted data), which the workbook marks N/A. Fig 6 (DCA panels
A/B/C/G/H + S6) and S18 are reproduced in-repo; the LGL panels 6D-F come from the
separate `github.com/morcoslab/LGL-VAE` (deposited coordinates) and are not part
of these tables. Both `source_data.xlsx` and `Source_Data.zip` are gitignored and
rebuilt by the commands above.

## Provenance

These 5 scripts were consolidated into this repository's `analysis/` layout
from the lab's internal development workspace, where they were originally
organized under a different, per-topic numbering. This repo renumbers those
as `analysis/s09_s10_epistatic_order/` …
`analysis/s17_dose_response_low_count/` (see `analysis/README.md`); the path
rewrites in these scripts track that renumbering one-to-one. No analysis
logic changed in the port — only path resolution.
