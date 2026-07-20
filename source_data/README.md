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
  scripts/            5 builder scripts (ported from the lab's internal revision workspace)
  derived/            per-figure CSVs written by 4 of the 5 scripts (gitignored; regenerable)
```

`derived/` is not committed. It is scratch output — every file in it is
regenerated deterministically from data already in this repository (see
"Reproducing" below). Nothing under `source_data/` other than `derived/`
is gitignored: the scripts and this README are tracked.

## Scripts

| Script | Produces | Reads |
|---|---|---|
| `build_parquet_source_data.py` | Fig 1E (ridgeline distribution summary + histograms), Fig 4E-F (measured vs. predicted fit summary), Fig S12 (R² by order across the concentration grid) | `data/processed/Epistasis_Combined.parquet`, `*_auc_long_df.parquet`, `analysis/s11_s12_concentration_grid/data/regression_r2_by_order.csv` (cross-check) |
| `build_fig3_fig4_source_data.py` | Fig 3A-D (landscape signatures: overlay points, sequence logos, fitness-vs-order), Fig 4A-B (single-mutant replicates), Fig 4C-D (single/double-mutant fitness matrix) | `Epistasis_Combined.parquet`, `*_auc_long_df.parquet`, `data/known_variants/encoded_variants.csv` |
| `build_fig5_source_data.py` | Fig 5C (mean\|SHAP\| per mutation) + seed-robustness table | `Epistasis_Combined.parquet` (refits the published LightGBM + SHAP pipeline deterministically, seed 42; also verifies the attributions are not an artifact of the original notebooks' unseeded split) |
| `build_s7s8_source_data.py` | Fig S7 (monoculture IC50 vs. mean AUC-fitness), Fig S8 (per-concentration AUC-fitness vs. IC50) | `validation/src/processed/{20260407,20260307}/xref_expanded_df.parquet` |
| `build_source_data.py` | Assembles `source_data.xlsx`: one sheet per display item, pulling the CSVs above plus the pre-computed tables in `analysis/s09_s10_epistatic_order/` through `analysis/s17_dose_response_low_count/` | Everything above, via `derived/` and `analysis/*/data/*.csv` |

Each script is self-contained: it walks up from its own location to find the
repository root (the directory containing `data/processed/Epistasis_Combined.parquet`)
and resolves every input relative to that root — no machine-specific paths.

## Reproducing

From the repository root, with the environment set up (`uv sync`):

```bash
# 1. Generate the per-figure CSVs into source_data/derived/
uv run python source_data/scripts/build_parquet_source_data.py
uv run python source_data/scripts/build_fig3_fig4_source_data.py
uv run python source_data/scripts/build_fig5_source_data.py
uv run python source_data/scripts/build_s7s8_source_data.py

# 2. Assemble the workbook
uv run python source_data/scripts/build_source_data.py
```

Step 1 has been re-run end to end in this repository as part of the port; all
four scripts complete cleanly, and `build_fig5_source_data.py`'s printed
seed-robustness check and `build_s7s8_source_data.py`'s printed Pearson r/n
match the values already documented in `build_source_data.py`'s per-sheet
notes (S7: AMP r=0.88 n=13, AZT r=0.80 n=18).

**Fig 2 is the one sheet these scripts do not populate.** Its node/edge
tables come from a separate, already-integrated generator,
`src/graph_analysis/graph_builder/extract_fig2_source_data.py`, which writes
`fig2_nodes.csv` / `fig2_edges.csv` into its own local `source_data/` folder
(not this one). To populate the "Fig2 landscape graphs" sheet, run that
script and copy (or symlink) its two output CSVs into `source_data/derived/`
before running the assembler.

## What is NOT done yet

`source_data.xlsx` is **not built or committed** by this port. Two sheets are
release-gated on outstanding contributor work:

- **Fig 6** (DCA logo, effective alphabet, Hamiltonian distributions, LGL
  coordinates) — pending F.M./Sophia/Alberto: we have requested the reproducible
  code (preferred — it also settles the Fig 6A/6B legend items) or the source
  values, with "available on request" only as a fallback if they decline.
- **S18** (peak-advantage box plots + neutral-threshold matrix) — Devin/Milo.

`build_source_data.py` stubs both sheets (plus Fig 1F, which has no plotted
data — it is a structural render / ChemDraw artwork) with a clear
"TO FILL (owner): ..." note rather than failing, so it can be run early to
sanity-check the rest of the workbook. The real assembly — and the Zenodo
release it feeds — waits for those two sheets.

## Provenance

These 5 scripts were ported from the lab's private revision workspace
(`rebuttal_response/editorial_revision/source_data/`), where they were
originally written against that workspace's `sprints/sN_*` layout. This repo
renumbers those as `analysis/s09_s10_epistatic_order/` …
`analysis/s17_dose_response_low_count/` (see `analysis/README.md`); the path
rewrites in these scripts track that renumbering one-to-one. No analysis
logic changed in the port — only path resolution.
