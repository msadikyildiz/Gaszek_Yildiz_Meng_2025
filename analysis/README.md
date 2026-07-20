# Revision analyses

Code, result tables and figures for the new analyses added during the Nature
Communications revision. Each sub-folder is self-contained and regenerates its
outputs from the deposited data in this repository
(`data/processed/Epistasis_Combined.parquet` and `data/raw/`); no external paths
are required.

## Contents

| Folder | Produces |
|---|---|
| `s09_s10_epistatic_order/` | Supplementary Figs **S9** (epistatic-order decomposition) and **S10** (drug-asymmetry at matched stringency) |
| `s11_s12_concentration_grid/` | Supplementary Figs **S11** (pairwise-epistasis heatmaps) and **S12** (measured-vs-predicted densities) across all concentrations; `library_invariance/` for the low-read robustness check |
| `s13_model_comparison/` | Supplementary Fig **S13** (linear / Lasso / decision-tree / LightGBM learning curves) |
| `s14_mutation_holdout/` | Supplementary Fig **S14** (Leave-One-Mutation-Out + Hamming-distance-stratified holdout); `classification_check.py` for the LOMO classification metrics |
| `s15_classification_metrics/` | Supplementary Fig **S15** (sensitivity / specificity / PPV / AUROC at three stringency thresholds; RMSD grounding) |
| `s16_replicate_reproducibility/` | Supplementary Fig **S16** (per-genotype %CV and replicate-pair correlations) |
| `s17_dose_response_low_count/` | Supplementary Fig **S17** (revised dose-response profiles with low-read flagging) |

Supplementary Fig **S18** (peak-absorption diagnostic and neutrality-cutoff
robustness) is produced by the graph-analysis code in `../src/graph_analysis/`.

## Running

From the repository root, with the environment set up (`uv sync`):

```bash
uv run python analysis/s13_model_comparison/analysis.py
# ... and likewise for each sub-folder
```

Each `analysis.py` locates the repository root automatically (it walks up to the
directory containing `data/processed/Epistasis_Combined.parquet`) and writes its
published Supplementary panel(s) into the shared `figures/supplementary/` gallery,
and its result tables and secondary/diagnostic panels into its own `data/` and
`figures/` sub-folders (the latter gitignored as regenerable scratch).
Per-folder `README.md` files describe the analysis and report the headline
numbers.
