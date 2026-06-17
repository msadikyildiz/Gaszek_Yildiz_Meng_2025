# Revision analyses

Code, result tables and figures for the new analyses added during the Nature
Communications revision. Each sub-folder is self-contained and regenerates its
outputs from the deposited data in this repository
(`data/processed/Epistasis_Combined.parquet` and `data/raw/`); no external paths
are required.

## Contents

| Folder | Produces | Addresses |
|---|---|---|
| `s4_epistasis_rewrite/` | Supplementary Figs **S9** (epistatic-order decomposition) and **S10** (drug-asymmetry at matched stringency) | R3 (epistasis rewrite) |
| `s7_concentration_grid/` | Supplementary Figs **S11** (pairwise-epistasis heatmaps) and **S12** (measured-vs-predicted densities) across all concentrations; `library_invariance/` for the low-read robustness check | R1, R3 |
| `s1_model_comparison/` | Supplementary Fig **S13** (linear / Lasso / decision-tree / LightGBM learning curves) | R3-2b |
| `s2_block_holdout/` | Supplementary Fig **S14** (Leave-One-Mutation-Out + Hamming-distance-stratified holdout); `classification_check.py` for the LOMO classification metrics | R4-ii |
| `s3_rmsd_justification/` | Supplementary Fig **S15** (sensitivity / specificity / PPV / AUROC at three stringency thresholds; RMSD grounding) | R3, R4 |
| `s5_replicate_correlations/` | Supplementary Fig **S16** (per-genotype %CV and replicate-pair correlations) | R1 |
| `s6_fig_s3_revision/` | Supplementary Fig **S17** (revised dose-response profiles with low-read flagging) | R1 |

Supplementary Fig **S18** (peak-absorption diagnostic and neutrality-cutoff
robustness) is produced by the graph-analysis code in `../src/graph_analysis/`.

## Running

From the repository root, with the environment set up (`uv sync`):

```bash
uv run python revision_analyses/s1_model_comparison/analysis.py
# ... and likewise for each sub-folder
```

Each `analysis.py` locates the repository root automatically (it walks up to the
directory containing `data/processed/Epistasis_Combined.parquet`) and writes its
figures and result tables back into its own `figures/` and `data/` sub-folders.
Per-folder `README.md` and `results.md` files describe the analysis and report
the headline numbers.
