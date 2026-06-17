# S5 — replicate reproducibility

Response to Reviewer #1 on manuscript line 163: quantify replicate-to-
replicate reproducibility of the AUC-Fitness metric, which the reviewer
flagged as unsupported by a correlation figure.

## Rerun

```
export MPLCONFIGDIR=/private/tmp/claude/mpl_cache FONTCONFIG_CACHE=/private/tmp/claude/fc_cache
mkdir -p $MPLCONFIGDIR $FONTCONFIG_CACHE
uv run python revision_analyses/s5_replicate_correlations/analysis.py
```

Runtime: ~30 s single-core. Polars + Matplotlib, reproducible (seed 20260420).

## Inputs

- `Gaszek_Yildiz_Meng_2025/data/raw/Ampicillin_auc_per_genotype.csv`
- `Gaszek_Yildiz_Meng_2025/data/raw/Aztreonam_auc_per_genotype.csv`
  Schema: `mut_profile_masked, <Drug> <conc> 1, <Drug> <conc> 2, <Drug> <conc> 3`.
  Cell values are already log10(AUC) (manuscript "fitness").

## Outputs

- `figures/replicate_scatter_amp.png` — 5x3 hex-bin grid of replicate-pair
  scatter plots for the five nonzero AMP concentrations. Identity line in
  drug colour; Pearson r, p-value, n annotated per panel.
- `figures/replicate_scatter_azt.png` — 7x3 hex-bin grid for the seven
  nonzero AZT concentrations.
- `figures/replicate_cv_amp.png` — per-concentration histograms of
  (left) %CV on raw AUC and (right) SD of log10(AUC). Dashed line at 10% CV
  marks the manuscript's claimed threshold; medians, IQR, and fraction of
  genotypes below 10% CV are annotated in-panel.
- `figures/replicate_cv_azt.png` — same layout for AZT.
- `data/per_genotype_replicate_stats.csv` — one row per (genotype, drug,
  concentration); columns: `genotype, drug, concentration, log_mean,
  log_sd, raw_mean, raw_sd, cv_percent` (985,674 rows).
- `data/summary_by_conc.csv` — 14 rows (drug x concentration including
  conc=0) with median/IQR/quantiles of %CV and log_sd, and fraction of
  genotypes with %CV < 10.
- `data/replicate_pair_correlations.csv` — 42 rows: one per unique
  replicate pair (1-2, 1-3, 2-3) per drug × concentration; columns:
  `drug, concentration, rep_i, rep_j, n, pearson_r, pearson_p`. The
  `pearson_p` column underflows to `0.0` at this sample size (double
  precision).
- `data/summary_by_conc_viable.csv` — 14 rows; same columns as
  `summary_by_conc.csv` but restricted to genotypes with mean
  log10(AUC) > 3.0 at the given condition. Confirms the "<10% SD"
  claim fails even in the viable subset.
- `results_table.csv` — 14-row headline verdict table: CV summary joined
  to mean Pearson r across pairs, with `median_cv_below_10pct` boolean.
- `results.md` — rebuttal response draft.

## Design choices

- **Conc=0 kept in tables, dropped from figures.** Zero-drug is a pure
  growth control and its dispersion does not speak to the AUC-fitness
  claim, but it is included in `summary_by_conc.csv` for completeness.
- **Two reproducibility metrics, reported side by side.**
  - %CV on raw AUC = 100 * SD(10^fitness) / mean(10^fitness). Matches the
    typical reader interpretation of "SD below 10%."
  - SD of log10(AUC) directly. Matches the manuscript units.
- **Three replicate pairs per concentration**, shown as 3 sub-panels rather
  than an average — the reviewer explicitly asked for scatter plots.
- **Hex-bin with log-count colour** scales cleanly across 65k points
  without over-plotting; identity line in drug colour.
- **Pearson is the right summary** (replicates should cluster around the
  identity line, not just "rank-correlate"). Spearman would over-state
  agreement in this regime.

## Status

- Figures and tables: complete, reproducible, correct.
- Manuscript text proposal: see `results.md`. Line 163 is replaced with
  honest, data-backed wording citing the Pearson r range and the
  log-SD dispersion directly. Denominators corrected, p-value wording
  fixed, and the viable-subset claim verified and saved to
  `data/summary_by_conc_viable.csv`.
