# s16 — replicate-to-replicate reproducibility of the AUC-fitness metric

Produces **Supplementary Figure S16**. Quantifies replicate-to-replicate
reproducibility of the AUC-fitness metric via per-genotype %CV and pairwise
replicate correlations, addressing the manuscript's reproducibility claim.

## Rerun

```
uv run python analysis/s16_replicate_reproducibility/analysis.py
```

Runtime: ~30 s single-core. Polars + Matplotlib, reproducible (seed 20260420).

## Inputs

- `data/raw/Ampicillin_auc_per_genotype.csv`
- `data/raw/Aztreonam_auc_per_genotype.csv`
  Schema: `mut_profile_masked, <Drug> <conc> 1, <Drug> <conc> 2, <Drug> <conc> 3`.
  Cell values are already log10(AUC) (manuscript "fitness").
- **Design filter.** Analysis is restricted to the 55,296 combinatorially-complete
  design genotypes (masked profiles present in `data/processed/Epistasis_Combined.parquet`);
  off-design barcode artefacts in the raw per-genotype table are excluded so the
  reproducibility statistics use the same population analysed in the paper. Per-concentration
  summaries report n = 55,293 — the three profiles without a finite replicate SD are dropped
  from the dispersion tables.

## Outputs

The replicate-scatter panels are the published Supplementary Fig. S16 and are
written directly to the shared gallery `figures/supplementary/`. The CV/SD
distribution panels (and the viable-subset variants) are secondary/diagnostic
and regenerate into this module's local `figures/` (gitignored, scratch).

- `figures/supplementary/figure_s16_amp.png` — 5x3 hex-bin grid of replicate-pair
  scatter plots for the five nonzero AMP concentrations. Identity line in
  drug colour; Pearson r, p-value, n annotated per panel. **Published
  Supplementary Fig. S16 (AMP)**.
- `figures/supplementary/figure_s16_azt.png` — 7x3 hex-bin grid for the seven
  nonzero AZT concentrations. **Published Supplementary Fig. S16 (AZT)**.
- `figures/replicate_cv_amp.png` — per-concentration histograms of
  (left) %CV on raw AUC and (right) SD of log10(AUC). Dashed line at 10% CV
  marks the manuscript's claimed threshold; medians, IQR, and fraction of
  genotypes below 10% CV are annotated in-panel. Secondary/diagnostic —
  regenerates locally (gitignored).
- `figures/replicate_cv_azt.png` — same layout for AZT. Secondary/diagnostic
  — regenerates locally (gitignored).
- `data/per_genotype_replicate_stats.csv` — one row per (genotype, drug,
  concentration); columns: `genotype, drug, concentration, log_mean,
  log_sd, raw_mean, raw_sd, cv_percent` (774,144 rows, design-filtered).
- `data/summary_by_conc.csv` — 14 rows (drug x concentration including
  conc=0) with median/IQR/quantiles of %CV and log_sd, and fraction of
  genotypes with %CV < 10.
- `data/replicate_pair_correlations.csv` — 42 rows: one per unique
  replicate pair (1-2, 1-3, 2-3) per drug × concentration; columns:
  `drug, concentration, rep_i, rep_j, n, pearson_r, pearson_p`. Most
  `pearson_p` values underflow to `0.0` at this sample size (double precision);
  the no-drug control (Aztreonam 0.0 µg/ml) has near-uniform fitness across
  genotypes, so its replicate correlation is ≈0 and its p is finite and even
  non-significant (max ≈0.10).
- `data/summary_by_conc_viable.csv` — 14 rows; same columns as
  `summary_by_conc.csv` but restricted to genotypes with mean
  log10(AUC) > 3.0 at the given condition. Confirms the "<10% SD"
  claim fails even in the viable subset.
- `results_table.csv` — 14-row headline verdict table: CV summary joined
  to mean Pearson r across pairs, with `median_cv_below_10pct` boolean.

## Design choices

- **Conc=0 kept in tables, dropped from figures.** Zero-drug is a pure
  growth control and its dispersion does not speak to the AUC-fitness
  claim, but it is included in `summary_by_conc.csv` for completeness.
- **Two reproducibility metrics, reported side by side.**
  - %CV on raw AUC = 100 * SD(10^fitness) / mean(10^fitness). Matches the
    typical reader interpretation of "SD below 10%."
  - SD of log10(AUC) directly. Matches the manuscript units.
- **Three replicate pairs per concentration**, shown as 3 sub-panels rather
  than an average, to display the scatter directly rather than collapsing
  it into a single averaged summary.
- **Hex-bin with log-count colour** scales cleanly across ~55k points
  without over-plotting; identity line in drug colour.
- **Pearson is the right summary** (replicates should cluster around the
  identity line, not just "rank-correlate"). Spearman would over-state
  agreement in this regime.

## Summary

- Reproducibility is quantified by 42 replicate-pair Pearson correlations
  (three pairs per drug × concentration) and per-genotype %CV, over the
  55,296 combinatorially-complete design genotypes (n = 55,293 with a finite
  replicate SD per condition). Pearson r ranges from -0.01 to 0.85 across
  pairs; most p-values underflow to 0.0 at this sample size; the largest (≈0.10)
  is at Aztreonam 0.0 µg/ml — the no-drug control — where fitness is near-uniform
  across genotypes, so replicate correlation is ≈0.
- Generated files: `figures/supplementary/figure_s16_amp.png`,
  `figures/supplementary/figure_s16_azt.png`,
  `data/per_genotype_replicate_stats.csv`, `data/summary_by_conc.csv`,
  `data/summary_by_conc_viable.csv`, `data/replicate_pair_correlations.csv`,
  and `results_table.csv`.
