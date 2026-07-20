# s15 — classification metrics for the 10%-trained fitness model

Produces **Supplementary Figure S15**. Translates the LightGBM regression error
(RMSD) into explicit classification performance at clinically-relevant resistance
thresholds, and benchmarks that error against an estimated replicate-noise floor.
This grounds the claim that a model trained on 10% of the landscape is sufficient
to screen for high-fitness (resistance-driving) variants.

## Rerun

```
uv run python analysis/s15_classification_metrics/analysis.py
```

Runtime: a few minutes single-core (LightGBM over 10 random 10/90 splits per drug).

## Inputs

- `data/processed/Epistasis_Combined.parquet` — measured fitness (`Fitness`,
  median log10(AUC)) and per-order regression predictions. No model is refit here
  beyond the LightGBM used for the classification read-out.

## Method

- LightGBM trained on 10% of the landscape, evaluated on the held-out 90%,
  averaged over 10 random splits, at AZT 36 µg/mL and AMP 781 µg/mL.
- **Noise-floor benchmark.** Per-replicate residuals around the per-genotype mean
  give σ_single; a bootstrap Gaussian model (20,000 draws) propagates this to the
  σ of the median-of-three target (the regression target). The model RMSD is
  reported as a multiple of this estimated floor.
- **Classification.** Predictions are thresholded at three resistance definitions
  — above WT, top 5% (95th percentile), top 1% (99th percentile) — and scored by
  sensitivity, specificity, PPV, NPV, AUROC and AUPR, pooled over the 10 seeds.

## Outputs

- `figures/supplementary/figure_s15_azt.png`, `figures/supplementary/figure_s15_amp.png` —
  three-panel figures (ROC, confusion matrix, measured-vs-predicted scatter),
  pooled across the 10 seeds. **Published Supplementary Fig. S15**, written
  directly to the shared gallery; this module has no secondary/diagnostic
  panel, so its local `figures/` dir is unused (gitignored).
- `results_table.csv` — per-threshold sens/spec/PPV/NPV/AUROC/AUPR with ± std,
  plus the `sigma_median3_estimated_floor` noise-floor column.

## Headline numbers

| Drug | Model RMSD | R² | Est. noise floor | RMSD ÷ floor |
|------|-----------:|---:|-----------------:|-------------:|
| AZT @ 36 µg/mL  | 0.310 ± 0.003 | 0.34 | 0.213 | 1.46× |
| AMP @ 781 µg/mL | 0.482 ± 0.002 | 0.53 | 0.350 | 1.38× |

Extreme-tail classification is strong where a screen would operate: at top 1% the
AZT model reaches sensitivity 0.15 / specificity 1.00 / PPV 0.99 / AUROC 0.99; at
top 5% the AMP model reaches sensitivity 0.57 / specificity 1.00 / AUROC 0.99. The
model cannot resolve fitness differences within the bulk distribution below the
≈0.2 log-unit measurement noise, so performance at the near-WT threshold is
intentionally reported alongside the tail thresholds.
