# S7 — concentration-grid extension of Figure 4 (and pairwise-epistasis panel of Figure 5)

Response to **Reviewer #1** (and parallel concerns from Reviewer #3) that the
manuscript's epistatic / predictive conclusions in Figures 4 and 5 are shown
at only one representative concentration per drug (AMP 781 µg/mL,
AZT 36 µg/mL) and need to be validated across the full measured concentration
range.

See `results.md` for the full rebuttal-response draft with reviewer
quotations and headline numbers.

## Rerun

```
export MPLCONFIGDIR=/private/tmp/claude/mpl_cache FONTCONFIG_CACHE=/private/tmp/claude/fc_cache
mkdir -p $MPLCONFIGDIR $FONTCONFIG_CACHE
uv run python revision_analyses/s7_concentration_grid/analysis.py
```

Runtime: ~6 s on CPU (data is pre-computed in the parquet — we only
aggregate and plot).

## Inputs

- `Gaszek_Yildiz_Meng_2025/data/processed/Epistasis_Combined.parquet`
  (774,144 rows = 55,296 genotypes × 2 drugs × 6–8 concentrations).
  Includes `Fitness` (log₁₀ AUC), `Biochemical Definition`
  (pairwise-epistasis measure), and pre-computed
  `Fitness_predicted for order K` for K ∈ {1, …, 13} produced by
  `src/utils/calculate_epistasis.py`.

## Outputs

### Figures

| File | Description |
| --- | --- |
| `figures/pairwise_heatmap_amp.png` | 19×19 pairwise-epistasis heatmap, one panel per AMP concentration (5 panels). |
| `figures/pairwise_heatmap_azt.png` | 19×19 pairwise-epistasis heatmap, one panel per AZT concentration (7 panels). |
| `figures/pred_vs_measured_amp.png` | Measured vs. predicted fitness density, rows = 5 AMP concentrations, columns = max. included epistatic order K ∈ 1..6. |
| `figures/pred_vs_measured_azt.png` | Measured vs. predicted fitness density, rows = 7 AZT concentrations, columns = max. included epistatic order K ∈ 1..6. |
| `figures/r2_vs_order.png` | Two-panel line plot: R² vs. maximum epistatic order, one curve per concentration, AMP (left) and AZT (right). |

### Data

| File | Shape | Content |
| --- | --- | --- |
| `data/pairwise_mean_fitness.csv` | 1,980 rows | drug, concentration, mut_i, mut_j, biochemical_epistasis (**note**: the authoritative content is the `biochemical_epistasis` column, not a mean-fitness value) |
| `data/regression_r2_by_order.csv` | 156 rows | drug, concentration, order, r2, rmsd, n |
| `results_table.csv` | 12 rows (2 drugs × 5 or 7 concs) | Wide R² table: r2_order_1, …, r2_order_13 |

## Design choices

- **Reuses manuscript's own partial-linear-regression predictions** stored in
  the parquet (`Fitness_predicted for order K`). Does *not* re-fit a separate
  model. This guarantees apples-to-apples comparison with main-text Fig. 4.
- **Concentrations** (zero-drug excluded):
  - AMP: 3.1, 12.2, 48.8, 195.0, 781.0 µg/mL (5 panels).
  - AZT: 0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0 µg/mL (7 panels).
- **Per-drug shared colour scales.** Pairwise heatmaps use symmetric
  RdBu_r with vlim = ±3 for AMP and ±4 for AZT, matching the
  manuscript's Fig. 5 E/F limits. Density panels use `Greys` (AMP) and
  `RdPu` (AZT), matching the manuscript's drug-colour convention.
- **Orders shown in density panels**: K ∈ {1, 2, 3, 4, 5, 6}. Main-text Fig.
  4 shows K ∈ {1, 2, 3}. We extend to K = 6 so the
  convergence to y = x is visible.
- **Orders in R² line plot**: K ∈ {1, …, 13}, the full range. R² reaches
  ≥ 0.98 by K = 10 and ≥ 0.99 by K = 11 at every concentration; K = 13
  asymptotes to near-exactness (R² between 0.99925 and 0.99998 depending
  on condition) — the partial linear-regression basis is over-determined
  relative to the 55,296-genotype library, so K = 13 is near-exact but
  not bit-exact.

## Related work in this response set

- `../s1_model_comparison/` — reviewer #3 comment 2b (model-class comparison).
- `../s2_block_holdout/` — reviewer #4 block-holdout / mutation-stratified
  validation.
- `../s3_rmsd_justification/` — reviewer #3/#4 sensitivity/specificity
  justification of the "RMSD ∼0.3 sufficient" claim.

## File layout

```
s7_concentration_grid/
├── README.md                                        (this file)
├── results.md                                       (draft rebuttal text)
├── analysis.py                                      (441 lines; <500 LOC)
├── results_table.csv                                (12 rows)
├── data/
│   ├── pairwise_mean_fitness.csv
│   └── regression_r2_by_order.csv
└── figures/
    ├── pairwise_heatmap_amp.png
    ├── pairwise_heatmap_azt.png
    ├── pred_vs_measured_amp.png
    ├── pred_vs_measured_azt.png
    └── r2_vs_order.png
```
