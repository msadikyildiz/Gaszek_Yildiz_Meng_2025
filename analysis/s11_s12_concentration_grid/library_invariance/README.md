# `library_invariance/` — AMP 781 library-invariance test

Tests whether the main-text Fig 4/Fig 5 conclusions at AMP 781 µg/mL are
robust to dropping the sequencing-floor fraction of the library.

## Rerun

```bash
uv run python analysis/s11_s12_concentration_grid/library_invariance/analysis.py
```

Takes about 2 min on a single CPU (LightGBM + SHAP dominate runtime).

## What the script produces

1. **`data/amp781_majority_low.csv`** — per-genotype majority-low flag
   computed from the raw read-counts table. 70,183 rows, with an
   `in_landscape` boolean identifying the 55,296 genotypes used in
   `Epistasis_Combined.parquet`. Raw-library majority-low fraction is
   32.15 %; landscape-subset rate is 14.89 %.

2. **Pairwise test** (Fig 4 C).
   - `data/pairwise_comparison.csv` — 19×19 mean-fitness matrix in long form
     (full vs clean vs Δ) for all unordered mutation pairs.
   - `figures/pairwise_full_vs_clean_amp781.png` — side-by-side heatmaps +
     difference panel with annotated Pearson r.

3. **R² vs K test** (Fig 4 G).
   - `data/r2_comparison.csv` — 13 epistatic orders × 2 subsets (full,
     clean). Uses the manuscript's own partial-linear-regression
     predictions (`Fitness_predicted for order K`) already in the parquet
     — we do not re-fit.
   - `figures/r2_vs_order_full_vs_clean_amp781.png`.

4. **LightGBM + SHAP test** (Fig 5).
   - `data/shap_comparison.csv` — per-mutation mean |SHAP| and ranks from
     two independently fit LightGBM models (identical manuscript
     hyperparameters: `n_estimators=100, lr=0.1, seed=20260420`).
   - `figures/shap_full_vs_clean_amp781.png` — paired bar chart + rank
     scatter.

5. **Summary** (`data/invariance_summary.csv`, 5 rows).

## Verdict

We evaluated three thresholds (pairwise r ≥ 0.95, |ΔR²| ≤ 0.02 at every K,
SHAP ρ ≥ 0.90). The pairwise mean-fitness Pearson r (0.858) falls below the
0.95 threshold; this is driven by the matrix cells shifting upward (median
+0.16 log₁₀ AUC, concentrated on L21P-containing pairs) when the
sequencing-floor tail is removed — a change at the level of raw means, not
of model fit or claimed epistatic pattern. The biochemical-epistasis view
(Fig 5E/F) is numerically unchanged (r = 1.0000); R² vs K is strictly higher
on the clean subset at every K; LightGBM is more accurate (RMSD 0.37 vs
0.46). We therefore soften the invariance claim at the level of raw means
while keeping the exclusion rule.

## File layout

```
library_invariance/
  analysis.py              # 4-step pipeline
  data/*.csv               # 5 CSVs
  figures/*.png            # 3 PNGs
  README.md                # this file
```
