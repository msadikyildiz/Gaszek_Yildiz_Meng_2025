# S1 — model-comparison analysis

Response to Reviewer #3 comment 2b: compare LightGBM against unregularised linear
regression, L1 (Lasso) regression, and a single decision tree; report predictive
accuracy as a function of training sample size.

## Rerun

```
export MPLCONFIGDIR=/tmp/claude/mpl_cache FONTCONFIG_CACHE=/tmp/claude/fc_cache
export LOKY_MAX_CPU_COUNT=8
mkdir -p $MPLCONFIGDIR $FONTCONFIG_CACHE
uv run python revision_analyses/s1_model_comparison/analysis.py
```

Runtime: ~1 minute on 8 cores.

## Inputs

- `Gaszek_Yildiz_Meng_2025/data/processed/Epistasis_Combined.parquet`
  (55,296 genotypes × 2 drugs × 6–8 concentrations; `Fitness` column is
  log10 AUC).

## Outputs

- `figures/model_comparison_amp.png` — AMP @ 781 µg/mL learning curves.
- `figures/model_comparison_azt.png` — AZT @ 36 µg/mL learning curves.
- `results_table.csv` — one row per (model, drug, fraction, seed): RMSD, R², runtime (396 rows; raw per-seed data).
- `data/learning_curves_summary.csv` — aggregated mean ± s.d. over seeds for every (model, drug, fraction) (132 rows; the spec-shaped table).
- `results.md` — rebuttal-response draft with quoted reviewer text + conclusions.

## Design choices

- **Fixed 5 % test set** across every model, training fraction, and seed
  (seed 20260419 permutes genotypes once). Training samples are drawn from the
  remaining 95 %. This is the only configuration in which learning-curve RMSDs
  from different training sizes are mutually comparable.
- **Training fractions** expressed as fractions of the 95 % training pool:
  {0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99}.
  At frac=0.99, n_train = 52,006 (~94 % of all 55,296 genotypes).
- **3 random subsample seeds** per training fraction; report mean ± s.d.
- **Features**. Bare additive: one binary column per non-WT substitution =
  18 features. Pairwise: all 18 × 17 / 2 = 153 pairwise AND-interactions
  appended = 171 features.
- **Models**.
  - `LinearRegression` (scikit-learn, no regularisation).
  - `LassoCV(cv=5, alphas=50, max_iter=20000)` — alpha chosen on training set.
  - `DecisionTreeRegressor(random_state=seed)` (default hyperparameters).
  - `LGBMRegressor(objective='regression', n_estimators=100, learning_rate=0.1,
    random_state=seed, verbose=-1)` — matches the manuscript.
- Linear and Lasso are fit on both bare and pairwise features; tree and
  LightGBM are fit only on bare features (they learn interactions internally).
