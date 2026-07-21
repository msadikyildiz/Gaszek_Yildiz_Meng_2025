# s13 — model comparison: LightGBM vs linear, Lasso, and decision-tree regression

Produces **Supplementary Figure S13**. Compares LightGBM against
unregularised linear regression, L1 (Lasso) regression, and a single decision
tree, reporting predictive accuracy as a function of training sample size.

## Rerun

```
export LOKY_MAX_CPU_COUNT=8
uv run python analysis/s13_model_comparison/analysis.py
```

Runtime: ~1 minute on 8 cores.

## Inputs

- `data/processed/Epistasis_Combined.parquet`
  (55,296 genotypes × 2 drugs × 6–8 concentrations; `Fitness` column is
  log10 AUC).

## Outputs

Both learning-curve panels are the published Supplementary Fig. S13 and are
written directly to the shared gallery `figures/supplementary/`; this module
has no secondary/diagnostic panel, so its local `figures/` dir is unused
(gitignored).

- `figures/supplementary/figure_s13_amp.png` — AMP @ 781 µg/mL learning curves. **Published Supplementary Fig. S13 (AMP)**.
- `figures/supplementary/figure_s13_azt.png` — AZT @ 36 µg/mL learning curves. **Published Supplementary Fig. S13 (AZT)**.
- `results_table.csv` — one row per (model, drug, fraction, seed): RMSD, R², runtime (396 rows; raw per-seed data).
- `data/learning_curves_summary.csv` — aggregated mean ± s.d. over seeds for every (model, drug, fraction) (132 rows; the spec-shaped table).

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
