# S2 — Block-holdout ML validation

Reviewer #4 point ii: block-holdout / mutation-stratified ML validation to
test whether the LightGBM fitness regressor learns transferable epistatic
rules rather than memorising local genotype neighbourhoods. See `results.md`
for the rebuttal-ready narrative.

## How to reproduce

```bash
# From repository root
export MPLCONFIGDIR=/tmp/claude/mpl_cache FONTCONFIG_CACHE=/tmp/claude/fc_cache
uv run python revision_analyses/s2_block_holdout/analysis.py
uv run python revision_analyses/s2_block_holdout/classification_check.py
```

`analysis.py` runs, for each of AMP @ 781 µg/mL and AZT @ 36 µg/mL:
- 18-fold leave-one-mutation-out (LOMO) cross-validation with LightGBM and
  a matched-N random-split baseline (also Lasso LOMO as a drug-class check).
- Distance-stratified Hamming holdout for D ∈ {2, …, 10} with two reference
  genotypes (wild type and a random deep mutant).
- Writes per-fold CSVs to `data/`, aggregated summary to `results_table.csv`,
  and two-panel publication PNGs to `figures/`.

`classification_check.py` adds clinical-classification diagnostics at
biologically motivated per-drug thresholds (AMP P25 loss-of-resistance; AZT
P90 gain-of-resistance). Output: `data/classification_results.csv`.

## Outputs

| File | Purpose |
|------|---------|
| `figures/block_holdout_amp.png` | Two-panel figure: LOMO bar chart + Hamming curve |
| `figures/block_holdout_azt.png` | Same for AZT |
| `data/lomo_results.csv`         | Per-mutation LightGBM + random + Lasso RMSD/R² |
| `data/hamming_results.csv`      | Per-D LightGBM + random RMSD/R² for both references |
| `data/classification_results.csv` | AUROC / balanced-accuracy under LOMO |
| `results_table.csv`             | Aggregated deliverable |
| `results.md`                    | Rebuttal response draft |

## Hyperparameters (held fixed)

Matches the manuscript's `src/04_epistasis_amp_regression.ipynb` and
`src/03_epistasis_azt_regression.ipynb`:

```python
lgbm.LGBMRegressor(
    objective='regression',
    n_estimators=100,
    learning_rate=0.1,
    random_state=42,
    verbose=-1,
)
```

Features: 18-dim one-hot of non-wild-type substitutions (same encoding as the
regression notebooks, minus the focal mutation column during LOMO folds).
