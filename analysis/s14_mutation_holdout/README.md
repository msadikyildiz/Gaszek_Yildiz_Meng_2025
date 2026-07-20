# s14 — mutation-stratified holdout validation (LOMO + Hamming distance)

Produces **Supplementary Figure S14**. Block-holdout / mutation-stratified
validation testing whether the LightGBM fitness regressor learns transferable
epistatic rules rather than memorising local genotype neighbourhoods. Added
during the Nature Communications revision.

## Rerun

```bash
# From repository root
uv run python analysis/s14_mutation_holdout/analysis.py
uv run python analysis/s14_mutation_holdout/classification_check.py
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

Both figures are the published Supplementary Fig. S14 and are written
directly to the shared gallery `figures/supplementary/`; this module has no
secondary/diagnostic panel, so its local `figures/` dir is unused
(gitignored).

| File | Purpose |
|------|---------|
| `figures/supplementary/figure_s14_amp.png` | Two-panel figure: LOMO bar chart + Hamming curve — **published Supplementary Fig. S14 (AMP)** |
| `figures/supplementary/figure_s14_azt.png` | Same for AZT — **published Supplementary Fig. S14 (AZT)** |
| `data/lomo_results.csv`         | Per-mutation LightGBM + random + Lasso RMSD/R² |
| `data/hamming_results.csv`      | Per-D LightGBM + random RMSD/R² for both references |
| `data/classification_results.csv` | AUROC / balanced-accuracy under LOMO |
| `results_table.csv`             | Aggregated deliverable |

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
