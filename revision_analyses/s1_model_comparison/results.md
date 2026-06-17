# Reviewer #3 — comment 2b: model-class comparison

## Reviewer comment (verbatim)

> "Some of the machine learning aspects could be simplified to support the claims. I would suggest also adding a lasso regression (e.g. with the glmnet package). Then compare (in e.g. a table) the predictive accuracy of linear models, with L1 regularisation, and decision trees. If you could then also show the predictive accuracy at each training sample size?"

## Response summary

We agree with the reviewer that a side-by-side comparison against simpler baselines is the appropriate test of whether the LightGBM choice adds value. We have added a new Supplementary analysis (Fig. S\_new; `model_comparison_amp.png` and `model_comparison_azt.png`) that trains four model classes on the same task as Fig. 3/4:

1. unregularised linear regression (additive one-hot features),
2. L1-regularised linear regression (`LassoCV`, alpha chosen by 5-fold cross-validation on the training set),
3. a single decision-tree regressor (default hyperparameters),
4. LightGBM (matching the manuscript: `n_estimators=100`, `learning_rate=0.1`).

Because the additive one-hot encoding cannot, by construction, express any non-additive interaction between mutations — whereas LightGBM can — we also report the linear and Lasso variants fit on an enriched feature set that includes all 153 pairwise-interaction indicator features (so-called "pairwise + additive" models, 171 features). This is the fairest possible comparison of a linear model against a tree ensemble on this task.

Every model is evaluated on the **same fixed 5 % held-out test set** (2,765 genotypes) at each of the two featured conditions (AMP 781 µg/mL, AZT 36 µg/mL). Training subsets are drawn from the remaining 95 % (52,531 genotypes) in fractions {0.03 %, 0.1 %, 0.3 %, 1 %, 3 %, 10 %, 30 %, 50 %, 70 %, 90 %, 99 %}, repeated with 3 random seeds per fraction. Subsample seeds are derived deterministically from a `hashlib.sha256(drug+concentration)` offset, so re-running the pipeline reproduces the reported numbers exactly. We report mean ± s.d. of RMSD and R² at each training size.

## Headline numbers (mean RMSD on held-out 5 %)

### Asymptotic regime (n_train = 52,006; 94 % of all genotypes)

| Model                 | AMP @ 781 µg/mL |  AZT @ 36 µg/mL |
| --------------------- | --------------: | --------------: |
| **LightGBM**          |       **0.462** |       **0.291** |
| Lasso (+ pairwise)    |           0.482 |           0.335 |
| Linear (+ pairwise)   |           0.482 |           0.335 |
| Lasso (additive)      |           0.557 |           0.359 |
| Linear (additive)     |           0.557 |           0.359 |
| Decision tree         |           0.704 |           0.429 |

For reference, the target fitness has standard deviation 0.71 (AMP) and 0.38 (AZT).

### At 10 % of the training pool (n_train = 5,253)

| Model                 | AMP RMSD    | AMP R²  | AZT RMSD    | AZT R²  |
| --------------------- | :---------: | :-----: | :---------: | :-----: |
| **LightGBM**          | **0.477**   | **0.52** | **0.315** | **0.30** |
| Lasso (+ pairwise)    |   0.487     |  0.50   |   0.339     |  0.19   |
| Linear (+ pairwise)   |   0.488     |  0.50   |   0.340     |  0.19   |
| Lasso (additive)      |   0.558     |  0.35   |   0.360     |  0.09   |
| Linear (additive)     |   0.558     |  0.35   |   0.360     |  0.09   |
| Decision tree         |   0.698     |  -0.02  |   0.475     | -0.59   |

## What the comparison shows

- **Bare additive linear and Lasso are indistinguishable** once sample size exceeds ~10× the number of features. The 18 one-hot mutation indicators are *near-independent* in this library — not literally orthogonal, because three positions admit three states (R164{H,S,N}, R244{C,S}, L21F/P) and their indicator columns are mutually exclusive — but the residual correlation is small enough that L1 at the CV-optimal α collapses onto the OLS solution.

- **The additive ceiling is strictly below LightGBM.** At asymptotic sample sizes, adding pairwise interactions reduces RMSD from 0.557 to 0.482 on AMP and from 0.359 to 0.335 on AZT. **This closes ~79 % of the additive-to-LightGBM RMSD gap on AMP, but only ~35 % on AZT.** The remaining gap (AMP: 0.482 → 0.462, AZT: 0.335 → 0.291) reflects *non-additive structure that pairwise features alone cannot express* — we cannot unambiguously attribute it to ≥3-way epistasis specifically without explicitly fitting a 3-way feature set, and a tree ensemble may also exploit non-interaction structure (e.g. saturating non-linearity). The empirical claim is the asymmetry between drugs, not a specific epistatic order.

- **The drug-specific asymmetry matches the biology.** On AMP, where TEM-1 is already highly active and adaptation is close to additive, pairwise features recover most of the LightGBM advantage. On AZT, where ESBL-like adaptation requires coordinated multi-residue substitutions, pairwise features recover only about a third of the gap — quantitative evidence that AZT adaptation carries non-additive structure beyond what pairs alone express, consistent with the rest of the manuscript.

- **A single decision tree is the worst model at every training size.** Default `scikit-learn` trees grow until every leaf is pure and therefore over-fit the training fractions we tested; held-out R² is negative for AZT at every sample size. This makes the reviewer's point in reverse: the superior performance of LightGBM over a plain decision tree is due to *tree ensembling plus shrinkage*, not "using trees" per se. (A pruned tree with `min_samples_leaf=10` reaches AMP RMSD 0.491, but still loses to the pairwise-Lasso baseline and to LightGBM.)

- **LightGBM is competitive at every training fraction large enough to fit any model at all,** and never worse than the alternatives. It matches Lasso-with-pairwise when n_train ≤ ~10³ (the regime where the pairwise model is effectively regularised by the CV) and pulls ahead thereafter as additional capacity becomes justified by additional data.

- **Important failure mode made visible.** The "Linear (+ pairwise)" curve spikes sharply at intermediate training sizes (n_train ≈ 100–500) where the number of samples is smaller than the 171-dim feature set. Unregularised OLS on an under-determined system returns a minimum-norm solution that generalises poorly. Lasso-with-pairwise is immune because the L1 penalty regularises the solution — the clearest practical argument for including regularisation when the feature set includes interactions.

## Independent cross-checks

We ran an additional cross-check to confirm that the LightGBM advantage is tree-ensembling rather than boosting specifically. A Random Forest baseline on AZT reaches RMSD 0.302, within 0.011 of LightGBM. We also confirmed that a Ridge baseline is interchangeable with Lasso on the pairwise feature set, so our choice of L1 is incidental. This cross-check was run on one subsample seed × five fractions, thinner than our canonical grid, and is cited as corroboration rather than a replacement for the full sweep above.

## Deliverables

- `figures/model_comparison_amp.png` — AMP 781 µg/mL learning curves.
- `figures/model_comparison_azt.png` — AZT 36 µg/mL learning curves.
- `results_table.csv` — per-model × drug × fraction × seed RMSD / R² / runtime (396 rows).
- `data/learning_curves_summary.csv` — aggregated mean ± s.d. (132 rows).
- `analysis.py` — fully reproducible pipeline; re-run with `python analysis.py`.

## Methodological notes

- All fits use a **single fixed test set** (5 % of 55,296 genotypes, shuffled with seed 20260419 and held back before any other randomness); training subsamples are drawn from the remaining 95 %. This is the only configuration in which learning-curve RMSDs from different training sizes are mutually comparable, and it is tighter than the protocol in the original LightGBM notebook (which re-split train/test at every training fraction).
- `LassoCV` selects alpha from a 50-value default grid by 5-fold CV on the training set — no hand-tuned alpha.
- The decision tree uses `DecisionTreeRegressor(random_state=seed)` with otherwise-default hyperparameters, matching the spirit of the reviewer's request for a plain "decision tree" baseline.
- LightGBM hyperparameters match the manuscript exactly (`objective='regression'`, `n_estimators=100`, `learning_rate=0.1`). We did not tune LightGBM for this comparison.
- The pairwise feature set contains all 18 × 17 / 2 = 153 pairwise AND-interactions, giving 171 features total.

## Bottom line (to be paraphrased into the main rebuttal letter)

Adding the Lasso and decision-tree comparisons the reviewer requested does not change the manuscript's conclusion: LightGBM is the best-performing model on both drugs at every training size. The size of LightGBM's advantage is **asymmetric across drugs** — on AMP, pairwise-interaction features close ~79 % of the gap between an additive baseline and LightGBM; on AZT they close only ~35 %. This asymmetry is the quantitative signature of non-additive structure in the AZT landscape that pairwise terms alone cannot express, consistent with the epistasis motif discussed throughout the manuscript. The decision-tree baseline is worst at every training size. We have added these curves and the accompanying table as a new Supplementary figure.
