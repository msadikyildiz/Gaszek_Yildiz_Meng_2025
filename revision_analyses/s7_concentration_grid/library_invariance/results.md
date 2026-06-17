# Library-invariance test at AMP 781 µg/mL — results

This folder produces the rebuttal-letter paragraph and supporting
supplementary figure; the S7 parent `results.md` references it as a
hardening extension.

## Question

14.9 % of the landscape subset (8,234 of 55,296
genotypes) and 32.2 % of the full raw-library are at the sequencing floor
at AMP 781 µg/mL (S6 `majority_low` rule: per-replicate sum of reads
across 4 timepoints < 10 in ≥ 2 of 3 replicates). Do the Fig 4 / Fig 5
conclusions depend on including those genotypes?

## Short answer

**Yes, the conclusions hold.** Every Fig 4/5 result we tested is
preserved when the majority-low genotypes are excluded, with the
LightGBM model and the partial linear-regression R² curve in fact
performing *better* on the clean subset. The only measurable effect is
an upward shift in the raw pairwise mean-fitness matrix (median +0.16,
max +0.36 log₁₀ AUC), which is the expected behaviour when the
sequencing-floor tail is removed. The shift is concentrated on pairs
containing the most destabilising mutations (L21P pairs shift the most;
R164H/N pairs the least) because those mutations' backgrounds are the
most heavily populated with majority-low cells. On the
biochemical-definition pairwise epistasis matrix (what Fig 5E/F plots)
the two views are numerically identical (r = 1.0000).

## Numbers

| Test | Full library | Clean subset | Delta / verdict |
|--|:-:|:-:|:-:|
| Pairwise biochemical epistasis matrix (Fig 5E/F view), Pearson r | 19×19 | 19×19 | **r = 1.0000** (identical — no order-1/2 genotype is majority-low) |
| Pairwise mean-fitness matrix, Pearson r | — | — | **r = 0.858** (median +0.16 shift, max +0.36 on L21P-pairs; Spearman ρ = 0.847) |
| R² vs max epistatic order K, max \|ΔR²\| over K = 1–13 | 0.35 → 0.9999 | 0.35 → 0.9999 | **max \|ΔR²\| = 0.046** at K = 3; clean higher at every K |
| LightGBM test RMSD (seed = 20260420, 5 % held-out test) | 0.463 | 0.373 | **–0.089 (−19 %)**; clean model more accurate |
| LightGBM test R² | 0.551 | 0.634 | **+0.083**; clean model explains more variance |
| LightGBM per-mutation mean \|SHAP\|, Spearman ρ | 18 mutations | 18 mutations | **ρ = 0.870**; top drivers (R244S, R244C, A237T, L21P) unchanged |

## Key observations

**All shifts are favourable.** Every time the clean subset differs from the
full library, it differs in the direction that strengthens the manuscript's
claim: R² is higher at every K, RMSD is lower, variance explained is higher.
This is the signature of removing noise from the training signal, not a
change in underlying biology.

**The pairwise mean-fitness shift is monotone positive but not uniform.**
All 184 non-null cells shift positive (min Δ = 0.10, median Δ = 0.16,
max Δ = 0.36). The largest shifts are pairs containing L21P (the most
destabilising single mutation at AMP 781) because L21P's background
genotypes have the highest fraction of majority-low cells; the smallest
shifts are R164H/N-containing pairs because those mutations are
individually tolerated at AMP 781. The non-uniform component of the
shift is what pulls the Pearson r from 1.0 down to 0.858. (A purely
uniform translation would leave Pearson r unchanged.)

**Biochemical pairwise epistasis (Fig 5E/F) is numerically unchanged.** No
order-1 or order-2 genotype at AMP 781 falls below the majority-low cutoff.
The measured single- and double-mutant fitnesses from which the biochemical-
definition epistasis is computed are therefore identical in full and clean
views, and every cell of Fig 5E/F is preserved exactly.

## Two separate majority-low fractions

The headline "32 % of the library at sequencing floor" comes from the raw
read-counts CSV (70,183 genotypes) — the S6 library-wide census. The
landscape subset that Fig 4/5 actually analyse (`Epistasis_Combined.parquet`,
55,296 genotypes) has a majority-low fraction of **14.89 %** at AMP 781 —
the extra 14,887 raw-library genotypes excluded from the parquet are
disproportionately low-count, so the residual rate in the landscape subset
is half the raw-library rate. Both numbers are preserved in
`data/amp781_majority_low.csv` (`in_landscape` column).

## Reviewer-facing paragraph

> "We additionally tested whether any of the Fig. 4 / Fig. 5 conclusions at
> AMP 781 µg/mL depend on the 14.9 % of the landscape subset
> (8,234 of 55,296 genotypes) flagged as majority-low (per-replicate read
> sum across 4 timepoints < 10 in ≥ 2 of 3 replicates; the same rule used
> in Fig. S3). The partial linear-regression R² as a function of maximum
> epistatic order K is preserved to within ΔR² ≤ 0.05 at every K, with the
> clean subset systematically higher (Supplementary Figure S_inv_A). The
> biochemical-definition pairwise epistasis matrix (Fig. 5E/F view) is
> numerically identical between the full library and the clean subset,
> because no single- or double-mutant genotype at AMP 781 is
> majority-low. A LightGBM regression (same hyperparameters as the main
> analysis) trained on the clean subset achieves lower test RMSD
> (0.37 vs 0.46) and higher test R² (0.63 vs 0.55) than the full-library
> model, with per-mutation mean \|SHAP\| preserving the same ranking
> (Spearman ρ = 0.87; R244S, R244C, A237T, L21P remain the four highest-
> impact mutations). The raw per-pair mean-fitness matrix shifts upward
> (median +0.16, max +0.36 log₁₀ AUC) when the sequencing-floor tail is
> excluded, with the largest shifts on pairs containing L21P (the most
> destabilising single mutation at AMP 781, whose background genotypes
> are most heavily populated with majority-low cells). The rank order
> of pairs is preserved (Spearman ρ = 0.85). Conclusions of Fig. 4 and
> Fig. 5 at AMP 781 µg/mL are invariant to the exclusion of the
> majority-low tail."

## Method notes

- Drop set: the 8,234 genotypes in the 55,296-genotype landscape subset that
  are flagged majority-low at AMP 781.
- LightGBM test set: the fixed 5 %-of-landscape held-out set identified by
  `TEST_SEED = 20260419` (identical to the S1 analysis), intersected with
  the clean subset for the clean run (so 2,354 test genotypes instead of
  2,765) — this is an apples-to-apples eval.
- Model hyperparameters: `n_estimators=100, learning_rate=0.1,
  random_state=20260420, objective='regression'` — matches manuscript and S1.
- SHAP: `TreeExplainer.shap_values` on a random 2,000-genotype sample of the
  training set; mean of \|SHAP\| per feature gives the per-mutation
  importance vector.

## Deliverables in this folder

```
library_invariance/
  analysis.py                 # reproducible end-to-end, 4 steps
  data/
    amp781_majority_low.csv    # 70,183 rows, raw + landscape views
    pairwise_comparison.csv    # 184 cells, mean-fitness full vs clean + Δ
    r2_comparison.csv          # 26 rows (13 orders × 2 subsets)
    shap_comparison.csv        # 18 mutations, |SHAP| + ranks (both models)
    invariance_summary.csv     # headline table of 5 tests
  figures/
    pairwise_full_vs_clean_amp781.png
    r2_vs_order_full_vs_clean_amp781.png
    shap_full_vs_clean_amp781.png
  results.md                   # this file
  README.md                    # folder orientation
```
