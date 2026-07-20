# S4 — Epistasis section rewrite and two new extended figures

Responds to Reviewer #3 (restructure the epistasis section to lead with
the overall linear-model results rather than with specific sequences/
mutations) and to Reviewer #1 comments R1-22 (comparison of prediction
methods; referred to S1) and R1-23 (aztreonam predictability and
comparison with other β-lactams).

## Reviewer comments (verbatim)

**Reviewer #3 (primary driver for S4):**

> "The description of the epistatic interactions a bit muddled.
> Probably starting with the overall results from the linear model(s),
> rather than talking about specific sequences/mutations would help
> improve this. Epistatic interactions sections will be extensively
> revised. New extended figures will be added."

**Reviewer #1 (R1-22):**

> "The authors describe two methods for predicting the fitness levels
> of genotypes: regression with increasing epistatic order and machine
> learning. However, they never discuss or compare the quality of the
> predictions, nor do they address the advantages or disadvantages of
> these methods for future use. Further discussion of the generality
> of this dataset and approach is necessary."

**Reviewer #1 (R1-23):**

> "Aztreonam predictability seems much lower than ampicillin,
> suggesting difficulty in precisely predicting adaptation to novel
> phenotypes ... How do results from aztreonam selection compare to
> those from other beta-lactams?"

## Response summary

We have (a) restructured the Epistasis section so that it now opens
with a model-level description of how the partial-linear-regression R²
grows with epistatic order before narrating residue-specific results,
(b) added two new extended figures (A and E) that make the model-level
opening concrete, and (c) inserted a 3-sentence bridge paragraph to
Supplementary Note S1 (model-class comparison) so that the reviewer's
predictive-method comparison concern is addressed in the same section.

R1-22 is primarily answered by S1 (`revision_analyses/s1_model_comparison/`);
the new Epistasis section refers to S1's headline numbers rather than
duplicating them. R1-23 is directly addressed by Extended Fig. E and the
"matched-stringency" analysis below.

The full rewritten section text lives in `draft.md` (~820 words, ready
for integration as a track-changes block).

## Headline numbers (from new Extended Figures A and E)

**Ext Fig. A — order decomposition.** For every one of the 12 measured
drug × concentration conditions the partial linear-regression R² grows
monotonically with included epistatic order K. Additive R² (K = 1) and
cumulative R² at selected K from `data/order_decomposition.csv`:

| Drug | Conc (µg/mL) | R²(K=1) | R²(K=2) | ΔR² at K=2 | R²(K=6) | K at which R² ≥ 0.95 |
|---|---:|---:|---:|---:|---:|---:|
| AMP |   3.1 | 0.241 | 0.351 | +0.110 | 0.611 | 10 |
| AMP |  12.2 | 0.551 | 0.667 | +0.116 | 0.861 | 9 |
| AMP |  48.8 | 0.641 | 0.720 | +0.079 | 0.898 | 9 |
| AMP | 195.0 | 0.556 | 0.705 | +0.149 | 0.873 | 9 |
| AMP | **781.0** *(main-text)* | 0.351 | 0.525 | +0.174 | 0.722 | 9 |
| AZT |  0.44 | 0.422 | 0.604 | +0.182 | 0.872 | 9 |
| AZT |  1.33 | 0.349 | 0.560 | +0.210 | 0.849 | 9 |
| AZT |  4.00 | 0.259 | 0.478 | +0.220 | 0.796 | 10 |
| AZT |  12.0 | 0.169 | 0.367 | +0.198 | 0.730 | 10 |
| AZT | **36.0** *(main-text)* | 0.100 | 0.250 | +0.150 | 0.654 | 10 |
| AZT | 108.0 | 0.040 | 0.128 | +0.088 | 0.600 | 10 |
| AZT | 324.0 | 0.063 | 0.099 | +0.037 | 0.532 | 11 |

R² crosses 0.95 at K = 9 – 11 at every condition (average K = 9.8) and
crosses 0.99 at K = 11 – 12. The full-order (K = 13) fit reaches R² =
0.999 – 0.99998 (figures above match the S7 numbers to machine
precision; max |ΔR²| = 0.0 between S4 and S7 regression tables).

**Ext Fig. E — drug asymmetry at matched selective stringency.** At the
main-text matched pair (AMP 781 µg/mL, AZT 36 µg/mL — both in the
above-WT-MIC selective window per S7) the drug asymmetry is:

| Metric | AMP 781 | AZT 36 | ratio |
|---|---:|---:|---:|
| additive R² (K=1) | 0.351 | 0.100 | 3.5× |
| pairwise R² (K=2) | 0.525 | 0.250 | 2.1× |
| order-6 R² | 0.722 | 0.654 | 1.1× |
| order-13 R² | 0.99990 | 0.99952 | 1.0× |
| additive RMSD (log₁₀ AUC) | 0.569 | 0.362 | 1.6× |
| measured fitness σ | 0.71 | 0.38 | 1.9× |
| median \|biochem. pairwise ε\| | 0.14 | 0.29 | 0.48× |
| max \|biochem. pairwise ε\| | 1.98 | 2.68 | 0.74× |

At an alternative matched pair (AMP 195 µg/mL vs AZT 12 µg/mL — each ~3×
WT-MIC) the additive asymmetry remains 3.3× (AMP 0.556 vs AZT 0.169).
The asymmetry is therefore not an artefact of concentration choice: the
aztreonam landscape is genuinely more rugged than the ampicillin
landscape at matched stringency. This is the quantitative answer to
R1-23 within our own data; the comparison with VIM-2 (Chen, Fowler &
Tokuriki 2022) cited in `draft.md` provides cross-β-lactamase context.

## Where R1-22 is answered

The Epistasis section draft now contains a 3-sentence bridge paragraph
to Supplementary Note S1. S1 shows:

- LightGBM beats every baseline at both drugs and every training size
- Pairwise-feature Lasso closes ≈79 % of the additive-to-LightGBM RMSD
  gap on AMP but only ≈35 % on AZT
- The drug-specific asymmetry matches the order-decomposition here:
  where pairwise terms suffice (AMP), Lasso-with-pairwise nearly catches
  LightGBM; where they don't (AZT), LightGBM retains a meaningful
  advantage.

We do not duplicate S1's content — the rewritten Epistasis section cites
S1 by cross-reference and contributes the *linear* half of the
prediction-methodology story that S1 deliberately left under-emphasised.

## Figures

| File | Description |
| --- | --- |
| `figures/ext_order_decomposition_amp.png` | Ext Fig. A (AMP only) — cumulative R² line + ΔR² bars, 5 concentrations. |
| `figures/ext_order_decomposition_azt.png` | Ext Fig. A (AZT only) — cumulative R² line + ΔR² bars, 7 concentrations. |
| `figures/ext_order_decomposition_combined.png` | Ext Fig. A combined AMP + AZT, stacked — primary publication version. |
| `figures/ext_drug_asymmetry.png` | Ext Fig. E — 2×2 panel: (a) R² vs K, (b) RMSD vs K, (c) fitness-σ comparison, (d) pairwise \|ε\| histogram. Primary and secondary matched pairs annotated. |

## Data files

| File | Rows | Content |
| --- | ---: | --- |
| `data/order_decomposition.csv` | 156 | drug × conc × order (1–13) × R², ΔR², RMSD, n |
| `data/matched_stringency_summary.csv` | 52 | primary + secondary matched pairs × drug × order × R², RMSD |
| `data/source_data.xlsx` | 6 sheets | Ext Fig A-AMP, Ext Fig A-AZT, Ext Fig E R²/RMSD, Ext Fig E fitness dist, Ext Fig E pairwise-ε, matched-stringency summary |
| `results_table.csv` | 4 | Compact wide summary: matched-stringency R² and RMSD at K = 1, 2, 3, 6, 13 |

## Methodological notes

- All R², RMSD, and ΔR² values are derived from the manuscript's own
  pre-computed `Fitness_predicted for order K` columns in
  `Gaszek_Yildiz_Meng_2025/data/processed/Epistasis_Combined.parquet`.
  No new model is fitted.
- R² is Pearson correlation squared (matching the main-text Fig. 4
  convention in `src/05_epistasis_figures.ipynb`). ΔR² is the first
  difference along K, with R²(K = 0) := 0.
- Cross-check with `revision_analyses/s7_concentration_grid/data/regression_r2_by_order.csv`:
  max |ΔR²(S4) − R²(S7)| = 0.0 across all 156 rows (machine-exact).
- Pairwise-epistasis magnitudes read from S7's
  `data/pairwise_mean_fitness.csv` (column `biochemical_epistasis`), the
  same column used to build Fig. 5E/F in the main text.
- Seed: 20260420 (though no stochastic step is performed — reported for
  reproducibility). `mpl.rcParams['pdf.fonttype'] = 42` so figures can
  be converted to PDF without font embedding issues.

## Bottom line for the rebuttal letter

> "We have restructured the Epistasis section to open with the
>  partial-linear-regression R² as a function of included epistatic
>  order (new Extended Fig. A), then describe the concentration- and
>  drug-specific *shape* of that curve before narrating residue-
>  specific results. The new Extended Fig. E directly quantifies the
>  drug asymmetry the reviewers flagged: at matched selective
>  stringency the additive R² is 0.35 for ampicillin and 0.10 for
>  aztreonam (a 3.5-fold asymmetry that persists at an alternative
>  matched pair). The section now also refers forward to Supplementary
>  Note S1 for a full comparison of regression- and machine-learning-
>  based prediction accuracy, addressing Reviewer #1's request for
>  explicit method comparison."
