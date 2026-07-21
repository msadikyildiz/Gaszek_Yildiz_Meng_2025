# s09/s10 — epistatic-order decomposition and drug asymmetry

Produces **Supplementary Figures S9** and **S10**. S9 decomposes the
manuscript's partial linear-regression R² by maximum included epistatic order
K, showing how much of the fitness-landscape variance is captured as
pairwise, three-way, and higher-order terms are added. S10 compares AMP vs
AZT predictability (R², RMSD, pairwise-epistasis magnitude) at matched
selective stringency, testing whether the apparent drug asymmetry reflects
genuine epistatic-structure differences rather than the choice of
concentration.

## Rerun

```
uv run python analysis/s09_s10_epistatic_order/analysis.py
```

Runtime: ~10 s on CPU. No GPU, no model refit — the parquet already
stores every `Fitness_predicted for order K` column we need.

## Inputs

- `data/processed/Epistasis_Combined.parquet`
  — 774,144 rows, 50 columns; the manuscript's master parquet with
  measured `Fitness` and `Fitness_predicted for order K` for
  K ∈ {1, …, 13}.
- `analysis/s11_s12_concentration_grid/data/regression_r2_by_order.csv`
  — used for cross-check only (S9/S10 R² must match S11/S12 R² to machine
  precision). This dependency is optional; if the file is missing the
  analysis still runs without the assertion.
- `analysis/s11_s12_concentration_grid/data/pairwise_mean_fitness.csv`
  — used to derive the pairwise-epistasis magnitude distributions
  shown in Supplementary Fig. S10 panel (d).

## Outputs

The two published Supplementary panels (S9, S10) are written directly to the
shared gallery `figures/supplementary/`. The single-drug decomposition panels
are secondary/diagnostic and regenerate into this module's local `figures/`
(gitignored, scratch).

### Figures

| File | Description |
| --- | --- |
| `figures/s09_order_decomposition_amp.png` | Supplementary Fig. S9 (AMP only). Cumulative R² line + ΔR² bars per concentration. Secondary/diagnostic — regenerates locally (gitignored). |
| `figures/s09_order_decomposition_azt.png` | Supplementary Fig. S9 (AZT only). Secondary/diagnostic — regenerates locally (gitignored). |
| `figures/supplementary/figure_s09.png` | Supplementary Fig. S9 combined — **published Supplementary Fig. S9**. |
| `figures/supplementary/figure_s10.png` | Supplementary Fig. S10, 2×2 panels: R², RMSD, fitness distribution, pairwise-ε histogram at matched stringency — **published Supplementary Fig. S10**. |

### Data

| File | Shape | Content |
| --- | --- | --- |
| `data/order_decomposition.csv` | 156 rows | drug × concentration × order × R², ΔR², cumulative R², RMSD, n |
| `data/matched_stringency_summary.csv` | 52 rows | primary + secondary matched-stringency pairs × drug × order × R², RMSD |
| `data/source_data.xlsx` | 6 sheets | one sheet per figure panel + matched-stringency summary |
| `results_table.csv` | 4 rows | compact summary of matched-stringency R² and RMSD at K ∈ {1, 2, 3, 6, 13} |

## Design choices

- **Reuses manuscript predictions.** All R² and RMSD values are read
  from `Fitness_predicted for order K` in the parquet; no refit.
- **ΔR² convention.** R²(0) := 0, so ΔR²(K = 1) = R²(K = 1). This makes
  the full stack of ΔR² bars sum to the cumulative R² at any K.
- **Matched-stringency pairs (Supplementary Fig. S10).** Primary pair is the main-
  text concentration per drug (AMP 781 µg/mL, AZT 36 µg/mL — the
  just-above-WT-MIC selective point). Secondary pair is AMP 195 µg/mL
  vs AZT 12 µg/mL (each ~3× WT-MIC) — included to demonstrate that the
  drug asymmetry is not an artefact of the primary concentration
  choice.
- **Drug-colour convention.** AMP = greyscale ramp, AZT = RdPu ramp —
  matches manuscript Figs 4, 5 and S7.
- **Fonts.** `pdf.fonttype = 42` so PNGs can be regenerated as
  searchable PDFs at manuscript integration.
- **LOC.** `analysis.py` is 497 non-blank non-comment non-docstring
  lines (project rule: ≤ 500).

## Related analyses

- `../s13_model_comparison/` — model-class comparison (linear / Lasso /
  decision tree / LightGBM). Referenced by the epistatic-order discussion
  for the predictive-method-comparison bridge.
- `../s11_s12_concentration_grid/` — concentration-grid extension. Its
  R² and pairwise-epistasis CSVs are direct inputs here.
- `../s14_mutation_holdout/` — mutation-stratified holdout validation.
- `../s15_classification_metrics/` — sensitivity/specificity/RMSD
  classification metrics.

## File layout

```
s09_s10_epistatic_order/
├── README.md                      (this file)
├── analysis.py
├── results_table.csv              (4 rows)
├── data/
│   ├── order_decomposition.csv
│   ├── matched_stringency_summary.csv
│   └── source_data.xlsx
└── figures/
    ├── s09_order_decomposition_amp.png
    └── s09_order_decomposition_azt.png
```
