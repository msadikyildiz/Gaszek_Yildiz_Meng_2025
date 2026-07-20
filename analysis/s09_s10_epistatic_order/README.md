# S4 — Epistasis section rewrite and two new extended figures

Addresses **Reviewer #3** (restructure the Epistasis section to lead with
the overall linear-model results) and reinforces answers to **Reviewer
#1 R1-22** (prediction-method comparison, primarily handled in S1) and
**Reviewer #1 R1-23** (aztreonam vs other β-lactams).

See `results.md` for the full rebuttal-format response and
`draft.md` for the rewritten Epistasis-section text (~820 words) ready
for integration.

## Rerun

```
export MPLCONFIGDIR=/private/tmp/claude/mpl_cache FONTCONFIG_CACHE=/private/tmp/claude/fc_cache
mkdir -p $MPLCONFIGDIR $FONTCONFIG_CACHE
uv run python revision_analyses/s4_epistasis_rewrite/analysis.py
```

Runtime: ~10 s on CPU. No GPU, no model refit — the parquet already
stores every `Fitness_predicted for order K` column we need.

## Inputs

- `Gaszek_Yildiz_Meng_2025/data/processed/Epistasis_Combined.parquet`
  — 774,144 rows, 50 columns; the manuscript's master parquet with
  measured `Fitness` and `Fitness_predicted for order K` for
  K ∈ {1, …, 13}.
- `revision_analyses/s7_concentration_grid/data/regression_r2_by_order.csv`
  — used for cross-check only (S4 R² must match S7 R² to machine
  precision). This dependency is optional; if the S7 file is missing the
  analysis still runs without the assertion.
- `revision_analyses/s7_concentration_grid/data/pairwise_mean_fitness.csv`
  — used to derive the pairwise-epistasis magnitude distributions
  shown in Ext Fig. E panel (d).

## Outputs

### Figures

| File | Description |
| --- | --- |
| `figures/ext_order_decomposition_amp.png` | Ext Fig. A (AMP only). Cumulative R² line + ΔR² bars per concentration. |
| `figures/ext_order_decomposition_azt.png` | Ext Fig. A (AZT only). |
| `figures/ext_order_decomposition_combined.png` | Ext Fig. A combined — primary publication panel. |
| `figures/ext_drug_asymmetry.png` | Ext Fig. E, 2×2 panels: R², RMSD, fitness distribution, pairwise-ε histogram at matched stringency. |

### Data

| File | Shape | Content |
| --- | --- | --- |
| `data/order_decomposition.csv` | 156 rows | drug × concentration × order × R², ΔR², cumulative R², RMSD, n |
| `data/matched_stringency_summary.csv` | 52 rows | primary + secondary matched-stringency pairs × drug × order × R², RMSD |
| `data/source_data.xlsx` | 6 sheets | one sheet per Extended Figure panel + matched-stringency summary |
| `results_table.csv` | 4 rows | compact summary of matched-stringency R² and RMSD at K ∈ {1, 2, 3, 6, 13} |

## Design choices

- **Reuses manuscript predictions.** All R² and RMSD values are read
  from `Fitness_predicted for order K` in the parquet; no refit.
- **ΔR² convention.** R²(0) := 0, so ΔR²(K = 1) = R²(K = 1). This makes
  the full stack of ΔR² bars sum to the cumulative R² at any K.
- **Matched-stringency pairs (Ext Fig. E).** Primary pair is the main-
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

## Related work in this response set

- `../s1_model_comparison/` — Reviewer #3 2b (model-class comparison).
  The rewritten Epistasis section references S1 for the predictive-
  method-comparison bridge paragraph.
- `../s7_concentration_grid/` — Reviewer #1 (concentration-grid
  extension). S7's R² and pairwise-epistasis CSVs are direct inputs
  to S4.
- `../s2_block_holdout/` — Reviewer #4 block-holdout / mutation-
  stratified validation.
- `../s3_rmsd_justification/` — Reviewer #3/#4 sensitivity/specificity
  justification.

## File layout

```
s4_epistasis_rewrite/
├── README.md                      (this file)
├── draft.md                       (rewritten Epistasis section, ~820 words)
├── results.md                     (rebuttal-format response + quoted comments)
├── analysis.py
├── results_table.csv              (4 rows)
├── data/
│   ├── order_decomposition.csv
│   ├── matched_stringency_summary.csv
│   └── source_data.xlsx
└── figures/
    ├── ext_order_decomposition_amp.png
    ├── ext_order_decomposition_azt.png
    ├── ext_order_decomposition_combined.png
    └── ext_drug_asymmetry.png
```
