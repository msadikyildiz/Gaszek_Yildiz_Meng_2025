# Supplementary Figures S7 / S8 — monoculture IC50 validates pooled AUC-fitness

Caption context: independent monoculture dose-response validation of the
pooled-fitness landscape. **S7** — log10(IC50) vs. mean AUC-fitness, one panel
per drug, Pearson r + p annotated. **S8** — the same comparison broken out per
selection concentration (5 AMP panels, 7 AZT panels).

Headline numbers (`validation/S7_S8_provenance.md` derives the same statistics
independently from the same parquets via a minimal snippet):

| Drug | Batch    | n  | Pearson r | p        |
|------|----------|----|-----------|----------|
| AMP  | 20260407 | 13 | 0.885     | 5.9e-05  |
| AZT  | 20260307 | 18 | 0.803     | 6.1e-05  |

## Inputs

`analysis.py` reads the processed cross-reference tables
`validation/src/processed/<batch>/xref_expanded_df.parquet` (columns: drug,
variant, genotype_13, mean_fitness, log10_ic50, fitness_<concentration>), which
regenerate from the raw plate-reader XLSX in `validation/data/` via
`validation/src/run_all.py` (see `validation/README.md`). Nothing under
`validation/` is modified; the script only reads the processed parquets. Palette,
fonts, and the open-circle WT/DD convention match the manuscript's other
extended figures (drug-primary colours, DejaVu Sans).

## Contents

```
si_figures/s07_s08_ic50/
  README.md              (this file)
  analysis.py            repo-relative, runnable
  figures/               output of analysis.py (convenience copies)
    ext_ic50_vs_auc.png            same image as figures/supplementary/figure_s07.png
    ext_ic50_per_conc_amp.png      same image as figures/supplementary/figure_s08_amp.png
    ext_ic50_per_conc_azt.png      same image as figures/supplementary/figure_s08_azt.png
```

## Method

- **S7**: one scatter panel per drug, log10(IC50) (y) vs mean AUC-fitness (x),
  with a dashed regression line and annotated Pearson r, p, n. Rows with a
  defined `genotype_13` are included; DD (dead control) has a null
  `mean_fitness` and is excluded from the correlation but still plotted.
- **S8**: a grid of per-selection-concentration panels (AMP batch 20260407, AZT
  batch 20260307), each log10(IC50) vs the AUC-fitness at that concentration,
  coloured by a per-concentration ramp, with a per-panel Pearson r.
- Label placement uses a self-contained iterative label-repulsion helper
  (`_nudge_labels_2d`), so the script has no import dependency on the
  `validation/` package internals.

## Outputs

Written directly to the published Supplementary location:
`figures/supplementary/figure_s07.png`, `figure_s08_amp.png`,
`figure_s08_azt.png`, plus a convenience copy of each under this folder's
`figures/`.

## Reproduction

```bash
uv run python si_figures/s07_s08_ic50/analysis.py
```

Runs in the repo's `uv` environment (numpy, polars, matplotlib, scipy; no new
dependencies). Requires `validation/src/processed/<batch>/xref_expanded_df.parquet`;
if missing, regenerate first with `uv run python validation/src/run_all.py`.
Reproduces the headline statistics:

```
AMP: r = 0.885, p = 5.869e-05, n = 13
AZT: r = 0.803, p = 6.127e-05, n = 18
AMP per-concentration r: 0.02 (3.1, n.s.) -> 0.67, 0.81, 0.91, 0.88 (781, ***)
AZT per-concentration r: 0.61-0.80 (all *** or **)
```

## Caveat

DD (catalytically dead control) has a null `mean_fitness` in the processed
xref (excluded from the Pearson correlation) but is still plotted as an open
circle, matching the manuscript's TEM-1$^{dead}$ reference-point convention.
