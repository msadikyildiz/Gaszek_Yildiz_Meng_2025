# Supplementary Figures S7 / S8 — monoculture IC50 validates pooled AUC-fitness

Caption context: independent monoculture dose-response validation of the
pooled-fitness landscape (Reviewer 1/4 request). **S7** — log10(IC50) vs. mean
AUC-fitness, one panel per drug, Pearson r + p annotated. **S8** — the same
comparison broken out per selection concentration (5 AMP panels, 7 AZT
panels).

Headline numbers (see `validation/S7_S8_provenance.md` for an independent
derivation of the same numbers from a different, minimal snippet):

| Drug | Batch    | n  | Pearson r | p        |
|------|----------|----|-----------|----------|
| AMP  | 20260407 | 13 | 0.885     | 5.9e-05  |
| AZT  | 20260307 | 18 | 0.803     | 6.1e-05  |

## Provenance

Ported from this lab's private rebuttal workspace,
`rebuttal_response/sprints/validation_figures_replot/analysis.py` (not part of
this public repository). That script was itself a manuscript-palette replot
of Ilona/Deniz's original `validation/src/fig_c_mic_vs_fitness.py` scatter
grid (Altair, per-variant qualitative colours, full 4-metric x N-column grid)
— restyled to the drug-primary palette / DejaVu Sans / open-circle-WT-DD
convention used throughout the manuscript's other extended figures, and
narrowed to just the two published panels (IC50 vs. mean fitness, and IC50
vs. per-concentration fitness).

The private script read pre-computed processed parquets from
`contributors/deniz_validation_experiments/src/processed/<batch>/xref_expanded_df.parquet`
(a private-workspace copy of the same processed validation data). This public
repository has the same table, same schema, at
`validation/src/processed/<batch>/xref_expanded_df.parquet`, which regenerates
from the raw plate-reader XLSX in `validation/data/` via
`validation/src/run_all.py` (see `validation/README.md`). No data was
re-fetched or re-derived for this port — only the plotting logic.

## What was imported

```
si_figures/s07_s08_ic50/
  README.md              (this file)
  analysis.py             new, repo-relative, runnable (see below)
  figures/                output of running analysis.py (convenience copies)
    ext_ic50_vs_auc.png            same image as figures/supplementary/figure_s07.png
    ext_ic50_per_conc_amp.png      same image as figures/supplementary/figure_s08_amp.png
    ext_ic50_per_conc_azt.png      same image as figures/supplementary/figure_s08_azt.png
```

## Path/logic ports (all flagged)

- `DENIZ_PROCESSED = PROJECT / "contributors" / "deniz_validation_experiments" / "src" / "processed"`
  -> `PROCESSED = REPO / "validation" / "src" / "processed"`, located via the
  same repo-root auto-discovery `_repo_root()` helper used by
  `analysis/s09_s10_epistatic_order/analysis.py` and
  `si_figures/s03_dose_response/analysis.py` (walks up from this file to the
  directory containing `data/processed/Epistasis_Combined.parquet`).
- Output target changed from the private sprint's own `figures/` scratch
  directory to writing **directly** into `figures/supplementary/figure_s07.png`,
  `figure_s08_amp.png`, `figure_s08_azt.png` — the actual published Supplementary
  location in this repo. A convenience copy under this directory's local
  `figures/` (using the private sprint's original filenames) is also written,
  for provenance diffing.
- Label-nudging (`_nudge_labels_2d`) stays inlined/self-contained rather than
  importing `validation/src/config_theme.nudge_labels`, so this SI script has
  no import-path dependency on the `validation/` package internals (same
  reasoning the private sprint used, now simply carried forward).
- **Dropped**: the private sprint's `build_before_after()` — a side-by-side
  QA comparison against Ilona's pre-replot `ic50_auc_dotplot.png` diagnostic
  PNG. Not ported, because (1) it is not one of the three published
  Supplementary panels, and (2) its "before" input is an internal QA image,
  not a tracked/regenerable artifact of this repository. Nothing else was
  removed — all styling, filtering, statistics, and per-panel logic for S7/S8
  are byte-for-byte the same code as the private sprint.

## Reproducing

```bash
uv run python si_figures/s07_s08_ic50/analysis.py
```

Runs in this repo's existing `uv` environment (numpy, polars, matplotlib,
scipy — all already in `pyproject.toml`; no new dependencies). Depends on
`validation/src/processed/<batch>/xref_expanded_df.parquet` existing on disk;
if missing, regenerate first with `uv run python validation/src/run_all.py`.

## Reproducibility status — VERIFIED, runs end-to-end

Ran via `uv run python si_figures/s07_s08_ic50/analysis.py` from the repo
root. Reproduces the headline statistics exactly:

```
AMP: r = 0.885, p = 5.869e-05, n = 13
AZT: r = 0.803, p = 6.127e-05, n = 18
AMP per-concentration r: 0.02 (3.1, n.s.) -> 0.67, 0.81, 0.91, 0.88 (781, ***)
AZT per-concentration r: 0.61-0.80 (all *** or **)
```

matching both `validation/S7_S8_provenance.md` and the private sprint's own
`results.md` narrative ("climbs from r = 0.02 at 3.1 µg/mL (n.s.) to r = 0.91
at 195 µg/mL", "AZT ... range 0.61-0.80") exactly.

Pixel comparison against the placeholders this code replaced (which were
themselves byte-for-byte copies of the private sprint's output, confirmed by
MD5 before this port):

- `figure_s07.png` — **byte-identical** (MD5 match) to the private sprint's
  `ext_ic50_vs_auc.png`.
- `figure_s08_amp.png` / `figure_s08_azt.png` — not byte-identical (differ by
  6-7 px out of ~4160-4165 px, i.e. ~0.15%, in the `bbox_inches="tight"`-computed
  canvas size), but a direct visual side-by-side confirms identical panel
  layout, data points, regression lines, r/p annotations, axis ranges, and
  colour ramps. The few-pixel delta is consistent with matplotlib/freetype
  text-metric drift between the environment that produced the placeholder
  (private sprint, Apr 2026) and this repo's current `.venv` (Jul 2026), not a
  logic or data difference — the printed statistics above are identical to
  the documented values in both cases.

**`figures/supplementary/figure_s07.png`, `figure_s08_amp.png`, and
`figure_s08_azt.png` now come directly from this code.**

## Caveats

- DD (catalytically dead control) has a null `mean_fitness` in the processed
  xref (excluded from the Pearson correlation) but is still plotted as an
  open circle, matching the manuscript's TEM-1$^{dead}$ reference-point
  convention — same behaviour as the private sprint, unchanged here.
