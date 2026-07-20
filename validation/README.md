# Monoculture IC50/MIC validation experiments

Independent wet-lab validation of the pooled AUC-fitness landscape: monoculture
dose-response of representative TEM-1 variants against ampicillin and aztreonam,
measured on a plate reader. These data underlie **Supplementary Figures S7 and
S8** and the correlation between pooled AUC-fitness and monoculture IC50
(Pearson r = 0.88 for AMP, r = 0.80 for AZT).

Panel sizes: **AMP n = 13** (batch 20260407), **AZT n = 18** (batch 20260307);
see `S7_S8_provenance.md` for the exact variant/batch provenance.

## Layout

```
data/                 raw plate-reader measurements + protocols
  20260220/ 20260307/ 20260407/   per-batch OD600 plate XLSX + variants.json + README
  *.heic                          lab-notebook pages (wet-lab protocol record)
  barcode_coverage_report.tsv     colony barcode-to-variant coverage
experiments.json      batch + drug + dilution configuration
variants.json         variant identity table
src/                  analysis pipeline (raw -> processed -> figures)
  process_data/       plate parsing, background subtraction, AUC, MIC/IC50 fits,
                      cross-reference to the pooled fitness landscape
  fig_*.py            figure scripts (fig_c_mic_vs_fitness.py = S7/S8)
  export_csv.py       per-variant IC50/MIC summary tables
  run_all.py          full pipeline driver
  processed/<batch>/  processed per-variant tables (IC50, MIC, AUC, dose-response)
plate_reader/         plate-reader I/O utilities (utils.py) + config notebook
```

## Reproducing

From the repository root, with the environment set up (`uv sync`):

```bash
uv run python validation_experiments/src/run_all.py            # all batches
uv run python validation_experiments/src/fig_c_mic_vs_fitness.py   # S7/S8 only
uv run python validation_experiments/src/export_csv.py            # per-variant tables
```

The pipeline reads the raw plate XLSX in `data/`, writes processed per-variant
tables into `src/processed/<batch>/`, and regenerates figures into
`validation_experiments/figures/` (gitignored; rebuilt on run). The
cross-reference step joins monoculture metrics to the pooled fitness landscape
in `../data/processed/Epistasis_Combined.parquet` (located automatically).
