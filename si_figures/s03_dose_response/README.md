# Supplementary Figure S3 - dose-response profiles for a representative subset

Caption: "Dose-response fitness profiles for TEM-1 beta-lactamase variants. (A)
Fitness of a representative subset of the TEM-1CML library across increasing
[antibiotic] concentration ..."

## Overview

`analysis.py` reproduces Supplementary Figure S3 from the deposited
genotype-AUC tables. For each drug it selects the 18 highest-AUC genotypes,
adds the spiked wild-type (TEM-1^WT) and dead-control (TEM-1^dead) references,
draws the full library as a faint grey background, and writes the combined
two-panel figure to `figures/supplementary/figure_s03.png`.

`source_notebook/04-the-most-abundant-genotypes.ipynb` is the source notebook
defining the 18-genotype selection (cell 18, "Select top 18
genotypes by auc of aucs across all conditions / add spiked wt and dead mutant
for comparison"), and its embedded cell-20 output lists the same 18 genotype
codes, in the same order, per drug, as the published Figure S3 legend.

## Contents

```
si_figures/s03_dose_response/
  README.md                                    (this file)
  analysis.py                                  repo-relative, runnable
  source_notebook/
    04-the-most-abundant-genotypes.ipynb       source notebook defining the 18-genotype selection
  figures/                                      output of analysis.py
    dose_response_ampicillin.png
    dose_response_aztreonam.png
    figure_s03_combined.png                    same figure as figures/supplementary/figure_s03.png
```

## Inputs

`analysis.py` reads the deposited genotype-AUC tables — one row per genotype,
one column per `drug concentration replicate`, sorted descending by total AUC:

| File | Bytes | SHA-256 |
|---|---|---|
| `data/raw/Ampicillin_auc_per_genotype.csv` | 21,968,691 | `a9e4f51e67521333285445c0ae73482db789f712bfeb6e9d414f3ff32ce92e59` |
| `data/raw/Aztreonam_auc_per_genotype.csv`  | 28,665,545 | `210af43e4381f5262768727b77e31fb8ce6475e9f7ef240e73c49dbb42572e8c` |

(SHA-256 checksums are provided so a downloaded copy can be verified.)

These tables hold the same per-genotype AUC values written by the source
notebook's cell 17. The deposited CSVs and the notebook's own intermediate have
equal byte counts (21,968,691 / 28,665,545 bytes); byte identity was not
verified. The notebook derives those tables from raw per-sample barcode-count
parquet files and a PacBio barcode-to-genotype matchbook — multi-GB sequencing
intermediates that are not part of the public data deposit (see
`DATA_README.md`) — so `analysis.py` starts from the deposited CSVs directly and
needs no barcode data.

## Methods

- Selection: the top 18 genotypes by total AUC (`df.index[:18]`) plus the WT
  (`.............`) and dead (`XXXXXXXXXXXXX`) rows, mirroring the source
  notebook's cell 18.
- Each panel plots a 1000-genotype random background sample (grey), the 18
  selected genotypes on a `tab20b` qualitative palette, and the WT / dead
  references drawn last in fixed black / blue.
- The background sample is drawn with a fixed seed (`SEED = 20240723`) for
  reproducibility.

## Reproduction

```bash
uv run python si_figures/s03_dose_response/analysis.py
```

Runs in the repo's `uv` environment with no additional dependencies. Writes
`figures/supplementary/figure_s03.png` and per-drug copies under `figures/`.
The selected genotypes (identical codes and rank order per drug), axis ranges,
and tick values (Ampicillin: 0, 3.1, 12.2, 48.8, 195, 781; Aztreonam: 0, 0.44,
1.33, 4, 12, 36, 108, 324) match the published Figure S3.

## Limitation

The per-genotype colour assignment on the `tab20b` palette is not guaranteed to
match the published figure exactly; the reproduced quantities are the set of 18
genotypes, their rank order, and the WT / dead highlighting convention.
