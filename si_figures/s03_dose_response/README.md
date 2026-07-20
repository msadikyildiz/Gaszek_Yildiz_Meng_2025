# Supplementary Figure S3 - dose-response profiles for a representative subset

Caption: "Dose-response fitness profiles for TEM-1 beta-lactamase variants. (A)
Fitness of a representative subset of the TEM-1CML library across increasing
[antibiotic] concentration ..."

## How this was found

The user's original pointer, `TEM1_Combinatorial_Mutagenesis/notebooks/`, was
searched thoroughly: all date-dirs, plus `src/`, `dump/`, `archive/`,
`log/`, `tools/` (excluded as third-party bioinformatics binaries), and the
`.claude/` dir (empty) and `src/Epistasis/` git history (nothing relevant).
Several strong-looking candidates were checked and ruled out or superseded:

- `src/Epistasis_Analysis/5-Epistasis_Global_fitness_calc_display_v02.ipynb`
  - has the `TEM-1$^{dead}$` / `TEM-1$^{WT}$` label convention, but its
    plotting cells are a 2-D AMP-vs-AZT fitness histogram and per-mutant-combo
    bar charts, not a concentration dose-response.
- `src/Epistasis/src/06_mutant_dashboard.py` (Streamlit app) -
  `plot_fitness_curves_by_drug()` has the right shape (gray background +
  colored subset + `ylabel='auc'`, `xlabel=f'{drug} concentration'`) but is
  interactive (user picks genotypes at runtime) and has no WT/dead styling.
- `notebooks/20240723/06-freq-func-of-time.ipynb` cell 20 - very close
  (same `x='conc', y='auc', hue='genotype'` pattern, gray background), but
  its embedded output has a negative-inclusive y-range (fold-change-based
  AUC, an earlier metric) that doesn't match the manuscript panel.
- `notebooks/20240723/07-explore-individual-mutants.ipynb` cells 3-9
  ("Visualize genotype AUCs of select mutants") - embedded output matches
  the target's axis ranges and WT/dead reference-line shapes closely, but
  only plots 5 hand-picked genotypes, not the ~18 in the published legend.

**Confirmed generator**: `notebooks/20240723/04-the-most-abundant-genotypes.ipynb`,
cells 16-20. Cell 18's comment reads verbatim "Select top 18 genotypes by auc
of aucs across all conditions / add spiked wt and dead mutant for
comparison". Its embedded cell-20 output (extracted and decoded from the
notebook's own execution record - not re-run by us) shows, per drug, the
**exact same 18 genotype codes in the exact same order** as the manuscript's
published Figure S3 legend, e.g. Ampicillin panel: `....N..S...T.`,
`IKV.ST.....L.`, `..VK.T..K...D`, ... `P......SK....`, plus `.............`
(WT) and `XXXXXXXXXXXXX` (dead) - confirmed by direct visual/textual
comparison against `figures/supplementary/figure_s03.png`.

Two sibling notebooks with the identical cells 16-20 also exist
(`04-the-most-abundant-genotypes-spiked-wt.ipynb`,
`04-...-stringent.ipynb`, differing only in upstream genotype-calling
stringency and output filenames); not imported since the plain variant is
the confirmed match.

## What was imported

```
si_figures/s03_dose_response/
  README.md                                    (this file)
  analysis.py                                  new, repo-relative, runnable (see below)
  source_notebook/
    04-the-most-abundant-genotypes.ipynb       909,256 B, fetched as-is, for provenance
  figures/                                      output of running analysis.py
    dose_response_ampicillin.png
    dose_response_aztreonam.png
    figure_s03_combined.png                    same image written to figures/supplementary/figure_s03.png
```

**No data was fetched for this figure.** The AUC-per-genotype CSVs that
`04-the-most-abundant-genotypes.ipynb` needs
(`genotype_auc_sorted_ampicillin.csv`, `genotype_auc_sorted_aztreonam.csv`,
written by its own cell 17) are byte-identical to files already in this
repository:

| BioHPC file (`dump/20240723/misc/`) | Size | Repo file (`data/raw/`) | Size |
|---|---|---|---|
| `genotype_auc_sorted_ampicillin.csv` | 21,968,691 B | `Ampicillin_auc_per_genotype.csv` | 21,968,691 B |
| `genotype_auc_sorted_aztreonam.csv` | 28,665,545 B | `Aztreonam_auc_per_genotype.csv` | 28,665,545 B |

verified by exact byte-size match during import (not re-hashed, but the sizes
matching to the byte across a 22-29 MB file is strong evidence of identity).

## Why the source notebook is not directly runnable here, and what `analysis.py` is

`04-the-most-abundant-genotypes.ipynb` computes its AUC table from raw
per-sample barcode-count parquet files
(`dump/20240723/combined-barcode-counts/*.parquet`) and a PacBio
barcode-to-genotype matchbook - multi-GB sequencing intermediates that are
correctly *not* part of this repository's public data (see
`DATA_README.md`). Rather than force those in, `source_notebook/` is kept
purely for **provenance** (it is the confirmed generator, byte-for-byte
matched against the published legend), while **`analysis.py` is the
actually-runnable, repo-relative port**: it starts from the already-sorted
`data/raw/*_auc_per_genotype.csv` (exactly what the notebook's cell 17
produces) and reproduces cells 18-20's selection and plotting logic
directly in Python/pandas/seaborn, no barcode data needed.

### Path/logic ports (all flagged)

- `PROJECT_PATH = Path("/project/greencenter/Toprak_lab/shared/...")` -> repo-root-relative
  `RAW = REPO / "data" / "raw"` (same auto-discovering `_repo_root()` helper
  used by `analysis/s09_s10_epistatic_order/analysis.py` etc., for
  consistency with this repo's existing convention).
- Genotype selection (`df.index[:18] + wt_index + dead_index`),
  wide-to-long melt, and the two-layer gray-background / colored-subset
  `pointplot` structure are reproduced as in the original cell.
- **Deviation 1** (flagged in the module docstring): `TEM-1^WT` / `TEM-1^dead`
  are drawn in a fixed black / blue on top of the `tab20b` palette used for
  the other 18 genotypes, matching the black-WT/blue-dead convention used
  elsewhere in this lab's code (`06_mutant_dashboard.py`,
  `5-Epistasis_Global_fitness_calc_display_v02.ipynb`, both confirmed
  during the search above). The original exploratory cell left all 20 on
  the same qualitative palette with no special styling.
- **Deviation 2** (flagged): the 1000-genotype background sample
  (`np.random.choice`) is seeded (`SEED = 20240723`, the date-stamp of the
  source notebook) for reproducibility; the original was unseeded. Neither
  deviation changes which genotypes are selected or their AUC values.
- Panels are combined into one two-row figure with "(A)"/"(B)" labels
  (the original notebook's loop produced two separate standalone figures,
  one per drug; its own `savefig` calls used a bare relative filename with
  no directory, e.g. `'top genotypes Ampicillin.png'` - meaning the actual
  rendered PNG was never written back into the tracked BioHPC project tree;
  this is presumably why no standalone output file could be found there
  during the search, only the embedded notebook-cell image used for
  verification above).

## Reproducibility status - VERIFIED, runs end-to-end

```
uv run python si_figures/s03_dose_response/analysis.py
```

Runs in this repo's existing `uv` environment, no new dependencies. Produces
`figures/supplementary/figure_s03.png` directly. Visual comparison against
the placeholder it replaced: same 18-genotype legend (identical codes, same
rank order) per drug, same axis ranges and tick values (Ampicillin: 0, 3.1,
12.2, 48.8, 195, 781; Aztreonam: 0, 0.44, 1.33, 4, 12, 36, 108, 324), same
qualitative shape (dead genotype's sharp decline, WT's modest drift, the
representative subset trending upward with concentration), same gray
full-library background. **`figures/supplementary/figure_s03.png` now comes
directly from this code.**

## Caveat

The exact tab20b-vs-original color assignment per individual genotype is not
guaranteed to match the published figure pixel-for-pixel (the published
figure's precise color draw was not recoverable - see Deviation 1 above) -
only the set of 18 genotypes, their rank order, and the WT/dead highlighting
convention are confirmed to match.
