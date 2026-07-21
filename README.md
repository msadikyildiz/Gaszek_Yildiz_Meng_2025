# Higher-order epistasis in TEM-1 β-lactamase

Code and data accompanying:

> Gaszek IK\*, Yildiz MS\*, Meng Z\*, de la Paz JA, Alvarez SM, Sezer D, Morcos F, Lin M, Toprak E.
> **Higher-order epistasis drives evolutionary unpredictability toward new antibiotic resistance.**
> *Nature Communications* (in press). \*Equal contribution.

We profiled the fitness of a combinatorial mutant library of TEM-1 β-lactamase — 13
resistance-associated positions, ~55,000 genotypes — under ampicillin (AMP) and
aztreonam (AZT) selection across a concentration gradient, quantified epistasis by
three independent definitions, and modeled the resulting fitness landscapes. This
repository holds the raw and processed data, the analysis pipeline, and the notebooks
and analysis modules that regenerate the epistasis calculations, regression models, and
the data-driven figures in the paper. A few panels — schematics (BioRender), graph
layouts (Gephi), and collaborator-supplied evolutionary statistics — are produced with
external tools and are flagged as such in `figures/manifest.csv`.

---

## Repository structure

```
data/
  raw/            # Read counts + AUC per genotype (AMP/AZT) and sample metadata (CSV)
  processed/      # Tidy AUC tables (wide/long) and Epistasis_Combined.parquet
  known_variants/ # Clinical / characterized-variant encodings (encoded_variants.csv)
src/
  01_data_preprocessing.ipynb       # Raw CSV → tidy Parquet (mean / median / std / CV)
  02_epistasis_pipeline.py          # Biochemical, ensemble, regression epistasis → Epistasis_Combined.parquet
  03_epistasis_azt_regression.ipynb # LightGBM + SHAP, AZT (Fig 5B/5D)
  04_epistasis_amp_regression.ipynb # LightGBM + SHAP, AMP (Fig 5A/5C)
  05_epistasis_figures.ipynb        # Main-figure panels (Figs 3-4) + SI panel S5
  utils/                            # W/H/V/G/X matrix generation, epistasis methods, config
  graph_analysis/                   # Landscape-graph + distribution analysis; graph_builder/ = Fig 2; s18_peak_robustness/ = Fig S18 (Meng)
  figure_scripts/                   # Standalone figure scripts (Fig 3A; known-variants overlay)
  evolutionary_statistics/          # Morcos-lab DCA analysis = Figure 6 A/B/C/G/H (LGL D/E/F via morcoslab/LGL-VAE)
analysis/                 # Analyses backing Supplementary Figs S9–S17 (one self-contained module each)
validation/               # Monoculture IC50/MIC pipeline (feeds Supplementary Figs S7, S8)
si_figures/               # Reproducing code for SI figures built from external data (S2, S3, S7/S8)
figures/
  main/                   # Main-text Figures 1–6
  supplementary/          # Supplementary Figures S1–S18
  diagnostics/            # Non-manuscript regression + exploratory renders
  manifest.csv            # figure → generator map (see figures/README.md)
source_data/              # Per-figure Source Data generators (source_data/scripts/)
env/                      # Conda environment (epistasis_env.yml) + Singularity definition
pyproject.toml            # uv project (Python 3.12); uv.lock pins the full environment
DATA_README.md            # Column-level description of every data file
CITATION.cff / LICENSE    # Citation metadata; GPL-3.0
```

> `archive/`, `cache/`, `.venv/`, and per-module generator scratch (e.g. `analysis/**/figures/`,
> `si_figures/**/figures/`, `source_data/derived/`) are git-ignored and not part of the distribution.

---

## Quick start

### uv (recommended)

Requires [uv](https://docs.astral.sh/uv/getting-started/installation/).

```bash
uv sync                # creates .venv with Python 3.12 + all pinned deps
uv run jupyter lab     # launch the notebooks
```

On Linux with CUDA the CUDA 12.4 `torch` wheel is selected automatically; on macOS
(Apple Silicon) the MPS-capable wheel is used. Selection is handled in `pyproject.toml`.

### conda (legacy)

```bash
conda env create -f env/epistasis_env.yml
conda activate epistasis_env
```

---

## Reproducing the analysis

1. **Preprocess** — run `src/01_data_preprocessing.ipynb` to turn the raw CSVs into the
   tidy Parquet tables under `data/processed/`.

2. **Compute epistasis** — `python src/02_epistasis_pipeline.py`, which writes
   `data/processed/Epistasis_Combined.parquet` (biochemical definition, ensemble
   averaging, and per-order linear regression).

   > **Hardware.** The full epistasis computation is GPU-bound. A CUDA GPU with ≥40 GB
   > memory is recommended; Apple Silicon (MPS) with comparable memory works; CPU-only
   > runs but can take weeks. The pre-computed `Epistasis_Combined.parquet` is shipped
   > so the downstream steps run without re-deriving it.

3. **Model fitness** — run `src/03_epistasis_azt_regression.ipynb` and
   `src/04_epistasis_amp_regression.ipynb` for LightGBM cross-validation, learning
   curves, permutation importance, and SHAP attributions.

4. **Generate figures** — main-text panels come from `src/05_epistasis_figures.ipynb`
   (and `src/03`/`src/04` for Fig 5). Supplementary Figures are produced by the
   self-contained modules under `analysis/` (S9–S17), `validation/` (S7/S8), and
   `si_figures/` (S2, S3). See `figures/README.md` + `figures/manifest.csv` for the full
   map, and `source_data/` for the per-figure Source Data generators.

---

## Data

Genotypes are encoded as a 13-character string over the 13 mutated positions
(`.` = wild-type, an amino-acid letter = that substitution, `X` = non-functional / dead).
Positions, in protein numbering, are 19, 37, 67, 102, 162, 180, 235, 236, 237, 241,
261, 271, 272 (Ambler: 21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276).

`Epistasis_Combined.parquet` is the analysis table: one row per genotype × drug ×
concentration, carrying fitness (median AUC), its error, and epistasis under each of the
three definitions. See **[DATA_README.md](DATA_README.md)** for the full column-level
description of every raw and processed file.

---

## Data and code availability

- **Raw sequencing reads** — barcode-count Illumina libraries and the PacBio barcode-to-genotype
  matchbook — are deposited in the NCBI Sequence Read Archive; the accession is given in the
  article's Data Availability statement.
- **Processed data and analysis code** are in this repository, archived at Zenodo:
  [10.5281/zenodo.21442350](https://doi.org/10.5281/zenodo.21442350) (concept DOI, resolves to the latest version).
- **Source data** underlying the figures are provided with the published article.

---

## Figures

`figures/` collects every published figure in one place:

| Folder | Contents |
|---|---|
| `main/` | Main-text Figures 1–6 |
| `supplementary/` | Supplementary Figures S1–S18 (`figure_sNN[_amp/azt].png`) |
| `diagnostics/` | Non-manuscript renders: regression learning-curve/SHAP/permutation panels and exploratory distributions |

`figures/manifest.csv` maps every figure and panel to the code that generates it, and
`figures/README.md` explains how each is reproduced. Most regenerate from this repository's
code — `src/05_epistasis_figures.ipynb` (Figs 3–4, S5), `src/03`/`src/04` (Fig 5 SHAP),
`src/figure_scripts/figure_3a.py` (Fig 3A), `src/graph_analysis/summary-stats.ipynb`
(Fig 1E global-fitness ridgelines), the `analysis/` modules (S9–S17), and the `si_figures/`
folders (S2, S3, S7/S8). A few are produced with external tools (Gephi for the Fig 2 / S4
graph layouts, BioRender for the Fig 1 schematics) or by collaborators (Fig 6 DCA — code folded
in at `src/evolutionary_statistics/`, needs an external MSA + `py-mfdca`; S6; S18),
and Fig 1F and S1 have no in-repo generator; these cannot be
regenerated by this repository's code alone — noted in the manifest.

Supplementary Figures whose generating code needed data outside the main pipeline are
provided under `si_figures/` for this public archive:

| Folder | Produces |
|---|---|
| `si_figures/s02_mic/` | Supplementary Fig **S2** (MIC panel for TEM-1-CML and select variants across 9 β-lactams). Uses the lab's `plategig` fitting package — an **optional external dependency** pinned to its `biohpc` commit and requiring **Python ≥ 3.12** (`uv sync --extra si-figures` under a 3.12+ interpreter). This folder reproduces the figure's underlying data. |
| `si_figures/s03_dose_response/` | Supplementary Fig **S3** (dose-response profiles for a representative genotype subset). `analysis.py` reads `data/raw/*_auc_per_genotype.csv` directly. |
| `si_figures/s07_s08_ic50/` | Supplementary Figs **S7/S8** (monoculture IC50 vs AUC-fitness). Reads `validation/src/processed/`. |

See each folder's README for provenance and reproducibility notes.

---

## License

Released under the GNU General Public License v3.0. See [LICENSE](LICENSE).
