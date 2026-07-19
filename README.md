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
that regenerate every epistasis calculation, regression model, and figure in the paper.

---

## Repository structure

```
data/
  raw/            # Read counts + AUC per genotype (AMP/AZT) and sample metadata (CSV)
  processed/      # Tidy AUC tables (wide/long) and Epistasis_Combined.parquet
src/
  01_data_preprocessing.ipynb       # Raw CSV → tidy Parquet (mean / median / std / CV)
  02_epistasis_pipeline.py          # Biochemical, ensemble, regression epistasis → Epistasis_Combined.parquet
  03_epistasis_azt_regression.ipynb # LightGBM + SHAP, AZT
  04_epistasis_amp_regression.ipynb # LightGBM + SHAP, AMP
  05_epistasis_figures.ipynb        # Publication figures
  utils/                            # W/H/V/G/X matrix generation, epistasis methods, config
  graph_analysis/                   # Fitness-landscape graphs + distribution analysis (Fig 1, 2, S4, S18)
figures/                  # Publication figures, named by manuscript panel (see "Figures" below)
revision_analyses/        # Code reproducing the revision's new Supplementary Figures (S9–S18)
validation_experiments/   # Monoculture IC50/MIC validation data + pipeline (Supplementary Figs S7, S8)
env/                      # Conda environment (epistasis_env.yml) + Singularity definition
pyproject.toml            # uv project (Python 3.12); uv.lock pins the full environment
DATA_README.md            # Column-level description of every data file
LICENSE                   # GPL-3.0
```

> `archive/`, `cache/`, and `.venv/` are git-ignored and not part of the distribution.

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

4. **Generate figures** — run `src/05_epistasis_figures.ipynb`.

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

Files in `figures/` are named by their manuscript identity — `Figure <N><panel>. <description>.png`
for main figures 1–6 and the supplementary panels (S1, S4, S5). Supporting material lives
in subfolders:

| Folder | Contents |
|---|---|
| `amp_regression/`, `azt_regression/` | Regression diagnostics: learning curves, permutation importance, SHAP heatmaps, 2-D prediction histograms, and learning-curve statistics (CSV). |
| `known_variants/` | Clinical / characterized-variant encodings and the notebook that places them on the landscape. |
| `Figure 3A/` | Standalone script for the Figure 3A AMP-vs-AZT scatter. |
| `_not_in_manuscript/` | Notebook renders with no manuscript panel, plus `_pre_s9_s13/` — a frozen pre-revision figure snapshot under the earlier numbering. See the folder's README. |

All figures are regenerable from `src/05_epistasis_figures.ipynb`.

---

## License

Released under the GNU General Public License v3.0. See [LICENSE](LICENSE).
