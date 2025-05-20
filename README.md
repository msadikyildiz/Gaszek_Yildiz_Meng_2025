# TEM1 Combinatorial Mutagenesis Epistasis Analysis

This repository accompanies our study on combinatorial mutagenesis of TEM-1 β-lactamase under Ampicillin (AMP) and Aztreonam (AZT) selection. It provides raw and processed data, analysis pipelines, Jupyter notebooks, and scripts to reproduce epistasis calculations, regression modeling, and figures from the paper.

---

## 📁 Repository structure

```
data/
  raw/            # Raw CSV data (read counts, AUC, metadata)
  processed/      # Pre-processed Parquet tables and combined epistasis data
src/
  01_data_preprocessing.ipynb    # Clean and summarize raw data
  02_epistasis_pipeline.py       # Orchestrate epistasis calculations
  03_epistasis_azt_regression.ipynb  # Regression and SHAP for AZT
  04_epistasis_amp_regression.ipynb  # Regression and SHAP for AMP
  05_epistasis_figures.ipynb      # Generate publication figures
  utils/                          # Core scripts for W-matrix generation & epistasis methods
figures/                          # Output figures organized by analysis
archive/                          # Historical notebooks and experiments
env/                              # Conda environment definition
data_readme.md                    # Detailed description of data files
README.md                         # Project overview and instructions
```

---

## 🚀 Quick start

1. **Install environment**:
   ```bash
   conda env create -f env/pytorch_cuda_env.def
   conda activate pytorch_cuda_env
   ```
2. **Preprocess data** *(if needed)*:
   - Run `src/01_data_preprocessing.ipynb` to generate Parquet files under `data/processed/`.
3. **Calculate epistasis**:
   ```bash
   python src/02_epistasis_pipeline.py
   ```
   Outputs `Epistasis_Combined.parquet` in `data/processed/`.
   
   > ⚠️ **Hardware Requirements**: Epistasis calculations require substantial computational resources:
   > - NVIDIA CUDA-enabled GPU with 40GB+ memory (recommended)
   > - Apple Silicon MPS with similar memory capacity
   > - CPU-only mode is supported but may take weeks to complete
   >
   > For convenience, pre-computed epistasis matrices are provided in `data/processed/Epistasis_Combined.parquet`

4. **Run regression analyses**:
   - Launch `03_epistasis_azt_regression.ipynb` and `04_epistasis_amp_regression.ipynb` for LightGBM modeling, learning curves, and SHAP analyses.
5. **Generate figures**:
   - Execute `05_epistasis_figures.ipynb` to reproduce main publication plots.

---

## 🛠️ Key components

- **Data preprocessing**: Converts raw CSVs to tidy Parquet; computes mean, median, std, and CV.
- **Epistasis pipeline**: Uses Polars and PyTorch to compute biochemical, ensemble, and regression-based epistasis.
- **Regression modeling**: LightGBM cross-validation, learning curves, permutation importance, and SHAP explanations for feature importance.
- **Figure generation**: Custom plotting utilities for 2D histograms, heatmaps, and publication-quality layouts.

---

## 📖 Citation

Please cite our paper when using this code:
> **[Gas], [Journal], [Year].** Higher-Order Epistasis Drives the Emergence of Extended-Spectrum β-Lactamases Under Novel β-Lactam Selection.

---

## 📜 License

This project is licensed under the GNU General Public License v3.0. See `LICENSE` for details.
