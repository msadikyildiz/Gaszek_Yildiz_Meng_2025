# Figure 6 (+ Supplementary Figure S6) — evolutionary statistics (DCA / family analysis)

Direct-coupling-analysis (DCA) of the Pfam **PF13354** β-lactamase family behind
**Figure 6** and **Supplementary Figure S6**, contributed by the Morcos lab
(F. Morcos, J. A. de la Paz, S. Alvarez, UT Dallas).

- **`Figure6_Analysis.ipynb`** — the analysis notebook from J. A. de la Paz.
- **`reproduce_fig6.py`** — a faithful, repo-relative transcription of that notebook that
  runs on our shipped data and writes the Source Data. This is the in-repo entry point.

## Panel coverage

`reproduce_fig6.py` regenerates the DCA panels and their Source Data:

| Panel | What it is | Source Data |
|---|---|---|
| **6A** | PF13354 family sequence logo at the nine family positions | `source_data/fig6a_family_logo_information.csv` |
| **6B** | contextual effective alphabet (family average ± SD vs TEM-1) | `source_data/fig6b_effective_alphabet.csv` |
| **6C** | family-Hamiltonian distribution by dataset (all / top-10k / top-2.5k common) | `source_data/full/fig6c_hamiltonian_by_dataset.csv` |
| **6G** | family-Hamiltonian distribution, AMP 781 µg ml⁻¹, top-2.5k | `source_data/full/fig6g_hamiltonian_amp.csv` |
| **6H** | family-Hamiltonian distribution, AZT 36 µg ml⁻¹, top-2.5k | `source_data/full/fig6h_hamiltonian_azt.csv` |
| **S6** | family Hamiltonian vs fitness, per concentration (Spearman ρ) | `source_data/figS6_spearman_stats.csv` + `full/figS6_hamiltonian_vs_fitness_{amp,azt}.csv` |

**Panels 6D / 6E / 6F (Latent Generative Landscape) are NOT produced here.** They come
from the separate LGL-VAE model at <https://github.com/morcoslab/LGL-VAE> (J. A. de la Paz;
**MIT** licensed), archived on Dryad (**doi:10.5061/dryad.51c59zwbn**). The VAE is
stochastic, so exact 6D/E/F coordinates require the deposited model, not a re-run. The
notebook writes the LGL input MSAs (`AMP/AZT_MSA_wrtPF13354_*_headers.fasta`); S. Alvarez
assembles 6D/E/F into the final figure.

## Reproduced in-repo

```bash
uv sync --extra evo-stats          # installs py-mfdca (pinned) + pyhmmer + logomaker
export MPLCONFIGDIR=$PWD/.mpl_cache
python src/evolutionary_statistics/reproduce_fig6.py     # ~5 min (mean-field DCA + Hamiltonians)
```

It (1) filters PF13354 to ≤5 % gaps, (2) back-maps the nine family positions with an
`hmmscan` of the TEM-1 reference (using the committed `MSAs/TEM1_263AA.scan` when the
HMMER binary is absent — the reference-to-profile alignment is fixed), (3) fits the
mean-field family model, and (4) scores the family and the 55,293 intended designs,
writing the plotted values to `source_data/`.

**Verified** against the assembled figure: the panel-6B effective alphabet matches to two
decimals at all nine positions (e.g. position 104 family-mean 1.86 / TEM-1 bar 3.39;
position 240 family-mean 2.41 / TEM-1 bar 3.76). Run constants: `M = 27242` filtered family
sequences, `N = 214`, `Meff = 8377`, TEM-1 Hamiltonian −1838.24.

> **Note (n).** The 6B family average is over **n = 27242** filtered family sequences — the
> count this pipeline produces, and the value carried by the figure and the Source Data.

## Figure 6B reporting details

- **6A shows 9 of the 13 library positions** — `['69','104','164','182','237','238','240','244','265']`
  (Ambler). Positions 21, 39, 275, 276 fall outside the PF13354 mature-domain alignment (21 is
  the signal-peptide residue).
- **6B error bars are the standard deviation** across the family sequences.
- **6B `n`** is the number of filtered family sequences (27242; see the legend note above).

## Dependencies (optional extra `evo-stats`)

- **`py-mfdca`** — the `dca` package (`from dca.dca_class import dca`),
  <https://github.com/utdal/py-mfdca> (UT Dallas). Pinned to commit `24bb3adf` (v1.6.0),
  installable on Python ≥ 3.12. (The `v1.0.0` release tag is **not** used: it predates
  `compute_EffAlphabet` and requires Python < 3.10.)
- **`pyhmmer`** — bundles HMMER; no separate binary needed. The notebook's `subprocess`
  `hmmscan` path is only exercised when a system HMMER is on `PATH`.
- `logomaker`, `biopython`, `numpy`, `pandas`, `seaborn`, `matplotlib`, `scipy`.

## Inputs

- **`MSAs/PF13354_ga.fasta.gz`** — the Pfam PF13354 alignment (104 745 sequences, gzipped;
  `reproduce_fig6.py` decompresses it on first run). Provenance (J. A. de la Paz): the
  PF13354 HMM profile from <https://www.ebi.ac.uk/interpro/entry/pfam/PF13354/>, `hmmsearch`
  against UniProt TrEMBL + Swiss-Prot. See `MSAs/README.md`.
- `../../data/processed/{amp,azt}_auc_wide_df.parquet` — our AUC-fitness tables.
- `../../data/processed/TEM1-combinatorial-mutagenesis-intended.csv` — the intended-genotype list.

## Source Data layout

Small summary tables (6A logo, 6B effective alphabet, S6 Spearman, DCA constants) are
committed here. The large per-genotype tables (6C/6G/6H Hamiltonians, S6 Hamiltonian-vs-fitness;
~28 MB) are regenerated into `source_data/full/` (git-ignored) and packaged into the release
Source Data archive.
