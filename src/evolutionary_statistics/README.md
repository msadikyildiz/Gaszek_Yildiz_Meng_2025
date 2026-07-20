# Figure 6 — evolutionary statistics (DCA / family analysis)

Direct-coupling-analysis (DCA) code behind **Figure 6**, contributed by the Morcos
lab (F. Morcos, A. de la Paz, S. Alvarez, UT Dallas). `Figure6_evolutionary_statistics.ipynb`
is the notebook **as received** from A. de la Paz — it is internally titled *"Figure 8
and Supplementary Figure 8"* because it was written against the original-submission
numbering. **Original Figure 8 = revision Figure 6.**

## Panel coverage

This notebook reproduces the DCA half of Figure 6:

| Figure 6 panel | Notebook output (`_8x` = original numbering) | What it is |
|---|---|---|
| **6A** | `FamilyLogo_8a.png/pdf` | PF13354 family sequence logo (logomaker) |
| **6B** | `TEM1_EffAlph_fam_ga_SpecificPositions_8b.png/pdf` | effective alphabet at the target positions |
| **6C** | `Ridge_Histogram_H_fam_by_Dataset_8C.png/pdf` | ridge histogram of the family Hamiltonian by dataset |
| **6G** | `H_fam_Distribution_withTop2.5k_781AMP_G.png/pdf` | family-Hamiltonian distribution, AMP 781 |
| **6H** | `H_fam_Distribution_withTop2.5k_36AZT_H.png/pdf` | family-Hamiltonian distribution, AZT 36 |

**Panels 6D / 6E / 6F (LGL coordinates) are NOT produced here** — they come from the
separate LGL-VAE model at <https://github.com/morcoslab/LGL-VAE> (A. de la Paz).

## What this settles for the legend / reporting checklist

- **6A shows 9 of the 13 library positions**, not all 13. The panel labels exactly
  `['69','104','164','182','237','238','240','244','265']` (Ambler); positions 21, 39,
  275 and 276 are absent (21 is the signal-peptide residue; the other three fall outside
  the PF13354 mature-domain alignment).
- **6B error bars are the standard deviation** across the family sequences
  (`yerr = OutputEffAlph_fam[:, positions].std(axis=0)`).
- **6B `n` is the effective number of sequences** `Meff`, printed by the notebook
  (`familyModel.Meff`, alongside `M` and `N`). It is emitted on running the notebook;
  it is not hard-coded here.

## Dependencies

- **`py-mfdca`** — the `dca` package (`from dca.dca_class import dca`), published at
  <https://github.com/utdal/py-mfdca> (UT Dallas). Not yet declared in this repo's
  `pyproject.toml` (see "Status" below).
- `logomaker`, `biopython` (`Bio`), `numpy`, `pandas`, `seaborn`, `matplotlib` — already
  in this repo's environment.
- **HMMER** (`hmmsearch` / `hmmalign`, called via `subprocess`) — a system tool, not a
  Python package.

## Inputs

Relative to the notebook's run directory (`MSADir=MSAs/`, `UTSWDir=data/`, `CleanFigDir=figures/`):

- **`MSAs/PF13354_ga.fasta`** — the Pfam PF13354 alignment. **NOT YET IN THIS REPO** — see
  `MSAs/README.md`. Provenance (A. de la Paz): the PF13354 HMM profile from
  <https://www.ebi.ac.uk/interpro/entry/pfam/PF13354/logo/>, `hmmsearch` against the
  UniProt TrEMBL + Swiss-Prot databases (<https://www.uniprot.org/help/downloads>).
- `data/amp_auc_wide_df.parquet`, `data/azt_auc_wide_df.parquet` — our AUC-fitness tables
  (present in this repo at **`data/processed/`**; the notebook's `UTSWDir` must be pointed
  there when run).
- `data/TEM1-combinatorial-mutagenesis-intended.csv` — the intended-genotype list.
- Intermediate MSAs (`TEM1_263AA.fasta`, `TEM1_single_MSA.fasta`,
  `AMP_MSA_wrtPF13354.fasta`, `AZT_MSA_wrtPF13354.fasta`) — built by the notebook from the
  above during a run.

## Status — folded in as received; not yet runnable end-to-end

The notebook is committed verbatim for provenance. Three things are needed before it runs
top-to-bottom in this repo and its Source Data numbers (6B effective alphabet; 6C/6G/6H
Hamiltonians) can be extracted:

1. **`MSAs/PF13354_ga.fasta`** — obtain from A. de la Paz (the OneDrive attachment on the
   2025-11-26 email) or regenerate from the InterPro/UniProt provenance above.
2. **`py-mfdca`** installed (then pin it as an optional `pyproject` extra, mirroring
   `plategig`, and re-lock).
3. Point `UTSWDir` at `../../data/processed/` (a one-line path port).

The **LGL panels (6D/E/F)** additionally require the `LGL-VAE` repo and its inputs.

Once (1)-(3) land, this becomes a fully reproducible in-repo generator for Figure 6A/B/C/G/H
(like `si_figures/s02_mic/`), and the Source Data sheet for Figure 6 can be populated.
