# MSAs — sequence alignments for the Figure 6 DCA analysis

Inputs to `../reproduce_fig6.py` (and `../Figure6_Analysis.ipynb`).

- **`PF13354_ga.fasta.gz`** — the Pfam PF13354 alignment, 104 745 sequences (gzipped,
  ~11 MB; `reproduce_fig6.py` decompresses it to `PF13354_ga.fasta` on first run).
  - Source: J. A. de la Paz (Morcos lab).
  - Reproducible provenance: the PF13354 HMM profile from
    <https://www.ebi.ac.uk/interpro/entry/pfam/PF13354/>, run through `hmmsearch` against
    the UniProt TrEMBL + Swiss-Prot databases (<https://www.uniprot.org/help/downloads>).
- **`PF13354.hmm`** — the PF13354 profile HMM (for the `hmmscan` back-mapping step).
- **`TEM1_263AA.fasta`** — the TEM-1 mature-domain reference sequence.
- **`TEM1_263AA.scan`** — cached `hmmscan` output of `TEM1_263AA.fasta` against `PF13354.hmm`.
  Deterministic; used to back-map the nine family positions when a system HMMER binary is
  not on `PATH`. `TEM1_263AA.align` is the parsed alignment derived from it.

Intermediate alignments the run produces from the above (`PF13354_ga_filtered_5pct.fasta`,
`TEM1_single_MSA.fasta`, `AMP/AZT_MSA_wrtPF13354*.fasta`) are git-ignored; the LGL input
MSAs (`*_headers.fasta`) feed the separate LGL-VAE model for panels 6D/E/F.
