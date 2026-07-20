# MSAs — sequence alignments for the Figure 6 DCA analysis

The Figure 6 notebook (`../Figure6_evolutionary_statistics.ipynb`) reads its Pfam
alignment from this directory:

- **`PF13354_ga.fasta`** — **not yet committed.** Drop it here to run the notebook.
  - Source: A. de la Paz (Morcos lab), OneDrive attachment on the 2025-11-26 email.
  - Reproducible provenance: the PF13354 HMM profile from
    <https://www.ebi.ac.uk/interpro/entry/pfam/PF13354/logo/>, run through `hmmsearch`
    against the UniProt TrEMBL + Swiss-Prot databases
    (<https://www.uniprot.org/help/downloads>).

The intermediate alignments the notebook uses (`TEM1_263AA.fasta`,
`TEM1_single_MSA.fasta`, `AMP_MSA_wrtPF13354.fasta`, `AZT_MSA_wrtPF13354.fasta`) are
generated from `PF13354_ga.fasta` plus our genotype data during a run; they do not need
to be supplied separately.
