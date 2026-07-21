# Figures not in the manuscript

These are renders from `src/05_epistasis_figures.ipynb` that do **not** map to
any numbered panel in the manuscript. They are regenerable from the notebook and
kept for reference only. The canonical published figures live under
`figures/main/` and `figures/supplementary/`.

| File | What it shows | Why it has no panel |
|---|---|---|
| `double_mutants_pairwise_epistasis_AMP.png` | 19×19 pairwise epistasis-*value* heatmap, AMP | Figure 4C delivers the pairwise epistasis of the Fig 4 title via the double-mutant *fitness* heatmap; the dedicated epistasis heatmaps across all concentrations are Supplementary Figure S11 (`analysis/s11_s12_concentration_grid`). |
| `double_mutants_pairwise_epistasis_AZT.png` | Same, AZT | Same as above. |
| `biochemical_definition_global_AZT_vs_num_mutations.png` | Global (biochemical-definition) epistasis vs number of mutated residues, AZT | Not a cited panel in any manuscript figure. |
| `fitness_distribution_AZT_AMP.png` | AMP-vs-AZT global-fitness histogram with WT/dead reference lines | Figure 3 presents the AMP-vs-AZT relationship as the scatter (Figure 3A), not as a histogram. |
