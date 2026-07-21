# s17 — dose-response profiles with low-count flagging (Supplementary Fig. S17)

Produces **Supplementary Figure S17**: dose-response profiles for a curated
set of highlighted genotypes, with read-count-based flagging of low-count
(extinction) replicates.

## Rerun

```
uv run python analysis/s17_dose_response_low_count/analysis.py
```

Runtime: ~45 s.

## Inputs

- `data/raw/Ampicillin_auc_per_genotype.csv`,
  `Aztreonam_auc_per_genotype.csv` — per-genotype log10(AUC) across 3
  biological replicates at 6 (AMP) or 8 (AZT) concentrations. Rows = 70,183
  (AMP) / 70,572 (AZT) genotypes. Schema: `mut_profile_masked, <Drug> <conc> <rep>`.
- `data/raw/Ampicillin_read_counts_per_genotype.csv`,
  `Aztreonam_read_counts_per_genotype.csv` — per-genotype Illumina read counts
  at every drug × conc × replicate × timepoint sample. Used to flag
  low-count (extinction) events. Column schema: `<DRUG>_<conc>_<rep>_<tp>_S<n>`.
- `data/raw/metadata.csv` — maps sample column names
  to drug / concentration / replicate / timepoint. Concentration read as a
  string because "3,100" is European decimal for 3.1; parsed with
  `str.replace(",", ".").cast(Float64)`.
- `data/known_variants/encoded_variants.csv` — 63
  clinical TEM alleles with their 13-letter encoded sequences. Converted
  from residue strings to the library's masked encoding via
  `"".join('.' if seq[i] == WT[i] else seq[i])`.

## Outputs

The per-drug panels are the published Supplementary Fig. S17 and are written
directly to the shared gallery `figures/supplementary/`. The combined
two-panel figure is secondary/diagnostic and regenerates into this module's
local `figures/` (gitignored, scratch).

- `figures/fig_s17.png` — combined two-panel local preview.
  Secondary/diagnostic — regenerates locally (gitignored).
- `figures/supplementary/figure_s17_amp.png`, `figures/supplementary/figure_s17_azt.png`
  — per-drug high-resolution panels (kept in case individual embedding is
  preferred). **Published Supplementary Fig. S17**.
- `data/highlighted_variants.csv` — 107 rows: one per highlighted genotype
  with its category (WT / dead / single / clinical / top-composite), label,
  Ambler position (singles only) and ranking (composites only).
- `data/per_variant_fitness.csv` — 1,498 rows: highlighted-genotype × drug ×
  concentration stats including per-replicate AUC values, per-replicate
  `sum_reads` across 4 timepoints, low-count flags and mean / SD (both raw
  and "clean" = non-flagged-reps-only).
- `results_table.csv` — library-wide extinction census: 14 rows (one per
  drug × conc) with the fraction of genotypes whose read counts fall below
  the threshold in any / majority / all replicates.

## Highlighted-variant selection

| Category                        |  n |
|---------------------------------|---:|
| WT (TEM-1, `.............`)     |  1 |
| TEM-1_dead (`XXXXXXXXXXXXX`)    |  1 |
| Single substitutions (one-hot)  | 18 |
| Clinical TEM alleles            | 47 |
| Top-20 composites at max AMP    | 20 |
| Top-20 composites at max AZT    | 20 |
| **Total**                       |**107**|

- **Singles**: one line per non-WT residue option — 13 positions with
  1–3 non-WT choices each = 18 distinct single-substitution genotypes. Each
  single is colored by its Ambler position (shared across the 13-color
  palette `POSITION_CMAP`).
- **Clinical**: 63 clinical TEM alleles from the Lahey database, encoded in
  `encoded_variants.csv`. After removing entries that
  collide with WT (TEM-1 itself) or with one of the 18 single-mutant
  highlights, 47 unique compound clinical alleles remain. Drawn in yellow,
  matching the yellow-triangle coloring of clinical isolates in Fig. 3A of
  the main manuscript.
- **Top composites**: per drug, the 20 genotypes with the highest mean
  log10(AUC) at the highest concentration in the library (AMP 781 µg/mL,
  AZT 324 µg/mL), ≥ 2 non-WT positions, not already highlighted. Intended
  to expose the "winners" of the landscape. AMP composites drawn red,
  AZT composites drawn blue.

## Low-count (extinction) flagging

A genotype × concentration × replicate point is **low-count** if the sum of
Illumina read counts across its four timepoints (3h, 6h, 9h, 12h) is less
than **10**. The threshold was chosen before any plotting, from the clear
shoulder in the library-wide read-count distribution. The library-wide
fraction of low-count
replicates varies between 17% (AMP 0 µg/mL, baseline) and 34% (AMP 781
µg/mL), with the sharp increase at the highest AMP concentration reflecting
the drug-induced extinction discussed under LOD band, below.

## LOD band

Shaded region at-or-below the TEM-1_dead curve per concentration. A genotype
with fitness inside this band is statistically indistinguishable from a
catalytically dead enzyme, marking the empirical limit of detection of the
assay. On AMP, the LOD drops sharply from ~3.2 (0 µg/mL) to ~0.7 (781 µg/mL)
as drug-induced killing shrinks non-resistant populations. On AZT, the LOD
drops only modestly (~3.2 → 1.9) because AZT is bacteriostatic: dead-enzyme
cells still contribute some OD signal.

## Design choices

| Choice | Value | Rationale |
|--------|-------|-----------|
| Threshold | 10 reads | Library-wide distribution shows a clear shoulder; fixes ~20% baseline floor across most conditions (= sequencer noise + library under-sampling), spike at AMP 781 diagnoses extinction. Per-conc frequencies are reported in `results_table.csv`. |
| Top-composite criterion | mean log10(AUC) at highest conc × 3 reps | Picks adaptive "winners" at the selection edge. Alternative (AUC of AUCs across concs) would overweight baseline growth and erase the selection effect. |
| Composite-count | 20 per drug | Dense enough to show the upper-envelope; sparse enough not to overwhelm. |
| Clinical overlap handling | WT + singles removed from clinical set | Each highlight appears exactly once. |
| 0-conc placement on log-x | proxy tick at `lowest_nonzero / 3` | Real "0" breaks log scale; dashed grey line separates proxy from first real conc. |
| LOD visualisation | shaded region below TEM-1_dead | Reflects that the LOD is conc-dependent, not a single horizontal threshold. |

## Reproducibility

- `SEED = 20260420`.
- All ranking operations are deterministic (stable sort keyed on
  `mean_high_conc` descending); no random subsampling.
- `pdf.fonttype = 42` / `ps.fonttype = 42` so the figure exports with
  editable text.
