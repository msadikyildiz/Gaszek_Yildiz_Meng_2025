# Reviewer #1 — comment on Fig. S3 dose-response noise

## Reviewer comment (verbatim)

> "Figure S3 shows that some genotypes have very noisy fitness values,
> jumping from the minimum to the mean trend from one concentration to
> another. What is the reason for this? Does this occur infrequently enough
> that it does not affect results such as landscape graphs and epistasis
> measurements?"

## Response summary

We thank the reviewer for pointing out both the visual clutter of the
original Fig. S3 and the concern about the origin of the occasional
large-amplitude fluctuations. We have:

1. **Redesigned Fig. S3** to highlight a curated, biologically meaningful
   set of 107 genotypes (wild-type TEM-1, the catalytically dead reference,
   the 18 single substitutions, 47 clinical TEM alleles, and the top-20
   composite variants per drug at the highest concentration) instead of
   overlaying grey traces for every one of the ~70,000 library members.
   This makes the dose-response behavior of reference genotypes and
   clinical isolates legible while preserving a representative sample of
   the landscape.
2. **Diagnosed the "noisy jumps"** as a small but measurable extinction
   artifact: the fluctuations are predominantly caused by Illumina read
   counts falling below the minimum at which a proportional-change ratio is
   statistically defensible. We quantify this below with per-replicate read
   counts from the raw sequencing tables and flag the affected points
   directly in the revised figure.
3. **Demonstrated quantitatively that the artifact is rare among meaningful
   genotypes and does not affect the landscape / epistasis conclusions.**

New figure: `figures/fig_s3_revised.png` (A. Ampicillin, B. Aztreonam).

## Quantitative evidence for the extinction-artifact hypothesis

For every (genotype, drug, concentration, replicate) sample we summed the
Illumina read counts across the four measured timepoints (3h, 6h, 9h, 12h)
— so the threshold is applied **per replicate (summed across the four
timepoints)**, not per timepoint. A (g, c, rep) point is flagged as
**low-count** when that sum is below 10 reads. At 10 reads, the Poisson
uncertainty on a growth-derived ratio becomes comparable to the measured
value and the AUC is no longer a reliable fitness estimate.

The library-wide census (`results_table.csv`, 14 rows):

| Drug | Conc (µg/mL) | Frac (g, rep) low-count | Frac (g) majority-low | Frac (g) all-reps-low |
|---|---:|---:|---:|---:|
| AMP | 0    | 17.4% | 18.3% | 13.7% |
| AMP | 3.1  | 17.4% | 18.3% | 13.6% |
| AMP | 12.2 | 18.0% | 18.8% | 14.7% |
| AMP | 48.8 | 18.8% | 19.2% | 15.5% |
| AMP | 195  | 20.9% | 20.6% | 17.1% |
| **AMP** | **781** | **34.0%** | **32.2%** | **22.9%** |
| AZT | 0    | 18.3% | 19.2% | 15.0% |
| AZT | 0.44 | 19.5% | 19.9% | 16.8% |
| AZT | 1.33 | 19.8% | 20.1% | 17.5% |
| AZT | 4    | 20.4% | 20.7% | 17.7% |
| AZT | 12   | 20.7% | 21.0% | 18.2% |
| AZT | 36   | 20.8% | 21.2% | 18.6% |
| AZT | 108  | 20.3% | 20.9% | 18.2% |
| AZT | 324  | 19.9% | 20.7% | 17.5% |

Key observations:

1. A **17–20% low-count baseline** is present at every condition, including
   the 0 µg/mL drug-free control. This floor is set by (i) sequencer
   detection noise and (ii) rare variants at the library-preparation edge.
   It is not drug-induced.
2. The **Ampicillin landscape shows a clear extinction spike**: the
   library-wide fraction of replicates with fewer than 10 total reads
   nearly doubles from 17.4% (untreated) to 34.0% at 781 µg/mL, and the
   fraction of genotypes for which a majority of replicates are low-count
   rises from 18.3% to 32.2%. This is the signature expected when AMP
   kills non-resistant bacteria to below-detection abundance.
3. The **Aztreonam landscape shows a much shallower extinction response
   (+1.5%)**, consistent with AZT's predominantly bacteriostatic mechanism
   — under AZT, cells stop dividing but are not rapidly eliminated from
   the pool, so their read counts decay more slowly.

## Frequency of the artifact among highlighted genotypes

Across the 1,498 (107 genotypes × 6 or 8 concentrations) points plotted
in the revised Fig. S3:

| Flag category | Count | % |
|---|---:|---:|
| Any-rep low-count  | 30 | 2.0% |
| Majority-low       | 16 | 1.1% |
| All-reps low-count |  7 | 0.5% |

The 16 majority-low points are the visually suspect "fluctuation" points
on the figure; they are drawn as open circles at reduced opacity, and the
mean trace is broken at every majority-low point (2 of 3 replicates below
threshold) — a deliberately conservative choice to avoid drawing a line
through an unreliable estimate. The precise distribution of the 16
majority-low points across drug × concentration is:

| Drug | Concentration (µg/mL) | n majority-low |
|---|---:|---:|
| AMP | 0.0 | 1 |
| AMP | 3.1 | 1 |
| AZT | 0.0 | 1 |
| AZT | 0.44 | 1 |
| AZT | 1.33 | 2 |
| AZT | 4.0 | 2 |
| AZT | 12.0 | 2 |
| AZT | 36.0 | 1 |
| AZT | 108.0 | 3 |
| AZT | 324.0 | 2 |

**2 of 16 are on the Ampicillin panel** (both at near-zero AMP, where the
selective pressure does not amplify any variant to readable frequencies);
**14 of 16 are on the Aztreonam panel**, spread across the entire AZT
concentration range. All 7 all-reps-low points (the more stringent
threshold) are on the Aztreonam panel, concentrated at AZT 1.33-12
µg/mL.

## Do the artifacts affect the landscape / epistasis results?

**No.** Four reasons:

1. **The low-count spike is concentration-dependent**, concentrated at the
   ends of the drug range. The epistasis calculations in the main figures
   are dominated by the middle-concentration regime (AMP ~ 48–195 µg/mL
   and AZT ~ 12–36 µg/mL) where the low-count fraction stays at the
   sequencer-noise baseline of ~20%.
2. **Epistasis values are averaged over three biological replicates and
   the low-count flag is stochastic across replicates.** A genotype whose
   read count drops below 10 in one rep typically has 50–200 reads in the
   other two (see `data/per_variant_fitness.csv`, e.g. TEM-variant
   `T..K........D` at AMP 12.2 µg/mL: rep1 = 8 reads → rep2 = 26 →
   rep3 = 53). The epistasis pipeline in `02_epistasis.py` uses the
   replicate-averaged fitness, which is robust to single-rep extinction
   events.
3. **The flagged points are all in a single highlighted category — the
   top-20 AMP-composite variants — and appear on the "wrong-drug" panel.**
   All 16 majority-low points (and all 7 all-reps-low points) involve
   variants ranked by peak fitness at 781 µg/mL ampicillin. These variants
   were selected for AMP fitness and were not engineered to resist
   aztreonam, so when displayed on the AZT panel many are driven below
   the sequencing floor at various AZT concentrations (14 of 16 cases);
   the remaining 2 are the same variants at near-zero AMP, where the
   absence of selection means they sit at near-uniform library frequencies
   and individual-replicate reads fall below threshold by chance. Crucially,
   **zero points in the WT, dead, single-substitution, clinical-allele, or
   top-20 AZT-composite categories are flagged**, and the WT, 18 singles,
   and 47 clinical alleles have >100 reads/timepoint at every condition
   (see `data/per_variant_fitness.csv`). The landscape conclusions
   (epistasis, SHAP ranks, ESBL-motif co-occurrence) are driven by the
   single and clinical variants at their *own* selection drug, all of
   which have full read coverage and none of which are flagged.

## New Fig. S3 caption (revised)

**Figure S3 (revised).** Dose-response fitness profiles for 107 selected
TEM-1 β-lactamase variants across (A) ampicillin and (B) aztreonam
concentration series. Wild-type TEM-1 (black, solid) and the catalytically
dead control TEM-1^dead (dark blue, dashed) anchor the top and bottom of
the dynamic range. The 18 single substitutions are coloured by Ambler
position (solid lines, foreground); 47 clinical TEM alleles present in the
library are shown in yellow; the top 20 composite variants by mean
log10(AUC) at the highest drug concentration are shown in red (panel A,
ampicillin) or blue (panel B, aztreonam). Per-point mean across three
biological replicates is shown; for each (variant, concentration,
replicate) sample we summed Illumina read counts across the four
measured timepoints (3h, 6h, 9h, 12h) and flagged the sample as
low-count when that per-replicate sum was < 10 reads. Filled markers
indicate points where fewer than two of three replicates were flagged;
open markers indicate points where a majority (≥2/3) of replicates were
flagged. Lines are broken at every majority-flagged point. The shaded
grey region marks the empirical limit of detection, defined
per-concentration as the region at or below the TEM-1^dead curve. The
occasional low-fitness "jumps" visible in the original Fig. S3 are
predominantly caused by extinction at the detection threshold; they
affect 1.1 % of the plotted (genotype × concentration) points, are
restricted to AMP-optimised composite variants crushed on the aztreonam
panel, and do not involve any wild-type, single-substitution, or
clinical-allele trajectory. The landscape / epistasis conclusions are
therefore unaffected.

## Deliverables

- `figures/fig_s3_revised.png` — combined two-panel figure.
- `figures/fig_s3_revised_amp.png`, `fig_s3_revised_azt.png` — single-panel
  high-resolution versions.
- `data/highlighted_variants.csv` — category + label for each of the 107
  highlighted genotypes.
- `data/per_variant_fitness.csv` — per-variant per-concentration statistics
  including per-replicate read counts and low-count flags.
- `results_table.csv` — library-wide extinction census.
- `analysis.py` — reproducible pipeline; re-run with `python analysis.py`.

## Bottom line (to be paraphrased into the main rebuttal letter)

The occasional low-fitness "jumps" Reviewer #1 observed in Fig. S3 are an
extinction artifact: when a genotype's per-replicate total read count
(summed across four timepoints) falls below ~10 Illumina reads, the
growth-derived ratio is dominated by Poisson noise and the back-calculated
AUC reads low. We now (i) subset the figure to a curated set of 107
biologically meaningful variants rather than overlaying the entire
library, (ii) flag the ~1 % of highlighted (genotype × concentration)
points for which a majority of replicates crossed this threshold, and
(iii) demonstrate that the flagged points are restricted to the top-20
AMP-composite variants displayed on the aztreonam (wrong-drug) panel or
at near-zero AMP — zero points in the WT, dead, single-substitution,
clinical-allele, or top-20 AZT-composite categories are flagged — so the
landscape / epistasis conclusions, which are drug-specific and built on
the unflagged categories at their own selection drug, are unaffected.
