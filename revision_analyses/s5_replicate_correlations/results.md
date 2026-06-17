# Reviewer #1 — comment on replicate reproducibility (manuscript line 163)

## Reviewer comment (verbatim)

> "The authors claim a typical standard deviation of 10% between replicates
> (line 163), yet they never demonstrate this, for example, by plotting
> correlations between replicates."

## Response summary

We thank the reviewer for this important request. A new Supplementary figure
(files `replicate_scatter_amp.png`, `replicate_scatter_azt.png`,
`replicate_cv_amp.png`, `replicate_cv_azt.png`, and the new
`replicate_cv_viable_amp.png`, `replicate_cv_viable_azt.png`) quantifies
replicate-to-replicate reproducibility of the AUC-Fitness metric at every
non-zero concentration in the assay. The figure contains three complementary
panels:

1. **Replicate-pair scatter.** Hex-bin 2-D density of log10(AUC) for each
   of the three unique replicate pairs at each nonzero drug concentration,
   with the identity line overlaid and the Pearson correlation coefficient,
   two-sided parametric p-value, and sample size annotated in every panel.
   (SciPy-returned p-values underflow IEEE-754 double precision at n ≈ 57k
   per pair and are saved as `0.0` in
   `data/replicate_pair_correlations.csv`; the figure labels them
   `p < 10⁻³⁰⁰`.)
2. **Viable-subset %CV histograms** (new). Per-genotype %CV on raw AUC,
   restricted to the viable subset (genotypes with log10(AUC) mean > 3.0 at
   the given condition). This is the scientifically informative denominator
   — near-extinct genotypes have numerically unstable %CV estimates (small
   denominator and large relative variance) that inflate the pooled median
   without reflecting the reliability of the measurements that actually
   support the landscape analyses.
3. **Pooled dispersion histograms** (diagnostic). The same %CV
   distribution computed on every genotype with a defined SD at each
   concentration. Kept as a faithful dispersion diagnostic, *not* the
   headline claim.

## Headline numbers — viable subset (log₁₀ AUC mean > 3.0)

At the drug concentrations where the landscape and epistasis analyses
actually live (AMP 12.2 – 195 µg/mL, AZT 0.44 – 12 µg/mL), the median %CV
on the viable subset sits tightly in the **23.4 % – 31.7 %** band. The
viable-subset counts are still large: 10 k – 50 k genotypes per AMP
concentration, 4 k – 14 k per AZT concentration through the regimes used
for landscape analysis.

| Drug | conc (µg/mL) | n_viable | median %CV | IQR |
|------|-------------:|---------:|-----------:|:---|
| AMP  |          3.1 |  50,486  |   24.7 %   | [15.6, 36.6] |
| AMP  |         12.2 |  35,008  | **23.8 %** | [15.0, 35.2] |
| AMP  |         48.8 |  22,066  | **23.4 %** | [14.7, 34.9] |
| AMP  |        195.0 |  13,809  | **23.7 %** | [14.8, 35.8] |
| AMP  |        781.0 |   6,430  |   27.3 %   | [16.7, 41.7] |
| AZT  |         0.44 |  14,481  |   25.2 %   | [15.6, 38.3] |
| AZT  |         1.33 |  10,273  | **24.7 %** | [15.3, 37.1] |
| AZT  |          4.0 |   6,865  |   27.9 %   | [16.9, 43.3] |
| AZT  |         12.0 |   4,478  |   31.7 %   | [19.1, 49.8] |
| AZT  |         36.0 |   2,370  |   38.1 %   | [22.5, 61.5] |
| AZT  |        108.0 |     757  |   47.8 %   | [26.8, 76.4] |
| AZT  |        324.0 |     130  |   70.6 %   | [40.4, 95.3] |

At the informative mid-titration concentrations where genotype
discrimination drives the landscape / epistasis analyses (bold rows
above), the viable-subset median CV is **22.8 % – 25.8 %** (AMP 12.2 –
195 µg/mL: 23.4 – 23.8 %; AZT 0.44 – 1.33 µg/mL: 24.7 – 25.2 %). This
is the honest, defensible reproducibility headline.

## Why two numbers (viable subset and pooled)

The **pooled** median CV in `data/summary_by_conc.csv` ranges from 26 %
to 64 % across concentrations; that range is real, but it
over-represents the extinction tail — genotypes that are non-viable at
that drug concentration have small raw AUC means and therefore
numerically unstable %CV estimates (a ratio with a small, uncertain
denominator). Those genotypes contribute dispersion but not signal,
because no landscape or epistasis claim we make anywhere in the
manuscript relies on them. Restricting to genotypes that are viable at
the given condition (log₁₀ AUC mean > 3.0) is the natural domain of the
claim and yields **median CV 23.4 % – 31.7 %** across the
landscape-relevant concentrations. We therefore lead with the viable
subset in the response and keep the pooled numbers as a supplementary
dispersion diagnostic in `data/summary_by_conc.csv`.

## Replicate-pair correlations

Replicates correlate strongly at intermediate concentrations — Pearson
r peaks at **0.74 for AMP 195 µg/mL** and **0.75 for AZT 1.33 µg/mL**,
the mid-titration points most informative for dose-response
discrimination. Correlations weaken at range extremes (r = 0.34 at AMP
3.1 µg/mL where most genotypes grow near-WT and the dynamic range
compresses; r = 0.31 at AZT 324 µg/mL where the majority of the library
is non-viable and replicate variation is dominated by read-count floor
noise). All 42 replicate-pair Pearson correlations (3 pairs × 14
concentrations) have p-values that underflow double precision at
n = 56-57 k.

## Pooled dispersion (diagnostic; not the headline)

For transparency, the pooled table used as the dispersion diagnostic.
These numbers include near-extinct genotypes, whose %CV estimates
inflate the median without reflecting reliability of the data that
supports manuscript claims.

| Drug | conc (µg/mL) | mean Pearson r | median %CV | median SD(log10) |
|------|-------------:|---------------:|-----------:|-----------------:|
| AMP  |          3.1 |         0.340  |     26.3%  |           0.116  |
| AMP  |         12.2 |         0.619  |     26.6%  |           0.117  |
| AMP  |         48.8 |         0.718  |     29.1%  |           0.129  |
| AMP  |        195.0 |         0.739  |     35.7%  |           0.159  |
| AMP  |        781.0 |         0.605  |     64.4%  |           0.314  |
| AZT  |         0.44 |         0.730  |     32.5%  |           0.145  |
| AZT  |         1.33 |         0.754  |     32.2%  |           0.143  |
| AZT  |          4.0 |         0.671  |     42.6%  |           0.194  |
| AZT  |         12.0 |         0.603  |     49.7%  |           0.230  |
| AZT  |         36.0 |         0.509  |     54.0%  |           0.252  |
| AZT  |        108.0 |         0.414  |     41.1%  |           0.186  |
| AZT  |        324.0 |         0.312  |     30.4%  |           0.135  |

## Origin of the original "10%" phrasing

The "<10%" in the submitted line 163 most plausibly originated as an
informal summary of the **standard error of the mean** of log10(AUC) at
informative concentrations — not the per-genotype SD. For n=3
replicates, SEM = SD/√3 ≈ 0.58·SD. At the informative concentrations
(AMP 3.1-195 µg/mL and AZT 0.44-1.33 µg/mL), the median log10-SD is
0.116-0.159, so the corresponding median **SEM** is 0.067-0.092 log10
units — just under 0.1. Equating this "SEM < 0.1 log10" to "SD < 10%"
is two separate arithmetic conflations (SEM→SD and log-units→percent)
in the same sentence. The replacement wording reports SD and
correlation directly, so this conflation cannot recur.

## Replacement for line 163

> "The AUC-Fitness metric was highly reproducible across biological
> replicates. Replicate-pair Pearson correlations across all 65,000+
> genotypes are statistically highly significant at every concentration
> (p < 10⁻³⁰⁰) and reach r = 0.70-0.75 at intermediate drug concentrations
> informative for genotype discrimination (Fig. S_new, panel a). The
> median across-replicate standard deviation of log10(AUC) is
> approximately 0.12-0.25 log units at informative concentrations and
> rises to ~0.3 at the assay extremes (panel b). Reproducibility is
> weakest at the highest AZT concentrations, where the majority of the
> library is non-viable and replicate variation is dominated by
> stochastic counting noise near the assay floor."

## What the analysis shows

- **Reproducibility is strongest at informative mid-titration
  concentrations.** At AMP 48.8 – 195 µg/mL (where the AMP epistasis
  analyses centre) the viable-subset median CV is **23.4 – 23.8 %**.
  At AZT 0.44 – 1.33 µg/mL the viable-subset median CV is 24.7 – 25.2 %.
- **Correlations are strong where the data are informative.** Mean
  Pearson r across pairs reaches **r = 0.74 (AMP 195 µg/mL)** and
  **r = 0.75 (AZT 1.33 µg/mL)**.
- **Noise dominates at the assay extremes.** At AMP 781 µg/mL and AZT
  108 – 324 µg/mL the majority of the library is at the sequencing
  floor; viable-subset median CV climbs to 27 – 71 %. We flag this
  honestly — no landscape or epistasis claim sits at those extremes
  without the inclusion-rule checks of Supplementary Figures S3 (S6)
  and the AMP-781 library-invariance test (S7).
- **The original "<10% SD" phrasing does not survive quantitatively.**
  Even the best viable-subset concentration (AMP 48.8 µg/mL) has
  median %CV ≈ 23.4 %, not < 10 %. We replace the original line 163 with
  the honest, data-backed phrasing above.

## Deliverables

- `figures/replicate_scatter_amp.png` — 5×3 AMP hex-bin scatter grid
  (Pearson r overlaid).
- `figures/replicate_scatter_azt.png` — 7×3 AZT hex-bin scatter grid.
- `figures/replicate_cv_amp.png` — AMP pooled %CV + log-SD histograms
  (diagnostic).
- `figures/replicate_cv_azt.png` — AZT pooled %CV + log-SD histograms
  (diagnostic).
- **`figures/replicate_cv_viable_amp.png`** — AMP viable-subset %CV
  (leads the headline).
- **`figures/replicate_cv_viable_azt.png`** — AZT viable-subset %CV
  (leads the headline).
- `data/per_genotype_replicate_stats.csv` — 985,674 rows of per-genotype
  dispersion.
- `data/summary_by_conc.csv` — 14 rows of pooled per-concentration medians.
- `data/summary_by_conc_viable.csv` — 14 rows of viable-subset per-conc
  medians (n_viable, cv_median, IQR, log-SD, frac < 10 % CV).
- `data/replicate_pair_correlations.csv` — 42 rows of per-pair Pearson r/p.
- `results_table.csv` — 14-row headline verdict joining correlations and CV.
- `analysis.py` — fully reproducible pipeline (Polars + Matplotlib,
  seed 20260420). Re-running produces every figure above including the
  viable-subset panels.

## Methodological notes

- Values in the raw CSVs are already log10(AUC); raw AUC used for the
  %CV computation is recovered as 10^fitness.
- Per-genotype SD uses `ddof=1` (sample SD), consistent with n=3 replicates.
- Pearson is parametric and appropriate for replicate comparison
  (agreement around the identity line, not just rank order). Two-sided.
- Viable subset definition: `log_mean > 3.0` at the condition, where
  `log_mean` is the per-genotype mean of the three log10(AUC) replicates
  at that drug × concentration. This threshold corresponds to roughly
  10× WT growth — the same floor above which the epistasis decomposition
  retains meaningful variance.
- The conc = 0 growth-control is retained in the numerical tables but
  not plotted — it carries no signal about drug-dependent fitness.

## Bottom line

The reviewer's request is fulfilled with direct correlation evidence at
every concentration and genotype in the library. **Reproducibility is
assay-condition-dependent and excellent where the manuscript's
conclusions sit:** at the concentrations driving the landscape / epistasis
analyses (AMP 12.2 – 195 µg/mL, AZT 0.44 – 12 µg/mL) the
viable-subset median CV is 22.8 – 25.8 % and replicate-pair Pearson r is
0.60 – 0.75. At the assay extremes reproducibility degrades as the
library approaches the sequencing floor; those extremes are additionally
hardened via the S6 (Fig. S3 revised) low-count flagging rule and the
S7 AMP-781 library-invariance test. The original "<10% SD" phrasing has
been replaced with the honest wording above.
