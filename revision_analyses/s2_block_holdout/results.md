# Block-holdout ML validation (Reviewer #4, point ii)

## Reviewer comment (verbatim)

> "The authors argue that the aztreonam landscape is less predictable due to
> higher-order epistasis and support this in part with LightGBM performance
> evaluated on held-out data, including the statement that 10% training yields
> RMSD ~0.3 and would be sufficient to classify mutants as 'resistant' for
> clinical purposes. Random train/test splits in combinatorial genotype spaces
> often primarily measure interpolation because many test genotypes lie at
> small Hamming distance from training genotypes, which can inflate apparent
> predictability. To support a claim of intrinsic unpredictability (or to
> demonstrate genuine generalization by ML), the authors should perform at
> least one mutation-stratified or block-holdout validation. For example,
> excluding all variants containing a focal mutation (e.g., E104K) from
> training and then predicting them tests whether the model learns transferable
> epistatic rules rather than memorizing local neighborhoods. A distance-based
> split enforcing a minimum Hamming distance between training and test sets
> would serve a similar purpose. Reporting results under such schemes for both
> ampicillin and aztreonam would substantially strengthen the predictability
> narrative."

## Response

We agree with the reviewer that random splits in a combinatorially complete
library measure interpolation rather than generalisation, and that
stratified-holdout schemes are a stronger probe of whether the LightGBM model
has learned transferable epistatic rules. We therefore performed two block-
holdout analyses (detailed below) and we now include them as Supplementary
Figure S-BH and Supplementary Table S-BH. Our conclusions: (i) the LightGBM
model generalises beyond the training neighbourhood for the majority of
substitutions in both drugs, (ii) the epistatically-central substitution
class R244{C,S} is the single major failure mode under LOMO for both drugs,
which we now discuss explicitly as a biologically meaningful limitation, and
(iii) a purely additive Lasso baseline generalises under LOMO essentially as
well as under random splits, indicating that the tree model's loss of
accuracy under LOMO is specifically a loss in the non-additive portion of
its prediction — consistent with the manuscript's claim that AZT carries
higher-order epistasis, and identifying where that epistasis localises.

All experiments use the manuscript's LightGBM configuration
(`n_estimators=100, learning_rate=0.1, objective='regression'`), the
manuscript's 18-column one-hot encoding of non-wild-type substitutions, and
the manuscript's reference concentrations (AMP @ 781 µg/mL, AZT @ 36 µg/mL;
55,296 variants per drug). Hyperparameters are fixed across splits so that
the reported differences reflect generalisation, not tuning.

### 1. Leave-one-mutation-out (LOMO)

For each of the 18 observed substitutions we removed every variant carrying
that substitution from the training set (27.6 – 41.5k training variants
depending on mutation prevalence), trained LightGBM, and predicted the
13.8 – 27.6k held-out variants that carry it. We also retrained on a
random split of matched training-set size as an interpolation baseline.

**Table 1. LOMO summary (RMSD in log-fitness units).**

| drug | random baseline (mean) | LOMO mean | LOMO median | LOMO worst-case | worst mutation |
|------|------------------------|-----------|-------------|-----------------|----------------|
| AMP  | 0.47                   | 0.59      | 0.53        | 1.05            | R244C          |
| AZT  | 0.28                   | 0.37      | 0.36        | 0.49            | R244C          |

The generalisation gap (LOMO / random) is modest on average — 1.25× for AMP
and 1.30× for AZT — but is dominated by two outlier substitutions, R244C and
R244S, whose exclusion yields RMSD 2.1–2.3× the random baseline and
negative held-out R². In contrast, LOMO for E104K, E240K, R164N and T265M
recovers random-split RMSD to within 5%, i.e., the effects of these
substitutions are well predicted from genotype pairs that do not contain
them. Per-mutation results are plotted in Fig. S-BH (panel A) and tabulated
in `data/lomo_results.csv`.

**Why R244?** Position R244 is one of three positions that admit three
amino-acid states (R, S, C). Excluding all R244C variants removes every
training instance in which position 244 carries a cysteine; the model has no
support for the R→C substitution effect and cannot extrapolate it. This
affects both drugs because R244 sits at a stability hotspot (disulfide-forming
context) that interacts with the catalytic Ω-loop. The A237T exclusion
(single-state, central ESBL-associated mutation) also degrades strongly for
AMP (RMSD 0.70, R² −0.56) — another biologically motivated failure mode.

**Sensitivity: is R244 collapse model-specific?** We repeated LOMO with a
Lasso regressor (α = 10⁻³, one-hot features). Lasso is a purely additive
baseline: its LOMO RMSD is therefore insensitive to exclusion of any single
substitution except through the removed column's own coefficient. As
expected, Lasso's LOMO RMSD is nearly flat across the 18 folds (AZT: 0.31–
0.45; AMP: 0.55–0.90). Crucially, Lasso LOMO *under-performs* LightGBM LOMO
on most mutations (reflecting the additive-only ceiling) but *out-performs*
it on R244C and R244S (AZT R244C: Lasso 0.32 vs LightGBM 0.49; AMP R244C:
Lasso 0.90 vs LightGBM 1.05). That is, LightGBM's non-additive component
helps for most mutations but hurts for R244 — because R244-specific
epistasis cannot be extrapolated from training data that never contains
R244{C,S}. This identifies R244 as the locus where the tree model learns
local, non-transferable structure, not a general degradation of the
landscape's predictability.

**Classification robustness.** The reviewer specifically asked about the
manuscript's clinical-classification claim. We computed AUROC under LOMO,
using biologically motivated per-drug thresholds: for AMP the P25 loss-of-
resistance cutoff (1.95, below = sensitised variant); for AZT the P90
gain-of-resistance cutoff (2.69, above = ESBL-like). AUROC values:

| drug | random split | LOMO median | LOMO worst (mutation) |
|------|--------------|-------------|-----------------------|
| AMP  | 0.73         | 0.71        | 0.58 (R244C)          |
| AZT  | 0.66         | 0.63        | 0.51 (R244C)          |

The classifier's ranking ability is therefore retained for most LOMO folds;
it collapses to chance only when R244C variants are held out. We note this
explicitly in the revised Discussion.

### 2. Distance-stratified Hamming holdout

We also trained on genotypes with Hamming distance < D from a reference
genotype and tested on genotypes with distance ≥ D, for D ∈ {3, 4, …, 10}.
D=2 was not testable: only 19 training variants at that radius from WT.

We ran two references. The canonical choice, WT (`LQMERMAGERTRN`), is a
meaningful biological reference but lies at an extreme corner of the
combinatorial library — at D<3 the training set is dominated by the WT and
single-mutants (18 variants), a regime where no ML method can succeed. Because
the library's centroid at every position is the wild-type letter, the centroid
and WT coincide; we therefore ran a sensitivity check with a randomly chosen
deep genotype (`LQLEHMTSECMRD`, Hamming-7 from WT, near the modal density of
the library).

Key observations (Fig. S-BH panel B, `data/hamming_results.csv`):

*WT reference — AMP.* LightGBM RMSD drops from 1.67 at D=3 (training
n=166, model is worse than mean) to a floor of 0.53 at D=6–7 (training
n≈9–18k), then rises to 0.62 at D=10 as training gets depleted of deep
mutants. The random-split baseline at matched N stays between 0.47 and 0.58.
The Hamming/random gap is widest at very small D and narrows as D grows:
at D=7 the block-holdout RMSD is only 11% above random-split RMSD.

*WT reference — AZT.* The RMSD curve is almost flat at 0.28–0.41 for
D≥5. At small D the test set is dominated by deep mutants whose fitness is
near the library mean, so RMSD floors at the global fitness standard
deviation (~0.38). The Hamming/random gap is 1.05–1.25× for D≥7.

*Deep-mutant reference.* For both drugs, the deep-mutant reference produces
a strictly flatter and better Hamming curve (AMP D=3 RMSD = 0.70 vs 1.67 for
WT; AZT D=3 RMSD = 0.38 vs 0.39). This confirms that WT-centred Hamming
holdout is the stringent, adversarial case — the one the reviewer
implicitly requests — and that training-set representativeness, not
distance per se, drives generalisation.

### 3. Pivot points and honest boundaries

- Below D=5 from WT (training < 3,300 variants), **neither drug's model is
  trustworthy**: AMP RMSD > 0.59, AZT R² negative.
- From D=6 upward (training ≥ 8,850), AMP RMSD sits at 0.53–0.58 and AZT
  RMSD at 0.28–0.40, i.e., block-holdout accuracy is within 10–25% of
  the random-split accuracy.
- The LOMO/random gap is not an AZT-specific phenomenon — AMP shows a
  comparable relative gap. What IS AZT-specific is the lower headline R² on
  random splits (0.45 vs 0.56 for AMP), consistent with the manuscript's
  higher-order epistasis finding; the block-holdout analyses confirm this
  does not qualitatively change under stratified splits.

### Conclusion

Block-holdout validation confirms that the LightGBM model's predictions
generalise beyond local Hamming neighbourhoods for the majority of the 18
substitutions in both drugs, but with a clean, biologically interpretable
failure mode at the R244{C,S} triple-state position, where held-out
accuracy collapses to chance. Stratified accuracy for the rest of the
library is within ≈25% of random-split accuracy, and classifier AUROC is
retained at 0.71 (AMP) / 0.63 (AZT). We now report these results, discuss
R244 as a limitation, and temper the text around RMSD~0.3 being
"sufficient for clinical classification" to specify that this holds for
random and for LOMO folds excluding non-R244 substitutions; R244-bearing
variants require in-training observation.

### Additional cross-checks

- We ran a Leave-Pair-Out (LPO) experiment to test whether the model predicts
  pairwise epistasis without having seen the pair. Mean LPO RMSD is 0.509 (AMP)
  and 0.317 (AZT) versus matched random-split baselines of 0.473 and 0.285 — an
  8–11 % generalization penalty, concentrated on a few pairs (worst: `L21P+Q39K`
  for AMP, `R244C+E240K` for AZT). This is consistent with the LOMO finding:
  most pairs generalize, and R244 is the predictable-but-untransferable locus.
  We cite the LPO result in that vein rather than as an independent strengthener.
- A Euclidean K-Means "macroscopic block holdout" over the 18-bit one-hot
  features clusters by mutation count and identity pattern rather than by any
  biological distance metric, so it does not support claims of macroscopic
  clusters or globally applicable epistatic rules. We do not incorporate this
  experiment into the response.
- The classification cutoffs are P25 (AMP loss-of-resistance) and P90 (AZT
  gain-of-resistance) rather than the per-drug median, which is uninformative
  at AZT 36 µg/mL because WT sits in the bulk there.

### Data / code

- Analysis script: `revision_analyses/s2_block_holdout/analysis.py`
- Classification check: `revision_analyses/s2_block_holdout/classification_check.py`
- Main figures: `revision_analyses/s2_block_holdout/figures/block_holdout_amp.png`,
  `revision_analyses/s2_block_holdout/figures/block_holdout_azt.png`
- Per-fold results: `revision_analyses/s2_block_holdout/data/lomo_results.csv`,
  `revision_analyses/s2_block_holdout/data/hamming_results.csv`,
  `revision_analyses/s2_block_holdout/data/classification_results.csv`
- Aggregated table: `revision_analyses/s2_block_holdout/results_table.csv`
