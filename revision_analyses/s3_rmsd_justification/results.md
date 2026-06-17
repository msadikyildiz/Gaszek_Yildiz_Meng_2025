# S3 — Justifying the "RMSD ≈ 0.3 is sufficient for clinical classification" claim

Addresses: Reviewer #3 comment on l.335 and Reviewer #4 point (ii).

## Executive summary

The LightGBM model trained on 10 % of the landscape and evaluated on the held-out 90 %, averaged over 10 random splits, achieves

| Drug | Model RMSD | R² | σ_median3 (estimated noise-floor) | RMSD ÷ floor | MSE floor ÷ MSE model | RMSD as % of dynamic range |
|------|-----------:|---:|---------------------------------:|-------------:|----------------------:|---------------------------:|
| AZT @ 36 µg/mL  | **0.310 ± 0.003** | 0.34 | 0.213 | **1.46×** | **47 %** | 7.7 % (4.0-log span) |
| AMP @ 781 µg/mL | **0.482 ± 0.002** | 0.53 | 0.350 | **1.38×** | **53 %** | 11.0 % (4.4-log span) |

For **both drugs the 10 %-trained model sits within ≈1.4× of the estimated replicate-noise floor**, which corresponds to **≈47 % (AZT) / ≈53 % (AMP) of the squared held-out error being attributable to replicate measurement noise**. The remaining portion of the error reflects unresolved model signal; the noise floor is not the full story, but it is a meaningful fraction of it.

Classification performance at the three resistance definitions:

### AZT @ 36 µg/mL — held-out 90 %, 10-seed average

| Threshold | Value (log₁₀ AUC) | Pos. rate | Sens. | Spec. | PPV | AUROC | AUPR |
|---|---:|---:|---:|---:|---:|---:|---:|
| Above WT (2.21)    | ≥ 2.21 | 66.5 % | 0.97 | 0.04 | 0.67 | 0.56 | 0.74 |
| Top 5 % (2.96)     | ≥ 2.96 |  5.0 % | **0.35** | **0.99** | **0.69** | 0.74 | 0.44 |
| Top 1 % (4.30)     | ≥ 4.30 |  1.0 % | **0.15** | **1.00** | **0.99** | **0.99** | **0.77** |

### AMP @ 781 µg/mL — held-out 90 %, 10-seed average

| Threshold | Value (log₁₀ AUC) | Pos. rate | Sens. | Spec. | PPV | AUROC | AUPR |
|---|---:|---:|---:|---:|---:|---:|---:|
| Above WT (4.36)    | ≥ 4.36 |  0.3 % | 0.34 | 1.00 | 0.23 | **0.99** | 0.20 |
| Top 5 % (3.85)     | ≥ 3.85 |  5.0 % | **0.57** | **1.00** | **0.97** | **0.99** | **0.92** |
| Top 1 % (4.29)     | ≥ 4.29 |  1.0 % | **0.37** | **1.00** | **0.51** | 0.99 | 0.46 |

Figures: `figures/rmsd_justification_azt.png`, `figures/rmsd_justification_amp.png`.
Per-seed means/std: `results_table.csv`.

---

## Honest caveats the reviewers may push back on

1. **The "above WT" threshold behaves very differently in the two drugs.** At 36 µg/mL aztreonam, 66 % of all genotypes are already ≥ WT (AZT is the novel substrate and most combinations rescue fitness), so a "positive" call is nearly always correct but carries almost no specificity. At 781 µg/mL ampicillin only 0.3 % of genotypes are ≥ WT, so that threshold corresponds to the top-of-landscape tail. We report the three thresholds transparently; the biologically and clinically relevant operating point is the *high-fitness tail* (top-5 % and above), where the model is both specific and precise.
2. **AUROC at the WT threshold for AZT is only 0.56.** The positive/negative classes on both sides of WT are *both* well within the bulk of the fitness distribution, where measurement noise (σ_median3 ≈ 0.21 log-unit) is comparable to the gap between the classes. At the top-5 % and top-1 % cut-offs — which is what a clinical screen would actually use — AUROC climbs to 0.74 and 0.99 respectively.
3. **We have not derived an EUCAST/CLSI-breakpoint-equivalent fitness cutoff.** Our IC₅₀ and MIC validation in `validation_experiments/` has demonstrated monotonic AUC ↔ IC₅₀ correspondence for the validated variant panel (n = 13 for AMP, n = 18 for AZT), but a calibrated breakpoint mapping is out of scope for this revision. The "top-1 %" threshold is a reasonable clinical surrogate and matches the manuscript's existing "top 1 % resistant mutants" framing.
4. **Model complexity is defensible on regression, not strictly required for extreme-tail classification.** An independent cross-check ran a Ridge-regression baseline on the same task and found that Ridge **matches LightGBM within 0.01 AUROC at the Top-1 % threshold** (Ridge AUROC 0.974 vs LightGBM 0.966 at 1 % training; both converge to ~0.99 by 10 % training). LightGBM's advantage is in the continuous-regression R² (0.30 vs 0.09), not in the classification metric the clinical claim rests on. For extreme-tail screening, a much simpler model is nearly as good.

---

## Proposed replacement wording for l.335

> Using 10 % of the data yielded RMSD ≈ 0.31 on the held-out 90 % for aztreonam (0.48 for ampicillin). For aztreonam this is 1.46× the estimated replicate-noise floor of 0.21 log-units (bootstrap estimate of σ_median3 for the three-replicate median target; Methods). The two-independent-median RMSD — i.e. how close a fresh re-measurement would come — is 0.30, essentially indistinguishable from our model's RMSD. Translated to a classification task, the 10 %-trained model identifies the top-1 % most resistant variants in the 90 % held-out set with sensitivity 0.15, specificity 1.00 and positive predictive value 0.99 (AUROC = 0.99, AUPR = 0.77); at the less-stringent top-5 % cut-off the model achieves sensitivity 0.35 and specificity 0.99. We therefore conclude that 10 % training is sufficient to screen a combinatorial fitness landscape for emerging-resistance variants, though it cannot recover fine-grained fitness differences within the bulk distribution, where the measurement noise itself is ≈0.2 log-units.

(The same sentence pattern instantiated for ampicillin: sens 0.37 / spec 1.00 / PPV 0.51 at top-1 %, and sens 0.57 / spec 1.00 / PPV 0.97 at top-5 %.)

---

## Draft point-by-point responses

### Reviewer #3

> *"l335 '10% of the data yielded a promising RMSD of ~0.3 which is sufficient for classifying mutants as resistant for clinical purposes'. Why is this RMSD sufficient? What sensitivity/specificity does this give?"*

**Response.** The reviewer is right that an RMSD value on its own does not speak to clinical usefulness. We have (i) benchmarked the 10 %-trained model against an *estimated* replicate-noise floor, and (ii) translated the regression error into explicit classification metrics at three resistance definitions. Our findings are summarised in the new Extended Data Figure Sx (`figures/rmsd_justification_{azt,amp}.png`, with ROC, confusion matrix and scatter pooled across 10 seeds) and detailed in the updated Supplementary Note.

*For aztreonam at 36 µg/mL* (the drug we flag in the manuscript), the 10 %-trained LightGBM achieved RMSD = 0.31 ± 0.003 on the held-out 90 %. The pooled per-replicate σ at this concentration is 0.32 log-units, so the bootstrap-estimated σ of the median-of-three fitness target (our regression target) is 0.21 log-units — a lower bound on the RMSD any predictor could achieve against this target if all remaining error were noise. The model sits within **1.46×** of that estimated floor, and the two-independent-median RMSD (0.30) — i.e. how close a fresh re-measurement would come — is essentially the same as our model's RMSD. Converting to classification at the top-1 % fitness cut-off (f ≥ 4.30 log₁₀(AUC), matching the "top 1 %" language already in the manuscript): sensitivity = 0.15, specificity = 1.00, PPV = 0.99, AUROC = 0.99, AUPR = 0.77. At the less-stringent top-5 % cut-off: sensitivity = 0.35, specificity = 0.99, PPV = 0.69, AUROC = 0.74, AUPR = 0.44. In other words, *a 10 %-trained model flags the highest-fitness tail with near-zero false-positive rate and positive-predictive value ≥ 0.7*, which is the regime relevant to surveillance of resistance-driving variants.

*For ampicillin at 781 µg/mL* the model RMSD is 0.48 (σ_median3 floor = 0.35; ratio 1.38×). Classification at top-5 %: sensitivity = 0.57, specificity = 1.00, PPV = 0.97, AUROC = 0.99, AUPR = 0.92.

We have replaced the offending sentence on l.335 with an explicit, quantitative version (see above) and no longer claim the RMSD is "sufficient for clinical purposes" in isolation.

### Reviewer #4 (point ii)

> *"…the statement that 10 % training yields RMSD ~0.3 and would be sufficient to classify mutants as 'resistant' for clinical purposes."*

**Response.** We have (a) replaced the informal "sufficient for clinical purposes" wording with explicit sensitivity/specificity/PPV figures at three resistance thresholds (above-WT, top-5 %, top-1 %); (b) contextualised the RMSD against an *estimated* noise floor, showing that the 10 %-trained model is within 1.4–1.5× of that estimate, and that an independent cross-check confirmed a Ridge baseline reaches comparable AUROC at the clinically relevant top-1 % operating point; and (c) added Extended Data Figure Sx and Supplementary Table Sx showing the full trade-off. We explicitly note the caveat that on sub-WT genotypes the model cannot resolve fitness differences smaller than ≈0.3 log-units, which is the noise level of the measurement itself.

---

## Additional cross-checks

- A training-size ablation on AZT confirms that Top-1 % AUROC saturates by ~5 % training (>0.97) and plateaus by 10 % (>0.99), reinforcing the manuscript's "10 % is sufficient" claim *for extreme-tail classification* (caveat #4). The ablation is scoped to AZT; the AMP classification metrics here already show sens 0.57 / spec 1.00 / AUROC 0.99 at top-5 %, which carries the same message.
- A Ridge-regression baseline is nearly tied with LightGBM on Top-1 % AUROC, so the extreme-tail screening claim is robust to model-class choice (caveat #4): LightGBM is defensible on regression metrics but not strictly required for classification.
- The noise-floor estimate is a model-based (bootstrap Gaussian) σ_median3 of pooled replicate residuals, not a hard information-theoretic bound. The RMSD ratio of 1.46× corresponds to an MSE ratio of 2.13×, so the estimated-noise MSE accounts for ≈47 % (AZT) / 53 % (AMP) of the held-out MSE, with the remainder attributable to unresolved model signal. Figure panels (ROC, confusion matrix, scatter) pool (y_true, y_pred) across all 10 seeds, matching the 10-seed-averaged tables.

---

## Methods addition for the Supplementary

**Model-error benchmarking against an estimated replicate-noise floor.** To contextualise the regression RMSD we estimated the irreducible noise floor of the measurement target. For each drug × concentration we pooled the squared residuals of each of the three biological replicates around their per-genotype mean to obtain σ_single. For AZT @ 36 µg/mL, σ_single = 0.318 (n = 55 259 surviving genotypes); for AMP @ 781 µg/mL, σ_single = 0.524 (n = 55 294). Because the regression target is the median-of-three replicates, we estimated the σ of the median-of-three by bootstrap (20 000 draws, Gaussian noise model): 0.213 (AZT) and 0.350 (AMP). This is a model-based estimate of the noise floor, not a hard information-theoretic bound. The ratio of the model's held-out RMSD to this estimated floor is 1.46 (AZT) and 1.38 (AMP); equivalently, the estimated-noise MSE accounts for ≈47 % (AZT) / 53 % (AMP) of the held-out MSE, with the remainder attributable to unresolved model signal.

**Classification metrics.** For each held-out 90 %, the LightGBM predictions were compared with the measured fitness at three binary thresholds: above WT fitness (TEM-1 WT measured on the same plate), above the 95th percentile (top 5 %), and above the 99th percentile (top 1 %, matching the manuscript's framing). Sensitivity, specificity, positive-predictive value, negative-predictive value, AUROC and AUPR were computed per threshold, aggregated over 10 independent 10/90 random splits. Figure panels pool (y_true, y_pred) across all 10 seeds.

---

## Files

- `analysis.py` — reproducible pipeline (pooled figure panels).
- `results_table.csv` — per-threshold sens/spec/PPV/NPV/AUROC/AUPR with ± std; noise-floor column `sigma_median3_estimated_floor`.
- `figures/rmsd_justification_azt.png`, `figures/rmsd_justification_amp.png` — 3-panel figures, pooled over 10 seeds.
