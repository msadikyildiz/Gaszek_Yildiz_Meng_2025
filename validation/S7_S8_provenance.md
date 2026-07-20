# S7 / S8 provenance — monoculture IC50 vs AUC-fitness

The exact data behind the "n = 13 (AMP) / n = 18 (AZT), Pearson r = 0.88 (AMP) /
0.80 (AZT)" correlation reported in Supplementary Figures S7 and S8.

## The numbers and where they come from

Each headline number is the `mean_fitness × log10_ic50` panel of
`src/fig_c_mic_vs_fitness.py` — monoculture AUC-based IC50 vs mean pooled
AUC-fitness, after the figure's own filter (`drug == D & genotype_13 is not
null`, then `drop_nulls` on the two columns).

| Drug | Batch    | n  | Pearson r | p       |
|------|----------|----|-----------|---------|
| AMP  | 20260407 | 13 | +0.885    | < 0.001 |
| AZT  | 20260307 | 18 | +0.803    | < 0.001 |

These are the two largest single-drug batches. The result is timepoint-invariant:
the same n and r hold at the default, `_12h`, and `_24h` snapshots (the IC50 fit
pools timepoints; only the OD-derived metrics move with timepoint).

## Why these batches

Batch 20260220 ran both drugs but with a smaller shared panel (n = 12 each, AMP
r = 0.59, AZT r = 0.86); we report the later, larger per-drug batches. The
pipeline is batch-scoped (`config.set_batch()` repoints `PROCESSED_DIR`), so
there is no single cross-batch pooled object — a "combined" AMP/AZT scatter is a
two-panel layout, not a single pooled fit.

## Reproduce

From the repository root:

```bash
uv run python - <<'PY'
import polars as pl
from scipy import stats
for drug, b in [("AMP", "20260407"), ("AZT", "20260307")]:
    d = pl.read_parquet(f"validation_experiments/src/processed/{b}/xref_expanded_df.parquet")
    d = (d.filter((pl.col("drug") == drug) & pl.col("genotype_13").is_not_null())
           .select(["mean_fitness", "log10_ic50"]).drop_nulls())
    r, p = stats.pearsonr(d["mean_fitness"], d["log10_ic50"])
    print(drug, b, "n =", len(d), "r =", round(r, 3), "p =", round(p, 4))
PY
```
