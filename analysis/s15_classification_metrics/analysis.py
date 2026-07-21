"""RMSD justification analysis for Gaszek-Yildiz-Meng 2026 (Supplementary Figure S15:
resistance-classification metrics).

Reproduces the 10%-trained LightGBM model (matching notebooks 03/04) for
AZT @ 36 ug/mL and AMP @ 781 ug/mL, then computes:

  (1) sensitivity / specificity / PPV / NPV / AUROC / AUPR on the held-out
      90% at three resistance-threshold definitions:
        WT     - fitness above TEM-1(WT) @ that concentration
        TOP1   - fitness in top 1% (99th percentile, manuscript framing)
        TOP5   - fitness in top 5% (less stringent comparator)
  (2) replicate-pair RMSD (single-rep-vs-single-rep) and single-rep-vs-
      median-of-other-two (noise floor on the target the model is trained
      to predict)
  (3) Fitness dynamic range → express RMSD as fraction

Outputs:
  figures/rmsd_justification_{amp,azt}.png
  results_table.csv
  data/*.parquet (for inspection)

Uses Gaszek_Yildiz_Meng_2025/.venv/bin/python (Polars, scikit-learn, LightGBM).
"""
from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))
os.environ.setdefault("FONTCONFIG_CACHE", os.path.join(tempfile.gettempdir(), "fontconfig"))

import numpy as np
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import lightgbm as lgbm
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    mean_squared_error,
    r2_score,
    roc_curve,
    auc,
    precision_recall_curve,
    average_precision_score,
    confusion_matrix,
)

# ─────────────────────────────────────────────────────────────────────────────
# Configuration

def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
DATA_DIR = REPO / "data" / "processed"
OUT_DIR = Path(__file__).resolve().parent
FIG_DIR = OUT_DIR / "figures"
SUPP = REPO / "figures" / "supplementary"
FIG_DIR.mkdir(exist_ok=True, parents=True)
(OUT_DIR / "data").mkdir(exist_ok=True, parents=True)
SUPP.mkdir(exist_ok=True, parents=True)

# Ambler numbering – matches notebooks 03/04
SIGN_POS = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]

# Wild-type sequence at the 13 Ambler positions (from config.py INTENDED first AA)
WT_SEQ = "LQMERMAGERTRN"

# Features used in the published LightGBM (drops the WT one-hot per position)
FEATURE_COLS = [
    "L21P", "Q39K", "M69L", "M69V", "E104K",
    "R164H", "R164N", "R164S", "M182T", "A237T",
    "G238S", "E240K", "R244C", "R244S", "T265M",
    "R275L", "R275Q", "N276D",
]
COLUMN_MAPPING = {
    "21_P": "L21P", "39_K": "Q39K", "69_L": "M69L", "69_V": "M69V",
    "104_K": "E104K", "164_H": "R164H", "164_N": "R164N", "164_S": "R164S",
    "182_T": "M182T", "237_T": "A237T", "238_S": "G238S", "240_K": "E240K",
    "244_C": "R244C", "244_S": "R244S", "265_M": "T265M", "275_L": "R275L",
    "275_Q": "R275Q", "276_D": "N276D",
}

# Config for each drug – exact concentration used in the original notebooks
DRUG_CONFIG = {
    "AMP": {"drug_name": "AMP", "conc": 781.0, "long_file": "amp_auc_long_df.parquet",
            "display": "Ampicillin (781 µg/mL)", "color": "#6E7FAE"},
    "AZT": {"drug_name": "AZT", "conc": 36.0, "long_file": "azt_auc_long_df.parquet",
            "display": "Aztreonam (36 µg/mL)", "color": "#C15B9C"},
}

RANDOM_STATE_TRAIN = 42  # fixed for reproducibility; notebooks use random seed
TRAIN_SIZE = 0.10
N_REPEATS = 10            # how many seeds to average metrics over (matches notebooks)

# ─────────────────────────────────────────────────────────────────────────────
# Utilities


def dots_to_aa(dot_genotype: str, wt: str = WT_SEQ) -> str:
    """Convert long_df genotype ('.' = WT) → Epistasis_Combined genotype (real AA)."""
    return "".join(wt[i] if c == "." else c for i, c in enumerate(dot_genotype))


def build_feature_matrix(epi_df_pd: pd.DataFrame) -> pd.DataFrame:
    """Match the exact feature matrix produced by the published notebooks."""
    epi_df_pd = epi_df_pd.copy()
    epi_df_pd["genotype_split"] = epi_df_pd["Genotype"].apply(list)
    genotype_df = pd.DataFrame(
        epi_df_pd["genotype_split"].tolist(), index=epi_df_pd.index
    )
    genotype_df.columns = [f"{SIGN_POS[i]}" for i in genotype_df.columns]
    df = pd.concat([epi_df_pd, genotype_df], axis=1)
    df = df.drop(["Genotype", "genotype_split"], axis=1)
    df_enc = pd.get_dummies(df, columns=genotype_df.columns.tolist())

    cols_needed = [f"{pos}_{aa}" for pos, aa in [
        (21, "P"), (39, "K"), (69, "L"), (69, "V"), (104, "K"),
        (164, "H"), (164, "N"), (164, "S"), (182, "T"), (237, "T"),
        (238, "S"), (240, "K"), (244, "C"), (244, "S"), (265, "M"),
        (275, "L"), (275, "Q"), (276, "D"),
    ]]
    for c in cols_needed:
        if c not in df_enc.columns:
            df_enc[c] = 0
    df_enc = df_enc[cols_needed + ["Fitness"]]
    df_enc = df_enc.rename(columns=COLUMN_MAPPING)
    return df_enc


def train_lightgbm(
    X: pd.DataFrame, y: pd.Series, train_size: float, seed: int
) -> tuple[lgbm.LGBMRegressor, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce the notebook's 10%/90% LightGBM training."""
    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y, train_size=train_size, random_state=seed
    )
    model = lgbm.LGBMRegressor(
        objective="regression",
        n_estimators=100,
        learning_rate=0.1,
        random_state=42,
        verbose=-1,
    )
    model.fit(X_tr, y_tr)
    y_pred = model.predict(X_te)
    return model, X_te, y_te.values, y_pred, X_tr.index.values


def classification_metrics(y_true: np.ndarray, y_pred: np.ndarray, threshold: float) -> dict:
    """Compute sens/spec/PPV/NPV/AUROC/AUPR at a single threshold."""
    tl = (y_true >= threshold).astype(int)
    pl = (y_pred >= threshold).astype(int)
    n_pos = int(tl.sum())
    n_neg = int(len(tl) - n_pos)

    if n_pos == 0 or n_neg == 0:
        return {"n_pos": n_pos, "n_neg": n_neg}

    tn, fp, fn, tp = confusion_matrix(tl, pl, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    ppv = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan

    # Continuous-score metrics use the raw prediction as the score
    fpr, tpr, _ = roc_curve(tl, y_pred)
    auroc = auc(fpr, tpr)
    aupr = average_precision_score(tl, y_pred)

    return {
        "threshold": float(threshold),
        "n_pos": n_pos, "n_neg": n_neg,
        "positive_rate": n_pos / len(tl),
        "tp": int(tp), "fp": int(fp), "fn": int(fn), "tn": int(tn),
        "sensitivity": float(sens),
        "specificity": float(spec),
        "precision_ppv": float(ppv),
        "npv": float(npv),
        "auroc": float(auroc),
        "aupr": float(aupr),
        "fpr_curve": fpr, "tpr_curve": tpr,
    }


def noise_floor(long_df_conc: pl.DataFrame) -> dict:
    """Replicate-level noise metrics.

    The target the LightGBM predicts is the median of 3 biological replicates.
    Several "noise floor" numbers are informative:

      * `pair_rmsd_mean` — RMSD between any two single replicates.  Directly
        measurable; approximately √2·σ_single.  Useful comparator for
        "replicate-level reproducibility".
      * `sigma_single_est` — pooled residual-from-mean σ per single replicate.
      * `sigma_median3_est` — empirical σ of the median of 3 reps.  Equals the
        *estimated noise-floor* RMSD: the minimum RMSD that any predictor, no matter
        how powerful, could achieve against the median-of-3 target, because
        this much noise separates the measured target from the unobserved true
        fitness.
      * `rmsd_two_medians` — RMSD between two independent median-of-3
        measurements of the same genotype ≈ √2·σ_median3.  Practical upper
        bound on the "replicate-level" benchmark — a predictor that hit this
        RMSD would be indistinguishable from re-running the whole experiment.
    """
    r1 = long_df_conc["replicate1"].to_numpy()
    r2 = long_df_conc["replicate2"].to_numpy()
    r3 = long_df_conc["replicate3"].to_numpy()
    stacked = np.stack([r1, r2, r3])
    keep = ~np.any(np.isnan(stacked), axis=0)
    stacked = stacked[:, keep]

    pairs = [(0, 1, "r1_r2"), (0, 2, "r1_r3"), (1, 2, "r2_r3")]
    pair_rmsds = {
        name: float(np.sqrt(np.mean((stacked[i] - stacked[j]) ** 2)))
        for i, j, name in pairs
    }
    pair_mean = float(np.mean(list(pair_rmsds.values())))

    # pooled σ_single from residual to the 3-rep mean (unbiased for Gaussian)
    mu = stacked.mean(axis=0)
    resid = stacked - mu
    sigma_single = float(np.sqrt(np.mean(resid ** 2) * 3.0 / 2.0))

    # bootstrap σ_median3 and RMSD between two independent median-of-3 draws
    rng = np.random.default_rng(0)
    n_bs = 20_000
    truth = rng.choice(mu, size=n_bs, replace=True)
    draws_a = truth[:, None] + rng.normal(0, sigma_single, size=(n_bs, 3))
    draws_b = truth[:, None] + rng.normal(0, sigma_single, size=(n_bs, 3))
    med_a = np.median(draws_a, axis=1)
    med_b = np.median(draws_b, axis=1)
    sigma_median3 = float(np.std(med_a - truth, ddof=1))
    rmsd_two_medians = float(np.sqrt(np.mean((med_a - med_b) ** 2)))

    return {
        "pair_rmsds": pair_rmsds,
        "pair_rmsd_mean": pair_mean,
        "sigma_single_est": sigma_single,
        "sigma_median3_est": sigma_median3,  # estimated noise-floor model RMSD floor
        "rmsd_two_medians": rmsd_two_medians,  # replicate-level benchmark
    }


# ─────────────────────────────────────────────────────────────────────────────
# Main per-drug pipeline


def run_drug(drug_key: str) -> dict:
    cfg = DRUG_CONFIG[drug_key]
    print(f"\n=== {drug_key} @ {cfg['conc']} ug/mL ===")

    # Load data -------------------------------------------------------------
    epi = pl.read_parquet(DATA_DIR / "Epistasis_Combined.parquet")
    epi_sub = epi.filter(
        (pl.col("Drug") == cfg["drug_name"]) & (pl.col("Concentration") == cfg["conc"])
    ).to_pandas()
    print(f"  epi subset shape: {epi_sub.shape}")

    # Drop any nan-fitness genotypes (shouldn't be any, but be safe)
    epi_sub = epi_sub.dropna(subset=["Fitness"]).reset_index(drop=True)

    # WT fitness (at this concentration) — genotype == WT string in real AA
    wt_row = epi_sub[epi_sub["Genotype"] == WT_SEQ]
    if len(wt_row) == 0:
        print(f"  WARNING: WT genotype not found in epi table")
        wt_fitness = np.nan
    else:
        wt_fitness = float(wt_row["Fitness"].values[0])
    print(f"  WT fitness: {wt_fitness:.3f}")

    df_enc = build_feature_matrix(epi_sub)
    X = df_enc.drop("Fitness", axis=1)
    y = df_enc["Fitness"]

    # Multi-seed training — aggregate predictions for all 90% test partitions
    print(f"  Training LightGBM with 10%/90% split × {N_REPEATS} seeds ...")
    rmsd_per_seed = []
    r2_per_seed = []
    # Pool predictions across all seeds for figure panels (scatter, confusion,
    # ROC) so plots reflect the same evidence as the numeric summary rather
    # than a single arbitrary split. Per-seed rows are kept for the tables.
    pooled_y_te_list = []
    pooled_y_pred_list = []
    headline = None
    for seed_idx in range(N_REPEATS):
        seed = 1000 + seed_idx * 7
        model, X_te, y_te, y_pred, train_idx = train_lightgbm(X, y, TRAIN_SIZE, seed)
        rmsd = math.sqrt(mean_squared_error(y_te, y_pred))
        r2 = r2_score(y_te, y_pred)
        rmsd_per_seed.append(rmsd)
        r2_per_seed.append(r2)
        pooled_y_te_list.append(np.asarray(y_te))
        pooled_y_pred_list.append(np.asarray(y_pred))
        if seed_idx == N_REPEATS - 1:
            headline = dict(y_te=y_te, y_pred=y_pred, X_te=X_te, model=model)
    pooled_y_te = np.concatenate(pooled_y_te_list)
    pooled_y_pred = np.concatenate(pooled_y_pred_list)
    rmsd_mean = float(np.mean(rmsd_per_seed))
    rmsd_std = float(np.std(rmsd_per_seed))
    r2_mean = float(np.mean(r2_per_seed))
    r2_std = float(np.std(r2_per_seed))
    print(f"  RMSD = {rmsd_mean:.4f} ± {rmsd_std:.4f}  (range "
          f"{min(rmsd_per_seed):.4f}-{max(rmsd_per_seed):.4f})")
    print(f"  R²   = {r2_mean:.4f} ± {r2_std:.4f}")

    # Build thresholds -----------------------------------------------------
    y_all = y.values
    thresholds = {
        "WT": {"value": wt_fitness, "label": f"Above WT (fitness ≥ {wt_fitness:.2f})"},
        "TOP5": {
            "value": float(np.percentile(y_all, 95)),
            "label": "Top 5% (95th pctile)",
        },
        "TOP1": {
            "value": float(np.percentile(y_all, 99)),
            "label": "Top 1% (99th pctile)",
        },
    }
    print("  Thresholds:")
    for k, v in thresholds.items():
        print(f"    {k}: {v['value']:.3f}  ({v['label']})")

    y_te = headline["y_te"]
    y_pred = headline["y_pred"]

    # Classification metrics per threshold (aggregate over all seeds too) --
    all_metrics = {}
    for key, thr in thresholds.items():
        rows = []
        for seed_idx in range(N_REPEATS):
            seed = 1000 + seed_idx * 7
            _, _, y_te_k, y_pred_k, _ = train_lightgbm(X, y, TRAIN_SIZE, seed)
            m = classification_metrics(y_te_k, y_pred_k, thr["value"])
            if "sensitivity" not in m:
                continue
            rows.append(m)
        # aggregate
        agg = {}
        for field in ["sensitivity", "specificity", "precision_ppv",
                      "npv", "auroc", "aupr"]:
            vals = [r[field] for r in rows]
            agg[f"{field}_mean"] = float(np.mean(vals))
            agg[f"{field}_std"] = float(np.std(vals))
        agg["threshold_value"] = thr["value"]
        agg["threshold_label"] = thr["label"]
        agg["headline"] = rows[-1]  # for plotting
        agg["n_positive_mean"] = float(np.mean([r["n_pos"] for r in rows]))
        agg["positive_rate_mean"] = float(np.mean([r["positive_rate"] for r in rows]))
        all_metrics[key] = agg
        print(f"    [{key}]  sens={agg['sensitivity_mean']:.3f}±{agg['sensitivity_std']:.3f} "
              f"spec={agg['specificity_mean']:.3f}±{agg['specificity_std']:.3f} "
              f"PPV={agg['precision_ppv_mean']:.3f}  "
              f"AUROC={agg['auroc_mean']:.3f}  AUPR={agg['aupr_mean']:.3f}")

    # Replicate noise floor -----------------------------------------------
    long_df = pl.read_parquet(DATA_DIR / cfg["long_file"])
    long_sub = long_df.filter(pl.col("concentration") == cfg["conc"])
    nf = noise_floor(long_sub)
    print(f"  σ_single       = {nf['sigma_single_est']:.4f}")
    print(f"  σ_median3      = {nf['sigma_median3_est']:.4f}   (estimated noise-floor model RMSD floor)")
    print(f"  pair RMSD      = {nf['pair_rmsd_mean']:.4f}      (single-rep vs single-rep)")
    print(f"  two-median RMSD= {nf['rmsd_two_medians']:.4f}     (replicate-level benchmark)")

    # Dynamic range --------------------------------------------------------
    y_range = float(y_all.max() - y_all.min())
    y_iqr_range = float(np.percentile(y_all, 99) - np.percentile(y_all, 1))
    print(f"  Dynamic range: max-min = {y_range:.3f}, p99-p1 = {y_iqr_range:.3f}")
    print(f"  RMSD as fraction of range: {rmsd_mean/y_range*100:.2f}% (full), "
          f"{rmsd_mean/y_iqr_range*100:.2f}% (p99-p1)")

    return {
        "drug": drug_key,
        "cfg": cfg,
        "n_genotypes": len(y_all),
        "wt_fitness": wt_fitness,
        "rmsd_mean": rmsd_mean, "rmsd_std": rmsd_std,
        "r2_mean": r2_mean, "r2_std": r2_std,
        "thresholds": thresholds,
        "metrics": all_metrics,
        "noise_floor": nf,
        "dynamic_range": y_range,
        "iqr_range": y_iqr_range,
        "headline": headline,
        "pooled_y_te": pooled_y_te,
        "pooled_y_pred": pooled_y_pred,
        "y_all": y_all,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Plotting


def plot_drug(res: dict) -> Path:
    drug = res["drug"]
    cfg = res["cfg"]
    accent = cfg["color"]
    rmsd = res["rmsd_mean"]
    sigma_med = res["noise_floor"]["sigma_median3_est"]
    two_med_rmsd = res["noise_floor"]["rmsd_two_medians"]
    pair_rmsd = res["noise_floor"]["pair_rmsd_mean"]
    y_all = res["y_all"]

    fig = plt.figure(figsize=(13, 4.2), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.42)

    # Pool predictions across all seeds so ROC / confusion / scatter reflect
    # the same evidence as the multi-seed numeric summary.
    pooled_y_te = res["pooled_y_te"]
    pooled_y_pred = res["pooled_y_pred"]

    from sklearn.metrics import roc_curve as _roc
    # ── Panel A: ROC curves across threshold definitions (pooled, 10 seeds)
    ax_a = fig.add_subplot(gs[0, 0])
    curve_colors = {"WT": "#222222", "TOP5": accent, "TOP1": "#E07B00"}
    for key in ["WT", "TOP5", "TOP1"]:
        thr_val_a = res["thresholds"][key]["value"]
        y_true_bin = (pooled_y_te >= thr_val_a).astype(int)
        if y_true_bin.sum() == 0 or y_true_bin.sum() == len(y_true_bin):
            continue
        fpr_c, tpr_c, _ = _roc(y_true_bin, pooled_y_pred)
        label = (f"{res['thresholds'][key]['label']}\n"
                 f"AUROC={res['metrics'][key]['auroc_mean']:.3f}  "
                 f"AUPR={res['metrics'][key]['aupr_mean']:.3f}")
        ax_a.plot(fpr_c, tpr_c, color=curve_colors[key], lw=1.8, label=label)
    ax_a.plot([0, 1], [0, 1], color="#AAAAAA", lw=0.8, ls="--")
    ax_a.set_xlabel("False positive rate", fontsize=10)
    ax_a.set_ylabel("True positive rate", fontsize=10)
    ax_a.set_title(f"a  ROC — {drug} @ {cfg['conc']:g} µg/mL  (pooled over 10 seeds)",
                   fontsize=9, loc="left", fontweight="bold")
    ax_a.legend(fontsize=7, loc="lower right", frameon=False)
    ax_a.set_aspect("equal")
    for s in ax_a.spines.values():
        s.set_linewidth(0.6)

    # ── Panel B: confusion matrix at TOP5 (summed counts across all seeds)
    ax_b = fig.add_subplot(gs[0, 1])
    key_cm = "TOP5"
    thr_val = res["thresholds"][key_cm]["value"]
    yt_bin = (pooled_y_te >= thr_val).astype(int)
    yp_bin = (pooled_y_pred >= thr_val).astype(int)
    tp = int(((yt_bin == 1) & (yp_bin == 1)).sum())
    tn = int(((yt_bin == 0) & (yp_bin == 0)).sum())
    fp = int(((yt_bin == 0) & (yp_bin == 1)).sum())
    fn = int(((yt_bin == 1) & (yp_bin == 0)).sum())
    mat = np.array([[tn, fp], [fn, tp]])
    im = ax_b.imshow(mat, cmap="RdPu" if drug == "AZT" else "Blues",
                     norm=LogNorm(vmin=max(1, mat.min())))
    for i in range(2):
        for j in range(2):
            val = mat[i, j]
            ax_b.text(j, i, f"{val:,}", ha="center", va="center",
                      fontsize=11,
                      color="white" if val > mat.max() / 4 else "black",
                      fontweight="bold")
    ax_b.set_xticks([0, 1])
    ax_b.set_xticklabels([f"Pred. < p95", f"Pred. ≥ p95"])
    ax_b.set_yticks([0, 1])
    ax_b.set_yticklabels([f"True < p95", f"True ≥ p95"])
    sens = res["metrics"][key_cm]["sensitivity_mean"]
    spec = res["metrics"][key_cm]["specificity_mean"]
    ppv = res["metrics"][key_cm]["precision_ppv_mean"]
    ax_b.set_title(
        f"b  Confusion @ Top-5% (fitness ≥ {thr_val:.2f})\n"
        f"   sens={sens:.2f}  spec={spec:.2f}  PPV={ppv:.2f}",
        fontsize=9, loc="left", fontweight="bold")
    for s in ax_b.spines.values():
        s.set_linewidth(0.6)

    # ── Panel C: predicted-vs-measured scatter + RMSD + noise floor band
    ax_c = fig.add_subplot(gs[0, 2])
    y_te = pooled_y_te; y_pred = pooled_y_pred
    fmin = min(float(y_te.min()), float(y_pred.min())) - 0.2
    fmax = max(float(y_te.max()), float(y_pred.max())) + 0.2
    bins = np.linspace(fmin, fmax, 50)
    h, xe, ye, img = ax_c.hist2d(
        y_te, y_pred, bins=[bins, bins],
        cmap="RdPu" if drug == "AZT" else "Blues", norm=LogNorm())
    ax_c.plot([fmin, fmax], [fmin, fmax], color="#333333", lw=1.0, ls="--")
    # ± RMSD band
    xs = np.linspace(fmin, fmax, 100)
    ax_c.fill_between(xs, xs - rmsd, xs + rmsd, color=accent, alpha=0.12,
                      label=f"± Model RMSD = {rmsd:.3f}", zorder=1)
    # ± estimated noise-floor (σ_median3) band — minimum achievable
    ax_c.fill_between(xs, xs - sigma_med, xs + sigma_med, color="#555555",
                      alpha=0.22,
                      label=f"± σ_median3 (estimated noise-floor floor) = {sigma_med:.3f}",
                      zorder=1)
    ax_c.set_xlim(fmin, fmax); ax_c.set_ylim(fmin, fmax)
    ax_c.set_xlabel("Measured fitness  log₁₀(AUC)", fontsize=10)
    ax_c.set_ylabel("Predicted fitness  log₁₀(AUC)", fontsize=10)
    ax_c.set_title(
        f"c  Predicted vs measured  (n={len(y_te):,} pooled over 10 seeds, R²={res['r2_mean']:.2f})",
        fontsize=10, loc="left", fontweight="bold")
    # annotate
    ax_c.text(0.03, 0.97,
              f"Model RMSD     = {rmsd:.3f}\n"
              f"σ_median3 floor = {sigma_med:.3f}   "
              f"({rmsd/sigma_med:.2f}× floor)\n"
              f"Two-median RMSD = {two_med_rmsd:.3f}\n"
              f"Pair (single-rep) RMSD = {pair_rmsd:.3f}\n"
              f"Dynamic range = {res['dynamic_range']:.2f}   "
              f"RMSD / range = {rmsd/res['dynamic_range']*100:.1f}%",
              transform=ax_c.transAxes, fontsize=7.5, va="top", ha="left",
              family="monospace",
              bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
    ax_c.legend(fontsize=7, loc="lower right", frameon=False)
    ax_c.set_aspect("equal")
    for s in ax_c.spines.values():
        s.set_linewidth(0.6)

    fig.tight_layout()
    out = SUPP / f"figure_s15_{drug.lower()}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  saved {out}")
    plt.close(fig)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Main


def main():
    results = {}
    for drug in ["AZT", "AMP"]:
        results[drug] = run_drug(drug)
        plot_drug(results[drug])

    # Results table
    rows = []
    for drug, res in results.items():
        cfg = res["cfg"]
        for key in ["WT", "TOP5", "TOP1"]:
            m = res["metrics"][key]
            rows.append({
                "drug": drug,
                "concentration_ug_mL": cfg["conc"],
                "threshold_name": key,
                "threshold_value": m["threshold_value"],
                "threshold_label": m["threshold_label"],
                "n_positive_mean": m["n_positive_mean"],
                "positive_rate_mean": m["positive_rate_mean"],
                "sensitivity_mean": m["sensitivity_mean"],
                "sensitivity_std": m["sensitivity_std"],
                "specificity_mean": m["specificity_mean"],
                "specificity_std": m["specificity_std"],
                "precision_ppv_mean": m["precision_ppv_mean"],
                "precision_ppv_std": m["precision_ppv_std"],
                "npv_mean": m["npv_mean"],
                "npv_std": m["npv_std"],
                "auroc_mean": m["auroc_mean"],
                "auroc_std": m["auroc_std"],
                "aupr_mean": m["aupr_mean"],
                "aupr_std": m["aupr_std"],
                "model_rmsd_mean": res["rmsd_mean"],
                "model_rmsd_std": res["rmsd_std"],
                "r2_mean": res["r2_mean"],
                "wt_fitness": res["wt_fitness"],
                "replicate_pair_rmsd": res["noise_floor"]["pair_rmsd_mean"],
                "sigma_single_est": res["noise_floor"]["sigma_single_est"],
                "sigma_median3_estimated_floor": res["noise_floor"]["sigma_median3_est"],
                "rmsd_two_independent_medians": res["noise_floor"]["rmsd_two_medians"],
                "dynamic_range": res["dynamic_range"],
                "iqr_range_p99_p1": res["iqr_range"],
                "rmsd_fraction_of_range_pct": res["rmsd_mean"]
                    / res["dynamic_range"] * 100,
            })
    out_table = pd.DataFrame(rows)
    out_path = OUT_DIR / "results_table.csv"
    out_table.to_csv(out_path, index=False)
    print(f"\nWrote results table → {out_path}")

    # Drop numeric summary for quick inspection
    summary = out_table[[
        "drug", "threshold_name", "threshold_value",
        "sensitivity_mean", "specificity_mean",
        "precision_ppv_mean", "auroc_mean", "aupr_mean",
        "model_rmsd_mean", "sigma_median3_estimated_floor",
        "rmsd_two_independent_medians", "dynamic_range",
    ]].round(3)
    print("\nSummary:")
    print(summary.to_string(index=False))
    return results


if __name__ == "__main__":
    main()
