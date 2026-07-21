"""
S1 — Four-model comparison: linear, Lasso, single decision tree, and LightGBM.

Compares four regression model classes on the TEM-1 fitness landscape:
  1. Unregularized linear regression (bare additive features).
  2. L1-regularized linear regression (LassoCV, bare additive features).
  3. Single decision tree (DecisionTreeRegressor, default hyperparameters).
  4. LightGBM (matching manuscript hyperparameters).

To make the comparison fair (LightGBM captures interactions, the bare
additive model cannot), we *additionally* fit the linear and Lasso variants
on one-hot + pairwise-interaction features (171 features).

Learning curve: RMSD and R^2 vs training fraction. Test set is held out once
(5 % of all 55,296 genotypes at a given drug x concentration) and kept FIXED
across every training fraction and seed. Training fractions are expressed as
fractions of the remaining 95 % training pool.

Two drug x concentration conditions (featured in manuscript main figures):
  - AMP at 781 ug/mL
  - AZT at  36 ug/mL

Outputs:
  figures/model_comparison_amp.png
  figures/model_comparison_azt.png
  results_table.csv                   (all model x drug x fraction x seed rows)
  data/learning_curves_summary.csv    (aggregated mean +/- std)
"""

from __future__ import annotations

import hashlib
import math
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
import lightgbm as lgbm
from sklearn.linear_model import LinearRegression, LassoCV
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error, r2_score

# --- paths -------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / "data" / "processed" / "Epistasis_Combined.parquet"
FIGDIR = HERE / "figures"
DATADIR = HERE / "data"
SUPP = REPO / "figures" / "supplementary"
FIGDIR.mkdir(exist_ok=True, parents=True)
DATADIR.mkdir(exist_ok=True, parents=True)
SUPP.mkdir(exist_ok=True, parents=True)

# --- constants ---------------------------------------------------------------
# The fitness parquet encodes each of the 13 positions with its wild-type letter
# when unmutated and with the substitution letter otherwise.
AMBLER_POS = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
WT_LETTER = ["L", "Q", "M", "E", "R", "M", "A", "G", "E", "R", "T", "R", "N"]
ALT_LETTERS = [
    ["P"], ["K"], ["L", "V"], ["K"], ["H", "N", "S"],
    ["T"], ["T"], ["S"], ["K"], ["C", "S"],
    ["M"], ["L", "Q"], ["D"],
]  # 18 substitutions in total

CONDITIONS = [
    ("AMP", 781.0),
    ("AZT", 36.0),
]

# Training fractions as fractions of the 95% training pool.
# 0.99 ~= 99% of pool = ~94% of full data; enforces max ~5% held out.
TRAIN_FRACTIONS = [0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
TEST_FRAC = 0.05
N_SEEDS = 3
TEST_SEED = 20260419            # fixed across all runs — keeps test set identical
SUBSAMPLE_SEEDS = [1001, 1002, 1003]

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


# --- feature construction ----------------------------------------------------
def mutation_labels() -> list[str]:
    labs: list[str] = []
    for pos, wt, alts in zip(AMBLER_POS, WT_LETTER, ALT_LETTERS):
        for a in alts:
            labs.append(f"{wt}{pos}{a}")
    return labs


def build_feature_matrix(df: pl.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """One-hot encode mutations (WT = reference => 18 binary columns)."""
    gts = df["Genotype"].to_numpy()
    y = df["Fitness"].to_numpy().astype(np.float64)
    feats: list[np.ndarray] = []
    labs: list[str] = []
    for i, (wt, alts) in enumerate(zip(WT_LETTER, ALT_LETTERS)):
        col = np.array([g[i] for g in gts])
        for a in alts:
            feats.append((col == a).astype(np.float32))
            labs.append(f"{wt}{AMBLER_POS[i]}{a}")
    X = np.stack(feats, axis=1)
    return X, y, labs


def add_pairwise_interactions(X: np.ndarray, labs: list[str]) -> tuple[np.ndarray, list[str]]:
    """Append pairwise AND-interactions (x_i * x_j) for all i<j."""
    n, d = X.shape
    idx = list(combinations(range(d), 2))
    pw = np.empty((n, len(idx)), dtype=np.float32)
    pw_labs: list[str] = []
    for k, (i, j) in enumerate(idx):
        pw[:, k] = X[:, i] * X[:, j]
        pw_labs.append(f"{labs[i]}*{labs[j]}")
    return np.concatenate([X, pw], axis=1), labs + pw_labs


# --- model definitions -------------------------------------------------------
def fit_linear(Xtr, ytr, Xte):
    m = LinearRegression()
    m.fit(Xtr, ytr)
    return m.predict(Xte)


def fit_lasso(Xtr, ytr, Xte):
    # LassoCV picks alpha by 5-fold CV on the training set.
    # n_jobs=-1 parallelises across CV folds.
    # Cast to float64 — float32 binary features trigger spurious Gram-matrix
    # precompute mismatches on sklearn >=1.7.
    Xtr64 = np.asarray(Xtr, dtype=np.float64)
    Xte64 = np.asarray(Xte, dtype=np.float64)
    m = LassoCV(cv=5, alphas=50, max_iter=20000, n_jobs=-1, random_state=0)
    m.fit(Xtr64, ytr)
    return m.predict(Xte64)


def fit_tree(Xtr, ytr, Xte, seed: int):
    m = DecisionTreeRegressor(random_state=seed)
    m.fit(Xtr, ytr)
    return m.predict(Xte)


def fit_lgbm(Xtr, ytr, Xte, seed: int):
    m = lgbm.LGBMRegressor(
        objective="regression",
        n_estimators=100,
        learning_rate=0.1,
        random_state=seed,
        verbose=-1,
    )
    m.fit(Xtr, ytr)
    return m.predict(Xte)


# --- core experiment ---------------------------------------------------------
def run_condition(drug: str, conc: float) -> pl.DataFrame:
    df = pl.read_parquet(PARQUET)
    df = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == conc))
    print(f"[{drug} {conc}] rows: {df.height}")

    X_bare, y, labs = build_feature_matrix(df)
    X_pw, _ = add_pairwise_interactions(X_bare, labs)

    n = X_bare.shape[0]
    rng = np.random.default_rng(TEST_SEED)
    perm = rng.permutation(n)
    n_test = int(round(n * TEST_FRAC))
    test_idx = perm[:n_test]
    pool_idx = perm[n_test:]          # 95 % training pool, fixed across everything
    print(f"  fixed test n={len(test_idx)}, training pool n={len(pool_idx)}")

    X_bare_te, y_te = X_bare[test_idx], y[test_idx]
    X_pw_te = X_pw[test_idx]

    rows: list[dict] = []

    # Configuration for each model family: (name, feature_set, fit_fn_name)
    # fit_fn_name picks the right call signature below
    model_specs = [
        ("Linear (additive)", "bare", "linear"),
        ("Lasso (additive)", "bare", "lasso"),
        ("Decision tree", "bare", "tree"),
        ("LightGBM", "bare", "lgbm"),
        ("Linear (+ pairwise)", "pw", "linear"),
        ("Lasso (+ pairwise)", "pw", "lasso"),
    ]

    drug_offset = int.from_bytes(
        hashlib.sha256(f"{drug}_{conc}".encode()).digest()[:4], "little"
    ) % 97
    for frac in TRAIN_FRACTIONS:
        n_train = max(2, int(round(frac * len(pool_idx))))
        for seed in SUBSAMPLE_SEEDS:
            sub_rng = np.random.default_rng(seed + drug_offset)
            sub = sub_rng.choice(pool_idx, size=n_train, replace=False)

            for mname, fset, fn in model_specs:
                Xtr_full = X_bare if fset == "bare" else X_pw
                Xte = X_bare_te if fset == "bare" else X_pw_te
                Xtr = Xtr_full[sub]
                ytr = y[sub]

                # Skip linear (unregularised) when pairwise features outnumber
                # training samples and the system is severely ill-conditioned.
                # We still fit and report NaN if numerically undefined.
                t0 = time.perf_counter()
                try:
                    if fn == "linear":
                        yhat = fit_linear(Xtr, ytr, Xte)
                    elif fn == "lasso":
                        # LassoCV needs cv=5 folds; fall back to smaller cv when training is tiny
                        if n_train < 10:
                            # Fit a plain Lasso at alpha ~= 0 is pointless; skip with NaN
                            yhat = np.full_like(y_te, np.nan, dtype=np.float64)
                        else:
                            yhat = fit_lasso(Xtr, ytr, Xte)
                    elif fn == "tree":
                        yhat = fit_tree(Xtr, ytr, Xte, seed)
                    elif fn == "lgbm":
                        yhat = fit_lgbm(Xtr, ytr, Xte, seed)
                    else:
                        raise ValueError(fn)
                except Exception as e:
                    print(f"    !! {mname} frac={frac} seed={seed} failed: {e}")
                    yhat = np.full_like(y_te, np.nan, dtype=np.float64)
                runtime = time.perf_counter() - t0

                if np.isfinite(yhat).all():
                    rmsd = math.sqrt(mean_squared_error(y_te, yhat))
                    r2 = r2_score(y_te, yhat)
                else:
                    rmsd = float("nan")
                    r2 = float("nan")

                rows.append({
                    "drug": drug,
                    "concentration": conc,
                    "model": mname,
                    "feature_set": fset,
                    "train_fraction_of_pool": frac,
                    "n_train": n_train,
                    "n_test": len(y_te),
                    "seed": seed,
                    "rmsd": rmsd,
                    "r2": r2,
                    "runtime_s": runtime,
                })
            print(f"  frac={frac:.4f} seed={seed} done (n_train={n_train})")

    return pl.DataFrame(rows)


# --- aggregate + plot --------------------------------------------------------
MODEL_ORDER = [
    "Linear (additive)",
    "Lasso (additive)",
    "Decision tree",
    "LightGBM",
    "Linear (+ pairwise)",
    "Lasso (+ pairwise)",
]
MODEL_COLORS = {
    "Linear (additive)": "#7a7a7a",
    "Lasso (additive)": "#4a6fa5",
    "Decision tree": "#d17a22",
    "LightGBM": "#b02a2a",
    "Linear (+ pairwise)": "#2b8a3e",
    "Lasso (+ pairwise)": "#6f42c1",
}
MODEL_LINESTYLE = {
    "Linear (additive)": "-",
    "Lasso (additive)": "-",
    "Decision tree": "-",
    "LightGBM": "-",
    "Linear (+ pairwise)": "--",
    "Lasso (+ pairwise)": "--",
}


def summarise(results: pl.DataFrame) -> pl.DataFrame:
    return (
        results.group_by(["drug", "concentration", "model", "feature_set", "train_fraction_of_pool", "n_train"])
        .agg([
            pl.col("rmsd").mean().alias("rmsd_mean"),
            pl.col("rmsd").std().alias("rmsd_std"),
            pl.col("r2").mean().alias("r2_mean"),
            pl.col("r2").std().alias("r2_std"),
            pl.col("runtime_s").mean().alias("runtime_mean"),
        ])
        .sort(["drug", "model", "train_fraction_of_pool"])
    )


def plot_curve(summary: pl.DataFrame, drug: str, conc: float, out_path: Path) -> None:
    sub = summary.filter(pl.col("drug") == drug)
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), dpi=150)
    ax_rmsd, ax_r2 = axes

    # Determine a sensible RMSD cap for visibility — the Linear(+pw) spike
    # when features > samples can reach ~3 and crushes everything else.
    rmsd_finite = sub["rmsd_mean"].drop_nulls().to_numpy()
    ymax_rmsd = float(np.nanpercentile(rmsd_finite, 92)) * 1.15

    for mname in MODEL_ORDER:
        d = sub.filter(pl.col("model") == mname).sort("n_train")
        if d.is_empty():
            continue
        ntr = d["n_train"].to_numpy()
        rm = d["rmsd_mean"].to_numpy()
        rs = d["rmsd_std"].to_numpy()
        r2 = d["r2_mean"].to_numpy()
        r2s = d["r2_std"].to_numpy()
        ax_rmsd.errorbar(
            ntr, rm, yerr=rs, label=mname,
            color=MODEL_COLORS[mname], linestyle=MODEL_LINESTYLE[mname],
            marker="o", markersize=4, linewidth=1.6, capsize=2, alpha=0.95,
        )
        ax_r2.errorbar(
            ntr, r2, yerr=r2s, label=mname,
            color=MODEL_COLORS[mname], linestyle=MODEL_LINESTYLE[mname],
            marker="o", markersize=4, linewidth=1.6, capsize=2, alpha=0.95,
        )

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Training samples (log scale)")
        ax.grid(True, which="both", alpha=0.25, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax_rmsd.set_ylabel("Test-set RMSD  (log$_{10}$ AUC units)")
    ax_rmsd.set_ylim(bottom=0)
    ax_rmsd.set_ylim(top=ymax_rmsd)
    # Annotate truncation if any point exceeds the cap
    if np.nanmax(rmsd_finite) > ymax_rmsd:
        ax_rmsd.text(
            0.98, 0.97,
            f"axis capped at {ymax_rmsd:.1f};\nunderdetermined OLS+pairwise\nspikes beyond this",
            transform=ax_rmsd.transAxes, ha="right", va="top",
            fontsize=7, color="#555555",
        )

    ax_r2.set_ylabel("Test-set R$^{2}$")
    # R^2 can go below 0 for the unpruned decision tree — use a floor of -0.2
    ax_r2.set_ylim(-0.2, 1.02)
    ax_r2.axhline(0.0, color="#888888", linewidth=0.6, linestyle=":")

    fig.suptitle(
        f"{drug} at {conc:g} µg/mL  —  learning curves (mean ± s.d. over {N_SEEDS} seeds, fixed 5% test set)",
        fontsize=11, y=0.98,
    )
    handles, labels_ = ax_rmsd.get_legend_handles_labels()
    fig.legend(
        handles, labels_, loc="center left", bbox_to_anchor=(0.855, 0.5),
        frameon=False, fontsize=9, title="Model", title_fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 0.84, 0.95))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# --- main --------------------------------------------------------------------
def main():
    all_rows: list[pl.DataFrame] = []
    for drug, conc in CONDITIONS:
        t0 = time.perf_counter()
        res = run_condition(drug, conc)
        print(f"[{drug} {conc}] total time: {time.perf_counter() - t0:.1f} s")
        all_rows.append(res)

    results = pl.concat(all_rows)
    results.write_csv(HERE / "results_table.csv")
    print(f"wrote {HERE / 'results_table.csv'}  ({results.height} rows)")

    summary = summarise(results)
    summary.write_csv(DATADIR / "learning_curves_summary.csv")

    for drug, conc in CONDITIONS:
        out = SUPP / f"figure_s13_{drug.lower()}.png"
        plot_curve(summary, drug, conc, out)


if __name__ == "__main__":
    main()
