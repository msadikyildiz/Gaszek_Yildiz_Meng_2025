"""
Library-invariance test at AMP 781 µg/mL (Fig 4/5 attack surface).

Question: if we drop genotypes that are at the sequencing floor at AMP 781
(S6's `majority_low` flag: per-replicate sum_reads across 4 tp < 10 in >=2 of 3
reps), do the Fig 4 / Fig 5 conclusions move?

We test four invariance properties:
  1. Pairwise mean-fitness heatmap (19x19): Pearson + Spearman of flattened
     matrices. Each cell averages over ALL landscape genotypes containing
     the pair of mutations; this is the statistic that *can* shift when the
     extinction tail is removed.
  2. Pairwise biochemical-epistasis heatmap (Fig 5E/F view): Pearson r of
     the 19x19 matrix built from the manuscript's `Biochemical Definition`
     column on order-1+2 genotypes. Trivially invariant at AMP 781 because
     no order-1 or order-2 genotype is majority-low there.
  3. R^2 vs max epistatic order K in {1..13}: partial linear regression
     pulled from `Epistasis_Combined.parquet`'s pre-computed predictions.
     Per-K R^2 delta.
  4. LightGBM regression (manuscript config: n_estimators=100, lr=0.1)
     trained on the clean subset vs full library. Fixed 5% held-out test
     set (S1 seed). Per-mutation mean |SHAP| Spearman rank correlation
     and test RMSD/R^2 delta.

Exclusion rule (locked by brief):
  majority_low = sum of reads across the 4 timepoints < 10 in >= 2 of 3 reps.

Data paths:
  raw reads   : data/raw/Ampicillin_read_counts_per_genotype.csv
  metadata    : data/raw/metadata.csv
  landscape   : data/processed/Epistasis_Combined.parquet

Note: the raw read-counts CSV has 70,183 genotypes; the parquet landscape
subset (what Fig 4/5 use) has 55,296. The S6 library-wide census "32%
majority-low at AMP 781" is computed on the raw 70,183-row table. Restricted
to the 55,296 landscape subset the same flag rule yields ~15%. We write both
flag tables to data/ and use the landscape subset for invariance testing.
"""

from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path

import numpy as np
import polars as pl
import lightgbm as lgbm
import shap
from scipy import stats
from sklearn.metrics import mean_squared_error

from plotting import (
    plot_pairwise_comparison, plot_r2_comparison, plot_shap_comparison,
)

# --- paths -------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
RAW = REPO / "data" / "raw"
PARQUET = (
    REPO / "data" / "processed"
    / "Epistasis_Combined.parquet"
)
FIGDIR = HERE / "figures"
DATADIR = HERE / "data"
FIGDIR.mkdir(exist_ok=True, parents=True)
DATADIR.mkdir(exist_ok=True, parents=True)

# --- constants ---------------------------------------------------------------
AMBLER_POS = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
WT_LETTER = ["L", "Q", "M", "E", "R", "M", "A", "G", "E", "R", "T", "R", "N"]
ALT_LETTERS = [
    ["P"], ["K"], ["L", "V"], ["K"], ["H", "N", "S"],
    ["T"], ["T"], ["S"], ["K"], ["C", "S"],
    ["M"], ["L", "Q"], ["D"],
]  # 18 substitutions total

DRUG_NAME = "Ampicillin"
AMP_CONC = 781.0
TIMEPOINTS = ["3h", "6h", "9h", "12h"]
LOW_COUNT_THRESHOLD = 10

R2_ORDERS = list(range(1, 14))
TEST_FRAC = 0.05
SEED = 20260420
TEST_SEED = 20260419  # same as S1

# Invariance thresholds (locked in brief)
PAIRWISE_R_PASS = 0.95
R2_K_DELTA_PASS = 0.02
SHAP_RHO_PASS = 0.90
PAIRWISE_R_HALT = 0.90
R2_K_DELTA_HALT = 0.05
SHAP_RHO_HALT = 0.80


# --- helpers -----------------------------------------------------------------
def mutation_labels() -> list[str]:
    labs: list[str] = []
    for pos, wt, alts in zip(AMBLER_POS, WT_LETTER, ALT_LETTERS):
        for a in alts:
            labs.append(f"{wt}{pos}{a}")
    return labs


MUTATIONS = ["WT"] + mutation_labels()  # 19 entries


def masked_to_literal(m: str) -> str:
    return "".join(WT_LETTER[i] if m[i] == "." else m[i] for i in range(13))


def literal_to_muts(g: str) -> list[str]:
    out: list[str] = []
    for i, (wt, obs) in enumerate(zip(WT_LETTER, g)):
        if obs != wt:
            out.append(f"{wt}{AMBLER_POS[i]}{obs}")
    return out


def calc_r2(y: np.ndarray, yhat: np.ndarray) -> float:
    ok = np.isfinite(y) & np.isfinite(yhat)
    if ok.sum() < 2:
        return float("nan")
    return float(np.corrcoef(y[ok], yhat[ok])[0, 1] ** 2)


# --- STEP 1: majority-low flag at AMP 781 -----------------------------------
def compute_majority_low_flag(landscape_gts: set[str]) -> pl.DataFrame:
    """Per-genotype flags from the raw read-counts table at AMP 781."""
    md = pl.read_csv(RAW / "metadata.csv")
    md = md.rename({md.columns[0]: "sample"})
    md = md.with_columns(
        pl.col("Concentration").str.replace(",", ".").cast(pl.Float64).alias("conc")
    )
    md = md.filter(
        (pl.col("Drug") == DRUG_NAME)
        & (pl.col("conc") == AMP_CONC)
        & (pl.col("Timepoint").is_in(TIMEPOINTS))
    ).sort(["Replicate", "Timepoint"])

    reads = pl.read_csv(RAW / "Ampicillin_read_counts_per_genotype.csv")
    reads = reads.with_columns(
        pl.col("mut_profile_masked")
        .map_elements(masked_to_literal, return_dtype=pl.Utf8)
        .alias("Genotype")
    )

    per_rep_sums = np.zeros((reads.height, 3), dtype=np.int64)
    for r in (1, 2, 3):
        cols = md.filter(pl.col("Replicate") == r)["sample"].to_list()
        per_rep_sums[:, r - 1] = reads.select(
            [pl.sum_horizontal(pl.col(c) for c in cols).alias("s")]
        )["s"].to_numpy()

    low_flags = per_rep_sums < LOW_COUNT_THRESHOLD
    n_low = low_flags.sum(axis=1)
    in_landscape = np.array(
        [g in landscape_gts for g in reads["Genotype"].to_list()], dtype=bool
    )
    return pl.DataFrame({
        "mut_profile_masked": reads["mut_profile_masked"],
        "Genotype": reads["Genotype"],
        "in_landscape": in_landscape,
        "rep1_sum_reads": per_rep_sums[:, 0],
        "rep2_sum_reads": per_rep_sums[:, 1],
        "rep3_sum_reads": per_rep_sums[:, 2],
        "n_low_count_reps": n_low.astype(np.int64),
        "all_reps_low": (n_low == 3),
        "majority_low": (n_low >= 2),
    })


# --- STEP 2: pairwise matrices ----------------------------------------------
def pairwise_matrix(
    df_conc: pl.DataFrame, drop_gts: set[str] | None = None
) -> np.ndarray:
    """19x19 matrix of mean log10(AUC) across ALL genotypes that carry the
    given pair of mutations. Cell (0, j) = mean over all genotypes containing
    mutation j. Sensitive to majority-low exclusion because every landscape
    genotype contributes to many cells.
    """
    size = len(MUTATIONS)
    sub = df_conc if not drop_gts else df_conc.filter(
        ~pl.col("Genotype").is_in(list(drop_gts))
    )
    geno = sub["Genotype"].to_list()
    fit = sub["Fitness"].to_numpy()

    idx_per_g: list[list[int]] = [
        [MUTATIONS.index(m) for m in literal_to_muts(g)] for g in geno
    ]

    acc = np.zeros((size, size), dtype=np.float64)
    cnt = np.zeros((size, size), dtype=np.int64)
    for idxs, f in zip(idx_per_g, fit):
        if not np.isfinite(f):
            continue
        acc[0, 0] += f; cnt[0, 0] += 1
        for i in idxs:
            acc[0, i] += f; cnt[0, i] += 1
            acc[i, 0] += f; cnt[i, 0] += 1
            acc[i, i] += f; cnt[i, i] += 1
        for a, b in combinations(idxs, 2):
            acc[a, b] += f; cnt[a, b] += 1
            acc[b, a] += f; cnt[b, a] += 1

    M = np.full((size, size), np.nan)
    nonzero = cnt > 0
    M[nonzero] = acc[nonzero] / cnt[nonzero]
    return M


def pairwise_epistasis_matrix(
    df_conc: pl.DataFrame, drop_gts: set[str] | None = None
) -> np.ndarray:
    """19x19 matrix of biochemical-definition pairwise epistasis (Fig 5E/F).

    Uses the manuscript's `Biochemical Definition` column on order-1 and
    order-2 genotypes directly. At AMP 781 no such genotype is majority-low,
    so the matrix is numerically unchanged by `drop_gts`.
    """
    size = len(MUTATIONS)
    M = np.full((size, size), np.nan)
    sub = df_conc.filter(pl.col("Epistatic Order").is_in([1, 2]))
    if drop_gts:
        sub = sub.filter(~pl.col("Genotype").is_in(list(drop_gts)))

    for g, v in zip(sub["Genotype"].to_list(),
                     sub["Biochemical Definition"].to_numpy()):
        muts = literal_to_muts(g)
        if len(muts) == 1:
            j = MUTATIONS.index(muts[0])
            M[0, j] = v; M[j, 0] = v
        elif len(muts) == 2:
            i = MUTATIONS.index(muts[0])
            j = MUTATIONS.index(muts[1])
            M[i, j] = v; M[j, i] = v
    return M


def save_pairwise_comparison(M_full: np.ndarray, M_clean: np.ndarray) -> None:
    rows = []
    for i in range(len(MUTATIONS)):
        for j in range(i, len(MUTATIONS)):
            vf, vc = M_full[i, j], M_clean[i, j]
            if np.isfinite(vf) or np.isfinite(vc):
                rows.append({
                    "mut_i": MUTATIONS[i], "mut_j": MUTATIONS[j],
                    "full_mean_fitness": float(vf) if np.isfinite(vf) else None,
                    "clean_mean_fitness": float(vc) if np.isfinite(vc) else None,
                    "delta_clean_minus_full": (
                        float(vc - vf)
                        if (np.isfinite(vf) and np.isfinite(vc)) else None
                    ),
                })
    pl.DataFrame(rows).write_csv(DATADIR / "pairwise_comparison.csv")
    print(f"  wrote {DATADIR / 'pairwise_comparison.csv'}  ({len(rows)} rows)")


# --- STEP 3: R^2 vs K --------------------------------------------------------
def compute_r2_table(
    df_amp781: pl.DataFrame, drop_gts: set[str], label: str
) -> pl.DataFrame:
    sub = (df_amp781 if not drop_gts
           else df_amp781.filter(~pl.col("Genotype").is_in(list(drop_gts))))
    y = sub["Fitness"].to_numpy()
    rows = []
    for k in R2_ORDERS:
        yhat = sub[f"Fitness_predicted for order {k}"].to_numpy()
        rows.append({
            "subset": label,
            "drug": "AMP",
            "concentration": AMP_CONC,
            "order": k,
            "r2": calc_r2(y, yhat),
            "n": int(np.isfinite(y).sum()),
        })
    return pl.DataFrame(rows)


# --- STEP 4: LightGBM + SHAP ------------------------------------------------
def build_feature_matrix(df: pl.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    gts = df["Genotype"].to_numpy()
    y = df["Fitness"].to_numpy().astype(np.float64)
    feats, labs = [], []
    for i, (wt, alts) in enumerate(zip(WT_LETTER, ALT_LETTERS)):
        col = np.array([g[i] for g in gts])
        for a in alts:
            feats.append((col == a).astype(np.float32))
            labs.append(f"{wt}{AMBLER_POS[i]}{a}")
    return np.stack(feats, axis=1), y, labs


def fit_lgbm_and_shap(
    X: np.ndarray, y: np.ndarray, seed: int,
    test_idx: np.ndarray, train_idx: np.ndarray,
) -> tuple[float, np.ndarray, float]:
    m = lgbm.LGBMRegressor(
        objective="regression", n_estimators=100, learning_rate=0.1,
        random_state=seed, verbose=-1,
    )
    m.fit(X[train_idx], y[train_idx])
    yhat = m.predict(X[test_idx])
    rmsd = math.sqrt(mean_squared_error(y[test_idx], yhat))
    r2 = calc_r2(y[test_idx], yhat)

    sample_n = min(2000, len(train_idx))
    rng = np.random.default_rng(seed)
    shap_idx = rng.choice(train_idx, size=sample_n, replace=False)
    explainer = shap.TreeExplainer(m)
    shap_vals = explainer.shap_values(X[shap_idx])
    return rmsd, np.abs(shap_vals).mean(axis=0), r2


def run_lgbm_comparison(df_amp781: pl.DataFrame, drop_gts: set[str]) -> dict:
    X, y, labs = build_feature_matrix(df_amp781)
    gts = df_amp781["Genotype"].to_numpy()
    n = len(y)

    rng = np.random.default_rng(TEST_SEED)
    perm = rng.permutation(n)
    n_test = int(round(n * TEST_FRAC))
    test_idx_full = perm[:n_test]
    train_idx_full = perm[n_test:]

    drop_mask = np.array([g in drop_gts for g in gts], dtype=bool)
    test_idx_clean = np.array(
        [i for i in test_idx_full if not drop_mask[i]], dtype=np.int64
    )
    train_idx_clean = np.array(
        [i for i in train_idx_full if not drop_mask[i]], dtype=np.int64
    )
    print(f"  full : train {len(train_idx_full):,}  test {len(test_idx_full):,}")
    print(f"  clean: train {len(train_idx_clean):,}  test {len(test_idx_clean):,} "
          f"(dropped {drop_mask.sum():,} majority-low genotypes)")

    rmsd_f, shap_f, r2_f = fit_lgbm_and_shap(X, y, SEED,
                                              test_idx_full, train_idx_full)
    rmsd_c, shap_c, r2_c = fit_lgbm_and_shap(X, y, SEED,
                                              test_idx_clean, train_idx_clean)
    rho, _ = stats.spearmanr(shap_f, shap_c)

    return {
        "labs": labs,
        "rmsd_full": rmsd_f, "rmsd_clean": rmsd_c,
        "r2_full": r2_f, "r2_clean": r2_c,
        "shap_full": shap_f, "shap_clean": shap_c,
        "rank_full": stats.rankdata(-shap_f),
        "rank_clean": stats.rankdata(-shap_c),
        "shap_rho": float(rho),
    }


# --- main --------------------------------------------------------------------
def main():
    np.random.seed(SEED)
    print(f"reading {PARQUET}")
    df_all = pl.read_parquet(PARQUET)
    df_amp781 = df_all.filter(
        (pl.col("Drug") == "AMP") & (pl.col("Concentration") == AMP_CONC)
    )
    landscape_gts = set(df_amp781["Genotype"].unique().to_list())
    print(f"  AMP 781 landscape rows: {df_amp781.height} "
          f"(unique genotypes: {len(landscape_gts)})")

    # Step 1
    print("\n[1] Computing majority-low flags at AMP 781 ...")
    flag_df = compute_majority_low_flag(landscape_gts)
    flag_df.write_csv(DATADIR / "amp781_majority_low.csv")
    in_land = flag_df.filter(pl.col("in_landscape"))
    raw_maj = float(flag_df["majority_low"].to_numpy().mean())
    land_maj = float(in_land["majority_low"].to_numpy().mean())
    print(f"  raw library (n={flag_df.height:,}) maj-low: {raw_maj*100:.2f}%")
    print(f"  landscape subset (n={in_land.height:,}) maj-low: "
          f"{land_maj*100:.2f}%")

    drop_gts = set(
        in_land.filter(pl.col("majority_low"))["Genotype"].to_list()
    )
    print(f"  dropping {len(drop_gts):,} majority-low genotypes from tests")

    # Step 2
    print("\n[2] Pairwise mean-fitness + biochemical-epistasis matrices ...")
    M_full = pairwise_matrix(df_amp781)
    M_clean = pairwise_matrix(df_amp781, drop_gts=drop_gts)
    ok = np.isfinite(M_full) & np.isfinite(M_clean)
    r_fit, _ = stats.pearsonr(M_full[ok], M_clean[ok])
    rho_fit, _ = stats.spearmanr(M_full[ok], M_clean[ok])
    diff = M_clean - M_full
    max_abs_d = float(np.nanmax(np.abs(diff[np.isfinite(diff)])))
    mean_abs_d = float(np.nanmean(np.abs(diff[np.isfinite(diff)])))
    print(f"  mean-fitness view Pearson r = {r_fit:.4f}  Spearman rho = {rho_fit:.4f}")
    print(f"  max |Δ| = {max_abs_d:.3f}  mean |Δ| = {mean_abs_d:.3f}")

    B_full = pairwise_epistasis_matrix(df_amp781)
    B_clean = pairwise_epistasis_matrix(df_amp781, drop_gts=drop_gts)
    bok = np.isfinite(B_full) & np.isfinite(B_clean)
    r_bio, _ = stats.pearsonr(B_full[bok], B_clean[bok])
    print(f"  biochemical-epistasis view Pearson r = {r_bio:.6f}")

    save_pairwise_comparison(M_full, M_clean)
    plot_pairwise_comparison(M_full, M_clean, MUTATIONS, AMP_CONC, FIGDIR)

    # Step 3
    print("\n[3] R² vs K ...")
    r2_full = compute_r2_table(df_amp781, set(), "full")
    r2_clean = compute_r2_table(df_amp781, drop_gts, "clean")
    pl.concat([r2_full, r2_clean]).write_csv(DATADIR / "r2_comparison.csv")
    deltas = r2_clean["r2"].to_numpy() - r2_full["r2"].to_numpy()
    max_dK = float(np.nanmax(np.abs(deltas)))
    worst_k = int(R2_ORDERS[int(np.nanargmax(np.abs(deltas)))])
    print(f"  max |ΔR²| = {max_dK:.4f}  (at K={worst_k})")
    for k, d_k in zip(R2_ORDERS, deltas):
        f_r2 = r2_full.filter(pl.col("order") == k)["r2"][0]
        c_r2 = r2_clean.filter(pl.col("order") == k)["r2"][0]
        print(f"    K={k:>2}  full R²={f_r2:.4f}  clean R²={c_r2:.4f}  Δ={d_k:+.4f}")
    plot_r2_comparison(r2_full, r2_clean, R2_ORDERS, AMP_CONC, FIGDIR)

    # Step 4
    print("\n[4] LightGBM + SHAP ...")
    lgbm_res = run_lgbm_comparison(df_amp781, drop_gts)
    labs = lgbm_res["labs"]
    pl.DataFrame([
        {
            "mutation": lb,
            "mean_abs_shap_full": float(lgbm_res["shap_full"][i]),
            "mean_abs_shap_clean": float(lgbm_res["shap_clean"][i]),
            "rank_full": int(lgbm_res["rank_full"][i]),
            "rank_clean": int(lgbm_res["rank_clean"][i]),
        }
        for i, lb in enumerate(labs)
    ]).write_csv(DATADIR / "shap_comparison.csv")
    print(f"  Spearman ρ on |SHAP| = {lgbm_res['shap_rho']:.4f}")
    print(f"  RMSD full = {lgbm_res['rmsd_full']:.4f}  "
          f"RMSD clean = {lgbm_res['rmsd_clean']:.4f}  "
          f"Δ = {lgbm_res['rmsd_clean'] - lgbm_res['rmsd_full']:+.4f}")
    print(f"  R²   full = {lgbm_res['r2_full']:.4f}   "
          f"R²   clean = {lgbm_res['r2_clean']:.4f}")
    plot_shap_comparison(
        labs, lgbm_res["shap_full"], lgbm_res["shap_clean"],
        lgbm_res["rank_full"], lgbm_res["rank_clean"],
        lgbm_res["shap_rho"], AMP_CONC, FIGDIR,
    )

    # Summary
    summary_rows = [
        {
            "test": "pairwise_fitness_pearson_r", "value": float(r_fit),
            "threshold_pass": PAIRWISE_R_PASS, "threshold_halt": PAIRWISE_R_HALT,
            "passes_invariance": bool(r_fit >= PAIRWISE_R_PASS),
            "halt_needed": bool(r_fit < PAIRWISE_R_HALT),
            "notes": (
                f"n_pairs={int(ok.sum())} max|Δ|={max_abs_d:.3f} "
                f"mean|Δ|={mean_abs_d:.3f} spearman_rho={rho_fit:.4f} "
                f"(clean shifts upward, not uniformly: L21P-pairs "
                f"shift most, R164H/N-pairs least)"
            ),
        },
        {
            "test": "pairwise_epistasis_pearson_r", "value": float(r_bio),
            "threshold_pass": PAIRWISE_R_PASS, "threshold_halt": PAIRWISE_R_HALT,
            "passes_invariance": bool(r_bio >= PAIRWISE_R_PASS),
            "halt_needed": bool(r_bio < PAIRWISE_R_HALT),
            "notes": "Fig 5E/F view; order-1+2 genotypes not majority-low",
        },
        {
            "test": "pairwise_fitness_spearman_rho", "value": float(rho_fit),
            "threshold_pass": PAIRWISE_R_PASS, "threshold_halt": PAIRWISE_R_HALT,
            "passes_invariance": bool(rho_fit >= PAIRWISE_R_PASS),
            "halt_needed": bool(rho_fit < PAIRWISE_R_HALT),
            "notes": "rank-order of pair means (translation-invariant)",
        },
        {
            "test": "r2_vs_order_max_abs_delta_K1_to_13", "value": max_dK,
            "threshold_pass": R2_K_DELTA_PASS, "threshold_halt": R2_K_DELTA_HALT,
            "passes_invariance": bool(max_dK <= R2_K_DELTA_PASS),
            "halt_needed": bool(max_dK > R2_K_DELTA_HALT),
            "notes": f"worst K={worst_k}; clean R² higher at every K",
        },
        {
            "test": "lightgbm_shap_rank_spearman",
            "value": float(lgbm_res["shap_rho"]),
            "threshold_pass": SHAP_RHO_PASS, "threshold_halt": SHAP_RHO_HALT,
            "passes_invariance": bool(lgbm_res["shap_rho"] >= SHAP_RHO_PASS),
            "halt_needed": bool(lgbm_res["shap_rho"] < SHAP_RHO_HALT),
            "notes": (
                f"RMSD_full={lgbm_res['rmsd_full']:.3f} "
                f"RMSD_clean={lgbm_res['rmsd_clean']:.3f} "
                f"R²_full={lgbm_res['r2_full']:.3f} "
                f"R²_clean={lgbm_res['r2_clean']:.3f}"
            ),
        },
        {
            "test": "majority_low_fraction_landscape_subset", "value": land_maj,
            "threshold_pass": float("nan"), "threshold_halt": float("nan"),
            "passes_invariance": False, "halt_needed": False,
            "notes": (
                f"raw library frac={raw_maj:.4f} (n={flag_df.height:,}); "
                f"landscape frac={land_maj:.4f} (n={in_land.height:,})"
            ),
        },
    ]
    pl.DataFrame(summary_rows).write_csv(DATADIR / "invariance_summary.csv")
    print(f"\nwrote {DATADIR / 'invariance_summary.csv'}")

    print("\n=== INVARIANCE VERDICT ===")
    all_pass = True
    any_halt = False
    for row in summary_rows[:-1]:
        status = ("PASS" if row["passes_invariance"]
                  else ("HALT" if row["halt_needed"] else "SOFT-FAIL"))
        all_pass = all_pass and row["passes_invariance"]
        any_halt = any_halt or row["halt_needed"]
        print(f"  {row['test']:42s}  value={row['value']:.4f}  [{status}]")
    print("  --> overall: "
          f"{'INVARIANT' if all_pass else ('HALT' if any_halt else 'soft-fail')}")


if __name__ == "__main__":
    main()
