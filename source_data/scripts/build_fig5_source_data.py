#!/usr/bin/env python3
"""Fixed-seed Source Data for Figure 5 (SHAP interpretation of TEM-1 fitness).

Figure 5 SHAP values come from the src/04 (AMP) and src/03 (AZT) regression pipelines. This
script recomputes that pipeline with a fixed seed so the SHAP arrays behind each panel are
reproducible, and reports that the per-mutation attributions are stable across split seeds
(per-mutation mean|SHAP| Spearman rho and max drift across seeds are written out).

Pipeline (identical to the notebooks): Epistasis_Combined at the reference dose (AMP 781 /
AZT 36) -> 18 one-hot mutation-presence features -> LGBMRegressor(n_estimators=100,
learning_rate=0.1, random_state=42) trained on a 10% split -> shap.TreeExplainer over the 90%
test set. Panels: 5A/5B heatmap = top-553 test genotypes by predicted fitness (18 x 553);
beeswarm = all test genotypes; 5C = mean|SHAP| per mutation/position.
"""
import csv
from pathlib import Path

import numpy as np
import polars as pl
import lightgbm as lgbm
import shap
from scipy import stats
from sklearn.model_selection import train_test_split

HERE = Path(__file__).resolve().parent


def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))


REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / "data/processed/Epistasis_Combined.parquet"
OUT = HERE.parent / "derived"
OUT.mkdir(exist_ok=True)

REF = {"AMP": 781.0, "AZT": 36.0}
SEED = 42                       # fixed split+model seed for the shipped Source Data
STAB_SEEDS = [42, 7, 123, 2024]  # split seeds for the robustness check
WT_SEQ = "LQMERMAGERTRN"
# (feature name, string index, alt residue, Ambler position) — the 18 notebook features.
FEATURES = [
    ("L21P", 0, "P", 21), ("Q39K", 1, "K", 39), ("M69L", 2, "L", 69), ("M69V", 2, "V", 69),
    ("E104K", 3, "K", 104), ("R164H", 4, "H", 164), ("R164N", 4, "N", 164), ("R164S", 4, "S", 164),
    ("M182T", 5, "T", 182), ("A237T", 6, "T", 237), ("G238S", 7, "S", 238), ("E240K", 8, "K", 240),
    ("R244C", 9, "C", 244), ("R244S", 9, "S", 244), ("T265M", 10, "M", 265), ("R275L", 11, "L", 275),
    ("R275Q", 11, "Q", 275), ("N276D", 12, "D", 276),
]
NAMES = [f[0] for f in FEATURES]


def write_csv(name, header, rows):
    with open(OUT / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {name}  ({len(rows)} rows)")


def design_matrix(df, drug):
    sub = df.filter((pl.col("Drug") == drug) & (pl.col("Concentration") == REF[drug]))
    genos = sub["Genotype"].to_list()
    y = np.asarray(sub["Fitness"].to_list(), dtype=float)
    X = np.zeros((len(genos), len(FEATURES)), dtype=np.int8)
    for j, (_, idx, alt, _) in enumerate(FEATURES):
        X[:, j] = [1 if g[idx] == alt else 0 for g in genos]
    return np.array(genos), X, y


def fit_shap(X, y, seed):
    Xtr, Xte, ytr, yte, itr, ite = train_test_split(
        X, y, np.arange(len(y)), train_size=0.1, random_state=seed)
    model = lgbm.LGBMRegressor(objective="regression", n_estimators=100,
                               learning_rate=0.1, random_state=42, verbose=-1)
    model.fit(Xtr, ytr)
    ypred = model.predict(Xte)
    sv = shap.TreeExplainer(model)(Xte).values          # (n_test, 18)
    return ite, yte, ypred, sv


def build():
    df = pl.read_parquet(PARQUET)
    stab_rows, c5_rows = [], []
    for drug in ("AMP", "AZT"):
        genos, X, y = design_matrix(df, drug)
        ite, yte, ypred, sv = fit_shap(X, y, SEED)
        gte = genos[ite]

        # 5A/5B: top-553 test genotypes by predicted fitness, most -> least resistant
        top = np.argsort(ypred)[-553:][::-1]
        rows = [[r + 1, gte[i], round(float(ypred[i]), 4), round(float(yte[i]), 4)]
                + [round(float(sv[i, j]), 5) for j in range(len(FEATURES))]
                for r, i in enumerate(top)]
        write_csv(f"fig5{'A' if drug=='AMP' else 'B'}_shap_top553_{drug.lower()}.csv",
                  ["rank", "genotype_13", "predicted_fitness", "measured_fitness"] + NAMES, rows)

        # beeswarm master: all test-set SHAP (companion)
        allrows = [[gte[i], round(float(ypred[i]), 4), round(float(yte[i]), 4)]
                   + [round(float(sv[i, j]), 5) for j in range(len(FEATURES))]
                   for i in range(len(gte))]
        write_csv(f"fig5_shap_all_test_{drug.lower()}.csv",
                  ["genotype_13", "predicted_fitness", "measured_fitness"] + NAMES, allrows)

        # 5C: mean|SHAP| and mean SHAP per mutation
        mabs = np.abs(sv).mean(axis=0)
        msign = sv.mean(axis=0)
        for j, (nm, _, _, amb) in enumerate(FEATURES):
            c5_rows.append([drug, nm, amb, round(float(mabs[j]), 5), round(float(msign[j]), 5)])

        # robustness: per-mutation mean|SHAP| across seeds
        per_seed = {SEED: mabs}
        for s in STAB_SEEDS:
            if s == SEED:
                continue
            _, _, _, sv_s = fit_shap(X, y, s)
            per_seed[s] = np.abs(sv_s).mean(axis=0)
        base = per_seed[SEED]
        rhos, drifts = [], []
        for s in STAB_SEEDS:
            if s == SEED:
                continue
            rho = stats.spearmanr(base, per_seed[s]).statistic
            drift = float(np.max(np.abs(base - per_seed[s])))
            rhos.append(rho); drifts.append(drift)
        print(f"  {drug} seed-robustness: Spearman rho(mean|SHAP|) vs seed42 = "
              f"{['%.4f'%r for r in rhos]}, max drift = {max(drifts):.4f}")
        for j, (nm, _, _, amb) in enumerate(FEATURES):
            stab_rows.append([drug, nm, amb] + [round(float(per_seed[s][j]), 5) for s in STAB_SEEDS])

    write_csv("fig5C_mean_abs_shap.csv",
              ["drug", "mutation", "ambler_position", "mean_abs_shap", "mean_shap"], c5_rows)
    write_csv("fig5_shap_seed_stability.csv",
              ["drug", "mutation", "ambler_position"] + [f"mean_abs_shap_seed{s}" for s in STAB_SEEDS],
              stab_rows)


if __name__ == "__main__":
    build()
