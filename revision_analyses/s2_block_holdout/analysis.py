"""
Block-holdout ML validation for Gaszek, Yildiz, Meng (2026) rebuttal (Reviewer #4 ii).

Runs three train/test regimes on the 55,296-genotype TEM-1 landscape:
  1. Leave-one-mutation-out (LOMO): exclude every variant carrying a focal
     substitution, train LightGBM on the remainder, evaluate on held-out set.
  2. Hamming-distance holdout: train on variants within Hamming radius D of a
     reference genotype, test on those outside.  References: WT (canonical),
     and a random deep-mutant (sensitivity check).
  3. Matched-size random split baseline to quantify the interpolation gap.

LightGBM config matches the manuscript (n_estimators=100, learning_rate=0.1,
objective='regression').  Features: 18-column one-hot of non-WT substitutions.

Outputs
-------
  figures/block_holdout_amp.png
  figures/block_holdout_azt.png
  data/lomo_results.csv
  data/hamming_results.csv
  data/random_baseline.csv
  results_table.csv             (aggregated deliverable)
"""
from __future__ import annotations
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/claude/mpl_cache')
os.environ.setdefault('FONTCONFIG_CACHE', '/tmp/claude/fc_cache')

import numpy as np
import pandas as pd
import lightgbm as lgbm
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ paths
def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / 'data' / 'processed' / 'Epistasis_Combined.parquet'
OUT = Path(__file__).resolve().parent
FIG = OUT / 'figures'
DAT = OUT / 'data'
FIG.mkdir(parents=True, exist_ok=True)
DAT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ schema
# Ambler-numbered mutation labels (paper convention)
MUTATION_LABELS = [
    'L21P', 'Q39K', 'M69L', 'M69V', 'E104K',
    'R164H', 'R164N', 'R164S', 'M182T', 'A237T',
    'G238S', 'E240K', 'R244C', 'R244S', 'T265M',
    'R275L', 'R275Q', 'N276D',
]
# Position (0-based in 13-char genotype) and mutant letter for each
POS_LETTER = [
    (0, 'P'), (1, 'K'), (2, 'L'), (2, 'V'), (3, 'K'),
    (4, 'H'), (4, 'N'), (4, 'S'), (5, 'T'), (6, 'T'),
    (7, 'S'), (8, 'K'), (9, 'C'), (9, 'S'), (10, 'M'),
    (11, 'L'), (11, 'Q'), (12, 'D'),
]
WT = 'LQMERMAGERTRN'

DRUG_SPEC = {
    'AMP': ('AMP', 781.0),
    'AZT': ('AZT', 36.0),
}

LGBM_KW = dict(
    objective='regression',
    n_estimators=100,
    learning_rate=0.1,
    random_state=42,
    verbose=-1,
    n_jobs=-1,
)

# ------------------------------------------------------------------ helpers
def load_drug(drug: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return one-hot X (N x 18), y (N,), and Hamming-from-WT d (N,)."""
    tag, conc = DRUG_SPEC[drug]
    df = pd.read_parquet(PARQUET)
    df = df[(df.Drug == tag) & (df.Concentration == conc)].reset_index(drop=True)
    genotypes = df.Genotype.values
    y = df.Fitness.values.astype(np.float64)

    # Stack genotypes to (N,13) char matrix
    G = np.stack([np.frombuffer(g.encode(), dtype='S1') for g in genotypes])
    # One-hot: X[:, k] = 1 if position POS_LETTER[k][0] has letter POS_LETTER[k][1]
    X = np.zeros((len(genotypes), len(POS_LETTER)), dtype=np.float32)
    for k, (p, letter) in enumerate(POS_LETTER):
        X[:, k] = (G[:, p] == letter.encode()).astype(np.float32)

    WT_bytes = np.frombuffer(WT.encode(), dtype='S1')
    d = (G != WT_bytes).sum(axis=1).astype(np.int32)
    return X, y, d, genotypes


def hamming_from(genotypes: np.ndarray, ref: str) -> np.ndarray:
    G = np.stack([np.frombuffer(g.encode(), dtype='S1') for g in genotypes])
    R = np.frombuffer(ref.encode(), dtype='S1')
    return (G != R).sum(axis=1).astype(np.int32)


def fit_eval(X_tr, y_tr, X_te, y_te, model_fn=None):
    model = model_fn() if model_fn else lgbm.LGBMRegressor(**LGBM_KW)
    model.fit(X_tr, y_tr)
    pred = model.predict(X_te)
    rmsd = math.sqrt(mean_squared_error(y_te, pred))
    r2 = r2_score(y_te, pred)
    return rmsd, r2


# ------------------------------------------------------------------ LOMO
def run_lomo(X, y, drug, include_lasso=True):
    rows = []
    for k, label in enumerate(MUTATION_LABELS):
        mask_test = X[:, k] == 1  # variants carrying the focal mutation
        mask_train = ~mask_test
        n_test = int(mask_test.sum())
        n_train = int(mask_train.sum())

        # Drop the focal column so the model can't see it; it's constant per side.
        keep_cols = [j for j in range(X.shape[1]) if j != k]
        X_tr = X[mask_train][:, keep_cols]
        X_te = X[mask_test][:, keep_cols]
        y_tr = y[mask_train]
        y_te = y[mask_test]

        # LightGBM
        rmsd, r2 = fit_eval(X_tr, y_tr, X_te, y_te)

        # Matched-N random baseline: randomly sample n_train from the whole set,
        # test on the complement (size n_test-ish -- we just use disjoint split).
        X_rtr, X_rte, y_rtr, y_rte = train_test_split(
            X, y, train_size=n_train, random_state=42, shuffle=True,
        )
        rmsd_rand, r2_rand = fit_eval(X_rtr, y_rtr, X_rte, y_rte)

        row = dict(
            drug=drug, mutation=label, position=POS_LETTER[k][0],
            n_train=n_train, n_test=n_test,
            rmsd_lomo=rmsd, r2_lomo=r2,
            rmsd_random=rmsd_rand, r2_random=r2_rand,
        )
        if include_lasso:
            rmsd_l, r2_l = fit_eval(
                X_tr, y_tr, X_te, y_te,
                model_fn=lambda: Lasso(alpha=1e-3, max_iter=20000),
            )
            row.update(rmsd_lomo_lasso=rmsd_l, r2_lomo_lasso=r2_l)
        rows.append(row)
        print(f'  [{drug} LOMO] {label}: n_train={n_train} n_test={n_test} '
              f'RMSD={rmsd:.3f} (rand={rmsd_rand:.3f}) R²={r2:.3f}')
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ Hamming
def run_hamming(X, y, d, drug, label='WT'):
    rows = []
    for D in range(2, 11):
        mask_train = d < D
        mask_test = d >= D
        n_train = int(mask_train.sum())
        n_test = int(mask_test.sum())
        if n_train < 50 or n_test < 50:
            rows.append(dict(drug=drug, reference=label, D=D,
                             n_train=n_train, n_test=n_test,
                             rmsd=np.nan, r2=np.nan,
                             rmsd_random=np.nan, r2_random=np.nan))
            print(f'  [{drug} Hamming ref={label}] D={D}: SKIP (n_train={n_train}, n_test={n_test})')
            continue

        rmsd, r2 = fit_eval(X[mask_train], y[mask_train], X[mask_test], y[mask_test])

        # Matched random baseline
        X_rtr, X_rte, y_rtr, y_rte = train_test_split(
            X, y, train_size=n_train, random_state=42, shuffle=True,
        )
        rmsd_rand, r2_rand = fit_eval(X_rtr, y_rtr, X_rte, y_rte)
        rows.append(dict(drug=drug, reference=label, D=D,
                         n_train=n_train, n_test=n_test,
                         rmsd=rmsd, r2=r2,
                         rmsd_random=rmsd_rand, r2_random=r2_rand))
        print(f'  [{drug} Hamming ref={label}] D={D}: n_train={n_train} n_test={n_test} '
              f'RMSD={rmsd:.3f} (rand={rmsd_rand:.3f}) R²={r2:.3f}')
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ main
def main():
    np.random.seed(42)
    print('=== Loading drugs ===')
    data = {}
    for drug in ('AMP', 'AZT'):
        X, y, d, genotypes = load_drug(drug)
        data[drug] = dict(X=X, y=y, d=d, genotypes=genotypes)
        print(f'{drug}: X={X.shape}, y mean={y.mean():.3f}, d median={np.median(d)}')

    # Pick a random deep genotype for the sensitivity-check reference.
    # Seed-based so it is reproducible; we want a genotype ~7 mutations from WT
    # (near the mode of the Hamming distribution).
    genotypes_amp = data['AMP']['genotypes']
    d_wt = data['AMP']['d']
    pool = np.where(d_wt == 7)[0]
    rng = np.random.default_rng(42)
    ref_idx = rng.choice(pool)
    REF_DEEP = genotypes_amp[ref_idx]
    print(f'random-deep reference genotype = {REF_DEEP} (Hamming-to-WT = 7)')

    # ---- LOMO ----
    print('\n=== LOMO (leave-one-mutation-out) ===')
    lomo_rows = []
    for drug in ('AMP', 'AZT'):
        d = data[drug]
        print(f'-- {drug} --')
        part = run_lomo(d['X'], d['y'], drug, include_lasso=True)
        lomo_rows.append(part)
    lomo = pd.concat(lomo_rows, ignore_index=True)
    lomo.to_csv(DAT / 'lomo_results.csv', index=False)

    # ---- Hamming ----
    print('\n=== Hamming-distance holdout ===')
    ham_rows = []
    for drug in ('AMP', 'AZT'):
        dd = data[drug]
        print(f'-- {drug}, reference=WT --')
        ham_rows.append(run_hamming(dd['X'], dd['y'], dd['d'], drug, label='WT'))
        # Sensitivity: deep reference
        d_deep = hamming_from(dd['genotypes'], REF_DEEP)
        print(f'-- {drug}, reference=random deep ({REF_DEEP}) --')
        ham_rows.append(run_hamming(dd['X'], dd['y'], d_deep, drug,
                                    label=f'deep[{REF_DEEP}]'))
    ham = pd.concat(ham_rows, ignore_index=True)
    ham.to_csv(DAT / 'hamming_results.csv', index=False)

    # ---- Aggregated summary ----
    summary_rows = []
    for drug in ('AMP', 'AZT'):
        sub = lomo[lomo.drug == drug]
        summary_rows.append(dict(
            scope='LOMO_summary', drug=drug, detail='LightGBM',
            mean_rmsd=sub.rmsd_lomo.mean(), median_rmsd=sub.rmsd_lomo.median(),
            max_rmsd=sub.rmsd_lomo.max(),
            worst_mutation=sub.loc[sub.rmsd_lomo.idxmax(), 'mutation'],
            mean_random_rmsd=sub.rmsd_random.mean(),
        ))
        if 'rmsd_lomo_lasso' in sub.columns:
            summary_rows.append(dict(
                scope='LOMO_summary', drug=drug, detail='Lasso',
                mean_rmsd=sub.rmsd_lomo_lasso.mean(),
                median_rmsd=sub.rmsd_lomo_lasso.median(),
                max_rmsd=sub.rmsd_lomo_lasso.max(),
                worst_mutation=sub.loc[sub.rmsd_lomo_lasso.idxmax(), 'mutation'],
                mean_random_rmsd=sub.rmsd_random.mean(),
            ))
    pd.concat([
        pd.DataFrame(summary_rows),
        lomo.assign(scope='LOMO_per_mutation'),
        ham.assign(scope='Hamming_holdout'),
    ], ignore_index=True).to_csv(OUT / 'results_table.csv', index=False)

    # ---- Figures ----
    for drug in ('AMP', 'AZT'):
        make_figure(drug, lomo[lomo.drug == drug], ham[ham.drug == drug])

    print('\nDone. Outputs in', OUT)


def make_figure(drug, lomo_df, ham_df):
    lomo_df = lomo_df.sort_values('rmsd_lomo', ascending=False).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), dpi=150,
                             gridspec_kw=dict(width_ratios=[1.35, 1]))

    # --- Panel A: per-mutation LOMO bars ---
    ax = axes[0]
    y_pos = np.arange(len(lomo_df))
    bar_h = 0.38
    ax.barh(y_pos - bar_h/2, lomo_df.rmsd_lomo.values, height=bar_h,
            color='#D55E00', edgecolor='black', linewidth=0.4, label='LOMO')
    ax.barh(y_pos + bar_h/2, lomo_df.rmsd_random.values, height=bar_h,
            color='#999999', edgecolor='black', linewidth=0.4,
            label='random split (matched N)')
    # Overlay Lasso LOMO if present
    if 'rmsd_lomo_lasso' in lomo_df.columns:
        ax.scatter(lomo_df.rmsd_lomo_lasso.values, y_pos - bar_h/2,
                   color='#0072B2', s=22, zorder=5, label='Lasso LOMO', marker='d')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(lomo_df.mutation.values, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('RMSD (held-out fitness)', fontsize=10)
    ax.set_title(f'{drug}: Leave-one-mutation-out (LOMO)', fontsize=11)
    ax.grid(axis='x', linestyle='--', alpha=0.35)
    ax.legend(fontsize=8, loc='lower right', frameon=True)
    # annotate held-out N, flushed to right edge to avoid clashing with
    # the Lasso-LOMO diamond markers that sit on the LOMO-bar row.
    x_right = max(lomo_df.rmsd_lomo.max(),
                  lomo_df.get('rmsd_lomo_lasso', lomo_df.rmsd_lomo).max())
    x_text = x_right * 1.04
    for i, row in lomo_df.iterrows():
        ax.text(x_text, i, f'n={int(row.n_test)}',
                va='center', ha='left', fontsize=6.5, color='#444')
    ax.set_xlim(0, x_right * 1.18)

    # --- Panel B: Hamming-distance curve (WT reference, with sensitivity) ---
    ax = axes[1]
    for ref, marker, color in [
        ('WT', 'o', '#D55E00'),
        (ham_df.reference[ham_df.reference.str.startswith('deep', na=False)].iloc[0]
            if (ham_df.reference.astype(str).str.startswith('deep')).any() else None,
         's', '#0072B2'),
    ]:
        if ref is None:
            continue
        sub = ham_df[ham_df.reference == ref].sort_values('D')
        sub = sub.dropna(subset=['rmsd'])
        ax.plot(sub.D.values, sub.rmsd.values, '-' + marker, color=color,
                label=f'block (ref={ref.split("[")[0]})', markersize=5, linewidth=1.5)
        ax.plot(sub.D.values, sub.rmsd_random.values, '--' + marker, color=color,
                alpha=0.55, label=f'random (ref={ref.split("[")[0]})',
                markersize=4, linewidth=1.0)

    ax.set_xlabel('Minimum Hamming distance D\n(train: d<D, test: d≥D)', fontsize=10)
    ax.set_ylabel('RMSD (held-out fitness)', fontsize=10)
    ax.set_title(f'{drug}: Distance-stratified holdout', fontsize=11)
    ax.grid(linestyle='--', alpha=0.35)
    ax.legend(fontsize=7, loc='best', frameon=True)

    plt.tight_layout()
    out_path = FIG / f'block_holdout_{drug.lower()}.png'
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out_path}')


if __name__ == '__main__':
    main()
