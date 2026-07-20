"""
Extra check: does the 'RMSD ~0.3 is sufficient to classify resistant/sensitive
for clinical purposes' claim survive the block-holdout regime?

For each LOMO fold and for the matched-N random baseline we compute binary
classification metrics at a biologically reasonable cutoff.  Thresholds are
per-drug percentiles (AMP: P25 loss-of-resistance; AZT: P90 gain-of-resistance,
see DRUG_THRESHOLDS below). These were chosen because the raw median is
uninformative for AZT (WT sits in the bulk; the resistant tail lives above P90).

Outputs: appended rows to results_table.csv, and two JSON snippets for the
results.md narrative.
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
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, roc_auc_score,
    matthews_corrcoef, confusion_matrix,
)

# paths + schema (copied minimally)
def _repo_root(_p):
    """Walk up from this file to the repo root (dir with data/processed/Epistasis_Combined.parquet)."""
    for _a in (_p, *_p.parents):
        if (_a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return _a
    raise FileNotFoundError("repo root not found from " + str(_p))
REPO = _repo_root(Path(__file__).resolve())
PARQUET = REPO / 'data' / 'processed' / 'Epistasis_Combined.parquet'
OUT = Path(__file__).resolve().parent

MUTATION_LABELS = [
    'L21P', 'Q39K', 'M69L', 'M69V', 'E104K',
    'R164H', 'R164N', 'R164S', 'M182T', 'A237T',
    'G238S', 'E240K', 'R244C', 'R244S', 'T265M',
    'R275L', 'R275Q', 'N276D',
]
POS_LETTER = [
    (0, 'P'), (1, 'K'), (2, 'L'), (2, 'V'), (3, 'K'),
    (4, 'H'), (4, 'N'), (4, 'S'), (5, 'T'), (6, 'T'),
    (7, 'S'), (8, 'K'), (9, 'C'), (9, 'S'), (10, 'M'),
    (11, 'L'), (11, 'Q'), (12, 'D'),
]
DRUG_SPEC = {'AMP': ('AMP', 781.0), 'AZT': ('AZT', 36.0)}
LGBM_KW = dict(objective='regression', n_estimators=100, learning_rate=0.1,
               random_state=42, verbose=-1, n_jobs=-1)


def load_drug(drug):
    tag, conc = DRUG_SPEC[drug]
    df = pd.read_parquet(PARQUET)
    df = df[(df.Drug == tag) & (df.Concentration == conc)].reset_index(drop=True)
    gs = df.Genotype.values
    y = df.Fitness.values.astype(np.float64)
    G = np.stack([np.frombuffer(g.encode(), dtype='S1') for g in gs])
    X = np.zeros((len(gs), len(POS_LETTER)), dtype=np.float32)
    for k, (p, letter) in enumerate(POS_LETTER):
        X[:, k] = (G[:, p] == letter.encode()).astype(np.float32)
    return X, y


def cls_metrics(y_true, y_pred, thr, direction='above'):
    """direction: 'above' -> class 1 = y>=thr; 'below' -> class 1 = y<=thr."""
    if direction == 'above':
        y_true_b = (y_true >= thr).astype(int)
        y_pred_b = (y_pred >= thr).astype(int)
        score = y_pred  # higher = more likely class 1
    else:
        y_true_b = (y_true <= thr).astype(int)
        y_pred_b = (y_pred <= thr).astype(int)
        score = -y_pred
    tn, fp, fn, tp = confusion_matrix(y_true_b, y_pred_b, labels=[0, 1]).ravel()
    auroc = (roc_auc_score(y_true_b, score)
             if len(np.unique(y_true_b)) > 1 else np.nan)
    return dict(
        accuracy=accuracy_score(y_true_b, y_pred_b),
        balanced_accuracy=balanced_accuracy_score(y_true_b, y_pred_b),
        mcc=matthews_corrcoef(y_true_b, y_pred_b),
        auroc=auroc,
        sensitivity=tp / max(1, tp + fn),
        specificity=tn / max(1, tn + fp),
        threshold=thr,
        class1_rate_true=float(y_true_b.mean()),
    )


DRUG_THRESHOLDS = {
    # AMP: WT is fully resistant (fitness 4.36, above P99); the clinically
    # meaningful classification is "loss-of-resistance" variants.  Use P25
    # as the cutoff (below = sensitive / destabilized).
    'AMP': ('loss_of_resistance_P25', 25),
    # AZT: WT (2.21) is near the bulk; the clinically meaningful label is
    # "gain-of-function / ESBL-like resistance".  Use P90 as the cutoff.
    'AZT': ('gain_of_resistance_P90', 90),
}


def main():
    rows = []
    for drug in ('AMP', 'AZT'):
        X, y = load_drug(drug)
        label, pct = DRUG_THRESHOLDS[drug]
        thr = float(np.percentile(y, pct))
        print(f'\n== {drug}: threshold = P{pct} fitness = {thr:.3f} ({label}) ==')
        # random baseline (90/10 split like manuscript)
        Xtr, Xte, ytr, yte = train_test_split(X, y, train_size=0.9, random_state=42)
        m = lgbm.LGBMRegressor(**LGBM_KW); m.fit(Xtr, ytr); yp = m.predict(Xte)
        direction = 'below' if drug == 'AMP' else 'above'
        metrics = cls_metrics(yte, yp, thr, direction=direction)
        metrics.update(drug=drug, scheme='random_split_90_10', mutation='(all)',
                       n_test=len(yte), rmsd=math.sqrt(np.mean((yte - yp) ** 2)))
        rows.append(metrics); print('  random 90/10:', metrics)

        # LOMO each mutation
        for k, label in enumerate(MUTATION_LABELS):
            mask_te = X[:, k] == 1
            keep = [j for j in range(X.shape[1]) if j != k]
            Xtr = X[~mask_te][:, keep]; ytr = y[~mask_te]
            Xte = X[mask_te][:, keep]; yte = y[mask_te]
            m = lgbm.LGBMRegressor(**LGBM_KW); m.fit(Xtr, ytr); yp = m.predict(Xte)
            direction = 'below' if drug == 'AMP' else 'above'
            metrics = cls_metrics(yte, yp, thr, direction=direction)
            metrics.update(drug=drug, scheme='LOMO', mutation=label,
                           n_test=len(yte), rmsd=math.sqrt(np.mean((yte - yp) ** 2)))
            rows.append(metrics)
            print(f'  LOMO {label}: acc={metrics["accuracy"]:.3f} '
                  f'bal={metrics["balanced_accuracy"]:.3f} '
                  f'MCC={metrics["mcc"]:.3f} AUROC={metrics["auroc"]:.3f}')

    out = pd.DataFrame(rows)
    out.to_csv(OUT / 'data' / 'classification_results.csv', index=False)

    # Summarize
    print('\n== Summary (balanced accuracy) ==')
    for drug in ('AMP', 'AZT'):
        lomo = out[(out.drug == drug) & (out.scheme == 'LOMO')]
        rand = out[(out.drug == drug) & (out.scheme == 'random_split_90_10')].iloc[0]
        print(f'{drug}: random={rand["balanced_accuracy"]:.3f} | '
              f'LOMO median={lomo["balanced_accuracy"].median():.3f} '
              f'min={lomo["balanced_accuracy"].min():.3f} '
              f'worst={lomo.loc[lomo["balanced_accuracy"].idxmin(), "mutation"]}')


if __name__ == '__main__':
    main()
