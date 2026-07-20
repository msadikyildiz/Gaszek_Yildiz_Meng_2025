"""Plot helpers for the AMP 781 library-invariance test."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"


def plot_pairwise_comparison(
    M_full: np.ndarray, M_clean: np.ndarray,
    mutations: list[str], amp_conc: float, figdir: Path,
) -> Path:
    """Side-by-side mean-fitness heatmaps + difference panel."""
    diff = M_clean - M_full
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.6), dpi=150)
    vmax = float(np.nanmax(np.concatenate([
        M_full[np.isfinite(M_full)], M_clean[np.isfinite(M_clean)]
    ])))
    dv = float(np.nanmax(np.abs(diff[np.isfinite(diff)])))
    dv = max(dv, 1e-3)

    labels = ["none" if m == "WT" else m for m in mutations]
    panels = [
        (axes[0], M_full, "full library", "viridis", 0, vmax),
        (axes[1], M_clean, "non-majority-low subset", "viridis", 0, vmax),
        (axes[2], diff, "clean − full", "RdBu_r", -dv, dv),
    ]
    for ax, M, ttl, cmap, lo, hi in panels:
        im = ax.imshow(M, cmap=cmap, vmin=lo, vmax=hi, origin="upper")
        ax.set_title(ttl, fontsize=10)
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=5.4)
        ax.set_yticklabels(labels, fontsize=5.4)
        ax.tick_params(length=1.5, pad=1)
        fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)

    ok = np.isfinite(M_full) & np.isfinite(M_clean)
    r, _ = stats.pearsonr(M_full[ok], M_clean[ok])
    axes[2].text(
        0.02, 0.98, f"Pearson r = {r:.4f}",
        transform=axes[2].transAxes, ha="left", va="top",
        fontsize=9, color="#222",
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2),
    )

    fig.suptitle(
        f"AMP {amp_conc:g} µg/mL — pairwise mean-fitness matrices (full vs clean)",
        fontsize=11, y=1.02,
    )
    out = figdir / "pairwise_full_vs_clean_amp781.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    return out


def plot_r2_comparison(
    full: pl.DataFrame, clean: pl.DataFrame, r2_orders: list[int],
    amp_conc: float, figdir: Path,
) -> Path:
    fig, ax = plt.subplots(1, 1, figsize=(6.0, 3.8), dpi=150)
    f = full.sort("order")
    c = clean.sort("order")
    ax.plot(
        f["order"].to_numpy(), f["r2"].to_numpy(),
        marker="o", markersize=5, linewidth=1.8,
        color="#222", label=f"full library (n={f['n'][0]:,})",
    )
    ax.plot(
        c["order"].to_numpy(), c["r2"].to_numpy(),
        marker="s", markersize=5, linewidth=1.8,
        color="#b02a2a", label=f"non-majority-low (n={c['n'][0]:,})",
    )
    ax.set_xticks(r2_orders)
    ax.set_xlim(0.7, 13.3)
    ax.set_ylim(0.0, 1.02)
    ax.set_xlabel("highest-order epistatic terms included")
    ax.set_ylabel("R² (measured vs. predicted fitness)")
    ax.grid(True, alpha=0.25, linewidth=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="lower right", fontsize=8, frameon=False)

    max_delta = float(np.nanmax(np.abs(
        c["r2"].to_numpy() - f["r2"].to_numpy()
    )))
    ax.set_title(
        f"AMP {amp_conc:g} µg/mL — R² vs order  (max |Δ| = {max_delta:.3f})",
        fontsize=10,
    )

    out = figdir / "r2_vs_order_full_vs_clean_amp781.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    return out


def plot_shap_comparison(
    labs: list[str],
    shap_full: np.ndarray, shap_clean: np.ndarray,
    rank_full: np.ndarray, rank_clean: np.ndarray,
    shap_rho: float, amp_conc: float, figdir: Path,
) -> Path:
    order = np.argsort(-shap_full)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), dpi=150)

    ax = axes[0]
    x = np.arange(len(labs))
    w = 0.42
    ax.bar(x - w/2, shap_full[order], width=w, color="#222",
           label="full library")
    ax.bar(x + w/2, shap_clean[order], width=w, color="#b02a2a",
           label="non-majority-low")
    ax.set_xticks(x)
    ax.set_xticklabels([labs[i] for i in order], rotation=75, fontsize=7)
    ax.set_ylabel("mean |SHAP|")
    ax.set_title("mean |SHAP| per mutation (ordered by full)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = axes[1]
    ax.scatter(rank_full, rank_clean, s=28, color="#4a6fa5", edgecolor="white")
    for i, lb in enumerate(labs):
        ax.annotate(lb, (rank_full[i], rank_clean[i]), fontsize=6,
                    textcoords="offset points", xytext=(3, 3))
    mx = len(labs) + 1
    ax.plot([0, mx], [0, mx], "k--", linewidth=0.7, alpha=0.5)
    ax.set_xlabel("|SHAP| rank, full library")
    ax.set_ylabel("|SHAP| rank, non-majority-low")
    ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.25, linewidth=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"Spearman ρ = {shap_rho:.3f}", fontsize=10)

    fig.suptitle(
        f"AMP {amp_conc:g} µg/mL — LightGBM SHAP, full vs clean", fontsize=11,
    )
    out = figdir / "shap_full_vs_clean_amp781.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    return out
