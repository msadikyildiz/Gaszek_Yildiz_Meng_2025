"""Supplementary Figures S7/S8 - monoculture IC50 validates pooled AUC-fitness.

Produces the two Supplementary panels backing the wet-lab validation of the
pooled-fitness landscape: (S7) IC50 vs mean AUC-fitness scatter with Pearson
r + p, one panel per drug, and (S8) per-concentration panels for AMP (batch
20260407) and AZT (batch 20260307). Palette mirrors
analysis/s09_s10_epistatic_order/analysis.py and
si_figures/s03_dose_response/analysis.py -- AMP = greys, AZT = RdPu,
DejaVu Sans font, Illustrator-editable PDF fonts.

Provenance
----------
Ported from an earlier internal development version of this script (not part
of this public repository). That script reads pre-computed processed
parquets under
    contributors/deniz_validation_experiments/src/processed/<batch>/xref_expanded_df.parquet
The same processed tables live in this public repo at
    validation/src/processed/<batch>/xref_expanded_df.parquet
(identical schema: drug, variant, genotype_13, mean_fitness, log10_ic50,
fitness_<concentration> columns), and regenerate from the raw plate-reader
XLSX via `validation/src/run_all.py`. Nothing under `validation/` is modified
by this script; it only reads the already-processed parquets.

Numbers reproduced (cross-checked against validation/S7_S8_provenance.md,
which independently derives the same Pearson r/n from the same parquets via
a different, minimal snippet):
    AMP (batch 20260407): n = 13, Pearson r = 0.885, p = 5.9e-05
    AZT (batch 20260307): n = 18, Pearson r = 0.803, p = 6.1e-05

Scope note: the earlier version also built a `before_after_comparison.png`
(side-by-side against Ilona's pre-replot `ic50_auc_dotplot.png` diagnostic
PNG, used only for that version's own visual QA). That comparison is
intentionally NOT ported here -- it is not one of the three published
Supplementary panels, and its "before" input is an internal diagnostic image
that isn't part of this repository's tracked/regenerable outputs.

Outputs (written directly to the published Supplementary location):
    figures/supplementary/figure_s07.png       (S7: IC50 vs mean AUC-fitness)
    figures/supplementary/figure_s08_amp.png   (S8, AMP per-concentration)
    figures/supplementary/figure_s08_azt.png   (S8, AZT per-concentration)
A convenience copy of each (same content) is also written locally under
si_figures/s07_s08_ic50/figures/ using the earlier version's original
filenames, for provenance diffing.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy import stats

# --- paths -------------------------------------------------------------------
HERE = Path(__file__).resolve().parent


def _repo_root(p: Path) -> Path:
    """Walk up from this file to the repo root (dir with the manuscript parquet)."""
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


REPO = _repo_root(Path(__file__).resolve())
PROCESSED = REPO / "validation" / "src" / "processed"

FIGDIR = HERE / "figures"
SUPP = REPO / "figures" / "supplementary"
FIGDIR.mkdir(parents=True, exist_ok=True)
SUPP.mkdir(parents=True, exist_ok=True)

# --- constants ---------------------------------------------------------------
AMP_BATCH = "20260407"
AZT_BATCH = "20260307"

AMP_CONCS = [3.1, 12.2, 48.8, 195.0, 781.0]
AZT_CONCS = [0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]
DRUG_CONCS = {"AMP": AMP_CONCS, "AZT": AZT_CONCS}
DRUG_NAME = {"AMP": "Ampicillin", "AZT": "Aztreonam"}

# Manuscript palette (matches analysis/s09_s10_epistatic_order/analysis.py).
AMP_PRIMARY = "#2f2f2f"
AMP_SECONDARY = "#7a7a7a"
AZT_PRIMARY = "#c2185b"
AZT_SECONDARY = "#e898b9"
DRUG_PRIMARY = {"AMP": AMP_PRIMARY, "AZT": AZT_PRIMARY}
DRUG_SECONDARY = {"AMP": AMP_SECONDARY, "AZT": AZT_SECONDARY}
WT_EDGE = "#000000"
WT_FACE = "#ffffff"

# Concentration ramps (mirror the s09_s10 rewrite).
AMP_RAMP = plt.cm.Greys(np.linspace(0.35, 0.92, len(AMP_CONCS)))
AZT_RAMP = plt.cm.RdPu(np.linspace(0.35, 0.92, len(AZT_CONCS)))
CONC_COLOR = {("AMP", c): AMP_RAMP[i] for i, c in enumerate(AMP_CONCS)}
CONC_COLOR.update({("AZT", c): AZT_RAMP[i] for i, c in enumerate(AZT_CONCS)})

# --- matplotlib globals ------------------------------------------------------
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False


# --- helpers -----------------------------------------------------------------
def load_xref(drug: str) -> pl.DataFrame:
    """Load the xref_expanded_df for the latest batch of the given drug."""
    batch = AMP_BATCH if drug == "AMP" else AZT_BATCH
    path = PROCESSED / batch / "xref_expanded_df.parquet"
    return (
        pl.read_parquet(path)
        .filter(pl.col("drug") == drug)
        # Keep WT separately; DD can be noisy at bottom — keep it but colour softly.
    )


def pearson_r(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    """Pearson r, two-sided p, n on finite pairs."""
    ok = np.isfinite(x) & np.isfinite(y)
    x_ok, y_ok = x[ok], y[ok]
    if len(x_ok) < 3:
        return float("nan"), float("nan"), int(ok.sum())
    r, p = stats.pearsonr(x_ok, y_ok)
    return float(r), float(p), int(ok.sum())


def format_ic50_tick(val: float) -> str:
    """Log10(IC50) tick -> human readable ug/mL."""
    ic = 10 ** val
    if ic >= 1000:
        return f"{ic / 1000:.0f}k"
    if ic >= 10:
        return f"{ic:.0f}"
    if ic >= 1:
        return f"{ic:.1f}"
    if ic >= 0.1:
        return f"{ic:.2f}"
    return f"{ic:.3f}"


def sig_stars(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def marker_style_for(variant: str, drug: str) -> dict:
    """Return marker kwargs for a variant.

    WT / DD get open-circle (they are reference genotypes in manuscript).
    All other variants get filled circles in the drug's primary colour.
    """
    if variant == "WT":
        return dict(
            facecolors=WT_FACE, edgecolors=DRUG_PRIMARY[drug],
            linewidths=1.8, s=80, marker="o", zorder=5,
        )
    if variant == "DD":
        return dict(
            facecolors=WT_FACE, edgecolors=DRUG_PRIMARY[drug],
            linewidths=1.2, s=70, marker="o",
            linestyle="--", zorder=4,
        )
    return dict(
        facecolors=DRUG_PRIMARY[drug], edgecolors="#fafafa",
        linewidths=1.0, s=70, marker="o", zorder=4,
    )


def plot_variant_points(
    ax, df: pl.DataFrame, drug: str, x_col: str, y_col: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Scatter all variants, draw dashed edge for DD, open circle for WT.

    Also return x,y arrays for the regression / r calculation.
    """
    xs = df[x_col].to_numpy()
    ys = df[y_col].to_numpy()
    variants = df["variant"].to_list()

    for x, y, v in zip(xs, ys, variants):
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        kw = marker_style_for(v, drug).copy()
        # scatter doesn't accept linestyle directly for edges, split handling.
        ls = kw.pop("linestyle", "-")
        sc = ax.scatter([x], [y], **kw)
        if ls == "--":
            # Matplotlib scatter edge style: emulate with circle patch.
            sc.set_linestyle("--")

    return xs, ys


def _nudge_labels_2d(
    x: np.ndarray, y: np.ndarray, x_range: float, y_range: float,
    dx_init: float = 0.04, dy_init: float = -0.03,
    passes: int = 30, repel_radius: float = 0.11,
    point_repel_radius: float = 0.06,
) -> tuple[np.ndarray, np.ndarray]:
    """Iterative label repulsion (label-label + label-point + boundary clamp).

    Self-contained re-implementation of the label-nudging behaviour in
    validation/src/config_theme.py's nudge_labels, inlined here so this
    script has no import dependency on the validation/ package internals.
    """
    lx = x + dx_init * x_range
    ly = y + dy_init * y_range
    x_min, x_max = x.min() - 0.04 * x_range, x.max() + 0.04 * x_range
    y_min, y_max = y.min() - 0.04 * y_range, y.max() + 0.04 * y_range

    for _ in range(passes):
        # label-label
        for i in range(len(lx)):
            for j in range(i + 1, len(lx)):
                dxn = (lx[i] - lx[j]) / x_range
                dyn = (ly[i] - ly[j]) / y_range
                d = np.sqrt(dxn * dxn + dyn * dyn)
                if 0 < d < repel_radius:
                    f = (repel_radius - d) / 2
                    lx[i] += f * dxn / d * x_range
                    ly[i] += f * dyn / d * y_range
                    lx[j] -= f * dxn / d * x_range
                    ly[j] -= f * dyn / d * y_range
        # label-point
        for i in range(len(lx)):
            for px, py in zip(x, y):
                dxn = (lx[i] - px) / x_range
                dyn = (ly[i] - py) / y_range
                d = np.sqrt(dxn * dxn + dyn * dyn)
                if 0 < d < point_repel_radius:
                    f = (point_repel_radius - d) * 0.5
                    lx[i] += f * dxn / d * x_range
                    ly[i] += f * dyn / d * y_range
        lx = np.clip(lx, x_min, x_max)
        ly = np.clip(ly, y_min, y_max)
    return lx, ly


def add_variant_labels(
    ax, df: pl.DataFrame, x_col: str, y_col: str,
    skip: set[str] | None = None,
    nudge: bool = True,
) -> None:
    """Annotate points with short variant names; optionally nudge for overlap."""
    skip = skip or set()
    xs = df[x_col].to_numpy()
    ys = df[y_col].to_numpy()
    variants = df["variant"].to_list()

    keep = [
        (x, y, v) for x, y, v in zip(xs, ys, variants)
        if np.isfinite(x) and np.isfinite(y) and v not in skip
    ]
    if not keep:
        return
    xs_k = np.array([k[0] for k in keep])
    ys_k = np.array([k[1] for k in keep])
    vs_k = [k[2] for k in keep]

    x_range = float(xs_k.max() - xs_k.min()) or 1.0
    y_range = float(ys_k.max() - ys_k.min()) or 1.0

    if nudge:
        lx, ly = _nudge_labels_2d(xs_k.copy(), ys_k.copy(), x_range, y_range)
        # Connector lines from point to nudged label when offset is substantial.
        for x, y, v, lxi, lyi in zip(xs_k, ys_k, vs_k, lx, ly):
            dxn = (lxi - x) / x_range
            dyn = (lyi - y) / y_range
            if np.sqrt(dxn * dxn + dyn * dyn) > 0.04:
                ax.plot(
                    [x, lxi], [y, lyi],
                    color="#bdbdbd", linewidth=0.5, alpha=0.5, zorder=1,
                )
            ax.text(
                lxi, lyi, v,
                fontsize=7.4, color="#333333", ha="center", va="center",
                zorder=6,
            )
    else:
        dx = 0.02 * x_range
        dy = 0.025 * y_range
        for x, y, v in zip(xs_k, ys_k, vs_k):
            ax.annotate(
                v, xy=(x, y), xytext=(x + dx, y + dy),
                fontsize=7.2, color="#333333", ha="left", va="bottom",
                annotation_clip=False,
            )


def add_regression(
    ax, x: np.ndarray, y: np.ndarray, color: str, significant: bool,
) -> None:
    """Dashed regression line across the plot domain."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return
    x_ok, y_ok = x[ok], y[ok]
    slope, intercept = np.polyfit(x_ok, y_ok, 1)
    x_line = np.array([x_ok.min(), x_ok.max()])
    y_line = slope * x_line + intercept
    ax.plot(
        x_line, y_line,
        linestyle="--", linewidth=1.4 if significant else 1.0,
        color=color if significant else "#bdbdbd",
        alpha=0.8 if significant else 0.45,
        zorder=2,
    )


def format_log10_ic50_axis(ax, axis: str = "y") -> None:
    """Set ticks/labels as 'ug/mL' strings even though values are log10."""
    ticks_fn = ax.set_yticks if axis == "y" else ax.set_xticks
    label_fn = ax.set_yticklabels if axis == "y" else ax.set_xticklabels

    # Pick integer-ish log10 ticks in the current range.
    limits = ax.get_ylim() if axis == "y" else ax.get_xlim()
    lo, hi = float(limits[0]), float(limits[1])
    candidate_logs = np.arange(np.floor(lo) - 1, np.ceil(hi) + 1)
    ticks = [t for t in candidate_logs if lo <= t <= hi]
    if len(ticks) == 0:
        ticks = list(np.linspace(lo, hi, 5))
    ticks_fn(ticks)
    label_fn([format_ic50_tick(t) for t in ticks])


# --- Figure S7: IC50 vs mean AUC-fitness (two drug panels) -------------------
def plot_ic50_vs_auc(
    amp_df: pl.DataFrame, azt_df: pl.DataFrame, out_paths: list[Path],
) -> dict:
    """Two-panel scatter: log10(IC50) (y) vs mean AUC-fitness (x) for AMP and AZT.

    Pearson r + p, n annotated per panel. Regression line dashed.
    """
    fig = plt.figure(figsize=(9.5, 4.6), dpi=150)
    gs = fig.add_gridspec(1, 2, wspace=0.32)
    ax_amp = fig.add_subplot(gs[0, 0])
    ax_azt = fig.add_subplot(gs[0, 1])

    stats_out: dict = {}

    for ax, df, drug in [(ax_amp, amp_df, "AMP"), (ax_azt, azt_df, "AZT")]:
        # Keep only rows with both fitness and IC50 defined (so DD which has
        # null mean_fitness is excluded from correlation but still plotted).
        sub = df.filter(pl.col("genotype_13").is_not_null())

        x_all = sub["mean_fitness"].to_numpy()
        y_all = sub["log10_ic50"].to_numpy()

        r, p, n = pearson_r(x_all, y_all)
        stats_out[drug] = {"r": r, "p": p, "n": n}

        plot_variant_points(
            ax, sub, drug, x_col="mean_fitness", y_col="log10_ic50",
        )
        # Skip labelling WT/DD because their marker style is already a legend.
        add_variant_labels(
            ax, sub, x_col="mean_fitness", y_col="log10_ic50",
            skip={"WT", "DD"},
        )

        # Padding.
        x_min, x_max = float(np.nanmin(x_all)), float(np.nanmax(x_all))
        y_min, y_max = float(np.nanmin(y_all)), float(np.nanmax(y_all))
        x_pad = (x_max - x_min) * 0.10 + 0.05
        y_pad = (y_max - y_min) * 0.10 + 0.10
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)

        add_regression(
            ax, x_all, y_all, color=DRUG_PRIMARY[drug], significant=(p < 0.05),
        )

        # r,p text annotation.
        stars = sig_stars(p)
        txt = f"Pearson r = {r:.2f}{'  ' + stars if stars else ''}\n"
        if p < 1e-4:
            txt += f"p = {p:.1e},  n = {n}"
        else:
            txt += f"p = {p:.3g},  n = {n}"
        ax.text(
            0.98, 0.04, txt,
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9.5, color=DRUG_PRIMARY[drug], fontweight="bold",
        )

        ax.set_xlabel(
            "mean AUC-fitness across drug concentrations",
            fontsize=10.5,
        )
        ax.set_ylabel("monoculture IC$_{50}$ (µg/mL)", fontsize=10.5)
        ax.set_title(
            f"{DRUG_NAME[drug]}", fontsize=12, loc="left",
            color=DRUG_PRIMARY[drug], fontweight="bold", pad=6,
        )
        ax.grid(True, alpha=0.2, linewidth=0.5)

        # Reformat IC50 ticks to human-readable concentrations.
        format_log10_ic50_axis(ax, axis="y")

        # Minimal reference-genotype legend (matches Figure 3B convention).
        legend_handles = [
            plt.Line2D(
                [0], [0], marker="o", color="none",
                markerfacecolor="none", markeredgecolor=DRUG_PRIMARY[drug],
                markersize=8, markeredgewidth=1.6,
                label="TEM-1$^{WT}$",
            ),
            plt.Line2D(
                [0], [0], marker="o", color="none",
                markerfacecolor=DRUG_PRIMARY[drug], markeredgecolor="#fafafa",
                markersize=7, markeredgewidth=0.8,
                label="TEM-1 variant",
            ),
        ]
        ax.legend(
            handles=legend_handles, loc="upper left",
            fontsize=8.5, frameon=False, handlelength=1.2,
            borderpad=0.2, labelspacing=0.3,
        )

    fig.suptitle(
        "Monoculture IC$_{50}$ validates pooled AUC-fitness",
        fontsize=13, fontweight="bold", y=1.02, x=0.5,
    )

    for out_path in out_paths:
        fig.savefig(out_path, dpi=350, bbox_inches="tight")
        print(f"saved {out_path}")
    plt.close(fig)
    return stats_out


# --- Figure S8: per-concentration AUC-fitness vs IC50 (one drug per file) ----
def plot_ic50_per_conc(
    df: pl.DataFrame, drug: str, out_paths: list[Path],
) -> dict:
    """Grid of panels: one per selection concentration.

    y-axis: log10(IC50) (shared across panels)
    x-axis: AUC-fitness at that concentration
    Marker colour: concentration ramp (so the reader can quickly see which
                   concentration each panel reflects).
    """
    concs = DRUG_CONCS[drug]
    n = len(concs)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.0 * ncols, 3.2 * nrows),
        dpi=150, squeeze=False,
    )

    sub = df.filter(pl.col("genotype_13").is_not_null())

    # Common y (log10 IC50) range.
    y_all = sub["log10_ic50"].to_numpy()
    y_lo = float(np.nanmin(y_all))
    y_hi = float(np.nanmax(y_all))
    y_pad = (y_hi - y_lo) * 0.12 + 0.2

    # Find panel-agnostic x-axis range for consistency (optional — each panel
    # uses its own range to make small signal at low conc visible).
    per_panel_stats: dict = {}

    for i, conc in enumerate(concs):
        ax = axes[i // ncols, i % ncols]
        fit_col = f"fitness_{conc:g}"
        # Col name uses trailing .0 for some concentrations (e.g. fitness_195.0)
        # and matches what polars rendered; verify:
        if fit_col not in sub.columns:
            # fall back to actual string, which includes the ".0" or similar.
            # Build by scanning columns that start with "fitness_".
            matches = [c for c in sub.columns if c.startswith("fitness_")]
            num_map = {float(c.replace("fitness_", "")): c for c in matches}
            if conc in num_map:
                fit_col = num_map[conc]
            elif float(conc) in num_map:
                fit_col = num_map[float(conc)]
            else:
                raise KeyError(f"No fitness column for conc={conc}; have {matches}")

        x_all = sub[fit_col].to_numpy()
        r, p, nobs = pearson_r(x_all, y_all)
        per_panel_stats[conc] = {"r": r, "p": p, "n": nobs}

        # Scatter with per-panel conc colour (so the mini-bar legend can show
        # which panel corresponds to which concentration in the colour ramp).
        panel_color = CONC_COLOR[(drug, conc)]
        for x, y, v in zip(x_all, y_all, sub["variant"].to_list()):
            if not (np.isfinite(x) and np.isfinite(y)):
                continue
            if v == "WT":
                ax.scatter(
                    x, y, facecolors=WT_FACE, edgecolors=DRUG_PRIMARY[drug],
                    linewidths=1.8, s=60, marker="o", zorder=5,
                )
            elif v == "DD":
                ax.scatter(
                    x, y, facecolors=WT_FACE, edgecolors=DRUG_PRIMARY[drug],
                    linewidths=1.0, s=55, marker="o", zorder=4,
                )
            else:
                ax.scatter(
                    x, y, facecolors=panel_color, edgecolors="#fafafa",
                    linewidths=0.8, s=55, marker="o", zorder=4,
                )

        add_regression(
            ax, x_all, y_all,
            color=DRUG_PRIMARY[drug], significant=(p < 0.05),
        )

        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)
        x_min, x_max = float(np.nanmin(x_all)), float(np.nanmax(x_all))
        x_pad = (x_max - x_min) * 0.10 + 0.05
        ax.set_xlim(x_min - x_pad, x_max + x_pad)

        # Annotation.
        stars = sig_stars(p)
        txt = f"r = {r:.2f}{'  ' + stars if stars else ''}"
        ax.text(
            0.03, 0.97, txt,
            transform=ax.transAxes, ha="left", va="top",
            fontsize=9, fontweight="bold",
            color=DRUG_PRIMARY[drug] if p < 0.05 else "#777777",
        )

        ax.set_title(
            f"{conc:g} µg/mL",
            fontsize=10, loc="left", color=DRUG_PRIMARY[drug], fontweight="bold",
        )
        ax.grid(True, alpha=0.2, linewidth=0.5)
        ax.tick_params(axis="both", labelsize=8)

        # Only outer axes get labels.
        if i % ncols == 0:
            ax.set_ylabel(
                "monoculture IC$_{50}$ (µg/mL)", fontsize=9,
            )
            format_log10_ic50_axis(ax, axis="y")
        else:
            ax.set_yticklabels([])
        if i // ncols == nrows - 1:
            ax.set_xlabel("AUC-fitness at this conc.", fontsize=9)

    # Hide empty axes if any.
    for j in range(n, nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    # Title.
    fig.suptitle(
        f"{DRUG_NAME[drug]} — monoculture IC$_{{50}}$ vs pooled AUC-fitness "
        f"at each selection concentration",
        fontsize=12, fontweight="bold", y=0.99, x=0.5,
    )

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for out_path in out_paths:
        fig.savefig(out_path, dpi=350, bbox_inches="tight")
        print(f"saved {out_path}")
    plt.close(fig)
    return per_panel_stats


# --- main --------------------------------------------------------------------
def main() -> None:
    print(f"Loading processed data from {PROCESSED}...")
    amp_df = load_xref("AMP")
    azt_df = load_xref("AZT")
    print(f"  AMP: {amp_df.height} variants (batch {AMP_BATCH})")
    print(f"  AZT: {azt_df.height} variants (batch {AZT_BATCH})")

    print("\n1/3 - figure_s07.png")
    main_stats = plot_ic50_vs_auc(
        amp_df, azt_df,
        [SUPP / "figure_s07.png", FIGDIR / "ext_ic50_vs_auc.png"],
    )
    for d, s in main_stats.items():
        print(f"     {d}: r = {s['r']:.3f}, p = {s['p']:.3e}, n = {s['n']}")

    print("\n2/3 - figure_s08_amp.png")
    amp_panel_stats = plot_ic50_per_conc(
        amp_df, "AMP",
        [SUPP / "figure_s08_amp.png", FIGDIR / "ext_ic50_per_conc_amp.png"],
    )
    for c, s in amp_panel_stats.items():
        print(f"     AMP @ {c:g}: r = {s['r']:.3f}, p = {s['p']:.3g}")

    print("\n3/3 - figure_s08_azt.png")
    azt_panel_stats = plot_ic50_per_conc(
        azt_df, "AZT",
        [SUPP / "figure_s08_azt.png", FIGDIR / "ext_ic50_per_conc_azt.png"],
    )
    for c, s in azt_panel_stats.items():
        print(f"     AZT @ {c:g}: r = {s['r']:.3f}, p = {s['p']:.3g}")


if __name__ == "__main__":
    main()
