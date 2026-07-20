"""
S6 — Revised Fig. S3 for Reviewer #1 comment on dose-response line noise.

Reviewer #1: "Figure S3 shows that some genotypes have very noisy fitness values,
jumping from the minimum to the mean trend from one concentration to another.
What is the reason for this? Does this occur infrequently enough that it does
not affect results such as landscape graphs and epistasis measurements?"

Response strategy (agreed with PI):
  - Revise Fig S3 by replacing the all-library grey lines with a small curated
    set of highlighted variants:
      (1) wild type (LQMERMAGERTRN, black);
      (2) TEM-1_dead (XXXXXXXXXXXXX, blue-grey);
      (3) the 18 single-substitution mutants (one per non-WT column), colored
          by Ambler position;
      (4) the clinical TEM alleles in the library (yellow triangles on Fig 3A);
      (5) the top-N composite variants by mean fitness at the highest drug
          concentration (orange-red).
  - Annotate the extinction artifact by per-timepoint read counts. Any
    (genotype, drug, conc, replicate) with sum_reads < 10 across its four
    timepoints is flagged as low-count; plotted as an open marker with reduced
    alpha. The mean ± SD band is computed over non-flagged replicates only.
  - Shade the empirical AUC floor (the TEM-1_dead band) as the limit of
    detection (LOD).

Inputs (PROJECT/Gaszek_Yildiz_Meng_2025/data/raw):
  Ampicillin_auc_per_genotype.csv,  Aztreonam_auc_per_genotype.csv
  Ampicillin_read_counts_per_genotype.csv, Aztreonam_read_counts_per_genotype.csv
  metadata.csv    (sample -> Drug/Concentration/Replicate/Timepoint)
  data/known_variants/encoded_variants.csv  (63 clinical TEM alleles)

Outputs:
  figures/fig_s3_revised_amp.png, figures/fig_s3_revised_azt.png,
  figures/fig_s3_revised.png  (combined two-panel)
  data/highlighted_variants.csv
  data/per_variant_fitness.csv
  results_table.csv   (per-conc summary: n_flagged, fraction_flagged, LOD)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
CLINICAL_CSV = (
    REPO / "data" / "known_variants" / "encoded_variants.csv"
)
FIGDIR = HERE / "figures"
DATADIR = HERE / "data"
SUPP = REPO / "figures" / "supplementary"
FIGDIR.mkdir(exist_ok=True, parents=True)
DATADIR.mkdir(exist_ok=True, parents=True)
SUPP.mkdir(exist_ok=True, parents=True)

# --- constants ---------------------------------------------------------------
# The library encodes genotypes in "masked" form: '.' at a position where the
# residue matches wild-type, the substitution letter otherwise. Wild-type is
# therefore 13 dots; the catalytically dead control is 13 X's.
WT_LITERAL = "LQMERMAGERTRN"   # the raw residue string for wild-type TEM-1
WT_GENOTYPE = "." * 13         # "masked" encoding used throughout the library
DEAD_GENOTYPE = "XXXXXXXXXXXXX"
AMBLER_POS = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]
WT_LETTER = ["L", "Q", "M", "E", "R", "M", "A", "G", "E", "R", "T", "R", "N"]
ALT_LETTERS = [
    ["P"], ["K"], ["L", "V"], ["K"], ["H", "N", "S"],
    ["T"], ["T"], ["S"], ["K"], ["C", "S"],
    ["M"], ["L", "Q"], ["D"],
]   # 18 substitutions total

AMP_CONCS = [0.0, 3.1, 12.2, 48.8, 195.0, 781.0]
AZT_CONCS = [0.0, 0.44, 1.33, 4.0, 12.0, 36.0, 108.0, 324.0]
DRUG_CONCS = {"Ampicillin": AMP_CONCS, "Aztreonam": AZT_CONCS}

LOW_COUNT_THRESHOLD = 10                # sum_reads across 4 tp per rep
N_TOP_COMPOSITE = 20
TIMEPOINTS = ["3h", "6h", "9h", "12h"]
SEED = 20260420

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]


# --- IO ----------------------------------------------------------------------
def load_metadata() -> pl.DataFrame:
    """metadata.csv with concentration coerced to float (handles '3,100' -> 3.1)."""
    md = pl.read_csv(RAW / "metadata.csv")
    md = md.rename({md.columns[0]: "sample"})
    md = md.with_columns([
        pl.col("Concentration").str.replace(",", ".").cast(pl.Float64).alias("conc"),
    ])
    return md


def load_auc(drug: str) -> pl.DataFrame:
    """Per-genotype per-replicate log10(AUC) table."""
    path = RAW / ("Ampicillin_auc_per_genotype.csv"
                  if drug == "Ampicillin" else "Aztreonam_auc_per_genotype.csv")
    df = pl.read_csv(path)
    df = df.rename({df.columns[0]: "_row_index"})
    return df


def load_read_counts(drug: str) -> pl.DataFrame:
    path = RAW / ("Ampicillin_read_counts_per_genotype.csv"
                  if drug == "Ampicillin" else "Aztreonam_read_counts_per_genotype.csv")
    return pl.read_csv(path)


def samples_for(md: pl.DataFrame, drug: str, conc: float, rep: int) -> list[str]:
    """Return the 4 timepoint sample names for this drug/conc/rep."""
    sub = md.filter(
        (pl.col("Drug") == drug)
        & (pl.col("conc") == conc)
        & (pl.col("Replicate") == rep)
        & (pl.col("Timepoint").is_in(TIMEPOINTS))
    )
    return sub.sort("Timepoint")["sample"].to_list()


# --- highlighted variant selection -------------------------------------------
def literal_to_masked(seq: str) -> str:
    """Convert a 13-letter residue string to the masked library encoding."""
    return "".join("." if seq[i] == WT_LITERAL[i] else seq[i] for i in range(13))


def single_mutant_genotypes() -> pl.DataFrame:
    """Return the 18 single-substitution genotypes (exactly one non-WT position)
    in the library's masked encoding: '.' elsewhere, substitution at that
    single position.
    """
    rows = []
    for i, (pos, wt, alts) in enumerate(zip(AMBLER_POS, WT_LETTER, ALT_LETTERS)):
        for a in alts:
            g = ["."] * 13
            g[i] = a
            rows.append({
                "genotype": "".join(g),
                "category": "single",
                "label": f"{wt}{pos}{a}",
                "position": pos,
                "rank": None,
            })
    return pl.DataFrame(rows)


def load_clinical_variants(auc_gts: set[str], exclude: set[str]) -> pl.DataFrame:
    """Clinical TEM alleles that are present in the library, converted to the
    masked encoding. Alleles that collide with WT or with single-mutant
    highlights are excluded so each highlight appears exactly once.
    """
    df = pl.read_csv(CLINICAL_CSV)
    df = df.with_columns([
        pl.col("Encoded_Sequence")
        .map_elements(literal_to_masked, return_dtype=pl.Utf8)
        .alias("genotype"),
    ])
    df = df.filter(pl.col("genotype").is_in(list(auc_gts)))
    df = df.filter(~pl.col("genotype").is_in(list(exclude)))
    # Collapse duplicates (e.g. TEM-12 / TEM-97 with identical encoded sequences)
    df = df.group_by("genotype", maintain_order=True).agg([
        pl.col("Variant").str.concat("/").alias("label"),
    ])
    return df.select([
        pl.col("genotype"),
        pl.lit("clinical").alias("category"),
        pl.col("label"),
        pl.lit(None, dtype=pl.Int64).alias("position"),
        pl.lit(None, dtype=pl.Int64).alias("rank"),
    ])


def top_composite_genotypes(
    auc: pl.DataFrame, drug: str, conc_high: float, n: int,
    exclude: set[str],
) -> pl.DataFrame:
    """Top-N composite variants (>= 2 non-WT positions, not already highlighted)
    ranked by mean log10(AUC) at the highest drug concentration across 3 reps.
    """
    cols = [f"{drug} {conc_high} {r}" for r in (1, 2, 3)]
    mean_expr = (
        (pl.col(cols[0]) + pl.col(cols[1]) + pl.col(cols[2])) / 3.0
    ).alias("mean_high_conc")
    sub = auc.select(["mut_profile_masked", *cols]).with_columns(mean_expr)

    def count_muts(g: str) -> int:
        if "X" in g:
            return -1  # exclude the dead genotype explicitly
        return sum(1 for ch in g if ch != ".")

    sub = sub.with_columns([
        pl.col("mut_profile_masked")
        .map_elements(count_muts, return_dtype=pl.Int64).alias("n_mut"),
    ])
    sub = sub.filter(pl.col("n_mut") >= 2)
    sub = sub.filter(~pl.col("mut_profile_masked").is_in(list(exclude)))
    sub = sub.filter(pl.col("mean_high_conc").is_finite())
    sub = sub.sort("mean_high_conc", descending=True).head(n)
    return sub.select([
        pl.col("mut_profile_masked").alias("genotype"),
        pl.lit(f"top_composite_{drug[:3].lower()}").alias("category"),
        pl.col("mut_profile_masked").alias("label"),
        pl.lit(None, dtype=pl.Int64).alias("position"),
        pl.int_range(1, sub.height + 1).alias("rank"),
    ])


def build_highlighted_table(auc_amp: pl.DataFrame, auc_azt: pl.DataFrame) -> pl.DataFrame:
    amp_gts = set(auc_amp["mut_profile_masked"].to_list())
    azt_gts = set(auc_azt["mut_profile_masked"].to_list())
    shared = amp_gts & azt_gts

    wt = pl.DataFrame({
        "genotype": [WT_GENOTYPE], "category": ["WT"], "label": ["WT (TEM-1)"],
        "position": [None], "rank": [None],
    }, schema={"genotype": pl.Utf8, "category": pl.Utf8, "label": pl.Utf8,
               "position": pl.Int64, "rank": pl.Int64})
    dead = pl.DataFrame({
        "genotype": [DEAD_GENOTYPE], "category": ["dead"],
        "label": ["TEM-1$^{dead}$"],
        "position": [None], "rank": [None],
    }, schema={"genotype": pl.Utf8, "category": pl.Utf8, "label": pl.Utf8,
               "position": pl.Int64, "rank": pl.Int64})
    singles = single_mutant_genotypes()
    pre_exclude = (
        set(wt["genotype"].to_list())
        | set(dead["genotype"].to_list())
        | set(singles["genotype"].to_list())
    )
    clinical = load_clinical_variants(shared, exclude=pre_exclude)

    already = pre_exclude | set(clinical["genotype"].to_list())

    top_amp = top_composite_genotypes(auc_amp, "Ampicillin", 781.0,
                                       N_TOP_COMPOSITE, already)
    already |= set(top_amp["genotype"].to_list())
    top_azt = top_composite_genotypes(auc_azt, "Aztreonam", 324.0,
                                       N_TOP_COMPOSITE, already)

    return pl.concat([wt, dead, singles, clinical, top_amp, top_azt])


# --- per-variant fitness and low-count flags ---------------------------------
def compute_variant_fitness(
    auc: pl.DataFrame, reads: pl.DataFrame, md: pl.DataFrame,
    drug: str, highlighted: list[str],
) -> pl.DataFrame:
    """Per highlighted genotype x concentration: mean, SD across reps + per-rep
    sum_reads + low_count flags.
    Only genotypes present in both tables are returned.
    """
    auc_sub = auc.filter(pl.col("mut_profile_masked").is_in(highlighted))
    reads_sub = reads.filter(pl.col("mut_profile_masked").is_in(highlighted))

    # Join on genotype
    auc_map = {row["mut_profile_masked"]: row for row in auc_sub.iter_rows(named=True)}
    reads_map = {row["mut_profile_masked"]: row for row in reads_sub.iter_rows(named=True)}

    rows: list[dict] = []
    concs = DRUG_CONCS[drug]
    for gt in highlighted:
        if gt not in auc_map or gt not in reads_map:
            continue
        auc_row = auc_map[gt]
        reads_row = reads_map[gt]
        for c in concs:
            auc_cols = [f"{drug} {c} {r}" for r in (1, 2, 3)]
            aucs = np.array([auc_row.get(cc) for cc in auc_cols], dtype=np.float64)
            sums = np.zeros(3, dtype=np.int64)
            for r in (1, 2, 3):
                s_cols = samples_for(md, drug, c, r)
                sums[r - 1] = int(np.sum([
                    reads_row.get(s, 0) if reads_row.get(s, 0) is not None else 0
                    for s in s_cols
                ]))
            low_flag = (sums < LOW_COUNT_THRESHOLD)  # per-rep
            n_low = int(low_flag.sum())

            # valid AUCs: replicates with sums >= threshold and finite AUC
            valid = (~low_flag) & np.isfinite(aucs)
            if valid.any():
                mean_auc = float(np.mean(aucs[valid]))
                sd_auc = float(np.std(aucs[valid], ddof=1)) if valid.sum() > 1 else float("nan")
            else:
                mean_auc = float("nan")
                sd_auc = float("nan")

            mean_raw = float(np.nanmean(aucs))
            sd_raw = float(np.nanstd(aucs, ddof=1)) if np.isfinite(aucs).sum() > 1 else float("nan")

            rows.append({
                "drug": drug,
                "genotype": gt,
                "concentration": c,
                "rep1_auc": float(aucs[0]) if np.isfinite(aucs[0]) else None,
                "rep2_auc": float(aucs[1]) if np.isfinite(aucs[1]) else None,
                "rep3_auc": float(aucs[2]) if np.isfinite(aucs[2]) else None,
                "rep1_sum_reads": int(sums[0]),
                "rep2_sum_reads": int(sums[1]),
                "rep3_sum_reads": int(sums[2]),
                "rep1_low_count": bool(low_flag[0]),
                "rep2_low_count": bool(low_flag[1]),
                "rep3_low_count": bool(low_flag[2]),
                "n_low_count_reps": n_low,
                "all_reps_low": bool(n_low == 3),
                # Flag as "majority low-count" when 2 of 3 reps have insufficient
                # reads; mean is then driven by a single surviving rep and is
                # visually unreliable.
                "majority_low": bool(n_low >= 2),
                "mean_raw": mean_raw,
                "sd_raw": sd_raw,
                "mean_clean": mean_auc,  # excludes flagged reps
                "sd_clean": sd_auc,
            })
    return pl.DataFrame(rows)


# --- whole-library low-count census (for results_table + LOD) ----------------
def library_census(
    reads: pl.DataFrame, md: pl.DataFrame, drug: str,
) -> pl.DataFrame:
    """Per drug x concentration: fraction of (genotype, rep) with sum_reads < 10;
    total fraction all-reps-low (collapses genotype to a single-conc point
    that would be flagged on the plot)."""
    rows: list[dict] = []
    concs = DRUG_CONCS[drug]
    for c in concs:
        flag = np.zeros((reads.height, 3), dtype=bool)
        for r in (1, 2, 3):
            s_cols = samples_for(md, drug, c, r)
            arr = reads.select([
                pl.sum_horizontal(pl.col(cc) for cc in s_cols).alias("s"),
            ])["s"].to_numpy()
            flag[:, r - 1] = arr < LOW_COUNT_THRESHOLD
        rows.append({
            "drug": drug,
            "concentration": c,
            "n_genotypes": reads.height,
            "frac_rep_low_count": float(flag.mean()),
            "frac_all_reps_low": float(np.all(flag, axis=1).mean()),
            "frac_any_rep_low": float(np.any(flag, axis=1).mean()),
            "frac_majority_low": float((flag.sum(axis=1) >= 2).mean()),
        })
    return pl.DataFrame(rows)


# --- LOD estimate -------------------------------------------------------------
def compute_lod(variant_df: pl.DataFrame, drug: str) -> dict[float, float]:
    """Per-concentration LOD = median fitness of TEM-1_dead across reps (the
    non-functional control sets the empirical AUC floor).
    """
    lod: dict[float, float] = {}
    for c in DRUG_CONCS[drug]:
        row = variant_df.filter(
            (pl.col("drug") == drug)
            & (pl.col("genotype") == DEAD_GENOTYPE)
            & (pl.col("concentration") == c)
        )
        if row.height:
            v = row["mean_raw"][0]
            lod[c] = float(v) if v is not None else float("nan")
        else:
            lod[c] = float("nan")
    return lod


# --- plotting ----------------------------------------------------------------
POSITION_CMAP = {  # one color per Ambler position for single-mutant highlights
    21: "#8e44ad", 39: "#2980b9", 69: "#16a085", 104: "#27ae60",
    164: "#f39c12", 182: "#d35400", 237: "#c0392b", 238: "#7f8c8d",
    240: "#2c3e50", 244: "#e67e22", 265: "#1abc9c", 275: "#34495e",
    276: "#c9a227",
}
CATEGORY_STYLE = {
    "WT":       dict(color="#000000", linewidth=2.4, linestyle="-", zorder=12),
    "dead":     dict(color="#3a4a7a", linewidth=2.2, linestyle="--", zorder=11),
    "single":   dict(linewidth=1.6, linestyle="-", zorder=9, alpha=0.90),  # color set below
    "clinical": dict(color="#f5b700", linewidth=0.7, linestyle="-", zorder=5, alpha=0.55),
    "top_composite_amp": dict(color="#b02a2a", linewidth=0.7, linestyle="-", zorder=4, alpha=0.45),
    "top_composite_azt": dict(color="#1f6fb5", linewidth=0.7, linestyle="-", zorder=4, alpha=0.45),
}
LOW_COUNT_ALPHA = 0.25


def plot_variant(ax, concs, means, sds, low_flags, style: dict,
                 show_sd_band: bool = True) -> None:
    """Plot a single variant's dose-response trace.

    concs:      list of concentrations (x-axis; already offset from 0 for log).
    means:      mean fitness at each conc (np.nan where undefined).
    sds:        SD across clean replicates.
    low_flags:  list[bool] indicating "all reps flagged" at this conc.
    """
    concs = np.asarray(concs)
    means = np.asarray(means)
    sds = np.asarray(sds)
    low_flags = np.asarray(low_flags)

    # Break the line at fully-flagged points so the "noisy jump" is visually
    # disconnected from the reliable trend.
    x_draw = np.where(low_flags, np.nan, concs.astype(float))
    y_draw = np.where(low_flags, np.nan, means)

    kw = {k: v for k, v in style.items() if k in {"color", "linewidth",
                                                     "linestyle", "zorder",
                                                     "alpha"}}
    ax.plot(x_draw, y_draw, marker=None, **kw)

    if show_sd_band:
        lo = means - sds
        hi = means + sds
        # Same masking rule
        lo = np.where(low_flags, np.nan, lo)
        hi = np.where(low_flags, np.nan, hi)
        color = style.get("color", "#555555")
        alpha_band = min(0.18, style.get("alpha", 0.8) * 0.25)
        ax.fill_between(concs, lo, hi, color=color, alpha=alpha_band,
                         linewidth=0, zorder=style.get("zorder", 5) - 1)

    # Per-point markers; open marker if low-count-flagged, filled otherwise
    for x, y, flag in zip(concs, means, low_flags):
        if not np.isfinite(y):
            continue
        if flag:
            ax.plot(
                x, y, marker="o", markersize=3.6,
                markerfacecolor="none", markeredgecolor=style.get("color", "#555555"),
                markeredgewidth=0.9, alpha=LOW_COUNT_ALPHA,
                zorder=style.get("zorder", 5) - 1,
            )
        else:
            ax.plot(
                x, y, marker="o", markersize=3.6,
                color=style.get("color", "#555555"),
                alpha=min(1.0, style.get("alpha", 0.9) + 0.1),
                zorder=style.get("zorder", 5) + 1,
            )


def make_panel(ax, drug: str, per_variant: pl.DataFrame,
               highlighted: pl.DataFrame, lod: dict[float, float]) -> None:
    concs = DRUG_CONCS[drug]
    # For log-x: substitute 0 -> an artificial small tick just below lowest nonzero conc
    nonzero = [c for c in concs if c > 0]
    # Place the 0 point at half the lowest nonzero concentration
    x0_proxy = min(nonzero) / 3.0
    x_vals = [x0_proxy if c == 0 else c for c in concs]

    # LOD band: shade the region at-or-below the TEM-1_dead curve per
    # concentration. Fitness values in this band are indistinguishable from a
    # catalytically dead enzyme and mark the assay's limit of detection.
    lod_vals = np.array([lod.get(c, np.nan) for c in concs], dtype=np.float64)
    if np.any(np.isfinite(lod_vals)):
        # Use a floor one fitness-unit below the dead curve so the patch is
        # visible even when dead fitness is near zero.
        bottom = float(np.nanmin(lod_vals)) - 1.5
        ax.fill_between(
            x_vals, bottom, lod_vals,
            where=np.isfinite(lod_vals),
            color="#d0d0d0", alpha=0.40, linewidth=0, zorder=1,
            step="mid",
        )

    cat_to_style = {c: s.copy() for c, s in CATEGORY_STYLE.items()}
    drug_composite_cat = f"top_composite_{drug[:3].lower()}"

    # Build per-genotype row map for rapid plotting
    hl_rows = {r["genotype"]: r for r in highlighted.iter_rows(named=True)}

    # Stable draw order: background first, foreground last. Composites and
    # clinical lines are drawn translucent in the back so the single-mutant
    # family and the anchor traces (WT, dead) remain legible.
    draw_order = [drug_composite_cat, "clinical", "single", "dead", "WT"]
    for cat in draw_order:
        if cat == drug_composite_cat:
            subs = highlighted.filter(pl.col("category") == drug_composite_cat)
        else:
            subs = highlighted.filter(pl.col("category") == cat)

        for row in subs.iter_rows(named=True):
            gt = row["genotype"]
            vd = per_variant.filter(
                (pl.col("drug") == drug) & (pl.col("genotype") == gt)
            ).sort("concentration")
            if vd.is_empty():
                continue
            means = vd["mean_clean"].to_numpy()
            sds = vd["sd_clean"].to_numpy()
            # "Low-count" visualisation: either all 3 reps flagged (no valid
            # mean) OR a majority (>=2 of 3) reps flagged so the surviving
            # value drives the apparent mean.
            low_flags = vd["majority_low"].to_numpy()

            style = cat_to_style[cat].copy()
            if cat == "single":
                pos = row["position"]
                style["color"] = POSITION_CMAP.get(pos, "#888888")

            # For single / clinical / composite: hide SD band (too busy).
            show_sd = cat in ("WT", "dead")

            plot_variant(ax, x_vals, means, sds, low_flags, style,
                          show_sd_band=show_sd)

    ax.set_xscale("log")
    ax.set_xlabel(f"{drug} (µg/mL)", fontsize=10)
    ax.set_ylabel("Fitness (log$_{10}$ AUC)", fontsize=10)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # x tick labels: show 0 at the proxy position
    ax.set_xticks(x_vals)
    ax.set_xticklabels([f"{c:g}" for c in concs], fontsize=8, rotation=0)
    ax.set_xlim(x_vals[0] / 1.3, x_vals[-1] * 1.3)

    # "0 axis break" mark: dashed vertical gap between proxy and lowest nonzero
    ax.axvline(x_vals[0] * 1.7, color="#e0e0e0", linewidth=0.7,
                linestyle=":", zorder=0)


def build_legend_handles(drug: str) -> list:
    handles = [
        Line2D([0], [0], color="#000000", linewidth=2.2, linestyle="-",
                label="WT (TEM-1)"),
        Line2D([0], [0], color="#3a4a7a", linewidth=2.0, linestyle="--",
                label="TEM-1$^{dead}$ (LOD)"),
        Line2D([0], [0], color="#f5b700", linewidth=1.3,
                label="Clinical TEM alleles"),
        Line2D([0], [0], color="#666666", linewidth=1.3,
                label="Single substitutions (18)"),
    ]
    composite_color = "#b02a2a" if drug == "Ampicillin" else "#1f6fb5"
    handles.append(Line2D([0], [0], color=composite_color, linewidth=1.3,
                            label=f"Top-{N_TOP_COMPOSITE} composites at max conc."))
    handles.append(Line2D([0], [0], color="#888888", linewidth=0, marker="o",
                            markersize=5, markerfacecolor="none",
                            markeredgecolor="#888888",
                            label=(f"Low-count (>=2 of 3 reps < "
                                   f"{LOW_COUNT_THRESHOLD} reads)")))
    handles.append(Patch(facecolor="#d0d0d0", alpha=0.5,
                          label="LOD band (TEM-1$^{dead}$)"))
    return handles


def plot_combined(per_variant: pl.DataFrame, highlighted: pl.DataFrame,
                   lods: dict[str, dict[float, float]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 7.0), dpi=150)
    make_panel(axes[0], "Ampicillin", per_variant, highlighted, lods["Ampicillin"])
    make_panel(axes[1], "Aztreonam", per_variant, highlighted, lods["Aztreonam"])

    axes[0].set_title("A. Ampicillin", loc="left", fontsize=13, fontweight="bold")
    axes[1].set_title("B. Aztreonam", loc="left", fontsize=13, fontweight="bold")

    # Unified y-axis range across both drugs
    ylo, yhi = +np.inf, -np.inf
    for ax in axes:
        a, b = ax.get_ylim()
        ylo, yhi = min(ylo, a), max(yhi, b)
    for ax in axes:
        ax.set_ylim(ylo, yhi)

    handles = build_legend_handles("Ampicillin")
    fig.legend(
        handles=handles,
        loc="lower center", bbox_to_anchor=(0.5, 0.005),
        ncol=7, fontsize=8.5, frameon=False, handletextpad=0.5,
        columnspacing=1.1,
    )

    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    out = FIGDIR / "fig_s3_revised.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def plot_single_panel(drug: str, per_variant: pl.DataFrame,
                       highlighted: pl.DataFrame,
                       lod: dict[float, float]) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 5.6), dpi=150)
    make_panel(ax, drug, per_variant, highlighted, lod)
    handles = build_legend_handles(drug)
    ax.legend(handles=handles, loc="lower left",
               fontsize=7.5, frameon=False, ncol=1, handletextpad=0.6)
    tag = drug[:3].lower()
    fig.tight_layout()
    out = SUPP / f"figure_s17_{tag}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# --- main --------------------------------------------------------------------
def main():
    np.random.seed(SEED)
    md = load_metadata()
    auc_amp = load_auc("Ampicillin")
    auc_azt = load_auc("Aztreonam")
    rc_amp = load_read_counts("Ampicillin")
    rc_azt = load_read_counts("Aztreonam")

    print("Building highlighted-variant table...")
    highlighted = build_highlighted_table(auc_amp, auc_azt)
    highlighted.write_csv(DATADIR / "highlighted_variants.csv")
    print(f"  {highlighted.height} highlighted variants")
    print(highlighted.group_by("category").len().sort("category"))

    highlighted_gts = highlighted["genotype"].unique().to_list()

    print("Computing per-variant dose-response stats...")
    pv_amp = compute_variant_fitness(auc_amp, rc_amp, md, "Ampicillin",
                                       highlighted_gts)
    pv_azt = compute_variant_fitness(auc_azt, rc_azt, md, "Aztreonam",
                                       highlighted_gts)
    per_variant = pl.concat([pv_amp, pv_azt])
    per_variant.write_csv(DATADIR / "per_variant_fitness.csv")
    print(f"  wrote {per_variant.height:,} rows")

    print("Running library-wide low-count census...")
    cen_amp = library_census(rc_amp, md, "Ampicillin")
    cen_azt = library_census(rc_azt, md, "Aztreonam")
    census = pl.concat([cen_amp, cen_azt])
    census.write_csv(HERE / "results_table.csv")
    print(f"  wrote {census.height} rows to results_table.csv")

    print("Computing TEM-1_dead LOD per conc...")
    lods = {
        "Ampicillin": compute_lod(per_variant, "Ampicillin"),
        "Aztreonam": compute_lod(per_variant, "Aztreonam"),
    }
    for drug, d in lods.items():
        print(f"  {drug} LOD:")
        for c, v in d.items():
            print(f"    conc={c:>7g}  LOD={v:.3f}")

    print("Drawing figures...")
    plot_single_panel("Ampicillin", per_variant, highlighted, lods["Ampicillin"])
    plot_single_panel("Aztreonam", per_variant, highlighted, lods["Aztreonam"])
    plot_combined(per_variant, highlighted, lods)

    print("\n=== Library-wide extinction-artifact census ===")
    for row in census.iter_rows(named=True):
        print(
            f"  {row['drug']:11s}  conc={row['concentration']:>7g}  "
            f"frac_any_rep<10r={row['frac_any_rep_low']*100:5.1f}%  "
            f"frac_maj_rep<10r={row['frac_majority_low']*100:5.1f}%  "
            f"frac_all_reps<10r={row['frac_all_reps_low']*100:5.1f}%"
        )

    # Highlighted-variant flag rates (the numbers that matter for the figure)
    all_low = per_variant.filter(pl.col("all_reps_low"))
    maj_low = per_variant.filter(pl.col("majority_low"))
    hl_total = per_variant.height
    print(
        f"\n=== Highlighted variants ({highlighted.height}): "
        f"all-reps-low = {all_low.height}/{hl_total} ({all_low.height/hl_total*100:.2f}%); "
        f"majority-low = {maj_low.height}/{hl_total} ({maj_low.height/hl_total*100:.2f}%) ==="
    )
    by_drug_maj = maj_low.group_by("drug").len().sort("drug")
    print("majority-low by drug:")
    print(by_drug_maj)


if __name__ == "__main__":
    main()
