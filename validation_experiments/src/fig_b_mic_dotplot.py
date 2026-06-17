"""Figure B: MIC Dot Plot with Error Bars — faceted by drug."""

import pandas as pd
import polars as pl
import altair as alt
from scipy.stats import t as t_dist

from config import (
    load_processed, register_theme, FIGURES_DIR, experiments,
    FONT_MONO, SHORT_COLORS, VARIANT_CLUSTERS,
    make_cluster_legend, extract_short_name, VARIANT_ORDER,
)

# ── Variant x-positions: integer coords with +1 gap between clusters ──
_CLUSTER_LISTS = list(VARIANT_CLUSTERS.values())


def _build_xpos() -> dict[str, int]:
    """Map short variant name → integer x position, +1 gap between clusters."""
    pos = {}
    x = 0
    for i, members in enumerate(_CLUSTER_LISTS):
        if i > 0:
            x += 1  # gap
        for m in members:
            pos[m] = x
            x += 1
    return pos


XPOS = _build_xpos()

# Cluster boundary x-values (midpoint of gap)
_CLUSTER_BOUNDARIES: list[float] = []
_cursor = 0
for i, members in enumerate(_CLUSTER_LISTS):
    if i > 0:
        _CLUSTER_BOUNDARIES.append(_cursor - 0.5)
    _cursor += len(members)
    if i < len(_CLUSTER_LISTS) - 1:
        _cursor += 1  # account for gap


def _drug_panel(
    drug: str,
    summary_pd: pd.DataFrame,
    rep_pd: pd.DataFrame,
    y_col: str = "log10_mic",
    y_title: str = "MIC (µg/mL)",
) -> alt.LayerChart:
    """Single drug panel with quantitative x-axis, cluster gaps, variant whiskers."""
    drug_sum = summary_pd[summary_pd["drug"] == drug].copy()
    drug_rep = rep_pd[rep_pd["drug"] == drug].copy()

    # Map labels to short names and x positions
    drug_sum["short"] = drug_sum["label"].apply(extract_short_name)
    drug_sum["xpos"] = drug_sum["short"].map(XPOS)
    drug_rep["short"] = drug_rep["label"].apply(extract_short_name)
    drug_rep["xpos"] = drug_rep["short"].map(XPOS)

    # Variant color by short name
    color_domain = list(SHORT_COLORS.keys())
    color_range = [SHORT_COLORS[k] for k in color_domain]
    variant_color = alt.Color(
        "short:N",
        scale=alt.Scale(domain=color_domain, range=color_range),
        legend=None,
    )

    # Build labelExpr: map xpos → short name
    label_parts = [f"datum.value === {v} ? '{k}'" for k, v in XPOS.items()]
    label_expr = " : ".join(label_parts) + " : ''"

    x_enc = alt.X(
        "xpos:Q",
        title=None,
        scale=alt.Scale(domain=[-0.8, max(XPOS.values()) + 0.8]),
        axis=alt.Axis(
            values=list(XPOS.values()),
            labelExpr=label_expr,
            labelFont=FONT_MONO,
            labelFontSize=10,
            labelAngle=0,
            grid=False,
            tickSize=0,
        ),
    )

    y_axis_cfg = alt.Axis(
        labelExpr=(
            "pow(10, datum.value) >= 1000"
            " ? format(pow(10, datum.value) / 1000, '.1f') + 'k'"
            " : format(pow(10, datum.value), '.1f')"
        ),
    )

    # Cluster separator rules
    sep_df = pd.DataFrame({"bx": _CLUSTER_BOUNDARIES})
    separators = (
        alt.Chart(sep_df)
        .mark_rule(strokeDash=[3, 3], strokeWidth=0.8, color="#BDBDBD")
        .encode(x=alt.X("bx:Q", scale=alt.Scale(domain=[-0.8, max(XPOS.values()) + 0.8])))
    )

    whiskers = (
        alt.Chart(drug_sum)
        .mark_rule(strokeWidth=1.5)
        .encode(
            x=x_enc,
            y=alt.Y("ci_lo:Q", title=y_title, axis=y_axis_cfg),
            y2="ci_hi:Q",
            color=variant_color,
        )
    )

    means = (
        alt.Chart(drug_sum)
        .mark_point(size=80, filled=True)
        .encode(
            x=x_enc,
            y=alt.Y("mean_log10:Q", title=y_title, axis=y_axis_cfg),
            color=variant_color,
        )
    )

    # Split reps into censored (open circles) and uncensored (filled)
    has_censored = "is_censored" in drug_rep.columns
    rep_uncensored = drug_rep[~drug_rep["is_censored"]] if has_censored else drug_rep
    rep_censored = drug_rep[drug_rep["is_censored"]] if has_censored else drug_rep.iloc[:0]

    reps_filled = (
        alt.Chart(rep_uncensored)
        .mark_point(size=20, filled=True, opacity=0.35)
        .encode(
            x=x_enc,
            y=alt.Y(f"{y_col}:Q", title=y_title, axis=y_axis_cfg),
            color=variant_color,
        )
    )

    reps_open = (
        alt.Chart(rep_censored)
        .mark_point(size=20, filled=False, strokeWidth=1.2, opacity=0.5)
        .encode(
            x=x_enc,
            y=alt.Y(f"{y_col}:Q", title=y_title, axis=y_axis_cfg),
            color=variant_color,
        )
    )

    return (separators + whiskers + reps_filled + reps_open + means).properties(
        width=320, height=280, title=drug,
    )


def _build_dotplot(
    rep_df: pl.DataFrame,
    y_col: str,
    y_title: str,
    title: str,
    stem: str,
    suffix: str,
):
    """Build and save a MIC dotplot from a rep-level dataframe."""
    mic_summary = (
        rep_df
        .group_by("drug", "label")
        .agg(
            pl.col(y_col).mean().alias("mean_log10"),
            pl.col(y_col).std().alias("std_log10"),
            pl.col(y_col).count().alias("n"),
        )
        .with_columns(
            (pl.col("std_log10") / pl.col("n").sqrt()).alias("sem_log10"),
        )
    )

    summary_pd = mic_summary.to_pandas()
    summary_pd["ci95"] = summary_pd.apply(
        lambda r: (
            t_dist.ppf(0.975, r["n"] - 1) * r["sem_log10"]
            if r["n"] > 1 and r["sem_log10"] > 0
            else 0
        ),
        axis=1,
    )
    summary_pd["ci_lo"] = summary_pd["mean_log10"] - summary_pd["ci95"]
    summary_pd["ci_hi"] = summary_pd["mean_log10"] + summary_pd["ci95"]

    rep_pd = rep_df.to_pandas()

    panels = [
        _drug_panel(drug, summary_pd, rep_pd, y_col=y_col, y_title=y_title)
        for drug in experiments
    ]
    main = (
        alt.hconcat(*panels)
        .resolve_scale(y="shared")
        .properties(title=title)
    )

    fig = (
        alt.vconcat(main, make_cluster_legend())
        .properties(padding={"bottom": 10, "right": 10})
    )

    full_stem = f"{stem}{suffix}"
    fig.save(str(FIGURES_DIR / f"{full_stem}.png"), scale_factor=2)
    fig.save(str(FIGURES_DIR / f"{full_stem}.html"))
    print(f"Saved {full_stem}")
    return fig


def make_figure(suffix: str = ""):
    register_theme()

    # D1: MIC from AUC
    mic_rep_df = load_processed("mic_rep_df", suffix)
    _build_dotplot(
        mic_rep_df, y_col="log10_mic", y_title="AUC-based MIC (µg/mL)",
        title="AUC-based MIC by Variant (geometric mean ± 95% CI)",
        stem="mic_dotplot", suffix=suffix,
    )

    # D2: MIC from OD
    od_mic_rep_df = load_processed("od_mic_rep_df", suffix)
    _build_dotplot(
        od_mic_rep_df, y_col="log10_mic", y_title="FinalOD-based MIC (µg/mL)",
        title="FinalOD-based MIC by Variant (geometric mean ± 95% CI)",
        stem="od_mic_dotplot", suffix=suffix,
    )


if __name__ == "__main__":
    import sys
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    s = ""
    if "--suffix" in sys.argv:
        s = sys.argv[sys.argv.index("--suffix") + 1]
    make_figure(suffix=s)
