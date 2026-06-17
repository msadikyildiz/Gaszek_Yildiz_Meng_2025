"""Figure B1: AUC vs OD Metric Correlation — 2×2 grid (MIC/IC50 × AZT/AMP)."""

import numpy as np
import pandas as pd
import polars as pl
import altair as alt
from scipy import stats

from config import (
    load_processed, register_theme, FIGURES_DIR, experiments,
    VARIANT_ORDER, VARIANT_COLORS, FONT_BODY, FONT_MONO,
    make_cluster_legend, extract_short_name,
    nudge_labels, make_connectors,
)


def _sig_stars(p: float) -> str:
    if p < 0.001:
        return " ***"
    if p < 0.01:
        return " **"
    if p < 0.05:
        return " *"
    return " (n.s.)"


def _scatter(
    sub_pd: pd.DataFrame,
    x_col: str, y_col: str,
    x_label: str, y_label: str,
    color_enc: alt.Color,
    drug: str,
) -> alt.LayerChart:
    """Scatter: AUC-based vs OD-based metric with identity + regression."""
    x_arr = sub_pd[x_col].to_numpy().astype(float)
    y_arr = sub_pd[y_col].to_numpy().astype(float)
    sub_pd = sub_pd.copy()

    x_range = x_arr.max() - x_arr.min() or 1.0
    y_range = y_arr.max() - y_arr.min() or 1.0
    x_pad = max(x_range * 0.08, 0.05)
    y_pad = max(y_range * 0.08, 0.05)

    # Use shared domain for identity line
    all_vals = np.concatenate([x_arr, y_arr])
    lo = float(all_vals.min()) - max(x_pad, y_pad)
    hi = float(all_vals.max()) + max(x_pad, y_pad)
    x_scale = alt.Scale(domain=[lo, hi])
    y_scale = alt.Scale(domain=[lo, hi])

    r, p = stats.pearsonr(x_arr, y_arr)
    significant = p < 0.05

    # Nudge labels
    lx, ly = nudge_labels(
        x_arr.copy(), y_arr.copy(), x_range, y_range,
        all_pts_x=x_arr, all_pts_y=y_arr,
    )
    sub_pd["_lx"] = lx
    sub_pd["_ly"] = ly
    sub_pd["short_name"] = sub_pd["label"].apply(extract_short_name)

    # Points
    points = (
        alt.Chart(sub_pd)
        .mark_point(size=120, filled=True, opacity=0.85,
                    stroke="#FAFAFA", strokeWidth=1.5)
        .encode(
            x=alt.X(f"{x_col}:Q", title=x_label, scale=x_scale),
            y=alt.Y(f"{y_col}:Q", title=y_label, scale=y_scale),
            color=color_enc,
        )
    )

    labels = (
        alt.Chart(sub_pd)
        .mark_text(fontSize=8, font=FONT_MONO, align="left", color="#424242")
        .encode(
            x=alt.X("_lx:Q", scale=x_scale),
            y=alt.Y("_ly:Q", scale=y_scale),
            text="short_name:N",
        )
    )

    conn = make_connectors(
        sub_pd, x_col, y_col,
        x_range=x_range, y_range=y_range,
    )

    # Identity line
    id_df = pd.DataFrame({x_col: [lo, hi], y_col: [lo, hi]})
    id_line = (
        alt.Chart(id_df)
        .mark_line(strokeDash=[2, 2], strokeWidth=0.8, color="#BDBDBD")
        .encode(x=f"{x_col}:Q", y=f"{y_col}:Q")
    )

    # Regression line
    slope, intercept = np.polyfit(x_arr, y_arr, 1)
    x_line = np.array([x_arr.min(), x_arr.max()])
    y_line = slope * x_line + intercept
    reg_df = pd.DataFrame({x_col: x_line, y_col: y_line})
    reg = (
        alt.Chart(reg_df)
        .mark_line(
            strokeDash=[6, 3],
            strokeWidth=1.5 if significant else 1.0,
            opacity=0.5 if significant else 0.2,
            color="#424242" if significant else "#BDBDBD",
        )
        .encode(x=f"{x_col}:Q", y=f"{y_col}:Q")
    )

    # r,p annotation
    stars = _sig_stars(p)
    ann_color = "#1A1A1A" if significant else "#9E9E9E"
    ann_text = f"r={r:.2f}, p={p:.3f}{stars}"
    rp_data = pd.DataFrame({
        "ax": [hi - (hi - lo) * 0.02],
        "ay": [lo + (hi - lo) * 0.05],
        "t": [ann_text],
    })
    rp_ann = (
        alt.Chart(rp_data)
        .mark_text(align="right", baseline="bottom",
                   fontSize=10, font=FONT_MONO,
                   fontWeight="bold" if significant else "normal",
                   color=ann_color)
        .encode(x=alt.X("ax:Q", scale=x_scale),
                y=alt.Y("ay:Q", scale=y_scale), text="t:N")
    )

    layers = [id_line, reg]
    if conn is not None:
        layers.append(conn)
    layers += [points, labels, rp_ann]

    return alt.layer(*layers).properties(
        width=280, height=240, title=f"{drug}",
    )


def make_figure(suffix: str = ""):
    register_theme()
    xref = load_processed("xref_expanded_df", suffix)

    color_domain = [v for v in VARIANT_ORDER if v != "DD"]
    color_range = [VARIANT_COLORS[v] for v in color_domain]
    color_enc = alt.Color(
        "label:N",
        scale=alt.Scale(domain=color_domain, range=color_range),
        legend=None,
    )

    metric_pairs = [
        ("log10_mic", "log10_od_mic", "log₁₀(MIC_AUC)", "log₁₀(MIC_FinalOD)"),
        ("log10_ic50", "log10_od_ic50", "log₁₀(IC₅₀_AUC)", "log₁₀(IC₅₀_FinalOD)"),
    ]

    scatter_rows = []
    for x_col, y_col, x_label, y_label in metric_pairs:
        row_panels = []
        for drug in experiments:
            drug_data = xref.filter(
                (pl.col("drug") == drug) & pl.col("genotype_13").is_not_null(),
            )
            cols = ["label", x_col, y_col]
            sub = (
                drug_data.select([c for c in cols if c in drug_data.columns])
                .drop_nulls()
            )
            if len(sub) < 3:
                continue
            panel = _scatter(
                sub.to_pandas(), x_col, y_col, x_label, y_label,
                color_enc, drug,
            )
            row_panels.append(panel)
        if row_panels:
            scatter_rows.append(alt.hconcat(*row_panels))

    if not scatter_rows:
        print("WARNING: No data for metric correlation figure")
        return None

    fig = (
        alt.vconcat(*scatter_rows, make_cluster_legend(exclude={"DD"}))
        .resolve_scale(color="shared")
        .properties(
            title="AUC vs OD Metric Correlation",
            padding={"bottom": 10, "right": 10},
        )
    )

    stem = f"mic_metric_correlation{suffix}"
    fig.save(str(FIGURES_DIR / f"{stem}.png"), scale_factor=2)
    fig.save(str(FIGURES_DIR / f"{stem}.html"))
    print(f"Saved {stem}")
    return fig


if __name__ == "__main__":
    import sys
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    s = ""
    if "--suffix" in sys.argv:
        s = sys.argv[sys.argv.index("--suffix") + 1]
    make_figure(suffix=s)
