"""Figure C: MIC/IC50 vs Epistasis Fitness — full scatter grid."""

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
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _vl_field(name: str) -> str:
    """Escape field name for Vega-Lite path parsing (e.g., fitness_0.0)."""
    return name.replace(".", r"\.")


def _scatter_panel(
    sub_pd: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    color_enc: alt.Color,
) -> alt.LayerChart:
    """Single scatter panel with Pearson r/p subtitle, regression, connectors."""
    x_arr = sub_pd[x_col].to_numpy().astype(float)
    y_arr = sub_pd[y_col].to_numpy().astype(float)
    sub_pd = sub_pd.copy()

    # Explicit axis domains with 8% padding
    x_range = x_arr.max() - x_arr.min() or 1.0
    y_range = y_arr.max() - y_arr.min() or 1.0
    x_pad = max(x_range * 0.08, 0.05)
    y_pad = max(y_range * 0.08, 0.3)
    x_domain = [float(x_arr.min() - x_pad), float(x_arr.max() + x_pad)]
    y_domain = [float(y_arr.min() - y_pad), float(y_arr.max() + y_pad)]

    # Shared label nudging (label-label + label-point repulsion)
    lx, ly = nudge_labels(
        x_arr.copy(), y_arr.copy(), x_range, y_range,
        all_pts_x=x_arr, all_pts_y=y_arr,
    )
    sub_pd["_lx"] = lx
    sub_pd["_ly"] = ly

    r, p = stats.pearsonr(x_arr, y_arr)
    stars = _sig_stars(p)
    significant = p < 0.05

    slope, intercept = np.polyfit(x_arr, y_arr, 1)
    x_line = np.array([x_arr.min(), x_arr.max()])
    y_line = slope * x_line + intercept
    reg_df = pd.DataFrame({x_col: x_line, y_col: y_line})

    sub_pd["_short"] = sub_pd["label"].apply(extract_short_name)
    x_field = _vl_field(x_col)
    y_field = _vl_field(y_col)

    x_scale = alt.Scale(domain=x_domain)
    y_scale = alt.Scale(domain=y_domain)

    points = (
        alt.Chart(sub_pd)
        .mark_point(
            size=90, filled=True, opacity=0.85,
            stroke="#FAFAFA", strokeWidth=1.5,
        )
        .encode(
            x=alt.X(f"{x_field}:Q", title=x_label, scale=x_scale),
            y=alt.Y(f"{y_field}:Q", title=y_label, scale=y_scale),
            color=color_enc,
        )
    )

    pt_labels = (
        alt.Chart(sub_pd)
        .mark_text(fontSize=8, font=FONT_MONO, align="left", color="#424242")
        .encode(
            x=alt.X("_lx:Q", scale=x_scale),
            y=alt.Y("_ly:Q", scale=y_scale),
            text="_short:N",
        )
    )

    conn = make_connectors(
        sub_pd, x_col, y_col,
        x_range=x_range, y_range=y_range,
    )

    reg_line = (
        alt.Chart(reg_df)
        .mark_line(
            strokeDash=[6, 3],
            strokeWidth=2.0 if significant else 1.0,
            opacity=0.6 if significant else 0.2,
            color="#424242" if significant else "#BDBDBD",
        )
        .encode(
            x=alt.X(f"{x_field}:Q", scale=x_scale),
            y=alt.Y(f"{y_field}:Q", scale=y_scale),
        )
    )

    # r,p as panel subtitle (non-overlapping, above plot)
    sub_color = "#1A1A1A" if significant else "#9E9E9E"
    ann_text = f"r={r:.2f}, p={p:.3f}" + (f" {stars}" if stars else "")
    panel_title = alt.Title(
        text=x_label,
        subtitle=ann_text,
        subtitleColor=sub_color,
        subtitleFont=FONT_MONO,
        subtitleFontSize=9,
    )

    layers = [reg_line]
    if conn is not None:
        layers.append(conn)
    layers += [points, pt_labels]

    return alt.layer(*layers).properties(
        width=180, height=160, title=panel_title,
    )


def _build_grid(drug: str, xref_expanded_df: pl.DataFrame) -> alt.VConcatChart:
    """Build 3-row × N-col scatter grid for one drug."""
    color_domain = [v for v in VARIANT_ORDER if v != "DD"]
    color_range = [VARIANT_COLORS[v] for v in color_domain]
    color_enc = alt.Color(
        "label:N",
        scale=alt.Scale(domain=color_domain, range=color_range),
        legend=None,
    )

    drug_data = xref_expanded_df.filter(
        (pl.col("drug") == drug) & pl.col("genotype_13").is_not_null(),
    )

    fitness_cols = sorted(
        [
            c for c in drug_data.columns
            if c.startswith("fitness_")
            and drug_data[c].drop_nulls().len() > 0
        ],
        key=lambda s: float(s.replace("fitness_", "")),
    )
    x_cols = fitness_cols + ["mean_fitness", "PC1"]
    x_labels = {}
    for c in fitness_cols:
        conc = c.replace("fitness_", "")
        x_labels[c] = f"AUC-F@{conc}"
    x_labels["mean_fitness"] = "Mean AUC-F"
    x_labels["PC1"] = "SVD/PC1"

    y_metrics = [
        ("log10_mic", "log₁₀(MIC_AUC)"),
        ("log10_ic50", "log₁₀(IC₅₀_AUC)"),
        ("log10_od_mic", "log₁₀(MIC_FinalOD)"),
        ("log10_od_ic50", "log₁₀(IC₅₀_FinalOD)"),
    ]

    rows = []
    for y_col, y_label in y_metrics:
        row_panels = []
        for x_col in x_cols:
            if x_col not in drug_data.columns or y_col not in drug_data.columns:
                continue
            cols = ["label", x_col, y_col]
            sub = (
                drug_data.select([c for c in cols if c in drug_data.columns])
                .drop_nulls()
            )
            if len(sub) < 4:
                continue
            sub_pd = sub.to_pandas()
            panel = _scatter_panel(sub_pd, x_col, y_col, x_labels[x_col], y_label, color_enc)
            row_panels.append(panel)

        # Chunk into sub-rows of 4 for A4 readability
        for i in range(0, len(row_panels), 4):
            chunk = row_panels[i:i + 4]
            rows.append(alt.hconcat(*chunk).resolve_scale(color="shared"))

    if not rows:
        return None

    return (
        alt.vconcat(*rows, make_cluster_legend(exclude={"DD"}))
        .properties(
            title=f"{drug}: Experimental Metrics vs Epistasis Fitness",
            padding={"bottom": 10, "right": 10},
        )
    )


def make_figure(suffix: str = ""):
    register_theme()
    xref_expanded_df = load_processed("xref_expanded_df", suffix)

    for drug in experiments:
        fig = _build_grid(drug, xref_expanded_df)
        if fig is None:
            print(f"WARNING: No data for {drug} scatter grid")
            continue
        stem = f"mic_scatter_grid_{drug.lower()}{suffix}"
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
