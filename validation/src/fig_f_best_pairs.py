"""Figure F: Best-Pair Scatter Plots (cell 16)."""

import numpy as np
import pandas as pd
import polars as pl
import altair as alt

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


def _vl_field(name: str) -> str:
    """Escape field name for Vega-Lite path parsing (e.g., fitness_0.0)."""
    return name.replace(".", r"\.")


def make_figure(suffix: str = ""):
    register_theme()
    corr_results_df = load_processed("corr_results_df", suffix)
    xref_expanded_df = load_processed("xref_expanded_df", suffix)

    color_domain = [v for v in VARIANT_ORDER if v != "DD"]
    color_range = [VARIANT_COLORS[v] for v in color_domain]

    color_enc = alt.Color(
        "label:N",
        scale=alt.Scale(domain=color_domain, range=color_range),
        legend=None,
    )

    y_col_map = {
        "log₁₀(MIC_AUC)": "log10_mic",
        "log₁₀(IC₅₀_AUC)": "log10_ic50",
        "log₁₀(MIC_FinalOD)": "log10_od_mic",
        "log₁₀(IC₅₀_FinalOD)": "log10_od_ic50",
    }
    x_col_map = {"Mean AUC-F": "mean_fitness", "SVD/PC1": "PC1"}
    for col in xref_expanded_df.columns:
        if col.startswith("fitness_") and col not in ("fitness_0.0", "fitness_0"):
            conc_str = col.replace("fitness_", "")
            x_col_map[f"AUC-F@{conc_str}"] = col

    scatter_panels = []

    for y_label in y_col_map:
        y_col = y_col_map[y_label]
        row_panels = []

        for drug in experiments:
            sub_corr = corr_results_df.filter(
                (pl.col("drug") == drug) & (pl.col("y_metric") == y_label),
            ).sort(pl.col("r").abs(), descending=True)

            if len(sub_corr) == 0:
                continue

            best = sub_corr.row(0, named=True)
            x_label = best["x_metric"]
            x_col = x_col_map.get(x_label)
            r_val, p_val = best["r"], best["p"]

            if x_col is None:
                continue

            drug_data = xref_expanded_df.filter(
                (pl.col("drug") == drug)
                & pl.col("genotype_13").is_not_null(),
            )

            plot_cols = ["label", x_col, y_col]
            drug_pd = (
                drug_data
                .select([c for c in plot_cols if c in drug_data.columns])
                .drop_nulls()
                .to_pandas()
            )

            if len(drug_pd) < 3:
                continue

            drug_pd["short_name"] = drug_pd["label"].apply(extract_short_name)

            significant = p_val < 0.05
            x_arr = drug_pd[x_col].to_numpy().astype(float)
            y_arr = drug_pd[y_col].to_numpy().astype(float)
            slope, intercept = np.polyfit(x_arr, y_arr, 1)
            x_line = np.array([x_arr.min(), x_arr.max()])
            y_line = slope * x_line + intercept
            reg_df = pd.DataFrame({x_col: x_line, y_col: y_line})
            x_field = _vl_field(x_col)
            y_field = _vl_field(y_col)

            # Nudge label positions (shared helper)
            x_range = x_arr.max() - x_arr.min() or 1.0
            y_range = y_arr.max() - y_arr.min() or 1.0
            lx, ly = nudge_labels(
                x_arr.copy(), y_arr.copy(), x_range, y_range,
                all_pts_x=x_arr, all_pts_y=y_arr,
            )
            drug_pd["_lx"] = lx
            drug_pd["_ly"] = ly

            # Explicit domains with padding
            x_pad = max(x_range * 0.08, 0.05)
            y_pad = max(y_range * 0.08, 0.05)
            x_scale = alt.Scale(domain=[float(x_arr.min() - x_pad), float(x_arr.max() + x_pad)])
            y_scale = alt.Scale(domain=[float(y_arr.min() - y_pad), float(y_arr.max() + y_pad)])

            points = (
                alt.Chart(drug_pd)
                .mark_point(
                    size=120, filled=True, opacity=0.85,
                    stroke="#FAFAFA", strokeWidth=1.5,
                )
                .encode(
                    x=alt.X(f"{x_field}:Q", title=x_label, scale=x_scale),
                    y=alt.Y(f"{y_field}:Q", title=y_label, scale=y_scale),
                    color=color_enc,
                )
            )

            labels = (
                alt.Chart(drug_pd)
                .mark_text(
                    fontSize=8, font=FONT_MONO,
                    align="left", color="#424242",
                )
                .encode(
                    x=alt.X("_lx:Q", scale=x_scale),
                    y=alt.Y("_ly:Q", scale=y_scale),
                    text="short_name:N",
                )
            )

            conn = make_connectors(
                drug_pd, x_col, y_col,
                x_range=x_range, y_range=y_range,
            )

            reg = (
                alt.Chart(reg_df)
                .mark_line(
                    strokeDash=[6, 3],
                    strokeWidth=1.5 if significant else 1.0,
                    opacity=0.5 if significant else 0.2,
                    color="#424242" if significant else "#BDBDBD",
                )
                .encode(
                    x=alt.X(f"{x_field}:Q", scale=x_scale),
                    y=alt.Y(f"{y_field}:Q", scale=y_scale),
                )
            )

            # r,p as text annotation at top-right of chart
            stars = _sig_stars(p_val)
            ann_color = "#1A1A1A" if significant else "#9E9E9E"
            ann_fw = "bold" if significant else "normal"
            ann_text = f"r={r_val:.2f}, p={p_val:.3f}{stars}"
            rp_data = pd.DataFrame({
                "ax": [float(x_arr.max() + x_pad * 0.5)],
                "ay": [float(y_arr.max() + y_pad * 0.8)],
                "t": [ann_text],
            })
            rp_ann = (
                alt.Chart(rp_data)
                .mark_text(
                    align="right", baseline="top",
                    fontSize=10, font=FONT_MONO,
                    fontWeight=ann_fw, color=ann_color,
                )
                .encode(
                    x=alt.X("ax:Q", scale=x_scale),
                    y=alt.Y("ay:Q", scale=y_scale),
                    text="t:N",
                )
            )

            panel_title = f"{drug}: {x_label}"

            layers = [reg]
            if conn is not None:
                layers.append(conn)
            layers += [points, labels, rp_ann]

            panel = alt.layer(*layers).properties(
                width=280, height=240, title=panel_title,
            )
            row_panels.append(panel)

        if row_panels:
            scatter_panels.append(alt.hconcat(*row_panels))

    if scatter_panels:
        fig = (
            alt.vconcat(*scatter_panels, make_cluster_legend(exclude={"DD"}))
            .resolve_scale(color="shared")
            .properties(
                title="Best Correlating Metric Pairs",
                padding={"bottom": 10, "right": 10},
            )
        )

        stem = f"mic_best_correlations{suffix}"
        fig.save(str(FIGURES_DIR / f"{stem}.png"), scale_factor=2)
        fig.save(str(FIGURES_DIR / f"{stem}.html"))
        print(f"Saved {stem}")
        return fig
    else:
        print("WARNING: No scatter panels (insufficient data)")
        return None


if __name__ == "__main__":
    import sys
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    s = ""
    if "--suffix" in sys.argv:
        s = sys.argv[sys.argv.index("--suffix") + 1]
    make_figure(suffix=s)
