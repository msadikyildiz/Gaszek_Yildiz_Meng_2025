"""Figure A: Dose-Response Curves (cell 12)."""

import altair as alt

from config import (
    load_processed, register_theme, FIGURES_DIR,
    VARIANT_ORDER, VARIANT_COLORS, MIC_THRESHOLD,
    make_cluster_legend,
)


def make_figure(suffix: str = ""):
    register_theme()
    dr_df = load_processed("dr_df", suffix)

    dr_plot = dr_df.to_pandas()

    color_domain = list(VARIANT_ORDER)
    color_range = [VARIANT_COLORS[v] for v in color_domain]
    color_scale = alt.Scale(domain=color_domain, range=color_range)

    band = alt.Chart(dr_plot).mark_area(opacity=0.15).encode(
        x=alt.X(
            "concentration:Q",
            scale=alt.Scale(type="log"),
            title="Concentration (µg/mL)",
            axis=alt.Axis(grid=False),
        ),
        y=alt.Y("norm_auc_lo:Q", title="Normalized AUC",
                 axis=alt.Axis(grid=False)),
        y2="norm_auc_hi:Q",
        color=alt.Color("label:N", scale=color_scale, legend=None),
    )

    line = alt.Chart(dr_plot).mark_line(strokeWidth=1.8).encode(
        x=alt.X(
            "concentration:Q",
            scale=alt.Scale(type="log"),
            title="Concentration (µg/mL)",
            axis=alt.Axis(grid=False),
        ),
        y=alt.Y("norm_auc_mean:Q", title="Normalized AUC",
                 axis=alt.Axis(grid=False)),
        color=alt.Color("label:N", scale=color_scale, legend=None),
    )

    threshold = alt.Chart(dr_plot).mark_rule(
        strokeDash=[6, 4], strokeWidth=1, color="#616161",
    ).encode(y=alt.datum(MIC_THRESHOLD))

    main = (
        (band + line + threshold)
        .properties(width=300, height=250)
        .facet(column=alt.Column(
            "drug:N", title=None,
            header=alt.Header(labelFontSize=14),
        ))
        .resolve_scale(x="independent")
    )

    fig = (
        alt.vconcat(main, make_cluster_legend())
        .properties(
            title="Dose-Response Curves",
            padding={"bottom": 10, "right": 10},
        )
    )

    stem = f"mic_dose_response{suffix}"
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
