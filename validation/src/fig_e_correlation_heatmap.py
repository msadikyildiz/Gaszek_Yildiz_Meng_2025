"""Figure E: Correlation Heatmap (cell 15, plot half)."""

import altair as alt

from config import (
    load_processed, register_theme, FIGURES_DIR, FONT_MONO, experiments,
)


def make_figure(suffix: str = ""):
    register_theme()
    corr_results_df = load_processed("corr_results_df", suffix)
    corr_pd = corr_results_df.to_pandas()

    corr_pd["r_val_text"] = corr_pd["r"].apply(lambda r: f"{r:.2f}")
    corr_pd["stars_text"] = corr_pd["stars"]

    x_order_base = ["Mean AUC-F", "SVD/PC1", "SVD/PC2"]
    x_fitness_cols = sorted(
        [x for x in corr_pd["x_metric"].unique() if x.startswith("AUC-F@")],
        key=lambda s: float(s.replace("AUC-F@", "")),
    )
    x_order = [
        x for x in x_order_base if x in corr_pd["x_metric"].values
    ] + x_fitness_cols
    y_order = ["log₁₀(MIC_AUC)", "log₁₀(IC₅₀_AUC)", "log₁₀(MIC_FinalOD)", "log₁₀(IC₅₀_FinalOD)"]

    heatmap_panels = []
    for drug in experiments:
        sub = corr_pd[corr_pd["drug"] == drug]

        base = alt.Chart(sub).encode(
            x=alt.X(
                "x_metric:N", sort=x_order, title=None,
                axis=alt.Axis(labelAngle=-45, labelFontSize=10),
            ),
            y=alt.Y(
                "y_metric:N", sort=y_order, title=None,
                axis=alt.Axis(labelFontSize=11),
            ),
        )

        rect = base.mark_rect(stroke="#FAFAFA", strokeWidth=2).encode(
            color=alt.Color(
                "r:Q",
                scale=alt.Scale(
                    domain=[-1, 0, 1],
                    range=["#1565C0", "#FAFAFA", "#C62828"],
                ),
                legend=alt.Legend(title="Pearson r", orient="right"),
            ),
        )

        text_r = base.mark_text(
            fontSize=10, font=FONT_MONO, fontWeight="bold", dy=-5,
        ).encode(
            text="r_val_text:N",
            color=alt.condition(
                "abs(datum.r) > 0.5",
                alt.value("white"),
                alt.value("#1A1A1A"),
            ),
        )

        text_stars = base.mark_text(
            fontSize=8, font=FONT_MONO, dy=7,
        ).encode(
            text="stars_text:N",
            color=alt.condition(
                "abs(datum.r) > 0.5",
                alt.value("white"),
                alt.value("#1A1A1A"),
            ),
        )

        panel = (rect + text_r + text_stars).properties(
            width=alt.Step(55), height=150, title=drug,
        )
        heatmap_panels.append(panel)

    fig = (
        alt.vconcat(*heatmap_panels)
        .properties(
            title="Correlation: Experimental Metrics vs Epistasis Predictors",
            padding={"bottom": 10, "right": 10},
        )
    )

    stem = f"mic_correlation_heatmap{suffix}"
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
