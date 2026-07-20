"""Figure H: Library Fitness — estimated cell counts from pooled sequencing."""

import sys
import numpy as np
import polars as pl
import altair as alt

from config import (
    load_processed, register_theme, FIGURES_DIR, FONT_MONO, experiments,
    VARIANT_ORDER, VARIANT_COLORS,
    make_cluster_legend,
)

# Map viewer suffix → sequencing timepoint (available: 3, 6, 9, 12)
SUFFIX_TO_TIMEPOINT = {
    "_12h": 12.0,
    "_24h": 12.0,
    "": 12.0,
}


def make_figure(suffix: str = ""):
    register_theme()

    # Cell counts are computed once from sequencing — no suffix on parquet
    cc_df = load_processed("cell_counts_df")

    tp = SUFFIX_TO_TIMEPOINT.get(suffix, 12.0)

    variant_order_no_dd = [v for v in VARIANT_ORDER if v != "DD"]
    color_domain = list(variant_order_no_dd)
    color_range = [VARIANT_COLORS[v] for v in color_domain]

    for drug in experiments:
        sub = cc_df.filter(
            (pl.col("drug") == drug) & (pl.col("timepoint_h") == tp),
        )

        if len(sub) == 0:
            # Fall back to closest available timepoint
            available = sorted(
                cc_df.filter(pl.col("drug") == drug)["timepoint_h"]
                .unique().to_list()
            )
            tp_actual = min(available, key=lambda t: abs(t - tp))
            sub = cc_df.filter(
                (pl.col("drug") == drug) & (pl.col("timepoint_h") == tp_actual),
            )

        # Average over replicates
        agg = (
            sub.group_by(["variant", "label", "concentration"])
            .agg(
                pl.col("estimated_cell_count").mean().alias("cell_count_mean"),
                pl.col("estimated_cell_count").std().alias("cell_count_std"),
                pl.col("read_fraction").mean().alias("read_frac_mean"),
            )
            .filter(pl.col("concentration") > 0)
            .with_columns(
                (pl.col("cell_count_mean") * 1e8).alias("cells"),
            )
            .sort("variant", "concentration")
        ).to_pandas()

        if agg.empty:
            print(f"No data for {drug} at timepoint {tp}, skipping")
            continue

        x_min = np.log10(agg["concentration"].min()) - 0.2
        x_max = np.log10(agg["concentration"].max()) + 0.2
        agg["log_conc"] = np.log10(agg["concentration"])
        # Clamp cells to positive for log scale
        agg["cells"] = agg["cells"].clip(lower=1.0)

        lines = (
            alt.Chart(agg)
            .mark_line(strokeWidth=1.8)
            .encode(
                x=alt.X("log_conc:Q", title="log₁₀ Conc (µg/mL)",
                         scale=alt.Scale(domain=[x_min, x_max])),
                y=alt.Y("cells:Q",
                         title="Estimated Cells (1 OD₆₀₀ = 10⁸ cells)",
                         scale=alt.Scale(type="log")),
                color=alt.Color(
                    "label:N",
                    scale=alt.Scale(domain=color_domain, range=color_range),
                    legend=None,
                ),
            )
        )

        dots = (
            alt.Chart(agg)
            .mark_point(size=35, filled=True)
            .encode(
                x="log_conc:Q",
                y="cells:Q",
                color=alt.Color(
                    "label:N",
                    scale=alt.Scale(domain=color_domain, range=color_range),
                    legend=None,
                ),
            )
        )

        main = (
            (lines + dots)
            .properties(width=400, height=280)
        )

        fig = (
            alt.vconcat(main, make_cluster_legend(exclude={"DD"}))
            .properties(
                title=f"Library Fitness — {drug} (t = {tp:.0f}h)",
                padding={"bottom": 10, "right": 10},
            )
        )

        fname = f"library_fitness_{drug.lower()}{suffix}"
        fig.save(str(FIGURES_DIR / f"{fname}.png"), scale_factor=2)
        print(f"Saved {fname}")


if __name__ == "__main__":
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    s = ""
    if "--suffix" in sys.argv:
        s = sys.argv[sys.argv.index("--suffix") + 1]
    make_figure(suffix=s)
