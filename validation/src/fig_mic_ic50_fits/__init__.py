"""Unified MIC/IC50 Fits — growth grid + dose-response per variant."""

import polars as pl
import altair as alt

from config import (
    load_processed, register_theme, FIGURES_DIR, FONT_BODY, experiments,
    SHORT_COLORS, MUTATION_CUES, _NAME_ORDER,
)
from .charts import growth_grid, dose_response_hill

_ATLAS_ORDER = list(_NAME_ORDER)


def _build_row(gc_pd, dr_pd, variant, show_titles, od_dr_pd=None):
    """Single variant row: growth grid | AUC dose-response | OD dose-response."""
    color = SHORT_COLORS[variant]

    left = growth_grid(gc_pd, color, n_cols=4)
    auc_panel, _, _ = dose_response_hill(
        dr_pd, color, norm_col="norm_auc", y_title="Normalized AUC",
    )

    parts = [left, auc_panel]
    if od_dr_pd is not None and len(od_dr_pd) > 0:
        od_panel, _, _ = dose_response_hill(
            od_dr_pd, color, norm_col="norm_od", y_title="Normalized OD",
        )
        parts.append(od_panel)

    cue = MUTATION_CUES.get(variant, "")
    return (
        alt.hconcat(*parts)
        .properties(
            spacing=12,
            title=alt.Title(
                f"{variant} — {cue}",
                fontSize=11, font=FONT_BODY, anchor="start",
            ),
        )
    )


def make_figure(suffix=""):
    register_theme()
    gc_raw = load_processed("growth_curves", suffix)
    dr_raw = load_processed("dr_df", suffix)
    od_dr_raw = load_processed("od_dr_df", suffix)

    for drug in experiments:
        rows = []
        for i, variant in enumerate(_ATLAS_ORDER):
            gc_df = gc_raw.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            dr_df = dr_raw.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            od_dr_df = od_dr_raw.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            if gc_df.height == 0 or dr_df.height == 0:
                continue

            row = _build_row(
                gc_df.to_pandas(), dr_df.to_pandas(),
                variant, show_titles=(i == 0),
                od_dr_pd=od_dr_df.to_pandas(),
            )
            rows.append(row)

        fig = alt.vconcat(*rows).properties(spacing=18)

        stem = f"mic_ic50_fits_{drug.lower()}{suffix}"
        fig.save(str(FIGURES_DIR / f"{stem}.png"), scale_factor=2)
        fig.save(str(FIGURES_DIR / f"{stem}.html"))
        print(f"Saved {stem}")
