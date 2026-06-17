"""Correlation computation and growth curve caching (cell 15 data + new)."""

import warnings

import polars as pl
import numpy as np
from scipy import stats

from config import (
    experiments, variant_label, utl,
)


def _sig_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def compute_correlations(xref_expanded_df: pl.DataFrame) -> pl.DataFrame:
    corr_source = xref_expanded_df.filter(
        pl.col("genotype_13").is_not_null(),
    )

    Y_METRICS = {
        "log₁₀(MIC_AUC)": "log10_mic",
        "log₁₀(IC₅₀_AUC)": "log10_ic50",
        "log₁₀(MIC_FinalOD)": "log10_od_mic",
        "log₁₀(IC₅₀_FinalOD)": "log10_od_ic50",
    }

    corr_records = []
    for drug in experiments:
        drug_df = corr_source.filter(pl.col("drug") == drug)

        x_metrics = {"Mean AUC-F": "mean_fitness"}
        if "PC1" in drug_df.columns:
            x_metrics["SVD/PC1"] = "PC1"
        for col in drug_df.columns:
            if (
                col.startswith("fitness_")
                and col != "fitness_0.0"
                and col != "fitness_0"
            ):
                conc_str = col.replace("fitness_", "")
                x_metrics[f"AUC-F@{conc_str}"] = col

        for y_label, y_col in Y_METRICS.items():
            if y_col not in drug_df.columns:
                continue
            y_vals = drug_df[y_col].to_numpy().astype(float)
            y_valid = ~np.isnan(y_vals)

            for x_label, x_col in x_metrics.items():
                if x_col not in drug_df.columns:
                    continue
                x_vals = drug_df[x_col].to_numpy().astype(float)
                x_valid = ~np.isnan(x_vals)

                mask = y_valid & x_valid
                if mask.sum() < 4:
                    continue

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    r, p = stats.pearsonr(x_vals[mask], y_vals[mask])

                corr_records.append({
                    "drug": drug,
                    "y_metric": y_label,
                    "x_metric": x_label,
                    "r": float(r),
                    "p": float(p),
                    "n": int(mask.sum()),
                    "stars": _sig_stars(p),
                    "r_text": f"{r:.2f}{_sig_stars(p)}",
                })

    corr_results_df = pl.DataFrame(corr_records)
    print(f"corr_results_df: {corr_results_df.shape}")
    return corr_results_df


def cache_growth_curves(end_hour: float | None = None) -> pl.DataFrame:
    """Pre-extract smoothed bg-subtracted curves for fig_d."""
    SMOOTH_WINDOW = 6
    gc_rows = []

    for drug_abbrev in experiments:
        experiment = experiments[drug_abbrev]
        concs, _, _ = utl.form_concentrations(experiment)

        for plate_idx in range(len(experiment["plates"])):
            plate_variants = experiment["plates"][plate_idx]
            variant_rows_list = experiment["variants"]

            df_raw = utl.read_plate(experiment, plate_idx)
            df_proc = utl.subtract_bg_and_integrate(
                experiment, df_raw, drop_first=3,
            )
            t_hours_full = df_proc["t_hours"].to_numpy().astype(float)
            # Truncate to valid (monotonically increasing) data
            valid_len = len(t_hours_full)
            for _vi in range(1, len(t_hours_full)):
                if t_hours_full[_vi] <= t_hours_full[_vi - 1]:
                    valid_len = _vi
                    break
            t_hours_full = t_hours_full[:valid_len]
            if end_hour is not None:
                cut_idx = int(np.argmin(np.abs(t_hours_full - end_hour))) + 1
                t_hours = t_hours_full[:cut_idx]
            else:
                t_hours = t_hours_full

            for slot_idx, variant_name in enumerate(plate_variants):
                letters = variant_rows_list[slot_idx]
                for col_idx in range(11):
                    col_num = col_idx + 1
                    concentration = concs[col_idx]
                    for letter in letters:
                        well = f"{letter}{col_num}"
                        bgsub_col = f"{well}_bgsub"
                        if bgsub_col not in df_proc.columns:
                            continue
                        curve = df_proc[bgsub_col].to_numpy().astype(float)[:len(t_hours)]

                        if len(curve) >= SMOOTH_WINDOW:
                            kernel = np.ones(SMOOTH_WINDOW) / SMOOTH_WINDOW
                            smoothed = np.convolve(curve, kernel, mode="valid")
                            t_smooth = t_hours[SMOOTH_WINDOW - 1:]
                        else:
                            smoothed = curve
                            t_smooth = t_hours

                        for i in range(len(smoothed)):
                            gc_rows.append({
                                "drug": drug_abbrev,
                                "variant": variant_name,
                                "label": variant_label(variant_name),
                                "concentration": concentration,
                                "time_h": t_smooth[i],
                                "od_bgsub": smoothed[i],
                                "rep_id": f"{plate_idx}_{letter}",
                            })

    gc_df = pl.DataFrame(gc_rows)
    print(f"growth_curves: {gc_df.shape}")
    return gc_df
