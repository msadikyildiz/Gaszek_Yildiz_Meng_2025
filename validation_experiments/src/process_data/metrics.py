"""IC50 metric computation."""

import polars as pl
import numpy as np
from scipy.optimize import curve_fit

from config import (
    experiments, TARGET_HOURS, utl, variant_label, MIC_THRESHOLD, MAX_CONC,
)


def _hill_func(x, bottom, top, ic50, hill):
    return bottom + (top - bottom) / (1.0 + (x / ic50) ** hill)


def compute_ic50(
    end_hour: float | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    active_hours = [h for h in TARGET_HOURS if end_hour is None or h <= end_hour]
    if not active_hours:
        return pl.DataFrame(), pl.DataFrame()

    ic50_rows = []
    for drug_abbrev, experiment in experiments.items():
        concs, _, _ = utl.form_concentrations(experiment)
        for plate_idx in range(len(experiment["plates"])):
            plate_variants = experiment["plates"][plate_idx]
            variant_rows_list = experiment["variants"]

            df_raw = utl.read_plate(experiment, plate_idx)
            df_proc = utl.subtract_bg_and_integrate(
                experiment, df_raw, drop_first=3,
            )
            t_hours = df_proc["t_hours"].to_numpy().astype(float)
            # Truncate to valid (monotonically increasing) data
            for _vi in range(1, len(t_hours)):
                if t_hours[_vi] <= t_hours[_vi - 1]:
                    t_hours = t_hours[:_vi]
                    break

            for slot_idx, variant_name in enumerate(plate_variants):
                letters = variant_rows_list[slot_idx]
                for target_h in active_hours:
                    if target_h > t_hours[-1]:
                        continue
                    t_idx = int(np.argmin(np.abs(t_hours - target_h)))
                    for rep_idx, letter in enumerate(letters):
                        ctrl_col = f"{letter}11_bgsub"
                        if ctrl_col not in df_proc.columns:
                            continue
                        ctrl_od = float(df_proc[ctrl_col][t_idx])
                        if ctrl_od <= 0.01:
                            continue
                        od_vals, conc_vals = [], []
                        for col_idx in range(10):
                            col_num = col_idx + 1
                            concentration = concs[col_idx]
                            if concentration <= 0:
                                continue
                            well = f"{letter}{col_num}"
                            bgsub_col = f"{well}_bgsub"
                            if bgsub_col not in df_proc.columns:
                                continue
                            od = float(df_proc[bgsub_col][t_idx])
                            od_vals.append(od / ctrl_od)
                            conc_vals.append(concentration)

                        if len(conc_vals) < 4:
                            continue

                        conc_arr = np.array(conc_vals)
                        od_arr = np.clip(np.array(od_vals), 0, 2.0)
                        min_c, max_c = conc_arr.min(), conc_arr.max()
                        try:
                            popt, _ = curve_fit(
                                _hill_func, conc_arr, od_arr,
                                p0=[0.0, 1.0, np.sqrt(min_c * max_c), 1.0],
                                bounds=(
                                    [0.0, 0.5, min_c / 10, 0.1],
                                    [0.5, 1.5, max_c * 10, 10.0],
                                ),
                                maxfev=5000,
                            )
                            ic50_val = float(popt[2])
                        except (RuntimeError, ValueError):
                            ic50_val = np.nan

                        ic50_rows.append({
                            "drug": drug_abbrev,
                            "variant": variant_name,
                            "target_hour": target_h,
                            "replicate": rep_idx,
                            "well_letter": letter,
                            "ic50_ugml": ic50_val,
                        })

    ic50_df = pl.DataFrame(ic50_rows)

    # Select the timepoint with the lowest mean coefficient of variation (CV)
    # across variants. Rationale: IC50 fits are noisy at early timepoints
    # (low OD → poor Hill fit) and at late timepoints (saturation → ambiguous
    # inflection). Min-CV picks the most reproducible window.
    # Candidate pool is small (4 timepoints: 6, 12, 18, 24 h) and selection
    # is per-drug, so overfitting risk is minimal.
    ic50_best = {}
    for drug in experiments:
        best_cv, best_h = np.inf, None
        cv_table = []
        for th in active_hours:
            sub = ic50_df.filter(
                (pl.col("drug") == drug)
                & (pl.col("target_hour") == th)
                & pl.col("ic50_ugml").is_not_nan(),
            )
            if len(sub) == 0:
                cv_table.append((th, 0, np.nan))
                continue
            var_stats = sub.group_by("variant").agg(
                pl.col("ic50_ugml").mean().alias("mean"),
                pl.col("ic50_ugml").std().alias("std"),
            ).with_columns(
                (pl.col("std") / pl.col("mean")).alias("cv"),
            )
            mean_cv = var_stats["cv"].drop_nulls().mean()
            cv_table.append((th, len(sub), mean_cv))
            if mean_cv is not None and mean_cv < best_cv:
                best_cv = mean_cv
                best_h = th
        ic50_best[drug] = best_h
        print(f"  IC50 timepoint CV ({drug}): "
              + ", ".join(f"{h}h: CV={cv:.3f} (n={n})"
                          if cv is not None and not np.isnan(cv)
                          else f"{h}h: n/a"
                          for h, n, cv in cv_table)
              + f"  → selected {best_h}h")

    ic50_metric_records = []
    for drug in experiments:
        bh = ic50_best[drug]
        if bh is None:
            continue
        sub = ic50_df.filter(
            (pl.col("drug") == drug)
            & (pl.col("target_hour") == bh)
            & pl.col("ic50_ugml").is_not_nan(),
        )
        for (variant,), vgrp in sub.group_by("variant"):
            vals = vgrp["ic50_ugml"].to_numpy()
            ic50_metric_records.append({
                "drug": drug,
                "variant": variant,
                "ic50_mean": float(np.nanmean(vals)),
                "ic50_std": float(np.nanstd(vals)),
                "ic50_n": int(np.sum(~np.isnan(vals))),
                "best_timepoint_h": bh,
            })

    ic50_metric_df = pl.DataFrame(ic50_metric_records)
    print(f"ic50_metric_df: {ic50_metric_df.shape}")

    # Per-replicate IC50: filter to best timepoint per drug, exclude NaN fits
    ic50_rep_rows = []
    for drug in experiments:
        bh = ic50_best.get(drug)
        if bh is None:
            continue
        sub = ic50_df.filter(
            (pl.col("drug") == drug)
            & (pl.col("target_hour") == bh)
            & pl.col("ic50_ugml").is_not_nan(),
        )
        for row in sub.iter_rows(named=True):
            ic50_rep_rows.append({
                "drug": row["drug"],
                "variant": row["variant"],
                "rep_letter": row["well_letter"],
                "ic50_ugml": row["ic50_ugml"],
                "log10_ic50": float(np.log10(row["ic50_ugml"])),
                "label": variant_label(row["variant"]),
            })

    ic50_rep_df = pl.DataFrame(ic50_rep_rows)
    print(f"ic50_rep_df: {ic50_rep_df.shape}")
    return ic50_metric_df, ic50_rep_df


def compute_ic50_from_long(
    long_df: pl.DataFrame, metric_col: str = "final_auc",
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Per-replicate Hill fits using long_df + specified metric column."""
    tag = metric_col.replace("final_", "")
    ic50_rows = []

    for (drug, variant), grp in long_df.group_by(["drug", "variant"]):
        ctrl = grp.filter(pl.col("is_control"))
        ctrl_median = ctrl[metric_col].median()
        if ctrl_median is None or ctrl_median < 0.05:
            continue

        drug_wells = grp.filter(~pl.col("is_control")).with_columns(
            pl.col("well").str.slice(0, 1).alias("rep_letter"),
            (pl.col(metric_col) / ctrl_median).alias("norm_val"),
        )

        for (rep_letter,), rep_grp in drug_wells.group_by("rep_letter"):
            rep = rep_grp.sort("concentration")
            conc_arr = rep["concentration"].to_numpy()
            val_arr = np.clip(rep["norm_val"].to_numpy(), 0, 2.0)
            mask = conc_arr > 0
            if mask.sum() < 4:
                continue
            c = conc_arr[mask]
            v = val_arr[mask]
            min_c, max_c = float(c.min()), float(c.max())
            try:
                popt, _ = curve_fit(
                    _hill_func, c, v,
                    p0=[0.0, 1.0, np.sqrt(min_c * max_c), 1.0],
                    bounds=(
                        [0.0, 0.5, min_c / 10, 0.1],
                        [0.5, 1.5, max_c * 10, 10.0],
                    ),
                    maxfev=5000,
                )
                ic50_val = float(popt[2])
            except (RuntimeError, ValueError):
                ic50_val = np.nan

            ic50_rows.append({
                "drug": drug,
                "variant": variant,
                "rep_letter": rep_letter,
                "ic50_ugml": ic50_val,
            })

    if not ic50_rows:
        empty_metric = pl.DataFrame(schema={
            "drug": pl.Utf8, "variant": pl.Utf8,
            "ic50_mean": pl.Float64, "ic50_std": pl.Float64,
            "ic50_n": pl.Int64,
        })
        empty_rep = pl.DataFrame(schema={
            "drug": pl.Utf8, "variant": pl.Utf8,
            "rep_letter": pl.Utf8, "ic50_ugml": pl.Float64,
            "log10_ic50": pl.Float64, "label": pl.Utf8,
        })
        return empty_metric, empty_rep

    raw_df = pl.DataFrame(ic50_rows)

    # Aggregate per variant
    ic50_metric_records = []
    for (drug, variant), vgrp in raw_df.filter(
        pl.col("ic50_ugml").is_not_nan()
    ).group_by(["drug", "variant"]):
        vals = vgrp["ic50_ugml"].to_numpy()
        ic50_metric_records.append({
            "drug": drug, "variant": variant,
            "ic50_mean": float(np.nanmean(vals)),
            "ic50_std": float(np.nanstd(vals)),
            "ic50_n": int(np.sum(~np.isnan(vals))),
        })

    ic50_metric_df = pl.DataFrame(ic50_metric_records)

    ic50_rep_df = (
        raw_df
        .filter(pl.col("ic50_ugml").is_not_nan())
        .with_columns(
            pl.col("ic50_ugml").log(base=10).alias("log10_ic50"),
            pl.col("variant").map_elements(
                variant_label, return_dtype=pl.Utf8,
            ).alias("label"),
        )
    )

    print(f"ic50_from_long({tag}): metric={ic50_metric_df.shape}, rep={ic50_rep_df.shape}")
    return ic50_metric_df, ic50_rep_df
