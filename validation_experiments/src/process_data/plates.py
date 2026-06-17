"""Process plates and compute MIC (cells 4-5)."""

import polars as pl
import numpy as np

from config import (
    experiments, variant_label, MIC_THRESHOLD, utl,
)


def _valid_length(t_hours: np.ndarray) -> int:
    """Find last index where t_hours is strictly increasing (skip padding)."""
    for i in range(1, len(t_hours)):
        if t_hours[i] <= t_hours[i - 1]:
            return i
    return len(t_hours)


def process_plates(end_hour: float | None = None) -> pl.DataFrame:
    rows = []
    for drug_abbrev, experiment in experiments.items():
        concs, _, _ = utl.form_concentrations(experiment)
        n_plates = len(experiment["plates"])

        for plate_idx in range(n_plates):
            plate_variants = experiment["plates"][plate_idx]
            variant_rows_list = experiment["variants"]

            df_raw = utl.read_plate(experiment, plate_idx)
            df_proc = utl.subtract_bg_and_integrate(
                experiment, df_raw, drop_first=3,
            )

            t_hours = df_proc["t_hours"].to_numpy().astype(float)
            valid_len = _valid_length(t_hours)

            # Determine cutoff index for end_hour truncation
            if end_hour is not None:
                cut_idx = int(np.argmin(np.abs(t_hours[:valid_len] - end_hour)))
            else:
                cut_idx = valid_len - 1

            last5_start = max(0, cut_idx + 1 - 5)
            last5_slice = slice(last5_start, cut_idx + 1)

            for slot_idx, variant_name in enumerate(plate_variants):
                letters = variant_rows_list[slot_idx]
                for col_idx in range(11):
                    col_num = col_idx + 1
                    concentration = concs[col_idx]
                    is_control = concentration == 0.0

                    for rep_idx, letter in enumerate(letters):
                        well = f"{letter}{col_num}"
                        auc_col = f"{well}_auc"
                        bgsub_col = f"{well}_bgsub"
                        if auc_col not in df_proc.columns:
                            continue

                        final_auc = df_proc[auc_col][cut_idx]
                        final_od = df_proc[bgsub_col][last5_slice].mean()

                        rows.append({
                            "drug": drug_abbrev,
                            "plate": plate_idx + 1,
                            "variant": variant_name,
                            "concentration": concentration,
                            "is_control": is_control,
                            "replicate": rep_idx,
                            "well": well,
                            "final_auc": float(final_auc),
                            "final_od": float(final_od),
                        })

    long_df = pl.DataFrame(rows)
    print(f"long_df: {long_df.shape}")
    return long_df


def compute_mic(long_df: pl.DataFrame, metric_col: str = "final_auc"):
    norm_col = "norm_auc" if metric_col == "final_auc" else "norm_od"
    mic_records = []
    dr_records = []
    mic_rep_records = []

    for (drug, variant), grp in long_df.group_by(["drug", "variant"]):
        ctrl = grp.filter(pl.col("is_control"))
        control_median = ctrl[metric_col].median()
        if control_median is None or control_median < 0.05:
            continue

        drug_wells = grp.filter(~pl.col("is_control")).with_columns(
            (pl.col(metric_col) / control_median).alias(norm_col),
        )

        conc_stats = (
            drug_wells
            .group_by("concentration")
            .agg(
                pl.col(norm_col).mean().alias(f"{norm_col}_mean"),
                pl.col(norm_col).std().alias(f"{norm_col}_std"),
                pl.col(norm_col).count().alias("n_reps"),
            )
            .sort("concentration")
        )

        concs_arr = conc_stats["concentration"].to_numpy()
        means_arr = conc_stats[f"{norm_col}_mean"].to_numpy()
        stds_arr = conc_stats[f"{norm_col}_std"].to_numpy()

        for i in range(len(concs_arr)):
            if concs_arr[i] == 0.0:
                continue
            std_val = stds_arr[i]
            if std_val is None or np.isnan(std_val):
                std_val = 0.0
            dr_records.append({
                "drug": drug,
                "variant": variant,
                "concentration": concs_arr[i],
                f"{norm_col}_mean": means_arr[i],
                f"{norm_col}_std": std_val,
                f"{norm_col}_lo": means_arr[i] - std_val,
                f"{norm_col}_hi": means_arr[i] + std_val,
            })

        # Zero-drug control point (pseudo-conc one dilution step below min)
        ctrl_norm = ctrl[metric_col].to_numpy() / control_median
        sorted_dc = np.sort(concs_arr[concs_arr > 0])
        if len(sorted_dc) >= 2:
            pseudo_conc = sorted_dc[0] ** 2 / sorted_dc[1]
        else:
            pseudo_conc = sorted_dc[0] / 3.0
        ctrl_mean = float(ctrl_norm.mean())
        ctrl_std = float(ctrl_norm.std()) if len(ctrl_norm) > 1 else 0.0
        dr_records.append({
            "drug": drug, "variant": variant,
            "concentration": float(pseudo_conc),
            f"{norm_col}_mean": ctrl_mean,
            f"{norm_col}_std": ctrl_std,
            f"{norm_col}_lo": ctrl_mean - ctrl_std,
            f"{norm_col}_hi": ctrl_mean + ctrl_std,
        })

        mask = concs_arr > 0
        concs_drug = concs_arr[mask]
        means_drug = means_arr[mask]
        if len(concs_drug) > 1:
            is_censored = bool(means_drug.min() >= MIC_THRESHOLD)
            mic_val = np.interp(
                MIC_THRESHOLD,
                means_drug[::-1],
                concs_drug[::-1],
            )
        else:
            mic_val = np.nan
            is_censored = True

        mic_records.append({
            "drug": drug,
            "variant": variant,
            "mic_ugml": float(mic_val),
            "is_censored": is_censored,
            f"control_median_{metric_col}": float(control_median),
        })

        # Per-replicate MIC
        drug_wells2 = grp.filter(~pl.col("is_control")).with_columns(
            pl.col("well").str.slice(0, 1).alias("rep_letter"),
            (pl.col(metric_col) / control_median).alias(norm_col),
        )
        for (rep_letter,), rep_grp in drug_wells2.group_by("rep_letter"):
            rep_conc = rep_grp.sort("concentration")
            rc = rep_conc["concentration"].to_numpy()
            rv = rep_conc[norm_col].to_numpy()
            m = rc > 0
            if m.sum() > 1:
                rep_censored = bool(rv[m].min() >= MIC_THRESHOLD)
                mic_r = np.interp(MIC_THRESHOLD, rv[m][::-1], rc[m][::-1])
            else:
                mic_r = np.nan
                rep_censored = True
            mic_rep_records.append({
                "drug": drug, "variant": variant,
                "rep_letter": rep_letter, "mic_ugml": float(mic_r),
                "is_censored": rep_censored,
            })

    mic_df = pl.DataFrame(mic_records)
    dr_df = pl.DataFrame(dr_records).with_columns(
        pl.col("variant").map_elements(
            variant_label, return_dtype=pl.Utf8,
        ).alias("label"),
    )
    mic_rep_df = pl.DataFrame(mic_rep_records).with_columns(
        pl.col("mic_ugml").clip(lower_bound=1e-3).log(base=10).alias("log10_mic"),
        pl.col("variant").map_elements(
            variant_label, return_dtype=pl.Utf8,
        ).alias("label"),
    )
    mic_df = mic_df.with_columns(
        pl.col("variant").map_elements(
            variant_label, return_dtype=pl.Utf8,
        ).alias("label"),
    )

    tag = metric_col.replace("final_", "")
    print(f"mic_df({tag}): {mic_df.shape}, dr_df: {dr_df.shape}, mic_rep_df: {mic_rep_df.shape}")
    return mic_df, dr_df, mic_rep_df
