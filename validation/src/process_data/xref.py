"""Cross-reference and expanded cross-reference tables (cells 8-10)."""

import polars as pl
import numpy as np

from config import (
    variants_meta, GENO_LOOKUP, variant_label,
    DRUG_MAP, MAX_CONC,
)


def cross_reference(mic_df: pl.DataFrame, epi_df: pl.DataFrame):
    xref_records = []
    for row in mic_df.iter_rows(named=True):
        drug = row["drug"]
        variant = row["variant"]
        genotype = GENO_LOOKUP.get(variant)

        if genotype is None:
            xref_records.append({
                **row,
                "genotype_13": None,
                "n_mutations": len(variants_meta[variant]["mutations"]),
                "fitness_baseline": None,
                "fitness_max_conc": None,
                "resistance_index": None,
            })
            continue

        fit_rows = epi_df.filter(
            (pl.col("Genotype") == genotype)
            & (pl.col("Drug") == DRUG_MAP[drug]),
        ).sort("Concentration")

        fitness_baseline = fit_rows.filter(
            pl.col("Concentration") == 0.0,
        )["Fitness"].item()
        fitness_max = fit_rows.filter(
            pl.col("Concentration") == MAX_CONC[drug],
        )["Fitness"].item()
        resistance_index = (
            fitness_max / fitness_baseline if fitness_baseline != 0 else None
        )

        xref_records.append({
            **row,
            "genotype_13": genotype,
            "n_mutations": len(variants_meta[variant]["mutations"]),
            "fitness_baseline": float(fitness_baseline),
            "fitness_max_conc": float(fitness_max),
            "resistance_index": (
                float(resistance_index) if resistance_index is not None else None
            ),
        })

    xref_df = pl.DataFrame(xref_records).with_columns(
        pl.col("variant").map_elements(
            variant_label, return_dtype=pl.Utf8,
        ).alias("label"),
    )
    print(f"xref_df: {xref_df.shape}")
    return xref_df


def expand_xref(
    mic_df: pl.DataFrame,
    ic50_metric_df: pl.DataFrame,
    pca_df: pl.DataFrame,
    epi_df: pl.DataFrame,
    od_mic_df: pl.DataFrame | None = None,
    auc_ic50_metric_df: pl.DataFrame | None = None,
    od_ic50_metric_df: pl.DataFrame | None = None,
) -> pl.DataFrame:
    EPI_CONCS = {
        "AZT": epi_df.filter(
            pl.col("Drug") == "AZT",
        )["Concentration"].unique().sort().to_list(),
        "AMP": epi_df.filter(
            pl.col("Drug") == "AMP",
        )["Concentration"].unique().sort().to_list(),
    }

    records = []
    for row in mic_df.iter_rows(named=True):
        drug = row["drug"]
        variant = row["variant"]
        genotype = GENO_LOOKUP.get(variant)

        rec = {
            "drug": drug,
            "variant": variant,
            "label": variant_label(variant),
            "genotype_13": genotype,
            "n_mutations": len(variants_meta[variant]["mutations"]),
            "mic_ugml": row["mic_ugml"],
            "mic_is_censored": row.get("is_censored", False),
            "log10_mic": float(np.log10(max(row["mic_ugml"], 1e-3))),
        }

        # IC50 (OD-based from original pipeline)
        ic50_row = ic50_metric_df.filter(
            (pl.col("drug") == drug) & (pl.col("variant") == variant),
        )
        if len(ic50_row) > 0 and ic50_row["ic50_mean"][0] is not None:
            ic50_val = float(ic50_row["ic50_mean"].item())
            rec["ic50_mean"] = ic50_val
            rec["log10_ic50"] = float(np.log10(max(ic50_val, 1e-3)))
        else:
            rec["ic50_mean"] = None
            rec["log10_ic50"] = None

        # OD-based MIC
        if od_mic_df is not None:
            od_row = od_mic_df.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            if len(od_row) > 0:
                od_mic_val = float(od_row["mic_ugml"].item())
                rec["od_mic_ugml"] = od_mic_val
                rec["log10_od_mic"] = float(np.log10(max(od_mic_val, 1e-3)))
            else:
                rec["od_mic_ugml"] = None
                rec["log10_od_mic"] = None
        else:
            rec["od_mic_ugml"] = None
            rec["log10_od_mic"] = None

        # AUC-Hill IC50
        if auc_ic50_metric_df is not None:
            ah_row = auc_ic50_metric_df.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            if len(ah_row) > 0 and ah_row["ic50_mean"][0] is not None:
                ah_val = float(ah_row["ic50_mean"].item())
                rec["auc_hill_ic50"] = ah_val
                rec["log10_auc_hill_ic50"] = float(np.log10(max(ah_val, 1e-3)))
            else:
                rec["auc_hill_ic50"] = None
                rec["log10_auc_hill_ic50"] = None
        else:
            rec["auc_hill_ic50"] = None
            rec["log10_auc_hill_ic50"] = None

        # OD-based IC50
        if od_ic50_metric_df is not None:
            oi_row = od_ic50_metric_df.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            if len(oi_row) > 0 and oi_row["ic50_mean"][0] is not None:
                oi_val = float(oi_row["ic50_mean"].item())
                rec["od_ic50_mean"] = oi_val
                rec["log10_od_ic50"] = float(np.log10(max(oi_val, 1e-3)))
            else:
                rec["od_ic50_mean"] = None
                rec["log10_od_ic50"] = None
        else:
            rec["od_ic50_mean"] = None
            rec["log10_od_ic50"] = None

        # PCA
        if len(pca_df) > 0:
            pca_row = pca_df.filter(
                (pl.col("drug") == drug) & (pl.col("variant") == variant),
            )
            for pc in ["PC1", "PC2", "PC3"]:
                rec[pc] = (
                    float(pca_row[pc].item()) if len(pca_row) > 0 else None
                )
        else:
            for pc in ["PC1", "PC2", "PC3"]:
                rec[pc] = None

        # Fitness at all concentrations
        if genotype is not None:
            fit_rows = epi_df.filter(
                (pl.col("Genotype") == genotype)
                & (pl.col("Drug") == DRUG_MAP[drug]),
            ).sort("Concentration")

            fitness_nonzero = []
            for conc in EPI_CONCS[drug]:
                fit_val = fit_rows.filter(pl.col("Concentration") == conc)
                if len(fit_val) > 0:
                    f = float(fit_val["Fitness"].item())
                    rec[f"fitness_{conc}"] = f
                    if conc > 0:
                        fitness_nonzero.append(f)
                else:
                    rec[f"fitness_{conc}"] = None

            rec["mean_fitness"] = (
                float(np.mean(fitness_nonzero)) if fitness_nonzero else None
            )
            f0 = rec.get("fitness_0.0") or rec.get("fitness_0")
            f_max = rec.get(f"fitness_{EPI_CONCS[drug][-1]}")
            rec["resistance_index"] = (
                float(f_max / f0) if f0 and f_max and f0 != 0 else None
            )
        else:
            for conc in EPI_CONCS[drug]:
                rec[f"fitness_{conc}"] = None
            rec["mean_fitness"] = None
            rec["resistance_index"] = None

        records.append(rec)

    xref_expanded_df = pl.DataFrame(records)
    print(f"xref_expanded_df: {xref_expanded_df.shape}")
    return xref_expanded_df
