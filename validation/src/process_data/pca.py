"""PCA from sequencing read counts (cell 9)."""

import polars as pl
import numpy as np

from config import GENO_LOOKUP, WT_GENOTYPE, EPISTASIS_ROOT


def compute_pca() -> pl.DataFrame:
    meta_df = pl.read_csv(EPISTASIS_ROOT / "data" / "raw" / "metadata.csv")
    meta_df = meta_df.with_columns(
        pl.col("Concentration").cast(pl.Utf8).str.replace_all(",", "")
        .cast(pl.Float64).alias("Concentration_f"),
    )
    meta_df = meta_df.filter(
        pl.col("Drug").is_in(["Aztreonam", "Ampicillin"])
        & (pl.col("Timepoint") != "0h")
        & ~pl.col("Sample Name").str.contains("3_1mg"),
    ).with_columns(
        pl.when(pl.col("Drug") == "Aztreonam").then(pl.lit("AZT"))
        .when(pl.col("Drug") == "Ampicillin").then(pl.lit("AMP"))
        .alias("drug_abbrev"),
        pl.col("Timepoint").str.replace("h", "").cast(pl.Float64).alias("timepoint_h"),
    )

    geno_to_dot = {}
    for name, geno in GENO_LOOKUP.items():
        dots = "".join(
            "." if g == w else g for g, w in zip(geno, WT_GENOTYPE)
        )
        geno_to_dot[name] = dots
    dot_to_name = {v: k for k, v in geno_to_dot.items()}
    target_dots = list(geno_to_dot.values())

    drug_file_map = {
        "AZT": "Aztreonam_read_counts_per_genotype.csv",
        "AMP": "Ampicillin_read_counts_per_genotype.csv",
    }

    pca_records = {}
    for drug_abbrev, fname in drug_file_map.items():
        rc_df = pl.read_csv(EPISTASIS_ROOT / "data" / "raw" / fname)
        rc_sub = rc_df.filter(
            pl.col("mut_profile_masked").is_in(target_dots),
        )

        drug_meta = meta_df.filter(pl.col("drug_abbrev") == drug_abbrev)
        ut_meta = drug_meta.filter(pl.col("Concentration_f") == 0.0)
        treated_meta = drug_meta.filter(pl.col("Concentration_f") > 0.0)

        # CPM normalisation: compute per-sample total reads from ALL genotypes
        # (not just target 8) so library-size estimate is robust.
        sample_cols = [c for c in rc_df.columns if c != "mut_profile_masked"]
        sample_totals = {}
        for sc in sample_cols:
            total = rc_df[sc].sum()
            sample_totals[sc] = total if total > 0 else 1.0

        ut_by_tp = {}
        for tp in ut_meta["timepoint_h"].unique().to_list():
            tp_samples = ut_meta.filter(
                pl.col("timepoint_h") == tp,
            )["Sample Name"].to_list()
            valid = [s for s in tp_samples if s in rc_sub.columns]
            if valid:
                # CPM-normalise each UT sample before averaging
                cpm_exprs = [
                    (pl.col(s) / sample_totals[s] * 1e6).alias(s)
                    for s in valid
                ]
                ut_by_tp[tp] = rc_sub.select(
                    "mut_profile_masked", *valid,
                ).with_columns(cpm_exprs).with_columns(
                    pl.mean_horizontal(*valid).alias(f"ut_{tp}h"),
                ).select("mut_profile_masked", f"ut_{tp}h")

        feature_cols = []
        treated_sorted = list(
            treated_meta.sort(
                "Concentration_f", "timepoint_h", "Replicate",
            ).iter_rows(named=True)
        )
        for row in treated_sorted:
            sample = row["Sample Name"]
            tp = row["timepoint_h"]
            if sample in rc_sub.columns and tp in ut_by_tp:
                feature_cols.append(sample)

        if not feature_cols:
            continue

        genotypes = rc_sub["mut_profile_masked"].to_list()
        # CPM-normalise treated samples
        raw_matrix = rc_sub.select(feature_cols).to_numpy().astype(float)
        for j_col, col_name in enumerate(feature_cols):
            raw_matrix[:, j_col] = raw_matrix[:, j_col] / sample_totals[col_name] * 1e6

        ut_matrix = np.zeros_like(raw_matrix)
        j = 0
        for row in treated_sorted:
            sample = row["Sample Name"]
            tp = row["timepoint_h"]
            if sample not in rc_sub.columns or tp not in ut_by_tp:
                continue
            ut_col = ut_by_tp[tp].select(
                f"ut_{tp}h",
            ).to_numpy().flatten().astype(float)
            ut_matrix[:, j] = ut_col
            j += 1

        pseudocount = 1.0
        lfc = np.log2(
            (raw_matrix + pseudocount) / (ut_matrix + pseudocount),
        )

        ct_keys = []
        for row in treated_sorted:
            if row["Sample Name"] in rc_sub.columns and row["timepoint_h"] in ut_by_tp:
                ct_keys.append((row["Concentration_f"], row["timepoint_h"]))

        unique_ct = sorted(set(ct_keys))
        avg_lfc = np.zeros((len(genotypes), len(unique_ct)))
        for k, (c, t) in enumerate(unique_ct):
            indices = [i for i, ct in enumerate(ct_keys) if ct == (c, t)]
            avg_lfc[:, k] = lfc[:, indices].mean(axis=1)

        col_mean = avg_lfc.mean(axis=0)
        col_std = avg_lfc.std(axis=0)
        col_std[col_std == 0] = 1.0
        Z = (avg_lfc - col_mean) / col_std

        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        variance_explained = (S ** 2) / (S ** 2).sum()
        scores = U * S

        # Sign-correct PC1: positive = enriched (higher fitness under drug)
        # Correlate PC1 with row-mean LFC; flip if anti-correlated.
        mean_lfc = avg_lfc.mean(axis=1)
        if np.corrcoef(scores[:, 0], mean_lfc)[0, 1] < 0:
            scores[:, 0] *= -1

        pca_records[drug_abbrev] = pl.DataFrame({
            "variant": [dot_to_name.get(g, g) for g in genotypes],
            "drug": [drug_abbrev] * len(genotypes),
            "PC1": scores[:, 0].tolist(),
            "PC2": scores[:, 1].tolist(),
            "PC3": scores[:, 2].tolist(),
            "var_PC1": [float(variance_explained[0])] * len(genotypes),
            "var_PC2": [float(variance_explained[1])] * len(genotypes),
            "var_PC3": [float(variance_explained[2])] * len(genotypes),
        })

    pca_df = (
        pl.concat(list(pca_records.values()))
        if pca_records
        else pl.DataFrame()
    )
    print(f"pca_df: {pca_df.shape}")
    return pca_df
