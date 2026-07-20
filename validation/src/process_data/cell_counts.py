"""Estimated cell counts from pooled library sequencing.

Computes  estimated_cell_count = OD600 × (variant_reads / total_sample_reads)
for each genotyped variant across drugs, concentrations, timepoints, replicates.
"""

import polars as pl

from config import GENO_LOOKUP, WT_GENOTYPE, EPISTASIS_ROOT, variant_label


def compute_cell_counts() -> pl.DataFrame:
    # ── Load bulk OD600 metadata ──
    meta_df = pl.read_csv(EPISTASIS_ROOT / "data" / "raw" / "metadata.csv")
    meta_df = meta_df.with_columns(
        pl.col("Concentration").cast(pl.Utf8).str.replace_all(",", "")
        .cast(pl.Float64).alias("conc_f"),
    ).filter(
        pl.col("Drug").is_in(["Aztreonam", "Ampicillin"]),
    ).with_columns(
        pl.when(pl.col("Drug") == "Aztreonam").then(pl.lit("AZT"))
        .when(pl.col("Drug") == "Ampicillin").then(pl.lit("AMP"))
        .alias("drug"),
        pl.col("Timepoint").str.replace("h", "").cast(pl.Float64).alias("timepoint_h"),
    )

    # ── Dot-notation lookup (same logic as pca.py) ──
    geno_to_dot: dict[str, str] = {}
    for name, geno in GENO_LOOKUP.items():
        dots = "".join("." if g == w else g for g, w in zip(geno, WT_GENOTYPE))
        geno_to_dot[name] = dots
    dot_to_name = {v: k for k, v in geno_to_dot.items()}
    target_dots = list(geno_to_dot.values())

    drug_file_map = {
        "AZT": "Aztreonam_read_counts_per_genotype.csv",
        "AMP": "Ampicillin_read_counts_per_genotype.csv",
    }

    records: list[dict] = []
    for drug_abbrev, fname in drug_file_map.items():
        rc_df = pl.read_csv(EPISTASIS_ROOT / "data" / "raw" / fname)
        rc_sub = rc_df.filter(
            pl.col("mut_profile_masked").is_in(target_dots),
        )

        drug_meta = meta_df.filter(pl.col("drug") == drug_abbrev)

        # Each sample column = one (drug, conc, timepoint, replicate) combo
        sample_cols = [c for c in rc_sub.columns if c != "mut_profile_masked"]
        meta_lookup = {
            row["Sample Name"]: row
            for row in drug_meta.iter_rows(named=True)
        }

        for sample in sample_cols:
            if sample not in meta_lookup:
                continue
            row = meta_lookup[sample]
            od600 = row["OD600"]
            conc = row["conc_f"]
            tp_h = row["timepoint_h"]
            replicate = row["Replicate"]

            # Total reads across ALL genotypes in this sample
            total_reads = float(rc_df[sample].sum())
            if total_reads <= 0:
                continue

            # Per-variant reads & cell count
            for dot, count_val in zip(
                rc_sub["mut_profile_masked"].to_list(),
                rc_sub[sample].to_list(),
            ):
                variant_name = dot_to_name.get(dot)
                if variant_name is None:
                    continue
                read_frac = float(count_val) / total_reads
                records.append({
                    "drug": drug_abbrev,
                    "variant": variant_name,
                    "label": variant_label(variant_name),
                    "concentration": conc,
                    "timepoint_h": tp_h,
                    "replicate": replicate,
                    "od600": od600,
                    "read_fraction": read_frac,
                    "estimated_cell_count": od600 * read_frac,
                })

    cc_df = pl.DataFrame(records)
    print(f"cell_counts_df: {cc_df.shape}")
    return cc_df
