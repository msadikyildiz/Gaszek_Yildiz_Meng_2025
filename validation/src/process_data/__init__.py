"""Process raw plates -> intermediate parquets for all downstream figures."""

import polars as pl

from config import load_epi_df, PROCESSED_DIR, GENO_LOOKUP

from .plates import process_plates, compute_mic
from .metrics import compute_ic50, compute_ic50_from_long
from .xref import cross_reference, expand_xref
from .pca import compute_pca
from .cell_counts import compute_cell_counts
from .correlations import compute_correlations, cache_growth_curves


def _has_epistasis_variants() -> bool:
    """Check if any batch variants exist in the epistasis dataset."""
    target_genotypes = list(GENO_LOOKUP.values())
    if not target_genotypes:
        return False
    try:
        epi_genos = set(
            pl.scan_parquet(
                str(PROCESSED_DIR).replace(
                    str(PROCESSED_DIR), ""
                ) or "dummy"
            ).select("Genotype").unique().collect()["Genotype"].to_list()
        )
    except Exception:
        pass
    # Simpler: just check if EPISTASIS_PARQUET has our genotypes
    from config import EPISTASIS_PARQUET
    try:
        hits = (
            pl.scan_parquet(EPISTASIS_PARQUET)
            .filter(pl.col("Genotype").is_in(target_genotypes))
            .select("Genotype")
            .head(1)
            .collect()
        )
        return len(hits) > 0
    except Exception:
        return False


def main(end_hour: float | None = None):
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    suffix = f"_{int(end_hour)}h" if end_hour is not None else ""

    long_df = process_plates(end_hour=end_hour)
    long_df.write_parquet(PROCESSED_DIR / f"long_df{suffix}.parquet")

    # AUC-based MIC (existing)
    mic_df, dr_df, mic_rep_df = compute_mic(long_df, "final_auc")
    mic_df.write_parquet(PROCESSED_DIR / f"mic_df{suffix}.parquet")
    dr_df.write_parquet(PROCESSED_DIR / f"dr_df{suffix}.parquet")
    mic_rep_df.write_parquet(PROCESSED_DIR / f"mic_rep_df{suffix}.parquet")

    # OD-based MIC (new)
    od_mic_df, od_dr_df, od_mic_rep_df = compute_mic(long_df, "final_od")
    od_mic_df.write_parquet(PROCESSED_DIR / f"od_mic_df{suffix}.parquet")
    od_dr_df.write_parquet(PROCESSED_DIR / f"od_dr_df{suffix}.parquet")
    od_mic_rep_df.write_parquet(PROCESSED_DIR / f"od_mic_rep_df{suffix}.parquet")

    # OD-snapshot IC50 (existing pipeline)
    ic50_metric_df, ic50_rep_df = compute_ic50(end_hour=end_hour)
    ic50_metric_df.write_parquet(PROCESSED_DIR / f"ic50_metric_df{suffix}.parquet")
    ic50_rep_df.write_parquet(PROCESSED_DIR / f"ic50_rep_df{suffix}.parquet")

    # AUC-Hill IC50 (new)
    auc_ic50_metric_df, auc_ic50_rep_df = compute_ic50_from_long(long_df, "final_auc")
    auc_ic50_metric_df.write_parquet(PROCESSED_DIR / f"auc_ic50_metric_df{suffix}.parquet")
    auc_ic50_rep_df.write_parquet(PROCESSED_DIR / f"auc_ic50_rep_df{suffix}.parquet")

    # OD-Hill IC50 (new)
    od_ic50_metric_df, od_ic50_rep_df = compute_ic50_from_long(long_df, "final_od")
    od_ic50_metric_df.write_parquet(PROCESSED_DIR / f"od_ic50_metric_df{suffix}.parquet")
    od_ic50_rep_df.write_parquet(PROCESSED_DIR / f"od_ic50_rep_df{suffix}.parquet")

    # Epistasis cross-ref, PCA, cell counts — only if variants match
    has_epi = _has_epistasis_variants()

    if has_epi:
        epi_df = load_epi_df()

        xref_df = cross_reference(mic_df, epi_df)
        xref_df.write_parquet(PROCESSED_DIR / f"xref_df{suffix}.parquet")

        # PCA: computed once from read counts, independent of end_hour
        pca_path = PROCESSED_DIR / "pca_df.parquet"
        if pca_path.exists() and suffix:
            pca_df = pl.read_parquet(pca_path)
        else:
            pca_df = compute_pca()
            pca_df.write_parquet(pca_path)

        # Cell counts: computed once from read counts, independent of end_hour
        cc_path = PROCESSED_DIR / "cell_counts_df.parquet"
        if cc_path.exists() and suffix:
            pass  # reuse cached
        else:
            cell_counts_df = compute_cell_counts()
            cell_counts_df.write_parquet(cc_path)

        xref_expanded_df = expand_xref(
            mic_df, ic50_metric_df, pca_df, epi_df,
            od_mic_df=od_mic_df,
            auc_ic50_metric_df=auc_ic50_metric_df,
            od_ic50_metric_df=od_ic50_metric_df,
        )
        xref_expanded_df.write_parquet(PROCESSED_DIR / f"xref_expanded_df{suffix}.parquet")

        corr_results_df = compute_correlations(xref_expanded_df)
        corr_results_df.write_parquet(PROCESSED_DIR / f"corr_results_df{suffix}.parquet")
    else:
        print("  Skipping epistasis cross-ref, PCA, cell counts (no matching variants)")

    gc_df = cache_growth_curves(end_hour=end_hour)
    gc_df.write_parquet(PROCESSED_DIR / f"growth_curves{suffix}.parquet")

    print(f"\nAll parquets written to: {PROCESSED_DIR} (suffix={suffix!r})")
