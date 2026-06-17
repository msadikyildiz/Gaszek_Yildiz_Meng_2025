"""Portable data loaders for the graph-analysis notebooks in this directory.

The notebooks were originally written against a cluster checkout that imported
`src.preprocess` / `src.pair_table_global`. These functions reproduce the inputs
those modules provided, backed entirely by this repository's canonical processed
data so the notebooks run directly from the repo.

- ``preprocess_data`` returns the per-drug long/wide AUC tables from
  ``data/processed/{amp,azt}_auc_long_df.parquet``.
- ``calculate_normalized_fitness`` returns per-genotype "global fitness",
  defined as the local median fitness at the representative concentration used
  throughout the paper (AMP 781 ug/ml, AZT 36 ug/ml). At that concentration the
  local median equals the ``Fitness`` column of ``Epistasis_Combined``
  (verified identical), and unlike the epistasis slice it retains the
  dead-mutant baseline that several plots reference.
- ``load_global_epistasis`` returns the per-drug global epistasis tables,
  sliced from ``Epistasis_Combined`` exactly as 05_epistasis_figures.ipynb
  builds ``Epistasis_Combined_{AMP,AZT}_auc_10``.
"""

from pathlib import Path

import polars as pl

# amp_select_conc / azt_select_conc from 05_epistasis_figures.ipynb
AMP_SELECT_CONC = 781.0
AZT_SELECT_CONC = 36.0


def find_repo_root(start=None):
    here = Path(start or __file__).resolve()
    for parent in (here, *here.parents):
        if (parent / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return parent
    raise FileNotFoundError(
        "could not locate repo root (data/processed/Epistasis_Combined.parquet)"
    )


REPO_ROOT = find_repo_root()
REPO_DATA = REPO_ROOT / "data"
PROC = REPO_DATA / "processed"


def _patch_joypy():
    """joypy 0.2.6 indexes the result of pandas' ``flatten_axes``; modern pandas
    returns a generator, which breaks ``_axes[i]`` in ``joypy.joyplot``. Wrap it
    to materialize a list so the ridgeline/joyplots render. No-op if joypy or the
    helper is absent. Imported by the notebooks before any joyplot call.
    """
    try:
        import functools
        import importlib
        # joypy/__init__ rebinds the name `joypy.joyplot` to the function, so the
        # submodule must be fetched via importlib to reach its `_flatten`.
        _jj = importlib.import_module("joypy.joyplot")
    except Exception:
        return
    orig = getattr(_jj, "_flatten", None)
    if orig is None or getattr(orig, "_coerced_to_list", False):
        return

    @functools.wraps(orig)
    def _flatten_to_list(*a, **k):
        return list(orig(*a, **k))

    _flatten_to_list._coerced_to_list = True
    _jj._flatten = _flatten_to_list


_patch_joypy()


def preprocess_data(amp_path=None, azt_path=None, amp_pairs_path=None,
                    azt_pairs_path=None, clean_nulls_flag=True):
    """Return ``{drug: {original, long, pairs}}`` from the repo's canonical tables.

    Path arguments are accepted for call-signature compatibility and ignored;
    data is always read from ``data/processed/``. ``pairs`` is unused by the
    notebooks and returned as ``None``.
    """
    def load(drug):
        long = pl.read_parquet(PROC / f"{drug}_auc_long_df.parquet")
        wide = pl.read_parquet(PROC / f"{drug}_auc_wide_df.parquet")
        return {"original": wide, "long": long, "pairs": None}

    return {"amp": load("amp"), "azt": load("azt")}


def calculate_normalized_fitness(long_df):
    """Per-genotype global fitness = local median fitness at the representative
    concentration (AMP 781 / AZT 36), inferred from the drug's concentration set.
    """
    concs = set(long_df["concentration"].unique().to_list())
    select_conc = AMP_SELECT_CONC if AMP_SELECT_CONC in concs else AZT_SELECT_CONC
    sub = long_df.filter(pl.col("concentration") == select_conc)
    if "median" in sub.columns:
        fitness = pl.col("median")
    else:
        fitness = pl.concat_list("replicate1", "replicate2", "replicate3").list.median()
    return sub.select("mutant_profile", fitness.alias("normalized_fitness"))


def load_global_epistasis():
    """Return ``(combined, amp_slice, azt_slice)`` global epistasis tables,
    sliced from ``Epistasis_Combined`` at the representative concentration.
    """
    ep = pl.read_parquet(PROC / "Epistasis_Combined.parquet")
    amp = ep.filter((pl.col("Drug") == "AMP") & (pl.col("Concentration") == AMP_SELECT_CONC))
    azt = ep.filter((pl.col("Drug") == "AZT") & (pl.col("Concentration") == AZT_SELECT_CONC))
    return ep, amp, azt
