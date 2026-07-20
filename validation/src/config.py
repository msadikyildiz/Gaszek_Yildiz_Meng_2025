"""Shared paths, constants, lookups, and Altair theme for MIC analysis."""

import sys
import json
from pathlib import Path

import polars as pl
import numpy as np
import altair as alt

from config_theme import (  # noqa: F401 — re-export for backward compat
    FONT_BODY, FONT_MONO, register_theme,
    build_cluster_col, nudge_labels, make_connectors,
)

# ── Paths ──
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
PROCESSED_DIR = Path(__file__).resolve().parent / "processed"
FIGURES_DIR = ROOT / "figures"
# validation/ lives inside the distribution repo, whose root holds
# data/processed/Epistasis_Combined.parquet. Walk up to find it.
def _epistasis_root(start: Path) -> Path:
    for anc in (start, *start.parents):
        if (anc / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return anc
    return ROOT.parent  # fallback
EPISTASIS_ROOT = _epistasis_root(ROOT)
EPISTASIS_PARQUET = EPISTASIS_ROOT / "data" / "processed" / "Epistasis_Combined.parquet"

# ── Plate-reader utilities ──
sys.path.insert(0, str(ROOT / "plate_reader"))
import utils as utl  # noqa: E402

# ── All experiments (batch-tagged) ──
with open(ROOT / "experiments.json") as _f:
    experiments_raw = json.load(_f)["experiments"]

# ── Batch-varying state (set by set_batch()) ──
variants_meta: dict = {}
VARIANT_NAMES: list[str] = []
GENO_LOOKUP: dict[str, str] = {}
NAME_LOOKUP: dict[str, str] = {}
WT_GENOTYPE: str = ""
BG_EXCLUDE: dict[str, list[str]] = {}
experiments: dict = {}

_NAME_ORDER: list[str] = []
MUTATION_CUES: dict[str, str] = {}
_RAW_COLORS: dict[str, str] = {}
VARIANT_COLORS: dict[str, str] = {}
SHORT_COLORS: dict[str, str] = {}
VARIANT_ORDER: list[str] = []
VARIANT_CLUSTERS: dict[str, list[str]] = {}

CURRENT_BATCH: str = ""


# ── Experiment adaptation ──
def adapt_experiment(exp_raw: dict) -> dict:
    """Convert experiments.json entry to Deniz's utils dict format."""
    drug_abbrev = exp_raw["drug_abbrev"]
    return {
        "file_name": str(DATA_DIR / exp_raw["file_name"]),
        "drug": exp_raw["drug"],
        "largest_conc": exp_raw["largest_conc_ugml"],
        "dilution": exp_raw["dilution_fold"],
        "plates": [[p["variant_1"], p["variant_2"]] for p in exp_raw["plates"]],
        "variants": [
            exp_raw["plate_layout"]["variant_1_rows"],
            exp_raw["plate_layout"]["variant_2_rows"],
        ],
        "background": [
            exp_raw["plate_layout"]["background_rows"],
            list(range(2, 13)),
            BG_EXCLUDE.get(drug_abbrev, []),
        ],
    }


# ── Short-label Vega-Lite expression ──
SHORT_LABEL_EXPR = (
    "indexof(datum.value, '(') >= 0"
    " ? substring(datum.value, indexof(datum.value, '(') + 1,"
    " indexof(datum.value, ')'))"
    " : datum.value"
)


def extract_short_name(label: str) -> str:
    """Pull short name from label like '..K.NT.A.K.T.N (v1.1)' -> 'v1.1'."""
    if "(" in label and label.endswith(")"):
        return label.rsplit("(", 1)[1].rstrip(")")
    return label


def variant_label(name: str) -> str:
    """Dot-notation label: matching -> '.', mismatches keep letter, + (short_name)."""
    if name == "DD":
        return "DD"
    geno = GENO_LOOKUP.get(name)
    if geno is None:
        return name
    if name == "WT":
        return f"{geno} ({name})"
    dots = "".join("." if g == w else g for g, w in zip(geno, WT_GENOTYPE))
    return f"{dots} ({name})"


# ── Constants ──
MIC_THRESHOLD = 0.10
TARGET_HOURS = [6, 12, 18, 24]
END_HOURS = [12, 24]
DRUG_MAP = {"AZT": "AZT", "AMP": "AMP"}
MAX_CONC = {"AZT": 324.0, "AMP": 781.0}

# ── Epistasis data (lazy load) ──
_epi_df_cache = None


def load_epi_df() -> pl.DataFrame:
    """Lazy-load epistasis parquet, filtered to our variants."""
    global _epi_df_cache
    if _epi_df_cache is not None:
        return _epi_df_cache
    target_genotypes = list(GENO_LOOKUP.values())
    _epi_df_cache = (
        pl.scan_parquet(EPISTASIS_PARQUET)
        .filter(pl.col("Genotype").is_in(target_genotypes))
        .select("Genotype", "Epistatic Order", "Drug", "Concentration", "Fitness")
        .collect()
    )
    return _epi_df_cache


def _compute_cue(name: str, meta: dict, no_q39k_base: set, plus_q39k_base: set) -> str:
    """Short string showing how a variant differs from its cluster base."""
    muts = set(meta["mutations"])
    if name == "WT":
        return "wild-type"
    if name == "DD":
        return "S70A + E166A"
    # Colony isolates: anything not WT/DD/v-series
    import re
    if not re.match(r'^v\d', name):
        return f"{len(muts)} mut"
    base = plus_q39k_base if "Q39K" in muts else no_q39k_base
    if muts == base:
        return "base"
    added = sorted(muts - base)
    removed = sorted(base - muts)
    parts = [f"-{m}" for m in removed] + [f"+{m}" for m in added]
    return " ".join(parts)


def make_cluster_legend(
    exclude: set[str] | None = None, cols_per_row: int = 2,
) -> alt.VConcatChart:
    """Cluster-grouped legend with padded tabular text and matched row heights."""
    exclude = exclude or set()
    entries: list[tuple[str, list[str]]] = []
    for cname, members in VARIANT_CLUSTERS.items():
        filtered = [m for m in members if m not in exclude]
        if filtered:
            entries.append((cname, filtered))

    row_chunks: list[alt.HConcatChart] = []
    for i in range(0, len(entries), cols_per_row):
        chunk = entries[i:i + cols_per_row]
        cols = [
            build_cluster_col(
                cname, ms,
                short_colors=SHORT_COLORS,
                geno_lookup=GENO_LOOKUP,
                wt_genotype=WT_GENOTYPE,
                mutation_cues=MUTATION_CUES,
            )
            for cname, ms in chunk
        ]
        row_chunks.append(alt.hconcat(*cols).properties(spacing=16))
    return alt.vconcat(*row_chunks).properties(spacing=6)


# ── Batch management ──

def set_batch(batch_id: str) -> None:
    """Load batch-specific variants, experiments, colors, and output dirs."""
    global variants_meta, VARIANT_NAMES, GENO_LOOKUP, NAME_LOOKUP, WT_GENOTYPE
    global BG_EXCLUDE, experiments, _NAME_ORDER, MUTATION_CUES
    global _RAW_COLORS, VARIANT_COLORS, SHORT_COLORS, VARIANT_ORDER, VARIANT_CLUSTERS
    global PROCESSED_DIR, FIGURES_DIR, CURRENT_BATCH, _epi_df_cache

    CURRENT_BATCH = batch_id

    # Load batch-specific variants.json
    batch_variants_path = DATA_DIR / batch_id / "variants.json"
    with open(batch_variants_path) as f:
        batch_data = json.load(f)

    variants_meta = batch_data["variants"]
    VARIANT_NAMES = list(variants_meta.keys())

    GENO_LOOKUP = {
        name: v["genotype_13"]
        for name, v in variants_meta.items()
        if v["genotype_13"] is not None
    }
    NAME_LOOKUP = {v: k for k, v in GENO_LOOKUP.items()}
    WT_GENOTYPE = variants_meta["WT"]["genotype_13"]

    # Display metadata from batch variants.json
    _NAME_ORDER = batch_data["display_order"]
    _RAW_COLORS = batch_data["colors"]
    VARIANT_CLUSTERS = batch_data["clusters"]
    BG_EXCLUDE = batch_data.get("bg_exclude", {})

    # Mutation cues
    no_q39k_base_names = [n for n in _NAME_ORDER if n not in ("WT", "DD") and not n.startswith("Col") and "Q39K" not in variants_meta.get(n, {}).get("mutations", [])]
    plus_q39k_base_names = [n for n in _NAME_ORDER if n not in ("WT", "DD") and not n.startswith("Col") and "Q39K" in variants_meta.get(n, {}).get("mutations", [])]
    no_q39k_base = set(variants_meta[no_q39k_base_names[0]]["mutations"]) if no_q39k_base_names else set()
    plus_q39k_base = set(variants_meta[plus_q39k_base_names[0]]["mutations"]) if plus_q39k_base_names else set()
    MUTATION_CUES = {
        name: _compute_cue(name, variants_meta[name], no_q39k_base, plus_q39k_base)
        for name in _NAME_ORDER
        if name in variants_meta
    }

    VARIANT_COLORS = {variant_label(k): v for k, v in _RAW_COLORS.items()}
    SHORT_COLORS = dict(_RAW_COLORS)
    VARIANT_ORDER = [variant_label(v) for v in _NAME_ORDER]

    # Filter experiments to this batch
    batch_exps = [e for e in experiments_raw if e["batch"] == batch_id]
    experiments = {e["drug_abbrev"]: adapt_experiment(e) for e in batch_exps}

    # Output dirs scoped to batch
    PROCESSED_DIR = Path(__file__).resolve().parent / "processed" / batch_id
    FIGURES_DIR = ROOT / "figures" / batch_id

    # Invalidate epistasis cache (different variants may be in play)
    _epi_df_cache = None


# ── Helpers ──
def load_processed(name: str, suffix: str = "") -> pl.DataFrame:
    """Read a parquet from the processed cache."""
    return pl.read_parquet(PROCESSED_DIR / f"{name}{suffix}.parquet")


# ── Auto-detect batch from CLI at import time ──
_default_batch = "20260220"
if "--batch" in sys.argv:
    _default_batch = sys.argv[sys.argv.index("--batch") + 1]
set_batch(_default_batch)
