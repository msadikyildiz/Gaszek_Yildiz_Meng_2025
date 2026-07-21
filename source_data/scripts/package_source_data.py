#!/usr/bin/env python3
"""Package the release Source Data archive: source_data.xlsx plus every bulk
per-figure CSV, zipped as Source_Data.zip.

The workbook (build_source_data.py) embeds the small plotted tables directly and
points, in each sheet's note, to the large per-genotype companion tables. This
script gathers the workbook and all of those CSVs — from the three locations
that hold them — into one archive:

  source_data/derived/                                  Fig 1-5, S7, S8, S12 (+ Fig2 nodes/edges)
  src/evolutionary_statistics/source_data/{,full/}      Fig 6, S6 (DCA)
  src/graph_analysis/s18_peak_robustness/source_data/   S18 (graph)

Run the builders and the workbook assembler first (see source_data/README.md),
then this script.
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path


def _repo_root(p: Path) -> Path:
    for a in (p, *p.parents):
        if (a / "data" / "processed" / "Epistasis_Combined.parquet").exists():
            return a
    raise FileNotFoundError("repo root not found from " + str(p))


REPO = _repo_root(Path(__file__).resolve())
XLSX = REPO / "source_data" / "source_data.xlsx"
OUT = REPO / "source_data" / "Source_Data.zip"
CSV_DIRS = [
    REPO / "source_data" / "derived",
    REPO / "src" / "evolutionary_statistics" / "source_data",
    REPO / "src" / "evolutionary_statistics" / "source_data" / "full",
    REPO / "src" / "graph_analysis" / "s18_peak_robustness" / "source_data",
]


def main() -> None:
    if not XLSX.exists():
        sys.exit(f"missing {XLSX} — run source_data/scripts/build_source_data.py first")
    csvs: dict[str, Path] = {}
    for d in CSV_DIRS:
        for p in sorted(d.glob("*.csv")) if d.exists() else []:
            if p.name in csvs:
                sys.exit(f"duplicate table name across source dirs: {p.name}")
            csvs[p.name] = p
    if not csvs:
        sys.exit("no CSV tables found — run the source-data builders first")

    OUT.unlink(missing_ok=True)
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        z.write(XLSX, "source_data.xlsx")
        for name in sorted(csvs):
            z.write(csvs[name], f"tables/{name}")
    size_mb = OUT.stat().st_size / 1e6
    print(f"wrote {OUT}  ({len(csvs)} tables + workbook, {size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
