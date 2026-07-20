"""Run full MIC analysis pipeline: process -> figures -> export.

Supports multi-timepoint batch: generates suffixed parquets + PNGs
for each END_HOUR, plus the unsuffixed (full-length) baseline.

Usage:
    python run_all.py                  # run all batches
    python run_all.py --batch 20260220 # run one batch
"""

import subprocess
import sys
from pathlib import Path

from config import END_HOURS, experiments_raw

FIGURE_SCRIPTS = [
    "fig_a_dose_response.py",
    "fig_b1_metric_correlation.py",
    "fig_b_mic_dotplot.py",
    "fig_b2_ic50_dotplot.py",
    "fig_c_mic_vs_fitness.py",
    "fig_e_correlation_heatmap.py",
    "fig_f_best_pairs.py",
    "fig_h_library_fitness.py",
    ["-m", "fig_mic_ic50_fits"],
]


def _run(args: list[str], label: str, cwd: str) -> bool:
    print(f"\n{'='*60}")
    print(f"Running {label}")
    print(f"{'='*60}")
    result = subprocess.run([sys.executable, *args], cwd=cwd)
    if result.returncode != 0:
        print(f"FAILED: {label} (exit code {result.returncode})")
        return False
    return True


def main():
    src_dir = str(Path(__file__).resolve().parent)

    # Determine which batches to run
    all_batches = sorted(set(e["batch"] for e in experiments_raw))
    if "--batch" in sys.argv:
        idx = sys.argv.index("--batch")
        batches_to_run = [sys.argv[idx + 1]]
    else:
        batches_to_run = all_batches

    failures = []

    for batch_id in batches_to_run:
        batch_flag = ["--batch", batch_id]
        print(f"\n{'#'*60}")
        print(f"# BATCH: {batch_id}")
        print(f"{'#'*60}")

        # 1. Full-length (unsuffixed) pipeline
        if not _run(["-m", "process_data"] + batch_flag, f"process_data (full) [{batch_id}]", src_dir):
            failures.append(f"process_data (full) [{batch_id}]")
        for entry in FIGURE_SCRIPTS:
            args = [entry] if isinstance(entry, str) else list(entry)
            label = entry if isinstance(entry, str) else " ".join(entry)
            if not _run(args + batch_flag, f"{label} [{batch_id}]", src_dir):
                failures.append(f"{label} [{batch_id}]")
        _run(["export_csv.py"] + batch_flag, f"export_csv.py [{batch_id}]", src_dir)

        # 2. Per-timepoint suffixed runs
        for hour in END_HOURS:
            suffix = f"_{int(hour)}h"
            plabel = f"process_data {suffix} [{batch_id}]"
            if not _run(["-m", "process_data", "--end-hour", str(hour)] + batch_flag, plabel, src_dir):
                failures.append(plabel)
                continue
            for entry in FIGURE_SCRIPTS:
                args = [entry] if isinstance(entry, str) else list(entry)
                label = entry if isinstance(entry, str) else " ".join(entry)
                tag = f"{label} {suffix} [{batch_id}]"
                if not _run(args + ["--suffix", suffix] + batch_flag, tag, src_dir):
                    failures.append(tag)

    steps_per_batch = 1 + len(FIGURE_SCRIPTS) + len(END_HOURS) * (1 + len(FIGURE_SCRIPTS))
    total = len(batches_to_run) * steps_per_batch
    print(f"\n{'='*60}")
    if failures:
        print(f"FAILURES ({len(failures)}/{total}): {failures}")
        sys.exit(1)
    else:
        print(f"All {total} steps across {len(batches_to_run)} batch(es) completed successfully")


if __name__ == "__main__":
    main()
