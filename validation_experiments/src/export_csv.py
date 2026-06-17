"""Export summary tables as CSV."""

from config import load_processed, PROCESSED_DIR


def main():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    try:
        xref_df = load_processed("xref_df")
        out = PROCESSED_DIR / "mic_summary.csv"
        xref_df.write_csv(out)
        print(f"Saved: {out}")
    except FileNotFoundError:
        print("  xref_df not found (no epistasis match) — skipping mic_summary.csv")

    try:
        xref_expanded_df = load_processed("xref_expanded_df")
        out2 = PROCESSED_DIR / "mic_expanded_summary.csv"
        xref_expanded_df.write_csv(out2)
        print(f"Saved: {out2}")
    except FileNotFoundError:
        print("  xref_expanded_df not found — skipping mic_expanded_summary.csv")

    try:
        corr_results_df = load_processed("corr_results_df")
        out3 = PROCESSED_DIR / "metric_correlations.csv"
        corr_results_df.write_csv(out3)
        print(f"Saved: {out3}")
    except FileNotFoundError:
        print("  corr_results_df not found — skipping metric_correlations.csv")


if __name__ == "__main__":
    main()
