#!/usr/bin/env python3
"""Parallel graph builder for global peak analysis.

Builds Aztreonam graphs for specific concentrations and neutral thresholds
using multiprocessing to leverage all available CPU cores.

This script is optimized for the global peak robustness analysis, building only
the graphs needed:
- Aztreonam only (not Ampicillin)
- 4 specific concentrations: 12.0, 36.0, 108.0, 324.0
- 31 thresholds: 0.15 to 0.45 (step 0.01)
- Total: 124 graphs

Performance: ~15-20 minutes on 8 cores vs. 60-90 minutes for sequential approach.
With --save-images: ~20-25 minutes (adds ~1-2 seconds per graph for image export).
"""

import argparse
import os
import sys
from multiprocessing import Pool, cpu_count

import networkx as nx

from fitness_landscape_graph.build_graph import convert_graph_attributes_to_json
from fitness_landscape_graph.graph_builder import GraphBuilder
from fitness_landscape_graph.preprocess import preprocess_data


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Build Aztreonam graphs in parallel for global peak analysis"
    )
    parser.add_argument(
        "--base-path",
        type=str,
        required=True,
        help="Base path to project root directory",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save the output graph files",
    )
    parser.add_argument(
        "--neutral-thresholds",
        type=float,
        nargs=3,
        metavar=("START", "END", "STEP"),
        default=[0.15, 0.45, 0.01],
        help="Neutral threshold range (start, end, step). Default: 0.15 0.45 0.01",
    )
    parser.add_argument(
        "--concentrations",
        type=float,
        nargs="+",
        default=[12.0, 36.0, 108.0, 324.0],
        help="Concentrations to build graphs for. Default: 12.0 36.0 108.0 324.0",
    )
    parser.add_argument(
        "--tiny-initial-threshold",
        type=float,
        default=0.02,
        help="Threshold for merging tiny initial edges. Default: 0.02",
    )
    parser.add_argument(
        "--large-edge-threshold",
        type=float,
        default=5.5,
        help="Threshold for identifying potential forbidden pairs. Default: 5.5",
    )
    parser.add_argument(
        "--num-forbidden-pairs",
        type=int,
        default=1,
        help="Number of forbidden pairs to select. Default: 1",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: min(cpu_count, 8))",
    )
    parser.add_argument(
        "--save-images",
        action="store_true",
        default=False,
        help="Save visualization images (PNG) alongside GraphML files",
    )
    parser.add_argument(
        "--image-format",
        type=str,
        default="png",
        choices=["png", "pdf", "svg"],
        help="Image format for visualizations. Default: png",
    )
    return parser.parse_args()


def build_single_graph(args):
    """Build one graph for a specific threshold/concentration and optionally save image.

    Args:
        args: Tuple of (concentration, threshold, azt_long_df, azt_pairs_df,
                       output_path, params, save_image, image_path, image_format)

    Returns:
        Tuple of (concentration, threshold, success, info)
        - success: True if graph built successfully, False otherwise
        - info: output_path if success=True, error message if success=False
    """
    (
        concentration,
        threshold,
        azt_long_df,
        azt_pairs_df,
        output_path,
        params,
        save_image,
        image_path,
        image_format,
    ) = args

    try:
        # Create GraphBuilder
        builder = GraphBuilder(long_table=azt_long_df, pairwise_table=azt_pairs_df)

        # Build graph
        graph = builder.build_graph(
            concentration=concentration,
            neutral_threshold=threshold,
            fitness_col="median",
            fitness_diff_col="median_diff",
            use_parallel_peak_merge=True,  # Internal parallelism for peak merging
            rename_by_avg_fitness=True,
            merge_peaks=False,
            **params,  # large_edge_threshold, tiny_initial_threshold, etc.
        )

        # Convert dict attributes to JSON and save
        graph = convert_graph_attributes_to_json(graph)
        nx.write_graphml(graph, output_path)

        # Save visualization image if requested
        if save_image and image_path:
            from fitness_landscape_graph.graph_analyzer import GraphAnalyzer, VisConfig

            # Create analyzer and generate plot
            analyzer = GraphAnalyzer(output_path)  # Load from saved GraphML
            config = VisConfig(
                figure_width=1200,
                figure_height=800,
            )
            fig = analyzer.plot(config=config)

            # Save image (requires kaleido package)
            fig.write_image(image_path, format=image_format, scale=2)

        return (concentration, threshold, True, output_path)

    except Exception as e:
        return (concentration, threshold, False, str(e))


def main():
    """Main execution function."""
    args = parse_args()

    # Print configuration
    print("=" * 80)
    print("Parallel Graph Builder for Global Peak Analysis")
    print("=" * 80)
    print(f"Base path: {args.base_path}")
    print(f"Output directory: {args.output_dir}")
    print(
        f"Neutral thresholds: {args.neutral_thresholds[0]} to "
        f"{args.neutral_thresholds[1]} (step {args.neutral_thresholds[2]})"
    )
    print(f"Concentrations: {args.concentrations}")
    print(f"Tiny initial threshold: {args.tiny_initial_threshold}")
    print(f"Large edge threshold: {args.large_edge_threshold}")
    print(f"Num forbidden pairs: {args.num_forbidden_pairs}")
    print(f"Save images: {args.save_images}")
    if args.save_images:
        print(f"Image format: {args.image_format}")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate threshold list
    start, end, step = args.neutral_thresholds
    # Use round to avoid floating point precision issues
    num_steps = int(round((end - start) / step)) + 1
    thresholds = [round(start + i * step, 2) for i in range(num_steps)]
    print(
        f"Generated {len(thresholds)} thresholds: {thresholds[:3]}...{thresholds[-3:]}"
    )

    # Calculate total graphs
    total_graphs = len(thresholds) * len(args.concentrations)
    print(f"Total graphs to build: {total_graphs}")

    # Determine number of workers
    if args.num_workers is None:
        num_workers = min(cpu_count(), 8)
    else:
        num_workers = args.num_workers
    print(f"Using {num_workers} parallel workers")
    print("=" * 80)

    # File paths
    base_path = args.base_path
    amp_path = f"{base_path}/data/raw/combined-auc/genotype_auc_sorted_ampicillin.csv"
    azt_path = f"{base_path}/data/raw/combined-auc/genotype_auc_sorted_aztreonam.csv"
    amp_pairs_path = f"{base_path}/data/processed/amp_pairs.csv"
    azt_pairs_path = f"{base_path}/data/processed/azt_pairs.csv"

    # Load and preprocess data (once, in main process)
    print("\nLoading and preprocessing Aztreonam data...")
    processed_data = preprocess_data(
        amp_path,
        azt_path,
        amp_pairs_path,
        azt_pairs_path,
        clean_nulls_flag=True,
    )

    # Extract only Aztreonam data
    azt_long_df = processed_data["azt"]["long"]
    azt_pairs_df = processed_data["azt"]["pairs"]
    print(f"Loaded {len(azt_long_df)} rows from long table")
    print(f"Loaded {len(azt_pairs_df)} pairs from pairs table")

    # Fixed parameters (same for all graphs)
    fixed_params = {
        "tiny_initial_threshold": args.tiny_initial_threshold,
        "large_edge_threshold": args.large_edge_threshold,
        "num_forbidden_pairs": args.num_forbidden_pairs,
    }

    # Create images subfolder if needed
    images_dir = None
    if args.save_images:
        images_dir = os.path.join(args.output_dir, "images")
        os.makedirs(images_dir, exist_ok=True)
        print(f"\nImages will be saved to: {images_dir}")

    # Generate all task arguments
    print("\nGenerating task list...")
    tasks = []
    for threshold in thresholds:
        for concentration in args.concentrations:
            # Format filename with underscores instead of dots
            conc_str = str(concentration).replace(".", "_")
            thresh_str = f"{threshold:.2f}".replace(".", "_")
            output_path = os.path.join(
                args.output_dir, f"azt_c{conc_str}_t{thresh_str}.graphml"
            )

            # Image path in separate subfolder
            image_path = None
            if args.save_images:
                image_path = os.path.join(
                    images_dir, f"azt_c{conc_str}_t{thresh_str}.{args.image_format}"
                )

            tasks.append(
                (
                    concentration,
                    threshold,
                    azt_long_df,
                    azt_pairs_df,  # Shared data
                    output_path,
                    fixed_params,
                    args.save_images,
                    image_path,
                    args.image_format,
                )
            )

    print(f"Generated {len(tasks)} tasks")

    # Build graphs in parallel
    print("\nBuilding graphs in parallel...")
    print("(This may take 15-20 minutes...)")
    print("-" * 80)

    with Pool(processes=num_workers) as pool:
        results = pool.map(build_single_graph, tasks)

    # Report results
    print("\n" + "=" * 80)
    print("Build Results")
    print("=" * 80)

    successes = []
    failures = []

    for concentration, threshold, success, info in results:
        if success:
            successes.append((concentration, threshold, info))
        else:
            failures.append((concentration, threshold, info))
            print(f"✗ FAILED: c={concentration}, t={threshold}: {info}")

    # Summary
    print("\n" + "=" * 80)
    print(f"Summary: {len(successes)} succeeded, {len(failures)} failed")
    print("=" * 80)

    if failures:
        print("\nFailed graphs:")
        for concentration, threshold, error in failures:
            print(f"  c={concentration}, t={threshold}: {error}")
        sys.exit(1)
    else:
        print("\n✓ All graphs built successfully!")
        print(f"\nOutput directory: {args.output_dir}")
        print(f"Total files: {len(successes)}")


if __name__ == "__main__":
    main()
