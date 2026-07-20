# This script is the main script that use GraphBuilder to build the fitness landscape graph.
import argparse
import json
import os

import networkx as nx
import polars as pl

from fitness_landscape_graph.graph_builder import GraphBuilder
from fitness_landscape_graph.pair_table_global import calculate_normalized_fitness
from fitness_landscape_graph.preprocess import preprocess_data


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build graphs from antibiotic resistance data"
    )
    parser.add_argument(
        "--base-path",
        type=str,
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        help="Directory to save the output graph files",
    )
    parser.add_argument(
        "--neutral-threshold",
        type=float,
        default=0.4,
        help="Threshold for merging neutral edges (for concentration-specific graphs)",
    )
    parser.add_argument(
        "--tiny-initial-threshold",
        type=float,
        default=0.02,
        help="Threshold for merging tiny initial edges",
    )
    parser.add_argument(
        "--large-edge-threshold",
        type=float,
        default=5.5,
        help="Threshold for identifying potential forbidden pairs",
    )
    parser.add_argument(
        "--num-forbidden-pairs",
        type=int,
        default=None,
        help="Number of forbidden pairs to select (if None, use all pairs above threshold)",
    )
    return parser.parse_args()


def convert_graph_attributes_to_json(graph):
    """Convert dictionary attributes in a NetworkX graph to JSON strings.

    Args:
        graph (nx.Graph): NetworkX graph whose attributes need to be converted

    Returns:
        nx.Graph: Graph with dictionary attributes converted to JSON strings
    """
    # Convert node attributes
    for node, data in graph.nodes(data=True):
        for key, value in data.items():
            if isinstance(value, dict):
                graph.nodes[node][key] = json.dumps(value)

    # Convert edge attributes
    for u, v, data in graph.edges(data=True):
        for key, value in data.items():
            if isinstance(value, dict):
                graph.edges[u, v][key] = json.dumps(value)

    return graph


def main():
    args = parse_args()
    print(args)

    # Define list of thresholds for global graphs
    global_thresholds = [0.4]
    print(f"Building global graphs with thresholds: {global_thresholds}")

    base_path = args.base_path

    # File paths
    amp_path = f"{base_path}/data/raw/combined-auc/genotype_auc_sorted_ampicillin.csv"
    azt_path = f"{base_path}/data/raw/combined-auc/genotype_auc_sorted_aztreonam.csv"
    amp_pairs_path = f"{base_path}/data/processed/amp_pairs.csv"
    azt_pairs_path = f"{base_path}/data/processed/azt_pairs.csv"
    amp_global_pairs_path = f"{base_path}/data/processed/amp_global_pairs.csv"
    azt_global_pairs_path = f"{base_path}/data/processed/azt_global_pairs.csv"

    # Load and preprocess data
    processed_data = preprocess_data(
        amp_path, azt_path, amp_pairs_path, azt_pairs_path, clean_nulls_flag=True
    )

    # Access the processed dataframes
    amp_long_df = processed_data["amp"]["long"]
    amp_pairs_df = processed_data["amp"]["pairs"]

    azt_long_df = processed_data["azt"]["long"]
    azt_pairs_df = processed_data["azt"]["pairs"]

    amp_global_fitness_df = calculate_normalized_fitness(amp_long_df)
    azt_global_fitness_df = calculate_normalized_fitness(azt_long_df)
    amp_global_pairs_df = pl.read_csv(amp_global_pairs_path)
    azt_global_pairs_df = pl.read_csv(azt_global_pairs_path)

    # Get the concentrations for ampicillin and aztreonam
    amp_concs = amp_long_df["concentration"].unique().sort()
    azt_concs = azt_long_df["concentration"].unique().sort()

    # Dictionary to save all the graphs
    all_graphs = {"amp": {}, "azt": {}}

    # Build global graphs for each threshold
    for threshold in global_thresholds:
        # Create threshold-specific output directory
        threshold_str = f"{threshold:.2f}".replace(".", "_")
        threshold_output_dir = os.path.join(
            args.output_dir, f"global_threshold_{threshold_str}"
        )
        os.makedirs(threshold_output_dir, exist_ok=True)

        # Build ampicillin global graph
        graph_builder = GraphBuilder(
            long_table=amp_global_fitness_df, pairwise_table=amp_global_pairs_df
        )
        amp_global_graph = graph_builder.build_graph(
            concentration=None,
            neutral_threshold=threshold,  # Use the current threshold from list
            fitness_col="fitness",
            fitness_diff_col="fitness_diff",
            use_parallel_peak_merge=True,
            rename_by_avg_fitness=True,
            merge_peaks=False,
            large_edge_threshold=args.large_edge_threshold,
            tiny_initial_threshold=args.tiny_initial_threshold,
            num_forbidden_pairs=args.num_forbidden_pairs,
        )

        # Build aztreonam global graph
        graph_builder = GraphBuilder(
            long_table=azt_global_fitness_df, pairwise_table=azt_global_pairs_df
        )
        azt_global_graph = graph_builder.build_graph(
            concentration=None,
            neutral_threshold=threshold,  # Use the current threshold from list
            fitness_col="fitness",
            fitness_diff_col="fitness_diff",
            use_parallel_peak_merge=True,
            rename_by_avg_fitness=True,
            merge_peaks=False,
            large_edge_threshold=args.large_edge_threshold,
            tiny_initial_threshold=args.tiny_initial_threshold,
            num_forbidden_pairs=args.num_forbidden_pairs,
        )

        # Convert and save graphs
        amp_global_graph = convert_graph_attributes_to_json(amp_global_graph)
        azt_global_graph = convert_graph_attributes_to_json(azt_global_graph)

        nx.write_graphml(
            amp_global_graph, f"{threshold_output_dir}/amp_global_t{threshold}.graphml"
        )
        nx.write_graphml(
            azt_global_graph, f"{threshold_output_dir}/azt_global_t{threshold}.graphml"
        )

        print(
            f"Saved global graphs with neutral threshold {threshold} to {threshold_output_dir}"
        )

    # Process each concentration for both antibiotics
    for conc in amp_concs:
        graph_builder = GraphBuilder(
            long_table=amp_long_df, pairwise_table=amp_pairs_df
        )
        all_graphs["amp"][conc] = graph_builder.build_graph(
            concentration=conc,
            neutral_threshold=args.neutral_threshold,
            fitness_col="median",
            fitness_diff_col="median_diff",
            use_parallel_peak_merge=True,
            rename_by_avg_fitness=True,
            tiny_initial_threshold=args.tiny_initial_threshold,
            merge_peaks=False,
            large_edge_threshold=args.large_edge_threshold,
            num_forbidden_pairs=args.num_forbidden_pairs,
        )

    for conc in azt_concs:
        graph_builder = GraphBuilder(
            long_table=azt_long_df, pairwise_table=azt_pairs_df
        )
        all_graphs["azt"][conc] = graph_builder.build_graph(
            concentration=conc,
            neutral_threshold=args.neutral_threshold,
            fitness_col="median",
            fitness_diff_col="median_diff",
            use_parallel_peak_merge=True,
            rename_by_avg_fitness=True,
            tiny_initial_threshold=args.tiny_initial_threshold,
            merge_peaks=False,
            large_edge_threshold=args.large_edge_threshold,
            num_forbidden_pairs=args.num_forbidden_pairs,
        )

    for antibiotic, conc_graphs in all_graphs.items():
        for conc, graph in conc_graphs.items():
            # Replace dots with underscores in the concentration value
            conc_str = str(conc).replace(".", "_")
            output_path = os.path.join(
                args.output_dir, f"{antibiotic}_{conc_str}.graphml"
            )

            graph = convert_graph_attributes_to_json(graph)

            nx.write_graphml(graph, output_path)
            print(
                f"Saved graph for {antibiotic} at concentration {conc} to {output_path}"
            )


if __name__ == "__main__":
    main()
