# This script contains the code that generate the global fitness table and pairs table for the fitness landscape graph.

import argparse
from collections import defaultdict

import numpy as np
import polars as pl

from fitness_landscape_graph.preprocess import preprocess_data


def calculate_normalized_fitness(long_df):
    """Calculate normalized fitness from a long-format dataframe containing replicate measurements.

    Normalizes fitness relative to wildtype (all dots) and dead mutant (all X's).

    Parameters:
    -----------
    long_df : polars.DataFrame
        DataFrame in long format with columns:
        - mutant_profile: str, mutation profile
        - concentration: float, drug concentration
        - replicate1, replicate2, replicate3: float, replicate measurements
        - median: float, median of replicates

    Returns:
    --------
    polars.DataFrame
        DataFrame with columns:
        - mutant_profile: str
        - auc1, auc2, auc3: float, area under curve for each replicate
        - fitness: float, mean fitness across replicates
        - error: float, standard error
        - normalized_fitness: float, normalized fitness relative to wildtype
    """
    # Define reference sequences
    dead_mutant = "XXXXXXXXXXXXX"
    wildtype_mutant = "............."

    # Calculate AUC for each replicate
    auc_df = (
        long_df.group_by("mutant_profile")
        .agg(
            pl.col("replicate1"),
            pl.col("replicate2"),
            pl.col("replicate3"),
            pl.col("concentration"),
        )
        .with_columns(
            [
                # Calculate AUC for each replicate
                pl.struct(["replicate1", "concentration"])
                .map_elements(
                    lambda x: np.log10(
                        np.trapezoid(
                            y=np.power(10, np.array(x["replicate1"])),
                            x=np.array(x["concentration"]),
                        )
                    )
                    if not np.all(np.isnan(x["replicate1"]))
                    else np.nan,
                    return_dtype=pl.Float64,
                )
                .alias("auc1"),
                pl.struct(["replicate2", "concentration"])
                .map_elements(
                    lambda x: np.log10(
                        np.trapezoid(
                            y=np.power(10, np.array(x["replicate2"])),
                            x=np.array(x["concentration"]),
                        )
                    )
                    if not np.all(np.isnan(x["replicate2"]))
                    else np.nan,
                    return_dtype=pl.Float64,
                )
                .alias("auc2"),
                pl.struct(["replicate3", "concentration"])
                .map_elements(
                    lambda x: np.log10(
                        np.trapezoid(
                            y=np.power(10, np.array(x["replicate3"])),
                            x=np.array(x["concentration"]),
                        )
                    )
                    if not np.all(np.isnan(x["replicate3"]))
                    else np.nan,
                    return_dtype=pl.Float64,
                )
                .alias("auc3"),
            ]
        )
        # Calculate statistics
        .with_columns(
            [
                pl.concat_list(["auc1", "auc2", "auc3"])
                .map_elements(
                    lambda x: np.nanmean(x), return_dtype=pl.Float64
                )  # consistent with Erdal's
                .alias("fitness"),
                pl.concat_list(["auc1", "auc2", "auc3"])
                .map_elements(lambda x: np.nanstd(x), return_dtype=pl.Float64)
                .alias("error"),
            ]
        )
        .select(["mutant_profile", "auc1", "auc2", "auc3", "fitness", "error"])
    )

    # Get fitness values for reference sequences
    f_dead = auc_df.filter(pl.col("mutant_profile") == dead_mutant)["fitness"].item()
    f_wt = auc_df.filter(pl.col("mutant_profile") == wildtype_mutant)["fitness"].item()

    # Normalize fitness relative to wildtype and dead mutant
    global_fitness_df = auc_df.with_columns(
        ((pl.col("fitness") - f_dead) / (f_wt - f_dead)).alias("normalized_fitness")
    )

    return global_fitness_df


def get_pairs(global_fitness_df):
    """Generate pairs of mutants that differ by one character and compute fitness differences.

    Args:
    global_fitness_df (pl.DataFrame): Input dataframe with normalized fitness data

    Returns:
    pl.DataFrame: Dataframe with mutant pairs and their fitness differences
    """

    def string_difference(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=True))

    mutant_profiles = global_fitness_df["mutant_profile"].unique().to_list()

    # Generate mutant pairs
    char_diff_dict = defaultdict(set)
    for i, mutant in enumerate(mutant_profiles):
        for j in range(len(mutant)):
            mutated = mutant[:j] + "?" + mutant[j + 1 :]
            char_diff_dict[mutated].add(i)

    mutant_pairs = set()
    for indices in char_diff_dict.values():
        for i in indices:
            for j in indices:
                if i < j:
                    m1, m2 = mutant_profiles[i], mutant_profiles[j]
                    if string_difference(m1, m2) == 1:
                        mutant_pairs.add((m1, m2))

    pairs = []
    for m1, m2 in mutant_pairs:
        m1_fitness = global_fitness_df.filter(pl.col("mutant_profile") == m1)[
            "fitness"
        ].item()
        m2_fitness = global_fitness_df.filter(pl.col("mutant_profile") == m2)[
            "fitness"
        ].item()
        fitness_diff = m1_fitness - m2_fitness
        pairs.append(
            {"mutant_profile1": m1, "mutant_profile2": m2, "fitness_diff": fitness_diff}
        )

    return pl.DataFrame(pairs)


def main():
    parser = argparse.ArgumentParser(
        description="Generate global fitness table and pairs table for fitness landscape graph."
    )
    parser.add_argument(
        "--base-path",
        type=str,
        help="Path to the project directory",
    )
    args = parser.parse_args()
    base_path = args.base_path

    # File paths
    amp_path = f"{base_path}/data/raw/combined-auc/genotype_auc_sorted_ampicillin.csv"
    azt_path = f"{base_path}/data/raw/combined-auc/genotype_auc_sorted_aztreonam.csv"
    amp_pairs_path = f"{base_path}/data/processed/amp_pairs.csv"
    azt_pairs_path = f"{base_path}/data/processed/azt_pairs.csv"

    # Load and preprocess data
    processed_data = preprocess_data(
        amp_path, azt_path, amp_pairs_path, azt_pairs_path, clean_nulls_flag=True
    )

    # Access the processed dataframes
    amp_long_df = processed_data["amp"]["long"]
    azt_long_df = processed_data["azt"]["long"]

    # Calculate normalized fitness
    amp_global_fitness_df = calculate_normalized_fitness(amp_long_df)
    azt_global_fitness_df = calculate_normalized_fitness(azt_long_df)

    # Get mutant pairs
    amp_pairs = get_pairs(amp_global_fitness_df)
    print("ampicillin pairs generated")
    azt_pairs = get_pairs(azt_global_fitness_df)
    print("aztreonam pairs generated")

    # Save results
    amp_global_fitness_df.sort("mutant_profile").write_csv(
        f"{base_path}/exps/exp-01-reproduce-results/outputs/reproduce-pairs/amp_global_fitness.csv"
    )
    azt_global_fitness_df.sort("mutant_profile").write_csv(
        f"{base_path}/exps/exp-01-reproduce-results/outputs/reproduce-pairs/azt_global_fitness.csv"
    )

    amp_pairs.write_csv(
        f"{base_path}/exps/exp-01-reproduce-results/outputs/reproduce-pairs/amp_global_pairs.csv"
    )
    azt_pairs.write_csv(
        f"{base_path}/exps/exp-01-reproduce-results/outputs/reproduce-pairs/azt_global_pairs.csv"
    )
    print("results saved")


if __name__ == "__main__":
    main()
