"""Preprocesses experimental data on antibiotic resistance.

Main function:
- preprocess_data: Reads raw data, filters for intended mutants, cleans nulls,
                   converts to long format, and processes mutant pairs.

Helper functions:
- clean_nulls
- convert_long_format
- get_mutant_pairs

Uses polars for dataframes, numpy for computations, and itertools for combinations.
"""

from itertools import product

import numpy as np
import polars as pl


def preprocess_data(
    amp_path, azt_path, amp_pairs_path, azt_pairs_path, clean_nulls_flag
):
    """Preprocess the data from input files and return processed dataframes.

    Args:
    amp_path (str): Path to ampicillin data file
    azt_path (str): Path to aztreonam data file
    amp_pairs_path (str): Path to ampicillin pairs file
    azt_pairs_path (str): Path to aztreonam pairs file
    clean_nulls_flag (bool): Whether to clean null values from the data (default: True)

    Returns:
    dict: A dictionary containing processed dataframes for both drugs
        - 'drug': {
            'original': original dataframe, filtered for intended mutants.
            'long': long format dataframe,
            'pairs': pairs dataframe
        }
    """
    # Read raw data from sequencing
    amp_df = pl.read_csv(amp_path)
    azt_df = pl.read_csv(azt_path)

    # Define intended mutations
    intended = {
        19: [".", "P"],
        37: [".", "K"],
        67: [".", "L", "V"],
        102: [".", "K"],
        162: [".", "S", "H", "N"],
        180: [".", "T"],
        235: [".", "T"],
        236: [".", "S"],
        237: [".", "K"],
        241: [".", "S", "C"],
        261: [".", "M"],
        271: [".", "L", "Q"],
        272: [".", "D"],
    }

    # we only include intended mutants in the analysis
    mutation_combinations = product(*[intended[pos] for pos in sorted(intended.keys())])
    mutation_combinations = ["".join(mut) for mut in mutation_combinations]
    possible_intended_mutants = np.array(mutation_combinations)
    wildtype = "............."
    dead_mutant = "XXXXXXXXXXXXX"
    possible_intended_mutants = np.append(
        possible_intended_mutants, [wildtype, dead_mutant], axis=0
    )
    # Filter data for intended mutations
    amp_df = amp_df.filter(
        amp_df["mut_profile_masked"].is_in(possible_intended_mutants)
    )
    azt_df = azt_df.filter(
        azt_df["mut_profile_masked"].is_in(possible_intended_mutants)
    )

    # Clean null values if flag is True
    if clean_nulls_flag:
        amp_df = clean_nulls(amp_df, num_rep_threshold=2)
        azt_df = clean_nulls(azt_df, num_rep_threshold=2)

    # Define concentrations
    amp_concs = [
        "Ampicillin 0.0",
        "Ampicillin 3.1",
        "Ampicillin 12.2",
        "Ampicillin 48.8",
        "Ampicillin 195.0",
        "Ampicillin 781.0",
    ]
    azt_concs = [
        "Aztreonam 0.0",
        "Aztreonam 0.44",
        "Aztreonam 1.33",
        "Aztreonam 4.0",
        "Aztreonam 12.0",
        "Aztreonam 36.0",
        "Aztreonam 108.0",
        "Aztreonam 324.0",
    ]

    # Convert to long format
    amp_long_df = convert_long_format(amp_df, amp_concs)
    azt_long_df = convert_long_format(azt_df, azt_concs)

    # Add median column
    amp_long_df = amp_long_df.with_columns(
        pl.struct(["replicate1", "replicate2", "replicate3"])
        .map_elements(
            lambda s: np.nanmedian([s["replicate1"], s["replicate2"], s["replicate3"]]),
            return_dtype=pl.Float64,
        )
        .alias("median")
    )
    azt_long_df = azt_long_df.with_columns(
        pl.struct(["replicate1", "replicate2", "replicate3"])
        .map_elements(
            lambda s: np.nanmedian([s["replicate1"], s["replicate2"], s["replicate3"]]),
            return_dtype=pl.Float64,
        )
        .alias("median")
    )

    # Read pairs data
    amp_pairs_df = pl.read_csv(amp_pairs_path)
    azt_pairs_df = pl.read_csv(azt_pairs_path)

    # Create a dictionary to store all dataframes
    processed_data = {
        "amp": {"original": amp_df, "long": amp_long_df, "pairs": amp_pairs_df},
        "azt": {"original": azt_df, "long": azt_long_df, "pairs": azt_pairs_df},
    }

    return processed_data


def clean_nulls(df, num_rep_threshold=2):
    """Remove rows with too many missing replicates from the dataframe.

    Args:
    df (pl.DataFrame): Input dataframe
    num_rep_threshold (int): Number of replicates missing to be considered too many

    Returns:
    pl.DataFrame: Cleaned dataframe with rows removed that have more than num_rep_threshold/3 of replicates missing
    """
    df_with_null_counts = df.with_columns(
        pl.sum_horizontal(pl.all().is_null()).alias("null_count_per_row")
    )
    threshold = (
        (df.shape[1] - 3) * num_rep_threshold / 3
    )  # delete rows with more than num_rep_threshold replicates missing
    filtered_df = df_with_null_counts.filter(pl.col("null_count_per_row") < threshold)

    # remove the null count column
    filtered_df = filtered_df.drop("null_count_per_row")
    return filtered_df


def convert_long_format(df, concs):
    """Convert the dataframe from wide to long format.

    In this context, converting to long format means:
    - Each mutant profile, concentration, and replicate measurement becomes its own row
    - Instead of having multiple columns for different concentrations and replicates,
      we have single columns for 'concentration' and 'replicate1', 'replicate2', 'replicate3'
    - This format makes it easier to perform analyses and visualizations across different
      concentrations and replicates

    The resulting long format dataframe will have the following columns:
    - mutant_profile: The genetic profile of the mutant
    - concentration: The concentration of the antibiotic
    - replicate1, replicate2, replicate3: The measurements for each replicate

    Args:
    df (pl.DataFrame): Input dataframe in wide format
    concs (list): List of concentration values

    Returns:
    pl.DataFrame: Dataframe in long format with columns for mutant_profile, concentration, and replicates
    """
    # add replicate to concs
    rep_concs = [f"{conc} {i + 1}" for conc in concs for i in range(3)]
    data_array = df[rep_concs].to_numpy()

    # Get the mutant profiles
    mutant_profiles = df["mut_profile_masked"].to_list()

    # Reshape the array
    n_mutants, n_columns = data_array.shape
    n_concs = len(concs)
    n_replicates = n_columns // n_concs

    reshaped_data = data_array.reshape(n_mutants, n_concs, n_replicates)

    # Create lists for the long format data
    long_mutants, long_concs, long_rep1, long_rep2, long_rep3 = [], [], [], [], []

    for i, mutant in enumerate(mutant_profiles):
        for j, conc in enumerate(concs):
            long_mutants.append(mutant)

            # only save the concentration number
            conc = conc.split(" ")[1]
            long_concs.append(float(conc))
            long_rep1.append(reshaped_data[i, j, 0])
            long_rep2.append(reshaped_data[i, j, 1] if n_replicates > 1 else None)
            long_rep3.append(reshaped_data[i, j, 2] if n_replicates > 2 else None)

    # Create the long format dataframe
    long_df = pl.DataFrame(
        {
            "mutant_profile": long_mutants,
            "concentration": long_concs,
            "replicate1": long_rep1,
            "replicate2": long_rep2,
            "replicate3": long_rep3,
        }
    )

    return long_df
