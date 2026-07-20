import time
from collections import defaultdict
from itertools import product

import numpy as np
import polars as pl
from scipy import stats


def clean_nulls(df):
    df_with_null_counts = df.with_columns(
        pl.sum_horizontal(pl.all().is_null()).alias("null_count_per_row")
    )
    threshold = (
        (df.shape[1] - 3) * 2 / 3
    )  # delete rows with more than 2 replicates missing
    filtered_df = df_with_null_counts.filter(pl.col("null_count_per_row") < threshold)

    # remove the null count column
    filtered_df = filtered_df.drop("null_count_per_row")
    return filtered_df


def convert_long_format(df, concs):
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
    long_mutants = []
    long_concs = []
    long_rep1 = []
    long_rep2 = []
    long_rep3 = []

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


def get_mutant_pairs(long_df: pl.DataFrame) -> pl.DataFrame:
    concs = long_df["concentration"].unique().to_numpy()
    mutant_profiles = long_df["mutant_profile"].unique().to_numpy()

    def string_difference(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=True))

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

    print("the number of pairs is: ", len(mutant_pairs))

    n_pairs = len(mutant_pairs)
    n_concs = len(concs)
    total_rows = n_pairs * n_concs

    results = {
        "mutant_profile1": [None] * total_rows,
        "mutant_profile2": [None] * total_rows,
        "concentration": [None] * total_rows,
        "t_stat": [None] * total_rows,
        "p_value": [None] * total_rows,
        "median_diff": [None] * total_rows,
        "profile1_rep1": [None] * total_rows,
        "profile1_rep2": [None] * total_rows,
        "profile1_rep3": [None] * total_rows,
        "profile2_rep1": [None] * total_rows,
        "profile2_rep2": [None] * total_rows,
        "profile2_rep3": [None] * total_rows,
    }

    idx = 0
    for conc in concs:
        conc_df = long_df.filter(pl.col("concentration") == conc)

        for m1, m2 in mutant_pairs:
            m1_df = conc_df.filter(pl.col("mutant_profile") == m1)
            m2_df = conc_df.filter(pl.col("mutant_profile") == m2)

            rep1 = m1_df.select(["replicate1", "replicate2", "replicate3"]).to_numpy()[
                0
            ]
            rep2 = m2_df.select(["replicate1", "replicate2", "replicate3"]).to_numpy()[
                0
            ]

            t_stat_val, p_val = stats.ttest_ind(rep1, rep2, nan_policy="omit")
            median_diff_val = np.nanmedian(rep1) - np.nanmedian(rep2)

            results["mutant_profile1"][idx] = str(m1)
            results["mutant_profile2"][idx] = str(m2)
            results["concentration"][idx] = float(conc)
            results["t_stat"][idx] = float(t_stat_val)
            results["p_value"][idx] = float(p_val)
            results["median_diff"][idx] = float(median_diff_val)
            results["profile1_rep1"][idx] = float(rep1[0])
            results["profile1_rep2"][idx] = float(rep1[1])
            results["profile1_rep3"][idx] = float(rep1[2])
            results["profile2_rep1"][idx] = float(rep2[0])
            results["profile2_rep2"][idx] = float(rep2[1])
            results["profile2_rep3"][idx] = float(rep2[2])
            idx += 1

    # Remove unused pre-allocated rows
    for key in results:
        results[key] = results[key][:idx]

    return pl.DataFrame(results)


amp_path = "/work/greencenter/s439821/fitness-landscape-graph/data/raw/combined-auc/genotype_auc_sorted_ampicillin.csv"
azt_path = "/work/greencenter/s439821/fitness-landscape-graph/data/raw/combined-auc/genotype_auc_sorted_aztreonam.csv"

amp_df = pl.read_csv(amp_path)
azt_df = pl.read_csv(azt_path)

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

mutants_id_arr = amp_df["mut_profile_masked"].to_numpy()
mutants_id_arr = mutants_id_arr.astype(str)

mutation_combinations = product(*[intended[pos] for pos in sorted(intended.keys())])

# convert the mutation combinations to strings
mutation_combinations = ["".join(mut) for mut in mutation_combinations]
possible_intended_mutants = np.array(mutation_combinations)

# clean the data so that only rows that have the intended mutations are kept
amp_df = amp_df.filter(amp_df["mut_profile_masked"].is_in(possible_intended_mutants))
azt_df = azt_df.filter(azt_df["mut_profile_masked"].is_in(possible_intended_mutants))

amp_df = clean_nulls(amp_df)
azt_df = clean_nulls(azt_df)

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

amp_long_df = convert_long_format(amp_df, amp_concs)
azt_long_df = convert_long_format(azt_df, azt_concs)

start = time.time()
amp_pairs = get_mutant_pairs(amp_long_df)
end = time.time()
print("amp time: ", end - start)

start = time.time()
azt_pairs = get_mutant_pairs(azt_long_df)
end = time.time()
print("azt time: ", end - start)

save_path = "/work/greencenter/s439821/fitness-landscape-graph/exps/exp-01-reproduce-results/outputs/reproduce-pairs"
amp_pairs.write_csv(f"{save_path}/amp_pairs.csv")
azt_pairs.write_csv(f"{save_path}/azt_pairs.csv")
print("done!")
