# This script contains the code to make weblogos and calculate logo strings for the fitness landscape graph.
import json

import logomaker as lm
import matplotlib.pyplot as plt

sign_pos = [19, 37, 67, 102, 162, 180, 235, 236, 237, 241, 261, 271, 272]

intended = {
    19: ["L", "P"],
    37: ["Q", "K"],
    67: ["M", "L", "V"],
    102: ["E", "K"],
    162: ["R", "S", "H", "N"],
    180: ["M", "T"],
    235: ["A", "T"],
    236: ["G", "S"],
    237: ["E", "K"],
    241: ["R", "S", "C"],
    261: ["T", "M"],
    271: ["R", "L", "Q"],
    272: ["N", "D"],
}

all_mutations = {
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


def convert_profile_to_intended(mutant_profile, all_mutations, intended):
    # Create mapping from position index to wildtype amino acid
    wt_map = {pos: intended[pos][0] for pos in intended}

    # Convert string to list for easier manipulation
    profile_list = list(mutant_profile)

    # Get positions in order
    positions = sorted(intended.keys())

    # For each position in the profile
    for i, pos in enumerate(positions):
        # If it's a dot, replace with wildtype amino acid
        if profile_list[i] == ".":
            profile_list[i] = wt_map[pos]

    return "".join(profile_list)


def convert_profile_to_dots(mutant_profile, intended):
    # Create mapping from position index to wildtype amino acid
    wt_map = {pos: intended[pos][0] for pos in intended}

    # Convert string to list for easier manipulation
    profile_list = list(mutant_profile)

    # Get positions in order
    positions = sorted(intended.keys())

    # For each position in the profile
    for i, pos in enumerate(positions):
        # If it matches wildtype, replace with dot
        if profile_list[i] == wt_map[pos]:
            profile_list[i] = "."

    return "".join(profile_list)


def convert_dots_to_intended_sequences(mutant_seqs, all_mutations, intended):
    """Convert sequences from dot notation to intended sequences.

    Args:
        mutant_seqs: List of sequences in dot notation
        all_mutations: Dict of allowed mutations
        intended: Dict mapping positions to allowed mutations

    Returns:
        converted_seqs: List of converted sequences
        wildtype: String of wildtype sequence
    """
    # Convert each sequence
    converted_seqs = []
    for seq in mutant_seqs:
        converted = convert_profile_to_intended(seq, all_mutations, intended)
        converted_seqs.append(converted)

    # Extract wildtype sequence
    wildtype = ""
    for pos in sorted(intended.keys()):
        wildtype += intended[pos][0]

    # # Print summary info
    # print("Intended sequence positions:")
    # for pos in sorted(intended.keys()):
    #     print(f"Position {pos}: {intended[pos][0]} (wildtype)")
    #     print(f"         Allowed mutations: {', '.join(intended[pos][1:])}")

    # print("\nSample of dot notation -> intended sequences:")
    # for orig, conv in zip(mutant_seqs[:5], converted_seqs[:5]):
    #     print(f"{orig} -> {conv}")

    # print(f"Wildtype sequence: {wildtype}")

    return converted_seqs, wildtype


def counts_to_string(counts_mat, threshold_ratio=0.8, wildtype=None, use_dots=False):
    """Convert count matrix to string, showing multiple characters when frequencies are similar.

    Matches the format: L[KQ]MKNMAGKR[TM]RN -> .[.K].KN...K.[.M]..

    Parameters:
    counts_mat: pandas DataFrame of counts
    threshold_ratio: float between 0 and 1, characters with frequency >= max_freq * threshold_ratio will be included
    wildtype: str, wildtype sequence (optional)
    use_dots: bool, if True, converts characters matching wildtype to dots
    """
    logo_chars = []

    for pos, pos_counts in counts_mat.iterrows():
        max_count = pos_counts.max()
        threshold = max_count * threshold_ratio

        # Get all characters above threshold
        significant_chars = pos_counts[pos_counts >= threshold].index.tolist()
        significant_chars.sort(key=lambda x: pos_counts[x], reverse=True)

        if len(significant_chars) == 1:
            char = significant_chars[0]
            if use_dots and wildtype and pos < len(wildtype) and char == wildtype[pos]:
                logo_chars.append(".")
            else:
                logo_chars.append(char)
        else:
            if use_dots and wildtype and pos < len(wildtype):
                # Convert each character to dot if it matches wildtype
                converted_chars = []
                for char in significant_chars:
                    if char == wildtype[pos]:
                        converted_chars.append(".")
                    else:
                        converted_chars.append(char)
                logo_chars.append(f"[{''.join(converted_chars)}]")
            else:
                logo_chars.append(f"[{''.join(significant_chars)}]")

    return "".join(logo_chars)


def mutant_dict_to_string(
    mutant_dict, all_mutations, intended, threshold_ratio=0.6, use_dots=True
):
    """Convert a mutant dictionary to a string representation from the count matrix.

    Args:
        mutant_dict: Dictionary of mutant sequences, which is a node attribute in the graph
        all_mutations: Dict of allowed mutations
        intended: Dict mapping positions to allowed mutations
        threshold_ratio: float between 0 and 1, characters with frequency >= max_freq * threshold_ratio will be included
        use_dots: bool, if True, converts characters matching wildtype to dots

    Returns:
        logo_string: String representation of the count matrix
    """
    # Convert mutant dictionary to list of sequences
    mutant_seqs = list(mutant_dict.keys())

    # Convert sequences from dot notation to intended sequences
    converted_seqs, wildtype = convert_dots_to_intended_sequences(
        mutant_seqs, all_mutations, intended
    )

    # Create count matrix from converted sequences
    counts_mat = lm.alignment_to_matrix(converted_seqs, characters_to_ignore="-")

    # Convert count matrix to string
    logo_string = counts_to_string(
        counts_mat,
        threshold_ratio=threshold_ratio,
        wildtype=wildtype,
        use_dots=use_dots,
    )

    return logo_string


def style_logo_colors(
    logo, sequence_length, wildtype, intended, wt_color="grey", mut_color="red"
):
    """Style logo colors with wildtype and mutation colors.

    Args:
        logo: Logomaker Logo object to style
        sequence_length: Length of the sequences in the logo
        wildtype: String of wildtype sequence
        intended: Dict mapping positions to allowed mutations
        wt_color: Color to use for wildtype residues (default: grey)
        mut_color: Color to use for mutation residues (default: red)
    """
    for pos in range(sequence_length):
        wt = wildtype[pos]
        mutations = list(intended.values())[pos][1:]  # Get allowed mutations

        # Only style wildtype if it exists in the logo
        if wt in logo.glyph_df.columns:
            logo.style_single_glyph(p=pos, c=wt, color=wt_color)

        # Only style mutations that exist in the logo
        for mut in mutations:
            if mut in logo.glyph_df.columns:
                logo.style_single_glyph(p=pos, c=mut, color=mut_color)


def plot_logo(
    mutant_seqs,
    all_mutations,
    intended,
    threshold_ratio=0.6,
    figsize=(10, 2),
    save_path=None,
):
    """Plot sequence logo from mutant dictionary.

    Args:
        mutant_seqs: Mutant sequences to plot
        all_mutations: Dict of allowed mutations at each position
        intended: Dict mapping positions to allowed mutations
        threshold_ratio: float between 0 and 1, characters with frequency >= max_freq * threshold_ratio will be included
        figsize: Tuple of (width, height) for figure size
        save_path: Optional path to save figure. If None, figure is displayed but not saved.

    Returns:
        matplotlib.figure.Figure: The generated logo figure
    """
    # Convert sequences from dot notation to intended sequences
    converted_seqs, wildtype = convert_dots_to_intended_sequences(
        mutant_seqs, all_mutations, intended
    )

    # Create count matrix from converted sequences
    counts_mat = lm.alignment_to_matrix(converted_seqs, characters_to_ignore="-")

    # Create and style logo
    fig, ax = plt.subplots(figsize=figsize)
    logo = lm.Logo(counts_mat, ax=ax, vpad=0.05)

    # Style colors
    style_logo_colors(logo, len(wildtype), wildtype, intended)

    # Define Ambler positions
    sign_pos_ambler = [21, 39, 69, 104, 164, 182, 237, 238, 240, 244, 265, 275, 276]

    # Set x-axis ticks and labels to Ambler positions
    ax.set_xticks(range(len(sign_pos_ambler)))
    ax.set_xticklabels(sign_pos_ambler)

    # Customize appearance
    ax.set_xlabel("Position")
    ax.set_ylabel("Frequency")

    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=1000)
        plt.close()

    return fig


if __name__ == "__main__":
    mutant_dict_file = (
        "/work/greencenter/s439821/TEM1CML/data/processed/node_group_mutants.json"
    )
    with open(mutant_dict_file) as f:
        mutant_dict = json.load(f)

    logo_string = mutant_dict_to_string(
        mutant_dict, all_mutations, intended, threshold_ratio=0.6, use_dots=True
    )
    print(logo_string)

    filename = "amp_3_1_big.png"
    save_path = (
        f"/work/greencenter/s439821/TEM1CML/data/outputs/figures/logo/graph_{filename}"
    )
    mutant_seqs = list(mutant_dict.keys())
    fig = plot_logo(
        mutant_seqs,
        all_mutations,
        intended,
        threshold_ratio=0.6,
        figsize=(8, 2),
        save_path=save_path,
    )
    plt.show()
