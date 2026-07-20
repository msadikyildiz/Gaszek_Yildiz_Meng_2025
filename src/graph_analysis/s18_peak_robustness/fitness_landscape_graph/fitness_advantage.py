"""Fitness advantage analysis for fitness landscape graphs.

Provides the FitnessAdvantageAnalyzer class to compute fitness advantages
of a group of genotypes (e.g., a peak supernode) over external neighbors
across multiple drug concentrations.
"""

import numpy as np
import polars as pl


class FitnessAdvantageAnalyzer:
    """Analyze fitness advantages of genotype groups over external neighbors.

    Uses pre-computed fitness values from the preprocessing pipeline
    (long_df.median) and pair topology from pairs_df. Vectorized NumPy
    operations for efficient 2-mutation neighbor finding.
    """

    def __init__(self, pairs_df: pl.DataFrame, long_df: pl.DataFrame):
        """Initialize the analyzer.

        Args:
            pairs_df: Pairwise mutation table with columns:
                mutant_profile1, mutant_profile2, concentration, median_diff.
            long_df: Long-format fitness table with columns:
                mutant_profile, concentration, median.
        """
        self.pairs_df = pairs_df
        self.long_df = long_df

    def compute_fitness_advantage(
        self,
        group_genotypes: set[str],
        concentrations: list[float] | None = None,
        max_distance: int = 2,
    ) -> pl.DataFrame:
        """Compute fitness advantage of group members over external neighbors.

        Args:
            group_genotypes: Set of genotype strings in the target group.
            concentrations: Concentrations to analyze. If None, uses all.
            max_distance: Maximum mutation distance (1 or 2).

        Returns:
            DataFrame with columns: concentration, group_member, neighbor,
                distance, group_member_fitness, neighbor_fitness, fitness_diff.
        """
        if concentrations is None:
            concentrations = sorted(
                self.long_df["concentration"].unique().to_list()
            )

        empty_schema = {
            "group_member": pl.Utf8,
            "neighbor": pl.Utf8,
            "distance": pl.Int32,
            "group_member_fitness": pl.Float64,
            "neighbor_fitness": pl.Float64,
            "fitness_diff": pl.Float64,
            "concentration": pl.Float64,
        }

        all_dfs = []
        for conc in concentrations:
            conc_pairs = self.pairs_df.filter(pl.col("concentration") == conc)
            conc_long = self.long_df.filter(pl.col("concentration") == conc)

            dfs = []

            if max_distance >= 1:
                df_1mut = self._find_1mut_neighbors(
                    conc_pairs, conc_long, group_genotypes
                )
                if df_1mut.height > 0:
                    dfs.append(df_1mut)

            if max_distance >= 2:
                df_2mut = self._find_2mut_neighbors(conc_long, group_genotypes)
                if df_2mut.height > 0:
                    dfs.append(df_2mut)

            if dfs:
                df = pl.concat(dfs).with_columns(
                    pl.lit(conc).alias("concentration")
                )
                all_dfs.append(df)

        if not all_dfs:
            return pl.DataFrame(schema=empty_schema)

        return pl.concat(all_dfs)

    def _find_1mut_neighbors(
        self,
        conc_pairs: pl.DataFrame,
        conc_long: pl.DataFrame,
        group_genotypes: set[str],
    ) -> pl.DataFrame:
        """Find 1-mutation neighbors using the pair table.

        Uses median_diff from pairs_df directly (no recomputation).
        Joins with conc_long for individual fitness values.

        Args:
            conc_pairs: Pair table filtered to one concentration.
            conc_long: Long table filtered to one concentration.
            group_genotypes: Set of genotype strings in the group.

        Returns:
            DataFrame with columns: group_member, neighbor, distance,
                group_member_fitness, neighbor_fitness, fitness_diff.
        """
        # Case 1: group member is profile1, neighbor is profile2
        # median_diff = median(profile1) - median(profile2) = gm - nb
        p1_in = conc_pairs.filter(
            pl.col("mutant_profile1").is_in(group_genotypes)
            & ~pl.col("mutant_profile2").is_in(group_genotypes)
        ).select(
            group_member=pl.col("mutant_profile1"),
            neighbor=pl.col("mutant_profile2"),
            fitness_diff=pl.col("median_diff"),
        )

        # Case 2: group member is profile2, neighbor is profile1
        # median_diff = median(profile1) - median(profile2) = nb - gm
        # So fitness_diff = -median_diff
        p2_in = conc_pairs.filter(
            pl.col("mutant_profile2").is_in(group_genotypes)
            & ~pl.col("mutant_profile1").is_in(group_genotypes)
        ).select(
            group_member=pl.col("mutant_profile2"),
            neighbor=pl.col("mutant_profile1"),
            fitness_diff=-pl.col("median_diff"),
        )

        combined = pl.concat([p1_in, p2_in])

        if combined.height == 0:
            return pl.DataFrame(
                schema={
                    "group_member": pl.Utf8,
                    "neighbor": pl.Utf8,
                    "distance": pl.Int32,
                    "group_member_fitness": pl.Float64,
                    "neighbor_fitness": pl.Float64,
                    "fitness_diff": pl.Float64,
                }
            )

        # Join with long_df to get individual fitness values
        fitness_lookup = conc_long.select("mutant_profile", "median")

        combined = (
            combined.join(
                fitness_lookup.rename(
                    {"mutant_profile": "group_member", "median": "group_member_fitness"}
                ),
                on="group_member",
                how="left",
            )
            .join(
                fitness_lookup.rename(
                    {"mutant_profile": "neighbor", "median": "neighbor_fitness"}
                ),
                on="neighbor",
                how="left",
            )
            .with_columns(pl.lit(1).cast(pl.Int32).alias("distance"))
            .select(
                "group_member",
                "neighbor",
                "distance",
                "group_member_fitness",
                "neighbor_fitness",
                "fitness_diff",
            )
        )

        return combined

    def _find_2mut_neighbors(
        self,
        conc_long: pl.DataFrame,
        group_genotypes: set[str],
    ) -> pl.DataFrame:
        """Find 2-mutation neighbors using vectorized Hamming distance.

        Args:
            conc_long: Long table filtered to one concentration.
            group_genotypes: Set of genotype strings in the group.

        Returns:
            DataFrame with columns: group_member, neighbor, distance,
                group_member_fitness, neighbor_fitness, fitness_diff.
        """
        all_genotypes = set(conc_long["mutant_profile"].unique().to_list())
        external_genotypes = all_genotypes - group_genotypes
        # Only keep group members that exist in this concentration
        group_at_conc = group_genotypes & all_genotypes

        group_list = list(group_at_conc)
        external_list = list(external_genotypes)

        if not group_list or not external_list:
            return pl.DataFrame(
                schema={
                    "group_member": pl.Utf8,
                    "neighbor": pl.Utf8,
                    "distance": pl.Int32,
                    "group_member_fitness": pl.Float64,
                    "neighbor_fitness": pl.Float64,
                    "fitness_diff": pl.Float64,
                }
            )

        pairs_at_dist_2 = self._compute_hamming_distances_vectorized(
            group_list, external_list, target_distance=2
        )

        if not pairs_at_dist_2:
            return pl.DataFrame(
                schema={
                    "group_member": pl.Utf8,
                    "neighbor": pl.Utf8,
                    "distance": pl.Int32,
                    "group_member_fitness": pl.Float64,
                    "neighbor_fitness": pl.Float64,
                    "fitness_diff": pl.Float64,
                }
            )

        # Build DataFrame from Hamming pairs and join with long_df for fitness
        pairs_df = pl.DataFrame(
            {
                "group_member": [gm for gm, _ in pairs_at_dist_2],
                "neighbor": [ext for _, ext in pairs_at_dist_2],
            }
        )

        fitness_lookup = conc_long.select("mutant_profile", "median")

        result = (
            pairs_df.join(
                fitness_lookup.rename(
                    {"mutant_profile": "group_member", "median": "group_member_fitness"}
                ),
                on="group_member",
                how="left",
            )
            .join(
                fitness_lookup.rename(
                    {"mutant_profile": "neighbor", "median": "neighbor_fitness"}
                ),
                on="neighbor",
                how="left",
            )
            .with_columns(
                (pl.col("group_member_fitness") - pl.col("neighbor_fitness")).alias(
                    "fitness_diff"
                ),
                pl.lit(2).cast(pl.Int32).alias("distance"),
            )
            .drop_nulls(subset=["group_member_fitness", "neighbor_fitness"])
            .select(
                "group_member",
                "neighbor",
                "distance",
                "group_member_fitness",
                "neighbor_fitness",
                "fitness_diff",
            )
        )

        return result

    @staticmethod
    def _compute_hamming_distances_vectorized(
        group_genotypes: list[str],
        external_genotypes: list[str],
        target_distance: int,
    ) -> list[tuple[str, str]]:
        """Compute Hamming distances using NumPy vectorization.

        Algorithm:
            1. Convert strings to 2D character arrays: [G, L] and [E, L]
            2. Broadcast comparison: [G, 1, L] != [1, E, L] → [G, E, L]
            3. Sum across positions: [G, E, L] → [G, E]
            4. Filter pairs at target distance

        Args:
            group_genotypes: List of genotype strings in the group.
            external_genotypes: List of external genotype strings.
            target_distance: Desired Hamming distance (e.g., 2).

        Returns:
            List of (group_genotype, external_genotype) pairs at target distance.
        """
        if not group_genotypes or not external_genotypes:
            return []

        group_arr = np.array([list(g) for g in group_genotypes])  # [G, L]
        ext_arr = np.array([list(e) for e in external_genotypes])  # [E, L]

        # Broadcast comparison: [G, 1, L] != [1, E, L] → [G, E, L]
        diff_matrix = group_arr[:, np.newaxis, :] != ext_arr[np.newaxis, :, :]

        # Sum to get Hamming distances: [G, E]
        hamming_matrix = diff_matrix.sum(axis=2)

        i_indices, j_indices = np.where(hamming_matrix == target_distance)

        return [
            (group_genotypes[i], external_genotypes[j])
            for i, j in zip(i_indices, j_indices, strict=True)
        ]
