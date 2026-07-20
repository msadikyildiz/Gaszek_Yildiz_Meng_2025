# This script contains the code to build the fitness landscape graph.
import contextlib
import logging
import math

# Parallel BFS for peak merging
import multiprocessing as mp
from functools import partial
from multiprocessing import get_context
from typing import Any

import networkx as nx
import numpy as np
import polars as pl

# set multiprocessing start method
# Try to set 'spawn' start method for safety in some platforms
with contextlib.suppress(RuntimeError):
    mp.set_start_method("spawn")

from fitness_landscape_graph.make_logo import (
    all_mutations,
    intended,
    mutant_dict_to_string,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GraphBuilder:
    def __init__(
        self,
        long_table: pl.DataFrame,
        pairwise_table: pl.DataFrame,
        wildtype_label: str = ".............",
    ):
        self.long_table = long_table
        self.pairwise_table = pairwise_table
        self.wildtype_label = wildtype_label

        self.initial_graph: nx.MultiDiGraph | None = None
        self.neutral_merged_graph: nx.MultiDiGraph | None = None
        self.peak_graph: nx.DiGraph | None = None

    def build_raw_graph(
        self,
        concentration: float | None = None,
        fitness_col: str = "fitness",
        fitness_diff_col: str = "fitness_diff",
        dead_mutant_label: str = "XXXXXXXXXXXXX",
        remove_dead_mutant: bool = True,
    ) -> nx.MultiDiGraph:
        """Initial raw graph: all mutants are nodes, all edges are one mutation away between mutants.

        If mutants have the same fitness, we add two edges in both directions.
        If mutants have different fitness, we add one edge in the direction of higher fitness.

        Supports global fitness graph by setting concentration to None.
        """
        if concentration is not None:
            table_nodes = self.long_table.filter(
                pl.col("concentration") == concentration
            )
            table_edges = self.pairwise_table.filter(
                pl.col("concentration") == concentration
            )
        else:
            # global fitness graph
            table_nodes = self.long_table
            table_edges = self.pairwise_table

        graph_init = nx.MultiDiGraph()

        # Add nodes
        for row in table_nodes.iter_rows(named=True):
            node_label = row["mutant_profile"]
            node_fitness = row[fitness_col]
            graph_init.add_node(
                node_label,
                fitness=node_fitness,
                group_size=1,
                group_mutants={node_label: node_fitness},
                contain_wildtype=int(node_label == self.wildtype_label),
            )

        # Add edges with weight = exp(|fitness_diff|)
        for row in table_edges.iter_rows(named=True):
            mut1 = row["mutant_profile1"]
            mut2 = row["mutant_profile2"]
            diff_val = row[fitness_diff_col]
            weight_val = math.exp(abs(diff_val))

            # If diff_val > 0 => mut2 -> mut1 (mut1 has higher fitness)
            # If diff_val < 0 => mut1 -> mut2 (mut2 has higher fitness)
            # If diff_val == 0 => direction can be arbitrary or skip
            if diff_val > 0:
                graph_init.add_edge(mut2, mut1, weight=weight_val, count=1)
            elif diff_val < 0:
                graph_init.add_edge(mut1, mut2, weight=weight_val, count=1)
            else:
                graph_init.add_edge(mut1, mut2, weight=1.0, count=1)
                graph_init.add_edge(mut2, mut1, weight=1.0, count=1)
                # logger.info(f"detected one 0 diff edge between {mut1} and {mut2}")

        # Optionally remove dead mutant
        if remove_dead_mutant and dead_mutant_label in graph_init:
            graph_init.remove_node(dead_mutant_label)

        return graph_init

    def merge_neutral_edges(
        self,
        input_graph: nx.MultiDiGraph,
        neutral_threshold: float,
        forbidden_pairs: set[frozenset] | None = None,
    ) -> nx.MultiDiGraph:
        """If edge weight <= neutral_threshold, merge the two nodes.

        Forbidden pairs are not merged.
        Aggregate edge weights as the net flow.
        """
        ############################################################
        # 1) Initialize parent map for Union-Find
        ############################################################
        parent_map: dict[str, str] = {}
        for node in input_graph.nodes():
            parent_map[node] = node

        def find_set(x: str) -> str:
            """Find the root of the set containing x."""
            if parent_map[x] != x:
                parent_map[x] = find_set(parent_map[x])
            return parent_map[x]

        def union_set(a: str, b: str):
            """Union b into a."""
            root_a, root_b = find_set(a), find_set(b)
            if root_a != root_b:
                parent_map[root_b] = root_a

        # Optional helper to check if merging u and v would unify a forbidden pair.
        # This checks if, after merging sets containing u and v, any forbidden pair
        # would end up in the same set.
        def would_merge_forbidden(u: str, v: str) -> bool:
            """Check if merging u and v would cause forbidden pair violations.

            1. Unify any forbidden pair directly
            2. Cause multiple forbidden pairs to collapse into the same edge
            3. Indirectly unify nodes that are part of different forbidden pairs
            """
            if not forbidden_pairs:
                return False  # No constraints to worry about

            ru, rv = find_set(u), find_set(v)
            if ru == rv:
                return False  # Already in the same set

            # Track all nodes involved in forbidden pairs with ru or rv
            ru_forbidden_partners = set()
            rv_forbidden_partners = set()

            for fset in forbidden_pairs:
                if len(fset) != 2:
                    logger.info(
                        f"skipping forbidden pair {fset} because it's not a pair"
                    )
                    continue

                node_a, node_b = list(fset)
                root_a, root_b = find_set(node_a), find_set(node_b)

                # Check if either representative is involved in this forbidden pair
                if root_a == ru:
                    ru_forbidden_partners.add(root_b)
                elif root_b == ru:
                    ru_forbidden_partners.add(root_a)

                if root_a == rv:
                    rv_forbidden_partners.add(root_b)
                elif root_b == rv:
                    rv_forbidden_partners.add(root_a)

            # If there's any overlap in forbidden partners, merging would collapse
            # multiple forbidden pairs into one edge
            if ru_forbidden_partners & rv_forbidden_partners:
                # logger.info(f"Prevented merging {u} and {v} to avoid collapsing multiple forbidden pairs")
                return True

            # Also check the original conditions for direct unification
            for fset in forbidden_pairs:
                if len(fset) != 2:
                    logger.info(
                        f"skipping forbidden pair {fset} because it's not a pair"
                    )
                    continue
                node_a, node_b = list(fset)
                root_a, root_b = find_set(node_a), find_set(node_b)
                if root_a in (ru, rv) and root_b in (ru, rv) and root_a != root_b:
                    return True

            return False

        ############################################################
        # 2) Union nodes if weight <= neutral_threshold, skipping forbidden merges
        ############################################################
        sorted_edges = sorted(
            input_graph.edges(data=True), key=lambda e: e[2]["weight"]
        )

        skipped_merges = set()  # Track pairs that weren't merged due to being forbidden
        for u, v, data in sorted_edges:
            w = data["weight"]
            if w > neutral_threshold:
                break

            # If this specific pair is forbidden, skip immediately
            if forbidden_pairs and frozenset({u, v}) in forbidden_pairs:
                skipped_merges.add(frozenset({u, v}))
                # logger.info(f"Skipped merging forbidden pair: {u} and {v}")
                continue

            # Otherwise, check if merging them would unify any forbidden pair indirectly
            if would_merge_forbidden(u, v):
                skipped_merges.add(frozenset({u, v}))
                # logger.info(f"Skipped merging {u} and {v} to prevent unifying forbidden pairs")
                continue

            # If it's safe to merge, union
            union_set(u, v)

        logger.info(
            f"Total pairs prevented from merging due to forbidden constraints: {len(skipped_merges)}"
        )

        ############################################################
        # 3) Build representative -> members
        # here representative nodes are results from Union-Find, so now rep nodes don't mean anything.
        ############################################################
        rep_members_map: dict[str, set[str]] = {}
        for node in input_graph.nodes():
            rep = find_set(node)
            rep_members_map.setdefault(rep, set()).add(node)

        ############################################################
        # 4) Create a new MultiDiGraph with super-nodes
        ############################################################
        merged_graph = nx.MultiDiGraph()
        for rep, members in rep_members_map.items():
            all_mutants = {}
            any_wildtype = 0
            total_size = 0
            fitness_values = []
            for m in members:
                node_data = input_graph.nodes[m]
                all_mutants.update(node_data["group_mutants"])
                any_wildtype |= node_data["contain_wildtype"]
                total_size += node_data["group_size"]
                fitness_values.append(node_data["fitness"])

            # For representative node's 'fitness', pick min or another rule
            # we will also revise this later, see rename_nodes_closest_to_average function.
            rep_fitness = min(fitness_values) if fitness_values else 0.0
            merged_graph.add_node(
                rep,
                fitness=rep_fitness,
                group_size=total_size,
                group_mutants=all_mutants,
                contain_wildtype=int(any_wildtype),
            )

        # Add logging for forbidden edges in input graph
        forbidden_edge_count = 0
        if forbidden_pairs:
            for u, v, _ in input_graph.edges(data=True):
                if frozenset({u, v}) in forbidden_pairs:
                    forbidden_edge_count += 1
        logger.info(f"Input graph has {forbidden_edge_count} forbidden edges")

        ############################################################
        # 5) Aggregate edges by net flow
        ############################################################
        edge_map: dict[tuple[str, str], dict[str, float]] = {}
        # Track forbidden pairs between representatives
        forbidden_rep_pairs = set()

        for u, v, data in input_graph.edges(data=True):
            ru, rv = find_set(u), find_set(v)
            if ru == rv:
                continue
            # If this edge is forbidden, track it at the representative level
            if forbidden_pairs and frozenset({u, v}) in forbidden_pairs:
                forbidden_rep_pairs.add(frozenset({ru, rv}))
                logger.info(f"Forbidden edge between {ru} and {rv}")

            key = (ru, rv)
            if key not in edge_map:
                edge_map[key] = {"count": 0, "total_weight": 0.0}
            edge_map[key]["count"] += data["count"]
            edge_map[key]["total_weight"] += data["weight"]

        processed_pairs = set()
        forbidden_merged_count = 0
        for (ru, rv), forward_data in edge_map.items():
            if (ru, rv) in processed_pairs or (rv, ru) in processed_pairs:
                continue

            fw_weight = forward_data["total_weight"]
            fw_count = forward_data["count"]
            reverse_data = edge_map.get((rv, ru), {"total_weight": 0.0, "count": 0})
            bw_weight = reverse_data["total_weight"]
            bw_count = reverse_data["count"]

            net_flow = fw_weight - bw_weight
            total_count = fw_count + bw_count

            # Check if this representative pair is forbidden
            is_forbidden = frozenset({ru, rv}) in forbidden_rep_pairs
            if is_forbidden:
                forbidden_merged_count += 1

            if net_flow > 0:
                merged_graph.add_edge(
                    ru,
                    rv,
                    weight=net_flow,
                    count=total_count,
                    forbidden=int(is_forbidden),
                )
            elif net_flow < 0:
                merged_graph.add_edge(
                    rv,
                    ru,
                    weight=abs(net_flow),
                    count=total_count,
                    forbidden=int(is_forbidden),
                )

            processed_pairs.add((ru, rv))
            processed_pairs.add((rv, ru))

        logger.info(f"Merged graph has {forbidden_merged_count} forbidden edges")
        return merged_graph

    @staticmethod
    def detect_peak_nodes(graph_in: nx.MultiDiGraph) -> list[str]:
        """A 'peak' node is one with out_degree=0 (no outgoing edges).

        Also requires in_degree>0.
        """
        return [
            n
            for n in graph_in.nodes()
            if graph_in.out_degree(n) == 0 and graph_in.in_degree(n) > 0
        ]

    @staticmethod
    def _gather_ancestors_worker(
        predecessor_map: dict[str, list[str]], peak: str
    ) -> tuple[str, set[str]]:
        """Parallel BFS: gather all predecessors (ancestors) of 'peak'."""
        visited = set([peak])
        ancestors = set()
        queue = [peak]

        while queue:
            current = queue.pop(0)
            for pred in predecessor_map.get(current, []):
                if pred not in visited:
                    visited.add(pred)
                    ancestors.add(pred)
                    queue.append(pred)

        # Return the union of peak with its ancestors
        return peak, ancestors.union({peak})

    @staticmethod
    def _merge_single_peak_group(
        peak_data: tuple[str, set[str], dict[str, dict[str, Any]], set[str]],
    ) -> tuple[str, dict[str, Any]]:
        """Merge all nodes in 'group_set' into 'peak', excluding 'connection_nodes'."""
        peak_node, group_set, node_data_map, connection_nodes = peak_data
        exclusive_members = group_set - connection_nodes

        combined_mutants = {}
        combined_size = 0
        combined_wildtype = 0

        for member in exclusive_members:
            member_data = node_data_map[member]
            combined_mutants.update(member_data["group_mutants"])
            combined_size += member_data["group_size"]
            combined_wildtype |= member_data["contain_wildtype"]

        peak_info = node_data_map[peak_node]
        # keep the peak's current fitness for the new super-node
        final_fitness = peak_info["fitness"]

        return peak_node, {
            "group_size": combined_size,
            "group_mutants": combined_mutants,
            "contain_wildtype": combined_wildtype,
            "is_peak": 1,
            "fitness": final_fitness,
        }

    def merge_peaks(
        self, input_graph: nx.MultiDiGraph, use_parallel: bool = True
    ) -> nx.DiGraph:
        """Peak nodes are nodes with out_degree=0 and in_degree>0.

        Connection nodes are nodes that appear in multiple peak groups.
        This function merges ancestors of each peak node, and keeps connection nodes as is.

        Output is a DiGraph as we no longer need multi-edges.
        """
        merged_graph = nx.DiGraph()

        ############################################################
        # 1) Identify all peaks
        ############################################################
        peaks = self.detect_peak_nodes(input_graph)
        logger.info(f"Found {len(peaks)} peaks.")

        ############################################################
        # 2) Build {node: [predecessors]} for BFS
        ############################################################
        predecessor_map = {}
        for n in input_graph.nodes():
            predecessor_map[n] = list(input_graph.predecessors(n))

        ############################################################
        # 3) Gather each peak's ancestor set
        ############################################################
        if use_parallel:
            with get_context("spawn").Pool() as pool:
                peak_anc_map = dict(
                    pool.map(
                        partial(self._gather_ancestors_worker, predecessor_map), peaks
                    )
                )
        else:
            peak_anc_map = {}
            for peak_node in peaks:
                peak_anc_map[peak_node], ancestors = self._gather_ancestors_worker(
                    predecessor_map, peak_node
                )
                # NOTE: The code above is slightly repetitive if single-threaded.
                # A simpler approach:
                #   pk, anc = self._gather_ancestors_worker(predecessor_map, peak_node)
                #   peak_anc_map[pk] = anc
                # but let's keep the structure consistent.

        ############################################################
        # 4) Identify "connection nodes" that appear in multiple sets
        ############################################################
        all_nodes_in_peaks = set().union(*peak_anc_map.values())
        connection_nodes = {
            x
            for x in all_nodes_in_peaks
            if sum(x in anc for anc in peak_anc_map.values()) > 1
        }

        ############################################################
        # 5) Prepare node_data map for easy merges
        ############################################################
        node_data_map: dict[str, dict[str, Any]] = {}
        for n in input_graph.nodes():
            node_data_map[n] = {
                "group_mutants": input_graph.nodes[n]["group_mutants"],
                "group_size": input_graph.nodes[n]["group_size"],
                "contain_wildtype": input_graph.nodes[n]["contain_wildtype"],
                "fitness": input_graph.nodes[n]["fitness"],
            }

        ############################################################
        # 6) Merge each peak group (excluding connection nodes)
        ############################################################
        peak_data_list = []
        for pk in peaks:
            group_set = peak_anc_map[pk]  # set of ancestors + pk
            peak_data_list.append((pk, group_set, node_data_map, connection_nodes))

        if use_parallel:
            with get_context("spawn").Pool() as pool:
                peak_supernodes = dict(
                    pool.map(self._merge_single_peak_group, peak_data_list)
                )
        else:
            peak_supernodes = {}
            for pd in peak_data_list:
                new_pk, info = self._merge_single_peak_group(pd)
                peak_supernodes[new_pk] = info

        ############################################################
        # 7) Add new peak super-nodes
        ############################################################
        for pk_node, pk_info in peak_supernodes.items():
            merged_graph.add_node(pk_node, **pk_info)

        ############################################################
        # 8) Add connection nodes individually
        ############################################################
        for cn in connection_nodes:
            cn_data = node_data_map[cn]
            merged_graph.add_node(
                cn,
                group_size=cn_data["group_size"],
                group_mutants=cn_data["group_mutants"],
                contain_wildtype=cn_data["contain_wildtype"],
                is_peak=0,
                fitness=cn_data["fitness"],
            )

        ############################################################
        # 9) Build old->new representative map
        ############################################################
        old_to_new_rep: dict[str, str] = {}
        for pk_node, anc_set in peak_anc_map.items():
            # exclude connection nodes
            for member in anc_set - connection_nodes:
                old_to_new_rep[member] = pk_node
        # each connection node maps to itself
        for cn in connection_nodes:
            old_to_new_rep[cn] = cn

        ############################################################
        # 10) Net-flow edges among final super-nodes
        ############################################################
        edge_buffer: dict[tuple[str, str], dict[str, float]] = {}
        for u, v, data in input_graph.edges(data=True):
            u_rep = old_to_new_rep[u]
            v_rep = old_to_new_rep[v]
            if u_rep == v_rep:
                continue
            key = (u_rep, v_rep)
            if key not in edge_buffer:
                edge_buffer[key] = {"count": 0.0, "total_weight": 0.0}
            edge_buffer[key]["count"] += data["count"]
            edge_buffer[key]["total_weight"] += data["weight"]

        processed_pairs = set()
        for (ru, rv), fwd_data in edge_buffer.items():
            if (ru, rv) in processed_pairs or (rv, ru) in processed_pairs:
                continue

            fw_weight = fwd_data["total_weight"]
            fw_count = fwd_data["count"]
            rev_data = edge_buffer.get((rv, ru), {"total_weight": 0.0, "count": 0.0})
            bw_weight = rev_data["total_weight"]
            bw_count = rev_data["count"]

            net_flow = fw_weight - bw_weight
            total_count = fw_count + bw_count

            if net_flow > 0:
                merged_graph.add_edge(ru, rv, weight=net_flow, count=total_count)
            elif net_flow < 0:
                merged_graph.add_edge(rv, ru, weight=abs(net_flow), count=total_count)

            processed_pairs.add((ru, rv))
            processed_pairs.add((rv, ru))

        # Optionally mark edges as "deterministic" if the target is a peak
        for source, target in merged_graph.edges():
            merged_graph[source][target]["deterministic"] = int(
                merged_graph.nodes[target].get("is_peak", 0)
            )

        return merged_graph

    def mark_peaks_and_edges(self, input_graph: nx.MultiDiGraph) -> nx.DiGraph:
        """Mark peak nodes and deterministic edges without merging ancestors.

        Similar to merge_peaks, but without merging.
        """
        marked_graph = nx.DiGraph()

        for node, data in input_graph.nodes(data=True):
            marked_graph.add_node(node, **data)
        for u, v, data in input_graph.edges(data=True):
            if not marked_graph.has_edge(u, v):
                marked_graph.add_edge(u, v, **data)

        peaks = self.detect_peak_nodes(input_graph)

        for node in marked_graph.nodes():
            marked_graph.nodes[node]["is_peak"] = int(node in peaks)

        predecessor_map = {}
        for n in input_graph.nodes():
            predecessor_map[n] = list(input_graph.predecessors(n))

        peak_anc_map = {}
        for peak_node in peaks:
            peak, ancestors = self._gather_ancestors_worker(predecessor_map, peak_node)
            peak_anc_map[peak] = ancestors

        all_nodes_in_peaks = set().union(*peak_anc_map.values())
        connection_nodes = {
            x
            for x in all_nodes_in_peaks
            if sum(x in anc for anc in peak_anc_map.values()) > 1
        }

        # Mark edges as deterministic if target is in exactly one peak group
        for u, v in marked_graph.edges():
            # If target is in a peak group but not a connection node, mark as deterministic
            marked_graph[u][v]["deterministic"] = int(
                v in all_nodes_in_peaks and v not in connection_nodes
            )

        return marked_graph

    ############################################################
    # utility functions for graph post-processing
    ############################################################
    def rename_nodes_closest_to_average(self, final_graph: nx.DiGraph) -> None:
        """Rename each node based on the mutant whose fitness is closest to the group average.

        Sets that as the node's "fitness".
        """
        rename_map = {}
        for node in final_graph.nodes():
            group_muts = final_graph.nodes[node].get("group_mutants", {})
            if not group_muts:
                continue
            avg_fit = sum(group_muts.values()) / len(group_muts)
            # find the mutant with fitness closest to avg
            best_mutant = min(group_muts, key=lambda m: abs(group_muts[m] - avg_fit))
            new_fit = group_muts[best_mutant]
            final_graph.nodes[node]["fitness"] = new_fit
            if best_mutant != node:
                rename_map[node] = best_mutant

        if rename_map:
            nx.relabel_nodes(final_graph, rename_map, copy=False)
            logger.info(
                f"Renamed {len(rename_map)} nodes to their closest-to-average mutants."
            )

    def update_peak_nodes_max_fitness(self):
        """Updates peak nodes' fitness to max fitness in their group and renames nodes accordingly."""
        nodes_to_rename = {}
        for node in self.peak_graph.nodes():
            if self.peak_graph.nodes[node].get("is_peak") != 1:
                continue

            group_mutants = self.peak_graph.nodes[node].get("group_mutants", {})
            if not group_mutants:
                logger.error(f"Peak node '{node}' has empty group_mutants")
                continue

            try:
                max_mutant = max(group_mutants, key=group_mutants.get)
                max_fitness = group_mutants[max_mutant]
                self.peak_graph.nodes[node]["fitness"] = max_fitness

                if max_mutant != node:
                    nodes_to_rename[node] = max_mutant
            except ValueError:
                logger.error(f"Peak node '{node}' has no mutants")
                continue

        if nodes_to_rename:
            self.peak_graph = nx.relabel_nodes(
                self.peak_graph, nodes_to_rename, copy=False
            )
            logger.info(f"Renamed {len(nodes_to_rename)} peak nodes")

    def get_forbidden_pairs(
        self,
        graph: nx.MultiDiGraph,
        num_forbidden_pairs: int | None = None,
        large_edge_threshold: float | None = None,
    ) -> set[frozenset]:
        """Identify forbidden pairs based on edge weights.

        First find all pairs with weight >= large_edge_threshold.
        If num_forbidden_pairs is specified, select the top k pairs from these.
        """
        # First identify all pairs above threshold
        threshold_pairs: dict[frozenset, float] = {}
        for u, v, data in graph.edges(data=True):
            w = data["weight"]
            if large_edge_threshold is not None and w >= large_edge_threshold:
                pair = frozenset({u, v})
                current_max_weight = threshold_pairs.get(pair, -1.0)
                threshold_pairs[pair] = max(current_max_weight, w)

        # Log all pairs above threshold
        sorted_threshold_pairs = sorted(
            threshold_pairs.items(), key=lambda x: x[1], reverse=True
        )
        logger.info(
            f"Found {len(sorted_threshold_pairs)} pairs with weight >= {large_edge_threshold}:"
        )
        for pair, weight in sorted_threshold_pairs:
            logger.info(f"  {sorted(pair)}: {weight:.3f}")

        # If num_forbidden_pairs specified, select top k pairs
        if num_forbidden_pairs is not None and num_forbidden_pairs > 0:
            selected_pairs = sorted_threshold_pairs[:num_forbidden_pairs]
            forbidden_pairs = {pair for pair, _ in selected_pairs}
            logger.info(f"\nSelected top {len(forbidden_pairs)} pairs as forbidden:")
            for pair, weight in selected_pairs:
                logger.info(f"  {sorted(pair)}: {weight:.3f}")
        else:
            forbidden_pairs = set(threshold_pairs.keys())

        return forbidden_pairs

    def build_graph(
        self,
        concentration: float | None = None,
        neutral_threshold: float = 0.4,
        fitness_col: str = "fitness",
        fitness_diff_col: str = "fitness_diff",
        use_parallel_peak_merge: bool = True,
        rename_by_avg_fitness: bool = True,
        large_edge_threshold: float | None = 5.5,
        num_forbidden_pairs: int | None = None,
        tiny_initial_threshold: float = 0.02,
        merge_peaks: bool = True,
        logo_threshold_ratio: float = 0.6,
    ) -> nx.DiGraph:
        logger.info(
            f"Starting build_graph, concentration={concentration}, "
            f"neutral_threshold={neutral_threshold}, tiny_initial_threshold={tiny_initial_threshold}, "
            f"large_edge_threshold={large_edge_threshold}, num_forbidden_pairs={num_forbidden_pairs}"
        )

        neutral_threshold = np.exp(neutral_threshold)
        tiny_initial_threshold = np.exp(tiny_initial_threshold)
        large_edge_threshold = np.exp(large_edge_threshold)

        logger.info(
            f"After exp -> "
            f"neutral_threshold={neutral_threshold}, "
            f"tiny_initial_threshold={tiny_initial_threshold}, "
            f"large_edge_threshold={large_edge_threshold}"
        )

        ############################################################
        # Step 1: Build raw graph
        ############################################################
        self.initial_graph = self.build_raw_graph(
            concentration=concentration,
            fitness_col=fitness_col,
            fitness_diff_col=fitness_diff_col,
        )
        logger.info(
            f"Initial raw graph -> "
            f"{self.initial_graph.number_of_nodes()} nodes, "
            f"{self.initial_graph.number_of_edges()} edges, "
            f"total edge weight={sum(data['weight'] for _, _, data in self.initial_graph.edges(data=True)):.2f}"
        )

        ############################################################
        # Step 2: Tiny neutral merging to find forbidden pairs
        ############################################################
        tiny_merged_graph = self.merge_neutral_edges(
            input_graph=self.initial_graph,
            neutral_threshold=tiny_initial_threshold,
            forbidden_pairs=None,
        )
        logger.info(
            f"After tiny merge -> "
            f"{tiny_merged_graph.number_of_nodes()} nodes, "
            f"{tiny_merged_graph.number_of_edges()} edges, "
            f"total edge weight={sum(data['weight'] for _, _, data in tiny_merged_graph.edges(data=True)):.2f}"
        )

        forbidden_pairs = self.get_forbidden_pairs(
            graph=tiny_merged_graph,
            num_forbidden_pairs=num_forbidden_pairs,
            large_edge_threshold=large_edge_threshold,
        )
        if forbidden_pairs:
            logger.info(
                f"Proceeding with {len(forbidden_pairs)} forbidden pairs for neutral merging."
            )

        ############################################################
        # Step 3: Neutral merging with forbidden pairs
        ############################################################
        self.neutral_merged_graph = self.merge_neutral_edges(
            input_graph=tiny_merged_graph,
            neutral_threshold=neutral_threshold,
            forbidden_pairs=forbidden_pairs,
        )
        logger.info(
            f"After merge neutral edges -> "
            f"{self.neutral_merged_graph.number_of_nodes()} nodes, "
            f"{self.neutral_merged_graph.number_of_edges()} edges, "
            f"total edge weight={sum(data['weight'] for _, _, data in self.neutral_merged_graph.edges(data=True)):.2f}"
        )

        ############################################################
        # Step 4: Merge peaks or just mark them
        ############################################################
        if merge_peaks:
            self.peak_graph = self.merge_peaks(
                input_graph=self.neutral_merged_graph,
                use_parallel=use_parallel_peak_merge,
            )
        else:
            self.peak_graph = self.mark_peaks_and_edges(
                input_graph=self.neutral_merged_graph
            )

        logger.info(
            f"After peak processing -> "
            f"{self.peak_graph.number_of_nodes()} nodes, "
            f"{self.peak_graph.number_of_edges()} edges, "
            f"total edge weight={sum(data['weight'] for _, _, data in self.peak_graph.edges(data=True)):.2f}"
        )

        ############################################################
        # Step 5: Rename and add logo_string
        ############################################################
        if rename_by_avg_fitness:
            self.rename_nodes_closest_to_average(self.peak_graph)

        self.update_peak_nodes_max_fitness()

        # Step 5: Add 'logo_string' attribute to each node
        for node in self.peak_graph.nodes():
            group_muts = self.peak_graph.nodes[node].get("group_mutants", {})
            logo_str = mutant_dict_to_string(
                group_muts,
                all_mutations,
                intended,
                threshold_ratio=logo_threshold_ratio,
                use_dots=True,
            )
            self.peak_graph.nodes[node]["logo_string"] = logo_str
            # Also store a 'fitness[z]' for clarity
            self.peak_graph.nodes[node]["fitness[z]"] = self.peak_graph.nodes[node][
                "fitness"
            ]

        return self.peak_graph
