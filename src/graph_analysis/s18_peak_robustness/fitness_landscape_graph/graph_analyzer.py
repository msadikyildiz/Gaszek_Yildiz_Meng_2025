"""Interactive graph analysis tool for fitness landscape graphs.

Provides a GraphAnalyzer class that loads GraphML files and supports:
- Summary of macroscopic structure (peaks, top connection nodes)
- Interactive plotly visualization using Hamming distance layout
- Node inspection and genotype search

The VisConfig dataclass allows fine-tuning of node colors, edge scaling,
arrow styles, and x-axis mode for optimal visualization.
"""

import json
import math
from collections.abc import Callable
from contextlib import suppress
from dataclasses import dataclass, field
from typing import Literal

import networkx as nx
import numpy as np
import plotly.graph_objects as go


@dataclass
class VisConfig:
    """Configuration for fitness landscape graph visualization."""

    # Node settings
    node_colors: dict[str, str] = field(
        default_factory=lambda: {
            "peak": "rgba(220,60,60,0.8)",
            "wildtype": "rgba(230,160,40,0.9)",
            "connection": "rgba(70,130,210,0.8)",
        }
    )
    node_size_base: float = 8.0
    node_size_scale: float = 5.0  # For math.log1p(group_size) * scale

    # Edge settings
    edge_styles: dict[str, dict[str, str]] = field(
        default_factory=lambda: {
            "forbidden": {
                "color": "rgba(40,180,40,0.7)",
                "dash": "dash",
                "name": "Forbidden edge",
            },
            "deterministic": {
                "color": "rgba(220,60,60,0.6)",
                "dash": "solid",
                "name": "Deterministic edge",
            },
            "other": {
                "color": "rgba(70,130,210,0.35)",
                "dash": "solid",
                "name": "Other edge",
            },
        }
    )
    edge_width_scaling: Callable[[float], float] = lambda w: min(
        8.0, max(1.0, math.log1p(w) * 0.5)
    )  # Exposed for experimentation

    # Arrow settings (optional)
    enable_arrows: bool = True
    arrow_size: int = 6  # Increased from 4
    arrow_symbol: str = "triangle-up"

    # Hamming layout settings
    x_mode: Literal["node", "average"] = "average"
    # "node": Use node ID's Hamming distance from wildtype
    # "average": Use average Hamming distance of group_mutants

    # Other
    highlight_border_width: int = 3
    hover_template: str | None = None  # Custom hover text if needed
    figure_width: int = 1400
    figure_height: int = 1000


class GraphAnalyzer:
    """Analyze and visualize a fitness landscape graph from a GraphML file."""

    def __init__(self, graphml_path: str) -> None:
        self.path = graphml_path
        self.graph = self._load_graph(graphml_path)
        self._peaks: list[tuple[str, dict]] = []
        self._connections: list[tuple[str, dict]] = []
        self._classify_nodes()

    @staticmethod
    def _load_graph(graphml_path: str) -> nx.DiGraph:
        """Load GraphML and deserialize JSON-encoded attributes."""
        graph = nx.read_graphml(graphml_path)
        for node in graph.nodes():
            data = graph.nodes[node]
            if isinstance(data.get("group_mutants"), str):
                data["group_mutants"] = json.loads(data["group_mutants"])
            for key in ("fitness", "group_size", "is_peak", "contain_wildtype"):
                if key in data and isinstance(data[key], str):
                    with suppress(ValueError):
                        data[key] = float(data[key])
            if "group_size" in data:
                data["group_size"] = int(data["group_size"])
            if "is_peak" in data:
                data["is_peak"] = int(data["is_peak"])
            if "contain_wildtype" in data:
                data["contain_wildtype"] = int(data["contain_wildtype"])
        for u, v in graph.edges():
            edata = graph[u][v]
            for key in ("weight",):
                if key in edata and isinstance(edata[key], str):
                    edata[key] = float(edata[key])
            for key in ("deterministic", "forbidden", "count"):
                if key in edata and isinstance(edata[key], str):
                    edata[key] = int(float(edata[key]))
        return graph

    @staticmethod
    def _hamming_distance(a: str, b: str) -> int:
        """Compute Hamming distance between equal-length strings."""
        if len(a) != len(b):
            raise ValueError(
                f"Unequal lengths for Hamming distance: {len(a)} != {len(b)}"
            )
        return sum(1 for ca, cb in zip(a, b, strict=True) if ca != cb)

    def _wildtype_profile(self) -> str:
        """Return the wildtype genotype profile string."""
        wildtype_default = "............."
        wildtype_node = None
        for n, d in self.graph.nodes(data=True):
            if d.get("contain_wildtype") == 1:
                wildtype_node = n
                break
        if wildtype_node is None:
            return wildtype_default
        group_mutants = self.graph.nodes[wildtype_node].get("group_mutants", {})
        if isinstance(group_mutants, dict) and wildtype_default in group_mutants:
            return wildtype_default
        if isinstance(wildtype_node, str) and len(wildtype_node) == len(
            wildtype_default
        ):
            return wildtype_node
        return wildtype_default

    def _classify_nodes(self) -> None:
        """Sort nodes into peaks, wildtype, and connections.

        Wildtype nodes (contain_wildtype == 1) are always classified as
        wildtype, even if they are also peaks.
        """
        self._peaks = []
        self._wildtype = []
        self._connections = []
        for n, d in self.graph.nodes(data=True):
            if d.get("contain_wildtype") == 1:
                self._wildtype.append((n, d))
            elif d.get("is_peak") == 1:
                self._peaks.append((n, d))
            else:
                self._connections.append((n, d))
        self._peaks.sort(key=lambda x: x[1]["fitness"], reverse=True)
        self._wildtype.sort(key=lambda x: x[1]["fitness"], reverse=True)
        self._connections.sort(key=lambda x: x[1]["group_size"], reverse=True)

    def summary(self, top_k: int = 10) -> None:
        """Print overview of the graph structure.

        Args:
            top_k: Number of top connection nodes to display.
        """
        g = self.graph
        print(f"Graph: {self.path}")
        print(f"Nodes: {g.number_of_nodes()}  Edges: {g.number_of_edges()}")
        print(f"Peaks: {len(self._peaks)}  Connection nodes: {len(self._connections)}")
        print()

        print("=== Peak nodes (sorted by fitness) ===")
        print(f"  {'Node':<22s} {'Fitness':>8s} {'Size':>6s}  Logo")
        print(f"  {'-' * 22} {'-' * 8} {'-' * 6}  {'-' * 30}")
        for n, d in self._peaks:
            print(
                f"  {n:<22s} {d['fitness']:>8.3f} {d['group_size']:>6d}  "
                f"{d.get('logo_string', '')}"
            )

        print()
        print(f"=== Top {top_k} connection nodes (sorted by group_size) ===")
        print(f"  {'Node':<22s} {'Fitness':>8s} {'Size':>6s}  Logo")
        print(f"  {'-' * 22} {'-' * 8} {'-' * 6}  {'-' * 30}")
        for n, d in self._connections[:top_k]:
            print(
                f"  {n:<22s} {d['fitness']:>8.3f} {d['group_size']:>6d}  "
                f"{d.get('logo_string', '')}"
            )

    # ===================================================================
    # Layout computation
    # ===================================================================

    def _compute_layout(
        self, g: nx.DiGraph, config: VisConfig
    ) -> dict[str, tuple[float, float]]:
        """Compute node positions using Hamming distance for x, fitness for y.

        Args:
            g: The graph to layout.
            config: Visualization configuration.

        Returns:
            Dictionary mapping node labels to (x, y) positions.
        """
        wildtype_profile = self._wildtype_profile()

        pos: dict[str, tuple[float, float]] = {}
        for n, d in g.nodes(data=True):
            if d.get("contain_wildtype") == 1:
                x_val = 0.0
            else:
                if config.x_mode == "node":
                    if isinstance(n, str) and len(n) == len(wildtype_profile):
                        x_val = float(self._hamming_distance(n, wildtype_profile))
                    else:
                        x_val = 0.0
                elif config.x_mode == "average":
                    group_mutants = d.get("group_mutants", {})
                    if isinstance(group_mutants, dict) and group_mutants:
                        distances = [
                            self._hamming_distance(m, wildtype_profile)
                            for m in group_mutants
                            if isinstance(m, str) and len(m) == len(wildtype_profile)
                        ]
                        x_val = float(np.mean(distances)) if distances else 0.0
                    else:
                        x_val = 0.0
                else:
                    raise ValueError("Unknown x_mode. Use 'node' or 'average'.")

            pos[n] = (x_val, float(d["fitness"]))

        return pos

    @staticmethod
    def _resolve_overlaps(
        pos: dict[str, tuple[float, float]],
        node_sizes: dict[str, float],
        iterations: int = 50,
        scale: float = 0.02,
    ) -> dict[str, tuple[float, float]]:
        """Push overlapping nodes apart along x-axis, keeping y (fitness) fixed.

        Args:
            pos: Node positions {node: (x, y)}.
            node_sizes: Node marker sizes {node: size_px}.
            iterations: Number of repulsion passes.
            scale: Coordinate-space scale factor for marker sizes.
        """
        nodes = list(pos)
        if len(nodes) <= 1:
            return pos
        xs = np.array([pos[n][0] for n in nodes], dtype=np.float64)
        ys = np.array([pos[n][1] for n in nodes], dtype=np.float64)
        sizes = np.array([node_sizes.get(n, 8.0) for n in nodes], dtype=np.float64)

        for _ in range(iterations):
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    dx = xs[j] - xs[i]
                    dy = ys[j] - ys[i]
                    d = math.sqrt(dx * dx + dy * dy)
                    min_dist = (sizes[i] + sizes[j]) * scale
                    if d < min_dist and d > 1e-9:
                        push = (min_dist - d) / d * 0.3
                        xs[i] -= push * dx
                        xs[j] += push * dx

        return {n: (float(xs[i]), float(ys[i])) for i, n in enumerate(nodes)}

    # ===================================================================
    # Visualization helpers
    # ===================================================================

    def _add_edge_traces(
        self, fig: go.Figure, pos: dict[str, tuple[float, float]], config: VisConfig
    ) -> None:
        """Add edge traces to the figure."""
        g = self.graph
        # --- Edge traces (3 categories) ---
        edge_groups: dict[str, list[tuple[str, str, dict]]] = {
            "forbidden": [],
            "deterministic": [],
            "other": [],
        }
        for u, v in g.edges():
            edata = g[u][v]
            if edata.get("forbidden", 0) == 1:
                edge_groups["forbidden"].append((u, v, edata))
            elif edata.get("deterministic", 0) == 1:
                edge_groups["deterministic"].append((u, v, edata))
            else:
                edge_groups["other"].append((u, v, edata))

        # Draw edges per category
        arrow_x, arrow_y, arrow_colors = [], [], []
        hover_x, hover_y, hover_text = [], [], []

        for cat, edges in edge_groups.items():
            if not edges:
                continue
            style = config.edge_styles[cat]
            ex, ey = [], []
            for u, v, edata in edges:
                x0, y0 = pos[u]
                x1, y1 = pos[v]
                w = edata.get("weight", 1.0)
                ex.extend([x0, x1, None])
                ey.extend([y0, y1, None])

                if config.enable_arrows:
                    # Arrow marker at ~80% along edge
                    ax = x0 + 0.8 * (x1 - x0)
                    ay = y0 + 0.8 * (y1 - y0)
                    arrow_x.append(ax)
                    arrow_y.append(ay)
                    arrow_colors.append(style["color"])

                # Hover at midpoint
                mx = (x0 + x1) / 2
                my = (y0 + y1) / 2
                hover_x.append(mx)
                hover_y.append(my)
                hover_text.append(
                    f"<b>{u} → {v}</b><br>"
                    f"Weight: {w:.2f}<br>"
                    f"Deterministic: {edata.get('deterministic', '?')}<br>"
                    f"Forbidden: {edata.get('forbidden', 0)}"
                )

            # Use median width for the category trace
            all_weights = [e[2].get("weight", 1.0) for e in edges]
            median_w = float(np.median(all_weights))
            cat_width = config.edge_width_scaling(median_w)

            fig.add_trace(
                go.Scatter(
                    x=ex,
                    y=ey,
                    mode="lines",
                    name=style["name"],
                    line={
                        "width": cat_width,
                        "color": style["color"],
                        "dash": style["dash"],
                    },
                    hoverinfo="none",
                    showlegend=True,
                )
            )

        # Arrow markers
        if config.enable_arrows and arrow_x:
            fig.add_trace(
                go.Scatter(
                    x=arrow_x,
                    y=arrow_y,
                    mode="markers",
                    marker={
                        "symbol": config.arrow_symbol,
                        "size": config.arrow_size,
                        "color": arrow_colors,
                        "opacity": 0.6,
                    },
                    hoverinfo="none",
                    showlegend=False,
                )
            )

        # Edge hover (invisible midpoint markers)
        if hover_x:
            fig.add_trace(
                go.Scatter(
                    x=hover_x,
                    y=hover_y,
                    mode="markers",
                    marker={"size": 8, "color": "rgba(0,0,0,0)", "opacity": 0},
                    text=hover_text,
                    hoverinfo="text",
                    showlegend=False,
                )
            )

    def _add_node_traces(
        self,
        fig: go.Figure,
        pos: dict[str, tuple[float, float]],
        node_sizes: dict[str, float],
        highlight_set: set[str],
        config: VisConfig,
    ) -> None:
        """Add node traces to the figure."""
        g = self.graph

        conn_nodes = [n for n, _ in self._connections]
        wt_nodes = [n for n, _ in self._wildtype]
        peak_nodes = [n for n, _ in self._peaks]

        def _make_node_trace(
            nodes: list[str], color: str, name: str, symbol: str = "circle"
        ) -> go.Scatter:
            x_vals, y_vals, sizes, hovers, borders = [], [], [], [], []
            for n in nodes:
                d = g.nodes[n]
                x_vals.append(pos[n][0])
                y_vals.append(pos[n][1])
                sizes.append(node_sizes[n])
                hover = (
                    f"<b>{n}</b><br>"
                    f"Fitness: {d['fitness']:.3f}<br>"
                    f"Group size: {d['group_size']}<br>"
                    f"Peak: {'Yes' if d.get('is_peak') == 1 else 'No'}<br>"
                    f"WT: {'Yes' if d.get('contain_wildtype') == 1 else 'No'}<br>"
                    f"Logo: {d.get('logo_string', 'N/A')}"
                )
                hovers.append(hover)
                borders.append("green" if n in highlight_set else "rgba(0,0,0,0.3)")
            border_width = [
                config.highlight_border_width if n in highlight_set else 1
                for n in nodes
            ]
            return go.Scatter(
                x=x_vals,
                y=y_vals,
                mode="markers",
                name=name,
                marker={
                    "size": sizes,
                    "color": color,
                    "symbol": symbol,
                    "line": {"width": border_width, "color": borders},
                    "opacity": 0.85,
                },
                text=hovers,
                hoverinfo="text",
            )

        fig.add_trace(
            _make_node_trace(conn_nodes, config.node_colors["connection"], "Connection")
        )
        fig.add_trace(
            _make_node_trace(wt_nodes, config.node_colors["wildtype"], "Wildtype")
        )
        fig.add_trace(_make_node_trace(peak_nodes, config.node_colors["peak"], "Peak"))

    # ===================================================================
    # Public API methods
    # ===================================================================

    def plot(
        self,
        highlight_nodes: list[str] | None = None,
        width: int | None = None,
        height: int | None = None,
        config: VisConfig | None = None,
    ) -> go.Figure:
        """Create interactive plotly visualization of the graph.

        Uses Hamming distance from wildtype for x-axis positioning.

        Args:
            highlight_nodes: Nodes to highlight with a green border.
            width: Figure width in pixels.
            height: Figure height in pixels.
            config: Visualization configuration. Defaults to VisConfig().

        Returns:
            Plotly Figure object.
        """
        if config is None:
            config = VisConfig()
        if width is None:
            width = config.figure_width
        if height is None:
            height = config.figure_height
        g = self.graph
        highlight_set = set(highlight_nodes or [])

        # --- Compute layout ---
        pos = self._compute_layout(g, config)

        # --- Compute node sizes (needed for overlap resolution) ---
        node_sizes: dict[str, float] = {}
        for n, d in g.nodes(data=True):
            node_sizes[n] = max(
                config.node_size_base,
                config.node_size_scale * math.log1p(d["group_size"]),
            )

        # --- Resolve overlaps ---
        pos = self._resolve_overlaps(pos, node_sizes)

        fig = go.Figure()

        self._add_edge_traces(fig, pos, config)
        self._add_node_traces(fig, pos, node_sizes, highlight_set, config)

        fig.update_layout(
            width=width,
            height=height,
            showlegend=True,
            xaxis={"visible": True, "title": "Hamming distance from WT"},
            yaxis={"title": "Fitness", "side": "left"},
            plot_bgcolor="white",
            hovermode="closest",
            margin={"l": 60, "r": 20, "t": 40, "b": 40},
        )

        return fig

    def inspect(self, node: str) -> None:
        """Print detailed information about a node.

        Args:
            node: Node label to inspect.
        """
        if node not in self.graph:
            print(f"Node '{node}' not found in graph.")
            return

        d = self.graph.nodes[node]
        print(f"Node: {node}")
        print(f"  Fitness:          {d['fitness']:.4f}")
        print(f"  Group size:       {d['group_size']}")
        print(f"  Is peak:          {'Yes' if d.get('is_peak') == 1 else 'No'}")
        print(f"  Contains WT:      {'Yes' if d.get('contain_wildtype') else 'No'}")
        print(f"  Logo:             {d.get('logo_string', 'N/A')}")

        # Predecessors (nodes pointing TO this node)
        preds = list(self.graph.predecessors(node))
        if preds:
            print(f"  Predecessors ({len(preds)}):")
            for p in sorted(
                preds, key=lambda n: self.graph.nodes[n]["fitness"], reverse=True
            ):
                w = self.graph[p][node].get("weight", "?")
                pd = self.graph.nodes[p]
                print(
                    f"    {p:<22s}  fit={pd['fitness']:.3f}  "
                    f"size={pd['group_size']:>5d}  w={w}"
                )

        # Successors (nodes this node points TO)
        succs = list(self.graph.successors(node))
        if succs:
            print(f"  Successors ({len(succs)}):")
            for s in sorted(
                succs, key=lambda n: self.graph.nodes[n]["fitness"], reverse=True
            ):
                w = self.graph[node][s].get("weight", "?")
                sd = self.graph.nodes[s]
                print(
                    f"    {s:<22s}  fit={sd['fitness']:.3f}  "
                    f"size={sd['group_size']:>5d}  w={w}"
                )

        if not preds and not succs:
            print("  No neighbors (isolated node).")

    def get_group_mutants(self, node: str) -> dict[str, float]:
        """Return the group_mutants dict for a node.

        Args:
            node: Node label.

        Returns:
            Dict mapping genotype strings to fitness values.
        """
        if node not in self.graph:
            raise KeyError(f"Node '{node}' not found in graph.")
        return self.graph.nodes[node].get("group_mutants", {})

    def find_genotype(self, genotype: str) -> str | None:
        """Find which node contains a specific genotype string.

        Args:
            genotype: A 13-character mutant profile string.

        Returns:
            The node label containing the genotype, or None.
        """
        for node, data in self.graph.nodes(data=True):
            if genotype in data.get("group_mutants", {}):
                return node
        return None

    def get_node_genotypes(self, node_id: str) -> set[str]:
        """Get all genotype strings from a supernode.

        Args:
            node_id: Node identifier in the graph.

        Returns:
            Set of genotype strings in the node's group_mutants.

        Raises:
            KeyError: If node_id not found in graph.
        """
        if node_id not in self.graph:
            raise KeyError(f"Node '{node_id}' not found in graph.")
        return set(self.graph.nodes[node_id].get("group_mutants", {}).keys())

    def get_nodes_genotypes(self, node_ids: list[str]) -> set[str]:
        """Get all genotype strings from multiple supernodes.

        Args:
            node_ids: List of node identifiers in the graph.

        Returns:
            Set of genotype strings from all specified nodes (union).
        """
        all_genotypes = set()
        for node_id in node_ids:
            all_genotypes.update(self.get_node_genotypes(node_id))
        return all_genotypes

    def get_peak_genotypes(self, rank: int = 0) -> set[str]:
        """Get all genotype strings from a peak supernode.

        Convenience method for common case. For custom selection,
        use get_node_genotypes() or get_nodes_genotypes() instead.

        Args:
            rank: 0 for highest-fitness peak, 1 for second, etc.

        Returns:
            Set of genotype strings in the peak's group_mutants.
        """
        peaks = [
            (n, d) for n, d in self.graph.nodes(data=True) if d.get("is_peak") == 1
        ]
        peaks.sort(key=lambda x: x[1]["fitness"], reverse=True)
        if rank >= len(peaks):
            return set()
        node_id, _ = peaks[rank]
        return self.get_node_genotypes(node_id)

    def has_global_peak(
        self, min_group_size: int = 32
    ) -> tuple[bool, str | None, dict]:
        """Detect if graph has a global peak.

        A global peak must:
        1. Be a peak node (is_peak=1)
        2. Have group_size > min_group_size
        3. Have the highest fitness in the graph

        Args:
            min_group_size: Minimum group size to qualify as global peak.

        Returns:
            Tuple of (has_global_peak, peak_node_id, peak_info) where:
            - has_global_peak: True if a global peak exists
            - peak_node_id: Node ID of the global peak (None if doesn't exist)
            - peak_info: Dict with 'fitness' and 'group_size' (empty if doesn't exist)
        """
        # Find the highest fitness node in the entire graph
        max_fitness = -float("inf")
        max_fitness_node = None
        for n, d in self.graph.nodes(data=True):
            if d["fitness"] > max_fitness:
                max_fitness = d["fitness"]
                max_fitness_node = n

        if max_fitness_node is None:
            return False, None, {}

        # Check if the highest fitness node is a large peak
        node_data = self.graph.nodes[max_fitness_node]
        is_peak = node_data.get("is_peak", 0) == 1
        group_size = node_data.get("group_size", 0)
        is_large = group_size > min_group_size

        if is_peak and is_large:
            peak_info = {
                "fitness": node_data["fitness"],
                "group_size": group_size,
            }
            return True, max_fitness_node, peak_info

        return False, None, {}

    @staticmethod
    def hamming_distance(a: str, b: str) -> int:
        """Compute Hamming distance between equal-length strings.

        Public static method for external use (replaces _hamming_distance).

        Args:
            a: First string.
            b: Second string.

        Returns:
            Number of positions where characters differ.

        Raises:
            ValueError: If strings have unequal lengths.
        """
        if len(a) != len(b):
            raise ValueError(f"Unequal lengths: {len(a)} != {len(b)}")
        return sum(1 for ca, cb in zip(a, b, strict=True) if ca != cb)
