"""Altair theme, cluster legend, label nudging, and connectors."""

import numpy as np
import altair as alt

# ── Fonts ──
FONT_BODY = "Avenir Next"
FONT_MONO = "Menlo"


def register_theme():
    """Register and enable the TEM-1 Altair theme."""

    @alt.theme.register("tem1", enable=True)
    def tem1_theme():
        return alt.theme.ThemeConfig({
            "config": {
                "font": FONT_BODY,
                "background": "#FAFAFA",
                "view": {"stroke": None},
                "autosize": {"type": "pad"},
                "axis": {
                    "labelFont": FONT_BODY,
                    "titleFont": FONT_BODY,
                    "labelFontSize": 11,
                    "titleFontSize": 13,
                    "gridColor": "#E8E8E8",
                    "domainColor": "#424242",
                    "tickColor": "#424242",
                },
                "header": {
                    "labelFont": FONT_BODY,
                    "titleFont": FONT_BODY,
                    "labelFontSize": 12,
                    "titleFontSize": 14,
                },
                "legend": {
                    "labelFont": FONT_MONO,
                    "titleFont": FONT_BODY,
                    "labelFontSize": 11,
                    "titleFontSize": 12,
                },
                "title": {
                    "font": FONT_BODY,
                    "fontSize": 15,
                    "anchor": "start",
                },
            }
        })

    alt.data_transformers.disable_max_rows()


_ROW_PX = 16  # pixels per legend item


def build_cluster_col(
    cluster_name: str,
    members: list[str],
    *,
    short_colors: dict[str, str],
    geno_lookup: dict[str, str],
    wt_genotype: str,
    mutation_cues: dict[str, str],
) -> alt.LayerChart:
    """Single cluster column for the legend with padded tabular text."""
    import pandas as pd

    _NAME_W = 4
    _GENO_W = 13

    geno_strs: list[str] = []
    for m in members:
        g = geno_lookup.get(m)
        if g is None:
            geno_strs.append("")
        elif m == "WT":
            geno_strs.append(g)
        else:
            geno_strs.append(
                "".join("." if c == w else c for c, w in zip(g, wt_genotype))
            )

    legend_text: list[str] = []
    for m, gs in zip(members, geno_strs):
        cue = mutation_cues.get(m, "")
        gs_padded = gs.ljust(_GENO_W) if gs else " " * _GENO_W
        cue_part = f"  {cue}" if cue and cue != "base" else ""
        legend_text.append(f"{m:<{_NAME_W}s}  {gs_padded}{cue_part}")

    n = len(members)
    height = max(n * _ROW_PX, 20)

    items = pd.DataFrame({
        "name": legend_text,
        "color": [short_colors[m] for m in members],
        "y_px": [i * _ROW_PX + _ROW_PX // 2 for i in range(n)],
        "x": [0] * n,
    })

    dots = (
        alt.Chart(items)
        .mark_point(size=45, filled=True)
        .encode(
            x=alt.X("x:Q", axis=None, scale=alt.Scale(domain=[0, 250])),
            y=alt.Y("y_px:Q", axis=None, scale=alt.Scale(domain=[0, height])),
            color=alt.Color("color:N", scale=None),
        )
    )
    labels = (
        alt.Chart(items)
        .mark_text(
            align="left", dx=10, fontSize=10, font=FONT_MONO,
            color="#424242",
        )
        .encode(
            x=alt.X("x:Q", axis=None, scale=alt.Scale(domain=[0, 250])),
            y=alt.Y("y_px:Q", axis=None, scale=alt.Scale(domain=[0, height])),
            text="name:N",
        )
    )
    return (
        (dots + labels)
        .properties(
            width=250, height=height,
            title=alt.Title(
                cluster_name, fontSize=11, font=FONT_BODY,
                anchor="start", offset=4,
            ),
        )
    )


# ── Label nudging + connectors ──

def nudge_labels(
    x: np.ndarray,
    y: np.ndarray,
    x_range: float,
    y_range: float,
    all_pts_x: np.ndarray | None = None,
    all_pts_y: np.ndarray | None = None,
    dx_init: float = 0.04,
    dy_init: float = -0.03,
    passes: int = 12,
    repel_radius: float = 0.06,
    point_repel_radius: float = 0.04,
    pad: float = 0.02,
) -> tuple[np.ndarray, np.ndarray]:
    """Iterative label repulsion (label-label + label-point + boundary clamp)."""
    lx = x + dx_init * x_range
    ly = y + dy_init * y_range
    if all_pts_x is None:
        all_pts_x = x
    if all_pts_y is None:
        all_pts_y = y

    x_min = x.min() - pad * x_range
    x_max = x.max() + pad * x_range
    y_min = y.min() - pad * y_range
    y_max = y.max() + pad * y_range

    for _ in range(passes):
        # Label-label repulsion
        for i in range(len(lx)):
            for j in range(i + 1, len(lx)):
                dx = (lx[i] - lx[j]) / x_range
                dy = (ly[i] - ly[j]) / y_range
                dist = np.sqrt(dx**2 + dy**2)
                if 0 < dist < repel_radius:
                    force = (repel_radius - dist) / 2
                    fx = force * dx / dist * x_range
                    fy = force * dy / dist * y_range
                    lx[i] += fx
                    ly[i] += fy
                    lx[j] -= fx
                    ly[j] -= fy

        # Label-point repulsion (labels pushed away from ALL dots)
        for i in range(len(lx)):
            for px, py in zip(all_pts_x, all_pts_y):
                dx = (lx[i] - px) / x_range
                dy = (ly[i] - py) / y_range
                dist = np.sqrt(dx**2 + dy**2)
                if 0 < dist < point_repel_radius:
                    force = (point_repel_radius - dist) * 0.5
                    lx[i] += force * dx / dist * x_range
                    ly[i] += force * dy / dist * y_range

        # Boundary clamp
        lx = np.clip(lx, x_min, x_max)
        ly = np.clip(ly, y_min, y_max)

    return lx, ly


def make_connectors(
    data: "pd.DataFrame",
    x_col: str,
    y_col: str,
    lx_col: str = "_lx",
    ly_col: str = "_ly",
    min_dist: float = 0.02,
    x_range: float = 1.0,
    y_range: float = 1.0,
) -> alt.Chart | None:
    """Thin grey connector lines from data points to nudged label positions."""
    import pandas as _pd

    rows = []
    for _, row in data.iterrows():
        dx = (row[lx_col] - row[x_col]) / x_range
        dy = (row[ly_col] - row[y_col]) / y_range
        if np.sqrt(dx**2 + dy**2) < min_dist:
            continue
        rows.append({"x": row[x_col], "y": row[y_col], "x2": row[lx_col], "y2": row[ly_col]})

    if not rows:
        return None

    conn_df = _pd.DataFrame(rows)
    return (
        alt.Chart(conn_df)
        .mark_rule(strokeWidth=0.5, color="#BDBDBD", opacity=0.5)
        .encode(x="x:Q", y="y:Q", x2="x2:Q", y2="y2:Q")
    )
