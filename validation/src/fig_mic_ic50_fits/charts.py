"""Chart generators for the unified MIC/IC50 fits figure."""

import numpy as np
import pandas as pd
import altair as alt
from scipy.optimize import curve_fit, brentq

from config import FONT_BODY, FONT_MONO, MIC_THRESHOLD

GREY = "#616161"
ACCENT = "#E63946"
IC50_COLOR = "#2E7D32"
_AX = {"labelFontSize": 8, "titleFontSize": 9}


def _fmt_conc(val: float) -> str:
    """Format a concentration value with appropriate precision."""
    if val >= 10:
        return f"{val:.0f}"
    if val >= 1:
        return f"{val:.1f}"
    if val >= 0.1:
        return f"{val:.2f}"
    return f"{val:.2g}"

# ── Hill model ──


def _hill(x, bottom, top, ic50, hill):
    return bottom + (top - bottom) / (1.0 + (x / ic50) ** hill)


# ── Growth grid: small-multiples of all concentrations ──


def growth_grid(gc_df: pd.DataFrame, color: str, n_cols: int = 4):
    """5×2 grid of small panels, one per concentration, all reps shown.

    Concentration label is placed inside each panel (top-left text overlay).
    """
    concs = sorted(gc_df["concentration"].unique())

    x_max = float(gc_df["time_h"].max()) * 1.05
    y_max = max(float(gc_df["od_bgsub"].max()) * 1.15, 0.01)

    panels = []
    for conc in concs:
        sub = gc_df[gc_df["concentration"] == conc].copy()
        if conc == 0:
            conc_str = "ctrl"
        elif conc >= 1:
            conc_str = f"{conc:.0f}"
        else:
            conc_str = f"{conc:.1f}"

        lines = (
            alt.Chart(sub)
            .mark_line(strokeWidth=0.8, opacity=0.5, color=color)
            .encode(
                x=alt.X(
                    "time_h:Q", axis=alt.Axis(**_AX, tickCount=3, grid=False),
                    scale=alt.Scale(domain=[0, x_max]), title=None,
                ),
                y=alt.Y(
                    "od_bgsub:Q", axis=alt.Axis(**_AX, tickCount=3),
                    scale=alt.Scale(domain=[0, y_max]), title=None,
                ),
                detail="rep_id:N",
            )
        )
        # Concentration label inside panel (top-left)
        lbl = (
            alt.Chart(pd.DataFrame({"x": [x_max * 0.04], "y": [y_max * 0.92],
                                     "t": [conc_str]}))
            .mark_text(
                align="left", fontSize=7, font=FONT_MONO,
                color="#757575", fontWeight="bold",
            )
            .encode(x="x:Q", y="y:Q", text="t:N")
        )

        chart = (
            (lines + lbl)
            .properties(width=68, height=48)
        )
        panels.append(chart)

    rows = []
    for i in range(0, len(panels), n_cols):
        rows.append(alt.hconcat(*panels[i : i + n_cols]).properties(spacing=2))
    return alt.vconcat(*rows).properties(spacing=2)


# ── Dose-response + Hill fit panel ──


def dose_response_hill(
    dr_df: pd.DataFrame,
    color: str,
    norm_col: str = "norm_auc",
    y_title: str = "Normalized AUC",
):
    """Dose-response scatter + Hill fit + MIC90/IC50 annotations.

    Parameters
    ----------
    dr_df : pandas DataFrame filtered to one variant+drug.
            Must have: concentration, {norm_col}_mean, {norm_col}_std,
            {norm_col}_lo, {norm_col}_hi
    color : hex color for data points
    norm_col : column prefix (norm_auc or norm_od)
    y_title : y-axis label

    Returns
    -------
    (chart, ic50_fit, mic_display)
    """
    mean_col = f"{norm_col}_mean"
    lo_col = f"{norm_col}_lo"
    hi_col = f"{norm_col}_hi"

    sub = dr_df[dr_df["concentration"] > 0].sort_values("concentration").copy()
    c_arr = sub["concentration"].values
    log_range = np.log10(c_arr)
    x_lo = float(log_range.min()) - 0.2
    x_hi = float(log_range.max()) + 0.2
    sub["log_conc"] = np.log10(sub["concentration"])

    scatter = (
        alt.Chart(sub)
        .mark_point(size=25, filled=True, color=color)
        .encode(
            x=alt.X(
                "log_conc:Q", title="log₁₀ Conc (µg/mL)",
                axis=alt.Axis(**_AX),
                scale=alt.Scale(domain=[x_lo, x_hi]),
            ),
            y=alt.Y(
                f"{mean_col}:Q", title=y_title,
                axis=alt.Axis(**_AX),
                scale=alt.Scale(domain=[-0.05, 1.25]),
            ),
        )
    )
    sub["lo"], sub["hi"] = sub[lo_col], sub[hi_col]
    errbar = (
        alt.Chart(sub)
        .mark_errorbar(color=color)
        .encode(x="log_conc:Q", y=alt.Y("lo:Q", title=""), y2="hi:Q")
    )
    layers = [errbar, scatter]

    y_arr = np.clip(sub[mean_col].values, 0, 2.0)
    min_c, max_c = float(c_arr.min()), float(c_arr.max())
    ic50_fit, mic_fit, mic_display = None, None, None

    try:
        popt, _ = curve_fit(
            _hill, c_arr, y_arr,
            p0=[0.0, 1.0, np.sqrt(min_c * max_c), 1.0],
            bounds=([0, 0.5, min_c / 10, 0.1], [0.5, 1.5, max_c * 10, 10]),
            maxfev=5000,
        )
        ic50_fit = float(popt[2])
        fit_x = np.logspace(np.log10(min_c), np.log10(max_c), 100)
        fit_df = pd.DataFrame({"lc": np.log10(fit_x), "y": _hill(fit_x, *popt)})
        layers.append(
            alt.Chart(fit_df)
            .mark_line(strokeWidth=1.2, strokeDash=[4, 3], color=GREY)
            .encode(x="lc:Q", y="y:Q")
        )
        if float(popt[0]) <= MIC_THRESHOLD:
            try:
                mic_fit = brentq(
                    lambda x: _hill(x, *popt) - MIC_THRESHOLD, min_c, max_c,
                )
                mic_display = _fmt_conc(mic_fit)
            except ValueError:
                mic_display = f"> {_fmt_conc(max_c)}"
        else:
            mic_display = f"> {_fmt_conc(max_c)}"
    except (RuntimeError, ValueError):
        if len(c_arr) > 1 and y_arr.min() < MIC_THRESHOLD:
            mic_interp = float(np.interp(
                MIC_THRESHOLD, y_arr[::-1], c_arr[::-1],
            ))
            mic_display = f"~{_fmt_conc(mic_interp)}"
        else:
            mic_display = f"> {_fmt_conc(max_c)}"

    # MIC threshold line
    layers.append(
        alt.Chart(pd.DataFrame({"y": [MIC_THRESHOLD]}))
        .mark_rule(strokeDash=[4, 3], strokeWidth=0.8, color=ACCENT)
        .encode(y="y:Q")
    )

    # IC50 horizontal line at 0.5
    if ic50_fit is not None:
        layers.append(
            alt.Chart(pd.DataFrame({"y": [0.5]}))
            .mark_rule(strokeDash=[6, 4], strokeWidth=0.6, color=IC50_COLOR)
            .encode(y="y:Q")
        )

    # MIC vertical rule
    if mic_fit is not None:
        log_mic = np.log10(mic_fit)
        if x_lo <= log_mic <= x_hi:
            layers.append(
                alt.Chart(pd.DataFrame({"x": [log_mic]}))
                .mark_rule(strokeDash=[3, 3], strokeWidth=1.0, color=ACCENT)
                .encode(x="x:Q")
            )

    # IC50 vertical rule
    if ic50_fit is not None:
        log_ic50 = np.log10(ic50_fit)
        if x_lo <= log_ic50 <= x_hi:
            layers.append(
                alt.Chart(pd.DataFrame({"x": [log_ic50]}))
                .mark_rule(strokeDash=[6, 4], strokeWidth=1.0, color=IC50_COLOR)
                .encode(x="x:Q")
            )

    # Text labels — top-right whitespace
    label_rows = []
    if mic_display is not None:
        label_rows.append({
            "x": x_hi - 0.15, "y": 1.20,
            "text": f"MIC₉₀ = {mic_display}",
            "color": ACCENT,
        })
    if ic50_fit is not None:
        label_rows.append({
            "x": x_hi - 0.15, "y": 1.08,
            "text": f"IC₅₀ = {_fmt_conc(ic50_fit)}",
            "color": IC50_COLOR,
        })
    if label_rows:
        lbl_df = pd.DataFrame(label_rows)
        layers.append(
            alt.Chart(lbl_df)
            .mark_text(
                align="right", fontSize=9, font=FONT_MONO,
            )
            .encode(
                x="x:Q", y="y:Q", text="text:N",
                color=alt.Color("color:N", scale=None),
            )
        )

    chart = (
        alt.layer(*layers)
        .properties(width=170, height=120)
    )
    return chart, ic50_fit, mic_display
