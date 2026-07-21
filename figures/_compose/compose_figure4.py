"""Compose manuscript Figure 4 from the eight rendered panels.

Layout (Nature double-column, 180 mm):
  row 1:  (A) AMP single mutants | (B) AZT single mutants
  row 2:  (C) AMP heat-map | (D) AZT heat-map | [ (G) R² / (H) ρ stacked ]
  row 3:  (E) AMP measured-vs-predicted density grid  (full width)
  row 4:  (F) AZT measured-vs-predicted density grid  (full width)
Panel letters are bold parenthesised Helvetica Neue, matching Figure 1.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

import figstyle

FIG = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "Figure_4_composite.png"

P = {
    # Panels A/B show per-replicate distributions (n = 3), consistent with
    # the Figure 4A/B caption.
    "A": "_compose/panel_A_dist.png",
    "B": "_compose/panel_B_dist.png",
    "C": "main/Figure 4C. Double mutants fitness AMP.png",
    "D": "main/Figure 4D. Double mutants fitness AZT.png",
    "E": "main/Figure 4E. AMP density grid.png",
    "F": "main/Figure 4F. AZT density grid.png",
    "G": "_compose/panel_G.png",
    "H": "_compose/panel_H.png",
}
img = {k: Image.open(FIG / v).convert("RGB") for k, v in P.items()}
asp = {k: im.width / im.height for k, im in img.items()}

# --- layout (mm) ---
W = 180.0
M_L, M_R, M_T, M_B = 6.0, 2.5, 3.0, 2.5
GX, GY = 3.5, 4.0          # gaps between panels / rows
GY_IN = 4.0                # gap between stacked G and H
cw = W - M_L - M_R          # content width

half = (cw - GX) / 2.0
# row 2: solve column widths so the stacked (G over H) height == heat-map height
_k = 1.0 / asp["G"] + 1.0 / asp["H"]
cd_w = (cw - 2 * GX + GY_IN / _k) / (2 + 1.0 / (asp["C"] * _k))
right_w = (cd_w / asp["C"] - GY_IN) / _k

h1 = half / asp["A"]
row2_h = cd_w / asp["C"]
h3 = cw / asp["E"]
h4 = cw / asp["F"]
total_h = M_T + M_B + h1 + row2_h + h3 + h4 + GY * 3

figstyle.apply()
fig = plt.figure(figsize=(W * figstyle.MM, total_h * figstyle.MM))
LET = dict(fontsize=10, fontweight="bold", family=figstyle.FONT, va="bottom",
           ha="left", color="#000000")


def place(letter, im, left, top, w, h):
    bottom = top - h
    ax = fig.add_axes([left / W, bottom / total_h, w / W, h / total_h])
    ax.imshow(im)
    ax.axis("off")
    fig.text(max((left - 6.0) / W, 0.002), top / total_h, f"({letter})", **LET)


y = total_h - M_T
# row 1
place("A", img["A"], M_L, y, half, h1)
place("B", img["B"], M_L + half + GX, y, half, h1)
y -= h1 + GY
# row 2
place("C", img["C"], M_L, y, cd_w, row2_h)
place("D", img["D"], M_L + cd_w + GX, y, cd_w, row2_h)
rx = M_L + 2 * cd_w + 2 * GX
g_h = right_w / asp["G"]
h_h = right_w / asp["H"]
place("G", img["G"], rx, y, right_w, g_h)
place("H", img["H"], rx, y - g_h - GY_IN, right_w, h_h)
y -= row2_h + GY
# row 3
place("E", img["E"], M_L, y, cw, h3)
y -= h3 + GY
# row 4
place("F", img["F"], M_L, y, cw, h4)

fig.savefig(OUT, dpi=600)
print(f"wrote {OUT}  ({W:.0f} x {total_h:.1f} mm @600dpi)")
