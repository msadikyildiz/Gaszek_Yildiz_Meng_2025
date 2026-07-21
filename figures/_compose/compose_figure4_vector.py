"""Compose manuscript Figure 4 as a single vector PDF from the eight vector
panels emitted by ``emit_fig4_vector.py``.

Reuses the exact 180 mm layout math of ``compose_figure4.py`` (raster composer)
but places each panel with pymupdf ``show_pdf_page`` so text and marks stay
vector. Panel letters are drawn as vector Helvetica-Bold text.

Layout (Nature double-column, 180 mm):
  row 1:  (A) AMP single mutants | (B) AZT single mutants
  row 2:  (C) AMP heat-map | (D) AZT heat-map | [ (G) R2 / (H) rho stacked ]
  row 3:  (E) AMP measured-vs-predicted density grid  (full width)
  row 4:  (F) AZT measured-vs-predicted density grid  (full width)
"""
from __future__ import annotations

from pathlib import Path

import fitz  # pymupdf

HERE = Path(__file__).resolve().parent
VEC = HERE / "vector"
OUT = HERE / "Figure_4.pdf"

PT = 72.0 / 25.4  # mm -> points

src = {k: fitz.open(VEC / f"panel_{k}.pdf") for k in "ABCDEFGH"}
asp = {k: d[0].rect.width / d[0].rect.height for k, d in src.items()}

# --- layout (mm), identical to compose_figure4.py ---
W = 180.0
M_L, M_R, M_T, M_B = 6.0, 2.5, 3.0, 2.5
GX, GY = 3.5, 4.0
GY_IN = 4.0
cw = W - M_L - M_R

half = (cw - GX) / 2.0
_k = 1.0 / asp["G"] + 1.0 / asp["H"]
cd_w = (cw - 2 * GX + GY_IN / _k) / (2 + 1.0 / (asp["C"] * _k))
right_w = (cd_w / asp["C"] - GY_IN) / _k

h1 = half / asp["A"]
row2_h = cd_w / asp["C"]
h3 = cw / asp["E"]
h4 = cw / asp["F"]
total_h = M_T + M_B + h1 + row2_h + h3 + h4 + GY * 3

doc = fitz.open()
page = doc.new_page(width=W * PT, height=total_h * PT)


def place(letter: str, left: float, top: float, w: float, h: float) -> None:
    """`top` is mm from the figure BOTTOM to the panel top (matplotlib origin)."""
    y_top = total_h - top          # mm from page top to panel top
    rect = fitz.Rect(left * PT, y_top * PT, (left + w) * PT, (y_top + h) * PT)
    page.show_pdf_page(rect, src[letter], 0)
    page.insert_text(
        fitz.Point(max(left - 5.5, 1.0) * PT, (y_top * PT) + 9.5),
        f"({letter})", fontname="hebo", fontsize=10, color=(0, 0, 0))


y = total_h - M_T
# row 1
place("A", M_L, y, half, h1)
place("B", M_L + half + GX, y, half, h1)
y -= h1 + GY
# row 2
place("C", M_L, y, cd_w, row2_h)
place("D", M_L + cd_w + GX, y, cd_w, row2_h)
rx = M_L + 2 * cd_w + 2 * GX
g_h = right_w / asp["G"]
h_h = right_w / asp["H"]
place("G", rx, y, right_w, g_h)
place("H", rx, y - g_h - GY_IN, right_w, h_h)
y -= row2_h + GY
# row 3
place("E", M_L, y, cw, h3)
y -= h3 + GY
# row 4
place("F", M_L, y, cw, h4)

doc.save(OUT, deflate=True, garbage=4)
print(f"wrote {OUT.name}  ({W:.0f} x {total_h:.1f} mm)")
