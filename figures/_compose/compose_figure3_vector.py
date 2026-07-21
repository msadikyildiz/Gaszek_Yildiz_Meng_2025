"""Compose manuscript Figure 3 as a single vector PDF from the panels emitted
by ``emit_fig3_vector.py``, placed with pymupdf ``show_pdf_page`` so text and
line-art stay vector (dense density/scatter layers are rasterised in place).

Vertical stack matching the manuscript figure:
  (A) AMP-vs-AZT fitness density, 2x3 concentration grid   (full width)
  (B) top-1% sequence logos: Ampicillin | Aztreonam        (two half-width)
  (C) fitness in AMP vs number of mutated residues         (full width)
  (D) fitness in AZT vs number of mutated residues         (full width)
"""
from __future__ import annotations

from pathlib import Path

import fitz  # pymupdf

HERE = Path(__file__).resolve().parent
VEC = HERE / "vector"
OUT = HERE / "Figure_3.pdf"
PT = 72.0 / 25.4

src = {k: fitz.open(VEC / f"fig3_{k}.pdf") for k in ("A", "Bamp", "Bazt", "C", "D")}
asp = {k: d[0].rect.width / d[0].rect.height for k, d in src.items()}

# --- layout (mm) ---
W = 180.0
M_L, M_R, M_T, M_B = 5.0, 3.0, 3.0, 3.0
GY = 5.0          # gap between panel rows
GX_B = 4.0        # gap between the two logos
cw = W - M_L - M_R

hA = cw / asp["A"]
bw = (cw - GX_B) / 2.0
hB = max(bw / asp["Bamp"], bw / asp["Bazt"])
hC = cw / asp["C"]
hD = cw / asp["D"]
total_h = M_T + M_B + hA + hB + hC + hD + GY * 3

# Fit within the Nature type area (247 mm) by uniform scale (no distortion).
MAX_H = 247.0
SCALE = min(1.0, MAX_H / total_h)
PTS = PT * SCALE

doc = fitz.open()
page = doc.new_page(width=W * PTS, height=total_h * PTS)


def place(key: str, left: float, top_mm: float, w: float, h: float) -> None:
    """top_mm = mm from page TOP to the panel top edge."""
    rect = fitz.Rect(left * PTS, top_mm * PTS, (left + w) * PTS, (top_mm + h) * PTS)
    page.show_pdf_page(rect, src[key], 0)


def letter(txt: str, left_mm: float, top_mm: float) -> None:
    page.insert_text(fitz.Point(max(left_mm, 1.0) * PTS, top_mm * PTS + 9.0 * SCALE),
                     txt, fontname="hebo", fontsize=10 * SCALE, color=(0, 0, 0))


y = M_T
place("A", M_L, y, cw, hA)
letter("(A)", 1.0, y)
y += hA + GY
place("Bamp", M_L, y, bw, hB)
place("Bazt", M_L + bw + GX_B, y, bw, hB)
letter("(B)", 1.0, y)
y += hB + GY
place("C", M_L, y, cw, hC)
letter("(C)", 1.0, y)
y += hC + GY
place("D", M_L, y, cw, hD)
letter("(D)", 1.0, y)

doc.save(OUT, deflate=True, garbage=4)
print(f"wrote {OUT.name}  ({W * SCALE:.1f} x {total_h * SCALE:.1f} mm, scale={SCALE:.3f})")
