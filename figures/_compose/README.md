# Vector figure composition (`figures/_compose/`)

Scripts that rebuild the main-text vector figures (Figure 3 and Figure 4) from
the processed data in `data/processed/`. Each figure is a two-step pipeline: an
`emit_*` script renders the individual panels as vector PDFs into `vector/`, and
a `compose_*` script assembles those panels into the final figure PDF.

| Script | Role |
|---|---|
| `emit_fig3_vector.py` | Emit the Figure 3 panels (A, B sequence logos, C, D) as vector PDFs into `vector/` |
| `compose_figure3_vector.py` | Assemble the emitted Figure 3 panels into the composed figure PDF |
| `emit_fig4_vector.py` | Emit the Figure 4 panels as vector PDFs into `vector/` |
| `compose_figure4_vector.py` | Assemble the emitted Figure 4 panels into the composed figure PDF |
| `figstyle.py` | Shared matplotlib typography and style (prefers Helvetica Neue, falls back to Helvetica/Arial/DejaVu Sans) |
| `compose_figure4.py` | Assemble the raster Figure 4 panels into the composed figure (uses `regen_gh.py` for panels G/H) |
| `regen_gh.py` | Regenerate Figure 4 panels G (R²) and H (ρ) as compact PNGs for the raster composite |

## Reproducing the vector figures

From the repository root, with the `figures` extra installed:

```bash
uv sync --extra figures
uv run python figures/_compose/emit_fig3_vector.py  &&  uv run python figures/_compose/compose_figure3_vector.py
uv run python figures/_compose/emit_fig4_vector.py  &&  uv run python figures/_compose/compose_figure4_vector.py
```

Run each `emit_*` step before its matching `compose_*` step: the compose script
reads the per-panel PDFs that the emit step writes into `vector/`.
