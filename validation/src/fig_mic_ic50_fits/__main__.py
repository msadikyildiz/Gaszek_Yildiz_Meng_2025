"""CLI: python -m fig_mic_ic50_fits [--suffix _12h]"""

import sys

from config import FIGURES_DIR
from . import make_figure

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

suffix = ""
if "--suffix" in sys.argv:
    suffix = sys.argv[sys.argv.index("--suffix") + 1]

make_figure(suffix=suffix)
