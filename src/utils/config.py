# Configuration module for TEM1 epistasis analysis pipeline
from pathlib import Path

# Project directory setup
#
# Historically the project assumed an absolute path on the Toprak lab
# server.  When running the pipeline outside of that environment the
# absolute path does not exist and all scripts fail to locate data and
# cache directories.  To make the code portable we derive the project
# root relative to this configuration file.

# `config.py` lives in `src/utils/`, therefore two ``parents`` levels
# up from this file points to the repository root.
PROJECT_PATH = Path(__file__).resolve().parents[2]
RUN_ID = "full"
CACHE_DIR = PROJECT_PATH / "cache" / RUN_ID

# Mutation settings
# INTENDED: Maps positions to amino acid variants (original aa first)
INTENDED = {
    19: ['L','P'], 37: ['Q','K'], 67: ['M','L','V'],
    102: ['E','K'], 162: ['R','S','H','N'],
    180: ['M','T'], 235: ['A','T'], 236: ['G','S'],
    237: ['E','K'], 241: ['R','S','C'], 261: ['T','M'],
    271: ['R','L','Q'], 272:['N','D']
}

# ALL_MUTATIONS: Maps positions to variants with '.' representing wildtype
ALL_MUTATIONS = {
    19: ['.','P'], 37: ['.','K'], 67: ['.','L','V'],
    102: ['.','K'], 162: ['.','S','H','N'],
    180: ['.','T'], 235: ['.','T'], 236: ['.','S'],
    237: ['.','K'], 241: ['.','S','C'], 261: ['.','M'],
    271: ['.','L','Q'], 272:['.','D']
}

def get_mutation_options(intended=INTENDED):
    """Returns list with number of possible states at each position"""
    return [len(aas) for pos, aas in intended.items()]

# Ensure cache directory exists
CACHE_DIR.mkdir(parents=True, exist_ok=True)
