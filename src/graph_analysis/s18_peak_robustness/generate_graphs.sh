#!/bin/bash

# run this script to generate the graphs used for fitness advantage analysis and global peak presence analysis.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." &> /dev/null && pwd )"
OUTPUT_DIR="$SCRIPT_DIR/outputs/global-peak-robustness"

python -m fitness_landscape_graph.build_graphs_parallel \
    --base-path "$PROJECT_ROOT" \
    --output-dir "$OUTPUT_DIR" \
    --neutral-thresholds 0.15 0.45 0.01 \
    --concentrations 0.0 0.44 1.33 4.0 12.0 36.0 108.0 324.0
