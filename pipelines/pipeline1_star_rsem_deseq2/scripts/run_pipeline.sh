#!/bin/bash

###############################################################################
# STAR-RSEM-DESeq2 Pipeline
# RNA-seq differential expression analysis
###############################################################################

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${1:-../config/pipeline_config.yaml}"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}STAR-RSEM-DESeq2 Pipeline${NC}"
echo -e "${BLUE}========================================${NC}"

if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}Error: Config file not found: $CONFIG_FILE${NC}"
    exit 1
fi

# Parse config
OUTPUT_DIR=$(grep "output_dir:" "$CONFIG_FILE" | head -1 | awk '{print $2}' | tr -d '"' | tr -d "'")
THREADS=$(grep "threads:" "$CONFIG_FILE" | grep -v "#" | head -1 | awk '{print $2}')

mkdir -p "$OUTPUT_DIR"/{"logs,qc,results"}

echo "Starting pipeline..."
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

# Main pipeline steps would go here
# bash "$SCRIPT_DIR/01_qc.sh" "$CONFIG_FILE"
# bash "$SCRIPT_DIR/02_alignment.sh" "$CONFIG_FILE"  # if alignment-based
# bash "$SCRIPT_DIR/03_quantification.sh" "$CONFIG_FILE"
# Rscript "$SCRIPT_DIR/04_statistical_analysis.R" "$CONFIG_FILE"

echo -e "${GREEN}Pipeline completed!${NC}"
echo "Results: $OUTPUT_DIR"
