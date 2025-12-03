#!/bin/bash

###############################################################################
# 01_run_all_pipelines.sh
# Run all 8 RAPTOR RNA-seq pipelines on the same dataset
# Integrates with raptor.benchmark module
#
# Usage: bash 01_run_all_pipelines.sh <input_dir> <output_dir> <reference_dir> [pipelines]
#
# Author: Ayeh Bolouki
# Organization: RAPTOR Project
# License: MIT
###############################################################################

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0;33m' # No Color

# Functions
print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${CYAN}ℹ $1${NC}"
}

# Check arguments
if [ $# -lt 3 ]; then
    echo -e "${RED}Error: Missing required arguments${NC}"
    echo ""
    echo "Usage: $0 <input_dir> <output_dir> <reference_dir> [pipeline_ids]"
    echo ""
    echo "Arguments:"
    echo "  input_dir     : Directory containing FASTQ files and sample sheet"
    echo "  output_dir    : Base output directory for all pipeline results"
    echo "  reference_dir : Directory containing reference genome and annotation"
    echo "  pipeline_ids  : Optional comma-separated list (e.g., '1,3,4' or 'all')"
    echo ""
    echo "Examples:"
    echo "  # Run all 8 pipelines"
    echo "  $0 data/fastq results/benchmark references/human_grch38"
    echo ""
    echo "  # Run only pipelines 1, 3, and 4"
    echo "  $0 data/fastq results/benchmark references/human_grch38 1,3,4"
    echo ""
    echo "  # Run specific pipelines explicitly"
    echo "  $0 data/fastq results/benchmark references/human_grch38 all"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_BASE="$2"
REF_DIR="$3"
PIPELINE_IDS="${4:-all}"

# Validate directories
if [ ! -d "$INPUT_DIR" ]; then
    print_error "Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -d "$REF_DIR" ]; then
    print_error "Reference directory not found: $REF_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_BASE"

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RAPTOR_ROOT="$(dirname "$SCRIPT_DIR")"
PIPELINES_DIR="$RAPTOR_ROOT/pipelines"

# Check if pipelines directory exists
if [ ! -d "$PIPELINES_DIR" ]; then
    print_error "Pipelines directory not found: $PIPELINES_DIR"
    print_info "Expected location: $PIPELINES_DIR"
    exit 1
fi

# Log file
LOG_FILE="$OUTPUT_BASE/pipeline_comparison_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

print_header "RAPTOR Pipeline Comparison"
echo "Date: $(date)"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_BASE"
echo "Reference: $REF_DIR"
echo "Log: $LOG_FILE"
echo ""

# Find sample sheet
SAMPLE_SHEET=$(find "$INPUT_DIR" -name "*.csv" -o -name "sample*.txt" | head -1)
if [ -z "$SAMPLE_SHEET" ]; then
    print_warning "Sample sheet not found in $INPUT_DIR"
    print_info "Pipelines will need to auto-detect samples"
else
    print_info "Sample sheet: $SAMPLE_SHEET"
fi

# Count total threads available
TOTAL_THREADS=$(nproc 2>/dev/null || echo 8)
THREADS_PER_PIPELINE=$((TOTAL_THREADS / 2))
if [ $THREADS_PER_PIPELINE -lt 1 ]; then
    THREADS_PER_PIPELINE=1
fi
print_info "Available threads: $TOTAL_THREADS (using $THREADS_PER_PIPELINE per pipeline)"
echo ""

# Parse pipeline IDs
if [ "$PIPELINE_IDS" = "all" ]; then
    SELECTED_PIPELINES=(1 2 3 4 5 6 7 8)
else
    IFS=',' read -ra SELECTED_PIPELINES <<< "$PIPELINE_IDS"
fi

print_info "Pipelines to run: ${SELECTED_PIPELINES[*]}"
echo ""

# Pipeline information
declare -A PIPELINE_NAMES
PIPELINE_NAMES[1]="STAR-RSEM-DESeq2"
PIPELINE_NAMES[2]="HISAT2-StringTie-Ballgown"
PIPELINE_NAMES[3]="Salmon-edgeR"
PIPELINE_NAMES[4]="Kallisto-sleuth"
PIPELINE_NAMES[5]="STAR-HTSeq-limma"
PIPELINE_NAMES[6]="STAR-featureCounts-NOISeq"
PIPELINE_NAMES[7]="Bowtie2-RSEM-EBSeq"
PIPELINE_NAMES[8]="HISAT2-Cufflinks-Cuffdiff"

declare -A PIPELINE_DIRS
PIPELINE_DIRS[1]="pipeline1_star_rsem_deseq2"
PIPELINE_DIRS[2]="pipeline2_hisat2_stringtie_ballgown"
PIPELINE_DIRS[3]="pipeline3_salmon_edger"
PIPELINE_DIRS[4]="pipeline4_kallisto_sleuth"
PIPELINE_DIRS[5]="pipeline5_star_htseq_limma"
PIPELINE_DIRS[6]="pipeline6_star_featurecounts_noiseq"
PIPELINE_DIRS[7]="pipeline7_bowtie2_rsem_ebseq"
PIPELINE_DIRS[8]="pipeline8_hisat2_cufflinks_cuffdiff"

# Track results
declare -a SUCCESS_PIPELINES=()
declare -a FAILED_PIPELINES=()
START_TIME=$(date +%s)

# Run each selected pipeline
for PIPELINE_ID in "${SELECTED_PIPELINES[@]}"; do
    PIPELINE_NAME="${PIPELINE_NAMES[$PIPELINE_ID]}"
    PIPELINE_DIR_NAME="${PIPELINE_DIRS[$PIPELINE_ID]}"
    
    print_header "Pipeline $PIPELINE_ID: $PIPELINE_NAME"
    
    PIPELINE_DIR="$PIPELINES_DIR/$PIPELINE_DIR_NAME"
    PIPELINE_OUTPUT="$OUTPUT_BASE/$PIPELINE_DIR_NAME"
    
    # Check if pipeline exists
    if [ ! -d "$PIPELINE_DIR" ]; then
        print_warning "Pipeline directory not found: $PIPELINE_DIR"
        FAILED_PIPELINES+=("$PIPELINE_NAME (not found)")
        echo ""
        continue
    fi
    
    # Check if config exists
    CONFIG_FILE="$PIPELINE_DIR/config/pipeline_config.yaml"
    if [ ! -f "$CONFIG_FILE" ]; then
        print_warning "Config file not found: $CONFIG_FILE"
        FAILED_PIPELINES+=("$PIPELINE_NAME (no config)")
        echo ""
        continue
    fi
    
    # Create pipeline output directory
    mkdir -p "$PIPELINE_OUTPUT"
    
    # Create runtime config with updated paths
    TEMP_CONFIG="$PIPELINE_OUTPUT/runtime_config.yaml"
    cp "$CONFIG_FILE" "$TEMP_CONFIG"
    
    # Update paths in config (using sed)
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        sed -i '' "s|input_dir:.*|input_dir: \"$INPUT_DIR\"|g" "$TEMP_CONFIG"
        sed -i '' "s|output_dir:.*|output_dir: \"$PIPELINE_OUTPUT\"|g" "$TEMP_CONFIG"
        sed -i '' "s|threads:.*|threads: $THREADS_PER_PIPELINE|g" "$TEMP_CONFIG"
        sed -i '' "s|reference_genome:.*|reference_genome: \"$REF_DIR/genome.fa\"|g" "$TEMP_CONFIG"
        sed -i '' "s|annotation_gtf:.*|annotation_gtf: \"$REF_DIR/annotation.gtf\"|g" "$TEMP_CONFIG"
        sed -i '' "s|transcriptome_fasta:.*|transcriptome_fasta: \"$REF_DIR/transcriptome.fa\"|g" "$TEMP_CONFIG"
        
        if [ -n "$SAMPLE_SHEET" ]; then
            sed -i '' "s|sample_sheet:.*|sample_sheet: \"$SAMPLE_SHEET\"|g" "$TEMP_CONFIG"
        fi
    else
        # Linux
        sed -i "s|input_dir:.*|input_dir: \"$INPUT_DIR\"|g" "$TEMP_CONFIG"
        sed -i "s|output_dir:.*|output_dir: \"$PIPELINE_OUTPUT\"|g" "$TEMP_CONFIG"
        sed -i "s|threads:.*|threads: $THREADS_PER_PIPELINE|g" "$TEMP_CONFIG"
        sed -i "s|reference_genome:.*|reference_genome: \"$REF_DIR/genome.fa\"|g" "$TEMP_CONFIG"
        sed -i "s|annotation_gtf:.*|annotation_gtf: \"$REF_DIR/annotation.gtf\"|g" "$TEMP_CONFIG"
        sed -i "s|transcriptome_fasta:.*|transcriptome_fasta: \"$REF_DIR/transcriptome.fa\"|g" "$TEMP_CONFIG"
        
        if [ -n "$SAMPLE_SHEET" ]; then
            sed -i "s|sample_sheet:.*|sample_sheet: \"$SAMPLE_SHEET\"|g" "$TEMP_CONFIG"
        fi
    fi
    
    print_info "Running $PIPELINE_NAME..."
    PIPELINE_START=$(date +%s)
    
    # Check for run script
    RUN_SCRIPT="$PIPELINE_DIR/scripts/run_pipeline.sh"
    if [ ! -f "$RUN_SCRIPT" ]; then
        print_warning "Run script not found: $RUN_SCRIPT"
        FAILED_PIPELINES+=("$PIPELINE_NAME (no run script)")
        echo ""
        continue
    fi
    
    # Make script executable
    chmod +x "$RUN_SCRIPT"
    
    # Run pipeline
    if bash "$RUN_SCRIPT" "$TEMP_CONFIG" 2>&1; then
        PIPELINE_END=$(date +%s)
        PIPELINE_TIME=$((PIPELINE_END - PIPELINE_START))
        print_success "$PIPELINE_NAME completed in ${PIPELINE_TIME}s"
        SUCCESS_PIPELINES+=("$PIPELINE_NAME")
        
        # Save timing info
        echo "$PIPELINE_TIME" > "$PIPELINE_OUTPUT/.runtime"
    else
        print_error "$PIPELINE_NAME failed"
        FAILED_PIPELINES+=("$PIPELINE_NAME (runtime error)")
    fi
    
    echo ""
done

# Summary
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

print_header "Pipeline Comparison Summary"
echo "Total runtime: ${TOTAL_TIME}s ($((TOTAL_TIME / 60)) minutes)"
echo ""

echo -e "${GREEN}Successful pipelines (${#SUCCESS_PIPELINES[@]}/${#SELECTED_PIPELINES[@]}):${NC}"
if [ ${#SUCCESS_PIPELINES[@]} -eq 0 ]; then
    echo "  None"
else
    for pipeline in "${SUCCESS_PIPELINES[@]}"; do
        echo "  ✓ $pipeline"
    done
fi
echo ""

echo -e "${RED}Failed pipelines (${#FAILED_PIPELINES[@]}/${#SELECTED_PIPELINES[@]}):${NC}"
if [ ${#FAILED_PIPELINES[@]} -eq 0 ]; then
    echo "  None"
else
    for pipeline in "${FAILED_PIPELINES[@]}"; do
        echo "  ✗ $pipeline"
    done
fi
echo ""

# Create results summary
SUMMARY_FILE="$OUTPUT_BASE/comparison_summary.txt"
cat > "$SUMMARY_FILE" << EOF
RAPTOR Pipeline Comparison Summary
===================================

Date: $(date)
Input: $INPUT_DIR
Output: $OUTPUT_BASE
Reference: $REF_DIR

Pipelines Run: ${#SELECTED_PIPELINES[@]}
Successful: ${#SUCCESS_PIPELINES[@]}
Failed: ${#FAILED_PIPELINES[@]}
Total runtime: ${TOTAL_TIME}s

Successful pipelines:
EOF

for pipeline in "${SUCCESS_PIPELINES[@]}"; do
    echo "  - $pipeline" >> "$SUMMARY_FILE"
done

echo "" >> "$SUMMARY_FILE"
echo "Failed pipelines:" >> "$SUMMARY_FILE"

for pipeline in "${FAILED_PIPELINES[@]}"; do
    echo "  - $pipeline" >> "$SUMMARY_FILE"
done

print_success "Summary saved to: $SUMMARY_FILE"
print_info "Full log saved to: $LOG_FILE"

echo ""
print_info "Next steps:"
echo "  1. Review individual pipeline outputs in: $OUTPUT_BASE"
echo "  2. Run: Rscript scripts/03_compare_results.R $OUTPUT_BASE"
echo "  3. Run: Rscript scripts/04_visualize_comparison.R $OUTPUT_BASE"
echo ""

if [ ${#FAILED_PIPELINES[@]} -eq 0 ]; then
    print_success "All pipelines completed successfully! 🦖🎉"
    exit 0
else
    print_warning "Some pipelines failed. Check logs for details."
    exit 1
fi
