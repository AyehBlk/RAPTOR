#!/bin/bash

###############################################################################
# RAPTOR Quick Salmon Pipeline
# Fast pseudo-alignment for QC and profiling
#
# Author: Ayeh Bolouki
# Email: ayehbolouki1988@gmail.com
# Version: 2.2.0
###############################################################################

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# =============================================================================
# Constants - Architecture Compliant
# =============================================================================
DEFAULT_OUTPUT_DIR="results/quick_counts"
OUTPUT_COUNTS_FILE="quick_gene_counts.csv"
OUTPUT_TPM_FILE="quick_tpm.csv"

# =============================================================================
# Color Output
# =============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# Functions
# =============================================================================

print_banner() {
    echo -e "${BLUE}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                🦖 RAPTOR Quick Salmon                        ║"
    echo "║         Fast Pseudo-alignment for QC & Profiling             ║"
    echo "║                      Version 2.2.0                           ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

check_salmon() {
    if ! command -v salmon &> /dev/null; then
        log_error "Salmon not found in PATH"
        log_info "Install with: conda install -c bioconda salmon"
        exit 1
    fi
    
    local version=$(salmon --version 2>&1 | grep -oP '\d+\.\d+\.\d+' | head -1)
    log_info "Salmon version: $version"
}

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "RAPTOR Quick Salmon Pipeline - Module 1: Quantify"
    echo ""
    echo "Options:"
    echo "  -s, --sample-sheet FILE Sample sheet CSV (required if no --input)"
    echo "  -i, --input DIR         FASTQ directory (required if no --sample-sheet)"
    echo "  -x, --index DIR         Salmon index directory (required)"
    echo "  -g, --gene-map FILE     Transcript to gene mapping (optional)"
    echo "  -o, --output DIR        Output directory (default: $DEFAULT_OUTPUT_DIR)"
    echo "  -p, --threads INT       Number of threads (default: 8)"
    echo "  -l, --lib-type TYPE     Library type (default: A=auto)"
    echo "  --gc-bias               Enable GC bias correction (default: on)"
    echo "  --no-gc-bias            Disable GC bias correction"
    echo "  --validate-mappings     Enable mapping validation (default: on)"
    echo "  --keep-quant            Keep individual quant files"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Output Files:"
    echo "  $OUTPUT_COUNTS_FILE   - Gene-level count matrix (for M2-M4)"
    echo "  $OUTPUT_TPM_FILE       - TPM normalized matrix"
    echo "  sample_info.csv        - Sample metadata and QC metrics"
    echo ""
    echo "Examples:"
    echo "  # With sample sheet"
    echo "  $0 -s samples.csv -x salmon_index/ -o results/quick_counts"
    echo ""
    echo "  # Auto-detect from directory"
    echo "  $0 -i fastq/ -x salmon_index/ -g tx2gene.csv"
    echo ""
    echo "  # With gene-level aggregation"
    echo "  $0 -s samples.csv -x salmon_index/ -g tx2gene.csv"
    echo ""
}

# =============================================================================
# Parse Arguments
# =============================================================================

SAMPLE_SHEET=""
INPUT_DIR=""
INDEX_DIR=""
GENE_MAP=""
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
THREADS=8
LIB_TYPE="A"
GC_BIAS="--gcBias"
VALIDATE="--validateMappings"
KEEP_QUANT=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample-sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -x|--index)
            INDEX_DIR="$2"
            shift 2
            ;;
        -g|--gene-map)
            GENE_MAP="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--threads)
            THREADS="$2"
            shift 2
            ;;
        -l|--lib-type)
            LIB_TYPE="$2"
            shift 2
            ;;
        --gc-bias)
            GC_BIAS="--gcBias"
            shift
            ;;
        --no-gc-bias)
            GC_BIAS=""
            shift
            ;;
        --validate-mappings)
            VALIDATE="--validateMappings"
            shift
            ;;
        --keep-quant)
            KEEP_QUANT=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# =============================================================================
# Main Pipeline
# =============================================================================

print_banner

# Check Salmon
check_salmon

# Validate inputs
if [[ -z "$SAMPLE_SHEET" && -z "$INPUT_DIR" ]]; then
    log_error "Either --sample-sheet or --input is required"
    usage
    exit 1
fi

if [[ -z "$INDEX_DIR" ]]; then
    log_error "--index is required"
    exit 1
fi

if [[ ! -d "$INDEX_DIR" ]]; then
    log_error "Salmon index not found: $INDEX_DIR"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"/{salmon_quant,logs}

LOG_FILE="$OUTPUT_DIR/logs/pipeline_$(date '+%Y%m%d_%H%M%S').log"

log_info "Output directory: $OUTPUT_DIR"
log_info "Log file: $LOG_FILE"

# =============================================================================
# Step 1: Parse Sample Sheet or Auto-Detect
# =============================================================================

log_info "Step 1: Processing sample information..."

SAMPLE_INFO="$OUTPUT_DIR/sample_info.csv"

if [[ -n "$SAMPLE_SHEET" ]]; then
    log_info "Using sample sheet: $SAMPLE_SHEET"
    cp "$SAMPLE_SHEET" "$SAMPLE_INFO"
else
    log_info "Auto-detecting from directory: $INPUT_DIR"
    python3 "$SCRIPT_DIR/detect_samples.py" \
        --input "$INPUT_DIR" \
        --output "$SAMPLE_INFO"
fi

# Count samples
N_SAMPLES=$(tail -n +2 "$SAMPLE_INFO" | wc -l)
log_info "Found $N_SAMPLES samples"

# =============================================================================
# Step 2: Run Salmon Quantification
# =============================================================================

log_info "Step 2: Running Salmon quantification..."

# Process each sample
SAMPLE_NUM=0
while IFS=, read -r sample_id condition batch fastq_r1 fastq_r2; do
    # Skip header
    if [[ "$sample_id" == "sample_id" ]]; then
        continue
    fi
    
    SAMPLE_NUM=$((SAMPLE_NUM + 1))
    QUANT_DIR="$OUTPUT_DIR/salmon_quant/$sample_id"
    
    log_info "[$SAMPLE_NUM/$N_SAMPLES] Processing: $sample_id"
    
    # Prepare Salmon command
    SALMON_CMD="salmon quant -i $INDEX_DIR -l $LIB_TYPE -o $QUANT_DIR -p $THREADS $GC_BIAS $VALIDATE"
    
    # Add gene map if provided
    if [[ -n "$GENE_MAP" && -f "$GENE_MAP" ]]; then
        SALMON_CMD="$SALMON_CMD -g $GENE_MAP"
    fi
    
    # Check if paired or single
    if [[ -n "$fastq_r2" && "$fastq_r2" != "" && "$fastq_r2" != "NA" ]]; then
        # Paired-end
        log_info "  Mode: Paired-end"
        SALMON_CMD="$SALMON_CMD -1 $fastq_r1 -2 $fastq_r2"
    else
        # Single-end
        log_info "  Mode: Single-end"
        SALMON_CMD="$SALMON_CMD -r $fastq_r1"
    fi
    
    # Run Salmon
    log_info "  Running Salmon..."
    eval "$SALMON_CMD" >> "$LOG_FILE" 2>&1
    
    if [[ $? -eq 0 ]]; then
        log_info "  ✓ Complete"
    else
        log_error "  ✗ Failed - check log: $LOG_FILE"
        exit 1
    fi
    
done < "$SAMPLE_INFO"

# =============================================================================
# Step 3: Combine Counts
# =============================================================================

log_info "Step 3: Combining count matrices..."

if [[ -n "$GENE_MAP" && -f "$GENE_MAP" ]]; then
    python3 "$SCRIPT_DIR/combine_counts.py" \
        --quant-dir "$OUTPUT_DIR/salmon_quant" \
        --tx2gene "$GENE_MAP" \
        --output "$OUTPUT_DIR" \
        --sample-info "$SAMPLE_INFO" \
        --tool salmon
else
    log_warn "No gene map provided - skipping gene-level aggregation"
    log_info "Output will be at transcript level"
fi

# =============================================================================
# Step 4: Generate Summary Report
# =============================================================================

log_info "Step 4: Generating summary report..."

python3 "$SCRIPT_DIR/generate_report.py" \
    --quant-dir "$OUTPUT_DIR/salmon_quant" \
    --counts "$OUTPUT_DIR/$OUTPUT_COUNTS_FILE" \
    --output "$OUTPUT_DIR/qc_report.html" \
    --tool salmon

# =============================================================================
# Step 5: Cleanup (optional)
# =============================================================================

if [[ "$KEEP_QUANT" == false ]]; then
    log_info "Cleaning up intermediate files..."
    rm -rf "$OUTPUT_DIR/salmon_quant"
fi

# =============================================================================
# Done
# =============================================================================

echo ""
echo -e "${GREEN}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}                    Pipeline Complete!                          ${NC}"
echo -e "${GREEN}═══════════════════════════════════════════════════════════════${NC}"
echo ""
echo "Output files:"
echo "  • Count matrix (gene):  $OUTPUT_DIR/$OUTPUT_COUNTS_FILE"
echo "  • Count matrix (TPM):   $OUTPUT_DIR/$OUTPUT_TPM_FILE"
echo "  • Sample info:          $OUTPUT_DIR/sample_info.csv"
echo "  • QC report:            $OUTPUT_DIR/qc_report.html"
echo ""
echo "Next steps (RAPTOR Workflow):"
echo "  1. Module 2 - Quality Assessment:"
echo "     raptor qc --counts $OUTPUT_DIR/$OUTPUT_COUNTS_FILE"
echo ""
echo "  2. Module 3 - Data Profiling:"
echo "     raptor profile --counts $OUTPUT_DIR/$OUTPUT_COUNTS_FILE"
echo ""
echo "  3. Module 4 - Get Recommendation:"
echo "     raptor recommend"
echo ""
