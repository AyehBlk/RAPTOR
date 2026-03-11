#!/bin/bash

###############################################################################
# RAPTOR Quick Kallisto Pipeline
# Ultra-fast pseudo-alignment for QC and profiling
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
PURPLE='\033[0;35m'
NC='\033[0m' # No Color

# =============================================================================
# Functions
# =============================================================================

print_banner() {
    echo -e "${PURPLE}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                🦖 RAPTOR Quick Kallisto                      ║"
    echo "║       Ultra-Fast Pseudo-alignment for QC & Profiling         ║"
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

check_kallisto() {
    if ! command -v kallisto &> /dev/null; then
        log_error "Kallisto not found in PATH"
        log_info "Install with: conda install -c bioconda kallisto"
        exit 1
    fi
    
    local version=$(kallisto version 2>&1 | grep -oP '\d+\.\d+\.\d+' | head -1)
    log_info "Kallisto version: $version"
}

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "RAPTOR Quick Kallisto Pipeline - Module 1: Quantify"
    echo ""
    echo "Options:"
    echo "  -s, --sample-sheet FILE Sample sheet CSV (required if no --input)"
    echo "  -i, --input DIR         FASTQ directory (required if no --sample-sheet)"
    echo "  -x, --index FILE        Kallisto index .idx file (required)"
    echo "  -t, --tx2gene FILE      Transcript to gene mapping (required)"
    echo "  -o, --output DIR        Output directory (default: $DEFAULT_OUTPUT_DIR)"
    echo "  -p, --threads INT       Number of threads (default: 8)"
    echo "  -r, --read-type TYPE    Read type: auto, paired, single (default: auto)"
    echo "  -l, --frag-length INT   Fragment length for single-end (default: 200)"
    echo "  -d, --frag-sd INT       Fragment length SD for single-end (default: 20)"
    echo "  --keep-quant            Keep individual quant files"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "⚠️  IMPORTANT: For single-end reads, you MUST provide --frag-length and --frag-sd"
    echo "   Kallisto cannot estimate fragment length from single-end data!"
    echo ""
    echo "Output Files:"
    echo "  $OUTPUT_COUNTS_FILE   - Gene-level count matrix (for M2-M4)"
    echo "  $OUTPUT_TPM_FILE       - TPM normalized matrix"
    echo "  sample_info.csv        - Sample metadata and QC metrics"
    echo ""
    echo "Examples:"
    echo "  # Paired-end with sample sheet"
    echo "  $0 -s samples.csv -x transcripts.idx -t tx2gene.tsv"
    echo ""
    echo "  # Auto-detect from directory"
    echo "  $0 -i fastq/ -x transcripts.idx -t tx2gene.tsv"
    echo ""
    echo "  # Single-end (MUST specify fragment length!)"
    echo "  $0 -i fastq/ -x transcripts.idx -t tx2gene.tsv -r single -l 200 -d 20"
    echo ""
}

# =============================================================================
# Parse Arguments
# =============================================================================

SAMPLE_SHEET=""
INPUT_DIR=""
INDEX_FILE=""
TX2GENE=""
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
THREADS=8
READ_TYPE="auto"
FRAG_LENGTH=200
FRAG_SD=20
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
            INDEX_FILE="$2"
            shift 2
            ;;
        -t|--tx2gene)
            TX2GENE="$2"
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
        -r|--read-type)
            READ_TYPE="$2"
            shift 2
            ;;
        -l|--frag-length)
            FRAG_LENGTH="$2"
            shift 2
            ;;
        -d|--frag-sd)
            FRAG_SD="$2"
            shift 2
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

# Check Kallisto
check_kallisto

# Validate inputs
if [[ -z "$SAMPLE_SHEET" && -z "$INPUT_DIR" ]]; then
    log_error "Either --sample-sheet or --input is required"
    usage
    exit 1
fi

if [[ -z "$INDEX_FILE" ]]; then
    log_error "--index is required"
    exit 1
fi

if [[ ! -f "$INDEX_FILE" ]]; then
    log_error "Kallisto index not found: $INDEX_FILE"
    exit 1
fi

if [[ -z "$TX2GENE" ]]; then
    log_error "--tx2gene is required for gene-level counts"
    exit 1
fi

if [[ ! -f "$TX2GENE" ]]; then
    log_error "tx2gene file not found: $TX2GENE"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"/{kallisto_quant,logs}

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
        --output "$SAMPLE_INFO" \
        --read-type "$READ_TYPE"
fi

# Count samples
N_SAMPLES=$(tail -n +2 "$SAMPLE_INFO" | wc -l)
log_info "Found $N_SAMPLES samples"

# Determine read type if auto
if [[ "$READ_TYPE" == "auto" ]]; then
    # Check if fastq_r2 column has values
    HAS_R2=$(awk -F',' 'NR>1 && $5!="" {print "yes"; exit}' "$SAMPLE_INFO")
    if [[ "$HAS_R2" == "yes" ]]; then
        READ_TYPE="paired"
    else
        READ_TYPE="single"
    fi
    log_info "Detected read type: $READ_TYPE"
fi

# =============================================================================
# Step 2: Run Kallisto Quantification
# =============================================================================

log_info "Step 2: Running Kallisto quantification..."

# Process each sample
SAMPLE_NUM=0
while IFS=, read -r sample_id condition batch fastq_r1 fastq_r2; do
    # Skip header
    if [[ "$sample_id" == "sample_id" ]]; then
        continue
    fi
    
    SAMPLE_NUM=$((SAMPLE_NUM + 1))
    QUANT_DIR="$OUTPUT_DIR/kallisto_quant/$sample_id"
    
    log_info "[$SAMPLE_NUM/$N_SAMPLES] Processing: $sample_id"
    
    # Prepare Kallisto command
    KALLISTO_CMD="kallisto quant -i $INDEX_FILE -o $QUANT_DIR -t $THREADS --plaintext"
    
    # Check if paired or single
    if [[ "$READ_TYPE" == "paired" && -n "$fastq_r2" && "$fastq_r2" != "" && "$fastq_r2" != "NA" ]]; then
        # Paired-end
        log_info "  Mode: Paired-end"
        KALLISTO_CMD="$KALLISTO_CMD $fastq_r1 $fastq_r2"
    else
        # Single-end - MUST specify fragment length
        log_info "  Mode: Single-end (fragment length: $FRAG_LENGTH ± $FRAG_SD)"
        KALLISTO_CMD="$KALLISTO_CMD --single -l $FRAG_LENGTH -s $FRAG_SD $fastq_r1"
    fi
    
    # Run Kallisto
    log_info "  Running Kallisto..."
    eval "$KALLISTO_CMD" >> "$LOG_FILE" 2>&1
    
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

python3 "$SCRIPT_DIR/combine_counts.py" \
    --quant-dir "$OUTPUT_DIR/kallisto_quant" \
    --tx2gene "$TX2GENE" \
    --output "$OUTPUT_DIR" \
    --sample-info "$SAMPLE_INFO" \
    --tool kallisto

# =============================================================================
# Step 4: Generate Summary Report
# =============================================================================

log_info "Step 4: Generating summary report..."

python3 "$SCRIPT_DIR/generate_report.py" \
    --quant-dir "$OUTPUT_DIR/kallisto_quant" \
    --counts "$OUTPUT_DIR/$OUTPUT_COUNTS_FILE" \
    --output "$OUTPUT_DIR/qc_report.html" \
    --tool kallisto

# =============================================================================
# Step 5: Cleanup (optional)
# =============================================================================

if [[ "$KEEP_QUANT" == false ]]; then
    log_info "Cleaning up intermediate files..."
    rm -rf "$OUTPUT_DIR/kallisto_quant"
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
