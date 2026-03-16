#!/bin/bash
# =============================================================================
# RAPTOR v2.2.0 - Kallisto Quantification Script
# =============================================================================
#
# Enhanced bash script with v2.2.0 features:
# - Input validation with clear error messages
# - Better error handling
# - Progress indicators
# - Quality checks
#
# Usage: ./run_kallisto.sh <index> <r1> <r2|-> <output_dir> [threads] [fragment_length] [fragment_sd]
#
# Arguments:
#   index           : Path to Kallisto index (.idx file)
#   r1              : Path to R1 FASTQ file
#   r2              : Path to R2 FASTQ file (or "-" for single-end)
#   output_dir      : Output directory for quantification
#   threads         : Number of threads (optional, default: 8)
#   fragment_length : Mean fragment length for single-end (optional)
#   fragment_sd     : Fragment length SD for single-end (optional)
#
# Examples:
#   # Paired-end
#   ./run_kallisto.sh kallisto.idx sample_R1.fq.gz sample_R2.fq.gz output/ 16
#
#   # Single-end (REQUIRES fragment length)
#   ./run_kallisto.sh kallisto.idx sample.fq.gz - output/ 8 200 20
#
# Author: Ayeh Bolouki
# Version: 2.2.0
# =============================================================================

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# =============================================================================
# COLOR CODES FOR OUTPUT
# =============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# LOGGING FUNCTIONS (v2.2.0)
# =============================================================================
log_info() {
    echo -e "${GREEN}✓${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
}

log_error() {
    echo -e "${RED}✗${NC} $1" >&2
}

log_step() {
    echo -e "${BLUE}▶${NC} $1"
}

# =============================================================================
# INPUT VALIDATION (v2.2.0)
# =============================================================================
validate_inputs() {
    local errors=0
    
    # Check argument count
    if [ $# -lt 4 ]; then
        log_error "Insufficient arguments"
        echo "Usage: $0 <index> <r1> <r2|-> <output_dir> [threads] [fragment_length] [fragment_sd]"
        echo ""
        echo "Examples:"
        echo "  Paired-end: $0 kallisto.idx R1.fq.gz R2.fq.gz output/ 16"
        echo "  Single-end: $0 kallisto.idx sample.fq.gz - output/ 8 200 20"
        exit 1
    fi
    
    # Check Kallisto index exists
    if [ ! -f "$INDEX" ]; then
        log_error "Kallisto index not found: $INDEX"
        log_error "  Build index with: kallisto index -i index.idx transcripts.fa"
        ((errors++))
    else
        log_info "Kallisto index validated: $INDEX"
    fi
    
    # Check R1 exists
    if [ ! -f "$R1" ]; then
        log_error "R1 file not found: $R1"
        ((errors++))
    else
        log_info "R1 file found: $R1"
    fi
    
    # Check R2 if paired-end
    if [ -n "$R2" ] && [ "$R2" != "-" ]; then
        if [ ! -f "$R2" ]; then
            log_error "R2 file not found: $R2"
            ((errors++))
        else
            log_info "R2 file found: $R2 (paired-end mode)"
        fi
    else
        log_info "Single-end mode"
        
        # Single-end REQUIRES fragment length
        if [ -z "$FRAG_LEN" ] || [ -z "$FRAG_SD" ]; then
            log_error "Single-end mode requires fragment length parameters!"
            log_error "  Usage: $0 <index> <r1> - <output> <threads> <fragment_length> <fragment_sd>"
            log_error "  Example: $0 kallisto.idx sample.fq.gz - output/ 8 200 20"
            log_error ""
            log_error "  Typical values:"
            log_error "    RNA-seq: length=200, sd=20"
            log_error "    Small RNA: length=50, sd=10"
            ((errors++))
        else
            log_info "Fragment length: ${FRAG_LEN} ± ${FRAG_SD}"
        fi
    fi
    
    # Check output directory is writable
    OUTPUT_PARENT=$(dirname "$OUTPUT")
    if [ ! -d "$OUTPUT_PARENT" ]; then
        log_warn "Parent directory doesn't exist, creating: $OUTPUT_PARENT"
        mkdir -p "$OUTPUT_PARENT" || {
            log_error "Cannot create parent directory: $OUTPUT_PARENT"
            ((errors++))
        }
    fi
    
    # Check Kallisto is installed
    if ! command -v kallisto &> /dev/null; then
        log_error "Kallisto not found in PATH"
        log_error "  Install with: conda install -c bioconda kallisto"
        log_error "  Or load module: module load kallisto"
        ((errors++))
    else
        KALLISTO_VERSION=$(kallisto version 2>&1)
        log_info "Kallisto found: $KALLISTO_VERSION"
    fi
    
    # Exit if any validation errors
    if [ $errors -gt 0 ]; then
        log_error "Validation failed with $errors error(s)"
        exit 1
    fi
    
    log_info "All inputs validated successfully"
}

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================
INDEX="$1"
R1="$2"
R2="${3:-}"
OUTPUT="$4"
THREADS="${5:-8}"
FRAG_LEN="${6:-}"
FRAG_SD="${7:-}"

# =============================================================================
# DISPLAY PARAMETERS
# =============================================================================
echo ""
echo "🦖 RAPTOR v2.2.0 - Kallisto Quantification"
echo "═══════════════════════════════════════════════════════════════"
log_step "Parameters:"
echo "  Index:   $INDEX"
echo "  R1:      $R1"
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    echo "  R2:      $R2"
else
    echo "  R2:      (single-end)"
    if [ -n "$FRAG_LEN" ]; then
        echo "  Fragment: ${FRAG_LEN} ± ${FRAG_SD}"
    fi
fi
echo "  Output:  $OUTPUT"
echo "  Threads: $THREADS"
echo ""

# =============================================================================
# VALIDATE INPUTS (v2.2.0)
# =============================================================================
log_step "Validating inputs..."
validate_inputs "$@"
echo ""

# =============================================================================
# CREATE OUTPUT DIRECTORY
# =============================================================================
log_step "Creating output directory..."
mkdir -p "$OUTPUT"
log_info "Output directory: $OUTPUT"
echo ""

# =============================================================================
# BUILD KALLISTO COMMAND
# =============================================================================
log_step "Building Kallisto command..."

CMD="kallisto quant -i $INDEX -o $OUTPUT -t $THREADS"

# Add input files
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    # Paired-end
    CMD="$CMD $R1 $R2"
    log_info "Paired-end quantification"
else
    # Single-end (requires fragment length)
    CMD="$CMD --single -l $FRAG_LEN -s $FRAG_SD $R1"
    log_info "Single-end quantification (fragment: ${FRAG_LEN}±${FRAG_SD})"
fi

echo ""
log_info "Command: $CMD"
echo ""

# =============================================================================
# RUN KALLISTO
# =============================================================================
log_step "Running Kallisto quantification..."
echo ""

# Run with timing
START_TIME=$(date +%s)

if $CMD 2>&1; then
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo ""
    log_info "Kallisto quantification completed successfully"
    log_info "Elapsed time: ${ELAPSED}s"
    
    # Check output files
    ABUNDANCE_FILE="$OUTPUT/abundance.tsv"
    RUN_INFO_FILE="$OUTPUT/run_info.json"
    
    if [ -f "$ABUNDANCE_FILE" ]; then
        NUM_TRANSCRIPTS=$(tail -n +2 "$ABUNDANCE_FILE" | wc -l)
        log_info "Quantified $NUM_TRANSCRIPTS transcripts"
        
        FILE_SIZE=$(du -h "$ABUNDANCE_FILE" | cut -f1)
        log_info "Output file: $ABUNDANCE_FILE ($FILE_SIZE)"
    else
        log_error "Output file not created: $ABUNDANCE_FILE"
        exit 1
    fi
    
else
    EXIT_CODE=$?
    echo ""
    log_error "Kallisto quantification failed (exit code: $EXIT_CODE)"
    log_error "Common issues:"
    log_error "  • Index version mismatch (rebuild index with same Kallisto version)"
    log_error "  • Corrupted FASTQ files"
    log_error "  • Wrong fragment length for single-end"
    log_error "  • Index built from wrong organism"
    log_error "  • Insufficient memory (Kallisto needs ~4GB)"
    exit $EXIT_CODE
fi

# =============================================================================
# QUALITY CHECKS (v2.2.0)
# =============================================================================
log_step "Running quality checks..."

# Parse run_info.json if available
if [ -f "$RUN_INFO_FILE" ]; then
    # Extract metrics using python or jq if available
    if command -v python3 &> /dev/null; then
        python3 -c "
import json
try:
    with open('$RUN_INFO_FILE') as f:
        data = json.load(f)
    n_processed = data.get('n_processed', 0)
    n_pseudoaligned = data.get('n_pseudoaligned', 0)
    if n_processed > 0:
        rate = (n_pseudoaligned / n_processed) * 100
        print(f'  Processed: {n_processed:,} reads')
        print(f'  Pseudoaligned: {n_pseudoaligned:,} reads')
        print(f'  Rate: {rate:.2f}%')
        if rate < 50:
            print('  ⚠️  WARNING: Low pseudoalignment rate (<50%)')
        elif rate < 70:
            print('  ⚠️  Moderate pseudoalignment rate (50-70%)')
        else:
            print('  ✓ Good pseudoalignment rate (≥70%)')
except Exception as e:
    print(f'  Could not parse run_info.json: {e}')
"
    else
        log_info "run_info.json found (use Python/jq to parse metrics)"
    fi
else
    log_warn "run_info.json not found"
fi

# Check abundance.tsv has data
if [ -f "$ABUNDANCE_FILE" ]; then
    NUM_LINES=$(wc -l < "$ABUNDANCE_FILE")
    if [ $NUM_LINES -lt 2 ]; then
        log_warn "abundance.tsv has no data (empty quantification)"
    else
        log_info "abundance.tsv validated: $NUM_LINES lines"
    fi
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "═══════════════════════════════════════════════════════════════"
log_info "Kallisto quantification completed"
log_info "Output files:"
echo "  • abundance.tsv  : Transcript abundance estimates"
echo "  • abundance.h5   : HDF5 format (for sleuth)"
echo "  • run_info.json  : Run metadata and QC metrics"
echo ""
log_info "Next steps:"
echo "  1. Aggregate samples: raptor pipeline kallisto ..."
echo "  2. Or use abundance.tsv directly in R/Python"
echo "  3. For gene-level: provide --tx2gene mapping"
echo ""

exit 0
