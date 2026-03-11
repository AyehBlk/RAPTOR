#!/bin/bash
# =============================================================================
# RAPTOR v2.2.0 - HISAT2 Alignment Script
# =============================================================================
#
# Enhanced bash script with v2.2.0 features:
# - Input validation with clear error messages
# - Better error handling
# - Progress indicators
# - Quality checks
#
# Usage: ./run_hisat2.sh <index> <r1> <r2> <output.sam> [threads]
#
# Arguments:
#   index    : Path to HISAT2 index (base name, not directory)
#   r1       : Path to R1 FASTQ file
#   r2       : Path to R2 FASTQ file (or "-" for single-end)
#   output   : Output SAM file path
#   threads  : Number of threads (optional, default: 8)
#
# Examples:
#   # Paired-end
#   ./run_hisat2.sh genome_index sample_R1.fq.gz sample_R2.fq.gz output.sam 16
#
#   # Single-end
#   ./run_hisat2.sh genome_index sample.fq.gz - output.sam 8
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
        echo "Usage: $0 <index> <r1> <r2|-> <output.sam> [threads]"
        exit 1
    fi
    
    # Check HISAT2 index exists
    if [ ! -f "${INDEX}.1.ht2" ] && [ ! -f "${INDEX}.1.ht2l" ]; then
        log_error "HISAT2 index not found: ${INDEX}"
        log_error "  Expected files: ${INDEX}.*.ht2 or ${INDEX}.*.ht2l"
        log_error "  Build index with: hisat2-build genome.fa genome_index"
        ((errors++))
    else
        log_info "HISAT2 index validated: ${INDEX}"
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
    fi
    
    # Check output directory is writable
    OUTPUT_DIR=$(dirname "$OUTPUT")
    if [ ! -d "$OUTPUT_DIR" ]; then
        log_warn "Output directory doesn't exist, creating: $OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR" || {
            log_error "Cannot create output directory: $OUTPUT_DIR"
            ((errors++))
        }
    fi
    
    if [ ! -w "$OUTPUT_DIR" ]; then
        log_error "Output directory not writable: $OUTPUT_DIR"
        ((errors++))
    fi
    
    # Check HISAT2 is installed
    if ! command -v hisat2 &> /dev/null; then
        log_error "HISAT2 not found in PATH"
        log_error "  Install with: conda install -c bioconda hisat2"
        log_error "  Or load module: module load hisat2"
        ((errors++))
    else
        HISAT2_VERSION=$(hisat2 --version 2>&1 | head -n1)
        log_info "HISAT2 found: $HISAT2_VERSION"
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

# =============================================================================
# DISPLAY PARAMETERS
# =============================================================================
echo ""
echo "🦖 RAPTOR v2.2.0 - HISAT2 Alignment"
echo "═══════════════════════════════════════════════════════════════"
log_step "Parameters:"
echo "  Index:   $INDEX"
echo "  R1:      $R1"
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    echo "  R2:      $R2"
else
    echo "  R2:      (single-end)"
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
# BUILD HISAT2 COMMAND
# =============================================================================
log_step "Building HISAT2 command..."

CMD="hisat2 -x $INDEX -p $THREADS --min-intronlen 20 --max-intronlen 500000"

# Add input files
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    CMD="$CMD -1 $R1 -2 $R2"
    log_info "Paired-end alignment"
else
    CMD="$CMD -U $R1"
    log_info "Single-end alignment"
fi

# Add output
CMD="$CMD -S $OUTPUT"

echo ""
log_info "Command: $CMD"
echo ""

# =============================================================================
# RUN HISAT2
# =============================================================================
log_step "Running HISAT2 alignment..."
echo ""

# Run with timing
START_TIME=$(date +%s)

if $CMD 2>&1; then
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo ""
    log_info "HISAT2 alignment completed successfully"
    log_info "Elapsed time: ${ELAPSED}s"
    
    # Check output file size
    if [ -f "$OUTPUT" ]; then
        OUTPUT_SIZE=$(du -h "$OUTPUT" | cut -f1)
        log_info "Output file: $OUTPUT ($OUTPUT_SIZE)"
        
        # Parse alignment stats from stderr (captured in log)
        # This is a basic check - real stats would come from stderr
        if [ -s "$OUTPUT" ]; then
            NUM_LINES=$(wc -l < "$OUTPUT")
            log_info "SAM file has $NUM_LINES lines"
        else
            log_warn "Output SAM file is empty"
            exit 1
        fi
    else
        log_error "Output file not created: $OUTPUT"
        exit 1
    fi
    
else
    EXIT_CODE=$?
    echo ""
    log_error "HISAT2 alignment failed (exit code: $EXIT_CODE)"
    log_error "Common issues:"
    log_error "  • Index version mismatch (rebuild index if needed)"
    log_error "  • Corrupted FASTQ files"
    log_error "  • Insufficient memory (HISAT2 needs ~16GB)"
    log_error "  • Wrong index path or missing index files"
    exit $EXIT_CODE
fi

# =============================================================================
# QUALITY CHECKS (v2.2.0)
# =============================================================================
log_step "Running quality checks..."

# Check if SAM file has aligned reads
NUM_MAPPED=$(grep -v "^@" "$OUTPUT" | grep -v "^\*" | wc -l || true)
TOTAL_READS=$(grep -v "^@" "$OUTPUT" | wc -l || true)

if [ $TOTAL_READS -gt 0 ]; then
    MAPPING_RATE=$(awk "BEGIN {printf \"%.2f\", ($NUM_MAPPED/$TOTAL_READS)*100}")
    log_info "Alignment rate: ${MAPPING_RATE}%"
    
    # Warn if low alignment rate
    if (( $(echo "$MAPPING_RATE < 50" | bc -l) )); then
        log_warn "Low alignment rate (<50%)"
        log_warn "  Check if index matches your organism"
        log_warn "  Check FASTQ quality"
    elif (( $(echo "$MAPPING_RATE < 70" | bc -l) )); then
        log_warn "Moderate alignment rate (50-70%)"
        log_warn "  This may be acceptable depending on your data"
    else
        log_info "Good alignment rate (≥70%)"
    fi
else
    log_warn "Cannot calculate alignment rate (empty SAM file)"
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "═══════════════════════════════════════════════════════════════"
log_info "HISAT2 alignment completed"
log_info "Next steps:"
echo "  1. Convert SAM to BAM: samtools view -bS $OUTPUT -o output.bam"
echo "  2. Sort BAM: samtools sort -o sorted.bam output.bam"
echo "  3. Index BAM: samtools index sorted.bam"
echo "  4. Run featureCounts for quantification"
echo ""

exit 0
