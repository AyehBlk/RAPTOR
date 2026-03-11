#!/bin/bash
# =============================================================================
# RAPTOR v2.2.0 - Salmon Quantification Script
# =============================================================================
#
# Enhanced bash script with v2.2.0 features:
# - Input validation with clear error messages
# - Better error handling
# - Progress indicators
# - Quality checks
#
# Usage: ./run_salmon.sh <index> <r1> <r2|-> <output_dir> [threads] [library_type]
#
# Arguments:
#   index         : Path to Salmon index directory
#   r1            : Path to R1 FASTQ file
#   r2            : Path to R2 FASTQ file (or "-" for single-end)
#   output_dir    : Output directory for quantification
#   threads       : Number of threads (optional, default: 8)
#   library_type  : Library type (optional, default: A for auto-detect)
#
# Examples:
#   # Paired-end with auto-detect
#   ./run_salmon.sh salmon_index sample_R1.fq.gz sample_R2.fq.gz output/ 16
#
#   # Paired-end with ISR (dUTP stranded)
#   ./run_salmon.sh salmon_index sample_R1.fq.gz sample_R2.fq.gz output/ 16 ISR
#
#   # Single-end
#   ./run_salmon.sh salmon_index sample.fq.gz - output/ 8 A
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
        echo "Usage: $0 <index> <r1> <r2|-> <output_dir> [threads] [library_type]"
        echo ""
        echo "Examples:"
        echo "  Paired-end: $0 salmon_index R1.fq.gz R2.fq.gz output/ 16"
        echo "  Single-end: $0 salmon_index sample.fq.gz - output/ 8"
        exit 1
    fi
    
    # Check Salmon index exists
    if [ ! -d "$INDEX" ]; then
        log_error "Salmon index directory not found: $INDEX"
        log_error "  Build index with: salmon index -t transcripts.fa -i salmon_index"
        ((errors++))
    else
        # Check for index files
        if [ ! -f "$INDEX/info.json" ]; then
            log_error "Invalid Salmon index (missing info.json): $INDEX"
            log_error "  Rebuild index with: salmon index -t transcripts.fa -i salmon_index"
            ((errors++))
        else
            log_info "Salmon index validated: $INDEX"
        fi
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
    
    # Validate library type
    VALID_TYPES=("A" "IU" "ISF" "ISR" "MU" "MSF" "MSR" "OU" "OSF" "OSR" "SF" "SR" "U")
    if [[ ! " ${VALID_TYPES[@]} " =~ " ${LIBRARY_TYPE} " ]]; then
        log_error "Invalid library type: $LIBRARY_TYPE"
        log_error "  Valid types: ${VALID_TYPES[*]}"
        log_error "  Common options:"
        log_error "    A   = auto-detect (recommended)"
        log_error "    ISR = inward, stranded reverse (dUTP)"
        log_error "    ISF = inward, stranded forward (ligation)"
        log_error "    IU  = inward, unstranded"
        ((errors++))
    else
        log_info "Library type: $LIBRARY_TYPE"
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
    
    # Check Salmon is installed
    if ! command -v salmon &> /dev/null; then
        log_error "Salmon not found in PATH"
        log_error "  Install with: conda install -c bioconda salmon"
        log_error "  Or load module: module load salmon"
        ((errors++))
    else
        SALMON_VERSION=$(salmon --version 2>&1)
        log_info "Salmon found: $SALMON_VERSION"
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
LIBRARY_TYPE="${6:-A}"

# =============================================================================
# DISPLAY PARAMETERS
# =============================================================================
echo ""
echo "🦖 RAPTOR v2.2.0 - Salmon Quantification"
echo "═══════════════════════════════════════════════════════════════"
log_step "Parameters:"
echo "  Index:        $INDEX"
echo "  R1:           $R1"
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    echo "  R2:           $R2"
else
    echo "  R2:           (single-end)"
fi
echo "  Output:       $OUTPUT"
echo "  Threads:      $THREADS"
echo "  Library type: $LIBRARY_TYPE"
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
# BUILD SALMON COMMAND
# =============================================================================
log_step "Building Salmon command..."

CMD="salmon quant -i $INDEX -l $LIBRARY_TYPE -o $OUTPUT -p $THREADS"

# Add recommended flags
CMD="$CMD --validateMappings --gcBias --seqBias"

# Add input files
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    # Paired-end
    CMD="$CMD -1 $R1 -2 $R2"
    log_info "Paired-end quantification"
else
    # Single-end
    CMD="$CMD -r $R1"
    log_info "Single-end quantification"
fi

echo ""
log_info "Command: $CMD"
echo ""

# =============================================================================
# RUN SALMON
# =============================================================================
log_step "Running Salmon quantification..."
echo ""

# Run with timing
START_TIME=$(date +%s)

if $CMD 2>&1; then
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo ""
    log_info "Salmon quantification completed successfully"
    log_info "Elapsed time: ${ELAPSED}s"
    
    # Check output files
    QUANT_FILE="$OUTPUT/quant.sf"
    CMD_INFO="$OUTPUT/cmd_info.json"
    LIB_FORMAT="$OUTPUT/lib_format_counts.json"
    
    if [ -f "$QUANT_FILE" ]; then
        NUM_TRANSCRIPTS=$(tail -n +2 "$QUANT_FILE" | wc -l)
        log_info "Quantified $NUM_TRANSCRIPTS transcripts"
        
        FILE_SIZE=$(du -h "$QUANT_FILE" | cut -f1)
        log_info "Output file: $QUANT_FILE ($FILE_SIZE)"
    else
        log_error "Output file not created: $QUANT_FILE"
        exit 1
    fi
    
else
    EXIT_CODE=$?
    echo ""
    log_error "Salmon quantification failed (exit code: $EXIT_CODE)"
    log_error "Common issues:"
    log_error "  • Index version mismatch (rebuild index with same Salmon version)"
    log_error "  • Corrupted FASTQ files"
    log_error "  • Index built from wrong organism"
    log_error "  • Wrong library type (try -l A for auto-detect)"
    log_error "  • Insufficient memory (Salmon needs ~8GB)"
    exit $EXIT_CODE
fi

# =============================================================================
# QUALITY CHECKS (v2.2.0)
# =============================================================================
log_step "Running quality checks..."

# Parse library format counts if available
if [ -f "$LIB_FORMAT" ]; then
    if command -v python3 &> /dev/null; then
        python3 -c "
import json
try:
    with open('$LIB_FORMAT') as f:
        data = json.load(f)
    
    # Show library format detection
    if 'expected_format' in data:
        print('  Library format detected: ' + data['expected_format'])
    
    # Show read counts by orientation
    if 'read_files' in data:
        for read_file in data['read_files']:
            if 'num_compatible_fragments' in read_file:
                print(f\"  Compatible fragments: {read_file['num_compatible_fragments']:,}\")
            if 'num_assigned_fragments' in read_file:
                print(f\"  Assigned fragments: {read_file['num_assigned_fragments']:,}\")
except Exception as e:
    print(f'  Could not parse lib_format_counts.json: {e}')
"
    else
        log_info "lib_format_counts.json found (use Python to parse)"
    fi
fi

# Parse mapping rate from log if available
SALMON_LOG="$OUTPUT/logs/salmon_quant.log"
if [ -f "$SALMON_LOG" ]; then
    if grep -q "Mapping rate" "$SALMON_LOG"; then
        MAPPING_RATE=$(grep "Mapping rate" "$SALMON_LOG" | head -n1 | awk '{print $NF}')
        log_info "Mapping rate: $MAPPING_RATE"
        
        # Extract numeric value for comparison
        RATE_NUM=$(echo "$MAPPING_RATE" | tr -d '%')
        if (( $(echo "$RATE_NUM < 50" | bc -l 2>/dev/null || echo 0) )); then
            log_warn "Low mapping rate (<50%)"
            log_warn "  Check if index matches organism"
            log_warn "  Check FASTQ quality"
            log_warn "  Try different library type"
        elif (( $(echo "$RATE_NUM < 70" | bc -l 2>/dev/null || echo 0) )); then
            log_warn "Moderate mapping rate (50-70%)"
            log_warn "  This may be acceptable depending on your data"
        else
            log_info "Good mapping rate (≥70%)"
        fi
    fi
fi

# Check quant.sf has data
if [ -f "$QUANT_FILE" ]; then
    NUM_LINES=$(wc -l < "$QUANT_FILE")
    if [ $NUM_LINES -lt 2 ]; then
        log_warn "quant.sf has no data (empty quantification)"
    else
        log_info "quant.sf validated: $NUM_LINES lines"
    fi
fi

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "═══════════════════════════════════════════════════════════════"
log_info "Salmon quantification completed"
log_info "Output files:"
echo "  • quant.sf             : Transcript abundance estimates"
echo "  • quant.genes.sf       : Gene-level (if --geneMap provided)"
echo "  • aux_info/            : Auxiliary information"
echo "  • logs/                : Log files"
echo "  • lib_format_counts.json : Library format detection"
echo ""
log_info "Next steps:"
echo "  1. Aggregate samples: raptor pipeline salmon ..."
echo "  2. Or use quant.sf directly in R/Python"
echo "  3. For gene-level: provide --tx2gene mapping"
echo "  4. Use with sleuth/DESeq2 for differential expression"
echo ""

exit 0
