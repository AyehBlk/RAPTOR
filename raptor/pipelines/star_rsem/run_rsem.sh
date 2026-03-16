#!/bin/bash
# =============================================================================
# RAPTOR v2.2.0 - RSEM Quantification Script
# =============================================================================
#
# Enhanced bash script with v2.2.0 features:
# - Input validation with clear error messages
# - Better error handling
# - Progress indicators
# - Quality checks
#
# Usage: ./run_rsem.sh <rsem_ref> <r1> <r2|-> <output_prefix> [threads] [strandedness]
#
# Arguments:
#   rsem_ref         : Path to RSEM reference prefix
#   r1               : Path to R1 FASTQ file
#   r2               : Path to R2 FASTQ file (or "-" for single-end)
#   output_prefix    : Output prefix for RSEM results
#   threads          : Number of threads (optional, default: 8)
#   strandedness     : none, forward, or reverse (optional, default: none)
#
# Examples:
#   # Paired-end unstranded
#   ./run_rsem.sh rsem_ref/genome sample_R1.fq.gz sample_R2.fq.gz output/sample 16 none
#
#   # Paired-end reverse stranded (dUTP)
#   ./run_rsem.sh rsem_ref/genome sample_R1.fq.gz sample_R2.fq.gz output/sample 16 reverse
#
#   # Single-end
#   ./run_rsem.sh rsem_ref/genome sample.fq.gz - output/sample 8 none
#
# Author: Ayeh Bolouki
# Version: 2.2.0
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${GREEN}✓${NC} $1"; }
log_warn() { echo -e "${YELLOW}⚠${NC} $1"; }
log_error() { echo -e "${RED}✗${NC} $1" >&2; }
log_step() { echo -e "${BLUE}▶${NC} $1"; }

validate_inputs() {
    local errors=0
    
    if [ $# -lt 4 ]; then
        log_error "Insufficient arguments"
        echo "Usage: $0 <rsem_ref> <r1> <r2|-> <output_prefix> [threads] [strandedness]"
        echo ""
        echo "Examples:"
        echo "  Paired-end: $0 rsem_ref/genome R1.fq.gz R2.fq.gz output/sample 16 reverse"
        echo "  Single-end: $0 rsem_ref/genome sample.fq.gz - output/sample 8 none"
        exit 1
    fi
    
    # Check RSEM reference exists
    if [ ! -f "${REF}.grp" ] && [ ! -f "${REF}.ti" ]; then
        log_error "RSEM reference not found: $REF"
        log_error "  Expected files: ${REF}.grp, ${REF}.ti, ${REF}.transcripts.fa"
        log_error "  Build reference with:"
        log_error "    rsem-prepare-reference --gtf annotation.gtf genome.fa rsem_ref/genome"
        ((errors++))
    else
        log_info "RSEM reference validated: $REF"
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
    
    # Validate strandedness
    case "$STRAND" in
        none|forward|reverse)
            log_info "Strandedness: $STRAND"
            ;;
        *)
            log_error "Invalid strandedness: $STRAND"
            log_error "  Must be: none, forward, or reverse"
            log_error "  none    = unstranded (most common)"
            log_error "  forward = first read forward strand"
            log_error "  reverse = first read reverse strand (dUTP)"
            ((errors++))
            ;;
    esac
    
    # Check rsem-calculate-expression installed
    if ! command -v rsem-calculate-expression &> /dev/null; then
        log_error "RSEM not found in PATH"
        log_error "  Install: conda install -c bioconda rsem"
        ((errors++))
    else
        RSEM_VERSION=$(rsem-calculate-expression --version 2>&1 | head -n1 || echo "unknown")
        log_info "RSEM found: $RSEM_VERSION"
    fi
    
    if [ $errors -gt 0 ]; then
        log_error "Validation failed with $errors error(s)"
        exit 1
    fi
}

REF="$1"
R1="$2"
R2="${3:-}"
PREFIX="$4"
THREADS="${5:-8}"
STRAND="${6:-none}"

echo ""
echo "🦖 RAPTOR v2.2.0 - RSEM Quantification"
echo "═══════════════════════════════════════════════════════════════"
log_step "Parameters:"
echo "  Reference: $REF"
echo "  R1:        $R1"
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    echo "  R2:        $R2"
else
    echo "  R2:        (single-end)"
fi
echo "  Prefix:    $PREFIX"
echo "  Threads:   $THREADS"
echo "  Strandedness: $STRAND"
echo ""

log_step "Validating inputs..."
validate_inputs "$@"
echo ""

mkdir -p "$(dirname $PREFIX)"

log_step "Building RSEM command..."

# Base command
CMD="rsem-calculate-expression -p $THREADS --star --star-output-genome-bam \
    --estimate-rspd --seed 12345"

# Determine forward probability based on strandedness
case "$STRAND" in
    forward)
        PROB="1.0"
        log_info "Forward strand (forward-prob=1.0)"
        ;;
    reverse)
        PROB="0.0"
        log_info "Reverse strand / dUTP (forward-prob=0.0)"
        ;;
    *)
        PROB="0.5"
        log_info "Unstranded (forward-prob=0.5)"
        ;;
esac

CMD="$CMD --forward-prob $PROB"

# Paired-end or single-end
if [ -n "$R2" ] && [ "$R2" != "-" ]; then
    CMD="$CMD --paired-end $R1 $R2"
    log_info "Paired-end quantification"
else
    CMD="$CMD $R1"
    log_info "Single-end quantification"
    log_warn "Consider providing --fragment-length-mean and --fragment-length-sd"
fi

# Add reference and output prefix
CMD="$CMD $REF $PREFIX"

echo ""
log_info "Command: $CMD"
echo ""

log_step "Running RSEM quantification (this may take a while)..."
START_TIME=$(date +%s)

if $CMD 2>&1; then
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    echo ""
    log_info "RSEM quantification completed"
    log_info "Elapsed time: ${ELAPSED}s"
    
    # Check outputs
    GENE_RESULTS="${PREFIX}.genes.results"
    ISOFORM_RESULTS="${PREFIX}.isoforms.results"
    BAM_FILE="${PREFIX}.STAR.genome.bam"
    
    if [ -f "$GENE_RESULTS" ]; then
        NUM_GENES=$(tail -n +2 "$GENE_RESULTS" | wc -l)
        log_info "Quantified $NUM_GENES genes"
        
        SIZE=$(du -h "$GENE_RESULTS" | cut -f1)
        log_info "Gene results: $GENE_RESULTS ($SIZE)"
    else
        log_error "Gene results not found: $GENE_RESULTS"
        exit 1
    fi
    
    if [ -f "$ISOFORM_RESULTS" ]; then
        NUM_ISOFORMS=$(tail -n +2 "$ISOFORM_RESULTS" | wc -l)
        log_info "Quantified $NUM_ISOFORMS isoforms"
        
        SIZE=$(du -h "$ISOFORM_RESULTS" | cut -f1)
        log_info "Isoform results: $ISOFORM_RESULTS ($SIZE)"
    fi
    
    if [ -f "$BAM_FILE" ]; then
        SIZE=$(du -h "$BAM_FILE" | cut -f1)
        log_info "BAM file: $BAM_FILE ($SIZE)"
    fi
    
    # Parse alignment rate from .cnt file
    CNT_FILE="${PREFIX}.stat/${PREFIX##*/}.cnt"
    if [ -f "$CNT_FILE" ]; then
        if command -v python3 &> /dev/null; then
            python3 << PYEND
import sys
try:
    with open("$CNT_FILE") as f:
        lines = f.readlines()
        if len(lines) >= 4:
            n_total = int(lines[0].strip())
            n_aligned = int(lines[1].strip())
            if n_total > 0:
                rate = (n_aligned / n_total) * 100
                print(f"  Alignment rate: {rate:.2f}%")
                if rate < 50:
                    print(f"  ⚠ Low alignment rate (<50%)")
                    print(f"    Check if reference matches organism")
                    print(f"    Check strandedness setting")
except Exception as e:
    sys.stderr.write(f"Could not parse alignment rate: {e}\n")
PYEND
        fi
    fi
    
else
    EXIT_CODE=$?
    echo ""
    log_error "RSEM quantification failed (exit code: $EXIT_CODE)"
    log_error "Common issues:"
    log_error "  • RSEM reference mismatch (rebuild with same version)"
    log_error "  • Corrupted FASTQ files"
    log_error "  • Insufficient memory (need 32GB+)"
    log_error "  • Wrong strandedness setting"
    log_error "  • STAR version mismatch with RSEM reference"
    exit $EXIT_CODE
fi

echo ""
echo "═══════════════════════════════════════════════════════════════"
log_info "RSEM quantification completed"
log_info "Output files:"
echo "  • ${PREFIX}.genes.results       : Gene-level counts"
echo "  • ${PREFIX}.isoforms.results    : Transcript-level counts"
echo "  • ${PREFIX}.STAR.genome.bam     : Genome alignment BAM"
echo "  • ${PREFIX}.transcript.bam      : Transcript alignment BAM"
echo "  • ${PREFIX}.stat/               : Statistics directory"
echo ""
log_info "Next steps:"
echo "  1. Aggregate samples into count matrix"
echo "  2. Use with DESeq2/edgeR for gene-level DE"
echo "  3. Use with sleuth/EBSeq for isoform-level DE"
echo ""

exit 0
