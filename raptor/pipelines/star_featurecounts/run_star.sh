#!/bin/bash
# =============================================================================
# RAPTOR v2.2.0 - STAR Alignment Script
# =============================================================================
#
# Enhanced bash script with v2.2.0 features:
# - Input validation with clear error messages
# - Better error handling
# - Progress indicators
# - Quality checks
#
# Usage: ./run_star.sh <index> <r1> <r2|-> <output_prefix> [threads]
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

log_info() { echo -e "${GREEN}✓${NC} \$1"; }
log_warn() { echo -e "${YELLOW}⚠${NC} \$1"; }
log_error() { echo -e "${RED}✗${NC} \$1" >&2; }
log_step() { echo -e "${BLUE}▶${NC} \$1"; }

validate_inputs() {
    local errors=0
    
    if [ \$# -lt 4 ]; then
        log_error "Insufficient arguments"
        echo "Usage: \$0 <index> <r1> <r2|-> <output_prefix> [threads]"
        exit 1
    fi
    
    # Check STAR index
    if [ ! -d "\$INDEX" ]; then
        log_error "STAR index directory not found: \$INDEX"
        ((errors++))
    elif [ ! -f "\$INDEX/genomeParameters.txt" ]; then
        log_error "Invalid STAR index (missing genomeParameters.txt)"
        ((errors++))
    else
        log_info "STAR index validated: \$INDEX"
    fi
    
    # Check FASTQ files
    if [ ! -f "\$R1" ]; then
        log_error "R1 file not found: \$R1"
        ((errors++))
    else
        log_info "R1 file found: \$R1"
    fi
    
    if [ -n "\$R2" ] && [ "\$R2" != "-" ]; then
        if [ ! -f "\$R2" ]; then
            log_error "R2 file not found: \$R2"
            ((errors++))
        else
            log_info "R2 file found: \$R2 (paired-end)"
        fi
    else
        log_info "Single-end mode"
    fi
    
    # Check STAR installed
    if ! command -v STAR &> /dev/null; then
        log_error "STAR not found in PATH"
        log_error "  Install: conda install -c bioconda star"
        ((errors++))
    else
        STAR_VERSION=\$(STAR --version)
        log_info "STAR found: \$STAR_VERSION"
    fi
    
    if [ \$errors -gt 0 ]; then
        log_error "Validation failed with \$errors error(s)"
        exit 1
    fi
}

INDEX="\$1"
R1="\$2"
R2="\${3:-}"
PREFIX="\$4"
THREADS="\${5:-8}"

echo ""
echo "🦖 RAPTOR v2.2.0 - STAR Alignment"
echo "═══════════════════════════════════════════════════════════════"
log_step "Parameters:"
echo "  Index:  \$INDEX"
echo "  R1:     \$R1"
if [ -n "\$R2" ] && [ "\$R2" != "-" ]; then
    echo "  R2:     \$R2"
fi
echo "  Prefix: \$PREFIX"
echo "  Threads: \$THREADS"
echo ""

log_step "Validating inputs..."
validate_inputs "\$@"
echo ""

mkdir -p "\$(dirname \$PREFIX)"

log_step "Building STAR command..."
CMD="STAR --runThreadN \$THREADS --genomeDir \$INDEX --outFileNamePrefix \${PREFIX}_ \\
    --twopassMode Basic --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \\
    --outFilterMultimapNmax 20 --outFilterMismatchNmax 10"

# Handle gzipped input
if [[ "\$R1" == *.gz ]]; then
    CMD="\$CMD --readFilesCommand zcat"
fi

# Input files
if [ -n "\$R2" ] && [ "\$R2" != "-" ]; then
    CMD="\$CMD --readFilesIn \$R1 \$R2"
    log_info "Paired-end alignment"
else
    CMD="\$CMD --readFilesIn \$R1"
    log_info "Single-end alignment"
fi

echo ""
log_info "Command: \$CMD"
echo ""

log_step "Running STAR alignment..."
START_TIME=\$(date +%s)

if \$CMD 2>&1; then
    END_TIME=\$(date +%s)
    ELAPSED=\$((END_TIME - START_TIME))
    echo ""
    log_info "STAR alignment completed"
    log_info "Elapsed time: \${ELAPSED}s"
    
    # Check output
    BAM_FILE="\${PREFIX}_Aligned.sortedByCoord.out.bam"
    if [ -f "\$BAM_FILE" ]; then
        SIZE=\$(du -h "\$BAM_FILE" | cut -f1)
        log_info "BAM file: \$BAM_FILE (\$SIZE)"
    fi
    
    # Parse log
    LOG_FILE="\${PREFIX}_Log.final.out"
    if [ -f "\$LOG_FILE" ]; then
        UNIQUE_RATE=\$(grep "Uniquely mapped reads %" "\$LOG_FILE" | awk '{print \$NF}')
        log_info "Unique mapping rate: \$UNIQUE_RATE"
    fi
else
    EXIT_CODE=\$?
    echo ""
    log_error "STAR alignment failed (exit code: \$EXIT_CODE)"
    log_error "Common issues:"
    log_error "  • Index version mismatch"
    log_error "  • Insufficient memory (need 32GB+)"
    log_error "  • Corrupted FASTQ files"
    log_error "  • Disk space full"
    exit \$EXIT_CODE
fi

echo ""
log_info "Next: Run featureCounts on BAM file"
exit 0
