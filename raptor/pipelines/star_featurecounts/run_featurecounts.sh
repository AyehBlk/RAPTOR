#!/bin/bash
# =============================================================================
# RAPTOR v2.2.0 - featureCounts Quantification Script
# =============================================================================
#
# Enhanced bash script with v2.2.0 features:
# - Input validation with clear error messages
# - Better error handling
# - Progress indicators
# - Quality checks
#
# Usage: ./run_featurecounts.sh <gtf> <bam> <output> [threads] [strand] [paired]
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
    
    if [ \$# -lt 3 ]; then
        log_error "Insufficient arguments"
        echo "Usage: \$0 <gtf> <bam> <output> [threads] [strand] [paired]"
        exit 1
    fi
    
    # Check GTF
    if [ ! -f "\$GTF" ]; then
        log_error "GTF file not found: \$GTF"
        ((errors++))
    else
        log_info "GTF file found: \$GTF"
    fi
    
    # Check BAM
    if [ ! -f "\$BAM" ]; then
        log_error "BAM file not found: \$BAM"
        ((errors++))
    else
        log_info "BAM file found: \$BAM"
    fi
    
    # Validate strand
    if [ "\$STRAND" != "0" ] && [ "\$STRAND" != "1" ] && [ "\$STRAND" != "2" ]; then
        log_error "Invalid strand: \$STRAND (must be 0, 1, or 2)"
        ((errors++))
    else
        case \$STRAND in
            0) log_info "Strandedness: unstranded" ;;
            1) log_info "Strandedness: stranded (forward)" ;;
            2) log_info "Strandedness: stranded (reverse/dUTP)" ;;
        esac
    fi
    
    # Check featureCounts installed
    if ! command -v featureCounts &> /dev/null; then
        log_error "featureCounts not found in PATH"
        log_error "  Install: conda install -c bioconda subread"
        ((errors++))
    else
        FC_VERSION=\$(featureCounts -v 2>&1 | head -n1)
        log_info "featureCounts found: \$FC_VERSION"
    fi
    
    if [ \$errors -gt 0 ]; then
        log_error "Validation failed with \$errors error(s)"
        exit 1
    fi
}

GTF="\$1"
BAM="\$2"
OUTPUT="\$3"
THREADS="\${4:-8}"
STRAND="\${5:-0}"
PAIRED="\${6:-yes}"

echo ""
echo "🦖 RAPTOR v2.2.0 - featureCounts"
echo "═══════════════════════════════════════════════════════════════"
log_step "Parameters:"
echo "  GTF:     \$GTF"
echo "  BAM:     \$BAM"
echo "  Output:  \$OUTPUT"
echo "  Threads: \$THREADS"
echo "  Strand:  \$STRAND"
echo "  Paired:  \$PAIRED"
echo ""

log_step "Validating inputs..."
validate_inputs "\$@"
echo ""

log_step "Building featureCounts command..."
CMD="featureCounts -a \$GTF -o \$OUTPUT -t exon -g gene_id -T \$THREADS -s \$STRAND"

if [ "\$PAIRED" = "yes" ]; then
    CMD="\$CMD -p --countReadPairs"
    log_info "Counting read pairs (paired-end)"
else
    log_info "Counting reads (single-end)"
fi

CMD="\$CMD \$BAM"

echo ""
log_info "Command: \$CMD"
echo ""

log_step "Running featureCounts..."
START_TIME=\$(date +%s)

if \$CMD 2>&1; then
    END_TIME=\$(date +%s)
    ELAPSED=\$((END_TIME - START_TIME))
    echo ""
    log_info "featureCounts completed"
    log_info "Elapsed time: \${ELAPSED}s"
    
    # Check output
    if [ -f "\$OUTPUT" ]; then
        NUM_GENES=\$(tail -n +3 "\$OUTPUT" | wc -l)
        log_info "Quantified \$NUM_GENES genes"
        
        # Parse summary
        SUMMARY="\${OUTPUT}.summary"
        if [ -f "\$SUMMARY" ]; then
            ASSIGNED=\$(grep "Assigned" "\$SUMMARY" | awk '{print \$2}')
            TOTAL=\$(awk 'NR>1 {sum+=\$2} END {print sum}' "\$SUMMARY")
            if [ \$TOTAL -gt 0 ]; then
                RATE=\$(awk "BEGIN {printf \"%.2f\", (\$ASSIGNED/\$TOTAL)*100}")
                log_info "Assignment rate: \${RATE}%"
                
                if (( \$(echo "\$RATE < 50" | bc -l) )); then
                    log_warn "Low assignment rate (<50%)"
                    log_warn "  Check if GTF matches genome"
                    log_warn "  Check strandedness setting"
                fi
            fi
        fi
    else
        log_error "Output file not created"
        exit 1
    fi
else
    EXIT_CODE=\$?
    echo ""
    log_error "featureCounts failed (exit code: \$EXIT_CODE)"
    log_error "Common issues:"
    log_error "  • GTF and BAM chromosome names don't match"
    log_error "  • Wrong strandedness setting"
    log_error "  • Corrupted BAM file"
    exit \$EXIT_CODE
fi

echo ""
log_info "Output: \$OUTPUT"
exit 0
