#!/bin/bash
# =============================================================================
# RAPTOR Quick Profile Script - Version 2.1.0
# =============================================================================
# Fast profiling and ML-based recommendation workflow for RNA-seq data
# Features: Quality assessment, ML recommendations, resource estimation
# Usage: ./quick_profile.sh <counts.csv> [metadata.csv]
# =============================================================================

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

print_info() { echo -e "${BLUE}➜ $1${NC}"; }
print_success() { echo -e "${GREEN}✓ $1${NC}"; }
print_error() { echo -e "${RED}✗ $1${NC}"; }
print_warning() { echo -e "${YELLOW}⚠ $1${NC}"; }

# Help message
show_help() {
    cat << EOF
RAPTOR Quick Profile - Fast RNA-seq Data Profiling (v2.1.0)

Features:
  • Advanced quality assessment
  • ML-based pipeline recommendations  
  • Resource requirement estimation
  • Interactive HTML reports

Usage:
    ./quick_profile.sh <counts.csv> [metadata.csv] [options]

Arguments:
    counts.csv      Count matrix (genes x samples)
    metadata.csv    Sample metadata (optional)

Options:
    -o, --output DIR        Output directory (default: quick_profile_results/)
    -t, --threads N         Number of threads (default: 4)
    -h, --help              Show this help message

Examples:
    # Basic profiling
    ./quick_profile.sh my_counts.csv

    # With metadata
    ./quick_profile.sh my_counts.csv my_metadata.csv

    # Custom output location
    ./quick_profile.sh counts.csv -o my_results/

    # With more threads
    ./quick_profile.sh counts.csv -t 8

Requirements:
    - RAPTOR installed (pip install raptor-rnaseq)
    - Count matrix in CSV format
    - Python 3.8+

EOF
}

# Parse arguments
COUNTS_FILE=""
METADATA_FILE=""
OUTPUT_DIR="quick_profile_results"
THREADS=4

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            if [ -z "$COUNTS_FILE" ]; then
                COUNTS_FILE="$1"
            elif [ -z "$METADATA_FILE" ]; then
                METADATA_FILE="$1"
            else
                print_error "Unknown argument: $1"
                show_help
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate inputs
if [ -z "$COUNTS_FILE" ]; then
    print_error "Count matrix file is required"
    show_help
    exit 1
fi

if [ ! -f "$COUNTS_FILE" ]; then
    print_error "Count file not found: $COUNTS_FILE"
    exit 1
fi

if [ -n "$METADATA_FILE" ] && [ ! -f "$METADATA_FILE" ]; then
    print_error "Metadata file not found: $METADATA_FILE"
    exit 1
fi

# Check RAPTOR installation
if ! command -v raptor &> /dev/null; then
    print_error "RAPTOR is not installed"
    echo "Install with: pip install raptor-rnaseq"
    exit 1
fi

# Main workflow
echo ""
echo "╔═══════════════════════════════════════════╗"
echo "║   RAPTOR QUICK PROFILE WORKFLOW v2.1.0    ║"
echo "║   ML Recommendations & Quality Assessment ║"
echo "╚═══════════════════════════════════════════╝"
echo ""

# Step 1: Validate data
print_info "Step 1/4: Validating input data..."

python3 << EOF
import pandas as pd
import sys

try:
    counts = pd.read_csv('${COUNTS_FILE}', index_col=0)
    print(f"✓ Count matrix loaded: {counts.shape[0]} genes × {counts.shape[1]} samples")
    
    # Basic checks
    if counts.isnull().any().any():
        print("✗ ERROR: Count matrix contains missing values")
        sys.exit(1)
    
    if (counts < 0).any().any():
        print("✗ ERROR: Count matrix contains negative values")
        sys.exit(1)
    
    if not counts.dtypes.apply(lambda x: x.kind in 'iuf').all():
        print("⚠ WARNING: Count matrix contains non-numeric values")
    
    # Check if looks like normalized data
    if counts.max().max() < 100:
        print("⚠ WARNING: Values are very small. Are these normalized counts?")
        print("  RAPTOR requires RAW counts, not FPKM/TPM/normalized values")
    
    print(f"Library sizes: {counts.sum(axis=0).min():,.0f} - {counts.sum(axis=0).max():,.0f}")
    
except Exception as e:
    print(f"✗ ERROR validating count matrix: {e}")
    sys.exit(1)

# Validate metadata if provided
if '${METADATA_FILE}':
    try:
        meta = pd.read_csv('${METADATA_FILE}')
        print(f"✓ Metadata loaded: {meta.shape[0]} samples")
        
        # Check if sample names match
        count_samples = set(counts.columns)
        meta_samples = set(meta['sample']) if 'sample' in meta.columns else set(meta.iloc[:, 0])
        
        if not count_samples.issubset(meta_samples):
            missing = count_samples - meta_samples
            print(f"⚠ WARNING: {len(missing)} samples in counts not found in metadata")
            
    except Exception as e:
        print(f"✗ ERROR validating metadata: {e}")
        sys.exit(1)

print("\n✓ Data validation complete")
EOF

if [ $? -ne 0 ]; then
    print_error "Data validation failed"
    exit 1
fi

# Step 2: Run profiling
print_info "Step 2/4: Profiling RNA-seq data characteristics..."

PROFILE_CMD="raptor profile --counts ${COUNTS_FILE}"

if [ -n "$METADATA_FILE" ]; then
    PROFILE_CMD="$PROFILE_CMD --metadata ${METADATA_FILE}"
fi

PROFILE_CMD="$PROFILE_CMD --output ${OUTPUT_DIR} --threads ${THREADS} --verbose"

echo "Command: $PROFILE_CMD"
eval $PROFILE_CMD

if [ $? -ne 0 ]; then
    print_error "Profiling failed"
    exit 1
fi

print_success "Profiling complete!"

# Step 3: Display key results
print_info "Step 3/4: Extracting key metrics..."

python3 << EOF
import json
import os

results_file = '${OUTPUT_DIR}/recommendations.json'

if os.path.exists(results_file):
    with open(results_file) as f:
        data = json.load(f)
    
    profile = data.get('profile', {})
    
    print("\n" + "="*60)
    print("DATA CHARACTERISTICS")
    print("="*60)
    print(f"Samples:          {profile.get('n_samples', 'N/A')}")
    print(f"Genes:            {profile.get('n_genes', 'N/A'):,}")
    print(f"BCV:              {profile.get('bcv', 0):.3f} ({profile.get('bcv_category', 'N/A')})")
    print(f"Mean depth:       {profile.get('mean_depth', 0):,.0f} ({profile.get('depth_category', 'N/A')})")
    print(f"Zero inflation:   {profile.get('zero_inflation', 0):.1%}")
    print(f"Library size CV:  {profile.get('library_size_cv', 0):.2f}")
    
    recommendations = data.get('recommendations', [])
    
    if recommendations:
        print("\n" + "="*60)
        print("TOP 3 RECOMMENDED PIPELINES")
        print("="*60)
        
        for i, rec in enumerate(recommendations[:3], 1):
            print(f"\n{i}. {rec.get('pipeline_name', 'N/A')} (Pipeline {rec.get('pipeline_id', 'N/A')})")
            print(f"   Score:           {rec.get('score', 0):.2f}")
            print(f"   Confidence:      {rec.get('confidence', 'N/A')}")
            print(f"   Expected time:   {rec.get('expected_runtime', 0):.1f} hours")
            print(f"   Expected memory: {rec.get('expected_memory', 0):.0f} GB")
            print(f"   Reasoning:       {rec.get('reasoning', 'N/A')[:80]}...")
    
    print("\n" + "="*60)
    
else:
    print("⚠ Results file not found")
EOF

# Step 4: Generate summary
print_info "Step 4/4: Generating summary..."

cat > "${OUTPUT_DIR}/QUICK_SUMMARY.txt" << EOF
RAPTOR Quick Profile Summary
============================

Input Data:
-----------
Count matrix:  ${COUNTS_FILE}
Metadata:      ${METADATA_FILE:-Not provided}
Output dir:    ${OUTPUT_DIR}
Date:          $(date)

Analysis completed successfully!

Next Steps:
-----------
1. Open HTML report:
   ${OUTPUT_DIR}/raptor_profile_report.html

2. Run recommended pipeline:
   raptor run --pipeline <ID> --data your_fastq_dir/

3. Run full benchmark (if needed):
   raptor compare --data your_fastq_dir/ --output benchmark_results/

4. Read the recommendations and choose your pipeline

For more information:
- Documentation: https://github.com/AyehBlk/RAPTOR/tree/main/docs
- Tutorials: https://github.com/AyehBlk/RAPTOR/tree/main/docs/tutorials
- Support: ayehbolouki1988@gmail.com
EOF

print_success "Summary saved to: ${OUTPUT_DIR}/QUICK_SUMMARY.txt"

# Final message
echo ""
echo "╔═══════════════════════════════════════════╗"
echo "║          QUICK PROFILE COMPLETE!          ║"
echo "╚═══════════════════════════════════════════╝"
echo ""
print_info "Results saved to: ${OUTPUT_DIR}/"
print_info "To view the report:"
echo "  xdg-open ${OUTPUT_DIR}/raptor_profile_report.html  # Linux"
echo "  open ${OUTPUT_DIR}/raptor_profile_report.html      # macOS"
echo ""
print_success "Profile complete in $(date +%H:%M:%S)! 🦖"
