#!/bin/bash
# =============================================================================
# RAPTOR Demo Script
# =============================================================================
# Comprehensive demonstration of RAPTOR's main features
# Author: Ayeh Bolouki
# License: MIT
# =============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print section headers
print_header() {
    echo -e "\n${BLUE}================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================${NC}\n"
}

# Function to print success messages
print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# Function to print info messages
print_info() {
    echo -e "${YELLOW}➜ $1${NC}"
}

# Function to print error messages
print_error() {
    echo -e "${RED}✗ $1${NC}"
}

# Check if RAPTOR is installed
check_installation() {
    print_header "Checking RAPTOR Installation"
    
    if command -v raptor &> /dev/null; then
        print_success "RAPTOR is installed"
        raptor --version
    else
        print_error "RAPTOR is not installed or not in PATH"
        echo "Please install RAPTOR first:"
        echo "  pip install raptor-rnaseq"
        exit 1
    fi
}

# Setup demo directory
setup_demo_dir() {
    print_header "Setting Up Demo Directory"
    
    DEMO_DIR="raptor_demo_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$DEMO_DIR"
    cd "$DEMO_DIR"
    
    print_success "Created demo directory: $DEMO_DIR"
}

# Demo 1: Generate simulated data
demo_simulate() {
    print_header "Demo 1: Generating Simulated RNA-seq Data"
    
    print_info "Creating small simulated dataset..."
    raptor simulate \
        --output simulated_small/ \
        --size small \
        --n-genes 1000 \
        --n-samples 6 \
        --n-de 200 \
        --seed 42
    
    print_success "Simulated data created in simulated_small/"
    
    # Show what was created
    echo ""
    echo "Generated files:"
    ls -lh simulated_small/
    
    echo ""
    echo "Count matrix preview:"
    head -n 5 simulated_small/counts.csv
    
    echo ""
    echo "Metadata preview:"
    cat simulated_small/metadata.csv
}

# Demo 2: Data profiling
demo_profile() {
    print_header "Demo 2: Profiling RNA-seq Data"
    
    print_info "Analyzing data characteristics..."
    raptor profile \
        --counts simulated_small/counts.csv \
        --metadata simulated_small/metadata.csv \
        --output profile_results/ \
        --verbose
    
    print_success "Profiling complete!"
    
    echo ""
    echo "Generated reports:"
    ls -lh profile_results/
}

# Demo 3: Python API usage
demo_python_api() {
    print_header "Demo 3: Using RAPTOR Python API"
    
    print_info "Creating Python script example..."
    
    cat > api_example.py << 'EOF'
#!/usr/bin/env python3
"""Example of using RAPTOR Python API"""

import pandas as pd
from raptor import RNAseqDataProfiler, PipelineRecommender

# Load data
print("Loading data...")
counts = pd.read_csv('simulated_small/counts.csv', index_col=0)
metadata = pd.read_csv('simulated_small/metadata.csv')

print(f"Data shape: {counts.shape}")
print(f"Samples: {counts.shape[1]}, Genes: {counts.shape[0]}")

# Profile the data
print("\nProfiling data...")
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# Display key metrics
print("\n" + "="*50)
print("DATA PROFILE RESULTS")
print("="*50)
print(f"BCV: {profile['bcv']:.3f} ({profile['bcv_category']})")
print(f"Mean depth: {profile['mean_depth']:,.0f}")
print(f"Zero inflation: {profile['zero_inflation']:.1%}")

# Get recommendations
print("\nGetting pipeline recommendations...")
recommender = PipelineRecommender()
recommendations = recommender.recommend(profile, n=3)

# Display recommendations
print("\n" + "="*50)
print("TOP 3 PIPELINE RECOMMENDATIONS")
print("="*50)
for i, rec in enumerate(recommendations, 1):
    print(f"\n{i}. {rec['pipeline_name']}")
    print(f"   Score: {rec['score']:.2f}")
    print(f"   Expected runtime: {rec['expected_runtime']:.1f}h")
    print(f"   Expected memory: {rec['expected_memory']:.0f}GB")
EOF

    chmod +x api_example.py
    print_success "Python script created: api_example.py"
    
    if command -v python3 &> /dev/null; then
        print_info "Running Python API example..."
        python3 api_example.py || print_error "Python execution failed (RAPTOR module may not be installed)"
    fi
}

# Summary
show_summary() {
    print_header "Demo Complete - Summary"
    
    echo "This demo covered:"
    echo "  ✓ Data simulation"
    echo "  ✓ Data profiling"
    echo "  ✓ Pipeline recommendation"
    echo "  ✓ Python API usage"
    echo ""
    echo "All results are in: $(pwd)"
    echo ""
    echo "Next steps:"
    echo "  1. Open the HTML report to see visualizations"
    echo "  2. Try quick_profile.sh for faster profiling"
    echo "  3. Try full_benchmark.sh for comprehensive comparison"
    echo "  4. Read tutorials: docs/tutorials/"
    echo ""
    print_success "Thank you for using RAPTOR! 🦖"
}

# Main execution
main() {
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║          RAPTOR COMPREHENSIVE DEMO                       ║"
    echo "║   RNA-seq Analysis Pipeline Testing & Optimization       ║"
    echo "╚══════════════════════════════════════════════════════════╝"
    
    check_installation
    setup_demo_dir
    demo_simulate
    demo_profile
    demo_python_api
    show_summary
}

# Run main function
main
