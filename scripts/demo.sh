#!/bin/bash
# =============================================================================
# RAPTOR Demo Script - Version 2.1.1
# =============================================================================
# Comprehensive demonstration of RAPTOR's main features including:
# - ML-based pipeline recommendations
# - Interactive dashboard
# - Advanced quality assessment
# - Real-time resource monitoring
# - 🆕 Adaptive Threshold Optimizer (ATO)
#
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
    echo -e "${GREEN}✔ $1${NC}"
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

# Demo 2: Data profiling with ML recommendations
demo_profile() {
    print_header "Demo 2: Profiling RNA-seq Data with ML Recommendations"
    
    print_info "Analyzing data characteristics and getting ML-based recommendations..."
    raptor profile \
        --counts simulated_small/counts.csv \
        --metadata simulated_small/metadata.csv \
        --output profile_results/ \
        --use-ml \
        --verbose
    
    print_success "Profiling complete with ML recommendations!"
    
    echo ""
    echo "Generated reports:"
    ls -lh profile_results/
    
    echo ""
    echo "Check the recommendations.json file for ML-based pipeline suggestions"
}

# Demo 3: Threshold Optimization (NEW in v2.1.1)
demo_threshold_optimizer() {
    print_header "Demo 3: Adaptive Threshold Optimizer (ATO) - NEW in v2.1.1!"
    
    print_info "Demonstrating data-driven threshold optimization..."
    
    # Check if DE results exist, otherwise create demo
    if [ -f "simulated_small/de_results.csv" ]; then
        DE_FILE="simulated_small/de_results.csv"
    else
        print_info "Creating demo DE results..."
        
        # Create a simple DE results file for demonstration
        python3 << 'EOF'
import numpy as np
import pandas as pd

np.random.seed(42)
n_genes = 1000
n_de = 150

# Null genes
null_logfc = np.random.normal(0, 0.2, n_genes - n_de)
null_pval = np.random.uniform(0.05, 1, n_genes - n_de)

# DE genes
de_logfc = np.concatenate([
    np.random.normal(1.5, 0.5, n_de // 2),
    np.random.normal(-1.5, 0.5, n_de - n_de // 2)
])
de_pval = np.random.exponential(0.001, n_de)
de_pval = np.clip(de_pval, 1e-300, 0.05)

df = pd.DataFrame({
    'log2FoldChange': np.concatenate([null_logfc, de_logfc]),
    'pvalue': np.concatenate([null_pval, de_pval]),
    'baseMean': np.random.exponential(1000, n_genes)
})
df.index = [f'Gene_{i}' for i in range(n_genes)]
df.to_csv('demo_de_results.csv')
print(f"Created demo_de_results.csv with {n_genes} genes ({n_de} DE)")
EOF
        DE_FILE="demo_de_results.csv"
    fi
    
    print_info "Running threshold optimization..."
    raptor optimize-thresholds \
        --input "$DE_FILE" \
        --goal balanced \
        --output ato_results/ \
        --verbose
    
    print_success "Threshold optimization complete!"
    
    echo ""
    echo "Generated files:"
    ls -lh ato_results/
    
    # Show comparison of goals
    print_info "Comparing different analysis goals..."
    raptor optimize-thresholds \
        --input "$DE_FILE" \
        --compare-goals \
        --output ato_comparison/
    
    print_success "Goal comparison complete!"
    
    echo ""
    echo "ATO Results:"
    if [ -f "ato_results/ato_results.json" ]; then
        python3 -c "import json; d=json.load(open('ato_results/ato_results.json')); print(f\"  LogFC cutoff: {d['optimization_results']['logfc_cutoff']:.3f}\"); print(f\"  Padj cutoff: {d['optimization_results']['padj_cutoff']}\"); print(f\"  Traditional DE: {d['gene_counts']['traditional']}\"); print(f\"  Optimized DE: {d['gene_counts']['optimized']}\")"
    fi
}

# Demo 4: Python API usage
demo_python_api() {
    print_header "Demo 4: Using RAPTOR Python API"
    
    print_info "Creating Python script example..."
    
    cat > api_example.py << 'EOF'
#!/usr/bin/env python3
"""Example of using RAPTOR Python API including ATO (v2.1.1)"""

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

# NEW in v2.1.1: Threshold Optimization
print("\n" + "="*50)
print("ADAPTIVE THRESHOLD OPTIMIZER (v2.1.1)")
print("="*50)

try:
    from raptor.threshold_optimizer import optimize_thresholds
    
    # Load DE results if available
    import os
    if os.path.exists('demo_de_results.csv'):
        de_results = pd.read_csv('demo_de_results.csv', index_col=0)
        
        result = optimize_thresholds(de_results, goal='balanced')
        
        print(f"Optimized |logFC| cutoff: {result.logfc_cutoff:.3f}")
        print(f"Optimized padj cutoff: {result.padj_cutoff}")
        print(f"Traditional DE genes: {result.n_significant_traditional}")
        print(f"Optimized DE genes: {result.n_significant_optimized}")
        print(f"\nMethods text for publication:")
        print(result.methods_text[:200] + "...")
    else:
        print("No DE results available. Run DE analysis first.")
        
except ImportError:
    print("ATO module not available in this installation")
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
    echo "  ✔ Data simulation"
    echo "  ✔ Data profiling with ML recommendations"
    echo "  ✔ Pipeline recommendation using machine learning"
    echo "  ✔ 🆕 Adaptive Threshold Optimizer (ATO)"
    echo "  ✔ Python API usage"
    echo ""
    echo "All results are in: $(pwd)"
    echo ""
    echo "🆕 New in v2.1.1:"
    echo "  • Adaptive Threshold Optimizer (ATO)"
    echo "  • Data-driven logFC and p-value thresholds"
    echo "  • Publication-ready methods text generation"
    echo "  • Goal-based optimization (discovery/balanced/validation)"
    echo ""
    echo "Previous features (v2.1.0):"
    echo "  • ML-based pipeline recommendations"
    echo "  • Interactive dashboard"
    echo "  • Advanced quality assessment"
    echo "  • Real-time resource monitoring"
    echo "  • Ensemble analysis methods"
    echo ""
    echo "Next steps:"
    echo "  1. Open the HTML report to see visualizations"
    echo "  2. Try quick_profile.sh for faster profiling"
    echo "  3. Try full_benchmark.sh for comprehensive comparison"
    echo "  4. Run: python 11_threshold_optimizer.py --demo (ATO demo)"
    echo "  5. Run: python example_threshold_optimizer.py (ATO examples)"
    echo "  6. Read tutorials: docs/tutorials/"
    echo ""
    print_success "Thank you for using RAPTOR v2.1.1! 🦖"
}

# Main execution
main() {
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║      RAPTOR COMPREHENSIVE DEMO - v2.1.1                  ║"
    echo "║   RNA-seq Analysis Pipeline Testing & Optimization       ║"
    echo "║   🆕 Now with Adaptive Threshold Optimizer!              ║"
    echo "╚══════════════════════════════════════════════════════════╝"
    
    check_installation
    setup_demo_dir
    demo_simulate
    demo_profile
    demo_threshold_optimizer  # NEW in v2.1.1
    demo_python_api
    show_summary
}

# Run main function
main
