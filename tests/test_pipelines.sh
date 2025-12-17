#!/bin/bash
# =============================================================================
# RAPTOR v2.1.1 Pipeline Integration Tests
# =============================================================================
# Integration tests for all 8 RNA-seq analysis pipelines
# Tests installation, configuration, and basic functionality
# 🆕 Now includes ATO availability check
# =============================================================================

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'

print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_success() { echo -e "${GREEN}✓ $1${NC}"; }
print_info() { echo -e "${CYAN}➜ $1${NC}"; }
print_warning() { echo -e "${YELLOW}⚠ $1${NC}"; }
print_error() { echo -e "${RED}✗ $1${NC}"; }

# Test counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_SKIPPED=0

# Test result function
test_result() {
    if [ $1 -eq 0 ]; then
        print_success "$2"
        ((TESTS_PASSED++))
    else
        print_error "$2"
        ((TESTS_FAILED++))
    fi
}

# Main test header
print_header "RAPTOR v2.1.1 PIPELINE INTEGRATION TESTS"

echo "Test Suite: Pipeline Installation & Configuration"
echo "Date: $(date)"
echo "🆕 Includes Adaptive Threshold Optimizer (ATO) checks"
echo ""

# TEST 1: RAPTOR Installation
print_header "Test 1: RAPTOR Installation"

if command -v raptor &> /dev/null; then
    VERSION=$(raptor --version 2>&1 | head -n 1)
    print_success "RAPTOR is installed: $VERSION"
    ((TESTS_PASSED++))
    
    # Check if version is 2.1.1
    if echo "$VERSION" | grep -q "2.1.1"; then
        print_success "Version is 2.1.1"
        ((TESTS_PASSED++))
    else
        print_warning "Expected version 2.1.1"
        ((TESTS_SKIPPED++))
    fi
else
    print_error "RAPTOR is not installed"
    ((TESTS_FAILED++))
    print_info "Install with: pip install raptor-rnaseq"
    exit 1
fi

# TEST 2: Python Module Imports
print_header "Test 2: Python Module Imports"

python3 -c "from raptor.profiler import RNAseqDataProfiler" 2>/dev/null
test_result $? "raptor.profiler"

python3 -c "from raptor.recommender import PipelineRecommender" 2>/dev/null
test_result $? "raptor.recommender"

python3 -c "from raptor.ml_recommender import MLPipelineRecommender" 2>/dev/null
test_result $? "raptor.ml_recommender"

# 🆕 v2.1.1: Test ATO import
python3 -c "from raptor.threshold_optimizer import AdaptiveThresholdOptimizer, optimize_thresholds" 2>/dev/null
test_result $? "raptor.threshold_optimizer (v2.1.1)"

# TEST 3: Required Tools Check
print_header "Test 3: Bioinformatics Tools"

TOOLS="samtools fastqc"
for tool in $TOOLS; do
    if command -v $tool &> /dev/null; then
        test_result 0 "$tool is installed"
    else
        print_warning "$tool NOT installed (optional)"
        ((TESTS_SKIPPED++))
    fi
done

# Alignment tools (any one is sufficient)
print_info "Checking alignment tools (at least one required)..."
ALIGNMENT_TOOLS="STAR salmon kallisto hisat2"
ALIGNMENT_FOUND=0
for tool in $ALIGNMENT_TOOLS; do
    if command -v $tool &> /dev/null; then
        print_success "$tool is installed"
        ((ALIGNMENT_FOUND++))
    else
        print_warning "$tool NOT installed"
    fi
done

if [ $ALIGNMENT_FOUND -gt 0 ]; then
    test_result 0 "At least one alignment tool available"
else
    test_result 1 "No alignment tools found"
fi

# TEST 4: Python Dependencies
print_header "Test 4: Python Dependencies"

PYTHON_DEPS="numpy pandas scipy matplotlib seaborn sklearn"
for dep in $PYTHON_DEPS; do
    python3 -c "import $dep" 2>/dev/null
    if [ $? -eq 0 ]; then
        test_result 0 "Python: $dep"
    else
        test_result 1 "Python: $dep NOT installed"
    fi
done

# TEST 5: ATO Quick Functionality Test (v2.1.1)
print_header "Test 5: ATO Quick Test (v2.1.1)"

python3 << 'EOF'
import numpy as np
import pandas as pd
from raptor.threshold_optimizer import optimize_thresholds

# Generate test data
np.random.seed(42)
df = pd.DataFrame({
    'log2FoldChange': np.concatenate([
        np.random.normal(0, 0.2, 400),
        np.random.normal(2, 0.5, 50),
        np.random.normal(-2, 0.5, 50)
    ]),
    'pvalue': np.concatenate([
        np.random.uniform(0.1, 1, 400),
        np.random.exponential(0.001, 100)
    ])
})

# Run ATO
result = optimize_thresholds(df, goal='balanced', verbose=False)

# Verify result
assert result.logfc_cutoff > 0
assert result.n_significant_optimized >= 0
print(f"ATO test passed: |logFC| > {result.logfc_cutoff:.3f}, {result.n_significant_optimized} DE genes")
EOF

test_result $? "ATO functionality test"

# TEST 6: Pipeline Directories
print_header "Test 6: Pipeline Structure"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RAPTOR_ROOT="$(dirname "$SCRIPT_DIR")"

if [ -d "$RAPTOR_ROOT/pipelines" ]; then
    test_result 0 "Pipelines directory exists"
    
    # Count pipelines
    PIPELINE_COUNT=$(ls -d "$RAPTOR_ROOT/pipelines/pipeline"* 2>/dev/null | wc -l)
    if [ $PIPELINE_COUNT -ge 8 ]; then
        test_result 0 "All 8 pipelines present"
    else
        print_warning "Found $PIPELINE_COUNT pipelines (expected 8)"
        ((TESTS_SKIPPED++))
    fi
else
    print_warning "Pipelines directory not found at $RAPTOR_ROOT/pipelines"
    ((TESTS_SKIPPED++))
fi

# TEST SUMMARY
print_header "TEST SUMMARY"

echo ""
echo "Test Results:"
echo "  ✓ Passed:  $TESTS_PASSED"
echo "  ✗ Failed:  $TESTS_FAILED"
echo "  ⚠ Skipped: $TESTS_SKIPPED"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}═══════════════════════════════════════════════════════════${NC}"
    print_success "All critical tests passed! RAPTOR v2.1.1 is ready! 🦖"
    echo -e "${GREEN}═══════════════════════════════════════════════════════════${NC}"
    exit 0
else
    echo -e "${RED}═══════════════════════════════════════════════════════════${NC}"
    print_error "Some tests failed. Check errors above."
    echo -e "${RED}═══════════════════════════════════════════════════════════${NC}"
    exit 1
fi
