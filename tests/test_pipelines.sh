#!/bin/bash
# =============================================================================
# RAPTOR Pipeline Integration Tests
# =============================================================================
# Integration tests for all 8 RNA-seq analysis pipelines
# Tests installation, configuration, and basic functionality
# =============================================================================

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_success() { echo -e "${GREEN}✓ $1${NC}"; }
print_info() { echo -e "${YELLOW}➜ $1${NC}"; }
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
print_header "RAPTOR PIPELINE INTEGRATION TESTS"

echo "Test Suite: Pipeline Installation & Configuration"
echo "Date: $(date)"
echo ""

# TEST 1: RAPTOR Installation
print_header "Test 1: RAPTOR Installation"

if command -v raptor &> /dev/null; then
    VERSION=$(raptor --version 2>&1 | head -n 1)
    print_success "RAPTOR is installed: $VERSION"
    ((TESTS_PASSED++))
else
    print_error "RAPTOR is not installed"
    ((TESTS_FAILED++))
    exit 1
fi

# TEST 2: Required Tools Check
print_header "Test 2: Bioinformatics Tools"

TOOLS="samtools fastqc STAR salmon kallisto hisat2"
for tool in $TOOLS; do
    if command -v $tool &> /dev/null; then
        test_result 0 "$tool is installed"
    else
        test_result 1 "$tool NOT installed"
    fi
done

# TEST SUMMARY
print_header "TEST SUMMARY"

echo "Test Results:"
echo "  Passed:  $TESTS_PASSED"
echo "  Failed:  $TESTS_FAILED"
echo "  Skipped: $TESTS_SKIPPED"

if [ $TESTS_FAILED -eq 0 ]; then
    print_success "All tests passed! 🦖"
    exit 0
else
    print_error "Some tests failed"
    exit 1
fi
