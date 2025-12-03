#!/bin/bash
# =============================================================================
# RAPTOR Full Benchmark Script - Version 2.1.0
# =============================================================================
# Comprehensive pipeline benchmarking workflow with:
# - Real-time resource monitoring
# - Parallel pipeline execution
# - ML-enhanced result analysis
# Compares multiple RNA-seq analysis pipelines on your data
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
RAPTOR Full Benchmark - Comprehensive Pipeline Comparison

Usage:
    ./full_benchmark.sh [options]

Options:
    -d, --data DIR          FASTQ data directory (required)
    -o, --output DIR        Output directory (default: benchmark_results/)
    -p, --pipelines LIST    Pipeline IDs to benchmark (default: 1,3,4,5)
    -t, --threads N         Number of threads (default: 8)
    -m, --memory SIZE       Memory limit (default: 32G)
    --mode MODE             Benchmark mode: quick, standard, full (default: standard)
    --parallel N            Run N pipelines in parallel (default: 1)
    --ground-truth FILE     Ground truth file for accuracy evaluation
    -h, --help              Show this help message

Benchmark Modes:
    quick       - 2-3 fastest pipelines, ~1-2 hours
    standard    - 4 recommended pipelines, ~4-6 hours
    full        - All 8 pipelines, ~12-24 hours

Pipeline Options:
    1  - STAR-RSEM-DESeq2 (highest accuracy)
    2  - HISAT2-StringTie-Ballgown (novel transcripts)
    3  - Salmon-edgeR (recommended, balanced)
    4  - Kallisto-Sleuth (fastest)
    5  - STAR-HTSeq-limma (complex designs)
    6  - STAR-featureCounts-NOISeq (small samples)
    7  - Bowtie2-RSEM-EBSeq (memory efficient)
    8  - HISAT2-Cufflinks-Cuffdiff (legacy)

Examples:
    # Quick benchmark with 3 pipelines
    ./full_benchmark.sh -d data/ --mode quick

    # Standard benchmark with custom pipelines
    ./full_benchmark.sh -d data/ -p 1,3,4,5 -t 16

    # Full benchmark with all 8 pipelines
    ./full_benchmark.sh -d data/ --mode full -t 32 -m 64G

    # Parallel execution
    ./full_benchmark.sh -d data/ -p 1,3,4,5 --parallel 2

    # With ground truth for accuracy evaluation
    ./full_benchmark.sh -d data/ --ground-truth truth.csv

Requirements:
    - RAPTOR installed
    - All bioinformatics tools installed (STAR, Salmon, etc.)
    - Reference genome indices built
    - Sufficient computational resources

Notes:
    - Full benchmarking is resource-intensive
    - Estimated time varies with data size
    - Results include runtime, memory, and accuracy metrics
    - Use --mode quick for testing your workflow first

EOF
}

# Default values
DATA_DIR=""
OUTPUT_DIR="benchmark_results"
PIPELINES="1,3,4,5"
THREADS=8
MEMORY="32G"
MODE="standard"
PARALLEL=1
GROUND_TRUTH=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -d|--data)
            DATA_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--pipelines)
            PIPELINES="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --ground-truth)
            GROUND_TRUTH="$2"
            shift 2
            ;;
        *)
            print_error "Unknown argument: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate inputs
if [ -z "$DATA_DIR" ]; then
    print_error "Data directory is required (-d/--data)"
    show_help
    exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
    print_error "Data directory not found: $DATA_DIR"
    exit 1
fi

# Check RAPTOR installation
if ! command -v raptor &> /dev/null; then
    print_error "RAPTOR is not installed"
    echo "Install with: pip install raptor-rnaseq"
    exit 1
fi

# Adjust pipelines based on mode
case $MODE in
    quick)
        PIPELINES="3,4"
        print_info "Quick mode: benchmarking pipelines 3,4"
        ;;
    standard)
        PIPELINES="1,3,4,5"
        print_info "Standard mode: benchmarking pipelines 1,3,4,5"
        ;;
    full)
        PIPELINES="1,2,3,4,5,6,7,8"
        print_info "Full mode: benchmarking all 8 pipelines"
        ;;
esac

# Header
echo ""
echo "╔═══════════════════════════════════════════════════════════╗"
echo "║      RAPTOR FULL BENCHMARK WORKFLOW v2.1.0                ║"
echo "║      Comprehensive RNA-seq Pipeline Comparison            ║"
echo "╚═══════════════════════════════════════════════════════════╝"
echo ""

# Configuration summary
echo "Configuration:"
echo "  Data directory:    $DATA_DIR"
echo "  Output directory:  $OUTPUT_DIR"
echo "  Pipelines:         $PIPELINES"
echo "  Threads:           $THREADS"
echo "  Memory:            $MEMORY"
echo "  Mode:              $MODE"
echo "  Parallel jobs:     $PARALLEL"
if [ -n "$GROUND_TRUTH" ]; then
    echo "  Ground truth:      $GROUND_TRUTH"
fi
echo ""

# Estimate runtime
print_info "Estimating benchmark runtime..."

python3 << EOF
pipelines = [int(x) for x in '$PIPELINES'.split(',')]
n_pipelines = len(pipelines)
parallel = int('$PARALLEL')

# Rough estimates (hours per pipeline for medium dataset)
times = {1: 3.5, 2: 2.0, 3: 0.5, 4: 0.3, 5: 3.0, 6: 3.0, 7: 2.5, 8: 4.0}

total_time = sum(times.get(p, 2.0) for p in pipelines)

if parallel > 1:
    estimated_time = total_time / parallel
else:
    estimated_time = total_time

print(f"\nEstimated runtime: {estimated_time:.1f} hours")
print(f"  ({n_pipelines} pipelines, {parallel} parallel job(s))")

if estimated_time > 12:
    print(f"\n⚠ WARNING: This will take a significant amount of time!")
    print(f"  Consider running overnight or using --mode quick for testing")
EOF

# User confirmation
read -p "Continue with benchmark? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_info "Benchmark cancelled"
    exit 0
fi

# Step 1: Pre-flight checks
print_info "Step 1/5: Pre-flight checks..."

# Check for FASTQ files
FASTQ_COUNT=$(find "$DATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" | wc -l)

if [ "$FASTQ_COUNT" -eq 0 ]; then
    print_warning "No FASTQ files found in $DATA_DIR"
    print_info "Looking for .fastq or .fq files..."
    FASTQ_COUNT=$(find "$DATA_DIR" -name "*.fastq" -o -name "*.fq" | wc -l)
    
    if [ "$FASTQ_COUNT" -eq 0 ]; then
        print_error "No FASTQ files found!"
        exit 1
    fi
fi

print_success "Found $FASTQ_COUNT FASTQ files"

# Check available resources
print_info "Checking system resources..."

TOTAL_MEM=$(free -g | awk '/^Mem:/{print $2}')
AVAIL_MEM=$(free -g | awk '/^Mem:/{print $7}')
CPU_COUNT=$(nproc)

echo "  Total memory:     ${TOTAL_MEM}GB"
echo "  Available memory: ${AVAIL_MEM}GB"
echo "  CPU cores:        $CPU_COUNT"

# Warn if resources are tight
MEM_NUM=$(echo $MEMORY | sed 's/G//')
if [ "$AVAIL_MEM" -lt "$MEM_NUM" ]; then
    print_warning "Available memory ($AVAIL_MEM GB) is less than requested ($MEMORY)"
fi

if [ "$THREADS" -gt "$CPU_COUNT" ]; then
    print_warning "Requested threads ($THREADS) exceeds available cores ($CPU_COUNT)"
fi

# Step 2: Create benchmark configuration
print_info "Step 2/5: Creating benchmark configuration..."

mkdir -p "$OUTPUT_DIR"

cat > "${OUTPUT_DIR}/benchmark_config.yaml" << EOF
# RAPTOR Benchmark Configuration
# Generated: $(date)

benchmarking:
  default_pipelines: [$(echo $PIPELINES | tr ',' ' ')]
  mode: $MODE
  
  metrics:
    runtime: true
    memory_usage: true
    cpu_usage: true
    accuracy: $([ -n "$GROUND_TRUTH" ] && echo "true" || echo "false")
    concordance: true
  
  ground_truth:
    available: $([ -n "$GROUND_TRUTH" ] && echo "true" || echo "false")
    file: $GROUND_TRUTH
    fdr_threshold: 0.05
    fc_threshold: 1.0
  
  parallel_pipelines: $([ "$PARALLEL" -gt 1 ] && echo "true" || echo "false")
  max_parallel: $PARALLEL
  timeout_hours: 48
  
  tracking:
    save_intermediate_results: true
    track_tool_versions: true
    save_command_logs: true

resources:
  default_threads: $THREADS
  default_memory_gb: $(echo $MEMORY | sed 's/G//')
  max_threads: $THREADS
  monitor_resources: true
EOF

print_success "Configuration saved to: ${OUTPUT_DIR}/benchmark_config.yaml"

# Step 3: Run benchmark
print_info "Step 3/5: Running benchmark (this will take a while)..."

START_TIME=$(date +%s)

# Build the command
BENCHMARK_CMD="raptor compare"
BENCHMARK_CMD="$BENCHMARK_CMD --data $DATA_DIR"
BENCHMARK_CMD="$BENCHMARK_CMD --output $OUTPUT_DIR"
BENCHMARK_CMD="$BENCHMARK_CMD --pipelines $PIPELINES"
BENCHMARK_CMD="$BENCHMARK_CMD --threads $THREADS"
BENCHMARK_CMD="$BENCHMARK_CMD --memory $MEMORY"
BENCHMARK_CMD="$BENCHMARK_CMD --config ${OUTPUT_DIR}/benchmark_config.yaml"
BENCHMARK_CMD="$BENCHMARK_CMD --mode $MODE"

if [ "$PARALLEL" -gt 1 ]; then
    BENCHMARK_CMD="$BENCHMARK_CMD --parallel --max-jobs $PARALLEL"
fi

if [ -n "$GROUND_TRUTH" ]; then
    BENCHMARK_CMD="$BENCHMARK_CMD --ground-truth $GROUND_TRUTH"
fi

BENCHMARK_CMD="$BENCHMARK_CMD --verbose"

echo "Running: $BENCHMARK_CMD"
echo ""

# Execute benchmark
eval $BENCHMARK_CMD 2>&1 | tee "${OUTPUT_DIR}/benchmark.log"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(( (ELAPSED % 3600) / 60 ))

print_success "Benchmark completed in ${HOURS}h ${MINUTES}m!"

# Step 4: Generate comparison report
print_info "Step 4/5: Generating comparison report..."

raptor report \
    --results "$OUTPUT_DIR" \
    --output "${OUTPUT_DIR}/benchmark_comparison.html" \
    --include-all

print_success "Report generated: ${OUTPUT_DIR}/benchmark_comparison.html"

# Step 5: Create summary
print_info "Step 5/5: Creating benchmark summary..."

python3 - "$OUTPUT_DIR" << 'PYEOF'
import json
import os
import sys

output_dir = sys.argv[1]
results_file = f"{output_dir}/benchmark_results.json"

if not os.path.exists(results_file):
    print("⚠ Results file not found")
    sys.exit(0)

with open(results_file) as f:
    results = json.load(f)

print("\n" + "="*70)
print("BENCHMARK RESULTS SUMMARY")
print("="*70)

# Sort pipelines by score
pipelines = results.get('pipelines', [])
pipelines_sorted = sorted(pipelines, key=lambda x: x.get('overall_score', 0), reverse=True)

print("\nPipeline Rankings:")
print("-"*70)
print(f"{'Rank':<6} {'Pipeline':<30} {'Runtime':<12} {'Memory':<10} {'Score':<8}")
print("-"*70)

for i, p in enumerate(pipelines_sorted, 1):
    name = p.get('name', 'N/A')
    runtime = p.get('runtime_hours', 0)
    memory = p.get('peak_memory_gb', 0)
    score = p.get('overall_score', 0)
    
    medals = {1: "🥇", 2: "🥈", 3: "🥉"}
    rank_symbol = medals.get(i, f"  {i}")
    
    print(f"{rank_symbol:<6} {name:<30} {runtime:>6.1f}h     {memory:>5.1f}GB    {score:.2f}")

print("\nKey Findings:")
print("-"*70)

# Best for accuracy
if any(p.get('accuracy_metrics') for p in pipelines):
    best_acc = max(pipelines, key=lambda x: x.get('accuracy_metrics', {}).get('f1_score', 0))
    print(f"✓ Best accuracy:  {best_acc.get('name', 'N/A')}")

# Fastest
fastest = min(pipelines, key=lambda x: x.get('runtime_hours', float('inf')))
print(f"✓ Fastest:        {fastest.get('name', 'N/A')} ({fastest.get('runtime_hours', 0):.1f}h)")

# Most memory efficient
mem_efficient = min(pipelines, key=lambda x: x.get('peak_memory_gb', float('inf')))
print(f"✓ Memory efficient: {mem_efficient.get('name', 'N/A')} ({mem_efficient.get('peak_memory_gb', 0):.1f}GB)")

# Best overall
best_overall = pipelines_sorted[0]
print(f"✓ Best overall:   {best_overall.get('name', 'N/A')}")

print("\n" + "="*70)

PYEOF

# Create text summary
cat > "${OUTPUT_DIR}/BENCHMARK_SUMMARY.txt" << EOF
RAPTOR Full Benchmark Summary
=============================

Benchmark Configuration:
-----------------------
Data directory:     $DATA_DIR
Pipelines tested:   $PIPELINES
Benchmark mode:     $MODE
Total runtime:      ${HOURS}h ${MINUTES}m
Threads used:       $THREADS
Memory limit:       $MEMORY
Parallel jobs:      $PARALLEL

Results:
--------
Detailed results:   ${OUTPUT_DIR}/benchmark_results.json
Comparison report:  ${OUTPUT_DIR}/benchmark_comparison.html
Full logs:          ${OUTPUT_DIR}/benchmark.log

Next Steps:
-----------
1. Open the HTML report to see visualizations
2. Review pipeline rankings and trade-offs
3. Choose optimal pipeline for your needs
4. Document your choice for reproducibility

For publication, include:
- Pipeline name and version
- Tool versions (in ${OUTPUT_DIR}/tool_versions.txt)
- Configuration file (${OUTPUT_DIR}/benchmark_config.yaml)
- Key metrics from benchmark results

Generated: $(date)
EOF

print_success "Summary saved to: ${OUTPUT_DIR}/BENCHMARK_SUMMARY.txt"

# Final message
echo ""
echo "╔═══════════════════════════════════════════════════════════╗"
echo "║            BENCHMARK COMPLETE!                            ║"
echo "╚═══════════════════════════════════════════════════════════╝"
echo ""
print_info "Total runtime: ${HOURS}h ${MINUTES}m"
print_info "Results saved to: ${OUTPUT_DIR}/"
echo ""
print_info "To view results:"
echo "  xdg-open ${OUTPUT_DIR}/benchmark_comparison.html  # Linux"
echo "  open ${OUTPUT_DIR}/benchmark_comparison.html      # macOS"
echo ""
print_info "Summary file:"
echo "  cat ${OUTPUT_DIR}/BENCHMARK_SUMMARY.txt"
echo ""
print_success "Benchmark complete! 🦖"
