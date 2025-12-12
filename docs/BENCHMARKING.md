# RAPTOR v2.1.0 Benchmarking Guide

Compare all 8 pipelines on your data with comprehensive metrics, ML-powered ranking, and real-time monitoring.

## Quick Start
```bash
# Basic benchmarking with ML ranking
raptor compare --data fastq/ --output results/ --use-ml

# With resource monitoring
raptor compare --data fastq/ \
  --pipelines 1,3,4 \
  --threads 16 \
  --memory 64G \
  --monitor-resources

# With quality checks
raptor compare --data fastq/ \
  --output results/ \
  --quality-check \
  --track-intermediate
```

## What's New in v2.1.0

⭐ **ML-Powered Ranking** - Automatic ranking based on ML models  
⭐ **Resource Monitoring** - Real-time CPU/memory tracking  
⭐ **Quality Assessment** - Automatic QC before/after analysis  
⭐ **Interactive Dashboard** - View results in web interface  
⭐ **Ensemble Integration** - Combine results automatically  
⭐ **Checkpoint System** - Resume interrupted benchmarks  

## What Gets Benchmarked

**Performance Metrics:**
- **Runtime**: Wallclock time for each pipeline
- **Memory**: Peak RAM usage (now tracked in real-time!)
- **CPU**: Core utilization (live monitoring)
- **Disk I/O**: Read/write speeds (NEW!)

**Quality Metrics:**
- **Accuracy**: If ground truth available
- **Concordance**: Agreement between pipelines
- **Sensitivity**: True positive rate
- **Precision**: False discovery rate

**NEW in v2.1.0:**
- **ML Score**: Machine learning prediction
- **Quality Score**: Data quality assessment
- **Resource Efficiency**: Performance per resource unit
- **Confidence Intervals**: Statistical validation

## Input Data

**FASTQ files:**
```
fastq/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz  
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

**Or simulated data:**
```bash
raptor simulate --output sim_data/ --size medium
raptor compare \
  --data sim_data/fastq/ \
  --ground-truth sim_data/ground_truth.csv \
  --use-ml
```

## Understanding Results

The benchmark generates:

**v2.0.0 outputs:**
- `benchmark_results.json`: Raw metrics
- `comparison_report.html`: Interactive visualizations
- `pipeline_*/`: Individual pipeline outputs

**NEW v2.1.0 outputs:**
- `ml_rankings.json`: ML-based rankings
- `resource_usage.csv`: Detailed resource tracking
- `quality_assessment.json`: QC results
- `ensemble_results/`: Combined analysis
- `interactive_dashboard/`: Launch with `raptor dashboard`

## Comparison Metrics

**Accuracy Metrics** (with ground truth):
- Sensitivity (True Positive Rate)
- Specificity (True Negative Rate)  
- Precision, F1 Score
- Concordance with truth
- **NEW**: ML confidence score

**Performance Metrics**:
- Runtime (hours)
- Peak memory (GB)
- **NEW**: Real-time CPU tracking
- **NEW**: Disk I/O statistics
- **NEW**: Resource efficiency score

**Quality Metrics (NEW):**
- Data quality score (0-100)
- Contamination check
- Batch effect detection
- Outlier identification

## Advanced Options

```bash
# Quick benchmark with ML (NEW)
raptor compare --data fastq/ --quick --use-ml

# Full benchmark with all features (NEW)
raptor compare --data fastq/ \
  --full \
  --use-ml \
  --monitor-resources \
  --quality-check \
  --ensemble

# Parallel execution with monitoring
raptor compare --data fastq/ \
  --parallel \
  --max-jobs 4 \
  --monitor-resources

# Resume interrupted benchmark (NEW)
raptor compare --data fastq/ --resume

# Custom config
raptor compare --data fastq/ --config my_config.yaml
```

## ML-Powered Ranking (NEW)

```bash
# Enable ML ranking
raptor compare --data fastq/ --use-ml

# ML will rank pipelines based on:
# - Historical success on similar data
# - Your data characteristics
# - Resource constraints
# - Expected accuracy
```

**ML Ranking Output:**
```
ML-Powered Pipeline Rankings
═══════════════════════════

 Pipeline 3: Salmon-edgeR
   ML Score: 0.89 (High Confidence)
   Expected Accuracy: 0.88
   Expected Runtime: 22 min
   
 Pipeline 1: STAR-RSEM-DESeq2
   ML Score: 0.85 (High Confidence)
   Expected Accuracy: 0.92
   Expected Runtime: 3.5 hours
   
 Pipeline 4: Kallisto-Sleuth
   ML Score: 0.82 (Medium Confidence)
   Expected Accuracy: 0.83
   Expected Runtime: 15 min
```

## Real-Time Monitoring (NEW)

Watch your benchmark in real-time:

```bash
# In terminal 1: Start benchmark with monitoring
raptor compare --data fastq/ --monitor-resources

# In terminal 2: Watch live stats
raptor monitor --live --benchmark results/

# Or use interactive dashboard
raptor dashboard --results results/
```

**Live monitoring shows:**
- Current pipeline running
- CPU/Memory usage graphs
- Estimated time remaining
- Bottleneck detection
- Resource alerts

## Quality Assessment (NEW)

Automatic quality checks before and after:

```bash
raptor compare --data fastq/ --quality-check
```

**Quality checks include:**
- Library size distribution
- Gene detection rate
- Zero inflation
- Batch effects
- Contamination
- Outlier samples

**QC Report:**
```
Quality Assessment Summary
═══════════════════════════

Pre-Analysis QC:
✅ Library sizes: Good (CV = 12%)
✅ Gene detection: Excellent (18,234 genes)
⚠️  Batch effect: Minor (PVR = 8%)
✅ No contamination detected

Post-Analysis Validation:
✅ All pipelines completed successfully
✅ Results consistent across pipelines
✅ No quality flags raised
```

## Interpreting Results

**Traditional Ranking** (v2.0.0):
1. Accuracy (if ground truth available)
2. Speed-accuracy tradeoff
3. Resource efficiency

**NEW ML-Enhanced Ranking** (v2.1.0):
1. ML prediction score (0-1)
2. Historical success rate
3. Confidence level
4. Resource efficiency
5. Quality compatibility

**Use the results** to:
- Validate ML recommendations
- Understand trade-offs
- Choose optimal pipeline
- Publish comprehensive comparisons
- **NEW**: Train custom ML models
- **NEW**: Generate ensemble results

## Interactive Dashboard (NEW)

View results in web interface:

```bash
# Launch dashboard after benchmarking
raptor dashboard --results results/

# Opens at http://localhost:8501
```

**Dashboard features:**
- Side-by-side pipeline comparison
- Interactive plots
- Resource usage graphs
- Parameter exploration
- Report generation

## Ensemble Analysis Integration (NEW)

Automatically combine results:

```bash
# Run benchmark with ensemble
raptor compare --data fastq/ --ensemble

# Or combine after benchmarking
raptor ensemble \
  --results results/ \
  --method vote \
  --min-agreement 2
```

**Ensemble gives you:**
- High-confidence genes (found by multiple pipelines)
- Consensus results
- Robustness validation

## Python API

```python
from raptor import PipelineBenchmark

# v2.1.0 enhanced API
benchmark = PipelineBenchmark(
    data_dir='fastq/',
    output_dir='results/',
    threads=8,
    monitor_resources=True,  # NEW
    use_ml_ranking=True,      # NEW
    quality_check=True        # NEW
)

# Run with monitoring
results = benchmark.run_pipelines(
    [1, 3, 4],
    track_intermediate=True,  # NEW
    enable_checkpoints=True   # NEW
)

# Get ML rankings (NEW)
ranked = benchmark.rank_results_ml(results)

# Save comprehensive results
benchmark.save_results(
    results,
    include_resources=True,   # NEW
    include_quality=True      # NEW
)
```

## Resume Interrupted Benchmarks (NEW)

```bash
# If benchmark is interrupted
raptor compare --data fastq/ --resume

# Or specify checkpoint
raptor compare --data fastq/ --resume-from checkpoint_20251119.pkl
```

**Checkpoints save:**
- Completed pipelines
- Intermediate results
- Resource usage so far
- Quality assessments

## Batch Processing (NEW)

Benchmark multiple datasets:

```bash
# Process multiple experiments
for dataset in exp1 exp2 exp3; do
  raptor compare \
    --data ${dataset}/fastq/ \
    --output results_${dataset}/ \
    --use-ml \
    --monitor-resources
done

# Compare across datasets
raptor meta-compare \
  --benchmarks results_*/\
  --output meta_comparison/
```

## Resource Prediction (NEW)

Estimate requirements before running:

```bash
# Predict resource needs
raptor predict-resources \
  --data fastq/ \
  --pipelines 1,3,4

# Output:
Estimated Requirements:
├─ Pipeline 1: 45 GB RAM, 3.5 hours, 32 cores
├─ Pipeline 3: 12 GB RAM, 22 min, 16 cores
└─ Pipeline 4: 6 GB RAM, 15 min, 8 cores

Total estimated time (sequential): 4.1 hours
Total estimated time (parallel): 3.5 hours
```

## Tips

**v2.0.0 tips still apply:**
- Start with quick benchmark (2-3 pipelines)
- Use simulated data for testing
- Run overnight for large datasets
- Save intermediate results
- Compare to profile recommendation

**NEW v2.1.0 tips:**
- ⭐ Always enable ML ranking (`--use-ml`)
- ⭐ Monitor resources to identify bottlenecks
- ⭐ Run quality checks to catch issues early
- ⭐ Use checkpoints for long benchmarks
- ⭐ View results in interactive dashboard
- ⭐ Consider ensemble for critical projects
- ⭐ Let ML guide which pipelines to benchmark

## Common Workflows

### Workflow 1: Quick Validation (NEW)

```bash
# Get ML recommendation
raptor profile --counts counts.csv --use-ml

# Benchmark top 3 ML recommendations
raptor compare \
  --data fastq/ \
  --pipelines <top-3-from-ML> \
  --quick \
  --monitor-resources

# View in dashboard
raptor dashboard --results results/
```

**Time: 2-4 hours**

### Workflow 2: Comprehensive Validation

```bash
# Full benchmark with all features
raptor compare \
  --data fastq/ \
  --full \
  --use-ml \
  --quality-check \
  --monitor-resources \
  --ensemble

# Generate publication report
raptor report \
  --results results/ \
  --type comprehensive \
  --output benchmark_report.html
```

**Time: 12-24 hours**

### Workflow 3: Critical Analysis (NEW)

```bash
# 1. Quality check first
raptor qc --counts counts.csv --metadata metadata.csv

# 2. ML-guided benchmark
raptor compare \
  --data fastq/ \
  --use-ml \
  --quality-check \
  --monitor-resources

# 3. Ensemble for maximum confidence
raptor ensemble \
  --results results/ \
  --method weighted \
  --min-agreement 3

# 4. Interactive exploration
raptor dashboard --results results/
```

**Time: 16-30 hours** (but highest confidence!)

## Troubleshooting

**Problem: Benchmark taking too long**

```bash
# Solution 1: Use ML to select fewer pipelines
raptor profile --counts counts.csv --use-ml
# Then benchmark only top 3

# Solution 2: Enable parallel execution
raptor compare --data fastq/ --parallel --max-jobs 4

# Solution 3: Use spot instances in cloud
raptor cloud compare --data fastq/ --spot-instances
```

**Problem: Out of memory**

```bash
# Check predicted requirements first (NEW)
raptor predict-resources --data fastq/ --pipelines 1,3,4

# Monitor in real-time
raptor compare --data fastq/ --monitor-resources

# Dashboard will show memory alerts
raptor dashboard --results results/
```

**Problem: Interrupted benchmark**

```bash
# Resume from checkpoint (NEW)
raptor compare --data fastq/ --resume

# Check what was completed
raptor benchmark-status --results results/
```

## See Also

- [ENSEMBLE.md](ENSEMBLE.md) - Combining pipeline results
- [RESOURCE_MONITORING.md](RESOURCE_MONITORING.md) - Resource tracking
- [QUALITY_ASSESSMENT.md](QUALITY_ASSESSMENT.md) - Data QC
- [ML_TRAINING.md](ML_TRAINING.md) - ML recommendations
- [DASHBOARD.md](DASHBOARD.md) - Interactive visualization
- [CLOUD_DEPLOYMENT.md](CLOUD_DEPLOYMENT.md) - Cloud benchmarking

For detailed metrics and interpretation, see the full documentation.

---

**RAPTOR v2.1.0 Benchmarking**  
**Author**: Ayeh Bolouki  
**License**: MIT

*"Benchmark smarter, not longer!"* ⚡🦖
