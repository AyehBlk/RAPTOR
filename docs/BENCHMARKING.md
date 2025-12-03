# RAPTOR Benchmarking Guide

Compare all 8 pipelines on your data with comprehensive metrics.

## Quick Start
```bash
# Basic benchmarking
raptor compare --data fastq/ --output results/

# With custom settings
raptor compare --data fastq/ \
  --pipelines 1,3,4 \
  --threads 16 \
  --memory 64G
```

## What Gets Benchmarked

- **Runtime**: Wallclock time for each pipeline
- **Memory**: Peak RAM usage
- **CPU**: Core utilization
- **Accuracy**: If ground truth available
- **Concordance**: Agreement between pipelines

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
raptor compare --data sim_data/fastq/ --ground-truth sim_data/ground_truth.csv
```

## Understanding Results

The benchmark generates:
- `benchmark_results.json`: Raw metrics
- `comparison_report.html`: Interactive visualizations
- `pipeline_*/`: Individual pipeline outputs

## Comparison Metrics

**Accuracy Metrics** (with ground truth):
- Sensitivity (True Positive Rate)
- Specificity (True Negative Rate)  
- Precision, F1 Score
- Concordance with truth

**Performance Metrics**:
- Runtime (hours)
- Peak memory (GB)
- CPU efficiency
- Disk usage

## Advanced Options

```bash
# Quick benchmark (2 fastest pipelines)
raptor compare --data fastq/ --quick

# Full benchmark (all 8 pipelines)
raptor compare --data fastq/ --full

# Parallel execution
raptor compare --data fastq/ --parallel --max-jobs 4

# Custom config
raptor compare --data fastq/ --config my_config.yaml
```

## Interpreting Results

**Pipeline Rankings** based on:
1. Accuracy (if ground truth available)
2. Speed-accuracy tradeoff
3. Resource efficiency

**Use the results** to:
- Validate profile recommendations
- Understand trade-offs
- Choose optimal pipeline
- Publish comprehensive comparisons

## Python API
```python
from raptor import PipelineBenchmark

benchmark = PipelineBenchmark(
    data_dir='fastq/',
    output_dir='results/',
    threads=8
)

results = benchmark.run_pipelines([1, 3, 4])
benchmark.save_results(results)
```

## Tips

- Start with quick benchmark (2-3 pipelines)
- Use simulated data for testing
- Run overnight for large datasets
- Save intermediate results
- Compare to profile recommendation

See full benchmarking guide for detailed metrics and interpretation.

