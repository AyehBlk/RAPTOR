# Tutorial 2: Running Full Pipeline Benchmarking

**Level**: Intermediate  
**Time**: 3-6 hours (mostly automated)  
**Goal**: Compare all 8 pipelines on your data and choose the best one

---

## What You'll Learn

- How to run comprehensive pipeline benchmarking
- How to compare accuracy, speed, and resource usage
- How to interpret benchmark results
- How to make evidence-based pipeline decisions

---

## Prerequisites

- Completed [Tutorial 1: Getting Started](tutorial_01_getting_started.md)
- All bioinformatics tools installed
- FASTQ files or sufficient computational resources
- Time (benchmarking takes several hours)

---

## When to Use Full Benchmarking

### Use Benchmarking When:

✅ **Publishing results** - Need rigorous validation  
✅ **Method comparison** - Evaluating different approaches  
✅ **New data type** - Working with unusual experimental design  
✅ **High stakes** - Results inform important decisions  
✅ **Sufficient resources** - Have time and computational power  

### Skip Benchmarking When:

❌ **Exploratory analysis** - Use quick profiling instead  
❌ **Standard dataset** - Trust the recommendation  
❌ **Limited resources** - Benchmarking is resource-intensive  
❌ **Time pressure** - Profiling gives quick guidance  

---

## Step 1: Prepare Your Data

### Option A: Use Simulated Data (For Learning)

```bash
# Create benchmark directory
mkdir raptor_benchmark
cd raptor_benchmark

# Generate larger simulated dataset
raptor simulate --output sim_data/ --size medium

# This creates:
# - 5000 genes
# - 12 samples (6 per condition)
# - 1000 DE genes
# - Known ground truth
```

### Option B: Use Your Real FASTQ Files

```bash
# Organize your FASTQ files
data/
├── Sample1_R1.fastq.gz
├── Sample1_R2.fastq.gz
├── Sample2_R1.fastq.gz
├── Sample2_R2.fastq.gz
├── ...
└── metadata.csv
```

**Metadata format:**
```csv
sample,condition,replicate,fastq_r1,fastq_r2
Sample1,Control,1,Sample1_R1.fastq.gz,Sample1_R2.fastq.gz
Sample2,Control,2,Sample2_R1.fastq.gz,Sample2_R2.fastq.gz
Sample3,Treatment,1,Sample3_R1.fastq.gz,Sample3_R2.fastq.gz
Sample4,Treatment,2,Sample4_R1.fastq.gz,Sample4_R2.fastq.gz
```

---

## Step 2: Quick Benchmark (3 Fast Pipelines)

Before running all 8 pipelines, let's try a quick benchmark with the 3 fastest:

```bash
# Run quick benchmark (Pipelines 3, 4, 6)
raptor compare \
  --data sim_data/ \
  --output quick_benchmark/ \
  --pipelines 3 4 6 \
  --threads 8 \
  --mode quick
```

**Time**: ~30-60 minutes  
**Pipelines tested**:
- Pipeline 3: Salmon-edgeR (recommended for most)
- Pipeline 4: Kallisto-Sleuth (fastest)
- Pipeline 6: STAR-featureCounts-NOISeq (small sample friendly)

### Monitoring Progress

```bash
# In another terminal, watch progress
watch -n 5 'tail -n 20 quick_benchmark/benchmark.log'

# Or check current status
raptor status --benchmark quick_benchmark/
```

---

## Step 3: Review Quick Benchmark Results

```bash
# Generate comparison report
raptor report \
  --results quick_benchmark/ \
  --output quick_benchmark_report.html

# Open report
xdg-open quick_benchmark_report.html
```

### What to Look For:

**1. Runtime Comparison**
```
Pipeline | Runtime  | Speedup
---------|----------|--------
P4       | 15 min   | 4.0×
P3       | 22 min   | 2.7×
P6       | 60 min   | 1.0× (baseline)
```

**2. Accuracy Comparison (if ground truth available)**
```
Pipeline | Sensitivity | Precision | F1-Score
---------|-------------|-----------|----------
P3       | 0.89        | 0.87      | 0.88
P6       | 0.85        | 0.82      | 0.83
P4       | 0.82        | 0.84      | 0.83
```

**3. Resource Usage**
```
Pipeline | Peak Memory | Disk Space
---------|-------------|------------
P4       | 6 GB        | 15 GB
P3       | 12 GB       | 25 GB
P6       | 38 GB       | 45 GB
```

**Decision Point**: If one pipeline clearly dominates, you might not need to run all 8!

---

## Step 4: Full Benchmark (All 8 Pipelines)

If you need comprehensive comparison:

```bash
# Run complete benchmark
raptor compare \
  --data sim_data/ \
  --output full_benchmark/ \
  --pipelines all \
  --threads 16 \
  --memory 64G \
  --mode full \
  --save-intermediate
```

### Configuration Options:

```bash
# Customize resource allocation
raptor compare \
  --data sim_data/ \
  --output full_benchmark/ \
  --config custom_benchmark.yaml \
  --parallel 4  # Run 4 pipelines simultaneously
```

**Custom config example** (`custom_benchmark.yaml`):
```yaml
benchmarking:
  default_pipelines: [1, 2, 3, 4, 5, 6, 7, 8]
  mode: full
  metrics:
    runtime: true
    memory_usage: true
    accuracy: true
    concordance: true
  parallel_pipelines: true
  max_parallel: 4
  timeout_hours: 24
```

### Estimated Times:

| Dataset Size | Time (Sequential) | Time (Parallel) |
|--------------|-------------------|-----------------|
| Small (6 samples, 10M reads) | 3-4 hours | 1-2 hours |
| Medium (12 samples, 25M reads) | 8-12 hours | 3-4 hours |
| Large (24 samples, 50M reads) | 24-36 hours | 8-12 hours |

---

## Step 5: Monitor Long-Running Benchmark

### Real-time Monitoring:

```bash
# Check overall progress
raptor status --benchmark full_benchmark/

# Watch log file
tail -f full_benchmark/benchmark.log

# Check resource usage
htop  # or top
```

### What Each Pipeline Does:

```
✓ Pipeline 1 [=====>             ] 35% - STAR alignment
✓ Pipeline 2 [===========>       ] 75% - StringTie assembly  
✓ Pipeline 3 [=================] 100% - Complete
• Pipeline 4 [                  ] 0% - Queued
• Pipeline 5 [                  ] 0% - Queued
...
```

### Email Notifications (Optional):

```bash
# Get email when complete
raptor compare \
  --data sim_data/ \
  --output full_benchmark/ \
  --notify ayeh@example.com
```

---

## Step 6: Analyze Comprehensive Results

### Generate Comparison Report:

```bash
raptor report \
  --results full_benchmark/ \
  --output comprehensive_report.html \
  --include-all
```

The report includes:

### 1. Performance Summary Table

| Pipeline | Runtime | Memory | Accuracy (F1) | Rank |
|----------|---------|--------|---------------|------|
| P3: Salmon-edgeR | 22 min | 12 GB | 0.88 | 🥇 1 |
| P1: STAR-RSEM-DESeq2 | 3.5 h | 45 GB | 0.92 | 🥈 2 |
| P4: Kallisto-Sleuth | 15 min | 6 GB | 0.83 | 🥉 3 |
| P5: STAR-HTSeq-limma | 3.2 h | 42 GB | 0.89 | 4 |
| P2: HISAT2-StringTie | 2.1 h | 28 GB | 0.82 | 5 |
| P6: STAR-featureCounts-NOISeq | 3.0 h | 38 GB | 0.80 | 6 |
| P7: Bowtie2-RSEM-EBSeq | 2.8 h | 32 GB | 0.81 | 7 |
| P8: HISAT2-Cufflinks | 4.2 h | 35 GB | 0.75 | 8 |

### 2. Visualization Plots

**Speed vs. Accuracy Trade-off:**
```
Accuracy (F1)
    ↑
1.0 |     P1
0.9 |   P5   P3
0.8 |  P2 P6 P7
0.7 | P8     P4
    └───────────────→ Runtime (log scale)
      10m  1h   3h
```

**Resource Efficiency:**
```
Memory Usage (GB)
    ↑
50  | P1
40  | P5 P6
30  | P2 P7 P8
20  |
10  | P3
5   | P4
    └───────────────→ Accuracy (F1)
    0.7  0.8  0.9  1.0
```

### 3. Concordance Analysis

How well do pipelines agree?

```
Pipeline Pair        | Overlap | Jaccard
---------------------|---------|--------
P1 ↔ P3             | 87%     | 0.78
P1 ↔ P5             | 91%     | 0.84
P3 ↔ P4             | 82%     | 0.71
```

**High concordance** = Different methods agree  
**Low concordance** = Methods finding different genes

### 4. Detailed Performance Metrics

For each pipeline:
- **Sensitivity** (True Positive Rate)
- **Specificity** (True Negative Rate)
- **Precision** (Positive Predictive Value)
- **F1-Score** (Harmonic mean)
- **False Discovery Rate**
- **Runtime breakdown** by step
- **Memory profile** over time

---

## Step 7: Make Your Decision

### Decision Framework:

#### If Accuracy is Priority:
→ Choose **Pipeline 1** (STAR-RSEM-DESeq2)
- Highest accuracy
- Best for publications
- Worth the extra time

#### If Speed is Priority:
→ Choose **Pipeline 4** (Kallisto-Sleuth)
- Fastest option
- Good for exploratory analysis
- Sufficient for many uses

#### If Balance is Priority:
→ Choose **Pipeline 3** (Salmon-edgeR)
- Best speed/accuracy trade-off
- Recommended for most users
- Good resource efficiency

#### If Small Samples (n<6):
→ Choose **Pipeline 6** (STAR-featureCounts-NOISeq)
- Non-parametric statistics
- Better for limited replicates
- Handles high variation

---

## Step 8: Validate Your Choice

### Cross-validation:

```bash
# Run chosen pipeline on subset of data
raptor run \
  --pipeline 3 \
  --data sim_data/ \
  --output validation/ \
  --cross-validate 5  # 5-fold CV
```

### Robustness testing:

```bash
# Test with different parameters
raptor run \
  --pipeline 3 \
  --data sim_data/ \
  --output robust_test/ \
  --param-sweep "fdr:[0.01,0.05,0.1]" \
  --param-sweep "fc:[0.5,1.0,1.5]"
```

---

## Step 9: Document Your Analysis

### Generate Methods Section:

```bash
# Auto-generate methods text for your paper
raptor report \
  --results full_benchmark/ \
  --export-methods methods.txt
```

**Example output:**
```
RNA-seq differential expression analysis was performed using 
RAPTOR v2.0.0 (Bolouki, 2025). Eight different analysis pipelines 
were benchmarked on the dataset. Pipeline 3 (Salmon v1.10.0 with 
edgeR v3.40.0) was selected based on superior performance 
(F1-score: 0.88, runtime: 22 minutes). Raw reads were pseudo-aligned 
to the reference transcriptome using Salmon with --validateMappings 
and 30 bootstrap iterations. Gene-level quantification was performed 
using tximport. Differential expression was assessed using edgeR 
with TMM normalization and quasi-likelihood F-test. Genes with 
FDR < 0.05 and |log2FC| > 1 were considered differentially expressed.
```

### Save Configuration:

```bash
# Save exact configuration for reproducibility
raptor config export --output my_analysis_config.yaml
```

---

## Best Practices

### Do's:
✅ Profile data first before benchmarking  
✅ Use simulated data to test workflow  
✅ Save intermediate results  
✅ Document tool versions  
✅ Keep raw benchmark data  

### Don'ts:
❌ Run all 8 pipelines without quick test first  
❌ Benchmark on full dataset immediately  
❌ Ignore resource limitations  
❌ Delete intermediate files too early  
❌ Forget to document parameters  

---

## Troubleshooting

### Problem: Pipeline failed mid-run

```bash
# Check which pipeline failed
cat full_benchmark/benchmark.log | grep ERROR

# Resume from checkpoint
raptor compare \
  --data sim_data/ \
  --output full_benchmark/ \
  --resume
```

### Problem: Out of disk space

```bash
# Clean intermediate files
raptor clean --benchmark full_benchmark/ --keep-results

# Or run with minimal storage
raptor compare \
  --data sim_data/ \
  --output full_benchmark/ \
  --disk-efficient
```

### Problem: Taking too long

```bash
# Run subset of pipelines
raptor compare \
  --data sim_data/ \
  --output quick_test/ \
  --pipelines 1 3 4  # Just top 3

# Or use smaller data subset
raptor compare \
  --data sim_data/ \
  --output test/ \
  --subsample 5000  # Use 5000 genes only
```

---

## Advanced Topics

### Custom Pipeline Addition:

```bash
# Add your own pipeline to benchmark
raptor add-pipeline \
  --name "My_Custom_Pipeline" \
  --script my_pipeline.sh \
  --config my_pipeline_config.yaml
```

### Batch Processing:

```bash
# Benchmark multiple datasets
for dataset in dataset1 dataset2 dataset3; do
  raptor compare \
    --data $dataset/ \
    --output benchmark_$dataset/ \
    --pipelines all
done

# Compare across datasets
raptor meta-analysis \
  --benchmarks benchmark_*/\
  --output meta_comparison.html
```

---

## Summary

You've learned to:
- ✅ Run quick benchmarks for rapid comparison
- ✅ Execute comprehensive 8-pipeline benchmarking
- ✅ Monitor long-running analyses
- ✅ Interpret complex comparison reports
- ✅ Make evidence-based pipeline decisions
- ✅ Document your analysis for reproducibility

---

## Next Steps

1. **Try Tutorial 3**: [Understanding Pipeline Differences](tutorial_03_pipeline_details.md)
2. **Read**: [BENCHMARKING.md](../BENCHMARKING.md) for more details
3. **Explore**: [PIPELINES.md](../PIPELINES.md) for pipeline specifics
4. **Apply**: Use benchmarking on your real data!

---

**Tutorial by Ayeh Bolouki**  
University of Namur, Belgium  
For RAPTOR v2.0.0
