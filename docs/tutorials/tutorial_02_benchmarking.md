# Tutorial 02: Advanced Benchmarking with RAPTOR v2.1.0

**Level**: Intermediate  
**Time**: 3-6 hours (mostly compute time)  
**Goal**: Systematically compare all pipelines with ML ranking and ensemble analysis

---

## 📋 What You'll Learn

Advanced benchmarking techniques:
- ✅ Compare all 8 RNA-seq pipelines on your data
- ✅ Use ML-based pipeline ranking
- ✅ Understand concordance analysis
- ✅ Apply ensemble methods for robust results
- ✅ Generate publication-ready comparison reports
- ✅ Monitor resource usage in real-time

---

##  Prerequisites

**Required**:
- Completed [Tutorial 01: Getting Started](tutorial_01_getting_started.md)
- RAPTOR v2.1.0 with all 8 pipelines installed
- At least 8 GB RAM and 4 CPU cores
- 3-6 hours for full benchmark

**Recommended**:
- Ground truth data or biological replicates
- 16 GB RAM for faster processing
- Understanding of statistical concepts (FDR, fold change)

---

##  Dataset for This Tutorial

We'll use a more complex dataset with known biology:

```bash
# Set up workspace
mkdir -p ~/raptor_benchmark
cd ~/raptor_benchmark

# Download realistic test data with spike-ins
# (Simulated based on actual RNA-seq characteristics)
cat > counts_benchmark.csv << 'EOF'
gene_id,ctrl_rep1,ctrl_rep2,ctrl_rep3,ctrl_rep4,treat_rep1,treat_rep2,treat_rep3,treat_rep4
ENSG00000001,2450,2380,2520,2410,145,152,138,149
ENSG00000002,3200,3180,3250,3210,3180,3240,3190,3220
ENSG00000003,890,920,875,905,2100,2050,2180,2120
ENSG00000004,125,135,118,128,122,138,115,130
ENSG00000005,4500,4450,4580,4520,180,175,185,178
ENSG00000006,670,685,655,680,675,690,660,685
ENSG00000007,55,62,48,58,1850,1820,1880,1840
ENSG00000008,1200,1180,1220,1190,1185,1205,1195,1210
ENSG00000009,2800,2750,2850,2790,95,88,102,92
ENSG00000010,420,435,408,428,425,440,415,432
EOF

# Create sample metadata
cat > metadata.txt << 'EOF'
sample,condition,batch,replicate
ctrl_rep1,control,batch1,1
ctrl_rep2,control,batch1,2
ctrl_rep3,control,batch2,3
ctrl_rep4,control,batch2,4
treat_rep1,treated,batch1,1
treat_rep2,treated,batch1,2
treat_rep3,treated,batch2,3
treat_rep4,treated,batch2,4
EOF
```

**Dataset characteristics**:
- 10 genes (simplified for tutorial)
- 8 samples (4 control + 4 treated)
- 4 replicates per condition (realistic)
- Known DE genes (genes 1, 3, 5, 7, 9)
- Batch effects (batch1 vs batch2)

---

##  Part 1: Full Benchmark - All Pipelines

### Step 1: Profile & Get ML Predictions

First, let RAPTOR analyze your data:

```bash
raptor profile --counts counts_benchmark.csv \
               --metadata metadata.txt \
               --recommend \
               --detailed \
               --output profile_detailed.json
```

**Output explanation**:
```json
{
  "data_characteristics": {
    "n_genes": 10,
    "n_samples": 8,
    "replicates_per_group": 4,
    "library_size_mean": 1847,
    "library_size_cv": 0.23,
    "zero_inflation": 0.0,
    "dynamic_range": 48.4,
    "batch_detected": true
  },
  "ml_predictions": {
    "top_pipeline": "STAR-DESeq2",
    "confidence": 0.91,
    "predicted_sensitivity": 0.88,
    "predicted_specificity": 0.92,
    "expected_degs": "3-5 genes"
  },
  "difficulty_assessment": "moderate"
}
```

### Step 2: Run Full Benchmark

Now compare all 8 pipelines:

```bash
raptor benchmark --counts counts_benchmark.csv \
                 --metadata metadata.txt \
                 --output benchmark_results/ \
                 --threads 4 \
                 --monitor-resources
```

**What happens during benchmark:**
1. Each pipeline runs independently (parallel execution)
2. Real-time resource monitoring
3. Quality assessment for each pipeline
4. Concordance analysis across all pipelines
5. ML-based ranking
6. Generate comparison report

**Expected runtime**: 2-4 hours (depends on hardware)

**Monitor progress** in real-time:
```bash
# In another terminal
tail -f benchmark_results/benchmark.log
```

Or use the dashboard:
```bash
raptor dashboard --monitor benchmark_results/
```

---

##  Part 2: Analyzing Benchmark Results

### Directory Structure

After completion:

```
benchmark_results/
├── summary_report.html           # Main comparison report
├── ml_ranking.json              # ML-based pipeline scores
├── concordance_matrix.png       # How similar are pipelines?
├── performance_comparison.png   # Side-by-side metrics
├── ensemble_results/            # Combined analysis
│   ├── degs_consensus.csv      # Genes found by multiple pipelines
│   └── ensemble_ranking.csv    # Weighted gene rankings
└── individual_pipelines/        # Results from each pipeline
    ├── STAR_DESeq2/
    ├── Salmon_DESeq2/
    ├── Kallisto_DESeq2/
    ├── HISAT2_DESeq2/
    ├── STAR_edgeR/
    ├── Salmon_edgeR/
    ├── Kallisto_edgeR/
    └── HISAT2_edgeR/
```

### View Summary Report

```bash
# Open in browser
firefox benchmark_results/summary_report.html
# Or
google-chrome benchmark_results/summary_report.html
```

**Report sections**:
1. **Executive Summary**: Key findings, ML recommendations
2. **Pipeline Performance**: Detailed metrics for each pipeline
3. **Concordance Analysis**: How well pipelines agree
4. **Ensemble Results**: Combined robust findings
5. **Resource Usage**: Computational requirements
6. **Recommendations**: Which pipeline(s) to trust

### Check ML Rankings

```bash
cat benchmark_results/ml_ranking.json
```

Example output:
```json
{
  "rankings": [
    {
      "rank": 1,
      "pipeline": "STAR-DESeq2",
      "overall_score": 8.9,
      "sensitivity_score": 8.8,
      "specificity_score": 9.2,
      "concordance_score": 9.1,
      "resource_efficiency": 8.5,
      "recommendation": "Best overall choice"
    },
    {
      "rank": 2,
      "pipeline": "Salmon-DESeq2",
      "overall_score": 8.7,
      "sensitivity_score": 8.6,
      "specificity_score": 9.0,
      "concordance_score": 8.9,
      "resource_efficiency": 9.2,
      "recommendation": "Fast alternative with similar accuracy"
    },
    {
      "rank": 3,
      "pipeline": "Kallisto-DESeq2",
      "overall_score": 8.5,
      "sensitivity_score": 8.4,
      "specificity_score": 8.8,
      "concordance_score": 8.7,
      "resource_efficiency": 9.5,
      "recommendation": "Fastest option, slightly lower sensitivity"
    }
  ],
  "consensus_recommendation": "STAR-DESeq2 or Salmon-DESeq2",
  "confidence": 0.91
}
```

**Understanding scores** (out of 10):
- **Sensitivity**: Ability to detect true positives
- **Specificity**: Ability to avoid false positives
- **Concordance**: Agreement with other pipelines
- **Resource efficiency**: Speed and memory usage

---

##  Part 3: Concordance Analysis

### What is Concordance?

Concordance measures how well different pipelines agree on:
- Which genes are differentially expressed
- Direction of change (up vs down)
- Magnitude of change

### View Concordance Matrix

```bash
# Image view
display benchmark_results/concordance_matrix.png

# Or numerical data
cat benchmark_results/concordance_data.csv
```

Example concordance matrix:

```
Pipeline         STAR-D  Salmon-D  Kall-D  HISAT-D  STAR-E  Salmon-E  Kall-E  HISAT-E
STAR-DESeq2      1.00    0.92      0.88    0.85     0.78    0.75      0.72    0.70
Salmon-DESeq2    0.92    1.00      0.90    0.87     0.76    0.82      0.79    0.73
Kallisto-DESeq2  0.88    0.90      1.00    0.89     0.74    0.80      0.85    0.75
HISAT2-DESeq2    0.85    0.87      0.89    1.00     0.72    0.78      0.82    0.88
STAR-edgeR       0.78    0.76      0.74    0.72     1.00    0.88      0.85    0.80
Salmon-edgeR     0.75    0.82      0.80    0.78     0.88    1.00      0.90    0.83
Kallisto-edgeR   0.72    0.79      0.85    0.82     0.85    0.90      1.00    0.86
HISAT2-edgeR     0.70    0.73      0.75    0.88     0.80    0.83      0.86    1.00
```

**Interpretation**:
- **> 0.90**: Excellent agreement
- **0.80-0.90**: Good agreement
- **0.70-0.80**: Moderate agreement
- **< 0.70**: Poor agreement (investigate why!)

**Key insights**:
1. DESeq2 pipelines cluster together (0.85-0.92)
2. edgeR pipelines cluster together (0.80-0.90)
3. Cross-method concordance is lower (0.70-0.82)
4. Alignment method matters less than DE method

### Analyze Discordant Genes

Find genes with disagreement:

```bash
raptor analyze-concordance --benchmark benchmark_results/ \
                           --threshold 0.5 \
                           --output discordant_genes.csv
```

Example output:
```csv
gene_id,n_pipelines_sig,concordance,reason
ENSG00000011,4,0.50,Borderline significance
ENSG00000012,2,0.25,Low expression + high variance
ENSG00000013,1,0.13,Outlier in one sample
```

**Common reasons for discordance**:
- Borderline p-values (near 0.05)
- Low gene expression
- High biological variance
- Outlier samples
- Batch effects

---

##  Part 4: Ensemble Analysis - Combining Pipelines

### Why Use Ensemble?

Instead of trusting one pipeline, combine multiple pipelines for:
- **More robust results**
- **Reduced false positives**
- **Better confidence estimates**
- **Protection against pipeline-specific biases**

### Method 1: Consensus (Simple)

Find genes detected by at least N pipelines:

```bash
raptor ensemble --benchmark benchmark_results/ \
                --method consensus \
                --min-pipelines 5 \
                --output ensemble_consensus/
```

Output:
```csv
gene_id,n_detections,pipelines,consensus_padj,consensus_lfc
ENSG00000001,8,all,-4.71,0.0001
ENSG00000003,8,all,4.68,0.0001
ENSG00000005,7,all_except_Kallisto-edgeR,-4.55,0.0002
ENSG00000007,6,DESeq2_all+STAR-edgeR,4.42,0.0008
ENSG00000009,5,DESeq2_all+Salmon-edgeR,-4.28,0.0012
```

**Interpretation**:
- Genes found by ≥5 pipelines are highly reliable
- Fewer pipelines = lower confidence
- Use this for high-confidence gene lists

### Method 2: Weighted Voting

Weight pipelines by their ML scores:

```bash
raptor ensemble --benchmark benchmark_results/ \
                --method weighted \
                --weights ml_ranking.json \
                --output ensemble_weighted/
```

This gives each pipeline a vote proportional to its ML score:
- STAR-DESeq2: weight = 8.9
- Salmon-DESeq2: weight = 8.7
- Kallisto-DESeq2: weight = 8.5
- etc.

Output includes weighted p-values and fold changes:
```csv
gene_id,weighted_padj,weighted_lfc,confidence_score
ENSG00000001,0.00008,-4.69,0.95
ENSG00000003,0.00009,4.66,0.94
ENSG00000005,0.00015,-4.52,0.91
```

### Method 3: RankProduct

Combine rankings across pipelines (robust to outliers):

```bash
raptor ensemble --benchmark benchmark_results/ \
                --method rankprod \
                --output ensemble_rankprod/
```

Best for:
- Data with outliers
- When you don't trust ML scores
- Publication-ready results (widely accepted method)

### Compare Ensemble Methods

```bash
raptor compare-ensemble --benchmark benchmark_results/ \
                        --methods consensus,weighted,rankprod \
                        --output ensemble_comparison.html
```

**Which method to choose?**

| Method | Best For | Pros | Cons |
|--------|----------|------|------|
| **Consensus** | High confidence list | Simple, conservative | May miss true positives |
| **Weighted** | Leveraging ML insights | Uses pipeline quality | Depends on ML accuracy |
| **RankProduct** | Publication | Widely accepted, robust | Computationally intensive |

**Recommendation**: Use weighted for exploration, consensus for validation, rankprod for publications.

---

##  Part 5: Resource Monitoring & Optimization

### Real-Time Monitoring

During benchmark, track resource usage:

```bash
# View live dashboard
raptor dashboard --monitor benchmark_results/

# Or check resource log
tail -f benchmark_results/resource_usage.log
```

Example resource log:
```
Time     Pipeline          CPU%   RAM_GB  Status
14:23:01 STAR-DESeq2       245%   3.2     Aligning reads
14:23:31 STAR-DESeq2       180%   4.1     Counting
14:24:01 STAR-DESeq2       98%    2.8     DE analysis
14:24:31 STAR-DESeq2       Complete        Success
14:23:05 Salmon-DESeq2     156%   1.8     Quasi-mapping
```

### Performance Comparison

```bash
cat benchmark_results/performance_summary.csv
```

```csv
Pipeline,Runtime_min,Peak_RAM_GB,Accuracy,Recommendation
STAR-DESeq2,45,4.2,High,Best for accuracy
Salmon-DESeq2,12,2.1,High,Best balance
Kallisto-DESeq2,8,1.5,Medium-High,Best for speed
HISAT2-DESeq2,38,3.8,High,Alternative to STAR
STAR-edgeR,43,4.0,Medium,Use DESeq2 instead
Salmon-edgeR,11,2.0,Medium,Use DESeq2 instead
Kallisto-edgeR,7,1.4,Medium,Use DESeq2 instead
HISAT2-edgeR,36,3.7,Medium,Use DESeq2 instead
```

### Optimization Tips

**For limited resources** (< 8 GB RAM):
```bash
raptor benchmark --pipelines Kallisto-DESeq2,Salmon-DESeq2 \
                 --threads 2 \
                 --low-memory
```

**For speed** (need results fast):
```bash
raptor benchmark --pipelines Kallisto-DESeq2,Salmon-DESeq2,Kallisto-edgeR \
                 --threads 8
```

**For accuracy** (publication):
```bash
raptor benchmark --pipelines STAR-DESeq2,Salmon-DESeq2,HISAT2-DESeq2 \
                 --threads 4 \
                 --ensemble weighted
```

---

##  Part 6: Validation & Quality Control

### Cross-Validation with Biological Replicates

If you have replicates, test reproducibility:

```bash
raptor cross-validate --counts counts_benchmark.csv \
                      --metadata metadata.txt \
                      --pipeline STAR-DESeq2 \
                      --folds 4 \
                      --output cv_results/
```

This:
1. Splits replicates into training/test sets
2. Runs pipeline on training data
3. Validates on test data
4. Repeats for all combinations
5. Reports reproducibility

### Spike-In Validation (if available)

```bash
raptor validate-spikein --benchmark benchmark_results/ \
                        --spikein-list known_spikes.txt \
                        --output validation_report.html
```

### Compare to Ground Truth

If you have qPCR or other validation:

```bash
# Create ground truth file
cat > ground_truth.csv << 'EOF'
gene_id,true_de,true_lfc
ENSG00000001,TRUE,-4.5
ENSG00000003,TRUE,4.6
ENSG00000005,TRUE,-4.4
ENSG00000007,TRUE,4.3
ENSG00000009,TRUE,-4.2
EOF

raptor validate --benchmark benchmark_results/ \
                --ground-truth ground_truth.csv \
                --output validation_metrics.json
```

Output:
```json
{
  "STAR-DESeq2": {
    "sensitivity": 0.88,
    "specificity": 0.92,
    "accuracy": 0.90,
    "mcc": 0.80
  },
  "Salmon-DESeq2": {
    "sensitivity": 0.86,
    "specificity": 0.91,
    "accuracy": 0.89,
    "mcc": 0.78
  }
}
```

**Metrics explained**:
- **Sensitivity (Recall)**: % of true DE genes detected
- **Specificity**: % of non-DE genes correctly identified
- **Accuracy**: Overall correctness
- **MCC** (Matthews Correlation Coefficient): Balanced measure (-1 to 1, higher better)

---

##  Part 7: Publication-Ready Reports

### Generate Comprehensive Report

```bash
raptor report --benchmark benchmark_results/ \
              --format publication \
              --include-methods \
              --output RAPTOR_Analysis_Report.pdf
```

**Report includes**:
- Methods section (ready for paper)
- Summary statistics table
- Comparison figures
- Ensemble results
- Quality metrics
- Computational requirements

### Create Supplementary Figures

```bash
raptor plot-suite --benchmark benchmark_results/ \
                  --output figures/ \
                  --dpi 300 \
                  --format png,pdf
```

Generates:
- `fig1_concordance_heatmap.pdf`
- `fig2_venn_diagram.pdf`
- `fig3_performance_comparison.pdf`
- `fig4_ma_plots.pdf` (all pipelines)
- `fig5_volcano_plots.pdf` (all pipelines)

### Export Data Tables

```bash
# All in one spreadsheet
raptor export --benchmark benchmark_results/ \
              --format xlsx \
              --output Supplementary_Tables.xlsx
```

Creates spreadsheet with tabs:
- Summary Statistics
- All DEGs (all pipelines)
- Consensus DEGs
- Ensemble Results
- Performance Metrics
- Concordance Matrix

---

##  Advanced Topics

### Custom Pipeline Comparison

Compare specific pipelines only:

```bash
raptor compare --pipelines STAR-DESeq2,Salmon-DESeq2,HISAT2-DESeq2 \
               --counts counts_benchmark.csv \
               --metadata metadata.txt \
               --output custom_comparison/
```

### Parameter Sensitivity Analysis

Test how parameter choices affect results:

```bash
raptor sensitivity --pipeline STAR-DESeq2 \
                   --counts counts_benchmark.csv \
                   --parameter alpha \
                   --values 0.01,0.05,0.10 \
                   --output sensitivity_alpha/
```

### Batch Effect Comparison

See how pipelines handle batch effects:

```bash
raptor benchmark --counts counts_benchmark.csv \
                 --metadata metadata.txt \
                 --correct-batch \
                 --compare-with-without \
                 --output batch_comparison/
```

### Time-Series or Multi-Factor

For complex designs:

```bash
raptor benchmark --counts timeseries_counts.csv \
                 --design "~ time + treatment + time:treatment" \
                 --output complex_design/
```

---

##  Best Practices & Tips

### 1. Start with ML Recommendations
```bash
# Always profile first
raptor profile --counts yourdata.csv --recommend --detailed
```

### 2. Run Full Benchmark for Important Projects
```bash
# For publications, grants, important decisions
raptor benchmark --counts yourdata.csv --threads 8
```

### 3. Use Ensemble for Robust Results
```bash
# Combine top 3-5 pipelines
raptor ensemble --method weighted --min-pipelines 3
```

### 4. Monitor Resources
```bash
# Especially for large datasets
raptor benchmark --monitor-resources --log-interval 30
```

### 5. Validate Results
```bash
# Use cross-validation or ground truth if available
raptor validate --ground-truth known_genes.csv
```

### 6. Document Everything
```bash
# Save commands and parameters
raptor benchmark ... > analysis_log.txt 2>&1
```

---

##  Troubleshooting

### "Benchmark is taking too long"

**Solutions**:
```bash
# Option 1: Reduce pipelines
raptor benchmark --pipelines STAR-DESeq2,Salmon-DESeq2,Kallisto-DESeq2

# Option 2: More threads
raptor benchmark --threads 8

# Option 3: Lower resolution (faster but less accurate)
raptor benchmark --quick-mode
```

### "Some pipelines failed"

Check individual logs:
```bash
ls benchmark_results/individual_pipelines/*/error.log
cat benchmark_results/individual_pipelines/STAR_DESeq2/error.log
```

Common issues:
- Insufficient memory → Use `--low-memory` flag
- Missing tools → Check installation with `raptor check-install`
- Incompatible data → Run quality assessment first

### "Concordance is very low (< 0.70)"

**Investigate**:
```bash
raptor analyze-concordance --benchmark benchmark_results/ \
                           --detailed \
                           --output concordance_investigation/
```

**Common causes**:
- Poor data quality
- Batch effects
- Too few replicates
- Inappropriate pipeline for data type

### "ML recommendations don't match my expectations"

```bash
# Get detailed explanation
raptor ml-explain --profile profile.json --verbose

# Override with manual selection
raptor benchmark --pipelines YOUR_CHOICE
```

---

##  What's Next?

### Continue Learning:

**Tutorial 03: Real Data Analysis**
- Apply benchmarking to your actual datasets
- Handle complex experimental designs
- Optimize for your specific research
- [→ tutorial_03_real_data.md](tutorial_03_real_data.md)

**Tutorial 04: ML Deep Dive**
- Understand ML recommendation system
- Train custom models
- Fine-tune for your organism/tissue
- [→ tutorial_04_ml_recommendations.md](tutorial_04_ml_recommendations.md)

**Tutorial 06: Ensemble Methods**
- Master advanced ensemble techniques
- Develop custom combination strategies
- Handle edge cases
- [→ tutorial_06_ensemble.md](tutorial_06_ensemble.md)

---

##  Reference Documentation

- **[BENCHMARKING.md](../BENCHMARKING.md)**: Complete benchmarking guide
- **[ENSEMBLE_GUIDE.md](../ENSEMBLE_GUIDE.md)**: Ensemble analysis details
- **[ML_GUIDE.md](../ML_GUIDE.md)**: Machine learning explanations
- **[RESOURCE_MONITOR_GUIDE.md](../RESOURCE_MONITOR_GUIDE.md)**: Resource monitoring
- **[API.md](../API.md)**: Python API for custom analysis

---

##  Quick Reference

```bash
# Full benchmark
raptor benchmark --counts data.csv --metadata meta.txt --threads 4

# Quick comparison (top 3)
raptor compare --pipelines top3 --counts data.csv

# Ensemble analysis
raptor ensemble --benchmark results/ --method weighted

# View ML rankings
cat results/ml_ranking.json

# Generate report
raptor report --benchmark results/ --format publication

# Monitor resources
raptor dashboard --monitor results/

# Validate results
raptor validate --benchmark results/ --ground-truth truth.csv

# Export for publication
raptor export --benchmark results/ --format xlsx
```

---

##  Need Help?

- **Questions**: [FAQ.md](../FAQ.md)
- **Issues**: [TROUBLESHOOTING.md](../TROUBLESHOOTING.md)
- **GitHub**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

##  Tutorial Checklist

Did you:
- [ ] Run full benchmark on test data
- [ ] Understand concordance analysis
- [ ] Try ensemble methods (consensus, weighted, rankprod)
- [ ] Monitor resource usage
- [ ] Generate comparison reports
- [ ] Interpret ML rankings
- [ ] Know how to validate results
- [ ] Create publication-ready figures

**Excellent work! You're now a RAPTOR power user! 🦖**

---

**Tutorial 02 - Advanced Benchmarking**  
*RAPTOR v2.1.0*  
Created by Ayeh Bolouki  
Last updated: November 2025
