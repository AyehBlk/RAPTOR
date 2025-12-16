# Tutorial 03: Real Data Analysis with RAPTOR v2.1.1

**Level**: Intermediate  
**Time**: 2-4 hours  
**Goal**: Analyze your own RNA-seq data with quality assessment, ML guidance, and production-ready workflows

---

##  What You'll Learn

Real-world analysis techniques:
- ✅ Prepare and quality-check your own data
- ✅ Use ML recommendations for YOUR specific dataset
- ✅ Handle complex experimental designs
- ✅ Monitor and optimize resource usage
- ✅ Generate publication-ready results
- ✅ Troubleshoot common real-data issues
- ✅ Validate and interpret results confidently
- ✅ **NEW in v2.1.1**: Apply data-driven threshold optimization

---

## 🆕 What's New in v2.1.1

**Adaptive Threshold Optimizer (ATO)** for publication-ready analysis:

```bash
# After running your pipeline, optimize thresholds
raptor optimize-thresholds --results degs_all.csv --goal balanced --output ato_results/
```

ATO replaces arbitrary |logFC| > 1 with data-driven thresholds and generates **publication methods text** automatically. See [Tutorial 07](tutorial_07_threshold_optimizer.md) for details.

---

##  Prerequisites

**Required**:
- Completed [Tutorial 01](tutorial_01_getting_started.md) and [Tutorial 02](tutorial_02_benchmarking.md)
- RAPTOR v2.1.1 fully installed
- Your own RNA-seq data (count matrix or FASTQ files)
- Understanding of your experimental design

**Your Data Needs**:
- Count matrix (genes × samples) OR
- FASTQ files + reference genome
- Sample metadata (conditions, batches, etc.)
- Biological question in mind

---

##  Part 1: Data Preparation & Quality Assessment

### Step 1: Organize Your Data

Create a project directory:

```bash
# Replace with your project name
PROJECT="MyRNAseq_2025"
mkdir -p ~/$PROJECT/{raw_data,quality_reports,analysis_results}
cd ~/$PROJECT
```

### Step 2: Prepare Count Matrix

If you have a count matrix:

```bash
# Your file should look like:
# gene_id,sample1,sample2,sample3,...
# ENSG00000000001,245,267,198,...
# ENSG00000000002,89,102,95,...

# Copy your data
cp /path/to/your/counts.csv raw_data/
```

**Important formatting rules**:
- First column: Gene IDs (Ensembl, Symbol, or custom)
- Header row with sample names
- No missing values (use 0 for undetected genes)
- CSV or TSV format

### Step 3: Create Metadata File

**Critical for analysis!** Describe your samples:

```bash
cat > raw_data/metadata.txt << 'EOF'
sample_id,condition,batch,sequencing_date,rna_quality,notes
sample1,control,batch1,2025-01-15,RIN_8.2,good_quality
sample2,control,batch1,2025-01-15,RIN_8.5,good_quality
sample3,control,batch2,2025-01-22,RIN_8.1,good_quality
sample4,treated,batch1,2025-01-15,RIN_7.8,slightly_degraded
sample5,treated,batch1,2025-01-15,RIN_8.3,good_quality
sample6,treated,batch2,2025-01-22,RIN_8.0,good_quality
EOF
```

**Include all relevant factors**:
- Biological conditions (required)
- Batch information (if applicable)
- Technical factors (sequencing run, lane, date)
- Quality metrics (RIN score, library prep method)
- Any known issues or outliers

### Step 4: Initial Quality Assessment

**This is NEW in v2.1.1 and CRITICAL for real data!**

```bash
raptor assess-quality --counts raw_data/counts.csv \
                      --metadata raw_data/metadata.txt \
                      --output quality_reports/initial_assessment.html \
                      --detailed
```

**Open the report**:
```bash
firefox quality_reports/initial_assessment.html
# Or copy to your local machine and open
```

**What to look for**:

✅ **Pass** (proceed with analysis):
- Library size: > 1M reads per sample
- Gene detection: > 15,000 genes detected
- Sample correlation: > 0.85 within groups
- No major outliers in PCA
- Low zero inflation (< 70%)

⚠️ **Warning** (proceed with caution):
- Library size: 500K - 1M reads
- Gene detection: 10K - 15K genes
- Sample correlation: 0.70 - 0.85
- Potential batch effects detected
- Moderate zero inflation (70-85%)

❌ **Fail** (investigate before proceeding):
- Library size: < 500K reads
- Gene detection: < 10K genes
- Sample correlation: < 0.70
- Clear outliers or quality issues
- High zero inflation (> 85%)

### Step 5: Handle Quality Issues

#### Issue: Low Library Size

**Filter low-quality samples**:
```bash
raptor filter-samples --counts raw_data/counts.csv \
                      --metadata raw_data/metadata.txt \
                      --min-library-size 500000 \
                      --output raw_data/counts_filtered.csv
```

#### Issue: Batch Effects Detected

**Check severity**:
```bash
raptor detect-batch --counts raw_data/counts.csv \
                    --metadata raw_data/metadata.txt \
                    --output quality_reports/batch_analysis.html
```

If severe, plan to correct during analysis (see Part 3).

#### Issue: Outlier Samples

**Identify and decide**:
```bash
raptor detect-outliers --counts raw_data/counts.csv \
                       --metadata raw_data/metadata.txt \
                       --output quality_reports/outliers.csv
```

Review each outlier:
- Biological reason? (Keep, note in metadata)
- Technical issue? (Consider removing)
- Borderline? (Run analysis with and without)

---

##  Part 2: ML-Guided Pipeline Selection

### Step 1: Profile Your Data

```bash
raptor profile --counts raw_data/counts_filtered.csv \
               --metadata raw_data/metadata.txt \
               --detailed \
               --recommend \
               --output analysis_results/data_profile.json \
               --verbose
```

**Example output**:
```json
{
  "data_summary": {
    "n_genes": 58302,
    "n_samples": 12,
    "n_conditions": 2,
    "replicates_per_condition": 6,
    "experimental_design": "two_group_comparison"
  },
  "quality_metrics": {
    "library_size_mean": 15234567,
    "library_size_cv": 0.18,
    "gene_detection_rate": 0.72,
    "zero_inflation": 0.28,
    "dynamic_range": 126.5,
    "quality_score": 8.2
  },
  "complexity_assessment": {
    "biological_complexity": "moderate",
    "technical_noise": "low",
    "batch_effect_severity": "mild",
    "analysis_difficulty": "easy"
  },
  "ml_recommendation": {
    "top_pipeline": "STAR-DESeq2",
    "confidence": 0.89,
    "reasoning": [
      "High read count suitable for alignment-based methods",
      "Sufficient replicates (n=6) for DESeq2's dispersion estimation",
      "Low technical noise favors STAR's accuracy",
      "Mild batch effects manageable with DESeq2 design matrix"
    ],
    "alternatives": [
      {
        "pipeline": "Salmon-DESeq2",
        "score": 8.6,
        "note": "Faster, nearly equivalent accuracy"
      },
      {
        "pipeline": "HISAT2-DESeq2",
        "score": 8.3,
        "note": "Alternative aligner, similar results expected"
      }
    ],
    "not_recommended": [
      {
        "pipeline": "Kallisto-edgeR",
        "reason": "Your data favors DESeq2's methods"
      }
    ]
  },
  "expected_results": {
    "estimated_degs": "500-800 genes at FDR<0.05",
    "expected_sensitivity": 0.87,
    "expected_specificity": 0.91,
    "confidence_interval": "±5%"
  },
  "recommendations": {
    "primary_analysis": "Run STAR-DESeq2",
    "validation": "Compare with Salmon-DESeq2",
    "ensemble": "Consider ensemble if results differ",
    "batch_correction": "Include batch in design matrix"
  }
}
```

### Step 2: Understand the Recommendation

**Key questions to ask**:

1. **Does the confidence make sense?**
   - High (> 0.85): Strong evidence, proceed confidently
   - Medium (0.70-0.85): Good choice, consider comparing top 2
   - Low (< 0.70): Difficult dataset, compare top 3-5

2. **Do the reasons match your biology?**
   - Check if ML considers your specific tissue/organism
   - Verify it understands your experimental design
   - Ensure batch effect assessment is accurate

3. **Are alternatives close in score?**
   - Score difference < 0.5: Consider running both
   - Score difference > 1.0: Clear winner, use top choice

### Step 3: Decide on Strategy

**Strategy A: High Confidence (> 0.85)**
```bash
# Trust ML, run single pipeline
raptor run --pipeline STAR-DESeq2 \
           --counts raw_data/counts_filtered.csv \
           --metadata raw_data/metadata.txt \
           --output analysis_results/STAR_DESeq2/
```

**Strategy B: Medium Confidence (0.70-0.85)**
```bash
# Compare top 2
raptor compare --pipelines STAR-DESeq2,Salmon-DESeq2 \
               --counts raw_data/counts_filtered.csv \
               --metadata raw_data/metadata.txt \
               --output analysis_results/comparison/
```

**Strategy C: Low Confidence (< 0.70) or Important Project**
```bash
# Full ensemble analysis
raptor ensemble --counts raw_data/counts_filtered.csv \
                --metadata raw_data/metadata.txt \
                --pipelines top5 \
                --method weighted \
                --output analysis_results/ensemble/
```

---

##  Part 3: Running the Analysis

### Basic Two-Group Comparison

Most common case: Control vs. Treated

```bash
raptor run --pipeline STAR-DESeq2 \
           --counts raw_data/counts_filtered.csv \
           --metadata raw_data/metadata.txt \
           --condition-column condition \
           --control-level control \
           --treatment-level treated \
           --alpha 0.05 \
           --lfc-threshold 1.0 \
           --output analysis_results/main_analysis/ \
           --threads 4 \
           --monitor-resources
```

**Parameters explained**:
- `--alpha 0.05`: FDR cutoff (adjust for stringency)
- `--lfc-threshold 1.0`: Log2 fold change ≥ 1 (2-fold change)
- `--threads 4`: Use 4 CPU cores
- `--monitor-resources`: Track CPU/RAM usage

### Complex Design: Batch Correction

If batch effects are present:

```bash
raptor run --pipeline STAR-DESeq2 \
           --counts raw_data/counts_filtered.csv \
           --metadata raw_data/metadata.txt \
           --design "~ batch + condition" \
           --contrast condition_treated_vs_control \
           --output analysis_results/batch_corrected/ \
           --correct-batch \
           --batch-column batch
```

**The design formula**:
- `~ batch + condition`: Control for batch, test condition effect
- Order matters: covariates before variable of interest
- DESeq2 automatically handles this

### Multi-Factor Design

Example: Time × Treatment interaction

```bash
raptor run --pipeline STAR-DESeq2 \
           --counts raw_data/counts_filtered.csv \
           --metadata raw_data/metadata.txt \
           --design "~ timepoint + treatment + timepoint:treatment" \
           --output analysis_results/interaction_analysis/
```

Contrasts to test:
1. Main effect of treatment: `treatment_effect`
2. Main effect of time: `timepoint_effect`
3. Interaction: `timepoint:treatment`

### Paired Samples

Example: Before/After in same patients

```bash
raptor run --pipeline STAR-DESeq2 \
           --counts raw_data/counts_filtered.csv \
           --metadata raw_data/metadata.txt \
           --design "~ patient + condition" \
           --paired \
           --pair-column patient \
           --output analysis_results/paired_analysis/
```

### Monitor Analysis Progress

**Real-time monitoring (NEW in v2.1.1)**:

```bash
# In another terminal
raptor monitor --analysis analysis_results/main_analysis/ \
               --refresh 10
```

Or use the dashboard:
```bash
raptor dashboard --monitor analysis_results/main_analysis/
```

**What you'll see**:
- Current analysis step (alignment, counting, DE analysis)
- CPU and memory usage
- Estimated time remaining
- Warnings or errors

---

##  Part 4: Interpreting Results

### Directory Structure

```
analysis_results/main_analysis/
├── degs_all.csv                 # All genes with statistics
├── degs_significant.csv         # FDR < 0.05 & |LFC| > 1
├── normalized_counts.csv        # For further analysis
├── quality_metrics.json         # Analysis QC
├── plots/
│   ├── ma_plot.png             # M vs A plot
│   ├── volcano_plot.png        # -log10(p) vs log2FC
│   ├── pca_plot.png            # Sample clustering
│   ├── heatmap_top100.png      # Top DEGs heatmap
│   ├── dispersion_plot.png     # DESeq2 diagnostic
│   └── sample_distances.png    # Sample correlation
├── reports/
│   ├── analysis_report.html    # Comprehensive report
│   └── qc_report.html          # Quality control
└── logs/
    ├── raptor.log              # Analysis log
    └── resource_usage.log      # CPU/RAM tracking
```

### Check Analysis Quality

```bash
cat analysis_results/main_analysis/quality_metrics.json
```

**Look for**:
```json
{
  "dispersion_fit": "good",
  "outlier_genes": 23,
  "outlier_samples": 0,
  "model_convergence": "all_converged",
  "size_factors": [0.89, 1.12, 0.95, 1.05, 0.98, 1.01],
  "cook_distance_outliers": 2,
  "warnings": [],
  "overall_quality": "pass"
}
```

✅ **Good signs**:
- dispersion_fit: "good"
- outlier_samples: 0
- model_convergence: "all_converged"
- size_factors close to 1.0
- Few Cook's distance outliers

❌ **Red flags**:
- dispersion_fit: "poor"
- Many outlier samples
- Convergence issues
- Extreme size factors (< 0.5 or > 2.0)

### Examine DEGs

```bash
# Count significant genes
wc -l analysis_results/main_analysis/degs_significant.csv

# Preview top genes
head -20 analysis_results/main_analysis/degs_significant.csv
```

**Example output**:
```csv
gene_id,gene_name,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
ENSG00000130600,H19,25634.2,-5.67,0.23,-24.65,1.2e-134,2.1e-130
ENSG00000157764,BRAF,8932.4,4.82,0.19,25.37,3.4e-142,5.8e-138
ENSG00000171862,PTEN,12456.8,-3.94,0.17,-23.18,8.9e-119,1.5e-115
```

**Columns explained**:
- **baseMean**: Average expression across all samples
- **log2FoldChange**: Effect size (positive = upregulated)
- **lfcSE**: Standard error of LFC
- **stat**: Test statistic
- **pvalue**: Unadjusted p-value
- **padj**: FDR-adjusted p-value (Benjamini-Hochberg)

### Visualize Results

**MA Plot** - Overall view:
```bash
display analysis_results/main_analysis/plots/ma_plot.png
```
- X-axis: Average expression (log2)
- Y-axis: Log2 fold change
- **Red**: Significant DEGs
- **Blue**: Non-significant

**Volcano Plot** - Significance:
```bash
display analysis_results/main_analysis/plots/volcano_plot.png
```
- X-axis: Log2 fold change
- Y-axis: -log10(adjusted p-value)
- **Upper corners**: Most significant DEGs

**PCA Plot** - Sample clustering:
```bash
display analysis_results/main_analysis/plots/pca_plot.png
```
- Should see clustering by condition
- Check for batch effects (unwanted clustering)
- Identify potential outliers

### Compare to Expected Results

Recall ML prediction:
```bash
grep "estimated_degs" analysis_results/data_profile.json
# Predicted: 500-800 genes
```

Count actual DEGs:
```bash
wc -l analysis_results/main_analysis/degs_significant.csv
# If you got ~600 genes: ✅ matches prediction
# If you got ~50 genes: ⚠️ fewer than expected, investigate
# If you got ~5000 genes: ⚠️ more than expected, check thresholds
```

**If results differ significantly from predictions**:
1. Check data quality again
2. Verify experimental design is correct
3. Consider ensemble analysis
4. Review ML assumptions

### Step 5: Optimize Thresholds with ATO (NEW in v2.1.1)

After getting results, optimize your significance thresholds:

```bash
raptor optimize-thresholds --results analysis_results/main_analysis/degs_all.csv \
                           --goal balanced \
                           --output analysis_results/threshold_optimization/
```

**Why use ATO?**
- Replace arbitrary |logFC| > 1 with data-driven threshold
- Generate publication-ready methods text
- Justify your thresholds to reviewers

**View ATO results**:
```bash
cat analysis_results/threshold_optimization/ato_summary.json
```

```json
{
  "optimal_logfc_threshold": 0.847,
  "pi0_estimate": 0.782,
  "fdr_threshold": 0.05,
  "n_significant": 623,
  "analysis_goal": "balanced",
  "methods_text": "Differential expression significance..."
}
```

**Apply ATO threshold**:
```bash
raptor apply-threshold --results analysis_results/main_analysis/degs_all.csv \
                       --logfc-threshold 0.847 \
                       --output analysis_results/ato_filtered_degs.csv
```

---

##  Part 5: Advanced Analysis & Validation

### Ensemble Analysis for Robust Results

If:
- ML confidence was < 0.80
- Results seem unexpected
- This is for publication
- You want extra confidence

Run ensemble:

```bash
raptor ensemble --counts raw_data/counts_filtered.csv \
                --metadata raw_data/metadata.txt \
                --pipelines STAR-DESeq2,Salmon-DESeq2,HISAT2-DESeq2 \
                --method weighted \
                --output analysis_results/ensemble/ \
                --threads 8
```

**Compare results**:
```bash
raptor compare-results --single analysis_results/main_analysis/degs_significant.csv \
                       --ensemble analysis_results/ensemble/degs_consensus.csv \
                       --output analysis_results/comparison_report.html
```

**Interpretation**:
- **High overlap (> 80%)**: Single pipeline is reliable
- **Moderate overlap (60-80%)**: Ensemble adds value
- **Low overlap (< 60%)**: Challenging dataset, trust ensemble

### Functional Enrichment Analysis

Find biological meaning in your DEGs:

```bash
# Prepare gene list
cut -f1 analysis_results/main_analysis/degs_significant.csv > deglist.txt

# Option 1: Use built-in enrichment (if available)
raptor enrich --genes deglist.txt \
              --organism human \
              --output enrichment_results/

# Option 2: Export for external tools
raptor export-for-enrichment --degs analysis_results/main_analysis/degs_significant.csv \
                             --format gsea,david,ipa \
                             --output enrichment_input/
```

**Then use**:
- DAVID: https://david.ncifcrf.gov/
- GSEA: https://www.gsea-msigdb.org/
- Enrichr: https://maayanlab.cloud/Enrichr/
- clusterProfiler (R package)

### Gene-Level Validation

Check specific genes of interest:

```bash
raptor plot-gene --counts raw_data/counts_filtered.csv \
                 --metadata raw_data/metadata.txt \
                 --gene BRAF \
                 --output plots/BRAF_expression.png
```

Creates boxplot showing expression across samples.

### Cross-Validation

Test reproducibility by subsampling:

```bash
raptor cross-validate --counts raw_data/counts_filtered.csv \
                      --metadata raw_data/metadata.txt \
                      --pipeline STAR-DESeq2 \
                      --folds 5 \
                      --output cv_results/
```

**Good reproducibility**: > 85% of DEGs found across all folds

---

##  Part 6: Resource Optimization

### Monitor Resource Usage

Check what happened during analysis:

```bash
raptor resource-report --analysis analysis_results/main_analysis/ \
                       --output resource_analysis.html
```

**Example metrics**:
```
Pipeline: STAR-DESeq2
Total runtime: 2h 15m
Peak RAM: 8.2 GB
Average CPU: 340%
Disk I/O: 156 GB read, 42 GB written

Breakdown by step:
  - STAR alignment: 1h 45m, 7.8 GB RAM
  - featureCounts: 15m, 2.1 GB RAM
  - DESeq2 analysis: 15m, 3.5 GB RAM
```

### Optimize for Your System

**If you have limited RAM (< 8 GB)**:
```bash
raptor run --pipeline Salmon-DESeq2 \  # More memory efficient
           --low-memory \
           --threads 2
```

**If you want speed**:
```bash
raptor run --pipeline Kallisto-DESeq2 \  # Fastest
           --threads 8 \
           --quick-mode
```

**If you want accuracy** (publication):
```bash
raptor run --pipeline STAR-DESeq2 \  # Most accurate
           --threads 4 \
           --high-quality
```

### Batch Processing Multiple Datasets

For multiple experiments:

```bash
# Create batch configuration
cat > batch_config.yaml << 'EOF'
datasets:
  - name: experiment1
    counts: data/exp1_counts.csv
    metadata: data/exp1_metadata.txt
  - name: experiment2
    counts: data/exp2_counts.csv
    metadata: data/exp2_metadata.txt
  - name: experiment3
    counts: data/exp3_counts.csv
    metadata: data/exp3_metadata.txt

pipeline: STAR-DESeq2
threads: 4
output_dir: batch_results/
EOF

# Run batch analysis
raptor batch --config batch_config.yaml
```

---

##  Part 7: Publication-Ready Outputs

### Generate Comprehensive Report

```bash
raptor report --analysis analysis_results/main_analysis/ \
              --format publication \
              --include-methods \
              --include-qc \
              --output Final_Report.pdf
```

**Report sections**:
1. Executive Summary
2. Methods (ready for manuscript)
3. Quality Control
4. Results Summary
5. Key Findings
6. Statistical Details
7. All Figures
8. Supplementary Tables

### Create Figure Panel

```bash
raptor create-figure-panel --analysis analysis_results/main_analysis/ \
                           --layout "2x2" \
                           --panels "pca,volcano,ma,heatmap" \
                           --dpi 300 \
                           --output Figure1_DEG_Analysis.pdf
```

### Export Results Tables

**For supplementary materials**:
```bash
raptor export --analysis analysis_results/main_analysis/ \
              --format xlsx \
              --include-normalized-counts \
              --include-metadata \
              --output Supplementary_Data.xlsx
```

Creates Excel file with multiple sheets:
- Summary statistics
- All DEGs
- Significant DEGs
- Normalized counts
- Sample metadata
- Analysis parameters

### Methods Section

```bash
raptor generate-methods --analysis analysis_results/main_analysis/ \
                        --include-ato \
                        --output methods_section.txt
```

**Example output** (with ATO):
```
RNA-seq Differential Expression Analysis

Quality control and differential expression analysis were performed using 
RAPTOR v2.1.1 (Bolouki et al., 2025). Raw count data were quality-filtered 
using thresholds of minimum library size >500,000 reads and gene detection 
rate >70%. The STAR-DESeq2 pipeline was selected based on machine learning 
recommendations (confidence: 89%). 

Reads were aligned to the [species] genome ([version]) using STAR v2.7.10b 
with default parameters. Gene-level counts were quantified using featureCounts 
v2.0.3. Differential expression analysis was performed using DESeq2 v1.38.0.

Significance thresholds were determined using the Adaptive Threshold Optimizer 
(ATO) from RAPTOR v2.1.1. The proportion of true null hypotheses (π₀) was 
estimated at 0.782 using the Storey method. Based on the 'balanced' analysis 
goal, an optimal |log₂FC| threshold of 0.847 was determined. This data-driven 
threshold, combined with an FDR threshold of 0.050, yielded 623 significant 
differentially expressed genes.

[Analysis completed with n=X genes detected, n=Y significantly 
differentially expressed]
```

### ATO Methods Text (NEW in v2.1.1)

If you used threshold optimization, get the methods text separately:

```bash
cat analysis_results/threshold_optimization/methods_text.txt
# Copy directly to your paper's methods section
```

---

##  Troubleshooting Real Data Issues

### Issue 1: "No Significant DEGs Found"

**Possible causes**:
1. Small sample size (n<3 per group)
2. High biological variance
3. Small effect sizes
4. Inappropriate thresholds

**Solutions**:
```bash
# Lower thresholds
raptor run --pipeline STAR-DESeq2 \
           --alpha 0.10 \  # More lenient FDR
           --lfc-threshold 0.5 \  # Smaller fold change
           ...

# Check power analysis
raptor power-analysis --counts data.csv \
                      --metadata metadata.txt \
                      --output power_report.html
```

### Issue 2: "Too Many DEGs (>50% of genes)"

**Possible causes**:
1. Poor normalization
2. Major batch effects
3. Sample mislabeling
4. Wrong control group

**Solutions**:
```bash
# Verify sample labels
head metadata.txt

# Check normalization
raptor plot-normalization --counts data.csv --output norm_check.png

# Try batch correction
raptor run --design "~ batch + condition" --correct-batch ...
```

### Issue 3: "Samples Don't Cluster by Condition in PCA"

**Possible causes**:
1. Batch effects dominate
2. High biological variance
3. Mislabeled samples
4. Other confounders

**Solutions**:
```bash
# Investigate batch effects
raptor detect-batch --detailed ...

# Try regression of unwanted variation
raptor run --remove-unwanted-variation \
           --n-factors 2 \
           ...

# Check sample metadata carefully
```

### Issue 4: "Pipeline Failed: Out of Memory"

**Solutions**:
```bash
# Use memory-efficient pipeline
raptor run --pipeline Salmon-DESeq2 --low-memory

# Reduce threads
raptor run --threads 2

# Filter lowly expressed genes first
raptor filter-genes --counts data.csv \
                    --min-count 10 \
                    --min-samples 3 \
                    --output counts_filtered.csv
```

### Issue 5: "Results Differ from Previous Analysis"

**Check**:
1. RAPTOR version
2. Tool versions (STAR, DESeq2, etc.)
3. Filtering thresholds
4. Normalization method
5. Statistical thresholds

**Reproduce previous analysis**:
```bash
raptor run --pipeline STAR-DESeq2 \
           --version-lock previous_analysis.json \
           ...
```

---

##  Best Practices Summary

### Before Analysis
1. ✅ Quality-assess your data
2. ✅ Check for batch effects
3. ✅ Verify sample metadata
4. ✅ Get ML recommendations
5. ✅ Plan your analysis strategy

### During Analysis
1. ✅ Monitor resources
2. ✅ Check intermediate outputs
3. ✅ Save all parameters
4. ✅ Document decisions

### After Analysis
1. ✅ Validate key findings
2. ✅ Check QC metrics
3. ✅ Optimize thresholds with ATO (v2.1.1)
4. ✅ Compare to expectations
5. ✅ Perform functional enrichment
6. ✅ Prepare publication materials with methods text

---

##  What's Next?

**Tutorial 04: ML Deep Dive**
- Understand ML recommendations in depth
- Train custom models for your organism
- [→ tutorial_04_ml_recommendations.md](tutorial_04_ml_recommendations.md)

**Tutorial 05: Dashboard Mastery**
- Build custom workflows
- Share results with collaborators
- [→ tutorial_05_dashboard.md](tutorial_05_dashboard.md)

**Tutorial 06: Ensemble Analysis**
- Master advanced ensemble techniques
- Handle challenging datasets
- [→ tutorial_06_ensemble.md](tutorial_06_ensemble.md)

**Tutorial 07: Threshold Optimization** 🎯 (NEW)
- Master data-driven threshold selection
- Generate publication methods text
- [→ tutorial_07_threshold_optimizer.md](tutorial_07_threshold_optimizer.md)

---

##  Reference Documentation

- **[PROFILE_RECOMMEND.md](../PROFILE_RECOMMEND.md)**: Data profiling guide
- **[QUALITY_GUIDE.md](../QUALITY_GUIDE.md)**: Quality assessment details
- **[ML_GUIDE.md](../ML_GUIDE.md)**: ML recommendations explained
- **[THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md)**: ATO complete guide
- **[RESOURCE_MONITOR_GUIDE.md](../RESOURCE_MONITOR_GUIDE.md)**: Resource optimization
- **[TROUBLESHOOTING.md](../TROUBLESHOOTING.md)**: Common issues
- **[API.md](../API.md)**: Python API reference

---

##  Quick Command Reference

```bash
# Quality assessment
raptor assess-quality --counts data.csv --metadata meta.txt

# Get recommendations
raptor profile --counts data.csv --recommend

# Run analysis
raptor run --pipeline RECOMMENDED --counts data.csv --metadata meta.txt

# Optimize thresholds (NEW in v2.1.1)
raptor optimize-thresholds --results degs_all.csv --goal balanced

# Monitor progress
raptor dashboard --monitor results/

# Generate report (with ATO methods)
raptor report --analysis results/ --format publication --include-ato

# Ensemble analysis
raptor ensemble --counts data.csv --method weighted

# Export results
raptor export --analysis results/ --format xlsx
```

---

##  Need Help?

- **Documentation**: See [FAQ.md](../FAQ.md) and [TROUBLESHOOTING.md](../TROUBLESHOOTING.md)
- **GitHub Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com
- **Cite RAPTOR**: [Citation info in main README]

---

##  Tutorial Complete!

You should now be able to:
- [ ] Quality-assess real RNA-seq data
- [ ] Use ML recommendations effectively
- [ ] Handle complex experimental designs
- [ ] Monitor and optimize resource usage
- [ ] Interpret results confidently
- [ ] Optimize thresholds with ATO (v2.1.1)
- [ ] Generate publication-ready outputs with methods text
- [ ] Troubleshoot common issues
- [ ] Validate and ensemble results

**You're ready for real research! Go make discoveries! 🦖**

---

**Tutorial 03 - Real Data Analysis**  
*RAPTOR v2.1.1*  
Created by Ayeh Bolouki  
Last updated: December 2025
