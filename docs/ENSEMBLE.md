# 🎭 RAPTOR v2.1.1 Ensemble Analysis Guide

**Combine Multiple Pipelines for Robust Results**

Ensemble analysis combines results from multiple RNA-seq pipelines to generate more robust and reliable differential expression results.

---

## 🆕 What's New in v2.1.1

### 🎯 Adaptive Threshold Optimizer (ATO) Integration

- **Uniform Thresholds** - Apply data-driven thresholds across all pipelines
- **Ensemble-Optimized Cutoffs** - Thresholds optimized for combined results
- **Publication Methods Text** - Auto-generated for ensemble analyses

```python
from raptor.threshold_optimizer import optimize_thresholds

# After ensemble combination, optimize thresholds
result = optimize_thresholds(ensemble_results, goal='balanced')
print(f"Uniform threshold for all pipelines: |logFC| > {result.logfc_threshold:.2f}")
```

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Why Use Ensemble Analysis](#why-use-ensemble-analysis)
3. [Quick Start](#quick-start)
4. [Ensemble Methods](#ensemble-methods)
5. [🎯 ATO Integration](#ato-integration) - **NEW!**
6. [Running Ensemble Analysis](#running-ensemble-analysis)
7. [Interpreting Results](#interpreting-results)
8. [Advanced Options](#advanced-options)
9. [Best Practices](#best-practices)
10. [Case Studies](#case-studies)
11. [Troubleshooting](#troubleshooting)

---

## 🌟 Overview

### What is Ensemble Analysis?

Ensemble analysis runs multiple RNA-seq pipelines on the same data and combines their results to produce a **consensus gene list** that is more reliable than any single pipeline.

```
Pipeline 1 (STAR-RSEM-DESeq2)  ]
Pipeline 3 (Salmon-edgeR)       ] → Ensemble → 🎯 ATO → Robust
Pipeline 4 (Kallisto-DESeq2)    ]   Analysis   Thresholds  Results
```

### Key Benefits

✅ **More Robust** - Reduces pipeline-specific biases  
✅ **Higher Confidence** - Genes found by multiple methods  
✅ **Publication-Ready** - Demonstrates rigorous validation  
✅ **🎯 Optimized Thresholds** - Data-driven cutoffs (NEW in v2.1.1)  
✅ **Catches Errors** - Outlier detection  

---

## 🤔 Why Use Ensemble Analysis

### When to Use Ensemble

✅ **High-Stakes Decisions** - Clinical applications, drug targets  
✅ **Publication Requirements** - Reviewers request validation  
✅ **Novel Experiments** - Unusual data, new conditions  
✅ **Conflicting Results** - Different pipelines disagree  
✅ **Need Defensible Thresholds** - Use ATO on combined results (NEW)

### When Single Pipeline is Fine

❌ **Exploratory Analysis** - Initial data examination  
❌ **Time Constraints** - Quick turnaround needed  
❌ **Limited Resources** - Computational constraints  
❌ **Standard Data** - Well-characterized experiments  

---

## ⚡ Quick Start

### Basic Ensemble with ATO (NEW Recommended Workflow)

```bash
# 1. Run 3-pipeline ensemble
raptor ensemble \
  --data data/fastq/ \
  --pipelines 1,3,4 \
  --output ensemble_results/

# 2. Optimize thresholds with ATO (NEW!)
python -c "
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd
df = pd.read_csv('ensemble_results/consensus_genes.csv')
result = optimize_thresholds(df, goal='balanced')
print(f'Uniform threshold: |logFC| > {result.logfc_threshold:.2f}')
result.results_df.to_csv('ensemble_results/optimized_consensus.csv')
"
```

### Full Ensemble

```bash
# Run all 8 pipelines (thorough)
raptor ensemble \
  --data data/fastq/ \
  --pipelines all \
  --use-ato \
  --ato-goal balanced \
  --output full_ensemble/
```

### From Existing Results

```bash
# Combine pre-computed results
raptor ensemble combine \
  --results pipeline1_results/,pipeline3_results/,pipeline4_results/ \
  --use-ato \
  --output ensemble_combined/
```

---

## 🔧 Ensemble Methods

### 1. Weighted Average (Default)

Each pipeline gets a weight based on its performance. Results are averaged proportionally.

```python
weights = {
    'Pipeline 1': 0.40,  # Highest accuracy
    'Pipeline 3': 0.35,  # Fast and accurate
    'Pipeline 4': 0.25   # Ultra-fast
}
```

**Best for:** Most situations, balances all inputs

### 2. Rank Product

Combines ranks rather than raw values. Less sensitive to outliers.

**Best for:** Datasets with very different scales

### 3. Vote Counting

A gene is DE if majority of pipelines agree.

**Best for:** Conservative analysis, high specificity

### 4. Robust Rank Aggregation (RRA)

Statistical method to combine ranked lists.

**Best for:** Rigorous statistical combination

### 5. Meta-Analysis

Statistical meta-analysis combining effect sizes and p-values.

**Best for:** Publication-quality combined statistics

---

## 🎯 ATO Integration (NEW in v2.1.1)

### Why Use ATO with Ensemble?

**Problem:** Different pipelines may use different thresholds
```
Pipeline 1: |logFC| > 1.0, FDR < 0.05 → 1,234 genes
Pipeline 3: |logFC| > 0.58, FDR < 0.05 → 1,567 genes  
Pipeline 4: |logFC| > 1.0, FDR < 0.01 → 892 genes
```

**Solution:** ATO provides uniform, data-driven thresholds
```python
from raptor.threshold_optimizer import optimize_thresholds

# Apply to combined results
result = optimize_thresholds(ensemble_combined, goal='balanced')
# Now all pipelines use: |logFC| > 0.73, FDR < 0.05
```

### Configuration with ATO

```yaml
# ensemble_config.yaml (v2.1.1)
ensemble:
  pipelines:
    - 1  # STAR-RSEM-DESeq2
    - 3  # Salmon-edgeR
    - 4  # Kallisto-DESeq2
  
  method: weighted_average
  min_pipelines: 2

  # NEW: ATO settings
  threshold_optimizer:
    enabled: true
    goal: "balanced"           # discovery, balanced, validation
    uniform_thresholds: true   # Apply same threshold to all pipelines
    generate_methods_text: true
```

### ATO Workflow Options

#### Option 1: Optimize After Combination (Recommended)

```python
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

# 1. Run ensemble with default thresholds
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(pipeline_results, method='weighted_average')

# 2. Apply ATO to combined results
result = optimize_thresholds(
    ensemble['combined_de'],
    logfc_col='log2FC_consensus',
    pvalue_col='FDR_consensus',
    goal='balanced'
)

# 3. Use optimized thresholds
significant = result.results_df[result.results_df['significant']]
print(f"High-confidence genes: {len(significant)}")
print(f"\nMethods:\n{result.methods_text}")
```

#### Option 2: Uniform Thresholds Across Pipelines

```python
from raptor.threshold_optimizer import optimize_thresholds

# First, determine optimal threshold from combined data
combined_logfc = pd.concat([p1['logFC'], p3['logFC'], p4['logFC']])
combined_pval = pd.concat([p1['pvalue'], p3['pvalue'], p4['pvalue']])
combined = pd.DataFrame({'logFC': combined_logfc, 'pvalue': combined_pval})

result = optimize_thresholds(combined, goal='balanced')
threshold = result.logfc_threshold

# Apply uniform threshold to all pipelines
for pipeline_df in [p1, p3, p4]:
    pipeline_df['significant'] = (
        (abs(pipeline_df['logFC']) > threshold) & 
        (pipeline_df['FDR'] < 0.05)
    )
```

#### Option 3: Integrated Ensemble + ATO

```bash
# Single command with ATO
raptor ensemble \
  --data data/ \
  --pipelines 1,3,4 \
  --use-ato \
  --ato-goal balanced \
  --output results/
```

### ATO Methods Text for Ensemble

```python
result = optimize_thresholds(ensemble_combined, goal='balanced')
print(result.methods_text)
```

**Example Output:**
```
"Ensemble differential expression analysis was performed using 
three complementary pipelines (STAR-RSEM-DESeq2, Salmon-edgeR, 
Kallisto-DESeq2) with results combined via weighted averaging. 
Significance thresholds were determined using the Adaptive 
Threshold Optimizer (ATO) with the 'balanced' analysis goal. 
The proportion of true null hypotheses (π₀) was estimated at 
0.83 using Storey's spline method. A uniform adjusted p-value 
threshold of 0.05 (Benjamini-Hochberg FDR correction) and log₂ 
fold change threshold of 0.71 were applied across all pipelines, 
identifying 1,156 consensus differentially expressed genes 
(agreement score ≥ 0.67)."
```

---

## 🔬 Running Ensemble Analysis

### Configuration with ATO

```yaml
# ensemble_config_v2.1.1.yaml
ensemble:
  pipelines:
    - 1
    - 3
    - 4
  
  method: weighted_average
  
  weights:
    pipeline_1: 0.40
    pipeline_3: 0.35
    pipeline_4: 0.25
  
  min_pipelines: 2
  
  # NEW: Threshold Optimizer settings
  threshold_optimizer:
    enabled: true
    goal: "balanced"
    padj_method: "BH"
    logfc_method: "auto"
    uniform_thresholds: true
  
  detect_outliers: true
  generate_plots: true
  generate_report: true
```

### Run with Configuration

```bash
raptor ensemble \
  --config ensemble_config_v2.1.1.yaml \
  --data data/fastq/ \
  --output ensemble_analysis/
```

---

## 📊 Interpreting Results

### Ensemble Output Structure (v2.1.1)

```
ensemble_results/
├── consensus_genes.csv         # Final DE gene list
├── optimized_consensus.csv     # ATO-optimized results (NEW!)
├── threshold_report.txt        # ATO thresholds used (NEW!)
├── methods_text.txt            # Publication methods (NEW!)
├── pipeline_agreement.csv      # Agreement metrics
├── individual_results/         # Each pipeline's output
├── plots/
│   ├── venn_diagram.png
│   ├── volcano_optimized.png   # With ATO thresholds (NEW!)
│   └── concordance_heatmap.png
└── ensemble_report.html        # Interactive report
```

### Understanding optimized_consensus.csv (NEW)

```csv
gene_id,log2FC_consensus,FDR_consensus,n_pipelines,agreement_score,significant,direction
ENSG00000111640,3.45,1.2e-08,3,1.00,True,up
ENSG00000087086,3.21,2.3e-07,3,1.00,True,up
ENSG00000148584,-3.12,1.5e-07,3,1.00,True,down
ENSG00000183878,-0.65,3.2e-06,2,0.67,False,down  # Below ATO threshold
```

The `significant` column uses ATO-optimized thresholds!

---

## 🔧 Advanced Options

### Custom Weights with ATO

```python
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

analyzer = EnsembleAnalyzer()

# Combine with custom weights
ensemble = analyzer.combine(results, method='weighted_average', weights={
    1: 0.5,   # Trust STAR-RSEM most
    3: 0.3,
    4: 0.2
})

# Optimize thresholds
result = optimize_thresholds(ensemble['combined'], goal='validation')  # Stringent
```

### Sensitivity Analysis with ATO

```python
from raptor.ensemble import sensitivity_analysis
from raptor.threshold_optimizer import optimize_thresholds

# Test different ATO goals
results = {}
for goal in ['discovery', 'balanced', 'validation']:
    threshold_result = optimize_thresholds(ensemble_data, goal=goal)
    results[goal] = {
        'threshold': threshold_result.logfc_threshold,
        'n_genes': threshold_result.n_significant
    }
    
print("Goal\t\tlogFC\tGenes")
for goal, r in results.items():
    print(f"{goal}\t{r['threshold']:.2f}\t{r['n_genes']}")
```

---

## ✅ Best Practices

### Pipeline Selection

✅ **Use 3-5 pipelines** - Good balance  
✅ **Include diverse methods** - Both alignment and pseudo-alignment  
✅ **Mix statistical approaches** - Different DE methods  
✅ **Apply ATO for uniform thresholds** (NEW in v2.1.1)

**Recommended Combinations:**

**Standard (3 pipelines) with ATO:**
- Pipeline 1 (STAR-RSEM-DESeq2) - Gold standard
- Pipeline 3 (Salmon-edgeR) - Fast and accurate
- Pipeline 4 (Kallisto-DESeq2) - Ultra-fast
- **+ ATO goal='balanced'**

### Reporting in Publications (Updated for v2.1.1)

**Methods Section:**
```
Differential expression analysis was performed using ensemble 
analysis of three complementary pipelines: STAR-RSEM-DESeq2, 
Salmon-edgeR, and Kallisto-DESeq2. Results were combined using 
weighted averaging. Significance thresholds were determined using 
the Adaptive Threshold Optimizer (ATO) with the 'balanced' goal, 
yielding an optimized log₂ fold change threshold of 0.71 based on 
MAD estimation from the combined data distribution. Genes were 
considered differentially expressed if they met FDR < 0.05 and 
|log₂FC| > 0.71 in at least 2 of 3 pipelines (agreement score ≥ 0.67).
```

---

## 📚 Case Studies

### Case Study 1: Cancer vs Normal with ATO

**Setup:**
- 24 samples (12 cancer, 12 normal)
- High biological variation (BCV = 0.68)
- Clinical implications - need confidence

**Approach:**
```bash
raptor ensemble \
  --data cancer_study/ \
  --pipelines 1,3,5 \
  --use-ato \
  --ato-goal validation \  # Stringent for clinical
  --min-pipelines 3 \
  --output cancer_ensemble/
```

**Results:**
- ATO threshold: |logFC| > 0.89 (data-driven, not arbitrary)
- 1,523 high-confidence DE genes
- Publication-ready methods text generated ✅

### Case Study 2: Exploratory with Discovery Goal

**Setup:**
- Preliminary experiment
- Want to capture more candidates

**Approach:**
```python
result = optimize_thresholds(ensemble_combined, goal='discovery')
# More permissive threshold: |logFC| > 0.52
# Captures more potential hits for validation
```

---

## 🔧 Troubleshooting

### Issue: Different Pipelines Give Very Different Gene Counts

**Old Problem:**
```
Pipeline 1: 1,234 genes (|logFC| > 1)
Pipeline 3: 1,876 genes (|logFC| > 0.58)
Pipeline 4: 921 genes (|logFC| > 1, FDR < 0.01)
```

**v2.1.1 Solution - Use ATO:**
```python
from raptor.threshold_optimizer import optimize_thresholds

# Determine uniform threshold
result = optimize_thresholds(combined_data, goal='balanced')

# Apply to all pipelines
uniform_threshold = result.logfc_threshold
# Now: All pipelines use |logFC| > 0.73
```

### Issue: Too Few Consensus Genes

**Solutions:**
```bash
# 1. Use ATO with discovery goal
raptor ensemble --use-ato --ato-goal discovery

# 2. Lower agreement threshold
raptor ensemble --min-pipelines 2

# 3. Combine both approaches
```

### Issue: Ensemble Takes Too Long

**Solutions:**
```bash
# 1. Use fast pipelines
raptor ensemble --pipelines 3,4  # Pseudo-alignment only

# 2. Run ATO separately (fast)
# ATO takes <1 second regardless of ensemble size
```

---

## 📋 Summary

Ensemble analysis in v2.1.1 provides:
- ✅ **Robust Results** - Reduces false positives
- ✅ **Higher Confidence** - Multiple methods agree
- ✅ **🎯 Optimized Thresholds** - Data-driven via ATO (NEW!)
- ✅ **Publication Quality** - Auto-generated methods text
- ✅ **Uniform Analysis** - Same thresholds across pipelines

**Recommended workflow:**
1. Run ensemble with 3-5 pipelines
2. Apply ATO to combined results
3. Use optimized thresholds for final gene list
4. Include ATO methods text in publication

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.1  
**License:** MIT

---

*"Many pipelines, optimized thresholds, one truth!"* 🎯🦖
