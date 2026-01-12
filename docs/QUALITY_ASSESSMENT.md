# 🦖 RAPTOR v2.2.0 Quality Assessment Guide

**Comprehensive Data Quality Control and Assessment**

Automatically assess the quality of your RNA-seq data to ensure reliable results and identify potential issues before analysis.

---

## 🆕 What's New in v2.2.0

### Advanced Outlier Detection
- **6 detection methods** with consensus voting
- PCA + Mahalanobis, Isolation Forest, LOF, Elliptic Envelope, Correlation-based, Library Size
- Configurable consensus threshold
- Detailed per-method results

### PLS-DA Analysis
- **Supervised classification** for group separation
- **Cross-validation** with Q² calculation
- **Permutation testing** for model significance
- Component optimization

### VIP Scores
- **Variable Importance in Projection** for feature selection
- Biomarker discovery workflow
- Configurable significance thresholds
- Export to CSV

---

## 📑 Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Quality Metrics](#quality-metrics)
4. [Quality Scoring](#quality-scoring)
5. [Advanced Outlier Detection](#advanced-outlier-detection) ← **NEW**
6. [PLS-DA Analysis](#pls-da-analysis) ← **NEW**
7. [VIP Scores](#vip-scores) ← **NEW**
8. [Batch Effect Analysis](#batch-effect-analysis)
9. [Sample Quality](#sample-quality)
10. [Interpretation Guide](#interpretation-guide)
11. [Troubleshooting](#troubleshooting)
12. [Best Practices](#best-practices)

---

## 📋 Overview

### What is Quality Assessment?

Quality assessment evaluates your RNA-seq data to identify:

✅ **Low-quality samples** - Technical failures  
✅ **Outliers** - Unusual samples (6 detection methods)  
✅ **Batch effects** - Systematic biases  
✅ **Library problems** - Prep issues  
✅ **Biological signal** - Signal strength  
✅ **Biomarkers** - Key features via PLS-DA/VIP  

### Quality Assessment Workflow

```
Your Data → QC Check → Outlier Detection → Clean Data → PLS-DA/VIP → Reliable Results
            ↓              ↓                              ↓
         Issues?      Remove outliers            Biomarker candidates
```

---

## ⚡ Quick Start

### Basic Quality Assessment

```python
from raptor.data_quality_assessment import quick_quality_check

# Quick check with visualization
report = quick_quality_check(counts, metadata, output_file='quality.png')
print(f"Quality Score: {report['overall']['score']:.1f}/100")
```

### Advanced Outlier Detection (NEW)

```python
from raptor.data_quality_assessment import detect_outliers_quick

# Detect outliers using 6 methods
result = detect_outliers_quick(counts, consensus_threshold=3)
print(f"Outliers: {result.outlier_samples}")
```

### PLS-DA with VIP (NEW)

```python
from raptor.data_quality_assessment import run_plsda_quick

# Run PLS-DA and get VIP scores
result = run_plsda_quick(counts, metadata, group_column='condition')
biomarkers = result.get_significant_features(vip_threshold=1.5)
```

### Command Line

```bash
# Basic quality check
raptor qc --counts data/counts.csv --output qc_report/

# With outlier detection
raptor qc --counts data/counts.csv --outliers --threshold 3

# With PLS-DA
raptor plsda --counts data/counts.csv --metadata meta.csv --group condition
```

---

## 📊 Quality Metrics

### 1. Library Size Distribution

**What it measures:** Total counts per sample

```
Library Sizes (Million Reads)
    ↑
30M │     ▂▄█▆▃
25M │   ▃█████▆▂
20M │  ▅████████▄
15M │ ▇███████████▇
10M │██████████████
    └───────────────→ Samples

Mean: 24.5M | CV: 12.3%
```

**Interpretation:**
| CV | Status | Action |
|----|--------|--------|
| < 20% | ✅ Good | Proceed |
| 20-30% | ⚠️ Acceptable | Monitor |
| > 30% | ❌ Poor | Investigate |

---

### 2. Gene Detection Rate

**What it measures:** Number of expressed genes per sample

| Organism | Expected Genes | Status |
|----------|---------------|--------|
| Human | 15,000-20,000 | ✅ Good |
| Mouse | 13,000-18,000 | ✅ Good |
| Any | < 10,000 | ❌ Low complexity |

---

### 3. Zero Inflation

**What it measures:** Percentage of zero counts

| Zero % | Status | Notes |
|--------|--------|-------|
| 30-50% | ✅ Normal | Typical RNA-seq |
| < 30% | ⚠️ Unusual | Very deep sequencing? |
| > 60% | ❌ Poor | Low coverage |

---

### 4. Biological Coefficient of Variation (BCV)

**What it measures:** Biological variability between replicates

| BCV | Category | Typical Source |
|-----|----------|----------------|
| < 0.2 | Low | Cell lines |
| 0.2-0.6 | Medium | Most experiments |
| > 0.6 | High | Clinical samples |

---

## 🎯 Quality Scoring

### Overall Quality Score (0-100)

```
Component Weights:
├─ Library Sizes:       15%
├─ Gene Detection:      20%
├─ Outlier Detection:   15%
├─ Variance Structure:  15%
├─ Batch Effects:       20%
└─ Biological Signal:   15%
```

### Score Interpretation

| Score | Status | Action |
|-------|--------|--------|
| 90-100 | 🟢 Excellent | Proceed with confidence |
| 80-89 | 🟢 Good | Proceed, note minor issues |
| 70-79 | 🟡 Acceptable | Address flagged issues |
| 60-69 | 🟡 Marginal | Fix problems before analysis |
| < 60 | 🔴 Poor | Investigate, consider re-sequencing |

---

## 🔍 Advanced Outlier Detection

### NEW in v2.2.0

Detect outliers using **6 independent methods** with consensus voting.

### Methods

| Method | Description | Best For |
|--------|-------------|----------|
| **PCA + Mahalanobis** | Distance from center in PC space | Global outliers |
| **Isolation Forest** | Tree-based anomaly detection | Complex patterns |
| **Local Outlier Factor** | Local density comparison | Local outliers |
| **Elliptic Envelope** | Robust covariance estimation | Gaussian data |
| **Correlation-based** | Low correlation with other samples | Technical failures |
| **Library Size** | Z-score of total counts | Sequencing failures |

### Usage

```python
from raptor.data_quality_assessment import DataQualityAssessor

assessor = DataQualityAssessor(counts, metadata)

# Detect outliers (sample must be flagged by 3+ methods)
result = assessor.detect_outliers_advanced(
    methods=['pca_mahalanobis', 'isolation_forest', 'lof', 
             'elliptic_envelope', 'correlation', 'library_size'],
    consensus_threshold=3,
    contamination=0.1
)

# View results
print(result.summary())
print(f"Outliers: {result.outlier_samples}")

# See which methods flagged each sample
for sample, score in result.outlier_scores.items():
    if score >= 3:
        print(f"{sample}: flagged by {score}/6 methods")

# Visualize
assessor.plot_outliers_advanced(result, 'outliers.png')
```

### Output

```
============================================================
ADVANCED OUTLIER DETECTION RESULTS
============================================================
Consensus threshold: 3 methods
Total outliers: 2 (8.3%)

Outlier samples:
  • Sample_7 (flagged by 5 methods)
  • Sample_12 (flagged by 4 methods)

Detection by method:
  • PCA + Mahalanobis: 3 outliers
  • Isolation Forest: 2 outliers
  • Local Outlier Factor: 3 outliers
  • Elliptic Envelope: 2 outliers
  • Correlation-based: 2 outliers
  • Library Size: 1 outliers
============================================================
```

### Consensus Threshold Guide

| Threshold | Stringency | Use Case |
|-----------|------------|----------|
| 2 | Low | Exploratory, catch borderline |
| 3 | Medium | **Recommended default** |
| 4 | High | Conservative, clear outliers only |
| 5-6 | Very High | Only extreme outliers |

### Remove Outliers

```python
# Get clean data
clean_counts = counts.drop(columns=result.outlier_samples)
clean_metadata = metadata[~metadata['sample'].isin(result.outlier_samples)]

print(f"Removed {len(result.outlier_samples)} outliers")
print(f"Remaining: {clean_counts.shape[1]} samples")
```

---

## 📊 PLS-DA Analysis

### NEW in v2.2.0

**Partial Least Squares Discriminant Analysis** for supervised classification and biomarker discovery.

### What is PLS-DA?

PLS-DA finds the directions (components) in your data that best separate your experimental groups. It's useful for:

- **Classification**: Predicting group membership
- **Feature Selection**: Identifying important genes (via VIP)
- **Visualization**: Seeing group separation
- **Biomarker Discovery**: Finding discriminating features

### Usage

```python
from raptor.data_quality_assessment import DataQualityAssessor

assessor = DataQualityAssessor(counts, metadata)

# Run PLS-DA
result = assessor.run_plsda(
    group_column='condition',    # Column with group labels
    n_components=2,              # Number of components
    run_permutation=True,        # Test significance
    n_permutations=100           # Number of permutations
)

# View results
print(result.summary())

# Visualize
assessor.plot_plsda(result, 'plsda_results.png')
```

### Key Metrics

| Metric | Description | Good Value |
|--------|-------------|------------|
| **R²** | Variance explained | > 0.7 |
| **Q²** | Cross-validated R² | > 0.5 |
| **Accuracy** | Classification accuracy | > 80% |
| **CV Accuracy** | Cross-validated accuracy | > 70% |
| **Permutation p** | Model significance | < 0.05 |

### Output

```
============================================================
PLS-DA ANALYSIS RESULTS
============================================================

MODEL PERFORMANCE:
  • Components: 2
  • R² (explained variance): 0.847
  • Q² (cross-validated): 0.723
  • Accuracy: 91.7%
  • CV Accuracy: 83.3%

PERMUTATION TEST:
  • p-value: 0.0099
  • Status: ✓ Significant

VIP SCORES (Variable Importance in Projection):
  • Features with VIP > 1: 127

Top 10 features by VIP:
   1. GENE00234: 2.847
   2. GENE00891: 2.634
   3. GENE00456: 2.512
   ...
============================================================
```

### Optimize Components

```python
# Find optimal number of components
optimal_n = assessor.optimize_plsda_components(
    group_column='condition',
    max_components=10
)
print(f"Optimal components: {optimal_n}")

# Use optimal
result = assessor.run_plsda(n_components=optimal_n)
```

### Interpreting Q²

| Q² Value | Interpretation |
|----------|----------------|
| > 0.9 | Excellent predictive ability |
| 0.7 - 0.9 | Good predictive ability |
| 0.5 - 0.7 | Moderate predictive ability |
| 0.3 - 0.5 | Weak predictive ability |
| < 0.3 | Poor - model may be overfitting |

---

## 🎯 VIP Scores

### NEW in v2.2.0

**Variable Importance in Projection (VIP)** quantifies each feature's contribution to the PLS-DA model.

### Interpretation

| VIP Score | Significance |
|-----------|--------------|
| > 1.5 | ⭐⭐⭐ Highly significant |
| > 1.0 | ⭐⭐ Significant |
| 0.8 - 1.0 | ⭐ Moderately important |
| < 0.8 | Not significant |

### Get Significant Features

```python
# Get all VIP scores
all_vip = result.vip_scores  # pandas Series

# Get significant features (VIP > 1)
sig_features = result.get_significant_features(vip_threshold=1.0)

# Get highly significant features (VIP > 1.5)
high_sig = result.get_significant_features(vip_threshold=1.5)

print(f"VIP > 1.0: {len(sig_features)} features")
print(f"VIP > 1.5: {len(high_sig)} features")
```

### Quick VIP Selection

```python
from raptor.data_quality_assessment import get_vip_features

# One-liner for biomarker candidates
biomarkers = get_vip_features(
    counts, 
    metadata, 
    group_column='condition',
    vip_threshold=1.5
)

print(biomarkers.head(10))
```

### Export VIP Scores

```python
# Export to DataFrame
vip_df = result.to_dataframe()
vip_df.to_csv('vip_scores.csv', index=False)

# Output format:
# feature,vip_score,significant
# GENE001,2.847,True
# GENE002,2.634,True
# ...
```

### VIP for Biomarker Discovery

```python
# Complete biomarker discovery workflow
assessor = DataQualityAssessor(counts, metadata)

# 1. Run PLS-DA with permutation test
result = assessor.run_plsda(
    group_column='condition',
    n_components=2,
    run_permutation=True
)

# 2. Check if model is significant
if result.permutation_pvalue < 0.05:
    print("✓ Model is significant")
    
    # 3. Get biomarker candidates
    biomarkers = result.get_significant_features(vip_threshold=1.5)
    
    print(f"Found {len(biomarkers)} biomarker candidates:")
    for gene, vip in biomarkers.head(10).items():
        print(f"  {gene}: VIP = {vip:.3f}")
    
    # 4. Export
    biomarkers.to_csv('biomarker_candidates.csv')
else:
    print("✗ Model not significant - groups may not be separable")
```

---

## 🔬 Batch Effect Analysis

### Detecting Batch Effects

```python
assessor = DataQualityAssessor(counts, metadata)
report = assessor.assess_quality()

batch_info = report['components']['batch_effects']
print(f"Batch detected: {batch_info['batch_detected']}")
print(f"Batch variable: {batch_info['batch_variable']}")
print(f"Recommendation: {batch_info['recommendation']}")
```

### Batch Effect Severity

| Strength | Level | Action |
|----------|-------|--------|
| < 2 | Minor | Optional correction |
| 2-5 | Moderate | Include in model |
| 5-10 | Strong | Use ComBat |
| > 10 | Severe | Consider re-design |

---

## 👤 Sample Quality

### Individual Sample Assessment

```python
# Get per-sample quality metrics
lib_sizes = counts.sum(axis=0)
zero_pct = (counts == 0).sum(axis=0) / counts.shape[0] * 100

for sample in counts.columns:
    print(f"{sample}:")
    print(f"  Library size: {lib_sizes[sample]:,.0f}")
    print(f"  Zero %: {zero_pct[sample]:.1f}%")
```

### Sample Flags

| Issue | Threshold | Action |
|-------|-----------|--------|
| Low library size | < 30% of median | Exclude |
| High zeros | > 70% | Investigate |
| Low correlation | < 0.7 with group | Check |
| Multiple method outlier | ≥ 3 methods | Remove |

---

## 📖 Interpretation Guide

### Decision Tree

```
Start
  │
  ├─ Quality Score < 60?
  │   └─ Yes → Investigate issues, consider re-sequencing
  │
  ├─ Outliers detected?
  │   └─ Yes → Remove outliers, re-assess
  │
  ├─ Batch effects?
  │   └─ Yes → Include batch in model or correct
  │
  ├─ PLS-DA Q² < 0.5?
  │   └─ Yes → Groups may not be well separated
  │
  └─ Proceed with analysis
```

### Common Patterns

#### Pattern 1: Good Quality Data ✅
```
Quality: 85/100
Outliers: 0
Batch: None
PLS-DA Q²: 0.72
→ Proceed with any pipeline
```

#### Pattern 2: Outliers Present ⚠️
```
Quality: 72/100
Outliers: 2 samples (consensus 4+)
→ Remove outliers, re-assess
```

#### Pattern 3: Batch Effects ⚠️
```
Quality: 78/100
Batch: Strong (strength=8.5)
→ Use ~ batch + condition model
```

#### Pattern 4: Poor Separation ❌
```
Quality: 80/100
PLS-DA Q²: 0.25
Permutation p: 0.34
→ Groups not well separated, weak biological signal
```

---

## 🔧 Troubleshooting

### Issue: Too Many Outliers

```python
# Try higher consensus threshold
result = assessor.detect_outliers_advanced(consensus_threshold=4)

# Or adjust contamination
result = assessor.detect_outliers_advanced(contamination=0.05)
```

### Issue: PLS-DA Won't Converge

```python
# Reduce components
result = assessor.run_plsda(n_components=1)

# Or filter low-variance genes first
variance = counts.var(axis=1)
high_var = counts.loc[variance > variance.quantile(0.5)]
```

### Issue: VIP All Near 1

This suggests features contribute equally - possibly:
- No clear group separation
- Too many components
- Data normalization issues

```python
# Try fewer components
result = assessor.run_plsda(n_components=1)

# Or check Q² - if low, groups aren't separable
```

---

## ✅ Best Practices

### Quality Assessment Workflow

1. **Always run QC first**
   ```python
   report = quick_quality_check(counts, metadata)
   ```

2. **Check for outliers**
   ```python
   outliers = assessor.detect_outliers_advanced(consensus_threshold=3)
   ```

3. **Remove outliers before analysis**
   ```python
   clean_counts = counts.drop(columns=outliers.outlier_samples)
   ```

4. **Use PLS-DA for biomarkers**
   ```python
   result = assessor.run_plsda(run_permutation=True)
   biomarkers = result.get_significant_features(vip_threshold=1.5)
   ```

### For Publications

Always report:
- Overall quality score
- Number of outliers removed
- PLS-DA metrics (R², Q², permutation p-value)
- VIP threshold used

**Example methods text:**
```
RNA-seq data quality was assessed using RAPTOR v2.2.0. The overall 
quality score was 85/100. Two samples were identified as outliers 
by consensus of 4/6 detection methods (PCA-Mahalanobis, Isolation 
Forest, LOF, Elliptic Envelope, correlation-based, library size) 
and excluded. PLS-DA analysis (R²=0.85, Q²=0.72, permutation 
p<0.01) identified 45 features with VIP>1.5 as potential biomarkers.
```

---

## 📚 API Reference

### Classes

| Class | Description |
|-------|-------------|
| `DataQualityAssessor` | Main quality assessment class |
| `OutlierResult` | Container for outlier detection results |
| `PLSDAResult` | Container for PLS-DA results |

### Functions

| Function | Description |
|----------|-------------|
| `quick_quality_check()` | One-line quality assessment |
| `detect_outliers_quick()` | One-line outlier detection |
| `run_plsda_quick()` | One-line PLS-DA analysis |
| `get_vip_features()` | One-line VIP feature selection |

### DataQualityAssessor Methods

| Method | Description |
|--------|-------------|
| `assess_quality()` | Full quality assessment |
| `detect_outliers_advanced()` | 6-method outlier detection |
| `run_plsda()` | PLS-DA with VIP |
| `optimize_plsda_components()` | Find optimal components |
| `plot_quality_report()` | Quality visualization |
| `plot_outliers_advanced()` | Outlier visualization |
| `plot_plsda()` | PLS-DA visualization |

---

## 🎉 Summary

RAPTOR v2.2.0 Quality Assessment provides:

- ✅ **Comprehensive quality scoring** - 6 component metrics
- ✅ **Advanced outlier detection** - 6 methods with consensus
- ✅ **PLS-DA analysis** - Supervised classification
- ✅ **VIP scores** - Feature importance for biomarker discovery
- ✅ **Batch effect detection** - Identify systematic biases
- ✅ **Publication-ready visualizations** - Professional plots
- ✅ **Easy-to-use API** - Quick functions and detailed classes

**Good QC is the foundation of good analysis!** 🦖

---

**Author:** Ayeh Bolouki  
**Version:** 2.2.0  
**Email:** ayehbolouki1988@gmail.com  
**License:** MIT

---

*"Garbage in, garbage out - QC saves the day!"* 🦖
