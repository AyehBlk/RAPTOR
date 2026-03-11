# RAPTOR User Guide
**Complete Tutorial for RNA-seq Analysis Pipeline Testing and Optimization**

**Version:** 2.2.0  
**Author:** Ayeh Bolouki  
**Date:** March 2026

---

## 📚 **Table of Contents**

1. [Introduction](#1-introduction)
2. [Installation](#2-installation)
3. [Quick Start](#3-quick-start)
4. [Module-by-Module Guide](#4-module-by-module-guide)
5. [Pipeline Guide](#5-pipeline-guide)
6. [Dashboard Guide](#6-dashboard-guide)
7. [Troubleshooting](#7-troubleshooting)
8. [FAQ](#8-faq)

---

## 1. Introduction

### 1.1 What is RAPTOR?

RAPTOR (**R**NA-seq **A**nalysis **P**ipeline **T**esting and **O**ptimization **R**esource) is a comprehensive framework that makes sophisticated RNA-seq differential expression analysis accessible to researchers regardless of computational background.

**The Problem RAPTOR Solves:**

Traditional RNA-seq analysis requires you to make many decisions:
- Which quantification pipeline? (STAR? Salmon? Kallisto?)
- Which DE method? (DESeq2? edgeR? limma?)
- What thresholds? (FDR < 0.05? log2FC > 1?)
- How to combine results from multiple methods?

**These decisions are often based on:**
- ❌ What others in your lab use
- ❌ What you find in a tutorial
- ❌ Trial and error

**RAPTOR provides:**
- ✅ **Data-driven recommendations** based on your dataset characteristics
- ✅ **Automated optimization** of analysis parameters
- ✅ **Ensemble methods** to combine results from multiple approaches
- ✅ **Publication-quality** outputs with detailed QC

### 1.2 Key Features

**🔬 6 Production Pipelines**
- Salmon (recommended for most uses)
- Kallisto (fastest)
- STAR + featureCounts
- STAR + RSEM
- STAR + Salmon (unique: BAM + bootstraps)
- HISAT2 + featureCounts

**📊 32-Feature Data Profiling**
- Automatically characterizes your dataset
- Extracts key metrics (BCV, sample size, sparsity, etc.)
- Powers ML-based recommendations

**🤖 ML-Powered Recommendations**
- Random Forest classifier
- Trained on diverse datasets
- Recommends optimal pipeline and parameters

**⚙️ 4 Optimization Methods**
- Ground Truth (simulated data)
- FDR Control (Storey's method)
- Stability (cross-validation)
- Reproducibility (independent cohorts)

**🎯 5 Ensemble Methods**
- Fisher's Method
- Brown's Method
- Robust Rank Aggregation
- Voting Consensus
- Weighted Ensemble

**📈 Comprehensive QC**
- 6 outlier detection methods
- Consensus-based reporting
- Batch effect detection

### 1.3 Who Should Use RAPTOR?

**Biologists:**
- No need to be an expert in statistics or programming
- Clear recommendations for pipeline selection
- Interactive dashboard for exploration

**Bioinformaticians:**
- Systematic benchmarking framework
- Reproducible analysis workflows
- Ensemble approaches for robust results

**Computational Researchers:**
- Extensible architecture
- Python API for custom workflows
- Integration with existing pipelines

### 1.4 System Requirements

**Minimum:**
- OS: Linux, macOS, Windows (WSL2 for full features)
- Python: 3.8-3.12
- RAM: 4 GB (for Python package only)
- Disk: 500 MB (Python only)

**Recommended:**
- RAM: 16-32 GB (for alignment-based pipelines)
- Disk: 50-100 GB (for bioinformatics tools)
- CPU: 8+ cores

---

## 2. Installation

### 2.1 Method 1: PyPI (Quickest)

**For most users - just want the Python package:**

```bash
# Basic installation
pip install raptor-rnaseq

# Verify installation
raptor --version
python -c "import raptor; print(raptor.__version__)"
```

**Expected output:**
```
RAPTOR version 2.2.0
2.2.0
```

**With all optional features:**
```bash
pip install raptor-rnaseq[all]
```

**Development installation:**
```bash
pip install raptor-rnaseq[dev]
```

---

### 2.2 Method 2: Conda (Recommended for Bioinformatics)

**Two options:**

#### **Option A: Core Environment (Python only, 5-10 minutes)**

```bash
# Download environment file
wget https://raw.githubusercontent.com/AyehBlk/RAPTOR/main/environment.yml

# Create environment
conda env create -f environment.yml

# Activate
conda activate raptor

# Verify
raptor --version
```

**Size:** ~500 MB  
**Includes:** Python package only  
**Best for:** Development, testing, or if you already have bioinformatics tools

#### **Option B: Full Environment (Complete suite, 30-60 minutes)**

```bash
# Download environment file
wget https://raw.githubusercontent.com/AyehBlk/RAPTOR/main/environment-full.yml

# Create environment (grab coffee - this takes a while!)
conda env create -f environment-full.yml

# Activate
conda activate raptor-full

# Verify everything
raptor --version
STAR --version
salmon --version
kallisto version
```

**Size:** ~5-8 GB  
**Includes:** Python package + STAR, Salmon, Kallisto, HISAT2, SAMtools, R, DESeq2, edgeR, limma  
**Best for:** Complete analysis workstation, new users, workshops

**See [docs/CONDA_ENVIRONMENTS.md](docs/CONDA_ENVIRONMENTS.md) for detailed comparison.**

---

### 2.3 Method 3: From Source (For Developers)

```bash
# Clone repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Install in editable mode
pip install -e .

# Or with development tools
pip install -e .[dev]

# Verify
raptor --version

# Run tests
pytest tests/
```

---

### 2.4 Post-Installation Setup

**Verify installation:**

```bash
# Check CLI
raptor --help

# Check Python API
python -c "from raptor import ensemble_fisher; print('✅ API OK')"
```

---

## 3. Quick Start

**Complete RNA-seq analysis in 5 minutes!**

### 3.1 Prerequisites

You need:
- Count matrix (genes × samples)
- Metadata file (sample information)

**Example data format:**

**counts.csv:**
```csv
gene_id,sample_1,sample_2,sample_3,sample_4,sample_5,sample_6
ENSG00000000003,1523,1872,1456,2103,2987,2456
ENSG00000000005,0,0,0,0,0,0
ENSG00000000419,456,523,412,678,789,654
...
```

**metadata.csv:**
```csv
sample_id,condition,batch,replicate
sample_1,control,1,1
sample_2,control,1,2
sample_3,control,2,3
sample_4,treatment,1,1
sample_5,treatment,2,2
sample_6,treatment,2,3
```

### 3.2 5-Minute Workflow

```bash
# Step 1: Quality Control (30 seconds)
raptor qc \
    --counts counts.csv \
    --metadata metadata.csv \
    --output results/qc

# Step 2: Profile Your Data (1 minute)
raptor profile \
    --counts counts.csv \
    --metadata metadata.csv \
    --group-column condition \
    --output results/profile

# Step 3: Get Recommendation (10 seconds)
raptor recommend \
    --profile results/profile/data_profile.json \
    --method ml \
    --output results/recommendation

# Step 4: View recommendation
cat results/recommendation/recommendation.txt
```

**Expected output:**
```
╔═══════════════════════════════════════════╗
║     RAPTOR Pipeline Recommendation        ║
╠═══════════════════════════════════════════╣

  Recommended Pipeline: salmon
  Confidence: 0.89
  
  Reasoning:
    - BCV: 0.234 (moderate biological variation)
    - Sample size: 6 (adequate for parametric methods)
    - Sparsity: 42% (typical for bulk RNA-seq)
    - Need for isoforms: Yes (based on research questions)
  
  Alternative Options:
    1. kallisto (confidence: 0.78) - Faster, similar results
    2. star_featurecounts (confidence: 0.65) - Gene-level only
    
╚═══════════════════════════════════════════╝
```

### 3.3 Python Quick Start

```python
from raptor import (
    quick_quality_check,
    profile_data_quick,
    recommend_pipeline
)
import pandas as pd

# Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Step 1: QC
qc_report = quick_quality_check(
    counts='counts.csv',
    metadata='metadata.csv',
    output_dir='results/qc'
)
print(f"✅ QC complete. Outliers: {qc_report.outliers}")

# Step 2: Profile
profile = profile_data_quick(
    counts=counts,
    metadata=metadata,
    group_column='condition',
    output_dir='results/profile'
)
print(f"✅ Profile complete. BCV: {profile.features['bcv']:.3f}")

# Step 3: Recommend
recommendation = recommend_pipeline(
    profile_file='results/profile/data_profile.json',
    method='ml'
)
print(f"✅ Recommended pipeline: {recommendation.pipeline_name}")
print(f"   Confidence: {recommendation.confidence:.2f}")
```

---

## 4. Module-by-Module Guide

### 4.1 Module 2: Quality Assessment

**Purpose:** Identify outlier samples and assess data quality

#### **CLI Usage**

```bash
raptor qc \
    --counts counts.csv \
    --metadata metadata.csv \
    --output results/qc \
    --normalization tpm \
    --consensus-threshold 3 \
    --plot \
    --plot-output qc_plots.pdf
```

**Parameters:**
- `--counts`: Count matrix file (required)
- `--metadata`: Sample metadata (optional but recommended)
- `--output`: Output directory (default: results/qc)
- `--normalization`: tpm, fpkm, cpm, or counts (default: tpm)
- `--consensus-threshold`: Min methods to flag outlier (default: 3)
- `--plot`: Generate diagnostic plots
- `--plot-output`: Plot output file (default: qc_plots.pdf)

#### **Python API**

```python
from raptor import DataQualityAssessor, quick_quality_check

# Quick way
report = quick_quality_check(
    counts='counts.csv',
    metadata='metadata.csv',
    output_dir='qc_results/'
)

# Advanced way
assessor = DataQualityAssessor(
    normalization='tpm',
    consensus_threshold=3
)

report = assessor.assess_quality(
    counts=counts_df,
    metadata=metadata_df
)

# Results
print(f"Outliers: {report.outliers}")
print(f"Quality issues: {report.quality_issues}")
print(f"Recommendations: {report.recommendations}")

# Access metrics
metrics_df = report.metrics
```

#### **Understanding Results**

**Outlier Detection Methods (6):**
1. **MAD (Median Absolute Deviation)** - Statistical outliers
2. **Isolation Forest** - ML-based anomaly detection
3. **LOF (Local Outlier Factor)** - Density-based
4. **PCA-based** - Multivariate outliers
5. **Clustering-based** - Distance from cluster centers
6. **Statistical Distance** - Mahalanobis distance

**Sample flagged if ≥3 methods agree (consensus)**

**Output files:**
```
results/qc/
├── quality_report.json        # Complete QC report
├── outlier_samples.csv        # List of outliers
├── qc_metrics.csv             # All QC metrics
└── qc_plots.pdf               # Diagnostic plots
```

**Example interpretation:**
```python
if len(report.outliers) > 0:
    print(f"⚠️ Found {len(report.outliers)} outlier samples")
    print("Review these samples and consider removing them")
else:
    print("✅ No outliers detected - data quality looks good")
```

**📖 Full documentation:** [docs/MODULE_2_Quality_Assessment.md](docs/MODULE_2_Quality_Assessment.md)

---

### 4.2 Module 3: Data Profiler

**Purpose:** Extract 32 statistical features to characterize your dataset

#### **The Most Important Module for Recommendations!**

Module 3 calculates the **BCV (Biological Coefficient of Variation)** - the single most important metric for pipeline selection.

#### **CLI Usage**

```bash
raptor profile \
    --counts counts.csv \
    --metadata metadata.csv \
    --group-column condition \
    --output results/profile \
    --verbose
```

**Parameters:**
- `--counts`: Count matrix (required)
- `--metadata`: Sample metadata (required for group-based features)
- `--group-column`: Column in metadata defining groups (default: condition)
- `--min-count`: Minimum count threshold (default: 1)
- `--output`: Output directory
- `--verbose`: Show detailed progress

#### **Python API**

```python
from raptor import RNAseqDataProfiler, profile_data_quick

# Quick way
profile = profile_data_quick(
    counts='counts.csv',
    metadata='metadata.csv',
    group_column='condition',
    output_dir='profile_results/'
)

# Advanced way
profiler = RNAseqDataProfiler(
    min_count=1,
    verbose=True
)

profile = profiler.profile(
    counts=counts_df,
    metadata=metadata_df,
    group_column='condition'
)

# Access features
print(f"BCV: {profile.features['bcv']:.3f}")
print(f"Sample size: {profile.features['n_samples']}")
print(f"Sparsity: {profile.features['sparsity']:.1%}")

# Get all 32 features
all_features = profile.features
```

#### **Understanding the 32 Features**

**8 Categories:**

1. **Sample Characteristics (5 features)**
   - n_samples, n_groups, samples_per_group, min_group_size, group_balance

2. **Library Size (4 features)**
   - mean_library_size, cv_library_size, median_library_size, library_size_range

3. **Gene Detection (3 features)**
   - n_genes, detection_rate, reliable_detection_rate

4. **Expression Distribution (4 features)**
   - mean_expression, variance_expression, skewness, kurtosis

5. **Dispersion ⭐ MOST IMPORTANT (5 features)**
   - **bcv** ← Key metric!
   - common_dispersion
   - mean_dispersion
   - median_dispersion
   - max_dispersion

6. **Sparsity (3 features)**
   - sparsity, zero_inflation, low_count_proportion

7. **Count Distribution (4 features)**
   - zero_proportion, low_count_proportion, medium_count_proportion, high_count_proportion

8. **Mean-Variance Relationship (3 features)**
   - mv_slope, mv_r_squared, poisson_fit

#### **Understanding BCV**

**BCV (Biological Coefficient of Variation) = sqrt(dispersion)**

| BCV | Interpretation | Typical For | Best Method |
|-----|----------------|-------------|-------------|
| < 0.2 | Low variation | Technical replicates, cell lines | limma-voom, Wilcoxon |
| 0.2-0.4 | Moderate (typical) | Most biological experiments | DESeq2, edgeR, limma-voom |
| 0.4-0.6 | High variation | Human studies, heterogeneous samples | edgeR, edgeR_robust |
| > 0.6 | Very high | Check data quality first | edgeR_robust |

**Example:**
```python
bcv = profile.features['bcv']
if bcv < 0.2:
    print("Low biological variation - methods will agree")
elif bcv < 0.4:
    print("Typical biological variation - standard methods work")
elif bcv < 0.6:
    print("High biological variation - need robust methods")
else:
    print("⚠️ Very high BCV - check data quality")
```

**📖 Full documentation:** [docs/MODULE_3_Data_Profiling.md](docs/MODULE_3_Data_Profiling.md)  
**📖 Quick reference:** [docs/MODULE_3_QUICK_REFERENCE.md](docs/MODULE_3_QUICK_REFERENCE.md)

---

### 4.3 Module 4: Pipeline Recommender

**Purpose:** Get data-driven pipeline recommendation

#### **Two Methods:**

1. **Rule-based** - Expert rules based on BCV, sample size, etc.
2. **ML-based** - Random Forest trained on diverse datasets

#### **CLI Usage**

```bash
# ML-based (recommended)
raptor recommend \
    --profile results/profile/data_profile.json \
    --method ml \
    --output results/recommendation \
    --verbose-explanation

# Rule-based
raptor recommend \
    --profile results/profile/data_profile.json \
    --method rule-based

# Both methods
raptor recommend \
    --profile results/profile/data_profile.json \
    --method both
```

#### **Python API**

```python
from raptor import recommend_pipeline, MLRecommender

# Quick way
recommendation = recommend_pipeline(
    profile_file='profile.json',
    method='ml'
)

# ML-based (detailed)
ml_recommender = MLRecommender()
ml_recommendation = ml_recommender.recommend(profile_object)

print(f"Pipeline: {ml_recommendation.pipeline_name}")
print(f"Confidence: {ml_recommendation.confidence:.2f}")
print(f"Feature importance: {ml_recommendation.feature_importance}")

# Rule-based
from raptor import PipelineRecommender
recommender = PipelineRecommender()
rule_recommendation = recommender.recommend(profile_object)
```

#### **Understanding Recommendations**

**Recommendation object contains:**
- `pipeline_name`: Recommended pipeline
- `confidence`: Model confidence (0-1)
- `reasoning`: Why this pipeline was chosen
- `alternatives`: Other good options
- `feature_importance`: Which features influenced decision

**Example output:**
```
Recommended: salmon
Confidence: 0.89

Top 3 influential features:
  1. BCV (0.234) → moderate variation
  2. Sample size (6) → adequate for parametric
  3. Sparsity (0.42) → typical bulk RNA-seq

Alternatives:
  - kallisto (confidence: 0.78) - faster, similar
  - star_featurecounts (confidence: 0.65) - gene-level only
```

**📖 Full documentation:** [docs/MODULE_4_Pipeline_Recommender.md](docs/MODULE_4_Pipeline_Recommender.md)

---

### 4.4 Module 7: DE Import

**Purpose:** Import DE results from any tool into standardized format

#### **Supported Formats:**
- DESeq2
- edgeR
- limma-voom
- Wilcoxon
- Custom (any tool)

#### **CLI Usage**

```bash
# Import DESeq2 results
raptor import-de \
    --input deseq2_results.csv \
    --method deseq2 \
    --output results/de_imported

# Import edgeR
raptor import-de \
    --input edger_results.csv \
    --method edger

# Import limma
raptor import-de \
    --input limma_results.csv \
    --method limma

# Auto-detect
raptor import-de \
    --input results.csv \
    --method auto

# Compare multiple results
raptor compare-de \
    deseq2.csv edger.csv limma.csv \
    --output results/comparison \
    --venn-diagram
```

#### **Python API**

```python
from raptor import (
    import_deseq2,
    import_edger,
    import_limma,
    compare_de_results
)

# Import DESeq2
deseq2_result = import_deseq2(
    'deseq2_results.csv',
    gene_column='gene_id',
    pvalue_column='pvalue',
    padj_column='padj',
    lfc_column='log2FoldChange'
)

# Import edgeR (auto-detects columns)
edger_result = import_edger('edger_results.csv')

# Import limma
limma_result = import_limma('limma_results.csv')

# Compare all three
comparison = compare_de_results(
    deseq2=deseq2_result,
    edger=edger_result,
    limma=limma_result,
    threshold=0.05
)

print(f"Overlap: {comparison.overlap_matrix}")
print(f"Consensus genes: {len(comparison.consensus_genes)}")
```

#### **Standardized Format**

All imports produce DEResult object with:
- `gene_id`: Gene identifier
- `log2FoldChange`: Effect size
- `pvalue`: P-value
- `padj`: Adjusted p-value (FDR)
- `baseMean`: Mean expression (or AveExpr for limma)

This allows easy comparison across methods!

**📖 Full documentation:** [docs/MODULE_7_DE_Import.md](docs/MODULE_7_DE_Import.md)

---

### 4.5 Module 8: Parameter Optimization

**Purpose:** Optimize FDR and log2FC thresholds

#### **4 Methods:**

1. **Ground Truth** - Have simulated/spike-in data with known true positives
2. **FDR Control** - Want to achieve specific FDR (e.g., 0.05)
3. **Stability** - Maximize stability across cross-validation folds
4. **Reproducibility** - Maximize agreement with independent cohort

#### **CLI Usage**

```bash
# Method 1: Ground Truth
raptor optimize \
    --de-result de_results.csv \
    --method ground-truth \
    --ground-truth true_positives.csv \
    --output results/optimization

# Method 2: FDR Control
raptor optimize \
    --de-result de_results.csv \
    --method fdr-control \
    --fdr-target 0.05 \
    --output results/optimization

# Method 3: Stability
raptor optimize \
    --de-result de_results.csv \
    --method stability \
    --counts counts.csv \
    --metadata metadata.csv \
    --output results/optimization

# Method 4: Reproducibility
raptor optimize \
    --de-result de_results.csv \
    --method reproducibility \
    --counts counts_cohort1.csv \
    --cohort2 counts_cohort2.csv \
    --output results/optimization
```

#### **Python API**

```python
from raptor import (
    optimize_with_ground_truth,
    optimize_with_fdr_control,
    optimize_with_stability,
    optimize_with_reproducibility
)

# Method 1: Ground Truth
result = optimize_with_ground_truth(
    de_result=de_result_object,
    ground_truth=true_positives,
    output_dir='optimization/ground_truth'
)

# Method 2: FDR Control
result = optimize_with_fdr_control(
    de_result=de_result_object,
    fdr_target=0.05,
    output_dir='optimization/fdr'
)

# Method 3: Stability
result = optimize_with_stability(
    counts=counts_df,
    metadata=metadata_df,
    output_dir='optimization/stability'
)

# Method 4: Reproducibility
result = optimize_with_reproducibility(
    counts=counts_cohort1,
    metadata=metadata1,
    cohort2=counts_cohort2,
    output_dir='optimization/reproducibility'
)

# Results
print(f"Optimal FDR threshold: {result.optimal_threshold['padj']}")
print(f"Optimal log2FC threshold: {result.optimal_threshold['lfc']}")
print(f"Performance metrics: {result.metrics}")
```

#### **Understanding Results**

**OptimizationResult contains:**
- `optimal_threshold`: Best threshold values
- `metrics`: Performance at optimal threshold
- `search_history`: Optimization trajectory
- `recommendations`: Suggested thresholds

**Example:**
```
Optimal Thresholds:
  FDR (padj): 0.05
  log2FC: 0.58

Performance:
  F1 Score: 0.847
  Precision: 0.901
  Recall: 0.798
  FDR (actual): 0.049

Number of DE genes: 1,247
```

**📖 Full documentation:** [docs/MODULE_8_Parameter_Optimization.md](docs/MODULE_8_Parameter_Optimization.md)

---

### 4.6 Module 9: Ensemble Analysis

**Purpose:** Combine DE results from multiple methods for robust conclusions

#### **5 Ensemble Methods:**

1. **Fisher's Method** - Classic p-value combination (assumes independence)
2. **Brown's Method** - Accounts for correlation between methods
3. **RRA (Robust Rank Aggregation)** - Rank-based, robust to outliers
4. **Voting** - Gene significant in ≥K methods
5. **Weighted** - Use optimization weights from Module 8

#### **CLI Usage**

```bash
# Single method (Fisher's)
raptor ensemble \
    --methods fisher \
    --deseq2 deseq2.csv \
    --edger edger.csv \
    --limma limma.csv \
    --output results/ensemble \
    --threshold 0.05

# Compare all 5 methods
raptor ensemble-compare \
    --deseq2 deseq2.csv \
    --edger edger.csv \
    --limma limma.csv \
    --output results/ensemble_comparison \
    --threshold 0.05
```

#### **Python API**

```python
from raptor import (
    ensemble_fisher,
    ensemble_brown,
    ensemble_rra,
    ensemble_voting,
    ensemble_weighted
)

# Prepare results dictionary
de_results = {
    'deseq2': deseq2_result,
    'edger': edger_result,
    'limma': limma_result
}

# Method 1: Fisher's
fisher_result = ensemble_fisher(
    de_results_dict=de_results,
    significance_threshold=0.05
)

# Method 2: Brown's (correlation-aware)
brown_result = ensemble_brown(
    de_results_dict=de_results,
    significance_threshold=0.05
)

# Method 3: RRA
rra_result = ensemble_rra(
    de_results_dict=de_results,
    significance_threshold=0.05
)

# Method 4: Voting
voting_result = ensemble_voting(
    de_results_dict=de_results,
    min_methods=2,  # Gene must be in ≥2 methods
    significance_threshold=0.05
)

# Method 5: Weighted
weighted_result = ensemble_weighted(
    de_results_dict=de_results,
    weights={'deseq2': 0.4, 'edger': 0.3, 'limma': 0.3},
    significance_threshold=0.05
)

# Results
print(f"Consensus genes: {len(fisher_result.consensus_genes)}")
print(f"Direction consistency: {fisher_result.direction_consistency}")
```

#### **Choosing an Ensemble Method**

| Method | When to Use | Pros | Cons |
|--------|-------------|------|------|
| **Fisher's** | Methods are independent | Simple, well-established | Assumes independence |
| **Brown's** | Methods are correlated (DESeq2+edgeR) | Accounts for correlation | More complex |
| **RRA** | Effect sizes vary widely | Robust to outliers | Loses magnitude info |
| **Voting** | Want conservative consensus | Very robust | Can be too conservative |
| **Weighted** | You trust one method more | Uses optimization weights | Need to determine weights |

**Recommendation:** Start with **Brown's method** (accounts for correlation between DESeq2 and edgeR).

#### **Understanding Results**

**EnsembleResult contains:**
- `consensus_genes`: Genes passing ensemble threshold
- `combined_pvalues`: Combined p-value per gene
- `meta_lfc`: Meta-analysis log fold change
- `direction_consistency`: Direction agreement across methods
- `method_agreement`: Per-gene method agreement

**Example:**
```
Ensemble Analysis (Brown's Method)
===================================
Input methods: DESeq2, edgeR, limma (3)
Total genes analyzed: 20,145

Results:
  Consensus DE genes: 1,847
  Direction consistent: 1,823 (98.7%)
  
Agreement:
  All 3 methods: 1,456 genes (79%)
  2 of 3 methods: 391 genes (21%)
  
Correlation between methods:
  DESeq2-edgeR: 0.92
  DESeq2-limma: 0.85
  edgeR-limma: 0.88
```

**📖 Full documentation:** [docs/MODULE_9_Ensemble_Analysis.md](docs/MODULE_9_Ensemble_Analysis.md)

---

## 5. Pipeline Guide

### 5.1 Overview of 6 Pipelines

| Pipeline | Memory | Time | Produces | Best For |
|----------|--------|------|----------|----------|
| **salmon** ⭐ | 8 GB | 10-20 min | genes + isoforms + bootstraps | **Recommended** |
| **kallisto** | 4 GB | 5-10 min | genes + isoforms + bootstraps | Speed priority |
| **star_featurecounts** | 32 GB | 40-70 min | BAM + genes | Gene-level publication |
| **star_rsem** | 32 GB | 60-120 min | BAM + genes + isoforms | Isoform quantification |
| **star_salmon** | 32 GB | 50-90 min | BAM + genes + isoforms + bootstraps | Need BAM + bootstraps |
| **hisat2_featurecounts** | 16 GB | 30-60 min | BAM + genes | Low memory systems |

### 5.2 Running Pipelines

#### **General Workflow:**

```bash
# 1. List available pipelines
raptor pipeline list

# 2. Run specific pipeline
raptor pipeline run \
    --name salmon \
    --samples samples.csv \
    --genome-index /path/to/salmon_index \
    --output results/salmon \
    --threads 8
```

#### **Sample Sheet Format:**

**samples.csv:**
```csv
sample_id,fastq_1,fastq_2,condition
sample_1,data/sample_1_R1.fq.gz,data/sample_1_R2.fq.gz,control
sample_2,data/sample_2_R1.fq.gz,data/sample_2_R2.fq.gz,control
sample_3,data/sample_3_R1.fq.gz,data/sample_3_R2.fq.gz,treatment
```

### 5.3 Pipeline-Specific Examples

#### **Salmon (Recommended)**

```bash
raptor pipeline run \
    --name salmon \
    --samples samples.csv \
    --genome-index salmon_index/ \
    --output results/salmon \
    --threads 8 \
    --bootstraps 100  # For uncertainty quantification
```

**Outputs:**
```
results/salmon/
├── counts/
│   ├── gene_counts.csv
│   ├── transcript_counts.csv
│   └── tx2gene.tsv
├── quant/
│   ├── sample_1/
│   ├── sample_2/
│   └── sample_3/
└── qc/
    └── salmon_qc_report.html
```

#### **Kallisto (Fastest)**

```bash
raptor pipeline run \
    --name kallisto \
    --samples samples.csv \
    --genome-index kallisto_index/ \
    --output results/kallisto \
    --threads 8 \
    --bootstraps 100
```

#### **STAR + featureCounts (Gene-level)**

```bash
raptor pipeline run \
    --name star_featurecounts \
    --samples samples.csv \
    --genome-index STAR_index/ \
    --annotation genes.gtf \
    --output results/star_fc \
    --threads 16  # STAR needs more memory/threads
```

**Outputs BAM files** for visualization in IGV!

### 5.4 Pipeline Comparison

After running multiple pipelines:

```python
from raptor.pipelines import compare_pipelines

comparison = compare_pipelines([
    'results/salmon/counts/gene_counts.csv',
    'results/kallisto/counts/gene_counts.csv',
    'results/star_fc/counts/gene_counts.csv'
])

print(comparison.correlation_matrix)
print(comparison.agreement_stats)
```

---

## 6. Dashboard Guide

### 6.1 Launching the Dashboard

```bash
# Activate environment with streamlit
conda activate raptor

# Launch dashboard
streamlit run raptor/dashboard/app.py

# Or with custom port
streamlit run raptor/dashboard/app.py --server.port 8080
```

Opens in browser at `http://localhost:8501`

### 6.2 Dashboard Pages

**Home/Welcome**
- Overview of RAPTOR
- Quick links
- System status

**📊 Quality Assessment**
- Upload count matrix + metadata
- Run 6 outlier detection methods
- Interactive QC plots
- Download report

**📈 Data Profiler**
- Upload data
- View 32 features
- Visualize BCV, sparsity, etc.
- Compare to reference datasets

**🤖 Recommender**
- Upload profile
- Get ML and rule-based recommendations
- View confidence scores
- See feature importance

**📥 DE Import**
- Upload DE results from any tool
- Standardize format
- Compare multiple methods
- Venn diagrams

**⚙️ Optimization**
- Select optimization method
- Upload data
- Run optimization
- Visualize search trajectory

**🎯 Ensemble**
- Upload multiple DE results
- Select ensemble method
- View consensus genes
- Direction consistency checks
- Download combined results

**📝 Reports**
- Generate publication-quality reports
- Export to PDF/HTML
- Include QC, profiling, recommendations
- Custom branding

**⚙️ Settings**
- Configure analysis parameters
- Set default thresholds
- Manage projects

**📊 Visualization**
- Interactive plots
- MA plots, volcano plots
- Heatmaps
- PCA

### 6.3 Dashboard Workflow Example

```
1. Upload Data (Quality page)
   → Upload counts.csv, metadata.csv
   → Click "Run QC"
   → Review outliers

2. Profile Data (Profiler page)
   → Automatically uses QC'd data
   → Click "Profile Dataset"
   → View 32 features, BCV

3. Get Recommendation (Recommender page)
   → Automatically uses profile
   → Click "Get Recommendation"
   → View recommended pipeline

4. Compare DE Results (Ensemble page)
   → Upload DESeq2, edgeR, limma results
   → Select "Brown's method"
   → Click "Run Ensemble"
   → Download consensus genes

5. Generate Report (Reports page)
   → Select modules to include
   → Click "Generate Report"
   → Download PDF
```

---

## 7. Troubleshooting

### 7.1 Installation Issues

#### **Problem:** `pip install raptor-rnaseq` fails

**Solutions:**
```bash
# Update pip
pip install --upgrade pip

# Try with --no-cache
pip install --no-cache-dir raptor-rnaseq

# Install specific version
pip install raptor-rnaseq==2.2.0
```

#### **Problem:** Conda environment creation stuck

**Solutions:**
```bash
# Use mamba (faster)
conda install -c conda-forge mamba
mamba env create -f environment.yml

# Or install in stages
conda create -n raptor python=3.10
conda activate raptor
pip install raptor-rnaseq
```

#### **Problem:** Import errors

```python
ImportError: No module named 'raptor'
```

**Solutions:**
```bash
# Verify installation
pip list | grep raptor

# Reinstall
pip uninstall raptor-rnaseq
pip install raptor-rnaseq

# Check Python version
python --version  # Must be 3.8-3.12
```

### 7.2 Data Issues

#### **Problem:** "Too few genes" error

**Solution:**
```python
# Check your count matrix
import pandas as pd
counts = pd.read_csv('counts.csv', index_col=0)
print(f"Genes: {len(counts)}")  # Should be >100

# Remove low-count genes
counts_filtered = counts[counts.sum(axis=1) > 10]
```

#### **Problem:** "Sample mismatch" error

**Solution:**
```python
# Check sample IDs match
counts_samples = set(counts.columns)
metadata_samples = set(metadata['sample_id'])

missing = counts_samples - metadata_samples
print(f"Missing in metadata: {missing}")

extra = metadata_samples - counts_samples
print(f"Missing in counts: {extra}")
```

#### **Problem:** High BCV (>0.6)

**Check:**
1. Are there batch effects? (Module 2 QC)
2. Are there outliers? (Module 2 QC)
3. Is data quality good? (Check FastQC reports)
4. Are samples actually heterogeneous? (This might be real biology!)

### 7.3 Analysis Issues

#### **Problem:** No consensus genes in ensemble

**Possible causes:**
- Methods disagree too much → Check correlation
- Threshold too strict → Try 0.1 instead of 0.05
- Different normalizations → Re-run with same normalization

**Solution:**
```python
# Check method agreement
comparison = compare_de_results(deseq2, edger, limma)
print(comparison.correlation_matrix)

# If correlation < 0.7, methods disagree significantly
# Try voting with min_methods=2 instead of 3
```

#### **Problem:** Optimization gives extreme thresholds

**Possible causes:**
- Insufficient data (need ≥50 validated genes)
- Poor ground truth
- Unstable results

**Solution:**
```python
# Check optimization trajectory
print(result.search_history)

# If unstable, use FDR control method instead
result = optimize_with_fdr_control(
    de_result=de_result,
    fdr_target=0.05
)
```

### 7.4 Performance Issues

#### **Problem:** Pipeline very slow

**Solutions:**
- Use more threads: `--threads 16`
- Use faster pipeline: kallisto instead of STAR
- Check disk I/O (is data on slow drive?)
- Use SSD for temporary files

#### **Problem:** Out of memory

**Solutions:**
- Use lower-memory pipeline (kallisto, salmon)
- Increase system RAM
- Process samples in batches
- Use HPC cluster

### 7.5 Output Issues

#### **Problem:** Can't find output files

**Check:**
```bash
# List output directory
ls -lR results/

# Check for error logs
cat results/*/logs/*.err

# Verify permissions
ls -l results/
```

#### **Problem:** Output format wrong

**For standard format:**
```python
from raptor import import_de_result

# Auto-detect format
result = import_de_result(
    'my_results.csv',
    method='auto'
)

# Check what was detected
print(result.method_detected)
```

---

## 8. FAQ

### 8.1 General Questions

**Q: What is RAPTOR?**
A: RAPTOR is a comprehensive framework for RNA-seq analysis that provides ML-powered pipeline recommendations, threshold optimization, and ensemble methods to combine results from multiple DE tools.

**Q: Do I need to know programming?**
A: No! RAPTOR has:
- Command-line interface (just copy/paste commands)
- Interactive dashboard (point-and-click)
- Python API (for advanced users)

**Q: What data do I need?**
A: Minimum:
- Count matrix (genes × samples)
- Metadata file (sample information)

**Q: Can I use my existing DE results?**
A: Yes! Module 7 imports results from DESeq2, edgeR, limma, or any custom format.

**Q: Is RAPTOR free?**
A: Yes! MIT license, completely free and open-source.

### 8.2 Installation Questions

**Q: Which installation method should I use?**
A:
- **pip:** Fastest, Python package only
- **conda core:** Python + optional tools
- **conda full:** Everything included (recommended for beginners)
- **source:** For developers

**Q: Do I need bioinformatics tools installed?**
A: Only if you want to run Module 5 (quantification pipelines). Modules 2-4, 7-9 work with count matrices.

**Q: Can I use RAPTOR on Windows?**
A: Python package (Modules 2-4, 7-9): Yes
Bioinformatics tools (Module 5): Use WSL2 or Docker

**Q: How much disk space do I need?**
A: 
- Python only: ~500 MB
- Full environment: ~5-8 GB
- Analysis data: Varies (plan for 50-100 GB)

### 8.3 Analysis Questions

**Q: Which pipeline should I use?**
A: Use Module 4 to get a data-driven recommendation! Generally:
- **salmon** for most uses (recommended)
- **kallisto** if speed is priority
- **STAR-based** if you need BAM files

**Q: Which DE method is best?**
A: No single "best" method. That's why RAPTOR provides ensemble analysis! Use:
- **DESeq2** for general use
- **edgeR** for small samples or high variation
- **limma-voom** for large samples
- **Ensemble** for robust results

**Q: What is BCV and why does it matter?**
A: BCV (Biological Coefficient of Variation) = sqrt(dispersion)
- Measures biological variability
- Key metric for method selection
- BCV = 0.2-0.4 is typical
- See Module 3 docs for details

**Q: How do I optimize thresholds?**
A: Module 8 provides 4 methods:
- **FDR Control** (most common) - achieve target FDR
- **Ground Truth** - if you have simulated data
- **Stability** - maximize reproducibility
- **Reproducibility** - if you have replicates

**Q: Should I use ensemble methods?**
A: Yes, especially for important findings! Ensemble methods:
- Increase confidence in results
- Reduce false positives
- Find robust consensus genes
- Publication reviewers love them

### 8.4 Technical Questions

**Q: What are the 32 profiling features?**
A: See Module 3 docs. Key categories:
- Sample characteristics
- Library size
- Dispersion (BCV)
- Sparsity
- Expression distribution

**Q: How does the ML recommender work?**
A: Random Forest classifier trained on:
- 32 profiling features
- Diverse datasets
- Validated pipeline performance
- Returns confidence score + reasoning

**Q: Can I add my own pipelines?**
A: Yes! RAPTOR is extensible:
```python
from raptor.pipelines import BasePipeline

class MyPipeline(BasePipeline):
    # Implement required methods
    pass
```

**Q: Can I use custom DE methods?**
A: Yes! Import with Module 7:
```python
from raptor import import_de_result

result = import_de_result(
    'custom_results.csv',
    method='custom',
    gene_column='gene',
    pvalue_column='p',
    lfc_column='logFC'
)
```

### 8.5 Performance Questions

**Q: How long does analysis take?**
A:
- QC (Module 2): 1-5 minutes
- Profiling (Module 3): 1-2 minutes
- Recommendation (Module 4): <10 seconds
- Quantification (Module 5): 5-120 minutes (pipeline-dependent)
- Optimization (Module 8): 5-30 minutes
- Ensemble (Module 9): <1 minute

**Q: How much memory do I need?**
A:
- Modules 2-4, 7-9: 4-8 GB RAM
- Salmon/Kallisto: 8-16 GB RAM
- STAR-based: 32-64 GB RAM

**Q: Can I run on HPC cluster?**
A: Yes! RAPTOR CLI commands work on any system. Use SLURM/PBS job scripts for pipelines.

### 8.6 Results Questions

**Q: How do I interpret the profile?**
A: Key metrics:
- **BCV < 0.4:** Good, typical variation
- **Min group size ≥ 3:** Adequate samples
- **Sparsity 30-60%:** Normal for bulk RNA-seq
- See Module 3 Quick Reference

**Q: What if I have outliers?**
A: Module 2 detects outliers. Options:
1. Remove them (if technical issue)
2. Keep them (if biological variation)
3. Use robust methods (edgeR_robust)

**Q: How do I choose ensemble method?**
A: 
- **Brown's:** Best for correlated methods (DESeq2+edgeR)
- **Fisher's:** If methods independent
- **Voting:** Most conservative
- **RRA:** Robust to outliers
- **Weighted:** If you trust one method more

**Q: What is a good number of DE genes?**
A: Depends on biology! Typical ranges:
- 100-1000: Small effect
- 1000-3000: Moderate effect
- 3000+: Large effect
- 10000+: Check for batch effects!

### 8.7 Compatibility Questions

**Q: What file formats are supported?**
A:
- **Counts:** CSV, TSV, Excel
- **Metadata:** CSV, TSV
- **DE results:** CSV, TSV, Excel
- **Annotations:** GTF, GFF

**Q: Can I use with 10X Genomics data?**
A: Yes, but:
- Single-cell needs specialized methods (Seurat, Scanpy)
- RAPTOR optimized for bulk RNA-seq
- Can use for pseudo-bulk analysis

**Q: Does RAPTOR work with non-model organisms?**
A: Yes! Modules 2-4, 7-9 work with any organism.
Module 5 (pipelines) needs:
- Transcriptome/genome FASTA
- Annotation GTF (for gene-level)

**Q: Can I analyze microarray data?**
A: Not directly. RAPTOR is designed for count data (RNA-seq, etc.).

### 8.8 Publication Questions

**Q: How do I cite RAPTOR?**
A: See [CITATION.cff](CITATION.cff) or README for BibTeX.

**Q: Can I use RAPTOR for publication?**
A: Yes! RAPTOR is designed for publication-quality analysis:
- Validated methods
- Reproducible workflows
- Detailed QC reports
- Publication-ready figures

**Q: What should I include in Methods section?**
A: Example:
```
RNA-seq data quality was assessed using RAPTOR v2.2.0 
(Bolouki, 2026). Dataset profiling identified moderate 
biological variation (BCV=0.234), and machine learning 
recommendation suggested Salmon quantification. Differential 
expression was performed using DESeq2, edgeR, and limma-voom, 
with consensus genes identified via Brown's ensemble method 
(accounting for method correlation). Optimal thresholds 
(FDR < 0.05, |log2FC| > 0.58) were determined via FDR control 
optimization.
```

**Q: Is there a preprint/paper?**
A: Coming soon! Check GitHub for updates.

---

## 📚 Additional Resources

### Complete Module Documentation
- [MODULE_1_Quick_Quantification.md](docs/MODULE_1_Quick_Quantification.md)
- [MODULE_2_Quality_Assessment.md](docs/MODULE_2_Quality_Assessment.md)
- [MODULE_3_Data_Profiling.md](docs/MODULE_3_Data_Profiling.md)
- [MODULE_3_QUICK_REFERENCE.md](docs/MODULE_3_QUICK_REFERENCE.md)
- [MODULE_4_Pipeline_Recommender.md](docs/MODULE_4_Pipeline_Recommender.md)
- [MODULE_7_DE_Import.md](docs/MODULE_7_DE_Import.md)
- [MODULE_8_Parameter_Optimization.md](docs/MODULE_8_Parameter_Optimization.md)
- [MODULE_9_Ensemble_Analysis.md](docs/MODULE_9_Ensemble_Analysis.md)
- [CONDA_ENVIRONMENTS.md](docs/CONDA_ENVIRONMENTS.md)

### Example Scripts
- [examples/02_quality_assessment.py](examples/02_quality_assessment.py)
- [examples/03_data_profiler.py](examples/03_data_profiler.py)
- [examples/04_recommender.py](examples/04_recommender.py)
- [examples/07_DE_Import.py](examples/07_DE_Import.py)
- [examples/08_Parameter_Optimization.py](examples/08_Parameter_Optimization.py)
- [examples/09_Ensemble_Analysis.py](examples/09_Ensemble_Analysis.py)

### Online Resources
- **GitHub:** https://github.com/AyehBlk/RAPTOR
- **Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Discussions:** https://github.com/AyehBlk/RAPTOR/discussions
- **PyPI:** https://pypi.org/project/raptor-rnaseq/

### Get Help
- 📧 Email: ayehbolouki1988@gmail.com
- 🐛 Report bugs: [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- 💬 Ask questions: [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)

---

**Happy analyzing! 🦖**

**Last updated:** March 2026  
**Version:** 2.2.0
