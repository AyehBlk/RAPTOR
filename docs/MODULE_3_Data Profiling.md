# RAPTOR Module 3: Data Profiling & Feature Extraction

**Version 2.2.0** | **Stage 3** | **Module M3**

---

## Table of Contents

1. [Overview](#overview)
2. [Philosophy & Design](#philosophy--design)
3. [Quick Start](#quick-start)
4. [The 32 Features](#the-32-features)
5. [Scientific Basis](#scientific-basis)
6. [Usage Examples](#usage-examples)
7. [Interpreting Profiles](#interpreting-profiles)
8. [Integration with RAPTOR Workflow](#integration-with-raptor-workflow)
9. [Feature Categories](#feature-categories)
10. [Mathematical Formulations](#mathematical-formulations)
11. [Troubleshooting](#troubleshooting)
12. [Best Practices](#best-practices)
13. [Technical Details](#technical-details)

---

## Overview

### What is Module 3?

Module 3 performs **comprehensive statistical profiling** of RNA-seq count matrices to extract features that are predictive of optimal differential expression pipeline performance. It characterizes your data across multiple dimensions to enable intelligent, ML-based tool recommendations.

### Key Characteristics

- 📊 **32 Features Extracted**: Organized into 8 categories
- 🔬 **Literature-Based**: Features derived from RNA-seq benchmarking studies
- 🤖 **ML-Ready**: Output formatted for machine learning recommendation
- 🎯 **Predictive**: Features correlated with pipeline performance
- 📈 **Comprehensive**: Covers sample size, dispersion, sparsity, quality
- ⚡ **Fast**: Profiles 50k genes × 100 samples in < 1 minute

### What Module 3 Does

✅ Extract sample and group characteristics  
✅ Calculate library size metrics  
✅ Assess gene detection patterns  
✅ Analyze expression distribution  
✅ **Estimate dispersion (CRITICAL for pipeline selection)**  
✅ Measure sparsity and zero-inflation  
✅ Profile count distribution  
✅ Analyze mean-variance relationship  
✅ Generate feature vector for ML recommendation  
✅ Produce human-readable summary

### What Module 3 Does NOT Do

❌ Recommend pipelines (that's Module 4)  
❌ Perform differential expression  
❌ Quality control (that's Module 2)  
❌ Modify or normalize data

---

## Philosophy & Design

### The Profiling Strategy

RAPTOR Module 3 employs a **feature-based characterization approach** grounded in RNA-seq benchmarking literature.

#### Why Profiling Matters

Different DE pipelines excel under different data conditions:

| Condition | Best Pipeline | Why |
|-----------|--------------|-----|
| High dispersion (BCV > 0.4) | edgeR | Robust to biological variability |
| Low dispersion (BCV < 0.2) | limma-voom | Efficient with homogeneous data |
| Small samples (n < 8) | DESeq2, edgeR | Parametric shrinkage helps |
| Large samples (n ≥ 8) | Wilcoxon, limma-voom | Non-parametric or empirical Bayes |
| Many outliers | edgeR_robust | Downweights outliers |
| Low counts | edgeR | More sensitive detection |

**Module 3 extracts the features that distinguish these conditions.**

### Design Principles

1. **Literature-Grounded**: Every feature is motivated by benchmarking studies
2. **Comprehensive**: 32 features cover all aspects affecting pipeline performance
3. **Interpretable**: Features have clear biological/statistical meaning
4. **ML-Compatible**: Standardized feature vector for recommendation models
5. **Fast Computation**: Optimized for large datasets
6. **Robust Estimation**: Handles edge cases and missing data

### The Three-Stage Flow

```
Module 2: QC
    ↓
    Quality passed → Module 3: Profile
    ↓
    Extract 32 features → Module 4: Recommend
    ↓
    Select optimal pipeline
```

---

## Quick Start

### Minimal Working Example

```bash
# Profile data after QC
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# With metadata for group-based features
raptor profile \
    --counts results/quick_counts/quick_gene_counts.csv \
    --metadata sample_metadata.csv \
    --group-column condition

# Custom output location
raptor profile \
    --counts counts.csv \
    --metadata metadata.csv \
    --output my_profile.json
```

### What Happens

1. **Validation** (1 second): Checks count matrix format
2. **Sample Profiling** (2 seconds): Group sizes, balance
3. **Library Profiling** (2 seconds): Size distribution
4. **Expression Profiling** (5 seconds): Detection rates
5. **Dispersion Estimation** (15 seconds): **CRITICAL** - BCV calculation
6. **Sparsity Analysis** (3 seconds): Zero patterns
7. **Count Distribution** (3 seconds): Low/medium/high counts
8. **Mean-Variance** (5 seconds): Relationship fitting
9. **Feature Vector** (1 second): Build 32-feature vector

**Total**: ~37 seconds for typical dataset

### Output Location

```
results/profile/
├── data_profile.json          # Complete profile with all 32 features
├── profile_summary.txt        # Human-readable summary
└── feature_vector.csv         # Feature vector for ML (optional)
```

---

## The 32 Features

Module 3 extracts **32 numerical features** organized into 8 categories:

### Feature Summary Table

| Category | # Features | Purpose | Examples |
|----------|-----------|---------|----------|
| **Sample Characteristics** | 5 | Sample size, balance, groups | n_samples, min_group_size |
| **Library Size** | 4 | Sequencing depth variation | library_size_cv, range |
| **Gene Detection** | 3 | Expression coverage | detection_rate, reliable_rate |
| **Expression Distribution** | 4 | Overall expression patterns | mean, variance, skewness |
| **Dispersion** | 5 | **CRITICAL** - Biological variation | BCV, common_dispersion |
| **Sparsity** | 3 | Zero counts and dropout | sparsity, zero_inflation |
| **Count Distribution** | 4 | Count magnitude patterns | low/medium/high proportions |
| **Mean-Variance** | 3 | Technical model fit | slope, R², Poisson fit |
| **Quality** | 1 | Overall data quality | quality_score |

**Total**: 32 features

---

## Scientific Basis

### Key Literature

Module 3 features are based on these benchmarking studies:

#### 1. **Dispersion & BCV** (MOST CRITICAL)

**Robinson et al. (2010)** - *Bioinformatics*  
"edgeR: a Bioconductor package for differential expression analysis"
- Introduced biological coefficient of variation (BCV)
- BCV = sqrt(dispersion) represents biological variability
- **High BCV (>0.4)**: edgeR handles well
- **Low BCV (<0.2)**: limma-voom more efficient

**Love et al. (2014)** - *Genome Biology*  
"Moderated estimation of fold change and dispersion for RNA-seq data"
- DESeq2 uses dispersion shrinkage
- Performance depends on dispersion structure
- **Stable dispersion**: DESeq2 excels
- **Variable dispersion**: edgeR more robust

#### 2. **Sample Size Thresholds**

**Li et al. (2022)** - *Genome Biology*  
"A comprehensive evaluation of differential gene expression methods"
- **n < 8 per group**: Parametric methods (DESeq2, edgeR) better
- **n ≥ 8 per group**: Wilcoxon test has better FDR control
- **n > 20 per group**: limma-voom or Wilcoxon recommended

**Soneson & Delorenzi (2013)** - *BMC Bioinformatics*  
"A comparison of methods for differential expression analysis"
- DESeq2/edgeR need ≥3 replicates for reliable estimates
- **n = 2**: Insufficient for variance estimation
- **n = 3**: Minimum for parametric methods
- **n ≥ 5**: Good performance across methods

#### 3. **Library Size Variation**

**Evans et al. (2018)** - *BMC Genomics*  
"Selecting normalization methods for RNA-seq"
- High CV in library sizes → normalization critical
- **CV > 0.5**: Use CPM or TMM
- **Range > 5x**: Careful normalization needed
- DESeq2 median-of-ratios vs edgeR TMM performance

#### 4. **Count Distribution**

**Law et al. (2014)** - *Genome Biology*  
"voom: precision weights unlock linear model analysis tools"
- **Low counts (<10)**: edgeR more sensitive
- **Medium/high counts**: All methods comparable
- Zero-inflation: Negative binomial handles better

#### 5. **Outlier Sensitivity**

**Zhou et al. (2014)** - *Genome Biology*  
"Robustly detecting differential expression in RNA sequencing data"
- Small samples + outliers: edgeR_robust
- Large samples + outliers: limma-voom with robust weights
- Cook's distance in DESeq2 for outlier detection

### Mathematical Models

#### Negative Binomial Distribution

RNA-seq counts follow:
```
Y ~ NB(μ, φ)
E(Y) = μ
Var(Y) = μ + φμ²
```

Where:
- `μ` = mean count
- `φ` = dispersion parameter
- `BCV = sqrt(φ)` = biological coefficient of variation

**Interpretation**:
- φ = 0 → Poisson (no overdispersion)
- φ = 0.01 → BCV = 10% (low biological variation)
- φ = 0.04 → BCV = 20% (moderate)
- φ = 0.16 → BCV = 40% (high)

---

## Usage Examples

### Example 1: Basic Profiling

```bash
# Minimal usage
raptor profile --counts counts.csv
```

**Output**:
```
╔════════════════════════════════════════════════════════════════╗
║              🦖 RAPTOR DATA PROFILE                            ║
╠════════════════════════════════════════════════════════════════╣

  📊 DIMENSIONS
  ────────────────────────────────────────
    Genes:        20,145
    Samples:           60
    Groups:             2
    Balance:         1.00

  📚 LIBRARY SIZE
  ────────────────────────────────────────
    Mean:      15,234,567 reads
    Median:    14,987,342 reads
    CV:             0.234
    Range:          2.3x
    Category:   moderate

  🧬 EXPRESSION
  ────────────────────────────────────────
    Detected genes:     16,234 (80.6%)
    Reliable (≥10):     12,456 (61.8%)
    Mean (log2):          6.32
    Sparsity:           45.2%

  📈 DISPERSION (Critical for Pipeline Selection)
  ────────────────────────────────────────
    Common φ:          0.0289
    BCV:               0.170 (17.0% biological variation)
    BCV Category:  moderate
    Overdispersion:    2.34x

  📊 COUNT DISTRIBUTION
  ────────────────────────────────────────
    Zero:              45.2%
    Low (1-9):         22.3%
    Medium (10-99):    18.7%
    High (100-999):    10.2%
    Very High (≥1000):  3.6%
    Dynamic Range:     14.2 log2 units

  📐 MEAN-VARIANCE RELATIONSHIP
  ────────────────────────────────────────
    Slope:             1.623 (Poisson=1, NB≈1.5-2)
    R²:                0.892

╚════════════════════════════════════════════════════════════════╝
```

### Example 2: With Metadata (Group-Based Features)

```bash
raptor profile \
    --counts counts.csv \
    --metadata metadata.csv \
    --group-column treatment
```

**Metadata Format** (`metadata.csv`):
```csv
sample_id,treatment,batch
Sample1,Control,Batch1
Sample2,Control,Batch1
Sample3,Control,Batch2
Sample4,Drug,Batch1
Sample5,Drug,Batch2
```

**Benefits**:
- Calculates group-specific metrics
- Assesses sample balance
- Enables biological signal estimation
- Better dispersion estimation with groups

### Example 3: Python API Usage

```python
from raptor.profiler import RNAseqDataProfiler, profile_data_quick
import pandas as pd

# Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Quick profiling
profile = profile_data_quick(counts, metadata, group_column='condition')

# Print summary
print(profile.summary())

# Access specific features
print(f"BCV: {profile.bcv:.3f}")
print(f"Dispersion category: {profile.bcv_category}")
print(f"Min group size: {profile.min_group_size}")
print(f"Sparsity: {profile.sparsity:.1%}")

# Get recommendation-relevant features
rec_features = profile.get_recommendation_features()
print(f"Has sufficient replicates: {rec_features['has_sufficient_replicates']}")
print(f"Overdispersion ratio: {rec_features['overdispersion_ratio']:.2f}")

# Access ML feature vector
feature_vector = profile.feature_vector
feature_names = profile.feature_names
print(f"\nFeature vector shape: {feature_vector.shape}")
print(f"Number of features: {len(feature_names)}")

# Save profile
with open('profile.json', 'w') as f:
    f.write(profile.to_json())
```

### Example 4: Advanced - Custom Profiler

```python
from raptor.profiler import RNAseqDataProfiler

# Create profiler with custom settings
profiler = RNAseqDataProfiler(
    counts=counts,
    metadata=metadata,
    group_column='condition',
    min_count_threshold=5  # Higher threshold for "expressed"
)

# Run profiling
profile = profiler.run_full_profile()

# Access detailed metrics
print(f"Library sizes: {profile.library_sizes}")
print(f"Group sizes: {profile.group_sizes}")
print(f"Design complexity: {profile.design_complexity}")

# Dispersion details
print(f"\nDispersion Analysis:")
print(f"  Common: {profile.common_dispersion:.4f}")
print(f"  Mean: {profile.dispersion_mean:.4f}")
print(f"  Median: {profile.dispersion_median:.4f}")
print(f"  IQR: {profile.dispersion_iqr:.4f}")

# Mean-variance relationship
print(f"\nMean-Variance Relationship:")
print(f"  Slope: {profile.mean_var_slope:.3f}")
print(f"  Intercept: {profile.mean_var_intercept:.3f}")
print(f"  R²: {profile.mean_var_r_squared:.3f}")
print(f"  Poisson fit: {profile.poisson_fit_score:.3f}")
```

### Example 5: Quick Characteristics

```python
from raptor.profiler import get_key_characteristics

# Fast assessment without full profiling
chars = get_key_characteristics(counts)

print(f"Samples: {chars['n_samples']}")
print(f"Genes: {chars['n_genes']}")
print(f"Mean library size: {chars['mean_library_size']:,.0f}")
print(f"Sparsity: {chars['sparsity']:.1%}")
print(f"Estimated BCV: {chars['bcv']:.3f}")
```

---

## Interpreting Profiles

### Critical Features for Pipeline Selection

#### 1. BCV (Biological Coefficient of Variation)

**Most Important Feature for Pipeline Selection**

```python
BCV = sqrt(common_dispersion)
```

**Interpretation**:

| BCV | Category | Biological Variation | Best Pipeline |
|-----|----------|---------------------|---------------|
| < 0.1 | Very low | Highly reproducible | limma-voom |
| 0.1-0.2 | Low | Homogeneous samples | limma-voom, DESeq2 |
| 0.2-0.4 | Moderate | Typical RNA-seq | DESeq2, edgeR |
| 0.4-0.6 | High | Variable biology | edgeR, edgeR_robust |
| > 0.6 | Very high | Heterogeneous | edgeR_robust |

**Examples**:
- **Cell lines, controlled conditions**: BCV = 0.10-0.15
- **Inbred mice, same age/sex**: BCV = 0.15-0.25
- **Human samples, matched**: BCV = 0.25-0.35
- **Human samples, diverse**: BCV = 0.35-0.50
- **Clinical samples, heterogeneous**: BCV > 0.50

**What Affects BCV**:
- ✅ Biological variability (genetics, environment)
- ✅ Sample heterogeneity
- ✅ Experimental control
- ❌ NOT sequencing depth
- ❌ NOT library size

#### 2. Sample Size (Min Group Size)

**Second Most Important**

| Min Group Size | Recommendation | Why |
|----------------|----------------|-----|
| < 3 | **STOP** | Insufficient for variance estimation |
| 3-4 | DESeq2, edgeR | Parametric shrinkage helps |
| 5-7 | DESeq2, edgeR, limma-voom | Good for parametric |
| 8-12 | DESeq2, edgeR, limma, Wilcoxon | All methods work |
| 13-20 | limma-voom, Wilcoxon | Non-parametric competitive |
| > 20 | limma-voom, Wilcoxon | Non-parametric preferred |

**Key Threshold**: n = 8 per group
- Below: Use parametric (DESeq2, edgeR)
- Above: Can use non-parametric (Wilcoxon)

#### 3. Library Size CV

**Affects Normalization Choice**

| CV | Category | Action |
|----|----------|--------|
| < 0.2 | Low | Simple normalization OK |
| 0.2-0.5 | Moderate | Standard (CPM, TMM, RLE) |
| 0.5-1.0 | High | Critical - use CPM/TMM |
| > 1.0 | Very high | Check for failed libraries |

**Example**:
```
Library sizes: 5M, 8M, 12M, 15M
Mean: 10M
CV: 0.42 (moderate)
→ Standard normalization sufficient
```

#### 4. Sparsity

**Affects Low-Count Gene Detection**

| Sparsity | Interpretation | Impact |
|----------|----------------|--------|
| < 30% | Low | Easy to analyze |
| 30-60% | Moderate | Typical RNA-seq |
| 60-80% | High | Need sensitive methods |
| > 80% | Very high | Single-cell or targeted |

**High sparsity (>60%)**:
- edgeR more sensitive to low counts
- Consider filtering very lowly expressed genes
- Zero-inflation may be issue

#### 5. Overdispersion Ratio

**Median(Variance / Mean)**

| Ratio | Interpretation | Model Fit |
|-------|----------------|-----------|
| ≈ 1 | No overdispersion | Poisson adequate |
| 1-3 | Moderate | Negative binomial |
| 3-10 | High | NB, may need robustness |
| > 10 | Very high | Check data quality |

**Typical**: 2-4x for RNA-seq

#### 6. Sample Balance

**Min group size / Max group size**

| Balance | Interpretation | Action |
|---------|----------------|--------|
| 1.0 | Perfect | Optimal power |
| 0.8-0.99 | Good | Minor imbalance OK |
| 0.5-0.79 | Moderate | Consider resampling |
| < 0.5 | Poor | Unbalanced design issues |

**Example**:
```
Group A: 10 samples
Group B: 10 samples
Balance: 10/10 = 1.0 (perfect)

Group A: 5 samples
Group B: 10 samples
Balance: 5/10 = 0.5 (moderate imbalance)
```

### Decision Trees Based on Profile

#### Tree 1: Pipeline Selection by BCV and Sample Size

```
BCV Category?
├─ Low (BCV < 0.2)
│   └─ n ≥ 8? 
│       ├─ Yes → limma-voom or Wilcoxon
│       └─ No → DESeq2 or limma-voom
│
├─ Moderate (BCV 0.2-0.4)
│   └─ n < 8? 
│       ├─ Yes → DESeq2 or edgeR
│       └─ No → DESeq2, edgeR, limma-voom
│
└─ High (BCV > 0.4)
    └─ Outliers present?
        ├─ Yes → edgeR_robust
        └─ No → edgeR
```

#### Tree 2: Normalization by Library Size

```
Library Size CV?
├─ CV < 0.2 → Simple log2 or CPM
├─ CV 0.2-0.5 → CPM or TMM (standard)
├─ CV 0.5-1.0 → TMM or RLE (careful)
└─ CV > 1.0 → Check for failed libraries
```

#### Tree 3: Filtering by Sparsity

```
Sparsity?
├─ < 30% → Minimal filtering (keep mean > 1)
├─ 30-60% → Standard (keep mean > 5-10)
├─ 60-80% → Aggressive (keep mean > 10)
└─ > 80% → Very aggressive or specialized methods
```

---

## Integration with RAPTOR Workflow

### The Complete Pipeline

```
Module 1: Quick Quantification
    ↓
    quick_gene_counts.csv
    ↓
Module 2: Quality Assessment
    ↓
    Quality score, outliers detected
    ↓
    Decision: Quality OK?
    ├─ No → Fix issues, re-run QC
    └─ Yes ↓
Module 3: Data Profiling (THIS MODULE)
    ↓
    32 features extracted → data_profile.json
    ↓
Module 4: Pipeline Recommendation
    ↓
    Optimal pipeline selected based on profile
    ↓
Module 5: Production Quantification (if needed)
    ↓
Modules 6-7: Differential Expression
```

### Handoff to Module 4

**What Module 3 Provides to Module 4**:

```json
{
  "recommendation_features": {
    "bcv": 0.234,
    "bcv_category": "moderate",
    "min_group_size": 5,
    "sample_balance": 0.833,
    "has_sufficient_replicates": true,
    "library_size_cv": 0.312,
    "low_count_proportion": 0.223,
    "sparsity": 0.452,
    "has_outliers": false,
    "has_batch_effect": true,
    "quality_score": 87.3
  },
  "feature_vector": [2.807, 4.332, 2, 5, 0.833, ...]
}
```

**Module 4 Uses These Features To**:
1. Select optimal DE pipeline (DESeq2, edgeR, limma-voom, Wilcoxon)
2. Configure pipeline parameters
3. Choose normalization method
4. Set filtering thresholds
5. Decide on robustness options

### Typical Workflow Commands

```bash
# Step 1: QC (Module 2)
raptor qc --counts results/quick_counts/quick_gene_counts.csv \
          --metadata metadata.csv

# Step 2: Profile (Module 3)
raptor profile --counts results/quick_counts/quick_gene_counts.csv \
               --metadata metadata.csv \
               --group-column treatment

# Step 3: Get Recommendation (Module 4)
raptor recommend --profile results/profile/data_profile.json

# Or in one command (if QC passed)
raptor qc --counts counts.csv --metadata metadata.csv && \
raptor profile --counts counts.csv --metadata metadata.csv && \
raptor recommend --profile results/profile/data_profile.json
```

---

## Feature Categories

### Category 1: Sample Characteristics (5 features)

#### Features:
1. **log2_n_samples**: Log₂ of total samples
2. **log2_n_genes**: Log₂ of total genes  
3. **n_groups**: Number of experimental groups
4. **min_group_size**: Smallest group size
5. **sample_balance**: Min/max group size ratio

#### Why These Matter:
- **Sample size** determines statistical power
- **Group balance** affects false positive rate
- **Number of groups** indicates design complexity
- Small samples need parametric shrinkage (DESeq2, edgeR)
- Large samples can use empirical methods (limma-voom)

#### Example:
```json
{
  "n_samples": 60,
  "n_groups": 2,
  "min_group_size": 28,
  "max_group_size": 32,
  "sample_balance": 0.875  // Good balance
}
```

### Category 2: Library Size (4 features)

#### Features:
1. **library_size_cv**: Coefficient of variation
2. **log2_library_size_range**: Log₂(max/min ratio)
3. **log2_library_size_mean**: Log₂ of mean size
4. **library_size_iqr_norm**: IQR normalized by mean

#### Why These Matter:
- High CV → normalization critical
- Large range → potential failed libraries
- Mean size indicates sequencing depth
- IQR shows distribution spread

#### Example:
```json
{
  "library_size_mean": 15234567,
  "library_size_cv": 0.234,  // Moderate variation
  "library_size_range": 2.3,  // Modest range
  "sequencing_depth_category": "standard"
}
```

### Category 3: Gene Detection (3 features)

#### Features:
1. **detection_rate**: Proportion of genes detected
2. **reliable_detection_rate**: Proportion with mean ≥ 10
3. **genes_per_sample_cv**: Variation in genes detected per sample

#### Why These Matter:
- Detection rate indicates data quality
- Reliable detection affects number of testable genes
- CV shows sample consistency

#### Example:
```json
{
  "n_total_genes": 20145,
  "n_expressed_genes": 16234,  // 80.6%
  "n_reliably_expressed": 12456,  // 61.8%
  "detection_rate": 0.806,
  "reliable_detection_rate": 0.618
}
```

### Category 4: Expression Distribution (4 features)

#### Features:
1. **expression_mean**: Mean log₂ expression
2. **expression_variance**: Variance of log₂ expression
3. **expression_skewness**: Distribution skewness
4. **expression_iqr**: Interquartile range

#### Why These Matter:
- Mean shows overall expression level
- Variance indicates heterogeneity
- Skewness shows distribution shape
- IQR robust to outliers

### Category 5: Dispersion (5 features) ⭐ CRITICAL

#### Features:
1. **bcv**: Biological coefficient of variation (√φ)
2. **common_dispersion**: Shared dispersion estimate (φ)
3. **dispersion_trend_slope**: Mean-dependent dispersion trend
4. **overdispersion_ratio**: Variance/mean ratio (log₂)
5. **dispersion_iqr**: IQR of gene dispersions

#### Why These Matter:
**MOST IMPORTANT FOR PIPELINE SELECTION**

- **BCV** directly determines which pipeline works best
- Common dispersion used by DESeq2/edgeR
- Trend slope shows if dispersion is mean-dependent
- Overdispersion ratio validates NB model
- IQR shows dispersion heterogeneity

#### Mathematical Basis:
```
Var(Y) = μ + φμ²
BCV = √φ

For a gene with mean count = 100:
  BCV = 0.2 → Var = 100 + 0.04×100² = 500
  BCV = 0.4 → Var = 100 + 0.16×100² = 1700
```

**Decision Rules**:
```
BCV < 0.2  → limma-voom (low variation)
BCV 0.2-0.4 → DESeq2, edgeR (moderate)
BCV > 0.4  → edgeR, edgeR_robust (high variation)
```

### Category 6: Sparsity (3 features)

#### Features:
1. **sparsity**: Overall proportion of zeros
2. **zero_inflation_index**: Excess zeros vs Poisson
3. **dropout_rate**: Per-sample average zero rate

#### Why These Matter:
- High sparsity → need sensitive methods (edgeR)
- Zero-inflation indicates model assumptions
- Dropout shows technical vs biological zeros

#### Example:
```json
{
  "sparsity": 0.452,  // 45.2% zeros
  "zero_inflation_index": 0.12,  // Slightly inflated
  "is_zero_inflated": true,
  "sparsity_category": "moderate"
}
```

### Category 7: Count Distribution (4 features)

#### Features:
1. **low_count_proportion**: Proportion of counts 1-9
2. **medium_count_proportion**: Proportion 10-99
3. **high_count_proportion**: Proportion 100-999
4. **dynamic_range**: Log₂(max/min) of nonzero counts

#### Why These Matter:
- Low counts → edgeR more sensitive
- High counts → all methods similar
- Dynamic range → normalization complexity

#### Example:
```json
{
  "zero_proportion": 0.452,
  "low_count_proportion": 0.223,  // 22.3% counts 1-9
  "medium_count_proportion": 0.187,
  "high_count_proportion": 0.102,
  "very_high_count_proportion": 0.036,
  "dynamic_range": 14.2  // log2 units
}
```

### Category 8: Mean-Variance (3 features)

#### Features:
1. **mean_var_slope**: Slope of log(var) ~ log(mean)
2. **mean_var_r_squared**: Fit quality
3. **poisson_fit_score**: How close to Poisson (slope=1)

#### Why These Matter:
- Slope = 1 → Poisson (rare for RNA-seq)
- Slope ≈ 1.5-2 → Negative binomial (typical)
- R² shows if relationship is consistent
- Validates statistical model assumptions

#### Expected Values:
```
Poisson: slope = 1
RNA-seq: slope ≈ 1.5-2.0
Overdispersed: slope > 2
```

---

## Mathematical Formulations

### Dispersion Estimation (Method of Moments)

**For each gene i**:
```
μᵢ = mean count across samples
σᵢ² = variance across samples

Dispersion: φᵢ = (σᵢ² - μᵢ) / μᵢ²
```

**Common dispersion** (trimmed mean):
```
φ_common = trimmed_mean(φ₁, φ₂, ..., φₙ)
          (trim top/bottom 10%)

BCV = √φ_common
```

### Mean-Variance Relationship

**Fit model**:
```
log₂(σᵢ²) = α + β·log₂(μᵢ) + εᵢ
```

**Interpretation**:
- β = 1: Var ∝ μ (Poisson)
- β ≈ 1.5-2: Typical RNA-seq (NB)
- β > 2: High overdispersion

**R² Quality**:
- R² > 0.8: Good fit
- R² < 0.5: Heterogeneous variance

### Zero-Inflation Index

**Expected zeros under Poisson**:
```
E[zeros] = Σᵢ exp(-μᵢ) × n_samples
```

**Observed zeros**:
```
O[zeros] = count(observations == 0)
```

**Zero-inflation index**:
```
ZI = (O - E) / E

ZI ≈ 0: As expected
ZI > 0.1: Zero-inflated
ZI < -0.1: Zero-deflated (rare)
```

### Library Size Metrics

**Coefficient of variation**:
```
CV = σ / μ
   = std(library_sizes) / mean(library_sizes)
```

**Range**:
```
Range = max(library_sizes) / min(library_sizes)
```

**Normalized IQR**:
```
IQR_norm = IQR(library_sizes) / mean(library_sizes)
```

---

## Troubleshooting

### Issue 1: "Too Few Samples for Dispersion"

**Error**:
```
ValidationError: Cannot estimate dispersion with < 2 samples per group
```

**Cause**: Insufficient replication

**Solution**:
```bash
# Check sample sizes
cat metadata.csv | cut -d',' -f2 | sort | uniq -c

# Need ≥2 per group for variance
# Recommended: ≥3 per group
```

**Workaround**:
- Profile without groups: `raptor profile --counts counts.csv`
- Merge groups if biologically acceptable
- Get more samples if possible

### Issue 2: "Extreme BCV Value"

**Symptom**: BCV > 1.0 or BCV = NaN

**Causes**:
1. Very heterogeneous samples
2. Low-quality data
3. Mixed sample types
4. Calculation error

**Diagnosis**:
```python
# Check dispersion distribution
print(f"Dispersion range: {profile.dispersion_mean:.4f}")
print(f"Common dispersion: {profile.common_dispersion:.4f}")
print(f"BCV: {profile.bcv:.3f}")

# Check for outliers
if profile.has_outliers:
    print("Outliers detected - may affect BCV")
```

**Solutions**:
1. **Run QC first**: Remove outliers
2. **Check sample mixing**: Ensure homogeneous groups
3. **Filter genes**: Remove very lowly expressed
4. **Check metadata**: Verify group assignments

### Issue 3: "Sparsity Too High"

**Symptom**: Sparsity > 80%

**Interpretation**:
- Normal for single-cell RNA-seq
- Problematic for bulk RNA-seq

**Solutions**:
```bash
# For bulk RNA-seq with high sparsity:
# 1. Check library sizes
raptor profile --counts counts.csv | grep "Library Size"

# 2. Filter very lowly expressed genes first
# Keep genes with mean count ≥ 5-10

# 3. Check if data is actually single-cell
# If yes, use specialized tools (not RAPTOR)
```

### Issue 4: "Poor Mean-Variance Fit"

**Symptom**: R² < 0.5

**Causes**:
- Heterogeneous variance structure
- Different gene classes (housekeeping vs regulated)
- Batch effects
- Data quality issues

**Diagnosis**:
```python
print(f"Mean-variance R²: {profile.mean_var_r_squared:.3f}")
print(f"Slope: {profile.mean_var_slope:.3f}")
print(f"Has batch effect: {profile.has_batch_effect}")
```

**Solutions**:
1. Run QC (Module 2) to check quality
2. Correct batch effects before profiling
3. Filter extreme outlier genes
4. May indicate need for robust methods

### Issue 5: "Library Size CV Very High"

**Symptom**: CV > 1.0

**Interpretation**:
- Likely failed libraries present
- OR intentionally different depths

**Diagnosis**:
```python
print(f"Library sizes: {profile.library_sizes}")
print(f"Min: {profile.library_size_min:,.0f}")
print(f"Max: {profile.library_size_max:,.0f}")
print(f"Range: {profile.library_size_range:.1f}x")
```

**Solutions**:
```bash
# Identify failed libraries (< 1M reads)
# Remove or re-sequence them
# Then re-run profiling
```

### Issue 6: Memory Error

**Symptom**:
```
MemoryError: Unable to allocate array
```

**Cause**: Very large dataset (100k+ genes, 1000+ samples)

**Solutions**:
```python
# 1. Profile on filtered data
filtered_counts = counts[counts.mean(axis=1) > 5]
profile = profile_data_quick(filtered_counts, metadata)

# 2. Use quick characteristics instead
from raptor.profiler import get_key_characteristics
chars = get_key_characteristics(counts)  # Much faster, less memory
```

---

## Best Practices

### 1. Always Profile After QC

❌ **Don't**:
```bash
# Skip QC and go straight to profiling
raptor profile --counts counts.csv
```

✅ **Do**:
```bash
# QC first to ensure quality
raptor qc --counts counts.csv --metadata metadata.csv

# Then profile
raptor profile --counts counts.csv --metadata metadata.csv
```

**Why**: Outliers and quality issues distort feature estimates, especially BCV.

### 2. Include Metadata

❌ **Don't**:
```bash
raptor profile --counts counts.csv
# Missing group-based features
```

✅ **Do**:
```bash
raptor profile --counts counts.csv --metadata metadata.csv
# Gets sample balance, group sizes, better dispersion
```

### 3. Use Correct Group Column

❌ **Don't**:
```bash
# Using batch as groups
raptor profile --counts counts.csv --metadata metadata.csv --group-column batch
```

✅ **Do**:
```bash
# Use biological groups (condition, treatment)
raptor profile --counts counts.csv --metadata metadata.csv --group-column condition
```

### 4. Check Profile Before Recommendation

✅ **Do**:
```python
profile = profile_data_quick(counts, metadata)
print(profile.summary())

# Sanity checks:
assert profile.min_group_size >= 3, "Need ≥3 replicates"
assert profile.quality_score >= 60, "Quality too low"
assert 0 < profile.bcv < 1.0, "BCV unrealistic"
```

### 5. Save Profiles for Reproducibility

✅ **Do**:
```python
# Save profile
with open('profile_20240309.json', 'w') as f:
    f.write(profile.to_json())

# Version control
git add profile_20240309.json
git commit -m "Data profile for analysis 2024-03-09"
```

### 6. Document Key Features

✅ **Do**:
```markdown
# Analysis Notes

**Data Profile (2024-03-09)**
- BCV: 0.234 (moderate biological variation)
- Sample size: 30 per group (adequate)
- Sparsity: 45% (typical)
- Recommended: DESeq2 or edgeR
```

### 7. Profile Different Subsets

```python
# Profile full dataset
full_profile = profile_data_quick(counts, metadata)

# Profile high-quality subset (after filtering)
expressed = counts[counts.mean(axis=1) >= 10]
filtered_profile = profile_data_quick(expressed, metadata)

# Compare
print(f"Full BCV: {full_profile.bcv:.3f}")
print(f"Filtered BCV: {filtered_profile.bcv:.3f}")
```

---

## Technical Details

### Feature Vector Construction

The 32-feature vector is standardized for ML:

```python
feature_vector = [
    # Sample (5)
    log2(n_samples + 1),           # Compress range
    log2(n_genes + 1),             # Compress range
    n_groups,                      # Raw count
    min_group_size,                # Raw count
    sample_balance,                # [0, 1]
    
    # Library (4)
    library_size_cv,               # [0, ∞)
    log2(library_size_range + 1),  # Compress range
    log2(library_size_mean + 1),   # Compress range
    library_size_iqr_norm,         # Normalized
    
    # Detection (3)
    detection_rate,                # [0, 1]
    reliable_detection_rate,       # [0, 1]
    genes_per_sample_cv,           # [0, ∞)
    
    # Expression (4)
    expression_mean,               # log2 scale
    expression_variance,           # Raw
    expression_skewness,           # [-∞, ∞]
    expression_iqr,                # Raw
    
    # Dispersion (5) - CRITICAL
    bcv,                          # [0, ∞)
    common_dispersion,            # [0, ∞)
    dispersion_trend_slope,       # Raw
    log2(overdispersion_ratio+1), # Compress range
    dispersion_iqr,               # Raw
    
    # Sparsity (3)
    sparsity,                     # [0, 1]
    clip(zero_inflation, -5, 5),  # Bounded
    dropout_rate,                 # [0, 1]
    
    # Counts (4)
    low_count_proportion,         # [0, 1]
    medium_count_proportion,      # [0, 1]
    high_count_proportion,        # [0, 1]
    log2(dynamic_range + 1),      # Compress range
    
    # Mean-var (3)
    mean_var_slope,               # Typically 1-2
    mean_var_r_squared,           # [0, 1]
    poisson_fit_score,            # [0, 1]
    
    # Quality (1)
    quality_score / 100.0         # Normalize [0, 1]
]
```

**Transformations Applied**:
- **log2(x + 1)**: For count-like features (compress large values)
- **Normalization**: Divide by mean/total to get proportions
- **Clipping**: Bound outliers to reasonable ranges
- **NaN handling**: Replace with 0.0 or appropriate default

### Computation Complexity

| Operation | Time Complexity | Space | Notes |
|-----------|----------------|-------|-------|
| Library sizes | O(n × m) | O(n) | Fast |
| Gene detection | O(n × m) | O(n) | Fast |
| Expression dist | O(n × m) | O(n × m) | Moderate |
| **Dispersion** | **O(n × m²)** | **O(n)** | **Slowest** |
| Sparsity | O(n × m) | O(1) | Fast |
| Mean-variance | O(n) | O(n) | Fast |

Where:
- n = number of genes
- m = number of samples

**Bottleneck**: Dispersion estimation (requires variance calculation per gene)

**Optimization**:
- Uses NumPy vectorization
- Trims extremes (10% each end)
- Skips very lowly expressed genes (mean < 1)

### Performance Benchmarks

| Dataset Size | Time | Memory | Notes |
|--------------|------|--------|-------|
| 10k genes × 10 samples | 3 sec | < 100 MB | Tiny |
| 20k genes × 50 samples | 30 sec | 500 MB | Typical |
| 50k genes × 100 samples | 2 min | 2 GB | Large |
| 100k genes × 500 samples | 15 min | 10 GB | Very large |

**Tested On**: MacBook Pro, 16GB RAM, M1 chip

---

## Summary

### Key Takeaways

1. **Module 3 extracts 32 features** that characterize RNA-seq data across 8 dimensions

2. **BCV is the most critical feature** for pipeline selection:
   - Low BCV → limma-voom
   - Moderate BCV → DESeq2, edgeR
   - High BCV → edgeR, edgeR_robust

3. **Sample size matters**:
   - n < 8: Parametric methods (DESeq2, edgeR)
   - n ≥ 8: Non-parametric competitive (Wilcoxon)

4. **Always run QC (Module 2) before profiling** to ensure clean feature estimates

5. **Profile output feeds Module 4** for ML-based pipeline recommendation

6. **Features are literature-grounded** from RNA-seq benchmarking studies

### Workflow Checklist

- [ ] Run Module 2 QC first
- [ ] Quality score ≥ 60
- [ ] Remove outliers if detected
- [ ] Run Module 3 profiling with metadata
- [ ] Check profile summary
- [ ] Verify BCV is reasonable (0.1-0.6)
- [ ] Confirm sample sizes adequate (≥3 per group)
- [ ] Save profile JSON
- [ ] Proceed to Module 4 recommendation

---

**Document Version**: 2.2.0  
**Last Updated**: March 2026  
**Author**: Ayeh Bolouki  
**License**: MIT
