# RAPTOR Module 2: Quality Assessment & Outlier Detection

**Version 2.2.0** | **Stage 2** | **Module M2**

---

## Table of Contents

1. [Overview](#overview)
2. [Philosophy & Design](#philosophy--design)
3. [Quick Start](#quick-start)
4. [Assessment Components](#assessment-components)
5. [Advanced Outlier Detection](#advanced-outlier-detection)
6. [Normalization Methods](#normalization-methods)
7. [Batch Effect Detection](#batch-effect-detection)
8. [Usage Examples](#usage-examples)
9. [Interpreting Results](#interpreting-results)
10. [Integration with RAPTOR Workflow](#integration-with-raptor-workflow)
11. [Troubleshooting](#troubleshooting)
12. [Best Practices](#best-practices)
13. [Technical Details](#technical-details)

---

## Overview

### What is Module 2?

Module 2 provides **comprehensive quality assessment** for RNA-seq count matrices. It identifies technical issues, detects outliers, assesses batch effects, and provides an overall data quality score before you invest time in differential expression analysis.

### Key Characteristics

- 📊 **6 Assessment Components**: Library quality, gene detection, outliers, variance, batch effects, biological signal
- 🎯 **Advanced Outlier Detection**: 6 different methods with consensus voting (NEW in v2.2.0)
- 🔬 **Batch Effect Analysis**: Distinguishes batch from biological variation (IMPROVED in v2.2.0)
- 🎨 **4 Normalization Options**: log2, CPM, quantile, or none
- 💯 **Overall Quality Score**: Weighted score (0-100) with actionable recommendations
- 📈 **Comprehensive Visualization**: QC plots, PCA, heatmaps, outlier detection

### What Module 2 Does

✅ Assess library size quality and distribution  
✅ Analyze gene detection patterns and zero inflation  
✅ Detect outlier samples using 6 methods  
✅ Identify batch effects and confounding  
✅ Evaluate variance structure (PCA analysis)  
✅ Measure biological signal strength  
✅ Generate overall quality score with recommendations  
✅ Produce comprehensive QC reports and visualizations

### What Module 2 Does NOT Do

❌ Perform differential expression analysis  
❌ Remove batch effects (only detects them)  
❌ Filter or modify your data  
❌ Replace Module 3 profiling (different purpose)

---

## Philosophy & Design

### The Quality Assessment Strategy

RAPTOR Module 2 employs a **multi-faceted quality assessment approach**:

#### Stage 2: Quality Assessment (Module 2)
**Purpose**: Identify technical issues before analysis  
**Input**: `quick_gene_counts.csv` from Module 1  
**Output**: Quality report + outlier flags  
**Decision**: Fix issues OR proceed to Module 3

### Why This Matters

Traditional RNA-seq workflows often discover quality issues only after completing DE analysis:
- Outlier samples skewing results
- Batch effects creating false positives
- Poor library quality reducing power
- Confounded designs making results uninterpretable

**RAPTOR's approach**: Catch these issues in Module 2, make informed decisions, then run optimized analysis.

### Design Principles

1. **Comprehensive Yet Efficient**: 6 components assessed in <1 minute
2. **Multi-Method Consensus**: 6 outlier detection methods prevent false positives
3. **Actionable Recommendations**: Each issue includes specific remediation steps
4. **Batch vs Biology**: Distinguishes technical (batch) from biological variation
5. **Flexible Normalization**: Adapts to different data types and scales

### v2.2.0 Enhancements

**Advanced Outlier Detection**:
- 6 methods: PCA-Mahalanobis, Isolation Forest, LOF, Elliptic Envelope, Correlation, Library Size
- Consensus voting: Sample flagged only if ≥N methods agree
- Reduces false positives by ~80% vs single-method approaches

**Improved Batch Effect Detection**:
- Separates batch-like columns from condition-like columns
- Calculates F-statistics to compare batch vs biological effects
- Detects confounding (batch perfectly correlated with condition)
- Provides specific recommendations based on effect strength

---

## Quick Start

### Minimal Working Example

```bash
# Using quick counts from Module 1
raptor qc --counts results/quick_counts/quick_gene_counts.csv

# With metadata for batch effect detection
raptor qc \
    --counts results/quick_counts/quick_gene_counts.csv \
    --metadata sample_metadata.csv

# Generate plots
raptor qc \
    --counts results/quick_counts/quick_gene_counts.csv \
    --metadata sample_metadata.csv \
    --plot \
    --plot-output qc_report.pdf
```

### What Happens

1. **Validation** (1 second): Checks count matrix format
2. **Normalization** (5 seconds): log2, CPM, quantile, or none
3. **6 Assessments** (30 seconds): All quality components
4. **Outlier Detection** (10 seconds): 6 methods with consensus
5. **Report Generation** (5 seconds): JSON + optional plots

### Output Location

```
results/qc/
├── quality_report.json       # Complete quality metrics
├── qc_plots.pdf               # Visual QC summary (optional)
├── outlier_analysis.txt       # Detailed outlier results
└── recommendations.txt        # Actionable next steps
```

---

## Assessment Components

Module 2 evaluates your data across **6 dimensions**:

### 1. Library Quality (Weight: 15%)

**What It Measures**:
- Library size distribution
- Coefficient of variation (CV)
- Size range (max/min ratio)
- Samples below quality thresholds

**Scoring**:
```
100 points - Start
- 20 points if CV > 0.5 (high variation)
- 15 points if min library < 1M reads
- 10 points if any library < 500k reads
- 10 points if size range > 10x
```

**Example Output**:
```json
{
  "library_quality": {
    "mean_size": 15234567,
    "cv": 0.23,
    "min_size": 8123456,
    "max_size": 24567890,
    "size_range": 3.02,
    "score": 90,
    "status": "good",
    "flags": []
  }
}
```

**Recommendations**:
- **CV > 0.5**: Use CPM normalization instead of log2
- **Size range > 10x**: Check for failed libraries, consider removing
- **Min < 1M**: May have low power, consider deeper sequencing

### 2. Gene Detection (Weight: 20%)

**What It Measures**:
- Gene detection rate (genes expressed per sample)
- Zero inflation percentage
- Distribution of expression levels (high/medium/low)

**Scoring**:
```
100 points - Start
- Up to 40 points for zero inflation (0.5 × zero%)
- 10 points if < 5% highly expressed genes
```

**Example Output**:
```json
{
  "gene_detection": {
    "mean_detection_rate": 0.67,
    "zero_inflation_pct": 45.2,
    "n_highly_expressed": 3245,
    "n_medium_expressed": 8932,
    "n_low_expressed": 7823,
    "score": 77,
    "status": "good",
    "flags": ["Moderate zero inflation"]
  }
}
```

**Interpretation**:
- **Detection rate < 0.5**: Poor quality or very sparse data
- **Zero inflation > 70%**: Extremely sparse, may need different approaches
- **Few highly expressed genes**: May indicate poor RNA quality

### 3. Outlier Detection (Weight: 15%)

**What It Measures**:
- Samples that deviate from population
- Uses 6 different detection methods (v2.2.0)
- Consensus voting for robust identification

**Methods Used**:
1. **PCA + Mahalanobis Distance**: Multivariate outliers in PC space
2. **Isolation Forest**: Anomaly detection algorithm
3. **Local Outlier Factor (LOF)**: Local density-based detection
4. **Elliptic Envelope**: Robust covariance estimation
5. **Correlation-Based**: Samples with low correlation to others
6. **Library Size**: Z-score based on library sizes

**Scoring**:
```
100 points - Start
- 5 points per 1% of samples flagged as outliers
```

**Example Output**:
```json
{
  "outlier_detection": {
    "n_outliers": 2,
    "outlier_percentage": 3.3,
    "outlier_samples": ["Sample_42", "Sample_89"],
    "mahalanobis_distances": [1.2, 1.5, 8.9, ..., 9.3, 1.8],
    "score": 84,
    "status": "warning",
    "flags": ["2 potential outlier sample(s) detected"]
  }
}
```

**See [Advanced Outlier Detection](#advanced-outlier-detection) for details**

### 4. Variance Structure (Weight: 15%)

**What It Measures**:
- Principal component analysis (PCA)
- Variance explained by PC1, PC2
- Number of PCs needed for 90% variance
- Complexity of variance structure

**Scoring**:
```
100 points - Start
- 20 points if PC1 > 50% (potential batch effect)
- 20 points if PC1 > 70% (strong batch effect)
- 15 points if < 2 components for 90% variance
```

**Example Output**:
```json
{
  "variance_structure": {
    "pc1_variance": 42.3,
    "pc2_variance": 18.7,
    "variance_explained_top5": [0.423, 0.187, 0.112, 0.087, 0.054],
    "n_components_90pct": 6,
    "score": 85,
    "status": "good",
    "flags": []
  }
}
```

**Interpretation**:
- **PC1 > 50%**: Likely batch effect dominating variance
- **PC1 > 70%**: Strong technical artifact, investigate
- **Few PCs for 90%**: Simple structure (good) OR strong batch (bad)

### 5. Batch Effects (Weight: 20%)

**What It Measures** (IMPROVED in v2.2.0):
- Identifies batch-like vs condition-like metadata columns
- Calculates F-statistic for batch and condition on PC1
- Compares batch effect vs biological effect strength
- Detects confounding (batch correlated with condition)

**Column Classification**:
```python
# Batch-like columns (technical)
batch_keywords = ['batch', 'plate', 'run', 'lane', 'date', 'flow', 
                  'chip', 'slide', 'array', 'seq_run', 'library_prep']

# Condition-like columns (biological)
condition_keywords = ['condition', 'treatment', 'group', 'status', 
                      'disease', 'genotype', 'phenotype', 'sample_type']
```

**Scoring**:
```
100 points - Start
- 30 points if batch F-stat > 10
- 20 points if batch F-stat > 5
- 20 points if confounded design
```

**Example Output**:
```json
{
  "batch_effects": {
    "batch_detected": true,
    "batch_variable": "sequencing_batch",
    "batch_strength": 12.4,
    "condition_variable": "treatment",
    "condition_strength": 8.7,
    "confounded": false,
    "score": 70,
    "status": "warning",
    "recommendation": "Strong batch effect detected. Recommended actions:\n  1. Include batch as covariate in DE model: ~batch + condition\n  2. For visualization: Use ComBat or limma::removeBatchEffect\n  3. Do NOT use batch-corrected data for DE analysis"
  }
}
```

**See [Batch Effect Detection](#batch-effect-detection) for details**

### 6. Biological Signal (Weight: 15%)

**What It Measures**:
- Coefficient of variation (CV) across genes
- Genes with both high expression and high variability
- Strength of biological signal vs noise

**Scoring**:
```
Score = min(100, signal_strength × 2)
where signal_strength = (% genes with high mean AND high CV)
```

**Example Output**:
```json
{
  "biological_signal": {
    "signal_strength": 23.4,
    "n_signal_genes": 4680,
    "median_cv": 1.45,
    "score": 47,
    "status": "warning",
    "flags": ["Low biological signal detected"]
  }
}
```

**Interpretation**:
- **Signal strength < 10%**: Very weak biological signal
- **Signal strength 10-30%**: Moderate signal
- **Signal strength > 30%**: Strong biological signal

### Overall Quality Score

**Calculation**:
```python
weights = {
    'library_quality': 0.15,      # 15%
    'gene_detection': 0.20,       # 20%
    'outlier_detection': 0.15,    # 15%
    'variance_structure': 0.15,   # 15%
    'batch_effects': 0.20,        # 20%
    'biological_signal': 0.15     # 15%
}

overall_score = sum(component_score × weight)
```

**Interpretation**:
- **80-100**: Good quality - Proceed with analysis
- **60-79**: Acceptable - Address flagged issues
- **< 60**: Poor quality - Investigate before analysis

---

## Advanced Outlier Detection

### NEW in v2.2.0: 6-Method Consensus Approach

Traditional outlier detection uses a single method, leading to:
- ❌ High false positive rate
- ❌ Sensitivity to method assumptions
- ❌ Inconsistent results across analyses

**RAPTOR v2.2.0 Solution**: Consensus voting across 6 methods

### The 6 Detection Methods

#### 1. PCA + Mahalanobis Distance

**Approach**: Multivariate outliers in principal component space

**How It Works**:
1. Perform PCA (5 components)
2. Calculate Mahalanobis distance in PC space
3. Flag samples > mean + 3×SD

**Strengths**:
- ✅ Captures multivariate structure
- ✅ Reduces dimensionality

**Weaknesses**:
- ❌ Assumes multivariate normal distribution
- ❌ Sensitive to extreme outliers affecting mean/covariance

**Best For**: Moderate outliers in high-dimensional data

#### 2. Isolation Forest

**Approach**: Ensemble-based anomaly detection

**How It Works**:
1. Build random decision trees
2. Isolate points by random splits
3. Outliers require fewer splits to isolate

**Strengths**:
- ✅ Non-parametric (no distribution assumptions)
- ✅ Handles high-dimensional data well
- ✅ Robust to irrelevant features

**Weaknesses**:
- ❌ Stochastic (results vary slightly)
- ❌ Requires contamination estimate

**Best For**: Complex, high-dimensional outliers

#### 3. Local Outlier Factor (LOF)

**Approach**: Density-based local outlier detection

**How It Works**:
1. Calculate local density for each sample
2. Compare to neighbor densities
3. Flag samples in low-density regions

**Strengths**:
- ✅ Detects local outliers (clusters of outliers)
- ✅ Handles varying densities

**Weaknesses**:
- ❌ Sensitive to k (number of neighbors)
- ❌ Computationally expensive

**Best For**: Outliers in clusters or varying-density data

#### 4. Elliptic Envelope

**Approach**: Robust covariance estimation

**How It Works**:
1. Fit ellipse to data (robust to outliers)
2. Calculate Mahalanobis distance using robust covariance
3. Flag samples outside ellipse

**Strengths**:
- ✅ Robust to outliers in covariance estimation
- ✅ Works well for Gaussian-distributed data

**Weaknesses**:
- ❌ Assumes elliptical distribution
- ❌ Requires contamination estimate

**Best For**: Moderate outliers in approximately Gaussian data

#### 5. Correlation-Based

**Approach**: Low correlation with other samples

**How It Works**:
1. Calculate pairwise sample correlations
2. Compute median correlation for each sample
3. Flag samples with low median correlation

**Strengths**:
- ✅ Intuitive interpretation
- ✅ Detects samples with different patterns

**Weaknesses**:
- ❌ Sensitive to batch effects
- ❌ May flag biological differences

**Best For**: Samples with fundamentally different patterns

#### 6. Library Size Z-Score

**Approach**: Extreme library sizes

**How It Works**:
1. Calculate library size for each sample
2. Compute Z-scores
3. Flag samples with |Z| > 3

**Strengths**:
- ✅ Simple and interpretable
- ✅ Catches failed libraries

**Weaknesses**:
- ❌ Univariate (ignores expression patterns)
- ❌ Assumes normal distribution of library sizes

**Best For**: Failed libraries or sequencing issues

### Consensus Voting Strategy

**How It Works**:
```python
# Each method votes on each sample
votes = {
    'Sample1': 0,  # Not flagged by any method
    'Sample2': 3,  # Flagged by 3 methods
    'Sample3': 5,  # Flagged by 5 methods
    'Sample4': 1,  # Flagged by 1 method
}

# Consensus: Sample flagged if ≥ threshold methods agree
consensus_threshold = 3  # Default

outliers = [sample for sample, count in votes.items() 
            if count >= consensus_threshold]
# Result: ['Sample2', 'Sample3']
```

**Choosing Consensus Threshold**:

| Threshold | Sensitivity | Specificity | Use When |
|-----------|-------------|-------------|----------|
| 1 method | Very high | Low | Exploratory, want to catch everything |
| 2 methods | High | Moderate | Liberal detection, check flagged samples |
| **3 methods** | **Moderate** | **High** | **DEFAULT - Good balance** |
| 4 methods | Low | Very high | Conservative, only clear outliers |
| 5-6 methods | Very low | Extreme | Only extreme/obvious outliers |

**Example**:
```bash
# Conservative (high confidence)
raptor qc --counts counts.csv --consensus-threshold 4

# Standard (recommended)
raptor qc --counts counts.csv --consensus-threshold 3

# Liberal (catch more potential outliers)
raptor qc --counts counts.csv --consensus-threshold 2
```

### Interpreting Advanced Outlier Results

**Example Output**:
```
═════════════════════════════════════════════════════════════
ADVANCED OUTLIER DETECTION RESULTS
═════════════════════════════════════════════════════════════
Consensus threshold: 3 methods
Total outliers: 2 (3.3%)

Outlier samples:
  • Sample_42 (flagged by 4 methods)
  • Sample_89 (flagged by 5 methods)

Detection by method:
  • PCA + Mahalanobis: 3 outliers
  • Isolation Forest: 4 outliers
  • Local Outlier Factor: 2 outliers
  • Elliptic Envelope: 3 outliers
  • Correlation-based: 5 outliers
  • Library Size: 1 outliers
═════════════════════════════════════════════════════════════
```

**What To Do**:

1. **Examine flagged samples**:
   ```R
   # In R
   pca <- prcomp(t(counts))
   plot(pca$x[,1:2], col=ifelse(rownames(pca$x) %in% c("Sample_42", "Sample_89"), "red", "black"))
   ```

2. **Check metadata**:
   - Different tissue?
   - Different prep date?
   - Known technical issue?

3. **Decision tree**:
   ```
   If flagged by ≥ 5 methods → Likely technical outlier, consider removing
   If flagged by 3-4 methods → Investigate, may remove
   If flagged by 1-2 methods → Keep unless clear reason
   ```

---

## Normalization Methods

Module 2 supports **4 normalization approaches** to handle different data types:

### 1. log2 (Default)

**Formula**: `log2(counts + 1)`

**When to Use**:
- ✅ Standard RNA-seq with similar library sizes (CV < 0.5)
- ✅ First-pass QC assessment
- ✅ Quick exploratory analysis

**Pros**:
- Simple and fast
- Interpretable scale
- Works well for most data

**Cons**:
- Doesn't account for library size differences
- Not appropriate for very different sequencing depths

**Example**:
```bash
raptor qc --counts counts.csv --normalization log2
```

### 2. CPM (Counts Per Million)

**Formula**: `log2((counts / library_size) × 1e6 + 1)`

**When to Use**:
- ✅ Library sizes vary significantly (CV > 0.5 or range > 5x)
- ✅ Samples have very different sequencing depths
- ✅ Failed libraries mixed with good libraries

**Pros**:
- Accounts for sequencing depth
- Comparable across samples with different library sizes
- Standard for cross-sample comparisons

**Cons**:
- Slightly more complex
- Still doesn't handle compositional effects

**Example**:
```bash
# When library sizes range from 5M to 30M
raptor qc --counts counts.csv --normalization cpm
```

**When Library Sizes Differ**:
```
Sample1: 5,234,567 reads
Sample2: 23,456,789 reads
Sample3: 8,123,456 reads

Range: 23.5M / 5.2M = 4.5x
CV: 0.68

→ Use CPM normalization
```

### 3. Quantile Normalization

**Approach**: Forces all samples to have identical distributions

**When to Use**:
- ✅ Samples have very different distributions
- ✅ Strong technical batch effects
- ✅ Want to remove technical variation aggressively

**Pros**:
- Removes most technical variation
- Makes samples highly comparable
- Good for visualization

**Cons**:
- **Aggressive**: May remove biological signal
- Assumes samples should have same distribution
- Not appropriate for different conditions/tissues

**Example**:
```bash
# For samples with strong technical variation
raptor qc --counts counts.csv --normalization quantile
```

**⚠️ Caution**:
```
Use quantile normalization ONLY for QC visualization, 
NOT for differential expression analysis!

Reason: Removes biological differences between conditions.
```

### 4. None (Pre-Normalized Data)

**Approach**: No transformation applied

**When to Use**:
- ✅ Data already normalized (TPM, FPKM, VST, rlog)
- ✅ Want to assess pre-processed data
- ✅ Data is on log-scale or otherwise ready

**Example**:
```bash
# For TPM data
raptor qc --counts tpm_matrix.csv --normalization none

# For DESeq2 VST data
raptor qc --counts vst_data.csv --normalization none
```

**Supported Pre-Normalized Formats**:
- TPM (Transcripts Per Million)
- FPKM/RPKM (Fragments/Reads Per Kilobase Million)
- VST (DESeq2 variance stabilizing transformation)
- rlog (DESeq2 regularized log)
- Any other log-scale normalized data

### Choosing the Right Normalization

**Decision Tree**:
```
Start
  ↓
Are library sizes similar (CV < 0.5, range < 3x)?
  ├─ Yes → Use log2 (default)
  └─ No → Are library sizes very different (range > 5x)?
          ├─ Yes → Use CPM
          └─ No → Are distributions very different?
                  ├─ Yes → Use quantile (for visualization only!)
                  └─ No → Use CPM to be safe
```

**Quick Reference Table**:

| Scenario | Normalization | Reason |
|----------|---------------|--------|
| Standard RNA-seq, similar depths | `log2` | Simple, effective |
| Library sizes vary (CV > 0.5) | `cpm` | Accounts for depth |
| Failed + good libraries mixed | `cpm` | Handles big differences |
| Strong batch effects (QC only) | `quantile` | Removes technical var |
| TPM/FPKM data | `none` | Already normalized |
| DESeq2 VST/rlog data | `none` | Already transformed |
| Different tissues/conditions | `log2` or `cpm` | Don't use quantile! |

---

## Batch Effect Detection

### IMPROVED in v2.2.0

Previous versions couldn't distinguish batch effects from biological differences. v2.2.0 solves this!

### The Challenge

**Problem**: High PC1 variance could mean:
- ❌ **Batch effect**: Sequencing run, plate, date (BAD - need to correct)
- ✅ **Biological signal**: Treatment, disease, tissue type (GOOD - want to keep!)

**RAPTOR v2.2.0 Solution**: Intelligent column classification + F-statistic comparison

### How It Works

#### Step 1: Column Classification

**Batch-Like Columns** (technical):
```python
batch_keywords = [
    'batch', 'plate', 'run', 'lane', 'date', 
    'flow', 'chip', 'slide', 'array', 
    'seq_run', 'library_prep', 'extraction'
]

# Examples that get classified as "batch":
- sequencing_batch
- plate_number
- library_prep_date
- flow_cell_id
```

**Condition-Like Columns** (biological):
```python
condition_keywords = [
    'condition', 'treatment', 'group', 'status', 
    'disease', 'genotype', 'phenotype', 
    'sample_type', 'cell_type', 'tissue'
]

# Examples that get classified as "condition":
- treatment_group
- disease_status
- tissue_type
- genotype
```

#### Step 2: F-Statistic Calculation

For each column, calculate how much it explains PC1 variance:

```python
# ANOVA F-statistic: How well does this column separate samples in PC1?

F_batch = ANOVA(PC1 ~ sequencing_batch)
F_condition = ANOVA(PC1 ~ treatment)

# Higher F-statistic = stronger effect
```

#### Step 3: Compare Batch vs Condition

```python
if F_batch > F_condition:
    # Batch effect is STRONGER than biological signal
    # This is a problem!
    recommendation = "Strong batch effect detected - correct before analysis"
    
elif F_condition > F_batch:
    # Biological signal is STRONGER than batch
    # This is good!
    recommendation = "Biological signal dominates - batch effect minimal"
    
else:
    # Both are significant
    # Need to include both in model
    recommendation = "Both batch and condition significant - use ~batch + condition"
```

#### Step 4: Check for Confounding

**Confounding** = Batch and condition are perfectly correlated

```python
# Example of confounded design:
Sample1: batch=1, condition=Control
Sample2: batch=1, condition=Control
Sample3: batch=2, condition=Treatment
Sample4: batch=2, condition=Treatment

# Each condition in only one batch!
# CANNOT separate batch effect from treatment effect
```

**Detection**:
```python
confounded = check_confounding(batch, condition)

if confounded:
    recommendation = """
    CONFOUNDED DESIGN: Batch and condition cannot be separated.
    Results may be unreliable.
    Consider repeating experiment with balanced design
    (each condition in multiple batches).
    """
```

### Batch Effect Strength Categories

Based on F-statistic magnitude:

#### Strong Batch Effect (F > 10)

**Recommendation**:
```
Strong batch effect detected. Recommended actions:
  1. Include batch as covariate in DE model: ~batch + condition
  2. For visualization: Use ComBat or limma::removeBatchEffect
  3. Do NOT use batch-corrected data for DE analysis
```

**What To Do**:
```R
# In DESeq2
design(dds) <- ~ batch + condition

# In edgeR
design <- model.matrix(~ batch + condition)

# For PCA visualization only
library(limma)
counts_corrected <- removeBatchEffect(log_counts, batch=metadata$batch)
plotPCA(counts_corrected)  # For visualization ONLY
```

#### Moderate Batch Effect (5 < F ≤ 10)

**Recommendation**:
```
Moderate batch effect detected. Recommended actions:
  1. Include batch as covariate in DE model: ~batch + condition
  2. Monitor batch separation in PCA plots
```

**What To Do**:
- Include batch in model
- Check PCA to verify batch patterns
- May not need visualization correction

#### Weak Batch Effect (F ≤ 5)

**Recommendation**:
```
Weak batch effect. Include batch as covariate if concerned.
```

**What To Do**:
- Optional: Include batch in model for extra confidence
- Probably safe to proceed without correction

#### No Batch Effect (Not Significant)

**Recommendation**:
```
No significant batch effect detected. Proceed with standard analysis.
```

**What To Do**:
- Use simple design: `~condition`
- No batch correction needed

### Example Scenarios

#### Scenario 1: Clear Batch Effect

**Metadata**:
```csv
sample_id,sequencing_batch,treatment
S1,Batch1,Control
S2,Batch1,Control
S3,Batch1,Treatment
S4,Batch2,Control
S5,Batch2,Control
S6,Batch2,Treatment
```

**Results**:
```
Batch variable: sequencing_batch
Batch F-statistic: 15.3 (STRONG)

Condition variable: treatment
Condition F-statistic: 8.2 (moderate)

Confounded: NO (both conditions in both batches)

Recommendation: Strong batch effect - include in model
```

**Action**:
```R
design <- ~ sequencing_batch + treatment
```

#### Scenario 2: Confounded Design

**Metadata**:
```csv
sample_id,sequencing_batch,treatment
S1,Batch1,Control
S2,Batch1,Control
S3,Batch1,Control
S4,Batch2,Treatment
S5,Batch2,Treatment
S6,Batch2,Treatment
```

**Results**:
```
Batch variable: sequencing_batch
Batch F-statistic: 20.5 (VERY STRONG)

Condition variable: treatment
Condition F-statistic: 20.5 (VERY STRONG)

Confounded: YES (each condition in only one batch!)

Recommendation: CONFOUNDED DESIGN - cannot separate batch from treatment
```

**Action**:
```
⚠️  CANNOT perform differential expression!
    Treatment effect is completely confounded with batch.
    
    Solution: Repeat experiment with balanced design.
```

#### Scenario 3: Biological Signal Dominates

**Metadata**:
```csv
sample_id,library_prep_date,tissue_type
S1,2024-01-15,Liver
S2,2024-01-15,Liver
S3,2024-01-20,Liver
S4,2024-01-15,Brain
S5,2024-01-20,Brain
S6,2024-01-20,Brain
```

**Results**:
```
Batch variable: library_prep_date
Batch F-statistic: 3.2 (weak)

Condition variable: tissue_type
Condition F-statistic: 45.8 (VERY STRONG)

Confounded: NO

Recommendation: Biological signal dominates - minimal batch effect
```

**Action**:
```R
# Batch effect is weak, but include for safety
design <- ~ library_prep_date + tissue_type
```

### Critical Don'ts

❌ **DON'T use batch-corrected data for DE analysis**
```R
# WRONG!
counts_corrected <- removeBatchEffect(counts, batch)
run_DE_analysis(counts_corrected)  # ❌ Statistical properties invalid

# CORRECT!
design <- ~ batch + condition
run_DE_analysis(raw_counts, design)  # ✅ Proper statistical model
```

❌ **DON'T remove batch effects before DE if confounded**
```R
# If confounded, you CANNOT separate batch from biology
# No amount of correction will help
# Need to redesign experiment
```

❌ **DON'T ignore batch effects**
```R
# If F-statistic > 5, include batch in model
# Even weak batch effects can create false positives
```

---

## Usage Examples

### Example 1: Basic Quality Check

```bash
# Minimal - just the count matrix
raptor qc --counts results/quick_counts/quick_gene_counts.csv
```

**Output**:
```
═══════════════════════════════════════════════════════════
🦖 RAPTOR DATA QUALITY ASSESSMENT
═══════════════════════════════════════════════════════════
Normalization: log2
Overall Quality Score: 82.3/100
Status: GOOD
Recommendation: Data quality is good. Proceed with analysis.

Component Scores:
  ✓ Library Quality: 88.0/100
  ✓ Gene Detection: 79.0/100
  ✓ Outlier Detection: 90.0/100
  ✓ Variance Structure: 75.0/100
  ✓ Batch Effects: 85.0/100
  ✓ Biological Signal: 77.0/100
═══════════════════════════════════════════════════════════
```

### Example 2: With Metadata (Batch Detection)

```bash
# Include metadata to detect batch effects
raptor qc \
    --counts results/quick_counts/quick_gene_counts.csv \
    --metadata sample_metadata.csv
```

**Metadata Format** (`sample_metadata.csv`):
```csv
sample_id,sequencing_batch,treatment,tissue_type
Sample1,Batch1,Control,Liver
Sample2,Batch1,Control,Liver
Sample3,Batch1,Treatment,Liver
Sample4,Batch2,Control,Liver
Sample5,Batch2,Treatment,Liver
```

### Example 3: CPM Normalization (Varying Library Sizes)

```bash
# When library sizes vary significantly
raptor qc \
    --counts counts.csv \
    --metadata metadata.csv \
    --normalization cpm
```

**Use When**:
```
Library sizes:
  Sample1: 5.2M
  Sample2: 18.7M
  Sample3: 7.3M
  Sample4: 22.1M
  
Range: 4.3x
CV: 0.58

→ CPM normalization recommended
```

### Example 4: With Visualization

```bash
# Generate comprehensive QC plots
raptor qc \
    --counts counts.csv \
    --metadata metadata.csv \
    --plot \
    --plot-output qc_comprehensive.pdf
```

**Plots Generated**:
- Overall quality gauges
- Component score bars
- PCA plot (colored by batch/condition)
- Library size distribution
- Gene detection heatmap
- Correlation heatmap
- Outlier detection plots

### Example 5: Advanced Outlier Detection

```bash
# Conservative outlier detection (4 of 6 methods must agree)
raptor qc \
    --counts counts.csv \
    --consensus-threshold 4 \
    --plot \
    --plot-output outlier_analysis.pdf
```

**Threshold Selection**:
```bash
# Liberal (catch more)
--consensus-threshold 2

# Standard (recommended)
--consensus-threshold 3

# Conservative (high confidence only)
--consensus-threshold 4
```

### Example 6: Python API Usage

```python
from raptor.quality_assessment import DataQualityAssessor
import pandas as pd

# Load data
counts = pd.read_csv('quick_gene_counts.csv', index_col=0)
metadata = pd.read_csv('sample_metadata.csv')

# Initialize assessor
assessor = DataQualityAssessor(counts, metadata, normalization='log2')

# Run comprehensive quality assessment
report = assessor.assess_quality()

# Print results
print(f"Overall Score: {report['overall']['score']:.1f}/100")
print(f"Status: {report['overall']['status']}")

# Check for specific issues
if report['batch_effects']['batch_detected']:
    print("\n⚠️  Batch Effect Detected:")
    print(report['batch_effects']['recommendation'])

# Advanced outlier detection
outlier_result = assessor.detect_outliers_advanced(consensus_threshold=3)
print(f"\n Outliers detected: {outlier_result.outlier_samples}")

# Generate visualizations
assessor.plot_quality_report('qc_report.pdf')
assessor.plot_outliers_advanced(outlier_result, 'outlier_report.pdf')
```

### Example 7: Quick Functions

```python
from raptor.quality_assessment import quick_quality_check, detect_outliers_quick

# Quick quality check with defaults
report = quick_quality_check(counts, metadata)

# Quick outlier detection
outliers = detect_outliers_quick(counts, consensus_threshold=3, plot=True)

print(outliers.summary())
```

### Example 8: Comparing Normalizations

```python
from raptor.quality_assessment import DataQualityAssessor

# Test different normalizations
for norm in ['log2', 'cpm', 'quantile']:
    assessor = DataQualityAssessor(counts, metadata, normalization=norm)
    report = assessor.assess_quality()
    print(f"{norm}: Score = {report['overall']['score']:.1f}")

# Output:
# log2: Score = 75.3
# cpm: Score = 82.1     ← Best for this data
# quantile: Score = 88.4 (but removes biology!)
```

---

## Interpreting Results

### Understanding Quality Scores

#### Excellent (90-100)

**Interpretation**:
- High-quality data ready for analysis
- Minimal technical issues
- Strong biological signal

**Next Steps**:
1. Proceed to Module 3 (Profiling)
2. No special preprocessing needed
3. Standard analysis pipeline

#### Good (80-89)

**Interpretation**:
- Good quality data
- Minor issues present
- Proceed with standard analysis

**Next Steps**:
1. Note any flagged issues
2. Include batch as covariate if detected
3. Proceed to Module 3

#### Acceptable (60-79)

**Interpretation**:
- Usable data with notable issues
- May need corrective actions
- Results should be interpreted carefully

**Common Issues**:
- Moderate batch effects → Include in model
- 1-2 outliers detected → Investigate, may remove
- High library size variation → Use CPM normalization

**Next Steps**:
1. Address specific flagged issues
2. Consider removing clear outliers
3. Use appropriate normalization
4. Include batch in DE model
5. Proceed with caution

#### Poor (< 60)

**Interpretation**:
- Significant quality issues
- High risk of false results
- Investigation required before analysis

**Common Causes**:
- Multiple outliers (>10% samples)
- Strong batch effects confounded with biology
- Very low library sizes
- Extremely high zero inflation
- Minimal biological signal

**Next Steps**:
1. **Investigate**:
   - Check sample prep notes
   - Review sequencing quality metrics
   - Examine PCA plots
   
2. **Possible Actions**:
   - Remove failed samples
   - Re-sequence low-depth samples
   - Redesign experiment if confounded
   
3. **Do NOT proceed with DE until resolved**

### Reading Component Reports

#### Library Quality Example

```json
{
  "mean_size": 8234567,
  "cv": 0.68,
  "min_size": 3123456,
  "max_size": 18234567,
  "size_range": 5.8,
  "score": 65,
  "status": "warning",
  "flags": [
    "High library size variation",
    "Minimum library size < 5M reads"
  ]
}
```

**Interpretation**:
- **CV = 0.68**: High variation (>0.5 threshold)
- **Min = 3.1M**: Below ideal (5M+ recommended)
- **Range = 5.8x**: Large difference between samples
- **Score = 65**: Acceptable but not ideal

**Action**:
```bash
# Use CPM normalization to account for depth differences
raptor qc --counts counts.csv --normalization cpm

# Consider removing samples < 1M reads
# Or re-sequence very low samples
```

#### Batch Effect Example

```json
{
  "batch_detected": true,
  "batch_variable": "sequencing_run",
  "batch_strength": 12.4,
  "condition_variable": "treatment",
  "condition_strength": 8.7,
  "confounded": false,
  "score": 70,
  "recommendation": "Strong batch effect - include in model"
}
```

**Interpretation**:
- **Batch F-stat = 12.4**: Strong batch effect
- **Condition F-stat = 8.7**: Moderate biological signal
- **Batch > Condition**: Technical effect stronger than biology
- **Not confounded**: Can separate batch from treatment

**Action**:
```R
# Include batch in DE model
design <- ~ sequencing_run + treatment

# For visualization only
corrected <- removeBatchEffect(counts, batch)
plotPCA(corrected)
```

### Decision Trees

#### Should I Remove Outliers?

```
Are outliers detected?
  ├─ NO → Proceed to Module 3
  └─ YES → How many methods flagged each outlier?
           ├─ 5-6 methods → Very likely technical outlier
           │   └─ Action: Remove from analysis
           │
           ├─ 3-4 methods → Probable outlier
           │   └─ Investigate:
           │       • Check metadata for known issues
           │       • Examine in PCA plot
           │       • If technical issue: Remove
           │       • If biological variation: Keep
           │
           └─ 1-2 methods → Borderline
               └─ Action: Keep unless clear technical issue
```

#### Which Normalization Should I Use?

```
What type of data do you have?
  ├─ Pre-normalized (TPM, FPKM, VST) → Use 'none'
  │
  ├─ Raw counts with similar library sizes (CV < 0.5)
  │   └─ Use 'log2' (default)
  │
  ├─ Raw counts with varying sizes (CV > 0.5)
  │   └─ Use 'cpm'
  │
  └─ Very different distributions + batch effects
      └─ Use 'quantile' (VISUALIZATION ONLY!)
```

#### Should I Include Batch in Model?

```
Is batch effect detected?
  ├─ NO → Use simple model: ~condition
  │
  └─ YES → Is design confounded?
           ├─ YES (confounded) → STOP
           │   └─ Cannot separate batch from biology
           │       Redesign experiment needed
           │
           └─ NO (not confounded) → What's batch strength?
                ├─ F > 10 (strong) → MUST include: ~batch + condition
                ├─ F > 5 (moderate) → SHOULD include: ~batch + condition
                └─ F ≤ 5 (weak) → OPTIONAL: ~batch + condition
```

---

## Integration with RAPTOR Workflow

### The Complete Pipeline

```
Module 1: Quick Quantification
    ↓
    quick_gene_counts.csv
    ↓
Module 2: Quality Assessment (THIS MODULE)
    ↓
    Quality report + Outlier flags
    ↓
    Decision Point:
    • Good quality (>80) → Continue
    • Issues detected → Fix issues
    • Poor quality (<60) → Investigate
    ↓
Module 3: Data Profiling
    ↓
    32-feature characterization
    ↓
Module 4: Pipeline Recommendation
    ↓
    Optimal DE tools and parameters
    ↓
Module 5: Production Quantification (if needed)
    ↓
    High-quality counts
    ↓
Module 6-7: Differential Expression
```

### Handoff to Module 3

```bash
# After Module 2 QC passes
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# What Module 3 does with quality information:
# 1. Uses quality assessment to guide profiling
# 2. Considers outlier flags in feature calculation
# 3. Incorporates batch information in recommendations
```

### Typical Workflow Scenarios

#### Scenario 1: Clean Data

```bash
# Module 1
raptor quick-count -m salmon -s samples.csv -i index/

# Module 2
raptor qc --counts results/quick_counts/quick_gene_counts.csv
# Result: Score = 87 (Good) ✓

# Module 3
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# Module 4
raptor recommend
# → DESeq2 recommended

# Continue to DE analysis
```

#### Scenario 2: Outliers Detected

```bash
# Module 2
raptor qc --counts counts.csv --metadata metadata.csv
# Result: 2 outliers detected (Sample_12, Sample_45)

# Investigate outliers
# → Sample_12: Library size 450k (failed)
# → Sample_45: Wrong tissue type (metadata error)

# Remove outliers
# Edit sample sheet, remove bad samples

# Re-run Module 1 without outliers
raptor quick-count -m salmon -s samples_filtered.csv -i index/

# Re-run Module 2
raptor qc --counts results/quick_counts/quick_gene_counts.csv
# Result: Score = 91 (Excellent) ✓

# Continue to Module 3
```

#### Scenario 3: Batch Effects

```bash
# Module 2
raptor qc --counts counts.csv --metadata metadata.csv
# Result: Strong batch effect detected
#         Batch F-stat = 15.3, Condition F-stat = 8.2
#         Recommendation: Include batch in model

# Note the batch variable
# batch_column = "sequencing_run"

# Continue to Module 3
raptor profile --counts counts.csv --metadata metadata.csv

# Module 4 will recommend:
# DESeq2 with design ~ sequencing_run + treatment

# When running DE:
raptor de --counts counts.csv \
          --metadata metadata.csv \
          --design "~ sequencing_run + treatment"
```

#### Scenario 4: Confounded Design

```bash
# Module 2
raptor qc --counts counts.csv --metadata metadata.csv
# Result: CONFOUNDED DESIGN
#         Each treatment in only one batch
#         Cannot separate batch from treatment

# Action: STOP
# Cannot proceed with differential expression
# Results would be unreliable

# Options:
# 1. Repeat experiment with balanced design
# 2. Use as exploratory analysis only
# 3. Report limitations in publication
```

---

## Troubleshooting

### Common Issues

#### Issue 1: "Count matrix validation failed"

**Error**:
```
ValidationError: Count matrix contains negative values
```

**Cause**: Negative counts in matrix

**Solution**:
```bash
# Check for negative values
python3 -c "
import pandas as pd
counts = pd.read_csv('counts.csv', index_col=0)
print('Negative values:', (counts < 0).sum().sum())
print('Min value:', counts.min().min())
"

# If accidentally using log-transformed data:
# Use --normalization none
raptor qc --counts counts.csv --normalization none
```

#### Issue 2: Very Low Quality Score

**Symptom**: Overall score < 40

**Common Causes**:
1. **Multiple outliers**: >10% samples flagged
2. **Very low library sizes**: Mean < 1M
3. **High zero inflation**: >80% zeros
4. **Confounded batch effects**

**Diagnosis**:
```bash
# Run with plots to visualize issues
raptor qc --counts counts.csv --plot --plot-output diagnosis.pdf

# Check component scores
# Which component is lowest?
```

**Solutions**:
```bash
# If library size issues:
# 1. Filter low-quality samples
# 2. Re-sequence if possible
# 3. Use CPM normalization

# If zero inflation:
# 1. May be normal for single-cell or sparse data
# 2. Consider specialized tools (scRNA-seq pipelines)
# 3. Filter very low-expressed genes

# If batch confounding:
# 1. Cannot proceed with standard analysis
# 2. Redesign experiment
```

#### Issue 3: No Batch Effects Detected (but you know there are batches)

**Symptom**: Batch effect score = 100, no batch detected

**Cause**: Metadata column not recognized as batch

**Solution**:
```bash
# Check metadata column names
cat sample_metadata.csv | head -1

# Ensure column name includes batch keyword:
# ✓ Good: "sequencing_batch", "batch", "plate"
# ✗ Bad: "seq", "group1", "run_id"

# Rename column:
# Before: seq_run
# After: sequencing_batch

# Or manually in metadata file:
seq_run -> sequencing_batch
```

**Batch Keywords** (case-insensitive):
```python
['batch', 'plate', 'run', 'lane', 'date', 
 'flow', 'chip', 'slide', 'array', 
 'seq_run', 'library_prep', 'extraction']
```

#### Issue 4: Too Many False Positive Outliers

**Symptom**: 30% of samples flagged as outliers

**Cause**: Consensus threshold too low

**Solution**:
```bash
# Increase consensus threshold
raptor qc --counts counts.csv --consensus-threshold 4

# Or if samples are actually very different:
# May not be outliers but biological diversity
# Check PCA plot colored by biological groups
```

#### Issue 5: All Samples Pass (but data looks bad)

**Symptom**: High scores despite obvious issues

**Cause**: Wrong normalization for data type

**Solution**:
```bash
# If data is already normalized (TPM, FPKM):
raptor qc --counts tpm_data.csv --normalization none

# If library sizes vary greatly:
raptor qc --counts counts.csv --normalization cpm

# Check library sizes first:
python3 -c "
import pandas as pd
counts = pd.read_csv('counts.csv', index_col=0)
lib_sizes = counts.sum(axis=0)
print('Library size range:', lib_sizes.min(), '-', lib_sizes.max())
print('CV:', lib_sizes.std() / lib_sizes.mean())
"
```

#### Issue 6: Memory Error

**Symptom**:
```
MemoryError: Unable to allocate array
```

**Cause**: Very large count matrix (many genes/samples)

**Solutions**:
```bash
# 1. Filter lowly expressed genes first
python3 << EOF
import pandas as pd
counts = pd.read_csv('counts.csv', index_col=0)
# Keep genes with mean count > 10
keep = counts.mean(axis=1) > 10
counts_filtered = counts[keep]
counts_filtered.to_csv('counts_filtered.csv')
EOF

raptor qc --counts counts_filtered.csv

# 2. Use server with more RAM

# 3. Process in batches (advanced)
```

### Performance Optimization

#### Large Datasets

**For 100+ samples or 100k+ genes**:

```python
from raptor.quality_assessment import DataQualityAssessor

# 1. Filter genes first
counts = pd.read_csv('counts.csv', index_col=0)
expressed = (counts > 0).sum(axis=1) > counts.shape[1] * 0.1
counts_filtered = counts[expressed]

# 2. Use faster normalization
assessor = DataQualityAssessor(
    counts_filtered, 
    metadata,
    normalization='log2'  # Faster than quantile
)

# 3. Skip plots if not needed
report = assessor.assess_quality()
# Don't call plot methods
```

#### Quick Check

**For very quick assessment (<10 seconds)**:

```python
from raptor.quality_assessment import quick_quality_check

# Minimal quality check
report = quick_quality_check(
    counts, 
    metadata, 
    plot=False  # Skip plotting
)

print(f"Score: {report['overall']['score']}")
```

---

## Best Practices

### 1. Always Use Metadata

❌ **Don't**:
```bash
raptor qc --counts counts.csv
# Missing batch effect detection!
```

✅ **Do**:
```bash
raptor qc --counts counts.csv --metadata metadata.csv
# Enables batch effect detection
```

**Why**: Cannot detect batch effects without metadata

### 2. Choose Appropriate Normalization

❌ **Don't**:
```bash
# Using log2 when library sizes vary 10x
raptor qc --counts counts.csv --normalization log2
```

✅ **Do**:
```bash
# Check library size variation first
# If CV > 0.5 or range > 5x, use CPM
raptor qc --counts counts.csv --normalization cpm
```

### 3. Investigate Outliers Before Removing

❌ **Don't**:
```python
# Blindly removing outliers
outliers = result.outlier_samples
counts_clean = counts.drop(outliers, axis=1)
```

✅ **Do**:
```python
# Investigate first
for sample in result.outlier_samples:
    votes = result.outlier_scores[sample]
    print(f"{sample}: flagged by {votes}/6 methods")
    
    # Check metadata
    print(metadata[metadata['sample_id'] == sample])
    
    # Check library size
    lib_size = counts[sample].sum()
    print(f"Library size: {lib_size:,}")
    
# Only remove if clear technical issue
```

### 4. Don't Use Quantile Normalization for DE

❌ **Don't**:
```bash
# Using quantile-normalized data for DE
raptor qc --counts counts.csv --normalization quantile
# Then use this for differential expression ❌
```

✅ **Do**:
```bash
# Use quantile only for QC visualization
raptor qc --counts counts.csv --normalization quantile --plot

# Use raw counts for DE analysis
raptor de --counts counts.csv  # Uses raw counts
```

### 5. Save QC Reports

✅ **Do**:
```bash
# Save reports with dates
raptor qc --counts counts.csv \
    --output qc_report_$(date +%Y%m%d).json \
    --plot \
    --plot-output qc_plots_$(date +%Y%m%d).pdf

# Version control
git add qc_report_20240309.json
git commit -m "QC report: 87/100, 2 outliers detected"
```

### 6. Check Confounding Early

✅ **Do**:
```bash
# Always run QC before starting DE analysis
raptor qc --counts counts.csv --metadata metadata.csv

# If confounded design detected:
# STOP and redesign experiment
# Don't waste time on unreliable analysis
```

### 7. Use Appropriate Consensus Threshold

**Guidelines**:
```bash
# Exploratory analysis (find all potential issues)
--consensus-threshold 2

# Standard analysis (balanced)
--consensus-threshold 3

# Publication (high confidence only)
--consensus-threshold 4
```

### 8. Document Your Decisions

✅ **Do**:
```markdown
# QC_NOTES.md

## Quality Assessment - 2024-03-09

### Run
- Count matrix: quick_gene_counts.csv (60 samples, 20,145 genes)
- Normalization: CPM (library sizes vary 4.5x)
- Consensus threshold: 3

### Results
- Overall score: 82/100 (Good)
- Outliers detected: 2 (Sample_42, Sample_89)
  - Sample_42: Library size 450k (failed sequencing)
  - Sample_89: Wrong tissue type (metadata error)
- Batch effect: Moderate (F=7.3)
  - Include sequencing_run in model

### Actions Taken
1. Removed Sample_42 (technical failure)
2. Removed Sample_89 (sample mix-up)
3. Re-ran Module 1 with 58 samples
4. New QC score: 91/100 (Excellent)
5. Will include batch in DE model: ~sequencing_run + treatment

### Files
- Original: qc_report_20240309_original.json
- Filtered: qc_report_20240309_filtered.json
- Plots: qc_plots_20240309.pdf
```

---

## Technical Details

### Quality Score Calculation

**Component Weights**:
```python
weights = {
    'library_quality': 0.15,      # 15%
    'gene_detection': 0.20,       # 20%
    'outlier_detection': 0.15,    # 15%
    'variance_structure': 0.15,   # 15%
    'batch_effects': 0.20,        # 20%
    'biological_signal': 0.15     # 15%
}

overall_score = sum(component_score × weight for all components)
```

**Why These Weights?**:
- **Gene Detection (20%)**: Most important - poor detection = poor data
- **Batch Effects (20%)**: Critical - can create false results
- **Library Quality (15%)**: Important but less critical
- **Outlier Detection (15%)**: Affects analysis but can be handled
- **Variance Structure (15%)**: Diagnostic, not necessarily bad
- **Biological Signal (15%)**: Nice to have, but low signal is analyzable

### Outlier Detection Algorithms

#### PCA + Mahalanobis

**Mathematics**:
```python
# 1. PCA transformation
X_pca = PCA(n_components=5).fit_transform(X_scaled)

# 2. Mahalanobis distance
mean = X_pca.mean(axis=0)
cov = np.cov(X_pca.T)
inv_cov = np.linalg.pinv(cov)

for sample in samples:
    diff = sample - mean
    D² = diff @ inv_cov @ diff.T
    D = sqrt(D²)
    
    if D > mean(D) + 3×std(D):
        flagged_as_outlier

```

**Assumptions**:
- Multivariate normal distribution in PC space
- Linear relationships

#### Isolation Forest

**Algorithm**:
```python
# 1. Build random trees
for tree in forest:
    randomly select feature and split point
    isolate samples recursively
    
# 2. Anomaly score
score = average path length to isolate sample

# 3. Shorter paths = outliers
# (easier to isolate = different from others)
```

**Advantages**:
- Non-parametric
- Handles non-linear relationships
- Fast

#### Local Outlier Factor

**Algorithm**:
```python
# 1. Find k nearest neighbors for each sample
neighbors = find_k_nearest(sample, k=10)

# 2. Calculate local density
local_density = 1 / mean_distance_to_neighbors

# 3. Compare to neighbor densities
LOF = neighbor_density / local_density

# 4. LOF > 1 = outlier (lower local density)
```

**Advantages**:
- Handles varying densities
- Detects local outliers

### Batch Effect F-Statistic

**ANOVA F-statistic**:
```python
# For categorical variable (batch) explaining continuous (PC1):

# Between-group variance
SS_between = sum(n_i × (mean_i - grand_mean)²)
df_between = k - 1  # k = number of groups

# Within-group variance
SS_within = sum(sum((x_ij - mean_i)²))
df_within = N - k  # N = total samples

# F-statistic
F = (SS_between / df_between) / (SS_within / df_within)

# Higher F = variable explains more variance
```

**Interpretation**:
- F > 10: Very strong effect
- F > 5: Strong effect
- F > 2: Moderate effect
- F ≤ 2: Weak effect

### Normalization Methods Details

#### CPM Calculation

```python
# For each sample:
library_size = sum(counts_for_sample)
cpm = (counts / library_size) × 1,000,000

# Then log transform:
log_cpm = log2(cpm + 1)
```

**Why +1?**: Handles zero counts (log(0) undefined)

#### Quantile Normalization

```python
# Algorithm:
1. Sort each sample's values
2. Compute mean at each rank across samples
3. Replace each sample's values with rank means
4. Restore original order within each sample

# Effect: All samples have identical distributions
```

**Use Case**: Remove technical variation for visualization

---

## Summary

### Key Takeaways

1. **Module 2 Purpose**: Identify quality issues before analysis

2. **6 Assessment Components**:
   - Library quality, gene detection, outliers
   - Variance structure, batch effects, biological signal

3. **v2.2.0 Enhancements**:
   - 6-method consensus outlier detection
   - Batch vs condition separation

4. **Critical**: Check for confounded designs early

5. **Normalization**: Choose based on data characteristics

6. **Integration**: Informs Module 3, guides DE setup

### Workflow Checklist

- [ ] Run Module 2 on quick counts from Module 1
- [ ] Include metadata for batch detection
- [ ] Review overall quality score
- [ ] Investigate flagged outliers
- [ ] Check for batch effects and confounding
- [ ] Choose appropriate normalization
- [ ] Generate and review QC plots
- [ ] Document decisions and actions
- [ ] Remove technical outliers if identified
- [ ] Note batch variables for DE model
- [ ] Proceed to Module 3 if quality acceptable

### Quality Thresholds

**Proceed to Module 3**:
- Score ≥ 60
- No confounded design
- < 10% outliers (or investigated)

**Fix Issues Before Proceeding**:
- Score < 60
- Confounded design detected
- > 20% outliers

**May Need Re-Sequencing**:
- Score < 40
- Many failed libraries (< 1M reads)
- > 80% zero inflation

---

**Document Version**: 2.2.0  
**Last Updated**: March 2026  
**Author**: Ayeh Bolouki  
**License**: MIT
