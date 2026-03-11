# RAPTOR Module 9: Ensemble Analysis - Comprehensive User Guide

**Version:** 2.2.0  
**Author:** Ayeh Bolouki  
**Last Updated:** January 2026  

---

## Table of Contents

1. [Overview](#overview)
2. [Scientific Background](#scientific-background)
3. [Available Ensemble Methods](#available-ensemble-methods)
4. [Installation & Setup](#installation--setup)
5. [CLI Usage](#cli-usage)
6. [Python API Usage](#python-api-usage)
7. [Input Requirements](#input-requirements)
8. [Output Files](#output-files)
9. [Method Selection Guide](#method-selection-guide)
10. [Complete Examples](#complete-examples)
11. [Parameter Tuning](#parameter-tuning)
12. [Troubleshooting](#troubleshooting)
13. [Integration with Other Modules](#integration-with-other-modules)
14. [Scientific References](#scientific-references)

---

## Overview

Module 9 provides comprehensive ensemble analysis capabilities for combining differential expression (DE) results from multiple methods (DESeq2, edgeR, limma, Wilcoxon) to create high-confidence consensus gene lists.

### Why Use Ensemble Analysis?

**Problem:** Different DE methods often produce different gene lists due to:
- Different statistical models and assumptions
- Different normalization approaches
- Different handling of outliers and low counts
- Different sensitivity to biological variation

**Solution:** Ensemble methods combine results to:
- ✅ Identify genes with consistent evidence across methods
- ✅ Reduce method-specific false positives
- ✅ Create prioritized gene lists for validation
- ✅ Provide statistical confidence measures
- ✅ Enable publication-quality consensus results

### Key Features

- 🔬 **5 Ensemble Methods**: Fisher's, Brown's, RRA, Voting, Weighted
- 📊 **Direction Consistency Checking**: Ensures genes have consistent regulation direction
- 🎯 **Tiered Confidence Levels**: Very High, High, Moderate confidence genes
- 📈 **Statistical Rigor**: P-value combination, rank aggregation, meta-analysis
- 🔄 **Flexible Integration**: Works with Module 7 DE results
- 💾 **Comprehensive Outputs**: CSV files, statistics, visualizations

---

## Scientific Background

### The Ensemble Principle

Ensemble analysis in bioinformatics combines evidence from multiple statistical tests to achieve:

1. **Higher Specificity**: Genes detected by multiple methods are less likely to be false positives
2. **Better Reproducibility**: Consensus genes show more consistent results across studies
3. **Robust Rankings**: Combined rankings are less sensitive to outliers
4. **Statistical Power**: Meta-analytical approaches increase detection power

### When to Use Each Method

| Method | Best For | Sensitivity | Specificity | Computational Cost |
|--------|----------|-------------|-------------|-------------------|
| **Fisher's** | Exploratory analysis, maximum gene discovery | Very High | Moderate | Low |
| **Brown's** | Correlated methods (same normalization) | Very High | Moderate | Low |
| **RRA** | Robust ranking, balanced results | High | High | Medium |
| **Voting** | High-confidence validation targets | Moderate | Very High | Low |
| **Weighted** | Quality-based, leverages Module 8 weights | High | High | Low |

---

## Available Ensemble Methods

### 1. Fisher's Method (P-value Combination)

**Principle:** Combines p-values using Fisher's combined probability test

**Formula:**
```
χ² = -2 Σ ln(p_i)  where i = 1 to k methods
Combined p-value from chi-squared distribution with 2k degrees of freedom
```

**Best For:**
- Maximum sensitivity
- Exploratory analysis
- Large-scale gene discovery

**Advantages:**
- Simple, well-established method
- Maximum gene detection
- Fast computation

**Limitations:**
- Assumes independence between methods
- May have more false positives

**Reference:** Fisher, R.A. (1925). Statistical Methods for Research Workers.

---

### 2. Brown's Method (Correlation-Aware)

**Principle:** Extension of Fisher's method that accounts for correlation between tests

**Formula:**
```
Scale factor κ accounts for correlation between p-values
χ² = κ × (-2 Σ ln(p_i))
```

**Best For:**
- Methods using same data/normalization
- When p-values are correlated
- More accurate p-value combination

**Advantages:**
- Accounts for method correlation
- More accurate than Fisher's when methods are related
- Still computationally efficient

**Limitations:**
- Requires correlation estimation
- Slightly more complex than Fisher's

**Reference:** Brown, M.B. (1975). A method for combining non-independent tests. Biometrics.

---

### 3. Robust Rank Aggregation (RRA)

**Principle:** Combines ranked gene lists using order statistics

**Algorithm:**
1. Rank genes by p-value in each method
2. For each gene, compute probability of observing ranks this high
3. Correct for multiple testing
4. Select genes below significance threshold

**Best For:**
- Balanced sensitivity/specificity
- Robust to outliers
- Publication-quality results

**Advantages:**
- Doesn't require p-value assumptions
- Robust to extreme values
- Well-suited for meta-analysis

**Limitations:**
- Slightly more computational cost
- May miss some true positives

**Reference:** Kolde, R., et al. (2012). Robust rank aggregation. Bioinformatics.

---

### 4. Simple Voting (Count-Based)

**Principle:** Genes detected by ≥ min_methods are included

**Algorithm:**
1. Filter each method by significance threshold
2. Count how many methods detect each gene
3. Include genes detected by ≥ min_methods

**Best For:**
- High-confidence validation targets
- Strict quality control
- Publication main results

**Advantages:**
- Extremely simple to interpret
- Very high confidence
- Flexible threshold (≥2, ≥3, etc.)

**Limitations:**
- May be too conservative (fewer genes)
- Binary decision (no ranking within consensus)

---

### 5. Weighted Ensemble (Performance-Based)

**Principle:** Weight each method by quality metrics from Module 8

**Algorithm:**
1. Load method weights from optimization results
2. Compute weighted score for each gene
3. Include genes above min_score threshold

**Best For:**
- Leveraging Module 8 optimization
- Quality-based prioritization
- Advanced users

**Advantages:**
- Incorporates method performance
- Adaptive to dataset characteristics
- Can use custom weights

**Limitations:**
- Requires Module 8 results or manual weights
- More complex to interpret

---

## Installation & Setup

### Prerequisites

```bash
# Python 3.8+
python --version

# Required packages
pip install numpy pandas scipy --break-system-packages
```

### Verify Installation

```bash
# Check if module is available
python -c "from raptor.ensemble import ensemble_fisher; print('✓ Module 9 available')"

# Check CLI commands
raptor ensemble --help
raptor ensemble-compare --help
```

### Import Structure

```python
# Import specific functions
from raptor.ensemble import (
    ensemble_fisher,      # Fisher's method
    ensemble_brown,       # Brown's method  
    ensemble_rra,         # Robust Rank Aggregation
    ensemble_voting,      # Simple voting
    ensemble_weighted     # Weighted ensemble
)

# Import result class
from raptor.ensemble import EnsembleResult

# Import lower-level functions
from raptor.ensemble import (
    fishers_method,
    browns_method,
    robust_rank_aggregation,
    check_direction_consistency
)
```

---

## CLI Usage

### Command 1: `raptor ensemble`

Run one or more ensemble methods on DE results.

#### Basic Syntax

```bash
raptor ensemble \
    --methods METHOD [METHOD ...] \
    --deseq2 PATH \
    --edger PATH \
    --limma PATH \
    [--wilcoxon PATH] \
    --output DIR \
    [OPTIONS]
```

#### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--methods, -m` | Ensemble method(s) to use | `fisher`, `brown`, `rra` |
| `--output, -o` | Output directory | `results/ensemble/` |
| **At least 2 of:** | | |
| `--deseq2` | DESeq2 result file (.pkl) | `deseq2_result.pkl` |
| `--edger` | edgeR result file (.pkl) | `edger_result.pkl` |
| `--limma` | limma result file (.pkl) | `limma_result.pkl` |
| `--wilcoxon` | Wilcoxon result file (.pkl) | `wilcoxon_result.pkl` |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--threshold` | 0.05 | Significance threshold |
| `--min-methods` | 2 | Minimum methods for voting/weighted |

#### Examples

**Example 1: Fisher's method only**
```bash
raptor ensemble \
    -m fisher \
    --deseq2 output/de_results/deseq2_result.pkl \
    --edger output/de_results/edger_result.pkl \
    --limma output/de_results/limma_result.pkl \
    -o results/ensemble_fisher/
```

**Example 2: Multiple methods**
```bash
raptor ensemble \
    -m fisher -m brown -m rra \
    --deseq2 deseq2_result.pkl \
    --edger edger_result.pkl \
    --limma limma_result.pkl \
    -o results/ensemble_all/
```

**Example 3: Custom threshold**
```bash
raptor ensemble \
    -m rra \
    --deseq2 deseq2_result.pkl \
    --edger edger_result.pkl \
    --threshold 0.01 \
    -o results/ensemble_strict/
```

**Example 4: Four methods**
```bash
raptor ensemble \
    -m fisher -m rra \
    --deseq2 deseq2_result.pkl \
    --edger edger_result.pkl \
    --limma limma_result.pkl \
    --wilcoxon wilcoxon_result.pkl \
    -o results/ensemble_4methods/
```

---

### Command 2: `raptor ensemble-compare`

Compare all 5 ensemble methods side-by-side.

#### Basic Syntax

```bash
raptor ensemble-compare \
    --deseq2 PATH \
    --edger PATH \
    --limma PATH \
    --output DIR \
    [OPTIONS]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--deseq2` | DESeq2 result file (.pkl) |
| `--edger` | edgeR result file (.pkl) |
| `--limma` | limma result file (.pkl) |
| `--output, -o` | Output directory |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--threshold` | 0.05 | Significance threshold |

#### Example

```bash
raptor ensemble-compare \
    --deseq2 output/de_results/deseq2_result.pkl \
    --edger output/de_results/edger_result.pkl \
    --limma output/de_results/limma_result.pkl \
    -o results/ensemble_comparison/ \
    --threshold 0.05
```

#### Output

Creates subdirectories for each method plus a comparison summary:
```
results/ensemble_comparison/
├── fisher/
│   ├── consensus_genes.csv
│   ├── combined_pvalues.csv
│   └── statistics.json
├── brown/
├── rra/
├── voting/
├── weighted/
└── comparison_summary.csv
```

---

## Python API Usage

### Basic Usage Pattern

```python
import pickle
from raptor.ensemble import ensemble_fisher, ensemble_brown, ensemble_rra

# 1. Load DE results from Module 7
with open('deseq2_result.pkl', 'rb') as f:
    deseq2_result = pickle.load(f)
with open('edger_result.pkl', 'rb') as f:
    edger_result = pickle.load(f)
with open('limma_result.pkl', 'rb') as f:
    limma_result = pickle.load(f)

# 2. Organize into dictionary
de_results = {
    'DESeq2': deseq2_result,
    'edgeR': edger_result,
    'limma': limma_result
}

# 3. Run ensemble analysis
result = ensemble_fisher(
    de_results=de_results,
    significance_threshold=0.05,
    check_direction=True,
    output_dir='output/fisher/'
)

# 4. Access results
print(f"Found {result.n_consensus_genes} consensus genes")
print(result.summary())

# 5. Get consensus gene list
consensus_genes = result.consensus_genes
print(consensus_genes.head())

# 6. Save results
result.save('output/fisher/')
```

### Function Signatures

#### 1. `ensemble_fisher()` / `ensemble_brown()`

```python
def ensemble_pvalue_combination(
    de_results: Dict[str, DEResult],
    method: str = 'fisher',           # 'fisher' or 'brown'
    use_padj: bool = False,           # Use adjusted p-values?
    significance_threshold: float = 0.05,
    check_direction: bool = True,
    direction_threshold: float = 1.0,  # 1.0 = all must agree
    output_dir: Optional[Path] = None
) -> EnsembleResult
```

**Parameters:**
- `de_results`: Dictionary mapping method names to DEResult objects
- `method`: 'fisher' or 'brown'
- `use_padj`: If True, use adjusted p-values; if False, use raw p-values (RECOMMENDED: False)
- `significance_threshold`: Combined p-value threshold (default: 0.05)
- `check_direction`: Filter for direction consistency (default: True)
- `direction_threshold`: Fraction of methods that must agree on direction (default: 1.0)
- `output_dir`: Where to save results (optional)

**Returns:**
- `EnsembleResult` object with consensus genes, statistics, and metadata

---

#### 2. `ensemble_rra()`

```python
def ensemble_rra(
    de_results: Dict[str, DEResult],
    rank_by: str = 'pvalue',          # 'pvalue' or 'padj'
    use_padj: bool = False,
    significance_threshold: float = 0.05,
    check_direction: bool = True,
    direction_threshold: float = 1.0,
    output_dir: Optional[Path] = None
) -> EnsembleResult
```

**Parameters:**
- `de_results`: Dictionary mapping method names to DEResult objects
- `rank_by`: What to rank by ('pvalue' recommended)
- `use_padj`: If True, rank by padj; if False, by pvalue
- `significance_threshold`: RRA p-value threshold
- `check_direction`: Filter for direction consistency
- `direction_threshold`: Fraction agreement required
- `output_dir`: Where to save results

**Returns:**
- `EnsembleResult` object with ranked genes and consensus list

---

#### 3. `ensemble_voting()`

```python
def ensemble_voting(
    de_results: Dict[str, DEResult],
    min_methods: int = 2,              # Minimum methods required
    filters: Dict[str, float] = None,  # e.g., {'padj': 0.05, 'lfc': 0.0}
    check_direction: bool = True,
    direction_threshold: float = 1.0,
    output_dir: Optional[Path] = None
) -> EnsembleResult
```

**Parameters:**
- `de_results`: Dictionary mapping method names to DEResult objects
- `min_methods`: Minimum number of methods that must detect gene (default: 2)
- `filters`: Significance filters to apply per method
- `check_direction`: Filter for direction consistency
- `direction_threshold`: Fraction agreement required
- `output_dir`: Where to save results

**Returns:**
- `EnsembleResult` object with voted consensus genes

---

#### 4. `ensemble_weighted()`

```python
def ensemble_weighted(
    de_results: Dict[str, DEResult],
    weights: Optional[Dict[str, float]] = None,  # Method weights
    min_score: float = 1.5,           # Minimum weighted score
    filters: Dict[str, float] = None,
    check_direction: bool = True,
    output_dir: Optional[Path] = None
) -> EnsembleResult
```

**Parameters:**
- `de_results`: Dictionary mapping method names to DEResult objects
- `weights`: Dictionary of method weights (default: equal weights)
- `min_score`: Minimum weighted score threshold
- `filters`: Significance filters per method
- `check_direction`: Filter for direction consistency
- `output_dir`: Where to save results

**Returns:**
- `EnsembleResult` object with weighted consensus genes

---

### EnsembleResult Class

The result object contains all outputs and methods:

```python
class EnsembleResult:
    # Main outputs
    consensus_genes: pd.DataFrame        # Final consensus gene list
    n_consensus_genes: int               # Number of consensus genes
    
    # Method details
    ensemble_method: str                 # Which method was used
    n_methods: int                       # Number of input methods
    method_names: List[str]              # Names of methods
    
    # Direction analysis
    direction_consistency: pd.DataFrame  # Per-gene direction agreement
    n_direction_inconsistent: int        # Genes with inconsistent direction
    
    # Statistics
    method_statistics: Dict              # Per-method and overall stats
    
    # Method-specific outputs (optional)
    combined_pvalues: pd.DataFrame       # For Fisher/Brown
    ranked_genes: pd.DataFrame           # For RRA
    
    # Metadata
    parameters: Dict                     # Parameters used
    timestamp: str                       # When analysis was run
    
    # Methods
    def summary() -> str                 # Print summary
    def save(output_dir: Path)           # Save all results
    @classmethod
    def load(output_dir: Path) -> EnsembleResult  # Load saved results
```

#### Accessing Results

```python
# Get consensus genes
genes_df = result.consensus_genes
print(f"Columns: {genes_df.columns.tolist()}")

# Check statistics
stats = result.method_statistics
print(f"Total genes: {stats['overall']['n_total_genes']}")
print(f"Consensus genes: {stats['overall']['n_consensus_genes']}")

# Direction consistency
direction_df = result.direction_consistency
inconsistent = direction_df[~direction_df['is_consistent']]
print(f"Direction inconsistent: {len(inconsistent)}")

# Method-specific data
if result.combined_pvalues is not None:
    # Fisher/Brown results
    pval_df = result.combined_pvalues
    sig_genes = pval_df[pval_df['combined_padj'] < 0.05]
    
if result.ranked_genes is not None:
    # RRA results
    ranked_df = result.ranked_genes
    top_genes = ranked_df.head(100)
```

---

## Input Requirements

### DE Result Files

Must be pickle files (.pkl) from Module 7 containing `DEResult` objects.

**Required columns in DEResult.data:**
- `gene_id`: Gene identifier
- `pvalue`: Raw p-value
- `padj`: Adjusted p-value (FDR)
- `log2FoldChange`: Log2 fold change

**Expected structure:**
```python
de_result = pickle.load(f)
de_result.data  # pandas DataFrame with results
de_result.n_significant  # Number of significant genes
de_result.method  # Method name
```

### Minimum Requirements

- **At least 2 DE methods** are required for ensemble analysis
- **Recommended:** 3+ methods for robust consensus
- **Optimal:** 4-5 methods for comprehensive analysis

### File Format Validation

```python
# Check if files are valid
import pickle

with open('deseq2_result.pkl', 'rb') as f:
    result = pickle.load(f)
    
# Should have these attributes
assert hasattr(result, 'data')
assert hasattr(result, 'n_significant')
assert 'gene_id' in result.data.columns
assert 'pvalue' in result.data.columns
assert 'padj' in result.data.columns
assert 'log2FoldChange' in result.data.columns

print("✓ File is valid")
```

---

## Output Files

### Directory Structure

```
output_directory/
├── consensus_genes.csv              # Main output - consensus gene list
├── direction_consistency.csv        # Direction agreement per gene
├── statistics.json                  # Summary statistics
├── summary.txt                      # Human-readable summary
├── combined_pvalues.csv            # (Fisher/Brown only)
└── ranked_genes.csv                # (RRA only)
```

### File Descriptions

#### 1. `consensus_genes.csv`

**Description:** Final consensus gene list with all details

**Columns:**
- `gene_id`: Gene identifier
- `direction`: Consensus direction ('up' or 'down')
- `direction_agreement`: Fraction of methods agreeing
- `is_consistent`: Boolean - passed direction filter
- `n_methods_detected`: Number of methods detecting gene
- `methods_detected`: Comma-separated list of methods
- Method-specific columns (pvalue, padj, lfc per method)

**Example:**
```csv
gene_id,direction,direction_agreement,is_consistent,n_methods_detected,methods_detected,DESeq2_pvalue,edgeR_pvalue,limma_pvalue
ENSG00000001,up,1.0,True,3,"DESeq2,edgeR,limma",0.0001,0.0002,0.0003
ENSG00000002,down,1.0,True,3,"DESeq2,edgeR,limma",0.001,0.002,0.003
```

---

#### 2. `direction_consistency.csv`

**Description:** Per-gene direction consistency across methods

**Columns:**
- `gene_id`: Gene identifier
- `direction`: Consensus direction
- `direction_agreement`: Fraction of methods agreeing (0.0-1.0)
- `is_consistent`: Boolean
- `n_up`: Number of methods calling upregulated
- `n_down`: Number of methods calling downregulated
- `n_methods`: Total methods detecting gene

**Example:**
```csv
gene_id,direction,direction_agreement,is_consistent,n_up,n_down,n_methods
ENSG00000001,up,1.0,True,3,0,3
ENSG00000002,down,1.0,True,0,3,3
ENSG00000003,up,0.67,False,2,1,3
```

---

#### 3. `statistics.json`

**Description:** Comprehensive statistics in JSON format

**Structure:**
```json
{
  "overall": {
    "ensemble_method": "pvalue_combination_fisher",
    "n_methods": 3,
    "method_names": ["DESeq2", "edgeR", "limma"],
    "n_total_genes": 15000,
    "n_consensus_genes": 450,
    "n_upregulated": 250,
    "n_downregulated": 200,
    "n_direction_consistent": 440,
    "n_inconsistent": 10,
    "significance_threshold": 0.05,
    "timestamp": "2026-01-22T10:30:00",
    "raptor_version": "2.2.0"
  },
  "DESeq2": {
    "n_genes": 14500,
    "n_in_consensus": 445
  },
  "edgeR": {
    "n_genes": 14800,
    "n_in_consensus": 448
  },
  "limma": {
    "n_genes": 14200,
    "n_in_consensus": 440
  }
}
```

---

#### 4. `combined_pvalues.csv` (Fisher/Brown only)

**Description:** Combined p-values for all genes

**Columns:**
- `gene_id`: Gene identifier
- `combined_pvalue`: Combined raw p-value
- `combined_padj`: Combined adjusted p-value (FDR)
- `n_methods`: Number of methods with p-values
- Individual method p-values

**Example:**
```csv
gene_id,combined_pvalue,combined_padj,n_methods,DESeq2_pvalue,edgeR_pvalue,limma_pvalue
ENSG00000001,1.2e-08,3.4e-05,3,0.0001,0.0002,0.0003
ENSG00000002,3.5e-07,5.1e-04,3,0.001,0.002,0.003
```

---

#### 5. `ranked_genes.csv` (RRA only)

**Description:** Genes ranked by RRA score

**Columns:**
- `gene_id`: Gene identifier
- `rra_score`: RRA score (lower = better)
- `rra_pvalue`: RRA p-value
- `n_methods`: Number of methods detecting gene
- Individual method ranks

**Example:**
```csv
gene_id,rra_score,rra_pvalue,n_methods,DESeq2_rank,edgeR_rank,limma_rank
ENSG00000001,0.001,1.2e-05,3,1,2,1
ENSG00000002,0.002,3.4e-05,3,3,1,2
```

---

## Method Selection Guide

### Decision Tree

```
START: Choose ensemble method

Do you want MAXIMUM gene discovery?
├─ YES → Use Fisher's or Brown's method
│   └─ Are your methods correlated (same normalization)?
│       ├─ YES → Brown's method
│       └─ NO → Fisher's method
│
└─ NO → Do you want HIGH CONFIDENCE genes?
    ├─ YES → Use Voting
    │   └─ How many methods do you have?
    │       ├─ 3 methods → Use min_methods=3 (unanimous)
    │       └─ 4+ methods → Use min_methods=3 or 2
    │
    └─ NO → Want BALANCED approach?
        ├─ YES → Use RRA
        └─ NO → Do you have Module 8 optimization results?
            ├─ YES → Use Weighted ensemble
            └─ NO → Use RRA or Voting
```

### Use Case Recommendations

| Research Goal | Recommended Method | Typical Results |
|---------------|-------------------|----------------|
| Exploratory analysis, hypothesis generation | Fisher's or Brown's | 400-600 genes |
| Publication main results | Voting (≥3 methods) | 150-300 genes |
| Validation experiment planning | Voting (≥3 methods) | 150-300 genes |
| Balanced discovery + confidence | RRA | 300-500 genes |
| Quality-aware selection | Weighted | 250-400 genes |
| Meta-analysis across studies | Brown's | 400-600 genes |
| High-throughput validation | RRA or Weighted | 300-500 genes |

### Parameter Recommendations

| Scenario | Method | Threshold | Direction | Min Methods |
|----------|--------|-----------|-----------|-------------|
| **Strict (low FDR)** | Voting | 0.01 | 1.0 | 3 |
| **Standard** | RRA | 0.05 | 1.0 | 2 |
| **Exploratory** | Fisher's | 0.05 | 0.67 | 2 |
| **Very strict** | Voting | 0.01 | 1.0 | 4 |
| **Balanced** | Brown's | 0.05 | 1.0 | 2 |

---

## Complete Examples

### Example 1: Quick Start (Fisher's Method)

```python
import pickle
from raptor.ensemble import ensemble_fisher
from pathlib import Path

# Load DE results
with open('deseq2_result.pkl', 'rb') as f:
    deseq2 = pickle.load(f)
with open('edger_result.pkl', 'rb') as f:
    edger = pickle.load(f)
with open('limma_result.pkl', 'rb') as f:
    limma = pickle.load(f)

de_results = {
    'DESeq2': deseq2,
    'edgeR': edger,
    'limma': limma
}

# Run Fisher's method
result = ensemble_fisher(
    de_results=de_results,
    significance_threshold=0.05,
    output_dir='output/fisher/'
)

# View results
print(result.summary())
print(f"\nConsensus genes: {result.n_consensus_genes}")

# Get top genes
top_genes = result.consensus_genes.head(20)
print("\nTop 20 genes:")
print(top_genes[['gene_id', 'direction', 'combined_padj']])
```

---

### Example 2: Compare All Methods

```python
import pickle
from raptor.ensemble import (
    ensemble_fisher,
    ensemble_brown,
    ensemble_rra,
    ensemble_voting,
    ensemble_weighted
)

# Load results (same as Example 1)
de_results = {...}

# Run all methods
fisher_result = ensemble_fisher(de_results, significance_threshold=0.05)
brown_result = ensemble_brown(de_results, significance_threshold=0.05)
rra_result = ensemble_rra(de_results, significance_threshold=0.05)
voting_result = ensemble_voting(de_results, min_methods=2)
weighted_result = ensemble_weighted(de_results, min_score=1.5)

# Compare results
methods = {
    'Fisher': fisher_result,
    'Brown': brown_result,
    'RRA': rra_result,
    'Voting': voting_result,
    'Weighted': weighted_result
}

print("Method Comparison:")
print("-" * 50)
for name, result in methods.items():
    print(f"{name:12} {result.n_consensus_genes:6} genes")
```

---

### Example 3: Creating Tiered Confidence Levels

```python
import pickle
import pandas as pd
from raptor.ensemble import (
    ensemble_fisher,
    ensemble_rra,
    ensemble_voting
)

# Load results
de_results = {...}

# Run multiple methods
fisher = ensemble_fisher(de_results)
rra = ensemble_rra(de_results)
voting_2 = ensemble_voting(de_results, min_methods=2)
voting_3 = ensemble_voting(de_results, min_methods=3)

# Get gene sets
fisher_genes = set(fisher.consensus_genes['gene_id'])
rra_genes = set(rra.consensus_genes['gene_id'])
voting_2_genes = set(voting_2.consensus_genes['gene_id'])
voting_3_genes = set(voting_3.consensus_genes['gene_id'])

# Create tiers
tier1 = voting_3_genes  # Very high confidence
tier2 = voting_2_genes - voting_3_genes  # High confidence
tier3 = rra_genes - voting_2_genes  # Moderate confidence
tier4 = fisher_genes - rra_genes  # Exploratory

print(f"Tier 1 (Very High): {len(tier1)} genes - unanimous")
print(f"Tier 2 (High):      {len(tier2)} genes - majority")
print(f"Tier 3 (Moderate):  {len(tier3)} genes - RRA only")
print(f"Tier 4 (Explore):   {len(tier4)} genes - Fisher only")

# Save tiered list
tiered_df = pd.DataFrame({
    'gene_id': list(tier1) + list(tier2) + list(tier3),
    'tier': [1]*len(tier1) + [2]*len(tier2) + [3]*len(tier3),
    'confidence': ['Very High']*len(tier1) + ['High']*len(tier2) + ['Moderate']*len(tier3)
})

tiered_df.to_csv('tiered_genes.csv', index=False)
```

---

### Example 4: Weighted Ensemble with Module 8

```python
import pickle
import pandas as pd
from raptor.ensemble import ensemble_weighted

# Load DE results
de_results = {...}

# Load weights from Module 8 optimization
weights_df = pd.read_csv('output/module8/ensemble_weights.csv')
weights = dict(zip(weights_df['Method'], weights_df['Weight']))

print("Using optimized weights:")
for method, weight in weights.items():
    print(f"  {method}: {weight:.4f}")

# Run weighted ensemble
result = ensemble_weighted(
    de_results=de_results,
    weights=weights,
    min_score=1.5,  # Adjust based on weights
    check_direction=True,
    output_dir='output/weighted/'
)

print(f"\nWeighted consensus: {result.n_consensus_genes} genes")
```

---

### Example 5: Strict Validation Targets

```python
from raptor.ensemble import ensemble_voting

# Very strict parameters for validation
result = ensemble_voting(
    de_results=de_results,
    min_methods=3,  # All 3 must agree
    filters={
        'padj': 0.01,  # Very strict p-value
        'lfc': 1.5     # High fold change
    },
    check_direction=True,
    direction_threshold=1.0,  # Perfect agreement
    output_dir='output/validation_targets/'
)

print(f"Validation targets: {result.n_consensus_genes} genes")

# Get top candidates
candidates = result.consensus_genes.copy()

# Add additional filters
candidates = candidates[
    abs(candidates['meta_lfc']) > 2.0  # Very high fold change
]

print(f"High-priority candidates: {len(candidates)} genes")

# Save for qPCR validation
top_20 = candidates.nsmallest(20, 'combined_padj')
top_20[['gene_id', 'direction', 'meta_lfc', 'combined_padj']].to_csv(
    'qpcr_validation_list.csv',
    index=False
)
```

---

### Example 6: CLI Workflow

```bash
# Step 1: Run Module 7 to get DE results
raptor import-de \
    --deseq2 de_results/deseq2_results.csv \
    --edger de_results/edger_results.csv \
    --limma de_results/limma_results.csv \
    -o output/de_imported/

# Step 2: Run ensemble analysis (single method)
raptor ensemble \
    -m fisher \
    --deseq2 output/de_imported/deseq2_result.pkl \
    --edger output/de_imported/edger_result.pkl \
    --limma output/de_imported/limma_result.pkl \
    -o output/ensemble_fisher/

# Step 3: Compare all methods
raptor ensemble-compare \
    --deseq2 output/de_imported/deseq2_result.pkl \
    --edger output/de_imported/edger_result.pkl \
    --limma output/de_imported/limma_result.pkl \
    -o output/ensemble_compare/

# Step 4: Check results
ls output/ensemble_compare/
cat output/ensemble_compare/comparison_summary.csv
```

---

## Parameter Tuning

### Significance Threshold

**Default:** 0.05

**When to adjust:**
- **More strict (0.01)**: Fewer, higher-confidence genes
- **More lenient (0.1)**: More genes, exploratory analysis

```python
# Strict
result = ensemble_fisher(de_results, significance_threshold=0.01)

# Standard
result = ensemble_fisher(de_results, significance_threshold=0.05)

# Lenient
result = ensemble_fisher(de_results, significance_threshold=0.1)
```

---

### Direction Threshold

**Default:** 1.0 (all methods must agree)

**When to adjust:**
- **1.0**: Perfect agreement (strict)
- **0.75**: 3/4 methods agree
- **0.67**: 2/3 majority
- **0.5**: Simple majority

```python
# Perfect agreement (strict)
result = ensemble_fisher(de_results, direction_threshold=1.0)

# 2/3 majority (balanced)
result = ensemble_fisher(de_results, direction_threshold=0.67)

# Simple majority (lenient)
result = ensemble_fisher(de_results, direction_threshold=0.5)
```

**Recommendation:**
- 2-3 methods: Use 1.0 (perfect agreement)
- 4-5 methods: Can use 0.75 or 0.67
- Never go below 0.5

---

### Min Methods (Voting)

**Default:** 2

**When to adjust:**
- **3**: Very high confidence (strict)
- **2**: Balanced
- **1**: Not recommended (defeats ensemble purpose)

```python
# High confidence (3/3 methods)
result = ensemble_voting(de_results, min_methods=3)

# Balanced (2/3 methods)
result = ensemble_voting(de_results, min_methods=2)
```

---

### Check Direction

**Default:** True

**When to turn off:**
- Generally **not recommended**
- May be appropriate for exploratory analysis only

```python
# Standard (recommended)
result = ensemble_fisher(de_results, check_direction=True)

# Exploratory only (not recommended)
result = ensemble_fisher(de_results, check_direction=False)
```

---

### Use Adjusted P-values

**Default:** False (use raw p-values)

**Important:** For p-value combination methods (Fisher/Brown):
- **Always use raw p-values** (use_padj=False)
- Adjusted p-values have already been corrected
- Combining corrected p-values leads to over-correction

```python
# CORRECT (recommended)
result = ensemble_fisher(de_results, use_padj=False)

# INCORRECT (over-correction)
result = ensemble_fisher(de_results, use_padj=True)  # Don't do this!
```

For RRA:
- Either raw or adjusted is acceptable
- Raw p-values recommended for consistency

---

## Troubleshooting

### Common Issues

#### Issue 1: "Need at least 2 DE result files"

**Cause:** Only provided 1 DE result file

**Solution:**
```bash
# Wrong
raptor ensemble -m fisher --deseq2 d.pkl -o out/

# Correct
raptor ensemble -m fisher --deseq2 d.pkl --edger e.pkl -o out/
```

---

#### Issue 2: "Module 9 not available"

**Cause:** Import error, module not installed correctly

**Solution:**
```bash
# Check installation
python -c "from raptor.ensemble import ensemble_fisher"

# If error, reinstall
cd /path/to/raptor
pip install -e . --break-system-packages
```

---

#### Issue 3: "No genes in consensus"

**Cause:** Threshold too strict or direction filtering too stringent

**Solution:**
```python
# Relax parameters
result = ensemble_fisher(
    de_results,
    significance_threshold=0.1,  # More lenient
    direction_threshold=0.67     # 2/3 majority
)

# Check individual method results
for name, result in de_results.items():
    n_sig = (result.data['padj'] < 0.05).sum()
    print(f"{name}: {n_sig} significant genes")
```

---

#### Issue 4: "Direction inconsistent genes rejected"

**Cause:** Many genes have conflicting regulation direction

**Possible reasons:**
1. Methods using different normalizations
2. Different statistical models
3. Actual biology (context-dependent regulation)

**Solution:**
```python
# Check direction consistency
direction_df = result.direction_consistency
inconsistent = direction_df[~direction_df['is_consistent']]

print(f"Inconsistent genes: {len(inconsistent)}")
print("\nExamples:")
print(inconsistent.head()[['gene_id', 'n_up', 'n_down', 'direction_agreement']])

# Option 1: Relax threshold
result = ensemble_fisher(de_results, direction_threshold=0.67)

# Option 2: Investigate specific genes
gene = 'ENSG00000001'
for method, res in de_results.items():
    gene_data = res.data[res.data['gene_id'] == gene]
    if not gene_data.empty:
        lfc = gene_data.iloc[0]['log2FoldChange']
        padj = gene_data.iloc[0]['padj']
        print(f"{method}: LFC={lfc:.2f}, padj={padj:.4f}")
```

---

#### Issue 5: "Different gene counts than expected"

**Cause:** Understanding what each method reports

**Explanation:**
```python
# Fisher/Brown: Reports genes with combined_padj < threshold
# RRA: Reports genes with rra_pvalue < threshold
# Voting: Reports genes detected by ≥ min_methods
# Weighted: Reports genes with score ≥ min_score

# They will give different counts!
fisher_n = 450   # Maximum sensitivity
rra_n = 380      # Balanced
voting_n = 250   # High confidence
```

---

### Performance Issues

#### Issue: Slow processing

**Cause:** Large number of genes or methods

**Solution:**
```python
# Profile execution time
import time

start = time.time()
result = ensemble_fisher(de_results)
elapsed = time.time() - start

print(f"Execution time: {elapsed:.2f} seconds")

# For very large datasets (>50,000 genes):
# - Filter genes before ensemble (e.g., padj < 0.5)
# - Use RRA instead of Fisher/Brown
# - Process on HPC if available
```

---

### Data Quality Issues

#### Issue: Too many/few consensus genes

**Diagnosis:**
```python
# Check individual methods
for name, result in de_results.items():
    n_genes = len(result.data)
    n_sig_005 = (result.data['padj'] < 0.05).sum()
    n_sig_001 = (result.data['padj'] < 0.01).sum()
    
    print(f"\n{name}:")
    print(f"  Total genes: {n_genes}")
    print(f"  Significant (0.05): {n_sig_005}")
    print(f"  Significant (0.01): {n_sig_001}")
```

**Solutions:**

If too many genes (>1000):
```python
# Use stricter threshold
result = ensemble_fisher(de_results, significance_threshold=0.01)

# Or use voting instead
result = ensemble_voting(de_results, min_methods=3)

# Or add fold change filter
result = ensemble_voting(
    de_results,
    min_methods=2,
    filters={'padj': 0.05, 'lfc': 1.0}
)
```

If too few genes (<50):
```python
# Use more lenient threshold
result = ensemble_fisher(de_results, significance_threshold=0.1)

# Or reduce direction requirement
result = ensemble_fisher(de_results, direction_threshold=0.67)

# Or use Fisher instead of voting
result = ensemble_fisher(de_results)  # More sensitive
```

---

## Integration with Other Modules

### From Module 7 (DE Import)

**Module 7 → Module 9 workflow:**

```bash
# Step 1: Import DE results (Module 7)
raptor import-de \
    --deseq2 results/deseq2_results.csv \
    --edger results/edger_results.csv \
    --limma results/limma_results.csv \
    -o output/module7/

# Step 2: Ensemble analysis (Module 9)
raptor ensemble \
    -m fisher -m rra \
    --deseq2 output/module7/deseq2_result.pkl \
    --edger output/module7/edger_result.pkl \
    --limma output/module7/limma_result.pkl \
    -o output/module9/
```

**Python:**
```python
# Load Module 7 results
import pickle

de_results = {}
for method in ['deseq2', 'edger', 'limma']:
    with open(f'output/module7/{method}_result.pkl', 'rb') as f:
        de_results[method.capitalize()] = pickle.load(f)

# Run Module 9
from raptor.ensemble import ensemble_fisher
result = ensemble_fisher(de_results, output_dir='output/module9/')
```

---

### To Module 8 (Parameter Optimization)

**Module 8 optimization can inform Module 9:**

```python
# Load Module 8 optimization results
import pandas as pd

weights_df = pd.read_csv('output/module8/ensemble_weights.csv')
weights = dict(zip(weights_df['Method'], weights_df['Weight']))

# Use in weighted ensemble
from raptor.ensemble import ensemble_weighted

result = ensemble_weighted(
    de_results=de_results,
    weights=weights,  # From Module 8
    min_score=1.5,
    output_dir='output/module9_weighted/'
)
```

---

### Downstream Analysis

**Module 9 → Enrichment analysis:**

```python
# Get consensus genes for enrichment
consensus_genes = result.consensus_genes

# Export for enrichment tools
gene_list = consensus_genes['gene_id'].tolist()

# Save for enrichR, GSEA, etc.
with open('consensus_genes_for_enrichment.txt', 'w') as f:
    f.write('\n'.join(gene_list))

# Separate by direction
up_genes = consensus_genes[consensus_genes['direction'] == 'up']['gene_id']
down_genes = consensus_genes[consensus_genes['direction'] == 'down']['gene_id']

up_genes.to_csv('upregulated_genes.txt', index=False, header=False)
down_genes.to_csv('downregulated_genes.txt', index=False, header=False)
```

---

## Scientific References

### Primary Methods

1. **Fisher's Method**
   - Fisher, R.A. (1925). Statistical Methods for Research Workers. Oliver and Boyd, Edinburgh.
   - DOI: Classic statistical reference

2. **Brown's Method**
   - Brown, M.B. (1975). A method for combining non-independent, one-sided tests of significance. Biometrics, 31(4), 987-992.
   - DOI: 10.2307/2529826

3. **Robust Rank Aggregation**
   - Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012). Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics, 28(4), 573-580.
   - DOI: 10.1093/bioinformatics/btr709

### Meta-Analysis Reviews

4. **Hong & Breitling (2008)**
   - Hong, F., & Breitling, R. (2008). A comparison of meta-analysis methods for detecting differentially expressed genes in microarray experiments. Bioinformatics, 24(3), 374-382.
   - DOI: 10.1093/bioinformatics/btm620

5. **Ramasamy et al. (2008)**
   - Ramasamy, A., Mondry, A., Holmes, C.C., & Altman, D.G. (2008). Key issues in conducting a meta-analysis of gene expression microarray datasets. PLoS Medicine, 5(9), e184.
   - DOI: 10.1371/journal.pmed.0050184

### Ensemble Applications in Genomics

6. **Yang et al. (2020)**
   - Yang, X., et al. (2020). Concepts of artificial intelligence for computer-assisted drug discovery. Chemical Reviews, 119(18), 10520-10594.
   - Application of ensemble methods in computational biology

7. **Tseng et al. (2012)**
   - Tseng, G.C., Ghosh, D., & Feingold, E. (2012). Comprehensive literature review and statistical considerations for microarray meta-analysis. Nucleic Acids Research, 40(9), 3785-3799.
   - DOI: 10.1093/nar/gkr1265

---

## Quick Reference Card

### Essential Commands

```bash
# Run single method
raptor ensemble -m fisher --deseq2 d.pkl --edger e.pkl -o out/

# Compare all methods
raptor ensemble-compare --deseq2 d.pkl --edger e.pkl --limma l.pkl -o out/

# Custom threshold
raptor ensemble -m rra --deseq2 d.pkl --edger e.pkl --threshold 0.01 -o out/
```

### Essential Python

```python
from raptor.ensemble import ensemble_fisher, ensemble_rra, ensemble_voting

# Load results
de_results = {...}

# Quick analysis
result = ensemble_fisher(de_results)
print(f"Found {result.n_consensus_genes} genes")
result.save('output/')
```

### Default Parameters

| Parameter | Default | Alternative |
|-----------|---------|-------------|
| significance_threshold | 0.05 | 0.01, 0.1 |
| direction_threshold | 1.0 | 0.67, 0.75 |
| min_methods | 2 | 3, 4 |
| use_padj | False | True (not recommended) |
| check_direction | True | False |

### Method Decision Quick Guide

- **Exploratory** → Fisher's or Brown's
- **Validation targets** → Voting (min_methods=3)
- **Balanced** → RRA
- **Quality-aware** → Weighted (requires Module 8)

---

## Version History

**v2.2.0 (January 2026)**
- Complete implementation of all 5 ensemble methods
- Direction consistency checking
- Comprehensive CLI integration
- Full Python API
- Extensive documentation

---

## Support

For questions, issues, or feature requests:
- **Email:** ayehbolouki1988@gmail.com
- **GitHub:** github.com/ayehbolouki/RAPTOR
- **Documentation:** docs/MODULE9_ENSEMBLE_ANALYSIS.md

---

## Appendix: Statistical Details

### Fisher's Method Formula

```
Test statistic: X² = -2 Σ ln(p_i)  for i = 1 to k

Under H₀: X² ~ χ²(2k)

Combined p-value = P(X² ≥ observed | H₀)
```

### Brown's Method Correction

```
Scale factor: κ = E[X²] / (2k)

Corrected statistic: X²_c = X² / κ

Where E[X²] accounts for covariance between p-values
```

### RRA Score Calculation

```
For gene g with ranks r₁, r₂, ..., rₖ in k methods:

ρ_min = min(r₁/n₁, r₂/n₂, ..., rₖ/nₖ)

RRA score = Beta(ρ_min | 1, k)

Where Beta is the beta distribution CDF
```

---

**End of Guide**

This comprehensive guide covers all aspects of Module 9: Ensemble Analysis. For additional examples and use cases, see the example scripts in the `examples/` directory.
