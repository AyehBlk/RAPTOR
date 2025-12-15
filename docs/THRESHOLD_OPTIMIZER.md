# Adaptive Threshold Optimizer (ATO)

## Overview

The **Adaptive Threshold Optimizer (ATO)** is a data-driven module for optimizing significance thresholds in RNA-seq differential expression analysis. It replaces arbitrary cutoffs (|logFC| > 1, padj < 0.05) with scientifically justified, dataset-specific values.

## The Problem

Most RNA-seq analyses use arbitrary thresholds:
- **|log2FoldChange| > 1**: Why 1? Why not 0.5 or 1.5?
- **padj < 0.05**: Convention, not optimized for your data
- **BH adjustment**: May not be optimal for your correlation structure

These choices are rarely justified and can miss real biological signals or include false positives.

## The Solution

ATO analyzes your DE results and recommends:
1. **Optimal p-value adjustment method** based on correlation structure and π₀
2. **Data-driven logFC cutoff** based on null distribution
3. **Evidence-based reasoning** for all recommendations

---

## Quick Start

### Installation

The ATO module is included in RAPTOR v2.1.1+:

```bash
pip install raptor-rnaseq>=2.1.1
```

### Basic Usage

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load your DESeq2/edgeR/limma results
df = pd.read_csv('deseq2_results.csv', index_col=0)

# Optimize thresholds
result = optimize_thresholds(df, goal='discovery')

# View results
print(result.summary())

# Get significant genes with optimized thresholds
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
ato = AdaptiveThresholdOptimizer(df)
result = ato.optimize()
significant_genes = ato.get_significant_genes()
```

### Dashboard

Access ATO through the RAPTOR dashboard:

```bash
streamlit run dashboard.py
```

Navigate to **🎯 Threshold Optimizer** in the sidebar.

---

## Input Requirements

ATO accepts DE results from any tool (DESeq2, edgeR, limma-voom, etc.):

### Required Columns

| Column | Aliases | Description |
|--------|---------|-------------|
| `log2FoldChange` | `logFC`, `log2FC`, `lfc` | Log2 fold change |
| `pvalue` | `PValue`, `P.Value`, `pval` | Raw p-values |

### Optional Columns

| Column | Aliases | Description |
|--------|---------|-------------|
| `padj` | `adj.P.Val`, `FDR`, `qvalue` | Pre-calculated adjusted p-values |
| `baseMean` | `AveExpr`, `logCPM` | Mean expression |
| `lfcSE` | `SE` | LogFC standard error |

---

## Analysis Goals

Set `goal` parameter based on your analysis purpose:

### `goal='discovery'` (Default)
- **Use case**: Exploratory analysis, hypothesis generation
- **Error control**: FDR (False Discovery Rate)
- **Behavior**: More lenient, maximizes true positives
- **Typical output**: More DE genes

### `goal='validation'`
- **Use case**: Confirming known targets, clinical validation
- **Error control**: FWER (Family-Wise Error Rate)
- **Behavior**: Conservative, minimizes false positives
- **Typical output**: Fewer DE genes, higher confidence

### `goal='balanced'`
- **Use case**: Publication, balanced analysis
- **Error control**: Standard FDR with moderate stringency
- **Behavior**: Balance between sensitivity and specificity

---

## P-Value Adjustment Methods

ATO supports and auto-selects from these methods:

### FDR Control Methods

#### Benjamini-Hochberg (BH)
- **When used**: Standard choice, valid for independent or PRDS tests
- **Formula**: `p_adj[i] = p[i] × n / rank[i]`
- **Reference**: Benjamini & Hochberg (1995)

#### Benjamini-Yekutieli (BY)
- **When used**: Negative correlations detected between genes
- **Formula**: `p_adj = BH × c(n)` where `c(n) = Σ(1/i) ≈ ln(n) + 0.577`
- **Reference**: Benjamini & Yekutieli (2001)

#### Storey's q-value
- **When used**: Many DE genes detected (π₀ < 0.8)
- **Advantage**: Estimates π₀ for ~20% more power than BH
- **Reference**: Storey (2002), Storey & Tibshirani (2003)

### FWER Control Methods

#### Holm (Step-Down)
- **When used**: `goal='validation'`
- **Advantage**: More powerful than Bonferroni, valid for any dependence
- **Reference**: Holm (1979)

#### Hochberg (Step-Up)
- **When used**: PRDS structure detected
- **Advantage**: More powerful than Holm
- **Reference**: Hochberg (1988)

#### Bonferroni
- **When used**: Very conservative validation
- **Formula**: `p_adj = p × n`
- **Note**: Most conservative, rarely recommended

---

## LogFC Cutoff Methods

### MAD-Based (Default)
Uses Median Absolute Deviation for robust null estimation:

```
1. Identify null genes (padj > 0.5)
2. Calculate MAD of logFC in null set
3. robust_σ = 1.4826 × MAD
4. cutoff = k × robust_σ (k = 2.5)
```

**Pros**: Robust to outliers, fast
**Reference**: Standard robust statistics

### Mixture Model
Fits 2-component Gaussian mixture:

```
Component 1 (Null): μ ≈ 0, small σ
Component 2 (DE): spread away from 0

Cutoff = point where P(null) < 0.05
```

**Pros**: Models biological reality
**Cons**: May not converge for weak signals

### Power-Based
Calculates minimum detectable effect:

```
min_logFC = (z_α + z_β) × SE(logFC)

Where:
- z_α = 1.96 (for α = 0.05, two-sided)
- z_β = 0.84 (for 80% power)
- SE estimated from data
```

**Pros**: Statistically principled
**Requires**: Sample size information

### Percentile-Based
Simple 95th percentile of null distribution:

```
cutoff = 95th percentile of |logFC| for genes with padj > 0.5
```

**Pros**: Simple, interpretable
**Cons**: Sensitive to null set definition

### Consensus (Auto)
Runs all methods and takes median:

```python
cutoff = median(MAD, Mixture, Power, Percentile)
```

**Recommended** for most analyses.

---

## Key Statistics

### π₀ (Pi-Zero)
The estimated proportion of true null hypotheses:

- **π₀ = 1.0**: All genes are null (no DE genes)
- **π₀ = 0.5**: Half of genes are truly DE
- **π₀ = 0.1**: 90% of genes are DE (rare)

**Interpretation**:
- π₀ > 0.9: Weak signal, use standard BH
- π₀ < 0.8: Strong signal, consider q-value
- π₀ < 0.5: Very strong signal, check data quality

### Impact Metrics
ATO reports:
- **n_significant_optimized**: DE genes with optimized thresholds
- **n_significant_traditional**: DE genes with |logFC|>1, padj<0.05
- **Difference**: Genes you would have missed/incorrectly included

---

## API Reference

### AdaptiveThresholdOptimizer

```python
class AdaptiveThresholdOptimizer:
    def __init__(
        self,
        df: pd.DataFrame,
        goal: str = 'discovery',  # 'discovery', 'validation', 'balanced'
        verbose: bool = True
    )
    
    def optimize(
        self,
        logfc_method: str = 'auto',  # 'auto', 'mad', 'mixture', 'power', 'percentile'
        n1: int = 3,  # Sample size group 1
        n2: int = 3   # Sample size group 2
    ) -> ThresholdResult
    
    def get_significant_genes(
        self,
        logfc_cutoff: float = None,
        padj_cutoff: float = None,
        use_optimized: bool = True
    ) -> pd.DataFrame
    
    def compare_thresholds(
        self,
        logfc_values: List[float] = [0.5, 1.0, 1.5, 2.0],
        padj_values: List[float] = [0.01, 0.05, 0.1]
    ) -> pd.DataFrame
    
    def get_adjustment_comparison(self) -> pd.DataFrame
```

### ThresholdResult

```python
@dataclass
class ThresholdResult:
    logfc_cutoff: float
    padj_cutoff: float
    padj_method: str
    logfc_method: str
    logfc_reasoning: str
    padj_reasoning: str
    n_significant_optimized: int
    n_significant_traditional: int
    pi0_estimate: float
    metrics: Dict[str, Any]
    
    def summary(self) -> str
```

### Convenience Function

```python
def optimize_thresholds(
    df: pd.DataFrame,
    goal: str = 'discovery',
    logfc_method: str = 'auto',
    n1: int = 3,
    n2: int = 3,
    verbose: bool = True
) -> ThresholdResult
```

---

## Visualization Functions

```python
from raptor.threshold_optimizer import (
    plot_optimization_summary,
    plot_volcano,
    plot_logfc_distribution,
    plot_pvalue_distribution,
    plot_threshold_comparison,
    plot_adjustment_comparison
)

# Comprehensive 4-panel summary
fig = plot_optimization_summary(result, df, save_path='summary.png')

# Volcano plot with optimized thresholds
fig = plot_volcano(df, logfc_cutoff=0.5, padj_cutoff=0.05)

# Threshold comparison heatmap
fig = plot_threshold_comparison(df, optimized_logfc=0.5, optimized_padj=0.05)
```

---

## Best Practices

### 1. Always Report Your Thresholds
```
"DE genes were identified using data-driven thresholds 
(|log2FC| > 0.45, q-value < 0.05) determined by the 
Adaptive Threshold Optimizer (RAPTOR v2.1.1), which 
estimated π₀ = 0.82 and recommended the Storey q-value 
method for p-value adjustment."
```

### 2. Compare with Traditional Thresholds
Always report both optimized and traditional results for transparency.

### 3. Consider Your Goal
- Exploratory → `goal='discovery'`
- Validation → `goal='validation'`
- Publication → `goal='balanced'`

### 4. Provide Sample Sizes
For accurate power-based cutoffs:
```python
result = ato.optimize(n1=4, n2=4)  # 4 samples per group
```

### 5. Examine the Diagnostics
Always check:
- π₀ estimate (reasonable? 0.5-0.95)
- LogFC distribution (bimodal expected for good data)
- P-value histogram (should be uniform + spike near 0)

---

## Troubleshooting

### "Not enough null genes"
**Cause**: padj > 0.5 yields < 100 genes
**Solution**: Strong signal in your data, results still valid

### "Mixture model did not converge"
**Cause**: Weak signal or unusual distribution
**Solution**: Use `logfc_method='mad'` instead

### "π₀ estimate is 1.0"
**Cause**: No significant genes detected
**Solution**: Check data quality, consider less stringent analysis

### "Negative logFC cutoff"
**Cause**: Numerical issue with MAD
**Solution**: Use `logfc_method='percentile'`

---

## Scientific References

### P-Value Adjustment
1. Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS-B* 57:289-300.

2. Benjamini Y, Yekutieli D (2001). The control of the false discovery rate in multiple testing under dependency. *Ann Stat* 29:1165-1188.

3. Holm S (1979). A simple sequentially rejective multiple test procedure. *Scand J Stat* 6:65-70.

4. Hochberg Y (1988). A sharper Bonferroni procedure for multiple tests of significance. *Biometrika* 75:800-802.

5. Storey JD (2002). A direct approach to false discovery rates. *JRSS-B* 64:479-498.

6. Storey JD, Tibshirani R (2003). Statistical significance for genomewide studies. *PNAS* 100:9440-9445.

### Effect Size Thresholds
7. McCarthy DJ, Smyth GK (2009). Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25:765-771.

---


---

## Citation

If you use ATO in your research, please cite:

```bibtex
@software{raptor2024,
  author = {Bolouki, Ayeh},
  title = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year = {2025},
  url = {https://github.com/AyehBlk/RAPTOR}
}
```
