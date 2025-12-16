# Tutorial 07: Threshold Optimization with ATO 

**Level**: Intermediate  
**Time**: 45-60 minutes  
**Goal**: Master data-driven threshold selection for differential expression analysis

---

## 🆕 New in RAPTOR v2.1.1!

This tutorial introduces the **Adaptive Threshold Optimizer (ATO)** - a data-driven approach to selecting significance thresholds that replaces arbitrary cutoffs like |logFC| > 1.

---

##  What You'll Learn

- ✅ Why arbitrary thresholds (|logFC| > 1) are problematic
- ✅ How ATO calculates optimal thresholds from YOUR data
- ✅ Choosing the right analysis goal (discovery/balanced/validation)
- ✅ Understanding p-value adjustment methods
- ✅ LogFC optimization techniques
- ✅ Generating publication-ready methods text
- ✅ Using ATO in the dashboard
- ✅ Applying ATO to ensemble analysis

---

## Prerequisites

- Basic understanding of differential expression analysis
- RAPTOR v2.1.1 installed
- DE results from any tool (DESeq2, edgeR, limma, etc.)

**No prior RAPTOR experience required!** This tutorial works with any DE results.

---

## Why Threshold Optimization?

### The Problem with Arbitrary Thresholds

Every RNA-seq publication uses thresholds like:
- "Genes with |log₂FC| > 1 and FDR < 0.05 were considered significant"

**But why |logFC| > 1?** 🤔

```
Researcher A: "I use |logFC| > 1 because... everyone does?"
Researcher B: "I use |logFC| > 0.58 because it's a 1.5-fold change"
Researcher C: "I use |logFC| > 2 because I want fewer genes"
Reviewer: "Why didn't you use |logFC| > 1.5?"
```

**The truth:** These are arbitrary choices that don't consider YOUR data!

### The Solution: Data-Driven Thresholds

```python
from raptor.threshold_optimizer import optimize_thresholds

# Let your DATA decide the threshold
result = optimize_thresholds(de_results, goal='balanced')

print(f"Optimal threshold: |logFC| > {result.logfc_threshold:.2f}")
# Output: Optimal threshold: |logFC| > 0.73

# This threshold is calculated from YOUR data distribution!
```

**Benefits:**
- ✅ Defensible in peer review
- ✅ Tailored to your specific dataset
- ✅ Publication-ready methods text included
- ✅ Statistically justified

---

## Step 1: Your First Threshold Optimization

### Load Your DE Results

```python
import pandas as pd
from raptor.threshold_optimizer import optimize_thresholds

# Load DE results from any tool
# DESeq2:
df = pd.read_csv('deseq2_results.csv')
# Columns: gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

# Or edgeR:
df = pd.read_csv('edger_results.csv')
# Columns: gene_id, logFC, logCPM, F, PValue, FDR

# Or limma:
df = pd.read_csv('limma_results.csv')
# Columns: gene_id, logFC, AveExpr, t, P.Value, adj.P.Val
```

### Run Basic Optimization

```python
# For DESeq2 results
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',  # Column name for log fold change
    pvalue_col='pvalue',          # Column name for p-values
    goal='balanced'               # Analysis goal
)

# View results
print(f"Optimal |logFC| threshold: {result.logfc_threshold:.3f}")
print(f"P-value threshold: {result.pvalue_threshold}")
print(f"Adjusted p-value threshold: {result.padj_threshold}")
print(f"Significant genes: {result.n_significant}")
print(f"Up-regulated: {result.n_up}")
print(f"Down-regulated: {result.n_down}")
```

**Example Output:**
```
Optimal |logFC| threshold: 0.731
P-value threshold: 0.05
Adjusted p-value threshold: 0.05
Significant genes: 1,247
Up-regulated: 623
Down-regulated: 624
```

### Get Publication Methods Text

```python
# This is the magic! 
print(result.methods_text)
```

**Output (copy directly to your paper!):**
```
Differential expression analysis was performed using DESeq2 v1.38.0. 
Significance thresholds were determined using the Adaptive Threshold 
Optimizer (ATO) with the 'balanced' analysis goal. The proportion of 
true null hypotheses (π₀) was estimated at 0.847 using Storey's spline 
method. An adjusted p-value threshold of 0.05 (Benjamini-Hochberg FDR 
correction) was applied. The log₂ fold change threshold was determined 
to be 0.731 using the MAD-based method, calculated from the data 
distribution rather than an arbitrary cutoff. This identified 1,247 
differentially expressed genes (623 up-regulated, 624 down-regulated).
```

**Compare to traditional methods text:**
```
"Genes with |log2FC| > 1 and FDR < 0.05 were considered significant."
```

Which would YOU trust more in peer review? 

---

## Step 2: Understanding Analysis Goals

ATO provides three analysis goals that balance sensitivity and specificity:

### Discovery Goal (Most Permissive)

```python
result = optimize_thresholds(df, goal='discovery')
```

**Use when:**
- Exploratory analysis
- Pilot studies
- Planning validation experiments
- Don't want to miss potential hits

**Effect:**
- Lower logFC threshold
- More genes pass filter
- Higher sensitivity, slightly lower specificity

### Balanced Goal (Recommended Default)

```python
result = optimize_thresholds(df, goal='balanced')
```

**Use when:**
- Standard publication
- Most typical analyses
- Want good balance of sensitivity/specificity

**Effect:**
- Standard FDR control
- Optimized logFC threshold
- Balanced trade-off

### Validation Goal (Most Stringent)

```python
result = optimize_thresholds(df, goal='validation')
```

**Use when:**
- Clinical applications
- Confirmation studies
- Expensive follow-up experiments
- Need high confidence

**Effect:**
- Higher logFC threshold
- Fewer genes pass filter
- Higher specificity, lower sensitivity

### Compare All Goals

```python
print("Goal\t\t|logFC|\tGenes")
print("-" * 35)
for goal in ['discovery', 'balanced', 'validation']:
    r = optimize_thresholds(df, goal=goal)
    print(f"{goal}\t\t{r.logfc_threshold:.2f}\t{r.n_significant}")
```

**Example Output:**
```
Goal            |logFC| Genes
-----------------------------------
discovery       0.52    1,876
balanced        0.73    1,247
validation      0.94    892
```

**Decision guide:**
- Exploratory? → Discovery
- Publication? → Balanced (default)
- Clinical/expensive follow-up? → Validation

---

## Step 3: P-Value Adjustment Methods

### Available Methods

```python
# Benjamini-Hochberg (default, recommended)
result = optimize_thresholds(df, padj_method='BH')

# Benjamini-Yekutieli (more conservative, for dependent tests)
result = optimize_thresholds(df, padj_method='BY')

# Storey's q-value (estimates π₀)
result = optimize_thresholds(df, padj_method='storey')

# Holm (step-down, family-wise error rate)
result = optimize_thresholds(df, padj_method='holm')

# Bonferroni (most conservative)
result = optimize_thresholds(df, padj_method='bonferroni')
```

### When to Use Each

| Method | Use When | Conservativeness |
|--------|----------|------------------|
| `BH` | Default choice, independent tests | Moderate |
| `storey` | Want π₀ estimation, large studies | Adaptive |
| `BY` | Correlated tests, unknown dependency | Conservative |
| `holm` | Need FWER control | Very conservative |
| `bonferroni` | Very few genes, strict control | Most conservative |

### Example: Comparing Methods

```python
methods = ['BH', 'storey', 'BY', 'holm', 'bonferroni']
print("Method\t\tSignificant Genes")
print("-" * 35)
for method in methods:
    r = optimize_thresholds(df, padj_method=method, goal='balanced')
    print(f"{method}\t\t{r.n_significant}")
```

---

## Step 4: LogFC Optimization Methods

### Available Methods

```python
# Auto (recommended - selects best method for your data)
result = optimize_thresholds(df, logfc_method='auto')

# MAD-based (robust to outliers)
result = optimize_thresholds(df, logfc_method='mad')

# Mixture model (models null + DE distributions)
result = optimize_thresholds(df, logfc_method='mixture')

# Power analysis (estimates detectable effect size)
result = optimize_thresholds(df, logfc_method='power')

# Percentile-based (simple, robust)
result = optimize_thresholds(df, logfc_method='percentile')
```

### Understanding Each Method

**MAD (Median Absolute Deviation):**
- Uses median-based statistics
- Robust to extreme values
- Good default choice
- Formula: `threshold = median + k × MAD`

**Mixture Model:**
- Fits two distributions (null + DE)
- More sophisticated
- Requires enough DE genes
- Good for large datasets

**Power Analysis:**
- Based on statistical power
- Considers your sample size
- Good when effect sizes matter
- Requires sample size information

**Percentile:**
- Simple approach
- Uses distribution percentiles
- Very robust
- Good for unusual distributions

### Auto Selection

The `auto` method intelligently selects:
```python
result = optimize_thresholds(df, logfc_method='auto')
print(f"Auto selected: {result.logfc_method_used}")
```

---

## Step 5: Working with Different DE Tools

### DESeq2

```python
# DESeq2 column names
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    padj_col='padj',  # Optional: use pre-computed adjusted p-values
    goal='balanced'
)
```

### edgeR

```python
# edgeR column names
result = optimize_thresholds(
    df,
    logfc_col='logFC',
    pvalue_col='PValue',
    padj_col='FDR',
    goal='balanced'
)
```

### limma

```python
# limma column names
result = optimize_thresholds(
    df,
    logfc_col='logFC',
    pvalue_col='P.Value',
    padj_col='adj.P.Val',
    goal='balanced'
)
```

### NOISeq (Special Case)

```python
# NOISeq uses probability, not p-values
# Convert probability to p-value-like metric
df['pvalue_like'] = 1 - df['prob']

result = optimize_thresholds(
    df,
    logfc_col='log2FC',
    pvalue_col='pvalue_like',
    goal='balanced'
)
```

---

## Step 6: Using ATO in the Dashboard

### Launch Dashboard

```bash
raptor dashboard
```

### Navigate to ATO Page

1. Open browser to `http://localhost:8501`
2. Click **"🎯 Threshold Optimizer"** in sidebar
3. Upload your DE results (CSV file)
4. Select your column names from dropdowns
5. Choose analysis goal
6. Click **"Optimize Thresholds"**

### Dashboard ATO Features

```
┌────────────────────────────────────────────────────┐
│  🎯 Threshold Optimizer                            │
├────────────────────────────────────────────────────┤
│                                                     │
│  📁 Upload DE Results: [Choose File]               │
│                                                     │
│  Column Mapping:                                   │
│  LogFC column: [log2FoldChange ▼]                 │
│  P-value column: [pvalue ▼]                       │
│                                                     │
│  Analysis Goal: ○ Discovery ● Balanced ○ Validation│
│                                                     │
│  [🎯 Optimize Thresholds]                          │
│                                                     │
├────────────────────────────────────────────────────┤
│  Results:                                           │
│  ✅ Optimal |logFC|: 0.731                         │
│  ✅ P-value threshold: 0.05                        │
│  ✅ Significant genes: 1,247                       │
│                                                     │
│  [📊 View Volcano Plot] [📝 Get Methods Text]      │
│  [💾 Download Results]                              │
└────────────────────────────────────────────────────┘
```

### Download Options

- **Optimized Results CSV**: DE results with `significant` column
- **Methods Text**: Copy-paste for your paper
- **Volcano Plot**: Publication-ready figure with optimized thresholds

---

## Step 7: ATO + Ensemble Analysis

### Problem: Different Thresholds Across Pipelines

```
Pipeline 1: |logFC| > 1.0 → 1,234 genes
Pipeline 3: |logFC| > 0.58 → 1,567 genes
Pipeline 5: |logFC| > 1.5 → 892 genes

How do you compare? 🤔
```

### Solution: Uniform Thresholds with ATO

```python
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

# 1. Combine results from all pipelines
analyzer = EnsembleAnalyzer()
combined = analyzer.combine(pipeline_results, method='weighted_average')

# 2. Determine optimal threshold from combined data
result = optimize_thresholds(
    combined['combined_de'],
    logfc_col='log2FC_consensus',
    pvalue_col='FDR_consensus',
    goal='balanced'
)

# 3. Apply uniform threshold
threshold = result.logfc_threshold
print(f"Uniform threshold for all pipelines: |logFC| > {threshold:.2f}")

# Now all comparisons use the same threshold!
```

### Ensemble + ATO Workflow

```bash
# Single command
raptor ensemble \
  --data data/ \
  --pipelines 1,3,5 \
  --use-ato \
  --ato-goal balanced \
  --output results/
```

---

## Step 8: Visualizing Optimized Results

### Volcano Plot with Optimized Thresholds

```python
from raptor.threshold_optimizer import optimize_thresholds, create_volcano_plot

# Optimize
result = optimize_thresholds(df, goal='balanced')

# Create publication-ready volcano plot
fig = create_volcano_plot(
    result.results_df,
    logfc_threshold=result.logfc_threshold,
    pvalue_threshold=result.padj_threshold,
    title="Differential Expression (ATO-Optimized Thresholds)"
)
fig.savefig('volcano_optimized.png', dpi=300)
```

### Compare Traditional vs Optimized

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Traditional arbitrary threshold
ax1 = create_volcano_plot(df, logfc_threshold=1.0, ax=axes[0])
axes[0].set_title('Traditional (|logFC| > 1.0)\n1,089 genes')

# ATO-optimized threshold
ax2 = create_volcano_plot(df, logfc_threshold=result.logfc_threshold, ax=axes[1])
axes[1].set_title(f'ATO-Optimized (|logFC| > {result.logfc_threshold:.2f})\n{result.n_significant} genes')

plt.savefig('threshold_comparison.png', dpi=300)
```

---

## Step 9: Saving and Exporting Results

### Save Optimized Results

```python
# Save full results with significance column
result.results_df.to_csv('de_results_optimized.csv', index=False)

# Save only significant genes
significant = result.results_df[result.results_df['significant']]
significant.to_csv('significant_genes.csv', index=False)

# Save methods text
with open('methods_text.txt', 'w') as f:
    f.write(result.methods_text)

# Save threshold summary
summary = {
    'logfc_threshold': result.logfc_threshold,
    'pvalue_threshold': result.pvalue_threshold,
    'padj_threshold': result.padj_threshold,
    'n_significant': result.n_significant,
    'n_up': result.n_up,
    'n_down': result.n_down,
    'pi0': result.pi0,
    'goal': 'balanced',
    'logfc_method': result.logfc_method_used,
    'padj_method': 'BH'
}

import json
with open('threshold_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
```

---

## Step 10: Complete Workflow Example

### Publication-Ready Analysis

```python
import pandas as pd
from raptor.threshold_optimizer import optimize_thresholds, create_volcano_plot

# 1. Load your DE results
df = pd.read_csv('deseq2_results.csv')
print(f"Loaded {len(df)} genes")

# 2. Optimize thresholds
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced',
    padj_method='BH',
    logfc_method='auto'
)

# 3. Print summary
print("\n" + "="*50)
print("🎯 THRESHOLD OPTIMIZATION RESULTS")
print("="*50)
print(f"Optimal |logFC| threshold: {result.logfc_threshold:.3f}")
print(f"P-value adjustment: Benjamini-Hochberg")
print(f"Adjusted p-value threshold: {result.padj_threshold}")
print(f"π₀ estimate: {result.pi0:.3f}")
print(f"\nSignificant genes: {result.n_significant}")
print(f"  - Up-regulated: {result.n_up}")
print(f"  - Down-regulated: {result.n_down}")

# 4. Generate publication materials
print("\n" + "="*50)
print("📝 METHODS TEXT (copy to your paper):")
print("="*50)
print(result.methods_text)

# 5. Save results
result.results_df.to_csv('final_de_results.csv', index=False)
print("\n✅ Results saved to: final_de_results.csv")

# 6. Create volcano plot
fig = create_volcano_plot(result.results_df, 
                          logfc_threshold=result.logfc_threshold,
                          pvalue_threshold=result.padj_threshold)
fig.savefig('Figure1_volcano.png', dpi=300, bbox_inches='tight')
print("✅ Volcano plot saved to: Figure1_volcano.png")
```

---

## Troubleshooting

### Problem: "Column not found"

```python
# Check your column names
print(df.columns.tolist())

# Use the correct names
result = optimize_thresholds(df, logfc_col='YOUR_COLUMN_NAME')
```

### Problem: "π₀ estimation failed"

This happens with unusual p-value distributions:

```python
# Try different method
result = optimize_thresholds(df, pi0_method='histogram')

# Or use default π₀
result = optimize_thresholds(df, pi0_method='fixed', pi0_value=0.9)
```

### Problem: Extreme threshold values

```python
# Check if threshold seems wrong
if result.logfc_threshold > 5 or result.logfc_threshold < 0.1:
    print("Warning: Unusual threshold - check data quality")
    
    # Try different logfc method
    result = optimize_thresholds(df, logfc_method='percentile')
```

### Problem: No significant genes

```python
if result.n_significant == 0:
    # Try discovery goal
    result = optimize_thresholds(df, goal='discovery')
    
    # Or check your data
    print(f"Min p-value: {df['pvalue'].min()}")
    print(f"LogFC range: {df['log2FoldChange'].min():.2f} to {df['log2FoldChange'].max():.2f}")
```

---

## Best Practices

### Do's ✅
- Always report the threshold AND how it was determined
- Include ATO methods text in publications
- Compare goals when unsure which to use
- Use uniform thresholds for ensemble analysis
- Save all parameters for reproducibility

### Don'ts ❌
- Don't use arbitrary thresholds without justification
- Don't change thresholds after seeing results (p-hacking!)
- Don't ignore the π₀ estimate - it tells you about your data
- Don't use validation goal for exploratory studies

---

## Summary

You've learned to:
- ✅ **Replace arbitrary thresholds** with data-driven ones
- ✅ **Choose appropriate analysis goals** for your study
- ✅ **Select p-value adjustment methods** based on your needs
- ✅ **Understand logFC optimization** techniques
- ✅ **Generate publication-ready methods text**
- ✅ **Use ATO in the dashboard**
- ✅ **Apply uniform thresholds to ensemble analysis**

---

## Next Steps

1. **Apply to your data**: Run ATO on your DE results today!
2. **Read**: [THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md) - Complete ATO guide
3. **Try**: [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md) - ATO + Ensemble
4. **Explore**: Dashboard ATO page for no-code threshold optimization

---

## Citation

When using ATO in your research:

```
Significance thresholds were determined using the Adaptive Threshold 
Optimizer (ATO) from RAPTOR v2.1.1 (Bolouki, 2025) with the 
'[goal]' analysis goal. The proportion of true null hypotheses (π₀) 
was estimated at [value] using [method]. An adjusted p-value threshold 
of [value] ([adjustment method] correction) and log₂ fold change 
threshold of [value] ([logfc method]) were applied, identifying 
[N] differentially expressed genes.
```

**Or simply use `result.methods_text` - it's generated automatically!**

---

**Tutorial by Ayeh Bolouki**  
For RAPTOR v2.1.1

*"Let your data decide the thresholds!"* 🎯🦖
