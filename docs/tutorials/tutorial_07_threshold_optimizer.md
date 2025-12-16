# Tutorial 07: Adaptive Threshold Optimization 🎯

**Level**: Intermediate | **Time**: 45 min | **Goal**: Master data-driven threshold selection

---

## 🆕 NEW in RAPTOR v2.1.1

Replace arbitrary |logFC| > 1 with statistically-grounded, data-driven thresholds.

---

## What You'll Learn

- ✅ Why arbitrary thresholds are problematic
- ✅ How ATO determines optimal cutoffs
- ✅ Choosing analysis goals (discovery/balanced/validation)
- ✅ Generating publication methods text
- ✅ Using ATO with ensemble analysis
- ✅ Using ATO in the dashboard

---

## The Problem: Arbitrary Thresholds

```python
# Why 1.0? Why 0.05? 🤷
significant = results[(abs(results['logFC']) > 1.0) & (results['padj'] < 0.05)]
```

**Issues:**
- |logFC| > 1 is arbitrary - why not 0.5 or 1.5?
- Different datasets need different thresholds
- Reviewers increasingly ask for justification

---

## The Solution: ATO

```python
from raptor.threshold_optimizer import optimize_thresholds

result = optimize_thresholds(de_results, goal='balanced')
# ATO analyzes YOUR data to determine optimal threshold!
```

---

## Step 1: Basic Usage

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

de_results = pd.read_csv('deseq2_results.csv')

result = optimize_thresholds(
    de_results,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)

print(f"Optimal |logFC|: {result.logfc_threshold:.3f}")
print(f"Significant genes: {result.n_significant}")
```

---

## Step 2: Analysis Goals

###  Discovery
```python
result = optimize_thresholds(de_results, goal='discovery')
# More permissive - catch more candidates
```

###  Balanced (Recommended)
```python
result = optimize_thresholds(de_results, goal='balanced')
# Standard FDR control
```

###  Validation
```python
result = optimize_thresholds(de_results, goal='validation')
# Stringent - fewer false positives
```

---

## Step 3: Publication Methods Text

```python
result = optimize_thresholds(de_results, goal='balanced')
print(result.methods_text)  # Copy to your paper!
```

**Output:**
```
Differential expression significance thresholds were determined 
using the Adaptive Threshold Optimizer (ATO) from RAPTOR v2.1.1. 
The proportion of true null hypotheses (π₀) was estimated at 0.847...
```

---

## Step 4: Different DE Tools

### DESeq2
```python
result = optimize_thresholds(df, logfc_col='log2FoldChange', pvalue_col='pvalue')
```

### edgeR
```python
result = optimize_thresholds(df, logfc_col='logFC', pvalue_col='PValue')
```

### limma
```python
result = optimize_thresholds(df, logfc_col='logFC', pvalue_col='P.Value')
```

---

## Step 5: Dashboard ATO

```bash
raptor dashboard
# Navigate to "🎯 Threshold Optimizer" tab
# No coding required!
```

---

## Step 6: ATO with Ensemble

```python
# Apply uniform threshold across pipelines
combined = pd.concat([p1_results, p2_results, p3_results])
result = optimize_thresholds(combined, goal='balanced')
uniform_threshold = result.logfc_threshold

# Apply to all pipelines
for df in [p1_results, p2_results, p3_results]:
    df['significant'] = (abs(df['logFC']) > uniform_threshold) & (df['padj'] < 0.05)
```

Or use:
```bash
raptor ensemble --pipelines 1,3,4 --use-ato --ato-goal balanced
```

---

## Best Practices

### Do's ✅
✅ Use ATO for every publication  
✅ Include methods text  
✅ Match goal to research needs  
✅ Save thresholds for reproducibility  

### Don'ts ❌
❌ Use arbitrary |logFC| > 1 without justification  
❌ Mix ATO and arbitrary in same analysis  
❌ Ignore π₀ estimate  

---

## Troubleshooting

**Column not found:**
```python
print(df.columns.tolist())  # Check column names
```

**π₀ estimation failed:**
```python
result = optimize_thresholds(df, pi0_method='histogram')
```

---

## Summary

- ✅ Replace arbitrary thresholds with data-driven values
- ✅ Generate publication methods text automatically
- ✅ Use goals: discovery/balanced/validation
- ✅ Works with DESeq2, edgeR, limma
- ✅ Dashboard ATO tab for no-code usage

---

## Next Steps

- **[THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md)** - Complete reference
- **[Tutorial 06](tutorial_06_ensemble.md)** - ATO + ensemble

---

*RAPTOR v2.1.1 | Ayeh Bolouki | December 2025*

*"Your data deserves better than arbitrary thresholds!"* 🎯
