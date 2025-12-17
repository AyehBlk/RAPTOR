# Tutorial 07: Adaptive Threshold Optimizer (ATO) 🎯

**Level**: Intermediate  
**Time**: 45 minutes  
**Goal**: Master data-driven threshold optimization for publication-ready analysis

---

## 🆕 NEW in RAPTOR v2.1.1

The Adaptive Threshold Optimizer (ATO) replaces arbitrary thresholds with **statistically-grounded, data-driven cutoffs** for differential expression analysis.

---

##  What You'll Learn

- ✅ Why arbitrary thresholds are problematic
- ✅ How ATO determines optimal cutoffs
- ✅ Using different analysis goals (discovery/balanced/validation)
- ✅ Generating publication-ready methods text
- ✅ Applying ATO with different DE tools (DESeq2, edgeR, limma, NOISeq)
- ✅ Using ATO in the dashboard (no coding!)
- ✅ Combining ATO with ensemble analysis

---

##  Prerequisites

**Required**:
- Completed [Tutorial 01: Getting Started](tutorial_01_getting_started.md)
- RAPTOR v2.1.1 installed
- DE results from any pipeline (DESeq2, edgeR, limma, etc.)

**Verify ATO installation**:
```bash
raptor --version
# Should show: RAPTOR v2.1.1

python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

---

## 🤔 The Problem: Arbitrary Thresholds

### Traditional Approach

```python
# The classic filtering criteria
significant = results[
    (abs(results['log2FoldChange']) > 1.0) &  # Why 1.0?
    (results['padj'] < 0.05)                   # Why 0.05?
]
```

### Why This is Problematic

| Issue | Impact |
|-------|--------|
| **|logFC| > 1 is arbitrary** | Why not 0.5? Or 1.5? Or 0.8? |
| **One size doesn't fit all** | Different datasets need different thresholds |
| **No justification** | Reviewers ask: "Why did you choose 1.0?" |
| **Inconsistent literature** | Some papers use 1.0, others 1.5, others 0.585 |

### Real-World Example

```python
# Dataset A (high effect sizes):
# With |logFC| > 1.0 → 500 genes
# With |logFC| > 1.5 → 480 genes (almost the same!)
# Optimal might be 1.2

# Dataset B (subtle changes):
# With |logFC| > 1.0 → 50 genes
# With |logFC| > 0.5 → 800 genes (huge difference!)
# Optimal might be 0.6
```

---

## ✅ The Solution: Adaptive Threshold Optimizer

### How ATO Works

```
Your DE Results
      ↓
1. Estimate π₀ (proportion of true nulls)
      ↓
2. Analyze p-value distribution
      ↓
3. Calculate data-driven logFC threshold
      ↓
4. Apply appropriate FDR control
      ↓
Statistically-justified threshold + Methods text
```

### Key Statistical Concepts

**π₀ (pi-zero)**: The estimated proportion of genes that are NOT differentially expressed.
- High π₀ (e.g., 0.95): Most genes are null → more stringent threshold
- Low π₀ (e.g., 0.70): Many DE genes → can be more permissive

**Data-driven threshold**: Based on YOUR data's effect size distribution, not arbitrary convention.

---

##  Step 1: Basic ATO Usage

### Command Line

```bash
raptor optimize-thresholds \
    --results deseq2_results.csv \
    --logfc-col log2FoldChange \
    --pvalue-col pvalue \
    --goal balanced \
    --output threshold_results/
```

### Python API

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load your DE results
de_results = pd.read_csv('deseq2_results.csv')

# Run ATO
result = optimize_thresholds(
    de_results,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)

# View results
print(f"Optimal |logFC| threshold: {result.logfc_threshold:.3f}")
print(f"π₀ estimate: {result.pi0:.3f}")
print(f"Significant genes: {result.n_significant}")
print(f"\nMethods text:\n{result.methods_text}")
```

**Example output**:
```
Optimal |logFC| threshold: 0.847
π₀ estimate: 0.782
Significant genes: 623

Methods text:
Differential expression significance thresholds were determined 
using the Adaptive Threshold Optimizer (ATO) from RAPTOR v2.1.1...
```

---

## 🎯 Step 2: Analysis Goals

ATO supports three analysis goals to match your research needs:

###  Discovery Mode

**Use when**: Exploratory research, generating hypotheses, don't want to miss candidates

```python
result = optimize_thresholds(de_results, goal='discovery')
```

**Characteristics**:
- More permissive threshold
- Higher sensitivity (fewer false negatives)
- May include more false positives
- Good for: pilot studies, hypothesis generation

###  Balanced Mode (Recommended)

**Use when**: Standard analysis, most publications, general research

```python
result = optimize_thresholds(de_results, goal='balanced')
```

**Characteristics**:
- Standard FDR control
- Balance between sensitivity and specificity
- Most commonly appropriate choice
- Good for: most research, publications

###  Validation Mode

**Use when**: Clinical research, important decisions, need high confidence

```python
result = optimize_thresholds(de_results, goal='validation')
```

**Characteristics**:
- Stringent threshold
- Higher specificity (fewer false positives)
- May miss some true positives
- Good for: clinical studies, validation experiments

### Compare Goals

```python
for goal in ['discovery', 'balanced', 'validation']:
    result = optimize_thresholds(de_results, goal=goal)
    print(f"{goal:12} → |logFC| > {result.logfc_threshold:.3f} → {result.n_significant} genes")
```

**Example output**:
```
discovery    → |logFC| > 0.523 → 892 genes
balanced     → |logFC| > 0.847 → 623 genes
validation   → |logFC| > 1.156 → 412 genes
```

---

##  Step 3: Publication Methods Text

ATO automatically generates publication-ready methods text:

```python
result = optimize_thresholds(de_results, goal='balanced')

# Save methods text to file
with open('ato_methods.txt', 'w') as f:
    f.write(result.methods_text)
```

**Example methods text**:
```
Differential expression significance thresholds were determined 
using the Adaptive Threshold Optimizer (ATO) from RAPTOR v2.1.1. 
The proportion of true null hypotheses (π₀) was estimated at 0.782 
using the Storey method. Based on the 'balanced' analysis goal, 
an optimal |log₂FC| threshold of 0.847 was determined. This data-
driven threshold, combined with an FDR threshold of 0.050, yielded 
623 significant differentially expressed genes.
```

**Copy this directly to your paper!** No more "we used |logFC| > 1 because... reasons."

---

##  Step 4: Working with Different DE Tools

### DESeq2

```python
# DESeq2 column names
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue'
)
```

### edgeR

```python
# edgeR column names
result = optimize_thresholds(
    df,
    logfc_col='logFC',
    pvalue_col='PValue'
)
```

### limma

```python
# limma column names
result = optimize_thresholds(
    df,
    logfc_col='logFC',
    pvalue_col='P.Value'
)
```

### NOISeq

```python
# NOISeq uses probability instead of p-value
result = optimize_thresholds(
    df,
    logfc_col='log2FC',
    pvalue_col='prob',
    noiseq_mode=True  # Handles probability conversion
)
```

### Auto-Detect Columns

```python
# Let ATO guess the column names
result = optimize_thresholds(df)  # Auto-detects common column patterns
```

---

##  Step 5: Dashboard ATO (No Coding!)

For users who prefer graphical interfaces:

### Launch Dashboard

```bash
raptor dashboard
```

### Navigate to Threshold Optimizer

1. Click on **🎯 Threshold Optimizer** tab

2. **Upload DE Results**: Click "Choose File" → select your CSV

3. **Configure Columns**:
   - logFC column: Select from dropdown
   - p-value column: Select from dropdown

4. **Select Analysis Goal**:
   - Discovery / Balanced / Validation

5. **Click "Optimize Thresholds"**

6. **View Results**:
   - Optimal threshold
   - Number of significant genes
   - Methods text preview

7. **Download**:
   - Results CSV (with significance column)
   - Methods text (for your paper)

### Dashboard Screenshot

```
┌────────────────────────────────────────────────────────┐
│  🎯 Threshold Optimizer                                │
├────────────────────────────────────────────────────────┤
│                                                        │
│  📤 Upload: [deseq2_results.csv]                       │
│                                                        │
│  📊 Columns:                                           │
│     logFC:  [log2FoldChange ▼]                        │
│     p-value: [pvalue ▼]                               │
│                                                        │
│  🎯 Goal: ○ Discovery  ● Balanced  ○ Validation       │
│                                                        │
│  [🚀 Optimize]                                         │
│                                                        │
│  ════════════════════════════════════════════════════ │
│                                                        │
│  ✅ Results:                                           │
│     Optimal |logFC|: 0.847                            │
│     π₀: 0.782                                          │
│     Significant genes: 623                            │
│                                                        │
│  📝 Methods Text: [Preview] [Copy] [Download]         │
│                                                        │
│  [📥 Download Results CSV]                             │
│                                                        │
└────────────────────────────────────────────────────────┘
```

---

##  Step 6: ATO with Ensemble Analysis

For maximum robustness, combine ATO with ensemble analysis:

### Why ATO + Ensemble?

- **Fair comparison**: All pipelines use the same threshold
- **No pipeline advantage**: Threshold doesn't favor any method
- **Combined methods text**: Documents both approaches

### Command Line

```bash
raptor ensemble \
    --counts data.csv \
    --metadata meta.txt \
    --pipelines STAR-DESeq2,Salmon-DESeq2,Kallisto-DESeq2 \
    --use-ato \
    --ato-goal balanced \
    --output ensemble_results/
```

### Python API

```python
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Combine results from all pipelines
all_results = pd.concat([
    pd.read_csv('star_deseq2.csv'),
    pd.read_csv('salmon_deseq2.csv'),
    pd.read_csv('kallisto_deseq2.csv')
])

# Get uniform threshold from combined data
ato_result = optimize_thresholds(all_results, goal='balanced')
uniform_threshold = ato_result.logfc_threshold

# Run ensemble with uniform threshold
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(
    results_dir='results/',
    logfc_threshold=uniform_threshold
)

# Combined methods text
print(f"Pipeline comparison methods:\n{ensemble.methods_text}")
print(f"\nThreshold methods:\n{ato_result.methods_text}")
```

---

##  Best Practices

### ✅ Do's

1. **Always use ATO for publications**
   - Provides justification for thresholds
   - Auto-generated methods text
   - Data-driven approach

2. **Match goal to research stage**
   - Discovery → exploratory
   - Balanced → most research
   - Validation → clinical/important

3. **Document your choice**
   - Include methods text in paper
   - State which goal you used
   - Report π₀ estimate

4. **Use with ensemble for robustness**
   - Uniform threshold across pipelines
   - Fair comparison
   - Combined documentation

### ❌ Don'ts

1. **Don't mix ATO and arbitrary thresholds**
   - Be consistent within analysis
   - Don't use |logFC| > 1 for some, ATO for others

2. **Don't ignore π₀ estimate**
   - Very high π₀ (> 0.95): Few DE genes expected
   - Very low π₀ (< 0.50): Check data quality

3. **Don't skip the methods text**
   - It's free documentation
   - Reviewers appreciate it

---

##  Troubleshooting

### "Column not found"

```python
# Check your column names
print(df.columns.tolist())

# Specify exact names
result = optimize_thresholds(
    df,
    logfc_col='YOUR_EXACT_LOGFC_COLUMN',
    pvalue_col='YOUR_EXACT_PVALUE_COLUMN'
)
```

### "π₀ estimation failed"

```python
# Try alternative estimation method
result = optimize_thresholds(df, pi0_method='histogram')

# Or bootstrap method
result = optimize_thresholds(df, pi0_method='bootstrap')
```

### "Threshold seems too permissive/stringent"

```python
# Try different goals
for goal in ['discovery', 'balanced', 'validation']:
    result = optimize_thresholds(df, goal=goal)
    print(f"{goal}: {result.n_significant} genes")
    
# Choose the one that matches your research needs
```

### "Very few significant genes"

This may be biological! But check:
1. Data quality (run Tutorial 03 first)
2. Try 'discovery' goal
3. Check if π₀ is very high (few true DE genes)

---

##  Reference Documentation

- **[THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md)**: Complete ATO guide
- **[API.md](../API.md)**: Python API reference
- **[ENSEMBLE_GUIDE.md](../ENSEMBLE_GUIDE.md)**: ATO with ensemble
- **[FAQ.md](../FAQ.md)**: Threshold-related questions

---

##  Quick Reference

```bash
# Command line - basic
raptor optimize-thresholds --results degs.csv --goal balanced

# Command line - specify columns
raptor optimize-thresholds --results degs.csv \
    --logfc-col log2FoldChange --pvalue-col pvalue --goal balanced

# With ensemble
raptor ensemble --counts data.csv --use-ato --ato-goal balanced

# Python - basic
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(df, goal='balanced')
print(result.logfc_threshold)
print(result.methods_text)

# Dashboard - no code
raptor dashboard  # → 🎯 Threshold Optimizer tab
```

---

##  Tutorial Checklist

Did you:
- [ ] Understand why arbitrary thresholds are problematic
- [ ] Run basic ATO on test data
- [ ] Try different analysis goals
- [ ] Generate methods text for a publication
- [ ] Try ATO with your DE tool (DESeq2/edgeR/limma)
- [ ] Use the dashboard Threshold Optimizer tab
- [ ] Combine ATO with ensemble analysis
- [ ] Know where to find more help

**Congratulations! Your thresholds are now data-driven! 🎯🦖**

---

##  What's Next?

- **[Tutorial 01](tutorial_01_getting_started.md)**: Getting Started basics
- **[Tutorial 03](tutorial_03_real_data.md)**: Complete real-data workflow with ATO
- **[Tutorial 06](tutorial_06_ensemble.md)**: Advanced ensemble + ATO
- **[THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md)**: Deep dive into ATO

---

**Tutorial 07 - Adaptive Threshold Optimizer**  
*RAPTOR v2.1.1*  
Created by Ayeh Bolouki  
Last updated: December 2025

*"Your data deserves better than arbitrary thresholds!"* 🎯
