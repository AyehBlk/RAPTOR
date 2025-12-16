# ❓ RAPTOR v2.1.1 Frequently Asked Questions (FAQ)

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

Your questions answered about RAPTOR v2.1.1!

---

##  Table of Contents

1. [General Questions](#general-questions)
2. [Installation & Setup](#installation--setup)
3. [ML Recommendations](#ml-recommendations)
4. [🎯 Threshold Selection](#threshold-selection) - **NEW!**
5. [Dashboard](#dashboard)
6. [Pipeline Selection](#pipeline-selection)
7. [Resource Requirements](#resource-requirements)
8. [Cloud Computing](#cloud-computing)
9. [Results & Output](#results--output)
10. [Troubleshooting](#troubleshooting)
11. [Advanced Usage](#advanced-usage)

---

##  General Questions

### What is RAPTOR?

RAPTOR (RNA-seq Analysis Pipeline Testing and Optimization Resource) is a comprehensive framework for testing, comparing, and optimizing RNA-seq analysis pipelines. Version 2.1.1 adds **Adaptive Threshold Optimizer (ATO)** for data-driven significance thresholds.

**Key Features:**
- 🤖 ML-based pipeline recommendations (85-90% accuracy)
- 🎯 **Adaptive Threshold Optimizer** (NEW in v2.1.1!)
- 📊 Interactive web dashboard
- 📈 Real-time resource monitoring
- 🎭 Ensemble analysis across pipelines
- ☁️ Cloud computing support

---

### What's new in v2.1.1?

**🎯 Adaptive Threshold Optimizer (ATO)**

Stop using arbitrary thresholds! ATO provides:
- Data-driven logFC threshold selection
- Multiple p-value adjustment methods
- π₀ estimation for true null proportion
- Publication-ready methods text
- Interactive dashboard page

```python
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(de_results, goal='discovery')
print(f"Optimal |logFC| > {result.logfc_threshold:.2f}")
```

[See complete changelog](CHANGELOG.md)

---

### Is RAPTOR free?

**Yes!** RAPTOR is completely free and open-source under the MIT License.

- ✅ Free to use
- ✅ Free to modify
- ✅ Free for commercial use
- ✅ No registration required

---

### How do I cite RAPTOR?

```
Bolouki, A. (2025). RAPTOR v2.1.1: RNA-seq Analysis Pipeline 
Testing and Optimization Resource with Adaptive Threshold 
Optimization. GitHub. https://github.com/AyehBlk/RAPTOR
```

---

##  Installation & Setup

### What are the system requirements?

**Minimum:**
- Python 3.8+
- 8 GB RAM
- 50 GB free disk space

**Recommended:**
- Python 3.9+
- 32 GB RAM
- 500 GB free disk space

[Detailed requirements](INSTALLATION.md)

---

### How do I verify v2.1.1 installation?

```bash
# Check version
python -c "import raptor; print(raptor.__version__)"
# Should show: 2.1.1

# Verify Threshold Optimizer
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

---

### How long does installation take?

**Quick install (pip):** 5-10 minutes
```bash
pip install raptor-rnaseq
```

---

## 🤖 ML Recommendations

### How accurate are ML recommendations?

**Very accurate!**

- ✅ 85-90% accuracy overall
- ✅ 95%+ for common organisms
- ✅ 80-85% for novel conditions

---

### Can I trust ML over my expertise?

**Use ML as a guide, not gospel!**

Best practice: Use ML suggestion as starting point, then apply your expertise!

**NEW in v2.1.1:** After ML recommends a pipeline, use ATO for data-driven thresholds:
```python
# 1. Get ML recommendation
rec = recommender.recommend(profile)

# 2. Run pipeline
# ...

# 3. Optimize thresholds with ATO
result = optimize_thresholds(de_results, goal='balanced')
```

---

## 🎯 Threshold Selection (NEW!)

### Why shouldn't I use |logFC| > 1 and padj < 0.05?

**These are arbitrary!** The common thresholds have no scientific basis for YOUR data:

❌ **Problems with arbitrary thresholds:**
- One size doesn't fit all datasets
- May miss important genes (too stringent)
- May include noise (too permissive)
- Hard to justify to reviewers

✅ **ATO solution:**
- Data-driven threshold selection
- Optimized for YOUR specific data
- Scientifically justified
- Publication-ready methods text

```python
# Instead of arbitrary thresholds:
# bad: fdr < 0.05 & |logFC| > 1

# Use ATO:
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(df, goal='balanced')
# good: fdr < 0.05 & |logFC| > 0.73 (data-driven!)
```

---

### What is ATO?

**ATO (Adaptive Threshold Optimizer)** determines optimal significance thresholds based on YOUR data characteristics:

**What it does:**
1. Estimates π₀ (proportion of true nulls)
2. Calculates optimal logFC threshold
3. Applies appropriate p-value adjustment
4. Generates publication methods text

**Example:**
```python
from raptor.threshold_optimizer import optimize_thresholds

result = optimize_thresholds(de_results, goal='discovery')
print(f"Optimal logFC: {result.logfc_threshold:.2f}")
print(f"Significant genes: {result.n_significant}")
print(result.methods_text)  # Copy to paper!
```

---

### What analysis goal should I use?

| Goal | Use Case | Threshold Behavior |
|------|----------|-------------------|
| **discovery** | Exploratory analysis, hypothesis generation | More permissive, maximize sensitivity |
| **balanced** | Standard publication, typical analysis | Balanced sensitivity/specificity |
| **validation** | Clinical, confirmation studies | Stringent, maximize specificity |

```python
# Exploratory - find more genes
result = optimize_thresholds(df, goal='discovery')

# Standard publication
result = optimize_thresholds(df, goal='balanced')

# Confirmation study
result = optimize_thresholds(df, goal='validation')
```

---

### What p-value adjustment method should I use?

**Quick guide:**

| Method | When to Use |
|--------|-------------|
| **BH** (default) | Most analyses - good balance |
| **BY** | When tests are correlated |
| **storey** | Want more power with π₀ estimation |
| **holm/bonferroni** | Need FWER control (very stringent) |

```python
# Standard (recommended)
result = optimize_thresholds(df, padj_method='BH')

# More power
result = optimize_thresholds(df, padj_method='storey')

# Very stringent
result = optimize_thresholds(df, padj_method='bonferroni')
```

---

### What logFC method should I use?

**Use 'auto' (default)** - it takes consensus of all methods.

If you need specific behavior:

| Method | Best For |
|--------|----------|
| **auto** | Most cases (recommended) |
| **mad** | Data with outliers |
| **mixture** | Clean bimodal distribution |
| **power** | When you know desired power |
| **percentile** | Simple, distribution-based |

---

### How do I get methods text for my paper?

ATO automatically generates publication-ready text:

```python
result = optimize_thresholds(df, goal='balanced')

# Get methods text
print(result.methods_text)
```

**Example output:**
```
"Differential expression significance thresholds were determined 
using the Adaptive Threshold Optimizer (ATO) with the 'balanced' 
analysis goal. The proportion of true null hypotheses (π₀) was 
estimated at 0.85 using Storey's spline method. An adjusted 
p-value threshold of 0.05 (Benjamini-Hochberg FDR correction) 
and log₂ fold change threshold of 0.73 (determined by MAD-based 
robust estimation) were applied, identifying 1,247 differentially 
expressed genes (623 upregulated, 624 downregulated)."
```

Just copy this into your Methods section!

---

### Can I use ATO with any DE tool?

**Yes!** ATO works with output from:

- ✅ DESeq2
- ✅ edgeR
- ✅ limma
- ✅ NOISeq
- ✅ EBSeq
- ✅ Any tool with logFC and p-value columns

Just specify the correct column names:

```python
# DESeq2
result = optimize_thresholds(df, logfc_col='log2FoldChange', pvalue_col='pvalue')

# edgeR
result = optimize_thresholds(df, logfc_col='logFC', pvalue_col='PValue')

# limma
result = optimize_thresholds(df, logfc_col='logFC', pvalue_col='P.Value')
```

---

### Where is ATO in the dashboard?

Click **"🎯 Threshold Optimizer"** in the sidebar!

**Dashboard workflow:**
1. Launch: `raptor dashboard`
2. Click "🎯 Threshold Optimizer"
3. Upload your DE results CSV
4. Select analysis goal
5. Click "Optimize Thresholds"
6. Download results and methods text

---

## 📊 Dashboard

### Do I need coding to use the dashboard?

**No coding required!** Including the new Threshold Optimizer:

```bash
raptor dashboard
# Open http://localhost:8501
# Click "🎯 Threshold Optimizer"
# Upload, click, download!
```

---

### Is there a Threshold Optimizer page?

**Yes!** New in v2.1.1:

1. Launch dashboard: `raptor dashboard`
2. Click "🎯 Threshold Optimizer" in sidebar
3. Upload DE results
4. Get optimized thresholds
5. Download results and methods text

---

## 🔬 Pipeline Selection

### Which pipeline should I use?

**Use ML recommendations + ATO for best results:**

```bash
# 1. Get ML recommendation
raptor ml recommend --config config.yaml

# 2. Run recommended pipeline
raptor run --pipeline 3 ...

# 3. Optimize thresholds (NEW)
python -c "
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd
df = pd.read_csv('results.csv')
result = optimize_thresholds(df, goal='balanced')
print(f'Use |logFC| > {result.logfc_threshold:.2f}')
"
```

---

### Do all pipelines work with ATO?

**Yes!** ATO is pipeline-agnostic. It works with any DE results that have:
- Log fold change column
- P-value column

All 8 RAPTOR pipelines are compatible.

---

## 💾 Resource Requirements

### How much memory does ATO need?

**Very little!** ATO is lightweight:
- Runtime: <1 second for typical datasets
- Memory: ~10 MB for 20,000 genes
- No GPU required

---

## ☁️ Cloud Computing

### Can I use ATO in cloud deployments?

**Yes!** ATO works in all environments:
- Local machine
- HPC clusters
- AWS/GCP/Azure
- Docker containers

Container image updated to v2.1.1 with ATO included.

---

## 📈 Results & Output

### What outputs does ATO generate?

1. **ThresholdResult object** with:
   - `logfc_threshold` - Optimal |logFC| cutoff
   - `pvalue_threshold` - P-value cutoff
   - `n_significant` - Significant gene count
   - `methods_text` - Publication paragraph
   - `results_df` - Full results with flags

2. **Visualizations** (in dashboard):
   - Volcano plot with thresholds
   - P-value distribution
   - LogFC distribution
   - Threshold comparison heatmap

---

### How do I export ATO results?

```python
result = optimize_thresholds(df, goal='balanced')

# Save results
result.results_df.to_csv('optimized_results.csv')

# Save methods text
with open('methods.txt', 'w') as f:
    f.write(result.methods_text)
```

---

## 🔧 Troubleshooting

### ATO not found - what do I do?

```bash
# Update RAPTOR
pip install --upgrade raptor-rnaseq

# Verify
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('OK')"
```

---

### ATO gives strange thresholds

**Check your data:**

```python
# Verify columns exist
print(df.columns.tolist())

# Check data range
print(df['log2FoldChange'].describe())
print(df['pvalue'].describe())

# Try different goal
result = optimize_thresholds(df, goal='discovery')  # More permissive
```

---

### No significant genes found

**Try:**
1. Use 'discovery' goal (more permissive)
2. Check your data has real signal
3. Verify p-values aren't all 1.0

```python
# More permissive
result = optimize_thresholds(df, goal='discovery')

# Check data
print(f"P-values < 0.05: {(df['pvalue'] < 0.05).sum()}")
```

---

## 🔬 Advanced Usage

### Can I customize ATO behavior?

**Yes!** Full control over all parameters:

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer

ato = AdaptiveThresholdOptimizer(df, 'log2FoldChange', 'pvalue')

result = ato.optimize(
    goal='balanced',
    padj_method='storey',      # Custom p-value adjustment
    logfc_method='mad',         # Custom logFC method
    fdr_threshold=0.01          # Stricter FDR
)
```

---

### Can I use ATO in ensemble analysis?

**Yes!** Use ATO for uniform thresholds across pipelines:

```yaml
# config.yaml
threshold_optimizer:
  enabled: true
  ensemble_mode:
    uniform_thresholds: true
```

```python
# Or in Python
from raptor.threshold_optimizer import optimize_thresholds

# After ensemble combination
combined = ensemble_analyzer.combine(results)
thresholds = optimize_thresholds(combined, goal='balanced')
```

---

## 📚 Additional Resources

### Where can I learn more about ATO?

- [THRESHOLD_OPTIMIZER.md](THRESHOLD_OPTIMIZER.md) - Comprehensive guide
- [DASHBOARD.md](DASHBOARD.md) - Dashboard with ATO page
- [API.md](API.md) - Python API reference

---

## 🚀 Quick Reference

### Common Commands (Updated for v2.1.1)

```bash
# Install
pip install raptor-rnaseq

# Verify ATO
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('OK')"

# Get ML recommendation
raptor ml recommend --config config.yaml

# Run analysis
raptor analyze --config config.yaml

# Start dashboard (includes ATO page!)
raptor dashboard

# Quick threshold optimization
python -c "
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd
df = pd.read_csv('results.csv')
result = optimize_thresholds(df, goal='balanced')
print(f'|logFC| > {result.logfc_threshold:.2f}')
print(result.methods_text)
"
```

---

## 📞 Contact

**Author:** Ayeh Bolouki  
**Email:** ayehbolouki1988@gmail.com  
**GitHub:** https://github.com/AyehBlk/RAPTOR  
**Version:** 2.1.1  
**License:** MIT

---

**Still have questions?**  
Don't hesitate to reach out! We're here to help. 💙

*"The only stupid question is the one not asked!"* 🦖
