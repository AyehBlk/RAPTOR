# RAPTOR Migration Guide: v2.1.0 → v2.1.1

**Complete guide for upgrading to RAPTOR v2.1.1**

---

## 🌟 Overview

RAPTOR v2.1.1 is a **feature release** that adds the **Adaptive Threshold Optimizer (ATO)** while maintaining **100% backward compatibility** with v2.1.0 and v2.0.0.

### 🎯 What's New in v2.1.1

**One Major Feature:**
1. 🎯 **Adaptive Threshold Optimizer (ATO)** - Data-driven threshold selection

**Key Benefits:**
- Replace arbitrary |logFC| > 1 with optimized values
- Multiple p-value adjustment methods
- π₀ estimation for true null proportion
- Auto-generated publication methods text
- Dashboard integration

---

## ✅ Is Migration Right for You?

### Migrate to v2.1.1 if:
- ✅ You want data-driven thresholds
- ✅ You're preparing manuscripts (need methods text)
- ✅ You want defensible threshold choices
- ✅ You use ensemble analysis
- ✅ You want the dashboard ATO page

### Stay on v2.1.0 if:
- ⚠️ You're in the middle of an analysis (wait until complete)
- ⚠️ Your current thresholds are validated
- ⚠️ You need exact result reproducibility

**Good news:** v2.1.1 is **100% backward compatible!**

---

## 📦 Installation

### Option 1: Upgrade Existing Installation (Recommended)

```bash
# If installed via pip
pip install --upgrade raptor-rnaseq

# Verify version
raptor --version  # Should show v2.1.1

# Verify ATO
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

### Option 2: Fresh Installation

```bash
# Create new environment
conda create -n raptor-v2.1.1 python=3.9
conda activate raptor-v2.1.1

# Install RAPTOR v2.1.1
pip install raptor-rnaseq

# Verify
raptor --version
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

### Option 3: Keep Both Versions

```bash
# Keep v2.1.0 environment
conda activate raptor-v2.1.0

# Create separate v2.1.1 environment
conda create -n raptor-v2.1.1 python=3.9
conda activate raptor-v2.1.1
pip install raptor-rnaseq

# Switch as needed
conda activate raptor-v2.1.0  # For reproducibility
conda activate raptor-v2.1.1  # For new ATO features
```

---

## 📚 New Dependencies

### No New Required Dependencies!

ATO uses existing scipy, numpy, and pandas.

### Optional Enhancement

```bash
# For advanced π₀ estimation (optional)
pip install statsmodels>=0.14.0
```

---

## ⚙️ Configuration Files

### v2.1.0 Config Files Still Work!

Your existing configs are **100% compatible**:

```bash
raptor profile \
  --counts data.csv \
  --config my_v2.1.0_config.yaml  # ✓ Still works!
```

### New v2.1.1 Config Options

Add these sections to enable ATO:

```yaml
# Add to your existing config.yaml

# NEW: Threshold Optimizer settings
threshold_optimizer:
  enabled: true
  goal: "balanced"           # discovery, balanced, or validation
  default_padj_method: "BH"  # BH, BY, storey, holm, bonferroni
  default_logfc_method: "auto"  # auto, mad, mixture, power, percentile

# Update statistics section
statistics:
  use_adaptive_thresholds: true  # NEW: Use ATO instead of fixed thresholds
  fdr_threshold: 0.05
  log2fc_threshold: 1.0  # Fallback if ATO disabled
```

### Example Config Files

```bash
# Copy example config
cp config/examples/config_with_ato.yaml my_config.yaml
```

---

## 🔄 Command Changes

### Good News: No Breaking Changes!

All v2.1.0 commands work identically:

```bash
# These all work exactly the same
raptor profile --counts data.csv --use-ml
raptor run --pipeline 3 --data fastq/
raptor ensemble --pipelines 1,3,4
raptor dashboard
```

### New Commands in v2.1.1

```bash
# NEW: Optimize thresholds from command line
raptor optimize-thresholds \
  --input de_results.csv \
  --goal balanced \
  --output optimized_results.csv

# NEW: Ensemble with ATO
raptor ensemble \
  --data data/ \
  --pipelines 1,3,4 \
  --use-ato \
  --ato-goal balanced
```

---

## 🔬 Migration Scenarios

### Scenario 1: Basic User (Threshold Optimization)

**v2.1.0 workflow:**
```bash
raptor run --pipeline 3 --data fastq/ --output results/
# Then manually filter: |logFC| > 1, FDR < 0.05 (arbitrary!)
```

**v2.1.1 workflow (enhanced):**
```bash
# Run pipeline (same)
raptor run --pipeline 3 --data fastq/ --output results/

# NEW: Optimize thresholds
python -c "
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd
df = pd.read_csv('results/de_results.csv')
result = optimize_thresholds(df, goal='balanced')
print(f'Optimal threshold: |logFC| > {result.logfc_threshold:.2f}')
print(result.methods_text)  # Copy to paper!
result.results_df.to_csv('results/optimized_results.csv')
"
```

**Migration effort:** 5 minutes to add ATO step

### Scenario 2: Dashboard User

**v2.1.0 workflow:**
```bash
raptor dashboard
# Use ML recommendations, view results
```

**v2.1.1 workflow (enhanced):**
```bash
raptor dashboard
# Same as before, PLUS:
# Click "🎯 Threshold Optimizer" in sidebar (NEW!)
# Upload DE results, get optimized thresholds
# Download methods text for publication
```

**Migration effort:** None (new page appears automatically)

### Scenario 3: Ensemble User

**v2.1.0 workflow:**
```bash
raptor ensemble --pipelines 1,3,4 --output results/
# Results use fixed thresholds per pipeline
```

**v2.1.1 workflow (enhanced):**
```bash
raptor ensemble \
  --pipelines 1,3,4 \
  --use-ato \
  --ato-goal balanced \
  --output results/
# Results use uniform, data-driven thresholds!
```

**Migration effort:** Add `--use-ato` flag

### Scenario 4: Python API User

**v2.1.0 code:**
```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler

profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile()

recommender = MLPipelineRecommender()
recommendation = recommender.recommend(profile)
```

**v2.1.1 code (backward compatible + enhanced):**
```python
# Your v2.1.0 code still works!
from raptor import MLPipelineRecommender, RNAseqDataProfiler

profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile()

recommender = MLPipelineRecommender()
recommendation = recommender.recommend(profile)

# NEW: After running pipeline, optimize thresholds
from raptor.threshold_optimizer import optimize_thresholds

de_results = pd.read_csv('pipeline_results.csv')
result = optimize_thresholds(de_results, goal='balanced')

print(f"Optimal |logFC|: {result.logfc_threshold:.3f}")
print(f"Significant genes: {result.n_significant}")
print(f"\nMethods:\n{result.methods_text}")
```

**Migration effort:** Add ATO import and calls where needed

---

## 📊 Results Compatibility

### Will My Results Change?

**Short answer:** No, unless you enable ATO.

| Feature | v2.1.0 Result | v2.1.1 Result | Notes |
|---------|---------------|---------------|-------|
| **Pipeline execution** | Identical | Identical | Same |
| **ML recommendations** | Identical | Identical | Same |
| **Ensemble analysis** | Identical | Identical | Unless `--use-ato` |
| **Dashboard** | Same | Same + ATO page | New feature |
| **Thresholds** | Fixed | **Data-driven with ATO** | Opt-in |

### Reproducing v2.1.0 Results

```bash
# Disable ATO (not needed - it's off by default)
raptor run --pipeline 3 --output results/

# Or explicitly in config
threshold_optimizer:
  enabled: false
```

---

## 🚨 Potential Issues & Solutions

### Issue 1: ATO Not Found

**Error:**
```python
ModuleNotFoundError: No module named 'raptor.threshold_optimizer'
```

**Solution:**
```bash
pip install --upgrade raptor-rnaseq
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('OK')"
```

### Issue 2: Column Not Found

**Error:**
```python
KeyError: 'log2FoldChange' not found
```

**Solution:**
```python
# Check your column names
print(df.columns.tolist())

# Use correct column names
result = optimize_thresholds(
    df,
    logfc_col='logFC',      # edgeR uses 'logFC'
    pvalue_col='PValue'     # edgeR uses 'PValue'
)
```

### Issue 3: Dashboard ATO Page Missing

**Solution:**
```bash
# Update RAPTOR
pip install --upgrade raptor-rnaseq

# Restart dashboard
raptor dashboard

# Look for "🎯 Threshold Optimizer" in sidebar
```

### Issue 4: π₀ Estimation Warning

**Warning:**
```
Warning: π₀ estimation failed, using default 0.9
```

**Solution:**
```python
# Usually means insufficient data or unusual p-value distribution
# The analysis will still work with default π₀

# Or try different method
result = optimize_thresholds(df, goal='balanced', pi0_method='histogram')
```

---

## 📚 Updated Documentation

### New Documentation Files

1. **[THRESHOLD_OPTIMIZER.md](THRESHOLD_OPTIMIZER.md)** - Comprehensive ATO guide

### Updated Documentation

- **INSTALLATION.md** - ATO verification steps
- **QUICK_START.md** - ATO quick example
- **DASHBOARD.md** - ATO page documentation
- **ENSEMBLE.md** - ATO integration
- **API.md** - threshold_optimizer module
- **FAQ.md** - Threshold selection Q&A
- **TROUBLESHOOTING.md** - ATO issues
- **PIPELINES.md** - ATO compatibility
- **DOCUMENTATION_INDEX.md** - ATO links

---

## 📋 Migration Checklist

### Before Migrating

- [ ] Read this migration guide
- [ ] Back up any in-progress analyses
- [ ] Note current threshold choices (for comparison)

### Installation

- [ ] Upgrade: `pip install --upgrade raptor-rnaseq`
- [ ] Verify version: `raptor --version` → 2.1.1
- [ ] Verify ATO: `python -c "from raptor.threshold_optimizer import optimize_thresholds; print('OK')"`

### Learning ATO

- [ ] Read [THRESHOLD_OPTIMIZER.md](THRESHOLD_OPTIMIZER.md)
- [ ] Try ATO on existing DE results
- [ ] Compare ATO thresholds to your usual choices
- [ ] Try dashboard ATO page

### Adoption

- [ ] Update scripts to include ATO (optional)
- [ ] Update config files with ATO settings (optional)
- [ ] Update documentation/SOPs
- [ ] Share methods text template with team

---

## 📈 Feature Comparison

### v2.1.0 vs v2.1.1

| Feature | v2.1.0 | v2.1.1 | Notes |
|---------|--------|--------|-------|
| **All v2.1.0 features** | ✅ | ✅ | Identical |
| **ML recommendations** | ✅ | ✅ | Same |
| **Dashboard** | ✅ | ✅ Enhanced | + ATO page |
| **Ensemble analysis** | ✅ | ✅ Enhanced | + ATO integration |
| **🎯 Threshold Optimizer** | ❌ | ✅ | **NEW** |
| **Publication methods text** | ❌ | ✅ | **NEW** |
| **π₀ estimation** | ❌ | ✅ | **NEW** |
| **Data-driven logFC** | ❌ | ✅ | **NEW** |

---

## 🎯 Quick Start with ATO

### Immediate Use (No Config Changes)

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load any DE results
df = pd.read_csv('deseq2_results.csv')

# Optimize
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)

# Use results
print(f"Optimal |logFC| threshold: {result.logfc_threshold:.3f}")
print(f"Significant genes: {result.n_significant}")

# Copy to paper
print(result.methods_text)

# Save
result.results_df.to_csv('optimized_results.csv')
```

### Dashboard Use

```bash
raptor dashboard
# Click "🎯 Threshold Optimizer"
# Upload, optimize, download!
```

---

## 📞 Getting Help

### Resources

1. **Documentation:** [THRESHOLD_OPTIMIZER.md](THRESHOLD_OPTIMIZER.md)
2. **FAQ:** [FAQ.md](FAQ.md) - Threshold selection section
3. **Troubleshooting:** [TROUBLESHOOTING.md](TROUBLESHOOTING.md)

### Community Support

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Email:** ayehbolouki1988@gmail.com

### Reporting Issues

```markdown
**Issue Title:** v2.1.1 Migration: [Brief Description]

**Information:**
- Previous version: v2.1.0
- New version: v2.1.1
- Feature: Threshold Optimizer

**Problem:**
[Describe what's not working]

**Error messages:**
```
[Paste error]
```
```

---

## 🗓️ Version Timeline

- **v2.0.0** (Jan 2025) - Initial release
- **v2.1.0** (Nov 2025) - ML, Dashboard, Ensemble
- **v2.1.1** (Dec 2025) - 🎯 Adaptive Threshold Optimizer ← **YOU ARE HERE**
- **v2.2.0** (Planned Q1 2026) - Single-cell support

---

## ✅ Conclusion

### Key Takeaways

✅ **v2.1.1 is backward compatible** - All v2.1.0 workflows work unchanged

✅ **ATO is opt-in** - Your existing analyses won't change

✅ **Easy to adopt** - One import, one function call

✅ **Publication-ready** - Auto-generated methods text

✅ **Dashboard integration** - No coding required

### Recommended Migration Timeline

- **Immediate:** If preparing manuscript
- **This week:** For active projects
- **When convenient:** For established workflows

---

**Welcome to RAPTOR v2.1.1 with Adaptive Threshold Optimizer!** 🎯🦖

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.1  
**Date:** December 2025  
**License:** MIT

---

*"Data-driven thresholds - because your data deserves better than arbitrary cutoffs."* 🎯
