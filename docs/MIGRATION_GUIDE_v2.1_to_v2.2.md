# Migration Guide: v2.1.2 → v2.2.0

**RAPTOR Version Upgrade Guide**

This guide helps you migrate from RAPTOR v2.1.2 to v2.2.0, which includes significant improvements and some breaking changes.

---

## 📋 **Overview**

### **What Changed?**

- ✅ **Module 8 renamed and expanded** (Threshold Optimizer → Parameter Optimization)
- ✅ **NEW Module 9** (Ensemble Analysis)
- ✅ **Pipeline consolidation** (8 → 6 validated pipelines)
- ✅ **Module 1 removed** (was placeholder only)
- ✅ **Expanded documentation**
- ✅ **Broader Python support** (3.8-3.12)

### **Migration Difficulty: LOW** 🟢

Most changes involve:
- Updating import statements
- Using new function names
- Exploring new features

**Estimated time:** 15-30 minutes

---

## 🎯 **Quick Migration Checklist**

- [ ] Update RAPTOR: `pip install --upgrade raptor-rnaseq`
- [ ] Update imports for Module 8
- [ ] Review pipeline structure if using Module 5
- [ ] Update CLI commands if using old threshold optimizer
- [ ] Explore Module 9 (Ensemble Analysis) - NEW!
- [ ] Update conda environment if applicable
- [ ] Run tests to verify migration

---

## 🔄 **Breaking Changes & Fixes**

### **1. Module 8: Threshold Optimizer → Parameter Optimization**

#### **What Changed:**

Module 8 was renamed from "Threshold Optimizer" to "Parameter Optimization" and expanded from 1 method (ATO) to 4 methods.

#### **Migration Steps:**

**Old Code (v2.1.2):**
```python
from raptor.threshold_optimizer import ATO

# Create optimizer
optimizer = ATO()

# Run optimization
result = optimizer.optimize(
    de_results=de_result,
    target_fdr=0.05
)

# Get optimal thresholds
optimal_fdr = result['optimal_fdr']
optimal_lfc = result['optimal_lfc']
```

**New Code (v2.2.0):**
```python
from raptor import optimize_with_fdr_control

# Run optimization (more direct)
result = optimize_with_fdr_control(
    de_result=de_result,
    fdr_target=0.05,
    output_dir='optimization/fdr'
)

# Get optimal thresholds
optimal_fdr = result.optimal_threshold['padj']
optimal_lfc = result.optimal_threshold['lfc']
```

#### **CLI Migration:**

**Old (v2.1.2):**
```bash
# ATO command (if existed)
raptor threshold-optimize --results results.csv --target-fdr 0.05
```

**New (v2.2.0):**
```bash
raptor optimize \
    --de-result results.csv \
    --method fdr-control \
    --fdr-target 0.05 \
    --output optimization/
```

#### **Key Differences:**

| Aspect | v2.1.2 | v2.2.0 |
|--------|--------|--------|
| **Import** | `from raptor.threshold_optimizer import ATO` | `from raptor import optimize_with_fdr_control` |
| **Class** | `ATO()` object | Function-based API |
| **Methods** | 1 (ATO) | 4 (ground-truth, FDR, stability, reproducibility) |
| **Result** | Dictionary | `OptimizationResult` object |
| **CLI** | `threshold-optimize` | `optimize` |

---

### **2. Pipeline Structure Reorganization**

#### **What Changed:**

Consolidated from 8 pipelines to 6 validated pipelines with cleaner organization.

#### **Old Structure (v2.1.2):**

```
pipelines/
├── pipeline1_star_rsem_deseq2/
├── pipeline2_hisat2_stringtie_ballgown/  ← REMOVED
├── pipeline3_salmon_edger/
├── pipeline4_kallisto_sleuth/            ← REMOVED
├── pipeline5_star_htseq_limma/
├── pipeline6_star_featurecounts_noiseq/  ← REMOVED
├── pipeline7_bowtie2_rsem_ebseq/         ← REMOVED
└── pipeline8_hisat2_cufflinks_cuffdiff/  ← REMOVED
```

#### **New Structure (v2.2.0):**

```
raptor/pipelines/
├── salmon/                     # Recommended
├── kallisto/                   # Fastest
├── star_featurecounts/         # Gene-level with BAM
├── star_rsem/                  # Gene + isoform with BAM
├── star_salmon/                # BAM + bootstraps
└── hisat2_featurecounts/       # Low memory
```

#### **Migration Steps:**

**If you used removed pipelines:**

| Old Pipeline | Migration Path |
|--------------|----------------|
| **pipeline2** (HISAT2+StringTie+Ballgown) | Use `hisat2_featurecounts` + Module 7 import |
| **pipeline4** (Kallisto+Sleuth) | Use `kallisto` + Module 7 import |
| **pipeline6** (STAR+featureCounts+NOISeq) | Use `star_featurecounts` + Module 7 import |
| **pipeline7** (Bowtie2+RSEM+EBSeq) | Use `star_rsem` + Module 7 import |
| **pipeline8** (HISAT2+Cufflinks+Cuffdiff) | Use `hisat2_featurecounts` + Module 7 import |

**Example Migration:**

**Old (v2.1.2) - Pipeline 4:**
```bash
# Used Kallisto + Sleuth directly
python pipelines/pipeline4_kallisto_sleuth/run.py
```

**New (v2.2.0) - Kallisto + Module 7:**
```bash
# Step 1: Run Kallisto quantification
raptor pipeline run --name kallisto --samples samples.csv

# Step 2: Run Sleuth in R separately

# Step 3: Import Sleuth results
raptor import-de --input sleuth_results.csv --method limma
```

---

### **3. Module 1 Removed**

#### **What Changed:**

Module 1 (Quick Count) was removed as it was a placeholder and never implemented.

#### **Migration Steps:**

**If you referenced Module 1:**

Module 1 didn't have any functionality, so no code changes needed. Just note that:
- Module numbering now starts at 2
- All subsequent modules keep their numbers (2-4, 7-9)

**Old documentation references:**
```
Module 1: Quick Count      ← REMOVED
Module 2: Quality Assessment
Module 3: Data Profiler
...
```

**New documentation:**
```
Module 2: Quality Assessment
Module 3: Data Profiler
Module 4: Recommender
...
```

---

### **4. Python Version Support**

#### **What Changed:**

Expanded Python support from just 3.10 to 3.8-3.12.

#### **Migration Steps:**

**If you're using Python 3.10:** No changes needed ✅

**If you want to use Python 3.8, 3.9, 3.11, or 3.12:**

```bash
# Create new environment with desired Python version
conda create -n raptor-new python=3.11
conda activate raptor-new
pip install raptor-rnaseq==2.2.0
```

---

## ✨ **New Features to Explore**

### **1. Module 9: Ensemble Analysis (NEW!)**

Combine DE results from multiple methods for robust consensus.

**Quick Start:**

```python
from raptor import ensemble_brown, import_deseq2, import_edger, import_limma

# Import results from different tools
deseq2 = import_deseq2('deseq2_results.csv')
edger = import_edger('edger_results.csv')
limma = import_limma('limma_results.csv')

# Combine using Brown's method (correlation-aware)
result = ensemble_brown(
    de_results={
        'deseq2': deseq2,
        'edger': edger,
        'limma': limma
    },
    significance_threshold=0.05
)

# View consensus genes
print(f"Consensus DE genes: {len(result.consensus_genes)}")
print(f"Direction consistency: {result.direction_consistency.mean():.1%}")

# Save results
result.save('ensemble_results/')
```

**5 Methods Available:**
1. `ensemble_fisher()` - Fisher's method
2. `ensemble_brown()` - Brown's method (recommended)
3. `ensemble_rra()` - Robust Rank Aggregation
4. `ensemble_voting()` - Voting consensus
5. `ensemble_weighted()` - Weighted ensemble

**CLI:**
```bash
# Single method
raptor ensemble \
    --methods brown \
    --deseq2 deseq2.csv \
    --edger edger.csv \
    --limma limma.csv \
    --output ensemble/

# Compare all methods
raptor ensemble-compare \
    --deseq2 deseq2.csv \
    --edger edger.csv \
    --limma limma.csv \
    --output ensemble_comparison/
```

---

### **2. Module 8: Additional Optimization Methods**

Beyond FDR control, you now have 3 more optimization approaches.

**Method 1: Ground Truth (for simulated data)**
```python
from raptor import optimize_with_ground_truth

result = optimize_with_ground_truth(
    de_result=de_result,
    ground_truth=true_positive_genes,  # Known DE genes
    output_dir='optimization/ground_truth'
)
```

**Method 2: Stability (cross-validation)**
```python
from raptor import optimize_with_stability

result = optimize_with_stability(
    counts=counts,
    metadata=metadata,
    n_folds=5,  # 5-fold CV
    output_dir='optimization/stability'
)
```

**Method 3: Reproducibility (independent cohorts)**
```python
from raptor import optimize_with_reproducibility

result = optimize_with_reproducibility(
    counts=counts_cohort1,
    metadata=metadata1,
    cohort2=counts_cohort2,
    output_dir='optimization/reproducibility'
)
```

---

### **3. Enhanced Module 3: 32-Feature Profiler**

Module 3 now extracts 32 features for ML-based recommendations.

**Key New Metric: BCV (Biological Coefficient of Variation)**

```python
from raptor import profile_data_quick

profile = profile_data_quick(
    counts='counts.csv',
    metadata='metadata.csv',
    group_column='condition'
)

# Access key metrics
print(f"BCV: {profile.bcv:.3f}")
print(f"BCV category: {profile.bcv_category}")
print(f"Sample size: {profile.n_samples}")
print(f"Sparsity: {profile.sparsity:.1%}")

# Interpretation
if profile.bcv < 0.2:
    print("Low variation → Use limma-voom")
elif profile.bcv < 0.4:
    print("Moderate variation → Use DESeq2 or edgeR")
else:
    print("High variation → Use edgeR_robust")
```

**All 32 features available:**
```python
# Get all features
features = profile.features
print(f"Total features: {len(features)}")

# Access feature vector for ML
feature_names = profile.feature_names
feature_values = profile.feature_vector
```

---

### **4. Comprehensive Documentation**

New documentation files available:

```
docs/
├── RAPTOR_USER_GUIDE.md          # Complete tutorial
├── RAPTOR_API_DOCUMENTATION.md   # Python API reference
├── RAPTOR_QUICK_REFERENCE.md     # Command cheat sheet
├── CONDA_ENVIRONMENTS.md         # Conda setup
├── MODULE_8_Parameter_Optimization.md  # NEW
├── MODULE_9_Ensemble_Analysis.md       # NEW
└── ...
```

**Read the USER_GUIDE for comprehensive tutorials!**

---

### **5. Conda Environment Support**

Two new conda environment files:

**Core environment (fast):**
```bash
conda env create -f environment.yml
conda activate raptor
# ~500 MB, 5-10 minutes
```

**Full environment (complete):**
```bash
conda env create -f environment-full.yml
conda activate raptor-full
# ~5-8 GB, 30-60 minutes
# Includes STAR, Salmon, Kallisto, R packages
```

---

## 🔍 **Testing Your Migration**

### **1. Verify Installation**

```bash
# Check version
raptor --version
# Should show: RAPTOR version 2.2.0

python -c "import raptor; print(raptor.__version__)"
# Should show: 2.2.0
```

### **2. Test Module 8 (Parameter Optimization)**

```python
import pandas as pd
from raptor import optimize_with_fdr_control

# Load test DE results
de_result = pd.read_csv('test_de_results.csv')

# Run optimization
result = optimize_with_fdr_control(
    de_result=de_result,
    fdr_target=0.05
)

# Verify results
assert result.optimal_threshold is not None
print("✅ Module 8 working!")
```

### **3. Test Module 9 (Ensemble Analysis)**

```python
from raptor import ensemble_brown, import_deseq2, import_edger

# Load test results
deseq2 = import_deseq2('test_deseq2.csv')
edger = import_edger('test_edger.csv')

# Run ensemble
result = ensemble_brown(
    de_results={'deseq2': deseq2, 'edger': edger}
)

# Verify results
assert len(result.consensus_genes) > 0
print("✅ Module 9 working!")
```

### **4. Run Full Test Suite**

```bash
# If you have tests/
pytest tests/

# Or run module-specific tests
python -m pytest tests/test_parameter_optimization.py
python -m pytest tests/test_ensemble.py
```

---

## 📚 **Additional Resources**

### **Documentation**

- [RAPTOR_USER_GUIDE.md](RAPTOR_USER_GUIDE.md) - Complete tutorial
- [RAPTOR_QUICK_REFERENCE.md](RAPTOR_QUICK_REFERENCE.md) - Command cheat sheet
- [RAPTOR_API_DOCUMENTATION.md](RAPTOR_API_DOCUMENTATION.md) - API reference
- [CHANGELOG.md](../CHANGELOG.md) - Detailed change log

### **Module Documentation**

- [MODULE_8_Parameter_Optimization.md](MODULE_8_Parameter_Optimization.md)
- [MODULE_9_Ensemble_Analysis.md](MODULE_9_Ensemble_Analysis.md)
- [MODULE_3_Data_Profiling.md](MODULE_3_Data_Profiling.md)

### **Examples**

Check the `examples/` directory for updated example scripts:
- `08_Parameter_Optimization.py`
- `09_Ensemble_Analysis.py`

---

## 🆘 **Troubleshooting**

### **Issue: Import errors after upgrade**

```python
ImportError: cannot import name 'ATO' from 'raptor.threshold_optimizer'
```

**Solution:**
```python
# Old
from raptor.threshold_optimizer import ATO

# New
from raptor import optimize_with_fdr_control
```

---

### **Issue: Old pipeline structure not found**

```bash
Error: pipeline2_hisat2_stringtie_ballgown not found
```

**Solution:** Use new pipeline structure
```bash
# Old
python pipelines/pipeline2_hisat2_stringtie_ballgown/run.py

# New
raptor pipeline run --name hisat2_featurecounts --samples samples.csv
```

---

### **Issue: Version conflict**

```bash
ERROR: raptor-rnaseq 2.1.2 has requirement numpy>=1.21.0, but you have numpy 1.19.0
```

**Solution:** Update dependencies
```bash
pip install --upgrade raptor-rnaseq
# This will update all dependencies to compatible versions
```

---

### **Issue: Conda environment conflicts**

**Solution:** Create fresh environment
```bash
# Remove old environment
conda env remove -n raptor

# Create new v2.2.0 environment
conda env create -f environment.yml
conda activate raptor
```

---

## ✅ **Migration Complete!**

Once you've:
- ✅ Updated RAPTOR to v2.2.0
- ✅ Updated your code (imports, function calls)
- ✅ Tested core functionality
- ✅ Explored new features (Module 9!)

**You're ready to use RAPTOR v2.2.0!** 🎉

---

## 🚀 **Next Steps**

1. **Explore Module 9** - Try ensemble analysis with your DE results
2. **Use new optimization methods** - Test ground-truth or stability optimization
3. **Read the USER_GUIDE** - Comprehensive tutorial for all features
4. **Update your workflows** - Incorporate new features
5. **Share feedback** - Open issues or discussions on GitHub

---

## 💬 **Need Help?**

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Discussions:** https://github.com/AyehBlk/RAPTOR/discussions
- **Email:** ayehbolouki1988@gmail.com

---

**Migration Guide Version:** 1.0  
**Last Updated:** March 2026  
**RAPTOR Version:** 2.2.0
