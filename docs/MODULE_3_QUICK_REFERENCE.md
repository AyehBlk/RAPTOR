# RAPTOR Module 3: Data Profiling - Quick Reference

## ⚡ Quick Summary

Module 3 extracts **32 statistical features** from RNA-seq data to enable intelligent pipeline recommendation.

---

## 🎯 The 32 Features (8 Categories)

| Category | Features | Key Metrics |
|----------|----------|-------------|
| **Sample Characteristics** | 5 | n_samples, min_group_size, balance |
| **Library Size** | 4 | mean, CV, range |
| **Gene Detection** | 3 | detection_rate, reliable_rate |
| **Expression Distribution** | 4 | mean, variance, skewness |
| **Dispersion** ⭐ | 5 | **BCV, common_dispersion** |
| **Sparsity** | 3 | sparsity, zero_inflation |
| **Count Distribution** | 4 | low/medium/high proportions |
| **Mean-Variance** | 3 | slope, R², Poisson fit |

---

## 🔑 Most Critical Features

### 1. BCV (Biological Coefficient of Variation)

**THE MOST IMPORTANT FEATURE**

```
BCV = sqrt(common_dispersion)
```

| BCV | Interpretation | Best Pipeline |
|-----|----------------|---------------|
| < 0.2 | Low variation | limma-voom |
| 0.2-0.4 | Moderate (typical) | DESeq2, edgeR |
| > 0.4 | High variation | edgeR, edgeR_robust |

### 2. Min Group Size

| Size | Recommendation |
|------|----------------|
| < 3 | Insufficient - need more samples |
| 3-7 | DESeq2, edgeR (parametric) |
| 8-12 | All methods work |
| > 12 | limma-voom, Wilcoxon (non-parametric) |

### 3. Library Size CV

| CV | Action |
|----|--------|
| < 0.2 | Simple normalization |
| 0.2-0.5 | Standard (CPM, TMM) |
| > 0.5 | Careful normalization needed |

---

## 📋 Quick Start

```bash
# Basic profiling
raptor profile --counts counts.csv

# With metadata (recommended)
raptor profile --counts counts.csv \
               --metadata metadata.csv \
               --group-column condition

# Custom output
raptor profile -c counts.csv -m metadata.csv -o my_profile.json
```

---

## 🐍 Python API

```python
from raptor.profiler import profile_data_quick
import pandas as pd

# Load and profile
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')
profile = profile_data_quick(counts, metadata, group_column='treatment')

# Key metrics
print(f"BCV: {profile.bcv:.3f}")
print(f"Category: {profile.bcv_category}")
print(f"Min group size: {profile.min_group_size}")
print(f"Sparsity: {profile.sparsity:.1%}")

# Get recommendation features
rec_features = profile.get_recommendation_features()

# Save
with open('profile.json', 'w') as f:
    f.write(profile.to_json())
```

---

## 📊 Understanding Output

### Profile Summary

```
╔════════════════════════════════════════════════════════════╗
║              🦖 RAPTOR DATA PROFILE                        ║
╠════════════════════════════════════════════════════════════╣

  📊 DIMENSIONS
    Genes:        20,145
    Samples:           60
    Groups:             2
    Balance:         1.00

  📈 DISPERSION (Critical!)
    Common φ:      0.0289
    BCV:           0.170 (17.0% biological variation)
    Category:   moderate
    
  📊 COUNT DISTRIBUTION
    Zero:          45.2%
    Low (1-9):     22.3%
    Sparsity:   moderate

╚════════════════════════════════════════════════════════════╝
```

### What To Look For

✅ **Good Profile**:
- BCV: 0.1-0.4 (typical range)
- Min group size ≥ 3
- Sparsity: 30-60%
- Sample balance > 0.8

⚠️ **Potential Issues**:
- BCV > 0.6 (very high variation)
- Min group size < 3 (insufficient)
- Sparsity > 80% (unusual for bulk)
- Balance < 0.5 (unbalanced)

---

## 🎯 Pipeline Selection Guide

### Based on BCV + Sample Size

```
if BCV < 0.2:
    if n >= 8:
        → limma-voom or Wilcoxon
    else:
        → DESeq2 or limma-voom

elif BCV < 0.4:
    if n >= 8:
        → DESeq2, edgeR, or limma-voom
    else:
        → DESeq2 or edgeR

else:  # BCV >= 0.4
    if has_outliers:
        → edgeR_robust
    else:
        → edgeR
```

---

## 🔍 Interpreting Specific Features

### Sparsity

| Sparsity | Type | Typical For |
|----------|------|-------------|
| < 30% | Low | High-depth bulk RNA-seq |
| 30-60% | Moderate | Standard bulk RNA-seq |
| 60-80% | High | Low-depth or targeted |
| > 80% | Very high | Single-cell (use specialized tools) |

### Overdispersion Ratio

| Ratio | Interpretation |
|-------|----------------|
| ≈ 1 | Poisson (rare for RNA-seq) |
| 2-4 | Typical negative binomial |
| > 10 | Check data quality |

### Library Size Categories

| Category | Mean Library Size |
|----------|------------------|
| Shallow | < 5M reads |
| Standard | 5-30M reads |
| Deep | > 30M reads |

---

## 📈 Workflow Integration

```
Module 2: QC
    ↓
    Quality OK?
    ↓ Yes
Module 3: Profile (THIS MODULE)
    ↓
    Extract 32 features
    ↓
Module 4: Recommend
    ↓
    Select pipeline based on BCV + sample size
```

---

## ⚙️ Common Commands

```bash
# Profile with verbose output
raptor profile -c counts.csv -m metadata.csv --verbose

# Specify group column
raptor profile -c counts.csv -m metadata.csv -g treatment

# Different output location
raptor profile -c counts.csv -o results/my_profile.json

# Complete workflow
raptor qc -c counts.csv -m metadata.csv && \
raptor profile -c counts.csv -m metadata.csv && \
raptor recommend --profile results/profile/data_profile.json
```

---

## 🐛 Troubleshooting

### Issue: "BCV is NaN or extreme"
**Solution**: Check for outliers in Module 2 QC

### Issue: "Too few samples per group"
**Solution**: Need ≥2 samples per group (≥3 recommended)

### Issue: "Very high sparsity (>80%)"
**Check**: Is this single-cell data? Use specialized tools

### Issue: "Library size CV > 1.0"
**Check**: Failed libraries present? Remove them

---

## 💾 Output Files

```
results/profile/
├── data_profile.json          # Complete profile (all features)
├── profile_summary.txt        # Human-readable summary
└── feature_vector.csv         # 32-feature vector for ML
```

---

## 📚 Key Literature

- **Robinson et al. (2010)**: edgeR & BCV concept
- **Love et al. (2014)**: DESeq2 dispersion shrinkage
- **Li et al. (2022)**: Sample size thresholds (n=8)
- **Law et al. (2014)**: limma-voom mean-variance
- **Soneson & Delorenzi (2013)**: Method comparison

---

## ✅ Best Practices Checklist

- [ ] Always run QC (Module 2) before profiling
- [ ] Include metadata for group-based features
- [ ] Check BCV is in reasonable range (0.1-0.6)
- [ ] Verify sample sizes adequate (≥3 per group)
- [ ] Review profile summary before recommendation
- [ ] Save profile JSON for reproducibility
- [ ] Document BCV and key features

---

## 🎓 Key Concepts

### BCV (Biological Coefficient of Variation)
```
BCV = √(dispersion)

Represents biological variability as a percentage.
Example: BCV = 0.2 means 20% biological variation.
```

### Dispersion
```
Variance = Mean + Dispersion × Mean²

NB model parameter representing overdispersion.
```

### Why BCV Matters
Different pipelines handle biological variation differently:
- **Low BCV**: Empirical methods work (limma-voom)
- **High BCV**: Need robust variance estimation (edgeR)

---

## 🚀 Quick Decision Guide

**Want to know**: Which pipeline should I use?

**Quick answer**: 
1. Check `profile.bcv`
2. Check `profile.min_group_size`
3. Follow decision tree above

**Want to know**: Is my data good quality?

**Quick answer**:
1. Check QC score from Module 2
2. Check `profile.quality_score`
3. Look for outliers, batch effects

**Want to know**: Do I have enough samples?

**Quick answer**:
1. Check `profile.min_group_size`
2. Need ≥3 (minimum), ≥5 (good), ≥8 (excellent)

---

**Version**: 2.2.0  
**Module**: M3 (Stage 3)  
**Status**: Production Ready
