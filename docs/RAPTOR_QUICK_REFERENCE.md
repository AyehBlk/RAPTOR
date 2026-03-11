# RAPTOR Quick Reference Cheat Sheet

**Print this for quick command reference!**

---

## ⚡ Installation

```bash
# PyPI (quick)
pip install raptor-rnaseq

# Conda (with tools)
conda env create -f environment-full.yml
conda activate raptor-full

# Check version
raptor --version
```

---

## 🔥 Most Common Commands

### **Complete Workflow**
```bash
# 1. Quality Check
raptor qc -c counts.csv -m metadata.csv -o qc/

# 2. Profile Data
raptor profile -c counts.csv -m metadata.csv -o profile/

# 3. Get Recommendation
raptor recommend -p profile/profile.json -o recommend/

# 4. Optimize Thresholds
raptor optimize -d de_results.csv -m fdr-control --fdr-target 0.05

# 5. Ensemble Analysis
raptor ensemble-compare --deseq2 d.csv --edger e.csv --limma l.csv -o ensemble/
```

---

## 📋 Module Quick Reference

### **Module 2: Quality Assessment**
```bash
# Basic
raptor qc -c counts.csv -m metadata.csv

# With plots
raptor qc -c counts.csv -m metadata.csv --plot --plot-output qc.pdf

# Python
from raptor import quick_quality_check
report = quick_quality_check('counts.csv', 'metadata.csv')
print(f"Outliers: {report.outliers}")
```

**Detects:** Outliers using 6 methods (MAD, Isolation Forest, LOF, PCA, Clustering, Statistical)

---

### **Module 3: Data Profiler**
```bash
# Basic
raptor profile -c counts.csv -m metadata.csv -g condition

# Python
from raptor import profile_data_quick
profile = profile_data_quick('counts.csv', 'metadata.csv', group_column='condition')
print(f"BCV: {profile.features['bcv']:.3f}")  # Most important!
```

**Extracts:** 32 features including BCV (key metric for pipeline selection)

**BCV Guide:**
- < 0.2: Low variation → limma-voom
- 0.2-0.4: Moderate (typical) → DESeq2, edgeR
- > 0.4: High variation → edgeR_robust

---

### **Module 4: Recommender**
```bash
# ML-based (recommended)
raptor recommend -p profile.json -m ml

# Both methods
raptor recommend -p profile.json -m both --verbose-explanation

# Python
from raptor import recommend_pipeline
rec = recommend_pipeline(profile_file='profile.json', method='ml')
print(f"Use: {rec.pipeline_name} (confidence: {rec.confidence:.2f})")
```

**Returns:** Pipeline recommendation with confidence score and reasoning

---

### **Module 7: DE Import**
```bash
# Import DESeq2
raptor import-de -i deseq2.csv -m deseq2 -o imported/

# Import edgeR
raptor import-de -i edger.csv -m edger

# Import limma
raptor import-de -i limma.csv -m limma

# Compare all
raptor compare-de deseq2.csv edger.csv limma.csv -o comparison/ --venn-diagram

# Python
from raptor import import_deseq2, compare_de_results
deseq2 = import_deseq2('deseq2.csv')
edger = import_edger('edger.csv')
comparison = compare_de_results(deseq2=deseq2, edger=edger)
```

**Standardizes:** Any DE format into common schema

---

### **Module 8: Parameter Optimization**
```bash
# Ground Truth (simulated data)
raptor optimize -d de.csv -m ground-truth -g truth.csv

# FDR Control (most common)
raptor optimize -d de.csv -m fdr-control --fdr-target 0.05

# Stability (cross-validation)
raptor optimize -d de.csv -m stability --counts counts.csv --metadata metadata.csv

# Reproducibility (independent cohorts)
raptor optimize -d de.csv -m reproducibility --counts c1.csv --cohort2 c2.csv

# Python
from raptor import optimize_with_fdr_control
result = optimize_with_fdr_control(de_result, fdr_target=0.05)
print(f"Optimal thresholds: FDR={result.optimal_threshold['padj']}, log2FC={result.optimal_threshold['lfc']}")
```

**Optimizes:** FDR and log2FC thresholds using 4 validated methods

---

### **Module 9: Ensemble Analysis**
```bash
# Single method (Fisher's)
raptor ensemble -m fisher --deseq2 d.csv --edger e.csv --limma l.csv -o ens/

# Compare all 5 methods
raptor ensemble-compare --deseq2 d.csv --edger e.csv --limma l.csv -o comparison/

# Python
from raptor import ensemble_fisher, ensemble_brown, ensemble_rra

# Fisher's Method
fisher = ensemble_fisher({
    'deseq2': deseq2_result,
    'edger': edger_result,
    'limma': limma_result
})

# Brown's Method (recommended - accounts for correlation)
brown = ensemble_brown({
    'deseq2': deseq2_result,
    'edger': edger_result,
    'limma': limma_result
})

print(f"Consensus genes: {len(brown.consensus_genes)}")
```

**Combines:** Multiple DE methods using 5 statistical approaches

---

## 🔧 Common Parameters

| Parameter | Flag | Description | Example |
|-----------|------|-------------|---------|
| **Counts** | `-c` | Count matrix | `-c counts.csv` |
| **Metadata** | `-m` | Sample metadata | `-m metadata.csv` |
| **Output** | `-o` | Output directory | `-o results/` |
| **Group** | `-g` | Group column | `-g condition` |
| **Threshold** | `--threshold` | Significance threshold | `--threshold 0.05` |
| **Method** | `-m` | Analysis method | `-m ml` |

---

## 📊 Pipeline Quick Reference

| Pipeline | Command | Memory | Time | Best For |
|----------|---------|--------|------|----------|
| **salmon** ⭐ | `raptor pipeline run -n salmon` | 8 GB | 10-20 min | **Recommended** |
| **kallisto** | `raptor pipeline run -n kallisto` | 4 GB | 5-10 min | Speed |
| **star_featurecounts** | `raptor pipeline run -n star_featurecounts` | 32 GB | 40-70 min | Publication genes |
| **star_rsem** | `raptor pipeline run -n star_rsem` | 32 GB | 60-120 min | Isoforms |

```bash
# List available
raptor pipeline list

# Run salmon
raptor pipeline run -n salmon -s samples.csv -g salmon_index/ -o results/ -t 8
```

---

## 🐍 Python API Quick Reference

### **Import Everything**
```python
from raptor import (
    # Module 2: QC
    quick_quality_check, DataQualityAssessor,
    
    # Module 3: Profiler
    profile_data_quick, RNAseqDataProfiler,
    
    # Module 4: Recommender
    recommend_pipeline, MLRecommender,
    
    # Module 7: DE Import
    import_deseq2, import_edger, import_limma, compare_de_results,
    
    # Module 8: Optimization
    optimize_with_ground_truth, optimize_with_fdr_control,
    optimize_with_stability, optimize_with_reproducibility,
    
    # Module 9: Ensemble
    ensemble_fisher, ensemble_brown, ensemble_rra,
    ensemble_voting, ensemble_weighted
)
```

### **Minimal Example**
```python
from raptor import quick_quality_check, profile_data_quick, ensemble_fisher

# QC
qc = quick_quality_check('counts.csv', 'metadata.csv')

# Profile
profile = profile_data_quick('counts.csv', 'metadata.csv', group_column='condition')

# Ensemble
result = ensemble_fisher({
    'deseq2': deseq2_result,
    'edger': edger_result
})
```

---

## 📁 File Formats

### **Count Matrix (counts.csv)**
```csv
gene_id,sample_1,sample_2,sample_3,sample_4
ENSG00001,1523,1872,2103,2987
ENSG00002,456,523,678,789
```
- Rows = genes
- Columns = samples
- Values = integer counts

### **Metadata (metadata.csv)**
```csv
sample_id,condition,batch,replicate
sample_1,control,1,1
sample_2,control,1,2
sample_3,treatment,1,1
sample_4,treatment,2,2
```
- Must have `sample_id` column
- Must match count matrix columns
- Include grouping variable (condition)

### **Sample Sheet for Pipelines (samples.csv)**
```csv
sample_id,fastq_1,fastq_2,condition
sample_1,s1_R1.fq.gz,s1_R2.fq.gz,control
sample_2,s2_R1.fq.gz,s2_R2.fq.gz,treatment
```

---

## ⚡ Quick Decision Flowchart

```
Do you have count matrix?
  ├─ YES → Start with Module 2 (QC)
  └─ NO  → Run Module 5 (Pipeline) or use quick-count

Is data quality good?
  ├─ YES → Continue to Module 3 (Profile)
  └─ NO  → Remove outliers, check batch effects

What pipeline should I use?
  → Module 4 (Recommender) - Use ML method!

Have DE results from DESeq2/edgeR/limma?
  → Module 7 (Import) → Module 9 (Ensemble)

Need to optimize thresholds?
  → Module 8 (FDR Control method most common)

Want robust results?
  → Module 9 (Ensemble - Brown's method recommended)
```

---

## 🆘 Common Errors & Fixes

### **Error: "Too few genes"**
```bash
# Check count matrix
head -n 5 counts.csv
wc -l counts.csv  # Should have >100 genes

# Filter low counts
python -c "
import pandas as pd
counts = pd.read_csv('counts.csv', index_col=0)
filtered = counts[counts.sum(axis=1) > 10]
filtered.to_csv('counts_filtered.csv')
"
```

### **Error: "Sample mismatch"**
```bash
# Check sample IDs
head -n 1 counts.csv
cut -d',' -f1 metadata.csv

# They must match exactly!
```

### **Error: "Module not found"**
```bash
# Reinstall
pip uninstall raptor-rnaseq
pip install raptor-rnaseq

# Verify
python -c "import raptor; print(raptor.__version__)"
```

---

## 📚 Get Help

```bash
# General help
raptor --help

# Module-specific help
raptor qc --help
raptor profile --help
raptor recommend --help
raptor optimize --help
raptor ensemble --help

# Check installation
raptor validate-installation

# Version
raptor --version
```

**Online:**
- 📖 Docs: [docs/](https://github.com/AyehBlk/RAPTOR/tree/main/docs)
- 🐛 Issues: [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)
- 📧 Email: ayehbolouki1988@gmail.com

---

## 🎯 Workflow Templates

### **Template 1: Quick QC Only**
```bash
raptor qc -c counts.csv -m metadata.csv
```

### **Template 2: QC + Profile**
```bash
raptor qc -c counts.csv -m metadata.csv
raptor profile -c counts.csv -m metadata.csv
```

### **Template 3: Complete Analysis**
```bash
raptor qc -c counts.csv -m metadata.csv -o qc/
raptor profile -c counts.csv -m metadata.csv -o profile/
raptor recommend -p profile/profile.json -o recommend/
raptor optimize -d de_results.csv -m fdr-control -o opt/
raptor ensemble-compare --deseq2 d.csv --edger e.csv --limma l.csv -o ens/
```

### **Template 4: Ensemble Only**
```bash
raptor import-de -i deseq2.csv -m deseq2 -o imported/
raptor import-de -i edger.csv -m edger -o imported/
raptor ensemble-compare \
    --deseq2 imported/deseq2_standardized.csv \
    --edger imported/edger_standardized.csv \
    -o ensemble_results/
```

---

## 💡 Pro Tips

1. **Always run QC first** - outliers will mess up everything downstream
2. **BCV is key** - Check profile.bcv to understand your data
3. **Use ML recommender** - It's trained on real data
4. **Brown's method for ensemble** - Accounts for correlation between DESeq2/edgeR
5. **FDR control for optimization** - Most practical method
6. **Keep intermediate files** - Useful for troubleshooting
7. **Document parameters** - Save commands in script for reproducibility

---

**Version:** 2.2.0  
**Updated:** March 2026

**Save this file for quick reference!**
