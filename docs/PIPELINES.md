# RAPTOR v2.1.1 Pipelines Reference

Deep dive into all 8 RNA-seq analysis pipelines with ML-powered selection guidance and **ATO compatibility**.

## 🆕 What's New in v2.1.1

- **🎯 ATO Compatibility** - All pipelines work with Adaptive Threshold Optimizer
- **Output Format Reference** - Column names for each pipeline's DE output
- **Threshold Optimization** - Data-driven thresholds for all pipelines

---

## Pipeline Overview

| ID | Name | Type | Speed | Accuracy | ATO Compatible |
|----|------|------|-------|----------|----------------|
| 1 | STAR-RSEM-DESeq2 | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚫⚫ | ✅ Yes |
| 2 | HISAT2-StringTie-Ballgown | Alignment | ⚫⚫⚫⚪⚪ | ⚫⚫⚫⚫⚪ | ✅ Yes |
| 3 | Salmon-edgeR | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚫⚪ | ✅ Yes ⭐ |
| 4 | Kallisto-DESeq2 | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚪⚪ | ✅ Yes |
| 5 | STAR-featureCounts-limma | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚫⚪ | ✅ Yes |
| 6 | Salmon-NOISeq | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚪⚪ | ⚠️ Special* |
| 7 | Bowtie2-RSEM-EBSeq | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚪⚪ | ⚠️ Special* |
| 8 | HISAT2-Cufflinks-Cuffdiff | Alignment | ⚫⚪⚪⚪⚪ | ⚫⚫⚪⚪⚪ | ✅ Yes |

*Special handling required - see pipeline section

---

## 🎯 ATO Integration (NEW in v2.1.1)

### All Pipelines Support ATO

After running any pipeline, use ATO for data-driven thresholds:

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load DE results from any pipeline
df = pd.read_csv('pipeline_results.csv')

# Optimize thresholds
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',  # Adjust per pipeline
    pvalue_col='pvalue',
    goal='balanced'
)

print(f"Optimal |logFC| > {result.logfc_threshold:.2f}")
print(f"Significant genes: {result.n_significant}")
```

### Pipeline Output Formats for ATO

| Pipeline | LogFC Column | P-value Column | Notes |
|----------|-------------|----------------|-------|
| 1 (DESeq2) | `log2FoldChange` | `pvalue` | Standard |
| 2 (Ballgown) | `fc` | `pval` | Check column names |
| 3 (edgeR) | `logFC` | `PValue` | Case-sensitive |
| 4 (DESeq2) | `log2FoldChange` | `pvalue` | Standard |
| 5 (limma) | `logFC` | `P.Value` | Note the dot |
| 6 (NOISeq) | `M` | `prob` | Probability, not p-value* |
| 7 (EBSeq) | `PostFC` | `PPDE` | Posterior probability* |
| 8 (Cuffdiff) | `log2(fold_change)` | `p_value` | Check exact name |

*NOISeq and EBSeq use different statistical frameworks - see pipeline sections

---

## Pipeline 1: STAR-RSEM-DESeq2

**Gold Standard - Highest Accuracy**

### Components
- **Alignment**: STAR (splice-aware)
- **Quantification**: RSEM (EM algorithm)
- **Statistics**: DESeq2 (negative binomial)

### Best For
- Publication-quality results
- Small sample sizes (n<6)
- Complex experimental designs
- When accuracy is paramount

### 🎯 ATO Integration (NEW)

```python
# DESeq2 output columns
from raptor.threshold_optimizer import optimize_thresholds

df = pd.read_csv('deseq2_results.csv')
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)
```

### Running with ATO

```bash
# 1. Run pipeline
raptor run --pipeline 1 --data fastq/ --output results/

# 2. Optimize thresholds
python -c "
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd
df = pd.read_csv('results/deseq2_results.csv')
result = optimize_thresholds(df, goal='balanced')
print(f'Use |logFC| > {result.logfc_threshold:.2f}')
result.results_df.to_csv('results/optimized_results.csv')
"
```

### ML Insights
- ML Success Rate: 92%
- Recommended for: n < 10 samples, high BCV
- ATO works excellently with DESeq2 output

---

## Pipeline 2: HISAT2-StringTie-Ballgown

**Transcript Assembly & Novel Discovery**

### Components
- **Alignment**: HISAT2
- **Assembly**: StringTie
- **Statistics**: Ballgown

### Best For
- Novel transcript discovery
- Isoform-level analysis
- Non-model organisms

### 🎯 ATO Integration (NEW)

```python
# Ballgown output columns may vary
df = pd.read_csv('ballgown_results.csv')

# Check column names first
print(df.columns.tolist())

# Common Ballgown columns
result = optimize_thresholds(
    df,
    logfc_col='fc',      # or 'log2fc'
    pvalue_col='pval',   # or 'qval'
    goal='discovery'
)
```

---

## Pipeline 3: Salmon-edgeR ⭐ ML TOP CHOICE

**Best Balance - Fast & Accurate**

### Components
- **Quantification**: Salmon (quasi-mapping)
- **Statistics**: edgeR (quasi-likelihood)

### Best For
- Most RNA-seq experiments
- Large datasets (>20 samples)
- Quick turnaround needed
- **Best ATO compatibility** (well-calibrated p-values)

### 🎯 ATO Integration (NEW) ⭐

```python
# edgeR output columns
df = pd.read_csv('edger_results.csv')
result = optimize_thresholds(
    df,
    logfc_col='logFC',      # Note: different from DESeq2
    pvalue_col='PValue',    # Note: capital P
    goal='balanced'
)

# edgeR p-values work excellently with ATO
print(f"π₀ estimate: {result.pi0:.3f}")
print(f"Optimal |logFC|: {result.logfc_threshold:.3f}")
```

### Why ATO Works Best with edgeR

- Well-calibrated p-values
- Good π₀ estimation
- Clean logFC distribution
- Most tested combination

### Running with ATO

```bash
raptor run --pipeline 3 --data fastq/ --use-ato --ato-goal balanced
```

---

## Pipeline 4: Kallisto-DESeq2

**Ultra-Fast - Large Studies**

### Components
- **Quantification**: Kallisto
- **Statistics**: DESeq2

### Best For
- Very large datasets (>50 samples)
- Exploratory analysis
- Minimal resources

### 🎯 ATO Integration (NEW)

```python
# Same as Pipeline 1 (DESeq2 output)
df = pd.read_csv('kallisto_deseq2_results.csv')
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='discovery'  # Good for exploratory
)
```

---

## Pipeline 5: STAR-featureCounts-limma

**Flexible Modeling - Complex Designs**

### Components
- **Alignment**: STAR
- **Counting**: featureCounts
- **Statistics**: limma-voom

### Best For
- Complex experimental designs
- Multi-factor analysis
- Batch correction needed

### 🎯 ATO Integration (NEW)

```python
# limma output columns
df = pd.read_csv('limma_results.csv')
result = optimize_thresholds(
    df,
    logfc_col='logFC',
    pvalue_col='P.Value',  # Note the dot!
    goal='balanced'
)
```

### Note on limma P-values

limma uses moderated t-statistics. ATO works well but:
- P-values may be more conservative
- Consider 'discovery' goal for more sensitivity

---

## Pipeline 6: Salmon-NOISeq

**Non-Parametric - Small Samples**

### Components
- **Quantification**: Salmon
- **Statistics**: NOISeq

### Best For
- Very small samples (n=2-3)
- No replicates

### 🎯 ATO Integration (SPECIAL)

⚠️ **NOISeq uses probability, not p-values!**

```python
# NOISeq output is different
df = pd.read_csv('noiseq_results.csv')

# NOISeq uses probability (prob), not p-value
# Higher prob = more likely DE (opposite of p-value!)

# Convert probability to pseudo p-value
df['pseudo_pvalue'] = 1 - df['prob']

result = optimize_thresholds(
    df,
    logfc_col='M',              # NOISeq uses 'M' for log-ratio
    pvalue_col='pseudo_pvalue',
    goal='balanced'
)
```

### Alternative: Use NOISeq's Native Threshold

NOISeq recommends prob > 0.8 or prob > 0.9. You can also:

```python
# Use NOISeq's native approach
significant = df[df['prob'] > 0.9]
```

---

## Pipeline 7: Bowtie2-RSEM-EBSeq

**Empirical Bayes - Isoform Analysis**

### Components
- **Alignment**: Bowtie2
- **Quantification**: RSEM
- **Statistics**: EBSeq

### Best For
- Isoform-level analysis
- Two-condition comparisons

### 🎯 ATO Integration (SPECIAL)

⚠️ **EBSeq uses posterior probability, not p-values!**

```python
# EBSeq output columns
df = pd.read_csv('ebseq_results.csv')

# EBSeq uses PPDE (Posterior Probability of being DE)
# Convert to pseudo p-value
df['pseudo_pvalue'] = 1 - df['PPDE']

result = optimize_thresholds(
    df,
    logfc_col='PostFC',          # Posterior fold change
    pvalue_col='pseudo_pvalue',
    goal='balanced'
)
```

### Alternative: Use EBSeq's Native Threshold

EBSeq recommends PPDE > 0.95 (FDR = 0.05):

```python
significant = df[df['PPDE'] > 0.95]
```

---

## Pipeline 8: HISAT2-Cufflinks-Cuffdiff

**Legacy Pipeline**

### Components
- **Alignment**: HISAT2
- **Assembly**: Cufflinks
- **Statistics**: Cuffdiff

### Best For
- Reproducing legacy analyses
- **Not recommended for new analyses**

### 🎯 ATO Integration (NEW)

```python
# Cuffdiff output columns
df = pd.read_csv('cuffdiff_results.csv')
result = optimize_thresholds(
    df,
    logfc_col='log2(fold_change)',  # Check exact name
    pvalue_col='p_value',
    goal='balanced'
)
```

### Note

ML rarely recommends this pipeline. Consider Pipeline 2 for transcript discovery or Pipeline 3 for standard DE.

---

## 🎯 ATO Best Practices by Pipeline

### Recommended Configurations

| Pipeline | ATO Goal | Notes |
|----------|----------|-------|
| 1 (STAR-RSEM-DESeq2) | balanced | Best for publications |
| 2 (HISAT2-StringTie) | discovery | More permissive for novel transcripts |
| 3 (Salmon-edgeR) ⭐ | balanced | Best ATO compatibility |
| 4 (Kallisto-DESeq2) | discovery | Good for large exploratory studies |
| 5 (STAR-limma) | balanced | Works well with complex designs |
| 6 (Salmon-NOISeq) | - | Use native prob threshold |
| 7 (Bowtie2-EBSeq) | - | Use native PPDE threshold |
| 8 (HISAT2-Cuffdiff) | balanced | Legacy support |

### Quick Reference Code

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Pipeline 1 & 4 (DESeq2)
result = optimize_thresholds(df, 'log2FoldChange', 'pvalue', goal='balanced')

# Pipeline 3 (edgeR)
result = optimize_thresholds(df, 'logFC', 'PValue', goal='balanced')

# Pipeline 5 (limma)
result = optimize_thresholds(df, 'logFC', 'P.Value', goal='balanced')

# Pipeline 2 & 8 (Ballgown/Cuffdiff) - check columns first
print(df.columns.tolist())
```

---

## ML-Powered Pipeline Selection with ATO

### Complete Workflow

```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# 1. Profile data
profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile()

# 2. Get ML recommendation
recommender = MLPipelineRecommender()
rec = recommender.recommend(profile)[0]
print(f"Recommended: {rec['pipeline_name']} ({rec['confidence']})")

# 3. Run recommended pipeline
# raptor run --pipeline {rec['pipeline_id']} ...

# 4. Optimize thresholds with ATO
de_results = pd.read_csv('de_results.csv')

# Get column names for this pipeline
if rec['pipeline_id'] in [1, 4]:
    logfc_col, pval_col = 'log2FoldChange', 'pvalue'
elif rec['pipeline_id'] == 3:
    logfc_col, pval_col = 'logFC', 'PValue'
elif rec['pipeline_id'] == 5:
    logfc_col, pval_col = 'logFC', 'P.Value'

result = optimize_thresholds(de_results, logfc_col, pval_col, goal='balanced')

print(f"\n🎯 Data-Driven Thresholds:")
print(f"   |logFC| > {result.logfc_threshold:.3f}")
print(f"   Significant: {result.n_significant} genes")
print(f"\nMethods text:\n{result.methods_text}")
```

---

## Ensemble Analysis with ATO

For critical analyses, combine pipelines with uniform ATO thresholds:

```yaml
# config.yaml
ensemble:
  enabled: true
  pipelines: [1, 3, 5]
  
threshold_optimizer:
  enabled: true
  goal: "balanced"
  ensemble_mode:
    uniform_thresholds: true  # Same thresholds for all pipelines
```

```python
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

# Combine results
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(pipeline_results, method='vote')

# Apply uniform ATO thresholds
result = optimize_thresholds(ensemble['combined'], goal='balanced')
```

---

## Decision Guide with ATO

### Quick Decision Tree

```
Start Here
    ↓
Need highest accuracy?
├─ YES → Pipeline 1 (STAR-RSEM-DESeq2) + ATO goal='balanced'
└─ NO → Continue
    ↓
Have >50 samples?
├─ YES → Pipeline 4 (Kallisto-DESeq2) + ATO goal='discovery'
└─ NO → Continue
    ↓
Need novel transcripts?
├─ YES → Pipeline 2 (HISAT2-StringTie) + ATO goal='discovery'
└─ NO → Continue
    ↓
Have <4 samples?
├─ YES → Pipeline 6 (NOISeq) + Use native prob threshold
└─ NO → Pipeline 3 (Salmon-edgeR) ⭐ + ATO goal='balanced'

ALL PIPELINES: Use ATO for data-driven thresholds!
```

---

## Configuration Reference

```yaml
# pipelines.yaml (v2.1.1)

pipelines:
  - id: 3
    name: "Salmon-edgeR"
    
    # NEW: ATO configuration
    threshold_optimizer_support:
      output_format: "edger"
      columns:
        logfc: "logFC"
        pvalue: "PValue"
        fdr: "FDR"
      ato_compatible: true
      recommended_goal: "balanced"
```

---

**RAPTOR v2.1.1 Pipelines**  
**Author**: Ayeh Bolouki  
**License**: MIT

*"Eight pipelines, one ML brain, optimized thresholds!"* 🎯🦖
