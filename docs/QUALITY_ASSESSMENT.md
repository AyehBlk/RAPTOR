# 🔍 RAPTOR v2.1.0 Quality Assessment Guide

**Comprehensive Data Quality Control and Assessment**

Automatically assess the quality of your RNA-seq data to ensure reliable results and identify potential issues before analysis.

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Quality Metrics](#quality-metrics)
4. [Quality Scoring](#quality-scoring)
5. [Contamination Detection](#contamination-detection)
6. [Batch Effect Analysis](#batch-effect-analysis)
7. [Sample Quality](#sample-quality)
8. [Interpretation Guide](#interpretation-guide)
9. [Troubleshooting Quality Issues](#troubleshooting-quality-issues)
10. [Best Practices](#best-practices)

---

## 🎯 Overview

### What is Quality Assessment?

Quality assessment evaluates your RNA-seq data to identify:

✅ **Low-quality samples** - Technical failures  
✅ **Contamination** - Non-target organisms  
✅ **Batch effects** - Systematic biases  
✅ **Outliers** - Unusual samples  
✅ **Library problems** - Prep issues  
✅ **Sequencing quality** - Read quality  

### Why Quality Assessment?

**Before assessment:**
```
Your Data → Analysis → Results
                       ↓
                    ??? Quality ???
```

**With quality assessment:**
```
Your Data → QC Check → Clean Data → Analysis → Reliable Results
            ↓
         Issues Found
            ↓
         Fix Problems
```

---

## ⚡ Quick Start

### Basic Quality Assessment

```bash
# Assess count matrix
raptor qc \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --output qc_report/
```

**Output:**
```
Quality Assessment Complete!

Overall Quality Score: 87/100 🟢 GOOD

✅ Library sizes: Good (CV = 12%)
✅ Gene detection: Excellent (18,234 genes)
✅ Zero inflation: Normal (42%)
⚠️  Batch effects: Minor (detected 1 batch)
✅ No contamination detected
⚠️  1 potential outlier: Sample_7

Recommendation: Remove Sample_7, include batch in model

Full report: qc_report/quality_report.html
```

### From FASTQ Files

```bash
# Assess raw FASTQ quality
raptor qc \
  --fastq data/fastq/ \
  --output qc_fastq_report/
```

### View Interactive Report

```bash
# Open report
xdg-open qc_report/quality_report.html

# Or via dashboard
raptor dashboard
# Click "Quality Assessment" tab
```

---

## 📊 Quality Metrics

### 1. Library Size Distribution

**What it measures:** Total counts per sample

```
Library Sizes (Million Reads)
    ↑
30M │     ▂▄█▆▃
25M │   ▃█████▆▂
20M │  ▅████████▄
15M │ ▇███████████▇
10M │██████████████
    └───────────────→ Samples

Mean: 24.5M
Median: 25.1M
CV: 12.3%
Range: 18.2M - 28.9M
```

**Interpretation:**
- **CV < 20%**: ✅ Good - Consistent libraries
- **CV 20-30%**: ⚠️ Acceptable - Some variation
- **CV > 30%**: ❌ Poor - Highly variable, investigate
- **Very low samples**: ❌ Failed libraries

**Quality Score Calculation:**
```python
if cv < 0.15:
    score = 100
elif cv < 0.20:
    score = 90
elif cv < 0.30:
    score = 75
else:
    score = max(50, 100 - cv * 100)
```

---

### 2. Gene Detection Rate

**What it measures:** Number of expressed genes per sample

```
Genes Detected per Sample
    ↑
20K │  ▂▄█▇▆▅▄▃
15K │▄██████████▇
10K │█████████████
 5K │█████████████
    └───────────────→ Samples

Average: 18,234 genes
Expected: 15,000-20,000 (human)
Samples below threshold: 1 (Sample_7)
```

**Interpretation:**
- **Human: 15,000-20,000**: ✅ Good
- **Mouse: 13,000-18,000**: ✅ Good
- **< 10,000**: ❌ Poor - Low complexity
- **Very high (>25,000)**: ⚠️ Check for contamination

**Quality Score:**
```python
if organism == "human":
    optimal = 17500
    tolerance = 2500
else:
    optimal = 15000
    tolerance = 2000

deviation = abs(detected - optimal)
score = max(0, 100 - (deviation / tolerance * 100))
```

---

### 3. Zero Inflation

**What it measures:** Percentage of zero counts

```
Zero Count Percentage per Sample

Sample_1:  [████████████░░░░] 42%
Sample_2:  [████████████░░░░] 41%
Sample_3:  [█████████████░░░] 45%
Sample_4:  [████████████░░░░] 43%
...

Average: 42.3%
Expected: 30-50% (typical RNA-seq)
Outliers: None
```

**Interpretation:**
- **30-50%**: ✅ Normal for bulk RNA-seq
- **< 30%**: ⚠️ Unusual (very deep sequencing?)
- **> 60%**: ❌ Poor - Low coverage or quality
- **> 80%**: ❌ Critical - Major issues

**Quality Score:**
```python
if 0.30 <= zero_rate <= 0.50:
    score = 100
elif 0.25 <= zero_rate <= 0.60:
    score = 85
elif 0.20 <= zero_rate <= 0.70:
    score = 70
else:
    score = max(40, 100 - abs(zero_rate - 0.40) * 200)
```

---

### 4. Biological Coefficient of Variation (BCV)

**What it measures:** Biological variability between replicates

```
BCV Analysis

Within Groups:
Control:     BCV = 0.38 ✅
Treatment:   BCV = 0.42 ✅

Overall:     BCV = 0.40 ✅

Category: Medium variation
Expected for: Typical experiments
```

**Interpretation:**
- **< 0.2**: Low variation (cell lines, controlled)
- **0.2-0.6**: Medium variation (most experiments)
- **> 0.6**: High variation (clinical, heterogeneous)
- **> 1.0**: Very high (investigate!)

**Not necessarily bad - depends on experiment!**

---

### 5. Sample Correlation

**What it measures:** How similar samples are

```
Sample Correlation Matrix

         S1   S2   S3   S4   S5   S6
    S1 [1.0  0.95 0.93 0.45 0.42 0.48]
    S2 [0.95 1.0  0.94 0.43 0.44 0.46]
    S3 [0.93 0.94 1.0  0.41 0.43 0.45]
    S4 [0.45 0.43 0.41 1.0  0.96 0.94]
    S5 [0.42 0.44 0.43 0.96 1.0  0.95]
    S6 [0.48 0.46 0.45 0.94 0.95 1.0 ]

Within-group: 0.94 ✅ High
Between-group: 0.44 ✅ Low (as expected)
```

**Interpretation:**
- **Within-group > 0.85**: ✅ Good replicates
- **Within-group 0.70-0.85**: ⚠️ Acceptable
- **Within-group < 0.70**: ❌ Poor - Check samples
- **Between-group < 0.60**: ✅ Groups distinct
- **Between-group > 0.85**: ❌ Groups not different!

---

### 6. PCA Clustering

**What it measures:** Whether samples cluster by condition

```
Principal Component Analysis

    PC2 (15%)
    ↑
    │    ●Control
    │  ●●●
    │●●●
    │      ▲▲▲Treatment
    │    ▲▲▲
    └─────────────→ PC1 (68%)

✅ Samples cluster by condition
❌ No batch effect visible
```

**Interpretation:**
- **Clustering by condition**: ✅ Good - Real biology
- **Clustering by batch**: ❌ Bad - Technical artifact
- **No clustering**: ❌ Weak signal or bad samples
- **Outliers**: ⚠️ Check individual samples

---

## 🎯 Quality Scoring

### Overall Quality Score (0-100)

```
Quality Score Calculation
═════════════════════════

Component Scores:
├─ Library Sizes:        92/100 (12% CV) ⭐
├─ Gene Detection:       90/100 (18,234 genes) ⭐
├─ Zero Inflation:       88/100 (42%) ⭐
├─ Sample Correlation:   95/100 (high within-group) ⭐
├─ BCV:                  85/100 (0.40 - medium) ⭐
├─ Batch Effects:        75/100 (minor detected) ⚠️
├─ Outliers:             80/100 (1 flagged) ⚠️
└─ Contamination:       100/100 (none detected) ⭐

Weights:
├─ Library Sizes:       15%
├─ Gene Detection:      20%
├─ Zero Inflation:      10%
├─ Correlation:         15%
├─ BCV:                 10%
├─ Batch Effects:       15%
├─ Outliers:            10%
└─ Contamination:        5%

OVERALL SCORE: 87/100 🟢 GOOD
```

### Score Interpretation

```
🟢 90-100: EXCELLENT
   ├─ High-quality data
   ├─ Proceed with confidence
   └─ All pipelines should work well

🟢 80-89: GOOD
   ├─ Good quality data
   ├─ Minor issues present
   └─ Proceed with awareness of issues

🟡 70-79: ACCEPTABLE
   ├─ Usable but not ideal
   ├─ Consider addressing issues
   └─ May affect sensitivity

🟡 60-69: MARGINAL
   ├─ Significant issues present
   ├─ Address problems before analysis
   └─ Results may be unreliable

🔴 <60: POOR
   ├─ Major quality problems
   ├─ Do not proceed without fixes
   └─ Consider re-sequencing
```

---

## 🦠 Contamination Detection

### What Gets Checked

```
Contamination Screen
═══════════════════

1. Species Verification
   Expected: Homo sapiens
   Detected: 98.5% H. sapiens ✅
            1.2% unmapped
            0.3% other

2. rRNA Contamination
   rRNA reads: 2.1% ✅ (acceptable)
   Expected: <5%

3. Adapter Contamination
   Adapter content: 0.8% ✅
   Expected: <1%

4. Vector Contamination
   Vector sequences: Not detected ✅

5. Cross-sample Contamination
   Sample barcode mixing: <0.1% ✅

Overall: No significant contamination ✅
```

### Detailed Contamination Analysis

```bash
raptor qc --contamination-detailed \
  --fastq data/fastq/ \
  --output contamination_report/
```

**Checks for:**
1. **Other organisms**
   - Bacteria
   - Fungi
   - Viruses
   - Other eukaryotes

2. **Technical contaminants**
   - Adapter sequences
   - Primer sequences
   - Vector sequences
   - PhiX (sequencing control)

3. **Sample mix-up**
   - Barcode bleeding
   - Index hopping
   - Cross-contamination

**Example output:**
```
Detailed Contamination Report
════════════════════════════

Sample: Sample_3
─────────────────

Organism Composition:
├─ Homo sapiens:     97.2% ✅
├─ Unmapped:          2.1%
├─ Bacteria:          0.5% ⚠️  (E. coli detected)
├─ Fungi:             0.1%
└─ Other:             0.1%

⚠️  WARNING: Bacterial contamination detected
   Source: Likely environmental (E. coli)
   Impact: Minimal (<1% of reads)
   Action: Optional - can filter or ignore

Technical Contaminants:
├─ Adapters:          0.3% ✅
├─ Primers:           0.0% ✅
├─ Vectors:           0.0% ✅
└─ PhiX:              0.0% ✅

Recommendation: Proceed, contamination is minimal
```

---

## 📉 Batch Effect Analysis

### Detecting Batch Effects

```bash
raptor qc --check-batch-effects \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --batch-column sequencing_run
```

**Output:**
```
Batch Effect Analysis
═══════════════════

Batch Variable: sequencing_run
Batches: Batch1, Batch2, Batch3

PCA Analysis:
PC1 (68%): Separates conditions ✅
PC2 (15%): Separates batches ⚠️

Batch Effect Strength: 15% (Moderate)

Visualization:
    PC2 (15%)
    ↑
    │  ●Batch1  ▲Batch2  ■Batch3
    │ ●●●  ▲▲▲  ■■■
    │●●●  ▲▲▲  ■■■
    │●●  ▲▲  ■■
    └─────────────→ PC1 (68%)

Statistical Tests:
├─ ANOVA on PC1: p < 0.001 (condition effect ✅)
├─ ANOVA on PC2: p = 0.023 (batch effect ⚠️)
└─ PVR (Proportion Variance Removed): 15%

Recommendation:
⚠️  Batch effects detected but manageable
   Action: Include batch in your model
   Model: ~ batch + condition
```

### Batch Effect Severity

```
Severity Levels:
═══════════════

PVR < 5%:  MINOR ✅
├─ Minimal impact
└─ Optional to correct

PVR 5-15%: MODERATE ⚠️
├─ Noticeable but manageable
├─ Include in model
└─ Model: ~ batch + condition

PVR 15-30%: STRONG ⚠️⚠️
├─ Significant impact
├─ Must correct
└─ Use ComBat or similar

PVR > 30%: SEVERE ❌
├─ Dominates signal
├─ May obscure biology
└─ Consider re-design
```

---

## 🎨 Sample Quality

### Individual Sample Assessment

```
Sample Quality Report
═══════════════════

Sample: Sample_1 ✅ PASS
├─ Library size: 25.2M (✅ good)
├─ Genes detected: 18,456 (✅ excellent)
├─ Zero inflation: 41% (✅ normal)
├─ Correlation with group: 0.94 (✅ high)
├─ Outlier score: 0.23 (✅ not outlier)
└─ Quality score: 93/100 🟢

Sample: Sample_7 ⚠️ FLAG
├─ Library size: 12.1M (❌ low - 50% of median)
├─ Genes detected: 11,234 (❌ low)
├─ Zero inflation: 68% (❌ high)
├─ Correlation with group: 0.62 (❌ low)
├─ Outlier score: 3.45 (❌ strong outlier)
└─ Quality score: 48/100 🔴

Recommendation:
⚠️  Sample_7 should be excluded from analysis
   Reason: Multiple quality flags
   Impact: May introduce noise and reduce power
```

### Outlier Detection Methods

**1. Statistical Outliers (Z-score)**
```python
outlier_score = (sample_value - mean) / std
if abs(outlier_score) > 3:
    flag_as_outlier()
```

**2. Distance-based (Mahalanobis)**
```python
distance = mahalanobis(sample, group_center, covariance)
if distance > threshold:
    flag_as_outlier()
```

**3. PCA-based**
```python
# Samples far from their group in PCA space
if pca_distance > 3 * std:
    flag_as_outlier()
```

**4. Correlation-based**
```python
# Low correlation with other samples in group
if median_correlation < 0.70:
    flag_as_outlier()
```

---

## 📖 Interpretation Guide

### Common Quality Patterns

#### Pattern 1: Perfect Data ✅

```
Quality Score: 95/100

Characteristics:
├─ Library sizes: CV < 15%
├─ Gene detection: Optimal range
├─ Zero inflation: 35-45%
├─ High within-group correlation (>0.90)
├─ Clear PCA separation
├─ No batch effects
└─ No outliers

Action: Proceed with any pipeline
```

#### Pattern 2: Minor Issues ⚠️

```
Quality Score: 82/100

Issues:
├─ Slight library size variation (CV = 22%)
└─ One potential outlier

Action: 
├─ Check outlier sample
├─ Consider robust methods
└─ Proceed with caution
```

#### Pattern 3: Batch Effects ⚠️

```
Quality Score: 76/100

Issues:
└─ Moderate batch effect (PVR = 18%)

Action:
├─ Include batch in model: ~ batch + condition
├─ Or use batch correction (ComBat)
└─ Verify results are robust
```

#### Pattern 4: Poor Quality ❌

```
Quality Score: 52/100

Issues:
├─ High library size variation (CV = 45%)
├─ Low gene detection (<12,000)
├─ Multiple outliers (3/12 samples)
└─ High zero inflation (>65%)

Action:
❌ Do NOT proceed
├─ Investigate failed samples
├─ Consider re-sequencing
└─ Consult with sequencing facility
```

---

## 🔧 Troubleshooting Quality Issues

### Issue 1: High Library Size Variation

**Diagnosis:**
```bash
raptor qc --diagnose library-sizes \
  --counts counts.csv
```

**Possible causes:**
1. Library prep variability
2. Sequencing depth differences
3. RNA quality differences
4. Technical failures

**Solutions:**
```bash
# Option 1: Normalize more aggressively
raptor analyze --normalize TMM

# Option 2: Exclude extreme samples
raptor analyze --filter-samples \
  --min-library-size 5000000

# Option 3: Use robust methods
raptor analyze --pipeline 6  # NOISeq handles variation well
```

---

### Issue 2: Low Gene Detection

**Diagnosis:**
```bash
raptor qc --diagnose gene-detection \
  --counts counts.csv
```

**Possible causes:**
1. Low sequencing depth
2. rRNA contamination
3. Poor RNA quality
4. Wrong reference genome

**Solutions:**
```bash
# Check sequencing depth
raptor qc --check-depth

# Check for rRNA
raptor qc --check-rrna

# If depth is issue:
# → Re-sequence or pool replicates

# If rRNA contamination:
# → Use rRNA-depleted library prep next time
# → Analyze as-is but note limitation
```

---

### Issue 3: Batch Effects

**Solutions:**
```bash
# Option 1: Include in model
raptor analyze --design "~ batch + condition"

# Option 2: ComBat correction
raptor analyze --batch-correction combat

# Option 3: Use limma (handles batches well)
raptor analyze --pipeline 5

# Verify correction worked
raptor qc --check-batch-correction \
  --original counts.csv \
  --corrected corrected_counts.csv
```

---

### Issue 4: Contamination

**For minor contamination (<5%):**
```bash
# Usually safe to proceed
raptor analyze --counts counts.csv

# Optional: filter contamination
raptor qc --remove-contamination \
  --output cleaned_counts.csv
```

**For major contamination (>10%):**
```bash
# Do NOT proceed
# Contact sequencing facility
# Consider re-sequencing
```

---

## 📚 Best Practices

### Before Sequencing

✅ **Plan for quality**
- Include QC samples
- Balance batches
- Randomize sample order
- Include spike-ins (optional)

✅ **Document everything**
- Sample metadata
- Batch information
- Library prep protocol
- Sequencing parameters

### After Sequencing

✅ **Always run QC first**
```bash
# Before any analysis
raptor qc --counts counts.csv --metadata metadata.csv
```

✅ **Review QC report carefully**
- Check overall score
- Read all warnings
- Understand each metric

✅ **Address issues before analysis**
- Remove bad samples
- Correct batch effects
- Filter contamination

### During Analysis

✅ **Include QC info in model**
```bash
# If batch effects detected
raptor analyze --design "~ batch + condition"

# If outliers present
raptor analyze --robust
```

✅ **Use robust methods for poor quality**
```bash
raptor analyze --pipeline 6  # NOISeq
# or
raptor analyze --pipeline 5  # limma with robust=TRUE
```

### Reporting

✅ **Always report QC metrics**
```
Materials and Methods:
"RNA-seq data quality was assessed using RAPTOR v2.1.0.
 Overall quality score was 87/100. One sample (Sample_7)
 was excluded due to low library size (12.1M vs median 25.2M)
 and low gene detection (11,234 genes). Minor batch effects
 (PVR = 8%) were accounted for by including batch as a
 covariate in the differential expression model."
```

---

## 🎉 Summary

Quality assessment provides:
- ✅ **Early problem detection** - Before wasting time
- ✅ **Objective metrics** - Not just visual inspection
- ✅ **Automated flagging** - Catches issues you might miss
- ✅ **Contamination detection** - Ensures data purity
- ✅ **Batch effect identification** - Prevent false discoveries
- ✅ **Sample QC** - Identify failures
- ✅ **Confidence in results** - Know your data is good

**Good QC is the foundation of good analysis!** 🔍✨

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**Version:** 2.1.0  
**License:** MIT

---

*"Garbage in, garbage out - QC saves the day!"* 🔍🎯
