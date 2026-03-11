# 🦖 RAPTOR Complete Workflow - Modules 6-10 Overview

**Version:** 2.2.0  
**Date:** January 19, 2026

---

## 🎯 RAPTOR ARCHITECTURE OVERVIEW

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         RAPTOR PIPELINE ARCHITECTURE                     │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  STAGE 1: QUANTIFICATION                                                │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M1: Quick Count (Kallisto/Salmon)                                │  │
│  │     Input: FASTQ files                                           │  │
│  │     Output: quick_gene_counts.csv                                │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  STAGE 2: QUALITY & PROFILING                                           │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M2: Quality Assessment                                           │  │
│  │     • Outlier detection                                          │  │
│  │     • QC metrics                                                 │  │
│  │     Output: outliers.txt, qc_report.html                        │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M3: Data Profiling (32 features)                                │  │
│  │     • Sample size, depth, variance                               │  │
│  │     • Batch effects, library complexity                          │  │
│  │     Output: profile.json                                         │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M4: Pipeline Recommender (ML-based)                             │  │
│  │     • Best pipeline for your data                                │  │
│  │     • Optimal DE parameters                                      │  │
│  │     Output: recommendation.yaml                                  │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M5: Production Pipeline (7 options)                             │  │
│  │     • HISAT2/STAR/Kallisto/Salmon + featureCounts/RSEM         │  │
│  │     Output: gene_counts.csv (final counts)                       │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                                                          │
├══════════════════════════════════════════════════════════════════════════┤
│                                                                          │
│  STAGE 3: DIFFERENTIAL EXPRESSION (EXTERNAL R MODULE) 🔴                │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M6: DE Analysis (Run OUTSIDE RAPTOR)                            │  │
│  │     ┌────────────────────────────────────────────────────────┐  │  │
│  │     │ Option 1: DESeq2 (recommended)                         │  │  │
│  │     │ Option 2: edgeR (small samples)                        │  │  │
│  │     │ Option 3: limma-voom (complex designs)                 │  │  │
│  │     └────────────────────────────────────────────────────────┘  │  │
│  │                                                                  │  │
│  │     Command: Rscript run_deseq2.R \                             │  │
│  │                --counts gene_counts.csv \                        │  │
│  │                --metadata metadata.csv \                         │  │
│  │                --config recommendation.yaml                      │  │
│  │                                                                  │  │
│  │     Output: de_results.csv (standardized format)                │  │
│  │             de_summary.json                                      │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M7: Import DE Results (BACK TO RAPTOR)                          │  │
│  │     • Validate format                                            │  │
│  │     • Load to database                                           │  │
│  │     Command: raptor import-de --de-file de_results.csv          │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                                                          │
├══════════════════════════════════════════════════════════════════════════┤
│                                                                          │
│  STAGE 4: OPTIMIZATION & ENSEMBLE                                       │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M8: Parameter Optimization                                       │  │
│  │     • Optimize FDR threshold                                     │  │
│  │     • Optimize LFC cutoff                                        │  │
│  │     • Maximize reproducibility                                   │  │
│  │     Output: optimized_params.yaml                                │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M9: Ensemble Analysis                                            │  │
│  │     • Combine multiple DE methods                                │  │
│  │     • Consensus gene lists                                       │  │
│  │     • Weighted voting                                            │  │
│  │     Output: ensemble_results.csv                                 │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                ↓                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ M10: Biomarker Discovery                                         │  │
│  │     • Rank genes by biomarker potential                          │  │
│  │     • Multi-criteria scoring                                     │  │
│  │     • Pathway enrichment                                         │  │
│  │     Output: biomarkers.csv, ranked_genes.csv                     │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## 📋 MODULE DETAILS

### Module 6: Differential Expression Analysis (External) 🔴

**Type:** External R Module  
**Language:** R/Bioconductor  
**Status:** v2.2.0 Complete

**Purpose:**
Statistical differential expression analysis using industry-standard R packages.

**Why External?**
- ✅ Uses validated R/Bioconductor tools (DESeq2, edgeR, limma)
- ✅ Full control over statistical parameters
- ✅ Access to all package features
- ✅ Familiar workflow for bioinformaticians

**Input:**
- Count matrix (from M5)
- Sample metadata
- RAPTOR recommendations (from M4)

**Analysis Options:**
1. **DESeq2** - General use, recommended for most experiments
2. **edgeR** - Small sample sizes, fast analysis
3. **limma-voom** - Complex designs, maximum flexibility

**Output (Standardized):**
```
results/de_analysis/
├── de_results.csv        # Full results (all genes)
├── de_significant.csv    # Significant genes only
├── de_summary.json       # Analysis metadata
└── de_plots.pdf          # QC plots (optional)
```

**How to Run:**
```bash
# Step 1: Install R packages (first time only)
Rscript install_packages.R

# Step 2: Run analysis
Rscript run_deseq2.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_analysis \
    --config recommendation.yaml \
    --plots

# Step 3: Results ready for import!
```

**Key Features:**
- ✅ Batch correction support
- ✅ LFC shrinkage (DESeq2)
- ✅ Multiple normalization methods
- ✅ Comprehensive QC plots
- ✅ RAPTOR parameter integration
- ✅ Standardized output format

---

### Module 7: Import DE Results

**Type:** Python Pipeline  
**Status:** Future Implementation

**Purpose:**
Import external DE results back into RAPTOR for downstream analysis.

**Input:**
- de_results.csv (from M6)
- de_summary.json (metadata)

**Process:**
1. **Validate Format**
   - Check required columns
   - Verify data types
   - Check value ranges

2. **Parse Results**
   - Load gene IDs
   - Extract statistics
   - Parse metadata

3. **Store in Database**
   - RAPTOR internal format
   - Ready for M8-M10

**Output:**
- DE results in RAPTOR database
- Import confirmation
- Validation report

**Command:**
```bash
raptor import-de \
    --de-file results/de_analysis/de_results.csv \
    --summary results/de_analysis/de_summary.json \
    --method deseq2
```

**Validation Checks:**
- ✅ Required columns present
- ✅ Gene IDs match reference
- ✅ No missing critical values
- ✅ FDR values valid (0-1)
- ✅ Metadata consistent

---

### Module 8: Parameter Optimization

**Type:** Python Pipeline  
**Status:** Future Implementation

**Purpose:**
Optimize statistical thresholds (FDR, LFC) using data-driven approaches.

**Input:**
- DE results (from M7)
- Count matrix
- Sample metadata

**Optimization Strategies:**

1. **Reproducibility-Based**
   - Maximize overlap between resampled datasets
   - Find stable threshold

2. **Sensitivity-Specificity Balance**
   - ROC curve analysis (if truth available)
   - Maximize F1 score

3. **Biological Coherence**
   - Pathway enrichment consistency
   - Gene ontology coherence

4. **Effect Size Distribution**
   - Model null and alternative distributions
   - Find optimal separation

**Output:**
```yaml
optimized_parameters:
  fdr_threshold: 0.03
  lfc_threshold: 0.58
  min_count: 12
  
optimization_metrics:
  reproducibility: 0.92
  estimated_fdr: 0.028
  power: 0.87
  
recommendations:
  - FDR 0.05 too lenient (estimated FDR: 0.08)
  - Recommended: FDR 0.03 (estimated FDR: 0.028)
  - LFC threshold improves specificity
```

**Command:**
```bash
raptor optimize-params \
    --de-results results/de_analysis/de_results.csv \
    --optimize fdr,lfc \
    --criterion reproducibility \
    --resamples 100
```

**Benefits:**
- ✅ Data-driven thresholds
- ✅ Improved reproducibility
- ✅ Balanced sensitivity/specificity
- ✅ Reduced false positives

---

### Module 9: Ensemble Analysis

**Type:** Python Pipeline  
**Status:** Future Implementation

**Purpose:**
Combine multiple DE methods to create robust, consensus gene lists.

**Input:**
- Multiple DE results (DESeq2, edgeR, limma)
- Weighting scheme

**Ensemble Strategies:**

1. **Consensus (Intersection)**
   - Genes significant in ALL methods
   - Most conservative
   - Highest confidence

2. **Union**
   - Genes significant in ANY method
   - Most liberal
   - Highest sensitivity

3. **Majority Vote**
   - Genes significant in ≥N methods
   - Balanced approach

4. **Weighted Voting**
   - Weight by method confidence
   - Combine p-values/rankings
   - Advanced meta-analysis

5. **Rank Aggregation**
   - Combine rankings from each method
   - Robust order statistics

**Output:**
```
results/ensemble/
├── ensemble_results.csv      # Combined results
├── method_comparison.pdf     # Venn diagrams
├── concordance_matrix.csv    # Method agreement
└── consensus_genes.txt       # High-confidence list
```

**Command:**
```bash
raptor ensemble-analysis \
    --deseq2 results/de_deseq2/de_results.csv \
    --edger results/de_edger/de_results.csv \
    --limma results/de_limma/de_results.csv \
    --strategy consensus \
    --weights 0.4,0.3,0.3
```

**Benefits:**
- ✅ Reduced false positives
- ✅ Increased confidence
- ✅ Method-independent results
- ✅ Robust to method-specific biases

**Analysis Outputs:**
- Consensus gene list
- Venn diagrams (2-way, 3-way)
- Correlation of log fold changes
- Method concordance metrics

---

### Module 10: Biomarker Discovery

**Type:** Python Pipeline  
**Status:** Future Implementation

**Purpose:**
Identify and rank potential biomarkers from differential expression results.

**Input:**
- DE results (from M7 or M9)
- Expression matrix
- Sample metadata
- Clinical outcomes (optional)

**Biomarker Scoring Criteria:**

1. **Statistical Significance**
   - Low adjusted p-value
   - High test statistic

2. **Effect Size**
   - Large fold change
   - High discriminative power

3. **Expression Level**
   - Adequate expression
   - Detectable in target tissue

4. **Consistency**
   - Stable across batches
   - Reproducible in resamples

5. **Biological Relevance**
   - Known disease association
   - Pathway membership
   - Literature evidence

6. **Clinical Utility** (if data available)
   - Correlation with outcomes
   - Predictive power
   - Diagnostic accuracy

**Multi-Criteria Ranking:**
```python
biomarker_score = (
    w1 * significance_score +
    w2 * effect_size_score +
    w3 * expression_score +
    w4 * consistency_score +
    w5 * biological_relevance_score +
    w6 * clinical_utility_score
)
```

**Output:**
```
results/biomarkers/
├── ranked_genes.csv           # All genes ranked
├── top_biomarkers.csv         # Top N candidates
├── biomarker_panels.txt       # Multi-gene panels
├── enrichment_analysis.csv    # Pathway enrichment
├── expression_heatmap.pdf     # Top genes
└── biomarker_report.html      # Interactive report
```

**Command:**
```bash
raptor discover-biomarkers \
    --de-results results/ensemble/ensemble_results.csv \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --rank-by significance,effect-size,consistency \
    --top-n 50 \
    --generate-panels
```

**Biomarker Panel Generation:**
- Single-gene biomarkers (top N)
- Multi-gene signatures (2-10 genes)
- Pathway-based panels
- Machine learning classifiers

**Output Analysis:**
- ROC curves (if outcomes available)
- Expression boxplots
- Correlation matrices
- Pathway enrichment
- Literature mining

**Benefits:**
- ✅ Prioritized candidate list
- ✅ Multi-criteria evaluation
- ✅ Ready for validation
- ✅ Publication-ready figures

---

## 🔄 COMPLETE WORKFLOW EXAMPLE

### Step-by-Step Execution

```bash
# ============================================================================
# STAGE 1-2: Quantification, QC, Profiling (RAPTOR Python)
# ============================================================================

# Module 1: Quick quantification
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/ \
    -o results/quick_counts

# Module 2: Quality assessment
raptor qc \
    --counts results/quick_counts/quick_gene_counts.csv \
    --metadata data/metadata.csv

# Module 3: Data profiling
raptor profile \
    --counts results/quick_counts/quick_gene_counts.csv

# Module 4: Get recommendation
raptor recommend

# Module 5: Production pipeline (recommended method)
raptor run-pipeline \
    --pipeline salmon \
    --sample-sheet samples.csv \
    --index salmon_index/ \
    --output results/production

# ============================================================================
# STAGE 3: Differential Expression (EXTERNAL R)
# ============================================================================

# Module 6: Install R packages (first time only)
Rscript install_packages.R --test

# Module 6: Run DESeq2 analysis
Rscript run_deseq2.R \
    --counts results/production/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_deseq2 \
    --config recommendation.yaml \
    --condition treatment \
    --reference Control \
    --fdr 0.05 \
    --plots

# Optional: Run additional methods for ensemble
Rscript run_edger.R \
    --counts results/production/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_edger \
    --config recommendation.yaml

Rscript run_limma.R \
    --counts results/production/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_limma \
    --config recommendation.yaml

# ============================================================================
# STAGE 4: Import & Downstream Analysis (BACK TO RAPTOR Python)
# ============================================================================

# Module 7: Import DE results
raptor import-de \
    --de-file results/de_deseq2/de_results.csv \
    --summary results/de_deseq2/de_summary.json \
    --method deseq2

# Module 8: Optimize parameters
raptor optimize-params \
    --de-results results/de_deseq2/de_results.csv \
    --optimize fdr,lfc \
    --criterion reproducibility

# Module 9: Ensemble analysis (if multiple methods run)
raptor ensemble-analysis \
    --deseq2 results/de_deseq2/de_results.csv \
    --edger results/de_edger/de_results.csv \
    --limma results/de_limma/de_results.csv \
    --strategy consensus

# Module 10: Discover biomarkers
raptor discover-biomarkers \
    --de-results results/ensemble/ensemble_results.csv \
    --counts results/production/gene_counts.csv \
    --metadata data/metadata.csv \
    --rank-by significance,effect-size,consistency \
    --top-n 50
```

---

## 📊 DATA FLOW

```
FASTQ files
    ↓
[M1] Quick Count → quick_gene_counts.csv
    ↓
[M2] QC → outliers.txt, qc_report.html
    ↓
[M3] Profile → profile.json (32 features)
    ↓
[M4] Recommend → recommendation.yaml
    ↓
[M5] Production Pipeline → gene_counts.csv (final)
    ↓
═══════════════════════════════════════════════════════
    ↓
[M6] R Analysis → de_results.csv (EXTERNAL)
    ↓
[M7] Import → RAPTOR database
═══════════════════════════════════════════════════════
    ↓
[M8] Optimize → optimized_params.yaml
    ↓
[M9] Ensemble → ensemble_results.csv
    ↓
[M10] Biomarkers → ranked_genes.csv, biomarker_panels.txt
    ↓
Publication-ready results
```

---

## 📝 IMPLEMENTATION STATUS

| Module | Status | Language | Priority |
|--------|--------|----------|----------|
| M1: Quick Count | ✅ Complete | Python | - |
| M2: QC | ✅ Complete | Python | - |
| M3: Profile | ✅ Complete | Python | - |
| M4: Recommend | ✅ Complete | Python | - |
| M5: Production | ✅ Complete | Python | - |
| **M6: DE Analysis** | **✅ Complete** | **R** | **Done** |
| M7: Import | 🔄 Future | Python | High |
| M8: Optimize | 🔄 Future | Python | Medium |
| M9: Ensemble | 🔄 Future | Python | Medium |
| M10: Biomarkers | 🔄 Future | Python | Low |

**Legend:**
- ✅ Complete and tested
- 🔄 Future implementation
- ⚠️ In progress

---

## 🎯 MODULE 6 FOCUS

Since Module 6 is now complete, users can:

1. ✅ **Run DE analysis** using provided R scripts
2. ✅ **Use RAPTOR recommendations** (from M4)
3. ✅ **Generate standardized output** for M7 import
4. ✅ **Choose best method** (DESeq2/edgeR/limma)
5. ✅ **Create QC plots** for publication

**Files Provided:**
- ✅ run_deseq2.R (DESeq2 analysis)
- ✅ run_edger.R (edgeR analysis)
- ✅ run_limma.R (limma-voom analysis)
- ✅ install_packages.R (package installation)
- ✅ run_de_analysis.py (Python wrapper)
- ✅ MODULE6_README.md (complete documentation)

**Future Development:**
Modules 7-10 will be implemented to provide:
- Seamless import of M6 results
- Parameter optimization
- Method consensus
- Biomarker discovery

---

## 📧 CONTACT

**Author:** Ayeh Bolouki  
**Email:** ayehbolouki1988@gmail.com  
**Version:** 2.2.0

**Making free science for everybody around the world** 🌍
