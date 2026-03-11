# RAPTOR Example Scripts: Complete User Guide (Modules 1-9)

**Comprehensive Guide for the Complete RAPTOR RNA-seq Analysis Workflow**

This guide provides detailed instructions for using all RAPTOR example scripts from fast quantification to biomarker discovery.

---

## 📋 Table of Contents

### Overview
1. [RAPTOR Workflow Overview](#raptor-workflow-overview)
2. [Prerequisites & Installation](#prerequisites--installation)
3. [Quick Start: Complete Workflow](#quick-start-complete-workflow)

### Module Documentation
4. [Module 1: Quick Quantification](#module-1-quick-quantification)
5. [Module 2: Quality Assessment](#module-2-quality-assessment)
6. [Module 3: Data Profiling](#module-3-data-profiling)
7. [Module 4: Pipeline Recommendation](#module-4-pipeline-recommendation)
8. [Module 5: Production Pipeline](#module-5-production-pipeline)
9. [Module 6: R-based DE Analysis](#module-6-r-based-de-analysis)
10. [Module 7: DE Results Import](#module-7-de-results-import)
11. [Module 8: Parameter Optimization](#module-8-parameter-optimization)
12. [Module 9: Ensemble Analysis](#module-9-ensemble-analysis)
13. [Module 10: Biomarker Discovery](#module-10-biomarker-discovery)

### Reference
14. [File Formats & Data Structures](#file-formats--data-structures)
15. [Troubleshooting Guide](#troubleshooting-guide)
16. [Tips & Best Practices](#tips--best-practices)
17. [FAQ](#faq)

---

## RAPTOR Workflow Overview

### The Complete Pipeline

RAPTOR follows a **4-stage workflow** from raw FASTQ files to biomarker discovery:

```
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 1: FAST PROFILING                      │
│  FASTQ → Counts → QC → Profile → Recommendation (5-15 min)     │
└─────────────────────────────────────────────────────────────────┘
   M1: Quick Count    →  quick_gene_counts.csv
   M2: Quality QC     →  qc_results.json, counts_clean.csv
   M3: Data Profile   →  profile.json (32 features)
   M4: Recommend      →  recommendation.json (DESeq2/edgeR/limma)

┌─────────────────────────────────────────────────────────────────┐
│              STAGE 2: PRODUCTION QUANTIFICATION                 │
│  FASTQ → High-Quality Counts (30-120 min per sample)           │
└─────────────────────────────────────────────────────────────────┘
   M5: Production Pipeline  →  gene_counts.csv, transcripts.csv

┌─────────────────────────────────────────────────────────────────┐
│                   STAGE 3: DE ANALYSIS                          │
│  Counts → DE Results → Optimized → Consensus (R + Python)      │
└─────────────────────────────────────────────────────────────────┘
   M6: R DE Analysis  →  deseq2_results.csv, edger_results.csv
   M7: DE Import      →  de_result.pkl (standardized)
   M8: Optimize       →  optimized_parameters.json
   M9: Ensemble       →  consensus_genes.csv

┌─────────────────────────────────────────────────────────────────┐
│                  STAGE 4: BIOMARKER DISCOVERY                   │
│  Consensus Genes → Validated Biomarkers                        │
└─────────────────────────────────────────────────────────────────┘
   M10: Biomarkers    →  biomarker_panel.csv, validation_report.json
```

### Module Dependencies

```
M1 (Quick Count)
  └─→ M2 (Quality QC)
       └─→ M3 (Data Profile)
            └─→ M4 (Recommend)
                 └─→ M5 (Production Pipeline)
                      └─→ M6 (R DE Analysis)
                           └─→ M7 (DE Import)
                                ├─→ M8 (Parameter Optimization)
                                └─→ M9 (Ensemble Analysis)
                                     └─→ M10 (Biomarker Discovery)
```

### Key Design Principles

1. **Fast Profiling First** (M1-M4): Get recommendations in <15 minutes
2. **Demo Mode Available**: All modules work without real data
3. **Production Ready**: Full quantification when ready (M5)
4. **Multi-Method Consensus**: Combine results for high confidence (M7-M9)
4. **Validated Outputs**: Publication-quality results

---

## Prerequisites & Installation

### System Requirements

**Minimum:**
- Python 3.9+
- 8 GB RAM
- 20 GB disk space

**Recommended:**
- Python 3.10+
- 16 GB RAM
- 100 GB disk space (for real data)

### Python Dependencies

```bash
# Minimum (for demo mode)
pip install numpy pandas scipy

# Full installation (from RAPTOR root)
pip install -e .

# Optional for visualization
pip install matplotlib seaborn plotly
```

### External Tools (Optional for M5)

```bash
# Salmon (recommended)
conda install -c bioconda salmon

# Or download: https://github.com/COMBINE-lab/salmon

# Kallisto (alternative)
conda install -c bioconda kallisto

# STAR (for BAM alignment)
conda install -c bioconda star
```

### Verify Installation

```bash
# Check Python version
python --version  # Should be 3.9+

# Check RAPTOR installation
python -c "import raptor; print(raptor.__version__)"  # Should print 2.2.0

# Check external tools (if installed)
salmon --version
kallisto version
STAR --version
```

---

## Quick Start: Complete Workflow

### Demo Mode (No Real Data Needed!)

Try the complete workflow with simulated data:

```bash
# Navigate to examples directory
cd RAPTOR/examples/

# Stage 1: Fast Profiling (5 minutes)
python 01_quick_count.py --demo              # M1: Generate counts
python 02_quality_assessment.py --demo       # M2: Quality QC
python 03_data_profiler.py --demo           # M3: Extract features
python 04_recommender.py --demo             # M4: Get recommendation

# Stage 2: Production Pipeline (skipped in demo)
# python 05_production_pipeline.py --demo   # M5: Full quantification

# Stage 3: DE Analysis
python 07_DE_Import.py --demo               # M7: Import DE results
python 08_Parameter_Optimization.py --demo  # M8: Optimize parameters
python 09_Ensemble_Analysis.py --demo       # M9: Consensus genes

# Check all outputs
ls -R results/
```

### Real Data Workflow

```bash
# Stage 1: Fast Profiling
python 01_quick_count.py \
  --fastq-dir data/fastqs/ \
  --index data/salmon_index/

python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv

python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv

python 04_recommender.py \
  --profile results/profile.json

# Stage 2: Production (if recommended)
python 05_production_pipeline.py \
  --fastq-dir data/fastqs/ \
  --pipeline salmon \
  --index data/salmon_index/

# Stage 3: DE Analysis (after R analysis)
python 07_DE_Import.py \
  --de-file results/deseq2_results.csv \
  --pipeline deseq2

python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method fdr_control

python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method all
```

---

## Module 1: Quick Quantification

### Purpose

Fast RNA-seq quantification for profiling and QC. Uses lightweight pseudoalignment (Salmon/Kallisto) to generate quick count matrices in 5-15 minutes.

**Use this for:**
- Initial data assessment
- Quality control before full processing
- Getting quick recommendations
- Testing your workflow

**Don't use this for:**
- Final publication results (use M5 instead)
- When you need BAM files

### Quick Start

```bash
# Demo mode (no FASTQ files needed)
python 01_quick_count.py --demo

# With real FASTQ files
python 01_quick_count.py \
  --fastq-dir data/fastqs/ \
  --index data/salmon_index/ \
  --pipeline salmon
```

### Input Requirements

**Required:**
- FASTQ files (paired-end or single-end)
- Salmon/Kallisto index

**Optional:**
- Sample metadata CSV
- Custom reference annotation

### Output Files

```
results/quick_counts/
├── quick_gene_counts.csv       # Gene-level count matrix
├── quick_transcript_counts.csv # Transcript-level counts
├── sample_info.csv             # Sample metadata
├── quantification_stats.json   # QC metrics
└── logs/                       # Quantification logs
```

### Common Usage Patterns

#### Pattern 1: Basic Quantification

```bash
python 01_quick_count.py \
  --fastq-dir data/fastqs/ \
  --index data/salmon_index/
```

#### Pattern 2: With Metadata

```bash
python 01_quick_count.py \
  --fastq-dir data/fastqs/ \
  --index data/salmon_index/ \
  --metadata data/sample_metadata.csv
```

#### Pattern 3: Custom Output Location

```bash
python 01_quick_count.py \
  --fastq-dir data/fastqs/ \
  --index data/salmon_index/ \
  --output results/my_counts/
```

### Key Parameters

```
--fastq-dir PATH        Directory containing FASTQ files
--index PATH            Salmon/Kallisto index path
--pipeline {salmon,kallisto}  Quantification tool (default: salmon)
--threads N             Number of CPU threads (default: 4)
--output DIR            Output directory
--metadata FILE         Sample metadata CSV
--demo                  Run in demo mode
```

### What Happens Internally

1. **Sample Discovery**: Finds FASTQ files automatically
2. **Pseudoalignment**: Maps reads to transcriptome (5-15 min/sample)
3. **Quantification**: Estimates transcript abundances
4. **Gene Summarization**: Aggregates to gene-level
5. **QC Metrics**: Calculates mapping rates, library sizes

### Troubleshooting

**Issue: "No FASTQ files found"**
```bash
# Check your directory
ls data/fastqs/

# FASTQ files should end with .fastq, .fq, .fastq.gz, or .fq.gz
# Format: sample_1.fastq.gz, sample_2.fastq.gz (paired-end)
```

**Issue: "Index not found"**
```bash
# Build Salmon index first
salmon index -t transcripts.fa -i salmon_index

# Or download pre-built
wget https://refgenomes.databio.org/v3/assets/splash/...
```

**Issue: "Low mapping rate"**
- Check: Are you using the correct reference organism?
- Check: Is the index compatible with your FASTQ files?
- Check: Are FASTQ files corrupted?

### Next Steps

```bash
# After M1, run M2 for quality assessment
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv
```

---

## Module 2: Quality Assessment

### Purpose

Comprehensive quality control with 6-method outlier detection, batch effect assessment, and overall quality scoring (0-100).

**Use this for:**
- Identifying low-quality samples
- Detecting outliers before DE analysis
- Assessing batch effects
- Quality scoring for publication

**Key Features:**
- 6 outlier detection methods (PCA, Mahalanobis, IsolationForest, etc.)
- Library quality assessment
- Batch effect detection
- Overall quality score (0-100)

### Quick Start

```bash
# Demo mode
python 02_quality_assessment.py --demo

# With count matrix from M1
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv

# With metadata
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv \
  --metadata results/quick_counts/sample_info.csv
```

### Input Requirements

**Required:**
- Count matrix CSV (genes × samples)

**Optional:**
- Sample metadata CSV

### Output Files

```
results/qc/
├── qc_results.json          # Quality assessment report
├── counts_clean.csv         # Clean count matrix (outliers removed)
├── outlier_report.csv       # Detailed outlier analysis
├── quality_plots/           # QC visualizations
│   ├── pca_plot.png
│   ├── library_sizes.png
│   └── correlation_heatmap.png
└── batch_effects.json       # Batch effect analysis
```

### Common Usage Patterns

#### Pattern 1: Basic QC

```bash
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv
```

#### Pattern 2: With Outlier Removal

```bash
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv \
  --remove-outliers \
  --outlier-threshold 2.5
```

#### Pattern 3: Custom Normalization

```bash
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv \
  --normalization log2  # or cpm, quantile
```

### Key Parameters

```
--counts FILE           Count matrix CSV
--metadata FILE         Sample metadata CSV
--normalization METHOD  log2, cpm, quantile, none (default: log2)
--remove-outliers       Remove detected outliers
--outlier-threshold N   Z-score threshold (default: 3.0)
--output DIR            Output directory
--demo                  Run in demo mode
```

### Understanding Quality Scores

**Overall Quality Score (0-100):**
- **90-100**: Excellent - Ready for publication
- **80-89**: Good - Suitable for most analyses
- **70-79**: Fair - Acceptable with caveats
- **60-69**: Poor - Consider filtering samples
- **<60**: Very Poor - Significant issues

**Component Scores:**
- Library Quality (30%): Library size consistency
- Gene Detection (30%): Number of expressed genes
- Technical Quality (20%): Outliers, technical noise
- Biological Signal (20%): Signal-to-noise ratio

### What Happens Internally

1. **Input Validation**: Checks count matrix format and structure
2. **Library Quality**: Assesses library size distribution
3. **Gene Detection**: Counts detected genes per sample
4. **Outlier Detection**: 6 independent methods
5. **Batch Effects**: Detects systematic biases
6. **Biological Signal**: Estimates signal-to-noise
7. **Overall Scoring**: Combines components → 0-100 score

### Troubleshooting

**Issue: "Quality score is low (<70)"**
```bash
# Check the issues in qc_results.json
cat results/qc/qc_results.json | grep -A 10 "issues"

# Common issues:
# - Low library sizes → May need more sequencing
# - Many outliers → Check sample quality
# - Batch effects → Use batch correction (limma, DESeq2)
```

**Issue: "Too many outliers detected"**
```bash
# Relax the threshold
python 02_quality_assessment.py \
  --counts results/quick_counts/quick_gene_counts.csv \
  --outlier-threshold 3.5  # More lenient

# Or review individual samples
cat results/qc/outlier_report.csv
```

### Next Steps

```bash
# After M2, run M3 for data profiling
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv
```

---

## Module 3: Data Profiling

### Purpose

Extract 32 statistical features for machine learning-based pipeline recommendation. Profiles sample characteristics, library metrics, expression patterns, and variance structure.

**Use this for:**
- Getting ML-based pipeline recommendations
- Understanding your dataset characteristics
- Comparing datasets quantitatively

**32 Features Include:**
- Sample characteristics (n_samples, n_groups, group sizes)
- Library metrics (mean, CV, range)
- Expression patterns (sparsity, zero fraction)
- Dispersion estimates (BCV, common dispersion)
- Quality indicators (outliers, batch effects)

### Quick Start

```bash
# Demo mode
python 03_data_profiler.py --demo

# With clean counts from M2
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv

# With metadata
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv \
  --metadata results/quick_counts/sample_info.csv
```

### Input Requirements

**Required:**
- Count matrix CSV (genes × samples)

**Optional:**
- Sample metadata CSV with group/condition column

### Output Files

```
results/
├── profile.json              # 32 statistical features
└── profiling_report.txt      # Human-readable summary
```

### Common Usage Patterns

#### Pattern 1: Basic Profiling

```bash
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv
```

#### Pattern 2: With Group Information

```bash
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv \
  --metadata data/sample_metadata.csv \
  --group-column condition
```

#### Pattern 3: Custom Output

```bash
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv \
  --output results/my_profile/
```

### Key Parameters

```
--counts FILE          Count matrix CSV
--metadata FILE        Sample metadata CSV
--group-column NAME    Column for grouping (default: condition)
--output DIR           Output directory
--demo                 Run in demo mode
```

### Understanding the 32 Features

**Sample Characteristics:**
- `n_samples`: Total number of samples
- `n_genes`: Total number of genes
- `n_groups`: Number of experimental groups
- `min_group_size`: Smallest group size
- `max_group_size`: Largest group size

**Library Metrics:**
- `library_size_mean`: Average library size
- `library_size_cv`: Library size coefficient of variation
- `library_size_range`: Max/min library size ratio

**Expression Patterns:**
- `zero_fraction`: Proportion of zero counts
- `low_count_proportion`: Proportion of low counts (<10)
- `detection_rate`: Proportion of detected genes

**Dispersion (CRITICAL for DE method selection):**
- `bcv`: Biological coefficient of variation
- `bcv_category`: low (<0.2), moderate (0.2-0.4), high (>0.4)
- `common_dispersion`: Shared dispersion across genes

**Quality Indicators:**
- `has_outliers`: Boolean outlier flag
- `has_batch_effect`: Boolean batch effect flag
- `design_complexity`: simple/moderate/complex

### What Happens Internally

1. **Input Validation**: Checks matrix and metadata
2. **Sample Stats**: Counts samples, groups, sizes
3. **Library Analysis**: Calculates size metrics
4. **Expression Profiling**: Analyzes count distribution
5. **Dispersion Estimation**: Computes BCV and common dispersion
6. **Quality Checks**: Detects outliers, batch effects
7. **Feature Extraction**: Exports all 32 features to JSON

### Troubleshooting

**Issue: "Small sample size warning"**
```bash
# Profile still works, but recommendations may be limited
# Consider:
# - Pooling groups if scientifically appropriate
# - Using more conservative DE methods (edgeR)
```

**Issue: "High BCV detected"**
```bash
# BCV > 0.4 indicates high biological variability
# Recommended methods: edgeR or edgeR_robust
# See profile.json for bcv value
cat results/profile.json | grep "bcv"
```

### Next Steps

```bash
# After M3, run M4 for pipeline recommendation
python 04_recommender.py \
  --profile results/profile.json
```

---

## Module 4: Pipeline Recommendation

### Purpose

ML-based recommendation for optimal DE analysis pipeline (DESeq2, edgeR, limma-voom, or Wilcoxon). Uses data profile features to score and rank methods.

**Recommends based on:**
- Sample size and group sizes
- Biological coefficient of variation (BCV)
- Library characteristics
- Presence of outliers or batch effects
- Overall data quality

**Output:**
- Primary recommendation with score
- Alternative methods with scores
- Detailed reasoning
- Pipeline-specific R code templates

### Quick Start

```bash
# Demo mode
python 04_recommender.py --demo

# With profile from M3
python 04_recommender.py \
  --profile results/profile.json

# Get full comparison
python 04_recommender.py \
  --profile results/profile.json \
  --show-all-scores
```

### Input Requirements

**Required:**
- Profile JSON from Module 3 (or data profile object)

### Output Files

```
results/recommendation/
├── recommendation.json          # Full recommendation details
├── recommendation_report.txt    # Human-readable report
├── method_comparison.csv        # All methods scored
└── r_code_templates/            # Pipeline-specific R scripts
    ├── deseq2_template.R
    ├── edger_template.R
    └── limma_template.R
```

### Common Usage Patterns

#### Pattern 1: Basic Recommendation

```bash
python 04_recommender.py \
  --profile results/profile.json
```

#### Pattern 2: See All Scores

```bash
python 04_recommender.py \
  --profile results/profile.json \
  --show-all-scores
```

#### Pattern 3: Export R Code

```bash
python 04_recommender.py \
  --profile results/profile.json \
  --export-r-code
```

### Key Parameters

```
--profile FILE         Profile JSON from M3
--show-all-scores      Display all method scores
--export-r-code        Generate R code templates
--output DIR           Output directory
--demo                 Run in demo mode
```

### Understanding Recommendations

**Scoring Algorithm:**

1. **Sample Size (MOST CRITICAL):**
   - `min_group_size < 3` → DESeq2 (+30 bonus)
   - `3-8` → DESeq2/edgeR
   - `8-20` → limma
   - `≥20` → Wilcoxon

2. **BCV (Biological Variation):**
   - `< 0.2` → limma preferred
   - `0.2-0.4` → DESeq2/edgeR
   - `> 0.4` → edgeR_robust

3. **Low Count Proportion:**
   - `> 30%` → edgeR (+10)

4. **Outliers:**
   - Present → edgeR_robust (+20)

5. **Batch Effects:**
   - Present → DESeq2/limma (+10)

**Method Profiles:**

| Method | Best For | Strengths | Weaknesses |
|--------|----------|-----------|------------|
| **DESeq2** | Small samples (n<3 per group) | Robust with small n, handles batch effects | Slower for large datasets |
| **edgeR** | Moderate BCV (0.2-0.4) | Fast, handles low counts well | Less robust to outliers |
| **edgeR_robust** | High BCV, outliers | Robust to outliers and high variability | More conservative |
| **limma-voom** | Large samples (n≥8), batch effects | Fast, flexible, batch correction | Needs adequate replicates |
| **Wilcoxon** | Large samples (n≥20), non-parametric | No distribution assumptions | Less powerful than parametric methods |

### What Happens Internally

1. **Load Profile**: Reads 32 features from JSON
2. **Feature Validation**: Checks critical features present
3. **Score Calculation**: Applies scoring rules to each method
4. **Ranking**: Orders methods by total score
5. **Reasoning Generation**: Explains recommendation
6. **R Code Generation**: Creates method-specific templates

### Troubleshooting

**Issue: "Profile missing critical features"**
```bash
# Re-run M3 to regenerate complete profile
python 03_data_profiler.py \
  --counts results/qc/counts_clean.csv \
  --metadata data/sample_metadata.csv
```

**Issue: "Multiple methods have similar scores"**
```bash
# This is normal! Check reasoning for each
cat results/recommendation/recommendation_report.txt

# Consider running both top methods and using ensemble (M9)
```

### Next Steps

```bash
# After M4, you have two paths:

# Path A: Production quantification (recommended)
python 05_production_pipeline.py \
  --fastq-dir data/fastqs/ \
  --pipeline star_salmon

# Path B: Skip to DE analysis if counts are adequate
# Use recommended method in R (see r_code_templates/)
```

---

## Module 5: Production Pipeline

### Purpose

Full production-quality RNA-seq quantification with multiple pipeline options. Provides high-quality gene and transcript counts for final DE analysis.

**Use this when:**
- You need publication-quality results
- You need BAM files for visualization
- Quick counts (M1) are insufficient
- You have time for full processing (30-120 min/sample)

**6 Pipeline Options:**
1. `salmon` - Fast, gene+transcript, bootstraps ⭐
2. `kallisto` - Fastest, gene+transcript, bootstraps
3. `star_featurecounts` - BAM files, gene-level only
4. `hisat2_featurecounts` - BAM files, lower memory
5. `star_rsem` - BAM files, gene+transcript, slower
6. `star_salmon` - BAM files, gene+transcript, bootstraps 🌟

### Quick Start

```bash
# Demo mode
python 05_production_pipeline.py --demo

# With real FASTQ files (salmon - fastest)
python 05_production_pipeline.py \
  --fastq-dir data/fastqs/ \
  --pipeline salmon \
  --index data/salmon_index/

# With STAR (produces BAM files)
python 05_production_pipeline.py \
  --fastq-dir data/fastqs/ \
  --pipeline star_featurecounts \
  --genome-index data/star_index/ \
  --annotation data/genes.gtf
```

### Input Requirements

**Required:**
- FASTQ files (paired-end or single-end)
- Pipeline-specific index/references

**For Salmon/Kallisto:**
- Transcriptome index

**For STAR/HISAT2:**
- Genome index
- Gene annotation (GTF)

### Output Files

```
results/production_pipeline/
├── gene_counts.csv              # Gene-level count matrix
├── transcript_counts.csv        # Transcript-level counts (if applicable)
├── sample_results/              # Per-sample outputs
│   ├── sample1/
│   │   ├── quant.sf            # Salmon/Kallisto quantification
│   │   ├── aligned.bam         # BAM file (STAR/HISAT2 pipelines)
│   │   └── logs/
│   └── sample2/
├── pipeline_stats.json          # Overall statistics
└── qc_report.html              # MultiQC report (if available)
```

### Common Usage Patterns

#### Pattern 1: Salmon (Recommended - Fast & Accurate)

```bash
python 05_production_pipeline.py \
  --fastq-dir data/fastqs/ \
  --pipeline salmon \
  --index data/salmon_index/ \
  --threads 8
```

#### Pattern 2: STAR + Salmon (BAM + Transcripts)

```bash
python 05_production_pipeline.py \
  --fastq-dir data/fastqs/ \
  --pipeline star_salmon \
  --genome-index data/star_index/ \
  --annotation data/genes.gtf \
  --threads 12
```

#### Pattern 3: Reuse Quick Counts (Instant!)

```bash
python 05_production_pipeline.py \
  --use-quick-counts results/quick_counts/ \
  --pipeline salmon
```

### Key Parameters

```
--fastq-dir DIR         FASTQ files directory
--pipeline PIPELINE     salmon, kallisto, star_featurecounts, 
                       hisat2_featurecounts, star_rsem, star_salmon
--index DIR            Salmon/Kallisto index
--genome-index DIR     STAR/HISAT2 genome index
--annotation FILE      GTF annotation file
--threads N            CPU threads (default: 4)
--output DIR           Output directory
--use-quick-counts     Reuse M1 counts (instant)
--demo                 Run in demo mode
```

### Pipeline Comparison

| Pipeline | BAM | Gene | Transcript | Bootstraps | Memory | Time | Best For |
|----------|-----|------|------------|------------|--------|------|----------|
| salmon | ❌ | ✅ | ✅ | ✅ | 8GB | 10-20min | Speed + accuracy ⭐ |
| kallisto | ❌ | ✅ | ✅ | ✅ | 4GB | 5-10min | Maximum speed |
| star_featurecounts | ✅ | ✅ | ❌ | ❌ | 32GB | 40-70min | BAM files needed |
| hisat2_featurecounts | ✅ | ✅ | ❌ | ❌ | 16GB | 30-60min | Lower memory |
| star_rsem | ✅ | ✅ | ✅ | ❌ | 32GB | 60-120min | Gene + transcript + BAM |
| star_salmon | ✅ | ✅ | ✅ | ✅ | 32GB | 50-90min | Everything! 🌟 |

### What Happens Internally

**Salmon Pipeline:**
1. Read validation
2. Pseudoalignment to transcriptome
3. EM algorithm for abundance estimation
4. Bootstrap sampling for uncertainty
5. Gene-level aggregation

**STAR Pipeline:**
1. Read trimming (optional)
2. Genome alignment
3. BAM file generation
4. Feature counting (featureCounts/RSEM/Salmon)
5. Quality metrics

### Troubleshooting

**Issue: "Pipeline not found"**
```bash
# Install required tool
conda install -c bioconda salmon  # or star, hisat2, kallisto
```

**Issue: "Out of memory"**
```bash
# Use lower-memory pipeline
python 05_production_pipeline.py \
  --pipeline kallisto  # Only 4GB needed
  
# Or HISAT2 instead of STAR
python 05_production_pipeline.py \
  --pipeline hisat2_featurecounts  # 16GB vs 32GB for STAR
```

**Issue: "Too slow"**
```bash
# Use more threads
python 05_production_pipeline.py \
  --pipeline salmon \
  --threads 16  # Increase from default 4

# Or use faster pipeline
python 05_production_pipeline.py \
  --pipeline kallisto  # Fastest option
```

### Next Steps

```bash
# After M5, run DE analysis in R (M6)
# Then import results (M7)
python 07_DE_Import.py \
  --de-file results/deseq2_results.csv \
  --pipeline deseq2
```

---

## Module 6: R-based DE Analysis

### Purpose

Perform differential expression analysis in R using the recommended method from M4 (DESeq2, edgeR, or limma).

**Note:** This module consists of R scripts, not Python. See `examples/R/` directory.

### Available R Scripts

```
examples/R/
├── deseq2_analysis.R           # DESeq2 workflow
├── edger_analysis.R            # edgeR workflow
├── limma_voom_analysis.R       # limma-voom workflow
└── helper_functions.R          # Shared functions
```

### Quick Start

```bash
# In R console
source("examples/R/deseq2_analysis.R")

# Or command line
Rscript examples/R/deseq2_analysis.R \
  --counts results/production_pipeline/gene_counts.csv \
  --metadata data/sample_metadata.csv \
  --output results/de_analysis/
```

### Input Requirements

**Required:**
- Gene count matrix (from M5 or M1)
- Sample metadata with condition/group column

### Output Files

```
results/de_analysis/
├── deseq2_results.csv          # Full DE results
├── deseq2_significant.csv      # Significant genes only
├── normalized_counts.csv       # Normalized count matrix
└── plots/
    ├── ma_plot.png
    ├── volcano_plot.png
    └── heatmap.png
```

### Common Usage Patterns

#### Pattern 1: DESeq2 (Recommended for Small Samples)

```r
library(DESeq2)

# Load data
counts <- read.csv("results/production_pipeline/gene_counts.csv", row.names=1)
metadata <- read.csv("data/sample_metadata.csv")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Run analysis
dds <- DESeq(dds)
results <- results(dds, alpha=0.05)

# Save results
write.csv(as.data.frame(results), "results/deseq2_results.csv")
```

#### Pattern 2: edgeR (Recommended for Moderate BCV)

```r
library(edgeR)

# Load data
counts <- read.csv("results/production_pipeline/gene_counts.csv", row.names=1)
metadata <- read.csv("data/sample_metadata.csv")

# Create DGEList
dge <- DGEList(counts=counts, group=metadata$condition)

# Normalization
dge <- calcNormFactors(dge)

# Dispersion estimation
design <- model.matrix(~condition, data=metadata)
dge <- estimateDisp(dge, design)

# DE testing
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef=2)

# Save results
write.csv(topTags(qlf, n=Inf), "results/edger_results.csv")
```

#### Pattern 3: limma-voom (Recommended for Large Samples)

```r
library(limma)
library(edgeR)

# Load data
counts <- read.csv("results/production_pipeline/gene_counts.csv", row.names=1)
metadata <- read.csv("data/sample_metadata.csv")

# Create DGEList and normalize
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~condition, data=metadata)

# Voom transformation
v <- voom(dge, design, plot=FALSE)

# Fit model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Save results
results <- topTable(fit, coef=2, number=Inf)
write.csv(results, "results/limma_results.csv")
```

### Next Steps

```bash
# After R analysis (M6), import results to RAPTOR (M7)
python 07_DE_Import.py \
  --de-file results/de_analysis/deseq2_results.csv \
  --pipeline deseq2
```

---

## Module 7: DE Results Import

### Purpose

Import and standardize DE results from R (DESeq2, edgeR, limma) into RAPTOR's unified DEResult format. Enables downstream optimization (M8) and ensemble analysis (M9).

**Key Features:**
- Auto-detects DE pipeline (DESeq2/edgeR/limma)
- Standardizes column names
- Validates results
- Computes summary statistics
- Saves as pickle for M8/M9

### Quick Start

```bash
# Demo mode
python 07_DE_Import.py --demo

# Import DESeq2 results
python 07_DE_Import.py \
  --de-file results/de_analysis/deseq2_results.csv \
  --pipeline auto  # Auto-detects from column names

# Import edgeR results
python 07_DE_Import.py \
  --de-file results/de_analysis/edger_results.csv \
  --pipeline edger
```

### Input Requirements

**Required:**
- DE results CSV from R (DESeq2, edgeR, or limma)

**Expected Columns:**

**DESeq2:**
- Gene IDs (row names or 'gene_id' column)
- `log2FoldChange` or `lfc` or `log2FC`
- `pvalue` or `pval`
- `padj` or `FDR`

**edgeR:**
- `logFC`
- `PValue`
- `FDR`

**limma:**
- `logFC`
- `P.Value`
- `adj.P.Val`

### Output Files

```
results/de_imported/
├── de_result.pkl               # DEResult object (for M8/M9)
├── de_standardized.csv         # Standardized column names
├── de_significant.csv          # Significant genes only
├── de_summary.json             # Summary statistics
└── import_report.txt           # Import details
```

### Common Usage Patterns

#### Pattern 1: Auto-Detection (Recommended)

```bash
python 07_DE_Import.py \
  --de-file results/deseq2_results.csv \
  --pipeline auto
```

#### Pattern 2: Specify Pipeline

```bash
python 07_DE_Import.py \
  --de-file results/edger_results.csv \
  --pipeline edger \
  --output results/de_imported/edger/
```

#### Pattern 3: Custom Thresholds

```bash
python 07_DE_Import.py \
  --de-file results/limma_results.csv \
  --pipeline limma \
  --fdr-threshold 0.01 \
  --lfc-threshold 1.0
```

### Key Parameters

```
--de-file FILE          DE results CSV from R
--pipeline PIPELINE     auto, deseq2, edger, limma (default: auto)
--fdr-threshold FLOAT   FDR cutoff (default: 0.05)
--lfc-threshold FLOAT   Log2 fold change cutoff (default: 0.0)
--output DIR            Output directory
--demo                  Run in demo mode
```

### DEResult Structure

The standardized DEResult object has:

```python
@dataclass
class DEResult:
    results_df: pd.DataFrame     # Standardized results
    pipeline: str                # 'DESEQ2', 'EDGER', 'LIMMA-VOOM'
    parameters: Dict             # FDR/LFC thresholds
    metadata: Dict               # Import info
    
    # Properties
    n_genes: int                 # Total genes
    n_significant: int           # Significant genes
    n_up: int                    # Upregulated
    n_down: int                  # Downregulated
    
    # Methods
    def get_top_genes(n=50)
    def filter_by_threshold(fdr, lfc)
    def summary()
```

**Standardized Columns:**
- `gene_id` (index)
- `log2_fold_change`
- `p_value`
- `adjusted_p_value`
- `base_mean` (if available)

### What Happens Internally

1. **File Loading**: Reads CSV, handles various formats
2. **Pipeline Detection**: Auto-detects from column names
3. **Column Mapping**: Maps to standardized names
4. **Validation**: Checks required columns present
5. **Gene ID Handling**: Ensures gene_id is index
6. **Statistics**: Computes n_genes, n_significant, etc.
7. **Export**: Saves DEResult pickle + CSV files

### Troubleshooting

**Issue: "Could not detect pipeline"**
```bash
# Specify manually
python 07_DE_Import.py \
  --de-file results/my_results.csv \
  --pipeline deseq2  # or edger, limma
```

**Issue: "Missing required columns"**
```bash
# Check what columns you have
head results/my_results.csv

# Required for DESeq2: log2FoldChange, pvalue, padj
# Required for edgeR: logFC, PValue, FDR
# Required for limma: logFC, P.Value, adj.P.Val
```

**Issue: "Gene IDs not found"**
```bash
# Make sure gene IDs are either:
# 1. Row names (no column name)
# 2. Column named 'gene_id', 'GeneID', 'Gene', or 'ID'

# Fix in R before exporting:
write.csv(results, "results.csv", row.names=TRUE)
```

### Next Steps

```bash
# After M7, you can run:

# M8: Parameter optimization
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method fdr_control

# M9: Ensemble analysis (if you have multiple methods)
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method all
```

---

## Module 8: Parameter Optimization

### Purpose

Optimize FDR (alpha) and log fold change thresholds for DE analysis using statistical methods. No validation data needed!

**4 Optimization Methods:**
1. **FDR Control** - Statistical FDR estimation ⭐ (No validation needed)
2. **Ground Truth** - With validated genes (Gold standard)
3. **Stability** - Bootstrap-based (Needs count matrix)
4. **Reproducibility** - Cross-cohort (Needs 2 independent datasets)

**Optimizes:**
- `alpha`: FDR/p-value threshold (0.001 - 0.20)
- `lfc_threshold`: Minimum log2 fold change (0.0 - 2.0)

### Quick Start

```bash
# Demo mode (all 4 methods)
python 08_Parameter_Optimization.py --demo --method all

# FDR Control (recommended - no validation needed!)
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method fdr_control

# With validated genes (gold standard)
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method ground_truth \
  --ground-truth data/validated_genes.csv
```

### Input Requirements

**For FDR Control (Method 1):**
- DE results from M7 (de_result.pkl)

**For Ground Truth (Method 2):**
- DE results + validated genes CSV

**For Stability (Method 3):**
- DE results + count matrix + sample metadata

**For Reproducibility (Method 4):**
- Two independent DE results

### Output Files

```
results/parameter_optimization/
├── fdr_control/
│   └── optimization_results.json
├── ground_truth/
│   └── optimization_results.json
├── stability/
│   └── optimization_results.json
├── reproducibility/
│   └── optimization_results.json
├── method_comparison.csv          # If multiple methods used
└── consensus_parameters.json      # Median across methods
```

### Common Usage Patterns

#### Pattern 1: FDR Control (Recommended - Works for Everyone!)

```bash
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method fdr_control \
  --strategy grid \
  --grid-points 5
```

#### Pattern 2: Ground Truth (With Validated Genes)

```bash
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method ground_truth \
  --ground-truth data/validated_genes.csv

# validated_genes.csv format:
# gene_id
# ENSG00000000001
# ENSG00000000002
# ...
```

#### Pattern 3: Stability (Bootstrap-Based)

```bash
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method stability \
  --counts results/production_pipeline/gene_counts.csv \
  --metadata data/sample_metadata.csv \
  --n-bootstrap 100
```

#### Pattern 4: All Methods (Most Robust)

```bash
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method all \
  --ground-truth data/validated_genes.csv \
  --counts results/production_pipeline/gene_counts.csv \
  --metadata data/sample_metadata.csv
```

### Key Parameters

```
--de-result FILE        DE result from M7 (.pkl)
--method METHOD         fdr_control, ground_truth, stability, 
                       reproducibility, or all
--strategy STRATEGY     grid, random, differential_evolution
--grid-points N         Grid points per parameter (default: 5)
--ground-truth FILE     Validated genes CSV
--counts FILE           Count matrix CSV (for stability)
--metadata FILE         Sample metadata (for stability)
--cohort2 FILE          Second cohort DE results (for reproducibility)
--n-bootstrap N         Bootstrap iterations (default: 100)
--output DIR            Output directory
--demo                  Run in demo mode
```

### Method Selection Guide

**Decision Tree:**

```
Do you have validated genes from qPCR/literature?
  └─ YES → Use ground_truth (gold standard)

Do you have 2 independent cohorts?
  └─ YES → Use reproducibility (highest confidence)

Do you have the original count matrix?
  └─ YES → Use stability (bootstrap-based)

None of the above?
  └─ Use fdr_control (works for everyone!) ⭐
```

### Search Strategies

**Grid Search (Default - Recommended):**
- Tests all combinations: grid_points²
- 5 points = 25 combinations (5×5)
- Guaranteed to find optimum
- Time: 1-5 minutes

**Random Search:**
- Monte Carlo sampling
- Faster for large spaces
- Good approximation

**Differential Evolution:**
- Population-based global optimization
- Best for complex optimization landscapes

### What Happens Internally

**FDR Control Method:**
1. Loads DE results (DataFrame)
2. Estimates π₀ (proportion of true nulls)
3. Tests parameter combinations
4. Optimizes for FDR control accuracy
5. Returns best alpha and LFC threshold

**Ground Truth Method:**
1. Loads DE results and validated genes
2. Tests parameter combinations
3. Calculates precision, recall, F1 score
4. Optimizes for maximum F1 score
5. Returns best parameters

### Understanding Results

```json
{
  "best_parameters": {
    "alpha": 0.0350,
    "lfc_threshold": 0.75
  },
  "best_score": 0.8934,
  "optimization_method": "fdr_control",
  "search_strategy": "grid",
  "n_iterations": 25
}
```

**Use these parameters in your R analysis:**
```r
# In R
alpha_optimized <- 0.0350
lfc_threshold <- 0.75

# Filter results
sig_genes <- results(dds, alpha=alpha_optimized) %>%
  filter(padj < alpha_optimized, abs(log2FoldChange) > lfc_threshold)
```

### Troubleshooting

**Issue: "Not enough validated genes"**
```bash
# Need at least 10 validated genes for ground_truth
# Consider using fdr_control instead
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl \
  --method fdr_control
```

**Issue: "Optimization takes too long"**
```bash
# Reduce grid points
python 08_Parameter_Optimization.py \
  --method fdr_control \
  --grid-points 3  # 3×3 = 9 combinations

# Or use random search
python 08_Parameter_Optimization.py \
  --method fdr_control \
  --strategy random \
  --n-iterations 20
```

### Next Steps

```bash
# After M8, apply parameters in R or continue to M9
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method all
```

---

## Module 9: Ensemble Analysis

### Purpose

Combine DE results from multiple methods (DESeq2, edgeR, limma) to create high-confidence consensus gene lists using statistical ensemble techniques.

**5 Ensemble Methods:**
1. **Voting** - Count-based (high confidence) ⭐
2. **Fisher's** - P-value combination (maximum sensitivity)
3. **Brown's** - Correlation-aware combination
4. **RRA** - Robust Rank Aggregation (handles outliers)
5. **Weighted** - Performance-weighted (with validation data)

**IMPORTANT:** Includes automatic adapter to bridge DEResult (M7) with ensemble module.

### Quick Start

```bash
# Demo mode (all 5 methods)
python 09_Ensemble_Analysis.py --demo --method all

# Voting ensemble (high confidence)
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/deseq2_result.pkl \
               results/de_imported/edger_result.pkl \
               results/de_imported/limma_result.pkl \
  --method voting \
  --min-methods 3

# Fisher's method (maximum sensitivity)
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method fisher
```

### Input Requirements

**Required:**
- 2+ DE results from M7 (.pkl files)

**Note:** The script automatically adapts DEResult objects for the ensemble module. No manual conversion needed!

### Output Files

```
results/ensemble_analysis/
├── voting/
│   └── consensus_genes.csv
├── fisher/
│   └── consensus_genes.csv
├── brown/
│   └── consensus_genes.csv
├── rra/
│   └── consensus_genes.csv
├── weighted/
│   └── consensus_genes.csv
├── method_comparison.csv          # Compares all methods
└── ensemble_summary.json          # Overall statistics
```

### Common Usage Patterns

#### Pattern 1: Voting (High Confidence)

```bash
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/deseq2_result.pkl \
               results/de_imported/edger_result.pkl \
               results/de_imported/limma_result.pkl \
  --method voting \
  --min-methods 3  # All 3 must agree
```

#### Pattern 2: Fisher's Method (Exploratory)

```bash
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method fisher \
  --significance 0.05
```

#### Pattern 3: All Methods (Most Comprehensive)

```bash
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method all
```

#### Pattern 4: Weighted (With Performance Data)

```bash
python 09_Ensemble_Analysis.py \
  --de-results results/de_imported/*.pkl \
  --method weighted \
  --weights '{"DESeq2": 1.0, "edgeR": 0.8, "limma": 0.9}'
```

### Key Parameters

```
--de-results FILES      DE result files from M7 (.pkl)
--method METHOD         voting, fisher, brown, rra, weighted, or all
--min-methods N         Minimum methods for voting (default: 2)
--significance FLOAT    Significance threshold (default: 0.05)
--weights JSON          Method weights (for weighted ensemble)
--output DIR            Output directory
--demo                  Run in demo mode
```

### Method Comparison

| Method | Stringency | Best For | Requires |
|--------|-----------|----------|----------|
| **Voting** | HIGH | Publication-quality, high confidence | 2+ methods |
| **Fisher's** | LOW | Exploratory, finding weak signals | 2+ methods |
| **Brown's** | MEDIUM | Same data/normalization | 2+ methods |
| **RRA** | MEDIUM | Robust to outliers, balanced | 2+ methods |
| **Weighted** | VARIABLE | Custom weighting by performance | Validation data |

### Method Selection Guide

```
Want maximum confidence?
  └─ Use Voting (min_methods = all)

Want to find weak signals?
  └─ Use Fisher's (most sensitive)

Methods use same data?
  └─ Use Brown's (accounts for correlation)

Want balanced approach?
  └─ Use RRA (robust to outliers)

Have performance data?
  └─ Use Weighted (customizable)
```

### What Happens Internally

**CRITICAL - Automatic Adapter:**

The script includes `adapt_deresult_for_ensemble()` which automatically:
1. Converts `.results_df` → `.data`
2. Renames `adjusted_p_value` → `padj`
3. Renames `log2_fold_change` → `log2FoldChange`
4. Renames `p_value` → `pvalue`
5. Converts gene_id from index to column

**You don't need to do anything!** Just pass the .pkl files.

**Voting Method:**
1. Loads and adapts DEResult objects
2. Counts how many methods detect each gene
3. Requires ≥min_methods for consensus
4. Checks direction consistency
5. Outputs high-confidence genes

**Fisher's Method:**
1. Loads and adapts DEResult objects
2. Combines p-values using Fisher's method
3. Calculates meta-analytic effect sizes
4. Applies significance threshold
5. Outputs consensus with combined p-values

### Understanding Results

**Voting Output:**
```csv
gene_id,n_votes,direction,direction_agreement,methods_detected
ENSG00000000001,3,up,1.0,"DESeq2,edgeR,limma"
ENSG00000000002,3,down,1.0,"DESeq2,edgeR,limma"
ENSG00000000003,2,up,1.0,"DESeq2,edgeR"
```

**Fisher's Output:**
```csv
gene_id,combined_pvalue,combined_padj,direction,meta_lfc
ENSG00000000001,1.2e-15,3.5e-12,up,2.34
ENSG00000000002,2.3e-12,4.6e-09,down,-1.87
```

### Troubleshooting

**Issue: "AttributeError: 'DEResult' has no attribute 'data'"**
```bash
# This should NOT happen with improved script!
# The adapter handles this automatically

# If you see this error, you're using the old script
# Use: 09_Ensemble_Analysis_IMPROVED.py
```

**Issue: "Only 1 DE result provided"**
```bash
# Ensemble needs at least 2 methods
# Run another DE method in R, then import (M7)

# Or use single-method parameter optimization (M8)
python 08_Parameter_Optimization.py \
  --de-result results/de_imported/de_result.pkl
```

**Issue: "Methods disagree on direction"**
```bash
# This is flagged in direction_consistency column
# Review individual methods carefully
# Consider biological plausibility

# Or use more stringent voting
python 09_Ensemble_Analysis.py \
  --method voting \
  --min-methods 3  # All must agree
```

### Next Steps

```bash
# After M9, use consensus genes for:

# 1. Publication (use voting results)
# 2. Biomarker discovery (M10)
# 3. Pathway analysis
# 4. Validation experiments

# Example: Get top genes for validation
head -n 20 results/ensemble_analysis/voting/consensus_genes.csv
```

---

## Module 10: Biomarker Discovery

### Purpose

Identify and validate potential biomarker panels from consensus genes using machine learning and statistical validation.

**Note:** Module 10 scripts are under development. Check `examples/` directory for latest versions.

### Expected Features

- Feature selection from consensus genes
- Classification model training
- Cross-validation
- Performance metrics (AUC, sensitivity, specificity)
- Biomarker panel optimization
- Validation recommendations

---

## File Formats & Data Structures

### Count Matrix Format

**Gene × Sample Matrix (CSV):**

```csv
,Sample1,Sample2,Sample3,Sample4
ENSG00000000001,1523,1845,1632,1789
ENSG00000000002,234,189,256,198
ENSG00000000003,0,2,1,0
ENSG00000000004,8456,9123,8734,8967
```

**Requirements:**
- Row names: Gene IDs
- Column names: Sample IDs
- Values: Raw counts (integers)
- No missing values

### Sample Metadata Format

**CSV with Sample Information:**

```csv
sample_id,condition,batch,replicate
Sample1,Control,Batch1,1
Sample2,Control,Batch1,2
Sample3,Treatment,Batch2,1
Sample4,Treatment,Batch2,2
```

**Requirements:**
- `sample_id` column (must match count matrix columns)
- `condition` or `group` column (for grouping)
- Optional: `batch`, `replicate`, etc.

### DE Results Format

**DESeq2:**
```csv
gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
ENSG00000000001,1672.45,2.34,0.23,10.17,2.31e-24,1.45e-20
ENSG00000000002,213.67,-1.87,0.31,-6.03,1.63e-09,4.21e-07
```

**edgeR:**
```csv
gene_id,logFC,logCPM,F,PValue,FDR
ENSG00000000001,2.34,10.71,103.49,2.31e-24,1.45e-20
ENSG00000000002,-1.87,7.74,36.41,1.63e-09,4.21e-07
```

**limma:**
```csv
gene_id,logFC,AveExpr,t,P.Value,adj.P.Val,B
ENSG00000000001,2.34,10.71,10.17,2.31e-24,1.45e-20,40.23
ENSG00000000002,-1.87,7.74,-6.03,1.63e-09,4.21e-07,14.56
```

### Profile JSON Format

```json
{
  "timestamp": "2026-03-09T10:30:00",
  "raptor_version": "2.2.0",
  "module": "M3",
  "features": {
    "n_samples": 6,
    "n_genes": 20000,
    "min_group_size": 3,
    "library_size_mean": 4500000,
    "library_size_cv": 0.08,
    "bcv": 0.32,
    "bcv_category": "moderate",
    "zero_fraction": 0.15,
    "has_outliers": false,
    "has_batch_effect": false
  }
}
```

---

## Troubleshooting Guide

### Common Issues Across All Modules

#### Issue: "ModuleNotFoundError: No module named 'raptor'"

**Solution:**
```bash
# Install RAPTOR
cd RAPTOR/
pip install -e .

# Verify installation
python -c "import raptor; print(raptor.__version__)"
```

#### Issue: "FileNotFoundError: Count file not found"

**Solution:**
```bash
# Check current directory
pwd

# Check file exists
ls -l results/quick_counts/quick_gene_counts.csv

# Use absolute path
python 02_quality_assessment.py \
  --counts /full/path/to/quick_gene_counts.csv
```

#### Issue: "MemoryError" or "Killed"

**Solution:**
```bash
# Use smaller dataset for testing
# Or increase available RAM
# Or use chunked processing

# For M5, use lower-memory pipeline
python 05_production_pipeline.py \
  --pipeline kallisto  # Only 4GB needed
```

### Module-Specific Issues

#### M1: Quantification Issues

**Low mapping rate (<50%):**
- Check organism matches reference
- Verify FASTQ quality
- Try different index

**Slow processing:**
- Increase threads: `--threads 16`
- Use Kallisto instead of Salmon
- Check system resources

#### M2: Quality Issues

**Low quality scores (<70):**
- Review specific issues in qc_results.json
- May need more sequencing depth
- Check for technical problems

**Too many outliers:**
- Relax threshold: `--outlier-threshold 3.5`
- Review individual samples
- May indicate biological variation

#### M3: Profiling Issues

**Missing features:**
- Re-run with metadata
- Check count matrix format
- Ensure adequate sample size

#### M4: Recommendation Issues

**No clear recommendation:**
- Normal for edge cases
- Run top 2 methods in R
- Use ensemble (M9) to combine

#### M5: Pipeline Issues

**Tool not found:**
```bash
conda install -c bioconda salmon star kallisto
```

**Out of memory:**
- Use lower-memory pipeline
- Reduce threads
- Process samples individually

#### M7: Import Issues

**Column names not recognized:**
- Specify pipeline: `--pipeline deseq2`
- Check column names match expected
- Ensure gene IDs present

#### M8: Optimization Issues

**Takes too long:**
- Reduce grid points
- Use random search
- Use fdr_control method

#### M9: Ensemble Issues

**Methods disagree:**
- Normal biological variation
- Use voting with min_methods=3
- Check individual method quality

---

## Tips & Best Practices

### General Workflow

1. **Always start with demo mode** to learn
2. **Use fast profiling (M1-M4)** before full processing
3. **Run QC (M2)** to catch issues early
4. **Follow recommendations (M4)** - they're data-driven
5. **Use ensemble (M9)** for high-confidence results

### Data Organization

```
project/
├── data/
│   ├── fastqs/              # Raw FASTQ files
│   ├── references/          # Salmon/STAR indices
│   └── metadata/            # Sample information
├── results/
│   ├── quick_counts/        # M1 output
│   ├── qc/                  # M2 output
│   ├── profile.json         # M3 output
│   ├── recommendation/      # M4 output
│   ├── production_pipeline/ # M5 output
│   ├── de_analysis/         # M6 output (R)
│   ├── de_imported/         # M7 output
│   ├── parameter_optimization/ # M8 output
│   └── ensemble_analysis/   # M9 output
└── RAPTOR/
    └── examples/            # Example scripts
```

### Performance Optimization

**For Speed:**
- Use Kallisto (M1, M5)
- Increase threads
- Use SSD storage
- Process samples in parallel

**For Memory:**
- Use Kallisto instead of STAR
- Process samples sequentially
- Use smaller reference

**For Accuracy:**
- Use Salmon with bootstraps (M1, M5)
- Use STAR + Salmon for BAM + accuracy (M5)
- Run multiple DE methods → ensemble (M6-M9)

### Quality Control Checkpoints

**After M1:**
- Check mapping rates >70%
- Verify library sizes are reasonable
- Check for failed samples

**After M2:**
- Overall quality score >70
- No excessive outliers
- Batch effects identified

**After M3:**
- Profile has all 32 features
- Sample size adequate for DE
- BCV in reasonable range

**After M7:**
- Import successful for all methods
- n_significant > 0
- Direction consistency checked

**After M9:**
- Methods show reasonable agreement
- Consensus genes make biological sense
- Ready for validation

---

## FAQ

### General Questions

**Q: Do I need to run all modules?**

A: No! Common workflows:
- Quick QC: M1 → M2 → M3
- Full analysis: M1-M4 → M5 → M6 → M7 → M9
- With existing counts: M2 → M3 → M4
- With existing DE results: M7 → M8 → M9

**Q: Can I use my existing count matrix?**

A: Yes! Start at M2:
```bash
python 02_quality_assessment.py --counts my_counts.csv
```

**Q: Can I use my existing DE results?**

A: Yes! Start at M7:
```bash
python 07_DE_Import.py --de-file my_deseq2_results.csv
```

**Q: How long does the complete workflow take?**

A:
- M1-M4 (Fast profiling): 5-15 minutes
- M5 (Production): 30-120 min/sample
- M6 (R DE): 5-30 minutes
- M7-M9 (Import & Ensemble): 1-5 minutes

### Technical Questions

**Q: What's the difference between M1 and M5?**

A: 
- **M1**: Fast quantification for QC and profiling (5-15 min)
- **M5**: Production quantification for final analysis (30-120 min)

Use M1 for initial assessment, M5 for publication.

**Q: Should I use DESeq2, edgeR, or limma?**

A: Run M4 to get data-driven recommendation! Generally:
- Small samples (n<3): DESeq2
- Moderate BCV: edgeR
- Large samples (n≥8): limma

**Q: What if methods disagree in ensemble (M9)?**

A: Normal! Options:
1. Use voting with min_methods=3 (strictest)
2. Check individual method quality
3. Use biological knowledge to adjudicate
4. Consider all methods may be correct (biological complexity)

**Q: Do I need validation data for optimization (M8)?**

A: No! Use `fdr_control` method - works for everyone.

**Q: How many DE methods should I run?**

A: For ensemble (M9):
- Minimum: 2 methods
- Recommended: 3 methods (DESeq2 + edgeR + limma)
- Ideal: All applicable methods

### Practical Questions

**Q: My dataset has only 2 samples per group. Can I still use RAPTOR?**

A: Yes, but with limitations:
- M1-M4 work fine
- M4 will recommend DESeq2 (best for small n)
- Statistical power will be limited
- Consider validating top genes

**Q: Can I analyze single-cell RNA-seq data?**

A: RAPTOR is designed for bulk RNA-seq. For scRNA-seq, use specialized tools (Seurat, Scanpy).

**Q: What if I have batch effects?**

A: 
- M2 detects them
- M4 recommends appropriate methods (DESeq2/limma handle batches)
- Include batch in design matrix in R

**Q: Can I use RAPTOR for non-model organisms?**

A: Yes! You need:
- Transcriptome FASTA for Salmon/Kallisto
- Or genome + annotation for STAR
- Module 1-9 work independently of organism

---

## Getting Help

### Documentation

- **Full Documentation**: `docs/` directory
- **API Reference**: `docs/api/`
- **Tutorials**: `docs/tutorials/`

### Example Data

```bash
# Download example dataset
wget https://raptor.example.com/data/demo_dataset.tar.gz
tar -xzf demo_dataset.tar.gz
```

### Community

- **GitHub Issues**: https://github.com/ayehbolouki/RAPTOR/issues
- **Discussions**: https://github.com/ayehbolouki/RAPTOR/discussions
- **Email**: ayehbolouki1988@gmail.com

### Citation

If you use RAPTOR in your research, please cite:

```
Bolouki, A. (2026). RAPTOR v2.2.0: RNA-seq Analysis Pipeline Tool 
for Omics Research. GitHub. https://github.com/ayehbolouki/RAPTOR
```

---

## Quick Reference Card

### Complete Workflow Commands

```bash
# Stage 1: Fast Profiling (5-15 min)
python 01_quick_count.py --fastq-dir data/fastqs/ --index data/salmon_index/
python 02_quality_assessment.py --counts results/quick_counts/quick_gene_counts.csv
python 03_data_profiler.py --counts results/qc/counts_clean.csv
python 04_recommender.py --profile results/profile.json

# Stage 2: Production Quantification (30-120 min)
python 05_production_pipeline.py --fastq-dir data/fastqs/ --pipeline salmon

# Stage 3: DE Analysis
# [Run R analysis - M6]
python 07_DE_Import.py --de-file results/deseq2_results.csv --pipeline auto
python 08_Parameter_Optimization.py --de-result results/de_imported/de_result.pkl --method fdr_control
python 09_Ensemble_Analysis.py --de-results results/de_imported/*.pkl --method all
```

### Demo Mode Commands

```bash
python 01_quick_count.py --demo
python 02_quality_assessment.py --demo
python 03_data_profiler.py --demo
python 04_recommender.py --demo
python 07_DE_Import.py --demo
python 08_Parameter_Optimization.py --demo
python 09_Ensemble_Analysis.py --demo
```

### Common Options

```
--demo              Run in demo mode (all modules)
--output DIR        Custom output directory
--threads N         Number of CPU threads
--help              Show help message
```

---

**Making free science for everybody around the world 🌍**

*Version: 2.2.0*  
*Last Updated: March 2026*  
*Author: Ayeh Bolouki*
