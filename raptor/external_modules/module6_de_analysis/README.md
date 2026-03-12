# 🧬 RAPTOR Module 6 - Differential Expression Analysis (External Module)

**Version:** 2.2.0  
**Type:** External R Module  
**Purpose:** Statistical differential expression analysis using DESeq2, edgeR, or limma-voom

---

## 📋 TABLE OF CONTENTS

1. [Overview](#overview)
2. [Module 6 Workflow](#module-6-workflow)
3. [Installation](#installation)
4. [Running DE Analysis](#running-de-analysis)
5. [Importing Results to RAPTOR](#importing-results-to-raptor)
6. [Output Files](#output-files)
7. [Method Comparison](#method-comparison)
8. [Next Steps](#next-steps)
9. [Troubleshooting](#troubleshooting)

---

## 🎯 OVERVIEW

### What is Module 6?

Module 6 performs **statistical differential expression (DE) analysis** using industry-standard R/Bioconductor packages. This module is **external to RAPTOR's Python pipeline** - you run the R scripts independently, then import results back into RAPTOR for downstream analysis.

### Why External?

- ✅ Uses mature, validated R/Bioconductor tools (DESeq2, edgeR, limma)
- ✅ Full control over statistical parameters
- ✅ Access to all package features and updates
- ✅ Familiar workflow for bioinformaticians
- ✅ Results are standardized for RAPTOR import

### RAPTOR Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    RAPTOR WORKFLOW                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  M1-M5: Python Pipelines (Quantification → Recommendation)     │
│         ↓                                                       │
│  M6:    External R Analysis (YOU RUN THIS OUTSIDE RAPTOR) ← 🔴│
│         ↓                                                       │
│  M7:    Import DE Results (BACK TO RAPTOR)                     │
│         ↓                                                       │
│  M8-M10: Python Pipelines (Optimization → Biomarkers)          │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🔄 MODULE 6 WORKFLOW

### Step 1: Complete RAPTOR Modules 1-5

Run RAPTOR's quantification and recommendation pipeline:

```bash
# Module 1: Quick Quantification
raptor quick-count -m salmon -s samples.csv -i index/ -o results/quick_counts

# Module 2: Quality Assessment
raptor qc --counts results/quick_counts/quick_gene_counts.csv

# Module 3: Data Profiling
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# Module 4: Pipeline Recommendation
raptor recommend
```

**Output:** `recommendation.yaml` (contains suggested DE parameters)

### Step 2: Run External DE Analysis (Module 6)

Use the R scripts provided in this module to perform statistical analysis:

```bash
# Choose one method: DESeq2, edgeR, or limma-voom

# DESeq2 (recommended for most cases)
Rscript run_deseq2.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_analysis \
    --config recommendation.yaml \
    --plots

# OR edgeR (for smaller samples)
Rscript run_edger.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_analysis \
    --config recommendation.yaml

# OR limma-voom (for versatility)
Rscript run_limma.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_analysis \
    --config recommendation.yaml
```

**Output:** 
- `de_results.csv` (standardized for RAPTOR)
- `de_significant.csv` (significant genes only)
- `de_summary.json` (metadata)

### Step 3: Import Results to RAPTOR (Module 7)

Import the DE results back into RAPTOR for downstream analysis:

```bash
# Import DE results
raptor import-de \
    --de-file results/de_analysis/de_results.csv \
    --summary results/de_analysis/de_summary.json
```

**Output:** DE results integrated into RAPTOR database

### Step 4: Continue with Modules 8-10

Proceed with RAPTOR's downstream analysis:

```bash
# Module 8: Parameter Optimization
raptor optimize-params

# Module 9: Ensemble Analysis
raptor ensemble-analysis

# Module 10: Biomarker Discovery
raptor discover-biomarkers
```

---

## 📦 INSTALLATION

### Prerequisites

- **R** ≥ 4.0.0
- **Python** ≥ 3.8 (for RAPTOR)
- **RAPTOR** v2.2.0 installed

### Install R Packages

Run the installation script to install all required Bioconductor packages:

```bash
# Install all required packages
Rscript install_packages.R

# Test installation
Rscript install_packages.R --test
```

**Required Packages:**
- DESeq2 (differential expression)
- edgeR (differential expression)
- limma (differential expression)
- apeglm (LFC shrinkage for DESeq2)
- ashr (LFC shrinkage for DESeq2)
- optparse (command-line parsing)
- jsonlite (JSON output)
- yaml (config files)
- ggplot2 (plotting)

**Installation time:** ~10-15 minutes

---

## 🚀 RUNNING DE ANALYSIS

### Input Requirements

**1. Count Matrix** (CSV or TSV):
```
gene_id,Sample1,Sample2,Sample3,Sample4
ENSG00000000003,1234,2345,3456,4567
ENSG00000000005,567,678,789,890
...
```

**2. Metadata** (CSV):
```
sample_id,condition,batch
Sample1,Control,1
Sample2,Control,1
Sample3,Treatment,2
Sample4,Treatment,2
```

**3. RAPTOR Recommendation** (optional YAML):
```yaml
fdr_threshold: 0.05
lfc_threshold: 0.0
min_count: 10
normalization: TMM
recommended_method: DESeq2
```

### Method 1: DESeq2 (Recommended)

**Best for:** Most RNA-seq experiments, especially with biological replicates

```bash
Rscript run_deseq2.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_deseq2 \
    --condition treatment \
    --reference Control \
    --fdr 0.05 \
    --lfc 0 \
    --shrinkage apeglm \
    --plots
```

**Key Parameters:**
- `--shrinkage apeglm` - LFC shrinkage (apeglm, ashr, normal, none)
- `--fit-type parametric` - Dispersion fit (parametric, local, mean)
- `--batch batch_column` - Batch correction
- `--threads 8` - Parallel processing

### Method 2: edgeR

**Best for:** Small sample sizes, robust performance

```bash
Rscript run_edger.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_edger \
    --condition treatment \
    --reference Control \
    --fdr 0.05 \
    --normalization TMM \
    --method QLF \
    --robust
```

**Key Parameters:**
- `--normalization TMM` - Normalization (TMM, RLE, upperquartile)
- `--method QLF` - Test method (QLF, LRT, exact)
- `--robust` - Robust dispersion estimation

### Method 3: limma-voom

**Best for:** Versatile, handles complex designs well

```bash
Rscript run_limma.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_limma \
    --condition treatment \
    --reference Control \
    --fdr 0.05 \
    --voom-method voom \
    --robust
```

**Key Parameters:**
- `--voom-method voom` - Transformation (voom, voomWithQualityWeights)
- `--robust` - Robust empirical Bayes
- `--trend` - Mean-variance trend

### Using RAPTOR Recommendations

If you have a `recommendation.yaml` from Module 4:

```bash
Rscript run_deseq2.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_analysis \
    --config recommendation.yaml \
    --plots
```

RAPTOR will automatically use recommended parameters (FDR, LFC thresholds, etc.)

---

## 📥 IMPORTING RESULTS TO RAPTOR

After running DE analysis, import results back into RAPTOR:

### Import Command

```bash
raptor import-de \
    --de-file results/de_analysis/de_results.csv \
    --summary results/de_analysis/de_summary.json \
    --method deseq2
```

### What Gets Imported

The standardized DE results include:
- Gene IDs
- Log2 fold changes
- P-values and adjusted p-values
- Base mean expression
- Significance flags
- Direction (up/down)

### Validation

RAPTOR validates the imported data:
- ✅ Required columns present
- ✅ No missing critical values
- ✅ FDR values in valid range
- ✅ Gene IDs match reference

### Import Confirmation

```
📥 Importing DE results...
   Method: DESeq2
   Genes tested: 15,234
   Significant: 1,456 (9.6%)
   Upregulated: 823
   Downregulated: 633

✅ Import successful!
   Results saved to RAPTOR database
   Ready for Module 8: Parameter Optimization
```

---

## 📁 OUTPUT FILES

Each DE analysis produces standardized output:

### 1. de_results.csv (Full Results)

All genes tested with complete statistics:

```csv
gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,direction,is_significant
ENSG00000000003,1234.5,2.34,0.23,10.2,1.2e-24,3.4e-22,up,TRUE
ENSG00000000005,567.8,-1.45,0.34,-4.3,2.3e-5,0.003,down,TRUE
ENSG00000000419,8901.2,0.12,0.15,0.8,0.42,0.89,unchanged,FALSE
```

**Use for:** RAPTOR import (Module 7)

### 2. de_significant.csv (Significant Only)

Filtered to significant genes only:

```csv
gene_id,baseMean,log2FoldChange,padj,direction
ENSG00000000003,1234.5,2.34,3.4e-22,up
ENSG00000000005,567.8,-1.45,0.003,down
```

**Use for:** Quick review, gene lists, downstream analysis

### 3. de_summary.json (Metadata)

Analysis metadata for reproducibility:

```json
{
  "timestamp": "2026-01-19T14:30:00",
  "raptor_version": "2.2.0",
  "module": "M6",
  "pipeline": "DESeq2",
  "parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 0.0,
    "shrinkage": "apeglm"
  },
  "results": {
    "n_significant": 1456,
    "n_upregulated": 823,
    "n_downregulated": 633
  }
}
```

**Use for:** RAPTOR import validation, provenance tracking

### 4. de_plots.pdf (Optional QC)

Quality control plots:
- MA plot (log2FC vs mean expression)
- Dispersion estimates
- PCA (sample clustering)
- Volcano plot (significance vs fold change)

**Use for:** Quality assessment, publication figures

---

## 📊 METHOD COMPARISON

### When to Use Each Method

| Method | Best For | Sample Size | Speed | Flexibility |
|--------|----------|-------------|-------|-------------|
| **DESeq2** | General use, biological replicates | Any (≥3 per group) | Medium | High |
| **edgeR** | Small samples, TMM normalization | ≥2 per group | Fast | High |
| **limma-voom** | Complex designs, paired samples | ≥3 per group | Fastest | Highest |

### Recommendations

**Use DESeq2 if:**
- ✅ Standard RNA-seq experiment
- ✅ Biological replicates available
- ✅ Want LFC shrinkage
- ✅ Need robust analysis

**Use edgeR if:**
- ✅ Small sample size (2-3 per group)
- ✅ Need fast analysis
- ✅ Simple design
- ✅ Want exact tests

**Use limma-voom if:**
- ✅ Complex experimental design
- ✅ Multiple factors/covariates
- ✅ Paired samples
- ✅ Need speed and flexibility

**Not sure?** Use DESeq2 - it's the most widely used and robust.

### Method Concordance

In most cases, all three methods will identify similar significant genes:
- ~80-90% overlap for highly significant genes (padj < 0.01)
- ~70-80% overlap for moderately significant genes (padj < 0.05)
- Differences mainly in borderline cases

---

## 🔜 NEXT STEPS

After importing DE results (Module 7), continue with RAPTOR's downstream modules:

### Module 8: Parameter Optimization

Optimize statistical thresholds using your data:

```bash
raptor optimize-params \
    --de-results results/de_analysis/de_results.csv \
    --optimize fdr,lfc \
    --criterion reproducibility
```

**Purpose:**
- Find optimal FDR threshold
- Determine appropriate LFC cutoff
- Balance sensitivity vs specificity
- Maximize reproducibility

### Module 9: Ensemble Analysis

Combine multiple DE methods for robust gene lists:

```bash
raptor ensemble-analysis \
    --deseq2 results/de_deseq2/de_results.csv \
    --edger results/de_edger/de_results.csv \
    --limma results/de_limma/de_results.csv \
    --strategy consensus
```

**Purpose:**
- Combine results from multiple methods
- Identify consensus genes
- Reduce false positives
- Increase confidence

### Module 10: Biomarker Discovery

Identify potential biomarkers from DE genes:

```bash
raptor discover-biomarkers \
    --de-results results/de_analysis/de_results.csv \
    --strategy multi-criteria \
    --rank-by significance,effect-size,consistency
```

**Purpose:**
- Rank genes by biomarker potential
- Multi-criteria scoring
- Pathway enrichment
- Generate biomarker panels

---

## 🔧 TROUBLESHOOTING

### Installation Issues

**Problem:** BiocManager not found
```bash
# Solution: Install BiocManager first
R -e "install.packages('BiocManager')"
```

**Problem:** Package installation fails
```bash
# Solution: Update R and Bioconductor
R -e "BiocManager::install(version = '3.18', update = TRUE, ask = FALSE)"
```

### Analysis Issues

**Problem:** "Condition column not found"
```bash
# Solution: Check metadata column names
# Use --condition to specify correct column
Rscript run_deseq2.R ... --condition treatment_group
```

**Problem:** "Samples don't match"
```bash
# Solution: Ensure sample names match between counts and metadata
# Check for extra spaces, case sensitivity
```

**Problem:** "Not enough replicates"
```bash
# Solution: Need ≥2 replicates per condition for DESeq2/edgeR
# Use limma-voom for n=1 (not recommended)
```

**Problem:** "All p-values are NA"
```bash
# Solution: Likely due to low counts or filtering
# Try: --min-count 5 (instead of 10)
# Check: Are there enough reads?
```

### Import Issues

**Problem:** "de_results.csv not found"
```bash
# Solution: Verify output directory
# Check: Did R script complete successfully?
ls -lh results/de_analysis/
```

**Problem:** "Invalid format"
```bash
# Solution: Ensure you're using output from our R scripts
# They produce standardized format for RAPTOR
```

### R Script Errors

**Problem:** "Rscript: command not found"
```bash
# Solution: Add R to PATH or use full path
/usr/local/bin/Rscript run_deseq2.R ...
```

**Problem:** "Package X is not available"
```bash
# Solution: Install missing package
R -e "BiocManager::install('X')"
```

---

## 📚 ADDITIONAL RESOURCES

### Documentation
- **DESeq2:** https://bioconductor.org/packages/DESeq2/
- **edgeR:** https://bioconductor.org/packages/edgeR/
- **limma:** https://bioconductor.org/packages/limma/

### Tutorials
- DESeq2 vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/
- edgeR user guide: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/
- limma user guide: https://bioconductor.org/packages/release/bioc/vignettes/limma/

### Citation
If you use this module, please cite:
- **RAPTOR:** Bolouki, A. (2026). RAPTOR v2.2.0
- **DESeq2:** Love, M.I., Huber, W., Anders, S. (2014). Genome Biology
- **edgeR:** Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). Bioinformatics
- **limma:** Ritchie, M.E., et al. (2015). Nucleic Acids Research

---

## 📋 QUICK REFERENCE

### Complete Workflow

```bash
# 1. Run RAPTOR Modules 1-5
raptor quick-count -m salmon -s samples.csv -i index/ -o results/
raptor qc --counts results/quick_counts/quick_gene_counts.csv
raptor profile --counts results/quick_counts/quick_gene_counts.csv
raptor recommend

# 2. Install R packages (first time only)
Rscript install_packages.R

# 3. Run DE analysis (Module 6)
Rscript run_deseq2.R \
    --counts results/gene_counts.csv \
    --metadata data/metadata.csv \
    --output results/de_analysis \
    --config recommendation.yaml \
    --plots

# 4. Import results (Module 7)
raptor import-de \
    --de-file results/de_analysis/de_results.csv \
    --summary results/de_analysis/de_summary.json

# 5. Continue with Modules 8-10
raptor optimize-params
raptor ensemble-analysis
raptor discover-biomarkers
```

### Common Commands

```bash
# DESeq2 with defaults
Rscript run_deseq2.R -c counts.csv -m metadata.csv -o results/

# edgeR with TMM normalization
Rscript run_edger.R -c counts.csv -m metadata.csv --normalization TMM

# limma with quality weights
Rscript run_limma.R -c counts.csv -m metadata.csv --voom-method voomWithQualityWeights

# Check installed packages
Rscript install_packages.R --test

# Get help
Rscript run_deseq2.R --help
```

---

## 📧 SUPPORT

**Author:** Ayeh Bolouki  
**Email:** ayehbolouki1988@gmail.com  
**Version:** 2.2.0  
**License:** MIT

---

**Making free science for everybody around the world** 🌍
