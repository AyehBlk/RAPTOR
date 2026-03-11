# RAPTOR Module 7: DE Import - User Guide

**A Comprehensive Guide to Importing and Standardizing Differential Expression Results**

**Version**: 2.2.0  
**Last Updated**: March 9, 2026  
**Difficulty**: Beginner to Intermediate  
**Estimated Time**: 5-10 minutes per dataset

---

## 📚 Table of Contents

1. [Introduction](#introduction)
2. [Why Import DE Results?](#why-import-de-results)
3. [Supported DE Analysis Tools](#supported-de-analysis-tools)
4. [Getting Started](#getting-started)
5. [Input File Formats](#input-file-formats)
6. [Auto-Detection Feature](#auto-detection-feature)
7. [Step-by-Step Tutorials](#step-by-step-tutorials)
8. [Understanding the Output](#understanding-the-output)
9. [Working with DEResult Objects](#working-with-deresult-objects)
10. [Integration with Your Workflow](#integration-with-your-workflow)
11. [Best Practices](#best-practices)
12. [Troubleshooting](#troubleshooting)
13. [Frequently Asked Questions](#frequently-asked-questions)
14. [Real-World Examples](#real-world-examples)
15. [Technical Details](#technical-details)

---

## Introduction

### What is Module 7?

Module 7 (DE Import) is RAPTOR's **differential expression results importer and standardizer**. It takes DE results from popular R packages (DESeq2, edgeR, limma-voom, Wilcoxon) and converts them into a unified format that works seamlessly with RAPTOR's downstream modules.

**The Problem**: Different DE analysis tools produce results in different formats with different column names, making it difficult to:
- Compare results across methods
- Use multiple tools in ensemble analysis
- Apply downstream analysis consistently

**The Solution**: Module 7 automatically detects the analysis method, standardizes column names, calculates significance metrics, and prepares data for Modules 8-10.

### What Does Module 7 Do?

**In simple terms**: Module 7 is a "translator" that speaks 4 different DE analysis languages and converts them all to a common format.

```
Input:  DESeq2 results with "baseMean", "log2FoldChange", "padj"
        ↓
Module 7: Standardization + Validation
        ↓
Output: Unified format with "log2_fold_change", "adjusted_p_value", "is_significant"
```

### Key Features

✅ **4 Supported Tools**: DESeq2, edgeR, limma-voom, Wilcoxon  
✅ **Auto-Detection**: Automatically identifies which tool was used  
✅ **Smart Standardization**: Maps different column names to standard format  
✅ **Significance Calculation**: Adds `is_significant` and `direction` columns  
✅ **Quality Validation**: Checks for missing values, invalid ranges  
✅ **Multiple Outputs**: CSV files, JSON summary, pickle object  
✅ **Ready for Downstream**: Prepares data for Modules 8-10  

### Where Does Module 7 Fit?

```
Your RNA-seq Analysis Pipeline:

Step 1: Quality Control (Module 2)
Step 2: Data Profiling (Module 3)
Step 3: Pipeline Selection (Module 4)
         ↓
Step 4: Run DE Analysis in R
        - DESeq2, edgeR, limma-voom, or Wilcoxon
        - Export results to CSV
         ↓
Step 5: IMPORT DE RESULTS (Module 7) ⭐ YOU ARE HERE
         ↓
        Creates unified format for:
         ↓
Step 6: Parameter Optimization (Module 8)
Step 7: Ensemble Analysis (Module 9)
Step 8: Biomarker Discovery (Module 10)
```

---

## Why Import DE Results?

### The Column Name Problem

**Different tools, different names**:

| What it means | DESeq2 | edgeR | limma-voom | Wilcoxon |
|---------------|--------|-------|------------|----------|
| Fold change | `log2FoldChange` | `logFC` | `logFC` | `log2FC` |
| P-value | `pvalue` | `PValue` | `P.Value` | `pvalue` |
| Adjusted p-value | `padj` | `FDR` | `adj.P.Val` | `padj` |
| Mean expression | `baseMean` | `logCPM` | `AveExpr` | N/A |

**Without Module 7**: You'd need custom code for each tool  
**With Module 7**: One command works for all!

---

### The Standardization Solution

**Before** (inconsistent):
```csv
# DESeq2 format
gene_id,baseMean,log2FoldChange,pvalue,padj
GENE1,500.2,2.5,0.001,0.01

# edgeR format  
gene,logCPM,logFC,PValue,FDR
GENE1,5.2,2.5,0.001,0.01
```

**After** (standardized):
```csv
# Unified RAPTOR format (works for all tools!)
gene_id,log2_fold_change,p_value,adjusted_p_value,base_mean,is_significant,direction
GENE1,2.5,0.001,0.01,500.2,True,up
```

---

### Benefits of Importing

1. ✅ **Consistency**: Same format regardless of DE tool used
2. ✅ **Validation**: Automatic quality checks
3. ✅ **Enrichment**: Adds significance flags and direction
4. ✅ **Integration**: Works seamlessly with Modules 8-10
5. ✅ **Comparison**: Easy to compare different DE methods
6. ✅ **Documentation**: Tracks analysis metadata
7. ✅ **Reproducibility**: Standardized, versioned format

---

## Supported DE Analysis Tools

Module 7 supports the **4 most popular differential expression tools**:

### 1. DESeq2

**What it is**: The most widely used DE tool for RNA-seq

**Based on**: Negative binomial distribution  
**Best for**: General purpose RNA-seq analysis  
**Paper**: Love et al. (2014) Genome Biology

**Typical R workflow**:
```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData, ~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "deseq2_results.csv")
```

**RAPTOR import**:
```bash
raptor import-de -i deseq2_results.csv -m deseq2
```

---

### 2. edgeR

**What it is**: Popular tool from Robinson lab, uses empirical Bayes

**Based on**: Negative binomial GLM  
**Best for**: Small sample sizes, high biological variability  
**Paper**: Robinson et al. (2010) Bioinformatics

**Typical R workflow**:
```r
library(edgeR)
y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
res <- topTags(qlf, n=Inf)
write.csv(res$table, "edger_results.csv")
```

**RAPTOR import**:
```bash
raptor import-de -i edger_results.csv -m edger
```

---

### 3. limma-voom

**What it is**: Limma adapted for RNA-seq using voom transformation

**Based on**: Linear modeling with precision weights  
**Best for**: Large sample sizes (n>20), complex designs  
**Paper**: Law et al. (2014) Genome Biology

**Typical R workflow**:
```r
library(limma)
library(edgeR)
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~condition)
v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, n=Inf)
write.csv(res, "limma_results.csv")
```

**RAPTOR import**:
```bash
raptor import-de -i limma_results.csv -m limma
```

---

### 4. Wilcoxon Rank-Sum Test

**What it is**: Non-parametric test, no distribution assumptions

**Based on**: Rank statistics  
**Best for**: Non-normal data, n≥8 per group  
**Paper**: Li et al. (2022) Briefings in Bioinformatics

**Typical R workflow**:
```r
# Manual implementation or from tools
results <- data.frame(
    gene_id = rownames(counts),
    log2FC = log2FoldChanges,
    pvalue = wilcox_pvalues,
    padj = p.adjust(wilcox_pvalues, method="BH")
)
write.csv(results, "wilcoxon_results.csv")
```

**RAPTOR import**:
```bash
raptor import-de -i wilcoxon_results.csv -m wilcoxon
```

---

### Comparison of Methods

| Tool | Distribution | Sample Size | Complexity | Speed |
|------|-------------|-------------|------------|-------|
| **DESeq2** | Neg. Binomial | n≥3 | Medium | Medium |
| **edgeR** | Neg. Binomial | n≥3 | Medium | Fast |
| **limma-voom** | Normal (transformed) | n>20 | Low | Fastest |
| **Wilcoxon** | None (non-parametric) | n≥8 | Low | Fast |

**Module 7 works with all 4!**

---

## Getting Started

### Prerequisites

1. ✅ **DE analysis completed in R**
   - You've run DESeq2, edgeR, limma, or Wilcoxon
   - Results exported to CSV/TSV file

2. ✅ **Results file saved**
   - CSV or TSV format
   - Contains required columns (fold change, p-values)
   - Gene IDs in first column or in column named 'gene_id'

3. ✅ **RAPTOR installed**
   ```bash
   pip install raptor-rnaseq
   ```

### Quick Start (3 Steps)

**Step 1**: Check your DE results file
```bash
head deseq2_results.csv
```

**Step 2**: Import with Module 7
```bash
raptor import-de -i deseq2_results.csv -m deseq2
```

**Step 3**: Check the output
```bash
ls results/de_imported/
# You should see:
# - de_standardized.csv
# - de_significant.csv
# - de_summary.json
# - de_result.pkl
```

That's it! Your results are now in standardized RAPTOR format.

### Basic Command Structure

```bash
raptor import-de \
    -i INPUT_FILE.csv \          # Your DE results
    -m METHOD \                  # deseq2, edger, limma, or wilcoxon
    -o OUTPUT_DIR \              # Where to save (optional)
    --fdr-threshold 0.05 \       # Significance cutoff (optional)
    --lfc-threshold 0.0          # Fold-change cutoff (optional)
```

### Expected Output

After running, you'll see:
```
🦖 RAPTOR v2.2.0 - Import DE Results (Module 7)

📂 Loading DE results from: deseq2_results.csv
   Loaded 15,000 genes

   Auto-detected pipeline: DESeq2

🔄 Standardizing column names...
   ✓ Mapped: log2FoldChange → log2_fold_change
   ✓ Mapped: pvalue → p_value
   ✓ Mapped: padj → adjusted_p_value
   ✓ Mapped: baseMean → base_mean

📊 Calculating significance (FDR=0.05, LFC=0.0)...
   ✓ 1,850 significant genes (12.3%)

📁 Saving results to: results/de_imported/
   ✓ de_standardized.csv (15,000 genes)
   ✓ de_significant.csv (1,850 genes)
   ✓ de_summary.json
   ✓ de_result.pkl (for M8-M10)

✅ Import complete!
```

---

## Input File Formats

### General Requirements

**All input files must have**:

1. ✅ **Gene identifiers**
   - First column OR column named 'gene_id'
   - Can be: Ensembl IDs, gene symbols, transcript IDs

2. ✅ **Fold change column**
   - log2 fold change values
   - Can be named: log2FoldChange, logFC, log2FC, etc.

3. ✅ **P-value columns**
   - Raw p-value AND adjusted p-value (FDR)
   - Various names accepted (see below)

4. ✅ **File format**
   - CSV (comma-separated) or TSV (tab-separated)
   - Text file with .csv, .tsv, or .txt extension

---

### DESeq2 Format

**Required columns** (at least one name from each group):

| Data | Accepted Column Names |
|------|----------------------|
| **Gene ID** | `gene_id`, `gene`, `Gene` (or first column) |
| **Fold Change** | `log2FoldChange`, `log2FC`, `lfc`, `logFC` |
| **P-value** | `pvalue`, `pval`, `PValue`, `p.value` |
| **Adjusted P** | `padj`, `FDR`, `adj.P.Val`, `qvalue`, `q.value` |
| **Mean Expression** | `baseMean`, `AveExpr`, `logCPM`, `meanExpr` (optional) |
| **Statistic** | `stat`, `t`, `LR`, `F`, `z` (optional) |

**Example DESeq2 output**:
```csv
gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
ENSG00000141510,5234.5,2.45,0.234,10.47,1.23e-25,4.56e-22
ENSG00000012048,3456.2,-1.87,0.198,-9.44,3.89e-21,7.12e-18
ENSG00000139618,8901.3,1.92,0.187,10.27,8.92e-25,3.21e-21
```

**Minimal DESeq2 format** (will work):
```csv
gene_id,log2FoldChange,pvalue,padj
GENE1,2.45,1.23e-25,4.56e-22
GENE2,-1.87,3.89e-21,7.12e-18
```

---

### edgeR Format

**Required columns**:

| Data | Accepted Column Names |
|------|----------------------|
| **Gene ID** | `gene`, `gene_id`, `Gene` (or first column) |
| **Fold Change** | `logFC`, `log2FC` |
| **P-value** | `PValue`, `pvalue` |
| **Adjusted P** | `FDR`, `adj.P.Val` |
| **Expression** | `logCPM`, `AveExpr` (optional) |

**Example edgeR output**:
```csv
gene,logFC,logCPM,LR,PValue,FDR
GENE1,2.45,5.234,109.5,1.23e-25,4.56e-22
GENE2,-1.87,4.156,89.1,3.89e-21,7.12e-18
GENE3,1.92,6.789,105.4,8.92e-25,3.21e-21
```

**Minimal edgeR format**:
```csv
gene,logFC,PValue,FDR
GENE1,2.45,1.23e-25,4.56e-22
GENE2,-1.87,3.89e-21,7.12e-18
```

---

### limma-voom Format

**Required columns**:

| Data | Accepted Column Names |
|------|----------------------|
| **Gene ID** | `gene`, `gene_id`, `Gene` (or first column) |
| **Fold Change** | `logFC`, `log2FC` |
| **P-value** | `P.Value`, `pvalue` |
| **Adjusted P** | `adj.P.Val`, `FDR` |
| **Expression** | `AveExpr`, `logCPM` (optional) |
| **Statistic** | `t`, `F` (optional) |

**Example limma output**:
```csv
gene,logFC,AveExpr,t,P.Value,adj.P.Val,B
GENE1,2.45,8.234,10.47,1.23e-25,4.56e-22,48.5
GENE2,-1.87,7.156,-9.44,3.89e-21,7.12e-18,38.9
GENE3,1.92,9.789,10.27,8.92e-25,3.21e-21,45.2
```

**Minimal limma format**:
```csv
gene,logFC,P.Value,adj.P.Val
GENE1,2.45,1.23e-25,4.56e-22
GENE2,-1.87,3.89e-21,7.12e-18
```

---

### Wilcoxon Format

**Required columns**:

| Data | Accepted Column Names |
|------|----------------------|
| **Gene ID** | `gene_id`, `gene` (or first column) |
| **Fold Change** | `log2FC`, `logFC` |
| **P-value** | `pvalue`, `PValue` |
| **Adjusted P** | `padj`, `FDR` |

**Example Wilcoxon output**:
```csv
gene_id,log2FC,pvalue,padj
GENE1,2.45,1.23e-10,4.56e-8
GENE2,-1.87,3.89e-9,7.12e-7
GENE3,1.92,8.92e-10,3.21e-8
```

---

### What if My Columns Have Different Names?

**Don't worry!** Module 7 is smart and recognizes many column name variations:

**For fold change**, it accepts:
- `log2FoldChange` (DESeq2 style)
- `logFC` (edgeR/limma style)
- `log2FC` (Wilcoxon style)
- `lfc` (lowercase)
- `log_fold_change` (snake_case)

**For p-values**, it accepts:
- `pvalue`, `pval`, `PValue`, `p.value`, `p_value`

**For adjusted p-values**, it accepts:
- `padj` (DESeq2 style)
- `FDR` (edgeR style)
- `adj.P.Val` (limma style)
- `qvalue`, `q.value`

**If your columns still don't match**, you can:
1. Rename them manually in Excel/R before importing
2. Use the `--gene-id-column` option to specify gene column

---

## Auto-Detection Feature

### What is Auto-Detection?

Module 7 can **automatically detect** which DE tool you used based on column names - you don't need to specify `-m`!

**How to use**:
```bash
# Instead of specifying method
raptor import-de -i results.csv -m deseq2

# Let Module 7 detect it automatically
raptor import-de -i results.csv -m auto
```

### How Auto-Detection Works

Module 7 looks at your column names and matches them to patterns:

**DESeq2 detected if**:
- Has `baseMean` column
- Has `log2FoldChange` (DESeq2-specific naming)
- Has `padj` column

**edgeR detected if**:
- Has `logCPM` column
- Has `logFC` AND `PValue` (case-sensitive)
- Has `FDR` column

**limma detected if**:
- Has `AveExpr` column
- Has `P.Value` (with dot)
- Has `adj.P.Val` or `t` statistic

**Wilcoxon detected if**:
- Has `log2FC` (not logFC)
- Simple structure (fewer columns)
- Uses `pvalue` and `padj`

### Auto-Detection Example

```bash
$ raptor import-de -i myresults.csv -m auto

🦖 RAPTOR v2.2.0 - Import DE Results (Module 7)

📂 Loading DE results from: myresults.csv
   Loaded 15,000 genes

   Auto-detected pipeline: DESeq2 ✓
   
   Detected columns:
   - baseMean → DESeq2 signature
   - log2FoldChange → DESeq2 naming
   - padj → DESeq2 FDR column

🔄 Standardizing column names...
```

### When to Use Auto vs Manual

**Use Auto-Detection** (`-m auto`) when:
- ✅ You have standard output from DESeq2/edgeR/limma
- ✅ You're not sure which tool was used
- ✅ You're processing multiple files from different tools

**Specify Method Manually** (`-m deseq2`) when:
- ✅ You modified column names
- ✅ You have custom/non-standard format
- ✅ Auto-detection fails or is uncertain
- ✅ You want to be explicit (more reproducible)

---

## Step-by-Step Tutorials

### Tutorial 1: Import DESeq2 Results

**Scenario**: You ran DESeq2 in R and exported results to CSV.

#### Step 1: Verify Your DESeq2 Results File

```bash
# Check first few lines
head deseq2_results.csv
```

**Expected content**:
```csv
gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
ENSG00000141510,5234.5,2.45,0.234,10.47,1.23e-25,4.56e-22
ENSG00000012048,3456.2,-1.87,0.198,-9.44,3.89e-21,7.12e-18
...
```

#### Step 2: Run Import Command

**Basic import** (default settings):
```bash
raptor import-de \
    -i deseq2_results.csv \
    -m deseq2
```

**Custom settings**:
```bash
raptor import-de \
    -i deseq2_results.csv \
    -m deseq2 \
    -o results/my_deseq2_import \
    --fdr-threshold 0.01 \
    --lfc-threshold 1.0
```

#### Step 3: Examine the Output

```bash
# Check what was created
ls results/de_imported/

# View standardized results
head results/de_imported/de_standardized.csv

# View significant genes only
head results/de_imported/de_significant.csv

# Check summary
cat results/de_imported/de_summary.json
```

**You should see**:
```
results/de_imported/
├── de_standardized.csv   # All genes, standardized format
├── de_significant.csv    # Significant genes only
├── de_summary.json       # Summary statistics
└── de_result.pkl         # For Modules 8-10
```

#### Step 4: Verify Results

**Check summary**:
```bash
cat results/de_imported/de_summary.json
```

**Expected output**:
```json
{
  "n_genes": 15000,
  "n_significant": 1850,
  "n_up": 950,
  "n_down": 900,
  "proportion_significant": 0.123,
  "pipeline": "DESEQ2",
  "parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 0.0
  }
}
```

---

### Tutorial 2: Import edgeR Results

**Scenario**: You ran edgeR exact test or GLM and saved topTags output.

#### Step 1: Export from R

```r
# In R, after running edgeR analysis
library(edgeR)

# Get all results (not just top)
all_results <- topTags(qlf, n=Inf)

# Save to CSV
write.csv(all_results$table, "edger_results.csv")
```

#### Step 2: Import to RAPTOR

```bash
raptor import-de \
    -i edger_results.csv \
    -m edger \
    -o results/edger_import
```

**With strict thresholds**:
```bash
raptor import-de \
    -i edger_results.csv \
    -m edger \
    --fdr-threshold 0.01 \
    --lfc-threshold 0.5
```

#### Step 3: Compare with DESeq2

If you ran both DESeq2 and edgeR:
```bash
# Import both
raptor import-de -i deseq2_results.csv -m deseq2 -o results/deseq2
raptor import-de -i edger_results.csv -m edger -o results/edger

# Compare significant genes
wc -l results/deseq2/de_significant.csv
wc -l results/edger/de_significant.csv
```

---

### Tutorial 3: Import limma-voom Results

**Scenario**: Large RNA-seq study (n>20), used limma-voom for speed.

#### Step 1: Export from R

```r
# In R, after limma-voom analysis
library(limma)

# Get all genes
all_results <- topTable(fit, coef=2, n=Inf, sort.by="none")

# Save with row names as gene column
all_results$gene <- rownames(all_results)
write.csv(all_results, "limma_results.csv", row.names=FALSE)
```

#### Step 2: Import to RAPTOR

```bash
raptor import-de \
    -i limma_results.csv \
    -m limma \
    -o results/limma_import
```

**Auto-detection**:
```bash
# Let RAPTOR detect it's limma
raptor import-de -i limma_results.csv -m auto
```

---

### Tutorial 4: Import Wilcoxon Results

**Scenario**: Non-parametric test for non-normal data.

#### Step 1: Prepare Wilcoxon Results

```r
# Example Wilcoxon test implementation
results <- data.frame(
    gene_id = rownames(counts),
    log2FC = log2FoldChanges,  # Calculate from mean counts
    pvalue = wilcox_pvalues,    # From wilcox.test()
    padj = p.adjust(wilcox_pvalues, method="BH")
)

write.csv(results, "wilcoxon_results.csv", row.names=FALSE)
```

#### Step 2: Import to RAPTOR

```bash
raptor import-de \
    -i wilcoxon_results.csv \
    -m wilcoxon \
    -o results/wilcoxon_import
```

---

### Tutorial 5: Batch Import Multiple Comparisons

**Scenario**: You have multiple comparisons (Treatment1 vs Control, Treatment2 vs Control, etc.)

```bash
# Import each comparison separately
for COMP in treatment1_vs_ctrl treatment2_vs_ctrl treatment3_vs_ctrl; do
    raptor import-de \
        -i results/${COMP}_deseq2.csv \
        -m deseq2 \
        -o results/de_imported_${COMP}
done

# Now you have standardized results for all comparisons
ls results/de_imported_*/de_result.pkl
```

---

### Tutorial 6: Custom Gene ID Column

**Scenario**: Your gene IDs are in a column named "ensembl_id" instead of "gene_id"

```bash
raptor import-de \
    -i results.csv \
    -m deseq2 \
    --gene-id-column ensembl_id \
    -o results/de_imported
```

---

## Understanding the Output

### Output File Structure

After running Module 7, you get **4 files**:

```
results/de_imported/
├── de_standardized.csv      # All genes in standard format
├── de_significant.csv       # Significant genes only
├── de_summary.json          # Summary statistics
└── de_result.pkl            # DEResult object for M8-M10
```

---

### File 1: de_standardized.csv

**What it is**: All genes from your DE analysis in standardized RAPTOR format

**Structure**:
```csv
gene_id,log2_fold_change,p_value,adjusted_p_value,base_mean,statistic,is_significant,direction
GENE1,2.45,1.23e-25,4.56e-22,5234.5,10.47,True,up
GENE2,-1.87,3.89e-21,7.12e-18,3456.2,-9.44,True,down
GENE3,0.23,0.234,0.456,8901.3,1.12,False,unchanged
GENE4,1.92,8.92e-25,3.21e-21,4567.8,10.27,True,up
```

**Columns explained**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `gene_id` | string | Gene identifier | ENSG00000141510 |
| `log2_fold_change` | float | Log2 fold change | 2.45 |
| `p_value` | float | Raw p-value | 1.23e-25 |
| `adjusted_p_value` | float | FDR-adjusted p-value | 4.56e-22 |
| `base_mean` | float | Mean expression (optional) | 5234.5 |
| `statistic` | float | Test statistic (optional) | 10.47 |
| `is_significant` | boolean | Passes thresholds? | True |
| `direction` | string | 'up', 'down', 'unchanged' | up |

**When to use**:
- ✅ Full dataset for visualization (MA plots, volcano plots)
- ✅ Custom filtering with different thresholds
- ✅ Input for statistical tools (IHW, swfdr)
- ✅ Comprehensive gene lists for supplementary materials

---

### File 2: de_significant.csv

**What it is**: Filtered subset containing only significant genes

**Filter criteria**:
```python
is_significant = (
    (adjusted_p_value < fdr_threshold) AND
    (abs(log2_fold_change) > lfc_threshold)
)
```

**Default thresholds**:
- FDR threshold: 0.05 (5%)
- LFC threshold: 0.0 (no fold-change filter)

**Example**:
```csv
gene_id,log2_fold_change,p_value,adjusted_p_value,base_mean,is_significant,direction
GENE1,2.45,1.23e-25,4.56e-22,5234.5,True,up
GENE2,-1.87,3.89e-21,7.12e-18,3456.2,True,down
GENE4,1.92,8.92e-25,3.21e-21,4567.8,True,up
```

**When to use**:
- ✅ Pathway enrichment analysis (GO, KEGG, Reactome)
- ✅ Gene lists for follow-up validation (qRT-PCR)
- ✅ Input for network analysis
- ✅ Quick overview of DE genes

---

### File 3: de_summary.json

**What it is**: JSON file with summary statistics and metadata

**Full example**:
```json
{
  "n_genes": 15000,
  "n_significant": 1850,
  "n_up": 950,
  "n_down": 900,
  "proportion_significant": 0.123,
  "mean_abs_lfc_significant": 2.35,
  "median_adjusted_p_value_significant": 0.003,
  "parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 0.0,
    "source_pipeline": "DESEQ2"
  },
  "metadata": {
    "source_file": "deseq2_results.csv",
    "timestamp": "2026-03-09T14:30:00",
    "raptor_version": "2.2.0",
    "module": "M7"
  },
  "top_upregulated": [
    {"gene_id": "GENE1", "log2_fold_change": 5.67, "adjusted_p_value": 1.2e-45},
    {"gene_id": "GENE5", "log2_fold_change": 4.89, "adjusted_p_value": 3.4e-38}
  ],
  "top_downregulated": [
    {"gene_id": "GENE2", "log2_fold_change": -4.23, "adjusted_p_value": 2.1e-35},
    {"gene_id": "GENE8", "log2_fold_change": -3.98, "adjusted_p_value": 5.6e-32}
  ]
}
```

**When to use**:
- ✅ Quick quality check
- ✅ Comparing multiple DE analyses
- ✅ Tracking analysis provenance
- ✅ Reporting in manuscripts

**Access in Python**:
```python
import json

with open('de_summary.json') as f:
    summary = json.load(f)

print(f"Total genes: {summary['n_genes']}")
print(f"Significant: {summary['n_significant']} ({summary['proportion_significant']:.1%})")
print(f"Upregulated: {summary['n_up']}")
print(f"Downregulated: {summary['n_down']}")
```

---

### File 4: de_result.pkl

**What it is**: Python pickle file containing complete `DEResult` object

**Why it exists**: Seamless integration with Modules 8-10

**Contains**:
- Complete standardized DataFrame
- All metadata
- Analysis parameters
- Convenient methods for filtering, summarizing

**How to use**:
```python
from raptor.de_import import DEResult

# Load
de_result = DEResult.load('de_result.pkl')

# Access data
print(de_result.n_significant)
print(de_result.upregulated_genes[:10])  # Top 10 up genes

# Get summary
print(de_result.summary())

# Filter
filtered = de_result.filter_by_threshold(fdr=0.01, lfc=1.0)

# Pass to Module 8
from raptor.parameter_optimization import optimize_with_fdr_control
optimized = optimize_with_fdr_control(de_result.results_df)
```

**Don't delete this file!** It's needed for Modules 8, 9, and 10.

---

## Working with DEResult Objects

### What is a DEResult Object?

`DEResult` is a Python dataclass that wraps your DE results with convenient methods and properties.

**Structure**:
```python
@dataclass
class DEResult:
    results_df: pd.DataFrame    # Standardized results
    pipeline: str               # 'DESEQ2', 'EDGER', etc.
    parameters: Dict            # FDR/LFC thresholds
    metadata: Dict              # Source file, timestamp
```

---

### Loading a DEResult

**From pickle file**:
```python
from raptor.de_import import DEResult

de_result = DEResult.load('results/de_imported/de_result.pkl')
```

**From import function**:
```python
from raptor.de_import import import_de_results

de_result = import_de_results(
    de_file='deseq2_results.csv',
    pipeline='deseq2',
    fdr_threshold=0.05,
    lfc_threshold=0.0
)
```

---

### Properties (Quick Access)

**Basic counts**:
```python
print(f"Total genes: {de_result.n_genes}")
print(f"Significant: {de_result.n_significant}")
print(f"Upregulated: {de_result.n_up}")
print(f"Downregulated: {de_result.n_down}")
```

**Gene lists**:
```python
# Get all significant gene IDs
sig_genes = de_result.significant_genes  # List of gene IDs

# Get upregulated gene IDs
up_genes = de_result.upregulated_genes

# Get downregulated gene IDs
down_genes = de_result.downregulated_genes
```

**Example output**:
```python
>>> de_result.n_genes
15000
>>> de_result.n_significant
1850
>>> de_result.significant_genes[:5]
['GENE1', 'GENE2', 'GENE4', 'GENE7', 'GENE9']
```

---

### Methods (Actions)

**Get top genes**:
```python
# Top 50 by adjusted p-value (most significant)
top_genes = de_result.get_top_genes(n=50, by='adjusted_p_value')

# Top 20 by fold change (largest changes)
top_fc = de_result.get_top_genes(n=20, by='log2_fold_change')

# Top 10 significant genes only
top_sig = de_result.get_top_genes(n=10, significant_only=True)
```

**Filter by different thresholds**:
```python
# Create new DEResult with stricter thresholds
strict = de_result.filter_by_threshold(fdr=0.01, lfc=1.0)

print(f"Original: {de_result.n_significant} genes")
print(f"Strict: {strict.n_significant} genes")

# Original: 1850 genes
# Strict: 450 genes
```

**Calculate metrics**:
```python
metrics = de_result.calculate_metrics()

print(metrics)
# {
#   'n_genes': 15000,
#   'n_significant': 1850,
#   'proportion_significant': 0.123,
#   'mean_abs_lfc_significant': 2.35,
#   'median_p_value': 0.234,
#   ...
# }
```

**Get summary**:
```python
print(de_result.summary())
```

**Output**:
```
======================================================================
RAPTOR Module 7: DE Results Import
======================================================================
Pipeline: DESEQ2
Total Genes: 15,000
Significant Genes: 1,850 (12.3%)
Upregulated: 950 (6.3%)
Downregulated: 900 (6.0%)
Unchanged: 13,150 (87.7%)

Thresholds:
  FDR: 0.05
  LFC: 0.0

Top Upregulated Genes:
  GENE1: log2FC=5.67, FDR=1.2e-45
  GENE5: log2FC=4.89, FDR=3.4e-38
  GENE9: log2FC=4.56, FDR=7.8e-35

Top Downregulated Genes:
  GENE2: log2FC=-4.23, FDR=2.1e-35
  GENE8: log2FC=-3.98, FDR=5.6e-32
  GENE12: log2FC=-3.67, FDR=1.2e-28

Timestamp: 2026-03-09T14:30:00
======================================================================
```

**Save to new location**:
```python
de_result.save('backup/my_de_result.pkl')
```

---

### Working with the DataFrame

**Access full data**:
```python
# Get the DataFrame
df = de_result.results_df

# Standard pandas operations
print(df.head())
print(df.describe())

# Filter manually
sig = df[df['is_significant']]
high_fc = df[abs(df['log2_fold_change']) > 2]
```

**Export to other formats**:
```python
# Save to CSV
de_result.results_df.to_csv('my_de_results.csv')

# Save to Excel
de_result.results_df.to_excel('my_de_results.xlsx', index=False)

# Save to TSV
de_result.results_df.to_csv('my_de_results.tsv', sep='\t')
```

---

### Integration Examples

**Example 1: Pathway Enrichment**
```python
from raptor.de_import import DEResult

# Load results
de_result = DEResult.load('de_result.pkl')

# Get significant genes for enrichment
sig_genes = de_result.significant_genes

# Use with clusterProfiler in R
import rpy2.robjects as ro
ro.r.assign('gene_list', sig_genes)
ro.r('library(clusterProfiler)')
ro.r('ego <- enrichGO(gene=gene_list, OrgDb=org.Hs.eg.db, ont="BP")')
```

**Example 2: Volcano Plot**
```python
import matplotlib.pyplot as plt

df = de_result.results_df

plt.figure(figsize=(10, 6))
plt.scatter(
    df['log2_fold_change'],
    -np.log10(df['adjusted_p_value']),
    c=df['is_significant'],
    alpha=0.5
)
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(FDR)')
plt.title(f'Volcano Plot - {de_result.pipeline}')
plt.axhline(-np.log10(0.05), color='red', linestyle='--')
plt.axvline(0, color='black', linestyle='-', linewidth=0.5)
plt.show()
```

**Example 3: Heatmap of Top Genes**
```python
import seaborn as sns

# Get top 50 significant genes
top = de_result.get_top_genes(n=50, significant_only=True)
gene_ids = top.index.tolist()

# Assuming you have normalized counts
counts_subset = normalized_counts.loc[gene_ids]

sns.clustermap(counts_subset, z_score=0, cmap='RdBu_r', figsize=(12, 10))
plt.show()
```

---

## Integration with Your Workflow

### Before Module 7: R Analysis

**Step 1: Run DE analysis in R**

Choose your tool and run analysis:

**DESeq2**:
```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, metadata, ~condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "deseq2_results.csv")
```

**edgeR**:
```r
library(edgeR)
y <- DGEList(counts, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
write.csv(topTags(qlf, n=Inf)$table, "edger_results.csv")
```

**limma-voom**:
```r
library(limma)
library(edgeR)
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~condition)
v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
write.csv(topTable(fit, n=Inf, coef=2), "limma_results.csv")
```

---

### During Module 7: Import & Standardize

**Step 2: Import with Module 7**

```bash
raptor import-de \
    -i deseq2_results.csv \
    -m deseq2 \
    -o results/de_imported
```

**What happens**:
1. ✅ Validates file format
2. ✅ Standardizes column names
3. ✅ Calculates significance
4. ✅ Adds direction labels
5. ✅ Creates 4 output files
6. ✅ Generates summary statistics

---

### After Module 7: Downstream Analysis

**Option 1: Proceed to Module 8 (Parameter Optimization)**
```bash
raptor optimize \
    -d results/de_imported/de_result.pkl \
    -m fdr-control \
    -o results/optimization
```

**Option 2: Proceed to Module 9 (Ensemble Analysis)**

If you ran multiple DE tools:
```bash
# Import all methods
raptor import-de -i deseq2_results.csv -m deseq2 -o results/deseq2
raptor import-de -i edger_results.csv -m edger -o results/edger
raptor import-de -i limma_results.csv -m limma -o results/limma

# Ensemble analysis
raptor ensemble \
    --methods fisher brown rra \
    --deseq2 results/deseq2/de_result.pkl \
    --edger results/edger/de_result.pkl \
    --limma results/limma/de_result.pkl \
    -o results/ensemble
```

**Option 3: External Analysis**

Use standardized outputs with other tools:

**Pathway enrichment (R)**:
```r
library(clusterProfiler)

# Load RAPTOR output
genes <- read.csv("results/de_imported/de_significant.csv")
gene_list <- genes$gene_id

# GO enrichment
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH"
)
```

**GSEA (R)**:
```r
library(fgsea)

# Load all genes with fold changes
all_genes <- read.csv("results/de_imported/de_standardized.csv")

# Create ranked list
ranked <- setNames(
    all_genes$log2_fold_change,
    all_genes$gene_id
)
ranked <- sort(ranked, decreasing=TRUE)

# Run GSEA
fgseaRes <- fgsea(pathways, ranked, minSize=15, maxSize=500)
```

**Visualization (Python)**:
```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load results
df = pd.read_csv('results/de_imported/de_standardized.csv')

# MA plot
plt.scatter(df['base_mean'], df['log2_fold_change'], 
            c=df['is_significant'], alpha=0.3)
plt.xscale('log')
plt.xlabel('Mean Expression')
plt.ylabel('log2 Fold Change')
plt.show()

# Volcano plot
plt.scatter(df['log2_fold_change'], 
            -np.log10(df['adjusted_p_value']),
            c=df['is_significant'], alpha=0.3)
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(FDR)')
plt.show()
```

---

### Complete Workflow Example

```bash
# 1. Quality control
raptor qc --counts counts.csv --metadata metadata.csv -o results/qc

# 2. Data profiling
raptor profile --counts counts.csv --metadata metadata.csv -o results/profile

# 3. Get recommendation
raptor recommend --profile results/profile/data_profile.json -o results/recommend

# 4. Run DE analysis in R (DESeq2 recommended)
# [Run R script - produces deseq2_results.csv]

# 5. Import DE results ⭐ MODULE 7
raptor import-de -i deseq2_results.csv -m deseq2 -o results/de_imported

# 6. Optimize parameters
raptor optimize -d results/de_imported/de_result.pkl -m fdr-control -o results/opt

# 7. Ensemble (if multiple methods)
raptor ensemble --methods fisher --deseq2 results/opt/deg_genes.csv -o results/ensemble

# 8. Biomarker discovery
raptor biomarkers --de-result results/ensemble/result.pkl -o results/biomarkers
```

---

## Best Practices

### Before Importing

1. ✅ **Check your R output first**
   ```r
   # In R, verify you have results
   head(res)
   summary(res)
   
   # Check for NAs
   sum(is.na(res$padj))
   ```

2. ✅ **Export all genes, not just significant**
   ```r
   # GOOD: All genes
   write.csv(as.data.frame(res), "results.csv")
   
   # BAD: Only significant genes
   write.csv(subset(res, padj < 0.05), "results.csv")  # Don't do this!
   ```

3. ✅ **Include gene IDs as a column**
   ```r
   # GOOD: Gene IDs in column
   results_df <- as.data.frame(res)
   results_df$gene_id <- rownames(results_df)
   write.csv(results_df, "results.csv", row.names=FALSE)
   
   # OKAY: Gene IDs as row names (RAPTOR handles this)
   write.csv(as.data.frame(res), "results.csv")
   ```

4. ✅ **Remove quotes from gene names if present**
   ```r
   # Use quote=FALSE to avoid quotes around gene names
   write.csv(results, "results.csv", quote=FALSE)
   ```

---

### During Import

1. ✅ **Use consistent thresholds**
   ```bash
   # Use same thresholds for all comparisons
   FDR=0.05
   LFC=0.0
   
   raptor import-de -i treatment1.csv -m deseq2 --fdr-threshold $FDR --lfc-threshold $LFC
   raptor import-de -i treatment2.csv -m deseq2 --fdr-threshold $FDR --lfc-threshold $LFC
   ```

2. ✅ **Organize output directories**
   ```bash
   # Good organization
   raptor import-de -i deseq2_trt1_vs_ctrl.csv -m deseq2 -o results/deseq2/trt1_vs_ctrl
   raptor import-de -i deseq2_trt2_vs_ctrl.csv -m deseq2 -o results/deseq2/trt2_vs_ctrl
   
   # Keeps everything organized by method and comparison
   ```

3. ✅ **Document your imports**
   ```bash
   # Create a log
   echo "Imported DESeq2 results on $(date)" >> import_log.txt
   echo "  File: deseq2_results.csv" >> import_log.txt
   echo "  FDR: 0.05, LFC: 0.0" >> import_log.txt
   echo "  Significant genes: $(wc -l < results/de_imported/de_significant.csv)" >> import_log.txt
   ```

4. ✅ **Verify gene counts match**
   ```bash
   # Check number of genes matches R output
   R_GENES=$(wc -l < deseq2_results.csv)
   RAPTOR_GENES=$(tail -n +2 results/de_imported/de_standardized.csv | wc -l)
   
   if [ "$R_GENES" -eq "$RAPTOR_GENES" ]; then
       echo "✓ Gene counts match: $R_GENES genes"
   else
       echo "✗ Gene count mismatch! R: $R_GENES, RAPTOR: $RAPTOR_GENES"
   fi
   ```

---

### After Import

1. ✅ **Review the summary**
   ```bash
   cat results/de_imported/de_summary.json | jq .
   ```

2. ✅ **Spot-check a few genes**
   ```python
   # In Python
   import pandas as pd
   
   original = pd.read_csv('deseq2_results.csv')
   standardized = pd.read_csv('results/de_imported/de_standardized.csv')
   
   # Check gene GENE1
   print("Original:")
   print(original[original['gene_id'] == 'GENE1'])
   
   print("\nStandardized:")
   print(standardized[standardized['gene_id'] == 'GENE1'])
   
   # Values should match!
   ```

3. ✅ **Visualize to verify**
   ```python
   import matplotlib.pyplot as plt
   
   df = pd.read_csv('results/de_imported/de_standardized.csv')
   
   # Quick volcano plot
   plt.scatter(df['log2_fold_change'], 
                -np.log10(df['adjusted_p_value']),
                c=df['is_significant'], alpha=0.3)
   plt.xlabel('log2 Fold Change')
   plt.ylabel('-log10(FDR)')
   plt.title('Quick QC: Volcano Plot')
   plt.show()
   ```

4. ✅ **Keep original R files**
   ```bash
   # Don't delete original DE results!
   # Keep for:
   # - Verification
   # - Re-import with different parameters
   # - Reviewer requests
   
   mkdir -p results/original_de_outputs
   cp deseq2_results.csv results/original_de_outputs/
   ```

---

### Common Mistakes to Avoid

❌ **Don't**:
1. Export only significant genes from R
2. Use different thresholds for different comparisons
3. Modify gene IDs during export
4. Delete original R output files
5. Import without checking data quality first
6. Mix different gene ID types (Ensembl vs symbols)
7. Forget to specify output directory (overwrites previous)

✅ **Do**:
1. Export all genes (including non-significant)
2. Use consistent thresholds across comparisons
3. Keep gene IDs exactly as they appear in counts
4. Archive original outputs
5. Review summary statistics after import
6. Use same gene ID format throughout
7. Organize outputs by comparison/method

---

## Troubleshooting

### Issue 1: "Column not found" Error

**Error**:
```
ValidationError: Required column 'log2_fold_change' not found in input file
Available columns: ['gene', 'logFC', 'PValue', 'FDR']
```

**Cause**: Module can't find fold change column

**Solutions**:

1. **Check if you used the right method**:
   ```bash
   # Are you sure it's DESeq2? Looks like edgeR
   raptor import-de -i results.csv -m edger  # Try this instead
   ```

2. **Let auto-detection help**:
   ```bash
   raptor import-de -i results.csv -m auto
   ```

3. **Check your column names**:
   ```bash
   head -1 results.csv
   ```

4. **Rename columns in R before export**:
   ```r
   # Rename to standard names
   names(results) <- c('gene_id', 'log2_fold_change', 'p_value', 'adjusted_p_value')
   write.csv(results, "results.csv", row.names=FALSE)
   ```

---

### Issue 2: Gene ID Column Not Recognized

**Error**:
```
ValidationError: No gene ID column found. Expected 'gene_id' or first column.
```

**Cause**: Gene IDs in unusual column name

**Solution**:
```bash
# Specify the gene ID column explicitly
raptor import-de \
    -i results.csv \
    -m deseq2 \
    --gene-id-column ensembl_id
```

Or rename in R:
```r
results$gene_id <- results$ensembl_id
results$ensembl_id <- NULL
```

---

### Issue 3: Different Number of Genes

**Observation**:
```
R output: 20,000 genes
RAPTOR import: 15,000 genes
```

**Cause**: Module removed rows with missing values

**Solution**:

1. **Check for NAs in R**:
   ```r
   # How many genes have NA padj?
   sum(is.na(res$padj))
   
   # DESeq2 sometimes has NAs for:
   # - Low count genes (filtered)
   # - Outliers
   # - Genes with all zeros
   ```

2. **This is normal** - RAPTOR removes invalid rows for safety

3. **If you want to keep them**, fill NAs in R first:
   ```r
   # Fill NA padj with 1.0 (non-significant)
   res$padj[is.na(res$padj)] <- 1.0
   res$pvalue[is.na(res$pvalue)] <- 1.0
   ```

---

### Issue 4: All Genes Marked Non-Significant

**Observation**:
```
de_summary.json shows: "n_significant": 0
```

**Causes**:

1. **Too strict thresholds**:
   ```bash
   # You used very strict thresholds
   raptor import-de -i results.csv -m deseq2 --fdr-threshold 0.0001 --lfc-threshold 5.0
   ```
   
   **Solution**: Use more reasonable thresholds
   ```bash
   raptor import-de -i results.csv -m deseq2 --fdr-threshold 0.05 --lfc-threshold 0.0
   ```

2. **No significant genes in original analysis**:
   ```r
   # Check in R
   summary(res$padj < 0.05)
   # If all FALSE, you really have no significant genes
   ```
   
   **Solution**: This might be biology! Check:
   - Sample size adequate?
   - Biological effect size small?
   - High variability?

3. **Wrong p-value column**:
   Check that adjusted p-values are actually FDR, not raw p-values

---

### Issue 5: Direction Labels Wrong

**Observation**: Upregulated genes showing as "down" or vice versa

**Cause**: Log fold changes inverted (Control vs Treatment instead of Treatment vs Control)

**Solution**:

1. **Check in R**:
   ```r
   # What was your contrast?
   resultsNames(dds)
   res <- results(dds, contrast=c("condition", "treatment", "control"))
   # This gives: Treatment/Control (positive = upregulated in treatment)
   ```

2. **If inverted, fix in R**:
   ```r
   res$log2FoldChange <- -res$log2FoldChange
   write.csv(as.data.frame(res), "results.csv")
   ```

3. **Or note in documentation**: Direction is relative to comparison order

---

### Issue 6: Import is Slow

**Problem**: Import takes >5 minutes for 20,000 genes

**Solutions**:

1. **Check file size**:
   ```bash
   ls -lh results.csv
   # If >100MB, might be including extra columns
   ```

2. **Remove unnecessary columns in R**:
   ```r
   # Keep only essential columns
   results_minimal <- res[, c('log2FoldChange', 'pvalue', 'padj', 'baseMean')]
   write.csv(results_minimal, "results.csv")
   ```

3. **Use CSV not Excel**:
   ```bash
   # Don't export from R to Excel then to CSV
   # Export directly to CSV for speed
   ```

---

### Issue 7: "Validation Error: p-values > 1"

**Error**:
```
ValidationError: Found 150 p-values > 1.0
```

**Cause**: File has incorrect values

**Solution**:

1. **Check in R**:
   ```r
   summary(res$pvalue)
   summary(res$padj)
   # Should all be between 0 and 1
   ```

2. **Check for data corruption**:
   ```bash
   # Look for obviously wrong values
   awk -F',' '{print $4}' results.csv | sort -n | tail -20
   ```

3. **Re-export from R**:
   ```r
   # Clean export
   write.csv(as.data.frame(res), "results_clean.csv", quote=FALSE)
   ```

---

## Frequently Asked Questions

### General Questions

**Q1: Do I need to use Module 7?**

A: Yes, if you want to use Modules 8-10 (Parameter Optimization, Ensemble Analysis, Biomarker Discovery). Module 7 standardizes your DE results into the format these modules expect.

**Q2: Can I use Module 7 with other DE tools not listed?**

A: If your tool outputs similar columns (gene IDs, log fold changes, p-values, FDR), you can try. Use `-m auto` or create a CSV that matches one of the supported formats.

**Q3: What if I already have standardized results?**

A: You can still use Module 7 to:
- Add significance flags (`is_significant`, `direction`)
- Create the `DEResult` object for downstream modules
- Generate summary statistics
- Validate data quality

---

### File Format Questions

**Q4: Does the gene ID have to be in a specific format?**

A: No! Module 7 works with:
- Ensembl IDs (ENSG00000141510)
- Gene symbols (TP53, BRCA1)
- RefSeq IDs (NM_000546)
- Any unique gene identifier

Just be consistent across your analysis.

**Q5: What if my CSV uses semicolons instead of commas?**

A: Convert it first:
```bash
# Convert semicolon to comma
sed 's/;/,/g' results.csv > results_comma.csv
raptor import-de -i results_comma.csv -m deseq2
```

**Q6: Can I use Excel files (.xlsx)?**

A: No, must be CSV or TSV. Convert in Excel:
- File → Save As → CSV (Comma delimited)

Or in R:
```r
library(readxl)
library(writexl)
data <- read_excel("results.xlsx")
write.csv(data, "results.csv", row.names=FALSE)
```

---

### Method-Specific Questions

**Q7: What's the difference between DESeq2 and edgeR import?**

A: Just the column name mapping:
- DESeq2: `log2FoldChange`, `pvalue`, `padj`, `baseMean`
- edgeR: `logFC`, `PValue`, `FDR`, `logCPM`

The standardized output is identical.

**Q8: When should I use Wilcoxon instead of DESeq2/edgeR?**

A: Use Wilcoxon when:
- Data is non-normal
- You have n≥8 samples per group
- You want a non-parametric test

But you still import it the same way!

**Q9: Can I import results from multiple tools and compare?**

A: Yes! That's exactly what ensemble analysis (Module 9) is for:
```bash
raptor import-de -i deseq2.csv -m deseq2 -o results/deseq2
raptor import-de -i edger.csv -m edger -o results/edger
raptor import-de -i limma.csv -m limma -o results/limma

# Then ensemble
raptor ensemble --deseq2 results/deseq2/de_result.pkl \
                --edger results/edger/de_result.pkl \
                --limma results/limma/de_result.pkl
```

---

### Threshold Questions

**Q10: What FDR threshold should I use?**

A: Common choices:
- **0.05 (5%)**: Standard for most studies
- **0.01 (1%)**: Strict, for follow-up validation
- **0.10 (10%)**: Lenient, for exploratory analysis

Default is 0.05, but you can optimize later with Module 8!

**Q11: Should I filter by fold change?**

A: Depends on your goals:
- **LFC = 0**: Include all statistically significant changes (discovery)
- **LFC = 0.5**: ~1.4-fold change minimum (moderate filter)
- **LFC = 1.0**: 2-fold change minimum (strict, common standard)

Default is 0.0 (no fold-change filter).

**Q12: Can I change thresholds after import?**

A: Yes! Two ways:

1. Re-import with different thresholds:
   ```bash
   raptor import-de -i results.csv -m deseq2 --fdr-threshold 0.01 --lfc-threshold 1.0
   ```

2. Use `filter_by_threshold()` in Python:
   ```python
   de_result = DEResult.load('de_result.pkl')
   strict = de_result.filter_by_threshold(fdr=0.01, lfc=1.0)
   ```

---

### Output Questions

**Q13: What's the difference between de_standardized.csv and de_significant.csv?**

A:
- **de_standardized.csv**: ALL genes (15,000)
- **de_significant.csv**: Only significant genes (1,850)

Use standardized for visualization, significant for downstream analysis.

**Q14: What is de_result.pkl and why do I need it?**

A: It's a Python object containing all your results + methods for analysis. You need it for:
- Module 8 (Parameter Optimization)
- Module 9 (Ensemble Analysis)
- Module 10 (Biomarker Discovery)

**Q15: Can I delete the pkl file if I have the CSV?**

A: Not recommended! The pkl file has:
- Convenient methods
- Metadata tracking
- Integration with downstream modules

CSV is for human reading, pkl is for RAPTOR processing.

---

### Integration Questions

**Q16: How do I use imported results with clusterProfiler?**

A:
```r
library(clusterProfiler)

# Load RAPTOR output
genes <- read.csv("results/de_imported/de_significant.csv")

# Extract gene IDs
gene_list <- genes$gene_id

# GO enrichment
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = "BP"
)
```

**Q17: Can I use this with IPA or Ingenuity?**

A: Yes! Use de_significant.csv:
```bash
# IPA needs specific format:
# ID, log2FC, p-value

# Create IPA format
awk -F',' 'NR==1 {print "ID,logFC,pval"} NR>1 {print $1","$2","$4}' \
    results/de_imported/de_significant.csv > for_ipa.csv
```

**Q18: How do I visualize results?**

A: Many options:

**In R**:
```r
library(ggplot2)

df <- read.csv("results/de_imported/de_standardized.csv")

# Volcano plot
ggplot(df, aes(x=log2_fold_change, y=-log10(adjusted_p_value), color=is_significant)) +
    geom_point(alpha=0.5) +
    theme_minimal()
```

**In Python**:
```python
import seaborn as sns

df = pd.read_csv('results/de_imported/de_standardized.csv')

sns.scatterplot(data=df, x='log2_fold_change', y='adjusted_p_value', 
                hue='is_significant', alpha=0.5)
```

---

### Workflow Questions

**Q19: In what order should I run RAPTOR modules?**

A:
```
1. QC (Module 2)
2. Profiler (Module 3)
3. Recommender (Module 4)
   ↓ [Run DE in R]
4. Import (Module 7) ← YOU ARE HERE
5. Optimize (Module 8)
6. Ensemble (Module 9)
7. Biomarkers (Module 10)
```

**Q20: Can I skip modules?**

A: Yes, but:
- Modules 2-4 help you choose the right DE method
- Module 7 is required for Modules 8-10
- Module 8 helps find optimal thresholds
- Modules 9-10 are optional advanced analysis

Minimum workflow: Module 7 → Module 8 → Module 9

---

## Real-World Examples

### Example 1: Breast Cancer Study

**Study Design**:
- 80 breast cancer samples
- ER+ vs ER- comparison
- DESeq2 analysis

**R Analysis**:
```r
library(DESeq2)

# Load data
counts <- read.csv("breast_cancer_counts.csv", row.names=1)
metadata <- read.csv("breast_cancer_metadata.csv")

# DESeq2
dds <- DESeqDataSetFromMatrix(counts, metadata, ~ER_status)
dds <- DESeq(dds)
res <- results(dds, contrast=c("ER_status", "positive", "negative"))

# Export
write.csv(as.data.frame(res), "deseq2_er_positive_vs_negative.csv")
```

**RAPTOR Import**:
```bash
raptor import-de \
    -i deseq2_er_positive_vs_negative.csv \
    -m deseq2 \
    -o results/breast_cancer_de \
    --fdr-threshold 0.05 \
    --lfc-threshold 0.0
```

**Results**:
```
📊 Summary:
   Total genes: 18,500
   Significant: 3,250 (17.6%)
   Upregulated (ER+): 1,680
   Downregulated (ER+): 1,570
   
Top upregulated in ER+:
   ESR1: log2FC=8.9, FDR=1.2e-180
   PGR: log2FC=6.7, FDR=3.4e-145
   GATA3: log2FC=5.4, FDR=8.9e-120
```

**Next Steps**:
```bash
# Optimize parameters
raptor optimize -d results/breast_cancer_de/de_result.pkl -m fdr-control

# Pathway enrichment
# Use de_significant.csv with clusterProfiler
```

---

### Example 2: Drug Response Study

**Study Design**:
- 24 cell lines
- Drug-sensitive vs drug-resistant
- Ran both DESeq2 AND edgeR for comparison

**R Analysis**:
```r
# DESeq2
dds <- DESeq(dds)
res_deseq2 <- results(dds)
write.csv(as.data.frame(res_deseq2), "deseq2_sensitive_vs_resistant.csv")

# edgeR
y <- DGEList(counts, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
write.csv(topTags(qlf, n=Inf)$table, "edger_sensitive_vs_resistant.csv")
```

**RAPTOR Import**:
```bash
# Import both
raptor import-de -i deseq2_sensitive_vs_resistant.csv -m deseq2 -o results/deseq2
raptor import-de -i edger_sensitive_vs_resistant.csv -m edger -o results/edger

# Compare
echo "DESeq2 significant:" $(tail -n +2 results/deseq2/de_significant.csv | wc -l)
echo "edgeR significant:" $(tail -n +2 results/edger/de_significant.csv | wc -l)

# Output:
# DESeq2 significant: 450
# edgeR significant: 520
```

**Comparison Analysis**:
```python
import pandas as pd

deseq2 = pd.read_csv('results/deseq2/de_significant.csv')
edger = pd.read_csv('results/edger/de_significant.csv')

# Find overlap
deseq2_genes = set(deseq2['gene_id'])
edger_genes = set(edger['gene_id'])

overlap = deseq2_genes & edger_genes
print(f"Overlap: {len(overlap)} genes")
print(f"DESeq2 only: {len(deseq2_genes - edger_genes)}")
print(f"edgeR only: {len(edger_genes - deseq2_genes)}")

# Output:
# Overlap: 380 genes (high-confidence!)
# DESeq2 only: 70 genes
# edgeR only: 140 genes
```

---

### Example 3: Time Course Study

**Study Design**:
- 5 time points (0h, 6h, 12h, 24h, 48h)
- 4 replicates each
- Compare each time to baseline (0h)

**R Analysis**:
```r
# Run DESeq2 for each comparison
for (time in c("6h", "12h", "24h", "48h")) {
    dds <- DESeqDataSetFromMatrix(counts, colData, ~timepoint)
    dds <- dds[, dds$timepoint %in% c("0h", time)]
    dds$timepoint <- droplevels(dds$timepoint)
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("timepoint", time, "0h"))
    
    write.csv(as.data.frame(res), 
              paste0("deseq2_", time, "_vs_0h.csv"))
}
```

**RAPTOR Import** (all comparisons):
```bash
#!/bin/bash

# Import all time points
for TIME in 6h 12h 24h 48h; do
    raptor import-de \
        -i deseq2_${TIME}_vs_0h.csv \
        -m deseq2 \
        -o results/timecourse_${TIME}_vs_0h \
        --fdr-threshold 0.05
done

# Summarize
echo "Time-course DE gene counts:"
for TIME in 6h 12h 24h 48h; do
    COUNT=$(tail -n +2 results/timecourse_${TIME}_vs_0h/de_significant.csv | wc -l)
    echo "  ${TIME}: ${COUNT} genes"
done
```

**Results**:
```
Time-course DE gene counts:
  6h: 145 genes
  12h: 892 genes
  24h: 2150 genes
  48h: 3420 genes

Progressive increase in DE genes over time!
```

---

### Example 4: Multi-Platform Comparison

**Study Design**:
- Same samples analyzed with 3 methods
- DESeq2, edgeR, limma-voom
- Goal: Find consensus genes

**RAPTOR Workflow**:
```bash
# Import all three
raptor import-de -i deseq2_results.csv -m deseq2 -o results/deseq2
raptor import-de -i edger_results.csv -m edger -o results/edger
raptor import-de -i limma_results.csv -m limma -o results/limma

# Extract significant genes
cut -d',' -f1 results/deseq2/de_significant.csv | tail -n +2 > deseq2_genes.txt
cut -d',' -f1 results/edger/de_significant.csv | tail -n +2 > edger_genes.txt
cut -d',' -f1 results/limma/de_significant.csv | tail -n +2 > limma_genes.txt

# Find consensus (in all 3)
cat deseq2_genes.txt edger_genes.txt limma_genes.txt | \
    sort | uniq -c | awk '$1==3 {print $2}' > consensus_genes.txt

echo "Consensus genes (all 3 methods): $(wc -l < consensus_genes.txt)"
```

**Or use Module 9 Ensemble**:
```bash
raptor ensemble \
    --methods fisher \
    --deseq2 results/deseq2/de_result.pkl \
    --edger results/edger/de_result.pkl \
    --limma results/limma/de_result.pkl \
    -o results/ensemble
```

---

## Technical Details

### Column Name Mapping

**Complete mapping table**:

| Standard Name | DESeq2 | edgeR | limma | Wilcoxon |
|---------------|--------|-------|-------|----------|
| `log2_fold_change` | log2FoldChange | logFC | logFC | log2FC |
| `p_value` | pvalue | PValue | P.Value | pvalue |
| `adjusted_p_value` | padj | FDR | adj.P.Val | padj |
| `base_mean` | baseMean | logCPM | AveExpr | - |
| `statistic` | stat | LR | t | - |

### Significance Calculation

**Algorithm**:
```python
def calculate_significance(df, fdr_threshold, lfc_threshold):
    # Significance flag
    df['is_significant'] = (
        (df['adjusted_p_value'] < fdr_threshold) &
        (abs(df['log2_fold_change']) > lfc_threshold)
    )
    
    # Direction
    df['direction'] = 'unchanged'
    df.loc[
        (df['is_significant']) & (df['log2_fold_change'] > lfc_threshold),
        'direction'
    ] = 'up'
    df.loc[
        (df['is_significant']) & (df['log2_fold_change'] < -lfc_threshold),
        'direction'
    ] = 'down'
    
    return df
```

### Auto-Detection Logic

**Priority**:
1. Check for DESeq2-specific columns (baseMean + log2FoldChange)
2. Check for edgeR-specific columns (logCPM + logFC + PValue)
3. Check for limma-specific columns (AveExpr + P.Value + adj.P.Val)
4. Default to Wilcoxon if simple structure

**Code**:
```python
def detect_pipeline(df):
    columns = set(df.columns)
    
    # DESeq2
    if 'baseMean' in columns and 'log2FoldChange' in columns:
        return 'deseq2'
    
    # edgeR
    if 'logCPM' in columns and 'logFC' in columns and 'PValue' in columns:
        return 'edger'
    
    # limma
    if 'AveExpr' in columns and 'P.Value' in columns:
        return 'limma'
    
    # Wilcoxon
    return 'wilcoxon'
```

### Data Validation

**Checks performed**:
1. ✅ Required columns present
2. ✅ No duplicate gene IDs
3. ✅ P-values in range [0, 1]
4. ✅ No infinite fold changes
5. ✅ At least 100 genes (warning if fewer)
6. ✅ Reasonable fold change distribution

### Performance

**Typical processing times**:
- 10,000 genes: ~2 seconds
- 20,000 genes: ~3 seconds
- 50,000 genes: ~8 seconds

**Memory usage**:
- ~50 MB for 20,000 genes
- Scales linearly with gene count

---

## Conclusion

Module 7 (DE Import) is your **gateway to RAPTOR's advanced analysis modules**. It:

✅ **Standardizes** results from 4 popular DE tools  
✅ **Validates** data quality automatically  
✅ **Enriches** with significance flags and direction  
✅ **Prepares** data for Modules 8-10  
✅ **Simplifies** downstream analysis  

### Quick Recap

1. **Run DE analysis in R** (DESeq2, edgeR, limma, Wilcoxon)
2. **Export to CSV** (all genes, with gene IDs)
3. **Import with Module 7**:
   ```bash
   raptor import-de -i results.csv -m deseq2
   ```
4. **Use standardized outputs** for downstream analysis

### What's Next?

After importing:

**Option 1**: Optimize parameters (Module 8)
```bash
raptor optimize -d de_result.pkl -m fdr-control
```

**Option 2**: Ensemble analysis (Module 9)
```bash
raptor ensemble --methods fisher --deseq2 de_result.pkl
```

**Option 3**: External tools
- Pathway enrichment (clusterProfiler, GSEA)
- Visualization (R, Python)
- Network analysis (STRING, Cytoscape)

---

**Thank you for using RAPTOR Module 7!**

For questions or issues:
- **GitHub**: https://github.com/AyehBlk/RAPTOR
- **Email**: ayehbolouki1988@gmail.com
- **Documentation**: Full docs at RAPTOR website

**Happy analyzing!** 🦖📊

---

**Version**: 2.2.0  
**Last Updated**: March 9, 2026  
**Author**: Ayeh Bolouki  
**License**: MIT
