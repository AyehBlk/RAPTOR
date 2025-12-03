# Tutorial 3: Working with Your Own RNA-seq Data

**Level**: Intermediate  
**Time**: 1-2 hours  
**Goal**: Apply RAPTOR to your real RNA-seq experiment

---

## What You'll Learn

- How to prepare your data for RAPTOR
- How to handle different data formats
- How to troubleshoot common data issues
- How to integrate RAPTOR into your workflow

---

## Prerequisites

- Completed [Tutorial 1: Getting Started](tutorial_01_getting_started.md)
- Your own RNA-seq data (FASTQ or count matrix)
- Metadata about your samples
- Reference genome/transcriptome (if starting from FASTQ)

---

## Understanding Your Starting Point

RAPTOR can work with data at different stages:

```
Your Data Stage              RAPTOR Can Help
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FASTQ files                  ✓ Full pipeline (align → quantify → stats)
↓
BAM files                    ✓ From quantification onward
↓
Count matrix                 ✓ Profile & recommend pipeline
↓
Already analyzed             ✓ Compare your results with benchmarks
```

This tutorial covers the **most common case**: You have a count matrix.

---

## Step 1: Organize Your Data

### Create a Project Directory

```bash
mkdir my_rnaseq_project
cd my_rnaseq_project

# Create subdirectories
mkdir data/
mkdir results/
mkdir config/
mkdir reports/
```

### Expected Structure

```
my_rnaseq_project/
├── data/
│   ├── counts.csv           # Your count matrix
│   ├── metadata.csv         # Sample information
│   └── reference/           # (Optional) Reference files
├── results/                 # RAPTOR outputs go here
├── config/                  # Custom configurations
└── reports/                 # Generated reports
```

---

## Step 2: Prepare Your Count Matrix

### Format Requirements

RAPTOR expects a CSV file with:
- **Rows**: Genes
- **Columns**: Samples
- **First column**: Gene IDs
- **No missing values**

#### Example Format

```csv
gene_id,Sample_1,Sample_2,Sample_3,Sample_4,Sample_5,Sample_6
ENSG00000000003,523,612,498,1250,1180,1340
ENSG00000000005,89,95,102,210,185,198
ENSG00000000419,2341,2567,2234,2401,2389,2456
ENSG00000000457,156,178,145,334,298,312
ENSG00000000460,1023,1156,987,1045,1123,1078
...
```

### Converting from Common Formats

#### From RSEM/Salmon Output

```r
# In R
library(tximport)

# List sample files
files <- list.files("quant_results/", pattern="quant.sf", full.names=TRUE)
names(files) <- c("Sample_1", "Sample_2", ...)  # Your sample names

# Import and aggregate to gene level
tx2gene <- read.csv("tx2gene.csv")  # Transcript-to-gene mapping
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)

# Export count matrix
counts <- txi$counts
write.csv(counts, "data/counts.csv", row.names=TRUE)
```

#### From HTSeq/featureCounts Output

```python
# In Python
import pandas as pd

# Read individual count files
sample_files = {
    'Sample_1': 'htseq_Sample_1.txt',
    'Sample_2': 'htseq_Sample_2.txt',
    # ... more samples
}

# Combine into matrix
dfs = []
for sample, file in sample_files.items():
    df = pd.read_csv(file, sep='\t', names=['gene_id', sample], 
                     skiprows=5, index_col=0)  # Skip HTSeq summary rows
    dfs.append(df)

counts = pd.concat(dfs, axis=1)
counts.to_csv('data/counts.csv')
```

#### From DESeq2/edgeR Objects

```r
# In R - from DESeq2
library(DESeq2)
dds <- readRDS("my_deseq2_object.rds")
counts <- counts(dds, normalized=FALSE)  # Get raw counts
write.csv(counts, "data/counts.csv")

# In R - from edgeR
library(edgeR)
y <- readRDS("my_edger_object.rds")
counts <- y$counts
write.csv(counts, "data/counts.csv")
```

### Data Validation

```bash
# Quick check of your count matrix
head -n 5 data/counts.csv
wc -l data/counts.csv  # Count lines (genes)
```

**Common issues to check:**
- ✓ First row contains sample names
- ✓ First column contains gene IDs
- ✓ All values are integers (or can be rounded)
- ✓ No NA or missing values
- ✓ Reasonable value ranges (not already normalized)

---

## Step 3: Prepare Your Metadata

### Metadata Requirements

Essential columns:
- `sample`: Sample identifier (must match count matrix columns)
- `condition`: Experimental condition (e.g., Control, Treatment)

Optional but recommended:
- `replicate`: Biological replicate number
- `batch`: Batch information (if applicable)
- `pairing`: For paired samples
- Any other covariates

#### Example Metadata

```csv
sample,condition,replicate,batch,sex,age
Sample_1,Control,1,Batch1,F,45
Sample_2,Control,2,Batch1,M,52
Sample_3,Control,3,Batch2,F,48
Sample_4,Treatment,1,Batch1,M,51
Sample_5,Treatment,2,Batch2,F,47
Sample_6,Treatment,3,Batch2,M,49
```

### Metadata Validation

```python
# Check metadata
import pandas as pd

meta = pd.read_csv('data/metadata.csv')

# Verify sample names match count matrix
counts = pd.read_csv('data/counts.csv', index_col=0)
print("Samples in counts:", list(counts.columns))
print("Samples in metadata:", list(meta['sample']))

# Check for required columns
required = ['sample', 'condition']
missing = set(required) - set(meta.columns)
if missing:
    print(f"WARNING: Missing columns: {missing}")
else:
    print("✓ All required columns present")

# Check for balanced design
print("\nSamples per condition:")
print(meta['condition'].value_counts())
```

---

## Step 4: Initial Data Exploration

Before running RAPTOR, let's understand your data:

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load data
counts = pd.read_csv('data/counts.csv', index_col=0)
metadata = pd.read_csv('data/metadata.csv')

# Basic statistics
print("Dataset Overview:")
print(f"  Genes: {counts.shape[0]:,}")
print(f"  Samples: {counts.shape[1]}")
print(f"  Total counts: {counts.sum().sum():,.0f}")

# Library sizes
lib_sizes = counts.sum(axis=0)
print(f"\nLibrary Sizes:")
print(f"  Mean: {lib_sizes.mean():,.0f}")
print(f"  Range: {lib_sizes.min():,.0f} - {lib_sizes.max():,.0f}")
print(f"  CV: {lib_sizes.std() / lib_sizes.mean():.2%}")

# Gene detection
genes_detected = (counts > 0).sum(axis=1)
print(f"\nGene Detection:")
print(f"  Avg samples per gene: {genes_detected.mean():.1f}")
print(f"  Genes in all samples: {(genes_detected == counts.shape[1]).sum()}")

# Zero inflation
zero_pct = (counts == 0).sum().sum() / counts.size
print(f"\nZero Inflation: {zero_pct:.1%}")
```

### Expected Ranges

| Metric | Typical Range | Your Value | Status |
|--------|---------------|------------|--------|
| Library Size | 10M - 50M | ? | ? |
| CV of Library Size | 0.1 - 0.3 | ? | ? |
| Zero Inflation | 30% - 60% | ? | ? |
| Genes Detected | 15,000 - 20,000 | ? | ? |

---

## Step 5: Run RAPTOR Profiling

Now let's see what RAPTOR recommends:

```bash
raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --output results/profile/ \
  --threads 8
```

### Understanding the Output

RAPTOR creates:
```
results/profile/
├── raptor_profile_report.html    # Main report
├── data_characteristics.csv      # Detailed metrics
├── plots/
│   ├── library_sizes.png
│   ├── pca.png
│   ├── correlation_heatmap.png
│   └── bcv_plot.png
└── recommendations.json          # Machine-readable results
```

---

## Step 6: Interpret Results for Your Data

### Case Study 1: Small Pilot Study

**Your data:**
- 6 samples (3 vs 3)
- BCV: 0.35 (moderate variation)
- Library size CV: 0.12 (good)

**RAPTOR recommendation:**
```
🥇 Pipeline 3 (Salmon-edgeR)
   - Good for small samples
   - Robust to moderate variation
   - Fast turnaround
```

### Case Study 2: Large Clinical Cohort

**Your data:**
- 48 samples (24 vs 24)
- BCV: 0.65 (high variation - expected for human samples)
- Some batch effects detected

**RAPTOR recommendation:**
```
🥇 Pipeline 1 (STAR-RSEM-DESeq2)
   - Handles high variation well
   - Built-in batch effect handling
   - Worth the computational cost
```

### Case Study 3: Very Limited Replication

**Your data:**
- 4 samples (2 vs 2)
- BCV: Cannot be reliably estimated

**RAPTOR recommendation:**
```
🥇 Pipeline 6 (STAR-featureCounts-NOISeq)
   - Non-parametric statistics
   - Designed for n<3
   - Conservative but appropriate
```

---

## Step 7: Address Common Data Issues

### Issue 1: Library Size Variation

**Problem**: CV > 0.3, samples have very different total counts

**Solution**:
```bash
# RAPTOR will handle this with normalization
# But you can pre-filter low-quality samples

raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --filter-low-quality \
  --min-lib-size 5000000  # 5M reads minimum
```

### Issue 2: Batch Effects

**Problem**: Samples cluster by batch, not condition

**Solution**:
```bash
# Include batch in metadata
# RAPTOR will recommend appropriate methods

raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --batch-column batch \
  --output results/profile/
```

**In the report**, look for:
- PCA colored by batch
- Recommendation will favor pipelines with batch correction (DESeq2, limma)

### Issue 3: Outlier Samples

**Problem**: One or more samples are very different

**Solution**:
```python
# Identify outliers
from scipy.stats import zscore

lib_sizes = counts.sum(axis=0)
z_scores = zscore(lib_sizes)

outliers = lib_sizes[abs(z_scores) > 3]
print("Potential outliers:", outliers.index.tolist())
```

Then either:
```bash
# Remove outliers
raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --exclude-samples Sample_7,Sample_12
```

Or:
```bash
# Use robust methods
# RAPTOR will recommend NOISeq or limma with robust=TRUE
```

### Issue 4: Low Gene Detection

**Problem**: < 10,000 genes detected

**Possible causes**:
- Poor sequencing depth
- rRNA contamination
- Wrong reference genome
- Wrong organism

**Solution**:
```bash
# Check data quality first
raptor qc --counts data/counts.csv

# If depth is issue, consider:
# - Resequencing
# - Pooling biological replicates
# - Using more sensitive methods (RAPTOR will recommend these)
```

---

## Step 8: Customize Analysis

### Custom Configuration

Create `config/my_config.yaml`:

```yaml
# Adjust for your specific needs
profiling:
  bcv_thresholds:
    low: 0.15      # Stricter for cell lines
    medium: 0.35
    high: 0.55
    
statistics:
  fdr_threshold: 0.01  # More stringent
  log2fc_threshold: 1.5  # Larger fold change

recommendation:
  weights:
    accuracy: 0.6     # Prioritize accuracy over speed
    speed: 0.1
    memory: 0.1
    data_compatibility: 0.2
```

Use it:
```bash
raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --config config/my_config.yaml \
  --output results/custom_profile/
```

---

## Step 9: Integration with Existing Workflow

### Scenario: You already have DESeq2 results

```bash
# Compare RAPTOR recommendation with your current approach
raptor profile --counts data/counts.csv --metadata data/metadata.csv

# Then benchmark to see if another pipeline would be better
raptor compare \
  --data data/ \
  --output benchmark/ \
  --pipelines 1 3 5  # Include your current method
```

### Scenario: You want to try multiple thresholds

```bash
# Generate recommendations for different stringency levels
for fdr in 0.01 0.05 0.10; do
  raptor profile \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    --fdr $fdr \
    --output results/fdr_$fdr/
done

# Compare results
raptor compare-profiles --profiles results/fdr_*/
```

---

## Step 10: Generate Analysis Report

### For Your Lab Notebook

```bash
raptor report \
  --results results/profile/ \
  --format notebook \
  --output reports/analysis_notebook.html
```

### For Your Manuscript

```bash
raptor report \
  --results results/profile/ \
  --export-methods \
  --output reports/methods_section.txt
```

**Example output:**
```
RNA-seq data from [YOUR EXPERIMENT] consisting of [N] samples 
(condition A: n=[X], condition B: n=[Y]) were analyzed using 
RAPTOR v2.0.0 (Bolouki, 2025). Data profiling revealed [moderate/high] 
biological coefficient of variation (BCV = [VALUE]) and [good/variable] 
library size distribution (CV = [VALUE]). Based on these characteristics, 
Pipeline [X] ([TOOLS]) was selected for differential expression analysis...
```

---

## Troubleshooting Real Data

### Problem: Recommendation seems wrong

**Investigate:**
```bash
# Re-run with verbose output
raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --verbose \
  --output results/profile_verbose/

# Check detailed metrics
cat results/profile_verbose/data_characteristics.csv
```

**Common reasons:**
- Mislabeled metadata
- Wrong normalization applied before RAPTOR
- Outlier samples affecting metrics

### Problem: "Data quality warning" message

**This means:**
- Library sizes very uneven
- Too many zeros
- Possible sequencing failure

**Actions:**
1. Check QC plots in report
2. Identify problem samples
3. Consider excluding or re-sequencing

### Problem: Multiple conditions (>2 groups)

```bash
# RAPTOR handles this automatically
# Just make sure metadata reflects all groups

raptor profile \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --design "~ condition"  # or more complex: "~ batch + condition"
```

---

## Best Practices Checklist

Before running RAPTOR on your data:

- [ ] Count matrix has genes as rows, samples as columns
- [ ] All values are raw counts (not normalized)
- [ ] No missing values (NA)
- [ ] Metadata sample IDs exactly match count matrix columns
- [ ] At least 2 replicates per condition (3+ strongly recommended)
- [ ] Library sizes reasonably similar (CV < 0.3)
- [ ] Outliers identified and handled
- [ ] Batch information included if applicable
- [ ] Reference genome/transcriptome version documented

---

## Summary

You've learned to:
- ✅ Prepare real RNA-seq data for RAPTOR
- ✅ Convert from various input formats
- ✅ Validate data quality
- ✅ Address common data issues
- ✅ Interpret recommendations in context
- ✅ Customize analysis for your needs
- ✅ Integrate RAPTOR into existing workflows

---

## Next Steps

1. **Run the recommended pipeline** on your full dataset
2. **Compare with your current method** (if applicable)
3. **Share results** with your collaborators
4. **Cite RAPTOR** in your publications

---

## Getting Help

**Still having issues?**

1. Check [FAQ.md](../FAQ.md)
2. Check [TROUBLESHOOTING.md](../TROUBLESHOOTING.md)
3. Post on [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)
4. Email: ayehbolouki1988@gmail.com

---

**Tutorial by Ayeh Bolouki**  
University of Namur, Belgium  
For RAPTOR v2.0.0
