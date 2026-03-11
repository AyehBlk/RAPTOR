# RAPTOR Module 1: Quick Quantification

**Version 2.2.0** | **Stage 1** | **Module M1**

---

## Table of Contents

1. [Overview](#overview)
2. [Philosophy & Design](#philosophy--design)
3. [Quick Start](#quick-start)
4. [Tool Selection: Kallisto vs Salmon](#tool-selection-kallisto-vs-salmon)
5. [Input Requirements](#input-requirements)
6. [Configuration](#configuration)
7. [Validation System](#validation-system)
8. [Output Files](#output-files)
9. [Usage Examples](#usage-examples)
10. [Critical Parameters](#critical-parameters)
11. [Integration with RAPTOR Workflow](#integration-with-raptor-workflow)
12. [Troubleshooting](#troubleshooting)
13. [Best Practices](#best-practices)
14. [Technical Details](#technical-details)

---

## Overview

### What is Module 1?

Module 1 provides **ultra-fast transcript quantification** for quality control and data profiling. It's the entry point of the RAPTOR workflow, designed to give you rapid insights into your RNA-seq data before committing to computationally intensive downstream analyses.

### Key Characteristics

- ⚡ **Speed-Optimized**: Pseudo-alignment instead of full alignment
- 🎯 **Purpose-Built**: QC and profiling, NOT final differential expression
- 🔧 **Dual Implementation**: Both Kallisto and Salmon pipelines
- ✅ **Comprehensive Validation**: 8-10 pre-flight checks prevent common errors
- 📊 **Architecture-Compliant**: Outputs integrate seamlessly with Modules 2-4

### What Module 1 Does

✅ Generate gene-level count matrices in minutes  
✅ Provide quantification QC metrics  
✅ Enable rapid data profiling (Module 3)  
✅ Support pipeline recommendation (Module 4)  
✅ Detect outliers and quality issues (Module 2)

### What Module 1 Does NOT Do

❌ Replace production-grade quantification (use Module 5)  
❌ Generate publication-quality DE analysis  
❌ Perform transcript-level analysis with bootstraps  
❌ Include extensive bias modeling (kept minimal for speed)

---

## Philosophy & Design

### The Two-Stage Quantification Strategy

RAPTOR employs a unique **two-stage quantification approach**:

#### Stage 1: Quick Quantification (Module 1)
**Purpose**: Rapid exploration and decision-making  
**Tools**: Kallisto or Salmon (speed-optimized)  
**Output**: `quick_gene_counts.csv`  
**Use Case**: QC, profiling, pipeline selection

#### Stage 2: Production Quantification (Module 5)
**Purpose**: Final analysis and publication  
**Tools**: Full STAR/RSEM or Salmon with bootstraps  
**Output**: `production_counts.csv`  
**Use Case**: Differential expression, manuscript figures

### Why This Matters

Traditional RNA-seq workflows force you to commit hours or days to quantification before discovering:
- Batch effects requiring re-processing
- Outlier samples needing removal
- Unsuitable pipeline for your data characteristics

**RAPTOR's approach**: Invest 30 minutes in Module 1, make informed decisions, then run optimized production analysis.

### Design Principles

1. **Speed First**: Sacrifices precision for rapid turnaround
2. **Fail Fast**: Comprehensive validation catches errors early
3. **Distinct Outputs**: Clear separation from production counts
4. **Tool Flexibility**: Choose Kallisto or Salmon based on needs
5. **Workflow Integration**: Seamless handoff to Modules 2-4

---

## Quick Start

### Minimal Working Example (Paired-End)

```bash
# Using Kallisto
raptor quick-count -m kallisto \
    -s samples.csv \
    -i kallisto_index.idx \
    -g tx2gene.csv

# Using Salmon  
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/ \
    -g tx2gene.csv
```

### What Happens

1. **Validation** (5-10 seconds): Checks all inputs
2. **Quantification** (2-10 min): Pseudo-aligns reads
3. **Aggregation** (30 sec): Combines to gene-level matrix
4. **QC Report** (5 sec): Generates HTML summary

### Output Location

```
results/quick_counts/
├── quick_gene_counts.csv    # Main output for M2-M4
├── quick_tpm.csv             # Normalized expression
├── sample_info.csv           # Sample metadata + QC
├── qc_report.html            # Visual QC summary
└── count_summary.txt         # Statistics
```

---

## Tool Selection: Kallisto vs Salmon

### When to Use Kallisto

✅ **Maximum Speed Required**
- Large cohorts (100+ samples)
- Exploratory analysis on HPC with time limits
- Testing multiple parameters quickly

✅ **Memory-Constrained Environments**
- Typical usage: 4-8 GB RAM
- Minimal index size (~2-4 GB for human)

✅ **Simple Data**
- Standard RNA-seq without complex biases
- Well-annotated organisms with clean transcriptomes

### When to Use Salmon

✅ **Better Accuracy Needed** (even for QC)
- GC bias correction built-in
- Mapping validation reduces false positives
- Better handling of sequence-specific biases

✅ **Single-End Reads**
- Auto-estimates fragment length distribution
- More robust than Kallisto for SE data

✅ **Complex Samples**
- Degraded RNA
- Samples with known technical artifacts
- Non-model organisms with incomplete annotations

### Performance Comparison

| Feature | Kallisto | Salmon |
|---------|----------|--------|
| **Speed** | ⚡⚡⚡ Fastest | ⚡⚡ Fast |
| **Memory** | 4-8 GB | 8-16 GB |
| **Accuracy (QC)** | ✓ Good | ✓✓ Better |
| **GC Bias Correction** | ❌ No | ✅ Yes |
| **Mapping Validation** | ❌ No | ✅ Yes |
| **Single-End Support** | ⚠️ Manual params | ✅ Auto-detect |
| **Index Size** | ~2 GB (human) | ~4 GB (human) |
| **Bootstraps** | ✓ Supported | ✓ Supported |

### RAPTOR Recommendation

**Default Choice**: **Salmon**
- Better accuracy without significant speed penalty
- More robust across diverse data types
- Handles edge cases better

**Use Kallisto When**:
- Speed is absolutely critical (>200 samples)
- Memory is severely limited (<8 GB)
- You have clean, high-quality paired-end data

---

## Input Requirements

### 1. Sample Sheet (Required)

**Format**: CSV file with specific columns

**Required Columns**:
```csv
sample_id,condition,batch,fastq_r1,fastq_r2
Sample1,Control,,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
Sample2,Control,,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
Sample3,Treatment,,/path/to/Sample3_R1.fastq.gz,/path/to/Sample3_R2.fastq.gz
```

**Column Descriptions**:
- `sample_id`: Unique identifier (no spaces or special chars)
- `condition`: Experimental group (can be empty, filled in Module 2)
- `batch`: Batch identifier (can be empty, filled in Module 2)
- `fastq_r1`: Full path to Read 1 FASTQ file
- `fastq_r2`: Full path to Read 2 FASTQ file (empty for single-end)

**Auto-Generation**:
```bash
# Automatically detect and create sample sheet
python3 detect_samples.py \
    --input /path/to/fastq_directory/ \
    --output samples.csv
```

### 2. Quantification Index (Required)

#### Kallisto Index

**File**: Single `.idx` file  
**Size**: ~2-4 GB for human transcriptome

**Build Command**:
```bash
# Download transcriptome
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Build index
kallisto index \
    -i kallisto_index.idx \
    Homo_sapiens.GRCh38.cdna.all.fa.gz
```

#### Salmon Index

**File**: Directory containing index files  
**Size**: ~4-8 GB for human transcriptome

**Build Command**:
```bash
# Basic index
salmon index \
    -t Homo_sapiens.GRCh38.cdna.all.fa.gz \
    -i salmon_index

# With decoys (recommended for accuracy)
salmon index \
    -t gentrome.fa.gz \
    -d decoys.txt \
    -i salmon_index
```

### 3. Transcript-to-Gene Mapping (Recommended)

**Purpose**: Aggregate transcript-level counts to gene-level

**Formats Supported**:
- TSV/CSV: Two columns (transcript_id, gene_id)
- GTF/GFF: Extracts from annotation file

**Example TSV**:
```tsv
transcript_id	gene_id
ENST00000456328	ENSG00000223972
ENST00000450305	ENSG00000223972
ENST00000488147	ENSG00000227232
```

**From GTF** (automatic):
```bash
# Module 1 can parse GTF directly
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/ \
    -g Homo_sapiens.GRCh38.109.gtf
```

**Without tx2gene**:
- Output will be at transcript level
- Not suitable for Modules 2-4 (which expect genes)
- Use only for transcript-specific analysis

### 4. FASTQ Files (Required)

**Formats Supported**:
- `.fastq.gz` (compressed, recommended)
- `.fq.gz` (compressed)
- `.fastq` (uncompressed)
- `.fq` (uncompressed)

**Naming Conventions** (auto-detected):
- Illumina: `Sample_S1_L001_R1_001.fastq.gz`
- Generic: `Sample_R1.fastq.gz`, `Sample_1.fastq.gz`

**Quality Requirements**:
- Minimum: Phred 20
- Recommended: Phred 30+
- Adapters trimmed (not required but recommended)

---

## Configuration

### Kallisto Configuration (`config.yaml`)

```yaml
# =============================================================================
# Quick Kallisto Configuration
# =============================================================================

kallisto:
  # Strand-specific options
  strand: ""  # Options: "", "rf-stranded", "fr-stranded"
  
  # Bootstrap samples (0 = disabled for speed)
  num_bootstraps: 0
  
  # Random seed for reproducibility
  seed: 42

# CRITICAL for single-end reads
single_end:
  fragment_length: 200  # REQUIRED: Mean fragment length
  fragment_sd: 20       # REQUIRED: Standard deviation

resources:
  threads: 8
  memory_gb: 8
  parallel_samples: false

output:
  dir: "results/quick_counts"
  keep_quant_files: false  # Delete intermediate files
```

### Salmon Configuration (`config.yaml`)

```yaml
# =============================================================================
# Quick Salmon Configuration
# =============================================================================

salmon:
  # Library type detection
  lib_type: "A"  # A = auto-detect (recommended)
  
  # Accuracy features
  validate_mappings: true  # Improves accuracy
  gc_bias: true           # GC bias correction
  seq_bias: false         # Sequence bias (slower)
  
  # Bootstrap samples (0 = disabled for speed)
  num_bootstraps: 0

single_end:
  fragment_length_mean: 200
  fragment_length_sd: 80

resources:
  threads: 8
  memory_gb: 16
  parallel_samples: false

output:
  dir: "results/quick_counts"
  keep_quant_files: false
```

### Key Configuration Decisions

#### Bootstraps: Why Disabled?

**Module 1 (Quick)**:
```yaml
num_bootstraps: 0  # Speed optimized
```
- QC and profiling don't need uncertainty estimates
- Saves 10-20x computation time
- Uncertainty not used in Modules 2-4

**Module 5 (Production)**:
```yaml
num_bootstraps: 100  # Precision for DE analysis
```
- Required for tools like Sleuth, DESeq2
- Enables confidence intervals
- Publication-quality results

#### Thread Allocation

**Sweet Spot**: 8 threads per sample
- Diminishing returns after 16 threads
- Better to run samples in parallel with 8 threads each

**Scaling**:
```yaml
# For large cohorts
threads: 4-8
parallel_samples: true  # Run 4-8 samples simultaneously

# For deep samples
threads: 16-32
parallel_samples: false  # Focus resources on one sample
```

---

## Validation System

### Version 2.2.0 Enhancements

Module 1 now includes **comprehensive pre-flight validation** to catch errors before wasting computation time.

### Validation Checks

#### Kallisto Pipeline (10 Checks)

1. ✅ **Sample sheet exists**
   - Error: File not found or path typo
   - Fix: Verify file path

2. ✅ **Sample sheet format**
   - Error: Not a .csv file
   - Fix: Convert to CSV format

3. ✅ **FASTQ files exist**
   - Error: Paths in sample sheet don't exist
   - Fix: Update paths or check file locations

4. ✅ **Index file exists**
   - Error: Index not found
   - Fix: Build index or verify path

5. ✅ **Index extension**
   - Error: File is not `.idx`
   - Fix: Build proper Kallisto index

6. ✅ **Threads range**
   - Error: Threads < 1 or > 128
   - Fix: Use 1-128 threads

7. ✅ **Fragment length (single-end)**
   - Error: Not provided or out of range
   - Fix: Add `--fragment-length 200`

8. ✅ **Fragment SD (single-end)**
   - Error: Not provided or out of range
   - Fix: Add `--fragment-sd 20`

9. ✅ **Kallisto installed**
   - Error: Tool not in PATH
   - Fix: `conda install -c bioconda kallisto`

10. ✅ **Bootstraps range**
    - Error: Negative or too large
    - Fix: Use 0-1000

#### Salmon Pipeline (8-10 Checks)

1. ✅ **Sample sheet validation** (same as Kallisto)
2. ✅ **FASTQ files exist** (same as Kallisto)
3. ✅ **Index directory exists**
4. ✅ **Index validity** (contains `versionInfo.json`)
5. ✅ **Library type** (one of 13 valid types)
6. ✅ **Threads range** (1-128)
7. ✅ **Fragment length** (50-1000 bp)
8. ✅ **Fragment SD** (10-500 bp)
9. ✅ **Salmon installed**
10. ✅ **Memory allocation** (4-256 GB)

### Validation Error Messages

**Example: Single-End Fragment Length Error (Kallisto)**

```
❌ ERROR: Fragment length required for single-end reads

CRITICAL: Kallisto CANNOT estimate fragment length from single-end data!

You MUST provide both:
  --fragment-length <mean>  (e.g., 200)
  --fragment-sd <sd>        (e.g., 20)

How to estimate:
  • Bioanalyzer/TapeStation: Direct measurement
  • Picard CollectInsertSizeMetrics: From paired-end BAM
  • Default: 200 ± 20 bp (typical Illumina)

Solution:
  raptor quick-count -m kallisto \
      -s samples.csv \
      -i index.idx \
      --fragment-length 200 \
      --fragment-sd 20
```

**Example: Invalid Index Error**

```
❌ ERROR: Kallisto index must have .idx extension

Provided file: transcripts.fa
Expected: transcripts.idx

Build a new index:
  kallisto index -i transcripts.idx transcripts.fa.gz

Or verify you provided the correct file path.
```

---

## Output Files

### Primary Outputs

#### 1. `quick_gene_counts.csv`

**Purpose**: Gene-level count matrix for Modules 2-4

**Format**:
```csv
gene_id,Sample1,Sample2,Sample3,Sample4
ENSG00000223972,15,12,18,14
ENSG00000227232,245,198,267,221
ENSG00000243485,0,0,1,0
```

**Characteristics**:
- Rows: Genes (typically 20,000-60,000)
- Columns: Samples
- Values: Integer counts (rounded from estimated counts)
- Index: Gene IDs

**Usage**:
```bash
# Module 2: Quality Control
raptor qc --counts results/quick_counts/quick_gene_counts.csv

# Module 3: Data Profiling
raptor profile --counts results/quick_counts/quick_gene_counts.csv
```

#### 2. `quick_tpm.csv`

**Purpose**: TPM-normalized expression for visualization

**Format**: Same as counts, but with TPM values (floats)

**When to Use**:
- Cross-sample comparisons
- Heatmaps and clustering
- Expression level filtering

**When NOT to Use**:
- Differential expression (use raw counts)
- Statistical modeling (use raw counts)

#### 3. `sample_info.csv`

**Purpose**: Sample metadata enriched with QC metrics

**Format**:
```csv
sample_id,condition,batch,fastq_r1,fastq_r2,n_processed,n_mapped,mapping_rate,library_size
Sample1,Control,,/path/R1.fq.gz,/path/R2.fq.gz,25000000,23500000,94.0,23234567
Sample2,Control,,/path/R1.fq.gz,/path/R2.fq.gz,28000000,26320000,94.0,26109876
```

**Additional Columns Added**:
- `n_processed`: Total reads processed
- `n_mapped`: Successfully mapped reads
- `mapping_rate`: Percentage mapped
- `library_size`: Total gene-level counts

#### 4. `qc_report.html`

**Purpose**: Visual QC summary

**Contents**:
- Summary statistics cards
- Per-sample mapping rates
- Library size distribution
- Quality flags (Good/Warning/Poor)
- Next steps in RAPTOR workflow

**View**:
```bash
# Open in browser
open results/quick_counts/qc_report.html
```

### Intermediate Files (Optional)

**Location**: `results/quick_counts/{kallisto|salmon}_quant/`

**Per-Sample Directories**:
```
kallisto_quant/
├── Sample1/
│   ├── abundance.tsv       # Transcript-level counts
│   └── run_info.json       # QC metrics
├── Sample2/
│   ├── abundance.tsv
│   └── run_info.json
```

**Retention**:
```yaml
output:
  keep_quant_files: true   # Keep for debugging
  keep_quant_files: false  # Delete to save space (default)
```

**Disk Space**:
- ~100-500 MB per sample
- For 50 samples: 5-25 GB
- Recommendation: Delete after confirming success

---

## Usage Examples

### Example 1: Basic Paired-End (Kallisto)

```bash
raptor quick-count -m kallisto \
    -s samples.csv \
    -i /refs/kallisto_index.idx \
    -g /refs/tx2gene.csv \
    -o results/quick_counts
```

**Timeline**:
- Validation: 5 seconds
- Quantification: 3-5 minutes (20 samples)
- Aggregation: 10 seconds
- Total: ~5 minutes

### Example 2: Basic Paired-End (Salmon)

```bash
raptor quick-count -m salmon \
    -s samples.csv \
    -i /refs/salmon_index/ \
    -g /refs/tx2gene.csv \
    -o results/quick_counts
```

**Timeline**:
- Validation: 5 seconds
- Quantification: 5-8 minutes (20 samples)
- Aggregation: 10 seconds
- Total: ~8 minutes

### Example 3: Single-End with Kallisto

```bash
# CRITICAL: Must provide fragment length!
raptor quick-count -m kallisto \
    -s samples_single_end.csv \
    -i /refs/kallisto_index.idx \
    -g /refs/tx2gene.csv \
    --fragment-length 200 \
    --fragment-sd 20
```

**Fragment Length Estimation**:
```bash
# From Bioanalyzer CSV
fragment_length=$(awk -F',' 'NR>1 {sum+=$2; n++} END {print int(sum/n)}' bioanalyzer.csv)

# Manual specification
raptor quick-count -m kallisto ... \
    --fragment-length $fragment_length \
    --fragment-sd 20
```

### Example 4: Auto-Detect from Directory

```bash
# Step 1: Auto-detect FASTQ files
python3 detect_samples.py \
    --input /data/fastq_files/ \
    --output auto_samples.csv

# Step 2: Review and edit (add conditions)
nano auto_samples.csv

# Step 3: Run quantification
raptor quick-count -m salmon \
    -s auto_samples.csv \
    -i /refs/salmon_index/ \
    -g /refs/tx2gene.csv
```

### Example 5: Using Configuration File

```bash
# Custom config
cat > my_config.yaml << EOF
salmon:
  lib_type: "ISR"  # dUTP stranded
  gc_bias: true
  validate_mappings: true
  
resources:
  threads: 16
  memory_gb: 32
EOF

# Run with config
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/ \
    -g tx2gene.csv \
    --config my_config.yaml
```

### Example 6: Large Cohort (100+ samples)

```bash
# Optimize for throughput
raptor quick-count -m kallisto \
    -s large_cohort.csv \
    -i kallisto_index.idx \
    -g tx2gene.csv \
    --threads 8 \
    --parallel-samples 8  # Run 8 samples at once
```

**Resource Planning**:
- 8 threads × 8 samples = 64 cores
- 8 GB × 8 samples = 64 GB RAM
- Timeline: ~20-30 minutes for 100 samples

### Example 7: Stranded RNA-seq

```bash
# Kallisto: rf-stranded (dUTP)
raptor quick-count -m kallisto \
    -s samples.csv \
    -i kallisto_index.idx \
    --strand rf-stranded

# Salmon: ISR (dUTP)
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/ \
    --lib-type ISR
```

**Library Type Reference**:
```
Kallisto          Salmon      Protocol
---------         ------      --------
(default)         A           Auto-detect
rf-stranded       ISR         dUTP (Illumina TruSeq stranded)
fr-stranded       ISF         Ligation-based (fr-firststrand)
```

### Example 8: Using Bash Script

```bash
# Direct script usage
cd raptor/pipelines/quick_kallisto/
./scripts/run_quick.sh \
    -s /data/samples.csv \
    -x /refs/kallisto_index.idx \
    -t /refs/tx2gene.csv \
    -p 16
```

---

## Critical Parameters

### Fragment Length (Single-End Only)

**Why Critical**:
- Kallisto cannot infer from single-end data
- Affects effective length calculations
- Impacts count accuracy

**How to Determine**:

1. **Bioanalyzer/TapeStation** (Best):
   ```
   Mean: 250 bp
   SD: 30 bp
   ```

2. **From Paired-End Pilot**:
   ```bash
   # Run Picard on paired-end sample
   java -jar picard.jar CollectInsertSizeMetrics \
       I=sample.bam \
       O=metrics.txt \
       H=histogram.pdf
   
   # Extract median
   grep "MEDIAN_INSERT_SIZE" metrics.txt
   ```

3. **Literature/Protocol** (Fallback):
   - Standard TruSeq: 200 ± 20 bp
   - Small RNA: 140 ± 15 bp
   - Long-fragment: 300 ± 50 bp

**Impact of Errors**:
```
True: 200 bp
Specified: 250 bp
Impact: ~10-15% count bias for short genes

True: 200 bp  
Specified: 150 bp
Impact: ~20-25% count bias for long genes
```

**Recommendation**:
- If uncertain, use 200 ± 20 (conservative default)
- Better to be approximately right than precisely wrong
- Module 1 tolerates some imprecision (QC purposes)

### Thread Allocation

**Optimal Settings**:

```bash
# Standard workflow (balanced)
--threads 8

# Speed priority (diminishing returns)
--threads 16

# Memory-constrained
--threads 4
```

**Performance Scaling**:
```
Threads    Time (20 samples)    Efficiency
-------    -----------------    ----------
1          40 min               100%
4          12 min               83%
8          6 min                83%
16         4 min                62%
32         3 min                42%
```

**Best Practice**:
- Use 8 threads for most workflows
- Use 4 threads if memory limited
- Use 16 threads only for very large samples

### Library Type (Salmon)

**Auto-Detect (Recommended)**:
```bash
--lib-type A  # Let Salmon figure it out
```

**When to Specify**:
```bash
# dUTP stranded (most common)
--lib-type ISR

# You're certain and want to enforce
--lib-type ISF

# Error if wrong type detected
```

**Valid Types**:
```
Type    Description                         Use Case
----    -----------                         --------
A       Auto-detect                         Default (recommended)
IU      Inward, unstranded                  Standard unstranded
ISR     Inward, read1 sense                 dUTP (TruSeq stranded)
ISF     Inward, read1 antisense             Ligation (fr-firststrand)
MU      Matching, unstranded                Mate-pair
MSR     Matching, read1 sense               Mate-pair stranded
MSF     Matching, read1 antisense           Mate-pair stranded
OU      Outward, unstranded                 Mate-pair
OSR     Outward, read1 sense                Mate-pair stranded
OSF     Outward, read1 antisense            Mate-pair stranded
SF      Read1 sense                         Single-end stranded
SR      Read1 antisense                     Single-end stranded
U       Unstranded                          Single-end unstranded
```

### GC Bias Correction

**Salmon Only** (Kallisto doesn't support):

```yaml
gc_bias: true   # Default for Module 1
```

**When to Enable**:
- ✅ Always (minimal performance cost)
- ✅ PCR-amplified libraries
- ✅ Samples with GC-rich/poor regions

**When to Disable**:
- ❌ When speed is absolutely critical
- ❌ Paired with `--numGibbsSamples` (incompatible)

**Performance Impact**:
- Adds ~10-15% runtime
- Improves count accuracy 2-5%

---

## Integration with RAPTOR Workflow

### The Complete Pipeline

```
Module 1: Quantify (THIS MODULE)
    ↓
    quick_gene_counts.csv
    ↓
Module 2: Quality Control
    ↓
    Outlier detection
    Batch effect visualization
    Sample filtering
    ↓
Module 3: Data Profiling  
    ↓
    32-feature characterization
    Complexity metrics
    Data type classification
    ↓
Module 4: Recommendation
    ↓
    Optimal pipeline selection
    Parameter suggestions
    ↓
Module 5: Production Quantification
    ↓
    High-quality counts for DE
```

### Handoff to Module 2

```bash
# After Module 1 completes
raptor qc --counts results/quick_counts/quick_gene_counts.csv

# What Module 2 does with quick counts:
# 1. PCA outlier detection
# 2. Sample correlation analysis
# 3. Batch effect visualization
# 4. Filtering recommendations
```

### Handoff to Module 3

```bash
# Profile your data
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# What Module 3 analyzes:
# - Sparsity (zeros in matrix)
# - Dynamic range
# - Mean-variance relationship
# - Gene detection patterns
# - 32 total features
```

### Handoff to Module 4

```bash
# Get pipeline recommendation
raptor recommend

# What Module 4 recommends based on Module 1 + 3:
# - Best DE tool (DESeq2, edgeR, limma-voom)
# - Whether to use blocking
# - Batch correction strategy
# - Filtering thresholds
```

### Why Not Use Quick Counts for DE?

**Technical Reasons**:

1. **No Bootstraps**
   - Can't estimate uncertainty
   - Tools like Sleuth need this

2. **Speed Optimizations**
   - Simplified bias modeling
   - Less rigorous mapping validation

3. **Intended Use**
   - Quick counts: "Is my data good?"
   - Production counts: "What's differentially expressed?"

**Practical Example**:

```bash
# ❌ WRONG: Using quick counts for publication
raptor quick-count -m kallisto ...
raptor de --counts results/quick_counts/quick_gene_counts.csv

# ✅ CORRECT: Use Module 5 for DE
raptor quick-count -m kallisto ...  # Explore data
raptor recommend                     # Get optimized params
raptor production --pipeline star    # Final quantification
raptor de --counts results/production/counts.csv
```

### File Naming Convention

**Critical Design Decision**:

```
Quick Counts:       quick_gene_counts.csv
Production Counts:  production_counts.csv
                    gene_counts.csv
                    counts_final.csv
```

**Why Different Names**:
- Prevents accidental use of quick counts in DE
- Makes file purpose immediately clear
- Enables both to coexist in same directory

---

## Troubleshooting

### Common Errors

#### 1. Fragment Length Error (Kallisto, Single-End)

**Error Message**:
```
❌ ERROR: Fragment length required for single-end reads
CRITICAL: Kallisto CANNOT estimate fragment length!
```

**Cause**: Single-end data without fragment length parameters

**Solution**:
```bash
# Add required parameters
raptor quick-count -m kallisto \
    -s samples.csv \
    -i index.idx \
    --fragment-length 200 \
    --fragment-sd 20
```

**How to Estimate**:
- Check Bioanalyzer report: "Mean" and "Std Dev"
- Use previous paired-end run with same library prep
- Default: 200 ± 20 bp for standard Illumina

#### 2. Low Mapping Rate

**Symptom**: QC report shows <50% mapping rate

**Possible Causes**:

1. **Wrong Index Organism**:
   ```
   Data: Mouse samples
   Index: Human transcriptome
   Result: ~5-10% mapping
   ```
   
   **Fix**: Rebuild index with correct organism

2. **Poor Quality Reads**:
   ```bash
   # Check FASTQ quality
   fastqc sample_R1.fastq.gz
   
   # Look for:
   # - Per base quality < 20
   # - Adapter contamination
   # - Over-represented sequences
   ```
   
   **Fix**: Trim adapters and low-quality bases

3. **Contamination**:
   ```bash
   # Screen for contaminants
   fastq_screen --aligner bowtie2 sample_R1.fastq.gz
   ```
   
   **Fix**: Filter out contaminant reads

#### 3. Index Not Found

**Error Message**:
```
❌ ERROR: Kallisto index not found: transcripts.idx
```

**Causes**:
- Typo in path
- File moved or deleted
- Wrong file type provided

**Solution**:
```bash
# Verify file exists
ls -lh /refs/kallisto_index.idx

# Check it's actually a Kallisto index
file /refs/kallisto_index.idx
# Should show: data

# Rebuild if necessary
kallisto index -i kallisto_index.idx transcripts.fa.gz
```

#### 4. Out of Memory

**Symptom**:
```
Killed
-bash: line 1: 12345 Killed  kallisto quant ...
```

**Cause**: Insufficient RAM

**Solutions**:

1. **Reduce Threads**:
   ```bash
   # Fewer threads = less memory per sample
   --threads 4  # Instead of 8
   ```

2. **Switch to Kallisto**:
   ```bash
   # Kallisto uses less memory than Salmon
   raptor quick-count -m kallisto ...
   ```

3. **Process Sequentially**:
   ```yaml
   resources:
     parallel_samples: false  # One at a time
   ```

4. **Use HPC**:
   ```bash
   # Request more memory
   #SBATCH --mem=32G
   ```

#### 5. FASTQ Files Not Found

**Error Message**:
```
❌ ERROR: FASTQ file not found: /path/to/Sample1_R1.fastq.gz
```

**Causes**:
- Absolute vs relative paths
- Files moved after sample sheet creation
- Network drive unmounted

**Solution**:
```bash
# Use absolute paths in sample sheet
realpath sample_R1.fastq.gz

# Or create symbolic links
ln -s /mnt/data/fastq/*.fastq.gz ./
```

#### 6. Tx2Gene Mapping Issues

**Symptom**: Many unmapped transcripts

**Error in Log**:
```
Warning: 15,234 transcripts not found in tx2gene mapping
```

**Causes**:
1. Index and tx2gene from different sources
2. Version mismatch (e.g., Ensembl 109 vs 110)
3. Wrong organism

**Solution**:
```bash
# Verify compatibility
head -5 tx2gene.csv
# ENST00000456328,ENSG00000223972

zgrep ">" transcripts.fa.gz | head -5
# >ENST00000456328.2 ...

# Transcript IDs should match (ignore version numbers)
```

**Fix**:
```bash
# Extract tx2gene from same GTF used for index
python3 - << 'EOF'
import re
import sys

tx2gene = {}
with open('annotation.gtf', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] == 'transcript':
            tx_match = re.search(r'transcript_id "([^"]+)"', fields[8])
            gene_match = re.search(r'gene_id "([^"]+)"', fields[8])
            if tx_match and gene_match:
                tx = tx_match.group(1).split('.')[0]  # Remove version
                gene = gene_match.group(1).split('.')[0]
                tx2gene[tx] = gene

with open('tx2gene.csv', 'w') as out:
    out.write('transcript_id,gene_id\n')
    for tx, gene in tx2gene.items():
        out.write(f'{tx},{gene}\n')
EOF
```

#### 7. Pipeline Hangs

**Symptom**: No progress for >30 minutes

**Causes**:
1. Network FASTQ files (slow read)
2. Large file decompression
3. Insufficient memory causing swapping

**Diagnosis**:
```bash
# Check if process is running
ps aux | grep kallisto

# Check I/O wait
top
# Look for high %wa (I/O wait)

# Check memory
free -h
# Look for low available memory
```

**Solutions**:
```bash
# Copy FASTQ locally first
cp /network/path/*.fastq.gz /local/tmp/
raptor quick-count ... -s samples_local.csv

# Or use faster local storage
raptor quick-count ... -o /tmp/results/
```

### Performance Optimization

#### Baseline Performance

**Kallisto** (8 threads, human transcriptome):
```
Sample size    Time per sample
-----------    ---------------
10M reads      30 seconds
50M reads      2 minutes
100M reads     4 minutes
```

**Salmon** (8 threads, human transcriptome):
```
Sample size    Time per sample
-----------    ---------------
10M reads      1 minute
50M reads      4 minutes
100M reads     8 minutes
```

#### Optimization Strategies

1. **Use Kallisto for Large Cohorts**:
   ```bash
   # 100 samples, 50M reads each
   # Kallisto: ~3 hours
   # Salmon: ~7 hours
   raptor quick-count -m kallisto ...
   ```

2. **Parallel Processing**:
   ```bash
   # If you have 64 cores
   raptor quick-count -m kallisto \
       --threads 8 \
       --parallel-samples 8
   ```

3. **Local Storage**:
   ```bash
   # Copy to NVMe/SSD
   cp -r /slow/storage/fastq/ /fast/nvme/
   raptor quick-count ... -i /fast/nvme/fastq/
   ```

4. **Predecompress FASTQ** (if running multiple times):
   ```bash
   # Decompress once
   gunzip *.fastq.gz
   
   # Faster subsequent runs
   raptor quick-count -m kallisto -s samples.csv ...
   ```

---

## Best Practices

### 1. Always Use tx2gene Mapping

❌ **Bad**:
```bash
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/
# Output: transcript-level counts (not useful for Modules 2-4)
```

✅ **Good**:
```bash
raptor quick-count -m salmon \
    -s samples.csv \
    -i salmon_index/ \
    -g tx2gene.csv
# Output: gene-level counts (works with Modules 2-4)
```

### 2. Verify Sample Sheet Before Running

```bash
# Check format
head -5 samples.csv

# Verify files exist
awk -F',' 'NR>1 {print $4}' samples.csv | while read f; do
    [ -f "$f" ] || echo "Missing: $f"
done

# Check for duplicates
awk -F',' 'NR>1 {print $1}' samples.csv | sort | uniq -d
```

### 3. Use Consistent Naming

```bash
# ✅ Good: Systematic naming
Sample001_Control_Rep1
Sample002_Control_Rep2
Sample003_Treatment_Rep1

# ❌ Bad: Inconsistent naming  
John_experiment_1
sample2_FINAL
Sample_03_v2_updated
```

### 4. Document Fragment Length Estimates

```bash
# Create a README with your experiment
cat > experiment_notes.md << EOF
# RNA-seq Experiment Notes

## Library Prep
- Protocol: TruSeq Stranded mRNA
- Fragment selection: 200-300 bp
- Bioanalyzer: Mean 250 bp, SD 25 bp

## Quantification Parameters
- Fragment length: 250
- Fragment SD: 25
- Tool: Kallisto v0.48.0

## Date: 2024-03-09
EOF
```

### 5. Save QC Report

```bash
# Archive QC reports with date
cp results/quick_counts/qc_report.html \
   results/qc_reports/qc_$(date +%Y%m%d).html
```

### 6. Use Configuration Files

```bash
# Create reusable configs
cat > config_standard.yaml << EOF
salmon:
  lib_type: "A"
  gc_bias: true
  validate_mappings: true
  
resources:
  threads: 8
  memory_gb: 16
EOF

# Reference in runs
raptor quick-count -m salmon \
    --config config_standard.yaml \
    -s samples.csv -i index/
```

### 7. Test with Subset First

```bash
# Test with 2-3 samples first
head -4 samples.csv > samples_test.csv

raptor quick-count -m salmon \
    -s samples_test.csv \
    -i index/ -g tx2gene.csv

# Verify success, then run full cohort
raptor quick-count -m salmon \
    -s samples.csv \
    -i index/ -g tx2gene.csv
```

### 8. Monitor Resource Usage

```bash
# During run, check resources
htop  # Interactive
top   # Simple

# Or log usage
while true; do
    date >> resource_usage.log
    free -h >> resource_usage.log
    ps aux | grep kallisto >> resource_usage.log
    sleep 60
done &
```

---

## Technical Details

### Pseudo-Alignment vs Traditional Alignment

#### Traditional Alignment (STAR, HISAT2)

**Process**:
1. Build full genome index (30-50 GB)
2. Align reads to genome positions
3. Resolve splice junctions
4. Generate BAM files (large)
5. Count reads per gene

**Pros**:
- ✅ Precise genomic locations
- ✅ Can detect novel transcripts
- ✅ Enables variant calling

**Cons**:
- ❌ Slow (1-3 hours per sample)
- ❌ Large disk usage (50-100 GB per sample)
- ❌ Complex processing

#### Pseudo-Alignment (Kallisto, Salmon)

**Process**:
1. Build transcriptome index (2-8 GB)
2. Match k-mers to transcript sequences
3. Assign reads probabilistically
4. Estimate counts directly

**Pros**:
- ✅ Fast (2-10 minutes per sample)
- ✅ Small disk usage (500 MB per sample)
- ✅ Direct quantification

**Cons**:
- ❌ No genomic coordinates
- ❌ Can't detect novel transcripts
- ❌ No variant calling

**Module 1 Choice**: Pseudo-alignment
- Speed enables rapid iteration
- Sufficient accuracy for QC/profiling
- Genomic coordinates not needed for count matrices

### Count Aggregation Method

**Transcript → Gene Aggregation**:

```python
# Simplified version of what Module 1 does

# 1. Load transcript counts
transcript_counts = {
    'ENST00000456328': 15.4,
    'ENST00000450305': 8.2,
    'ENST00000488147': 22.1
}

# 2. Map to genes
tx2gene = {
    'ENST00000456328': 'ENSG00000223972',
    'ENST00000450305': 'ENSG00000223972',
    'ENST00000488147': 'ENSG00000227232'
}

# 3. Aggregate by summing
gene_counts = {}
for tx_id, count in transcript_counts.items():
    gene_id = tx2gene[tx_id]
    if gene_id not in gene_counts:
        gene_counts[gene_id] = 0
    gene_counts[gene_id] += count

# 4. Round to integers
gene_counts = {g: round(c) for g, c in gene_counts.items()}

# Result:
# ENSG00000223972: 24  (15.4 + 8.2)
# ENSG00000227232: 22  (22.1)
```

**Why Summing**:
- ✅ Preserves total count mass
- ✅ Standard approach (tximport, DESeq2)
- ✅ Compatible with downstream tools

**Alternatives** (not used):
- ❌ Averaging: Loses count information
- ❌ Median: Biased towards low-expressed isoforms
- ❌ Maximum: Ignores minor isoforms

### TPM Calculation

**Transcripts Per Million (TPM)**:

```python
# TPM calculation (simplified)

# 1. Normalize by feature length
rpm_per_kb = count / (effective_length / 1000)

# 2. Scale to per-million
tpm = rpm_per_kb / (sum(rpm_per_kb) / 1e6)

# Example:
# Gene A: 1000 counts, 2000 bp effective length
# Gene B: 500 counts, 1000 bp effective length
# 
# Gene A RPK: 1000 / 2 = 500
# Gene B RPK: 500 / 1 = 500
# Sum RPK: 1000
# 
# Gene A TPM: 500 / (1000/1e6) = 500,000
# Gene B TPM: 500 / (1000/1e6) = 500,000
```

**Properties**:
- ✅ Accounts for gene length
- ✅ Comparable across samples
- ✅ Sums to 1 million per sample

**Use Cases**:
- ✅ Heatmaps
- ✅ Cross-sample visualization
- ✅ Expression filtering

**NOT for**:
- ❌ Differential expression (use raw counts)
- ❌ Statistical testing (use raw counts)

### Effective Length

**Definition**: Expected length of fragments from a transcript

**Calculation** (Kallisto):
```
effective_length = transcript_length - mean_fragment_length + 1
```

**Why It Matters**:
- Longer transcripts generate more fragments
- Must normalize counts for length
- Affects TPM calculation

**Example**:
```
Transcript: 3000 bp
Fragment length: 200 bp
Effective length: 3000 - 200 + 1 = 2801 bp
```

**Gene-Level Aggregation**:
```python
# Module 1 uses mean effective length across transcripts
gene_eff_length = mean([
    tx1_eff_length,
    tx2_eff_length,
    tx3_eff_length
])
```

---

## Version History

### 2.2.0 (Current)

**Major Changes**:
- ✅ Comprehensive validation system (10 checks)
- ✅ Enhanced error messages with solutions
- ✅ Architecture-compliant file naming
- ✅ Validation schemas in configuration
- ✅ Improved documentation

### 2.1.0

- ✅ Dual implementation (Kallisto + Salmon)
- ✅ Auto-detection of FASTQ files
- ✅ HTML QC report generation

### 2.0.0

- ✅ Initial modular design
- ✅ Separation from production pipelines

---

## Summary

### Key Takeaways

1. **Purpose**: Module 1 is for rapid QC and profiling, NOT final DE analysis

2. **Tool Choice**: 
   - Default: Salmon (better accuracy)
   - Alternative: Kallisto (maximum speed)

3. **Critical**: Single-end Kallisto requires fragment length parameters

4. **Output**: `quick_gene_counts.csv` feeds into Modules 2-4

5. **Validation**: v2.2.0 catches errors before wasting time

6. **Integration**: Designed to inform Module 5 production quantification

### Workflow Reminder

```bash
# 1. Quick quantification (Module 1) - THIS MODULE
raptor quick-count -m salmon -s samples.csv -i index/ -g tx2gene.csv

# 2. Quality control (Module 2)
raptor qc --counts results/quick_counts/quick_gene_counts.csv

# 3. Data profiling (Module 3)
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# 4. Get recommendation (Module 4)
raptor recommend

# 5. Production quantification (Module 5)
raptor production --pipeline <recommended>

# 6. Differential expression
raptor de --counts results/production/counts.csv
```

### Getting Help

**Documentation**:
- This file: Module 1 comprehensive guide
- `README.md`: Quick reference
- `config.yaml`: Parameter explanations

**Support**:
- GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
- Email: ayehbolouki1988@gmail.com

**Citation**:
```
Bolouki, A. (2024). RAPTOR: RNA-seq Analysis Pipeline Testing and 
Optimization Resource. GitHub. https://github.com/AyehBlk/RAPTOR
```

---

**Document Version**: 2.2.0  
**Last Updated**: March 2026  
**Author**: Ayeh Bolouki  
**License**: MIT
