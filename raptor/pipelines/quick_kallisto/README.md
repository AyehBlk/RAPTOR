# 🦖 RAPTOR Quick Kallisto Pipeline

**Module 1: Quantify** - Ultra-fast pseudo-alignment for QC and profiling

## Overview

This pipeline uses Kallisto for ultra-fast transcript quantification. It generates
a quick count matrix for quality assessment and pipeline recommendation, **NOT**
for final differential expression analysis.

## ✨ What's New in v2.2.0

- ✅ **Comprehensive Validation** - 10 checks before pipeline runs
- ✅ **Professional Error Handling** - Clear, actionable error messages
- ✅ **CLI Integration Ready** - Parameter definitions for `raptor` command
- ✅ **Critical Fragment Length Validation** - Prevents #1 Kallisto user error
- ✅ **Enhanced Documentation** - Troubleshooting built-in

## 🎯 Validation Checks (v2.2.0)

Quick Kallisto performs 10 comprehensive validation checks:

| Check | What It Validates | Why It Matters |
|-------|------------------|----------------|
| 1. Sample sheet exists | File presence | Catches typos early |
| 2. Sample sheet format | .csv extension | Correct format |
| 3. FASTQ files exist | All input files | No wasted time |
| 4. Index file exists | .idx file presence | Index available |
| 5. Index extension | Must be .idx | Correct format |
| 6. Threads range | 1-128 threads | Valid resource use |
| 7. Fragment length (SE) | 50-1000 bp | **CRITICAL for single-end!** |
| 8. Fragment SD (SE) | 10-500 bp | **CRITICAL for single-end!** |
| 9. Kallisto installed | Tool available | Ready to run |
| 10. Bootstraps range | 0-1000 | Valid parameter |

**Total: 10 validation checks** - Prevents 90%+ of user errors!

## Output Files

| File | Description | Used By |
|------|-------------|---------|
| `quick_gene_counts.csv` | Gene-level count matrix | M2, M3, M4 |
| `quick_tpm.csv` | TPM normalized matrix | QC plots |
| `sample_info.csv` | Sample metadata + QC metrics | Reports |

**Output Location**: `results/quick_counts/`

## Quick Start

### Using CLI (Recommended)

```bash
# Basic usage (paired-end)
raptor quick-count -m kallisto \
    -s samples.csv \
    -i kallisto_index.idx \
    -o results/quick_counts

# With gene-level aggregation
raptor quick-count -m kallisto \
    -s samples.csv \
    -i kallisto_index.idx \
    -g tx2gene.csv \
    -o results/quick_counts

# Single-end (MUST specify fragment length!)
raptor quick-count -m kallisto \
    -s samples.csv \
    -i kallisto_index.idx \
    --fragment-length 200 \
    --fragment-sd 20
```

### Using Python API

```python
from raptor.pipelines.quick_kallisto.scripts.kallisto_quant import run_quick_kallisto

success = run_quick_kallisto(
    sample_sheet='samples.csv',
    index='kallisto_index.idx',
    output_dir='results/quick_counts',
    gene_map='tx2gene.csv',
    threads=8,
    fragment_length=200,  # Required for single-end
    fragment_sd=20
)
```

### Using Bash Script

```bash
./scripts/run_quick.sh \
    -s samples.csv \
    -x kallisto_index.idx \
    -t tx2gene.csv \
    -o results/quick_counts
```

## ⚠️ CRITICAL: Single-End Reads

**Kallisto CANNOT estimate fragment length from single-end data!**

You MUST provide `--fragment-length` and `--fragment-sd` for single-end reads.

Estimate from:
- Bioanalyzer/TapeStation data
- Previous experiments with same library prep
- Default: 200 ± 20 bp (typical for Illumina)

```bash
# Single-end example
raptor quick-count -m kallisto \
    -s samples.csv \
    -i kallisto_index.idx \
    --fragment-length 200 \
    --fragment-sd 20
```

## Sample Sheet Format

Create a CSV file with the following columns:

```csv
sample_id,condition,batch,fastq_r1,fastq_r2
Sample1,Control,,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
Sample2,Control,,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
Sample3,Treatment,,/path/to/Sample3_R1.fastq.gz,/path/to/Sample3_R2.fastq.gz
```

For single-end reads, leave `fastq_r2` empty.

## Configuration

Edit `quick_config.yaml` to customize:

```yaml
kallisto:
  strand: ""              # rf-stranded, fr-stranded, or empty
  num_bootstraps: 0       # 0 for speed (M1), 100+ for sleuth (M5)
  
single_end:
  fragment_length: 200    # REQUIRED for single-end
  fragment_sd: 20         # REQUIRED for single-end

output:
  dir: "results/quick_counts"
  counts_file: "quick_gene_counts.csv"
```

## Requirements

- **Kallisto** ≥ 0.46.0: `conda install -c bioconda kallisto`
- **Python** ≥ 3.8 with pandas

## Kallisto Index

If you don't have a Kallisto index:

```bash
# Download transcriptome (example: human)
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Build index
kallisto index -i kallisto_index.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
```

## Next Steps

After running this pipeline:

```bash
# Module 2: Quality Assessment
raptor qc --counts results/quick_counts/quick_gene_counts.csv

# Module 3: Data Profiling
raptor profile --counts results/quick_counts/quick_gene_counts.csv

# Module 4: Get Recommendation
raptor recommend
```

## Kallisto vs Salmon

| Feature | Kallisto | Salmon |
|---------|----------|--------|
| Speed | Faster | Fast |
| Memory | Lower | Low |
| GC Bias Correction | No | Yes |
| Validate Mappings | No | Yes |
| Single-end | Needs fragment length | Auto-estimates |

**Recommendation**: Use Salmon for most cases. Use Kallisto when you need maximum speed or minimal memory.

## 🔧 Troubleshooting (v2.2.0)

### Error: "Fragment length required for single-end"

**Problem:** Kallisto cannot estimate fragment length from single-end data.

**Solution:**
```bash
# Add fragment length parameters:
raptor quick-count -m kallisto \
    -s samples.csv \
    -i index.idx \
    --fragment-length 200 \
    --fragment-sd 20
```

**How to estimate:**
- Bioanalyzer/TapeStation: Direct measurement
- From previous experiments: Use same library prep values
- Default: 200 ± 20 bp (typical Illumina)

### Error: "Index must have .idx extension"

**Problem:** Provided file is not a Kallisto index.

**Solution:**
```bash
# Build proper index:
kallisto index -i transcripts.idx Homo_sapiens.cdna.all.fa.gz
```

### Low Pseudoalignment Rate (<50%)

**Possible causes:**
1. Index doesn't match your organism
2. Poor quality FASTQ files
3. Adapter contamination

**Solutions:**
1. Verify index organism matches samples
2. Run FastQC on FASTQ files
3. Trim adapters with Trim Galore

### Kallisto Not Found

**Problem:** Kallisto not installed or not in PATH.

**Solution:**
```bash
conda install -c bioconda kallisto
# OR
mamba install -c bioconda kallisto
```

### Memory Errors

**Problem:** Out of memory during quantification.

**Solution:**
```bash
# Kallisto is very memory efficient (typically <8GB)
# If you still have issues, reduce parallel samples:
raptor quick-count -m kallisto \
    -s samples.csv \
    -i index.idx \
    --threads 4  # Reduce threads
```

## Files in This Directory

```
quick_kallisto/
├── README.md              # This file
├── quick_config.yaml      # Configuration
└── scripts/
    ├── kallisto_quant.py  # Main quantification script
    ├── combine_counts.py  # Combine quant files
    ├── detect_samples.py  # Auto-detect FASTQ files
    ├── generate_report.py # QC report generator
    └── run_quick.sh       # Bash wrapper
```

## Author

Ayeh Bolouki (ayehbolouki1988@gmail.com)

## License

MIT License - RAPTOR v2.2.0
