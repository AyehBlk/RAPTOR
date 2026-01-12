# Quick Kallisto Pipeline

## Overview

**Ultra-fast pseudo-alignment quantification for QC and profiling purposes.**

This pipeline runs Kallisto quantification on all samples and generates a count matrix that can be used for:
- Quality Assessment (outlier detection, batch effects)
- Data Profiling (BCV, library size stats)
- ML Pipeline Recommendation

### Output Files

```
output/
├── counts.csv         # Gene/transcript count matrix
├── tpm.csv            # TPM normalized matrix
├── sample_info.csv    # Sample QC metrics (pseudoalignment rate, etc.)
└── quick_kallisto.log # Pipeline log
```

---

## Quick Start

### 1. Prepare Sample Sheet

Create a CSV file with your samples:

**Paired-end:**
```csv
sample_id,condition,fastq_r1,fastq_r2
Sample1,Control,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
Sample2,Control,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
Sample3,Treatment,/path/to/Sample3_R1.fastq.gz,/path/to/Sample3_R2.fastq.gz
Sample4,Treatment,/path/to/Sample4_R1.fastq.gz,/path/to/Sample4_R2.fastq.gz
```

**Single-end:**
```csv
sample_id,condition,fastq
Sample1,Control,/path/to/Sample1.fastq.gz
Sample2,Control,/path/to/Sample2.fastq.gz
Sample3,Treatment,/path/to/Sample3.fastq.gz
Sample4,Treatment,/path/to/Sample4.fastq.gz
```

### 2. Run Pipeline

**Via RAPTOR CLI (recommended):**
```bash
# Paired-end
raptor quick-count \
  --method kallisto \
  --sample-sheet samples.csv \
  --index /path/to/kallisto_index.idx \
  --output quick_counts/ \
  --threads 8

# Single-end (MUST specify fragment length!)
raptor quick-count \
  --method kallisto \
  --sample-sheet samples.csv \
  --index /path/to/kallisto_index.idx \
  --output quick_counts/ \
  --fragment-length 200 \
  --fragment-sd 20 \
  --threads 8
```

**Via Python:**
```python
from quick_kallisto.scripts.kallisto_quant import run_quick_kallisto

run_quick_kallisto(
    sample_sheet='samples.csv',
    index='/path/to/kallisto_index.idx',
    output_dir='quick_counts/',
    threads=8,
    fragment_length=200,  # Required for single-end
    fragment_sd=20        # Required for single-end
)
```

### 3. Next Steps

After generating counts, proceed with:

```bash
# Quality Assessment
raptor qc --counts quick_counts/counts.csv --output qc_report/

# Profile and get ML recommendation
raptor profile --counts quick_counts/counts.csv --use-ml
```

---

## ⚠️ Important: Single-End Reads

**Kallisto REQUIRES fragment length parameters for single-end data!**

Unlike Salmon, Kallisto cannot estimate fragment length from single-end reads. You MUST provide:
- `--fragment-length`: Mean fragment length (typically 150-300)
- `--fragment-sd`: Standard deviation (typically 20-50)

**How to estimate fragment length:**
1. **Bioanalyzer/TapeStation**: Check your library QC report
2. **Previous runs**: Look at insert size distribution
3. **Defaults**: Use 200 ± 20 for standard libraries

```bash
# Single-end with fragment length
raptor quick-count \
  --method kallisto \
  --sample-sheet samples.csv \
  --index kallisto_index.idx \
  --fragment-length 200 \
  --fragment-sd 20 \
  --output counts/
```

---

## Kallisto Index

You need a Kallisto index before running. Create one with:

```bash
# Download transcriptome
# Human: https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/
# Mouse: https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/

# Build index
kallisto index -i kallisto_index.idx transcriptome.fa
```

Or use pre-built indices from [refgenie](http://refgenie.databio.org/).

---

## Gene-Level Counts

By default, Kallisto outputs transcript-level counts. For gene-level counts, provide a tx2gene mapping:

```bash
raptor quick-count \
  --method kallisto \
  --sample-sheet samples.csv \
  --index kallisto_index.idx \
  --gene-map tx2gene.csv \
  --output quick_counts/
```

**tx2gene.csv format:**
```csv
transcript_id,gene_id
ENST00000456328,ENSG00000223972
ENST00000450305,ENSG00000223972
ENST00000488147,ENSG00000227232
...
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--sample-sheet` | Sample sheet CSV file | Required |
| `--index` | Kallisto index file (.idx) | Required |
| `--output` | Output directory | Required |
| `--threads` | Number of threads | 8 |
| `--gene-map` | tx2gene mapping file | None |
| `--fragment-length` | Fragment length mean (single-end) | 200 |
| `--fragment-sd` | Fragment length SD (single-end) | 20 |
| `--keep-quant` | Keep individual quant files | False |

### Strand-Specific Options

| Option | Description |
|--------|-------------|
| (none) | Unstranded (default) |
| `--rf-stranded` | Reverse-forward stranded (e.g., TruSeq) |
| `--fr-stranded` | Forward-reverse stranded |

---

## Output Files

### counts.csv

Gene/transcript count matrix (estimated counts):
```csv
,Sample1,Sample2,Sample3,Sample4
ENSG00000223972,145,162,89,95
ENSG00000227232,523,498,612,589
ENSG00000243485,0,0,2,1
...
```

### tpm.csv

TPM normalized expression:
```csv
,Sample1,Sample2,Sample3,Sample4
ENSG00000223972,2.45,2.62,1.89,1.95
ENSG00000227232,8.23,7.98,9.12,8.89
...
```

### sample_info.csv

Sample QC metrics:
```csv
sample_id,success,num_reads,num_pseudoaligned,pseudoalign_rate,read_type,condition,batch
Sample1,True,25000000,21875000,87.5,paired,Control,
Sample2,True,28000000,24976000,89.2,paired,Control,
Sample3,True,22000000,18722000,85.1,paired,Treatment,
Sample4,True,26000000,22958000,88.3,paired,Treatment,
```

---

## Salmon vs Kallisto: When to Use Which?

| Aspect | Salmon | Kallisto |
|--------|--------|----------|
| **Speed** | Fast | Faster |
| **Memory** | Low-Medium | Low |
| **Accuracy** | Slightly better | Very good |
| **GC bias correction** | ✅ Built-in | ❌ No |
| **Single-end** | Easy (auto-estimates) | Requires fragment length |
| **Bootstrap** | Supported | Supported |

**Recommendation:**
- **Use Salmon** if you have single-end data (easier)
- **Use Kallisto** if you need maximum speed
- **Both work well** for QC purposes

---

## Troubleshooting

### "Kallisto index not found"

Ensure you've built the index:
```bash
ls -la /path/to/kallisto_index.idx
# Should be a single .idx file
```

### "Single-end requires fragment length"

For single-end data, you MUST specify:
```bash
raptor quick-count --method kallisto \
  --fragment-length 200 --fragment-sd 20 ...
```

### "Sample sheet validation failed"

Check that:
1. All FASTQ files exist at specified paths
2. File extensions are correct (.fastq, .fq, .fastq.gz, .fq.gz)
3. No duplicate sample IDs
4. For paired-end: both R1 and R2 columns are present

### Low pseudoalignment rate

Low rate (<50%) may indicate:
- Wrong reference transcriptome (species mismatch)
- Contamination
- Low quality reads

---

## Citation

If you use this pipeline, please cite:

**Kallisto:**
```
Bray, N.L., et al. (2016). Near-optimal probabilistic RNA-seq 
quantification. Nature Biotechnology, 34(5), 525-527.
```

**RAPTOR:**
```
Bolouki, A. (2025). RAPTOR: RNA-seq Analysis Pipeline Testing 
and Optimization Resource. https://github.com/AyehBlk/RAPTOR
```

---

**Author:** Ayeh Bolouki  
**Version:** 2.2.0  
**License:** MIT
