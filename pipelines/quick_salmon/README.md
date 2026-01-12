# Quick Salmon Pipeline

## Overview

**Fast pseudo-alignment quantification for QC and profiling purposes.**

This pipeline runs Salmon quantification on all samples and generates a count matrix that can be used for:
- Quality Assessment (outlier detection, batch effects)
- Data Profiling (BCV, library size stats)
- ML Pipeline Recommendation

### Output Files

```
output/
├── counts.csv         # Gene/transcript count matrix
├── tpm.csv            # TPM normalized matrix
├── sample_info.csv    # Sample QC metrics (mapping rate, etc.)
└── quick_salmon.log   # Pipeline log
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
raptor quick-count \
  --method salmon \
  --sample-sheet samples.csv \
  --index /path/to/salmon_index \
  --output quick_counts/ \
  --threads 8
```

**Via Python:**
```python
from quick_salmon.scripts.salmon_quant import run_quick_salmon

run_quick_salmon(
    sample_sheet='samples.csv',
    index='/path/to/salmon_index',
    output_dir='quick_counts/',
    threads=8
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

## Salmon Index

You need a Salmon index before running. Create one with:

```bash
# Download transcriptome
# Human: https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/
# Mouse: https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/

# Build index
salmon index -t transcriptome.fa -i salmon_index -k 31
```

Or use pre-built indices from [refgenie](http://refgenie.databio.org/).

---

## Gene-Level Counts

By default, Salmon outputs transcript-level counts. For gene-level counts, provide a tx2gene mapping:

```bash
raptor quick-count \
  --method salmon \
  --sample-sheet samples.csv \
  --index salmon_index \
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

Generate tx2gene from GTF:
```bash
# Using awk
awk -F'\t' '$3=="transcript" {
    match($9, /gene_id "([^"]+)"/, g);
    match($9, /transcript_id "([^"]+)"/, t);
    print t[1]","g[1]
}' annotation.gtf > tx2gene.csv
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--sample-sheet` | Sample sheet CSV file | Required |
| `--index` | Salmon index directory | Required |
| `--output` | Output directory | Required |
| `--threads` | Number of threads | 8 |
| `--gene-map` | tx2gene mapping file | None |
| `--lib-type` | Library type (A=auto) | A |
| `--keep-quant` | Keep individual quant files | False |

### Library Types

| Type | Description |
|------|-------------|
| `A` | Automatic detection (recommended) |
| `ISR` | Inward, stranded, reverse (TruSeq) |
| `ISF` | Inward, stranded, forward |
| `IU` | Inward, unstranded |
| `SR` | Stranded, reverse |
| `SF` | Stranded, forward |
| `U` | Unstranded |

---

## Output Files

### counts.csv

Gene/transcript count matrix:
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
sample_id,success,num_reads,mapping_rate,read_type,condition,batch
Sample1,True,25000000,87.5,paired,Control,
Sample2,True,28000000,89.2,paired,Control,
Sample3,True,22000000,85.1,paired,Treatment,
Sample4,True,26000000,88.3,paired,Treatment,
```

---

## Troubleshooting

### "Salmon index not found"

Ensure you've built the index and the path is correct:
```bash
ls -la /path/to/salmon_index/
# Should contain: versionInfo.json, hash.bin, etc.
```

### "Sample sheet validation failed"

Check that:
1. All FASTQ files exist at specified paths
2. File extensions are correct (.fastq, .fq, .fastq.gz, .fq.gz)
3. No duplicate sample IDs
4. For paired-end: both R1 and R2 columns are present

### Low mapping rate

Low mapping rate (<50%) may indicate:
- Wrong reference transcriptome (species mismatch)
- Contamination
- Low quality reads
- Wrong library type setting

Try:
```bash
raptor quick-count --method salmon --lib-type ISR ...  # If TruSeq
```

---

## Citation

If you use this pipeline, please cite:

**Salmon:**
```
Patro, R., et al. (2017). Salmon provides fast and bias-aware 
quantification of transcript expression. Nature Methods, 14(4), 417-419.
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
