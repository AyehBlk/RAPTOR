# RAPTOR Pipelines Reference

Deep dive into all 8 RNA-seq analysis pipelines.

## Pipeline Overview

| ID | Name | Type | Speed | Accuracy | Memory |
|----|------|------|-------|----------|--------|
| 1 | STAR-RSEM-DESeq2 | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚫⚫ | High |
| 2 | HISAT2-StringTie-Ballgown | Alignment | ⚫⚫⚫⚪⚪ | ⚫⚫⚫⚫⚪ | Medium |
| 3 | Salmon-edgeR | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚫⚪ | Low |
| 4 | Kallisto-Sleuth | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚪⚪ | Low |
| 5 | STAR-HTSeq-limma | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚫⚪ | High |
| 6 | STAR-featureCounts-NOISeq | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚪⚪ | High |
| 7 | Bowtie2-RSEM-EBSeq | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚪⚪ | Medium |
| 8 | HISAT2-Cufflinks-Cuffdiff | Alignment | ⚫⚪⚪⚪⚪ | ⚫⚫⚪⚪⚪ | Medium |

## Pipeline 1: STAR-RSEM-DESeq2

**Gold Standard - Highest Accuracy**

### Components
- **Alignment**: STAR (splice-aware)
- **Quantification**: RSEM (EM algorithm)
- **Statistics**: DESeq2 (negative binomial)

### Best For
- Publication-quality results
- Small sample sizes (n<6)
- Complex experimental designs
- When accuracy is paramount

### Performance
- Runtime: ~2-6 hours (typical dataset)
- Memory: 32-48 GB
- Accuracy: 95%

### Running
```bash
raptor run --pipeline 1 \
  --data fastq/ \
  --reference /path/to/star_index \
  --annotation genes.gtf \
  --output pipeline1_results/
```

## Pipeline 2: HISAT2-StringTie-Ballgown

**Transcript Assembly & Novel Discovery**

### Components
- **Alignment**: HISAT2
- **Assembly**: StringTie
- **Statistics**: Ballgown

### Best For
- Novel transcript discovery
- Isoform-level analysis
- Non-model organisms
- When reference incomplete

### Performance
- Runtime: ~1-4 hours
- Memory: 16-24 GB
- Accuracy: 88%

## Pipeline 3: Salmon-edgeR ⭐ RECOMMENDED

**Best Balance - Fast & Accurate**

### Components
- **Quantification**: Salmon (quasi-mapping)
- **Statistics**: edgeR (quasi-likelihood)

### Best For
- Most RNA-seq experiments
- Large datasets (>20 samples)
- Quick turnaround needed
- Good balance of all metrics

### Performance
- Runtime: ~0.5-2 hours
- Memory: 8-16 GB
- Accuracy: 90%

### Running
```bash
raptor run --pipeline 3 \
  --data fastq/ \
  --transcriptome /path/to/salmon_index \
  --output pipeline3_results/
```

## Pipeline 4: Kallisto-Sleuth

**Ultra-Fast - Large Studies**

### Components
- **Quantification**: Kallisto
- **Statistics**: Sleuth (bootstrap-based)

### Best For
- Very large datasets (>50 samples)
- Exploratory analysis
- Minimal resources
- Speed is critical

### Performance
- Runtime: ~0.3-1 hour
- Memory: 4-8 GB
- Accuracy: 88%

## Pipeline 5: STAR-HTSeq-limma-voom

**Flexible Modeling - Complex Designs**

### Components
- **Alignment**: STAR
- **Counting**: HTSeq
- **Statistics**: limma-voom

### Best For
- Complex experimental designs
- Multi-factor analysis
- Batch correction needed
- Repeated measures

### Performance
- Runtime: ~2-7 hours
- Memory: 32-40 GB
- Accuracy: 92%

## Pipeline 6: STAR-featureCounts-NOISeq

**Non-Parametric - Small Samples**

### Components
- **Alignment**: STAR
- **Counting**: featureCounts
- **Statistics**: NOISeq

### Best For
- Very small samples (n=2-3)
- No replicates
- Non-normal distributions

### Performance
- Runtime: ~2-7 hours
- Memory: 32-36 GB
- Accuracy: 85%

## Pipeline 7: Bowtie2-RSEM-EBSeq

**Memory-Efficient Alternative**

### Components
- **Alignment**: Bowtie2
- **Quantification**: RSEM
- **Statistics**: EBSeq

### Best For
- Moderate resource environments
- Isoform-level analysis
- Two-condition comparisons

### Performance
- Runtime: ~3-8 hours
- Memory: 16-24 GB
- Accuracy: 87%

## Pipeline 8: HISAT2-Cufflinks-Cuffdiff

**Legacy Pipeline**

### Components
- **Alignment**: HISAT2
- **Assembly**: Cufflinks
- **Statistics**: Cuffdiff

### Best For
- Reproducing legacy analyses
- When Cufflinks ecosystem required

### Performance
- Runtime: ~4-12 hours
- Memory: 20-32 GB
- Accuracy: 82%

**Note**: Newer methods preferred for new projects

## Choosing a Pipeline

Use the decision guide:

1. **Need highest accuracy?** → Pipeline 1
2. **Large dataset?** → Pipeline 3 or 4
3. **Novel transcripts?** → Pipeline 2
4. **Complex design?** → Pipeline 5
5. **Small samples?** → Pipeline 6
6. **Not sure?** → Use `raptor profile`

See configuration details in `config/pipelines.yaml`

