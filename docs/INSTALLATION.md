# RAPTOR Installation Guide

Complete installation guide for RAPTOR v2.0.0

## Quick Installation

### Option 1: Conda (Recommended)
```bash
# Clone repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Create environment with all tools
conda env create -f environment.yml
conda activate raptor

# Verify
raptor --version
```

### Option 2: pip
```bash
pip install raptor-rnaseq
# Then install bioinformatics tools separately
```

## System Requirements
- **Python**: 3.8+
- **RAM**: 16GB minimum, 32GB recommended  
- **Storage**: 50GB+ (100GB+ with references)
- **OS**: Linux (recommended), macOS, Windows WSL2

## Bioinformatics Tools

All tools required for 8 pipelines:
```bash
conda install -c bioconda \
  star hisat2 bowtie2 \
  rsem salmon kallisto \
  stringtie htseq subread \
  samtools fastqc multiqc
```

**R/Bioconductor packages:**
```r
BiocManager::install(c("DESeq2", "edgeR", "limma",
                       "tximport", "ballgown", "NOISeq",  
                       "EBSeq", "sleuth"))
```

## Verification
```bash
raptor --version
raptor check-tools
raptor demo  # Run quick test
```

## Reference Data

Download reference genome and build indices for:
- STAR, HISAT2, Bowtie2 (alignment tools)
- Salmon, Kallisto (quantification tools)

See full guide: https://github.com/AyehBlk/RAPTOR/docs/INSTALLATION.md

## Troubleshooting

**Tool not found:** Add to PATH or reinstall via conda
**Import errors:** `pip install --force-reinstall raptor-rnaseq`
**Memory errors:** Reduce threads/memory in config

Contact: ayehbolouki1988@gmail.com
