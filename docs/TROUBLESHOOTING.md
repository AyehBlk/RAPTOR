# RAPTOR Troubleshooting Guide

Solutions to common problems and errors.

## Installation Issues

### Problem: "Command not found: raptor"

**Symptoms**: `raptor: command not found` after installation

**Solutions**:
```bash
# Solution 1: Check if installed
pip list | grep raptor

# Solution 2: Add to PATH
export PATH=$PATH:~/.local/bin
echo 'export PATH=$PATH:~/.local/bin' >> ~/.bashrc

# Solution 3: Reinstall
pip install --force-reinstall raptor-rnaseq

# Solution 4: Use python -m
python -m raptor.cli --help
```

### Problem: Import errors for Python packages

**Symptoms**: `ModuleNotFoundError: No module named 'pandas'`

**Solutions**:
```bash
# Reinstall with dependencies
pip install --force-reinstall raptor-rnaseq[all]

# Or install missing packages
pip install numpy pandas scipy matplotlib seaborn scikit-learn
```

### Problem: Conda environment conflicts

**Symptoms**: Dependency conflicts during conda install

**Solutions**:
```bash
# Remove and recreate
conda deactivate
conda env remove -n raptor
conda env create -f environment.yml

# Or use mamba (faster solver)
mamba env create -f environment.yml
```

## Tool Installation Issues

### Problem: "Tool not found: STAR"

**Symptoms**: `STAR: command not found`

**Solutions**:
```bash
# Check if installed
which STAR

# Install via conda
conda install -c bioconda star

# Add to PATH if installed
export PATH=$PATH:/path/to/STAR/bin

# Verify
STAR --version
```

### Problem: R package installation fails

**Symptoms**: Error installing Bioconductor packages

**Solutions**:
```r
# Update BiocManager
install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Install specific package
BiocManager::install("DESeq2", force = TRUE)

# Check dependencies
BiocManager::valid()
```

## Runtime Issues

### Problem: "Killed" error during analysis

**Symptoms**: Process terminates with "Killed"

**Cause**: Insufficient memory (OOM)

**Solutions**:
```bash
# Reduce resources
raptor profile --counts data.csv --threads 4 --memory 16G

# Use memory-efficient pipeline
raptor profile --counts data.csv --fast

# Check available memory
free -h

# Monitor during run
watch -n 1 free -h
```

### Problem: Very slow performance

**Symptoms**: Analysis takes much longer than expected

**Solutions**:
```bash
# Increase threads
raptor compare --data fastq/ --threads 16

# Use SSD, not HDD
# Move data to SSD: cp -r data/ /ssd/data/

# Check CPU usage
top

# Use faster pipeline
raptor run --pipeline 3  # Salmon-edgeR
```

### Problem: Disk space errors

**Symptoms**: "No space left on device"

**Solutions**:
```bash
# Check space
df -h

# Clean up intermediate files
raptor clean --output-dir results/

# Use compression
raptor compare --compress-results

# Change output location
raptor compare --output /path/with/space/
```

## Data Issues

### Problem: "Invalid count matrix"

**Symptoms**: ValueError during profiling

**Solutions**:
```python
# Check your data
import pandas as pd
counts = pd.read_csv('counts.csv', index_col=0)

# Verify format
print(counts.shape)  # Should be genes × samples
print(counts.dtypes)  # Should be numeric
print((counts < 0).any())  # Should be False

# Fix negative values
counts[counts < 0] = 0

# Fix non-numeric
counts = counts.apply(pd.to_numeric, errors='coerce')
```

### Problem: Sample/metadata mismatch

**Symptoms**: Warning about missing samples

**Solutions**:
```python
# Check names
count_samples = set(counts.columns)
meta_samples = set(metadata['sample'])

# Find mismatches
print("In counts only:", count_samples - meta_samples)
print("In metadata only:", meta_samples - count_samples)

# Fix names (if just case/whitespace)
counts.columns = counts.columns.str.strip().str.lower()
metadata['sample'] = metadata['sample'].str.strip().str.lower()
```

### Problem: Too many zeros

**Symptoms**: Warning about high zero-inflation

**Solutions**:
```python
# Filter low-count genes
min_count = 10
min_samples = 2
keep = (counts >= min_count).sum(axis=1) >= min_samples
counts_filtered = counts.loc[keep]

# Or use CPM filtering
cpm = counts / counts.sum(axis=0) * 1e6
keep = (cpm >= 1).sum(axis=1) >= min_samples
counts_filtered = counts.loc[keep]
```

## Pipeline-Specific Issues

### Problem: STAR alignment fails

**Symptoms**: STAR terminates with error

**Common causes and solutions**:

**Issue**: Genome index not found
```bash
# Build STAR index
STAR --runMode genomeGenerate \
  --genomeDir star_index/ \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile genes.gtf \
  --runThreadN 8
```

**Issue**: Out of memory
```bash
# Reduce memory with --limitGenomeGenerateRAM
STAR --runMode genomeGenerate \
  --limitGenomeGenerateRAM 31000000000 \
  ...
```

### Problem: Salmon quantification fails

**Symptoms**: Salmon exits with error

**Solutions**:
```bash
# Check index
salmon index -t transcripts.fa -i salmon_index

# Check FASTQ format
head sample_R1.fastq

# Try with --validateMappings
salmon quant -i salmon_index \
  --validateMappings \
  -l A -1 R1.fq -2 R2.fq \
  -o output
```

### Problem: DESeq2 fails in R

**Symptoms**: Error in DESeq2 analysis

**Solutions**:
```r
# Check count matrix
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)

# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run with error catching
tryCatch({
    dds <- DESeq(dds)
}, error = function(e) {
    print(paste("Error:", e$message))
})
```

## Configuration Issues

### Problem: Config file not loading

**Symptoms**: `FileNotFoundError: config.yaml not found`

**Solutions**:
```bash
# Check path
ls config/config.yaml

# Use absolute path
raptor profile --config /full/path/to/config.yaml

# Or copy default
cp /path/to/raptor/config/config.yaml ./my_config.yaml
```

### Problem: Invalid configuration values

**Symptoms**: Warning about config validation

**Solutions**:
```python
# Validate config
from raptor.utils import load_config
config = load_config('my_config.yaml')

# Check critical values
assert 0 < config['statistics']['fdr_threshold'] <= 1
assert config['resources']['default_threads'] > 0
```

## Output & Results Issues

### Problem: Empty output directory

**Symptoms**: No results generated

**Solutions**:
```bash
# Check for errors in log
cat raptor.log

# Run with verbose output
raptor compare --data fastq/ --verbose --log-level DEBUG

# Check permissions
ls -la output_dir/
chmod 755 output_dir/
```

### Problem: Cannot open HTML report

**Symptoms**: Report doesn't display correctly

**Solutions**:
```bash
# Check file exists
ls -lh report.html

# Try different browser
# Chrome, Firefox, Safari

# Check for JavaScript errors
# Open browser console (F12)

# Regenerate report
raptor profile --counts data.csv --output new_report.html
```

## Permission Issues

### Problem: Permission denied errors

**Symptoms**: `PermissionError: [Errno 13]`

**Solutions**:
```bash
# Check file permissions
ls -l file.csv

# Make readable
chmod 644 file.csv

# Make directory writable
chmod 755 output_dir/

# Check ownership
chown $USER:$USER file.csv
```

## Network Issues

### Problem: Cannot download references

**Symptoms**: Download fails or times out

**Solutions**:
```bash
# Use wget with retry
wget -c --tries=10 URL

# Use aria2 (faster, parallel)
aria2c -x 16 URL

# Use mirror sites
# GENCODE: multiple mirrors available

# Check proxy settings
echo $http_proxy
```

## Getting More Help

If your issue isn't listed:

1. **Check logs**: `cat raptor.log`
2. **Search issues**: GitHub Issues page
3. **Ask community**: GitHub Discussions
4. **Report bug**: Open new issue with:
   - RAPTOR version
   - Operating system
   - Complete error message
   - Minimal reproducible example

5. **Email**: ayehbolouki1988@gmail.com

## Debug Mode

Enable debug mode for detailed output:

```bash
# CLI
raptor profile --counts data.csv --debug --log-level DEBUG

# Python
import logging
logging.basicConfig(level=logging.DEBUG)

from raptor import RNAseqDataProfiler
profiler = RNAseqDataProfiler(counts, metadata, verbose=True)
```

## System Information

Collect system info for bug reports:

```bash
# RAPTOR version
raptor --version

# Python version
python --version

# Tool versions
raptor check-tools

# System info
uname -a
free -h
df -h

# Environment
conda list  # if using conda
pip list  # if using pip
```

---

**Still stuck?** Open an issue on GitHub with details!

