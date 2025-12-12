# 🦖 RAPTOR Installation Guide

Complete installation instructions for RAPTOR (RNA-seq Analysis Pipeline Testing and Optimization Resource).

---

## System Requirements

### Minimum Requirements
- **Python**: 3.8 or higher
- **RAM**: 8 GB
- **Disk**: 10 GB free space
- **OS**: Linux, macOS, or Windows (with WSL)

### Recommended
- **Python**: 3.10+
- **RAM**: 16 GB
- **Disk**: 50 GB free space
- **CPU**: 4+ cores

### For Full Pipeline Execution
- **R**: 4.0 or higher
- **Bioconductor packages**: DESeq2, edgeR, limma

---

## Installation Methods

### Method 1: Install from PyPI (Recommended) ⭐

The simplest way to install RAPTOR:

```bash
pip install raptor-rnaseq
```

#### With Optional Features

```bash
# With interactive dashboard
pip install raptor-rnaseq[dashboard]

# With ML features
pip install raptor-rnaseq[ml]

# With advanced features (Bayesian optimization, etc.)
pip install raptor-rnaseq[advanced]

# With development tools
pip install raptor-rnaseq[dev]

# With ALL features (recommended for full experience)
pip install raptor-rnaseq[all]
```

#### Verify Installation

```bash
python -c "import raptor; print(f'RAPTOR {raptor.__version__} installed successfully!')"
```

---

### Method 2: Install from GitHub

For the latest development version:

```bash
# Clone repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Install in development mode
pip install -e .

# Or install dependencies manually
pip install -r requirements.txt
```

---

### Method 3: Conda Environment

For isolated environment with all dependencies:

```bash
# Create environment from file
conda env create -f environment.yml
conda activate raptor

# Or create manually
conda create -n raptor python=3.10
conda activate raptor
pip install raptor-rnaseq[all]
```

---

## Platform-Specific Instructions

### Linux (Ubuntu/Debian)

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y python3-pip python3-dev

# Install RAPTOR
pip install raptor-rnaseq[all]
```

### macOS

```bash
# Using Homebrew
brew install python

# Install RAPTOR
pip3 install raptor-rnaseq[all]
```

### Windows

**Option A: Windows Subsystem for Linux (Recommended)**
```bash
# In WSL terminal
pip install raptor-rnaseq[all]
```

**Option B: Native Windows**
```cmd
# In Command Prompt or PowerShell
py -m pip install raptor-rnaseq[all]
```

---

## Installing R Dependencies (Optional)

Required only if you want to run the actual RNA-seq pipelines:

### Automatic Installation

```bash
Rscript scripts/install_r_packages.R
```

### Manual Installation

```r
# In R console
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2",
    "edgeR", 
    "limma",
    "tximport",
    "GenomicFeatures"
))
```

---

## Installing RNA-seq Tools (Optional)

Required only for running full pipelines:

### Salmon
```bash
conda install -c bioconda salmon
# or download from https://github.com/COMBINE-lab/salmon/releases
```

### Kallisto
```bash
conda install -c bioconda kallisto
# or download from https://pachterlab.github.io/kallisto/download
```

### STAR
```bash
conda install -c bioconda star
# or download from https://github.com/alexdobin/STAR/releases
```

---

## Verify Full Installation

Run the verification script:

```bash
python install.py
```

Or manually verify each component:

```python
# Python
import raptor
print(f"RAPTOR version: {raptor.__version__}")

# Core modules
from raptor import RNAseqDataProfiler, MLPipelineRecommender
from raptor.data_quality_assessment import DataQualityAssessor
from raptor.ensemble_analysis import EnsembleAnalyzer

print("✅ All modules imported successfully!")
```

---

## Quick Test

After installation, test with:

```bash
# Launch dashboard
python -c "import raptor" && echo "✅ RAPTOR ready!"

# Or launch the dashboard
raptor --help
```

---

## Troubleshooting

### Common Issues

#### "ModuleNotFoundError: No module named 'raptor'"
```bash
pip install raptor-rnaseq
```

#### "pip: command not found"
```bash
python -m pip install raptor-rnaseq
# or on Windows
py -m pip install raptor-rnaseq
```

#### Permission errors
```bash
pip install --user raptor-rnaseq
```

#### Outdated pip
```bash
pip install --upgrade pip
pip install raptor-rnaseq
```

#### Conflicts with existing packages
```bash
# Create fresh virtual environment
python -m venv raptor_env
source raptor_env/bin/activate  # Linux/macOS
# or
raptor_env\Scripts\activate  # Windows

pip install raptor-rnaseq[all]
```

---

## Updating RAPTOR

### From PyPI
```bash
pip install --upgrade raptor-rnaseq
```

### From GitHub
```bash
cd RAPTOR
git pull
pip install -e .
```

---

## Uninstalling

```bash
pip uninstall raptor-rnaseq
```

---

## Next Steps

After installation:

1. **Launch the dashboard**: `python launch_dashboard.py`
2. **Read the Quick Start**: [QUICK_START.md](QUICK_START.md)
3. **Try an example**: [examples/](examples/)
4. **Read the docs**: [docs/](docs/)

---

## Support

- **PyPI**: https://pypi.org/project/raptor-rnaseq/
- **GitHub**: https://github.com/AyehBlk/RAPTOR
- **Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

**Author:** Ayeh Bolouki  
**License:** MIT
