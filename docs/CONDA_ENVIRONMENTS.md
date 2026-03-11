# Conda Environment Files - Complete Guide

**RAPTOR v2.2.0**  
**Date:** January 2026  
**Author:** Ayeh Bolouki

---

## 📋 **Overview: Two Environment Files**

| File | Size | Time | Use Case |
|------|------|------|----------|
| **environment.yml** | ~500 MB | 5-10 min | Python package only |
| **environment-full.yml** | ~5-8 GB | 30-60 min | Complete bioinformatics suite |

---

## 🎯 **Which File Should You Use?**

### **Use `environment.yml` if:**
- ✅ You only need RAPTOR Python package
- ✅ You'll install bioinformatics tools separately
- ✅ You want fast installation
- ✅ You already have STAR, Salmon, etc. installed
- ✅ You're developing RAPTOR code

**Installation:**
```bash
conda env create -f environment.yml
conda activate raptor
```

---

### **Use `environment-full.yml` if:**
- ✅ You want everything in one environment
- ✅ You need all bioinformatics tools
- ✅ You're setting up a new analysis workstation
- ✅ You want reproducible complete environment
- ✅ You don't mind long installation time

**Installation:**
```bash
conda env create -f environment-full.yml
conda activate raptor-full
```

---

## 📊 **Changes from v2.1.1 → v2.2.0**

### **REMOVED from Required:**

1. **bayesian-optimization** - Moved to optional (pip section)
2. **jinja2** - Moved to optional (only for PDF reports)
3. **markdown** - Moved to optional (only for PDF reports)
4. **psutil** - Moved to optional (only for monitoring)
5. **colorlog** - Moved to optional (only for colored logs)

**Impact:** ~200 MB smaller, 2-3 minutes faster

---

### **ADDED to Required:**

1. **toml>=0.10.0** - NEW in v2.2.0
   - Needed for TOML configuration files
   - Used by dashboard and app settings

---

### **VERSION UPDATES:**

Aligned with updated `requirements.txt`:

| Package | v2.1.1 | v2.2.0 | Reason |
|---------|--------|--------|--------|
| Python | =3.10 | >=3.8,<3.13 | Flexibility |
| numpy | >=1.21.0 | >=1.19.0 | Broader compatibility |
| pandas | >=1.3.0 | >=1.1.0 | Python 3.8 support |
| scipy | >=1.7.0 | >=1.5.0 | Compatibility |
| scikit-learn | >=1.0.0 | >=0.24.0 | Python 3.8 |
| statsmodels | >=0.13.0 | >=0.12.0 | Stability |
| matplotlib | >=3.4.0 | >=3.3.0 | Compatibility |
| streamlit | >=1.28.0 | >=1.20.0 | Broader support |
| plotly | >=5.14.0 | >=5.0.0 | Compatibility |
| click | >=8.0.0 | >=7.0 | Compatibility |
| tqdm | >=4.62.0 | >=4.50.0 | Compatibility |
| pytest | >=7.0.0 | >=6.0.0 | Compatibility |
| black | >=22.0 | >=21.0 | Compatibility |
| sphinx | >=4.0.0 | >=3.5.0 | Compatibility |

**Why lower?** Ensures compatibility with Python 3.8 and reduces conflicts

---

### **REORGANIZATION:**

Better organized by function:

```yaml
# OLD (v2.1.1) - Mixed organization
- numpy
- streamlit  
- star
- pytest

# NEW (v2.2.0) - Grouped by purpose
# Core Scientific Computing
- numpy, pandas, scipy

# Visualization  
- matplotlib, seaborn, plotly

# Optional: Dashboard
- streamlit, altair

# Optional: Development
- pytest, black, flake8
```

---

## 📁 **File Structure Comparison**

### **environment.yml (Core)**

```yaml
name: raptor
channels:
  - conda-forge
  - defaults

dependencies:
  # Python 3.8-3.12
  - python>=3.8,<3.13
  
  # Core packages (14 total)
  - numpy, pandas, scipy
  - matplotlib, seaborn, plotly
  - statsmodels
  - pyyaml, toml
  - click, colorama, tqdm
  
  # Optional (commented out by default)
  # - streamlit, altair (dashboard)
  # - pytest, black, flake8 (dev)
  
  - pip:
    - raptor-rnaseq>=2.2.0
```

**Size:** ~500 MB  
**Time:** 5-10 minutes  
**Packages:** ~50

---

### **environment-full.yml (Complete)**

```yaml
name: raptor-full
channels:
  - conda-forge
  - bioconda  # ← Includes bioinformatics tools
  - defaults

dependencies:
  # Python 3.8-3.12
  - python>=3.8,<3.13
  
  # All core packages
  - numpy, pandas, scipy...
  
  # All optional features
  - streamlit, altair
  - pytest, black, flake8
  - sphinx, sphinx_rtd_theme
  
  # Bioinformatics tools (NOT in environment.yml)
  - star, hisat2, bowtie2
  - salmon, kallisto, rsem
  - samtools, bedtools, sra-tools
  
  # R + Bioconductor
  - r-base
  - bioconductor-deseq2
  - bioconductor-edger
  - bioconductor-limma
  
  - pip:
    - raptor-rnaseq[all]>=2.2.0
    - bayesian-optimization
    - jinja2, markdown
```

**Size:** ~5-8 GB  
**Time:** 30-60 minutes  
**Packages:** ~200+

---

## 🚀 **Quick Start Examples**

### **Example 1: Python Developer**

```bash
# Use core environment
conda env create -f environment.yml
conda activate raptor

# Install RAPTOR in editable mode
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -e .[dev]

# Start coding!
```

---

### **Example 2: Bioinformatics Researcher**

```bash
# Use full environment (grab coffee, takes 30-60 min)
conda env create -f environment-full.yml
conda activate raptor-full

# Verify everything
raptor --version
STAR --version
salmon --version

# Ready to analyze data!
```

---

### **Example 3: Dashboard User**

```bash
# Use core environment
conda env create -f environment.yml
conda activate raptor

# Uncomment dashboard lines in environment.yml first, then:
conda env update -f environment.yml

# Launch dashboard
streamlit run raptor/dashboard/app.py
```

---

## 🔧 **Customizing Environment Files**

### **Add Dashboard Support:**

In `environment.yml`, uncomment:
```yaml
# ==========================================================================
# Optional: Dashboard (comment out if not needed)
# ==========================================================================
- streamlit>=1.20.0
- altair>=4.0.0
```

Then:
```bash
conda env update -f environment.yml
```

---

### **Add Development Tools:**

In `environment.yml`, uncomment:
```yaml
# ==========================================================================
# Optional: Development Tools (comment out for production)
# ==========================================================================
- pytest>=6.0.0
- black>=21.0
- flake8>=3.8.0
```

---

### **Remove Bioinformatics Tools:**

In `environment-full.yml`, comment out:
```yaml
# ==========================================================================
# MODULE 5: Alignment Tools (Production Pipelines)
# ==========================================================================
# - star>=2.7.10a
# - hisat2>=2.2.1
```

---

## 📊 **Installation Comparison**

| Aspect | environment.yml | environment-full.yml |
|--------|----------------|----------------------|
| **Python** | ✅ | ✅ |
| **RAPTOR package** | ✅ | ✅ |
| **Dashboard** | Optional | ✅ |
| **Dev tools** | Optional | ✅ |
| **STAR aligner** | ❌ | ✅ |
| **Salmon** | ❌ | ✅ |
| **Kallisto** | ❌ | ✅ |
| **R + DESeq2** | ❌ | ✅ |
| **Install time** | 5-10 min | 30-60 min |
| **Disk space** | ~500 MB | ~5-8 GB |

---

## 🐛 **Troubleshooting**

### **Problem: Installation Takes Forever**

**Solution:** Use `environment.yml` instead of `environment-full.yml`

```bash
# Fast installation
conda env create -f environment.yml

# Install bioinformatics tools separately as needed
conda install -c bioconda star salmon kallisto
```

---

### **Problem: Conda Solver Takes Long**

**Solution:** Use mamba (faster conda solver)

```bash
# Install mamba
conda install -c conda-forge mamba

# Use mamba instead of conda
mamba env create -f environment.yml
```

---

### **Problem: Conflicts During Installation**

**Solution:** Create with specific Python version

```bash
conda env create -f environment.yml python=3.10
```

Or install in stages:
```bash
# Stage 1: Core
conda create -n raptor python=3.10
conda activate raptor
pip install raptor-rnaseq

# Stage 2: Tools (if needed)
conda install -c bioconda star salmon
```

---

### **Problem: Windows Installation Fails**

**Solution:** Some tools not available on Windows

**Options:**
1. Use `environment.yml` (Python only)
2. Use WSL2 (Windows Subsystem for Linux)
3. Use Docker

```bash
# WSL2
wsl
conda env create -f environment-full.yml

# Docker
docker pull ayehblk/raptor:2.2.0
```

---

## ✅ **Verification After Installation**

### **Test Core Installation:**

```bash
conda activate raptor

# Test Python
python --version

# Test RAPTOR
raptor --version
python -c "import raptor; print(raptor.__version__)"

# Test modules
python -c "from raptor import ensemble_fisher; print('✅')"
```

---

### **Test Full Installation:**

```bash
conda activate raptor-full

# Test RAPTOR
raptor --version

# Test bioinformatics tools
STAR --version
salmon --version
kallisto version
samtools --version

# Test R packages
Rscript -e "library(DESeq2); print('DESeq2 OK')"
Rscript -e "library(edgeR); print('edgeR OK')"
```

---

## 🎯 **Best Practices**

### **For Development:**
```bash
# Use core environment + editable install
conda env create -f environment.yml
conda activate raptor
pip install -e .[dev]
```

### **For Production:**
```bash
# Use core environment + stable install
conda env create -f environment.yml
conda activate raptor
pip install raptor-rnaseq==2.2.0
```

### **For Teaching/Workshops:**
```bash
# Use full environment (everything included)
conda env create -f environment-full.yml
conda activate raptor-full
```

### **For HPC/Cluster:**
```bash
# Use module system + core environment
module load star salmon kallisto
conda env create -f environment.yml
```

---

## 📚 **Summary**

### **Key Changes v2.1.1 → v2.2.0:**
- ✅ Split into two files (core vs full)
- ✅ Removed 5 packages from required
- ✅ Added `toml` package
- ✅ Updated version numbers
- ✅ Better organization
- ✅ Faster installation

### **Which to Use:**
- **Most users:** `environment.yml` (fast, minimal)
- **Complete setup:** `environment-full.yml` (slow, comprehensive)
- **Windows users:** `environment.yml` only (or use WSL2/Docker)

### **Migration from v2.1.1:**
```bash
# Remove old environment
conda env remove -n raptor

# Create new environment
conda env create -f environment.yml

# Test
conda activate raptor
raptor --version
```

---

## 📝 **Action Items**

1. **Choose your file:**
   - Need bioinformatics tools? → `environment-full.yml`
   - Python only? → `environment.yml`

2. **Create environment:**
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate and test:**
   ```bash
   conda activate raptor
   raptor --version
   ```

4. **Done!** ✅

---

**Questions?** Check the troubleshooting section or open an issue on GitHub!
