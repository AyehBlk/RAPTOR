# RAPTOR Tutorials

Welcome to the RAPTOR tutorial series! These hands-on guides will help you learn RAPTOR from basics to advanced usage.

---

##  Tutorial Overview

### Tutorial 1: Getting Started with RAPTOR
**File**: [tutorial_01_getting_started.md](tutorial_01_getting_started.md)  
**Level**: Beginner  
**Time**: 30 minutes  
**Prerequisites**: RAPTOR installed

**What you'll learn:**
- Basic RAPTOR commands
- Data profiling workflow
- How to interpret recommendations
- Using simulated data

**Start here if:** You're new to RAPTOR or want a quick introduction.

---

### Tutorial 2: Running Full Pipeline Benchmarking
**File**: [tutorial_02_benchmarking.md](tutorial_02_benchmarking.md)  
**Level**: Intermediate  
**Time**: 3-6 hours (mostly automated)  
**Prerequisites**: Tutorial 1 completed, all tools installed

**What you'll learn:**
- Comprehensive pipeline comparison
- Interpreting benchmark results
- Making evidence-based decisions
- Documentation for publications

**Start here if:** You need rigorous validation or are publishing results.

---

### Tutorial 3: Working with Your Own RNA-seq Data
**File**: [tutorial_03_real_data.md](tutorial_03_real_data.md)  
**Level**: Intermediate  
**Time**: 1-2 hours  
**Prerequisites**: Tutorial 1 completed, your own data

**What you'll learn:**
- Data preparation and formatting
- Handling common data issues
- Customizing analysis
- Integration with existing workflows

**Start here if:** You have your own RNA-seq data ready to analyze.

---

##  Learning Paths

### Path 1: Quick Start (1 hour)
For users who want to get started immediately:
1. **Tutorial 1** → Learn basics with simulated data
2. Jump to your analysis!

### Path 2: Comprehensive Learning (4-6 hours)
For users who want deep understanding:
1. **Tutorial 1** → Learn basics
2. **Tutorial 2** → Understand benchmarking
3. **Tutorial 3** → Apply to real data

### Path 3: Real Data Focus (2-3 hours)
For users with data ready:
1. **Tutorial 1** (skim) → Understand concepts
2. **Tutorial 3** (detailed) → Apply to your data
3. **Tutorial 2** (if needed) → Validate choices

---

##  Tutorial Features

Each tutorial includes:
- ✅ **Clear learning objectives**
- ✅ **Step-by-step instructions**
- ✅ **Expected outputs** with examples
- ✅ **Troubleshooting** for common issues
- ✅ **Best practices** and tips
- ✅ **Real-world examples**

---

##  Before You Start

### System Requirements

Minimum for tutorials:
- 8 GB RAM (16 GB recommended)
- 20 GB free disk space
- 4 CPU cores

For Tutorial 2 (benchmarking):
- 32 GB RAM (64 GB recommended)
- 100 GB free disk space
- 8+ CPU cores

### Software Requirements

**Required:**
- Python 3.8+
- RAPTOR v2.0.0+

**For full pipeline runs:**
- All bioinformatics tools (see [INSTALLATION.md](../INSTALLATION.md))
- R with Bioconductor packages

### Installation

If you haven't installed RAPTOR yet:

```bash
# Quick install
pip install raptor-rnaseq

# With conda (includes all tools)
conda env create -f https://raw.githubusercontent.com/AyehBlk/RAPTOR/main/environment.yml
conda activate raptor
```

See [INSTALLATION.md](../INSTALLATION.md) for detailed instructions.

---

##  Additional Resources

### Documentation

- [README.md](../../README.md) - Project overview
- [INSTALLATION.md](../INSTALLATION.md) - Setup guide
- [PROFILE_RECOMMEND.md](../PROFILE_RECOMMEND.md) - Profiling details
- [BENCHMARKING.md](../BENCHMARKING.md) - Benchmarking guide
- [PIPELINES.md](../PIPELINES.md) - Pipeline descriptions
- [API.md](../API.md) - Python API reference
- [FAQ.md](../FAQ.md) - Common questions
- [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) - Problem solving

### Community Support

- **GitHub Discussions**: https://github.com/AyehBlk/RAPTOR/discussions
- **Issue Tracker**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

##  Quick Reference

### Common Commands

```bash
# Generate test data
raptor simulate --output test_data/ --size small

# Profile data
raptor profile --counts counts.csv --metadata metadata.csv

# Quick benchmark (3 pipelines)
raptor compare --data data/ --pipelines 3 4 6 --mode quick

# Full benchmark (all 8 pipelines)
raptor compare --data data/ --pipelines all --mode full

# Generate report
raptor report --results results/ --output report.html

# Get help
raptor --help
raptor profile --help
```

### File Formats

**Count Matrix (counts.csv):**
```csv
gene_id,Sample1,Sample2,Sample3
Gene1,100,150,120
Gene2,50,45,60
```

**Metadata (metadata.csv):**
```csv
sample,condition,replicate
Sample1,Control,1
Sample2,Control,2
Sample3,Treatment,1
```

---

##  Citation

If you use these tutorials in your research or teaching:

```bibtex
@software{raptor_tutorials2025,
  author = {Ayeh Bolouki},
  title = {RAPTOR Tutorials: Hands-on Guides for RNA-seq Pipeline Selection},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/AyehBlk/RAPTOR/tree/main/docs/tutorials}
}
```

---

##  Ready to Start?

Jump into **[Tutorial 1: Getting Started](tutorial_01_getting_started.md)** now!

---

**Tutorials created by Ayeh Bolouki**  
University of Namur, Belgium  
For RAPTOR v2.0.0  
MIT License

Last updated: January 2025
