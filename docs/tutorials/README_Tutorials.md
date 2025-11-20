# RAPTOR v2.1.0 Tutorials

Welcome to the RAPTOR tutorial series! These hands-on guides will help you master RAPTOR from basics to advanced ML-powered analysis.

---

## 🎓 Tutorial Overview

### **Beginner Tutorials** (Start Here!)

#### Tutorial 1: Getting Started with RAPTOR
**File**: [tutorial_01_getting_started.md](tutorial_01_getting_started.md)  
**Level**: Beginner | **Time**: 30 minutes  
**Prerequisites**: RAPTOR installed

**What you'll learn:**
- Basic RAPTOR commands
- Data profiling workflow
- How to interpret recommendations
- Using simulated data

**Start here if:** You're new to RAPTOR

---

#### Tutorial 5: Using the Interactive Dashboard ⭐ NEW!
**File**: [tutorial_05_dashboard.md](tutorial_05_dashboard.md)  
**Level**: Beginner | **Time**: 30 minutes  
**Prerequisites**: Tutorial 1, web browser

**What you'll learn:**
- Launching and navigating the dashboard
- Real-time resource monitoring
- Interactive data exploration
- Parameter tuning without re-running

**Start here if:** You prefer visual/interactive tools

---

### **Intermediate Tutorials**

#### Tutorial 2: Running Full Pipeline Benchmarking
**File**: [tutorial_02_benchmarking.md](tutorial_02_benchmarking.md)  
**Level**: Intermediate | **Time**: 3-6 hours (automated)  
**Prerequisites**: Tutorial 1, all tools installed

**What you'll learn:**
- Comprehensive pipeline comparison
- Interpreting benchmark results
- Making evidence-based decisions
- Publication-ready validation

**Start here if:** You need rigorous validation

---

#### Tutorial 3: Working with Your Own RNA-seq Data
**File**: [tutorial_03_real_data.md](tutorial_03_real_data.md)  
**Level**: Intermediate | **Time**: 1-2 hours  
**Prerequisites**: Tutorial 1, your own data

**What you'll learn:**
- Data preparation and formatting
- Handling common data issues
- Customizing analysis
- Integration with existing workflows

**Start here if:** You have real data ready

---

#### Tutorial 4: ML-Powered Pipeline Recommendations ⭐ NEW!
**File**: [tutorial_04_ml_recommendations.md](tutorial_04_ml_recommendations.md)  
**Level**: Intermediate | **Time**: 45 minutes  
**Prerequisites**: Tutorial 1, Python basics

**What you'll learn:**
- How RAPTOR's ML system works
- Interpreting ML confidence scores
- When to trust ML vs rule-based recommendations
- Training custom models on your data

**Start here if:** You want to understand ML recommendations

---

### **Advanced Tutorials**

#### Tutorial 6: Ensemble Analysis ⭐ NEW!
**File**: [tutorial_06_ensemble.md](tutorial_06_ensemble.md)  
**Level**: Advanced | **Time**: 2-3 hours (automated)  
**Prerequisites**: Tutorials 1-2, multiple pipelines

**What you'll learn:**
- Combining results from multiple pipelines
- Identifying high-confidence genes
- Handling discordant results
- Maximum confidence analysis

**Start here if:** You need critical validation

---

## 🎯 Learning Paths

### Path 1: Quick Start (1-2 hours)
**For users who want immediate results:**

1. **Tutorial 1** (30 min) → Learn basics with simulated data
2. **Tutorial 5** (30 min) → Use interactive dashboard
3. Jump to your analysis!

**Best for:** Exploratory analysis, getting familiar with RAPTOR

---

### Path 2: ML-Powered Analysis (2-3 hours)
**For users who want smart recommendations:**

1. **Tutorial 1** (30 min) → Understand basics
2. **Tutorial 4** (45 min) → Master ML recommendations
3. **Tutorial 5** (30 min) → Interactive exploration
4. **Tutorial 3** (1 hour) → Apply to your data

**Best for:** Most users, modern RNA-seq analysis

---

### Path 3: Comprehensive Validation (6-8 hours)
**For users publishing or making critical decisions:**

1. **Tutorial 1** (30 min) → Learn basics
2. **Tutorial 2** (4 hours) → Full benchmarking
3. **Tutorial 6** (3 hours) → Ensemble analysis
4. **Tutorial 3** (1 hour) → Apply to your data

**Best for:** Publications, drug development, critical research

---

### Path 4: Real Data Focus (2-3 hours)
**For users with data ready to analyze:**

1. **Tutorial 1** (skim, 15 min) → Understand concepts
2. **Tutorial 3** (detailed, 90 min) → Prepare and analyze your data
3. **Tutorial 4** (45 min) → Get ML recommendations
4. **Tutorial 5** (30 min) → Explore results interactively

**Best for:** Working scientists with data in hand

---

## 🆕 What's New in v2.1.0 Tutorials

### New Tutorials Added:
- ⭐ **Tutorial 4**: ML-Powered Recommendations
- ⭐ **Tutorial 5**: Interactive Dashboard
- ⭐ **Tutorial 6**: Ensemble Analysis

### Updated Content:
- Tutorial 1-3 updated with v2.1.0 features
- Added ML recommendation examples
- Included dashboard integration
- New quality assessment sections

### New Features Covered:
- 🤖 Machine learning recommendations
- 📊 Interactive web dashboard
- ⚡ Real-time resource monitoring
- 🎯 Parameter optimization
- ✅ Quality assessment
- 🔄 Ensemble analysis

---

## 📚 Tutorial Features

Each tutorial includes:
- ✅ **Clear learning objectives**
- ✅ **Step-by-step instructions**
- ✅ **Expected outputs** with examples
- ✅ **Troubleshooting** for common issues
- ✅ **Best practices** and tips
- ✅ **Real-world examples**
- ✅ **Python and command-line examples**

---

## 💻 System Requirements

### Minimum for Basic Tutorials (1, 3, 4, 5):
- 8 GB RAM (16 GB recommended)
- 20 GB free disk space
- 4 CPU cores
- Web browser (for Tutorial 5)

### Required for Benchmarking (Tutorial 2):
- 32 GB RAM (64 GB recommended)
- 100 GB free disk space
- 8+ CPU cores

### Required for Ensemble (Tutorial 6):
- 32 GB RAM (64 GB recommended)
- 150 GB free disk space
- 16+ CPU cores

---

## 🔧 Before You Start

### Software Requirements

**Required for all tutorials:**
- Python 3.8+
- RAPTOR v2.1.0+

**For full pipeline runs (Tutorials 2, 6):**
- All bioinformatics tools (see [INSTALLATION.md](../INSTALLATION.md))
- R with Bioconductor packages

### Installation

If you haven't installed RAPTOR yet:

```bash
# Quick install
pip install raptor-rnaseq>=2.1.0

# With conda (includes all tools)
conda env create -f https://raw.githubusercontent.com/AyehBlk/RAPTOR/main/environment.yml
conda activate raptor

# Verify installation
raptor --version  # Should show 2.1.0 or higher
```

See [INSTALLATION.md](../INSTALLATION.md) for detailed instructions.

---

## 📖 Additional Resources

### v2.1.0 Documentation

- [README.md](../../README.md) - Project overview
- [INSTALLATION.md](../INSTALLATION.md) - Setup guide
- [UNDERSTANDING_ML.md](../UNDERSTANDING_ML.md) - ML concepts explained
- [ML_TRAINING.md](../ML_TRAINING.md) - Advanced ML guide
- [DASHBOARD.md](../DASHBOARD.md) - Dashboard complete guide
- [ENSEMBLE.md](../ENSEMBLE.md) - Ensemble analysis guide
- [RESOURCE_MONITORING.md](../RESOURCE_MONITORING.md) - Resource tracking
- [PARAMETER_OPTIMIZATION.md](../PARAMETER_OPTIMIZATION.md) - Parameter tuning
- [QUALITY_ASSESSMENT.md](../QUALITY_ASSESSMENT.md) - Data QC guide
- [CLOUD_DEPLOYMENT.md](../CLOUD_DEPLOYMENT.md) - AWS/GCP/Azure
- [API.md](../API.md) - Python API reference v2.1.0
- [FAQ.md](../FAQ.md) - Common questions
- [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) - Problem solving

### Community Support

- **GitHub Discussions**: https://github.com/AyehBlk/RAPTOR/discussions
- **Issue Tracker**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

## 🚀 Quick Reference

### Common Commands (v2.1.0)

```bash
# Generate test data
raptor simulate --output test_data/ --size small

# Profile data with ML
raptor profile --counts counts.csv --metadata metadata.csv --use-ml

# Launch dashboard
raptor dashboard --results results/

# Quick benchmark (3 pipelines)
raptor compare --data data/ --pipelines 3 4 6 --mode quick

# Ensemble analysis
raptor ensemble --results results/ --method vote --min-agreement 2

# Monitor resources
raptor monitor --live

# Optimize parameters
raptor optimize --counts counts.csv --method bayesian

# Quality assessment
raptor qc --counts counts.csv --metadata metadata.csv

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

## 🎓 Tutorial Completion Certificate

After completing all tutorials, you'll be able to:

- ✅ Use RAPTOR for RNA-seq analysis (Tutorial 1)
- ✅ Run comprehensive benchmarking (Tutorial 2)
- ✅ Analyze your own real data (Tutorial 3)
- ✅ Leverage ML-powered recommendations (Tutorial 4)
- ✅ Use the interactive dashboard (Tutorial 5)
- ✅ Perform ensemble analysis (Tutorial 6)
- ✅ Understand quality assessment
- ✅ Monitor computational resources
- ✅ Optimize analysis parameters
- ✅ Generate publication-ready reports

**Congratulations - you're now a RAPTOR expert!** 🦖🎓

---

## 📊 Estimated Time Investment

```
Quick Start Path:           1-2 hours
ML-Powered Path:            2-3 hours
Comprehensive Path:         6-8 hours
Real Data Path:             2-3 hours

All Tutorials:              10-15 hours total
```

**Worth it?** Absolutely! Proper RNA-seq analysis can take weeks without RAPTOR.

---

## 🏆 Learning Goals by Tutorial

### Tutorial 1: Foundation
- Understand RAPTOR basics
- Profile RNA-seq data
- Interpret recommendations
- Generate reports

### Tutorial 2: Validation
- Run comprehensive benchmarks
- Compare pipeline performance
- Make evidence-based decisions
- Document for publications

### Tutorial 3: Application
- Prepare real data
- Handle data issues
- Customize analyses
- Integrate workflows

### Tutorial 4: Intelligence (NEW!)
- Understand ML recommendations
- Interpret confidence scores
- Use similar projects
- Train custom models

### Tutorial 5: Visualization (NEW!)
- Navigate dashboard
- Explore data interactively
- Monitor in real-time
- Generate reports visually

### Tutorial 6: Robustness (NEW!)
- Combine multiple pipelines
- Identify high-confidence genes
- Handle discordance
- Maximum validation

---

## 🎯 Which Tutorial Should You Start With?

### "I'm completely new to RAPTOR"
→ **Start with Tutorial 1**

### "I want the easiest way to use RAPTOR"
→ **Start with Tutorial 5 (Dashboard)**

### "I have data and need results fast"
→ **Start with Tutorial 3**

### "I want to understand the ML system"
→ **Start with Tutorial 4**

### "I need maximum confidence in my results"
→ **Complete Tutorials 1, 2, and 6**

### "I'm publishing high-impact research"
→ **Complete ALL tutorials**

---

## 💡 Tips for Success

### Do's:
✅ Follow tutorials in order if you're new  
✅ Use simulated data first to learn  
✅ Take notes as you progress  
✅ Try each command yourself  
✅ Experiment with parameters  
✅ Ask questions in discussions  

### Don'ts:
❌ Skip Tutorial 1 if you're new  
❌ Use real data before understanding basics  
❌ Copy-paste without understanding  
❌ Get discouraged by errors (they're learning opportunities!)  
❌ Rush through advanced tutorials  
❌ Ignore troubleshooting sections  

---

## 📞 Getting Help

**Stuck on a tutorial?**

1. Check the **Troubleshooting** section in that tutorial
2. Review [TROUBLESHOOTING.md](../TROUBLESHOOTING.md)
3. Search [FAQ.md](../FAQ.md)
4. Ask on [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)
5. Open an [Issue](https://github.com/AyehBlk/RAPTOR/issues) with:
   - Tutorial name and step number
   - Error message (complete)
   - Your system (OS, RAM, etc.)
   - RAPTOR version (`raptor --version`)

**For tutorial improvements:**
- Suggestions welcome!
- Pull requests appreciated
- Feedback valued

---

## 🎓 After the Tutorials

### Next Steps:
1. Apply RAPTOR to your research
2. Explore advanced features in documentation
3. Train custom ML models with your data
4. Contribute to RAPTOR development
5. Share your success stories!

### Advanced Topics:
- Custom pipeline development
- Integration with workflow managers
- Cloud deployment at scale
- Custom ML model architectures
- Batch processing multiple projects

---

## 📜 Citation

If you use these tutorials in your research or teaching:

```bibtex
@software{raptor_tutorials2025,
  author = {Ayeh Bolouki},
  title = {RAPTOR v2.1.0 Tutorials: ML-Powered RNA-seq Pipeline Selection},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/AyehBlk/RAPTOR/tree/main/docs/tutorials},
  version = {2.1.0}
}
```

---

## 🎉 Ready to Start?

### Beginners:
Jump into **[Tutorial 1: Getting Started](tutorial_01_getting_started.md)** now!

### Want Interactive Experience:
Try **[Tutorial 5: Interactive Dashboard](tutorial_05_dashboard.md)** now!

### Have Real Data:
Go to **[Tutorial 3: Working with Real Data](tutorial_03_real_data.md)** now!

---

**Happy Learning!** 🦖📚

**Tutorials created by Ayeh Bolouki**  
University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
For RAPTOR v2.1.0  
MIT License

Last updated: November 2025

---

*"From novice to expert in 6 tutorials!"* 🎓✨
