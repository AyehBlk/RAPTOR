# RAPTOR v2.1.0 Tutorials

**Welcome to RAPTOR Tutorials!** 🦖

Step-by-step guides to master RNA-seq pipeline selection with ML-powered recommendations, ensemble analysis, and interactive dashboards.

---

##  Tutorial Overview

| Tutorial | Level | Time | Focus |
|----------|-------|------|-------|
| [Tutorial 01: Getting Started](tutorial_01_getting_started.md) | Beginner | 30-45 min | ML recommendations, dashboard basics |
| [Tutorial 02: Advanced Benchmarking](tutorial_02_benchmarking.md) | Intermediate | 3-6 hours | Compare pipelines, ensemble analysis |
| [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md) | Intermediate | 2-4 hours | Quality assessment, production workflows |
| [Tutorial 04: ML Recommendations](tutorial_04_ml_recommendations.md) | Advanced | 1-2 hours | Understanding ML, custom models |
| [Tutorial 05: Interactive Dashboard](tutorial_05_dashboard.md) | Intermediate | 1 hour | Web interface, visual analysis |
| [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md) | Advanced | 2-3 hours | Combining pipelines, robust results |

---

##  Learning Paths

###  **New to RAPTOR?**

**Start here**:
1. [Tutorial 01: Getting Started](tutorial_01_getting_started.md)
   - Install and configure RAPTOR
   - Get your first ML recommendations
   - Explore the dashboard
   - Run a simple analysis

2. [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md)
   - Prepare your own data
   - Quality assessment
   - Run production analysis
   - Interpret results

**You'll learn**: Basic operations, quality checking, and how to analyze your data confidently.

---

### **Wet-Lab Scientists (No Coding Required)**

**Recommended path**:
1. [Tutorial 01: Getting Started](tutorial_01_getting_started.md) → Dashboard basics
2. [Tutorial 05: Interactive Dashboard](tutorial_05_dashboard.md) → Master web interface
3. [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md) → Analyze real data

**You'll learn**: Point-and-click RNA-seq analysis with no command line needed.

---

###  **Computational Biologists**

**Power user path**:
1. [Tutorial 01: Getting Started](tutorial_01_getting_started.md) → Quick overview
2. [Tutorial 02: Advanced Benchmarking](tutorial_02_benchmarking.md) → Compare all pipelines
3. [Tutorial 04: ML Recommendations](tutorial_04_ml_recommendations.md) → Understand ML deeply
4. [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md) → Advanced techniques

**You'll learn**: Full command-line control, ensemble strategies, and ML customization.

---

###  **For Publications & Critical Decisions**

**Rigorous analysis path**:
1. [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md) → Quality assessment
2. [Tutorial 02: Advanced Benchmarking](tutorial_02_benchmarking.md) → Systematic comparison
3. [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md) → Robust results
4. [Tutorial 04: ML Recommendations](tutorial_04_ml_recommendations.md) → Validate ML choices

**You'll learn**: Publication-ready workflows with validation and ensemble analysis.

---

##  Quick Start by Use Case

### "I have a count matrix and need DEGs quickly"

→ **[Tutorial 01](tutorial_01_getting_started.md)**: Part 1-3 (15 minutes)
```bash
raptor profile --counts data.csv --recommend
raptor run --pipeline RECOMMENDED --counts data.csv --groups ...
```

---

### "I need to compare multiple pipelines"

→ **[Tutorial 02](tutorial_02_benchmarking.md)**: Part 1-2 (3-4 hours)
```bash
raptor benchmark --counts data.csv --metadata meta.txt
```

---

### "I want to use the web interface"

→ **[Tutorial 05](tutorial_05_dashboard.md)**: Complete tutorial (1 hour)
```bash
raptor dashboard
# Open browser, upload data, click buttons!
```

---

### "I'm not sure which pipeline to trust"

→ **[Tutorial 06](tutorial_06_ensemble.md)**: Part 1-3 (2 hours)
```bash
raptor ensemble --counts data.csv --method weighted
```

---

### "My data has quality issues"

→ **[Tutorial 03](tutorial_03_real_data.md)**: Part 1 (30 minutes)
```bash
raptor assess-quality --counts data.csv --detailed
```

---

## 📖 Tutorial Details

### Tutorial 01: Getting Started ⭐ **START HERE**

**What you'll learn**:
- ML-powered pipeline recommendations
- Interactive dashboard basics
- Data profiling and quality checks
- Running your first analysis
- Understanding results

**Best for**: Everyone! This is the foundation.

**Key concepts**: Profiling, ML recommendations, basic workflows

[→ Start Tutorial 01](tutorial_01_getting_started.md)

---

### Tutorial 02: Advanced Benchmarking

**What you'll learn**:
- Compare all 8 pipelines systematically
- Concordance analysis
- ML-based pipeline ranking
- Resource monitoring
- Ensemble analysis introduction

**Best for**: Users who want comprehensive pipeline comparison

**Prerequisites**: Tutorial 01

**Key concepts**: Benchmarking, concordance, ensemble methods

[→ Start Tutorial 02](tutorial_02_benchmarking.md)

---

### Tutorial 03: Real Data Analysis

**What you'll learn**:
- Quality assessment workflows
- Handling complex experimental designs
- Batch effect detection and correction
- Production-ready analysis
- Publication outputs

**Best for**: Analyzing your own RNA-seq data

**Prerequisites**: Tutorial 01

**Key concepts**: Quality control, real-world workflows, troubleshooting

[→ Start Tutorial 03](tutorial_03_real_data.md)

---

### Tutorial 04: ML Recommendations (NEW in v2.1.0)

**What you'll learn**:
- How ML recommendations work
- Training custom models
- Fine-tuning for your organism/tissue
- Interpreting ML confidence scores
- When to trust (or override) ML

**Best for**: Advanced users wanting to understand or customize ML

**Prerequisites**: Tutorials 01-03

**Key concepts**: Machine learning, model training, customization

[→ Start Tutorial 04](tutorial_04_ml_recommendations.md)

---

### Tutorial 05: Interactive Dashboard (NEW in v2.1.0)

**What you'll learn**:
- Complete dashboard guide
- No-code data analysis
- Visual parameter tuning
- Real-time monitoring
- Sharing results with collaborators

**Best for**: Users who prefer graphical interfaces

**Prerequisites**: Tutorial 01

**Key concepts**: Web interface, visual analysis, collaboration

[→ Start Tutorial 05](tutorial_05_dashboard.md)

---

### Tutorial 06: Ensemble Methods (NEW in v2.1.0)

**What you'll learn**:
- Advanced ensemble strategies
- Consensus vs. weighted voting
- Rank-product methods
- Handling discordant results
- Publication-ready ensemble analysis

**Best for**: Critical projects, publications, challenging datasets

**Prerequisites**: Tutorials 01-02

**Key concepts**: Ensemble analysis, robustness, validation

[→ Start Tutorial 06](tutorial_06_ensemble.md)

---

##  Skills You'll Gain

### After Tutorial 01
- ✅ Run RAPTOR with ML recommendations
- ✅ Use the interactive dashboard
- ✅ Profile data characteristics
- ✅ Execute basic differential expression analysis
- ✅ Interpret standard outputs

### After Tutorial 02
- ✅ Benchmark multiple pipelines
- ✅ Analyze concordance
- ✅ Use ensemble methods
- ✅ Monitor computational resources
- ✅ Generate comparison reports

### After Tutorial 03
- ✅ Quality-assess real RNA-seq data
- ✅ Handle complex experimental designs
- ✅ Correct for batch effects
- ✅ Troubleshoot common issues
- ✅ Create publication-ready outputs

### After Tutorials 04-06
- ✅ Understand ML recommendations deeply
- ✅ Train custom ML models
- ✅ Master the web dashboard
- ✅ Apply advanced ensemble techniques
- ✅ Validate and optimize workflows

---

##  Tutorial Datasets

All tutorials use realistic test datasets:

### Simple Test Data (Tutorial 01)
- 5-10 genes
- 6 samples (3+3)
- Clear differential expression
- **Purpose**: Learn basics quickly

### Benchmark Data (Tutorial 02)
- ~50 genes
- 8 samples (4+4)
- Known true/false positives
- Batch effects included
- **Purpose**: Compare pipelines systematically

### Real-World Simulation (Tutorial 03)
- Full transcriptome (~60,000 genes)
- 12+ samples
- Complex design (batches, covariates)
- Quality issues to handle
- **Purpose**: Practice with realistic data

**Download all datasets**:
```bash
raptor download-tutorial-data --output ~/raptor_tutorials/
```

---

##  Tips for Success

### 1. Follow Tutorials in Order
Each tutorial builds on previous concepts. Skip ahead only if you're experienced.

### 2. Use Your Own Data Early
After Tutorial 01, try Tutorial 03 with your actual data for immediate value.

### 3. Don't Skip Quality Assessment
Tutorial 03 Part 1 is crucial - most analysis problems start with data quality.

### 4. Experiment with the Dashboard
Tutorial 05 makes RAPTOR accessible to everyone - explore freely!

### 5. Trust But Verify
ML recommendations are powerful, but Tutorial 02 shows you how to validate them.

### 6. Join the Community
Ask questions, share experiences:
- GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
- Email: ayehbolouki1988@gmail.com

---

##  Tutorial Setup

### Prerequisites

**All tutorials require**:
- RAPTOR v2.1.0 installed ([INSTALLATION.md](../INSTALLATION.md))
- Basic command-line knowledge (except Tutorial 05)
- RNA-seq data or tutorial datasets

**Recommended**:
- 8+ GB RAM
- 4+ CPU cores
- Web browser (for dashboard tutorials)
- Python 3.9+ (for advanced tutorials)

### Quick Setup Check

```bash
# Verify installation
raptor --version
# Should show: RAPTOR v2.1.0

# Check installed pipelines
raptor check-install

# Download tutorial data
mkdir -p ~/raptor_tutorials
cd ~/raptor_tutorials
raptor download-tutorial-data

# You're ready!
```

---

##  Additional Resources

### Core Documentation
- **[README.md](../README.md)**: RAPTOR overview
- **[INSTALLATION.md](../INSTALLATION.md)**: Setup guide
- **[PIPELINES.md](../PIPELINES.md)**: Pipeline descriptions
- **[FAQ.md](../FAQ.md)**: Common questions

### New v2.1.0 Features
- **[ML_GUIDE.md](../ML_GUIDE.md)**: ML recommendations explained
- **[ENSEMBLE_GUIDE.md](../ENSEMBLE_GUIDE.md)**: Ensemble methods
- **[DASHBOARD_GUIDE.md](../DASHBOARD_GUIDE.md)**: Web interface
- **[QUALITY_GUIDE.md](../QUALITY_GUIDE.md)**: Quality assessment
- **[RESOURCE_MONITOR_GUIDE.md](../RESOURCE_MONITOR_GUIDE.md)**: Resource optimization

### Technical Documentation
- **[API.md](../API.md)**: Python API reference
- **[BENCHMARKING.md](../BENCHMARKING.md)**: Benchmarking details
- **[TROUBLESHOOTING.md](../TROUBLESHOOTING.md)**: Problem solving

---

##  Tutorial Goals

By completing all tutorials, you will:

✅ **Confidently analyze** RNA-seq data with best practices  
✅ **Understand** when to trust ML vs. when to validate manually  
✅ **Apply** ensemble methods for robust, publication-ready results  
✅ **Optimize** analyses for your computational resources  
✅ **Troubleshoot** common RNA-seq analysis problems  
✅ **Generate** publication-quality figures and reports  
✅ **Customize** workflows for your specific needs  

---

##  Getting Help

### During Tutorials
- Each tutorial has a "Troubleshooting" section
- Check [TROUBLESHOOTING.md](../TROUBLESHOOTING.md)
- Review [FAQ.md](../FAQ.md)

### Still Stuck?
1. **GitHub Issues**: https://github.com/AyehBlk/RAPTOR/issues
   - Search existing issues first
   - Provide: RAPTOR version, command used, error message
   
2. **Email**: ayehbolouki1988@gmail.com
   - Include tutorial name and step number
   - Attach relevant files if possible

3. **Documentation**: 
   - Main docs folder: See [DOCS_README.md](../DOCS_README.md)
   - Feature guides: Check guides for your specific question

---

##  Tutorial Feedback

We'd love to hear from you:
- Found a bug? → [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- Tutorial unclear? → Email us with suggestions
- Have an idea? → Open a discussion on GitHub
- Success story? → Share it with us!

---

##  Ready to Start?

**First time with RAPTOR?**  
→ [Tutorial 01: Getting Started](tutorial_01_getting_started.md)

**Have some experience?**  
→ Choose your learning path above

**Need specific help?**  
→ Jump to the relevant tutorial

**Want the dashboard?**  
→ [Tutorial 05: Interactive Dashboard](tutorial_05_dashboard.md)

---

##  Tutorial Version History

### v2.1.0 (Current)
-  Added Tutorial 04: ML Recommendations
-  Added Tutorial 05: Interactive Dashboard  
-  Added Tutorial 06: Ensemble Methods
-  Updated Tutorial 01-03 for v2.1.0 features
-  Added quality assessment to all tutorials
-  Improved all examples and figures

### v2.0.0
- Initial tutorial release
- Tutorials 01-03 (Getting Started, Benchmarking, Real Data)

---

## 🦖 About RAPTOR

**RAPTOR** (RNA-seq Analysis Pipeline Testing and Optimization Resource) helps you:
- Choose the best RNA-seq pipeline for YOUR data
- Validate results with ensemble analysis
- Ensure reproducible, publication-ready research

**Mission**: *Making free science for everybody around the world* 🌍

---

**Happy Learning!🦖**

*Created by Ayeh Bolouki*  
*Last updated: November 2025*

---

## Quick Links

- 📘 [Tutorial 01: Getting Started](tutorial_01_getting_started.md)
- 📗 [Tutorial 02: Advanced Benchmarking](tutorial_02_benchmarking.md)
- 📙 [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md)
- 📕 [Tutorial 04: ML Recommendations](tutorial_04_ml_recommendations.md)
- 📔 [Tutorial 05: Interactive Dashboard](tutorial_05_dashboard.md)
- 📓 [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md)
- 🏠 [Main Documentation](../DOCS_README.md)
- 💻 [GitHub Repository](https://github.com/AyehBlk/RAPTOR)
