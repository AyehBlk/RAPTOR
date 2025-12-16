#  RAPTOR v2.1.1 Documentation Index

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

Welcome to the complete documentation for RAPTOR v2.1.1! This index will help you find exactly what you need.

---

## 🆕 What's New in v2.1.1

### Adaptive Threshold Optimizer (ATO)
A new feature for data-driven threshold selection in differential expression analysis.

**Quick Start:**
```python
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(de_results, goal='discovery')
print(f"Optimal thresholds: |logFC| > {result.logfc_threshold:.2f}")
```

**New Documentation:** [THRESHOLD_OPTIMIZER.md](#threshold-optimizer-guide)

---

##  Quick Navigation

###  Getting Started
- [Installation Guide](#installation-guide) - Set up RAPTOR
- [Quick Start](#quick-start) - Get running in 5 minutes
- [User Guide](#user-guide) - Complete tutorials
- [FAQ](#faq) - Frequently asked questions

###  Core Documentation
- [README](#readme) - Project overview
- [Configuration Guide](#configuration-guide) - All settings explained
- [API Reference](#api-reference) - Python API documentation

###  New in v2.1.1
- [🎯 Threshold Optimizer Guide](#threshold-optimizer-guide) - **NEW!** Data-driven thresholds

###  Advanced Features
- [ML Training Guide](#ml-training-guide) - Custom model training
- [Dashboard Guide](#dashboard-guide) - Interactive interface (updated with ATO!)
- [Cloud Deployment](#cloud-deployment) - AWS/GCP/Azure setup
- [Ensemble Analysis](#ensemble-analysis) - Multi-pipeline consensus

### 🆘 Help & Support
- [Troubleshooting](#troubleshooting) - Common issues & solutions
- [FAQ](#faq) - Quick answers
- [Contributing](#contributing) - Join the project
- [Changelog](#changelog) - What's new in v2.1.1

---

##  Documentation Files

### Installation Guide
**File:** `INSTALLATION.md` | **Size:** 6 KB  
**Purpose:** Complete installation instructions for all platforms

**Contents:**
- System requirements (minimum & recommended)
- Quick installation (pip, conda, Docker)
- Platform-specific instructions
- **ATO verification steps** (v2.1.1)
- Troubleshooting common issues

**Who needs this:**
- ✅ First-time users
- ✅ System administrators
- ✅ Anyone setting up RAPTOR

[**View INSTALLATION.md →**](INSTALLATION.md)

---

### Quick Start Guide
**File:** `QUICK_START.md` | **Size:** 10 KB  
**Purpose:** Get running in 5 minutes

**Contents:**
- Installation (2 minutes)
- First run options (Dashboard, CLI, Python API)
- **Threshold optimization quick example** (v2.1.1)
- Basic usage patterns
- Common use cases

**Who needs this:**
- ✅ Everyone! Start here
- ✅ New users
- ✅ Anyone wanting quick results

[**View QUICK_START.md →**](QUICK_START.md)

---

### 🎯 Threshold Optimizer Guide (NEW!)
**File:** `THRESHOLD_OPTIMIZER.md` | **Size:** 20 KB  
**Purpose:** Comprehensive guide to the Adaptive Threshold Optimizer

**Contents:**
- Why use ATO (vs arbitrary thresholds)
- Analysis goals (discovery, balanced, validation)
- P-value adjustment methods (BH, BY, Storey, Holm, Bonferroni)
- LogFC optimization methods (MAD, mixture, power, percentile)
- π₀ estimation theory and methods
- Python API usage
- Dashboard usage
- Publication methods text
- Best practices

**Who needs this:**
- ✅ Anyone doing differential expression analysis
- ✅ Researchers preparing publications
- ✅ Users wanting data-driven thresholds

**Read if:**
- You're tired of arbitrary thresholds
- Reviewers ask "why |logFC| > 1?"
- You want publication-ready methods text
- You want to understand your data better

[**View THRESHOLD_OPTIMIZER.md →**](THRESHOLD_OPTIMIZER.md)

---

### Dashboard Guide
**File:** `DASHBOARD.md` | **Size:** 28 KB  
**Purpose:** Complete guide to the interactive web interface

**Contents:**
- Quick start
- All dashboard features
- Data upload
- Data profiling
- ML recommendations
- **🎯 Threshold Optimizer page** (v2.1.1)
- Results visualization
- Multi-user setup
- Troubleshooting

**Who needs this:**
- ✅ Users who prefer visual interfaces
- ✅ Teams setting up shared dashboards
- ✅ Anyone avoiding command-line

[**View DASHBOARD.md →**](DASHBOARD.md)

---

### Changelog
**File:** `CHANGELOG.md` | **Size:** 14 KB  
**Purpose:** Complete history of changes and updates

**Contents:**
- **v2.1.1 release notes** (current)
  - Adaptive Threshold Optimizer
  - Dashboard updates
  - Configuration updates
- v2.1.0 release notes
- v2.0.0 release notes
- Migration guides

**Who needs this:**
- ✅ Users updating from previous versions
- ✅ Anyone wanting to know what's new
- ✅ Contributors tracking changes

[**View CHANGELOG.md →**](CHANGELOG.md)

---

### Complete README
**File:** `COMPLETE_README.md` | **Size:** 15 KB  
**Purpose:** Master guide to all RAPTOR features

**Contents:**
- All features overview
- **v2.1.1 ATO highlights**
- Quick start options
- Performance comparisons
- Documentation links
- Troubleshooting
- Roadmap

**Who needs this:**
- ✅ New users wanting overview
- ✅ Anyone evaluating RAPTOR
- ✅ Documentation starting point

[**View COMPLETE_README.md →**](COMPLETE_README.md)

---

### Troubleshooting Guide
**File:** `TROUBLESHOOTING.md` | **Size:** 17 KB  
**Purpose:** Diagnose and fix common problems

**Contents:**
- Installation issues
- Configuration problems
- ML recommendation issues
- **Threshold Optimizer issues** (v2.1.1)
- Dashboard problems
- Pipeline execution errors
- Getting help

**Who needs this:**
- ✅ Users encountering errors
- ✅ Anyone debugging issues

[**View TROUBLESHOOTING.md →**](TROUBLESHOOTING.md)

---

### FAQ
**File:** `FAQ.md` | **Size:** 19 KB  
**Purpose:** Quick answers to common questions

**Contents:**
- General questions
- Installation & setup
- ML recommendations
- **Threshold selection** (v2.1.1)
- Dashboard usage
- Pipeline selection
- Cloud computing

**Who needs this:**
- ✅ Everyone! Quick answers
- ✅ New users

[**View FAQ.md →**](FAQ.md)

---

##  Documentation by User Type

###  Complete Beginners

**Start here:**
1. [FAQ](#faq) - Get familiar with RAPTOR
2. [Installation Guide](#installation-guide) - Set up RAPTOR
3. [Quick Start](#quick-start) - Run first analysis
4. [Dashboard Guide](#dashboard-guide) - Use visual interface

**Common questions:**
- "What is RAPTOR?" → [FAQ](#faq)
- "How do I install it?" → [Installation Guide](#installation-guide)
- "What thresholds should I use?" → [Threshold Optimizer Guide](#threshold-optimizer-guide)

---

###  Researchers

**Essential reading:**
1. [Quick Start](#quick-start) - Get running fast
2. [Threshold Optimizer Guide](#threshold-optimizer-guide) - **Essential for publications!**
3. [Dashboard Guide](#dashboard-guide) - Visual workflows
4. [FAQ](#faq) - Quick answers

**Common tasks:**
- Run basic analysis → [Quick Start](#quick-start)
- Get optimal thresholds → [Threshold Optimizer Guide](#threshold-optimizer-guide)
- Generate methods text → [Threshold Optimizer Guide](#threshold-optimizer-guide)

---

###  Bioinformaticians

**Focus on:**
1. [API Reference](#api-reference) - Python API
2. [Threshold Optimizer Guide](#threshold-optimizer-guide) - ATO API
3. [Configuration Guide](#configuration-guide) - Advanced settings
4. [ML Training Guide](#ml-training-guide) - Custom models

**Common needs:**
- Integrate ATO with existing pipeline → [Threshold Optimizer Guide](#threshold-optimizer-guide)
- Custom analysis → [API Reference](#api-reference)
- Optimize parameters → [Configuration Guide](#configuration-guide)

---

###  Core Facility Managers

**Important docs:**
1. [Installation Guide](#installation-guide) - Setup on servers
2. [Dashboard Guide](#dashboard-guide) - Multi-user interface
3. [Cloud Deployment](#cloud-deployment) - Scale for many users
4. [Troubleshooting](#troubleshooting) - Support users

---

##  Find Documentation by Topic

### Threshold Optimization (NEW in v2.1.1)
- [Threshold Optimizer Guide](#threshold-optimizer-guide) - Complete guide
- [Dashboard Guide: ATO Page](#dashboard-guide) - Visual interface
- [Quick Start: ATO Example](#quick-start) - Quick example
- [FAQ: Threshold Selection](#faq) - Common questions

### Installation & Setup
- [Installation Guide](#installation-guide) - Complete setup
- [Quick Start](#quick-start) - Fast setup
- [Troubleshooting](#troubleshooting) - Fix problems

### Machine Learning
- [ML Training Guide](#ml-training-guide) - Custom models
- [FAQ: ML](#faq) - Common questions
- [Dashboard Guide](#dashboard-guide) - ML in dashboard

### Dashboard
- [Dashboard Guide](#dashboard-guide) - Full guide
- [FAQ: Dashboard](#faq) - Common questions
- [Troubleshooting](#troubleshooting) - Fix issues

### Configuration
- [Configuration Guide](#configuration-guide) - All options
- [FAQ](#faq) - Common questions
- `config.yaml` examples - In configs directory

### Pipelines
- [Pipeline Guide](#pipeline-guide) - Pipeline details
- [FAQ: Pipeline Selection](#faq) - Which to use?
- [Ensemble Analysis](#ensemble-analysis) - Combine pipelines

### Cloud Computing
- [Cloud Deployment Guide](#cloud-deployment) - Full setup
- [FAQ: Cloud](#faq) - Costs, safety
- [Installation: Cloud](#installation-guide) - Cloud installation

---

##  Documentation Statistics

| Type | Files | Total Size | Purpose |
|------|-------|------------|---------|
| Getting Started | 2 | 16 KB | Installation, quick start |
| User Guides | 4 | 75 KB | Dashboard, ATO, ML, ensemble |
| Reference | 3 | 45 KB | API, config, pipelines |
| Support | 3 | 50 KB | FAQ, troubleshooting, contributing |
| **Total** | **12** | **~186 KB** | Complete documentation |

---

##  Quick Start Paths

### Path 1: "I just want optimal thresholds"
1. [Installation Guide](#installation-guide) - Install RAPTOR
2. [Threshold Optimizer Guide](#threshold-optimizer-guide) - Use ATO
3. Done! Copy methods text to paper

**Time:** 15 minutes

---

### Path 2: "I want visual interface"
1. [Installation Guide](#installation-guide) - Install RAPTOR
2. [Dashboard Guide](#dashboard-guide) - Launch dashboard
3. Use 🎯 Threshold Optimizer page

**Time:** 20 minutes

---

### Path 3: "I need ML recommendations"
1. [Installation Guide](#installation-guide) - Install RAPTOR
2. [Quick Start](#quick-start) - ML overview
3. [Dashboard Guide](#dashboard-guide) - Use ML in dashboard

**Time:** 30 minutes

---

### Path 4: "I need everything"
1. [Installation Guide](#installation-guide) - Full setup
2. [Quick Start](#quick-start) - Overview
3. [Dashboard Guide](#dashboard-guide) - Visual interface
4. [Threshold Optimizer Guide](#threshold-optimizer-guide) - Optimal thresholds
5. [ML Training Guide](#ml-training-guide) - Custom models

**Time:** 2-3 hours

---

## ✅ Documentation Checklist

Before running your first analysis:

- [ ] Read [Installation Guide](#installation-guide)
- [ ] Install RAPTOR successfully
- [ ] Verify ATO: `python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅')"`
- [ ] Read [Quick Start](#quick-start)
- [ ] Try the dashboard: `raptor dashboard`
- [ ] Try Threshold Optimizer (🎯 page in dashboard)
- [ ] Bookmark this index

---

##  All Documentation Files

### Core Documentation (v2.1.1)
1. ✅ **INSTALLATION.md** - Installation guide
2. ✅ **QUICK_START.md** - Get running fast
3. ✅ **COMPLETE_README.md** - Master overview
4. ✅ **CHANGELOG.md** - Version history
5. ✅ **DASHBOARD.md** - Dashboard guide (with ATO)
6. ✅ **DOCUMENTATION_INDEX.md** - This file
7. ✅ **THRESHOLD_OPTIMIZER.md** - 🆕 ATO guide

### Additional Documentation
- FAQ.md - Frequently asked questions
- TROUBLESHOOTING.md - Problem solving
- CONTRIBUTING.md - Contribution guide
- API.md - Python API reference
- ML_TRAINING.md - ML guide
- ENSEMBLE.md - Ensemble analysis
- CLOUD_DEPLOYMENT.md - Cloud guide
- LICENSE - Legal terms

---

##  Online Resources

### Official
- **GitHub**: https://github.com/AyehBlk/RAPTOR
- **PyPI**: https://pypi.org/project/raptor-rnaseq/
- **Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Discussions**: https://github.com/AyehBlk/RAPTOR/discussions

### Community
- **Email**: ayehbolouki1988@gmail.com

---

## 💙 Thank You!

We hope this documentation helps you make the most of RAPTOR v2.1.1!

**Questions? Suggestions?**
- Email: ayehbolouki1988@gmail.com
- GitHub: https://github.com/AyehBlk/RAPTOR

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.1  
**License:** MIT

---

*"Good documentation is like a good friend - always there when you need it!"* 
