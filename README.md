<p align="center">
  <img src="https://img.shields.io/badge/🦖_RAPTOR-v2.2.0-brightgreen?style=for-the-badge" alt="RAPTOR v2.2.0"/>
</p>

<h1 align="center">RAPTOR</h1>
<h3 align="center">RNA-seq Analysis Pipeline Testing and Optimization Resource</h3>

<p align="center">
  <strong>Making free science for everybody around the world 🌍</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/raptor-rnaseq/"><img src="https://img.shields.io/pypi/v/raptor-rnaseq.svg?style=flat&logo=pypi&logoColor=white" alt="PyPI version"/></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.8--3.12-blue.svg" alt="Python 3.8-3.12"/></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="MIT License"/></a>
  <a href="https://doi.org/10.5281/zenodo.17607161"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17607161.svg" alt="DOI"/></a>
  <a href="https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.0"><img src="https://img.shields.io/badge/Release-v2.2.0-orange.svg" alt="Release v2.2.0"/></a>
</p>

<p align="center">
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-features">Features</a> •
  <a href="#-installation">Installation</a> •
  <a href="#-architecture">Architecture</a> •
  <a href="#-documentation">Documentation</a> •
  <a href="#-pipelines">Pipelines</a> •
  <a href="#-citation">Citation</a>
</p>

---

## 🦖 What is RAPTOR?

**RAPTOR** is a comprehensive framework for RNA-seq analysis that makes sophisticated differential expression workflows accessible to everyone. Stop wondering which pipeline to use or what thresholds to set—RAPTOR provides **ML-powered recommendations** and **ensemble methods** for robust, reproducible results.

### Why RAPTOR?

| Challenge | RAPTOR Solution |
|-----------|-----------------|
| Which pipeline should I use? | ✅ **ML recommendations** based on 32 dataset features |
| Which DE method (DESeq2/edgeR/limma)? | ✅ **Ensemble analysis** combines all methods |
| What thresholds should I use? | ✅ **4 optimization methods** for data-driven cutoffs |
| Is my data quality good enough? | ✅ **6 outlier detection methods** with consensus |
| How do I know results are reliable? | ✅ **Ensemble consensus** with direction checking |
| What if methods disagree? | ✅ **Brown's method** accounts for correlation |

---

## ✨ Features

<table>
<tr>
<td width="50%">

### 🎯 Ensemble Analysis (NEW!)
- 5 statistical combination methods
- Fisher's, Brown's, RRA, Voting, Weighted
- Direction consistency checking
- Meta-analysis fold changes
- 33% fewer false positives

### ⚙️ Parameter Optimization (NEW!)
- 4 validated optimization methods
- Ground truth, FDR control, Stability, Reproducibility
- Automated threshold selection
- Performance metrics tracking
- Publication-ready results

### 📊 32-Feature Data Profiling
- BCV (Biological Coefficient of Variation)
- Sample characteristics & balance
- Dispersion patterns
- Sparsity analysis
- ML-ready feature vectors

</td>
<td width="50%">

### 🤖 ML-Powered Recommendations
- Random Forest classifier
- 32-feature profiling
- Confidence scoring
- Alternative suggestions
- Feature importance analysis

### 🔬 6 Production Pipelines
- Salmon ⭐ (recommended)
- Kallisto (fastest)
- STAR + featureCounts
- STAR + RSEM
- STAR + Salmon (unique!)
- HISAT2 + featureCounts

### 📈 Quality Assessment
- 6 outlier detection methods
- Consensus-based reporting
- Batch effect detection
- Actionable recommendations

</td>
</tr>
</table>

### 🎨 Interactive Dashboard

- Web-based interface (no coding!)
- Real-time visualizations
- Drag-and-drop data upload
- One-click ensemble analysis
- Export publication-ready reports

---

## 🚀 Quick Start

### Option 1: Interactive Dashboard (Recommended)

```bash
# Install
pip install raptor-rnaseq

# Launch dashboard
streamlit run raptor/dashboard/app.py

# Opens at http://localhost:8501
# Upload data → Profile → Get recommendation → Run ensemble → Done!
```

### Option 2: Command Line

```bash
# 1. Quality check
raptor qc --counts counts.csv --metadata metadata.csv

# 2. Profile your data
raptor profile --counts counts.csv --metadata metadata.csv --group-column condition

# 3. Get ML recommendation
raptor recommend --profile profile.json --method ml

# 4. Import DE results from different methods
raptor import-de --input deseq2.csv --method deseq2
raptor import-de --input edger.csv --method edger
raptor import-de --input limma.csv --method limma

# 5. Optimize thresholds (NEW!)
raptor optimize --de-result de_results.csv --method fdr-control --fdr-target 0.05

# 6. Ensemble analysis - combine all methods (NEW!)
raptor ensemble-compare --deseq2 deseq2.csv --edger edger.csv --limma limma.csv
```

### Option 3: Python API

```python
from raptor import (
    quick_quality_check,
    profile_data_quick,
    recommend_pipeline,
    optimize_with_fdr_control,
    ensemble_brown
)

# 1. Quality check
qc_report = quick_quality_check('counts.csv', 'metadata.csv')
print(f"Outliers: {qc_report.outliers}")

# 2. Profile data (32 features extracted)
profile = profile_data_quick('counts.csv', 'metadata.csv', group_column='condition')
print(f"BCV: {profile.bcv:.3f} ({profile.bcv_category})")

# 3. Get ML recommendation
recommendation = recommend_pipeline(profile_file='profile.json', method='ml')
print(f"Recommended: {recommendation.pipeline_name} (confidence: {recommendation.confidence:.2f})")

# 4. After running DE analysis, optimize thresholds (NEW!)
result = optimize_with_fdr_control(de_result, fdr_target=0.05)
print(f"Optimal thresholds: {result.optimal_threshold}")

# 5. Ensemble analysis - combine DESeq2, edgeR, limma (NEW!)
consensus = ensemble_brown({
    'deseq2': deseq2_result,
    'edger': edger_result,
    'limma': limma_result
})
print(f"Consensus DE genes: {len(consensus.consensus_genes)}")
```

---

## 📦 Installation

### Requirements

- **Python**: 3.8 - 3.12
- **R**: 4.0+ (optional, for Module 6 DE analysis)
- **RAM**: 4GB minimum (16GB recommended for pipelines)
- **Disk**: 500MB (Python package) / 5-8GB (with bioinformatics tools)

### Install from PyPI (Recommended)

```bash
# Basic installation
pip install raptor-rnaseq

# With dashboard support
pip install raptor-rnaseq[dashboard]

# With all features
pip install raptor-rnaseq[all]

# Development installation
pip install raptor-rnaseq[dev]
```

### Conda Installation

**Core environment** (Python only, ~500MB, 5-10 min):
```bash
conda env create -f environment.yml
conda activate raptor
```

**Full environment** (with STAR, Salmon, Kallisto, R, ~5-8GB, 30-60 min):
```bash
conda env create -f environment-full.yml
conda activate raptor-full
```

See [docs/CONDA_ENVIRONMENTS.md](docs/CONDA_ENVIRONMENTS.md) for detailed comparison.

### Install from Source

```bash
# Clone repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Install in editable mode
pip install -e .

# Or with development tools
pip install -e .[dev]

# Verify installation
raptor --version
pytest tests/
```

---

## 🏗️ Architecture

RAPTOR is organized into **9 modules** spanning 4 analysis stages:

```
┌─────────────────────────────────────────────────────────────┐
│                    RAPTOR v2.2.0                            │
│         RNA-seq Analysis Pipeline Framework                 │
└─────────────────────────────────────────────────────────────┘

Stage 1: Data Preparation & QC
├── Module 1: Quick Quantification (Salmon/Kallisto)
├── Module 2: Quality Assessment (6 outlier methods)
└── Module 3: Data Profiling (32 features)

Stage 2: Pipeline Selection
├── Module 4: ML Recommender (Random Forest)
└── Module 5: Production Pipelines (6 methods)
         ├── Salmon ⭐ (recommended)
         ├── Kallisto (fastest)
         ├── STAR + featureCounts
         ├── STAR + RSEM
         ├── STAR + Salmon (unique: BAM + bootstraps)
         └── HISAT2 + featureCounts

Stage 3: Differential Expression
├── Module 6: DE Analysis (R: DESeq2, edgeR, limma)
└── Module 7: DE Import (standardize any format)

Stage 4: Advanced Analysis ⭐ NEW in v2.2.0
├── Module 8: Parameter Optimization (4 methods)
│    ├── Ground Truth Optimization
│    ├── FDR Control Optimization
│    ├── Stability Optimization
│    └── Reproducibility Optimization
└── Module 9: Ensemble Analysis (5 methods)
     ├── Fisher's Method
     ├── Brown's Method
     ├── Robust Rank Aggregation
     ├── Voting Consensus
     └── Weighted Ensemble
```

---

## 🧬 Pipelines

RAPTOR supports **6 production RNA-seq quantification pipelines**:

| Pipeline | Memory | Time | Produces | Best For | Recommended |
|----------|--------|------|----------|----------|-------------|
| **salmon** | 8 GB | 10-20 min | genes + isoforms + bootstraps | Standard DE analysis | ⭐ **YES** |
| **kallisto** | 4 GB | 5-10 min | genes + isoforms + bootstraps | Speed priority | ✓ |
| **star_featurecounts** | 32 GB | 40-70 min | BAM + genes | Gene-level publication | ✓ |
| **star_rsem** | 32 GB | 60-120 min | BAM + genes + isoforms | Isoform analysis | ✓ |
| **star_salmon** | 32 GB | 50-90 min | BAM + genes + isoforms + bootstraps | Unique: BAM + bootstraps | ✓ |
| **hisat2_featurecounts** | 16 GB | 30-60 min | BAM + genes | Low memory systems | ✓ |

⭐ **Salmon is recommended** for most use cases due to optimal speed/accuracy balance and bootstrap support.

### Pipeline Features

**All pipelines support:**
- ✅ Paired-end and single-end reads
- ✅ Automatic parameter optimization
- ✅ QC report generation
- ✅ Multi-threading
- ✅ Sample sheet-based workflows

**Pipeline selection:**
```bash
# List available pipelines
raptor pipeline list

# Get detailed info
raptor pipeline run --name salmon --help

# Run with ML recommendation
raptor recommend --profile profile.json --method ml
# Recommended: salmon (confidence: 0.89)

raptor pipeline run --name salmon --samples samples.csv --index salmon_index/
```

---

## 🗂️ Repository Structure

```
RAPTOR/
├── raptor/                         # Core Python package
│   ├── __init__.py                 # Package initialization (v2.2.0)
│   ├── cli.py                      # Command-line interface (11 commands)
│   ├── quality_assessment.py       # Module 2: QC (6 methods)
│   ├── profiler.py                 # Module 3: Profiling (32 features)
│   ├── recommender.py              # Module 4: Rule-based
│   ├── ml_recommender.py           # Module 4: ML-based
│   ├── de_import.py                # Module 7: DE import
│   ├── parameter_optimization.py   # Module 8: Optimization ⭐ NEW
│   ├── ensemble.py                 # Module 9: Ensemble ⭐ NEW
│   ├── simulation.py               # Simulation tools
│   │
│   ├── pipelines/                  # Module 5: Production pipelines
│   │   ├── base.py
│   │   ├── salmon/
│   │   ├── kallisto/
│   │   ├── star_featurecounts/
│   │   ├── star_rsem/
│   │   ├── star_salmon/
│   │   └── hisat2_featurecounts/
│   │
│   ├── external_modules/           # Module 6: R integration
│   │   └── module6_de_analysis/
│   │       └── r_scripts/          # DESeq2, edgeR, limma
│   │
│   ├── dashboard/                  # Interactive Streamlit app
│   │   ├── app.py
│   │   ├── pages/                  # 9 dashboard pages
│   │   ├── components/
│   │   └── utils/
│   │
│   └── utils/                      # Utilities
│       ├── validation.py
│       ├── errors.py
│       └── sample_sheet.py
│
├── docs/                           # Documentation
│   ├── MODULE_1_Quick_Quantification.md
│   ├── MODULE_2_Quality_Assessment.md
│   ├── MODULE_3_Data_Profiling.md
│   ├── MODULE_3_QUICK_REFERENCE.md
│   ├── MODULE_4_Pipeline_Recommender.md
│   ├── MODULE_7_DE_Import.md
│   ├── MODULE_8_Parameter_Optimization.md      ⭐ NEW
│   ├── MODULE_9_Ensemble_Analysis.md           ⭐ NEW
│   ├── CONDA_ENVIRONMENTS.md
│   ├── RAPTOR_QUICK_REFERENCE.md               # Cheat sheet
│   └── RAPTOR_API_DOCUMENTATION.md             # Python API
│
├── examples/                       # Example scripts
│   ├── 02_quality_assessment.py
│   ├── 03_data_profiler.py
│   ├── 04_recommender.py
│   ├── 07_DE_Import.py
│   ├── 08_Parameter_Optimization.py            ⭐ NEW
│   └── 09_Ensemble_Analysis.py                 ⭐ NEW
│
├── tests/                          # Test suite (85%+ coverage)
│   ├── test_profiler.py
│   ├── test_quality_assessment.py
│   ├── test_parameter_optimization.py          ⭐ NEW
│   ├── test_ensemble.py                        ⭐ NEW
│   └── ...
│
├── templates/                      # Sample sheets
│   ├── sample_sheet_paired.csv
│   └── sample_sheet_single.csv
│
├── .github/                        # GitHub templates
│   └── ISSUE_TEMPLATE/
│
├── setup.py                        # Package setup
├── requirements.txt                # Python dependencies
├── environment.yml                 # Conda environment (core)
├── environment-full.yml            # Conda environment (complete)
├── CITATION.cff                    # Citation metadata
├── CHANGELOG.md                    # Version history
├── CONTRIBUTING.md                 # Contribution guidelines
└── LICENSE                         # MIT License
```

---

## 📚 Documentation

### Getting Started
| Document | Description |
|----------|-------------|
| [Quick Start](#-quick-start) | 5-minute quick start guide |
| [Installation](#-installation) | Detailed installation instructions |
| [CONDA_ENVIRONMENTS.md](docs/CONDA_ENVIRONMENTS.md) | Conda setup (core vs full) |

### Core Features (v2.2.0)
| Document | Description |
|----------|-------------|
| [MODULE_2_Quality_Assessment.md](docs/MODULE_2_Quality_Assessment.md) | QC with 6 outlier methods |
| [MODULE_3_Data_Profiling.md](docs/MODULE_3_Data_Profiling.md) | 32-feature profiling |
| [MODULE_3_QUICK_REFERENCE.md](docs/MODULE_3_QUICK_REFERENCE.md) | Profiling cheat sheet |
| [MODULE_4_Pipeline_Recommender.md](docs/MODULE_4_Pipeline_Recommender.md) | ML recommendations |
| [MODULE_7_DE_Import.md](docs/MODULE_7_DE_Import.md) | Import & standardize DE results |
| [MODULE_8_Parameter_Optimization.md](docs/MODULE_8_Parameter_Optimization.md) | ⭐ 4 optimization methods |
| [MODULE_9_Ensemble_Analysis.md](docs/MODULE_9_Ensemble_Analysis.md) | ⭐ 5 ensemble methods |

### Reference
| Document | Description |
|----------|-------------|
| [RAPTOR_QUICK_REFERENCE.md](docs/RAPTOR_QUICK_REFERENCE.md) | Command cheat sheet |
| [RAPTOR_API_DOCUMENTATION.md](docs/RAPTOR_API_DOCUMENTATION.md) | Complete Python API |
| [examples/](examples/) | Example scripts for all modules |
| [CHANGELOG.md](CHANGELOG.md) | Version history |

---

## 💡 Usage Examples

### Example 1: Complete Workflow (v2.2.0)

```python
from raptor import (
    quick_quality_check,
    profile_data_quick,
    recommend_pipeline,
    import_deseq2,
    import_edger,
    import_limma,
    optimize_with_fdr_control,
    ensemble_brown
)

# 1. Quality Check
print("Step 1: Quality Assessment...")
qc_report = quick_quality_check('counts.csv', 'metadata.csv')
if len(qc_report.outliers) > 0:
    print(f"⚠️ Warning: {len(qc_report.outliers)} outliers detected")
else:
    print("✅ No outliers detected")

# 2. Profile Data (32 features)
print("\nStep 2: Data Profiling...")
profile = profile_data_quick('counts.csv', 'metadata.csv', group_column='condition')
print(f"  BCV: {profile.bcv:.3f} ({profile.bcv_category})")
print(f"  Sample size: {profile.n_samples}")

# 3. Get ML Recommendation
print("\nStep 3: ML Recommendation...")
rec = recommend_pipeline(profile_file='results/profile/data_profile.json', method='ml')
print(f"  Recommended: {rec.pipeline_name} (confidence: {rec.confidence:.2f})")

# 4. [Run recommended pipeline, then DE analysis in R]

# 5. Import DE Results
print("\nStep 4: Import DE Results...")
deseq2 = import_deseq2('deseq2_results.csv')
edger = import_edger('edger_results.csv')
limma = import_limma('limma_results.csv')

# 6. Optimize Thresholds (NEW!)
print("\nStep 5: Optimize Thresholds...")
opt_result = optimize_with_fdr_control(deseq2, fdr_target=0.05)
print(f"  Optimal FDR: {opt_result.optimal_threshold['padj']:.3f}")
print(f"  Optimal |logFC|: {opt_result.optimal_threshold['lfc']:.3f}")

# 7. Ensemble Analysis (NEW!)
print("\nStep 6: Ensemble Analysis (Brown's Method)...")
consensus = ensemble_brown({
    'deseq2': deseq2,
    'edger': edger,
    'limma': limma
})
print(f"  Consensus genes: {len(consensus.consensus_genes)}")
print(f"  Direction consistency: {consensus.direction_consistency.mean():.1%}")

# 8. Export Results
consensus.to_csv('consensus_genes.csv')
print("\n✅ Analysis complete!")
```

### Example 2: Ensemble Analysis Only

```python
from raptor import import_de_result, ensemble_fisher, ensemble_brown, ensemble_rra

# Import results from different tools
deseq2 = import_de_result('deseq2_results.csv', method='deseq2')
edger = import_de_result('edger_results.csv', method='edger')
limma = import_de_result('limma_results.csv', method='limma')

# Try multiple ensemble methods
results = {}

# Fisher's Method (classic)
results['fisher'] = ensemble_fisher({'deseq2': deseq2, 'edger': edger, 'limma': limma})

# Brown's Method (recommended - accounts for correlation)
results['brown'] = ensemble_brown({'deseq2': deseq2, 'edger': edger, 'limma': limma})

# Robust Rank Aggregation
results['rra'] = ensemble_rra({'deseq2': deseq2, 'edger': edger, 'limma': limma})

# Compare results
for method, result in results.items():
    print(f"{method}: {len(result.consensus_genes)} consensus genes")

# Use Brown's method (best for correlated methods)
final_result = results['brown']
final_result.to_csv('final_consensus.csv')
```

### Example 3: CLI Workflow

```bash
#!/bin/bash
# Complete RAPTOR v2.2.0 workflow using CLI

# Step 1: QC
raptor qc --counts counts.csv --metadata metadata.csv --output qc_results/

# Step 2: Profile
raptor profile --counts counts.csv --metadata metadata.csv --group-column condition

# Step 3: Recommend
raptor recommend --profile profile.json --method ml

# Step 4: Import DE results
raptor import-de --input deseq2_results.csv --method deseq2 --output imported/
raptor import-de --input edger_results.csv --method edger --output imported/
raptor import-de --input limma_results.csv --method limma --output imported/

# Step 5: Optimize thresholds (NEW!)
raptor optimize --de-result imported/deseq2.csv --method fdr-control --fdr-target 0.05

# Step 6: Ensemble analysis (NEW!)
raptor ensemble-compare \
    --deseq2 imported/deseq2.csv \
    --edger imported/edger.csv \
    --limma imported/limma.csv \
    --output ensemble_results/

echo "✅ Complete! Check ensemble_results/ for consensus genes."
```

---

## 📊 Performance

### Module Performance

| Module | Time | Memory | Key Output |
|--------|------|--------|------------|
| **Module 2: QC** | 1-5 min | 4 GB | 6 methods consensus |
| **Module 3: Profiler** | 1-2 min | 4 GB | 32 features + BCV |
| **Module 4: Recommender** | <10 sec | <1 GB | ML recommendation |
| **Module 8: Optimization** | 5-30 min | 4 GB | Optimal thresholds |
| **Module 9: Ensemble** | <1 min | 2 GB | Consensus genes |

### Ensemble Analysis Benefits

| Metric | Single Method | Ensemble (Brown's) |
|--------|---------------|-------------------|
| False Positive Rate | Higher | **33% lower** |
| Reproducibility | Variable | **Higher** |
| Confidence | Method-specific | **Consensus-based** |
| Publication Impact | Good | **Better** |

---

## 🤝 Contributing

We welcome contributions! RAPTOR is open-source and aims to make free science accessible to everyone.

```bash
# Fork and clone
git clone https://github.com/YOUR_USERNAME/RAPTOR.git
cd RAPTOR

# Create feature branch
git checkout -b feature/amazing-feature

# Make changes and test
pytest tests/

# Submit pull request
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

### Ways to Contribute

- 🐛 **Report bugs** via [Issues](https://github.com/AyehBlk/RAPTOR/issues)
- ✨ **Request features**
- 📝 **Improve documentation**
- 🔧 **Submit pull requests**
- 💡 **Share use cases** and feedback
- ⭐ **Star the repository**

---

## 📖 Citation

If you use RAPTOR in your research, please cite:

```bibtex
@software{bolouki2026raptor,
  author       = {Bolouki, Ayeh},
  title        = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year         = {2026},
  version      = {2.2.0},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17607161},
  url          = {https://github.com/AyehBlk/RAPTOR}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17607161.svg)](https://doi.org/10.5281/zenodo.17607161)

---

## 📜 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

```
MIT License
Copyright (c) 2026 Ayeh Bolouki
```

---

## 📧 Contact

**Ayeh Bolouki**

- 🏛️ **GIGA, University of Liège, Belgium**
- 📧 **Email:** ayehbolouki1988@gmail.com
- 🐙 **GitHub:** [@AyehBlk](https://github.com/AyehBlk)
- 🔬 **Research:** Computational Biology, Bioinformatics, Multi-omics Analysis

### Support

- **📖 Documentation:** [docs/](docs/)
- **🐛 Issues:** [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- **💬 Discussions:** [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)
- **📧 Email:** ayehbolouki1988@gmail.com

---

## 🙏 Acknowledgments

RAPTOR builds on the excellent work of the RNA-seq community:

- **Bioconductor** community for the R package ecosystem
- **DESeq2** (Love et al., 2014) - Differential expression analysis
- **edgeR** (Robinson et al., 2010) - Empirical analysis of DGE
- **limma** (Ritchie et al., 2015) - Linear models for microarray and RNA-seq
- **Salmon** (Patro et al., 2017) - Wicked-fast transcript quantification
- **Kallisto** (Bray et al., 2016) - Near-optimal probabilistic RNA-seq quantification
- **STAR** (Dobin et al., 2013) - Ultrafast universal RNA-seq aligner
- All users who provided feedback and suggestions

---

<p align="center">
  <strong>⭐ Star this repository if you find RAPTOR useful!</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/github/stars/AyehBlk/RAPTOR?style=social" alt="GitHub Stars"/>
  <img src="https://img.shields.io/github/forks/AyehBlk/RAPTOR?style=social" alt="GitHub Forks"/>
</p>

<p align="center">
  <em>RAPTOR v2.2.0 - Making pipeline selection evidence-based, not guesswork 🦖</em>
</p>

<p align="center">
  <em>Making free science for everybody around the world 🌍</em>
</p>
