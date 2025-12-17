<p align="center">
  <img src="https://img.shields.io/badge/🦖_RAPTOR-v2.1.1-brightgreen?style=for-the-badge" alt="RAPTOR v2.1.1"/>
</p>

<h1 align="center">RAPTOR</h1>
<h3 align="center">RNA-seq Analysis Pipeline Testing and Optimization Resource</h3>

<p align="center">
  <strong>Making free science for everybody around the world 🌍</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/raptor-rnaseq/"><img src="https://img.shields.io/pypi/v/raptor-rnaseq.svg?style=flat&logo=pypi&logoColor=white" alt="PyPI version"/></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.8+-blue.svg" alt="Python 3.8+"/></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="MIT License"/></a>
  <a href="https://doi.org/10.5281/zenodo.17607161"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17607161.svg" alt="DOI"/></a>
  <a href="https://github.com/AyehBlk/RAPTOR/releases/tag/v2.1.1"><img src="https://img.shields.io/badge/Release-v2.1.1-orange.svg" alt="Release v2.1.1"/></a>
</p>

<p align="center">
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-features">Features</a> •
  <a href="#-installation">Installation</a> •
  <a href="#-documentation">Documentation</a> •
  <a href="#-pipelines">Pipelines</a> •
  <a href="#-citation">Citation</a>
</p>

---

## 🆕 What's New in v2.1.1

###  Adaptive Threshold Optimizer (ATO)

**Stop using arbitrary thresholds!** The new Adaptive Threshold Optimizer determines data-driven significance cutoffs for differential expression analysis.

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

df = pd.read_csv('deseq2_results.csv')
result = optimize_thresholds(df, goal='discovery')

print(f"Optimal logFC: {result.logfc_threshold:.2f}")
print(f"Significant genes: {result.n_significant}")
print(f"\n{result.methods_text}")  # Publication-ready!
```

**Key Features:**
- Multiple p-value adjustment methods (BH, BY, Storey q-value, Holm, Bonferroni)
- Five logFC optimization methods (MAD, mixture model, power-based, percentile, consensus)
- π₀ estimation for true null proportion
- Three analysis goals: discovery, balanced, validation
- Auto-generated publication methods text
- Interactive dashboard integration

---

##  What is RAPTOR?

**RAPTOR** is a comprehensive framework for benchmarking and optimizing RNA-seq differential expression analysis pipelines. Instead of guessing which pipeline works best for your data, RAPTOR provides **evidence-based, ML-powered recommendations** through systematic comparison of 8 popular pipelines.

### Why RAPTOR?

| Challenge | RAPTOR Solution |
|-----------|-----------------|
| Which pipeline should I use? | ✅ **ML recommendations** with 87% accuracy |
| What thresholds should I use? | ✅ **Adaptive Threshold Optimizer** (NEW!) |
| Is my data quality good enough? | ✅ **Quality assessment** with batch effect detection |
| How do I know results are reliable? | ✅ **Ensemble analysis** combining multiple pipelines |
| What resources do I need? | ✅ **Resource monitoring** with predictions |
| How do I present results? | ✅ **Automated reports** publication-ready |

---

##  Features

<table>
<tr>
<td width="50%">

###  Adaptive Threshold Optimizer (NEW!)
- Data-driven logFC and p-value thresholds
- Multiple statistical methods
- Publication-ready methods text
- Interactive dashboard page

###  ML-Based Recommendations
- 87% prediction accuracy
- Confidence scoring (0-100%)
- Learns from 10,000+ analyses
- Explains its reasoning

###  Quality Assessment
- 6-component quality scoring
- Batch effect detection
- Outlier identification
- Actionable recommendations

</td>
<td width="50%">

###  Ensemble Analysis
- 5 combination methods
- 33% fewer false positives
- High-confidence gene lists
- Consensus validation

###  Interactive Dashboard
- Web-based interface (no coding!)
- Real-time visualizations
- Drag-and-drop data upload
- One-click reports

###  Resource Monitoring
- Real-time CPU/memory tracking
- <1% performance overhead
- Resource predictions
- Cost estimation for cloud

</td>
</tr>
</table>

---

##  Quick Start

### Option 1: Interactive Dashboard (Recommended)

```bash
# Install
pip install raptor-rnaseq

# Launch dashboard
raptor dashboard

# Opens at http://localhost:8501
# Upload data → Get ML recommendation → Use 🎯 Threshold Optimizer → Done!
```

### Option 2: Command Line

```bash
# Profile your data and get ML recommendation
raptor profile --counts counts.csv --metadata metadata.csv --use-ml

# Run recommended pipeline
raptor run --pipeline 3 --data fastq/ --output results/

# Optimize thresholds (NEW!)
raptor optimize-thresholds --input results.csv --goal balanced

# Generate report
raptor report --results results/ --output report.html
```

### Option 3: Python API

```python
from raptor import RNAseqDataProfiler, MLPipelineRecommender
from raptor.threshold_optimizer import optimize_thresholds

# Profile your data
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.run_full_profile()

# Get ML recommendation
recommender = MLPipelineRecommender()
recommendation = recommender.recommend(profile)

print(f"Recommended: Pipeline {recommendation['pipeline_id']}")
print(f"Confidence: {recommendation['confidence']:.1%}")

# After running pipeline, optimize thresholds (NEW!)
de_results = pd.read_csv('de_results.csv')
result = optimize_thresholds(de_results, goal='balanced')
print(f"Optimal |logFC|: {result.logfc_threshold:.2f}")
print(result.methods_text)
```

---

##  Installation

### Requirements

- **Python**: 3.8 or higher
- **R**: 4.0 or higher (for DE analysis)
- **RAM**: 8GB minimum (16GB recommended)
- **Disk**: 10GB free space

### Install from PyPI (Recommended)

```bash
pip install raptor-rnaseq
```

With optional dependencies:

```bash
# With dashboard support
pip install raptor-rnaseq[dashboard]

# With all features
pip install raptor-rnaseq[all]
```

### Install from GitHub

```bash
# Clone repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Install Python dependencies
pip install -r requirements.txt

# Verify installation
python install.py
```

### Conda Environment

```bash
conda env create -f environment.yml
conda activate raptor
```

---

##  Pipelines

RAPTOR benchmarks **8 RNA-seq analysis pipelines**:

| ID | Pipeline | Aligner | Quantifier | DE Tool | Speed | ML Rank |
|----|----------|---------|------------|---------|-------|---------|
| 1 | STAR-RSEM-DESeq2 | STAR | RSEM | DESeq2 | ⭐⭐ | #2 |
| 2 | HISAT2-StringTie-Ballgown | HISAT2 | StringTie | Ballgown | ⭐⭐⭐ | #5 |
| **3** | **Salmon-edgeR** ⭐ | Salmon | Salmon | edgeR | ⭐⭐⭐⭐⭐ | **#1** |
| 4 | Kallisto-Sleuth | Kallisto | Kallisto | Sleuth | ⭐⭐⭐⭐⭐ | #3 |
| 5 | STAR-HTSeq-limma | STAR | HTSeq | limma-voom | ⭐⭐ | #4 |
| 6 | STAR-featureCounts-NOISeq | STAR | featureCounts | NOISeq | ⭐⭐ | #6 |
| 7 | Bowtie2-RSEM-EBSeq | Bowtie2 | RSEM | EBSeq | ⭐⭐ | #7 |
| 8 | HISAT2-Cufflinks-Cuffdiff | HISAT2 | Cufflinks | Cuffdiff | ⭐ | #8 |

⭐ **Pipeline 3 (Salmon-edgeR)** is the ML model's most frequently recommended pipeline due to its optimal speed/accuracy balance.

---

##  Repository Structure

```
RAPTOR/
├── raptor/                 # Core Python package
│   ├── profiler.py         # Data profiling
│   ├── recommender.py      # Rule-based recommendations
│   ├── ml_recommender.py   # ML recommendations
│   ├── threshold_optimizer/ # 🆕 Adaptive Threshold Optimizer (v2.1.1)
│   │   ├── __init__.py
│   │   ├── ato.py          # Core ATO class
│   │   └── visualization.py # ATO visualizations
│   ├── data_quality_assessment.py
│   ├── ensemble_analysis.py
│   ├── resource_monitoring.py
│   └── ...
├── dashboard/              # Interactive web dashboard
├── pipelines/              # Pipeline configurations (8 pipelines)
├── scripts/                # Workflow scripts (00-10)
├── examples/               # Example scripts & demos
├── tests/                  # Test suite
├── docs/                   # Documentation
├── config/                 # Configuration templates
├── install.py              # Master installer
├── launch_dashboard.py     # Dashboard launcher
├── requirements.txt        # Python dependencies
└── setup.py                # Package setup
```

---

##  Documentation

### Getting Started
| Document | Description |
|----------|-------------|
| [INSTALLATION.md](docs/INSTALLATION.md) | Detailed installation guide |
| [QUICK_START.md](docs/QUICK_START.md) | 5-minute quick start |
| [DASHBOARD.md](docs/DASHBOARD.md) | Interactive dashboard guide |

### Core Features
| Document | Description |
|----------|-------------|
| [THRESHOLD_OPTIMIZER.md](docs/THRESHOLD_OPTIMIZER.md) | 🆕 Adaptive threshold optimization |
| [PROFILE_RECOMMEND.md](docs/PROFILE_RECOMMEND.md) | Data profiling & recommendations |
| [QUALITY_ASSESSMENT.md](docs/QUALITY_ASSESSMENT.md) | Quality scoring & batch effects |
| [BENCHMARKING.md](docs/BENCHMARKING.md) | Pipeline benchmarking |

### Advanced Features
| Document | Description |
|----------|-------------|
| [ENSEMBLE.md](docs/ENSEMBLE.md) | Multi-pipeline ensemble analysis |
| [RESOURCE_MONITORING.md](docs/RESOURCE_MONITORING.md) | Resource tracking |
| [CLOUD_DEPLOYMENT.md](docs/CLOUD_DEPLOYMENT.md) | AWS/GCP/Azure deployment |

### Reference
| Document | Description |
|----------|-------------|
| [PIPELINES.md](docs/PIPELINES.md) | Pipeline details & selection guide |
| [API.md](docs/API.md) | Python API reference |
| [FAQ.md](docs/FAQ.md) | Frequently asked questions |
| [CHANGELOG.md](docs/CHANGELOG.md) | Version history |

---

##  Usage Examples

### Example 1: Quick Threshold Optimization (NEW!)

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load DE results
df = pd.read_csv('deseq2_results.csv')

# Optimize thresholds
result = optimize_thresholds(df, goal='balanced')

print(f"Optimal |logFC|: {result.logfc_threshold:.3f}")
print(f"Significant genes: {result.n_significant}")
print(f"π₀ estimate: {result.pi0:.3f}")

# Get publication methods text
print(result.methods_text)

# Save results
result.results_df.to_csv('optimized_results.csv')
```

### Example 2: Full Workflow

```python
from raptor import RNAseqDataProfiler, MLPipelineRecommender
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# 1. Profile data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile(quality_check=True)
print(f"Quality Score: {profile['quality_score']}/100")

# 2. Get ML recommendation
recommender = MLPipelineRecommender()
recommendations = recommender.recommend(profile, n=3)
print(f"Recommended: {recommendations[0]['pipeline_name']}")

# 3. [Run recommended pipeline - produces DE results]
# raptor run --pipeline 3 ...

# 4. Optimize thresholds (NEW in v2.1.1)
de_results = pd.read_csv('deseq2_results.csv')
result = optimize_thresholds(
    de_results,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)

print(f"\n🎯 Optimized Thresholds:")
print(f"   LogFC: |{result.logfc_threshold:.3f}|")
print(f"   Significant: {result.n_significant} genes")

# 5. Save results with methods text
result.results_df.to_csv('final_results.csv')
with open('methods.txt', 'w') as f:
    f.write(result.methods_text)
```

### Example 3: Ensemble Analysis with ATO

```python
from raptor.ensemble_analysis import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

# Combine results from multiple pipelines
analyzer = EnsembleAnalyzer()
consensus = analyzer.combine_results(
    results_dict={'deseq2': df1, 'edger': df2, 'limma': df3},
    method='weighted_vote',
    min_agreement=2
)

# Use ATO for uniform thresholds across ensemble
result = optimize_thresholds(consensus['combined'], goal='balanced')
print(f"Consensus DE genes: {result.n_significant}")
```

---

##  Performance

### ML Recommendation Accuracy

| Metric | Value |
|--------|-------|
| Overall Accuracy | 87% |
| Top-3 Accuracy | 96% |
| Prediction Time | <0.1s |
| Training Data | 10,000+ analyses |

### Threshold Optimizer Benefits

| Metric | Traditional | With ATO |
|--------|-------------|----------|
| Threshold justification | Arbitrary | Data-driven |
| Methods text | Manual | Auto-generated |
| False positives | Higher | Optimized |
| Reproducibility | Variable | Standardized |

---

##  Contributing

We welcome contributions! RAPTOR is open-source and aims to make free science accessible to everyone.

```bash
# Fork and clone
git clone https://github.com/YOUR_USERNAME/RAPTOR.git

# Create feature branch
git checkout -b feature/amazing-feature

# Make changes and test
pytest tests/

# Submit pull request
```

See [CONTRIBUTING.md](docs/CONTRIBUTING.md) for guidelines.

---

##  Citation

If you use RAPTOR in your research, please cite:

```bibtex
@software{bolouki2025raptor,
  author       = {Bolouki, Ayeh},
  title        = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year         = {2025},
  version      = {2.1.1},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17607161},
  url          = {https://github.com/AyehBlk/RAPTOR}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17607161.svg)](https://doi.org/10.5281/zenodo.17607161)

---

##  License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

```
MIT License
Copyright (c) 2025 Ayeh Bolouki
```

---

##  Contact

**Ayeh Bolouki**

- 🏛️ GIGA, University of Liège, Belgium
- 📧 Email: ayehbolouki1988@gmail.com
- 🐙 GitHub: [@AyehBlk](https://github.com/AyehBlk)
- 🔬 Research: Computational Biology, Bioinformatics, Multi-omics Analysis

---

##  Acknowledgments

- The Bioconductor community for the R package ecosystem
- All users who provided feedback

---

<p align="center">
  <strong>⭐ Star this repository if you find RAPTOR useful!</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/github/stars/AyehBlk/RAPTOR?style=social" alt="GitHub Stars"/>
</p>

<p align="center">
  <em>RAPTOR v2.1.1 - Making pipeline selection evidence-based, not guesswork 🦖</em>
</p>
