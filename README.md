<p align="center">
  <img src="https://img.shields.io/badge/🦖_RAPTOR-v2.2.2-brightgreen?style=for-the-badge" alt="RAPTOR v2.2.2"/>
</p>

<h1 align="center">RAPTOR</h1>
<h3 align="center">RNA-seq Analysis Pipeline Testing and Optimization Resource</h3>

<p align="center">
  <strong>Making free science for everybody around the world</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/raptor-rnaseq/"><img src="https://img.shields.io/pypi/v/raptor-rnaseq.svg?style=flat&logo=pypi&logoColor=white" alt="PyPI version"/></a>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/Python-3.8--3.12-blue.svg" alt="Python 3.8-3.12"/></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="MIT License"/></a>
  <a href="https://doi.org/10.5281/zenodo.17607161"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17607161.svg" alt="DOI"/></a>
  <a href="https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.2"><img src="https://img.shields.io/badge/Release-v2.2.2-orange.svg" alt="Release v2.2.2"/></a>
</p>

<p align="center">
  <a href="#-quick-start">Quick Start</a> •
  <a href="#modules">Modules</a> •
  <a href="#installation">Installation</a> •
  <a href="#architecture">Architecture</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#pipelines">Pipelines</a> •
  <a href="#citation">Citation</a>
</p>

---

## What is RAPTOR?

**RAPTOR** is an open-source Python framework that covers the full RNA-seq analysis workflow — from finding and downloading datasets, through quality assessment and differential expression, to ensemble analysis and biomarker discovery.

It's built for researchers who want to spend time on biology, not on figuring out which pipeline to use, what thresholds to set, or how to combine results from different tools. Everything runs from a visual dashboard or the command line.

**What you get:**

- Search and download RNA-seq datasets from GEO, SRA, TCGA, and ArrayExpress directly from a visual dashboard — no coding required. GEO and SRA are fully functional; TCGA and ArrayExpress are under active development.
- Upload your own count matrices and work interactively with sample metadata — add columns, rename groups, exclude samples, assign batches. Then pool multiple studies into one dataset with automatic gene ID harmonization and batch correction. RAPTOR checks whether your datasets are actually poolable and flags problems like library size differences, low gene overlap, or batch effects before you proceed.
- Profile your data across 32 features (BCV, dispersion, sparsity, sample balance, and more) and get pipeline recommendations. RAPTOR offers two recommendation modes: an ML-based approach using a Random Forest classifier, and a rule-based approach for simpler guidance.
- Run differential expression with DESeq2, edgeR, and limma, then combine results into a consensus through ensemble analysis. Five combination methods (Fisher's, Brown's, RRA, Voting, Weighted) reduce false positives by about 33% compared to any single method.
- Optimize your significance thresholds using data-driven methods instead of arbitrary cutoffs. Four approaches available: Ground Truth (calibrate against known DE genes if you have them), FDR Control, Stability (consistent across subsamples), and Reproducibility (replicates across datasets). Works with or without prior data, but prior data like ground truth gives better calibration.
- Validate your findings across independent cohorts for robust biomarker candidates. Module 10 (Biomarker Discovery, under development — target mid-April 2026) will support multiple discovery scenarios: single-cohort candidate identification, cross-cohort validation by pooling your data with public studies, consistency ranking across studies, effect size filtering, and pathway enrichment analysis. The goal is to surface biomarker candidates that are not just statistically significant in one experiment, but reproducible across independent datasets from different labs and patient populations.

---

## Modules

Everything below is accessible through a visual dashboard — no coding required. 10 pages covering the full workflow from data acquisition to visualization. Launch with `python -m streamlit run raptor/dashboard/app.py`.

### Data Acquisition (Module 6b)

Your analysis starts with data. RAPTOR connects to public repositories so you can search, preview, and download RNA-seq datasets without leaving the dashboard.

| Repository | Status | What you get |
|------------|--------|-------------|
| **GEO** (NCBI) | Ready | Search 200,000+ datasets, download processed count matrices |
| **SRA** (NCBI) | Ready | Explore run tables, auto-detect linked GEO studies, generate FASTQ download commands |
| **TCGA** (NCI) | In development | 33 cancer types via GDC API, multi-omic (gene expression, miRNA, methylation, CNV, RPPA) |
| **ArrayExpress** (EBI) | In development | European studies via BioStudies API |

You can also upload your own count matrix — from your own experiment or from collaborators at your institute. RAPTOR treats uploaded data the same way as public data.

**Interactive metadata editing:** Before pooling, you can work directly with sample metadata in the dashboard. Add new columns (condition, timepoint, tissue), rename groups, assign batch labels, exclude samples you don't want to include — all from a visual table editor with color-coded datasets.

**Dataset pooling:** Select 2 or more datasets and RAPTOR merges them into one. It harmonizes gene IDs across studies (Ensembl, Symbol, Entrez via MyGene.info), applies batch correction (ComBat, quantile normalization, or median ratio), and runs statistical checks to tell you whether the datasets are actually reliable enough to pool. If they're not, it flags the problems — library size differences, low gene overlap, batch clustering in PCA — and suggests how to fix them.

This is the foundation for cross-cohort biomarker validation — a gene that is consistently differentially expressed across your data and independent public cohorts is a far more credible biomarker candidate than one found in a single study.

### Quality Assessment (Module 2)

Six outlier detection methods with consensus-based reporting. Detects batch effects, library size issues, and sample quality problems. Gives you actionable recommendations — not just flags.

### Data Profiling (Module 3)

Extracts 32 features from your dataset: BCV, sample balance, dispersion patterns, sparsity, and more. These features feed into the ML recommender to help you choose the right pipeline.

### Pipeline Recommender (Module 4)

A Random Forest classifier trained on dataset features tells you which quantification pipeline to use and how confident it is. Also provides rule-based recommendations as a simpler alternative.

### Production Pipelines (Module 5)

Six quantification pipelines, all supporting paired-end and single-end reads:

| Pipeline | Speed | Memory | Best for |
|----------|-------|--------|----------|
| Salmon | Fast | 8 GB | Standard DE analysis (recommended) |
| Kallisto | Fastest | 4 GB | Speed priority |
| STAR + featureCounts | Moderate | 32 GB | Gene-level publication |
| STAR + RSEM | Slow | 32 GB | Isoform analysis |
| STAR + Salmon | Moderate | 32 GB | BAM + bootstraps |
| HISAT2 + featureCounts | Moderate | 16 GB | Lower memory systems |

### Differential Expression (Module 6 + Module 7)

Python wrappers for R-based DE analysis (DESeq2, edgeR, limma). Module 7 imports and standardizes results from any DE tool into a common format for downstream analysis.

### Parameter Optimization (Module 8)

Four methods for finding optimal thresholds instead of using arbitrary cutoffs:

- **Ground Truth** — when you have known DE genes to calibrate against
- **FDR Control** — data-driven FDR threshold selection
- **Stability** — thresholds that give consistent results across subsamples
- **Reproducibility** — thresholds that replicate across independent datasets

### Ensemble Analysis (Module 9)

Combine results from DESeq2, edgeR, and limma into a single consensus. Five statistical methods: Fisher's, Brown's (accounts for method correlation), Robust Rank Aggregation, Voting, and Weighted Ensemble. Reduces false positives by about 33% compared to any single method.

### Biomarker Discovery (Module 10)

Under development — target mid-April 2026. This module consumes pooled datasets from Module 6b and applies multiple discovery strategies depending on your research scenario:

- **Single-cohort discovery** — identify candidate biomarkers from your own study using differential expression, effect size filtering, and statistical confidence
- **Cross-cohort validation** — pool your data with 2-3 public studies of the same condition and find genes that are consistently differentially expressed across all cohorts
- **Consistency ranking** — rank candidates by how reproducible they are across independent datasets, not just by p-value in a single experiment
- **Effect size filtering** — prioritize candidates with biologically meaningful fold changes, not just statistically significant ones
- **Pathway enrichment** — connect individual gene candidates to biological pathways and functional categories

The key idea: a biomarker that survives batch correction and remains significant across 3+ independent cohorts from different labs and patient populations is far more credible than one found in a single study. RAPTOR's Data Acquisition and Pooling modules make this cross-cohort workflow accessible without custom integration code.

---

## Current Version

**GitHub: v2.2.2** — includes Data Acquisition module, all fixes, and 621 tests passing.

**PyPI: v2.2.1** — the PyPI release will be updated after all modules (including TCGA, ArrayExpress, and Biomarker Discovery) are complete and tested.

To get the latest version, install from GitHub:

```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -e .
pip install streamlit GEOparse biopython mygene
```

We are actively testing the Data Acquisition module with researchers. If you want to help, see the [Beta Testing Guide](BETA_TESTING_GUIDE.md) and report issues at [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues).

---

## Quick Start

### Option 1: Interactive Dashboard (Recommended)

**Install with dashboard support:**
```bash
pip install raptor-rnaseq[dashboard]
```

**Launch the dashboard** — choose the method that matches your setup:

| Scenario | Command |
|----------|---------|
| **pip install** (recommended) | `python -m raptor.launch_dashboard` |
| **Cloned from GitHub** | `python launch_dashboard.py` |
| **Direct streamlit** | `python -m streamlit run raptor/dashboard/app.py` |
| **Inside a virtual environment** | Activate venv first, then any command above |

The dashboard opens at **http://localhost:8501**. Search & download data → Profile → Get recommendation → Run ensemble → Done!

**Detailed instructions by setup:**

<details>
<summary><strong>A. Installed via pip (most users)</strong></summary>

```bash
# Install
pip install raptor-rnaseq[dashboard]

# Launch
python -m raptor.launch_dashboard
```
</details>

<details>
<summary><strong>B. Cloned from GitHub (developers)</strong></summary>

```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -e .[dashboard]

# Launch from repo root
python launch_dashboard.py

# Or directly
python -m streamlit run raptor/dashboard/app.py
```
</details>

<details>
<summary><strong>C. Using a virtual environment</strong></summary>

```bash
# Create and activate venv
python -m venv .venv

# Windows
.venv\Scripts\Activate

# Linux/macOS
source .venv/bin/activate

# Install and launch
pip install raptor-rnaseq[dashboard]
python -m raptor.launch_dashboard
```
</details>

<details>
<summary><strong>D. Conda environment</strong></summary>

```bash
conda env create -f environment.yml
conda activate raptor
pip install streamlit altair
python -m streamlit run raptor/dashboard/app.py
```
</details>

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
    ensemble_brown,
    GEOConnector,        # Data Acquisition
    GeneIDMapper,        # Gene ID conversion
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

## Installation

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

# Or with all features (dashboard + acquisition + dev tools)
pip install -e .[all]

# Or install acquisition connectors separately
pip install GEOparse biopython mygene

# Verify installation
raptor --version
pytest tests/
```

---

## Architecture

RAPTOR is organized into **10 modules** spanning 5 analysis stages:

```
┌─────────────────────────────────────────────────────────────┐
│                    RAPTOR v2.2.2                            │
│         RNA-seq Analysis Pipeline Framework                 │
└─────────────────────────────────────────────────────────────┘

Stage 1: Data Acquisition & Preparation
├── Module 6b: Data Acquisition (NEW in v2.2.2)
│    ├── GEO Connector (ready)
│    ├── SRA Connector (ready)
│    ├── TCGA Connector (in progress)
│    ├── ArrayExpress Connector (in progress)
│    ├── Gene ID Mapper (Ensembl/Symbol/Entrez)
│    ├── Pooling Engine (batch correction)
│    └── Upload your own data
├── Module 1: Quick Quantification (Salmon/Kallisto)
├── Module 2: Quality Assessment (6 outlier methods)
└── Module 3: Data Profiling (32 features)

Stage 2: Pipeline Selection
├── Module 4: ML Recommender (Random Forest)
└── Module 5: Production Pipelines (6 methods)
         ├── Salmon (recommended)
         ├── Kallisto (fastest)
         ├── STAR + featureCounts
         ├── STAR + RSEM
         ├── STAR + Salmon
         └── HISAT2 + featureCounts

Stage 3: Differential Expression
├── Module 6: DE Analysis (R: DESeq2, edgeR, limma)
└── Module 7: DE Import (standardize any format)

Stage 4: Advanced Analysis
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

Stage 5: Biomarker Discovery (under development — target mid-April 2026)
└── Module 10: Biomarker Discovery
     ├── Cross-cohort validation
     ├── Candidate ranking by consistency
     └── Pathway enrichment
```

---

## Pipelines

RAPTOR supports **6 production RNA-seq quantification pipelines**:

| Pipeline | Memory | Time | Produces | Best For | Recommended |
|----------|--------|------|----------|----------|-------------|
| **salmon** | 8 GB | 10-20 min | genes + isoforms + bootstraps | Standard DE analysis | Recommended |
| **kallisto** | 4 GB | 5-10 min | genes + isoforms + bootstraps | Speed priority | Yes |
| **star_featurecounts** | 32 GB | 40-70 min | BAM + genes | Gene-level publication | Yes |
| **star_rsem** | 32 GB | 60-120 min | BAM + genes + isoforms | Isoform analysis | Yes |
| **star_salmon** | 32 GB | 50-90 min | BAM + genes + isoforms + bootstraps | BAM + bootstraps | Yes |
| **hisat2_featurecounts** | 16 GB | 30-60 min | BAM + genes | Low memory systems | Yes |

**Salmon is recommended** for most use cases due to optimal speed/accuracy balance and bootstrap support.

### Pipeline Features

**All pipelines support:**
- Paired-end and single-end reads
- Automatic parameter optimization
- QC report generation
- Multi-threading
- Sample sheet-based workflows

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

## Repository Structure

```
RAPTOR/
├── raptor/                         # Core Python package
│   ├── __init__.py                 # Package initialization (v2.2.2)
│   ├── cli.py                      # Command-line interface (55 commands)
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
│   ├── external_modules/           # Module 6: R integration + Data Acquisition
│   │   ├── module6_de_analysis/
│   │   │   └── r_scripts/          # DESeq2, edgeR, limma
│   │   └── acquisition/            # Module 6b: Data Acquisition (NEW in v2.2.2)
│   │       ├── geo.py              # GEO connector (ready)
│   │       ├── sra.py              # SRA connector (ready)
│   │       ├── tcga.py             # TCGA connector (in progress)
│   │       ├── arrayexpress.py     # ArrayExpress connector (in progress)
│   │       ├── gene_mapping.py     # Gene ID conversion
│   │       ├── pooling.py          # Dataset pooling + batch correction
│   │       ├── datasets.py         # AcquiredDataset, PooledDataset
│   │       ├── cache.py            # Disk caching system
│   │       ├── catalog.py          # Dataset registry
│   │       └── base.py             # BaseConnector, SearchResult
│   │
│   ├── dashboard/                  # Interactive Streamlit app
│   │   ├── app.py
│   │   ├── pages/                  # 10 dashboard pages
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
├── tests/                          # Test suite (621 tests passing)
│   ├── test_acquisition.py                        (105 tests, NEW)
│   ├── test_profiler.py
│   ├── test_quality_assessment.py
│   ├── test_parameter_optimization.py
│   ├── test_ensemble.py
│   ├── test_cli_comprehensive.py                  (55 commands)
│   └── ...
│
├── templates/                      # Sample sheets
│   ├── sample_sheet_paired.csv
│   └── sample_sheet_single.csv
│
├── .github/                        # GitHub templates
│   └── ISSUE_TEMPLATE/
│
├── .streamlit/                     # Dashboard theme configuration
│   └── config.toml
│
├── setup.py                        # Package setup
├── MANIFEST.in                     # Package file inclusions
├── requirements.txt                # Python dependencies
├── environment.yml                 # Conda environment (core)
├── environment-full.yml            # Conda environment (complete)
├── check_raptor.py                 # Diagnostic suite (16 checks)
├── BETA_TESTING_GUIDE.md           # Beta testing scenarios
├── CITATION.cff                    # Citation metadata
├── CHANGELOG.md                    # Version history
├── CONTRIBUTING.md                 # Contribution guidelines
├── .gitignore                      # Git ignore rules
└── LICENSE                         # MIT License
```

---

## Documentation

### Getting Started
| Document | Description |
|----------|-------------|
| [Quick Start](#-quick-start) | 5-minute quick start guide |
| [Installation](#-installation) | Detailed installation instructions |
| [CONDA_ENVIRONMENTS.md](docs/CONDA_ENVIRONMENTS.md) | Conda setup (core vs full) |

### Core Features
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

## Usage Examples

### Example 1: Complete Workflow

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
# Complete RAPTOR workflow using CLI

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

### Example 4: Data Acquisition (Python API)

```python
from raptor.external_modules.acquisition import (
    GEOConnector, SRAConnector, GeneIDMapper, PoolingEngine
)

# Search GEO for datasets
geo = GEOConnector()
results = geo.search("pancreatic cancer RNA-seq human", max_results=10)
for r in results:
    print(f"{r.accession}: {r.title} ({r.n_samples} samples)")

# Download a dataset
dataset = geo.download("GSE12345")
print(f"Downloaded: {dataset.n_genes} genes x {dataset.n_samples} samples")

# Search SRA and find linked GEO studies
sra = SRAConnector()
sra_results = sra.search("PS19 tauopathy mouse RNA-seq")
run_table = sra.get_run_table(sra_results[0].accession)
linked_gse = sra.find_linked_gse(run_table)

# Convert gene IDs between systems
mapper = GeneIDMapper("Homo sapiens")
id_type = mapper.detect_id_type(dataset.counts_df.index.tolist())
converted = mapper.convert_index(dataset.counts_df, from_type=id_type, to_type="symbol")

# Pool multiple datasets with batch correction
engine = PoolingEngine(target_gene_id="symbol", species="Homo sapiens")
pooled = engine.merge([dataset1, dataset2], method="inner")
print(f"Pooled: {pooled.n_genes} genes, {pooled.n_samples} samples, {pooled.n_studies} studies")
```

---

## Performance

### Module Performance

| Module | Time | Memory | Key Output |
|--------|------|--------|------------|
| **Module 6b: Data Acquisition** | Varies | 4 GB | Downloaded datasets, pooled counts |
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

## Contributing

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

- **Test the Data Acquisition module** — see [Beta Testing Guide](BETA_TESTING_GUIDE.md)
- **Report bugs** via [Issues](https://github.com/AyehBlk/RAPTOR/issues)
- **Request features** via [Issues](https://github.com/AyehBlk/RAPTOR/issues)
- **Improve documentation**
- **Submit pull requests**
- **Share use cases** and feedback
- **Star the repository**

### Beta Testing (Active)

We are actively testing the Data Acquisition module (GEO and SRA search, dataset pooling, quality check). Your feedback on real-world datasets directly improves the tool.

- [Beta Testing Guide](BETA_TESTING_GUIDE.md) — 9 scenarios to try
- [Report an Issue](https://github.com/AyehBlk/RAPTOR/issues) — bug reports and feature requests
- [Project Board](https://github.com/users/AyehBlk/projects/5) — track progress on reported issues

---

## Citation

If you use RAPTOR in your research, please cite:

```bibtex
@software{bolouki2026raptor,
  author       = {Bolouki, Ayeh},
  title        = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year         = {2026},
  version      = {2.2.2},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17607161},
  url          = {https://github.com/AyehBlk/RAPTOR}
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17607161.svg)](https://doi.org/10.5281/zenodo.17607161)

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

```
MIT License
Copyright (c) 2026 Ayeh Bolouki
```

---

## Contact

**Ayeh Bolouki**

- **GIGA, University of Liege, Belgium**
- **Email:** ayehbolouki1988@gmail.com
- **GitHub:** [@AyehBlk](https://github.com/AyehBlk)
- **Research:** Computational Biology, Bioinformatics, Multi-omics Analysis

### Support

- **Documentation:** [docs/](docs/)
- **Issues:** [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- **Discussions:** [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)
- **Email:** ayehbolouki1988@gmail.com

---

## Known Issues & Troubleshooting

### v2.2.2 Known Issues

| Issue | Status | Workaround |
|-------|--------|------------|
| `raptor dashboard` CLI command fails with `ModuleNotFoundError` | Fix planned | Use `python -m streamlit run raptor/dashboard/app.py` |
| TCGA connector incomplete | In development | Use GEO for processed counts |
| ArrayExpress connector incomplete | In development | Use GEO as alternative |
| PyPI version is still 2.2.1 | By design | Install from GitHub for v2.2.2 (see below) |
| ML Recommender requires pre-trained model file | By design | Use rule-based recommender (`--method rule`) |

### Install v2.2.2 from GitHub

PyPI will be updated after all modules are complete and tested. For now:

```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS / Linux
source .venv/bin/activate

pip install -e .
pip install streamlit GEOparse biopython mygene
python -m streamlit run raptor/dashboard/app.py
```

### Common Errors

<details>
<summary><strong>"streamlit: not recognized" or "command not found"</strong></summary>

Streamlit is not installed or not on your PATH.

```bash
# Fix: install with dashboard extras
pip install raptor-rnaseq[dashboard]

# Or install streamlit directly
pip install streamlit

# If using a virtual environment, make sure it's activated first
# Windows: .venv\Scripts\Activate
# Linux/macOS: source .venv/bin/activate
```
</details>

<details>
<summary><strong>"ModuleNotFoundError: No module named 'raptor'"</strong></summary>

This can happen when running `raptor dashboard` from the CLI.

```bash
# Workaround: run streamlit directly
python -m streamlit run raptor/dashboard/app.py
```
</details>

<details>
<summary><strong>"raptor --version" shows wrong version</strong></summary>

If you installed from PyPI, it will show 2.2.1 (the latest PyPI release). For v2.2.2, install from GitHub.

```bash
# Check the version
python -c "import raptor; print(raptor.__version__)"
# Should print: 2.2.2 (GitHub) or 2.2.1 (PyPI)
```
</details>

<details>
<summary><strong>Dashboard pages show errors with no data loaded</strong></summary>

The dashboard requires data to be uploaded or loaded into session state before most pages will work. Start by:

1. Go to **Data Acquisition** to search and download datasets from GEO or SRA
2. Or go to the **Quality Assessment** page and upload a count matrix CSV and metadata CSV
3. Or go to **Import DE** and upload DESeq2/edgeR/limma result files
4. Pages will populate once data is available in the session
</details>

<details>
<summary><strong>Import errors on Windows (path issues)</strong></summary>

If you see path-related errors on Windows:

```bash
# Make sure you're not inside a folder named 'raptor' — this shadows the package
cd C:\Users\YourName
python -c "import raptor; print(raptor.__version__)"

# If using the source code, install in editable mode
cd C:\path\to\RAPTOR
pip install -e .
```
</details>

### Reporting Bugs

Found a bug? We have structured templates to make reporting easy:

1. Check [existing issues](https://github.com/AyehBlk/RAPTOR/issues) first
2. Create a [Bug Report](https://github.com/AyehBlk/RAPTOR/issues/new?template=bug_report.md) or [Feature Request](https://github.com/AyehBlk/RAPTOR/issues/new?template=feature_request.md)
3. Include your RAPTOR version (`python -c "import raptor; print(raptor.__version__)"`) and the dataset accession if applicable

---

## Acknowledgments

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
  <strong>Star this repository if you find RAPTOR useful!</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/github/stars/AyehBlk/RAPTOR?style=social" alt="GitHub Stars"/>
  <img src="https://img.shields.io/github/forks/AyehBlk/RAPTOR?style=social" alt="GitHub Forks"/>
</p>

<p align="center">
  <em>RAPTOR v2.2.2 — Making pipeline selection evidence-based, not guesswork</em>
</p>

<p align="center">
  <em>Making free science for everybody around the world</em>
</p>