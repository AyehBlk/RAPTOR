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

- Search and download RNA-seq datasets from GEO, SRA, and TCGA directly from a visual dashboard — no coding required. GEO and SRA are production-ready; TCGA is functional with multi-omic support (tested with 5 cancer cohorts); ArrayExpress is under active development.
- Upload your own count matrices and work interactively with sample metadata — add columns, rename groups, exclude samples, assign batches. Then pool multiple studies into one dataset with automatic gene ID harmonization and batch correction. RAPTOR checks whether your datasets are actually poolable and flags problems like library size differences, low gene overlap, or batch effects before you proceed.
- Profile your data across 32 features (BCV, dispersion, sparsity, sample balance, and more) and get pipeline recommendations. RAPTOR offers two recommendation modes: an ML-based approach using a Random Forest classifier, and a rule-based approach for simpler guidance.
- Run differential expression with DESeq2, edgeR, and limma, then combine results into a consensus through ensemble analysis. Five combination methods (Fisher's, Brown's, RRA, Voting, Weighted) reduce false positives by about 33% compared to any single method.
- Optimize your significance thresholds using data-driven methods instead of arbitrary cutoffs. Four approaches available: Ground Truth (calibrate against known DE genes if you have them), FDR Control, Stability (consistent across subsamples), and Reproducibility (replicates across datasets). Works with or without prior data, but prior data like ground truth gives better calibration.
- Discover gene-panel biomarkers with full nested cross-validation, multi-method feature selection, and clinical-utility metrics. Module 10 (Biomarker Discovery) ships in v2.2.2 with eight feature-selection methods (univariate filter, Elastic Net, LASSO, Boruta, mRMR, RFE, SHAP, WGCNA), four classifiers (Logistic Regression, Random Forest, SVM, XGBoost), and the Ambroise & McLachlan (2002 PNAS) leakage fix that runs feature selection inside each outer CV fold. Panel size is auto-detected by the kneedle algorithm with consensus pinning across folds. The pipeline reports out-of-fold clinical metrics (Youden's threshold, sensitivity, specificity, PPV/NPV at user-supplied prevalence, decision curve analysis), the Nogueira (2018) panel stability index Phi with bootstrap CI, and an optimism diagnostic comparing CV AUC to training-set AUC. A tiered small-cohort advisory and a multi-seed verification button surface the run-to-run variance that is genuinely present at small n. Six biomarker intent modes (diagnostic, prognostic, predictive, monitoring, exploratory, translational) configure the analysis end-to-end; diagnostic and exploratory are fully implemented, the others are scaffolded for v2.3.0.

---

## Modules

Everything below is accessible through a visual dashboard — no coding required. 10 pages covering the full workflow from data acquisition to visualization. Launch with `python -m streamlit run raptor/dashboard/app.py`.

### Data Acquisition (Module 6b)

Your analysis starts with data. RAPTOR connects to public repositories so you can search, preview, and download RNA-seq datasets without leaving the dashboard.

| Repository | Status | What you get |
|------------|--------|-------------|
| **GEO** (NCBI) | Ready | Search 200,000+ datasets, download processed count matrices |
| **SRA** (NCBI) | Ready | Explore run tables, auto-detect linked GEO studies, generate FASTQ download commands |
| **TCGA** (NCI) | Functional | 33 cancer types via GDC API, multi-omic (gene expression, miRNA, methylation, CNV, RPPA, mutations). Tested with 5 cancer cohorts |
| **ArrayExpress** (EBI) | Functional | European studies via BioStudies API |

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

End-to-end gene-panel biomarker discovery with nested cross-validation, panel-size auto-detection, classification evaluation, clinical-utility metrics, and biological annotation. Accessible via both the dashboard (page 08) and the CLI (`raptor biomarker`).

**Pipeline structure.** The pipeline runs feature selection, panel-size optimization, and classifier fitting **inside** each outer cross-validation fold, following Ambroise & McLachlan (2002 PNAS 99:6562). When feature selection sees the labels of validation samples, the resulting CV AUC is optimistically biased; pushing every label-using step inside the training fold restores honest generalization estimates. A separate "final-panel pass" then re-runs the same procedure on all data to produce the deployable panel. The OOF AUC tells you how well the procedure generalizes; the final panel is what you actually use.

**Feature selection.** Eight methods are available, with the default set (univariate filter, Elastic Net, RFE) chosen for fold-safe behavior:

- **Univariate filter** (default) — per-fold Welch's t-test or Mann-Whitney U
- **Elastic Net / LASSO** — penalized logistic regression
- **Boruta** — Random Forest with shadow features (Kursa & Rudnicki 2010)
- **mRMR** — minimum-redundancy maximum-relevance (Ding & Peng 2005)
- **RFE** — recursive feature elimination with Random Forest backbone
- **SHAP** — Shapley-value feature importance (Lundberg & Lee 2017)
- **WGCNA** — co-expression hub genes (Langfelder & Horvath 2008; PyWGCNA)

Method results are aggregated into a consensus ranking weighted 70% by mean rank across methods and 30% by selection frequency. The consensus is then post-processed by significance calibration: genes whose univariate p-value falls below a threshold (default alpha = 0.05) keep full weight; non-significant genes are demoted by a 0.5x multiplicative shrinkage. Forward selection draws its top candidates from the calibrated ordering.

**Panel-size auto-detection.** The panel-size-vs-AUC curve is fit with the kneedle algorithm (Satopaa et al. 2011), which detects the maximum-curvature point on a polynomially smoothed normalized curve. Three failure modes of the legacy first-drop heuristic are handled correctly: non-monotone curves, saturated curves (kneedle falls back to argmax), and late knees past panel size 30. Three strategies are exposed: `kneedle` (default), `argmax`, and `first_drop` (legacy). By default, panel-size detection runs once on the full data and all CV folds are pinned to the same K (`panel_size_strategy='consensus'`). Per-fold detection is available as opt-in.

**Classification.** Four classifiers are evaluated by nested CV with deterministic tiebreak (LR > SVM > RF > XGBoost when AUCs and F1 coincide within numerical tolerance). Held-out predictions accumulate into honest out-of-fold arrays used by all downstream clinical metrics.

**Clinical metrics.** Computed on out-of-fold predictions (not in-sample fits):

- **Youden's optimal threshold** — sensitivity + specificity - 1 (Youden 1950)
- **Sensitivity, specificity** at the optimal threshold
- **PPV/NPV at user-supplied prevalence** — Bayes' theorem-adjusted predictive values
- **Bootstrap 95% CI on AUC** — non-parametric, 2000 resamples
- **Decision curve analysis** — net benefit across threshold probabilities (Vickers & Elkin 2006)
- **Net reclassification improvement** — when a baseline model is supplied (Pencina et al. 2008)

**Panel stability.** Nogueira (2018 JMLR 18:174) Phi is computed across the per-fold panels, the only similarity-based stability measure that satisfies all five of Nogueira's axiomatic properties including correction-for-chance. Reported with bootstrap 95% CI and benchmark labels (Phi >= 0.75 excellent, 0.4-0.75 intermediate, < 0.4 poor).

**Optimism diagnostic.** Reports CV AUC vs training-set AUC and the gap between them. A large gap flags overfitting; the dashboard surfaces it as a colored banner.

**Direction patterns and ratio biomarkers.** Direction patterns separate UP-regulated and DOWN-regulated panel members and support cross-cohort concordance checking. Ratio biomarkers run a generalized Top Scoring Pair search (Geman et al. 2004), evaluating all pairwise gene ratios for self-normalizing biomarkers that cancel batch effects, library-size differences, and platform variation.

**Biological annotation.** MyGene.info gene info, KEGG/Reactome/GO pathway enrichment via gseapy (or manual Fisher's exact test if gseapy is unavailable), Europe PMC literature mining with disease-context filtering, STRING protein-protein interaction network with PPI-enrichment p-value.

**Small-cohort tooling.** At small sample sizes (n < 50, especially n < 20), single-seed CV results carry substantial run-to-run variance. The dashboard surfaces this in two ways: a tiered advisory at the top of the Clinical Metrics tab (strong below n=20, moderate below n=50), and a multi-seed verification button on the Panel tab that re-runs discovery at additional seeds and produces a side-by-side truth table of stable vs variable diagnostics.

**Intent system.** Six study-design modes auto-configure analysis: `diagnostic` and `exploratory` are fully implemented; `prognostic`, `predictive`, `monitoring`, and `translational` are scaffolded for future releases.

---

## Current Version

**v2.2.2** (April 2026) — composite release adding Module 6b (Data Acquisition) and Module 10 (Biomarker Discovery). Available on both PyPI and GitHub.

```bash
pip install raptor-rnaseq[dashboard]
```

Highlights of this release:

- **Module 6b: Data Acquisition.** Search and download RNA-seq datasets from GEO, SRA, and TCGA from a visual dashboard. Pool multiple studies with batch correction (ComBat, quantile, median ratio) and gene-ID harmonization (Ensembl/Symbol/Entrez via MyGene.info).
- **Module 10: Biomarker Discovery.** End-to-end gene-panel discovery with nested CV, kneedle panel-size detection, four classifiers, out-of-fold clinical metrics, panel stability index (Nogueira Phi), optimism diagnostic, biological annotation, and a tiered small-cohort advisory. Six biomarker intent modes; diagnostic and exploratory fully implemented.
- **CLI extension.** `raptor biomarker` and `raptor biomarker-validate` for the discovery and validation flows. Three new flags (`--auto-panel-strategy`, `--panel-sensitivity`, `--panel-size-strategy`) expose kneedle and consensus-pinning options.
- **Polish.** Feature-selection-inside-CV leakage fix (Ambroise & McLachlan 2002 PNAS), Nogueira (2018 JMLR) panel stability, deterministic best-classifier tiebreak, multi-seed verification button on the dashboard.

See [CHANGELOG.md](CHANGELOG.md) for the complete v2.2.2 release notes.

We are actively testing both the Data Acquisition and Biomarker Discovery modules with researchers. If you want to help, see the [Beta Testing Guide](BETA_TESTING_GUIDE.md) and report issues at [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues).

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

# 5. Optimize thresholds
raptor optimize --de-result de_results.csv --method fdr-control --fdr-target 0.05

# 6. Ensemble analysis - combine all methods
raptor ensemble-compare --deseq2 deseq2.csv --edger edger.csv --limma limma.csv

# 7. Biomarker discovery (with diagnostic intent and 5% prevalence)
raptor biomarker \
    --counts counts.csv --metadata metadata.csv \
    --group-column condition \
    --intent diagnostic --prevalence 0.05 \
    --output results/biomarkers/

# 8. Validate the panel on an independent cohort
raptor biomarker-validate \
    --panel results/biomarkers/biomarker_panel.csv \
    --counts validation_counts.csv \
    --metadata validation_metadata.csv
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
from raptor.biomarker_discovery import discover_biomarkers

# 1. Quality check
qc_report = quick_quality_check('counts.csv', 'metadata.csv')
print(f"Outliers: {qc_report.outliers}")

# 2. Profile data (32 features extracted)
profile = profile_data_quick('counts.csv', 'metadata.csv', group_column='condition')
print(f"BCV: {profile.bcv:.3f} ({profile.bcv_category})")

# 3. Get ML recommendation
recommendation = recommend_pipeline(profile_file='profile.json', method='ml')
print(f"Recommended: {recommendation.pipeline_name} (confidence: {recommendation.confidence:.2f})")

# 4. After running DE analysis, optimize thresholds
result = optimize_with_fdr_control(de_result, fdr_target=0.05)
print(f"Optimal thresholds: {result.optimal_threshold}")

# 5. Ensemble analysis - combine DESeq2, edgeR, limma
consensus = ensemble_brown({
    'deseq2': deseq2_result,
    'edger': edger_result,
    'limma': limma_result,
})
print(f"Consensus DE genes: {len(consensus.consensus_genes)}")

# 6. Biomarker discovery with diagnostic intent
biomarker_result = discover_biomarkers(
    counts='counts.csv',
    metadata='metadata.csv',
    group_column='condition',
    intent='diagnostic',
    prevalence=0.05,
    output_dir='results/biomarkers/',
)
print(f"Panel: {biomarker_result.panel}")
best = biomarker_result.classification_results[biomarker_result.best_classifier]
print(f"OOF AUC: {best.auc:.3f}")
if biomarker_result.panel_stability is not None:
    print(biomarker_result.panel_stability.summary())
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
+-------------------------------------------------------------+
|                    RAPTOR v2.2.2                            |
|         RNA-seq Analysis Pipeline Framework                 |
+-------------------------------------------------------------+

Stage 1: Data Acquisition & Preparation
+-- Module 6b: Data Acquisition
|    +-- GEO Connector (production-ready)
|    +-- SRA Connector (production-ready)
|    +-- TCGA Connector (functional, multi-omic: gene expression,
|    |                   miRNA, methylation, CNV, RPPA)
|    +-- ArrayExpress Connector (in development)
|    +-- Gene ID Mapper (Ensembl/Symbol/Entrez via MyGene.info)
|    +-- Pooling Engine (ComBat / quantile / median-ratio batch correction)
|    +-- Data Catalog and Cache Manager
+-- Module 1: Quick Quantification (Salmon/Kallisto)
+-- Module 2: Quality Assessment (6 outlier methods)
+-- Module 3: Data Profiling (32 features)

Stage 2: Pipeline Selection
+-- Module 4: ML Recommender (Random Forest + rule-based fallback)
+-- Module 5: Production Pipelines (6 methods)
     +-- Salmon (recommended)
     +-- Kallisto (fastest)
     +-- STAR + featureCounts
     +-- STAR + RSEM
     +-- STAR + Salmon
     +-- HISAT2 + featureCounts

Stage 3: Differential Expression
+-- Module 6: DE Analysis (R: DESeq2, edgeR, limma)
+-- Module 7: DE Import (standardize any format)

Stage 4: Advanced DE Analysis
+-- Module 8: Parameter Optimization (4 methods)
|    +-- Ground Truth Optimization
|    +-- FDR Control Optimization
|    +-- Stability Optimization
|    +-- Reproducibility Optimization
+-- Module 9: Ensemble Analysis (5 methods)
     +-- Fisher's Method
     +-- Brown's Method (recommended)
     +-- Robust Rank Aggregation
     +-- Voting Consensus
     +-- Weighted Ensemble

Stage 5: Biomarker Discovery
+-- Module 10: Biomarker Discovery
     +-- Feature selection (8 methods, fold-safe defaults)
     |    +-- Univariate filter (default; Welch's t-test, Mann-Whitney U)
     |    +-- Elastic Net / LASSO
     |    +-- Boruta, mRMR, RFE, SHAP, WGCNA (when libraries installed)
     +-- Consensus ranking with significance calibration
     +-- Panel-size auto-detection (kneedle algorithm)
     |    +-- Consensus pinning across CV folds (default)
     |    +-- Per-fold detection (opt-in)
     +-- Classification (4 classifiers, deterministic tiebreak)
     |    +-- Logistic Regression, Random Forest, SVM, XGBoost
     +-- Pipeline-CV with Ambroise & McLachlan (2002) leakage fix
     +-- Out-of-fold clinical metrics
     |    +-- Youden's threshold, sensitivity, specificity
     |    +-- PPV/NPV at user-supplied prevalence
     |    +-- Decision curve analysis (Vickers & Elkin 2006)
     |    +-- Bootstrap 95% CI on AUC
     +-- Panel stability (Nogueira 2018 Phi)
     +-- Optimism diagnostic (CV vs training AUC)
     +-- Direction patterns (UP/DOWN, cross-cohort concordance)
     +-- Ratio biomarkers (Geman et al. 2004 TSP)
     +-- Biological annotation (MyGene, KEGG/Reactome/GO, Europe PMC, STRING)
     +-- Six biomarker intent modes (diagnostic/exploratory implemented;
                                     prognostic/predictive/monitoring/
                                     translational scaffolded for v2.3.0)
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
│   ├── cli.py                      # Command-line interface
│   ├── quality_assessment.py       # Module 2: QC (6 methods)
│   ├── profiler.py                 # Module 3: Profiling (32 features)
│   ├── recommender.py              # Module 4: Rule-based
│   ├── ml_recommender.py           # Module 4: ML-based
│   ├── de_import.py                # Module 7: DE import
│   ├── parameter_optimization.py   # Module 8: Optimization
│   ├── ensemble.py                 # Module 9: Ensemble
│   ├── simulation.py               # Simulation tools
│   │
│   ├── biomarker_discovery/        # Module 10: Biomarker Discovery (NEW in v2.2.2)
│   │   ├── __init__.py             # Public API re-exports
│   │   ├── core.py                 # discover_biomarkers, FeatureSelector,
│   │   │                           # PanelOptimizer, ClassifierEvaluator,
│   │   │                           # SurvivalAnalyzer, BiologicalAnnotator
│   │   ├── intent.py               # BiomarkerIntent (6 study-design modes)
│   │   ├── signature_score.py      # Signature score per patient
│   │   ├── direction_patterns.py   # UP/DOWN pattern + cross-cohort check
│   │   ├── clinical_metrics.py     # Youden, PPV/NPV, DCA, bootstrap CI
│   │   ├── ratio_biomarkers.py     # Self-normalizing ratio search (TSP)
│   │   ├── univariate_de.py        # Per-fold Welch / Mann-Whitney filter
│   │   └── enhanced.py             # EnhancedBiomarkerResult integration
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
│   │       ├── geo.py              # GEO connector (production-ready)
│   │       ├── sra.py              # SRA connector (production-ready)
│   │       ├── tcga.py             # TCGA connector (functional, multi-omic)
│   │       ├── arrayexpress.py     # ArrayExpress connector (in development)
│   │       ├── gene_mapping.py     # Gene ID conversion (MyGene.info)
│   │       ├── pooling.py          # Dataset pooling + batch correction
│   │       ├── datasets.py         # AcquiredDataset, PooledDataset
│   │       ├── cache.py            # Disk caching system
│   │       ├── catalog.py          # Dataset registry
│   │       └── base.py             # BaseConnector, SearchResult
│   │
│   ├── dashboard/                  # Interactive Streamlit app
│   │   ├── app.py
│   │   ├── pages/                  # Dashboard pages (10 pages)
│   │   │   ├── 01___Data_Acquisition.py        # NEW in v2.2.2
│   │   │   ├── 02_Quality_Assessment.py
│   │   │   ├── 03_Data_Profiler.py
│   │   │   ├── 04_Pipeline_Recommender.py
│   │   │   ├── 05_DE_Import_Compare.py
│   │   │   ├── 06_Optimization.py
│   │   │   ├── 07_Ensemble_Analysis.py
│   │   │   ├── 08_Biomarker_Discovery.py       # NEW in v2.2.2
│   │   │   ├── 09_Reports.py
│   │   │   ├── 10_Settings.py
│   │   │   └── 11_Visualization.py
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
│   ├── MODULE_8_Parameter_Optimization.md
│   ├── MODULE_9_Ensemble_Analysis.md
│   ├── MODULE_10_Biomarker_Discovery.md         # NEW in v2.2.2
│   ├── CONDA_ENVIRONMENTS.md
│   ├── RAPTOR_QUICK_REFERENCE.md                # Cheat sheet
│   └── RAPTOR_API_DOCUMENTATION.md              # Python API
│
├── examples/                       # Example scripts
│   ├── 02_quality_assessment.py
│   ├── 03_data_profiler.py
│   ├── 04_recommender.py
│   ├── 07_DE_Import.py
│   ├── 08_Parameter_Optimization.py
│   ├── 09_Ensemble_Analysis.py
│   └── 10_Biomarker_Discovery.py                # NEW in v2.2.2
│
├── tests/                          # Comprehensive test suite (700+ tests)
│   ├── test_acquisition.py                        # Module 6b
│   ├── test_biomarker_discovery.py                # Module 10 core
│   ├── test_biomarker_intent.py                   # Module 10 intent system
│   ├── test_biomarker_clinical_metrics.py         # Module 10 clinical
│   ├── test_biomarker_signature_score.py          # Module 10 signature
│   ├── test_biomarker_direction_patterns.py       # Module 10 direction
│   ├── test_biomarker_ratio.py                    # Module 10 ratios
│   ├── test_panel_size_kneedle.py                 # Module 10 kneedle
│   ├── test_cli_kneedle.py                        # CLI kneedle flags
│   ├── test_profiler.py
│   ├── test_quality_assessment.py
│   ├── test_parameter_optimization.py
│   ├── test_ensemble.py
│   ├── test_cli_comprehensive.py
│   └── ...
│
├── templates/                      # Sample sheets
│   ├── sample_sheet_paired.csv
│   └── sample_sheet_single.csv
│
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
| [MODULE_8_Parameter_Optimization.md](docs/MODULE_8_Parameter_Optimization.md) | 4 optimization methods |
| [MODULE_9_Ensemble_Analysis.md](docs/MODULE_9_Ensemble_Analysis.md) | 5 ensemble methods |
| [MODULE_10_Biomarker_Discovery.md](docs/MODULE_10_Biomarker_Discovery.md) | Nested-CV biomarker discovery, kneedle panel sizing, clinical metrics |

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
    print(f"Warning: {len(qc_report.outliers)} outliers detected")
else:
    print("No outliers detected")

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

# 6. Optimize Thresholds
print("\nStep 5: Optimize Thresholds...")
opt_result = optimize_with_fdr_control(deseq2, fdr_target=0.05)
print(f"  Optimal FDR: {opt_result.optimal_threshold['padj']:.3f}")
print(f"  Optimal |logFC|: {opt_result.optimal_threshold['lfc']:.3f}")

# 7. Ensemble Analysis
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
print("\nAnalysis complete!")
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

# Step 5: Optimize thresholds
raptor optimize --de-result imported/deseq2.csv --method fdr-control --fdr-target 0.05

# Step 6: Ensemble analysis
raptor ensemble-compare \
    --deseq2 imported/deseq2.csv \
    --edger imported/edger.csv \
    --limma imported/limma.csv \
    --output ensemble_results/

echo "Complete! Check ensemble_results/ for consensus genes."
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

### Example 5: Biomarker Discovery (Python API)

```python
from raptor.biomarker_discovery import discover_biomarkers

# Diagnostic biomarker discovery with default settings
# (kneedle panel-size detection, consensus pinning across folds,
#  univariate filter + Elastic Net + RFE feature selection,
#  4 classifiers, OOF clinical metrics, panel stability, annotation)
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    group_column="condition",
    intent="diagnostic",
    prevalence=0.05,                    # adjust to target population
    output_dir="results/biomarkers",
)

# The deployable panel
print(f"Panel ({len(result.panel)} genes): {result.panel}")
print(f"Best classifier: {result.best_classifier}")

# Honest out-of-fold performance
best = result.classification_results[result.best_classifier]
print(f"OOF AUC: {best.auc:.3f}")
if result.optimism_diagnostic:
    print(result.optimism_diagnostic.summary())

# Panel stability (Nogueira Phi)
if result.panel_stability is not None:
    print(result.panel_stability.summary())
    # Example output: "Phi=0.87 [0.74, 0.94] (excellent agreement, n=5 folds)"

# Clinical metrics (out-of-fold)
if result.clinical_metrics:
    cm = result.clinical_metrics
    print(f"Youden's threshold: {cm['youdens']['threshold']:.3f}")
    print(f"Sensitivity: {cm['youdens']['sensitivity']:.3f}")
    print(f"Specificity: {cm['youdens']['specificity']:.3f}")
    print(f"PPV at 5%: {cm['ppv_npv']['ppv']:.3f}")
    print(f"NPV at 5%: {cm['ppv_npv']['npv']:.3f}")
    print(f"AUC 95% CI: [{cm['bootstrap_ci']['ci_lower']:.3f}, "
          f"{cm['bootstrap_ci']['ci_upper']:.3f}]")
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
| **Module 10: Biomarker Discovery** | 2-15 min | 2-8 GB | Gene panel + OOF clinical metrics + Phi |

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

We are actively testing the Data Acquisition module (GEO, SRA, and TCGA search, dataset pooling, quality check). Your feedback on real-world datasets directly improves the tool.

- [Beta Testing Guide](BETA_TESTING_GUIDE.md) — 12 scenarios to try
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
| `raptor dashboard` CLI command fails with `ModuleNotFoundError` | Fix planned for v2.3.0 | Use `python -m raptor.launch_dashboard` or `python -m streamlit run raptor/dashboard/app.py` |
| TCGA advanced filtering | Functional with multi-omic support | Project, data type, sample type, and tumor/normal selection work; further filter granularity planned for v2.3.0 |
| ArrayExpress connector partially implemented | In development | Use GEO as alternative |
| GEO/SRA connectors require optional packages | By design | `pip install GEOparse biopython` — connectors fail gracefully without them |
| Gene ID conversion requires `mygene` | By design | `pip install mygene` — needed only for ID conversion feature |
| ComBat batch correction requires `combat` | By design | `pip install combat` — falls back to median-centering without it |
| ML Recommender requires pre-trained model file | By design | Use rule-based recommender (`--method rule`) |
| Module 10 prognostic intent not yet wired | Planned for v2.3.0 | Use `--intent diagnostic` or `--intent exploratory`; see [MODULE_10_Biomarker_Discovery.md](docs/MODULE_10_Biomarker_Discovery.md) for the intent system status table |
| Boruta, mRMR, SHAP, WGCNA, lifelines require optional packages | By design | `pip install raptor-rnaseq[biomarker]` to install all M10 optional dependencies; M10 falls back to the available subset of methods |

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

If you have multiple Python environments, an older RAPTOR install from a different environment may be on your PATH.

```bash
# Check the version
python -c "import raptor; print(raptor.__version__)"
# Should print: 2.2.2

# If you see an older version, upgrade:
pip install --upgrade raptor-rnaseq

# Or, if using a virtual environment, make sure it's activated:
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
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
