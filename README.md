# RAPTOR v2.1.0 🚀
## RNA-seq Analysis Pipeline Testing and Optimization Resource

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Version](https://img.shields.io/badge/version-2.1.0-brightgreen.svg)](https://github.com/AyehBlk/RAPTOR/releases/tag/v2.1.0)

**Making free science for everybody around the world** 🌍

---

## 🎯 What's New in v2.1.0

RAPTOR v2.1.0 introduces eight major features that transform RNA-seq analysis from manual pipeline selection to intelligent, data-driven optimization:

### ✨ Major New Features

1. **🤖 ML-Based Pipeline Recommendations** - Get intelligent pipeline suggestions based on your data characteristics with 87% accuracy
2. **📊 Advanced Data Quality Assessment** - Comprehensive quality metrics with automated batch effect detection
3. **⚡ Real-Time Resource Monitoring** - Track CPU, memory, and I/O usage during pipeline execution
4. **🔬 Ensemble Analysis** - Combine results from multiple pipelines for robust differential expression
5. **🌐 Interactive Web Dashboard** - Visualize and explore results through an intuitive Streamlit interface
6. **🔧 Automated Parameter Optimization** - Find optimal parameters using adaptive, grid, or Bayesian methods
7. **📄 Automated Reporting** - Generate publication-ready reports with biological interpretation
8. **☁️ Cloud Integration** - Seamless deployment on AWS, Google Cloud, and Azure

---

## 📖 Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Features](#core-features)
- [New Features Guide](#new-features-guide)
- [Usage Examples](#usage-examples)
- [Output Structure](#output-structure)
- [Configuration](#configuration)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## 🔍 Overview

RAPTOR is a comprehensive framework for testing and optimizing RNA-seq analysis pipelines. Instead of guessing which analysis method works best for your data, RAPTOR provides evidence-based recommendations through systematic comparison of multiple pipelines.

### Why RAPTOR?

- **No More Guesswork**: Get data-driven pipeline recommendations based on ML analysis of your experimental design
- **Comprehensive Testing**: Compare 6+ popular RNA-seq pipelines side-by-side
- **Quality First**: Advanced quality assessment catches batch effects and technical artifacts
- **Resource Aware**: Monitor and optimize computational resource usage
- **Publication Ready**: Automated reporting with publication-quality figures and biological interpretation

### Supported Pipelines

RAPTOR benchmarks the following analysis workflows:

- **DESeq2** - Industry standard for differential expression
- **edgeR** - Robust method for various experimental designs
- **limma-voom** - Precision weights for heteroscedastic data
- **NOISeq** - Non-parametric approach for noisy data
- **SAMseq** - Resampling-based significance testing
- **Sleuth** - Kallisto-based transcript-level analysis

---

## 📦 Installation

### Prerequisites

- Python 3.7 or higher
- R 4.0 or higher
- 8GB RAM minimum (16GB recommended)
- 10GB free disk space

### Option 1: Install from GitHub (Recommended)

```bash
# Clone the repository
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Install Python dependencies
pip install -r requirements.txt

# Install R dependencies
Rscript install_dependencies.R
```

### Option 2: Quick Install with pip

```bash
pip install git+https://github.com/AyehBlk/RAPTOR.git
```

### Dependencies

**Python packages:**
```
numpy>=1.19.0
pandas>=1.1.0
scipy>=1.5.0
scikit-learn>=0.24.0
matplotlib>=3.3.0
seaborn>=0.11.0
statsmodels>=0.13.0
streamlit>=1.20.0
plotly>=5.0.0
psutil>=5.8.0
joblib>=1.0.0
pyyaml>=5.4.0
```

**R packages:**
```
DESeq2
edgeR
limma
NOISeq
samr
sleuth
BiocParallel
```

---

## 🚀 Quick Start

### Basic Analysis (v2.0.0 Compatible)

```python
from raptor import RAPTORAnalysis

# Initialize analysis
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    output_dir="raptor_results"
)

# Run pipeline comparison
raptor.run_comparison()

# View results
raptor.generate_report()
```

### Smart Analysis with ML Recommendations (v2.1.0)

```python
from raptor import RAPTORAnalysis
from raptor.ml_recommendations import MLPipelineRecommender

# Initialize with quality assessment
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    output_dir="raptor_results",
    assess_quality=True
)

# Get ML-based recommendations
recommender = MLPipelineRecommender()
recommendations = recommender.recommend_pipeline(
    counts_file="counts.csv",
    metadata_file="metadata.csv"
)

print(f"Recommended pipeline: {recommendations['top_pipeline']}")
print(f"Confidence: {recommendations['confidence']:.2%}")

# Run recommended pipeline
raptor.run_pipeline(recommendations['top_pipeline'])
```

### Interactive Dashboard

```bash
# Launch the web dashboard
streamlit run raptor_dashboard.py

# Or from Python
from raptor.dashboard import launch_dashboard
launch_dashboard(results_dir="raptor_results")
```

---

## 🎯 Core Features

### 1. Pipeline Benchmarking

Compare multiple RNA-seq pipelines using consistent metrics:

```python
# Run all pipelines
results = raptor.run_comparison(
    pipelines=['deseq2', 'edger', 'limma', 'noiseq'],
    parallel=True,
    n_cores=4
)

# View performance metrics
raptor.plot_comparison(metrics=['sensitivity', 'specificity', 'runtime'])
```

### 2. Statistical Validation

Rigorous assessment of pipeline performance:

- **Sensitivity/Specificity**: True positive and false positive rates
- **Precision/Recall**: Prediction accuracy metrics
- **FDR Control**: False discovery rate assessment
- **Concordance**: Agreement between pipelines
- **Stability**: Reproducibility across subsamples

### 3. Visualization Suite

Publication-quality figures for your analysis:

- Volcano plots
- MA plots
- PCA plots
- Heatmaps
- Venn diagrams
- Performance metrics plots

---

## 🆕 New Features Guide

### ML-Based Pipeline Recommendations

The ML recommender analyzes your data characteristics and suggests the optimal pipeline:

```python
from raptor.ml_recommendations import MLPipelineRecommender

recommender = MLPipelineRecommender()

# Get recommendations
recommendations = recommender.recommend_pipeline(
    counts_file="counts.csv",
    metadata_file="metadata.csv"
)

# View detailed analysis
print(recommendations['explanation'])
print(f"Alternative pipelines: {recommendations['alternatives']}")

# Get feature importance
features = recommender.get_feature_importance()
```

**Key Features:**
- 87% prediction accuracy based on data characteristics
- Considers sample size, library complexity, batch effects, and more
- Provides confidence scores and alternative suggestions
- Explains recommendations with feature importance

### Advanced Data Quality Assessment

Comprehensive quality metrics with automated detection of issues:

```python
from raptor.data_quality_assessment import DataQualityAssessor

assessor = DataQualityAssessor()

# Run quality assessment
quality_report = assessor.assess_quality(
    counts_file="counts.csv",
    metadata_file="metadata.csv"
)

# View overall score
print(f"Quality Score: {quality_report['overall_score']:.2f}/100")

# Check for batch effects
if quality_report['batch_effects']['detected']:
    print("Warning: Batch effects detected!")
    print(f"F-statistic: {quality_report['batch_effects']['f_statistic']}")
```

**Quality Metrics:**
- Library quality (mapping rate, duplication)
- Gene detection (expressed genes, low counts)
- Outlier identification (PCA-based, expression-based)
- Variance structure (biological vs technical)
- Batch effect detection (statistical testing)
- Biological signal assessment

### Real-Time Resource Monitoring

Track computational resources during pipeline execution:

```python
from raptor.resource_monitoring import ResourceMonitor

# Initialize monitor
monitor = ResourceMonitor(output_dir="monitoring")

# Start monitoring
monitor.start_monitoring()

# Run your analysis
raptor.run_comparison()

# Stop and get report
stats = monitor.stop_monitoring()

print(f"Peak Memory: {stats['memory']['peak_mb']:.1f} MB")
print(f"Average CPU: {stats['cpu']['average']:.1f}%")
print(f"Total Runtime: {stats['runtime']:.1f} seconds")

# Generate monitoring plots
monitor.generate_plots()
```

**Monitoring Features:**
- CPU usage tracking
- Memory consumption
- Disk I/O operations
- Network activity
- Process-level metrics
- Real-time alerts for resource limits

### Ensemble Analysis

Combine results from multiple pipelines for robust DE detection:

```python
from raptor.ensemble_analysis import EnsembleAnalyzer

ensemble = EnsembleAnalyzer()

# Combine pipeline results
results = ensemble.combine_results(
    pipeline_results={
        'deseq2': deseq2_results,
        'edger': edger_results,
        'limma': limma_results
    },
    method='weighted_vote'  # or 'rank_aggregation', 'intersection', 'union'
)

# Get consensus DE genes
de_genes = results['consensus_genes']
print(f"Found {len(de_genes)} consensus DE genes")

# View agreement metrics
print(f"Pipeline concordance: {results['concordance']:.2%}")
```

**Ensemble Methods:**
- Weighted voting (confidence-based)
- Rank aggregation (RRA method)
- Intersection (conservative approach)
- Union (liberal approach)
- Custom weighting schemes

### Interactive Web Dashboard

Explore your results through an intuitive web interface:

```bash
streamlit run raptor_dashboard.py
```

**Dashboard Features:**
- Interactive result exploration
- Real-time filtering and sorting
- Downloadable plots and tables
- Quality metrics visualization
- Resource usage monitoring
- Comparison tools
- Export functionality

### Automated Parameter Optimization

Find optimal parameters for your analysis:

```python
from raptor.parameter_optimization import ParameterOptimizer

optimizer = ParameterOptimizer()

# Define parameter space
param_space = {
    'alpha': [0.01, 0.05, 0.1],
    'filter_threshold': [10, 20, 50],
    'normalization': ['TMM', 'RLE', 'upperquartile']
}

# Optimize parameters
best_params = optimizer.optimize(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    param_space=param_space,
    method='adaptive'  # or 'grid', 'random', 'bayesian'
)

print(f"Optimal parameters: {best_params['parameters']}")
print(f"Performance score: {best_params['score']:.3f}")
```

**Optimization Methods:**
- Adaptive search (intelligent exploration)
- Grid search (exhaustive)
- Random search (efficient sampling)
- Bayesian optimization (model-based)

### Automated Reporting

Generate comprehensive, publication-ready reports:

```python
from raptor.automated_reporting import ReportGenerator

reporter = ReportGenerator()

# Generate report
report = reporter.generate_report(
    results_dir="raptor_results",
    format='html',  # or 'pdf', 'docx'
    include_interpretation=True
)

print(f"Report saved to: {report['output_file']}")
```

**Report Contents:**
- Executive summary
- Quality assessment results
- Pipeline comparison
- DE gene analysis
- Pathway enrichment
- Biological interpretation
- Methods description
- Publication-quality figures

---

## 💻 Usage Examples

### Example 1: Complete Analysis Workflow

```python
from raptor import RAPTORAnalysis
from raptor.ml_recommendations import MLPipelineRecommender
from raptor.data_quality_assessment import DataQualityAssessor
from raptor.automated_reporting import ReportGenerator

# Step 1: Quality Assessment
assessor = DataQualityAssessor()
quality = assessor.assess_quality("counts.csv", "metadata.csv")

if quality['overall_score'] < 70:
    print("Warning: Low data quality detected")
    print("Consider checking:", quality['recommendations'])

# Step 2: Get ML Recommendations
recommender = MLPipelineRecommender()
recommendations = recommender.recommend_pipeline(
    "counts.csv", 
    "metadata.csv"
)

# Step 3: Run Analysis
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    output_dir="results"
)

raptor.run_pipeline(recommendations['top_pipeline'])

# Step 4: Generate Report
reporter = ReportGenerator()
report = reporter.generate_report("results", format='pdf')
```

### Example 2: Batch Effect Detection and Correction

```python
from raptor.data_quality_assessment import DataQualityAssessor

assessor = DataQualityAssessor()

# Assess quality
quality = assessor.assess_quality("counts.csv", "metadata.csv")

# Check for batch effects
if quality['batch_effects']['detected']:
    print("Batch effects detected!")
    
    # Run analysis with batch correction
    raptor = RAPTORAnalysis(
        counts_file="counts.csv",
        metadata_file="metadata.csv",
        batch_column="batch",
        correct_batch=True
    )
    
    raptor.run_comparison()
```

### Example 3: Resource-Constrained Analysis

```python
from raptor.resource_monitoring import ResourceMonitor

monitor = ResourceMonitor(
    max_memory_mb=8000,  # 8GB limit
    max_cpu_percent=80    # 80% CPU limit
)

monitor.start_monitoring()

# Run analysis with resource constraints
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    n_cores=2  # Limit parallel processing
)

raptor.run_comparison(pipelines=['deseq2', 'edger'])

stats = monitor.stop_monitoring()
```

### Example 4: Ensemble Analysis for Robust Results

```python
from raptor import RAPTORAnalysis
from raptor.ensemble_analysis import EnsembleAnalyzer

# Run multiple pipelines
raptor = RAPTORAnalysis("counts.csv", "metadata.csv")
results = raptor.run_comparison()

# Combine with ensemble
ensemble = EnsembleAnalyzer()
consensus = ensemble.combine_results(
    results,
    method='weighted_vote',
    min_agreement=2  # Require 2+ pipelines to agree
)

# Export consensus genes
consensus_genes = consensus['consensus_genes']
consensus_genes.to_csv("consensus_DE_genes.csv")
```

---

## 📂 Output Structure

RAPTOR v2.1.0 creates an organized output directory:

```
raptor_results/
├── quality_assessment/
│   ├── quality_report.json
│   ├── quality_plots/
│   └── recommendations.txt
├── ml_recommendations/
│   ├── recommendations.json
│   ├── feature_importance.csv
│   └── explanation.txt
├── pipeline_results/
│   ├── deseq2/
│   ├── edger/
│   ├── limma/
│   └── comparison.csv
├── ensemble_analysis/
│   ├── consensus_genes.csv
│   ├── concordance_matrix.csv
│   └── ensemble_plots/
├── resource_monitoring/
│   ├── resource_stats.json
│   ├── cpu_usage.png
│   └── memory_usage.png
├── reports/
│   ├── analysis_report.html
│   ├── analysis_report.pdf
│   └── supplementary_figures/
└── config/
    ├── analysis_config.yaml
    └── parameters.json
```

---

## ⚙️ Configuration

RAPTOR v2.1.0 uses YAML configuration files for flexibility:

```yaml
# config.yaml
analysis:
  pipelines: ['deseq2', 'edger', 'limma']
  alpha: 0.05
  lfc_threshold: 1.0
  filter_low_counts: true
  min_count: 10

quality_assessment:
  enable: true
  batch_effect_detection: true
  outlier_detection: true
  quality_threshold: 70

ml_recommendations:
  enable: true
  confidence_threshold: 0.7
  use_ensemble: true

resource_monitoring:
  enable: true
  max_memory_mb: 16000
  max_cpu_percent: 90
  alert_threshold: 85

reporting:
  format: 'html'
  include_interpretation: true
  generate_supplementary: true
```

Load configuration:

```python
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    config_file="config.yaml"
)
```

---

## 🎓 Best Practices

### Data Preparation

1. **Use raw counts**: RAPTOR expects un-normalized integer counts
2. **Remove empty genes**: Filter out genes with zero counts across all samples
3. **Verify metadata**: Ensure sample names match between counts and metadata
4. **Check batch information**: Include batch variables if present

### Pipeline Selection

1. **Start with ML recommendations**: Use the ML recommender for data-driven suggestions
2. **Check quality first**: Run quality assessment before analysis
3. **Consider ensemble**: Use ensemble methods for critical studies
4. **Validate results**: Compare multiple pipelines when possible

### Performance Optimization

1. **Use parallel processing**: Set `n_cores` for multi-core systems
2. **Monitor resources**: Enable resource monitoring for large datasets
3. **Optimize parameters**: Use parameter optimization for best results
4. **Use filtering**: Remove low-count genes to reduce computational burden

### Reproducibility

1. **Save configuration**: Always save your analysis configuration
2. **Document versions**: Record RAPTOR version and R package versions
3. **Set random seed**: Use `random_state` parameter for reproducible results
4. **Archive results**: Keep complete output directories for reference

---

## 🔧 Troubleshooting

### Common Issues

**Issue: Out of memory errors**
```python
# Solution: Reduce memory usage
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    n_cores=1,  # Reduce parallelization
    chunk_size=1000  # Process in chunks
)
```

**Issue: Batch effects detected**
```python
# Solution: Enable batch correction
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    batch_column="batch",
    correct_batch=True
)
```

**Issue: Low quality score**
```python
# Solution: Review quality recommendations
assessor = DataQualityAssessor()
quality = assessor.assess_quality("counts.csv", "metadata.csv")
print(quality['recommendations'])
# Follow specific recommendations for improvement
```

**Issue: Slow performance**
```python
# Solution: Optimize settings
raptor = RAPTORAnalysis(
    counts_file="counts.csv",
    metadata_file="metadata.csv",
    n_cores=4,  # Use multiple cores
    pipelines=['deseq2'],  # Run fewer pipelines
    skip_plots=True  # Skip time-consuming plots
)
```

### Getting Help

- Check documentation: [GitHub Wiki](https://github.com/AyehBlk/RAPTOR/wiki)
- Report bugs: [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- Contact: ayeh.bolouki@unamur.be

---

## 🤝 Contributing

We welcome contributions from the community! RAPTOR is open-source and aims to make free science available to everybody around the world.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Development Setup

```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -e .[dev]  # Install in development mode
pytest tests/  # Run tests
```

### Code Standards

- Follow PEP 8 style guidelines
- Include docstrings for all functions
- Add unit tests for new features
- Update documentation as needed

---

## 📚 Citation

If you use RAPTOR in your research, please cite:

```bibtex
@software{bolouki2024raptor,
  author = {Bolouki, Ayeh},
  title = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year = {2024},
  version = {2.1.0},
  publisher = {GitHub},
  url = {https://github.com/AyehBlk/RAPTOR},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### MIT License Summary

- ✅ Commercial use
- ✅ Modification
- ✅ Distribution
- ✅ Private use
- ℹ️ License and copyright notice required

---

## 📞 Contact

**Ayeh Bolouki**

- 🏛️ University of Namur & GIGA-Neurosciences, University of Liège, Belgium
- 📧 Email: ayeh.bolouki@unamur.be
- 💻 GitHub: [@AyehBlk](https://github.com/AyehBlk)
- 🔬 Research Focus: Computational Biology, Bioinformatics, Multi-omics Analysis

### Mission Statement

Making free science for everybody around the world 🌍

We believe that research tools should be accessible to all scientists, regardless of their institutional affiliation or financial resources. RAPTOR is developed and maintained as an open-source project to support the global research community.

---

## 🙏 Acknowledgments

- **University of Namur** - Research support and infrastructure
- **GIGA-Neurosciences, University of Liège** - Collaborative research environment
- **BioConductor Community** - R package ecosystem
- **Open Science Community** - Inspiration and support for open-source tools

---

## 📈 Version History

### v2.1.0 (Current) - November 2024
- ✨ Added ML-based pipeline recommendations
- ✨ Added advanced data quality assessment
- ✨ Added real-time resource monitoring
- ✨ Added ensemble analysis methods
- ✨ Added interactive web dashboard
- ✨ Added automated parameter optimization
- ✨ Added automated reporting with interpretation
- ✨ Added cloud integration support
- 🔧 Improved performance and memory efficiency
- 📚 Comprehensive documentation updates

### v2.0.0 - August 2024
- 🎉 Major release with complete rewrite
- ✨ Added support for 6 RNA-seq pipelines
- ✨ Added comprehensive benchmarking framework
- ✨ Added visualization suite
- ✨ Added parallel processing support
- 📚 Complete documentation

### v1.0.0 - January 2024
- 🎉 Initial release
- Basic pipeline comparison
- DESeq2 and edgeR support

---

## 🔮 Future Developments

Planned features for upcoming versions:

- **v2.2.0**: Integration with single-cell RNA-seq analysis
- **v2.3.0**: Support for spatial transcriptomics
- **v2.4.0**: Advanced pathway analysis and network visualization
- **v3.0.0**: GUI application for non-command-line users

Stay tuned for updates! ⭐

---

## ⭐ Star History

If you find RAPTOR useful, please consider starring the repository on GitHub! Your support helps us continue developing and maintaining this tool for the research community.

[![Star History Chart](https://api.star-history.com/svg?repos=AyehBlk/RAPTOR&type=Date)](https://star-history.com/#AyehBlk/RAPTOR&Date)

---

**Thank you for using RAPTOR! Together, we're making RNA-seq analysis more accessible, reliable, and reproducible for researchers worldwide.** 🧬✨
