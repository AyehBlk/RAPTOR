# 🦖 RAPTOR Quick Start Guide

## What's Included

This package contains a complete **Machine Learning-based Pipeline Recommendation System** for RAPTOR:

### Core Features

- **ML-Based Recommendations** - Intelligent pipeline selection with 87% accuracy
- **🎯 Threshold Optimizer** - Data-driven significance thresholds (**NEW in v2.1.1!**)
- **Interactive Dashboard** - Web-based interface (no coding required!)
- **Quality Assessment** - Comprehensive data quality scoring
- **Ensemble Analysis** - Combine multiple pipelines for robust results
- **Resource Monitoring** - Track CPU, memory, and runtime
- **Automated Reporting** - Publication-ready reports

---

## 🆕 What's New in v2.1.1

### Adaptive Threshold Optimizer (ATO)

Stop using arbitrary thresholds! ATO determines data-driven cutoffs:

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load your DE results
df = pd.read_csv('deseq2_results.csv')

# Get optimized thresholds
result = optimize_thresholds(df, goal='discovery')

print(f"Optimal logFC: {result.logfc_threshold:.2f}")
print(f"Significant genes: {result.n_significant}")
print(f"\nMethods text for publication:\n{result.methods_text}")
```

---

## Installation (2 Minutes)

### Option 1: Install from PyPI (Recommended)

```bash
pip install raptor-rnaseq
```

With all features:
```bash
pip install raptor-rnaseq[all]
```

### Option 2: Install from GitHub

```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -r requirements.txt
```

### Verify Installation

```bash
python -c "import raptor; print(raptor.__version__)"
```

Expected output:
```
2.1.1
```

### Verify Threshold Optimizer

```bash
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

---

## First Run (2 Minutes)

### Option A: Interactive Dashboard (Easiest)

```bash
raptor dashboard
```

Then open http://localhost:8501 in your browser. Upload your data and get instant recommendations!

**NEW in v2.1.1:** Click "🎯 Threshold Optimizer" in the sidebar to optimize DE thresholds!

### Option B: Command Line

```bash
# Get ML recommendation for your data
raptor profile --counts your_data.csv --use-ml
```

### Option C: Python API

```python
from raptor import RNAseqDataProfiler, MLPipelineRecommender
import pandas as pd

# Load your data
counts = pd.read_csv('counts.csv', index_col=0)

# Profile data
profiler = RNAseqDataProfiler(counts)
profile = profiler.run_full_profile()

# Get ML recommendation
recommender = MLPipelineRecommender()
recommendation = recommender.recommend(profile)

print(f"Recommended: Pipeline {recommendation['pipeline_id']}")
print(f"Confidence: {recommendation['confidence']:.1%}")
```

---

## Quick Threshold Optimization (NEW in v2.1.1!)

### From Python

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
import pandas as pd

# Load DE results from any pipeline
df = pd.read_csv('deseq2_results.csv')

# Create optimizer
ato = AdaptiveThresholdOptimizer(
    df, 
    logfc_col='log2FoldChange', 
    pvalue_col='pvalue'
)

# Optimize for discovery (exploratory analysis)
result = ato.optimize(goal='discovery')

# Or for validation (stringent analysis)
result = ato.optimize(goal='validation')

# Get results
print(f"Optimal logFC threshold: {result.logfc_threshold:.3f}")
print(f"Optimal p-value threshold: {result.pvalue_threshold}")
print(f"Significant genes: {result.n_significant}")
print(f"π₀ estimate: {result.pi0:.3f}")

# Get publication-ready methods text
print(result.methods_text)

# Save optimized results
result.results_df.to_csv('optimized_results.csv')
```

### From Dashboard

1. Launch: `raptor dashboard`
2. Click "🎯 Threshold Optimizer" in sidebar
3. Upload your DE results (CSV)
4. Select analysis goal (Discovery/Balanced/Validation)
5. Click "Optimize Thresholds"
6. Download results and methods text

---

## Basic Usage

### 1. Get Pipeline Recommendation

```bash
# Quick recommendation
raptor profile --counts data.csv --use-ml

# Output:
# 🦖 RECOMMENDED: Pipeline 3 (Salmon-edgeR)
# Confidence: 89%
# Reason: Optimal for your sample size (n=12) and moderate BCV (0.35)
```

### 2. Run Quality Assessment

```python
from raptor.data_quality_assessment import DataQualityAssessor

assessor = DataQualityAssessor(counts, metadata)
report = assessor.assess_quality()

print(f"Quality Score: {report['overall_score']}/100")
```

### 3. Optimize Thresholds (NEW!)

```python
from raptor.threshold_optimizer import optimize_thresholds

# Quick optimization
result = optimize_thresholds(de_results, goal='balanced')
print(f"Use |logFC| > {result.logfc_threshold:.2f}, padj < {result.pvalue_threshold}")
```

### 4. Ensemble Analysis

```python
from raptor.ensemble_analysis import EnsembleAnalyzer

analyzer = EnsembleAnalyzer()
consensus = analyzer.combine_results(
    results_dict={'deseq2': df1, 'edger': df2, 'limma': df3},
    method='weighted_vote'
)

print(f"Consensus DE genes: {len(consensus['de_genes'])}")
```

### 5. Generate Report

```bash
raptor report --results results/ --output report.html
```

---

## Quick Command Reference

```bash
# Install
pip install raptor-rnaseq

# Launch dashboard
raptor dashboard

# Get recommendation
raptor profile --counts data.csv --use-ml

# Run pipeline
raptor run --pipeline 3 --data fastq/ --output results/

# Generate report
raptor report --results results/ --output report.html

# Show help
raptor --help
```

---

## Analysis Goals for Threshold Optimizer

| Goal | Use Case | Behavior |
|------|----------|----------|
| **discovery** | Exploratory analysis | More permissive, maximize sensitivity |
| **balanced** | Standard publication | Balance sensitivity/specificity |
| **validation** | Clinical/confirmation | Stringent, maximize specificity |

```python
# Examples
result = optimize_thresholds(df, goal='discovery')   # More genes
result = optimize_thresholds(df, goal='balanced')    # Standard
result = optimize_thresholds(df, goal='validation')  # Fewer, high-confidence genes
```

---

## Directory Structure After Installation

```
your_project/
├── data/
│   ├── counts.csv          # Your count matrix
│   └── metadata.csv        # Sample metadata
├── results/                # Pipeline outputs
│   └── deseq2_results.csv  # DE results (use with ATO!)
├── reports/                # Generated reports
└── figures/                # Visualizations
```

---

## Common Use Cases

### Use Case 1: Quick Recommendation

```bash
raptor profile --counts data.csv --use-ml
```

### Use Case 2: Full Workflow with Threshold Optimization

```bash
# 1. Profile and get recommendation
raptor profile --counts data.csv --use-ml

# 2. Run the recommended pipeline
raptor run --pipeline 3 --data fastq/ --output results/

# 3. Optimize thresholds for the results
python -c "
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd
df = pd.read_csv('results/de_results.csv')
result = optimize_thresholds(df, goal='balanced')
print(f'Optimal thresholds: |logFC| > {result.logfc_threshold:.2f}, padj < {result.pvalue_threshold}')
result.results_df.to_csv('results/optimized_de_results.csv')
"

# 4. Generate report
raptor report --results results/
```

### Use Case 3: Batch Processing

```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler
from raptor.threshold_optimizer import optimize_thresholds
from pathlib import Path

recommender = MLPipelineRecommender()

for data_file in Path('datasets/').glob('*.csv'):
    counts = pd.read_csv(data_file, index_col=0)
    profiler = RNAseqDataProfiler(counts)
    profile = profiler.run_full_profile()
    
    rec = recommender.recommend(profile)
    print(f"{data_file.name}: Pipeline {rec['pipeline_id']}")
```

---

## Troubleshooting

### "Module not found"
```bash
pip install raptor-rnaseq[all]
```

### "threshold_optimizer not found"
```bash
pip install --upgrade raptor-rnaseq
```

### "Model not found"
```bash
# Train a model first
python example_ml_workflow.py --n-datasets 100
```

### "Low confidence predictions"
- Check data quality with `DataQualityAssessor`
- Ensure sufficient sample size (n ≥ 6)
- Consider ensemble analysis

---

## Need Help?

- **Documentation**: https://github.com/AyehBlk/RAPTOR/tree/main/docs
- **PyPI**: https://pypi.org/project/raptor-rnaseq/
- **Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

## Success Criteria

✓ Installation completes without errors
✓ `import raptor` works
✓ Dashboard launches at http://localhost:8501
✓ ML recommendation runs in <0.1 seconds
✓ Threshold Optimizer available (**v2.1.1**)
✓ Confidence scores provided

---

**Ready to start? Run:**
```bash
pip install raptor-rnaseq
python -c "import raptor; from raptor.threshold_optimizer import optimize_thresholds; print('✅ Ready!')"
```

🦖 Happy analyzing!

---

**Version:** 2.1.1  
**Author:** Ayeh Bolouki  
**License:** MIT
