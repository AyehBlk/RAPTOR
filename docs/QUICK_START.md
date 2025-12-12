# 🦖 RAPTOR Quick Start Guide

## What's Included

This package contains a complete **Machine Learning-based Pipeline Recommendation System** for RAPTOR:

### Core Features

- **ML-Based Recommendations** - Intelligent pipeline selection with 87% accuracy
- **Interactive Dashboard** - Web-based interface (no coding required!)
- **Quality Assessment** - Comprehensive data quality scoring
- **Ensemble Analysis** - Combine multiple pipelines for robust results
- **Resource Monitoring** - Track CPU, memory, and runtime
- **Automated Reporting** - Publication-ready reports

## Installation (2 Minutes)

### Option 1: Install from PyPI (Recommended)

```bash
pip install raptor-rnaseq
```

With all features:
```bash
pip install raptor-rnaseq[all]
```

With specific features:
```bash
pip install raptor-rnaseq[dashboard]   # Web dashboard
pip install raptor-rnaseq[ml]          # ML features
pip install raptor-rnaseq[advanced]    # Advanced features
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
2.1.0
```

Or run the test suite:
```bash
python test_ml_system.py
```

## First Run (2 Minutes)

### Option A: Interactive Dashboard (Easiest)

```bash
python launch_dashboard.py
```

Then open http://localhost:8501 in your browser. Upload your data and get instant recommendations!

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

### 3. Ensemble Analysis

```python
from raptor.ensemble_analysis import EnsembleAnalyzer

analyzer = EnsembleAnalyzer()
consensus = analyzer.combine_results(
    results_dict={'deseq2': df1, 'edger': df2, 'limma': df3},
    method='weighted_vote'
)

print(f"Consensus DE genes: {len(consensus['de_genes'])}")
```

### 4. Generate Report

```bash
raptor report --results results/ --output report.html
```

## Quick Command Reference

```bash
# Install
pip install raptor-rnaseq

# Launch dashboard
python launch_dashboard.py

# Get recommendation
raptor profile --counts data.csv --use-ml

# Run pipeline
raptor run --pipeline 3 --data fastq/ --output results/

# Generate report
raptor report --results results/ --output report.html

# Show help
raptor --help
```

## Directory Structure After Installation

```
your_project/
├── data/
│   ├── counts.csv          # Your count matrix
│   └── metadata.csv        # Sample metadata
├── results/                # Pipeline outputs
├── reports/                # Generated reports
└── figures/                # Visualizations
```

## Common Use Cases

### Use Case 1: Quick Recommendation

```bash
raptor profile --counts data.csv --use-ml
```

### Use Case 2: Full Workflow

```bash
# 1. Profile and get recommendation
raptor profile --counts data.csv --use-ml

# 2. Run the recommended pipeline
raptor run --pipeline 3 --data fastq/ --output results/

# 3. Generate report
raptor report --results results/
```

### Use Case 3: Batch Processing

```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler
from pathlib import Path

recommender = MLPipelineRecommender()

for data_file in Path('datasets/').glob('*.csv'):
    counts = pd.read_csv(data_file, index_col=0)
    profiler = RNAseqDataProfiler(counts)
    profile = profiler.run_full_profile()
    
    rec = recommender.recommend(profile)
    print(f"{data_file.name}: Pipeline {rec['pipeline_id']}")
```

## Troubleshooting

### "Module not found"
```bash
pip install raptor-rnaseq[all]
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

## Need Help?

- **Documentation**: https://github.com/AyehBlk/RAPTOR/tree/main/docs
- **PyPI**: https://pypi.org/project/raptor-rnaseq/
- **Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

## Success Criteria

✓ Installation completes without errors
✓ `import raptor` works
✓ Dashboard launches at http://localhost:8501
✓ ML recommendation runs in <0.1 seconds
✓ Confidence scores provided

---

**Ready to start? Run:**
```bash
pip install raptor-rnaseq
python -c "import raptor"
```

🦖 Happy analyzing!
