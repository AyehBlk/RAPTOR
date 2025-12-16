# 🦖 RAPTOR Examples

Example scripts demonstrating RAPTOR v2.1.1 features.

## 🆕 New in v2.1.1

**Adaptive Threshold Optimizer (ATO)** - Data-driven threshold selection for DE analysis:

```bash
python example_threshold_optimizer.py
```

---

## Official Example Workflows

| Script | Description |
|--------|-------------|
| `example_ml_workflow.py` | Complete ML training workflow (data generation → training → evaluation) |
| `example_quality_assessment.py` | Four comprehensive quality assessment examples |
| `example_threshold_optimizer.py` | **NEW!** ATO demonstration (7 examples of threshold optimization) |

## Feature Demonstrations (with --demo mode)

All these scripts support `--demo` mode for testing without RAPTOR installed:

| Script | Feature | Demo Command |
|--------|---------|--------------|
| `05_ml_recommendation.py` | ML-based pipeline recommendation | `python 05_ml_recommendation.py --demo` |
| `06_quality_assessment.py` | Data quality scoring | `python 06_quality_assessment.py --demo` |
| `07_ensemble_analysis.py` | Multi-pipeline consensus | `python 07_ensemble_analysis.py --demo` |
| `08_automated_report.py` | Report generation | `python 08_automated_report.py --demo` |
| `09_resource_monitoring.py` | Resource tracking | `python 09_resource_monitoring.py --demo` |
| `10_parameter_optimization.py` | Parameter tuning | `python 10_parameter_optimization.py --demo` |

##  Threshold Optimizer Examples

The new `example_threshold_optimizer.py` includes 7 examples:

1. **Basic Usage** - Quick optimization with convenience function
2. **Full Control** - AdaptiveThresholdOptimizer class with custom settings
3. **Threshold Comparison** - Compare different logFC/padj combinations
4. **Adjustment Methods** - Compare BH, BY, Holm, Bonferroni, qvalue
5. **Visualization** - Generate publication-ready plots
6. **Real Data** - Working with DESeq2/edgeR/limma output
7. **Analysis Goals** - Compare discovery/balanced/validation modes

```bash
# Run all ATO examples
python example_threshold_optimizer.py

# Or import in your own code
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(de_results, goal='balanced')
print(result.methods_text)  # Ready for your paper!
```

## Quick Start

```bash
# Test without RAPTOR installed
python 05_ml_recommendation.py --demo
python 06_quality_assessment.py --demo

# ATO example (runs with synthetic data)
python example_threshold_optimizer.py

# With real data (RAPTOR required)
python 05_ml_recommendation.py --counts counts.csv --metadata metadata.csv
python 06_quality_assessment.py --counts counts.csv --plot
```

## Installation

```bash
pip install raptor-rnaseq[ml]
```

## Example Output Structure

```
examples/
├── README.md                           # This file
├── example_ml_workflow.py              # ML training workflow
├── example_quality_assessment.py       # Quality assessment examples
├── example_threshold_optimizer.py      # 🆕 ATO examples
├── 05_ml_recommendation.py             # ML recommendation demo
├── 06_quality_assessment.py            # Quality scoring demo
├── 07_ensemble_analysis.py             # Ensemble analysis demo
├── 08_automated_report.py              # Report generation demo
├── 09_resource_monitoring.py           # Resource monitoring demo
└── 10_parameter_optimization.py        # Parameter optimization demo
```

## ATO Quick Reference

```python
from raptor.threshold_optimizer import optimize_thresholds

# Load your DE results
import pandas as pd
df = pd.read_csv('deseq2_results.csv')

# Optimize thresholds
result = optimize_thresholds(
    df,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'  # or 'discovery' or 'validation'
)

# Get results
print(f"Optimal |logFC|: {result.logfc_threshold:.3f}")
print(f"Significant genes: {result.n_significant}")

# Get publication methods text
print(result.methods_text)
```

---

*RAPTOR v2.1.1 - Making free science for everybody around the world 🌍*

*Author: Ayeh Bolouki*
