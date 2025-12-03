# 🦖 RAPTOR Examples

Example scripts demonstrating RAPTOR v2.1.0 features.

## Official Example Workflows

| Script | Description |
|--------|-------------|
| `example_ml_workflow.py` | Complete ML training workflow (data generation → training → evaluation) |
| `example_quality_assessment.py` | Four comprehensive quality assessment examples |

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

## Quick Start

```bash
# Test without RAPTOR installed
python 05_ml_recommendation.py --demo
python 06_quality_assessment.py --demo

# With real data (RAPTOR required)
python 05_ml_recommendation.py --counts counts.csv --metadata metadata.csv
python 06_quality_assessment.py --counts counts.csv --plot
```

## Installation

```bash
pip install raptor-rnaseq[ml]
```

---
*RAPTOR v2.1.0 - Making free science for everybody around the world 🌍*
