# 🦖 RAPTOR Scripts

Complete workflow scripts for RNA-seq pipeline benchmarking (v2.1.0).

## Numbered Workflow (00-10)

| Script | Language | Description |
|--------|----------|-------------|
| `00_simulate_data.R` | R | Generate simulated RNA-seq data with known DE genes |
| `01_run_all_pipelines.sh` | Bash | Run all 8 pipelines on FASTQ data |
| `01_run_all_pipelines_python.py` | Python | Run pipelines using RAPTOR benchmark module |
| `02_profile_data.py` | Python | Profile count data and get recommendations |
| `03_compare_results.R` | R | Compare DE results across pipelines |
| `04_visualize_comparison.R` | R | Create comparison visualizations |
| `05_ml_recommendation.py` | Python | ML-based pipeline recommendation |
| `06_quality_assessment.py` | Python | Comprehensive data quality scoring |
| `07_ensemble_analysis.py` | Python | Combine results from multiple pipelines |
| `08_automated_report.py` | Python | Generate publication-ready reports |
| `09_resource_monitoring.py` | Python | Track CPU/memory/disk usage |
| `10_parameter_optimization.py` | Python | Optimize analysis parameters |

## Utility Scripts

| Script | Description |
|--------|-------------|
| `demo.sh` | Interactive RAPTOR demonstration |
| `quick_profile.sh` | Fast profiling workflow |
| `full_benchmark.sh` | Complete pipeline benchmarking |

## Quick Start

```bash
# Run complete workflow
Rscript 00_simulate_data.R -o sim_data/
python 02_profile_data.py sim_data/counts.csv
Rscript 03_compare_results.R results/
Rscript 04_visualize_comparison.R results/

# Or use demo mode (no data required)
python 05_ml_recommendation.py --demo
python 06_quality_assessment.py --demo
```

## Requirements

**R packages:** polyester, Biostrings, rtracklayer, optparse, ggplot2, dplyr, tidyr, pheatmap, ggVennDiagram

**Python:** RAPTOR (`pip install raptor-rnaseq[ml]`)

---
*RAPTOR v2.1.0 - Making free science for everybody around the world 🌍*
