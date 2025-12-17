# 🦖 RAPTOR Scripts

Complete workflow scripts for RNA-seq pipeline benchmarking (v2.1.1).

## 🆕 New in v2.1.1

**Adaptive Threshold Optimizer (ATO)** - Data-driven threshold selection:

```bash
# Optimize thresholds for your DE results
python 11_threshold_optimizer.py --input deseq2_results.csv --goal balanced

# Compare different analysis goals
python 11_threshold_optimizer.py --input results.csv --compare-goals

# Demo mode (no data required)
python 11_threshold_optimizer.py --demo
```

---

## Numbered Workflow (00-11)

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
| `11_threshold_optimizer.py` | Python | **🆕 Adaptive threshold optimization (ATO)** |

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

# NEW: Optimize thresholds for DE results
python 11_threshold_optimizer.py --input deseq2_results.csv

# Or use demo mode (no data required)
python 05_ml_recommendation.py --demo
python 06_quality_assessment.py --demo
python 11_threshold_optimizer.py --demo
```

##  Threshold Optimizer Usage

The new `11_threshold_optimizer.py` script provides data-driven threshold selection:

```bash
# Basic usage with DESeq2 output
python 11_threshold_optimizer.py --input deseq2_results.csv

# Different analysis goals
python 11_threshold_optimizer.py --input results.csv --goal discovery   # Liberal
python 11_threshold_optimizer.py --input results.csv --goal balanced    # Default
python 11_threshold_optimizer.py --input results.csv --goal validation  # Conservative

# Compare all goals side-by-side
python 11_threshold_optimizer.py --input results.csv --compare-goals

# With visualization
python 11_threshold_optimizer.py --input results.csv --plot --output ato_results/

# Custom column names (for edgeR output)
python 11_threshold_optimizer.py --input edger.csv --logfc-col logFC --pvalue-col PValue
```

### Supported Input Formats

- **DESeq2**: `log2FoldChange`, `pvalue`
- **edgeR**: `logFC`, `PValue`
- **limma**: `logFC`, `P.Value`
- **NOISeq**: `log2FC`, `prob`

### Output

- `ato_results.json` - Optimization results with methods text
- `significant_genes_ato.csv` - Genes passing optimized thresholds
- `ato_optimization_summary.png` - Visualization (with `--plot`)

## Requirements

**R packages:** polyester, Biostrings, rtracklayer, optparse, ggplot2, dplyr, tidyr, pheatmap, ggVennDiagram

**Python:** RAPTOR (`pip install raptor-rnaseq[ml]`)

## Complete Workflow Example

```bash
# 1. Simulate data (or use your own)
Rscript 00_simulate_data.R -o sim_data/ --n-genes 10000 --n-de 1500

# 2. Profile your data
python 02_profile_data.py sim_data/counts.csv -o profile_results/

# 3. Get ML-based pipeline recommendation
python 05_ml_recommendation.py --counts sim_data/counts.csv --output recommendations/

# 4. Run pipelines (using recommended ones)
bash 01_run_all_pipelines.sh sim_data/ results/ references/ 1,3,5

# 5. Compare DE results
Rscript 03_compare_results.R results/

# 6. Visualize comparison
Rscript 04_visualize_comparison.R results/

# 7. Ensemble analysis
python 07_ensemble_analysis.py --input results/ --output ensemble/

# 8. Optimize thresholds (NEW in v2.1.1!)
python 11_threshold_optimizer.py --input results/deseq2/de_results.csv --goal balanced

# 9. Generate final report
python 08_automated_report.py --input results/ --output final_report/
```

---

*RAPTOR v2.1.1 - Making free science for everybody around the world 🌍*

*Author: Ayeh Bolouki*
