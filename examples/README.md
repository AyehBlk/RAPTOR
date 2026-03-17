# RAPTOR Example Scripts
## Complete User & Developer Guide (v2.2.0)

**Author:** Ayeh Bolouki  
**Package:** raptor-rnaseq (PyPI)  
**Updated:** 17 March 2025  

---

## Overview

This folder contains example scripts demonstrating the complete RAPTOR RNA-seq analysis workflow. Each script is self-contained with demo mode and can also accept real data.

### Script Inventory

| Script | Module | What It Does | Needs Tools? |
|--------|--------|-------------|--------------|
| `01_quick_count.py` | M1 | FASTQ → quick counts (Salmon/Kallisto) | Yes (Salmon or Kallisto) |
| `02_quality_assessment.py` | M2 | Count matrix → QC report (0-100 score) | No |
| `03_data_profiler.py` | M3 | Count matrix → 32 statistical features | No |
| `04_recommender.py` | M4 | Profile → pipeline recommendation | No |
| `05_production_pipeline.py` | M5 | FASTQ → production counts (STAR/HISAT2) | Yes (STAR, HISAT2, etc.) |
| `07_DE_Import.py` | M7 | R DE results → standardized format | No |
| `08_Parameter_Optimization.py` | M8 | DE results → optimized thresholds | No |
| `09_Ensemble_Analysis.py` | M9 | Multiple DE results → consensus genes | No |
| `10_qc_scenarios.py` | M2 | Stress test: 4 QC scenarios | No |
| `11_recommender_scenarios.py` | M4 | Stress test: 5 recommender scenarios | No |

> **Note:** M6 (R-based DE analysis) is performed in R, not Python.  
> **Note:** M10 (Biomarker Discovery) is planned for v3.0.0.

---

## Quick Start

### Prerequisites

```bash
# Install RAPTOR
pip install raptor-rnaseq

# Or editable install (for development)
cd RAPTOR-main
pip install -e .

# Verify
raptor --version
python -c "import raptor; print(raptor.__version__)"
```

### Windows Users

Set UTF-8 encoding before running any example (Windows terminals use cp1252 by default which cannot render Unicode characters in the scripts):

```powershell
$env:PYTHONIOENCODING = "utf-8"
```

### Run the Full Chain (Demo Mode)

No real data needed — each script generates simulated data internally:

```bash
# Stage 1: Fast Profiling (M2 → M3 → M4)
python examples/02_quality_assessment.py                    # Creates results/qc/
python examples/03_data_profiler.py                         # Creates results/profile.json
python examples/04_recommender.py --profile results/profile.json  # Creates results/recommendation.json

# Stage 3: DE Analysis (M7 → M8 → M9)
python examples/07_DE_Import.py --demo                      # Creates results/de_imported/
python examples/08_Parameter_Optimization.py --demo         # Creates results/parameter_optimization/
python examples/09_Ensemble_Analysis.py --demo              # Creates results/ensemble_analysis/
```

### Run with Real Data

```bash
# M2: QC with your count matrix
python examples/02_quality_assessment.py --counts my_counts.csv --metadata my_metadata.csv

# M3: Profile your data
python examples/03_data_profiler.py --counts results/qc/counts_clean.csv

# M4: Get recommendation based on your profile
python examples/04_recommender.py --profile results/profile.json

# M7: Import your R analysis results
python examples/07_DE_Import.py --de-file deseq2_results.csv --pipeline deseq2

# M8: Optimize parameters using M7 output
python examples/08_Parameter_Optimization.py --de-result results/de_imported/de_result.pkl

# M9: Combine multiple DE methods
python examples/09_Ensemble_Analysis.py --de-results results/de_imported/de_result.pkl
```

---

## Workflow Architecture

```
STAGE 1: Fast Profiling (Python only, ~5 min)
═══════════════════════════════════════════════
  M1: FASTQ files ──→ quick_gene_counts.csv     (requires Salmon/Kallisto)
           │
  M2: counts.csv ──→ qc_results.json             02_quality_assessment.py
           │           counts_clean.csv
  M3: counts.csv ──→ profile.json (32 features)  03_data_profiler.py
           │
  M4: profile.json ──→ recommendation.json        04_recommender.py
           │             (DESeq2/edgeR/limma)


STAGE 2: Production Quantification (requires tools, ~1-2 hrs)
═══════════════════════════════════════════════════════════════
  M5: FASTQ files ──→ gene_counts.csv            05_production_pipeline.py
                       tx_counts.csv


STAGE 3: DE Analysis (R + Python, ~30 min)
══════════════════════════════════════════
  M6: counts ──→ deseq2_results.csv               (run in R)
                  edger_results.csv
                  limma_results.csv
           │
  M7: R output ──→ de_standardized.csv             07_DE_Import.py
           │        de_result.pkl
  M8: de_result ──→ optimized_parameters.json      08_Parameter_Optimization.py
           │
  M9: multiple ──→ consensus_genes.csv             09_Ensemble_Analysis.py
       methods      method_comparison.csv
```

### Module Output Chain

Each module's output feeds into the next:

```
M2 output (counts_clean.csv) → M3 input
M3 output (profile.json)     → M4 input
M7 output (de_result.pkl)    → M8 input
M7 output (de_result.pkl)    → M9 input (multiple files)
```

---

## Script Details

### 02_quality_assessment.py (Module 2)

**Purpose:** Comprehensive QC scoring with 6 assessment components.

```bash
# Default (uses M1 output)
python examples/02_quality_assessment.py

# Custom input
python examples/02_quality_assessment.py --counts my_counts.csv --metadata my_meta.csv --output results/my_qc

# CLI equivalent
raptor qc --counts my_counts.csv --metadata my_meta.csv
```

**Output:** `results/qc/`
- `qc_results.json` — Full QC report with scores per component
- `counts_clean.csv` — Cleaned count matrix (outliers removed)

**QC Components:**
- Library Quality — library size variation, sequencing depth
- Gene Detection — detection rate, zero inflation
- Outlier Detection — Mahalanobis distance on PCA
- Variance Structure — PCA variance explained
- Batch Effects — batch-condition confounding (requires metadata)
- Biological Signal — DE signal strength

**Quality Score Interpretation:**
- 90-100: Excellent — proceed with confidence
- 70-89: Good — minor issues, safe to proceed
- 50-69: Acceptable — review warnings before proceeding
- <50: Poor — address issues before DE analysis

---

### 03_data_profiler.py (Module 3)

**Purpose:** Extract 32 statistical features for ML-based pipeline recommendation.

```bash
# Uses M2 output
python examples/03_data_profiler.py --counts results/qc/counts_clean.csv

# CLI equivalent
raptor profile --counts results/qc/counts_clean.csv
```

**Output:** `results/profile.json` — 32 features including:
- Sample characteristics (n_samples, n_genes, n_groups, balance)
- Library size metrics (mean, CV, range, depth category)
- Gene detection (detection rate, reliable expression)
- Expression distribution (mean, variance, skewness on log2 scale)
- Dispersion (BCV, common dispersion, overdispersion — **critical for M4**)
- Sparsity (zero proportion, zero inflation index)
- Mean-variance relationship (slope, Poisson/NB fit scores)
- Quality indicators (outliers, batch effects)

---

### 04_recommender.py (Module 4)

**Purpose:** Recommend the optimal DE analysis pipeline based on data characteristics.

```bash
# From profile (recommended)
python examples/04_recommender.py --profile results/profile.json

# From counts (auto-profiles first)
python examples/04_recommender.py --counts results/qc/counts_clean.csv

# CLI equivalent
raptor recommend --profile results/profile.json
```

**Output:** `results/recommendation.json` — Includes:
- Primary pipeline (DESeq2, edgeR, limma-voom, or Wilcoxon)
- Alternative pipeline
- Confidence scores (0-100) for all 5 pipelines
- Decision factors and warnings
- R code snippets for the recommended pipeline

**How the recommender decides:**

| Data Characteristic | Favors |
|---|---|
| Small samples (n < 8) | DESeq2 (best shrinkage) |
| Large samples (n > 20) | limma-voom (fast, powerful) |
| High dispersion / sparsity | edgeR (robust NB model) |
| Outlier samples detected | edgeR_robust |
| Very large samples (n > 50) | Wilcoxon (non-parametric) |

---

### 07_DE_Import.py (Module 7)

**Purpose:** Import and standardize DE results from R analysis tools.

```bash
# Demo mode
python examples/07_DE_Import.py --demo

# Import DESeq2 results
python examples/07_DE_Import.py --de-file deseq2_results.csv --pipeline deseq2

# Auto-detect pipeline format
python examples/07_DE_Import.py --de-file my_de_results.csv

# CLI equivalent
raptor import-de --de-file deseq2_results.csv
```

**Supported R formats:**

| Pipeline | Expected Columns |
|---|---|
| DESeq2 | gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj |
| edgeR | gene_id, logFC, logCPM, LR, PValue, FDR |
| limma-voom | gene_id, logFC, AveExpr, t, P.Value, adj.P.Val, B |
| Wilcoxon | gene_id, logFC, pvalue, padj |

**Output:** `results/de_imported/`
- `de_standardized.csv` — All genes, standardized column names
- `de_significant.csv` — Significant genes only (FDR < 0.05)
- `de_summary.json` — Summary statistics
- `de_result.pkl` — Serialized DEResult object for M8/M9

---

### 08_Parameter_Optimization.py (Module 8)

**Purpose:** Find optimal FDR and log2FC thresholds for DE analysis.

```bash
# Demo mode
python examples/08_Parameter_Optimization.py --demo

# From M7 output
python examples/08_Parameter_Optimization.py --de-result results/de_imported/de_result.pkl

# CLI equivalent
raptor optimize --de-result results/de_imported/de_result.pkl
```

**Optimization methods:**
- **FDR Control** — Statistical FDR estimation (Storey & Tibshirani, 2003)
- **Ground Truth** — Uses validated gene lists (if available)
- **Stability** — Cross-validation stability of gene lists
- **Reproducibility** — Agreement between cohorts

**Output:** `results/parameter_optimization/`
- `optimization_results.json` — Optimal alpha and LFC thresholds

**Example result:** alpha=0.03, LFC=0.25 means use `results(dds, alpha=0.03, lfcThreshold=0.25)` in R.

---

### 09_Ensemble_Analysis.py (Module 9)

**Purpose:** Combine DE results from multiple methods into high-confidence consensus gene lists.

```bash
# Demo mode
python examples/09_Ensemble_Analysis.py --demo

# CLI equivalent
raptor ensemble --de-results results/de_imported/de_result.pkl
```

**5 Ensemble methods:**

| Method | Approach | Best For |
|---|---|---|
| **Voting** | Gene must pass ≥2 methods | High-confidence publication lists |
| **Fisher** | Combine p-values (Fisher, 1925) | Maximum sensitivity / exploration |
| **Brown** | Correlation-aware p-value combination | Conservative, accounts for method similarity |
| **RRA** | Robust Rank Aggregation (Kolde, 2012) | Robust to outlier methods |
| **Weighted** | Performance-weighted combination | When you have validation data |

**Output:** `results/ensemble_analysis/`
- `voting/consensus_genes.csv`
- `fisher/consensus_genes.csv`
- `brown/consensus_genes.csv`
- `rra/consensus_genes.csv`
- `weighted/consensus_genes.csv`
- `method_comparison.csv` — Side-by-side comparison
- `ensemble_summary.json`

---

### 10_qc_scenarios.py (Stress Test)

**Purpose:** Verify RAPTOR QC handles different RNA-seq data quality problems.

```bash
python examples/10_qc_scenarios.py
```

**4 Scenarios tested:**

| Scenario | What It Simulates | Expected Detection |
|---|---|---|
| Gold Standard | 8 samples, 200 DE genes, balanced batches | Score > 85, all green |
| Outlier + Library Issues | 1 sample at 10%, library sizes vary 5x | Library quality warning |
| Confounded Batch | Batch1 = Controls, Batch2 = Treatments | Confounding warning |
| Sparse Data | 70% zeros, 4 samples | Zero inflation, low signal |

---

### 11_recommender_scenarios.py (Stress Test)

**Purpose:** Verify the recommender adapts to different data characteristics.

```bash
python examples/11_recommender_scenarios.py
```

**5 Scenarios tested:**

| Scenario | Expected Primary | Why |
|---|---|---|
| Small samples (n=3) | DESeq2 | Best shrinkage for few samples |
| Large samples (n=20) | limma-voom | Fast, powerful at large n |
| High sparsity (60% zeros) | edgeR | Handles sparse data best |
| Outlier present | edgeR_robust | Robust to outlier samples |
| Confounded batch | Any + warning | Must warn about confounding |

---

## CLI Equivalents

Every example script has a CLI equivalent using `raptor` commands:

```bash
# Module 2: Quality Assessment
raptor qc --counts counts.csv --metadata metadata.csv --output results/qc

# Module 3: Data Profiling
raptor profile --counts counts.csv --output results

# Module 4: Pipeline Recommendation
raptor recommend --counts counts.csv
raptor recommend --profile results/profile.json

# Module 5: Pipeline Management
raptor pipeline list
raptor pipeline run --pipeline salmon --sample-sheet samples.csv --index /path/to/index

# Module 7: DE Import
raptor import-de --de-file deseq2_results.csv --pipeline deseq2

# Module 8: Parameter Optimization
raptor optimize --de-result results/de_imported/de_result.pkl

# Module 9: Ensemble Analysis
raptor ensemble --de-results results/de_imported/de_result.pkl

# General
raptor --version
raptor --help
```

---

## Input File Formats

### Count Matrix (CSV)

```csv
,Sample1,Sample2,Sample3,Sample4
ENSG00000000001,150,230,180,210
ENSG00000000002,45,62,38,55
ENSG00000000003,0,3,1,2
```

- First column: gene IDs (Ensembl, gene symbol, or any identifier)
- Remaining columns: raw integer counts per sample
- No negative values
- Minimum: 10 genes, 2 samples

### Metadata (CSV)

```csv
sample_id,condition,batch
Sample1,Control,Batch1
Sample2,Control,Batch1
Sample3,Treatment,Batch2
Sample4,Treatment,Batch2
```

- `sample_id`: must match column names in count matrix
- `condition`: experimental groups (required for M4 recommendation)
- `batch`: batch information (optional, used for batch effect detection)

### DE Results from R (CSV)

See [Module 7 section](#07_de_importpy-module-7) for expected column formats per pipeline.

---

## Output Directory Structure

After running all examples:

```
results/
├── quick_counts/              # M1 output
│   ├── quick_gene_counts.csv
│   ├── quick_tpm.csv
│   └── sample_info.csv
├── qc/                        # M2 output
│   ├── qc_results.json
│   └── counts_clean.csv
├── profile.json               # M3 output
├── recommendation.json        # M4 output
├── de_imported/               # M7 output
│   ├── de_standardized.csv
│   ├── de_significant.csv
│   ├── de_summary.json
│   └── de_result.pkl
├── parameter_optimization/    # M8 output
│   └── fdr_control/
│       └── optimization_results.json
└── ensemble_analysis/         # M9 output
    ├── voting/consensus_genes.csv
    ├── fisher/consensus_genes.csv
    ├── brown/consensus_genes.csv
    ├── rra/consensus_genes.csv
    ├── weighted/consensus_genes.csv
    ├── method_comparison.csv
    └── ensemble_summary.json
```

---

## Troubleshooting

### Windows: UnicodeEncodeError (cp1252)

```
UnicodeEncodeError: 'charmap' codec can't encode characters
```

**Fix:** Set UTF-8 encoding before running:
```powershell
$env:PYTHONIOENCODING = "utf-8"
```

### Module Not Found / Import Error

```
ModuleNotFoundError: No module named 'raptor'
```

**Fix:** Install in editable mode:
```bash
cd RAPTOR-main
pip install -e .
```

### Count File Not Found

```
Count file not found: results/quick_counts/quick_gene_counts.csv
```

**Fix:** Either run M1 first, or create simulated data:
```python
import pandas as pd, numpy as np, os
os.makedirs('results/quick_counts', exist_ok=True)
np.random.seed(42)
counts = pd.DataFrame(
    np.random.negative_binomial(10, 0.05, (1000, 6)),
    index=[f'Gene_{i}' for i in range(1000)],
    columns=[f'Sample_{i}' for i in range(1, 7)]
)
counts.to_csv('results/quick_counts/quick_gene_counts.csv')
```

### --demo Flag

Scripts M7, M8, M9 require `--demo` flag when running without real DE results:
```bash
python examples/07_DE_Import.py --demo
python examples/08_Parameter_Optimization.py --demo
python examples/09_Ensemble_Analysis.py --demo
```

---

## For Developers

### Adding a New Example Script

1. Name it `XX_description.py` where XX is the module number
2. Include demo mode (`--demo` flag) that generates simulated data
3. Use argparse with `--help` documentation
4. Import from `raptor` package (with try/except fallback for demo mode)
5. Save output to `results/` following the directory convention
6. Print a completion summary with next steps
7. Test on both Windows (PowerShell) and Linux

### Testing All Examples

```bash
# Set encoding (Windows only)
$env:PYTHONIOENCODING = "utf-8"

# Run all examples
python examples/02_quality_assessment.py
python examples/03_data_profiler.py
python examples/04_recommender.py --profile results/profile.json
python examples/07_DE_Import.py --demo
python examples/08_Parameter_Optimization.py --demo
python examples/09_Ensemble_Analysis.py --demo

# Run stress tests
python examples/10_qc_scenarios.py
python examples/11_recommender_scenarios.py
```

### Common Example Script Bugs

These patterns caused bugs in v2.2.0 example scripts — avoid in new scripts:

| Bug Pattern | Example | Fix |
|---|---|---|
| Wrong dataclass field names | `EnsembleResult(method=...)` | Check `ClassName.__dataclass_fields__.keys()` |
| `to_dict()` on a dict | `result.to_dict() if hasattr(...)` | Check `isinstance(result, dict)` first |
| Wrong attribute names in display | `getattr(profile, 'lib_size_mean')` | Check actual DataProfile fields |
| numpy int64 in json.dump | `json.dump(data)` crashes | Always use `json.dump(data, default=str)` |
| Windows path in assertion | `'results/quick_counts' in str(path)` | Use `'quick_counts' in str(path)` |

---

## Version History

- **v2.2.0** (March 2025): Initial release with all example scripts
- **v2.2.0-fix** (March 2025): Fixed display bugs in 02, 03, 08, 09 example scripts
- Scripts 10, 11 added for stress testing QC and recommender modules

---

*RAPTOR: Making free science for everybody around the world*
