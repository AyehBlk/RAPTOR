# 🦖 RAPTOR Tests - v2.1.1

Test scripts for RAPTOR v2.1.1.

**🆕 New in v2.1.1:** Adaptive Threshold Optimizer (ATO) tests added!

---

##  Test Scripts

| Script | Type | Description |
|--------|------|-------------|
| `test_raptor_v2_1_1.py` | Integration | Comprehensive v2.1.1 tests (modules, version, dependencies, CLI, **ATO**) |
| `test_ml_system.py` | Unit | Tests ML recommender + **ATO** components |
| `test_recommender.py` | Unit | Pytest tests for PipelineRecommender + **ATO** classes |
| `test_threshold_optimizer.py` | Unit | **🆕 Comprehensive ATO tests** (pytest) |
| `test_pipelines.sh` | Integration | Pipeline installation and tool availability tests |

---

##  Quick Start

### Run Main Test (Recommended First)

```bash
python test_raptor_v2_1_1.py
```

Expected output:
```
 All tests passed! RAPTOR v2.1.1 is ready to use.
```

### Run All Tests

```bash
# Main comprehensive test
python test_raptor_v2_1_1.py

# ML system tests (includes ATO)
python test_ml_system.py

# Pytest unit tests
pytest test_recommender.py -v
pytest test_threshold_optimizer.py -v

# Shell integration tests
chmod +x test_pipelines.sh
./test_pipelines.sh
```

---

## 🆕 New in v2.1.1: ATO Tests

The Adaptive Threshold Optimizer (ATO) is a major new feature in v2.1.1. Test coverage includes:

### `test_threshold_optimizer.py` (Comprehensive pytest suite)

| Test Class | What It Tests |
|------------|---------------|
| `TestBasicFunctionality` | Initialization, column standardization, optimize() |
| `TestPValueAdjustment` | BH, BY, Holm, Bonferroni, Storey q-value |
| `TestLogFCOptimization` | MAD, mixture, power, percentile methods |
| `TestAnalysisGoals` | discovery, balanced, validation modes |
| `TestThresholdResult` | Result attributes, summary(), values |
| `TestHelperMethods` | get_significant_genes(), compare_thresholds() |
| `TestEdgeCases` | Small datasets, NaN handling, edge cases |

### Quick ATO Test

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load your DE results
df = pd.read_csv('deseq2_results.csv')

# Run ATO
result = optimize_thresholds(df, goal='balanced')

# Check results
print(f"LogFC cutoff: {result.logfc_cutoff}")
print(f"Padj cutoff: {result.padj_cutoff}")
print(f"DE genes: {result.n_significant_optimized}")
```

---

##  Test Coverage

### Module Imports (16 modules tested)

**Core v2.0.0:**
- `raptor.profiler` - RNAseqDataProfiler
- `raptor.recommender` - PipelineRecommender
- `raptor.benchmark` - PipelineBenchmark
- `raptor.simulate` - DataSimulator
- `raptor.report` - ReportGenerator
- `raptor.utils` - ensure_dir

**v2.1.0:**
- `raptor.ml_recommender` - MLPipelineRecommender
- `raptor.data_quality_assessment` - DataQualityAssessor
- `raptor.automated_reporting` - AutomatedReporter
- `raptor.ensemble_analysis` - EnsembleAnalyzer
- `raptor.parameter_optimization` - ParameterOptimizer
- `raptor.resource_monitoring` - ResourceMonitor
- `raptor.synthetic_benchmarks` - SyntheticBenchmarkGenerator

**🆕 v2.1.1:**
- `raptor.threshold_optimizer` - AdaptiveThresholdOptimizer, ThresholdResult, optimize_thresholds

### Dependencies Tested

**Required:**
- numpy, pandas, scipy, matplotlib, seaborn

**Optional:**
- scikit-learn (ML features)
- statsmodels (ATO π₀ estimation)
- boto3 (AWS)
- google-cloud-storage (GCP)

### Bioinformatics Tools

- samtools, fastqc
- STAR, Salmon, Kallisto, HISAT2

---

##  Understanding Test Output

### Successful Test

```
TEST 6: Threshold Optimizer (ATO) - NEW in v2.1.1
======================================================================

✓ ATO modules imported
✓ Generated test data: 1000 genes
✓ optimize_thresholds() executed
  - LogFC cutoff: 0.847
  - Padj cutoff: 0.05
  - Significant genes: 142
✓ Result is ThresholdResult type
✓ summary() method works
✓ Goal 'discovery': |logFC| > 0.623, padj < 0.1
✓ Goal 'balanced': |logFC| > 0.847, padj < 0.05
✓ Goal 'validation': |logFC| > 1.124, padj < 0.01
```

### Failed Test

```
✗ raptor.threshold_optimizer - Import Error: No module named 'raptor.threshold_optimizer'
```

If this happens, update RAPTOR: `pip install --upgrade raptor-rnaseq`

---

##  File Structure

```
tests/
├── README.md                     # This file
├── test_raptor_v2_1_1.py        # Main comprehensive tests
├── test_ml_system.py            # ML + ATO system tests
├── test_recommender.py          # Pytest: Recommender + ATO
├── test_threshold_optimizer.py  # 🆕 Pytest: Full ATO tests
└── test_pipelines.sh            # Shell integration tests
```


---

##  Requirements

```bash
# Install RAPTOR with ML features
pip install raptor-rnaseq[ml]

# Install test dependencies
pip install pytest

# Verify installation
python -c "import raptor; print(raptor.__version__)"
# Should output: 2.1.1
```

---

##  Troubleshooting

### "Module not found" Error

```bash
pip install --upgrade raptor-rnaseq
```

### "pytest not found" Error

```bash
pip install pytest
```

### ATO Tests Fail

Check that scipy is installed:
```bash
pip install scipy
```

### Version Mismatch

```bash
# Check current version
python -c "import raptor; print(raptor.__version__)"

# Update to 2.1.1
pip install raptor-rnaseq==2.1.1
```

---

*RAPTOR v2.1.1 - Making free science for everybody around the world 🌍*

*Author: Ayeh Bolouki*
