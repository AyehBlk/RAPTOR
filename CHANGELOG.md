# Changelog

All notable changes to RAPTOR will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.2.1] - 2026-03-18

### Dashboard Verified & Enhanced

All 9 Streamlit dashboard pages functionally tested with real RNA-seq data.

### Fixed
- **16 dashboard bugs** across Import DE, Ensemble, Optimization, Reports, and Quality pages
- **Column name mapping** between Module 7 (standardized: `log2_fold_change`, `p_value`) and Modules 8/9 (expected: `log2FoldChange`, `pvalue`)
- **Session state keys** aligned across all dashboard pages (`m7_de_results`, `m2_qc_report`, `m8_opt_result`)
- **QC report page** reads from correct nested dict keys (`m2_qc_report['overall']['score']`)
- **DE report** safely computes `n_up`/`n_down` from `results_df['direction']` column
- **Ensemble page** checks for `gene_id` column before `reset_index()`, removes invalid `filters` kwarg
- **Optimization page** reads `alpha` and `lfc_threshold` keys (not `lfc`/`log2fc`)
- **PCA plot** uses `log2(counts+1)` + `StandardScaler` instead of broken `assessor.data`
- **Correlation heatmap** applies `log2(counts+1)` before `.corr()` (was showing all 1.00)
- **Plotly deprecation** `titlefont` replaced with `title=dict(text=..., font=dict(size=...))`

### Enhanced
- **Quality Assessment page** expanded to 7 visualization tabs: PCA (2D/3D), sample correlation (3 methods), expression distribution, RLE plot, sample dendrogram, mean-variance/BCV
- **Optimization page** integrates all 4 scientific methods from Module 8 (FDR Control, Ground Truth, Stability, Reproducibility)
- **Visualization page** completely rewritten with 12 plot types and 7 gene expression styles (box, violin, beeswarm, raincloud, lollipop, heatmap bar, forest plot)
- **Professional styling** across all 9 pages (emoji cleanup, CSS cards, clean typography)
- **Recommender page** verified with rule-based recommendation (DESeq2 primary, 95% confidence)

### Added
- `raptor/cli.py` — all 55 CLI commands tested, version strings updated
- `tests/test_cli_comprehensive.py` — comprehensive CLI test suite
- `.gitignore` — proper Python/build/IDE ignore rules
- `.github/ISSUE_TEMPLATE/` — bug report, feature request, PR templates
- `.streamlit/config.toml` — dashboard theme configuration (RAPTOR green palette)
- `raptor/launch_dashboard.py` — pip-friendly dashboard launcher (`python -m raptor.launch_dashboard`)

### Changed
- Version bumped to 2.2.1 in `__init__.py`, `setup.py`, `cli.py` (15 locations)
- `use_column_width` replaced with `use_container_width` (Streamlit 1.28+ deprecation)
- Sidebar session state key changed from `m7_results` to `m7_de_results`

### Known Issues
- `raptor dashboard` CLI command broken (workaround: `python -m raptor.launch_dashboard`)
- `raptor --version` shows 2.2.0 (cosmetic; `cli.py` was built before version bump fix in source)
- ML Recommender inactive without pre-trained model file (rule-based works)

### Test Suite
- **413 passed, 0 failed, 50 expected skips** (unchanged from v2.2.0 baseline)

---

## [2.2.0] - 2026-03-12

### ⚠️ **Breaking Changes**

#### **Major Architecture Reorganization**

RAPTOR v2.2.0 introduces a significant restructuring to improve maintainability, add new features, and provide a cleaner API. Users upgrading from v2.1.x should review the [Migration Guide](docs/MIGRATION_GUIDE_v2.1_to_v2.2.md).

**Module Changes:**
- **Module 1 (Quick Count)** - REMOVED (was placeholder only)
- **Module 8 (Threshold Optimizer)** - RENAMED to "Parameter Optimization" with expanded functionality
  - Old: `raptor.threshold_optimizer.ATO`
  - New: `raptor.parameter_optimization` with 4 optimization methods
- **Module 9 (Ensemble Analysis)** - NEW module for combining DE results

**Pipeline Reorganization:**
- Consolidated from 8 pipelines to 6 validated, production-ready pipelines
- New structure: `raptor/pipelines/` with subdirectories per pipeline
- Removed experimental pipelines (2, 4, 6, 7, 8)
- Retained 6 core pipelines:
  1. Salmon (recommended)
  2. Kallisto (fastest)
  3. STAR + featureCounts (gene-level)
  4. STAR + RSEM (gene + isoform)
  5. STAR + Salmon (alignment + pseudo-alignment)
  6. HISAT2 + featureCounts (low memory)

**API Changes:**
```python
# OLD (v2.1.2)
from raptor.threshold_optimizer import ATO
optimizer = ATO()

# NEW (v2.2.0)
from raptor import optimize_with_fdr_control, optimize_with_ground_truth
```

---

### ✨ **New Features**

#### **Module 9: Ensemble Analysis**

Combine differential expression results from multiple methods (DESeq2, edgeR, limma) for robust, consensus gene lists.

**5 Ensemble Methods:**
1. **Fisher's Method** - Classic p-value combination
2. **Brown's Method** - Correlation-aware (recommended)
3. **Robust Rank Aggregation (RRA)** - Rank-based, outlier-robust
4. **Voting Consensus** - Conservative agreement-based
5. **Weighted Ensemble** - Use optimization-derived weights

**Usage:**
```python
from raptor import ensemble_brown

result = ensemble_brown(
    de_results={
        'deseq2': deseq2_result,
        'edger': edger_result,
        'limma': limma_result
    },
    significance_threshold=0.05
)
```

**CLI:**
```bash
raptor ensemble --methods brown --deseq2 deseq2.csv --edger edger.csv --limma limma.csv
raptor ensemble-compare --deseq2 deseq2.csv --edger edger.csv --limma limma.csv
```

---

#### **Module 8: Expanded Parameter Optimization**

Renamed from "Threshold Optimizer" and expanded from 1 to 4 optimization approaches.

**4 Optimization Methods:**

1. **Ground Truth Optimization** - Maximize F1 with known positives
2. **FDR Control Optimization** - Achieve target FDR (Storey's π₀)
3. **Stability-Based Optimization** - Cross-validation stability
4. **Reproducibility-Based Optimization** - Independent cohort agreement

**Usage:**
```python
from raptor import optimize_with_fdr_control

result = optimize_with_fdr_control(
    de_result=de_result,
    fdr_target=0.05
)
```

**CLI:**
```bash
raptor optimize --de-result results.csv --method fdr-control --fdr-target 0.05
```

---

#### **Module 3: Enhanced Data Profiling**

**32-Feature Data Profiler** for ML-based recommendations:

- 8 feature categories (sample, library, detection, expression, dispersion, sparsity, counts, mean-variance)
- **Key Metric: BCV (Biological Coefficient of Variation)**
  - BCV < 0.2: Low variation → limma-voom
  - BCV 0.2-0.4: Moderate → DESeq2, edgeR
  - BCV > 0.4: High → edgeR_robust

**Usage:**
```python
from raptor import profile_data_quick

profile = profile_data_quick('counts.csv', 'metadata.csv')
print(f"BCV: {profile.bcv:.3f}")
```

---

#### **Comprehensive Documentation**

**New Documentation Files:**
- `RAPTOR_USER_GUIDE.md` - Complete tutorial (~2000 lines)
- `RAPTOR_API_DOCUMENTATION.md` - Full Python API reference
- `RAPTOR_QUICK_REFERENCE.md` - Printable cheat sheet
- `CONDA_ENVIRONMENTS.md` - Conda setup guide
- `DOCUMENTATION_INDEX.md` - Documentation structure guide

**Updated Module Documentation:**
- MODULE_8_Parameter_Optimization.md (NEW)
- MODULE_9_Ensemble_Analysis.md (NEW)
- All existing module docs updated

---

#### **Conda Environment Support**

**Two Environment Options:**

1. **environment.yml** - Core (~500 MB, 5-10 min)
2. **environment-full.yml** - Complete suite (~5-8 GB, 30-60 min)

```bash
# Core
conda env create -f environment.yml
conda activate raptor

# Full (includes STAR, Salmon, Kallisto, R packages)
conda env create -f environment-full.yml
conda activate raptor-full
```

---

### 🔧 **Improvements**

#### **Code Quality**

- Custom exception hierarchy with actionable error messages
- Comprehensive input validation with type hints
- Expanded test suite with integration tests
- Better edge case handling

#### **Performance**

- Optimized dispersion estimation
- Faster p-value combination
- Improved memory efficiency

#### **User Experience**

- Consistent API across modules
- Informative progress messages
- Better error diagnostics
- Improved CLI help messages

#### **Dashboard Updates**

- Reorganized from 11 to 9 functional pages
- Removed Module 1 placeholder
- Updated module numbering
- Improved navigation and error handling

---

### 🐛 **Bug Fixes**

- Fixed dispersion estimation for small samples
- Corrected BCV calculation edge cases
- Resolved dashboard import errors
- Fixed metadata validation edge cases
- Corrected ensemble p-value combination for tied values
- Fixed pipeline initialization issues

---

### 🔄 **Deprecated**

- `raptor.threshold_optimizer.ATO` → Use `raptor.optimize_with_fdr_control()`
- Old pipeline structure (8 folders) → Use `raptor/pipelines/`
- Module 1 (Quick Count) → Removed (was placeholder)

---

### 🗑️ **Removed**

- **Module 1** - Quick Count (placeholder, never implemented)
- **Pipeline 2** - HISAT2 + StringTie + Ballgown (experimental)
- **Pipeline 4** - Kallisto + Sleuth
- **Pipeline 6** - STAR + featureCounts + NOISeq
- **Pipeline 7** - Bowtie2 + RSEM + EBSeq
- **Pipeline 8** - HISAT2 + Cufflinks + Cuffdiff

---

### 📦 **Dependencies**

#### **Updated**

- Python: `3.8-3.12` (expanded from `3.10`)
- numpy: `>=1.19.0` (was `>=1.21.0`)
- pandas: `>=1.1.0` (was `>=1.3.0`)
- scipy: `>=1.5.0` (was `>=1.7.0`)
- scikit-learn: `>=0.24.0` (was `>=1.0.0`)
- streamlit: `>=1.20.0` (was `>=1.28.0`)

**Reason:** Broader Python compatibility (3.8-3.12)

#### **New**

- `toml>=0.10.0` - TOML config support

#### **Optional** (install with `pip install raptor-rnaseq[all]`)

- bayesian-optimization, jinja2, markdown, psutil, colorlog

---

### 🎯 **Migration**

See [MIGRATION_GUIDE.md](docs/MIGRATION_GUIDE_v2.1_to_v2.2.md) for details.

**Quick Checklist:**
- [ ] Update imports: `threshold_optimizer` → `parameter_optimization`
- [ ] Review pipeline structure changes
- [ ] Explore Module 9 (Ensemble Analysis)
- [ ] Update conda environment
- [ ] Review Module 8 expansion

---

## [2.1.2] - 2025-01-XX

### Added
- Automated Threshold Optimization (ATO)
- ML-based pipeline recommendation
- Quality assessment with 6 outlier methods
- Data profiling

### Fixed
- Pipeline compatibility issues
- Documentation inconsistencies

---

## [2.1.1] - 2024-12-XX

### Added
- Initial ML recommender
- Basic profiler

### Fixed
- Windows installation issues
- Dependency conflicts

---

## [2.1.0] - 2024-11-XX

### Added
- First stable release
- 8 RNA-seq pipelines
- Basic recommendation
- CLI interface

---

## [Unreleased]

### Planned for v2.2.2
- Fix `raptor dashboard` CLI command
- Read version from `raptor.__version__` in `cli.py` (eliminate hardcoded strings)
- Update CITATION.cff DOI after Zenodo release

### Planned for v2.3.0
- Module 10: Biomarker Discovery
- Sphinx documentation website
- ReadTheDocs hosting
- GitHub Actions CI/CD
- Conda-forge package
- Docker containerization

---

**[2.2.1]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.1  
**[2.2.0]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.0  
**[2.1.2]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.1.2
