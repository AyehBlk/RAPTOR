# Changelog

All notable changes to RAPTOR will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.2.2] - 2026-04-28

This release adds two major subsystems to RAPTOR.

**Module 6b — Data Acquisition** lands first (work completed mid-April):
search, download, and pool RNA-seq datasets from public repositories
(GEO, SRA, TCGA, ArrayExpress) with automatic gene-ID conversion, batch
correction, and quality-check tooling. Dashboard-only.

**Module 10 — Biomarker Discovery** lands on top: nested-CV-based
biomarker panel selection with multi-method feature selection,
panel-size auto-detection via the kneedle algorithm with consensus
pinning across CV folds, classification evaluation, out-of-fold
clinical metrics, biological annotation, and honest small-cohort
reporting. Exposed via both the `raptor biomarker` CLI command and the
Biomarker Discovery dashboard page.

### Added — Module 10: Biomarker Discovery

#### Pipeline structure

- Multi-method feature selection: LASSO (elastic net), Boruta, mRMR,
  Recursive Feature Elimination (RFE), SHAP, WGCNA. Genes ranked by
  consensus across selected methods.
- Significance calibration of the consensus ranking: down-weights
  genes whose univariate two-sample test (Welch's t-test or
  Mann-Whitney U) returned a non-significant p-value before they
  enter the panel.
- Panel-size optimization with kneedle algorithm (see below) and
  consensus pinning across CV folds.
- Classification evaluation in nested CV across four classifiers
  (Random Forest, SVM, XGBoost, Logistic Regression) with
  deterministic tiebreak when AUCs are tied.
- Biological annotation via MyGene.info (gene symbols, descriptions,
  pathway membership), STRING PPI queries, and PubMed literature
  mining with disease-context filtering.
- Publication-ready report generation including ranked genes,
  classification performance, panel curve, and biomarker_result.pkl.

#### Leakage fix (Ambroise & McLachlan 2002 PNAS)

Feature selection now runs inside each outer CV fold rather than on
all data prior to CV. The previous external-CV-based feature
selection was giving optimistically biased AUC estimates. Per-fold
selection eliminates this leakage. Tests verify the per-fold
behavior.

#### Panel stability (Nogueira 2018 JMLR)

Computes the Nogueira stability index Phi across the per-fold
panels, with bootstrap 95% CI and benchmark labels (excellent /
intermediate / poor). Surfaces visually on the Panel tab of the
dashboard.

#### Honest performance reporting

- Optimism diagnostic: cross-validated AUC vs. training-data AUC
  with the gap reported. Large gap flags overfitting.
- Bootstrap 95% CI on the cross-validated AUC.
- Out-of-fold clinical metrics: Youden's optimal threshold,
  sensitivity, specificity, PPV/NPV at user-supplied prevalence,
  decision curve analysis. All computed on out-of-fold predictions
  rather than in-sample fits.

#### Direction patterns

- Separates UP-regulated and DOWN-regulated biomarkers in the
  selected panel for biological interpretation.

#### Kneedle panel-size auto-detection

- Replaces the legacy first-drop heuristic with the kneedle
  algorithm (Satopaa et al. 2011), with polynomial smoothing for
  noise robustness on small-n curves and argmax fallback for
  saturated curves with no detectable knee.
- Three strategies exposed: `kneedle` (default), `argmax`,
  `first_drop` (legacy, retained for backward compatibility).
- Per-fold size variability controlled via a discovery-level
  consensus scout: when `panel_size_strategy='consensus'` (default),
  panel-size detection runs once on the full data and all outer CV
  folds are pinned to that K. The `per_fold` strategy (each outer
  fold detects independently) is retained as an opt-in.
- Consensus made the default after a sanity check on real fixtures
  showed per-fold detection drove a 0.313 Phi drop on a 60-sample
  cohort under seed 42 (per-fold sizes (3,3,3,12,3) inflated apparent
  panel-stability variance without changing the discovery-level K).

#### Small-cohort tooling (dashboard-facing)

- Tiered advisory banner at the top of the Clinical Metrics tab:
  strong advisory below n=20, moderate advisory below n=50, no
  advisory at n>=50.
- "Verify across multiple seeds" button on the Panel tab: re-runs
  discovery at additional random seeds and reports which diagnostics
  are stable across seeds (trust these) vs. which vary (treat with
  caution). Surfaces single-seed-CV unreliability at small n.

#### Biomarker intent classification

- Six intent modes: diagnostic, prognostic, predictive, monitoring,
  exploratory, translational. Each configures study design,
  validation strategy, required metadata columns, minimum sample
  count, and required output metrics.
- `--intent` CLI flag and intent dropdown on the dashboard
  auto-configure enhanced analyses (signature score, direction
  patterns, clinical metrics, ratio biomarkers).

#### Ratio biomarkers

- Generalized Top Scoring Pair search (Geman et al. 2004): tests all
  pairwise gene ratios, ranks by AUC, returns the top-k discriminating
  pairs. Self-normalizing: batch effects, library-size differences,
  and platform variation cancel within ratios.

#### Files added

- `raptor/biomarker_discovery/` — 7-file subpackage:
  `core.py`, `intent.py`, `signature_score.py`,
  `direction_patterns.py`, `clinical_metrics.py`,
  `ratio_biomarkers.py`, `enhanced.py`, `univariate_de.py`.
- `raptor/dashboard/pages/08_🧬_Biomarker_Discovery.py` — full
  dashboard page with 10 tabs (Panel, Ranked Genes, Classification,
  Panel Curve, Signature Score, Direction Pattern, Clinical Metrics,
  Ratio Biomarkers, Annotations, Downloads).
- `docs/MODULE_10_Biomarker_Discovery.md` — module documentation.

#### CLI additions

- `raptor biomarker` — full discovery command with options for
  feature-selection methods, panel size, species, disease term,
  annotation toggles.
- `raptor biomarker-survival` — survival-context biomarker discovery
  (stub for prognostic intent, full implementation deferred).
- `raptor biomarker-validate` — independent-cohort validation of a
  discovered panel.
- New flags on `raptor biomarker`:
  `--intent [diagnostic|exploratory]`,
  `--prevalence FLOAT` (PPV/NPV calculation),
  `--auto-panel-strategy [kneedle|argmax|first_drop]`,
  `--panel-sensitivity FLOAT`,
  `--panel-size-strategy [per_fold|consensus]`.

#### Dashboard pages renumbered

- Module 10 (Biomarker Discovery) inserted as page 08. Reports,
  Settings, and Visualization shifted down by one position to
  accommodate. Page contents unchanged; only filenames renumbered.

#### Tests added

- `tests/test_biomarker_pipeline_cv.py` — feature-selection-inside-CV.
- `tests/test_biomarker_optimism_diagnostic_m2.py` — CV-vs-training
  AUC + bootstrap CI.
- `tests/test_biomarker_m4_calibration.py` — significance calibration.
- `tests/test_biomarker_rfe_ranking.py` — recursive feature elimination.
- `tests/test_biomarker_clinical_metrics_oof.py` — out-of-fold
  clinical metrics.
- `tests/test_biomarker_best_classifier_selection.py` — tiebreak
  determinism.
- `tests/test_dashboard_panel_stability.py` — Phi color/interpretation,
  small-cohort advisory tiers, multi-seed verification helpers.
- `tests/test_panel_size_kneedle.py` — kneedle algorithm on synthetic
  curves (clean concave, saturated, monotone, short), forward
  selection integration, dashboard text-scan tests, consensus-pinned
  label regression.
- `tests/test_cli_kneedle.py` — CLI flag advertisement, invalid value
  rejection, end-to-end selection_method verification, explicit
  --panel-size override.
- `tests/test_m10_cli_scenarios.py` — end-to-end CLI scenarios across
  multi-method, small-data, upstream integration, output completeness,
  classification quality, validation workflow.
- `tests/test_biomarker_intent.py` — intent mode validation.
- `tests/test_biomarker_discovery.py` — overall pipeline integration.
- `tests/generate_edge_case_data.py` and
  `tests/generate_realistic_test_data.py` — fixture generators.

### Added — Module 6b: Data Acquisition

Complete data acquisition subpackage for searching, downloading, and
pooling datasets from public genomics repositories. Accessible via the
dashboard (Data Acquisition page) — no CLI commands needed.

- **GEO Connector** — search GEO via NCBI Entrez, download expression
  matrices, extract processing metadata (aligner, assembly,
  quantification, library selection), fuzzy sample alignment with
  greedy optimal assignment.
- **SRA Connector** — search SRA via ENA `read_run` endpoint,
  study-level grouping, rich run tables, cross-platform FASTQ
  download scripts (PowerShell `.ps1` / Bash `.sh`),
  `find_linked_gse` method for GSM→GSE lookup.
- **TCGA Connector** — GDC API integration with multi-omic support
  (gene expression, miRNA, methylation, CNV, RPPA), data-type-aware
  dashboard labels, `tar.gz` recovery for truncated downloads.
- **ArrayExpress Connector** — BioStudies API integration (in
  development).
- **Gene ID Mapper** — detect ID type (Ensembl/Symbol/Entrez),
  convert via MyGene.info, mean aggregation for duplicates, failed
  gene download.
- **Pooling Engine** — inner/outer gene merging, batch correction
  (ComBat, quantile normalization, median ratio), study-label
  tracking.
- **Data Library** — repository-aware display separating GEO
  expression datasets from SRA run tables, gene ID conversion UI,
  sample metadata editor with batch auto-population.
- **Quality Check tab** — library size distributions, PCA, sample
  correlations, expression distributions, RLE plots for pooled data.
- **Cache Manager** — parquet-based caching for fast dataset reload.
- **Data Catalog** — track downloaded/uploaded datasets across
  sessions.
- `raptor/external_modules/acquisition/` — 10 source files:
  `__init__.py`, `base.py`, `datasets.py`, `cache.py`, `catalog.py`,
  `geo.py`, `tcga.py`, `arrayexpress.py`, `sra.py`,
  `gene_mapping.py`, `pooling.py`.
- `raptor/dashboard/pages/01___Data_Acquisition.py` — full dashboard
  page (~5200 lines) with 7 tabs.
- `tests/test_acquisition.py` — 105 tests including 17 mocked Entrez
  tests.
- `BETA_TESTING_GUIDE.md` — testing scenarios for beta testers.

### Changed

- `raptor/external_modules/__init__.py` — imports acquisition
  subpackage.
- `raptor/dashboard/app.py` — sidebar navigation includes Data
  Acquisition and Biomarker Discovery, version bumped to 2.2.2,
  workflow diagram updated.
- `raptor/launch_dashboard.py` — version reading made dynamic via
  `raptor.__version__` (no more hardcoded version strings); banner
  flexes for any version length; Module 6b acquisition check added.
- `setup.py` — `requests` and `pyarrow` added to core deps for
  Module 6b; `kneed>=0.8.0` added for Module 10 panel-size detection;
  `acquisition` extras group (GEOparse, biopython, mygene, combat).
- `requirements.txt` — acquisition + biomarker dependency sections
  added.
- `environment.yml` / `environment-full.yml` — acquisition + biomarker
  dependencies added.
- `check_raptor.py` — acquisition checks across all 15 diagnostic
  sections, smoke tests for offline functionality.
- `tests/test_cli_comprehensive.py` — version assertion now reads
  `raptor.__version__` dynamically (was hardcoded `2.2.0`).

### Fixed

- **TCGA multi-omic downloads** — miRNA, methylation, CNV, and RPPA
  data types now supported with expanded `gene_id_type` validation,
  consistent `_sample_label()`, data-type-aware labels, and
  `EOFError` catch for truncated `tar.gz` files.
- **Speculative version references** in `raptor/biomarker_discovery/core.py`
  — LOOCV deprecation docstring and DeprecationWarning string had
  been written referencing a non-existent `v2.2.3`; now correctly
  reference `v2.2.2` where the deprecation actually lands.
- **MANIFEST.in stale comment** — `tar.gz` example in usage notes
  referenced `raptor-2.2.1.tar.gz`; updated to `raptor-2.2.2.tar.gz`.

### Deprecated

- `validation='loocv'` in `discover_biomarkers()` — Ambroise &
  McLachlan (2002 PNAS) explicitly recommend k-fold over leave-one-out.
  LOOCV has higher variance with no bias advantage once the
  external-CV leakage fix is applied. Triggers a `DeprecationWarning`;
  use `validation='nested_cv'` with `n_folds=5` or `n_folds=10`.
  Removal scheduled for a future release.

### Known Issues

- TCGA and ArrayExpress connectors are partially implemented (GEO
  and SRA are production-ready).
- Gene ID conversion requires optional `mygene` package.
- ComBat batch correction requires optional `combat` package (falls
  back to median-centering).
- GEO download can be slow for very large datasets (>100 samples).
- `raptor dashboard` CLI command still broken (use
  `python -m raptor.launch_dashboard`).
- Prognostic intent on `raptor biomarker-survival` is stub-level
  (full implementation deferred).
- Dashboard pages Visualization, Reports, and Recommender still
  pending exhaustive testing.

### Test Suite

- **Module 6b**: 105 tests, 0 failures, 5 expected
  optional-dependency warnings. Diagnostic suite (`check_raptor.py`)
  scores 105/110.
- **Module 10**: 89 net-new tests across 11 test files. Full M10 CLI
  + dashboard suite passes 88/89 + 14/14 (kneedle CLI) + 21/21
  (kneedle algorithm + dashboard).

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

### Planned for v2.3.0
- Prognostic intent (full Cox regression-based survival biomarkers
  in `raptor biomarker-survival`)
- Batch confounding propagation in recommender
- GitHub Actions CI/CD (automated testing on Linux/macOS/Windows
  for every push)
- Sphinx documentation site hosted on ReadTheDocs
- Docker image published to Docker Hub
  (`ayehblk/raptor:<version>`)
- Conda-forge package distribution
- Real-data validation suite

---

**[2.2.2]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.2  
**[2.2.1]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.1  
**[2.2.0]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.2.0  
**[2.1.2]:** https://github.com/AyehBlk/RAPTOR/releases/tag/v2.1.2