# Changelog

All notable changes to RAPTOR will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.2.0] - 2026-01-XX (Planned)

### Added
- **Module 9:** Complete ensemble analysis implementation (5 methods)
  - Fisher's method for p-value combination
  - Brown's method (correlation-aware)
  - Robust Rank Aggregation (RRA)
  - Voting-based consensus
  - Weighted ensemble with Module 8 optimization weights
- **Module 8:** Parameter optimization with 4 methods
  - Ground truth-based optimization
  - FDR control optimization
  - Stability-based optimization
  - Reproducibility-based optimization
- **Pipelines:** All 6 production pipelines fully validated
  - HISAT2 + featureCounts (16GB, low memory)
  - Kallisto (4GB, fastest)
  - Salmon (8GB, recommended)
  - STAR + featureCounts (32GB, gene-level gold standard)
  - STAR + RSEM (32GB, isoform-level gold standard)
  - STAR + Salmon (32GB, unique: BAM + bootstraps)
- **Dashboard:** Streamlit-based interactive interface
  - 9 pages covering Modules 2-4, 7-9
  - Real-time visualization
  - Settings and preferences
- **Documentation:** Comprehensive guides for all modules
- **CITATION.cff:** Academic citation support

### Changed
- Improved validation system with detailed error messages
- Enhanced CLI with better help messages
- Optimized memory usage in ensemble methods
- Better progress reporting with tqdm

### Fixed
- Dashboard folder name typo (`dashbord` → `dashboard`)
- Pipeline `__init__.py` filename issues
- Import path consistency across modules
- Memory leaks in large dataset processing

---

## [2.1.1] - 2025-12-30

### Added
- Initial PyPI release
- Basic module structure (Modules 2-4, 7)
- Core pipeline implementations
- Command-line interface
- Example scripts

### Fixed
- Import errors in CLI
- Missing dependencies in setup.py
- Documentation formatting

---

## [2.0.0] - 2025-11-XX

### Added
- Complete restructure to modular architecture
- Module 2: Quality Assessment
- Module 3: Data Profiler  
- Module 4: Pipeline Recommender (Rule-based + ML)
- Module 7: DE Import

---

## [1.0.0] - 2025-XX-XX

### Added
- Initial release
- Basic RNA-seq analysis functionality

---

## Version Number Scheme

RAPTOR uses semantic versioning: `MAJOR.MINOR.PATCH`

- **MAJOR:** Incompatible API changes
- **MINOR:** New features, backward compatible
- **PATCH:** Bug fixes, backward compatible

---

**Legend:**
- `Added` - New features
- `Changed` - Changes in existing functionality
- `Deprecated` - Soon-to-be removed features
- `Removed` - Removed features
- `Fixed` - Bug fixes
- `Security` - Security fixes
