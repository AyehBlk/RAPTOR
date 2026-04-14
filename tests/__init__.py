"""
RAPTOR v2.2.0 - Test Suite

Comprehensive test suite covering all RAPTOR modules.

Run all tests:
    pytest tests/ -v

Run specific module tests:
    pytest tests/test_quality_assessment.py -v
    pytest tests/test_pipelines.py -v
    pytest tests/test_de_import.py -v
    pytest tests/test_parameter_optimization.py -v

Run with coverage:
    pytest tests/ --cov=raptor --cov-report=html

Run tests by marker:
    pytest tests/ -m "not slow"
    pytest tests/ -m "integration"

Test modules:
- test_detect_samples.py       : M1 - FASTQ detection
- test_combine_counts.py       : M1 - Count matrix combination
- test_generate_report.py      : M1 - QC report generation
- test_quant_pipelines.py      : M1 - Quick Salmon/Kallisto pipelines
- test_quality_assessment.py   : M2 - Quality assessment & outliers
- test_profiler.py             : M3 - Data profiling (32 features)
- test_recommender.py          : M4 - Pipeline recommendation (rule/ML)
- test_pipelines.py            : M5 - Production pipelines
- test_base_pipeline.py        : M5 - Base pipeline infrastructure
- test_pipelines_cli.py        : M5 - CLI integration tests
- test_de_import.py            : M7 - DE result import & standardization
- test_parameter_optimization.py : M8 - Parameter optimization
- test_acquisition.py           : M6b - Data acquisition & pooling

Module 5 (Production Pipelines) Tests:
======================================
The Module 5 test suite covers:

1. Base Infrastructure (test_base_pipeline.py):
   - SampleSheet parsing and validation
   - FASTQ auto-detection
   - Tool dependency checking
   - Base pipeline abstract class
   - Docker/HPC hybrid dependencies

2. Pipeline Classes (test_pipelines.py):
   - SalmonPipeline
   - KallistoPipeline
   - StarFeatureCountsPipeline
   - Hisat2FeatureCountsPipeline
   - StarRsemPipeline
   - StarSalmonPipeline
   - Pipeline registry functions

3. CLI Integration (test_pipelines_cli.py):
   - raptor pipeline list
   - raptor pipeline info
   - raptor pipeline salmon
   - raptor pipeline kallisto
   - raptor pipeline star-featurecounts
   - raptor pipeline hisat2-featurecounts
   - raptor pipeline star-rsem
   - raptor pipeline star-salmon
   - raptor pipeline run
   - raptor pipeline use-quantify

Module 7 (DE Import) Tests:
===========================
The Module 7 test suite covers:

1. Core Functionality (test_de_import.py):
   - DEResult class and all methods
   - Import from 4 DE methods (DESeq2, edgeR, limma, Wilcoxon)
   - Column mapping and standardization
   - Pipeline auto-detection
   - Threshold filtering
   - Serialization (pickle, JSON)
   - Performance metrics calculation
   - Integration with Module 6 outputs
   - CLI compatibility

Module 8 (Parameter Optimization) Tests:
========================================
The Module 8 test suite covers:

1. Core Functionality (test_parameter_optimization.py):
   - ParameterSpace and OptimizationResult classes
   - Ground truth optimization (with validated genes)
   - FDR control optimization (statistical method)
   - Search strategies (grid, random, differential evolution)
   - Parameter bounds validation
   - Convergence history tracking
   - Module 7 → Module 8 integration
   - CLI compatibility

Module 6b (Data Acquisition) Tests:
====================================
The Module 6b test suite covers (159 tests):

1. Data Containers (test_acquisition.py):
   - AcquiredDataset: creation, validation, integrity, filtering,
     subsetting, save/load (pkl + parquet), CSV export, from_csv
   - PooledDataset: creation, studies, LOSO splits, save/load

2. Cache & Catalog:
   - CacheManager: hybrid mode, save/load/delete, CSV export, clear
   - DataCatalog: register, search, filter, persistence

3. Connectors:
   - BaseConnector via MockConnector: search, download, auto-cache
   - GEO, TCGA, ArrayExpress, SRA constructors (offline)
   - ArrayExpConnector: _parse_organism, _parse_study_type,
     _simplify_platform, _detect_gene_ids, _parse_study_attrs,
     _parse_subsections, get_study_info, get_sample_types (28 tests)
   - SRAConnector: GSM extraction, GSE lookup, download_api (17 tests)
   - SearchResult formatting

4. Utilities:
   - GeneIDMapper: auto-detection, species aliases, validation
   - PoolingEngine: inner/outer join, 3 batch corrections,
     sample conflict resolution, min genes threshold

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# Test module imports for convenience
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Module information
MODULE_INFO = {
    'test_suite_version': '2.2.2',
    'modules_covered': ['M1', 'M2', 'M3', 'M4', 'M5', 'M6b', 'M7', 'M8'],
    'total_test_files': 13,
    'test_categories': [
        'unit',
        'integration', 
        'cli',
        'parametrized'
    ]
}

# Test file mapping
TEST_FILES = {
    # Module 1: Quick Count
    'M1': [
        'test_detect_samples.py',
        'test_combine_counts.py',
        'test_generate_report.py',
        'test_quant_pipelines.py'
    ],
    # Module 2: Quality Assessment
    'M2': [
        'test_quality_assessment.py'
    ],
    # Module 3: Data Profiling
    'M3': [
        'test_profiler.py'
    ],
    # Module 4: Pipeline Recommendation
    'M4': [
        'test_recommender.py'
    ],
    # Module 5: Production Pipelines
    'M5': [
        'test_pipelines.py',
        'test_base_pipeline.py',
        'test_pipelines_cli.py'
    ],
    # Module 6b: Data Acquisition
    'M6b': [
        'test_acquisition.py'
    ],
    # Module 7: DE Import
    'M7': [
        'test_de_import.py'
    ],
    # Module 8: Parameter Optimization
    'M8': [
        'test_parameter_optimization.py'
    ]
}


def get_test_files_for_module(module_id: str) -> list:
    """Get test files for a specific module.
    
    Parameters
    ----------
    module_id : str
        Module identifier (e.g., 'M1', 'M5', 'M7', 'M8')
    
    Returns
    -------
    list
        List of test file names for the module
    """
    return TEST_FILES.get(module_id, [])


def get_all_test_files() -> list:
    """Get all test files.
    
    Returns
    -------
    list
        List of all test file names
    """
    all_files = []
    for files in TEST_FILES.values():
        all_files.extend(files)
    return all_files


# Pytest markers
pytest_markers = """
# Custom markers for RAPTOR tests
# Add to pytest.ini or pyproject.toml:

[tool.pytest.ini_options]
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    "integration: marks tests as integration tests",
    "cli: marks tests as CLI integration tests",
    "m1: marks tests for Module 1 (Quick Count)",
    "m2: marks tests for Module 2 (Quality Assessment)",
    "m3: marks tests for Module 3 (Data Profiling)",
    "m4: marks tests for Module 4 (Pipeline Recommendation)",
    "m5: marks tests for Module 5 (Production Pipelines)",
    "m6b: marks tests for Module 6b (Data Acquisition)",
    "m7: marks tests for Module 7 (DE Import)",
    "m8: marks tests for Module 8 (Parameter Optimization)",
    "requires_docker: marks tests that require Docker",
    "requires_tools: marks tests that require external tools (salmon, kallisto, etc.)"
]
"""