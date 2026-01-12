"""
RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource

A comprehensive framework for benchmarking RNA-seq pipelines and getting
intelligent, data-driven pipeline recommendations.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0

Example Usage:
    >>> from raptor import RNAseqDataProfiler, PipelineRecommender
    >>> profiler = RNAseqDataProfiler(counts_df)
    >>> profile = profiler.run_full_profile()
    >>> recommender = PipelineRecommender(profile)
    >>> recommendation = recommender.get_recommendation()

    >>> from raptor import SampleSheet
    >>> sheet = SampleSheet('samples.csv')
    >>> for sample in sheet:
    ...     print(sample.sample_id, sample.fastq_r1)
"""

__version__ = "2.2.0"
__author__ = "Ayeh Bolouki"
__email__ = "ayehbolouki1988@gmail.com"

# Core analysis modules
from .profiler import RNAseqDataProfiler
from .recommender import PipelineRecommender

# Sample sheet handling (for pipeline execution)
from .sample_sheet import (
    Sample,
    SampleSheet,
    auto_detect_samples,
    create_sample_sheet_template
)

# Optional: Import other modules as they're added
# from .benchmark import PipelineBenchmark
# from .simulate import DataSimulator
# from .report import ReportGenerator

__all__ = [
    # Version info
    '__version__',
    '__author__',
    '__email__',
    
    # Core analysis
    'RNAseqDataProfiler',
    'PipelineRecommender',
    
    # Sample sheet handling
    'Sample',
    'SampleSheet',
    'auto_detect_samples',
    'create_sample_sheet_template',
]
