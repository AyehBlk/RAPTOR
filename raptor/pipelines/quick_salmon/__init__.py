"""
RAPTOR Quick Salmon Pipeline - Module 1: Quantify

Fast quasi-mapping with Salmon for QC and profiling.
NOT for final differential expression analysis.

Salmon parameters optimized for speed (Module 1):
- Library type: auto-detect ('A')
- GC bias correction: enabled
- Validate mappings: enabled
- Bootstraps: 0 (speed over precision)

Output: results/quick_counts/quick_gene_counts.csv

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

from pathlib import Path
import warnings

__version__ = "2.2.0"
__module__ = "M1"
__stage__ = 1

# Pipeline metadata
PIPELINE_INFO = {
    "id": "quick_salmon",
    "name": "Quick Salmon Quantification",
    "version": __version__,
    "module": __module__,
    "stage": __stage__,
    "description": "Fast quasi-mapping with Salmon for QC and profiling",
    "output_dir": "results/quick_counts",
    "output_files": [
        "quick_gene_counts.csv",
        "quick_tpm.csv",
        "sample_info.csv"
    ],
    "notes": [
        "Uses Salmon auto-detect library type (A) for flexibility",
        "GC bias correction enabled for accuracy",
        "Bootstraps disabled (0) for maximum speed in Module 1",
        "Use Quick Salmon for QC only, not for final DE analysis"
    ]
}

# Import main functions with error handling - v2.2.0
# NOTE: No scripts/ subfolder - all files in quick_salmon/ root
try:
    from .salmon_quant import (
        QuickSalmonPipeline,
        run_quick_salmon,
        QUICK_SALMON_CLI_PARAMS,
        VALID_LIBRARY_TYPES,
        PIPELINE_NAME,
        PIPELINE_VERSION,
        PIPELINE_MODULE,
        PIPELINE_STAGE
    )
    
    __all__ = [
        'QuickSalmonPipeline',
        'run_quick_salmon',
        'QUICK_SALMON_CLI_PARAMS',
        'VALID_LIBRARY_TYPES',
        'PIPELINE_INFO',
        'PIPELINE_NAME',
        'PIPELINE_VERSION',
        '__version__'
    ]
    
except ImportError as e:
    warnings.warn(
        f"Could not import Quick Salmon components: {e}. "
        f"Some functionality may be limited."
    )
    
    __all__ = ['PIPELINE_INFO', '__version__']


# Convenience accessor functions - v2.2.0

def get_cli_params():
    """Get CLI parameter definitions for Quick Salmon."""
    try:
        return QUICK_SALMON_CLI_PARAMS
    except NameError:
        return {}


def get_pipeline_info():
    """Get pipeline information dictionary."""
    return PIPELINE_INFO


def get_version():
    """Get pipeline version."""
    return __version__


def get_pipeline_class():
    """Get the pipeline class."""
    try:
        return QuickSalmonPipeline
    except NameError:
        return None
