"""
RAPTOR Quick Kallisto Pipeline - Module 1: Quantify

Ultra-fast pseudo-alignment with Kallisto for QC and profiling.
NOT for final differential expression analysis.

CRITICAL: For single-end reads, you MUST provide fragment_length and fragment_sd!
Kallisto CANNOT estimate fragment length from single-end data.

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
    "id": "quick_kallisto",
    "name": "Quick Kallisto Quantification",
    "version": __version__,
    "module": __module__,
    "stage": __stage__,
    "description": "Ultra-fast pseudo-alignment for QC and profiling",
    "output_dir": "results/quick_counts",
    "output_files": [
        "quick_gene_counts.csv",
        "quick_tpm.csv",
        "sample_info.csv"
    ],
    "notes": [
        "For single-end reads: MUST provide --fragment-length and --fragment-sd",
        "Kallisto CANNOT estimate fragment length from single-end data",
        "Use Quick Kallisto for QC only, not for final DE analysis"
    ]
}

# Import main functions with error handling - v2.2.0
# NOTE: No scripts/ subfolder - all files in quick_kallisto/ root
try:
    from .kallisto_quant import (
        QuickKallistoPipeline,
        run_quick_kallisto,
        QUICK_KALLISTO_CLI_PARAMS,
        PIPELINE_NAME,
        PIPELINE_VERSION,
        PIPELINE_MODULE,
        PIPELINE_STAGE
    )
    
    __all__ = [
        'QuickKallistoPipeline',
        'run_quick_kallisto',
        'QUICK_KALLISTO_CLI_PARAMS',
        'PIPELINE_INFO',
        'PIPELINE_NAME',
        'PIPELINE_VERSION',
        '__version__'
    ]
    
except ImportError as e:
    warnings.warn(
        f"Could not import Quick Kallisto components: {e}. "
        f"Some functionality may be limited."
    )
    
    __all__ = ['PIPELINE_INFO', '__version__']


# Convenience accessor functions - v2.2.0

def get_cli_params():
    """Get CLI parameter definitions for Quick Kallisto."""
    try:
        return QUICK_KALLISTO_CLI_PARAMS
    except NameError:
        return {}


def get_pipeline_info():
    """Get pipeline information dictionary."""
    return PIPELINE_INFO


def get_version():
    """Get pipeline version."""
    return __version__
