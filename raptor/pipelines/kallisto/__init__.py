"""
RAPTOR v2.2.0 - Kallisto Pipeline Package

This package provides Kallisto pseudo-alignment for ultra-fast quantification.
No BAM files produced, ideal for quick differential expression analysis.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# =============================================================================
# PUBLIC API - v2.2.0
# =============================================================================

# Import main pipeline class
from .pipeline import (
    KallistoPipeline,
    KALLISTO_CLI_PARAMS,
)

# =============================================================================
# MODULE METADATA
# =============================================================================

MODULE_INFO = {
    'id': 'M5',
    'name': 'kallisto',
    'stage': 2,
    'description': 'Kallisto pseudo-alignment (ultra-fast, no BAM)',
    'version': '2.2.0',
    'requires': ['kallisto>=0.48.0'],
    'input': 'FASTQ files + sample sheet',
    'output': 'gene_counts.csv, TPM matrix',
    'typical_time': '15-60 minutes',
    'memory_usage': '~4-8GB (very low)',
}

PIPELINE_INFO = {
    'name': 'kallisto',
    'class': KallistoPipeline,
    'cli_params': KALLISTO_CLI_PARAMS,
    'description': 'Kallisto pseudo-alignment quantification',
    'long_description': (
        'Ultra-fast pseudo-alignment using Kallisto. Uses exact k-mer matching '
        'for transcript quantification without alignment. No BAM files produced. '
        'Lower memory usage than Salmon (~4GB vs 8GB), but no GC bias correction. '
        'Ideal for quick differential expression analysis.'
    ),
    'features': [
        'Ultra-fast quantification',
        'Very low memory usage (~4GB)',
        'Transcript-level output',
        'Optional bootstrap for uncertainty',
        'No BAM files (saves disk space)',
    ],
    'use_cases': [
        'Quick differential expression',
        'High-throughput projects (many samples)',
        'Limited RAM systems (<8GB)',
        'When speed is critical',
        'Transcript isoform analysis',
    ],
    'docker_image': 'zlskidmore/kallisto:0.50.0',
    'version': '2.2.0',
}

# =============================================================================
# PUBLIC API EXPORTS
# =============================================================================

__all__ = [
    # Main pipeline class
    'KallistoPipeline',
    
    # CLI parameters
    'KALLISTO_CLI_PARAMS',
    
    # Metadata
    'MODULE_INFO',
    'PIPELINE_INFO',
]

# =============================================================================
# VERSION INFO
# =============================================================================

__version__ = '2.2.0'
__author__ = 'Ayeh Bolouki'
__email__ = 'ayehbolouki1988@gmail.com'
