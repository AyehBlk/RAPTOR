"""
RAPTOR v2.2.0 - Salmon Pipeline Package

This package provides Salmon pseudo-alignment for ultra-fast quantification.
No BAM files produced, ideal for differential expression analysis.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# =============================================================================
# PUBLIC API - v2.2.0
# =============================================================================

# Import main pipeline class
from .pipeline import (
    SalmonPipeline,
    SALMON_CLI_PARAMS,
)

# =============================================================================
# MODULE METADATA
# =============================================================================

MODULE_INFO = {
    'id': 'M5',
    'name': 'salmon',
    'stage': 2,
    'description': 'Salmon pseudo-alignment (fastest, no BAM)',
    'version': '2.2.0',
    'requires': ['salmon>=1.9.0'],
    'input': 'FASTQ files + sample sheet',
    'output': 'gene_counts.csv, TPM matrix',
    'typical_time': '30-90 minutes',
    'memory_usage': '~8GB (moderate)',
}

PIPELINE_INFO = {
    'name': 'salmon',
    'class': SalmonPipeline,
    'cli_params': SALMON_CLI_PARAMS,
    'description': 'Salmon pseudo-alignment quantification',
    'long_description': (
        'Ultra-fast pseudo-alignment using Salmon. Uses quasi-mapping '
        'for transcript quantification without alignment. Includes GC bias '
        'and sequence bias correction. No BAM files produced. '
        'Moderate memory usage (~8GB). '
        'Ideal for fast differential expression analysis with bias correction.'
    ),
    'features': [
        'Ultra-fast quantification',
        'GC bias correction',
        'Sequence bias correction',
        'Transcript-level output',
        'Optional bootstrap for uncertainty',
        'No BAM files (saves disk space)',
        'Selective alignment (--validateMappings)',
    ],
    'use_cases': [
        'Fast differential expression',
        'High-throughput projects (many samples)',
        'When bias correction is important',
        'When speed and accuracy balance is needed',
        'Transcript isoform analysis',
        'Standard RNA-seq quantification',
    ],
    'advantages': [
        'Faster than Kallisto with bias correction',
        'More accurate than Kallisto (GC/seq bias)',
        'Lower memory than STAR/HISAT2',
        'Well-maintained and widely used',
    ],
    'docker_image': 'combinelab/salmon:1.10.0',
    'version': '2.2.0',
}

# =============================================================================
# PUBLIC API EXPORTS
# =============================================================================

__all__ = [
    # Main pipeline class
    'SalmonPipeline',
    
    # CLI parameters
    'SALMON_CLI_PARAMS',
    
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
