"""
RAPTOR v2.2.0 - STAR + RSEM Pipeline Package

This package provides STAR alignment + RSEM for isoform-level quantification.
Gold standard for accurate transcript and gene expression estimation.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# =============================================================================
# PUBLIC API - v2.2.0
# =============================================================================

# Import main pipeline class
from .pipeline import (
    StarRsemPipeline,
    STAR_RSEM_CLI_PARAMS,
)

# =============================================================================
# MODULE METADATA
# =============================================================================

MODULE_INFO = {
    'id': 'M6',
    'name': 'star_rsem',
    'stage': 2,
    'description': 'STAR + RSEM (gold standard, isoform-level)',
    'version': '2.2.0',
    'requires': ['STAR>=2.7.0', 'RSEM>=1.3.0', 'samtools>=1.10'],
    'input': 'FASTQ files + sample sheet',
    'output': 'Gene counts + Transcript counts + BAM files',
    'typical_time': '4-12 hours',
    'memory_usage': '~32GB (high)',
}

PIPELINE_INFO = {
    'name': 'star_rsem',
    'class': StarRsemPipeline,
    'cli_params': STAR_RSEM_CLI_PARAMS,
    'description': 'STAR + RSEM isoform quantification',
    'long_description': (
        'Gold standard for isoform-level quantification using STAR alignment '
        'and RSEM expectation-maximization. Produces both gene-level and '
        'transcript-level abundance estimates. Uses probabilistic allocation '
        'of multi-mapping reads for accurate isoform quantification. '
        'High memory requirement (~32GB). Most accurate for isoform analysis.'
    ),
    'features': [
        'Isoform-level quantification',
        'EM-based probabilistic counting',
        'Gene and transcript abundance',
        'STAR alignment (BAM files)',
        'Multi-mapping read handling',
        'Credibility intervals (optional)',
        'Most accurate isoform estimates',
    ],
    'use_cases': [
        'Isoform-level differential expression',
        'Alternative splicing analysis',
        'When accurate transcript counts needed',
        'Publication-quality isoform analysis',
        'Multi-mapping read resolution',
        'Transcript-level studies',
    ],
    'advantages': [
        'Gold standard for isoforms',
        'EM algorithm for multi-mappers',
        'Both gene and transcript counts',
        'BAM files for visualization',
        'Widely validated and trusted',
    ],
    'disadvantages': [
        'High memory usage (~32GB)',
        'Slower than pseudo-aligners',
        'More complex than featureCounts',
        'CI calculation very slow',
    ],
    'docker_image': 'quay.io/biocontainers/rsem:1.3.3',
    'version': '2.2.0',
}

# =============================================================================
# PUBLIC API EXPORTS
# =============================================================================

__all__ = [
    # Main pipeline class
    'StarRsemPipeline',
    
    # CLI parameters
    'STAR_RSEM_CLI_PARAMS',
    
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
