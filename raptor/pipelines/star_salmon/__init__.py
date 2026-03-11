"""
RAPTOR v2.2.0 - STAR + Salmon Pipeline Package

This package provides STAR alignment + Salmon quantification.
Hybrid pipeline combining best of both worlds: BAM files + bootstrap uncertainty.

⚠️  REQUIRES TWO INDEXES:
    1. STAR genome index (--index)
    2. Salmon transcriptome index (--salmon-index)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# =============================================================================
# PUBLIC API - v2.2.0
# =============================================================================

# Import main pipeline class
from .pipeline import (
    StarSalmonPipeline,
    STAR_SALMON_CLI_PARAMS,
)

# =============================================================================
# MODULE METADATA
# =============================================================================

MODULE_INFO = {
    'id': 'M7',
    'name': 'star_salmon',
    'stage': 2,
    'description': 'STAR + Salmon (BAM + bootstraps)',
    'version': '2.2.0',
    'requires': ['STAR>=2.7.0', 'Salmon>=1.9.0', 'samtools>=1.10'],
    'input': 'FASTQ files + sample sheet',
    'output': 'Gene counts + TPM + BAM files',
    'typical_time': '4-12 hours',
    'memory_usage': '~32GB (high)',
    'special_requirements': 'TWO indexes (STAR genome + Salmon transcriptome)',
}

PIPELINE_INFO = {
    'name': 'star_salmon',
    'class': StarSalmonPipeline,
    'cli_params': STAR_SALMON_CLI_PARAMS,
    'description': 'STAR alignment + Salmon quantification (hybrid)',
    'long_description': (
        'Hybrid pipeline combining STAR alignment with Salmon quantification. '
        'STAR produces genome BAM and transcriptome BAM. Salmon quantifies from '
        'the transcriptome BAM with bootstrap uncertainty. Best of both worlds: '
        'you get BAM files for visualization AND bootstrap uncertainty for DE. '
        'Requires TWO indexes: STAR genome index and Salmon transcriptome index. '
        'High memory requirement (~32GB).'
    ),
    'features': [
        'BAM files from STAR',
        'Bootstrap uncertainty from Salmon',
        'Bias correction (GC + sequence)',
        'Two-pass mode for novel junctions',
        'Best of both worlds',
        'Useful for sleuth analysis',
    ],
    'use_cases': [
        'Need BAM files AND uncertainty estimates',
        'sleuth differential expression',
        'When you need both alignment and bootstraps',
        'Combining genome visualization with DE',
        'Publication with both IGV and uncertainty',
    ],
    'advantages': [
        'BAM files for visualization',
        'Bootstrap uncertainty for DE',
        'Combines STAR and Salmon benefits',
        'Bias correction from Salmon',
        'Widely compatible',
    ],
    'disadvantages': [
        'Requires TWO indexes (setup complexity)',
        'High memory usage (~32GB)',
        'Slower than pure Salmon',
        'More disk space (BAM files)',
        'Two-stage process',
    ],
    'docker_image': 'combinelab/salmon:1.10.0',
    'version': '2.2.0',
}

# =============================================================================
# PUBLIC API EXPORTS
# =============================================================================

__all__ = [
    # Main pipeline class
    'StarSalmonPipeline',
    
    # CLI parameters
    'STAR_SALMON_CLI_PARAMS',
    
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
