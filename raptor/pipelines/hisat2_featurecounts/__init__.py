"""
RAPTOR v2.2.0 - HISAT2 + featureCounts Pipeline Package

This package provides HISAT2 alignment + featureCounts quantification pipeline.
Lower memory alternative to STAR, ideal for systems with limited RAM.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# =============================================================================
# PUBLIC API - v2.2.0
# =============================================================================

# Import main pipeline class
from .pipeline import (
    Hisat2FeatureCountsPipeline,
    HISAT2_FC_CLI_PARAMS,
)

# =============================================================================
# MODULE METADATA
# =============================================================================

MODULE_INFO = {
    'id': 'M5',
    'name': 'hisat2_featurecounts',
    'stage': 2,
    'description': 'HISAT2 alignment + featureCounts (low memory)',
    'version': '2.2.0',
    'requires': ['hisat2>=2.2.0', 'samtools>=1.10', 'subread>=2.0.0'],
    'input': 'FASTQ files + sample sheet',
    'output': 'gene_counts.csv, BAM files',
    'typical_time': '4-8 hours',
    'memory_usage': '~16GB (lower than STAR)',
}

PIPELINE_INFO = {
    'name': 'hisat2_featurecounts',
    'class': Hisat2FeatureCountsPipeline,
    'cli_params': HISAT2_FC_CLI_PARAMS,
    'description': 'HISAT2 alignment + featureCounts quantification',
    'long_description': (
        'Production pipeline using HISAT2 for spliced alignment and '
        'featureCounts for gene quantification. Lower memory footprint '
        'than STAR (typically 16GB vs 32GB), making it ideal for systems '
        'with limited RAM or for processing many samples.'
    ),
    'features': [
        'Lower memory usage (~16GB)',
        'Fast spliced alignment',
        'Gene-level quantification',
        'BAM files for downstream analysis',
        'Compatible with limited-RAM systems',
    ],
    'use_cases': [
        'Systems with <32GB RAM',
        'Processing many samples in parallel',
        'When BAM files are needed',
        'Standard differential expression',
    ],
    'docker_image': 'quay.io/biocontainers/hisat2:2.2.1',
    'version': '2.2.0',
}

# =============================================================================
# PUBLIC API EXPORTS
# =============================================================================

__all__ = [
    # Main pipeline class
    'Hisat2FeatureCountsPipeline',
    
    # CLI parameters
    'HISAT2_FC_CLI_PARAMS',
    
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
