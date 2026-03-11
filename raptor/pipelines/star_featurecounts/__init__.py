"""
RAPTOR v2.2.0 - STAR + featureCounts Pipeline Package

This package provides STAR alignment + featureCounts for gene quantification.
Standard alignment-based workflow producing BAM files and gene counts.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

# =============================================================================
# PUBLIC API - v2.2.0
# =============================================================================

# Import main pipeline class
from .pipeline import (
    StarFeatureCountsPipeline,
    STAR_FC_CLI_PARAMS,
)

# =============================================================================
# MODULE METADATA
# =============================================================================

MODULE_INFO = {
    'id': 'M5',
    'name': 'star_featurecounts',
    'stage': 2,
    'description': 'STAR alignment + featureCounts (standard, produces BAM)',
    'version': '2.2.0',
    'requires': ['STAR>=2.7.0', 'subread>=2.0.0', 'samtools>=1.10'],
    'input': 'FASTQ files + sample sheet',
    'output': 'BAM files + gene_counts.csv',
    'typical_time': '4-12 hours',
    'memory_usage': '~32GB (high)',
}

PIPELINE_INFO = {
    'name': 'star_featurecounts',
    'class': StarFeatureCountsPipeline,
    'cli_params': STAR_FC_CLI_PARAMS,
    'description': 'STAR alignment + featureCounts quantification',
    'long_description': (
        'Standard alignment-based workflow using STAR for spliced alignment '
        'and featureCounts for gene-level quantification. Produces BAM files '
        'for downstream analysis (visualization, variant calling, etc.). '
        'Two-pass mode for better splice junction detection. '
        'High memory requirement (~32GB). Most accurate for gene quantification.'
    ),
    'features': [
        'Full genome alignment',
        'BAM files for downstream analysis',
        'Gene-level quantification',
        'Two-pass mode for sensitivity',
        'Splice junction detection',
        'High accuracy',
        'IGV visualization ready',
    ],
    'use_cases': [
        'Standard RNA-seq analysis',
        'When BAM files are needed',
        'Visualization in IGV/UCSC',
        'Variant calling from RNA-seq',
        'Novel isoform discovery',
        'Publication-quality analysis',
        'When accuracy is critical',
    ],
    'advantages': [
        'Most accurate alignment',
        'BAM files for visualization',
        'Gold standard for RNA-seq',
        'Well-established pipeline',
        'Widely used and trusted',
    ],
    'disadvantages': [
        'High memory usage (~32GB)',
        'Slower than pseudo-aligners',
        'Requires more disk space',
        'More complex setup',
    ],
    'docker_image': 'quay.io/biocontainers/star:2.7.10b',
    'version': '2.2.0',
}

# =============================================================================
# PUBLIC API EXPORTS
# =============================================================================

__all__ = [
    # Main pipeline class
    'StarFeatureCountsPipeline',
    
    # CLI parameters
    'STAR_FC_CLI_PARAMS',
    
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
