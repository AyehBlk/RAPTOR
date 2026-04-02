"""
RAPTOR External Modules

Contains external analysis modules implemented in other languages (R, etc.)
and data acquisition tools for public genomics repositories.

Available Modules:
    - module6_de_analysis: Differential expression analysis (R-based)
    - acquisition: Data acquisition from GEO, TCGA, ArrayExpress, SRA
"""

__version__ = "2.2.2"

# Existing: R-based DE analysis
from . import module6_de_analysis

# New: Data acquisition subpackage
try:
    from . import acquisition
    from .acquisition import (
        AcquiredDataset,
        PooledDataset,
        CacheManager,
        DataCatalog,
        GEOConnector,
        TCGAConnector,
        ArrayExpConnector,
        SRAConnector,
        GeneIDMapper,
        PoolingEngine,
    )
    _ACQUISITION_AVAILABLE = True
except ImportError as e:
    import warnings
    warnings.warn(f"Acquisition subpackage not fully available: {e}")
    _ACQUISITION_AVAILABLE = False
    acquisition = None

__all__ = [
    'module6_de_analysis',
    'acquisition',
    'AcquiredDataset',
    'PooledDataset',
    'CacheManager',
    'DataCatalog',
    'GEOConnector',
    'TCGAConnector',
    'ArrayExpConnector',
    'SRAConnector',
    'GeneIDMapper',
    'PoolingEngine',
]
