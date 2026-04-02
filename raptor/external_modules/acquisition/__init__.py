"""
RAPTOR v2.2.2 - Module 6: Data Acquisition

Connect to public genomics repositories (GEO, TCGA, ArrayExpress, SRA),
download datasets, pool across studies, and feed into Module 10 for
biomarker discovery.

Part of Module 6 (External Connections) alongside module6_de_analysis.

Quick Start
-----------
>>> from raptor.external_modules.acquisition import GEOConnector
>>> geo = GEOConnector()
>>> results = geo.search("pancreatic cancer RNA-seq human")
>>> dataset = geo.download("GSE12345")
>>> print(dataset.summary())

>>> from raptor.external_modules.acquisition import AcquiredDataset
>>> my_data = AcquiredDataset.from_csv("my_counts.csv", "my_metadata.csv")

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

__version__ = "2.2.2"

# Data containers (always available)
from .datasets import (
    AcquiredDataset,
    PooledDataset,
    SUPPORTED_REPOSITORIES,
    SUPPORTED_DATA_TYPES,
    SUPPORTED_ORGANISMS,
)

# Cache manager (always available)
from .base import SearchResult

from .cache import (
    CacheManager,
    DEFAULT_CACHE_DIR,
)

# Connectors and engines (imported as available)
# These will be enabled as each file is built

_CATALOG_AVAILABLE = False
_GEO_AVAILABLE = False
_TCGA_AVAILABLE = False
_ARRAYEXPRESS_AVAILABLE = False
_SRA_AVAILABLE = False
_GENE_MAPPING_AVAILABLE = False
_POOLING_AVAILABLE = False

try:
    from .catalog import DataCatalog
    _CATALOG_AVAILABLE = True
except ImportError:
    DataCatalog = None

try:
    from .geo import GEOConnector
    _GEO_AVAILABLE = True
except ImportError:
    GEOConnector = None

try:
    from .tcga import TCGAConnector
    _TCGA_AVAILABLE = True
except ImportError:
    TCGAConnector = None

try:
    from .arrayexpress import ArrayExpConnector
    _ARRAYEXPRESS_AVAILABLE = True
except ImportError:
    ArrayExpConnector = None

try:
    from .sra import SRAConnector
    _SRA_AVAILABLE = True
except ImportError:
    SRAConnector = None

try:
    from .gene_mapping import GeneIDMapper
    _GENE_MAPPING_AVAILABLE = True
except ImportError:
    GeneIDMapper = None

try:
    from .pooling import PoolingEngine
    _POOLING_AVAILABLE = True
except ImportError:
    PoolingEngine = None


def get_available_components() -> dict:
    """Get availability status of all acquisition components."""
    return {
        'datasets': True,
        'cache': True,
        'catalog': _CATALOG_AVAILABLE,
        'geo_connector': _GEO_AVAILABLE,
        'tcga_connector': _TCGA_AVAILABLE,
        'arrayexpress_connector': _ARRAYEXPRESS_AVAILABLE,
        'sra_connector': _SRA_AVAILABLE,
        'gene_mapping': _GENE_MAPPING_AVAILABLE,
        'pooling_engine': _POOLING_AVAILABLE,
    }


__all__ = [
    # Version
    '__version__',

    # Data containers
    'AcquiredDataset',
    'PooledDataset',
    'SUPPORTED_REPOSITORIES',
    'SUPPORTED_DATA_TYPES',
    'SUPPORTED_ORGANISMS',

    # Cache
    'CacheManager',
    'DEFAULT_CACHE_DIR',

    # Catalog
    'DataCatalog',

    # Connectors
    'GEOConnector',
    'TCGAConnector',
    'ArrayExpConnector',
    'SRAConnector',

    # Utilities
    'GeneIDMapper',
    'PoolingEngine',
    'SearchResult',

    # Info
    'get_available_components',
]
