"""
RAPTOR External Modules

Contains external analysis modules implemented in other languages (R, etc.)

Available Modules:
    - module6_de_analysis: Differential expression analysis (R-based)
"""

__version__ = "2.2.0"

# Import the module so it's actually accessible
from . import module6_de_analysis

__all__ = ['module6_de_analysis']
