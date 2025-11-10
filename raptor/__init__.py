"""
RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource

A comprehensive benchmarking framework for RNA-seq differential expression analysis
pipelines with intelligent, data-driven pipeline recommendations.

Author: Ayeh Bolouki
Affiliation: University of Namur, Belgium
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
License: MIT
"""

# Version information
__version__ = '2.0.0'
__author__ = 'Ayeh Bolouki'
__email__ = 'ayehbolouki1988@gmail.com'
__license__ = 'MIT'
__url__ = 'https://github.com/AyehBlk/RAPTOR'

# Package metadata
__all__ = [
    'RNAseqDataProfiler',
    'PipelineRecommender',
    'PipelineBenchmark',
    'DataSimulator',
    'ReportGenerator',
    '__version__',
]

# Import main classes for easy access
try:
    from raptor.profiler import RNAseqDataProfiler
    from raptor.recommender import PipelineRecommender
    from raptor.benchmark import PipelineBenchmark
    from raptor.simulate import DataSimulator
    from raptor.report import ReportGenerator
except ImportError as e:
    # Handle missing dependencies gracefully during installation
    import warnings
    warnings.warn(
        f"Some RAPTOR components could not be imported: {e}. "
        "This is normal during installation. If you see this after "
        "installation, please ensure all dependencies are installed.",
        ImportWarning
    )

# Package-level configuration
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Create package logger
logger = logging.getLogger(__name__)

# Welcome message (only shown once)
def _show_welcome():
    """Display welcome message on first import."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║                     🦖 RAPTOR v2.0.0                        ║
    ║   RNA-seq Analysis Pipeline Testing & Optimization Resource ║
    ║                                                              ║
    ║          Making pipeline selection evidence-based,           ║
    ║                      not guesswork.                          ║
    ║                                                              ║
    ║              Created by Ayeh Bolouki                         ║
    ║            University of Namur, Belgium                      ║
    ╚══════════════════════════════════════════════════════════════╝
    
    Quick Start:
    • raptor profile --counts data.csv    # Get recommendation
    • raptor compare --data fastq/        # Benchmark pipelines
    • raptor --help                       # See all commands
    
    Documentation: https://github.com/AyehBlk/RAPTOR
    """)

# Show welcome message only once per session
_WELCOME_SHOWN = False
if not _WELCOME_SHOWN:
    _show_welcome()
    _WELCOME_SHOWN = True
