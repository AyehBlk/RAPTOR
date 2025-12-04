"""
RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource

A comprehensive benchmarking framework for RNA-seq differential expression analysis
pipelines with intelligent, data-driven pipeline recommendations powered by machine learning.

Version 2.1.0 introduces ML-based recommendations, advanced quality assessment,
automated reporting, and parameter optimization capabilities.

Author: Ayeh Bolouki
Affiliation: University of Namur & GIGA-Neurosciences, University of Liège, Belgium
Email: ayeh.bolouki@unamur.be
License: MIT
"""

# Version information
__version__ = '2.1.0'
__author__ = 'Ayeh Bolouki'
__email__ = 'ayeh.bolouki@unamur.be'
__license__ = 'MIT'
__url__ = 'https://github.com/AyehBlk/RAPTOR'

# Package metadata
__all__ = [
    # Core classes (v2.0.0)
    'RNAseqDataProfiler',
    'PipelineRecommender',
    'PipelineBenchmark',
    'DataSimulator',
    'ReportGenerator',
    # New classes (v2.1.0)
    'MLPipelineRecommender',
    'DataQualityAssessor',
    'ParameterOptimizer',
    'AutomatedReportGenerator',
    # Version info
    '__version__',
]

# Import main classes for easy access
try:
    # Core v2.0.0 classes
    from raptor.profiler import RNAseqDataProfiler
    from raptor.recommender import PipelineRecommender
    from raptor.benchmark import PipelineBenchmark
    from raptor.simulate import DataSimulator
    from raptor.report import ReportGenerator
    
    # New v2.1.0 classes
    from raptor.ml_recommender import MLPipelineRecommender
    from raptor.data_quality_assessment import DataQualityAssessor
    from raptor.parameter_optimization import ParameterOptimizer
    from raptor.automated_reporting import AutomatedReportGenerator
    
except ImportError as e:
    # Handle missing dependencies gracefully during installation
    import warnings
    warnings.warn(
        f"Some RAPTOR components could not be imported: {e}. "
        "This is normal during installation. If you see this after "
        "installation, please ensure all dependencies are installed with: "
        "pip install -r requirements.txt",
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

# Welcome message (only shown once per session)
_WELCOME_SHOWN = False

def _show_welcome():
    """Display welcome message on first import."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║                     🦖 RAPTOR v2.1.0                        ║
    ║   RNA-seq Analysis Pipeline Testing & Optimization Resource ║
    ║                                                              ║
    ║          🤖 NOW WITH ML-POWERED RECOMMENDATIONS!            ║
    ║          📊 ADVANCED QUALITY ASSESSMENT                     ║
    ║          📄 AUTOMATED REPORTING                             ║
    ║                                                              ║
    ║              Created by Ayeh Bolouki                         ║
    ║         University of Namur & GIGA-Neurosciences            ║
    ║            University of Liège, Belgium                      ║
    ╚══════════════════════════════════════════════════════════════╝
    
    Quick Start:
    • python scripts/02_profile_data.py counts.csv    # Get ML recommendation
    • bash scripts/01_run_all_pipelines.sh data/      # Benchmark pipelines
    • python launch_dashboard.py                       # Launch web dashboard
    • raptor --help                                    # See all commands
    
    Documentation: https://github.com/AyehBlk/RAPTOR
    Making free science for everybody around the world 🌍
    """)

# Show welcome message
if not _WELCOME_SHOWN:
    try:
        _show_welcome()
        _WELCOME_SHOWN = True
    except:
        pass  # Suppress errors in non-interactive environments
