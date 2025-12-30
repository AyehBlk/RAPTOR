"""
RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource

A comprehensive benchmarking framework for RNA-seq differential expression analysis
pipelines with intelligent, data-driven pipeline recommendations powered by machine learning.

Version 2.1.2 is a hotfix release that fixes Python 3.8-3.11 compatibility
for the Adaptive Threshold Optimizer (ATO) feature introduced in v2.1.1.

Install from PyPI: pip install raptor-rnaseq

Author: Ayeh Bolouki
License: MIT
"""

# Version information
__version__ = '2.1.2'
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
    # ML classes (v2.1.0)
    'MLPipelineRecommender',
    'DataQualityAssessor',
    'ParameterOptimizer',
    'AutomatedReportGenerator',
    # Threshold Optimizer (v2.1.1) - NEW
    'AdaptiveThresholdOptimizer',
    'ThresholdResult',
    'optimize_thresholds',
    'THRESHOLD_OPTIMIZER_AVAILABLE',
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
        "pip install raptor-rnaseq[all]",
        ImportWarning
    )

# ==============================================================================
# Threshold Optimizer Module (v2.1.1) - NEW
# ==============================================================================

try:
    from raptor.threshold_optimizer import (
        # Main class
        AdaptiveThresholdOptimizer,
        
        # Result container
        ThresholdResult,
        
        # Convenience function
        optimize_thresholds,
        
        # Visualization functions
        plot_optimization_summary,
        plot_volcano,
        plot_logfc_distribution,
        plot_pvalue_distribution,
        plot_threshold_comparison,
        plot_adjustment_comparison
    )
    THRESHOLD_OPTIMIZER_AVAILABLE = True
except ImportError as e:
    import warnings
    warnings.warn(
        f"Threshold Optimizer module not available: {e}. "
        "Some features may be limited.",
        ImportWarning
    )
    THRESHOLD_OPTIMIZER_AVAILABLE = False
    
    # Create placeholder classes for graceful degradation
    AdaptiveThresholdOptimizer = None
    ThresholdResult = None
    optimize_thresholds = None

# ==============================================================================
# Package-level configuration
# ==============================================================================

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
    ╔══════════════════════════════════════════════════════════════════╗
    ║                     🦖 RAPTOR v2.1.2                             ║
    ║   RNA-seq Analysis Pipeline Testing & Optimization Resource      ║
    ║                                                                  ║
    ║          🤖 ML-POWERED RECOMMENDATIONS                           ║
    ║          📊 ADVANCED QUALITY ASSESSMENT                          ║
    ║          🎯 ADAPTIVE THRESHOLD OPTIMIZER (NEW!)                  ║
    ║          📄 AUTOMATED REPORTING                                  ║
    ║                                                                  ║
    ║              Created by Ayeh Bolouki                             ║
    ║             University of Liège, Belgium                         ║
    ║                                                                  ║
    ╚══════════════════════════════════════════════════════════════════╝
    
    What's New in v2.1.2:
    • 🎯 Adaptive Threshold Optimizer (ATO) for DE analysis
      - Data-driven logFC and padj threshold selection
      - Multiple p-value adjustment methods (BH, BY, q-value, Holm...)
      - π₀ estimation for true null proportion
      - Interactive dashboard integration
    
    Quick Start:
    • pip install raptor-rnaseq                        # Install from PyPI
    • raptor profile --counts data.csv --use-ml       # Get ML recommendation
    • python launch_dashboard.py                       # Launch web dashboard
    • raptor --help                                    # See all commands
    
    Threshold Optimizer:
    • from raptor import optimize_thresholds
    • result = optimize_thresholds(de_results, goal='discovery')
    • print(result.summary())
    
    PyPI: https://pypi.org/project/raptor-rnaseq/
    GitHub: https://github.com/AyehBlk/RAPTOR
    Making free science for everybody around the world 🌍
    """)

# Show welcome message
if not _WELCOME_SHOWN:
    try:
        _show_welcome()
        _WELCOME_SHOWN = True
    except:
        pass  # Suppress errors in non-interactive environments


# ==============================================================================
# Convenience functions
# ==============================================================================

def get_version():
    """Return the current RAPTOR version."""
    return __version__


def check_installation():
    """Check which RAPTOR components are available."""
    components = {
        'core': True,
        'ml_recommender': False,
        'threshold_optimizer': THRESHOLD_OPTIMIZER_AVAILABLE,
        'dashboard': False,
    }
    
    try:
        from raptor.ml_recommender import MLPipelineRecommender
        components['ml_recommender'] = True
    except ImportError:
        pass
    
    try:
        import streamlit
        components['dashboard'] = True
    except ImportError:
        pass
    
    return components


def print_status():
    """Print installation status."""
    components = check_installation()
    print("\n🦖 RAPTOR Installation Status:")
    print("=" * 40)
    for name, available in components.items():
        status = "✅" if available else "❌"
        print(f"  {status} {name}")
    print("=" * 40)
    print(f"  Version: {__version__}")
    print()
