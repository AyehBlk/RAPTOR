"""
RAPTOR v2.2.0 - RNA-seq Analysis Pipeline Testing and Optimization Resource

A comprehensive framework for RNA-seq quality control, analysis,
and pipeline optimization with ML-based recommendations.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

__version__ = "2.2.0"
__author__ = "Ayeh Bolouki"
__email__ = "ayehbolouki1988@gmail.com"

# =============================================================================
# UTILS IMPORTS (Foundation - Always Load First)
# =============================================================================

try:
    # Validation functions
    from .utils.validation import (
        # Basic validation
        validate_count_matrix,
        validate_metadata,
        validate_group_column,
        validate_file_path,
        validate_directory_path,
        validate_numeric_range,
        validate_positive_integer,
        validate_probability,
        
        # Module 8: Parameter optimization validation
        validate_de_result,
        validate_ground_truth,
        validate_fdr_target,
        validate_bootstrap_iterations,
        MIN_VALIDATED_GENES_ABSOLUTE,
        MIN_VALIDATED_GENES_RECOMMENDED,
        IDEAL_VALIDATED_GENES,
        
        # Module 9: Ensemble analysis validation
        validate_de_results_dict,
        validate_ensemble_method,
        validate_significance_threshold,
        validate_direction_threshold,
        validate_weights_dict,
        validate_min_methods,
        validate_filters_dict,
        MIN_METHODS_ENSEMBLE,
        RECOMMENDED_METHODS_ENSEMBLE,
        IDEAL_METHODS_ENSEMBLE,
    )

    # Error classes and handlers
    from .utils.errors import (
        RAPTORError,
        ValidationError,
        DependencyError,
        PipelineError,
        
        # Module 8 errors
        OptimizationError,
        GroundTruthError,
        InsufficientDataError,
        
        # Module 9 errors
        EnsembleError,
        MethodMismatchError,
        DirectionInconsistencyError,
        CombinationFailedError,
        
        # Error handlers
        handle_errors,
        check_file_exists,
        validate_output_writable,
    )
    
    _UTILS_AVAILABLE = True

except ImportError as e:
    import warnings
    warnings.warn(f"Failed to import utils: {e}. Core functionality may be impaired.")
    _UTILS_AVAILABLE = False
    
    # Create dummy placeholders
    class ValidationError(Exception):
        pass
    class RAPTORError(Exception):
        pass

# =============================================================================
# CORE MODULE IMPORTS (Always Available)
# =============================================================================

# Simulation Module
from .simulation import (
    SimulationConfig,
    SimulationResult,
    RNAseqSimulator,
    simulate_rnaseq,
    simulate_diverse_datasets,
    simulate_small_sample_high_dispersion,
    simulate_large_sample_low_dispersion,
    simulate_medium_balanced,
    simulate_with_outliers,
)

# Module 3: Data Profiler
from .profiler import (
    DataProfile,
    RNAseqDataProfiler,
    profile_data_quick,
)

# Module 2: Quality Assessment
from .quality_assessment import (
    DataQualityAssessor,
    QualityReport,
    quick_quality_check,
)

# Module 4: Pipeline Recommender (Rule-based)
from .recommender import (
    PipelineRecommender,
    Recommendation,
    recommend_pipeline,
)

# Module 4: ML-Based Recommender
from .ml_recommender import (
    MLRecommender,
    TrainingConfig,
    MLRecommendation,
    train_recommender,
    load_recommender,
)

# Synthetic Benchmarks
from .synthetic_benchmarks import (
    SyntheticBenchmarkGenerator,
    generate_training_data,
)

# =============================================================================
# OPTIONAL MODULE IMPORTS (May Not Be Available)
# =============================================================================

# Module 5: Pipeline Base Classes
try:
    from .pipelines import (
        get_pipeline,
        list_pipelines,
        BasePipeline,
    )
    _PIPELINES_AVAILABLE = True
except ImportError:
    get_pipeline = None
    list_pipelines = None
    BasePipeline = None
    _PIPELINES_AVAILABLE = False

# Module 7: DE Import
try:
    from .de_import import (
        DEResult,
        DEImporter,
        import_deseq2,
        import_edger,
        import_limma,
        import_wilcoxon,
        import_de_result,
        compare_de_results,
        merge_de_results,
    )
    _DE_IMPORT_AVAILABLE = True
except ImportError:
    DEResult = None
    DEImporter = None
    import_deseq2 = None
    import_edger = None
    import_limma = None
    import_wilcoxon = None
    import_de_result = None
    compare_de_results = None
    merge_de_results = None
    _DE_IMPORT_AVAILABLE = False

# Module 8: Parameter Optimization
try:
    from .parameter_optimization import (
        # Core classes
        ParameterSpace,
        OptimizationResult,
        ParameterOptimizer,
        
        # Optimizer classes (4 methods)
        GroundTruthOptimizer,
        FDRControlOptimizer,
        StabilityOptimizer,
        ReproducibilityOptimizer,
        
        # Convenience functions (4 methods)
        optimize_with_ground_truth,
        optimize_with_fdr_control,
        optimize_with_stability,
        optimize_with_reproducibility,
    )
    _PARAM_OPTIMIZATION_AVAILABLE = True
except ImportError:
    ParameterSpace = None
    OptimizationResult = None
    ParameterOptimizer = None
    GroundTruthOptimizer = None
    FDRControlOptimizer = None
    StabilityOptimizer = None
    ReproducibilityOptimizer = None
    optimize_with_ground_truth = None
    optimize_with_fdr_control = None
    optimize_with_stability = None
    optimize_with_reproducibility = None
    _PARAM_OPTIMIZATION_AVAILABLE = False

# Module 9: Ensemble Analysis
try:
    from .ensemble import (
        # Result class
        EnsembleResult,
        
        # Main ensemble functions (5 methods)
        ensemble_fisher,
        ensemble_brown,
        ensemble_rra,
        ensemble_voting,
        ensemble_weighted,
        
        # Unified functions
        ensemble_pvalue_combination,
        
        # Lower-level functions
        fishers_method,
        browns_method,
        robust_rank_aggregation,
        check_direction_consistency,
        get_consensus_direction,
        calculate_direction_consistency_table,
        combine_pvalues_across_methods,
        calculate_meta_lfc,
    )
    _ENSEMBLE_AVAILABLE = True
except ImportError:
    EnsembleResult = None
    ensemble_fisher = None
    ensemble_brown = None
    ensemble_rra = None
    ensemble_voting = None
    ensemble_weighted = None
    ensemble_pvalue_combination = None
    fishers_method = None
    browns_method = None
    robust_rank_aggregation = None
    check_direction_consistency = None
    get_consensus_direction = None
    calculate_direction_consistency_table = None
    combine_pvalues_across_methods = None
    calculate_meta_lfc = None
    _ENSEMBLE_AVAILABLE = False

# =============================================================================
# PUBLIC API (__all__)
# =============================================================================

__all__ = [
    # Version info
    '__version__',
    '__author__',
    '__email__',
    
    # =============================================================================
    # UTILS (Validation & Error Handling)
    # =============================================================================
    
    # Basic validation functions
    'validate_count_matrix',
    'validate_metadata',
    'validate_group_column',
    'validate_file_path',
    'validate_directory_path',
    'validate_numeric_range',
    'validate_positive_integer',
    'validate_probability',
    
    # Module 8: Parameter optimization validation
    'validate_de_result',
    'validate_ground_truth',
    'validate_fdr_target',
    'validate_bootstrap_iterations',
    'MIN_VALIDATED_GENES_ABSOLUTE',
    'MIN_VALIDATED_GENES_RECOMMENDED',
    'IDEAL_VALIDATED_GENES',
    
    # Module 9: Ensemble analysis validation
    'validate_de_results_dict',
    'validate_ensemble_method',
    'validate_significance_threshold',
    'validate_direction_threshold',
    'validate_weights_dict',
    'validate_min_methods',
    'validate_filters_dict',
    'MIN_METHODS_ENSEMBLE',
    'RECOMMENDED_METHODS_ENSEMBLE',
    'IDEAL_METHODS_ENSEMBLE',
    
    # Error classes
    'RAPTORError',
    'ValidationError',
    'DependencyError',
    'PipelineError',
    'OptimizationError',
    'GroundTruthError',
    'InsufficientDataError',
    'EnsembleError',
    'MethodMismatchError',
    'DirectionInconsistencyError',
    'CombinationFailedError',
    
    # Error handlers
    'handle_errors',
    'check_file_exists',
    'validate_output_writable',
    
    # =============================================================================
    # CORE MODULES (Always Available)
    # =============================================================================
    
    # Simulation
    'SimulationConfig',
    'SimulationResult',
    'RNAseqSimulator',
    'simulate_rnaseq',
    'simulate_diverse_datasets',
    'simulate_small_sample_high_dispersion',
    'simulate_large_sample_low_dispersion',
    'simulate_medium_balanced',
    'simulate_with_outliers',
    
    # Module 3: Data Profiler
    'DataProfile',
    'RNAseqDataProfiler',
    'profile_data_quick',
    
    # Module 2: Quality Assessment
    'DataQualityAssessor',
    'QualityReport',
    'quick_quality_check',
    
    # Module 4: Recommender (Rule-based)
    'PipelineRecommender',
    'Recommendation',
    'recommend_pipeline',
    
    # Module 4: ML Recommender
    'MLRecommender',
    'TrainingConfig',
    'MLRecommendation',
    'train_recommender',
    'load_recommender',
    
    # Synthetic Benchmarks
    'SyntheticBenchmarkGenerator',
    'generate_training_data',
    
    # =============================================================================
    # OPTIONAL MODULES (May Not Be Available)
    # =============================================================================
    
    # Module 5: Pipelines
    'get_pipeline',
    'list_pipelines',
    'BasePipeline',
    
    # Module 7: DE Import
    'DEResult',
    'DEImporter',
    'import_deseq2',
    'import_edger',
    'import_limma',
    'import_wilcoxon',
    'import_de_result',
    'compare_de_results',
    'merge_de_results',
    
    # Module 8: Parameter Optimization
    'ParameterSpace',
    'OptimizationResult',
    'ParameterOptimizer',
    'GroundTruthOptimizer',
    'FDRControlOptimizer',
    'StabilityOptimizer',
    'ReproducibilityOptimizer',
    'optimize_with_ground_truth',
    'optimize_with_fdr_control',
    'optimize_with_stability',
    'optimize_with_reproducibility',
    
    # Module 9: Ensemble Analysis
    'EnsembleResult',
    'ensemble_fisher',
    'ensemble_brown',
    'ensemble_rra',
    'ensemble_voting',
    'ensemble_weighted',
    'ensemble_pvalue_combination',
    'fishers_method',
    'browns_method',
    'robust_rank_aggregation',
    'check_direction_consistency',
    'get_consensus_direction',
    'calculate_direction_consistency_table',
    'combine_pvalues_across_methods',
    'calculate_meta_lfc',
]

# =============================================================================
# PACKAGE INFORMATION FUNCTIONS
# =============================================================================

def get_version():
    """Get RAPTOR version string."""
    return __version__


def get_info():
    """
    Get package information.
    
    Returns
    -------
    dict
        Package metadata
    """
    return {
        'name': 'RAPTOR',
        'version': __version__,
        'author': __author__,
        'email': __email__,
        'description': 'RNA-seq Analysis Pipeline Testing and Optimization Resource',
        'url': 'https://github.com/AyehBlk/RAPTOR',
    }


def get_available_modules():
    """
    Get list of available modules.
    
    Returns
    -------
    dict
        Module availability status
    """
    modules = {
        'utils': _UTILS_AVAILABLE,
        'simulation': True,
        'profiler': True,  # Module 3
        'quality_assessment': True,  # Module 2
        'recommender': True,  # Module 4 (rule-based)
        'ml_recommender': True,  # Module 4 (ML)
        'synthetic_benchmarks': True,
        'pipelines': _PIPELINES_AVAILABLE,  # Module 5
        'de_import': _DE_IMPORT_AVAILABLE,  # Module 7
        'parameter_optimization': _PARAM_OPTIMIZATION_AVAILABLE,  # Module 8
        'ensemble': _ENSEMBLE_AVAILABLE,  # Module 9
    }
    return modules


def validate_installation():
    """
    Comprehensive installation validation.
    
    Checks:
    1. All __all__ exports exist
    2. All modules are importable
    3. Dependencies are available
    
    Returns
    -------
    dict
        Validation report with keys:
        - version: RAPTOR version
        - exports_valid: bool, whether all exports exist
        - modules_available: dict of module availability
        - missing_exports: list of missing items from __all__
        - issues: list of detected issues
    
    Examples
    --------
    >>> import raptor
    >>> report = raptor.validate_installation()
    >>> if report['exports_valid']:
    ...     print("✅ Installation valid")
    ... else:
    ...     print(f"❌ Missing exports: {report['missing_exports']}")
    """
    import sys
    
    report = {
        'version': __version__,
        'exports_valid': False,
        'modules_available': {},
        'missing_exports': [],
        'issues': []
    }
    
    # Check exports
    current_module = sys.modules[__name__]
    
    missing_exports = []
    for name in __all__:
        if not hasattr(current_module, name):
            missing_exports.append(name)
    
    report['exports_valid'] = len(missing_exports) == 0
    report['missing_exports'] = missing_exports
    
    if missing_exports:
        report['issues'].append(
            f"{len(missing_exports)} exports in __all__ don't exist: {', '.join(missing_exports[:5])}"
            f"{'...' if len(missing_exports) > 5 else ''}"
        )
    
    # Check modules
    modules = get_available_modules()
    report['modules_available'] = modules
    
    unavailable = [k for k, v in modules.items() if not v]
    if unavailable:
        report['issues'].append(
            f"{len(unavailable)} modules unavailable: {', '.join(unavailable)}"
        )
    
    return report


# =============================================================================
# INTERNAL VALIDATION (Development Mode)
# =============================================================================

def _validate_exports():
    """
    Validate that all items in __all__ actually exist.
    
    This catches typos, missing imports, or import failures.
    Run during development with: RAPTOR_VALIDATE_EXPORTS=true
    
    Returns
    -------
    bool
        True if all exports exist, False otherwise
    """
    import sys
    current_module = sys.modules[__name__]
    
    missing = []
    for name in __all__:
        if not hasattr(current_module, name):
            missing.append(name)
    
    if missing:
        import warnings
        warnings.warn(
            f"❌ __all__ contains {len(missing)} non-existent exports:\n"
            f"   {', '.join(missing[:10])}\n"
            f"   {'...' if len(missing) > 10 else ''}\n"
            f"   This indicates import errors or typos in __init__.py"
        )
        return False
    
    return True


# =============================================================================
# LOGGING SETUP
# =============================================================================

import logging
import os

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Run validation in development mode
if os.environ.get('RAPTOR_VALIDATE_EXPORTS', 'false').lower() == 'true':
    if _validate_exports():
        logger.info("✅ All exports validated successfully")
    else:
        logger.warning("⚠️  Export validation failed - see warnings above")

# Log version on import
logger.info(f"RAPTOR v{__version__} loaded successfully")

# Log available modules
modules = get_available_modules()
available = [k for k, v in modules.items() if v]
unavailable = [k for k, v in modules.items() if not v]

logger.debug(f"Available modules: {', '.join(available)}")
if unavailable:
    logger.debug(f"Unavailable modules: {', '.join(unavailable)}")

# Add helper functions to __all__
__all__.extend([
    'get_version',
    'get_info',
    'get_available_modules',
    'validate_installation',
])
