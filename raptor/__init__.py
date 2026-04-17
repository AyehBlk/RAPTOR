"""
RAPTOR v2.2.2 - RNA-seq Analysis Pipeline Testing and Optimization Resource

A comprehensive framework for RNA-seq quality control, analysis,
and pipeline optimization with ML-based recommendations.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
"""

__version__ = "2.2.2"
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
    MLPipelineRecommender,
    TrainingConfig,
    MLRecommendation,
    train_ml_recommender,
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
        import_de_results,
    )
    _DE_IMPORT_AVAILABLE = True
except ImportError:
    DEResult = None
    import_de_results = None
    _DE_IMPORT_AVAILABLE = False

# Module 8: Parameter Optimization
try:
    from .parameter_optimization import (
        ParameterSpace,
        OptimizationResult,
        ParameterOptimizer,
        GroundTruthOptimizer,
        FDRControlOptimizer,
        StabilityOptimizer,
        ReproducibilityOptimizer,
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
        EnsembleResult,
        ensemble_fisher,
        ensemble_brown,
        ensemble_rra,
        ensemble_voting,
        ensemble_weighted,
        ensemble_pvalue_combination,
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

# Module 10: Biomarker Discovery
try:
    from .biomarker_discovery import (
        BiomarkerResult,
        FeatureSelectionResult,
        ClassificationResult,
        PanelOptimizationResult,
        SurvivalResult,
        AnnotationResult,
        FeatureSelector,
        ClassifierEvaluator,
        PanelOptimizer,
        SurvivalAnalyzer,
        BiologicalAnnotator,
        discover_biomarkers,
        discover_survival_biomarkers,
        validate_biomarkers,
        get_dependencies_status,
        # M10 enhancement modules
        BiomarkerIntent,
        VALID_INTENTS,
        SignatureScore,
        build_signature_score,
        DirectionPattern,
        build_direction_pattern,
        ppv_npv_at_prevalence,
        bootstrap_ci,
        youdens_optimal_threshold,
        decision_curve_analysis,
        net_reclassification_improvement,
        RatioBiomarkerSearcher,
        RatioPair,
        RatioSearchResult,
        apply_ratios,
        build_ratio_features,
        EnhancedBiomarkerResult,
        enhance_biomarker_result,
    )
    _BIOMARKER_AVAILABLE = True
except ImportError:
    BiomarkerResult = None
    FeatureSelectionResult = None
    ClassificationResult = None
    PanelOptimizationResult = None
    SurvivalResult = None
    AnnotationResult = None
    FeatureSelector = None
    ClassifierEvaluator = None
    PanelOptimizer = None
    SurvivalAnalyzer = None
    BiologicalAnnotator = None
    discover_biomarkers = None
    discover_survival_biomarkers = None
    validate_biomarkers = None
    get_dependencies_status = None
    BiomarkerIntent = None
    VALID_INTENTS = None
    SignatureScore = None
    build_signature_score = None
    DirectionPattern = None
    build_direction_pattern = None
    ppv_npv_at_prevalence = None
    bootstrap_ci = None
    youdens_optimal_threshold = None
    decision_curve_analysis = None
    net_reclassification_improvement = None
    RatioBiomarkerSearcher = None
    RatioPair = None
    RatioSearchResult = None
    apply_ratios = None
    build_ratio_features = None
    EnhancedBiomarkerResult = None
    enhance_biomarker_result = None
    _BIOMARKER_AVAILABLE = False

# Module 6b: Data Acquisition
try:
    from .external_modules.acquisition import (
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
except ImportError:
    AcquiredDataset = None
    PooledDataset = None
    CacheManager = None
    DataCatalog = None
    GEOConnector = None
    TCGAConnector = None
    ArrayExpConnector = None
    SRAConnector = None
    GeneIDMapper = None
    PoolingEngine = None
    _ACQUISITION_AVAILABLE = False

# =============================================================================
# PUBLIC API (__all__)
# =============================================================================

__all__ = [
    # Version info
    '__version__',
    '__author__',
    '__email__',
    
    # UTILS - Validation
    'validate_count_matrix',
    'validate_metadata',
    'validate_group_column',
    'validate_file_path',
    'validate_directory_path',
    'validate_numeric_range',
    'validate_positive_integer',
    'validate_probability',
    'validate_de_result',
    'validate_ground_truth',
    'validate_fdr_target',
    'validate_bootstrap_iterations',
    'MIN_VALIDATED_GENES_ABSOLUTE',
    'MIN_VALIDATED_GENES_RECOMMENDED',
    'IDEAL_VALIDATED_GENES',
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
    
    # UTILS - Errors
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
    'handle_errors',
    'check_file_exists',
    'validate_output_writable',
    
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
    'quick_quality_check',
    
    # Module 4: Recommender
    'PipelineRecommender',
    'Recommendation',
    'recommend_pipeline',
    'MLPipelineRecommender',
    'TrainingConfig',
    'MLRecommendation',
    'train_ml_recommender',
    
    # Synthetic Benchmarks
    'SyntheticBenchmarkGenerator',
    'generate_training_data',
    
    # Module 5: Pipelines (optional)
    'get_pipeline',
    'list_pipelines',
    'BasePipeline',
    
    # Module 7: DE Import (optional)
    'DEResult',
    'import_de_results',
    
    # Module 8: Parameter Optimization (optional)
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
    
    # Module 9: Ensemble Analysis (optional)
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
    
    # Module 10: Biomarker Discovery (optional)
    'BiomarkerResult',
    'FeatureSelectionResult',
    'ClassificationResult',
    'PanelOptimizationResult',
    'SurvivalResult',
    'AnnotationResult',
    'FeatureSelector',
    'ClassifierEvaluator',
    'PanelOptimizer',
    'SurvivalAnalyzer',
    'BiologicalAnnotator',
    'discover_biomarkers',
    'discover_survival_biomarkers',
    'validate_biomarkers',
    'get_dependencies_status',
    'BiomarkerIntent',
    'VALID_INTENTS',
    'SignatureScore',
    'build_signature_score',
    'DirectionPattern',
    'build_direction_pattern',
    'ppv_npv_at_prevalence',
    'bootstrap_ci',
    'youdens_optimal_threshold',
    'decision_curve_analysis',
    'net_reclassification_improvement',
    'RatioBiomarkerSearcher',
    'RatioPair',
    'RatioSearchResult',
    'apply_ratios',
    'build_ratio_features',
    'EnhancedBiomarkerResult',
    'enhance_biomarker_result',
    'BiomarkerIntent',
    'VALID_INTENTS',
    'SignatureScore',
    'build_signature_score',
    'DirectionPattern',
    'build_direction_pattern',
    'ppv_npv_at_prevalence',
    'bootstrap_ci',
    'youdens_optimal_threshold',
    'decision_curve_analysis',
    'net_reclassification_improvement',
    'RatioBiomarkerSearcher',
    'RatioPair',
    'RatioSearchResult',
    'apply_ratios',
    'build_ratio_features',
    'EnhancedBiomarkerResult',
    'enhance_biomarker_result',
    'BiomarkerIntent',
    'VALID_INTENTS',
    'SignatureScore',
    'build_signature_score',
    'DirectionPattern',
    'build_direction_pattern',
    'ppv_npv_at_prevalence',
    'bootstrap_ci',
    'youdens_optimal_threshold',
    'decision_curve_analysis',
    'net_reclassification_improvement',
    'RatioBiomarkerSearcher',
    'RatioPair',
    'RatioSearchResult',
    'apply_ratios',
    'build_ratio_features',
    'EnhancedBiomarkerResult',
    'enhance_biomarker_result',
    
    # Module 6b: Data Acquisition (optional)
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

# =============================================================================
# PACKAGE INFORMATION FUNCTIONS
# =============================================================================

def get_version():
    """Get RAPTOR version string."""
    return __version__


def get_info():
    """Get package information."""
    return {
        'name': 'RAPTOR',
        'version': __version__,
        'author': __author__,
        'email': __email__,
        'description': 'RNA-seq Analysis Pipeline Testing and Optimization Resource',
        'url': 'https://github.com/AyehBlk/RAPTOR',
    }


def get_available_modules():
    """Get list of available modules."""
    return {
        'utils': _UTILS_AVAILABLE,
        'simulation': True,
        'profiler': True,
        'quality_assessment': True,
        'recommender': True,
        'ml_recommender': True,
        'synthetic_benchmarks': True,
        'pipelines': _PIPELINES_AVAILABLE,
        'de_import': _DE_IMPORT_AVAILABLE,
        'parameter_optimization': _PARAM_OPTIMIZATION_AVAILABLE,
        'ensemble': _ENSEMBLE_AVAILABLE,
        'biomarker_discovery': _BIOMARKER_AVAILABLE,
        'acquisition': _ACQUISITION_AVAILABLE,
    }


def validate_installation():
    """Comprehensive installation validation."""
    import sys
    
    report = {
        'version': __version__,
        'exports_valid': False,
        'modules_available': {},
        'missing_exports': [],
        'issues': []
    }
    
    current_module = sys.modules[__name__]
    missing_exports = [n for n in __all__ if not hasattr(current_module, n)]
    
    report['exports_valid'] = len(missing_exports) == 0
    report['missing_exports'] = missing_exports
    
    if missing_exports:
        report['issues'].append(
            f"{len(missing_exports)} exports missing: {', '.join(missing_exports[:5])}"
        )
    
    modules = get_available_modules()
    report['modules_available'] = modules
    
    unavailable = [k for k, v in modules.items() if not v]
    if unavailable:
        report['issues'].append(
            f"{len(unavailable)} modules unavailable: {', '.join(unavailable)}"
        )
    
    return report


# =============================================================================
# LOGGING SETUP
# =============================================================================

import logging
import os

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

if os.environ.get('RAPTOR_VALIDATE_EXPORTS', 'false').lower() == 'true':
    import sys as _sys
    _current = _sys.modules[__name__]
    _missing = [n for n in __all__ if not hasattr(_current, n)]
    if _missing:
        import warnings
        warnings.warn(f"__all__ has {len(_missing)} missing exports: {', '.join(_missing[:10])}")

logger.info(f"RAPTOR v{__version__} loaded successfully")

__all__.extend([
    'get_version',
    'get_info',
    'get_available_modules',
    'validate_installation',
])