"""
RAPTOR Utils Package

Validation and error handling utilities for RAPTOR.
"""

from .validation import (
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
    
    # Module 1: Pipeline validation (v2.2.0 additions)
    validate_choice,
    check_directory_exists,
    validate_salmon_index,
    validate_kallisto_index,
    validate_fastq_file,
    validate_threads,
)

from .errors import (
    RAPTORError,
    ValidationError,
    DependencyError,
    PipelineError,
    
    # Module 8: Parameter optimization errors
    OptimizationError,
    GroundTruthError,
    InsufficientDataError,
    
    # Module 9: Ensemble analysis errors
    EnsembleError,
    MethodMismatchError,
    DirectionInconsistencyError,
    CombinationFailedError,
    
    # Error handling
    handle_errors,
    check_file_exists,
    validate_output_writable,
)

# Module 1: Sample sheet utilities (v2.2.0 additions)
from .sample_sheet import (
    Sample,
    SampleSheet,
    load_sample_sheet,
)

__all__ = [
    # Validation functions
    'validate_count_matrix',
    'validate_metadata',
    'validate_group_column',
    'validate_file_path',
    'validate_directory_path',
    'validate_numeric_range',
    'validate_positive_integer',
    'validate_probability',
    
    # Module 1: Pipeline validation (v2.2.0 additions)
    'validate_choice',
    'check_directory_exists',
    'validate_salmon_index',
    'validate_kallisto_index',
    'validate_fastq_file',
    'validate_threads',
    
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
    
    # Error classes and handlers
    'RAPTORError',
    'ValidationError',
    'DependencyError',
    'PipelineError',
    
    # Module 8: Parameter optimization errors
    'OptimizationError',
    'GroundTruthError',
    'InsufficientDataError',
    
    # Module 9: Ensemble analysis errors
    'EnsembleError',
    'MethodMismatchError',
    'DirectionInconsistencyError',
    'CombinationFailedError',
    
    # Error handling functions
    'handle_errors',
    'check_file_exists',
    'validate_output_writable',
    
    # Module 1: Sample sheet utilities (v2.2.0 additions)
    'Sample',
    'SampleSheet',
    'load_sample_sheet',
]
