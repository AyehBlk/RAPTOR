"""
RAPTOR v2.2.0 - Input Validation Utilities

Comprehensive validation functions for all RAPTOR inputs.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, List, Union, Tuple
import logging

logger = logging.getLogger(__name__)


# =============================================================================
# Count Matrix Validation
# =============================================================================

def validate_count_matrix(counts: pd.DataFrame, 
                          allow_negative: bool = False,
                          min_genes: int = 100,
                          min_samples: int = 2) -> bool:
    """
    Validate count matrix format and content.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes × samples)
    allow_negative : bool
        Allow negative values (for normalized data)
    min_genes : int
        Minimum number of genes required
    min_samples : int
        Minimum number of samples required
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If validation fails
    """
    # Check it's a DataFrame
    if not isinstance(counts, pd.DataFrame):
        raise TypeError(f"counts must be pd.DataFrame, got {type(counts)}")
    
    # Check not empty
    if counts.empty:
        raise ValueError("Count matrix is empty")
    
    # Check dimensions
    n_genes, n_samples = counts.shape
    if n_genes < min_genes:
        raise ValueError(
            f"Too few genes: {n_genes} < {min_genes}. "
            f"Is the matrix transposed?"
        )
    
    if n_samples < min_samples:
        raise ValueError(
            f"Too few samples: {n_samples} < {min_samples}"
        )
    
    # Check for negative values
    if not allow_negative and (counts < 0).any().any():
        raise ValueError(
            "Count matrix contains negative values. "
            "Raw counts must be non-negative."
        )
    
    # Check for unique gene names
    if not counts.index.is_unique:
        duplicates = counts.index[counts.index.duplicated()].unique()
        raise ValueError(
            f"Duplicate gene names found: {list(duplicates[:5])}"
            f"{' ...' if len(duplicates) > 5 else ''}"
        )
    
    # Check for unique sample names
    if not counts.columns.is_unique:
        duplicates = counts.columns[counts.columns.duplicated()].unique()
        raise ValueError(
            f"Duplicate sample names found: {list(duplicates[:5])}"
            f"{' ...' if len(duplicates) > 5 else ''}"
        )
    
    # Check for NaN values
    if counts.isna().any().any():
        n_nan = counts.isna().sum().sum()
        raise ValueError(
            f"Count matrix contains {n_nan} NaN values. "
            f"Please remove or impute missing values."
        )
    
    # Check for infinite values
    if np.isinf(counts.values).any():
        raise ValueError("Count matrix contains infinite values")
    
    # Warn if all zeros in a sample
    zero_samples = counts.columns[counts.sum(axis=0) == 0]
    if len(zero_samples) > 0:
        logger.warning(
            f"Samples with all zeros: {list(zero_samples)}. "
            f"Consider removing these samples."
        )
    
    # Warn if all zeros for a gene
    zero_genes = counts.index[counts.sum(axis=1) == 0]
    if len(zero_genes) > 100:
        logger.warning(
            f"{len(zero_genes)} genes have all zeros. "
            f"Consider filtering."
        )
    
    logger.info(
        f"✓ Count matrix valid: {n_genes:,} genes × {n_samples} samples"
    )
    
    return True


def validate_metadata(metadata: pd.DataFrame,
                      counts: Optional[pd.DataFrame] = None,
                      required_columns: Optional[List[str]] = None) -> bool:
    """
    Validate metadata format and content.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Sample metadata
    counts : pd.DataFrame, optional
        Count matrix to validate against
    required_columns : list, optional
        Required column names
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If validation fails
    """
    # Check it's a DataFrame
    if not isinstance(metadata, pd.DataFrame):
        raise TypeError(f"metadata must be pd.DataFrame, got {type(metadata)}")
    
    # Check not empty
    if metadata.empty:
        raise ValueError("Metadata is empty")
    
    # Check for sample_id column
    if 'sample_id' not in metadata.columns:
        raise ValueError(
            "Metadata must have 'sample_id' column. "
            f"Available columns: {list(metadata.columns)}"
        )
    
    # Check required columns
    if required_columns:
        missing = set(required_columns) - set(metadata.columns)
        if missing:
            raise ValueError(
                f"Missing required columns in metadata: {missing}"
            )
    
    # Check for duplicate sample IDs
    if metadata['sample_id'].duplicated().any():
        duplicates = metadata['sample_id'][metadata['sample_id'].duplicated()]
        raise ValueError(
            f"Duplicate sample IDs in metadata: {list(duplicates.unique())}"
        )
    
    # If counts provided, validate sample matching
    if counts is not None:
        count_samples = set(counts.columns)
        meta_samples = set(metadata['sample_id'])
        
        # Check all count samples have metadata
        missing_meta = count_samples - meta_samples
        if missing_meta:
            raise ValueError(
                f"Samples in counts missing from metadata: {missing_meta}"
            )
        
        # Warn if metadata has extra samples
        extra_meta = meta_samples - count_samples
        if extra_meta:
            logger.warning(
                f"Samples in metadata not in counts (will be ignored): "
                f"{extra_meta}"
            )
    
    logger.info(f"✓ Metadata valid: {len(metadata)} samples")
    
    return True


def validate_group_column(metadata: pd.DataFrame,
                         group_column: str,
                         min_groups: int = 2,
                         min_samples_per_group: int = 2) -> bool:
    """
    Validate experimental group column.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Sample metadata
    group_column : str
        Column name for groups/conditions
    min_groups : int
        Minimum number of groups
    min_samples_per_group : int
        Minimum samples per group
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If validation fails
    """
    # Check column exists
    if group_column not in metadata.columns:
        raise ValueError(
            f"Group column '{group_column}' not found in metadata. "
            f"Available columns: {list(metadata.columns)}"
        )
    
    # Get group sizes
    group_sizes = metadata[group_column].value_counts()
    
    # Check minimum number of groups
    if len(group_sizes) < min_groups:
        raise ValueError(
            f"Too few groups: {len(group_sizes)} < {min_groups}. "
            f"Need at least {min_groups} groups for comparison."
        )
    
    # Check minimum samples per group
    small_groups = group_sizes[group_sizes < min_samples_per_group]
    if len(small_groups) > 0:
        raise ValueError(
            f"Groups with <{min_samples_per_group} samples: "
            f"{dict(small_groups)}. "
            f"Each group needs ≥{min_samples_per_group} samples."
        )
    
    logger.info(
        f"✓ Group column '{group_column}' valid: "
        f"{len(group_sizes)} groups ({dict(group_sizes)})"
    )
    
    return True


# =============================================================================
# File Path Validation
# =============================================================================

def validate_file_path(filepath: Union[str, Path],
                       must_exist: bool = True,
                       file_type: Optional[str] = None) -> Path:
    """
    Validate file path.
    
    Parameters
    ----------
    filepath : str or Path
        File path to validate
    must_exist : bool
        If True, file must exist
    file_type : str, optional
        Expected file extension (e.g., '.csv', '.tsv')
    
    Returns
    -------
    Path
        Validated Path object
    
    Raises
    ------
    FileNotFoundError
        If file doesn't exist and must_exist=True
    ValueError
        If file type doesn't match
    """
    path = Path(filepath)
    
    if must_exist and not path.exists():
        raise FileNotFoundError(
            f"File not found: {filepath}"
        )
    
    if file_type and path.suffix.lower() != file_type.lower():
        raise ValueError(
            f"Expected {file_type} file, got {path.suffix}: {filepath}"
        )
    
    return path


def validate_directory_path(dirpath: Union[str, Path],
                           must_exist: bool = True,
                           create_if_missing: bool = False) -> Path:
    """
    Validate directory path.
    
    Parameters
    ----------
    dirpath : str or Path
        Directory path to validate
    must_exist : bool
        If True, directory must exist
    create_if_missing : bool
        If True, create directory if it doesn't exist
    
    Returns
    -------
    Path
        Validated Path object
    
    Raises
    ------
    NotADirectoryError
        If path exists but is not a directory
    FileNotFoundError
        If directory doesn't exist and must_exist=True
    """
    path = Path(dirpath)
    
    if path.exists() and not path.is_dir():
        raise NotADirectoryError(
            f"Path exists but is not a directory: {dirpath}"
        )
    
    if not path.exists():
        if must_exist and not create_if_missing:
            raise FileNotFoundError(
                f"Directory not found: {dirpath}"
            )
        elif create_if_missing:
            path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created directory: {dirpath}")
    
    return path


# =============================================================================
# Numerical Parameter Validation
# =============================================================================

def validate_numeric_range(value: Union[int, float],
                          name: str,
                          min_val: Optional[Union[int, float]] = None,
                          max_val: Optional[Union[int, float]] = None,
                          allow_none: bool = False) -> Union[int, float, None]:
    """
    Validate numerical parameter is within range.
    
    Parameters
    ----------
    value : int or float
        Value to validate
    name : str
        Parameter name (for error messages)
    min_val : int or float, optional
        Minimum allowed value (inclusive)
    max_val : int or float, optional
        Maximum allowed value (inclusive)
    allow_none : bool
        If True, None is acceptable
    
    Returns
    -------
    int, float, or None
        Validated value
    
    Raises
    ------
    ValueError
        If value is out of range
    """
    if value is None:
        if allow_none:
            return None
        else:
            raise ValueError(f"{name} cannot be None")
    
    if min_val is not None and value < min_val:
        raise ValueError(
            f"{name} must be >= {min_val}, got {value}"
        )
    
    if max_val is not None and value > max_val:
        raise ValueError(
            f"{name} must be <= {max_val}, got {value}"
        )
    
    return value


def validate_probability(value: float, name: str) -> float:
    """
    Validate probability value (0-1).
    
    Parameters
    ----------
    value : float
        Probability to validate
    name : str
        Parameter name
    
    Returns
    -------
    float
        Validated probability
    
    Raises
    ------
    ValueError
        If not in range [0, 1]
    """
    return validate_numeric_range(value, name, min_val=0.0, max_val=1.0)


def validate_positive_integer(value: int, name: str) -> int:
    """
    Validate positive integer.
    
    Parameters
    ----------
    value : int
        Integer to validate
    name : str
        Parameter name
    
    Returns
    -------
    int
        Validated integer
    
    Raises
    ------
    ValueError
        If not positive
    """
    return int(validate_numeric_range(value, name, min_val=1))


# =============================================================================
# Sample Sheet Validation
# =============================================================================

def validate_sample_sheet(sample_sheet: pd.DataFrame,
                         paired_end: bool = True) -> bool:
    """
    Validate sample sheet format.
    
    Parameters
    ----------
    sample_sheet : pd.DataFrame
        Sample sheet to validate
    paired_end : bool
        If True, expect both R1 and R2 columns
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If validation fails
    """
    required_cols = ['sample_id', 'fastq_r1']
    if paired_end:
        required_cols.append('fastq_r2')
    
    missing = set(required_cols) - set(sample_sheet.columns)
    if missing:
        raise ValueError(
            f"Sample sheet missing required columns: {missing}"
        )
    
    # Check for duplicate sample IDs
    if sample_sheet['sample_id'].duplicated().any():
        dups = sample_sheet['sample_id'][sample_sheet['sample_id'].duplicated()]
        raise ValueError(
            f"Duplicate sample IDs in sample sheet: {list(dups.unique())}"
        )
    
    # Check FASTQ files exist (warn, don't fail)
    for col in ['fastq_r1', 'fastq_r2']:
        if col in sample_sheet.columns:
            missing_files = []
            for fq in sample_sheet[col]:
                if pd.notna(fq) and fq != '' and not Path(fq).exists():
                    missing_files.append(fq)
            
            if missing_files:
                logger.warning(
                    f"{len(missing_files)} FASTQ files not found. "
                    f"First 3: {missing_files[:3]}"
                )
    
    logger.info(
        f"✓ Sample sheet valid: {len(sample_sheet)} samples"
    )
    
    return True


# =============================================================================
# Quick Validation Functions
# =============================================================================

def quick_validate_counts_and_metadata(
    counts: pd.DataFrame,
    metadata: Optional[pd.DataFrame] = None,
    group_column: str = 'condition'
) -> Tuple[bool, List[str]]:
    """
    Quick validation of counts and metadata.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix
    metadata : pd.DataFrame, optional
        Sample metadata
    group_column : str
        Group column name
    
    Returns
    -------
    tuple
        (is_valid, warnings_list)
    """
    warnings = []
    
    try:
        validate_count_matrix(counts)
    except Exception as e:
        return False, [f"Count matrix invalid: {e}"]
    
    if metadata is not None:
        try:
            validate_metadata(metadata, counts)
            validate_group_column(metadata, group_column)
        except Exception as e:
            return False, [f"Metadata invalid: {e}"]
    
    return True, warnings


# =============================================================================
# Module 8: Parameter Optimization Validation
# =============================================================================

# Validation constants for ground truth optimization
MIN_VALIDATED_GENES_ABSOLUTE = 10  # Absolute minimum
MIN_VALIDATED_GENES_RECOMMENDED = 20  # Recommended minimum
IDEAL_VALIDATED_GENES = 50  # Ideal for robust optimization


def validate_de_result(de_result: pd.DataFrame,
                       required_columns: Optional[List[str]] = None) -> bool:
    """
    Validate differential expression result format.
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results from DESeq2, edgeR, or limma
    required_columns : list, optional
        Required columns (default: ['log2FoldChange', 'pvalue'])
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If validation fails
    """
    # Check it's a DataFrame
    if not isinstance(de_result, pd.DataFrame):
        raise TypeError(f"de_result must be pd.DataFrame, got {type(de_result)}")
    
    # Check not empty
    if de_result.empty:
        raise ValueError("DE result is empty")
    
    # Default required columns
    if required_columns is None:
        required_columns = ['log2FoldChange', 'pvalue']
    
    # Check required columns
    missing = set(required_columns) - set(de_result.columns)
    if missing:
        raise ValueError(
            f"DE result missing required columns: {missing}. "
            f"Available columns: {list(de_result.columns)}"
        )
    
    # Check for gene identifiers
    if 'gene_id' not in de_result.columns and de_result.index.name is None:
        raise ValueError(
            "DE result must have 'gene_id' column or named index"
        )
    
    # Check for NaN values
    n_na_lfc = de_result['log2FoldChange'].isna().sum()
    n_na_pval = de_result['pvalue'].isna().sum()
    
    if n_na_lfc > 0:
        logger.warning(f"{n_na_lfc} genes have NaN log2FoldChange")
    if n_na_pval > 0:
        logger.warning(f"{n_na_pval} genes have NaN pvalue")
    
    # Check for infinite values
    if np.isinf(de_result['log2FoldChange']).any():
        logger.warning("DE result contains infinite log2FoldChange values")
    
    logger.info(f"✓ DE result valid: {len(de_result)} genes")
    
    return True


def validate_ground_truth(ground_truth: pd.DataFrame,
                         check_minimum: bool = True,
                         allow_borderline: bool = True) -> Tuple[bool, str]:
    """
    Validate ground truth gene list for parameter optimization.
    
    Parameters
    ----------
    ground_truth : pd.DataFrame
        Validated genes with 'gene_id' column
    check_minimum : bool
        If True, enforce minimum 10 genes
    allow_borderline : bool
        If True, allow 10-19 genes with warning
    
    Returns
    -------
    tuple
        (is_valid, message)
    
    Raises
    ------
    ValueError
        If validation fails
    """
    # Check it's a DataFrame
    if not isinstance(ground_truth, pd.DataFrame):
        raise TypeError(f"ground_truth must be pd.DataFrame, got {type(ground_truth)}")
    
    # Check not empty
    if ground_truth.empty:
        raise ValueError("Ground truth is empty")
    
    # Check for gene_id column
    if 'gene_id' not in ground_truth.columns:
        raise ValueError(
            "Ground truth must have 'gene_id' column. "
            f"Available columns: {list(ground_truth.columns)}"
        )
    
    n_validated = len(ground_truth)
    
    # Check minimum genes
    if check_minimum:
        if n_validated < MIN_VALIDATED_GENES_ABSOLUTE:
            raise ValueError(
                f"Only {n_validated} validated genes provided. "
                f"Minimum {MIN_VALIDATED_GENES_ABSOLUTE} required for optimization.\n"
                f"\n"
                f"RECOMMENDATION: Need at least {MIN_VALIDATED_GENES_RECOMMENDED} genes "
                f"for reliable optimization.\n"
                f"\n"
                f"Validation sources:\n"
                f"  - qRT-PCR confirmation (~$75/gene)\n"
                f"  - Western blot validation\n"
                f"  - Independent RNA-seq cohort\n"
                f"  - Literature-validated genes (same tissue/condition)\n"
                f"\n"
                f"Cost estimate for {MIN_VALIDATED_GENES_RECOMMENDED} genes: ~$2,700\n"
                f"\n"
                f"If validation not possible:\n"
                f"  - Skip Module 8 (Parameter Optimization)\n"
                f"  - Use standard thresholds (alpha=0.05, lfc=1.0)\n"
                f"  - Try alternative methods (FDR control, stability, reproducibility)"
            )
    
    # Provide feedback based on gene count
    if n_validated >= IDEAL_VALIDATED_GENES:
        message = f"✓ Excellent: {n_validated} validated genes (high-quality optimization expected)"
        logger.info(message)
        return True, message
    
    elif n_validated >= MIN_VALIDATED_GENES_RECOMMENDED:
        message = f"✓ Good: {n_validated} validated genes (recommend ≥{IDEAL_VALIDATED_GENES} for ideal optimization)"
        logger.info(message)
        return True, message
    
    elif n_validated >= MIN_VALIDATED_GENES_ABSOLUTE:
        message = (
            f"⚠️  WARNING: Only {n_validated} validated genes\n"
            f"This is below the recommended minimum of {MIN_VALIDATED_GENES_RECOMMENDED}.\n"
            f"Optimization may be unreliable.\n"
            f"\n"
            f"Recommendations:\n"
            f"  BEST: Validate {MIN_VALIDATED_GENES_RECOMMENDED - n_validated} more genes "
            f"(cost: ~${(MIN_VALIDATED_GENES_RECOMMENDED - n_validated) * 75})\n"
            f"  ALTERNATIVE: Use optimization results cautiously\n"
            f"  ALTERNATIVE: Try FDR control or stability methods instead"
        )
        
        if allow_borderline:
            logger.warning(message)
            return True, message
        else:
            raise ValueError(message)
    
    # Should not reach here if check_minimum=True
    return True, f"✓ {n_validated} validated genes"


def validate_fdr_target(target_fdr: float) -> bool:
    """
    Validate target FDR value for FDR control optimization.
    
    Parameters
    ----------
    target_fdr : float
        Target false discovery rate
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If not in valid range
    """
    if not 0 < target_fdr < 1:
        raise ValueError(
            f"target_fdr must be between 0 and 1, got {target_fdr}. "
            f"Common values: 0.01 (1%), 0.05 (5%), 0.10 (10%)"
        )
    
    if target_fdr < 0.001:
        logger.warning(
            f"target_fdr={target_fdr} is very stringent. "
            f"May result in very few significant genes."
        )
    
    if target_fdr > 0.20:
        logger.warning(
            f"target_fdr={target_fdr} is very relaxed. "
            f"May result in many false positives."
        )
    
    return True


def validate_bootstrap_iterations(n_bootstrap: int,
                                  minimum: int = 50,
                                  recommended: int = 100) -> bool:
    """
    Validate number of bootstrap iterations for stability optimization.
    
    Parameters
    ----------
    n_bootstrap : int
        Number of bootstrap iterations
    minimum : int
        Minimum acceptable
    recommended : int
        Recommended value
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValueError
        If below minimum
    """
    if n_bootstrap < minimum:
        raise ValueError(
            f"n_bootstrap must be ≥{minimum}, got {n_bootstrap}. "
            f"Recommend {recommended} for reliable results."
        )
    
    if n_bootstrap < recommended:
        logger.warning(
            f"n_bootstrap={n_bootstrap} is below recommended {recommended}. "
            f"Consider increasing for more reliable results."
        )
    
    if n_bootstrap > 500:
        logger.warning(
            f"n_bootstrap={n_bootstrap} is very high. "
            f"Optimization will be slow. Consider using 100-200."
        )
    
    return True


if __name__ == '__main__':
    print("RAPTOR v2.2.0 - Validation Utilities")
    print("=" * 60)
    print("\nThis module provides validation functions for RAPTOR inputs.")
    print("\nImport and use as:")
    print("  from raptor.utils.validation import validate_count_matrix")
    print("  validate_count_matrix(counts_df)")
# =============================================================================
# Module 9: Ensemble Analysis Validation
# =============================================================================

def validate_de_results_dict(de_results: Dict[str, Any],
                             min_methods: int = 1) -> bool:
    """
    Validate dictionary of DEResult objects for ensemble analysis.
    
    Parameters
    ----------
    de_results : dict
        Dictionary mapping method names to DEResult objects
    min_methods : int
        Minimum number of methods required
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If validation fails
    """
    from .errors import ValidationError, InsufficientDataError
    
    # Check not empty
    if not de_results:
        raise ValidationError(
            "de_results cannot be empty. "
            "Provide at least one DE method result."
        )
    
    # Check minimum methods
    n_methods = len(de_results)
    if n_methods < min_methods:
        raise InsufficientDataError(
            f"Insufficient methods: {n_methods} < {min_methods}. "
            f"Ensemble analysis requires at least {min_methods} method(s)."
        )
    
    # Check all are dictionaries with data attribute or DataFrames
    for method_name, result in de_results.items():
        if not isinstance(method_name, str):
            raise ValidationError(
                f"Method names must be strings, got {type(method_name)} for {method_name}"
            )
        
        # Check has 'data' attribute (DEResult) or is DataFrame
        if hasattr(result, 'data'):
            data = result.data
        elif isinstance(result, pd.DataFrame):
            data = result
        else:
            raise ValidationError(
                f"Method '{method_name}': result must have 'data' attribute or be DataFrame, "
                f"got {type(result)}"
            )
        
        # Check required columns
        required_cols = ['gene_id', 'log2FoldChange', 'pvalue']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            raise ValidationError(
                f"Method '{method_name}': missing required columns: {missing_cols}. "
                f"Required: {required_cols}"
            )
        
        # Check not empty
        if len(data) == 0:
            raise InsufficientDataError(
                f"Method '{method_name}': DE results are empty. "
                f"At least some genes required for ensemble analysis."
            )
        
        # Check for padj if using adjusted p-values
        if 'padj' not in data.columns:
            logger.warning(
                f"Method '{method_name}': 'padj' column not found. "
                f"Will not be able to use adjusted p-values."
            )
    
    logger.info(f"✅ Validated {n_methods} DE method results")
    return True


def validate_ensemble_method(method: str,
                             allowed_methods: Optional[List[str]] = None) -> bool:
    """
    Validate ensemble method name.
    
    Parameters
    ----------
    method : str
        Ensemble method name
    allowed_methods : list, optional
        Allowed method names
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If method invalid
    """
    from .errors import ValidationError
    
    if allowed_methods is None:
        allowed_methods = ['fisher', 'brown', 'rra', 'voting', 'weighted']
    
    if not isinstance(method, str):
        raise ValidationError(
            f"method must be string, got {type(method)}"
        )
    
    if method not in allowed_methods:
        raise ValidationError(
            f"Unknown ensemble method: '{method}'. "
            f"Allowed methods: {allowed_methods}"
        )
    
    return True


def validate_significance_threshold(threshold: float,
                                    min_value: float = 0.0,
                                    max_value: float = 1.0) -> bool:
    """
    Validate significance threshold (p-value or FDR).
    
    Parameters
    ----------
    threshold : float
        Significance threshold
    min_value : float
        Minimum allowed value (exclusive)
    max_value : float
        Maximum allowed value (exclusive)
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If threshold invalid
    """
    from .errors import ValidationError
    
    if not isinstance(threshold, (int, float)):
        raise ValidationError(
            f"threshold must be numeric, got {type(threshold)}"
        )
    
    if not (min_value < threshold < max_value):
        raise ValidationError(
            f"threshold must be in ({min_value}, {max_value}), got {threshold}"
        )
    
    # Warn if unusual
    if threshold > 0.1:
        logger.warning(
            f"⚠️  Significance threshold is high: {threshold}. "
            f"Common values: 0.01, 0.05"
        )
    
    return True


def validate_direction_threshold(threshold: float) -> bool:
    """
    Validate direction consistency threshold.
    
    Parameters
    ----------
    threshold : float
        Direction agreement threshold (0-1)
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If threshold invalid
    """
    from .errors import ValidationError
    
    if not isinstance(threshold, (int, float)):
        raise ValidationError(
            f"direction_threshold must be numeric, got {type(threshold)}"
        )
    
    if not (0 < threshold <= 1):
        raise ValidationError(
            f"direction_threshold must be in (0, 1], got {threshold}"
        )
    
    # Inform about threshold meaning
    if threshold == 1.0:
        logger.info("ℹ️  Direction threshold = 1.0: all methods must agree (strict)")
    elif threshold >= 0.67:
        logger.info(f"ℹ️  Direction threshold = {threshold}: 2/3 majority required")
    elif threshold >= 0.5:
        logger.info(f"ℹ️  Direction threshold = {threshold}: simple majority required")
    else:
        logger.warning(
            f"⚠️  Direction threshold = {threshold}: very permissive! "
            f"Genes with inconsistent direction will be included."
        )
    
    return True


def validate_weights_dict(weights: Dict[str, float],
                         method_names: List[str],
                         require_all: bool = False) -> bool:
    """
    Validate weights dictionary for weighted ensemble.
    
    Parameters
    ----------
    weights : dict
        Dictionary mapping method names to weights
    method_names : list
        List of all method names
    require_all : bool
        If True, require weights for all methods
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If weights invalid
    """
    from .errors import ValidationError
    
    if not isinstance(weights, dict):
        raise ValidationError(
            f"weights must be dictionary, got {type(weights)}"
        )
    
    # Check weights for known methods only
    unknown_methods = set(weights.keys()) - set(method_names)
    if unknown_methods:
        raise ValidationError(
            f"Unknown methods in weights: {unknown_methods}. "
            f"Available methods: {method_names}"
        )
    
    # Check all methods have weights if required
    if require_all:
        missing_methods = set(method_names) - set(weights.keys())
        if missing_methods:
            raise ValidationError(
                f"Missing weights for methods: {missing_methods}. "
                f"All methods must have weights: {method_names}"
            )
    
    # Check all weights are positive
    for method, weight in weights.items():
        if not isinstance(weight, (int, float)):
            raise ValidationError(
                f"Weight for '{method}' must be numeric, got {type(weight)}"
            )
        
        if weight < 0:
            raise ValidationError(
                f"Weight for '{method}' must be non-negative, got {weight}"
            )
        
        if weight == 0:
            logger.warning(
                f"⚠️  Weight for '{method}' is zero - this method will be ignored"
            )
    
    # Check sum of weights
    total_weight = sum(weights.values())
    if total_weight == 0:
        raise ValidationError(
            "Sum of weights is zero. At least one weight must be positive."
        )
    
    # Normalize if needed
    if abs(total_weight - 1.0) > 0.01:
        logger.info(
            f"ℹ️  Weights sum to {total_weight:.3f}, not 1.0. "
            f"Will be normalized automatically."
        )
    
    return True


def validate_min_methods(min_methods: int,
                        n_methods: int) -> bool:
    """
    Validate minimum methods parameter for voting.
    
    Parameters
    ----------
    min_methods : int
        Minimum number of methods required to detect gene
    n_methods : int
        Total number of methods available
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If min_methods invalid
    """
    from .errors import ValidationError
    
    if not isinstance(min_methods, int):
        raise ValidationError(
            f"min_methods must be integer, got {type(min_methods)}"
        )
    
    if min_methods < 1:
        raise ValidationError(
            f"min_methods must be ≥ 1, got {min_methods}"
        )
    
    if min_methods > n_methods:
        raise ValidationError(
            f"min_methods ({min_methods}) cannot exceed total methods ({n_methods})"
        )
    
    # Inform about stringency
    if min_methods == 1:
        logger.warning(
            "⚠️  min_methods = 1: very permissive, genes detected by ANY method included"
        )
    elif min_methods == n_methods:
        logger.info(
            f"ℹ️  min_methods = {n_methods}: very stringent, only unanimous genes included"
        )
    else:
        fraction = min_methods / n_methods
        if fraction >= 0.67:
            logger.info(
                f"ℹ️  min_methods = {min_methods}/{n_methods}: stringent, majority required"
            )
        else:
            logger.info(
                f"ℹ️  min_methods = {min_methods}/{n_methods}: moderate stringency"
            )
    
    return True


def validate_filters_dict(filters: Dict[str, float]) -> bool:
    """
    Validate filters dictionary for voting/weighted ensemble.
    
    Parameters
    ----------
    filters : dict
        Dictionary of filters (e.g., {'padj': 0.05, 'lfc': 1.0})
    
    Returns
    -------
    bool
        True if valid
    
    Raises
    ------
    ValidationError
        If filters invalid
    """
    from .errors import ValidationError
    
    if not isinstance(filters, dict):
        raise ValidationError(
            f"filters must be dictionary, got {type(filters)}"
        )
    
    allowed_filters = ['pvalue', 'padj', 'lfc', 'log2FoldChange', 'baseMean']
    
    for filter_name, filter_value in filters.items():
        # Check filter name
        if filter_name not in allowed_filters:
            logger.warning(
                f"⚠️  Unknown filter: '{filter_name}'. "
                f"Common filters: {allowed_filters}"
            )
        
        # Check filter value
        if not isinstance(filter_value, (int, float)):
            raise ValidationError(
                f"Filter value for '{filter_name}' must be numeric, "
                f"got {type(filter_value)}"
            )
        
        # Check p-value filters
        if filter_name in ['pvalue', 'padj']:
            if not (0 < filter_value < 1):
                raise ValidationError(
                    f"Filter '{filter_name}' must be in (0, 1), got {filter_value}"
                )
        
        # Check LFC filters
        if filter_name in ['lfc', 'log2FoldChange']:
            if filter_value < 0:
                logger.warning(
                    f"⚠️  Negative LFC filter: {filter_value}. "
                    f"Will filter on absolute LFC."
                )
    
    return True


# Constants for ensemble analysis
MIN_METHODS_ENSEMBLE = 2  # Minimum methods for meaningful ensemble
RECOMMENDED_METHODS_ENSEMBLE = 3  # Recommended for robust results
IDEAL_METHODS_ENSEMBLE = 4  # Ideal for comprehensive analysis
# =============================================================================
# Additional Validation Functions for Module 1 (Quick Salmon Pipeline)
# Add these to the END of validation.py (after line 1269, before constants)
# =============================================================================

"""
These functions are required by Module 1 (Quick Salmon) and other pipelines.
They provide general-purpose validation for common bioinformatics inputs.
"""

import os


def validate_choice(
    value: Any,
    choices: List[Any],
    param_name: str
) -> Any:
    """
    Validate that value is one of allowed choices.
    
    Commonly used for:
    - Library types (Salmon/Kallisto)
    - Strand orientation
    - Algorithm options
    
    Parameters
    ----------
    value : Any
        Value to validate
    choices : List[Any]
        List of valid choices
    param_name : str
        Parameter name for error message
    
    Returns
    -------
    Any
        Validated value
    
    Raises
    ------
    ValidationError
        If value not in choices
    
    Examples
    --------
    >>> validate_choice('A', ['A', 'IU', 'ISR'], 'lib_type')
    'A'
    
    >>> validate_choice('XYZ', ['A', 'IU', 'ISR'], 'lib_type')
    ValidationError: Invalid lib_type: XYZ. Must be one of: A, IU, ISR
    """
    from .errors import ValidationError
    
    if value not in choices:
        raise ValidationError(
            f"Invalid {param_name}: {value}\n"
            f"Must be one of: {', '.join(map(str, choices))}"
        )
    return value


def check_file_exists(
    filepath: Union[str, Path],
    param_name: str = "file"
) -> Path:
    """
    Check that file exists and is readable.
    
    Parameters
    ----------
    filepath : Union[str, Path]
        File path to check
    param_name : str
        Parameter name for error message
    
    Returns
    -------
    Path
        Validated file path
    
    Raises
    ------
    ValidationError
        If file doesn't exist or is not a file
    
    Examples
    --------
    >>> check_file_exists('samples.csv')
    Path('samples.csv')
    
    >>> check_file_exists('missing.csv')
    ValidationError: file not found: missing.csv
    """
    from .errors import ValidationError
    
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise ValidationError(
            f"{param_name} not found: {filepath}"
        )
    
    if not filepath.is_file():
        raise ValidationError(
            f"{param_name} is not a file: {filepath}"
        )
    
    return filepath


def check_directory_exists(
    dirpath: Union[str, Path],
    param_name: str = "directory"
) -> Path:
    """
    Check that directory exists.
    
    Parameters
    ----------
    dirpath : Union[str, Path]
        Directory path to check
    param_name : str
        Parameter name for error message
    
    Returns
    -------
    Path
        Validated directory path
    
    Raises
    ------
    ValidationError
        If directory doesn't exist or is not a directory
    
    Examples
    --------
    >>> check_directory_exists('salmon_index/')
    Path('salmon_index')
    
    >>> check_directory_exists('missing_dir/')
    ValidationError: directory not found: missing_dir/
    """
    from .errors import ValidationError
    
    dirpath = Path(dirpath)
    
    if not dirpath.exists():
        raise ValidationError(
            f"{param_name} not found: {dirpath}"
        )
    
    if not dirpath.is_dir():
        raise ValidationError(
            f"{param_name} is not a directory: {dirpath}"
        )
    
    return dirpath


def validate_salmon_index(
    index_path: Union[str, Path]
) -> Path:
    """
    Validate Salmon index directory.
    
    Checks that:
    1. Directory exists
    2. Contains versionInfo.json (indicates valid Salmon index)
    
    Parameters
    ----------
    index_path : Union[str, Path]
        Path to Salmon index directory
    
    Returns
    -------
    Path
        Validated index path
    
    Raises
    ------
    ValidationError
        If index is invalid or corrupted
    
    Examples
    --------
    >>> validate_salmon_index('salmon_index/')
    Path('salmon_index')
    
    Scientific Background
    ---------------------
    Salmon uses a quasi-mapping index that must contain:
    - versionInfo.json: Index metadata and Salmon version
    - pos.bin, ctable.bin: Index data structures
    
    Reference: Patro et al. 2017, Nature Methods
    """
    from .errors import ValidationError
    
    index_path = check_directory_exists(index_path, "Salmon index")
    
    # Check for versionInfo.json
    version_file = index_path / 'versionInfo.json'
    if not version_file.exists():
        raise ValidationError(
            f"Invalid Salmon index: {index_path}\n"
            f"Missing versionInfo.json - index may be corrupted or incomplete.\n"
            f"\n"
            f"Troubleshooting:\n"
            f"1. Check if index was built successfully\n"
            f"2. Rebuild with: salmon index -t transcripts.fa -i {index_path}\n"
            f"3. Ensure Salmon version matches index version"
        )
    
    logger.debug(f"✓ Valid Salmon index: {index_path}")
    return index_path


def validate_kallisto_index(
    index_path: Union[str, Path]
) -> Path:
    """
    Validate Kallisto index file.
    
    Parameters
    ----------
    index_path : Union[str, Path]
        Path to Kallisto index file (.idx)
    
    Returns
    -------
    Path
        Validated index path
    
    Raises
    ------
    ValidationError
        If index is invalid
    
    Examples
    --------
    >>> validate_kallisto_index('kallisto.idx')
    Path('kallisto.idx')
    
    Scientific Background
    ---------------------
    Kallisto uses a transcriptome de Bruijn graph index stored in .idx format.
    
    Reference: Bray et al. 2016, Nature Biotechnology
    """
    from .errors import ValidationError
    
    index_path = check_file_exists(index_path, "Kallisto index")
    
    # Kallisto index files should have .idx extension
    if index_path.suffix != '.idx':
        logger.warning(
            f"⚠️  Kallisto index typically has .idx extension: {index_path}\n"
            f"   If index is valid, you can ignore this warning."
        )
    
    # Check file size (Kallisto indexes are typically 5-15 GB)
    size_bytes = index_path.stat().st_size
    size_gb = size_bytes / (1024**3)
    
    if size_gb < 0.1:
        raise ValidationError(
            f"Kallisto index file is very small ({size_gb:.2f} GB): {index_path}\n"
            f"Index may be corrupted or incomplete.\n"
            f"Rebuild with: kallisto index -i {index_path} transcripts.fa"
        )
    
    logger.debug(f"✓ Valid Kallisto index: {index_path} ({size_gb:.2f} GB)")
    return index_path


def validate_fastq_file(
    fastq_path: Union[str, Path],
    param_name: str = "FASTQ file"
) -> Path:
    """
    Validate FASTQ file.
    
    Checks that:
    1. File exists
    2. Has valid FASTQ extension (.fastq, .fq, .fastq.gz, .fq.gz)
    
    Parameters
    ----------
    fastq_path : Union[str, Path]
        Path to FASTQ file
    param_name : str
        Parameter name for error message
    
    Returns
    -------
    Path
        Validated FASTQ path
    
    Raises
    ------
    ValidationError
        If FASTQ file is invalid
    
    Examples
    --------
    >>> validate_fastq_file('sample_R1.fastq.gz')
    Path('sample_R1.fastq.gz')
    
    >>> validate_fastq_file('sample.txt')
    ValidationError: FASTQ file has invalid extension: sample.txt
    """
    from .errors import ValidationError
    
    fastq_path = check_file_exists(fastq_path, param_name)
    
    # Check extension
    valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    has_valid_ext = any(str(fastq_path).endswith(ext) for ext in valid_extensions)
    
    if not has_valid_ext:
        raise ValidationError(
            f"{param_name} has invalid extension: {fastq_path}\n"
            f"Valid extensions: {', '.join(valid_extensions)}\n"
            f"\n"
            f"Hint: FASTQ files should be named like:\n"
            f"  - sample_R1.fastq.gz (compressed, recommended)\n"
            f"  - sample_R1.fastq (uncompressed)"
        )
    
    # Check file size (warn if very small or very large)
    size_bytes = fastq_path.stat().st_size
    size_mb = size_bytes / (1024**2)
    
    if size_mb < 1:
        logger.warning(
            f"⚠️  Very small FASTQ file ({size_mb:.2f} MB): {fastq_path}\n"
            f"   File may be empty or contain very few reads."
        )
    elif size_mb > 50000:  # > 50 GB
        logger.warning(
            f"⚠️  Very large FASTQ file ({size_mb/1024:.1f} GB): {fastq_path}\n"
            f"   Processing may require substantial resources."
        )
    
    return fastq_path


def validate_threads(
    threads: int,
    min_threads: int = 1,
    max_threads: int = 128
) -> int:
    """
    Validate number of threads.
    
    Parameters
    ----------
    threads : int
        Number of threads to use
    min_threads : int
        Minimum allowed threads (default: 1)
    max_threads : int
        Maximum allowed threads (default: 128)
    
    Returns
    -------
    int
        Validated thread count
    
    Raises
    ------
    ValidationError
        If threads out of valid range
    
    Examples
    --------
    >>> validate_threads(8)
    8
    
    >>> validate_threads(0)
    ValidationError: threads must be positive integer, got 0
    
    Notes
    -----
    Thread recommendations:
    - Salmon: 8-16 threads optimal, diminishing returns after 16
    - Kallisto: 4-8 threads optimal
    - STAR: 16-32 threads optimal
    - Generally: More threads helps, but check system availability
    """
    # First check it's a positive integer
    threads = validate_positive_integer(threads, 'threads')
    
    from .errors import ValidationError
    
    if threads < min_threads:
        raise ValidationError(
            f"threads must be at least {min_threads}, got {threads}"
        )
    
    if threads > max_threads:
        logger.warning(
            f"⚠️  threads ({threads}) exceeds recommended maximum ({max_threads}).\n"
            f"   This may not improve performance and could cause issues."
        )
    
    # Warn about system limits
    try:
        import multiprocessing
        available_cpus = multiprocessing.cpu_count()
        
        if threads > available_cpus:
            logger.warning(
                f"⚠️  threads ({threads}) exceeds available CPUs ({available_cpus}).\n"
                f"   Performance may be degraded by oversubscription."
            )
    except:
        pass  # Can't determine CPU count, skip warning
    
    return threads


# =============================================================================
# Keep existing constants at the end
# =============================================================================

# Constants for ensemble analysis (already in file)
# MIN_METHODS_ENSEMBLE = 2
# RECOMMENDED_METHODS_ENSEMBLE = 3
# IDEAL_METHODS_ENSEMBLE = 4
