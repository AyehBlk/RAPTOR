#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Data Profiler Module

Comprehensive RNA-seq data profiling for ML-based pipeline recommendation.

This module extracts statistical features from count matrices that are
predictive of optimal DE analysis pipeline performance.

============================================================================
LITERATURE BASIS FOR FEATURE SELECTION
============================================================================

The features extracted are based on extensive benchmarking studies:

1. SAMPLE SIZE & REPLICATION
   - Soneson & Delorenzi (2013): DESeq2/edgeR need ≥3 replicates
   - Law et al. (2014): limma-voom requires sufficient replication
   - Li et al. (2022): Wilcoxon test better for n > 8 per condition

2. DISPERSION / BIOLOGICAL COEFFICIENT OF VARIATION
   - Love et al. (2014): DESeq2 uses shrinkage for stable dispersion
   - Robinson et al. (2010): edgeR moderates toward common dispersion
   - BCV = sqrt(dispersion) represents biological variation
   - High BCV (>0.4): edgeR handles well
   - Low BCV (<0.2): limma-voom appropriate

3. LIBRARY SIZE VARIATION
   - High CV in library sizes → CPM normalization important
   - DESeq2 uses median-of-ratios (RLE)
   - edgeR uses TMM
   - Large variation (>5x) needs careful normalization

4. COUNT DISTRIBUTION
   - Low counts (<10): edgeR more sensitive
   - Medium/high counts: all methods comparable
   - Zero-inflation: negative binomial handles better

5. OUTLIERS
   - Small samples + outliers: edgeR_robust
   - Large samples + outliers: limma-voom robust
   - DESeq2 uses Cook's distance for outlier detection

6. SAMPLE SIZE THRESHOLDS (from Li et al. 2022, Genome Biology)
   - n < 8 per condition: parametric methods (DESeq2, edgeR)
   - n ≥ 8 per condition: Wilcoxon test has better FDR control
   - n > 20 per condition: limma-voom or Wilcoxon recommended

============================================================================
MATHEMATICAL FORMULATIONS
============================================================================

NEGATIVE BINOMIAL MODEL:
    Y ~ NB(μ, φ)
    E(Y) = μ
    Var(Y) = μ + φμ²
    
    where φ is the dispersion parameter
    BCV = sqrt(φ) = biological coefficient of variation

DISPERSION ESTIMATION (Method of Moments):
    φ̂ = (Var - μ) / μ²
    
MEAN-VARIANCE RELATIONSHIP:
    log(Var) = α + β·log(μ)
    
    - Poisson: β = 1
    - Overdispersed: β > 1
    - Typical RNA-seq: β ≈ 1.5-2

ZERO-INFLATION INDEX:
    ZI = (observed_zeros - expected_zeros) / expected_zeros
    expected_zeros = Σ exp(-μᵢ) × n_samples

============================================================================

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com  
Version: 2.2.0
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field, asdict
import warnings
import json
warnings.filterwarnings('ignore')

import logging


# RAPTOR v2.2.0 imports
from raptor.utils.validation import (
    validate_count_matrix,
    validate_metadata,
    validate_group_column,
    validate_file_path,
    validate_directory_path,
    validate_numeric_range,
    validate_positive_integer,
    validate_probability
)

from raptor.utils.errors import (
    handle_errors,
    ValidationError,
    DependencyError,
    check_file_exists,
    validate_output_writable
)

logger = logging.getLogger(__name__)


# =============================================================================
# CLI PARAMETERS & MODULE METADATA - v2.2.0
# =============================================================================

PROFILER_CLI_PARAMS = {
    'counts': {
        'flag': '--counts',
        'short': '-c',
        'required': True,
        'type': 'file',
        'help': 'Count matrix CSV file (genes x samples)',
        'validation': 'Must be CSV with gene IDs as rows, sample IDs as columns'
    },
    'metadata': {
        'flag': '--metadata',
        'short': '-m',
        'type': 'file',
        'help': 'Sample metadata CSV (optional, enables group-based features)',
        'validation': 'Must contain group/condition column matching sample IDs'
    },
    'group_column': {
        'flag': '--group-column',
        'short': '-g',
        'default': 'condition',
        'type': 'str',
        'help': 'Column name in metadata for grouping samples (default: condition)'
    },
    'output': {
        'flag': '--output',
        'short': '-o',
        'default': 'data_profile.json',
        'type': 'file',
        'help': 'Output profile JSON file'
    },
    'min_count': {
        'flag': '--min-count',
        'default': 1,
        'type': 'int',
        'range': [0, 1000],
        'help': 'Minimum count threshold for considering a gene expressed (default: 1)'
    },
    'verbose': {
        'flag': '--verbose',
        'short': '-v',
        'action': 'store_true',
        'help': 'Print detailed profiling information'
    }
}

# Module metadata
MODULE_NAME = "Data Profiler"
MODULE_VERSION = "2.2.0"
MODULE_ID = "M3"
MODULE_STAGE = 3

__all__ = [
    # Main class
    'RNAseqDataProfiler',
    
    # Data class
    'DataProfile',
    
    # Convenience functions
    'profile_data_quick',
    'get_key_characteristics',
    
    # CLI integration
    'PROFILER_CLI_PARAMS',
    'MODULE_NAME',
    'MODULE_VERSION',
    'MODULE_ID',
    'MODULE_STAGE'
]


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class DataProfile:
    """
    Comprehensive RNA-seq data profile.
    
    Contains all extracted features organized by category for
    ML-based pipeline recommendation.
    
    Attributes are grouped into:
    - Sample characteristics
    - Library size metrics
    - Gene detection metrics
    - Expression distribution
    - Dispersion estimates
    - Sparsity metrics
    - Count distribution
    - Mean-variance relationship
    - Quality indicators
    """
    
    # -------------------------------------------------------------------------
    # SAMPLE CHARACTERISTICS
    # -------------------------------------------------------------------------
    n_samples: int = 0
    n_genes: int = 0
    n_groups: int = 2
    group_sizes: List[int] = field(default_factory=list)
    min_group_size: int = 0
    max_group_size: int = 0
    mean_group_size: float = 0.0
    sample_balance: float = 1.0  # min/max ratio (1.0 = perfectly balanced)
    has_replicates: bool = True  # ≥2 per group
    has_sufficient_replicates: bool = True  # ≥3 per group
    design_complexity: str = 'simple'  # simple, moderate, complex
    
    # -------------------------------------------------------------------------
    # LIBRARY SIZE METRICS
    # -------------------------------------------------------------------------
    library_sizes: np.ndarray = field(default_factory=lambda: np.array([]))
    library_size_mean: float = 0.0
    library_size_median: float = 0.0
    library_size_std: float = 0.0
    library_size_cv: float = 0.0  # Coefficient of variation
    library_size_range: float = 0.0  # max/min ratio
    library_size_iqr: float = 0.0
    library_size_min: float = 0.0
    library_size_max: float = 0.0
    library_size_category: str = 'moderate'  # low, moderate, high, very_high
    sequencing_depth_category: str = 'standard'  # shallow, standard, deep
    
    # -------------------------------------------------------------------------
    # GENE DETECTION METRICS
    # -------------------------------------------------------------------------
    n_total_genes: int = 0
    n_expressed_genes: int = 0  # Any expression
    n_reliably_expressed: int = 0  # Mean count ≥ 10
    n_highly_expressed: int = 0  # Mean count ≥ 100
    detection_rate: float = 0.0  # Expressed / Total
    reliable_detection_rate: float = 0.0
    genes_per_sample_mean: float = 0.0
    genes_per_sample_std: float = 0.0
    genes_per_sample_cv: float = 0.0
    
    # -------------------------------------------------------------------------
    # EXPRESSION DISTRIBUTION (log2 scale)
    # -------------------------------------------------------------------------
    expression_mean: float = 0.0
    expression_median: float = 0.0
    expression_std: float = 0.0
    expression_variance: float = 0.0
    expression_skewness: float = 0.0
    expression_kurtosis: float = 0.0
    expression_iqr: float = 0.0
    expression_range: float = 0.0
    expression_percentiles: Dict[int, float] = field(default_factory=dict)
    
    # -------------------------------------------------------------------------
    # DISPERSION ESTIMATES (CRITICAL FOR PIPELINE SELECTION)
    # -------------------------------------------------------------------------
    # Gene-wise dispersion statistics
    dispersion_mean: float = 0.0
    dispersion_median: float = 0.0
    dispersion_std: float = 0.0
    dispersion_iqr: float = 0.0
    
    # Common/shared dispersion (DESeq2/edgeR style)
    common_dispersion: float = 0.0  # Trimmed mean
    trended_dispersion_intercept: float = 0.0
    trended_dispersion_slope: float = 0.0
    
    # Biological Coefficient of Variation
    bcv: float = 0.0  # sqrt(common_dispersion)
    bcv_category: str = 'moderate'  # low (<0.2), moderate (0.2-0.4), high (>0.4)
    
    # Overdispersion metrics
    overdispersion_ratio: float = 0.0  # Median(Var/Mean)
    overdispersion_category: str = 'moderate'
    
    # -------------------------------------------------------------------------
    # SPARSITY METRICS
    # -------------------------------------------------------------------------
    sparsity: float = 0.0  # Overall proportion of zeros
    sparsity_category: str = 'moderate'  # low (<0.3), moderate (0.3-0.6), high (>0.6)
    zero_inflation_index: float = 0.0  # Excess zeros vs Poisson expectation
    is_zero_inflated: bool = False
    dropout_rate: float = 0.0  # Per-sample average zero rate
    
    # -------------------------------------------------------------------------
    # COUNT DISTRIBUTION
    # -------------------------------------------------------------------------
    # Proportion of counts in each category
    zero_proportion: float = 0.0
    low_count_proportion: float = 0.0  # 1-9
    medium_count_proportion: float = 0.0  # 10-99
    high_count_proportion: float = 0.0  # 100-999
    very_high_count_proportion: float = 0.0  # ≥1000
    
    # Distribution characteristics
    count_mean: float = 0.0
    count_median: float = 0.0
    count_max: float = 0.0
    dynamic_range: float = 0.0  # log2(max/min nonzero)
    dynamic_range_category: str = 'moderate'
    
    # -------------------------------------------------------------------------
    # MEAN-VARIANCE RELATIONSHIP
    # -------------------------------------------------------------------------
    mean_var_slope: float = 0.0  # Slope of log(var) ~ log(mean)
    mean_var_intercept: float = 0.0
    mean_var_r_squared: float = 0.0  # Fit quality
    poisson_fit_score: float = 0.0  # How close to Poisson (slope=1)
    nb_fit_score: float = 0.0  # How well NB fits
    
    # -------------------------------------------------------------------------
    # QUALITY INDICATORS (from QC module if available)
    # -------------------------------------------------------------------------
    quality_score: float = 100.0
    has_outliers: bool = False
    n_outliers: int = 0
    outlier_severity: str = 'none'  # none, mild, moderate, severe
    has_batch_effect: bool = False
    batch_strength: float = 0.0
    batch_confounded: bool = False
    
    # -------------------------------------------------------------------------
    # DERIVED METRICS FOR ML
    # -------------------------------------------------------------------------
    feature_vector: np.ndarray = field(default_factory=lambda: np.array([]))
    feature_names: List[str] = field(default_factory=list)
    
    # -------------------------------------------------------------------------
    # METADATA
    # -------------------------------------------------------------------------
    profiler_version: str = '2.2.0'
    profile_timestamp: str = ''
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary, handling numpy arrays."""
        result = {}
        for key, value in asdict(self).items():
            if isinstance(value, np.ndarray):
                result[key] = value.tolist()
            else:
                result[key] = value
        return result
    
    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "",
            "╔" + "═" * 68 + "╗",
            "║" + "  🦖 RAPTOR DATA PROFILE".center(68) + "║",
            "╠" + "═" * 68 + "╣",
            "",
            "  📊 DIMENSIONS",
            "  " + "─" * 40,
            f"    Genes:    {self.n_genes:>10,}",
            f"    Samples:  {self.n_samples:>10}",
            f"    Groups:   {self.n_groups:>10}",
            f"    Balance:  {self.sample_balance:>10.2f}",
            "",
            "  📚 LIBRARY SIZE",
            "  " + "─" * 40,
            f"    Mean:     {self.library_size_mean:>12,.0f} reads",
            f"    Median:   {self.library_size_median:>12,.0f} reads",
            f"    CV:       {self.library_size_cv:>12.3f}",
            f"    Range:    {self.library_size_range:>12.1f}x",
            f"    Category: {self.library_size_category:>12}",
            "",
            "  🧬 EXPRESSION",
            "  " + "─" * 40,
            f"    Detected genes:   {self.n_expressed_genes:>8,} ({self.detection_rate:.1%})",
            f"    Reliable (≥10):   {self.n_reliably_expressed:>8,} ({self.reliable_detection_rate:.1%})",
            f"    Mean (log2):      {self.expression_mean:>8.2f}",
            f"    Sparsity:         {self.sparsity:>8.1%}",
            "",
            "  📈 DISPERSION (Critical for Pipeline Selection)",
            "  " + "─" * 40,
            f"    Common φ:         {self.common_dispersion:>8.4f}",
            f"    BCV:              {self.bcv:>8.3f} ({self.bcv*100:.1f}% biological variation)",
            f"    BCV Category:     {self.bcv_category:>8}",
            f"    Overdispersion:   {self.overdispersion_ratio:>8.2f}x",
            "",
            "  📊 COUNT DISTRIBUTION",
            "  " + "─" * 40,
            f"    Zero:             {self.zero_proportion:>8.1%}",
            f"    Low (1-9):        {self.low_count_proportion:>8.1%}",
            f"    Medium (10-99):   {self.medium_count_proportion:>8.1%}",
            f"    High (100-999):   {self.high_count_proportion:>8.1%}",
            f"    Very High (≥1000):{self.very_high_count_proportion:>8.1%}",
            f"    Dynamic Range:    {self.dynamic_range:>8.1f} log2 units",
            "",
            "  📐 MEAN-VARIANCE RELATIONSHIP",
            "  " + "─" * 40,
            f"    Slope:            {self.mean_var_slope:>8.3f} (Poisson=1, NB≈1.5-2)",
            f"    R²:               {self.mean_var_r_squared:>8.3f}",
            "",
            "╚" + "═" * 68 + "╝",
            ""
        ]
        return "\n".join(lines)
    
    def get_recommendation_features(self) -> Dict[str, Any]:
        """Get features most relevant for pipeline recommendation."""
        return {
            # Sample characteristics
            'n_samples': self.n_samples,
            'n_groups': self.n_groups,
            'min_group_size': self.min_group_size,
            'sample_balance': self.sample_balance,
            'has_sufficient_replicates': self.has_sufficient_replicates,
            
            # Dispersion (most critical)
            'bcv': self.bcv,
            'bcv_category': self.bcv_category,
            'common_dispersion': self.common_dispersion,
            'overdispersion_ratio': self.overdispersion_ratio,
            
            # Library size
            'library_size_cv': self.library_size_cv,
            'library_size_range': self.library_size_range,
            
            # Count distribution
            'low_count_proportion': self.low_count_proportion,
            'sparsity': self.sparsity,
            'dynamic_range': self.dynamic_range,
            
            # Mean-variance
            'mean_var_slope': self.mean_var_slope,
            
            # Quality
            'has_outliers': self.has_outliers,
            'outlier_severity': self.outlier_severity,
            'has_batch_effect': self.has_batch_effect,
            'quality_score': self.quality_score
        }


# =============================================================================
# MAIN PROFILER CLASS
# =============================================================================

class RNAseqDataProfiler:
    """
    Comprehensive RNA-seq data profiler for pipeline recommendation.
    
    Extracts statistical features that are predictive of differential
    expression pipeline performance, based on benchmark literature.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix with genes as rows and samples as columns.
        Raw counts (integers) are expected.
    metadata : pd.DataFrame, optional
        Sample metadata with at least a group/condition column.
    group_column : str, default='condition'
        Column name in metadata containing group labels.
    min_count_threshold : int, default=1
        Minimum count to consider a gene "expressed".
    
    Examples
    --------
    >>> # Basic usage
    >>> profiler = RNAseqDataProfiler(counts, metadata)
    >>> profile = profiler.run_full_profile()
    >>> print(profile.summary())
    
    >>> # Get recommendation features
    >>> features = profile.get_recommendation_features()
    >>> print(f"BCV: {features['bcv']:.3f}")
    >>> print(f"Min group size: {features['min_group_size']}")
    
    >>> # Access feature vector for ML
    >>> X = profile.feature_vector
    """
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame = None,
        group_column: str = 'condition',
        min_count_threshold: int = 1
    ):
        """Initialize profiler with count data."""
        
        # Validate inputs
        validate_count_matrix(counts, allow_negative=False, min_genes=10, min_samples=2)
        
        if metadata is not None:
            validate_metadata(metadata, counts)
            if group_column in metadata.columns:
                validate_group_column(metadata, group_column, min_groups=1, min_samples_per_group=1)
        
        if min_count_threshold < 0:
            raise ValidationError(f"min_count_threshold must be >= 0, got {min_count_threshold}")
        # Store parameters
        self.group_column = group_column
        self.min_count_threshold = min_count_threshold
        
        # Process counts
        self.counts = counts.copy()
        self.counts = self.counts.apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Convert to numpy for faster computation
        self.count_matrix = self.counts.values.astype(np.float64)
        
        # Dimensions
        self.n_genes, self.n_samples = self.count_matrix.shape
        self.gene_names = list(self.counts.index)
        self.sample_names = list(self.counts.columns)
        
        # Calculate library sizes
        self.library_sizes = self.count_matrix.sum(axis=0)
        
        # Log-transformed counts
        self.log_counts = np.log2(self.count_matrix + 1)
        
        # Process metadata
        self.metadata = metadata
        self.groups = None
        self.n_groups = 1
        self.group_labels = None
        
        if metadata is not None and group_column in metadata.columns:
            # Align metadata with samples
            if 'sample' in metadata.columns:
                meta_aligned = metadata.set_index('sample').loc[self.sample_names]
                self.groups = meta_aligned[group_column].values
            else:
                self.groups = metadata[group_column].values[:self.n_samples]
            
            self.group_labels, self.group_indices = np.unique(
                self.groups, return_inverse=True
            )
            self.n_groups = len(self.group_labels)
        
        logger.info(f"Profiler initialized: {self.n_genes:,} genes × {self.n_samples} samples")
    
    def run_full_profile(self) -> DataProfile:
        """
        Run comprehensive data profiling.
        
        Extracts all features relevant for pipeline recommendation.
        
        Returns
        -------
        DataProfile
            Complete profile with all extracted features.
        """
        logger.info("Starting full data profiling...")
        
        # Initialize profile
        profile = DataProfile()
        profile.profile_timestamp = pd.Timestamp.now().isoformat()
        
        # Run all profiling methods
        self._profile_samples(profile)
        self._profile_library_sizes(profile)
        self._profile_gene_detection(profile)
        self._profile_expression_distribution(profile)
        self._profile_dispersion(profile)
        self._profile_sparsity(profile)
        self._profile_count_distribution(profile)
        self._profile_mean_variance(profile)
        
        # Build feature vector for ML
        self._build_feature_vector(profile)
        
        logger.info("Data profiling complete")
        return profile
    
    # =========================================================================
    # PROFILING METHODS
    # =========================================================================
    
    def _profile_samples(self, profile: DataProfile):
        """Profile sample and group characteristics."""
        logger.debug("Profiling sample characteristics...")
        
        profile.n_samples = self.n_samples
        profile.n_genes = self.n_genes
        profile.n_total_genes = self.n_genes
        profile.n_groups = self.n_groups
        
        if self.groups is not None:
            _, group_counts = np.unique(self.groups, return_counts=True)
            profile.group_sizes = group_counts.tolist()
            profile.min_group_size = int(group_counts.min())
            profile.max_group_size = int(group_counts.max())
            profile.mean_group_size = float(group_counts.mean())
            profile.sample_balance = profile.min_group_size / profile.max_group_size
            profile.has_replicates = profile.min_group_size >= 2
            profile.has_sufficient_replicates = profile.min_group_size >= 3
        else:
            profile.group_sizes = [self.n_samples]
            profile.min_group_size = self.n_samples
            profile.max_group_size = self.n_samples
            profile.mean_group_size = float(self.n_samples)
            profile.sample_balance = 1.0
            profile.has_replicates = self.n_samples >= 2
            profile.has_sufficient_replicates = self.n_samples >= 3
        
        # Design complexity
        if self.n_groups <= 2 and profile.sample_balance >= 0.8:
            profile.design_complexity = 'simple'
        elif self.n_groups <= 4 and profile.sample_balance >= 0.5:
            profile.design_complexity = 'moderate'
        else:
            profile.design_complexity = 'complex'
    
    def _profile_library_sizes(self, profile: DataProfile):
        """Profile library size characteristics."""
        logger.debug("Profiling library sizes...")
        
        lib = self.library_sizes
        
        profile.library_sizes = lib
        profile.library_size_mean = float(np.mean(lib))
        profile.library_size_median = float(np.median(lib))
        profile.library_size_std = float(np.std(lib))
        profile.library_size_min = float(np.min(lib))
        profile.library_size_max = float(np.max(lib))
        
        # Coefficient of variation
        if profile.library_size_mean > 0:
            profile.library_size_cv = float(np.std(lib) / np.mean(lib))
        
        # Range (fold difference)
        if profile.library_size_min > 0:
            profile.library_size_range = float(np.max(lib) / np.min(lib))
        else:
            profile.library_size_range = float('inf')
        
        # IQR
        q75, q25 = np.percentile(lib, [75, 25])
        profile.library_size_iqr = float(q75 - q25)
        
        # Categorize library size variation
        if profile.library_size_cv < 0.1:
            profile.library_size_category = 'low'
        elif profile.library_size_cv < 0.3:
            profile.library_size_category = 'moderate'
        elif profile.library_size_cv < 0.5:
            profile.library_size_category = 'high'
        else:
            profile.library_size_category = 'very_high'
        
        # Sequencing depth category
        median_millions = profile.library_size_median / 1e6
        if median_millions < 5:
            profile.sequencing_depth_category = 'shallow'
        elif median_millions < 30:
            profile.sequencing_depth_category = 'standard'
        else:
            profile.sequencing_depth_category = 'deep'
    
    def _profile_gene_detection(self, profile: DataProfile):
        """Profile gene detection rates."""
        logger.debug("Profiling gene detection...")
        
        counts = self.count_matrix
        
        # Gene-level statistics
        gene_sums = counts.sum(axis=1)
        gene_means = counts.mean(axis=1)
        
        # Detection counts
        profile.n_expressed_genes = int((gene_sums > 0).sum())
        profile.n_reliably_expressed = int((gene_means >= 10).sum())
        profile.n_highly_expressed = int((gene_means >= 100).sum())
        
        # Detection rates
        profile.detection_rate = profile.n_expressed_genes / self.n_genes
        profile.reliable_detection_rate = profile.n_reliably_expressed / self.n_genes
        
        # Genes per sample
        genes_per_sample = (counts > 0).sum(axis=0)
        profile.genes_per_sample_mean = float(np.mean(genes_per_sample))
        profile.genes_per_sample_std = float(np.std(genes_per_sample))
        if profile.genes_per_sample_mean > 0:
            profile.genes_per_sample_cv = float(
                profile.genes_per_sample_std / profile.genes_per_sample_mean
            )
    
    def _profile_expression_distribution(self, profile: DataProfile):
        """Profile expression level distribution."""
        logger.debug("Profiling expression distribution...")
        
        # Use log2(counts + 1) for non-zero entries
        log_vals = self.log_counts.flatten()
        nonzero_log = log_vals[log_vals > 0]
        
        if len(nonzero_log) > 0:
            profile.expression_mean = float(np.mean(nonzero_log))
            profile.expression_median = float(np.median(nonzero_log))
            profile.expression_std = float(np.std(nonzero_log))
            profile.expression_variance = float(np.var(nonzero_log))
            
            # Shape statistics
            try:
                profile.expression_skewness = float(stats.skew(nonzero_log))
                profile.expression_kurtosis = float(stats.kurtosis(nonzero_log))
            except:
                pass
            
            # IQR and range
            q75, q25 = np.percentile(nonzero_log, [75, 25])
            profile.expression_iqr = float(q75 - q25)
            profile.expression_range = float(nonzero_log.max() - nonzero_log.min())
            
            # Percentiles
            for p in [5, 10, 25, 50, 75, 90, 95]:
                profile.expression_percentiles[p] = float(np.percentile(nonzero_log, p))
    
    def _profile_dispersion(self, profile: DataProfile):
        """
        Estimate dispersion parameters.
        
        This is the most critical feature for pipeline selection.
        
        The negative binomial dispersion φ relates variance to mean:
            Var(Y) = μ + φ·μ²
        
        BCV (Biological Coefficient of Variation) = sqrt(φ)
        
        Method: Uses method-of-moments estimation with trimming
        for robustness, similar to DESeq2/edgeR approaches.
        """
        logger.debug("Profiling dispersion...")
        
        counts = self.count_matrix
        
        # Gene-wise statistics
        gene_means = counts.mean(axis=1)
        gene_vars = counts.var(axis=1, ddof=1)
        
        # Filter to reliably expressed genes (mean ≥ 5)
        # This is similar to filterByExpr in edgeR
        expr_mask = gene_means >= 5
        
        if expr_mask.sum() < 50:
            logger.warning("Too few expressed genes for reliable dispersion estimation")
            # Use all genes with mean ≥ 1
            expr_mask = gene_means >= 1
        
        if expr_mask.sum() < 10:
            logger.error("Insufficient data for dispersion estimation")
            return
        
        # Filtered statistics
        means = gene_means[expr_mask]
        variances = gene_vars[expr_mask]
        
        # Method of moments dispersion estimation
        # φ = (Var - μ) / μ²
        # Only valid where variance > mean (overdispersed)
        valid = (variances > means) & (means > 0) & (means**2 > 0)
        
        if valid.sum() < 10:
            # Use all genes if too few pass the filter
            valid = means > 0
        
        dispersions = np.zeros(len(means))
        dispersions[valid] = (variances[valid] - means[valid]) / (means[valid] ** 2)
        
        # Clip to reasonable range
        dispersions = np.clip(dispersions, 0.001, 10)
        
        # Gene-wise dispersion statistics
        profile.dispersion_mean = float(np.mean(dispersions[dispersions > 0]))
        profile.dispersion_median = float(np.median(dispersions[dispersions > 0]))
        profile.dispersion_std = float(np.std(dispersions[dispersions > 0]))
        
        q75, q25 = np.percentile(dispersions[dispersions > 0], [75, 25])
        profile.dispersion_iqr = float(q75 - q25)
        
        # Common dispersion (trimmed mean for robustness)
        sorted_disp = np.sort(dispersions[dispersions > 0])
        n_disp = len(sorted_disp)
        trim_frac = 0.1  # 10% trimming on each side
        trim_n = int(n_disp * trim_frac)
        
        if n_disp > 2 * trim_n + 1:
            trimmed = sorted_disp[trim_n:-trim_n]
            profile.common_dispersion = float(np.mean(trimmed))
        else:
            profile.common_dispersion = float(np.median(sorted_disp))
        
        # BCV = sqrt(common dispersion)
        profile.bcv = float(np.sqrt(profile.common_dispersion))
        
        # Categorize BCV
        if profile.bcv < 0.2:
            profile.bcv_category = 'low'
        elif profile.bcv < 0.4:
            profile.bcv_category = 'moderate'
        else:
            profile.bcv_category = 'high'
        
        # Dispersion-mean trend (log-log scale)
        log_means = np.log10(means[valid] + 1)
        log_disp = np.log10(dispersions[valid] + 0.001)
        
        try:
            slope, intercept, r_value, _, _ = stats.linregress(log_means, log_disp)
            profile.trended_dispersion_slope = float(slope)
            profile.trended_dispersion_intercept = float(intercept)
        except:
            pass
        
        # Overdispersion ratio: median(var/mean)
        var_mean_ratio = np.zeros(len(means))
        nonzero = means > 0
        var_mean_ratio[nonzero] = variances[nonzero] / means[nonzero]
        
        profile.overdispersion_ratio = float(np.median(var_mean_ratio[var_mean_ratio > 0]))
        
        # Categorize overdispersion
        if profile.overdispersion_ratio < 2:
            profile.overdispersion_category = 'low'
        elif profile.overdispersion_ratio < 10:
            profile.overdispersion_category = 'moderate'
        else:
            profile.overdispersion_category = 'high'
    
    def _profile_sparsity(self, profile: DataProfile):
        """Profile sparsity and zero-inflation."""
        logger.debug("Profiling sparsity...")
        
        counts = self.count_matrix
        total_entries = counts.size
        
        # Overall sparsity
        n_zeros = (counts == 0).sum()
        profile.sparsity = float(n_zeros / total_entries)
        profile.zero_proportion = profile.sparsity
        
        # Categorize sparsity
        if profile.sparsity < 0.3:
            profile.sparsity_category = 'low'
        elif profile.sparsity < 0.6:
            profile.sparsity_category = 'moderate'
        else:
            profile.sparsity_category = 'high'
        
        # Zero-inflation index
        # Compare observed zeros to expected under Poisson
        gene_means = counts.mean(axis=1)
        # Expected zeros under Poisson: sum of P(Y=0) = exp(-λ) for each gene
        expected_zeros = np.sum(np.exp(-gene_means) * self.n_samples)
        
        if expected_zeros > 0:
            profile.zero_inflation_index = float(
                (n_zeros - expected_zeros) / expected_zeros
            )
        
        # Is zero-inflated? (more than 50% excess zeros)
        profile.is_zero_inflated = profile.zero_inflation_index > 0.5
        
        # Dropout rate (per-sample zero proportion)
        zeros_per_sample = (counts == 0).sum(axis=0)
        profile.dropout_rate = float(np.mean(zeros_per_sample / self.n_genes))
    
    def _profile_count_distribution(self, profile: DataProfile):
        """Profile count value distribution."""
        logger.debug("Profiling count distribution...")
        
        counts = self.count_matrix.flatten()
        total = len(counts)
        nonzero = counts[counts > 0]
        
        # Count categories (as proportion of all entries)
        profile.zero_proportion = float((counts == 0).sum() / total)
        profile.low_count_proportion = float(((counts >= 1) & (counts < 10)).sum() / total)
        profile.medium_count_proportion = float(((counts >= 10) & (counts < 100)).sum() / total)
        profile.high_count_proportion = float(((counts >= 100) & (counts < 1000)).sum() / total)
        profile.very_high_count_proportion = float((counts >= 1000).sum() / total)
        
        # Basic statistics
        profile.count_mean = float(np.mean(nonzero)) if len(nonzero) > 0 else 0
        profile.count_median = float(np.median(nonzero)) if len(nonzero) > 0 else 0
        profile.count_max = float(np.max(counts))
        
        # Dynamic range
        if len(nonzero) > 1 and nonzero.min() > 0:
            profile.dynamic_range = float(np.log2(nonzero.max() / nonzero.min()))
        
        # Categorize dynamic range
        if profile.dynamic_range < 10:
            profile.dynamic_range_category = 'low'
        elif profile.dynamic_range < 20:
            profile.dynamic_range_category = 'moderate'
        else:
            profile.dynamic_range_category = 'high'
    
    def _profile_mean_variance(self, profile: DataProfile):
        """
        Profile mean-variance relationship.
        
        This is fundamental to RNA-seq modeling:
        - Poisson: Var = μ (slope = 1 on log-log)
        - Negative Binomial: Var = μ + φμ² (slope > 1)
        - Over-dispersed: slope typically 1.5-2
        """
        logger.debug("Profiling mean-variance relationship...")
        
        counts = self.count_matrix
        
        # Gene-wise statistics
        gene_means = counts.mean(axis=1)
        gene_vars = counts.var(axis=1, ddof=1)
        
        # Filter to expressed genes
        mask = (gene_means > 1) & (gene_vars > 0)
        
        if mask.sum() < 50:
            logger.warning("Too few genes for mean-variance analysis")
            return
        
        means = gene_means[mask]
        variances = gene_vars[mask]
        
        # Log-log regression: log(var) = α + β·log(mean)
        log_means = np.log10(means)
        log_vars = np.log10(variances)
        
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                log_means, log_vars
            )
            
            profile.mean_var_slope = float(slope)
            profile.mean_var_intercept = float(intercept)
            profile.mean_var_r_squared = float(r_value ** 2)
            
            # Poisson fit score: how close is slope to 1?
            # Score of 1 means perfect Poisson, decays as slope deviates
            profile.poisson_fit_score = float(np.exp(-abs(slope - 1)))
            
            # NB fit score: based on R² and slope in expected range
            if 1.3 <= slope <= 2.2:
                profile.nb_fit_score = float(r_value ** 2)
            else:
                profile.nb_fit_score = float(r_value ** 2 * 0.8)
                
        except Exception as e:
            logger.warning(f"Mean-variance regression failed: {e}")
    
    def _build_feature_vector(self, profile: DataProfile):
        """Build feature vector for ML pipeline recommendation."""
        logger.debug("Building feature vector...")
        
        # Define feature names and values
        feature_specs = [
            # Sample characteristics (5 features)
            ('log2_n_samples', np.log2(profile.n_samples + 1)),
            ('log2_n_genes', np.log2(profile.n_genes + 1)),
            ('n_groups', profile.n_groups),
            ('min_group_size', profile.min_group_size),
            ('sample_balance', profile.sample_balance),
            
            # Library size (4 features)
            ('library_size_cv', profile.library_size_cv),
            ('log2_library_size_range', np.log2(profile.library_size_range + 1) if profile.library_size_range != float('inf') else 10),
            ('log2_library_size_mean', np.log2(profile.library_size_mean + 1)),
            ('library_size_iqr_norm', profile.library_size_iqr / (profile.library_size_mean + 1)),
            
            # Gene detection (3 features)
            ('detection_rate', profile.detection_rate),
            ('reliable_detection_rate', profile.reliable_detection_rate),
            ('genes_per_sample_cv', profile.genes_per_sample_cv),
            
            # Expression distribution (4 features)
            ('expression_mean', profile.expression_mean),
            ('expression_variance', profile.expression_variance),
            ('expression_skewness', profile.expression_skewness),
            ('expression_iqr', profile.expression_iqr),
            
            # Dispersion (CRITICAL - 5 features)
            ('bcv', profile.bcv),
            ('common_dispersion', profile.common_dispersion),
            ('dispersion_trend_slope', profile.trended_dispersion_slope),
            ('overdispersion_ratio', np.log2(profile.overdispersion_ratio + 1)),
            ('dispersion_iqr', profile.dispersion_iqr),
            
            # Sparsity (3 features)
            ('sparsity', profile.sparsity),
            ('zero_inflation_index', np.clip(profile.zero_inflation_index, -5, 5)),
            ('dropout_rate', profile.dropout_rate),
            
            # Count distribution (4 features)
            ('low_count_proportion', profile.low_count_proportion),
            ('medium_count_proportion', profile.medium_count_proportion),
            ('high_count_proportion', profile.high_count_proportion),
            ('dynamic_range', np.log2(profile.dynamic_range + 1)),
            
            # Mean-variance (3 features)
            ('mean_var_slope', profile.mean_var_slope),
            ('mean_var_r_squared', profile.mean_var_r_squared),
            ('poisson_fit_score', profile.poisson_fit_score),
            
            # Quality (1 feature)
            ('quality_score', profile.quality_score / 100.0),
        ]
        
        # Extract names and values
        names = [f[0] for f in feature_specs]
        values = np.array([f[1] for f in feature_specs], dtype=np.float32)
        
        # Handle NaN/Inf
        values = np.nan_to_num(values, nan=0.0, posinf=10.0, neginf=-10.0)
        
        profile.feature_names = names
        profile.feature_vector = values


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def profile_data_quick(
    counts: pd.DataFrame,
    metadata: pd.DataFrame = None,
    group_column: str = 'condition'
) -> DataProfile:
    """
    Quick data profiling with default settings.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes × samples)
    metadata : pd.DataFrame, optional
        Sample metadata
    group_column : str
        Column name for groups
    
    Returns
    -------
    DataProfile
        Complete data profile
    
    Examples
    --------
    >>> profile = profile_data_quick(counts, metadata)
    >>> print(f"BCV: {profile.bcv:.3f}")
    >>> print(f"Sample size per group: {profile.min_group_size}")
    """
    profiler = RNAseqDataProfiler(counts, metadata, group_column)
    return profiler.run_full_profile()


def get_key_characteristics(counts: pd.DataFrame) -> Dict[str, Any]:
    """
    Get key data characteristics without full profiling.
    
    Faster than full profiling, good for quick assessment.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix
    
    Returns
    -------
    Dict
        Key characteristics
    """
    counts_array = counts.values.astype(float)
    n_genes, n_samples = counts_array.shape
    
    lib_sizes = counts_array.sum(axis=0)
    gene_means = counts_array.mean(axis=1)
    gene_vars = counts_array.var(axis=1, ddof=1)
    
    # Quick dispersion estimate
    expressed = gene_means >= 5
    if expressed.sum() > 0:
        means = gene_means[expressed]
        variances = gene_vars[expressed]
        valid = (variances > means) & (means > 0)
        if valid.sum() > 0:
            dispersions = (variances[valid] - means[valid]) / (means[valid] ** 2)
            dispersions = np.clip(dispersions, 0.001, 10)
            bcv = float(np.sqrt(np.median(dispersions)))
        else:
            bcv = 0.3  # Default
    else:
        bcv = 0.3
    
    nonzero = counts_array[counts_array > 0]
    
    return {
        'n_genes': n_genes,
        'n_samples': n_samples,
        'library_size_mean': float(lib_sizes.mean()),
        'library_size_cv': float(lib_sizes.std() / lib_sizes.mean()) if lib_sizes.mean() > 0 else 0,
        'library_size_range': float(lib_sizes.max() / lib_sizes.min()) if lib_sizes.min() > 0 else float('inf'),
        'sparsity': float((counts_array == 0).sum() / counts_array.size),
        'bcv': bcv,
        'median_count': float(np.median(nonzero)) if len(nonzero) > 0 else 0,
        'dynamic_range': float(np.log2(nonzero.max() / nonzero.min())) if len(nonzero) > 1 and nonzero.min() > 0 else 0,
        'detection_rate': float((gene_means > 0).sum() / n_genes)
    }


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("""
╔═══════════════════════════════════════════════════════════════════════════╗
║              🦖 RAPTOR Data Profiler Module v2.2.0                        ║
╠═══════════════════════════════════════════════════════════════════════════╣
║                                                                           ║
║  Comprehensive RNA-seq data profiling for ML pipeline recommendation.    ║
║                                                                           ║
║  FEATURES EXTRACTED (32 total):                                          ║
║  ─────────────────────────────                                           ║
║  Sample characteristics:     5 features                                  ║
║  Library size metrics:       4 features                                  ║
║  Gene detection:            3 features                                  ║
║  Expression distribution:    4 features                                  ║
║  Dispersion (CRITICAL):      5 features                                  ║
║  Sparsity metrics:          3 features                                  ║
║  Count distribution:         4 features                                  ║
║  Mean-variance:             3 features                                  ║
║  Quality:                   1 feature                                   ║
║                                                                           ║
║  USAGE:                                                                   ║
║  ──────                                                                   ║
║  from raptor.profiler import RNAseqDataProfiler, profile_data_quick      ║
║                                                                           ║
║  profiler = RNAseqDataProfiler(counts, metadata, 'condition')            ║
║  profile = profiler.run_full_profile()                                   ║
║  print(profile.summary())                                                ║
║                                                                           ║
║  # Quick access                                                          ║
║  profile = profile_data_quick(counts, metadata)                          ║
║  features = profile.get_recommendation_features()                        ║
║  X = profile.feature_vector  # For ML model                              ║
║                                                                           ║
║  LITERATURE BASIS:                                                        ║
║  ─────────────────                                                        ║
║  • Love et al. (2014) - DESeq2 paper                                     ║
║  • Robinson et al. (2010) - edgeR paper                                  ║
║  • Law et al. (2014) - voom paper                                        ║
║  • Li et al. (2022) - Pipeline comparison, Genome Biology                ║
║  • Soneson & Delorenzi (2013) - Method comparison                        ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
    """)
