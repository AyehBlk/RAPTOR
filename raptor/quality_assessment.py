#!/usr/bin/env python3

"""
Advanced Data Quality Assessment Module

Provides comprehensive quality metrics including:
- Batch effect detection
- Overall quality scoring
- Sample clustering analysis
- Advanced Outlier detection (6 methods)
- Technical variation assessment
- Biological signal strength
- Normalization options (log2, CPM, quantile, none)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import warnings
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

QUALITY_ASSESSMENT_CLI_PARAMS = {
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
        'help': 'Sample metadata CSV (optional, enables batch effect detection)',
        'validation': 'Must contain sample_id column matching count matrix'
    },
    'output': {
        'flag': '--output',
        'short': '-o',
        'default': 'quality_report.json',
        'type': 'file',
        'help': 'Output report JSON file'
    },
    'normalization': {
        'flag': '--normalization',
        'short': '-n',
        'default': 'log2',
        'choices': ['log2', 'cpm', 'quantile', 'none'],
        'help': 'Normalization method (log2=default, cpm=varying library sizes, quantile=aggressive, none=pre-normalized)'
    },
    'consensus_threshold': {
        'flag': '--consensus-threshold',
        'short': '-t',
        'default': 3,
        'type': 'int',
        'range': [1, 6],
        'help': 'Outlier consensus threshold: how many of 6 methods must flag a sample (default: 3)'
    },
    'plot': {
        'flag': '--plot',
        'action': 'store_true',
        'help': 'Generate QC visualization plots'
    },
    'plot_output': {
        'flag': '--plot-output',
        'default': 'qc_plots.pdf',
        'type': 'file',
        'help': 'Output file for plots (PDF)'
    }
}

# Module metadata
MODULE_NAME = "Quality Assessment"
MODULE_VERSION = "2.2.0"
MODULE_ID = "M2"
MODULE_STAGE = 2

__all__ = [
    # Main class
    'DataQualityAssessor',
    
    # Convenience functions
    'quick_quality_check',
    'detect_outliers_quick',
    
    # Data classes
    'OutlierResult',
    'QualityReport',
    'BatchEffectResult',
    
    # CLI integration
    'QUALITY_ASSESSMENT_CLI_PARAMS',
    'MODULE_NAME',
    'MODULE_VERSION',
    'MODULE_ID',
    'MODULE_STAGE'
]


# =============================================================================
# DATA CLASSES FOR RESULTS
# =============================================================================

@dataclass
class OutlierResult:
    """Container for advanced outlier detection results."""
    outlier_samples: List[str]
    outlier_scores: Dict[str, int]  # Sample -> number of methods flagging it
    method_results: Dict[str, List[str]]  # Method -> list of outliers
    consensus_threshold: int
    n_outliers: int
    outlier_percentage: float
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 60,
            "ADVANCED OUTLIER DETECTION RESULTS",
            "=" * 60,
            f"Consensus threshold: {self.consensus_threshold} methods",
            f"Total outliers: {self.n_outliers} ({self.outlier_percentage:.1f}%)",
            "",
        ]
        
        if self.outlier_samples:
            lines.append("Outlier samples:")
            for sample in self.outlier_samples:
                score = self.outlier_scores.get(sample, 0)
                lines.append(f"  • {sample} (flagged by {score} methods)")
        else:
            lines.append("No consensus outliers detected ✓")
        
        lines.extend(["", "Detection by method:"])
        for method, samples in self.method_results.items():
            lines.append(f"  • {method}: {len(samples)} outliers")
        
        lines.append("=" * 60)
        return "\n".join(lines)


# =============================================================================
# MAIN CLASS
# =============================================================================

class DataQualityAssessor:
    """
    Comprehensive data quality assessment for RNA-seq data.
    
    Performs:
    - Library quality assessment
    - Gene detection analysis
    - Advanced outlier detection (6 methods)  # UPGRADED v2.2.0
    - Batch effect detection (improved in v2.2.0)
    - Variance structure analysis
    - Biological signal assessment
    - Overall quality scoring (0-100)
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes x samples)
    metadata : pd.DataFrame, optional
        Sample metadata with batch/condition information
    normalization : str, default='log2'
        Normalization method: 'log2', 'cpm', 'quantile', 'none'
    
    Examples
    --------
    >>> assessor = DataQualityAssessor(counts, metadata)
    >>> quality_report = assessor.assess_quality()
    >>> print(f"Overall Quality Score: {quality_report['overall']['score']:.1f}/100")
    >>> 
    >>> # Advanced outlier detection (6 methods)
    >>> outlier_result = assessor.detect_outliers_advanced(consensus_threshold=3)
    >>> print(outlier_result.summary())
    >>> 
    >>> # Check batch effect recommendations
    >>> if quality_report['batch_effects']['batch_detected']:
    ...     print(quality_report['batch_effects']['recommendation'])
    """
    
    # Valid normalization methods
    NORMALIZATION_METHODS = ['log2', 'cpm', 'quantile', 'none']
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self, 
        counts: pd.DataFrame, 
        metadata: pd.DataFrame = None,
        normalization: str = 'log2'
    ):
        """
        Initialize quality assessor.
        
        Parameters
        ----------
        counts : pd.DataFrame
            Count matrix (genes x samples). Should be raw counts for most
            normalization methods, or pre-normalized data if normalization='none'.
        metadata : pd.DataFrame, optional
            Sample metadata with batch/condition information.
        normalization : str, default='log2'
            Normalization method to apply:
            
            - 'log2': log2(counts + 1) - DEFAULT
                Simple and effective for most RNA-seq QC.
                Good when library sizes are relatively similar.
            
            - 'cpm': Counts Per Million + log2
                CPM = (counts / library_size) * 1e6, then log2(CPM + 1)
                Use when library sizes vary significantly (>5x).
                Accounts for sequencing depth differences.
            
            - 'quantile': Quantile normalization + log2
                Forces all samples to have the same distribution.
                More aggressive normalization, removes technical variation.
                Use when samples have very different distributions.
            
            - 'none': No transformation
                Use if data is already normalized (TPM, FPKM, VST, rlog).
                Assumes data is on log-scale or otherwise ready for analysis.
        
        Examples
        --------
        >>> # Standard usage with raw counts
        >>> assessor = DataQualityAssessor(counts, metadata)
        
        >>> # When library sizes vary a lot
        >>> assessor = DataQualityAssessor(counts, metadata, normalization='cpm')
        
        >>> # With pre-normalized TPM data
        >>> assessor = DataQualityAssessor(tpm_data, metadata, normalization='none')
        """
        # Validate inputs
        validate_count_matrix(counts, allow_negative=False, min_genes=10, min_samples=2)
        
        if metadata is not None:
            validate_metadata(metadata, counts)
        
        # Validate normalization method
        if normalization not in self.NORMALIZATION_METHODS:
            raise ValueError(
                f"Invalid normalization method: '{normalization}'. "
                f"Valid options: {self.NORMALIZATION_METHODS}"
            )
        
        self.counts = counts
        self.metadata = metadata
        self.normalization = normalization
        self.n_genes, self.n_samples = counts.shape
        self.sample_names = counts.columns.tolist()
        self.feature_names = counts.index.tolist()
        
        # Prepare normalized data for analyses
        self._prepare_data()
        
        # Results storage
        self.results = {}
        
        logger.info(
            f"Initialized quality assessor for {self.n_samples} samples, "
            f"{self.n_genes} genes (normalization: {normalization})"
        )
    
    def _prepare_data(self):
        """
        Prepare normalized and scaled data matrices.
        
        Applies the selected normalization method and scales data for
        PCA and other multivariate analyses.
        """
        if self.normalization == 'none':
            # User provides pre-normalized data (TPM, FPKM, VST, etc.)
            # Assume it's already on appropriate scale
            self.log_counts = self.counts.copy()
            logger.info("Using pre-normalized data (no transformation applied)")
            
        elif self.normalization == 'cpm':
            # Counts Per Million + log2
            # Accounts for library size differences
            lib_sizes = self.counts.sum(axis=0)
            cpm = self.counts.div(lib_sizes, axis=1) * 1e6
            self.log_counts = np.log2(cpm + 1)
            logger.info(
                f"Applied CPM normalization (library sizes: "
                f"{lib_sizes.min():,.0f} - {lib_sizes.max():,.0f})"
            )
            
        elif self.normalization == 'quantile':
            # Quantile normalization + log2
            # Forces same distribution across all samples
            self.log_counts = self._quantile_normalize(self.counts)
            logger.info("Applied quantile normalization")
            
        else:  # 'log2' (default)
            # Simple log2 transformation
            self.log_counts = np.log2(self.counts + 1)
            logger.info("Applied log2(counts + 1) transformation")
        
        # Scaled data (samples x features) for PCA/outlier detection
        self.scaler = StandardScaler()
        self.data_scaled = self.scaler.fit_transform(self.log_counts.T)
    
    def _quantile_normalize(self, counts: pd.DataFrame) -> pd.DataFrame:
        """
        Perform quantile normalization.
        
        Forces all samples to have the same distribution by:
        1. Ranking values within each sample
        2. Replacing with mean value at each rank across samples
        
        Parameters
        ----------
        counts : pd.DataFrame
            Raw count matrix (genes x samples)
        
        Returns
        -------
        pd.DataFrame
            Quantile normalized and log2 transformed data
        """
        # Log transform first to handle zeros and skewness
        log_data = np.log2(counts + 1)
        
        # Sort each column and get mean for each rank position
        sorted_data = np.sort(log_data.values, axis=0)
        rank_means = sorted_data.mean(axis=1)
        
        # Get ranks for each column (handling ties with average)
        ranks = np.zeros_like(log_data.values)
        for i in range(log_data.shape[1]):
            ranks[:, i] = stats.rankdata(log_data.values[:, i], method='ordinal') - 1
        
        # Replace values with rank means
        normalized = np.zeros_like(log_data.values)
        for i in range(log_data.shape[1]):
            normalized[:, i] = rank_means[ranks[:, i].astype(int)]
        
        return pd.DataFrame(
            normalized,
            index=counts.index,
            columns=counts.columns
        )
    
    # =========================================================================
    # MAIN QUALITY ASSESSMENT
    # =========================================================================
    
    @handle_errors(exit_on_error=False)
    def assess_quality(self) -> Dict:
        """
        Perform comprehensive quality assessment.
        
        Returns
        -------
        dict
            Complete quality report with all metrics
        """
        logger.info("Starting comprehensive quality assessment...")
        
        # Run all assessments
        self._assess_library_quality()
        self._assess_gene_detection()
        self._detect_outliers()
        self._assess_variance_structure()
        self._detect_batch_effects()
        self._assess_biological_signal()
        self._calculate_overall_score()
        
        # Compile report
        report = self._compile_report()
        
        logger.info("Quality assessment complete")
        return report
    
    def _assess_library_quality(self):
        """Assess library size and distribution quality."""
        logger.info("Assessing library quality...")
        
        lib_sizes = self.counts.sum(axis=0)
        
        mean_size = lib_sizes.mean()
        cv = lib_sizes.std() / mean_size if mean_size > 0 else 0
        min_size = lib_sizes.min()
        max_size = lib_sizes.max()
        
        flags = []
        if cv > 0.5:
            flags.append("High library size variation")
        if min_size < mean_size * 0.3:
            flags.append("Some libraries very small")
        if max_size > mean_size * 3:
            flags.append("Some libraries very large")
        
        score = 100
        score -= min(30, cv * 50)
        score -= min(20, (max_size / min_size - 1) * 10) if min_size > 0 else 20
        
        self.results['library_quality'] = {
            'mean_size': float(mean_size),
            'cv': float(cv),
            'min_size': float(min_size),
            'max_size': float(max_size),
            'size_range': float(max_size / min_size) if min_size > 0 else float('inf'),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _assess_gene_detection(self):
        """Assess gene detection and expression distribution."""
        logger.info("Assessing gene detection...")
        
        zero_counts = (self.counts == 0).sum(axis=1)
        detection_rate = 1 - (zero_counts / self.n_samples)
        
        mean_expr = self.counts.mean(axis=1)
        
        n_high = (mean_expr > mean_expr.quantile(0.75)).sum()
        n_medium = ((mean_expr > mean_expr.quantile(0.25)) & 
                    (mean_expr <= mean_expr.quantile(0.75))).sum()
        n_low = (mean_expr <= mean_expr.quantile(0.25)).sum()
        
        flags = []
        zero_pct = (zero_counts.mean() / self.n_samples) * 100
        if zero_pct > 70:
            flags.append("Very high zero inflation")
        elif zero_pct > 50:
            flags.append("High zero inflation")
        
        score = 100
        score -= min(40, zero_pct * 0.5)
        if n_high < self.n_genes * 0.05:
            score -= 10
        
        self.results['gene_detection'] = {
            'mean_detection_rate': float(detection_rate.mean()),
            'zero_inflation_pct': float(zero_pct),
            'n_highly_expressed': int(n_high),
            'n_medium_expressed': int(n_medium),
            'n_low_expressed': int(n_low),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _detect_outliers(self):
        """Basic outlier detection (for backward compatibility)."""
        logger.info("Detecting outlier samples (basic)...")
        
        # PCA-based outlier detection
        pca = PCA(n_components=min(5, self.n_samples - 1))
        pca_scores = pca.fit_transform(self.data_scaled)
        
        # Mahalanobis distance
        mean = pca_scores.mean(axis=0)
        cov = np.cov(pca_scores.T)
        inv_cov = np.linalg.pinv(cov)
        
        mahal_dist = []
        for i in range(len(pca_scores)):
            diff = pca_scores[i] - mean
            dist = np.sqrt(diff @ inv_cov @ diff.T)
            mahal_dist.append(dist)
        
        mahal_dist = np.array(mahal_dist)
        threshold = mahal_dist.mean() + 3 * mahal_dist.std()
        outliers_pca = mahal_dist > threshold
        
        # Library size outliers
        lib_sizes = self.counts.sum(axis=0)
        z_scores = np.abs(stats.zscore(lib_sizes))
        outliers_libsize = z_scores > 3
        
        # Combine
        outlier_scores = outliers_pca.astype(int) + outliers_libsize.astype(int)
        outlier_samples = self.counts.columns[outlier_scores >= 2].tolist()
        
        n_outliers = len(outlier_samples)
        outlier_pct = (n_outliers / self.n_samples) * 100
        score = 100 - min(50, outlier_pct * 5)
        
        flags = []
        if n_outliers > 0:
            flags.append(f"{n_outliers} potential outlier sample(s) detected")
        
        self.results['outlier_detection'] = {
            'n_outliers': n_outliers,
            'outlier_percentage': float(outlier_pct),
            'outlier_samples': outlier_samples,
            'mahalanobis_distances': mahal_dist.tolist(),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if n_outliers == 0 else 'warning' if n_outliers <= 2 else 'poor'
        }
    
    def _assess_variance_structure(self):
        """Assess variance structure and components."""
        logger.info("Assessing variance structure...")
        
        pca = PCA()
        pca.fit(self.data_scaled)
        
        variance_explained = pca.explained_variance_ratio_
        pc1_var = variance_explained[0] * 100
        pc2_var = variance_explained[1] * 100 if len(variance_explained) > 1 else 0
        
        cum_var = np.cumsum(variance_explained) * 100
        n_components_90 = np.where(cum_var >= 90)[0][0] + 1 if any(cum_var >= 90) else len(variance_explained)
        
        flags = []
        if pc1_var > 50:
            flags.append("First PC explains >50% variance (potential batch effect)")
        if n_components_90 < 3:
            flags.append("Very few components needed for 90% variance")
        
        score = 100
        if pc1_var > 50:
            score -= 20
        if pc1_var > 70:
            score -= 20
        if n_components_90 < 2:
            score -= 15
        
        self.results['variance_structure'] = {
            'pc1_variance': float(pc1_var),
            'pc2_variance': float(pc2_var),
            'variance_explained_top5': variance_explained[:5].tolist(),
            'n_components_90pct': int(n_components_90),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _detect_batch_effects(self):
        """
        Detect potential batch effects with improved logic.
        
        Strategy:
        1. Identify likely batch columns (batch, plate, run, date, etc.)
        2. Identify likely condition columns (condition, treatment, group, etc.)
        3. Calculate F-statistic for each on PC1
        4. Compare batch vs condition effects
        5. Check for confounding (batch correlated with condition)
        """
        logger.info("Detecting batch effects...")
        
        pca = PCA(n_components=min(3, self.n_samples - 1))
        pca_coords = pca.fit_transform(self.data_scaled)
        variance_explained = pca.explained_variance_ratio_
        
        batch_detected = False
        batch_strength = 0
        condition_strength = 0
        flags = []
        batch_variable = None
        condition_variable = None
        confounded = False
        
        # Keywords to identify column types
        batch_keywords = ['batch', 'plate', 'run', 'lane', 'date', 'flow', 'chip', 
                          'slide', 'array', 'seq_run', 'library_prep', 'extraction']
        condition_keywords = ['condition', 'treatment', 'group', 'status', 'disease',
                              'genotype', 'phenotype', 'sample_type', 'cell_type', 'tissue']
        
        if self.metadata is not None and len(self.metadata) > 0:
            categorical_cols = self.metadata.select_dtypes(include=['object', 'string', 'category']).columns.tolist()
            
            # Separate columns into batch-like and condition-like
            batch_cols = []
            condition_cols = []
            other_cols = []
            
            for col in categorical_cols:
                col_lower = col.lower()
                if any(kw in col_lower for kw in batch_keywords):
                    batch_cols.append(col)
                elif any(kw in col_lower for kw in condition_keywords):
                    condition_cols.append(col)
                else:
                    other_cols.append(col)
            
            # Calculate F-statistics for all columns
            all_f_stats = {}
            
            for col in categorical_cols:
                groups = self.metadata[col].values
                unique_groups = np.unique(groups)
                
                if 1 < len(unique_groups) < self.n_samples:
                    f_stat = self._calculate_f_statistic(pca_coords[:, 0], groups)
                    if f_stat is not None:
                        all_f_stats[col] = f_stat
            
            # Find strongest batch effect
            batch_f_stats = {k: v for k, v in all_f_stats.items() if k in batch_cols}
            if batch_f_stats:
                batch_variable = max(batch_f_stats.items(), key=lambda x: x[1])[0]
                batch_strength = batch_f_stats[batch_variable]
            
            # Find strongest condition effect
            condition_f_stats = {k: v for k, v in all_f_stats.items() if k in condition_cols}
            if condition_f_stats:
                condition_variable = max(condition_f_stats.items(), key=lambda x: x[1])[0]
                condition_strength = condition_f_stats[condition_variable]
            
            # If no explicit batch/condition columns found, use other columns
            if not batch_variable and not condition_variable and other_cols:
                # Assume first column is condition, second is batch (if exists)
                if other_cols:
                    condition_variable = other_cols[0]
                    condition_strength = all_f_stats.get(condition_variable, 0)
                if len(other_cols) > 1:
                    batch_variable = other_cols[1]
                    batch_strength = all_f_stats.get(batch_variable, 0)
            
            # Check for confounding
            if batch_variable and condition_variable:
                confounded = self._check_confounding(
                    self.metadata[batch_variable].values,
                    self.metadata[condition_variable].values
                )
            
            # Determine if batch effect is problematic
            # Batch is problematic if: batch F-stat > condition F-stat AND batch F-stat > threshold
            if batch_strength > 5:
                if condition_strength > 0 and batch_strength > condition_strength:
                    batch_detected = True
                    flags.append(f"Batch effect in '{batch_variable}' (F={batch_strength:.1f}) exceeds condition effect (F={condition_strength:.1f})")
                elif condition_strength == 0:
                    batch_detected = True
                    flags.append(f"Strong batch effect detected in '{batch_variable}' (F={batch_strength:.1f})")
            elif batch_strength > 2:
                flags.append(f"Moderate batch effect possible in '{batch_variable}' (F={batch_strength:.1f})")
            
            # Add condition info
            if condition_strength > 2:
                flags.append(f"Biological signal detected in '{condition_variable}' (F={condition_strength:.1f})")
            
            # Add confounding warning
            if confounded:
                flags.append(f"⚠ WARNING: '{batch_variable}' and '{condition_variable}' appear confounded!")
                flags.append("  Cannot separate batch effect from biological effect")
        else:
            flags.append("No metadata provided - cannot check for batch effects")
        
        # Calculate score
        score = 100
        if batch_detected:
            score -= min(40, batch_strength * 3)
        if confounded:
            score -= 20
        
        self.results['batch_effects'] = {
            'batch_detected': batch_detected,
            'batch_variable': batch_variable,
            'batch_strength': float(batch_strength) if batch_strength > 0 else None,
            'condition_variable': condition_variable,
            'condition_strength': float(condition_strength) if condition_strength > 0 else None,
            'confounded': confounded,
            'pca_coordinates': pca_coords.tolist(),
            'variance_explained_pc1': float(variance_explained[0]),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 80 else 'warning' if score > 60 else 'poor',
            'recommendation': self._get_batch_recommendation(batch_detected, batch_strength, confounded)
        }
    
    def _calculate_f_statistic(self, data: np.ndarray, groups: np.ndarray) -> Optional[float]:
        """Calculate F-statistic for ANOVA."""
        try:
            unique_groups = np.unique(groups)
            if len(unique_groups) < 2:
                return None
            
            between_var = 0
            within_var = 0
            overall_mean = data.mean()
            
            for group in unique_groups:
                group_mask = groups == group
                group_data = data[group_mask]
                
                if len(group_data) > 0:
                    group_mean = group_data.mean()
                    between_var += len(group_data) * (group_mean - overall_mean) ** 2
                    within_var += ((group_data - group_mean) ** 2).sum()
            
            df_between = len(unique_groups) - 1
            df_within = len(data) - len(unique_groups)
            
            if within_var > 0 and df_within > 0:
                f_stat = (between_var / df_between) / (within_var / df_within)
                return float(f_stat)
            return None
        except:
            return None
    
    def _check_confounding(self, batch: np.ndarray, condition: np.ndarray) -> bool:
        """Check if batch and condition are confounded."""
        try:
            # Create contingency table
            unique_batch = np.unique(batch)
            unique_condition = np.unique(condition)
            
            if len(unique_batch) < 2 or len(unique_condition) < 2:
                return False
            
            # Check if each condition is in only one batch (perfect confounding)
            for cond in unique_condition:
                cond_mask = condition == cond
                batches_for_cond = np.unique(batch[cond_mask])
                if len(batches_for_cond) == 1:
                    # This condition is only in one batch - potential confounding
                    return True
            
            return False
        except:
            return False
    
    def _get_batch_recommendation(self, detected: bool, strength: float, confounded: bool = False) -> str:
        """Get recommendation for batch effect handling."""
        if confounded:
            return ("CONFOUNDED DESIGN: Batch and condition cannot be separated. "
                    "Results may be unreliable. Consider repeating experiment with "
                    "balanced design (each condition in multiple batches).")
        if not detected:
            return "No significant batch effect detected. Proceed with standard analysis."
        if strength > 10:
            return ("Strong batch effect detected. Recommended actions:\n"
                    "  1. Include batch as covariate in DE model: ~batch + condition\n"
                    "  2. For visualization: Use ComBat or limma::removeBatchEffect\n"
                    "  3. Do NOT use batch-corrected data for DE analysis")
        elif strength > 5:
            return ("Moderate batch effect detected. Recommended actions:\n"
                    "  1. Include batch as covariate in DE model: ~batch + condition\n"
                    "  2. Monitor batch separation in PCA plots")
        else:
            return "Weak batch effect. Include batch as covariate if concerned."
    
    def _assess_biological_signal(self):
        """Assess strength of biological signal."""
        logger.info("Assessing biological signal strength...")
        
        means = self.counts.mean(axis=1)
        stds = self.counts.std(axis=1)
        
        cv = np.zeros(len(means))
        nonzero = means > 0
        cv[nonzero] = stds[nonzero] / means[nonzero]
        
        median_cv = np.median(cv[cv > 0]) if (cv > 0).any() else 0
        
        high_mean = means > means.quantile(0.5)
        high_var = cv > median_cv
        signal_genes = high_mean & high_var
        n_signal = signal_genes.sum()
        
        signal_strength = (n_signal / self.n_genes) * 100
        
        flags = []
        if signal_strength < 10:
            flags.append("Low biological signal detected")
        elif signal_strength > 30:
            flags.append("Strong biological signal detected")
        
        score = min(100, signal_strength * 2)
        
        self.results['biological_signal'] = {
            'signal_strength': float(signal_strength),
            'n_signal_genes': int(n_signal),
            'median_cv': float(median_cv),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _calculate_overall_score(self):
        """Calculate weighted overall quality score."""
        weights = {
            'library_quality': 0.15,
            'gene_detection': 0.20,
            'outlier_detection': 0.15,
            'variance_structure': 0.15,
            'batch_effects': 0.20,
            'biological_signal': 0.15
        }
        
        total_score = 0
        for component, weight in weights.items():
            if component in self.results:
                total_score += self.results[component]['score'] * weight
        
        # Determine status
        if total_score >= 80:
            status = 'good'
            recommendation = "Data quality is good. Proceed with analysis."
        elif total_score >= 60:
            status = 'acceptable'
            recommendation = "Data quality is acceptable. Address flagged issues."
        else:
            status = 'poor'
            recommendation = "Data quality is poor. Investigate issues before analysis."
        
        self.results['overall'] = {
            'score': float(total_score),
            'status': status,
            'recommendation': recommendation
        }
    
    def _compile_report(self) -> Dict:
        """Compile all results into final report."""
        return {
            'overall': self.results.get('overall', {}),
            'components': {
                'library_quality': self.results.get('library_quality', {}),
                'gene_detection': self.results.get('gene_detection', {}),
                'outlier_detection': self.results.get('outlier_detection', {}),
                'variance_structure': self.results.get('variance_structure', {}),
                'batch_effects': self.results.get('batch_effects', {}),
                'biological_signal': self.results.get('biological_signal', {})
            },
            'summary': self._generate_summary()
        }
    
    def _generate_summary(self) -> str:
        """Generate text summary of quality assessment."""
        lines = [
            f"Overall Quality Score: {self.results['overall']['score']:.1f}/100",
            f"Status: {self.results['overall']['status'].upper()}",
            f"Recommendation: {self.results['overall']['recommendation']}",
            "",
            "Component Scores:"
        ]
        
        for name in ['library_quality', 'gene_detection', 'outlier_detection',
                     'variance_structure', 'batch_effects', 'biological_signal']:
            if name in self.results:
                score = self.results[name]['score']
                status = self.results[name]['status']
                icon = "✓" if status == 'good' else "⚠" if status == 'warning' else "✗"
                lines.append(f"  {icon} {name.replace('_', ' ').title()}: {score:.1f}/100")
        
        return "\n".join(lines)
    
    # =========================================================================
    # ADVANCED OUTLIER DETECTION
    # =========================================================================
    
    def detect_outliers_advanced(
        self,
        methods: List[str] = None,
        consensus_threshold: int = 3,
        contamination: float = 0.1
    ) -> OutlierResult:
        """
        Advanced outlier detection using multiple methods.
        
        NEW in v2.2.0: Uses consensus of 6 different methods for robust
        outlier identification.
        
        Parameters
        ----------
        methods : list of str, optional
            Methods to use. Options:
            - 'pca_mahalanobis': PCA + Mahalanobis distance
            - 'isolation_forest': Isolation Forest algorithm
            - 'lof': Local Outlier Factor
            - 'elliptic_envelope': Robust covariance
            - 'correlation': Low correlation with other samples
            - 'library_size': Z-score of library sizes
            Default: all methods
        consensus_threshold : int
            Minimum number of methods that must flag a sample
        contamination : float
            Expected proportion of outliers (0-0.5)
        
        Returns
        -------
        OutlierResult
            Results with outlier samples, scores, and method details
        
        Examples
        --------
        >>> result = assessor.detect_outliers_advanced(consensus_threshold=3)
        >>> print(result.summary())
        >>> print(f"Outliers: {result.outlier_samples}")
        """
        if methods is None:
            methods = ['pca_mahalanobis', 'isolation_forest', 'lof',
                      'elliptic_envelope', 'correlation', 'library_size']
        
        logger.info(f"Running advanced outlier detection with {len(methods)} methods...")
        
        method_results = {}
        outlier_votes = {sample: 0 for sample in self.sample_names}
        
        # Run each method
        if 'pca_mahalanobis' in methods:
            outliers = self._outlier_pca_mahalanobis()
            method_results['PCA + Mahalanobis'] = outliers
            for s in outliers:
                outlier_votes[s] += 1
        
        if 'isolation_forest' in methods:
            outliers = self._outlier_isolation_forest(contamination)
            method_results['Isolation Forest'] = outliers
            for s in outliers:
                outlier_votes[s] += 1
        
        if 'lof' in methods:
            outliers = self._outlier_lof(contamination)
            method_results['Local Outlier Factor'] = outliers
            for s in outliers:
                outlier_votes[s] += 1
        
        if 'elliptic_envelope' in methods:
            outliers = self._outlier_elliptic_envelope(contamination)
            method_results['Elliptic Envelope'] = outliers
            for s in outliers:
                outlier_votes[s] += 1
        
        if 'correlation' in methods:
            outliers = self._outlier_correlation()
            method_results['Correlation-based'] = outliers
            for s in outliers:
                outlier_votes[s] += 1
        
        if 'library_size' in methods:
            outliers = self._outlier_library_size()
            method_results['Library Size'] = outliers
            for s in outliers:
                outlier_votes[s] += 1
        
        # Consensus outliers
        consensus_outliers = [
            sample for sample, votes in outlier_votes.items()
            if votes >= consensus_threshold
        ]
        
        result = OutlierResult(
            outlier_samples=consensus_outliers,
            outlier_scores=outlier_votes,
            method_results=method_results,
            consensus_threshold=consensus_threshold,
            n_outliers=len(consensus_outliers),
            outlier_percentage=(len(consensus_outliers) / self.n_samples) * 100
        )
        
        logger.info(f"Advanced outlier detection complete: {result.n_outliers} consensus outliers")
        return result
    
    def _outlier_pca_mahalanobis(self, n_components: int = 5) -> List[str]:
        """PCA-based outlier detection using Mahalanobis distance."""
        n_comp = min(n_components, self.n_samples - 1, self.n_genes)
        pca = PCA(n_components=n_comp)
        scores = pca.fit_transform(self.data_scaled)
        
        mean = scores.mean(axis=0)
        cov = np.cov(scores.T)
        inv_cov = np.linalg.pinv(cov)
        
        distances = []
        for i in range(len(scores)):
            diff = scores[i] - mean
            dist = np.sqrt(diff @ inv_cov @ diff.T)
            distances.append(dist)
        
        distances = np.array(distances)
        threshold = np.percentile(distances, 95)
        
        return [self.sample_names[i] for i, d in enumerate(distances) if d > threshold]
    
    def _outlier_isolation_forest(self, contamination: float) -> List[str]:
        """Isolation Forest outlier detection."""
        clf = IsolationForest(contamination=contamination, random_state=42, n_jobs=-1)
        predictions = clf.fit_predict(self.data_scaled)
        return [self.sample_names[i] for i, p in enumerate(predictions) if p == -1]
    
    def _outlier_lof(self, contamination: float) -> List[str]:
        """Local Outlier Factor detection."""
        n_neighbors = min(20, self.n_samples - 1)
        clf = LocalOutlierFactor(n_neighbors=n_neighbors, contamination=contamination)
        predictions = clf.fit_predict(self.data_scaled)
        return [self.sample_names[i] for i, p in enumerate(predictions) if p == -1]
    
    def _outlier_elliptic_envelope(self, contamination: float) -> List[str]:
        """Robust covariance-based outlier detection."""
        try:
            # Reduce dimensionality if needed
            if self.n_genes > self.n_samples:
                pca = PCA(n_components=min(10, self.n_samples - 1))
                data_reduced = pca.fit_transform(self.data_scaled)
            else:
                data_reduced = self.data_scaled
            
            clf = EllipticEnvelope(contamination=contamination, random_state=42)
            predictions = clf.fit_predict(data_reduced)
            return [self.sample_names[i] for i, p in enumerate(predictions) if p == -1]
        except Exception as e:
            logger.warning(f"Elliptic envelope failed: {e}")
            return []
    
    def _outlier_correlation(self) -> List[str]:
        """Detect samples with low average correlation to others."""
        corr_matrix = np.corrcoef(self.data_scaled)
        
        avg_corr = []
        for i in range(self.n_samples):
            other_corrs = [corr_matrix[i, j] for j in range(self.n_samples) if i != j]
            avg_corr.append(np.mean(other_corrs))
        
        avg_corr = np.array(avg_corr)
        median_corr = np.median(avg_corr)
        mad = np.median(np.abs(avg_corr - median_corr))
        threshold = median_corr - 3 * 1.4826 * mad
        
        return [self.sample_names[i] for i, c in enumerate(avg_corr) if c < threshold]
    
    def _outlier_library_size(self) -> List[str]:
        """Detect outliers based on library size."""
        lib_sizes = self.counts.sum(axis=0)
        z_scores = np.abs(stats.zscore(lib_sizes))
        return [self.sample_names[i] for i, z in enumerate(z_scores) if z > 3]
    
    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    
    def plot_quality_report(self, output_file: str = None):
        """Generate comprehensive quality visualization."""
        import matplotlib.pyplot as plt
        
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 1. Overall score gauge
        ax1 = fig.add_subplot(gs[0, 0])
        self._plot_score_gauge(ax1, self.results['overall']['score'], 'Overall Quality')
        
        # 2. Component scores
        ax2 = fig.add_subplot(gs[0, 1:])
        self._plot_component_scores(ax2)
        
        # 3. Library size distribution
        ax3 = fig.add_subplot(gs[1, 0])
        lib_sizes = self.counts.sum(axis=0)
        ax3.hist(lib_sizes / 1e6, bins=20, color='steelblue', edgecolor='black')
        ax3.set_xlabel('Library Size (millions)')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Library Size Distribution', fontweight='bold')
        ax3.axvline(lib_sizes.mean() / 1e6, color='red', linestyle='--', label='Mean')
        ax3.legend()
        
        # 4. PCA plot
        ax4 = fig.add_subplot(gs[1, 1])
        pca_coords = np.array(self.results['batch_effects']['pca_coordinates'])
        
        if self.results['outlier_detection']['outlier_samples']:
            outliers = self.results['outlier_detection']['outlier_samples']
            outlier_mask = np.array([s in outliers for s in self.sample_names])
            
            ax4.scatter(pca_coords[~outlier_mask, 0], pca_coords[~outlier_mask, 1],
                       c='steelblue', s=50, alpha=0.6, label='Normal')
            ax4.scatter(pca_coords[outlier_mask, 0], pca_coords[outlier_mask, 1],
                       c='red', s=50, alpha=0.8, label='Outlier')
            ax4.legend()
        else:
            ax4.scatter(pca_coords[:, 0], pca_coords[:, 1], c='steelblue', s=50, alpha=0.6)
        
        ax4.set_xlabel('PC1')
        ax4.set_ylabel('PC2')
        ax4.set_title('PCA Plot', fontweight='bold')
        
        # 5. Variance explained
        ax5 = fig.add_subplot(gs[1, 2])
        var_exp = self.results['variance_structure']['variance_explained_top5']
        ax5.bar(range(1, len(var_exp) + 1), np.array(var_exp) * 100, color='coral')
        ax5.set_xlabel('Principal Component')
        ax5.set_ylabel('Variance Explained (%)')
        ax5.set_title('PCA Variance', fontweight='bold')
        
        # 6. Gene detection
        ax6 = fig.add_subplot(gs[2, 0])
        detection = self.results['gene_detection']
        categories = ['High', 'Medium', 'Low']
        values = [detection['n_highly_expressed'],
                 detection['n_medium_expressed'],
                 detection['n_low_expressed']]
        colors = ['darkgreen', 'orange', 'lightcoral']
        ax6.pie(values, labels=categories, colors=colors, autopct='%1.1f%%', startangle=90)
        ax6.set_title('Gene Expression Levels', fontweight='bold')
        
        # 7. Quality summary
        ax7 = fig.add_subplot(gs[2, 1:])
        ax7.axis('off')
        
        summary_text = self._generate_summary()
        ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes,
                fontsize=9, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.suptitle('RAPTOR Data Quality Assessment Report v2.2.0',
                    fontsize=16, fontweight='bold', y=0.98)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Quality report saved to {output_file}")
        
        return fig
    
    def plot_outliers_advanced(self, result: OutlierResult, output_file: str = None):
        """
        Visualize advanced outlier detection results.
        
        NEW in v2.2.0
        """
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # 1. PCA with outliers
        ax1 = axes[0, 0]
        pca = PCA(n_components=2)
        scores = pca.fit_transform(self.data_scaled)
        
        outlier_mask = np.array([s in result.outlier_samples for s in self.sample_names])
        
        ax1.scatter(scores[~outlier_mask, 0], scores[~outlier_mask, 1],
                   c='steelblue', s=80, alpha=0.7, label='Normal', edgecolors='black')
        ax1.scatter(scores[outlier_mask, 0], scores[outlier_mask, 1],
                   c='red', s=120, alpha=0.9, label='Outlier', marker='X', edgecolors='black')
        
        for i, sample in enumerate(self.sample_names):
            if sample in result.outlier_samples:
                ax1.annotate(sample, (scores[i, 0], scores[i, 1]), fontsize=8,
                           xytext=(5, 5), textcoords='offset points')
        
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
        ax1.set_title('PCA with Outliers Highlighted', fontweight='bold')
        ax1.legend()
        
        # 2. Outlier votes
        ax2 = axes[0, 1]
        votes = list(result.outlier_scores.values())
        colors = ['red' if v >= result.consensus_threshold else 'steelblue' for v in votes]
        
        ax2.bar(range(len(votes)), votes, color=colors, edgecolor='black')
        ax2.axhline(result.consensus_threshold, color='black', linestyle='--',
                   label=f'Threshold ({result.consensus_threshold})')
        ax2.set_xlabel('Sample Index')
        ax2.set_ylabel('Methods Flagging as Outlier')
        ax2.set_title('Outlier Votes per Sample', fontweight='bold')
        ax2.legend()
        
        # 3. Method comparison
        ax3 = axes[1, 0]
        method_names = list(result.method_results.keys())
        method_counts = [len(v) for v in result.method_results.values()]
        
        ax3.barh(method_names, method_counts, color='coral', edgecolor='black')
        ax3.set_xlabel('Number of Outliers')
        ax3.set_title('Outliers by Method', fontweight='bold')
        
        # 4. Library size
        ax4 = axes[1, 1]
        lib_sizes = self.counts.sum(axis=0) / 1e6
        
        ax4.hist(lib_sizes[~outlier_mask], bins=15, alpha=0.7, label='Normal',
                color='steelblue', edgecolor='black')
        if outlier_mask.any():
            ax4.hist(lib_sizes[outlier_mask], bins=5, alpha=0.8, label='Outliers',
                    color='red', edgecolor='black')
        ax4.set_xlabel('Library Size (millions)')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Library Size Distribution', fontweight='bold')
        ax4.legend()
        
        plt.suptitle('RAPTOR Advanced Outlier Detection', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            logger.info(f"Outlier plot saved to {output_file}")
        
        return fig
    
    def _plot_score_gauge(self, ax, score: float, title: str):
        """Plot score as gauge chart."""
        if score >= 80:
            color = 'green'
        elif score >= 60:
            color = 'orange'
        else:
            color = 'red'
        
        theta = np.linspace(0, np.pi, 100)
        r = np.ones(100)
        
        ax.plot(theta, r, 'k-', linewidth=2)
        ax.fill_between(theta, 0, r, alpha=0.1, color='gray')
        
        score_theta = (score / 100) * np.pi
        ax.plot([0, score_theta], [0, 1], color=color, linewidth=4)
        ax.plot(score_theta, 1, 'o', color=color, markersize=15)
        
        ax.text(np.pi/2, 0.5, f'{score:.1f}', ha='center', va='center',
               fontsize=24, fontweight='bold')
        ax.text(np.pi/2, 0.2, title, ha='center', va='center',
               fontsize=12, fontweight='bold')
        
        ax.set_xlim([0, np.pi])
        ax.set_ylim([0, 1.2])
        ax.axis('off')
    
    def _plot_component_scores(self, ax):
        """Plot component scores as horizontal bar chart."""
        components = ['Library\nQuality', 'Gene\nDetection', 'Outlier\nDetection',
                     'Variance\nStructure', 'Batch\nEffects', 'Biological\nSignal']
        
        scores = [
            self.results['library_quality']['score'],
            self.results['gene_detection']['score'],
            self.results['outlier_detection']['score'],
            self.results['variance_structure']['score'],
            self.results['batch_effects']['score'],
            self.results['biological_signal']['score']
        ]
        
        colors = ['green' if s >= 80 else 'orange' if s >= 60 else 'red' for s in scores]
        
        y_pos = np.arange(len(components))
        bars = ax.barh(y_pos, scores, color=colors, alpha=0.7, edgecolor='black')
        
        for i, (bar, score) in enumerate(zip(bars, scores)):
            ax.text(score + 2, i, f'{score:.0f}', va='center', fontweight='bold')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(components, fontsize=9)
        ax.set_xlabel('Score')
        ax.set_title('Component Scores', fontweight='bold')
        ax.set_xlim([0, 105])
        ax.axvline(80, color='green', linestyle='--', alpha=0.3)
        ax.axvline(60, color='orange', linestyle='--', alpha=0.3)
        ax.grid(axis='x', alpha=0.3)


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def quick_quality_check(
    counts: pd.DataFrame,
    metadata: pd.DataFrame = None,
    normalization: str = 'log2',
    plot: bool = True,
    output_file: str = 'quality_report.png'
) -> Dict:
    """
    Quick quality assessment with visualization.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes x samples)
    metadata : pd.DataFrame, optional
        Sample metadata
    normalization : str, default='log2'
        Normalization method: 'log2', 'cpm', 'quantile', or 'none'
    plot : bool
        Generate visualization
    output_file : str
        Output file for plot
    
    Returns
    -------
    dict
        Quality report
    
    Examples
    --------
    >>> report = quick_quality_check(counts, metadata)
    >>> print(f"Quality Score: {report['overall']['score']:.1f}")
    
    >>> # With CPM normalization for varying library sizes
    >>> report = quick_quality_check(counts, metadata, normalization='cpm')
    """
    assessor = DataQualityAssessor(counts, metadata, normalization=normalization)
    report = assessor.assess_quality()
    
    if plot:
        assessor.plot_quality_report(output_file)
    
    print("\n" + "=" * 70)
    print("🦖 RAPTOR DATA QUALITY ASSESSMENT")
    print("=" * 70)
    print(f"Normalization: {normalization}")
    print(report['summary'])
    print("=" * 70 + "\n")
    
    return report


def detect_outliers_quick(
    counts: pd.DataFrame,
    consensus_threshold: int = 3,
    normalization: str = 'log2',
    plot: bool = True,
    output_file: str = 'outliers.png'
) -> OutlierResult:
    """
    Quick advanced outlier detection.
    
    NEW in v2.2.0
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix
    consensus_threshold : int
        Methods that must agree
    normalization : str, default='log2'
        Normalization method: 'log2', 'cpm', 'quantile', or 'none'
    plot : bool
        Generate visualization
    output_file : str
        Output file
    
    Returns
    -------
    OutlierResult
        Outlier detection results
    
    Examples
    --------
    >>> result = detect_outliers_quick(counts, consensus_threshold=3)
    >>> print(f"Outliers: {result.outlier_samples}")
    """
    assessor = DataQualityAssessor(counts, normalization=normalization)
    result = assessor.detect_outliers_advanced(consensus_threshold=consensus_threshold)
    
    print(result.summary())
    
    if plot:
        assessor.plot_outliers_advanced(result, output_file)
    
    return result


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("""
╔═══════════════════════════════════════════════════════════════════════╗
║     🦖 RAPTOR Advanced Data Quality Assessment Module v2.2.0          ║
╠═══════════════════════════════════════════════════════════════════════╣
║                                                                       ║
║  FEATURES:                                                            ║
║  ─────────                                                            ║
║  ✓ Library quality assessment                                         ║
║  ✓ Gene detection analysis                                            ║
║  ✓ Advanced outlier detection (6 methods)                             ║
║  ✓ Variance structure analysis                                        ║
║  ✓ Batch effect detection (improved - distinguishes batch vs cond.)   ║
║  ✓ Biological signal assessment                                       ║
║  ✓ Overall quality scoring (0-100)                                    ║
║  ✓ Normalization options (log2, CPM, quantile, none)                  ║
║  ✓ Comprehensive visualization                                        ║
║                                                                       ║
╠═══════════════════════════════════════════════════════════════════════╣
║  USAGE:                                                               ║
║  ──────                                                               ║
║  from raptor.quality_assessment import (                         ║
║      DataQualityAssessor,                                             ║
║      quick_quality_check,                                             ║
║      detect_outliers_quick                                            ║
║  )                                                                    ║
║                                                                       ║
║  # Basic quality check                                                ║
║  report = quick_quality_check(counts, metadata)                       ║
║                                                                       ║
║  # With CPM normalization (for varying library sizes)                 ║
║  report = quick_quality_check(counts, metadata, normalization='cpm')  ║
║                                                                       ║
║  # Advanced outlier detection (6 methods)                             ║
║  outliers = detect_outliers_quick(counts, consensus_threshold=3)      ║
║                                                                       ║
║  # Check batch effect recommendations                                 ║
║  if report['batch_effects']['batch_detected']:                        ║
║      print(report['batch_effects']['recommendation'])                 ║
║                                                                       ║
╚═══════════════════════════════════════════════════════════════════════╝
    """)
