"""
RNA-seq Data Profiler

Analyzes RNA-seq count data to extract statistical characteristics used for
pipeline recommendation. Calculates metrics related to library size, count
distribution, biological variation, and experimental design.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class RNAseqDataProfiler:
    """
    Profile RNA-seq count matrix to extract data characteristics.
    
    Analyzes library sizes, count distributions, biological variation,
    and experimental design to provide metrics for pipeline recommendation.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix with genes as rows and samples as columns
    metadata : pd.DataFrame, optional
        Sample metadata with at least a 'condition' column
    
    Attributes
    ----------
    counts : pd.DataFrame
        Input count matrix
    metadata : pd.DataFrame
        Sample metadata
    n_genes : int
        Number of genes
    n_samples : int
        Number of samples
    
    Examples
    --------
    >>> import pandas as pd
    >>> counts = pd.read_csv('counts.csv', index_col=0)
    >>> profiler = RNAseqDataProfiler(counts)
    >>> profile = profiler.run_full_profile()
    >>> print(profile['library_stats']['cv'])
    0.25
    """
    
    def __init__(self, counts: pd.DataFrame, metadata: Optional[pd.DataFrame] = None):
        """Initialize profiler with count data."""
        self.counts = counts
        self.metadata = metadata
        self.n_genes, self.n_samples = counts.shape
        
        logger.info(f"Initialized profiler: {self.n_genes} genes × {self.n_samples} samples")
        
        # Validate input
        self._validate_input()
    
    def _validate_input(self):
        """Validate count matrix format."""
        # Check for negative values
        if (self.counts < 0).any().any():
            raise ValueError("Count matrix contains negative values")
        
        # Check for non-numeric values
        if not np.issubdtype(self.counts.values.dtype, np.number):
            raise ValueError("Count matrix must contain only numeric values")
        
        # Warn about low sample size
        if self.n_samples < 3:
            logger.warning(f"Low sample size detected (n={self.n_samples}). "
                          "Results may be unreliable.")
        
        # Warn about very small dataset
        if self.n_genes < 100:
            logger.warning(f"Very few genes detected (n={self.n_genes}). "
                          "This may not be a full transcriptome.")
    
    def run_full_profile(self) -> Dict:
        """
        Run complete data profiling.
        
        Returns
        -------
        dict
            Profile containing all calculated metrics organized by category
        
        Examples
        --------
        >>> profiler = RNAseqDataProfiler(counts)
        >>> profile = profiler.run_full_profile()
        >>> print(profile.keys())
        dict_keys(['library_stats', 'count_distribution', 'biological_variation', 
                   'design', 'sequencing', 'summary'])
        """
        logger.info("Running full profile analysis...")
        
        profile = {
            'library_stats': self.calculate_library_stats(),
            'count_distribution': self.analyze_count_distribution(),
            'biological_variation': self.estimate_biological_variation(),
            'design': self.analyze_experimental_design(),
            'sequencing': self.assess_sequencing_quality(),
        }
        
        # Add summary
        profile['summary'] = self._generate_summary(profile)
        
        logger.info("Profile analysis complete")
        return profile
    
    # =========================================================================
    # Library Size Statistics
    # =========================================================================
    
    def calculate_library_stats(self) -> Dict:
        """
        Calculate library size statistics.
        
        Returns
        -------
        dict
            Library statistics including mean, CV, and range
        
        Notes
        -----
        Library size variation is important for normalization method selection.
        High variation (CV > 0.5) suggests need for robust normalization like
        DESeq2's median-of-ratios.
        """
        library_sizes = self.counts.sum(axis=0)
        
        mean_size = library_sizes.mean()
        std_size = library_sizes.std()
        cv = std_size / mean_size if mean_size > 0 else 0
        
        return {
            'mean': float(mean_size),
            'std': float(std_size),
            'cv': float(cv),
            'min': float(library_sizes.min()),
            'max': float(library_sizes.max()),
            'range_fold': float(library_sizes.max() / library_sizes.min()) if library_sizes.min() > 0 else np.inf,
            'sizes': library_sizes.tolist()
        }
    
    # =========================================================================
    # Count Distribution Analysis
    # =========================================================================
    
    def analyze_count_distribution(self) -> Dict:
        """
        Analyze count distribution characteristics.
        
        Returns
        -------
        dict
            Count distribution metrics including zero-inflation and low counts
        
        Notes
        -----
        High zero-inflation (>60%) indicates challenging data that may benefit
        from methods designed for zero-inflated data or more stringent filtering.
        """
        # Zero-inflation
        total_counts = self.counts.size
        zero_counts = (self.counts == 0).sum().sum()
        zero_pct = (zero_counts / total_counts) * 100
        
        # Low-count genes
        gene_means = self.counts.mean(axis=1)
        low_count_genes = (gene_means < 10).sum()
        low_count_pct = (low_count_genes / self.n_genes) * 100
        
        # Expression range
        nonzero_counts = self.counts.values[self.counts.values > 0]
        if len(nonzero_counts) > 0:
            log2_range = np.log2(nonzero_counts.max()) - np.log2(nonzero_counts.min() + 1)
        else:
            log2_range = 0
        
        # Count distribution shape
        all_counts = self.counts.values.flatten()
        all_counts = all_counts[all_counts > 0]  # Exclude zeros
        if len(all_counts) > 0:
            skewness = float(stats.skew(np.log1p(all_counts)))
            kurtosis = float(stats.kurtosis(np.log1p(all_counts)))
        else:
            skewness = 0
            kurtosis = 0
        
        return {
            'zero_pct': float(zero_pct),
            'zero_count': int(zero_counts),
            'low_count_genes': int(low_count_genes),
            'low_count_pct': float(low_count_pct),
            'log2_range': float(log2_range),
            'skewness': skewness,
            'kurtosis': kurtosis,
            'median_gene_expression': float(gene_means.median()),
            'mean_gene_expression': float(gene_means.mean())
        }
    
    # =========================================================================
    # Biological Variation Estimation
    # =========================================================================
    
    def estimate_biological_variation(self) -> Dict:
        """
        Estimate biological coefficient of variation (BCV).
        
        Returns
        -------
        dict
            BCV and overdispersion metrics
        
        Notes
        -----
        BCV represents expected variation between biological replicates.
        Typical values:
        - Cell lines: BCV ~ 0.1-0.2
        - Genetically identical organisms: BCV ~ 0.2-0.4
        - Human samples: BCV ~ 0.4-0.6
        - Very heterogeneous: BCV > 0.6
        """
        # Calculate per-gene mean and variance
        gene_means = self.counts.mean(axis=1)
        gene_vars = self.counts.var(axis=1)
        
        # Remove genes with zero or very low expression
        mask = gene_means > 5
        gene_means_filtered = gene_means[mask]
        gene_vars_filtered = gene_vars[mask]
        
        if len(gene_means_filtered) == 0:
            logger.warning("No genes with sufficient expression for BCV calculation")
            return {
                'bcv': 0.4,  # Default moderate value
                'overdispersion': 1.0,
                'mean_variance_trend': 'insufficient_data'
            }
        
        # Calculate BCV (biological coefficient of variation)
        # BCV^2 = (variance - mean) / mean^2 for overdispersed count data
        bcv_squared = (gene_vars_filtered - gene_means_filtered) / (gene_means_filtered ** 2)
        bcv_squared = bcv_squared[bcv_squared > 0]  # Keep only positive values
        
        if len(bcv_squared) > 0:
            bcv = float(np.sqrt(np.median(bcv_squared)))
        else:
            bcv = 0.4  # Default
        
        # Overdispersion index (variance / mean)
        overdispersion = float(np.median(gene_vars_filtered / gene_means_filtered))
        
        # Mean-variance trend
        if len(gene_means_filtered) > 10:
            # Fit trend: log(var) ~ log(mean)
            log_mean = np.log10(gene_means_filtered + 1)
            log_var = np.log10(gene_vars_filtered + 1)
            slope, intercept = np.polyfit(log_mean, log_var, 1)
            trend = 'linear' if 0.8 < slope < 1.2 else 'nonlinear'
        else:
            trend = 'insufficient_data'
        
        return {
            'bcv': bcv,
            'overdispersion': overdispersion,
            'mean_variance_trend': trend,
            'n_genes_analyzed': int(len(gene_means_filtered))
        }
    
    # =========================================================================
    # Experimental Design Analysis
    # =========================================================================
    
    def analyze_experimental_design(self) -> Dict:
        """
        Analyze experimental design characteristics.
        
        Returns
        -------
        dict
            Design characteristics including sample size and replication
        
        Notes
        -----
        Sample size impacts statistical power and method choice:
        - n < 3: Very limited, requires robust methods
        - 3 ≤ n < 6: Moderate, standard methods work
        - 6 ≤ n < 10: Good power
        - n ≥ 10: High power, can use faster methods
        """
        design = {
            'n_samples': self.n_samples,
            'n_genes': self.n_genes,
        }
        
        if self.metadata is not None and 'condition' in self.metadata.columns:
            # Analyze conditions
            conditions = self.metadata['condition'].value_counts()
            design['n_conditions'] = len(conditions)
            design['samples_per_condition'] = conditions.to_dict()
            design['min_replicates'] = int(conditions.min())
            design['max_replicates'] = int(conditions.max())
            design['balanced'] = len(conditions.unique()) == 1
        else:
            # Assume two-group comparison
            design['n_conditions'] = 2
            design['min_replicates'] = self.n_samples // 2
            design['max_replicates'] = self.n_samples // 2
            design['balanced'] = self.n_samples % 2 == 0
        
        # Power assessment
        if design['min_replicates'] < 3:
            design['power'] = 'low'
        elif design['min_replicates'] < 6:
            design['power'] = 'moderate'
        elif design['min_replicates'] < 10:
            design['power'] = 'good'
        else:
            design['power'] = 'high'
        
        return design
    
    # =========================================================================
    # Sequencing Quality Assessment
    # =========================================================================
    
    def assess_sequencing_quality(self) -> Dict:
        """
        Assess sequencing depth and quality indicators.
        
        Returns
        -------
        dict
            Sequencing quality metrics
        
        Notes
        -----
        Sequencing depth categories:
        - Low: < 10M reads per sample
        - Medium: 10-25M reads
        - High: 25-50M reads
        - Very high: > 50M reads
        """
        library_sizes = self.counts.sum(axis=0)
        mean_depth = library_sizes.mean()
        
        # Categorize depth
        if mean_depth < 10e6:
            depth_category = 'low'
        elif mean_depth < 25e6:
            depth_category = 'medium'
        elif mean_depth < 50e6:
            depth_category = 'high'
        else:
            depth_category = 'very_high'
        
        # Genes detected
        genes_detected = (self.counts.sum(axis=1) > 0).sum()
        detection_rate = genes_detected / self.n_genes
        
        # Saturation estimate (genes detected vs depth)
        # More genes should be detected with more depth
        if detection_rate > 0.7:
            saturation = 'good'
        elif detection_rate > 0.5:
            saturation = 'moderate'
        else:
            saturation = 'poor'
        
        return {
            'mean_depth': float(mean_depth),
            'depth_category': depth_category,
            'genes_detected': int(genes_detected),
            'detection_rate': float(detection_rate),
            'saturation': saturation
        }
    
    # =========================================================================
    # Summary Generation
    # =========================================================================
    
    def _generate_summary(self, profile: Dict) -> Dict:
        """
        Generate human-readable summary of profile.
        
        Parameters
        ----------
        profile : dict
            Complete profile data
        
        Returns
        -------
        dict
            Summary with difficulty assessment and key characteristics
        """
        # Assess data difficulty
        difficulty_factors = []
        
        # High library variation
        if profile['library_stats']['cv'] > 0.5:
            difficulty_factors.append('high_library_variation')
        
        # High zero-inflation
        if profile['count_distribution']['zero_pct'] > 60:
            difficulty_factors.append('high_zero_inflation')
        
        # Low replication
        if profile['design']['min_replicates'] < 3:
            difficulty_factors.append('low_replication')
        
        # High biological variation
        if profile['biological_variation']['bcv'] > 0.6:
            difficulty_factors.append('high_biological_variation')
        
        # Low depth
        if profile['sequencing']['depth_category'] == 'low':
            difficulty_factors.append('low_sequencing_depth')
        
        # Overall difficulty
        if len(difficulty_factors) == 0:
            difficulty = 'easy'
        elif len(difficulty_factors) <= 2:
            difficulty = 'moderate'
        else:
            difficulty = 'challenging'
        
        return {
            'difficulty': difficulty,
            'difficulty_factors': difficulty_factors,
            'n_difficulty_factors': len(difficulty_factors),
            'recommended_approach': self._suggest_approach(difficulty, profile)
        }
    
    def _suggest_approach(self, difficulty: str, profile: Dict) -> str:
        """Suggest general analysis approach based on difficulty."""
        if difficulty == 'easy':
            return 'Standard pipelines appropriate. Fast methods likely sufficient.'
        elif difficulty == 'moderate':
            return 'Standard pipelines recommended. Consider robust normalization.'
        else:
            return 'Use robust methods with careful normalization and filtering.'


# =============================================================================
# Convenience Functions
# =============================================================================

def profile_from_file(counts_file: str, metadata_file: Optional[str] = None) -> Dict:
    """
    Profile RNA-seq data directly from files.
    
    Parameters
    ----------
    counts_file : str
        Path to count matrix file (CSV or TSV)
    metadata_file : str, optional
        Path to metadata file
    
    Returns
    -------
    dict
        Complete profile data
    
    Examples
    --------
    >>> profile = profile_from_file('counts.csv', 'metadata.csv')
    >>> print(profile['summary']['difficulty'])
    moderate
    """
    # Load counts
    counts = pd.read_csv(counts_file, index_col=0)
    
    # Load metadata if provided
    metadata = None
    if metadata_file:
        metadata = pd.read_csv(metadata_file)
    
    # Create profiler and run
    profiler = RNAseqDataProfiler(counts, metadata)
    return profiler.run_full_profile()


if __name__ == '__main__':
    # Example usage
    print("RAPTOR Data Profiler")
    print("====================")
    print("\nExample usage:")
    print("  from raptor.profiler import RNAseqDataProfiler")
    print("  profiler = RNAseqDataProfiler(counts_df)")
    print("  profile = profiler.run_full_profile()")
