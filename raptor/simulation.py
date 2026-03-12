"""
RNA-seq Data Simulation Module

Generates synthetic RNA-seq count data with KNOWN ground truth differentially
expressed genes. This enables proper evaluation of DE pipelines.

Based on negative binomial model used by DESeq2/edgeR.

Mathematical Model:
    Y_ij ~ NegBinom(μ_ij, φ_i)
    
    Where:
    - Y_ij: count for gene i in sample j
    - μ_ij: mean expression (baseline × size_factor × fold_change)
    - φ_i: gene-specific dispersion
    
    Dispersion-mean relationship:
    φ_i = α_0 + α_1/μ_i  (typical for RNA-seq)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union
import logging
from scipy import stats


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
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Main functions
    'simulate_rnaseq',
    'simulate_diverse_datasets',
    'add_batch_effects',
    
    # Data class
    'SimulationResult',
    
    # Helper functions
    'generate_base_expression',
    'add_differential_expression',
    'generate_counts_from_expression'
]


@dataclass
class SimulationConfig:
    """Configuration for RNA-seq simulation.
    
    Parameters
    ----------
    n_genes : int
        Total number of genes to simulate
    n_samples_per_group : int
        Samples per condition (total = 2 × this)
    de_fraction : float
        Fraction of genes that are DE (0-1)
    dispersion_mean : float
        Mean dispersion (φ). BCV = √φ
        Typical: 0.1 (low), 0.2 (moderate), 0.4+ (high)
    dispersion_sd : float
        Standard deviation of dispersion across genes
    baseline_mean : float
        Mean of log-normal baseline expression
    baseline_sd : float
        SD of log-normal baseline expression
    fold_change_mean : float
        Mean absolute log2 fold change for DE genes
    fold_change_sd : float
        SD of log2 fold changes
    library_size_mean : float
        Mean library size (total counts per sample)
    library_size_cv : float
        Coefficient of variation for library sizes
    dropout_rate : float
        Fraction of counts to set to zero (technical dropout)
    outlier_fraction : float
        Fraction of samples with outlier expression
    seed : int
        Random seed for reproducibility
    """
    n_genes: int = 10000
    n_samples_per_group: int = 4
    de_fraction: float = 0.10
    dispersion_mean: float = 0.2
    dispersion_sd: float = 0.1
    baseline_mean: float = 6.0
    baseline_sd: float = 2.0
    fold_change_mean: float = 1.5
    fold_change_sd: float = 0.5
    library_size_mean: float = 1e7
    library_size_cv: float = 0.2
    dropout_rate: float = 0.0
    outlier_fraction: float = 0.0
    seed: Optional[int] = None


@dataclass
class SimulationResult:
    """Result of RNA-seq simulation.
    
    Attributes
    ----------
    counts : pd.DataFrame
        Simulated count matrix (genes × samples)
    metadata : pd.DataFrame
        Sample metadata with condition labels
    ground_truth : pd.DataFrame
        True DE status and fold changes for each gene
    config : SimulationConfig
        Configuration used for simulation
    """
    counts: pd.DataFrame
    metadata: pd.DataFrame
    ground_truth: pd.DataFrame
    config: SimulationConfig
    
    def __post_init__(self):
        """Validate SimulationResult after initialization."""
        if self.counts is None or len(self.counts) == 0:
            raise ValueError("SimulationResult requires non-empty counts DataFrame")


    def de_genes(self) -> List[str]:
        """List of true DE gene names."""
        return self.ground_truth[self.ground_truth['is_de']].index.tolist()
    
    @property
    def non_de_genes(self) -> List[str]:
        """List of true non-DE gene names."""
        return self.ground_truth[~self.ground_truth['is_de']].index.tolist()
    
    @property
    def n_de(self) -> int:
        """Number of DE genes."""
        return self.ground_truth['is_de'].sum()
    
    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=== Simulation Summary ===",
            f"Genes: {len(self.counts)}",
            f"Samples: {len(self.counts.columns)} ({self.config.n_samples_per_group} per group)",
            f"DE genes: {self.n_de} ({100*self.n_de/len(self.counts):.1f}%)",
            f"Dispersion (BCV): {np.sqrt(self.config.dispersion_mean):.3f}",
            f"Mean library size: {self.counts.sum().mean():,.0f}",
        ]
        return "\n".join(lines)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return {
            'n_genes': len(self.counts),
            'n_samples': len(self.counts.columns),
            'n_de_genes': self.n_de,
            'de_fraction': self.n_de / len(self.counts),
            'dispersion_mean': self.config.dispersion_mean,
            'bcv': np.sqrt(self.config.dispersion_mean),
            'library_size_mean': float(self.counts.sum().mean()),
            'library_size_cv': float(self.counts.sum().std() / self.counts.sum().mean()),
        }


class RNAseqSimulator:
    """
    Simulate RNA-seq count data with known ground truth.
    
    Uses negative binomial model consistent with DESeq2/edgeR assumptions.
    Generates data with controlled:
    - Sample size
    - Dispersion (biological variation)
    - DE gene fraction and effect sizes
    - Library size variation
    - Optional: dropout, outliers
    
    Parameters
    ----------
    config : SimulationConfig
        Simulation parameters
    
    Examples
    --------
    >>> config = SimulationConfig(
    ...     n_genes=10000,
    ...     n_samples_per_group=4,
    ...     de_fraction=0.10,
    ...     dispersion_mean=0.3
    ... )
    >>> simulator = RNAseqSimulator(config)
    >>> result = simulator.simulate()
    >>> print(f"Generated {result.n_de} DE genes")
    """
    
    @handle_errors(exit_on_error=False)
    def __init__(self, config: Optional[SimulationConfig] = None):
        """Initialize simulator with configuration."""
        self.config = config or SimulationConfig()
        
        if self.config.seed is not None:
            np.random.seed(self.config.seed)
        
        logger.info(f"Initialized simulator: {self.config.n_genes} genes, "
                   f"{self.config.n_samples_per_group} samples/group")
    
    @handle_errors(exit_on_error=False)
    def simulate(self) -> SimulationResult:
        """
        Run simulation and return results.
        
        Returns
        -------
        SimulationResult
            Contains counts, metadata, ground truth, and config
        """
        logger.info("Starting RNA-seq simulation...")
        
        n_samples = 2 * self.config.n_samples_per_group
        n_de = int(self.config.n_genes * self.config.de_fraction)
        
        # Step 1: Generate gene-level parameters
        baseline_expr, dispersions = self._generate_gene_params()
        
        # Step 2: Assign DE status and fold changes
        is_de, fold_changes = self._assign_de_genes(n_de)
        
        # Step 3: Generate library size factors
        size_factors = self._generate_size_factors(n_samples)
        
        # Step 4: Generate counts using negative binomial
        counts = self._generate_counts(
            baseline_expr, dispersions, fold_changes, size_factors
        )
        
        # Step 5: Apply dropout if specified
        if self.config.dropout_rate > 0:
            counts = self._apply_dropout(counts)
        
        # Step 6: Add outliers if specified
        if self.config.outlier_fraction > 0:
            counts = self._add_outliers(counts)
        
        # Create DataFrames
        gene_names = [f"Gene_{i:05d}" for i in range(self.config.n_genes)]
        sample_names = [f"Sample_{i:02d}" for i in range(n_samples)]
        
        counts_df = pd.DataFrame(
            counts,
            index=gene_names,
            columns=sample_names
        )
        
        # Create metadata
        conditions = ['Control'] * self.config.n_samples_per_group + \
                    ['Treatment'] * self.config.n_samples_per_group
        metadata_df = pd.DataFrame({
            'sample_id': sample_names,
            'condition': conditions,
            'group': [0] * self.config.n_samples_per_group + \
                    [1] * self.config.n_samples_per_group
        })
        
        # Create ground truth
        ground_truth_df = pd.DataFrame({
            'is_de': is_de,
            'log2_fold_change': fold_changes,
            'baseline_expression': baseline_expr,
            'dispersion': dispersions
        }, index=gene_names)
        
        result = SimulationResult(
            counts=counts_df,
            metadata=metadata_df,
            ground_truth=ground_truth_df,
            config=self.config
        )
        
        logger.info(f"Simulation complete: {n_de} DE genes, "
                   f"mean library size {counts_df.sum().mean():,.0f}")
        
        return result
    
    def _generate_gene_params(self) -> Tuple[np.ndarray, np.ndarray]:
        """Generate baseline expression and dispersion for each gene."""
        
        # Baseline expression: log-normal distribution
        # This gives realistic range from low to highly expressed
        log_baseline = np.random.normal(
            self.config.baseline_mean,
            self.config.baseline_sd,
            self.config.n_genes
        )
        baseline_expr = np.exp(log_baseline)
        
        # Dispersion: gamma distribution (always positive, right-skewed)
        # Shape parameter controls the spread
        shape = (self.config.dispersion_mean / self.config.dispersion_sd) ** 2
        scale = self.config.dispersion_sd ** 2 / self.config.dispersion_mean
        
        dispersions = np.random.gamma(shape, scale, self.config.n_genes)
        
        # Add mean-dispersion trend (typical for RNA-seq)
        # Higher dispersion for lowly expressed genes
        trend_factor = 0.1 / (baseline_expr / baseline_expr.mean() + 0.1)
        dispersions = dispersions * (1 + 0.5 * trend_factor)
        
        # Ensure minimum dispersion
        dispersions = np.maximum(dispersions, 0.01)
        
        return baseline_expr, dispersions
    
    def _assign_de_genes(self, n_de: int) -> Tuple[np.ndarray, np.ndarray]:
        """Assign DE status and fold changes."""
        
        # Random selection of DE genes
        is_de = np.zeros(self.config.n_genes, dtype=bool)
        de_indices = np.random.choice(
            self.config.n_genes, n_de, replace=False
        )
        is_de[de_indices] = True
        
        # Generate fold changes (mixture of up and down)
        fold_changes = np.zeros(self.config.n_genes)
        
        # Log2 fold changes from normal distribution
        lfc = np.abs(np.random.normal(
            self.config.fold_change_mean,
            self.config.fold_change_sd,
            n_de
        ))
        
        # Randomly assign direction (up or down)
        directions = np.random.choice([-1, 1], n_de)
        fold_changes[de_indices] = lfc * directions
        
        return is_de, fold_changes
    
    def _generate_size_factors(self, n_samples: int) -> np.ndarray:
        """Generate library size factors for each sample."""
        
        # Log-normal size factors
        log_factors = np.random.normal(
            0,
            self.config.library_size_cv,
            n_samples
        )
        size_factors = np.exp(log_factors)
        
        # Scale to achieve target mean library size
        size_factors = size_factors * self.config.library_size_mean / \
                      (size_factors.mean() * self.config.n_genes * 100)
        
        return size_factors
    
    def _generate_counts(
        self,
        baseline_expr: np.ndarray,
        dispersions: np.ndarray,
        fold_changes: np.ndarray,
        size_factors: np.ndarray
    ) -> np.ndarray:
        """Generate counts using negative binomial distribution."""
        
        n_samples = len(size_factors)
        n_per_group = n_samples // 2
        counts = np.zeros((self.config.n_genes, n_samples))
        
        for j in range(n_samples):
            # Determine if treatment group
            is_treatment = j >= n_per_group
            
            for i in range(self.config.n_genes):
                # Calculate mean for this gene-sample combination
                mu = baseline_expr[i] * size_factors[j]
                
                # Apply fold change for treatment samples
                if is_treatment and fold_changes[i] != 0:
                    mu = mu * (2 ** fold_changes[i])
                
                # Negative binomial parameters
                # NB parameterization: mean=mu, var=mu + mu^2 * dispersion
                phi = dispersions[i]
                
                # Convert to scipy's parameterization
                # scipy uses: n (number of successes), p (probability)
                # n = 1/phi, p = 1/(1 + mu*phi)
                n_param = 1 / phi
                p_param = n_param / (n_param + mu)
                
                # Generate count
                counts[i, j] = np.random.negative_binomial(n_param, p_param)
        
        return counts.astype(int)
    
    def _apply_dropout(self, counts: np.ndarray) -> np.ndarray:
        """Apply technical dropout (set random counts to zero)."""
        
        # Dropout probability inversely related to expression
        # More dropout for lowly expressed genes
        mean_expr = counts.mean(axis=1, keepdims=True)
        dropout_prob = self.config.dropout_rate * np.exp(-mean_expr / 100)
        
        # Apply dropout
        mask = np.random.random(counts.shape) < dropout_prob
        counts[mask] = 0
        
        return counts
    
    def _add_outliers(self, counts: np.ndarray) -> np.ndarray:
        """Add outlier samples with extreme expression."""
        
        n_outlier_samples = max(1, int(counts.shape[1] * self.config.outlier_fraction))
        outlier_indices = np.random.choice(
            counts.shape[1], n_outlier_samples, replace=False
        )
        
        for idx in outlier_indices:
            # Multiply random subset of genes by large factor
            n_outlier_genes = int(0.05 * self.config.n_genes)
            gene_indices = np.random.choice(
                self.config.n_genes, n_outlier_genes, replace=False
            )
            
            # Random multiplier (2-10x)
            multipliers = np.random.uniform(2, 10, n_outlier_genes)
            counts[gene_indices, idx] = counts[gene_indices, idx] * multipliers
        
        return counts.astype(int)


# =============================================================================
# Convenience Functions
# =============================================================================

def simulate_rnaseq(
    n_genes: int = 10000,
    n_samples_per_group: int = 4,
    de_fraction: float = 0.10,
    dispersion: float = 0.2,
    seed: Optional[int] = None
) -> SimulationResult:
    """
    Quick simulation with common parameters.
    
    Parameters
    ----------
    n_genes : int
        Number of genes (default: 10000)
    n_samples_per_group : int
        Samples per condition (default: 4)
    de_fraction : float
        Fraction of DE genes (default: 0.10)
    dispersion : float
        Mean dispersion, BCV = √dispersion (default: 0.2)
    seed : int, optional
        Random seed
    
    Returns
    -------
    SimulationResult
        Simulation results with counts, metadata, ground truth
    
    Examples
    --------
    >>> result = simulate_rnaseq(n_samples_per_group=6, dispersion=0.3)
    >>> print(f"DE genes: {result.n_de}")
    """
    config = SimulationConfig(
        n_genes=n_genes,
        n_samples_per_group=n_samples_per_group,
        de_fraction=de_fraction,
        dispersion_mean=dispersion,
        seed=seed
    )
    
    simulator = RNAseqSimulator(config)
    return simulator.simulate()


def simulate_diverse_datasets(
    n_datasets: int = 100,
    seed: int = 42
) -> List[SimulationResult]:
    """
    Generate diverse datasets for ML training.
    
    Creates datasets spanning the full range of RNA-seq characteristics:
    - Sample sizes: 3-30 per group
    - Dispersions: 0.1-0.6 (BCV: 0.3-0.77)
    - DE fractions: 5-20%
    - Library size CV: 0.1-0.4
    
    Parameters
    ----------
    n_datasets : int
        Number of datasets to generate
    seed : int
        Random seed for reproducibility
    
    Returns
    -------
    List[SimulationResult]
        List of simulation results
    """
    np.random.seed(seed)
    
    datasets = []
    
    # Parameter ranges
    sample_sizes = [3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30]
    dispersions = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]
    de_fractions = [0.05, 0.10, 0.15, 0.20]
    lib_cvs = [0.1, 0.2, 0.3, 0.4]
    
    for i in range(n_datasets):
        config = SimulationConfig(
            n_genes=10000,
            n_samples_per_group=np.random.choice(sample_sizes),
            de_fraction=np.random.choice(de_fractions),
            dispersion_mean=np.random.choice(dispersions),
            library_size_cv=np.random.choice(lib_cvs),
            seed=seed + i
        )
        
        simulator = RNAseqSimulator(config)
        result = simulator.simulate()
        datasets.append(result)
        
        if (i + 1) % 20 == 0:
            logger.info(f"Generated {i + 1}/{n_datasets} datasets")
    
    return datasets


# =============================================================================
# Preset Scenarios
# =============================================================================

def simulate_small_sample_high_dispersion(seed: int = None) -> SimulationResult:
    """Scenario: Small samples, high biological variation → edgeR preferred."""
    config = SimulationConfig(
        n_genes=10000,
        n_samples_per_group=3,
        de_fraction=0.10,
        dispersion_mean=0.5,  # BCV ≈ 0.71
        seed=seed
    )
    return RNAseqSimulator(config).simulate()


def simulate_large_sample_low_dispersion(seed: int = None) -> SimulationResult:
    """Scenario: Large samples, low variation → limma-voom/Wilcoxon preferred."""
    config = SimulationConfig(
        n_genes=10000,
        n_samples_per_group=25,
        de_fraction=0.10,
        dispersion_mean=0.1,  # BCV ≈ 0.32
        seed=seed
    )
    return RNAseqSimulator(config).simulate()


def simulate_medium_balanced(seed: int = None) -> SimulationResult:
    """Scenario: Medium samples, moderate variation → DESeq2/edgeR both good."""
    config = SimulationConfig(
        n_genes=10000,
        n_samples_per_group=6,
        de_fraction=0.10,
        dispersion_mean=0.2,  # BCV ≈ 0.45
        seed=seed
    )
    return RNAseqSimulator(config).simulate()


def simulate_with_outliers(seed: int = None) -> SimulationResult:
    """Scenario: Dataset with outlier samples → edgeR_robust preferred."""
    config = SimulationConfig(
        n_genes=10000,
        n_samples_per_group=6,
        de_fraction=0.10,
        dispersion_mean=0.2,
        outlier_fraction=0.15,  # 15% of samples are outliers
        seed=seed
    )
    return RNAseqSimulator(config).simulate()


if __name__ == '__main__':
    print("🦖 RAPTOR RNA-seq Simulation Module")
    print("=" * 50)
    
    # Demo simulation
    print("\nGenerating example dataset...")
    result = simulate_rnaseq(
        n_samples_per_group=4,
        dispersion=0.3,
        seed=42
    )
    
    print(result.summary())
    
    print("\n--- Count Matrix Preview ---")
    print(result.counts.iloc[:5, :5])
    
    print("\n--- Ground Truth Preview ---")
    de_genes = result.ground_truth[result.ground_truth['is_de']].head()
    print(de_genes)
