"""
RNA-seq Data Simulator

Generates simulated RNA-seq data with known ground truth for benchmarking
and validation purposes.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Dict
import logging

logger = logging.getLogger(__name__)


class DataSimulator:
    """
    Generate simulated RNA-seq data.
    
    Creates synthetic count data with known differentially expressed genes
    for benchmarking pipeline performance.
    
    Parameters
    ----------
    n_genes : int
        Total number of genes
    n_samples : int
        Total number of samples
    n_de : int
        Number of differentially expressed genes
    fold_changes : list of float
        Fold changes for DE genes (e.g., [0.5, 2.0])
    seed : int, optional
        Random seed for reproducibility
    
    Examples
    --------
    >>> simulator = DataSimulator(
    ...     n_genes=2000,
    ...     n_samples=6,
    ...     n_de=400,
    ...     fold_changes=[0.5, 2.0]
    ... )
    >>> simulator.generate_data('simulated_data/')
    """
    
    def __init__(
        self,
        n_genes: int = 2000,
        n_samples: int = 6,
        n_de: int = 400,
        fold_changes: List[float] = [0.5, 2.0],
        seed: Optional[int] = None
    ):
        """Initialize data simulator."""
        self.n_genes = n_genes
        self.n_samples = n_samples
        self.n_de = n_de
        self.fold_changes = fold_changes
        self.seed = seed
        
        if seed is not None:
            np.random.seed(seed)
        
        logger.info(f"Initialized simulator: {n_genes} genes, {n_samples} samples, "
                   f"{n_de} DE genes")
    
    def generate_data(self, output_dir: str) -> Dict:
        """
        Generate simulated data and save to files.
        
        Parameters
        ----------
        output_dir : str
            Output directory for simulated data
        
        Returns
        -------
        dict
            Information about generated data including ground truth
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        logger.info("Generating simulated RNA-seq data...")
        
        # Generate counts
        counts, ground_truth = self._simulate_counts()
        
        # Create metadata
        metadata = self._create_metadata()
        
        # Save files
        counts_file = output_path / 'counts.csv'
        counts.to_csv(counts_file)
        logger.info(f"Saved counts: {counts_file}")
        
        metadata_file = output_path / 'metadata.csv'
        metadata.to_csv(metadata_file, index=False)
        logger.info(f"Saved metadata: {metadata_file}")
        
        ground_truth_file = output_path / 'ground_truth.csv'
        ground_truth.to_csv(ground_truth_file)
        logger.info(f"Saved ground truth: {ground_truth_file}")
        
        # Create summary
        summary = {
            'n_genes': self.n_genes,
            'n_samples': self.n_samples,
            'n_de': self.n_de,
            'fold_changes': self.fold_changes,
            'files': {
                'counts': str(counts_file),
                'metadata': str(metadata_file),
                'ground_truth': str(ground_truth_file)
            }
        }
        
        return summary
    
    def _simulate_counts(self) -> tuple:
        """
        Simulate count matrix with DE genes.
        
        Returns
        -------
        tuple
            (counts_dataframe, ground_truth_dataframe)
        """
        # Number of samples per condition
        n_per_condition = self.n_samples // 2
        
        # Base expression levels (log2 scale)
        # Simulate realistic gene expression distribution
        base_expression = np.random.gamma(shape=2, scale=100, size=self.n_genes)
        
        # Select DE genes
        de_indices = np.random.choice(self.n_genes, self.n_de, replace=False)
        de_mask = np.zeros(self.n_genes, dtype=bool)
        de_mask[de_indices] = True
        
        # Assign fold changes to DE genes
        fc_assignments = np.random.choice(self.fold_changes, self.n_de)
        
        # Initialize counts matrix
        counts = np.zeros((self.n_genes, self.n_samples))
        
        # Generate counts for condition 1 (control)
        for i in range(n_per_condition):
            # Add biological variation (negative binomial)
            size_param = 10  # Dispersion parameter
            counts[:, i] = np.random.negative_binomial(
                size_param,
                size_param / (size_param + base_expression)
            )
        
        # Generate counts for condition 2 (treatment)
        treatment_expression = base_expression.copy()
        treatment_expression[de_mask] *= fc_assignments
        
        for i in range(n_per_condition, self.n_samples):
            size_param = 10
            counts[:, i] = np.random.negative_binomial(
                size_param,
                size_param / (size_param + treatment_expression)
            )
        
        # Create DataFrame
        gene_names = [f'Gene_{i+1:05d}' for i in range(self.n_genes)]
        sample_names = [f'Sample_{i+1}' for i in range(self.n_samples)]
        
        counts_df = pd.DataFrame(
            counts.astype(int),
            index=gene_names,
            columns=sample_names
        )
        
        # Create ground truth
        ground_truth_df = pd.DataFrame({
            'gene': gene_names,
            'is_de': de_mask,
            'fold_change': np.ones(self.n_genes)
        })
        ground_truth_df.loc[de_mask, 'fold_change'] = fc_assignments
        ground_truth_df.set_index('gene', inplace=True)
        
        logger.info(f"Generated {self.n_genes} genes with {self.n_de} DE genes")
        
        return counts_df, ground_truth_df
    
    def _create_metadata(self) -> pd.DataFrame:
        """
        Create sample metadata.
        
        Returns
        -------
        pd.DataFrame
            Sample metadata
        """
        n_per_condition = self.n_samples // 2
        
        conditions = ['Control'] * n_per_condition + ['Treatment'] * n_per_condition
        sample_names = [f'Sample_{i+1}' for i in range(self.n_samples)]
        
        metadata = pd.DataFrame({
            'sample': sample_names,
            'condition': conditions,
            'replicate': list(range(1, n_per_condition + 1)) * 2
        })
        
        return metadata


# =============================================================================
# Convenience Functions
# =============================================================================

def quick_simulate(output_dir: str, size: str = 'small') -> Dict:
    """
    Quick simulation with preset parameters.
    
    Parameters
    ----------
    output_dir : str
        Output directory
    size : str
        Dataset size: 'small', 'medium', or 'large'
    
    Returns
    -------
    dict
        Simulation summary
    
    Examples
    --------
    >>> summary = quick_simulate('test_data/', size='small')
    >>> print(summary['n_genes'])
    1000
    """
    if size == 'small':
        params = {
            'n_genes': 1000,
            'n_samples': 6,
            'n_de': 200,
            'fold_changes': [0.5, 2.0]
        }
    elif size == 'medium':
        params = {
            'n_genes': 5000,
            'n_samples': 12,
            'n_de': 1000,
            'fold_changes': [0.5, 2.0, 4.0]
        }
    elif size == 'large':
        params = {
            'n_genes': 20000,
            'n_samples': 24,
            'n_de': 4000,
            'fold_changes': [0.25, 0.5, 2.0, 4.0]
        }
    else:
        raise ValueError(f"Size must be 'small', 'medium', or 'large', got '{size}'")
    
    simulator = DataSimulator(**params, seed=42)
    return simulator.generate_data(output_dir)


if __name__ == '__main__':
    print("RAPTOR Data Simulator")
    print("====================")
    print("\nGenerate simulated RNA-seq data for testing.")
    print("Use via CLI: raptor simulate --output sim_data/")
