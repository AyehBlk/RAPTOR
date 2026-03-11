"""
Synthetic Benchmark Data Generator for DE Pipeline Selection

Generates realistic synthetic benchmark data for training the ML-based
pipeline recommender. Creates data profiles and simulated pipeline performance
based on published benchmarking studies.

Key Literature:
- Love et al. (2014) - DESeq2: Genome Biology 15:550
- Robinson et al. (2010) - edgeR: Bioinformatics 26:139-140
- Law et al. (2014) - limma-voom: Genome Biology 15:R29
- Li et al. (2022) - Large sample comparison: Genome Biology 23:79
- Schurch et al. (2016) - Sample size: RNA 22:839-851

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
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
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Main functions
    'generate_benchmark_datasets',
    'create_benchmark_scenarios',
    'evaluate_pipeline_performance',
    
    # Data classes
    'BenchmarkDataset',
    'BenchmarkScenario',
    'PipelinePerformance',
    
    # Helper functions
    'run_pipeline_on_benchmark',
    'calculate_performance_metrics'
]


# =============================================================================
# Pipeline Performance Models (Literature-Based)
# =============================================================================

@dataclass
class PipelinePerformanceModel:
    """
    Literature-based performance model for each DE method.
    
    Encodes known strengths and weaknesses from benchmarking studies.
    """
    name: str
    base_f1: float
    min_replicates: int
    optimal_sample_range: Tuple[int, int]
    handles_high_bcv: bool
    handles_low_counts: bool
    handles_outliers: bool
    handles_large_samples: bool
    speed_factor: float  # 1.0 = baseline, higher = faster
    
    def calculate_f1(
        self,
        n_samples: int,
        min_group_size: int,
        bcv: float,
        low_count_prop: float,
        outlier_severity: float,
        sparsity: float
    ) -> float:
        """
        Calculate expected F1 score based on data characteristics.
        
        This implements the literature-based performance expectations.
        """
        f1 = self.base_f1
        
        # Sample size effects
        if min_group_size < self.min_replicates:
            f1 -= 0.20  # Cannot run properly
        elif min_group_size < 3:
            f1 -= 0.10  # Underpowered
        elif min_group_size >= self.optimal_sample_range[0]:
            if min_group_size <= self.optimal_sample_range[1]:
                f1 += 0.05  # Optimal range
            elif self.handles_large_samples:
                f1 += 0.03  # Handles large samples well
            else:
                f1 -= 0.08  # FDR inflation for large samples
        
        # BCV effects (biological coefficient of variation)
        if bcv > 0.4:
            if self.handles_high_bcv:
                f1 += 0.05
            else:
                f1 -= 0.08
        elif bcv < 0.2:
            if self.name == 'limma-voom':
                f1 += 0.05  # limma excels with low BCV
        
        # Low count effects
        if low_count_prop > 0.3:
            if self.handles_low_counts:
                f1 += 0.05
            else:
                f1 -= 0.05
        
        # Outlier effects
        if outlier_severity > 0.3:
            if self.handles_outliers:
                f1 += 0.08
            else:
                f1 -= 0.10
        
        # Sparsity effects
        if sparsity > 0.6:
            f1 -= 0.05  # All methods struggle with extreme sparsity
        
        # Add realistic noise
        f1 += np.random.normal(0, 0.02)
        
        # Clamp to realistic range
        return max(0.40, min(0.95, f1))


# Define performance models for each DE method
PIPELINE_MODELS = {
    'DESeq2': PipelinePerformanceModel(
        name='DESeq2',
        base_f1=0.82,
        min_replicates=2,
        optimal_sample_range=(3, 30),
        handles_high_bcv=True,
        handles_low_counts=False,
        handles_outliers=True,  # Cook's distance
        handles_large_samples=False,  # FDR inflation per Li et al.
        speed_factor=0.5
    ),
    'edgeR': PipelinePerformanceModel(
        name='edgeR',
        base_f1=0.81,
        min_replicates=2,
        optimal_sample_range=(3, 30),
        handles_high_bcv=True,  # Good dispersion estimation
        handles_low_counts=True,  # Best for low counts
        handles_outliers=False,
        handles_large_samples=False,  # FDR inflation per Li et al.
        speed_factor=1.0
    ),
    'limma-voom': PipelinePerformanceModel(
        name='limma-voom',
        base_f1=0.80,
        min_replicates=3,
        optimal_sample_range=(8, 100),
        handles_high_bcv=False,
        handles_low_counts=False,
        handles_outliers=True,  # robust=TRUE option
        handles_large_samples=True,  # Scales well
        speed_factor=2.0
    ),
    'Wilcoxon': PipelinePerformanceModel(
        name='Wilcoxon',
        base_f1=0.75,
        min_replicates=8,
        optimal_sample_range=(8, 200),
        handles_high_bcv=True,  # Non-parametric
        handles_low_counts=False,
        handles_outliers=True,  # Non-parametric = robust
        handles_large_samples=True,  # Best FDR control per Li et al.
        speed_factor=1.5
    ),
    'edgeR_robust': PipelinePerformanceModel(
        name='edgeR_robust',
        base_f1=0.79,
        min_replicates=2,
        optimal_sample_range=(3, 30),
        handles_high_bcv=True,
        handles_low_counts=True,
        handles_outliers=True,  # Specifically designed for outliers
        handles_large_samples=False,
        speed_factor=0.8
    )
}


# =============================================================================
# Synthetic Data Profile Generator
# =============================================================================

class SyntheticProfileGenerator:
    """
    Generate realistic RNA-seq data profiles for training.
    
    Creates diverse profiles covering the full range of experimental
    conditions encountered in real RNA-seq studies.
    """
    
    @handle_errors(exit_on_error=False)
    def __init__(self, seed: int = 42):
        """Initialize generator with random seed."""
        self.seed = seed
        np.random.seed(seed)
        
    def generate_profile(self) -> Dict:
        """
        Generate a single realistic data profile.
        
        Returns profile dictionary compatible with the 32-feature profiler.
        """
        # Sample design (weighted toward common configurations)
        n_samples = np.random.choice(
            [4, 6, 8, 10, 12, 16, 20, 30, 50],
            p=[0.10, 0.20, 0.20, 0.15, 0.12, 0.10, 0.07, 0.04, 0.02]
        )
        n_groups = np.random.choice([2, 3, 4], p=[0.70, 0.20, 0.10])
        
        # Group sizes (may be unbalanced)
        if np.random.random() < 0.7:
            # Balanced design
            group_sizes = [n_samples // n_groups] * n_groups
        else:
            # Unbalanced design
            group_sizes = np.random.multinomial(
                n_samples, 
                [1/n_groups] * n_groups
            ).tolist()
        
        min_group_size = min(group_sizes)
        max_group_size = max(group_sizes)
        sample_balance = min_group_size / max_group_size if max_group_size > 0 else 1.0
        
        # Gene count
        n_genes = np.random.choice(
            [10000, 15000, 20000, 25000, 30000],
            p=[0.10, 0.25, 0.35, 0.20, 0.10]
        )
        
        # Library size characteristics
        library_size_mean = np.random.lognormal(mean=np.log(20e6), sigma=0.4)
        library_size_cv = np.random.beta(2, 10)  # Usually low CV
        
        # Dispersion/BCV (varies by sample type)
        # Cell lines: 0.1-0.2, Tissues: 0.2-0.4, Human samples: 0.4-0.8
        sample_type = np.random.choice(['cell_line', 'tissue', 'human'], p=[0.2, 0.4, 0.4])
        if sample_type == 'cell_line':
            bcv = np.random.uniform(0.08, 0.25)
        elif sample_type == 'tissue':
            bcv = np.random.uniform(0.20, 0.50)
        else:
            bcv = np.random.uniform(0.35, 0.80)
        
        common_dispersion = bcv ** 2
        
        # Sparsity and zero-inflation
        sparsity = np.random.beta(3, 4)  # 0-1, tends toward middle
        zero_inflation_index = np.random.exponential(0.3)  # Usually low
        
        # Count distribution
        low_count_proportion = np.random.beta(2, 5) * 0.5  # 0-50%
        
        # Outlier presence
        has_outliers = np.random.random() < 0.25
        if has_outliers:
            outlier_severity = np.random.choice(
                ['mild', 'moderate', 'severe'],
                p=[0.5, 0.35, 0.15]
            )
            n_outliers = np.random.randint(1, max(2, n_samples // 4))
        else:
            outlier_severity = 'none'
            n_outliers = 0
        
        # Expression distribution characteristics
        expression_mean = np.log2(library_size_mean / n_genes)
        expression_variance = np.random.uniform(8, 16)
        
        # Detection rates
        detection_rate = 1 - sparsity * 0.5
        reliable_detection_rate = detection_rate * np.random.uniform(0.6, 0.9)
        
        # Mean-variance relationship
        mean_var_slope = np.random.uniform(1.3, 2.0)  # NB typically 1.5-2
        mean_var_r_squared = np.random.uniform(0.85, 0.98)
        
        # Build profile dictionary matching DataProfile structure
        profile = {
            # Sample characteristics
            'n_samples': int(n_samples),
            'n_genes': int(n_genes),
            'n_groups': int(n_groups),
            'min_group_size': int(min_group_size),
            'sample_balance': float(sample_balance),
            
            # Library size metrics
            'library_size_mean': float(library_size_mean),
            'library_size_cv': float(library_size_cv),
            'library_size_range': float(library_size_mean * library_size_cv * 4),
            'library_size_iqr': float(library_size_mean * library_size_cv * 1.5),
            
            # Gene detection
            'detection_rate': float(detection_rate),
            'reliable_detection_rate': float(reliable_detection_rate),
            'genes_per_sample_cv': float(np.random.uniform(0.02, 0.15)),
            
            # Expression distribution
            'expression_mean': float(expression_mean),
            'expression_variance': float(expression_variance),
            'expression_skewness': float(np.random.uniform(0.5, 2.0)),
            'expression_iqr': float(np.sqrt(expression_variance) * 1.35),
            
            # Dispersion estimates (CRITICAL)
            'bcv': float(bcv),
            'common_dispersion': float(common_dispersion),
            'trended_dispersion_slope': float(np.random.uniform(-0.3, 0.1)),
            'overdispersion_ratio': float(1 + common_dispersion * np.random.uniform(0.5, 2.0)),
            'dispersion_iqr': float(common_dispersion * np.random.uniform(0.3, 0.8)),
            
            # Sparsity metrics
            'sparsity': float(sparsity),
            'zero_inflation_index': float(zero_inflation_index),
            'dropout_rate': float(sparsity * np.random.uniform(0.8, 1.2)),
            
            # Count distribution
            'low_count_proportion': float(low_count_proportion),
            'medium_count_proportion': float(np.random.uniform(0.2, 0.4)),
            'high_count_proportion': float(1 - low_count_proportion - 0.3),
            'dynamic_range': float(np.random.uniform(14, 22)),
            
            # Mean-variance relationship
            'mean_var_slope': float(mean_var_slope),
            'mean_var_r_squared': float(mean_var_r_squared),
            'poisson_fit_score': float(1.0 / mean_var_slope),
            
            # Quality
            'quality_score': float(np.random.uniform(60, 95)),
            
            # Metadata (for simulation tracking)
            '_sample_type': sample_type,
            '_outlier_severity': outlier_severity,
            '_n_outliers': int(n_outliers),
            '_group_sizes': group_sizes
        }
        
        return profile


# =============================================================================
# Benchmark Results Generator
# =============================================================================

class SyntheticBenchmarkGenerator:
    """
    Generate synthetic benchmark datasets for ML training.
    
    Creates data profiles and simulated pipeline performance results
    based on literature-informed performance models.
    
    Parameters
    ----------
    n_datasets : int
        Number of synthetic datasets to generate
    seed : int
        Random seed for reproducibility
    
    Examples
    --------
    >>> generator = SyntheticBenchmarkGenerator(n_datasets=200, seed=42)
    >>> summary = generator.generate_benchmarks('training_data/')
    >>> print(f"Generated {summary['n_datasets']} datasets")
    """
    
    def __init__(self, n_datasets: int = 200, seed: int = 42):
        """Initialize generator."""
        self.n_datasets = n_datasets
        self.seed = seed
        self.profile_generator = SyntheticProfileGenerator(seed=seed)
        
        np.random.seed(seed)
        logger.info(f"Initialized SyntheticBenchmarkGenerator: {n_datasets} datasets")
    
    def generate_benchmarks(self, output_dir: str) -> Dict:
        """
        Generate complete benchmark dataset for ML training.
        
        Parameters
        ----------
        output_dir : str
            Output directory for benchmark data
        
        Returns
        -------
        dict
            Summary including pipeline distribution and file paths
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Generating {self.n_datasets} synthetic benchmarks...")
        
        # Track pipeline distribution
        pipeline_distribution = {name: 0 for name in PIPELINE_MODELS.keys()}
        all_profiles = []
        all_results = []
        
        for i in range(self.n_datasets):
            # Generate profile
            profile = self.profile_generator.generate_profile()
            
            # Generate benchmark results
            benchmark_results = self._simulate_pipeline_performance(profile)
            
            # Determine best pipeline
            best_pipeline = max(
                benchmark_results.items(),
                key=lambda x: x[1]['f1_score']
            )[0]
            
            pipeline_distribution[best_pipeline] += 1
            
            # Store for batch saving
            profile['_dataset_id'] = i
            profile['_best_pipeline'] = best_pipeline
            all_profiles.append(profile)
            
            benchmark_results['_dataset_id'] = i
            benchmark_results['_best_pipeline'] = best_pipeline
            all_results.append(benchmark_results)
            
            # Progress logging
            if (i + 1) % 50 == 0:
                logger.info(f"Generated {i + 1}/{self.n_datasets} datasets")
        
        # Save profiles as JSON
        profiles_file = output_path / 'profiles.json'
        with open(profiles_file, 'w') as f:
            json.dump(all_profiles, f, indent=2)
        
        # Save results as JSON
        results_file = output_path / 'benchmark_results.json'
        with open(results_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        # Create feature matrix for ML training
        X, y, feature_names = self._create_training_data(all_profiles, all_results)
        
        # Save as CSV for easy inspection
        features_df = pd.DataFrame(X, columns=feature_names)
        features_df['best_pipeline'] = y
        features_df.to_csv(output_path / 'training_data.csv', index=False)
        
        # Save feature names
        with open(output_path / 'feature_names.json', 'w') as f:
            json.dump(feature_names, f)
        
        # Summary
        summary = {
            'n_datasets': self.n_datasets,
            'output_dir': str(output_path),
            'pipeline_distribution': pipeline_distribution,
            'n_features': len(feature_names),
            'feature_names': feature_names,
            'seed': self.seed,
            'files': {
                'profiles': str(profiles_file),
                'results': str(results_file),
                'training_data': str(output_path / 'training_data.csv'),
                'feature_names': str(output_path / 'feature_names.json')
            }
        }
        
        # Save summary
        with open(output_path / 'generation_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Generation complete: {self.n_datasets} datasets")
        logger.info(f"Pipeline distribution: {pipeline_distribution}")
        
        return summary
    
    def _simulate_pipeline_performance(self, profile: Dict) -> Dict:
        """
        Simulate performance of each pipeline on the given data profile.
        
        Uses literature-based performance models to generate realistic
        F1 scores for each DE method.
        """
        results = {}
        
        # Extract key characteristics
        n_samples = profile['n_samples']
        min_group_size = profile['min_group_size']
        bcv = profile['bcv']
        low_count_prop = profile['low_count_proportion']
        sparsity = profile['sparsity']
        
        # Map outlier severity to numeric
        outlier_map = {'none': 0, 'mild': 0.2, 'moderate': 0.5, 'severe': 0.8}
        outlier_severity = outlier_map.get(profile.get('_outlier_severity', 'none'), 0)
        
        # Calculate F1 for each pipeline
        for name, model in PIPELINE_MODELS.items():
            f1_score = model.calculate_f1(
                n_samples=n_samples,
                min_group_size=min_group_size,
                bcv=bcv,
                low_count_prop=low_count_prop,
                outlier_severity=outlier_severity,
                sparsity=sparsity
            )
            
            # Generate other metrics with correlated noise
            precision = f1_score + np.random.normal(0, 0.03)
            recall = f1_score + np.random.normal(0, 0.03)
            
            # Runtime (scaled by data size and pipeline speed)
            base_runtime = 60 * n_samples / model.speed_factor
            runtime = base_runtime * (1 + np.random.uniform(-0.2, 0.2))
            
            results[name] = {
                'f1_score': float(f1_score),
                'precision': float(max(0.3, min(1.0, precision))),
                'recall': float(max(0.3, min(1.0, recall))),
                'runtime_seconds': float(runtime),
                'can_run': min_group_size >= model.min_replicates
            }
        
        return results
    
    def _create_training_data(
        self,
        profiles: List[Dict],
        results: List[Dict]
    ) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        Create feature matrix and labels for ML training.
        
        Returns
        -------
        X : np.ndarray
            Feature matrix (n_samples, n_features)
        y : list
            Best pipeline labels
        feature_names : list
            Names of features
        """
        # Define features to extract (matching profiler output)
        feature_keys = [
            'n_samples', 'n_genes', 'n_groups', 'min_group_size', 'sample_balance',
            'library_size_mean', 'library_size_cv', 'library_size_range', 'library_size_iqr',
            'detection_rate', 'reliable_detection_rate', 'genes_per_sample_cv',
            'expression_mean', 'expression_variance', 'expression_skewness', 'expression_iqr',
            'bcv', 'common_dispersion', 'trended_dispersion_slope', 
            'overdispersion_ratio', 'dispersion_iqr',
            'sparsity', 'zero_inflation_index', 'dropout_rate',
            'low_count_proportion', 'medium_count_proportion', 
            'high_count_proportion', 'dynamic_range',
            'mean_var_slope', 'mean_var_r_squared', 'poisson_fit_score',
            'quality_score'
        ]
        
        # Build feature matrix
        X = []
        y = []
        
        for profile, result in zip(profiles, results):
            # Extract features
            features = [profile.get(key, 0) for key in feature_keys]
            X.append(features)
            
            # Get best pipeline
            y.append(result['_best_pipeline'])
        
        return np.array(X), y, feature_keys


# =============================================================================
# Convenience Functions
# =============================================================================

def generate_training_data(
    n_datasets: int = 200,
    output_dir: str = 'ml_training_data',
    seed: int = 42
) -> Dict:
    """
    Convenience function to generate training data for ML recommender.
    
    Parameters
    ----------
    n_datasets : int
        Number of datasets to generate (default: 200)
    output_dir : str
        Output directory (default: 'ml_training_data')
    seed : int
        Random seed (default: 42)
    
    Returns
    -------
    dict
        Generation summary
    
    Examples
    --------
    >>> summary = generate_training_data(n_datasets=200)
    >>> print(f"Generated {summary['n_datasets']} datasets")
    >>> print(f"Pipeline distribution: {summary['pipeline_distribution']}")
    """
    # Validate parameters
    validate_positive_integer(n_datasets, 'n_datasets')
    
    if n_datasets < 10:
        raise ValidationError(
            f"n_datasets must be >= 10 for meaningful training data, got {n_datasets}"
        )
    
    # Validate output directory
    if output_dir:
        output_path = validate_directory_path(output_dir, must_exist=False, create_if_missing=True)
    
    generator = SyntheticBenchmarkGenerator(n_datasets=n_datasets, seed=seed)
    summary = generator.generate_benchmarks(output_dir)
    
    print("\n" + "=" * 60)
    print("🦖 RAPTOR Synthetic Benchmark Generation Complete")
    print("=" * 60)
    print(f"\nDatasets generated: {summary['n_datasets']}")
    print(f"Features: {summary['n_features']}")
    print(f"Output directory: {summary['output_dir']}")
    print(f"\nPipeline Distribution (Best Pipeline Counts):")
    print("-" * 40)
    
    for pipeline, count in sorted(
        summary['pipeline_distribution'].items(),
        key=lambda x: -x[1]
    ):
        pct = (count / summary['n_datasets']) * 100
        bar = "█" * int(pct / 2)
        print(f"  {pipeline:15s}: {count:4d} ({pct:5.1f}%) {bar}")
    
    print("\nFiles created:")
    for name, path in summary['files'].items():
        print(f"  • {name}: {path}")
    
    return summary


if __name__ == '__main__':
    print("🦖 RAPTOR Synthetic Benchmark Generator v2.2.0")
    print("=" * 50)
    print("\nGenerates realistic synthetic benchmark data for")
    print("training the ML-based pipeline recommender.")
    print("\nSupported DE Methods:")
    for name, model in PIPELINE_MODELS.items():
        print(f"  • {name}: base F1={model.base_f1:.2f}, "
              f"min_rep={model.min_replicates}")
    print("\nUsage:")
    print("  from raptor.synthetic_benchmarks import generate_training_data")
    print("  summary = generate_training_data(n_datasets=200)")
