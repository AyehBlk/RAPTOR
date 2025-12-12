"""
Synthetic Benchmark Data Generator

Generates realistic synthetic benchmark data for training and testing
the ML-based pipeline recommender.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List
import logging

logger = logging.getLogger(__name__)


class SyntheticBenchmarkGenerator:
    """
    Generate synthetic but realistic benchmark data.
    
    Creates data profiles and corresponding pipeline performance results
    based on known patterns and relationships in RNA-seq analysis.
    
    Parameters
    ----------
    n_datasets : int
        Number of synthetic datasets to generate
    seed : int, optional
        Random seed
    
    Examples
    --------
    >>> generator = SyntheticBenchmarkGenerator(n_datasets=100, seed=42)
    >>> generator.generate_benchmarks('training_data/')
    """
    
    def __init__(self, n_datasets: int = 100, seed: int = 42):
        """Initialize generator."""
        self.n_datasets = n_datasets
        self.seed = seed
        np.random.seed(seed)
        
        logger.info(f"Initialized generator for {n_datasets} datasets")
    
    def generate_benchmarks(self, output_dir: str) -> Dict:
        """
        Generate complete benchmark dataset.
        
        Parameters
        ----------
        output_dir : str
            Output directory for benchmark data
        
        Returns
        -------
        dict
            Summary of generated data
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Generating {self.n_datasets} synthetic benchmarks...")
        
        datasets_generated = 0
        pipeline_distribution = {i: 0 for i in range(1, 9)}
        
        for i in range(self.n_datasets):
            dataset_dir = output_path / f'dataset_{i:04d}'
            dataset_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate data profile
            profile = self._generate_profile()
            
            # Generate benchmark results based on profile
            benchmark_results = self._generate_benchmark_results(profile)
            
            # Determine best pipeline
            best_pipeline = max(
                benchmark_results.items(),
                key=lambda x: x[1].get('metrics', {}).get('f1_score', 0)
            )[0]
            
            pipeline_distribution[int(best_pipeline)] += 1
            
            # Save files
            with open(dataset_dir / 'data_profile.json', 'w') as f:
                json.dump(profile, f, indent=2)
            
            with open(dataset_dir / 'benchmark_results.json', 'w') as f:
                json.dump(benchmark_results, f, indent=2)
            
            datasets_generated += 1
            
            if (i + 1) % 20 == 0:
                logger.info(f"Generated {i + 1}/{self.n_datasets} datasets")
        
        summary = {
            'n_datasets': datasets_generated,
            'output_dir': str(output_path),
            'pipeline_distribution': pipeline_distribution,
            'seed': self.seed
        }
        
        # Save summary
        with open(output_path / 'generation_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Generation complete: {datasets_generated} datasets")
        logger.info(f"Pipeline distribution: {pipeline_distribution}")
        
        return summary
    
    def _generate_profile(self) -> Dict:
        """Generate realistic data profile."""
        
        # Sample size (biased towards common sizes)
        n_samples = np.random.choice(
            [3, 4, 6, 8, 12, 16, 24],
            p=[0.15, 0.15, 0.25, 0.20, 0.15, 0.05, 0.05]
        )
        
        # Gene count
        n_genes = np.random.choice(
            [15000, 20000, 25000, 30000],
            p=[0.2, 0.4, 0.3, 0.1]
        )
        
        # Library size (realistic distribution)
        mean_lib_size = np.random.lognormal(mean=np.log(20e6), sigma=0.3)
        lib_size_cv = np.random.beta(2, 8)  # Usually low CV
        
        # Zero inflation
        zero_pct = np.random.beta(6, 2) * 70  # 0-70%
        
        # Biological variation
        bcv = np.random.beta(2, 5) * 0.8  # 0-0.8
        
        # Sequencing depth
        reads_per_gene = mean_lib_size / n_genes
        if reads_per_gene < 500:
            depth_category = 'low'
        elif reads_per_gene < 1500:
            depth_category = 'medium'
        elif reads_per_gene < 3000:
            depth_category = 'high'
        else:
            depth_category = 'very_high'
        
        # Data complexity
        complexity_score = 50 + np.random.normal(0, 20)
        complexity_score = max(0, min(100, complexity_score))
        
        profile = {
            'design': {
                'n_samples': int(n_samples),
                'n_genes': int(n_genes),
                'n_conditions': 2,
                'samples_per_condition': int(n_samples / 2),
                'is_paired': bool(np.random.random() > 0.7)
            },
            'library_stats': {
                'mean': float(mean_lib_size),
                'median': float(mean_lib_size * 0.95),
                'cv': float(lib_size_cv),
                'range': float(mean_lib_size * lib_size_cv * 4),
                'skewness': float(np.random.normal(0, 0.5))
            },
            'count_distribution': {
                'zero_pct': float(zero_pct),
                'low_count_pct': float(zero_pct + np.random.uniform(10, 20)),
                'mean': float(mean_lib_size / n_genes),
                'median': float(mean_lib_size / n_genes * 0.7),
                'variance': float((mean_lib_size / n_genes) ** 2 * (1 + bcv ** 2))
            },
            'expression_distribution': {
                'high_expr_genes': int(n_genes * 0.1),
                'medium_expr_genes': int(n_genes * 0.3),
                'low_expr_genes': int(n_genes * 0.6),
                'dynamic_range': float(np.random.uniform(6, 12))
            },
            'biological_variation': {
                'bcv': float(bcv),
                'dispersion_mean': float(bcv ** 2),
                'dispersion_trend': float(np.random.uniform(-0.5, 0.5)),
                'outlier_genes': int(n_genes * np.random.uniform(0.01, 0.05))
            },
            'sequencing': {
                'total_reads': float(mean_lib_size * n_samples),
                'reads_per_gene': float(reads_per_gene),
                'depth_category': depth_category
            },
            'complexity': {
                'score': float(complexity_score),
                'noise_level': float(bcv * 2),
                'signal_strength': float(1.0 - bcv)
            },
            'summary': {
                'difficulty': self._categorize_difficulty(bcv, n_samples, zero_pct),
                'recommended_approach': 'standard'  # Will be determined by ML
            }
        }
        
        return profile
    
    def _generate_benchmark_results(self, profile: Dict) -> Dict:
        """
        Generate realistic benchmark results based on data characteristics.
        
        Different pipelines perform better under different conditions.
        """
        results = {}
        
        # Extract key characteristics
        n_samples = profile['design']['n_samples']
        bcv = profile['biological_variation']['bcv']
        depth = profile['sequencing']['depth_category']
        zero_pct = profile['count_distribution']['zero_pct']
        
        # Base performance for each pipeline
        base_performance = {
            1: 0.85,  # STAR-RSEM-DESeq2: generally excellent
            2: 0.78,  # HISAT2-StringTie-Ballgown: good
            3: 0.82,  # Salmon-edgeR: fast and accurate
            4: 0.80,  # Kallisto-Sleuth: very fast
            5: 0.83,  # STAR-HTSeq-limma: robust
            6: 0.75,  # NOISeq: for no replicates
            7: 0.77,  # Bowtie2-RSEM-EBSeq: older
            8: 0.72   # HISAT2-Cufflinks: less popular now
        }
        
        for pipeline_id in range(1, 9):
            base_f1 = base_performance[pipeline_id]
            
            # Adjust based on data characteristics
            f1_score = base_f1
            
            # Sample size effects
            if n_samples < 3 and pipeline_id == 6:
                f1_score += 0.10  # NOISeq excels here
            elif n_samples < 3 and pipeline_id in [1, 3, 5]:
                f1_score -= 0.15  # Others struggle
            elif n_samples >= 6 and pipeline_id in [1, 3, 4, 5]:
                f1_score += 0.05  # Benefit from replicates
            
            # BCV effects
            if bcv > 0.5 and pipeline_id in [1, 5]:
                f1_score += 0.05  # Handle variation well
            elif bcv > 0.5 and pipeline_id in [4, 8]:
                f1_score -= 0.08  # Struggle with high variation
            
            # Depth effects
            if depth == 'low' and pipeline_id in [3, 4]:
                f1_score += 0.05  # Efficient pseudo-aligners
            elif depth == 'low' and pipeline_id in [1, 2]:
                f1_score -= 0.05  # Heavy aligners less efficient
            
            # Zero inflation effects
            if zero_pct > 60 and pipeline_id in [3, 4]:
                f1_score -= 0.05
            
            # Add realistic noise
            f1_score += np.random.normal(0, 0.03)
            f1_score = max(0.4, min(0.95, f1_score))
            
            # Generate other metrics
            accuracy = f1_score + np.random.uniform(-0.05, 0.05)
            precision = f1_score + np.random.uniform(-0.03, 0.03)
            recall = f1_score + np.random.uniform(-0.03, 0.03)
            
            # Runtime (varies by pipeline and data size)
            base_runtimes = {
                1: 3600, 2: 2400, 3: 600, 4: 300,
                5: 3000, 6: 2800, 7: 3200, 8: 2600
            }
            
            runtime = base_runtimes[pipeline_id]
            runtime *= (n_samples / 6)  # Scale with samples
            runtime *= (1 + np.random.uniform(-0.2, 0.2))  # Add variation
            
            # Memory usage
            base_memory = {
                1: 32, 2: 24, 3: 16, 4: 8,
                5: 28, 6: 26, 7: 30, 8: 22
            }
            
            memory_gb = base_memory[pipeline_id]
            memory_gb *= (1 + np.random.uniform(-0.15, 0.15))
            
            results[pipeline_id] = {
                'status': 'success',
                'pipeline': self._get_pipeline_name(pipeline_id),
                'metrics': {
                    'f1_score': float(f1_score),
                    'accuracy': float(max(0, min(1, accuracy))),
                    'precision': float(max(0, min(1, precision))),
                    'recall': float(max(0, min(1, recall)))
                },
                'runtime': float(runtime),
                'memory_gb': float(memory_gb),
                'n_de_genes_detected': int(np.random.poisson(500))
            }
        
        return results
    
    @staticmethod
    def _categorize_difficulty(bcv: float, n_samples: int, zero_pct: float) -> str:
        """Categorize dataset difficulty."""
        score = 0
        
        if bcv > 0.5:
            score += 2
        elif bcv > 0.3:
            score += 1
        
        if n_samples < 3:
            score += 2
        elif n_samples < 6:
            score += 1
        
        if zero_pct > 70:
            score += 2
        elif zero_pct > 60:
            score += 1
        
        if score >= 5:
            return 'very_hard'
        elif score >= 3:
            return 'hard'
        elif score >= 1:
            return 'medium'
        else:
            return 'easy'
    
    @staticmethod
    def _get_pipeline_name(pipeline_id: int) -> str:
        """Get pipeline name."""
        names = {
            1: "STAR-RSEM-DESeq2",
            2: "HISAT2-StringTie-Ballgown",
            3: "Salmon-edgeR",
            4: "Kallisto-Sleuth",
            5: "STAR-HTSeq-limma-voom",
            6: "STAR-featureCounts-NOISeq",
            7: "Bowtie2-RSEM-EBSeq",
            8: "HISAT2-Cufflinks-Cuffdiff"
        }
        return names.get(pipeline_id, f"Pipeline {pipeline_id}")


def generate_training_data(
    n_datasets: int = 200,
    output_dir: str = 'ml_training_data',
    seed: int = 42
) -> Dict:
    """
    Convenience function to generate training data.
    
    Parameters
    ----------
    n_datasets : int
        Number of datasets to generate
    output_dir : str
        Output directory
    seed : int
        Random seed
    
    Returns
    -------
    dict
        Generation summary
    
    Examples
    --------
    >>> summary = generate_training_data(n_datasets=200)
    >>> print(f"Generated {summary['n_datasets']} datasets")
    """
    generator = SyntheticBenchmarkGenerator(n_datasets=n_datasets, seed=seed)
    summary = generator.generate_benchmarks(output_dir)
    
    print("\n=== Synthetic Benchmark Generation ===")
    print(f"Datasets: {summary['n_datasets']}")
    print(f"Output: {summary['output_dir']}")
    print(f"\nPipeline Distribution:")
    for pid, count in summary['pipeline_distribution'].items():
        pname = SyntheticBenchmarkGenerator._get_pipeline_name(pid)
        print(f"  Pipeline {pid} ({pname}): {count}")
    
    return summary


if __name__ == '__main__':
    print("RAPTOR Synthetic Benchmark Generator")
    print("====================================")
    print("\nGenerates realistic synthetic benchmark data for ML training.")
    print("\nUsage:")
    print("  from raptor.synthetic_benchmarks import generate_training_data")
    print("  summary = generate_training_data(n_datasets=200)")
