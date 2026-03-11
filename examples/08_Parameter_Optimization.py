#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Parameter Optimization (Module 8)

Demonstrates comprehensive parameter optimization for differential expression analysis
using multiple scientifically validated methods.

This is Module 8 of the RAPTOR workflow (Stage 3: DE Analysis):
  M7: Import DE Results → DEResult object
  M8: Parameter Optimization (THIS SCRIPT) → Optimized parameters
  M9: Ensemble Analysis
  M10: Biomarker Discovery

Optimization Methods:
  1. Ground Truth - With validated genes (gold standard)
  2. FDR Control - Statistical FDR estimation (no validation needed)
  3. Stability - Bootstrap-based (no validation needed)
  4. Reproducibility - Independent cohort validation

Input: DE results from Module 7 (de_result.pkl)
Output: results/parameter_optimization/
    - optimization_results.json
    - best_parameters.json
    - convergence_plot.png
    - parameter_space_plot.png

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import argparse
import sys
import json
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONSTANTS - Architecture Compliant (v2.2.0)
# =============================================================================
DEFAULT_OUTPUT_DIR = "results/parameter_optimization"
DEFAULT_STRATEGY = "grid"
DEFAULT_GRID_POINTS = 5
DEFAULT_N_BOOTSTRAP = 100

# Check for dependencies
try:
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: numpy and pandas are required")
    print("Install with: pip install numpy pandas")
    sys.exit(1)

# RAPTOR imports
RAPTOR_AVAILABLE = True
try:
    from raptor.de_import import DEResult
    from raptor.parameter_optimization import (
        optimize_with_ground_truth,
        optimize_with_fdr_control,
        optimize_with_stability,
        optimize_with_reproducibility,
        OptimizationResult
    )
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR not installed. Running in demo mode only.")
    print("Install RAPTOR with: pip install -e .")
    
    # Create dummy classes for demo mode
    class DEResult:
        def __init__(self, results_df, pipeline, parameters, metadata):
            self.results_df = results_df
            self.pipeline = pipeline
            self.parameters = parameters
            self.metadata = metadata
            self.n_genes = len(results_df)
        
        @classmethod
        def load(cls, path):
            """Demo load - generates fake data"""
            return None
    
    class OptimizationResult:
        def __init__(self, best_parameters, best_score, method, strategy):
            self.best_parameters = best_parameters
            self.best_score = best_score
            self.optimization_method = method
            self.search_strategy = strategy
            self.n_iterations = 25
            
        def summary(self):
            return f"""
╔══════════════════════════════════════════════════════════╗
║  🦖 RAPTOR Optimization Results                          ║
╠══════════════════════════════════════════════════════════╣

  Method: {self.optimization_method}
  Strategy: {self.search_strategy}
  
  Best Parameters:
    • alpha: {self.best_parameters['alpha']:.4f}
    • lfc_threshold: {self.best_parameters['lfc_threshold']:.2f}
  
  Best Score: {self.best_score:.4f}
  Iterations: {self.n_iterations}

╚══════════════════════════════════════════════════════════╝
"""


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║    🦖 RAPTOR v2.2.0 - Parameter Optimization (Module 8)      ║
    ║                                                              ║
    ║    Optimize FDR & LFC Thresholds for DE Analysis            ║
    ║    ✅ 4 Methods | 3 Search Strategies                        ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_de_results(n_genes=15000, seed=42):
    """Generate demo DE results for testing."""
    np.random.seed(seed)
    
    # Generate gene IDs
    gene_ids = [f'ENSG{i+1:011d}' for i in range(n_genes)]
    
    # Generate realistic DE results
    base_mean = np.random.gamma(shape=2, scale=100, size=n_genes)
    log2_fold_change = np.random.normal(0, 1.5, n_genes)
    
    # P-values (with some truly DE genes)
    n_de = int(n_genes * 0.1)  # 10% truly DE
    p_values = np.ones(n_genes)
    p_values[:n_de] = np.random.beta(0.5, 10, n_de)
    p_values[n_de:] = np.random.uniform(0.1, 1.0, n_genes - n_de)
    
    # Shuffle
    indices = np.random.permutation(n_genes)
    p_values = p_values[indices]
    
    # Adjusted p-values
    adjusted_p_values = np.minimum(p_values * 10, 1.0)
    
    # Create DataFrame
    df = pd.DataFrame({
        'log2FoldChange': log2_fold_change,
        'pvalue': p_values,
        'padj': adjusted_p_values,
        'baseMean': base_mean
    }, index=gene_ids)
    df.index.name = 'gene_id'
    
    return df


def generate_demo_ground_truth(de_df, n_validated=50, seed=42):
    """Generate demo validated genes."""
    np.random.seed(seed)
    
    # Select truly DE genes (low p-values)
    sig_genes = de_df[de_df['padj'] < 0.01].index.tolist()
    
    if len(sig_genes) >= n_validated:
        validated = np.random.choice(sig_genes, n_validated, replace=False)
    else:
        validated = sig_genes
    
    ground_truth_df = pd.DataFrame({
        'gene_id': validated
    })
    
    return ground_truth_df


def generate_demo_counts(de_df, n_samples=12, seed=42):
    """Generate demo count matrix."""
    np.random.seed(seed)
    
    # Generate counts for each gene
    counts = np.zeros((len(de_df), n_samples))
    
    for i, base_mean in enumerate(de_df['baseMean']):
        # Negative binomial counts
        mu = base_mean
        dispersion = 0.1
        counts[i, :] = np.random.negative_binomial(
            n=1/dispersion,
            p=1/(1 + mu*dispersion),
            size=n_samples
        )
    
    counts_df = pd.DataFrame(
        counts,
        index=de_df.index,
        columns=[f'Sample_{i+1}' for i in range(n_samples)]
    )
    
    return counts_df


def generate_demo_metadata(n_samples=12):
    """Generate demo sample metadata."""
    metadata = pd.DataFrame({
        'sample_id': [f'Sample_{i+1}' for i in range(n_samples)],
        'condition': ['Control'] * (n_samples // 2) + ['Treatment'] * (n_samples // 2),
        'batch': ['Batch1'] * (n_samples // 3) + ['Batch2'] * (n_samples // 3) + ['Batch3'] * (n_samples // 3)
    })
    
    return metadata


def run_ground_truth_optimization(de_df, ground_truth_df, output_dir, 
                                   strategy='grid', grid_points=5, demo=False):
    """
    Run ground truth-based optimization.
    
    Best for: When you have validated DEGs from literature/qPCR
    Requirements: 10-15 genes minimum, 20-30 recommended, 50+ ideal
    """
    print("\n" + "="*70)
    print("  METHOD 1: Ground Truth Optimization")
    print("="*70)
    print("\n  Uses validated genes to find optimal parameters")
    print("  Gold standard method - highest confidence")
    print()
    
    output_path = Path(output_dir) / 'ground_truth'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • DE results: {len(de_df):,} genes")
        print(f"     • Validated genes: {len(ground_truth_df)} genes")
        print(f"     • Strategy: {strategy}")
        print(f"     • Grid points: {grid_points} ({grid_points**2} combinations)")
        print()
        
        # Simulate optimization
        best_params = {
            'alpha': 0.05,
            'lfc_threshold': 0.5
        }
        best_score = 0.87
        
        result = OptimizationResult(
            best_parameters=best_params,
            best_score=best_score,
            method='ground_truth',
            strategy=strategy
        )
        
        print(result.summary())
        
        # Save results
        with open(output_path / 'optimization_results.json', 'w') as f:
            json.dump({
                'method': 'ground_truth',
                'best_parameters': best_params,
                'best_score': best_score,
                'strategy': strategy,
                'n_validated_genes': len(ground_truth_df),
                'timestamp': datetime.now().isoformat()
            }, f, indent=2)
        
        return result
    
    # Real optimization
    print("  🚀 Running optimization with RAPTOR...")
    print(f"     • DE results: {len(de_df):,} genes")
    print(f"     • Validated genes: {len(ground_truth_df)} genes")
    print(f"     • Strategy: {strategy}")
    
    try:
        result = optimize_with_ground_truth(
            de_result=de_df,
            ground_truth=ground_truth_df,
            strategy=strategy,
            grid_points=grid_points,
            output_dir=str(output_path)
        )
        
        print(result.summary())
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during optimization: {e}")
        return None


def run_fdr_control_optimization(de_df, output_dir, strategy='grid', 
                                 grid_points=5, demo=False):
    """
    Run FDR control optimization.
    
    Best for: When you DON'T have validated genes
    Method: Statistical FDR estimation (Storey & Tibshirani 2003)
    Requirements: 1000+ genes for reliable FDR estimation
    """
    print("\n" + "="*70)
    print("  METHOD 2: FDR Control Optimization")
    print("="*70)
    print("\n  Statistical FDR estimation - NO validation needed!")
    print("  Based on: Storey & Tibshirani (2003) PNAS")
    print()
    
    output_path = Path(output_dir) / 'fdr_control'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • DE results: {len(de_df):,} genes")
        print(f"     • Strategy: {strategy}")
        print(f"     • Grid points: {grid_points}")
        print()
        
        # Simulate optimization
        best_params = {
            'alpha': 0.03,
            'lfc_threshold': 0.25
        }
        best_score = 0.92
        
        result = OptimizationResult(
            best_parameters=best_params,
            best_score=best_score,
            method='fdr_control',
            strategy=strategy
        )
        
        print(result.summary())
        
        # Save results
        with open(output_path / 'optimization_results.json', 'w') as f:
            json.dump({
                'method': 'fdr_control',
                'best_parameters': best_params,
                'best_score': best_score,
                'strategy': strategy,
                'timestamp': datetime.now().isoformat()
            }, f, indent=2)
        
        return result
    
    # Real optimization
    print("  🚀 Running optimization with RAPTOR...")
    print(f"     • DE results: {len(de_df):,} genes")
    
    try:
        result = optimize_with_fdr_control(
            de_result=de_df,
            strategy=strategy,
            grid_points=grid_points,
            output_dir=str(output_path)
        )
        
        print(result.summary())
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during optimization: {e}")
        return None


def run_stability_optimization(de_df, counts_df, metadata_df, output_dir,
                               n_bootstrap=100, strategy='grid', 
                               grid_points=5, demo=False):
    """
    Run stability-based optimization.
    
    Best for: When you want robust parameter selection
    Method: Bootstrap stability (Meinshausen & Bühlmann 2010)
    Requirements: Original count matrix and sample metadata
    """
    print("\n" + "="*70)
    print("  METHOD 3: Stability Optimization")
    print("="*70)
    print("\n  Bootstrap-based stability - NO validation needed!")
    print("  Based on: Meinshausen & Bühlmann (2010) JRSS-B")
    print()
    
    output_path = Path(output_dir) / 'stability'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • DE results: {len(de_df):,} genes")
        print(f"     • Count matrix: {counts_df.shape}")
        print(f"     • Samples: {len(metadata_df)}")
        print(f"     • Bootstrap iterations: {n_bootstrap}")
        print(f"     • Strategy: {strategy}")
        print()
        
        # Simulate optimization
        best_params = {
            'alpha': 0.04,
            'lfc_threshold': 0.75
        }
        best_score = 0.89
        
        result = OptimizationResult(
            best_parameters=best_params,
            best_score=best_score,
            method='stability',
            strategy=strategy
        )
        
        print(result.summary())
        
        # Save results
        with open(output_path / 'optimization_results.json', 'w') as f:
            json.dump({
                'method': 'stability',
                'best_parameters': best_params,
                'best_score': best_score,
                'strategy': strategy,
                'n_bootstrap': n_bootstrap,
                'timestamp': datetime.now().isoformat()
            }, f, indent=2)
        
        return result
    
    # Real optimization
    print("  🚀 Running optimization with RAPTOR...")
    print(f"     • DE results: {len(de_df):,} genes")
    print(f"     • Bootstrap iterations: {n_bootstrap}")
    
    try:
        result = optimize_with_stability(
            de_result=de_df,
            counts=counts_df,
            metadata=metadata_df,
            n_bootstrap=n_bootstrap,
            strategy=strategy,
            grid_points=grid_points,
            output_dir=str(output_path)
        )
        
        print(result.summary())
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during optimization: {e}")
        return None


def run_reproducibility_optimization(cohort1_df, cohort2_df, output_dir,
                                     strategy='grid', grid_points=5, demo=False):
    """
    Run reproducibility-based optimization.
    
    Best for: When you have independent cohorts
    Method: Cross-cohort validation
    Requirements: Two independent DE result datasets
    """
    print("\n" + "="*70)
    print("  METHOD 4: Reproducibility Optimization")
    print("="*70)
    print("\n  Independent cohort validation")
    print("  Gold standard for parameter robustness")
    print()
    
    output_path = Path(output_dir) / 'reproducibility'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • Cohort 1: {len(cohort1_df):,} genes")
        print(f"     • Cohort 2: {len(cohort2_df):,} genes")
        print(f"     • Strategy: {strategy}")
        print()
        
        # Simulate optimization
        best_params = {
            'alpha': 0.045,
            'lfc_threshold': 0.6
        }
        best_score = 0.85
        
        result = OptimizationResult(
            best_parameters=best_params,
            best_score=best_score,
            method='reproducibility',
            strategy=strategy
        )
        
        print(result.summary())
        
        # Save results
        with open(output_path / 'optimization_results.json', 'w') as f:
            json.dump({
                'method': 'reproducibility',
                'best_parameters': best_params,
                'best_score': best_score,
                'strategy': strategy,
                'timestamp': datetime.now().isoformat()
            }, f, indent=2)
        
        return result
    
    # Real optimization
    print("  🚀 Running optimization with RAPTOR...")
    
    try:
        result = optimize_with_reproducibility(
            de_result_cohort1=cohort1_df,
            de_result_cohort2=cohort2_df,
            strategy=strategy,
            grid_points=grid_points,
            output_dir=str(output_path)
        )
        
        print(result.summary())
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during optimization: {e}")
        return None


def compare_methods(results_dict, output_dir):
    """Compare optimization results from different methods."""
    print("\n" + "="*70)
    print("  COMPARISON: All Methods")
    print("="*70)
    print()
    
    comparison_data = []
    
    for method_name, result in results_dict.items():
        if result is None:
            continue
        
        comparison_data.append({
            'Method': method_name,
            'Alpha': result.best_parameters['alpha'],
            'LFC Threshold': result.best_parameters['lfc_threshold'],
            'Score': result.best_score,
            'Strategy': result.search_strategy
        })
    
    if not comparison_data:
        print("  No results to compare")
        return
    
    # Create comparison table
    comparison_df = pd.DataFrame(comparison_data)
    
    print("  Parameter Comparison:")
    print("  " + "-"*66)
    print(f"  {'Method':<20} {'Alpha':<10} {'LFC':<10} {'Score':<10} {'Strategy':<10}")
    print("  " + "-"*66)
    
    for _, row in comparison_df.iterrows():
        print(f"  {row['Method']:<20} {row['Alpha']:<10.4f} {row['LFC Threshold']:<10.2f} "
              f"{row['Score']:<10.4f} {row['Strategy']:<10}")
    
    print("  " + "-"*66)
    print()
    
    # Find consensus
    alphas = [r['Alpha'] for r in comparison_data]
    lfcs = [r['LFC Threshold'] for r in comparison_data]
    
    print("  📊 Consensus Parameters:")
    print(f"     • Alpha range: {min(alphas):.4f} - {max(alphas):.4f}")
    print(f"     • Alpha median: {np.median(alphas):.4f}")
    print(f"     • LFC range: {min(lfcs):.2f} - {max(lfcs):.2f}")
    print(f"     • LFC median: {np.median(lfcs):.2f}")
    print()
    
    # Save comparison
    output_path = Path(output_dir)
    comparison_df.to_csv(output_path / 'method_comparison.csv', index=False)
    
    consensus = {
        'alpha_median': float(np.median(alphas)),
        'alpha_range': [float(min(alphas)), float(max(alphas))],
        'lfc_median': float(np.median(lfcs)),
        'lfc_range': [float(min(lfcs)), float(max(lfcs))],
        'methods_compared': len(comparison_data),
        'timestamp': datetime.now().isoformat()
    }
    
    with open(output_path / 'consensus_parameters.json', 'w') as f:
        json.dump(consensus, f, indent=2)
    
    print(f"  ✓ Saved: {output_path / 'method_comparison.csv'}")
    print(f"  ✓ Saved: {output_path / 'consensus_parameters.json'}")


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Parameter Optimization (Module 8)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (no data needed)
  python 08_Parameter_Optimization.py --demo
  
  # Ground truth optimization (with validated genes)
  python 08_Parameter_Optimization.py \\
      --de-result results/de_imported/de_result.pkl \\
      --method ground_truth \\
      --ground-truth data/validated_genes.csv
  
  # FDR control (no validation needed)
  python 08_Parameter_Optimization.py \\
      --de-result results/de_imported/de_result.pkl \\
      --method fdr_control
  
  # Stability optimization (needs counts)
  python 08_Parameter_Optimization.py \\
      --de-result results/de_imported/de_result.pkl \\
      --method stability \\
      --counts data/counts.csv \\
      --metadata data/metadata.csv
  
  # All methods
  python 08_Parameter_Optimization.py \\
      --de-result results/de_imported/de_result.pkl \\
      --method all \\
      --ground-truth data/validated_genes.csv \\
      --counts data/counts.csv \\
      --metadata data/metadata.csv

Optimization Methods:
  
  1. Ground Truth (Validated Genes)
     • Best for: When you have qPCR/literature validation
     • Requires: 10-50 validated DE genes
     • Metric: F1 score
  
  2. FDR Control (Statistical)
     • Best for: No validation available
     • Requires: Nothing! (uses statistics)
     • Metric: FDR control accuracy
  
  3. Stability (Bootstrap)
     • Best for: Robust parameter selection
     • Requires: Count matrix + metadata
     • Metric: Selection stability
  
  4. Reproducibility (Independent Cohorts)
     • Best for: Cross-study validation
     • Requires: Two independent DE results
     • Metric: Cross-cohort overlap

Search Strategies:
  - grid: Exhaustive search (guaranteed optimum)
  - random: Monte Carlo sampling (faster)
  - differential_evolution: Population-based global optimization

Output Files:
  results/parameter_optimization/
    ├── ground_truth/optimization_results.json
    ├── fdr_control/optimization_results.json
    ├── stability/optimization_results.json
    ├── reproducibility/optimization_results.json
    ├── method_comparison.csv
    └── consensus_parameters.json
        """
    )
    
    parser.add_argument('--de-result', '-d',
                       help='DE result file from Module 7 (.pkl)')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--method', '-m',
                       choices=['ground_truth', 'fdr_control', 'stability', 
                               'reproducibility', 'all'],
                       default='fdr_control',
                       help='Optimization method (default: fdr_control)')
    parser.add_argument('--strategy', '-s',
                       choices=['grid', 'random', 'differential_evolution'],
                       default=DEFAULT_STRATEGY,
                       help=f'Search strategy (default: {DEFAULT_STRATEGY})')
    parser.add_argument('--grid-points', type=int, default=DEFAULT_GRID_POINTS,
                       help=f'Grid points per parameter (default: {DEFAULT_GRID_POINTS})')
    parser.add_argument('--ground-truth',
                       help='Validated genes CSV (for ground_truth method)')
    parser.add_argument('--counts',
                       help='Count matrix CSV (for stability method)')
    parser.add_argument('--metadata',
                       help='Sample metadata CSV (for stability method)')
    parser.add_argument('--cohort2',
                       help='Second cohort DE results (for reproducibility method)')
    parser.add_argument('--n-bootstrap', type=int, default=DEFAULT_N_BOOTSTRAP,
                       help=f'Bootstrap iterations (default: {DEFAULT_N_BOOTSTRAP})')
    parser.add_argument('--demo', action='store_true',
                       help='Run in demo mode with simulated data')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs for real run
    if not args.demo and not args.de_result:
        print("ERROR: --de-result is required (or use --demo)")
        parser.print_help()
        sys.exit(1)
    
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # =========================================================================
    # Load or Generate Data
    # =========================================================================
    
    if args.demo or not RAPTOR_AVAILABLE:
        print("\n🎮 DEMO MODE - Generating simulated data...")
        print("─" * 60)
        
        # Generate demo DE results
        print("  Generating DE results (15,000 genes)...")
        de_df = generate_demo_de_results(n_genes=15000)
        print(f"  ✓ Generated: {len(de_df):,} genes")
        
        # Generate ground truth
        print("  Generating validated genes (50 genes)...")
        ground_truth_df = generate_demo_ground_truth(de_df, n_validated=50)
        print(f"  ✓ Generated: {len(ground_truth_df)} validated genes")
        
        # Generate counts and metadata
        print("  Generating count matrix and metadata...")
        counts_df = generate_demo_counts(de_df, n_samples=12)
        metadata_df = generate_demo_metadata(n_samples=12)
        print(f"  ✓ Generated: {counts_df.shape} count matrix")
        print(f"  ✓ Generated: {len(metadata_df)} samples")
        
        # Generate second cohort for reproducibility
        cohort2_df = generate_demo_de_results(n_genes=15000, seed=123)
        print(f"  ✓ Generated: Cohort 2 ({len(cohort2_df):,} genes)")
        
    else:
        print("\n🚀 LOADING REAL DATA")
        print("─" * 60)
        
        # Load DE results
        print(f"  Loading DE results from: {args.de_result}")
        try:
            de_result_obj = DEResult.load(args.de_result)
            de_df = de_result_obj.results_df
            print(f"  ✓ Loaded: {len(de_df):,} genes")
        except Exception as e:
            print(f"  ❌ Error loading DE results: {e}")
            sys.exit(1)
        
        # Load ground truth if provided
        ground_truth_df = None
        if args.ground_truth:
            print(f"  Loading validated genes from: {args.ground_truth}")
            try:
                ground_truth_df = pd.read_csv(args.ground_truth)
                print(f"  ✓ Loaded: {len(ground_truth_df)} validated genes")
            except Exception as e:
                print(f"  ❌ Error loading ground truth: {e}")
        
        # Load counts if provided
        counts_df = None
        if args.counts:
            print(f"  Loading count matrix from: {args.counts}")
            try:
                counts_df = pd.read_csv(args.counts, index_col=0)
                print(f"  ✓ Loaded: {counts_df.shape} count matrix")
            except Exception as e:
                print(f"  ❌ Error loading counts: {e}")
        
        # Load metadata if provided
        metadata_df = None
        if args.metadata:
            print(f"  Loading metadata from: {args.metadata}")
            try:
                metadata_df = pd.read_csv(args.metadata)
                print(f"  ✓ Loaded: {len(metadata_df)} samples")
            except Exception as e:
                print(f"  ❌ Error loading metadata: {e}")
        
        # Load cohort 2 if provided
        cohort2_df = None
        if args.cohort2:
            print(f"  Loading cohort 2 from: {args.cohort2}")
            try:
                cohort2_result = DEResult.load(args.cohort2)
                cohort2_df = cohort2_result.results_df
                print(f"  ✓ Loaded: {len(cohort2_df):,} genes")
            except Exception as e:
                print(f"  ❌ Error loading cohort 2: {e}")
    
    # =========================================================================
    # Run Optimization
    # =========================================================================
    
    results = {}
    
    if args.method == 'ground_truth' or args.method == 'all':
        if ground_truth_df is not None:
            result = run_ground_truth_optimization(
                de_df, ground_truth_df, args.output,
                strategy=args.strategy, grid_points=args.grid_points,
                demo=args.demo
            )
            if result:
                results['Ground Truth'] = result
        else:
            print("\n⚠️  Ground truth method requires --ground-truth")
    
    if args.method == 'fdr_control' or args.method == 'all':
        result = run_fdr_control_optimization(
            de_df, args.output,
            strategy=args.strategy, grid_points=args.grid_points,
            demo=args.demo
        )
        if result:
            results['FDR Control'] = result
    
    if args.method == 'stability' or args.method == 'all':
        if counts_df is not None and metadata_df is not None:
            result = run_stability_optimization(
                de_df, counts_df, metadata_df, args.output,
                n_bootstrap=args.n_bootstrap, strategy=args.strategy,
                grid_points=args.grid_points, demo=args.demo
            )
            if result:
                results['Stability'] = result
        else:
            print("\n⚠️  Stability method requires --counts and --metadata")
    
    if args.method == 'reproducibility' or args.method == 'all':
        if cohort2_df is not None:
            result = run_reproducibility_optimization(
                de_df, cohort2_df, args.output,
                strategy=args.strategy, grid_points=args.grid_points,
                demo=args.demo
            )
            if result:
                results['Reproducibility'] = result
        else:
            print("\n⚠️  Reproducibility method requires --cohort2")
    
    # =========================================================================
    # Compare Methods
    # =========================================================================
    
    if len(results) > 1:
        compare_methods(results, args.output)
    
    # =========================================================================
    # Final Summary
    # =========================================================================
    
    print("\n" + "="*70)
    print("  ✅ MODULE 8 (PARAMETER OPTIMIZATION) COMPLETE!")
    print("="*70)
    
    print(f"\n  📂 Output Directory: {output_path}")
    
    print(f"\n  📊 Optimization Results:")
    for method_name, result in results.items():
        print(f"     • {method_name}:")
        print(f"       - Alpha: {result.best_parameters['alpha']:.4f}")
        print(f"       - LFC: {result.best_parameters['lfc_threshold']:.2f}")
        print(f"       - Score: {result.best_score:.4f}")
    
    if len(results) > 0:
        print(f"\n  📄 Output Files:")
        for method_name in results.keys():
            method_dir = method_name.lower().replace(' ', '_')
            print(f"     • {method_dir}/optimization_results.json")
        
        if len(results) > 1:
            print(f"     • method_comparison.csv")
            print(f"     • consensus_parameters.json")
    
    print(f"\n  🔜 Next Steps:")
    print(f"\n     Module 9 - Ensemble Analysis:")
    print(f"     python 09_Ensemble_Analysis.py \\")
    print(f"         --de-results results/de_imported/*.pkl")
    print(f"")
    print(f"     Apply optimized parameters in R:")
    print(f"     # Use alpha={results[list(results.keys())[0]].best_parameters['alpha']:.4f} "
          f"and LFC={results[list(results.keys())[0]].best_parameters['lfc_threshold']:.2f}")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
