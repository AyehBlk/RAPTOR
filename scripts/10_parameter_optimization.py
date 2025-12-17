#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Example Script: Parameter Optimization

Demonstrates automated parameter tuning with:
- Adaptive optimization based on data characteristics
- Grid search optimization
- Random search with intelligent sampling
- Bayesian optimization (optional)
- Cross-validation based evaluation
- Parameter importance analysis

Author: Ayeh Bolouki
License: MIT
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

# Check for dependencies
try:
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: numpy and pandas are required")
    print("Install with: pip install numpy pandas")
    sys.exit(1)

# RAPTOR imports with fallback
RAPTOR_AVAILABLE = True
try:
    from raptor.parameter_optimization import ParameterOptimizer, ParameterSpace, optimize_pipeline_parameters
    from raptor.profiler import RNAseqDataProfiler
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - Parameter Optimization             ║
    ║                                                              ║
    ║   Automated Parameter Tuning for RNA-seq Analysis           ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_profile():
    """Generate demonstration data profile."""
    np.random.seed(42)
    
    return {
        'design': {
            'n_samples': 12,
            'n_genes': 20000,
            'n_conditions': 2,
            'samples_per_condition': 6,
            'is_paired': False,
            'min_replicates': 6
        },
        'library_stats': {
            'mean': 25000000.0,
            'cv': 0.18
        },
        'count_distribution': {
            'zero_pct': 42.5,
            'low_count_pct': 58.2,
            'mean': 1250.0
        },
        'biological_variation': {
            'bcv': 0.35,
            'dispersion_mean': 0.12
        },
        'sequencing': {
            'depth_category': 'high',
            'mean_depth': 25000000.0
        },
        'complexity': {
            'score': 65.0
        },
        'summary': {
            'difficulty': 'moderate'
        }
    }


def simple_grid_search(profile, n_iterations=20):
    """Simple grid search implementation for demo mode."""
    np.random.seed(42)
    
    # Define parameter grid
    param_grid = {
        'min_count_threshold': [5, 10, 15, 20],
        'min_samples_expressing': [1, 2, 3],
        'normalization_method': ['TMM', 'RLE', 'upperquartile'],
        'fdr_threshold': [0.01, 0.05, 0.1],
        'log2fc_threshold': [0.0, 0.5, 1.0],
        'dispersion_type': ['parametric', 'local']
    }
    
    # Data characteristics affect optimal parameters
    n_samples = profile['design']['n_samples']
    bcv = profile['biological_variation']['bcv']
    zero_pct = profile['count_distribution']['zero_pct']
    
    # Simulate optimization
    optimization_history = []
    best_score = -np.inf
    best_params = None
    
    for iteration in range(n_iterations):
        # Generate random parameter combination
        params = {
            'min_count_threshold': np.random.choice(param_grid['min_count_threshold']),
            'min_samples_expressing': np.random.choice(param_grid['min_samples_expressing']),
            'normalization_method': np.random.choice(param_grid['normalization_method']),
            'fdr_threshold': np.random.choice(param_grid['fdr_threshold']),
            'log2fc_threshold': np.random.choice(param_grid['log2fc_threshold']),
            'dispersion_type': np.random.choice(param_grid['dispersion_type'])
        }
        
        # Calculate simulated score based on data characteristics
        score = 0.5  # Base score
        
        # Higher min_count is better for high zero inflation
        if zero_pct > 50 and params['min_count_threshold'] >= 15:
            score += 0.1
        elif zero_pct < 40 and params['min_count_threshold'] <= 10:
            score += 0.08
        
        # More samples expressing for lower sample counts
        if n_samples < 6 and params['min_samples_expressing'] == 1:
            score += 0.05
        elif n_samples >= 10 and params['min_samples_expressing'] >= 2:
            score += 0.08
        
        # TMM is generally good, RLE for high CV
        if params['normalization_method'] == 'TMM':
            score += 0.05
        if profile['library_stats']['cv'] > 0.3 and params['normalization_method'] == 'RLE':
            score += 0.05
        
        # Parametric dispersion for moderate BCV
        if bcv < 0.4 and params['dispersion_type'] == 'parametric':
            score += 0.05
        elif bcv >= 0.4 and params['dispersion_type'] == 'local':
            score += 0.05
        
        # Add noise
        score += np.random.normal(0, 0.03)
        score = max(0.3, min(0.95, score))
        
        optimization_history.append({
            'iteration': iteration,
            'parameters': params.copy(),
            'score': score
        })
        
        if score > best_score:
            best_score = score
            best_params = params.copy()
    
    # Calculate parameter importance
    importance = {}
    for param_name in param_grid.keys():
        param_values = []
        scores = []
        
        for entry in optimization_history:
            param_values.append(entry['parameters'][param_name])
            scores.append(entry['score'])
        
        # For categorical, use variance of scores per category
        unique_values = set(param_values)
        category_scores = {v: [] for v in unique_values}
        for val, score in zip(param_values, scores):
            category_scores[val].append(score)
        
        variances = [np.var(s) for s in category_scores.values() if len(s) > 1]
        importance[param_name] = np.mean(variances) if variances else 0.1
    
    # Normalize importance
    total = sum(importance.values())
    if total > 0:
        importance = {k: v / total for k, v in importance.items()}
    
    return {
        'best_parameters': best_params,
        'optimization_score': best_score,
        'parameter_importance': dict(sorted(importance.items(), key=lambda x: x[1], reverse=True)),
        'optimization_history': optimization_history,
        'strategy': 'grid',
        'n_iterations': n_iterations
    }


def display_optimization_results(results, profile):
    """Display optimization results with formatting."""
    
    print("\n" + "="*70)
    print("  🦖 PARAMETER OPTIMIZATION RESULTS")
    print("="*70)
    
    # Strategy and score
    print(f"\n  📊 Optimization Summary:")
    print(f"     Strategy:    {results['strategy']}")
    print(f"     Iterations:  {results['n_iterations']}")
    print(f"     Best Score:  {results['optimization_score']:.3f}")
    
    # Best parameters
    print(f"\n  🎯 Optimal Parameters:")
    print("  " + "-"*66)
    
    for param, value in results['best_parameters'].items():
        print(f"     {param:30s}: {value}")
    
    # Parameter importance
    print(f"\n  📈 Parameter Importance:")
    print("  " + "-"*66)
    
    for param, importance in list(results['parameter_importance'].items())[:6]:
        bar_len = int(importance * 40)
        bar = '█' * bar_len + '░' * (40 - bar_len)
        print(f"     {param:25s} [{bar}] {importance:.3f}")
    
    # Recommendations based on profile
    print(f"\n  💡 Parameter Recommendations:")
    
    best = results['best_parameters']
    
    if best.get('min_count_threshold', 10) >= 15:
        print(f"     • High min_count_threshold ({best.get('min_count_threshold')}) recommended due to zero inflation")
    
    if best.get('normalization_method') == 'TMM':
        print(f"     • TMM normalization recommended for this dataset")
    elif best.get('normalization_method') == 'RLE':
        print(f"     • RLE normalization recommended for better handling of outliers")
    
    if best.get('dispersion_type') == 'parametric':
        print(f"     • Parametric dispersion works well with moderate biological variation")
    else:
        print(f"     • Local dispersion fitting recommended for high biological variation")
    
    if best.get('fdr_threshold', 0.05) != 0.05:
        print(f"     • Consider FDR threshold of {best.get('fdr_threshold')} for this analysis")
    
    # Optimization progress
    history = results['optimization_history']
    if len(history) >= 5:
        early_scores = [h['score'] for h in history[:5]]
        late_scores = [h['score'] for h in history[-5:]]
        improvement = np.mean(late_scores) - np.mean(early_scores)
        
        if improvement > 0.05:
            print(f"\n  ✓ Significant optimization improvement ({improvement:.1%})")
        else:
            print(f"\n  ✓ Parameters converged (stable optimization)")


def run_parameter_optimization(counts_file=None, metadata_file=None, strategy='adaptive',
                               n_trials=20, cv_folds=5, demo=False):
    """Run parameter optimization."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated optimization...")
        profile = generate_demo_profile()
        
        print(f"\n📊 Demo profile:")
        print(f"   Samples: {profile['design']['n_samples']}")
        print(f"   Genes: {profile['design']['n_genes']:,}")
        print(f"   BCV: {profile['biological_variation']['bcv']:.3f}")
        print(f"   Difficulty: {profile['summary']['difficulty']}")
        
        print(f"\n🔍 Running {strategy} optimization ({n_trials} trials)...")
        
        # Show progress
        for i in range(n_trials):
            progress = (i + 1) / n_trials
            bar_len = int(progress * 40)
            bar = '█' * bar_len + '░' * (40 - bar_len)
            print(f"\r   Progress: [{bar}] {i+1}/{n_trials}", end='', flush=True)
        print()  # New line
        
        results = simple_grid_search(profile, n_iterations=n_trials)
        
    else:
        # Load data
        print(f"\n📂 Loading data from: {counts_file}")
        counts = pd.read_csv(counts_file, index_col=0)
        print(f"   Loaded: {counts.shape[0]} genes × {counts.shape[1]} samples")
        
        metadata = None
        if metadata_file:
            metadata = pd.read_csv(metadata_file)
            print(f"   Metadata: {len(metadata)} samples")
        
        # Profile data
        print("\n🔍 Profiling data characteristics...")
        profiler = RNAseqDataProfiler(counts, metadata)
        profile = profiler.run_full_profile()
        
        # Run optimization
        print(f"\n⚙️ Running {strategy} optimization ({n_trials} trials)...")
        
        optimizer = ParameterOptimizer(profile, optimization_strategy=strategy)
        results = optimizer.optimize_parameters(n_iterations=n_trials, cv_folds=cv_folds)
        
        # Generate plot
        optimizer.plot_optimization_progress('optimization_progress.png')
        print(f"\n📊 Progress plot saved: optimization_progress.png")
    
    # Display results
    display_optimization_results(results, profile)
    
    # Prepare output
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'strategy': strategy,
        'n_trials': n_trials,
        'cv_folds': cv_folds,
        'best_parameters': results['best_parameters'],
        'optimization_score': results['optimization_score'],
        'parameter_importance': results['parameter_importance'],
        'profile_summary': {
            'n_samples': profile['design']['n_samples'],
            'n_genes': profile['design']['n_genes'],
            'bcv': profile['biological_variation']['bcv'],
            'difficulty': profile['summary']['difficulty']
        }
    }
    
    return output, results


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Parameter Optimization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Optimize parameters for a dataset
  python 10_parameter_optimization.py --counts counts.csv
  
  # Use grid search strategy
  python 10_parameter_optimization.py --counts counts.csv --method grid
  
  # More optimization trials
  python 10_parameter_optimization.py --counts counts.csv --n-trials 50
  
  # Change cross-validation folds
  python 10_parameter_optimization.py --counts counts.csv --cv-folds 10
  
  # Demo mode (no data required)
  python 10_parameter_optimization.py --demo

Optimization Strategies:
  adaptive  - Data-driven adaptive optimization (default)
  grid      - Exhaustive grid search
  random    - Random search with intelligent sampling
  bayesian  - Bayesian optimization (requires additional dependencies)

Parameters Optimized:
  - min_count_threshold: Minimum count for gene filtering (1-50)
  - min_samples_expressing: Minimum samples with expression (1-n/2)
  - normalization_method: TMM, RLE, upperquartile, none
  - fdr_threshold: False discovery rate (0.001-0.2)
  - log2fc_threshold: Fold change threshold (0-2)
  - dispersion_type: parametric, local, mean
  - batch_correction_method: none, combat, limma, sva
        """
    )
    
    parser.add_argument('--counts', help='Count matrix CSV file (genes × samples)')
    parser.add_argument('--metadata', help='Sample metadata CSV file')
    parser.add_argument('--method', choices=['adaptive', 'grid', 'random', 'bayesian'],
                        default='adaptive', help='Optimization strategy')
    parser.add_argument('--n-trials', type=int, default=20, 
                        help='Number of optimization trials')
    parser.add_argument('--cv-folds', type=int, default=5,
                        help='Cross-validation folds')
    parser.add_argument('--output', '-o', help='Output JSON file')
    parser.add_argument('--output-params', help='Output best parameters JSON')
    parser.add_argument('--demo', action='store_true', help='Run in demo mode')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs
    if not args.demo and not args.counts:
        print("ERROR: Either --counts or --demo is required")
        parser.print_help()
        sys.exit(1)
    
    if args.counts and not Path(args.counts).exists():
        print(f"ERROR: Count file not found: {args.counts}")
        sys.exit(1)
    
    # Run optimization
    output, results = run_parameter_optimization(
        counts_file=args.counts,
        metadata_file=args.metadata,
        strategy=args.method,
        n_trials=args.n_trials,
        cv_folds=args.cv_folds,
        demo=args.demo
    )
    
    # Save output
    output_file = args.output or 'optimization_results.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\n💾 Results saved to: {output_file}")
    
    # Save best parameters separately if requested
    if args.output_params or args.demo:
        params_file = args.output_params or 'best_parameters.json'
        with open(params_file, 'w') as f:
            json.dump(results['best_parameters'], f, indent=2)
        print(f"💾 Best parameters saved to: {params_file}")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
