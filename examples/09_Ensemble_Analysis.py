#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Ensemble Analysis (Module 9)

Combines differential expression results from multiple methods to create
high-confidence consensus gene lists using statistically rigorous ensemble
techniques.

This is Module 9 of the RAPTOR workflow (Stage 3: DE Analysis):
  M7: Import DE Results → DEResult objects
  M8: Parameter Optimization → Optimized parameters
  M9: Ensemble Analysis (THIS SCRIPT) → Consensus genes
  M10: Biomarker Discovery

Ensemble Methods:
  1. Voting - Simple count-based (high confidence)
  2. Weighted - Performance-weighted (when validation available)
  3. Fisher's - P-value combination (maximum sensitivity)
  4. Brown's - Correlation-aware combination
  5. RRA - Robust Rank Aggregation (handles outliers)

Input: Multiple DE results from Module 7 (de_result.pkl files)
Output: results/ensemble_analysis/
    - consensus_genes.csv
    - ensemble_statistics.json
    - method_comparison.csv

CRITICAL: This module includes an adapter to bridge DEResult (Module 7)
with ensemble.py expectations. See adapt_deresult_for_ensemble() function.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import argparse
import sys
import json
from pathlib import Path
from datetime import datetime
from types import SimpleNamespace

# =============================================================================
# CONSTANTS - Architecture Compliant (v2.2.0)
# =============================================================================
DEFAULT_OUTPUT_DIR = "results/ensemble_analysis"
DEFAULT_MIN_METHODS = 2
DEFAULT_SIGNIFICANCE = 0.05

# Check for dependencies
try:
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: numpy and pandas are required")
    print("Install with: pip install numpy pandas scipy")
    sys.exit(1)

# RAPTOR imports
RAPTOR_AVAILABLE = True
try:
    from raptor.de_import import DEResult
    from raptor.ensemble import (
        ensemble_voting,
        ensemble_weighted,
        ensemble_pvalue_combination,
        ensemble_rra,
        ensemble_fisher,
        ensemble_brown,
        EnsembleResult
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
            self.n_significant = len(results_df[results_df['adjusted_p_value'] < 0.05])
        
        @classmethod
        def load(cls, path):
            return None
    
    class EnsembleResult:
        def __init__(self, consensus_genes, n_consensus, method, n_methods, method_names):
            self.consensus_genes = consensus_genes
            self.n_consensus_genes = n_consensus
            self.ensemble_method = method
            self.n_methods = n_methods
            self.method_names = method_names
            self.method_statistics = {}


# =============================================================================
# ADAPTER FUNCTION - CRITICAL FOR DERESULT COMPATIBILITY
# =============================================================================

def adapt_deresult_for_ensemble(de_result):
    """
    Adapt DEResult from Module 7 to format expected by ensemble.py.
    
    CRITICAL: This adapter is necessary because:
    - DEResult uses .results_df, ensemble.py expects .data
    - DEResult uses standardized names, ensemble.py expects R names
    - DEResult has gene_id as index, ensemble.py expects it as column
    
    Parameters
    ----------
    de_result : DEResult
        DEResult object from Module 7
    
    Returns
    -------
    SimpleNamespace
        Adapted object with .data attribute and R column names
    
    Conversions
    -----------
    Attribute: .results_df → .data
    Column: adjusted_p_value → padj
    Column: log2_fold_change → log2FoldChange
    Column: p_value → pvalue
    Column: base_mean → baseMean (if present)
    Index: gene_id → gene_id column
    """
    # Extract DataFrame
    df = de_result.results_df.copy()
    
    # Convert gene_id from index to column
    df = df.reset_index()
    
    # Rename columns to what ensemble.py expects (R naming convention)
    column_mapping = {
        'adjusted_p_value': 'padj',
        'log2_fold_change': 'log2FoldChange',
        'p_value': 'pvalue',
        'base_mean': 'baseMean'
    }
    
    df = df.rename(columns=column_mapping)
    
    # Create adapter object with .data attribute
    adapter = SimpleNamespace(
        data=df,
        pipeline=de_result.pipeline,
        parameters=de_result.parameters,
        metadata=de_result.metadata
    )
    
    return adapter


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║    🦖 RAPTOR v2.2.0 - Ensemble Analysis (Module 9)           ║
    ║                                                              ║
    ║    Combine Multiple DE Methods → High-Confidence Genes      ║
    ║    ✅ 5 Methods | Statistical Rigor | Publication Quality    ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_de_result(method_name, n_genes=15000, seed=None):
    """Generate demo DE results for testing."""
    if seed is not None:
        np.random.seed(seed)
    
    # Generate gene IDs
    gene_ids = [f'ENSG{i+1:011d}' for i in range(n_genes)]
    
    # Generate realistic DE results
    # Method-specific biases to make ensemble interesting
    bias_factor = {'DESeq2': 1.0, 'edgeR': 1.2, 'limma': 0.9}.get(method_name, 1.0)
    
    base_mean = np.random.gamma(shape=2, scale=100, size=n_genes)
    log2_fold_change = np.random.normal(0, 1.5 * bias_factor, n_genes)
    
    # P-values with some truly DE genes
    n_de = int(n_genes * 0.08)  # 8% truly DE
    p_values = np.ones(n_genes)
    p_values[:n_de] = np.random.beta(0.5, 10, n_de)
    p_values[n_de:] = np.random.uniform(0.1, 1.0, n_genes - n_de)
    
    # Shuffle
    indices = np.random.permutation(n_genes)
    p_values = p_values[indices]
    
    # Adjusted p-values
    adjusted_p_values = np.minimum(p_values * 10, 1.0)
    
    # Create DataFrame with standardized names (DEResult format)
    df = pd.DataFrame({
        'log2_fold_change': log2_fold_change,
        'p_value': p_values,
        'adjusted_p_value': adjusted_p_values,
        'base_mean': base_mean
    }, index=gene_ids)
    df.index.name = 'gene_id'
    
    # Create DEResult
    de_result = DEResult(
        results_df=df,
        pipeline=method_name.upper(),
        parameters={'fdr_threshold': 0.05, 'lfc_threshold': 0.0},
        metadata={'source': 'demo', 'timestamp': datetime.now().isoformat()}
    )
    
    return de_result


def run_voting_ensemble(de_results_adapted, output_dir, min_methods=2, demo=False):
    """
    Run voting ensemble analysis.
    
    Best for: High-confidence gene lists
    Method: Count how many methods detect each gene
    """
    print("\n" + "="*70)
    print("  METHOD 1: Voting Ensemble")
    print("="*70)
    print("\n  Simple count-based voting - HIGH CONFIDENCE")
    print(f"  Require detection by ≥{min_methods} methods")
    print()
    
    output_path = Path(output_dir) / 'voting'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • Methods: {list(de_results_adapted.keys())}")
        print(f"     • Min methods: {min_methods}")
        print()
        
        # Simulate voting
        consensus_genes = pd.DataFrame({
            'gene_id': [f'ENSG{i:011d}' for i in range(1, 251)],
            'n_votes': np.random.randint(min_methods, len(de_results_adapted)+1, 250),
            'direction': np.random.choice(['up', 'down'], 250),
            'direction_agreement': np.random.uniform(0.8, 1.0, 250)
        })
        
        print(f"  ✓ Voting found {len(consensus_genes)} consensus genes")
        print(f"     • All {len(de_results_adapted)} methods: {(consensus_genes['n_votes'] == len(de_results_adapted)).sum()}")
        print(f"     • {len(de_results_adapted)-1} methods: {(consensus_genes['n_votes'] == len(de_results_adapted)-1).sum()}")
        print()
        
        # Save
        consensus_genes.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return consensus_genes
    
    # Real ensemble
    print("  🚀 Running voting ensemble with RAPTOR...")
    
    try:
        result = ensemble_voting(
            de_results=de_results_adapted,
            min_methods=min_methods,
            filters={'padj': 0.05},  # Pre-filter for significance
            check_direction=True,
            direction_threshold=1.0
        )
        
        print(f"  ✓ Voting found {len(result)} consensus genes")
        
        # Save
        result.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during voting: {e}")
        return None


def run_fisher_ensemble(de_results_adapted, output_dir, demo=False):
    """
    Run Fisher's method for p-value combination.
    
    Best for: Maximum sensitivity, exploratory analysis
    Method: Combines p-values using Fisher's method
    """
    print("\n" + "="*70)
    print("  METHOD 2: Fisher's Method")
    print("="*70)
    print("\n  P-value combination - MAXIMUM SENSITIVITY")
    print("  Based on: Fisher (1925)")
    print()
    
    output_path = Path(output_dir) / 'fisher'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • Methods: {list(de_results_adapted.keys())}")
        print()
        
        # Simulate Fisher's method
        consensus_genes = pd.DataFrame({
            'gene_id': [f'ENSG{i:011d}' for i in range(1, 401)],
            'combined_pvalue': np.random.beta(0.5, 10, 400),
            'combined_padj': np.random.beta(1, 5, 400),
            'direction': np.random.choice(['up', 'down'], 400),
            'meta_lfc': np.random.normal(0, 2, 400)
        })
        
        n_sig = (consensus_genes['combined_padj'] < 0.05).sum()
        
        print(f"  ✓ Fisher's method found {n_sig} consensus genes")
        print(f"     • Combined p-value < 0.05: {(consensus_genes['combined_pvalue'] < 0.05).sum()}")
        print()
        
        # Save
        consensus_genes.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return EnsembleResult(
            consensus_genes=consensus_genes[consensus_genes['combined_padj'] < 0.05],
            n_consensus=n_sig,
            method='fisher',
            n_methods=len(de_results_adapted),
            method_names=list(de_results_adapted.keys())
        )
    
    # Real ensemble
    print("  🚀 Running Fisher's method with RAPTOR...")
    print("     ⚠️  Using raw p-values (not adjusted)")
    
    try:
        result = ensemble_fisher(
            de_results=de_results_adapted,
            use_padj=False,  # CRITICAL: Use raw p-values
            significance_threshold=0.05,
            check_direction=True,
            direction_threshold=1.0,
            output_dir=str(output_path)
        )
        
        print(f"  ✓ Fisher's method found {result.n_consensus_genes} consensus genes")
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during Fisher's method: {e}")
        return None


def run_brown_ensemble(de_results_adapted, output_dir, demo=False):
    """
    Run Brown's method for correlation-aware p-value combination.
    
    Best for: When methods use same data/normalization
    Method: Accounts for correlation between methods
    """
    print("\n" + "="*70)
    print("  METHOD 3: Brown's Method")
    print("="*70)
    print("\n  Correlation-aware combination")
    print("  Based on: Brown (1975) - Accounts for method correlation")
    print()
    
    output_path = Path(output_dir) / 'brown'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • Methods: {list(de_results_adapted.keys())}")
        print()
        
        # Simulate Brown's method
        consensus_genes = pd.DataFrame({
            'gene_id': [f'ENSG{i:011d}' for i in range(1, 351)],
            'combined_pvalue': np.random.beta(0.5, 10, 350),
            'combined_padj': np.random.beta(1, 5, 350),
            'direction': np.random.choice(['up', 'down'], 350),
            'meta_lfc': np.random.normal(0, 2, 350)
        })
        
        n_sig = (consensus_genes['combined_padj'] < 0.05).sum()
        
        print(f"  ✓ Brown's method found {n_sig} consensus genes")
        print(f"     • More conservative than Fisher's (accounts for correlation)")
        print()
        
        # Save
        consensus_genes.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return EnsembleResult(
            consensus_genes=consensus_genes[consensus_genes['combined_padj'] < 0.05],
            n_consensus=n_sig,
            method='brown',
            n_methods=len(de_results_adapted),
            method_names=list(de_results_adapted.keys())
        )
    
    # Real ensemble
    print("  🚀 Running Brown's method with RAPTOR...")
    
    try:
        result = ensemble_brown(
            de_results=de_results_adapted,
            use_padj=False,
            significance_threshold=0.05,
            check_direction=True,
            direction_threshold=1.0,
            output_dir=str(output_path)
        )
        
        print(f"  ✓ Brown's method found {result.n_consensus_genes} consensus genes")
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during Brown's method: {e}")
        return None


def run_rra_ensemble(de_results_adapted, output_dir, demo=False):
    """
    Run Robust Rank Aggregation.
    
    Best for: Robust ranking, handling outliers
    Method: Combines ranked gene lists using order statistics
    """
    print("\n" + "="*70)
    print("  METHOD 4: Robust Rank Aggregation (RRA)")
    print("="*70)
    print("\n  Rank-based combination - HANDLES OUTLIERS")
    print("  Based on: Kolde et al. (2012)")
    print()
    
    output_path = Path(output_dir) / 'rra'
    output_path.mkdir(parents=True, exist_ok=True)
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • Methods: {list(de_results_adapted.keys())}")
        print()
        
        # Simulate RRA
        consensus_genes = pd.DataFrame({
            'gene_id': [f'ENSG{i:011d}' for i in range(1, 301)],
            'rra_score': np.random.beta(0.5, 10, 300),
            'rra_padj': np.random.beta(1, 5, 300),
            'direction': np.random.choice(['up', 'down'], 300),
            'mean_rank': np.random.uniform(1, 1000, 300)
        })
        
        n_sig = (consensus_genes['rra_padj'] < 0.05).sum()
        
        print(f"  ✓ RRA found {n_sig} consensus genes")
        print(f"     • Robust to outliers and method disagreements")
        print()
        
        # Save
        consensus_genes.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return EnsembleResult(
            consensus_genes=consensus_genes[consensus_genes['rra_padj'] < 0.05],
            n_consensus=n_sig,
            method='rra',
            n_methods=len(de_results_adapted),
            method_names=list(de_results_adapted.keys())
        )
    
    # Real ensemble
    print("  🚀 Running RRA with RAPTOR...")
    
    try:
        result = ensemble_rra(
            de_results=de_results_adapted,
            rank_by='pvalue',
            use_padj=False,
            significance_threshold=0.05,
            check_direction=True,
            output_dir=str(output_path)
        )
        
        print(f"  ✓ RRA found {result.n_consensus_genes} consensus genes")
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during RRA: {e}")
        return None


def run_weighted_ensemble(de_results_adapted, output_dir, weights=None, demo=False):
    """
    Run weighted ensemble (when performance data available).
    
    Best for: When you have validation/performance data
    Method: Weight methods by their performance
    """
    print("\n" + "="*70)
    print("  METHOD 5: Weighted Ensemble")
    print("="*70)
    print("\n  Performance-weighted combination")
    print("  Best when you have validation data")
    print()
    
    output_path = Path(output_dir) / 'weighted'
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Default weights if not provided
    if weights is None:
        weights = {name: 1.0 for name in de_results_adapted.keys()}
    
    if demo or not RAPTOR_AVAILABLE:
        print("  🎮 Running in DEMO mode...")
        print(f"     • Methods: {list(de_results_adapted.keys())}")
        print(f"     • Weights: {weights}")
        print()
        
        # Simulate weighted ensemble
        consensus_genes = pd.DataFrame({
            'gene_id': [f'ENSG{i:011d}' for i in range(1, 281)],
            'weighted_score': np.random.uniform(0.5, 1.0, 280),
            'direction': np.random.choice(['up', 'down'], 280),
            'total_weight': np.random.uniform(2.0, 3.0, 280)
        })
        
        print(f"  ✓ Weighted ensemble found {len(consensus_genes)} consensus genes")
        print()
        
        # Save
        consensus_genes.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return consensus_genes
    
    # Real ensemble
    print("  🚀 Running weighted ensemble with RAPTOR...")
    
    try:
        result = ensemble_weighted(
            de_results=de_results_adapted,
            weights=weights,
            min_weight=0.5,
            check_direction=True,
            direction_threshold=1.0
        )
        
        print(f"  ✓ Weighted ensemble found {len(result)} consensus genes")
        
        # Save
        result.to_csv(output_path / 'consensus_genes.csv', index=False)
        
        return result
        
    except Exception as e:
        print(f"\n  ❌ Error during weighted ensemble: {e}")
        return None


def compare_ensemble_methods(results_dict, output_dir):
    """Compare results from different ensemble methods."""
    print("\n" + "="*70)
    print("  COMPARISON: All Ensemble Methods")
    print("="*70)
    print()
    
    comparison_data = []
    
    for method_name, result in results_dict.items():
        if result is None:
            continue
        
        # Handle different result types
        if hasattr(result, 'n_consensus_genes'):
            n_genes = result.n_consensus_genes
        elif isinstance(result, pd.DataFrame):
            n_genes = len(result)
        else:
            n_genes = 0
        
        comparison_data.append({
            'Method': method_name,
            'Consensus Genes': n_genes,
            'Stringency': {
                'Voting': 'High',
                'Fisher': 'Low',
                'Brown': 'Medium',
                'RRA': 'Medium',
                'Weighted': 'Variable'
            }.get(method_name, 'Unknown')
        })
    
    if not comparison_data:
        print("  No results to compare")
        return
    
    # Create comparison table
    comparison_df = pd.DataFrame(comparison_data)
    
    print("  Ensemble Method Comparison:")
    print("  " + "-"*50)
    print(f"  {'Method':<15} {'Genes':<15} {'Stringency':<15}")
    print("  " + "-"*50)
    
    for _, row in comparison_df.iterrows():
        print(f"  {row['Method']:<15} {row['Consensus Genes']:<15} {row['Stringency']:<15}")
    
    print("  " + "-"*50)
    print()
    
    # Recommendations
    print("  📊 Method Selection Guide:")
    print()
    print("     • VOTING: High-confidence genes (strict)")
    print("     • FISHER: Maximum sensitivity (exploratory)")
    print("     • BROWN: Correlation-aware (conservative)")
    print("     • RRA: Robust to outliers (balanced)")
    print("     • WEIGHTED: Use when you have performance data")
    print()
    
    # Save comparison
    output_path = Path(output_dir)
    comparison_df.to_csv(output_path / 'method_comparison.csv', index=False)
    
    with open(output_path / 'ensemble_summary.json', 'w') as f:
        json.dump({
            'methods_compared': len(comparison_data),
            'results': comparison_data,
            'timestamp': datetime.now().isoformat()
        }, f, indent=2)
    
    print(f"  ✓ Saved: {output_path / 'method_comparison.csv'}")
    print(f"  ✓ Saved: {output_path / 'ensemble_summary.json'}")


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Ensemble Analysis (Module 9)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (no data needed)
  python 09_Ensemble_Analysis.py --demo
  
  # Voting ensemble (high confidence)
  python 09_Ensemble_Analysis.py \\
      --de-results results/de_imported/deseq2_result.pkl \\
                   results/de_imported/edger_result.pkl \\
                   results/de_imported/limma_result.pkl \\
      --method voting \\
      --min-methods 3
  
  # Fisher's method (maximum sensitivity)
  python 09_Ensemble_Analysis.py \\
      --de-results results/de_imported/*.pkl \\
      --method fisher
  
  # All methods
  python 09_Ensemble_Analysis.py \\
      --de-results results/de_imported/*.pkl \\
      --method all

Ensemble Methods:
  
  1. Voting (High Confidence)
     • Simple count-based
     • Requires detection by ≥N methods
     • Most conservative
  
  2. Fisher's Method (Maximum Sensitivity)
     • P-value combination
     • Detects weak signals
     • Good for exploratory analysis
  
  3. Brown's Method (Correlation-Aware)
     • Accounts for method correlation
     • More conservative than Fisher's
     • Best when methods use same data
  
  4. RRA (Robust Ranking)
     • Rank aggregation
     • Handles outliers well
     • Balanced approach
  
  5. Weighted (Performance-Based)
     • Weight by method performance
     • Requires validation data
     • Highly customizable

Output Files:
  results/ensemble_analysis/
    ├── voting/consensus_genes.csv
    ├── fisher/consensus_genes.csv
    ├── brown/consensus_genes.csv
    ├── rra/consensus_genes.csv
    ├── weighted/consensus_genes.csv
    ├── method_comparison.csv
    └── ensemble_summary.json
        """
    )
    
    parser.add_argument('--de-results', '-d', nargs='+',
                       help='DE result files from Module 7 (.pkl files)')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--method', '-m',
                       choices=['voting', 'fisher', 'brown', 'rra', 'weighted', 'all'],
                       default='all',
                       help='Ensemble method (default: all)')
    parser.add_argument('--min-methods', type=int, default=DEFAULT_MIN_METHODS,
                       help=f'Minimum methods for voting (default: {DEFAULT_MIN_METHODS})')
    parser.add_argument('--significance', type=float, default=DEFAULT_SIGNIFICANCE,
                       help=f'Significance threshold (default: {DEFAULT_SIGNIFICANCE})')
    parser.add_argument('--weights',
                       help='Method weights (JSON format: {"DESeq2": 1.0, "edgeR": 0.8})')
    parser.add_argument('--demo', action='store_true',
                       help='Run in demo mode with simulated data')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs for real run
    if not args.demo and not args.de_results:
        print("ERROR: --de-results is required (or use --demo)")
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
        
        # Generate demo DE results for 3 methods
        print("  Generating DE results for 3 methods...")
        deseq2_result = generate_demo_de_result('DESeq2', seed=42)
        edger_result = generate_demo_de_result('edgeR', seed=123)
        limma_result = generate_demo_de_result('limma', seed=456)
        
        de_results = {
            'DESeq2': deseq2_result,
            'edgeR': edger_result,
            'limma': limma_result
        }
        
        print(f"  ✓ Generated: 3 methods, {len(deseq2_result.results_df):,} genes each")
        
        # Show overlap
        for name, result in de_results.items():
            print(f"     • {name}: {result.n_significant} significant genes")
        
    else:
        print("\n🚀 LOADING REAL DATA")
        print("─" * 60)
        
        # Load DE results
        de_results = {}
        
        for de_file in args.de_results:
            print(f"  Loading: {de_file}")
            try:
                de_result = DEResult.load(de_file)
                
                # Extract method name from filename or use pipeline
                method_name = de_result.pipeline
                de_results[method_name] = de_result
                
                print(f"  ✓ Loaded: {method_name} ({de_result.n_genes:,} genes, "
                      f"{de_result.n_significant} significant)")
                
            except Exception as e:
                print(f"  ❌ Error loading {de_file}: {e}")
        
        if len(de_results) < 2:
            print("\nERROR: Need at least 2 DE results for ensemble analysis")
            sys.exit(1)
    
    # =========================================================================
    # Adapt DEResults for Ensemble Module
    # =========================================================================
    
    print("\n📦 ADAPTING DE RESULTS FOR ENSEMBLE MODULE")
    print("─" * 60)
    print("  Converting DEResult format → ensemble.py format...")
    print("  (DEResult uses .results_df, ensemble.py expects .data)")
    print()
    
    de_results_adapted = {}
    
    for method_name, de_result in de_results.items():
        adapted = adapt_deresult_for_ensemble(de_result)
        de_results_adapted[method_name] = adapted
        print(f"  ✓ Adapted: {method_name}")
    
    print()
    print(f"  ✓ Ready: {len(de_results_adapted)} methods adapted for ensemble")
    
    # =========================================================================
    # Run Ensemble Methods
    # =========================================================================
    
    results = {}
    
    if args.method == 'voting' or args.method == 'all':
        result = run_voting_ensemble(
            de_results_adapted, args.output,
            min_methods=args.min_methods,
            demo=args.demo
        )
        if result is not None:
            results['Voting'] = result
    
    if args.method == 'fisher' or args.method == 'all':
        result = run_fisher_ensemble(
            de_results_adapted, args.output,
            demo=args.demo
        )
        if result is not None:
            results['Fisher'] = result
    
    if args.method == 'brown' or args.method == 'all':
        result = run_brown_ensemble(
            de_results_adapted, args.output,
            demo=args.demo
        )
        if result is not None:
            results['Brown'] = result
    
    if args.method == 'rra' or args.method == 'all':
        result = run_rra_ensemble(
            de_results_adapted, args.output,
            demo=args.demo
        )
        if result is not None:
            results['RRA'] = result
    
    if args.method == 'weighted' or args.method == 'all':
        # Parse weights if provided
        weights = None
        if args.weights:
            try:
                weights = json.loads(args.weights)
            except:
                print("WARNING: Could not parse weights, using defaults")
        
        result = run_weighted_ensemble(
            de_results_adapted, args.output,
            weights=weights,
            demo=args.demo
        )
        if result is not None:
            results['Weighted'] = result
    
    # =========================================================================
    # Compare Methods
    # =========================================================================
    
    if len(results) > 1:
        compare_ensemble_methods(results, args.output)
    
    # =========================================================================
    # Final Summary
    # =========================================================================
    
    print("\n" + "="*70)
    print("  ✅ MODULE 9 (ENSEMBLE ANALYSIS) COMPLETE!")
    print("="*70)
    
    print(f"\n  📂 Output Directory: {output_path}")
    
    if len(results) > 0:
        print(f"\n  📊 Ensemble Results:")
        for method_name, result in results.items():
            if hasattr(result, 'n_consensus_genes'):
                n_genes = result.n_consensus_genes
            elif isinstance(result, pd.DataFrame):
                n_genes = len(result)
            else:
                n_genes = 0
            print(f"     • {method_name}: {n_genes} consensus genes")
    
    print(f"\n  📄 Output Files:")
    for method_name in results.keys():
        method_dir = method_name.lower()
        print(f"     • {method_dir}/consensus_genes.csv")
    
    if len(results) > 1:
        print(f"     • method_comparison.csv")
        print(f"     • ensemble_summary.json")
    
    print(f"\n  🔜 Next Steps:")
    print(f"\n     Module 10 - Biomarker Discovery:")
    print(f"     Use consensus genes for biomarker identification")
    print(f"")
    print(f"     Or apply to your R analysis:")
    print(f"     # Use high-confidence genes from voting/RRA")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
