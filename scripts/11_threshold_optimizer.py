#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Script: Adaptive Threshold Optimizer (ATO)

Data-driven threshold optimization for differential expression analysis.
Replaces arbitrary cutoffs (|logFC| > 1, padj < 0.05) with statistically
justified thresholds based on your data's characteristics.

Features:
- π₀ estimation for null proportion
- Multiple logFC methods (MAD, IQR, mixture model)
- Goal-based optimization (discovery/balanced/validation)
- Publication-ready methods text generation
- Integration with ensemble analysis

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
    from raptor.threshold_optimizer import (
        AdaptiveThresholdOptimizer,
        optimize_thresholds,
        plot_optimization_summary
    )
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - Threshold Optimizer                 ║
    ║                                                              ║
    ║   🆕 Adaptive Threshold Optimizer (ATO)                      ║
    ║   Data-Driven Thresholds for DE Analysis                     ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_data(n_genes=10000, de_proportion=0.15, seed=42):
    """Generate synthetic DE data for demonstration."""
    np.random.seed(seed)
    
    n_null = int(n_genes * (1 - de_proportion))
    n_de = n_genes - n_null
    
    # Null genes: centered at 0
    null_logfc = np.random.normal(0, 0.2, n_null)
    null_pval = np.random.uniform(0.05, 1, n_null)
    
    # DE genes: shifted away from 0
    de_logfc = np.concatenate([
        np.random.normal(1.5, 0.5, n_de // 2),  # Upregulated
        np.random.normal(-1.5, 0.5, n_de - n_de // 2)  # Downregulated
    ])
    de_pval = np.random.exponential(0.001, n_de)
    de_pval = np.clip(de_pval, 1e-300, 0.05)
    
    df = pd.DataFrame({
        'gene_id': [f'Gene_{i}' for i in range(n_genes)],
        'log2FoldChange': np.concatenate([null_logfc, de_logfc]),
        'pvalue': np.concatenate([null_pval, de_pval]),
        'baseMean': np.random.exponential(1000, n_genes),
        'lfcSE': np.abs(np.random.normal(0.1, 0.05, n_genes))
    })
    df = df.set_index('gene_id')
    df['stat'] = df['log2FoldChange'] / df['lfcSE']
    
    return df


def generate_demo_result(goal='balanced'):
    """Generate demonstration ATO result."""
    
    goal_params = {
        'discovery': {
            'logfc_cutoff': 0.485,
            'padj_cutoff': 0.10,
            'padj_method': 'BH',
            'n_significant': 1842
        },
        'balanced': {
            'logfc_cutoff': 0.632,
            'padj_cutoff': 0.05,
            'padj_method': 'BH',
            'n_significant': 1523
        },
        'validation': {
            'logfc_cutoff': 0.891,
            'padj_cutoff': 0.01,
            'padj_method': 'BY',
            'n_significant': 1156
        }
    }
    
    params = goal_params.get(goal, goal_params['balanced'])
    
    result = {
        'logfc_cutoff': params['logfc_cutoff'],
        'padj_cutoff': params['padj_cutoff'],
        'padj_method': params['padj_method'],
        'goal': goal,
        'pi0_estimate': 0.847,
        'logfc_method': 'MAD',
        'n_significant_traditional': 1485,
        'n_significant_optimized': params['n_significant'],
        'methods_text': f"""Statistical thresholds were determined using the Adaptive Threshold 
Optimizer (ATO) from RAPTOR v2.1.1. The null proportion (π₀) was estimated 
at 0.847 using the Storey method. LogFC threshold was calculated using the 
MAD-based approach, yielding |log₂FC| > {params['logfc_cutoff']:.3f}. P-values were 
adjusted using the {params['padj_method']} method with α = {params['padj_cutoff']}. Analysis goal 
was set to '{goal}' mode. This data-driven approach identified 
{params['n_significant']} differentially expressed genes."""
    }
    
    return result


def display_result(result, df=None):
    """Display ATO result with formatting."""
    print("\n" + "="*70)
    print("  🎯 ADAPTIVE THRESHOLD OPTIMIZER RESULTS")
    print("="*70)
    
    # Thresholds
    print("\n  📊 Optimized Thresholds:")
    print("  " + "-"*66)
    print(f"     |log₂FC| threshold:  {result['logfc_cutoff']:.3f}")
    print(f"     Adjusted p-value:    {result['padj_cutoff']}")
    print(f"     Adjustment method:   {result['padj_method']}")
    print(f"     Analysis goal:       {result['goal']}")
    
    # Statistics
    print("\n  📈 Optimization Statistics:")
    print("  " + "-"*66)
    print(f"     π₀ estimate:         {result['pi0_estimate']:.3f} ({result['pi0_estimate']*100:.1f}% null genes)")
    print(f"     LogFC method:        {result['logfc_method']}")
    
    # Comparison
    print("\n  🔄 Traditional vs Optimized:")
    print("  " + "-"*66)
    trad = result['n_significant_traditional']
    opt = result['n_significant_optimized']
    diff = opt - trad
    pct = (diff / trad * 100) if trad > 0 else 0
    
    print(f"     Traditional (|logFC|>1, padj<0.05):  {trad:,} genes")
    print(f"     Optimized thresholds:                {opt:,} genes")
    
    if diff > 0:
        print(f"     Difference:                          +{diff:,} ({pct:+.1f}%)")
    else:
        print(f"     Difference:                          {diff:,} ({pct:+.1f}%)")
    
    # Visual comparison bar
    max_genes = max(trad, opt)
    trad_bar = int(trad / max_genes * 40) if max_genes > 0 else 0
    opt_bar = int(opt / max_genes * 40) if max_genes > 0 else 0
    
    print("\n     Traditional: [" + "█" * trad_bar + "░" * (40 - trad_bar) + f"] {trad:,}")
    print("     Optimized:   [" + "█" * opt_bar + "░" * (40 - opt_bar) + f"] {opt:,}")
    
    # Methods text
    print("\n  📝 Publication Methods Text:")
    print("  " + "-"*66)
    methods = result['methods_text']
    # Word wrap at ~60 chars
    words = methods.split()
    line = "     "
    for word in words:
        if len(line) + len(word) > 70:
            print(line)
            line = "     " + word + " "
        else:
            line += word + " "
    if line.strip():
        print(line)


def run_threshold_optimization(de_file, goal='balanced', logfc_col='log2FoldChange',
                                pvalue_col='pvalue', output_dir=None, plot=False,
                                demo=False):
    """Run adaptive threshold optimization."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated data...")
        df = generate_demo_data()
        result = generate_demo_result(goal)
        
        print(f"\n📊 Demo data generated:")
        print(f"   Genes: {len(df):,}")
        print(f"   True DE proportion: ~15%")
        
    else:
        # Load DE results
        print(f"\n📂 Loading DE results from: {de_file}")
        
        if de_file.endswith('.csv'):
            df = pd.read_csv(de_file, index_col=0)
        elif de_file.endswith('.tsv') or de_file.endswith('.txt'):
            df = pd.read_csv(de_file, sep='\t', index_col=0)
        else:
            df = pd.read_csv(de_file, index_col=0)
        
        print(f"   Loaded: {len(df):,} genes")
        
        # Check required columns
        if logfc_col not in df.columns:
            # Try common alternatives
            alt_names = ['log2FoldChange', 'logFC', 'log2FC', 'lfc', 'LogFC']
            found = False
            for alt in alt_names:
                if alt in df.columns:
                    logfc_col = alt
                    found = True
                    print(f"   Using '{alt}' as logFC column")
                    break
            if not found:
                print(f"   ERROR: logFC column '{logfc_col}' not found")
                print(f"   Available columns: {list(df.columns)}")
                sys.exit(1)
        
        if pvalue_col not in df.columns:
            alt_names = ['pvalue', 'PValue', 'pval', 'p.value', 'P.Value']
            found = False
            for alt in alt_names:
                if alt in df.columns:
                    pvalue_col = alt
                    found = True
                    print(f"   Using '{alt}' as p-value column")
                    break
            if not found:
                print(f"   ERROR: p-value column '{pvalue_col}' not found")
                print(f"   Available columns: {list(df.columns)}")
                sys.exit(1)
        
        # Run optimization
        print(f"\n🔬 Running threshold optimization (goal: {goal})...")
        
        ato = AdaptiveThresholdOptimizer(
            df=df,
            logfc_col=logfc_col,
            pvalue_col=pvalue_col,
            goal=goal,
            verbose=True
        )
        
        result_obj = ato.optimize()
        
        # Convert to dict for display
        result = {
            'logfc_cutoff': result_obj.logfc_cutoff,
            'padj_cutoff': result_obj.padj_cutoff,
            'padj_method': result_obj.padj_method,
            'goal': goal,
            'pi0_estimate': result_obj.pi0_estimate,
            'logfc_method': result_obj.logfc_method,
            'n_significant_traditional': result_obj.n_significant_traditional,
            'n_significant_optimized': result_obj.n_significant_optimized,
            'methods_text': result_obj.methods_text
        }
        
        # Generate plots if requested
        if plot and output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            print(f"\n📊 Generating visualization...")
            plot_file = output_path / 'ato_optimization_summary.png'
            plot_optimization_summary(
                result=result_obj,
                df=ato.df,
                save_path=plot_file
            )
            print(f"   Saved: {plot_file}")
        
        # Get significant genes
        sig_genes = ato.get_significant_genes()
        
        # Save significant genes if output dir specified
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            
            sig_file = output_path / 'significant_genes_ato.csv'
            sig_genes.to_csv(sig_file)
            print(f"   Saved significant genes: {sig_file}")
    
    # Display results
    display_result(result, df)
    
    # Prepare output
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'input_file': str(de_file) if not demo else 'demo_data',
        'goal': goal,
        'optimization_results': {
            'logfc_cutoff': result['logfc_cutoff'],
            'padj_cutoff': result['padj_cutoff'],
            'padj_method': result['padj_method'],
            'pi0_estimate': result['pi0_estimate'],
            'logfc_method': result['logfc_method']
        },
        'gene_counts': {
            'traditional': result['n_significant_traditional'],
            'optimized': result['n_significant_optimized']
        },
        'methods_text': result['methods_text']
    }
    
    return output


def compare_goals(de_file=None, demo=False):
    """Compare different analysis goals."""
    print("\n" + "="*70)
    print("  📊 COMPARING ANALYSIS GOALS")
    print("="*70)
    
    goals = ['discovery', 'balanced', 'validation']
    results = {}
    
    for goal in goals:
        if demo or not RAPTOR_AVAILABLE:
            result = generate_demo_result(goal)
        else:
            df = pd.read_csv(de_file, index_col=0)
            ato = AdaptiveThresholdOptimizer(df=df, goal=goal, verbose=False)
            result_obj = ato.optimize()
            result = {
                'logfc_cutoff': result_obj.logfc_cutoff,
                'padj_cutoff': result_obj.padj_cutoff,
                'padj_method': result_obj.padj_method,
                'n_significant': result_obj.n_significant_optimized
            }
        
        results[goal] = result
    
    # Display comparison table
    print("\n  Goal         |log₂FC|   padj     Method    DE Genes")
    print("  " + "-"*60)
    
    for goal in goals:
        r = results[goal]
        n_sig = r.get('n_significant_optimized', r.get('n_significant', 0))
        print(f"  {goal:12s}  {r['logfc_cutoff']:>7.3f}   {r['padj_cutoff']:<8}  {r['padj_method']:8s}  {n_sig:,}")
    
    print("\n  💡 Interpretation:")
    print("     • Discovery:   Liberal thresholds → More candidates for validation")
    print("     • Balanced:    Trade-off between sensitivity and specificity")
    print("     • Validation:  Conservative thresholds → High-confidence hits")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Adaptive Threshold Optimizer (ATO)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with DESeq2 results
  python 11_threshold_optimizer.py --input deseq2_results.csv

  # Specify analysis goal
  python 11_threshold_optimizer.py --input results.csv --goal discovery

  # Compare all goals
  python 11_threshold_optimizer.py --input results.csv --compare-goals

  # Generate visualization
  python 11_threshold_optimizer.py --input results.csv --plot --output ato_results/

  # Demo mode (no data required)
  python 11_threshold_optimizer.py --demo

  # With custom column names (edgeR output)
  python 11_threshold_optimizer.py --input edger.csv --logfc-col logFC --pvalue-col PValue

Analysis Goals:
  discovery   - Liberal thresholds, maximize candidate detection
  balanced    - Balance sensitivity and specificity (default)
  validation  - Conservative thresholds, high-confidence hits

Supported Input Formats:
  • DESeq2 output (log2FoldChange, pvalue)
  • edgeR output (logFC, PValue)
  • limma output (logFC, P.Value)
  • NOISeq output (log2FC, prob)
  • Any CSV/TSV with logFC and p-value columns
        """
    )
    
    parser.add_argument('--input', '-i', help='DE results file (CSV/TSV)')
    parser.add_argument('--goal', '-g', choices=['discovery', 'balanced', 'validation'],
                        default='balanced', help='Analysis goal (default: balanced)')
    parser.add_argument('--logfc-col', default='log2FoldChange',
                        help='Column name for log fold change')
    parser.add_argument('--pvalue-col', default='pvalue',
                        help='Column name for p-value')
    parser.add_argument('--output', '-o', help='Output directory')
    parser.add_argument('--plot', action='store_true', help='Generate visualization')
    parser.add_argument('--compare-goals', action='store_true',
                        help='Compare all analysis goals')
    parser.add_argument('--demo', action='store_true', help='Run in demo mode')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs
    if not args.demo and not args.input:
        print("ERROR: Either --input or --demo is required")
        parser.print_help()
        sys.exit(1)
    
    if args.input and not args.demo and not Path(args.input).exists():
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)
    
    # Compare goals mode
    if args.compare_goals:
        results = compare_goals(args.input, demo=args.demo)
        
        if args.output:
            output_path = Path(args.output)
            output_path.mkdir(parents=True, exist_ok=True)
            
            with open(output_path / 'goal_comparison.json', 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            print(f"\n💾 Comparison saved to: {output_path / 'goal_comparison.json'}")
    
    else:
        # Run single optimization
        results = run_threshold_optimization(
            de_file=args.input,
            goal=args.goal,
            logfc_col=args.logfc_col,
            pvalue_col=args.pvalue_col,
            output_dir=args.output,
            plot=args.plot,
            demo=args.demo
        )
        
        # Save output
        output_file = Path(args.output) / 'ato_results.json' if args.output else 'ato_results.json'
        if args.output:
            Path(args.output).mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"\n💾 Results saved to: {output_file}")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
