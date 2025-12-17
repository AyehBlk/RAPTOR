#!/usr/bin/env python3
"""
Example: Using the Adaptive Threshold Optimizer (ATO)

This script demonstrates how to use the ATO module for
data-driven threshold optimization in DE analysis.

Author: Ayeh Bolouki
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Import ATO
from raptor.threshold_optimizer import (
    AdaptiveThresholdOptimizer,
    optimize_thresholds,
    plot_optimization_summary,
    plot_volcano,
    plot_threshold_comparison
)


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
        'log2FoldChange': np.concatenate([null_logfc, de_logfc]),
        'pvalue': np.concatenate([null_pval, de_pval]),
        'baseMean': np.random.exponential(1000, n_genes),
        'lfcSE': np.abs(np.random.normal(0.1, 0.05, n_genes))
    })
    df.index = [f'Gene_{i}' for i in range(n_genes)]
    
    # Calculate stat column
    df['stat'] = df['log2FoldChange'] / df['lfcSE']
    
    return df


def example_basic_usage():
    """Example 1: Basic usage with convenience function."""
    print("\n" + "="*60)
    print("EXAMPLE 1: Basic Usage")
    print("="*60)
    
    # Load or generate data
    df = generate_demo_data()
    print(f"Generated {len(df)} synthetic genes")
    
    # Quick optimization
    result = optimize_thresholds(df, goal='discovery')
    
    # Print summary
    print(result.summary())
    
    return result, df


def example_full_control():
    """Example 2: Full control with AdaptiveThresholdOptimizer class."""
    print("\n" + "="*60)
    print("EXAMPLE 2: Full Control")
    print("="*60)
    
    # Load data
    df = generate_demo_data(n_genes=15000, de_proportion=0.1)
    
    # Create optimizer with custom settings
    ato = AdaptiveThresholdOptimizer(
        df=df,
        goal='balanced',
        verbose=True
    )
    
    # Run optimization with specific method
    result = ato.optimize(
        logfc_method='mad',  # Use only MAD method
        n1=4,  # 4 samples per group
        n2=4
    )
    
    print(result.summary())
    
    # Get significant genes
    sig_genes = ato.get_significant_genes()
    print(f"\nTop 10 significant genes:")
    print(sig_genes.sort_values('pvalue').head(10)[['logfc', 'pvalue', 'padj_optimized']])
    
    return ato, result


def example_comparison():
    """Example 3: Compare different thresholds."""
    print("\n" + "="*60)
    print("EXAMPLE 3: Threshold Comparison")
    print("="*60)
    
    df = generate_demo_data()
    ato = AdaptiveThresholdOptimizer(df, verbose=False)
    ato.optimize()
    
    # Compare across thresholds
    comparison = ato.compare_thresholds(
        logfc_values=[0.25, 0.5, 0.75, 1.0, 1.5, 2.0],
        padj_values=[0.001, 0.01, 0.05, 0.1]
    )
    
    print("\nDE genes by threshold combination:")
    print(comparison.to_string(index=False))
    
    return comparison


def example_adjustment_methods():
    """Example 4: Compare p-value adjustment methods."""
    print("\n" + "="*60)
    print("EXAMPLE 4: P-value Adjustment Methods")
    print("="*60)
    
    df = generate_demo_data()
    ato = AdaptiveThresholdOptimizer(df, verbose=False)
    
    # Get comparison of all methods
    adj_comparison = ato.get_adjustment_comparison()
    
    # Summary
    print("\nSignificant genes (padj < 0.05) by method:")
    for method in ['BH', 'BY', 'Holm', 'Hochberg', 'Bonferroni', 'qvalue']:
        count = (adj_comparison[method] < 0.05).sum()
        print(f"  {method}: {count}")
    
    return adj_comparison


def example_visualization():
    """Example 5: Generate visualizations."""
    print("\n" + "="*60)
    print("EXAMPLE 5: Visualization")
    print("="*60)
    
    df = generate_demo_data()
    ato = AdaptiveThresholdOptimizer(df, verbose=False)
    result = ato.optimize()
    
    # Create output directory
    output_dir = Path('ato_output')
    output_dir.mkdir(exist_ok=True)
    
    # 1. Comprehensive summary plot
    fig1 = plot_optimization_summary(
        result=result,
        df=ato.df,
        save_path=output_dir / 'optimization_summary.png'
    )
    print(f"Saved: {output_dir / 'optimization_summary.png'}")
    plt.close(fig1)
    
    # 2. Volcano plot
    fig2 = plot_volcano(
        df=ato.df,
        logfc_col='logfc',
        padj_col='padj_optimized',
        logfc_cutoff=result.logfc_cutoff,
        padj_cutoff=result.padj_cutoff,
        save_path=output_dir / 'volcano_plot.png'
    )
    print(f"Saved: {output_dir / 'volcano_plot.png'}")
    plt.close(fig2)
    
    # 3. Threshold comparison heatmap
    fig3 = plot_threshold_comparison(
        df=ato.df,
        logfc_col='logfc',
        padj_col='padj_optimized',
        optimized_logfc=result.logfc_cutoff,
        optimized_padj=result.padj_cutoff,
        save_path=output_dir / 'threshold_heatmap.png'
    )
    print(f"Saved: {output_dir / 'threshold_heatmap.png'}")
    plt.close(fig3)
    
    print(f"\nAll plots saved to {output_dir}/")
    
    return output_dir


def example_real_data():
    """Example 6: Working with real DESeq2 output."""
    print("\n" + "="*60)
    print("EXAMPLE 6: Real DESeq2 Data")
    print("="*60)
    
    # Example of loading real DESeq2 output
    example_path = Path('deseq2_results.csv')
    
    if example_path.exists():
        df = pd.read_csv(example_path, index_col=0)
        print(f"Loaded {len(df)} genes from {example_path}")
    else:
        print("No real data found, using synthetic data")
        print("To use real data, save your DESeq2 output as 'deseq2_results.csv'")
        df = generate_demo_data()
    
    # Optimize
    result = optimize_thresholds(df, goal='discovery')
    print(result.summary())
    
    # Export results
    ato = AdaptiveThresholdOptimizer(df, verbose=False)
    ato.optimize()
    
    sig_genes = ato.get_significant_genes()
    sig_genes.to_csv('significant_genes_optimized.csv')
    print(f"\nExported {len(sig_genes)} significant genes to significant_genes_optimized.csv")
    
    return result


def example_different_goals():
    """Example 7: Compare different analysis goals."""
    print("\n" + "="*60)
    print("EXAMPLE 7: Analysis Goals Comparison")
    print("="*60)
    
    df = generate_demo_data()
    
    goals = ['discovery', 'balanced', 'validation']
    results = {}
    
    for goal in goals:
        result = optimize_thresholds(df, goal=goal, verbose=False)
        results[goal] = result
        print(f"\n{goal.upper()}:")
        print(f"  LogFC cutoff: {result.logfc_cutoff:.3f}")
        print(f"  Padj method: {result.padj_method}")
        print(f"  DE genes: {result.n_significant_optimized}")
    
    return results


def main():
    """Run all examples."""
    print("="*60)
    print("ADAPTIVE THRESHOLD OPTIMIZER - EXAMPLES")
    print("="*60)
    
    # Run examples
    example_basic_usage()
    example_full_control()
    example_comparison()
    example_adjustment_methods()
    example_visualization()
    example_different_goals()
    
    print("\n" + "="*60)
    print("All examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
