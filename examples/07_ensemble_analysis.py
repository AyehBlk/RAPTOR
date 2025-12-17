#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Example Script: Ensemble Analysis

Demonstrates multi-pipeline consensus differential expression analysis:
- Weighted voting across pipelines
- Rank aggregation methods
- Intersection (conservative) approach
- Union (liberal) approach
- Pipeline concordance analysis

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
    from raptor.ensemble_analysis import EnsembleAnalyzer
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - Ensemble Analysis                  ║
    ║                                                              ║
    ║   Multi-Pipeline Consensus Differential Expression          ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_pipeline_results(n_genes=1000, seed=42):
    """Generate demonstration pipeline results for testing."""
    np.random.seed(seed)
    
    gene_names = [f'Gene_{i+1:05d}' for i in range(n_genes)]
    
    # Define true DE genes (20% of genes)
    n_de = n_genes // 5
    true_de = set(gene_names[:n_de])
    
    pipeline_results = {}
    
    # Pipeline 1: DESeq2 - high accuracy
    deseq2_de = set(np.random.choice(list(true_de), int(n_de * 0.85), replace=False))
    deseq2_fp = set(np.random.choice(gene_names[n_de:], int(n_de * 0.08), replace=False))
    deseq2_all = deseq2_de | deseq2_fp
    
    pipeline_results['DESeq2'] = pd.DataFrame({
        'log2FoldChange': [np.random.normal(2 if g in true_de else 0, 0.5) for g in gene_names],
        'pvalue': [0.001 if g in deseq2_all else np.random.uniform(0.1, 1) for g in gene_names],
        'padj': [0.005 if g in deseq2_all else np.random.uniform(0.1, 1) for g in gene_names]
    }, index=gene_names)
    
    # Pipeline 2: edgeR - slightly different results
    edger_de = set(np.random.choice(list(true_de), int(n_de * 0.82), replace=False))
    edger_fp = set(np.random.choice(gene_names[n_de:], int(n_de * 0.10), replace=False))
    edger_all = edger_de | edger_fp
    
    pipeline_results['edgeR'] = pd.DataFrame({
        'log2FoldChange': [np.random.normal(2 if g in true_de else 0, 0.5) for g in gene_names],
        'pvalue': [0.001 if g in edger_all else np.random.uniform(0.1, 1) for g in gene_names],
        'padj': [0.005 if g in edger_all else np.random.uniform(0.1, 1) for g in gene_names]
    }, index=gene_names)
    
    # Pipeline 3: limma-voom
    limma_de = set(np.random.choice(list(true_de), int(n_de * 0.80), replace=False))
    limma_fp = set(np.random.choice(gene_names[n_de:], int(n_de * 0.12), replace=False))
    limma_all = limma_de | limma_fp
    
    pipeline_results['limma'] = pd.DataFrame({
        'log2FoldChange': [np.random.normal(2 if g in true_de else 0, 0.5) for g in gene_names],
        'pvalue': [0.001 if g in limma_all else np.random.uniform(0.1, 1) for g in gene_names],
        'padj': [0.005 if g in limma_all else np.random.uniform(0.1, 1) for g in gene_names]
    }, index=gene_names)
    
    return pipeline_results, true_de


def simple_ensemble_analysis(pipeline_results, method='weighted_vote', min_pipelines=2, 
                              alpha=0.05, lfc_threshold=0.0):
    """Simple ensemble analysis implementation for demo mode."""
    
    # Get DE genes from each pipeline
    de_gene_sets = {}
    for pipeline, df in pipeline_results.items():
        if 'padj' in df.columns:
            significant = (df['padj'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
        else:
            significant = (df['pvalue'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
        de_gene_sets[pipeline] = set(df.index[significant])
    
    # Calculate consensus based on method
    if method == 'intersection':
        # Genes significant in ALL pipelines
        consensus_genes = set.intersection(*de_gene_sets.values())
    elif method == 'union':
        # Genes significant in ANY pipeline
        consensus_genes = set.union(*de_gene_sets.values())
    else:  # weighted_vote or rank_aggregation
        # Count votes for each gene
        all_genes = set.union(*de_gene_sets.values())
        gene_votes = {}
        for gene in all_genes:
            votes = sum(1 for de_set in de_gene_sets.values() if gene in de_set)
            if votes >= min_pipelines:
                gene_votes[gene] = votes
        consensus_genes = set(gene_votes.keys())
    
    # Calculate concordance
    from itertools import combinations
    concordances = []
    for p1, p2 in combinations(de_gene_sets.keys(), 2):
        s1, s2 = de_gene_sets[p1], de_gene_sets[p2]
        if len(s1.union(s2)) > 0:
            jaccard = len(s1.intersection(s2)) / len(s1.union(s2))
            concordances.append(jaccard)
    mean_concordance = np.mean(concordances) if concordances else 0.0
    
    # Create consensus DataFrame
    consensus_list = []
    for gene in consensus_genes:
        lfc_values = []
        pval_values = []
        n_pipelines = 0
        
        for pipeline, df in pipeline_results.items():
            if gene in de_gene_sets[pipeline]:
                lfc_values.append(df.loc[gene, 'log2FoldChange'])
                pval_values.append(df.loc[gene, 'pvalue'])
                n_pipelines += 1
        
        consensus_list.append({
            'gene_id': gene,
            'n_pipelines': n_pipelines,
            'mean_log2FC': np.mean(lfc_values),
            'std_log2FC': np.std(lfc_values) if len(lfc_values) > 1 else 0,
            'min_pvalue': min(pval_values)
        })
    
    consensus_df = pd.DataFrame(consensus_list)
    if len(consensus_df) > 0:
        consensus_df = consensus_df.sort_values('n_pipelines', ascending=False)
    
    return {
        'consensus_genes': consensus_df,
        'n_consensus': len(consensus_genes),
        'method': method,
        'concordance': mean_concordance,
        'n_pipelines': len(pipeline_results),
        'pipeline_names': list(pipeline_results.keys()),
        'per_pipeline_counts': {p: len(s) for p, s in de_gene_sets.items()},
        'full_agreement_genes': len([g for g in consensus_genes 
                                     if all(g in s for s in de_gene_sets.values())])
    }


def display_ensemble_results(results, pipeline_results):
    """Display ensemble analysis results with formatting."""
    
    print("\n" + "="*70)
    print("  🦖 ENSEMBLE ANALYSIS RESULTS")
    print("="*70)
    
    # Method summary
    print(f"\n  📊 Analysis Configuration:")
    print(f"     Method:     {results['method']}")
    print(f"     Pipelines:  {results['n_pipelines']} ({', '.join(results['pipeline_names'])})")
    
    # Per-pipeline DE genes
    print(f"\n  📈 DE Genes per Pipeline:")
    for pipeline, count in results['per_pipeline_counts'].items():
        bar_len = min(40, count // 5)
        bar = '█' * bar_len
        print(f"     {pipeline:15s} [{bar:40s}] {count}")
    
    # Concordance
    concordance = results['concordance']
    bar_len = int(concordance * 40)
    bar = '█' * bar_len + '░' * (40 - bar_len)
    print(f"\n  🤝 Pipeline Concordance:")
    print(f"     Jaccard Similarity: [{bar}] {concordance:.1%}")
    
    if concordance > 0.7:
        print(f"     → High agreement between pipelines")
    elif concordance > 0.5:
        print(f"     → Moderate agreement between pipelines")
    else:
        print(f"     → Low agreement - consider investigating differences")
    
    # Consensus genes
    print(f"\n  🎯 Consensus Results:")
    print(f"     Total consensus genes:     {results['n_consensus']}")
    print(f"     Full agreement genes:      {results['full_agreement_genes']}")
    
    # False positive reduction estimate
    max_pipeline = max(results['per_pipeline_counts'].values())
    if max_pipeline > 0:
        reduction = 1 - (results['n_consensus'] / max_pipeline)
        print(f"     Potential FP reduction:    {reduction:.1%}")
    
    # Top consensus genes
    if len(results['consensus_genes']) > 0:
        print(f"\n  🔝 Top 10 Consensus Genes:")
        print("  " + "-"*66)
        print(f"     {'Gene':15s} {'Pipelines':10s} {'Mean LFC':12s} {'Std LFC':10s}")
        print("  " + "-"*66)
        
        top_genes = results['consensus_genes'].head(10)
        for _, row in top_genes.iterrows():
            print(f"     {row['gene_id']:15s} {row['n_pipelines']:>5d}      "
                  f"{row['mean_log2FC']:>8.2f}    {row['std_log2FC']:>8.3f}")


def run_ensemble_analysis(results_dir=None, pipeline_files=None, method='weighted_vote',
                          min_pipelines=2, alpha=0.05, lfc_threshold=0.0, demo=False):
    """Run ensemble analysis on pipeline results."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated pipeline results...")
        pipeline_results, true_de = generate_demo_pipeline_results()
        
        print(f"\n📊 Demo pipeline results generated:")
        for name, df in pipeline_results.items():
            n_sig = (df['padj'] < alpha).sum()
            print(f"   {name}: {n_sig} significant genes")
        
        results = simple_ensemble_analysis(
            pipeline_results, 
            method=method,
            min_pipelines=min_pipelines,
            alpha=alpha,
            lfc_threshold=lfc_threshold
        )
        
    else:
        # Load pipeline results
        pipeline_results = {}
        
        if results_dir:
            results_path = Path(results_dir)
            for csv_file in results_path.glob('*.csv'):
                pipeline_name = csv_file.stem
                df = pd.read_csv(csv_file, index_col=0)
                pipeline_results[pipeline_name] = df
                print(f"   Loaded: {pipeline_name} ({len(df)} genes)")
        
        elif pipeline_files:
            for filepath in pipeline_files:
                pipeline_name = Path(filepath).stem
                df = pd.read_csv(filepath, index_col=0)
                pipeline_results[pipeline_name] = df
                print(f"   Loaded: {pipeline_name} ({len(df)} genes)")
        
        if len(pipeline_results) < 2:
            print("ERROR: At least 2 pipeline results are required")
            sys.exit(1)
        
        # Run ensemble analysis
        print(f"\n🔍 Running ensemble analysis ({method})...")
        
        ensemble = EnsembleAnalyzer(min_pipelines=min_pipelines)
        results = ensemble.combine_results(
            pipeline_results,
            method=method,
            alpha=alpha,
            lfc_threshold=lfc_threshold
        )
    
    # Display results
    display_ensemble_results(results, pipeline_results)
    
    # Prepare output
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'method': method,
        'parameters': {
            'min_pipelines': min_pipelines,
            'alpha': alpha,
            'lfc_threshold': lfc_threshold
        },
        'summary': {
            'n_consensus_genes': results['n_consensus'],
            'concordance': results['concordance'],
            'full_agreement_genes': results['full_agreement_genes'],
            'per_pipeline_counts': results['per_pipeline_counts']
        }
    }
    
    return output, results


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Ensemble Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze pipeline results in a directory
  python 07_ensemble_analysis.py --results-dir pipeline_results/
  
  # Specify individual pipeline files
  python 07_ensemble_analysis.py --pipelines deseq2.csv edger.csv limma.csv
  
  # Use intersection method (conservative)
  python 07_ensemble_analysis.py --results-dir results/ --method intersection
  
  # Require at least 3 pipelines to agree
  python 07_ensemble_analysis.py --results-dir results/ --min-pipelines 3
  
  # Demo mode (no data required)
  python 07_ensemble_analysis.py --demo

Available Methods:
  weighted_vote    - Weighted voting based on pipeline confidence (default)
  rank_aggregation - Combine gene rankings across pipelines
  intersection     - Genes significant in ALL pipelines (conservative)
  union            - Genes significant in ANY pipeline (liberal)
        """
    )
    
    parser.add_argument('--results-dir', help='Directory containing pipeline result CSVs')
    parser.add_argument('--pipelines', nargs='+', help='Individual pipeline result files')
    parser.add_argument('--method', choices=['weighted_vote', 'rank_aggregation', 
                                              'intersection', 'union'],
                        default='weighted_vote', help='Ensemble combination method')
    parser.add_argument('--min-pipelines', type=int, default=2,
                        help='Minimum pipelines that must agree')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance threshold')
    parser.add_argument('--lfc-threshold', type=float, default=0.0, 
                        help='Log2 fold-change threshold')
    parser.add_argument('--output', '-o', help='Output JSON file')
    parser.add_argument('--output-genes', help='Output consensus genes CSV')
    parser.add_argument('--demo', action='store_true', help='Run in demo mode')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs
    if not args.demo and not args.results_dir and not args.pipelines:
        print("ERROR: Either --results-dir, --pipelines, or --demo is required")
        parser.print_help()
        sys.exit(1)
    
    # Run analysis
    output, results = run_ensemble_analysis(
        results_dir=args.results_dir,
        pipeline_files=args.pipelines,
        method=args.method,
        min_pipelines=args.min_pipelines,
        alpha=args.alpha,
        lfc_threshold=args.lfc_threshold,
        demo=args.demo
    )
    
    # Save output
    output_file = args.output or 'ensemble_results.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\n💾 Summary saved to: {output_file}")
    
    # Save consensus genes if requested
    if args.output_genes or args.demo:
        genes_file = args.output_genes or 'consensus_genes.csv'
        if len(results['consensus_genes']) > 0:
            results['consensus_genes'].to_csv(genes_file, index=False)
            print(f"💾 Consensus genes saved to: {genes_file}")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
