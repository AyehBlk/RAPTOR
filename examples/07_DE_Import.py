#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: DE Import (Module 7)

Demonstrates importing differential expression results from R analysis tools
(DESeq2, edgeR, limma-voom, Wilcoxon) into standardized RAPTOR format.

This is Module 7 of the RAPTOR workflow (Stage 3: DE Analysis):
  M1-M5: Quantification & Profiling
  M6: External R Analysis (DESeq2/edgeR/limma) → de_results.csv
  M7: Import DE Results (THIS SCRIPT) → DEResult object
  M8: Parameter Optimization
  M9: Ensemble Analysis
  M10: Biomarker Discovery

Input: DE results CSV file from R
Output: results/de_imported/
    - de_standardized.csv  (standardized results)
    - de_significant.csv   (significant genes only)
    - de_summary.json      (summary statistics)
    - de_result.pkl        (DEResult object for M8-M10)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONSTANTS - Architecture Compliant (v2.2.0)
# =============================================================================
DEFAULT_OUTPUT_DIR = "results/de_imported"
DEFAULT_FDR = 0.05
DEFAULT_LFC = 0.0

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
    from raptor.de_import import (
        import_de_results,
        DEResult
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
            self.n_significant = len(results_df[results_df['is_significant']])
            self.n_up = len(results_df[(results_df['is_significant']) & (results_df['direction'] == 'up')])
            self.n_down = len(results_df[(results_df['is_significant']) & (results_df['direction'] == 'down')])
        
        def summary(self):
            return f"""
╔══════════════════════════════════════════════════════════╗
║  🦖 RAPTOR DE Results Summary                            ║
╠══════════════════════════════════════════════════════════╣

  Pipeline: {self.pipeline}
  Total genes tested: {self.n_genes:,}

  Significant genes: {self.n_significant:,} ({100*self.n_significant/self.n_genes:.1f}%)
    ↑ Upregulated:   {self.n_up:,}
    ↓ Downregulated: {self.n_down:,}

  Thresholds:
    FDR: {self.parameters.get('fdr_threshold', 0.05)}
    LFC: {self.parameters.get('lfc_threshold', 0.0)}

╚══════════════════════════════════════════════════════════╝
"""
        
        def get_top_genes(self, n=10, by='adjusted_p_value', significant_only=False):
            df = self.results_df
            if significant_only:
                df = df[df['is_significant']]
            if by == 'adjusted_p_value':
                return df.nsmallest(n, 'adjusted_p_value')
            return df.head(n)
    
    def import_de_results(de_file, output_dir=DEFAULT_OUTPUT_DIR, 
                         pipeline='auto', fdr_threshold=0.05, 
                         lfc_threshold=0.0, gene_id_column=None):
        """Demo mode import_de_results"""
        print("Running in demo mode...")
        return None


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.2.0 - DE Import (Module 7)                ║
    ║                                                              ║
    ║   Import DE Results from R Analysis Tools                   ║
    ║   ✅ DESeq2 | edgeR | limma-voom | Wilcoxon                  ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_de_results(pipeline='DESeq2', n_genes=15000, seed=42):
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
    p_values[:n_de] = np.random.beta(0.5, 10, n_de)  # Low p-values for DE genes
    p_values[n_de:] = np.random.uniform(0.1, 1.0, n_genes - n_de)  # Higher p-values
    
    # Shuffle
    indices = np.random.permutation(n_genes)
    p_values = p_values[indices]
    
    # Adjusted p-values (FDR correction simulation)
    adjusted_p_values = np.minimum(p_values * 10, 1.0)
    
    # Create DataFrame
    if pipeline.upper() == 'DESEQ2':
        df = pd.DataFrame({
            'gene_id': gene_ids,
            'baseMean': base_mean,
            'log2FoldChange': log2_fold_change,
            'lfcSE': np.abs(log2_fold_change) * 0.3,
            'stat': log2_fold_change / 0.3,
            'pvalue': p_values,
            'padj': adjusted_p_values
        })
    elif pipeline.upper() == 'EDGER':
        df = pd.DataFrame({
            'gene_id': gene_ids,
            'logFC': log2_fold_change,
            'logCPM': np.log2(base_mean + 1),
            'LR': np.random.chisquare(df=1, size=n_genes) * 2,
            'PValue': p_values,
            'FDR': adjusted_p_values
        })
    elif pipeline.upper() == 'LIMMA':
        df = pd.DataFrame({
            'gene_id': gene_ids,
            'logFC': log2_fold_change,
            'AveExpr': np.log2(base_mean + 1),
            't': log2_fold_change / 0.3,
            'P.Value': p_values,
            'adj.P.Val': adjusted_p_values,
            'B': np.random.normal(0, 3, n_genes)
        })
    else:
        df = pd.DataFrame({
            'gene_id': gene_ids,
            'log2FoldChange': log2_fold_change,
            'pvalue': p_values,
            'padj': adjusted_p_values
        })
    
    return df


def display_import_summary(de_result):
    """Display import summary."""
    if de_result is None:
        print("  Demo mode - import would be performed here")
        return
    
    print(de_result.summary())
    
    # Top significant genes
    if de_result.n_significant > 0:
        print("\n📊 Top 10 Most Significant Genes:")
        print("  " + "="*70)
        top_genes = de_result.get_top_genes(n=10, by='adjusted_p_value', significant_only=True)
        
        for idx, (gene_id, row) in enumerate(top_genes.iterrows(), 1):
            direction = "↑" if row['direction'] == 'up' else "↓"
            print(f"  {idx:2d}. {gene_id:<15} {direction} "
                  f"LFC={row['log2_fold_change']:>6.2f}  "
                  f"FDR={row['adjusted_p_value']:.2e}")
        print("  " + "="*70)


def run_import(de_file, output_dir=None, pipeline='auto', 
               fdr_threshold=DEFAULT_FDR, lfc_threshold=DEFAULT_LFC,
               gene_id_column=None, demo=False):
    """
    Run DE results import (Module 7).
    
    Parameters
    ----------
    de_file : str
        Path to DE results CSV from R
    output_dir : str, optional
        Output directory
    pipeline : str
        Pipeline name: 'auto', 'deseq2', 'edger', 'limma', 'wilcoxon'
    fdr_threshold : float
        FDR threshold for significance
    lfc_threshold : float
        Log2FC threshold for significance
    gene_id_column : str, optional
        Column containing gene IDs (auto-detect if None)
    demo : bool
        Run in demo mode
    
    Returns
    -------
    DEResult or None
        Imported DE result object
    """
    output_path = Path(output_dir or DEFAULT_OUTPUT_DIR)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # =========================================================================
    # Demo Mode
    # =========================================================================
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode...")
        print("─" * 60)
        
        # Generate demo data
        print("  Generating demo DE results (DESeq2)...")
        demo_df = generate_demo_de_results(pipeline='DESeq2', n_genes=15000)
        
        # Save demo file
        demo_file = output_path / 'demo_de_results.csv'
        demo_df.to_csv(demo_file, index=False)
        print(f"  ✓ Saved demo file: {demo_file}")
        
        # Simulate standardization
        print("\n  Standardizing columns...")
        standardized_df = demo_df.copy()
        standardized_df = standardized_df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value',
            'baseMean': 'base_mean',
            'lfcSE': 'se'
        })
        standardized_df = standardized_df.set_index('gene_id')
        
        # Calculate significance
        print(f"  Calculating significance (FDR={fdr_threshold}, LFC={lfc_threshold})...")
        standardized_df['is_significant'] = (
            (standardized_df['adjusted_p_value'] < fdr_threshold) &
            (standardized_df['log2_fold_change'].abs() > lfc_threshold)
        )
        standardized_df['direction'] = 'unchanged'
        mask_up = (standardized_df['is_significant']) & (standardized_df['log2_fold_change'] > 0)
        mask_down = (standardized_df['is_significant']) & (standardized_df['log2_fold_change'] < 0)
        standardized_df.loc[mask_up, 'direction'] = 'up'
        standardized_df.loc[mask_down, 'direction'] = 'down'
        
        # Create DEResult object
        de_result = DEResult(
            results_df=standardized_df,
            pipeline='DESEQ2',
            parameters={
                'fdr_threshold': fdr_threshold,
                'lfc_threshold': lfc_threshold,
                'source_pipeline': 'deseq2'
            },
            metadata={
                'source_file': str(demo_file),
                'timestamp': datetime.now().isoformat(),
                'raptor_version': '2.2.0',
                'module': 'M7',
                'mode': 'demo'
            }
        )
        
        # Display summary
        display_import_summary(de_result)
        
        # Save outputs
        print("\n📁 Saving results to:", output_path)
        standardized_df.to_csv(output_path / 'de_standardized.csv')
        print("  ✓ de_standardized.csv")
        
        sig_df = standardized_df[standardized_df['is_significant']]
        sig_df.to_csv(output_path / 'de_significant.csv')
        print(f"  ✓ de_significant.csv ({len(sig_df):,} genes)")
        
        # Summary JSON
        import json
        summary = {
            'timestamp': datetime.now().isoformat(),
            'raptor_version': '2.2.0',
            'module': 'M7',
            'mode': 'demo',
            'pipeline': 'DESEQ2',
            'n_genes': de_result.n_genes,
            'n_significant': de_result.n_significant,
            'n_up': de_result.n_up,
            'n_down': de_result.n_down,
            'fdr_threshold': fdr_threshold,
            'lfc_threshold': lfc_threshold
        }
        with open(output_path / 'de_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        print("  ✓ de_summary.json")
        
        return de_result
    
    # =========================================================================
    # Real Import with RAPTOR
    # =========================================================================
    print("\n🚀 IMPORTING DE RESULTS")
    print("─" * 60)
    
    print(f"  Input file: {de_file}")
    print(f"  Pipeline: {pipeline}")
    print(f"  FDR threshold: {fdr_threshold}")
    print(f"  LFC threshold: {lfc_threshold}")
    
    try:
        de_result = import_de_results(
            de_file=de_file,
            output_dir=str(output_path),
            pipeline=pipeline,
            fdr_threshold=fdr_threshold,
            lfc_threshold=lfc_threshold,
            gene_id_column=gene_id_column
        )
        
        display_import_summary(de_result)
        
        return de_result
        
    except FileNotFoundError:
        print(f"\n❌ DE results file not found: {de_file}")
        print(f"   Current directory: {Path.cwd()}")
        print(f"   Expected file: {Path(de_file).absolute()}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error importing DE results: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 DE Import (Module 7)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (no data needed)
  python 07_DE_Import.py --demo
  
  # Import DESeq2 results (auto-detect)
  python 07_DE_Import.py --de-file data/deseq2_results.csv
  
  # Import with specific pipeline
  python 07_DE_Import.py --de-file data/de_results.csv --pipeline deseq2
  
  # Import with custom thresholds
  python 07_DE_Import.py \\
      --de-file data/de_results.csv \\
      --fdr-threshold 0.01 \\
      --lfc-threshold 1.0
  
  # Import edgeR results
  python 07_DE_Import.py --de-file data/edger_results.csv --pipeline edger
  
  # Import limma results
  python 07_DE_Import.py --de-file data/limma_results.csv --pipeline limma

CLI Equivalent:
  raptor import-de --de-file data/deseq2_results.csv

Workflow:
  M6: Run R analysis (DESeq2/edgeR/limma) → de_results.csv
  M7: THIS SCRIPT → de_result.pkl
  M8: raptor optimize --de-result results/de_imported/de_result.pkl
  M9: raptor ensemble --de-results [multiple_files]

Expected DE Results Format (from R):
  
  DESeq2:
    gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
  
  edgeR:
    gene_id, logFC, logCPM, LR, PValue, FDR
  
  limma-voom:
    gene_id, logFC, AveExpr, t, P.Value, adj.P.Val, B

Output Files:
  results/de_imported/
    ├── de_standardized.csv   (all genes, standardized columns)
    ├── de_significant.csv    (significant genes only)
    ├── de_summary.json       (summary statistics)
    └── de_result.pkl         (DEResult object for M8-M10)
        """
    )
    
    parser.add_argument('--de-file', '-f',
                       help='DE results CSV file from R analysis')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--pipeline', '-p',
                       choices=['auto', 'deseq2', 'edger', 'limma', 'wilcoxon'],
                       default='auto',
                       help='Pipeline name (auto-detect if not specified)')
    parser.add_argument('--fdr-threshold', type=float, default=DEFAULT_FDR,
                       help=f'FDR threshold for significance (default: {DEFAULT_FDR})')
    parser.add_argument('--lfc-threshold', type=float, default=DEFAULT_LFC,
                       help=f'Log2FC threshold for significance (default: {DEFAULT_LFC})')
    parser.add_argument('--gene-id-column', type=str,
                       help='Column containing gene IDs (auto-detect if not specified)')
    parser.add_argument('--demo', action='store_true',
                       help='Run in demo mode with simulated data')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs for real run
    if not args.demo and not args.de_file:
        print("ERROR: --de-file is required (or use --demo)")
        parser.print_help()
        sys.exit(1)
    
    # Run import
    de_result = run_import(
        de_file=args.de_file,
        output_dir=args.output,
        pipeline=args.pipeline,
        fdr_threshold=args.fdr_threshold,
        lfc_threshold=args.lfc_threshold,
        gene_id_column=args.gene_id_column,
        demo=args.demo
    )
    
    # Final summary
    print("\n" + "="*70)
    print("  ✅ MODULE 7 (DE IMPORT) COMPLETE!")
    print("="*70)
    
    output_dir = Path(args.output)
    
    print(f"\n  📂 Output Directory: {output_dir}")
    
    if de_result:
        print(f"\n  📊 Results:")
        print(f"     • Total genes: {de_result.n_genes:,}")
        print(f"     • Significant: {de_result.n_significant:,} "
              f"({100*de_result.n_significant/de_result.n_genes:.1f}%)")
        print(f"     • Upregulated: {de_result.n_up:,}")
        print(f"     • Downregulated: {de_result.n_down:,}")
    
    print(f"\n  📄 Output Files:")
    print(f"     • de_standardized.csv  - All genes, standardized format")
    print(f"     • de_significant.csv   - Significant genes only")
    print(f"     • de_summary.json      - Summary statistics")
    if RAPTOR_AVAILABLE:
        print(f"     • de_result.pkl        - DEResult object (for M8-M10)")
    
    print(f"\n  🔜 Next Steps:")
    print(f"\n     Module 8 - Parameter Optimization:")
    print(f"     python 08_Parameter_Optimization.py \\")
    print(f"         --de-result {output_dir}/de_result.pkl")
    print(f"")
    print(f"     Module 9 - Ensemble Analysis:")
    print(f"     python 09_Ensemble_Analysis.py \\")
    print(f"         --de-results {output_dir}/de_result.pkl [more_files]")
    print(f"")
    print(f"     Or use the CLI:")
    print(f"     raptor optimize --de-result {output_dir}/de_result.pkl")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
