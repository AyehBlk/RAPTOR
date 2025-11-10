#!/usr/bin/env python3

"""
02_profile_data.py
Profile RNA-seq count data using RAPTOR's profiler module and get pipeline recommendations

This script integrates with raptor.profiler.RNAseqDataProfiler and 
raptor.recommender.PipelineRecommender to provide intelligent pipeline selection.

Usage: python 02_profile_data.py <count_matrix.csv> [options]

Author: Ayeh Bolouki
Organization: RAPTOR Project
License: MIT
"""

import sys
import argparse
import json
from pathlib import Path
import pandas as pd

# Try to import RAPTOR modules
try:
    from raptor.profiler import RNAseqDataProfiler
    from raptor.recommender import PipelineRecommender
    RAPTOR_AVAILABLE = True
except ImportError:
    RAPTOR_AVAILABLE = False
    print("Warning: RAPTOR package not found. Install with: pip install raptor-rnaseq")
    print("Falling back to standalone mode...\n")


def profile_with_raptor(counts_file, metadata_file=None, output_file='data_profile.json'):
    """
    Profile data using RAPTOR's RNAseqDataProfiler.
    
    Parameters
    ----------
    counts_file : str
        Path to count matrix CSV
    metadata_file : str, optional
        Path to metadata CSV
    output_file : str
        Output JSON file path
        
    Returns
    -------
    dict
        Profile data
    """
    print("=== RAPTOR Data Profiling ===\n")
    
    # Load data
    print(f"Loading count matrix: {counts_file}")
    counts = pd.read_csv(counts_file, index_col=0)
    print(f"  Loaded: {counts.shape[0]} genes × {counts.shape[1]} samples\n")
    
    metadata = None
    if metadata_file:
        print(f"Loading metadata: {metadata_file}")
        metadata = pd.read_csv(metadata_file)
        print(f"  Loaded: {len(metadata)} samples\n")
    
    # Create profiler
    profiler = RNAseqDataProfiler(counts, metadata)
    
    # Run full profile
    print("Running profile analysis...")
    profile = profiler.run_full_profile()
    
    print("\n=== Profile Summary ===")
    print(f"Samples: {profile['design']['n_samples']}")
    print(f"Genes: {profile['design']['n_genes']}")
    print(f"Mean library size: {profile['library_stats']['mean']:,.0f}")
    print(f"Library size CV: {profile['library_stats']['cv']:.2f}")
    print(f"Zero percentage: {profile['count_distribution']['zero_pct']:.1f}%")
    print(f"Biological variation (BCV): {profile['biological_variation']['bcv']:.3f}")
    print(f"Sequencing depth: {profile['sequencing']['depth_category']}")
    print(f"Data difficulty: {profile['summary']['difficulty']}")
    print(f"Recommended approach: {profile['summary']['recommended_approach']}")
    
    # Get pipeline recommendation
    print("\n=== Pipeline Recommendations ===")
    try:
        recommender = PipelineRecommender(profile)
        recommendation = recommender.recommend()
        
        print(f"\n🦖 RECOMMENDED: Pipeline {recommendation['pipeline_id']} - {recommendation['pipeline_name']}")
        print(f"\nReasoning:")
        for reason in recommendation['reasons']:
            print(f"  • {reason}")
        
        if 'alternatives' in recommendation:
            print(f"\nAlternative options:")
            for alt in recommendation['alternatives']:
                print(f"  • Pipeline {alt['pipeline_id']}: {alt['pipeline_name']}")
                print(f"    {alt['reason']}")
        
        # Add recommendation to profile
        profile['recommendation'] = recommendation
        
    except Exception as e:
        print(f"Could not generate recommendation: {e}")
        print("Profile data still available for manual review.")
    
    # Save profile
    output_path = Path(output_file)
    with open(output_path, 'w') as f:
        json.dump(profile, f, indent=2)
    
    print(f"\n✓ Profile saved to: {output_path}")
    
    return profile


def profile_standalone(counts_file, output_file='data_profile.json'):
    """
    Fallback profiling without RAPTOR package.
    
    Provides basic profiling when RAPTOR package is not installed.
    """
    import numpy as np
    from scipy import stats
    
    print("=== Standalone Data Profiling ===\n")
    print("Note: For full features, install RAPTOR: pip install raptor-rnaseq\n")
    
    # Load data
    print(f"Loading count matrix: {counts_file}")
    counts = pd.read_csv(counts_file, index_col=0)
    n_genes, n_samples = counts.shape
    print(f"  Loaded: {n_genes} genes × {n_samples} samples\n")
    
    print("Calculating basic statistics...")
    
    # Basic statistics
    library_sizes = counts.sum(axis=0)
    mean_lib_size = library_sizes.mean()
    cv_lib_size = library_sizes.std() / mean_lib_size
    
    # Zero inflation
    zero_pct = 100 * (counts == 0).sum().sum() / counts.size
    
    # Gene expression
    gene_means = counts.mean(axis=1)
    low_expr_genes = (gene_means < 10).sum()
    
    profile = {
        'n_samples': n_samples,
        'n_genes': n_genes,
        'mean_library_size': float(mean_lib_size),
        'library_size_cv': float(cv_lib_size),
        'zero_percentage': float(zero_pct),
        'low_expression_genes': int(low_expr_genes),
        'median_expression': float(gene_means.median()),
    }
    
    # Simple recommendations
    recommendations = []
    if n_samples < 3:
        recommendations.append("Pipeline 6 (NOISeq) - designed for no/few replicates")
        recommendations.append("Pipeline 7 (EBSeq) - optimized for small samples")
    elif n_samples >= 6:
        recommendations.append("Pipeline 1 (STAR-RSEM-DESeq2) - gold standard")
        recommendations.append("Pipeline 3 (Salmon-edgeR) - fast and accurate")
    else:
        recommendations.append("Pipeline 1 (STAR-RSEM-DESeq2) - robust for moderate samples")
    
    if mean_lib_size < 10e6:
        recommendations.append("Pipeline 3 (Salmon) or Pipeline 4 (Kallisto) - efficient for lower coverage")
    
    profile['recommendations'] = recommendations
    
    # Display
    print("\n=== Profile Summary ===")
    print(f"Samples: {n_samples}")
    print(f"Genes: {n_genes}")
    print(f"Mean library size: {mean_lib_size:,.0f}")
    print(f"Library size CV: {cv_lib_size:.2f}")
    print(f"Zero percentage: {zero_pct:.1f}%")
    print(f"Low expression genes: {low_expr_genes}")
    
    print("\n=== Basic Recommendations ===")
    for i, rec in enumerate(recommendations, 1):
        print(f"{i}. {rec}")
    
    # Save
    output_path = Path(output_file)
    with open(output_path, 'w') as f:
        json.dump(profile, f, indent=2)
    
    print(f"\n✓ Profile saved to: {output_path}")
    print("\nFor detailed profiling and recommendations, install RAPTOR:")
    print("  pip install raptor-rnaseq")
    
    return profile


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Profile RNA-seq count data for RAPTOR pipeline selection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic profiling
  python 02_profile_data.py counts.csv
  
  # With metadata for improved recommendations
  python 02_profile_data.py counts.csv --metadata metadata.csv
  
  # Custom output file
  python 02_profile_data.py counts.csv -o my_profile.json

For more information: https://github.com/AyehBlk/RAPTOR
        """
    )
    
    parser.add_argument('count_matrix', 
                       help='Path to count matrix CSV file (genes x samples)')
    parser.add_argument('-m', '--metadata', 
                       help='Path to metadata CSV file (optional)')
    parser.add_argument('-o', '--output', 
                       default='data_profile.json',
                       help='Output JSON file for profile (default: data_profile.json)')
    parser.add_argument('-v', '--verbose', 
                       action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Check if count file exists
    if not Path(args.count_matrix).exists():
        print(f"Error: Count matrix file not found: {args.count_matrix}")
        sys.exit(1)
    
    # Profile data
    try:
        if RAPTOR_AVAILABLE:
            profile = profile_with_raptor(
                args.count_matrix,
                args.metadata,
                args.output
            )
        else:
            profile = profile_standalone(
                args.count_matrix,
                args.output
            )
        
        print("\n✓ Data profiling complete!")
        
    except Exception as e:
        print(f"\nError during profiling: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
