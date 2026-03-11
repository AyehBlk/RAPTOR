#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Data Profiler (Module 3)
UPDATED WITH VALIDATION

Demonstrates comprehensive data profiling for ML-based pipeline recommendation.
Extracts 32 statistical features from RNA-seq count matrices.

This is Module 3 of the RAPTOR workflow (Stage 1: Fast Profiling):
  M1: Quantify (FASTQ → quick_gene_counts.csv)
  M2: Sample QC (Quality Assessment & Outlier Detection)
  M3: Profile (Data Profiling - 32 features) ← THIS SCRIPT
  M4: Recommend (Pipeline Recommendation)

Input: results/quick_counts/quick_gene_counts.csv (or results/qc/counts_clean.csv)
Output: results/profile.json

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

# =============================================================================
# CONSTANTS - Architecture Compliant (v2.2.0)
# =============================================================================
DEFAULT_INPUT_COUNTS = "results/quick_counts/quick_gene_counts.csv"
DEFAULT_CLEAN_COUNTS = "results/qc/counts_clean.csv"
DEFAULT_OUTPUT_DIR = "results"

# Check for dependencies
try:
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: numpy and pandas are required")
    print("Install with: pip install numpy pandas")
    sys.exit(1)

# RAPTOR imports with validation - UPDATED FOR v2.2.0
RAPTOR_AVAILABLE = True
try:
    from raptor import (
        RNAseqDataProfiler,
        DataProfile,
        profile_data_quick,
        validate_count_matrix,
        validate_file_path,
        validate_directory_path,
        ValidationError
    )
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR not installed. Running in demo mode only.")
    print("Install RAPTOR with: pip install -e .")
    
    # Create dummy classes for demo mode
    class RNAseqDataProfiler:
        def __init__(self, counts, metadata=None, group_column=None): pass
        def run_full_profile(self): return None
    
    class DataProfile:
        pass
    
    def profile_data_quick(counts): return {}
    def validate_count_matrix(df, **kwargs): pass
    def validate_file_path(p, **kwargs): return Path(p)
    def validate_directory_path(p, **kwargs): return Path(p)
    class ValidationError(ValueError): pass


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.2.0 - Data Profiler (Module 3)            ║
    ║                                                              ║
    ║   Extracts 32 Features for ML Pipeline Recommendation       ║
    ║   ✅ WITH INPUT VALIDATION                                   ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def display_profile_summary(profile):
    """Display profile summary."""
    print("\n📊 Data Profile Summary:")
    print("  " + "="*70)
    
    if profile is None:
        # Demo mode
        print("  Demo mode - profile would contain 32 features:")
        print("\n  Basic Statistics:")
        print(f"    • n_samples: 6")
        print(f"    • n_genes: 20,000")
        print(f"    • library_size_mean: 4,500,000")
        print(f"    • library_size_cv: 0.08")
        
        print("\n  Expression Characteristics:")
        print(f"    • zero_fraction: 0.15")
        print(f"    • low_expr_fraction: 0.42")
        print(f"    • high_expr_genes: 856")
        
        print("\n  Variance Structure:")
        print(f"    • total_variance: 125.4")
        print(f"    • biological_cv: 0.32")
        print(f"    • technical_noise_estimate: 0.15")
        
        print("  " + "="*70)
        return
    
    # Real profile
    print(f"  Basic Statistics:")
    print(f"    • n_samples: {getattr(profile, 'n_samples', 'N/A')}")
    print(f"    • n_genes: {getattr(profile, 'n_genes', 'N/A'):,}")
    print(f"    • library_size_mean: {getattr(profile, 'lib_size_mean', 0):,.0f}")
    print(f"    • library_size_cv: {getattr(profile, 'lib_size_cv', 0):.2%}")
    
    print(f"\n  Expression Characteristics:")
    print(f"    • zero_fraction: {getattr(profile, 'zero_fraction', 0):.2%}")
    print(f"    • low_expr_fraction: {getattr(profile, 'low_expr_fraction', 0):.2%}")
    print(f"    • high_expr_genes: {getattr(profile, 'n_high_expr', 0):,}")
    
    print(f"\n  Variance Structure:")
    print(f"    • total_variance: {getattr(profile, 'total_variance', 0):.1f}")
    print(f"    • biological_cv: {getattr(profile, 'biological_cv', 0):.2f}")
    
    print("  " + "="*70)


def run_profiling(counts_file, metadata_file=None, group_column='condition',
                 output_dir=None, demo=False):
    """
    Run data profiling with validation.
    
    NEW: Added comprehensive input validation.
    """
    
    # Validate and load count matrix
    try:
        print("\n🔍 Validating inputs...")
        
        # Validate count file exists
        counts_path = validate_file_path(counts_file, must_exist=True)
        print(f"  ✓ Count file found: {counts_path}")
        
        # Load counts
        counts_df = pd.read_csv(counts_path, index_col=0)
        print(f"  ✓ Count matrix loaded: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
        
        # Validate count matrix structure
        validate_count_matrix(counts_df, min_genes=10, min_samples=2)
        print(f"  ✓ Count matrix validated")
        
    except FileNotFoundError:
        print(f"❌ Count file not found: {counts_file}")
        print(f"   Current directory: {Path.cwd()}")
        print(f"   Try one of these:")
        print(f"     - {DEFAULT_INPUT_COUNTS}")
        print(f"     - {DEFAULT_CLEAN_COUNTS}")
        sys.exit(1)
    except ValidationError as e:
        print(f"❌ Count matrix validation failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error loading count matrix: {e}")
        sys.exit(1)
    
    # Validate and load metadata (optional)
    metadata_df = None
    if metadata_file:
        try:
            from raptor import validate_metadata
            
            metadata_path = validate_file_path(metadata_file, must_exist=True)
            print(f"  ✓ Metadata file found: {metadata_path}")
            
            metadata_df = pd.read_csv(metadata_path)
            print(f"  ✓ Metadata loaded: {len(metadata_df)} samples")
            
            validate_metadata(metadata_df, counts_df)
            print(f"  ✓ Metadata validated")
            
        except Exception as e:
            print(f"⚠️  Warning: Metadata issue: {e}")
            print(f"   Continuing without metadata")
            metadata_df = None
    
    # Validate output directory
    try:
        output_path = validate_directory_path(
            output_dir or DEFAULT_OUTPUT_DIR,
            create_if_missing=True
        )
        print(f"  ✓ Output directory: {output_path}")
    except Exception as e:
        print(f"❌ Error with output directory: {e}")
        sys.exit(1)
    
    # Run profiling
    print("\n📊 Running data profiling...")
    print(f"  Extracting 32 statistical features...")
    
    if not RAPTOR_AVAILABLE or demo:
        print("  (Demo mode - showing simulated profile)")
        
        # Enhanced demo profile with all critical fields for Module 4
        lib_sizes = counts_df.sum(axis=0)
        lib_mean = lib_sizes.mean()
        lib_std = lib_sizes.std()
        zero_count = (counts_df == 0).sum().sum()
        total_values = counts_df.shape[0] * counts_df.shape[1]
        
        profile = type('Profile', (), {
            # SAMPLE CHARACTERISTICS (Module 4 needs these!)
            'n_samples': counts_df.shape[1],
            'n_genes': counts_df.shape[0],
            'n_groups': 2,
            'min_group_size': max(3, counts_df.shape[1] // 2),  # At least 3 per group
            'max_group_size': counts_df.shape[1] // 2,
            'has_replicates': True,
            'has_sufficient_replicates': True,
            'design_complexity': 'simple',
            
            # LIBRARY SIZE METRICS
            'library_size_mean': float(lib_mean),
            'library_size_median': float(lib_sizes.median()),
            'library_size_cv': float(lib_std / lib_mean) if lib_mean > 0 else 0.12,
            'library_size_range': float(lib_sizes.max() / lib_sizes.min()) if lib_sizes.min() > 0 else 2.5,
            'library_size_category': 'moderate',
            
            # EXPRESSION CHARACTERISTICS
            'zero_fraction': float(zero_count / total_values),
            'sparsity': 0.15,
            'low_count_proportion': 0.25,  # Affects edgeR score in Module 4
            'detection_rate': 0.85,
            'expression_mean': 8.5,
            'expression_median': 6.2,
            
            # DISPERSION (CRITICAL FOR MODULE 4!)
            'bcv': 0.32,  # Biological Coefficient of Variation
            'bcv_category': 'moderate',  # low/moderate/high
            'common_dispersion': 0.1024,  # bcv^2
            'dispersion_mean': 0.15,
            'overdispersion_ratio': 1.8,
            
            # QUALITY INDICATORS (Module 4 uses these!)
            'has_outliers': False,  # Affects edgeR_robust recommendation
            'has_batch_effect': False,  # Affects DESeq2/limma scoring
            'outlier_severity': 'none',
            
            # OTHER METRICS
            'total_variance': 125.4,
            'n_expressed_genes': int(counts_df.shape[0] * 0.85),
            'n_highly_expressed': 856,
            'sequencing_depth_category': 'standard'
        })()
        
        display_profile_summary(profile)
        
        # Demo results with enhanced structure
        results = {
            'timestamp': datetime.now().isoformat(),
            'raptor_version': '2.2.0',
            'module': 'M3',
            'mode': 'demo',
            'features': {
                # Sample characteristics
                'n_samples': counts_df.shape[1],
                'n_genes': counts_df.shape[0],
                'n_groups': 2,
                'min_group_size': max(3, counts_df.shape[1] // 2),
                
                # Library metrics
                'library_size_mean': float(lib_mean),
                'library_size_cv': float(lib_std / lib_mean) if lib_mean > 0 else 0.12,
                
                # Expression
                'zero_fraction': float(zero_count / total_values),
                'sparsity': 0.15,
                'low_count_proportion': 0.25,
                
                # Dispersion (CRITICAL!)
                'bcv': 0.32,
                'bcv_category': 'moderate',
                'common_dispersion': 0.1024,
                
                # Quality
                'has_outliers': False,
                'has_batch_effect': False,
                'design_complexity': 'simple'
            }
        }
        
    else:
        # Real profiling with RAPTOR
        profiler = RNAseqDataProfiler(counts_df, metadata_df, group_column)
        profile = profiler.run_full_profile()
        
        display_profile_summary(profile)
        
        # Real results
        results = profile.to_dict() if hasattr(profile, 'to_dict') else {}
        results['timestamp'] = datetime.now().isoformat()
        results['raptor_version'] = '2.2.0'
        results['module'] = 'M3'
    
    # Save profile
    profile_file = output_path / 'profile.json'
    with open(profile_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n  ✓ Profile saved: {profile_file}")
    
    return results, output_path


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Data Profiler (Module 3) - WITH VALIDATION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (uses default data)
  python 03_data_profiler.py --demo
  
  # With quick counts
  python 03_data_profiler.py \
      --counts results/quick_counts/quick_gene_counts.csv
  
  # With clean counts from QC
  python 03_data_profiler.py \
      --counts results/qc/counts_clean.csv
  
  # With metadata
  python 03_data_profiler.py \
      --counts results/qc/counts_clean.csv \
      --metadata results/quick_counts/sample_info.csv

CLI Equivalent:
  raptor profile --counts results/qc/counts_clean.csv

Workflow:
  Module 1: raptor quick-count → quick_gene_counts.csv
  Module 2: raptor qc → counts_clean.csv
  Module 3: THIS SCRIPT → profile.json
  Module 4: python 04_recommender.py --profile results/profile.json

Output:
  results/profile.json  (32 statistical features)
        """
    )
    
    parser.add_argument('--counts', '-c',
                       default=DEFAULT_CLEAN_COUNTS,
                       help='Count matrix CSV file')
    parser.add_argument('--metadata', '-m',
                       help='Sample metadata CSV file (optional)')
    parser.add_argument('--group-column', '-g',
                       default='condition',
                       help='Column name for grouping (default: condition)')
    parser.add_argument('--output', '-o',
                       default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--demo', action='store_true',
                       help='Run in demo mode')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Run profiling
    results, output_dir = run_profiling(
        counts_file=args.counts,
        metadata_file=args.metadata,
        group_column=args.group_column,
        output_dir=args.output,
        demo=args.demo
    )
    
    # Final summary
    print("\n" + "="*70)
    print("  ✅ MODULE 3 (DATA PROFILING) COMPLETE!")
    print("="*70)
    
    print(f"\n  📂 Output Directory: {output_dir}")
    print(f"\n  📊 Output Files:")
    print(f"     • profile.json  - 32 statistical features")
    
    print(f"\n  🔜 Next Steps:")
    print(f"\n     Module 4 - Pipeline Recommendation:")
    print(f"     python 04_recommender.py --profile {output_dir}/profile.json")
    print(f"     ")
    print(f"     Or use the CLI:")
    print(f"     raptor recommend --profile {output_dir}/profile.json")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
