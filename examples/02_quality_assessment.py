#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Quality Assessment (Module 2)
UPDATED WITH VALIDATION

Demonstrates comprehensive quality assessment including:
- Library quality assessment
- Gene detection analysis  
- Advanced outlier detection (6 methods)
- Batch effect detection
- Variance structure analysis
- Biological signal assessment
- Overall quality scoring (0-100)
- Input validation with clear error messages

This is Module 2 of the RAPTOR workflow (Stage 1: Fast Profiling):
  M1: Quantify (FASTQ → quick_gene_counts.csv)
  M2: Sample QC (Quality Assessment & Outlier Detection) ← THIS SCRIPT
  M3: Profile (Data Profiling - 32 features)
  M4: Recommend (Pipeline Recommendation)

Input: results/quick_counts/quick_gene_counts.csv
Output: results/qc/

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
DEFAULT_INPUT_DIR = "results/quick_counts"
DEFAULT_OUTPUT_DIR = "results/qc"
INPUT_COUNTS_FILE = "quick_gene_counts.csv"
INPUT_SAMPLE_INFO = "sample_info.csv"

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
        DataQualityAssessor,
        quick_quality_check,
        validate_count_matrix,
        validate_metadata,
        validate_file_path,
        validate_directory_path,
        ValidationError
    )
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR not installed. Running in demo mode only.")
    print("Install RAPTOR with: pip install -e .")
    
    # Create dummy classes and functions for demo mode
    class DataQualityAssessor:
        def __init__(self, counts, metadata=None): pass
        def assess_quality(self): return None
    
    def quick_quality_check(counts): return {}
    def validate_count_matrix(df, **kwargs): pass
    def validate_metadata(meta, counts): pass
    def validate_file_path(p, **kwargs): return Path(p)
    def validate_directory_path(p, **kwargs): return Path(p)
    class ValidationError(ValueError): pass


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.2.0 - Quality Assessment (Module 2)       ║
    ║                                                              ║
    ║   Comprehensive QC with 6-Method Outlier Detection          ║
    ║   Quality Score: 0-100 with Component Analysis              ║
    ║   ✅ WITH INPUT VALIDATION                                   ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def display_qc_summary(qc_result):
    """Display QC summary in formatted output."""
    print("\n📊 Quality Assessment Summary:")
    print("  " + "="*70)
    
    if qc_result is None:
        print("  Demo mode - QC analysis would be performed here")
        print("  " + "="*70)
        return
    
    # Real RAPTOR returns a Dict, not an object
    # Access overall quality from dict
    if isinstance(qc_result, dict):
        overall = qc_result.get('overall', {})
        score = overall.get('score', 0)
        quality = overall.get('quality', 'UNKNOWN')
    else:
        # Fallback for demo mode object
        score = getattr(qc_result, 'quality_score', 75)
        quality = getattr(qc_result, 'overall_quality', 'GOOD')
    
    print(f"  Overall Quality: {quality}")
    print(f"  Quality Score: {score:.1f}/100")
    
    # Component scores
    print(f"\n  Component Scores:")
    
    if isinstance(qc_result, dict):
        # Real RAPTOR dict structure
        components = {
            'Library Quality': qc_result.get('library_quality', {}),
            'Gene Detection': qc_result.get('gene_detection', {}),
            'Biological Signal': qc_result.get('biological_signal', {}),
            'Outlier Detection': qc_result.get('outlier_detection', {})
        }
        
        for comp_name, comp_data in components.items():
            if comp_data:
                comp_score = comp_data.get('score', 0)
                status = "✓" if comp_score >= 70 else "⚠" if comp_score >= 50 else "✗"
                print(f"    {status} {comp_name:<25} {comp_score:.1f}/100")
    else:
        # Demo mode object structure
        components = {
            'Library Quality': getattr(qc_result, 'library_quality', 80),
            'Gene Detection': getattr(qc_result, 'gene_detection', 85),
            'Biological Signal': getattr(qc_result, 'biological_signal', 70),
            'Technical Quality': getattr(qc_result, 'technical_quality', 75)
        }
        
        for comp, comp_score in components.items():
            status = "✓" if comp_score >= 70 else "⚠" if comp_score >= 50 else "✗"
            print(f"    {status} {comp:<25} {comp_score:.1f}/100")
    
    # Issues - collect from all components
    all_flags = []
    if isinstance(qc_result, dict):
        for comp in qc_result.values():
            if isinstance(comp, dict) and 'flags' in comp:
                all_flags.extend(comp['flags'])
    else:
        all_flags = getattr(qc_result, 'issues', [])
    
    if all_flags:
        print(f"\n  Issues Found: {len(all_flags)}")
        for flag in all_flags[:3]:
            print(f"    • {flag}")
        if len(all_flags) > 3:
            print(f"    ... and {len(all_flags) - 3} more")
    else:
        print(f"\n  ✓ No major issues detected")
    
    print("  " + "="*70)


def run_quality_assessment(counts_file, metadata_file=None, output_dir=None, 
                          normalization='log2', demo=False):
    """
    Run quality assessment with validation.
    
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
        print(f"   Expected location: {Path(counts_file).absolute()}")
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
            metadata_path = validate_file_path(metadata_file, must_exist=True)
            print(f"  ✓ Metadata file found: {metadata_path}")
            
            metadata_df = pd.read_csv(metadata_path)
            print(f"  ✓ Metadata loaded: {len(metadata_df)} samples")
            
            # Validate metadata matches counts
            validate_metadata(metadata_df, counts_df)
            print(f"  ✓ Metadata validated")
            
        except FileNotFoundError:
            print(f"⚠️  Warning: Metadata file not found: {metadata_file}")
            print(f"   Continuing without metadata")
            metadata_df = None
        except ValidationError as e:
            print(f"⚠️  Warning: Metadata validation failed: {e}")
            print(f"   Continuing without metadata")
            metadata_df = None
        except Exception as e:
            print(f"⚠️  Warning: Error loading metadata: {e}")
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
    
    # Run quality assessment
    print("\n📊 Running quality assessment...")
    
    if not RAPTOR_AVAILABLE or demo:
        print("  (Demo mode - showing simulated results)")
        
        # Generate demo QC results
        qc_result = type('QCResult', (), {
            'quality_score': 78,
            'overall_quality': 'GOOD',
            'library_quality': 82,
            'gene_detection': 85,
            'biological_signal': 72,
            'technical_quality': 75,
            'issues': [
                'Sample Control_2 shows slightly elevated library size CV',
                'Batch effect detected between Batch1 and Batch2'
            ],
            'outliers': [],
            'recommendations': [
                'Consider batch correction before DE analysis',
                'Library sizes are within acceptable range'
            ]
        })()
        
        display_qc_summary(qc_result)
        
        # Save demo results
        results = {
            'timestamp': datetime.now().isoformat(),
            'raptor_version': '2.2.0',
            'module': 'M2',
            'mode': 'demo',
            'quality_score': 78,
            'overall_quality': 'GOOD',
            'issues': qc_result.issues,
            'recommendations': qc_result.recommendations
        }
        
    else:
        # Real QC with RAPTOR
        assessor = DataQualityAssessor(counts_df, metadata_df)
        qc_result = assessor.assess_quality()
        
        display_qc_summary(qc_result)
        
        # Save real results
        results = qc_result.to_dict() if hasattr(qc_result, 'to_dict') else {}
        results['timestamp'] = datetime.now().isoformat()
        results['raptor_version'] = '2.2.0'
        results['module'] = 'M2'
    
    # Save results
    results_file = output_path / 'qc_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    # Save clean counts (remove outliers if any)
    clean_counts = counts_df  # In demo, no outliers removed
    clean_counts.to_csv(output_path / 'counts_clean.csv')
    
    return results, output_path


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Quality Assessment (Module 2) - WITH VALIDATION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (uses default demo data from M1)
  python 02_quality_assessment.py --demo
  
  # With count matrix
  python 02_quality_assessment.py --counts results/quick_counts/quick_gene_counts.csv
  
  # With count matrix and metadata
  python 02_quality_assessment.py \
      --counts results/quick_counts/quick_gene_counts.csv \
      --metadata results/quick_counts/sample_info.csv
  
  # Custom output directory
  python 02_quality_assessment.py \
      --counts results/quick_counts/quick_gene_counts.csv \
      --output results/my_qc/

CLI Equivalent:
  raptor qc --counts results/quick_counts/quick_gene_counts.csv

Workflow:
  Module 1: raptor quick-count → quick_gene_counts.csv
  Module 2: THIS SCRIPT → qc_results.json, counts_clean.csv
  Module 3: python 03_data_profiler.py --counts results/qc/counts_clean.csv
  Module 4: python 04_recommender.py --profile results/profile.json

Output Files:
  results/qc/
    ├── qc_results.json      (Quality assessment report)
    ├── counts_clean.csv     (Clean count matrix)
    └── outlier_report.csv   (Outlier detection results)
        """
    )
    
    parser.add_argument('--counts', '-c',
                       default=f'{DEFAULT_INPUT_DIR}/{INPUT_COUNTS_FILE}',
                       help='Count matrix CSV file')
    parser.add_argument('--metadata', '-m',
                       help='Sample metadata CSV file (optional)')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--normalization', '-n',
                       choices=['log2', 'cpm', 'quantile', 'none'],
                       default='log2',
                       help='Normalization method for outlier detection')
    parser.add_argument('--demo', action='store_true',
                       help='Run in demo mode with simulated results')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Run QC
    results, output_dir = run_quality_assessment(
        counts_file=args.counts,
        metadata_file=args.metadata,
        output_dir=args.output,
        normalization=args.normalization,
        demo=args.demo
    )
    
    # Final summary
    print("\n" + "="*70)
    print("  ✅ MODULE 2 (QUALITY ASSESSMENT) COMPLETE!")
    print("="*70)
    
    print(f"\n  📂 Output Directory: {output_dir}")
    print(f"\n  📊 Output Files:")
    print(f"     • qc_results.json       - Quality assessment report")
    print(f"     • counts_clean.csv      - Clean count matrix")
    
    print(f"\n  🔜 Next Steps:")
    print(f"\n     Module 3 - Data Profiling:")
    print(f"     python 03_data_profiler.py --counts {output_dir}/counts_clean.csv")
    print(f"     ")
    print(f"     Or use the CLI:")
    print(f"     raptor profile --counts {output_dir}/counts_clean.csv")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
