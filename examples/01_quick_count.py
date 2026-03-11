#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Quick Quantification (Module 1)
UPDATED WITH VALIDATION

Demonstrates fast RNA-seq quantification with:
- Salmon pseudo-alignment
- Kallisto pseudo-alignment
- Automatic sample sheet creation
- Count matrix generation
- Quality metrics summary
- Input validation with clear error messages

This is Module 1 of the RAPTOR workflow (Stage 1: Fast Profiling):
  M1: Quantify (FASTQ → quick_gene_counts.csv) ← THIS SCRIPT
  M2: Sample QC (Quality Assessment & Outlier Detection)
  M3: Profile (Data Profiling - 32 features)
  M4: Recommend (Pipeline Recommendation)

Output Location: results/quick_counts/
Output Files:
  - quick_gene_counts.csv  (gene-level count matrix)
  - quick_tpm.csv          (TPM normalized matrix)
  - sample_info.csv        (sample metadata)

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
DEFAULT_OUTPUT_DIR = "results/quick_counts"
OUTPUT_COUNTS_FILE = "quick_gene_counts.csv"
OUTPUT_TPM_FILE = "quick_tpm.csv"
OUTPUT_SAMPLE_INFO = "sample_info.csv"

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
        validate_file_path,
        validate_directory_path,
        validate_positive_integer,
        ValidationError
    )
    # Try to import pipeline-specific modules
    try:
        from raptor.pipelines import SampleSheet, auto_detect_samples
    except ImportError:
        # Fallback if pipelines not fully implemented yet
        SampleSheet = None
        auto_detect_samples = None
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR not installed. Running in demo mode only.")
    print("Install RAPTOR with: pip install -e .")
    
    # Create dummy validation functions for demo mode
    def validate_file_path(p, **kwargs): return Path(p)
    def validate_directory_path(p, **kwargs): return Path(p)
    def validate_positive_integer(v, name): 
        if v < 1: raise ValueError(f"{name} must be positive")
    class ValidationError(ValueError): pass


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.2.0 - Quick Quantification (Module 1)     ║
    ║                                                              ║
    ║   Fast RNA-seq Quantification with Salmon or Kallisto       ║
    ║   FASTQ → Count Matrix in minutes                           ║
    ║   ✅ WITH INPUT VALIDATION                                   ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def print_workflow():
    """Print the RAPTOR workflow diagram."""
    print("""
    ┌─────────────────────────────────────────────────────────────────┐
    │                 RAPTOR v2.2.0 WORKFLOW                          │
    ├─────────────────────────────────────────────────────────────────┤
    │                                                                 │
    │  STAGE 1: Fast Profiling (M1-M4)                               │
    │  ═══════════════════════════════                               │
    │                                                                 │
    │  ┌──────────┐     ┌──────────┐     ┌─────────────────────┐     │
    │  │  FASTQ   │ ──► │  M1:     │ ──► │ quick_gene_counts   │     │
    │  │  files   │     │ Quantify │     │      .csv           │     │
    │  │          │     │          │     │ results/quick_counts│     │
    │  └──────────┘     └──────────┘     └─────────────────────┘     │
    │       ▲               │ ◄── YOU ARE HERE                       │
    │       │               ▼                                         │
    │  ┌──────────┐     ┌──────────┐     ┌──────────┐                │
    │  │  sample  │     │  M2:     │ ──► │  clean   │                │
    │  │  sheet   │     │ Sample   │     │  counts  │                │
    │  │   .csv   │     │   QC     │     │          │                │
    │  └──────────┘     └──────────┘     └──────────┘                │
    │                       │                                         │
    │                       ▼                                         │
    │                   ┌──────────┐                                  │
    │                   │  M3:     │                                  │
    │                   │ Profile  │                                  │
    │                   │(32 feat) │                                  │
    │                   └──────────┘                                  │
    │                       │                                         │
    │                       ▼                                         │
    │                   ┌──────────┐                                  │
    │                   │  M4:     │                                  │
    │                   │Recommend │                                  │
    │                   │          │                                  │
    │                   └──────────┘                                  │
    │                       │                                         │
    │                       ▼                                         │
    │  STAGE 2: Production Pipeline (M5)                             │
    │  ═════════════════════════════════                             │
    │                                                                 │
    │  STAGE 3: DE Analysis (M6-M10)                                 │
    │  ══════════════════════════════                                │
    │                                                                 │
    └─────────────────────────────────────────────────────────────────┘
    """)


def generate_demo_fastq_structure():
    """Generate demonstration FASTQ file structure."""
    return {
        'samples': [
            {'sample_id': 'Control_1', 'fastq_r1': 'Control_1_R1.fastq.gz', 'fastq_r2': 'Control_1_R2.fastq.gz'},
            {'sample_id': 'Control_2', 'fastq_r1': 'Control_2_R1.fastq.gz', 'fastq_r2': 'Control_2_R2.fastq.gz'},
            {'sample_id': 'Control_3', 'fastq_r1': 'Control_3_R1.fastq.gz', 'fastq_r2': 'Control_3_R2.fastq.gz'},
            {'sample_id': 'Treatment_1', 'fastq_r1': 'Treatment_1_R1.fastq.gz', 'fastq_r2': 'Treatment_1_R2.fastq.gz'},
            {'sample_id': 'Treatment_2', 'fastq_r1': 'Treatment_2_R1.fastq.gz', 'fastq_r2': 'Treatment_2_R2.fastq.gz'},
            {'sample_id': 'Treatment_3', 'fastq_r1': 'Treatment_3_R1.fastq.gz', 'fastq_r2': 'Treatment_3_R2.fastq.gz'},
        ],
        'read_type': 'paired-end',
        'n_samples': 6
    }


def generate_demo_counts(n_samples=6, n_genes=20000, seed=42):
    """Generate demonstration count matrix."""
    np.random.seed(seed)
    
    # Generate realistic RNA-seq counts using negative binomial
    base_expr = np.random.gamma(shape=2, scale=100, size=n_genes)
    
    counts = np.zeros((n_genes, n_samples))
    for i in range(n_samples):
        size_param = 10
        counts[:, i] = np.random.negative_binomial(
            size_param, 
            size_param / (size_param + base_expr)
        )
    
    # Gene names (mix of protein-coding, lncRNA, and pseudogenes)
    gene_names = []
    for i in range(n_genes):
        if i < 15000:
            gene_names.append(f'ENSG{i+1:011d}')
        elif i < 18000:
            gene_names.append(f'ENSG{i+1:011d}_lncRNA')
        else:
            gene_names.append(f'ENSG{i+1:011d}_pseudogene')
    
    sample_names = ['Control_1', 'Control_2', 'Control_3', 
                    'Treatment_1', 'Treatment_2', 'Treatment_3']
    
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=gene_names,
        columns=sample_names[:n_samples]
    )
    
    return counts_df


def generate_demo_tpm(counts_df, seed=42):
    """Generate demonstration TPM matrix from counts."""
    np.random.seed(seed)
    
    # Simulate gene lengths (1000-10000 bp)
    gene_lengths = np.random.uniform(1000, 10000, size=len(counts_df))
    
    # Calculate RPK (reads per kilobase)
    rpk = counts_df.div(gene_lengths / 1000, axis=0)
    
    # Calculate TPM (transcripts per million)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    
    return tpm


def generate_demo_quant_stats():
    """Generate demonstration quantification statistics."""
    return {
        'samples': {
            'Control_1': {
                'num_processed': 25000000,
                'num_mapped': 23500000,
                'percent_mapped': 94.0,
                'num_unique': 22000000,
                'num_multimapped': 1500000
            },
            'Control_2': {
                'num_processed': 28000000,
                'num_mapped': 26320000,
                'percent_mapped': 94.0,
                'num_unique': 24640000,
                'num_multimapped': 1680000
            },
            'Control_3': {
                'num_processed': 26000000,
                'num_mapped': 24440000,
                'percent_mapped': 94.0,
                'num_unique': 22880000,
                'num_multimapped': 1560000
            },
            'Treatment_1': {
                'num_processed': 27000000,
                'num_mapped': 25380000,
                'percent_mapped': 94.0,
                'num_unique': 23760000,
                'num_multimapped': 1620000
            },
            'Treatment_2': {
                'num_processed': 29000000,
                'num_mapped': 27260000,
                'percent_mapped': 94.0,
                'num_unique': 25520000,
                'num_multimapped': 1740000
            },
            'Treatment_3': {
                'num_processed': 26500000,
                'num_mapped': 24910000,
                'percent_mapped': 94.0,
                'num_unique': 23320000,
                'num_multimapped': 1590000
            }
        },
        'overall': {
            'total_processed': 161500000,
            'total_mapped': 151810000,
            'avg_percent_mapped': 94.0
        }
    }


def display_sample_sheet(samples, read_type):
    """Display sample sheet in a formatted table."""
    print("\n📋 Sample Sheet:")
    print("  " + "="*70)
    print(f"  {'Sample ID':<20} {'Read Type':<15} {'FASTQ Files'}")
    print("  " + "-"*70)
    
    for sample in samples[:5]:  # Show first 5
        sample_id = sample['sample_id']
        if read_type == 'paired-end':
            fastq_info = f"{sample['fastq_r1']}, {sample['fastq_r2']}"
        else:
            fastq_info = sample.get('fastq_r1', sample.get('fastq', 'N/A'))
        print(f"  {sample_id:<20} {read_type:<15} {fastq_info}")
    
    if len(samples) > 5:
        print(f"  ... and {len(samples) - 5} more samples")
    
    print("  " + "="*70)


def display_count_matrix_summary(counts_df):
    """Display count matrix summary statistics."""
    print("\n📊 Count Matrix Summary:")
    print("  " + "="*70)
    print(f"  Dimensions: {counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples")
    print(f"  Total counts: {counts_df.sum().sum():,.0f}")
    print(f"\n  Library Sizes (total counts per sample):")
    
    lib_sizes = counts_df.sum(axis=0)
    for sample in lib_sizes.index:
        print(f"    {sample:<20} {lib_sizes[sample]:>12,.0f}")
    
    print(f"\n  Library size range: {lib_sizes.min():,.0f} - {lib_sizes.max():,.0f}")
    print(f"  Library size CV: {(lib_sizes.std() / lib_sizes.mean()):.2%}")
    print("  " + "="*70)


def display_quality_check(counts_df):
    """Display quality check results."""
    print("\n✅ Quality Checks:")
    print("  " + "="*70)
    
    # Check 1: Library sizes
    lib_sizes = counts_df.sum(axis=0)
    lib_cv = lib_sizes.std() / lib_sizes.mean()
    
    if lib_cv < 0.15:
        status = "✓ PASS"
    elif lib_cv < 0.30:
        status = "⚠ WARNING"
    else:
        status = "✗ FAIL"
    print(f"  Library Size Consistency: {status}")
    print(f"    CV = {lib_cv:.2%} (target: <15%)")
    
    # Check 2: Gene detection
    genes_detected_per_sample = (counts_df > 0).sum(axis=0)
    avg_genes = genes_detected_per_sample.mean()
    
    if avg_genes > 15000:
        status = "✓ PASS"
    elif avg_genes > 10000:
        status = "⚠ WARNING"
    else:
        status = "✗ FAIL"
    print(f"\n  Gene Detection: {status}")
    print(f"    Average genes detected: {avg_genes:,.0f} (target: >15,000)")
    
    # Check 3: Zero-count genes
    zero_genes = (counts_df == 0).all(axis=1).sum()
    zero_pct = (zero_genes / len(counts_df)) * 100
    
    if zero_pct < 5:
        status = "✓ PASS"
    elif zero_pct < 15:
        status = "⚠ WARNING"
    else:
        status = "✗ FAIL"
    print(f"\n  Zero-Count Genes: {status}")
    print(f"    {zero_genes:,} genes ({zero_pct:.1f}%) with zero counts (target: <5%)")
    
    print("  " + "="*70)


def run_quick_count(fastq_dir=None, sample_sheet=None, index=None, 
                   output=None, method='salmon', threads=8, demo=False):
    """
    Run quick quantification with validation.
    
    NEW: Added input validation with clear error messages.
    """
    
    # Demo mode
    if demo:
        print("\n🎬 Running in DEMO mode (simulated data)")
        print("  " + "="*70)
        
        # Generate demo data
        fastq_info = generate_demo_fastq_structure()
        display_sample_sheet(fastq_info['samples'], fastq_info['read_type'])
        
        # Generate demo counts
        counts_df = generate_demo_counts(n_samples=6, n_genes=20000)
        tpm_df = generate_demo_tpm(counts_df)
        
        # Output directory with validation
        try:
            output_dir = validate_directory_path(output or DEFAULT_OUTPUT_DIR, create_if_missing=True)
        except Exception as e:
            print(f"❌ Error creating output directory: {e}")
            sys.exit(1)
        
        # Save demo results
        counts_df.to_csv(output_dir / OUTPUT_COUNTS_FILE)
        tpm_df.to_csv(output_dir / OUTPUT_TPM_FILE)
        
        # Sample info
        sample_info = pd.DataFrame({
            'sample_id': counts_df.columns,
            'condition': ['Control'] * 3 + ['Treatment'] * 3,
            'batch': ['Batch1'] * 6
        })
        sample_info.to_csv(output_dir / OUTPUT_SAMPLE_INFO, index=False)
        
        display_count_matrix_summary(counts_df)
        display_quality_check(counts_df)
        
        results = {
            'timestamp': datetime.now().isoformat(),
            'raptor_version': '2.2.0',
            'module': 'M1',
            'stage': 1,
            'mode': 'demo',
            'method': method,
            'data_info': {
                'n_samples': 6,
                'n_genes': 20000,
                'read_type': 'paired-end'
            },
            'output_files': {
                'counts': str(output_dir / OUTPUT_COUNTS_FILE),
                'tpm': str(output_dir / OUTPUT_TPM_FILE),
                'sample_info': str(output_dir / OUTPUT_SAMPLE_INFO)
            }
        }
        
        return results, output_dir
    
    # Real data mode - REQUIRES RAPTOR
    if not RAPTOR_AVAILABLE or SampleSheet is None:
        print("❌ ERROR: RAPTOR pipelines not available for real data")
        print("   Install RAPTOR with: pip install -e .")
        print("   Or use --demo flag for demonstration")
        sys.exit(1)
    
    # Validate inputs - NEW VALIDATION
    try:
        # Validate sample sheet
        if sample_sheet:
            sample_sheet_path = validate_file_path(sample_sheet, must_exist=True)
            print(f"✓ Sample sheet validated: {sample_sheet_path}")
        elif fastq_dir:
            fastq_dir_path = validate_directory_path(fastq_dir, must_exist=True)
            print(f"✓ FASTQ directory validated: {fastq_dir_path}")
        else:
            raise ValidationError("Either --sample-sheet or --fastq-dir required")
        
        # Validate index
        index_path = validate_file_path(index, must_exist=True)
        print(f"✓ Index validated: {index_path}")
        
        # Validate threads
        validate_positive_integer(threads, 'threads')
        if threads > 64:
            print(f"⚠️  Warning: threads={threads} exceeds recommended maximum (64)")
        print(f"✓ Threads validated: {threads}")
        
        # Validate/create output directory
        output_dir = validate_directory_path(output or DEFAULT_OUTPUT_DIR, create_if_missing=True)
        print(f"✓ Output directory validated: {output_dir}")
        
    except FileNotFoundError as e:
        print(f"❌ File not found: {e}")
        print(f"   Current directory: {Path.cwd()}")
        sys.exit(1)
    except ValidationError as e:
        print(f"❌ Validation error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error during validation: {e}")
        sys.exit(1)
    
    # Continue with real quantification...
    print(f"\n🚀 Running {method.upper()} quantification...")
    print(f"   This would run the actual pipeline in production")
    print(f"   Use --demo to see simulated results")
    
    sys.exit(0)


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Quick Quantification (Module 1) - WITH VALIDATION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (no data required)
  python 01_quick_count.py --demo
  
  # Demo with Kallisto
  python 01_quick_count.py --demo --method kallisto
  
  # Real data with Salmon
  python 01_quick_count.py --sample-sheet samples.csv --index salmon_index/
  
  # Real data with Kallisto
  python 01_quick_count.py --sample-sheet samples.csv --index transcripts.idx --method kallisto
  
  # Auto-create sample sheet from FASTQ directory
  python 01_quick_count.py --fastq-dir data/fastq/ --index salmon_index/

CLI Equivalent:
  raptor quick-count -m salmon -s samples.csv -i salmon_index/

Workflow (Stage 1: Fast Profiling):
  Module 1: quick-count (this script) → quick_gene_counts.csv
  Module 2: raptor qc --counts results/quick_counts/quick_gene_counts.csv
  Module 3: raptor profile --counts results/quick_counts/quick_gene_counts.csv
  Module 4: raptor recommend

Output Location:
  results/quick_counts/
    ├── quick_gene_counts.csv   (gene-level count matrix)
    ├── quick_tpm.csv           (TPM normalized matrix)
    └── sample_info.csv         (sample metadata)

Required Files:
  - Sample sheet CSV (or FASTQ directory)
  - Salmon index directory OR Kallisto index file
  
Sample Sheet Format:
  sample_id,condition,batch,fastq_r1,fastq_r2
  Control_1,Control,Batch1,/path/to/Control_1_R1.fastq.gz,/path/to/Control_1_R2.fastq.gz
  ...
        """
    )
    
    # Input options
    parser.add_argument('--sample-sheet', '-s', help='Sample sheet CSV file')
    parser.add_argument('--fastq-dir', help='Directory containing FASTQ files (auto-detect samples)')
    parser.add_argument('--index', '-i', help='Salmon index directory or Kallisto index file')
    
    # Method selection
    parser.add_argument('--method', '-m', choices=['salmon', 'kallisto'], default='salmon',
                        help='Quantification method (default: salmon)')
    
    # Output options
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    
    # Performance options
    parser.add_argument('--threads', '-t', type=int, default=8,
                        help='Number of threads (default: 8)')
    
    # Demo mode
    parser.add_argument('--demo', action='store_true',
                        help='Run in demo mode with simulated data')
    
    # Show workflow
    parser.add_argument('--show-workflow', action='store_true',
                        help='Show RAPTOR workflow diagram')
    
    args = parser.parse_args()
    
    print_banner()
    
    if args.show_workflow:
        print_workflow()
        sys.exit(0)
    
    # Validate inputs
    if not args.demo:
        if not args.sample_sheet and not args.fastq_dir:
            print("❌ ERROR: Either --sample-sheet, --fastq-dir, or --demo is required")
            parser.print_help()
            sys.exit(1)
        
        if not args.index:
            print("❌ ERROR: --index is required for real data")
            print("   Or use --demo flag for demonstration")
            sys.exit(1)
    
    # Run quantification
    results, output_dir = run_quick_count(
        fastq_dir=args.fastq_dir,
        sample_sheet=args.sample_sheet,
        index=args.index,
        output=args.output,
        method=args.method,
        threads=args.threads,
        demo=args.demo
    )
    
    # Save results
    results_file = output_dir / 'quick_count_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    # Final summary
    print("\n" + "="*70)
    print("  ✅ MODULE 1 (QUANTIFY) COMPLETE!")
    print("="*70)
    
    print(f"\n  📂 Output Directory: {output_dir}")
    print(f"\n  📊 Output Files:")
    print(f"     • {OUTPUT_COUNTS_FILE:<25} - Count matrix ({results['data_info']['n_genes']:,} genes × {results['data_info']['n_samples']} samples)")
    print(f"     • {OUTPUT_TPM_FILE:<25} - TPM normalized matrix")
    print(f"     • {OUTPUT_SAMPLE_INFO:<25} - Sample metadata")
    print(f"     • quick_count_results.json  - Quantification stats")
    
    print(f"\n  🔜 Next Steps (Continue RAPTOR Workflow):")
    print(f"\n     Module 2 - Quality Assessment:")
    print(f"     python 02_quality_assessment_UPDATED.py --counts {output_dir}/{OUTPUT_COUNTS_FILE}")
    print(f"     ")
    print(f"     Or use the CLI:")
    print(f"     raptor qc --counts {output_dir}/{OUTPUT_COUNTS_FILE}")
    
    print(f"\n     Module 3 - Data Profiling:")
    print(f"     raptor profile --counts {output_dir}/{OUTPUT_COUNTS_FILE}")
    
    print(f"\n     Module 4 - Pipeline Recommendation:")
    print(f"     raptor recommend")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
