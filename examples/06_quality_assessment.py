#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Example Script: Advanced Data Quality Assessment

Demonstrates comprehensive data quality assessment with:
- Library quality scoring
- Gene detection analysis
- Outlier detection (multiple methods)
- Variance structure analysis
- Batch effect detection
- Biological signal assessment
- Overall quality scoring (0-100)

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
    from raptor.data_quality_assessment import DataQualityAssessor, quick_quality_check
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - Data Quality Assessment            ║
    ║                                                              ║
    ║   Comprehensive Quality Scoring & Batch Effect Detection    ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_data(n_samples=12, n_genes=5000, with_batch=True, seed=42):
    """Generate demonstration data with optional batch effects."""
    np.random.seed(seed)
    
    # Generate base expression
    base_expr = np.random.gamma(shape=2, scale=100, size=n_genes)
    
    # Generate counts
    counts = np.zeros((n_genes, n_samples))
    for i in range(n_samples):
        size_param = 10
        counts[:, i] = np.random.negative_binomial(
            size_param, 
            size_param / (size_param + base_expr)
        )
    
    # Add batch effect if requested
    if with_batch:
        batch_size = n_samples // 2
        batch_effect = np.random.normal(1.5, 0.3, n_genes)
        counts[:, :batch_size] = counts[:, :batch_size] * batch_effect[:, np.newaxis]
    
    # Create DataFrame
    gene_names = [f'Gene_{i+1:05d}' for i in range(n_genes)]
    sample_names = [f'Sample_{i+1}' for i in range(n_samples)]
    
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=gene_names,
        columns=sample_names
    )
    
    # Create metadata
    conditions = ['Control'] * (n_samples // 2) + ['Treatment'] * (n_samples // 2)
    batches = ['Batch1'] * (n_samples // 2) + ['Batch2'] * (n_samples // 2)
    
    metadata_df = pd.DataFrame({
        'sample': sample_names,
        'condition': conditions,
        'batch': batches
    })
    
    return counts_df, metadata_df


def generate_demo_quality_report():
    """Generate demonstration quality report."""
    return {
        'overall': {
            'score': 72.5,
            'status': 'good',
            'grade': 'B'
        },
        'library_quality': {
            'mean_size': 25000000.0,
            'cv': 0.22,
            'min_size': 18000000.0,
            'max_size': 32000000.0,
            'score': 85.0,
            'status': 'good',
            'flags': []
        },
        'gene_detection': {
            'mean_detection_rate': 0.72,
            'zero_inflation_pct': 45.2,
            'n_highly_expressed': 500,
            'n_medium_expressed': 1500,
            'n_low_expressed': 3000,
            'score': 75.0,
            'status': 'good',
            'flags': []
        },
        'outlier_detection': {
            'n_outliers': 1,
            'outlier_percentage': 8.3,
            'outlier_samples': ['Sample_7'],
            'mahalanobis_threshold': 3.5,
            'score': 70.0,
            'status': 'warning',
            'flags': ['1 potential outlier sample(s) detected']
        },
        'variance_structure': {
            'variance_explained_pc1': 0.35,
            'variance_explained_top5': [0.35, 0.18, 0.12, 0.08, 0.05],
            'total_variance_top5': 0.78,
            'score': 80.0,
            'status': 'good',
            'flags': []
        },
        'batch_effects': {
            'batch_detected': True,
            'batch_variable': 'batch',
            'batch_variance_explained': 0.25,
            'f_statistic': 15.2,
            'p_value': 0.0001,
            'pca_coordinates': [[1.2, 0.5], [-1.1, 0.3]],  # Simplified
            'score': 55.0,
            'status': 'warning',
            'flags': ['Significant batch effect detected (p < 0.001)']
        },
        'biological_signal': {
            'condition_separation': 0.65,
            'effect_size': 1.2,
            'signal_to_noise': 2.5,
            'score': 75.0,
            'status': 'good',
            'flags': []
        },
        'summary': """
RAPTOR Data Quality Assessment Summary
======================================
Overall Score: 72.5/100 (Good)

✓ Library Quality:     85.0  [████████░░]  Good
✓ Gene Detection:      75.0  [███████░░░]  Good
⚠ Outlier Detection:   70.0  [███████░░░]  Warning
✓ Variance Structure:  80.0  [████████░░]  Good
⚠ Batch Effects:       55.0  [█████░░░░░]  Warning
✓ Biological Signal:   75.0  [███████░░░]  Good

Recommendations:
1. Consider batch correction before differential expression analysis
2. Investigate Sample_7 for potential quality issues
3. Current quality is suitable for standard analysis pipelines
        """,
        'recommendations': [
            'Consider batch correction before differential expression analysis',
            'Investigate Sample_7 for potential quality issues',
            'Current quality is suitable for standard analysis pipelines'
        ]
    }


def display_quality_score(score, name, width=40):
    """Display a quality score with visual bar."""
    filled = int(score / 100 * width)
    bar = '█' * filled + '░' * (width - filled)
    
    if score >= 80:
        status = '✓'
        status_text = 'Good'
    elif score >= 60:
        status = '⚠'
        status_text = 'Warning'
    else:
        status = '✗'
        status_text = 'Poor'
    
    print(f"  {status} {name:20s} [{bar}] {score:5.1f}  {status_text}")


def display_quality_report(report):
    """Display quality report with formatting."""
    overall = report['overall']
    
    print("\n" + "="*70)
    print("  🦖 DATA QUALITY ASSESSMENT REPORT")
    print("="*70)
    
    # Overall score gauge
    score = overall['score']
    if score >= 80:
        grade = 'A'
        color_desc = 'Excellent'
    elif score >= 70:
        grade = 'B'
        color_desc = 'Good'
    elif score >= 60:
        grade = 'C'
        color_desc = 'Acceptable'
    elif score >= 50:
        grade = 'D'
        color_desc = 'Poor'
    else:
        grade = 'F'
        color_desc = 'Very Poor'
    
    print(f"\n  📊 OVERALL QUALITY SCORE: {score:.1f}/100 ({color_desc}) [Grade: {grade}]")
    
    # Component scores
    print("\n  Component Scores:")
    print("  " + "-"*66)
    
    display_quality_score(report['library_quality']['score'], 'Library Quality')
    display_quality_score(report['gene_detection']['score'], 'Gene Detection')
    display_quality_score(report['outlier_detection']['score'], 'Outlier Detection')
    display_quality_score(report['variance_structure']['score'], 'Variance Structure')
    display_quality_score(report['batch_effects']['score'], 'Batch Effects')
    display_quality_score(report['biological_signal']['score'], 'Biological Signal')
    
    # Flags and warnings
    all_flags = []
    for component in ['library_quality', 'gene_detection', 'outlier_detection', 
                      'variance_structure', 'batch_effects', 'biological_signal']:
        flags = report.get(component, {}).get('flags', [])
        all_flags.extend(flags)
    
    if all_flags:
        print("\n  ⚠️  Issues Detected:")
        for flag in all_flags:
            print(f"     • {flag}")
    
    # Recommendations
    if 'recommendations' in report:
        print("\n  📋 Recommendations:")
        for i, rec in enumerate(report['recommendations'], 1):
            print(f"     {i}. {rec}")
    
    # Detailed metrics
    print("\n  📈 Detailed Metrics:")
    print("  " + "-"*66)
    
    lib = report['library_quality']
    print(f"     Library Size: {lib['mean_size']/1e6:.1f}M reads (CV: {lib['cv']:.2f})")
    
    gene = report['gene_detection']
    print(f"     Zero Inflation: {gene['zero_inflation_pct']:.1f}%")
    print(f"     Gene Detection Rate: {gene['mean_detection_rate']:.1%}")
    
    outlier = report['outlier_detection']
    if outlier['n_outliers'] > 0:
        print(f"     Outliers: {outlier['n_outliers']} samples ({', '.join(outlier['outlier_samples'])})")
    else:
        print(f"     Outliers: None detected")
    
    batch = report['batch_effects']
    if batch['batch_detected']:
        print(f"     Batch Effect: Detected (variance explained: {batch['batch_variance_explained']:.1%})")
    else:
        print(f"     Batch Effect: Not detected")


def run_quality_assessment(counts_file, metadata_file=None, plot=False, 
                           no_batch=False, demo=False):
    """Run comprehensive quality assessment."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated data...")
        counts, metadata = generate_demo_data(with_batch=not no_batch)
        report = generate_demo_quality_report()
        
        print(f"\n📊 Demo data generated:")
        print(f"   Counts: {counts.shape[0]} genes × {counts.shape[1]} samples")
        print(f"   Metadata: {len(metadata)} samples with batch information")
        
    else:
        # Load data
        print(f"\n📂 Loading data from: {counts_file}")
        counts = pd.read_csv(counts_file, index_col=0)
        print(f"   Loaded: {counts.shape[0]} genes × {counts.shape[1]} samples")
        
        metadata = None
        if metadata_file:
            metadata = pd.read_csv(metadata_file)
            print(f"   Metadata: {len(metadata)} samples")
        
        # Run quality assessment
        print("\n🔍 Running comprehensive quality assessment...")
        
        assessor = DataQualityAssessor(counts, metadata)
        report = assessor.assess_quality()
        
        # Generate plot if requested
        if plot:
            plot_file = 'quality_report.png'
            print(f"\n📊 Generating quality visualization: {plot_file}")
            assessor.plot_quality_report(plot_file)
    
    # Display report
    display_quality_report(report)
    
    # Prepare output
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'data_info': {
            'n_genes': counts.shape[0] if not demo else 5000,
            'n_samples': counts.shape[1] if not demo else 12
        },
        'quality_report': report
    }
    
    return output


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Data Quality Assessment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic quality assessment
  python 06_quality_assessment.py --counts counts.csv
  
  # With metadata for batch effect detection
  python 06_quality_assessment.py --counts counts.csv --metadata samples.csv
  
  # Generate visualization
  python 06_quality_assessment.py --counts counts.csv --metadata samples.csv --plot
  
  # Skip batch effect detection
  python 06_quality_assessment.py --counts counts.csv --no-batch
  
  # Demo mode (no data required)
  python 06_quality_assessment.py --demo
        """
    )
    
    parser.add_argument('--counts', help='Count matrix CSV file (genes × samples)')
    parser.add_argument('--metadata', help='Sample metadata CSV file')
    parser.add_argument('--plot', action='store_true', help='Generate quality visualization')
    parser.add_argument('--no-batch', action='store_true', help='Skip batch effect detection')
    parser.add_argument('--output', '-o', help='Output JSON file')
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
    
    # Run assessment
    results = run_quality_assessment(
        counts_file=args.counts,
        metadata_file=args.metadata,
        plot=args.plot,
        no_batch=args.no_batch,
        demo=args.demo
    )
    
    # Save output
    output_file = args.output or 'quality_assessment_results.json'
    with open(output_file, 'w') as f:
        # Remove non-serializable items
        clean_results = json.loads(json.dumps(results, default=str))
        json.dump(clean_results, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_file}")
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
