#!/usr/bin/env python3

"""
Example: Advanced Data Quality Assessment

Demonstrates comprehensive quality assessment including:
- Batch effect detection
- Quality scoring
- Outlier detection
- Visualization

Author: Ayeh Bolouki
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from data_quality_assessment import DataQualityAssessor, quick_quality_check


def example_1_basic_usage():
    """Example 1: Basic quality assessment."""
    print("\n" + "="*70)
    print("EXAMPLE 1: Basic Quality Assessment")
    print("="*70 + "\n")
    
    # Generate sample data
    np.random.seed(42)
    n_genes = 1000
    n_samples = 8
    
    # Simulate counts with realistic distribution
    means = np.random.lognormal(mean=6, sigma=2, size=n_genes)
    counts = np.random.negative_binomial(n=10, p=0.1, size=(n_genes, n_samples))
    
    counts_df = pd.DataFrame(
        counts,
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"Sample{i+1}" for i in range(n_samples)]
    )
    
    print(f"Generated sample data: {n_genes} genes × {n_samples} samples\n")
    
    # Quick quality check
    report = quick_quality_check(counts_df, plot=True, output_file='example1_quality.png')
    
    print("\n✅ Assessment complete!")
    print(f"   Overall Score: {report['overall']['score']:.1f}/100")
    print(f"   Status: {report['overall']['status'].upper()}")
    print(f"   Report saved: example1_quality.png")


def example_2_with_metadata():
    """Example 2: Quality assessment with batch detection."""
    print("\n" + "="*70)
    print("EXAMPLE 2: Quality Assessment with Batch Detection")
    print("="*70 + "\n")
    
    # Generate data with batch effect
    np.random.seed(42)
    n_genes = 1000
    n_samples = 12
    
    # Create two batches
    batch1_samples = 6
    batch2_samples = 6
    
    # Batch 1: Normal distribution
    batch1 = np.random.negative_binomial(n=10, p=0.1, size=(n_genes, batch1_samples))
    
    # Batch 2: Shifted distribution (batch effect)
    batch2 = np.random.negative_binomial(n=8, p=0.15, size=(n_genes, batch2_samples))
    batch2 = batch2 * 1.5  # Systematic shift
    
    counts = np.hstack([batch1, batch2])
    
    counts_df = pd.DataFrame(
        counts,
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"Sample{i+1}" for i in range(n_samples)]
    )
    
    # Create metadata with batch information
    metadata = pd.DataFrame({
        'sample': counts_df.columns,
        'batch': ['Batch1'] * batch1_samples + ['Batch2'] * batch2_samples,
        'condition': ['Control', 'Treatment'] * 6
    })
    
    print("Generated data with batch effect:")
    print(f"  - {n_genes} genes × {n_samples} samples")
    print(f"  - Batch 1: {batch1_samples} samples")
    print(f"  - Batch 2: {batch2_samples} samples (with systematic shift)\n")
    
    # Assess quality with batch detection
    assessor = DataQualityAssessor(counts_df, metadata)
    report = assessor.assess_quality()
    
    # Print results
    print("\n" + "-"*70)
    print("QUALITY ASSESSMENT RESULTS")
    print("-"*70)
    
    print(f"\nOverall Score: {report['overall']['score']:.1f}/100")
    print(f"Status: {report['overall']['status'].upper()}")
    print(f"Recommendation: {report['overall']['recommendation']}")
    
    # Batch effect details
    batch_info = report['components']['batch_effects']
    print(f"\n🔍 Batch Effect Detection:")
    print(f"   Detected: {'YES' if batch_info['batch_detected'] else 'NO'}")
    if batch_info['batch_detected']:
        print(f"   Variable: {batch_info['batch_variable']}")
        print(f"   Strength: {batch_info['batch_strength']:.2f}")
        print(f"   Recommendation: {batch_info['recommendation']}")
    
    # Generate visualization
    assessor.plot_quality_report('example2_quality_with_batch.png')
    print(f"\n✅ Visualization saved: example2_quality_with_batch.png")


def example_3_poor_quality_data():
    """Example 3: Detecting poor quality data."""
    print("\n" + "="*70)
    print("EXAMPLE 3: Poor Quality Data Detection")
    print("="*70 + "\n")
    
    np.random.seed(42)
    n_genes = 1000
    n_samples = 10
    
    # Simulate poor quality data
    counts = np.random.negative_binomial(n=5, p=0.3, size=(n_genes, n_samples))
    
    # Add problems:
    # 1. Very uneven library sizes
    lib_size_multipliers = np.array([0.3, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5, 2.0, 3.0])
    counts = counts * lib_size_multipliers
    
    # 2. High zero inflation
    zero_mask = np.random.random((n_genes, n_samples)) < 0.6
    counts[zero_mask] = 0
    
    # 3. Add outlier sample
    counts[:, 9] = counts[:, 9] * 5  # Make last sample an outlier
    
    counts_df = pd.DataFrame(
        counts,
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"Sample{i+1}" for i in range(n_samples)]
    )
    
    print("Generated poor quality data with:")
    print("  - Highly uneven library sizes")
    print("  - 60% zero inflation")
    print("  - One outlier sample\n")
    
    # Assess quality
    report = quick_quality_check(counts_df, plot=True, output_file='example3_poor_quality.png')
    
    print(f"\n⚠️  Quality Issues Detected:")
    
    for component_name, component in report['components'].items():
        if component['flags']:
            print(f"\n   {component_name.replace('_', ' ').title()}:")
            for flag in component['flags']:
                print(f"     • {flag}")
    
    print(f"\n✅ Report saved: example3_poor_quality.png")


def example_4_compare_datasets():
    """Example 4: Compare quality across multiple datasets."""
    print("\n" + "="*70)
    print("EXAMPLE 4: Compare Multiple Datasets")
    print("="*70 + "\n")
    
    np.random.seed(42)
    n_genes = 1000
    
    datasets = {
        'High Quality': {
            'n_samples': 12,
            'lib_size_cv': 0.1,
            'zero_pct': 0.2
        },
        'Medium Quality': {
            'n_samples': 8,
            'lib_size_cv': 0.3,
            'zero_pct': 0.4
        },
        'Low Quality': {
            'n_samples': 6,
            'lib_size_cv': 0.6,
            'zero_pct': 0.7
        }
    }
    
    results = {}
    
    for name, params in datasets.items():
        # Generate data
        counts = np.random.negative_binomial(n=10, p=0.1, size=(n_genes, params['n_samples']))
        
        # Apply library size variation
        lib_multipliers = np.random.lognormal(0, params['lib_size_cv'], params['n_samples'])
        counts = counts * lib_multipliers
        
        # Apply zero inflation
        zero_mask = np.random.random((n_genes, params['n_samples'])) < params['zero_pct']
        counts[zero_mask] = 0
        
        counts_df = pd.DataFrame(
            counts,
            index=[f"GENE{i:05d}" for i in range(n_genes)],
            columns=[f"Sample{i+1}" for i in range(params['n_samples'])]
        )
        
        # Assess (no plot for comparison)
        assessor = DataQualityAssessor(counts_df)
        report = assessor.assess_quality()
        
        results[name] = report
        
        print(f"\n{name}:")
        print(f"  Overall Score: {report['overall']['score']:.1f}/100")
        print(f"  Status: {report['overall']['status']}")
        
        # Show component scores
        print(f"  Component Scores:")
        for comp_name, comp_data in report['components'].items():
            print(f"    - {comp_name.replace('_', ' ').title()}: {comp_data['score']:.1f}/100")
    
    # Create comparison plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Overall scores
    names = list(results.keys())
    scores = [results[name]['overall']['score'] for name in names]
    colors = ['green' if s >= 80 else 'orange' if s >= 60 else 'red' for s in scores]
    
    axes[0].bar(names, scores, color=colors, alpha=0.7, edgecolor='black')
    axes[0].set_ylabel('Overall Score', fontsize=12)
    axes[0].set_title('Overall Quality Scores', fontsize=14, fontweight='bold')
    axes[0].set_ylim([0, 100])
    axes[0].axhline(80, color='green', linestyle='--', alpha=0.3, label='Good')
    axes[0].axhline(60, color='orange', linestyle='--', alpha=0.3, label='Acceptable')
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)
    
    # Plot 2: Component comparison
    components = ['library_quality', 'gene_detection', 'outlier_detection', 
                  'biological_signal']
    comp_labels = [c.replace('_', '\n').title() for c in components]
    
    x = np.arange(len(components))
    width = 0.25
    
    for i, name in enumerate(names):
        comp_scores = [results[name]['components'][c]['score'] for c in components]
        axes[1].bar(x + i*width, comp_scores, width, label=name, alpha=0.7)
    
    axes[1].set_xlabel('Component', fontsize=12)
    axes[1].set_ylabel('Score', fontsize=12)
    axes[1].set_title('Component Score Comparison', fontsize=14, fontweight='bold')
    axes[1].set_xticks(x + width)
    axes[1].set_xticklabels(comp_labels, fontsize=9)
    axes[1].legend()
    axes[1].set_ylim([0, 100])
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example4_comparison.png', dpi=300, bbox_inches='tight')
    
    print(f"\n\n✅ Comparison plot saved: example4_comparison.png")


def main():
    """Run all examples."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║     🦖 RAPTOR Advanced Data Quality Assessment Examples     ║
    ╚══════════════════════════════════════════════════════════════╝
    """)
    
    try:
        # Example 1: Basic usage
        example_1_basic_usage()
        
        # Example 2: With batch detection
        example_2_with_metadata()
        
        # Example 3: Poor quality data
        example_3_poor_quality_data()
        
        # Example 4: Compare datasets
        example_4_compare_datasets()
        
        print("\n" + "="*70)
        print("✅ ALL EXAMPLES COMPLETE!")
        print("="*70)
        print("\nGenerated files:")
        print("  - example1_quality.png")
        print("  - example2_quality_with_batch.png")
        print("  - example3_poor_quality.png")
        print("  - example4_comparison.png")
        print("\n" + "="*70 + "\n")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
