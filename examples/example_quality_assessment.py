#!/usr/bin/env python3

"""
Advanced Data Quality Assessment Examples - RAPTOR v2.1.0

Demonstrates comprehensive quality assessment including:
- Basic quality scoring
- Batch effect detection
- Outlier identification  
- Poor quality data detection
- Multi-dataset comparison

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.1.0
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import RAPTOR quality assessment modules
try:
    from raptor.data_quality_assessment import DataQualityAssessor, quick_quality_check
    QA_AVAILABLE = True
except ImportError as e:
    print(f"Error: RAPTOR quality assessment module not found")
    print(f"Details: {e}")
    print("\nPlease install RAPTOR:")
    print("  pip install raptor-rnaseq")
    sys.exit(1)


def print_header(text):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70 + "\n")


def print_success(text):
    """Print success message."""
    print(f"✓ {text}")


def print_info(text):
    """Print info message."""
    print(f"→ {text}")


def example_1_basic_usage():
    """Example 1: Basic quality assessment with visualization."""
    print_header("EXAMPLE 1: Basic Quality Assessment")
    
    print_info("Generating sample RNA-seq count data...")
    
    # Generate sample data with realistic distribution
    np.random.seed(42)
    n_genes = 1000
    n_samples = 8
    
    # Simulate counts with log-normal distribution (realistic for RNA-seq)
    means = np.random.lognormal(mean=6, sigma=2, size=n_genes)
    counts = np.zeros((n_genes, n_samples))
    
    for i in range(n_genes):
        # Use negative binomial for overdispersion
        counts[i, :] = np.random.negative_binomial(n=10, p=10/(10+means[i]), size=n_samples)
    
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"Sample{i+1}" for i in range(n_samples)]
    )
    
    print_success(f"Generated: {n_genes} genes × {n_samples} samples")
    print(f"\nLibrary sizes: {counts_df.sum(axis=0).min():,.0f} - "
          f"{counts_df.sum(axis=0).max():,.0f}")
    print(f"Zero proportion: {(counts_df == 0).sum().sum() / counts_df.size:.1%}")
    
    # Run quick quality check
    print_info("\nRunning quality assessment...")
    
    report = quick_quality_check(
        counts_df, 
        plot=True, 
        output_file='example1_quality.png'
    )
    
    # Display results
    print("\n" + "─" * 70)
    print("QUALITY ASSESSMENT RESULTS")
    print("─" * 70)
    print(f"Overall Score:  {report['overall']['score']:.1f}/100")
    print(f"Status:         {report['overall']['status'].upper()}")
    print(f"Recommendation: {report['overall']['recommendation']}")
    
    # Component scores
    print("\nComponent Scores:")
    for comp_name, comp_data in report['components'].items():
        status_icon = "✓" if comp_data['score'] >= 70 else "⚠" if comp_data['score'] >= 50 else "✗"
        print(f"  {status_icon} {comp_name.replace('_', ' ').title():20s}: "
              f"{comp_data['score']:5.1f}/100")
    
    # Flags
    if any(comp['flags'] for comp in report['components'].values()):
        print("\nWarnings:")
        for comp_name, comp_data in report['components'].items():
            if comp_data['flags']:
                for flag in comp_data['flags']:
                    print(f"  ⚠ {flag}")
    
    print_success("\nVisualization saved: example1_quality.png")


def example_2_with_metadata():
    """Example 2: Quality assessment with batch effect detection."""
    print_header("EXAMPLE 2: Batch Effect Detection")
    
    print_info("Generating data with batch effect...")
    
    # Generate data with systematic batch effect
    np.random.seed(42)
    n_genes = 1000
    batch1_samples = 6
    batch2_samples = 6
    n_samples = batch1_samples + batch2_samples
    
    # Batch 1: Normal distribution
    means1 = np.random.lognormal(mean=6, sigma=2, size=n_genes)
    batch1 = np.zeros((n_genes, batch1_samples))
    for i in range(n_genes):
        batch1[i, :] = np.random.negative_binomial(
            n=10, p=10/(10+means1[i]), size=batch1_samples
        )
    
    # Batch 2: Shifted distribution (batch effect)
    means2 = means1 * 1.3  # Systematic shift
    batch2 = np.zeros((n_genes, batch2_samples))
    for i in range(n_genes):
        batch2[i, :] = np.random.negative_binomial(
            n=8, p=8/(8+means2[i]), size=batch2_samples
        )
    
    # Combine batches
    counts = np.hstack([batch1, batch2])
    
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"Sample{i+1}" for i in range(n_samples)]
    )
    
    # Create metadata with batch information
    metadata = pd.DataFrame({
        'sample': counts_df.columns,
        'batch': ['Batch1'] * batch1_samples + ['Batch2'] * batch2_samples,
        'condition': (['Control', 'Treatment'] * 3) + (['Control', 'Treatment'] * 3)
    })
    
    print_success("Generated data with batch structure:")
    print(f"  • Batch 1: {batch1_samples} samples")
    print(f"  • Batch 2: {batch2_samples} samples (with 30% systematic shift)")
    
    # Assess quality with batch detection
    print_info("\nRunning quality assessment with batch detection...")
    
    assessor = DataQualityAssessor(counts_df, metadata)
    report = assessor.assess_quality()
    
    # Display results
    print("\n" + "─" * 70)
    print("QUALITY ASSESSMENT RESULTS")
    print("─" * 70)
    
    print(f"\nOverall Score:  {report['overall']['score']:.1f}/100")
    print(f"Status:         {report['overall']['status'].upper()}")
    print(f"Recommendation: {report['overall']['recommendation']}")
    
    # Batch effect details
    batch_info = report['components'].get('batch_effects', {})
    print(f"\n🔍 Batch Effect Detection:")
    print(f"   Detected:       {'YES' if batch_info.get('batch_detected') else 'NO'}")
    
    if batch_info.get('batch_detected'):
        print(f"   Variable:       {batch_info.get('batch_variable', 'N/A')}")
        print(f"   Strength:       {batch_info.get('batch_strength', 0):.2f}")
        print(f"   Recommendation: {batch_info.get('recommendation', 'N/A')}")
        
        if batch_info.get('flags'):
            print(f"\n   Warnings:")
            for flag in batch_info['flags']:
                print(f"     • {flag}")
    
    # Component scores
    print("\nComponent Scores:")
    for comp_name, comp_data in report['components'].items():
        status_icon = "✓" if comp_data['score'] >= 70 else "⚠" if comp_data['score'] >= 50 else "✗"
        print(f"  {status_icon} {comp_name.replace('_', ' ').title():20s}: "
              f"{comp_data['score']:5.1f}/100")
    
    # Generate visualization
    print_info("\nGenerating comprehensive visualization...")
    assessor.plot_quality_report('example2_quality_with_batch.png')
    
    print_success("Visualization saved: example2_quality_with_batch.png")


def example_3_poor_quality_data():
    """Example 3: Detecting poor quality data with multiple issues."""
    print_header("EXAMPLE 3: Poor Quality Data Detection")
    
    print_info("Generating problematic dataset with multiple quality issues...")
    
    np.random.seed(42)
    n_genes = 1000
    n_samples = 10
    
    # Simulate poor quality data
    means = np.random.lognormal(mean=5, sigma=2, size=n_genes)
    counts = np.zeros((n_genes, n_samples))
    
    for i in range(n_genes):
        counts[i, :] = np.random.negative_binomial(
            n=5, p=5/(5+means[i]), size=n_samples
        )
    
    # Add problems:
    # 1. Very uneven library sizes
    lib_size_multipliers = np.array([0.3, 0.5, 0.8, 1.0, 1.0, 1.0, 1.2, 1.5, 2.0, 3.5])
    counts = counts * lib_size_multipliers
    
    # 2. High zero inflation
    zero_mask = np.random.random((n_genes, n_samples)) < 0.65
    counts[zero_mask] = 0
    
    # 3. Add outlier sample (last sample)
    counts[:, 9] = counts[:, 9] * 4.5
    
    # 4. Add some extreme outlier genes
    outlier_genes = np.random.choice(n_genes, size=20, replace=False)
    counts[outlier_genes, :] = counts[outlier_genes, :] * 10
    
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"Sample{i+1}" for i in range(n_samples)]
    )
    
    print_success("Generated problematic data with:")
    print("  • Highly uneven library sizes (3.5× range)")
    print("  • 65% zero inflation")
    print("  • One outlier sample (Sample10)")
    print("  • Extreme outlier genes")
    
    lib_sizes = counts_df.sum(axis=0)
    print(f"\nLibrary size range: {lib_sizes.min():,.0f} - {lib_sizes.max():,.0f}")
    print(f"Library size CV: {lib_sizes.std() / lib_sizes.mean():.2f}")
    print(f"Zero proportion: {(counts_df == 0).sum().sum() / counts_df.size:.1%}")
    
    # Assess quality
    print_info("\nRunning quality assessment...")
    
    report = quick_quality_check(
        counts_df, 
        plot=True, 
        output_file='example3_poor_quality.png'
    )
    
    # Display results
    print("\n" + "─" * 70)
    print("QUALITY ISSUES DETECTED")
    print("─" * 70)
    
    print(f"\nOverall Score:  {report['overall']['score']:.1f}/100")
    print(f"Status:         {report['overall']['status'].upper()}")
    
    # List all issues
    print("\n⚠️  Detected Issues:")
    issue_count = 0
    
    for comp_name, comp_data in report['components'].items():
        if comp_data['flags']:
            print(f"\n   {comp_name.replace('_', ' ').title()}:")
            for flag in comp_data['flags']:
                print(f"     • {flag}")
                issue_count += 1
    
    if issue_count == 0:
        print("   No critical issues detected.")
    else:
        print(f"\n   Total issues found: {issue_count}")
    
    # Recommendations
    if report['overall']['recommendation']:
        print(f"\n💡 Recommendation:")
        print(f"   {report['overall']['recommendation']}")
    
    print_success("\nVisualization saved: example3_poor_quality.png")


def example_4_compare_datasets():
    """Example 4: Compare quality across multiple datasets."""
    print_header("EXAMPLE 4: Multi-Dataset Quality Comparison")
    
    print_info("Generating three datasets with different quality levels...")
    
    np.random.seed(42)
    n_genes = 1000
    
    # Define dataset parameters
    datasets = {
        'High Quality': {
            'n_samples': 12,
            'lib_size_cv': 0.1,
            'zero_pct': 0.15,
            'mean_depth': 8.0
        },
        'Medium Quality': {
            'n_samples': 8,
            'lib_size_cv': 0.35,
            'zero_pct': 0.40,
            'mean_depth': 6.0
        },
        'Low Quality': {
            'n_samples': 6,
            'lib_size_cv': 0.70,
            'zero_pct': 0.65,
            'mean_depth': 4.0
        }
    }
    
    results = {}
    
    print()
    for name, params in datasets.items():
        print(f"Generating '{name}' dataset...")
        
        # Generate data
        means = np.random.lognormal(mean=params['mean_depth'], sigma=1.5, size=n_genes)
        counts = np.zeros((n_genes, params['n_samples']))
        
        for i in range(n_genes):
            counts[i, :] = np.random.negative_binomial(
                n=10, p=10/(10+means[i]), size=params['n_samples']
            )
        
        # Apply library size variation
        lib_multipliers = np.random.lognormal(0, params['lib_size_cv'], params['n_samples'])
        counts = counts * lib_multipliers
        
        # Apply zero inflation
        zero_mask = np.random.random((n_genes, params['n_samples'])) < params['zero_pct']
        counts[zero_mask] = 0
        
        counts_df = pd.DataFrame(
            counts.astype(int),
            index=[f"GENE{i:05d}" for i in range(n_genes)],
            columns=[f"Sample{i+1}" for i in range(params['n_samples'])]
        )
        
        # Assess (without plotting for comparison)
        assessor = DataQualityAssessor(counts_df)
        report = assessor.assess_quality()
        
        results[name] = report
    
    print_success("All datasets generated and assessed")
    
    # Display comparison
    print("\n" + "─" * 70)
    print("QUALITY COMPARISON")
    print("─" * 70)
    
    for name, report in results.items():
        score = report['overall']['score']
        status = report['overall']['status']
        
        # Status icon
        if status == 'good':
            icon = "✓"
        elif status == 'acceptable':
            icon = "⚠"
        else:
            icon = "✗"
        
        print(f"\n{icon} {name}:")
        print(f"   Overall Score: {score:.1f}/100 ({status.upper()})")
        
        # Component scores
        print(f"   Component Scores:")
        for comp_name, comp_data in report['components'].items():
            print(f"     • {comp_name.replace('_', ' ').title():20s}: "
                  f"{comp_data['score']:5.1f}/100")
    
    # Create comparison visualization
    print_info("\nCreating comparison visualization...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Overall scores
    names = list(results.keys())
    scores = [results[name]['overall']['score'] for name in names]
    colors = ['green' if s >= 80 else 'orange' if s >= 60 else 'red' for s in scores]
    
    bars = axes[0].bar(names, scores, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    axes[0].set_ylabel('Overall Score', fontsize=12)
    axes[0].set_title('Overall Quality Scores', fontsize=14, fontweight='bold')
    axes[0].set_ylim([0, 100])
    axes[0].axhline(80, color='green', linestyle='--', alpha=0.3, label='Good (≥80)')
    axes[0].axhline(60, color='orange', linestyle='--', alpha=0.3, label='Acceptable (≥60)')
    axes[0].axhline(50, color='red', linestyle='--', alpha=0.3, label='Poor (<60)')
    axes[0].legend(fontsize=9)
    axes[0].grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, score in zip(bars, scores):
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{score:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Plot 2: Component comparison
    components = ['library_quality', 'gene_detection', 'outlier_detection', 
                  'biological_signal']
    comp_labels = [c.replace('_', '\n').title() for c in components]
    
    x = np.arange(len(components))
    width = 0.25
    
    for i, name in enumerate(names):
        comp_scores = [results[name]['components'][c]['score'] for c in components]
        offset = (i - 1) * width
        axes[1].bar(x + offset, comp_scores, width, label=name, alpha=0.7, edgecolor='black')
    
    axes[1].set_xlabel('Component', fontsize=12)
    axes[1].set_ylabel('Score', fontsize=12)
    axes[1].set_title('Component Score Comparison', fontsize=14, fontweight='bold')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(comp_labels, fontsize=9)
    axes[1].legend(fontsize=10)
    axes[1].set_ylim([0, 100])
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example4_comparison.png', dpi=300, bbox_inches='tight')
    
    print_success("Comparison plot saved: example4_comparison.png")


def main():
    """Run all quality assessment examples."""
    
    print("""
╔═══════════════════════════════════════════════════════════════════╗
║     🦖 RAPTOR Advanced Data Quality Assessment Examples          ║
║        Version 2.1.0                                              ║
╚═══════════════════════════════════════════════════════════════════╝
    """)
    
    print("This script demonstrates RAPTOR's comprehensive quality assessment")
    print("capabilities including batch detection, outlier identification,")
    print("and multi-dataset comparison.\n")
    
    try:
        # Run all examples
        example_1_basic_usage()
        example_2_with_metadata()
        example_3_poor_quality_data()
        example_4_compare_datasets()
        
        # Success summary
        print_header("✅ ALL EXAMPLES COMPLETE!")
        
        print("Generated files:")
        print("  • example1_quality.png              - Basic quality assessment")
        print("  • example2_quality_with_batch.png   - Batch effect detection")
        print("  • example3_poor_quality.png         - Poor quality data")
        print("  • example4_comparison.png           - Multi-dataset comparison")
        
        print("\nKey takeaways:")
        print("  ✓ Quality assessment provides comprehensive scoring")
        print("  ✓ Batch effects can be automatically detected")
        print("  ✓ Multiple quality issues are flagged with recommendations")
        print("  ✓ Datasets can be compared for quality control")
        
        print("\nNext steps:")
        print("  • Apply quality assessment to your own data")
        print("  • Use metadata for batch effect detection")
        print("  • Integrate with RAPTOR profiling workflow")
        print("  • Address quality issues before analysis")
        
        print("\n" + "=" * 70)
        print("Thank you for using RAPTOR! 🦖")
        print("For support: ayehbolouki1988@gmail.com")
        print("=" * 70 + "\n")
        
    except KeyboardInterrupt:
        print("\n\nExamples interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
