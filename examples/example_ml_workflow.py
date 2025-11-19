#!/usr/bin/env python3

"""
ML Recommender Complete Example

Demonstrates the full workflow:
1. Generate synthetic training data
2. Train ML recommender
3. Evaluate performance
4. Make predictions on new data

Author: Ayeh Bolouki
Affiliation: University of Namur, Belgium
Email: ayehbolouki1988@gmail.com
"""

import sys
import argparse
from pathlib import Path
import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Import RAPTOR ML modules
try:
    from ml_recommender import MLPipelineRecommender, FeatureExtractor, train_recommender
    from synthetic_benchmarks import generate_training_data
    ML_AVAILABLE = True
except ImportError:
    print("Error: ML modules not found in current directory")
    print("Make sure ml_recommender.py and synthetic_benchmarks.py are present")
    sys.exit(1)


def print_header(text):
    """Print section header."""
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70 + "\n")


def step1_generate_data(n_datasets=200, output_dir='ml_training_data'):
    """Step 1: Generate synthetic training data."""
    print_header("STEP 1: Generate Synthetic Training Data")
    
    print(f"Generating {n_datasets} synthetic benchmark datasets...")
    print("This simulates diverse RNA-seq datasets with known optimal pipelines.\n")
    
    summary = generate_training_data(
        n_datasets=n_datasets,
        output_dir=output_dir,
        seed=42
    )
    
    print(f"\n✓ Generated {summary['n_datasets']} datasets")
    print(f"✓ Saved to: {summary['output_dir']}")
    
    return summary


def step2_train_model(benchmark_dir, model_type='random_forest', output_dir='models'):
    """Step 2: Train ML recommender."""
    print_header("STEP 2: Train ML Recommender")
    
    print(f"Training {model_type} model on benchmark data...")
    print(f"Reading from: {benchmark_dir}\n")
    
    recommender = MLPipelineRecommender(model_type=model_type)
    
    # Train
    results = recommender.train_from_benchmarks(
        benchmark_dir=benchmark_dir,
        performance_metric='f1_score'
    )
    
    # Display results
    print("\n--- Training Results ---")
    print(f"Model: {results['model_type']}")
    print(f"Training samples: {results['n_samples']}")
    print(f"Features: {results['n_features']}")
    print(f"Train accuracy: {results['train_score']:.3f}")
    print(f"Test accuracy: {results['test_score']:.3f}")
    print(f"Cross-validation: {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
    
    # Show top features
    if results['feature_importance']:
        print("\n--- Top 10 Most Important Features ---")
        for i, feat in enumerate(results['feature_importance'][:10], 1):
            print(f"{i:2d}. {feat['feature']:30s} {feat['importance']:.4f}")
    
    # Save model
    recommender.save_model(output_dir)
    print(f"\n✓ Model saved to: {output_dir}")
    
    return recommender, results


def step3_visualize_performance(results, output_dir='figures'):
    """Step 3: Visualize model performance."""
    print_header("STEP 3: Visualize Performance")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 1. Confusion Matrix
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    cm = np.array(results['confusion_matrix'])
    
    # Normalize confusion matrix
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    # Plot
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0], cbar=True)
    axes[0].set_title('Confusion Matrix (Counts)', fontsize=14, fontweight='bold')
    axes[0].set_xlabel('Predicted Pipeline', fontsize=12)
    axes[0].set_ylabel('True Pipeline', fontsize=12)
    
    sns.heatmap(cm_norm, annot=True, fmt='.2f', cmap='Greens', ax=axes[1], cbar=True)
    axes[1].set_title('Confusion Matrix (Normalized)', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Predicted Pipeline', fontsize=12)
    axes[1].set_ylabel('True Pipeline', fontsize=12)
    
    plt.tight_layout()
    cm_file = output_path / 'confusion_matrix.png'
    plt.savefig(cm_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved confusion matrix: {cm_file}")
    plt.close()
    
    # 2. Feature Importance
    if results['feature_importance']:
        feat_df = pd.DataFrame(results['feature_importance'][:15])
        
        fig, ax = plt.subplots(figsize=(10, 8))
        bars = ax.barh(feat_df['feature'], feat_df['importance'], color='steelblue')
        
        # Color top 5 differently
        for i in range(min(5, len(bars))):
            bars[i].set_color('darkred')
        
        ax.set_xlabel('Importance', fontsize=12)
        ax.set_title('Top 15 Feature Importances', fontsize=14, fontweight='bold')
        ax.invert_yaxis()
        
        plt.tight_layout()
        feat_file = output_path / 'feature_importance.png'
        plt.savefig(feat_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved feature importance: {feat_file}")
        plt.close()
    
    # 3. Classification Report
    report = results['classification_report']
    
    # Extract per-class metrics
    classes = [str(i) for i in range(1, 9)]
    metrics_data = []
    
    for cls in classes:
        if cls in report:
            metrics_data.append({
                'Pipeline': int(cls),
                'Precision': report[cls]['precision'],
                'Recall': report[cls]['recall'],
                'F1-Score': report[cls]['f1-score']
            })
    
    if metrics_data:
        metrics_df = pd.DataFrame(metrics_data)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        x = np.arange(len(metrics_df))
        width = 0.25
        
        ax.bar(x - width, metrics_df['Precision'], width, label='Precision', color='skyblue')
        ax.bar(x, metrics_df['Recall'], width, label='Recall', color='lightcoral')
        ax.bar(x + width, metrics_df['F1-Score'], width, label='F1-Score', color='lightgreen')
        
        ax.set_xlabel('Pipeline ID', fontsize=12)
        ax.set_ylabel('Score', fontsize=12)
        ax.set_title('Per-Pipeline Performance Metrics', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_df['Pipeline'])
        ax.legend()
        ax.set_ylim([0, 1.1])
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        metrics_file = output_path / 'pipeline_metrics.png'
        plt.savefig(metrics_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved pipeline metrics: {metrics_file}")
        plt.close()
    
    print(f"\n✓ All figures saved to: {output_dir}")


def step4_test_predictions(recommender, n_test_cases=5):
    """Step 4: Test predictions on new synthetic data."""
    print_header("STEP 4: Test Predictions on New Data")
    
    from synthetic_benchmarks import SyntheticBenchmarkGenerator
    
    print(f"Generating {n_test_cases} test cases...\n")
    
    # Generate test profiles
    generator = SyntheticBenchmarkGenerator(n_datasets=n_test_cases, seed=999)
    
    for i in range(n_test_cases):
        profile = generator._generate_profile()
        
        print(f"--- Test Case {i+1} ---")
        print(f"Samples: {profile['design']['n_samples']}")
        print(f"Genes: {profile['design']['n_genes']}")
        print(f"BCV: {profile['biological_variation']['bcv']:.3f}")
        print(f"Depth: {profile['sequencing']['depth_category']}")
        print(f"Zero%: {profile['count_distribution']['zero_pct']:.1f}%")
        
        # Get recommendation
        recommendation = recommender.recommend(profile, top_k=3)
        
        print(f"\n🦖 RECOMMENDATION:")
        print(f"   Pipeline {recommendation['pipeline_id']}: {recommendation['pipeline_name']}")
        print(f"   Confidence: {recommendation['confidence']:.1%}")
        print(f"   Reasons:")
        for reason in recommendation['reasons']:
            print(f"     • {reason}")
        
        if recommendation['alternatives']:
            print(f"\n   Alternatives:")
            for alt in recommendation['alternatives']:
                print(f"     • Pipeline {alt['pipeline_id']} ({alt['confidence']:.1%}): {alt['pipeline_name']}")
        
        print()


def step5_compare_models(benchmark_dir):
    """Step 5: Compare Random Forest vs Gradient Boosting."""
    print_header("STEP 5: Compare Model Types")
    
    models = ['random_forest', 'gradient_boosting']
    results_comparison = {}
    
    for model_type in models:
        print(f"\nTraining {model_type}...")
        recommender = MLPipelineRecommender(model_type=model_type)
        results = recommender.train_from_benchmarks(benchmark_dir)
        
        results_comparison[model_type] = {
            'test_accuracy': results['test_score'],
            'cv_mean': results['cv_mean'],
            'cv_std': results['cv_std']
        }
        
        print(f"  Test Accuracy: {results['test_score']:.3f}")
        print(f"  CV Score: {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
    
    # Visualize comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(models))
    width = 0.35
    
    test_scores = [results_comparison[m]['test_accuracy'] for m in models]
    cv_scores = [results_comparison[m]['cv_mean'] for m in models]
    cv_errs = [results_comparison[m]['cv_std'] for m in models]
    
    ax.bar(x - width/2, test_scores, width, label='Test Accuracy', color='steelblue')
    ax.bar(x + width/2, cv_scores, width, label='CV Accuracy', 
           yerr=cv_errs, capsize=5, color='coral')
    
    ax.set_xlabel('Model Type', fontsize=12)
    ax.set_ylabel('Accuracy', fontsize=12)
    ax.set_title('Model Comparison', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace('_', ' ').title() for m in models])
    ax.legend()
    ax.set_ylim([0, 1.1])
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/model_comparison.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved comparison: figures/model_comparison.png")
    plt.close()
    
    # Determine winner
    winner = max(results_comparison.items(), key=lambda x: x[1]['test_accuracy'])
    print(f"\n🏆 Best Model: {winner[0].replace('_', ' ').title()}")
    print(f"   Test Accuracy: {winner[1]['test_accuracy']:.3f}")


def main():
    """Main workflow."""
    parser = argparse.ArgumentParser(
        description='Complete ML Recommender Workflow Demo',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--n-datasets', type=int, default=200,
                       help='Number of synthetic datasets to generate (default: 200)')
    parser.add_argument('--model-type', choices=['random_forest', 'gradient_boosting'],
                       default='random_forest',
                       help='Model type to train (default: random_forest)')
    parser.add_argument('--skip-generation', action='store_true',
                       help='Skip data generation (use existing data)')
    parser.add_argument('--compare-models', action='store_true',
                       help='Compare Random Forest vs Gradient Boosting')
    parser.add_argument('--data-dir', default='ml_training_data',
                       help='Directory for training data')
    parser.add_argument('--model-dir', default='models',
                       help='Directory for saved models')
    
    args = parser.parse_args()
    
    print("\n" + "=" * 70)
    print("🦖 RAPTOR ML RECOMMENDER - COMPLETE WORKFLOW")
    print("=" * 70)
    
    try:
        # Step 1: Generate data (unless skipped)
        if not args.skip_generation:
            summary = step1_generate_data(
                n_datasets=args.n_datasets,
                output_dir=args.data_dir
            )
        else:
            print_header("STEP 1: Using Existing Data")
            print(f"Using data from: {args.data_dir}")
        
        # Step 2: Train model
        recommender, results = step2_train_model(
            benchmark_dir=args.data_dir,
            model_type=args.model_type,
            output_dir=args.model_dir
        )
        
        # Step 3: Visualize
        step3_visualize_performance(results, output_dir='figures')
        
        # Step 4: Test predictions
        step4_test_predictions(recommender, n_test_cases=5)
        
        # Step 5: Compare models (optional)
        if args.compare_models:
            step5_compare_models(args.data_dir)
        
        print_header("✓ WORKFLOW COMPLETE!")
        print("\nWhat you've accomplished:")
        print("  ✓ Generated synthetic training data")
        print("  ✓ Trained ML pipeline recommender")
        print("  ✓ Evaluated model performance")
        print("  ✓ Visualized results")
        print("  ✓ Tested predictions on new data")
        
        print("\nNext steps:")
        print("  • Use the trained model on real RNA-seq data")
        print("  • Integrate with RAPTOR's profiler module")
        print("  • Collect real benchmark data to improve the model")
        print("  • Deploy as part of RAPTOR CLI")
        
        print(f"\nModel saved at: {args.model_dir}")
        print(f"Figures saved at: figures/")
        print("\n" + "=" * 70 + "\n")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
