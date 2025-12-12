#!/usr/bin/env python3

"""
ML Recommender Complete Example - RAPTOR v2.1.0

Demonstrates the full machine learning workflow:
1. Generate synthetic training data
2. Train ML recommender
3. Evaluate performance
4. Make predictions on new data
5. Compare different models

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.1.0
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
    from raptor.ml_recommender import MLPipelineRecommender, FeatureExtractor, train_recommender
    from raptor.synthetic_benchmarks import generate_training_data, SyntheticBenchmarkGenerator
    ML_AVAILABLE = True
except ImportError as e:
    print(f"Error: RAPTOR ML modules not found")
    print(f"Details: {e}")
    print("\nPlease install RAPTOR with ML support:")
    print("  pip install raptor-rnaseq[ml]")
    print("\nOr install dependencies manually:")
    print("  pip install scikit-learn>=1.0.0")
    sys.exit(1)


def print_header(text):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70 + "\n")


def print_success(text):
    """Print success message with checkmark."""
    print(f"✓ {text}")


def print_info(text):
    """Print info message."""
    print(f"→ {text}")


def print_error(text):
    """Print error message."""
    print(f"✗ ERROR: {text}", file=sys.stderr)


def step1_generate_data(n_datasets=200, output_dir='ml_training_data'):
    """
    Step 1: Generate synthetic training data.
    
    Creates simulated RNA-seq datasets with known optimal pipelines
    for training the machine learning recommender.
    """
    print_header("STEP 1: Generate Synthetic Training Data")
    
    print_info(f"Generating {n_datasets} synthetic benchmark datasets...")
    print("This simulates diverse RNA-seq datasets with known optimal pipelines.")
    print("Each dataset varies in:")
    print("  • Sample size (4-50 samples)")
    print("  • Number of genes (5,000-30,000)")
    print("  • Biological variation (BCV: 0.1-0.8)")
    print("  • Sequencing depth (low, medium, high)")
    print("  • Zero inflation and library size variation\n")
    
    try:
        summary = generate_training_data(
            n_datasets=n_datasets,
            output_dir=output_dir,
            seed=42
        )
        
        print_success(f"Generated {summary['n_datasets']} datasets")
        print_success(f"Saved to: {summary['output_dir']}")
        
        # Display distribution of optimal pipelines
        if 'pipeline_distribution' in summary:
            print("\nOptimal Pipeline Distribution:")
            for pipeline_id, count in sorted(summary['pipeline_distribution'].items()):
                pct = (count / summary['n_datasets']) * 100
                print(f"  Pipeline {pipeline_id}: {count:3d} datasets ({pct:5.1f}%)")
        
        return summary
        
    except Exception as e:
        print_error(f"Failed to generate training data: {e}")
        raise


def step2_train_model(benchmark_dir, model_type='random_forest', output_dir='models'):
    """
    Step 2: Train ML recommender.
    
    Trains a machine learning model to predict optimal pipelines
    based on dataset characteristics.
    """
    print_header("STEP 2: Train ML Recommender")
    
    print_info(f"Training {model_type.replace('_', ' ').title()} model...")
    print(f"Reading benchmark data from: {benchmark_dir}\n")
    
    try:
        recommender = MLPipelineRecommender(model_type=model_type)
        
        # Train the model
        results = recommender.train_from_benchmarks(
            benchmark_dir=benchmark_dir,
            performance_metric='f1_score'
        )
        
        # Display training results
        print("\n" + "-" * 70)
        print("TRAINING RESULTS")
        print("-" * 70)
        print(f"Model Type:           {results['model_type'].replace('_', ' ').title()}")
        print(f"Training Samples:     {results['n_samples']}")
        print(f"Features:             {results['n_features']}")
        print(f"Train Accuracy:       {results['train_score']:.3f}")
        print(f"Test Accuracy:        {results['test_score']:.3f}")
        print(f"Cross-Validation:     {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
        
        # Show top important features
        if results.get('feature_importance'):
            print("\n" + "-" * 70)
            print("TOP 10 MOST IMPORTANT FEATURES")
            print("-" * 70)
            for i, feat in enumerate(results['feature_importance'][:10], 1):
                bar_length = int(feat['importance'] * 40)
                bar = "█" * bar_length
                print(f"{i:2d}. {feat['feature']:30s} {bar} {feat['importance']:.4f}")
        
        # Classification report summary
        if 'classification_report' in results:
            report = results['classification_report']
            print("\n" + "-" * 70)
            print("CLASSIFICATION METRICS")
            print("-" * 70)
            print(f"Overall Accuracy:     {report.get('accuracy', 0):.3f}")
            print(f"Macro Avg F1:         {report.get('macro avg', {}).get('f1-score', 0):.3f}")
            print(f"Weighted Avg F1:      {report.get('weighted avg', {}).get('f1-score', 0):.3f}")
        
        # Save model
        model_path = recommender.save_model(output_dir)
        print_success(f"\nModel saved to: {output_dir}")
        
        return recommender, results
        
    except Exception as e:
        print_error(f"Training failed: {e}")
        raise


def step3_visualize_performance(results, output_dir='figures'):
    """
    Step 3: Visualize model performance.
    
    Creates comprehensive visualizations of model performance including
    confusion matrices, feature importance, and per-pipeline metrics.
    """
    print_header("STEP 3: Visualize Model Performance")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print_info("Generating performance visualizations...")
    
    try:
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 300
        
        # 1. Confusion Matrix
        if 'confusion_matrix' in results:
            print_info("Creating confusion matrix...")
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            
            cm = np.array(results['confusion_matrix'])
            
            # Normalize confusion matrix
            cm_norm = cm.astype('float') / (cm.sum(axis=1)[:, np.newaxis] + 1e-10)
            
            # Plot raw counts
            sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0], 
                       cbar_kws={'label': 'Count'}, square=True)
            axes[0].set_title('Confusion Matrix (Counts)', fontsize=14, fontweight='bold')
            axes[0].set_xlabel('Predicted Pipeline', fontsize=12)
            axes[0].set_ylabel('True Pipeline', fontsize=12)
            
            # Plot normalized
            sns.heatmap(cm_norm, annot=True, fmt='.2f', cmap='Greens', ax=axes[1], 
                       cbar_kws={'label': 'Proportion'}, square=True, vmin=0, vmax=1)
            axes[1].set_title('Confusion Matrix (Normalized)', fontsize=14, fontweight='bold')
            axes[1].set_xlabel('Predicted Pipeline', fontsize=12)
            axes[1].set_ylabel('True Pipeline', fontsize=12)
            
            plt.tight_layout()
            cm_file = output_path / 'confusion_matrix.png'
            plt.savefig(cm_file, dpi=300, bbox_inches='tight')
            print_success(f"Saved: {cm_file}")
            plt.close()
        
        # 2. Feature Importance
        if results.get('feature_importance'):
            print_info("Creating feature importance plot...")
            feat_df = pd.DataFrame(results['feature_importance'][:15])
            
            fig, ax = plt.subplots(figsize=(10, 8))
            bars = ax.barh(range(len(feat_df)), feat_df['importance'], color='steelblue')
            
            # Highlight top 5
            for i in range(min(5, len(bars))):
                bars[i].set_color('darkred')
            
            ax.set_yticks(range(len(feat_df)))
            ax.set_yticklabels(feat_df['feature'])
            ax.set_xlabel('Importance Score', fontsize=12)
            ax.set_title('Top 15 Feature Importances', fontsize=14, fontweight='bold')
            ax.invert_yaxis()
            
            # Add value labels
            for i, (idx, row) in enumerate(feat_df.iterrows()):
                ax.text(row['importance'], i, f" {row['importance']:.4f}", 
                       va='center', fontsize=9)
            
            plt.tight_layout()
            feat_file = output_path / 'feature_importance.png'
            plt.savefig(feat_file, dpi=300, bbox_inches='tight')
            print_success(f"Saved: {feat_file}")
            plt.close()
        
        # 3. Per-Pipeline Performance Metrics
        if 'classification_report' in results:
            print_info("Creating per-pipeline metrics...")
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
                        'F1-Score': report[cls]['f1-score'],
                        'Support': report[cls]['support']
                    })
            
            if metrics_data:
                metrics_df = pd.DataFrame(metrics_data)
                
                fig, ax = plt.subplots(figsize=(12, 6))
                x = np.arange(len(metrics_df))
                width = 0.25
                
                ax.bar(x - width, metrics_df['Precision'], width, 
                      label='Precision', color='skyblue', edgecolor='black')
                ax.bar(x, metrics_df['Recall'], width, 
                      label='Recall', color='lightcoral', edgecolor='black')
                ax.bar(x + width, metrics_df['F1-Score'], width, 
                      label='F1-Score', color='lightgreen', edgecolor='black')
                
                ax.set_xlabel('Pipeline ID', fontsize=12)
                ax.set_ylabel('Score', fontsize=12)
                ax.set_title('Per-Pipeline Performance Metrics', fontsize=14, fontweight='bold')
                ax.set_xticks(x)
                ax.set_xticklabels(metrics_df['Pipeline'])
                ax.legend(fontsize=11)
                ax.set_ylim([0, 1.1])
                ax.grid(axis='y', alpha=0.3)
                
                plt.tight_layout()
                metrics_file = output_path / 'pipeline_metrics.png'
                plt.savefig(metrics_file, dpi=300, bbox_inches='tight')
                print_success(f"Saved: {metrics_file}")
                plt.close()
        
        print_success(f"\nAll figures saved to: {output_dir}/")
        
    except Exception as e:
        print_error(f"Visualization failed: {e}")
        raise


def step4_test_predictions(recommender, n_test_cases=5):
    """
    Step 4: Test predictions on new synthetic data.
    
    Generates new test datasets and demonstrates the recommender
    making predictions with confidence scores.
    """
    print_header("STEP 4: Test Predictions on New Data")
    
    print_info(f"Generating {n_test_cases} test cases and making predictions...\n")
    
    try:
        generator = SyntheticBenchmarkGenerator(n_datasets=n_test_cases, seed=999)
        
        for i in range(n_test_cases):
            profile = generator._generate_profile()
            
            print(f"{'─' * 70}")
            print(f"TEST CASE {i+1}/{n_test_cases}")
            print(f"{'─' * 70}")
            
            # Display dataset characteristics
            print("Dataset Characteristics:")
            print(f"  • Samples:        {profile['design']['n_samples']}")
            print(f"  • Genes:          {profile['design']['n_genes']:,}")
            print(f"  • BCV:            {profile['biological_variation']['bcv']:.3f}")
            print(f"  • Depth:          {profile['sequencing']['depth_category']}")
            print(f"  • Zero%:          {profile['count_distribution']['zero_pct']:.1f}%")
            print(f"  • Library CV:     {profile['technical_variation']['library_size_cv']:.2f}")
            
            # Get recommendation
            recommendation = recommender.recommend(profile, top_k=3)
            
            # Display recommendation
            print(f"\n🦖 RECOMMENDATION:")
            print(f"   Pipeline {recommendation['pipeline_id']}: {recommendation['pipeline_name']}")
            print(f"   Confidence: {recommendation['confidence']:.1%}")
            
            if recommendation.get('reasons'):
                print(f"   Reasoning:")
                for reason in recommendation['reasons']:
                    print(f"     • {reason}")
            
            # Display alternatives
            if recommendation.get('alternatives'):
                print(f"\n   Alternative Pipelines:")
                for alt in recommendation['alternatives'][:2]:  # Show top 2
                    print(f"     • Pipeline {alt['pipeline_id']} ({alt['confidence']:.1%}): "
                          f"{alt['pipeline_name']}")
            
            print()
        
        print_success("All test predictions completed successfully")
        
    except Exception as e:
        print_error(f"Prediction testing failed: {e}")
        raise


def step5_compare_models(benchmark_dir):
    """
    Step 5: Compare Random Forest vs Gradient Boosting.
    
    Trains both model types and compares their performance
    to determine which works best for this task.
    """
    print_header("STEP 5: Compare Model Types")
    
    models = ['random_forest', 'gradient_boosting']
    results_comparison = {}
    
    print_info("Training and comparing both model types...\n")
    
    try:
        for model_type in models:
            print(f"{'─' * 70}")
            print(f"Training {model_type.replace('_', ' ').title()}...")
            print(f"{'─' * 70}\n")
            
            recommender = MLPipelineRecommender(model_type=model_type)
            results = recommender.train_from_benchmarks(benchmark_dir)
            
            results_comparison[model_type] = {
                'test_accuracy': results['test_score'],
                'cv_mean': results['cv_mean'],
                'cv_std': results['cv_std'],
                'train_time': results.get('train_time', 0)
            }
            
            print(f"  Test Accuracy:    {results['test_score']:.3f}")
            print(f"  CV Score:         {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
            if results.get('train_time'):
                print(f"  Training Time:    {results['train_time']:.2f}s")
            print()
        
        # Visualize comparison
        print_info("Creating comparison visualization...")
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Accuracy comparison
        x = np.arange(len(models))
        width = 0.35
        
        test_scores = [results_comparison[m]['test_accuracy'] for m in models]
        cv_scores = [results_comparison[m]['cv_mean'] for m in models]
        cv_errs = [results_comparison[m]['cv_std'] for m in models]
        
        axes[0].bar(x - width/2, test_scores, width, label='Test Accuracy', 
                   color='steelblue', edgecolor='black')
        axes[0].bar(x + width/2, cv_scores, width, label='CV Accuracy', 
                   yerr=cv_errs, capsize=5, color='coral', edgecolor='black')
        
        axes[0].set_xlabel('Model Type', fontsize=12)
        axes[0].set_ylabel('Accuracy', fontsize=12)
        axes[0].set_title('Model Accuracy Comparison', fontsize=14, fontweight='bold')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels([m.replace('_', ' ').title() for m in models])
        axes[0].legend(fontsize=11)
        axes[0].set_ylim([0, 1.1])
        axes[0].grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, (test_val, cv_val) in enumerate(zip(test_scores, cv_scores)):
            axes[0].text(i - width/2, test_val + 0.02, f'{test_val:.3f}', 
                        ha='center', va='bottom', fontsize=9)
            axes[0].text(i + width/2, cv_val + 0.02, f'{cv_val:.3f}', 
                        ha='center', va='bottom', fontsize=9)
        
        # Plot 2: Training time comparison (if available)
        train_times = [results_comparison[m].get('train_time', 0) for m in models]
        if any(train_times):
            axes[1].bar(x, train_times, color='mediumpurple', edgecolor='black')
            axes[1].set_xlabel('Model Type', fontsize=12)
            axes[1].set_ylabel('Training Time (seconds)', fontsize=12)
            axes[1].set_title('Training Time Comparison', fontsize=14, fontweight='bold')
            axes[1].set_xticks(x)
            axes[1].set_xticklabels([m.replace('_', ' ').title() for m in models])
            axes[1].grid(axis='y', alpha=0.3)
            
            # Add value labels
            for i, time_val in enumerate(train_times):
                if time_val > 0:
                    axes[1].text(i, time_val + 0.1, f'{time_val:.2f}s', 
                                ha='center', va='bottom', fontsize=9)
        else:
            axes[1].text(0.5, 0.5, 'Training time data\nnot available', 
                        ha='center', va='center', transform=axes[1].transAxes,
                        fontsize=12)
            axes[1].set_xticks([])
            axes[1].set_yticks([])
        
        plt.tight_layout()
        plt.savefig('figures/model_comparison.png', dpi=300, bbox_inches='tight')
        print_success("Saved comparison: figures/model_comparison.png")
        plt.close()
        
        # Determine and announce winner
        winner = max(results_comparison.items(), 
                    key=lambda x: x[1]['test_accuracy'])
        
        print(f"\n{'═' * 70}")
        print(f"🏆 BEST MODEL: {winner[0].replace('_', ' ').title()}")
        print(f"{'═' * 70}")
        print(f"Test Accuracy:    {winner[1]['test_accuracy']:.3f}")
        print(f"CV Score:         {winner[1]['cv_mean']:.3f} ± {winner[1]['cv_std']:.3f}")
        print(f"{'═' * 70}\n")
        
    except Exception as e:
        print_error(f"Model comparison failed: {e}")
        raise


def main():
    """Main workflow execution."""
    
    parser = argparse.ArgumentParser(
        description='Complete ML Recommender Workflow Demo for RAPTOR v2.1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic workflow with default settings
  python example_ml_workflow.py
  
  # Generate more training data
  python example_ml_workflow.py --n-datasets 500
  
  # Use existing data
  python example_ml_workflow.py --skip-generation --data-dir my_data/
  
  # Compare models
  python example_ml_workflow.py --compare-models
  
  # Use Gradient Boosting
  python example_ml_workflow.py --model-type gradient_boosting

For more information, visit:
  https://github.com/AyehBlk/RAPTOR
        """
    )
    
    parser.add_argument('--n-datasets', type=int, default=200,
                       help='Number of synthetic datasets to generate (default: 200)')
    parser.add_argument('--model-type', 
                       choices=['random_forest', 'gradient_boosting'],
                       default='random_forest',
                       help='Model type to train (default: random_forest)')
    parser.add_argument('--skip-generation', action='store_true',
                       help='Skip data generation and use existing data')
    parser.add_argument('--compare-models', action='store_true',
                       help='Compare Random Forest vs Gradient Boosting')
    parser.add_argument('--data-dir', default='ml_training_data',
                       help='Directory for training data (default: ml_training_data)')
    parser.add_argument('--model-dir', default='models',
                       help='Directory for saved models (default: models)')
    parser.add_argument('--figures-dir', default='figures',
                       help='Directory for figures (default: figures)')
    
    args = parser.parse_args()
    
    # Welcome message
    print("\n" + "=" * 70)
    print("🦖 RAPTOR ML RECOMMENDER - COMPLETE WORKFLOW")
    print("   RNA-seq Analysis Pipeline Testing & Optimization Resource")
    print("   Version 2.1.0")
    print("=" * 70)
    
    try:
        # Step 1: Generate or load data
        if not args.skip_generation:
            summary = step1_generate_data(
                n_datasets=args.n_datasets,
                output_dir=args.data_dir
            )
        else:
            print_header("STEP 1: Using Existing Data")
            print_info(f"Using data from: {args.data_dir}")
            
            # Verify data exists
            data_path = Path(args.data_dir)
            if not data_path.exists():
                print_error(f"Data directory not found: {args.data_dir}")
                print("Please generate data first or check the path.")
                sys.exit(1)
        
        # Step 2: Train model
        recommender, results = step2_train_model(
            benchmark_dir=args.data_dir,
            model_type=args.model_type,
            output_dir=args.model_dir
        )
        
        # Step 3: Visualize performance
        step3_visualize_performance(results, output_dir=args.figures_dir)
        
        # Step 4: Test predictions
        step4_test_predictions(recommender, n_test_cases=5)
        
        # Step 5: Compare models (optional)
        if args.compare_models:
            step5_compare_models(args.data_dir)
        
        # Success summary
        print_header("✅ WORKFLOW COMPLETE!")
        
        print("What you've accomplished:")
        print("  ✓ Generated/loaded synthetic training data")
        print("  ✓ Trained ML pipeline recommender")
        print("  ✓ Evaluated model performance")
        print("  ✓ Visualized results")
        print("  ✓ Tested predictions on new data")
        if args.compare_models:
            print("  ✓ Compared different model types")
        
        print("\nGenerated outputs:")
        print(f"  • Models:         {args.model_dir}/")
        print(f"  • Figures:        {args.figures_dir}/")
        print(f"  • Training data:  {args.data_dir}/")
        
        print("\nNext steps:")
        print("  • Use the trained model on real RNA-seq data")
        print("  • Integrate with RAPTOR's profiler module")
        print("  • Collect real benchmark data to improve the model")
        print("  • Deploy as part of RAPTOR CLI: raptor recommend")
        
        print("\n" + "=" * 70)
        print("Thank you for using RAPTOR! 🦖")
        print("For support: ayehbolouki1988@gmail.com")
        print("=" * 70 + "\n")
        
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print_error(f"Workflow failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
