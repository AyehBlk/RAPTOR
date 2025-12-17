#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Example Script: ML-Based Pipeline Recommendation

Demonstrates the machine learning-based pipeline recommendation system that
uses Random Forest or Gradient Boosting to predict optimal pipelines.

Features:
- Data profiling with 30+ features
- ML-based pipeline prediction
- Confidence scoring with explanations
- Feature importance visualization

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
    from raptor.ml_recommender import MLPipelineRecommender, FeatureExtractor
    from raptor.profiler import RNAseqDataProfiler
    from raptor.recommender import PipelineRecommender
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - ML Pipeline Recommendation         ║
    ║                                                              ║
    ║   Machine Learning-Powered Pipeline Selection                ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_profile():
    """Generate demonstration profile for testing."""
    np.random.seed(42)
    
    profile = {
        'design': {
            'n_samples': 12,
            'n_genes': 20000,
            'n_conditions': 2,
            'samples_per_condition': 6,
            'is_paired': False,
            'min_replicates': 6
        },
        'library_stats': {
            'mean': 25000000.0,
            'median': 24500000.0,
            'cv': 0.18,
            'range': 8000000.0,
            'skewness': 0.3
        },
        'count_distribution': {
            'zero_pct': 42.5,
            'low_count_pct': 58.2,
            'mean': 1250.0,
            'median': 450.0,
            'variance': 2500000.0
        },
        'expression_distribution': {
            'high_expr_genes': 2000,
            'medium_expr_genes': 6000,
            'low_expr_genes': 12000,
            'dynamic_range': 9.5
        },
        'biological_variation': {
            'bcv': 0.35,
            'dispersion_mean': 0.12,
            'dispersion_trend': -0.2,
            'outlier_genes': 150
        },
        'sequencing': {
            'total_reads': 300000000.0,
            'reads_per_gene': 1250.0,
            'depth_category': 'high',
            'mean_depth': 25000000.0
        },
        'complexity': {
            'score': 65.0,
            'noise_level': 0.7,
            'signal_strength': 0.65
        },
        'summary': {
            'difficulty': 'moderate',
            'recommended_approach': 'standard'
        }
    }
    
    return profile


def generate_demo_recommendation(profile):
    """Generate demonstration recommendation."""
    recommendations = [
        {
            'pipeline_id': 3,
            'pipeline_name': 'Salmon-edgeR',
            'confidence': 0.82,
            'score': 0.82,
            'reasons': [
                'High confidence recommendation (82%)',
                'Fast pseudo-alignment with robust statistics',
                'Efficient for high-depth sequencing',
                'Handles moderate biological variability well',
                'Well-suited for your sample size (12 samples)'
            ],
            'feature_contributions': [
                {'feature': 'bcv', 'importance': 0.18, 'value': 0.35},
                {'feature': 'n_samples', 'importance': 0.15, 'value': 12},
                {'feature': 'depth_category', 'importance': 0.12, 'value': 'high'},
                {'feature': 'zero_pct', 'importance': 0.10, 'value': 42.5}
            ]
        },
        {
            'pipeline_id': 1,
            'pipeline_name': 'STAR-RSEM-DESeq2',
            'confidence': 0.78,
            'score': 0.78,
            'reasons': [
                'Gold standard pipeline with excellent accuracy',
                'Well-suited for your sample size'
            ]
        },
        {
            'pipeline_id': 5,
            'pipeline_name': 'STAR-HTSeq-limma-voom',
            'confidence': 0.72,
            'score': 0.72,
            'reasons': [
                'Robust for complex experimental designs',
                'Handles biological variation well'
            ]
        }
    ]
    
    return recommendations


def display_recommendation(recommendation, rank=1):
    """Display a single recommendation with formatting."""
    print(f"\n{'='*60}")
    print(f"  RANK #{rank}: {recommendation['pipeline_name']} (Pipeline {recommendation['pipeline_id']})")
    print(f"{'='*60}")
    
    confidence = recommendation.get('confidence', recommendation.get('score', 0))
    
    # Confidence bar
    bar_length = 40
    filled = int(confidence * bar_length)
    bar = '█' * filled + '░' * (bar_length - filled)
    print(f"\n  Confidence: [{bar}] {confidence:.1%}")
    
    # Reasons
    print(f"\n  📋 Reasoning:")
    for reason in recommendation.get('reasons', []):
        print(f"     • {reason}")
    
    # Feature contributions (if available)
    if 'feature_contributions' in recommendation:
        print(f"\n  🔬 Key Factors:")
        for contrib in recommendation['feature_contributions'][:4]:
            importance = contrib['importance']
            bar_len = int(importance * 20)
            bar = '▓' * bar_len + '░' * (20 - bar_len)
            print(f"     {contrib['feature']:20s} [{bar}] {importance:.3f}")


def run_ml_recommendation(counts_file, metadata_file=None, model_type='random_forest', 
                          model_path=None, top_k=3, demo=False):
    """Run ML-based pipeline recommendation."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated data...")
        profile = generate_demo_profile()
        recommendations = generate_demo_recommendation(profile)
    else:
        # Load data
        print(f"\n📂 Loading data from: {counts_file}")
        counts = pd.read_csv(counts_file, index_col=0)
        print(f"   Loaded: {counts.shape[0]} genes × {counts.shape[1]} samples")
        
        metadata = None
        if metadata_file:
            metadata = pd.read_csv(metadata_file)
            print(f"   Metadata: {len(metadata)} samples")
        
        # Profile data
        print("\n🔍 Profiling data characteristics...")
        profiler = RNAseqDataProfiler(counts, metadata)
        profile = profiler.run_full_profile()
        
        # Get ML recommendation
        print(f"\n🤖 Running ML recommendation ({model_type})...")
        
        if model_path and Path(model_path).exists():
            recommender = MLPipelineRecommender(model_type=model_type)
            recommender.load_model(model_path)
        else:
            recommender = MLPipelineRecommender(model_type=model_type)
            print("   Note: Using untrained model. Train with benchmark data for better results.")
        
        # Get recommendation
        recommendation = recommender.recommend(profile, top_k=top_k)
        
        # Format as list for consistent display
        if isinstance(recommendation, dict):
            recommendations = [recommendation]
            if 'alternatives' in recommendation:
                recommendations.extend(recommendation['alternatives'])
        else:
            recommendations = recommendation[:top_k] if len(recommendation) > top_k else recommendation
    
    # Display results
    print("\n" + "="*60)
    print("  🦖 ML PIPELINE RECOMMENDATIONS")
    print("="*60)
    
    # Show profile summary
    print("\n📊 Data Profile Summary:")
    print(f"   Samples:          {profile['design']['n_samples']}")
    print(f"   Genes:            {profile['design']['n_genes']:,}")
    print(f"   BCV:              {profile['biological_variation']['bcv']:.3f}")
    print(f"   Sequencing depth: {profile['sequencing']['depth_category']}")
    print(f"   Zero inflation:   {profile['count_distribution']['zero_pct']:.1f}%")
    print(f"   Difficulty:       {profile['summary']['difficulty']}")
    
    # Show recommendations
    for i, rec in enumerate(recommendations[:top_k], 1):
        display_recommendation(rec, rank=i)
    
    # Prepare output
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'model_type': model_type,
        'profile_summary': {
            'n_samples': profile['design']['n_samples'],
            'n_genes': profile['design']['n_genes'],
            'bcv': profile['biological_variation']['bcv'],
            'depth_category': profile['sequencing']['depth_category'],
            'difficulty': profile['summary']['difficulty']
        },
        'recommendations': recommendations[:top_k],
        'full_profile': profile
    }
    
    return output


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR ML-Based Pipeline Recommendation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python 05_ml_recommendation.py --counts counts.csv
  
  # With metadata
  python 05_ml_recommendation.py --counts counts.csv --metadata samples.csv
  
  # Use Gradient Boosting model
  python 05_ml_recommendation.py --counts counts.csv --model gradient_boosting
  
  # Get top 5 recommendations
  python 05_ml_recommendation.py --counts counts.csv --top-k 5
  
  # Demo mode (no data required)
  python 05_ml_recommendation.py --demo
        """
    )
    
    parser.add_argument('--counts', help='Count matrix CSV file (genes × samples)')
    parser.add_argument('--metadata', help='Sample metadata CSV file')
    parser.add_argument('--model', choices=['random_forest', 'gradient_boosting'],
                        default='random_forest', help='ML model type')
    parser.add_argument('--model-path', help='Path to pre-trained model directory')
    parser.add_argument('--top-k', type=int, default=3, help='Number of recommendations')
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
    
    # Run recommendation
    results = run_ml_recommendation(
        counts_file=args.counts,
        metadata_file=args.metadata,
        model_type=args.model,
        model_path=args.model_path,
        top_k=args.top_k,
        demo=args.demo
    )
    
    # Save output
    output_file = args.output or 'ml_recommendation_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n💾 Results saved to: {output_file}")
    print("\n" + "="*60)
    print("  Making free science for everybody around the world 🌍")
    print("="*60 + "\n")


if __name__ == '__main__':
    main()
