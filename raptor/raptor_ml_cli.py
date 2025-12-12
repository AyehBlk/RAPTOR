#!/usr/bin/env python3

"""
Enhanced RAPTOR CLI with ML Recommender

Command-line interface for RAPTOR with integrated machine learning-based
pipeline recommendation.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import sys
import argparse
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║                     🦖 RAPTOR v2.1.0                        ║
    ║   RNA-seq Analysis Pipeline Testing & Optimization Resource ║
    ║                                                              ║
    ║          NOW WITH ML-POWERED RECOMMENDATIONS! 🤖            ║
    ║                                                              ║
    ║              Created by Ayeh Bolouki                         ║
    ║  University of Namur & GIGA-Neurosciences, University of    ║
    ║                    Liège, Belgium                            ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def cmd_profile(args):
    """Profile RNA-seq data and get recommendations."""
    print("\n=== Data Profiling ===\n")
    
    # Check if count file exists
    if not Path(args.counts).exists():
        logger.error(f"Count file not found: {args.counts}")
        sys.exit(1)
    
    # Import modules
    try:
        from raptor.profiler import RNAseqDataProfiler
        from raptor.recommender import PipelineRecommender
        
        # Try to import ML recommender
        try:
            from raptor.ml_recommender import MLPipelineRecommender
            ML_AVAILABLE = True
        except ImportError:
            ML_AVAILABLE = False
            logger.warning("ML recommender not available, using rule-based system")
        
    except ImportError as e:
        logger.error(f"RAPTOR modules not found: {e}")
        logger.info("Install with: pip install raptor-rnaseq")
        sys.exit(1)
    
    # Load data
    import pandas as pd
    
    logger.info(f"Loading count matrix: {args.counts}")
    counts = pd.read_csv(args.counts, index_col=0)
    
    metadata = None
    if args.metadata:
        logger.info(f"Loading metadata: {args.metadata}")
        metadata = pd.read_csv(args.metadata)
    
    # Profile data
    logger.info("Running data profiling...")
    profiler = RNAseqDataProfiler(counts, metadata)
    profile = profiler.run_full_profile()
    
    # Display profile summary
    print("\n--- Profile Summary ---")
    print(f"Samples: {profile['design']['n_samples']}")
    print(f"Genes: {profile['design']['n_genes']}")
    print(f"Mean library size: {profile['library_stats']['mean']:,.0f}")
    print(f"Library size CV: {profile['library_stats']['cv']:.2f}")
    print(f"Zero percentage: {profile['count_distribution']['zero_pct']:.1f}%")
    print(f"Biological variation (BCV): {profile['biological_variation']['bcv']:.3f}")
    print(f"Sequencing depth: {profile['sequencing']['depth_category']}")
    print(f"Data difficulty: {profile['summary']['difficulty']}")
    
    # Get recommendation
    print("\n--- Pipeline Recommendation ---\n")
    
    if ML_AVAILABLE and args.use_ml and Path(args.ml_model).exists():
        # Use ML recommender
        logger.info("Using ML-based recommendation")
        recommender = MLPipelineRecommender(model_type='random_forest')
        recommender.load_model(args.ml_model)
        recommendation = recommender.recommend(profile, top_k=3)
        
        print(f"🦖 ML RECOMMENDATION (Confidence: {recommendation['confidence']:.1%}):")
        print(f"   Pipeline {recommendation['pipeline_id']}: {recommendation['pipeline_name']}")
        print(f"\n   Reasons:")
        for reason in recommendation['reasons']:
            print(f"     • {reason}")
        
        if recommendation['feature_contributions']:
            print(f"\n   Key factors:")
            for contrib in recommendation['feature_contributions'][:3]:
                print(f"     • {contrib['feature']}: {contrib['importance']:.3f}")
        
        if recommendation['alternatives']:
            print(f"\n   Alternatives:")
            for alt in recommendation['alternatives']:
                print(f"     • Pipeline {alt['pipeline_id']} ({alt['confidence']:.1%})")
    
    else:
        # Use rule-based recommender
        logger.info("Using rule-based recommendation")
        recommender = PipelineRecommender(profile)
        recommendation = recommender.recommend()
        
        print(f"🦖 RECOMMENDED: Pipeline {recommendation['pipeline_id']}")
        print(f"   {recommendation['pipeline_name']}")
        print(f"\n   Reasons:")
        for reason in recommendation['reasons']:
            print(f"     • {reason}")
    
    # Save profile if requested
    if args.output:
        import json
        output_path = Path(args.output)
        
        # Add recommendation to profile
        profile['recommendation'] = recommendation
        
        with open(output_path, 'w') as f:
            json.dump(profile, f, indent=2)
        
        logger.info(f"Profile saved to: {output_path}")
    
    print()


def cmd_compare(args):
    """Run and compare multiple pipelines."""
    print("\n=== Pipeline Comparison ===\n")
    
    try:
        from raptor.benchmark import PipelineBenchmark
    except ImportError:
        logger.error("RAPTOR benchmark module not found")
        sys.exit(1)
    
    # Create benchmark
    benchmark = PipelineBenchmark(
        data_dir=args.data,
        output_dir=args.output,
        threads=args.threads,
        memory=args.memory,
        reference=args.reference
    )
    
    # Parse pipeline IDs
    if args.pipelines == 'all':
        pipeline_ids = list(range(1, 9))
    else:
        pipeline_ids = [int(x) for x in args.pipelines.split(',')]
    
    logger.info(f"Running pipelines: {pipeline_ids}")
    
    # Run pipelines
    results = benchmark.run_pipelines(pipeline_ids)
    
    # Save results
    benchmark.save_results(results)
    
    # Summary
    print("\n--- Comparison Summary ---")
    success_count = sum(1 for r in results.values() if r.get('status') == 'success')
    print(f"Successful: {success_count}/{len(results)}")
    
    for pipeline_id, result in results.items():
        status = "✓" if result.get('status') == 'success' else "✗"
        name = result.get('pipeline', f'Pipeline {pipeline_id}')
        runtime = result.get('runtime', 0)
        print(f"  {status} {name}: {runtime:.1f}s")
    
    print()


def cmd_train_ml(args):
    """Train ML recommender."""
    print("\n=== Train ML Recommender ===\n")
    
    try:
        from raptor.ml_recommender import MLPipelineRecommender
    except ImportError:
        logger.error("ML recommender module not found")
        logger.info("Ensure RAPTOR is installed: pip install raptor-rnaseq[ml]")
        sys.exit(1)
    
    # Check if benchmark directory exists
    if not Path(args.benchmark_dir).exists():
        logger.error(f"Benchmark directory not found: {args.benchmark_dir}")
        logger.info("Generate synthetic data with: raptor generate-data")
        sys.exit(1)
    
    # Create recommender
    logger.info(f"Training {args.model_type} model...")
    recommender = MLPipelineRecommender(model_type=args.model_type)
    
    # Train
    results = recommender.train_from_benchmarks(
        benchmark_dir=args.benchmark_dir,
        performance_metric=args.metric
    )
    
    # Display results
    print("\n--- Training Results ---")
    print(f"Model: {results['model_type']}")
    print(f"Samples: {results['n_samples']}")
    print(f"Test accuracy: {results['test_score']:.3f}")
    print(f"CV score: {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
    
    # Save model
    recommender.save_model(args.output)
    logger.info(f"Model saved to: {args.output}")
    
    print()


def cmd_generate_data(args):
    """Generate synthetic training data."""
    print("\n=== Generate Synthetic Data ===\n")
    
    try:
        from raptor.synthetic_benchmarks import generate_training_data
    except ImportError:
        logger.error("Synthetic benchmark module not found")
        logger.info("Ensure RAPTOR is installed: pip install raptor-rnaseq[ml]")
        sys.exit(1)
    
    logger.info(f"Generating {args.n_datasets} synthetic datasets...")
    
    summary = generate_training_data(
        n_datasets=args.n_datasets,
        output_dir=args.output,
        seed=args.seed
    )
    
    logger.info(f"Data saved to: {summary['output_dir']}")
    
    print()


def cmd_simulate(args):
    """Simulate RNA-seq data."""
    print("\n=== Simulate RNA-seq Data ===\n")
    
    try:
        from raptor.simulate import quick_simulate
    except ImportError:
        logger.error("RAPTOR simulate module not found")
        sys.exit(1)
    
    logger.info(f"Generating {args.size} dataset...")
    
    summary = quick_simulate(args.output, size=args.size)
    
    logger.info(f"Simulated data saved to: {args.output}")
    print(f"\nGenerated:")
    print(f"  • {summary['n_genes']} genes")
    print(f"  • {summary['n_samples']} samples")
    print(f"  • {summary['n_de']} DE genes")
    
    print()


def main():
    """Main CLI."""
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR - RNA-seq Analysis Pipeline Testing & Optimization Resource',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Profile data and get ML recommendation
  raptor profile --counts data.csv --use-ml --ml-model models/
  
  # Compare pipelines
  raptor compare --data fastq/ --output results/ --pipelines 1,3,4
  
  # Train ML recommender
  raptor train-ml --benchmark-dir ml_training_data/ --output models/
  
  # Generate synthetic training data
  raptor generate-data --n-datasets 200 --output ml_training_data/

For more information: https://github.com/AyehBlk/RAPTOR
        """
    )
    
    # Create subparsers
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Profile command
    profile_parser = subparsers.add_parser('profile', help='Profile data and get recommendations')
    profile_parser.add_argument('--counts', required=True, help='Count matrix CSV file')
    profile_parser.add_argument('--metadata', help='Sample metadata CSV file')
    profile_parser.add_argument('--output', '-o', help='Output JSON file')
    profile_parser.add_argument('--use-ml', action='store_true', 
                               help='Use ML-based recommendation')
    profile_parser.add_argument('--ml-model', default='models/',
                               help='Path to ML model directory')
    
    # Compare command
    compare_parser = subparsers.add_parser('compare', help='Compare pipelines')
    compare_parser.add_argument('--data', required=True, help='Input data directory')
    compare_parser.add_argument('--output', required=True, help='Output directory')
    compare_parser.add_argument('--reference', help='Reference genome directory')
    compare_parser.add_argument('--pipelines', default='all',
                               help='Pipeline IDs (e.g., "1,3,4" or "all")')
    compare_parser.add_argument('--threads', type=int, default=8, help='Number of threads')
    compare_parser.add_argument('--memory', default='32G', help='Memory limit')
    
    # Train ML command
    train_parser = subparsers.add_parser('train-ml', help='Train ML recommender')
    train_parser.add_argument('--benchmark-dir', required=True,
                             help='Directory with benchmark results')
    train_parser.add_argument('--output', default='models/', help='Output directory')
    train_parser.add_argument('--model-type', choices=['random_forest', 'gradient_boosting'],
                             default='random_forest', help='Model type')
    train_parser.add_argument('--metric', default='f1_score',
                             choices=['f1_score', 'accuracy', 'runtime', 'combined'],
                             help='Performance metric to optimize')
    
    # Generate data command
    generate_parser = subparsers.add_parser('generate-data',
                                           help='Generate synthetic training data')
    generate_parser.add_argument('--n-datasets', type=int, default=200,
                                help='Number of datasets to generate')
    generate_parser.add_argument('--output', default='ml_training_data/',
                                help='Output directory')
    generate_parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    # Simulate command
    simulate_parser = subparsers.add_parser('simulate', help='Simulate RNA-seq data')
    simulate_parser.add_argument('--output', required=True, help='Output directory')
    simulate_parser.add_argument('--size', choices=['small', 'medium', 'large'],
                                default='small', help='Dataset size')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Show banner
    print_banner()
    
    # Execute command
    if args.command == 'profile':
        cmd_profile(args)
    elif args.command == 'compare':
        cmd_compare(args)
    elif args.command == 'train-ml':
        cmd_train_ml(args)
    elif args.command == 'generate-data':
        cmd_generate_data(args)
    elif args.command == 'simulate':
        cmd_simulate(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
