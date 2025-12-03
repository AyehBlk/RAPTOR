#!/usr/bin/env python3

"""
01_run_all_pipelines_python.py
Python wrapper for running all pipelines using raptor.benchmark module

This script provides a pure Python alternative to the bash script,
directly using raptor.benchmark.PipelineBenchmark.

Usage: python 01_run_all_pipelines_python.py <data_dir> <output_dir> [options]

Author: Ayeh Bolouki
Organization: RAPTOR Project
License: MIT
"""

import sys
import argparse
from pathlib import Path

try:
    from raptor.benchmark import PipelineBenchmark
    RAPTOR_AVAILABLE = True
except ImportError:
    RAPTOR_AVAILABLE = False
    print("Error: RAPTOR package not found.")
    print("Install with: pip install raptor-rnaseq")
    sys.exit(1)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Run RAPTOR RNA-seq pipelines using benchmark module',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all 8 pipelines
  python 01_run_all_pipelines_python.py data/fastq results/benchmark --reference refs/
  
  # Run specific pipelines
  python 01_run_all_pipelines_python.py data/fastq results/ --pipelines 1 3 4
  
  # With custom resources
  python 01_run_all_pipelines_python.py data/fastq results/ --threads 16 --memory 64G

For more information: https://github.com/AyehBlk/RAPTOR
        """
    )
    
    parser.add_argument('data_dir',
                       help='Directory containing input data (FASTQ or counts)')
    parser.add_argument('output_dir',
                       help='Base output directory for results')
    parser.add_argument('-r', '--reference',
                       help='Path to reference genome/transcriptome directory')
    parser.add_argument('-p', '--pipelines',
                       type=int,
                       nargs='+',
                       default=[1, 2, 3, 4, 5, 6, 7, 8],
                       help='Pipeline IDs to run (default: all 8)')
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=8,
                       help='Number of threads (default: 8)')
    parser.add_argument('-m', '--memory',
                       default='32G',
                       help='Maximum memory (default: 32G)')
    parser.add_argument('-v', '--verbose',
                       action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.data_dir).exists():
        print(f"Error: Data directory not found: {args.data_dir}")
        sys.exit(1)
    
    # Create benchmark
    print("=== RAPTOR Pipeline Benchmark ===\n")
    print(f"Data: {args.data_dir}")
    print(f"Output: {args.output_dir}")
    print(f"Reference: {args.reference}")
    print(f"Pipelines: {args.pipelines}")
    print(f"Resources: {args.threads} threads, {args.memory} memory\n")
    
    benchmark = PipelineBenchmark(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        threads=args.threads,
        memory=args.memory,
        reference=args.reference
    )
    
    # Run pipelines
    print("Running pipelines...\n")
    results = benchmark.run_pipelines(args.pipelines)
    
    # Save results
    benchmark.save_results(results)
    
    # Summary
    print("\n=== Summary ===")
    success_count = sum(1 for r in results.values() if r.get('status') == 'success')
    print(f"Successful: {success_count}/{len(results)}")
    print(f"Failed: {len(results) - success_count}/{len(results)}")
    
    print("\nResults:")
    for pipeline_id, result in results.items():
        status = "✓" if result.get('status') == 'success' else "✗"
        pipeline_name = result.get('pipeline', f'Pipeline {pipeline_id}')
        runtime = result.get('runtime', 0)
        print(f"  {status} {pipeline_name}: {runtime:.1f}s")
    
    print(f"\nResults saved to: {args.output_dir}/benchmark_results.json")
    print("\nNext steps:")
    print(f"  1. Compare results: Rscript scripts/03_compare_results.R {args.output_dir}")
    print(f"  2. Visualize: Rscript scripts/04_visualize_comparison.R {args.output_dir}")
    
    if success_count == len(results):
        print("\n✓ All pipelines completed successfully! 🦖")
        return 0
    else:
        print("\n⚠ Some pipelines failed. Check logs for details.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
