"""
RAPTOR Command-Line Interface

Provides the `raptor` command with subcommands for profiling, benchmarking,
and generating reports.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import click
import sys
import logging
from pathlib import Path

from raptor import __version__
from raptor.profiler import RNAseqDataProfiler
from raptor.recommender import PipelineRecommender
from raptor.benchmark import PipelineBenchmark
from raptor.simulate import DataSimulator
from raptor.report import ReportGenerator

# Configure logging
logger = logging.getLogger(__name__)


# Main CLI group
@click.group()
@click.version_option(version=__version__, prog_name='RAPTOR')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.option('--quiet', '-q', is_flag=True, help='Suppress non-error output')
def main(verbose, quiet):
    """
    🦖 RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource
    
    A comprehensive framework for benchmarking RNA-seq pipelines and getting
    intelligent, data-driven pipeline recommendations.
    
    Examples:
    
        # Get pipeline recommendation
        raptor profile --counts data.csv --output report.html
        
        # Run full benchmarking
        raptor compare --data fastq/ --output results/
        
        # Generate report
        raptor report --results results/ --output report.html
    
    For more information: https://github.com/AyehBlk/RAPTOR
    """
    # Configure logging level
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    elif quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)


# =============================================================================
# Profile Command - Data Profiling and Recommendation
# =============================================================================

@main.command()
@click.option('--counts', '-c', type=click.Path(exists=True), required=True,
              help='Count matrix file (CSV or TSV)')
@click.option('--metadata', '-m', type=click.Path(exists=True),
              help='Sample metadata file (optional)')
@click.option('--output', '-o', type=click.Path(), default='raptor_profile.html',
              help='Output report file (HTML or JSON)')
@click.option('--priority', type=click.Choice(['accuracy', 'speed', 'memory', 'balanced']),
              default='balanced', help='Optimization priority')
@click.option('--format', type=click.Choice(['html', 'json', 'text']),
              default='html', help='Output format')
def profile(counts, metadata, output, priority, format):
    """
    Profile RNA-seq data and get pipeline recommendation.
    
    Analyzes count matrix characteristics and recommends the optimal pipeline
    based on data properties like library size variation, zero-inflation,
    biological variation, and sample size.
    
    Examples:
    
        # Basic profiling
        raptor profile --counts counts.csv
        
        # With metadata and specific priority
        raptor profile --counts counts.csv --metadata meta.csv --priority accuracy
        
        # JSON output for automation
        raptor profile --counts counts.csv --format json --output rec.json
    """
    click.echo("🔍 Profiling RNA-seq data...")
    
    try:
        # Load data
        import pandas as pd
        counts_df = pd.read_csv(counts, index_col=0)
        
        metadata_df = None
        if metadata:
            metadata_df = pd.read_csv(metadata)
        
        # Profile the data
        click.echo(f"📊 Analyzing {counts_df.shape[0]} genes across {counts_df.shape[1]} samples...")
        profiler = RNAseqDataProfiler(counts_df, metadata_df)
        profile_data = profiler.run_full_profile()
        
        # Get recommendation
        click.echo(f"🎯 Generating recommendation (priority: {priority})...")
        recommender = PipelineRecommender(profile_data)
        recommendation = recommender.get_recommendation(priority=priority)
        
        # Generate output
        if format == 'html':
            from raptor.report import generate_profile_report
            generate_profile_report(profile_data, recommendation, output)
            click.echo(f"✅ Profile report saved to: {output}")
        elif format == 'json':
            import json
            with open(output, 'w') as f:
                json.dump({
                    'profile': profile_data,
                    'recommendation': recommendation
                }, f, indent=2)
            click.echo(f"✅ Profile data saved to: {output}")
        else:  # text
            _print_recommendation(recommendation)
        
        click.echo(f"\n🦖 Recommended: {recommendation['primary']['pipeline_name']}")
        click.echo(f"   Score: {recommendation['primary']['score']:.1f}/200")
        
    except Exception as e:
        logger.error(f"Error during profiling: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Compare Command - Full Pipeline Benchmarking
# =============================================================================

@main.command()
@click.option('--data', '-d', type=click.Path(exists=True), required=True,
              help='Input data directory (FASTQ files or count matrix)')
@click.option('--output', '-o', type=click.Path(), default='raptor_benchmark',
              help='Output directory for results')
@click.option('--pipelines', '-p', default='all',
              help='Pipelines to run (comma-separated IDs or "all")')
@click.option('--threads', '-t', type=int, default=8,
              help='Number of threads to use')
@click.option('--memory', type=str, default='32G',
              help='Maximum memory (e.g., "32G")')
@click.option('--reference', type=click.Path(exists=True),
              help='Reference genome/transcriptome')
def compare(data, output, pipelines, threads, memory, reference):
    """
    Benchmark multiple RNA-seq pipelines on your data.
    
    Runs complete RNA-seq analysis using multiple pipelines and compares
    their performance in terms of accuracy, speed, and resource usage.
    
    Examples:
    
        # Benchmark all pipelines
        raptor compare --data fastq/ --output results/
        
        # Benchmark specific pipelines
        raptor compare --data fastq/ --pipelines 1,3,5 --threads 16
        
        # With custom reference
        raptor compare --data fastq/ --reference genome.fa --output results/
    """
    click.echo("⚖️  Starting pipeline benchmarking...")
    
    try:
        # Parse pipeline selection
        if pipelines == 'all':
            pipeline_ids = list(range(1, 9))  # All 8 pipelines
        else:
            pipeline_ids = [int(p.strip()) for p in pipelines.split(',')]
        
        click.echo(f"📋 Running {len(pipeline_ids)} pipelines: {pipeline_ids}")
        click.echo(f"💻 Using {threads} threads and {memory} memory")
        
        # Create benchmark object
        benchmark = PipelineBenchmark(
            data_dir=data,
            output_dir=output,
            threads=threads,
            memory=memory,
            reference=reference
        )
        
        # Run benchmarking
        results = benchmark.run_pipelines(pipeline_ids)
        
        # Save results
        benchmark.save_results(results)
        
        click.echo(f"✅ Benchmarking complete! Results saved to: {output}")
        click.echo(f"📊 Generate report with: raptor report --results {output}")
        
    except Exception as e:
        logger.error(f"Error during benchmarking: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Run Command - Run Specific Pipeline
# =============================================================================

@main.command()
@click.option('--pipeline', '-p', type=int, required=True,
              help='Pipeline ID (1-8)')
@click.option('--data', '-d', type=click.Path(exists=True), required=True,
              help='Input data directory (FASTQ files)')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output directory')
@click.option('--threads', '-t', type=int, default=8,
              help='Number of threads')
@click.option('--reference', type=click.Path(exists=True),
              help='Reference genome/transcriptome')
def run(pipeline, data, output, threads, reference):
    """
    Run a specific RNA-seq pipeline.
    
    Execute a single pipeline on your data. Useful after getting a
    recommendation from the profile command.
    
    Examples:
    
        # Run recommended pipeline
        raptor run --pipeline 3 --data fastq/ --output results/
        
        # With custom threads and reference
        raptor run --pipeline 1 --data fastq/ --threads 16 --reference genome.fa
    """
    click.echo(f"🚀 Running Pipeline {pipeline}...")
    
    try:
        benchmark = PipelineBenchmark(
            data_dir=data,
            output_dir=output,
            threads=threads,
            reference=reference
        )
        
        results = benchmark.run_single_pipeline(pipeline)
        
        click.echo(f"✅ Pipeline {pipeline} completed!")
        click.echo(f"📂 Results saved to: {output}")
        
    except Exception as e:
        logger.error(f"Error running pipeline: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Report Command - Generate Comparison Report
# =============================================================================

@main.command()
@click.option('--results', '-r', type=click.Path(exists=True), required=True,
              help='Benchmark results directory')
@click.option('--output', '-o', type=click.Path(), default='raptor_report.html',
              help='Output report file')
@click.option('--format', type=click.Choice(['html', 'pdf', 'markdown']),
              default='html', help='Report format')
def report(results, output, format):
    """
    Generate comparison report from benchmark results.
    
    Creates a comprehensive report with visualizations comparing pipeline
    performance across accuracy, speed, and resource usage metrics.
    
    Examples:
    
        # Generate HTML report
        raptor report --results benchmark_results/
        
        # Generate PDF report
        raptor report --results results/ --format pdf --output report.pdf
    """
    click.echo("📊 Generating comparison report...")
    
    try:
        reporter = ReportGenerator(results)
        reporter.generate_report(output_file=output, format=format)
        
        click.echo(f"✅ Report generated: {output}")
        
    except Exception as e:
        logger.error(f"Error generating report: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Simulate Command - Generate Simulated Data
# =============================================================================

@main.command()
@click.option('--n-genes', type=int, default=2000,
              help='Number of genes')
@click.option('--n-samples', type=int, default=6,
              help='Number of samples')
@click.option('--n-de', type=int, default=400,
              help='Number of DE genes')
@click.option('--fold-changes', default='0.5,2',
              help='Fold changes (comma-separated)')
@click.option('--output', '-o', type=click.Path(), default='simulated_data',
              help='Output directory')
def simulate(n_genes, n_samples, n_de, fold_changes, output):
    """
    Generate simulated RNA-seq data for testing.
    
    Creates synthetic RNA-seq data with known ground truth for benchmarking
    and validation purposes.
    
    Examples:
    
        # Basic simulation
        raptor simulate --output sim_data/
        
        # Custom parameters
        raptor simulate --n-genes 5000 --n-samples 12 --n-de 1000
    """
    click.echo("🧬 Generating simulated RNA-seq data...")
    
    try:
        # Parse fold changes
        fold_changes_list = [float(fc.strip()) for fc in fold_changes.split(',')]
        
        simulator = DataSimulator(
            n_genes=n_genes,
            n_samples=n_samples,
            n_de=n_de,
            fold_changes=fold_changes_list
        )
        
        simulator.generate_data(output)
        
        click.echo(f"✅ Simulated data generated: {output}")
        click.echo(f"📊 {n_genes} genes, {n_samples} samples, {n_de} DE genes")
        
    except Exception as e:
        logger.error(f"Error during simulation: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Demo Command - Quick Demonstration
# =============================================================================

@main.command()
def demo():
    """
    Run a quick demo of RAPTOR functionality.
    
    Generates small simulated dataset, profiles it, and shows recommendation.
    Perfect for testing installation and seeing RAPTOR in action.
    
    Example:
    
        raptor demo
    """
    click.echo("🎬 Running RAPTOR demo...")
    click.echo("=" * 60)
    
    try:
        import tempfile
        import pandas as pd
        import numpy as np
        
        # Generate small test dataset
        click.echo("\n1️⃣  Generating test data...")
        np.random.seed(42)
        n_genes, n_samples = 100, 6
        counts = pd.DataFrame(
            np.random.negative_binomial(10, 0.1, (n_genes, n_samples)),
            columns=[f'Sample_{i+1}' for i in range(n_samples)]
        )
        
        click.echo(f"   Created {n_genes} genes × {n_samples} samples")
        
        # Profile the data
        click.echo("\n2️⃣  Profiling data...")
        profiler = RNAseqDataProfiler(counts)
        profile_data = profiler.run_full_profile()
        
        click.echo(f"   Library size CV: {profile_data['library_stats']['cv']:.2f}")
        click.echo(f"   Zero inflation: {profile_data['count_distribution']['zero_pct']:.1f}%")
        click.echo(f"   BCV: {profile_data['biological_variation']['bcv']:.2f}")
        
        # Get recommendation
        click.echo("\n3️⃣  Getting recommendation...")
        recommender = PipelineRecommender(profile_data)
        recommendation = recommender.get_recommendation()
        
        click.echo(f"\n{'=' * 60}")
        click.echo(f"🎯 RECOMMENDED: {recommendation['primary']['pipeline_name']}")
        click.echo(f"   Score: {recommendation['primary']['score']:.1f}/200")
        click.echo(f"\n   Reasoning:")
        for reason in recommendation['primary']['reasoning'][:3]:
            click.echo(f"   ✅ {reason}")
        
        click.echo(f"\n{'=' * 60}")
        click.echo("✅ Demo complete!")
        click.echo("\nNext steps:")
        click.echo("  • Try: raptor profile --counts YOUR_DATA.csv")
        click.echo("  • Documentation: https://github.com/AyehBlk/RAPTOR")
        
    except Exception as e:
        logger.error(f"Error during demo: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Helper Functions
# =============================================================================

def _print_recommendation(recommendation):
    """Print recommendation in text format."""
    rec = recommendation['primary']
    click.echo("\n" + "=" * 60)
    click.echo(f"🎯 RECOMMENDED PIPELINE")
    click.echo("=" * 60)
    click.echo(f"\nPipeline: {rec['pipeline_name']}")
    click.echo(f"ID: {rec['pipeline_id']}")
    click.echo(f"Score: {rec['score']:.1f}/200")
    click.echo(f"\nReasoning:")
    for i, reason in enumerate(rec['reasoning'], 1):
        click.echo(f"  {i}. {reason}")
    
    if 'alternatives' in recommendation:
        click.echo(f"\n{'=' * 60}")
        click.echo("Alternative Options:")
        for alt in recommendation['alternatives'][:2]:
            click.echo(f"\n  • {alt['pipeline_name']} (Score: {alt['score']:.1f})")


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == '__main__':
    main()
