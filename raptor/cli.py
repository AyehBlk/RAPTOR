"""
RAPTOR Command-Line Interface

Provides the `raptor` command with subcommands for profiling, benchmarking,
pipeline execution, and quick quantification.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
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
from raptor.sample_sheet import SampleSheet, auto_detect_samples

# Configure logging
logger = logging.getLogger(__name__)


# =============================================================================
# Main CLI Group
# =============================================================================

@click.group()
@click.version_option(version=__version__, prog_name='RAPTOR')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.option('--quiet', '-q', is_flag=True, help='Suppress non-error output')
def main(verbose, quiet):
    """
    🦖 RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource
    
    A comprehensive framework for benchmarking RNA-seq pipelines and getting
    intelligent, data-driven pipeline recommendations.
    
    \b
    WORKFLOW:
    1. quick-count  : Generate count matrix from FASTQ (fast)
    2. qc           : Quality assessment and outlier detection
    3. profile      : Analyze data and get ML recommendation
    4. run          : Execute recommended pipeline
    
    \b
    EXAMPLES:
    
        # Quick count from FASTQ files
        raptor quick-count --method salmon --sample-sheet samples.csv --index idx/ --output counts/
        
        # Get pipeline recommendation
        raptor profile --counts counts.csv --output report.html
        
        # Run full benchmarking
        raptor compare --data fastq/ --output results/
    
    For more information: https://github.com/AyehBlk/RAPTOR
    """
    # Configure logging level
    if verbose:
        logging.basicConfig(level=logging.DEBUG, 
                          format='%(asctime)s - %(levelname)s - %(message)s')
    elif quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO,
                          format='%(asctime)s - %(levelname)s - %(message)s')


# =============================================================================
# Quick-Count Command (Step 1: FASTQ → Count Matrix)
# =============================================================================

@main.command('quick-count')
@click.option('--method', '-m', 
              type=click.Choice(['salmon', 'kallisto']),
              required=True,
              help='Quantification method (salmon or kallisto)')
@click.option('--sample-sheet', '-s',
              type=click.Path(exists=True),
              required=True,
              help='Sample sheet CSV file with FASTQ paths')
@click.option('--index', '-i',
              type=click.Path(exists=True),
              required=True,
              help='Path to Salmon index directory or Kallisto index file (.idx)')
@click.option('--output', '-o',
              type=click.Path(),
              required=True,
              help='Output directory for count matrix')
@click.option('--threads', '-t',
              type=int,
              default=8,
              help='Number of threads (default: 8)')
@click.option('--gene-map', '-g',
              type=click.Path(exists=True),
              default=None,
              help='Gene-to-transcript mapping file (tx2gene.csv) for gene-level counts')
@click.option('--lib-type', '-l',
              type=str,
              default='A',
              help='Salmon library type (default: A for auto-detect)')
@click.option('--fragment-length',
              type=int,
              default=200,
              help='Fragment length mean for single-end reads (Kallisto requires this)')
@click.option('--fragment-sd',
              type=int,
              default=20,
              help='Fragment length SD for single-end reads (Kallisto requires this)')
@click.option('--keep-quant',
              is_flag=True,
              default=False,
              help='Keep individual sample quantification files')
@click.option('--validate-only',
              is_flag=True,
              default=False,
              help='Only validate inputs without running quantification')
def quick_count(method, sample_sheet, index, output, threads, gene_map,
                lib_type, fragment_length, fragment_sd, keep_quant, validate_only):
    """
    Generate count matrix from FASTQ files (fast).
    
    Uses Salmon or Kallisto pseudo-alignment for quick quantification.
    Optimized for speed to assess data quality before full analysis.
    
    \b
    EXAMPLES:
    
        # Paired-end with Salmon (recommended)
        raptor quick-count \\
          --method salmon \\
          --sample-sheet samples.csv \\
          --index salmon_index/ \\
          --output quick_counts/
        
        # Single-end with Kallisto
        raptor quick-count \\
          --method kallisto \\
          --sample-sheet samples.csv \\
          --index kallisto.idx \\
          --fragment-length 200 \\
          --output quick_counts/
    
    \b
    SAMPLE SHEET FORMAT (paired-end):
        sample_id,condition,fastq_r1,fastq_r2
        Sample1,Control,Sample1_R1.fq.gz,Sample1_R2.fq.gz
    
    \b
    OUTPUT FILES:
        counts.csv      - Raw count matrix
        tpm.csv         - TPM normalized matrix
        sample_info.csv - Sample QC metrics
    """
    click.echo("")
    click.echo("=" * 60)
    click.echo(f"🦖 RAPTOR Quick-Count ({method.upper()})")
    click.echo("=" * 60)
    click.echo("")
    
    # Validate inputs
    click.echo("📋 Validating inputs...")
    
    # Check index
    index_path = Path(index)
    if method == 'salmon':
        if not index_path.is_dir():
            click.echo(f"❌ Salmon index should be a directory: {index}", err=True)
            sys.exit(1)
        version_file = index_path / "versionInfo.json"
        if not version_file.exists():
            click.echo(f"⚠️  Warning: Salmon index may not be valid (missing versionInfo.json)")
    elif method == 'kallisto':
        if not index_path.is_file():
            click.echo(f"❌ Kallisto index should be a .idx file: {index}", err=True)
            sys.exit(1)
    
    # Create output directory
    output_path = Path(output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load and validate sample sheet
    try:
        sheet = SampleSheet(sample_sheet)
        is_valid, errors = sheet.validate(check_files=True)
        
        if not is_valid:
            click.echo("❌ Sample sheet validation failed:", err=True)
            for err in errors:
                click.echo(f"   • {err}", err=True)
            sys.exit(1)
        
        click.echo(f"   ✓ Sample sheet: {len(sheet)} samples")
        click.echo(f"   ✓ Read type: {'paired-end' if sheet.is_paired else 'single-end'}")
        
        if method == 'kallisto' and not sheet.is_paired:
            click.echo(f"   ✓ Fragment length: {fragment_length} ± {fragment_sd}")
        
    except Exception as e:
        click.echo(f"❌ Error validating sample sheet: {e}", err=True)
        sys.exit(1)
    
    click.echo(f"   ✓ Index: {index}")
    click.echo(f"   ✓ Output: {output}")
    click.echo(f"   ✓ Threads: {threads}")
    
    if gene_map:
        click.echo(f"   ✓ Gene map: {gene_map} (gene-level counts)")
    
    click.echo("")
    
    if validate_only:
        click.echo("✅ Validation complete!")
        return
    
    # Run quantification
    click.echo(f"🚀 Running {method.upper()} quantification...")
    click.echo("")
    
    try:
        if method == 'salmon':
            from raptor.pipelines.quick_salmon.scripts.salmon_quant import run_quick_salmon
            success = run_quick_salmon(
                sample_sheet=sample_sheet,
                index=index,
                output_dir=output,
                threads=threads,
                gene_map=gene_map,
                lib_type=lib_type,
                keep_quant=keep_quant
            )
        else:  # kallisto
            from raptor.pipelines.quick_kallisto.scripts.kallisto_quant import run_quick_kallisto
            success = run_quick_kallisto(
                sample_sheet=sample_sheet,
                index=index,
                output_dir=output,
                threads=threads,
                gene_map=gene_map,
                fragment_length=fragment_length,
                fragment_sd=fragment_sd,
                keep_quant=keep_quant
            )
        
        if success:
            click.echo("")
            click.echo("=" * 60)
            click.echo("✅ Quick-count completed successfully!")
            click.echo("=" * 60)
            click.echo("")
            click.echo("📂 Output files:")
            click.echo(f"   • {output}/counts.csv")
            click.echo(f"   • {output}/tpm.csv")
            click.echo(f"   • {output}/sample_info.csv")
            click.echo("")
            click.echo("📊 Next steps:")
            click.echo(f"   raptor profile --counts {output}/counts.csv --use-ml")
        else:
            click.echo("⚠️  Completed with some failures. Check sample_info.csv")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Create Sample Sheet Command
# =============================================================================

@main.command('create-sample-sheet')
@click.option('--input', '-i',
              type=click.Path(exists=True),
              required=True,
              help='Directory containing FASTQ files')
@click.option('--output', '-o',
              type=click.Path(),
              default='samples.csv',
              help='Output sample sheet file (default: samples.csv)')
@click.option('--paired/--single',
              default=True,
              help='Paired-end or single-end reads (default: paired)')
def create_sample_sheet_cmd(input, output, paired):
    """
    Create a sample sheet from FASTQ files.
    
    Scans a directory for FASTQ files and creates a sample sheet template.
    
    \b
    EXAMPLES:
        raptor create-sample-sheet --input fastq/ --output samples.csv
    """
    click.echo(f"🔍 Scanning {input} for FASTQ files...")
    
    try:
        import pandas as pd
        
        samples, read_type = auto_detect_samples(input)
        
        if not samples:
            click.echo("❌ No FASTQ files found!", err=True)
            sys.exit(1)
        
        click.echo(f"   Found {len(samples)} samples ({read_type})")
        
        # Create dataframe
        if read_type == 'paired':
            data = [{
                'sample_id': s.sample_id,
                'condition': '',
                'batch': '',
                'fastq_r1': s.fastq_r1,
                'fastq_r2': s.fastq_r2
            } for s in samples]
        else:
            data = [{
                'sample_id': s.sample_id,
                'condition': '',
                'batch': '',
                'fastq': s.fastq_r1
            } for s in samples]
        
        df = pd.DataFrame(data)
        df.to_csv(output, index=False)
        
        click.echo(f"✅ Sample sheet created: {output}")
        click.echo("\n📝 Edit the file to add condition and batch columns")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Profile Command (Data Profiling and Recommendation)
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
    based on data properties.
    
    \b
    EXAMPLES:
        raptor profile --counts counts.csv
        raptor profile --counts counts.csv --priority accuracy --format json
    """
    click.echo("🔍 Profiling RNA-seq data...")
    
    try:
        import pandas as pd
        counts_df = pd.read_csv(counts, index_col=0)
        
        metadata_df = None
        if metadata:
            metadata_df = pd.read_csv(metadata)
        
        click.echo(f"📊 Analyzing {counts_df.shape[0]} genes × {counts_df.shape[1]} samples...")
        
        profiler = RNAseqDataProfiler(counts_df, metadata_df)
        profile_data = profiler.run_full_profile()
        
        click.echo(f"🎯 Generating recommendation (priority: {priority})...")
        recommender = PipelineRecommender(profile_data)
        recommendation = recommender.get_recommendation(priority=priority)
        
        # Generate output
        if format == 'html':
            from raptor.report import generate_profile_report
            generate_profile_report(profile_data, recommendation, output)
            click.echo(f"✅ Report saved to: {output}")
        elif format == 'json':
            import json
            with open(output, 'w') as f:
                json.dump({'profile': profile_data, 'recommendation': recommendation}, f, indent=2)
            click.echo(f"✅ Profile saved to: {output}")
        else:
            _print_recommendation(recommendation)
        
        click.echo(f"\n🦖 Recommended: {recommendation['primary']['pipeline_name']}")
        click.echo(f"   Score: {recommendation['primary']['score']:.1f}/200")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Compare Command (Full Pipeline Benchmarking)
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
    their performance.
    
    \b
    EXAMPLES:
        raptor compare --data fastq/ --output results/
        raptor compare --data fastq/ --pipelines 1,3,5 --threads 16
    """
    click.echo("⚖️  Starting pipeline benchmarking...")
    
    try:
        if pipelines == 'all':
            pipeline_ids = list(range(1, 9))
        else:
            pipeline_ids = [int(p.strip()) for p in pipelines.split(',')]
        
        click.echo(f"📋 Running {len(pipeline_ids)} pipelines: {pipeline_ids}")
        click.echo(f"💻 Using {threads} threads and {memory} memory")
        
        benchmark = PipelineBenchmark(
            data_dir=data,
            output_dir=output,
            threads=threads,
            memory=memory,
            reference=reference
        )
        
        results = benchmark.run_pipelines(pipeline_ids)
        benchmark.save_results(results)
        
        click.echo(f"✅ Benchmarking complete! Results saved to: {output}")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Run Command (Run Specific Pipeline)
# =============================================================================

@main.command()
@click.option('--pipeline', '-p', type=int, required=True,
              help='Pipeline ID (1-8)')
@click.option('--data', '-d', type=click.Path(exists=True), required=True,
              help='Input data directory')
@click.option('--output', '-o', type=click.Path(), default='pipeline_results',
              help='Output directory')
@click.option('--threads', '-t', type=int, default=8,
              help='Number of threads')
@click.option('--reference', type=click.Path(exists=True),
              help='Reference genome/transcriptome')
def run(pipeline, data, output, threads, reference):
    """
    Run a specific pipeline on your data.
    
    \b
    EXAMPLES:
        raptor run --pipeline 3 --data fastq/ --output results/
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
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Report Command (Generate Comparison Report)
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
    
    \b
    EXAMPLES:
        raptor report --results benchmark_results/
        raptor report --results results/ --format pdf
    """
    click.echo("📊 Generating comparison report...")
    
    try:
        reporter = ReportGenerator(results)
        reporter.generate_report(output_file=output, format=format)
        click.echo(f"✅ Report generated: {output}")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Simulate Command (Generate Simulated Data)
# =============================================================================

@main.command()
@click.option('--n-genes', type=int, default=2000, help='Number of genes')
@click.option('--n-samples', type=int, default=6, help='Number of samples')
@click.option('--n-de', type=int, default=400, help='Number of DE genes')
@click.option('--fold-changes', default='0.5,2', help='Fold changes (comma-separated)')
@click.option('--output', '-o', type=click.Path(), default='simulated_data',
              help='Output directory')
def simulate(n_genes, n_samples, n_de, fold_changes, output):
    """
    Generate simulated RNA-seq data for testing.
    
    \b
    EXAMPLES:
        raptor simulate --output sim_data/
        raptor simulate --n-genes 5000 --n-samples 12 --n-de 1000
    """
    click.echo("🧬 Generating simulated RNA-seq data...")
    
    try:
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
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Demo Command
# =============================================================================

@main.command()
def demo():
    """
    Run a quick demo of RAPTOR functionality.
    
    Generates a small simulated dataset, profiles it, and shows recommendation.
    """
    click.echo("🎬 Running RAPTOR demo...")
    click.echo("=" * 60)
    
    try:
        import pandas as pd
        import numpy as np
        
        # Generate test data
        click.echo("\n1️⃣  Generating test data...")
        np.random.seed(42)
        n_genes, n_samples = 100, 6
        counts = pd.DataFrame(
            np.random.negative_binomial(10, 0.1, (n_genes, n_samples)),
            columns=[f'Sample_{i+1}' for i in range(n_samples)]
        )
        click.echo(f"   Created {n_genes} genes × {n_samples} samples")
        
        # Profile
        click.echo("\n2️⃣  Profiling data...")
        profiler = RNAseqDataProfiler(counts)
        profile_data = profiler.run_full_profile()
        
        click.echo(f"   Library size CV: {profile_data['library_stats']['cv']:.2f}")
        click.echo(f"   Zero inflation: {profile_data['count_distribution']['zero_pct']:.1f}%")
        click.echo(f"   BCV: {profile_data['biological_variation']['bcv']:.2f}")
        
        # Recommend
        click.echo("\n3️⃣  Getting recommendation...")
        recommender = PipelineRecommender(profile_data)
        recommendation = recommender.get_recommendation()
        
        click.echo(f"\n{'=' * 60}")
        click.echo(f"🎯 RECOMMENDED: {recommendation['primary']['pipeline_name']}")
        click.echo(f"   Score: {recommendation['primary']['score']:.1f}/200")
        
        click.echo(f"\n{'=' * 60}")
        click.echo("✅ Demo complete!")
        click.echo("\nNext steps:")
        click.echo("  • raptor quick-count --method salmon --sample-sheet samples.csv --index idx/ --output counts/")
        click.echo("  • raptor profile --counts counts.csv")
        
    except Exception as e:
        logger.error(f"Error: {e}", exc_info=True)
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# Helper Functions
# =============================================================================

def _print_recommendation(recommendation):
    """Print recommendation in text format."""
    rec = recommendation['primary']
    click.echo("\n" + "=" * 60)
    click.echo("🎯 RECOMMENDED PIPELINE")
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
