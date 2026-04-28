"""
RAPTOR Command-Line Interface v2.2.2 - COMPLETE PRODUCTION VERSION

Complete CLI with all commands for Modules 1-9.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
"""

import click
import sys
import raptor
import json
import logging
from pathlib import Path
from typing import Optional

# =============================================================================
# VALIDATION IMPORTS
# =============================================================================

try:
    from raptor.utils.validation import (
        validate_file_path,
        validate_directory_path,
        validate_positive_integer,
        validate_probability,
    )
    from raptor.utils.errors import (
        handle_errors,
        ValidationError,
        check_file_exists,
        validate_output_writable,
    )
except ImportError:
    # Fallback for development
    try:
        from utils.validation import (
            validate_file_path,
            validate_directory_path,
            validate_positive_integer,
            validate_probability,
        )
        from utils.errors import (
            handle_errors,
            ValidationError,
            check_file_exists,
            validate_output_writable,
        )
    except ImportError:
        # Last resort - create dummy functions
        def validate_file_path(p, **kwargs): return Path(p)
        def validate_directory_path(p, **kwargs): return Path(p)
        def validate_positive_integer(v, name): 
            if v < 1: raise ValueError(f"{name} must be positive")
        def validate_probability(v, name): 
            if not 0 <= v <= 1: raise ValueError(f"{name} must be 0-1")
        def check_file_exists(p, msg): 
            if not Path(p).exists(): raise FileNotFoundError(msg)
        def validate_output_writable(p): pass
        class ValidationError(ValueError): pass

def validate_cli_file(filepath, description="File"):
    """Validate file exists and is readable."""
    try:
        path = validate_file_path(filepath, must_exist=True)
        return path
    except FileNotFoundError:
        click.echo(f"❌ {description} not found: {filepath}", err=True)
        click.echo(f"   Current directory: {Path.cwd()}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error with {description}: {e}", err=True)
        sys.exit(1)

def validate_cli_directory(dirpath, create=True, description="Directory"):
    """Validate/create output directory."""
    try:
        path = validate_directory_path(dirpath, create_if_missing=create)
        return path
    except Exception as e:
        click.echo(f"❌ Cannot access {description}: {e}", err=True)
        sys.exit(1)

logger = logging.getLogger(__name__)

# =============================================================================
# MAIN CLI GROUP
# =============================================================================

@click.group()
# Version read dynamically — never hardcode here
@click.version_option(version=raptor.__version__, prog_name='RAPTOR')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.option('--quiet', '-q', is_flag=True, help='Suppress non-error output')
def main(verbose, quiet):
    """
    🦖 RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource
    
    A comprehensive framework for RNA-seq quality control, analysis,
    and pipeline optimization with ML-based recommendations.
    
    \b
    WORKFLOW (v2.2.2 - Complete Architecture):
    ══════════════════════════════════════════════════════════════
    
    STAGE 1: Fast Profiling (M1-M4)
    ─────────────────────────────────
    raptor quick-count    → M1: Generate quick counts (Salmon/Kallisto)
    raptor qc             → M2: Quality assessment & outlier detection
    raptor profile        → M3: Data profiling (32 features)
    raptor recommend      → M4: Get pipeline recommendation
    
    STAGE 2: Production Pipeline (M5)
    ─────────────────────────────────
    raptor pipeline list  → List available pipelines
    raptor pipeline run   → Run production pipeline
    
    STAGE 3: DE Analysis (M6-M9)
    ─────────────────────────────────
    Rscript run_deseq2.R  → M6: DE Analysis (external R script)
    raptor import-de      → M7: Import DE results
    raptor compare-de     → M7: Compare DE methods
    raptor optimize       → M8: Parameter optimization
    raptor ensemble       → M9: Ensemble analysis
    
    STAGE 4: Biomarker Discovery (M10)
    ─────────────────────────────────
    raptor biomarker          → M10: Discover biomarker panel
    raptor biomarker-survival → M10D: Survival biomarkers
    raptor biomarker-validate → M10: Validate on independent cohort
    
    ══════════════════════════════════════════════════════════════
    
    \b
    QUICK START:
        raptor init my_project
        raptor quick-count -m salmon -s samples.csv -i index/
        raptor qc --counts results/quick_counts/quick_gene_counts.csv
        raptor profile --counts results/quick_counts/quick_gene_counts.csv
        raptor recommend
        raptor import-de --input de_results/ --method deseq2
        raptor optimize --de-result de_result.pkl --ground-truth truth.csv
        raptor ensemble --methods fisher --deseq2 d.pkl --edger e.pkl
        raptor biomarker -c counts.csv -m metadata.csv -g condition
    
    For more information: https://github.com/AyehBlk/RAPTOR
    """
    if verbose:
        logging.basicConfig(level=logging.DEBUG,
                          format='%(asctime)s - %(levelname)s - %(message)s')
    elif quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO,
                          format='%(asctime)s - %(levelname)s - %(message)s')


# =============================================================================
# UTILITY COMMANDS
# =============================================================================

@main.command('init')
@click.argument('project_name')
@click.option('--template', type=click.Choice(['basic', 'full']), default='basic',
              help='Project template (basic or full)')
def init(project_name, template):
    """
    Initialize a new RAPTOR project directory.
    
    Creates a standardized directory structure with configuration files,
    sample templates, and documentation.
    
    \b
    CREATES:
        project_name/
        ├── config/
        │   ├── raptor_config.yaml
        │   └── sample_sheet_template.csv
        ├── data/
        │   └── fastq/
        ├── results/
        │   ├── quick_counts/
        │   ├── qc/
        │   ├── profile/
        │   └── de_analysis/
        ├── scripts/
        └── README.md
    
    \b
    EXAMPLES:
        raptor init my_rnaseq_project
        raptor init my_project --template full
    """
    try:
        project_path = Path(project_name)
        
        if project_path.exists():
            click.echo(f"❌ Project '{project_name}' already exists!", err=True)
            sys.exit(1)
        
        click.echo(f"🦖 Initializing RAPTOR project: {project_name}")
        click.echo()
        
        # Create directory structure
        dirs = [
            'config',
            'data/fastq',
            'results/quick_counts',
            'results/qc',
            'results/profile',
            'results/de_analysis',
            'scripts',
        ]
        
        if template == 'full':
            dirs.extend([
                'results/pipeline',
                'results/optimization',
                'results/ensemble',
            ])
        
        for d in dirs:
            (project_path / d).mkdir(parents=True, exist_ok=True)
            click.echo(f"   ✓ Created {d}/")
        
        # Create sample sheet template
        sample_sheet = project_path / 'config' / 'sample_sheet_template.csv'
        sample_sheet.write_text(
            "sample_id,condition,batch,fastq_r1,fastq_r2\n"
            "Sample1,Control,1,data/fastq/Sample1_R1.fastq.gz,data/fastq/Sample1_R2.fastq.gz\n"
            "Sample2,Control,1,data/fastq/Sample2_R1.fastq.gz,data/fastq/Sample2_R2.fastq.gz\n"
            "Sample3,Treatment,1,data/fastq/Sample3_R1.fastq.gz,data/fastq/Sample3_R2.fastq.gz\n"
            "Sample4,Treatment,1,data/fastq/Sample4_R1.fastq.gz,data/fastq/Sample4_R2.fastq.gz\n"
        )
        click.echo(f"   ✓ Created sample_sheet_template.csv")
        
        # Create README
        readme = project_path / 'README.md'
        readme.write_text(f"""# {project_name}

RAPTOR v2.2.2 RNA-seq Analysis Project

## Quick Start

1. Edit sample sheet:
   ```bash
   nano config/sample_sheet_template.csv
   ```

2. Run RAPTOR workflow:
   ```bash
   cd {project_name}
   
   # Stage 1: Fast profiling
   raptor quick-count -s config/sample_sheet_template.csv -i /path/to/index
   raptor qc --counts results/quick_counts/quick_gene_counts.csv
   raptor profile --counts results/quick_counts/quick_gene_counts.csv
   raptor recommend
   
   # Stage 2: Production pipeline
   raptor pipeline run -s config/sample_sheet_template.csv
   
   # Stage 3: DE analysis
   raptor import-de --input de_results/ --method deseq2
   raptor optimize --de-result results/de_analysis/de_result.pkl
   ```

## Documentation

- RAPTOR GitHub: https://github.com/AyehBlk/RAPTOR
- Documentation: https://raptor.readthedocs.io

## Project Structure

- `config/`: Configuration files and sample sheets
- `data/`: Raw data (FASTQ files)
- `results/`: Analysis outputs
- `scripts/`: Custom analysis scripts
""")
        click.echo(f"   ✓ Created README.md")
        
        click.echo()
        click.echo(f"✅ Project initialized successfully!")
        click.echo()
        click.echo(f"📁 Next steps:")
        click.echo(f"   1. cd {project_name}")
        click.echo(f"   2. Edit config/sample_sheet_template.csv")
        click.echo(f"   3. Run: raptor quick-count -s config/sample_sheet_template.csv")
        
    except Exception as e:
        click.echo(f"❌ Error initializing project: {e}", err=True)
        sys.exit(1)


@main.command('validate-installation')
def validate_installation():
    """
    Validate RAPTOR installation and check module availability.
    
    Checks:
    - All core modules are importable
    - All CLI commands are registered
    - Dependencies are installed
    - Version compatibility
    """
    try:
        click.echo("🦖 RAPTOR Installation Validation")
        click.echo("=" * 60)
        click.echo()
        
        # Import RAPTOR
        try:
            import raptor
            click.echo(f"✓ RAPTOR v{raptor.__version__} imported successfully")
        except ImportError as e:
            click.echo(f"❌ Cannot import RAPTOR: {e}", err=True)
            sys.exit(1)
        
        # Check module availability
        click.echo()
        click.echo("📦 Module Availability:")
        modules = raptor.get_available_modules()
        
        for module_name, available in modules.items():
            status = "✓" if available else "✗"
            click.echo(f"   {status} {module_name}")
        
        # Validate exports
        click.echo()
        click.echo("🔍 Validating Exports:")
        report = raptor.validate_installation()
        
        if report['exports_valid']:
            click.echo(f"   ✓ All {len(raptor.__all__)} exports valid")
        else:
            click.echo(f"   ❌ {len(report['missing_exports'])} missing exports:")
            for export in report['missing_exports'][:5]:
                click.echo(f"      - {export}")
        
        # Check for issues
        click.echo()
        if report['issues']:
            click.echo("⚠️  Issues Found:")
            for issue in report['issues']:
                click.echo(f"   - {issue}")
        else:
            click.echo("✅ No issues found!")
        
        click.echo()
        click.echo("=" * 60)
        
        if report['exports_valid'] and not report['issues']:
            click.echo("✅ Installation is valid and ready to use!")
            return 0
        else:
            click.echo("⚠️  Some issues detected - please review above")
            return 1
            
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)


# =============================================================================
# MODULE 1: QUICK COUNT
# =============================================================================


@main.command('quick-count')
@click.option('--method', '-m', type=click.Choice(['salmon', 'kallisto']), 
              default='salmon', help='Quantification method (default: salmon)')
@click.option('--samples', '-s', type=click.Path(exists=True),
              help='Sample sheet CSV file')
@click.option('--fastq-dir', '-f', type=click.Path(exists=True),
              help='Directory containing FASTQ files (auto-detect samples)')
@click.option('--index', '-i', type=click.Path(exists=True), required=True,
              help='Salmon index directory or Kallisto index file (.idx)')
@click.option('--tx2gene', '-g', type=click.Path(exists=True),
              help='Transcript-to-gene mapping file (TSV/CSV/GTF)')
@click.option('--output', '-o', type=click.Path(), default='results/quick_counts',
              help='Output directory (default: results/quick_counts)')
@click.option('--threads', '-t', type=int, default=8,
              help='Number of threads per sample (default: 8)')
@click.option('--lib-type', '-l', default='A',
              help='Salmon library type (default: A = auto-detect)')
@click.option('--fragment-length', type=int, default=200,
              help='Fragment length for single-end reads (default: 200)')
@click.option('--fragment-sd', type=int, default=20,
              help='Fragment length SD for single-end reads (default: 20)')
@click.option('--bootstraps', '-b', type=int, default=0,
              help='Bootstrap samples for Kallisto (default: 0 = disabled)')
@click.option('--keep-quant', is_flag=True,
              help='Keep individual sample quant directories')
@click.option('--no-gc-bias', is_flag=True,
              help='Disable GC bias correction (Salmon only)')
@click.option('--no-validate-mappings', is_flag=True,
              help='Disable mapping validation (Salmon only)')
def quick_count(method, samples, fastq_dir, index, tx2gene, output, threads,
                lib_type, fragment_length, fragment_sd, bootstraps, keep_quant,
                no_gc_bias, no_validate_mappings):
    """
    Module 1: Quick quantification with Salmon or Kallisto.
    
    Generates fast pseudo-alignment counts for QC and profiling (M2-M4).
    NOT for final DE analysis - use Module 5 production pipelines instead.
    
    \b
    WORKFLOW:
        FASTQ → Salmon/Kallisto → Gene counts → M2 QC → M3 Profile → M4 Recommend
    
    \b
    REQUIRES ONE OF:
        --samples  (-s) : Sample sheet CSV with columns: sample_id, condition, batch, fastq_r1, fastq_r2
        --fastq-dir (-f): Directory to auto-detect FASTQ files
    
    \b
    OUTPUT FILES (in results/quick_counts/):
        • quick_gene_counts.csv    - Gene-level count matrix (for M2-M4)
        • quick_tpm.csv            - TPM normalized counts
        • sample_info.csv          - Sample metadata & QC metrics
        • qc_report.html           - Quality control report
    
    \b
    EXAMPLES:
        # Salmon with sample sheet (most common)
        raptor quick-count -m salmon -s samples.csv -i salmon_index/
        
        # Kallisto with gene mapping
        raptor quick-count -m kallisto -s samples.csv -i kallisto.idx -g tx2gene.csv
        
        # Auto-detect samples from FASTQ directory
        raptor quick-count -m salmon -f data/fastq/ -i salmon_index/
        
        # Single-end reads (must specify fragment length)
        raptor quick-count -s samples.csv -i index/ --fragment-length 200 --fragment-sd 20
        
        # Custom settings
        raptor quick-count -s samples.csv -i index/ -t 16 --keep-quant
    
    \b
    NEXT STEPS:
        raptor qc --counts results/quick_counts/quick_gene_counts.csv
        raptor profile --counts results/quick_counts/quick_gene_counts.csv
        raptor recommend
    """
    try:
        click.echo("\n🦖 RAPTOR v2.2.2 - Module 1: Quick Quantification")
        click.echo("=" * 70)
        
        # =====================================================================
        # STEP 1: VALIDATION
        # =====================================================================
        
        # Must have either samples or fastq-dir
        if not samples and not fastq_dir:
            click.echo("\n❌ ERROR: Must provide either --samples or --fastq-dir", err=True)
            click.echo("\nExamples:", err=True)
            click.echo("  raptor quick-count -s samples.csv -i index/", err=True)
            click.echo("  raptor quick-count -f data/fastq/ -i index/", err=True)
            sys.exit(1)
        
        # Validate index exists
        index_path = Path(index)
        if not index_path.exists():
            click.echo(f"\n❌ ERROR: Index not found: {index}", err=True)
            click.echo(f"   Current directory: {Path.cwd()}", err=True)
            sys.exit(1)
        
        # Validate output directory
        output_path = validate_cli_directory(output, create=True, description="Output directory")
        
        # Validate sample sheet if provided
        if samples:
            sample_sheet_path = validate_cli_file(samples, "Sample sheet")
        else:
            # Auto-detect samples from FASTQ directory
            fastq_path = validate_cli_directory(fastq_dir, create=False, description="FASTQ directory")
            
            click.echo("\n📂 Auto-detecting samples from FASTQ directory...")
            
            # Import auto-detection utility
            try:
                # Try importing from quick_salmon first
                from raptor.pipelines.quick_salmon.detect_samples import auto_detect_samples
            except ImportError:
                try:
                    # Fall back to quick_kallisto
                    from raptor.pipelines.quick_kallisto.detect_samples import auto_detect_samples
                except ImportError:
                    click.echo("❌ ERROR: Could not import sample detection utility", err=True)
                    click.echo("   Please provide a sample sheet with --samples instead", err=True)
                    sys.exit(1)
            
            # Auto-detect samples
            detected_samples, read_type = auto_detect_samples(str(fastq_path))
            
            if not detected_samples:
                click.echo("❌ ERROR: No FASTQ files found in directory", err=True)
                sys.exit(1)
            
            # Create sample sheet
            import pandas as pd
            sample_data = [{
                'sample_id': s.sample_id,
                'condition': '',
                'batch': '',
                'fastq_r1': s.fastq_r1,
                'fastq_r2': s.fastq_r2 or ''
            } for s in detected_samples]
            
            sample_sheet_path = output_path / 'auto_detected_samples.csv'
            pd.DataFrame(sample_data).to_csv(sample_sheet_path, index=False)
            
            click.echo(f"✓ Detected {len(detected_samples)} samples ({read_type}-end)")
            click.echo(f"✓ Sample sheet created: {sample_sheet_path}")
        
        # Validate tx2gene if provided
        if tx2gene:
            tx2gene_path = validate_cli_file(tx2gene, "Transcript-to-gene mapping")
        else:
            tx2gene_path = None
        
        # Validate threads
        if threads < 1 or threads > 128:
            click.echo(f"❌ ERROR: Threads must be between 1-128 (got {threads})", err=True)
            sys.exit(1)
        
        # =====================================================================
        # STEP 2: DISPLAY CONFIGURATION
        # =====================================================================
        
        click.echo(f"\n📋 Configuration:")
        click.echo(f"   Method:       {method.upper()}")
        click.echo(f"   Index:        {index_path}")
        click.echo(f"   Sample sheet: {sample_sheet_path}")
        if tx2gene_path:
            click.echo(f"   Gene mapping: {tx2gene_path}")
        click.echo(f"   Output:       {output_path}")
        click.echo(f"   Threads:      {threads}")
        
        if method == 'salmon':
            click.echo(f"   Library type: {lib_type}")
            click.echo(f"   GC bias:      {'Disabled' if no_gc_bias else 'Enabled'}")
            click.echo(f"   Validation:   {'Disabled' if no_validate_mappings else 'Enabled'}")
            click.echo(f"   Fragment:     {fragment_length} ± {fragment_sd} bp (for single-end)")
        elif method == 'kallisto':
            click.echo(f"   Bootstraps:   {bootstraps}")
            click.echo(f"   Fragment:     {fragment_length} ± {fragment_sd} bp (for single-end)")
            click.echo(f"   Keep quant:   {keep_quant}")
        
        # =====================================================================
        # STEP 3: IMPORT AND INITIALIZE PIPELINE
        # =====================================================================
        
        click.echo(f"\n🚀 Initializing {method.upper()} pipeline...")
        
        if method == 'salmon':
            # Import Salmon pipeline
            try:
                from raptor.pipelines.quick_salmon import QuickSalmonPipeline
            except ImportError as e:
                click.echo(f"\n❌ ERROR: Could not import Salmon pipeline", err=True)
                click.echo(f"   {e}", err=True)
                click.echo("\n💡 Solutions:", err=True)
                click.echo("   1. Install RAPTOR: pip install -e .", err=True)
                click.echo("   2. Check import path is correct", err=True)
                click.echo("   3. Verify quick_salmon module exists", err=True)
                sys.exit(1)
            
            # Initialize Salmon pipeline with exact parameter names
            try:
                pipeline = QuickSalmonPipeline(
                    sample_sheet=str(sample_sheet_path),
                    index=str(index_path),
                    output_dir=str(output_path),
                    threads=threads,
                    gene_map=str(tx2gene_path) if tx2gene_path else None,
                    lib_type=lib_type,
                    validate_mappings=not no_validate_mappings,
                    gc_bias=not no_gc_bias,
                    seq_bias=False,  # Not exposed in CLI (advanced)
                    fragment_length_mean=fragment_length,
                    fragment_length_sd=fragment_sd
                )
                click.echo("✓ Salmon pipeline initialized successfully")
                
            except ValidationError as e:
                click.echo(f"\n❌ VALIDATION ERROR:", err=True)
                click.echo(f"{e}", err=True)
                sys.exit(1)
            except Exception as e:
                click.echo(f"\n❌ ERROR: Pipeline initialization failed", err=True)
                click.echo(f"   {e}", err=True)
                if logger.isEnabledFor(logging.DEBUG):
                    import traceback
                    traceback.print_exc()
                sys.exit(1)
            
        elif method == 'kallisto':
            # Import Kallisto pipeline
            try:
                from raptor.pipelines.quick_kallisto import QuickKallistoPipeline
            except ImportError as e:
                click.echo(f"\n❌ ERROR: Could not import Kallisto pipeline", err=True)
                click.echo(f"   {e}", err=True)
                click.echo("\n💡 Solutions:", err=True)
                click.echo("   1. Install RAPTOR: pip install -e .", err=True)
                click.echo("   2. Check import path is correct", err=True)
                click.echo("   3. Verify quick_kallisto module exists", err=True)
                sys.exit(1)
            
            # Initialize Kallisto pipeline with exact parameter names
            try:
                pipeline = QuickKallistoPipeline(
                    sample_sheet=str(sample_sheet_path),
                    index=str(index_path),
                    output_dir=str(output_path),
                    threads=threads,
                    gene_map=str(tx2gene_path) if tx2gene_path else None,
                    config={
                        'fragment_length': fragment_length,
                        'fragment_sd': fragment_sd,
                        'num_bootstraps': bootstraps,
                        'keep_quant': keep_quant
                    }
                )
                click.echo("✓ Kallisto pipeline initialized successfully")
                
            except ValidationError as e:
                click.echo(f"\n❌ VALIDATION ERROR:", err=True)
                click.echo(f"{e}", err=True)
                sys.exit(1)
            except Exception as e:
                click.echo(f"\n❌ ERROR: Pipeline initialization failed", err=True)
                click.echo(f"   {e}", err=True)
                if logger.isEnabledFor(logging.DEBUG):
                    import traceback
                    traceback.print_exc()
                sys.exit(1)
        
        # =====================================================================
        # STEP 4: RUN QUANTIFICATION
        # =====================================================================
        
        click.echo(f"\n▶️  Running quantification...")
        click.echo("   This may take several minutes depending on sample count")
        click.echo()
        
        try:
            # Run the pipeline
            result = pipeline.run()
            
            click.echo("\n" + "=" * 70)
            click.echo("✅ QUANTIFICATION COMPLETE!")
            click.echo("=" * 70)
            
        except KeyboardInterrupt:
            click.echo("\n\n⚠️  Interrupted by user", err=True)
            click.echo("   Partial results may be available in output directory", err=True)
            sys.exit(130)
        except Exception as e:
            click.echo(f"\n❌ ERROR: Pipeline execution failed", err=True)
            click.echo(f"   {e}", err=True)
            click.echo(f"\n📝 Check log file: {output_path}/quick_{method}.log", err=True)
            if logger.isEnabledFor(logging.DEBUG):
                import traceback
                traceback.print_exc()
            sys.exit(1)
        
        # =====================================================================
        # STEP 5: DISPLAY RESULTS & NEXT STEPS
        # =====================================================================
        
        click.echo(f"\n📂 Output Directory: {output_path}")
        click.echo(f"\n📊 Output Files:")
        click.echo(f"   • quick_gene_counts.csv   - Gene-level count matrix")
        click.echo(f"   • quick_tpm.csv           - TPM normalized counts")
        click.echo(f"   • sample_info.csv         - Sample metadata")
        
        # Check if files exist and show quick stats
        counts_file = output_path / 'quick_gene_counts.csv'
        if counts_file.exists():
            try:
                import pandas as pd
                counts = pd.read_csv(counts_file, index_col=0)
                click.echo(f"\n📈 Quick Stats:")
                click.echo(f"   Genes:        {len(counts):,}")
                click.echo(f"   Samples:      {len(counts.columns)}")
                click.echo(f"   Total counts: {counts.sum().sum():,.0f}")
                
                # Show library sizes
                lib_sizes = counts.sum(axis=0)
                click.echo(f"   Mean library: {lib_sizes.mean():,.0f}")
                click.echo(f"   Library CV:   {(lib_sizes.std() / lib_sizes.mean() * 100):.1f}%")
            except Exception as e:
                # Don't fail if we can't read stats
                click.echo(f"   (Could not read stats: {e})")
        
        # Check for QC report
        qc_report = output_path / 'qc_report.html'
        if qc_report.exists():
            click.echo(f"   • qc_report.html          - Quality control report")
            click.echo(f"\n📚 View QC report: open {qc_report}")
        
        # =====================================================================
        # NEXT STEPS
        # =====================================================================
        
        click.echo(f"\n🔜 Next Steps (RAPTOR Workflow):")
        click.echo(f"\n   Module 2 - Quality Assessment & Outlier Detection:")
        click.echo(f"   raptor qc --counts {output_path}/quick_gene_counts.csv")
        click.echo(f"\n   Module 3 - Data Profiling (32 features):")
        click.echo(f"   raptor profile --counts {output_path}/quick_gene_counts.csv")
        click.echo(f"\n   Module 4 - Pipeline Recommendation:")
        click.echo(f"   raptor recommend")
        
        click.echo("\n" + "=" * 70)
        click.echo("🦖 Module 1 complete! Continue to Module 2 for QC analysis.")
        click.echo("=" * 70 + "\n")
        
    except ValidationError as e:
        click.echo(f"\n❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except KeyboardInterrupt:
        click.echo("\n\n⚠️  Interrupted by user", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"\n❌ Unexpected Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)

# =============================================================================
# MODULE 2: QUALITY ASSESSMENT
# =============================================================================

@main.command('qc')
@click.option('--counts', '-c', type=click.Path(exists=True), required=True,
              help='Count matrix CSV file (genes x samples)')
@click.option('--metadata', '-m', type=click.Path(exists=True),
              help='Sample metadata CSV (enables batch effect detection)')
@click.option('--output', '-o', type=click.Path(), default='results/qc',
              help='Output directory (default: results/qc)')
@click.option('--normalization', '-n', 
              type=click.Choice(['log2', 'cpm', 'quantile', 'none'], case_sensitive=False),
              default='log2',
              help='Normalization method: log2 (default), cpm (varying lib sizes), quantile (aggressive), none (pre-normalized)')
@click.option('--consensus-threshold', '-t', type=int, default=3,
              help='Outlier consensus threshold: 1-6 methods must agree (default: 3)')
@click.option('--plot', is_flag=True,
              help='Generate QC visualization plots')
@click.option('--plot-output', type=click.Path(), default='qc_plots.pdf',
              help='Output file for plots (default: qc_plots.pdf)')
def qc(counts, metadata, output, normalization, consensus_threshold, plot, plot_output):
    """
    Comprehensive quality assessment and outlier detection (Module 2).
    
    Performs 6-component quality assessment with advanced outlier detection
    using consensus voting across 6 different methods. Identifies technical
    issues, batch effects, and outlier samples before DE analysis.
    
    \b
    QC COMPONENTS (6 total):
        • Library Quality (15%): Size distribution, CV, ranges
        • Gene Detection (20%): Zero inflation, detection rates
        • Outlier Detection (15%): 6-method consensus voting
        • Variance Structure (15%): PCA analysis, PC1/PC2 variance
        • Batch Effects (20%): Intelligent batch vs condition separation
        • Biological Signal (15%): Signal strength assessment
    
    \b
    OUTLIER DETECTION METHODS (6):
        1. PCA + Mahalanobis Distance
        2. Isolation Forest
        3. Local Outlier Factor (LOF)
        4. Elliptic Envelope
        5. Correlation-based
        6. Library Size Z-score
    
    \b
    OUTPUT FILES:
        results/qc/
        ├── quality_report.json      # Complete quality metrics
        ├── outliers.txt             # Detailed outlier analysis
        ├── qc_summary.txt           # Human-readable summary
        ├── qc_plots.pdf             # QC visualizations (if --plot)
        └── outlier_qc_plots.pdf     # Outlier plots (if --plot)
    
    \b
    NORMALIZATION OPTIONS:
        log2     → log2(counts + 1) - DEFAULT for similar library sizes
        cpm      → CPM + log2 - for varying library sizes (CV > 0.5)
        quantile → Quantile normalization - for visualization only!
        none     → No transformation - for pre-normalized data (TPM, FPKM, VST)
    
    \b
    EXAMPLES:
        # Basic QC (most common)
        raptor qc --counts results/quick_counts/quick_gene_counts.csv
        
        # With metadata for batch detection
        raptor qc -c counts.csv -m metadata.csv
        
        # CPM normalization for varying library sizes
        raptor qc -c counts.csv -n cpm
        
        # Conservative outlier detection (4 of 6 methods must agree)
        raptor qc -c counts.csv --consensus-threshold 4
        
        # With visualization plots
        raptor qc -c counts.csv -m metadata.csv --plot
        
        # Complete analysis
        raptor qc -c counts.csv -m metadata.csv -n cpm -t 3 --plot
    """
    try:
        click.echo("╔═══════════════════════════════════════════════════════════════╗")
        click.echo("║     🦖 RAPTOR v2.2.2 - Module 2: Quality Assessment          ║")
        click.echo("╚═══════════════════════════════════════════════════════════════╝")
        click.echo()
        
        # Validate inputs
        counts_path = validate_cli_file(counts, "Count matrix")
        output_path = validate_cli_directory(output, create=True)
        
        if metadata:
            metadata_path = validate_cli_file(metadata, "Metadata")
        else:
            metadata_path = None
        
        # Validate normalization
        normalization = normalization.lower()
        if normalization not in ['log2', 'cpm', 'quantile', 'none']:
            click.echo(f"❌ Invalid normalization: {normalization}", err=True)
            click.echo("   Valid options: log2, cpm, quantile, none", err=True)
            sys.exit(1)
        
        # Validate consensus threshold
        if not 1 <= consensus_threshold <= 6:
            click.echo(f"❌ Consensus threshold must be 1-6, got {consensus_threshold}", err=True)
            sys.exit(1)
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Count matrix: {counts_path}")
        if metadata_path:
            click.echo(f"   Metadata: {metadata_path}")
        else:
            click.echo(f"   Metadata: None (batch detection disabled)")
        click.echo(f"   Output: {output_path}")
        click.echo(f"   Normalization: {normalization}")
        click.echo(f"   Consensus threshold: {consensus_threshold}/6 methods")
        click.echo(f"   Generate plots: {plot}")
        click.echo()
        
        # Import QC module
        try:
            from raptor.quality_assessment import DataQualityAssessor
            import pandas as pd
            import numpy as np
        except ImportError as e:
            click.echo(f"❌ Cannot import RAPTOR modules: {e}", err=True)
            click.echo("   Ensure RAPTOR is installed: pip install -e .", err=True)
            click.echo("   Or check your PYTHONPATH", err=True)
            sys.exit(1)
        
        # Load data
        click.echo("📂 Loading data...")
        try:
            counts_df = pd.read_csv(counts_path, index_col=0)
        except Exception as e:
            click.echo(f"❌ Error loading count matrix: {e}", err=True)
            sys.exit(1)
        
        if metadata_path:
            try:
                metadata_df = pd.read_csv(metadata_path)
            except Exception as e:
                click.echo(f"❌ Error loading metadata: {e}", err=True)
                sys.exit(1)
        else:
            metadata_df = None
        
        click.echo(f"   ✓ Loaded {counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples")
        click.echo()
        
        # Run QC
        click.echo("🔬 Running quality assessment...")
        click.echo(f"   Normalization: {normalization}")
        
        try:
            assessor = DataQualityAssessor(
                counts_df, 
                metadata_df, 
                normalization=normalization
            )
            report = assessor.assess_quality()
        except Exception as e:
            click.echo(f"❌ Error during quality assessment: {e}", err=True)
            import traceback
            traceback.print_exc()
            sys.exit(1)
        
        click.echo(f"   ✓ 6-component assessment complete")
        click.echo()
        
        # Advanced outlier detection
        click.echo("🎯 Running advanced outlier detection...")
        click.echo(f"   Using 6 methods with consensus threshold: {consensus_threshold}")
        
        try:
            outlier_result = assessor.detect_outliers_advanced(
                consensus_threshold=consensus_threshold
            )
        except Exception as e:
            click.echo(f"❌ Error during outlier detection: {e}", err=True)
            import traceback
            traceback.print_exc()
            sys.exit(1)
        
        click.echo(f"   ✓ Outlier detection complete")
        click.echo()
        
        # Display results
        click.echo("╔═══════════════════════════════════════════════════════════════╗")
        click.echo("║                  Quality Assessment Results                   ║")
        click.echo("╚═══════════════════════════════════════════════════════════════╝")
        click.echo()
        
        # Print summary (already formatted)
        click.echo(report['summary'])
        click.echo()
        
        # Outliers
        n_outliers = outlier_result.n_outliers
        click.echo("─" * 70)
        click.echo("Outlier Detection Results:")
        click.echo("─" * 70)
        
        if n_outliers > 0:
            click.echo(f"⚠️  {n_outliers} outlier sample(s) detected ({outlier_result.outlier_percentage:.1f}%)")
            click.echo(f"   Consensus threshold: {consensus_threshold} methods")
            click.echo()
            click.echo("   Flagged samples:")
            for sample in outlier_result.outlier_samples:
                score = outlier_result.outlier_scores[sample]
                click.echo(f"     • {sample} (flagged by {score}/6 methods)")
            
            click.echo()
            click.echo("   Detection by method:")
            for method, samples in outlier_result.method_results.items():
                click.echo(f"     • {method:25s}: {len(samples)} outliers")
        else:
            click.echo("✅ No consensus outliers detected")
            click.echo(f"   (Threshold: {consensus_threshold}/6 methods must agree)")
        
        click.echo()
        
        # Batch effects
        if metadata_df is not None:
            batch_info = report['components']['batch_effects']
            
            if batch_info['batch_detected']:
                click.echo("─" * 70)
                click.echo("Batch Effect Analysis:")
                click.echo("─" * 70)
                click.echo("⚠️  Batch Effect Detected")
                click.echo(f"   Batch variable: {batch_info.get('batch_variable', 'Unknown')}")
                click.echo(f"   Batch F-statistic: {batch_info.get('batch_strength', 0):.2f}")
                
                if batch_info.get('condition_variable'):
                    click.echo(f"   Condition variable: {batch_info['condition_variable']}")
                    click.echo(f"   Condition F-statistic: {batch_info.get('condition_strength', 0):.2f}")
                
                if batch_info.get('confounded', False):
                    click.echo()
                    click.echo("   ⚠️  CRITICAL: CONFOUNDED DESIGN DETECTED")
                    click.echo("   Batch and condition are perfectly correlated!")
                    click.echo("   Cannot separate batch effect from biological signal.")
                
                click.echo()
                click.echo("   Recommendation:")
                for line in batch_info['recommendation'].split('\n'):
                    if line.strip():
                        click.echo(f"     {line}")
                click.echo()
        
        click.echo("═" * 70)
        click.echo()
        
        # Save results
        click.echo("💾 Saving results...")
        
        # Save quality report
        import json
        report_file = output_path / 'quality_report.json'
        
        # Helper function to convert numpy types to Python types
        def convert_types(obj):
            if isinstance(obj, dict):
                return {k: convert_types(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_types(v) for v in obj]
            elif hasattr(obj, 'item'):  # numpy scalar types
                return obj.item()
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                return obj
        
        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                json.dump(convert_types(report), f, indent=2)
            click.echo(f"   ✓ Quality report: {report_file}")
        except Exception as e:
            click.echo(f"   ⚠️  Could not save quality report: {e}")
        
        # Save outlier results
        outlier_file = output_path / 'outliers.txt'
        try:
            with open(outlier_file, 'w', encoding='utf-8') as f:
                f.write(outlier_result.summary())
            click.echo(f"   ✓ Outlier report: {outlier_file}")
        except Exception as e:
            click.echo(f"   ⚠️  Could not save outlier report: {e}")
        
        # Save summary
        summary_file = output_path / 'qc_summary.txt'
        try:
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write("╔═══════════════════════════════════════════════════════════════╗\n")
                f.write("║     🦖 RAPTOR Quality Assessment Summary                      ║\n")
                f.write("╚═══════════════════════════════════════════════════════════════╝\n\n")
                f.write(f"Generated: {pd.Timestamp.now()}\n\n")
                f.write("QUALITY ASSESSMENT:\n")
                f.write("=" * 70 + "\n")
                f.write(report['summary'])
                f.write("\n\n")
                f.write("OUTLIER DETECTION:\n")
                f.write("=" * 70 + "\n")
                f.write(outlier_result.summary())
                
                if metadata_df is not None and report['components']['batch_effects']['batch_detected']:
                    batch_info = report['components']['batch_effects']
                    f.write("\n\n")
                    f.write("BATCH EFFECT ANALYSIS:\n")
                    f.write("=" * 70 + "\n")
                    f.write(f"Batch detected: Yes\n")
                    f.write(f"Batch variable: {batch_info.get('batch_variable', 'Unknown')}\n")
                    f.write(f"Batch strength: {batch_info.get('batch_strength', 0):.2f}\n")
                    if batch_info.get('confounded', False):
                        f.write("\n⚠️  CONFOUNDED DESIGN - Cannot separate batch from condition\n")
                    f.write(f"\nRecommendation:\n{batch_info['recommendation']}\n")
            
            click.echo(f"   ✓ Summary: {summary_file}")
        except Exception as e:
            click.echo(f"   ⚠️  Could not save summary: {e}")
        
        # Generate plots if requested
        if plot:
            click.echo()
            click.echo("📊 Generating visualization plots...")
            try:
                plot_path = output_path / plot_output
                assessor.plot_quality_report(str(plot_path))
                click.echo(f"   ✓ QC plots: {plot_path}")
                
                outlier_plot_path = output_path / f"outlier_{plot_output}"
                assessor.plot_outliers_advanced(outlier_result, str(outlier_plot_path))
                click.echo(f"   ✓ Outlier plots: {outlier_plot_path}")
            except ImportError:
                click.echo(f"   ⚠️  Plotting requires matplotlib - install with: pip install matplotlib")
            except Exception as e:
                click.echo(f"   ⚠️  Could not generate plots: {e}")
        
        click.echo()
        
        # Recommendations
        overall_score = report['overall']['score']
        overall_status = report['overall']['status']
        
        click.echo("╔═══════════════════════════════════════════════════════════════╗")
        click.echo("║                     Recommendations                           ║")
        click.echo("╚═══════════════════════════════════════════════════════════════╝")
        click.echo()
        
        if overall_score >= 80:
            click.echo("   ✅ Excellent quality - Ready for analysis")
        elif overall_score >= 60:
            click.echo("   ✓ Good quality - Address flagged issues")
        elif overall_score >= 40:
            click.echo("   ⚠️  Acceptable quality - Carefully review all issues")
        else:
            click.echo("   ❌ Poor quality - Investigation required before proceeding")
        
        click.echo()
        
        if n_outliers > 0:
            click.echo(f"   → {n_outliers} outlier sample(s) detected:")
            click.echo("     • Investigate flagged samples")
            click.echo("     • Consider removing if technical failures")
            click.echo("     • Re-run QC after filtering")
            click.echo()
        
        if metadata_df is not None and report['components']['batch_effects']['batch_detected']:
            batch_info = report['components']['batch_effects']
            
            if batch_info.get('confounded', False):
                click.echo("   ⚠️  CRITICAL: Confounded design detected")
                click.echo("     • Cannot proceed with standard DE analysis")
                click.echo("     • Batch and condition perfectly correlated")
                click.echo("     • SOLUTION: Redesign experiment with balanced design")
                click.echo()
            else:
                batch_var = batch_info.get('batch_variable', 'batch')
                click.echo(f"   → Batch effect detected in '{batch_var}':")
                click.echo(f"     • Include as covariate in DE model: ~{batch_var} + condition")
                click.echo("     • Use ComBat/removeBatchEffect for visualization only")
                click.echo("     • Do NOT use batch-corrected data for DE analysis")
                click.echo()
        
        click.echo("📝 Next Steps:")
        if overall_score >= 60 and n_outliers == 0:
            click.echo("   1. ✓ Quality passed - Proceed to Module 3")
            click.echo(f"   2. raptor profile --counts {counts_path}")
        elif n_outliers > 0:
            click.echo("   1. Review outlier samples in results/qc/outliers.txt")
            click.echo("   2. Remove technical outliers from count matrix")
            click.echo("   3. Re-run: raptor qc --counts <filtered_counts.csv>")
            click.echo("   4. If quality improves, proceed to Module 3")
        else:
            click.echo("   1. Review detailed report in results/qc/")
            click.echo("   2. Address flagged quality issues")
            click.echo("   3. Re-run QC to verify improvements")
            click.echo("   4. Once quality ≥60, proceed to Module 3")
        
        click.echo()
        click.echo("═" * 70)
        click.echo()
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Unexpected Error: {e}", err=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)



# =============================================================================
# MODULE 3: DATA PROFILER (FIXED v2.2.2)
# =============================================================================

@main.command('profile')
@click.option('--counts', '-c', type=click.Path(exists=True), required=True,
              help='Count matrix CSV file (genes × samples)')
@click.option('--metadata', '-m', type=click.Path(exists=True),
              help='Sample metadata CSV file (optional)')
@click.option('--group-column', '-g', type=str, default='condition',
              help='Column name for experimental groups (default: condition)')
@click.option('--min-count', type=int, default=1,
              help='Minimum count threshold for expressed genes (default: 1)')
@click.option('--output', '-o', type=click.Path(), default='results/profile',
              help='Output directory (default: results/profile)')
@click.option('--verbose', '-v', is_flag=True,
              help='Print detailed profiling information')
def profile(counts, metadata, group_column, min_count, output, verbose):
    """
    Profile RNA-seq data characteristics (Module 3).
    
    Extracts 32 statistical features that characterize the dataset and
    are predictive of optimal DE pipeline performance.
    
    \b
    KEY FEATURES EXTRACTED:
        • BCV (Biological Coefficient of Variation) - CRITICAL
        • Sample size and group balance
        • Library size distribution
        • Dispersion estimates
        • Count distribution patterns
        • Sparsity and zero-inflation
        • Mean-variance relationship
    
    \b
    OUTPUT:
        results/profile/
        ├── data_profile.json           # Complete profile
        ├── profile_summary.txt         # Human-readable summary
        └── recommendation_features.json # Key features for Module 4
    
    \b
    EXAMPLES:
        # Basic profiling
        raptor profile --counts counts.csv
        
        # With metadata for group-based features
        raptor profile -c counts.csv -m metadata.csv -g treatment
        
        # Custom threshold and verbose output
        raptor profile -c counts.csv --min-count 5 --verbose
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 3: Data Profiler")
        click.echo()
        
        # Validate inputs
        counts_path = validate_cli_file(counts, "Count matrix")
        output_path = validate_cli_directory(output, create=True)
        
        if metadata:
            metadata_path = validate_cli_file(metadata, "Metadata")
        else:
            metadata_path = None
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Count matrix: {counts_path}")
        if metadata_path:
            click.echo(f"   Metadata: {metadata_path}")
            click.echo(f"   Group column: {group_column}")
        click.echo(f"   Min count threshold: {min_count}")
        click.echo(f"   Output: {output_path}")
        click.echo()
        
        # FIXED: Correct import path from raptor.profiler
        try:
            from raptor.profiler import RNAseqDataProfiler
            import pandas as pd
            import json
        except ImportError as e:
            click.echo(f"❌ Cannot import RAPTOR modules: {e}", err=True)
            click.echo("Make sure RAPTOR is properly installed:", err=True)
            click.echo("  pip install -e .", err=True)
            sys.exit(1)
        
        # Load data
        click.echo("📂 Loading data...")
        try:
            counts_df = pd.read_csv(counts_path, index_col=0)
        except Exception as e:
            click.echo(f"❌ Error loading count matrix: {e}", err=True)
            sys.exit(1)
        
        if metadata_path:
            try:
                metadata_df = pd.read_csv(metadata_path)
            except Exception as e:
                click.echo(f"❌ Error loading metadata: {e}", err=True)
                sys.exit(1)
        else:
            metadata_df = None
        
        click.echo(f"   ✓ Loaded {counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples")
        click.echo()
        
        # Run profiling
        click.echo("🔬 Profiling dataset characteristics...")
        
        try:
            # FIXED: Include min_count_threshold parameter
            profiler = RNAseqDataProfiler(
                counts_df, 
                metadata_df, 
                group_column=group_column,
                min_count_threshold=min_count
            )
            
            # FIXED: Correct method call - run_full_profile(), not profile()
            profile_result = profiler.run_full_profile()
            
        except Exception as e:
            click.echo(f"❌ Error during profiling: {e}", err=True)
            if verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
        
        # FIXED: Correct attribute access - DataProfile is a dataclass
        click.echo(f"   ✓ Extracted {len(profile_result.feature_vector)} features")
        click.echo()
        
        # Display key features
        if verbose:
            # Print full summary
            click.echo(profile_result.summary())
        else:
            # Print compact key metrics with FIXED attribute names
            click.echo("📊 Key Data Characteristics:")
            click.echo("=" * 70)
            
            click.echo()
            click.echo("Sample Characteristics:")
            click.echo(f"   Samples: {profile_result.n_samples}")
            click.echo(f"   Genes: {profile_result.n_genes:,}")
            click.echo(f"   Groups: {profile_result.n_groups}")
            click.echo(f"   Min group size: {profile_result.min_group_size}")
            click.echo(f"   Sample balance: {profile_result.sample_balance:.2f}")
            
            click.echo()
            click.echo("Library Size:")
            click.echo(f"   Mean: {profile_result.library_size_mean:,.0f} reads")
            click.echo(f"   Median: {profile_result.library_size_median:,.0f} reads")
            click.echo(f"   CV: {profile_result.library_size_cv:.3f}")
            click.echo(f"   Range: {profile_result.library_size_range:.1f}x")
            click.echo(f"   Category: {profile_result.library_size_category}")
            
            click.echo()
            click.echo("Dispersion (CRITICAL for Pipeline Selection):")
            click.echo(f"   Common φ: {profile_result.common_dispersion:.4f}")
            click.echo(f"   BCV: {profile_result.bcv:.3f} ({profile_result.bcv*100:.1f}% biological variation)")
            click.echo(f"   BCV Category: {profile_result.bcv_category}")
            click.echo(f"   Overdispersion: {profile_result.overdispersion_ratio:.2f}x")
            
            click.echo()
            click.echo("Expression & Detection:")
            click.echo(f"   Detected genes: {profile_result.n_expressed_genes:,} ({profile_result.detection_rate:.1%})")
            click.echo(f"   Reliable (≥10): {profile_result.n_reliably_expressed:,} ({profile_result.reliable_detection_rate:.1%})")
            click.echo(f"   Mean (log2): {profile_result.expression_mean:.2f}")
            
            click.echo()
            click.echo("Sparsity & Count Distribution:")
            click.echo(f"   Sparsity: {profile_result.sparsity:.1%}")
            click.echo(f"   Zero inflation: {profile_result.zero_inflation_index:.3f}")
            click.echo(f"   Low count proportion: {profile_result.low_count_proportion:.1%}")
            
            click.echo()
            click.echo("Mean-Variance Relationship:")
            click.echo(f"   Slope: {profile_result.mean_var_slope:.3f} (Poisson=1, NB≈1.5-2)")
            click.echo(f"   R²: {profile_result.mean_var_r_squared:.3f}")
            
        click.echo()
        click.echo("=" * 70)
        click.echo()
        
        # FIXED: Correct save implementation using DataProfile methods
        click.echo("💾 Saving profile...")
        
        # Save complete profile JSON
        profile_file = output_path / 'data_profile.json'
        with open(profile_file, 'w', encoding='utf-8') as f:
            f.write(profile_result.to_json())
        click.echo(f"   ✓ Profile JSON: {profile_file}")
        
        # Save human-readable summary
        summary_file = output_path / 'profile_summary.txt'
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("╔═══════════════════════════════════════════════════════════════╗\n")
            f.write("║         🦖 RAPTOR v2.2.2 - Data Profile Summary               ║\n")
            f.write("╚═══════════════════════════════════════════════════════════════╝\n\n")
            f.write(f"Generated: {pd.Timestamp.now()}\n")
            f.write(f"Count matrix: {counts_path}\n")
            if metadata_path:
                f.write(f"Metadata: {metadata_path}\n")
                f.write(f"Group column: {group_column}\n")
            f.write(f"Min count threshold: {min_count}\n\n")
            f.write(profile_result.summary())
        click.echo(f"   ✓ Summary text: {summary_file}")
        
        # Save recommendation features for Module 4
        rec_features = profile_result.get_recommendation_features()
        rec_file = output_path / 'recommendation_features.json'
        with open(rec_file, 'w', encoding='utf-8') as f:
            json.dump(rec_features, f, indent=2)
        click.echo(f"   ✓ Recommendation features: {rec_file}")
        
        click.echo()
        
        # Next steps
        click.echo("📝 Next Steps:")
        click.echo(f"   1. Review profile: {output_path}/")
        click.echo(f"   2. Key metric - BCV: {profile_result.bcv:.3f} ({profile_result.bcv_category})")
        click.echo(f"   3. Get pipeline recommendation:")
        click.echo(f"      raptor recommend --profile {profile_file}")
        click.echo()
        
        # Quick interpretation
        if profile_result.bcv < 0.2:
            click.echo("💡 Low BCV → Consider limma-voom or DESeq2")
        elif profile_result.bcv < 0.4:
            click.echo("💡 Moderate BCV → DESeq2, edgeR, or limma-voom work well")
        else:
            click.echo("💡 High BCV → edgeR or edgeR_robust recommended")
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Unexpected error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# MODULE 4: PIPELINE RECOMMENDER
# =============================================================================

@main.command('recommend')
@click.option('--profile', '-p', type=click.Path(exists=True),
              help='Data profile JSON (if not provided, will look in results/profile/)')
@click.option('--method', '-m', type=click.Choice(['rule-based', 'ml', 'both']),
              default='rule-based', help='Recommendation method (default: rule-based)')
@click.option('--output', '-o', type=click.Path(), default='results/recommendation',
              help='Output directory')
@click.option('--verbose-explanation', is_flag=True,
              help='Show detailed explanation of recommendation')
def recommend(profile, method, output, verbose_explanation):
    """
    Get pipeline recommendation based on data characteristics (Module 4).
    
    Uses either rule-based logic, ML model, or both to recommend the optimal
    DE analysis pipeline (DESeq2, edgeR, limma-voom, or Wilcoxon) based on
    your dataset's characteristics.
    
    \b
    RECOMMENDATION METHODS:
        rule-based : Literature-based decision tree (DEFAULT)
        ml         : Random Forest classifier (requires training)
        both       : Combine both methods for confidence
    
    \b
    RECOMMENDED PIPELINES:
        DESeq2      : General-purpose, good for n < 8 per group
        edgeR       : High dispersion, small samples
        limma-voom  : Large samples (n > 20), low dispersion
        Wilcoxon    : Non-parametric, n ≥ 8 per group
        edgeR_robust: Outliers + small samples
    
    \b
    OUTPUT:
        results/recommendation/
        ├── recommendation.json      # Recommended pipeline + confidence
        ├── explanation.txt          # Detailed reasoning
        └── summary.txt              # Formatted summary
    
    \b
    EXAMPLES:
        # Automatic (uses results/profile/data_profile.json)
        raptor recommend
        
        # Explicit profile path
        raptor recommend --profile my_profile.json
        
        # ML-based (if model trained)
        raptor recommend --method ml
        
        # Both methods for validation
        raptor recommend --method both
        
        # Detailed explanation
        raptor recommend --verbose-explanation
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 4: Pipeline Recommender")
        click.echo()
        
        # Determine profile path
        if not profile:
            profile_path = Path('results/profile/data_profile.json')
            if not profile_path.exists():
                click.echo("❌ No profile found! Run 'raptor profile' first", err=True)
                click.echo(f"   Looking for: {profile_path}", err=True)
                sys.exit(1)
        else:
            profile_path = validate_cli_file(profile, "Data profile")
        
        output_path = validate_cli_directory(output, create=True)
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Data profile: {profile_path}")
        click.echo(f"   Method: {method}")
        click.echo(f"   Output: {output_path}")
        click.echo()
        
        # FIXED: Import from correct modules
        try:
            from raptor.recommender import PipelineRecommender
            from raptor.profiler import DataProfile
        except ImportError as e:
            click.echo(f"❌ Cannot import RAPTOR modules: {e}", err=True)
            click.echo("   Make sure RAPTOR is installed: pip install -e .", err=True)
            sys.exit(1)
        
        # Load profile
        click.echo("📂 Loading data profile...")
        with open(profile_path) as f:
            profile_data = json.load(f)
        
        # FIXED: PipelineRecommender accepts dict directly, no need for from_dict()
        # FIXED: Check feature_vector length instead of .features attribute
        n_features = len(profile_data.get('feature_vector', []))
        if n_features == 0:
            # Fallback: count non-None values in profile
            n_features = sum(1 for k, v in profile_data.items() 
                           if v is not None and k not in ['timestamp', 'version'])
        
        click.echo(f"   ✓ Loaded profile with {n_features} features")
        click.echo()
        
        # Get recommendations
        recommendations = {}
        
        # Convert dict to object with attribute access for recommender
        from types import SimpleNamespace
        profile_obj = SimpleNamespace(**profile_data)
        
        if method in ['rule-based', 'both']:
            click.echo("🧠 Rule-based recommendation...")
            
            # FIXED: Pass profile object to __init__, then call get_recommendation()
            rule_recommender = PipelineRecommender(profile_obj)
            rule_rec = rule_recommender.get_recommendation()
            recommendations['rule-based'] = rule_rec
            
            # FIXED: Use correct attribute names
            click.echo(f"   ✓ {rule_rec.primary_pipeline} (score: {rule_rec.primary_score:.0f}%)")
            click.echo()
        
        if method in ['ml', 'both']:
            click.echo("🤖 ML-based recommendation...")
            try:
                # Try to load ML recommender
                from raptor.ml_recommender import MLPipelineRecommender
                
                # Look for model file
                model_path = Path('models/ml_recommender.pkl')
                if not model_path.exists():
                    raise FileNotFoundError(f"ML model not found: {model_path}")
                
                ml_recommender = MLPipelineRecommender.load(str(model_path))
                ml_rec = ml_recommender.predict(profile_obj)
                recommendations['ml'] = ml_rec
                
                # FIXED: Use correct attribute names
                click.echo(f"   ✓ {ml_rec.primary_pipeline} (score: {ml_rec.primary_score:.0f}%)")
            except (ImportError, FileNotFoundError) as e:
                click.echo(f"   ⚠️  ML model not available: {e}")
                if method == 'ml':
                    click.echo("   Cannot proceed with ML-only method", err=True)
                    sys.exit(1)
                click.echo("   Continuing with rule-based only...")
            click.echo()
        
        # Determine final recommendation
        if len(recommendations) == 0:
            click.echo("❌ No recommendations generated!", err=True)
            sys.exit(1)
        
        # Display final recommendation
        click.echo("=" * 70)
        click.echo("📊 PIPELINE RECOMMENDATION")
        click.echo("=" * 70)
        click.echo()
        
        if method == 'both' and len(recommendations) == 2:
            rule_rec = recommendations['rule-based']
            ml_rec = recommendations['ml']
            
            # FIXED: Use correct attribute names for comparison
            if rule_rec.primary_pipeline == ml_rec.primary_pipeline:
                click.echo(f"✅ CONSENSUS: {rule_rec.primary_pipeline}")
                click.echo(f"   Both methods agree!")
                avg_score = (rule_rec.primary_score + ml_rec.primary_score) / 2
                click.echo(f"   Average confidence: {avg_score:.0f}%")
                final_rec = rule_rec
            else:
                click.echo(f"⚠️  DISAGREEMENT:")
                click.echo(f"   Rule-based: {rule_rec.primary_pipeline} ({rule_rec.primary_score:.0f}%)")
                click.echo(f"   ML-based:   {ml_rec.primary_pipeline} ({ml_rec.primary_score:.0f}%)")
                click.echo()
                
                # Use higher confidence
                if rule_rec.primary_score >= ml_rec.primary_score:
                    final_rec = rule_rec
                    click.echo(f"   → Recommending: {final_rec.primary_pipeline} (higher confidence)")
                else:
                    final_rec = ml_rec
                    click.echo(f"   → Recommending: {final_rec.primary_pipeline} (higher confidence)")
        else:
            # Single method
            final_rec = list(recommendations.values())[0]
            click.echo(f"✅ Recommended: {final_rec.primary_pipeline}")
            click.echo(f"   Confidence: {final_rec.primary_score:.0f}%")
            click.echo(f"   Reason: {final_rec.primary_reason}")
        
        click.echo()
        
        # ADDED: Display full summary
        click.echo(final_rec.summary())
        click.echo()
        
        # Show detailed explanation if requested
        if verbose_explanation and final_rec.decision_factors:
            click.echo("=" * 70)
            click.echo("📖 DETAILED EXPLANATION")
            click.echo("=" * 70)
            click.echo()
            click.echo("Decision Factors:")
            for factor, value in final_rec.decision_factors.items():
                if isinstance(value, float):
                    click.echo(f"  • {factor}: {value:.3f}")
                else:
                    click.echo(f"  • {factor}: {value}")
            click.echo()
            
            if final_rec.warnings:
                click.echo("Warnings:")
                for warning in final_rec.warnings:
                    click.echo(f"  ⚠️  {warning}")
                click.echo()
        
        # Save results
        click.echo("=" * 70)
        click.echo("💾 Saving recommendation...")
        click.echo()
        
        # FIXED: Save using to_dict() method
        rec_file = output_path / 'recommendation.json'
        with open(rec_file, 'w', encoding='utf-8') as f:
            result = {
                'method': method,
                'profile_file': str(profile_path),
                **final_rec.to_dict()  # Use the built-in to_dict() method
            }
            json.dump(result, f, indent=2)
        
        click.echo(f"   ✓ recommendation.json")
        
        # Save formatted summary
        summary_file = output_path / 'summary.txt'
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write(final_rec.summary())
        click.echo(f"   ✓ summary.txt")
        
        # Save explanation if verbose
        if verbose_explanation:
            explanation_file = output_path / 'explanation.txt'
            with open(explanation_file, 'w', encoding='utf-8') as f:
                f.write("RAPTOR Module 4 - Pipeline Recommendation\n")
                f.write("=" * 70 + "\n\n")
                f.write(f"Profile: {profile_path}\n")
                f.write(f"Method: {method}\n\n")
                f.write(final_rec.summary())
                f.write("\n\nDecision Factors:\n")
                for factor, value in final_rec.decision_factors.items():
                    f.write(f"  {factor}: {value}\n")
                if final_rec.warnings:
                    f.write("\nWarnings:\n")
                    for warning in final_rec.warnings:
                        f.write(f"  - {warning}\n")
            click.echo(f"   ✓ explanation.txt")
        
        click.echo()
        click.echo(f"   All files saved to: {output_path}/")
        click.echo()
        
        # Next steps
        click.echo("=" * 70)
        click.echo("📝 Next Steps")
        click.echo("=" * 70)
        click.echo()
        click.echo(f"1. Run {final_rec.primary_pipeline} DE analysis in R")
        click.echo()
        click.echo("2. After DE analysis, import results:")
        
        # FIXED: Format pipeline name correctly for import command
        pipeline_cmd = final_rec.primary_pipeline.lower().replace('-', '').replace('_', '')
        click.echo(f"   raptor import-de --method {pipeline_cmd} --results de_results.csv")
        click.echo()
        
        # Show alternative
        click.echo(f"Alternative: {final_rec.alternative_pipeline} ({final_rec.alternative_score:.0f}%)")
        if abs(final_rec.primary_score - final_rec.alternative_score) < 10:
            click.echo("   ℹ️  Scores are close - consider running both for ensemble analysis")
        click.echo()
        
    except KeyboardInterrupt:
        click.echo("\n\n⚠️  Operation cancelled by user", err=True)
        sys.exit(1)
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# MODULE 5: PRODUCTION PIPELINES
# =============================================================================

@main.group('pipeline')
def pipeline():
    """
    Production pipeline commands (Module 5).
    
    Run full RNA-seq analysis pipelines from FASTQ to counts.
    These generate publication-quality results for DE analysis.
    """
    pass


@pipeline.command('list')
def pipeline_list():
    """
    List available production pipelines.
    
    Shows all registered pipelines with descriptions and requirements.
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Available Pipelines")
        click.echo("=" * 60)
        click.echo()
        
        pipelines_info = [
            {
                'name': 'salmon',
                'description': 'Salmon quasi-mapping → tximport',
                'speed': 'Fast (~10 min/sample)',
                'accuracy': 'Good',
                'use_case': 'Quick analysis, transcript-level DE',
            },
            {
                'name': 'kallisto',
                'description': 'Kallisto pseudo-alignment → tximport',
                'speed': 'Very fast (~5 min/sample)',
                'accuracy': 'Good',
                'use_case': 'Ultra-fast quantification',
            },
            {
                'name': 'star-featurecounts',
                'description': 'STAR alignment → featureCounts',
                'speed': 'Slow (~30 min/sample)',
                'accuracy': 'Best',
                'use_case': 'Publication-quality, splice-aware',
            },
            {
                'name': 'hisat2-htseq',
                'description': 'HISAT2 alignment → HTSeq',
                'speed': 'Medium (~20 min/sample)',
                'accuracy': 'Good',
                'use_case': 'Alternative to STAR',
            },
        ]
        
        for i, pipeline_info in enumerate(pipelines_info, 1):
            click.echo(f"{i}. {pipeline_info['name']}")
            click.echo(f"   {pipeline_info['description']}")
            click.echo(f"   Speed: {pipeline_info['speed']}")
            click.echo(f"   Accuracy: {pipeline_info['accuracy']}")
            click.echo(f"   Use case: {pipeline_info['use_case']}")
            click.echo()
        
        click.echo("Usage:")
        click.echo("   raptor pipeline run --name <pipeline_name> -s samples.csv")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@pipeline.command('run')
@click.option('--name', '-n', required=True,
              type=click.Choice(['salmon', 'kallisto', 'star-featurecounts', 'hisat2-htseq']),
              help='Pipeline name')
@click.option('--samples', '-s', type=click.Path(exists=True), required=True,
              help='Sample sheet CSV')
@click.option('--genome-index', '-g', type=click.Path(exists=True),
              help='Genome index directory')
@click.option('--annotation', '-a', type=click.Path(exists=True),
              help='GTF/GFF annotation file')
@click.option('--output', '-o', type=click.Path(), default='results/pipeline',
              help='Output directory')
@click.option('--threads', '-t', type=int, default=8,
              help='Number of threads')
def pipeline_run(name, samples, genome_index, annotation, output, threads):
    """
    Run production pipeline from FASTQ to counts.
    
    \b
    EXAMPLES:
        # Salmon pipeline
        raptor pipeline run -n salmon -s samples.csv -g salmon_index/
        
        # STAR + featureCounts
        raptor pipeline run -n star-featurecounts -s samples.csv \\
            -g star_index/ -a genes.gtf
    """
    try:
        click.echo(f"🦖 RAPTOR v2.2.2 - Pipeline: {name}")
        click.echo()
        
        # Validate inputs
        sample_path = validate_cli_file(samples, "Sample sheet")
        output_path = validate_cli_directory(output, create=True)
        
        if genome_index:
            index_path = validate_cli_directory(genome_index, create=False, description="Genome index")
        
        if annotation:
            annot_path = validate_cli_file(annotation, "Annotation")
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Pipeline: {name}")
        click.echo(f"   Sample sheet: {sample_path}")
        if genome_index:
            click.echo(f"   Genome index: {index_path}")
        if annotation:
            click.echo(f"   Annotation: {annot_path}")
        click.echo(f"   Output: {output_path}")
        click.echo(f"   Threads: {threads}")
        click.echo()
        
        click.echo("⚠️  NOTE: Production pipelines require external tools:")
        
        if name == 'salmon':
            click.echo("   → Salmon must be installed")
            click.echo("   → R with tximport package")
        elif name == 'kallisto':
            click.echo("   → Kallisto must be installed")
            click.echo("   → R with tximport package")
        elif name == 'star-featurecounts':
            click.echo("   → STAR aligner must be installed")
            click.echo("   → featureCounts (subread package)")
        elif name == 'hisat2-htseq':
            click.echo("   → HISAT2 must be installed")
            click.echo("   → HTSeq-count")
        
        click.echo()
        click.echo("💡 This is a placeholder - integrate with your pipeline wrappers")
        click.echo()
        click.echo("📝 Next steps:")
        click.echo("   1. Run pipeline (integrate with Snakemake/Nextflow)")
        click.echo("   2. raptor import-de --input pipeline_results/")
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


# =============================================================================
# MODULE 7: DE IMPORT COMMANDS
# =============================================================================

@main.command('import-de')
@click.option('--input', '-i', type=click.Path(exists=True), required=True,
              help='Input directory or file with DE results')
@click.option('--method', '-m', 
              type=click.Choice(['deseq2', 'edger', 'limma', 'wilcoxon']),
              required=True,
              help='DE analysis method')
@click.option('--output', '-o', type=click.Path(), default='results/de_imported',
              help='Output directory for imported results')
@click.option('--fdr-threshold', type=float, default=0.05,
              help='FDR threshold for significance (default: 0.05)')
@click.option('--lfc-threshold', type=float, default=0.0,
              help='Log fold-change threshold (default: 0.0)')
def import_de(input, method, output, fdr_threshold, lfc_threshold):
    """
    Import differential expression results from various tools (Module 7).
    
    \b
    SUPPORTED METHODS:
        deseq2   : DESeq2 results (CSV/TSV)
        edger    : edgeR results (CSV/TSV)
        limma    : limma-voom results (CSV/TSV)
        wilcoxon : Wilcoxon rank-sum test results
    
    \b
    INPUT FORMATS:
        DESeq2   : baseMean, log2FoldChange, pvalue, padj
        edgeR    : logFC, PValue, FDR
        limma    : logFC, P.Value, adj.P.Val
        Wilcoxon : log2FC, pvalue, padj
    
    \b
    OUTPUT:
        results/de_imported/
        ├── de_result.pkl         # Serialized DEResult object
        ├── de_summary.json       # Summary statistics
        ├── significant_genes.csv # Filtered genes
        └── volcano_plot.png      # Visualization
    
    \b
    EXAMPLES:
        raptor import-de -i deseq2_results.csv -m deseq2
        raptor import-de -i edger/ -m edger --fdr-threshold 0.01
        raptor import-de -i limma_results.txt -m limma
    """
    try:
        # Validate inputs
        input_path = validate_cli_file(input, "DE results file/directory")
        output_path = validate_cli_directory(output, create=True, description="Output directory")
        
        # Validate thresholds
        if not 0 <= fdr_threshold <= 1:
            raise ValidationError(f"FDR threshold must be between 0 and 1, got {fdr_threshold}")
        
        if lfc_threshold < 0:
            raise ValidationError(f"LFC threshold must be non-negative, got {lfc_threshold}")
        
        click.echo(f"\n📊 Importing {method.upper()} results...")
        click.echo(f"   Input: {input_path}")
        click.echo(f"   Output: {output_path}")
        click.echo(f"   FDR threshold: {fdr_threshold}")
        click.echo(f"   LFC threshold: {lfc_threshold}")
        click.echo()
        
        # Import DE import module
        try:
            from raptor.de_import import import_de_results
            import pandas as pd
        except ImportError as e:
            click.echo(f"❌ Module 7 not available: {e}", err=True)
            sys.exit(1)
        
        # Load and import
        click.echo("📂 Importing DE results...")
        de_result = import_de_results(
            de_file=input_path,
            output_dir=output_path,
            pipeline=method,
            fdr_threshold=fdr_threshold,
            lfc_threshold=lfc_threshold
        )
        
        click.echo(f"   ✓ Imported {de_result.n_genes} genes")
        click.echo(f"   ✓ {de_result.n_significant} significant")
        click.echo()
        
        # Display summary
        click.echo("📊 DE Analysis Summary:")
        click.echo("=" * 60)
        click.echo(f"Pipeline: {de_result.pipeline}")
        click.echo(f"Total genes: {de_result.n_genes}")
        click.echo(f"Significant (FDR < {fdr_threshold}): {de_result.n_significant}")
        click.echo(f"Upregulated: {de_result.n_up}")
        click.echo(f"Downregulated: {de_result.n_down}")
        click.echo()
        
        # Files already saved by import_de_results()
        click.echo("✅ Import complete!")
        click.echo(f"   Results saved to: {output_path}/")
        click.echo(f"   - de_standardized.csv ({de_result.n_genes:,} genes)")
        click.echo(f"   - de_significant.csv ({de_result.n_significant:,} genes)")
        click.echo(f"   - de_summary.json")
        click.echo(f"   - de_result.pkl (for M8-M10)")
        click.echo()
        
        # Next steps
        click.echo("📝 Next steps:")
        click.echo("   1. Review results in results/de_imported/")
        click.echo(f"   2. raptor optimize --de-result {output_path}/de_result.pkl")
        click.echo("   3. raptor ensemble --methods fisher ...")
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command('compare-de')
@click.argument('de_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), default='results/de_comparison',
              help='Output directory')
@click.option('--venn-diagram', is_flag=True, help='Generate Venn diagram')
def compare_de(de_files, output, venn_diagram):
    """
    Compare multiple DE results (Module 7).
    
    Analyzes overlap between different DE methods or parameter sets.
    Generates comparison tables, Venn diagrams, and concordance metrics.
    
    \b
    EXAMPLES:
        raptor compare-de deseq2.pkl edger.pkl limma.pkl
        raptor compare-de method1.pkl method2.pkl --venn-diagram
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 7: DE Comparison")
        click.echo()
        
        # Validate inputs
        output_path = validate_cli_directory(output, create=True)
        
        if len(de_files) < 2:
            click.echo("❌ Need at least 2 DE result files to compare", err=True)
            sys.exit(1)
        
        click.echo(f"📋 Comparing {len(de_files)} DE results:")
        for i, f in enumerate(de_files, 1):
            click.echo(f"   {i}. {f}")
        click.echo()
        
        # Import comparison module
        # FIXED: compare_de_results doesn't exist - implement comparison directly
        try:
            from raptor.de_import import DEResult
            import pickle
        except ImportError as e:
            click.echo(f"❌ Module 7 not available: {e}", err=True)
            sys.exit(1)
        
        # Load all results
        click.echo("📂 Loading DE results...")
        de_results = {}
        for fpath in de_files:
            with open(fpath, 'rb') as f:
                de_result = pickle.load(f)
            name = Path(fpath).stem
            de_results[name] = de_result
            click.echo(f"   ✓ {name}: {de_result.n_significant} significant genes")
        click.echo()
        
        # Compare
        # FIXED: Implement comparison logic directly
        click.echo("🔬 Comparing results...")
        
        # Get significant gene sets from each result
        gene_sets = {
            name: set(result.significant_genes)
            for name, result in de_results.items()
        }
        
        # Calculate overlap and union
        all_sets = list(gene_sets.values())
        if len(all_sets) >= 2:
            overlap = set.intersection(*all_sets)
            union = set.union(*all_sets)
        else:
            overlap = all_sets[0] if all_sets else set()
            union = all_sets[0] if all_sets else set()
        
        # Calculate metrics
        overlap_count = len(overlap)
        union_count = len(union)
        jaccard_index = overlap_count / union_count if union_count > 0 else 0.0
        
        # Build comparison dictionary
        comparison = {
            'n_methods': len(gene_sets),
            'overlap_count': overlap_count,
            'overlap_genes': sorted(list(overlap)),
            'union_count': union_count,
            'union_genes': sorted(list(union)),
            'jaccard_index': jaccard_index,
            'per_method': {
                name: {
                    'n_significant': len(genes),
                    'genes': sorted(list(genes))
                }
                for name, genes in gene_sets.items()
            },
            'pairwise_overlaps': {}
        }
        
        # Calculate pairwise overlaps for detailed comparison
        method_names = list(gene_sets.keys())
        for i, name1 in enumerate(method_names):
            for name2 in method_names[i+1:]:
                overlap_pair = gene_sets[name1] & gene_sets[name2]
                union_pair = gene_sets[name1] | gene_sets[name2]
                jaccard_pair = len(overlap_pair) / len(union_pair) if union_pair else 0.0
                
                comparison['pairwise_overlaps'][f"{name1}_vs_{name2}"] = {
                    'overlap': len(overlap_pair),
                    'jaccard_index': jaccard_pair
                }
        
        click.echo()
        click.echo("📊 Comparison Results:")
        click.echo("=" * 60)
        click.echo(f"Methods compared: {comparison['n_methods']}")
        click.echo(f"Overall overlap: {comparison['overlap_count']} genes")
        click.echo(f"Overall union: {comparison['union_count']} genes")
        click.echo(f"Overall Jaccard index: {comparison['jaccard_index']:.3f}")
        click.echo()
        
        # Per-method breakdown
        click.echo("Per-method significant genes:")
        for name, info in comparison['per_method'].items():
            click.echo(f"  • {name}: {info['n_significant']} genes")
        click.echo()
        
        # Pairwise comparisons (if more than 2 methods)
        if len(comparison['pairwise_overlaps']) > 0:
            click.echo("Pairwise comparisons:")
            for pair_name, pair_info in comparison['pairwise_overlaps'].items():
                click.echo(f"  • {pair_name}:")
                click.echo(f"      Overlap: {pair_info['overlap']} genes")
                click.echo(f"      Jaccard: {pair_info['jaccard_index']:.3f}")
        click.echo()
        
        # Save
        click.echo("💾 Saving comparison...")
        comparison_file = output_path / 'comparison.json'
        
        # Create a summary without the full gene lists for readability
        summary = {
            'n_methods': comparison['n_methods'],
            'overlap_count': comparison['overlap_count'],
            'union_count': comparison['union_count'],
            'jaccard_index': comparison['jaccard_index'],
            'per_method': {
                name: info['n_significant']
                for name, info in comparison['per_method'].items()
            },
            'pairwise_overlaps': comparison['pairwise_overlaps']
        }
        
        # Save summary
        with open(comparison_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2)
        click.echo(f"   ✓ comparison.json (summary)")
        
        # Save detailed comparison with gene lists
        detailed_file = output_path / 'comparison_detailed.json'
        with open(detailed_file, 'w', encoding='utf-8') as f:
            json.dump(comparison, f, indent=2)
        click.echo(f"   ✓ comparison_detailed.json (with gene lists)")
        
        # Save overlap genes to CSV
        if comparison['overlap_count'] > 0:
            overlap_file = output_path / 'overlap_genes.txt'
            with open(overlap_file, 'w', encoding='utf-8') as f:
                f.write('\n'.join(comparison['overlap_genes']))
            click.echo(f"   ✓ overlap_genes.txt ({comparison['overlap_count']} genes)")
        
        click.echo()
        click.echo(f"   All files saved to: {output_path}/")
        click.echo()
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# MODULE 8: PARAMETER OPTIMIZATION
# =============================================================================

@main.command('optimize')
@click.option('--de-result', '-d', type=click.Path(exists=True), required=True,
              help='DE result file (.pkl)')
@click.option('--method', '-m', 
              type=click.Choice(['ground-truth', 'fdr-control', 'stability', 'reproducibility']),
              default='fdr-control',
              help='Optimization method')
@click.option('--ground-truth', '-g', type=click.Path(exists=True),
              help='Ground truth genes (CSV) - required for ground-truth method')
@click.option('--fdr-target', type=float, default=0.05,
              help='Target FDR for fdr-control method')
@click.option('--counts', type=click.Path(exists=True),
              help='Count matrix CSV - required for stability method')
@click.option('--metadata', type=click.Path(exists=True),
              help='Sample metadata CSV - required for stability method')
@click.option('--cohort2', type=click.Path(exists=True),
              help='Second cohort DE result (.pkl) - required for reproducibility method')
@click.option('--output', '-o', type=click.Path(), default='results/optimization',
              help='Output directory')
def optimize(de_result, method, ground_truth, fdr_target, counts, metadata, cohort2, output):
    """
    Optimize FDR and LFC thresholds (Module 8).
    
    \b
    OPTIMIZATION METHODS:
        ground-truth     : Maximize overlap with known genes
        fdr-control      : Control FDR at target level
        stability        : Bootstrap stability analysis
        reproducibility  : Cross-validation consistency
    
    \b
    OUTPUT:
        results/optimization/
        ├── optimal_parameters.json  # Best FDR/LFC thresholds
        ├── optimization_curve.png   # Parameter sweep plot
        └── validation_metrics.csv   # Performance metrics
    
    \b
    EXAMPLES:
        # With ground truth
        raptor optimize -d de_result.pkl -m ground-truth -g truth.csv
        
        # FDR control
        raptor optimize -d de_result.pkl -m fdr-control --fdr-target 0.01
        
        # Stability analysis
        raptor optimize -d de_result.pkl -m stability
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 8: Parameter Optimization")
        click.echo()
        
        # Validate inputs
        de_path = validate_cli_file(de_result, "DE result")
        output_path = validate_cli_directory(output, create=True)
        
        # Method-specific validation
        if method == 'ground-truth' and not ground_truth:
            click.echo("❌ --ground-truth required for ground-truth method", err=True)
            sys.exit(1)
        
        if method == 'stability' and (not counts or not metadata):
            click.echo("❌ --counts and --metadata required for stability method", err=True)
            sys.exit(1)
        
        if method == 'reproducibility' and not cohort2:
            click.echo("❌ --cohort2 required for reproducibility method", err=True)
            sys.exit(1)
        
        # Validate additional files
        if ground_truth:
            gt_path = validate_cli_file(ground_truth, "Ground truth")
        if counts:
            counts_path = validate_cli_file(counts, "Count matrix")
        if metadata:
            metadata_path = validate_cli_file(metadata, "Sample metadata")
        if cohort2:
            cohort2_path = validate_cli_file(cohort2, "Cohort 2 DE result")
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   DE result: {de_path}")
        click.echo(f"   Method: {method}")
        if ground_truth:
            click.echo(f"   Ground truth: {gt_path}")
        if method == 'fdr-control':
            click.echo(f"   FDR target: {fdr_target}")
        if counts:
            click.echo(f"   Counts: {counts_path}")
        if metadata:
            click.echo(f"   Metadata: {metadata_path}")
        if cohort2:
            click.echo(f"   Cohort 2: {cohort2_path}")
        click.echo(f"   Output: {output_path}")
        click.echo()
        
        # Import optimization module
        try:
            from raptor.parameter_optimization import (
                optimize_with_ground_truth,
                optimize_with_fdr_control,
                optimize_with_stability,
                optimize_with_reproducibility,
            )
            import pickle
            import pandas as pd
        except ImportError as e:
            click.echo(f"❌ Module 8 not available: {e}", err=True)
            sys.exit(1)
        
        # Load DE result
        click.echo("📂 Loading DE result...")
        with open(de_path, 'rb') as f:
            de_res = pickle.load(f)
        click.echo(f"   ✓ Loaded {de_res.n_genes} genes")
        click.echo()
        
        # Adapt column names: DEResult uses log2_fold_change/p_value/adjusted_p_value
        # but parameter_optimization expects log2FoldChange/pvalue/padj
        de_df = de_res.results_df.copy()
        col_remap = {
            'log2_fold_change': 'log2FoldChange',
            'p_value': 'pvalue',
            'adjusted_p_value': 'padj',
        }
        de_df.rename(columns={k: v for k, v in col_remap.items() if k in de_df.columns}, inplace=True)
        
        # Run optimization
        click.echo(f"🔬 Running {method} optimization...")
        
        if method == 'ground-truth':
            gt_df = pd.read_csv(gt_path)
            if 'gene_id' not in gt_df.columns:
                click.echo("❌ Ground truth must have 'gene_id' column", err=True)
                sys.exit(1)
            result = optimize_with_ground_truth(de_df, gt_df, output_dir=output_path)
            
        elif method == 'fdr-control':
            result = optimize_with_fdr_control(de_df, target_fdr=fdr_target, output_dir=output_path)
            
        elif method == 'stability':
            counts_df = pd.read_csv(counts_path, index_col=0)
            metadata_df = pd.read_csv(metadata_path)
            result = optimize_with_stability(de_df, counts_df, metadata_df, output_dir=output_path)
            
        elif method == 'reproducibility':
            with open(cohort2_path, 'rb') as f:
                de_res2 = pickle.load(f)
            de_df2 = de_res2.results_df.copy()
            de_df2.rename(columns={k: v for k, v in col_remap.items() if k in de_df2.columns}, inplace=True)
            result = optimize_with_reproducibility(de_df, de_df2, output_dir=output_path)
        
        click.echo(f"   ✓ Optimization complete")
        click.echo()
        
        # Display results
        click.echo("📊 Optimal Parameters:")
        click.echo("=" * 60)
        click.echo(f"Alpha (FDR) threshold: {result.best_parameters['alpha']:.4f}")
        click.echo(f"LFC threshold: {result.best_parameters['lfc_threshold']:.2f}")
        click.echo(f"Best Score ({result.metric}): {result.best_score:.4f}")
        click.echo(f"DEG genes at optimal: {result.n_deg_genes}")
        click.echo()
        
        # Files already saved by optimization functions
        click.echo("✅ Results saved:")
        click.echo(f"   {output_path}/optimized_params.yaml")
        click.echo(f"   {output_path}/deg_genes.csv ({result.n_deg_genes} genes)")
        click.echo(f"   {output_path}/optimization_result.pkl")
        click.echo(f"   {output_path}/convergence_history.json")
        click.echo()
        
        # Next steps
        click.echo("📝 Next steps:")
        click.echo(f"   1. Review optimal parameters in {output_path}/")
        click.echo(f"   2. Use optimized thresholds:")
        click.echo(f"      raptor filter-de --de-result {de_path} \\")
        click.echo(f"          --alpha {result.best_parameters['alpha']:.4f} \\")
        click.echo(f"          --lfc {result.best_parameters['lfc_threshold']:.2f}")
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# MODULE 9: ENSEMBLE ANALYSIS
# =============================================================================

@main.command('ensemble')
@click.option('--methods', '-m', required=True,
              type=click.Choice(['fisher', 'brown', 'rra']),
              multiple=True,
              help='Ensemble method(s) to use (voting/weighted temporarily disabled)')
@click.option('--deseq2', type=click.Path(exists=True),
              help='DESeq2 result file (.pkl)')
@click.option('--edger', type=click.Path(exists=True),
              help='edgeR result file (.pkl)')
@click.option('--limma', type=click.Path(exists=True),
              help='limma result file (.pkl)')
@click.option('--wilcoxon', type=click.Path(exists=True),
              help='Wilcoxon result file (.pkl)')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output directory')
@click.option('--threshold', type=float, default=0.05,
              help='Significance threshold')
@click.option('--min-methods', type=int, default=2,
              help='Minimum methods for voting/weighted (default: 2)')
def ensemble(methods, deseq2, edger, limma, wilcoxon, output, threshold, min_methods):
    """
    Ensemble analysis combining multiple DE methods (Module 9).
    
    \b
    ENSEMBLE METHODS (AVAILABLE):
        fisher   : Fisher's combined probability test (meta-analysis)
        brown    : Brown's method (accounts for correlation)
        rra      : Robust Rank Aggregation (rank-based)
    
    \b
    TEMPORARILY DISABLED:
        voting   : Majority voting (being upgraded)
        weighted : Weighted combination (being upgraded)
    
    \b
    REQUIRES:
        At least 2 DE result files (DESeq2, edgeR, limma, or Wilcoxon)
    
    \b
    OUTPUT:
        results/ensemble/
        ├── consensus_genes.csv          # Final gene list
        ├── ensemble_result.pkl          # Full results
        ├── statistics.json              # Summary stats
        └── direction_consistency.csv    # Direction agreement
    
    \b
    EXAMPLES:
        # Fisher's method
        raptor ensemble -m fisher --deseq2 d.pkl --edger e.pkl --limma l.pkl -o fisher/
        
        # Multiple methods
        raptor ensemble -m fisher -m brown -m rra \\
            --deseq2 d.pkl --edger e.pkl --limma l.pkl -o ensemble/
        
        # All three methods together
        raptor ensemble -m fisher -m brown -m rra \\
            --deseq2 d.pkl --edger e.pkl --limma l.pkl -o all_methods/
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 9: Ensemble Analysis")
        click.echo()
        
        # Collect DE results
        de_files = {}
        if deseq2:
            de_files['DESeq2'] = deseq2
        if edger:
            de_files['edgeR'] = edger
        if limma:
            de_files['limma'] = limma
        if wilcoxon:
            de_files['Wilcoxon'] = wilcoxon
        
        if len(de_files) < 2:
            click.echo("❌ Need at least 2 DE result files", err=True)
            sys.exit(1)
        
        # Validate output
        output_path = validate_cli_directory(output, create=True)
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Methods: {', '.join(methods)}")
        click.echo(f"   DE results: {len(de_files)}")
        for name, path in de_files.items():
            click.echo(f"      - {name}: {path}")
        click.echo(f"   Output: {output_path}")
        click.echo(f"   Threshold: {threshold}")
        click.echo()
        
        # Note about temporarily disabled methods
        click.echo("ℹ️  Note: voting and weighted methods temporarily disabled")
        click.echo("   (being upgraded to match fisher/brown/rra architecture)")
        click.echo("   Use fisher, brown, or rra for ensemble analysis")
        click.echo()
        
        # Import ensemble module
        try:
            from raptor.ensemble import (
                ensemble_fisher,
                ensemble_brown,
                ensemble_rra,
            )
            from raptor.de_import import DEResult
            import pickle
            import pandas as pd
        except ImportError as e:
            click.echo(f"❌ Module 9 not available: {e}", err=True)
            sys.exit(1)
        
        # Load DE results
        click.echo("📂 Loading DE results...")
        # Column mapping: DEResult uses snake_case, ensemble expects camelCase
        _col_remap = {
            'log2_fold_change': 'log2FoldChange',
            'p_value': 'pvalue',
            'adjusted_p_value': 'padj',
        }
        from types import SimpleNamespace
        de_results = {}
        for name, path in de_files.items():
            with open(path, 'rb') as f:
                loaded = pickle.load(f)
            # Extract DataFrame from DEResult
            if hasattr(loaded, 'results_df'):
                click.echo(f"   ✓ {name}: {loaded.n_significant} significant genes")
                df = loaded.results_df.copy()
            else:
                df = loaded if isinstance(loaded, pd.DataFrame) else pd.DataFrame(loaded)
                click.echo(f"   ✓ {name}: loaded")
            # Rename columns to what ensemble functions expect
            df.rename(columns={k: v for k, v in _col_remap.items() if k in df.columns}, inplace=True)
            # Ensure gene_id is a column (not just the index)
            if 'gene_id' not in df.columns and df.index.name == 'gene_id':
                df = df.reset_index()
            elif 'gene_id' not in df.columns:
                df['gene_id'] = df.index
            # Wrap in object with .data attribute (ensemble expects result.data['gene_id'])
            de_results[name] = SimpleNamespace(data=df)
        click.echo()
        
        # Run ensemble methods
        results = {}
        
        for method in methods:
            click.echo(f"🔬 Running {method} ensemble...")
            
            if method == 'fisher':
                result = ensemble_fisher(de_results, significance_threshold=threshold)
            elif method == 'brown':
                result = ensemble_brown(de_results, significance_threshold=threshold)
            elif method == 'rra':
                result = ensemble_rra(de_results, significance_threshold=threshold)
            
            results[method] = result
            click.echo(f"   ✓ {result.n_consensus_genes} consensus genes")
            click.echo()
        
        # Display results
        click.echo("📊 Ensemble Results:")
        click.echo("=" * 60)
        
        for method, result in results.items():
            click.echo(f"{method.upper()}:")
            click.echo(f"   Consensus genes: {result.n_consensus_genes}")
            
            stats = result.method_statistics.get('overall', {})
            if stats:
                click.echo(f"   Upregulated: {stats.get('n_upregulated', 0)}")
                click.echo(f"   Downregulated: {stats.get('n_downregulated', 0)}")
                click.echo(f"   Direction consistent: {stats.get('n_direction_consistent', 0)}")
                click.echo(f"   Inconsistent direction: {stats.get('n_inconsistent', 0)}")
            click.echo()
        
        # Save results
        click.echo("💾 Saving results...")
        
        for method, result in results.items():
            method_dir = output_path / method
            method_dir.mkdir(exist_ok=True)
            result.save(method_dir)
            click.echo(f"   ✓ {method} → {method_dir}/")
        
        click.echo()
        click.echo("✅ Ensemble analysis complete!")
        click.echo()
        click.echo("💡 Next steps:")
        click.echo("   • Use consensus genes for pathway enrichment")
        click.echo("   • Validate with known biology")
        click.echo("   • Check direction_consistency.csv for conflicts")
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command('ensemble-compare')
@click.option('--deseq2', type=click.Path(exists=True), required=True,
              help='DESeq2 result file (.pkl)')
@click.option('--edger', type=click.Path(exists=True), required=True,
              help='edgeR result file (.pkl)')
@click.option('--limma', type=click.Path(exists=True), required=True,
              help='limma result file (.pkl)')
@click.option('--output', '-o', type=click.Path(), required=True,
              help='Output directory')
@click.option('--threshold', type=float, default=0.05,
              help='Significance threshold')
def ensemble_compare(deseq2, edger, limma, output, threshold):
    """
    Compare all 5 ensemble methods side-by-side (Module 9).
    
    Runs Fisher, Brown, RRA, Voting, and Weighted ensemble methods
    and generates comprehensive comparison showing overlap, differences,
    and method-specific statistics.
    
    \b
    EXAMPLE:
        raptor ensemble-compare --deseq2 d.pkl --edger e.pkl --limma l.pkl \\
            --output comparison/ --threshold 0.05
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 9: Ensemble Method Comparison")
        click.echo()
        
        # Import ensemble module
        try:
            from raptor.ensemble import (
                ensemble_fisher,
                ensemble_brown,
                ensemble_rra,
                ensemble_voting,
                ensemble_weighted,
            )
            from raptor.de_import import DEResult
            import pickle
            import pandas as pd
        except ImportError as e:
            click.echo(f"❌ Module 9 not available: {e}", err=True)
            sys.exit(1)
        
        # Validate and load DE results
        click.echo("📂 Loading DE results...")
        
        # Column mapping: DEResult uses snake_case, ensemble expects camelCase
        _col_remap = {
            'log2_fold_change': 'log2FoldChange',
            'p_value': 'pvalue',
            'adjusted_p_value': 'padj',
        }
        from types import SimpleNamespace
        
        def _extract_de(obj):
            """Extract DataFrame from DEResult, rename columns, wrap with .data attribute."""
            df = obj.results_df.copy() if hasattr(obj, 'results_df') else obj
            df.rename(columns={k: v for k, v in _col_remap.items() if k in df.columns}, inplace=True)
            if 'gene_id' not in df.columns and df.index.name == 'gene_id':
                df = df.reset_index()
            elif 'gene_id' not in df.columns:
                df['gene_id'] = df.index
            return SimpleNamespace(data=df)
        
        with open(deseq2, 'rb') as f:
            deseq2_result = _extract_de(pickle.load(f))
        with open(edger, 'rb') as f:
            edger_result = _extract_de(pickle.load(f))
        with open(limma, 'rb') as f:
            limma_result = _extract_de(pickle.load(f))
        
        de_results = {
            'DESeq2': deseq2_result,
            'edgeR': edger_result,
            'limma': limma_result
        }
        
        click.echo(f"   ✓ Loaded 3 DE method results")
        click.echo()
        
        # Validate output
        output_path = validate_cli_directory(output, create=True)
        
        # Run all methods
        methods_results = {}
        
        click.echo("🔬 Running all ensemble methods...")
        click.echo()
        
        # 1. Fisher
        click.echo("   1/5 Fisher's method...")
        methods_results['fisher'] = ensemble_fisher(
            de_results=de_results,
            significance_threshold=threshold
        )
        
        # 2. Brown
        click.echo("   2/5 Brown's method...")
        methods_results['brown'] = ensemble_brown(
            de_results=de_results,
            significance_threshold=threshold
        )
        
        # 3. RRA
        click.echo("   3/5 RRA...")
        methods_results['rra'] = ensemble_rra(
            de_results=de_results,
            significance_threshold=threshold
        )
        
        # 4. Voting (≥2 methods)
        click.echo("   4/5 Voting (≥2 methods)...")
        methods_results['voting'] = ensemble_voting(
            de_results=de_results,
            min_methods=2
        )
        
        # 5. Weighted (equal weights)
        click.echo("   5/5 Weighted ensemble...")
        methods_results['weighted'] = ensemble_weighted(
            de_results=de_results,
        )
        
        click.echo()
        
        # Create comparison summary
        click.echo("📊 Comparison Results:")
        click.echo()
        click.echo(f"{'Method':<15} {'Genes':<10} {'Up':<10} {'Down':<10}")
        click.echo("=" * 50)
        
        for method_name, result in methods_results.items():
            # Handle both EnsembleResult objects and DataFrames
            if hasattr(result, 'n_consensus_genes'):
                n_genes = result.n_consensus_genes
                stats = result.method_statistics.get('overall', {})
                n_up = stats.get('n_upregulated', 0)
                n_down = stats.get('n_downregulated', 0)
            elif isinstance(result, pd.DataFrame):
                n_genes = len(result)
                n_up = int((result.get('log2FoldChange', result.get('log2_fold_change', pd.Series())) > 0).sum()) if len(result) > 0 else 0
                n_down = int((result.get('log2FoldChange', result.get('log2_fold_change', pd.Series())) < 0).sum()) if len(result) > 0 else 0
            else:
                n_genes = 0
                n_up = 0
                n_down = 0
            
            click.echo(f"{method_name.capitalize():<15} {n_genes:<10} {n_up:<10} {n_down:<10}")
        
        click.echo()
        
        # Save all results
        click.echo(f"💾 Saving results to {output_path}/...")
        
        for method_name, result in methods_results.items():
            method_dir = output_path / method_name
            method_dir.mkdir(exist_ok=True)
            if hasattr(result, 'save'):
                result.save(method_dir)
            elif isinstance(result, pd.DataFrame):
                result.to_csv(method_dir / 'consensus_genes.csv', index=False)
        
        # Save comparison summary
        summary_data = []
        for method_name, result in methods_results.items():
            if hasattr(result, 'n_consensus_genes'):
                stats = result.method_statistics.get('overall', {})
                summary_data.append({
                    'method': method_name,
                    'n_genes': result.n_consensus_genes,
                    'n_upregulated': stats.get('n_upregulated', 0),
                    'n_downregulated': stats.get('n_downregulated', 0),
                })
            else:
                n = len(result) if isinstance(result, pd.DataFrame) else 0
                summary_data.append({
                    'method': method_name,
                    'n_genes': n,
                    'n_upregulated': 0,
                    'n_downregulated': 0,
                })
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_path / 'comparison_summary.csv', index=False)
        
        click.echo(f"   ✅ All results saved")
        click.echo()
        click.echo(f"💡 Recommendation:")
        click.echo(f"   • Fisher/Brown: Maximum sensitivity (~500 genes)")
        click.echo(f"   • RRA: Robust ranking (~450 genes)")
        click.echo(f"   • Voting: Balanced confidence (~300 genes)")
        click.echo(f"   • Weighted: Quality-based (~350 genes)")
        click.echo()
        click.echo(f"   Choose based on your research goals:")
        click.echo(f"   - Exploratory → Fisher/Brown")
        click.echo(f"   - Validation → Voting (min_methods=3)")
        click.echo(f"   - Best of both → RRA or Weighted")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# MODULE 10: BIOMARKER DISCOVERY
# =============================================================================

@main.command('biomarker')
@click.option('--counts', '-c', type=click.Path(exists=True), required=True,
              help='Count matrix CSV (genes x samples)')
@click.option('--metadata', '-m', type=click.Path(exists=True), required=True,
              help='Sample metadata CSV')
@click.option('--group-column', '-g', default='condition',
              help='Column in metadata defining groups')
@click.option('--de-result', '-d', type=click.Path(exists=True),
              help='DE result pickle from M7 (.pkl)')
@click.option('--ensemble-result', '-e', type=click.Path(exists=True),
              help='Ensemble result pickle from M9 (.pkl)')
@click.option('--methods', type=str, multiple=True,
              help='Feature selection methods (elastic_net, boruta, mrmr, rfe, shap, wgcna)')
@click.option('--panel-size', type=int, default=None,
              help='Target panel size (auto if not specified)')
@click.option('--species', default='human',
              help='Species for annotation (human, mouse, rat)')
@click.option('--disease-term', default=None,
              help='Disease context for literature search')
@click.option('--no-annotate', is_flag=True,
              help='Skip biological annotation')
@click.option('--no-literature', is_flag=True,
              help='Skip literature mining')
@click.option('--no-ppi', is_flag=True,
              help='Skip STRING PPI query')
@click.option('--output', '-o', type=click.Path(), default='results/biomarkers',
              help='Output directory')
@click.option('--intent', type=click.Choice(
              ['diagnostic', 'exploratory'], case_sensitive=False),
              default=None,
              help='Biomarker intent — auto-configures enhanced analyses '
                   '(signature score, direction patterns, clinical metrics, '
                   'ratio biomarkers)')
@click.option('--prevalence', type=float, default=0.05,
              help='Disease prevalence for PPV/NPV calculation (default: 0.05)')
@click.option('--auto-panel-strategy',
              type=click.Choice(['kneedle', 'argmax', 'first_drop'],
                                case_sensitive=False),
              default='kneedle',
              help='Panel-size auto-detection strategy when --panel-size '
                   'is not specified. kneedle (default): Satopaa et al. '
                   '2011 with argmax fallback for saturated curves. '
                   'argmax: smallest size at maximum CV AUC. first_drop: '
                   'legacy pre-2.2.2 heuristic, retained for backward '
                   'compatibility.')
@click.option('--panel-sensitivity', type=float, default=1.0,
              help='Kneedle sensitivity parameter S (default: 1.0). Lower '
                   'values pick earlier knees, higher values pick later. '
                   'Ignored unless --auto-panel-strategy=kneedle.')
@click.option('--panel-size-strategy',
              type=click.Choice(['per_fold', 'consensus'],
                                case_sensitive=False),
              default='consensus',
              help='Where panel-size auto-detection runs. consensus '
                   '(default): detect once on full data, pin all CV folds '
                   'to that K. per_fold: each outer fold detects '
                   'independently. Consensus reduces Phi instability at '
                   'small n (Q5 sanity-check evidence).')
def biomarker(counts, metadata, group_column, de_result, ensemble_result,
              methods, panel_size, species, disease_term,
              no_annotate, no_literature, no_ppi, output, intent, prevalence,
              auto_panel_strategy, panel_sensitivity, panel_size_strategy):
    """
    Discover biomarker gene panel (Module 10).
    
    \b
    PIPELINE STAGES:
        1. Multi-method feature selection (LASSO, Boruta, mRMR, RFE, SHAP, WGCNA)
        2. Panel size optimization (kneedle algorithm + consensus pinning)
        3. Classification evaluation (nested CV with RF, SVM, XGBoost, LogReg)
        4. Biological annotation (MyGene, pathways, literature, PPI)
        5. Publication-ready report generation
    
    \b
    OUTPUT:
        results/biomarkers/
        ├── biomarker_panel.csv
        ├── ranked_genes.csv
        ├── classification_performance.csv
        ├── panel_curve.csv
        ├── biomarker_report.md
        ├── biomarker_result.pkl
        └── annotations/
    
    \b
    EXAMPLES:
        # Basic discovery
        raptor biomarker -c counts.csv -m metadata.csv -g condition
        
        # With M9 ensemble result
        raptor biomarker -c counts.csv -m metadata.csv -g condition \\
            -e ensemble_result.pkl --panel-size 10
        
        # Mouse data with disease context
        raptor biomarker -c counts.csv -m metadata.csv -g genotype \\
            --species mouse --disease-term "tauopathy"
        
        # Fast run (no annotation, no network queries)
        raptor biomarker -c counts.csv -m metadata.csv -g condition --no-annotate
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 10: Biomarker Discovery")
        click.echo()
        
        # Validate inputs
        counts_path = validate_cli_file(counts, "Count matrix")
        meta_path = validate_cli_file(metadata, "Metadata")
        output_path = validate_cli_directory(output, create=True)
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Counts: {counts_path}")
        click.echo(f"   Metadata: {meta_path}")
        click.echo(f"   Group column: {group_column}")
        if de_result:
            click.echo(f"   DE result: {de_result}")
        if ensemble_result:
            click.echo(f"   Ensemble result: {ensemble_result}")
        if panel_size:
            click.echo(f"   Target panel size: {panel_size}")
        click.echo(f"   Species: {species}")
        if disease_term:
            click.echo(f"   Disease term: {disease_term}")
        click.echo(f"   Annotation: {'off' if no_annotate else 'on'}")
        click.echo(f"   Auto panel strategy: {auto_panel_strategy}")
        click.echo(f"   Panel size strategy: {panel_size_strategy}")
        if intent:
            click.echo(f"   Intent: {intent}")
            click.echo(f"   Prevalence: {prevalence}")
        click.echo(f"   Output: {output_path}")
        click.echo()
        
        # Import M10
        try:
            from raptor.biomarker_discovery import discover_biomarkers
            import pickle
        except ImportError as e:
            click.echo(f"❌ Module 10 not available: {e}", err=True)
            click.echo("   Install dependencies: pip install scikit-learn", err=True)
            sys.exit(1)
        
        # Load optional upstream results
        de_res = None
        if de_result:
            de_path = validate_cli_file(de_result, "DE result")
            click.echo("📂 Loading DE result...")
            with open(de_path, 'rb') as f:
                de_res = pickle.load(f)
            click.echo(f"   ✓ Loaded {de_res.n_genes} genes")
        
        ens_res = None
        if ensemble_result:
            ens_path = validate_cli_file(ensemble_result, "Ensemble result")
            click.echo("📂 Loading ensemble result...")
            with open(ens_path, 'rb') as f:
                ens_res = pickle.load(f)
            click.echo(f"   ✓ Loaded ensemble result")
        
        click.echo()
        
        # Run biomarker discovery
        methods_list = list(methods) if methods else None
        
        result = discover_biomarkers(
            counts=str(counts_path),
            metadata=str(meta_path),
            group_column=group_column,
            de_result=de_res,
            ensemble_result=ens_res,
            methods=methods_list,
            target_panel_size=panel_size,
            species=species,
            disease_term=disease_term,
            annotate=not no_annotate,
            run_literature=not no_literature,
            run_ppi=not no_ppi,
            output_dir=str(output_path),
            intent=intent,
            prevalence=prevalence,
            auto_panel_strategy=auto_panel_strategy,
            panel_sensitivity=panel_sensitivity,
            panel_size_strategy=panel_size_strategy,
        )
        
        # Display summary
        click.echo()
        click.echo("📊 Results:")
        click.echo("=" * 60)
        click.echo(f"Panel size: {result.panel_size} genes")
        click.echo(f"Best classifier: {result.best_classifier}")
        
        best_res = result.classification_results.get(result.best_classifier)
        if best_res:
            click.echo(f"AUC: {best_res.auc:.3f}")
            click.echo(f"F1:  {best_res.f1:.3f}")
        
        if result.panel_optimization:
            click.echo(f"Panel AUC: {result.panel_optimization.optimal_auc:.3f}")
        
        # Enhanced analysis results (if intent was set)
        if intent and hasattr(result, 'clinical_metrics') and result.clinical_metrics:
            cm = result.clinical_metrics
            click.echo()
            click.echo(f"Enhanced Analysis ({intent}):")
            if 'youdens' in cm:
                click.echo(f"  Youden's J: {cm['youdens']['youdens_j']:.3f} "
                           f"(threshold={cm['youdens']['threshold']:.3f})")
            if 'bootstrap_ci' in cm:
                b = cm['bootstrap_ci']
                click.echo(f"  AUC 95% CI: [{b['ci_lower']:.3f}, {b['ci_upper']:.3f}]")
            if 'ppv_npv' in cm:
                p = cm['ppv_npv']
                click.echo(f"  PPV at {p['prevalence']:.1%} prevalence: {p['ppv']:.3f}")
                click.echo(f"  NPV at {p['prevalence']:.1%} prevalence: {p['npv']:.3f}")
        if intent and hasattr(result, 'ratio_result') and result.ratio_result:
            rr = result.ratio_result
            if rr.best_pair:
                click.echo(f"  Top ratio: {rr.best_pair.name} (AUC={rr.best_pair.auc:.3f})")
        
        click.echo()
        click.echo(f"Panel genes: {', '.join(result.panel[:10])}")
        if len(result.panel) > 10:
            click.echo(f"   ... and {len(result.panel) - 10} more")
        
        click.echo()
        click.echo("✅ Results saved:")
        click.echo(f"   {output_path}/biomarker_panel.csv")
        click.echo(f"   {output_path}/ranked_genes.csv")
        click.echo(f"   {output_path}/classification_performance.csv")
        click.echo(f"   {output_path}/biomarker_result.pkl")
        if not no_annotate:
            click.echo(f"   {output_path}/biomarker_report.md")
            click.echo(f"   {output_path}/annotations/")
        
        click.echo()
        click.echo("📝 Next steps:")
        click.echo(f"   1. Review panel: {output_path}/biomarker_panel.csv")
        click.echo(f"   2. Read report: {output_path}/biomarker_report.md")
        click.echo(f"   3. Validate on independent cohort:")
        click.echo(f"      raptor biomarker-validate \\")
        click.echo(f"          --panel {output_path}/biomarker_panel.csv \\")
        click.echo(f"          --counts validation_counts.csv \\")
        click.echo(f"          --metadata validation_metadata.csv")
        
    except ValidationError as e:
        click.echo(f"❌ Validation Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command('biomarker-survival')
@click.option('--counts', '-c', type=click.Path(exists=True), required=True,
              help='Count matrix CSV (genes x samples)')
@click.option('--clinical', type=click.Path(exists=True), required=True,
              help='Clinical data CSV with survival columns')
@click.option('--time-column', default='os_time',
              help='Survival time column name')
@click.option('--event-column', default='os_event',
              help='Event indicator column (1=event, 0=censored)')
@click.option('--de-result', '-d', type=click.Path(exists=True),
              help='DE result pickle to restrict candidates')
@click.option('--fdr-threshold', type=float, default=0.05,
              help='FDR threshold for Cox univariate screen')
@click.option('--output', '-o', type=click.Path(), default='results/survival_biomarkers',
              help='Output directory')
def biomarker_survival(counts, clinical, time_column, event_column,
                       de_result, fdr_threshold, output):
    """
    Discover prognostic biomarkers via survival analysis (Module 10D).
    
    \b
    WORKFLOW:
        1. Cox univariate screen (per-gene hazard ratios)
        2. FDR correction
        3. CoxNet panel selection (L1-penalized Cox regression)
        4. C-index evaluation
    
    \b
    REQUIRES: pip install lifelines
    
    \b
    OUTPUT:
        results/survival_biomarkers/
        ├── cox_univariate_screen.csv
        ├── survival_panel.csv
        └── survival_summary.json
    
    \b
    EXAMPLES:
        raptor biomarker-survival -c tcga_counts.csv \\
            --clinical tcga_clinical.csv \\
            --time-column OS.time --event-column OS
        
        raptor biomarker-survival -c counts.csv \\
            --clinical clinical.csv \\
            -d de_result.pkl --fdr-threshold 0.01
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 10D: Survival Biomarkers")
        click.echo()
        
        counts_path = validate_cli_file(counts, "Count matrix")
        clinical_path = validate_cli_file(clinical, "Clinical data")
        output_path = validate_cli_directory(output, create=True)
        
        click.echo(f"📋 Configuration:")
        click.echo(f"   Counts: {counts_path}")
        click.echo(f"   Clinical: {clinical_path}")
        click.echo(f"   Time column: {time_column}")
        click.echo(f"   Event column: {event_column}")
        click.echo(f"   FDR threshold: {fdr_threshold}")
        click.echo(f"   Output: {output_path}")
        click.echo()
        
        try:
            from raptor.biomarker_discovery import discover_survival_biomarkers
            import pickle
        except ImportError as e:
            click.echo(f"❌ Module 10 not available: {e}", err=True)
            click.echo("   Install: pip install scikit-learn lifelines", err=True)
            sys.exit(1)
        
        # Load optional DE genes
        de_genes = None
        if de_result:
            de_path = validate_cli_file(de_result, "DE result")
            with open(de_path, 'rb') as f:
                de_res = pickle.load(f)
            de_genes = de_res.significant_genes
            click.echo(f"   Using {len(de_genes)} DE genes as candidates")
        
        click.echo("🔬 Running survival analysis...")
        click.echo()
        
        result = discover_survival_biomarkers(
            counts=str(counts_path),
            clinical=str(clinical_path),
            time_column=time_column,
            event_column=event_column,
            de_genes=de_genes,
            fdr_threshold=fdr_threshold,
            output_dir=str(output_path),
        )
        
        click.echo()
        click.echo("📊 Results:")
        click.echo("=" * 60)
        click.echo(f"Significant prognostic genes: {len(result.significant_genes)}")
        click.echo(f"Panel genes: {len(result.panel_genes)}")
        click.echo(f"C-index: {result.c_index:.3f}")
        
        if result.panel_genes:
            click.echo(f"\nPanel: {', '.join(result.panel_genes[:10])}")
            if len(result.panel_genes) > 10:
                click.echo(f"   ... and {len(result.panel_genes) - 10} more")
        
        click.echo()
        click.echo("✅ Results saved:")
        click.echo(f"   {output_path}/cox_univariate_screen.csv")
        click.echo(f"   {output_path}/survival_panel.csv")
        click.echo(f"   {output_path}/survival_summary.json")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


@main.command('biomarker-validate')
@click.option('--panel', '-p', type=click.Path(exists=True), required=True,
              help='Biomarker panel CSV (from raptor biomarker)')
@click.option('--counts', '-c', type=click.Path(exists=True), required=True,
              help='Validation cohort count matrix CSV')
@click.option('--metadata', '-m', type=click.Path(exists=True), required=True,
              help='Validation cohort metadata CSV')
@click.option('--group-column', '-g', default='condition',
              help='Group column in metadata')
@click.option('--n-folds', type=int, default=5,
              help='Number of CV folds')
@click.option('--output', '-o', type=click.Path(), default='results/biomarker_validation',
              help='Output directory')
def biomarker_validate(panel, counts, metadata, group_column, n_folds, output):
    """
    Validate biomarker panel on independent cohort (Module 10).
    
    \b
    Tests a previously discovered panel on new data to assess
    whether the biomarkers replicate in an independent dataset.
    
    \b
    OUTPUT:
        results/biomarker_validation/
        ├── validation_performance.csv
        └── validation_summary.txt
    
    \b
    EXAMPLES:
        raptor biomarker-validate \\
            --panel results/biomarkers/biomarker_panel.csv \\
            --counts validation_counts.csv \\
            --metadata validation_metadata.csv \\
            --group-column condition
    """
    try:
        click.echo("🦖 RAPTOR v2.2.2 - Module 10: Biomarker Validation")
        click.echo()
        
        panel_path = validate_cli_file(panel, "Biomarker panel")
        counts_path = validate_cli_file(counts, "Validation counts")
        meta_path = validate_cli_file(metadata, "Validation metadata")
        output_path = validate_cli_directory(output, create=True)
        
        try:
            from raptor.biomarker_discovery import validate_biomarkers
            import pandas as pd
        except ImportError as e:
            click.echo(f"❌ Module 10 not available: {e}", err=True)
            sys.exit(1)
        
        # Load panel
        panel_df = pd.read_csv(panel_path)
        if 'gene_id' in panel_df.columns:
            panel_genes = panel_df['gene_id'].tolist()
        else:
            panel_genes = panel_df.iloc[:, 0].tolist()
        
        click.echo(f"📋 Panel: {len(panel_genes)} genes")
        click.echo(f"   Validation data: {counts_path}")
        click.echo(f"   Group column: {group_column}")
        click.echo()
        
        click.echo("🔬 Validating panel...")
        click.echo()
        
        results = validate_biomarkers(
            panel_genes=panel_genes,
            counts=str(counts_path),
            metadata=str(meta_path),
            group_column=group_column,
            n_folds=n_folds,
        )
        
        click.echo("📊 Validation Performance:")
        click.echo("=" * 60)
        
        rows = []
        for name, res in results.items():
            click.echo(f"   {name}: AUC={res.auc:.3f}, F1={res.f1:.3f}")
            rows.append({
                'classifier': name,
                'auc': res.auc,
                'accuracy': res.accuracy,
                'f1': res.f1,
                'sensitivity': res.sensitivity,
                'specificity': res.specificity,
            })
        
        import pandas as pd
        perf_df = pd.DataFrame(rows)
        perf_df.to_csv(output_path / 'validation_performance.csv', index=False)
        
        # Summary
        best = max(results.keys(), key=lambda k: results[k].auc)
        summary = (
            f"Validation of {len(panel_genes)}-gene panel\n"
            f"Best classifier: {best} (AUC={results[best].auc:.3f})\n"
        )
        with open(output_path / 'validation_summary.txt', 'w') as f:
            f.write(summary)
        
        click.echo()
        click.echo("✅ Results saved:")
        click.echo(f"   {output_path}/validation_performance.csv")
        click.echo(f"   {output_path}/validation_summary.txt")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        if logger.isEnabledFor(logging.DEBUG):
            import traceback
            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == '__main__':
    main()