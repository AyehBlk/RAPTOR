#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Production Pipeline (Module 5)

Demonstrates production-grade RNA-seq quantification pipelines.
User chooses based on needs: speed, BAM files, isoform-level, or bootstrap uncertainty.

This is Module 5 of the RAPTOR workflow (Stage 2: Production Pipeline):
  M1: Quantify (FASTQ → quick_gene_counts.csv)
  M2: Sample QC (Quality Assessment & Outlier Detection)
  M3: Profile (Data Profiling - 32 features)
  M4: Recommend (Pipeline Recommendation)
  M5: Pipeline (Production Pipeline) ← THIS SCRIPT

Input: Sample sheet CSV + Index (or --use-quantify to reuse M1 counts)
Output: results/production/
    - gene_counts.csv       (main output for DE)
    - tx_counts.csv         (transcript-level, if applicable)
    - tpm.csv               (TPM normalized)
    - bam/                   (if alignment-based)
    - pipeline_info.json    (run metadata)

Available Pipelines:
    - salmon              : Fast pseudo-alignment, no BAM
    - kallisto            : Ultra-fast pseudo-alignment, lowest memory
    - star_featurecounts  : Standard alignment, produces BAM
    - hisat2_featurecounts: Low-memory alignment alternative
    - star_rsem           : Gold standard for isoforms
    - star_salmon         : BAM + bootstraps (best of both worlds)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any

# =============================================================================
# CONSTANTS - Architecture Compliant (v2.2.0)
# =============================================================================
DEFAULT_QUICK_COUNTS_DIR = "results/quick_counts"
DEFAULT_OUTPUT_DIR = "results/production"
QUICK_COUNTS_FILE = "quick_gene_counts.csv"
OUTPUT_COUNTS_FILE = "gene_counts.csv"
OUTPUT_TX_FILE = "tx_counts.csv"
OUTPUT_TPM_FILE = "tpm.csv"
OUTPUT_INFO_FILE = "pipeline_info.json"

# Available pipelines with metadata
AVAILABLE_PIPELINES = {
    'salmon': {
        'name': 'Salmon',
        'description': 'Fast pseudo-alignment, no BAM output',
        'produces_bam': False,
        'gene_level': True,
        'isoform_level': True,
        'bootstrap': True,
        'memory_gb': 8,
        'speed': '⚡⚡⚡',
        'best_for': 'General DE analysis, speed-sensitive workflows'
    },
    'kallisto': {
        'name': 'Kallisto',
        'description': 'Ultra-fast pseudo-alignment, lowest memory',
        'produces_bam': False,
        'gene_level': True,
        'isoform_level': True,
        'bootstrap': True,
        'memory_gb': 4,
        'speed': '⚡⚡⚡⚡',
        'best_for': 'Large sample sizes, memory-constrained environments'
    },
    'star_featurecounts': {
        'name': 'STAR + featureCounts',
        'description': 'Standard alignment with BAM output',
        'produces_bam': True,
        'gene_level': True,
        'isoform_level': False,
        'bootstrap': False,
        'memory_gb': 32,
        'speed': '⚡⚡',
        'best_for': 'When BAM files needed, variant calling downstream'
    },
    'hisat2_featurecounts': {
        'name': 'HISAT2 + featureCounts',
        'description': 'Low-memory alignment alternative',
        'produces_bam': True,
        'gene_level': True,
        'isoform_level': False,
        'bootstrap': False,
        'memory_gb': 16,
        'speed': '⚡⚡',
        'best_for': 'BAM files with limited memory resources'
    },
    'star_rsem': {
        'name': 'STAR + RSEM',
        'description': 'Gold standard for isoform quantification',
        'produces_bam': True,
        'gene_level': True,
        'isoform_level': True,
        'bootstrap': False,
        'memory_gb': 32,
        'speed': '⚡',
        'best_for': 'Isoform-level DE, transcript switching analysis'
    },
    'star_salmon': {
        'name': 'STAR + Salmon',
        'description': 'BAM files + bootstrap uncertainty',
        'produces_bam': True,
        'gene_level': True,
        'isoform_level': True,
        'bootstrap': True,
        'memory_gb': 32,
        'speed': '⚡',
        'best_for': 'When you need both BAM AND uncertainty estimates'
    }
}

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
    from raptor.pipelines import get_pipeline, list_pipelines
    from raptor.pipelines.base import BasePipeline, SampleInfo, PipelineResult
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR pipeline modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.2.0 - Production Pipeline (Module 5)      ║
    ║                                                              ║
    ║   Production-Grade Quantification for DE Analysis            ║
    ║   6 Pipeline Options • Hybrid Dependency Handling            ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def print_workflow():
    """Print the RAPTOR workflow diagram."""
    print("""
    ┌─────────────────────────────────────────────────────────────────┐
    │                 RAPTOR v2.2.0 WORKFLOW                          │
    ├─────────────────────────────────────────────────────────────────┤
    │                                                                 │
    │  STAGE 1: Fast Profiling (M1-M4) ✓ COMPLETE                    │
    │  ═══════════════════════════════                               │
    │                                                                 │
    │  ┌──────────┐     ┌──────────┐     ┌──────────┐               │
    │  │  M1:     │ ──► │  M2:     │ ──► │  M3:     │               │
    │  │ Quantify │     │Sample QC │     │ Profile  │               │
    │  └──────────┘     └──────────┘     └──────────┘               │
    │                                          │                      │
    │                                          ▼                      │
    │                                    ┌──────────┐                 │
    │                                    │  M4:     │                 │
    │                                    │Recommend │                 │
    │                                    └──────────┘                 │
    │                                          │                      │
    │  ════════════════════════════════════════════════════════════  │
    │                                          │                      │
    │  STAGE 2: Production Pipeline (M5)       ▼                      │
    │  ══════════════════════════════    ┌──────────┐                │
    │                                    │  M5:     │ ◄── YOU ARE    │
    │                                    │ Pipeline │     HERE       │
    │                                    └──────────┘                 │
    │                                          │                      │
    │                                          ▼                      │
    │  ┌──────────────────────────────────────────────────────────┐  │
    │  │  OUTPUT FILES:                                           │  │
    │  │  • gene_counts.csv  (main input for DE analysis)        │  │
    │  │  • tx_counts.csv    (transcript-level, if applicable)   │  │
    │  │  • tpm.csv          (TPM normalized)                    │  │
    │  │  • bam/             (if alignment-based pipeline)       │  │
    │  └──────────────────────────────────────────────────────────┘  │
    │                                          │                      │
    │                                          ▼                      │
    │  STAGE 3: DE Analysis (M6-M10)                                 │
    │  ═════════════════════════════                                 │
    │                                                                 │
    └─────────────────────────────────────────────────────────────────┘
    
    MODULE 5: PRODUCTION PIPELINE
    ═════════════════════════════
    
    ┌─────────────────────────────────────────────────────────────────┐
    │  WHICH PIPELINE SHOULD I USE?                                   │
    │                                                                 │
    │  Need just gene counts for DE?                                  │
    │  ├── Yes → Use salmon or --use-quantify (fastest)              │
    │  └── No, need more →                                            │
    │      │                                                          │
    │      Need BAM files?                                            │
    │      ├── No → Use salmon or kallisto                           │
    │      └── Yes →                                                  │
    │          │                                                      │
    │          Need isoform-level counts?                             │
    │          ├── No → Use star_featurecounts or hisat2_featurecounts│
    │          └── Yes →                                              │
    │              │                                                  │
    │              Need bootstrap uncertainty?                        │
    │              ├── No → Use star_rsem (gold standard)            │
    │              └── Yes → Use star_salmon                         │
    │                                                                 │
    └─────────────────────────────────────────────────────────────────┘
    """)


def display_available_pipelines():
    """Display all available pipelines with details."""
    print("\n  ┌────────────────────────────────────────────────────────────────┐")
    print("  │  AVAILABLE PIPELINES                                           │")
    print("  ├────────────────────────────────────────────────────────────────┤")
    
    header = f"  │ {'Pipeline':<22} │ {'BAM':<3} │ {'Gene':<4} │ {'Iso':<3} │ {'Boot':<4} │ {'Speed':<6} │"
    print(header)
    print("  ├" + "─" * 64 + "┤")
    
    for pipeline_id, info in AVAILABLE_PIPELINES.items():
        bam = '✅' if info['produces_bam'] else '❌'
        gene = '✅' if info['gene_level'] else '❌'
        iso = '✅' if info['isoform_level'] else '❌'
        boot = '✅' if info['bootstrap'] else '❌'
        speed = info['speed']
        
        row = f"  │ {pipeline_id:<22} │ {bam:<3} │ {gene:<4} │ {iso:<3} │ {boot:<4} │ {speed:<6} │"
        print(row)
    
    print("  └────────────────────────────────────────────────────────────────┘")
    print("\n  Legend: BAM=Produces BAM files, Gene=Gene-level, Iso=Isoform-level, Boot=Bootstrap")


def display_pipeline_details(pipeline_id: str):
    """Display detailed information about a specific pipeline."""
    if pipeline_id not in AVAILABLE_PIPELINES:
        print(f"  ❌ Unknown pipeline: {pipeline_id}")
        return
    
    info = AVAILABLE_PIPELINES[pipeline_id]
    
    print(f"\n  ┌{'─' * 60}┐")
    print(f"  │  {info['name']:<56} │")
    print(f"  ├{'─' * 60}┤")
    print(f"  │  Description: {info['description']:<43} │")
    print(f"  │  Memory:      ~{info['memory_gb']} GB{' ' * 44}│")
    print(f"  │  Speed:       {info['speed']:<44} │")
    print(f"  │  Best for:    {info['best_for'][:44]:<44} │")
    print(f"  ├{'─' * 60}┤")
    print(f"  │  Features:                                                   │")
    print(f"  │    • Produces BAM:   {'Yes' if info['produces_bam'] else 'No':<36} │")
    print(f"  │    • Gene-level:     {'Yes' if info['gene_level'] else 'No':<36} │")
    print(f"  │    • Isoform-level:  {'Yes' if info['isoform_level'] else 'No':<36} │")
    print(f"  │    • Bootstrap:      {'Yes' if info['bootstrap'] else 'No':<36} │")
    print(f"  └{'─' * 60}┘")


def generate_demo_sample_sheet():
    """Generate demonstration sample sheet data."""
    return [
        {'sample_id': 'Control_1', 'condition': 'Control', 'batch': 'Batch1',
         'fastq_r1': '/path/to/Control_1_R1.fastq.gz', 'fastq_r2': '/path/to/Control_1_R2.fastq.gz'},
        {'sample_id': 'Control_2', 'condition': 'Control', 'batch': 'Batch1',
         'fastq_r1': '/path/to/Control_2_R1.fastq.gz', 'fastq_r2': '/path/to/Control_2_R2.fastq.gz'},
        {'sample_id': 'Control_3', 'condition': 'Control', 'batch': 'Batch2',
         'fastq_r1': '/path/to/Control_3_R1.fastq.gz', 'fastq_r2': '/path/to/Control_3_R2.fastq.gz'},
        {'sample_id': 'Treatment_1', 'condition': 'Treatment', 'batch': 'Batch1',
         'fastq_r1': '/path/to/Treatment_1_R1.fastq.gz', 'fastq_r2': '/path/to/Treatment_1_R2.fastq.gz'},
        {'sample_id': 'Treatment_2', 'condition': 'Treatment', 'batch': 'Batch2',
         'fastq_r1': '/path/to/Treatment_2_R1.fastq.gz', 'fastq_r2': '/path/to/Treatment_2_R2.fastq.gz'},
        {'sample_id': 'Treatment_3', 'condition': 'Treatment', 'batch': 'Batch2',
         'fastq_r1': '/path/to/Treatment_3_R1.fastq.gz', 'fastq_r2': '/path/to/Treatment_3_R2.fastq.gz'},
    ]


def generate_demo_counts(n_samples=6, n_genes=15000, seed=42):
    """Generate demonstration count matrix."""
    np.random.seed(seed)
    
    # Generate realistic RNA-seq counts using negative binomial
    base_expr = np.random.gamma(shape=2, scale=100, size=n_genes)
    
    counts = np.zeros((n_genes, n_samples))
    for i in range(n_samples):
        size_factor = np.random.uniform(0.8, 1.2)
        size_param = 10
        counts[:, i] = np.random.negative_binomial(
            size_param, 
            size_param / (size_param + base_expr * size_factor)
        )
    
    # Add some DE genes
    de_genes = np.random.choice(n_genes, 500, replace=False)
    half = n_samples // 2
    for gene in de_genes:
        fc = np.random.choice([-1, 1]) * np.random.uniform(0.5, 2)
        counts[gene, half:] *= (2 ** fc)
    
    gene_names = [f'ENSG{i+1:011d}' for i in range(n_genes)]
    sample_names = ['Control_1', 'Control_2', 'Control_3', 
                    'Treatment_1', 'Treatment_2', 'Treatment_3']
    
    counts_df = pd.DataFrame(
        counts.astype(int),
        index=gene_names,
        columns=sample_names[:n_samples]
    )
    
    return counts_df


def generate_demo_tpm(counts_df, seed=42):
    """Generate demonstration TPM matrix from counts."""
    np.random.seed(seed)
    
    gene_lengths = np.random.uniform(1000, 10000, size=len(counts_df))
    rpk = counts_df.div(gene_lengths / 1000, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    
    return tpm


def generate_demo_pipeline_result(pipeline_id: str, n_samples: int = 6):
    """Generate demonstration pipeline result."""
    
    class DemoPipelineResult:
        def __init__(self, pipeline_id, n_samples):
            self.pipeline = pipeline_id
            self.success = True
            self.n_samples = n_samples
            self.n_genes = 15000
            self.n_transcripts = 45000 if AVAILABLE_PIPELINES[pipeline_id]['isoform_level'] else None
            
            self.sample_results = []
            for i in range(n_samples):
                sample_id = f"{'Control' if i < n_samples//2 else 'Treatment'}_{(i % (n_samples//2)) + 1}"
                self.sample_results.append({
                    'sample_id': sample_id,
                    'success': True,
                    'mapping_rate': np.random.uniform(85, 95),
                    'num_processed': int(np.random.uniform(20e6, 30e6)),
                    'num_mapped': int(np.random.uniform(18e6, 28e6))
                })
            
            self.output_files = {
                'gene_counts': 'results/production/gene_counts.csv',
                'tpm': 'results/production/tpm.csv'
            }
            
            if AVAILABLE_PIPELINES[pipeline_id]['isoform_level']:
                self.output_files['tx_counts'] = 'results/production/tx_counts.csv'
            
            if AVAILABLE_PIPELINES[pipeline_id]['produces_bam']:
                self.output_files['bam_dir'] = 'results/production/bam/'
            
            self.runtime_seconds = np.random.uniform(600, 3600)
            self.memory_peak_gb = AVAILABLE_PIPELINES[pipeline_id]['memory_gb'] * np.random.uniform(0.8, 1.0)
        
        def summary(self):
            lines = []
            lines.append(f"\n  Pipeline: {AVAILABLE_PIPELINES[self.pipeline]['name']}")
            lines.append(f"  Status: {'✅ SUCCESS' if self.success else '❌ FAILED'}")
            lines.append(f"  Samples: {self.n_samples}")
            lines.append(f"  Genes: {self.n_genes:,}")
            if self.n_transcripts:
                lines.append(f"  Transcripts: {self.n_transcripts:,}")
            lines.append(f"  Runtime: {self.runtime_seconds/60:.1f} minutes")
            lines.append(f"  Peak Memory: {self.memory_peak_gb:.1f} GB")
            return '\n'.join(lines)
    
    return DemoPipelineResult(pipeline_id, n_samples)


def display_sample_results(result):
    """Display per-sample results."""
    print("\n  Per-Sample Results:")
    print("  ─" * 35)
    print(f"  {'Sample':<15} │ {'Status':<8} │ {'Mapping':<10} │ {'Reads':>12}")
    print("  ─" * 35)
    
    for sr in result.sample_results:
        status = '✅' if sr['success'] else '❌'
        mapping = f"{sr['mapping_rate']:.1f}%"
        reads = f"{sr['num_processed']:,}"
        print(f"  {sr['sample_id']:<15} │ {status:<8} │ {mapping:<10} │ {reads:>12}")
    
    print("  ─" * 35)
    
    avg_mapping = np.mean([sr['mapping_rate'] for sr in result.sample_results if sr['success']])
    total_reads = sum(sr['num_processed'] for sr in result.sample_results if sr['success'])
    print(f"  {'Average/Total':<15} │ {'─':<8} │ {avg_mapping:.1f}%{'':4} │ {total_reads:>12,}")


def display_output_files(result, output_dir):
    """Display output files."""
    print("\n  Output Files:")
    print("  ─" * 35)
    
    files = [
        (OUTPUT_COUNTS_FILE, 'Gene-level count matrix (main DE input)'),
        (OUTPUT_TPM_FILE, 'TPM normalized expression'),
    ]
    
    if result.n_transcripts:
        files.append((OUTPUT_TX_FILE, 'Transcript-level counts'))
    
    if AVAILABLE_PIPELINES[result.pipeline]['produces_bam']:
        files.append(('bam/', 'BAM files + indices'))
    
    files.append((OUTPUT_INFO_FILE, 'Pipeline run metadata'))
    
    for filename, description in files:
        print(f"  • {filename:<25} - {description}")


def run_pipeline(
    pipeline: str = 'salmon',
    sample_sheet: Optional[str] = None,
    index: Optional[str] = None,
    output: Optional[str] = None,
    threads: int = 8,
    use_quantify: bool = False,
    use_docker: bool = False,
    modules: Optional[str] = None,
    gtf: Optional[str] = None,
    salmon_index: Optional[str] = None,
    bootstraps: int = 0,
    demo: bool = False,
    **kwargs
) -> Dict[str, Any]:
    """
    Run production pipeline (Module 5).
    
    Parameters
    ----------
    pipeline : str
        Pipeline to run (salmon, kallisto, star_featurecounts, etc.)
    sample_sheet : str, optional
        Path to sample sheet CSV
    index : str, optional
        Path to index (Salmon/Kallisto index or STAR genome index)
    output : str, optional
        Output directory
    threads : int
        Number of threads
    use_quantify : bool
        Reuse Module 1 counts instead of running new pipeline
    use_docker : bool
        Use Docker containers for tools
    modules : str, optional
        HPC modules to load (comma-separated)
    gtf : str, optional
        GTF file (required for alignment-based pipelines)
    salmon_index : str, optional
        Salmon index (required for star_salmon pipeline)
    bootstraps : int
        Number of bootstraps (for Salmon/Kallisto)
    demo : bool
        Run in demo mode
    
    Returns
    -------
    Dict
        Results dictionary
    """
    output_dir = Path(output or DEFAULT_OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.2.0',
        'module': 'M5',
        'stage': 2,
        'pipeline': pipeline,
        'config': {
            'threads': threads,
            'use_docker': use_docker,
            'bootstraps': bootstraps
        }
    }
    
    # =========================================================================
    # Option A: Use-Quantify (Reuse M1 counts)
    # =========================================================================
    if use_quantify:
        print("\n🔄 USE-QUANTIFY MODE: Reusing Module 1 counts for production")
        print("─" * 60)
        
        quick_counts_dir = Path(DEFAULT_QUICK_COUNTS_DIR)
        quick_counts_file = quick_counts_dir / QUICK_COUNTS_FILE
        
        if demo or not quick_counts_file.exists():
            print("   📊 Demo: Generating simulated counts...")
            counts_df = generate_demo_counts()
            tpm_df = generate_demo_tpm(counts_df)
        else:
            print(f"   📂 Loading: {quick_counts_file}")
            counts_df = pd.read_csv(quick_counts_file, index_col=0)
            
            quick_tpm_file = quick_counts_dir / "quick_tpm.csv"
            if quick_tpm_file.exists():
                tpm_df = pd.read_csv(quick_tpm_file, index_col=0)
            else:
                tpm_df = generate_demo_tpm(counts_df)
        
        # Save as production counts
        counts_file = output_dir / OUTPUT_COUNTS_FILE
        tpm_file = output_dir / OUTPUT_TPM_FILE
        
        counts_df.to_csv(counts_file)
        tpm_df.to_csv(tpm_file)
        
        print(f"\n   ✅ Copied to production:")
        print(f"      • {counts_file}")
        print(f"      • {tpm_file}")
        
        results['mode'] = 'use_quantify'
        results['n_samples'] = counts_df.shape[1]
        results['n_genes'] = counts_df.shape[0]
        results['output_files'] = {
            'gene_counts': str(counts_file),
            'tpm': str(tpm_file)
        }
        
        # Save pipeline info
        info_file = output_dir / OUTPUT_INFO_FILE
        with open(info_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        return results
    
    # =========================================================================
    # Option B: Run Production Pipeline
    # =========================================================================
    print(f"\n🚀 RUNNING PRODUCTION PIPELINE: {AVAILABLE_PIPELINES[pipeline]['name']}")
    print("─" * 60)
    
    # Display pipeline info
    display_pipeline_details(pipeline)
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode...")
        
        # Generate demo results
        result = generate_demo_pipeline_result(pipeline)
        counts_df = generate_demo_counts()
        tpm_df = generate_demo_tpm(counts_df)
        
        # Display results
        print(result.summary())
        display_sample_results(result)
        
        # Save demo outputs
        counts_file = output_dir / OUTPUT_COUNTS_FILE
        tpm_file = output_dir / OUTPUT_TPM_FILE
        
        counts_df.to_csv(counts_file)
        tpm_df.to_csv(tpm_file)
        
        if result.n_transcripts:
            # Generate demo transcript counts
            tx_counts = generate_demo_counts(n_genes=45000)
            tx_counts.index = [f'ENST{i+1:011d}' for i in range(len(tx_counts))]
            tx_file = output_dir / OUTPUT_TX_FILE
            tx_counts.to_csv(tx_file)
        
        display_output_files(result, output_dir)
        
        results['mode'] = 'demo'
        results['n_samples'] = result.n_samples
        results['n_genes'] = result.n_genes
        results['n_transcripts'] = result.n_transcripts
        results['runtime_seconds'] = result.runtime_seconds
        results['output_files'] = result.output_files
        
    else:
        # Run actual pipeline
        print(f"   Sample sheet: {sample_sheet}")
        print(f"   Index: {index}")
        print(f"   Output: {output_dir}")
        print(f"   Threads: {threads}")
        
        # Get pipeline class
        PipelineClass = get_pipeline(pipeline)
        if not PipelineClass:
            print(f"❌ Pipeline '{pipeline}' not found")
            sys.exit(1)
        
        # Build pipeline kwargs
        pipeline_kwargs = {
            'output_dir': str(output_dir),
            'threads': threads,
            'use_docker': use_docker
        }
        
        if modules:
            pipeline_kwargs['modules'] = [m.strip() for m in modules.split(',')]
        
        if bootstraps > 0:
            pipeline_kwargs['bootstraps'] = bootstraps
        
        # Pipeline-specific parameters
        if pipeline in ['star_featurecounts', 'hisat2_featurecounts']:
            if not gtf:
                print("❌ --gtf is required for alignment-based pipelines")
                sys.exit(1)
            pipeline_kwargs['gtf'] = gtf
        
        if pipeline == 'star_salmon':
            if not salmon_index:
                print("❌ --salmon-index is required for star_salmon pipeline")
                sys.exit(1)
            pipeline_kwargs['salmon_index'] = salmon_index
        
        # Run pipeline
        try:
            pipeline_obj = PipelineClass(**pipeline_kwargs)
            result = pipeline_obj.run(sample_sheet, index)
            
            print(result.summary())
            
            if not result.success:
                print(f"❌ Pipeline failed")
                sys.exit(1)
            
            results['mode'] = 'production'
            results['success'] = result.success
            results['output_files'] = {
                'gene_counts': str(output_dir / OUTPUT_COUNTS_FILE)
            }
            
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
    
    # Save pipeline info
    info_file = output_dir / OUTPUT_INFO_FILE
    with open(info_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n   ✓ Saved: {info_file}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Production Pipeline (Module 5)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode (no data required)
  python 05_production_pipeline.py --demo
  
  # List available pipelines
  python 05_production_pipeline.py --list-pipelines
  
  # Reuse M1 counts (instant!)
  python 05_production_pipeline.py --use-quantify
  
  # Run Salmon pipeline
  python 05_production_pipeline.py --pipeline salmon -s samples.csv -i salmon_index/
  
  # Run STAR + featureCounts (needs GTF)
  python 05_production_pipeline.py --pipeline star_featurecounts -s samples.csv -i star_index/ --gtf genes.gtf
  
  # Run with Docker
  python 05_production_pipeline.py --pipeline salmon -s samples.csv -i index/ --use-docker
  
  # Run on HPC with modules
  python 05_production_pipeline.py --pipeline star_featurecounts -s samples.csv -i index/ --gtf genes.gtf --modules "STAR/2.7.10b,subread/2.0.3"

CLI Equivalent:
  raptor pipeline salmon -s samples.csv -i salmon_index/
  raptor pipeline star-featurecounts -s samples.csv -i star_index/ --gtf genes.gtf
  raptor pipeline use-quantify

Workflow:
  Stage 1 (M1-M4): Fast profiling → Recommendation ✓
  Stage 2 (M5): Production pipeline (this script) → gene_counts.csv
  Stage 3 (M6-M10): DE analysis → Results

Available Pipelines:
  salmon              : Fast pseudo-alignment, no BAM (⚡⚡⚡)
  kallisto            : Ultra-fast, lowest memory (⚡⚡⚡⚡)
  star_featurecounts  : Standard alignment with BAM (⚡⚡)
  hisat2_featurecounts: Low-memory alignment (⚡⚡)
  star_rsem           : Gold standard for isoforms (⚡)
  star_salmon         : BAM + bootstraps (⚡)

Output Location:
  results/production/
    ├── gene_counts.csv    (main input for DE)
    ├── tx_counts.csv      (transcript-level, if applicable)
    ├── tpm.csv            (TPM normalized)
    ├── bam/               (if alignment-based)
    └── pipeline_info.json (run metadata)
        """
    )
    
    # Pipeline selection
    parser.add_argument('--pipeline', '-p',
                        choices=list(AVAILABLE_PIPELINES.keys()),
                        default='salmon',
                        help='Pipeline to run (default: salmon)')
    parser.add_argument('--list-pipelines', action='store_true',
                        help='List available pipelines with details')
    parser.add_argument('--pipeline-info', type=str, metavar='NAME',
                        help='Show detailed info for a specific pipeline')
    
    # Input options
    parser.add_argument('--sample-sheet', '-s',
                        help='Sample sheet CSV file')
    parser.add_argument('--index', '-i',
                        help='Index path (Salmon/Kallisto index or STAR genome index)')
    
    # Output options
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    
    # Performance options
    parser.add_argument('--threads', '-t', type=int, default=8,
                        help='Number of threads (default: 8)')
    
    # Use-quantify option
    parser.add_argument('--use-quantify', action='store_true',
                        help='Reuse Module 1 counts as production (instant)')
    
    # Execution environment
    parser.add_argument('--use-docker', action='store_true',
                        help='Run tools in Docker containers')
    parser.add_argument('--modules', type=str,
                        help='HPC modules to load (comma-separated)')
    
    # Pipeline-specific options
    parser.add_argument('--gtf', '-g',
                        help='GTF file (required for alignment-based pipelines)')
    parser.add_argument('--salmon-index',
                        help='Salmon index (required for star_salmon pipeline)')
    parser.add_argument('--bootstraps', '-b', type=int, default=0,
                        help='Number of bootstraps (default: 0)')
    
    # Demo mode
    parser.add_argument('--demo', action='store_true',
                        help='Run in demo mode with simulated data')
    
    # Show workflow
    parser.add_argument('--show-workflow', action='store_true',
                        help='Show RAPTOR workflow diagram')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Handle info options
    if args.show_workflow:
        print_workflow()
        sys.exit(0)
    
    if args.list_pipelines:
        display_available_pipelines()
        sys.exit(0)
    
    if args.pipeline_info:
        display_pipeline_details(args.pipeline_info)
        sys.exit(0)
    
    # Validate inputs for real run
    if not args.demo and not args.use_quantify:
        if not args.sample_sheet:
            print("ERROR: --sample-sheet is required (or use --demo or --use-quantify)")
            parser.print_help()
            sys.exit(1)
        
        if not args.index:
            print("ERROR: --index is required (or use --demo or --use-quantify)")
            sys.exit(1)
    
    # Run pipeline
    results = run_pipeline(
        pipeline=args.pipeline,
        sample_sheet=args.sample_sheet,
        index=args.index,
        output=args.output,
        threads=args.threads,
        use_quantify=args.use_quantify,
        use_docker=args.use_docker,
        modules=args.modules,
        gtf=args.gtf,
        salmon_index=args.salmon_index,
        bootstraps=args.bootstraps,
        demo=args.demo
    )
    
    # Final summary
    print("\n" + "=" * 70)
    print("  ✅ MODULE 5 (PRODUCTION PIPELINE) COMPLETE!")
    print("=" * 70)
    
    output_dir = Path(args.output or DEFAULT_OUTPUT_DIR)
    
    print(f"\n  📂 Output Directory: {output_dir}")
    
    if 'n_genes' in results:
        print(f"\n  📊 Results:")
        print(f"     • Genes: {results['n_genes']:,}")
        if results.get('n_transcripts'):
            print(f"     • Transcripts: {results['n_transcripts']:,}")
        print(f"     • Samples: {results['n_samples']}")
    
    print(f"\n  📄 Output Files:")
    print(f"     • {OUTPUT_COUNTS_FILE:<25} - Main input for DE analysis")
    print(f"     • {OUTPUT_TPM_FILE:<25} - TPM normalized expression")
    if AVAILABLE_PIPELINES[args.pipeline]['isoform_level']:
        print(f"     • {OUTPUT_TX_FILE:<25} - Transcript-level counts")
    if AVAILABLE_PIPELINES[args.pipeline]['produces_bam']:
        print(f"     • bam/{'':23} - BAM files + indices")
    
    print(f"\n  🔜 Next Steps (Stage 3: DE Analysis):")
    print(f"\n     Run R script for differential expression:")
    print(f"     Rscript r_scripts/run_deseq2.R \\")
    print(f"         --counts {output_dir}/{OUTPUT_COUNTS_FILE} \\")
    print(f"         --metadata data/metadata.csv \\")
    print(f"         --output results/de_results.csv")
    print(f"\n     Then import results:")
    print(f"     raptor import-de --de-results results/de_results.csv")
    
    print("\n" + "=" * 70)
    print("  Making free science for everybody around the world 🌍")
    print("=" * 70 + "\n")


if __name__ == '__main__':
    main()
