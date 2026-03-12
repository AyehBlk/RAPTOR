#!/usr/bin/env python3

"""
Quick Kallisto Quantification Pipeline (Module 1: Quantify)

Runs Kallisto quantification on all samples from a sample sheet
and combines results into a count matrix for QC and profiling.

Output Files:
- results/quick_counts/quick_gene_counts.csv  (for Module 2-4)
- results/quick_counts/quick_tpm.csv
- results/quick_counts/sample_info.csv

IMPORTANT: For single-end reads, Kallisto REQUIRES fragment length
and standard deviation parameters. These must be estimated from
Bioanalyzer/TapeStation data or use reasonable defaults (~200, ~20).

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import os
import sys
import subprocess
import logging
import shutil
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

# Clean import from raptor.utils - v2.2.0
from raptor.utils.sample_sheet import SampleSheet, Sample

# v2.2.0 imports - Validation and error handling
from raptor.utils.errors import (
    handle_errors,
    ValidationError,
    PipelineError
)

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS - Architecture Compliant
# =============================================================================

# Default output directory (relative to project root)
DEFAULT_OUTPUT_DIR = "results/quick_counts"

# Output file names (distinct from Module 5 production counts)
OUTPUT_COUNTS_FILE = "quick_gene_counts.csv"
OUTPUT_TPM_FILE = "quick_tpm.csv"
OUTPUT_SAMPLE_INFO = "sample_info.csv"
OUTPUT_LOG_FILE = "quick_kallisto.log"

# =============================================================================
# CLI PARAMETERS - v2.2.0
# =============================================================================

QUICK_KALLISTO_CLI_PARAMS = {
    'sample_sheet': {
        'flag': '--sample-sheet',
        'short': '-s',
        'required': True,
        'type': 'file',
        'help': 'Sample sheet CSV file',
        'validation': 'Must be .csv file with required columns'
    },
    'index': {
        'flag': '--index',
        'short': '-i',
        'required': True,
        'type': 'file',
        'help': 'Kallisto index file (.idx)',
        'validation': 'Must be .idx file'
    },
    'output_dir': {
        'flag': '--output',
        'short': '-o',
        'default': 'results/quick_counts',
        'type': 'dir',
        'help': 'Output directory'
    },
    'threads': {
        'flag': '--threads',
        'short': '-t',
        'default': 8,
        'type': 'int',
        'range': [1, 128],
        'help': 'Number of threads per sample'
    },
    'fragment_length': {
        'flag': '--fragment-length',
        'short': '-l',
        'default': 200,
        'type': 'int',
        'range': [50, 1000],
        'help': 'Fragment length for single-end (REQUIRED!)',
        'validation': 'Required for single-end: 50-1000 bp'
    },
    'fragment_sd': {
        'flag': '--fragment-sd',
        'short': '-d',
        'default': 20,
        'type': 'int',
        'range': [10, 500],
        'help': 'Fragment length SD for single-end (REQUIRED!)',
        'validation': 'Required for single-end: 10-500 bp'
    },
    'bootstraps': {
        'flag': '--bootstraps',
        'short': '-b',
        'default': 0,
        'type': 'int',
        'range': [0, 1000],
        'help': 'Bootstrap samples (0=disabled)'
    },
    'gene_map': {
        'flag': '--gene-map',
        'short': '-g',
        'type': 'file',
        'help': 'Transcript-to-gene mapping'
    },
    'keep_quant': {
        'flag': '--keep-quant',
        'action': 'store_true',
        'help': 'Keep quant directories'
    }
}

# Pipeline metadata - v2.2.0
PIPELINE_NAME = "Quick Kallisto"
PIPELINE_VERSION = "2.2.0"
PIPELINE_MODULE = "M1"
PIPELINE_STAGE = 1


@dataclass
class KallistoResult:
    """Result from a single Kallisto quantification."""
    sample_id: str
    success: bool
    quant_dir: str
    num_reads: int = 0
    num_pseudoaligned: int = 0
    pseudoalign_rate: float = 0.0
    error_message: str = ""


class QuickKallistoPipeline:
    """
    Quick Kallisto quantification pipeline for QC purposes (Module 1).
    
    This is the FIRST stage of RAPTOR workflow:
    - Ultra-fast pseudo-alignment with Kallisto
    - Generates quick_gene_counts.csv for profiling (M2-M4)
    - NOT for final DE analysis (use Module 5 for production)
    
    Parameters
    ----------
    sample_sheet : str or Path
        Path to sample sheet CSV file
    index : str or Path
        Path to Kallisto index file (.idx)
    output_dir : str or Path
        Output directory (default: results/quick_counts)
    threads : int
        Number of threads per sample
    gene_map : str or Path, optional
        Path to tx2gene mapping file for gene-level aggregation
    config : dict, optional
        Additional configuration options
    
    Examples
    --------
    >>> pipeline = QuickKallistoPipeline(
    ...     sample_sheet='samples.csv',
    ...     index='kallisto_index.idx',
    ...     output_dir='results/quick_counts'
    ... )
    >>> pipeline.run()
    
    Output Files
    ------------
    - quick_gene_counts.csv : Gene-level count matrix (for M2-M4)
    - quick_tpm.csv : TPM normalized matrix
    - sample_info.csv : Sample QC metrics
    
    Notes
    -----
    For single-end reads, you MUST provide fragment_length and fragment_sd
    in the config. Kallisto cannot estimate these from single-end data.
    """
    
    @handle_errors  # NEW: v2.2.0 error handling
    def __init__(
        self,
        sample_sheet: str,
        index: str,
        output_dir: str = DEFAULT_OUTPUT_DIR,
        threads: int = 8,
        gene_map: Optional[str] = None,
        config: Optional[dict] = None
    ):
        self.sample_sheet_path = Path(sample_sheet)
        self.index = Path(index)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.gene_map = Path(gene_map) if gene_map else None
        self.config = config or {}
        
        # Extract config values - NEW: v2.2.0
        self.fragment_length = self.config.get('fragment_length', 200)
        self.fragment_sd = self.config.get('fragment_sd', 20)
        self.bootstraps = self.config.get('num_bootstraps', 0)
        
        # NEW: v2.2.0 validation
        self.validate_all()
        
        # Load sample sheet
        self.sample_sheet = SampleSheet(sample_sheet)
        
        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.quant_dir = self.output_dir / "kallisto_quant"
        self.quant_dir.mkdir(exist_ok=True)
        
        # Results storage
        self.results: List[KallistoResult] = []
        
        # Setup logging to file
        self._setup_logging()
    
    def _setup_logging(self):
        """Setup logging to both console and file."""
        log_file = self.output_dir / OUTPUT_LOG_FILE
        
        # File handler
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        logger.addHandler(fh)
    
    # =============================================================================
    # VALIDATION METHODS - v2.2.0
    # =============================================================================
    
    def validate_all(self) -> None:
        """Run all validation checks - v2.2.0."""
        self.validate_sample_sheet()
        self.validate_index()
        self.validate_threads()
        self.validate_fragment_length()
        self.validate_bootstraps()
        if self.gene_map:
            self.validate_gene_map()
    
    def validate_sample_sheet(self) -> None:
        """Validate sample sheet exists and has correct structure."""
        if not self.sample_sheet_path.exists():
            raise ValidationError(
                f"Sample sheet not found: {self.sample_sheet_path}\n"
                f"Please provide a valid CSV file with columns: sample_id, condition, batch, fastq_r1, fastq_r2"
            )
        
        if self.sample_sheet_path.suffix.lower() != '.csv':
            raise ValidationError(
                f"Sample sheet must be a .csv file, got: {self.sample_sheet_path.suffix}\n"
                f"File: {self.sample_sheet_path}"
            )
    
    def validate_index(self) -> None:
        """Validate Kallisto index file - v2.2.0."""
        if not self.index.exists():
            raise ValidationError(
                f"Kallisto index not found: {self.index}\n"
                f"Build index with: kallisto index -i index.idx transcripts.fa.gz\n"
                f"Download transcriptome from Ensembl or GENCODE"
            )
        
        if self.index.suffix != '.idx':
            raise ValidationError(
                f"Kallisto index must have .idx extension, got: {self.index.suffix}\n"
                f"File: {self.index}\n"
                f"Build index with: kallisto index -i transcripts.idx transcripts.fa.gz"
            )
    
    def validate_threads(self) -> None:
        """Validate thread count - v2.2.0."""
        if not (1 <= self.threads <= 128):
            raise ValidationError(
                f"Threads must be between 1-128, got: {self.threads}\n"
                f"Recommended: 4-8 threads for Kallisto (optimal performance)\n"
                f"Maximum: Use fewer threads than available CPUs"
            )
    
    def validate_fragment_length(self) -> None:
        """Validate fragment length for single-end - CRITICAL! v2.2.0."""
        # Will check again after loading sample sheet
        if self.fragment_length is not None:
            if not (50 <= self.fragment_length <= 1000):
                raise ValidationError(
                    f"Fragment length must be between 50-1000 bp, got: {self.fragment_length}\n"
                    f"Typical values:\n"
                    f"  - Illumina TruSeq: 200-300 bp\n"
                    f"  - Small RNA-seq: 50-100 bp\n"
                    f"  - Long fragment: 400-600 bp"
                )
        
        if self.fragment_sd is not None:
            if not (10 <= self.fragment_sd <= 500):
                raise ValidationError(
                    f"Fragment length SD must be between 10-500 bp, got: {self.fragment_sd}\n"
                    f"Typical values:\n"
                    f"  - Standard library prep: 20-30 bp\n"
                    f"  - Tight distribution: 10-20 bp\n"
                    f"  - Variable library: 40-80 bp"
                )
    
    def validate_bootstraps(self) -> None:
        """Validate bootstrap count - v2.2.0."""
        if not (0 <= self.bootstraps <= 1000):
            raise ValidationError(
                f"Bootstraps must be between 0-1000, got: {self.bootstraps}\n"
                f"Recommended values:\n"
                f"  - Module 1 (Quick): 0 (fastest, no uncertainty estimation)\n"
                f"  - Module 5 (Production): 100 (recommended for sleuth)\n"
                f"  - High precision: 200-500 (slower but better uncertainty)"
            )
    
    def validate_gene_map(self) -> None:
        """Validate gene map file if provided - v2.2.0."""
        if not self.gene_map.exists():
            raise ValidationError(
                f"Gene map file not found: {self.gene_map}\n"
                f"Transcript-to-gene mapping required for gene-level counts.\n"
                f"Supported formats: TSV, CSV, GTF\n"
                f"Example format (TSV):\n"
                f"  transcript_id\tgene_id\n"
                f"  ENST00000456328\tENSG00000223972"
            )
    
    def validate_single_end_requirements(self) -> None:
        """Validate single-end requirements after loading sample sheet - v2.2.0."""
        if not self.sample_sheet.is_paired:
            if self.fragment_length is None or self.fragment_sd is None:
                raise ValidationError(
                    f"❌ CRITICAL: Single-end reads detected but fragment length not provided!\n"
                    f"\n"
                    f"Kallisto CANNOT estimate fragment length from single-end data.\n"
                    f"You MUST provide both parameters:\n"
                    f"  --fragment-length <mean>  (e.g., 200)\n"
                    f"  --fragment-sd <sd>        (e.g., 20)\n"
                    f"\n"
                    f"How to estimate fragment length:\n"
                    f"  1. Bioanalyzer/TapeStation: Direct measurement (best)\n"
                    f"  2. From paired-end BAM: Use Picard CollectInsertSizeMetrics\n"
                    f"  3. Library prep protocol: Check manufacturer specifications\n"
                    f"  4. Default values: 200 ± 20 bp (typical Illumina TruSeq)\n"
                    f"\n"
                    f"Example command:\n"
                    f"  raptor quick-count -m kallisto -s samples.csv -i index.idx \\\n"
                    f"    --fragment-length 200 --fragment-sd 20"
                )
    
    def validate_kallisto_installed(self) -> None:
        """Validate Kallisto is installed - v2.2.0."""
        try:
            result = subprocess.run(
                ['kallisto', 'version'],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info(f"✓ Kallisto version: {result.stdout.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise ValidationError(
                f"Kallisto not found in PATH\n"
                f"\n"
                f"Install Kallisto with conda/mamba:\n"
                f"  conda install -c bioconda kallisto\n"
                f"  # OR\n"
                f"  mamba install -c bioconda kallisto\n"
                f"\n"
                f"Verify installation:\n"
                f"  kallisto version\n"
                f"\n"
                f"Required version: >= 0.46.0"
            )
    
    def _run_kallisto_quant(self, sample: Sample) -> KallistoResult:
        """
        Run Kallisto quantification for a single sample.
        
        Parameters
        ----------
        sample : Sample
            Sample object with FASTQ paths
        
        Returns
        -------
        KallistoResult
            Result object with status and metrics
        """
        sample_output = self.quant_dir / sample.sample_id
        
        # Build Kallisto command
        cmd = [
            "kallisto", "quant",
            "-i", str(self.index),
            "-o", str(sample_output),
            "-t", str(self.threads)
        ]
        
        # Add bootstrap if requested
        n_bootstrap = self.config.get('num_bootstraps', 0)
        if n_bootstrap > 0:
            cmd.extend(["-b", str(n_bootstrap)])
            seed = self.config.get('seed', 42)
            cmd.extend(["--seed", str(seed)])
        
        # Add strand-specific options
        strand = self.config.get('strand')
        if strand == 'rf-stranded':
            cmd.append('--rf-stranded')
        elif strand == 'fr-stranded':
            cmd.append('--fr-stranded')
        
        # Add input files based on read type
        if sample.is_paired:
            # Paired-end: just add the two files
            cmd.extend([sample.fastq_r1, sample.fastq_r2])
        else:
            # Single-end: MUST specify --single, -l (length), -s (sd)
            frag_len = self.config.get('fragment_length', 200)
            frag_sd = self.config.get('fragment_sd', 20)
            
            cmd.extend([
                "--single",
                "-l", str(frag_len),
                "-s", str(frag_sd),
                sample.fastq_r1
            ])
        
        # Add extra options
        extra = self.config.get('extra_options', '')
        if extra:
            cmd.extend(extra.split())
        
        logger.info(f"Running Kallisto for {sample.sample_id}...")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                logger.error(f"Kallisto failed for {sample.sample_id}: {result.stderr}")
                return KallistoResult(
                    sample_id=sample.sample_id,
                    success=False,
                    quant_dir=str(sample_output),
                    error_message=result.stderr[:500]
                )
            
            # Parse mapping stats
            num_reads, num_pseudo, rate = self._parse_kallisto_stats(sample_output)
            
            logger.info(f"✓ {sample.sample_id}: {num_reads:,} reads, {rate:.1f}% pseudoaligned")
            
            return KallistoResult(
                sample_id=sample.sample_id,
                success=True,
                quant_dir=str(sample_output),
                num_reads=num_reads,
                num_pseudoaligned=num_pseudo,
                pseudoalign_rate=rate
            )
            
        except subprocess.TimeoutExpired:
            logger.error(f"Kallisto timed out for {sample.sample_id}")
            return KallistoResult(
                sample_id=sample.sample_id,
                success=False,
                quant_dir=str(sample_output),
                error_message="Timeout after 1 hour"
            )
        except Exception as e:
            logger.error(f"Error running Kallisto for {sample.sample_id}: {e}")
            return KallistoResult(
                sample_id=sample.sample_id,
                success=False,
                quant_dir=str(sample_output),
                error_message=str(e)
            )
    
    def _parse_kallisto_stats(self, quant_dir: Path) -> Tuple[int, int, float]:
        """Parse Kallisto run_info.json for statistics."""
        run_file = quant_dir / "run_info.json"
        
        if not run_file.exists():
            return 0, 0, 0.0
        
        try:
            with open(run_file) as f:
                info = json.load(f)
            
            num_reads = info.get("n_processed", 0)
            num_pseudo = info.get("n_pseudoaligned", 0)
            rate = info.get("p_pseudoaligned", 0.0)
            
            return num_reads, num_pseudo, rate
        except Exception:
            return 0, 0, 0.0
    
    @handle_errors  # NEW: v2.2.0 error handling
    def run(self, parallel: bool = False) -> bool:
        """
        Run the Quick Kallisto pipeline - v2.2.0.
        
        Parameters
        ----------
        parallel : bool
            Run samples in parallel (experimental)
        
        Returns
        -------
        bool
            True if all samples succeeded
        """
        # NEW: v2.2.0 - Validate Kallisto installed
        self.validate_kallisto_installed()
        
        # NEW: v2.2.0 - Validate single-end requirements
        self.validate_single_end_requirements()
        
        logger.info("=" * 60)
        logger.info(f"🦖 RAPTOR {PIPELINE_NAME} (v{PIPELINE_VERSION})")
        logger.info("=" * 60)
        logger.info(f"Module: {PIPELINE_MODULE} | Stage: {PIPELINE_STAGE}")
        logger.info(f"Samples: {len(self.sample_sheet)}")
        logger.info(f"Output: {self.output_dir}")
        logger.info("")
        
        # Run quantification for each sample
        samples = list(self.sample_sheet.samples.values())
        
        if parallel and len(samples) > 1:
            # Parallel execution
            with ThreadPoolExecutor(max_workers=min(4, len(samples))) as executor:
                futures = {
                    executor.submit(self._run_kallisto_quant, s): s 
                    for s in samples
                }
                for future in as_completed(futures):
                    result = future.result()
                    self.results.append(result)
        else:
            # Sequential execution
            for i, sample in enumerate(samples, 1):
                logger.info(f"[{i}/{len(samples)}] Processing {sample.sample_id}")
                result = self._run_kallisto_quant(sample)
                self.results.append(result)
        
        # Check for failures
        failures = [r for r in self.results if not r.success]
        if failures:
            logger.warning(f"⚠ {len(failures)} samples failed:")
            for f in failures:
                logger.warning(f"  - {f.sample_id}: {f.error_message}")
        
        # Combine results into count matrix
        if any(r.success for r in self.results):
            self._combine_counts()
            self._save_sample_info()
            self._print_summary()
            return len(failures) == 0
        else:
            logger.error("All samples failed!")
            return False
    
    def _combine_counts(self):
        """Combine individual Kallisto abundance.tsv files into count matrix."""
        logger.info("Combining counts into matrix...")
        
        counts_dict = {}
        tpm_dict = {}
        
        for result in self.results:
            if not result.success:
                continue
            
            abundance_file = Path(result.quant_dir) / "abundance.tsv"
            
            if not abundance_file.exists():
                logger.warning(f"abundance.tsv not found for {result.sample_id}")
                continue
            
            # Read Kallisto output
            df = pd.read_csv(abundance_file, sep='\t')
            
            # Extract counts and TPM
            # Kallisto columns: target_id, length, eff_length, est_counts, tpm
            counts_dict[result.sample_id] = df.set_index('target_id')['est_counts']
            tpm_dict[result.sample_id] = df.set_index('target_id')['tpm']
        
        if not counts_dict:
            logger.error("No count data to combine!")
            return
        
        # Create DataFrames
        counts_df = pd.DataFrame(counts_dict)
        tpm_df = pd.DataFrame(tpm_dict)
        
        # Round counts to integers
        counts_df = counts_df.round().astype(int)
        
        # If gene map provided, aggregate to gene level
        if self.gene_map:
            counts_df, tpm_df = self._aggregate_to_gene_level(counts_df, tpm_df)
        
        # Save matrices with architecture-compliant names
        counts_file = self.output_dir / OUTPUT_COUNTS_FILE
        tpm_file = self.output_dir / OUTPUT_TPM_FILE
        
        counts_df.to_csv(counts_file)
        tpm_df.to_csv(tpm_file)
        
        # Also save as TSV for compatibility
        counts_df.to_csv(self.output_dir / "quick_gene_counts.tsv", sep='\t')
        
        logger.info(f"✓ Count matrix: {counts_file}")
        logger.info(f"  Shape: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
        logger.info(f"✓ TPM matrix: {tpm_file}")
    
    def _aggregate_to_gene_level(
        self, 
        counts_df: pd.DataFrame, 
        tpm_df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Aggregate transcript-level counts to gene level."""
        logger.info("Aggregating to gene level...")
        
        # Load gene map
        tx2gene = pd.read_csv(self.gene_map)
        
        # Expected columns: transcript_id, gene_id
        if 'transcript_id' not in tx2gene.columns or 'gene_id' not in tx2gene.columns:
            # Try common alternative column names
            if len(tx2gene.columns) >= 2:
                tx2gene.columns = ['transcript_id', 'gene_id'] + list(tx2gene.columns[2:])
            else:
                logger.warning("Could not parse gene map, keeping transcript level")
                return counts_df, tpm_df
        
        # Create mapping
        tx2gene_map = tx2gene.set_index('transcript_id')['gene_id'].to_dict()
        
        # Map transcripts to genes
        counts_df['gene_id'] = counts_df.index.map(tx2gene_map)
        tpm_df['gene_id'] = tpm_df.index.map(tx2gene_map)
        
        # Drop unmapped
        unmapped = counts_df['gene_id'].isna().sum()
        if unmapped > 0:
            logger.warning(f"  {unmapped} transcripts not mapped to genes")
        
        counts_df = counts_df.dropna(subset=['gene_id'])
        tpm_df = tpm_df.dropna(subset=['gene_id'])
        
        # Aggregate
        counts_gene = counts_df.groupby('gene_id').sum()
        tpm_gene = tpm_df.groupby('gene_id').sum()
        
        logger.info(f"  Aggregated to {len(counts_gene)} genes")
        
        return counts_gene, tpm_gene
    
    def _save_sample_info(self):
        """Save sample information and QC metrics."""
        info_data = []
        
        for result in self.results:
            sample = self.sample_sheet.get_sample(result.sample_id)
            
            info_data.append({
                'sample_id': result.sample_id,
                'success': result.success,
                'num_reads': result.num_reads,
                'num_pseudoaligned': result.num_pseudoaligned,
                'pseudoalign_rate': result.pseudoalign_rate,
                'read_type': 'paired' if (sample and sample.is_paired) else 'single',
                'condition': sample.condition if sample else None,
                'batch': sample.batch if sample else None
            })
        
        info_df = pd.DataFrame(info_data)
        info_file = self.output_dir / OUTPUT_SAMPLE_INFO
        info_df.to_csv(info_file, index=False)
        
        logger.info(f"✓ Sample info: {info_file}")
    
    def _print_summary(self):
        """Print pipeline summary."""
        successful = sum(1 for r in self.results if r.success)
        total_reads = sum(r.num_reads for r in self.results if r.success)
        avg_rate = sum(r.pseudoalign_rate for r in self.results if r.success) / max(successful, 1)
        
        logger.info("")
        logger.info("=" * 60)
        logger.info("📊 SUMMARY")
        logger.info("=" * 60)
        logger.info(f"  Samples processed: {successful}/{len(self.results)}")
        logger.info(f"  Total reads: {total_reads:,}")
        logger.info(f"  Average pseudoalignment rate: {avg_rate:.1f}%")
        logger.info("")
        logger.info("📂 OUTPUT FILES:")
        logger.info(f"  • {self.output_dir / OUTPUT_COUNTS_FILE}")
        logger.info(f"  • {self.output_dir / OUTPUT_TPM_FILE}")
        logger.info(f"  • {self.output_dir / OUTPUT_SAMPLE_INFO}")
        logger.info("")
        logger.info("📋 NEXT STEPS:")
        logger.info(f"  raptor qc --counts {self.output_dir / OUTPUT_COUNTS_FILE}")
        logger.info(f"  raptor profile --counts {self.output_dir / OUTPUT_COUNTS_FILE}")
        logger.info("=" * 60)
    
    def cleanup(self, keep_quant: bool = False):
        """Clean up intermediate files."""
        if not keep_quant and self.quant_dir.exists():
            shutil.rmtree(self.quant_dir)
            logger.info("Cleaned up Kallisto quant directories")


def run_quick_kallisto(
    sample_sheet: str,
    index: str,
    output_dir: str = DEFAULT_OUTPUT_DIR,
    threads: int = 8,
    gene_map: Optional[str] = None,
    fragment_length: int = 200,
    fragment_sd: int = 20,
    keep_quant: bool = False
) -> bool:
    """
    Convenience function to run Quick Kallisto pipeline (Module 1).
    
    Parameters
    ----------
    sample_sheet : str
        Path to sample sheet CSV
    index : str
        Path to Kallisto index (.idx file)
    output_dir : str
        Output directory (default: results/quick_counts)
    threads : int
        Number of threads
    gene_map : str, optional
        Path to tx2gene mapping
    fragment_length : int
        Fragment length mean (REQUIRED for single-end)
    fragment_sd : int
        Fragment length SD (REQUIRED for single-end)
    keep_quant : bool
        Keep individual quant files
    
    Returns
    -------
    bool
        True if successful
    
    Examples
    --------
    >>> from raptor.pipelines.quick_kallisto.scripts.kallisto_quant import run_quick_kallisto
    >>> success = run_quick_kallisto(
    ...     sample_sheet='samples.csv',
    ...     index='kallisto_index.idx',
    ...     output_dir='results/quick_counts',
    ...     threads=8
    ... )
    
    Output
    ------
    - results/quick_counts/quick_gene_counts.csv
    - results/quick_counts/quick_tpm.csv
    - results/quick_counts/sample_info.csv
    """
    config = {
        'fragment_length': fragment_length,
        'fragment_sd': fragment_sd,
        'num_bootstraps': 0
    }
    
    pipeline = QuickKallistoPipeline(
        sample_sheet=sample_sheet,
        index=index,
        output_dir=output_dir,
        threads=threads,
        gene_map=gene_map,
        config=config
    )
    
    success = pipeline.run()
    
    if not keep_quant:
        pipeline.cleanup()
    
    return success


# =============================================================================
# EXPORTS for CLI Integration - v2.2.0
# =============================================================================

__all__ = [
    'QuickKallistoPipeline',
    'run_quick_kallisto',
    'QUICK_KALLISTO_CLI_PARAMS',
    'PIPELINE_NAME',
    'PIPELINE_VERSION',
    'PIPELINE_MODULE',
    'PIPELINE_STAGE'
]


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Quick Kallisto Quantification (Module 1)'
    )
    parser.add_argument('--sample-sheet', '-s', required=True,
                        help='Sample sheet CSV file')
    parser.add_argument('--index', '-i', required=True,
                        help='Kallisto index file (.idx)')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--threads', '-t', type=int, default=8,
                        help='Number of threads')
    parser.add_argument('--gene-map', '-g',
                        help='Gene-to-transcript mapping file')
    parser.add_argument('--fragment-length', '-l', type=int, default=200,
                        help='Fragment length mean (for single-end)')
    parser.add_argument('--fragment-sd', type=int, default=20,
                        help='Fragment length SD (for single-end)')
    parser.add_argument('--keep-quant', action='store_true',
                        help='Keep individual quant files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose output')
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    success = run_quick_kallisto(
        sample_sheet=args.sample_sheet,
        index=args.index,
        output_dir=args.output,
        threads=args.threads,
        gene_map=args.gene_map,
        fragment_length=args.fragment_length,
        fragment_sd=args.fragment_sd,
        keep_quant=args.keep_quant
    )
    
    sys.exit(0 if success else 1)
