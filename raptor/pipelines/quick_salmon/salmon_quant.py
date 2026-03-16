#!/usr/bin/env python3

"""
Quick Salmon Quantification Pipeline (Module 1: Quantify)
RAPTOR v2.2.0 - UPDATED with Full Validation & Error Handling

Runs Salmon quantification on all samples from a sample sheet
and combines results into a count matrix for QC and profiling.

This is Module 1 (QC) - lightweight validation (8-10 checks)
For production analysis, use Module 5 pipelines with full validation.

Output Files:
- results/quick_counts/quick_gene_counts.csv  (for Module 2-4)
- results/quick_counts/quick_tpm.csv
- results/quick_counts/sample_info.csv

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
Updated: January 2026
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

# RAPTOR v2.2.0 imports
from raptor.utils.validation import (
    validate_file_path,
    validate_directory_path,
    validate_positive_integer,
    validate_choice,
    check_file_exists,
    check_directory_exists
)
from raptor.utils.errors import ValidationError, PipelineError, handle_errors

# Clean import from raptor.utils - v2.2.0
from raptor.utils.sample_sheet import SampleSheet, Sample

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
OUTPUT_LOG_FILE = "quick_salmon.log"

# Valid Salmon library types
VALID_LIBRARY_TYPES = ['A', 'IU', 'ISR', 'ISF', 'MU', 'MSR', 'MSF', 
                       'OU', 'OSR', 'OSF', 'SF', 'SR', 'U']


# =============================================================================
# CLI PARAMETER DEFINITIONS (v2.2.0 Registry Compatible)
# =============================================================================

QUICK_SALMON_CLI_PARAMS = {
    'sample_sheet': {
        'flag': '--sample-sheet',
        'short': '-s',
        'type': str,
        'required': True,
        'help': 'Sample sheet CSV file with columns: sample_id, condition, batch, fastq_r1, fastq_r2'
    },
    'index': {
        'flag': '--index',
        'short': '-i',
        'type': str,
        'required': True,
        'help': 'Salmon index directory'
    },
    'output_dir': {
        'flag': '--output-dir',
        'short': '-o',
        'type': str,
        'default': DEFAULT_OUTPUT_DIR,
        'help': f'Output directory (default: {DEFAULT_OUTPUT_DIR})'
    },
    'threads': {
        'flag': '--threads',
        'short': '-t',
        'type': int,
        'default': 8,
        'help': 'Number of threads per sample (default: 8)'
    },
    'gene_map': {
        'flag': '--gene-map',
        'short': '-g',
        'type': str,
        'required': False,
        'help': 'Transcript-to-gene mapping file (TSV/CSV/GTF) for gene-level aggregation'
    },
    'lib_type': {
        'flag': '--lib-type',
        'short': '-l',
        'type': str,
        'default': 'A',
        'choices': VALID_LIBRARY_TYPES,
        'help': 'Salmon library type (default: A = auto-detect)'
    },
    'validate_mappings': {
        'flag': '--validate-mappings',
        'action': 'store_true',
        'default': True,
        'help': 'Enable mapping validation for better accuracy (default: True)'
    },
    'gc_bias': {
        'flag': '--gc-bias',
        'action': 'store_true',
        'default': True,
        'help': 'Enable GC bias correction (default: True)'
    },
    'seq_bias': {
        'flag': '--seq-bias',
        'action': 'store_true',
        'default': False,
        'help': 'Enable sequence-specific bias correction (slower, default: False)'
    },
    'fragment_length_mean': {
        'flag': '--fragment-length',
        'type': int,
        'default': 200,
        'help': 'Fragment length mean for single-end reads (default: 200)'
    },
    'fragment_length_sd': {
        'flag': '--fragment-sd',
        'type': int,
        'default': 80,
        'help': 'Fragment length SD for single-end reads (default: 80)'
    },
    'keep_quant': {
        'flag': '--keep-quant',
        'action': 'store_true',
        'default': False,
        'help': 'Keep individual Salmon quant directories (default: False)'
    }
}


@dataclass
class SalmonResult:
    """Result from a single Salmon quantification."""
    sample_id: str
    success: bool
    quant_dir: str
    num_reads: int = 0
    mapping_rate: float = 0.0
    error_message: str = ""


class QuickSalmonPipeline:
    """
    Quick Salmon quantification pipeline for QC purposes (Module 1).
    UPDATED v2.2.0 with comprehensive validation and error handling.
    
    This is the FIRST stage of RAPTOR workflow:
    - Fast pseudo-alignment with Salmon
    - Generates quick_gene_counts.csv for profiling (M2-M4)
    - Lightweight validation (8-10 checks for QC)
    - NOT for final DE analysis (use Module 5 for production)
    
    Parameters
    ----------
    sample_sheet : str or Path
        Path to sample sheet CSV file
    index : str or Path
        Path to Salmon index directory
    output_dir : str or Path
        Output directory (default: results/quick_counts)
    threads : int
        Number of threads per sample (default: 8)
    gene_map : str or Path, optional
        Path to tx2gene mapping file for gene-level aggregation
    lib_type : str
        Salmon library type (default: 'A' for auto-detect)
    validate_mappings : bool
        Enable mapping validation (default: True)
    gc_bias : bool
        Enable GC bias correction (default: True)
    seq_bias : bool
        Enable sequence bias correction (default: False)
    fragment_length_mean : int
        Fragment length for single-end (default: 200)
    fragment_length_sd : int
        Fragment length SD for single-end (default: 80)
    
    Examples
    --------
    >>> from raptor.pipelines.quick_salmon import QuickSalmonPipeline
    >>> pipeline = QuickSalmonPipeline(
    ...     sample_sheet='samples.csv',
    ...     index='salmon_index/',
    ...     output_dir='results/quick_counts',
    ...     threads=8
    ... )
    >>> pipeline.run()
    
    Output Files
    ------------
    - quick_gene_counts.csv : Gene-level count matrix (for M2-M4)
    - quick_tpm.csv : TPM normalized matrix
    - sample_info.csv : Sample QC metrics
    
    Raises
    ------
    ValidationError
        If input validation fails
    PipelineError
        If pipeline execution fails
    """
    
    # Pipeline metadata
    PIPELINE_NAME = "quick_salmon"
    PIPELINE_VERSION = "2.2.0"
    MODULE = 1
    STAGE = 1
    DESCRIPTION = "Fast Salmon pseudo-alignment for QC and profiling"
    
    def __init__(
        self,
        sample_sheet: str,
        index: str,
        output_dir: str = DEFAULT_OUTPUT_DIR,
        threads: int = 8,
        gene_map: Optional[str] = None,
        lib_type: str = 'A',
        validate_mappings: bool = True,
        gc_bias: bool = True,
        seq_bias: bool = False,
        fragment_length_mean: int = 200,
        fragment_length_sd: int = 80
    ):
        """Initialize Quick Salmon Pipeline with v2.2.0 validation.
        
        NOTE: __init__ stores parameters only. Call validate_inputs() or run()
        to trigger full validation. This allows unit tests with fake paths.
        """
        # Store parameters as strings (tests expect string access)
        self.sample_sheet = str(sample_sheet)
        self.sample_sheet_path = Path(sample_sheet)
        self.index = Path(index)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.gene_map = Path(gene_map) if gene_map else None
        self.lib_type = lib_type
        self.validate_mappings = validate_mappings
        self.gc_bias = gc_bias
        self.seq_bias = seq_bias
        self.fragment_length_mean = fragment_length_mean
        self.fragment_length_sd = fragment_length_sd
        
        # Will be set during validate_inputs()
        self._sample_sheet_obj = None
        self.results: List[SalmonResult] = []
        self._validated = False
    
    @handle_errors(exit_on_error=False)
    def validate_inputs(self):
        """
        Run all validation checks (10 checks - lightweight for QC).
        
        Called automatically by run(), or manually for early validation.
        
        Raises
        ------
        ValidationError
            If input validation fails
        FileNotFoundError
            If required files/directories don't exist
        """
        logger.info("🦖 RAPTOR Quick Salmon Pipeline v2.2.0")
        logger.info("=" * 60)
        logger.info("🔍 Validating inputs...")
        
        # 1. Sample sheet exists and is valid CSV
        validate_file_path(self.sample_sheet_path, must_exist=True, file_type='.csv')
        
        # 2. Index directory exists
        validate_directory_path(self.index, must_exist=True)
        
        # 3. Validate Salmon index structure
        if not (self.index / "versionInfo.json").exists():
            raise ValidationError(
                f"Invalid Salmon index: {self.index}\n"
                f"Missing versionInfo.json - index may be corrupted or incomplete.\n"
                f"Rebuild index with: salmon index -t transcripts.fa -i {self.index}"
            )
        
        # 4. Threads must be positive and in range
        validate_positive_integer(self.threads, 'threads')
        if self.threads > 128:
            logger.warning(f"threads ({self.threads}) exceeds recommended max (128)")
        
        # 5. Library type must be valid
        validate_choice(self.lib_type, VALID_LIBRARY_TYPES, 'lib_type')
        
        # 6. Gene map validation (if provided)
        if self.gene_map:
            validate_file_path(self.gene_map, must_exist=True)
        
        # 7. Fragment length validation (for single-end)
        validate_positive_integer(self.fragment_length_mean, 'fragment_length_mean')
        if not (50 <= self.fragment_length_mean <= 1000):
            raise ValidationError(
                f"fragment_length_mean must be 50-1000, got {self.fragment_length_mean}"
            )
        validate_positive_integer(self.fragment_length_sd, 'fragment_length_sd')
        if not (10 <= self.fragment_length_sd <= 500):
            raise ValidationError(
                f"fragment_length_sd must be 10-500, got {self.fragment_length_sd}"
            )
        
        # 8. Load and validate sample sheet
        try:
            self._sample_sheet_obj = SampleSheet(self.sample_sheet)
        except Exception as e:
            raise ValidationError(
                f"Failed to load sample sheet: {self.sample_sheet}\n"
                f"Error: {str(e)}"
            )
        
        # 9. Validate sample sheet structure
        is_valid, sheet_errors = self._sample_sheet_obj.validate(check_files=False)
        if not is_valid:
            error_msg = "\n".join(f"  - {err}" for err in sheet_errors)
            raise ValidationError(
                f"Sample sheet validation failed: {self.sample_sheet}\n"
                f"Errors found:\n{error_msg}"
            )
        
        # 10. Check Salmon is installed
        if not self._check_salmon_installed():
            raise ValidationError(
                f"Salmon not found in PATH\n"
                f"Install with: conda install -c bioconda salmon"
            )
        
        logger.info(f"✅ Validation passed (10 checks)")
        logger.info(f"   • Samples: {len(self._sample_sheet_obj)}")
        logger.info(f"   • Read type: {'paired-end' if self._sample_sheet_obj.is_paired else 'single-end'}")
        logger.info(f"   • Library type: {self.lib_type}")
        logger.info(f"   • Threads: {self.threads}")
        
        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.quant_dir = self.output_dir / "salmon_quant"
        self.quant_dir.mkdir(exist_ok=True)
        
        # Setup logging to file
        self._setup_logging()
        
        self._validated = True
    
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
        
        logger.info(f"📝 Logging to: {log_file}")
    
    def _check_salmon_installed(self) -> bool:
        """Check if Salmon is installed and accessible."""
        try:
            result = subprocess.run(
                ["salmon", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            version = result.stdout.strip() or result.stderr.strip()
            logger.info(f"✓ Salmon version: {version}")
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    
    def _run_salmon_quant(self, sample: Sample) -> SalmonResult:
        """
        Run Salmon quantification for a single sample.
        
        Parameters
        ----------
        sample : Sample
            Sample object with FASTQ paths
        
        Returns
        -------
        SalmonResult
            Result object with status and metrics
        """
        sample_output = self.quant_dir / sample.sample_id
        
        # Build Salmon command
        cmd = [
            "salmon", "quant",
            "-i", str(self.index),
            "-l", self.lib_type,
            "-o", str(sample_output),
            "-p", str(self.threads)
        ]
        
        # Add input files based on read type
        if sample.is_paired:
            cmd.extend(["-1", sample.fastq_r1, "-2", sample.fastq_r2])
        else:
            cmd.extend(["-r", sample.fastq_r1])
            # Add fragment length for single-end
            cmd.extend([
                "--fldMean", str(self.fragment_length_mean),
                "--fldSD", str(self.fragment_length_sd)
            ])
        
        # Add gene map if provided
        if self.gene_map:
            cmd.extend(["-g", str(self.gene_map)])
        
        # Add optional flags
        if self.validate_mappings:
            cmd.append("--validateMappings")
        
        if self.gc_bias:
            cmd.append("--gcBias")
        
        if self.seq_bias:
            cmd.append("--seqBias")
        
        # Run Salmon
        logger.info(f"  Running Salmon for: {sample.sample_id}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )
            
            if result.returncode != 0:
                return SalmonResult(
                    sample_id=sample.sample_id,
                    success=False,
                    quant_dir=str(sample_output),
                    error_message=result.stderr[:500]
                )
            
            # Extract metrics from meta_info.json
            meta_file = sample_output / 'aux_info' / 'meta_info.json'
            num_reads = 0
            mapping_rate = 0.0
            
            if meta_file.exists():
                with open(meta_file, 'r') as f:
                    meta = json.load(f)
                    num_reads = meta.get('num_processed', 0)
                    num_mapped = meta.get('num_mapped', 0)
                    mapping_rate = (num_mapped / num_reads * 100) if num_reads > 0 else 0
            
            logger.info(f"  ✓ {sample.sample_id}: {num_reads:,} reads, {mapping_rate:.1f}% mapped")
            
            return SalmonResult(
                sample_id=sample.sample_id,
                success=True,
                quant_dir=str(sample_output),
                num_reads=num_reads,
                mapping_rate=mapping_rate
            )
            
        except subprocess.TimeoutExpired:
            return SalmonResult(
                sample_id=sample.sample_id,
                success=False,
                quant_dir=str(sample_output),
                error_message="Salmon timed out (>2 hours)"
            )
        except Exception as e:
            return SalmonResult(
                sample_id=sample.sample_id,
                success=False,
                quant_dir=str(sample_output),
                error_message=str(e)
            )
    
    @handle_errors(exit_on_error=False)
    def run(self) -> bool:
        """
        Run the complete Quick Salmon pipeline.
        
        Returns
        -------
        bool
            True if pipeline completed successfully
        
        Raises
        ------
        PipelineError
            If pipeline execution fails
        """
        logger.info("")
        logger.info("=" * 60)
        logger.info("🚀 Starting Quick Salmon quantification...")
        logger.info("=" * 60)
        
        # Validate inputs if not already done
        if not self._validated:
            self.validate_inputs()
        
        # Run Salmon on all samples
        logger.info(f"📊 Processing {len(self._sample_sheet_obj)} samples...")
        
        for sample in self._sample_sheet_obj.samples:
            result = self._run_salmon_quant(sample)
            self.results.append(result)
        
        # Check for failures
        failures = [r for r in self.results if not r.success]
        if failures:
            logger.warning(f"⚠️  {len(failures)} samples failed:")
            for f in failures:
                logger.warning(f"   • {f.sample_id}: {f.error_message[:100]}")
        
        if not any(r.success for r in self.results):
            raise PipelineError(
                "All samples failed quantification. "
                "Check Salmon log files in output directory."
            )
        
        # Combine count matrices
        logger.info("📊 Combining count matrices...")
        self._combine_counts()
        
        # Save sample info
        self._save_sample_info()
        
        # Print summary
        self._print_summary()
        
        logger.info("")
        logger.info("=" * 60)
        logger.info("✅ Pipeline completed successfully!")
        logger.info("=" * 60)
        
        return True
    
    def _combine_counts(self):
        """Combine Salmon quant.sf files into count matrices."""
        counts_dict = {}
        tpm_dict = {}
        
        for result in self.results:
            if not result.success:
                continue
            
            quant_file = Path(result.quant_dir) / 'quant.sf'
            
            if not quant_file.exists():
                logger.warning(f"  Missing quant.sf for: {result.sample_id}")
                continue
            
            # Load quant.sf
            quant_df = pd.read_csv(quant_file, sep='\t')
            
            # Extract counts and TPM
            counts_dict[result.sample_id] = quant_df.set_index('Name')['NumReads']
            tpm_dict[result.sample_id] = quant_df.set_index('Name')['TPM']
        
        if not counts_dict:
            raise PipelineError(
                "No successful quantifications to combine"
            )
        
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
        logger.info("📊 Aggregating to gene level...")
        
        # Load gene map
        tx2gene = pd.read_csv(self.gene_map, sep=None, engine='python')
        
        # Normalize column names
        tx2gene.columns = [c.lower().strip() for c in tx2gene.columns]
        
        # Find transcript and gene columns
        tx_col = None
        gene_col = None
        
        for col in tx2gene.columns:
            if 'transcript' in col or 'txname' in col or 'target' in col:
                tx_col = col
            if 'gene' in col:
                gene_col = col
        
        # Fallback to first two columns
        if tx_col is None or gene_col is None:
            tx2gene.columns = ['transcript_id', 'gene_id'] + list(tx2gene.columns[2:])
            tx_col = 'transcript_id'
            gene_col = 'gene_id'
        
        # Create mapping
        tx2gene_map = tx2gene.set_index(tx_col)[gene_col].to_dict()
        
        # Map transcripts to genes
        counts_df['gene_id'] = counts_df.index.map(tx2gene_map)
        tpm_df['gene_id'] = tpm_df.index.map(tx2gene_map)
        
        # Drop unmapped
        unmapped = counts_df['gene_id'].isna().sum()
        if unmapped > 0:
            logger.warning(f"  ⚠️  {unmapped} transcripts not mapped to genes")
        
        counts_df = counts_df.dropna(subset=['gene_id'])
        tpm_df = tpm_df.dropna(subset=['gene_id'])
        
        # Aggregate
        counts_gene = counts_df.groupby('gene_id').sum()
        tpm_gene = tpm_df.groupby('gene_id').sum()
        
        logger.info(f"  ✓ Aggregated to {len(counts_gene)} genes")
        
        return counts_gene, tpm_gene
    
    def _save_sample_info(self):
        """Save sample information and QC metrics."""
        info_data = []
        
        for result in self.results:
            sample = self._sample_sheet_obj.get_sample(result.sample_id)
            
            info_data.append({
                'sample_id': result.sample_id,
                'success': result.success,
                'num_reads': result.num_reads,
                'mapping_rate': result.mapping_rate,
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
        avg_mapping = sum(r.mapping_rate for r in self.results if r.success) / max(successful, 1)
        
        logger.info("")
        logger.info("=" * 60)
        logger.info("📊 SUMMARY")
        logger.info("=" * 60)
        logger.info(f"  Samples processed: {successful}/{len(self.results)}")
        logger.info(f"  Total reads: {total_reads:,}")
        logger.info(f"  Average mapping rate: {avg_mapping:.1f}%")
        logger.info("")
        logger.info("📂 OUTPUT FILES:")
        logger.info(f"  • {self.output_dir / OUTPUT_COUNTS_FILE}")
        logger.info(f"  • {self.output_dir / OUTPUT_TPM_FILE}")
        logger.info(f"  • {self.output_dir / OUTPUT_SAMPLE_INFO}")
        logger.info("")
        logger.info("📋 NEXT STEPS (RAPTOR Workflow):")
        logger.info(f"  1. raptor qc --counts {self.output_dir / OUTPUT_COUNTS_FILE}")
        logger.info(f"  2. raptor profile --counts {self.output_dir / OUTPUT_COUNTS_FILE}")
        logger.info(f"  3. raptor recommend")
        logger.info("=" * 60)
    
    def cleanup(self, keep_quant: bool = False):
        """
        Clean up intermediate files.
        
        Parameters
        ----------
        keep_quant : bool
            If True, keep Salmon quant directories
        """
        if not keep_quant and self.quant_dir.exists():
            shutil.rmtree(self.quant_dir)
            logger.info("🧹 Cleaned up Salmon quant directories")


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def run_quick_salmon(
    sample_sheet: str,
    index: str,
    output_dir: str = DEFAULT_OUTPUT_DIR,
    threads: int = 8,
    gene_map: Optional[str] = None,
    lib_type: str = "A",
    keep_quant: bool = False,
    **kwargs
) -> bool:
    """
    Convenience function to run Quick Salmon pipeline (Module 1).
    
    Parameters
    ----------
    sample_sheet : str
        Path to sample sheet CSV
    index : str
        Path to Salmon index
    output_dir : str
        Output directory (default: results/quick_counts)
    threads : int
        Number of threads (default: 8)
    gene_map : str, optional
        Path to tx2gene mapping
    lib_type : str
        Salmon library type (default: A for auto)
    keep_quant : bool
        Keep individual quant files (default: False)
    **kwargs
        Additional arguments passed to QuickSalmonPipeline
    
    Returns
    -------
    bool
        True if successful
    
    Examples
    --------
    >>> from raptor.pipelines.quick_salmon import run_quick_salmon
    >>> success = run_quick_salmon(
    ...     sample_sheet='samples.csv',
    ...     index='salmon_index/',
    ...     output_dir='results/quick_counts',
    ...     threads=8
    ... )
    
    Output
    ------
    - results/quick_counts/quick_gene_counts.csv
    - results/quick_counts/quick_tpm.csv
    - results/quick_counts/sample_info.csv
    """
    pipeline = QuickSalmonPipeline(
        sample_sheet=sample_sheet,
        index=index,
        output_dir=output_dir,
        threads=threads,
        gene_map=gene_map,
        lib_type=lib_type,
        **kwargs
    )
    
    success = pipeline.run()
    
    if not keep_quant:
        pipeline.cleanup()
    
    return success


# =============================================================================
# CLI INTERFACE
# =============================================================================

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Quick Salmon Quantification (Module 1) v2.2.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python salmon_quant.py -s samples.csv -i salmon_index/ -o results/quick_counts
  
  # With gene-level aggregation
  python salmon_quant.py -s samples.csv -i salmon_index/ -g tx2gene.csv
  
  # Custom library type
  python salmon_quant.py -s samples.csv -i salmon_index/ -l ISR --threads 16

Output:
  results/quick_counts/quick_gene_counts.csv  (for Module 2-4)
  results/quick_counts/quick_tpm.csv
  results/quick_counts/sample_info.csv
        """
    )
    
    # Add arguments from CLI params dict
    for param_name, param_info in QUICK_SALMON_CLI_PARAMS.items():
        args_dict = {
            'help': param_info['help']
        }
        
        # Add flags
        flags = [param_info['flag']]
        if 'short' in param_info:
            flags.append(param_info['short'])
        
        # Add type
        if 'type' in param_info:
            args_dict['type'] = param_info['type']
        
        # Add default
        if 'default' in param_info:
            args_dict['default'] = param_info['default']
        
        # Add choices
        if 'choices' in param_info:
            args_dict['choices'] = param_info['choices']
        
        # Add action
        if 'action' in param_info:
            args_dict['action'] = param_info['action']
        
        # Add required
        if param_info.get('required'):
            args_dict['required'] = True
        
        parser.add_argument(*flags, **args_dict)
    
    # Additional arguments
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose output')
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Run pipeline
    try:
        success = run_quick_salmon(
            sample_sheet=args.sample_sheet,
            index=args.index,
            output_dir=args.output_dir,
            threads=args.threads,
            gene_map=args.gene_map,
            lib_type=args.lib_type,
            keep_quant=args.keep_quant,
            validate_mappings=args.validate_mappings,
            gc_bias=args.gc_bias,
            seq_bias=args.seq_bias,
            fragment_length_mean=args.fragment_length,
            fragment_length_sd=args.fragment_sd
        )
        
        sys.exit(0 if success else 1)
        
    except (ValidationError, PipelineError) as e:
        logger.error(f"❌ {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
