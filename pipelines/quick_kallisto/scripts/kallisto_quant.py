#!/usr/bin/env python3

"""
Quick Kallisto Quantification Pipeline

Runs Kallisto quantification on all samples from a sample sheet
and combines results into a count matrix for QC and profiling.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0

IMPORTANT: For single-end reads, Kallisto REQUIRES fragment length
and standard deviation parameters. These must be estimated from
Bioanalyzer/TapeStation data or use reasonable defaults (~200, ~20).
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

# Clean import from raptor core
from raptor.sample_sheet import SampleSheet, Sample

logger = logging.getLogger(__name__)


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
    Quick Kallisto quantification pipeline for QC purposes.
    
    Runs Kallisto quant on all samples and combines into count matrix.
    
    Parameters
    ----------
    sample_sheet : str or Path
        Path to sample sheet CSV file
    index : str or Path
        Path to Kallisto index file (.idx)
    output_dir : str or Path
        Output directory
    threads : int
        Number of threads per sample
    gene_map : str or Path, optional
        Path to tx2gene mapping file
    config : dict, optional
        Additional configuration options
    
    Examples
    --------
    >>> pipeline = QuickKallistoPipeline(
    ...     sample_sheet='samples.csv',
    ...     index='kallisto_index.idx',
    ...     output_dir='counts/'
    ... )
    >>> pipeline.run()
    
    Notes
    -----
    For single-end reads, you MUST provide fragment_length and fragment_sd
    in the config. Kallisto cannot estimate these from single-end data.
    """
    
    def __init__(
        self,
        sample_sheet: str,
        index: str,
        output_dir: str,
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
        
        # Load sample sheet
        self.sample_sheet = SampleSheet(sample_sheet)
        
        # Validate inputs
        self._validate_inputs()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.quant_dir = self.output_dir / "kallisto_quant"
        self.quant_dir.mkdir(exist_ok=True)
        
        # Results storage
        self.results: List[KallistoResult] = []
    
    def _validate_inputs(self):
        """Validate all input files and directories exist."""
        errors = []
        
        # Check sample sheet
        if not self.sample_sheet_path.exists():
            errors.append(f"Sample sheet not found: {self.sample_sheet_path}")
        
        # Check index
        if not self.index.exists():
            errors.append(f"Kallisto index not found: {self.index}")
        elif not self.index.suffix == '.idx' and not self.index.is_file():
            logger.warning(f"Kallisto index may not be valid: {self.index}")
        
        # Check gene map if provided
        if self.gene_map and not self.gene_map.exists():
            errors.append(f"Gene map file not found: {self.gene_map}")
        
        # Validate sample sheet
        is_valid, sheet_errors = self.sample_sheet.validate(check_files=True)
        if not is_valid:
            errors.extend(sheet_errors)
        
        # Check single-end requirements
        if not self.sample_sheet.is_paired:
            frag_len = self.config.get('fragment_length')
            frag_sd = self.config.get('fragment_sd')
            
            if frag_len is None or frag_sd is None:
                errors.append(
                    "Single-end reads detected! Kallisto REQUIRES --fragment-length "
                    "and --fragment-sd parameters for single-end data. "
                    "Please provide 'fragment_length' and 'fragment_sd' in config."
                )
        
        if errors:
            for err in errors:
                logger.error(err)
            raise ValueError(f"Validation failed with {len(errors)} errors")
        
        logger.info(f"Validated {len(self.sample_sheet)} samples")
        logger.info(f"Read type: {'paired-end' if self.sample_sheet.is_paired else 'single-end'}")
    
    def _check_kallisto_installed(self) -> bool:
        """Check if Kallisto is installed and accessible."""
        try:
            result = subprocess.run(
                ["kallisto", "version"],
                capture_output=True,
                text=True
            )
            version = result.stdout.strip() or result.stderr.strip()
            logger.info(f"Kallisto version: {version}")
            return True
        except FileNotFoundError:
            logger.error("Kallisto not found. Please install: conda install -c bioconda kallisto")
            return False
    
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
            
            # Parse stats from run_info.json
            num_reads, num_pseudoaligned, rate = self._parse_kallisto_stats(sample_output)
            
            logger.info(f"✓ {sample.sample_id}: {num_reads:,} reads, {rate:.1f}% pseudoaligned")
            
            return KallistoResult(
                sample_id=sample.sample_id,
                success=True,
                quant_dir=str(sample_output),
                num_reads=num_reads,
                num_pseudoaligned=num_pseudoaligned,
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
        info_file = quant_dir / "run_info.json"
        
        if not info_file.exists():
            return 0, 0, 0.0
        
        try:
            with open(info_file) as f:
                info = json.load(f)
            
            num_reads = info.get("n_processed", 0)
            num_pseudoaligned = info.get("n_pseudoaligned", 0)
            
            rate = (num_pseudoaligned / num_reads * 100) if num_reads > 0 else 0.0
            
            return num_reads, num_pseudoaligned, rate
            
        except Exception as e:
            logger.warning(f"Could not parse Kallisto stats: {e}")
            return 0, 0, 0.0
    
    def run(self, parallel: bool = False) -> bool:
        """
        Run Kallisto quantification for all samples.
        
        Parameters
        ----------
        parallel : bool
            Run samples in parallel (experimental)
        
        Returns
        -------
        bool
            True if all samples succeeded
        """
        logger.info("=" * 60)
        logger.info("QUICK KALLISTO PIPELINE")
        logger.info("=" * 60)
        
        # Check Kallisto installation
        if not self._check_kallisto_installed():
            return False
        
        logger.info(f"Processing {len(self.sample_sheet)} samples...")
        logger.info(f"Output directory: {self.output_dir}")
        
        # Run quantification for each sample
        if parallel and len(self.sample_sheet) > 1:
            # Parallel execution (use with caution)
            with ThreadPoolExecutor(max_workers=2) as executor:
                futures = {
                    executor.submit(self._run_kallisto_quant, sample): sample
                    for sample in self.sample_sheet
                }
                for future in as_completed(futures):
                    result = future.result()
                    self.results.append(result)
        else:
            # Sequential execution (recommended)
            for sample in self.sample_sheet:
                result = self._run_kallisto_quant(sample)
                self.results.append(result)
        
        # Check for failures
        failures = [r for r in self.results if not r.success]
        if failures:
            logger.warning(f"{len(failures)} samples failed:")
            for f in failures:
                logger.warning(f"  - {f.sample_id}: {f.error_message}")
        
        # Combine results into count matrix
        if any(r.success for r in self.results):
            self._combine_counts()
            self._save_sample_info()
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
        
        # Save matrices
        counts_file = self.output_dir / "counts.csv"
        tpm_file = self.output_dir / "tpm.csv"
        
        counts_df.to_csv(counts_file)
        tpm_df.to_csv(tpm_file)
        
        logger.info(f"✓ Count matrix: {counts_file}")
        logger.info(f"  Shape: {counts_df.shape[0]} features × {counts_df.shape[1]} samples")
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
        counts_df = counts_df.dropna(subset=['gene_id'])
        tpm_df = tpm_df.dropna(subset=['gene_id'])
        
        # Aggregate
        counts_gene = counts_df.groupby('gene_id').sum()
        tpm_gene = tpm_df.groupby('gene_id').sum()
        
        logger.info(f"Aggregated {len(counts_df)} transcripts to {len(counts_gene)} genes")
        
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
        info_file = self.output_dir / "sample_info.csv"
        info_df.to_csv(info_file, index=False)
        
        logger.info(f"✓ Sample info: {info_file}")
    
    def cleanup(self, keep_quant: bool = False):
        """Clean up intermediate files."""
        if not keep_quant and self.quant_dir.exists():
            shutil.rmtree(self.quant_dir)
            logger.info("Cleaned up Kallisto quant directories")


def run_quick_kallisto(
    sample_sheet: str,
    index: str,
    output_dir: str,
    threads: int = 8,
    gene_map: Optional[str] = None,
    fragment_length: int = 200,
    fragment_sd: int = 20,
    keep_quant: bool = False
) -> bool:
    """
    Convenience function to run Quick Kallisto pipeline.
    
    Parameters
    ----------
    sample_sheet : str
        Path to sample sheet CSV
    index : str
        Path to Kallisto index (.idx file)
    output_dir : str
        Output directory
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
    ...     output_dir='counts/',
    ...     threads=8
    ... )
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


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Quick Kallisto Quantification for RAPTOR'
    )
    parser.add_argument('--sample-sheet', '-s', required=True,
                        help='Sample sheet CSV file')
    parser.add_argument('--index', '-i', required=True,
                        help='Kallisto index file (.idx)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
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
