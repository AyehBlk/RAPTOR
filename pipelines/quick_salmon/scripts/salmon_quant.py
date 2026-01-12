#!/usr/bin/env python3

"""
Quick Salmon Quantification Pipeline

Runs Salmon quantification on all samples from a sample sheet
and combines results into a count matrix for QC and profiling.

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

# Clean import from raptor core
from raptor.sample_sheet import SampleSheet, Sample

logger = logging.getLogger(__name__)


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
    Quick Salmon quantification pipeline for QC purposes.
    
    Runs Salmon quant on all samples and combines into count matrix.
    
    Parameters
    ----------
    sample_sheet : str or Path
        Path to sample sheet CSV file
    index : str or Path
        Path to Salmon index directory
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
    >>> pipeline = QuickSalmonPipeline(
    ...     sample_sheet='samples.csv',
    ...     index='salmon_index/',
    ...     output_dir='counts/'
    ... )
    >>> pipeline.run()
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
        self.quant_dir = self.output_dir / "salmon_quant"
        self.quant_dir.mkdir(exist_ok=True)
        
        # Results storage
        self.results: List[SalmonResult] = []
    
    def _validate_inputs(self):
        """Validate all input files and directories exist."""
        errors = []
        
        # Check sample sheet
        if not self.sample_sheet_path.exists():
            errors.append(f"Sample sheet not found: {self.sample_sheet_path}")
        
        # Check index
        if not self.index.exists():
            errors.append(f"Salmon index not found: {self.index}")
        elif not (self.index / "versionInfo.json").exists():
            errors.append(f"Invalid Salmon index (missing versionInfo.json): {self.index}")
        
        # Check gene map if provided
        if self.gene_map and not self.gene_map.exists():
            errors.append(f"Gene map file not found: {self.gene_map}")
        
        # Validate sample sheet
        is_valid, sheet_errors = self.sample_sheet.validate(check_files=True)
        if not is_valid:
            errors.extend(sheet_errors)
        
        if errors:
            for err in errors:
                logger.error(err)
            raise ValueError(f"Validation failed with {len(errors)} errors")
        
        logger.info(f"Validated {len(self.sample_sheet)} samples")
        logger.info(f"Read type: {'paired-end' if self.sample_sheet.is_paired else 'single-end'}")
    
    def _check_salmon_installed(self) -> bool:
        """Check if Salmon is installed and accessible."""
        try:
            result = subprocess.run(
                ["salmon", "--version"],
                capture_output=True,
                text=True
            )
            version = result.stdout.strip() or result.stderr.strip()
            logger.info(f"Salmon version: {version}")
            return True
        except FileNotFoundError:
            logger.error("Salmon not found. Please install: conda install -c bioconda salmon")
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
            "-l", self.config.get("lib_type", "A"),
            "-o", str(sample_output),
            "-p", str(self.threads)
        ]
        
        # Add input files based on read type
        if sample.is_paired:
            cmd.extend(["-1", sample.fastq_r1, "-2", sample.fastq_r2])
        else:
            cmd.extend(["-r", sample.fastq_r1])
            # Add fragment length for single-end (optional in Salmon but helps)
            frag_mean = self.config.get("fragment_length_mean", 200)
            frag_sd = self.config.get("fragment_length_sd", 80)
            cmd.extend(["--fldMean", str(frag_mean), "--fldSD", str(frag_sd)])
        
        # Add gene map if provided
        if self.gene_map:
            cmd.extend(["-g", str(self.gene_map)])
        
        # Add optional flags
        if self.config.get("validate_mappings", True):
            cmd.append("--validateMappings")
        
        if self.config.get("gc_bias", True):
            cmd.append("--gcBias")
        
        if self.config.get("seq_bias", False):
            cmd.append("--seqBias")
        
        # Add extra options
        extra = self.config.get("extra_options", "")
        if extra:
            cmd.extend(extra.split())
        
        logger.info(f"Running Salmon for {sample.sample_id}...")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                logger.error(f"Salmon failed for {sample.sample_id}: {result.stderr}")
                return SalmonResult(
                    sample_id=sample.sample_id,
                    success=False,
                    quant_dir=str(sample_output),
                    error_message=result.stderr[:500]
                )
            
            # Parse mapping stats
            num_reads, mapping_rate = self._parse_salmon_stats(sample_output)
            
            logger.info(f"✓ {sample.sample_id}: {num_reads:,} reads, {mapping_rate:.1f}% mapped")
            
            return SalmonResult(
                sample_id=sample.sample_id,
                success=True,
                quant_dir=str(sample_output),
                num_reads=num_reads,
                mapping_rate=mapping_rate
            )
            
        except subprocess.TimeoutExpired:
            logger.error(f"Salmon timed out for {sample.sample_id}")
            return SalmonResult(
                sample_id=sample.sample_id,
                success=False,
                quant_dir=str(sample_output),
                error_message="Timeout after 1 hour"
            )
        except Exception as e:
            logger.error(f"Error running Salmon for {sample.sample_id}: {e}")
            return SalmonResult(
                sample_id=sample.sample_id,
                success=False,
                quant_dir=str(sample_output),
                error_message=str(e)
            )
    
    def _parse_salmon_stats(self, quant_dir: Path) -> Tuple[int, float]:
        """Parse Salmon meta_info.json for statistics."""
        meta_file = quant_dir / "aux_info" / "meta_info.json"
        
        if not meta_file.exists():
            return 0, 0.0
        
        try:
            with open(meta_file) as f:
                meta = json.load(f)
            
            num_reads = meta.get("num_processed", 0)
            num_mapped = meta.get("num_mapped", 0)
            
            mapping_rate = (num_mapped / num_reads * 100) if num_reads > 0 else 0.0
            
            return num_reads, mapping_rate
            
        except Exception as e:
            logger.warning(f"Could not parse Salmon stats: {e}")
            return 0, 0.0
    
    def run(self, parallel: bool = False) -> bool:
        """
        Run Salmon quantification for all samples.
        
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
        logger.info("QUICK SALMON PIPELINE")
        logger.info("=" * 60)
        
        # Check Salmon installation
        if not self._check_salmon_installed():
            return False
        
        logger.info(f"Processing {len(self.sample_sheet)} samples...")
        logger.info(f"Output directory: {self.output_dir}")
        
        # Run quantification for each sample
        if parallel and len(self.sample_sheet) > 1:
            # Parallel execution (use with caution - high memory)
            with ThreadPoolExecutor(max_workers=2) as executor:
                futures = {
                    executor.submit(self._run_salmon_quant, sample): sample
                    for sample in self.sample_sheet
                }
                for future in as_completed(futures):
                    result = future.result()
                    self.results.append(result)
        else:
            # Sequential execution (recommended)
            for sample in self.sample_sheet:
                result = self._run_salmon_quant(sample)
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
        """Combine individual Salmon quant.sf files into count matrix."""
        logger.info("Combining counts into matrix...")
        
        counts_dict = {}
        tpm_dict = {}
        
        for result in self.results:
            if not result.success:
                continue
            
            quant_file = Path(result.quant_dir) / "quant.sf"
            
            if not quant_file.exists():
                logger.warning(f"quant.sf not found for {result.sample_id}")
                continue
            
            # Read Salmon output
            df = pd.read_csv(quant_file, sep='\t')
            
            # Extract counts and TPM
            counts_dict[result.sample_id] = df.set_index('Name')['NumReads']
            tpm_dict[result.sample_id] = df.set_index('Name')['TPM']
        
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
                'mapping_rate': result.mapping_rate,
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
            logger.info("Cleaned up Salmon quant directories")


def run_quick_salmon(
    sample_sheet: str,
    index: str,
    output_dir: str,
    threads: int = 8,
    gene_map: Optional[str] = None,
    lib_type: str = "A",
    keep_quant: bool = False
) -> bool:
    """
    Convenience function to run Quick Salmon pipeline.
    
    Parameters
    ----------
    sample_sheet : str
        Path to sample sheet CSV
    index : str
        Path to Salmon index
    output_dir : str
        Output directory
    threads : int
        Number of threads
    gene_map : str, optional
        Path to tx2gene mapping
    lib_type : str
        Salmon library type (default: A for auto)
    keep_quant : bool
        Keep individual quant files
    
    Returns
    -------
    bool
        True if successful
    
    Examples
    --------
    >>> from raptor.pipelines.quick_salmon.scripts.salmon_quant import run_quick_salmon
    >>> success = run_quick_salmon(
    ...     sample_sheet='samples.csv',
    ...     index='salmon_index/',
    ...     output_dir='counts/',
    ...     threads=8
    ... )
    """
    config = {
        'lib_type': lib_type,
        'validate_mappings': True,
        'gc_bias': True
    }
    
    pipeline = QuickSalmonPipeline(
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
        description='Quick Salmon Quantification for RAPTOR'
    )
    parser.add_argument('--sample-sheet', '-s', required=True,
                        help='Sample sheet CSV file')
    parser.add_argument('--index', '-i', required=True,
                        help='Salmon index directory')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--threads', '-t', type=int, default=8,
                        help='Number of threads')
    parser.add_argument('--gene-map', '-g',
                        help='Gene-to-transcript mapping file')
    parser.add_argument('--lib-type', '-l', default='A',
                        help='Salmon library type')
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
    
    success = run_quick_salmon(
        sample_sheet=args.sample_sheet,
        index=args.index,
        output_dir=args.output,
        threads=args.threads,
        gene_map=args.gene_map,
        lib_type=args.lib_type,
        keep_quant=args.keep_quant
    )
    
    sys.exit(0 if success else 1)
