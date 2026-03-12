#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Kallisto Pipeline
FULLY INTEGRATED WITH VALIDATION & ERROR HANDLING

Production pipeline using Kallisto for ultra-fast pseudo-alignment.
No BAM files produced.

Features:
- Ultra-fast pseudo-alignment
- Transcript-level quantification
- Optional bootstrap for uncertainty estimation
- Lower memory footprint than Salmon
- ✅ Full input validation with clear error messages
- ✅ Integrated error handling
- ✅ Compatible with CLI system

Requirements:
- Kallisto >= 0.48.0
- Kallisto index (built from transcriptome)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import pandas as pd

# ============================================================================
# RAPTOR v2.2.0 IMPORTS - VALIDATION & ERROR HANDLING
# ============================================================================
from raptor.utils.validation import (
    validate_file_path,
    validate_directory_path,
    validate_positive_integer,
    validate_numeric_range,
)
from raptor.utils.errors import (
    handle_errors,
    ValidationError,
    PipelineError,
    check_file_exists,
)

# Pipeline base classes
from ..base import (
    BasePipeline,
    ToolDependency,
    SampleInfo,
    PipelineResult
)

logger = logging.getLogger(__name__)


class KallistoPipeline(BasePipeline):
    """
    Kallisto pseudo-alignment pipeline.
    
    Kallisto uses exact k-mer matching for ultra-fast transcript quantification.
    Lower memory usage than Salmon, but no GC bias correction.
    
    🆕 v2.2.0 Features:
    - ✅ Full input validation
    - ✅ Clear error messages
    - ✅ Graceful error handling
    - ✅ Type checking
    
    Parameters
    ----------
    output_dir : str or Path
        Output directory
    threads : int
        Number of threads (default: 8, must be > 0)
    bootstraps : int
        Number of bootstraps for uncertainty estimation (default: 0, must be >= 0)
    fragment_length : int, optional
        Fragment length mean (REQUIRED for single-end reads, must be > 0)
    fragment_sd : int, optional
        Fragment length SD (REQUIRED for single-end reads, must be > 0)
    rf_stranded : bool
        Strand specific reads, first read reverse (default: False)
    fr_stranded : bool
        Strand specific reads, first read forward (default: False)
    seed : int
        Random seed for bootstrap sampling (default: 42)
    tx2gene : str or Path, optional
        Transcript-to-gene mapping file for gene-level aggregation
    
    Raises
    ------
    ValidationError
        If any input parameter fails validation
    FileNotFoundError
        If index file or tx2gene file doesn't exist
    PipelineError
        If pipeline execution fails
    
    Examples
    --------
    >>> # Paired-end (auto fragment length)
    >>> pipeline = KallistoPipeline(output_dir='results/', threads=16)
    >>> result = pipeline.run('samples.csv', 'kallisto.idx')
    >>> 
    >>> # Single-end (MUST specify fragment length)
    >>> pipeline = KallistoPipeline(
    ...     output_dir='results/',
    ...     threads=16,
    ...     fragment_length=200,
    ...     fragment_sd=20
    ... )
    """
    
    PIPELINE_NAME = "kallisto"
    PIPELINE_VERSION = "2.2.0"
    DESCRIPTION = "Kallisto pseudo-alignment (ultra-fast, no BAM)"
    DOCKER_IMAGE = "zlskidmore/kallisto:0.50.0"
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        output_dir: Union[str, Path],
        threads: int = 8,
        memory_gb: int = 16,  # Kallisto uses less memory
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = False,  # Kallisto doesn't produce BAM
        extra_args: Optional[Dict[str, str]] = None,
        # Kallisto-specific parameters
        bootstraps: int = 0,
        fragment_length: Optional[int] = None,
        fragment_sd: Optional[int] = None,
        rf_stranded: bool = False,
        fr_stranded: bool = False,
        plaintext: bool = False,
        fusion: bool = False,
        seed: int = 42,
        # Gene aggregation
        tx2gene: Optional[Union[str, Path]] = None
    ):
        """
        Initialize Kallisto pipeline with validation.
        
        🆕 v2.2.0: All inputs are validated with clear error messages.
        """
        # ====================================================================
        # STEP 1: VALIDATE ALL INPUTS (v2.2.0 Enhancement)
        # ====================================================================
        logger.info("🔍 Validating pipeline inputs...")
        
        # Validate output directory
        try:
            output_dir = validate_directory_path(
                output_dir,
                create_if_missing=True,
                description="Output directory"
            )
            logger.info(f"  ✓ Output directory: {output_dir}")
        except Exception as e:
            raise ValidationError(
                f"Cannot create output directory: {output_dir}\n"
                f"  Error: {e}\n"
                f"  Please check write permissions."
            ) from e
        
        # Validate numeric parameters
        try:
            validate_positive_integer(threads, 'threads')
            validate_positive_integer(memory_gb, 'memory_gb')
            logger.info(f"  ✓ Numeric parameters validated")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid numeric parameter: {e}\n"
                f"  All numeric values must be positive integers.\n"
                f"  Current values: threads={threads}, memory_gb={memory_gb}"
            ) from e
        
        # Validate bootstraps (can be 0)
        if bootstraps < 0:
            raise ValidationError(
                f"bootstraps must be >= 0, got {bootstraps}\n"
                f"  Use 0 to disable bootstraps (faster)\n"
                f"  Use 100+ for uncertainty estimation (slower)"
            )
        
        # Validate seed
        try:
            validate_positive_integer(seed, 'seed')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid seed: {e}\n"
                f"  Seed must be a positive integer for reproducibility"
            ) from e
        
        # Validate fragment length/SD if provided
        if fragment_length is not None:
            try:
                validate_positive_integer(fragment_length, 'fragment_length')
                logger.info(f"  ✓ Fragment length: {fragment_length}")
            except ValidationError as e:
                raise ValidationError(
                    f"Invalid fragment_length: {e}\n"
                    f"  Fragment length must be positive (typically 150-300)"
                ) from e
        
        if fragment_sd is not None:
            try:
                validate_positive_integer(fragment_sd, 'fragment_sd')
                logger.info(f"  ✓ Fragment SD: {fragment_sd}")
            except ValidationError as e:
                raise ValidationError(
                    f"Invalid fragment_sd: {e}\n"
                    f"  Fragment SD must be positive (typically 10-30)"
                ) from e
        
        # Validate fragment length parameters relationship
        if (fragment_length is None) != (fragment_sd is None):
            raise ValidationError(
                "Both fragment_length AND fragment_sd must be specified together\n"
                f"  Got: fragment_length={fragment_length}, fragment_sd={fragment_sd}\n"
                f"  For single-end reads, provide both parameters"
            )
        
        # Validate strand-specific options
        if rf_stranded and fr_stranded:
            raise ValidationError(
                "Cannot specify both --rf-stranded and --fr-stranded\n"
                f"  Choose one:\n"
                f"    --rf-stranded: First read reverse (dUTP protocol)\n"
                f"    --fr-stranded: First read forward (ligation protocol)"
            )
        
        # Validate tx2gene file if provided
        if tx2gene:
            try:
                self.tx2gene = validate_file_path(
                    tx2gene,
                    must_exist=True,
                    description="Transcript-to-gene mapping file"
                )
                logger.info(f"  ✓ tx2gene file: {self.tx2gene}")
            except FileNotFoundError as e:
                raise ValidationError(
                    f"tx2gene file not found: {tx2gene}\n"
                    f"  Create this file with format:\n"
                    f"    transcript_id<TAB>gene_id\n"
                    f"  Example:\n"
                    f"    ENST00000456328	ENSG00000223972\n"
                    f"    ENST00000450305	ENSG00000223972"
                ) from e
        else:
            self.tx2gene = None
        
        logger.info("✅ All inputs validated successfully")
        
        # ====================================================================
        # STEP 2: INITIALIZE BASE PIPELINE
        # ====================================================================
        super().__init__(
            output_dir=output_dir,
            threads=threads,
            memory_gb=memory_gb,
            use_docker=use_docker,
            docker_image=docker_image,
            modules=modules,
            keep_bam=keep_bam,
            extra_args=extra_args
        )
        
        # ====================================================================
        # STEP 3: STORE VALIDATED PARAMETERS
        # ====================================================================
        
        # Kallisto-specific parameters
        self.bootstraps = bootstraps
        self.fragment_length = fragment_length
        self.fragment_sd = fragment_sd
        self.rf_stranded = rf_stranded
        self.fr_stranded = fr_stranded
        self.plaintext = plaintext
        self.fusion = fusion
        self.seed = seed
        
        # Create validated directories
        self.quant_dir = self.output_dir / "quant"
        self.quant_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info("🦖 Kallisto pipeline initialized")
    
    def get_dependencies(self) -> List[ToolDependency]:
        """Return Kallisto dependency."""
        return [
            ToolDependency(
                name="Kallisto",
                command="kallisto",
                min_version="0.48.0",
                docker_image="zlskidmore/kallisto:0.50.0",
                conda_package="kallisto",
                hpc_module="kallisto",
                install_url="https://pachterlab.github.io/kallisto/",
                install_conda="conda install -c bioconda kallisto"
            )
        ]
    
    @handle_errors(exit_on_error=False)
    def run_sample(
        self,
        sample: SampleInfo,
        index_path: Path,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run Kallisto on a single sample with error handling.
        
        🆕 v2.2.0: Enhanced error handling with detailed messages.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to Kallisto index (.idx file)
        
        Returns
        -------
        Dict
            Result dictionary with success status and metrics
            
        Raises
        ------
        PipelineError
            If any step of the pipeline fails
        """
        # Validate index exists
        try:
            check_file_exists(
                index_path,
                f"Kallisto index not found: {index_path}\n"
                f"  Build index with: kallisto index -i index.idx transcripts.fa"
            )
        except FileNotFoundError as e:
            raise PipelineError(str(e)) from e
        
        # Check for single-end without fragment length
        if not sample.is_paired:
            if self.fragment_length is None or self.fragment_sd is None:
                raise PipelineError(
                    f"Single-end sample requires fragment length parameters\n"
                    f"  Sample: {sample.sample_id}\n"
                    f"  Solution: Provide --fragment-length and --fragment-sd\n"
                    f"  Example: --fragment-length 200 --fragment-sd 20\n"
                    f"  \n"
                    f"  Typical values:\n"
                    f"    RNA-seq: length=200, sd=20\n"
                    f"    Small RNA: length=50, sd=10"
                )
        
        sample_out = self.quant_dir / sample.sample_id
        sample_out.mkdir(exist_ok=True, parents=True)
        
        log_file = self.logs_dir / f"{sample.sample_id}_kallisto.log"
        
        # ============== BUILD COMMAND ==============
        cmd = [
            'kallisto', 'quant',
            '-i', str(index_path),
            '-o', str(sample_out),
            '-t', str(self.threads)
        ]
        
        # Add bootstraps
        if self.bootstraps > 0:
            cmd.extend(['-b', str(self.bootstraps)])
            cmd.extend(['--seed', str(self.seed)])
        
        # Add strand-specific options
        if self.rf_stranded:
            cmd.append('--rf-stranded')
        if self.fr_stranded:
            cmd.append('--fr-stranded')
        
        # Other options
        if self.plaintext:
            cmd.append('--plaintext')
        if self.fusion:
            cmd.append('--fusion')
        
        # Add input files
        if sample.is_paired:
            cmd.extend([str(sample.fastq_1), str(sample.fastq_2)])
            logger.info(f"  Running Kallisto (paired-end) for {sample.sample_id}")
        else:
            cmd.extend([
                '--single',
                '-l', str(self.fragment_length),
                '-s', str(self.fragment_sd),
                str(sample.fastq_1)
            ])
            logger.info(f"  Running Kallisto (single-end) for {sample.sample_id}")
        
        # Add extra args
        if self.extra_args:
            for key, value in self.extra_args.items():
                if value is not None and value != '':
                    cmd.extend([f'--{key}', str(value)])
                else:
                    cmd.append(f'--{key}')
        
        # ============== RUN KALLISTO ==============
        returncode, stdout, stderr = self.run_command(cmd, log_file)
        
        if returncode != 0:
            raise PipelineError(
                f"Kallisto failed for {sample.sample_id}\n"
                f"  Check log: {log_file}\n"
                f"  Common issues:\n"
                f"    - Index version mismatch (rebuild with same Kallisto version)\n"
                f"    - Corrupted FASTQ files\n"
                f"    - Wrong fragment length for single-end\n"
                f"    - Index built from wrong organism"
            )
        
        # ============== CHECK OUTPUT ==============
        abundance_file = sample_out / "abundance.tsv"
        if not abundance_file.exists():
            raise PipelineError(
                f"Kallisto output not found for {sample.sample_id}\n"
                f"  Expected: {abundance_file}\n"
                f"  This usually means Kallisto completed but produced no output.\n"
                f"  Check if input FASTQ files are empty or corrupted."
            )
        
        # Parse metrics from run_info.json
        metrics = self._parse_run_info(sample_out / "run_info.json")
        
        # Warn if low pseudoalignment rate
        if metrics.get('p_pseudoaligned'):
            if metrics['p_pseudoaligned'] < 50:
                logger.warning(
                    f"Low pseudoalignment rate for {sample.sample_id}: "
                    f"{metrics['p_pseudoaligned']:.1f}%"
                )
        
        return {
            'success': True,
            'sample_id': sample.sample_id,
            'counts_file': str(abundance_file),
            'quant_dir': str(sample_out),
            'log_file': str(log_file),
            'n_processed': metrics.get('n_processed'),
            'n_pseudoaligned': metrics.get('n_pseudoaligned'),
            'p_pseudoaligned': metrics.get('p_pseudoaligned')
        }
    
    def _parse_run_info(self, json_file: Path) -> Dict:
        """Parse run info from Kallisto output."""
        import json
        metrics = {}
        if json_file.exists():
            try:
                with open(json_file) as f:
                    data = json.load(f)
                    metrics['n_processed'] = data.get('n_processed')
                    metrics['n_pseudoaligned'] = data.get('n_pseudoaligned')
                    if metrics['n_processed'] and metrics['n_processed'] > 0:
                        metrics['p_pseudoaligned'] = (
                            metrics['n_pseudoaligned'] / metrics['n_processed'] * 100
                        )
            except Exception as e:
                logger.warning(f"Could not parse run_info.json: {e}")
        return metrics
    
    @handle_errors(exit_on_error=False)
    def combine_counts(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """
        Combine Kallisto abundance.tsv files into count matrix.
        
        🆕 v2.2.0: Enhanced error handling for combining counts.
        """
        if not sample_results:
            raise PipelineError("No successful samples to combine")
        
        counts_data = {}
        tpm_data = {}
        
        for result in sample_results:
            if not result.get('success'):
                logger.warning(f"Skipping failed sample: {result.get('sample_id')}")
                continue
                
            sample_id = result['sample_id']
            abundance_file = Path(result['counts_file'])
            
            if not abundance_file.exists():
                logger.warning(f"Abundance file not found for {sample_id}: {abundance_file}")
                continue
            
            try:
                # Read abundance.tsv
                df = pd.read_csv(abundance_file, sep='\t')
                
                # Extract counts (est_counts column)
                counts = df.set_index('target_id')['est_counts']
                counts_data[sample_id] = counts
                
                # Extract TPM
                tpm = df.set_index('target_id')['tpm']
                tpm_data[sample_id] = tpm
            except Exception as e:
                logger.warning(f"Could not read abundance file for {sample_id}: {e}")
                continue
        
        if not counts_data:
            raise PipelineError("No valid count data found in any samples")
        
        # Combine into matrices
        counts_df = pd.DataFrame(counts_data)
        tpm_df = pd.DataFrame(tpm_data)
        
        # Round counts to integers
        counts_df = counts_df.round().astype(int)
        
        # If tx2gene provided, aggregate to gene level
        if self.tx2gene and self.tx2gene.exists():
            logger.info("Aggregating transcript counts to gene level...")
            counts_df = self._aggregate_to_genes(counts_df, self.tx2gene)
            tpm_df = self._aggregate_to_genes(tpm_df, self.tx2gene)
        
        # Save counts
        try:
            counts_df.to_csv(output_file)
            logger.info(f"✅ Saved counts to {output_file}")
        except Exception as e:
            raise PipelineError(f"Could not save count matrix to {output_file}: {e}") from e
        
        # Store TPM for generate_tpm
        self._tpm_df = tpm_df
        
        return counts_df
    
    def _aggregate_to_genes(
        self,
        counts_df: pd.DataFrame,
        tx2gene_file: Path
    ) -> pd.DataFrame:
        """Aggregate transcript counts to gene level."""
        try:
            # Load tx2gene mapping
            tx2gene = pd.read_csv(tx2gene_file, sep='\t', header=None,
                                 names=['transcript_id', 'gene_id'])
            tx2gene_map = dict(zip(tx2gene['transcript_id'], tx2gene['gene_id']))
            
            # Map transcripts to genes
            counts_df['gene_id'] = counts_df.index.map(tx2gene_map)
            
            # Drop unmapped
            n_unmapped = counts_df['gene_id'].isna().sum()
            if n_unmapped > 0:
                logger.warning(f"Dropped {n_unmapped} unmapped transcripts")
            counts_df = counts_df.dropna(subset=['gene_id'])
            
            # Aggregate by gene
            gene_counts = counts_df.groupby('gene_id').sum()
            
            return gene_counts
        except Exception as e:
            logger.error(f"Error aggregating to genes: {e}")
            logger.warning("Returning transcript-level counts without aggregation")
            return counts_df.drop(columns=['gene_id'], errors='ignore')
    
    def generate_tpm(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """Generate TPM matrix."""
        if hasattr(self, '_tpm_df'):
            self._tpm_df.to_csv(output_file)
            logger.info(f"Saved TPM to {output_file}")
            return self._tpm_df
        
        # If not cached, recalculate
        tpm_data = {}
        for result in sample_results:
            if not result.get('success'):
                continue
            sample_id = result['sample_id']
            abundance_file = Path(result['counts_file'])
            
            try:
                df = pd.read_csv(abundance_file, sep='\t')
                tpm = df.set_index('target_id')['tpm']
                tpm_data[sample_id] = tpm
            except Exception:
                continue
        
        tpm_df = pd.DataFrame(tpm_data)
        tpm_df.to_csv(output_file)
        return tpm_df


# =============================================================================
# CLI PARAMETER DEFINITIONS (v2.2.0 Compatible)
# =============================================================================

KALLISTO_CLI_PARAMS = {
    'bootstraps': {
        'flag': '-b',
        'type': int,
        'default': 0,
        'help': 'Number of bootstrap samples (0=disabled)',
        'param_name': 'bootstraps'
    },
    'fragment-length': {
        'flag': '-l',
        'type': int,
        'default': None,
        'help': 'Fragment length mean (REQUIRED for single-end)',
        'param_name': 'fragment_length'
    },
    'fragment-sd': {
        'flag': '-s',
        'type': int,
        'default': None,
        'help': 'Fragment length SD (REQUIRED for single-end)',
        'param_name': 'fragment_sd'
    },
    'rf-stranded': {
        'flag': '--rf-stranded',
        'type': bool,
        'default': False,
        'help': 'Strand-specific, first read reverse',
        'param_name': 'rf_stranded'
    },
    'fr-stranded': {
        'flag': '--fr-stranded',
        'type': bool,
        'default': False,
        'help': 'Strand-specific, first read forward',
        'param_name': 'fr_stranded'
    },
    'seed': {
        'flag': '--seed',
        'type': int,
        'default': 42,
        'help': 'Random seed for bootstrap sampling',
        'param_name': 'seed'
    },
    'plaintext': {
        'flag': '--plaintext',
        'type': bool,
        'default': False,
        'help': 'Output plaintext instead of HDF5',
        'param_name': 'plaintext'
    },
    'tx2gene': {
        'flag': '--tx2gene',
        'type': str,
        'default': None,
        'help': 'Transcript-to-gene mapping file for gene-level aggregation',
        'param_name': 'tx2gene'
    }
}


if __name__ == '__main__':
    print(f"🦖 RAPTOR Kallisto Pipeline v{KallistoPipeline.PIPELINE_VERSION}")
    print("=" * 70)
    print("\n✅ v2.2.0 Features:")
    print("   • Full input validation with clear error messages")
    print("   • Integrated error handling (@handle_errors)")
    print("   • Compatible with CLI system")
    print("\n⚡ Ultra-fast pseudo-alignment")
    print("   Lower memory than Salmon, no BAM files")
    print("\n🔧 Available CLI parameters:")
    for param, info in KALLISTO_CLI_PARAMS.items():
        print(f"  --{param}: {info['help']} (default: {info['default']})")
    print("\n⚠️  WARNING: Single-end reads REQUIRE --fragment-length and --fragment-sd!")
