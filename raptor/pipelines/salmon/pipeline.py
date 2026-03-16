#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Salmon Pipeline
FULLY INTEGRATED WITH VALIDATION & ERROR HANDLING

Production pipeline using Salmon for pseudo-alignment quantification.
Fastest method, no BAM files produced.

Features:
- Pseudo-alignment (no genome alignment needed)
- Transcript-level quantification with gene aggregation
- Optional bootstrap for uncertainty estimation
- GC bias and sequence bias correction
- ✅ Full input validation with clear error messages
- ✅ Integrated error handling
- ✅ Compatible with CLI system

Requirements:
- Salmon >= 1.9.0
- Salmon index (built from transcriptome)

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
    validate_probability,
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


class SalmonPipeline(BasePipeline):
    """
    Salmon pseudo-alignment pipeline.
    
    Salmon uses quasi-mapping for ultra-fast transcript quantification.
    No BAM files are produced, making it the fastest option for gene-level DE.
    
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
    library_type : str
        Library type: 'A' (auto), 'ISF', 'ISR', 'IU', etc. (default: 'A')
    bootstraps : int
        Number of bootstraps for uncertainty estimation (default: 0, must be >= 0)
    gc_bias : bool
        Enable GC bias correction (default: True)
    seq_bias : bool
        Enable sequence bias correction (default: True)
    validate_mappings : bool
        Enable selective alignment for better accuracy (default: True)
    min_score_fraction : float
        Minimum alignment score as fraction of max (default: 0.65, range: 0-1)
    fragment_length : int, optional
        Fragment length mean for single-end (must be > 0 if provided)
    fragment_sd : int, optional
        Fragment length SD for single-end (must be > 0 if provided)
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
    >>> # Basic usage
    >>> pipeline = SalmonPipeline(output_dir='results/', threads=16)
    >>> result = pipeline.run('samples.csv', 'salmon_index/')
    >>> 
    >>> # With bootstraps for sleuth
    >>> pipeline = SalmonPipeline(
    ...     output_dir='results/',
    ...     threads=16,
    ...     bootstraps=100
    ... )
    >>> 
    >>> # Docker execution
    >>> pipeline = SalmonPipeline(
    ...     output_dir='results/',
    ...     use_docker=True,
    ...     docker_image='combinelab/salmon:1.10.0'
    ... )
    """
    
    PIPELINE_NAME = "salmon"
    PIPELINE_VERSION = "2.2.0"
    DESCRIPTION = "Salmon pseudo-alignment (fastest, no BAM)"
    DOCKER_IMAGE = "combinelab/salmon:1.10.0"
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        output_dir: Union[str, Path],
        threads: int = 8,
        memory_gb: int = 32,
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = False,  # Salmon doesn't produce BAM
        extra_args: Optional[Dict[str, str]] = None,
        # Salmon-specific parameters
        library_type: str = 'A',
        bootstraps: int = 0,
        gc_bias: bool = True,
        seq_bias: bool = True,
        pos_bias: bool = False,
        validate_mappings: bool = True,
        min_score_fraction: float = 0.65,
        consensus_slack: float = 0.35,
        range_factorization_bins: int = 4,
        num_gibbs_samples: int = 0,
        incompatible_prior: float = 0.0,
        write_unmapped_names: bool = False,
        no_length_correction: bool = False,
        no_effective_length_correction: bool = False,
        use_vbopt: bool = False,
        # Single-end specific
        fragment_length: Optional[int] = None,
        fragment_sd: Optional[int] = None,
        # Gene aggregation
        tx2gene: Optional[Union[str, Path]] = None
    ):
        """
        Initialize Salmon pipeline with validation.
        
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
        
        # Validate library type
        valid_library_types = [
            'A', 'IU', 'ISF', 'ISR', 'MU', 'MSF', 'MSR', 
            'OU', 'OSF', 'OSR', 'SF', 'SR', 'U'
        ]
        if library_type not in valid_library_types:
            raise ValidationError(
                f"Invalid library_type: '{library_type}'\n"
                f"  Must be one of: {', '.join(valid_library_types)}\n"
                f"  Common options:\n"
                f"    'A' = auto-detect (recommended)\n"
                f"    'IU' = inward, unstranded\n"
                f"    'ISF' = inward, stranded forward\n"
                f"    'ISR' = inward, stranded reverse (dUTP)"
            )
        
        # Validate bootstraps (can be 0)
        if bootstraps < 0:
            raise ValidationError(
                f"bootstraps must be >= 0, got {bootstraps}\n"
                f"  Use 0 to disable bootstraps (faster)\n"
                f"  Use 30-100 for standard uncertainty estimation\n"
                f"  Use 200+ for publication-quality"
            )
        
        # Validate min_score_fraction
        try:
            validate_probability(min_score_fraction, 'min_score_fraction')
            logger.info(f"  ✓ min_score_fraction: {min_score_fraction}")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid min_score_fraction: {e}\n"
                f"  Must be between 0 and 1 (typically 0.5-0.8)\n"
                f"  Default: 0.65"
            ) from e
        
        # Validate consensus_slack
        try:
            validate_probability(consensus_slack, 'consensus_slack')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid consensus_slack: {e}\n"
                f"  Must be between 0 and 1 (typically 0.2-0.5)\n"
                f"  Default: 0.35"
            ) from e
        
        # Validate range_factorization_bins
        try:
            validate_positive_integer(range_factorization_bins, 'range_factorization_bins')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid range_factorization_bins: {e}\n"
                f"  Must be a positive integer (typically 4)\n"
                f"  Default: 4"
            ) from e
        
        # Validate num_gibbs_samples
        if num_gibbs_samples < 0:
            raise ValidationError(
                f"num_gibbs_samples must be >= 0, got {num_gibbs_samples}\n"
                f"  Use 0 to disable (default)\n"
                f"  Use 1000+ for more accurate quantification"
            )
        
        # Validate incompatible_prior
        if incompatible_prior < 0 or incompatible_prior > 1:
            raise ValidationError(
                f"incompatible_prior must be between 0 and 1, got {incompatible_prior}\n"
                f"  Default: 0.0 (no prior)"
            )
        
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
        
        # Validate tx2gene file if provided
        if tx2gene:
            try:
                self.tx2gene = validate_file_path(
                    tx2gene,
                    must_exist=True,
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
        
        # Salmon-specific parameters
        self.library_type = library_type
        self.bootstraps = bootstraps
        self.gc_bias = gc_bias
        self.seq_bias = seq_bias
        self.pos_bias = pos_bias
        self.validate_mappings = validate_mappings
        self.min_score_fraction = min_score_fraction
        self.consensus_slack = consensus_slack
        self.range_factorization_bins = range_factorization_bins
        self.num_gibbs_samples = num_gibbs_samples
        self.incompatible_prior = incompatible_prior
        self.write_unmapped_names = write_unmapped_names
        self.no_length_correction = no_length_correction
        self.no_effective_length_correction = no_effective_length_correction
        self.use_vbopt = use_vbopt
        self.fragment_length = fragment_length
        self.fragment_sd = fragment_sd
        
        # Create validated directories
        self.quant_dir = self.output_dir / "quant"
        self.quant_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info("🦖 Salmon pipeline initialized")
    
    def get_dependencies(self) -> List[ToolDependency]:
        """Return Salmon dependency."""
        return [
            ToolDependency(
                name="Salmon",
                command="salmon",
                min_version="1.9.0",
                docker_image="combinelab/salmon:1.10.0",
                conda_package="salmon",
                hpc_module="salmon",
                install_url="https://github.com/COMBINE-lab/salmon",
                install_conda="conda install -c bioconda salmon"
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
        Run Salmon on a single sample with error handling.
        
        🆕 v2.2.0: Enhanced error handling with detailed messages.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to Salmon index
        
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
                f"Salmon index not found: {index_path}\n"
                f"  Build index with: salmon index -t transcripts.fa -i salmon_index"
            )
        except FileNotFoundError as e:
            raise PipelineError(str(e)) from e
        
        sample_out = self.quant_dir / sample.sample_id
        sample_out.mkdir(exist_ok=True, parents=True)
        
        log_file = self.logs_dir / f"{sample.sample_id}_salmon.log"
        
        # ============== BUILD COMMAND ==============
        cmd = [
            'salmon', 'quant',
            '-i', str(index_path),
            '-l', self.library_type,
            '-o', str(sample_out),
            '-p', str(self.threads)
        ]
        
        # Add input files
        if sample.is_paired:
            cmd.extend(['-1', str(sample.fastq_1), '-2', str(sample.fastq_2)])
            logger.info(f"  Running Salmon (paired-end) for {sample.sample_id}")
        else:
            cmd.extend(['-r', str(sample.fastq_1)])
            logger.info(f"  Running Salmon (single-end) for {sample.sample_id}")
            # Single-end may need fragment length
            if self.fragment_length:
                cmd.extend(['--fldMean', str(self.fragment_length)])
            if self.fragment_sd:
                cmd.extend(['--fldSD', str(self.fragment_sd)])
        
        # Add optional flags
        if self.validate_mappings:
            cmd.append('--validateMappings')
        
        if self.gc_bias:
            cmd.append('--gcBias')
        
        if self.seq_bias:
            cmd.append('--seqBias')
        
        if self.pos_bias:
            cmd.append('--posBias')
        
        if self.bootstraps > 0:
            cmd.extend(['--numBootstraps', str(self.bootstraps)])
        
        if self.num_gibbs_samples > 0:
            cmd.extend(['--numGibbsSamples', str(self.num_gibbs_samples)])
        
        if self.incompatible_prior > 0:
            cmd.extend(['--incompatPrior', str(self.incompatible_prior)])
        
        if self.min_score_fraction != 0.65:
            cmd.extend(['--minScoreFraction', str(self.min_score_fraction)])
        
        if self.consensus_slack != 0.35:
            cmd.extend(['--consensusSlack', str(self.consensus_slack)])
        
        if self.range_factorization_bins != 4:
            cmd.extend(['--rangeFactorizationBins', str(self.range_factorization_bins)])
        
        if self.write_unmapped_names:
            cmd.append('--writeUnmappedNames')
        
        if self.no_length_correction:
            cmd.append('--noLengthCorrection')
        
        if self.no_effective_length_correction:
            cmd.append('--noEffectiveLengthCorrection')
        
        if self.use_vbopt:
            cmd.append('--useVBOpt')
        
        # Add extra args
        if self.extra_args:
            for key, value in self.extra_args.items():
                if value is not None and value != '':
                    cmd.extend([f'--{key}', str(value)])
                else:
                    cmd.append(f'--{key}')
        
        # ============== RUN SALMON ==============
        returncode, stdout, stderr = self.run_command(cmd, log_file)
        
        if returncode != 0:
            raise PipelineError(
                f"Salmon failed for {sample.sample_id}\n"
                f"  Check log: {log_file}\n"
                f"  Common issues:\n"
                f"    - Index version mismatch (rebuild with same Salmon version)\n"
                f"    - Corrupted FASTQ files\n"
                f"    - Index built from wrong organism\n"
                f"    - Wrong library type (try -l A for auto-detect)\n"
                f"    - Insufficient memory ({self.memory_gb}GB allocated)"
            )
        
        # ============== CHECK OUTPUT ==============
        quant_file = sample_out / "quant.sf"
        if not quant_file.exists():
            raise PipelineError(
                f"Salmon output not found for {sample.sample_id}\n"
                f"  Expected: {quant_file}\n"
                f"  This usually means Salmon completed but produced no output.\n"
                f"  Check if input FASTQ files are empty or corrupted."
            )
        
        # Parse mapping rate from logs
        mapping_rate = self._parse_mapping_rate(sample_out / "logs" / "salmon_quant.log")
        
        # Warn if low mapping rate
        if mapping_rate and mapping_rate < 50:
            logger.warning(
                f"Low mapping rate for {sample.sample_id}: {mapping_rate:.1f}%"
            )
        
        return {
            'success': True,
            'sample_id': sample.sample_id,
            'counts_file': str(quant_file),
            'quant_dir': str(sample_out),
            'log_file': str(log_file),
            'mapping_rate': mapping_rate
        }
    
    def _parse_mapping_rate(self, log_file: Path) -> Optional[float]:
        """Parse mapping rate from Salmon log."""
        if not log_file.exists():
            return None
        
        try:
            with open(log_file) as f:
                for line in f:
                    if "Mapping rate" in line or "mapping rate" in line:
                        # Extract percentage
                        import re
                        match = re.search(r'(\d+\.?\d*)%', line)
                        if match:
                            return float(match.group(1))
        except Exception as e:
            logger.warning(f"Could not parse mapping rate: {e}")
        
        return None
    
    @handle_errors(exit_on_error=False)
    def combine_counts(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """
        Combine Salmon quant.sf files into count matrix.
        
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
            quant_file = Path(result['counts_file'])
            
            if not quant_file.exists():
                logger.warning(f"Quant file not found for {sample_id}: {quant_file}")
                continue
            
            try:
                # Read quant.sf
                df = pd.read_csv(quant_file, sep='\t')
                
                # Extract counts (NumReads column)
                counts = df.set_index('Name')['NumReads']
                counts_data[sample_id] = counts
                
                # Extract TPM
                tpm = df.set_index('Name')['TPM']
                tpm_data[sample_id] = tpm
            except Exception as e:
                logger.warning(f"Could not read quant file for {sample_id}: {e}")
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
            quant_file = Path(result['counts_file'])
            
            try:
                df = pd.read_csv(quant_file, sep='\t')
                tpm = df.set_index('Name')['TPM']
                tpm_data[sample_id] = tpm
            except Exception:
                continue
        
        tpm_df = pd.DataFrame(tpm_data)
        tpm_df.to_csv(output_file)
        return tpm_df


# =============================================================================
# CLI PARAMETER DEFINITIONS (v2.2.0 Compatible)
# =============================================================================

SALMON_CLI_PARAMS = {
    'library-type': {
        'flag': '-l',
        'type': str,
        'default': 'A',
        'choices': ['A', 'IU', 'ISF', 'ISR', 'MU', 'MSF', 'MSR', 'OU', 'OSF', 'OSR', 'SF', 'SR', 'U'],
        'help': 'Library type (A=auto-detect, ISR=dUTP)',
        'param_name': 'library_type'
    },
    'bootstraps': {
        'flag': '-b',
        'type': int,
        'default': 0,
        'help': 'Number of bootstrap samples (0=disabled, 100+ recommended)',
        'param_name': 'bootstraps'
    },
    'gc-bias': {
        'flag': '--gcBias',
        'type': bool,
        'default': True,
        'help': 'Enable GC bias correction',
        'param_name': 'gc_bias'
    },
    'seq-bias': {
        'flag': '--seqBias',
        'type': bool,
        'default': True,
        'help': 'Enable sequence bias correction',
        'param_name': 'seq_bias'
    },
    'pos-bias': {
        'flag': '--posBias',
        'type': bool,
        'default': False,
        'help': 'Enable positional bias correction',
        'param_name': 'pos_bias'
    },
    'validate-mappings': {
        'flag': '--validateMappings',
        'type': bool,
        'default': True,
        'help': 'Enable selective alignment (more accurate)',
        'param_name': 'validate_mappings'
    },
    'min-score-fraction': {
        'flag': '--minScoreFraction',
        'type': float,
        'default': 0.65,
        'help': 'Minimum alignment score fraction (0-1)',
        'param_name': 'min_score_fraction'
    },
    'fragment-length': {
        'flag': '--fldMean',
        'type': int,
        'default': None,
        'help': 'Fragment length mean (for single-end)',
        'param_name': 'fragment_length'
    },
    'fragment-sd': {
        'flag': '--fldSD',
        'type': int,
        'default': None,
        'help': 'Fragment length SD (for single-end)',
        'param_name': 'fragment_sd'
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
    print(f"🦖 RAPTOR Salmon Pipeline v{SalmonPipeline.PIPELINE_VERSION}")
    print("=" * 70)
    print("\n✅ v2.2.0 Features:")
    print("   • Full input validation with clear error messages")
    print("   • Integrated error handling (@handle_errors)")
    print("   • Compatible with CLI system")
    print("\n⚡ Fastest pseudo-alignment")
    print("   GC bias correction, no BAM files")
    print("\n🔧 Available CLI parameters:")
    for param, info in SALMON_CLI_PARAMS.items():
        print(f"  --{param}: {info['help']} (default: {info['default']})")
