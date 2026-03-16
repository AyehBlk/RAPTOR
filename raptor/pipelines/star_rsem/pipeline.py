#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - STAR + RSEM Pipeline
FULLY INTEGRATED WITH VALIDATION & ERROR HANDLING

Gold standard pipeline for isoform-level quantification.
Uses STAR for alignment and RSEM for expectation-maximization quantification.

Features:
- Accurate isoform-level quantification
- Gene and transcript counts
- STAR alignment with BAM output
- EM-based probabilistic counting
- ✅ Full input validation with clear error messages
- ✅ Integrated error handling
- ✅ Compatible with CLI system

Requirements:
- STAR >= 2.7.0
- RSEM >= 1.3.0
- RSEM reference prepared with rsem-prepare-reference

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
    SampleInfo
)

logger = logging.getLogger(__name__)


class StarRsemPipeline(BasePipeline):
    """
    STAR alignment + RSEM quantification pipeline.
    
    Gold standard for isoform-level quantification using expectation-maximization.
    Produces both gene-level and transcript-level counts.
    
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
    memory_gb : int
        Memory in GB (default: 32, RSEM needs high memory)
    strandedness : str
        Strand specificity: 'none', 'forward', 'reverse' (default: 'none')
    paired_end : bool
        Whether samples are paired-end (default: True)
    estimate_rspd : bool
        Estimate read start position distribution (default: True)
    calc_ci : bool
        Calculate 95% credibility intervals (default: False, slow!)
    ci_memory : int
        Memory for CI calculation in MB (default: 1024, must be > 0)
    seed : int
        Random seed for reproducibility (default: 12345, must be > 0)
    fragment_length_mean : float, optional
        Fragment length mean for single-end (must be > 0 if provided)
    fragment_length_sd : float, optional
        Fragment length SD for single-end (must be > 0 if provided)
    star_sjdb_overhang : int
        STAR splice junction overhang (default: 100, must be > 0)
    
    Raises
    ------
    ValidationError
        If any input parameter fails validation
    FileNotFoundError
        If RSEM reference doesn't exist
    PipelineError
        If pipeline execution fails
    
    Examples
    --------
    >>> # Basic usage
    >>> pipeline = StarRsemPipeline(
    ...     output_dir='results/',
    ...     threads=16,
    ...     strandedness='reverse'
    ... )
    >>> # Note: index_path should be RSEM reference prefix
    >>> result = pipeline.run('samples.csv', 'rsem_ref/genome')
    >>> 
    >>> # With credibility intervals
    >>> pipeline = StarRsemPipeline(
    ...     output_dir='results/',
    ...     calc_ci=True,
    ...     ci_memory=2048
    ... )
    """
    
    PIPELINE_NAME = "star_rsem"
    PIPELINE_VERSION = "2.2.0"
    DESCRIPTION = "STAR + RSEM (gold standard, isoform-level)"
    DOCKER_IMAGE = "quay.io/biocontainers/rsem:1.3.3"
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        output_dir: Union[str, Path],
        threads: int = 8,
        memory_gb: int = 32,
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = True,
        extra_args: Optional[Dict[str, str]] = None,
        # RSEM parameters
        strandedness: str = 'none',
        paired_end: bool = True,
        estimate_rspd: bool = True,
        calc_ci: bool = False,
        ci_memory: int = 1024,
        seed: int = 12345,
        fragment_length_mean: Optional[float] = None,
        fragment_length_sd: Optional[float] = None,
        no_bam_output: bool = False,
        output_genome_bam: bool = False,
        sort_bam_by_read_name: bool = False,
        sort_bam_by_coordinate: bool = True,
        sampling_for_bam: bool = False,
        # STAR pass-through
        star_output_genome_bam: bool = True,
        star_sjdb_overhang: int = 100
    ):
        """
        Initialize STAR + RSEM pipeline with validation.
        
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
            logger.info(f"  ✓ Resources: {threads} threads, {memory_gb}GB memory")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid numeric parameter: {e}\n"
                f"  All numeric values must be positive integers.\n"
                f"  Current values: threads={threads}, memory_gb={memory_gb}"
            ) from e
        
        # CRITICAL: Validate strandedness
        valid_strandedness = ['none', 'forward', 'reverse']
        if strandedness not in valid_strandedness:
            raise ValidationError(
                f"Invalid strandedness: '{strandedness}'\n"
                f"  Must be one of: {', '.join(valid_strandedness)}\n"
                f"  \n"
                f"  Strandedness guide:\n"
                f"    'none'    = unstranded (default)\n"
                f"    'forward' = first read maps to forward strand\n"
                f"    'reverse' = first read maps to reverse strand (dUTP)\n"
                f"  \n"
                f"  Tip: 'reverse' is most common for Illumina stranded"
            )
        logger.info(f"  ✓ Strandedness: {strandedness}")
        
        # Validate seed
        try:
            validate_positive_integer(seed, 'seed')
            logger.info(f"  ✓ Random seed: {seed}")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid seed: {e}\n"
                f"  Seed must be a positive integer for reproducibility"
            ) from e
        
        # Validate ci_memory if calc_ci is enabled
        if calc_ci:
            try:
                validate_positive_integer(ci_memory, 'ci_memory')
                logger.info(f"  ✓ CI calculation enabled with {ci_memory}MB memory")
                logger.warning(
                    "Credibility interval calculation is SLOW!\n"
                    "  Expect 2-5x longer runtime per sample."
                )
            except ValidationError as e:
                raise ValidationError(
                    f"Invalid ci_memory: {e}\n"
                    f"  CI memory must be a positive integer (in MB)\n"
                    f"  Typical values: 1024-4096 MB"
                ) from e
        
        # Validate fragment length parameters (for single-end)
        if fragment_length_mean is not None:
            if fragment_length_mean <= 0:
                raise ValidationError(
                    f"fragment_length_mean must be > 0, got {fragment_length_mean}\n"
                    f"  Typical values: 150-300 for single-end RNA-seq"
                )
            logger.info(f"  ✓ Fragment length mean: {fragment_length_mean}")
        
        if fragment_length_sd is not None:
            if fragment_length_sd <= 0:
                raise ValidationError(
                    f"fragment_length_sd must be > 0, got {fragment_length_sd}\n"
                    f"  Typical values: 20-50 for single-end RNA-seq"
                )
            logger.info(f"  ✓ Fragment length SD: {fragment_length_sd}")
        
        # Warn about single-end without fragment length
        if not paired_end:
            if fragment_length_mean is None or fragment_length_sd is None:
                logger.warning(
                    "Single-end mode without fragment length parameters!\n"
                    "  RSEM can estimate these, but providing them improves accuracy.\n"
                    "  Recommend: --fragment-length-mean 200 --fragment-length-sd 30"
                )
        
        # Validate STAR sjdb overhang
        try:
            validate_positive_integer(star_sjdb_overhang, 'star_sjdb_overhang')
            logger.info(f"  ✓ STAR sjdb overhang: {star_sjdb_overhang}")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid star_sjdb_overhang: {e}\n"
                f"  Must be a positive integer (typically 100 for 100bp reads)"
            ) from e
        
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
        
        # RSEM parameters
        self.strandedness = strandedness
        self.paired_end = paired_end
        self.estimate_rspd = estimate_rspd
        self.calc_ci = calc_ci
        self.ci_memory = ci_memory
        self.seed = seed
        self.fragment_length_mean = fragment_length_mean
        self.fragment_length_sd = fragment_length_sd
        self.no_bam_output = no_bam_output
        self.output_genome_bam = output_genome_bam
        self.sort_bam_by_read_name = sort_bam_by_read_name
        self.sort_bam_by_coordinate = sort_bam_by_coordinate
        self.sampling_for_bam = sampling_for_bam
        
        # STAR parameters
        self.star_output_genome_bam = star_output_genome_bam
        self.star_sjdb_overhang = star_sjdb_overhang
        
        # Create validated directories
        self.rsem_dir = self.output_dir / "rsem"
        self.rsem_dir.mkdir(exist_ok=True, parents=True)
        self.bam_dir = self.output_dir / "bam"
        self.bam_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info("🦖 STAR + RSEM pipeline initialized")
    
    def get_dependencies(self) -> List[ToolDependency]:
        """Return RSEM and STAR dependencies."""
        return [
            ToolDependency(
                name="RSEM",
                command="rsem-calculate-expression",
                min_version="1.3.0",
                docker_image="quay.io/biocontainers/rsem:1.3.3",
                conda_package="rsem",
                hpc_module="rsem",
                install_url="https://github.com/deweylab/RSEM",
                install_conda="conda install -c bioconda rsem"
            ),
            ToolDependency(
                name="STAR",
                command="STAR",
                min_version="2.7.0",
                docker_image="quay.io/biocontainers/star:2.7.10b",
                conda_package="star",
                hpc_module="STAR",
                install_url="https://github.com/alexdobin/STAR",
                install_conda="conda install -c bioconda star"
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
        Run STAR + RSEM on a single sample with error handling.
        
        🆕 v2.2.0: Enhanced error handling with detailed messages.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to RSEM reference prefix (built with rsem-prepare-reference)
        
        Returns
        -------
        Dict
            Result dictionary with success status and file paths
            
        Raises
        ------
        PipelineError
            If any step of the pipeline fails
        """
        # Validate RSEM reference exists
        # RSEM reference consists of multiple files with this prefix
        # Check for at least one expected file
        ref_files = list(Path(str(index_path).rsplit('.', 1)[0]).parent.glob(f"{index_path.name}*"))
        if not ref_files:
            raise PipelineError(
                f"RSEM reference not found: {index_path}\n"
                f"  Build reference with:\n"
                f"    rsem-prepare-reference --gtf annotation.gtf genome.fa rsem_ref/genome\n"
                f"  \n"
                f"  Expected files like:\n"
                f"    {index_path}.grp\n"
                f"    {index_path}.ti\n"
                f"    {index_path}.transcripts.fa"
            )
        
        sample_prefix = self.rsem_dir / sample.sample_id / sample.sample_id
        sample_prefix.parent.mkdir(exist_ok=True, parents=True)
        
        log_file = self.logs_dir / f"{sample.sample_id}_rsem.log"
        
        # ============== BUILD RSEM COMMAND ==============
        cmd = [
            'rsem-calculate-expression',
            '-p', str(self.threads),
            '--star',
            '--seed', str(self.seed)
        ]
        
        # STAR-specific options
        if self.star_output_genome_bam:
            cmd.append('--star-output-genome-bam')
        
        # Strandedness (CRITICAL for accuracy!)
        if self.strandedness == 'forward':
            cmd.append('--forward-prob')
            cmd.append('1')
        elif self.strandedness == 'reverse':
            cmd.append('--forward-prob')
            cmd.append('0')
        else:  # none
            cmd.append('--forward-prob')
            cmd.append('0.5')
        
        logger.info(f"  Running STAR+RSEM for {sample.sample_id}")
        logger.info(f"    Strandedness: {self.strandedness}")
        
        # Paired-end / single-end
        if sample.is_paired:
            cmd.append('--paired-end')
            logger.info(f"    Mode: Paired-end")
        else:
            logger.info(f"    Mode: Single-end")
            # Single-end may need fragment length
            if self.fragment_length_mean:
                cmd.extend(['--fragment-length-mean', str(self.fragment_length_mean)])
            if self.fragment_length_sd:
                cmd.extend(['--fragment-length-sd', str(self.fragment_length_sd)])
        
        # Additional options
        if self.estimate_rspd:
            cmd.append('--estimate-rspd')
        
        if self.calc_ci:
            cmd.append('--calc-ci')
            cmd.extend(['--ci-memory', str(self.ci_memory)])
            logger.info(f"    CI calculation enabled (slow!)")
        
        # BAM output options
        if self.no_bam_output:
            cmd.append('--no-bam-output')
        else:
            if self.output_genome_bam:
                cmd.append('--output-genome-bam')
            if self.sort_bam_by_coordinate:
                cmd.append('--sort-bam-by-coordinate')
            elif self.sort_bam_by_read_name:
                cmd.append('--sort-bam-by-read-name')
            if self.sampling_for_bam:
                cmd.append('--sampling-for-bam')
        
        # STAR options
        cmd.extend(['--star-sjdb-overhang', str(self.star_sjdb_overhang)])
        
        # Input files
        if sample.is_paired:
            cmd.extend([str(sample.fastq_1), str(sample.fastq_2)])
        else:
            cmd.append(str(sample.fastq_1))
        
        # Reference and output prefix
        cmd.extend([str(index_path), str(sample_prefix)])
        
        # ============== RUN RSEM ==============
        returncode, stdout, stderr = self.run_command(cmd, log_file)
        
        if returncode != 0:
            raise PipelineError(
                f"RSEM failed for {sample.sample_id}\n"
                f"  Check log: {log_file}\n"
                f"  Common issues:\n"
                f"    - RSEM reference mismatch (rebuild with same version)\n"
                f"    - Corrupted FASTQ files\n"
                f"    - Insufficient memory ({self.memory_gb}GB allocated, need 32GB+)\n"
                f"    - Wrong strandedness setting\n"
                f"    - STAR index embedded in RSEM ref doesn't match STAR version"
            )
        
        # ============== CHECK OUTPUTS ==============
        gene_results = Path(f"{sample_prefix}.genes.results")
        isoform_results = Path(f"{sample_prefix}.isoforms.results")
        
        if not gene_results.exists():
            raise PipelineError(
                f"Gene results not found for {sample.sample_id}\n"
                f"  Expected: {gene_results}\n"
                f"  RSEM completed but produced no output.\n"
                f"  Check log: {log_file}"
            )
        
        if not isoform_results.exists():
            logger.warning(
                f"Isoform results not found for {sample.sample_id}: {isoform_results}"
            )
        
        # Move BAM if generated
        bam_file = None
        genome_bam = Path(f"{sample_prefix}.STAR.genome.bam")
        
        if self.keep_bam and genome_bam.exists():
            final_bam = self.bam_dir / f"{sample.sample_id}.bam"
            genome_bam.rename(final_bam)
            bam_file = str(final_bam)
            # Index BAM
            self.run_command(['samtools', 'index', str(final_bam)])
        
        # Parse stats
        stats = self._parse_rsem_stats(Path(f"{sample_prefix}.stat"))
        
        # Warn if low alignment rate
        if stats and stats.get('alignment_rate'):
            if stats['alignment_rate'] < 50:
                logger.warning(
                    f"Low alignment rate for {sample.sample_id}: "
                    f"{stats['alignment_rate']:.1f}%"
                )
        
        return {
            'success': True,
            'sample_id': sample.sample_id,
            'gene_results': str(gene_results),
            'isoform_results': str(isoform_results),
            'counts_file': str(gene_results),  # For combine_counts
            'bam_file': bam_file,
            'log_file': str(log_file),
            'stats': stats
        }
    
    def _parse_rsem_stats(self, stat_dir: Path) -> Dict:
        """Parse RSEM statistics."""
        stats = {}
        cnt_file = stat_dir / f"{stat_dir.parent.name}.cnt"
        
        if not cnt_file.exists():
            return stats
        
        try:
            with open(cnt_file) as f:
                lines = f.readlines()
                if len(lines) >= 4:
                    # Parse alignment stats
                    n_aligned = int(lines[1].strip())
                    n_total = int(lines[0].strip())
                    if n_total > 0:
                        stats['alignment_rate'] = (n_aligned / n_total) * 100
                        stats['n_aligned'] = n_aligned
                        stats['n_total'] = n_total
        except Exception as e:
            logger.warning(f"Could not parse RSEM stats: {e}")
        
        return stats
    
    @handle_errors(exit_on_error=False)
    def combine_counts(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """
        Combine RSEM gene results into count matrix.
        
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
            gene_file = Path(result['gene_results'])
            
            if not gene_file.exists():
                logger.warning(f"Gene results not found for {sample_id}: {gene_file}")
                continue
            
            try:
                # Read RSEM gene results
                df = pd.read_csv(gene_file, sep='\t')
                
                # Extract expected counts (round to integers)
                counts = df.set_index('gene_id')['expected_count'].round().astype(int)
                counts_data[sample_id] = counts
                
                # Extract TPM
                tpm = df.set_index('gene_id')['TPM']
                tpm_data[sample_id] = tpm
            except Exception as e:
                logger.warning(f"Could not read gene results for {sample_id}: {e}")
                continue
        
        if not counts_data:
            raise PipelineError("No valid count data found in any samples")
        
        # Combine into matrices
        counts_df = pd.DataFrame(counts_data)
        tpm_df = pd.DataFrame(tpm_data)
        
        # Save counts
        try:
            counts_df.to_csv(output_file)
            logger.info(f"✅ Saved gene counts to {output_file}")
        except Exception as e:
            raise PipelineError(f"Could not save count matrix to {output_file}: {e}") from e
        
        # Store TPM for later use
        self._tpm_df = tpm_df
        
        return counts_df
    
    def generate_tpm(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """Generate TPM matrix from RSEM results."""
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
            gene_file = Path(result['gene_results'])
            
            try:
                df = pd.read_csv(gene_file, sep='\t')
                tpm = df.set_index('gene_id')['TPM']
                tpm_data[sample_id] = tpm
            except Exception:
                continue
        
        tpm_df = pd.DataFrame(tpm_data)
        tpm_df.to_csv(output_file)
        return tpm_df


# =============================================================================
# CLI PARAMETER DEFINITIONS (v2.2.0 Compatible)
# =============================================================================

STAR_RSEM_CLI_PARAMS = {
    'strandedness': {
        'flag': '--strandedness',
        'type': str,
        'default': 'none',
        'choices': ['none', 'forward', 'reverse'],
        'help': 'Strandedness: none, forward, reverse (dUTP)',
        'param_name': 'strandedness'
    },
    'estimate-rspd': {
        'flag': '--estimate-rspd',
        'type': bool,
        'default': True,
        'help': 'Estimate read start position distribution',
        'param_name': 'estimate_rspd'
    },
    'calc-ci': {
        'flag': '--calc-ci',
        'type': bool,
        'default': False,
        'help': 'Calculate 95%% credibility intervals (slow!)',
        'param_name': 'calc_ci'
    },
    'ci-memory': {
        'flag': '--ci-memory',
        'type': int,
        'default': 1024,
        'help': 'Memory for CI calculation (MB)',
        'param_name': 'ci_memory'
    },
    'seed': {
        'flag': '--seed',
        'type': int,
        'default': 12345,
        'help': 'Random seed for reproducibility',
        'param_name': 'seed'
    },
    'fragment-length-mean': {
        'flag': '--fragment-length-mean',
        'type': float,
        'default': None,
        'help': 'Fragment length mean (for single-end)',
        'param_name': 'fragment_length_mean'
    },
    'fragment-length-sd': {
        'flag': '--fragment-length-sd',
        'type': float,
        'default': None,
        'help': 'Fragment length SD (for single-end)',
        'param_name': 'fragment_length_sd'
    },
}


if __name__ == '__main__':
    print(f"🦖 RAPTOR STAR + RSEM Pipeline v{StarRsemPipeline.PIPELINE_VERSION}")
    print("=" * 70)
    print("\n✅ v2.2.0 Features:")
    print("   • Full input validation with clear error messages")
    print("   • Integrated error handling (@handle_errors)")
    print("   • Compatible with CLI system")
    print("\n⚡ Gold standard for isoform-level quantification")
    print("   EM-based probabilistic counting with STAR alignment")
    print("\n🔧 Available CLI parameters:")
    for param, info in STAR_RSEM_CLI_PARAMS.items():
        print(f"  --{param}: {info['help']} (default: {info['default']})")
    print("\n⚠️  High memory requirement: ~32GB")
