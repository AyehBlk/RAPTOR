#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - STAR + Salmon Pipeline
FULLY INTEGRATED WITH VALIDATION & ERROR HANDLING

Hybrid pipeline combining STAR alignment with Salmon quantification.
Produces BAM files AND supports bootstraps for uncertainty estimation.

Features:
- Full genome alignment with STAR (BAM output)
- Salmon quantification from BAM (bootstraps supported)
- Best of both worlds: BAM + bootstrap uncertainty
- ✅ Full input validation with clear error messages
- ✅ Integrated error handling
- ✅ Compatible with CLI system

Requirements:
- STAR >= 2.7.0
- Salmon >= 1.9.0
- STAR genome index + Salmon transcriptome index (BOTH required!)

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


class StarSalmonPipeline(BasePipeline):
    """
    STAR alignment + Salmon quantification pipeline.
    
    Hybrid pipeline combining STAR alignment (for BAM) with Salmon quantification 
    (for bootstraps). Use when you need both BAM files AND uncertainty estimates.
    
    🆕 v2.2.0 Features:
    - ✅ Full input validation
    - ✅ Clear error messages
    - ✅ Graceful error handling
    - ✅ Validates BOTH indexes (STAR + Salmon)
    
    Parameters
    ----------
    output_dir : str or Path
        Output directory
    salmon_index : str or Path
        Path to Salmon transcriptome index (REQUIRED in addition to STAR index!)
    threads : int
        Number of threads (default: 8, must be > 0)
    memory_gb : int
        Memory in GB (default: 32, STAR needs high memory)
    two_pass_mode : str
        STAR two-pass mode: 'None', 'Basic' (default: 'Basic')
    library_type : str
        Salmon library type (default: 'A' for auto-detect)
    bootstraps : int
        Number of bootstraps (default: 100, must be >= 0)
    gc_bias : bool
        Enable GC bias correction (default: True)
    seq_bias : bool
        Enable sequence bias correction (default: True)
    
    Raises
    ------
    ValidationError
        If any input parameter fails validation
    FileNotFoundError
        If STAR index or Salmon index doesn't exist
    PipelineError
        If pipeline execution fails
    
    Examples
    --------
    >>> # Basic usage - note BOTH indexes required!
    >>> pipeline = StarSalmonPipeline(
    ...     output_dir='results/',
    ...     salmon_index='salmon_index/',  # REQUIRED!
    ...     threads=16,
    ...     bootstraps=100
    ... )
    >>> result = pipeline.run('samples.csv', 'star_index/')
    >>> 
    >>> # With Docker
    >>> pipeline = StarSalmonPipeline(
    ...     output_dir='results/',
    ...     salmon_index='salmon_index/',
    ...     use_docker=True
    ... )
    """
    
    PIPELINE_NAME = "star_salmon"
    PIPELINE_VERSION = "2.2.0"
    DESCRIPTION = "STAR + Salmon (BAM + bootstraps)"
    DOCKER_IMAGE = "combinelab/salmon:1.10.0"
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        output_dir: Union[str, Path],
        salmon_index: Union[str, Path],  # REQUIRED!
        threads: int = 8,
        memory_gb: int = 32,
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = True,
        extra_args: Optional[Dict[str, str]] = None,
        # STAR parameters
        two_pass_mode: str = 'Basic',
        out_filter_multimap_nmax: int = 20,
        out_filter_mismatch_nmax: int = 10,
        align_intron_max: int = 1000000,
        # Salmon parameters
        library_type: str = 'A',
        bootstraps: int = 100,
        gc_bias: bool = True,
        seq_bias: bool = True,
        num_gibbs_samples: int = 0
    ):
        """
        Initialize STAR + Salmon hybrid pipeline with validation.
        
        🆕 v2.2.0: All inputs validated, including BOTH indexes.
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
        
        # CRITICAL: Validate Salmon index (REQUIRED for this hybrid pipeline!)
        try:
            self.salmon_index = validate_directory_path(
                salmon_index,
                must_exist=True,
                description="Salmon transcriptome index"
            )
            logger.info(f"  ✓ Salmon index: {self.salmon_index}")
        except Exception as e:
            raise ValidationError(
                f"Salmon index not found: {salmon_index}\n"
                f"  This hybrid pipeline requires BOTH:\n"
                f"    1. STAR genome index (provided with --index)\n"
                f"    2. Salmon transcriptome index (provide with --salmon-index)\n"
                f"  \n"
                f"  Build Salmon index with:\n"
                f"    salmon index -t transcripts.fa -i salmon_index\n"
                f"  \n"
                f"  Common locations:\n"
                f"    - /ref/salmon/\n"
                f"    - /data/indexes/salmon/"
            ) from e
        
        # Validate numeric parameters
        try:
            validate_positive_integer(threads, 'threads')
            validate_positive_integer(memory_gb, 'memory_gb')
            logger.info(f"  ✓ Resources: {threads} threads, {memory_gb}GB memory")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid numeric parameter: {e}\n"
                f"  All numeric values must be positive integers."
            ) from e
        
        # Validate STAR two-pass mode
        valid_two_pass = ['None', 'Basic']
        if two_pass_mode not in valid_two_pass:
            raise ValidationError(
                f"Invalid two_pass_mode: '{two_pass_mode}'\n"
                f"  Must be one of: {', '.join(valid_two_pass)}\n"
                f"  Recommended: 'Basic' for better sensitivity"
            )
        logger.info(f"  ✓ STAR two-pass mode: {two_pass_mode}")
        
        # Validate STAR numeric parameters
        try:
            validate_positive_integer(out_filter_multimap_nmax, 'out_filter_multimap_nmax')
            if out_filter_mismatch_nmax < 0:
                raise ValidationError("out_filter_mismatch_nmax must be >= 0")
            validate_positive_integer(align_intron_max, 'align_intron_max')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid STAR parameter: {e}"
            ) from e
        
        # CRITICAL: Validate Salmon library type (reuse from Salmon pipeline)
        valid_library_types = [
            'A',      # Auto-detect (recommended)
            'IU',     # Inward, unstranded
            'ISF',    # Inward, stranded forward
            'ISR',    # Inward, stranded reverse (dUTP - most common)
            'MU',     # Matching, unstranded
            'MSF',    # Matching, stranded forward
            'MSR',    # Matching, stranded reverse
            'OU',     # Outward, unstranded
            'OSF',    # Outward, stranded forward
            'OSR',    # Outward, stranded reverse
            'SF',     # Stranded forward (single-end)
            'SR',     # Stranded reverse (single-end)
            'U'       # Unstranded (single-end)
        ]
        
        if library_type not in valid_library_types:
            raise ValidationError(
                f"Invalid library_type: '{library_type}'\n"
                f"  Must be one of: {', '.join(valid_library_types)}\n"
                f"  \n"
                f"  Common options:\n"
                f"    'A'   = auto-detect (recommended)\n"
                f"    'IU'  = inward, unstranded\n"
                f"    'ISF' = inward, stranded forward\n"
                f"    'ISR' = inward, stranded reverse (dUTP - most common)\n"
                f"  \n"
                f"  See: https://salmon.readthedocs.io/en/latest/library_type.html"
            )
        logger.info(f"  ✓ Salmon library type: {library_type}")
        
        # Validate bootstraps
        if bootstraps < 0:
            raise ValidationError(
                f"bootstraps must be >= 0, got {bootstraps}\n"
                f"  Use 0 for no bootstraps (fastest)\n"
                f"  Use 30+ for uncertainty estimation\n"
                f"  Use 100+ for differential expression (recommended)"
            )
        logger.info(f"  ✓ Bootstraps: {bootstraps}")
        
        # Validate Gibbs samples
        if num_gibbs_samples < 0:
            raise ValidationError(
                f"num_gibbs_samples must be >= 0, got {num_gibbs_samples}"
            )
        
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
        
        # STAR parameters
        self.two_pass_mode = two_pass_mode
        self.out_filter_multimap_nmax = out_filter_multimap_nmax
        self.out_filter_mismatch_nmax = out_filter_mismatch_nmax
        self.align_intron_max = align_intron_max
        
        # Salmon parameters
        self.library_type = library_type
        self.bootstraps = bootstraps
        self.gc_bias = gc_bias
        self.seq_bias = seq_bias
        self.num_gibbs_samples = num_gibbs_samples
        
        # Create validated directories
        self.star_dir = self.output_dir / "star"
        self.star_dir.mkdir(exist_ok=True, parents=True)
        self.salmon_dir = self.output_dir / "salmon"
        self.salmon_dir.mkdir(exist_ok=True, parents=True)
        self.bam_dir = self.output_dir / "bam"
        self.bam_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info("🦖 STAR + Salmon hybrid pipeline initialized")
    
    def get_dependencies(self) -> List[ToolDependency]:
        """Return STAR, Salmon, and samtools dependencies."""
        return [
            ToolDependency(
                name="STAR",
                command="STAR",
                min_version="2.7.0",
                docker_image="quay.io/biocontainers/star:2.7.10b",
                conda_package="star",
                hpc_module="STAR",
                install_url="https://github.com/alexdobin/STAR",
                install_conda="conda install -c bioconda star"
            ),
            ToolDependency(
                name="Salmon",
                command="salmon",
                min_version="1.9.0",
                docker_image="combinelab/salmon:1.10.0",
                conda_package="salmon",
                hpc_module="salmon",
                install_url="https://github.com/COMBINE-lab/salmon",
                install_conda="conda install -c bioconda salmon"
            ),
            ToolDependency(
                name="samtools",
                command="samtools",
                min_version="1.10",
                docker_image="quay.io/biocontainers/samtools:1.17",
                conda_package="samtools",
                hpc_module="samtools",
                install_url="http://www.htslib.org/",
                install_conda="conda install -c bioconda samtools"
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
        Run STAR + Salmon on a single sample with error handling.
        
        🆕 v2.2.0: Enhanced two-stage validation and error messages.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to STAR genome index directory
        
        Returns
        -------
        Dict
            Result dictionary with success status and file paths
            
        Raises
        ------
        PipelineError
            If any step of the pipeline fails
        
        Notes
        -----
        This is a TWO-STAGE pipeline:
        1. STAR alignment → genome BAM + transcriptome BAM
        2. Salmon quantification from transcriptome BAM
        """
        # Validate STAR index
        try:
            check_file_exists(
                index_path,
                f"STAR index not found: {index_path}\n"
                f"  Build index with: STAR --runMode genomeGenerate ..."
            )
        except FileNotFoundError as e:
            raise PipelineError(str(e)) from e
        
        sample_star = self.star_dir / sample.sample_id
        sample_star.mkdir(exist_ok=True, parents=True)
        sample_salmon = self.salmon_dir / sample.sample_id
        sample_salmon.mkdir(exist_ok=True, parents=True)
        
        star_log = self.logs_dir / f"{sample.sample_id}_star.log"
        salmon_log = self.logs_dir / f"{sample.sample_id}_salmon.log"
        
        # ============== STEP 1: STAR Alignment ==============
        logger.info(f"  [1/2] Running STAR alignment for {sample.sample_id}")
        
        cmd_star = [
            'STAR',
            '--runThreadN', str(self.threads),
            '--genomeDir', str(index_path),
            '--outFileNamePrefix', str(sample_star) + '/',
            '--twopassMode', self.two_pass_mode,
            '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--outSAMunmapped', 'Within',
            '--quantMode', 'TranscriptomeSAM',  # CRITICAL: Output transcriptome BAM for Salmon!
            '--outFilterMultimapNmax', str(self.out_filter_multimap_nmax),
            '--outFilterMismatchNmax', str(self.out_filter_mismatch_nmax),
            '--alignIntronMax', str(self.align_intron_max)
        ]
        
        # Input files
        if sample.is_paired:
            cmd_star.extend(['--readFilesIn', str(sample.fastq_1), str(sample.fastq_2)])
            logger.info(f"    Paired-end mode")
        else:
            cmd_star.extend(['--readFilesIn', str(sample.fastq_1)])
            logger.info(f"    Single-end mode")
        
        # Handle gzipped files
        if str(sample.fastq_1).endswith('.gz'):
            cmd_star.extend(['--readFilesCommand', 'zcat'])
        
        # Run STAR
        returncode, stdout, stderr = self.run_command(cmd_star, star_log)
        
        if returncode != 0:
            raise PipelineError(
                f"STAR alignment failed for {sample.sample_id}\n"
                f"  Check log: {star_log}\n"
                f"  Common issues:\n"
                f"    - Index version mismatch\n"
                f"    - Corrupted FASTQ files\n"
                f"    - Insufficient memory ({self.memory_gb}GB allocated, need 32GB+)\n"
                f"    - Disk space full"
            )
        
        # Find transcriptome BAM (CRITICAL for Salmon!)
        tx_bam = sample_star / "Aligned.toTranscriptome.out.bam"
        if not tx_bam.exists():
            raise PipelineError(
                f"Transcriptome BAM not found for {sample.sample_id}\n"
                f"  Expected: {tx_bam}\n"
                f"  STAR must output transcriptome BAM for Salmon.\n"
                f"  Check: --quantMode TranscriptomeSAM was used.\n"
                f"  Check log: {star_log}"
            )
        
        # Move genome BAM to bam directory
        genome_bam = sample_star / "Aligned.sortedByCoord.out.bam"
        final_bam = None
        if self.keep_bam and genome_bam.exists():
            final_bam = self.bam_dir / f"{sample.sample_id}.bam"
            genome_bam.rename(final_bam)
            # Index BAM
            self.run_command(['samtools', 'index', str(final_bam)])
        
        # Parse STAR stats
        star_stats = self._parse_star_log(sample_star / "Log.final.out")
        
        # Warn if low mapping rate
        if star_stats.get('uniquely_mapped_reads_%'):
            mapping_rate = star_stats['uniquely_mapped_reads_%']
            if mapping_rate < 50:
                logger.warning(
                    f"Low STAR mapping rate for {sample.sample_id}: {mapping_rate:.1f}%"
                )
        
        # ============== STEP 2: Salmon Quantification ==============
        logger.info(f"  [2/2] Running Salmon quantification for {sample.sample_id}")
        
        cmd_salmon = [
            'salmon', 'quant',
            '-i', str(self.salmon_index),
            '-l', self.library_type,
            '-a', str(tx_bam),  # Read from transcriptome BAM!
            '-o', str(sample_salmon),
            '-p', str(self.threads)
        ]
        
        # Bootstraps
        if self.bootstraps > 0:
            cmd_salmon.extend(['--numBootstraps', str(self.bootstraps)])
            logger.info(f"    Bootstraps: {self.bootstraps}")
        
        # Gibbs samples
        if self.num_gibbs_samples > 0:
            cmd_salmon.extend(['--numGibbsSamples', str(self.num_gibbs_samples)])
        
        # Bias correction
        if self.gc_bias:
            cmd_salmon.append('--gcBias')
        if self.seq_bias:
            cmd_salmon.append('--seqBias')
        
        # Run Salmon
        returncode, stdout, stderr = self.run_command(cmd_salmon, salmon_log)
        
        if returncode != 0:
            raise PipelineError(
                f"Salmon quantification failed for {sample.sample_id}\n"
                f"  Check log: {salmon_log}\n"
                f"  Common issues:\n"
                f"    - Salmon index doesn't match transcripts in STAR index\n"
                f"    - Wrong library type (try -l A for auto-detect)\n"
                f"    - Corrupted transcriptome BAM\n"
                f"    - Salmon version mismatch with index"
            )
        
        # Check Salmon output
        quant_file = sample_salmon / "quant.sf"
        if not quant_file.exists():
            raise PipelineError(
                f"Salmon output not found for {sample.sample_id}\n"
                f"  Expected: {quant_file}\n"
                f"  Salmon completed but produced no output.\n"
                f"  Check log: {salmon_log}"
            )
        
        return {
            'success': True,
            'sample_id': sample.sample_id,
            'counts_file': str(quant_file),
            'quant_dir': str(sample_salmon),
            'bam_file': str(final_bam) if final_bam else None,
            'log_files': {
                'star': str(star_log),
                'salmon': str(salmon_log)
            },
            'star_stats': star_stats,
            'mapping_rate': star_stats.get('uniquely_mapped_reads_%')
        }
    
    def _parse_star_log(self, log_file: Path) -> Dict:
        """Parse STAR Log.final.out file."""
        stats = {}
        if not log_file.exists():
            return stats
        
        try:
            with open(log_file) as f:
                for line in f:
                    if 'Uniquely mapped reads %' in line:
                        stats['uniquely_mapped_reads_%'] = float(line.split('|')[1].strip().rstrip('%'))
                    elif 'Number of input reads' in line:
                        stats['n_input_reads'] = int(line.split('|')[1].strip())
        except Exception as e:
            logger.warning(f"Could not parse STAR log: {e}")
        
        return stats
    
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
                df = pd.read_csv(quant_file, sep='\t')
                counts = df.set_index('Name')['NumReads']
                counts_data[sample_id] = counts
                
                tpm = df.set_index('Name')['TPM']
                tpm_data[sample_id] = tpm
            except Exception as e:
                logger.warning(f"Could not read quant file for {sample_id}: {e}")
                continue
        
        if not counts_data:
            raise PipelineError("No valid count data found in any samples")
        
        # Combine into matrices
        counts_df = pd.DataFrame(counts_data)
        counts_df = counts_df.round().astype(int)
        
        # Save counts
        try:
            counts_df.to_csv(output_file)
            logger.info(f"✅ Saved counts to {output_file}")
        except Exception as e:
            raise PipelineError(f"Could not save count matrix to {output_file}: {e}") from e
        
        # Store TPM for later use
        self._tpm_df = pd.DataFrame(tpm_data)
        
        return counts_df
    
    def generate_tpm(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """Generate TPM matrix from Salmon results."""
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

STAR_SALMON_CLI_PARAMS = {
    'salmon-index': {
        'flag': '--salmon-index',
        'type': str,
        'required': True,
        'help': 'Salmon transcriptome index (REQUIRED!)',
        'param_name': 'salmon_index'
    },
    'two-pass-mode': {
        'flag': '--two-pass-mode',
        'type': str,
        'default': 'Basic',
        'choices': ['None', 'Basic'],
        'help': 'STAR two-pass mode',
        'param_name': 'two_pass_mode'
    },
    'library-type': {
        'flag': '--library-type',
        'type': str,
        'default': 'A',
        'help': 'Salmon library type (A=auto)',
        'param_name': 'library_type'
    },
    'bootstraps': {
        'flag': '--bootstraps',
        'type': int,
        'default': 100,
        'help': 'Number of bootstraps for uncertainty',
        'param_name': 'bootstraps'
    },
    'gc-bias': {
        'flag': '--gc-bias/--no-gc-bias',
        'type': bool,
        'default': True,
        'help': 'Enable GC bias correction',
        'param_name': 'gc_bias'
    },
}


if __name__ == '__main__':
    print(f"🦖 RAPTOR STAR + Salmon Pipeline v{StarSalmonPipeline.PIPELINE_VERSION}")
    print("=" * 70)
    print("\n✅ v2.2.0 Features:")
    print("   • Full input validation with clear error messages")
    print("   • Integrated error handling (@handle_errors)")
    print("   • Compatible with CLI system")
    print("\n🔥 Best of both worlds: BAM files + bootstrap uncertainty")
    print("   STAR alignment → Salmon quantification")
    print("\n⚠️  Requires TWO indexes:")
    print("   --index PATH         : STAR genome index")
    print("   --salmon-index PATH  : Salmon transcriptome index")
    print("\n🔧 Available CLI parameters:")
    for param, info in STAR_SALMON_CLI_PARAMS.items():
        required = " (REQUIRED!)" if info.get('required') else ""
        print(f"  --{param}: {info['help']}{required}")
