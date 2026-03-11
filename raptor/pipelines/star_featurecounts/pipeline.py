#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - STAR + featureCounts Pipeline
FULLY INTEGRATED WITH VALIDATION & ERROR HANDLING

Production pipeline using STAR alignment + featureCounts for quantification.
Standard alignment-based workflow, produces BAM files.

Features:
- Full genome alignment with STAR
- Gene-level quantification with featureCounts
- BAM files for downstream analysis (IGV, variants, etc.)
- Two-pass mode for better sensitivity
- ✅ Full input validation with clear error messages
- ✅ Integrated error handling
- ✅ Compatible with CLI system

Requirements:
- STAR >= 2.7.0
- featureCounts (subread >= 2.0.0)
- STAR index (built from genome + GTF)
- GTF annotation file

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import pandas as pd
import os

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
    SampleInfo
)

logger = logging.getLogger(__name__)


class StarFeatureCountsPipeline(BasePipeline):
    """
    STAR alignment + featureCounts quantification pipeline.
    
    Standard alignment-based workflow producing BAM files and gene counts.
    Recommended when BAM files are needed for downstream analysis.
    
    🆕 v2.2.0 Features:
    - ✅ Full input validation
    - ✅ Clear error messages
    - ✅ Graceful error handling
    - ✅ Type checking
    
    Parameters
    ----------
    output_dir : str or Path
        Output directory
    gtf : str or Path
        GTF annotation file (REQUIRED, must exist)
    threads : int
        Number of threads (default: 8, must be > 0)
    memory_gb : int
        Memory in GB (default: 32, STAR needs high memory)
    two_pass_mode : str
        STAR two-pass mode: 'None', 'Basic' (default: 'Basic')
    out_sam_type : str
        Output BAM type (default: 'BAM SortedByCoordinate')
    quant_mode : str
        STAR quantification mode (default: 'GeneCounts')
    out_filter_multimap_nmax : int
        Max number of loci read can map to (default: 20, must be > 0)
    out_filter_mismatch_nmax : int
        Max number of mismatches (default: 10, must be >= 0)
    align_intron_min : int
        Minimum intron length (default: 20, must be > 0)
    align_intron_max : int
        Maximum intron length (default: 1000000, must be > align_intron_min)
    feature_type : str
        featureCounts feature type (default: 'exon')
    attribute_type : str
        featureCounts attribute type (default: 'gene_id')
    strand_specific : int
        Strandedness: 0=unstranded, 1=stranded, 2=reverse (default: 0)
    
    Raises
    ------
    ValidationError
        If any input parameter fails validation
    FileNotFoundError
        If GTF file or index doesn't exist
    PipelineError
        If pipeline execution fails
    
    Examples
    --------
    >>> # Basic usage
    >>> pipeline = StarFeatureCountsPipeline(
    ...     output_dir='results/',
    ...     gtf='annotation.gtf',
    ...     threads=16
    ... )
    >>> result = pipeline.run('samples.csv', 'star_index/')
    >>> 
    >>> # With Docker
    >>> pipeline = StarFeatureCountsPipeline(
    ...     output_dir='results/',
    ...     gtf='annotation.gtf',
    ...     use_docker=True,
    ...     docker_image='quay.io/biocontainers/star:2.7.10b'
    ... )
    """
    
    PIPELINE_NAME = "star_featurecounts"
    PIPELINE_VERSION = "2.2.0"
    DESCRIPTION = "STAR alignment + featureCounts (standard, produces BAM)"
    DOCKER_IMAGE = "quay.io/biocontainers/star:2.7.10b"
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        output_dir: Union[str, Path],
        gtf: Union[str, Path],  # REQUIRED
        threads: int = 8,
        memory_gb: int = 32,
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = True,
        extra_args: Optional[Dict[str, str]] = None,
        # STAR parameters
        two_pass_mode: str = 'Basic',
        out_sam_type: str = 'BAM SortedByCoordinate',
        quant_mode: str = 'GeneCounts',
        out_filter_multimap_nmax: int = 20,
        out_filter_mismatch_nmax: int = 10,
        out_filter_mismatch_nover_read_lmax: float = 0.04,
        align_intron_min: int = 20,
        align_intron_max: int = 1000000,
        align_mates_gap_max: int = 1000000,
        align_sjdb_overhang_min: int = 1,
        out_sam_unmapped: str = 'Within',
        out_sam_attributes: str = 'Standard',
        out_bam_compression: int = 1,
        limit_bam_sort_ram: int = 0,
        chimeric_segment_min: int = 12,
        chimeric_junction_overhang_min: int = 12,
        # featureCounts parameters
        feature_type: str = 'exon',
        attribute_type: str = 'gene_id',
        count_multimapping: bool = False,
        count_primary: bool = False,
        count_fraction: bool = False,
        min_overlap: int = 1,
        strand_specific: int = 0,
        count_paired_end: bool = True,
        require_both_ends_mapped: bool = False,
        count_chimeric: bool = False
    ):
        """
        Initialize STAR + featureCounts pipeline with validation.
        
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
        
        # Validate GTF file (REQUIRED)
        try:
            self.gtf = validate_file_path(
                gtf,
                must_exist=True,
                extension='.gtf',
                description="GTF annotation file"
            )
            logger.info(f"  ✓ GTF file: {self.gtf}")
        except FileNotFoundError as e:
            raise ValidationError(
                f"GTF file not found: {gtf}\n"
                f"  Please provide a valid GTF annotation file.\n"
                f"  Example: --gtf gencode.v45.annotation.gtf\n"
                f"  \n"
                f"  Common locations:\n"
                f"    - /ref/gencode/\n"
                f"    - /data/annotations/\n"
                f"  \n"
                f"  Download from:\n"
                f"    - GENCODE: https://www.gencodegenes.org/\n"
                f"    - Ensembl: http://www.ensembl.org/"
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
        
        # Validate two-pass mode
        valid_two_pass = ['None', 'Basic']
        if two_pass_mode not in valid_two_pass:
            raise ValidationError(
                f"Invalid two_pass_mode: '{two_pass_mode}'\n"
                f"  Must be one of: {', '.join(valid_two_pass)}\n"
                f"  Recommended: 'Basic' for better sensitivity\n"
                f"  Use 'None' for faster alignment"
            )
        
        # Validate out_sam_type
        valid_sam_types = ['BAM SortedByCoordinate', 'BAM Unsorted', 'SAM']
        if out_sam_type not in valid_sam_types:
            raise ValidationError(
                f"Invalid out_sam_type: '{out_sam_type}'\n"
                f"  Must be one of: {', '.join(valid_sam_types)}\n"
                f"  Recommended: 'BAM SortedByCoordinate'"
            )
        
        # Validate quant_mode
        valid_quant_modes = ['GeneCounts', 'TranscriptomeSAM', '-']
        if quant_mode not in valid_quant_modes:
            raise ValidationError(
                f"Invalid quant_mode: '{quant_mode}'\n"
                f"  Must be one of: {', '.join(valid_quant_modes)}\n"
                f"  'GeneCounts': Generate gene counts\n"
                f"  '-': No quantification"
            )
        
        # Validate STAR numeric parameters
        try:
            validate_positive_integer(out_filter_multimap_nmax, 'out_filter_multimap_nmax')
            if out_filter_mismatch_nmax < 0:
                raise ValidationError("out_filter_mismatch_nmax must be >= 0")
            validate_positive_integer(align_intron_min, 'align_intron_min')
            validate_positive_integer(align_intron_max, 'align_intron_max')
            validate_positive_integer(align_sjdb_overhang_min, 'align_sjdb_overhang_min')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid STAR parameter: {e}\n"
                f"  Please check numeric values are positive"
            ) from e
        
        # Validate intron length relationship
        if align_intron_min >= align_intron_max:
            raise ValidationError(
                f"align_intron_min must be < align_intron_max\n"
                f"  Current: min={align_intron_min}, max={align_intron_max}\n"
                f"  Typical values: min=20, max=1000000"
            )
        
        # Validate mismatch rate
        try:
            validate_probability(
                out_filter_mismatch_nover_read_lmax,
                'out_filter_mismatch_nover_read_lmax'
            )
        except ValidationError as e:
            raise ValidationError(
                f"Invalid out_filter_mismatch_nover_read_lmax: {e}\n"
                f"  Must be between 0 and 1 (fraction of read length)\n"
                f"  Typical: 0.04 (4% of read length)"
            ) from e
        
        # Validate out_sam_unmapped
        valid_unmapped = ['Within', 'None']
        if out_sam_unmapped not in valid_unmapped:
            raise ValidationError(
                f"Invalid out_sam_unmapped: '{out_sam_unmapped}'\n"
                f"  Must be one of: {', '.join(valid_unmapped)}"
            )
        
        # Validate BAM compression
        if out_bam_compression < -1 or out_bam_compression > 10:
            raise ValidationError(
                f"out_bam_compression must be between -1 and 10, got {out_bam_compression}\n"
                f"  -1: uncompressed\n"
                f"  0: uncompressed\n"
                f"  1-9: compression level\n"
                f"  10: maximum compression"
            )
        
        # Validate chimeric parameters
        try:
            validate_positive_integer(chimeric_segment_min, 'chimeric_segment_min')
            validate_positive_integer(chimeric_junction_overhang_min, 'chimeric_junction_overhang_min')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid chimeric parameter: {e}\n"
                f"  Both values must be positive integers"
            ) from e
        
        # Validate featureCounts strand_specific
        if strand_specific not in [0, 1, 2]:
            raise ValidationError(
                f"strand_specific must be 0, 1, or 2, got {strand_specific}\n"
                f"  0 = unstranded\n"
                f"  1 = stranded\n"
                f"  2 = reverse stranded (dUTP)"
            )
        
        # Validate min_overlap
        try:
            validate_positive_integer(min_overlap, 'min_overlap')
        except ValidationError as e:
            raise ValidationError(
                f"Invalid min_overlap: {e}\n"
                f"  Must be a positive integer (typically 1)"
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
        
        # STAR parameters
        self.two_pass_mode = two_pass_mode
        self.out_sam_type = out_sam_type
        self.quant_mode = quant_mode
        self.out_filter_multimap_nmax = out_filter_multimap_nmax
        self.out_filter_mismatch_nmax = out_filter_mismatch_nmax
        self.out_filter_mismatch_nover_read_lmax = out_filter_mismatch_nover_read_lmax
        self.align_intron_min = align_intron_min
        self.align_intron_max = align_intron_max
        self.align_mates_gap_max = align_mates_gap_max
        self.align_sjdb_overhang_min = align_sjdb_overhang_min
        self.out_sam_unmapped = out_sam_unmapped
        self.out_sam_attributes = out_sam_attributes
        self.out_bam_compression = out_bam_compression
        self.limit_bam_sort_ram = limit_bam_sort_ram
        self.chimeric_segment_min = chimeric_segment_min
        self.chimeric_junction_overhang_min = chimeric_junction_overhang_min
        
        # featureCounts parameters
        self.feature_type = feature_type
        self.attribute_type = attribute_type
        self.count_multimapping = count_multimapping
        self.count_primary = count_primary
        self.count_fraction = count_fraction
        self.min_overlap = min_overlap
        self.strand_specific = strand_specific
        self.count_paired_end = count_paired_end
        self.require_both_ends_mapped = require_both_ends_mapped
        self.count_chimeric = count_chimeric
        
        # Create validated directories
        self.bam_dir = self.output_dir / "bam"
        self.bam_dir.mkdir(exist_ok=True, parents=True)
        self.star_dir = self.output_dir / "star"
        self.star_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info("🦖 STAR + featureCounts pipeline initialized")
    
    def get_dependencies(self) -> List[ToolDependency]:
        """Return STAR and featureCounts dependencies."""
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
                name="featureCounts",
                command="featureCounts",
                min_version="2.0.0",
                docker_image="quay.io/biocontainers/subread:2.0.3",
                conda_package="subread",
                hpc_module="subread",
                install_url="https://subread.sourceforge.net/",
                install_conda="conda install -c bioconda subread"
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
        Run STAR + featureCounts on a single sample with error handling.
        
        🆕 v2.2.0: Enhanced error handling with detailed messages.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to STAR index directory
        
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
                f"STAR index not found: {index_path}\n"
                f"  Build index with: STAR --runMode genomeGenerate ..."
            )
        except FileNotFoundError as e:
            raise PipelineError(str(e)) from e
        
        sample_star_dir = self.star_dir / sample.sample_id
        sample_star_dir.mkdir(exist_ok=True, parents=True)
        
        star_log = self.logs_dir / f"{sample.sample_id}_star.log"
        fc_log = self.logs_dir / f"{sample.sample_id}_featurecounts.log"
        
        # ============== STEP 1: STAR Alignment ==============
        logger.info(f"  Running STAR alignment for {sample.sample_id}")
        
        cmd_star = [
            'STAR',
            '--runThreadN', str(self.threads),
            '--genomeDir', str(index_path),
            '--outFileNamePrefix', str(sample_star_dir) + '/',
            '--twopassMode', self.two_pass_mode,
            '--outSAMtype', *self.out_sam_type.split(),
            '--quantMode', self.quant_mode,
            '--outFilterMultimapNmax', str(self.out_filter_multimap_nmax),
            '--outFilterMismatchNmax', str(self.out_filter_mismatch_nmax),
            '--outFilterMismatchNoverReadLmax', str(self.out_filter_mismatch_nover_read_lmax),
            '--alignIntronMin', str(self.align_intron_min),
            '--alignIntronMax', str(self.align_intron_max),
            '--alignMatesGapMax', str(self.align_mates_gap_max),
            '--alignSJDBoverhangMin', str(self.align_sjdb_overhang_min),
            '--outSAMunmapped', self.out_sam_unmapped,
            '--outSAMattributes', self.out_sam_attributes,
            '--outBAMcompression', str(self.out_bam_compression)
        ]
        
        # Add input files
        if sample.is_paired:
            cmd_star.extend(['--readFilesIn', str(sample.fastq_1), str(sample.fastq_2)])
            logger.info(f"    Paired-end mode")
        else:
            cmd_star.extend(['--readFilesIn', str(sample.fastq_1)])
            logger.info(f"    Single-end mode")
        
        # Handle gzipped files
        if str(sample.fastq_1).endswith('.gz'):
            cmd_star.extend(['--readFilesCommand', 'zcat'])
        
        # Memory limit for BAM sorting
        if self.limit_bam_sort_ram > 0:
            cmd_star.extend(['--limitBAMsortRAM', str(self.limit_bam_sort_ram)])
        
        # Run STAR
        returncode, stdout, stderr = self.run_command(cmd_star, star_log)
        
        if returncode != 0:
            raise PipelineError(
                f"STAR alignment failed for {sample.sample_id}\n"
                f"  Check log: {star_log}\n"
                f"  Common issues:\n"
                f"    - Index version mismatch (rebuild with same STAR version)\n"
                f"    - Corrupted FASTQ files\n"
                f"    - Insufficient memory ({self.memory_gb}GB allocated, need 32GB+)\n"
                f"    - Wrong index for organism\n"
                f"    - Disk space full"
            )
        
        # Find BAM file
        bam_file = sample_star_dir / "Aligned.sortedByCoord.out.bam"
        if not bam_file.exists():
            bam_file = sample_star_dir / "Aligned.out.bam"
        
        if not bam_file.exists():
            raise PipelineError(
                f"STAR BAM file not found for {sample.sample_id}\n"
                f"  Expected in: {sample_star_dir}\n"
                f"  STAR completed but produced no BAM output.\n"
                f"  Check log: {star_log}"
            )
        
        # Move BAM to bam directory
        final_bam = self.bam_dir / f"{sample.sample_id}.bam"
        bam_file.rename(final_bam)
        
        # Index BAM
        self._index_bam(final_bam)
        
        # Parse STAR stats
        star_stats = self._parse_star_log(sample_star_dir / "Log.final.out")
        
        # Warn if low mapping rate
        if star_stats.get('uniquely_mapped_percent'):
            if star_stats['uniquely_mapped_percent'] < 50:
                logger.warning(
                    f"Low unique mapping rate for {sample.sample_id}: "
                    f"{star_stats['uniquely_mapped_percent']:.1f}%"
                )
        
        # ============== STEP 2: featureCounts ==============
        logger.info(f"  Running featureCounts for {sample.sample_id}")
        
        counts_file = sample_star_dir / f"{sample.sample_id}_counts.txt"
        
        cmd_fc = [
            'featureCounts',
            '-a', str(self.gtf),
            '-o', str(counts_file),
            '-t', self.feature_type,
            '-g', self.attribute_type,
            '-T', str(self.threads),
            '--minOverlap', str(self.min_overlap),
            '-s', str(self.strand_specific)
        ]
        
        # Paired-end options
        if sample.is_paired:
            cmd_fc.append('-p')  # Count fragments
            if self.count_paired_end:
                cmd_fc.append('--countReadPairs')
            if self.require_both_ends_mapped:
                cmd_fc.append('-B')
            if self.count_chimeric:
                cmd_fc.append('-C')
        
        # Multimapping options
        if self.count_multimapping:
            cmd_fc.append('-M')
        if self.count_primary:
            cmd_fc.append('--primary')
        if self.count_fraction:
            cmd_fc.append('--fraction')
        
        # Add BAM file
        cmd_fc.append(str(final_bam))
        
        # Run featureCounts
        returncode, stdout, stderr = self.run_command(cmd_fc, fc_log)
        
        if returncode != 0:
            raise PipelineError(
                f"featureCounts failed for {sample.sample_id}\n"
                f"  Check log: {fc_log}\n"
                f"  Common issues:\n"
                f"    - GTF and BAM chromosome names don't match\n"
                f"    - Wrong feature_type (use 'exon' for genes)\n"
                f"    - Wrong attribute_type (use 'gene_id')\n"
                f"    - Corrupted BAM file"
            )
        
        # Check output exists
        if not counts_file.exists():
            raise PipelineError(
                f"featureCounts output not found for {sample.sample_id}\n"
                f"  Expected: {counts_file}\n"
                f"  Check log: {fc_log}"
            )
        
        # Parse featureCounts summary
        fc_stats = self._parse_featurecounts_summary(counts_file + ".summary")
        
        # Warn if low assignment rate
        if fc_stats.get('assignment_percent'):
            if fc_stats['assignment_percent'] < 50:
                logger.warning(
                    f"Low assignment rate for {sample.sample_id}: "
                    f"{fc_stats['assignment_percent']:.1f}%"
                )
        
        return {
            'success': True,
            'sample_id': sample.sample_id,
            'bam_file': str(final_bam),
            'counts_file': str(counts_file),
            'log_files': {
                'star': str(star_log),
                'featurecounts': str(fc_log)
            },
            'star_stats': star_stats,
            'featurecounts_stats': fc_stats
        }
    
    def _index_bam(self, bam_file: Path):
        """Index BAM file with samtools."""
        cmd = ['samtools', 'index', str(bam_file)]
        returncode, stdout, stderr = self.run_command(cmd)
        if returncode != 0:
            logger.warning(f"Failed to index BAM: {bam_file}")
    
    def _parse_star_log(self, log_file: Path) -> Dict:
        """Parse STAR Log.final.out file."""
        stats = {}
        if not log_file.exists():
            return stats
        
        try:
            with open(log_file) as f:
                for line in f:
                    if 'Uniquely mapped reads %' in line:
                        stats['uniquely_mapped_percent'] = float(line.split('|')[1].strip().rstrip('%'))
                    elif 'Number of input reads' in line:
                        stats['n_input_reads'] = int(line.split('|')[1].strip())
                    elif 'Uniquely mapped reads number' in line:
                        stats['n_uniquely_mapped'] = int(line.split('|')[1].strip())
        except Exception as e:
            logger.warning(f"Could not parse STAR log: {e}")
        
        return stats
    
    def _parse_featurecounts_summary(self, summary_file: Path) -> Dict:
        """Parse featureCounts summary file."""
        stats = {}
        if not summary_file.exists():
            return stats
        
        try:
            df = pd.read_csv(summary_file, sep='\t', index_col=0)
            col = df.columns[0]
            stats['assigned'] = int(df.loc['Assigned', col])
            stats['total'] = int(df[col].sum())
            if stats['total'] > 0:
                stats['assignment_percent'] = (stats['assigned'] / stats['total']) * 100
        except Exception as e:
            logger.warning(f"Could not parse featureCounts summary: {e}")
        
        return stats
    
    @handle_errors(exit_on_error=False)
    def combine_counts(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """
        Combine featureCounts files into count matrix.
        
        🆕 v2.2.0: Enhanced error handling for combining counts.
        """
        if not sample_results:
            raise PipelineError("No successful samples to combine")
        
        counts_data = {}
        
        for result in sample_results:
            if not result.get('success'):
                logger.warning(f"Skipping failed sample: {result.get('sample_id')}")
                continue
                
            sample_id = result['sample_id']
            counts_file = Path(result['counts_file'])
            
            if not counts_file.exists():
                logger.warning(f"Counts file not found for {sample_id}: {counts_file}")
                continue
            
            try:
                # Read featureCounts output (skip first line, use Geneid as index)
                df = pd.read_csv(counts_file, sep='\t', comment='#', index_col=0)
                # Get last column (counts)
                counts = df.iloc[:, -1]
                counts_data[sample_id] = counts
            except Exception as e:
                logger.warning(f"Could not read counts file for {sample_id}: {e}")
                continue
        
        if not counts_data:
            raise PipelineError("No valid count data found in any samples")
        
        # Combine into matrix
        counts_df = pd.DataFrame(counts_data)
        
        # Round to integers (should already be int)
        counts_df = counts_df.round().astype(int)
        
        # Save counts
        try:
            counts_df.to_csv(output_file)
            logger.info(f"✅ Saved counts to {output_file}")
        except Exception as e:
            raise PipelineError(f"Could not save count matrix to {output_file}: {e}") from e
        
        return counts_df


# =============================================================================
# CLI PARAMETER DEFINITIONS (v2.2.0 Compatible)
# =============================================================================

STAR_FC_CLI_PARAMS = {
    'gtf': {
        'flag': '-g',
        'type': str,
        'required': True,
        'help': 'GTF annotation file (REQUIRED)',
        'param_name': 'gtf'
    },
    'two-pass-mode': {
        'flag': '--twopassMode',
        'type': str,
        'default': 'Basic',
        'choices': ['None', 'Basic'],
        'help': 'STAR two-pass mode (Basic=better sensitivity)',
        'param_name': 'two_pass_mode'
    },
    'out-sam-type': {
        'flag': '--outSAMtype',
        'type': str,
        'default': 'BAM SortedByCoordinate',
        'help': 'Output type (BAM SortedByCoordinate recommended)',
        'param_name': 'out_sam_type'
    },
    'strand-specific': {
        'flag': '-s',
        'type': int,
        'default': 0,
        'choices': [0, 1, 2],
        'help': 'Strandedness: 0=unstranded, 1=stranded, 2=reverse',
        'param_name': 'strand_specific'
    },
    'feature-type': {
        'flag': '-t',
        'type': str,
        'default': 'exon',
        'help': 'featureCounts feature type (typically exon)',
        'param_name': 'feature_type'
    },
    'attribute-type': {
        'flag': '-g',
        'type': str,
        'default': 'gene_id',
        'help': 'featureCounts attribute type (typically gene_id)',
        'param_name': 'attribute_type'
    },
}


if __name__ == '__main__':
    print(f"🦖 RAPTOR STAR + featureCounts Pipeline v{StarFeatureCountsPipeline.PIPELINE_VERSION}")
    print("=" * 70)
    print("\n✅ v2.2.0 Features:")
    print("   • Full input validation with clear error messages")
    print("   • Integrated error handling (@handle_errors)")
    print("   • Compatible with CLI system")
    print("\n⚡ Standard alignment-based workflow")
    print("   Produces BAM files + gene counts")
    print("\n🔧 Available CLI parameters:")
    for param, info in STAR_FC_CLI_PARAMS.items():
        required = " (REQUIRED)" if info.get('required') else ""
        print(f"  --{param}: {info['help']}{required}")
    print("\n⚠️  High memory requirement: ~32GB")
