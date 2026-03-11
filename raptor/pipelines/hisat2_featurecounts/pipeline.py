#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - HISAT2 + featureCounts Pipeline
FULLY INTEGRATED WITH VALIDATION & ERROR HANDLING

Production pipeline using HISAT2 alignment + featureCounts.
Lower memory alternative to STAR.

Features:
- Fast spliced alignment with HISAT2
- Lower memory footprint than STAR
- Gene-level quantification with featureCounts
- BAM files for downstream analysis
- ✅ Full input validation with clear error messages
- ✅ Integrated error handling
- ✅ Compatible with CLI system

Requirements:
- HISAT2 >= 2.2.0
- featureCounts (subread >= 2.0.0)
- samtools >= 1.10
- HISAT2 index + GTF annotation

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


class Hisat2FeatureCountsPipeline(BasePipeline):
    """
    HISAT2 alignment + featureCounts quantification pipeline.
    
    Lower memory alternative to STAR. Good for systems with limited RAM
    or when processing many samples.
    
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
        GTF annotation file (REQUIRED)
    threads : int
        Number of threads (default: 8, must be > 0)
    memory_gb : int
        Memory in GB (default: 16, must be > 0)
    min_intron_len : int
        Minimum intron length (default: 20, must be > 0)
    max_intron_len : int
        Maximum intron length (default: 500000, must be > min_intron_len)
    rna_strandness : str
        Strand-specific protocol: '', 'RF', 'FR', 'F', 'R' (default: '')
    dta : bool
        Report alignments for transcript assembly (default: False)
    dta_cufflinks : bool
        Report alignments for Cufflinks (default: False)
    feature_type : str
        Feature type in GTF (default: 'exon')
    attribute_type : str
        GTF attribute for grouping (default: 'gene_id')
    count_multimapping : bool
        Count multi-mapping reads (default: False)
    strand_specific : int
        Strand specificity: 0=unstranded, 1=stranded, 2=reverse (default: 0)
    
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
    >>> pipeline = Hisat2FeatureCountsPipeline(
    ...     output_dir='results/',
    ...     gtf='annotation.gtf',
    ...     threads=16
    ... )
    >>> result = pipeline.run('samples.csv', 'hisat2_index/genome')
    """
    
    PIPELINE_NAME = "hisat2_featurecounts"
    PIPELINE_VERSION = "2.2.0"
    DESCRIPTION = "HISAT2 alignment + featureCounts (low memory)"
    DOCKER_IMAGE = "quay.io/biocontainers/hisat2:2.2.1"
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        output_dir: Union[str, Path],
        gtf: Union[str, Path],  # REQUIRED
        threads: int = 8,
        memory_gb: int = 16,  # HISAT2 uses less memory
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = True,
        extra_args: Optional[Dict[str, str]] = None,
        # HISAT2 parameters
        min_intron_len: int = 20,
        max_intron_len: int = 500000,
        rna_strandness: str = '',
        dta: bool = False,
        dta_cufflinks: bool = False,
        pen_cansplice: int = 0,
        pen_noncansplice: int = 12,
        pen_canintronlen: str = 'G,-8,1',
        pen_noncanintronlen: str = 'G,-8,1',
        no_softclip: bool = False,
        no_spliced_alignment: bool = False,
        novel_splicesite_outfile: Optional[str] = None,
        # featureCounts parameters
        feature_type: str = 'exon',
        attribute_type: str = 'gene_id',
        count_multimapping: bool = False,
        strand_specific: int = 0
    ):
        """
        Initialize HISAT2 + featureCounts pipeline with validation.
        
        🆕 v2.2.0: All inputs are validated with clear error messages.
        """
        # ====================================================================
        # STEP 1: VALIDATE ALL INPUTS (v2.2.0 Enhancement)
        # ====================================================================
        logger.info("🔍 Validating pipeline inputs...")
        
        # Validate GTF file (REQUIRED)
        try:
            self.gtf = validate_file_path(
                gtf,
                must_exist=True,
                extension='.gtf',
                description="GTF annotation file"
            )
            logger.info(f"  ✓ GTF file validated: {self.gtf}")
        except FileNotFoundError as e:
            raise ValidationError(
                f"GTF file not found: {gtf}\n"
                f"  Please provide a valid GTF annotation file.\n"
                f"  Example: --gtf gencode.v45.annotation.gtf"
            ) from e
        except ValidationError as e:
            raise ValidationError(
                f"Invalid GTF file: {e}\n"
                f"  GTF files must have .gtf extension.\n"
                f"  If you have a .gff file, convert it to GTF format first."
            ) from e
        
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
            validate_positive_integer(min_intron_len, 'min_intron_len')
            validate_positive_integer(max_intron_len, 'max_intron_len')
            logger.info(f"  ✓ Numeric parameters validated")
        except ValidationError as e:
            raise ValidationError(
                f"Invalid numeric parameter: {e}\n"
                f"  All numeric values must be positive integers.\n"
                f"  Current values: threads={threads}, memory_gb={memory_gb}"
            ) from e
        
        # Validate intron length relationship
        if min_intron_len >= max_intron_len:
            raise ValidationError(
                f"min_intron_len ({min_intron_len}) must be < max_intron_len ({max_intron_len})\n"
                f"  Typical values: min=20, max=500000"
            )
        
        # Validate strand specificity
        if rna_strandness not in ['', 'RF', 'FR', 'F', 'R']:
            raise ValidationError(
                f"Invalid rna_strandness: '{rna_strandness}'\n"
                f"  Must be one of: '', 'RF', 'FR', 'F', 'R'\n"
                f"  RF/FR for paired-end, F/R for single-end, '' for unstranded"
            )
        
        # Validate featureCounts strand parameter
        if strand_specific not in [0, 1, 2]:
            raise ValidationError(
                f"Invalid strand_specific: {strand_specific}\n"
                f"  Must be 0 (unstranded), 1 (stranded), or 2 (reverse stranded)"
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
        
        # HISAT2 parameters
        self.min_intron_len = min_intron_len
        self.max_intron_len = max_intron_len
        self.rna_strandness = rna_strandness
        self.dta = dta
        self.dta_cufflinks = dta_cufflinks
        self.pen_cansplice = pen_cansplice
        self.pen_noncansplice = pen_noncansplice
        self.pen_canintronlen = pen_canintronlen
        self.pen_noncanintronlen = pen_noncanintronlen
        self.no_softclip = no_softclip
        self.no_spliced_alignment = no_spliced_alignment
        self.novel_splicesite_outfile = novel_splicesite_outfile
        
        # featureCounts parameters
        self.feature_type = feature_type
        self.attribute_type = attribute_type
        self.count_multimapping = count_multimapping
        self.strand_specific = strand_specific
        
        # Create validated directories
        self.bam_dir = self.output_dir / "bam"
        self.bam_dir.mkdir(exist_ok=True, parents=True)
        self.hisat2_dir = self.output_dir / "hisat2"
        self.hisat2_dir.mkdir(exist_ok=True, parents=True)
        
        logger.info("🦖 HISAT2 + featureCounts pipeline initialized")
    
    def get_dependencies(self) -> List[ToolDependency]:
        """Return HISAT2, samtools, and featureCounts dependencies."""
        return [
            ToolDependency(
                name="HISAT2",
                command="hisat2",
                min_version="2.2.0",
                docker_image="quay.io/biocontainers/hisat2:2.2.1",
                conda_package="hisat2",
                hpc_module="hisat2",
                install_url="http://daehwankimlab.github.io/hisat2/",
                install_conda="conda install -c bioconda hisat2"
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
        Run HISAT2 + featureCounts on a single sample with error handling.
        
        🆕 v2.2.0: Enhanced error handling with detailed messages.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to HISAT2 index (base name, not directory)
        
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
                f"HISAT2 index not found: {index_path}\n"
                f"  Build index with: hisat2-build genome.fa genome_index"
            )
        except FileNotFoundError as e:
            raise PipelineError(str(e)) from e
        
        sample_dir = self.hisat2_dir / sample.sample_id
        sample_dir.mkdir(exist_ok=True, parents=True)
        
        hisat2_log = self.logs_dir / f"{sample.sample_id}_hisat2.log"
        fc_log = self.logs_dir / f"{sample.sample_id}_featurecounts.log"
        
        sam_file = sample_dir / f"{sample.sample_id}.sam"
        bam_unsorted = sample_dir / f"{sample.sample_id}.unsorted.bam"
        bam_sorted = self.bam_dir / f"{sample.sample_id}.bam"
        
        # ============== STEP 1: HISAT2 Alignment ==============
        logger.info(f"  Running HISAT2 alignment for {sample.sample_id}")
        
        cmd_hisat2 = [
            'hisat2',
            '-x', str(index_path),
            '-p', str(self.threads),
            '--min-intronlen', str(self.min_intron_len),
            '--max-intronlen', str(self.max_intron_len),
            '-S', str(sam_file)
        ]
        
        # Input files
        if sample.is_paired:
            cmd_hisat2.extend(['-1', str(sample.fastq_1), '-2', str(sample.fastq_2)])
        else:
            cmd_hisat2.extend(['-U', str(sample.fastq_1)])
        
        # Strand-specific
        if self.rna_strandness:
            cmd_hisat2.extend(['--rna-strandness', self.rna_strandness])
        
        # DTA options
        if self.dta:
            cmd_hisat2.append('--dta')
        if self.dta_cufflinks:
            cmd_hisat2.append('--dta-cufflinks')
        
        # Splicing parameters
        if self.pen_cansplice != 0:
            cmd_hisat2.extend(['--pen-cansplice', str(self.pen_cansplice)])
        if self.pen_noncansplice != 12:
            cmd_hisat2.extend(['--pen-noncansplice', str(self.pen_noncansplice)])
        cmd_hisat2.extend(['--pen-canintronlen', self.pen_canintronlen])
        cmd_hisat2.extend(['--pen-noncanintronlen', self.pen_noncanintronlen])
        
        if self.no_softclip:
            cmd_hisat2.append('--no-softclip')
        if self.no_spliced_alignment:
            cmd_hisat2.append('--no-spliced-alignment')
        
        if self.novel_splicesite_outfile:
            cmd_hisat2.extend(['--novel-splicesite-outfile', self.novel_splicesite_outfile])
        
        # Run HISAT2
        returncode, stdout, stderr = self.run_command(cmd_hisat2, hisat2_log)
        
        if returncode != 0 or not sam_file.exists():
            raise PipelineError(
                f"HISAT2 alignment failed for {sample.sample_id}\n"
                f"  Check log: {hisat2_log}\n"
                f"  Common issues:\n"
                f"    - Invalid HISAT2 index\n"
                f"    - Corrupted FASTQ files\n"
                f"    - Insufficient memory ({self.memory_gb}GB allocated)"
            )
        
        # Parse alignment stats
        alignment_stats = self._parse_hisat2_stderr(hisat2_log)
        
        # ============== STEP 2: SAM to BAM ==============
        logger.info(f"  Converting SAM to BAM for {sample.sample_id}")
        
        # Convert to BAM
        cmd_view = ['samtools', 'view', '-@ str(self.threads)', '-bS', str(sam_file), '-o', str(bam_unsorted)]
        returncode, stdout, stderr = self.run_command(cmd_view)
        
        if returncode != 0:
            raise PipelineError(
                f"SAM to BAM conversion failed for {sample.sample_id}\n"
                f"  This usually means samtools is not installed or the SAM file is corrupted"
            )
        
        # Sort BAM
        cmd_sort = ['samtools', 'sort', '-@', str(self.threads), '-o', str(bam_sorted), str(bam_unsorted)]
        returncode, stdout, stderr = self.run_command(cmd_sort)
        
        if returncode != 0 or not bam_sorted.exists():
            raise PipelineError(
                f"BAM sorting failed for {sample.sample_id}\n"
                f"  Check if you have enough disk space and memory"
            )
        
        # Index BAM
        cmd_index = ['samtools', 'index', str(bam_sorted)]
        self.run_command(cmd_index)
        
        # Clean up temporary files
        if sam_file.exists():
            sam_file.unlink()
        if bam_unsorted.exists():
            bam_unsorted.unlink()
        
        # ============== STEP 3: featureCounts ==============
        logger.info(f"  Running featureCounts for {sample.sample_id}")
        
        counts_file = sample_dir / f"{sample.sample_id}_counts.txt"
        
        cmd_fc = [
            'featureCounts',
            '-a', str(self.gtf),
            '-o', str(counts_file),
            '-t', self.feature_type,
            '-g', self.attribute_type,
            '-T', str(self.threads),
            '-s', str(self.strand_specific)
        ]
        
        # Paired-end options
        if sample.is_paired:
            cmd_fc.extend(['-p', '--countReadPairs'])
        
        # Multimapping
        if self.count_multimapping:
            cmd_fc.append('-M')
        
        # Input BAM
        cmd_fc.append(str(bam_sorted))
        
        # Run featureCounts
        returncode, stdout, stderr = self.run_command(cmd_fc, fc_log)
        
        if returncode != 0:
            raise PipelineError(
                f"featureCounts failed for {sample.sample_id}\n"
                f"  Check log: {fc_log}\n"
                f"  Common issues:\n"
                f"    - GTF format doesn't match BAM chromosome names\n"
                f"    - Wrong feature_type (current: {self.feature_type})\n"
                f"    - Wrong attribute_type (current: {self.attribute_type})"
            )
        
        # Parse featureCounts stats
        fc_stats = self._parse_fc_summary(Path(str(counts_file) + ".summary"))
        
        return {
            'success': True,
            'sample_id': sample.sample_id,
            'bam_file': str(bam_sorted),
            'counts_file': str(counts_file),
            'hisat2_log': str(hisat2_log),
            'fc_log': str(fc_log),
            'alignment_stats': alignment_stats,
            'fc_stats': fc_stats,
            'overall_alignment_rate': alignment_stats.get('overall_alignment_rate'),
            'assigned_rate': fc_stats.get('assigned_percent')
        }
    
    def _parse_hisat2_stderr(self, log_file: Path) -> Dict:
        """Parse HISAT2 alignment summary from log."""
        stats = {}
        if not log_file.exists():
            return stats
        
        try:
            with open(log_file) as f:
                content = f.read()
                
            # Parse key metrics
            import re
            
            # Total reads
            match = re.search(r'(\d+) reads', content)
            if match:
                stats['total_reads'] = int(match.group(1))
            
            # Overall alignment rate
            match = re.search(r'([\d.]+)% overall alignment rate', content)
            if match:
                stats['overall_alignment_rate'] = float(match.group(1))
            
            # Unique alignments
            match = re.search(r'(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time', content)
            if match:
                stats['unique_alignments'] = int(match.group(1))
                stats['unique_rate'] = float(match.group(2))
                
        except Exception as e:
            logger.warning(f"Could not parse HISAT2 log: {e}")
        
        return stats
    
    def _parse_fc_summary(self, summary_file: Path) -> Dict:
        """Parse featureCounts summary file."""
        stats = {}
        if not summary_file.exists():
            return stats
        
        try:
            df = pd.read_csv(summary_file, sep='\t', index_col=0)
            if len(df.columns) > 0:
                col = df.columns[0]
                total = df[col].sum()
                stats['total'] = total
                stats['assigned'] = df.loc['Assigned', col] if 'Assigned' in df.index else 0
                stats['assigned_percent'] = stats['assigned'] / total * 100 if total > 0 else 0
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
        Combine featureCounts output files into count matrix.
        
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
                # Read featureCounts output
                df = pd.read_csv(counts_file, sep='\t', comment='#')
                
                # Gene ID and counts columns
                gene_col = 'Geneid'
                count_col = df.columns[-1]
                
                counts = df.set_index(gene_col)[count_col]
                counts_data[sample_id] = counts
            except Exception as e:
                logger.warning(f"Could not read counts for {sample_id}: {e}")
                continue
        
        if not counts_data:
            raise PipelineError("No valid count data found in any samples")
        
        # Combine into matrix
        counts_df = pd.DataFrame(counts_data)
        counts_df = counts_df.fillna(0).astype(int)
        
        # Save
        try:
            counts_df.to_csv(output_file)
            logger.info(f"✅ Saved gene counts to {output_file}")
        except Exception as e:
            raise PipelineError(f"Could not save count matrix to {output_file}: {e}") from e
        
        return counts_df


# =============================================================================
# CLI PARAMETER DEFINITIONS (v2.2.0 Compatible)
# =============================================================================

HISAT2_FC_CLI_PARAMS = {
    # Required
    'gtf': {
        'flag': '--gtf',
        'type': str,
        'required': True,
        'help': 'GTF annotation file (REQUIRED)',
        'param_name': 'gtf'
    },
    
    # HISAT2 parameters
    'min-intron-len': {
        'flag': '--min-intron-len',
        'type': int,
        'default': 20,
        'help': 'Minimum intron length',
        'param_name': 'min_intron_len'
    },
    'max-intron-len': {
        'flag': '--max-intron-len',
        'type': int,
        'default': 500000,
        'help': 'Maximum intron length',
        'param_name': 'max_intron_len'
    },
    'rna-strandness': {
        'flag': '--rna-strandness',
        'type': str,
        'default': '',
        'choices': ['', 'RF', 'FR', 'F', 'R'],
        'help': 'Strand-specific protocol (RF/FR for paired, F/R for single)',
        'param_name': 'rna_strandness'
    },
    'dta': {
        'flag': '--dta',
        'type': bool,
        'default': False,
        'help': 'Report alignments for transcript assembly',
        'param_name': 'dta'
    },
    
    # featureCounts parameters
    'feature-type': {
        'flag': '--feature-type',
        'type': str,
        'default': 'exon',
        'help': 'Feature type in GTF',
        'param_name': 'feature_type'
    },
    'attribute-type': {
        'flag': '--attribute-type',
        'type': str,
        'default': 'gene_id',
        'help': 'GTF attribute for grouping',
        'param_name': 'attribute_type'
    },
    'strand-specific': {
        'flag': '--strand-specific',
        'type': int,
        'default': 0,
        'choices': [0, 1, 2],
        'help': 'Strand specificity: 0=unstranded, 1=stranded, 2=reverse',
        'param_name': 'strand_specific'
    }
}


if __name__ == '__main__':
    print(f"🦖 RAPTOR HISAT2+featureCounts Pipeline v{Hisat2FeatureCountsPipeline.PIPELINE_VERSION}")
    print("=" * 70)
    print("\n✅ v2.2.0 Features:")
    print("   • Full input validation with clear error messages")
    print("   • Integrated error handling (@handle_errors)")
    print("   • Compatible with CLI system")
    print("\n💡 Lower memory alternative to STAR")
    print("   Recommended when RAM is limited (<32 GB)")
    print("\n📋 Required parameters:")
    print("  --gtf PATH : GTF annotation file")
    print("\n🔧 Available CLI parameters:")
    for param, info in HISAT2_FC_CLI_PARAMS.items():
        if param != 'gtf':
            print(f"  --{param}: {info['help']} (default: {info.get('default', 'REQUIRED')})")
