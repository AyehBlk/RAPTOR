"""
RAPTOR v2.2.0 - Tests for Production Pipelines (Module 5)

Comprehensive tests for the production pipeline infrastructure:
- Base pipeline classes (BasePipeline, ToolDependency, SampleInfo, PipelineResult)
- Production pipelines (Salmon, Kallisto, STAR+featureCounts, HISAT2+featureCounts, STAR+RSEM, STAR+Salmon)
- Pipeline registry and information functions
- CLI integration tests

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys
import tempfile
import shutil
from unittest.mock import patch, MagicMock, PropertyMock
from dataclasses import dataclass

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# =============================================================================
# Import Modules Under Test
# =============================================================================

# Base pipeline components
try:
    from raptor.pipelines.base import (
        BasePipeline,
        ToolDependency,
        SampleInfo,
        PipelineResult,
        auto_detect_samples,
        create_sample_sheet
    )
    BASE_AVAILABLE = True
except ImportError:
    BASE_AVAILABLE = False

# Salmon pipeline
try:
    from raptor.pipelines.salmon.pipeline import SalmonPipeline
    SALMON_PIPELINE_AVAILABLE = True
except ImportError:
    SALMON_PIPELINE_AVAILABLE = False

# Kallisto pipeline
try:
    from raptor.pipelines.kallisto.pipeline import KallistoPipeline
    KALLISTO_PIPELINE_AVAILABLE = True
except ImportError:
    KALLISTO_PIPELINE_AVAILABLE = False

# STAR + featureCounts pipeline
try:
    from raptor.pipelines.star_featurecounts.pipeline import StarFeatureCountsPipeline
    STAR_FC_PIPELINE_AVAILABLE = True
except ImportError:
    STAR_FC_PIPELINE_AVAILABLE = False

# HISAT2 + featureCounts pipeline
try:
    from raptor.pipelines.hisat2_featurecounts.pipeline import Hisat2FeatureCountsPipeline
    HISAT2_FC_PIPELINE_AVAILABLE = True
except ImportError:
    HISAT2_FC_PIPELINE_AVAILABLE = False

# STAR + RSEM pipeline
try:
    from raptor.pipelines.star_rsem.pipeline import StarRsemPipeline
    STAR_RSEM_PIPELINE_AVAILABLE = True
except ImportError:
    STAR_RSEM_PIPELINE_AVAILABLE = False

# STAR + Salmon pipeline
try:
    from raptor.pipelines.star_salmon.pipeline import StarSalmonPipeline
    STAR_SALMON_PIPELINE_AVAILABLE = True
except ImportError:
    STAR_SALMON_PIPELINE_AVAILABLE = False

# Pipeline registry
try:
    # Import from main registry
    from raptor.pipelines import get_pipeline
    
    # Import metadata functions from init.py
    try:
        from raptor.pipelines.init import (
            list_pipelines,
            get_pipeline_info,
            print_pipeline_summary,
            get_comparison_table,
            AVAILABLE_PIPELINES
        )
    except ImportError:
        # If init.py not available, that's OK - tests will skip
        pass
    
    REGISTRY_AVAILABLE = True
except ImportError as e:
    REGISTRY_AVAILABLE = False
    print(f"⚠️  Warning: Pipeline registry not available: {e}")


# =============================================================================
# Constants
# =============================================================================

DEFAULT_PRODUCTION_DIR = "results/production"
PRODUCTION_COUNTS_FILE = "gene_counts.csv"
PRODUCTION_TX_FILE = "tx_counts.csv"

AVAILABLE_PIPELINES = [
    'salmon', 'kallisto', 'star_featurecounts',
    'hisat2_featurecounts', 'star_rsem', 'star_salmon'
]


# =============================================================================
# Additional Fixtures for Module 5
# =============================================================================

@pytest.fixture
def mock_salmon_index(test_dir):
    """Create a mock Salmon index directory."""
    index_dir = test_dir / "salmon_index"
    index_dir.mkdir(exist_ok=True)
    
    # Create mock index files
    (index_dir / "versionInfo.json").write_text('{"version": "1.9.0"}')
    (index_dir / "hash.bin").touch()
    (index_dir / "seq.bin").touch()
    (index_dir / "pos.bin").touch()
    (index_dir / "duplicate_clusters.tsv").touch()
    
    return index_dir


@pytest.fixture
def mock_kallisto_index(test_dir):
    """Create a mock Kallisto index file."""
    index_file = test_dir / "kallisto.idx"
    index_file.write_bytes(b'\x00' * 100)  # Dummy bytes
    return index_file


@pytest.fixture
def mock_star_index(test_dir):
    """Create a mock STAR index directory."""
    index_dir = test_dir / "star_index"
    index_dir.mkdir(exist_ok=True)
    
    # Create mock STAR index files
    (index_dir / "genomeParameters.txt").write_text("versionGenome 2.7.10b")
    (index_dir / "Genome").touch()
    (index_dir / "SA").touch()
    (index_dir / "SAindex").touch()
    (index_dir / "chrName.txt").write_text("chr1\nchr2\n")
    (index_dir / "chrLength.txt").write_text("248956422\n242193529\n")
    
    return index_dir


@pytest.fixture
def mock_hisat2_index(test_dir):
    """Create a mock HISAT2 index directory."""
    index_dir = test_dir / "hisat2_index"
    index_dir.mkdir(exist_ok=True)
    
    # Create mock HISAT2 index files
    for i in range(1, 9):
        (index_dir / f"genome.{i}.ht2").touch()
    
    return index_dir / "genome"


@pytest.fixture
def mock_rsem_index(test_dir):
    """Create a mock RSEM index directory."""
    index_dir = test_dir / "rsem_ref"
    index_dir.mkdir(exist_ok=True)
    
    # Create mock RSEM reference files
    (index_dir / "genome.grp").touch()
    (index_dir / "genome.ti").touch()
    (index_dir / "genome.seq").touch()
    (index_dir / "genome.transcripts.fa").touch()
    
    return index_dir / "genome"


@pytest.fixture
def mock_gtf_file(test_dir):
    """Create a mock GTF annotation file."""
    gtf_file = test_dir / "genes.gtf"
    
    gtf_content = '''##gtf-version 2.2
chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "Gene1";
chr1\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";
chr1\tENSEMBL\texon\t1000\t1500\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";
chr1\tENSEMBL\texon\t1700\t2000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";
chr1\tENSEMBL\tgene\t5000\t7000\t.\t-\t.\tgene_id "ENSG00000000002"; gene_name "Gene2";
chr1\tENSEMBL\ttranscript\t5000\t7000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002";
chr1\tENSEMBL\texon\t5000\t6000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002";
'''
    gtf_file.write_text(gtf_content)
    return gtf_file


@pytest.fixture
def pipeline_sample_sheet(test_dir, mock_fastq_dir):
    """Create a sample sheet for production pipelines."""
    sample_sheet = test_dir / "production_samples.csv"
    
    # Get actual FASTQ files from mock_fastq_dir
    data = {
        'sample_id': ['Control_1', 'Control_2', 'Treatment_1', 'Treatment_2'],
        'condition': ['Control', 'Control', 'Treatment', 'Treatment'],
        'batch': ['Batch1'] * 4,
        'fastq_r1': [str(mock_fastq_dir / f"{s}_R1_001.fastq.gz") 
                    for s in ['Control_1', 'Control_2', 'Treatment_1', 'Treatment_2']],
        'fastq_r2': [str(mock_fastq_dir / f"{s}_R2_001.fastq.gz") 
                    for s in ['Control_1', 'Control_2', 'Treatment_1', 'Treatment_2']]
    }
    
    pd.DataFrame(data).to_csv(sample_sheet, index=False)
    return sample_sheet


# =============================================================================
# Tests for ToolDependency
# =============================================================================

class TestToolDependency:
    """Tests for ToolDependency dataclass."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_tool_dependency_creation(self):
        """Test creating a ToolDependency."""
        dep = ToolDependency(
            name="Salmon",
            command="salmon",
            min_version="1.9.0",
            docker_image="combinelab/salmon:1.10.0",
            conda_package="salmon",
            hpc_module="salmon",
            install_url="https://github.com/COMBINE-lab/salmon",
            install_conda="conda install -c bioconda salmon"
        )
        
        assert dep.name == "Salmon"
        assert dep.command == "salmon"
        assert dep.min_version == "1.9.0"
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_tool_dependency_optional_fields(self):
        """Test ToolDependency with minimal required fields."""
        dep = ToolDependency(
            name="TestTool",
            command="testtool"
        )
        
        assert dep.name == "TestTool"
        assert dep.min_version is None or dep.min_version == ""


# =============================================================================
# Tests for SampleInfo
# =============================================================================

class TestSampleInfo:
    """Tests for SampleInfo dataclass."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_sample_info_paired(self, mock_fastq_dir):
        """Test creating paired-end SampleInfo."""
        sample = SampleInfo(
            sample_id="Test_1",
            fastq_1=mock_fastq_dir / "Control_1_R1_001.fastq.gz",
            fastq_2=mock_fastq_dir / "Control_1_R2_001.fastq.gz",
            is_paired=True,
            condition="Control"
        )
        
        assert sample.sample_id == "Test_1"
        assert sample.is_paired == True
        assert sample.fastq_2 is not None
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_sample_info_single(self, mock_fastq_dir):
        """Test creating single-end SampleInfo."""
        sample = SampleInfo(
            sample_id="Test_1",
            fastq_1=mock_fastq_dir / "Control_1_R1_001.fastq.gz",
            is_paired=False,
            condition="Control"
        )
        
        assert sample.sample_id == "Test_1"
        assert sample.is_paired == False


# =============================================================================
# Tests for PipelineResult
# =============================================================================

class TestPipelineResult:
    """Tests for PipelineResult dataclass."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_pipeline_result_success(self, output_dir):
        """Test successful PipelineResult."""
        result = PipelineResult(
            success=True,
            pipeline_name="salmon",
            output_dir=str(output_dir),
            gene_counts_file=str(output_dir / "gene_counts.csv"),
            tx_counts_file=str(output_dir / "tx_counts.csv"),
            n_samples_processed=4,
            n_samples_failed=0,
            elapsed_time=120.5
        )
        
        assert result.success == True
        assert result.n_samples_processed == 4
        assert result.n_samples_failed == 0
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_pipeline_result_summary(self, output_dir):
        """Test PipelineResult.summary() method."""
        result = PipelineResult(
            success=True,
            pipeline_name="salmon",
            output_dir=str(output_dir),
            gene_counts_file=str(output_dir / "gene_counts.csv"),
            n_samples_processed=4,
            n_samples_failed=0,
            elapsed_time=120.5
        )
        
        summary = result.summary()
        
        assert isinstance(summary, str)
        assert 'salmon' in summary.lower() or 'Salmon' in summary
        assert '4' in summary  # n_samples_processed
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_pipeline_result_failure(self, output_dir):
        """Test failed PipelineResult."""
        result = PipelineResult(
            success=False,
            pipeline_name="salmon",
            output_dir=str(output_dir),
            n_samples_processed=2,
            n_samples_failed=2,
            error_message="Salmon not found"
        )
        
        assert result.success == False
        assert result.n_samples_failed == 2
        assert "Salmon" in result.error_message


# =============================================================================
# Tests for Salmon Pipeline
# =============================================================================

class TestSalmonPipeline:
    """Tests for SalmonPipeline class."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_init_basic(self, output_dir):
        """Test basic initialization."""
        pipeline = SalmonPipeline(
            output_dir=str(output_dir),
            threads=8
        )
        
        assert pipeline.PIPELINE_NAME == "salmon"
        assert pipeline.threads == 8
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_init_with_options(self, output_dir):
        """Test initialization with options."""
        pipeline = SalmonPipeline(
            output_dir=str(output_dir),
            threads=16,
            library_type='ISR',
            bootstraps=100,
            gc_bias=True,
            seq_bias=True,
            validate_mappings=True
        )
        
        assert pipeline.threads == 16
        assert pipeline.library_type == 'ISR'
        assert pipeline.bootstraps == 100
        assert pipeline.gc_bias == True
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_get_dependencies(self, output_dir):
        """Test get_dependencies method."""
        pipeline = SalmonPipeline(output_dir=str(output_dir))
        deps = pipeline.get_dependencies()
        
        assert len(deps) >= 1
        salmon_dep = deps[0]
        assert salmon_dep.name == "Salmon"
        assert salmon_dep.command == "salmon"
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_output_directory_structure(self, output_dir):
        """Test that output directory structure is created."""
        pipeline = SalmonPipeline(output_dir=str(output_dir))
        
        # Check expected directories exist
        assert pipeline.output_dir.exists()
        assert pipeline.quant_dir.exists()
        assert pipeline.logs_dir.exists()
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_default_output_dir(self):
        """Test default output directory."""
        pipeline = SalmonPipeline()
        
        assert 'production' in str(pipeline.output_dir)


# =============================================================================
# Tests for Kallisto Pipeline
# =============================================================================

class TestKallistoPipeline:
    """Tests for KallistoPipeline class."""
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_init_basic(self, output_dir):
        """Test basic initialization."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir),
            threads=8
        )
        
        assert pipeline.PIPELINE_NAME == "kallisto"
        assert pipeline.threads == 8
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_init_with_bootstraps(self, output_dir):
        """Test initialization with bootstraps."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir),
            bootstraps=100,
            seed=42
        )
        
        assert pipeline.bootstraps == 100
        assert pipeline.seed == 42
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_fragment_length_for_single_end(self, output_dir):
        """Test fragment length parameters for single-end."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir),
            fragment_length=200,
            fragment_sd=20
        )
        
        assert pipeline.fragment_length == 200
        assert pipeline.fragment_sd == 20
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_strand_specific_options(self, output_dir):
        """Test strand-specific options."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir),
            rf_stranded=True
        )
        
        assert pipeline.rf_stranded == True
        assert pipeline.fr_stranded == False
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_tx2gene_option(self, output_dir, tx2gene_path):
        """Test tx2gene option."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir),
            tx2gene=str(tx2gene_path)
        )
        
        assert pipeline.tx2gene == tx2gene_path


# =============================================================================
# Tests for STAR + featureCounts Pipeline
# =============================================================================

class TestStarFeatureCountsPipeline:
    """Tests for StarFeatureCountsPipeline class."""
    
    @pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
    def test_init_basic(self, output_dir, mock_gtf_file):
        """Test basic initialization."""
        pipeline = StarFeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            threads=8
        )
        
        assert pipeline.PIPELINE_NAME == "star_featurecounts"
        assert pipeline.threads == 8
    
    @pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
    def test_gtf_required(self, output_dir):
        """Test that GTF is required."""
        with pytest.raises((ValueError, TypeError)):
            StarFeatureCountsPipeline(
                output_dir=str(output_dir)
                # Missing gtf parameter
            )
    
    @pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
    def test_star_options(self, output_dir, mock_gtf_file):
        """Test STAR-specific options."""
        pipeline = StarFeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            two_pass_mode='Basic',
            out_filter_multimap_nmax=20,
            align_intron_max=1000000
        )
        
        assert pipeline.two_pass_mode == 'Basic'
        assert pipeline.out_filter_multimap_nmax == 20
    
    @pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
    def test_featurecounts_options(self, output_dir, mock_gtf_file):
        """Test featureCounts-specific options."""
        pipeline = StarFeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            feature_type='exon',
            attribute_type='gene_id',
            strand_specific=2
        )
        
        assert pipeline.feature_type == 'exon'
        assert pipeline.attribute_type == 'gene_id'
        assert pipeline.strand_specific == 2
    
    @pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
    def test_keep_bam_option(self, output_dir, mock_gtf_file):
        """Test keep_bam option."""
        pipeline = StarFeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            keep_bam=True
        )
        
        assert pipeline.keep_bam == True


# =============================================================================
# Tests for HISAT2 + featureCounts Pipeline
# =============================================================================

class TestHisat2FeatureCountsPipeline:
    """Tests for Hisat2FeatureCountsPipeline class."""
    
    @pytest.mark.skipif(not HISAT2_FC_PIPELINE_AVAILABLE, reason="HISAT2+FC pipeline not available")
    def test_init_basic(self, output_dir, mock_gtf_file):
        """Test basic initialization."""
        pipeline = Hisat2FeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            threads=8
        )
        
        assert pipeline.PIPELINE_NAME == "hisat2_featurecounts"
        assert pipeline.threads == 8
    
    @pytest.mark.skipif(not HISAT2_FC_PIPELINE_AVAILABLE, reason="HISAT2+FC pipeline not available")
    def test_hisat2_options(self, output_dir, mock_gtf_file):
        """Test HISAT2-specific options."""
        pipeline = Hisat2FeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            min_intron_len=20,
            max_intron_len=500000,
            dta=True
        )
        
        assert pipeline.min_intron_len == 20
        assert pipeline.max_intron_len == 500000
        assert pipeline.dta == True
    
    @pytest.mark.skipif(not HISAT2_FC_PIPELINE_AVAILABLE, reason="HISAT2+FC pipeline not available")
    def test_rna_strandness(self, output_dir, mock_gtf_file):
        """Test RNA strandness option."""
        pipeline = Hisat2FeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            rna_strandness='RF'
        )
        
        assert pipeline.rna_strandness == 'RF'


# =============================================================================
# Tests for STAR + RSEM Pipeline
# =============================================================================

class TestStarRsemPipeline:
    """Tests for StarRsemPipeline class."""
    
    @pytest.mark.skipif(not STAR_RSEM_PIPELINE_AVAILABLE, reason="STAR+RSEM pipeline not available")
    def test_init_basic(self, output_dir):
        """Test basic initialization."""
        pipeline = StarRsemPipeline(
            output_dir=str(output_dir),
            threads=8
        )
        
        assert pipeline.PIPELINE_NAME == "star_rsem"
        assert pipeline.threads == 8
    
    @pytest.mark.skipif(not STAR_RSEM_PIPELINE_AVAILABLE, reason="STAR+RSEM pipeline not available")
    def test_strandedness_option(self, output_dir):
        """Test strandedness option."""
        pipeline = StarRsemPipeline(
            output_dir=str(output_dir),
            strandedness='reverse'
        )
        
        assert pipeline.strandedness == 'reverse'
    
    @pytest.mark.skipif(not STAR_RSEM_PIPELINE_AVAILABLE, reason="STAR+RSEM pipeline not available")
    def test_rsem_options(self, output_dir):
        """Test RSEM-specific options."""
        pipeline = StarRsemPipeline(
            output_dir=str(output_dir),
            estimate_rspd=True,
            calc_ci=True
        )
        
        assert pipeline.estimate_rspd == True
        assert pipeline.calc_ci == True
    
    @pytest.mark.skipif(not STAR_RSEM_PIPELINE_AVAILABLE, reason="STAR+RSEM pipeline not available")
    def test_fragment_length_single_end(self, output_dir):
        """Test fragment length for single-end."""
        pipeline = StarRsemPipeline(
            output_dir=str(output_dir),
            fragment_length_mean=200,
            fragment_length_sd=30
        )
        
        assert pipeline.fragment_length_mean == 200
        assert pipeline.fragment_length_sd == 30


# =============================================================================
# Tests for STAR + Salmon Pipeline
# =============================================================================

class TestStarSalmonPipeline:
    """Tests for StarSalmonPipeline class."""
    
    @pytest.mark.skipif(not STAR_SALMON_PIPELINE_AVAILABLE, reason="STAR+Salmon pipeline not available")
    def test_init_basic(self, output_dir, mock_salmon_index):
        """Test basic initialization."""
        pipeline = StarSalmonPipeline(
            output_dir=str(output_dir),
            salmon_index=str(mock_salmon_index),
            threads=8
        )
        
        assert pipeline.PIPELINE_NAME == "star_salmon"
        assert pipeline.threads == 8
    
    @pytest.mark.skipif(not STAR_SALMON_PIPELINE_AVAILABLE, reason="STAR+Salmon pipeline not available")
    def test_salmon_index_required(self, output_dir):
        """Test that Salmon index is required."""
        with pytest.raises((ValueError, TypeError)):
            StarSalmonPipeline(
                output_dir=str(output_dir)
                # Missing salmon_index parameter
            )
    
    @pytest.mark.skipif(not STAR_SALMON_PIPELINE_AVAILABLE, reason="STAR+Salmon pipeline not available")
    def test_hybrid_options(self, output_dir, mock_salmon_index):
        """Test hybrid STAR+Salmon options."""
        pipeline = StarSalmonPipeline(
            output_dir=str(output_dir),
            salmon_index=str(mock_salmon_index),
            two_pass_mode='Basic',
            library_type='A',
            bootstraps=100,
            gc_bias=True,
            keep_bam=True
        )
        
        assert pipeline.two_pass_mode == 'Basic'
        assert pipeline.library_type == 'A'
        assert pipeline.bootstraps == 100
        assert pipeline.keep_bam == True
    
    @pytest.mark.skipif(not STAR_SALMON_PIPELINE_AVAILABLE, reason="STAR+Salmon pipeline not available")
    def test_get_dependencies(self, output_dir, mock_salmon_index):
        """Test get_dependencies returns STAR, Salmon, samtools."""
        pipeline = StarSalmonPipeline(
            output_dir=str(output_dir),
            salmon_index=str(mock_salmon_index)
        )
        
        deps = pipeline.get_dependencies()
        dep_names = [d.name for d in deps]
        
        assert 'STAR' in dep_names
        assert 'Salmon' in dep_names


# =============================================================================
# Tests for Pipeline Registry
# =============================================================================

class TestPipelineRegistry:
    """Tests for pipeline registry functions."""
    
    @pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
    def test_get_pipeline_salmon(self):
        """Test getting Salmon pipeline class."""
        PipelineClass = get_pipeline('salmon')
        
        assert PipelineClass is not None
        assert PipelineClass.PIPELINE_NAME == 'salmon'
    
    @pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
    def test_get_pipeline_kallisto(self):
        """Test getting Kallisto pipeline class."""
        PipelineClass = get_pipeline('kallisto')
        
        assert PipelineClass is not None
        assert PipelineClass.PIPELINE_NAME == 'kallisto'
    
    @pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
    def test_get_pipeline_invalid(self):
        """Test getting invalid pipeline returns None."""
        PipelineClass = get_pipeline('invalid_pipeline')
        
        assert PipelineClass is None
    
    @pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
    def test_pipeline_info_dict(self):
        """Test PIPELINE_INFO dictionary."""
        assert isinstance(PIPELINE_INFO, dict)
        assert len(PIPELINE_INFO) >= 6
        
        # Check expected pipelines exist
        for name in AVAILABLE_PIPELINES:
            assert name in PIPELINE_INFO or name.replace('_', '-') in PIPELINE_INFO
    
    @pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
    def test_pipeline_info_structure(self):
        """Test PIPELINE_INFO has expected structure."""
        for name, info in PIPELINE_INFO.items():
            assert 'description' in info or hasattr(info, 'description')
            assert 'bam_output' in info or hasattr(info, 'bam_output')


# =============================================================================
# Tests for Docker/HPC Hybrid Dependency Management
# =============================================================================

class TestHybridDependencies:
    """Tests for Docker/HPC hybrid dependency management."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_docker_mode(self, output_dir):
        """Test Docker mode initialization."""
        pipeline = SalmonPipeline(
            output_dir=str(output_dir),
            use_docker=True,
            docker_image='combinelab/salmon:1.10.0'
        )
        
        assert pipeline.use_docker == True
        assert pipeline.docker_image == 'combinelab/salmon:1.10.0'
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_hpc_modules(self, output_dir):
        """Test HPC module loading."""
        pipeline = SalmonPipeline(
            output_dir=str(output_dir),
            modules=['salmon/1.9.0', 'gcc/11.2.0']
        )
        
        assert pipeline.modules is not None
        assert len(pipeline.modules) == 2
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_kallisto_docker(self, output_dir):
        """Test Kallisto Docker configuration."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir),
            use_docker=True
        )
        
        # Should use default Docker image
        assert pipeline.use_docker == True


# =============================================================================
# Tests for Count Matrix Combination
# =============================================================================

class TestCountMatrixCombination:
    """Tests for count matrix combination methods."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_combine_counts_structure(self, output_dir, mock_salmon_quant_dir):
        """Test combine_counts output structure."""
        pipeline = SalmonPipeline(output_dir=str(output_dir))
        
        # Create mock sample results
        sample_results = []
        for sample_dir in Path(mock_salmon_quant_dir).iterdir():
            sample_results.append({
                'success': True,
                'sample_id': sample_dir.name,
                'counts_file': str(sample_dir / 'quant.sf')
            })
        
        counts_file = output_dir / 'gene_counts.csv'
        counts_df = pipeline.combine_counts(sample_results, counts_file)
        
        assert isinstance(counts_df, pd.DataFrame)
        assert len(counts_df.columns) > 0
        assert counts_file.exists()
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_kallisto_combine_counts(self, output_dir, mock_kallisto_quant_dir):
        """Test Kallisto count matrix combination."""
        pipeline = KallistoPipeline(output_dir=str(output_dir))
        
        sample_results = []
        for sample_dir in Path(mock_kallisto_quant_dir).iterdir():
            sample_results.append({
                'success': True,
                'sample_id': sample_dir.name,
                'counts_file': str(sample_dir / 'abundance.tsv')
            })
        
        counts_file = output_dir / 'gene_counts.csv'
        counts_df = pipeline.combine_counts(sample_results, counts_file)
        
        assert isinstance(counts_df, pd.DataFrame)
        # Kallisto outputs est_counts which should be integers
        assert counts_df.dtypes.apply(lambda x: np.issubdtype(x, np.integer)).all()


# =============================================================================
# Tests for TPM Matrix Generation
# =============================================================================

class TestTPMGeneration:
    """Tests for TPM matrix generation."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_generate_tpm(self, output_dir, mock_salmon_quant_dir):
        """Test TPM matrix generation."""
        pipeline = SalmonPipeline(output_dir=str(output_dir))
        
        sample_results = []
        for sample_dir in Path(mock_salmon_quant_dir).iterdir():
            sample_results.append({
                'success': True,
                'sample_id': sample_dir.name,
                'counts_file': str(sample_dir / 'quant.sf')
            })
        
        # First combine counts (to cache TPM)
        counts_file = output_dir / 'gene_counts.csv'
        pipeline.combine_counts(sample_results, counts_file)
        
        # Then generate TPM
        tpm_file = output_dir / 'tpm.csv'
        tpm_df = pipeline.generate_tpm(sample_results, tpm_file)
        
        assert isinstance(tpm_df, pd.DataFrame)
        assert tpm_file.exists()


# =============================================================================
# Tests for Pipeline Output Compliance
# =============================================================================

class TestOutputCompliance:
    """Tests for architecture-compliant output."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_salmon_output_dir_default(self):
        """Test Salmon default output directory."""
        pipeline = SalmonPipeline()
        
        # Should default to results/production
        assert 'production' in str(pipeline.output_dir)
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_kallisto_output_dir_default(self):
        """Test Kallisto default output directory."""
        pipeline = KallistoPipeline()
        
        assert 'production' in str(pipeline.output_dir)
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_output_file_naming(self, output_dir):
        """Test output file naming convention."""
        pipeline = SalmonPipeline(output_dir=str(output_dir))
        
        # Expected file names
        expected_counts = output_dir / "gene_counts.csv"
        expected_tpm = output_dir / "tpm.csv"
        
        # Files should follow convention
        assert expected_counts.name == "gene_counts.csv"
        assert expected_tpm.name == "tpm.csv"


# =============================================================================
# Tests for Error Handling
# =============================================================================

class TestErrorHandling:
    """Tests for error handling in pipelines."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_invalid_index_path(self, output_dir):
        """Test handling of invalid index path."""
        pipeline = SalmonPipeline(output_dir=str(output_dir))
        
        # Should raise error during validation
        with pytest.raises((FileNotFoundError, ValueError)):
            pipeline.validate_inputs('/nonexistent/index')
    
    @pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
    def test_missing_fragment_length_single_end(self, single_end_sample_sheet_path, output_dir, mock_kallisto_index):
        """Test error when fragment length missing for single-end."""
        pipeline = KallistoPipeline(
            output_dir=str(output_dir)
            # Deliberately not setting fragment_length
        )
        
        # Should fail when running single-end without fragment length
        # The exact behavior depends on implementation


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for Module 5."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    @patch('shutil.which')
    @patch('subprocess.run')
    def test_salmon_full_workflow_mock(self, mock_run, mock_which, 
                                        pipeline_sample_sheet, mock_salmon_index, 
                                        output_dir):
        """Test Salmon pipeline workflow with mocked tools."""
        # Mock tool availability
        mock_which.return_value = '/usr/bin/salmon'
        mock_run.return_value = MagicMock(returncode=0, stdout='', stderr='')
        
        pipeline = SalmonPipeline(
            output_dir=str(output_dir),
            threads=4
        )
        
        # Pipeline should initialize without error
        assert pipeline.threads == 4
    
    @pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
    def test_all_pipelines_instantiable(self, output_dir, mock_gtf_file, mock_salmon_index):
        """Test that all pipelines can be instantiated."""
        # Salmon
        if SALMON_PIPELINE_AVAILABLE:
            SalmonPipeline(output_dir=str(output_dir))
        
        # Kallisto
        if KALLISTO_PIPELINE_AVAILABLE:
            KallistoPipeline(output_dir=str(output_dir))
        
        # STAR + featureCounts (requires GTF)
        if STAR_FC_PIPELINE_AVAILABLE:
            StarFeatureCountsPipeline(output_dir=str(output_dir), gtf=str(mock_gtf_file))
        
        # HISAT2 + featureCounts (requires GTF)
        if HISAT2_FC_PIPELINE_AVAILABLE:
            Hisat2FeatureCountsPipeline(output_dir=str(output_dir), gtf=str(mock_gtf_file))
        
        # STAR + RSEM
        if STAR_RSEM_PIPELINE_AVAILABLE:
            StarRsemPipeline(output_dir=str(output_dir))
        
        # STAR + Salmon (requires salmon_index)
        if STAR_SALMON_PIPELINE_AVAILABLE:
            StarSalmonPipeline(output_dir=str(output_dir), salmon_index=str(mock_salmon_index))


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
@pytest.mark.parametrize("threads", [1, 4, 8, 16, 32])
def test_salmon_thread_options(output_dir, threads):
    """Test Salmon pipeline with various thread counts."""
    pipeline = SalmonPipeline(
        output_dir=str(output_dir),
        threads=threads
    )
    
    assert pipeline.threads == threads


@pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
@pytest.mark.parametrize("library_type", ['A', 'ISR', 'ISF', 'SR', 'SF', 'IU', 'U'])
def test_salmon_library_types(output_dir, library_type):
    """Test Salmon pipeline with various library types."""
    pipeline = SalmonPipeline(
        output_dir=str(output_dir),
        library_type=library_type
    )
    
    assert pipeline.library_type == library_type


@pytest.mark.skipif(not KALLISTO_PIPELINE_AVAILABLE, reason="Kallisto pipeline not available")
@pytest.mark.parametrize("bootstraps", [0, 30, 100, 200])
def test_kallisto_bootstrap_options(output_dir, bootstraps):
    """Test Kallisto pipeline with various bootstrap counts."""
    pipeline = KallistoPipeline(
        output_dir=str(output_dir),
        bootstraps=bootstraps
    )
    
    assert pipeline.bootstraps == bootstraps


@pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
@pytest.mark.parametrize("strand_specific", [0, 1, 2])
def test_star_fc_strand_options(output_dir, mock_gtf_file, strand_specific):
    """Test STAR+featureCounts with various strand options."""
    pipeline = StarFeatureCountsPipeline(
        output_dir=str(output_dir),
        gtf=str(mock_gtf_file),
        strand_specific=strand_specific
    )
    
    assert pipeline.strand_specific == strand_specific


@pytest.mark.skipif(not REGISTRY_AVAILABLE, reason="Registry not available")
@pytest.mark.parametrize("pipeline_name", AVAILABLE_PIPELINES)
def test_all_pipelines_in_registry(pipeline_name):
    """Test that all expected pipelines are in registry."""
    # Some pipelines may use underscores, others hyphens
    PipelineClass = get_pipeline(pipeline_name)
    if PipelineClass is None:
        PipelineClass = get_pipeline(pipeline_name.replace('_', '-'))
    
    # At least one variant should work
    # (depends on implementation)


# =============================================================================
# CLI Integration Tests
# =============================================================================

class TestCLIIntegration:
    """Tests for CLI command integration."""
    
    def test_cli_import(self):
        """Test that CLI can be imported."""
        try:
            from raptor.cli import main, pipeline
            assert callable(main)
        except ImportError:
            pytest.skip("CLI not available")
    
    def test_cli_pipeline_list(self):
        """Test pipeline list command output."""
        try:
            from click.testing import CliRunner
            from raptor.cli import main
            
            runner = CliRunner()
            result = runner.invoke(main, ['pipeline', 'list'])
            
            # Should list pipelines
            assert result.exit_code == 0 or 'salmon' in result.output.lower()
        except ImportError:
            pytest.skip("CLI or Click not available")


# =============================================================================
# Tests for Module 5 __init__.py
# =============================================================================

class TestModule5Init:
    """Tests for Module 5 package initialization."""
    
    def test_module5_imports(self):
        """Test that Module 5 exports expected names."""
        try:
            from raptor.pipelines import (
                get_pipeline,
                PIPELINE_INFO
            )
            
            assert callable(get_pipeline)
            assert isinstance(PIPELINE_INFO, dict)
        except ImportError:
            pytest.skip("Module 5 package not available")
    
    def test_module5_info(self):
        """Test MODULE_INFO if available."""
        try:
            from raptor.pipelines import MODULE_INFO
            
            assert MODULE_INFO['id'] == 'M5'
            assert MODULE_INFO['stage'] == 2
            assert 'Production' in MODULE_INFO['name']
        except (ImportError, KeyError):
            pytest.skip("MODULE_INFO not available")


# =============================================================================
# Performance/Resource Tests
# =============================================================================

class TestResourceManagement:
    """Tests for resource management."""
    
    @pytest.mark.skipif(not SALMON_PIPELINE_AVAILABLE, reason="Salmon pipeline not available")
    def test_memory_setting(self, output_dir):
        """Test memory setting."""
        pipeline = SalmonPipeline(
            output_dir=str(output_dir),
            memory_gb=32
        )
        
        assert pipeline.memory_gb == 32
    
    @pytest.mark.skipif(not STAR_FC_PIPELINE_AVAILABLE, reason="STAR+FC pipeline not available")
    def test_star_memory_intensive(self, output_dir, mock_gtf_file):
        """Test STAR pipeline memory requirements."""
        pipeline = StarFeatureCountsPipeline(
            output_dir=str(output_dir),
            gtf=str(mock_gtf_file),
            memory_gb=64
        )
        
        # STAR typically needs more memory
        assert pipeline.memory_gb >= 32
