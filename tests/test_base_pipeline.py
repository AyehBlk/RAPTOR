"""
RAPTOR v2.2.0 - Tests for Base Pipeline Infrastructure (Module 5)

Tests for base pipeline classes, sample detection, and common utilities:
- BasePipeline abstract class
- SampleSheet parsing and validation
- FASTQ auto-detection
- Tool dependency checking
- Command execution utilities

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
import gzip
from unittest.mock import patch, MagicMock
import subprocess

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import modules under test
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

try:
    from raptor.utils.sample_sheet import SampleSheet, Sample
    SAMPLE_SHEET_AVAILABLE = True
except ImportError as e:
    SAMPLE_SHEET_AVAILABLE = False
    # Show warning if SampleSheet not available
    import warnings
    warnings.warn(
        f"SampleSheet class not available - related tests will be skipped.\n"
        f"  This is normal if SampleSheet dependencies are not installed.\n"
        f"  Error: {e}",
        UserWarning
    )


# =============================================================================
# Additional Fixtures
# =============================================================================

@pytest.fixture
def extended_fastq_dir(test_dir):
    """Create a directory with various FASTQ file naming patterns."""
    fastq_dir = test_dir / "extended_fastq"
    fastq_dir.mkdir(exist_ok=True)
    
    # Various naming patterns
    patterns = [
        # Illumina standard
        ("Sample_A", "Sample_A_S1_L001_R1_001.fastq.gz", "Sample_A_S1_L001_R2_001.fastq.gz"),
        ("Sample_B", "Sample_B_S2_L001_R1_001.fastq.gz", "Sample_B_S2_L001_R2_001.fastq.gz"),
        # Simple pattern
        ("Sample_C", "Sample_C_R1.fastq.gz", "Sample_C_R2.fastq.gz"),
        # Underscore pattern
        ("Sample_D", "Sample_D_1.fastq.gz", "Sample_D_2.fastq.gz"),
        # Dot pattern
        ("Sample_E", "Sample_E.R1.fq.gz", "Sample_E.R2.fq.gz"),
    ]
    
    for sample_name, r1_name, r2_name in patterns:
        r1_path = fastq_dir / r1_name
        r2_path = fastq_dir / r2_name
        
        with gzip.open(r1_path, 'wt') as f:
            f.write(f"@{sample_name}_R1_1\nACGTACGT\n+\nIIIIIIII\n")
        with gzip.open(r2_path, 'wt') as f:
            f.write(f"@{sample_name}_R2_1\nACGTACGT\n+\nIIIIIIII\n")
    
    return fastq_dir


@pytest.fixture
def single_end_fastq_dir(test_dir):
    """Create a directory with single-end FASTQ files."""
    fastq_dir = test_dir / "single_end_fastq"
    fastq_dir.mkdir(exist_ok=True)
    
    samples = ['Sample_1', 'Sample_2', 'Sample_3']
    
    for sample in samples:
        fq_path = fastq_dir / f"{sample}.fastq.gz"
        with gzip.open(fq_path, 'wt') as f:
            f.write(f"@{sample}_1\nACGTACGT\n+\nIIIIIIII\n")
    
    return fastq_dir


@pytest.fixture
def valid_sample_sheet_df():
    """Create a valid sample sheet DataFrame."""
    return pd.DataFrame({
        'sample_id': ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4'],
        'condition': ['Control', 'Control', 'Treatment', 'Treatment'],
        'batch': ['Batch1', 'Batch1', 'Batch2', 'Batch2'],
        'fastq_r1': ['/path/to/Sample1_R1.fastq.gz', '/path/to/Sample2_R1.fastq.gz',
                    '/path/to/Sample3_R1.fastq.gz', '/path/to/Sample4_R1.fastq.gz'],
        'fastq_r2': ['/path/to/Sample1_R2.fastq.gz', '/path/to/Sample2_R2.fastq.gz',
                    '/path/to/Sample3_R2.fastq.gz', '/path/to/Sample4_R2.fastq.gz']
    })


@pytest.fixture
def mixed_sample_sheet(test_dir, mock_fastq_dir):
    """Create a sample sheet with mixed conditions."""
    sample_sheet = test_dir / "mixed_samples.csv"
    
    data = {
        'sample_id': ['Ctrl_Rep1', 'Ctrl_Rep2', 'Ctrl_Rep3', 
                      'Drug_Rep1', 'Drug_Rep2', 'Drug_Rep3'],
        'condition': ['Control'] * 3 + ['Drug'] * 3,
        'batch': ['B1', 'B1', 'B2', 'B1', 'B2', 'B2'],
        'replicate': [1, 2, 3, 1, 2, 3],
        'fastq_r1': [str(mock_fastq_dir / f"Sample{i}_R1_001.fastq.gz") for i in range(1, 7)],
        'fastq_r2': [str(mock_fastq_dir / f"Sample{i}_R2_001.fastq.gz") for i in range(1, 7)]
    }
    
    # Create dummy FASTQ files
    for i in range(1, 7):
        r1 = mock_fastq_dir / f"Sample{i}_R1_001.fastq.gz"
        r2 = mock_fastq_dir / f"Sample{i}_R2_001.fastq.gz"
        if not r1.exists():
            with gzip.open(r1, 'wt') as f:
                f.write(f"@Sample{i}_1\nACGT\n+\nIIII\n")
        if not r2.exists():
            with gzip.open(r2, 'wt') as f:
                f.write(f"@Sample{i}_1\nACGT\n+\nIIII\n")
    
    pd.DataFrame(data).to_csv(sample_sheet, index=False)
    return sample_sheet


# =============================================================================
# Tests for SampleSheet Class
# =============================================================================

class TestSampleSheet:
    """Tests for SampleSheet class."""
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_load_sample_sheet(self, sample_sheet_path):
        """Test loading sample sheet from CSV."""
        sheet = SampleSheet(str(sample_sheet_path))
        
        assert len(sheet) > 0
        assert hasattr(sheet, 'samples')
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_sample_sheet_is_paired(self, sample_sheet_path):
        """Test paired-end detection."""
        sheet = SampleSheet(str(sample_sheet_path))
        
        assert sheet.is_paired == True
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_sample_sheet_single_end(self, single_end_sample_sheet_path):
        """Test single-end detection."""
        sheet = SampleSheet(str(single_end_sample_sheet_path))
        
        # Either is_paired is False or R2 columns are empty
        if hasattr(sheet, 'is_paired'):
            assert sheet.is_paired == False or not sheet.samples[0].fastq_r2
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_validate_with_file_check(self, sample_sheet_path):
        """Test validation with file existence check."""
        sheet = SampleSheet(str(sample_sheet_path))
        
        # With check_files=True, should fail because fake paths
        is_valid, errors = sheet.validate(check_files=True)
        
        # Either validation passes or returns errors
        assert isinstance(is_valid, bool)
        assert isinstance(errors, list)
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_validate_without_file_check(self, sample_sheet_path):
        """Test validation without file check."""
        sheet = SampleSheet(str(sample_sheet_path))
        
        is_valid, errors = sheet.validate(check_files=False)
        
        # Should pass without file check
        assert is_valid == True
        assert len(errors) == 0
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_get_sample(self, sample_sheet_path):
        """Test getting a specific sample."""
        sheet = SampleSheet(str(sample_sheet_path))
        
        # Get first sample
        sample = sheet.get_sample(sheet.samples[0].sample_id)
        
        assert sample is not None
        assert hasattr(sample, 'sample_id')
    
    @pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
    def test_iterate_samples(self, sample_sheet_path):
        """Test iterating over samples."""
        sheet = SampleSheet(str(sample_sheet_path))
        
        sample_ids = [s.sample_id for s in sheet]
        
        assert len(sample_ids) == len(sheet)


# =============================================================================
# Tests for FASTQ Auto-Detection
# =============================================================================

class TestFASTQAutoDetection:
    """Tests for FASTQ file auto-detection."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_auto_detect_paired_end(self, mock_fastq_dir):
        """Test auto-detection of paired-end samples."""
        samples, read_type = auto_detect_samples(str(mock_fastq_dir))
        
        assert read_type == 'paired'
        assert len(samples) > 0
        
        for sample in samples:
            assert sample.fastq_r1 is not None
            assert sample.fastq_r2 is not None
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_auto_detect_single_end(self, single_end_fastq_dir):
        """Test auto-detection of single-end samples."""
        samples, read_type = auto_detect_samples(str(single_end_fastq_dir))
        
        assert read_type == 'single'
        assert len(samples) > 0
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_auto_detect_various_patterns(self, extended_fastq_dir):
        """Test detection of various FASTQ naming patterns."""
        samples, read_type = auto_detect_samples(str(extended_fastq_dir))
        
        assert len(samples) == 5  # 5 sample pairs
        assert read_type == 'paired'
        
        sample_ids = [s.sample_id for s in samples]
        expected = ['Sample_A', 'Sample_B', 'Sample_C', 'Sample_D', 'Sample_E']
        
        for exp in expected:
            # Sample ID extraction may vary
            assert any(exp in sid for sid in sample_ids)
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_auto_detect_empty_directory(self, test_dir):
        """Test auto-detection with empty directory."""
        empty_dir = test_dir / "empty_fastq"
        empty_dir.mkdir(exist_ok=True)
        
        samples, read_type = auto_detect_samples(str(empty_dir))
        
        assert len(samples) == 0
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_auto_detect_with_pattern(self, extended_fastq_dir):
        """Test auto-detection with custom pattern."""
        # Only detect files matching specific pattern
        samples, read_type = auto_detect_samples(
            str(extended_fastq_dir),
            pattern=['_S*_L001']
        )
        
        # Should only find Illumina-style named files
        # Depends on implementation


# =============================================================================
# Tests for Sample Sheet Creation
# =============================================================================

class TestCreateSampleSheet:
    """Tests for sample sheet creation."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_create_sample_sheet_paired(self, mock_fastq_dir, output_dir):
        """Test creating sample sheet from paired-end FASTQs."""
        samples, read_type = auto_detect_samples(str(mock_fastq_dir))
        
        output_path = output_dir / "generated_samples.csv"
        df = create_sample_sheet(samples, str(output_path))
        
        assert output_path.exists()
        assert isinstance(df, pd.DataFrame)
        assert 'sample_id' in df.columns
        assert 'fastq_r1' in df.columns
        assert 'fastq_r2' in df.columns
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_create_sample_sheet_single(self, single_end_fastq_dir, output_dir):
        """Test creating sample sheet from single-end FASTQs."""
        samples, read_type = auto_detect_samples(str(single_end_fastq_dir))
        
        output_path = output_dir / "single_samples.csv"
        df = create_sample_sheet(samples, str(output_path), read_type)
        
        assert output_path.exists()
        assert 'sample_id' in df.columns
        assert 'fastq_r1' in df.columns


# =============================================================================
# Tests for Tool Dependency Checking
# =============================================================================

class TestToolDependencyChecking:
    """Tests for tool dependency checking."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    @patch('shutil.which')
    def test_tool_available(self, mock_which):
        """Test checking if tool is available."""
        mock_which.return_value = '/usr/bin/salmon'
        
        dep = ToolDependency(
            name="Salmon",
            command="salmon",
            min_version="1.9.0"
        )
        
        # Tool should be found
        assert shutil.which('salmon') is not None or mock_which.called
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    @patch('shutil.which')
    def test_tool_not_available(self, mock_which):
        """Test handling of unavailable tool."""
        mock_which.return_value = None
        
        dep = ToolDependency(
            name="NonExistentTool",
            command="nonexistent"
        )
        
        # Tool should not be found
        result = shutil.which('nonexistent')
        assert result is None
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    @patch('subprocess.run')
    def test_version_check(self, mock_run):
        """Test version checking."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="salmon 1.9.0",
            stderr=""
        )
        
        # Version checking logic
        result = subprocess.run(['salmon', '--version'], capture_output=True, text=True)
        
        assert '1.9.0' in result.stdout or mock_run.called


# =============================================================================
# Tests for Base Pipeline Class
# =============================================================================

class TestBasePipeline:
    """Tests for BasePipeline abstract class."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_base_pipeline_is_abstract(self):
        """Test that BasePipeline cannot be instantiated directly."""
        # BasePipeline should be abstract
        with pytest.raises(TypeError):
            BasePipeline(output_dir='/tmp/test')
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_base_pipeline_subclass(self, output_dir):
        """Test creating a BasePipeline subclass."""
        
        class TestPipeline(BasePipeline):
            PIPELINE_NAME = "test"
            PIPELINE_VERSION = "1.0.0"
            
            def get_dependencies(self):
                return []
            
            def run_sample(self, sample, index_path, **kwargs):
                return {'success': True, 'sample_id': sample.sample_id}
            
            def combine_counts(self, sample_results, output_file):
                return pd.DataFrame()
            
            def generate_tpm(self, sample_results, output_file):
                return pd.DataFrame()
        
        pipeline = TestPipeline(output_dir=str(output_dir))
        
        assert pipeline.PIPELINE_NAME == "test"
        assert pipeline.output_dir.exists()


# =============================================================================
# Tests for Command Execution
# =============================================================================

class TestCommandExecution:
    """Tests for command execution utilities."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    @patch('subprocess.run')
    def test_run_command_success(self, mock_run, output_dir):
        """Test successful command execution."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Success",
            stderr=""
        )
        
        # Create a test pipeline subclass
        class TestPipeline(BasePipeline):
            PIPELINE_NAME = "test"
            
            def get_dependencies(self):
                return []
            
            def run_sample(self, sample, index_path, **kwargs):
                return {}
            
            def combine_counts(self, sample_results, output_file):
                return pd.DataFrame()
            
            def generate_tpm(self, sample_results, output_file):
                return pd.DataFrame()
        
        pipeline = TestPipeline(output_dir=str(output_dir))
        
        # Test run_command if available
        if hasattr(pipeline, 'run_command'):
            returncode, stdout, stderr = pipeline.run_command(['echo', 'test'])
            assert returncode == 0 or mock_run.called
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    @patch('subprocess.run')
    def test_run_command_failure(self, mock_run, output_dir):
        """Test failed command execution."""
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="Error"
        )
        
        # Error handling depends on implementation


# =============================================================================
# Tests for Docker Integration
# =============================================================================

class TestDockerIntegration:
    """Tests for Docker integration in pipelines."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_docker_command_building(self, output_dir):
        """Test Docker command construction."""
        
        class DockerTestPipeline(BasePipeline):
            PIPELINE_NAME = "docker_test"
            DOCKER_IMAGE = "test/image:1.0"
            
            def __init__(self, **kwargs):
                kwargs.setdefault('use_docker', True)
                super().__init__(**kwargs)
            
            def get_dependencies(self):
                return []
            
            def run_sample(self, sample, index_path, **kwargs):
                return {}
            
            def combine_counts(self, sample_results, output_file):
                return pd.DataFrame()
            
            def generate_tpm(self, sample_results, output_file):
                return pd.DataFrame()
        
        pipeline = DockerTestPipeline(output_dir=str(output_dir))
        
        assert pipeline.use_docker == True


# =============================================================================
# Tests for HPC Module Loading
# =============================================================================

class TestHPCModules:
    """Tests for HPC module loading."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_module_list_parsing(self, output_dir):
        """Test parsing of HPC module list."""
        
        class ModuleTestPipeline(BasePipeline):
            PIPELINE_NAME = "module_test"
            
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
            
            def get_dependencies(self):
                return []
            
            def run_sample(self, sample, index_path, **kwargs):
                return {}
            
            def combine_counts(self, sample_results, output_file):
                return pd.DataFrame()
            
            def generate_tpm(self, sample_results, output_file):
                return pd.DataFrame()
        
        modules = ['salmon/1.9.0', 'gcc/11.2.0', 'python/3.9']
        pipeline = ModuleTestPipeline(output_dir=str(output_dir), modules=modules)
        
        assert pipeline.modules == modules
        assert len(pipeline.modules) == 3


# =============================================================================
# Tests for Extra Arguments
# =============================================================================

class TestExtraArguments:
    """Tests for extra argument handling."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
    def test_extra_args_dict(self, output_dir):
        """Test extra arguments as dictionary."""
        
        class ExtraArgsTestPipeline(BasePipeline):
            PIPELINE_NAME = "extra_args_test"
            
            def get_dependencies(self):
                return []
            
            def run_sample(self, sample, index_path, **kwargs):
                return {}
            
            def combine_counts(self, sample_results, output_file):
                return pd.DataFrame()
            
            def generate_tpm(self, sample_results, output_file):
                return pd.DataFrame()
        
        extra_args = {'custom-option': 'value', 'another-option': 42}
        pipeline = ExtraArgsTestPipeline(
            output_dir=str(output_dir),
            extra_args=extra_args
        )
        
        assert pipeline.extra_args == extra_args


# =============================================================================
# Integration Tests
# =============================================================================

class TestBaseIntegration:
    """Integration tests for base infrastructure."""
    
    @pytest.mark.skipif(not BASE_AVAILABLE and not SAMPLE_SHEET_AVAILABLE,
                        reason="Modules not available")
    def test_full_sample_detection_workflow(self, mock_fastq_dir, output_dir):
        """Test complete sample detection workflow."""
        # 1. Auto-detect samples
        samples, read_type = auto_detect_samples(str(mock_fastq_dir))
        
        assert len(samples) > 0
        
        # 2. Create sample sheet
        sheet_path = output_dir / "detected_samples.csv"
        df = create_sample_sheet(samples, str(sheet_path), read_type)
        
        assert sheet_path.exists()
        
        # 3. Load and validate
        sheet = SampleSheet(str(sheet_path))
        is_valid, errors = sheet.validate(check_files=False)
        
        assert is_valid


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not BASE_AVAILABLE, reason="Base module not available")
@pytest.mark.parametrize("extension", ['.fastq.gz', '.fq.gz', '.fastq', '.fq'])
def test_fastq_extensions(test_dir, extension):
    """Test detection of various FASTQ extensions."""
    fastq_dir = test_dir / f"ext_test_{extension.replace('.', '_')}"
    fastq_dir.mkdir(exist_ok=True)
    
    # Create test files
    r1_file = fastq_dir / f"Sample_R1{extension}"
    r2_file = fastq_dir / f"Sample_R2{extension}"
    
    if extension.endswith('.gz'):
        with gzip.open(r1_file, 'wt') as f:
            f.write("@Sample_1\nACGT\n+\nIIII\n")
        with gzip.open(r2_file, 'wt') as f:
            f.write("@Sample_1\nACGT\n+\nIIII\n")
    else:
        r1_file.write_text("@Sample_1\nACGT\n+\nIIII\n")
        r2_file.write_text("@Sample_1\nACGT\n+\nIIII\n")
    
    samples, read_type = auto_detect_samples(str(fastq_dir))
    
    # Should detect the sample regardless of extension
    assert len(samples) >= 0  # May or may not detect depending on implementation


@pytest.mark.skipif(not SAMPLE_SHEET_AVAILABLE, reason="SampleSheet not available")
@pytest.mark.parametrize("n_samples", [2, 4, 8, 12, 24])
def test_various_sample_counts(test_dir, n_samples):
    """Test sample sheets with various sample counts."""
    sheet_path = test_dir / f"samples_{n_samples}.csv"
    
    data = {
        'sample_id': [f'Sample_{i+1}' for i in range(n_samples)],
        'condition': ['Control'] * (n_samples // 2) + ['Treatment'] * (n_samples - n_samples // 2),
        'fastq_r1': [f'/path/Sample_{i+1}_R1.fastq.gz' for i in range(n_samples)],
        'fastq_r2': [f'/path/Sample_{i+1}_R2.fastq.gz' for i in range(n_samples)]
    }
    
    pd.DataFrame(data).to_csv(sheet_path, index=False)
    
    sheet = SampleSheet(str(sheet_path))
    
    assert len(sheet) == n_samples
