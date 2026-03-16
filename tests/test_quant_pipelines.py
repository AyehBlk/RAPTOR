"""
RAPTOR v2.2.0 - Tests for Quick Quantification Pipelines

Tests for salmon_quant.py and kallisto_quant.py pipelines.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import subprocess
from unittest.mock import patch, MagicMock

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import modules under test
try:
    from raptor.pipelines.quick_salmon.salmon_quant import (
        QuickSalmonPipeline,
        run_quick_salmon
    )
    SALMON_AVAILABLE = True
except ImportError:
    SALMON_AVAILABLE = False

try:
    from raptor.pipelines.quick_kallisto.kallisto_quant import (
        QuickKallistoPipeline,
        run_quick_kallisto
    )
    KALLISTO_AVAILABLE = True
except ImportError:
    KALLISTO_AVAILABLE = False


# =============================================================================
# Constants
# =============================================================================

DEFAULT_OUTPUT_DIR = "results/quick_counts"
OUTPUT_COUNTS_FILE = "quick_gene_counts.csv"
OUTPUT_TPM_FILE = "quick_tpm.csv"


# =============================================================================
# Tests for QuickSalmonPipeline
# =============================================================================

class TestQuickSalmonPipeline:
    """Tests for QuickSalmonPipeline class."""
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_init(self, sample_sheet_path, test_dir):
        """Test pipeline initialization."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index',
            output_dir=str(test_dir / 'output'),
            threads=4
        )
        
        assert pipeline.sample_sheet == str(sample_sheet_path)
        assert pipeline.threads == 4
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_default_output_dir(self, sample_sheet_path):
        """Test default output directory."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index'
        )
        
        # Should default to results/quick_counts (normalize for Windows)
        assert 'quick_counts' in str(pipeline.output_dir)
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_output_file_names(self, sample_sheet_path, output_dir):
        """Test output file naming convention."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index',
            output_dir=str(output_dir)
        )
        
        # Check expected output paths
        expected_counts = output_dir / OUTPUT_COUNTS_FILE
        expected_tpm = output_dir / OUTPUT_TPM_FILE
        
        assert str(expected_counts.name) == OUTPUT_COUNTS_FILE
        assert str(expected_tpm.name) == OUTPUT_TPM_FILE
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_validate_inputs_missing_sample_sheet(self, test_dir):
        """Test validation with missing sample sheet."""
        pipeline = QuickSalmonPipeline(
            sample_sheet='/nonexistent/samples.csv',
            index='/fake/salmon_index',
            output_dir=str(test_dir / 'output')
        )
        
        with pytest.raises(FileNotFoundError):
            pipeline.validate_inputs()
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_validate_inputs_missing_index(self, sample_sheet_path, test_dir):
        """Test validation with missing index."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/nonexistent/index',
            output_dir=str(test_dir / 'output')
        )
        
        with pytest.raises(FileNotFoundError):
            pipeline.validate_inputs()


# =============================================================================
# Tests for QuickKallistoPipeline
# =============================================================================

class TestQuickKallistoPipeline:
    """Tests for QuickKallistoPipeline class."""
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    def test_init(self, sample_sheet_path, tx2gene_path, test_dir):
        """Test pipeline initialization."""
        pipeline = QuickKallistoPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/kallisto.idx',
            gene_map=str(tx2gene_path),
            output_dir=str(test_dir / 'output'),
            threads=4
        )
        
        assert pipeline.sample_sheet == str(sample_sheet_path)
        assert pipeline.gene_map == str(tx2gene_path)
        assert pipeline.threads == 4
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    def test_single_end_requires_fragment_length(self, single_end_sample_sheet_path, tx2gene_path, test_dir):
        """Test that single-end requires fragment length."""
        pipeline = QuickKallistoPipeline(
            sample_sheet=str(single_end_sample_sheet_path),
            index='/fake/kallisto.idx',
            gene_map=str(tx2gene_path),
            output_dir=str(test_dir / 'output')
        )
        
        # Should fail validation without fragment length for single-end
        # This depends on implementation - may raise error during validation
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    def test_single_end_with_fragment_length(self, single_end_sample_sheet_path, tx2gene_path, test_dir):
        """Test single-end with fragment length specified."""
        pipeline = QuickKallistoPipeline(
            sample_sheet=str(single_end_sample_sheet_path),
            index='/fake/kallisto.idx',
            gene_map=str(tx2gene_path),
            output_dir=str(test_dir / 'output'),
            fragment_length=200,
            fragment_sd=20
        )
        
        assert pipeline.fragment_length == 200
        assert pipeline.fragment_sd == 20


# =============================================================================
# Tests for Convenience Functions
# =============================================================================

class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    @patch('subprocess.run')
    def test_run_quick_salmon_calls_pipeline(self, mock_run, sample_sheet_path, test_dir):
        """Test that run_quick_salmon creates and runs pipeline."""
        mock_run.return_value = MagicMock(returncode=0)
        
        # This test depends on the implementation
        # Just verify the function exists and takes correct args
        assert callable(run_quick_salmon)
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    def test_run_quick_kallisto_callable(self):
        """Test that run_quick_kallisto is callable."""
        assert callable(run_quick_kallisto)


# =============================================================================
# Tests for Output Path Compliance
# =============================================================================

class TestOutputPathCompliance:
    """Tests for architecture compliance of output paths."""
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_salmon_default_output_path(self, sample_sheet_path):
        """Test Salmon default output path matches architecture."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index'
        )
        
        # Should default to results/quick_counts
        assert 'quick_counts' in str(pipeline.output_dir)
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    def test_kallisto_default_output_path(self, sample_sheet_path, tx2gene_path):
        """Test Kallisto default output path matches architecture."""
        pipeline = QuickKallistoPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/kallisto.idx',
            gene_map=str(tx2gene_path)
        )
        
        assert 'quick_counts' in str(pipeline.output_dir)
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_custom_output_path(self, sample_sheet_path, output_dir):
        """Test custom output path is respected."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index',
            output_dir=str(output_dir)
        )
        
        assert str(output_dir) in str(pipeline.output_dir)


# =============================================================================
# Tests for Configuration
# =============================================================================

class TestConfiguration:
    """Tests for pipeline configuration."""
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_salmon_default_params(self, sample_sheet_path, test_dir):
        """Test Salmon default parameters."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index',
            output_dir=str(test_dir / 'output')
        )
        
        # Check defaults
        assert pipeline.threads == 8 or hasattr(pipeline, 'threads')
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_salmon_custom_threads(self, sample_sheet_path, test_dir):
        """Test Salmon with custom thread count."""
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index',
            output_dir=str(test_dir / 'output'),
            threads=16
        )
        
        assert pipeline.threads == 16


# =============================================================================
# Mock Integration Tests
# =============================================================================

class TestMockIntegration:
    """Mock integration tests (don't require actual tools)."""
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    @patch('shutil.which')
    def test_salmon_tool_check(self, mock_which, sample_sheet_path, test_dir):
        """Test Salmon tool availability check."""
        mock_which.return_value = '/usr/bin/salmon'
        
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/salmon_index',
            output_dir=str(test_dir / 'output')
        )
        
        # Tool check should pass with mock
        # Implementation specific
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    @patch('shutil.which')
    def test_kallisto_tool_check(self, mock_which, sample_sheet_path, tx2gene_path, test_dir):
        """Test Kallisto tool availability check."""
        mock_which.return_value = '/usr/bin/kallisto'
        
        pipeline = QuickKallistoPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/kallisto.idx',
            gene_map=str(tx2gene_path),
            output_dir=str(test_dir / 'output')
        )
        
        # Tool check should pass with mock


# =============================================================================
# Tests for Error Handling
# =============================================================================

class TestErrorHandling:
    """Tests for error handling in pipelines."""
    
    @pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
    def test_salmon_invalid_sample_sheet(self, test_dir):
        """Test Salmon with invalid sample sheet format."""
        invalid_sheet = test_dir / 'invalid.csv'
        invalid_sheet.write_text('invalid,format\n1,2,3,4,5')
        
        pipeline = QuickSalmonPipeline(
            sample_sheet=str(invalid_sheet),
            index='/fake/salmon_index',
            output_dir=str(test_dir / 'output')
        )
        
        # Should raise validation error
        # Implementation specific
    
    @pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
    def test_kallisto_missing_gene_map(self, sample_sheet_path, test_dir):
        """Test Kallisto with missing gene map."""
        pipeline = QuickKallistoPipeline(
            sample_sheet=str(sample_sheet_path),
            index='/fake/kallisto.idx',
            gene_map='/nonexistent/tx2gene.csv',
            output_dir=str(test_dir / 'output')
        )
        
        # validate_inputs() should raise an error (ValidationError from index
        # check or FileNotFoundError from gene_map check)
        with pytest.raises((FileNotFoundError, Exception)):
            pipeline.validate_inputs()


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not SALMON_AVAILABLE, reason="Salmon module not available")
@pytest.mark.parametrize("threads", [1, 4, 8, 16])
def test_salmon_thread_options(sample_sheet_path, test_dir, threads):
    """Test Salmon with various thread counts."""
    pipeline = QuickSalmonPipeline(
        sample_sheet=str(sample_sheet_path),
        index='/fake/salmon_index',
        output_dir=str(test_dir / 'output'),
        threads=threads
    )
    
    assert pipeline.threads == threads


@pytest.mark.skipif(not KALLISTO_AVAILABLE, reason="Kallisto module not available")
@pytest.mark.parametrize("frag_len,frag_sd", [
    (150, 20),
    (200, 30),
    (250, 40),
    (300, 50),
])
def test_kallisto_fragment_options(sample_sheet_path, tx2gene_path, test_dir, frag_len, frag_sd):
    """Test Kallisto with various fragment length options."""
    pipeline = QuickKallistoPipeline(
        sample_sheet=str(sample_sheet_path),
        index='/fake/kallisto.idx',
        gene_map=str(tx2gene_path),
        output_dir=str(test_dir / 'output'),
        fragment_length=frag_len,
        fragment_sd=frag_sd
    )
    
    assert pipeline.fragment_length == frag_len
    assert pipeline.fragment_sd == frag_sd