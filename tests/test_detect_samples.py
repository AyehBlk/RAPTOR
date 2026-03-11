"""
RAPTOR v2.2.0 - Tests for detect_samples.py

Tests for FASTQ auto-detection functionality.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import pytest
import pandas as pd
from pathlib import Path
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import module under test
try:
    from raptor.pipelines.quick_salmon.detect_samples import (
        FastqAutoDetector,
        auto_detect_samples,
        create_sample_sheet
    )
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False


# =============================================================================
# Tests for FastqAutoDetector
# =============================================================================

class TestFastqAutoDetector:
    """Tests for FastqAutoDetector class."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init(self, mock_fastq_dir):
        """Test detector initialization."""
        detector = FastqAutoDetector(str(mock_fastq_dir))
        assert detector.input_dir == mock_fastq_dir
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_find_fastq_files(self, mock_fastq_dir):
        """Test finding FASTQ files."""
        detector = FastqAutoDetector(str(mock_fastq_dir))
        files = detector.find_fastq_files()
        
        assert len(files) > 0
        assert all(f.endswith('.fastq.gz') for f in files)
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_classify_paired_end(self, mock_fastq_dir):
        """Test classification of paired-end files."""
        detector = FastqAutoDetector(str(mock_fastq_dir))
        files = detector.find_fastq_files()
        r1_files, r2_files, unknown = detector.classify_files(files)
        
        assert len(r1_files) > 0
        assert len(r2_files) > 0
        assert len(r1_files) == len(r2_files)
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_extract_sample_name(self, mock_fastq_dir):
        """Test sample name extraction from filenames."""
        detector = FastqAutoDetector(str(mock_fastq_dir))
        
        test_cases = [
            ('Control_1_R1_001.fastq.gz', 'Control_1'),
            ('Control_1_R2_001.fastq.gz', 'Control_1'),
            ('Sample_R1.fastq.gz', 'Sample'),
            ('Sample_1.fastq.gz', 'Sample'),
            ('Test_Sample.R1.fq.gz', 'Test_Sample'),
        ]
        
        for filename, expected in test_cases:
            result = detector.extract_sample_name(filename)
            assert result == expected, f"Failed for {filename}: got {result}"
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_detect_paired_end(self, mock_fastq_dir):
        """Test full detection of paired-end samples."""
        detector = FastqAutoDetector(str(mock_fastq_dir))
        samples, read_type = detector.detect()
        
        assert read_type == 'paired'
        assert len(samples) > 0
        
        for sample in samples:
            assert sample.sample_id is not None
            assert sample.fastq_r1 is not None
            assert sample.fastq_r2 is not None


class TestAutoDetectSamples:
    """Tests for auto_detect_samples function."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_auto_detect_paired(self, mock_fastq_dir):
        """Test auto-detection of paired-end samples."""
        samples, read_type = auto_detect_samples(str(mock_fastq_dir))
        
        assert read_type == 'paired'
        assert len(samples) > 0
    
    def test_empty_directory(self, test_dir):
        """Test handling of empty directory."""
        empty_dir = test_dir / "empty_fastq"
        empty_dir.mkdir(exist_ok=True)
        
        if IMPORTS_AVAILABLE:
            samples, read_type = auto_detect_samples(str(empty_dir))
            assert len(samples) == 0


class TestCreateSampleSheet:
    """Tests for sample sheet creation."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_create_sample_sheet(self, mock_fastq_dir, output_dir):
        """Test sample sheet CSV creation."""
        samples, read_type = auto_detect_samples(str(mock_fastq_dir))
        
        output_path = output_dir / "test_samples.csv"
        df = create_sample_sheet(samples, str(output_path), read_type)
        
        assert output_path.exists()
        assert len(df) == len(samples)
        assert 'sample_id' in df.columns
        assert 'fastq_r1' in df.columns
        assert 'fastq_r2' in df.columns


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    def test_nonexistent_directory(self):
        """Test handling of non-existent directory."""
        if IMPORTS_AVAILABLE:
            detector = FastqAutoDetector('/nonexistent/path')
            files = detector.find_fastq_files()
            assert len(files) == 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_mixed_file_types(self, test_dir):
        """Test handling of mixed file types in directory."""
        mixed_dir = test_dir / "mixed"
        mixed_dir.mkdir(exist_ok=True)
        
        # Create various file types
        (mixed_dir / "sample.txt").touch()
        (mixed_dir / "sample.fastq.gz").touch()
        (mixed_dir / "sample.bam").touch()
        
        detector = FastqAutoDetector(str(mixed_dir))
        files = detector.find_fastq_files()
        
        assert len(files) == 1
        assert files[0].endswith('.fastq.gz')


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
@pytest.mark.parametrize("filename,expected_sample", [
    ("Sample1_S1_L001_R1_001.fastq.gz", "Sample1_S1_L001"),
    ("Control_R1.fastq.gz", "Control"),
    ("Treatment.R1.fq.gz", "Treatment"),
    ("Sample_001_R1.fq", "Sample_001"),
    ("test_sample_1.fastq.gz", "test_sample"),
])
def test_sample_name_extraction(filename, expected_sample):
    """Parametrized test for sample name extraction."""
    from raptor.pipelines.quick_salmon.detect_samples import FastqAutoDetector
    
    detector = FastqAutoDetector("/fake/path")
    result = detector.extract_sample_name(filename)
    assert result == expected_sample
