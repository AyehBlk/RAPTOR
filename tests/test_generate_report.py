"""
RAPTOR v2.2.0 - Tests for generate_report.py

Tests for QC report generation.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import module under test
try:
    from raptor.pipelines.quick_salmon.scripts.generate_report import (
        load_salmon_meta,
        load_kallisto_run_info,
        collect_quant_stats,
        generate_html_report,
        OUTPUT_COUNTS_FILE
    )
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False


# =============================================================================
# Tests for Metadata Loading
# =============================================================================

class TestLoadMetadata:
    """Tests for loading tool-specific metadata."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_load_salmon_meta(self, mock_salmon_quant_dir):
        """Test loading Salmon meta_info.json."""
        sample_dirs = list(Path(mock_salmon_quant_dir).iterdir())
        sample_dir = sample_dirs[0]
        
        meta = load_salmon_meta(str(sample_dir))
        
        assert isinstance(meta, dict)
        assert 'salmon_version' in meta
        assert 'num_processed' in meta
        assert 'percent_mapped' in meta
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_load_kallisto_run_info(self, mock_kallisto_quant_dir):
        """Test loading Kallisto run_info.json."""
        sample_dirs = list(Path(mock_kallisto_quant_dir).iterdir())
        sample_dir = sample_dirs[0]
        
        run_info = load_kallisto_run_info(str(sample_dir))
        
        assert isinstance(run_info, dict)
        assert 'kallisto_version' in run_info
        assert 'n_processed' in run_info
        assert 'p_pseudoaligned' in run_info
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_missing_metadata_returns_empty(self, test_dir):
        """Test that missing metadata returns empty dict."""
        empty_dir = test_dir / "empty_sample"
        empty_dir.mkdir(exist_ok=True)
        
        meta = load_salmon_meta(str(empty_dir))
        assert meta == {}
        
        run_info = load_kallisto_run_info(str(empty_dir))
        assert run_info == {}


# =============================================================================
# Tests for Statistics Collection
# =============================================================================

class TestCollectQuantStats:
    """Tests for quantification statistics collection."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_collect_salmon_stats(self, mock_salmon_quant_dir):
        """Test collecting Salmon statistics."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        
        assert isinstance(stats_df, pd.DataFrame)
        assert len(stats_df) == 6  # 6 samples
        assert 'sample_id' in stats_df.columns
        assert 'n_processed' in stats_df.columns
        assert 'mapping_rate' in stats_df.columns
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_collect_kallisto_stats(self, mock_kallisto_quant_dir):
        """Test collecting Kallisto statistics."""
        stats_df = collect_quant_stats(str(mock_kallisto_quant_dir), tool='kallisto')
        
        assert isinstance(stats_df, pd.DataFrame)
        assert len(stats_df) == 6
        assert 'sample_id' in stats_df.columns
        assert 'mapping_rate' in stats_df.columns
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_stats_have_expected_values(self, mock_salmon_quant_dir):
        """Test that statistics have reasonable values."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        
        # All samples should have data
        assert stats_df['n_processed'].sum() > 0
        
        # Mapping rate should be between 0 and 100
        assert (stats_df['mapping_rate'] >= 0).all()
        assert (stats_df['mapping_rate'] <= 100).all()


# =============================================================================
# Tests for HTML Report Generation
# =============================================================================

class TestGenerateHtmlReport:
    """Tests for HTML report generation."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_generate_salmon_report(self, mock_salmon_quant_dir, sample_count_matrix, output_dir):
        """Test generating Salmon QC report."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='salmon')
        
        assert output_path.exists()
        
        # Check HTML content
        content = output_path.read_text()
        assert '<!DOCTYPE html>' in content
        assert 'RAPTOR' in content
        assert 'Salmon' in content
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_generate_kallisto_report(self, mock_kallisto_quant_dir, sample_count_matrix, output_dir):
        """Test generating Kallisto QC report."""
        stats_df = collect_quant_stats(str(mock_kallisto_quant_dir), tool='kallisto')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='kallisto')
        
        assert output_path.exists()
        content = output_path.read_text()
        assert 'Kallisto' in content
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_report_contains_sample_table(self, mock_salmon_quant_dir, sample_count_matrix, output_dir):
        """Test that report contains per-sample statistics table."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='salmon')
        
        content = output_path.read_text()
        
        # Check for sample names in table
        for sample in sample_count_matrix.columns:
            assert sample in content
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_report_contains_summary_stats(self, mock_salmon_quant_dir, sample_count_matrix, output_dir):
        """Test that report contains summary statistics."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='salmon')
        
        content = output_path.read_text()
        
        # Check for summary elements
        assert 'Samples' in content
        assert 'Genes' in content
        assert 'Library Size' in content


# =============================================================================
# Tests for Report Content Quality
# =============================================================================

class TestReportContentQuality:
    """Tests for report content quality and accuracy."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_report_has_correct_sample_count(self, mock_salmon_quant_dir, sample_count_matrix, output_dir):
        """Test that report shows correct sample count."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='salmon')
        
        content = output_path.read_text()
        
        # Should show 6 samples
        n_samples = len(sample_count_matrix.columns)
        assert str(n_samples) in content
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_report_has_next_steps(self, mock_salmon_quant_dir, sample_count_matrix, output_dir):
        """Test that report includes next steps."""
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='salmon')
        
        content = output_path.read_text()
        
        # Should include next step commands
        assert 'raptor qc' in content
        assert 'raptor profile' in content


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases in report generation."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_empty_stats_dataframe(self, sample_count_matrix, output_dir):
        """Test handling of empty statistics DataFrame."""
        empty_stats = pd.DataFrame()
        output_path = output_dir / 'qc_report.html'
        
        # Should handle gracefully
        generate_html_report(empty_stats, sample_count_matrix, str(output_path), tool='salmon')
        
        assert output_path.exists()
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_small_count_matrix(self, mock_salmon_quant_dir, output_dir):
        """Test report generation with small count matrix."""
        # Create small count matrix
        small_counts = pd.DataFrame({
            'Sample1': [10, 20, 30],
            'Sample2': [15, 25, 35]
        }, index=['Gene1', 'Gene2', 'Gene3'])
        
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        output_path = output_dir / 'qc_report.html'
        
        generate_html_report(stats_df, small_counts, str(output_path), tool='salmon')
        
        assert output_path.exists()


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for report generation workflow."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_full_report_workflow_salmon(self, mock_salmon_quant_dir, sample_count_matrix, output_dir):
        """Test complete report generation workflow for Salmon."""
        # Collect stats
        stats_df = collect_quant_stats(str(mock_salmon_quant_dir), tool='salmon')
        
        # Generate report
        output_path = output_dir / 'qc_report.html'
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='salmon')
        
        # Verify
        assert output_path.exists()
        assert output_path.stat().st_size > 1000  # Should have substantial content
        
        content = output_path.read_text()
        assert 'RAPTOR' in content
        assert 'Module 1' in content or 'Quick-Count' in content
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_full_report_workflow_kallisto(self, mock_kallisto_quant_dir, sample_count_matrix, output_dir):
        """Test complete report generation workflow for Kallisto."""
        stats_df = collect_quant_stats(str(mock_kallisto_quant_dir), tool='kallisto')
        
        output_path = output_dir / 'qc_report.html'
        generate_html_report(stats_df, sample_count_matrix, str(output_path), tool='kallisto')
        
        assert output_path.exists()
        content = output_path.read_text()
        assert 'Kallisto' in content
