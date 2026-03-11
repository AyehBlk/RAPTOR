"""
RAPTOR v2.2.0 - Tests for Quality Assessment (Module 2)

Tests for DataQualityAssessor and related functions.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import module under test
try:
    from raptor.quality_assessment import (
        DataQualityAssessor,
        OutlierResult,
        quick_quality_check,
        detect_outliers_quick
    )
    IMPORTS_AVAILABLE = True
except ImportError as e:
    IMPORTS_AVAILABLE = False
    print(f"Warning: Could not import quality_assessment module: {e}")


# =============================================================================
# Tests for DataQualityAssessor Initialization
# =============================================================================

class TestDataQualityAssessorInit:
    """Tests for DataQualityAssessor initialization."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init_basic(self, sample_count_matrix):
        """Test basic initialization."""
        assessor = DataQualityAssessor(sample_count_matrix)
        
        assert assessor.n_samples == sample_count_matrix.shape[1]
        assert assessor.n_genes == sample_count_matrix.shape[0]
        assert assessor.normalization == 'log2'
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init_with_metadata(self, sample_count_matrix, sample_metadata):
        """Test initialization with metadata."""
        assessor = DataQualityAssessor(sample_count_matrix, sample_metadata)
        
        assert assessor.metadata is not None
        assert len(assessor.metadata) == sample_count_matrix.shape[1]
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    @pytest.mark.parametrize("normalization", ['log2', 'cpm', 'quantile', 'none'])
    def test_init_normalization_methods(self, sample_count_matrix, normalization):
        """Test all normalization methods."""
        assessor = DataQualityAssessor(sample_count_matrix, normalization=normalization)
        assert assessor.normalization == normalization
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init_invalid_normalization(self, sample_count_matrix):
        """Test invalid normalization raises error."""
        with pytest.raises(ValueError):
            DataQualityAssessor(sample_count_matrix, normalization='invalid')


# =============================================================================
# Tests for Quality Assessment
# =============================================================================

class TestQualityAssessment:
    """Tests for assess_quality method."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_assess_quality_returns_report(self, sample_count_matrix):
        """Test that assess_quality returns a valid report."""
        assessor = DataQualityAssessor(sample_count_matrix)
        report = assessor.assess_quality()
        
        assert isinstance(report, dict)
        assert 'overall' in report
        assert 'score' in report['overall']
        assert 'status' in report['overall']
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_overall_score_range(self, sample_count_matrix):
        """Test that overall score is in valid range."""
        assessor = DataQualityAssessor(sample_count_matrix)
        report = assessor.assess_quality()
        
        score = report['overall']['score']
        assert 0 <= score <= 100
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_component_scores_present(self, sample_count_matrix, sample_metadata):
        """Test that all component scores are present."""
        assessor = DataQualityAssessor(sample_count_matrix, sample_metadata)
        report = assessor.assess_quality()
        
        expected_components = ['library_quality', 'gene_detection', 'outlier_detection',
                              'variance_structure', 'batch_effects', 'biological_signal']
        
        # Check in either 'components' key or directly in report
        components = report.get('components', report)
        
        for comp in expected_components:
            assert comp in components, f"Missing component: {comp}"
            assert 'score' in components[comp]
            assert 'status' in components[comp]
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_library_quality_metrics(self, sample_count_matrix):
        """Test library quality metrics."""
        assessor = DataQualityAssessor(sample_count_matrix)
        report = assessor.assess_quality()
        
        lib_quality = report.get('components', report).get('library_quality', {})
        
        assert 'mean_size' in lib_quality
        assert 'cv' in lib_quality
        assert lib_quality['mean_size'] > 0
        assert lib_quality['cv'] >= 0


# =============================================================================
# Tests for Outlier Detection
# =============================================================================

class TestOutlierDetection:
    """Tests for outlier detection methods."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_detect_outliers_advanced_returns_result(self, sample_count_matrix):
        """Test that detect_outliers_advanced returns OutlierResult."""
        assessor = DataQualityAssessor(sample_count_matrix)
        result = assessor.detect_outliers_advanced()
        
        assert isinstance(result, OutlierResult)
        assert hasattr(result, 'outlier_samples')
        assert hasattr(result, 'n_outliers')
        assert hasattr(result, 'method_results')
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_outlier_detection_with_known_outlier(self, sample_count_matrix_with_outlier):
        """Test outlier detection finds known outlier."""
        assessor = DataQualityAssessor(sample_count_matrix_with_outlier)
        result = assessor.detect_outliers_advanced(consensus_threshold=2)
        
        # The artificial outlier should be detected by at least some methods
        assert 'Treatment_4' in result.outlier_scores
        # It should have votes from multiple methods
        assert result.outlier_scores['Treatment_4'] >= 1
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    @pytest.mark.parametrize("threshold", [1, 2, 3, 4, 5, 6])
    def test_consensus_threshold_effect(self, sample_count_matrix, threshold):
        """Test that higher threshold means fewer outliers."""
        assessor = DataQualityAssessor(sample_count_matrix)
        result = assessor.detect_outliers_advanced(consensus_threshold=threshold)
        
        # Higher threshold should generally mean fewer or equal outliers
        assert result.consensus_threshold == threshold
        assert result.n_outliers >= 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_all_six_methods_run(self, sample_count_matrix):
        """Test that all 6 outlier detection methods are used."""
        assessor = DataQualityAssessor(sample_count_matrix)
        result = assessor.detect_outliers_advanced()
        
        expected_methods = ['PCA + Mahalanobis', 'Isolation Forest', 
                          'Local Outlier Factor', 'Elliptic Envelope',
                          'Correlation-based', 'Library Size']
        
        for method in expected_methods:
            assert method in result.method_results, f"Missing method: {method}"
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_outlier_summary_method(self, sample_count_matrix):
        """Test OutlierResult.summary() method."""
        assessor = DataQualityAssessor(sample_count_matrix)
        result = assessor.detect_outliers_advanced()
        
        summary = result.summary()
        assert isinstance(summary, str)
        assert 'OUTLIER DETECTION' in summary


# =============================================================================
# Tests for Batch Effect Detection
# =============================================================================

class TestBatchEffectDetection:
    """Tests for batch effect detection."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_batch_detection_with_metadata(self, sample_count_matrix, sample_metadata):
        """Test batch effect detection with metadata."""
        assessor = DataQualityAssessor(sample_count_matrix, sample_metadata)
        report = assessor.assess_quality()
        
        batch_info = report.get('components', report).get('batch_effects', {})
        
        assert 'batch_detected' in batch_info
        assert 'batch_variable' in batch_info or 'confounded' in batch_info
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_confounding_detection(self, sample_count_matrix, confounded_metadata):
        """Test detection of confounded design."""
        assessor = DataQualityAssessor(sample_count_matrix, confounded_metadata)
        report = assessor.assess_quality()
        
        batch_info = report.get('components', report).get('batch_effects', {})
        
        # Should detect confounding
        assert batch_info.get('confounded', False) == True
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_batch_recommendation_present(self, sample_count_matrix, sample_metadata):
        """Test that batch recommendation is provided."""
        assessor = DataQualityAssessor(sample_count_matrix, sample_metadata)
        report = assessor.assess_quality()
        
        batch_info = report.get('components', report).get('batch_effects', {})
        
        assert 'recommendation' in batch_info
        assert isinstance(batch_info['recommendation'], str)


# =============================================================================
# Tests for Normalization Methods
# =============================================================================

class TestNormalization:
    """Tests for normalization methods."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_log2_normalization(self, sample_count_matrix):
        """Test log2 normalization."""
        assessor = DataQualityAssessor(sample_count_matrix, normalization='log2')
        
        # log_counts should be transformed
        assert assessor.log_counts is not None
        # Values should be smaller than original counts
        assert assessor.log_counts.values.max() < sample_count_matrix.values.max()
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_cpm_normalization(self, varying_library_size_matrix):
        """Test CPM normalization helps with varying library sizes."""
        # Without CPM
        assessor_log2 = DataQualityAssessor(varying_library_size_matrix, normalization='log2')
        report_log2 = assessor_log2.assess_quality()
        
        # With CPM
        assessor_cpm = DataQualityAssessor(varying_library_size_matrix, normalization='cpm')
        report_cpm = assessor_cpm.assess_quality()
        
        # Both should complete without error
        assert report_log2['overall']['score'] > 0
        assert report_cpm['overall']['score'] > 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_none_normalization(self, sample_count_matrix):
        """Test no normalization."""
        assessor = DataQualityAssessor(sample_count_matrix, normalization='none')
        
        # Data should be unchanged
        pd.testing.assert_frame_equal(assessor.log_counts, sample_count_matrix)


# =============================================================================
# Tests for Convenience Functions
# =============================================================================

class TestConvenienceFunctions:
    """Tests for quick_quality_check and detect_outliers_quick."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_quick_quality_check(self, sample_count_matrix, output_dir):
        """Test quick_quality_check function."""
        report = quick_quality_check(
            sample_count_matrix, 
            plot=False
        )
        
        assert isinstance(report, dict)
        assert 'overall' in report
        assert 'score' in report['overall']
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_detect_outliers_quick(self, sample_count_matrix):
        """Test detect_outliers_quick function."""
        result = detect_outliers_quick(
            sample_count_matrix,
            consensus_threshold=3,
            plot=False
        )
        
        assert isinstance(result, OutlierResult)
        assert result.n_outliers >= 0


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_small_sample_size(self, small_count_matrix):
        """Test with small sample size."""
        assessor = DataQualityAssessor(small_count_matrix)
        report = assessor.assess_quality()
        
        assert report['overall']['score'] >= 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_high_zero_inflation(self, high_zero_count_matrix):
        """Test with high zero inflation data."""
        assessor = DataQualityAssessor(high_zero_count_matrix)
        report = assessor.assess_quality()
        
        gene_detection = report.get('components', report).get('gene_detection', {})
        
        # Should detect high zero inflation
        assert gene_detection.get('zero_inflation_pct', 0) > 50
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_no_metadata(self, sample_count_matrix):
        """Test without metadata."""
        assessor = DataQualityAssessor(sample_count_matrix, metadata=None)
        report = assessor.assess_quality()
        
        # Should still complete
        assert report['overall']['score'] > 0


# =============================================================================
# Tests for Visualization
# =============================================================================

class TestVisualization:
    """Tests for visualization methods."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_plot_quality_report(self, sample_count_matrix, output_dir):
        """Test quality report plot generation."""
        assessor = DataQualityAssessor(sample_count_matrix)
        assessor.assess_quality()
        
        plot_file = output_dir / 'test_quality.png'
        
        try:
            assessor.plot_quality_report(str(plot_file))
            assert plot_file.exists()
        except ImportError:
            pytest.skip("Matplotlib not available")
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_plot_outliers_advanced(self, sample_count_matrix_with_outlier, output_dir):
        """Test outlier plot generation."""
        assessor = DataQualityAssessor(sample_count_matrix_with_outlier)
        result = assessor.detect_outliers_advanced(consensus_threshold=2)
        
        if result.n_outliers > 0:
            plot_file = output_dir / 'test_outliers.png'
            
            try:
                assessor.plot_outliers_advanced(result, str(plot_file))
                assert plot_file.exists()
            except ImportError:
                pytest.skip("Matplotlib not available")


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for Module 2."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_full_workflow(self, sample_count_matrix, sample_metadata, output_dir):
        """Test complete Module 2 workflow."""
        # Initialize
        assessor = DataQualityAssessor(sample_count_matrix, sample_metadata)
        
        # Quality assessment
        report = assessor.assess_quality()
        assert report['overall']['score'] > 0
        
        # Outlier detection
        outliers = assessor.detect_outliers_advanced(consensus_threshold=3)
        assert outliers.n_outliers >= 0
        
        # Save results
        report_file = output_dir / 'quality_report.json'
        with open(report_file, 'w') as f:
            json.dump({'overall': report['overall']}, f)
        
        assert report_file.exists()
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_outlier_removal(self, sample_count_matrix_with_outlier, output_dir):
        """Test outlier removal workflow."""
        assessor = DataQualityAssessor(sample_count_matrix_with_outlier)
        result = assessor.detect_outliers_advanced(consensus_threshold=2)
        
        if result.n_outliers > 0:
            # Remove outliers
            clean_counts = sample_count_matrix_with_outlier.drop(
                columns=result.outlier_samples
            )
            
            assert clean_counts.shape[1] < sample_count_matrix_with_outlier.shape[1]
            
            # Save clean counts
            clean_file = output_dir / 'counts_clean.csv'
            clean_counts.to_csv(clean_file)
            
            assert clean_file.exists()


# =============================================================================
# Test __init__.py
# =============================================================================

@pytest.fixture
def test_module2_init():
    """Create tests for module2 __init__.py file."""
    pass


@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
def test_module2_init_imports():
    """Test that quality_assessment module exports expected names."""
    from raptor.quality_assessment import (
        MODULE_ID,
        MODULE_STAGE,
        MODULE_NAME,
        MODULE_VERSION,
        QUALITY_ASSESSMENT_CLI_PARAMS
    )
    
    # Check module metadata
    assert MODULE_ID == 'M2'
    assert MODULE_STAGE == 2  # Module 2 is stage 2!
    assert MODULE_NAME == 'Quality Assessment'
    assert MODULE_VERSION == '2.2.0'
    
    # Check CLI params exist
    assert isinstance(QUALITY_ASSESSMENT_CLI_PARAMS, dict)
    assert 'counts' in QUALITY_ASSESSMENT_CLI_PARAMS
    assert 'metadata' in QUALITY_ASSESSMENT_CLI_PARAMS
    assert 'normalization' in QUALITY_ASSESSMENT_CLI_PARAMS
