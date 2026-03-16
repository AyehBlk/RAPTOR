"""
RAPTOR v2.2.0 - Tests for Data Profiler (Module 3)

Tests for RNAseqDataProfiler and related functions.

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
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import module under test
try:
    from raptor.profiler import (
        RNAseqDataProfiler,
        DataProfile,
        profile_data_quick,
        get_key_characteristics
    )
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False


# =============================================================================
# Tests for RNAseqDataProfiler Initialization
# =============================================================================

class TestProfilerInit:
    """Tests for RNAseqDataProfiler initialization."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init_basic(self, sample_count_matrix):
        """Test basic initialization."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        
        assert profiler.n_samples == sample_count_matrix.shape[1]
        assert profiler.n_genes == sample_count_matrix.shape[0]
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init_with_metadata(self, sample_count_matrix, sample_metadata):
        """Test initialization with metadata."""
        profiler = RNAseqDataProfiler(
            sample_count_matrix, 
            sample_metadata, 
            group_column='condition'
        )
        
        assert profiler.n_groups == 2
        assert profiler.groups is not None
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_init_custom_group_column(self, sample_count_matrix, sample_metadata):
        """Test with custom group column."""
        profiler = RNAseqDataProfiler(
            sample_count_matrix,
            sample_metadata,
            group_column='batch'
        )
        
        assert profiler.n_groups >= 1


# =============================================================================
# Tests for Full Profiling
# =============================================================================

class TestFullProfiling:
    """Tests for run_full_profile method."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_run_full_profile_returns_dataprofile(self, sample_count_matrix):
        """Test that run_full_profile returns DataProfile."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert isinstance(profile, DataProfile)
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_has_all_features(self, sample_count_matrix):
        """Test that profile has 32 features."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert len(profile.feature_vector) == 32
        assert len(profile.feature_names) == 32
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_sample_characteristics(self, sample_count_matrix, sample_metadata):
        """Test sample characteristics extraction."""
        profiler = RNAseqDataProfiler(sample_count_matrix, sample_metadata)
        profile = profiler.run_full_profile()
        
        assert profile.n_samples == sample_count_matrix.shape[1]
        assert profile.n_genes == sample_count_matrix.shape[0]
        assert profile.n_groups == 2
        assert profile.min_group_size >= 1
        assert 0 <= profile.sample_balance <= 1
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_library_size_metrics(self, sample_count_matrix):
        """Test library size metrics extraction."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert profile.library_size_mean > 0
        assert profile.library_size_cv >= 0
        assert profile.library_size_range >= 1
        assert profile.library_size_category in ['low', 'moderate', 'high', 'very_high']
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_dispersion_metrics(self, sample_count_matrix):
        """Test dispersion metrics extraction (CRITICAL)."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert profile.bcv >= 0
        assert profile.bcv_category in ['low', 'moderate', 'high']
        assert profile.common_dispersion >= 0
        assert profile.overdispersion_ratio >= 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_sparsity_metrics(self, sample_count_matrix):
        """Test sparsity metrics extraction."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert 0 <= profile.sparsity <= 1
        assert profile.sparsity_category in ['low', 'moderate', 'high']
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_count_distribution(self, sample_count_matrix):
        """Test count distribution metrics."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        # Proportions should sum to ~1
        total = (profile.zero_proportion + 
                profile.low_count_proportion + 
                profile.medium_count_proportion + 
                profile.high_count_proportion + 
                profile.very_high_count_proportion)
        assert 0.99 <= total <= 1.01
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_mean_variance(self, sample_count_matrix):
        """Test mean-variance relationship metrics."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        # Slope should be > 1 for overdispersed data
        assert profile.mean_var_slope >= 0
        assert 0 <= profile.mean_var_r_squared <= 1


# =============================================================================
# Tests for BCV Categories
# =============================================================================

class TestBCVCategories:
    """Tests for BCV categorization."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_bcv_low_category(self, sample_count_matrix):
        """Test that low BCV is categorized correctly."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        if profile.bcv < 0.2:
            assert profile.bcv_category == 'low'
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_bcv_moderate_category(self, sample_count_matrix):
        """Test that moderate BCV is categorized correctly."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        if 0.2 <= profile.bcv < 0.4:
            assert profile.bcv_category == 'moderate'
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_bcv_high_category(self, sample_count_matrix):
        """Test that high BCV is categorized correctly."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        if profile.bcv >= 0.4:
            assert profile.bcv_category == 'high'


# =============================================================================
# Tests for DataProfile Methods
# =============================================================================

class TestDataProfileMethods:
    """Tests for DataProfile methods."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_to_dict(self, sample_count_matrix):
        """Test to_dict method."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        d = profile.to_dict()
        assert isinstance(d, dict)
        assert 'n_samples' in d
        assert 'bcv' in d
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_to_json(self, sample_count_matrix):
        """Test to_json method."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        j = profile.to_json()
        assert isinstance(j, str)
        
        # Should be valid JSON
        d = json.loads(j)
        assert 'n_samples' in d
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_summary(self, sample_count_matrix):
        """Test summary method."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        s = profile.summary()
        assert isinstance(s, str)
        assert 'RAPTOR' in s
        assert 'BCV' in s
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_get_recommendation_features(self, sample_count_matrix):
        """Test get_recommendation_features method."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        features = profile.get_recommendation_features()
        assert isinstance(features, dict)
        assert 'bcv' in features
        assert 'min_group_size' in features
        assert 'library_size_cv' in features


# =============================================================================
# Tests for Convenience Functions
# =============================================================================

class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_data_quick(self, sample_count_matrix, sample_metadata):
        """Test profile_data_quick function."""
        profile = profile_data_quick(sample_count_matrix, sample_metadata)
        
        assert isinstance(profile, DataProfile)
        assert profile.n_samples > 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_get_key_characteristics(self, sample_count_matrix):
        """Test get_key_characteristics function."""
        chars = get_key_characteristics(sample_count_matrix)
        
        assert isinstance(chars, dict)
        assert 'n_genes' in chars
        assert 'n_samples' in chars
        assert 'bcv' in chars
        assert 'sparsity' in chars


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_small_sample_size(self, small_count_matrix):
        """Test with small sample size."""
        profiler = RNAseqDataProfiler(small_count_matrix)
        profile = profiler.run_full_profile()
        
        assert profile.n_samples == small_count_matrix.shape[1]
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_high_zero_inflation(self, high_zero_count_matrix):
        """Test with high zero inflation."""
        profiler = RNAseqDataProfiler(high_zero_count_matrix)
        profile = profiler.run_full_profile()
        
        assert profile.sparsity > 0.5
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_varying_library_sizes(self, varying_library_size_matrix):
        """Test with varying library sizes."""
        profiler = RNAseqDataProfiler(varying_library_size_matrix)
        profile = profiler.run_full_profile()
        
        assert profile.library_size_cv > 0.1
        assert profile.library_size_range > 2
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_no_metadata(self, sample_count_matrix):
        """Test without metadata."""
        profiler = RNAseqDataProfiler(sample_count_matrix, metadata=None)
        profile = profiler.run_full_profile()
        
        assert profile.n_groups == 1
        assert profile.min_group_size == sample_count_matrix.shape[1]


# =============================================================================
# Tests for Feature Vector
# =============================================================================

class TestFeatureVector:
    """Tests for ML feature vector."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_feature_vector_length(self, sample_count_matrix):
        """Test feature vector has 32 elements."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert len(profile.feature_vector) == 32
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_feature_vector_no_nan(self, sample_count_matrix):
        """Test feature vector has no NaN values."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert not np.isnan(profile.feature_vector).any()
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_feature_vector_no_inf(self, sample_count_matrix):
        """Test feature vector has no infinite values."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert not np.isinf(profile.feature_vector).any()
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_feature_names_match_vector(self, sample_count_matrix):
        """Test feature names match vector length."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        assert len(profile.feature_names) == len(profile.feature_vector)


# =============================================================================
# Tests for Serialization
# =============================================================================

class TestSerialization:
    """Tests for profile serialization."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_json_serialization(self, sample_count_matrix, output_dir):
        """Test JSON serialization."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        # Save
        json_path = output_dir / 'profile.json'
        with open(json_path, 'w') as f:
            f.write(profile.to_json())
        
        # Load
        with open(json_path) as f:
            loaded = json.load(f)
        
        assert loaded['n_samples'] == profile.n_samples
        assert loaded['bcv'] == profile.bcv
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_dict_serialization(self, sample_count_matrix):
        """Test dictionary serialization."""
        profiler = RNAseqDataProfiler(sample_count_matrix)
        profile = profiler.run_full_profile()
        
        d = profile.to_dict()
        
        # All values should be JSON-serializable
        json.dumps(d)  # Should not raise


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for Module 3."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_full_workflow(self, sample_count_matrix, sample_metadata, output_dir):
        """Test complete Module 3 workflow."""
        # Profile
        profiler = RNAseqDataProfiler(sample_count_matrix, sample_metadata)
        profile = profiler.run_full_profile()
        
        # Save
        json_path = output_dir / 'profile.json'
        with open(json_path, 'w') as f:
            json.dump(profile.to_dict(), f, indent=2)
        
        # Verify
        assert json_path.exists()
        assert profile.feature_vector is not None
        assert len(profile.feature_vector) == 32
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_profile_from_file(self, saved_count_matrix, output_dir):
        """Test profiling from file."""
        counts_df = pd.read_csv(saved_count_matrix, index_col=0)
        
        profiler = RNAseqDataProfiler(counts_df)
        profile = profiler.run_full_profile()
        
        assert profile.n_samples > 0
        assert profile.bcv >= 0


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
@pytest.mark.parametrize("min_count_threshold", [1, 5, 10])
def test_min_count_threshold(sample_count_matrix, min_count_threshold):
    """Test different minimum count thresholds."""
    profiler = RNAseqDataProfiler(
        sample_count_matrix,
        min_count_threshold=min_count_threshold
    )
    profile = profiler.run_full_profile()
    
    assert profile.n_samples > 0


@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
@pytest.mark.parametrize("n_samples", [4, 6, 8, 10])
def test_different_sample_sizes(n_samples):
    """Test with different sample sizes."""
    np.random.seed(42)
    counts = pd.DataFrame(
        np.random.negative_binomial(10, 0.1, (100, n_samples)),
        columns=[f'Sample{i+1}' for i in range(n_samples)]
    )
    
    profiler = RNAseqDataProfiler(counts)
    profile = profiler.run_full_profile()
    
    assert profile.n_samples == n_samples
