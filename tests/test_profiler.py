#!/usr/bin/env python3
"""
RAPTOR Unit Tests - Data Profiler
==================================
Comprehensive tests for RNAseqDataProfiler class

Author: Ayeh Bolouki
License: MIT
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os

# Import RAPTOR classes
try:
    from raptor import RNAseqDataProfiler
    from raptor.utils import validate_count_matrix, validate_metadata
except ImportError:
    pytest.skip("RAPTOR not installed", allow_module_level=True)


class TestRNAseqDataProfiler:
    """Test suite for RNAseqDataProfiler"""
    
    @pytest.fixture
    def simple_counts(self):
        """Create simple count matrix for testing"""
        np.random.seed(42)
        genes = [f"Gene_{i}" for i in range(100)]
        samples = [f"Sample_{i}" for i in range(6)]
        
        # Generate realistic count data
        data = np.random.negative_binomial(n=10, p=0.3, size=(100, 6))
        counts = pd.DataFrame(data, index=genes, columns=samples)
        
        return counts
    
    @pytest.fixture
    def simple_metadata(self):
        """Create simple metadata for testing"""
        metadata = pd.DataFrame({
            'sample': [f"Sample_{i}" for i in range(6)],
            'condition': ['Control', 'Control', 'Control', 
                         'Treatment', 'Treatment', 'Treatment'],
            'replicate': [1, 2, 3, 1, 2, 3]
        })
        return metadata
    
    @pytest.fixture
    def complex_counts(self):
        """Create more complex count matrix"""
        np.random.seed(42)
        genes = [f"Gene_{i}" for i in range(5000)]
        samples = [f"Sample_{i}" for i in range(12)]
        
        # Simulate different expression patterns
        # Group 1: Low expression
        low_expr = np.random.negative_binomial(n=2, p=0.5, size=(2000, 12))
        # Group 2: Medium expression
        med_expr = np.random.negative_binomial(n=10, p=0.3, size=(2000, 12))
        # Group 3: High expression
        high_expr = np.random.negative_binomial(n=50, p=0.1, size=(1000, 12))
        
        data = np.vstack([low_expr, med_expr, high_expr])
        counts = pd.DataFrame(data, index=genes, columns=samples)
        
        return counts
    
    def test_profiler_initialization(self, simple_counts, simple_metadata):
        """Test profiler can be initialized"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        assert profiler is not None
        assert profiler.counts.shape == (100, 6)
        assert profiler.metadata.shape[0] == 6
    
    def test_profiler_without_metadata(self, simple_counts):
        """Test profiler works without metadata"""
        profiler = RNAseqDataProfiler(simple_counts, metadata=None)
        assert profiler is not None
        # Should create default metadata
        assert profiler.metadata is not None or profiler._has_metadata == False
    
    def test_profile_execution(self, simple_counts, simple_metadata):
        """Test profile() runs successfully"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert isinstance(profile, dict)
        assert 'n_samples' in profile
        assert 'n_genes' in profile
        assert 'bcv' in profile
        assert 'mean_depth' in profile
    
    def test_profile_sample_count(self, simple_counts, simple_metadata):
        """Test correct sample count"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert profile['n_samples'] == 6
    
    def test_profile_gene_count(self, simple_counts, simple_metadata):
        """Test correct gene count"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert profile['n_genes'] == 100
    
    def test_bcv_calculation(self, simple_counts, simple_metadata):
        """Test BCV is calculated correctly"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'bcv' in profile
        assert isinstance(profile['bcv'], (int, float))
        assert profile['bcv'] >= 0
        assert profile['bcv'] <= 2.0  # Reasonable range
    
    def test_bcv_categories(self, simple_counts, simple_metadata):
        """Test BCV categorization"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'bcv_category' in profile
        assert profile['bcv_category'] in ['low', 'medium', 'high', 'very_high']
    
    def test_depth_calculation(self, simple_counts, simple_metadata):
        """Test sequencing depth calculation"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'mean_depth' in profile
        assert profile['mean_depth'] > 0
        
        # Check it matches library sizes
        expected_mean = simple_counts.sum(axis=0).mean()
        assert abs(profile['mean_depth'] - expected_mean) < 1.0
    
    def test_depth_categories(self, simple_counts, simple_metadata):
        """Test depth categorization"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'depth_category' in profile
        assert profile['depth_category'] in ['very_low', 'low', 'medium', 'high', 'very_high']
    
    def test_zero_inflation(self, simple_counts, simple_metadata):
        """Test zero inflation calculation"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'zero_inflation' in profile
        assert 0 <= profile['zero_inflation'] <= 1.0
        
        # Check calculation
        expected_zeros = (simple_counts == 0).sum().sum() / simple_counts.size
        assert abs(profile['zero_inflation'] - expected_zeros) < 0.01
    
    def test_library_size_cv(self, simple_counts, simple_metadata):
        """Test library size coefficient of variation"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'library_size_cv' in profile
        assert profile['library_size_cv'] >= 0
        
        # Check calculation
        lib_sizes = simple_counts.sum(axis=0)
        expected_cv = lib_sizes.std() / lib_sizes.mean()
        assert abs(profile['library_size_cv'] - expected_cv) < 0.01
    
    def test_outlier_detection(self, simple_counts, simple_metadata):
        """Test outlier detection"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'outliers' in profile
        assert isinstance(profile['outliers'], list)
    
    def test_quality_flags(self, simple_counts, simple_metadata):
        """Test quality flag generation"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        assert 'quality_flags' in profile
        assert isinstance(profile['quality_flags'], list)
    
    def test_with_complex_data(self, complex_counts):
        """Test with more realistic complex dataset"""
        profiler = RNAseqDataProfiler(complex_counts, metadata=None)
        profile = profiler.profile()
        
        assert profile['n_genes'] == 5000
        assert profile['n_samples'] == 12
        assert profile['bcv'] > 0
        assert profile['mean_depth'] > 0
    
    def test_invalid_count_matrix_negative_values(self):
        """Test that negative values raise error"""
        counts = pd.DataFrame({
            'Sample1': [10, 20, -5],  # Negative value
            'Sample2': [15, 25, 10]
        }, index=['Gene1', 'Gene2', 'Gene3'])
        
        with pytest.raises((ValueError, AssertionError)):
            validate_count_matrix(counts)
    
    def test_invalid_count_matrix_missing_values(self):
        """Test that missing values raise error"""
        counts = pd.DataFrame({
            'Sample1': [10, np.nan, 30],  # Missing value
            'Sample2': [15, 25, 10]
        }, index=['Gene1', 'Gene2', 'Gene3'])
        
        with pytest.raises((ValueError, AssertionError)):
            validate_count_matrix(counts)
    
    def test_metadata_sample_mismatch(self, simple_counts):
        """Test detection of metadata-count mismatch"""
        wrong_metadata = pd.DataFrame({
            'sample': ['Wrong1', 'Wrong2', 'Wrong3'],
            'condition': ['A', 'B', 'C']
        })
        
        # Should raise warning or error
        with pytest.raises((ValueError, KeyError)):
            validate_metadata(wrong_metadata, simple_counts)
    
    def test_save_profile_results(self, simple_counts, simple_metadata, tmp_path):
        """Test saving profile results to file"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        profile = profiler.profile()
        
        # Save to temp file
        output_file = tmp_path / "profile_results.json"
        
        # Assuming profiler has a save method
        # profiler.save_profile(profile, output_file)
        
        import json
        with open(output_file, 'w') as f:
            json.dump(profile, f)
        
        assert output_file.exists()
        
        # Load and verify
        with open(output_file) as f:
            loaded = json.load(f)
        
        assert loaded['n_samples'] == profile['n_samples']
        assert loaded['n_genes'] == profile['n_genes']
    
    def test_high_zero_inflation_warning(self):
        """Test warning for high zero inflation"""
        # Create data with many zeros
        counts = pd.DataFrame(
            np.random.choice([0, 1, 2], size=(100, 6), p=[0.8, 0.15, 0.05]),
            index=[f"Gene_{i}" for i in range(100)],
            columns=[f"Sample_{i}" for i in range(6)]
        )
        
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        assert profile['zero_inflation'] > 0.5
        # Should generate quality flag
        assert len(profile['quality_flags']) > 0
    
    def test_unbalanced_library_sizes(self):
        """Test detection of highly variable library sizes"""
        np.random.seed(42)
        
        # Create unbalanced library sizes
        counts = pd.DataFrame({
            'Sample1': np.random.poisson(100, 100),    # Low
            'Sample2': np.random.poisson(100, 100),    # Low
            'Sample3': np.random.poisson(1000, 100),   # High
            'Sample4': np.random.poisson(1000, 100),   # High
        }, index=[f"Gene_{i}" for i in range(100)])
        
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        # Should have high CV
        assert profile['library_size_cv'] > 0.3
        # Should generate quality flag
        assert any('library' in flag.lower() for flag in profile['quality_flags'])
    
    def test_small_sample_size(self):
        """Test handling of very small sample sizes"""
        counts = pd.DataFrame({
            'Sample1': [100, 200, 150],
            'Sample2': [110, 210, 160],
        }, index=['Gene1', 'Gene2', 'Gene3'])
        
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        assert profile['n_samples'] == 2
        # Should generate warning
        assert any('sample' in flag.lower() for flag in profile['quality_flags'])


class TestProfilerEdgeCases:
    """Test edge cases and error handling"""
    
    def test_single_sample(self):
        """Test with single sample (should fail or warn)"""
        counts = pd.DataFrame({
            'Sample1': [100, 200, 150]
        }, index=['Gene1', 'Gene2', 'Gene3'])
        
        # Should either raise error or give strong warning
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        assert len(profile['quality_flags']) > 0
    
    def test_single_gene(self):
        """Test with single gene"""
        counts = pd.DataFrame({
            'Sample1': [100],
            'Sample2': [110],
            'Sample3': [120],
        }, index=['Gene1'])
        
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        assert profile['n_genes'] == 1
    
    def test_all_zeros(self):
        """Test with all zero counts"""
        counts = pd.DataFrame(
            np.zeros((10, 6)),
            index=[f"Gene_{i}" for i in range(10)],
            columns=[f"Sample_{i}" for i in range(6)]
        )
        
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        assert profile['zero_inflation'] == 1.0
        assert profile['mean_depth'] == 0
    
    def test_very_large_counts(self):
        """Test with very large count values"""
        counts = pd.DataFrame(
            np.random.poisson(1e6, size=(100, 6)),
            index=[f"Gene_{i}" for i in range(100)],
            columns=[f"Sample_{i}" for i in range(6)]
        )
        
        profiler = RNAseqDataProfiler(counts, metadata=None)
        profile = profiler.profile()
        
        assert profile['mean_depth'] > 1e7
        assert profile['depth_category'] == 'very_high'


class TestProfilerIntegration:
    """Integration tests with file I/O"""
    
    def test_load_from_csv(self, tmp_path, simple_counts, simple_metadata):
        """Test loading data from CSV files"""
        # Save to temp files
        counts_file = tmp_path / "counts.csv"
        metadata_file = tmp_path / "metadata.csv"
        
        simple_counts.to_csv(counts_file)
        simple_metadata.to_csv(metadata_file, index=False)
        
        # Load and profile
        counts = pd.read_csv(counts_file, index_col=0)
        metadata = pd.read_csv(metadata_file)
        
        profiler = RNAseqDataProfiler(counts, metadata)
        profile = profiler.profile()
        
        assert profile['n_samples'] == 6
        assert profile['n_genes'] == 100
    
    def test_profile_with_config(self, simple_counts, simple_metadata, tmp_path):
        """Test profiling with custom configuration"""
        # Create custom config
        config_file = tmp_path / "config.yaml"
        
        import yaml
        config = {
            'profiling': {
                'bcv_thresholds': {
                    'low': 0.15,
                    'medium': 0.35,
                    'high': 0.60
                }
            }
        }
        
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        
        # Profile with config
        profiler = RNAseqDataProfiler(
            simple_counts, 
            simple_metadata,
            config_file=str(config_file)
        )
        profile = profiler.profile()
        
        assert 'bcv' in profile
        assert 'bcv_category' in profile


# Pytest configuration
def pytest_configure(config):
    """Configure pytest"""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "--tb=short"])
