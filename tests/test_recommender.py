#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Unit Tests - Pipeline Recommender & Threshold Optimizer
======================================================================
Comprehensive pytest tests for PipelineRecommender and ATO classes

Author: Ayeh Bolouki
License: MIT
"""

import pytest
import pandas as pd
import numpy as np

# Import RAPTOR classes
try:
    from raptor import PipelineRecommender
    from raptor.profiler import RNAseqDataProfiler
    RAPTOR_AVAILABLE = True
except ImportError:
    RAPTOR_AVAILABLE = False
    pytest.skip("RAPTOR not installed", allow_module_level=True)


# ============================================================================
# Pipeline Recommender Tests
# ============================================================================

class TestPipelineRecommender:
    """Test suite for PipelineRecommender"""
    
    @pytest.fixture
    def sample_profile_small(self):
        """Create sample profile for small dataset"""
        return {
            'n_samples': 6,
            'n_genes': 2000,
            'bcv': 0.35,
            'bcv_category': 'medium',
            'mean_depth': 20000000,
            'depth_category': 'medium',
            'zero_inflation': 0.45,
            'library_size_cv': 0.15,
            'outliers': [],
            'quality_flags': []
        }
    
    @pytest.fixture
    def sample_profile_large(self):
        """Create sample profile for large dataset"""
        return {
            'n_samples': 48,
            'n_genes': 25000,
            'bcv': 0.65,
            'bcv_category': 'high',
            'mean_depth': 35000000,
            'depth_category': 'high',
            'zero_inflation': 0.38,
            'library_size_cv': 0.22,
            'outliers': [],
            'quality_flags': []
        }
    
    def test_recommender_initialization(self):
        """Test recommender can be initialized"""
        recommender = PipelineRecommender()
        assert recommender is not None
    
    def test_recommend_basic(self, sample_profile_small):
        """Test basic recommendation generation"""
        recommender = PipelineRecommender()
        recommendations = recommender.recommend(sample_profile_small)
        
        assert isinstance(recommendations, list)
        assert len(recommendations) > 0
        assert len(recommendations) <= 8
    
    def test_recommendation_structure(self, sample_profile_small):
        """Test structure of recommendation output"""
        recommender = PipelineRecommender()
        recommendations = recommender.recommend(sample_profile_small, n=3)
        
        assert len(recommendations) == 3
        
        for rec in recommendations:
            assert 'pipeline_id' in rec
            assert 'pipeline_name' in rec
            assert 'score' in rec
            assert 'reasoning' in rec
    
    def test_recommendations_sorted(self, sample_profile_small):
        """Test recommendations are sorted by score"""
        recommender = PipelineRecommender()
        recommendations = recommender.recommend(sample_profile_small, n=5)
        
        scores = [rec['score'] for rec in recommendations]
        assert scores == sorted(scores, reverse=True)
    
    def test_large_dataset_recommendation(self, sample_profile_large):
        """Test recommendation for large dataset"""
        recommender = PipelineRecommender()
        recommendations = recommender.recommend(sample_profile_large, n=1)
        
        top = recommendations[0]
        assert top['pipeline_id'] in [3, 4]  # Fast pipelines


# ============================================================================
# 🆕 Threshold Optimizer Tests (v2.1.1)
# ============================================================================

class TestThresholdOptimizer:
    """Test suite for Adaptive Threshold Optimizer (ATO) - NEW in v2.1.1"""
    
    @pytest.fixture
    def sample_de_results(self):
        """Generate sample DE results for testing"""
        np.random.seed(42)
        n_genes = 1000
        n_de = 150
        
        # Null genes
        null_logfc = np.random.normal(0, 0.2, n_genes - n_de)
        null_pval = np.random.uniform(0.1, 1, n_genes - n_de)
        
        # DE genes
        de_logfc = np.concatenate([
            np.random.normal(2, 0.5, n_de // 2),
            np.random.normal(-2, 0.5, n_de - n_de // 2)
        ])
        de_pval = np.random.exponential(0.001, n_de)
        de_pval = np.clip(de_pval, 1e-300, 0.05)
        
        df = pd.DataFrame({
            'log2FoldChange': np.concatenate([null_logfc, de_logfc]),
            'pvalue': np.concatenate([null_pval, de_pval])
        })
        df.index = [f'Gene_{i}' for i in range(n_genes)]
        
        return df
    
    @pytest.fixture
    def deseq2_like_df(self):
        """Create DataFrame mimicking DESeq2 output."""
        np.random.seed(123)
        n = 500
        
        df = pd.DataFrame({
            'baseMean': np.random.exponential(1000, n),
            'log2FoldChange': np.random.normal(0, 1, n),
            'lfcSE': np.abs(np.random.normal(0.2, 0.05, n)),
            'stat': np.random.normal(0, 2, n),
            'pvalue': np.random.uniform(0, 1, n),
            'padj': np.random.uniform(0, 1, n)
        })
        df.index = [f'ENSG{i:011d}' for i in range(n)]
        
        return df
    
    @pytest.fixture
    def edger_like_df(self):
        """Create DataFrame mimicking edgeR output."""
        np.random.seed(456)
        n = 500
        
        df = pd.DataFrame({
            'logFC': np.random.normal(0, 1, n),
            'logCPM': np.random.normal(5, 2, n),
            'PValue': np.random.uniform(0, 1, n),
            'FDR': np.random.uniform(0, 1, n)
        })
        df.index = [f'Gene{i}' for i in range(n)]
        
        return df
    
    def test_ato_import(self):
        """Test ATO module imports correctly"""
        try:
            from raptor.threshold_optimizer import (
                AdaptiveThresholdOptimizer,
                ThresholdResult,
                optimize_thresholds
            )
            assert True
        except ImportError:
            pytest.skip("ATO module not available")
    
    def test_optimize_thresholds_basic(self, sample_de_results):
        """Test basic threshold optimization"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds
        except ImportError:
            pytest.skip("ATO module not available")
        
        result = optimize_thresholds(sample_de_results, goal='balanced', verbose=False)
        
        assert result.logfc_cutoff > 0
        assert 0 < result.padj_cutoff <= 1
        assert result.n_significant_optimized >= 0
    
    def test_optimize_different_goals(self, sample_de_results):
        """Test different analysis goals"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds
        except ImportError:
            pytest.skip("ATO module not available")
        
        results = {}
        for goal in ['discovery', 'balanced', 'validation']:
            results[goal] = optimize_thresholds(
                sample_de_results, goal=goal, verbose=False
            )
        
        # Discovery should generally be most liberal
        assert results['discovery'].n_significant_optimized >= results['validation'].n_significant_optimized
    
    def test_result_type(self, sample_de_results):
        """Test that result is ThresholdResult type"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds, ThresholdResult
        except ImportError:
            pytest.skip("ATO module not available")
        
        result = optimize_thresholds(sample_de_results, verbose=False)
        assert isinstance(result, ThresholdResult)
    
    def test_summary_method(self, sample_de_results):
        """Test that summary() method works"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds
        except ImportError:
            pytest.skip("ATO module not available")
        
        result = optimize_thresholds(sample_de_results, verbose=False)
        summary = result.summary()
        
        assert summary is not None
        assert len(summary) > 50
        assert 'threshold' in summary.lower() or 'cutoff' in summary.lower()
    
    def test_ato_class_interface(self, sample_de_results):
        """Test AdaptiveThresholdOptimizer class"""
        try:
            from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            pytest.skip("ATO module not available")
        
        ato = AdaptiveThresholdOptimizer(
            sample_de_results, 
            goal='balanced', 
            verbose=False
        )
        result = ato.optimize()
        
        assert hasattr(result, 'logfc_cutoff')
        assert hasattr(result, 'padj_cutoff')
        assert hasattr(result, 'padj_method')
        assert hasattr(result, 'logfc_method')
        
        sig_genes = ato.get_significant_genes()
        assert isinstance(sig_genes, pd.DataFrame)
    
    def test_column_standardization_deseq2(self, deseq2_like_df):
        """Test column name standardization for DESeq2"""
        try:
            from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            pytest.skip("ATO module not available")
        
        ato = AdaptiveThresholdOptimizer(deseq2_like_df, verbose=False)
        assert 'logfc' in ato.df.columns
        assert 'pvalue' in ato.df.columns
    
    def test_column_standardization_edger(self, edger_like_df):
        """Test column name standardization for edgeR"""
        try:
            from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            pytest.skip("ATO module not available")
        
        ato = AdaptiveThresholdOptimizer(edger_like_df, verbose=False)
        assert 'logfc' in ato.df.columns
        assert 'pvalue' in ato.df.columns
    
    def test_invalid_goal(self, sample_de_results):
        """Test invalid goal raises error"""
        try:
            from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            pytest.skip("ATO module not available")
        
        with pytest.raises(ValueError):
            AdaptiveThresholdOptimizer(sample_de_results, goal='invalid', verbose=False)
    
    def test_missing_columns(self):
        """Test error on missing required columns"""
        try:
            from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            pytest.skip("ATO module not available")
        
        df = pd.DataFrame({'somecolumn': [1, 2, 3]})
        with pytest.raises(ValueError):
            AdaptiveThresholdOptimizer(df, verbose=False)
    
    def test_small_dataset(self):
        """Test with very small dataset"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds
        except ImportError:
            pytest.skip("ATO module not available")
        
        df = pd.DataFrame({
            'log2FoldChange': [0.5, -0.5, 2.0],
            'pvalue': [0.1, 0.2, 0.001]
        })
        
        result = optimize_thresholds(df, verbose=False)
        assert result is not None
    
    def test_with_nan_values(self):
        """Test handling of NaN values"""
        try:
            from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            pytest.skip("ATO module not available")
        
        df = pd.DataFrame({
            'log2FoldChange': [0.5, np.nan, 2.0, -1.0],
            'pvalue': [0.1, 0.2, np.nan, 0.001]
        })
        
        ato = AdaptiveThresholdOptimizer(df, verbose=False)
        assert len(ato.df) < 4  # Should have removed NaN rows
    
    def test_result_attributes(self, sample_de_results):
        """Test all result attributes are present"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds
        except ImportError:
            pytest.skip("ATO module not available")
        
        result = optimize_thresholds(sample_de_results, verbose=False)
        
        # Check all required attributes
        assert hasattr(result, 'logfc_cutoff')
        assert hasattr(result, 'padj_cutoff')
        assert hasattr(result, 'padj_method')
        assert hasattr(result, 'logfc_method')
        assert hasattr(result, 'logfc_reasoning')
        assert hasattr(result, 'padj_reasoning')
        assert hasattr(result, 'n_significant_optimized')
        assert hasattr(result, 'n_significant_traditional')
        assert hasattr(result, 'pi0_estimate')
    
    def test_result_values_reasonable(self, sample_de_results):
        """Test result values are reasonable"""
        try:
            from raptor.threshold_optimizer import optimize_thresholds
        except ImportError:
            pytest.skip("ATO module not available")
        
        result = optimize_thresholds(sample_de_results, verbose=False)
        
        assert 0 < result.logfc_cutoff < 10
        assert 0 < result.padj_cutoff <= 0.1
        assert 0 < result.pi0_estimate <= 1
        assert result.n_significant_optimized >= 0
        assert result.n_significant_traditional >= 0


# Pytest configuration
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
