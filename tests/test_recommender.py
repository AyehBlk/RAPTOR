#!/usr/bin/env python3
"""
RAPTOR Unit Tests - Pipeline Recommender
=========================================
Comprehensive tests for PipelineRecommender class

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
except ImportError:
    pytest.skip("RAPTOR not installed", allow_module_level=True)


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


# Pytest configuration
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
