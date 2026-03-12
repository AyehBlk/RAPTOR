"""
RAPTOR v2.2.0 - Tests for Pipeline Recommender (Module 4)

Tests for PipelineRecommender, MLPipelineRecommender, and related functions.

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
    from raptor.recommender import (
        PipelineRecommender,
        Recommendation,
        PipelineInfo,
        recommend_pipeline,
        get_pipeline_info,
        list_pipelines,
        PIPELINES
    )
    RECOMMENDER_AVAILABLE = True
except ImportError:
    RECOMMENDER_AVAILABLE = False

try:
    from raptor.ml_recommender import MLPipelineRecommender, MLRecommendation
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False

try:
    from raptor.profiler import RNAseqDataProfiler, DataProfile
    PROFILER_AVAILABLE = True
except ImportError:
    PROFILER_AVAILABLE = False

try:
    from raptor.simulation import simulate_rnaseq, SimulationResult
    SIMULATION_AVAILABLE = True
except ImportError:
    SIMULATION_AVAILABLE = False


# =============================================================================
# Fixtures for Module 4
# =============================================================================

@pytest.fixture
def demo_profile():
    """Create a demo profile for testing."""
    class DemoProfile:
        def __init__(self):
            self.n_samples = 8
            self.n_genes = 15000
            self.n_groups = 2
            self.min_group_size = 4
            self.max_group_size = 4
            self.sample_balance = 1.0
            self.has_replicates = True
            self.has_sufficient_replicates = True
            self.design_complexity = 'simple'
            
            self.bcv = 0.32
            self.bcv_category = 'moderate'
            self.common_dispersion = 0.10
            self.overdispersion_ratio = 2.5
            
            self.library_size_cv = 0.15
            self.library_size_range = 1.8
            
            self.sparsity = 0.45
            self.low_count_proportion = 0.12
            
            self.has_outliers = False
            self.outlier_severity = 'none'
            self.has_batch_effect = False
            self.batch_confounded = False
            self.quality_score = 78.0
    
    return DemoProfile()


@pytest.fixture
def small_sample_profile():
    """Profile with small sample size."""
    class SmallProfile:
        def __init__(self):
            self.n_samples = 4
            self.n_genes = 10000
            self.n_groups = 2
            self.min_group_size = 2
            self.max_group_size = 2
            self.sample_balance = 1.0
            self.has_replicates = True
            self.has_sufficient_replicates = False
            self.design_complexity = 'simple'
            
            self.bcv = 0.40
            self.bcv_category = 'high'
            self.common_dispersion = 0.16
            self.overdispersion_ratio = 3.0
            
            self.library_size_cv = 0.20
            self.library_size_range = 2.0
            
            self.sparsity = 0.50
            self.low_count_proportion = 0.20
            
            self.has_outliers = False
            self.outlier_severity = 'none'
            self.has_batch_effect = False
            self.batch_confounded = False
            self.quality_score = 65.0
    
    return SmallProfile()


@pytest.fixture
def large_sample_profile():
    """Profile with large sample size (n=25 per group)."""
    class LargeProfile:
        def __init__(self):
            self.n_samples = 50
            self.n_genes = 20000
            self.n_groups = 2
            self.min_group_size = 25
            self.max_group_size = 25
            self.sample_balance = 1.0
            self.has_replicates = True
            self.has_sufficient_replicates = True
            self.design_complexity = 'simple'
            
            self.bcv = 0.18
            self.bcv_category = 'low'
            self.common_dispersion = 0.03
            self.overdispersion_ratio = 1.5
            
            self.library_size_cv = 0.12
            self.library_size_range = 1.5
            
            self.sparsity = 0.35
            self.low_count_proportion = 0.08
            
            self.has_outliers = False
            self.outlier_severity = 'none'
            self.has_batch_effect = False
            self.batch_confounded = False
            self.quality_score = 90.0
    
    return LargeProfile()


@pytest.fixture
def outlier_profile():
    """Profile with outliers."""
    class OutlierProfile:
        def __init__(self):
            self.n_samples = 8
            self.n_genes = 15000
            self.n_groups = 2
            self.min_group_size = 4
            self.max_group_size = 4
            self.sample_balance = 1.0
            self.has_replicates = True
            self.has_sufficient_replicates = True
            self.design_complexity = 'simple'
            
            self.bcv = 0.35
            self.bcv_category = 'moderate'
            self.common_dispersion = 0.12
            self.overdispersion_ratio = 2.8
            
            self.library_size_cv = 0.25
            self.library_size_range = 3.0
            
            self.sparsity = 0.48
            self.low_count_proportion = 0.15
            
            self.has_outliers = True
            self.outlier_severity = 'moderate'
            self.has_batch_effect = False
            self.batch_confounded = False
            self.quality_score = 60.0
    
    return OutlierProfile()


@pytest.fixture
def batch_effect_profile():
    """Profile with batch effects."""
    class BatchProfile:
        def __init__(self):
            self.n_samples = 12
            self.n_genes = 18000
            self.n_groups = 2
            self.min_group_size = 6
            self.max_group_size = 6
            self.sample_balance = 1.0
            self.has_replicates = True
            self.has_sufficient_replicates = True
            self.design_complexity = 'moderate'
            
            self.bcv = 0.28
            self.bcv_category = 'moderate'
            self.common_dispersion = 0.08
            self.overdispersion_ratio = 2.2
            
            self.library_size_cv = 0.18
            self.library_size_range = 2.2
            
            self.sparsity = 0.42
            self.low_count_proportion = 0.11
            
            self.has_outliers = False
            self.outlier_severity = 'none'
            self.has_batch_effect = True
            self.batch_confounded = False
            self.quality_score = 72.0
    
    return BatchProfile()


# =============================================================================
# Tests for PipelineRecommender
# =============================================================================

class TestPipelineRecommenderInit:
    """Tests for PipelineRecommender initialization."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_init_basic(self, demo_profile):
        """Test basic initialization."""
        recommender = PipelineRecommender(demo_profile)
        assert recommender.profile is not None
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_init_creates_empty_scores(self, demo_profile):
        """Test that initialization creates empty scores."""
        recommender = PipelineRecommender(demo_profile)
        assert hasattr(recommender, 'scores')


class TestGetRecommendation:
    """Tests for get_recommendation method."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_returns_recommendation(self, demo_profile):
        """Test that get_recommendation returns Recommendation object."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        assert isinstance(rec, Recommendation)
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_recommendation_has_primary(self, demo_profile):
        """Test that recommendation has primary pipeline."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        assert rec.primary_pipeline in PIPELINES
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_recommendation_has_scores(self, demo_profile):
        """Test that recommendation has scores."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        assert 0 <= rec.primary_score <= 100
        assert 0 <= rec.alternative_score <= 100
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_recommendation_has_reason(self, demo_profile):
        """Test that recommendation has reason."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        assert len(rec.primary_reason) > 0
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_recommendation_has_r_code(self, demo_profile):
        """Test that recommendation includes R code."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        assert len(rec.r_code_primary) > 0
        assert 'library' in rec.r_code_primary


class TestRecommendationScenarios:
    """Tests for different data scenarios."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_small_samples_prefer_deseq2(self, small_sample_profile):
        """Test that small samples prefer DESeq2."""
        recommender = PipelineRecommender(small_sample_profile)
        rec = recommender.get_recommendation()
        # DESeq2 should score highly for small samples
        assert rec.all_scores['DESeq2'] >= rec.all_scores['Wilcoxon']
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_large_samples_allow_wilcoxon(self, large_sample_profile):
        """Test that large samples allow Wilcoxon."""
        recommender = PipelineRecommender(large_sample_profile)
        rec = recommender.get_recommendation()
        # Wilcoxon should be viable for large samples
        assert rec.all_scores['Wilcoxon'] > 50
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_outliers_prefer_robust(self, outlier_profile):
        """Test that outliers prefer edgeR_robust."""
        recommender = PipelineRecommender(outlier_profile)
        rec = recommender.get_recommendation()
        # edgeR_robust should score highly with outliers
        assert rec.all_scores['edgeR_robust'] >= rec.all_scores['edgeR'] - 10
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_batch_effects_generate_warning(self, batch_effect_profile):
        """Test that batch effects generate warning."""
        recommender = PipelineRecommender(batch_effect_profile)
        rec = recommender.get_recommendation()
        # Should have warning about batch effects
        assert any('batch' in w.lower() for w in rec.warnings)


class TestRecommendationMethods:
    """Tests for Recommendation object methods."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_summary_returns_string(self, demo_profile):
        """Test that summary() returns string."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        summary = rec.summary()
        assert isinstance(summary, str)
        assert 'RAPTOR' in summary
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_to_dict_returns_dict(self, demo_profile):
        """Test that to_dict() returns dictionary."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        d = rec.to_dict()
        assert isinstance(d, dict)
        assert 'primary_pipeline' in d
        assert 'all_scores' in d


# =============================================================================
# Tests for Convenience Functions
# =============================================================================

class TestConvenienceFunctions:
    """Tests for convenience functions."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_recommend_pipeline(self, demo_profile):
        """Test recommend_pipeline function."""
        rec = recommend_pipeline(demo_profile)
        assert isinstance(rec, Recommendation)
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_get_pipeline_info(self):
        """Test get_pipeline_info function."""
        info = get_pipeline_info('DESeq2')
        assert isinstance(info, PipelineInfo)
        assert info.name == 'DESeq2'
        assert len(info.strengths) > 0
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_get_pipeline_info_invalid(self):
        """Test get_pipeline_info with invalid name."""
        info = get_pipeline_info('InvalidPipeline')
        assert info is None
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_list_pipelines(self):
        """Test list_pipelines function."""
        pipelines = list_pipelines()
        assert isinstance(pipelines, list)
        assert 'DESeq2' in pipelines
        assert 'edgeR' in pipelines
        assert len(pipelines) >= 4


# =============================================================================
# Tests for PipelineInfo
# =============================================================================

class TestPipelineInfo:
    """Tests for pipeline information."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_all_pipelines_have_info(self):
        """Test that all pipelines have complete info."""
        for name in list_pipelines():
            info = get_pipeline_info(name)
            assert info is not None
            assert info.name == name
            assert info.min_replicates >= 1
            assert len(info.r_packages) > 0
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_deseq2_info(self):
        """Test DESeq2 pipeline info."""
        info = get_pipeline_info('DESeq2')
        assert 'DESeq2' in info.r_packages
        assert info.min_replicates == 2
        assert 'shrinkage' in info.description.lower() or 'negative binomial' in info.description.lower()
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_wilcoxon_requires_8(self):
        """Test that Wilcoxon requires 8 replicates."""
        info = get_pipeline_info('Wilcoxon')
        assert info.min_replicates == 8


# =============================================================================
# Tests for ML Recommender
# =============================================================================

class TestMLRecommender:
    """Tests for ML-based recommender."""
    
    @pytest.mark.skipif(not ML_AVAILABLE, reason="ML module not available")
    def test_init(self):
        """Test ML recommender initialization."""
        recommender = MLPipelineRecommender(model_type='random_forest')
        assert recommender.model_type == 'random_forest'
        assert not recommender.is_trained
    
    @pytest.mark.skipif(not ML_AVAILABLE, reason="ML module not available")
    def test_init_gradient_boosting(self):
        """Test gradient boosting initialization."""
        recommender = MLPipelineRecommender(model_type='gradient_boosting')
        assert recommender.model_type == 'gradient_boosting'
    
    @pytest.mark.skipif(not ML_AVAILABLE, reason="ML module not available")
    def test_recommend_requires_trained(self, demo_profile):
        """Test that recommend requires trained model."""
        recommender = MLPipelineRecommender()
        with pytest.raises(RuntimeError):
            recommender.recommend(demo_profile)


# =============================================================================
# Tests for Simulation (if available)
# =============================================================================

class TestSimulation:
    """Tests for RNA-seq simulation."""
    
    @pytest.mark.skipif(not SIMULATION_AVAILABLE, reason="Simulation not available")
    def test_simulate_rnaseq(self):
        """Test basic simulation."""
        result = simulate_rnaseq(
            n_genes=1000,
            n_samples_per_group=3,
            de_fraction=0.10,
            seed=42
        )
        
        assert isinstance(result, SimulationResult)
        assert result.counts.shape[0] == 1000
        assert result.counts.shape[1] == 6  # 2 groups × 3 samples
    
    @pytest.mark.skipif(not SIMULATION_AVAILABLE, reason="Simulation not available")
    def test_simulation_has_ground_truth(self):
        """Test that simulation has ground truth."""
        result = simulate_rnaseq(
            n_genes=1000,
            n_samples_per_group=3,
            de_fraction=0.10,
            seed=42
        )
        
        assert 'is_de' in result.ground_truth.columns
        assert result.n_de == int(1000 * 0.10)
    
    @pytest.mark.skipif(not SIMULATION_AVAILABLE, reason="Simulation not available")
    def test_simulation_metadata(self):
        """Test that simulation creates correct metadata."""
        result = simulate_rnaseq(
            n_genes=1000,
            n_samples_per_group=4,
            seed=42
        )
        
        assert 'condition' in result.metadata.columns
        assert len(result.metadata) == 8


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for Module 4."""
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE or not PROFILER_AVAILABLE,
                        reason="Modules not available")
    def test_full_workflow(self, sample_count_matrix, sample_metadata):
        """Test complete Module 4 workflow."""
        # Profile
        profiler = RNAseqDataProfiler(sample_count_matrix, sample_metadata)
        profile = profiler.run_full_profile()
        
        # Recommend
        recommender = PipelineRecommender(profile)
        rec = recommender.get_recommendation()
        
        # Verify
        assert rec.primary_pipeline in PIPELINES
        assert rec.primary_score > 0
    
    @pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
    def test_serialization(self, demo_profile, output_dir):
        """Test recommendation serialization."""
        recommender = PipelineRecommender(demo_profile)
        rec = recommender.get_recommendation()
        
        # Save as JSON
        rec_file = output_dir / 'recommendation.json'
        with open(rec_file, 'w') as f:
            json.dump(rec.to_dict(), f, indent=2)
        
        # Load and verify
        with open(rec_file) as f:
            loaded = json.load(f)
        
        assert loaded['primary_pipeline'] == rec.primary_pipeline
        assert loaded['primary_score'] == rec.primary_score


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
@pytest.mark.parametrize("pipeline_name", ['DESeq2', 'edgeR', 'limma-voom', 'Wilcoxon', 'edgeR_robust'])
def test_pipeline_info_exists(pipeline_name):
    """Test that info exists for all pipelines."""
    info = get_pipeline_info(pipeline_name)
    assert info is not None
    assert info.name == pipeline_name


@pytest.mark.skipif(not RECOMMENDER_AVAILABLE, reason="Module not available")
@pytest.mark.parametrize("min_group_size", [2, 3, 5, 8, 15, 30])
def test_different_sample_sizes(min_group_size):
    """Test recommendations for different sample sizes."""
    class TestProfile:
        def __init__(self, n):
            self.n_samples = n * 2
            self.n_genes = 10000
            self.n_groups = 2
            self.min_group_size = n
            self.max_group_size = n
            self.sample_balance = 1.0
            self.has_replicates = n >= 2
            self.has_sufficient_replicates = n >= 3
            self.design_complexity = 'simple'
            self.bcv = 0.3
            self.bcv_category = 'moderate'
            self.common_dispersion = 0.09
            self.overdispersion_ratio = 2.0
            self.library_size_cv = 0.15
            self.library_size_range = 1.8
            self.sparsity = 0.4
            self.low_count_proportion = 0.1
            self.has_outliers = False
            self.outlier_severity = 'none'
            self.has_batch_effect = False
            self.batch_confounded = False
            self.quality_score = 75.0
    
    profile = TestProfile(min_group_size)
    recommender = PipelineRecommender(profile)
    rec = recommender.get_recommendation()
    
    assert rec.primary_pipeline in PIPELINES
    
    # Wilcoxon should only be recommended for n >= 8
    if min_group_size < 8:
        assert rec.all_scores['Wilcoxon'] < rec.all_scores['DESeq2']
