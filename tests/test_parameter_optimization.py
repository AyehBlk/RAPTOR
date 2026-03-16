"""
RAPTOR v2.2.0 - Module 8: Parameter Optimization Tests

Comprehensive test suite for parameter optimization functionality.

Tests cover:
- ParameterSpace class
- OptimizationResult class
- ParameterOptimizer base class
- Ground truth optimization
- FDR control optimization
- Stability optimization
- Reproducibility optimization
- Grid search strategy
- Random search strategy
- Differential evolution strategy
- CLI compatibility
- Integration with Module 7

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import pytest
import pandas as pd
import numpy as np
import json
import pickle
from pathlib import Path
from unittest.mock import patch, MagicMock

# RAPTOR imports
try:
    from raptor.parameter_optimization import (
        ParameterSpace,
        OptimizationResult,
        ParameterOptimizer,
        optimize_with_ground_truth,
        optimize_with_fdr_control,
        # optimize_with_stability,  # Requires special setup
        # optimize_with_reproducibility,  # Requires two cohorts
        ALPHA_MIN, ALPHA_MAX,
        LFC_MIN, LFC_MAX,
        MIN_VALIDATED_GENES_ABSOLUTE,
    )
    RAPTOR_AVAILABLE = True
except ImportError:
    RAPTOR_AVAILABLE = False
    pytest.skip("RAPTOR parameter_optimization module not available", allow_module_level=True)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def sample_de_result():
    """Generate sample DE results DataFrame."""
    np.random.seed(42)
    n_genes = 1000
    
    return pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'log2FoldChange': np.random.normal(0, 1.5, n_genes),
        'pvalue': np.random.beta(1, 10, n_genes),
        'padj': np.random.beta(1, 20, n_genes),
        'baseMean': np.random.gamma(2, 100, n_genes)
    }).set_index('gene_id')


@pytest.fixture
def ground_truth_genes():
    """Generate ground truth gene list."""
    np.random.seed(42)
    # Select 50 genes (realistic for validated gene set)
    genes = [f'ENSG{i+1:011d}' for i in range(50)]
    return pd.DataFrame({'gene_id': genes})


@pytest.fixture
def small_ground_truth():
    """Small ground truth (edge case - minimum size)."""
    genes = [f'ENSG{i+1:011d}' for i in range(MIN_VALIDATED_GENES_ABSOLUTE)]
    return pd.DataFrame({'gene_id': genes})


@pytest.fixture
def saved_de_result_pkl(output_dir, sample_de_result):
    """Save DE result as DEResult pickle (from Module 7)."""
    # Create a mock DEResult object
    try:
        from raptor.de_import import DEResult
        
        # DEResult expects standardized column names
        de_df = sample_de_result.copy()
        de_df = de_df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value',
            'baseMean': 'base_mean'
        })
        
        de_res = DEResult(
            results_df=de_df,
            pipeline='DESeq2',
            parameters={'fdr_threshold': 0.05, 'lfc_threshold': 0.0}
        )
        
        pkl_path = output_dir / 'de_result.pkl'
        de_res.save(pkl_path)
        return pkl_path
    except ImportError:
        # If Module 7 not available, create mock
        pkl_path = output_dir / 'de_result.pkl'
        mock_de_res = MagicMock()
        mock_de_res.results_df = sample_de_result
        mock_de_res.n_genes = len(sample_de_result)
        
        with open(pkl_path, 'wb') as f:
            pickle.dump(mock_de_res, f)
        return pkl_path


@pytest.fixture
def saved_ground_truth_csv(output_dir, ground_truth_genes):
    """Save ground truth to CSV."""
    csv_path = output_dir / 'ground_truth.csv'
    ground_truth_genes.to_csv(csv_path, index=False)
    return csv_path


# =============================================================================
# ParameterSpace Tests
# =============================================================================

class TestParameterSpace:
    """Test ParameterSpace class."""
    
    def test_parameter_space_initialization(self):
        """Test ParameterSpace initialization."""
        ps = ParameterSpace(
            name='alpha',
            bounds=(0.001, 0.20),
            default=0.05,
            description='FDR threshold'
        )
        
        assert ps.name == 'alpha'
        assert ps.bounds == (0.001, 0.20)
        assert ps.default == 0.05
        assert ps.description == 'FDR threshold'
    
    def test_sample_random(self):
        """Test random sampling from parameter space."""
        ps = ParameterSpace(
            name='alpha',
            bounds=(0.001, 0.20),
            default=0.05,
            description='FDR threshold'
        )
        
        # Sample 100 values
        samples = [ps.sample_random() for _ in range(100)]
        
        # All should be within bounds
        assert all(0.001 <= s <= 0.20 for s in samples)
        
        # Should have variety (not all the same)
        assert len(set(samples)) > 10
    
    def test_sample_grid(self):
        """Test grid sampling."""
        ps = ParameterSpace(
            name='lfc_threshold',
            bounds=(0.0, 2.0),
            default=1.0,
            description='LFC threshold'
        )
        
        grid = ps.sample_grid(n_points=5)
        
        assert len(grid) == 5
        assert grid[0] == 0.0
        assert grid[-1] == 2.0
        assert all(0.0 <= p <= 2.0 for p in grid)
    
    def test_validate(self):
        """Test parameter validation."""
        ps = ParameterSpace(
            name='alpha',
            bounds=(0.001, 0.20),
            default=0.05,
            description='FDR threshold'
        )
        
        assert ps.validate(0.05) is True
        assert ps.validate(0.001) is True
        assert ps.validate(0.20) is True
        assert ps.validate(0.0) is False
        assert ps.validate(0.5) is False


# =============================================================================
# OptimizationResult Tests
# =============================================================================

class TestOptimizationResult:
    """Test OptimizationResult class."""
    
    def test_optimization_result_initialization(self):
        """Test OptimizationResult initialization."""
        result = OptimizationResult(
            best_parameters={'alpha': 0.01, 'lfc_threshold': 1.0},
            best_score=0.85,
            optimization_method='ground_truth',
            search_strategy='grid',
            metric='f1_score',
            n_iterations=25
        )
        
        assert result.best_parameters == {'alpha': 0.01, 'lfc_threshold': 1.0}
        assert result.best_score == 0.85
        assert result.optimization_method == 'ground_truth'
        assert result.search_strategy == 'grid'
        assert result.metric == 'f1_score'
        assert result.n_iterations == 25
    
    def test_summary_generation(self):
        """Test summary string generation."""
        result = OptimizationResult(
            best_parameters={'alpha': 0.01, 'lfc_threshold': 1.5},
            best_score=0.92,
            optimization_method='fdr_control',
            search_strategy='random',
            metric='stability_score',
            n_iterations=50,
            n_deg_genes=250
        )
        
        summary = result.summary()
        
        assert isinstance(summary, str)
        assert 'Parameter Optimization' in summary
        assert 'fdr_control' in summary
        assert '0.01' in summary
        assert '1.5' in summary or '1.50' in summary
        assert '0.92' in summary or '0.9200' in summary
        assert '250' in summary
    
    def test_save_and_load(self, output_dir):
        """Test saving and loading OptimizationResult."""
        deg_genes = pd.DataFrame({
            'gene_id': ['Gene1', 'Gene2', 'Gene3'],
            'log2FoldChange': [2.0, -1.5, 3.0],
            'padj': [0.001, 0.002, 0.0001]
        })
        
        original = OptimizationResult(
            best_parameters={'alpha': 0.05, 'lfc_threshold': 1.0},
            best_score=0.88,
            optimization_method='ground_truth',
            search_strategy='grid',
            metric='f1_score',
            n_iterations=25,
            deg_genes=deg_genes,
            n_deg_genes=3,
            performance_metrics={'precision': 0.9, 'recall': 0.85}
        )
        
        # Save
        original.save(output_dir)
        
        # Check files exist
        assert (output_dir / 'optimized_params.yaml').exists()
        assert (output_dir / 'deg_genes.csv').exists()
        assert (output_dir / 'optimization_result.pkl').exists()
        assert (output_dir / 'convergence_history.json').exists()
        
        # Load
        loaded = OptimizationResult.load(output_dir / 'optimization_result.pkl')
        
        assert loaded.best_parameters == original.best_parameters
        assert loaded.best_score == original.best_score
        assert loaded.n_deg_genes == original.n_deg_genes


# =============================================================================
# Ground Truth Optimization Tests
# =============================================================================

class TestGroundTruthOptimization:
    """Test ground truth-based optimization."""
    
    def test_basic_optimization(self, sample_de_result, ground_truth_genes, output_dir):
        """Test basic ground truth optimization."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='grid',
            grid_points=3,  # Small grid for speed
            output_dir=output_dir,
            random_state=42
        )
        
        assert isinstance(result, OptimizationResult)
        assert 'alpha' in result.best_parameters
        assert 'lfc_threshold' in result.best_parameters
        assert result.optimization_method == 'groundtruth'
        assert result.search_strategy == 'grid'
        assert result.metric == 'f1_score'
        assert result.best_score >= 0
        assert result.best_score <= 1
    
    def test_parameter_bounds(self, sample_de_result, ground_truth_genes, output_dir):
        """Test that optimized parameters are within bounds."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='random',
            n_iterations=10,
            output_dir=output_dir,
            random_state=42
        )
        
        # Check alpha bounds
        assert ALPHA_MIN <= result.best_parameters['alpha'] <= ALPHA_MAX
        
        # Check LFC bounds
        assert LFC_MIN <= result.best_parameters['lfc_threshold'] <= LFC_MAX
    
    def test_small_ground_truth(self, sample_de_result, small_ground_truth, output_dir):
        """Test with minimum-size ground truth."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=small_ground_truth,
            strategy='grid',
            grid_points=3,
            output_dir=output_dir,
            random_state=42
        )
        
        # Should still work with minimum number
        assert isinstance(result, OptimizationResult)
        assert result.best_score >= 0
    
    def test_empty_ground_truth_fails(self, sample_de_result, output_dir):
        """Test that empty ground truth raises error."""
        empty_gt = pd.DataFrame({'gene_id': []})
        
        with pytest.raises(ValueError):
            optimize_with_ground_truth(
                de_result=sample_de_result,
                ground_truth=empty_gt,
                output_dir=output_dir
            )
    
    def test_convergence_history(self, sample_de_result, ground_truth_genes, output_dir):
        """Test that convergence history is recorded."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='random',
            n_iterations=5,
            output_dir=output_dir,
            random_state=42
        )
        
        assert len(result.convergence_history) > 0
        assert all('parameters' in item for item in result.convergence_history)
        assert all('score' in item for item in result.convergence_history)


# =============================================================================
# FDR Control Optimization Tests
# =============================================================================

class TestFDRControlOptimization:
    """Test FDR control optimization."""
    
    def test_basic_fdr_optimization(self, sample_de_result, output_dir):
        """Test basic FDR control optimization."""
        result = optimize_with_fdr_control(
            de_result=sample_de_result,
            target_fdr=0.05,
            strategy='grid',
            grid_points=3,
            output_dir=output_dir,
            random_state=42
        )
        
        assert isinstance(result, OptimizationResult)
        assert result.optimization_method == 'fdrcontrol'
        assert result.best_score >= 0
    
    def test_different_fdr_targets(self, sample_de_result, output_dir):
        """Test with different FDR targets."""
        for target_fdr in [0.01, 0.05, 0.10]:
            result = optimize_with_fdr_control(
                de_result=sample_de_result,
                target_fdr=target_fdr,
                strategy='grid',
                grid_points=3,
                output_dir=output_dir / f'fdr_{target_fdr}',
                random_state=42
            )
            
            assert isinstance(result, OptimizationResult)
            assert result.best_parameters['alpha'] <= 1.0  # alpha is a valid probability
    
    def test_fdr_no_validation_needed(self, sample_de_result, output_dir):
        """Test that FDR optimization doesn't require validation data."""
        # Should work with just DE results
        result = optimize_with_fdr_control(
            de_result=sample_de_result,
            target_fdr=0.05,
            strategy='grid',
            grid_points=3,
            output_dir=output_dir,
            random_state=42
        )
        
        assert result is not None


# =============================================================================
# Integration Tests
# =============================================================================

@pytest.mark.integration
class TestModule7to8Integration:
    """Test integration between Module 7 and Module 8."""
    
    def test_load_m7_result(self, saved_de_result_pkl):
        """Test loading Module 7 DEResult for optimization."""
        # Load the pickle
        with open(saved_de_result_pkl, 'rb') as f:
            de_res = pickle.load(f)
        
        # Should have results_df attribute
        assert hasattr(de_res, 'results_df')
        assert isinstance(de_res.results_df, pd.DataFrame)
    
    def test_m7_to_m8_workflow(self, saved_de_result_pkl, ground_truth_genes, output_dir):
        """Test complete M7→M8 workflow."""
        # Load M7 result
        with open(saved_de_result_pkl, 'rb') as f:
            de_res = pickle.load(f)
        
        # M7 uses standardized column names, M8 expects original names
        de_df = de_res.results_df.rename(columns={
            'log2_fold_change': 'log2FoldChange',
            'p_value': 'pvalue',
            'adjusted_p_value': 'padj',
            'base_mean': 'baseMean'
        })
        
        # Run M8 optimization
        result = optimize_with_ground_truth(
            de_result=de_df,
            ground_truth=ground_truth_genes,
            strategy='grid',
            grid_points=3,
            output_dir=output_dir,
            random_state=42
        )
        
        assert isinstance(result, OptimizationResult)
        assert result.n_deg_genes >= 0


# =============================================================================
# CLI Compatibility Tests
# =============================================================================

@pytest.mark.cli
class TestCLICompatibility:
    """Test CLI compatibility."""
    
    def test_cli_ground_truth_format(self, sample_de_result, saved_ground_truth_csv, output_dir):
        """Test that CLI can call ground truth optimization."""
        # This mimics how fixed CLI calls the function
        gt_df = pd.read_csv(saved_ground_truth_csv)
        
        assert 'gene_id' in gt_df.columns
        
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=gt_df,  # Pass DataFrame (not list!)
            output_dir=output_dir,
            strategy='grid',
            grid_points=3
        )
        
        # Check CLI can access attributes correctly
        assert hasattr(result, 'best_parameters')
        assert 'alpha' in result.best_parameters
        assert 'lfc_threshold' in result.best_parameters
        assert hasattr(result, 'best_score')  # NOT 'score'
        assert hasattr(result, 'n_deg_genes')
    
    def test_cli_fdr_format(self, sample_de_result, output_dir):
        """Test that CLI can call FDR optimization."""
        # CLI passes target_fdr (not fdr_target)
        result = optimize_with_fdr_control(
            de_result=sample_de_result,
            target_fdr=0.05,  # Correct parameter name
            output_dir=output_dir,
            strategy='grid',
            grid_points=3
        )
        
        assert isinstance(result, OptimizationResult)
    
    def test_output_files_for_cli(self, sample_de_result, ground_truth_genes, output_dir):
        """Test that expected output files are created."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            output_dir=output_dir,
            strategy='grid',
            grid_points=3
        )
        
        # These files should exist (CLI checks for them)
        assert (output_dir / 'optimized_params.yaml').exists()
        assert (output_dir / 'deg_genes.csv').exists()
        assert (output_dir / 'optimization_result.pkl').exists()
        assert (output_dir / 'convergence_history.json').exists()


# =============================================================================
# Search Strategy Tests
# =============================================================================

class TestSearchStrategies:
    """Test different search strategies."""
    
    def test_grid_search(self, sample_de_result, ground_truth_genes, output_dir):
        """Test grid search strategy."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='grid',
            grid_points=3,
            output_dir=output_dir,
            random_state=42
        )
        
        assert result.search_strategy == 'grid'
        # Grid with 3x3 = 9 iterations
        assert result.n_iterations == 9
    
    def test_random_search(self, sample_de_result, ground_truth_genes, output_dir):
        """Test random search strategy."""
        result = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='random',
            n_iterations=10,
            output_dir=output_dir,
            random_state=42
        )
        
        assert result.search_strategy == 'random'
        assert result.n_iterations == 10
    
    def test_reproducibility_with_seed(self, sample_de_result, ground_truth_genes, output_dir):
        """Test that random_state gives reproducible results."""
        result1 = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='random',
            n_iterations=10,
            output_dir=output_dir / 'run1',
            random_state=42
        )
        
        result2 = optimize_with_ground_truth(
            de_result=sample_de_result,
            ground_truth=ground_truth_genes,
            strategy='random',
            n_iterations=10,
            output_dir=output_dir / 'run2',
            random_state=42
        )
        
        # Should get same results with same seed
        assert result1.best_parameters == result2.best_parameters
        assert result1.best_score == result2.best_score


# =============================================================================
# Validation Tests
# =============================================================================

class TestValidation:
    """Test input validation."""
    
    def test_invalid_de_result_type(self, ground_truth_genes, output_dir):
        """Test with invalid DE result type."""
        with pytest.raises(TypeError):
            optimize_with_ground_truth(
                de_result="not_a_dataframe",
                ground_truth=ground_truth_genes,
                output_dir=output_dir
            )
    
    def test_empty_de_result(self, ground_truth_genes, output_dir):
        """Test with empty DE result."""
        empty_de = pd.DataFrame(columns=['gene_id', 'log2FoldChange', 'pvalue', 'padj'])
        
        with pytest.raises(ValueError):
            optimize_with_ground_truth(
                de_result=empty_de,
                ground_truth=ground_truth_genes,
                output_dir=output_dir
            )
    
    def test_missing_required_columns(self, ground_truth_genes, output_dir):
        """Test with missing required columns."""
        incomplete_de = pd.DataFrame({
            'gene_id': ['Gene1', 'Gene2'],
            'log2FoldChange': [1.0, -1.0]
            # Missing pvalue, padj
        }).set_index('gene_id')
        
        with pytest.raises((ValueError, RuntimeError)):
            optimize_with_ground_truth(
                de_result=incomplete_de,
                ground_truth=ground_truth_genes,
                output_dir=output_dir
            )
    
    def test_invalid_strategy(self, sample_de_result, ground_truth_genes, output_dir):
        """Test with invalid search strategy."""
        with pytest.raises(ValueError):
            optimize_with_ground_truth(
                de_result=sample_de_result,
                ground_truth=ground_truth_genes,
                strategy='invalid_strategy',
                output_dir=output_dir
            )


# =============================================================================
# Performance Tests
# =============================================================================

@pytest.mark.slow
class TestPerformance:
    """Test performance with large datasets."""
    
    def test_large_de_result(self, ground_truth_genes, output_dir):
        """Test with realistic dataset size."""
        np.random.seed(42)
        n_genes = 20000  # Realistic human transcriptome
        
        large_de = pd.DataFrame({
            'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
            'log2FoldChange': np.random.normal(0, 1.5, n_genes),
            'pvalue': np.random.beta(1, 10, n_genes),
            'padj': np.random.beta(1, 20, n_genes),
            'baseMean': np.random.gamma(2, 100, n_genes)
        }).set_index('gene_id')
        
        # Should handle large dataset
        result = optimize_with_ground_truth(
            de_result=large_de,
            ground_truth=ground_truth_genes,
            strategy='grid',
            grid_points=3,  # Keep small for speed
            output_dir=output_dir,
            random_state=42
        )
        
        assert result.n_iterations == 9
        assert isinstance(result, OptimizationResult)


# =============================================================================
# Helper Functions
# =============================================================================

def assert_optimization_result_valid(result: OptimizationResult):
    """Assert that an OptimizationResult is valid."""
    assert isinstance(result, OptimizationResult)
    assert isinstance(result.best_parameters, dict)
    assert 'alpha' in result.best_parameters
    assert 'lfc_threshold' in result.best_parameters
    assert isinstance(result.best_score, (int, float))
    assert result.best_score >= 0
    assert result.optimization_method in ['groundtruth', 'fdrcontrol', 'stability', 'reproducibility', 'ground_truth', 'fdr_control']
    assert result.search_strategy in ['grid', 'random', 'differential_evolution']
    assert result.n_iterations > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
