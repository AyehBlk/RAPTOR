#!/usr/bin/env python3
"""
Unit Tests for Adaptive Threshold Optimizer (ATO)

Run with: pytest test_threshold_optimizer.py -v

Author: Ayeh Bolouki
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from raptor.threshold_optimizer import (
    AdaptiveThresholdOptimizer,
    ThresholdResult,
    optimize_thresholds
)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def sample_df():
    """Create sample DE results DataFrame."""
    np.random.seed(42)
    n = 1000
    
    df = pd.DataFrame({
        'log2FoldChange': np.concatenate([
            np.random.normal(0, 0.3, 850),  # Null genes
            np.random.normal(2, 0.5, 75),    # Upregulated
            np.random.normal(-2, 0.5, 75)    # Downregulated
        ]),
        'pvalue': np.concatenate([
            np.random.uniform(0.1, 1, 850),
            np.random.uniform(0.0001, 0.01, 150)
        ]),
        'baseMean': np.random.exponential(500, n),
        'lfcSE': np.abs(np.random.normal(0.1, 0.05, n))
    })
    df.index = [f'Gene_{i}' for i in range(n)]
    
    return df


@pytest.fixture
def deseq2_like_df():
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
def edger_like_df():
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


# ============================================================================
# Basic Functionality Tests
# ============================================================================

class TestBasicFunctionality:
    """Test basic ATO functionality."""
    
    def test_initialization(self, sample_df):
        """Test ATO initialization."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        assert ato is not None
        assert len(ato.df) <= len(sample_df)  # May remove NaN
    
    def test_column_standardization_deseq2(self, deseq2_like_df):
        """Test column name standardization for DESeq2."""
        ato = AdaptiveThresholdOptimizer(deseq2_like_df, verbose=False)
        assert 'logfc' in ato.df.columns
        assert 'pvalue' in ato.df.columns
    
    def test_column_standardization_edger(self, edger_like_df):
        """Test column name standardization for edgeR."""
        ato = AdaptiveThresholdOptimizer(edger_like_df, verbose=False)
        assert 'logfc' in ato.df.columns
        assert 'pvalue' in ato.df.columns
    
    def test_missing_required_columns(self):
        """Test error on missing required columns."""
        df = pd.DataFrame({'somecolumn': [1, 2, 3]})
        with pytest.raises(ValueError):
            AdaptiveThresholdOptimizer(df, verbose=False)
    
    def test_optimize_returns_result(self, sample_df):
        """Test that optimize() returns ThresholdResult."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        result = ato.optimize()
        assert isinstance(result, ThresholdResult)
    
    def test_convenience_function(self, sample_df):
        """Test optimize_thresholds convenience function."""
        result = optimize_thresholds(sample_df, verbose=False)
        assert isinstance(result, ThresholdResult)


# ============================================================================
# P-Value Adjustment Tests
# ============================================================================

class TestPValueAdjustment:
    """Test p-value adjustment methods."""
    
    def test_benjamini_hochberg(self, sample_df):
        """Test BH adjustment."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        pvals = sample_df['pvalue'].values
        adjusted = ato.benjamini_hochberg(pvals)
        
        assert len(adjusted) == len(pvals)
        assert all(adjusted >= pvals)  # Adjusted >= raw
        assert all(adjusted <= 1)
        assert all(adjusted >= 0)
    
    def test_benjamini_yekutieli(self, sample_df):
        """Test BY adjustment."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        pvals = sample_df['pvalue'].values
        adjusted = ato.benjamini_yekutieli(pvals)
        
        # BY should be more conservative than BH
        bh_adjusted = ato.benjamini_hochberg(pvals)
        assert all(adjusted >= bh_adjusted - 1e-10)  # Small tolerance for floating point
    
    def test_holm(self, sample_df):
        """Test Holm step-down."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        pvals = sample_df['pvalue'].values
        adjusted = ato.holm(pvals)
        
        assert len(adjusted) == len(pvals)
        assert all(adjusted <= 1)
    
    def test_bonferroni(self, sample_df):
        """Test Bonferroni correction."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        pvals = np.array([0.01, 0.02, 0.03])
        adjusted = ato.bonferroni(pvals)
        
        expected = np.array([0.03, 0.06, 0.09])
        np.testing.assert_array_almost_equal(adjusted, expected)
    
    def test_storey_qvalue(self, sample_df):
        """Test Storey q-value."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        pvals = sample_df['pvalue'].values
        qvals, pi0 = ato.storey_qvalue(pvals)
        
        assert len(qvals) == len(pvals)
        assert 0 < pi0 <= 1
        assert all(qvals <= 1)


# ============================================================================
# LogFC Optimization Tests
# ============================================================================

class TestLogFCOptimization:
    """Test logFC cutoff optimization methods."""
    
    def test_mad_method(self, sample_df):
        """Test MAD-based optimization."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        cutoff, reasoning = ato.optimize_logfc_mad()
        
        assert cutoff > 0
        assert isinstance(reasoning, str)
        assert len(reasoning) > 0
    
    def test_mixture_method(self, sample_df):
        """Test mixture model optimization."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        cutoff, reasoning = ato.optimize_logfc_mixture()
        
        assert cutoff > 0
        assert isinstance(reasoning, str)
    
    def test_power_method(self, sample_df):
        """Test power-based optimization."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        cutoff, reasoning = ato.optimize_logfc_power(n1=3, n2=3)
        
        assert cutoff > 0
        assert isinstance(reasoning, str)
    
    def test_percentile_method(self, sample_df):
        """Test percentile-based optimization."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        cutoff, reasoning = ato.optimize_logfc_percentile()
        
        assert cutoff > 0
        assert isinstance(reasoning, str)
    
    def test_auto_method(self, sample_df):
        """Test auto/consensus method."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        result = ato.optimize(logfc_method='auto')
        
        assert result.logfc_cutoff > 0
        assert 'consensus' in result.logfc_method.lower() or 'Consensus' in result.logfc_method


# ============================================================================
# Analysis Goal Tests
# ============================================================================

class TestAnalysisGoals:
    """Test different analysis goals."""
    
    def test_discovery_goal(self, sample_df):
        """Test discovery goal."""
        result = optimize_thresholds(sample_df, goal='discovery', verbose=False)
        assert result.padj_method in ['BH', 'qvalue', 'BY']
    
    def test_validation_goal(self, sample_df):
        """Test validation goal."""
        result = optimize_thresholds(sample_df, goal='validation', verbose=False)
        assert result.padj_method in ['holm', 'bonferroni', 'hochberg']
    
    def test_balanced_goal(self, sample_df):
        """Test balanced goal."""
        result = optimize_thresholds(sample_df, goal='balanced', verbose=False)
        assert result is not None
    
    def test_invalid_goal(self, sample_df):
        """Test invalid goal raises error."""
        with pytest.raises(ValueError):
            AdaptiveThresholdOptimizer(sample_df, goal='invalid', verbose=False)


# ============================================================================
# Result Object Tests
# ============================================================================

class TestThresholdResult:
    """Test ThresholdResult object."""
    
    def test_result_attributes(self, sample_df):
        """Test all result attributes are present."""
        result = optimize_thresholds(sample_df, verbose=False)
        
        assert hasattr(result, 'logfc_cutoff')
        assert hasattr(result, 'padj_cutoff')
        assert hasattr(result, 'padj_method')
        assert hasattr(result, 'logfc_method')
        assert hasattr(result, 'logfc_reasoning')
        assert hasattr(result, 'padj_reasoning')
        assert hasattr(result, 'n_significant_optimized')
        assert hasattr(result, 'n_significant_traditional')
        assert hasattr(result, 'pi0_estimate')
    
    def test_result_summary(self, sample_df):
        """Test summary() method."""
        result = optimize_thresholds(sample_df, verbose=False)
        summary = result.summary()
        
        assert isinstance(summary, str)
        assert len(summary) > 100
        assert 'THRESHOLD' in summary.upper()
    
    def test_result_values_reasonable(self, sample_df):
        """Test result values are reasonable."""
        result = optimize_thresholds(sample_df, verbose=False)
        
        assert 0 < result.logfc_cutoff < 10
        assert 0 < result.padj_cutoff <= 0.1
        assert 0 < result.pi0_estimate <= 1
        assert result.n_significant_optimized >= 0
        assert result.n_significant_traditional >= 0


# ============================================================================
# Helper Method Tests
# ============================================================================

class TestHelperMethods:
    """Test helper methods."""
    
    def test_get_significant_genes(self, sample_df):
        """Test get_significant_genes method."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        ato.optimize()
        
        sig_genes = ato.get_significant_genes()
        assert isinstance(sig_genes, pd.DataFrame)
        assert len(sig_genes) <= len(sample_df)
    
    def test_compare_thresholds(self, sample_df):
        """Test compare_thresholds method."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        ato.optimize()
        
        comparison = ato.compare_thresholds()
        assert isinstance(comparison, pd.DataFrame)
        assert 'logFC_cutoff' in comparison.columns
        assert 'padj_cutoff' in comparison.columns
        assert 'n_significant' in comparison.columns
    
    def test_adjustment_comparison(self, sample_df):
        """Test get_adjustment_comparison method."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        
        adj_comp = ato.get_adjustment_comparison()
        assert isinstance(adj_comp, pd.DataFrame)
        assert 'BH' in adj_comp.columns
        assert 'BY' in adj_comp.columns
        assert 'qvalue' in adj_comp.columns
    
    def test_estimate_pi0(self, sample_df):
        """Test π₀ estimation."""
        ato = AdaptiveThresholdOptimizer(sample_df, verbose=False)
        
        pi0_storey = ato.estimate_pi0(method='storey')
        pi0_pounds = ato.estimate_pi0(method='pounds')
        pi0_hist = ato.estimate_pi0(method='histogram')
        
        for pi0 in [pi0_storey, pi0_pounds, pi0_hist]:
            assert 0 < pi0 <= 1


# ============================================================================
# Edge Case Tests
# ============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_small_dataset(self):
        """Test with very small dataset."""
        df = pd.DataFrame({
            'log2FoldChange': [0.5, -0.5, 2.0],
            'pvalue': [0.1, 0.2, 0.001]
        })
        
        result = optimize_thresholds(df, verbose=False)
        assert result is not None
    
    def test_all_significant(self):
        """Test when all genes are significant."""
        df = pd.DataFrame({
            'log2FoldChange': np.random.normal(2, 0.5, 100),
            'pvalue': np.random.uniform(0.0001, 0.01, 100)
        })
        
        result = optimize_thresholds(df, verbose=False)
        assert result.pi0_estimate < 0.5  # Should detect low π₀
    
    def test_none_significant(self):
        """Test when no genes are significant."""
        df = pd.DataFrame({
            'log2FoldChange': np.random.normal(0, 0.1, 100),
            'pvalue': np.random.uniform(0.5, 1, 100)
        })
        
        result = optimize_thresholds(df, verbose=False)
        assert result.pi0_estimate > 0.8  # Should detect high π₀
    
    def test_with_nan_values(self):
        """Test handling of NaN values."""
        df = pd.DataFrame({
            'log2FoldChange': [0.5, np.nan, 2.0, -1.0],
            'pvalue': [0.1, 0.2, np.nan, 0.001]
        })
        
        ato = AdaptiveThresholdOptimizer(df, verbose=False)
        assert len(ato.df) < 4  # Should have removed NaN rows
    
    def test_with_zero_pvalues(self):
        """Test handling of zero p-values."""
        df = pd.DataFrame({
            'log2FoldChange': [0.5, 1.0, 2.0],
            'pvalue': [0.0, 0.1, 0.001]  # Contains 0
        })
        
        ato = AdaptiveThresholdOptimizer(df, verbose=False)
        # Should handle without error
        result = ato.optimize()
        assert result is not None


# ============================================================================
# Run tests
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
