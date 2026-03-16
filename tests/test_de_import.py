"""
RAPTOR v2.2.0 - Module 7: DE Import Tests

Comprehensive test suite for differential expression results import.

Tests cover:
- DEResult class functionality
- Import from all 4 methods (DESeq2, edgeR, limma, Wilcoxon)
- Column mapping and standardization
- Pipeline auto-detection
- Threshold filtering
- Serialization (pickle, JSON)
- Performance metrics calculation
- Integration with Module 6 outputs
- CLI compatibility

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
    from raptor.de_import import (
        DEResult,
        import_de_results,
        detect_pipeline,
        standardize_columns,
        calculate_significance,
        COLUMN_MAPPINGS,
        REQUIRED_COLUMNS
    )
    from raptor.utils.validation import ValidationError
    RAPTOR_AVAILABLE = True
except ImportError:
    RAPTOR_AVAILABLE = False
    pytest.skip("RAPTOR de_import module not available", allow_module_level=True)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def deseq2_results():
    """Generate mock DESeq2 results."""
    np.random.seed(42)
    n_genes = 1000
    
    return pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'baseMean': np.random.gamma(2, 100, n_genes),
        'log2FoldChange': np.random.normal(0, 1.5, n_genes),
        'lfcSE': np.random.uniform(0.1, 0.5, n_genes),
        'stat': np.random.normal(0, 3, n_genes),
        'pvalue': np.random.beta(1, 10, n_genes),
        'padj': np.random.beta(1, 20, n_genes)
    })


@pytest.fixture
def edger_results():
    """Generate mock edgeR results."""
    np.random.seed(42)
    n_genes = 1000
    
    return pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'logFC': np.random.normal(0, 1.5, n_genes),
        'logCPM': np.random.uniform(0, 10, n_genes),
        'LR': np.abs(np.random.normal(0, 5, n_genes)),
        'PValue': np.random.beta(1, 10, n_genes),
        'FDR': np.random.beta(1, 20, n_genes)
    })


@pytest.fixture
def limma_results():
    """Generate mock limma results."""
    np.random.seed(42)
    n_genes = 1000
    
    return pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'logFC': np.random.normal(0, 1.5, n_genes),
        'AveExpr': np.random.uniform(0, 10, n_genes),
        't': np.random.normal(0, 3, n_genes),
        'P.Value': np.random.beta(1, 10, n_genes),
        'adj.P.Val': np.random.beta(1, 20, n_genes),
        'B': np.random.normal(0, 2, n_genes)
    })


@pytest.fixture
def wilcoxon_results():
    """Generate mock Wilcoxon test results."""
    np.random.seed(42)
    n_genes = 1000
    
    return pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'log2FC': np.random.normal(0, 1.5, n_genes),
        'AveExpr': np.random.uniform(0, 10, n_genes),
        'W': np.abs(np.random.normal(0, 1000, n_genes)),
        'pvalue': np.random.beta(1, 10, n_genes),
        'padj': np.random.beta(1, 20, n_genes)
    })


@pytest.fixture
def deseq2_results_with_significance(deseq2_results):
    """DESeq2 results with pre-calculated significance."""
    df = deseq2_results.copy()
    df['is_significant'] = (df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 0)
    df['direction'] = np.where(
        df['log2FoldChange'] > 0, 'up',
        np.where(df['log2FoldChange'] < 0, 'down', 'unchanged')
    )
    df.loc[~df['is_significant'], 'direction'] = 'unchanged'
    return df


@pytest.fixture
def saved_deseq2_csv(output_dir, deseq2_results):
    """Save DESeq2 results to CSV."""
    csv_path = output_dir / 'deseq2_results.csv'
    deseq2_results.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def saved_edger_csv(output_dir, edger_results):
    """Save edgeR results to CSV."""
    csv_path = output_dir / 'edger_results.csv'
    edger_results.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def saved_limma_tsv(output_dir, limma_results):
    """Save limma results to TSV."""
    tsv_path = output_dir / 'limma_results.tsv'
    limma_results.to_csv(tsv_path, sep='\t', index=False)
    return tsv_path


@pytest.fixture
def ground_truth_de():
    """Generate ground truth DE status for testing metrics."""
    np.random.seed(42)
    n_genes = 1000
    
    # 100 true DE genes (10%)
    is_de = np.zeros(n_genes, dtype=bool)
    is_de[:100] = True
    np.random.shuffle(is_de)
    
    return pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'is_de': is_de,
        'true_lfc': np.where(is_de, np.random.normal(2, 0.5, n_genes), 0)
    }).set_index('gene_id')


@pytest.fixture
def module6_deseq2_output(output_dir):
    """
    Simulate Module 6 DESeq2 R script output.
    
    This is the exact format that run_deseq2.R produces,
    including the pre-calculated direction and is_significant columns.
    """
    np.random.seed(42)
    n_genes = 1000
    
    # Generate realistic DE results
    padj = np.random.beta(1, 20, n_genes)
    log2fc = np.random.normal(0, 1.5, n_genes)
    
    # Calculate significance (FDR < 0.05, LFC > 0)
    is_sig = (padj < 0.05) & (np.abs(log2fc) > 0)
    
    # Calculate direction
    direction = np.where(
        log2fc > 0, 'up',
        np.where(log2fc < 0, 'down', 'unchanged')
    )
    direction[~is_sig] = 'unchanged'
    
    df = pd.DataFrame({
        'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
        'baseMean': np.random.gamma(2, 100, n_genes),
        'log2FoldChange': log2fc,
        'lfcSE': np.random.uniform(0.1, 0.5, n_genes),
        'stat': np.random.normal(0, 3, n_genes),
        'pvalue': np.random.beta(1, 10, n_genes),
        'padj': padj,
        'direction': direction,
        'is_significant': is_sig
    })
    
    # Save as CSV (matching R script output)
    csv_path = output_dir / 'de_results.csv'
    df.to_csv(csv_path, index=False)
    
    # Also create summary JSON (as R script does)
    summary = {
        'timestamp': '2026-02-24T10:00:00',
        'raptor_version': '2.2.0',
        'module': 'M6',
        'pipeline': 'DESeq2',
        'parameters': {
            'fdr_threshold': 0.05,
            'lfc_threshold': 0.0
        },
        'results': {
            'n_significant': int(is_sig.sum()),
            'n_upregulated': int((direction == 'up').sum()),
            'n_downregulated': int((direction == 'down').sum())
        }
    }
    
    json_path = output_dir / 'de_summary.json'
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return csv_path


# =============================================================================
# DEResult Class Tests
# =============================================================================

class TestDEResultClass:
    """Test DEResult dataclass functionality."""
    
    def test_deresult_initialization(self, deseq2_results):
        """Test basic DEResult initialization."""
        df = deseq2_results.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2',
            parameters={'fdr_threshold': 0.05},
            metadata={'source': 'test'}
        )
        
        assert result.pipeline == 'DESeq2'
        assert result.parameters['fdr_threshold'] == 0.05
        assert result.metadata['source'] == 'test'
        assert result.results_df.index.name == 'gene_id'
    
    def test_deresult_properties(self, deseq2_results_with_significance):
        """Test DEResult computed properties."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2'
        )
        
        # Test basic properties
        assert result.n_genes == len(df)
        assert result.n_significant >= 0
        assert result.n_up >= 0
        assert result.n_down >= 0
        assert result.n_up + result.n_down <= result.n_significant
        
        # Test gene lists
        assert len(result.significant_genes) == result.n_significant
        assert len(result.upregulated_genes) == result.n_up
        assert len(result.downregulated_genes) == result.n_down
    
    def test_deresult_filter_by_threshold(self, deseq2_results_with_significance):
        """Test filtering by different thresholds."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2',
            parameters={'fdr_threshold': 0.05, 'lfc_threshold': 0.0}
        )
        
        original_sig = result.n_significant
        
        # Apply stricter threshold
        filtered = result.filter_by_threshold(fdr=0.01, lfc=1.0)
        
        # Should have fewer significant genes
        assert filtered.n_significant <= original_sig
        assert filtered.parameters['fdr_threshold'] == 0.01
        assert filtered.parameters['lfc_threshold'] == 1.0
    
    def test_deresult_get_top_genes(self, deseq2_results_with_significance):
        """Test getting top DE genes."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2'
        )
        
        # Get top 50 genes
        top_genes = result.get_top_genes(n=50)
        
        assert len(top_genes) == min(50, result.n_significant)
        
        # Should be sorted by adjusted p-value
        if len(top_genes) > 1:
            assert (top_genes['adjusted_p_value'].diff()[1:] >= 0).all()
    
    def test_deresult_get_genes_by_pattern(self, deseq2_results_with_significance):
        """Test gene filtering by regex pattern."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2'
        )
        
        # Get genes matching pattern
        pattern_genes = result.get_genes_by_pattern('ENSG00000000[0-9]{3}')
        
        assert all(pattern_genes.index.str.match('ENSG00000000[0-9]{3}'))
    
    def test_deresult_calculate_metrics(self, deseq2_results_with_significance, ground_truth_de):
        """Test performance metrics calculation."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2',
            parameters={'fdr_threshold': 0.05, 'lfc_threshold': 0.0}
        )
        
        # Basic metrics (no ground truth)
        metrics = result.calculate_metrics()
        
        assert 'n_genes' in metrics
        assert 'n_significant' in metrics
        assert 'n_up' in metrics
        assert 'n_down' in metrics
        
        # Metrics with ground truth
        metrics_gt = result.calculate_metrics(ground_truth=ground_truth_de)
        
        assert 'sensitivity' in metrics_gt
        assert 'specificity' in metrics_gt
        assert 'precision' in metrics_gt
        assert 'f1_score' in metrics_gt
        assert 'actual_fdr' in metrics_gt
        
        # Metrics should be between 0 and 1
        assert 0 <= metrics_gt['sensitivity'] <= 1
        assert 0 <= metrics_gt['specificity'] <= 1
        assert 0 <= metrics_gt['precision'] <= 1
        assert 0 <= metrics_gt['f1_score'] <= 1
    
    def test_deresult_compare_with(self, deseq2_results_with_significance, edger_results):
        """Test comparison between two DEResults."""
        # Create first DEResult
        df1 = deseq2_results_with_significance.set_index('gene_id')
        df1 = df1.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        result1 = DEResult(results_df=df1, pipeline='DESeq2')
        
        # Create second DEResult
        df2 = edger_results.set_index('gene_id')
        df2 = df2.rename(columns={
            'logFC': 'log2_fold_change',
            'PValue': 'p_value',
            'FDR': 'adjusted_p_value'
        })
        df2['is_significant'] = (df2['FDR'] < 0.05) & (df2['logFC'].abs() > 0)
        result2 = DEResult(results_df=df2, pipeline='edgeR')
        
        # Compare
        comparison = result1.compare_with(result2)
        
        assert 'overlap' in comparison
        assert 'only_in_first' in comparison
        assert 'only_in_second' in comparison
        assert 'correlation' in comparison
        
        # Overlap should be non-negative
        assert comparison['overlap'] >= 0
    
    def test_deresult_summary(self, deseq2_results_with_significance):
        """Test summary string generation."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(
            results_df=df,
            pipeline='DESeq2',
            parameters={'fdr_threshold': 0.05, 'lfc_threshold': 0.0}
        )
        
        summary = result.summary()
        
        assert isinstance(summary, str)
        assert 'RAPTOR' in summary
        assert 'DESeq2' in summary
        assert str(result.n_genes) in summary
        assert str(result.n_significant) in summary


# =============================================================================
# Column Mapping and Standardization Tests
# =============================================================================

class TestColumnMapping:
    """Test column mapping and standardization."""
    
    def test_detect_pipeline_deseq2(self, deseq2_results):
        """Test auto-detection of DESeq2 results."""
        detected = detect_pipeline(deseq2_results)
        assert detected == 'deseq2'
    
    def test_detect_pipeline_edger(self, edger_results):
        """Test auto-detection of edgeR results."""
        detected = detect_pipeline(edger_results)
        assert detected == 'edger'
    
    def test_detect_pipeline_limma(self, limma_results):
        """Test auto-detection of limma results."""
        detected = detect_pipeline(limma_results)
        assert detected == 'limma'
    
    def test_standardize_deseq2_columns(self, deseq2_results):
        """Test standardization of DESeq2 column names."""
        standardized = standardize_columns(deseq2_results, 'deseq2')
        
        # Check required columns exist
        for col in REQUIRED_COLUMNS:
            assert col in standardized.columns
        
        # Check index is gene_id
        assert standardized.index.name == 'gene_id'
        
        # Check column mapping worked
        assert 'log2_fold_change' in standardized.columns
        assert 'p_value' in standardized.columns
        assert 'adjusted_p_value' in standardized.columns
    
    def test_standardize_edger_columns(self, edger_results):
        """Test standardization of edgeR column names."""
        standardized = standardize_columns(edger_results, 'edger')
        
        for col in REQUIRED_COLUMNS:
            assert col in standardized.columns
        
        assert standardized.index.name == 'gene_id'
    
    def test_standardize_limma_columns(self, limma_results):
        """Test standardization of limma column names."""
        standardized = standardize_columns(limma_results, 'limma')
        
        for col in REQUIRED_COLUMNS:
            assert col in standardized.columns
        
        assert standardized.index.name == 'gene_id'
    
    def test_standardize_missing_gene_id(self, deseq2_results):
        """Test standardization with missing gene_id column."""
        df = deseq2_results.drop(columns=['gene_id'])
        
        # Should use first column or index
        standardized = standardize_columns(df, 'deseq2')
        assert standardized.index.name == 'gene_id'
    
    def test_standardize_custom_gene_column(self, deseq2_results):
        """Test standardization with custom gene ID column."""
        df = deseq2_results.rename(columns={'gene_id': 'ensembl_id'})
        
        standardized = standardize_columns(df, 'deseq2', gene_id_column='ensembl_id')
        assert standardized.index.name == 'gene_id'


# =============================================================================
# Significance Calculation Tests
# =============================================================================

class TestSignificanceCalculation:
    """Test significance and direction calculation."""
    
    def test_calculate_significance_default(self, deseq2_results):
        """Test significance calculation with default thresholds."""
        df = deseq2_results.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'padj': 'adjusted_p_value'
        })
        
        result = calculate_significance(df, fdr_threshold=0.05, lfc_threshold=0.0)
        
        assert 'is_significant' in result.columns
        assert 'direction' in result.columns
        
        # Check significance logic
        expected_sig = (result['adjusted_p_value'] < 0.05) & (result['log2_fold_change'].abs() > 0)
        assert (result['is_significant'] == expected_sig).all()
    
    def test_calculate_significance_strict(self, deseq2_results):
        """Test significance with strict thresholds."""
        df = deseq2_results.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'padj': 'adjusted_p_value'
        })
        
        result = calculate_significance(df, fdr_threshold=0.01, lfc_threshold=1.0)
        
        # Stricter thresholds should give fewer significant genes
        strict_sig = result['is_significant'].sum()
        
        result_loose = calculate_significance(df, fdr_threshold=0.05, lfc_threshold=0.0)
        loose_sig = result_loose['is_significant'].sum()
        
        assert strict_sig <= loose_sig
    
    def test_direction_calculation(self, deseq2_results):
        """Test direction assignment."""
        df = deseq2_results.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'padj': 'adjusted_p_value'
        })
        
        result = calculate_significance(df, fdr_threshold=0.05, lfc_threshold=0.0)
        
        # Check direction assignment
        sig_up = (result['is_significant']) & (result['log2_fold_change'] > 0)
        sig_down = (result['is_significant']) & (result['log2_fold_change'] < 0)
        
        assert (result.loc[sig_up, 'direction'] == 'up').all()
        assert (result.loc[sig_down, 'direction'] == 'down').all()
        assert (result.loc[~result['is_significant'], 'direction'] == 'unchanged').all()


# =============================================================================
# Import Function Tests
# =============================================================================

class TestImportFunction:
    """Test main import_de_results function."""
    
    def test_import_deseq2_basic(self, saved_deseq2_csv, output_dir):
        """Test basic DESeq2 import."""
        result = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'imported',
            pipeline='deseq2',
            fdr_threshold=0.05,
            lfc_threshold=0.0
        )
        
        assert isinstance(result, DEResult)
        assert result.pipeline == 'DESEQ2'
        assert result.n_genes > 0
        
        # Check outputs were created
        out_dir = output_dir / 'imported'
        assert (out_dir / 'de_standardized.csv').exists()
        assert (out_dir / 'de_significant.csv').exists()
        assert (out_dir / 'de_summary.json').exists()
        assert (out_dir / 'de_result.pkl').exists()
    
    def test_import_edger(self, saved_edger_csv, output_dir):
        """Test edgeR import."""
        result = import_de_results(
            de_file=saved_edger_csv,
            output_dir=output_dir / 'imported',
            pipeline='edger'
        )
        
        assert result.pipeline == 'EDGER'
        assert result.n_genes > 0
    
    def test_import_limma_tsv(self, saved_limma_tsv, output_dir):
        """Test limma import from TSV."""
        result = import_de_results(
            de_file=saved_limma_tsv,
            output_dir=output_dir / 'imported',
            pipeline='limma'
        )
        
        assert result.pipeline == 'LIMMA'
        assert result.n_genes > 0
    
    def test_import_auto_detect(self, saved_deseq2_csv, output_dir):
        """Test pipeline auto-detection."""
        result = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'imported',
            pipeline='auto'
        )
        
        assert result.pipeline in ['DESEQ2', 'EDGER', 'LIMMA', 'WILCOXON']
    
    def test_import_with_thresholds(self, saved_deseq2_csv, output_dir):
        """Test import with custom thresholds."""
        result = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'imported',
            pipeline='deseq2',
            fdr_threshold=0.01,
            lfc_threshold=1.0
        )
        
        assert result.parameters['fdr_threshold'] == 0.01
        assert result.parameters['lfc_threshold'] == 1.0
    
    def test_import_module6_output(self, module6_deseq2_output, output_dir):
        """Test importing actual Module 6 R script output."""
        result = import_de_results(
            de_file=module6_deseq2_output,
            output_dir=output_dir / 'imported',
            pipeline='auto'  # Should auto-detect DESeq2
        )
        
        assert result.pipeline == 'DESEQ2'
        assert result.n_genes > 0
        assert result.n_significant >= 0
        
        # Module 6 output already has significance columns
        # These should be preserved or recalculated consistently
        assert 'is_significant' in result.results_df.columns
        assert 'direction' in result.results_df.columns
    
    def test_import_creates_summary_json(self, saved_deseq2_csv, output_dir):
        """Test that import creates proper summary JSON."""
        result = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'imported',
            pipeline='deseq2'
        )
        
        summary_file = output_dir / 'imported' / 'de_summary.json'
        assert summary_file.exists()
        
        with open(summary_file) as f:
            summary = json.load(f)
        
        assert 'n_genes' in summary
        assert 'n_significant' in summary
        assert 'n_up' in summary
        assert 'n_down' in summary


# =============================================================================
# Serialization Tests
# =============================================================================

class TestSerialization:
    """Test DEResult save/load functionality."""
    
    def test_save_load_pickle(self, deseq2_results_with_significance, output_dir):
        """Test saving and loading as pickle."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        original = DEResult(
            results_df=df,
            pipeline='DESeq2',
            parameters={'fdr': 0.05},
            metadata={'test': 'value'}
        )
        
        # Save
        pkl_path = output_dir / 'result.pkl'
        original.save(pkl_path)
        assert pkl_path.exists()
        
        # Load
        loaded = DEResult.load(pkl_path)
        
        assert loaded.pipeline == original.pipeline
        assert loaded.parameters == original.parameters
        assert loaded.metadata == original.metadata
        assert loaded.n_genes == original.n_genes
        assert loaded.n_significant == original.n_significant
    
    def test_save_load_json(self, deseq2_results_with_significance, output_dir):
        """Test saving and loading as JSON."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        original = DEResult(
            results_df=df,
            pipeline='DESeq2',
            parameters={'fdr': 0.05}
        )
        
        # Save
        json_path = output_dir / 'result.json'
        original.save(json_path)
        assert json_path.exists()
        
        # Load
        loaded = DEResult.load(json_path)
        
        assert loaded.pipeline == original.pipeline
        assert loaded.n_genes == original.n_genes
    
    def test_save_invalid_extension(self, deseq2_results_with_significance, output_dir):
        """Test save with invalid extension."""
        df = deseq2_results_with_significance.set_index('gene_id')
        df = df.rename(columns={
            'log2FoldChange': 'log2_fold_change',
            'pvalue': 'p_value',
            'padj': 'adjusted_p_value'
        })
        
        result = DEResult(results_df=df, pipeline='DESeq2')
        
        with pytest.raises(ValidationError):
            result.save(output_dir / 'result.txt')


# =============================================================================
# Integration Tests
# =============================================================================

@pytest.mark.integration
class TestM6toM7Integration:
    """Test integration between Module 6 (R) and Module 7 (Python)."""
    
    def test_full_workflow(self, module6_deseq2_output, output_dir):
        """Test complete M6→M7 workflow."""
        # Import Module 6 output
        de_result = import_de_results(
            de_file=module6_deseq2_output,
            output_dir=output_dir / 'imported',
            pipeline='auto'
        )
        
        # Verify proper import
        assert de_result.pipeline == 'DESEQ2'
        assert de_result.n_genes > 0
        
        # Verify can be saved for M8-M10
        pkl_path = output_dir / 'for_m8.pkl'
        de_result.save(pkl_path)
        
        # Verify can be reloaded
        reloaded = DEResult.load(pkl_path)
        assert reloaded.n_genes == de_result.n_genes
    
    def test_multiple_methods(self, saved_deseq2_csv, saved_edger_csv, output_dir):
        """Test importing results from multiple methods."""
        # Import DESeq2
        deseq2 = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'deseq2',
            pipeline='deseq2'
        )
        
        # Import edgeR
        edger = import_de_results(
            de_file=saved_edger_csv,
            output_dir=output_dir / 'edger',
            pipeline='edger'
        )
        
        # Compare
        comparison = deseq2.compare_with(edger)
        
        assert 'overlap' in comparison
        assert comparison['overlap'] >= 0


# =============================================================================
# Validation Tests
# =============================================================================

class TestValidation:
    """Test input validation."""
    
    def test_invalid_file_path(self, output_dir):
        """Test with non-existent file."""
        with pytest.raises(Exception):  # Could be FileNotFoundError or ValidationError
            import_de_results(
                de_file=output_dir / 'nonexistent.csv',
                output_dir=output_dir / 'out',
                pipeline='deseq2'
            )
    
    def test_invalid_fdr_threshold(self, saved_deseq2_csv, output_dir):
        """Test with invalid FDR threshold."""
        with pytest.raises(ValidationError):
            import_de_results(
                de_file=saved_deseq2_csv,
                output_dir=output_dir / 'out',
                pipeline='deseq2',
                fdr_threshold=1.5  # Invalid: > 1
            )
    
    def test_missing_required_columns(self, output_dir):
        """Test with missing required columns."""
        # Create DataFrame missing required column
        df = pd.DataFrame({
            'gene_id': ['Gene1', 'Gene2'],
            'log2FoldChange': [1.0, -1.0]
            # Missing pvalue and padj
        })
        
        csv_path = output_dir / 'incomplete.csv'
        df.to_csv(csv_path, index=False)
        
        with pytest.raises(Exception):
            import_de_results(
                de_file=csv_path,
                output_dir=output_dir / 'out',
                pipeline='deseq2'
            )


# =============================================================================
# CLI Compatibility Tests
# =============================================================================

@pytest.mark.cli
class TestCLICompatibility:
    """Test that import works as expected by CLI."""
    
    def test_cli_parameter_format(self, saved_deseq2_csv, output_dir):
        """Test with parameters as CLI would call them."""
        # This mimics how the fixed CLI calls the function
        result = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'cli_test',
            pipeline='deseq2',  # CLI passes this as 'method' → 'pipeline'
            fdr_threshold=0.05,
            lfc_threshold=0.0
        )
        
        assert result.pipeline == 'DESEQ2'
        
        # Check CLI can access attributes correctly
        assert hasattr(result, 'pipeline')  # NOT 'method'
        assert hasattr(result, 'n_up')  # NOT 'n_upregulated'
        assert hasattr(result, 'n_down')  # NOT 'n_downregulated'
    
    def test_output_files_for_cli(self, saved_deseq2_csv, output_dir):
        """Test that all expected output files are created for CLI."""
        result = import_de_results(
            de_file=saved_deseq2_csv,
            output_dir=output_dir / 'cli_out',
            pipeline='deseq2'
        )
        
        out_dir = output_dir / 'cli_out'
        
        # These files should exist (CLI checks for them)
        assert (out_dir / 'de_standardized.csv').exists()
        assert (out_dir / 'de_significant.csv').exists()
        assert (out_dir / 'de_summary.json').exists()
        assert (out_dir / 'de_result.pkl').exists()


# =============================================================================
# Performance Tests
# =============================================================================

@pytest.mark.slow
class TestPerformance:
    """Test performance with large datasets."""
    
    def test_large_dataset(self, output_dir):
        """Test with large number of genes."""
        np.random.seed(42)
        n_genes = 50000  # Realistic human transcriptome size
        
        df = pd.DataFrame({
            'gene_id': [f'ENSG{i+1:011d}' for i in range(n_genes)],
            'baseMean': np.random.gamma(2, 100, n_genes),
            'log2FoldChange': np.random.normal(0, 1.5, n_genes),
            'lfcSE': np.random.uniform(0.1, 0.5, n_genes),
            'stat': np.random.normal(0, 3, n_genes),
            'pvalue': np.random.beta(1, 10, n_genes),
            'padj': np.random.beta(1, 20, n_genes)
        })
        
        csv_path = output_dir / 'large_dataset.csv'
        df.to_csv(csv_path, index=False)
        
        # Should handle large dataset without issues
        result = import_de_results(
            de_file=csv_path,
            output_dir=output_dir / 'large_out',
            pipeline='deseq2'
        )
        
        assert result.n_genes == n_genes


# =============================================================================
# Helper Functions
# =============================================================================

def assert_deresult_valid(result: DEResult):
    """Assert that a DEResult object is valid."""
    assert isinstance(result, DEResult)
    assert isinstance(result.results_df, pd.DataFrame)
    assert result.pipeline in ['DESEQ2', 'EDGER', 'LIMMA', 'WILCOXON']
    assert result.n_genes > 0
    assert result.n_significant >= 0
    assert result.n_up >= 0
    assert result.n_down >= 0
    assert result.n_up + result.n_down <= result.n_significant
    
    # Check required columns
    for col in REQUIRED_COLUMNS:
        assert col in result.results_df.columns
    
    # Check index
    assert result.results_df.index.name == 'gene_id'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
