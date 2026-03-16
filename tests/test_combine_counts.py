"""
RAPTOR v2.2.0 - Tests for combine_counts.py

Tests for count matrix combination and aggregation.

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
    from raptor.pipelines.quick_salmon.combine_counts import (
        load_tx2gene,
        load_abundance,
        aggregate_to_gene,
        combine_samples,
        save_matrices,
        OUTPUT_COUNTS_FILE,
        OUTPUT_TPM_FILE
    )
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False


# =============================================================================
# Tests for tx2gene Loading
# =============================================================================

class TestLoadTx2Gene:
    """Tests for tx2gene loading functions."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_load_tsv_tx2gene(self, tx2gene_path):
        """Test loading tx2gene from TSV file."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        assert isinstance(tx2gene, dict)
        assert len(tx2gene) > 0
        
        # Check mapping format
        for tx_id, gene_id in list(tx2gene.items())[:5]:
            assert tx_id.startswith('ENST')
            assert gene_id.startswith('ENSG')
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_load_gtf_tx2gene(self, gtf_tx2gene_path):
        """Test loading tx2gene from GTF file."""
        tx2gene = load_tx2gene(str(gtf_tx2gene_path))
        
        assert isinstance(tx2gene, dict)
        assert len(tx2gene) > 0
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_tx2gene_mapping_consistency(self, tx2gene_path):
        """Test that transcript-gene mapping is consistent."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        # Same transcript should always map to same gene
        for tx_id, gene_id in tx2gene.items():
            assert tx2gene[tx_id] == gene_id


# =============================================================================
# Tests for Abundance Loading
# =============================================================================

class TestLoadAbundance:
    """Tests for abundance file loading."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_load_salmon_abundance(self, mock_salmon_quant_dir):
        """Test loading Salmon quant.sf files."""
        sample_dirs = list(Path(mock_salmon_quant_dir).iterdir())
        sample_dir = sample_dirs[0]
        
        df = load_abundance(str(sample_dir), tool='salmon')
        
        assert isinstance(df, pd.DataFrame)
        assert 'Name' in df.columns
        assert 'NumReads' in df.columns
        assert 'TPM' in df.columns
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_load_kallisto_abundance(self, mock_kallisto_quant_dir):
        """Test loading Kallisto abundance.tsv files."""
        sample_dirs = list(Path(mock_kallisto_quant_dir).iterdir())
        sample_dir = sample_dirs[0]
        
        df = load_abundance(str(sample_dir), tool='kallisto')
        
        assert isinstance(df, pd.DataFrame)
        assert 'Name' in df.columns  # Renamed from target_id
        assert 'NumReads' in df.columns  # Renamed from est_counts


# =============================================================================
# Tests for Gene Aggregation
# =============================================================================

class TestAggregateToGene:
    """Tests for transcript to gene aggregation."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_aggregate_to_gene(self, mock_salmon_quant_dir, tx2gene_path):
        """Test aggregation of transcript counts to gene level."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        sample_dirs = list(Path(mock_salmon_quant_dir).iterdir())
        sample_dir = sample_dirs[0]
        abundance_df = load_abundance(str(sample_dir), tool='salmon')
        
        gene_counts = aggregate_to_gene(abundance_df, tx2gene, 'NumReads')
        
        assert isinstance(gene_counts, pd.Series)
        assert len(gene_counts) > 0
        assert gene_counts.index[0].startswith('ENSG')
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_aggregate_sums_correctly(self, tx2gene_path):
        """Test that aggregation sums counts correctly."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        # Create test abundance data
        abundance_df = pd.DataFrame({
            'Name': list(tx2gene.keys())[:100],
            'NumReads': [10] * 100,
            'TPM': [100] * 100
        })
        
        gene_counts = aggregate_to_gene(abundance_df, tx2gene, 'NumReads')
        
        # Total should be preserved
        assert gene_counts.sum() == 1000  # 100 transcripts × 10 reads


# =============================================================================
# Tests for Sample Combination
# =============================================================================

class TestCombineSamples:
    """Tests for combining multiple samples."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_combine_salmon_samples(self, mock_salmon_quant_dir, tx2gene_path):
        """Test combining Salmon samples."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        counts_df, tpm_df, length_df = combine_samples(
            str(mock_salmon_quant_dir),
            tx2gene,
            tool='salmon'
        )
        
        assert isinstance(counts_df, pd.DataFrame)
        assert isinstance(tpm_df, pd.DataFrame)
        assert len(counts_df.columns) == 6  # 6 samples
        assert counts_df.shape[0] == tpm_df.shape[0]
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_combine_kallisto_samples(self, mock_kallisto_quant_dir, tx2gene_path):
        """Test combining Kallisto samples."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        counts_df, tpm_df, length_df = combine_samples(
            str(mock_kallisto_quant_dir),
            tx2gene,
            tool='kallisto'
        )
        
        assert isinstance(counts_df, pd.DataFrame)
        assert len(counts_df.columns) == 6
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_counts_are_integers(self, mock_salmon_quant_dir, tx2gene_path):
        """Test that counts are rounded to integers."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        counts_df, _, _ = combine_samples(
            str(mock_salmon_quant_dir),
            tx2gene,
            tool='salmon'
        )
        
        # All values should be integers
        assert counts_df.dtypes.apply(lambda x: np.issubdtype(x, np.integer)).all()


# =============================================================================
# Tests for Matrix Saving
# =============================================================================

class TestSaveMatrices:
    """Tests for saving count matrices."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_save_matrices(self, sample_count_matrix, sample_tpm_matrix, output_dir):
        """Test saving count and TPM matrices."""
        length_df = pd.DataFrame()  # Empty for this test
        
        save_matrices(sample_count_matrix, sample_tpm_matrix, length_df, str(output_dir))
        
        # Check files exist
        assert (output_dir / OUTPUT_COUNTS_FILE).exists()
        assert (output_dir / OUTPUT_TPM_FILE).exists()
        assert (output_dir / 'count_summary.txt').exists()
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_saved_matrix_readable(self, sample_count_matrix, sample_tpm_matrix, output_dir):
        """Test that saved matrices can be read back."""
        length_df = pd.DataFrame()
        
        save_matrices(sample_count_matrix, sample_tpm_matrix, length_df, str(output_dir))
        
        # Read back
        loaded_counts = pd.read_csv(output_dir / OUTPUT_COUNTS_FILE, index_col=0)
        loaded_tpm = pd.read_csv(output_dir / OUTPUT_TPM_FILE, index_col=0)
        
        # Compare
        pd.testing.assert_frame_equal(sample_count_matrix, loaded_counts)


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_unmapped_transcripts_warning(self, mock_salmon_quant_dir, capsys):
        """Test warning for unmapped transcripts."""
        # Create incomplete tx2gene
        incomplete_tx2gene = {
            'ENST00000000001': 'ENSG00000000001',
            'ENST00000000002': 'ENSG00000000001'
        }
        
        counts_df, _, _ = combine_samples(
            str(mock_salmon_quant_dir),
            incomplete_tx2gene,
            tool='salmon'
        )
        
        # Should still work, just with fewer genes
        assert isinstance(counts_df, pd.DataFrame)
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_empty_quant_directory(self, test_dir, tx2gene_path):
        """Test handling of empty quant directory."""
        empty_dir = test_dir / "empty_quant"
        empty_dir.mkdir(exist_ok=True)
        
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        with pytest.raises(ValueError):
            combine_samples(str(empty_dir), tx2gene, tool='salmon')


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for full workflow."""
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_full_workflow_salmon(self, mock_salmon_quant_dir, tx2gene_path, output_dir):
        """Test complete workflow for Salmon."""
        # Load tx2gene
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        # Combine samples
        counts_df, tpm_df, length_df = combine_samples(
            str(mock_salmon_quant_dir),
            tx2gene,
            tool='salmon'
        )
        
        # Save matrices
        save_matrices(counts_df, tpm_df, length_df, str(output_dir))
        
        # Verify output
        assert (output_dir / OUTPUT_COUNTS_FILE).exists()
        
        loaded = pd.read_csv(output_dir / OUTPUT_COUNTS_FILE, index_col=0)
        assert loaded.shape[0] > 0
        assert loaded.shape[1] == 6
    
    @pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Module not available")
    def test_full_workflow_kallisto(self, mock_kallisto_quant_dir, tx2gene_path, output_dir):
        """Test complete workflow for Kallisto."""
        tx2gene = load_tx2gene(str(tx2gene_path))
        
        counts_df, tpm_df, length_df = combine_samples(
            str(mock_kallisto_quant_dir),
            tx2gene,
            tool='kallisto'
        )
        
        save_matrices(counts_df, tpm_df, length_df, str(output_dir))
        
        assert (output_dir / OUTPUT_COUNTS_FILE).exists()
