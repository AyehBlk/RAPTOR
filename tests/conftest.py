"""
RAPTOR v2.2.0 - Test Fixtures

Shared pytest fixtures for all RAPTOR tests.

Supports:
- Module 1: Quick Quantification
- Module 2: Quality Assessment
- Module 3: Data Profiling
- Module 4: Pipeline Recommendation
- (Future modules)

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import json
import gzip


# =============================================================================
# Constants
# =============================================================================

# Module 1 outputs
DEFAULT_QUICK_COUNTS_DIR = "results/quick_counts"
OUTPUT_COUNTS_FILE = "quick_gene_counts.csv"
OUTPUT_TPM_FILE = "quick_tpm.csv"
OUTPUT_SAMPLE_INFO = "sample_info.csv"

# Module 2 outputs
DEFAULT_QC_DIR = "results/qc"


# =============================================================================
# Base Fixtures
# =============================================================================

@pytest.fixture(scope="session")
def test_dir():
    """Create a temporary directory for test files."""
    tmp_dir = tempfile.mkdtemp(prefix="raptor_test_")
    yield Path(tmp_dir)
    # Cleanup after all tests
    shutil.rmtree(tmp_dir, ignore_errors=True)


@pytest.fixture(scope="function")
def output_dir(test_dir):
    """Create a clean output directory for each test."""
    out_dir = test_dir / "output" / f"test_{np.random.randint(10000)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    yield out_dir
    # Cleanup after test
    shutil.rmtree(out_dir, ignore_errors=True)


# =============================================================================
# Count Matrix Fixtures
# =============================================================================

@pytest.fixture
def sample_count_matrix():
    """Generate a sample count matrix for testing (8 samples)."""
    np.random.seed(42)
    
    n_genes = 1000
    n_samples = 8
    
    # Generate realistic counts
    base_expr = np.random.gamma(shape=2, scale=100, size=n_genes)
    counts = np.zeros((n_genes, n_samples))
    
    for i in range(n_samples):
        size_param = 10
        counts[:, i] = np.random.negative_binomial(
            size_param, 
            size_param / (size_param + base_expr)
        )
    
    gene_names = [f'ENSG{i+1:011d}' for i in range(n_genes)]
    sample_names = ['Control_1', 'Control_2', 'Control_3', 'Control_4',
                    'Treatment_1', 'Treatment_2', 'Treatment_3', 'Treatment_4']
    
    return pd.DataFrame(
        counts.astype(int),
        index=gene_names,
        columns=sample_names
    )


@pytest.fixture
def sample_count_matrix_6(sample_count_matrix):
    """6-sample count matrix (for M1 compatibility)."""
    return sample_count_matrix.iloc[:, :6].copy()


@pytest.fixture
def sample_count_matrix_with_outlier(sample_count_matrix):
    """Generate count matrix with an artificial outlier."""
    counts = sample_count_matrix.copy()
    # Make one sample an outlier (very different expression)
    counts['Treatment_4'] = counts['Treatment_4'] * 0.1
    return counts


@pytest.fixture
def small_count_matrix():
    """Generate a small count matrix for fast testing."""
    np.random.seed(42)
    
    n_genes = 100
    n_samples = 4
    
    counts = np.random.negative_binomial(10, 0.1, (n_genes, n_samples))
    
    return pd.DataFrame(
        counts,
        index=[f'Gene{i+1}' for i in range(n_genes)],
        columns=['Sample1', 'Sample2', 'Sample3', 'Sample4']
    )


@pytest.fixture
def high_zero_count_matrix():
    """Generate count matrix with high zero inflation."""
    np.random.seed(42)
    
    n_genes = 500
    n_samples = 6
    
    # Generate sparse counts (70% zeros)
    counts = np.random.negative_binomial(5, 0.5, (n_genes, n_samples))
    zero_mask = np.random.random((n_genes, n_samples)) < 0.7
    counts[zero_mask] = 0
    
    return pd.DataFrame(
        counts,
        index=[f'Gene{i+1}' for i in range(n_genes)],
        columns=[f'Sample{i+1}' for i in range(n_samples)]
    )


@pytest.fixture
def varying_library_size_matrix():
    """Generate count matrix with varying library sizes."""
    np.random.seed(42)
    
    n_genes = 500
    n_samples = 6
    
    # Base expression
    base_expr = np.random.gamma(2, 100, n_genes)
    
    # Varying size factors (10x range)
    size_factors = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    
    counts = np.zeros((n_genes, n_samples))
    for i, sf in enumerate(size_factors):
        counts[:, i] = np.random.negative_binomial(
            10, 10 / (10 + base_expr * sf), n_genes
        )
    
    return pd.DataFrame(
        counts.astype(int),
        index=[f'Gene{i+1}' for i in range(n_genes)],
        columns=[f'Sample{i+1}' for i in range(n_samples)]
    )


@pytest.fixture
def sample_tpm_matrix(sample_count_matrix):
    """Generate a sample TPM matrix from counts."""
    np.random.seed(42)
    gene_lengths = np.random.uniform(1000, 10000, size=len(sample_count_matrix))
    rpk = sample_count_matrix.div(gene_lengths / 1000, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    return tpm


# =============================================================================
# Metadata Fixtures
# =============================================================================

@pytest.fixture
def sample_metadata():
    """Generate sample metadata (8 samples)."""
    return pd.DataFrame({
        'sample_id': ['Control_1', 'Control_2', 'Control_3', 'Control_4',
                      'Treatment_1', 'Treatment_2', 'Treatment_3', 'Treatment_4'],
        'condition': ['Control'] * 4 + ['Treatment'] * 4,
        'batch': ['Batch1', 'Batch1', 'Batch2', 'Batch2',
                  'Batch1', 'Batch1', 'Batch2', 'Batch2']
    })


@pytest.fixture
def confounded_metadata():
    """Generate confounded metadata (batch = condition)."""
    return pd.DataFrame({
        'sample_id': ['Control_1', 'Control_2', 'Control_3', 'Control_4',
                      'Treatment_1', 'Treatment_2', 'Treatment_3', 'Treatment_4'],
        'condition': ['Control'] * 4 + ['Treatment'] * 4,
        'batch': ['Batch1'] * 4 + ['Batch2'] * 4  # Perfectly confounded!
    })


# =============================================================================
# Sample Sheet Fixtures
# =============================================================================

@pytest.fixture
def sample_sheet_path(test_dir):
    """Create a sample sheet CSV file (paired-end)."""
    sample_sheet = test_dir / "samples.csv"
    
    data = {
        'sample_id': ['Control_1', 'Control_2', 'Control_3',
                      'Treatment_1', 'Treatment_2', 'Treatment_3'],
        'condition': ['Control'] * 3 + ['Treatment'] * 3,
        'batch': ['Batch1'] * 6,
        'fastq_r1': [f'/fake/path/Sample{i}_R1.fastq.gz' for i in range(1, 7)],
        'fastq_r2': [f'/fake/path/Sample{i}_R2.fastq.gz' for i in range(1, 7)]
    }
    
    pd.DataFrame(data).to_csv(sample_sheet, index=False)
    return sample_sheet


@pytest.fixture
def single_end_sample_sheet_path(test_dir):
    """Create a single-end sample sheet CSV file."""
    sample_sheet = test_dir / "samples_single.csv"
    
    data = {
        'sample_id': ['Sample_1', 'Sample_2', 'Sample_3'],
        'condition': ['Control', 'Control', 'Treatment'],
        'batch': ['Batch1'] * 3,
        'fastq_r1': [f'/fake/path/Sample{i}.fastq.gz' for i in range(1, 4)],
        'fastq_r2': ['', '', '']
    }
    
    pd.DataFrame(data).to_csv(sample_sheet, index=False)
    return sample_sheet


# =============================================================================
# Transcript/Gene Mapping Fixtures
# =============================================================================

@pytest.fixture
def tx2gene_path(test_dir):
    """Create a tx2gene mapping file."""
    tx2gene_file = test_dir / "tx2gene.tsv"
    
    # Create mapping for 1000 transcripts to 500 genes
    data = {
        'transcript_id': [f'ENST{i+1:011d}' for i in range(1000)],
        'gene_id': [f'ENSG{(i//2)+1:011d}' for i in range(1000)]
    }
    
    pd.DataFrame(data).to_csv(tx2gene_file, sep='\t', index=False)
    return tx2gene_file


@pytest.fixture
def gtf_tx2gene_path(test_dir):
    """Create a GTF file for tx2gene extraction."""
    gtf_file = test_dir / "annotation.gtf"
    
    lines = ['##gtf-version 2.2']
    for i in range(100):
        gene_id = f'ENSG{i+1:011d}'
        for j in range(2):  # 2 transcripts per gene
            tx_id = f'ENST{i*2+j+1:011d}'
            attr = f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "Gene{i+1}";'
            lines.append(f'chr1\tENSEMBL\ttranscript\t{i*1000}\t{i*1000+999}\t.\t+\t.\t{attr}')
    
    with open(gtf_file, 'w') as f:
        f.write('\n'.join(lines))
    
    return gtf_file


# =============================================================================
# Mock FASTQ/Quant Fixtures
# =============================================================================

@pytest.fixture
def mock_fastq_dir(test_dir):
    """Create a directory with mock FASTQ files."""
    fastq_dir = test_dir / "fastq"
    fastq_dir.mkdir(exist_ok=True)
    
    # Create empty gzipped FASTQ files
    samples = ['Control_1', 'Control_2', 'Treatment_1', 'Treatment_2']
    
    for sample in samples:
        for read in ['R1', 'R2']:
            fq_path = fastq_dir / f"{sample}_{read}_001.fastq.gz"
            with gzip.open(fq_path, 'wt') as f:
                f.write(f"@{sample}_{read}_1\nACGT\n+\nIIII\n")
    
    return fastq_dir


@pytest.fixture
def mock_salmon_quant_dir(test_dir, sample_count_matrix_6):
    """Create mock Salmon quant output directory."""
    quant_dir = test_dir / "salmon_quant"
    
    for sample in sample_count_matrix_6.columns:
        sample_dir = quant_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        # Create quant.sf
        n_transcripts = 1000
        quant_data = pd.DataFrame({
            'Name': [f'ENST{i+1:011d}' for i in range(n_transcripts)],
            'Length': np.random.randint(500, 5000, n_transcripts),
            'EffectiveLength': np.random.randint(400, 4500, n_transcripts),
            'TPM': np.random.uniform(0, 100, n_transcripts),
            'NumReads': np.random.negative_binomial(10, 0.1, n_transcripts)
        })
        quant_data.to_csv(sample_dir / 'quant.sf', sep='\t', index=False)
        
        # Create aux_info/meta_info.json
        aux_dir = sample_dir / 'aux_info'
        aux_dir.mkdir(exist_ok=True)
        
        meta_info = {
            'salmon_version': '1.9.0',
            'num_processed': 25000000,
            'num_mapped': 23500000,
            'percent_mapped': 94.0
        }
        with open(aux_dir / 'meta_info.json', 'w') as f:
            json.dump(meta_info, f)
    
    return quant_dir


@pytest.fixture
def mock_kallisto_quant_dir(test_dir, sample_count_matrix_6):
    """Create mock Kallisto quant output directory."""
    quant_dir = test_dir / "kallisto_quant"
    
    for sample in sample_count_matrix_6.columns:
        sample_dir = quant_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        # Create abundance.tsv
        n_transcripts = 1000
        abundance_data = pd.DataFrame({
            'target_id': [f'ENST{i+1:011d}' for i in range(n_transcripts)],
            'length': np.random.randint(500, 5000, n_transcripts),
            'eff_length': np.random.randint(400, 4500, n_transcripts),
            'est_counts': np.random.negative_binomial(10, 0.1, n_transcripts),
            'tpm': np.random.uniform(0, 100, n_transcripts)
        })
        abundance_data.to_csv(sample_dir / 'abundance.tsv', sep='\t', index=False)
        
        # Create run_info.json
        run_info = {
            'n_processed': 25000000,
            'n_pseudoaligned': 23500000,
            'n_unique': 22000000,
            'p_pseudoaligned': 94.0,
            'p_unique': 88.0,
            'kallisto_version': '0.48.0'
        }
        with open(sample_dir / 'run_info.json', 'w') as f:
            json.dump(run_info, f)
    
    return quant_dir


# =============================================================================
# Saved File Fixtures
# =============================================================================

@pytest.fixture
def saved_count_matrix(output_dir, sample_count_matrix):
    """Save count matrix to file and return path."""
    counts_path = output_dir / OUTPUT_COUNTS_FILE
    sample_count_matrix.to_csv(counts_path)
    return counts_path


@pytest.fixture
def saved_tpm_matrix(output_dir, sample_tpm_matrix):
    """Save TPM matrix to file and return path."""
    tpm_path = output_dir / OUTPUT_TPM_FILE
    sample_tpm_matrix.to_csv(tpm_path)
    return tpm_path


@pytest.fixture
def saved_metadata(output_dir, sample_metadata):
    """Save metadata to file and return path."""
    meta_path = output_dir / "metadata.csv"
    sample_metadata.to_csv(meta_path, index=False)
    return meta_path


# =============================================================================
# Helper Functions - Module 1 (Quick Count)
# =============================================================================

def assert_count_matrix_valid(df: pd.DataFrame):
    """Assert that a count matrix is valid."""
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert len(df.columns) > 0
    assert df.values.min() >= 0  # No negative counts
    assert df.dtypes.apply(lambda x: np.issubdtype(x, np.integer)).all()


def assert_tpm_matrix_valid(df: pd.DataFrame):
    """Assert that a TPM matrix is valid."""
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert len(df.columns) > 0
    assert df.values.min() >= 0  # No negative TPM
    # TPM should sum to ~1M per sample
    for col in df.columns:
        assert 999000 < df[col].sum() < 1001000


def assert_output_files_exist(output_dir: Path):
    """Assert that all expected M1 output files exist."""
    assert (output_dir / OUTPUT_COUNTS_FILE).exists()
    assert (output_dir / OUTPUT_TPM_FILE).exists()
    assert (output_dir / OUTPUT_SAMPLE_INFO).exists()


# =============================================================================
# Helper Functions - Module 2 (Quality Assessment)
# =============================================================================

def assert_quality_report_valid(report: dict):
    """Assert that a quality report has expected structure."""
    assert 'overall' in report
    assert 'score' in report['overall']
    assert 'status' in report['overall']
    assert 0 <= report['overall']['score'] <= 100


def assert_outlier_result_valid(result):
    """Assert that outlier result has expected attributes."""
    assert hasattr(result, 'outlier_samples')
    assert hasattr(result, 'n_outliers')
    assert hasattr(result, 'method_results')
    assert hasattr(result, 'consensus_threshold')
    assert result.n_outliers >= 0
    assert result.n_outliers == len(result.outlier_samples)


# =============================================================================
# Module 6b: Data Acquisition Fixtures
# =============================================================================

try:
    from raptor.external_modules.acquisition import (
        AcquiredDataset,
        PooledDataset,
        CacheManager,
        DataCatalog,
        PoolingEngine,
    )
    _ACQUISITION_AVAILABLE = True
except ImportError:
    _ACQUISITION_AVAILABLE = False


@pytest.fixture
def acquisition_count_matrix():
    """NB-distributed count matrix for acquisition tests (500 genes x 8 samples)."""
    np.random.seed(42)
    counts = np.random.negative_binomial(n=5, p=0.3, size=(500, 8))
    df = pd.DataFrame(
        counts,
        index=[f'GENE_{i:04d}' for i in range(500)],
        columns=[f'Sample_{i}' for i in range(8)],
    )
    df.index.name = 'gene_id'
    return df


@pytest.fixture
def acquisition_metadata():
    """Sample metadata for acquisition tests."""
    meta = pd.DataFrame({
        'condition': ['control'] * 4 + ['treatment'] * 4,
        'batch': ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B'],
        'tissue': ['liver'] * 8,
    }, index=[f'Sample_{i}' for i in range(8)])
    meta.index.name = 'sample_id'
    return meta


@pytest.fixture
def geo_source_info():
    """Source info for a mock GEO dataset."""
    return {
        'repository': 'GEO',
        'accession': 'GSE99999',
        'organism': 'Homo sapiens',
        'platform': 'GPL16791',
        'data_type': 'raw_counts',
        'description': 'Test RNA-seq dataset for liver cancer',
    }


@pytest.fixture
def acquired_dataset(acquisition_count_matrix, acquisition_metadata, geo_source_info):
    """A complete AcquiredDataset for testing."""
    if not _ACQUISITION_AVAILABLE:
        pytest.skip("Acquisition subpackage not available")
    return AcquiredDataset(
        counts_df=acquisition_count_matrix,
        metadata=acquisition_metadata,
        source_info=geo_source_info,
        gene_id_type='symbol',
    )


@pytest.fixture
def acquired_dataset_2():
    """A second dataset with partial gene overlap for pooling tests."""
    if not _ACQUISITION_AVAILABLE:
        pytest.skip("Acquisition subpackage not available")
    np.random.seed(123)
    counts = pd.DataFrame(
        np.random.negative_binomial(n=4, p=0.25, size=(500, 6)),
        index=[f'GENE_{i:04d}' for i in range(200, 700)],
        columns=[f'SampleB_{i}' for i in range(6)],
    )
    counts.index.name = 'gene_id'
    meta = pd.DataFrame({
        'condition': ['control'] * 3 + ['treatment'] * 3,
        'batch': ['X', 'X', 'Y', 'X', 'Y', 'Y'],
    }, index=[f'SampleB_{i}' for i in range(6)])
    meta.index.name = 'sample_id'
    return AcquiredDataset(
        counts_df=counts,
        metadata=meta,
        source_info={
            'repository': 'GEO',
            'accession': 'GSE88888',
            'organism': 'Homo sapiens',
            'data_type': 'raw_counts',
        },
        gene_id_type='symbol',
    )


@pytest.fixture
def acquisition_cache(tmp_path):
    """CacheManager with a temporary directory."""
    if not _ACQUISITION_AVAILABLE:
        pytest.skip("Acquisition subpackage not available")
    return CacheManager(cache_dir=tmp_path / 'raptor_test_cache', enabled=True)


@pytest.fixture
def acquisition_catalog(tmp_path):
    """DataCatalog with a temporary directory."""
    if not _ACQUISITION_AVAILABLE:
        pytest.skip("Acquisition subpackage not available")
    return DataCatalog(tmp_path / 'raptor_test_catalog')


@pytest.fixture
def pooling_engine():
    """PoolingEngine for testing."""
    if not _ACQUISITION_AVAILABLE:
        pytest.skip("Acquisition subpackage not available")
    return PoolingEngine(target_gene_id='symbol', species='Homo sapiens')


# =============================================================================
# Helper Functions - Module 6b (Acquisition)
# =============================================================================

def assert_acquired_dataset_valid(ds):
    """Assert that an AcquiredDataset has expected structure."""
    assert ds.n_genes > 0
    assert ds.n_samples > 0
    assert ds.counts_df.index.name == 'gene_id'
    assert ds.source_info.get('repository') is not None
    report = ds.validate_integrity()
    assert report['valid'] is True


def assert_pooled_dataset_valid(pool):
    """Assert that a PooledDataset has expected structure."""
    assert pool.n_genes > 0
    assert pool.n_samples > 0
    assert pool.n_studies >= 2
    assert len(pool.studies) == pool.n_studies
    assert 'study' in pool.metadata.columns
