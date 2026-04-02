"""
RAPTOR v2.2.2 - Tests for Module 6b: Data Acquisition

Comprehensive test suite for the acquisition subpackage.
All tests run offline (no network) using mock data and mock API responses.

Test coverage:
    - datasets.py: AcquiredDataset, PooledDataset (creation, validation, save/load, export)
    - cache.py: CacheManager (save, load, list, delete, hybrid mode)
    - catalog.py: DataCatalog (register, search, filter, persistence)
    - base.py: BaseConnector, SearchResult (via MockConnector)
    - gene_mapping.py: GeneIDMapper (detection, same-type convert)
    - geo.py: GEOConnector (constructor, static methods)
    - tcga.py: TCGAConnector (constructor, project list)
    - arrayexpress.py: ArrayExpConnector (constructor)
    - sra.py: SRAConnector (constructor, GSM extraction, GSE lookup, download_api GSE population)
    - pooling.py: PoolingEngine (merge, batch correction, conflict resolution)

Run:
    pytest tests/test_acquisition.py -v
    pytest tests/test_acquisition.py -v --tb=short 2>&1 | Tee-Object test_results\\acquisition.txt

Author: Ayeh Bolouki
Version: 2.2.2
"""

import json
import pickle
import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def sample_counts():
    """Create a realistic simulated count matrix (NB-distributed)."""
    np.random.seed(42)
    n_genes = 500
    n_samples = 8

    # Negative binomial counts (realistic RNA-seq)
    counts = np.random.negative_binomial(n=5, p=0.3, size=(n_genes, n_samples))

    df = pd.DataFrame(
        counts,
        index=[f'GENE_{i:04d}' for i in range(n_genes)],
        columns=[f'Sample_{i}' for i in range(n_samples)],
    )
    df.index.name = 'gene_id'
    return df


@pytest.fixture
def sample_metadata():
    """Sample metadata matching sample_counts."""
    meta = pd.DataFrame({
        'condition': ['control'] * 4 + ['treatment'] * 4,
        'batch': ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B'],
        'tissue': ['liver'] * 8,
    }, index=[f'Sample_{i}' for i in range(8)])
    meta.index.name = 'sample_id'
    return meta


@pytest.fixture
def source_info_geo():
    """Source info for a GEO dataset."""
    return {
        'repository': 'GEO',
        'accession': 'GSE99999',
        'organism': 'Homo sapiens',
        'platform': 'GPL16791',
        'data_type': 'raw_counts',
        'description': 'Test RNA-seq dataset for liver cancer',
    }


@pytest.fixture
def source_info_tcga():
    """Source info for a TCGA dataset."""
    return {
        'repository': 'TCGA',
        'accession': 'TCGA-LIHC',
        'organism': 'Homo sapiens',
        'platform': 'Illumina HiSeq',
        'data_type': 'raw_counts',
        'description': 'Liver hepatocellular carcinoma',
    }


@pytest.fixture
def acquired_dataset(sample_counts, sample_metadata, source_info_geo):
    """A complete AcquiredDataset."""
    from raptor.external_modules.acquisition import AcquiredDataset
    return AcquiredDataset(
        counts_df=sample_counts,
        metadata=sample_metadata,
        source_info=source_info_geo,
        gene_id_type='symbol',
    )


@pytest.fixture
def acquired_dataset_2():
    """A second dataset with partial gene overlap for pooling tests."""
    from raptor.external_modules.acquisition import AcquiredDataset
    np.random.seed(123)

    # Genes 200-699: 300 overlap with dataset 1 (genes 200-499)
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


# =============================================================================
# TEST: AcquiredDataset
# =============================================================================

class TestAcquiredDataset:
    """Tests for AcquiredDataset data container."""

    def test_creation(self, acquired_dataset):
        """Basic creation and properties."""
        ds = acquired_dataset
        assert ds.n_genes == 500
        assert ds.n_samples == 8
        assert ds.accession == 'GSE99999'
        assert ds.repository == 'GEO'
        assert ds.organism == 'Homo sapiens'
        assert ds.data_type == 'raw_counts'
        assert ds.gene_id_type == 'symbol'

    def test_gene_and_sample_ids(self, acquired_dataset):
        """Gene and sample ID lists."""
        ds = acquired_dataset
        assert len(ds.gene_ids) == 500
        assert ds.gene_ids[0] == 'GENE_0000'
        assert len(ds.sample_ids) == 8
        assert ds.sample_ids[0] == 'Sample_0'

    def test_samples_in_metadata(self, acquired_dataset):
        """Metadata-count overlap check."""
        assert acquired_dataset.samples_in_metadata == 8

    def test_validation_empty_counts(self, sample_metadata):
        """Reject empty count matrix."""
        from raptor.external_modules.acquisition import AcquiredDataset
        with pytest.raises(Exception):
            AcquiredDataset(
                counts_df=pd.DataFrame(),
                metadata=sample_metadata,
            )

    def test_validation_invalid_gene_id_type(self, sample_counts, sample_metadata):
        """Reject invalid gene_id_type."""
        from raptor.external_modules.acquisition import AcquiredDataset
        with pytest.raises(Exception):
            AcquiredDataset(
                counts_df=sample_counts,
                metadata=sample_metadata,
                gene_id_type='invalid_type',
            )

    def test_integrity_check_valid(self, acquired_dataset):
        """Integrity check on valid dataset."""
        report = acquired_dataset.validate_integrity()
        assert report['valid'] is True
        assert len(report['errors']) == 0

    def test_integrity_check_missing_metadata(self, sample_counts):
        """Integrity warning for mismatched metadata."""
        from raptor.external_modules.acquisition import AcquiredDataset
        partial_meta = pd.DataFrame(
            {'condition': ['ctrl'] * 4},
            index=[f'Sample_{i}' for i in range(4)],  # only 4 of 8
        )
        partial_meta.index.name = 'sample_id'
        ds = AcquiredDataset(counts_df=sample_counts, metadata=partial_meta)
        report = ds.validate_integrity()
        assert len(report['warnings']) > 0

    def test_filter_genes(self, acquired_dataset):
        """Gene filtering by expression level."""
        filtered = acquired_dataset.filter_genes(
            min_total_count=50,
            min_samples_detected=4,
        )
        assert filtered.n_genes < acquired_dataset.n_genes
        assert filtered.n_samples == acquired_dataset.n_samples
        assert filtered.n_genes > 0

    def test_subset_samples(self, acquired_dataset):
        """Sample subsetting."""
        subset = acquired_dataset.subset_samples(['Sample_0', 'Sample_1'])
        assert subset.n_samples == 2
        assert subset.n_genes == acquired_dataset.n_genes

    def test_subset_samples_invalid(self, acquired_dataset):
        """Reject subset with no matching samples."""
        with pytest.raises(Exception):
            acquired_dataset.subset_samples(['NONEXISTENT_1', 'NONEXISTENT_2'])

    def test_save_load_pickle(self, acquired_dataset, tmp_path):
        """Pickle round-trip."""
        path = tmp_path / 'test_ds.pkl'
        acquired_dataset.save(path)
        assert path.exists()

        from raptor.external_modules.acquisition import AcquiredDataset
        loaded = AcquiredDataset.load(path)
        assert loaded.n_genes == acquired_dataset.n_genes
        assert loaded.n_samples == acquired_dataset.n_samples
        assert loaded.accession == acquired_dataset.accession

    def test_save_load_parquet(self, acquired_dataset, tmp_path):
        """Parquet round-trip."""
        path = tmp_path / 'test_ds.parquet'
        acquired_dataset.save(path)
        assert path.exists()

        from raptor.external_modules.acquisition import AcquiredDataset
        loaded = AcquiredDataset.load(path)
        assert loaded.n_genes == acquired_dataset.n_genes
        assert loaded.n_samples == acquired_dataset.n_samples

    def test_export_csv(self, acquired_dataset, tmp_path):
        """CSV export."""
        result = acquired_dataset.export_csv(tmp_path / 'export')
        assert result['counts_file'].exists()
        assert result['metadata_file'].exists()
        assert result['info_file'].exists()

        # Verify CSV is readable
        reloaded = pd.read_csv(result['counts_file'], index_col=0)
        assert reloaded.shape == (500, 8)

    def test_from_csv(self, acquired_dataset, tmp_path):
        """Create from user-provided CSV files."""
        # Save CSVs
        counts_path = tmp_path / 'counts.csv'
        meta_path = tmp_path / 'meta.csv'
        acquired_dataset.counts_df.to_csv(counts_path, encoding='utf-8')
        acquired_dataset.metadata.to_csv(meta_path, encoding='utf-8')

        from raptor.external_modules.acquisition import AcquiredDataset
        ds = AcquiredDataset.from_csv(
            counts_path, meta_path,
            organism='Mus musculus',
        )
        assert ds.n_genes == 500
        assert ds.repository == 'user'

    def test_summary(self, acquired_dataset):
        """Summary string generation."""
        summary = acquired_dataset.summary()
        assert 'GSE99999' in summary
        assert 'GEO' in summary
        assert '500' in summary

    def test_repr(self, acquired_dataset):
        """String representation."""
        r = repr(acquired_dataset)
        assert 'GSE99999' in r
        assert '500' in r

    def test_default_source_info(self, sample_counts, sample_metadata):
        """Default source_info fields are populated."""
        from raptor.external_modules.acquisition import AcquiredDataset
        ds = AcquiredDataset(counts_df=sample_counts, metadata=sample_metadata)
        assert ds.source_info['repository'] == 'unknown'
        assert 'download_date' in ds.source_info

    def test_save_unsupported_format(self, acquired_dataset, tmp_path):
        """Reject unsupported save format."""
        with pytest.raises(Exception):
            acquired_dataset.save(tmp_path / 'bad.xlsx')


# =============================================================================
# TEST: PooledDataset
# =============================================================================

class TestPooledDataset:
    """Tests for PooledDataset data container."""

    @pytest.fixture
    def pooled_dataset(self, sample_counts, sample_metadata):
        """A PooledDataset with 2 studies."""
        from raptor.external_modules.acquisition import PooledDataset
        labels = pd.Series(
            ['GSE111'] * 4 + ['GSE222'] * 4,
            index=sample_counts.columns,
        )
        return PooledDataset(
            counts_df=sample_counts,
            metadata=sample_metadata,
            study_labels=labels,
            pooling_info={'method': 'inner', 'batch_correction': 'none'},
        )

    def test_creation(self, pooled_dataset):
        """Basic properties."""
        assert pooled_dataset.n_genes == 500
        assert pooled_dataset.n_samples == 8
        assert pooled_dataset.n_studies == 2

    def test_studies_list(self, pooled_dataset):
        """Study accession list."""
        studies = pooled_dataset.studies
        assert 'GSE111' in studies
        assert 'GSE222' in studies

    def test_samples_per_study(self, pooled_dataset):
        """Per-study sample counts."""
        sps = pooled_dataset.samples_per_study
        assert sps['GSE111'] == 4
        assert sps['GSE222'] == 4

    def test_get_study_samples(self, pooled_dataset):
        """Get samples from a specific study."""
        samples = pooled_dataset.get_study_samples('GSE111')
        assert len(samples) == 4

    def test_leave_one_study_out(self, pooled_dataset):
        """LOSO cross-validation split."""
        train, test = pooled_dataset.leave_one_study_out('GSE111')
        assert train.shape[1] == 4  # GSE222
        assert test.shape[1] == 4   # GSE111
        assert train.shape[0] == test.shape[0] == 500

    def test_save_load_pickle(self, pooled_dataset, tmp_path):
        """Pickle round-trip."""
        path = tmp_path / 'pool.pkl'
        pooled_dataset.save(path)

        from raptor.external_modules.acquisition import PooledDataset
        loaded = PooledDataset.load(path)
        assert loaded.n_studies == 2
        assert loaded.n_genes == 500

    def test_save_load_parquet(self, pooled_dataset, tmp_path):
        """Parquet round-trip."""
        path = tmp_path / 'pool.parquet'
        pooled_dataset.save(path)

        from raptor.external_modules.acquisition import PooledDataset
        loaded = PooledDataset.load(path)
        assert loaded.n_studies == 2

    def test_summary(self, pooled_dataset):
        """Summary string."""
        s = pooled_dataset.summary()
        assert 'GSE111' in s
        assert 'GSE222' in s
        assert '2' in s

    def test_study_column_in_metadata(self, pooled_dataset):
        """Study column auto-added to metadata."""
        assert 'study' in pooled_dataset.metadata.columns


# =============================================================================
# TEST: CacheManager
# =============================================================================

class TestCacheManager:
    """Tests for CacheManager."""

    @pytest.fixture
    def cache(self, tmp_path):
        from raptor.external_modules.acquisition import CacheManager
        return CacheManager(cache_dir=tmp_path / 'cache', enabled=True)

    @pytest.fixture
    def disabled_cache(self, tmp_path):
        from raptor.external_modules.acquisition import CacheManager
        return CacheManager(cache_dir=tmp_path / 'nocache', enabled=False)

    def test_init_creates_dirs(self, cache, tmp_path):
        """Cache directory structure created on init."""
        cache_dir = tmp_path / 'cache'
        assert cache_dir.exists()
        assert (cache_dir / 'geo').exists()
        assert (cache_dir / 'tcga').exists()
        assert (cache_dir / 'pools').exists()

    def test_disabled_mode(self, disabled_cache):
        """Disabled cache returns None/False for all ops."""
        assert disabled_cache.is_cached('GEO', 'GSE12345') is False
        assert disabled_cache.save_dataset(
            pd.DataFrame(), pd.DataFrame(), {}
        ) is None

    def test_save_and_load_dataset(self, cache, sample_counts, sample_metadata, source_info_geo):
        """Save then load a dataset."""
        path = cache.save_dataset(
            sample_counts, sample_metadata, source_info_geo, 'symbol'
        )
        assert path is not None
        assert cache.is_cached('GEO', 'GSE99999') is True

        loaded = cache.load_dataset('GEO', 'GSE99999')
        assert loaded is not None
        assert loaded['counts_df'].shape == sample_counts.shape
        assert loaded['gene_id_type'] == 'symbol'

    def test_load_uncached(self, cache):
        """Loading uncached dataset returns None."""
        assert cache.load_dataset('GEO', 'NONEXISTENT') is None

    def test_delete_dataset(self, cache, sample_counts, sample_metadata, source_info_geo):
        """Delete a cached dataset."""
        cache.save_dataset(sample_counts, sample_metadata, source_info_geo)
        assert cache.is_cached('GEO', 'GSE99999') is True

        result = cache.delete_dataset('GEO', 'GSE99999')
        assert result is True
        assert cache.is_cached('GEO', 'GSE99999') is False

    def test_delete_nonexistent(self, cache):
        """Deleting non-existent dataset returns False."""
        assert cache.delete_dataset('GEO', 'FAKE') is False

    def test_list_cached_datasets(self, cache, sample_counts, sample_metadata, source_info_geo):
        """List cached datasets."""
        cache.save_dataset(sample_counts, sample_metadata, source_info_geo)
        datasets = cache.list_cached_datasets()
        assert len(datasets) == 1
        assert datasets[0]['accession'] == 'GSE99999'

    def test_save_and_load_pool(self, cache, sample_counts, sample_metadata):
        """Save and load a pooled dataset."""
        labels = pd.Series(['S1'] * 4 + ['S2'] * 4, index=sample_counts.columns)
        path = cache.save_pool(
            'test_pool', sample_counts, sample_metadata, labels,
            {'method': 'inner'}, [{'accession': 'GSE111'}],
        )
        assert path is not None
        assert cache.is_pool_cached('test_pool') is True

        loaded = cache.load_pool('test_pool')
        assert loaded is not None
        assert loaded['counts_df'].shape == sample_counts.shape

    def test_csv_export(self, cache, sample_counts, sample_metadata, source_info_geo, tmp_path):
        """Export cached dataset as CSV."""
        cache.save_dataset(sample_counts, sample_metadata, source_info_geo)
        result = cache.export_dataset_csv('GEO', 'GSE99999', tmp_path / 'csv_out')
        assert result is not None
        assert result['counts_file'].exists()

    def test_cache_summary(self, cache, sample_counts, sample_metadata, source_info_geo):
        """Summary string."""
        cache.save_dataset(sample_counts, sample_metadata, source_info_geo)
        s = cache.cache_summary()
        assert 'GSE99999' in s

    def test_clear_all(self, cache, sample_counts, sample_metadata, source_info_geo):
        """Clear all requires confirm=True."""
        cache.save_dataset(sample_counts, sample_metadata, source_info_geo)

        # Without confirm — no-op
        cache.clear_all(confirm=False)
        assert cache.is_cached('GEO', 'GSE99999') is True

        # With confirm
        cache.clear_all(confirm=True)
        assert cache.is_cached('GEO', 'GSE99999') is False


# =============================================================================
# TEST: DataCatalog
# =============================================================================

class TestDataCatalog:
    """Tests for DataCatalog."""

    @pytest.fixture
    def catalog(self, tmp_path):
        from raptor.external_modules.acquisition import DataCatalog
        return DataCatalog(tmp_path / 'catalog')

    def test_register_dataset(self, catalog):
        """Register and retrieve a dataset."""
        catalog.register_dataset(
            repository='GEO', accession='GSE12345',
            organism='Homo sapiens', n_genes=20000, n_samples=12,
            description='Test liver cancer study',
            tags=['liver', 'cancer'],
        )
        assert catalog.n_datasets == 1
        entry = catalog.get_dataset('GEO', 'GSE12345')
        assert entry is not None
        assert entry['organism'] == 'Homo sapiens'

    def test_register_pool(self, catalog):
        """Register and retrieve a pool."""
        catalog.register_pool(
            pool_name='liver_pool',
            n_studies=3, n_genes=15000, n_samples=100,
            studies=['GSE111', 'GSE222', 'GSE333'],
        )
        assert catalog.n_pools == 1
        entry = catalog.get_pool('liver_pool')
        assert entry['n_studies'] == 3

    def test_update_existing(self, catalog):
        """Re-registering updates, not duplicates."""
        catalog.register_dataset(repository='GEO', accession='GSE111', n_samples=10)
        catalog.register_dataset(repository='GEO', accession='GSE111', n_samples=20)
        assert catalog.n_datasets == 1
        assert catalog.get_dataset('GEO', 'GSE111')['n_samples'] == 20

    def test_unregister(self, catalog):
        """Unregister dataset and pool."""
        catalog.register_dataset(repository='GEO', accession='GSE111')
        assert catalog.unregister_dataset('GEO', 'GSE111') is True
        assert catalog.n_datasets == 0
        assert catalog.unregister_dataset('GEO', 'GSE111') is False

    def test_search(self, catalog):
        """Free-text search across datasets and pools."""
        catalog.register_dataset(
            repository='GEO', accession='GSE12345',
            organism='Homo sapiens', description='Liver cancer RNA-seq',
            tags=['liver', 'cancer'],
        )
        catalog.register_dataset(
            repository='GEO', accession='GSE67890',
            organism='Mus musculus', description='Mouse pancreas',
            tags=['pancreas'],
        )
        catalog.register_pool(
            pool_name='liver_pool', description='Pooled liver datasets',
        )

        results = catalog.search('liver')
        assert len(results) == 2  # GSE12345 + liver_pool

        results = catalog.search('mouse pancreas')
        assert len(results) == 1
        assert results[0]['accession'] == 'GSE67890'

    def test_filter_by_repository(self, catalog):
        """Filter datasets by repository."""
        catalog.register_dataset(repository='GEO', accession='GSE111')
        catalog.register_dataset(repository='TCGA', accession='TCGA-LIHC')

        geo = catalog.list_datasets(repository='GEO')
        assert len(geo) == 1
        assert geo[0]['accession'] == 'GSE111'

    def test_filter_by_organism(self, catalog):
        """Filter datasets by organism."""
        catalog.register_dataset(repository='GEO', accession='GSE111', organism='Homo sapiens')
        catalog.register_dataset(repository='GEO', accession='GSE222', organism='Mus musculus')

        human = catalog.list_datasets(organism='sapiens')
        assert len(human) == 1

    def test_persistence(self, tmp_path):
        """Catalog persists across instances."""
        from raptor.external_modules.acquisition import DataCatalog

        cat1 = DataCatalog(tmp_path / 'persist')
        cat1.register_dataset(repository='GEO', accession='GSE111')

        cat2 = DataCatalog(tmp_path / 'persist')
        assert cat2.n_datasets == 1
        assert cat2.get_dataset('GEO', 'GSE111') is not None

    def test_summary(self, catalog):
        """Summary string."""
        catalog.register_dataset(repository='GEO', accession='GSE111', organism='Homo sapiens')
        s = catalog.summary()
        assert 'GSE111' in s


# =============================================================================
# TEST: BaseConnector + SearchResult (via MockConnector)
# =============================================================================

class TestBaseConnector:
    """Tests for BaseConnector via a mock subclass."""

    @pytest.fixture
    def mock_connector(self, tmp_path, sample_counts, sample_metadata):
        from raptor.external_modules.acquisition.base import BaseConnector, SearchResult
        from raptor.external_modules.acquisition import AcquiredDataset

        counts = sample_counts
        meta = sample_metadata

        class MockConnector(BaseConnector):
            REPOSITORY_NAME = 'MockRepo'

            def _search_api(self, query, organism=None, max_results=20):
                return [
                    SearchResult('MOCK001', 'Dataset one', 'Homo sapiens', 6),
                    SearchResult('MOCK002', 'Dataset two', 'Mus musculus', 4),
                ]

            def _download_api(self, accession, data_type='raw_counts'):
                return AcquiredDataset(
                    counts_df=counts.copy(),
                    metadata=meta.copy(),
                    source_info={
                        'repository': 'MockRepo',
                        'accession': accession,
                        'organism': 'Homo sapiens',
                    },
                    gene_id_type='symbol',
                )

        return MockConnector(cache=True, cache_dir=tmp_path / 'mock_cache')

    def test_search(self, mock_connector):
        """Search returns SearchResult objects."""
        results = mock_connector.search('test query')
        assert len(results) == 2
        assert results[0].accession == 'MOCK001'
        assert results[0].organism == 'Homo sapiens'

    def test_download_and_cache(self, mock_connector):
        """Download caches automatically."""
        ds = mock_connector.download('MOCK001')
        assert ds.n_genes == 500
        assert ds.accession == 'MOCK001'

        # Second download should hit cache
        assert mock_connector.is_cached('MOCK001') is True
        ds2 = mock_connector.download('MOCK001')
        assert ds2.n_genes == 500

    def test_download_force(self, mock_connector):
        """Force re-download bypasses cache."""
        mock_connector.download('MOCK001')
        ds = mock_connector.download('MOCK001', force=True)
        assert ds.n_genes == 500

    def test_download_multiple(self, mock_connector):
        """Batch download."""
        datasets = mock_connector.download_multiple(['MOCK001', 'MOCK002'])
        assert len(datasets) == 2

    def test_list_cached(self, mock_connector):
        """List cached datasets for this connector."""
        mock_connector.download('MOCK001')
        cached = mock_connector.list_cached()
        assert len(cached) >= 1

    def test_format_search_results(self, mock_connector):
        """Search results table formatting."""
        results = mock_connector.search('test')
        table = mock_connector.format_search_results(results)
        assert 'MOCK001' in table
        assert 'MOCK002' in table

    def test_format_empty_results(self):
        """Empty results formatting."""
        from raptor.external_modules.acquisition.base import BaseConnector
        table = BaseConnector.format_search_results([])
        assert 'No results' in table


class TestSearchResult:
    """Tests for SearchResult."""

    def test_creation(self):
        from raptor.external_modules.acquisition.base import SearchResult
        sr = SearchResult('GSE12345', 'My study', 'Homo sapiens', 24)
        assert sr.accession == 'GSE12345'
        assert sr.n_samples == 24

    def test_to_dict(self):
        from raptor.external_modules.acquisition.base import SearchResult
        sr = SearchResult('GSE12345', 'My study', 'Homo sapiens', 24, extra_field='test')
        d = sr.to_dict()
        assert d['accession'] == 'GSE12345'
        assert d['extra_field'] == 'test'


# =============================================================================
# TEST: GeneIDMapper
# =============================================================================

class TestGeneIDMapper:
    """Tests for GeneIDMapper."""

    def test_creation(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        mapper = GeneIDMapper('Homo sapiens')
        assert mapper.species == 'Homo sapiens'
        assert mapper.taxid == 9606

    def test_species_aliases(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        assert GeneIDMapper('human').taxid == 9606
        assert GeneIDMapper('mouse').taxid == 10090
        assert GeneIDMapper('Mus musculus').taxid == 10090

    def test_unknown_species_defaults_human(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        mapper = GeneIDMapper('Unknown alien species')
        assert mapper.taxid == 9606

    def test_detect_ensembl(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        ids = ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000146648']
        assert GeneIDMapper.detect_id_type(ids) == 'ensembl'

    def test_detect_symbol(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        ids = ['TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS']
        assert GeneIDMapper.detect_id_type(ids) == 'symbol'

    def test_detect_entrez(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        ids = ['7157', '672', '1956', '4609', '3845']
        assert GeneIDMapper.detect_id_type(ids) == 'entrez'

    def test_detect_empty(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        assert GeneIDMapper.detect_id_type([]) == 'unknown'

    def test_same_type_conversion(self):
        """No API call needed for same-type conversion."""
        from raptor.external_modules.acquisition import GeneIDMapper
        mapper = GeneIDMapper('Homo sapiens')
        result = mapper.convert(['TP53', 'BRCA1'], 'symbol', 'symbol')
        assert result == {'TP53': 'TP53', 'BRCA1': 'BRCA1'}

    def test_clear_cache(self):
        from raptor.external_modules.acquisition import GeneIDMapper
        mapper = GeneIDMapper('Homo sapiens')
        mapper.convert(['TP53'], 'symbol', 'symbol')
        mapper.clear_cache()
        assert repr(mapper).endswith('cached_mappings=0)')

    def test_invalid_id_type(self):
        """Reject invalid ID types."""
        from raptor.external_modules.acquisition import GeneIDMapper
        mapper = GeneIDMapper('Homo sapiens')
        with pytest.raises(Exception):
            mapper.convert(['TP53'], 'symbol', 'invalid_type')

        with pytest.raises(Exception):
            mapper.convert(['TP53'], 'invalid_type', 'symbol')


# =============================================================================
# TEST: Connectors (constructors, no network)
# =============================================================================

class TestConnectorConstructors:
    """Test connector creation without network access."""

    def test_geo_constructor(self):
        from raptor.external_modules.acquisition import GEOConnector
        geo = GEOConnector(email='test@example.com', cache=False)
        assert geo.REPOSITORY_NAME == 'GEO'
        assert 'test@example.com' in repr(geo)

    def test_tcga_constructor(self):
        from raptor.external_modules.acquisition import TCGAConnector
        tcga = TCGAConnector(cache=False)
        assert tcga.REPOSITORY_NAME == 'TCGA'
        projects = tcga.list_projects()
        assert 'TCGA-BRCA' in projects
        assert 'TCGA-LIHC' in projects
        assert len(projects) >= 15

    def test_arrayexpress_constructor(self):
        from raptor.external_modules.acquisition import ArrayExpConnector
        ae = ArrayExpConnector(cache=False)
        assert ae.REPOSITORY_NAME == 'ArrayExpress'

    def test_sra_constructor(self):
        from raptor.external_modules.acquisition import SRAConnector
        sra = SRAConnector(email='test@example.com', cache=False)
        assert sra.REPOSITORY_NAME == 'SRA'

    def test_geo_detect_gene_id_type(self):
        """GEOConnector static gene ID detection."""
        from raptor.external_modules.acquisition.geo import GEOConnector

        ensembl_df = pd.DataFrame(
            {'S1': [1, 2]},
            index=['ENSG00000141510', 'ENSG00000012048'],
        )
        assert GEOConnector._detect_gene_id_type(ensembl_df) == 'ensembl'

        symbol_df = pd.DataFrame(
            {'S1': [1, 2]},
            index=['TP53', 'BRCA1'],
        )
        assert GEOConnector._detect_gene_id_type(symbol_df) == 'symbol'


# =============================================================================
# TEST: SRA GEO Cross-Referencing
# =============================================================================

class TestSRAGeoLookup:
    """Tests for SRA GSM extraction and GSE cross-referencing (all offline)."""

    @pytest.fixture
    def sra_connector(self):
        from raptor.external_modules.acquisition import SRAConnector
        return SRAConnector(email='test@example.com', cache=False)

    @pytest.fixture
    def run_table_with_gsm(self):
        """Run table where sample_alias contains GSM IDs."""
        return pd.DataFrame({
            'run_accession': ['SRR36229312', 'SRR36229313', 'SRR36229314', 'SRR36229315'],
            'sample_alias': ['GSM8716844', 'GSM8716843', 'GSM8716846', 'GSM8716845'],
            'sample_title': [
                'WT mice and TAM treatment (CfD-TAM)',
                'WT mice and Oil treatment (CfD-OIL)',
                'PS19 mice and TAM treatment (PCfD-TAM)',
                'PS19 mice and Oil treatment (PCfD-OIL)',
            ],
            'instrument_model': ['Illumina NovaSeq X Plus'] * 4,
            'library_layout': ['PAIRED'] * 4,
            'read_count': [602809299, 260825463, 218231416, 386964084],
            'scientific_name': ['Mus musculus'] * 4,
        })

    @pytest.fixture
    def run_table_no_gsm(self):
        """Run table without GSM IDs (pure SRA study)."""
        return pd.DataFrame({
            'run_accession': ['SRR100001', 'SRR100002'],
            'sample_alias': ['my_sample_1', 'my_sample_2'],
            'sample_title': ['Control replicate 1', 'Treatment replicate 1'],
            'instrument_model': ['Illumina HiSeq 2500'] * 2,
            'library_layout': ['PAIRED'] * 2,
            'read_count': [50000000, 48000000],
            'scientific_name': ['Homo sapiens'] * 2,
        })

    @pytest.fixture
    def run_table_gsm_in_title(self):
        """Run table where GSM appears in sample_title but not sample_alias."""
        return pd.DataFrame({
            'run_accession': ['SRR200001', 'SRR200002'],
            'sample_alias': ['SAMN12345', 'SAMN12346'],
            'sample_title': ['GSM7777777 brain tissue control', 'GSM7777778 brain tissue disease'],
            'instrument_model': ['Illumina NovaSeq 6000'] * 2,
            'library_layout': ['SINGLE'] * 2,
            'read_count': [30000000, 35000000],
            'scientific_name': ['Homo sapiens'] * 2,
        })

    # ----- _extract_gsm_ids tests -----

    def test_extract_gsm_from_alias(self, sra_connector, run_table_with_gsm):
        """Extract GSM IDs from sample_alias column."""
        gsm_ids = sra_connector._extract_gsm_ids(run_table_with_gsm)
        assert len(gsm_ids) == 4
        assert 'GSM8716844' in gsm_ids
        assert 'GSM8716843' in gsm_ids
        assert 'GSM8716846' in gsm_ids
        assert 'GSM8716845' in gsm_ids

    def test_extract_gsm_from_title(self, sra_connector, run_table_gsm_in_title):
        """Extract GSM IDs from sample_title when not in sample_alias."""
        gsm_ids = sra_connector._extract_gsm_ids(run_table_gsm_in_title)
        assert len(gsm_ids) == 2
        assert 'GSM7777777' in gsm_ids
        assert 'GSM7777778' in gsm_ids

    def test_extract_gsm_none_found(self, sra_connector, run_table_no_gsm):
        """Return empty list when no GSM IDs present."""
        gsm_ids = sra_connector._extract_gsm_ids(run_table_no_gsm)
        assert gsm_ids == []

    def test_extract_gsm_custom_columns(self, sra_connector, run_table_with_gsm):
        """Only search specified columns."""
        # Search only sample_title (which has no GSMs in this fixture)
        gsm_ids = sra_connector._extract_gsm_ids(
            run_table_with_gsm, columns=['sample_title']
        )
        assert gsm_ids == []

    def test_extract_gsm_missing_column(self, sra_connector):
        """Gracefully handle missing columns."""
        df = pd.DataFrame({'run_accession': ['SRR1']})
        gsm_ids = sra_connector._extract_gsm_ids(df)
        assert gsm_ids == []

    def test_extract_gsm_returns_sorted_unique(self, sra_connector):
        """Duplicate GSMs across rows should be deduplicated."""
        df = pd.DataFrame({
            'sample_alias': ['GSM100', 'GSM100', 'GSM200', 'GSM200'],
        })
        gsm_ids = sra_connector._extract_gsm_ids(df)
        assert gsm_ids == ['GSM100', 'GSM200']

    # ----- _gse_from_gsm tests (mocked Entrez) -----

    def test_gse_from_gsm_success(self, sra_connector):
        """Successfully look up GSE from GSM via mocked Entrez."""
        mock_entrez = MagicMock()

        # Mock esearch
        mock_search_handle = MagicMock()
        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [
            {'IdList': ['300099999']},  # esearch result
            [{'Accession': 'GSE306761', 'Id': '300099999'}],  # esummary result
        ]
        mock_search_handle.close = MagicMock()

        # Mock esummary
        mock_summary_handle = MagicMock()
        mock_entrez.esummary.return_value = mock_summary_handle
        mock_summary_handle.close = MagicMock()

        with patch(
            'raptor.external_modules.acquisition.sra._check_entrez',
            return_value=mock_entrez
        ):
            gse = sra_connector._gse_from_gsm('GSM8716844')
            assert gse == 'GSE306761'

    def test_gse_from_gsm_not_found(self, sra_connector):
        """Return empty string when GSM has no linked GSE."""
        mock_entrez = MagicMock()
        mock_search_handle = MagicMock()
        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.return_value = {'IdList': []}
        mock_search_handle.close = MagicMock()

        with patch(
            'raptor.external_modules.acquisition.sra._check_entrez',
            return_value=mock_entrez
        ):
            gse = sra_connector._gse_from_gsm('GSM9999999')
            assert gse == ''

    def test_gse_from_gsm_entrez_error(self, sra_connector):
        """Gracefully handle Entrez errors."""
        with patch(
            'raptor.external_modules.acquisition.sra._check_entrez',
            side_effect=Exception("Network error")
        ):
            gse = sra_connector._gse_from_gsm('GSM8716844')
            assert gse == ''

    # ----- find_linked_gse tests -----

    def test_find_linked_gse_with_run_table(self, sra_connector, run_table_with_gsm):
        """find_linked_gse uses provided run table, finds GSE."""
        with patch.object(sra_connector, '_gse_from_gsm', return_value='GSE306761'):
            gse = sra_connector.find_linked_gse('SRP555825', run_table=run_table_with_gsm)
            assert gse == 'GSE306761'

    def test_find_linked_gse_no_gsm(self, sra_connector, run_table_no_gsm):
        """Return empty when run table has no GSM IDs."""
        gse = sra_connector.find_linked_gse('SRP999999', run_table=run_table_no_gsm)
        assert gse == ''

    def test_find_linked_gse_fallback_second_gsm(self, sra_connector, run_table_with_gsm):
        """Try second GSM if first lookup fails."""
        call_count = [0]

        def mock_lookup(gsm_id):
            call_count[0] += 1
            if call_count[0] == 1:
                return ''  # First GSM fails
            return 'GSE306761'  # Second succeeds

        with patch.object(sra_connector, '_gse_from_gsm', side_effect=mock_lookup):
            gse = sra_connector.find_linked_gse('SRP555825', run_table=run_table_with_gsm)
            assert gse == 'GSE306761'
            assert call_count[0] == 2

    def test_find_linked_gse_fetches_run_table(self, sra_connector):
        """When no run_table provided, fetches it automatically."""
        mock_table = pd.DataFrame({
            'sample_alias': ['GSM100001'],
            'run_accession': ['SRR100001'],
        })

        with patch.object(sra_connector, 'get_run_table', return_value=mock_table):
            with patch.object(sra_connector, '_gse_from_gsm', return_value='GSE999999'):
                gse = sra_connector.find_linked_gse('SRP123456')
                assert gse == 'GSE999999'

    # ----- _download_api GSE population tests -----

    def test_download_api_populates_gse(self, sra_connector):
        """_download_api should store GSE and GSM IDs in source_info."""
        mock_table = pd.DataFrame({
            'run_accession': ['SRR100001', 'SRR100002'],
            'sample_alias': ['GSM100001', 'GSM100002'],
            'sample_title': ['Control', 'Treatment'],
            'instrument_model': ['Illumina NovaSeq 6000', 'Illumina NovaSeq 6000'],
            'library_layout': ['PAIRED', 'PAIRED'],
            'read_count': [50000000, 48000000],
            'scientific_name': ['Homo sapiens', 'Homo sapiens'],
        })

        with patch.object(sra_connector, 'get_run_table', return_value=mock_table):
            with patch.object(sra_connector, 'find_linked_gse', return_value='GSE123456'):
                ds = sra_connector._download_api('SRP123456')
                assert ds.source_info.get('gse') == 'GSE123456'
                assert 'GSM100001' in ds.source_info.get('gsm_ids', [])
                assert 'GSM100002' in ds.source_info.get('gsm_ids', [])
                assert ds.source_info.get('repository') == 'SRA'
                assert ds.source_info.get('n_runs') == 2

    def test_download_api_no_gse(self, sra_connector):
        """_download_api should work fine when no GSE is found."""
        mock_table = pd.DataFrame({
            'run_accession': ['SRR100001'],
            'sample_alias': ['my_sample_1'],
            'sample_title': ['Control'],
            'instrument_model': ['Illumina HiSeq 2500'],
            'library_layout': ['PAIRED'],
            'read_count': [50000000],
            'scientific_name': ['Homo sapiens'],
        })

        with patch.object(sra_connector, 'get_run_table', return_value=mock_table):
            ds = sra_connector._download_api('SRP999999')
            assert ds.source_info.get('gse') is None or ds.source_info.get('gse', '') == ''
            assert ds.source_info.get('gsm_ids') is None or ds.source_info.get('gsm_ids', []) == []
            assert ds.source_info.get('repository') == 'SRA'

    # ----- SearchResult extra fields tests -----

    def test_search_result_sra_extra_fields(self):
        """SRA SearchResult should carry GEO link fields."""
        from raptor.external_modules.acquisition.base import SearchResult
        result = SearchResult(
            accession='SRP555825',
            title='Test study',
            organism='Mus musculus',
            n_samples=4,
            platform='Illumina NovaSeq X Plus',
            description='',
            repository='SRA',
            gse='GSE306761',
            has_geo_link=True,
            geo_hint='GSE306761',
            n_runs=8,
        )
        assert result.extra.get('gse') == 'GSE306761'
        assert result.extra.get('has_geo_link') is True
        assert result.extra.get('n_runs') == 8

    def test_search_result_sra_no_gse(self):
        """SRA SearchResult without GSE should still work."""
        from raptor.external_modules.acquisition.base import SearchResult
        result = SearchResult(
            accession='SRP999999',
            title='Pure SRA study',
            organism='Homo sapiens',
            n_samples=10,
            platform='Illumina HiSeq 2500',
            description='',
            repository='SRA',
            gse='',
            has_geo_link=False,
        )
        assert result.extra.get('gse') == ''
        assert result.extra.get('has_geo_link') is False


# =============================================================================
# TEST: PoolingEngine
# =============================================================================

class TestPoolingEngine:
    """Tests for PoolingEngine."""

    @pytest.fixture
    def engine(self):
        from raptor.external_modules.acquisition import PoolingEngine
        return PoolingEngine(target_gene_id='symbol', species='Homo sapiens')

    def test_creation(self, engine):
        assert 'symbol' in repr(engine)
        assert 'Homo sapiens' in repr(engine)

    def test_invalid_target_type(self):
        from raptor.external_modules.acquisition import PoolingEngine
        with pytest.raises(Exception):
            PoolingEngine(target_gene_id='invalid_type')

    def test_merge_inner(self, engine, acquired_dataset, acquired_dataset_2):
        """Inner join keeps only common genes."""
        pool = engine.merge(
            [acquired_dataset, acquired_dataset_2],
            method='inner',
        )
        # GENE_0200 to GENE_0499 = 300 common genes
        assert pool.n_genes == 300
        assert pool.n_samples == 14  # 8 + 6
        assert pool.n_studies == 2
        assert 'GSE99999' in pool.studies
        assert 'GSE88888' in pool.studies

    def test_merge_outer(self, engine, acquired_dataset, acquired_dataset_2):
        """Outer join keeps all genes."""
        pool = engine.merge(
            [acquired_dataset, acquired_dataset_2],
            method='outer',
        )
        # 500 + 500 - 300 overlap = 700
        assert pool.n_genes == 700
        assert pool.n_samples == 14

    def test_merge_minimum_datasets(self, engine, acquired_dataset):
        """Reject single dataset."""
        with pytest.raises(Exception):
            engine.merge([acquired_dataset])

    def test_merge_invalid_method(self, engine, acquired_dataset, acquired_dataset_2):
        """Reject invalid join method."""
        with pytest.raises(Exception):
            engine.merge([acquired_dataset, acquired_dataset_2], method='left')

    def test_merge_invalid_batch_correction(self, engine, acquired_dataset, acquired_dataset_2):
        """Reject invalid batch correction."""
        with pytest.raises(Exception):
            engine.merge(
                [acquired_dataset, acquired_dataset_2],
                batch_correction='invalid_method',
            )

    def test_batch_correction_combat_fallback(self, engine, acquired_dataset, acquired_dataset_2):
        """ComBat falls back to median centering if pycombat not installed."""
        pool = engine.merge(
            [acquired_dataset, acquired_dataset_2],
            batch_correction='combat',
        )
        assert pool.n_genes > 0
        assert pool.pooling_info['batch_correction'] == 'combat'

    def test_batch_correction_quantile(self, engine, acquired_dataset, acquired_dataset_2):
        """Quantile normalization."""
        pool = engine.merge(
            [acquired_dataset, acquired_dataset_2],
            batch_correction='quantile',
        )
        assert pool.n_genes > 0

    def test_batch_correction_median_ratio(self, engine, acquired_dataset, acquired_dataset_2):
        """Median ratio normalization."""
        pool = engine.merge(
            [acquired_dataset, acquired_dataset_2],
            batch_correction='median_ratio',
        )
        assert pool.n_genes > 0

    def test_sample_conflict_resolution(self, engine):
        """Conflicting sample IDs get prefixed."""
        from raptor.external_modules.acquisition import AcquiredDataset
        np.random.seed(42)

        shared_cols = ['S1', 'S2', 'S3']

        ds1 = AcquiredDataset(
            counts_df=pd.DataFrame(
                np.random.poisson(10, (100, 3)),
                index=[f'G{i}' for i in range(100)],
                columns=shared_cols,
            ).rename_axis('gene_id'),
            metadata=pd.DataFrame(
                {'cond': ['a', 'b', 'c']}, index=shared_cols,
            ).rename_axis('sample_id'),
            source_info={'repository': 'GEO', 'accession': 'GSE111', 'organism': 'Homo sapiens'},
            gene_id_type='symbol',
        )

        ds2 = AcquiredDataset(
            counts_df=pd.DataFrame(
                np.random.poisson(12, (100, 3)),
                index=[f'G{i}' for i in range(100)],
                columns=shared_cols,  # SAME as ds1!
            ).rename_axis('gene_id'),
            metadata=pd.DataFrame(
                {'cond': ['x', 'y', 'z']}, index=shared_cols,
            ).rename_axis('sample_id'),
            source_info={'repository': 'GEO', 'accession': 'GSE222', 'organism': 'Homo sapiens'},
            gene_id_type='symbol',
        )

        pool = engine.merge([ds1, ds2], method='inner')
        assert pool.n_samples == 6  # 3 + 3, no name collision
        # Check that prefixes were added
        assert any('GSE111' in s for s in pool.counts_df.columns)
        assert any('GSE222' in s for s in pool.counts_df.columns)

    def test_leave_one_study_out_from_pool(self, engine, acquired_dataset, acquired_dataset_2):
        """LOSO validation on pooled result."""
        pool = engine.merge([acquired_dataset, acquired_dataset_2], method='inner')
        train, test = pool.leave_one_study_out('GSE99999')
        assert train.shape[1] == 6   # GSE88888 samples
        assert test.shape[1] == 8    # GSE99999 samples

    def test_pool_auto_name(self, engine, acquired_dataset, acquired_dataset_2):
        """Auto-generated pool name when none provided."""
        pool = engine.merge([acquired_dataset, acquired_dataset_2])
        assert pool.pooling_info.get('timestamp') is not None

    def test_min_common_genes_threshold(self, engine):
        """Reject pools with too few common genes."""
        from raptor.external_modules.acquisition import AcquiredDataset
        np.random.seed(42)

        # Two datasets with NO gene overlap
        ds1 = AcquiredDataset(
            counts_df=pd.DataFrame(
                np.random.poisson(10, (50, 3)),
                index=[f'GENEA_{i}' for i in range(50)],
                columns=['S1', 'S2', 'S3'],
            ).rename_axis('gene_id'),
            metadata=pd.DataFrame(index=['S1', 'S2', 'S3']).rename_axis('sample_id'),
            source_info={'repository': 'GEO', 'accession': 'DS1', 'organism': 'Homo sapiens'},
            gene_id_type='symbol',
        )
        ds2 = AcquiredDataset(
            counts_df=pd.DataFrame(
                np.random.poisson(10, (50, 3)),
                index=[f'GENEB_{i}' for i in range(50)],
                columns=['S4', 'S5', 'S6'],
            ).rename_axis('gene_id'),
            metadata=pd.DataFrame(index=['S4', 'S5', 'S6']).rename_axis('sample_id'),
            source_info={'repository': 'GEO', 'accession': 'DS2', 'organism': 'Homo sapiens'},
            gene_id_type='symbol',
        )

        with pytest.raises(Exception, match="common genes"):
            engine.merge([ds1, ds2], method='inner', min_common_genes=100)


# =============================================================================
# TEST: Module availability and imports
# =============================================================================

class TestModuleAvailability:
    """Test that everything imports correctly."""

    def test_get_available_components(self):
        from raptor.external_modules.acquisition import get_available_components
        components = get_available_components()
        assert components['datasets'] is True
        assert components['cache'] is True
        assert components['catalog'] is True
        assert components['pooling_engine'] is True

    def test_top_level_imports(self):
        """All public classes importable from top-level."""
        from raptor.external_modules.acquisition import (
            AcquiredDataset,
            PooledDataset,
            CacheManager,
            DataCatalog,
            GEOConnector,
            TCGAConnector,
            ArrayExpConnector,
            SRAConnector,
            GeneIDMapper,
            PoolingEngine,
        )
        assert AcquiredDataset is not None
        assert PooledDataset is not None
        assert CacheManager is not None
        assert GEOConnector is not None

    def test_constants_exported(self):
        from raptor.external_modules.acquisition import (
            SUPPORTED_REPOSITORIES,
            SUPPORTED_DATA_TYPES,
            SUPPORTED_ORGANISMS,
            DEFAULT_CACHE_DIR,
        )
        assert 'GEO' in SUPPORTED_REPOSITORIES
        assert 'TCGA' in SUPPORTED_REPOSITORIES
        assert 'raw_counts' in SUPPORTED_DATA_TYPES
