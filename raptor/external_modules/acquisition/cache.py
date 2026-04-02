"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - Cache Manager

Hybrid caching system for acquired datasets. Stores data on disk by default
(~/.raptor/data/) with Parquet as the primary format. Supports in-memory
mode and custom cache directories.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import json
import shutil
from pathlib import Path
from datetime import datetime
from typing import Optional, Union, Dict, Any

import pandas as pd

# RAPTOR imports
try:
    from raptor.utils.errors import ValidationError, RAPTORError
except ImportError:
    class ValidationError(Exception):
        def __init__(self, parameter: str, message: str = "", hint: str = ""):
            self.parameter = parameter
            self.message = message
            self.hint = hint
            super().__init__(f"{parameter}: {message}. {hint}")

    class RAPTORError(Exception):
        pass


logger = logging.getLogger(__name__)


# =============================================================================
# Constants
# =============================================================================

DEFAULT_CACHE_DIR = Path.home() / '.raptor' / 'data'

REPOSITORY_SUBDIRS = {
    'GEO': 'geo',
    'TCGA': 'tcga',
    'ArrayExpress': 'arrayexpress',
    'SRA': 'sra',
}

POOLS_SUBDIR = 'pools'


# =============================================================================
# CacheManager
# =============================================================================

class CacheManager:
    """
    Manages on-disk caching of acquired datasets.

    Handles all file I/O for the acquisition subpackage. Supports
    Parquet (primary, fast) and CSV (export, portable) formats.

    Parameters
    ----------
    cache_dir : str or Path, optional
        Root cache directory. Defaults to ~/.raptor/data/
    enabled : bool
        If False, operates in memory-only mode (no disk writes).

    Examples
    --------
    >>> cache = CacheManager()  # uses ~/.raptor/data/
    >>> cache = CacheManager(enabled=False)  # in-memory only
    >>> cache = CacheManager(cache_dir="/scratch/raptor_data")  # custom
    """

    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        enabled: bool = True,
    ):
        self.enabled = enabled
        self._cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR

        if self.enabled:
            self._init_cache_dir()

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------

    def _init_cache_dir(self) -> None:
        """Create the cache directory structure."""
        try:
            self._cache_dir.mkdir(parents=True, exist_ok=True)
            for subdir in REPOSITORY_SUBDIRS.values():
                (self._cache_dir / subdir).mkdir(exist_ok=True)
            (self._cache_dir / POOLS_SUBDIR).mkdir(exist_ok=True)
            logger.debug(f"Cache directory initialized: {self._cache_dir}")
        except PermissionError:
            logger.warning(
                f"Cannot create cache directory {self._cache_dir}. "
                f"Falling back to in-memory mode."
            )
            self.enabled = False

    @property
    def cache_dir(self) -> Path:
        """Root cache directory path."""
        return self._cache_dir

    # -------------------------------------------------------------------------
    # Path helpers
    # -------------------------------------------------------------------------

    def _repo_dir(self, repository: str) -> Path:
        """Get the subdirectory for a repository."""
        subdir = REPOSITORY_SUBDIRS.get(repository, repository.lower())
        return self._cache_dir / subdir

    def _dataset_dir(self, repository: str, accession: str) -> Path:
        """Get the directory for a specific dataset."""
        safe_accession = accession.replace('/', '_').replace('\\', '_')
        return self._repo_dir(repository) / safe_accession

    def _pool_dir(self, pool_name: str) -> Path:
        """Get the directory for a pooled dataset."""
        safe_name = pool_name.replace('/', '_').replace('\\', '_')
        return self._cache_dir / POOLS_SUBDIR / safe_name

    # -------------------------------------------------------------------------
    # Check / Exists
    # -------------------------------------------------------------------------

    def is_cached(self, repository: str, accession: str) -> bool:
        """
        Check if a dataset is already cached.

        Parameters
        ----------
        repository : str
            Repository name (GEO, TCGA, etc.)
        accession : str
            Dataset accession ID.

        Returns
        -------
        bool
        """
        if not self.enabled:
            return False

        dataset_dir = self._dataset_dir(repository, accession)
        manifest = dataset_dir / 'manifest.json'
        counts_file = dataset_dir / 'counts.parquet'

        return manifest.exists() and counts_file.exists()

    def is_pool_cached(self, pool_name: str) -> bool:
        """Check if a pooled dataset is cached."""
        if not self.enabled:
            return False

        pool_dir = self._pool_dir(pool_name)
        return (pool_dir / 'pooled_counts.parquet').exists()

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------

    def save_dataset(
        self,
        counts_df: pd.DataFrame,
        metadata: pd.DataFrame,
        source_info: Dict[str, Any],
        gene_id_type: str = 'symbol',
    ) -> Optional[Path]:
        """
        Save an acquired dataset to cache.

        Parameters
        ----------
        counts_df : pd.DataFrame
            Count matrix (genes x samples).
        metadata : pd.DataFrame
            Sample metadata.
        source_info : dict
            Provenance information (must include 'repository' and 'accession').
        gene_id_type : str
            Gene identifier type.

        Returns
        -------
        Path or None
            Path to cached dataset directory, or None if caching is disabled.
        """
        if not self.enabled:
            logger.debug("Cache disabled, skipping save")
            return None

        repository = source_info.get('repository', 'unknown')
        accession = source_info.get('accession', 'unknown')
        dataset_dir = self._dataset_dir(repository, accession)
        dataset_dir.mkdir(parents=True, exist_ok=True)

        # Save counts as Parquet (primary format)
        counts_path = dataset_dir / 'counts.parquet'
        counts_df.to_parquet(counts_path, engine='pyarrow')

        # Save metadata as Parquet
        meta_path = dataset_dir / 'metadata.parquet'
        metadata.to_parquet(meta_path, engine='pyarrow')

        # Save manifest
        manifest = {
            **source_info,
            'gene_id_type': gene_id_type,
            'n_genes': len(counts_df),
            'n_samples': len(counts_df.columns),
            'cached_date': datetime.now().isoformat(),
            'files': {
                'counts': 'counts.parquet',
                'metadata': 'metadata.parquet',
            },
        }
        manifest_path = dataset_dir / 'manifest.json'
        with open(manifest_path, 'w', encoding='utf-8') as f:
            json.dump(manifest, f, indent=2, default=str)

        logger.info(f"Cached {accession} → {dataset_dir}")
        return dataset_dir

    def save_pool(
        self,
        pool_name: str,
        counts_df: pd.DataFrame,
        metadata: pd.DataFrame,
        study_labels: pd.Series,
        pooling_info: Dict[str, Any],
        component_datasets: list,
    ) -> Optional[Path]:
        """
        Save a pooled dataset to cache.

        Parameters
        ----------
        pool_name : str
            Name for the pooled dataset.
        counts_df : pd.DataFrame
            Merged count matrix.
        metadata : pd.DataFrame
            Merged metadata.
        study_labels : pd.Series
            Sample-to-study mapping.
        pooling_info : dict
            Pooling parameters.
        component_datasets : list
            Info about constituent datasets.

        Returns
        -------
        Path or None
            Path to cached pool directory, or None if disabled.
        """
        if not self.enabled:
            return None

        pool_dir = self._pool_dir(pool_name)
        pool_dir.mkdir(parents=True, exist_ok=True)

        counts_df.to_parquet(pool_dir / 'pooled_counts.parquet', engine='pyarrow')
        metadata.to_parquet(pool_dir / 'harmonized_metadata.parquet', engine='pyarrow')

        pool_manifest = {
            'pool_name': pool_name,
            'n_studies': study_labels.nunique(),
            'n_genes': len(counts_df),
            'n_samples': len(counts_df.columns),
            'study_labels': study_labels.to_dict(),
            'pooling_info': pooling_info,
            'component_datasets': component_datasets,
            'cached_date': datetime.now().isoformat(),
        }
        with open(pool_dir / 'pool_manifest.json', 'w', encoding='utf-8') as f:
            json.dump(pool_manifest, f, indent=2, default=str)

        logger.info(f"Cached pool '{pool_name}' → {pool_dir}")
        return pool_dir

    # -------------------------------------------------------------------------
    # Load
    # -------------------------------------------------------------------------

    def load_dataset(
        self,
        repository: str,
        accession: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Load a cached dataset.

        Parameters
        ----------
        repository : str
            Repository name.
        accession : str
            Dataset accession ID.

        Returns
        -------
        dict or None
            Dictionary with keys: counts_df, metadata, source_info, gene_id_type.
            Returns None if not cached.
        """
        if not self.is_cached(repository, accession):
            return None

        dataset_dir = self._dataset_dir(repository, accession)

        counts_df = pd.read_parquet(dataset_dir / 'counts.parquet', engine='pyarrow')
        metadata = pd.read_parquet(dataset_dir / 'metadata.parquet', engine='pyarrow')

        with open(dataset_dir / 'manifest.json', 'r', encoding='utf-8') as f:
            manifest = json.load(f)

        gene_id_type = manifest.pop('gene_id_type', 'symbol')
        # Remove cache-specific keys from source_info
        source_info = {
            k: v for k, v in manifest.items()
            if k not in ('n_genes', 'n_samples', 'cached_date', 'files')
        }

        logger.info(f"Loaded {accession} from cache")
        return {
            'counts_df': counts_df,
            'metadata': metadata,
            'source_info': source_info,
            'gene_id_type': gene_id_type,
        }

    def load_pool(self, pool_name: str) -> Optional[Dict[str, Any]]:
        """
        Load a cached pooled dataset.

        Parameters
        ----------
        pool_name : str
            Pool name.

        Returns
        -------
        dict or None
            Dictionary with keys: counts_df, metadata, study_labels,
            pooling_info, component_datasets. Returns None if not cached.
        """
        if not self.is_pool_cached(pool_name):
            return None

        pool_dir = self._pool_dir(pool_name)

        counts_df = pd.read_parquet(pool_dir / 'pooled_counts.parquet', engine='pyarrow')
        metadata = pd.read_parquet(pool_dir / 'harmonized_metadata.parquet', engine='pyarrow')

        with open(pool_dir / 'pool_manifest.json', 'r', encoding='utf-8') as f:
            manifest = json.load(f)

        study_labels = pd.Series(manifest.get('study_labels', {}))

        logger.info(f"Loaded pool '{pool_name}' from cache")
        return {
            'counts_df': counts_df,
            'metadata': metadata,
            'study_labels': study_labels,
            'pooling_info': manifest.get('pooling_info', {}),
            'component_datasets': manifest.get('component_datasets', []),
        }

    # -------------------------------------------------------------------------
    # Export (CSV)
    # -------------------------------------------------------------------------

    def export_dataset_csv(
        self,
        repository: str,
        accession: str,
        output_dir: Union[str, Path],
    ) -> Optional[Dict[str, Path]]:
        """
        Export a cached dataset as CSV files.

        Parameters
        ----------
        repository : str
            Repository name.
        accession : str
            Dataset accession ID.
        output_dir : str or Path
            Directory for CSV output.

        Returns
        -------
        dict or None
            Paths to exported files, or None if not cached.
        """
        data = self.load_dataset(repository, accession)
        if data is None:
            logger.warning(f"{accession} not found in cache")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        prefix = accession.replace('/', '_')

        counts_path = output_dir / f"{prefix}_counts.csv"
        data['counts_df'].to_csv(counts_path, encoding='utf-8')

        meta_path = output_dir / f"{prefix}_metadata.csv"
        data['metadata'].to_csv(meta_path, encoding='utf-8')

        info_path = output_dir / f"{prefix}_info.json"
        with open(info_path, 'w', encoding='utf-8') as f:
            json.dump(data['source_info'], f, indent=2, default=str)

        logger.info(f"Exported {accession} as CSV → {output_dir}")
        return {
            'counts_file': counts_path,
            'metadata_file': meta_path,
            'info_file': info_path,
        }

    # -------------------------------------------------------------------------
    # Delete / Cleanup
    # -------------------------------------------------------------------------

    def delete_dataset(self, repository: str, accession: str) -> bool:
        """
        Remove a cached dataset.

        Parameters
        ----------
        repository : str
            Repository name.
        accession : str
            Dataset accession ID.

        Returns
        -------
        bool
            True if deleted, False if not found.
        """
        if not self.is_cached(repository, accession):
            return False

        dataset_dir = self._dataset_dir(repository, accession)
        shutil.rmtree(dataset_dir)
        logger.info(f"Deleted cached dataset: {accession}")
        return True

    def delete_pool(self, pool_name: str) -> bool:
        """Remove a cached pooled dataset."""
        if not self.is_pool_cached(pool_name):
            return False

        pool_dir = self._pool_dir(pool_name)
        shutil.rmtree(pool_dir)
        logger.info(f"Deleted cached pool: {pool_name}")
        return True

    def clear_all(self, confirm: bool = False) -> None:
        """
        Delete all cached data.

        Parameters
        ----------
        confirm : bool
            Must be True to proceed. Safety check.
        """
        if not confirm:
            logger.warning("clear_all() requires confirm=True")
            return

        if self._cache_dir.exists():
            shutil.rmtree(self._cache_dir)
            self._init_cache_dir()
            logger.info(f"Cleared all cached data from {self._cache_dir}")

    # -------------------------------------------------------------------------
    # Listing / Info
    # -------------------------------------------------------------------------

    def list_cached_datasets(self) -> list:
        """
        List all cached datasets.

        Returns
        -------
        list of dict
            Each dict has keys: repository, accession, n_genes, n_samples,
            cached_date, organism.
        """
        if not self.enabled or not self._cache_dir.exists():
            return []

        datasets = []
        for repo_name, subdir in REPOSITORY_SUBDIRS.items():
            repo_path = self._cache_dir / subdir
            if not repo_path.exists():
                continue
            for dataset_dir in sorted(repo_path.iterdir()):
                manifest_path = dataset_dir / 'manifest.json'
                if manifest_path.exists():
                    with open(manifest_path, 'r', encoding='utf-8') as f:
                        manifest = json.load(f)
                    datasets.append({
                        'repository': repo_name,
                        'accession': manifest.get('accession', dataset_dir.name),
                        'organism': manifest.get('organism', 'unknown'),
                        'n_genes': manifest.get('n_genes', 0),
                        'n_samples': manifest.get('n_samples', 0),
                        'data_type': manifest.get('data_type', 'unknown'),
                        'cached_date': manifest.get('cached_date', 'unknown'),
                    })

        return datasets

    def list_cached_pools(self) -> list:
        """
        List all cached pooled datasets.

        Returns
        -------
        list of dict
            Each dict has keys: pool_name, n_studies, n_genes, n_samples,
            cached_date.
        """
        if not self.enabled or not self._cache_dir.exists():
            return []

        pools = []
        pools_path = self._cache_dir / POOLS_SUBDIR
        if not pools_path.exists():
            return []

        for pool_dir in sorted(pools_path.iterdir()):
            manifest_path = pool_dir / 'pool_manifest.json'
            if manifest_path.exists():
                with open(manifest_path, 'r', encoding='utf-8') as f:
                    manifest = json.load(f)
                pools.append({
                    'pool_name': manifest.get('pool_name', pool_dir.name),
                    'n_studies': manifest.get('n_studies', 0),
                    'n_genes': manifest.get('n_genes', 0),
                    'n_samples': manifest.get('n_samples', 0),
                    'cached_date': manifest.get('cached_date', 'unknown'),
                })

        return pools

    def cache_summary(self) -> str:
        """Human-readable summary of cached data."""
        datasets = self.list_cached_datasets()
        pools = self.list_cached_pools()

        lines = [
            "",
            "=" * 62,
            "  RAPTOR Data Cache Summary",
            "=" * 62,
            "",
            f"  Cache directory: {self._cache_dir}",
            f"  Cache enabled:   {self.enabled}",
            "",
            f"  Datasets cached: {len(datasets)}",
        ]

        if datasets:
            for ds in datasets:
                lines.append(
                    f"    {ds['repository']}/{ds['accession']} "
                    f"({ds['n_genes']:,} genes, {ds['n_samples']} samples)"
                )

        lines.append(f"\n  Pools cached:    {len(pools)}")

        if pools:
            for pool in pools:
                lines.append(
                    f"    {pool['pool_name']} "
                    f"({pool['n_studies']} studies, {pool['n_samples']} samples)"
                )

        lines.extend(["", "=" * 62, ""])
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"CacheManager("
            f"dir='{self._cache_dir}', "
            f"enabled={self.enabled})"
        )
