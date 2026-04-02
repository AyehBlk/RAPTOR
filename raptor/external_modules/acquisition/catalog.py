"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - Data Catalog

Maintains a searchable index (catalog.json) of all acquired and pooled
datasets in the local cache. Used by the dashboard for browsing and by
Module 10 for listing available inputs.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Union

# RAPTOR imports
try:
    from raptor.utils.errors import ValidationError
except ImportError:
    class ValidationError(Exception):
        def __init__(self, parameter: str, message: str = "", hint: str = ""):
            self.parameter = parameter
            self.message = message
            self.hint = hint
            super().__init__(f"{parameter}: {message}. {hint}")


logger = logging.getLogger(__name__)


# =============================================================================
# DataCatalog
# =============================================================================

class DataCatalog:
    """
    Searchable index of all acquired and pooled datasets.

    The catalog is a lightweight JSON file that tracks what has been
    downloaded, when, and where it lives on disk. It avoids scanning
    the filesystem on every query.

    Parameters
    ----------
    cache_dir : str or Path
        Root cache directory (same as CacheManager's).
    auto_load : bool
        Load existing catalog.json on init if it exists.

    Examples
    --------
    >>> from raptor.external_modules.acquisition import DataCatalog
    >>> catalog = DataCatalog("~/.raptor/data")
    >>> catalog.search("liver cancer human")
    >>> catalog.list_datasets()
    >>> catalog.list_pools()
    """

    def __init__(
        self,
        cache_dir: Union[str, Path],
        auto_load: bool = True,
    ):
        self._cache_dir = Path(cache_dir).expanduser()
        self._catalog_path = self._cache_dir / 'catalog.json'
        self._datasets: List[Dict[str, Any]] = []
        self._pools: List[Dict[str, Any]] = []

        if auto_load and self._catalog_path.exists():
            self._load()

    # -------------------------------------------------------------------------
    # Persistence
    # -------------------------------------------------------------------------

    def _load(self) -> None:
        """Load catalog from disk."""
        try:
            with open(self._catalog_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            self._datasets = data.get('datasets', [])
            self._pools = data.get('pools', [])
            logger.debug(
                f"Catalog loaded: {len(self._datasets)} datasets, "
                f"{len(self._pools)} pools"
            )
        except (json.JSONDecodeError, KeyError) as e:
            logger.warning(f"Corrupt catalog.json, starting fresh: {e}")
            self._datasets = []
            self._pools = []

    def _save(self) -> None:
        """Write catalog to disk."""
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        data = {
            'version': '2.2.2',
            'updated': datetime.now().isoformat(),
            'datasets': self._datasets,
            'pools': self._pools,
        }
        with open(self._catalog_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, default=str)
        logger.debug("Catalog saved")

    # -------------------------------------------------------------------------
    # Register / Unregister
    # -------------------------------------------------------------------------

    def register_dataset(
        self,
        repository: str,
        accession: str,
        organism: str = 'unknown',
        data_type: str = 'raw_counts',
        n_genes: int = 0,
        n_samples: int = 0,
        platform: str = 'unknown',
        description: str = '',
        gene_id_type: str = 'symbol',
        tags: Optional[List[str]] = None,
        **extra,
    ) -> None:
        """
        Add or update a dataset entry in the catalog.

        Parameters
        ----------
        repository : str
            Source repository (GEO, TCGA, ArrayExpress, SRA).
        accession : str
            Dataset accession ID.
        organism : str
            Species name.
        data_type : str
            Type of data (raw_counts, normalized_counts, deg_table).
        n_genes : int
            Number of genes.
        n_samples : int
            Number of samples.
        platform : str
            Sequencing platform.
        description : str
            Brief dataset description.
        gene_id_type : str
            Gene identifier type.
        tags : list of str, optional
            User-defined tags for searching.
        **extra
            Any additional metadata fields.
        """
        # Remove existing entry if updating
        self._datasets = [
            d for d in self._datasets
            if not (d['repository'] == repository and d['accession'] == accession)
        ]

        entry = {
            'repository': repository,
            'accession': accession,
            'organism': organism,
            'data_type': data_type,
            'n_genes': n_genes,
            'n_samples': n_samples,
            'platform': platform,
            'description': description,
            'gene_id_type': gene_id_type,
            'tags': tags or [],
            'registered': datetime.now().isoformat(),
            **extra,
        }

        self._datasets.append(entry)
        self._save()
        logger.info(f"Registered {repository}/{accession} in catalog")

    def register_pool(
        self,
        pool_name: str,
        n_studies: int = 0,
        n_genes: int = 0,
        n_samples: int = 0,
        studies: Optional[List[str]] = None,
        method: str = 'unknown',
        batch_correction: str = 'none',
        description: str = '',
        tags: Optional[List[str]] = None,
        **extra,
    ) -> None:
        """
        Add or update a pool entry in the catalog.

        Parameters
        ----------
        pool_name : str
            Name of the pooled dataset.
        n_studies : int
            Number of studies merged.
        n_genes : int
            Number of genes in pooled matrix.
        n_samples : int
            Total samples across studies.
        studies : list of str, optional
            Accession IDs of component studies.
        method : str
            Pooling method used.
        batch_correction : str
            Batch correction method.
        description : str
            Brief pool description.
        tags : list of str, optional
            User-defined tags.
        **extra
            Additional metadata.
        """
        self._pools = [
            p for p in self._pools if p['pool_name'] != pool_name
        ]

        entry = {
            'pool_name': pool_name,
            'n_studies': n_studies,
            'n_genes': n_genes,
            'n_samples': n_samples,
            'studies': studies or [],
            'method': method,
            'batch_correction': batch_correction,
            'description': description,
            'tags': tags or [],
            'registered': datetime.now().isoformat(),
            **extra,
        }

        self._pools.append(entry)
        self._save()
        logger.info(f"Registered pool '{pool_name}' in catalog")

    def unregister_dataset(self, repository: str, accession: str) -> bool:
        """
        Remove a dataset from the catalog.

        Returns
        -------
        bool
            True if found and removed.
        """
        before = len(self._datasets)
        self._datasets = [
            d for d in self._datasets
            if not (d['repository'] == repository and d['accession'] == accession)
        ]
        removed = len(self._datasets) < before
        if removed:
            self._save()
            logger.info(f"Unregistered {repository}/{accession}")
        return removed

    def unregister_pool(self, pool_name: str) -> bool:
        """Remove a pool from the catalog."""
        before = len(self._pools)
        self._pools = [p for p in self._pools if p['pool_name'] != pool_name]
        removed = len(self._pools) < before
        if removed:
            self._save()
            logger.info(f"Unregistered pool '{pool_name}'")
        return removed

    # -------------------------------------------------------------------------
    # Query / Search
    # -------------------------------------------------------------------------

    def list_datasets(
        self,
        repository: Optional[str] = None,
        organism: Optional[str] = None,
        data_type: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        List datasets with optional filters.

        Parameters
        ----------
        repository : str, optional
            Filter by repository.
        organism : str, optional
            Filter by organism (case-insensitive partial match).
        data_type : str, optional
            Filter by data type.

        Returns
        -------
        list of dict
            Matching dataset entries.
        """
        results = self._datasets

        if repository:
            results = [d for d in results if d['repository'] == repository]
        if organism:
            org_lower = organism.lower()
            results = [
                d for d in results
                if org_lower in d.get('organism', '').lower()
            ]
        if data_type:
            results = [d for d in results if d.get('data_type') == data_type]

        return results

    def list_pools(self) -> List[Dict[str, Any]]:
        """List all pooled datasets."""
        return list(self._pools)

    def search(self, query: str) -> List[Dict[str, Any]]:
        """
        Free-text search across datasets and pools.

        Searches accession, organism, description, tags, and platform.
        Returns both datasets and pools with a 'type' field to distinguish.

        Parameters
        ----------
        query : str
            Search query (case-insensitive, space-separated terms
            are ANDed).

        Returns
        -------
        list of dict
            Matching entries with an added 'type' key ('dataset' or 'pool').
        """
        terms = query.lower().split()
        results = []

        for ds in self._datasets:
            searchable = ' '.join([
                ds.get('accession', ''),
                ds.get('organism', ''),
                ds.get('description', ''),
                ds.get('platform', ''),
                ds.get('repository', ''),
                ' '.join(ds.get('tags', [])),
            ]).lower()

            if all(term in searchable for term in terms):
                results.append({**ds, 'type': 'dataset'})

        for pool in self._pools:
            searchable = ' '.join([
                pool.get('pool_name', ''),
                pool.get('description', ''),
                ' '.join(pool.get('studies', [])),
                ' '.join(pool.get('tags', [])),
            ]).lower()

            if all(term in searchable for term in terms):
                results.append({**pool, 'type': 'pool'})

        return results

    def get_dataset(self, repository: str, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get a specific dataset entry.

        Parameters
        ----------
        repository : str
            Repository name.
        accession : str
            Dataset accession ID.

        Returns
        -------
        dict or None
        """
        for ds in self._datasets:
            if ds['repository'] == repository and ds['accession'] == accession:
                return ds
        return None

    def get_pool(self, pool_name: str) -> Optional[Dict[str, Any]]:
        """
        Get a specific pool entry.

        Parameters
        ----------
        pool_name : str
            Pool name.

        Returns
        -------
        dict or None
        """
        for pool in self._pools:
            if pool['pool_name'] == pool_name:
                return pool
        return None

    # -------------------------------------------------------------------------
    # Sync with filesystem
    # -------------------------------------------------------------------------

    def rebuild_from_cache(self, cache_manager) -> None:
        """
        Rebuild catalog by scanning the cache directory.

        Useful if catalog.json is lost or corrupted. Reads all
        manifest.json files from the cache.

        Parameters
        ----------
        cache_manager : CacheManager
            The cache manager to scan.
        """
        self._datasets = []
        self._pools = []

        for ds_info in cache_manager.list_cached_datasets():
            self.register_dataset(**ds_info)

        for pool_info in cache_manager.list_cached_pools():
            self.register_pool(**pool_info)

        logger.info(
            f"Rebuilt catalog: {len(self._datasets)} datasets, "
            f"{len(self._pools)} pools"
        )

    # -------------------------------------------------------------------------
    # Info
    # -------------------------------------------------------------------------

    @property
    def n_datasets(self) -> int:
        """Number of cataloged datasets."""
        return len(self._datasets)

    @property
    def n_pools(self) -> int:
        """Number of cataloged pools."""
        return len(self._pools)

    def summary(self) -> str:
        """Human-readable catalog summary."""
        lines = [
            "",
            "=" * 62,
            "  RAPTOR Data Catalog",
            "=" * 62,
            "",
            f"  Datasets: {self.n_datasets}",
        ]

        # Group by repository
        repos = {}
        for ds in self._datasets:
            repo = ds['repository']
            repos.setdefault(repo, []).append(ds)

        for repo, datasets in sorted(repos.items()):
            lines.append(f"    {repo}: {len(datasets)} datasets")
            for ds in datasets[:5]:
                lines.append(
                    f"      {ds['accession']} "
                    f"({ds.get('organism', '?')}, "
                    f"{ds.get('n_samples', '?')} samples)"
                )
            if len(datasets) > 5:
                lines.append(f"      ... and {len(datasets) - 5} more")

        lines.append(f"\n  Pools: {self.n_pools}")
        for pool in self._pools[:5]:
            lines.append(
                f"    {pool['pool_name']} "
                f"({pool.get('n_studies', '?')} studies, "
                f"{pool.get('n_samples', '?')} samples)"
            )

        lines.extend(["", "=" * 62, ""])
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"DataCatalog("
            f"datasets={self.n_datasets}, "
            f"pools={self.n_pools}, "
            f"path='{self._catalog_path}')"
        )
