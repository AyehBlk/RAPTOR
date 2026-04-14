"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - Base Connector

Abstract base class for all repository connectors. Defines the interface
(search, download, validate) and provides shared logic for caching,
rate limiting, and catalog registration.

Subclasses only need to implement the API-specific methods:
    _search_api(), _download_api(), _parse_metadata()

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

import pandas as pd

from .datasets import AcquiredDataset
from .cache import CacheManager, DEFAULT_CACHE_DIR
from .catalog import DataCatalog

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
# Search Result Container
# =============================================================================

class SearchResult:
    """
    A single search result from a repository query.

    Attributes
    ----------
    accession : str
        Dataset accession ID (e.g., GSE12345).
    title : str
        Dataset title.
    organism : str
        Species.
    n_samples : int
        Number of samples.
    platform : str
        Sequencing platform.
    description : str
        Brief description or abstract excerpt.
    repository : str
        Source repository.
    extra : dict
        Any additional repository-specific fields.
    """

    def __init__(
        self,
        accession: str,
        title: str = '',
        organism: str = 'unknown',
        n_samples: int = 0,
        platform: str = 'unknown',
        description: str = '',
        repository: str = 'unknown',
        **extra,
    ):
        self.accession = accession
        self.title = title
        self.organism = organism
        self.n_samples = n_samples
        self.platform = platform
        self.description = description
        self.repository = repository
        self.extra = extra

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'accession': self.accession,
            'title': self.title,
            'organism': self.organism,
            'n_samples': self.n_samples,
            'platform': self.platform,
            'description': self.description,
            'repository': self.repository,
            **self.extra,
        }

    def __repr__(self) -> str:
        return (
            f"SearchResult("
            f"accession='{self.accession}', "
            f"title='{self.title[:40]}...', "
            f"samples={self.n_samples})"
        )


# =============================================================================
# BaseConnector
# =============================================================================

class BaseConnector(ABC):
    """
    Abstract base class for repository connectors.

    Provides shared infrastructure for caching, catalog registration,
    and rate limiting. Subclasses implement repository-specific API
    calls.

    Parameters
    ----------
    cache : bool
        Enable disk caching. Default True.
    cache_dir : str or Path, optional
        Custom cache directory. Default ~/.raptor/data/
    rate_limit : float
        Minimum seconds between API requests. Default 0.34 (3 req/s).

    Examples
    --------
    Subclasses are used directly:

    >>> from raptor.external_modules.acquisition import GEOConnector
    >>> geo = GEOConnector()
    >>> results = geo.search("pancreatic cancer RNA-seq")
    >>> dataset = geo.download("GSE12345")
    """

    # Subclasses must set this
    REPOSITORY_NAME: str = 'unknown'

    def __init__(
        self,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
        rate_limit: float = 0.34,
    ):
        self._cache_manager = CacheManager(
            cache_dir=cache_dir,
            enabled=cache,
        )
        self._catalog = DataCatalog(
            cache_dir=cache_dir or DEFAULT_CACHE_DIR,
        )
        self._rate_limit = rate_limit
        self._last_request_time = 0.0

    # -------------------------------------------------------------------------
    # Rate limiting
    # -------------------------------------------------------------------------

    def _wait_for_rate_limit(self) -> None:
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._rate_limit:
            wait_time = self._rate_limit - elapsed
            logger.debug(f"Rate limit: waiting {wait_time:.2f}s")
            time.sleep(wait_time)
        self._last_request_time = time.time()

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def search(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 20,
        **kwargs,
    ) -> List[SearchResult]:
        """
        Search the repository for datasets.

        Parameters
        ----------
        query : str
            Search query (free text, keywords).
        organism : str, optional
            Filter by organism (e.g., 'Homo sapiens').
        max_results : int
            Maximum number of results to return.
        **kwargs
            Additional repository-specific options (e.g., dataset_type for GEO).

        Returns
        -------
        list of SearchResult
            Matching datasets.
        """
        logger.info(
            f"Searching {self.REPOSITORY_NAME} for: '{query}'"
            + (f" (organism={organism})" if organism else "")
        )

        self._wait_for_rate_limit()

        try:
            results = self._search_api(
                query=query,
                organism=organism,
                max_results=max_results,
                **kwargs,
            )
        except Exception as e:
            logger.error(f"Search failed: {e}")
            raise RAPTORError(f"{self.REPOSITORY_NAME} search failed: {e}")

        logger.info(f"Found {len(results)} results")
        return results

    def download(
        self,
        accession: str,
        data_type: str = 'raw_counts',
        force: bool = False,
        **kwargs,
    ) -> AcquiredDataset:
        """
        Download a dataset from the repository.

        Checks the local cache first. If cached and force=False,
        returns the cached version without hitting the API.

        Parameters
        ----------
        accession : str
            Dataset accession ID.
        data_type : str
            What to download: 'raw_counts', 'normalized_counts', 'deg_table'.
        force : bool
            Re-download even if cached.

        Returns
        -------
        AcquiredDataset
            The downloaded (or cached) dataset.
        """
        # Check cache first
        if not force and self._cache_manager.is_cached(self.REPOSITORY_NAME, accession):
            logger.info(f"{accession} found in cache, loading...")
            cached = self._cache_manager.load_dataset(
                self.REPOSITORY_NAME, accession
            )
            if cached is not None:
                return AcquiredDataset(
                    counts_df=cached['counts_df'],
                    metadata=cached['metadata'],
                    source_info=cached['source_info'],
                    gene_id_type=cached['gene_id_type'],
                )

        # Download from API
        logger.info(f"Downloading {accession} from {self.REPOSITORY_NAME}...")
        self._wait_for_rate_limit()

        try:
            dataset = self._download_api(
                accession=accession,
                data_type=data_type,
                **kwargs,
            )
        except Exception as e:
            logger.error(f"Download failed for {accession}: {e}")
            raise RAPTORError(
                f"Failed to download {accession} from "
                f"{self.REPOSITORY_NAME}: {e}"
            )

        # Validate the downloaded dataset
        integrity = dataset.validate_integrity()
        if not integrity['valid']:
            for error in integrity['errors']:
                logger.warning(f"Integrity issue: {error}")
            # Don't fail — real-world data often has mismatches.
            # Log warnings and continue with what we have.
            logger.warning(
                f"{accession}: integrity check found issues, "
                f"but data was downloaded. Proceeding with available data."
            )

        for warning in integrity['warnings']:
            logger.warning(f"Integrity warning: {warning}")

        # Cache the dataset
        self._cache_manager.save_dataset(
            counts_df=dataset.counts_df,
            metadata=dataset.metadata,
            source_info=dataset.source_info,
            gene_id_type=dataset.gene_id_type,
        )

        # Register in catalog
        self._catalog.register_dataset(
            repository=self.REPOSITORY_NAME,
            accession=accession,
            organism=dataset.organism,
            data_type=dataset.data_type,
            n_genes=dataset.n_genes,
            n_samples=dataset.n_samples,
            gene_id_type=dataset.gene_id_type,
            description=dataset.source_info.get('description', ''),
            platform=dataset.source_info.get('platform', 'unknown'),
        )

        logger.info(
            f"Downloaded {accession}: "
            f"{dataset.n_genes:,} genes, {dataset.n_samples} samples"
        )
        return dataset

    def download_multiple(
        self,
        accessions: List[str],
        data_type: str = 'raw_counts',
        force: bool = False,
    ) -> List[AcquiredDataset]:
        """
        Download multiple datasets.

        Parameters
        ----------
        accessions : list of str
            Dataset accession IDs.
        data_type : str
            Data type to download.
        force : bool
            Re-download even if cached.

        Returns
        -------
        list of AcquiredDataset
            Successfully downloaded datasets. Failures are logged
            and skipped.
        """
        datasets = []
        for i, accession in enumerate(accessions, 1):
            logger.info(f"[{i}/{len(accessions)}] {accession}")
            try:
                ds = self.download(accession, data_type=data_type, force=force)
                datasets.append(ds)
            except Exception as e:
                logger.error(f"Skipping {accession}: {e}")
                continue

        logger.info(
            f"Downloaded {len(datasets)}/{len(accessions)} datasets "
            f"successfully"
        )
        return datasets

    def info(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get metadata about a dataset without downloading it.

        Parameters
        ----------
        accession : str
            Dataset accession ID.

        Returns
        -------
        dict or None
            Dataset metadata, or None if not found.
        """
        # Check catalog first
        entry = self._catalog.get_dataset(self.REPOSITORY_NAME, accession)
        if entry:
            return entry

        # Query API
        self._wait_for_rate_limit()
        try:
            return self._info_api(accession)
        except Exception as e:
            logger.warning(f"Could not get info for {accession}: {e}")
            return None

    # -------------------------------------------------------------------------
    # Abstract methods (subclasses must implement)
    # -------------------------------------------------------------------------

    @abstractmethod
    def _search_api(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 20,
    ) -> List[SearchResult]:
        """
        Execute a search against the repository API.

        Parameters
        ----------
        query : str
            Search query.
        organism : str, optional
            Organism filter.
        max_results : int
            Max results.

        Returns
        -------
        list of SearchResult
        """
        ...

    @abstractmethod
    def _download_api(
        self,
        accession: str,
        data_type: str = 'raw_counts',
    ) -> AcquiredDataset:
        """
        Download a dataset from the repository API.

        Parameters
        ----------
        accession : str
            Accession ID.
        data_type : str
            Type of data to download.

        Returns
        -------
        AcquiredDataset
        """
        ...

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get metadata from the repository API.

        Optional override. Default returns None.

        Parameters
        ----------
        accession : str
            Accession ID.

        Returns
        -------
        dict or None
        """
        return None

    # -------------------------------------------------------------------------
    # Cache / Catalog access
    # -------------------------------------------------------------------------

    @property
    def cache(self) -> CacheManager:
        """Access the cache manager."""
        return self._cache_manager

    @property
    def catalog(self) -> DataCatalog:
        """Access the data catalog."""
        return self._catalog

    def list_cached(self) -> List[Dict[str, Any]]:
        """List datasets from this repository that are cached."""
        return self._catalog.list_datasets(repository=self.REPOSITORY_NAME)

    def is_cached(self, accession: str) -> bool:
        """Check if a specific dataset is cached."""
        return self._cache_manager.is_cached(self.REPOSITORY_NAME, accession)

    # -------------------------------------------------------------------------
    # Display
    # -------------------------------------------------------------------------

    @staticmethod
    def format_search_results(results: List[SearchResult]) -> str:
        """
        Format search results as a readable table.

        Parameters
        ----------
        results : list of SearchResult
            Search results to format.

        Returns
        -------
        str
            Formatted text table.
        """
        if not results:
            return "  No results found."

        lines = [
            "",
            f"  {'#':<4} {'Accession':<15} {'Samples':<9} {'Organism':<20} Title",
            "  " + "-" * 75,
        ]

        for i, r in enumerate(results, 1):
            title = r.title[:35] + '...' if len(r.title) > 38 else r.title
            organism = r.organism[:18] if len(r.organism) > 18 else r.organism
            lines.append(
                f"  {i:<4} {r.accession:<15} {r.n_samples:<9} "
                f"{organism:<20} {title}"
            )

        lines.append("")
        return "\n".join(lines)

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"{self.__class__.__name__}("
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )
