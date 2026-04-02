"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - ArrayExpress/BioStudies Connector

Connects to EMBL-EBI BioStudies (formerly ArrayExpress) for downloading
RNA-seq datasets. Uses the BioStudies REST API.

API docs: https://www.ebi.ac.uk/biostudies/api

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import io
import gzip
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Union

import numpy as np
import pandas as pd

from .base import BaseConnector, SearchResult
from .datasets import AcquiredDataset

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

BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1"
BIOSTUDIES_SEARCH = f"{BIOSTUDIES_API}/search"
BIOSTUDIES_STUDIES = f"{BIOSTUDIES_API}/studies"

# ArrayExpress file download base
AE_FILES_BASE = "https://www.ebi.ac.uk/biostudies/files"


def _check_requests():
    """Check if requests is available."""
    try:
        import requests
        return requests
    except ImportError:
        raise RAPTORError(
            "requests is required for ArrayExpress API access. "
            "Install with: pip install requests"
        )


# =============================================================================
# ArrayExpConnector
# =============================================================================

class ArrayExpConnector(BaseConnector):
    """
    Connector for EMBL-EBI BioStudies (ArrayExpress).

    Searches and downloads RNA-seq datasets from the European
    Bioinformatics Institute repositories.

    Parameters
    ----------
    cache : bool
        Enable disk caching.
    cache_dir : str or Path, optional
        Custom cache directory.

    Examples
    --------
    >>> ae = ArrayExpConnector()
    >>> results = ae.search("breast cancer RNA-seq", organism="Homo sapiens")
    >>> dataset = ae.download("E-MTAB-1234")
    >>> print(dataset.summary())
    """

    REPOSITORY_NAME = 'ArrayExpress'

    def __init__(
        self,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
    ):
        super().__init__(
            cache=cache,
            cache_dir=cache_dir,
            rate_limit=0.25,
        )

    # -------------------------------------------------------------------------
    # Search
    # -------------------------------------------------------------------------

    def _search_api(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 20,
        **kwargs,
    ) -> List[SearchResult]:
        """Search BioStudies for RNA-seq datasets."""
        requests = _check_requests()

        params = {
            'query': f'{query} AND type:rnaseq',
            'pageSize': max_results,
            'page': 1,
        }

        if organism:
            params['query'] += f' AND organism:"{organism}"'

        try:
            resp = requests.get(BIOSTUDIES_SEARCH, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(f"BioStudies search failed: {e}")

        results = []
        for hit in data.get('hits', []):
            accession = hit.get('accession', '')
            title = hit.get('title', '')
            author = hit.get('author', '')

            # Extract organism and sample count from attributes
            attrs = {
                a.get('name', ''): a.get('value', '')
                for a in hit.get('attributes', [])
            }

            org = attrs.get('organism', 'unknown')
            n_samples = 0

            # Try to extract sample count
            links = hit.get('links', 0)
            if isinstance(links, int):
                n_samples = links

            results.append(SearchResult(
                accession=accession,
                title=title,
                organism=org,
                n_samples=n_samples,
                platform=attrs.get('technology', 'RNA-seq'),
                description=f"Author: {author}" if author else '',
                repository='ArrayExpress',
                release_date=hit.get('releaseDate', ''),
            ))

        return results

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    def _download_api(
        self,
        accession: str,
        data_type: str = 'raw_counts',
    ) -> AcquiredDataset:
        """Download a dataset from BioStudies/ArrayExpress."""
        requests = _check_requests()

        logger.info(f"Fetching study info for {accession}...")

        # Step 1: Get study metadata
        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                f"{BIOSTUDIES_STUDIES}/{accession}",
                timeout=30,
            )
            resp.raise_for_status()
            study_data = resp.json()
        except Exception as e:
            raise RAPTORError(f"Failed to fetch study {accession}: {e}")

        # Step 2: Find downloadable files
        file_list = self._find_files(study_data, accession)

        # Step 3: Download and parse count matrix
        counts_df, metadata = self._download_and_parse(
            accession, file_list, requests
        )

        # Extract study-level metadata
        study_attrs = {}
        for section in study_data.get('section', {}).get('attributes', []):
            name = section.get('name', '')
            value = section.get('value', '')
            if name and value:
                study_attrs[name.lower()] = value

        source_info = {
            'repository': 'ArrayExpress',
            'accession': accession,
            'organism': study_attrs.get('organism', 'unknown'),
            'platform': study_attrs.get('technology type', 'RNA-seq'),
            'download_date': datetime.now().isoformat(),
            'data_type': data_type,
            'description': study_attrs.get('title', ''),
        }

        gene_id_type = self._detect_gene_ids(counts_df)

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info=source_info,
            gene_id_type=gene_id_type,
        )

    def _find_files(self, study_data: dict, accession: str) -> List[dict]:
        """
        Find downloadable data files from study metadata.

        Looks for processed count matrices or SDRF files.
        """
        files = []

        def _walk_sections(section):
            """Recursively walk study sections for file references."""
            if isinstance(section, dict):
                # Check for file tables
                for f in section.get('files', []):
                    if isinstance(f, list):
                        for item in f:
                            if isinstance(item, dict) and item.get('path'):
                                files.append(item)
                    elif isinstance(f, dict) and f.get('path'):
                        files.append(f)

                # Recurse into subsections
                for subsection in section.get('subsections', []):
                    _walk_sections(subsection)

            elif isinstance(section, list):
                for item in section:
                    _walk_sections(item)

        _walk_sections(study_data.get('section', {}))

        return files

    def _download_and_parse(
        self,
        accession: str,
        file_list: List[dict],
        requests,
    ) -> tuple:
        """
        Download and parse data files.

        Strategies:
        1. Look for processed count matrix files
        2. Look for SDRF (sample metadata) files
        3. Fall back to raw data if needed
        """
        counts_df = None
        metadata = pd.DataFrame()

        # Prioritize count matrix files
        count_files = [
            f for f in file_list
            if any(kw in f.get('path', '').lower()
                   for kw in ['count', 'expression', 'matrix', 'quantification'])
            and f.get('path', '').endswith(('.tsv', '.csv', '.txt', '.tsv.gz', '.csv.gz'))
        ]

        # Also look for SDRF files
        sdrf_files = [
            f for f in file_list
            if 'sdrf' in f.get('path', '').lower()
            and f.get('path', '').endswith(('.tsv', '.txt'))
        ]

        # Try to download count matrix
        for cf in count_files:
            file_path = cf.get('path', '')
            url = f"{AE_FILES_BASE}/{accession}/{file_path}"

            logger.info(f"Downloading count file: {file_path}")
            self._wait_for_rate_limit()

            try:
                resp = requests.get(url, timeout=120)
                resp.raise_for_status()
                content = resp.content

                # Handle gzipped files
                if file_path.endswith('.gz'):
                    content = gzip.decompress(content)

                sep = '\t' if file_path.endswith(('.tsv', '.tsv.gz', '.txt')) else ','
                df = pd.read_csv(io.BytesIO(content), sep=sep, index_col=0)

                if df.shape[0] > 100 and df.shape[1] >= 2:
                    counts_df = df.apply(pd.to_numeric, errors='coerce')
                    counts_df = counts_df.dropna(how='all')
                    counts_df.index.name = 'gene_id'
                    logger.info(
                        f"Parsed count matrix: "
                        f"{counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples"
                    )
                    break

            except Exception as e:
                logger.debug(f"Failed to parse {file_path}: {e}")
                continue

        # Try to download SDRF
        for sf in sdrf_files:
            file_path = sf.get('path', '')
            url = f"{AE_FILES_BASE}/{accession}/{file_path}"

            self._wait_for_rate_limit()
            try:
                resp = requests.get(url, timeout=60)
                resp.raise_for_status()
                metadata = pd.read_csv(io.BytesIO(resp.content), sep='\t')

                # Use the first column as index (usually source/sample name)
                if len(metadata.columns) > 0:
                    metadata = metadata.set_index(metadata.columns[0])
                metadata.index.name = 'sample_id'
                break
            except Exception as e:
                logger.debug(f"Failed to parse SDRF {file_path}: {e}")

        if counts_df is None:
            raise RAPTORError(
                f"Could not find or parse a count matrix for {accession}. "
                f"Found {len(file_list)} files, {len(count_files)} candidate "
                f"count files. The study may not have processed data available."
            )

        return counts_df, metadata

    @staticmethod
    def _detect_gene_ids(counts_df: pd.DataFrame) -> str:
        """Detect gene ID type from index."""
        sample = counts_df.index[:100].astype(str).tolist()

        ensembl = sum(1 for g in sample if g.startswith('ENS'))
        if ensembl > len(sample) * 0.5:
            return 'ensembl'

        numeric = sum(1 for g in sample if g.isdigit())
        if numeric > len(sample) * 0.5:
            return 'entrez'

        return 'symbol'

    # -------------------------------------------------------------------------
    # Info
    # -------------------------------------------------------------------------

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """Get study metadata without downloading data."""
        requests = _check_requests()

        try:
            resp = requests.get(
                f"{BIOSTUDIES_STUDIES}/{accession}",
                timeout=15,
            )
            if resp.status_code == 200:
                data = resp.json()
                attrs = {}
                for a in data.get('section', {}).get('attributes', []):
                    attrs[a.get('name', '').lower()] = a.get('value', '')

                return {
                    'accession': accession,
                    'title': attrs.get('title', ''),
                    'organism': attrs.get('organism', 'unknown'),
                    'description': attrs.get('description', ''),
                    'repository': 'ArrayExpress',
                }
        except Exception:
            pass

        return None

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"ArrayExpConnector("
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )
