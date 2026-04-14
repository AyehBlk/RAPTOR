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

# ArrayExpress is a *collection* within BioStudies.
# The collection-scoped endpoint returns only ArrayExpress datasets,
# avoiding the millions of unrelated BioStudies entries (imaging, etc.).
BIOSTUDIES_AE_SEARCH = f"{BIOSTUDIES_API}/arrayexpress/search"

# Study detail endpoint (works for any accession, not collection-scoped)
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

    # Known organisms for parsing the `content` blob
    _KNOWN_ORGANISMS = [
        'Homo sapiens', 'Mus musculus', 'Rattus norvegicus',
        'Drosophila melanogaster', 'Danio rerio',
        'Caenorhabditis elegans', 'Saccharomyces cerevisiae',
        'Arabidopsis thaliana', 'Sus scrofa', 'Gallus gallus',
        'Bos taurus', 'Ovis aries', 'Macaca mulatta',
    ]

    # Known study types (experiment types from ArrayExpress/EFO)
    _KNOWN_STUDY_TYPES = [
        'RNA-seq of coding RNA',
        'RNA-seq of non coding RNA',
        'RNA-seq of total RNA',
        'microRNA profiling by array',
        'transcription profiling by array',
        'transcription profiling by high throughput sequencing',
        'methylation profiling by array',
        'methylation profiling by high throughput sequencing',
        'comparative genomic hybridization by array',
        'ChIP-seq',
        'ATAC-seq',
    ]

    # Map dashboard data-type selections to search terms and study-type filters
    _DATA_TYPE_QUERY = {
        'RNA-seq': '"RNA-seq"',
        'Microarray': '"transcription profiling by array"',
        'Any': '',
    }
    _DATA_TYPE_ACCEPT = {
        'RNA-seq': {
            'rna-seq of coding rna', 'rna-seq of non coding rna',
            'rna-seq of total rna',
            'transcription profiling by high throughput sequencing',
        },
        'Microarray': {
            'transcription profiling by array',
            'microrna profiling by array',
        },
        'Any': None,  # accept everything
    }

    def _search_api(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 20,
        **kwargs,
    ) -> List[SearchResult]:
        """Search BioStudies ArrayExpress collection for datasets."""
        requests = _check_requests()

        dataset_type = kwargs.get('dataset_type', 'Any')

        # Build query: add data-type keywords so BioStudies returns
        # relevant results (e.g. "breast cancer" + '"RNA-seq"')
        search_query = query
        type_query = self._DATA_TYPE_QUERY.get(dataset_type, '')
        if type_query:
            search_query = f'{query} {type_query}'

        # Request more results than needed so we can post-filter
        accepted_types = self._DATA_TYPE_ACCEPT.get(dataset_type)
        fetch_size = max_results * 3 if accepted_types else max_results

        params = {
            'query': search_query,
            'pageSize': min(fetch_size, 200),
            'page': 1,
        }

        # Organism filter
        if organism:
            params['query'] += f' AND organism:"{organism}"'

        try:
            resp = requests.get(
                BIOSTUDIES_AE_SEARCH, params=params, timeout=30,
                headers={'Accept': 'application/json'},
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(f"BioStudies ArrayExpress search failed: {e}")

        results = []
        for hit in data.get('hits', []):
            accession = hit.get('accession', '')
            title = hit.get('title', '')
            author = hit.get('author', '')

            # BioStudies search hits do NOT have an 'attributes' array.
            # Organism and study type are embedded in a flat 'content' blob.
            content = hit.get('content', '')

            org = self._parse_organism(content)
            study_type = self._parse_study_type(content)

            # Post-filter: skip studies that don't match the selected data type
            if accepted_types is not None:
                if study_type.lower() not in accepted_types:
                    continue

            # The search API does NOT return sample counts.
            # 'links' is a link count (not samples), 'files' is file count.
            # Real sample count is only in the study detail (too slow to fetch
            # per-result). Show 0; users can click Preview Info for real count.
            n_samples = 0

            # Platform = simplified technology category (not the full study type)
            platform = self._simplify_platform(study_type)

            results.append(SearchResult(
                accession=accession,
                title=title,
                organism=org,
                n_samples=n_samples,
                platform=platform,
                description=f"Author: {author}" if author else '',
                repository='ArrayExpress',
                release_date=hit.get('release_date', ''),
                gds_type=study_type,  # detailed type for "Type" column
            ))

            # Stop once we have enough after filtering
            if len(results) >= max_results:
                break

        return results

    @classmethod
    def _parse_organism(cls, content: str) -> str:
        """Extract organism name from the BioStudies content blob."""
        if not content:
            return 'unknown'
        for org in cls._KNOWN_ORGANISMS:
            if org in content:
                return org
        return 'unknown'

    @classmethod
    def _parse_study_type(cls, content: str) -> str:
        """Extract study/experiment type from the BioStudies content blob."""
        if not content:
            return ''
        content_lower = content.lower()
        for stype in cls._KNOWN_STUDY_TYPES:
            if stype.lower() in content_lower:
                return stype
        return ''

    @staticmethod
    def _simplify_platform(study_type: str) -> str:
        """Convert detailed study type to a short technology label for the Platform column."""
        if not study_type:
            return 'unknown'
        st = study_type.lower()
        if 'rna-seq' in st or 'high throughput sequencing' in st:
            return 'RNA-Sequencing'
        elif 'chip-seq' in st:
            return 'ChIP-Sequencing'
        elif 'atac-seq' in st:
            return 'ATAC-Sequencing'
        elif 'methylation' in st and 'array' in st:
            return 'Methylation Array'
        elif 'hybridization' in st:
            return 'CGH Array'
        elif 'array' in st:
            return 'Microarray'
        elif 'sequencing' in st:
            return 'Sequencing'
        return study_type

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
            # List found files for diagnostic
            found_names = [f.get('path', '') for f in file_list[:10]]
            found_str = ', '.join(found_names) if found_names else 'none'

            raise RAPTORError(
                f"Could not find or parse a count matrix for {accession}. "
                f"Found {len(file_list)} files ({found_str}), "
                f"{len(count_files)} candidate count files.\n\n"
                f"Many ArrayExpress RNA-seq studies store only raw FASTQs on "
                f"ENA (European Nucleotide Archive) without processed count "
                f"matrices. To use this data:\n"
                f"1. Check ENA for raw reads: "
                f"https://www.ebi.ac.uk/ena/browser/view/{accession}\n"
                f"2. Download FASTQs and quantify with RAPTOR Module 1 "
                f"(Salmon/Kallisto)\n"
                f"3. Or search for the linked GEO accession (E-GEOD-* studies "
                f"often have processed data on GEO)"
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
    # Info / Preview / Sample Types
    # -------------------------------------------------------------------------

    def _fetch_study_detail(self, accession: str) -> Optional[dict]:
        """Fetch the full study JSON from BioStudies. Cached per session."""
        requests = _check_requests()
        try:
            resp = requests.get(
                f"{BIOSTUDIES_STUDIES}/{accession}",
                headers={'Accept': 'application/json'},
                timeout=20,
            )
            if resp.status_code == 200:
                return resp.json()
        except Exception:
            pass
        return None

    @staticmethod
    def _parse_study_attrs(study_data: dict) -> Dict[str, str]:
        """Extract key-value attributes from a BioStudies study JSON."""
        attrs = {}
        for a in study_data.get('section', {}).get('attributes', []):
            name = a.get('name', '')
            value = a.get('value', '')
            if name and value:
                attrs[name.lower()] = value
        return attrs

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """Get study metadata without downloading data."""
        study_data = self._fetch_study_detail(accession)
        if study_data is None:
            return None

        attrs = self._parse_study_attrs(study_data)
        subsection_info = self._parse_subsections(study_data)

        return {
            'accession': accession,
            'title': attrs.get('title', ''),
            'organism': attrs.get('organism', 'unknown'),
            'description': attrs.get('description', ''),
            'repository': 'ArrayExpress',
            'platform_id': subsection_info.get('technology', attrs.get('technology type', 'unknown')),
            'data_type': subsection_info.get('experimental_designs', 'unknown'),
            'n_samples': subsection_info.get('sample_count', 0),
        }

    @staticmethod
    def _parse_subsections(study_data: dict) -> Dict[str, Any]:
        """
        Parse the typed subsection blocks in a BioStudies study JSON.

        ArrayExpress studies have subsections like:
          - "Samples" → Sample count, Experimental Designs, Experimental Factors
          - "Assays and Data" → Technology, Assay by Molecule
          - "Publication" → Title, Authors, DOI
          - "Organization" → Name, Address
        """
        result = {
            'sample_count': 0,
            'experimental_designs': '',
            'experimental_factors': [],
            'technology': '',
            'assay_molecule': '',
            'pub_title': '',
            'pub_authors': '',
            'pub_doi': '',
        }

        section = study_data.get('section', {})
        subsections = section.get('subsections', [])

        for sub in subsections:
            items = [sub] if isinstance(sub, dict) else (sub if isinstance(sub, list) else [])
            for item in items:
                if not isinstance(item, dict):
                    continue
                sec_type = item.get('type', '')
                item_attrs = {}
                for a in item.get('attributes', []):
                    name = a.get('name', '')
                    value = a.get('value', '')
                    if name:
                        # Some attributes repeat (e.g., multiple Experimental Factors)
                        if name in item_attrs:
                            if isinstance(item_attrs[name], list):
                                item_attrs[name].append(value)
                            else:
                                item_attrs[name] = [item_attrs[name], value]
                        else:
                            item_attrs[name] = value

                if sec_type == 'Samples':
                    # Parse "Sample count: 38"
                    sc = item_attrs.get('Sample count', '0')
                    try:
                        result['sample_count'] = int(sc)
                    except (ValueError, TypeError):
                        pass
                    result['experimental_designs'] = item_attrs.get('Experimental Designs', '')
                    ef = item_attrs.get('Experimental Factors', [])
                    if isinstance(ef, str):
                        ef = [ef]
                    result['experimental_factors'] = ef

                elif sec_type == 'Assays and Data':
                    result['technology'] = item_attrs.get('Technology', '')
                    result['assay_molecule'] = item_attrs.get('Assay by Molecule', '')

                elif sec_type == 'Publication':
                    result['pub_title'] = item_attrs.get('Title', '')
                    result['pub_authors'] = item_attrs.get('Authors', '')
                    result['pub_doi'] = item_attrs.get('DOI', '')

        return result

    def get_study_info(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get rich study metadata for the dashboard Preview panel.

        Returns a dict with keys the dashboard expects:
        title, organism, platform_id, platform_name, data_type,
        n_samples, summary, samples (DataFrame), supplementary_files, etc.
        """
        study_data = self._fetch_study_detail(accession)
        if study_data is None:
            return None

        attrs = self._parse_study_attrs(study_data)
        sub_info = self._parse_subsections(study_data)

        # Collect file paths from the study JSON
        files_list = []

        def _walk_files(section):
            if isinstance(section, dict):
                for f in section.get('files', []):
                    if isinstance(f, list):
                        for item in f:
                            if isinstance(item, dict) and item.get('path'):
                                files_list.append(item['path'])
                    elif isinstance(f, dict) and f.get('path'):
                        files_list.append(f['path'])
                for sub in section.get('subsections', []):
                    _walk_files(sub)
            elif isinstance(section, list):
                for item in section:
                    _walk_files(item)

        _walk_files(study_data.get('section', {}))

        # Build experimental factors display
        factors = sub_info.get('experimental_factors', [])
        factors_str = ', '.join(factors) if factors else ''

        # Technology / platform
        technology = sub_info.get('technology', '') or attrs.get('technology type', '')
        exp_design = sub_info.get('experimental_designs', '')

        # Release date from top-level attributes
        release_date = ''
        for a in study_data.get('attributes', []):
            if a.get('name', '').lower() == 'releasedate':
                release_date = a.get('value', '')

        info = {
            'accession': accession,
            'title': attrs.get('title', ''),
            'organism': attrs.get('organism', 'unknown'),
            'description': attrs.get('description', ''),
            'summary': attrs.get('description', ''),
            'platform_id': technology if technology else 'unknown',
            'platform_name': exp_design,
            'data_type': exp_design if exp_design else 'unknown',
            'n_samples': sub_info.get('sample_count', 0),
            'repository': 'ArrayExpress',
            'samples': pd.DataFrame(),  # individual samples require SDRF
            'supplementary_files': files_list[:50],
            'submission_date': release_date,
            'last_update': '',
            # ArrayExpress-specific extras
            'experimental_factors': factors_str,
            'pub_title': sub_info.get('pub_title', ''),
            'pub_authors': sub_info.get('pub_authors', ''),
            'pubmed_ids': [sub_info['pub_doi']] if sub_info.get('pub_doi', '').isdigit() else [],
        }

        return info

    def get_sample_types(self, accession: str) -> Dict[str, Any]:
        """
        Get sample type breakdown for the dashboard "Show Sample Types" button.

        For ArrayExpress, individual sample data lives in the SDRF file, not
        the study JSON. We attempt to download and parse it. If unavailable,
        we return the experimental factors from the study metadata.
        """
        requests = _check_requests()

        # First, try to fetch the SDRF for real sample-level info
        study_data = self._fetch_study_detail(accession)
        if study_data is None:
            return {'type_counts': {}, 'unique_types': []}

        # Find SDRF files
        sdrf_paths = []

        def _find_sdrf(section):
            if isinstance(section, dict):
                for f in section.get('files', []):
                    if isinstance(f, list):
                        for item in f:
                            if isinstance(item, dict) and 'sdrf' in item.get('path', '').lower():
                                sdrf_paths.append(item['path'])
                    elif isinstance(f, dict) and 'sdrf' in f.get('path', '').lower():
                        sdrf_paths.append(f['path'])
                for sub in section.get('subsections', []):
                    _find_sdrf(sub)
            elif isinstance(section, list):
                for item in section:
                    _find_sdrf(item)

        _find_sdrf(study_data.get('section', {}))

        type_counts = {}

        # Try to download and parse the first SDRF
        if sdrf_paths:
            sdrf_path = sdrf_paths[0]
            url = f"{AE_FILES_BASE}/{accession}/{sdrf_path}"
            try:
                self._wait_for_rate_limit()
                resp = requests.get(url, timeout=30)
                if resp.status_code == 200:
                    sdrf_df = pd.read_csv(io.BytesIO(resp.content), sep='\t')

                    # Look for characteristic or factor columns
                    for col in sdrf_df.columns:
                        col_lower = col.lower()
                        if any(kw in col_lower for kw in [
                            'characteristics[disease',
                            'characteristics[organism part',
                            'characteristics[cell type',
                            'factor',
                            'characteristics[phenotype',
                            'characteristics[tissue',
                        ]):
                            type_counts = sdrf_df[col].value_counts().to_dict()
                            break

                    # If no specific column found, try any Characteristics column
                    if not type_counts:
                        char_cols = [c for c in sdrf_df.columns
                                     if 'characteristic' in c.lower() or 'factor' in c.lower()]
                        # Pick the one with fewest unique values (most likely a grouping)
                        best_col = None
                        best_n = float('inf')
                        for col in char_cols:
                            n = sdrf_df[col].nunique()
                            if 2 <= n <= 15 and n < best_n:
                                best_col = col
                                best_n = n
                        if best_col:
                            type_counts = sdrf_df[best_col].value_counts().to_dict()

            except Exception as e:
                logger.debug(f"SDRF download failed for {accession}: {e}")

        # Fallback: return experimental factors from study metadata
        if not type_counts:
            sub_info = self._parse_subsections(study_data)
            factors = sub_info.get('experimental_factors', [])
            if factors:
                type_counts = {f"Factor: {f}": 0 for f in factors}

        return {
            'type_counts': type_counts,
            'unique_types': list(type_counts.keys()),
        }

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"ArrayExpConnector("
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )