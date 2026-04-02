"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - GEO Connector

Connects to NCBI Gene Expression Omnibus (GEO) for searching,
downloading, and parsing RNA-seq datasets. Uses Entrez for search
and GEOparse for data retrieval.

Dependencies:
    pip install GEOparse biopython

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import os
import io
import gzip
import tempfile
import warnings
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


def _list_to_str(val, sep='; ') -> str:
    """Convert a GEO metadata value (str or list) to a string."""
    if isinstance(val, list):
        return sep.join(str(v) for v in val) if val else ''
    return str(val) if val else ''


# Common GPL platform lookup
GPL_PLATFORMS = {
    # Illumina RNA-seq — Mus musculus
    '24247': 'NovaSeq 6000',
    '21103': 'HiSeq 4000',
    '19057': 'NextSeq 500',
    '17021': 'HiSeq 2500',
    '13112': 'HiSeq 2000',
    '21290': 'HiSeq 3000',
    '34290': 'NovaSeq 6000',
    # Illumina RNA-seq — Homo sapiens
    '24676': 'NovaSeq 6000',
    '20301': 'HiSeq 4000',
    '18573': 'NextSeq 500',
    '16791': 'HiSeq 2500',
    '11154': 'HiSeq 2000',
    '15520': 'HiSeq 2000',
    '21493': 'HiSeq 2500',
    '23227': 'NovaSeq 6000',
    '25947': 'NovaSeq X',
    '26530': 'NovaSeq X Plus',
    '28913': 'NovaSeq X Plus',
    '19415': 'NextSeq 500',
    '20795': 'HiSeq X Ten',
    '21697': 'MiSeq',
    # Other organisms
    '21626': 'NextSeq 500',
    '22396': 'HiSeq 2500',
    '24718': 'NovaSeq 6000',
    '30173': 'NovaSeq 6000',
    # Microarray
    '570': 'Affy HG-U133+2',
    '571': 'Affy HG-U133A 2.0',
    '96': 'Affy HG-U133A',
    '1261': 'Agilent Mouse G4121A',
    '6885': 'Illumina MouseRef-8 v2',
    '6887': 'Illumina MouseWG-6 v2',
    '10558': 'Illumina HumanHT-12 v4',
    '6244': 'Affy HG-U133+2',
    '6246': 'Affy MoGene-1_0-st',
    '6480': 'Agilent Human GE 8x60K',
    '13912': 'Agilent Human GE v2 8x60K',
}


def _gpl_to_name(gpl_str: str) -> str:
    """Convert GPL number(s) to human-readable platform name(s)."""
    if not gpl_str or gpl_str == 'unknown':
        return 'unknown'

    # Handle multiple platforms (e.g., "24247;19057")
    gpls = [g.strip() for g in str(gpl_str).replace('GPL', '').split(';') if g.strip()]

    names = []
    for gpl in gpls:
        name = GPL_PLATFORMS.get(gpl, f'GPL{gpl}')
        names.append(name)

    return '; '.join(names)


# =============================================================================
# Dependency checks
# =============================================================================

def _check_geoparse():
    """Check if GEOparse is installed."""
    try:
        import GEOparse
        return GEOparse
    except ImportError:
        raise RAPTORError(
            "GEOparse is required for GEO data access. "
            "Install with: pip install GEOparse"
        )


def _check_entrez():
    """Check if Biopython Entrez is available."""
    try:
        from Bio import Entrez
        return Entrez
    except ImportError:
        raise RAPTORError(
            "Biopython is required for GEO search. "
            "Install with: pip install biopython"
        )


# =============================================================================
# GEOConnector
# =============================================================================

class GEOConnector(BaseConnector):
    """
    Connector for NCBI Gene Expression Omnibus (GEO).

    Searches GEO via NCBI Entrez and downloads datasets using GEOparse.
    Supports GSE (series) accessions and extracts count matrices and
    sample metadata.

    Parameters
    ----------
    email : str, optional
        Email for NCBI Entrez (required by NCBI policy).
        Defaults to NCBI_EMAIL environment variable.
    api_key : str, optional
        NCBI API key for higher rate limits (10 req/s vs 3 req/s).
        Defaults to NCBI_API_KEY environment variable.
    cache : bool
        Enable disk caching.
    cache_dir : str or Path, optional
        Custom cache directory.

    Examples
    --------
    >>> geo = GEOConnector(email="user@example.com")
    >>> results = geo.search("pancreatic cancer RNA-seq", organism="Homo sapiens")
    >>> for r in results[:5]:
    ...     print(f"{r.accession}: {r.title} ({r.n_samples} samples)")
    >>> dataset = geo.download("GSE12345")
    >>> print(dataset.summary())
    """

    REPOSITORY_NAME = 'GEO'

    def __init__(
        self,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
    ):
        # NCBI requires email; use env var as fallback
        self._email = email or os.environ.get('NCBI_EMAIL', '')
        self._api_key = api_key or os.environ.get('NCBI_API_KEY', None)

        # Rate limit: 3/s without API key, 10/s with
        rate_limit = 0.1 if self._api_key else 0.34

        super().__init__(
            cache=cache,
            cache_dir=cache_dir,
            rate_limit=rate_limit,
        )

        if not self._email:
            logger.warning(
                "No email set for NCBI Entrez. Set via email parameter "
                "or NCBI_EMAIL environment variable. NCBI requires this "
                "for API access."
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
        """Search GEO via NCBI Entrez."""
        Entrez = _check_entrez()
        Entrez.email = self._email
        if self._api_key:
            Entrez.api_key = self._api_key

        # Build search query
        search_term = f"{query}[Title/Abstract]"
        if organism:
            search_term += f" AND {organism}[Organism]"

        # Data type filter (default: RNA-seq only)
        dataset_type = kwargs.get('dataset_type', 'RNA-seq')
        TYPE_FILTERS = {
            'RNA-seq': '"expression profiling by high throughput sequencing"[DataSet Type]',
            'Microarray': '"expression profiling by array"[DataSet Type]',
            'Methylation': '"methylation profiling by high throughput sequencing"[DataSet Type]',
            'All expression': (
                '("expression profiling by high throughput sequencing"[DataSet Type] '
                'OR "expression profiling by array"[DataSet Type])'
            ),
            'Any': '',
        }
        type_filter = TYPE_FILTERS.get(dataset_type, TYPE_FILTERS['RNA-seq'])
        if type_filter:
            search_term += f' AND {type_filter}'

        logger.debug(f"Entrez query: {search_term}")

        # Search GDS database
        try:
            handle = Entrez.esearch(
                db='gds',
                term=search_term,
                retmax=max_results,
                sort='relevance',
            )
            search_results = Entrez.read(handle)
            handle.close()
        except Exception as e:
            raise RAPTORError(f"GEO search failed: {e}")

        id_list = search_results.get('IdList', [])
        if not id_list:
            logger.info("No results found")
            return []

        # Fetch summaries
        self._wait_for_rate_limit()
        try:
            handle = Entrez.esummary(db='gds', id=','.join(id_list))
            summaries = Entrez.read(handle)
            handle.close()
        except Exception as e:
            raise RAPTORError(f"Failed to fetch GEO summaries: {e}")

        results = []
        for summary in summaries:
            # Skip non-GSE entries
            accession = summary.get('Accession', '')
            if not accession.startswith('GSE'):
                continue

            # Extract experiment type
            gds_type = summary.get('gdsType', '')
            type_display = gds_type
            if 'high throughput sequencing' in gds_type.lower():
                type_display = 'RNA-seq'
            elif 'array' in gds_type.lower():
                type_display = 'Microarray'
            elif 'methylation' in gds_type.lower():
                type_display = 'Methylation'

            # Extract date (PDAT = publication date)
            pdat = str(summary.get('PDAT', ''))
            year = pdat[:4] if len(pdat) >= 4 else ''

            # PubMed IDs
            pubmed_ids = summary.get('PubMedIds', [])
            if isinstance(pubmed_ids, list):
                pubmed_ids = [str(p) for p in pubmed_ids if p]

            # FTP link indicates downloadable data
            ftp_link = summary.get('FTPLink', '')
            has_data = bool(ftp_link)

            # Extract a short design hint from the summary
            full_summary = summary.get('summary', '')
            # Take first sentence (up to 150 chars)
            design_hint = full_summary[:150]
            if '.' in design_hint:
                design_hint = design_hint[:design_hint.index('.') + 1]

            # Platform name lookup
            gpl_raw = summary.get('GPL', 'unknown')
            platform_name = _gpl_to_name(str(gpl_raw))

            results.append(SearchResult(
                accession=accession,
                title=summary.get('title', ''),
                organism=summary.get('taxon', 'unknown'),
                n_samples=int(summary.get('n_samples', 0)),
                platform=platform_name,
                description=full_summary[:300],
                repository='GEO',
                gds_type=type_display,
                gpl_id=str(gpl_raw),
                pubmed_ids=pubmed_ids,
                has_data=has_data,
                ftp_link=ftp_link,
                design_hint=design_hint,
            ))

        # Resolve unknown GPL IDs to platform names dynamically
        results = self._resolve_gpl_names(results)

        return results

    def _resolve_gpl_names(self, results: List[SearchResult]) -> List[SearchResult]:
        """
        Resolve GPL IDs that aren't in the static lookup by querying Entrez.
        One batch API call for all unknown GPLs.
        """
        Entrez = _check_entrez()

        # Collect all unknown GPL IDs
        unknown_gpls = set()
        for r in results:
            gpl_id = r.extra.get('gpl_id', '')
            for gpl in str(gpl_id).split(';'):
                gpl = gpl.strip()
                if gpl and gpl not in GPL_PLATFORMS:
                    unknown_gpls.add(gpl)

        if not unknown_gpls:
            return results

        # Batch fetch GPL info from Entrez
        gpl_names = {}
        try:
            self._wait_for_rate_limit()
            # Search for GPLs in the geo database
            gpl_terms = ' OR '.join(f'GPL{g}[Accession]' for g in unknown_gpls)
            handle = Entrez.esearch(
                db='geo',
                term=f'({gpl_terms}) AND GPL[EntryType]',
                retmax=len(unknown_gpls) + 5,
            )
            result = Entrez.read(handle)
            handle.close()

            geo_ids = result.get('IdList', [])
            if geo_ids:
                self._wait_for_rate_limit()
                handle = Entrez.esummary(db='geo', id=','.join(geo_ids))
                summaries = Entrez.read(handle)
                handle.close()

                for s in summaries:
                    acc = s.get('Accession', '')
                    title = s.get('title', '')
                    if acc.startswith('GPL') and title:
                        gpl_num = acc.replace('GPL', '')
                        # Clean up: "Illumina NovaSeq 6000 (Mus musculus)" → "NovaSeq 6000"
                        clean_name = title
                        # Remove manufacturer prefix if present
                        for prefix in ['Illumina ', 'Affymetrix ', 'Agilent ']:
                            if clean_name.startswith(prefix):
                                clean_name = clean_name[len(prefix):]
                        gpl_names[gpl_num] = clean_name
                        # Also cache in module-level dict for future use
                        GPL_PLATFORMS[gpl_num] = clean_name

        except Exception as e:
            logger.debug(f"GPL resolution failed: {e}")

        # Update results with resolved names
        if gpl_names:
            for r in results:
                gpl_id = r.extra.get('gpl_id', '')
                r.platform = _gpl_to_name(str(gpl_id))

        return results

    # -------------------------------------------------------------------------
    # Quick sample type lookup (fast, no SOFT download)
    # -------------------------------------------------------------------------

    def get_sample_types(self, accession: str) -> Dict[str, Any]:
        """
        Fetch sample titles for a GSE to show unique conditions/groups.

        Tries Entrez first (fast), falls back to GEOparse (reliable).

        Parameters
        ----------
        accession : str
            GEO accession (GSE...).

        Returns
        -------
        dict
            Keys: unique_types (list), type_counts (dict), n_samples (int)
        """
        import re

        all_titles = []

        # Try Entrez first (fast)
        try:
            Entrez = _check_entrez()
            Entrez.email = self._email
            if self._api_key:
                Entrez.api_key = self._api_key

            handle = Entrez.esearch(
                db='geo',
                term=f'{accession}[Accession]',
                retmax=500,
            )
            result = Entrez.read(handle)
            handle.close()

            gsm_ids = result.get('IdList', [])
            if gsm_ids:
                self._wait_for_rate_limit()
                handle = Entrez.esummary(db='geo', id=','.join(gsm_ids[:200]))
                summaries = Entrez.read(handle)
                handle.close()

                for s in summaries:
                    title = s.get('title', '')
                    # Only include GSM entries (not the GSE itself)
                    etype = s.get('entryType', '')
                    if title and (etype == 'GSM' or not etype):
                        all_titles.append(title)

        except Exception as e:
            logger.debug(f"Entrez sample lookup failed: {e}")

        # Fallback: use GEOparse (always works)
        if not all_titles:
            try:
                logger.info(f"Falling back to GEOparse for {accession} sample types...")
                info = self.get_study_info(accession)
                samples_df = info.get('samples')
                if samples_df is not None and 'Title' in samples_df.columns:
                    all_titles = samples_df['Title'].tolist()
            except Exception as e:
                logger.warning(f"GEOparse fallback also failed: {e}")

        if not all_titles:
            return {'unique_types': [], 'type_counts': {}, 'n_samples': 0}

        # Extract unique sample types by stripping replicate suffixes
        type_map = {}
        for title in all_titles:
            clean = re.sub(
                r'[,_\s]*(rep|replicate|sample|rep\.?\s*)\s*\d+\s*$',
                '', title, flags=re.IGNORECASE
            ).strip()
            clean = re.sub(r'[_\s]+\d+\s*$', '', clean).strip()
            if not clean:
                clean = title
            type_map.setdefault(clean, 0)
            type_map[clean] += 1

        sorted_types = sorted(type_map.items(), key=lambda x: -x[1])

        return {
            'unique_types': [t[0] for t in sorted_types],
            'type_counts': dict(sorted_types),
            'n_samples': len(all_titles),
        }

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    def _download_api(
        self,
        accession: str,
        data_type: str = 'raw_counts',
    ) -> AcquiredDataset:
        """Download a GSE dataset using GEOparse."""
        GEOparse = _check_geoparse()

        if not accession.startswith('GSE'):
            raise ValidationError('accession', f"Invalid GEO accession: {accession}", "GEOConnector supports GSE accessions (e.g., GSE12345)"
            )

        # Download to temp directory (GEOparse manages its own download)
        with tempfile.TemporaryDirectory() as tmpdir:
            logger.info(f"Downloading {accession} via GEOparse...")

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    gse = GEOparse.get_GEO(
                        geo=accession,
                        destdir=tmpdir,
                        silent=True,
                    )
                except Exception as e:
                    raise RAPTORError(
                        f"GEOparse failed to download {accession}: {e}"
                    )

            # Extract count matrix
            counts_df = self._extract_counts(gse, data_type)

            # Extract metadata from GSM objects
            metadata = self._extract_metadata(gse)

            # Align metadata with count matrix columns
            metadata = self._align_metadata(counts_df, metadata, gse)

            # Build source info
            source_info = {
                'repository': 'GEO',
                'accession': accession,
                'organism': self._extract_organism(gse),
                'platform': self._extract_platform(gse),
                'download_date': datetime.now().isoformat(),
                'data_type': data_type,
                'description': getattr(gse.metadata, 'get',
                    lambda k, d: gse.metadata.get(k, [d]) if isinstance(gse.metadata, dict) else d
                )('title', [''])[0] if hasattr(gse, 'metadata') else '',
                'n_supplementary_files': self._count_supplementary(gse),
            }

            # Detect gene ID type
            gene_id_type = self._detect_gene_id_type(counts_df)

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info=source_info,
            gene_id_type=gene_id_type,
        )

    # -------------------------------------------------------------------------
    # Data extraction helpers
    # -------------------------------------------------------------------------

    def _extract_counts(
        self,
        gse,
        data_type: str,
    ) -> pd.DataFrame:
        """
        Extract count matrix from a GSE object.

        Tries multiple strategies:
        1. Build from GSM tables
        2. GEOparse pivot table
        3. Merge individual sample tables
        4. Download supplementary count files (most reliable for RNA-seq)

        After extraction, auto-detects if matrix is transposed and fixes it.
        """
        # Strategy 1: Try to build from GSM tables
        counts = self._counts_from_gsms(gse)
        if counts is not None and not counts.empty:
            logger.info(
                f"Extracted counts from GSMs: "
                f"{counts.shape[0]:,} genes x {counts.shape[1]} samples"
            )
            return self._orient_counts(counts, gse)

        # Strategy 2: Try the pivot table method
        try:
            pivot = gse.pivot_samples('VALUE')
            if pivot is not None and not pivot.empty:
                pivot.index.name = 'gene_id'
                logger.info(
                    f"Extracted via pivot: "
                    f"{pivot.shape[0]:,} genes x {pivot.shape[1]} samples"
                )
                return self._orient_counts(pivot, gse)
        except Exception as e:
            logger.debug(f"Pivot method failed: {e}")

        # Strategy 3: Merge individual sample tables
        counts = self._counts_from_merge(gse)
        if counts is not None and not counts.empty:
            return self._orient_counts(counts, gse)

        # Strategy 4: Download supplementary count files
        accession = gse.name if hasattr(gse, 'name') else 'unknown'
        logger.info(f"Trying supplementary file download for {accession}...")
        counts = self._counts_from_supplementary(gse)
        if counts is not None and not counts.empty:
            logger.info(
                f"Extracted from supplementary files: "
                f"{counts.shape[0]:,} genes x {counts.shape[1]} samples"
            )
            return self._orient_counts(counts, gse)

        raise RAPTORError(
            f"Could not extract count matrix from {accession}. "
            f"Tried GSM tables, pivot, merge, and supplementary files. "
            f"This dataset may store data in a format RAPTOR cannot parse automatically. "
            f"Try downloading supplementary files from GEO manually."
        )

    def _orient_counts(self, counts: pd.DataFrame, gse) -> pd.DataFrame:
        """
        Auto-detect if a count matrix is transposed and fix it.

        Expected orientation: genes (rows) x samples (columns).
        Many GEO supplementary files store data as samples/cells x genes,
        which needs to be transposed.
        """
        n_rows, n_cols = counts.shape
        expected_samples = len(gse.gsms) if gse.gsms else 0

        # Heuristic 1: Compare dimensions with expected sample count
        if expected_samples > 0:
            row_match = abs(n_rows - expected_samples) <= max(2, expected_samples * 0.2)
            col_match = abs(n_cols - expected_samples) <= max(2, expected_samples * 0.2)

            if col_match and not row_match:
                return counts  # correct orientation
            if row_match and not col_match:
                logger.info(
                    f"Transposing: {n_rows} x {n_cols} -> {n_cols} x {n_rows} "
                    f"(rows match expected {expected_samples} samples)"
                )
                counts = counts.T
                counts.index.name = 'gene_id'
                return counts

        # Heuristic 2: Check if COLUMNS look like gene IDs
        # (regardless of dimension ratio — catches single-cell data)
        col_gene_score = self._gene_id_score(
            [str(c) for c in counts.columns[:200]]
        )
        row_gene_score = self._gene_id_score(
            [str(r) for r in counts.index[:200]]
        )

        if col_gene_score > 0.4 and row_gene_score < 0.2:
            # Columns are genes, rows are not → need to transpose
            logger.info(
                f"Transposing: {n_rows} x {n_cols} -> {n_cols} x {n_rows} "
                f"(columns are gene IDs, rows are not)"
            )
            # Warn if this looks like single-cell data
            if n_rows > 500 and n_cols > 10000:
                logger.warning(
                    f"This dataset may be single-cell RNA-seq "
                    f"({n_rows:,} barcodes x {n_cols:,} genes). "
                    f"RAPTOR is designed for bulk RNA-seq."
                )
            counts = counts.T
            counts.index.name = 'gene_id'
            return counts

        # Heuristic 3: Row names are GSM accessions
        row_sample = [str(r) for r in counts.index[:20]]
        gsm_rows = sum(1 for r in row_sample if r.startswith('GSM'))
        if gsm_rows > len(row_sample) * 0.5:
            logger.info(
                f"Transposing: {n_rows} x {n_cols} -> {n_cols} x {n_rows} "
                f"(rows are GSM accessions)"
            )
            counts = counts.T
            counts.index.name = 'gene_id'
            return counts

        return counts

    @staticmethod
    def _gene_id_score(names: list) -> float:
        """
        Score how likely a list of names are gene identifiers.

        Returns 0.0 (not genes) to 1.0 (definitely genes).
        Checks for Ensembl IDs, gene symbols, and Entrez IDs.
        """
        if not names:
            return 0.0

        import re
        gene_count = 0
        for name in names:
            name = str(name).strip()
            # Ensembl IDs
            if re.match(r'^ENS[A-Z]*G\d{11}', name):
                gene_count += 1
            # Entrez numeric IDs (pure digits, 1-9 digits)
            elif name.isdigit() and 1 <= len(name) <= 9:
                gene_count += 1
            # Gene symbols: alphabetic or alphanumeric, 2-15 chars
            # e.g., TP53, BRCA1, Xkr4, Gm1992, Sox17
            elif re.match(r'^[A-Za-z][A-Za-z0-9]{1,14}$', name):
                gene_count += 1

        return gene_count / len(names)

    def _counts_from_gsms(self, gse) -> Optional[pd.DataFrame]:
        """Build count matrix from individual GSM tables."""
        gsms = gse.gsms
        if not gsms:
            return None

        sample_dfs = {}
        for gsm_name, gsm in gsms.items():
            table = gsm.table
            if table is None or table.empty:
                continue

            # Find the gene ID column
            id_col = None
            for candidate in ['ID_REF', 'Gene', 'gene_id', 'GeneID']:
                if candidate in table.columns:
                    id_col = candidate
                    break

            if id_col is None and len(table.columns) >= 1:
                id_col = table.columns[0]

            # Find the value column
            val_col = None
            for candidate in ['VALUE', 'value', 'counts', 'Count',
                              'expected_count', 'TPM', 'FPKM']:
                if candidate in table.columns:
                    val_col = candidate
                    break

            if id_col is None or val_col is None:
                continue

            series = table.set_index(id_col)[val_col]
            series.name = gsm_name
            sample_dfs[gsm_name] = series

        if not sample_dfs:
            return None

        counts = pd.DataFrame(sample_dfs)
        counts.index.name = 'gene_id'

        # Convert to numeric
        counts = counts.apply(pd.to_numeric, errors='coerce')
        counts = counts.dropna(how='all')

        return counts

    def _counts_from_merge(self, gse) -> Optional[pd.DataFrame]:
        """Merge individual GSM tables as a fallback."""
        gsms = gse.gsms
        if not gsms:
            return None

        first_gsm = next(iter(gsms.values()))
        table = first_gsm.table
        if table is None or table.empty:
            return None

        # If table has more than 2 columns, it might already be a matrix
        if len(table.columns) > 3:
            # Assume first column is gene ID, rest are values
            gene_col = table.columns[0]
            df = table.set_index(gene_col)
            df = df.apply(pd.to_numeric, errors='coerce')
            df.index.name = 'gene_id'
            return df

        return None

    def _counts_from_supplementary(self, gse) -> Optional[pd.DataFrame]:
        """
        Download and parse supplementary count files from GEO.

        Most modern RNA-seq GEO studies store their count matrices as
        supplementary files (e.g., GSE12345_raw_counts.csv.gz). This
        method finds those files, downloads them, and parses them.
        """
        try:
            import requests
        except ImportError:
            logger.debug("requests not available for supplementary download")
            return None

        # Get supplementary file URLs from GSE metadata
        gse_meta = gse.metadata if hasattr(gse, 'metadata') else {}
        supp_files = gse_meta.get('supplementary_file', [])

        if isinstance(supp_files, str):
            supp_files = [supp_files]

        if not supp_files:
            logger.debug("No supplementary files listed in metadata")
            return None

        # Keywords that suggest a count matrix file
        count_keywords = [
            'count', 'counts', 'raw_count', 'gene_count',
            'readcount', 'read_count', 'expression', 'expr',
            'matrix', 'abundance', 'fpkm', 'tpm', 'rpkm',
        ]

        # Extensions we can parse
        parseable_ext = [
            '.csv.gz', '.tsv.gz', '.txt.gz',
            '.csv', '.tsv', '.txt',
            '.xlsx', '.xls',
        ]

        # Filter and rank supplementary files
        candidates = []
        for url in supp_files:
            url = url.strip()
            if not url:
                continue

            filename = url.split('/')[-1].lower()

            # Check extension
            is_parseable = any(filename.endswith(ext) for ext in parseable_ext)
            if not is_parseable:
                continue

            # Score by keyword match (higher = more likely count matrix)
            score = sum(1 for kw in count_keywords if kw in filename)

            # Boost score for count-specific names
            if 'count' in filename:
                score += 5
            if 'raw' in filename:
                score += 3

            # Skip obviously non-count files
            skip_keywords = ['readme', 'metadata', 'sample_info', 'design',
                             'annotation', 'filelist', '.bed', '.bw', '.bam']
            if any(kw in filename for kw in skip_keywords):
                continue

            candidates.append((score, url, filename))

        # Sort by score (best candidates first)
        candidates.sort(key=lambda x: -x[0])

        if not candidates:
            logger.debug("No parseable supplementary count files found")
            return None

        logger.info(f"Found {len(candidates)} candidate supplementary file(s)")

        # Try each candidate until one works
        for score, url, filename in candidates:
            logger.info(f"Trying supplementary file: {filename}")

            # Convert FTP to HTTPS if needed
            if url.startswith('ftp://'):
                url = url.replace('ftp://ftp.ncbi.nlm.nih.gov',
                                  'https://ftp.ncbi.nlm.nih.gov', 1)

            try:
                self._wait_for_rate_limit()
                resp = requests.get(url, timeout=120, stream=True)
                resp.raise_for_status()
                content = resp.content

                # Decompress if gzipped
                if filename.endswith('.gz'):
                    try:
                        content = gzip.decompress(content)
                    except (gzip.BadGzipFile, OSError):
                        pass  # wasn't actually gzipped

                # Detect separator
                if '.tsv' in filename or '.txt' in filename:
                    sep = '\t'
                else:
                    # Peek at content to auto-detect
                    first_line = content[:2000].decode('utf-8', errors='ignore').split('\n')[0]
                    sep = '\t' if '\t' in first_line else ','

                # Parse
                df = pd.read_csv(io.BytesIO(content), sep=sep, index_col=0)

                # Validate: should look like a gene x sample matrix
                # At least 100 genes and 2 samples, mostly numeric
                if df.shape[0] < 100 or df.shape[1] < 2:
                    logger.debug(
                        f"Skipping {filename}: too small "
                        f"({df.shape[0]} rows x {df.shape[1]} cols)"
                    )
                    continue

                # Try converting to numeric
                numeric_df = df.apply(pd.to_numeric, errors='coerce')
                numeric_pct = numeric_df.notna().sum().sum() / (df.shape[0] * df.shape[1])

                if numeric_pct < 0.5:
                    logger.debug(
                        f"Skipping {filename}: only {numeric_pct:.0%} numeric"
                    )
                    continue

                # Clean up
                # Drop rows that are entirely NaN
                numeric_df = numeric_df.dropna(how='all')
                # Drop COLUMNS that are entirely NaN (annotation cols like gene_name, description)
                numeric_df = numeric_df.dropna(axis=1, how='all')
                # Also drop columns where >90% of values are NaN (mostly text columns)
                col_nan_pct = numeric_df.isna().sum() / len(numeric_df)
                text_cols = col_nan_pct[col_nan_pct > 0.9].index
                if len(text_cols) > 0:
                    logger.info(f"Dropping non-numeric columns: {list(text_cols)}")
                    numeric_df = numeric_df.drop(columns=text_cols)
                numeric_df = numeric_df.fillna(0)
                numeric_df.index.name = 'gene_id'

                logger.info(
                    f"Successfully parsed {filename}: "
                    f"{numeric_df.shape[0]:,} genes x {numeric_df.shape[1]} samples"
                )
                return numeric_df

            except Exception as e:
                logger.debug(f"Failed to parse {filename}: {e}")
                continue

        return None

    def _align_metadata(
        self,
        counts_df: pd.DataFrame,
        metadata: pd.DataFrame,
        gse,
    ) -> pd.DataFrame:
        """
        Align metadata index with count matrix columns.

        When counts come from supplementary files, column names often
        don't match GSM IDs. This method tries multiple strategies:
        1. Direct match (columns are already GSM IDs)
        2. Map via GSM title/source_name fields
        3. Build fresh metadata from count columns
        """
        count_cols = set(counts_df.columns)
        meta_idx = set(metadata.index)

        # Strategy 1: Already aligned
        overlap = count_cols & meta_idx
        if len(overlap) > len(count_cols) * 0.5:
            logger.info(f"Metadata aligned: {len(overlap)}/{len(count_cols)} samples match")
            return metadata

        # Strategy 2: Map count columns to GSM IDs via title/source_name
        logger.info("Count columns don't match GSM IDs. Attempting mapping...")

        # Build lookup: various GSM fields -> GSM ID
        gsm_lookup = {}
        gsm_metadata_dict = {}
        for gsm_name, gsm in gse.gsms.items():
            gsm_meta = gsm.metadata if hasattr(gsm, 'metadata') else {}

            possible_names = [gsm_name]
            for field in ['title', 'source_name_ch1', 'description']:
                values = gsm_meta.get(field, [])
                if isinstance(values, list):
                    for v in values:
                        v = str(v).strip()
                        if v and v.lower() not in ('', 'na', 'none'):
                            possible_names.append(v)

            for name in possible_names:
                gsm_lookup[name] = gsm_name
                gsm_lookup[name.replace(' ', '_')] = gsm_name
                gsm_lookup[name.replace(' ', '-')] = gsm_name
                gsm_lookup[name.lower()] = gsm_name

            gsm_metadata_dict[gsm_name] = gsm_meta

        # Try to map each count column to a GSM ID
        col_to_gsm = {}
        for col in counts_df.columns:
            col_str = str(col)
            if col_str in gsm_lookup:
                col_to_gsm[col] = gsm_lookup[col_str]
            elif col_str.lower() in gsm_lookup:
                col_to_gsm[col] = gsm_lookup[col_str.lower()]
            else:
                for name, gsm_id in gsm_lookup.items():
                    if col_str in name or name in col_str:
                        col_to_gsm[col] = gsm_id
                        break

        mapped_pct = len(col_to_gsm) / max(len(counts_df.columns), 1)

        if mapped_pct > 0.5:
            logger.info(
                f"Mapped {len(col_to_gsm)}/{len(counts_df.columns)} "
                f"columns to GSM IDs"
            )
            counts_df.rename(columns=col_to_gsm, inplace=True)
            return metadata

        # Strategy 3: Build metadata from count columns directly
        logger.info(
            f"Could not map columns to GSM IDs ({len(col_to_gsm)} mapped). "
            f"Building metadata from count matrix columns."
        )

        gsm_list = list(gse.gsms.keys())

        meta_records = []
        for i, col in enumerate(counts_df.columns):
            record = {'sample_id': col, 'original_column': col}

            # Pull GSM metadata by position if available
            if i < len(gsm_list):
                gsm_name = gsm_list[i]
                gsm_meta = gsm_metadata_dict.get(gsm_name, {})
                record['gsm_id'] = gsm_name

                for field in ['title', 'source_name_ch1', 'organism_ch1']:
                    value = gsm_meta.get(field, [''])
                    if isinstance(value, list):
                        value = '; '.join(str(v) for v in value)
                    record[field] = value

                chars = gsm_meta.get('characteristics_ch1', [])
                if isinstance(chars, list):
                    for char in chars:
                        char = str(char)
                        if ':' in char:
                            k, v = char.split(':', 1)
                            record[k.strip().lower().replace(' ', '_')] = v.strip()

            meta_records.append(record)

        new_metadata = pd.DataFrame(meta_records)
        if 'sample_id' in new_metadata.columns:
            new_metadata = new_metadata.set_index('sample_id')
        new_metadata.index.name = 'sample_id'

        return new_metadata

    def _extract_metadata(self, gse) -> pd.DataFrame:
        """Extract sample metadata from GSE, including processing details."""
        meta_records = []

        for gsm_name, gsm in gse.gsms.items():
            record = {'sample_id': gsm_name}

            # Extract from GSM metadata
            gsm_meta = gsm.metadata if hasattr(gsm, 'metadata') else {}

            # Basic fields
            for key in ['title', 'source_name_ch1', 'organism_ch1',
                        'characteristics_ch1', 'description',
                        'platform_id', 'type']:
                value = gsm_meta.get(key, [''])
                if isinstance(value, list):
                    value = '; '.join(str(v) for v in value)
                record[key] = value

            # Processing / library fields (critical for pooling decisions)
            processing_keys = {
                'molecule_ch1': 'molecule',
                'extract_protocol_ch1': 'extraction_protocol',
                'library_strategy': 'library_strategy',
                'library_selection': 'library_selection',
                'library_source': 'library_source',
                'instrument_model': 'instrument',
                'data_processing': 'data_processing',
            }
            for geo_key, col_name in processing_keys.items():
                value = gsm_meta.get(geo_key, [''])
                if isinstance(value, list):
                    value = '; '.join(str(v) for v in value)
                record[col_name] = str(value).strip()

            # Parse assembly and quantification from data_processing
            dp = record.get('data_processing', '').lower()
            
            # Detect genome assembly
            assembly = 'unknown'
            for asm in ['grch38', 'grch37', 'hg38', 'hg19', 'hg18',
                        'mm39', 'mm10', 'mm9', 'grcm39', 'grcm38',
                        'rn7', 'rn6', 'danrer11', 'dm6']:
                if asm in dp:
                    assembly = asm
                    break
            record['assembly'] = assembly

            # Detect aligner
            aligner = 'unknown'
            for tool in ['star', 'hisat2', 'hisat', 'bowtie2', 'bowtie',
                         'bwa', 'tophat', 'salmon', 'kallisto', 'minimap2']:
                if tool in dp:
                    aligner = tool.upper() if tool in ['star', 'bwa'] else tool.capitalize()
                    break
            record['aligner'] = aligner

            # Detect quantification tool
            quant_tool = 'unknown'
            for tool in ['featurecounts', 'htseq', 'stringtie', 'ballgown',
                         'salmon', 'kallisto', 'rsem', 'cufflinks']:
                if tool in dp:
                    quant_tool = tool.capitalize()
                    break
            record['quantification'] = quant_tool

            # Clean molecule type
            mol = record.get('molecule', '').lower()
            if 'polya' in mol or 'poly-a' in mol or 'poly(a)' in mol:
                record['rna_type'] = 'polyA'
            elif 'total' in mol:
                record['rna_type'] = 'total RNA'
            else:
                record['rna_type'] = record.get('molecule', 'unknown')

            # Parse characteristics into separate columns
            chars = gsm_meta.get('characteristics_ch1', [])
            if isinstance(chars, list):
                for char in chars:
                    char = str(char)
                    if ':' in char:
                        k, v = char.split(':', 1)
                        record[k.strip().lower().replace(' ', '_')] = v.strip()

            meta_records.append(record)

        metadata = pd.DataFrame(meta_records)
        if 'sample_id' in metadata.columns:
            metadata = metadata.set_index('sample_id')

        metadata.index.name = 'sample_id'
        return metadata

    def _align_metadata(
        self,
        counts_df: pd.DataFrame,
        metadata: pd.DataFrame,
        gse,
    ) -> pd.DataFrame:
        """
        Align metadata index with count matrix columns.

        Supplementary files often use sample names (PS19_1, WT_1)
        while GSM metadata uses accession IDs (GSM1234567). This method
        tries multiple strategies to align them, always preserving
        the GSM accession as a column.
        """
        import re

        count_cols = list(counts_df.columns)
        count_cols_set = set(count_cols)
        meta_idx = set(metadata.index)

        # Always preserve GSM IDs as a column
        gsm_ids = list(metadata.index)

        # Build lookup dicts from GSM metadata
        gsm_titles = {}
        gsm_sources = {}
        for gsm_name, gsm in gse.gsms.items():
            gsm_meta = gsm.metadata if hasattr(gsm, 'metadata') else {}
            title = gsm_meta.get('title', [''])
            if isinstance(title, list):
                title = title[0] if title else ''
            gsm_titles[gsm_name] = str(title).strip()

            source = gsm_meta.get('source_name_ch1', [''])
            if isinstance(source, list):
                source = source[0] if source else ''
            gsm_sources[gsm_name] = str(source).strip()

        # --- Strategy 1: Already aligned (count cols = GSM IDs) ---
        overlap = count_cols_set & meta_idx
        if len(overlap) > len(count_cols) * 0.5:
            logger.debug("Metadata already aligned with count matrix")
            metadata['gsm_id'] = metadata.index
            return metadata

        # --- Strategy 2: Exact title match (col name = GSM title) ---
        title_to_gsm = {t: g for g, t in gsm_titles.items() if t}
        matched = {}
        for col in count_cols:
            if col in title_to_gsm:
                matched[col] = title_to_gsm[col]
            else:
                for title, gsm in title_to_gsm.items():
                    if col.lower() == title.lower():
                        matched[col] = gsm
                        break

        if len(matched) > len(count_cols) * 0.5:
            logger.info(f"Aligned {len(matched)}/{len(count_cols)} via title matching")
            return self._build_aligned_meta(metadata, matched, count_cols)

        # --- Strategy 3: Source name match ---
        source_to_gsm = {s: g for g, s in gsm_sources.items() if s}
        matched = {}
        for col in count_cols:
            if col in source_to_gsm:
                matched[col] = source_to_gsm[col]
            else:
                for source, gsm in source_to_gsm.items():
                    if col.lower() == source.lower():
                        matched[col] = gsm
                        break

        if len(matched) > len(count_cols) * 0.5:
            logger.info(f"Aligned {len(matched)}/{len(count_cols)} via source matching")
            return self._build_aligned_meta(metadata, matched, count_cols)

        # --- Strategy 4: Fuzzy substring matching ---
        # Column "PS19_1" should match title "PS19,10month,rep1"
        # Uses base name matching + replicate number matching
        # Only assigns matches when confident (base name found in title)

        def _extract_rep_num(name):
            """Extract trailing number: 'PS19_1' → 1, 'WT_3' → 3"""
            m = re.search(r'[\._\-]?(\d+)$', name)
            return int(m.group(1)) if m else None

        def _extract_rep_from_title(title):
            """Extract rep number: 'WT,10month,rep1' → 1"""
            m = re.search(r'rep(?:licate)?[\._\-\s]?(\d+)', title, re.IGNORECASE)
            if m:
                return int(m.group(1))
            m = re.search(r'[\._\-](\d+)\s*$', title)
            return int(m.group(1)) if m else None

        def _base_name(name):
            """Strip trailing number: 'PS19_1' → 'ps19', 'WT_3' → 'wt'"""
            return re.sub(r'[\._\-]?\d+$', '', name).lower().strip('_- ')

        # Compute all pairwise scores
        score_matrix = {}
        for col in count_cols:
            col_base = _base_name(col)
            col_rep = _extract_rep_num(col)
            score_matrix[col] = {}

            for gsm_name, title in gsm_titles.items():
                title_lower = title.lower()
                title_rep = _extract_rep_from_title(title)
                title_base = _base_name(title)

                score = 0

                # Core: base name must appear in title (or vice versa)
                base_match = False
                if col_base and (col_base in title_lower or title_base in col.lower()):
                    score += 20
                    base_match = True

                # Bonus: prefer titles with fewer extra tokens (closer match)
                # "PS19_1" should prefer "PS19,10month,rep1" over "PS19, Zbp1+-,10month,rep1"
                if base_match:
                    col_len = len(col_base)
                    title_stripped = re.sub(r'[\d,\s\-\.\+]+|rep\d*|10month', '', title_lower).strip()
                    len_diff = abs(len(title_stripped) - col_len)
                    # Bonus inversely proportional to length difference
                    score += max(0, 15 - len_diff)

                # Only add rep bonus if base matched
                if base_match and col_rep is not None and title_rep is not None:
                    if col_rep == title_rep:
                        score += 15
                    else:
                        score -= 10

                score_matrix[col][gsm_name] = score

        # Greedy assignment: match highest scores first
        matched = {}
        used_gsms = set()

        # Sort all (col, gsm) pairs by score descending
        all_pairs = []
        for col in count_cols:
            for gsm_name, score in score_matrix[col].items():
                if score > 0:
                    all_pairs.append((score, col, gsm_name))
        all_pairs.sort(reverse=True)

        for score, col, gsm_name in all_pairs:
            if col not in matched and gsm_name not in used_gsms:
                matched[col] = gsm_name
                used_gsms.add(gsm_name)

        if len(matched) > len(count_cols) * 0.5:
            logger.info(
                f"Aligned {len(matched)}/{len(count_cols)} via fuzzy token matching"
            )
            return self._build_aligned_meta(metadata, matched, count_cols)

        # --- Strategy 5: Positional (last resort, warn user) ---
        if len(count_cols) == len(metadata):
            logger.warning(
                f"Using positional alignment for {len(count_cols)} samples. "
                f"Verify that sample order matches!"
            )
            matched = {col: gsm for col, gsm in zip(count_cols, metadata.index)}
            return self._build_aligned_meta(metadata, matched, count_cols)

        # --- Strategy 6: Build basic metadata ---
        logger.warning(
            f"Could not align metadata. Creating basic metadata. "
            f"Count cols: {count_cols[:3]}..., GSMs: {gsm_ids[:3]}..."
        )
        basic_meta = pd.DataFrame(index=count_cols)
        basic_meta.index.name = 'sample_id'
        basic_meta['gsm_id'] = ''
        # Include GSM info as reference
        for i, gsm_name in enumerate(gsm_ids):
            if i < len(count_cols):
                basic_meta.loc[count_cols[i], 'gsm_id'] = gsm_name
                basic_meta.loc[count_cols[i], 'gsm_title'] = gsm_titles.get(gsm_name, '')

        return basic_meta

    def _build_aligned_meta(
        self,
        metadata: pd.DataFrame,
        col_to_gsm: dict,
        count_cols: list,
    ) -> pd.DataFrame:
        """
        Build aligned metadata from a column→GSM mapping.
        Always includes gsm_id as a column.
        """
        # For each count column, pull the matched GSM's metadata row
        rows = []
        for col in count_cols:
            gsm = col_to_gsm.get(col)
            if gsm and gsm in metadata.index:
                row = metadata.loc[gsm].to_dict()
                row['gsm_id'] = gsm
            else:
                row = {'gsm_id': gsm or ''}
            rows.append(row)

        new_meta = pd.DataFrame(rows, index=count_cols)
        new_meta.index.name = 'sample_id'

        return new_meta

    def _extract_organism(self, gse) -> str:
        """Extract organism from GSE metadata."""
        gse_meta = gse.metadata if hasattr(gse, 'metadata') else {}

        # Strategy 1: GSE-level platform_organism
        for key in ['platform_organism', 'organism', 'sample_organism']:
            organism = gse_meta.get(key, [])
            if isinstance(organism, list):
                organism = organism[0] if organism else ''
            organism = str(organism).strip()
            if organism and organism not in ('', 'unknown'):
                return organism

        # Strategy 2: First GSM organism_ch1 (most reliable)
        if gse.gsms:
            for gsm in gse.gsms.values():
                gsm_meta = gsm.metadata if hasattr(gsm, 'metadata') else {}
                org = gsm_meta.get('organism_ch1', [])
                if isinstance(org, list):
                    org = org[0] if org else ''
                org = str(org).strip()
                if org and org not in ('', 'unknown'):
                    return org

        # Strategy 3: GPL metadata
        if gse.gpls:
            for gpl in gse.gpls.values():
                gpl_meta = gpl.metadata if hasattr(gpl, 'metadata') else {}
                org = gpl_meta.get('organism', [])
                if isinstance(org, list):
                    org = org[0] if org else ''
                org = str(org).strip()
                if org and org not in ('', 'unknown'):
                    return org

        return 'unknown'

    def _extract_platform(self, gse) -> str:
        """Extract platform from GSE metadata."""
        gse_meta = gse.metadata if hasattr(gse, 'metadata') else {}
        platform = gse_meta.get('platform_id', ['unknown'])
        if isinstance(platform, list) and platform:
            return platform[0]
        return str(platform)

    @staticmethod
    def _count_supplementary(gse) -> int:
        """Count supplementary files."""
        gse_meta = gse.metadata if hasattr(gse, 'metadata') else {}
        supp = gse_meta.get('supplementary_file', [])
        if isinstance(supp, list):
            return len(supp)
        return 0

    @staticmethod
    def _detect_gene_id_type(counts_df: pd.DataFrame) -> str:
        """Heuristically detect gene ID type from the index."""
        sample = counts_df.index[:100].astype(str).tolist()

        ensembl = sum(1 for g in sample if str(g).startswith('ENS'))
        if ensembl > len(sample) * 0.5:
            return 'ensembl'

        numeric = sum(1 for g in sample if str(g).isdigit())
        if numeric > len(sample) * 0.5:
            return 'entrez'

        return 'symbol'

    # -------------------------------------------------------------------------
    # Study Info (preview before download)
    # -------------------------------------------------------------------------

    def get_study_info(self, accession: str) -> Dict[str, Any]:
        """
        Get comprehensive study information without downloading count data.

        Fetches the GSE SOFT file via GEOparse for full metadata including
        sample details, platform info, and supplementary file list.
        Much faster than a full download since it skips count extraction.

        Parameters
        ----------
        accession : str
            GEO accession (GSE...).

        Returns
        -------
        dict
            Study info with keys: accession, title, summary, organism,
            platform_id, platform_name, n_samples, submission_date,
            samples (DataFrame), supplementary_files, overall_design,
            pubmed_ids, contact, data_type.
        """
        GEOparse = _check_geoparse()

        if not accession.startswith('GSE'):
            raise ValidationError('accession',
                f"Expected GSE accession, got: {accession}",
                "Use GSE followed by numbers (e.g., GSE306761)"
            )

        logger.info(f"Fetching study info for {accession}...")

        with tempfile.TemporaryDirectory() as tmpdir:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    gse = GEOparse.get_GEO(
                        geo=accession,
                        destdir=tmpdir,
                        silent=True,
                    )
                except Exception as e:
                    raise RAPTORError(f"Failed to fetch {accession}: {e}")

        # Extract GSE-level metadata
        gse_meta = gse.metadata if hasattr(gse, 'metadata') else {}

        def _get(key, default=''):
            val = gse_meta.get(key, [default])
            if isinstance(val, list):
                return '; '.join(str(v) for v in val) if val else default
            return str(val)

        # Build sample table
        sample_records = []
        for gsm_name, gsm in gse.gsms.items():
            gsm_meta = gsm.metadata if hasattr(gsm, 'metadata') else {}

            record = {
                'GSM': gsm_name,
                'Title': _list_to_str(gsm_meta.get('title', [''])),
                'Source': _list_to_str(gsm_meta.get('source_name_ch1', [''])),
                'Organism': _list_to_str(gsm_meta.get('organism_ch1', [''])),
                'Platform': _list_to_str(gsm_meta.get('platform_id', [''])),
                'Molecule': _list_to_str(gsm_meta.get('molecule_ch1', [''])),
                'Instrument': _list_to_str(gsm_meta.get('instrument_model', [''])),
                'Library Selection': _list_to_str(gsm_meta.get('library_selection', [''])),
            }

            # Classify RNA type
            mol = record['Molecule'].lower()
            if 'polya' in mol or 'poly-a' in mol or 'poly(a)' in mol:
                record['RNA Type'] = 'polyA'
            elif 'total' in mol:
                record['RNA Type'] = 'total RNA'
            else:
                record['RNA Type'] = record['Molecule'] or 'unknown'

            # Parse data processing for assembly and tools
            dp = _list_to_str(gsm_meta.get('data_processing', [''])).lower()

            assembly = ''
            for asm in ['grch38', 'grch37', 'hg38', 'hg19',
                        'mm39', 'mm10', 'mm9', 'grcm39', 'grcm38',
                        'rn7', 'rn6']:
                if asm in dp:
                    assembly = asm
                    break
            record['Assembly'] = assembly

            aligner = ''
            for tool in ['star', 'hisat2', 'hisat', 'bowtie2',
                         'salmon', 'kallisto', 'tophat', 'bwa']:
                if tool in dp:
                    aligner = tool.upper() if tool in ['star', 'bwa'] else tool.capitalize()
                    break
            record['Aligner'] = aligner

            quant = ''
            for tool in ['featurecounts', 'htseq', 'stringtie',
                         'salmon', 'kallisto', 'rsem', 'cufflinks']:
                if tool in dp:
                    quant = tool.capitalize()
                    break
            record['Quantification'] = quant

            # Parse characteristics
            chars = gsm_meta.get('characteristics_ch1', [])
            if isinstance(chars, list):
                for char in chars:
                    char = str(char)
                    if ':' in char:
                        k, v = char.split(':', 1)
                        key = k.strip().replace(' ', '_').lower()
                        record[key] = v.strip()

            sample_records.append(record)

        samples_df = pd.DataFrame(sample_records)

        # Supplementary files
        supp_files = gse_meta.get('supplementary_file', [])
        if isinstance(supp_files, str):
            supp_files = [supp_files]

        # Platform info
        platform_id = _get('platform_id', 'unknown')
        platform_name = ''
        if gse.gpls:
            first_gpl = next(iter(gse.gpls.values()))
            gpl_meta = first_gpl.metadata if hasattr(first_gpl, 'metadata') else {}
            platform_name = _list_to_str(gpl_meta.get('title', ['']))

        # Detect data type from supplementary files
        data_type = 'unknown'
        supp_lower = ' '.join(str(f).lower() for f in supp_files)
        if any(kw in supp_lower for kw in ['count', 'raw_count', 'readcount']):
            data_type = 'raw_counts'
        elif any(kw in supp_lower for kw in ['fpkm', 'tpm', 'rpkm', 'normalized']):
            data_type = 'normalized'
        elif any(kw in supp_lower for kw in ['.fastq', '.sra', '.bam']):
            data_type = 'raw_reads'

        # Detect organism from multiple sources
        organism = _get('platform_organism', '')
        if not organism or organism == 'unknown':
            # Try from first sample
            if sample_records:
                organism = sample_records[0].get('Organism', '')
        if not organism or organism == 'unknown':
            # Try from GPL
            if gse.gpls:
                first_gpl = next(iter(gse.gpls.values()))
                gpl_meta = first_gpl.metadata if hasattr(first_gpl, 'metadata') else {}
                org_list = gpl_meta.get('organism', [''])
                organism = _list_to_str(org_list)
        if not organism:
            organism = 'unknown'

        # Fix encoding issues (smart quotes etc.)
        def _fix_encoding(text):
            if not isinstance(text, str):
                return text
            try:
                # Fix mojibake from UTF-8 interpreted as latin-1
                return text.encode('latin-1').decode('utf-8')
            except (UnicodeDecodeError, UnicodeEncodeError):
                return text

        title = _fix_encoding(_get('title'))
        summary = _fix_encoding(_get('summary'))
        overall_design = _fix_encoding(_get('overall_design'))

        return {
            'accession': accession,
            'title': title,
            'summary': summary,
            'organism': organism,
            'overall_design': overall_design,
            'platform_id': platform_id,
            'platform_name': platform_name,
            'n_samples': len(gse.gsms),
            'submission_date': _get('submission_date'),
            'last_update': _get('last_update_date'),
            'pubmed_ids': gse_meta.get('pubmed_id', []),
            'contact': _get('contact_name'),
            'data_type': data_type,
            'supplementary_files': supp_files,
            'samples': samples_df,
            'repository': 'GEO',
        }

    # -------------------------------------------------------------------------
    # Info (lightweight, Entrez only)
    # -------------------------------------------------------------------------

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """Get metadata for a GSE without downloading the full dataset."""
        Entrez = _check_entrez()
        Entrez.email = self._email
        if self._api_key:
            Entrez.api_key = self._api_key

        try:
            # Search for the specific accession
            handle = Entrez.esearch(
                db='gds',
                term=f'{accession}[Accession]',
                retmax=1,
            )
            result = Entrez.read(handle)
            handle.close()

            id_list = result.get('IdList', [])
            if not id_list:
                return None

            self._wait_for_rate_limit()
            handle = Entrez.esummary(db='gds', id=id_list[0])
            summary = Entrez.read(handle)
            handle.close()

            if summary:
                s = summary[0]
                return {
                    'accession': s.get('Accession', accession),
                    'title': s.get('title', ''),
                    'organism': s.get('taxon', 'unknown'),
                    'n_samples': int(s.get('n_samples', 0)),
                    'platform': s.get('GPL', 'unknown'),
                    'description': s.get('summary', ''),
                    'gds_type': s.get('gdsType', ''),
                    'repository': 'GEO',
                }

        except Exception as e:
            logger.warning(f"Info query failed for {accession}: {e}")

        return None

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"GEOConnector("
            f"email='{self._email[:20]}...', "
            f"api_key={'set' if self._api_key else 'none'}, "
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )
