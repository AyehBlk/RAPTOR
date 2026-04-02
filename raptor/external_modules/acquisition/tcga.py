"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - TCGA/GDC Connector

Connects to the Genomic Data Commons (GDC) REST API for downloading
TCGA and TARGET RNA-seq datasets. No authentication required for
open-access data.

GDC API docs: https://docs.gdc.cancer.gov/API/Users_Guide/

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import io
import gzip
import json
import tempfile
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

GDC_API_BASE = "https://api.gdc.cancer.gov"
GDC_FILES_ENDPOINT = f"{GDC_API_BASE}/files"
GDC_CASES_ENDPOINT = f"{GDC_API_BASE}/cases"
GDC_PROJECTS_ENDPOINT = f"{GDC_API_BASE}/projects"
GDC_DATA_ENDPOINT = f"{GDC_API_BASE}/data"

# TCGA project IDs (common ones)
TCGA_PROJECTS = {
    'TCGA-BRCA': 'Breast invasive carcinoma',
    'TCGA-LUAD': 'Lung adenocarcinoma',
    'TCGA-LUSC': 'Lung squamous cell carcinoma',
    'TCGA-COAD': 'Colon adenocarcinoma',
    'TCGA-LIHC': 'Liver hepatocellular carcinoma',
    'TCGA-PRAD': 'Prostate adenocarcinoma',
    'TCGA-KIRC': 'Kidney renal clear cell carcinoma',
    'TCGA-HNSC': 'Head and neck squamous cell carcinoma',
    'TCGA-THCA': 'Thyroid carcinoma',
    'TCGA-STAD': 'Stomach adenocarcinoma',
    'TCGA-BLCA': 'Bladder urothelial carcinoma',
    'TCGA-OV': 'Ovarian serous cystadenocarcinoma',
    'TCGA-GBM': 'Glioblastoma multiforme',
    'TCGA-PAAD': 'Pancreatic adenocarcinoma',
    'TCGA-SKCM': 'Skin cutaneous melanoma',
    'TCGA-UCEC': 'Uterine corpus endometrial carcinoma',
}


def _check_requests():
    """Check if requests is available."""
    try:
        import requests
        return requests
    except ImportError:
        raise RAPTORError(
            "requests is required for GDC API access. "
            "Install with: pip install requests"
        )


# =============================================================================
# TCGAConnector
# =============================================================================

class TCGAConnector(BaseConnector):
    """
    Connector for TCGA/TARGET data via the GDC REST API.

    Downloads open-access RNA-seq gene expression quantification files
    (HTSeq counts or STAR counts) from the Genomic Data Commons.

    Parameters
    ----------
    cache : bool
        Enable disk caching.
    cache_dir : str or Path, optional
        Custom cache directory.
    workflow_type : str
        Which quantification workflow to download.
        Options: 'STAR - Counts', 'HTSeq - Counts'.

    Examples
    --------
    >>> tcga = TCGAConnector()
    >>> results = tcga.search("liver hepatocellular")
    >>> dataset = tcga.download("TCGA-LIHC")
    >>> print(dataset.summary())

    >>> # List available TCGA projects
    >>> tcga.list_projects()
    """

    REPOSITORY_NAME = 'TCGA'

    def __init__(
        self,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
        workflow_type: str = 'STAR - Counts',
    ):
        super().__init__(
            cache=cache,
            cache_dir=cache_dir,
            rate_limit=0.2,  # GDC is generous but be polite
        )
        self._workflow_type = workflow_type

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
        """Search GDC projects matching the query."""
        requests = _check_requests()

        # Search projects endpoint
        params = {
            'filters': json.dumps({
                'op': 'and',
                'content': [
                    {
                        'op': 'in',
                        'content': {
                            'field': 'program.name',
                            'value': ['TCGA', 'TARGET'],
                        }
                    }
                ]
            }),
            'fields': (
                'project_id,name,primary_site,disease_type,'
                'summary.case_count,summary.file_count,'
                'summary.data_categories.data_category,'
                'summary.data_categories.file_count'
            ),
            'size': 100,
        }

        try:
            resp = requests.get(GDC_PROJECTS_ENDPOINT, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(f"GDC project search failed: {e}")

        query_lower = query.lower()
        results = []

        for hit in data.get('data', {}).get('hits', []):
            project_id = hit.get('project_id', '')
            name = hit.get('name', '')
            primary_site = hit.get('primary_site', ['unknown'])
            disease_type = hit.get('disease_type', ['unknown'])

            if isinstance(primary_site, list):
                primary_site = '; '.join(primary_site)
            if isinstance(disease_type, list):
                disease_type = '; '.join(disease_type)

            searchable = f"{project_id} {name} {primary_site} {disease_type}".lower()

            if query_lower in searchable or any(
                term in searchable for term in query_lower.split()
            ):
                case_count = hit.get('summary', {}).get('case_count', 0)
                results.append(SearchResult(
                    accession=project_id,
                    title=name,
                    organism='Homo sapiens',
                    n_samples=case_count,
                    platform='Illumina HiSeq/NovaSeq',
                    description=f"{disease_type} | Site: {primary_site}",
                    repository='TCGA',
                    primary_site=primary_site,
                    disease_type=disease_type,
                ))

        # Also check the static dictionary for known projects
        for proj_id, desc in TCGA_PROJECTS.items():
            if query_lower in f"{proj_id} {desc}".lower():
                if not any(r.accession == proj_id for r in results):
                    results.append(SearchResult(
                        accession=proj_id,
                        title=desc,
                        organism='Homo sapiens',
                        n_samples=0,
                        platform='Illumina HiSeq',
                        description=desc,
                        repository='TCGA',
                    ))

        return results[:max_results]

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    def _download_api(
        self,
        accession: str,
        data_type: str = 'raw_counts',
    ) -> AcquiredDataset:
        """Download RNA-seq counts for a TCGA project."""
        requests = _check_requests()

        logger.info(f"Querying GDC for {accession} RNA-seq files...")

        # Step 1: Find all RNA-seq count files for this project
        file_ids, file_metadata = self._query_files(accession, requests)

        if not file_ids:
            raise RAPTORError(
                f"No RNA-seq count files found for {accession}. "
                f"Check the project ID or try a different workflow_type."
            )

        logger.info(f"Found {len(file_ids)} count files for {accession}")

        # Step 2: Download files (GDC supports batch download)
        counts_df = self._download_counts(file_ids, file_metadata, requests)

        # Step 3: Build metadata from case info
        metadata = self._build_metadata(file_metadata)

        source_info = {
            'repository': 'TCGA',
            'accession': accession,
            'organism': 'Homo sapiens',
            'platform': 'Illumina HiSeq/NovaSeq',
            'download_date': datetime.now().isoformat(),
            'data_type': data_type,
            'workflow_type': self._workflow_type,
            'description': TCGA_PROJECTS.get(accession, accession),
            'n_files_downloaded': len(file_ids),
        }

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info=source_info,
            gene_id_type='ensembl',
        )

    def _query_files(
        self,
        project_id: str,
        requests,
    ) -> tuple:
        """Query GDC for RNA-seq quantification files."""
        filters = {
            'op': 'and',
            'content': [
                {
                    'op': '=',
                    'content': {
                        'field': 'cases.project.project_id',
                        'value': project_id,
                    }
                },
                {
                    'op': '=',
                    'content': {
                        'field': 'data_category',
                        'value': 'Transcriptome Profiling',
                    }
                },
                {
                    'op': '=',
                    'content': {
                        'field': 'data_type',
                        'value': 'Gene Expression Quantification',
                    }
                },
                {
                    'op': '=',
                    'content': {
                        'field': 'analysis.workflow_type',
                        'value': self._workflow_type,
                    }
                },
                {
                    'op': '=',
                    'content': {
                        'field': 'access',
                        'value': 'open',
                    }
                },
            ]
        }

        params = {
            'filters': json.dumps(filters),
            'fields': (
                'file_id,file_name,cases.case_id,'
                'cases.submitter_id,cases.samples.sample_type,'
                'cases.samples.tissue_type,'
                'cases.demographic.gender,'
                'cases.demographic.vital_status,'
                'cases.diagnoses.primary_diagnosis,'
                'cases.diagnoses.tumor_stage'
            ),
            'size': 2000,
        }

        self._wait_for_rate_limit()
        try:
            resp = requests.get(GDC_FILES_ENDPOINT, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(f"GDC file query failed: {e}")

        hits = data.get('data', {}).get('hits', [])
        file_ids = []
        file_metadata = {}

        for hit in hits:
            fid = hit.get('file_id')
            if fid:
                file_ids.append(fid)

                # Extract case metadata
                cases = hit.get('cases', [{}])
                case = cases[0] if cases else {}
                submitter_id = case.get('submitter_id', fid[:12])

                samples = case.get('samples', [{}])
                sample = samples[0] if samples else {}

                demographics = case.get('demographic', {})
                diagnoses = case.get('diagnoses', [{}])
                diagnosis = diagnoses[0] if diagnoses else {}

                file_metadata[fid] = {
                    'case_id': case.get('case_id', ''),
                    'submitter_id': submitter_id,
                    'sample_type': sample.get('sample_type', 'unknown'),
                    'tissue_type': sample.get('tissue_type', 'unknown'),
                    'gender': demographics.get('gender', 'unknown') if isinstance(demographics, dict) else 'unknown',
                    'vital_status': demographics.get('vital_status', 'unknown') if isinstance(demographics, dict) else 'unknown',
                    'primary_diagnosis': diagnosis.get('primary_diagnosis', 'unknown') if isinstance(diagnosis, dict) else 'unknown',
                    'tumor_stage': diagnosis.get('tumor_stage', 'unknown') if isinstance(diagnosis, dict) else 'unknown',
                    'file_name': hit.get('file_name', ''),
                }

        return file_ids, file_metadata

    def _download_counts(
        self,
        file_ids: List[str],
        file_metadata: dict,
        requests,
    ) -> pd.DataFrame:
        """
        Download count files from GDC.

        Uses the GDC data endpoint for individual file downloads.
        Parses the TSV count files and merges into a matrix.
        """
        sample_counts = {}

        for i, fid in enumerate(file_ids):
            if i % 50 == 0 and i > 0:
                logger.info(f"  Downloaded {i}/{len(file_ids)} files...")

            self._wait_for_rate_limit()

            try:
                resp = requests.get(
                    f"{GDC_DATA_ENDPOINT}/{fid}",
                    timeout=60,
                )
                resp.raise_for_status()
            except Exception as e:
                logger.warning(f"Failed to download file {fid}: {e}")
                continue

            # Parse the count file
            content = resp.content

            # GDC may return gzipped content
            try:
                content = gzip.decompress(content)
            except (gzip.BadGzipFile, OSError):
                pass  # not gzipped

            try:
                df = pd.read_csv(
                    io.BytesIO(content),
                    sep='\t',
                    comment='#',
                    header=0,
                )
            except Exception as e:
                logger.warning(f"Failed to parse file {fid}: {e}")
                continue

            # Find gene ID and count columns
            gene_col = None
            count_col = None

            for col in df.columns:
                col_lower = col.lower()
                if 'gene_id' in col_lower or 'ensembl' in col_lower:
                    gene_col = col
                elif 'unstranded' in col_lower or 'count' in col_lower:
                    count_col = col

            if gene_col is None:
                gene_col = df.columns[0]
            if count_col is None:
                # Try common column names
                for candidate in ['unstranded', 'tpm_unstranded',
                                  'fpkm_unstranded', 'expected_count']:
                    if candidate in df.columns:
                        count_col = candidate
                        break
                if count_col is None and len(df.columns) >= 2:
                    count_col = df.columns[1]

            if gene_col is None or count_col is None:
                continue

            # Strip version numbers from Ensembl IDs (ENSG00000141510.17 → ENSG00000141510)
            genes = df[gene_col].astype(str).str.split('.').str[0]

            # Filter out non-gene rows (STAR counts have summary rows)
            gene_mask = genes.str.startswith('ENSG')

            submitter_id = file_metadata.get(fid, {}).get('submitter_id', fid[:12])
            series = pd.Series(
                df.loc[gene_mask, count_col].values,
                index=genes[gene_mask].values,
                name=submitter_id,
                dtype=float,
            )

            sample_counts[submitter_id] = series

        if not sample_counts:
            raise RAPTORError(
                "No count data could be extracted from downloaded files"
            )

        counts_df = pd.DataFrame(sample_counts)
        counts_df.index.name = 'gene_id'
        counts_df = counts_df.fillna(0)

        logger.info(
            f"Built count matrix: "
            f"{counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples"
        )

        return counts_df

    def _build_metadata(self, file_metadata: dict) -> pd.DataFrame:
        """Build sample metadata DataFrame from file metadata."""
        records = []
        seen_submitters = set()

        for fid, meta in file_metadata.items():
            submitter_id = meta.get('submitter_id', fid[:12])
            if submitter_id in seen_submitters:
                continue
            seen_submitters.add(submitter_id)

            records.append({
                'sample_id': submitter_id,
                'sample_type': meta.get('sample_type', 'unknown'),
                'tissue_type': meta.get('tissue_type', 'unknown'),
                'gender': meta.get('gender', 'unknown'),
                'vital_status': meta.get('vital_status', 'unknown'),
                'primary_diagnosis': meta.get('primary_diagnosis', 'unknown'),
                'tumor_stage': meta.get('tumor_stage', 'unknown'),
            })

        metadata = pd.DataFrame(records)
        if 'sample_id' in metadata.columns:
            metadata = metadata.set_index('sample_id')
        metadata.index.name = 'sample_id'

        return metadata

    # -------------------------------------------------------------------------
    # Info
    # -------------------------------------------------------------------------

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """Get project info from GDC."""
        requests = _check_requests()

        try:
            resp = requests.get(
                f"{GDC_PROJECTS_ENDPOINT}/{accession}",
                params={
                    'fields': 'name,primary_site,disease_type,summary.case_count'
                },
                timeout=15,
            )
            if resp.status_code == 200:
                data = resp.json().get('data', {})
                return {
                    'accession': accession,
                    'title': data.get('name', ''),
                    'organism': 'Homo sapiens',
                    'n_samples': data.get('summary', {}).get('case_count', 0),
                    'primary_site': data.get('primary_site', []),
                    'disease_type': data.get('disease_type', []),
                    'repository': 'TCGA',
                }
        except Exception as e:
            logger.debug(f"GDC info failed for {accession}: {e}")

        return None

    def list_projects(self) -> Dict[str, str]:
        """
        List known TCGA project IDs and descriptions.

        Returns
        -------
        dict
            Project ID → description mapping.
        """
        return dict(TCGA_PROJECTS)

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"TCGAConnector("
            f"workflow='{self._workflow_type}', "
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )
