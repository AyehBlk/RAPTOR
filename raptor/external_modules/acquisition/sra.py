"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - SRA Connector

Connects to NCBI Sequence Read Archive for searching and retrieving
RNA-seq run metadata. SRA provides raw reads rather than processed
counts, so this connector focuses on:

1. Searching and discovering SRA studies/runs
2. Downloading run metadata (sample info, library strategy, etc.)
3. Generating shell commands for fastq-dump/fasterq-dump
4. Importing pre-quantified count matrices from SRA-linked studies

For actual FASTQ download and quantification, use RAPTOR's Module 1
(Quick Quantification) pipelines after getting run accessions from here.

Dependencies:
    pip install biopython requests

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import os
import re
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Union

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

# EBI ENA API (provides structured SRA metadata)
ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# NCBI Entrez for searching SRA
SRA_DB = 'sra'

# Useful SRA fields from ENA
ENA_RUN_FIELDS = (
    'run_accession,experiment_accession,study_accession,'
    'sample_accession,instrument_model,library_strategy,'
    'library_source,library_layout,read_count,'
    'base_count,fastq_ftp,fastq_bytes,'
    'sample_alias,sample_title,scientific_name,'
    'study_title,experiment_title'
)


def _check_requests():
    try:
        import requests
        return requests
    except ImportError:
        raise RAPTORError(
            "requests is required for SRA access. "
            "Install with: pip install requests"
        )


def _check_entrez():
    try:
        from Bio import Entrez
        return Entrez
    except ImportError:
        raise RAPTORError(
            "Biopython is required for SRA search. "
            "Install with: pip install biopython"
        )


# =============================================================================
# SRAConnector
# =============================================================================

class SRAConnector(BaseConnector):
    """
    Connector for NCBI Sequence Read Archive (SRA).

    Searches for RNA-seq studies in SRA, retrieves run metadata, and
    generates download commands. Since SRA stores raw reads (FASTQ),
    not processed counts, this connector provides:

    - Study/run discovery and metadata retrieval
    - Run table export (for use with RAPTOR's quantification pipelines)
    - Shell command generation for fasterq-dump
    - Import of pre-computed count matrices if available

    Parameters
    ----------
    email : str, optional
        Email for NCBI Entrez.
    cache : bool
        Enable disk caching.
    cache_dir : str or Path, optional
        Custom cache directory.

    Examples
    --------
    >>> sra = SRAConnector(email="user@example.com")
    >>> results = sra.search("liver RNA-seq Homo sapiens")
    >>> run_table = sra.get_run_table("SRP123456")
    >>> print(run_table.head())

    >>> # Generate download commands
    >>> commands = sra.generate_download_commands("SRP123456", "/data/fastq")
    >>> for cmd in commands[:3]:
    ...     print(cmd)
    """

    REPOSITORY_NAME = 'SRA'

    def __init__(
        self,
        email: Optional[str] = None,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
    ):
        self._email = email or os.environ.get('NCBI_EMAIL', '')

        super().__init__(
            cache=cache,
            cache_dir=cache_dir,
            rate_limit=0.34,
        )

        if not self._email:
            logger.warning(
                "No email set for NCBI Entrez. Set via email parameter "
                "or NCBI_EMAIL environment variable."
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
        """Search SRA via ENA portal API with rich per-study metadata."""
        requests = _check_requests()

        ena_query = f'library_strategy="RNA-Seq" AND {query}'
        if organism:
            ena_query += f' AND tax_name("{organism}")'
        library_layout = kwargs.get('library_layout')
        if library_layout and library_layout != 'Any':
            ena_query += f' AND library_layout="{library_layout.upper()}"'

        params = {
            'result': 'read_run',
            'query': ena_query,
            'fields': (
                'run_accession,study_accession,study_title,'
                'scientific_name,instrument_model,'
                'library_layout,library_source,library_strategy,'
                'library_selection,read_count,base_count,'
                'sample_accession,sample_alias,sample_title'
            ),
            'format': 'json',
            'limit': min(max_results * 10, 1000),
        }

        try:
            resp = requests.get(ENA_SEARCH_URL, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            logger.debug(f"ENA search failed ({e}), trying Entrez...")
            return self._search_via_entrez(query, organism, max_results, **kwargs)

        if not data:
            return self._search_via_entrez(query, organism, max_results, **kwargs)

        # Group runs by study
        study_map = {}
        for run in data:
            study_acc = run.get('study_accession', '')
            if not study_acc:
                continue
            if study_acc not in study_map:
                study_map[study_acc] = {
                    'title': run.get('study_title', ''),
                    'organism': run.get('scientific_name', 'unknown'),
                    'runs': [], 'platforms': set(), 'layouts': set(),
                    'selections': set(), 'total_reads': 0, 'samples': set(),
                    'sample_aliases': set(),
                }
            info = study_map[study_acc]
            info['runs'].append(run.get('run_accession', ''))
            info['platforms'].add(run.get('instrument_model', ''))
            info['layouts'].add(run.get('library_layout', ''))
            info['selections'].add(run.get('library_selection', ''))
            info['samples'].add(run.get('sample_accession', ''))
            alias = run.get('sample_alias', '')
            if alias:
                info['sample_aliases'].add(alias)
            try:
                info['total_reads'] += int(run.get('read_count', 0))
            except (ValueError, TypeError):
                pass

        results = []
        for study_acc, info in list(study_map.items())[:max_results]:
            n_runs = len(set(info['runs']))
            avg_reads = info['total_reads'] // n_runs if n_runs > 0 else 0
            avg_str = f"{avg_reads/1e6:.1f}M" if avg_reads > 1e6 else f"{avg_reads/1e3:.0f}K" if avg_reads > 1e3 else str(avg_reads)

            # Check for GSM IDs in sample aliases (indicates GEO link)
            gsm_ids = sorted({
                m for alias in info['sample_aliases']
                for m in re.findall(r'GSM\d+', alias)
            })
            gse_hint = f"GSMs: {len(gsm_ids)}" if gsm_ids else ''

            results.append(SearchResult(
                accession=study_acc,
                title=info['title'][:100],
                organism=info['organism'],
                n_samples=len(info['samples']),
                platform='; '.join(p for p in info['platforms'] if p) or 'unknown',
                description='', repository='SRA',
                gse='',
                has_geo_link=bool(gsm_ids),
                geo_hint=gse_hint,
                n_runs=n_runs,
                library_layout='; '.join(l for l in info['layouts'] if l),
                library_selection='; '.join(s for s in info['selections'] if s),
                avg_reads=avg_str,
            ))
        return results

    def _search_via_entrez(self, query, organism, max_results, **kwargs):
        """Fallback search via NCBI Entrez with study-level grouping."""
        import re
        Entrez = _check_entrez()
        Entrez.email = self._email

        search_term = f'{query} AND "rna seq"[Strategy]'
        if organism:
            search_term += f' AND {organism}[Organism]'

        try:
            handle = Entrez.esearch(
                db=SRA_DB,
                term=search_term,
                retmax=min(max_results * 3, 500),  # fetch more, group later
                sort='relevance',
            )
            search_results = Entrez.read(handle)
            handle.close()
        except Exception as e:
            raise RAPTORError(f"SRA Entrez search failed: {e}")

        id_list = search_results.get('IdList', [])
        if not id_list:
            return []

        self._wait_for_rate_limit()
        try:
            # Fetch in batches of 100
            all_summaries = []
            for i in range(0, len(id_list), 100):
                batch = id_list[i:i + 100]
                handle = Entrez.esummary(db=SRA_DB, id=','.join(batch))
                summaries = Entrez.read(handle)
                handle.close()
                all_summaries.extend(summaries)
                if i + 100 < len(id_list):
                    self._wait_for_rate_limit()
        except Exception as e:
            raise RAPTORError(f"Failed to fetch SRA summaries: {e}")

        # Parse and group by study
        study_map = {}  # study_accession → info dict

        for summary in all_summaries:
            exp_xml = summary.get('ExpXml', '')
            runs_xml = summary.get('Runs', '')

            # Parse XML fields
            def _xml_tag(xml, tag):
                pattern = f'<{tag}>(.*?)</{tag}>'
                m = re.search(pattern, xml, re.DOTALL)
                return m.group(1).strip() if m else ''

            def _xml_attr(xml, tag, attr):
                pattern = f'<{tag}[^>]*{attr}="([^"]*)"'
                m = re.search(pattern, xml)
                return m.group(1) if m else ''

            title = _xml_tag(exp_xml, 'Title')
            organism_name = _xml_tag(exp_xml, 'ScientificName')
            platform = _xml_attr(exp_xml, 'Platform', 'instrument_model')
            study_acc = _xml_attr(exp_xml, 'Study', 'acc')
            study_name = _xml_attr(exp_xml, 'Study', 'name')
            library_layout = _xml_attr(exp_xml, 'LIBRARY_LAYOUT', '').strip('</ >')
            library_source = _xml_tag(exp_xml, 'LIBRARY_SOURCE')
            library_strategy = _xml_tag(exp_xml, 'LIBRARY_STRATEGY')

            # Extract GSE if available
            gse = ''
            ext_links = re.findall(r'GSE\d+', exp_xml)
            if ext_links:
                gse = ext_links[0]

            # Count runs
            run_accs = re.findall(r'Run acc="(SRR\d+)"', runs_xml)

            # Use study accession as grouping key, fall back to GSE
            group_key = study_acc or gse or summary.get('Accession', '')
            if not group_key:
                continue

            if group_key not in study_map:
                study_map[group_key] = {
                    'accession': study_acc,
                    'gse': gse,
                    'title': study_name or title,
                    'organism': organism_name or 'unknown',
                    'platform': platform or 'various',
                    'library_source': library_source,
                    'library_strategy': library_strategy,
                    'runs': [],
                    'samples': 0,
                }

            study_map[group_key]['runs'].extend(run_accs)
            study_map[group_key]['samples'] += 1

        # Build results
        results = []
        for key, info in list(study_map.items())[:max_results]:
            display_acc = info['accession']
            if info['gse']:
                display_acc = f"{info['accession']} ({info['gse']})" if info['accession'] else info['gse']

            results.append(SearchResult(
                accession=info['accession'] or info['gse'] or key,
                title=info['title'][:100],
                organism=info['organism'],
                n_samples=info['samples'],
                platform=info['platform'],
                description='',
                repository='SRA',
                gse=info['gse'],
                has_geo_link=bool(info['gse']),
                geo_hint=info['gse'] if info['gse'] else '',
                n_runs=len(set(info['runs'])),
                library_source=info['library_source'],
            ))

        return results

    # -------------------------------------------------------------------------
    # Run table (main SRA utility)
    # -------------------------------------------------------------------------

    def get_run_table(self, study_accession: str) -> pd.DataFrame:
        """
        Get the run table for an SRA study.

        This is the most useful SRA function — it returns a DataFrame
        with all runs, their metadata, and FASTQ download links.

        Parameters
        ----------
        study_accession : str
            SRA study accession (SRP*, ERP*, DRP*).

        Returns
        -------
        pd.DataFrame
            Run table with columns: run_accession, instrument_model,
            library_layout, read_count, sample_alias, scientific_name, etc.
        """
        requests = _check_requests()

        params = {
            'accession': study_accession,
            'result': 'read_run',
            'fields': ENA_RUN_FIELDS,
            'format': 'json',
            'limit': 5000,
        }

        self._wait_for_rate_limit()
        try:
            resp = requests.get(ENA_FILEREPORT_URL, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(
                f"Failed to get run table for {study_accession}: {e}"
            )

        if not data:
            raise RAPTORError(
                f"No runs found for {study_accession}"
            )

        run_table = pd.DataFrame(data)
        run_table.index.name = 'index'

        logger.info(
            f"Retrieved {len(run_table)} runs for {study_accession}"
        )

        return run_table

    # -------------------------------------------------------------------------
    # Download commands
    # -------------------------------------------------------------------------

    def generate_download_commands(
        self,
        study_accession: str,
        output_dir: str = './fastq',
        tool: str = 'fasterq-dump',
        threads: int = 4,
    ) -> List[str]:
        """
        Generate shell commands for downloading FASTQ files.

        Parameters
        ----------
        study_accession : str
            SRA study accession.
        output_dir : str
            Directory for FASTQ output.
        tool : str
            Download tool: 'fasterq-dump' or 'fastq-dump'.
        threads : int
            Number of threads for fasterq-dump.

        Returns
        -------
        list of str
            Shell commands, one per run.
        """
        run_table = self.get_run_table(study_accession)

        if 'run_accession' not in run_table.columns:
            raise RAPTORError("Run table missing run_accession column")

        commands = []
        for _, row in run_table.iterrows():
            run_acc = row['run_accession']

            if tool == 'fasterq-dump':
                cmd = (
                    f"fasterq-dump --split-files --threads {threads} "
                    f"--outdir {output_dir} {run_acc}"
                )
            else:
                cmd = (
                    f"fastq-dump --split-files --gzip "
                    f"--outdir {output_dir} {run_acc}"
                )
            commands.append(cmd)

        return commands

    def export_run_table(
        self,
        study_accession: str,
        output_path: Union[str, Path],
    ) -> Path:
        """
        Export run table as CSV for use with RAPTOR pipelines.

        Parameters
        ----------
        study_accession : str
            SRA study accession.
        output_path : str or Path
            Output CSV path.

        Returns
        -------
        Path
            Path to exported file.
        """
        run_table = self.get_run_table(study_accession)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        run_table.to_csv(output_path, index=False, encoding='utf-8')
        logger.info(f"Exported run table to {output_path}")
        return output_path

    def download_fastqs(
        self,
        study_accession: str,
        output_dir: Union[str, Path] = './fastq',
        tool: str = 'fasterq-dump',
        threads: int = 4,
        max_runs: Optional[int] = None,
        dry_run: bool = False,
    ) -> Dict[str, Any]:
        """
        Download FASTQ files for an SRA study.

        Executes fasterq-dump or fastq-dump for each run. Requires
        SRA Toolkit to be installed and on PATH.

        Parameters
        ----------
        study_accession : str
            SRA study accession.
        output_dir : str or Path
            Directory for FASTQ output.
        tool : str
            'fasterq-dump' (faster, uncompressed) or 'fastq-dump' (gzipped).
        threads : int
            Threads for fasterq-dump.
        max_runs : int, optional
            Limit number of runs to download. None = all.
        dry_run : bool
            If True, print commands without executing.

        Returns
        -------
        dict
            Report with keys: completed, failed, skipped, commands, output_dir
        """
        import shutil
        import subprocess

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Check tool availability
        tool_path = shutil.which(tool)
        if tool_path is None and not dry_run:
            raise RAPTORError(
                f"{tool} not found on PATH. Install SRA Toolkit:\n"
                f"  conda install -c bioconda sra-tools\n"
                f"  OR download from: https://github.com/ncbi/sra-tools"
            )

        # Get run table
        run_table = self.get_run_table(study_accession)
        if 'run_accession' not in run_table.columns:
            raise RAPTORError("Run table missing run_accession column")

        runs = run_table['run_accession'].tolist()
        if max_runs:
            runs = runs[:max_runs]

        logger.info(
            f"{'[DRY RUN] ' if dry_run else ''}"
            f"Downloading {len(runs)} runs to {output_dir}"
        )

        report = {
            'completed': [],
            'failed': [],
            'skipped': [],
            'commands': [],
            'output_dir': str(output_dir),
        }

        for i, run_acc in enumerate(runs, 1):
            # Check if already downloaded
            existing = list(output_dir.glob(f"{run_acc}*fastq*"))
            if existing:
                logger.info(f"[{i}/{len(runs)}] {run_acc}: already exists, skipping")
                report['skipped'].append(run_acc)
                continue

            # Build command
            if tool == 'fasterq-dump':
                cmd = [
                    tool, '--split-files',
                    '--threads', str(threads),
                    '--outdir', str(output_dir),
                    run_acc,
                ]
            else:
                cmd = [
                    tool, '--split-files', '--gzip',
                    '--outdir', str(output_dir),
                    run_acc,
                ]

            cmd_str = ' '.join(cmd)
            report['commands'].append(cmd_str)

            if dry_run:
                logger.info(f"[{i}/{len(runs)}] [DRY RUN] {cmd_str}")
                continue

            # Execute
            logger.info(f"[{i}/{len(runs)}] Downloading {run_acc}...")
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=3600,  # 1 hour per run
                )

                if result.returncode == 0:
                    report['completed'].append(run_acc)
                    logger.info(f"  Completed: {run_acc}")
                else:
                    report['failed'].append({
                        'run': run_acc,
                        'error': result.stderr[:200],
                    })
                    logger.warning(f"  Failed: {run_acc}: {result.stderr[:100]}")

            except subprocess.TimeoutExpired:
                report['failed'].append({
                    'run': run_acc,
                    'error': 'Timeout (>1 hour)',
                })
                logger.warning(f"  Timeout: {run_acc}")
            except Exception as e:
                report['failed'].append({
                    'run': run_acc,
                    'error': str(e),
                })
                logger.warning(f"  Error: {run_acc}: {e}")

        # Summary
        logger.info(
            f"Download complete: {len(report['completed'])} completed, "
            f"{len(report['failed'])} failed, {len(report['skipped'])} skipped"
        )

        return report

    # -------------------------------------------------------------------------
    # Download (AcquiredDataset)
    # -------------------------------------------------------------------------

    def _download_api(
        self,
        accession: str,
        data_type: str = 'raw_counts',
    ) -> AcquiredDataset:
        """
        Create an AcquiredDataset from SRA run metadata.

        Since SRA stores raw reads, not counts, this creates a metadata-only
        dataset with the run table as the 'counts' placeholder. For actual
        count data, use RAPTOR Module 1 pipelines after downloading FASTQs.
        """
        run_table = self.get_run_table(accession)

        # Build a metadata DataFrame
        metadata_cols = [
            'run_accession', 'sample_alias', 'sample_title',
            'scientific_name', 'instrument_model',
            'library_strategy', 'library_layout',
            'read_count', 'base_count',
        ]
        available_cols = [c for c in metadata_cols if c in run_table.columns]
        metadata = run_table[available_cols].copy()

        if 'run_accession' in metadata.columns:
            metadata = metadata.set_index('run_accession')
        metadata.index.name = 'sample_id'

        # Create a placeholder count matrix (run accessions x basic stats)
        count_data = {}
        if 'read_count' in run_table.columns:
            for _, row in run_table.iterrows():
                run_acc = row.get('run_accession', '')
                count_data[run_acc] = {'read_count': int(row.get('read_count', 0))}

        if count_data:
            counts_df = pd.DataFrame(count_data).T
            counts_df.index.name = 'gene_id'
        else:
            counts_df = pd.DataFrame(
                {'placeholder': 1},
                index=metadata.index,
            )
            counts_df.index.name = 'gene_id'

        organism = 'unknown'
        if 'scientific_name' in run_table.columns:
            organism = run_table['scientific_name'].mode().iloc[0] if not run_table['scientific_name'].empty else 'unknown'

        source_info = {
            'repository': 'SRA',
            'accession': accession,
            'organism': organism,
            'platform': run_table['instrument_model'].mode().iloc[0] if 'instrument_model' in run_table.columns and not run_table['instrument_model'].empty else 'unknown',
            'download_date': datetime.now().isoformat(),
            'data_type': 'run_metadata',
            'description': f'SRA study {accession} with {len(run_table)} runs',
            'n_runs': len(run_table),
            'note': (
                'SRA provides raw reads, not count matrices. '
                'Use RAPTOR Module 1 pipelines to quantify after '
                'downloading FASTQs.'
            ),
        }

        # Try to find linked GEO Series (GSE) from GSM IDs in run table
        gsm_ids = self._extract_gsm_ids(run_table)
        if gsm_ids:
            source_info['gsm_ids'] = gsm_ids
            try:
                gse = self.find_linked_gse(accession, run_table=run_table)
                if gse:
                    source_info['gse'] = gse
                    source_info['note'] = (
                        f'This SRA study has a linked GEO entry: {gse}. '
                        f'You can download processed counts from GEO '
                        f'instead of quantifying raw FASTQs.'
                    )
            except Exception as e:
                logger.debug(f"GSE lookup failed for {accession}: {e}")

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info=source_info,
            gene_id_type='unknown',
        )

    # -------------------------------------------------------------------------
    # Info
    # -------------------------------------------------------------------------

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """Get study info from ENA."""
        requests = _check_requests()

        params = {
            'accession': accession,
            'result': 'read_study',
            'fields': 'study_accession,study_title,scientific_name,sample_count',
            'format': 'json',
        }

        try:
            resp = requests.get(ENA_SEARCH_URL, params=params, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                if data:
                    hit = data[0]
                    return {
                        'accession': hit.get('study_accession', accession),
                        'title': hit.get('study_title', ''),
                        'organism': hit.get('scientific_name', 'unknown'),
                        'n_samples': int(hit.get('sample_count', 0)),
                        'repository': 'SRA',
                    }
        except Exception:
            pass

        return None

    # -------------------------------------------------------------------------
    # GEO cross-referencing
    # -------------------------------------------------------------------------

    @staticmethod
    def _extract_gsm_ids(
        run_table: pd.DataFrame,
        columns: Optional[List[str]] = None,
    ) -> List[str]:
        """
        Extract GSM accession IDs from a run table.

        Many SRA studies are also deposited in GEO, and the GSM IDs
        appear in the sample_alias or sample_title columns.

        Parameters
        ----------
        run_table : pd.DataFrame
            SRA run table.
        columns : list of str, optional
            Columns to search. Defaults to sample_alias, sample_accession,
            sample_title.

        Returns
        -------
        list of str
            Unique GSM accession IDs found.
        """
        if columns is None:
            columns = ['sample_alias', 'sample_accession', 'sample_title']

        gsm_ids = set()
        for col in columns:
            if col in run_table.columns:
                matches = run_table[col].astype(str).str.findall(r'GSM\d+')
                for match_list in matches:
                    gsm_ids.update(match_list)

        return sorted(gsm_ids)

    def _gse_from_gsm(self, gsm_id: str) -> str:
        """
        Look up the parent GSE for a GSM accession via NCBI Entrez.

        Parameters
        ----------
        gsm_id : str
            A GEO sample accession (e.g., 'GSM8716844').

        Returns
        -------
        str
            The parent GSE accession, or '' if not found.
        """
        try:
            Entrez = _check_entrez()
            Entrez.email = self._email

            self._wait_for_rate_limit()
            handle = Entrez.esearch(
                db='gds',
                term=f'{gsm_id}[ACCN]',
                retmax=5,
            )
            result = Entrez.read(handle)
            handle.close()

            id_list = result.get('IdList', [])
            if not id_list:
                return ''

            self._wait_for_rate_limit()
            handle = Entrez.esummary(db='gds', id=','.join(id_list[:5]))
            summaries = Entrez.read(handle)
            handle.close()

            for summary in summaries:
                accession = summary.get('Accession', '')
                if accession.startswith('GSE'):
                    return accession

        except Exception as e:
            logger.debug(f"GSE lookup from {gsm_id} failed: {e}")

        return ''

    def find_linked_gse(
        self,
        study_accession: str,
        run_table: Optional[pd.DataFrame] = None,
    ) -> str:
        """
        Find a linked GEO Series (GSE) for an SRA study.

        Many SRA studies have a companion GEO entry with processed
        count matrices. This method finds the GSE accession by:

        1. Checking the run table for GSM sample accessions
        2. Looking up the parent GSE via NCBI Entrez

        Parameters
        ----------
        study_accession : str
            SRA study accession (e.g., 'SRP555825').
        run_table : pd.DataFrame, optional
            Pre-fetched run table. If None, fetches it.

        Returns
        -------
        str
            The linked GSE accession, or '' if none found.

        Examples
        --------
        >>> sra = SRAConnector(email='user@example.com')
        >>> gse = sra.find_linked_gse('SRP555825')
        >>> print(gse)  # e.g., 'GSE306761'
        """
        if run_table is None:
            try:
                run_table = self.get_run_table(study_accession)
            except Exception as e:
                logger.debug(f"Could not fetch run table for {study_accession}: {e}")
                return ''

        # Extract GSM IDs from the run table
        gsm_ids = self._extract_gsm_ids(run_table)
        if not gsm_ids:
            logger.debug(f"No GSM IDs found in run table for {study_accession}")
            return ''

        logger.info(
            f"Found {len(gsm_ids)} GSM IDs in {study_accession}, "
            f"looking up parent GSE..."
        )

        # Try the first GSM — all should belong to the same GSE
        gse = self._gse_from_gsm(gsm_ids[0])

        if gse:
            logger.info(f"{study_accession} → {gse}")
        else:
            # Try a second GSM as fallback
            if len(gsm_ids) > 1:
                gse = self._gse_from_gsm(gsm_ids[1])
                if gse:
                    logger.info(f"{study_accession} → {gse} (2nd GSM)")

        return gse

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"SRAConnector("
            f"email='{self._email[:20]}...', "
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )
