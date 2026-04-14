"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - TCGA/GDC Connector

Comprehensive connector for the NCI Genomic Data Commons (GDC) REST API.
Covers TCGA (33 cancer types), TARGET (pediatric), and other GDC programs.

Supported Data Types
--------------------
- Gene Expression Quantification (STAR counts / TPM / FPKM)
- miRNA Expression Quantification
- Copy Number Variation (gene-level and segment-level)
- DNA Methylation (beta values)
- Somatic Mutations (open-access MAF)
- Clinical metadata (demographics, survival, staging, treatments, exposures)

Features
--------
- Search across all GDC projects by keyword, site, disease, or program
- Sample type filtering (Primary Tumor, Solid Tissue Normal, etc.)
- Batch download for large projects (tar.gz bundles via POST)
- Rich clinical metadata with survival data for Kaplan-Meier analysis
- TCGA barcode parsing and sample type classification
- Paired tumor/normal matching within the same patient
- Multi-project download (e.g., all lung cancers: LUAD + LUSC)
- Manifest generation for GDC Data Transfer Tool
- Gene-level mutation frequency via GDC analysis endpoints
- Quantification column selection (unstranded, TPM, FPKM, FPKM-UQ)
- Progress reporting for dashboard integration
- GDC API status check and data release version

No authentication required for open-access data.

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
import tarfile
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Union, Callable, Tuple

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
# Constants — GDC API Endpoints
# =============================================================================

GDC_API_BASE = "https://api.gdc.cancer.gov"
GDC_FILES_ENDPOINT = f"{GDC_API_BASE}/files"
GDC_CASES_ENDPOINT = f"{GDC_API_BASE}/cases"
GDC_PROJECTS_ENDPOINT = f"{GDC_API_BASE}/projects"
GDC_DATA_ENDPOINT = f"{GDC_API_BASE}/data"
GDC_STATUS_ENDPOINT = f"{GDC_API_BASE}/status"
GDC_MANIFEST_ENDPOINT = f"{GDC_API_BASE}/manifest"
GDC_SSM_ENDPOINT = f"{GDC_API_BASE}/ssms"
GDC_SSM_OCCUR_ENDPOINT = f"{GDC_API_BASE}/ssm_occurrences"
GDC_CNV_ENDPOINT = f"{GDC_API_BASE}/cnvs"
GDC_CNV_OCCUR_ENDPOINT = f"{GDC_API_BASE}/cnv_occurrences"
GDC_GENES_ENDPOINT = f"{GDC_API_BASE}/genes"
GDC_ANALYSIS_ENDPOINT = f"{GDC_API_BASE}/analysis"
GDC_SURVIVAL_ENDPOINT = f"{GDC_ANALYSIS_ENDPOINT}/survival"
GDC_TOP_MUTATED_ENDPOINT = f"{GDC_ANALYSIS_ENDPOINT}/top_mutated_genes_by_project"
GDC_TOP_CASES_ENDPOINT = f"{GDC_ANALYSIS_ENDPOINT}/top_cases_counts_by_genes"


# =============================================================================
# Constants — All TCGA projects (33 cancer types)
# =============================================================================

TCGA_PROJECTS = {
    'TCGA-ACC':  'Adrenocortical carcinoma',
    'TCGA-BLCA': 'Bladder urothelial carcinoma',
    'TCGA-BRCA': 'Breast invasive carcinoma',
    'TCGA-CESC': 'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
    'TCGA-CHOL': 'Cholangiocarcinoma',
    'TCGA-COAD': 'Colon adenocarcinoma',
    'TCGA-DLBC': 'Lymphoid neoplasm diffuse large B-cell lymphoma',
    'TCGA-ESCA': 'Esophageal carcinoma',
    'TCGA-GBM':  'Glioblastoma multiforme',
    'TCGA-HNSC': 'Head and neck squamous cell carcinoma',
    'TCGA-KICH': 'Kidney chromophobe',
    'TCGA-KIRC': 'Kidney renal clear cell carcinoma',
    'TCGA-KIRP': 'Kidney renal papillary cell carcinoma',
    'TCGA-LAML': 'Acute myeloid leukemia',
    'TCGA-LGG':  'Brain lower grade glioma',
    'TCGA-LIHC': 'Liver hepatocellular carcinoma',
    'TCGA-LUAD': 'Lung adenocarcinoma',
    'TCGA-LUSC': 'Lung squamous cell carcinoma',
    'TCGA-MESO': 'Mesothelioma',
    'TCGA-OV':   'Ovarian serous cystadenocarcinoma',
    'TCGA-PAAD': 'Pancreatic adenocarcinoma',
    'TCGA-PCPG': 'Pheochromocytoma and paraganglioma',
    'TCGA-PRAD': 'Prostate adenocarcinoma',
    'TCGA-READ': 'Rectum adenocarcinoma',
    'TCGA-SARC': 'Sarcoma',
    'TCGA-SKCM': 'Skin cutaneous melanoma',
    'TCGA-STAD': 'Stomach adenocarcinoma',
    'TCGA-TGCT': 'Testicular germ cell tumors',
    'TCGA-THCA': 'Thyroid carcinoma',
    'TCGA-THYM': 'Thymoma',
    'TCGA-UCEC': 'Uterine corpus endometrial carcinoma',
    'TCGA-UCS':  'Uterine carcinosarcoma',
    'TCGA-UVM':  'Uveal melanoma',
}

# =============================================================================
# Constants — TARGET projects (pediatric cancers)
# =============================================================================

TARGET_PROJECTS = {
    'TARGET-AML':    'Acute myeloid leukemia',
    'TARGET-ALL-P2': 'Acute lymphoblastic leukemia - Phase II',
    'TARGET-ALL-P3': 'Acute lymphoblastic leukemia - Phase III',
    'TARGET-NBL':    'Neuroblastoma',
    'TARGET-OS':     'Osteosarcoma',
    'TARGET-WT':     'Wilms tumor',
    'TARGET-RT':     'Rhabdoid tumors',
}

# Combined
ALL_GDC_PROJECTS = {**TCGA_PROJECTS, **TARGET_PROJECTS}


# =============================================================================
# Constants — TCGA barcode sample type codes
# =============================================================================

TCGA_SAMPLE_TYPES = {
    '01': 'Primary Solid Tumor',
    '02': 'Recurrent Solid Tumor',
    '03': 'Primary Blood Derived Cancer - Peripheral Blood',
    '04': 'Recurrent Blood Derived Cancer - Bone Marrow',
    '05': 'Additional - New Primary',
    '06': 'Metastatic',
    '07': 'Additional Metastatic',
    '08': 'Human Tumor Original Cells',
    '09': 'Primary Blood Derived Cancer - Bone Marrow',
    '10': 'Blood Derived Normal',
    '11': 'Solid Tissue Normal',
    '12': 'Buccal Cell Normal',
    '13': 'EBV Immortalized Normal',
    '14': 'Bone Marrow Normal',
    '20': 'Control Analyte',
    '40': 'Recurrent Blood Derived Cancer - Peripheral Blood',
    '50': 'Cell Lines',
    '60': 'Primary Xenograft Tissue',
    '61': 'Cell Line Derived Xenograft Tissue',
}

SAMPLE_TYPE_GROUPS = {
    'Tumor': [
        'Primary Tumor', 'Recurrent Tumor',
        'Primary Blood Derived Cancer - Peripheral Blood',
        'Primary Blood Derived Cancer - Bone Marrow',
        'Metastatic', 'Additional Metastatic',
        'Additional - New Primary',
        'Recurrent Blood Derived Cancer - Peripheral Blood',
        'Recurrent Blood Derived Cancer - Bone Marrow',
    ],
    'Normal': [
        'Solid Tissue Normal', 'Blood Derived Normal',
        'Buccal Cell Normal', 'EBV Immortalized Normal',
        'Bone Marrow Normal',
    ],
    'Primary Tumor': [
        'Primary Tumor',
        'Primary Blood Derived Cancer - Peripheral Blood',
        'Primary Blood Derived Cancer - Bone Marrow',
    ],
}


# =============================================================================
# Constants — Quantification columns in STAR-Counts files
# =============================================================================

QUANT_COLUMNS = {
    'unstranded':         'Raw counts (unstranded)',
    'stranded_first':     'Raw counts (stranded first-read)',
    'stranded_second':    'Raw counts (stranded second-read)',
    'tpm_unstranded':     'TPM (unstranded)',
    'fpkm_unstranded':    'FPKM (unstranded)',
    'fpkm_uq_unstranded': 'FPKM upper-quartile (unstranded)',
}


# =============================================================================
# Constants — GDC data categories and types for multi-omic queries
# =============================================================================

GDC_DATA_TYPES = {
    'gene_expression': {
        'data_category': 'Transcriptome Profiling',
        'data_type': 'Gene Expression Quantification',
        'workflow_type': 'STAR - Counts',
    },
    'mirna_expression': {
        'data_category': 'Transcriptome Profiling',
        'data_type': 'miRNA Expression Quantification',
    },
    'isoform_expression': {
        'data_category': 'Transcriptome Profiling',
        'data_type': 'Isoform Expression Quantification',
    },
    'cnv_segment': {
        'data_category': 'Copy Number Variation',
        'data_type': 'Copy Number Segment',
    },
    'cnv_masked': {
        'data_category': 'Copy Number Variation',
        'data_type': 'Masked Copy Number Segment',
    },
    'cnv_gene_level': {
        'data_category': 'Copy Number Variation',
        'data_type': 'Gene Level Copy Number',
    },
    'methylation': {
        'data_category': 'DNA Methylation',
        'data_type': 'Methylation Beta Value',
    },
    'somatic_mutation': {
        'data_category': 'Simple Nucleotide Variation',
        'data_type': 'Masked Somatic Mutation',
    },
    'protein_expression': {
        'data_category': 'Proteome Profiling',
        'data_type': 'Protein Expression Quantification',
    },
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
# TCGA Barcode Utilities
# =============================================================================

def parse_tcga_barcode(barcode: str) -> Dict[str, str]:
    """
    Parse a TCGA barcode into its component parts.

    TCGA barcodes follow the format:
        TCGA-XX-XXXX-SSV-PPP-AAAA-BB

    Fields:
        XX     = Tissue Source Site (e.g., 'BH' = Christiana Healthcare)
        XXXX   = Participant (4-char alphanumeric)
        SS     = Sample type (01=Primary Tumor, 11=Normal, etc.)
        V      = Vial (A, B, ...)
        PPP    = Portion + analyte
        AAAA   = Plate
        BB     = Center

    Parameters
    ----------
    barcode : str
        Full or partial TCGA barcode.

    Returns
    -------
    dict
        Parsed components.
    """
    parts = barcode.split('-')
    result = {
        'full_barcode': barcode,
        'project': parts[0] if len(parts) > 0 else '',
        'tss': parts[1] if len(parts) > 1 else '',
        'participant': parts[2] if len(parts) > 2 else '',
    }

    result['patient_barcode'] = '-'.join(parts[:3]) if len(parts) >= 3 else barcode

    if len(parts) > 3:
        sample_part = parts[3]
        result['sample_type_code'] = sample_part[:2]
        result['sample_type'] = TCGA_SAMPLE_TYPES.get(
            sample_part[:2], f'Unknown ({sample_part[:2]})'
        )
        result['vial'] = sample_part[2:] if len(sample_part) > 2 else ''
        result['sample_barcode'] = '-'.join(parts[:4])
    else:
        result['sample_type_code'] = ''
        result['sample_type'] = 'Unknown'
        result['vial'] = ''
        result['sample_barcode'] = barcode

    if len(parts) > 4:
        portion_part = parts[4]
        result['portion'] = portion_part[:2] if len(portion_part) >= 2 else portion_part
        result['analyte'] = portion_part[2:] if len(portion_part) > 2 else ''
    else:
        result['portion'] = ''
        result['analyte'] = ''

    result['plate'] = parts[5] if len(parts) > 5 else ''
    result['center'] = parts[6] if len(parts) > 6 else ''

    return result


def classify_sample_type(sample_type_str: str) -> str:
    """Classify a GDC sample_type string into 'Tumor', 'Normal', or 'Other'."""
    st = sample_type_str.lower()
    if any(kw in st for kw in ['tumor', 'cancer', 'metastatic']):
        return 'Tumor'
    elif 'normal' in st:
        return 'Normal'
    return 'Other'


def is_tumor_barcode(barcode: str) -> bool:
    """Check if a TCGA barcode represents a tumor sample (codes 01-09)."""
    parsed = parse_tcga_barcode(barcode)
    code = parsed.get('sample_type_code', '')
    return code.isdigit() and 1 <= int(code) <= 9


def is_normal_barcode(barcode: str) -> bool:
    """Check if a TCGA barcode represents a normal sample (codes 10-14)."""
    parsed = parse_tcga_barcode(barcode)
    code = parsed.get('sample_type_code', '')
    return code.isdigit() and 10 <= int(code) <= 14


# =============================================================================
# TCGAConnector
# =============================================================================

class TCGAConnector(BaseConnector):
    """
    Comprehensive connector for TCGA/TARGET data via the GDC REST API.

    Downloads open-access RNA-seq, miRNA, CNV, methylation, and mutation
    data from the Genomic Data Commons. Includes rich clinical metadata
    with survival, staging, treatments, and exposure data.

    Parameters
    ----------
    cache : bool
        Enable disk caching.
    cache_dir : str or Path, optional
        Custom cache directory.
    workflow_type : str
        Quantification workflow: 'STAR - Counts' or 'HTSeq - Counts'.
    quant_column : str
        Default column from STAR count files: 'unstranded', 'tpm_unstranded',
        'fpkm_unstranded', 'fpkm_uq_unstranded', 'stranded_first',
        'stranded_second'.

    Examples
    --------
    >>> tcga = TCGAConnector()
    >>> results = tcga.search("liver hepatocellular")
    >>> dataset = tcga.download("TCGA-LIHC")

    >>> # Preview before downloading
    >>> info = tcga.get_study_info("TCGA-BRCA")
    >>> types = tcga.get_sample_types("TCGA-BRCA")

    >>> # Download only tumor samples
    >>> dataset = tcga.download("TCGA-LIHC", sample_types=["Primary Tumor"])

    >>> # Multi-project download
    >>> lung = tcga.download_multiple_projects(["TCGA-LUAD", "TCGA-LUSC"])

    >>> # Paired tumor/normal
    >>> pairs = tcga.get_paired_samples("TCGA-LIHC")

    >>> # miRNA expression
    >>> mirna = tcga.download_mirna("TCGA-BRCA")

    >>> # Copy number variation
    >>> cnv = tcga.download_cnv("TCGA-BRCA")

    >>> # Somatic mutations
    >>> maf = tcga.download_mutations("TCGA-BRCA")

    >>> # Methylation
    >>> methyl = tcga.download_methylation("TCGA-BRCA")

    >>> # Protein expression (RPPA)
    >>> rppa = tcga.download_rppa("TCGA-BRCA")

    >>> # Survival data
    >>> surv = tcga.get_survival_data("TCGA-BRCA")

    >>> # Gene mutation frequencies
    >>> freq = tcga.get_gene_mutation_frequency("TCGA-BRCA")

    >>> # GDC transfer tool manifest
    >>> tcga.generate_manifest("TCGA-BRCA", output_path="manifest.tsv")
    """

    REPOSITORY_NAME = 'TCGA'

    def __init__(
        self,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
        workflow_type: str = 'STAR - Counts',
        quant_column: str = 'unstranded',
    ):
        super().__init__(
            cache=cache,
            cache_dir=cache_dir,
            rate_limit=0.2,
        )
        self._workflow_type = workflow_type
        self._quant_column = quant_column

        if quant_column not in QUANT_COLUMNS:
            raise ValidationError(
                'quant_column',
                f"Invalid quantification column: '{quant_column}'",
                f"Must be one of: {list(QUANT_COLUMNS.keys())}"
            )

    # =========================================================================
    # SEARCH
    # =========================================================================

    def _search_api(
        self,
        query: str,
        organism: Optional[str] = None,
        max_results: int = 20,
        **kwargs,
    ) -> List[SearchResult]:
        """Search GDC projects matching the query."""
        requests = _check_requests()

        programs = kwargs.get('programs', ['TCGA', 'TARGET'])
        if isinstance(programs, str):
            programs = [programs]

        params = {
            'filters': json.dumps({
                'op': 'and',
                'content': [{
                    'op': 'in',
                    'content': {
                        'field': 'program.name',
                        'value': programs,
                    }
                }]
            }),
            'fields': (
                'project_id,name,primary_site,disease_type,'
                'summary.case_count,summary.file_count,'
                'summary.data_categories.data_category,'
                'summary.data_categories.file_count,'
                'summary.experimental_strategies.experimental_strategy,'
                'summary.experimental_strategies.file_count'
            ),
            'size': 200,
        }

        try:
            resp = requests.get(GDC_PROJECTS_ENDPOINT, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(f"GDC project search failed: {e}")

        query_lower = query.lower()
        query_terms = query_lower.split()
        results = []

        for hit in data.get('data', {}).get('hits', []):
            project_id = hit.get('project_id', '')
            name = hit.get('name', '')
            primary_site = hit.get('primary_site', ['unknown'])
            disease_type = hit.get('disease_type', ['unknown'])

            if isinstance(primary_site, list):
                primary_site_str = '; '.join(primary_site)
            else:
                primary_site_str = str(primary_site)
            if isinstance(disease_type, list):
                disease_type_str = '; '.join(disease_type)
            else:
                disease_type_str = str(disease_type)

            searchable = f"{project_id} {name} {primary_site_str} {disease_type_str}".lower()

            if (query_lower in searchable
                    or all(term in searchable for term in query_terms)):
                case_count = hit.get('summary', {}).get('case_count', 0)

                strategies = hit.get('summary', {}).get('experimental_strategies', [])
                has_rnaseq = any(s.get('experimental_strategy') == 'RNA-Seq' for s in strategies)
                rnaseq_files = sum(s.get('file_count', 0) for s in strategies
                                   if s.get('experimental_strategy') == 'RNA-Seq')

                data_cats = hit.get('summary', {}).get('data_categories', [])
                available_data = [d['data_category'] for d in data_cats]

                results.append(SearchResult(
                    accession=project_id,
                    title=name,
                    organism='Homo sapiens',
                    n_samples=case_count,
                    platform='Illumina HiSeq/NovaSeq',
                    description=f"{disease_type_str} | Site: {primary_site_str}",
                    repository='TCGA',
                    primary_site=primary_site_str,
                    disease_type=disease_type_str,
                    has_rnaseq=has_rnaseq,
                    rnaseq_files=rnaseq_files,
                    available_data=available_data,
                ))

        # Fallback to static dictionaries
        for proj_dict in [TCGA_PROJECTS, TARGET_PROJECTS]:
            for proj_id, desc in proj_dict.items():
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

    # =========================================================================
    # STUDY INFO & SAMPLE TYPES (preview without downloading)
    # =========================================================================

    def get_study_info(self, accession: str) -> Dict[str, Any]:
        """
        Get comprehensive study information without downloading data.

        Returns
        -------
        dict
            Keys: accession, title, organism, n_samples, primary_site,
            disease_type, data_categories, experimental_strategies,
            sample_type_counts, clinical_summary, samples (DataFrame),
            n_cases, n_files, platform_id, platform_name, data_type.
        """
        requests = _check_requests()

        info = {
            'accession': accession,
            'title': ALL_GDC_PROJECTS.get(accession, accession),
            'organism': 'Homo sapiens',
            'repository': 'TCGA',
        }

        # ---- Project-level summary ----
        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                f"{GDC_PROJECTS_ENDPOINT}/{accession}",
                params={
                    'fields': (
                        'name,primary_site,disease_type,'
                        'summary.case_count,summary.file_count,'
                        'summary.data_categories.data_category,'
                        'summary.data_categories.file_count,'
                        'summary.experimental_strategies.experimental_strategy,'
                        'summary.experimental_strategies.file_count'
                    ),
                },
                timeout=30,
            )
            resp.raise_for_status()
            proj = resp.json().get('data', {})

            info['title'] = proj.get('name', info['title'])
            info['primary_site'] = proj.get('primary_site', [])
            info['disease_type'] = proj.get('disease_type', [])
            info['n_samples'] = proj.get('summary', {}).get('case_count', 0)
            info['n_files'] = proj.get('summary', {}).get('file_count', 0)

            info['data_categories'] = {
                d['data_category']: d['file_count']
                for d in proj.get('summary', {}).get('data_categories', [])
            }
            info['experimental_strategies'] = {
                s['experimental_strategy']: s['file_count']
                for s in proj.get('summary', {}).get('experimental_strategies', [])
            }
        except Exception as e:
            logger.warning(f"Project info query failed: {e}")

        # ---- Cases with clinical data ----
        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_CASES_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': '=',
                        'content': {
                            'field': 'project.project_id',
                            'value': accession,
                        }
                    }),
                    'fields': (
                        'submitter_id,'
                        'samples.sample_type,samples.tissue_type,'
                        'samples.submitter_id,'
                        'demographic.gender,demographic.race,'
                        'demographic.ethnicity,demographic.vital_status,'
                        'demographic.days_to_death,demographic.year_of_birth,'
                        'diagnoses.primary_diagnosis,diagnoses.tumor_stage,'
                        'diagnoses.age_at_diagnosis,'
                        'diagnoses.classification_of_tumor,'
                        'diagnoses.tissue_or_organ_of_origin'
                    ),
                    'size': 5000,
                },
                timeout=60,
            )
            resp.raise_for_status()
            cases_data = resp.json().get('data', {}).get('hits', [])
        except Exception as e:
            logger.warning(f"Cases query failed: {e}")
            cases_data = []

        # Build sample type counts and sample records
        sample_type_counts = {}
        sample_records = []

        for case in cases_data:
            sub_id = case.get('submitter_id', '')
            demo = case.get('demographic', {}) or {}
            if not isinstance(demo, dict):
                demo = {}
            diags = case.get('diagnoses', [{}])
            diag = diags[0] if diags and isinstance(diags[0], dict) else {}

            gender = demo.get('gender', 'unknown')
            vital = demo.get('vital_status', 'unknown')
            race = demo.get('race', 'unknown')
            dtd = demo.get('days_to_death')
            age = diag.get('age_at_diagnosis')

            for sample in (case.get('samples') or []):
                if not isinstance(sample, dict):
                    continue
                st = sample.get('sample_type', 'unknown')
                sample_type_counts[st] = sample_type_counts.get(st, 0) + 1
                sample_records.append({
                    'case_id': sub_id,
                    'sample_id': sample.get('submitter_id', ''),
                    'sample_type': st,
                    'tissue_type': sample.get('tissue_type', 'unknown'),
                    'gender': gender,
                    'vital_status': vital,
                    'race': race,
                    'primary_diagnosis': diag.get('primary_diagnosis', 'unknown'),
                    'tumor_stage': diag.get('tumor_stage', 'unknown'),
                    'age_at_diagnosis_days': age,
                    'days_to_death': dtd,
                })

        info['sample_type_counts'] = sample_type_counts
        info['n_cases'] = len(cases_data)

        # ---- Clinical summary statistics ----
        clinical_summary = {}

        genders = [c.get('demographic', {}).get('gender', 'unknown')
                   for c in cases_data if isinstance(c.get('demographic'), dict)]
        if genders:
            gc = {}
            for g in genders:
                gc[g] = gc.get(g, 0) + 1
            clinical_summary['gender'] = gc

        vitals = [c.get('demographic', {}).get('vital_status', 'unknown')
                  for c in cases_data if isinstance(c.get('demographic'), dict)]
        if vitals:
            vc = {}
            for v in vitals:
                vc[v] = vc.get(v, 0) + 1
            clinical_summary['vital_status'] = vc

        stages = []
        for c in cases_data:
            d = c.get('diagnoses', [])
            if d and isinstance(d[0], dict):
                s = d[0].get('tumor_stage', 'unknown')
                if s and s != 'not reported':
                    stages.append(s)
        if stages:
            sc = {}
            for s in stages:
                sc[s] = sc.get(s, 0) + 1
            clinical_summary['tumor_stage'] = sc

        ages_days = []
        for c in cases_data:
            d = c.get('diagnoses', [])
            if d and isinstance(d[0], dict):
                a = d[0].get('age_at_diagnosis')
                if a is not None:
                    ages_days.append(a)
        if ages_days:
            ages_y = [a / 365.25 for a in ages_days]
            clinical_summary['age_at_diagnosis'] = {
                'mean': round(np.mean(ages_y), 1),
                'median': round(np.median(ages_y), 1),
                'min': round(min(ages_y), 1),
                'max': round(max(ages_y), 1),
            }

        info['clinical_summary'] = clinical_summary
        info['samples'] = pd.DataFrame(sample_records) if sample_records else pd.DataFrame()
        if not info['samples'].empty and 'age_at_diagnosis_days' in info['samples'].columns:
            info['samples']['age_years'] = (
                pd.to_numeric(info['samples']['age_at_diagnosis_days'], errors='coerce') / 365.25
            ).round(1)

        info['data_type'] = 'RNA-seq (STAR counts)'
        info['platform_id'] = 'Illumina HiSeq/NovaSeq'
        info['platform_name'] = 'Illumina (multiple models)'
        return info

    def get_sample_types(self, accession: str) -> Dict[str, Any]:
        """
        Get sample type breakdown for a project using GDC facets.

        Returns
        -------
        dict
            type_counts, unique_types, n_samples, tumor_count, normal_count.
        """
        requests = _check_requests()
        self._wait_for_rate_limit()

        try:
            resp = requests.get(
                GDC_CASES_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': '=',
                        'content': {
                            'field': 'project.project_id',
                            'value': accession,
                        }
                    }),
                    'facets': 'samples.sample_type',
                    'size': 0,
                },
                timeout=30,
            )
            resp.raise_for_status()
            buckets = resp.json().get('data', {}).get(
                'aggregations', {}
            ).get('samples.sample_type', {}).get('buckets', [])

            type_counts = {
                b['key']: b['doc_count']
                for b in buckets if b.get('key') and b['key'] != '_missing'
            }

            tumor_count = sum(c for st, c in type_counts.items()
                              if classify_sample_type(st) == 'Tumor')
            normal_count = sum(c for st, c in type_counts.items()
                               if classify_sample_type(st) == 'Normal')

            return {
                'type_counts': type_counts,
                'unique_types': list(type_counts.keys()),
                'n_samples': sum(type_counts.values()),
                'tumor_count': tumor_count,
                'normal_count': normal_count,
            }
        except Exception as e:
            logger.warning(f"Sample types query failed: {e}")
            return {'type_counts': {}, 'unique_types': [], 'n_samples': 0,
                    'tumor_count': 0, 'normal_count': 0}

    def get_sample_data_availability(
        self,
        project_id: str,
    ) -> pd.DataFrame:
        """
        Get per-case data type availability for a project.

        Queries the GDC files endpoint to determine which cases have
        RNA-seq, miRNA, methylation, CNV, mutations, etc.

        Parameters
        ----------
        project_id : str

        Returns
        -------
        pd.DataFrame
            One row per case with columns: case_id, and boolean columns
            for each experimental strategy (RNA-Seq, miRNA-Seq,
            Methylation Array, WXS, etc.) plus sample counts.
        """
        requests = _check_requests()

        # Query files grouped by case and experimental strategy
        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_FILES_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': 'and',
                        'content': [
                            {'op': '=', 'content': {
                                'field': 'cases.project.project_id',
                                'value': project_id,
                            }},
                            {'op': '=', 'content': {
                                'field': 'access',
                                'value': 'open',
                            }},
                        ]
                    }),
                    'fields': (
                        'cases.submitter_id,'
                        'cases.samples.submitter_id,'
                        'cases.samples.sample_type,'
                        'experimental_strategy,'
                        'data_category,'
                        'data_type'
                    ),
                    'size': 10000,
                },
                timeout=60,
            )
            resp.raise_for_status()
            hits = resp.json().get('data', {}).get('hits', [])
        except Exception as e:
            logger.warning(f"Data availability query failed: {e}")
            return pd.DataFrame()

        # Build per-case, per-sample data availability
        # Structure: {case_id: {sample_id: {sample_type, strategies: set()}}}
        case_data = {}

        for hit in hits:
            strategy = hit.get('experimental_strategy', '')
            data_cat = hit.get('data_category', '')
            cases = hit.get('cases', [])

            for case in cases:
                if not isinstance(case, dict):
                    continue
                case_id = case.get('submitter_id', '')
                if not case_id:
                    continue

                if case_id not in case_data:
                    case_data[case_id] = {'strategies': set(), 'samples': {}}

                case_data[case_id]['strategies'].add(strategy if strategy else data_cat)

                for sample in (case.get('samples') or []):
                    if not isinstance(sample, dict):
                        continue
                    sid = sample.get('submitter_id', '')
                    if sid and sid not in case_data[case_id]['samples']:
                        case_data[case_id]['samples'][sid] = {
                            'sample_type': sample.get('sample_type', 'unknown'),
                        }

        if not case_data:
            return pd.DataFrame()

        # Collect all strategies (skip blanks)
        all_strategies = set()
        for cd in case_data.values():
            all_strategies.update(cd['strategies'])
        all_strategies.discard('')
        all_strategies.discard(None)
        all_strategies = sorted(all_strategies)

        # Build DataFrame
        rows = []
        for case_id, cd in sorted(case_data.items()):
            row = {'case_id': case_id}
            for strat in all_strategies:
                row[strat] = 'Yes' if strat in cd['strategies'] else ''
            row['n_samples'] = len(cd['samples'])

            # Sample type summary
            types = [s['sample_type'] for s in cd['samples'].values()]
            tumor_n = sum(1 for t in types if classify_sample_type(t) == 'Tumor')
            normal_n = sum(1 for t in types if classify_sample_type(t) == 'Normal')
            row['tumor_samples'] = tumor_n
            row['normal_samples'] = normal_n

            rows.append(row)

        avail_df = pd.DataFrame(rows)
        logger.info(
            f"{project_id}: data availability for {len(avail_df)} cases, "
            f"{len(all_strategies)} experimental strategies"
        )
        return avail_df

    # =========================================================================
    # DOWNLOAD — Gene Expression (primary)
    # =========================================================================

    def _download_api(
        self,
        accession: str,
        data_type: str = 'raw_counts',
        **kwargs,
    ) -> AcquiredDataset:
        """
        Download RNA-seq gene expression for a TCGA/TARGET project.

        Extra kwargs
        ------------
        sample_types : list of str
            Filter by sample type(s).
        quant_column : str
            Override default quantification column.
        batch_size : int
            Files per batch download (default 50).
        progress_callback : callable(current, total)
            Progress reporting for dashboard.
        """
        requests = _check_requests()

        sample_types = kwargs.get('sample_types', None)
        quant_col = kwargs.get('quant_column', self._quant_column)
        batch_size = kwargs.get('batch_size', 50)
        progress_cb = kwargs.get('progress_callback', None)

        logger.info(f"Querying GDC for {accession} RNA-seq files...")

        file_ids, file_metadata = self._query_files(
            accession, requests, sample_types=sample_types
        )

        if not file_ids:
            msg = f"No RNA-seq count files found for {accession}."
            if sample_types:
                msg += f" (filtered by: {sample_types})"
            raise RAPTORError(msg)

        logger.info(f"Found {len(file_ids)} count files for {accession}")

        if len(file_ids) > 10:
            counts_df = self._download_counts_batch(
                file_ids, file_metadata, requests,
                quant_col=quant_col, batch_size=batch_size,
                progress_callback=progress_cb,
            )
        else:
            counts_df = self._download_counts_individual(
                file_ids, file_metadata, requests,
                quant_col=quant_col, progress_callback=progress_cb,
            )

        metadata = self._fetch_clinical_metadata(accession, file_metadata, requests)

        if quant_col.startswith('tpm'):
            actual_type = 'tpm'
        elif quant_col.startswith('fpkm'):
            actual_type = 'fpkm'
        else:
            actual_type = 'raw_counts'

        source_info = {
            'repository': 'TCGA',
            'accession': accession,
            'organism': 'Homo sapiens',
            'platform': 'Illumina HiSeq/NovaSeq',
            'download_date': datetime.now().isoformat(),
            'data_type': actual_type,
            'workflow_type': self._workflow_type,
            'quant_column': quant_col,
            'description': ALL_GDC_PROJECTS.get(accession, accession),
            'n_files_downloaded': len(file_ids),
            'sample_types_filter': sample_types,
        }

        # Performance: save as Parquet for fast reload
        try:
            cache_name = f"{accession}_{quant_col}"
            if sample_types:
                cache_name += '_' + '_'.join(s.replace(' ', '') for s in sorted(sample_types))
            parquet_dir = Path.home() / '.raptor' / 'parquet_cache'
            parquet_dir.mkdir(parents=True, exist_ok=True)
            counts_df.to_parquet(
                parquet_dir / f"{cache_name}.parquet",
                engine='pyarrow', compression='snappy',
            )
            logger.info(f"Saved Parquet cache: {cache_name}.parquet")
        except Exception as e:
            logger.debug(f"Parquet cache save failed (non-critical): {e}")

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info=source_info,
            gene_id_type='ensembl',
        )

    # =========================================================================
    # DOWNLOAD — Multi-project
    # =========================================================================

    def download_multiple_projects(
        self,
        project_ids: List[str],
        sample_types: Optional[List[str]] = None,
        quant_column: Optional[str] = None,
        progress_callback: Optional[Callable] = None,
    ) -> AcquiredDataset:
        """
        Download and merge gene expression from multiple GDC projects.

        Useful for pan-cancer or related-cancer analyses (e.g., LUAD + LUSC
        for all lung cancers, or COAD + READ for colorectal).

        Parameters
        ----------
        project_ids : list of str
            GDC project IDs.
        sample_types : list of str, optional
            Filter by sample type(s) across all projects.
        quant_column : str, optional
            Override default quantification column.
        progress_callback : callable, optional
            Called with (project_index, total_projects).

        Returns
        -------
        AcquiredDataset
            Merged dataset with 'project' column in metadata.
        """
        all_counts = []
        all_metadata = []
        total_files = 0

        for i, pid in enumerate(project_ids):
            logger.info(f"[{i + 1}/{len(project_ids)}] Downloading {pid}...")

            try:
                ds = self.download(
                    pid,
                    sample_types=sample_types,
                    quant_column=quant_column or self._quant_column,
                    progress_callback=progress_callback,
                )

                # Tag metadata with project source
                meta = ds.metadata.copy()
                meta['project'] = pid
                all_metadata.append(meta)
                all_counts.append(ds.counts_df)
                total_files += ds.source_info.get('n_files_downloaded', 0)

            except Exception as e:
                logger.error(f"Failed to download {pid}: {e}")
                continue

            if progress_callback:
                progress_callback(i + 1, len(project_ids))

        if not all_counts:
            raise RAPTORError(
                f"No data downloaded from any project: {project_ids}"
            )

        # Merge on shared genes (inner join)
        merged_counts = all_counts[0]
        for df in all_counts[1:]:
            shared_genes = merged_counts.index.intersection(df.index)
            merged_counts = pd.concat(
                [merged_counts.loc[shared_genes], df.loc[shared_genes]],
                axis=1,
            )

        merged_meta = pd.concat(all_metadata)

        source_info = {
            'repository': 'TCGA',
            'accession': '+'.join(project_ids),
            'organism': 'Homo sapiens',
            'platform': 'Illumina HiSeq/NovaSeq',
            'download_date': datetime.now().isoformat(),
            'data_type': 'raw_counts',
            'description': f"Multi-project: {', '.join(project_ids)}",
            'n_files_downloaded': total_files,
            'projects': project_ids,
        }

        return AcquiredDataset(
            counts_df=merged_counts,
            metadata=merged_meta,
            source_info=source_info,
            gene_id_type='ensembl',
        )

    # =========================================================================
    # DOWNLOAD — miRNA Expression
    # =========================================================================

    def download_mirna(
        self,
        project_id: str,
        sample_types: Optional[List[str]] = None,
        progress_callback: Optional[Callable] = None,
    ) -> AcquiredDataset:
        """
        Download miRNA expression quantification for a project.

        Returns a matrix of miRNA IDs × samples with reads_per_million values.
        """
        requests = _check_requests()

        # Check cache
        cache_name = f"{project_id}_mirna"
        if sample_types:
            cache_name += '_' + '_'.join(s.replace(' ', '') for s in sorted(sample_types))
        cached_counts = self._load_cached(cache_name, fmt='parquet')
        cached_meta = self._load_cached(f"{cache_name}_meta", fmt='csv')
        if cached_counts is not None:
            logger.info(f"Loaded miRNA from cache: {cached_counts.shape}")
            if progress_callback:
                progress_callback(1, 1)
            return AcquiredDataset(
                counts_df=cached_counts,
                metadata=cached_meta if cached_meta is not None else pd.DataFrame(),
                source_info={
                    'repository': 'TCGA', 'accession': project_id,
                    'organism': 'Homo sapiens', 'data_type': 'mirna_rpm',
                    'description': f"miRNA expression — {ALL_GDC_PROJECTS.get(project_id, project_id)} (cached)",
                },
                gene_id_type='mirna',
            )

        file_ids, file_meta = self._query_generic_files(
            project_id, requests,
            data_category='Transcriptome Profiling',
            data_type='miRNA Expression Quantification',
            sample_types=sample_types,
        )

        if not file_ids:
            raise RAPTORError(f"No miRNA expression files for {project_id}")

        logger.info(f"Downloading {len(file_ids)} miRNA files for {project_id}")

        sample_data = {}
        for i, fid in enumerate(file_ids):
            content = self._download_raw_file(fid, requests)
            if content is None:
                continue

            try:
                df = pd.read_csv(io.BytesIO(content), sep='\t')
                mirna_col = next((c for c in df.columns if 'mirna_id' in c.lower()
                                  or 'miRNA_ID' in c), df.columns[0])
                rpm_col = next((c for c in df.columns if 'reads_per_million' in c.lower()),
                               None)
                if rpm_col is None:
                    rpm_col = df.columns[1] if len(df.columns) > 1 else None
                if rpm_col is None:
                    continue

                label = self._sample_label(file_meta.get(fid, {})) or fid[:12]
                sample_data[label] = pd.Series(
                    df[rpm_col].values, index=df[mirna_col].values,
                    name=label, dtype=float,
                )
            except Exception as e:
                logger.warning(f"Failed to parse miRNA file {fid}: {e}")

            if progress_callback:
                progress_callback(i + 1, len(file_ids))

        if not sample_data:
            raise RAPTORError("No miRNA data could be parsed")

        counts_df = pd.DataFrame(sample_data).fillna(0)
        counts_df.index.name = 'mirna_id'

        metadata = self._build_metadata(file_meta, keep_samples=set(counts_df.columns))
        self._cache_data(counts_df, cache_name, fmt='parquet')
        self._cache_data(metadata, f"{cache_name}_meta", fmt='csv')

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info={
                'repository': 'TCGA', 'accession': project_id,
                'organism': 'Homo sapiens', 'data_type': 'mirna_rpm',
                'download_date': datetime.now().isoformat(),
                'description': f"miRNA expression — {ALL_GDC_PROJECTS.get(project_id, project_id)}",
                'n_files_downloaded': len(sample_data),
            },
            gene_id_type='mirna',
        )

    # =========================================================================
    # DOWNLOAD — Copy Number Variation
    # =========================================================================

    def download_cnv(
        self,
        project_id: str,
        cnv_type: str = 'gene_level',
        sample_types: Optional[List[str]] = None,
        progress_callback: Optional[Callable] = None,
    ) -> AcquiredDataset:
        """
        Download copy number variation data.

        Parameters
        ----------
        project_id : str
        cnv_type : str
            'gene_level' — gene-level copy number scores (recommended)
            'segment' — segment-level (chr, start, end, num_probes, segment_mean)
            'masked' — masked somatic-only segments
        sample_types : list of str, optional
        progress_callback : callable, optional

        Returns
        -------
        AcquiredDataset
        """
        requests = _check_requests()

        # Check cache
        cache_name = f"{project_id}_cnv_{cnv_type}"
        if sample_types:
            cache_name += '_' + '_'.join(s.replace(' ', '') for s in sorted(sample_types))
        cached_counts = self._load_cached(cache_name, fmt='parquet')
        cached_meta = self._load_cached(f"{cache_name}_meta", fmt='csv')
        if cached_counts is not None:
            logger.info(f"Loaded CNV from cache: {cached_counts.shape}")
            if progress_callback:
                progress_callback(1, 1)
            return AcquiredDataset(
                counts_df=cached_counts,
                metadata=cached_meta if cached_meta is not None else pd.DataFrame(),
                source_info={
                    'repository': 'TCGA', 'accession': project_id,
                    'organism': 'Homo sapiens', 'data_type': f'cnv_{cnv_type}',
                    'description': f"CNV {cnv_type} — {project_id} (cached)",
                },
                gene_id_type='ensembl' if cnv_type == 'gene_level' else 'genomic_region',
            )

        type_map = {
            'gene_level': 'Gene Level Copy Number',
            'segment': 'Copy Number Segment',
            'masked': 'Masked Copy Number Segment',
        }
        if cnv_type not in type_map:
            raise ValidationError('cnv_type', f"Must be one of: {list(type_map.keys())}")

        file_ids, file_meta = self._query_generic_files(
            project_id, requests,
            data_category='Copy Number Variation',
            data_type=type_map[cnv_type],
            sample_types=sample_types,
        )

        if not file_ids:
            raise RAPTORError(f"No {cnv_type} CNV files for {project_id}")

        logger.info(f"Downloading {len(file_ids)} CNV ({cnv_type}) files for {project_id}")

        if cnv_type == 'gene_level':
            result = self._download_cnv_gene_level(
                file_ids, file_meta, project_id, requests, progress_callback
            )
        else:
            result = self._download_cnv_segments(
                file_ids, file_meta, project_id, requests, cnv_type, progress_callback
            )
        self._cache_data(result.counts_df, cache_name, fmt='parquet')
        if hasattr(result, 'metadata') and not result.metadata.empty:
            self._cache_data(result.metadata, f"{cache_name}_meta", fmt='csv')
        return result

    def _download_cnv_gene_level(self, file_ids, file_meta, project_id,
                                  requests, progress_callback):
        """Parse gene-level CNV into genes × samples matrix."""
        sample_data = {}
        for i, fid in enumerate(file_ids):
            content = self._download_raw_file(fid, requests)
            if content is None:
                continue
            try:
                df = pd.read_csv(io.BytesIO(content), sep='\t')
                gene_col = next((c for c in df.columns if 'gene_id' in c.lower()), df.columns[0])
                cn_col = next((c for c in df.columns if 'copy_number' in c.lower()), None)
                if cn_col is None:
                    cn_col = df.columns[-1]

                genes = df[gene_col].astype(str).str.split('.').str[0]
                label = self._sample_label(file_meta.get(fid, {})) or fid[:12]
                sample_data[label] = pd.Series(
                    df[cn_col].values, index=genes.values,
                    name=label, dtype=float,
                )
            except Exception as e:
                logger.warning(f"Failed to parse CNV file {fid}: {e}")
            if progress_callback:
                progress_callback(i + 1, len(file_ids))

        if not sample_data:
            raise RAPTORError("No gene-level CNV data could be parsed")

        counts_df = pd.DataFrame(sample_data).fillna(0)
        counts_df.index.name = 'gene_id'

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=self._build_metadata(file_meta, keep_samples=set(counts_df.columns)),
            source_info={
                'repository': 'TCGA', 'accession': project_id,
                'organism': 'Homo sapiens', 'data_type': 'cnv_gene_level',
                'download_date': datetime.now().isoformat(),
                'n_files_downloaded': len(sample_data),
            },
            gene_id_type='ensembl',
        )

    def _download_cnv_segments(self, file_ids, file_meta, project_id,
                                requests, cnv_type, progress_callback):
        """Parse segment-level CNV into a long-format DataFrame."""
        all_segments = []
        for i, fid in enumerate(file_ids):
            content = self._download_raw_file(fid, requests)
            if content is None:
                continue
            try:
                df = pd.read_csv(io.BytesIO(content), sep='\t')
                label = self._sample_label(file_meta.get(fid, {})) or fid[:12]
                df.insert(0, 'sample_id', label)
                all_segments.append(df)
            except Exception as e:
                logger.warning(f"Failed to parse CNV segment file {fid}: {e}")
            if progress_callback:
                progress_callback(i + 1, len(file_ids))

        if not all_segments:
            raise RAPTORError("No CNV segment data could be parsed")

        segments_df = pd.concat(all_segments, ignore_index=True)
        parsed_samples = set(segments_df['sample_id'].unique()) if 'sample_id' in segments_df.columns else set()

        return AcquiredDataset(
            counts_df=segments_df.set_index('sample_id') if 'sample_id' in segments_df.columns else segments_df,
            metadata=self._build_metadata(file_meta, keep_samples=parsed_samples),
            source_info={
                'repository': 'TCGA', 'accession': project_id,
                'organism': 'Homo sapiens', 'data_type': f'cnv_{cnv_type}',
                'download_date': datetime.now().isoformat(),
                'n_files_downloaded': len(all_segments),
            },
            gene_id_type='segment',
        )

    # =========================================================================
    # DOWNLOAD — DNA Methylation
    # =========================================================================

    def download_methylation(
        self,
        project_id: str,
        sample_types: Optional[List[str]] = None,
        max_files: int = 100,
        progress_callback: Optional[Callable] = None,
    ) -> AcquiredDataset:
        """
        Download DNA methylation beta values.

        Note: Methylation files are large (~480K probes per sample).
        Use max_files to limit for initial exploration.

        Parameters
        ----------
        project_id : str
        sample_types : list of str, optional
        max_files : int
            Maximum files to download (default 100).
        progress_callback : callable, optional

        Returns
        -------
        AcquiredDataset
            CpG probes × samples matrix of beta values.
        """
        requests = _check_requests()

        # Check cache
        cache_name = f"{project_id}_methylation"
        if sample_types:
            cache_name += '_' + '_'.join(s.replace(' ', '') for s in sorted(sample_types))
        cached_counts = self._load_cached(cache_name, fmt='parquet')
        cached_meta = self._load_cached(f"{cache_name}_meta", fmt='csv')
        if cached_counts is not None:
            logger.info(f"Loaded methylation from cache: {cached_counts.shape}")
            if progress_callback:
                progress_callback(1, 1)
            return AcquiredDataset(
                counts_df=cached_counts,
                metadata=cached_meta if cached_meta is not None else pd.DataFrame(),
                source_info={
                    'repository': 'TCGA', 'accession': project_id,
                    'organism': 'Homo sapiens', 'data_type': 'methylation_beta',
                    'description': f"Methylation — {project_id} (cached)",
                },
                gene_id_type='probe',
            )

        file_ids, file_meta = self._query_generic_files(
            project_id, requests,
            data_category='DNA Methylation',
            data_type='Methylation Beta Value',
            sample_types=sample_types,
        )

        if not file_ids:
            raise RAPTORError(f"No methylation files for {project_id}")

        if len(file_ids) > max_files:
            logger.warning(
                f"Limiting to {max_files} of {len(file_ids)} methylation files. "
                f"Set max_files= to download more."
            )
            file_ids = file_ids[:max_files]

        logger.info(f"Downloading {len(file_ids)} methylation files for {project_id}")

        sample_data = {}
        for i, fid in enumerate(file_ids):
            content = self._download_raw_file(fid, requests)
            if content is None:
                continue
            try:
                df = pd.read_csv(io.BytesIO(content), sep='\t')
                probe_col = next((c for c in df.columns
                                  if 'composite' in c.lower() or 'probe' in c.lower()),
                                 df.columns[0])
                beta_col = next((c for c in df.columns if 'beta' in c.lower()), None)
                if beta_col is None and len(df.columns) > 1:
                    beta_col = df.columns[1]
                if beta_col is None:
                    continue

                label = self._sample_label(file_meta.get(fid, {})) or fid[:12]
                sample_data[label] = pd.Series(
                    pd.to_numeric(df[beta_col], errors='coerce').values,
                    index=df[probe_col].values,
                    name=label, dtype=float,
                )
            except Exception as e:
                logger.warning(f"Failed to parse methylation file {fid}: {e}")
            if progress_callback:
                progress_callback(i + 1, len(file_ids))

        if not sample_data:
            raise RAPTORError("No methylation data could be parsed")

        counts_df = pd.DataFrame(sample_data)
        counts_df.index.name = 'probe_id'

        metadata = self._build_metadata(file_meta, keep_samples=set(counts_df.columns))
        self._cache_data(counts_df, cache_name, fmt='parquet')
        self._cache_data(metadata, f"{cache_name}_meta", fmt='csv')

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info={
                'repository': 'TCGA', 'accession': project_id,
                'organism': 'Homo sapiens', 'data_type': 'methylation_beta',
                'download_date': datetime.now().isoformat(),
                'n_files_downloaded': len(sample_data),
            },
            gene_id_type='probe',
        )

    # =========================================================================
    # DOWNLOAD — Somatic Mutations (MAF)
    # =========================================================================

    def download_mutations(
        self,
        project_id: str,
        progress_callback: Optional[Callable] = None,
    ) -> pd.DataFrame:
        """
        Download open-access MAF (Masked Somatic Mutation) files.

        Returns a combined DataFrame in standard MAF format with columns:
        Hugo_Symbol, Chromosome, Start_Position, End_Position,
        Variant_Classification, Variant_Type, Tumor_Sample_Barcode, etc.

        Parameters
        ----------
        project_id : str
        progress_callback : callable, optional

        Returns
        -------
        pd.DataFrame
            Combined MAF table from all available files.
        """
        requests = _check_requests()

        # Check cache first
        cache_name = f"{project_id}_mutations"
        cached = self._load_cached(cache_name, fmt='csv')
        if cached is not None:
            logger.info(f"Loaded {len(cached):,} mutations from cache")
            if progress_callback:
                progress_callback(1, 1)
            return cached

        file_ids, file_meta = self._query_generic_files(
            project_id, requests,
            data_category='Simple Nucleotide Variation',
            data_type='Masked Somatic Mutation',
        )

        if not file_ids:
            raise RAPTORError(f"No open-access MAF files for {project_id}")

        logger.info(f"Downloading {len(file_ids)} MAF files for {project_id}")

        all_mafs = []
        for i, fid in enumerate(file_ids):
            content = self._download_raw_file(fid, requests)
            if content is None:
                continue
            try:
                df = pd.read_csv(io.BytesIO(content), sep='\t', comment='#',
                                 low_memory=False)
                all_mafs.append(df)
            except Exception as e:
                logger.warning(f"Failed to parse MAF file {fid}: {e}")
            if progress_callback:
                progress_callback(i + 1, len(file_ids))

        if not all_mafs:
            raise RAPTORError("No MAF data could be parsed")

        maf_df = pd.concat(all_mafs, ignore_index=True)
        logger.info(
            f"Combined MAF: {len(maf_df):,} mutations across "
            f"{maf_df['Tumor_Sample_Barcode'].nunique() if 'Tumor_Sample_Barcode' in maf_df.columns else '?'} samples"
        )
        self._cache_data(maf_df, cache_name, fmt='csv')
        return maf_df

    # =========================================================================
    # DOWNLOAD — Protein Expression (RPPA)
    # =========================================================================

    def download_rppa(
        self,
        project_id: str,
        sample_types: Optional[List[str]] = None,
        progress_callback: Optional[Callable] = None,
    ) -> AcquiredDataset:
        """
        Download Reverse Phase Protein Array (RPPA) data.

        RPPA measures ~200 proteins and phosphoproteins per sample.
        Returns a matrix of protein/phosphoprotein IDs x samples.

        Parameters
        ----------
        project_id : str
        sample_types : list of str, optional
        progress_callback : callable, optional

        Returns
        -------
        AcquiredDataset
            Protein expression matrix (~200 proteins x samples).
        """
        requests = _check_requests()

        # Check cache
        cache_name = f"{project_id}_rppa"
        if sample_types:
            cache_name += '_' + '_'.join(s.replace(' ', '') for s in sorted(sample_types))
        cached_counts = self._load_cached(cache_name, fmt='parquet')
        cached_meta = self._load_cached(f"{cache_name}_meta", fmt='csv')
        if cached_counts is not None:
            logger.info(f"Loaded RPPA from cache: {cached_counts.shape}")
            if progress_callback:
                progress_callback(1, 1)
            return AcquiredDataset(
                counts_df=cached_counts,
                metadata=cached_meta if cached_meta is not None else pd.DataFrame(),
                source_info={
                    'repository': 'TCGA', 'accession': project_id,
                    'organism': 'Homo sapiens', 'data_type': 'protein_expression_rppa',
                    'description': f"RPPA — {project_id} (cached)",
                },
                gene_id_type='protein',
            )

        file_ids, file_meta = self._query_generic_files(
            project_id, requests,
            data_category='Proteome Profiling',
            data_type='Protein Expression Quantification',
            sample_types=sample_types,
        )

        if not file_ids:
            raise RAPTORError(f"No RPPA files for {project_id}")

        logger.info(f"Downloading {len(file_ids)} RPPA files for {project_id}")

        sample_data = {}
        for i, fid in enumerate(file_ids):
            content = self._download_raw_file(fid, requests)
            if content is None:
                continue
            try:
                df = pd.read_csv(io.BytesIO(content), sep='\t')

                # RPPA files have peptide_target and protein_expression columns
                protein_col = None
                value_col = None
                for col in df.columns:
                    cl = col.lower()
                    if 'peptide_target' in cl or 'protein' in cl and 'expression' not in cl:
                        protein_col = col
                    elif 'protein_expression' in cl or 'expression' in cl:
                        value_col = col

                if protein_col is None:
                    protein_col = df.columns[0]
                if value_col is None:
                    # Try last numeric column
                    for col in reversed(df.columns):
                        if pd.to_numeric(df[col], errors='coerce').notna().sum() > len(df) * 0.5:
                            value_col = col
                            break
                if value_col is None and len(df.columns) > 1:
                    value_col = df.columns[-1]
                if value_col is None:
                    continue

                label = self._sample_label(file_meta.get(fid, {})) or fid[:12]
                sample_data[label] = pd.Series(
                    pd.to_numeric(df[value_col], errors='coerce').values,
                    index=df[protein_col].values,
                    name=label, dtype=float,
                )
            except Exception as e:
                logger.warning(f"Failed to parse RPPA file {fid}: {e}")
            if progress_callback:
                progress_callback(i + 1, len(file_ids))

        if not sample_data:
            raise RAPTORError("No RPPA data could be parsed")

        counts_df = pd.DataFrame(sample_data).fillna(0)
        counts_df.index.name = 'protein_id'

        logger.info(
            f"RPPA matrix: {counts_df.shape[0]} proteins x {counts_df.shape[1]} samples"
        )

        metadata = self._build_metadata(file_meta, keep_samples=set(counts_df.columns))
        self._cache_data(counts_df, cache_name, fmt='parquet')
        self._cache_data(metadata, f"{cache_name}_meta", fmt='csv')

        return AcquiredDataset(
            counts_df=counts_df,
            metadata=metadata,
            source_info={
                'repository': 'TCGA', 'accession': project_id,
                'organism': 'Homo sapiens', 'data_type': 'protein_expression_rppa',
                'download_date': datetime.now().isoformat(),
                'description': f"RPPA protein expression — {ALL_GDC_PROJECTS.get(project_id, project_id)}",
                'n_files_downloaded': len(sample_data),
            },
            gene_id_type='protein',
        )

    # =========================================================================
    # PAIRED TUMOR/NORMAL MATCHING
    # =========================================================================

    def get_paired_samples(
        self,
        project_id: str,
    ) -> pd.DataFrame:
        """
        Find patients that have both tumor and normal samples.

        Essential for differential expression between matched pairs.

        Parameters
        ----------
        project_id : str

        Returns
        -------
        pd.DataFrame
            Columns: patient_barcode, tumor_sample, normal_sample,
            tumor_type, normal_type. Only patients with both are returned.
        """
        requests = _check_requests()

        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_CASES_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': '=',
                        'content': {
                            'field': 'project.project_id',
                            'value': project_id,
                        }
                    }),
                    'fields': (
                        'submitter_id,'
                        'samples.sample_type,'
                        'samples.submitter_id'
                    ),
                    'size': 5000,
                },
                timeout=60,
            )
            resp.raise_for_status()
            cases = resp.json().get('data', {}).get('hits', [])
        except Exception as e:
            raise RAPTORError(f"Paired samples query failed: {e}")

        pairs = []
        for case in cases:
            patient = case.get('submitter_id', '')
            samples = case.get('samples', [])
            if not samples:
                continue

            tumors = []
            normals = []
            for s in samples:
                if not isinstance(s, dict):
                    continue
                st = s.get('sample_type', '')
                sid = s.get('submitter_id', '')
                cat = classify_sample_type(st)
                if cat == 'Tumor':
                    tumors.append((sid, st))
                elif cat == 'Normal':
                    normals.append((sid, st))

            if tumors and normals:
                # Pair first tumor with first normal
                for t_id, t_type in tumors:
                    for n_id, n_type in normals:
                        pairs.append({
                            'patient_barcode': patient,
                            'tumor_sample': t_id,
                            'tumor_type': t_type,
                            'normal_sample': n_id,
                            'normal_type': n_type,
                        })

        pairs_df = pd.DataFrame(pairs)
        logger.info(
            f"{project_id}: {len(pairs_df)} tumor-normal pairs "
            f"from {pairs_df['patient_barcode'].nunique() if not pairs_df.empty else 0} patients"
        )
        return pairs_df

    # =========================================================================
    # SURVIVAL DATA
    # =========================================================================

    def get_survival_data(
        self,
        project_id: str,
    ) -> pd.DataFrame:
        """
        Get Kaplan-Meier-ready survival data for a project.

        Returns a table with time-to-event and censoring information
        ready for lifelines, R survival, or RAPTOR M10 analysis.

        Parameters
        ----------
        project_id : str

        Returns
        -------
        pd.DataFrame
            Columns: case_id, vital_status, days_to_death,
            days_to_last_follow_up, overall_survival_days,
            overall_survival_event (1=dead, 0=censored),
            gender, age_at_diagnosis_years, tumor_stage,
            primary_diagnosis.
        """
        requests = _check_requests()

        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_CASES_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': '=',
                        'content': {
                            'field': 'project.project_id',
                            'value': project_id,
                        }
                    }),
                    'fields': (
                        'submitter_id,'
                        'demographic.vital_status,'
                        'demographic.days_to_death,'
                        'demographic.gender,'
                        'diagnoses.days_to_last_follow_up,'
                        'diagnoses.age_at_diagnosis,'
                        'diagnoses.tumor_stage,'
                        'diagnoses.primary_diagnosis,'
                        'diagnoses.days_to_recurrence'
                    ),
                    'size': 5000,
                },
                timeout=60,
            )
            resp.raise_for_status()
            cases = resp.json().get('data', {}).get('hits', [])
        except Exception as e:
            raise RAPTORError(f"Survival query failed: {e}")

        records = []
        for case in cases:
            sub_id = case.get('submitter_id', '')
            demo = case.get('demographic', {})
            if not isinstance(demo, dict):
                demo = {}
            diags = case.get('diagnoses', [])
            diag = diags[0] if diags and isinstance(diags[0], dict) else {}

            vital = demo.get('vital_status', '')
            dtd = demo.get('days_to_death')
            dtlfu = diag.get('days_to_last_follow_up')
            dtr = diag.get('days_to_recurrence')
            age = diag.get('age_at_diagnosis')

            # Compute overall survival
            if vital == 'Dead' and dtd is not None:
                os_days = dtd
                os_event = 1
            elif dtlfu is not None:
                os_days = dtlfu
                os_event = 0
            else:
                os_days = None
                os_event = None

            records.append({
                'case_id': sub_id,
                'vital_status': vital,
                'days_to_death': dtd,
                'days_to_last_follow_up': dtlfu,
                'days_to_recurrence': dtr,
                'overall_survival_days': os_days,
                'overall_survival_event': os_event,
                'gender': demo.get('gender', ''),
                'age_at_diagnosis_years': round(age / 365.25, 1) if age else None,
                'tumor_stage': diag.get('tumor_stage', ''),
                'primary_diagnosis': diag.get('primary_diagnosis', ''),
            })

        surv_df = pd.DataFrame(records)
        n_events = (surv_df['overall_survival_event'] == 1).sum() if not surv_df.empty else 0
        logger.info(
            f"{project_id}: {len(surv_df)} patients, "
            f"{n_events} events (deaths)"
        )
        return surv_df

    # =========================================================================
    # GENE MUTATION FREQUENCY (via GDC analysis endpoint)
    # =========================================================================

    def get_gene_mutation_frequency(
        self,
        project_id: str,
        n_genes: int = 50,
    ) -> pd.DataFrame:
        """
        Get the most frequently mutated genes in a project.

        Uses the GDC analysis endpoint for aggregated mutation data.

        Parameters
        ----------
        project_id : str
        n_genes : int
            Number of top mutated genes to return.

        Returns
        -------
        pd.DataFrame
            Columns: gene_id, symbol, num_cases, frequency.
        """
        requests = _check_requests()
        self._wait_for_rate_limit()

        try:
            resp = requests.get(
                GDC_TOP_MUTATED_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': '=',
                        'content': {
                            'field': 'cases.project.project_id',
                            'value': project_id,
                        }
                    }),
                    'size': n_genes,
                },
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            raise RAPTORError(f"Mutation frequency query failed: {e}")

        records = []
        for hit in data.get('data', {}).get('hits', []):
            records.append({
                'gene_id': hit.get('gene_id', ''),
                'symbol': hit.get('symbol', ''),
                'name': hit.get('name', ''),
                'num_cases': hit.get('_score', 0),
                'cytoband': hit.get('cytoband', []),
                'biotype': hit.get('biotype', ''),
            })

        freq_df = pd.DataFrame(records)
        if not freq_df.empty and 'num_cases' in freq_df.columns:
            freq_df = freq_df.sort_values('num_cases', ascending=False)
        return freq_df

    # =========================================================================
    # MANIFEST GENERATION (for GDC Data Transfer Tool)
    # =========================================================================

    def generate_manifest(
        self,
        project_id: str,
        output_path: Optional[Union[str, Path]] = None,
        data_category: str = 'Transcriptome Profiling',
        data_type: str = 'Gene Expression Quantification',
        sample_types: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Generate a GDC download manifest for use with gdc-client.

        The manifest TSV file can be fed to the GDC Data Transfer Tool:
            gdc-client download -m manifest.tsv

        Parameters
        ----------
        project_id : str
        output_path : str or Path, optional
            If given, write the manifest to this TSV file.
        data_category : str
        data_type : str
        sample_types : list of str, optional

        Returns
        -------
        pd.DataFrame
            Manifest table with columns: id, filename, md5, size, state.
        """
        requests = _check_requests()

        filter_content = [
            {'op': '=', 'content': {'field': 'cases.project.project_id', 'value': project_id}},
            {'op': '=', 'content': {'field': 'data_category', 'value': data_category}},
            {'op': '=', 'content': {'field': 'data_type', 'value': data_type}},
            {'op': '=', 'content': {'field': 'access', 'value': 'open'}},
        ]
        if sample_types:
            filter_content.append({
                'op': 'in',
                'content': {'field': 'cases.samples.sample_type', 'value': sample_types}
            })

        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_FILES_ENDPOINT,
                params={
                    'filters': json.dumps({'op': 'and', 'content': filter_content}),
                    'fields': 'file_id,file_name,md5sum,file_size,state',
                    'size': 10000,
                    'return_type': 'manifest',
                },
                timeout=60,
            )
            resp.raise_for_status()

            manifest_df = pd.read_csv(io.StringIO(resp.text), sep='\t')

            if output_path:
                output_path = Path(output_path)
                manifest_df.to_csv(output_path, sep='\t', index=False)
                logger.info(f"Manifest written to {output_path} ({len(manifest_df)} files)")

            return manifest_df

        except Exception as e:
            raise RAPTORError(f"Manifest generation failed: {e}")

    # =========================================================================
    # PROJECT BROWSING
    # =========================================================================

    def list_projects(self, program: Optional[str] = None) -> Dict[str, str]:
        """List known project IDs and descriptions (static)."""
        if program == 'TCGA':
            return dict(TCGA_PROJECTS)
        elif program == 'TARGET':
            return dict(TARGET_PROJECTS)
        return dict(ALL_GDC_PROJECTS)

    def list_projects_live(self) -> List[Dict[str, Any]]:
        """Query GDC for all projects with current case counts."""
        requests = _check_requests()
        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_PROJECTS_ENDPOINT,
                params={
                    'fields': (
                        'project_id,name,primary_site,disease_type,'
                        'program.name,summary.case_count,'
                        'summary.data_categories.data_category,'
                        'summary.data_categories.file_count'
                    ),
                    'size': 500,
                },
                timeout=30,
            )
            resp.raise_for_status()
        except Exception as e:
            raise RAPTORError(f"GDC project listing failed: {e}")

        projects = []
        for hit in resp.json().get('data', {}).get('hits', []):
            data_cats = {
                d['data_category']: d['file_count']
                for d in hit.get('summary', {}).get('data_categories', [])
            }
            projects.append({
                'project_id': hit.get('project_id', ''),
                'name': hit.get('name', ''),
                'primary_site': hit.get('primary_site', []),
                'disease_type': hit.get('disease_type', []),
                'case_count': hit.get('summary', {}).get('case_count', 0),
                'program': hit.get('program', {}).get('name', ''),
                'data_categories': data_cats,
            })

        projects.sort(key=lambda p: (p.get('program', ''), p.get('project_id', '')))
        return projects

    def compare_projects(
        self,
        project_ids: List[str],
    ) -> pd.DataFrame:
        """
        Compare multiple projects side by side.

        Returns a DataFrame with case counts, available data types,
        and clinical summaries for each project.
        """
        rows = []
        for pid in project_ids:
            try:
                info = self.get_study_info(pid)
                cs = info.get('clinical_summary', {})
                stc = info.get('sample_type_counts', {})

                tumor_count = sum(c for st, c in stc.items()
                                  if classify_sample_type(st) == 'Tumor')
                normal_count = sum(c for st, c in stc.items()
                                   if classify_sample_type(st) == 'Normal')

                rows.append({
                    'project_id': pid,
                    'name': info.get('title', ''),
                    'cases': info.get('n_samples', 0),
                    'tumor_samples': tumor_count,
                    'normal_samples': normal_count,
                    'has_rnaseq': 'RNA-Seq' in info.get('experimental_strategies', {}),
                    'has_methylation': 'DNA Methylation' in info.get('data_categories', {}),
                    'has_cnv': 'Copy Number Variation' in info.get('data_categories', {}),
                    'has_mutations': 'Simple Nucleotide Variation' in info.get('data_categories', {}),
                    'median_age': cs.get('age_at_diagnosis', {}).get('median', ''),
                    'pct_male': round(
                        100 * cs.get('gender', {}).get('male', 0) /
                        max(sum(cs.get('gender', {}).values()), 1), 1
                    ) if cs.get('gender') else '',
                })
            except Exception as e:
                logger.warning(f"Failed to get info for {pid}: {e}")
                rows.append({'project_id': pid, 'name': 'Error', 'cases': 0})

        return pd.DataFrame(rows)

    # =========================================================================
    # UTILITY — API status, quant columns, etc.
    # =========================================================================

    @staticmethod
    def get_quant_columns() -> Dict[str, str]:
        """Available quantification columns from STAR-Counts files."""
        return dict(QUANT_COLUMNS)

    @staticmethod
    def get_sample_type_groups() -> Dict[str, List[str]]:
        """Predefined sample type groups for filtering."""
        return dict(SAMPLE_TYPE_GROUPS)

    @staticmethod
    def get_available_data_types() -> Dict[str, Dict[str, str]]:
        """All downloadable data types from GDC."""
        return dict(GDC_DATA_TYPES)

    def check_status(self) -> Dict[str, Any]:
        """Check GDC API health and data release version."""
        requests = _check_requests()
        try:
            resp = requests.get(GDC_STATUS_ENDPOINT, timeout=10)
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            return {'status': 'error', 'message': str(e)}

    # =========================================================================
    # INFO (lightweight, for BaseConnector interface)
    # =========================================================================

    def _info_api(self, accession: str) -> Optional[Dict[str, Any]]:
        """Lightweight project info for the catalog."""
        requests = _check_requests()
        try:
            resp = requests.get(
                f"{GDC_PROJECTS_ENDPOINT}/{accession}",
                params={'fields': 'name,primary_site,disease_type,summary.case_count'},
                timeout=15,
            )
            if resp.status_code == 200:
                data = resp.json().get('data', {})
                ps = data.get('primary_site', [])
                dt = data.get('disease_type', [])
                return {
                    'accession': accession,
                    'title': data.get('name', ''),
                    'organism': 'Homo sapiens',
                    'n_samples': data.get('summary', {}).get('case_count', 0),
                    'primary_site': ps if isinstance(ps, list) else [ps],
                    'disease_type': dt if isinstance(dt, list) else [dt],
                    'repository': 'TCGA',
                    'platform_id': 'Illumina HiSeq/NovaSeq',
                    'data_type': 'RNA-seq',
                }
        except Exception as e:
            logger.debug(f"GDC info failed for {accession}: {e}")
        return None

    # =========================================================================
    # INTERNAL — File queries and downloads
    # =========================================================================

    def _query_files(
        self,
        project_id: str,
        requests,
        sample_types: Optional[List[str]] = None,
    ) -> Tuple[List[str], Dict]:
        """Query GDC for RNA-seq quantification files with pagination."""
        return self._query_generic_files(
            project_id, requests,
            data_category='Transcriptome Profiling',
            data_type='Gene Expression Quantification',
            workflow_type=self._workflow_type,
            experimental_strategy='RNA-Seq',
            sample_types=sample_types,
        )

    def _query_generic_files(
        self,
        project_id: str,
        requests,
        data_category: str,
        data_type: str,
        workflow_type: Optional[str] = None,
        experimental_strategy: Optional[str] = None,
        sample_types: Optional[List[str]] = None,
    ) -> Tuple[List[str], Dict]:
        """
        Generic file query with pagination. Used by all download methods.
        """
        filter_content = [
            {'op': '=', 'content': {'field': 'cases.project.project_id', 'value': project_id}},
            {'op': '=', 'content': {'field': 'data_category', 'value': data_category}},
            {'op': '=', 'content': {'field': 'data_type', 'value': data_type}},
            {'op': '=', 'content': {'field': 'access', 'value': 'open'}},
        ]
        if workflow_type:
            filter_content.append(
                {'op': '=', 'content': {'field': 'analysis.workflow_type', 'value': workflow_type}}
            )
        if experimental_strategy:
            filter_content.append(
                {'op': '=', 'content': {'field': 'experimental_strategy', 'value': experimental_strategy}}
            )
        if sample_types:
            filter_content.append(
                {'op': 'in', 'content': {'field': 'cases.samples.sample_type', 'value': sample_types}}
            )

        params = {
            'filters': json.dumps({'op': 'and', 'content': filter_content}),
            'fields': (
                'file_id,file_name,file_size,'
                'cases.case_id,cases.submitter_id,'
                'cases.samples.sample_type,cases.samples.tissue_type,'
                'cases.samples.submitter_id,'
                'cases.demographic.gender,cases.demographic.vital_status,'
                'cases.demographic.days_to_death,cases.demographic.year_of_birth,'
                'cases.diagnoses.primary_diagnosis,cases.diagnoses.tumor_stage,'
                'cases.diagnoses.age_at_diagnosis,cases.diagnoses.morphology'
            ),
            'size': 5000,
        }

        all_hits = []
        page_from = 0

        while True:
            params['from'] = page_from
            self._wait_for_rate_limit()
            try:
                resp = requests.get(GDC_FILES_ENDPOINT, params=params, timeout=60)
                resp.raise_for_status()
                data = resp.json()
            except Exception as e:
                raise RAPTORError(f"GDC file query failed: {e}")

            hits = data.get('data', {}).get('hits', [])
            all_hits.extend(hits)

            total = data.get('data', {}).get('pagination', {}).get('total', 0)
            page_from += len(hits)
            if page_from >= total or not hits:
                break
            logger.info(f"  Paginating... {page_from}/{total} files queried")

        file_ids = []
        file_metadata = {}

        for hit in all_hits:
            fid = hit.get('file_id')
            if not fid:
                continue
            file_ids.append(fid)

            cases = hit.get('cases', [{}])
            case = cases[0] if cases else {}
            samples = case.get('samples', [{}])
            sample = samples[0] if samples and isinstance(samples[0], dict) else {}
            demo = case.get('demographic', {})
            diags = case.get('diagnoses', [{}])
            diag = diags[0] if diags and isinstance(diags[0], dict) else {}
            if not isinstance(demo, dict):
                demo = {}

            file_metadata[fid] = {
                'case_id': case.get('case_id', ''),
                'submitter_id': case.get('submitter_id', fid[:12]),
                'sample_type': sample.get('sample_type', 'unknown'),
                'tissue_type': sample.get('tissue_type', 'unknown'),
                'sample_submitter_id': sample.get('submitter_id', ''),
                'gender': demo.get('gender', 'unknown'),
                'vital_status': demo.get('vital_status', 'unknown'),
                'days_to_death': demo.get('days_to_death'),
                'year_of_birth': demo.get('year_of_birth'),
                'primary_diagnosis': diag.get('primary_diagnosis', 'unknown'),
                'tumor_stage': diag.get('tumor_stage', 'unknown'),
                'age_at_diagnosis': diag.get('age_at_diagnosis'),
                'morphology': diag.get('morphology', ''),
                'file_name': hit.get('file_name', ''),
                'file_size': hit.get('file_size', 0),
            }

        return file_ids, file_metadata

    def _cache_data(self, df: pd.DataFrame, cache_name: str, fmt: str = 'parquet'):
        """Cache any downloaded data to disk for fast reload."""
        try:
            cache_dir = Path.home() / '.raptor' / 'parquet_cache'
            cache_dir.mkdir(parents=True, exist_ok=True)
            if fmt == 'parquet':
                df.to_parquet(cache_dir / f"{cache_name}.parquet",
                              engine='pyarrow', compression='snappy')
            else:
                df.to_csv(cache_dir / f"{cache_name}.csv.gz",
                          compression='gzip', index=True)
            logger.info(f"Cached: {cache_name}")
        except Exception as e:
            logger.debug(f"Cache save failed (non-critical): {e}")

    def _load_cached(self, cache_name: str, fmt: str = 'parquet') -> Optional[pd.DataFrame]:
        """Load cached data if available."""
        cache_dir = Path.home() / '.raptor' / 'parquet_cache'
        try:
            if fmt == 'parquet':
                path = cache_dir / f"{cache_name}.parquet"
                if path.exists():
                    df = pd.read_parquet(path, engine='pyarrow')
                    logger.info(f"Loaded from cache: {cache_name}")
                    return df
            else:
                path = cache_dir / f"{cache_name}.csv.gz"
                if path.exists():
                    df = pd.read_csv(path, compression='gzip', index_col=0)
                    logger.info(f"Loaded from cache: {cache_name}")
                    return df
        except Exception as e:
            logger.debug(f"Cache load failed: {e}")
        return None

    def _download_raw_file(self, file_id: str, requests, max_retries: int = 3) -> Optional[bytes]:
        """Download a single file with retry logic and exponential backoff."""
        import time as _time

        for attempt in range(max_retries):
            self._wait_for_rate_limit()
            try:
                resp = requests.get(
                    f"{GDC_DATA_ENDPOINT}/{file_id}",
                    timeout=120 + (attempt * 60),  # increase timeout each retry
                )
                resp.raise_for_status()

                content = resp.content
                try:
                    content = gzip.decompress(content)
                except (gzip.BadGzipFile, OSError):
                    pass
                return content

            except Exception as e:
                if attempt < max_retries - 1:
                    wait = 2 ** (attempt + 1)  # 2s, 4s, 8s
                    logger.warning(
                        f"Download {file_id} failed (attempt {attempt + 1}/{max_retries}): {e}. "
                        f"Retrying in {wait}s..."
                    )
                    _time.sleep(wait)
                else:
                    logger.warning(f"Download {file_id} failed after {max_retries} attempts: {e}")
                    return None

    # =========================================================================
    # INTERNAL — Expression count downloads (batch + individual)
    # =========================================================================

    def _download_counts_batch(
        self, file_ids, file_metadata, requests,
        quant_col='unstranded', batch_size=50,
        progress_callback=None,
    ) -> pd.DataFrame:
        """
        Batch download via POST to /data with retry, adaptive batch size,
        and resume on partial failure.
        """
        import time as _time

        sample_counts = {}
        total = len(file_ids)
        downloaded = 0
        failed_ids = []
        current_batch_size = min(batch_size, max(10, total // 10))  # adaptive start

        # For very large projects, use smaller batches
        if total > 500:
            current_batch_size = 30
        elif total > 200:
            current_batch_size = 40

        batch_start = 0
        max_retries = 3

        while batch_start < total:
            # Check for cancel flag
            cancel_path = Path.home() / '.raptor' / '_cancel_download'
            if cancel_path.exists():
                cancel_path.unlink(missing_ok=True)
                logger.info(f"Download cancelled by user at {downloaded}/{total} files")
                break
            
            batch_ids = file_ids[batch_start:batch_start + current_batch_size]
            batch_end = min(batch_start + current_batch_size, total)

            logger.info(
                f"  Batch download: files {batch_start + 1}-{batch_end} "
                f"of {total} (batch_size={current_batch_size})..."
            )

            success = False
            for attempt in range(max_retries):
                self._wait_for_rate_limit()
                try:
                    resp = requests.post(
                        GDC_DATA_ENDPOINT,
                        json={"ids": batch_ids},
                        headers={"Content-Type": "application/json"},
                        timeout=300 + (attempt * 120),  # 300s, 420s, 540s
                    )
                    resp.raise_for_status()
                    success = True
                    break
                except Exception as e:
                    if attempt < max_retries - 1:
                        wait = 3 * (attempt + 1)
                        logger.warning(
                            f"Batch {batch_start + 1}-{batch_end} failed "
                            f"(attempt {attempt + 1}/{max_retries}): {e}. "
                            f"Retrying in {wait}s..."
                        )
                        _time.sleep(wait)

                        # Reduce batch size on retry
                        if current_batch_size > 10:
                            current_batch_size = max(10, current_batch_size // 2)
                            batch_ids = file_ids[batch_start:batch_start + current_batch_size]
                            batch_end = min(batch_start + current_batch_size, total)
                            logger.info(f"  Reduced batch size to {current_batch_size}")
                    else:
                        logger.warning(
                            f"Batch {batch_start + 1}-{batch_end} failed after "
                            f"{max_retries} attempts. Falling back to individual downloads."
                        )

            if success:
                # Parse tar.gz response
                parsed_in_batch = set()
                try:
                    with tarfile.open(fileobj=io.BytesIO(resp.content), mode='r:gz') as tar:
                        for member in tar.getmembers():
                            if not member.isfile():
                                continue
                            path_parts = member.name.split('/')
                            file_uuid = path_parts[0] if len(path_parts) >= 2 else None
                            if file_uuid not in file_metadata:
                                fname = path_parts[-1]
                                file_uuid = next(
                                    (fid for fid, m in file_metadata.items()
                                     if m.get('file_name') == fname), None
                                )
                                if file_uuid is None:
                                    continue

                            f = tar.extractfile(member)
                            if f is None:
                                continue
                            content = f.read()
                            try:
                                content = gzip.decompress(content)
                            except (gzip.BadGzipFile, OSError):
                                pass

                            series = self._parse_count_file(content, file_uuid, file_metadata, quant_col)
                            if series is not None:
                                sample_counts[series.name] = series
                                downloaded += 1
                                parsed_in_batch.add(file_uuid)

                except (tarfile.TarError, EOFError, OSError, gzip.BadGzipFile) as tar_err:
                    # Truncated response OR single-file response OR corrupt archive
                    if len(batch_ids) == 1 and not parsed_in_batch:
                        content = resp.content
                        try:
                            content = gzip.decompress(content)
                        except (gzip.BadGzipFile, OSError):
                            pass
                        series = self._parse_count_file(content, batch_ids[0], file_metadata, quant_col)
                        if series is not None:
                            sample_counts[series.name] = series
                            downloaded += 1
                    else:
                        # Fall back to individual downloads for unparsed files
                        remaining = [fid for fid in batch_ids if fid not in parsed_in_batch]
                        logger.warning(
                            f"Tar parse failed ({type(tar_err).__name__}): "
                            f"parsed {len(parsed_in_batch)}/{len(batch_ids)} files. "
                            f"Falling back to individual downloads for {len(remaining)} remaining."
                        )
                        for fid in remaining:
                            series = self._download_and_parse_single(
                                fid, file_metadata, requests, quant_col
                            )
                            if series is not None:
                                sample_counts[series.name] = series
                                downloaded += 1
            else:
                # All batch retries failed — download individually
                for fid in batch_ids:
                    series = self._download_and_parse_single(
                        fid, file_metadata, requests, quant_col
                    )
                    if series is not None:
                        sample_counts[series.name] = series
                        downloaded += 1
                    else:
                        failed_ids.append(fid)

                    if progress_callback:
                        progress_callback(downloaded, total)

            if progress_callback:
                progress_callback(downloaded, total)

            batch_start += len(batch_ids)

        # Retry failed files one more time
        if failed_ids:
            logger.info(f"  Retrying {len(failed_ids)} failed files individually...")
            for fid in failed_ids:
                series = self._download_and_parse_single(
                    fid, file_metadata, requests, quant_col
                )
                if series is not None:
                    sample_counts[series.name] = series
                    downloaded += 1

        if not sample_counts:
            raise RAPTORError("No count data could be extracted from downloads")

        counts_df = pd.DataFrame(sample_counts).fillna(0)
        counts_df.index.name = 'gene_id'

        n_failed = total - downloaded
        if n_failed > 0:
            logger.warning(f"  {n_failed}/{total} files could not be downloaded")
        logger.info(f"Built count matrix: {counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples")
        return counts_df

    def _download_counts_individual(
        self, file_ids, file_metadata, requests,
        quant_col='unstranded', progress_callback=None,
    ) -> pd.DataFrame:
        """Individual file downloads (for small projects or fallback)."""
        sample_counts = {}
        for i, fid in enumerate(file_ids):
            series = self._download_and_parse_single(fid, file_metadata, requests, quant_col)
            if series is not None:
                sample_counts[series.name] = series
            if progress_callback:
                progress_callback(i + 1, len(file_ids))
            if (i + 1) % 50 == 0:
                logger.info(f"  Downloaded {i + 1}/{len(file_ids)} files...")

        if not sample_counts:
            raise RAPTORError("No count data could be extracted from downloads")

        counts_df = pd.DataFrame(sample_counts).fillna(0)
        counts_df.index.name = 'gene_id'
        logger.info(f"Built count matrix: {counts_df.shape[0]:,} genes x {counts_df.shape[1]} samples")
        return counts_df

    def _download_and_parse_single(self, file_id, file_metadata, requests, quant_col):
        """Download + parse a single count file."""
        content = self._download_raw_file(file_id, requests)
        if content is None:
            return None
        return self._parse_count_file(content, file_id, file_metadata, quant_col)

    def _parse_count_file(self, content, file_id, file_metadata, quant_col):
        """Parse a STAR-Counts TSV → pd.Series with Ensembl gene IDs."""
        try:
            df = pd.read_csv(io.BytesIO(content), sep='\t', comment='#', header=0)
        except Exception as e:
            logger.warning(f"Failed to parse file {file_id}: {e}")
            return None

        gene_col = next((c for c in df.columns if 'gene_id' in c.lower()), df.columns[0])

        count_col = None
        if quant_col in df.columns:
            count_col = quant_col
        else:
            for col in df.columns:
                if quant_col.lower() in col.lower():
                    count_col = col
                    break
        if count_col is None:
            for candidate in ['unstranded', 'tpm_unstranded', 'fpkm_unstranded', 'expected_count']:
                if candidate in df.columns:
                    count_col = candidate
                    break
        if count_col is None and len(df.columns) >= 2:
            count_col = df.columns[1]
        if count_col is None:
            return None

        genes = df[gene_col].astype(str).str.split('.').str[0]
        gene_mask = genes.str.startswith('ENSG')

        meta = file_metadata.get(file_id, {})
        sample_sub = meta.get('sample_submitter_id', '')
        label = sample_sub if sample_sub and sample_sub != meta.get('submitter_id', '') else meta.get('submitter_id', file_id[:12])

        return pd.Series(
            df.loc[gene_mask, count_col].values,
            index=genes[gene_mask].values,
            name=label, dtype=float,
        )

    # =========================================================================
    # INTERNAL — Metadata building
    # =========================================================================

    @staticmethod
    def _sample_label(meta: dict) -> str:
        """
        Consistent sample label from file metadata.

        Uses sample_submitter_id (e.g. TCGA-W5-AA2R-01A) when available,
        falls back to submitter_id (case-level, e.g. TCGA-W5-AA2R).
        Must match the label used in _build_metadata so filtering works.
        """
        sub_id = meta.get('submitter_id', '')
        sample_sub = meta.get('sample_submitter_id', '')
        return sample_sub if sample_sub and sample_sub != sub_id else (sub_id or 'unknown')

    def _build_metadata(self, file_metadata: dict, keep_samples: set = None) -> pd.DataFrame:
        """Build sample metadata from file-level metadata.
        
        Parameters
        ----------
        file_metadata : dict
            File UUID → metadata dict from GDC query.
        keep_samples : set, optional
            If provided, only include samples whose label is in this set.
            Use this to align metadata with successfully downloaded data.
        """
        records = []
        seen = set()
        for fid, meta in file_metadata.items():
            label = self._sample_label(meta)
            if not label or label == 'unknown':
                label = fid[:12]
            if label in seen:
                continue
            if keep_samples is not None and label not in keep_samples:
                continue
            seen.add(label)

            sub_id = meta.get('submitter_id', fid[:12])
            st = meta.get('sample_type', 'unknown')
            barcode_info = parse_tcga_barcode(label) if label.startswith('TCGA-') else {}
            age = meta.get('age_at_diagnosis')

            records.append({
                'sample_id': label,
                'case_id': sub_id,
                'sample_type': st,
                'sample_category': classify_sample_type(st),
                'tissue_type': meta.get('tissue_type', 'unknown'),
                'gender': meta.get('gender', 'unknown'),
                'vital_status': meta.get('vital_status', 'unknown'),
                'days_to_death': meta.get('days_to_death'),
                'primary_diagnosis': meta.get('primary_diagnosis', 'unknown'),
                'tumor_stage': meta.get('tumor_stage', 'unknown'),
                'age_at_diagnosis_years': round(age / 365.25, 1) if age else None,
                'patient_barcode': barcode_info.get('patient_barcode', ''),
            })

        metadata = pd.DataFrame(records)
        if 'sample_id' in metadata.columns and not metadata.empty:
            metadata = metadata.set_index('sample_id')
        metadata.index.name = 'sample_id'
        return metadata

    def _fetch_clinical_metadata(self, project_id, file_metadata, requests):
        """Build metadata from file_metadata + enrich from /cases."""
        metadata = self._build_metadata(file_metadata)

        try:
            enriched = self._enrich_clinical(project_id, requests)
            if enriched is not None and not enriched.empty:
                metadata = metadata.reset_index()
                metadata = metadata.merge(
                    enriched, left_on='case_id', right_index=True,
                    how='left', suffixes=('', '_clinical'),
                )
                for col in [c for c in metadata.columns if c.endswith('_clinical')]:
                    base = col.replace('_clinical', '')
                    if base in metadata.columns:
                        mask = metadata[base].isna() | (metadata[base] == 'unknown')
                        metadata.loc[mask, base] = metadata.loc[mask, col]
                    metadata = metadata.drop(columns=[col])
                metadata = metadata.set_index('sample_id')
        except Exception as e:
            logger.debug(f"Clinical enrichment skipped: {e}")

        return metadata

    def _enrich_clinical(self, project_id, requests):
        """Fetch survival, staging, exposure data from /cases."""
        self._wait_for_rate_limit()
        try:
            resp = requests.get(
                GDC_CASES_ENDPOINT,
                params={
                    'filters': json.dumps({
                        'op': '=',
                        'content': {
                            'field': 'project.project_id',
                            'value': project_id,
                        }
                    }),
                    'fields': (
                        'submitter_id,'
                        'demographic.days_to_death,demographic.year_of_death,'
                        'diagnoses.days_to_last_follow_up,'
                        'diagnoses.days_to_recurrence,'
                        'diagnoses.site_of_resection_or_biopsy,'
                        'diagnoses.prior_malignancy,'
                        'diagnoses.ajcc_pathologic_stage,'
                        'diagnoses.ajcc_pathologic_t,'
                        'diagnoses.ajcc_pathologic_n,'
                        'diagnoses.ajcc_pathologic_m,'
                        'exposures.alcohol_history,'
                        'exposures.tobacco_smoking_status,'
                        'exposures.pack_years_smoked'
                    ),
                    'size': 5000,
                },
                timeout=60,
            )
            resp.raise_for_status()
            cases = resp.json().get('data', {}).get('hits', [])
        except Exception:
            return None

        records = []
        for case in cases:
            demo = case.get('demographic', {})
            if not isinstance(demo, dict):
                demo = {}
            diags = case.get('diagnoses', [])
            diag = diags[0] if diags and isinstance(diags[0], dict) else {}
            exps = case.get('exposures', [])
            exp = exps[0] if exps and isinstance(exps[0], dict) else {}

            records.append({
                'case_submitter_id': case.get('submitter_id', ''),
                'days_to_last_follow_up': diag.get('days_to_last_follow_up'),
                'days_to_recurrence': diag.get('days_to_recurrence'),
                'site_of_resection': diag.get('site_of_resection_or_biopsy', ''),
                'prior_malignancy': diag.get('prior_malignancy', ''),
                'ajcc_stage': diag.get('ajcc_pathologic_stage', ''),
                'ajcc_t': diag.get('ajcc_pathologic_t', ''),
                'ajcc_n': diag.get('ajcc_pathologic_n', ''),
                'ajcc_m': diag.get('ajcc_pathologic_m', ''),
                'alcohol_history': exp.get('alcohol_history', ''),
                'tobacco_smoking_status': exp.get('tobacco_smoking_status', ''),
                'pack_years_smoked': exp.get('pack_years_smoked'),
            })

        if records:
            df = pd.DataFrame(records).set_index('case_submitter_id')
            df.index.name = 'case_id'
            return df
        return None

    # =========================================================================
    # DISPLAY
    # =========================================================================

    def __repr__(self) -> str:
        cached = len(self.list_cached())
        return (
            f"TCGAConnector("
            f"workflow='{self._workflow_type}', "
            f"quant='{self._quant_column}', "
            f"cache={'on' if self._cache_manager.enabled else 'off'}, "
            f"cached={cached})"
        )
