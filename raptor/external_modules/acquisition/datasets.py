"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - Dataset Containers

Data container classes for acquired and pooled datasets from public repositories.
Follows the DEResult pattern from Module 7 for consistency across RAPTOR.

Workflow Position:
    M6 Acquisition: Repository API → AcquiredDataset
    M6 Pooling: Multiple AcquiredDatasets → PooledDataset
    M10: AcquiredDataset/PooledDataset → Biomarker Discovery

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import json
import pickle
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# RAPTOR imports
try:
    from raptor.utils.errors import ValidationError, RAPTORError
except ImportError:
    class ValidationError(Exception):
        """Validation error fallback."""
        def __init__(self, parameter: str, message: str = "", hint: str = ""):
            self.parameter = parameter
            self.message = message
            self.hint = hint
            super().__init__(f"{parameter}: {message}. {hint}")

    class RAPTORError(Exception):
        """Base RAPTOR error fallback."""
        pass


logger = logging.getLogger(__name__)


# =============================================================================
# Constants
# =============================================================================

SUPPORTED_REPOSITORIES = ['GEO', 'TCGA', 'ArrayExpress', 'SRA']

SUPPORTED_DATA_TYPES = ['raw_counts', 'normalized_counts', 'deg_table']

SUPPORTED_ORGANISMS = [
    'Homo sapiens', 'Mus musculus', 'Rattus norvegicus',
    'Drosophila melanogaster', 'Danio rerio', 'Caenorhabditis elegans',
]

REQUIRED_COUNTS_COLUMNS = []  # gene_id is the index; samples are columns

REQUIRED_DEG_COLUMNS = ['log2_fold_change', 'p_value', 'adjusted_p_value']


# =============================================================================
# AcquiredDataset
# =============================================================================

@dataclass
class AcquiredDataset:
    """
    A dataset acquired from a public repository.

    Core data container for Module 6 data acquisition. Represents a single
    dataset downloaded from GEO, TCGA, ArrayExpress, or SRA.

    Attributes
    ----------
    counts_df : pd.DataFrame
        Count matrix with gene IDs as index and sample IDs as columns.
        Can be raw counts, normalized counts, or a DEG table depending
        on data_type.

    metadata : pd.DataFrame
        Sample-level metadata with sample IDs as index and annotation
        columns (condition, tissue, batch, etc.).

    source_info : dict
        Provenance information including:
        - repository : str (GEO, TCGA, ArrayExpress, SRA)
        - accession : str (e.g., GSE12345, TCGA-LIHC)
        - organism : str
        - platform : str (sequencing platform)
        - download_date : str (ISO format)
        - data_type : str (raw_counts, normalized_counts, deg_table)
        - description : str

    gene_id_type : str
        Type of gene identifiers used: 'ensembl', 'symbol', 'entrez'

    Examples
    --------
    >>> from raptor.external_modules.acquisition import GEOConnector
    >>> geo = GEOConnector()
    >>> results = geo.search("pancreatic cancer RNA-seq")
    >>> dataset = geo.download("GSE12345")
    >>> print(f"Samples: {dataset.n_samples}, Genes: {dataset.n_genes}")
    >>> dataset.save("my_dataset.pkl")
    """

    counts_df: pd.DataFrame
    metadata: pd.DataFrame
    source_info: Dict[str, Any] = field(default_factory=dict)
    gene_id_type: str = 'symbol'

    def __post_init__(self):
        """Validate AcquiredDataset after initialization."""
        # Validate counts_df
        if not isinstance(self.counts_df, pd.DataFrame):
            raise ValidationError('counts_df', "counts_df must be a pandas DataFrame", "Pass a DataFrame with gene IDs as index and samples as columns"
            )

        if self.counts_df.empty:
            raise ValidationError('counts_df', "counts_df is empty", "DataFrame must contain at least one gene and one sample"
            )

        # Ensure index has a name
        if self.counts_df.index.name is None:
            self.counts_df.index.name = 'gene_id'

        # Validate metadata
        if not isinstance(self.metadata, pd.DataFrame):
            raise ValidationError('metadata', "metadata must be a pandas DataFrame", "Pass a DataFrame with sample IDs as index"
            )

        if self.metadata.index.name is None:
            self.metadata.index.name = 'sample_id'

        # Validate gene_id_type
        valid_id_types = [
            'ensembl', 'symbol', 'entrez', 'refseq', 'unknown',
            # Multi-omic types (TCGA miRNA, methylation, CNV, protein)
            'mirna', 'probe', 'segment', 'genomic_region', 'protein',
        ]
        if self.gene_id_type not in valid_id_types:
            raise ValidationError('gene_id_type', f"Invalid gene_id_type: {self.gene_id_type}", f"Must be one of: {valid_id_types}"
            )

        # Set default source_info fields
        defaults = {
            'repository': 'unknown',
            'accession': 'unknown',
            'organism': 'unknown',
            'platform': 'unknown',
            'download_date': datetime.now().isoformat(),
            'data_type': 'raw_counts',
            'description': '',
        }
        for key, value in defaults.items():
            self.source_info.setdefault(key, value)

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

    @property
    def n_genes(self) -> int:
        """Number of genes (rows) in the count matrix."""
        return len(self.counts_df)

    @property
    def n_samples(self) -> int:
        """Number of samples (columns) in the count matrix."""
        return len(self.counts_df.columns)

    @property
    def accession(self) -> str:
        """Repository accession ID."""
        return self.source_info.get('accession', 'unknown')

    @property
    def repository(self) -> str:
        """Source repository name."""
        return self.source_info.get('repository', 'unknown')

    @property
    def organism(self) -> str:
        """Organism of the dataset."""
        return self.source_info.get('organism', 'unknown')

    @property
    def data_type(self) -> str:
        """Type of data: raw_counts, normalized_counts, or deg_table."""
        return self.source_info.get('data_type', 'raw_counts')

    @property
    def gene_ids(self) -> List[str]:
        """List of gene IDs."""
        return self.counts_df.index.tolist()

    @property
    def sample_ids(self) -> List[str]:
        """List of sample IDs."""
        return self.counts_df.columns.tolist()

    @property
    def samples_in_metadata(self) -> int:
        """Number of samples that have metadata."""
        common = set(self.counts_df.columns) & set(self.metadata.index)
        return len(common)

    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------

    def validate_integrity(self) -> Dict[str, Any]:
        """
        Check dataset integrity and return a report.

        Adapts terminology and thresholds based on data type:
        methylation (probes, beta values), miRNA, CNV, protein, etc.

        Returns
        -------
        dict
            Integrity report with keys: valid, warnings, errors
        """
        report = {'valid': True, 'warnings': [], 'errors': []}

        # Data-type-aware labels
        _labels = {
            'mirna': ('miRNAs', 'expression matrix'),
            'probe': ('probes', 'beta value matrix'),
            'segment': ('segments', 'segment matrix'),
            'genomic_region': ('regions', 'CNV matrix'),
            'protein': ('proteins', 'protein expression matrix'),
        }
        feat_name, mat_name = _labels.get(
            self.gene_id_type, ('genes', 'count matrix')
        )

        # Check sample overlap between data and metadata
        count_samples = set(self.counts_df.columns)
        meta_samples = set(self.metadata.index)

        if not count_samples & meta_samples:
            report['valid'] = False
            report['errors'].append(
                f"No overlap between {mat_name} columns and metadata index"
            )
        elif count_samples != meta_samples:
            missing_meta = count_samples - meta_samples
            missing_counts = meta_samples - count_samples
            if missing_meta:
                report['warnings'].append(
                    f"{len(missing_meta)} samples in {mat_name} lack metadata"
                )
            if missing_counts:
                report['warnings'].append(
                    f"{len(missing_counts)} samples in metadata lack data"
                )

        # Check for all-zero features
        zero_feats = (self.counts_df.sum(axis=1) == 0).sum()
        if zero_feats > 0:
            pct = 100 * zero_feats / self.n_genes
            if self.gene_id_type == 'probe':
                # Methylation: zero beta = fully unmethylated, only warn if many
                if pct > 20:
                    report['warnings'].append(
                        f"{zero_feats} probes ({pct:.1f}%) have beta = 0 across all samples "
                        f"(may indicate failed probes or fully unmethylated CpGs)"
                    )
            else:
                report['warnings'].append(
                    f"{zero_feats} {feat_name} ({pct:.1f}%) have zero values across all samples"
                )

        # Check for NaN values
        nan_count = self.counts_df.isna().sum().sum()
        if nan_count > 0:
            if self.gene_id_type == 'probe':
                # Methylation NaN = failed probe detection, expected behavior
                total_cells = self.counts_df.shape[0] * self.counts_df.shape[1]
                nan_pct = 100 * nan_count / total_cells
                report['warnings'].append(
                    f"{nan_count:,} NaN values ({nan_pct:.1f}%) in {mat_name} "
                    f"(probes that failed detection — expected for methylation arrays)"
                )
            else:
                report['warnings'].append(
                    f"{nan_count:,} NaN values found in {mat_name}"
                )

        # Data-type-specific validation
        if self.data_type == 'raw_counts':
            # Integer counts should not be negative
            neg_count = (self.counts_df < 0).sum().sum()
            if neg_count > 0:
                report['warnings'].append(
                    f"{neg_count} negative values in raw {mat_name}"
                )
        elif self.gene_id_type == 'probe':
            # Methylation beta values should be 0–1
            df_clean = self.counts_df.dropna(how='all')
            out_of_range = ((df_clean < 0) | (df_clean > 1)).sum().sum()
            if out_of_range > 0:
                report['warnings'].append(
                    f"{out_of_range} values outside expected 0–1 beta range"
                )

        return report

    def filter_genes(
        self,
        min_total_count: int = 10,
        min_samples_detected: int = 2,
        min_count_per_sample: int = 1,
    ) -> 'AcquiredDataset':
        """
        Filter low-expression genes.

        Parameters
        ----------
        min_total_count : int
            Minimum total count across all samples.
        min_samples_detected : int
            Minimum number of samples where gene must be detected.
        min_count_per_sample : int
            Minimum count to consider a gene 'detected' in a sample.

        Returns
        -------
        AcquiredDataset
            New dataset with filtered genes.
        """
        mask_total = self.counts_df.sum(axis=1) >= min_total_count
        mask_detected = (self.counts_df >= min_count_per_sample).sum(axis=1) >= min_samples_detected

        filtered_df = self.counts_df.loc[mask_total & mask_detected]

        n_removed = self.n_genes - len(filtered_df)
        logger.info(
            f"Filtered {n_removed:,} genes "
            f"({self.n_genes:,} → {len(filtered_df):,})"
        )

        return AcquiredDataset(
            counts_df=filtered_df,
            metadata=self.metadata.copy(),
            source_info={**self.source_info, 'filtered': True},
            gene_id_type=self.gene_id_type,
        )

    def subset_samples(self, sample_ids: List[str]) -> 'AcquiredDataset':
        """
        Create a subset with only the specified samples.

        Parameters
        ----------
        sample_ids : list of str
            Sample IDs to keep.

        Returns
        -------
        AcquiredDataset
            New dataset with subset of samples.
        """
        valid_ids = [s for s in sample_ids if s in self.counts_df.columns]
        if not valid_ids:
            raise ValidationError('sample_ids', "None of the specified sample IDs found in dataset", f"Available samples: {self.sample_ids[:5]}..."
            )

        return AcquiredDataset(
            counts_df=self.counts_df[valid_ids],
            metadata=self.metadata.loc[self.metadata.index.isin(valid_ids)],
            source_info=self.source_info.copy(),
            gene_id_type=self.gene_id_type,
        )

    def export_csv(self, output_dir: Union[str, Path]) -> Dict[str, Path]:
        """
        Export dataset as CSV files.

        Parameters
        ----------
        output_dir : str or Path
            Output directory.

        Returns
        -------
        dict
            Paths to exported files: counts_file, metadata_file, info_file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        prefix = self.accession.replace('/', '_')

        counts_path = output_dir / f"{prefix}_counts.csv"
        self.counts_df.to_csv(counts_path, encoding='utf-8')

        meta_path = output_dir / f"{prefix}_metadata.csv"
        self.metadata.to_csv(meta_path, encoding='utf-8')

        info_path = output_dir / f"{prefix}_info.json"
        with open(info_path, 'w', encoding='utf-8') as f:
            json.dump(self.source_info, f, indent=2, default=str)

        logger.info(f"Exported {prefix} to {output_dir}")

        return {
            'counts_file': counts_path,
            'metadata_file': meta_path,
            'info_file': info_path,
        }

    def summary(self) -> str:
        """Generate human-readable summary."""
        integrity = self.validate_integrity()
        status = "OK" if integrity['valid'] else "ISSUES FOUND"
        warnings_count = len(integrity['warnings'])

        lines = [
            "",
            "=" * 62,
            "  RAPTOR Acquired Dataset Summary",
            "=" * 62,
            "",
            f"  Repository:  {self.repository}",
            f"  Accession:   {self.accession}",
            f"  Organism:    {self.organism}",
            f"  Data type:   {self.data_type}",
            f"  Gene IDs:    {self.gene_id_type}",
            "",
            f"  Genes:       {self.n_genes:,}",
            f"  Samples:     {self.n_samples:,}",
            f"  With metadata: {self.samples_in_metadata:,}",
            "",
            f"  Integrity:   {status} ({warnings_count} warnings)",
            "",
            "=" * 62,
            "",
        ]
        return "\n".join(lines)

    def save(self, path: Union[str, Path]) -> None:
        """
        Save AcquiredDataset to disk.

        Parameters
        ----------
        path : str or Path
            Output path (.pkl or .parquet)
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        if path.suffix in ['.pkl', '.pickle']:
            with open(path, 'wb') as f:
                pickle.dump(self, f)
        elif path.suffix == '.parquet':
            # Save counts as parquet, metadata and info alongside
            self.counts_df.to_parquet(path, engine='pyarrow')
            meta_path = path.with_suffix('.metadata.parquet')
            self.metadata.to_parquet(meta_path, engine='pyarrow')
            info_path = path.with_suffix('.info.json')
            with open(info_path, 'w', encoding='utf-8') as f:
                json.dump(
                    {**self.source_info, 'gene_id_type': self.gene_id_type},
                    f, indent=2, default=str
                )
        else:
            raise ValidationError('path', f"Unsupported file extension: {path.suffix}", "Use .pkl, .pickle, or .parquet"
            )

        logger.info(f"AcquiredDataset saved to {path}")

    @classmethod
    def load(cls, path: Union[str, Path]) -> 'AcquiredDataset':
        """
        Load AcquiredDataset from disk.

        Parameters
        ----------
        path : str or Path
            Path to saved dataset (.pkl or .parquet)

        Returns
        -------
        AcquiredDataset
        """
        path = Path(path)
        if not path.exists():
            raise ValidationError('path', f"File not found: {path}", "Check the file path"
            )

        if path.suffix in ['.pkl', '.pickle']:
            with open(path, 'rb') as f:
                return pickle.load(f)
        elif path.suffix == '.parquet':
            counts_df = pd.read_parquet(path, engine='pyarrow')
            meta_path = path.with_suffix('.metadata.parquet')
            metadata = pd.read_parquet(meta_path, engine='pyarrow') if meta_path.exists() else pd.DataFrame()
            info_path = path.with_suffix('.info.json')
            if info_path.exists():
                with open(info_path, 'r', encoding='utf-8') as f:
                    info = json.load(f)
                gene_id_type = info.pop('gene_id_type', 'symbol')
            else:
                info = {}
                gene_id_type = 'symbol'
            return cls(
                counts_df=counts_df,
                metadata=metadata,
                source_info=info,
                gene_id_type=gene_id_type,
            )
        else:
            raise ValidationError('path', f"Unsupported extension: {path.suffix}", "Use .pkl, .pickle, or .parquet"
            )

    @classmethod
    def from_csv(
        cls,
        counts_file: Union[str, Path],
        metadata_file: Optional[Union[str, Path]] = None,
        gene_id_column: Optional[str] = None,
        **source_kwargs,
    ) -> 'AcquiredDataset':
        """
        Create AcquiredDataset from user-provided CSV files.

        This is the entry point for users bringing their own data
        into Module 10 without going through a repository connector.

        Parameters
        ----------
        counts_file : str or Path
            Path to count matrix CSV (genes as rows, samples as columns).
        metadata_file : str or Path, optional
            Path to sample metadata CSV.
        gene_id_column : str, optional
            Column name for gene IDs. If None, uses the first column or index.
        **source_kwargs
            Additional source_info fields (organism, data_type, etc.)

        Returns
        -------
        AcquiredDataset
        """
        counts_file = Path(counts_file)
        if not counts_file.exists():
            raise ValidationError('counts_file', f"File not found: {counts_file}", "Provide a valid path to a CSV count matrix"
            )

        counts_df = pd.read_csv(counts_file)

        # Detect and set gene ID index
        if gene_id_column and gene_id_column in counts_df.columns:
            counts_df = counts_df.set_index(gene_id_column)
        elif counts_df.columns[0].lower() in ['gene_id', 'geneid', 'gene', 'id']:
            counts_df = counts_df.set_index(counts_df.columns[0])
        else:
            counts_df.index.name = 'gene_id'

        # Load metadata if provided
        if metadata_file:
            metadata_file = Path(metadata_file)
            metadata = pd.read_csv(metadata_file, index_col=0)
        else:
            metadata = pd.DataFrame(index=counts_df.columns)
            metadata.index.name = 'sample_id'

        source_info = {
            'repository': 'user',
            'accession': counts_file.stem,
            'data_type': 'raw_counts',
            'source_file': str(counts_file),
            **source_kwargs,
        }

        return cls(
            counts_df=counts_df,
            metadata=metadata,
            source_info=source_info,
            gene_id_type=source_kwargs.get('gene_id_type', 'symbol'),
        )

    def __repr__(self) -> str:
        return (
            f"AcquiredDataset("
            f"accession='{self.accession}', "
            f"genes={self.n_genes:,}, "
            f"samples={self.n_samples}, "
            f"repository='{self.repository}')"
        )


# =============================================================================
# PooledDataset
# =============================================================================

@dataclass
class PooledDataset:
    """
    A dataset created by merging multiple AcquiredDatasets.

    Extends the dataset concept with study-level tracking for cross-study
    validation in Module 10 biomarker discovery.

    Attributes
    ----------
    counts_df : pd.DataFrame
        Merged count matrix (genes x samples) after harmonization.

    metadata : pd.DataFrame
        Merged sample metadata with a 'study' column indicating origin.

    study_labels : pd.Series
        Maps each sample to its source study (accession ID).

    component_datasets : list of dict
        Summary info for each dataset that was pooled.

    pooling_info : dict
        How the pooling was performed:
        - method : str (e.g., 'inner_join', 'outer_join')
        - batch_correction : str (e.g., 'combat', 'none')
        - gene_id_type : str (unified gene ID type)
        - n_studies : int
        - timestamp : str

    Examples
    --------
    >>> from raptor.external_modules.acquisition import PoolingEngine
    >>> engine = PoolingEngine()
    >>> pool = engine.merge([dataset1, dataset2, dataset3])
    >>> print(f"Pooled: {pool.n_studies} studies, {pool.n_samples} samples")
    >>> pool.save("my_pool.pkl")
    """

    counts_df: pd.DataFrame
    metadata: pd.DataFrame
    study_labels: pd.Series
    component_datasets: List[Dict[str, Any]] = field(default_factory=list)
    pooling_info: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate PooledDataset after initialization."""
        if not isinstance(self.counts_df, pd.DataFrame) or self.counts_df.empty:
            raise ValidationError('counts_df', "counts_df must be a non-empty DataFrame", "Run PoolingEngine.merge() to create a PooledDataset"
            )

        if self.counts_df.index.name is None:
            self.counts_df.index.name = 'gene_id'

        if self.metadata.index.name is None:
            self.metadata.index.name = 'sample_id'

        # Ensure study column exists in metadata
        if 'study' not in self.metadata.columns and not self.study_labels.empty:
            self.metadata['study'] = self.study_labels

        # Set defaults
        self.pooling_info.setdefault('method', 'unknown')
        self.pooling_info.setdefault('batch_correction', 'none')
        self.pooling_info.setdefault('gene_id_type', 'symbol')
        self.pooling_info.setdefault('n_studies', self.n_studies)
        self.pooling_info.setdefault('timestamp', datetime.now().isoformat())

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

    @property
    def n_genes(self) -> int:
        """Number of genes in the pooled matrix."""
        return len(self.counts_df)

    @property
    def n_samples(self) -> int:
        """Total number of samples across all studies."""
        return len(self.counts_df.columns)

    @property
    def n_studies(self) -> int:
        """Number of distinct studies pooled."""
        return self.study_labels.nunique()

    @property
    def studies(self) -> List[str]:
        """List of study accession IDs."""
        return sorted(self.study_labels.unique().tolist())

    @property
    def samples_per_study(self) -> Dict[str, int]:
        """Number of samples from each study."""
        return self.study_labels.value_counts().to_dict()

    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------

    def get_study_samples(self, study: str) -> List[str]:
        """
        Get sample IDs belonging to a specific study.

        Parameters
        ----------
        study : str
            Study accession ID.

        Returns
        -------
        list of str
            Sample IDs from that study.
        """
        return self.study_labels[self.study_labels == study].index.tolist()

    def leave_one_study_out(self, study: str) -> tuple:
        """
        Split into train (all other studies) and test (one study).

        Used for cross-study validation in Module 10.

        Parameters
        ----------
        study : str
            Study accession ID to hold out.

        Returns
        -------
        tuple of (pd.DataFrame, pd.DataFrame)
            (train_counts, test_counts) — both are gene x sample DataFrames.
        """
        test_samples = self.get_study_samples(study)
        train_samples = [s for s in self.counts_df.columns if s not in test_samples]

        return self.counts_df[train_samples], self.counts_df[test_samples]

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "",
            "=" * 62,
            "  RAPTOR Pooled Dataset Summary",
            "=" * 62,
            "",
            f"  Studies pooled:    {self.n_studies}",
            f"  Total genes:       {self.n_genes:,}",
            f"  Total samples:     {self.n_samples:,}",
            "",
            f"  Pooling method:    {self.pooling_info.get('method', 'unknown')}",
            f"  Batch correction:  {self.pooling_info.get('batch_correction', 'none')}",
            f"  Gene ID type:      {self.pooling_info.get('gene_id_type', 'unknown')}",
            "",
            "  Studies:",
        ]

        for study, count in self.samples_per_study.items():
            lines.append(f"    {study}: {count} samples")

        lines.extend(["", "=" * 62, ""])
        return "\n".join(lines)

    def save(self, path: Union[str, Path]) -> None:
        """
        Save PooledDataset to disk.

        Parameters
        ----------
        path : str or Path
            Output path (.pkl or .parquet)
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        if path.suffix in ['.pkl', '.pickle']:
            with open(path, 'wb') as f:
                pickle.dump(self, f)
        elif path.suffix == '.parquet':
            self.counts_df.to_parquet(path, engine='pyarrow')
            meta_path = path.with_suffix('.metadata.parquet')
            self.metadata.to_parquet(meta_path, engine='pyarrow')
            info_path = path.with_suffix('.pool_info.json')
            with open(info_path, 'w', encoding='utf-8') as f:
                json.dump({
                    'component_datasets': self.component_datasets,
                    'pooling_info': self.pooling_info,
                    'study_labels': self.study_labels.to_dict(),
                }, f, indent=2, default=str)
        else:
            raise ValidationError('path', f"Unsupported extension: {path.suffix}", "Use .pkl, .pickle, or .parquet"
            )

        logger.info(f"PooledDataset saved to {path}")

    @classmethod
    def load(cls, path: Union[str, Path]) -> 'PooledDataset':
        """
        Load PooledDataset from disk.

        Parameters
        ----------
        path : str or Path
            Path to saved pool (.pkl or .parquet)

        Returns
        -------
        PooledDataset
        """
        path = Path(path)
        if not path.exists():
            raise ValidationError('path', f"File not found: {path}", "Check the file path"
            )

        if path.suffix in ['.pkl', '.pickle']:
            with open(path, 'rb') as f:
                return pickle.load(f)
        elif path.suffix == '.parquet':
            counts_df = pd.read_parquet(path, engine='pyarrow')
            meta_path = path.with_suffix('.metadata.parquet')
            metadata = pd.read_parquet(meta_path, engine='pyarrow') if meta_path.exists() else pd.DataFrame()
            info_path = path.with_suffix('.pool_info.json')
            if info_path.exists():
                with open(info_path, 'r', encoding='utf-8') as f:
                    info = json.load(f)
                study_labels = pd.Series(info.get('study_labels', {}))
                component_datasets = info.get('component_datasets', [])
                pooling_info = info.get('pooling_info', {})
            else:
                study_labels = pd.Series(dtype=str)
                component_datasets = []
                pooling_info = {}
            return cls(
                counts_df=counts_df,
                metadata=metadata,
                study_labels=study_labels,
                component_datasets=component_datasets,
                pooling_info=pooling_info,
            )
        else:
            raise ValidationError('path', f"Unsupported extension: {path.suffix}", "Use .pkl, .pickle, or .parquet"
            )

    def __repr__(self) -> str:
        return (
            f"PooledDataset("
            f"studies={self.n_studies}, "
            f"genes={self.n_genes:,}, "
            f"samples={self.n_samples})"
        )
