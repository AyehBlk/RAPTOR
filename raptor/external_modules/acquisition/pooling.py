"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - Pooling Engine

Merges multiple AcquiredDatasets into a single PooledDataset with
gene ID harmonization, batch effect correction, and metadata merging.
The output feeds directly into Module 10 for cross-study biomarker discovery.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
from datetime import datetime
from typing import Dict, List, Optional, Any, Union

import numpy as np
import pandas as pd

from .datasets import AcquiredDataset, PooledDataset
from .gene_mapping import GeneIDMapper

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

SUPPORTED_JOIN_METHODS = ['inner', 'outer']

SUPPORTED_BATCH_METHODS = ['none', 'combat', 'quantile', 'median_ratio']


# =============================================================================
# PoolingEngine
# =============================================================================

class PoolingEngine:
    """
    Merge multiple AcquiredDatasets into a PooledDataset.

    Handles gene ID harmonization, count matrix merging, batch effect
    correction, and metadata consolidation.

    Parameters
    ----------
    target_gene_id : str
        Gene ID type for harmonization: 'symbol', 'ensembl', 'entrez'.
    species : str
        Species for gene ID mapping.

    Examples
    --------
    >>> engine = PoolingEngine(target_gene_id='symbol')
    >>> pool = engine.merge(
    ...     datasets=[dataset1, dataset2, dataset3],
    ...     method='inner',
    ...     batch_correction='combat',
    ... )
    >>> print(pool.summary())
    >>> pool.save("my_pool.pkl")
    """

    def __init__(
        self,
        target_gene_id: str = 'symbol',
        species: str = 'Homo sapiens',
    ):
        valid_types = ['ensembl', 'symbol', 'entrez']
        if target_gene_id not in valid_types:
            raise ValidationError('target_gene_id', f"Invalid gene ID type: '{target_gene_id}'", f"Must be one of: {valid_types}"
            )

        self._target_gene_id = target_gene_id
        self._species = species
        self._mapper = GeneIDMapper(species=species)

    # -------------------------------------------------------------------------
    # Main merge
    # -------------------------------------------------------------------------

    def merge(
        self,
        datasets: List[AcquiredDataset],
        pool_name: Optional[str] = None,
        method: str = 'inner',
        batch_correction: str = 'none',
        min_common_genes: int = 100,
    ) -> PooledDataset:
        """
        Merge multiple datasets into a PooledDataset.

        Parameters
        ----------
        datasets : list of AcquiredDataset
            Datasets to merge (minimum 2).
        pool_name : str, optional
            Name for the pool. Auto-generated if not provided.
        method : str
            Join method: 'inner' (common genes only) or 'outer'
            (all genes, NaN-filled).
        batch_correction : str
            Batch correction: 'none', 'combat', 'quantile', 'median_ratio'.
        min_common_genes : int
            Minimum number of common genes required for merging.

        Returns
        -------
        PooledDataset
            Merged dataset ready for Module 10.
        """
        # Validate inputs
        if len(datasets) < 2:
            raise ValidationError('datasets', f"Need at least 2 datasets to pool, got {len(datasets)}", "Pass a list of AcquiredDataset objects"
            )

        if method not in SUPPORTED_JOIN_METHODS:
            raise ValidationError('method', f"Unsupported join method: '{method}'", f"Use one of: {SUPPORTED_JOIN_METHODS}"
            )

        if batch_correction not in SUPPORTED_BATCH_METHODS:
            raise ValidationError('batch_correction', f"Unsupported batch correction: '{batch_correction}'", f"Use one of: {SUPPORTED_BATCH_METHODS}"
            )

        # Auto-generate pool name
        if pool_name is None:
            accessions = '_'.join(ds.accession for ds in datasets[:3])
            pool_name = f"pool_{accessions}_{datetime.now():%Y%m%d}"

        logger.info(
            f"Pooling {len(datasets)} datasets "
            f"(method={method}, batch_correction={batch_correction})"
        )

        # Step 1: Harmonize gene IDs
        logger.info("Step 1/4: Harmonizing gene IDs...")
        harmonized = self._harmonize_gene_ids(datasets)

        # Step 2: Resolve sample ID conflicts
        logger.info("Step 2/4: Resolving sample IDs...")
        harmonized = self._resolve_sample_conflicts(harmonized)

        # Step 3: Merge count matrices
        logger.info("Step 3/4: Merging count matrices...")
        counts_df, study_labels = self._merge_counts(
            harmonized, method, min_common_genes
        )

        # Step 4: Apply batch correction
        if batch_correction != 'none':
            logger.info(f"Step 4/4: Applying {batch_correction} correction...")
            counts_df = self._apply_batch_correction(
                counts_df, study_labels, batch_correction
            )
        else:
            logger.info("Step 4/4: No batch correction applied")

        # Merge metadata
        metadata = self._merge_metadata(harmonized, study_labels)

        # Build component info
        component_datasets = [
            {
                'accession': ds.accession,
                'repository': ds.repository,
                'organism': ds.organism,
                'n_genes': ds.n_genes,
                'n_samples': ds.n_samples,
                'gene_id_type': ds.gene_id_type,
            }
            for ds in datasets
        ]

        pooling_info = {
            'method': method,
            'batch_correction': batch_correction,
            'gene_id_type': self._target_gene_id,
            'n_studies': len(datasets),
            'min_common_genes': min_common_genes,
            'timestamp': datetime.now().isoformat(),
            'species': self._species,
        }

        pool = PooledDataset(
            counts_df=counts_df,
            metadata=metadata,
            study_labels=study_labels,
            component_datasets=component_datasets,
            pooling_info=pooling_info,
        )

        logger.info(
            f"Pooling complete: {pool.n_genes:,} genes, "
            f"{pool.n_samples} samples from {pool.n_studies} studies"
        )

        return pool

    # -------------------------------------------------------------------------
    # Step 1: Gene ID harmonization
    # -------------------------------------------------------------------------

    def _harmonize_gene_ids(
        self,
        datasets: List[AcquiredDataset],
    ) -> List[AcquiredDataset]:
        """Convert all datasets to the target gene ID type."""
        needs_conversion = any(
            ds.gene_id_type != self._target_gene_id
            for ds in datasets
        )

        if not needs_conversion:
            logger.info(
                f"All datasets already use {self._target_gene_id} IDs"
            )
            return datasets

        return self._mapper.harmonize_to_common(
            datasets, target_type=self._target_gene_id
        )

    # -------------------------------------------------------------------------
    # Step 2: Sample ID conflict resolution
    # -------------------------------------------------------------------------

    @staticmethod
    def _resolve_sample_conflicts(
        datasets: List[AcquiredDataset],
    ) -> List[AcquiredDataset]:
        """
        Ensure unique sample IDs across all datasets.

        If two datasets share sample IDs, prefix with accession.
        """
        all_sample_ids = []
        for ds in datasets:
            all_sample_ids.extend(ds.sample_ids)

        if len(all_sample_ids) == len(set(all_sample_ids)):
            return datasets  # no conflicts

        logger.info("Sample ID conflicts detected, adding study prefixes...")

        resolved = []
        for ds in datasets:
            prefix = ds.accession.replace('-', '_')

            new_cols = {
                col: f"{prefix}__{col}" for col in ds.counts_df.columns
            }
            new_counts = ds.counts_df.rename(columns=new_cols)

            new_meta = ds.metadata.copy()
            new_meta.index = [
                f"{prefix}__{sid}" for sid in new_meta.index
            ]

            resolved.append(AcquiredDataset(
                counts_df=new_counts,
                metadata=new_meta,
                source_info=ds.source_info.copy(),
                gene_id_type=ds.gene_id_type,
            ))

        return resolved

    # -------------------------------------------------------------------------
    # Step 3: Merge count matrices
    # -------------------------------------------------------------------------

    def _merge_counts(
        self,
        datasets: List[AcquiredDataset],
        method: str,
        min_common_genes: int,
    ) -> tuple:
        """Merge count matrices and build study labels."""
        # Collect all DataFrames
        dfs = []
        labels = {}
        for ds in datasets:
            dfs.append(ds.counts_df)
            for sample in ds.counts_df.columns:
                labels[sample] = ds.accession

        study_labels = pd.Series(labels)

        # Check gene overlap
        gene_sets = [set(df.index) for df in dfs]
        common_genes = gene_sets[0]
        for gs in gene_sets[1:]:
            common_genes = common_genes & gs

        all_genes = set()
        for gs in gene_sets:
            all_genes = all_genes | gs

        logger.info(
            f"Gene overlap: {len(common_genes):,} common / "
            f"{len(all_genes):,} total across {len(datasets)} datasets"
        )

        if len(common_genes) < min_common_genes:
            raise ValidationError('datasets', (
                    f"Only {len(common_genes)} common genes found "
                    f"(minimum: {min_common_genes})"
                ), (
                    "Datasets may use different gene ID types or have "
                    "very different gene coverage. Try using 'outer' join "
                    "or check gene_id_type consistency."
                )
            )

        # Merge
        if method == 'inner':
            # Keep only common genes
            gene_list = sorted(common_genes)
            merged = pd.concat(
                [df.loc[df.index.isin(common_genes)] for df in dfs],
                axis=1,
            )
            # Ensure consistent gene order
            merged = merged.loc[merged.index.isin(gene_list)]

        elif method == 'outer':
            merged = pd.concat(dfs, axis=1)
            merged = merged.fillna(0)

        merged.index.name = 'gene_id'

        # Remove duplicate columns (shouldn't happen after conflict resolution)
        if merged.columns.duplicated().any():
            logger.warning("Removing duplicate sample columns")
            merged = merged.loc[:, ~merged.columns.duplicated()]

        logger.info(
            f"Merged matrix: {merged.shape[0]:,} genes x "
            f"{merged.shape[1]} samples"
        )

        return merged, study_labels

    # -------------------------------------------------------------------------
    # Step 4: Batch correction
    # -------------------------------------------------------------------------

    def _apply_batch_correction(
        self,
        counts_df: pd.DataFrame,
        study_labels: pd.Series,
        method: str,
    ) -> pd.DataFrame:
        """Apply batch effect correction."""
        if method == 'combat':
            return self._combat_correction(counts_df, study_labels)
        elif method == 'quantile':
            return self._quantile_normalization(counts_df)
        elif method == 'median_ratio':
            return self._median_ratio_normalization(counts_df, study_labels)
        else:
            return counts_df

    def _combat_correction(
        self,
        counts_df: pd.DataFrame,
        study_labels: pd.Series,
    ) -> pd.DataFrame:
        """
        Apply ComBat batch correction.

        Uses a pure-Python implementation if pycombat is available,
        otherwise falls back to a simpler median-centering approach.
        """
        try:
            from combat.pycombat import pycombat

            # pycombat expects genes x samples with batch as a list
            batch = study_labels.loc[counts_df.columns].tolist()

            logger.info("Applying ComBat batch correction (pycombat)...")
            corrected = pycombat(counts_df, batch)
            return corrected

        except ImportError:
            logger.warning(
                "pycombat not installed, using median-centering fallback. "
                "For full ComBat: pip install combat"
            )
            return self._median_centering(counts_df, study_labels)

    @staticmethod
    def _median_centering(
        counts_df: pd.DataFrame,
        study_labels: pd.Series,
    ) -> pd.DataFrame:
        """
        Simple median-centering per study.

        For each study, subtract the per-gene median so that each study
        has the same center. A lightweight batch correction alternative.
        """
        corrected = counts_df.copy()

        # Log-transform for centering (add pseudocount)
        log_counts = np.log2(corrected + 1)

        # Global median per gene
        global_median = log_counts.median(axis=1)

        for study in study_labels.unique():
            study_samples = study_labels[study_labels == study].index
            study_cols = [s for s in study_samples if s in log_counts.columns]

            if not study_cols:
                continue

            study_median = log_counts[study_cols].median(axis=1)
            shift = global_median - study_median

            for col in study_cols:
                log_counts[col] = log_counts[col] + shift

        # Back-transform
        corrected = (2 ** log_counts) - 1
        corrected = corrected.clip(lower=0)

        return corrected

    @staticmethod
    def _quantile_normalization(counts_df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply quantile normalization across all samples.

        Forces all samples to have the same distribution of values.
        """
        logger.info("Applying quantile normalization...")

        ranked = counts_df.rank(method='average')
        sorted_vals = np.sort(counts_df.values, axis=0)
        mean_sorted = sorted_vals.mean(axis=1)

        normalized = counts_df.copy()
        for col in normalized.columns:
            ranks = ranked[col].values
            # Map ranks to mean values
            normalized[col] = np.interp(
                ranks,
                np.arange(1, len(mean_sorted) + 1),
                mean_sorted,
            )

        return normalized

    @staticmethod
    def _median_ratio_normalization(
        counts_df: pd.DataFrame,
        study_labels: pd.Series,
    ) -> pd.DataFrame:
        """
        DESeq2-style median ratio normalization per study.

        Computes size factors within each study to make libraries
        comparable across studies.
        """
        logger.info("Applying median ratio normalization...")

        normalized = counts_df.copy()

        # Geometric mean per gene (reference)
        log_counts = np.log(counts_df.replace(0, np.nan))
        geo_mean = np.exp(log_counts.mean(axis=1))

        for study in study_labels.unique():
            study_samples = study_labels[study_labels == study].index
            study_cols = [s for s in study_samples if s in counts_df.columns]

            if not study_cols:
                continue

            for col in study_cols:
                ratios = counts_df[col] / geo_mean
                ratios = ratios.replace([np.inf, -np.inf], np.nan)
                size_factor = ratios.median(skipna=True)

                if size_factor > 0 and np.isfinite(size_factor):
                    normalized[col] = counts_df[col] / size_factor

        return normalized

    # -------------------------------------------------------------------------
    # Metadata merge
    # -------------------------------------------------------------------------

    @staticmethod
    def _merge_metadata(
        datasets: List[AcquiredDataset],
        study_labels: pd.Series,
    ) -> pd.DataFrame:
        """Merge metadata from all datasets."""
        meta_dfs = []

        for ds in datasets:
            meta = ds.metadata.copy()
            # Add study column
            meta['study'] = ds.accession
            meta['repository'] = ds.repository
            meta['organism'] = ds.organism
            meta_dfs.append(meta)

        if not meta_dfs:
            metadata = pd.DataFrame(index=study_labels.index)
            metadata['study'] = study_labels
            return metadata

        # Concat with outer join (different datasets may have different columns)
        metadata = pd.concat(meta_dfs, axis=0, sort=False)
        metadata.index.name = 'sample_id'

        return metadata

    # -------------------------------------------------------------------------
    # Display
    # -------------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"PoolingEngine("
            f"target='{self._target_gene_id}', "
            f"species='{self._species}')"
        )
