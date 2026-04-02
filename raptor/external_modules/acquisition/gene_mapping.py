"""
RAPTOR v2.2.2 - Module 6: Data Acquisition - Gene ID Mapper

Maps gene identifiers between Ensembl, HGNC symbols, and Entrez IDs.
Essential for pooling datasets from different repositories that use
different gene ID systems.

Uses MyGene.info as the backend — a free, high-performance gene
annotation service that requires no API key.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import warnings
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd
import numpy as np

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

VALID_ID_TYPES = ['ensembl', 'symbol', 'entrez']

# Scope fields for MyGene.info queries
MYGENE_SCOPES = {
    'ensembl': 'ensembl.gene',
    'symbol': 'symbol',
    'entrez': 'entrezgene',
}

MYGENE_FIELDS = 'symbol,ensembl.gene,entrezgene,name'

# Common species taxonomy IDs
SPECIES_TAXID = {
    'homo sapiens': 9606,
    'human': 9606,
    'mus musculus': 10090,
    'mouse': 10090,
    'rattus norvegicus': 10116,
    'rat': 10116,
    'drosophila melanogaster': 7227,
    'danio rerio': 7955,
    'caenorhabditis elegans': 6239,
}


# =============================================================================
# GeneIDMapper
# =============================================================================

class GeneIDMapper:
    """
    Convert gene identifiers between Ensembl, HGNC symbol, and Entrez ID.

    Uses MyGene.info for mapping. Results are cached in memory for the
    session to avoid redundant API calls.

    Parameters
    ----------
    species : str
        Species name (e.g., 'Homo sapiens', 'human', 'Mus musculus').

    Examples
    --------
    >>> mapper = GeneIDMapper('Homo sapiens')
    >>> mapped = mapper.convert(['TP53', 'BRCA1', 'EGFR'], 'symbol', 'ensembl')
    >>> print(mapped)
    {'TP53': 'ENSG00000141510', 'BRCA1': 'ENSG00000012048', ...}

    >>> df = mapper.convert_index(counts_df, from_type='ensembl', to_type='symbol')
    """

    def __init__(self, species: str = 'Homo sapiens'):
        self._species = species
        self._taxid = self._resolve_taxid(species)
        self._mg = None  # lazy load mygene
        self._cache: Dict[str, Dict[str, str]] = {}

    # -------------------------------------------------------------------------
    # Setup
    # -------------------------------------------------------------------------

    @staticmethod
    def _resolve_taxid(species: str) -> int:
        """Resolve species name to NCBI taxonomy ID."""
        taxid = SPECIES_TAXID.get(species.lower())
        if taxid is None:
            logger.warning(
                f"Unknown species '{species}', defaulting to human (9606). "
                f"Known species: {list(SPECIES_TAXID.keys())}"
            )
            return 9606
        return taxid

    def _get_mygene(self):
        """Lazy-load the mygene client."""
        if self._mg is None:
            try:
                import mygene
                self._mg = mygene.MyGeneInfo()
                logger.debug("MyGene.info client initialized")
            except ImportError:
                raise RAPTORError(
                    "mygene package is required for gene ID mapping. "
                    "Install it with: pip install mygene"
                )
        return self._mg

    # -------------------------------------------------------------------------
    # Core mapping
    # -------------------------------------------------------------------------

    def convert(
        self,
        gene_ids: List[str],
        from_type: str,
        to_type: str,
        drop_unmapped: bool = False,
    ) -> Dict[str, Optional[str]]:
        """
        Convert a list of gene IDs from one type to another.

        Parameters
        ----------
        gene_ids : list of str
            Input gene identifiers.
        from_type : str
            Source ID type: 'ensembl', 'symbol', or 'entrez'.
        to_type : str
            Target ID type: 'ensembl', 'symbol', or 'entrez'.
        drop_unmapped : bool
            If True, omit genes that couldn't be mapped.

        Returns
        -------
        dict
            Mapping from input ID to target ID (or None if unmapped).
        """
        self._validate_id_type(from_type, 'from_type')
        self._validate_id_type(to_type, 'to_type')

        if from_type == to_type:
            return {gid: gid for gid in gene_ids}

        # Check cache
        cache_key = f"{from_type}_to_{to_type}"
        if cache_key not in self._cache:
            self._cache[cache_key] = {}

        cached = self._cache[cache_key]
        uncached = [gid for gid in gene_ids if gid not in cached]

        if uncached:
            new_mappings = self._query_mygene(uncached, from_type, to_type)
            cached.update(new_mappings)

        if drop_unmapped:
            return {
                gid: cached.get(gid)
                for gid in gene_ids
                if cached.get(gid) is not None
            }

        return {gid: cached.get(gid) for gid in gene_ids}

    def _query_mygene(
        self,
        gene_ids: List[str],
        from_type: str,
        to_type: str,
    ) -> Dict[str, Optional[str]]:
        """
        Query MyGene.info for ID mappings.

        Parameters
        ----------
        gene_ids : list of str
            Gene IDs to query.
        from_type : str
            Source type.
        to_type : str
            Target type.

        Returns
        -------
        dict
            Mapping from input to target ID.
        """
        mg = self._get_mygene()
        scope = MYGENE_SCOPES[from_type]
        target_field = MYGENE_SCOPES[to_type]

        logger.info(
            f"Querying MyGene.info: {len(gene_ids)} genes "
            f"({from_type} → {to_type})"
        )

        # MyGene.info querymany handles batching internally
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            results = mg.querymany(
                gene_ids,
                scopes=scope,
                fields=target_field,
                species=self._taxid,
                returnall=True,
            )

        mapping = {}
        for hit in results.get('out', []):
            query_id = str(hit.get('query', ''))
            if hit.get('notfound', False):
                mapping[query_id] = None
                continue

            # Extract target value
            target_val = self._extract_field(hit, to_type)
            mapping[query_id] = target_val

        # Log mapping stats
        mapped = sum(1 for v in mapping.values() if v is not None)
        total = len(gene_ids)
        logger.info(
            f"Mapped {mapped}/{total} genes "
            f"({100 * mapped / total:.1f}%)"
        )

        if results.get('missing'):
            logger.debug(
                f"{len(results['missing'])} genes returned no results"
            )

        return mapping

    @staticmethod
    def _extract_field(hit: dict, id_type: str) -> Optional[str]:
        """Extract a gene ID field from a MyGene.info result."""
        if id_type == 'symbol':
            return hit.get('symbol')

        elif id_type == 'entrez':
            val = hit.get('entrezgene')
            return str(val) if val is not None else None

        elif id_type == 'ensembl':
            ensembl = hit.get('ensembl')
            if ensembl is None:
                return None
            # Can be a dict or a list of dicts
            if isinstance(ensembl, list):
                return ensembl[0].get('gene') if ensembl else None
            return ensembl.get('gene')

        return None

    # -------------------------------------------------------------------------
    # DataFrame operations
    # -------------------------------------------------------------------------

    def convert_index(
        self,
        df: pd.DataFrame,
        from_type: str,
        to_type: str,
        drop_unmapped: bool = True,
        aggregate: str = 'mean',
    ) -> pd.DataFrame:
        """
        Convert a DataFrame's index from one gene ID type to another.

        Handles duplicates (multiple source IDs mapping to the same
        target) by aggregating with the specified method.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with gene IDs as index.
        from_type : str
            Current index ID type.
        to_type : str
            Target index ID type.
        drop_unmapped : bool
            Drop genes that couldn't be mapped.
        aggregate : str
            How to handle duplicates: 'mean', 'sum', 'max', 'first'.

        Returns
        -------
        pd.DataFrame
            DataFrame with converted index.
        """
        gene_ids = df.index.astype(str).tolist()
        mapping = self.convert(gene_ids, from_type, to_type)

        # Build new index
        new_index = [mapping.get(gid) for gid in gene_ids]
        df_mapped = df.copy()
        df_mapped.index = new_index
        df_mapped.index.name = 'gene_id'

        # Drop unmapped
        if drop_unmapped:
            df_mapped = df_mapped.loc[df_mapped.index.notna()]

        # Handle duplicates
        if df_mapped.index.duplicated().any():
            n_dups = df_mapped.index.duplicated().sum()
            logger.info(
                f"Aggregating {n_dups} duplicate mappings using '{aggregate}'"
            )
            if aggregate == 'mean':
                df_mapped = df_mapped.groupby(level=0).mean()
            elif aggregate == 'sum':
                df_mapped = df_mapped.groupby(level=0).sum()
            elif aggregate == 'max':
                df_mapped = df_mapped.groupby(level=0).max()
            elif aggregate == 'first':
                df_mapped = df_mapped.loc[~df_mapped.index.duplicated(keep='first')]
            else:
                logger.warning(f"Unknown aggregate '{aggregate}', using 'mean'")
                df_mapped = df_mapped.groupby(level=0).mean()

        return df_mapped

    def harmonize_to_common(
        self,
        datasets: list,
        target_type: str = 'symbol',
    ) -> list:
        """
        Convert multiple datasets to a common gene ID type.

        Parameters
        ----------
        datasets : list of AcquiredDataset
            Datasets to harmonize.
        target_type : str
            Target gene ID type.

        Returns
        -------
        list of AcquiredDataset
            Datasets with harmonized gene IDs.
        """
        from .datasets import AcquiredDataset

        harmonized = []
        for ds in datasets:
            if ds.gene_id_type == target_type:
                harmonized.append(ds)
                logger.info(
                    f"{ds.accession}: already {target_type}, skipping conversion"
                )
                continue

            logger.info(
                f"{ds.accession}: converting {ds.gene_id_type} → {target_type}"
            )
            new_counts = self.convert_index(
                ds.counts_df,
                from_type=ds.gene_id_type,
                to_type=target_type,
            )

            new_ds = AcquiredDataset(
                counts_df=new_counts,
                metadata=ds.metadata.copy(),
                source_info={
                    **ds.source_info,
                    'original_gene_id_type': ds.gene_id_type,
                },
                gene_id_type=target_type,
            )
            harmonized.append(new_ds)

        return harmonized

    # -------------------------------------------------------------------------
    # Detection
    # -------------------------------------------------------------------------

    @staticmethod
    def detect_id_type(gene_ids: List[str]) -> str:
        """
        Heuristically detect the gene ID type from a list of IDs.

        Parameters
        ----------
        gene_ids : list of str
            Sample of gene identifiers.

        Returns
        -------
        str
            Detected type: 'ensembl', 'entrez', 'symbol', or 'unknown'.
        """
        if not gene_ids:
            return 'unknown'

        sample = [str(gid) for gid in gene_ids[:100]]

        # Ensembl: starts with ENS
        ensembl_count = sum(
            1 for gid in sample if gid.startswith('ENS')
        )
        if ensembl_count > len(sample) * 0.5:
            return 'ensembl'

        # Entrez: all numeric
        entrez_count = sum(1 for gid in sample if gid.isdigit())
        if entrez_count > len(sample) * 0.5:
            return 'entrez'

        # Symbol: alphabetic, typically short, mixed case
        alpha_count = sum(
            1 for gid in sample
            if gid[0].isalpha() and len(gid) < 20
        )
        if alpha_count > len(sample) * 0.5:
            return 'symbol'

        return 'unknown'

    # -------------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------------

    @staticmethod
    def _validate_id_type(id_type: str, param_name: str) -> None:
        """Validate that id_type is supported."""
        if id_type not in VALID_ID_TYPES:
            raise ValidationError(param_name, f"Invalid gene ID type: '{id_type}'", f"Must be one of: {VALID_ID_TYPES}"
            )

    def clear_cache(self) -> None:
        """Clear the in-memory mapping cache."""
        self._cache.clear()
        logger.debug("Gene mapping cache cleared")

    @property
    def species(self) -> str:
        """Current species."""
        return self._species

    @property
    def taxid(self) -> int:
        """NCBI taxonomy ID."""
        return self._taxid

    def __repr__(self) -> str:
        cached_mappings = sum(len(v) for v in self._cache.values())
        return (
            f"GeneIDMapper("
            f"species='{self._species}', "
            f"taxid={self._taxid}, "
            f"cached_mappings={cached_mappings})"
        )
