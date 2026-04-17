"""
Ratio Biomarker Module
======================

Gene-expression ratios (geneA / geneB) are powerful biomarkers because they
are **self-normalizing**: batch effects, library-size differences, and
platform variation cancel out when you divide one gene by another measured
in the same sample.

This module implements a generalized Top Scoring Pair (TSP) search
(Geman et al., 2004):

1. ``RatioBiomarkerSearcher``  — tests all (or filtered) pairwise gene
   ratios, ranks them by AUC, and returns the top-k most discriminating
   pairs.
2. ``apply_ratios``            — computes learned ratio features on new
   patient data.
3. ``build_ratio_features``    — convenience function that wraps the
   full search-then-apply workflow.

Usage
-----
>>> searcher = RatioBiomarkerSearcher(top_k=10, min_auc=0.80)
>>> results = searcher.search(expression, labels)
>>> results.best_pair
RatioPair(gene_a='MYC', gene_b='TP53', auc=0.94, direction='higher_in_disease')
>>> new_ratios = apply_ratios(new_expression, results.pairs)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class RatioPair:
    """A single gene-ratio biomarker: gene_a / gene_b."""

    gene_a: str
    gene_b: str
    auc: float
    direction: str          # 'higher_in_reference' or 'higher_in_baseline'
    mean_ratio_reference: float
    mean_ratio_baseline: float

    @property
    def name(self) -> str:
        return f"{self.gene_a}/{self.gene_b}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gene_a": self.gene_a,
            "gene_b": self.gene_b,
            "name": self.name,
            "auc": self.auc,
            "direction": self.direction,
            "mean_ratio_reference": self.mean_ratio_reference,
            "mean_ratio_baseline": self.mean_ratio_baseline,
        }


@dataclass
class RatioSearchResult:
    """Container for the full ratio biomarker search output."""

    pairs: List[RatioPair]
    reference_group: str
    baseline_group: str
    n_pairs_tested: int
    n_genes_input: int
    top_k: int
    min_auc: float

    @property
    def best_pair(self) -> Optional[RatioPair]:
        return self.pairs[0] if self.pairs else None

    @property
    def n_found(self) -> int:
        return len(self.pairs)

    def to_dataframe(self) -> pd.DataFrame:
        if not self.pairs:
            return pd.DataFrame(
                columns=["gene_a", "gene_b", "name", "auc", "direction",
                          "mean_ratio_reference", "mean_ratio_baseline"]
            )
        return pd.DataFrame([p.to_dict() for p in self.pairs])

    def summary(self) -> str:
        lines = [
            f"Ratio Biomarker Search Results",
            f"  Groups: {self.reference_group} vs {self.baseline_group}",
            f"  Genes input:    {self.n_genes_input}",
            f"  Pairs tested:   {self.n_pairs_tested}",
            f"  Pairs found:    {self.n_found} (AUC >= {self.min_auc})",
        ]
        if self.best_pair:
            lines.append(
                f"  Best pair:      {self.best_pair.name} "
                f"(AUC = {self.best_pair.auc:.4f})"
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main searcher class
# ---------------------------------------------------------------------------

@dataclass
class RatioBiomarkerSearcher:
    """
    Search all pairwise gene ratios for discriminating biomarkers.

    Parameters
    ----------
    top_k : int, default 20
        Maximum number of top-ranked ratio pairs to return.
    min_auc : float, default 0.70
        Minimum AUC for a ratio pair to be included in results.
    min_variance : float, default 1e-6
        Genes with variance below this across all samples are excluded
        (they would produce near-constant ratios).
    pseudocount : float, default 1.0
        Added to both numerator and denominator genes before dividing,
        to avoid division by zero and stabilize ratios for low-expressed
        genes.  Set to 0 to disable (not recommended).
    candidate_genes : list of str, optional
        If provided, only these genes are considered. Useful for
        restricting the search to DE genes or a prior gene set.
    """

    top_k: int = 20
    min_auc: float = 0.70
    min_variance: float = 1e-6
    pseudocount: float = 1.0
    candidate_genes: Optional[List[str]] = None

    def __post_init__(self) -> None:
        if self.top_k < 1:
            raise ValueError(f"top_k must be >= 1, got {self.top_k}.")
        if not 0.5 <= self.min_auc <= 1.0:
            raise ValueError(f"min_auc must be in [0.5, 1.0], got {self.min_auc}.")
        if self.pseudocount < 0:
            raise ValueError(f"pseudocount must be >= 0, got {self.pseudocount}.")

    def search(
        self,
        expression: pd.DataFrame,
        labels: Union[pd.Series, np.ndarray, list],
        reference_group: str,
        baseline_group: str,
    ) -> RatioSearchResult:
        """
        Test all pairwise gene ratios and return top discriminators.

        Parameters
        ----------
        expression : pd.DataFrame (samples × genes)
            Expression matrix. Values should be non-negative (e.g. TPM,
            CPM, or normalized counts). Log-transformed data works but
            ratios of logs are less interpretable than ratios of linear
            values.
        labels : Series, array, or list
            Group label per sample.
        reference_group : str
            The group of interest (e.g. 'disease').
        baseline_group : str
            The control group (e.g. 'healthy').

        Returns
        -------
        RatioSearchResult
        """
        labels = pd.Series(labels).reset_index(drop=True)
        expression = expression.reset_index(drop=True)

        if len(labels) != len(expression):
            raise ValueError(
                f"labels length ({len(labels)}) != expression rows "
                f"({len(expression)})."
            )

        for g in (reference_group, baseline_group):
            if g not in set(labels.unique()):
                raise ValueError(
                    f"Group {g!r} not in labels. "
                    f"Available: {sorted(labels.unique())}."
                )

        # --- select genes ---
        if self.candidate_genes is not None:
            available = [g for g in self.candidate_genes
                         if g in expression.columns]
            if len(available) < 2:
                raise ValueError(
                    f"Need at least 2 candidate genes in expression columns, "
                    f"found {len(available)}."
                )
            missing = set(self.candidate_genes) - set(available)
            if missing:
                warnings.warn(
                    f"{len(missing)} candidate gene(s) not found in "
                    f"expression columns and will be skipped.",
                    UserWarning,
                    stacklevel=2,
                )
            genes = available
        else:
            genes = list(expression.columns)

        # --- filter low-variance genes ---
        variances = expression[genes].var()
        genes = [g for g in genes if variances[g] > self.min_variance]

        if len(genes) < 2:
            raise ValueError(
                f"Need at least 2 genes after variance filtering, "
                f"got {len(genes)}."
            )

        n_genes_input = len(genes)

        # --- build binary label vector (1 = reference, 0 = baseline) ---
        mask = labels.isin([reference_group, baseline_group])
        y = (labels[mask] == reference_group).astype(int).to_numpy()
        expr_sub = expression.loc[mask, genes].to_numpy(dtype=float)

        # gene name lookup
        gene_idx = {name: i for i, name in enumerate(genes)}

        # --- test all pairs ---
        scored: List[Tuple[float, str, str, str, float, float]] = []
        n_pairs = 0

        for ga, gb in combinations(genes, 2):
            ia, ib = gene_idx[ga], gene_idx[gb]
            ratio = (expr_sub[:, ia] + self.pseudocount) / (
                expr_sub[:, ib] + self.pseudocount
            )

            # Skip if ratio is constant
            if np.std(ratio) < 1e-10:
                n_pairs += 1
                continue

            try:
                auc = float(roc_auc_score(y, ratio))
            except Exception:
                n_pairs += 1
                continue

            # AUC < 0.5 means the ratio discriminates in reverse
            if auc < 0.5:
                auc = 1 - auc
                direction = "higher_in_baseline"
            else:
                direction = "higher_in_reference"

            n_pairs += 1

            if auc >= self.min_auc:
                ref_mean = float(np.mean(ratio[y == 1]))
                base_mean = float(np.mean(ratio[y == 0]))
                scored.append((auc, ga, gb, direction, ref_mean, base_mean))

        # --- rank and truncate ---
        scored.sort(key=lambda x: x[0], reverse=True)
        top = scored[: self.top_k]

        pairs = [
            RatioPair(
                gene_a=ga, gene_b=gb, auc=auc,
                direction=direction,
                mean_ratio_reference=ref_mean,
                mean_ratio_baseline=base_mean,
            )
            for auc, ga, gb, direction, ref_mean, base_mean in top
        ]

        return RatioSearchResult(
            pairs=pairs,
            reference_group=reference_group,
            baseline_group=baseline_group,
            n_pairs_tested=n_pairs,
            n_genes_input=n_genes_input,
            top_k=self.top_k,
            min_auc=self.min_auc,
        )


# ---------------------------------------------------------------------------
# Apply learned ratios to new data
# ---------------------------------------------------------------------------

def apply_ratios(
    expression: pd.DataFrame,
    pairs: List[RatioPair],
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    """
    Compute ratio features on new data using previously discovered pairs.

    Parameters
    ----------
    expression : pd.DataFrame (samples × genes)
        New expression data.
    pairs : list of RatioPair
        Pairs from a previous ``RatioBiomarkerSearcher.search()`` call.
    pseudocount : float, default 1.0
        Same pseudocount used during the search.

    Returns
    -------
    pd.DataFrame (samples × n_pairs)
        Each column is named "geneA/geneB".
    """
    if not pairs:
        raise ValueError("No ratio pairs provided.")

    all_genes = set()
    for p in pairs:
        all_genes.add(p.gene_a)
        all_genes.add(p.gene_b)

    missing = all_genes - set(expression.columns)
    if missing:
        raise ValueError(
            f"Expression data is missing genes required by ratio pairs: "
            f"{sorted(missing)}"
        )

    result = pd.DataFrame(index=expression.index)
    for p in pairs:
        ratio = (
            (expression[p.gene_a].to_numpy(dtype=float) + pseudocount)
            / (expression[p.gene_b].to_numpy(dtype=float) + pseudocount)
        )
        result[p.name] = ratio

    return result


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def build_ratio_features(
    expression: pd.DataFrame,
    labels: Union[pd.Series, np.ndarray, list],
    reference_group: str,
    baseline_group: str,
    top_k: int = 20,
    min_auc: float = 0.70,
    pseudocount: float = 1.0,
    candidate_genes: Optional[List[str]] = None,
) -> Tuple[RatioSearchResult, pd.DataFrame]:
    """
    One-call wrapper: search for ratio biomarkers and build ratio features.

    Parameters
    ----------
    expression : pd.DataFrame (samples × genes)
    labels : group labels
    reference_group, baseline_group : str
    top_k, min_auc, pseudocount, candidate_genes : see RatioBiomarkerSearcher

    Returns
    -------
    (result, ratio_df) where result is a RatioSearchResult and ratio_df
    is the expression data transformed into ratio features.
    """
    searcher = RatioBiomarkerSearcher(
        top_k=top_k, min_auc=min_auc,
        pseudocount=pseudocount,
        candidate_genes=candidate_genes,
    )
    result = searcher.search(expression, labels, reference_group, baseline_group)

    if result.n_found == 0:
        raise ValueError(
            f"No ratio pairs found with AUC >= {min_auc}. "
            f"Try lowering min_auc."
        )

    # Apply to the full expression (all samples, not just the two groups)
    ratio_df = apply_ratios(expression, result.pairs, pseudocount=pseudocount)

    return result, ratio_df