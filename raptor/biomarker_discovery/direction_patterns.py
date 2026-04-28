"""
Direction Pattern Module
========================

A biomarker signature is not just "these genes matter" -- it's "these genes
go UP and these go DOWN" in disease. DirectionPattern captures that
structure: per-gene direction (UP/DOWN), magnitude (log2 fold change),
and confidence (-log10 p-value).

Once a pattern is learned from a discovery cohort it can be applied to:

1. A single patient -- how well does this patient's expression match the
   disease pattern? (concordance score in [-1, +1])
2. A second cohort's pattern -- do the directions replicate across studies?
   (cross-cohort direction check)

This is the workhorse for translational and multi-cohort biomarker work,
where direction consistency matters more than exact fold-change values.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Dict, List, Union

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DIRECTION_UP = "UP"
DIRECTION_DOWN = "DOWN"
VALID_DIRECTIONS = (DIRECTION_UP, DIRECTION_DOWN)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

@dataclass
class DirectionPattern:
    """
    An up/down expression pattern learned from a discovery cohort.

    Parameters
    ----------
    gene_directions : Dict[str, str]
        Gene symbol -> 'UP' or 'DOWN'.
    fold_changes : Dict[str, float]
        Gene symbol -> log2 fold change (reference - baseline).
    confidence : Dict[str, float]
        Gene symbol -> -log10(p-value). Higher means more confident.
    reference_group : str
        Name of the group of interest (e.g. 'disease', 'treated').
    baseline_group : str
        Name of the control group (e.g. 'healthy', 'untreated').
    baseline_means : Dict[str, float]
        Per-gene mean expression in the baseline group, used later for
        z-scoring new patients in :meth:`concordance`.
    baseline_stds : Dict[str, float]
        Per-gene std (ddof=1) in the baseline group.
    n_samples : int
        Total samples (reference + baseline) used to derive the pattern.
    test_used : str, default 't-test'
        The statistical test used to derive confidence values.
    """

    gene_directions: Dict[str, str]
    fold_changes: Dict[str, float]
    confidence: Dict[str, float]
    reference_group: str
    baseline_group: str
    baseline_means: Dict[str, float]
    baseline_stds: Dict[str, float]
    n_samples: int
    test_used: str = "t-test"

    # -- validation ----------------------------------------------------------

    def __post_init__(self) -> None:
        gene_set = set(self.gene_directions.keys())
        for name, d in (
            ("fold_changes", self.fold_changes),
            ("confidence", self.confidence),
            ("baseline_means", self.baseline_means),
            ("baseline_stds", self.baseline_stds),
        ):
            if set(d.keys()) != gene_set:
                raise ValueError(
                    f"Gene set mismatch: gene_directions has {len(gene_set)} "
                    f"gene(s) but {name} has {len(d)}. Every dict must cover "
                    f"exactly the same genes."
                )
        for gene, direction in self.gene_directions.items():
            if direction not in VALID_DIRECTIONS:
                raise ValueError(
                    f"Invalid direction {direction!r} for gene {gene!r}. "
                    f"Must be one of: {VALID_DIRECTIONS}."
                )

    # -- cheap properties ----------------------------------------------------

    @property
    def genes(self) -> List[str]:
        return list(self.gene_directions.keys())

    @property
    def n_genes(self) -> int:
        return len(self.gene_directions)

    @property
    def n_up(self) -> int:
        return sum(1 for d in self.gene_directions.values() if d == DIRECTION_UP)

    @property
    def n_down(self) -> int:
        return sum(1 for d in self.gene_directions.values() if d == DIRECTION_DOWN)

    # -- main analyses -------------------------------------------------------

    def concordance(
        self,
        patient_expression: Union[pd.Series, Dict[str, float], pd.DataFrame],
        weighted: bool = True,
    ) -> Union[float, pd.Series]:
        """
        Score how well a patient's expression matches this pattern.

        For each gene in the pattern:
          - z = (patient_value - baseline_mean) / baseline_std
          - expected sign: +1 if UP, -1 if DOWN
          - per-gene agreement = tanh(z) * expected_sign, squashed to [-1, +1]
            (tanh makes it robust to one outlier gene dominating the score)

        Then take the (confidence-weighted) mean across genes. Final score
        is in [-1, +1]:

          +1  pattern is fully realized in this patient
           0  no relationship to the pattern
          -1  opposite pattern

        Parameters
        ----------
        patient_expression : pd.Series, dict, or pd.DataFrame
            Single patient as Series or dict (gene -> value).
            Multiple patients as DataFrame (samples x genes).
        weighted : bool, default True
            Weight each gene's contribution by its confidence (-log10 p).
            If False, all pattern genes contribute equally.

        Returns
        -------
        float or pd.Series
            Single float if input was Series/dict; Series (one score per
            sample) if input was a DataFrame.
        """
        # -- normalize input to DataFrame for uniform handling --
        if isinstance(patient_expression, dict):
            patient_expression = pd.Series(patient_expression)
        is_single = isinstance(patient_expression, pd.Series)
        if is_single:
            patient_expression = patient_expression.to_frame().T

        # -- find gene overlap with the pattern --
        common = [g for g in self.genes if g in patient_expression.columns]
        if not common:
            raise ValueError(
                "No genes in patient data match this pattern. "
                f"Pattern expects (first 5): {self.genes[:5]}"
            )

        missing = set(self.genes) - set(common)
        if missing:
            warnings.warn(
                f"{len(missing)} of {self.n_genes} pattern genes are missing "
                f"from patient data; concordance computed on the remaining "
                f"{len(common)}.",
                UserWarning,
                stacklevel=2,
            )

        # -- vectorized computation --
        means = np.array([self.baseline_means[g] for g in common])
        stds = np.array([self.baseline_stds[g] for g in common])
        stds = np.where(stds > 1e-10, stds, 1e-10)  # guard against zero var

        expected_signs = np.array(
            [1.0 if self.gene_directions[g] == DIRECTION_UP else -1.0
             for g in common]
        )
        if weighted:
            weights = np.array([self.confidence[g] for g in common], dtype=float)
        else:
            weights = np.ones(len(common), dtype=float)

        # Fall back to uniform weights if all confidences collapse to 0.
        total_weight = weights.sum()
        if total_weight < 1e-10:
            weights = np.ones(len(common), dtype=float)
            total_weight = float(len(common))

        patient_vals = patient_expression[common].to_numpy(dtype=float)
        z = (patient_vals - means) / stds
        per_gene = np.tanh(z) * expected_signs          # samples x genes
        score = (per_gene * weights).sum(axis=1) / total_weight

        result = pd.Series(score, index=patient_expression.index,
                           name="concordance")
        return float(result.iloc[0]) if is_single else result

    def cross_cohort_check(
        self,
        other: "DirectionPattern",
    ) -> Dict[str, Any]:
        """
        Check direction agreement between this pattern and another's.

        Returns a report dict with:
          - n_common, n_agree, n_disagree
          - agreement_fraction   (n_agree / n_common)
          - binomial_p           (two-sided test vs chance p=0.5)
          - agreements           list of genes that agree
          - disagreements        list of (gene, this_dir, other_dir)
          - only_in_this, only_in_other
        """
        if not isinstance(other, DirectionPattern):
            raise TypeError(
                "cross_cohort_check requires another DirectionPattern, "
                f"got {type(other).__name__}."
            )

        this_genes = set(self.genes)
        other_genes = set(other.genes)
        common = sorted(this_genes & other_genes)
        only_in_this = sorted(this_genes - other_genes)
        only_in_other = sorted(other_genes - this_genes)

        agreements: List[str] = []
        disagreements: List[tuple] = []
        for gene in common:
            if self.gene_directions[gene] == other.gene_directions[gene]:
                agreements.append(gene)
            else:
                disagreements.append(
                    (gene, self.gene_directions[gene], other.gene_directions[gene])
                )

        n_common = len(common)
        n_agree = len(agreements)
        if n_common > 0:
            agreement_fraction = n_agree / n_common
            binom_p = float(stats.binomtest(n_agree, n_common, p=0.5).pvalue)
        else:
            agreement_fraction = float("nan")
            binom_p = float("nan")

        return {
            "n_common": n_common,
            "n_agree": n_agree,
            "n_disagree": n_common - n_agree,
            "agreement_fraction": agreement_fraction,
            "binomial_p": binom_p,
            "agreements": agreements,
            "disagreements": disagreements,
            "only_in_this": only_in_this,
            "only_in_other": only_in_other,
        }

    # -- inspection ----------------------------------------------------------

    def to_dataframe(self) -> pd.DataFrame:
        """Return the pattern as a gene-indexed DataFrame for inspection."""
        return pd.DataFrame({
            "direction": pd.Series(self.gene_directions),
            "log2FC": pd.Series(self.fold_changes),
            "neg_log10_p": pd.Series(self.confidence),
            "baseline_mean": pd.Series(self.baseline_means),
            "baseline_std": pd.Series(self.baseline_stds),
        })


# ---------------------------------------------------------------------------
# Factory function
# ---------------------------------------------------------------------------

def build_direction_pattern(
    expression: pd.DataFrame,
    labels: Union[pd.Series, np.ndarray, List],
    reference_group: str,
    baseline_group: str,
    p_threshold: float = 0.05,
    fc_threshold: float = 1.0,       # log2
    test: str = "t-test",
    min_nonzero: int = 3,
) -> DirectionPattern:
    """
    Learn a direction pattern from a discovery cohort.

    Parameters
    ----------
    expression : pd.DataFrame (samples x genes)
        Expression matrix. Typically log2-transformed normalized counts,
        because fold changes are computed as mean(ref) - mean(baseline).
    labels : Series, array, or list
        Group label per sample. Length must equal expression.shape[0].
    reference_group : str
        The group of interest (e.g. 'disease'). UP genes are elevated here.
    baseline_group : str
        The control group. Baseline means and stds are computed from this
        group for later patient z-scoring in DirectionPattern.concordance.
    p_threshold : float, default 0.05
        Per-gene p-value cutoff. Genes with p > threshold are dropped.
    fc_threshold : float, default 1.0
        Absolute log2 fold-change cutoff. Genes below this are dropped.
    test : {'t-test', 'mann-whitney'}, default 't-test'
    min_nonzero : int, default 3
        Genes with fewer nonzero values across both groups are skipped.

    Implementation notes
    --------------------
    As of M4 (2026-04), the per-gene DE computation is delegated to
    ``raptor.biomarker_discovery.univariate_de.compute_per_gene_de``.
    This function then applies the p_threshold / fc_threshold filter
    and assembles the surviving genes into a DirectionPattern. Both
    this function and M4's ``apply_significance_calibration`` consume
    the same underlying computation, so their p-values always agree
    to full float precision on the same (expression, labels) input.

    The test name mapping is:
        't-test'        -> compute_per_gene_de(test='welch')
        'mann-whitney'  -> compute_per_gene_de(test='mann_whitney')
    """
    from raptor.biomarker_discovery.univariate_de import compute_per_gene_de

    if test not in ("t-test", "mann-whitney"):
        raise ValueError(
            f"test must be 't-test' or 'mann-whitney', got {test!r}."
        )

    # Translate direction_patterns test-name convention to univariate_de convention.
    test_de = "welch" if test == "t-test" else "mann_whitney"

    # Delegate per-gene DE to the shared utility. This returns a row
    # per gene with p_value, neg_log10p, log2fc, mean_ref, mean_base,
    # n_ref, n_base, and skipped columns. Skipped genes have p=1.0 so
    # they're naturally dropped by the p_threshold filter below; the
    # filter logic therefore doesn't need to special-case them.
    de_table = compute_per_gene_de(
        expression=expression,
        labels=labels,
        reference_group=reference_group,
        baseline_group=baseline_group,
        test=test_de,
        min_nonzero=min_nonzero,
    )

    # Need per-gene baseline means/stds for the concordance computation
    # later. compute_per_gene_de gives us mean_base but not std_base,
    # so we compute std_base separately. This stays vectorized.
    labels_s = pd.Series(labels).reset_index(drop=True)
    expression_rs = expression.reset_index(drop=True)
    base_mask = (labels_s == baseline_group).to_numpy()
    X_base = expression_rs.to_numpy(dtype=float)[base_mask, :]
    std_base_all = np.std(X_base, axis=0, ddof=1)

    # Apply the DirectionPattern filter: p < threshold AND |log2fc| > threshold.
    # Skipped genes (p=1.0, log2fc=0.0) are dropped automatically by these conditions.
    keep = (
        (de_table["p_value"] < p_threshold) &
        (de_table["log2fc"].abs() >= fc_threshold)
    )
    kept_de = de_table.loc[keep]

    if kept_de.empty:
        raise ValueError(
            f"No genes passed thresholds (p < {p_threshold}, "
            f"|log2FC| > {fc_threshold}). Loosen thresholds or check input."
        )

    # Build the per-gene dicts expected by the DirectionPattern dataclass.
    # gene order is preserved from de_table (which preserves expression.columns).
    gene_directions: Dict[str, str] = {}
    fold_changes: Dict[str, float] = {}
    confidence: Dict[str, float] = {}
    baseline_means: Dict[str, float] = {}
    baseline_stds: Dict[str, float] = {}

    # Map gene_id -> std_base via positional index (de_table is in expression.columns order)
    gene_col_idx = {g: i for i, g in enumerate(expression.columns)}

    for gene, row in kept_de.iterrows():
        log2fc = float(row["log2fc"])
        gene_directions[gene] = DIRECTION_UP if log2fc > 0 else DIRECTION_DOWN
        fold_changes[gene] = log2fc
        confidence[gene] = float(row["neg_log10p"])
        baseline_means[gene] = float(row["mean_base"])
        baseline_stds[gene] = float(std_base_all[gene_col_idx[gene]])

    # n_samples = n_ref + n_base is constant across rows in de_table.
    n_total = int(kept_de["n_ref"].iloc[0] + kept_de["n_base"].iloc[0])

    return DirectionPattern(
        gene_directions=gene_directions,
        fold_changes=fold_changes,
        confidence=confidence,
        reference_group=reference_group,
        baseline_group=baseline_group,
        baseline_means=baseline_means,
        baseline_stds=baseline_stds,
        n_samples=n_total,
        test_used=test,
    )