"""
Univariate differential expression utility.

Computes per-gene two-sample tests (Welch's t or Mann-Whitney U) across
all genes of an expression matrix, returning a full per-gene p-value
table including genes that would normally be filtered out.

This module is the shared substrate for:
    * ``direction_patterns.build_direction_pattern`` — filters to
      p < alpha and |log2FC| > threshold, wraps in a DirectionPattern.
    * ``core.apply_significance_calibration`` (M4) — needs p-values
      for ALL genes, including non-significant ones, to down-weight
      their consensus score.
    * ``FeatureSelector.select_univariate_filter`` (M6) — alternative
      implementation of the same computation optimized for the hot
      per-fold pipeline-CV loop (kept separate to avoid cross-module
      import cost during 5x outer CV).

Design choices
--------------
* Welch's t-test is vectorized across genes via ``scipy.stats.ttest_ind``
  with ``axis=0``. A 20k-gene by 60-sample matrix completes in ~10 ms.
* Mann-Whitney U is looped (scipy's implementation is not vectorized
  across an axis). Tolerable cost: ~2 s on the same matrix.
* All-NaN / zero-variance / skipped genes get ``p_value=1.0`` and
  ``log2fc=0.0``. Consumers that filter on significance will therefore
  drop them without special-casing. Consumers that calibrate (M4) will
  treat them as non-significant (weight = 0.5 by default), which is
  the honest behavior for genes that carry no information.
* No FDR correction is applied here. Downstream consumers can apply it
  if they want; M4 explicitly does not, for reasons documented in the
  M4 scoping document.

References
----------
Welch, B. L. (1947). The generalization of "Student's" problem when
    several different population variances are involved. Biometrika,
    34(1/2), 28–35.
Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two
    random variables is stochastically larger than the other. Annals of
    Mathematical Statistics, 18(1), 50–60.
"""

from __future__ import annotations

import logging
from typing import List, Union

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public names
# ---------------------------------------------------------------------------

VALID_TESTS = ("welch", "mann_whitney")
"""Supported two-sample test names accepted by ``compute_per_gene_de``."""


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def compute_per_gene_de(
    expression: pd.DataFrame,
    labels: Union[pd.Series, np.ndarray, List],
    reference_group: str,
    baseline_group: str,
    test: str = "welch",
    min_nonzero: int = 3,
) -> pd.DataFrame:
    """
    Compute per-gene two-sample differential expression across all genes.

    Returns ONE row per gene in ``expression.columns``, including genes
    that fail the minimum-nonzero filter or have zero variance. Those
    "skipped" genes receive ``p_value=1.0`` and ``log2fc=0.0`` so that
    downstream consumers can either filter (direction_patterns style)
    or calibrate (M4 style) without having to special-case missing rows.

    Parameters
    ----------
    expression : pd.DataFrame (samples x genes)
        Expression matrix. Typically log2-transformed normalized counts,
        because fold changes are computed as ``mean(ref) - mean(base)``.
        Row order matches the ``labels`` argument.
    labels : Series, array, or list
        Group label per sample. Length must equal ``expression.shape[0]``.
    reference_group : str
        The group of interest (e.g. 'disease'). UP genes are elevated
        here relative to baseline.
    baseline_group : str
        The control group. log2fc = mean(ref) - mean(base).
    test : {'welch', 'mann_whitney'}, default 'welch'
        'welch'        Welch's t-test, two-sided, unequal variance.
                       Vectorized over all genes (~10 ms for 20k genes).
        'mann_whitney' Mann-Whitney U, two-sided. Looped per gene
                       (~2 s for 20k genes). Use for small-n or when
                       heavy outliers are expected.
    min_nonzero : int, default 3
        Minimum total non-zero values across both groups for a gene to
        be tested. Genes below this are marked ``skipped=True`` with
        ``p_value=1.0``, ``log2fc=0.0``. Prevents spurious tests on
        near-dead genes.

    Returns
    -------
    pd.DataFrame indexed by gene_id with columns:
        p_value     : float in [0, 1]
        neg_log10p  : -log10(max(p, 1e-300))
        log2fc      : mean(ref) - mean(base); 0.0 for skipped genes
        mean_ref    : float
        mean_base   : float
        n_ref       : int, sample count in reference group
        n_base      : int, sample count in baseline group
        skipped     : bool, True when gene failed min_nonzero or was
                      fully constant across both groups

    Raises
    ------
    ValueError
        If test is not recognized, labels/expression lengths mismatch,
        or either group has fewer than 2 samples.
    """
    # --- Validate inputs ---
    if test not in VALID_TESTS:
        raise ValueError(
            f"test must be one of {VALID_TESTS}, got {test!r}."
        )

    labels = pd.Series(labels).reset_index(drop=True)
    expression = expression.reset_index(drop=True)

    if len(labels) != len(expression):
        raise ValueError(
            f"labels length ({len(labels)}) does not match "
            f"expression rows ({len(expression)})."
        )

    groups_present = set(labels.unique())
    for g in (reference_group, baseline_group):
        if g not in groups_present:
            raise ValueError(
                f"Group {g!r} not found in labels. "
                f"Available: {sorted(groups_present)}."
            )

    ref_mask = (labels == reference_group).to_numpy()
    base_mask = (labels == baseline_group).to_numpy()
    n_ref = int(ref_mask.sum())
    n_base = int(base_mask.sum())

    if n_ref < 2 or n_base < 2:
        raise ValueError(
            f"Need at least 2 samples per group; got "
            f"{n_ref} for {reference_group!r} and {n_base} for {baseline_group!r}."
        )

    gene_ids = list(expression.columns)
    n_genes = len(gene_ids)

    X = expression.to_numpy(dtype=float)          # (n_samples, n_genes)
    X_ref = X[ref_mask, :]                        # (n_ref, n_genes)
    X_base = X[base_mask, :]                      # (n_base, n_genes)

    # --- Skipped-gene detection (vectorized) ---
    # Genes with fewer than min_nonzero total non-zero values, or zero
    # variance in BOTH groups simultaneously, are marked skipped.
    nonzero_total = (X_ref != 0).sum(axis=0) + (X_base != 0).sum(axis=0)
    ref_std = np.std(X_ref, axis=0)
    base_std = np.std(X_base, axis=0)
    both_constant = (ref_std == 0) & (base_std == 0)
    skipped_mask = (nonzero_total < min_nonzero) | both_constant

    # --- log2fc: always compute, even for skipped (harmless) ---
    mean_ref = np.mean(X_ref, axis=0)
    mean_base = np.mean(X_base, axis=0)
    log2fc = mean_ref - mean_base
    # For skipped genes we want log2fc=0.0 in the output so consumers
    # can filter on |log2fc| cleanly.
    log2fc_out = np.where(skipped_mask, 0.0, log2fc)

    # --- Two-sample test across all genes ---
    p_values = np.ones(n_genes, dtype=float)

    if test == "welch":
        # Vectorized Welch's t-test across all columns at once.
        # scipy.stats.ttest_ind with axis=0 handles the matrix; nan_policy
        # = 'omit' makes sure stray NaN cells don't propagate.
        tested_mask = ~skipped_mask
        if tested_mask.any():
            with np.errstate(all="ignore"):
                _, p_tested = sp_stats.ttest_ind(
                    X_ref[:, tested_mask],
                    X_base[:, tested_mask],
                    axis=0,
                    equal_var=False,
                    nan_policy="omit",
                )
            p_tested = np.asarray(p_tested, dtype=float).ravel()
            # Any remaining NaN (e.g. zero-variance on one side only) -> 1.0
            p_tested = np.where(np.isnan(p_tested), 1.0, p_tested)
            p_values[tested_mask] = p_tested

    elif test == "mann_whitney":
        # Mann-Whitney is not vectorized over axis in scipy; loop per gene.
        for j in range(n_genes):
            if skipped_mask[j]:
                continue
            try:
                with np.errstate(all="ignore"):
                    _, p = sp_stats.mannwhitneyu(
                        X_ref[:, j], X_base[:, j], alternative="two-sided"
                    )
                p_values[j] = float(p) if np.isfinite(p) else 1.0
            except ValueError:
                # All values tied across both groups -> no test defined
                p_values[j] = 1.0

    # --- Compile output ---
    # Floor p-values to avoid log10(0); 1e-300 is well below any real
    # two-sample p-value and matches direction_patterns.py convention.
    p_floor = np.maximum(p_values, 1e-300)
    neg_log10p = -np.log10(p_floor)

    out = pd.DataFrame(
        {
            "p_value": p_values,
            "neg_log10p": neg_log10p,
            "log2fc": log2fc_out,
            "mean_ref": mean_ref,
            "mean_base": mean_base,
            "n_ref": n_ref,
            "n_base": n_base,
            "skipped": skipped_mask,
        },
        index=pd.Index(gene_ids, name="gene_id"),
    )

    return out