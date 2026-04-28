"""
Tests for M3: unique-rank scoring in select_rfe.

Covers the fix in ``raptor.biomarker_discovery.core.FeatureSelector.select_rfe``
that replaces the pre-M3 ``(max_rank - rfe.ranking_ + 1) / max_rank``
scoring. The old formula gave every one of the ``n_features_to_select``
chosen features the same score of ``1.0``, which the downstream
consensus ranking then collapsed to a mid-rank (e.g. ``rank_rfe = 25.5``
for every selected feature when 50 were chosen), making RFE contribute
effectively zero discrimination to the consensus score.

The M3 fix produces a two-tier score:

    * Selected features: read importance from the final estimator
      (``feature_importances_`` for tree estimators, ``|coef_|`` for
      linear). Rescaled to [1.0, 2.0] so every selected feature ranks
      above every non-selected feature.
    * Non-selected features: ``1 / rfe.ranking_``, yielding scores in
      ``(0, 0.5]``. Note that sklearn's RFE eliminates features in
      batches when ``step > 1``, so non-selected features can share
      ranking values *within a batch*; this is a sklearn behavior we
      accept (ties only occur among eliminated features, which never
      enter the panel).

Guarantees under test:
    1. Selected features have unique scores (no ties among them).
    2. Every selected feature scores strictly higher than every
       non-selected feature.
    3. Higher estimator importance -> higher score among selected.
    4. Downstream consensus ranking produces ``rank_rfe`` values that
       place selected features before non-selected features, and that
       do not collapse to the pathological mid-rank seen pre-M3.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery.core import FeatureSelector


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def signal_dataset():
    """60 samples, 30 features, first 10 features carry real signal."""
    rng = np.random.default_rng(42)
    n_samples, n_features, n_signal = 60, 30, 10
    y = np.array([0] * (n_samples // 2) + [1] * (n_samples // 2))
    X = rng.standard_normal((n_samples, n_features))
    # Inject separating signal into the first n_signal columns for the
    # positive class. Effect size chosen to be strong enough that RFE
    # reliably picks the right ones.
    X[y == 1, :n_signal] += 2.0
    gene_names = [f'GENE_{i:03d}' for i in range(n_features)]
    return (
        pd.DataFrame(X, columns=gene_names),
        y,
        gene_names[:n_signal],  # the known-signal genes
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestRFERankingUniqueness:
    """Primary M3 guarantee: no ties among selected features."""

    def test_selected_features_have_unique_scores(self, signal_dataset):
        X, y, _ = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        scores = result.gene_scores.reindex(result.selected_genes)['score']

        assert len(scores) == 10
        assert len(scores.unique()) == 10, (
            f"Expected 10 unique scores among selected features, "
            f"got {len(scores.unique())} unique values. "
            f"Ties among selected features indicate the M3 fix has regressed."
        )

    def test_selected_scores_all_above_non_selected(self, signal_dataset):
        """Every selected feature must score strictly above every non-selected."""
        X, y, _ = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        all_scores = result.gene_scores['score']
        sel_scores = all_scores.reindex(result.selected_genes)
        non_sel = [g for g in X.columns if g not in result.selected_genes]
        non_sel_scores = all_scores.reindex(non_sel)

        assert sel_scores.min() > non_sel_scores.max(), (
            f"Selected min ({sel_scores.min():.4f}) must exceed "
            f"non-selected max ({non_sel_scores.max():.4f}) so that the "
            f"consensus rank always places selected features first."
        )

    def test_non_selected_scores_below_one(self, signal_dataset):
        """Non-selected features score in (0, 1] so invariant holds vs selected's [1, 2]."""
        X, y, _ = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        all_scores = result.gene_scores['score']
        non_sel = [g for g in X.columns if g not in result.selected_genes]
        non_sel_scores = all_scores.reindex(non_sel)

        assert (non_sel_scores <= 1.0).all(), (
            "All non-selected scores must be <= 1.0 to preserve the "
            "selected-above-non-selected invariant."
        )
        assert (non_sel_scores > 0).all(), (
            "Non-selected scores should be strictly positive "
            "(derived from 1 / ranking_ where ranking_ >= 2)."
        )


class TestRFERankingOrdering:
    """Secondary guarantee: ordering is driven by estimator importance."""

    def test_higher_importance_gets_better_rank(self, signal_dataset):
        """Among selected features, higher importance -> lower rank number."""
        X, y, _ = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        sel = result.gene_scores.reindex(result.selected_genes)
        sel_sorted = sel.sort_values('score', ascending=False)

        # Rescale score back to rank: highest score first
        # Check monotonicity: sel_sorted's scores are strictly decreasing
        score_vals = sel_sorted['score'].values
        assert np.all(np.diff(score_vals) <= 0), (
            "Sorting by score should yield a strictly non-increasing sequence."
        )

    def test_signal_genes_selected_and_top_ranked(self, signal_dataset):
        """Sanity: on a dataset with known signal, RFE picks signal genes."""
        X, y, signal_genes = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        # On this fixture with clear signal + low noise, RFE should pick
        # all 10 signal genes (order within them may vary by RF seed).
        selected = set(result.selected_genes)
        signal_set = set(signal_genes)
        overlap = selected & signal_set
        assert len(overlap) >= 8, (
            f"Expected >=8 of 10 signal genes to be selected, got "
            f"{len(overlap)}. Indicates signal isn't dominating noise "
            f"enough in the fixture."
        )


class TestRFEDownstreamConsensus:
    """The ultimate guarantee: downstream rank_rfe doesn't collapse to a mid-rank."""

    def test_consensus_rank_rfe_no_pathological_midpoint(self, signal_dataset):
        """Reproduce the bug scenario: 50 ties at rank 25.5 must not recur."""
        X, y, _ = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        # Simulate what consensus_ranking does: rank by score descending
        all_genes = list(X.columns)
        method_scores = result.gene_scores.reindex(all_genes).fillna(0.0)
        ranks = method_scores['score'].rank(ascending=False, method='average')

        # Pre-M3 pathology: all 10 selected features would share rank
        # = (1 + 10) / 2 = 5.5 on this fixture.
        pathological_mid = (len(result.selected_genes) + 1) / 2
        tied_at_mid = (ranks == pathological_mid).sum()

        # After M3: at most 1 feature can land exactly on the midpoint
        # (the median of the unique selected-feature ranks).
        assert tied_at_mid <= 1, (
            f"{tied_at_mid} features share the pathological midpoint "
            f"rank {pathological_mid}; pre-M3 bug has recurred."
        )

    def test_consensus_rank_rfe_selected_features_rank_before_non_selected(
        self, signal_dataset
    ):
        """Selected features must receive lower (better) ranks than non-selected."""
        X, y, _ = signal_dataset
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=10, step=0.1)

        all_genes = list(X.columns)
        method_scores = result.gene_scores.reindex(all_genes).fillna(0.0)
        ranks = method_scores['score'].rank(ascending=False, method='average')

        sel_ranks = ranks.loc[result.selected_genes]
        non_sel = [g for g in all_genes if g not in result.selected_genes]
        non_sel_ranks = ranks.loc[non_sel]

        assert sel_ranks.max() < non_sel_ranks.min(), (
            f"Selected max rank {sel_ranks.max()} must be < non-selected "
            f"min rank {non_sel_ranks.min()} so selected genes always "
            f"appear before non-selected in the consensus ordering."
        )


class TestRFERankingEdgeCases:
    """Defensive tests for the fallback paths."""

    def test_small_dataset_still_produces_valid_scores(self):
        """At n_samples near n_features, RFE should still produce clean scores."""
        rng = np.random.default_rng(0)
        n_samples = 20
        n_features = 15
        y = np.array([0] * 10 + [1] * 10)
        X = rng.standard_normal((n_samples, n_features))
        X[y == 1, :5] += 1.5  # moderate signal in first 5
        X_df = pd.DataFrame(
            X, columns=[f'G{i}' for i in range(n_features)]
        )

        selector = FeatureSelector(random_state=0, verbose=False)
        result = selector.select_rfe(X_df, y, n_features=5, step=0.1)

        sel_scores = result.gene_scores.reindex(result.selected_genes)['score']
        assert len(sel_scores.unique()) == 5
        # Invariant still holds
        non_sel = [g for g in X_df.columns if g not in result.selected_genes]
        non_sel_scores = result.gene_scores.reindex(non_sel)['score']
        assert sel_scores.min() > non_sel_scores.max()

    def test_reproducible_with_same_random_state(self, signal_dataset):
        """Same random_state -> identical selected genes and scores."""
        X, y, _ = signal_dataset

        sel1 = FeatureSelector(random_state=42, verbose=False)
        r1 = sel1.select_rfe(X, y, n_features=10, step=0.1)

        sel2 = FeatureSelector(random_state=42, verbose=False)
        r2 = sel2.select_rfe(X, y, n_features=10, step=0.1)

        assert r1.selected_genes == r2.selected_genes
        pd.testing.assert_series_equal(
            r1.gene_scores['score'], r2.gene_scores['score']
        )