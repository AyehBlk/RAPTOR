"""
Tests for M4: significance-calibrated consensus score.

Rationale
---------
Pre-M4, ``FeatureSelector.consensus_ranking`` produced ``consensus_score``
via 0.7*rank_score + 0.3*selection_count. A noise gene that got lucky
on all 5 methods could score ~1.0 identical to a real-signal gene with
the same profile: the formula had no significance floor. M4 adds
two-tier multiplicative shrinkage based on per-gene DE p-values:
weight=1.0 if p < alpha else 0.5. Kolde et al. 2012 (Robust Rank
Aggregation) is the theoretical precedent; Efron 2004 (empirical-Bayes
local FDR) is the shrinkage precedent. M4 uses a fixed 0.5 factor to
keep the interface at one knob (alpha).

Tests are organized into four groups per the M4 scoping document §10:

    10.1 compute_per_gene_de unit tests (4)
    10.2 apply_significance_calibration unit tests (5)
    10.3 Integration with discover_biomarkers (3)
    10.4 Interaction with existing modules (3)
"""

from __future__ import annotations

import warnings
from typing import List

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery import (
    apply_significance_calibration,
    compute_per_gene_de,
)
from raptor.biomarker_discovery.core import FeatureSelector


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def signal_fixture():
    """60 samples x 500 genes; first 10 carry 2.5 sigma signal."""
    rng = np.random.default_rng(42)
    n_samples, n_features, n_signal = 60, 500, 10
    y = np.array(['disease'] * 30 + ['healthy'] * 30)
    X = rng.standard_normal((n_samples, n_features))
    X[:30, :n_signal] += 2.5
    gene_names = [f'G{i:04d}' for i in range(n_features)]
    X_df = pd.DataFrame(X, columns=gene_names)
    return X_df, y, gene_names[:n_signal]


@pytest.fixture
def noise_fixture():
    """60 samples x 1000 genes; labels shuffled -> no signal."""
    rng = np.random.default_rng(0)
    n_samples, n_features = 60, 1000
    X = rng.standard_normal((n_samples, n_features))
    y = np.array(['a'] * 30 + ['b'] * 30)
    np.random.default_rng(0).shuffle(y)
    gene_names = [f'N{i:04d}' for i in range(n_features)]
    X_df = pd.DataFrame(X, columns=gene_names)
    return X_df, y


@pytest.fixture
def small_signal_fixture():
    """Small fixture for fast integration tests via discover_biomarkers.

    Returns (counts_df, metadata_df) in the shapes discover_biomarkers
    expects: counts = genes x samples (non-negative int), metadata
    indexed by sample_id with a 'condition' column.
    """
    rng = np.random.default_rng(42)
    n_samples, n_genes, n_signal = 40, 200, 10
    y_label = ['disease'] * 20 + ['healthy'] * 20

    X = rng.standard_normal((n_samples, n_genes)) * 0.8
    X[:20, :n_signal] += 2.5
    sample_ids = [f'S{i:03d}' for i in range(n_samples)]
    gene_ids = [f'GENE_{i:04d}' for i in range(n_genes)]

    # counts matrix as RAPTOR expects: genes x samples, non-negative
    counts = pd.DataFrame(
        np.clip(X.T * 100 + 500, 0, None).astype(int),
        index=gene_ids,
        columns=sample_ids,
    )
    metadata = pd.DataFrame(
        {'sample_id': sample_ids, 'condition': y_label}
    ).set_index('sample_id')
    return counts, metadata, gene_ids[:n_signal]


@pytest.fixture
def small_ranked_genes():
    """Build a synthetic ranked_genes DataFrame for unit-testing the
    calibration function without invoking the whole consensus_ranking
    pipeline."""
    return pd.DataFrame(
        {
            'consensus_rank': [1, 2, 3, 4, 5],
            'consensus_score': [0.95, 0.90, 0.80, 0.70, 0.60],
            'n_methods_selected': [3, 3, 2, 2, 1],
        },
        index=pd.Index(['G001', 'G002', 'G003', 'G004', 'G005'], name='gene_id'),
    )


# =============================================================================
# 10.1 compute_per_gene_de unit tests
# =============================================================================


class TestComputePerGeneDE:
    """Per-gene DE utility, the shared substrate for M4 and
    direction_patterns."""

    def test_welch_recovers_signal_genes(self, signal_fixture):
        """On a fixture with known signal in the first 10 genes,
        Welch's t-test should flag them with p_value < 1e-6."""
        X, y, signal_genes = signal_fixture
        out = compute_per_gene_de(
            X, y, reference_group='disease', baseline_group='healthy',
            test='welch',
        )
        assert out.shape == (500, 8)
        assert (out.loc[signal_genes, 'p_value'] < 1e-6).all()
        # Non-signal genes have a uniform distribution of p-values
        # with most far above the signal-gene range
        nonsig_p = out.loc[~out.index.isin(signal_genes), 'p_value']
        assert nonsig_p.median() > 0.1  # most non-signal genes non-significant

    def test_noise_produces_uniform_pvalues(self, noise_fixture):
        """On random-label data, the fraction of genes with p < 0.05
        should match the expected 5% uniform-under-null rate."""
        X, y = noise_fixture
        out = compute_per_gene_de(
            X, y, reference_group='a', baseline_group='b', test='welch',
        )
        pct_sig = float((out['p_value'] < 0.05).mean())
        # Allow wide 2-9% band for sampling variance at n=1000 genes
        assert 0.02 <= pct_sig <= 0.09, (
            f"Expected ~5% significant under null, got {pct_sig:.3f}. "
            f"If outside this band, the test statistic or filtering has "
            f"drifted."
        )

    def test_mann_whitney_agrees_on_signal(self, signal_fixture):
        """Mann-Whitney should also flag the signal genes, though with
        slightly different exact p-values than Welch."""
        X, y, signal_genes = signal_fixture
        out = compute_per_gene_de(
            X, y, reference_group='disease', baseline_group='healthy',
            test='mann_whitney',
        )
        # Every signal gene p < 1e-3 (Mann-Whitney more conservative)
        assert (out.loc[signal_genes, 'p_value'] < 1e-3).all()

    def test_skipped_genes_get_sensible_defaults(self, signal_fixture):
        """Constant-value or near-dead genes should be flagged
        skipped=True with p_value=1.0, log2fc=0.0."""
        X, y, _ = signal_fixture
        X = X.copy()
        X.iloc[:, 0] = 0.0          # constant zero
        X.iloc[:2, 1] = 0.5         # 2 non-zero out of 60 -> below min_nonzero=3
        X.iloc[2:, 1] = 0.0
        out = compute_per_gene_de(
            X, y, reference_group='disease', baseline_group='healthy',
        )
        # Constant-zero gene
        assert out.iloc[0]['skipped'] == True
        assert out.iloc[0]['p_value'] == 1.0
        assert out.iloc[0]['log2fc'] == 0.0
        # Near-dead gene
        assert out.iloc[1]['skipped'] == True
        assert out.iloc[1]['p_value'] == 1.0


# =============================================================================
# 10.2 apply_significance_calibration unit tests
# =============================================================================


class TestApplySignificanceCalibration:
    """The calibration post-processor."""

    def test_adds_expected_columns(self, signal_fixture, small_ranked_genes):
        """Output has p_value, weight, consensus_score_calibrated."""
        X, y, _ = signal_fixture
        # Align ranked index to X columns for this test
        ranked = small_ranked_genes.rename(index={
            'G001': 'G0000', 'G002': 'G0001', 'G003': 'G0002',
            'G004': 'G0003', 'G005': 'G0004',
        })
        ranked.index.name = 'gene_id'
        calibrated = apply_significance_calibration(
            ranked_genes=ranked, X=X, y=y,
            reference_group='disease', baseline_group='healthy',
        )
        assert 'p_value' in calibrated.columns
        assert 'weight' in calibrated.columns
        assert 'consensus_score_calibrated' in calibrated.columns
        # Original column preserved unchanged
        assert 'consensus_score' in calibrated.columns

    def test_significant_genes_get_full_weight(self, signal_fixture):
        """Signal genes with p < 0.05 get weight=1.0 and
        consensus_score_calibrated == consensus_score."""
        X, y, signal_genes = signal_fixture
        ranked = pd.DataFrame(
            {'consensus_rank': range(1, len(signal_genes) + 1),
             'consensus_score': [0.9] * len(signal_genes),
             'n_methods_selected': [3] * len(signal_genes)},
            index=pd.Index(signal_genes, name='gene_id'),
        )
        calibrated = apply_significance_calibration(
            ranked_genes=ranked, X=X, y=y,
            reference_group='disease', baseline_group='healthy',
        )
        assert (calibrated['weight'] == 1.0).all()
        # consensus_score_calibrated == consensus_score for these
        assert np.allclose(
            calibrated['consensus_score_calibrated'],
            calibrated['consensus_score'],
        )

    def test_nonsignificant_genes_get_half_weight(self, noise_fixture):
        """On noise data, most genes get weight=0.5 and
        consensus_score_calibrated == 0.5 * consensus_score."""
        X, y = noise_fixture
        # Use first 50 gene names
        gene_names = list(X.columns[:50])
        ranked = pd.DataFrame(
            {'consensus_rank': range(1, 51),
             'consensus_score': np.linspace(1.0, 0.5, 50),
             'n_methods_selected': [2] * 50},
            index=pd.Index(gene_names, name='gene_id'),
        )
        calibrated = apply_significance_calibration(
            ranked_genes=ranked, X=X, y=y,
            reference_group='a', baseline_group='b',
        )
        # Most noise genes have weight=0.5 (allow ~5% with weight=1.0 by chance)
        n_half = int((calibrated['weight'] == 0.5).sum())
        assert n_half >= 40, (
            f"Expected at least 40/50 noise genes down-weighted; got {n_half}"
        )
        # Their calibrated score is exactly 0.5 * raw
        half_mask = calibrated['weight'] == 0.5
        np.testing.assert_allclose(
            calibrated.loc[half_mask, 'consensus_score_calibrated'],
            calibrated.loc[half_mask, 'consensus_score'] * 0.5,
        )

    def test_consensus_rank_recomputed_on_calibrated(self):
        """A gene with high raw score but p > alpha should be outranked
        by a gene with mid raw score but p < alpha."""
        rng = np.random.default_rng(42)
        n = 30
        # gene A: pure noise, y-independent
        # gene B: carries strong signal
        noise = rng.standard_normal(n)
        y = np.array([0]*(n//2) + [1]*(n//2))
        signal = rng.standard_normal(n) + 3.0 * y

        X = pd.DataFrame({'A': noise, 'B': signal})
        # Hand-build ranked_genes: A ranks higher on raw consensus_score
        ranked = pd.DataFrame(
            {'consensus_rank': [1, 2],
             'consensus_score': [0.99, 0.60],     # A > B raw
             'n_methods_selected': [5, 3]},
            index=pd.Index(['A', 'B'], name='gene_id'),
        )
        calibrated = apply_significance_calibration(
            ranked_genes=ranked, X=X, y=y,
        )
        # A has no signal: p near 1, weight=0.5 -> calibrated=0.495
        # B has strong signal: p near 0, weight=1.0 -> calibrated=0.60
        # So after calibration B should outrank A
        assert calibrated.loc['B', 'consensus_rank'] == 1
        assert calibrated.loc['A', 'consensus_rank'] == 2
        assert calibrated.loc['B', 'consensus_score_calibrated'] > \
               calibrated.loc['A', 'consensus_score_calibrated']

    def test_noop_when_weight_nonsignificant_equals_one(self, signal_fixture):
        """weight_nonsignificant=1.0 leaves all scores unchanged."""
        X, y, _ = signal_fixture
        gene_names = list(X.columns[:20])
        ranked = pd.DataFrame(
            {'consensus_rank': range(1, 21),
             'consensus_score': np.linspace(1.0, 0.5, 20),
             'n_methods_selected': [2] * 20},
            index=pd.Index(gene_names, name='gene_id'),
        )
        calibrated = apply_significance_calibration(
            ranked_genes=ranked, X=X, y=y,
            reference_group='disease', baseline_group='healthy',
            weight_nonsignificant=1.0,
        )
        # All weights == 1.0, so calibrated == raw
        assert (calibrated['weight'] == 1.0).all()
        np.testing.assert_allclose(
            calibrated['consensus_score_calibrated'],
            calibrated['consensus_score'],
        )

    def test_skips_on_non_binary_labels(self, small_ranked_genes):
        """With 3-class or continuous y, function returns unchanged
        with a logged warning."""
        y_3class = np.array([0, 1, 2] * 20)
        X = pd.DataFrame(
            np.random.default_rng(0).standard_normal((60, 5)),
            columns=['G001', 'G002', 'G003', 'G004', 'G005'],
        )
        out = apply_significance_calibration(
            ranked_genes=small_ranked_genes, X=X, y=y_3class,
        )
        # Returned unchanged
        assert 'p_value' not in out.columns
        pd.testing.assert_frame_equal(out, small_ranked_genes)


# =============================================================================
# 10.3 Integration with discover_biomarkers
# =============================================================================


class TestDiscoverBiomarkersIntegration:
    """End-to-end: M4 produces calibrated ranked_genes via the public API."""

    def test_m4_columns_present_on_ranked_genes(self, small_signal_fixture):
        """discover_biomarkers with calibrate_consensus=True produces
        ranked_genes with the three M4 columns."""
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata, _ = small_signal_fixture

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m4_test_integ',
                random_state=42, verbose=False,
                calibrate_consensus=True, alpha=0.05,
            )
        assert 'p_value' in result.ranked_genes.columns
        assert 'weight' in result.ranked_genes.columns
        assert 'consensus_score_calibrated' in result.ranked_genes.columns
        # Parameters recorded
        assert result.parameters['calibrate_consensus'] == True
        assert result.parameters['alpha'] == 0.05

    def test_calibrate_consensus_false_disables_m4(self, small_signal_fixture):
        """calibrate_consensus=False: no M4 columns on ranked_genes."""
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata, _ = small_signal_fixture

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m4_test_integ_off',
                random_state=42, verbose=False,
                calibrate_consensus=False,
            )
        assert 'p_value' not in result.ranked_genes.columns
        assert 'weight' not in result.ranked_genes.columns
        assert 'consensus_score_calibrated' not in result.ranked_genes.columns
        assert result.parameters['calibrate_consensus'] == False

    def test_m4_preserves_m6_pipeline_cv(self, small_signal_fixture):
        """M4 is reporting-layer only; M6's panel_stability field must
        still be populated, OOF AUC should still be sensible."""
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata, _ = small_signal_fixture

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m4_test_m6_preserved',
                random_state=42, verbose=False,
                calibrate_consensus=True,
            )
        # M6 diagnostic still present
        assert result.panel_stability is not None
        assert len(result.panel_stability.per_fold_panels) == 5
        # OOF AUC on signal fixture should be high
        best_clf = result.classification_results[result.best_classifier]
        assert best_clf.oof_true is not None
        from sklearn.metrics import roc_auc_score
        pooled_auc = roc_auc_score(best_clf.oof_true, best_clf.oof_prob)
        assert pooled_auc >= 0.80


# =============================================================================
# 10.4 Interaction with existing modules
# =============================================================================


class TestInteractionWithExistingModules:
    """M4's prerequisite refactor must not regress direction_patterns;
    M4 must handle LOOCV path and edge cases cleanly."""

    def test_build_direction_pattern_still_works(self, signal_fixture):
        """After the refactor to delegate DE to compute_per_gene_de,
        build_direction_pattern returns the same kind of result."""
        from raptor.biomarker_discovery import build_direction_pattern
        X, y, signal_genes = signal_fixture
        dp = build_direction_pattern(
            expression=X, labels=y,
            reference_group='disease', baseline_group='healthy',
            p_threshold=0.05, fc_threshold=1.0,
        )
        # All 10 signal genes should survive p<0.05 AND |log2fc| > 1.0
        assert dp.n_genes >= 10
        for g in signal_genes:
            assert g in dp.gene_directions
            # Signal was added to disease samples, so direction = UP
            assert dp.gene_directions[g] == 'UP'

    def test_m4_runs_inside_pipeline_without_error(self, small_signal_fixture):
        """End-to-end check that M4 does not raise even under the full
        M6 pipeline-CV + M1 OOF + enhanced-analysis chain."""
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata, _ = small_signal_fixture

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m4_test_pipeline',
                random_state=42, verbose=False,
                calibrate_consensus=True, alpha=0.05,
            )
        # Basic sanity: result structure intact
        assert result.panel is not None
        assert len(result.panel) == 5
        assert 'consensus_score_calibrated' in result.ranked_genes.columns

    def test_m4_alpha_parameter_takes_effect(self, small_signal_fixture):
        """Lowering alpha below 0.05 should reduce the number of genes
        with weight=1.0 (fewer pass the stricter threshold)."""
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata, _ = small_signal_fixture

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r_loose = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m4_test_alpha_loose',
                random_state=42, verbose=False,
                calibrate_consensus=True, alpha=0.05,
            )
            r_strict = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m4_test_alpha_strict',
                random_state=42, verbose=False,
                calibrate_consensus=True, alpha=1e-10,
            )
        n_loose = int((r_loose.ranked_genes['weight'] == 1.0).sum())
        n_strict = int((r_strict.ranked_genes['weight'] == 1.0).sum())
        # strict alpha = 1e-10 should cut the pass-rate dramatically
        assert n_strict <= n_loose, (
            f"Stricter alpha should gate fewer genes: loose={n_loose}, "
            f"strict={n_strict}"
        )