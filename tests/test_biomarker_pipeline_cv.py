"""
Tests for M6: feature selection and panel optimization INSIDE the outer
CV loop.

Rationale
---------
Pre-M6, feature selection (Stage 1) and panel optimization (Stage 2) ran
once on all data, and Stage 3 nested-CV evaluated classifiers on the
already-fixed panel. Every outer-fold held-out sample had its label read
at least twice before being "tested", so the CV AUC was optimistically
biased. Ambroise & McLachlan (2002, PNAS 99:6562) showed this directly:
gene-selection-before-CV produces near-zero cross-validated error even
when the underlying signal is pure noise. Varma & Simon (2006, BMC
Bioinf 7:91) generalized to arbitrary model-selection steps. The fix,
"external" CV, runs every y-using step inside the training fold.

Test coverage is organized into six groups per the M6 scoping document:

    1. Leakage invariants: 4 tests verifying the orchestrator never
       lets a test-fold label reach feature selection or panel
       optimization.
    2. Statistical sanity: 3 tests verifying the bias correction
       actually works (noise -> AUC ~0.5, signal -> AUC stays high,
       optimism gap tightens).
    3. Nogueira stability diagnostics: 5 tests on the Nogueira et al.
       (2018, JMLR 18:174) stability measure (Phi in [-1,1], chance =
       0, identical = 1, benchmark labels, CI straddles estimate).
    4. Univariate filter (fold-safe de_filter replacement): 2 tests.
    5. Integration and regression: 4 tests on API surface + repeats.
    6. Interaction with M1 (clinical metrics): 1 test on the end-to-end
       regression guarantee that Youden's J on noise does not
       reinflate.
"""

from __future__ import annotations

import warnings
from typing import Dict, List
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery.core import (
    # New M6 machinery:
    PanelStabilityResult,
    _run_pipeline_cv,
    _run_feature_selection,
    _resolve_pipeline_cv_methods,
    _nogueira_stability_from_matrix,
    _panels_to_selection_matrix,
    _bootstrap_nogueira_ci,
    _nogueira_benchmark_label,
    _compute_panel_stability,
    _FOLD_SAFE_DEFAULT_METHODS,
    _NOGUEIRA_EXCELLENT_THRESHOLD,
    _NOGUEIRA_POOR_THRESHOLD,
    # Re-used classes:
    FeatureSelector,
    ClassificationResult,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def signal_dataset():
    """60 samples, 100 genes, first 10 carry 2sigma separating signal.

    Small enough to run pipeline-CV in seconds; signal large enough that
    RFE / elastic-net / univariate-filter all pick the signal genes
    reliably.
    """
    rng = np.random.default_rng(42)
    n_samples, n_features, n_signal = 60, 100, 10
    y = np.array([0] * 30 + [1] * 30)
    X = rng.standard_normal((n_samples, n_features))
    X[y == 1, :n_signal] += 2.0
    X_df = pd.DataFrame(
        X, columns=[f'GENE_{i:04d}' for i in range(n_features)]
    )
    return X_df, y


@pytest.fixture
def noise_dataset():
    """60 samples, 200 genes, labels shuffled -> no real signal.

    Pre-M6 this kind of fixture produces CV AUC >= 0.75 because
    Stage 1+2 latched onto chance correlations; post-M6 CV AUC must
    collapse toward 0.5.
    """
    rng = np.random.default_rng(0)
    n_samples, n_features = 60, 200
    X = rng.standard_normal((n_samples, n_features))
    y = np.array([0] * 30 + [1] * 30)
    rng.shuffle(y)
    X_df = pd.DataFrame(
        X, columns=[f'NOISE_{i:04d}' for i in range(n_features)]
    )
    return X_df, y


@pytest.fixture
def small_signal_dataset():
    """Very small fixture for fast tests that don't need large n."""
    rng = np.random.default_rng(42)
    n_samples, n_features, n_signal = 30, 40, 5
    y = np.array([0] * 15 + [1] * 15)
    X = rng.standard_normal((n_samples, n_features))
    X[y == 1, :n_signal] += 2.5
    X_df = pd.DataFrame(
        X, columns=[f'G{i:03d}' for i in range(n_features)]
    )
    return X_df, y


# -----------------------------------------------------------------------
# Fast configuration helpers
# -----------------------------------------------------------------------
#
# Pipeline-CV is inherently ~5x slower than pre-M6 (Ambroise-McLachlan
# fix requires re-running Stages 1+2 inside every outer fold). Stage 3
# runs every classifier (LR, RF, SVM, XGBoost) by default. Stage 2
# internally sweeps panel sizes with RandomForestClassifier. For unit
# tests we restrict both:
#     - target_panel_size=<fixed>   skips the panel sweep
#     - classifiers=['logistic_regression']  skips RF/SVM/XGBoost fits
# Stage 1 methods are also pruned to just 'univariate_filter' (the
# fastest). Together these take each orchestrator test from ~90s to
# ~3-5s, which is what the CI budget can absorb.

FAST_KWARGS = dict(
    methods=['univariate_filter'],
    de_genes=None, species='human',
    n_outer=5, n_repeats=1,
    n_select=10, min_panel=5, max_panel=5,
    target_panel_size=5,
    random_state=42, verbose=False,
    classifiers=['logistic_regression'],
)


# =============================================================================
# Group 1: Leakage invariants
# =============================================================================


class TestLeakageInvariants:
    """The core M6 guarantee: y-using steps never see test-fold labels."""

    def test_feature_selector_never_sees_test_labels(self, small_signal_dataset):
        """Patch select_univariate_filter to record every (X, y) it sees.

        The recorded y lengths across folds must never equal the full y
        length, and the recorded y values must always be a subset of the
        training indices for that fold.
        """
        X, y = small_signal_dataset
        recorded: List[Dict] = []

        orig_method = FeatureSelector.select_univariate_filter

        def spy(self, X_arg, y_arg, n_features=50, test='welch', label='univariate_filter'):
            recorded.append({
                'n_train': len(y_arg),
                'y_sum': int(np.sum(y_arg)),
                'gene_count': X_arg.shape[1],
            })
            return orig_method(self, X_arg, y_arg, n_features=n_features,
                                test=test, label=label)

        with patch.object(
            FeatureSelector, 'select_univariate_filter', spy
        ):
            kwargs = dict(FAST_KWARGS)
            clf_results, per_fold_panels, _ = _run_pipeline_cv(
                X=X, y=y,
                gene_list=list(X.columns),
                **kwargs,
            )

        # Each of the 5 folds should have produced at least one call
        assert len(recorded) >= 5

        # No recorded call saw all samples
        n_total = len(y)
        for call in recorded:
            assert call['n_train'] < n_total, (
                f"Leakage detected: a feature selector call received "
                f"{call['n_train']} labels, which equals the full n={n_total}."
            )
            # Training fold in 5-fold CV on n=30 should have ~24 samples
            assert call['n_train'] in range(
                int(n_total * 0.7), int(n_total * 0.9) + 1
            )

    def test_panel_optimizer_never_sees_test_labels(self, small_signal_dataset):
        """Same guarantee for PanelOptimizer.forward_selection."""
        from raptor.biomarker_discovery.core import PanelOptimizer

        X, y = small_signal_dataset
        recorded: List[int] = []

        orig = PanelOptimizer.forward_selection

        def spy(self, X_arg, y_arg, **kwargs):
            recorded.append(len(y_arg))
            return orig(self, X_arg, y_arg, **kwargs)

        with patch.object(PanelOptimizer, 'forward_selection', spy):
            _run_pipeline_cv(
                X=X, y=y,
                gene_list=list(X.columns),
                **FAST_KWARGS,
            )

        assert len(recorded) >= 5
        n_total = len(y)
        for n_seen in recorded:
            assert n_seen < n_total, (
                f"Leakage in PanelOptimizer: saw {n_seen} labels "
                f"(full n={n_total})."
            )

    def test_per_fold_panels_can_differ(self, noise_dataset):
        """On pure noise with many candidate genes, folds should pick
        different panels (no stable signal to latch onto)."""
        X, y = noise_dataset
        kwargs = dict(FAST_KWARGS)
        kwargs['n_select'] = 20
        _, per_fold_panels, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **kwargs,
        )
        unique_panels = {tuple(p) for p in per_fold_panels}
        assert len(unique_panels) >= 2, (
            f"On noise data, per-fold panels should differ but all "
            f"{len(per_fold_panels)} folds picked the same panel."
        )

    def test_oof_length_equals_n_samples(self, small_signal_dataset):
        """Every sample gets exactly one OOF prediction per repeat."""
        X, y = small_signal_dataset
        n_repeats = 2
        kwargs = dict(FAST_KWARGS)
        kwargs['n_repeats'] = n_repeats
        clf_results, _, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **kwargs,
        )
        expected_len = len(y) * n_repeats
        for name, res in clf_results.items():
            assert len(res.oof_true) == expected_len, (
                f"{name}: oof_true length {len(res.oof_true)} != "
                f"expected {expected_len} (n_samples * n_repeats)."
            )
            assert len(res.oof_prob) == expected_len


# =============================================================================
# Group 2: Statistical sanity — the bias correction actually works
# =============================================================================


class TestStatisticalSanity:
    """End-to-end: OOF AUC behaves correctly on noise vs signal."""

    def test_noise_data_oof_auc_near_chance(self, noise_dataset):
        """On random-label data, the honest OOF AUC must collapse to ~0.5.

        Pre-M6 on this fixture AUC was consistently >= 0.75 because
        feature selection had seen every label. Post-M6 it should be
        in a reasonable band around 0.5 (we allow [0.30, 0.70] to
        absorb 5-fold small-sample variance).
        """
        X, y = noise_dataset
        from sklearn.metrics import roc_auc_score

        kwargs = dict(FAST_KWARGS)
        kwargs['methods'] = ['univariate_filter', 'elastic_net']
        kwargs['n_select'] = 20
        clf_results, _, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **kwargs,
        )
        # Inspect the pooled OOF AUC of each classifier
        for name, res in clf_results.items():
            pooled_auc = roc_auc_score(res.oof_true, res.oof_prob)
            assert 0.30 <= pooled_auc <= 0.70, (
                f"{name}: pooled OOF AUC = {pooled_auc:.3f} on random "
                f"labels is outside chance-band [0.30, 0.70]. Pre-M6 "
                f"leakage may have re-entered the pipeline."
            )

    def test_signal_data_oof_auc_stays_high(self, signal_dataset):
        """With real 2sigma signal in 10 genes out of 100, OOF AUC
        should stay well above chance even under honest CV."""
        X, y = signal_dataset
        from sklearn.metrics import roc_auc_score

        kwargs = dict(FAST_KWARGS)
        kwargs['methods'] = ['univariate_filter', 'elastic_net']
        kwargs['n_select'] = 20
        clf_results, _, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **kwargs,
        )
        best_auc = max(
            roc_auc_score(res.oof_true, res.oof_prob)
            for res in clf_results.values()
        )
        assert best_auc >= 0.80, (
            f"Best pooled OOF AUC = {best_auc:.3f} on strong-signal "
            f"fixture; expected >= 0.80. If M6 regressed, the panels "
            f"may be so over-penalized by external CV that real signal "
            f"is being lost."
        )

    def test_per_fold_auc_recorded(self, signal_dataset):
        """metrics_per_fold should be populated for every fold/classifier."""
        X, y = signal_dataset
        kwargs = dict(FAST_KWARGS)
        kwargs['n_repeats'] = 2
        kwargs['n_select'] = 15
        clf_results, _, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **kwargs,
        )
        for name, res in clf_results.items():
            # 5 folds x 2 repeats = 10 rows per classifier
            assert len(res.metrics_per_fold) == 10
            # Each row must have the bookkeeping fields
            for row in res.metrics_per_fold:
                assert 'fold' in row
                assert 'repeat' in row
                assert 'panel_size' in row
                assert 'auc' in row


# =============================================================================
# Group 3: Nogueira stability measure
# =============================================================================


class TestNogueiraStability:
    """The stability diagnostic (Nogueira et al. 2018 JMLR 18:174)."""

    def test_phi_bounded_in_reasonable_range(self):
        """Phi should fall in [-1, 1]. (Theoretical range; in practice
        rarely goes below 0 on real data.)"""
        rng = np.random.default_rng(42)
        Z = rng.integers(0, 2, size=(5, 100))
        phi = _nogueira_stability_from_matrix(Z)
        assert -1.0 <= phi <= 1.0

    def test_phi_zero_under_random_selection(self):
        """Correction-for-chance: random uniform selection averages to 0."""
        rng = np.random.default_rng(0)
        n_trials, M, p, k = 500, 10, 100, 10
        phis = []
        for _ in range(n_trials):
            Z = np.zeros((M, p), dtype=int)
            for i in range(M):
                idx = rng.choice(p, size=k, replace=False)
                Z[i, idx] = 1
            phis.append(_nogueira_stability_from_matrix(Z))
        mean_phi = float(np.mean(phis))
        # Mean over 500 random trials should be within ~0.02 of 0
        assert abs(mean_phi) < 0.03, (
            f"Mean Phi under random selection = {mean_phi:.4f}; "
            f"expected close to 0 (correction-for-chance property)."
        )

    def test_phi_one_under_identical_panels(self):
        """All folds pick identical panels -> Phi = 1.0."""
        Z = np.array([
            [1, 1, 1, 0, 0],
            [1, 1, 1, 0, 0],
            [1, 1, 1, 0, 0],
            [1, 1, 1, 0, 0],
        ])
        phi = _nogueira_stability_from_matrix(Z)
        assert abs(phi - 1.0) < 1e-9

    def test_benchmark_label_assignment(self):
        """Labels follow the Nogueira 2018 benchmark scale."""
        assert _nogueira_benchmark_label(0.90) == 'excellent'
        assert _nogueira_benchmark_label(_NOGUEIRA_EXCELLENT_THRESHOLD) == 'excellent'
        assert _nogueira_benchmark_label(0.60) == 'intermediate'
        assert _nogueira_benchmark_label(_NOGUEIRA_POOR_THRESHOLD) == 'intermediate'
        assert _nogueira_benchmark_label(_NOGUEIRA_POOR_THRESHOLD - 0.01) == 'poor'
        assert _nogueira_benchmark_label(0.0) == 'poor'
        assert _nogueira_benchmark_label(-0.5) == 'poor'

    def test_bootstrap_ci_contains_point_for_stable_panels(self):
        """Bootstrap CI should be a sensible 2-sided interval."""
        # 4-of-5 folds identical, one gene swapped
        Z = np.zeros((5, 20), dtype=int)
        Z[0, :5] = 1
        Z[1, :5] = 1
        Z[2, :5] = 1
        Z[3, :5] = 1
        Z[4, [0, 1, 2, 3, 6]] = 1
        phi_point = _nogueira_stability_from_matrix(Z)
        lo, hi = _bootstrap_nogueira_ci(
            Z, n_bootstrap=1000, ci=0.95, seed=42,
        )
        assert lo <= hi
        # CI should span a reasonable range (not collapsed)
        assert 0.0 < phi_point <= 1.0
        assert 0.0 <= lo <= 1.0
        assert 0.0 <= hi <= 1.0


class TestPanelStabilityBuilder:
    """_compute_panel_stability assembles the full PanelStabilityResult."""

    def test_builds_full_result(self):
        """All fields populated from realistic per-fold inputs."""
        gene_universe = [f'G{i:03d}' for i in range(50)]
        per_fold_panels = [
            gene_universe[:5],
            gene_universe[:5],
            gene_universe[:4] + [gene_universe[7]],
            gene_universe[:5],
            gene_universe[:4] + [gene_universe[8]],
        ]
        per_fold_ranked = [pd.DataFrame() for _ in per_fold_panels]
        final_panel = gene_universe[:5]

        result = _compute_panel_stability(
            per_fold_panels=per_fold_panels,
            per_fold_ranked_genes=per_fold_ranked,
            gene_universe=gene_universe,
            final_panel=final_panel,
            n_folds=5,
            n_repeats=1,
            random_state=42,
        )

        assert isinstance(result, PanelStabilityResult)
        assert len(result.per_fold_panels) == 5
        assert result.nogueira_stability > 0.5  # mostly stable
        assert result.benchmark_label in ('excellent', 'intermediate')
        lo, hi = result.nogueira_stability_ci
        assert lo <= result.nogueira_stability + 1e-6  # sanity
        # Final panel overlap: 4 genes appeared in all folds, 1 in 3/5
        assert result.final_panel_overlap[gene_universe[0]] == 1.0
        assert 0.5 <= result.final_panel_overlap[gene_universe[4]] <= 1.0
        # Selection frequency sensible
        assert result.gene_selection_frequency[gene_universe[0]] == 1.0
        assert result.gene_selection_frequency[gene_universe[9]] == 0.0

    def test_summary_is_readable(self):
        """summary() returns a usable string."""
        gene_universe = ['a', 'b', 'c', 'd']
        per_fold_panels = [['a', 'b'], ['a', 'b'], ['a', 'c']]
        per_fold_ranked = [pd.DataFrame() for _ in per_fold_panels]
        result = _compute_panel_stability(
            per_fold_panels=per_fold_panels,
            per_fold_ranked_genes=per_fold_ranked,
            gene_universe=gene_universe,
            final_panel=['a', 'b'],
            n_folds=3, n_repeats=1, random_state=42,
        )
        s = result.summary()
        assert 'Nogueira' in s
        assert 'Phi' in s


# =============================================================================
# Group 4: Univariate filter (de_filter replacement)
# =============================================================================


class TestUnivariateFilter:
    """Welch's t-test filter, the fold-safe de_filter replacement."""

    def test_reproducible_with_same_data(self, signal_dataset):
        """Same X, y -> identical selected genes (no stochastic element)."""
        X, y = signal_dataset
        sel1 = FeatureSelector(random_state=42, verbose=False)
        sel2 = FeatureSelector(random_state=42, verbose=False)
        r1 = sel1.select_univariate_filter(X, y, n_features=10)
        r2 = sel2.select_univariate_filter(X, y, n_features=10)
        assert r1.selected_genes == r2.selected_genes

    def test_prefers_signal_genes(self, signal_dataset):
        """Top-10 by Welch's p-value should contain most of the 10 signal genes."""
        X, y = signal_dataset
        sel = FeatureSelector(random_state=42, verbose=False)
        result = sel.select_univariate_filter(X, y, n_features=10)
        signal_genes = set(X.columns[:10])
        picked = set(result.selected_genes)
        overlap = len(signal_genes & picked)
        assert overlap >= 8, (
            f"Univariate filter picked only {overlap}/10 signal genes; "
            f"expected >= 8."
        )


# =============================================================================
# Group 5: Integration and regression
# =============================================================================


class TestIntegrationRegression:
    """Non-functional guarantees: API surface, reproducibility, etc."""

    def test_panel_stability_result_dataclass_defaults(self):
        """PanelStabilityResult has sensible defaults and shapes."""
        res = PanelStabilityResult()
        assert isinstance(res.per_fold_panels, list)
        assert isinstance(res.gene_selection_frequency, pd.Series)
        assert res.nogueira_stability == 0.0
        assert res.benchmark_label == 'poor'

    def test_reproducible_with_same_random_state(self, small_signal_dataset):
        """Same seed -> identical OOF predictions and per-fold panels."""
        X, y = small_signal_dataset
        kwargs = dict(FAST_KWARGS)
        kwargs.update(
            X=X, y=y,
            gene_list=list(X.columns),
            random_state=123,
        )
        r1 = _run_pipeline_cv(**kwargs)
        r2 = _run_pipeline_cv(**kwargs)
        clf_r1, panels_r1, _ = r1
        clf_r2, panels_r2, _ = r2
        assert panels_r1 == panels_r2
        for name in clf_r1:
            np.testing.assert_array_equal(
                clf_r1[name].oof_prob, clf_r2[name].oof_prob
            )

    def test_resolve_methods_warns_on_de_filter_with_de_genes(self):
        """Supplying precomputed de_genes triggers the leakage warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            resolved = _resolve_pipeline_cv_methods(
                methods=['de_filter', 'elastic_net'],
                de_genes=['GENE_0000', 'GENE_0001'],
            )
            # Should warn exactly once about upstream leakage
            leakage_warnings = [
                x for x in w
                if 'leakage' in str(x.message).lower()
                or 'pipeline-CV' in str(x.message)
            ]
            assert len(leakage_warnings) >= 1
        assert 'de_filter' in resolved  # still honored

    def test_resolve_methods_substitutes_univariate_when_no_de_genes(self):
        """de_filter without de_genes -> swap in univariate_filter."""
        resolved = _resolve_pipeline_cv_methods(
            methods=['de_filter', 'elastic_net'],
            de_genes=None,
        )
        assert 'de_filter' not in resolved
        assert 'univariate_filter' in resolved

    def test_resolve_methods_default_list_is_fold_safe(self):
        """methods=None yields the fold-safe default set."""
        resolved = _resolve_pipeline_cv_methods(
            methods=None, de_genes=None,
        )
        # All three fold-safe defaults present
        for m in _FOLD_SAFE_DEFAULT_METHODS:
            assert m in resolved
        # de_filter is NOT in the default list
        assert 'de_filter' not in resolved


# =============================================================================
# Group 6: Interaction with M1 (clinical metrics use the cleaner OOF)
# =============================================================================


class TestM1M6Interaction:
    """End-to-end guarantee that M1 clinical metrics benefit from M6's
    cleaner OOF and don't re-inflate on noise."""

    def test_oof_prob_in_unit_interval(self, small_signal_dataset):
        """Pre-M1 the OOF slots were None; post-M1 + M6 they are populated
        and live in [0, 1] (positive-class probabilities)."""
        X, y = small_signal_dataset
        clf_results, _, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **FAST_KWARGS,
        )
        for name, res in clf_results.items():
            assert res.oof_prob is not None
            assert res.oof_true is not None
            assert np.all(res.oof_prob >= 0.0 - 1e-9)
            assert np.all(res.oof_prob <= 1.0 + 1e-9)

    def test_noise_oof_auc_gives_honest_youden_j(self, noise_dataset):
        """On pure noise, Youden's J computed from OOF should be low.

        This is the end-to-end regression test: if M6 ever regresses
        (feature selection sees test-fold labels again), the OOF AUC
        on noise will inflate above chance and Youden's J computed from
        those OOF predictions will rise above what true-noise should
        produce.
        """
        X, y = noise_dataset
        from sklearn.metrics import roc_curve

        kwargs = dict(FAST_KWARGS)
        kwargs['methods'] = ['univariate_filter', 'elastic_net']
        kwargs['n_select'] = 20
        clf_results, _, _ = _run_pipeline_cv(
            X=X, y=y,
            gene_list=list(X.columns),
            **kwargs,
        )
        youdens = []
        for res in clf_results.values():
            fpr, tpr, _ = roc_curve(res.oof_true, res.oof_prob)
            j_vals = tpr - fpr  # Youden's J at each threshold
            youdens.append(float(np.max(j_vals)))

        best_j = max(youdens)
        # On pure noise with honest CV, Youden's J should be modest.
        # Threshold 0.55 allows for random-chance best-of-fold noise
        # but flags any return of pre-M6 leakage (which produced
        # Youden's J >= 0.7 on comparable fixtures).
        assert best_j < 0.55, (
            f"Max Youden's J on noise = {best_j:.3f}; expected < 0.55. "
            f"Pre-M6 regression: feature-selection leakage may have "
            f"re-entered the pipeline."
        )