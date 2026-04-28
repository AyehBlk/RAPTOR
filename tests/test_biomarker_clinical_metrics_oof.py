"""
Tests for M1: honest clinical metrics via out-of-fold predictions.

These tests verify three intertwined guarantees:

1. ``evaluate_nested_cv`` and ``evaluate_loocv`` in ``core.py`` persist
   their out-of-fold (OOF) predictions on the returned
   ``ClassificationResult`` (via new fields ``oof_true`` and
   ``oof_prob``).
2. ``enhance_biomarker_result`` in ``enhanced.py`` prefers those OOF
   predictions when computing clinical metrics (Youden's threshold,
   bootstrap CI, PPV/NPV, DCA), setting ``oof_used=True`` on the
   returned ``clinical_metrics`` dict.
3. When OOF predictions are absent, the clinical-metrics block falls
   back to full-data predictions but attaches a ``warning`` note and
   sets ``oof_used=False`` so the user is informed rather than misled.

The critical statistical guarantee is tested via a shuffled-label
(null-signal) dataset: full-data predictions produce inflated Youden's
J (approaching 1.0) while OOF predictions produce near-chance Youden's J
(below 0.4 on the fixture). This is the difference M1 makes visible.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery.core import (
    ClassificationResult,
    ClassifierEvaluator,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def real_signal_dataset():
    """60 samples, 10 features; first 5 carry separating signal.

    Samples are labeled ``S000``..``S059`` (strings) so they align
    with metadata.sample_id in the integration-style tests below.
    """
    rng = np.random.default_rng(42)
    n_samples, n_features = 60, 10
    y = np.array([0] * 30 + [1] * 30)
    X = rng.standard_normal((n_samples, n_features))
    X[y == 1, :5] += 2.5  # strong signal in first 5 features
    sample_ids = [f'S{i:03d}' for i in range(n_samples)]
    return (
        pd.DataFrame(
            X,
            columns=[f'G{i}' for i in range(n_features)],
            index=sample_ids,
        ),
        y,
    )


@pytest.fixture
def null_signal_dataset():
    """60 samples, 20 features; labels shuffled relative to data.

    Used to test that OOF predictions collapse toward chance while
    full-data predictions overfit to the shuffled labels.
    """
    rng = np.random.default_rng(42)
    n_samples, n_features = 60, 20
    y = np.array([0] * 30 + [1] * 30)
    rng.shuffle(y)
    X = rng.standard_normal((n_samples, n_features))
    sample_ids = [f'S{i:03d}' for i in range(n_samples)]
    return (
        pd.DataFrame(
            X,
            columns=[f'G{i:02d}' for i in range(n_features)],
            index=sample_ids,
        ),
        y,
    )


@pytest.fixture
def small_dataset():
    """10 samples for LOOCV path testing."""
    rng = np.random.default_rng(7)
    n_samples, n_features = 10, 5
    y = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    X = rng.standard_normal((n_samples, n_features))
    X[y == 1, :3] += 2.0
    sample_ids = [f'S{i:03d}' for i in range(n_samples)]
    return (
        pd.DataFrame(
            X,
            columns=[f'G{i}' for i in range(n_features)],
            index=sample_ids,
        ),
        y,
    )


@pytest.fixture
def real_signal_counts_fixture():
    """Raw-count fixture for integration tests through discover_biomarkers.

    Produces a counts matrix (genes x samples) with a realistic mean in
    the hundreds so that RAPTOR's `_prepare_expression_data` low-count
    filter (gene_means >= min_count=10) keeps all genes. The first
    handful of genes carry a strong separating signal between the two
    groups; the rest are noise.

    Returns
    -------
    counts : pd.DataFrame (genes x samples)
        Raw counts, non-negative integers.
    metadata : pd.DataFrame
        Columns: sample_id, condition.
    sample_ids, y : for reference
    """
    rng = np.random.default_rng(42)
    n_samples, n_genes, n_signal = 60, 40, 8
    y = np.array([0] * 30 + [1] * 30)

    # Baseline counts around 500 per gene per sample; Negative Binomial-ish
    # variability to mimic real RNA-seq. Signal = multiplicative boost
    # in one group for the signal genes.
    base = rng.poisson(lam=500, size=(n_genes, n_samples)).astype(float)
    # Boost first n_signal genes in disease (y==1) samples
    for i in range(n_signal):
        base[i, y == 1] *= 4.0   # 4x upregulation on signal genes
    counts = base.astype(int)

    sample_ids = [f'S{i:03d}' for i in range(n_samples)]
    gene_ids = [f'G{i:03d}' for i in range(n_genes)]
    counts_df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)
    counts_df.index.name = 'gene_id'

    metadata = pd.DataFrame({
        'sample_id': sample_ids,
        'condition': ['healthy' if lbl == 0 else 'disease' for lbl in y],
    })

    return counts_df, metadata, sample_ids, y


# ---------------------------------------------------------------------------
# core.py: OOF persistence
# ---------------------------------------------------------------------------

class TestNestedCVPersistsOOF:
    """evaluate_nested_cv stores per-sample OOF predictions."""

    def test_oof_fields_populated(self, real_signal_dataset):
        X, y = real_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']
        assert r.oof_true is not None, "oof_true not persisted"
        assert r.oof_prob is not None, "oof_prob not persisted"

    def test_oof_length_equals_n_samples(self, real_signal_dataset):
        """Every sample is held out exactly once across the outer folds."""
        X, y = real_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']
        assert len(r.oof_true) == len(y)
        assert len(r.oof_prob) == len(y)

    def test_oof_true_matches_y_values(self, real_signal_dataset):
        """OOF true labels are a permutation of the original y."""
        X, y = real_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']
        # Sorted y and sorted oof_true should be identical: same class balance
        assert sorted(r.oof_true.tolist()) == sorted(y.tolist())

    def test_oof_prob_in_unit_interval(self, real_signal_dataset):
        X, y = real_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']
        assert (r.oof_prob >= 0.0).all()
        assert (r.oof_prob <= 1.0).all()

    def test_oof_persisted_across_all_classifiers(self, real_signal_dataset):
        X, y = real_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5,
            classifiers=['logistic_regression', 'random_forest', 'svm'],
        )
        for name in ['logistic_regression', 'random_forest', 'svm']:
            r = results[name]
            assert r.oof_true is not None, f"{name}: oof_true missing"
            assert r.oof_prob is not None, f"{name}: oof_prob missing"
            assert len(r.oof_prob) == len(y), f"{name}: wrong OOF length"


class TestLOOCVPersistsOOF:
    """evaluate_loocv stores per-sample OOF predictions (one per sample)."""

    def test_oof_fields_populated(self, small_dataset):
        X, y = small_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_loocv(
            X, y, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']
        assert r.oof_true is not None
        assert r.oof_prob is not None
        assert len(r.oof_true) == len(y)
        assert len(r.oof_prob) == len(y)


# ---------------------------------------------------------------------------
# enhanced.py: clinical metrics prefer OOF
# ---------------------------------------------------------------------------

class TestClinicalMetricsUseOOF:
    """The clinical-metrics block reads OOF fields when available."""

    def test_oof_used_flag_true_when_nested_cv_populates_oof(
        self, real_signal_counts_fixture,
    ):
        from raptor.biomarker_discovery import discover_biomarkers

        counts, metadata, _, _ = real_signal_counts_fixture

        result = discover_biomarkers(
            counts=counts,
            metadata=metadata,
            group_column='condition',
            reference_group='disease',
            baseline_group='healthy',
            intent='diagnostic',
            min_panel=3,
            max_panel=5,
            validation='nested_cv',
            n_folds=5,
            annotate=False,
            run_literature=False,
            run_ppi=False,
            random_state=42,
            verbose=False,
        )

        assert result.clinical_metrics is not None
        assert result.clinical_metrics.get('oof_used') is True, (
            "Clinical metrics should use OOF predictions when "
            "evaluate_nested_cv has populated them."
        )

    def test_oof_used_false_on_bare_fallback(self, real_signal_counts_fixture):
        """Manually clear OOF fields; enhanced.py should fall back and warn."""
        from raptor.biomarker_discovery import discover_biomarkers
        from raptor.biomarker_discovery.enhanced import enhance_biomarker_result

        counts, metadata, sample_ids, y = real_signal_counts_fixture

        result = discover_biomarkers(
            counts=counts,
            metadata=metadata,
            group_column='condition',
            reference_group='disease',
            baseline_group='healthy',
            intent='diagnostic',
            min_panel=3,
            max_panel=5,
            validation='nested_cv',
            n_folds=5,
            annotate=False,
            run_literature=False,
            run_ppi=False,
            random_state=42,
            verbose=False,
        )
        # Simulate an older saved result that lacks OOF fields by
        # clearing them, then re-run the enhanced analysis.
        best_name = result.best_classifier
        result.base_result.classification_results[best_name].oof_true = None
        result.base_result.classification_results[best_name].oof_prob = None

        # Rebuild expression matrix the way discover_biomarkers would
        # have done internally: log2-CPM of filtered counts, samples x
        # genes, aligned to the metadata order.
        lib_sizes = counts.sum(axis=0)
        cpm = counts.div(lib_sizes, axis=1) * 1e6
        log_cpm = np.log2(cpm + 1)
        # Filter low-count genes the same way _prepare_expression_data does
        kept = counts.mean(axis=1) >= 10
        expr = log_cpm.loc[kept].T  # samples x genes

        labels_arr = np.array(
            [1 if c == 'disease' else 0 for c in metadata['condition']],
            dtype=int,
        )
        enhanced = enhance_biomarker_result(
            base_result=result.base_result,
            expression=expr,
            labels=labels_arr,
            group_names=('healthy', 'disease'),
            intent='diagnostic',
            verbose=False,
        )
        clin = enhanced.clinical_metrics
        assert clin is not None
        assert clin.get('oof_used') is False
        assert 'warning' in clin
        assert 'training-data' in clin['warning'].lower()


# ---------------------------------------------------------------------------
# The statistical payoff: honest metrics on null signal
# ---------------------------------------------------------------------------

class TestHonestMetricsOnNullSignal:
    """The core paper claim: OOF drops clinical metrics toward chance on noise."""

    def test_oof_auc_far_below_training_auc_on_null_signal(
        self, null_signal_dataset,
    ):
        from sklearn.metrics import roc_auc_score

        X, y = null_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']

        # OOF AUC on noise: expected near 0.5
        oof_auc = roc_auc_score(r.oof_true, r.oof_prob)

        # Training AUC: apply the full-data model to full data
        train_prob = r.trained_model.predict_proba(X)[:, 1]
        train_auc = roc_auc_score(y, train_prob)

        # Training AUC should be substantially inflated relative to OOF
        assert train_auc > oof_auc, (
            f"Training AUC ({train_auc:.3f}) should exceed OOF AUC "
            f"({oof_auc:.3f}) on null-signal data; otherwise the test "
            f"fixture is not noisy enough."
        )
        # And OOF should be below 0.70 (generously: it's noise, should
        # be near 0.5 but small sample can jitter)
        assert oof_auc < 0.70, (
            f"OOF AUC on null-signal should be near chance, got {oof_auc:.3f}"
        )
        # Training AUC on this fixture should be clearly high
        assert train_auc > 0.80, (
            f"Training AUC should be clearly inflated on null signal "
            f"with n=60 p=20, got {train_auc:.3f}"
        )

    def test_youden_j_drops_with_oof_on_null_signal(self, null_signal_dataset):
        """The fix's clinical punchline: Youden's J on noise drops dramatically."""
        from sklearn.metrics import roc_curve

        def youden_j(y_true, y_prob):
            fpr, tpr, _ = roc_curve(y_true, y_prob)
            return float((tpr - fpr).max())

        X, y = null_signal_dataset
        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X, y, n_outer=5, classifiers=['logistic_regression'],
        )
        r = results['logistic_regression']

        oof_j = youden_j(r.oof_true, r.oof_prob)
        train_prob = r.trained_model.predict_proba(X)[:, 1]
        train_j = youden_j(y, train_prob)

        # On noise, training Youden typically saturates near 1.0;
        # OOF Youden should be substantially lower. A gap of 0.3+ is
        # expected on this fixture.
        gap = train_j - oof_j
        assert gap > 0.30, (
            f"Expected a Youden gap > 0.30 on null-signal data; got "
            f"train_j={train_j:.3f}, oof_j={oof_j:.3f}, gap={gap:.3f}. "
            f"If gap is small, M1 is not routing through OOF correctly."
        )


# ---------------------------------------------------------------------------
# Optimism gap diagnostic
# ---------------------------------------------------------------------------

class TestOptimismGapDiagnostic:
    """enhanced.py reports cv_auc_oof, training_auc, gap when OOF is used."""

    def test_optimism_gap_fields_present_with_oof(self, real_signal_counts_fixture):
        from raptor.biomarker_discovery import discover_biomarkers

        counts, metadata, _, _ = real_signal_counts_fixture
        result = discover_biomarkers(
            counts=counts,
            metadata=metadata,
            group_column='condition',
            reference_group='disease',
            baseline_group='healthy',
            intent='diagnostic',
            min_panel=3,
            max_panel=5,
            validation='nested_cv',
            n_folds=5,
            annotate=False,
            run_literature=False,
            run_ppi=False,
            random_state=42,
            verbose=False,
        )
        assert result.clinical_metrics is not None
        og = result.clinical_metrics.get('optimism_gap')
        assert og is not None, "optimism_gap block missing"
        assert 'cv_auc_oof' in og
        assert 'training_auc' in og
        assert 'gap' in og
        # training_auc >= cv_auc_oof on average (can equal at AUC=1.0 saturation)
        assert og['training_auc'] >= og['cv_auc_oof'] - 1e-9


# ---------------------------------------------------------------------------
# Backwards-compat: a ClassificationResult with no OOF still works
# ---------------------------------------------------------------------------

class TestBackwardsCompatibility:
    """Existing code that doesn't set oof_* fields should still work."""

    def test_classification_result_defaults_to_none_oof(self):
        r = ClassificationResult(model_name='test', auc=0.8, f1=0.7)
        assert r.oof_true is None
        assert r.oof_prob is None

    def test_classification_result_accepts_oof_kwargs(self):
        r = ClassificationResult(
            model_name='test',
            auc=0.8,
            oof_true=np.array([0, 1, 0, 1]),
            oof_prob=np.array([0.2, 0.8, 0.3, 0.7]),
        )
        assert r.oof_true is not None
        assert len(r.oof_prob) == 4