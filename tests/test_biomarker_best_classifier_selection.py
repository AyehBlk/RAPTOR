"""
Tests for M5: deterministic best-classifier selection with auditable tiebreak.

Covers ``raptor.biomarker_discovery.core._select_best_classifier``, which
replaces the previous ``max(clf_results.keys(), key=lambda k: ...auc)``
pattern that silently relied on Python's dict iteration order when AUCs
tied.

Rule precedence under test:
    1. Highest AUC.
    2. If tied on AUC (within auc_epsilon), highest F1.
    3. If tied on F1 (within f1_epsilon), interpretability preference
       order: LR > SVM > RF > XGBoost, then unknowns alphabetically.
"""

from __future__ import annotations

import pytest

from raptor.biomarker_discovery.core import (
    ClassificationResult,
    _select_best_classifier,
)


def _mk(name: str, auc: float, f1: float = 0.0) -> ClassificationResult:
    """Build a minimal ClassificationResult for tiebreak testing."""
    return ClassificationResult(model_name=name, auc=auc, f1=f1)


class TestSelectBestClassifier:
    """M5 tiebreak: AUC -> F1 -> interpretability preference."""

    def test_single_clear_winner_no_tiebreak(self):
        """Common case: one classifier has strictly higher AUC."""
        results = {
            'logistic_regression': _mk('logistic_regression', auc=0.80, f1=0.75),
            'random_forest':       _mk('random_forest',       auc=0.85, f1=0.80),
            'svm':                 _mk('svm',                 auc=0.70, f1=0.65),
        }
        assert _select_best_classifier(results) == 'random_forest'

    def test_auc_tie_resolved_by_f1(self):
        """Two classifiers tie on AUC; higher F1 wins."""
        results = {
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=0.90),
            'random_forest':       _mk('random_forest',       auc=1.000, f1=0.95),
            'svm':                 _mk('svm',                 auc=0.80,  f1=0.99),
        }
        # svm is excluded at the AUC step; RF beats LR on F1.
        assert _select_best_classifier(results) == 'random_forest'

    def test_full_tie_resolved_by_interpretability_preference(self):
        """Three-way tie on AUC and F1 -> LR wins per preference order.

        Dict insertion order deliberately does NOT start with LR, to
        confirm selection is by preference, not by iteration order.
        """
        results = {
            'random_forest':       _mk('random_forest',       auc=1.000, f1=1.000),
            'svm':                 _mk('svm',                 auc=1.000, f1=1.000),
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=1.000),
        }
        assert _select_best_classifier(results) == 'logistic_regression'

    def test_preference_order_svm_over_rf_when_lr_absent(self):
        """Preference order applies transitively when LR is excluded."""
        results = {
            'random_forest': _mk('random_forest', auc=1.000, f1=1.000),
            'svm':           _mk('svm',           auc=1.000, f1=1.000),
            'xgboost':       _mk('xgboost',       auc=1.000, f1=1.000),
        }
        assert _select_best_classifier(results) == 'svm'

    def test_preference_order_rf_over_xgb_when_lr_svm_absent(self):
        """RF beats XGBoost when only those two remain."""
        results = {
            'xgboost':       _mk('xgboost',       auc=1.000, f1=1.000),
            'random_forest': _mk('random_forest', auc=1.000, f1=1.000),
        }
        assert _select_best_classifier(results) == 'random_forest'

    def test_auc_epsilon_treats_near_equal_as_tied(self):
        """AUCs within 1e-6 are treated as tied, so F1 decides."""
        results = {
            'logistic_regression': _mk('logistic_regression', auc=1.000_000_0, f1=0.80),
            'random_forest':       _mk('random_forest',       auc=0.999_999_5, f1=0.95),
        }
        # AUC gap is 5e-7 < 1e-6 -> tied -> RF wins on F1.
        assert _select_best_classifier(results) == 'random_forest'

    def test_auc_epsilon_respects_non_tied_gap(self):
        """An AUC gap larger than epsilon is NOT treated as tied."""
        results = {
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=0.60),
            'random_forest':       _mk('random_forest',       auc=0.995, f1=0.99),
        }
        # AUC gap 5e-3 >> 1e-6 -> LR wins outright despite lower F1.
        assert _select_best_classifier(results) == 'logistic_regression'

    def test_f1_epsilon_treats_near_equal_as_tied(self):
        """F1s within 1e-6 fall through to interpretability preference."""
        results = {
            'random_forest':       _mk('random_forest',       auc=1.000, f1=1.000_000_0),
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=0.999_999_5),
        }
        # AUC tied; F1 gap 5e-7 < 1e-6 -> tied -> LR wins on preference.
        assert _select_best_classifier(results) == 'logistic_regression'

    def test_empty_clf_results_raises(self):
        """Empty input is a programming error, not a silent pick."""
        with pytest.raises(ValueError, match="empty"):
            _select_best_classifier({})

    def test_unknown_classifier_sorts_after_known_ones(self):
        """A name not in the preference tuple loses the preference tiebreak."""
        results = {
            'mystery_model':       _mk('mystery_model',       auc=1.000, f1=1.000),
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=1.000),
        }
        assert _select_best_classifier(results) == 'logistic_regression'

    def test_two_unknown_classifiers_break_alphabetically(self):
        """If only unknown names remain, alphabetical order is deterministic."""
        results = {
            'zeta_model':  _mk('zeta_model',  auc=1.000, f1=1.000),
            'alpha_model': _mk('alpha_model', auc=1.000, f1=1.000),
        }
        assert _select_best_classifier(results) == 'alpha_model'

    def test_tiebreak_fires_info_log_on_f1_resolution(self, caplog):
        """Log a line at INFO when F1 actually breaks a tie."""
        import logging
        caplog.set_level(logging.INFO, logger='raptor.biomarker_discovery.core')
        results = {
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=0.90),
            'random_forest':       _mk('random_forest',       auc=1.000, f1=0.95),
        }
        _select_best_classifier(results)
        assert any('resolved by F1' in rec.message for rec in caplog.records)

    def test_tiebreak_fires_info_log_on_preference_resolution(self, caplog):
        """Log a line at INFO when preference order breaks a tie."""
        import logging
        caplog.set_level(logging.INFO, logger='raptor.biomarker_discovery.core')
        results = {
            'random_forest':       _mk('random_forest',       auc=1.000, f1=1.000),
            'logistic_regression': _mk('logistic_regression', auc=1.000, f1=1.000),
        }
        _select_best_classifier(results)
        assert any(
            'resolved by interpretability preference' in rec.message
            for rec in caplog.records
        )

    def test_no_log_on_clear_winner(self, caplog):
        """Silent in the common case where nothing ties."""
        import logging
        caplog.set_level(logging.INFO, logger='raptor.biomarker_discovery.core')
        results = {
            'logistic_regression': _mk('logistic_regression', auc=0.80, f1=0.75),
            'random_forest':       _mk('random_forest',       auc=0.85, f1=0.80),
        }
        _select_best_classifier(results)
        assert not any(
            'tie' in rec.message.lower() for rec in caplog.records
        )