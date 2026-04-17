"""Tests for clinical_metrics module."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from raptor.biomarker_discovery.clinical_metrics import (
    ppv_npv_at_prevalence,
    bootstrap_ci,
    youdens_optimal_threshold,
    decision_curve_analysis,
    net_reclassification_improvement,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def binary_data() -> tuple:
    """100 samples with a decent classifier: AUC ~ 0.9."""
    rng = np.random.default_rng(42)
    n = 100
    y_true = np.array([0] * 50 + [1] * 50)
    # class-0 scores centered at 0.3, class-1 at 0.7
    y_score = np.concatenate([
        rng.normal(0.3, 0.15, 50),
        rng.normal(0.7, 0.15, 50),
    ])
    y_score = np.clip(y_score, 0, 1)
    return y_true, y_score


# ═══════════════════════════════════════════════════════════════════════════
# ppv_npv_at_prevalence
# ═══════════════════════════════════════════════════════════════════════════

class TestPPVNPV:

    def test_high_prevalence_high_ppv(self) -> None:
        result = ppv_npv_at_prevalence(0.90, 0.90, prevalence=0.50)
        assert result["ppv"] > 0.85

    def test_low_prevalence_drops_ppv(self) -> None:
        high = ppv_npv_at_prevalence(0.90, 0.90, prevalence=0.50)
        low = ppv_npv_at_prevalence(0.90, 0.90, prevalence=0.01)
        assert low["ppv"] < high["ppv"]
        assert low["ppv"] < 0.10  # classic Bayesian trap

    def test_low_prevalence_high_npv(self) -> None:
        result = ppv_npv_at_prevalence(0.90, 0.90, prevalence=0.01)
        assert result["npv"] > 0.99

    def test_perfect_test(self) -> None:
        result = ppv_npv_at_prevalence(1.0, 1.0, prevalence=0.05)
        assert abs(result["ppv"] - 1.0) < 1e-10
        assert abs(result["npv"] - 1.0) < 1e-10

    def test_returns_input_values(self) -> None:
        result = ppv_npv_at_prevalence(0.80, 0.75, prevalence=0.10)
        assert result["sensitivity"] == 0.80
        assert result["specificity"] == 0.75
        assert result["prevalence"] == 0.10

    def test_invalid_sensitivity_raises(self) -> None:
        with pytest.raises(ValueError, match="sensitivity"):
            ppv_npv_at_prevalence(0.0, 0.90, 0.10)

    def test_invalid_specificity_raises(self) -> None:
        with pytest.raises(ValueError, match="specificity"):
            ppv_npv_at_prevalence(0.90, 1.5, 0.10)

    def test_invalid_prevalence_raises(self) -> None:
        with pytest.raises(ValueError, match="prevalence"):
            ppv_npv_at_prevalence(0.90, 0.90, 0.0)

    def test_prevalence_at_boundary(self) -> None:
        with pytest.raises(ValueError, match="prevalence"):
            ppv_npv_at_prevalence(0.90, 0.90, 1.0)

    def test_known_calculation(self) -> None:
        # Manual: sens=0.9, spec=0.9, prev=0.5
        # PPV = 0.9*0.5 / (0.9*0.5 + 0.1*0.5) = 0.45/0.50 = 0.9
        result = ppv_npv_at_prevalence(0.9, 0.9, 0.5)
        assert abs(result["ppv"] - 0.9) < 1e-10


# ═══════════════════════════════════════════════════════════════════════════
# bootstrap_ci
# ═══════════════════════════════════════════════════════════════════════════

class TestBootstrapCI:

    def test_returns_expected_keys(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = bootstrap_ci(y_true, y_score)
        assert set(result.keys()) == {
            "point_estimate", "ci_lower", "ci_upper", "ci_level", "n_bootstrap"
        }

    def test_ci_contains_point_estimate(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = bootstrap_ci(y_true, y_score)
        assert result["ci_lower"] <= result["point_estimate"] <= result["ci_upper"]

    def test_ci_width_is_reasonable(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = bootstrap_ci(y_true, y_score)
        width = result["ci_upper"] - result["ci_lower"]
        assert 0.01 < width < 0.30  # not too narrow, not too wide for n=100

    def test_wider_ci_with_higher_confidence(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        ci90 = bootstrap_ci(y_true, y_score, ci=0.90)
        ci99 = bootstrap_ci(y_true, y_score, ci=0.99)
        width90 = ci90["ci_upper"] - ci90["ci_lower"]
        width99 = ci99["ci_upper"] - ci99["ci_lower"]
        assert width99 > width90

    def test_point_estimate_matches_metric(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = bootstrap_ci(y_true, y_score)
        direct_auc = roc_auc_score(y_true, y_score)
        assert abs(result["point_estimate"] - direct_auc) < 1e-10

    def test_custom_metric_function(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data

        def accuracy_at_50(yt, ys):
            pred = (ys >= 0.5).astype(int)
            return float(np.mean(pred == yt))

        result = bootstrap_ci(y_true, y_score, metric_fn=accuracy_at_50)
        assert 0.5 < result["point_estimate"] < 1.0

    def test_reproducible_with_same_seed(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        r1 = bootstrap_ci(y_true, y_score, seed=123)
        r2 = bootstrap_ci(y_true, y_score, seed=123)
        assert r1["ci_lower"] == r2["ci_lower"]
        assert r1["ci_upper"] == r2["ci_upper"]

    def test_different_seeds_give_different_results(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        r1 = bootstrap_ci(y_true, y_score, seed=1)
        r2 = bootstrap_ci(y_true, y_score, seed=999)
        assert r1["ci_lower"] != r2["ci_lower"]

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="length"):
            bootstrap_ci(np.array([0, 1]), np.array([0.5]))

    def test_invalid_ci_raises(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        with pytest.raises(ValueError, match="ci must be"):
            bootstrap_ci(y_true, y_score, ci=0.30)


# ═══════════════════════════════════════════════════════════════════════════
# youdens_optimal_threshold
# ═══════════════════════════════════════════════════════════════════════════

class TestYoudensThreshold:

    def test_returns_expected_keys(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = youdens_optimal_threshold(y_true, y_score)
        assert set(result.keys()) == {
            "threshold", "sensitivity", "specificity", "youdens_j", "auc"
        }

    def test_threshold_in_score_range(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = youdens_optimal_threshold(y_true, y_score)
        assert y_score.min() <= result["threshold"] <= y_score.max() + 1

    def test_j_is_positive_for_good_classifier(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = youdens_optimal_threshold(y_true, y_score)
        assert result["youdens_j"] > 0.5

    def test_sensitivity_specificity_in_range(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = youdens_optimal_threshold(y_true, y_score)
        assert 0 <= result["sensitivity"] <= 1
        assert 0 <= result["specificity"] <= 1

    def test_perfect_separation(self) -> None:
        y_true = np.array([0, 0, 0, 1, 1, 1])
        y_score = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
        result = youdens_optimal_threshold(y_true, y_score)
        assert result["youdens_j"] == 1.0
        assert result["sensitivity"] == 1.0
        assert result["specificity"] == 1.0

    def test_single_class_raises(self) -> None:
        with pytest.raises(ValueError, match="must contain both"):
            youdens_optimal_threshold(np.array([1, 1, 1]), np.array([0.5, 0.6, 0.7]))

    def test_auc_matches_sklearn(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        result = youdens_optimal_threshold(y_true, y_score)
        assert abs(result["auc"] - roc_auc_score(y_true, y_score)) < 1e-10


# ═══════════════════════════════════════════════════════════════════════════
# decision_curve_analysis
# ═══════════════════════════════════════════════════════════════════════════

class TestDecisionCurve:

    def test_returns_dataframe(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        df = decision_curve_analysis(y_true, y_score)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        df = decision_curve_analysis(y_true, y_score)
        assert set(df.columns) == {
            "threshold", "net_benefit_model",
            "net_benefit_treat_all", "net_benefit_treat_none",
        }

    def test_treat_none_always_zero(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        df = decision_curve_analysis(y_true, y_score)
        assert (df["net_benefit_treat_none"] == 0.0).all()

    def test_model_beats_treat_all_sometimes(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        df = decision_curve_analysis(y_true, y_score)
        advantage = df["net_benefit_model"] > df["net_benefit_treat_all"]
        assert advantage.any(), "Good model should beat treat-all at some thresholds"

    def test_custom_thresholds(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        th = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        df = decision_curve_analysis(y_true, y_score, thresholds=th)
        assert len(df) == 5

    def test_default_threshold_count(self, binary_data: tuple) -> None:
        y_true, y_score = binary_data
        df = decision_curve_analysis(y_true, y_score)
        assert len(df) == 99  # 0.01 to 0.99

    def test_single_class_raises(self) -> None:
        with pytest.raises(ValueError, match="must contain both"):
            decision_curve_analysis(np.array([0, 0, 0]), np.array([0.1, 0.2, 0.3]))


# ═══════════════════════════════════════════════════════════════════════════
# net_reclassification_improvement
# ═══════════════════════════════════════════════════════════════════════════

class TestNRI:

    def test_returns_expected_keys(self) -> None:
        y = np.array([1, 1, 0, 0, 1, 0])
        old = np.array([0.3, 0.4, 0.6, 0.3, 0.5, 0.4])
        new = np.array([0.6, 0.7, 0.1, 0.1, 0.8, 0.2])
        result = net_reclassification_improvement(y, old, new)
        expected_keys = {
            "nri", "nri_events", "nri_nonevents", "z_score", "p_value",
            "n_events", "n_nonevents", "events_up", "events_down",
            "nonevents_up", "nonevents_down",
        }
        assert set(result.keys()) == expected_keys

    def test_perfect_reclassification(self) -> None:
        # Old model: everyone at medium risk. New model: events high, non-events low.
        y = np.array([1, 1, 1, 0, 0, 0])
        old = np.array([0.35, 0.35, 0.35, 0.35, 0.35, 0.35])
        new = np.array([0.90, 0.90, 0.90, 0.05, 0.05, 0.05])
        result = net_reclassification_improvement(y, old, new)
        assert result["nri_events"] == 1.0      # all events moved up
        assert result["nri_nonevents"] == 1.0    # all non-events moved down
        assert result["nri"] == 2.0

    def test_no_change_nri_zero(self) -> None:
        y = np.array([1, 0, 1, 0])
        risk = np.array([0.3, 0.3, 0.7, 0.7])
        result = net_reclassification_improvement(y, risk, risk)
        assert result["nri"] == 0.0
        assert result["events_up"] == 0
        assert result["events_down"] == 0

    def test_harmful_biomarker_negative_nri(self) -> None:
        # New model moves events DOWN and non-events UP -> negative NRI
        y = np.array([1, 1, 1, 0, 0, 0])
        old = np.array([0.60, 0.60, 0.60, 0.10, 0.10, 0.10])
        new = np.array([0.10, 0.10, 0.10, 0.60, 0.60, 0.60])
        result = net_reclassification_improvement(y, old, new)
        assert result["nri"] < 0

    def test_custom_cutoffs(self) -> None:
        y = np.array([1, 0, 1, 0])
        old = np.array([0.15, 0.15, 0.45, 0.45])
        new = np.array([0.55, 0.05, 0.85, 0.05])
        result = net_reclassification_improvement(y, old, new,
                                                   cutoffs=[0.10, 0.40, 0.70])
        assert result["nri"] > 0

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            net_reclassification_improvement(
                np.array([0, 1]),
                np.array([0.5, 0.5]),
                np.array([0.5]),
            )

    def test_single_class_raises(self) -> None:
        with pytest.raises(ValueError, match="must contain both"):
            net_reclassification_improvement(
                np.array([1, 1, 1]),
                np.array([0.5, 0.5, 0.5]),
                np.array([0.6, 0.6, 0.6]),
            )

    def test_p_value_significant_for_good_reclassification(self) -> None:
        rng = np.random.default_rng(42)
        n = 200
        y = np.array([1] * 100 + [0] * 100)
        old = np.clip(rng.normal(0.5, 0.15, n), 0, 1)
        # New model: events get boosted, non-events get reduced
        new = old.copy()
        new[:100] += 0.25
        new[100:] -= 0.25
        new = np.clip(new, 0, 1)
        result = net_reclassification_improvement(y, old, new)
        assert result["nri"] > 0
        assert result["p_value"] < 0.05