"""
Unit tests for raptor.biomarker_discovery.signature_score

Covers:
    - SignatureScore construction and validation
    - build_signature_score() in diagnostic mode (logistic regression)
    - build_signature_score() in prognostic mode (Cox, if lifelines available)
    - .score() applies learned weights and z-score normalization correctly
    - .stratify() produces correct risk group assignments
    - .to_dict() serialization
    - Error handling: missing genes, invalid modes, missing survival data,
      bad label types, wrong risk group counts
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery import (
    SignatureScore,
    build_signature_score,
)
from raptor.biomarker_discovery.signature_score import (
    VALID_MODES,
    _LIFELINES_AVAILABLE,
)


# ---------------------------------------------------------------------------
# Synthetic data fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_diagnostic_data():
    """
    40 samples x 6 genes. Three 'panel' genes are informative
    (higher in cases), three 'noise' genes are random.
    """
    rng = np.random.RandomState(42)
    n_per_class = 20
    n_samples = n_per_class * 2

    panel = ["GENE_A", "GENE_B", "GENE_C"]
    noise = ["NOISE_1", "NOISE_2", "NOISE_3"]
    genes = panel + noise

    # Controls: mean 5 for all
    controls = rng.normal(5.0, 1.0, size=(n_per_class, len(genes)))
    # Cases: panel genes shifted +2, noise unchanged
    cases = rng.normal(5.0, 1.0, size=(n_per_class, len(genes)))
    cases[:, :3] += 2.0

    X = np.vstack([controls, cases])
    y = np.array([0] * n_per_class + [1] * n_per_class)

    sample_names = [f"sample_{i:02d}" for i in range(n_samples)]
    df = pd.DataFrame(X, index=sample_names, columns=genes)
    y_series = pd.Series(y, index=sample_names, name="label")

    return df, y_series, panel


@pytest.fixture
def synthetic_prognostic_data():
    """
    60 samples x 5 genes. Panel gene expression correlates with hazard.
    """
    rng = np.random.RandomState(7)
    n = 60
    genes = ["GENE_A", "GENE_B", "GENE_C", "NOISE_1", "NOISE_2"]
    X = rng.normal(0.0, 1.0, size=(n, len(genes)))

    # Linear predictor from panel genes only
    true_coef = np.array([0.8, 0.6, 0.4, 0.0, 0.0])
    lp = X @ true_coef

    # Exponential survival time inversely proportional to hazard
    baseline = 10.0
    time = baseline * np.exp(-lp) * rng.exponential(1.0, size=n)
    # 70% observed events, 30% censored
    event = rng.binomial(1, 0.7, size=n)

    sample_names = [f"sample_{i:02d}" for i in range(n)]
    df = pd.DataFrame(X, index=sample_names, columns=genes)
    time_s = pd.Series(time, index=sample_names, name="time")
    event_s = pd.Series(event, index=sample_names, name="event")
    panel = ["GENE_A", "GENE_B", "GENE_C"]

    return df, time_s, event_s, panel


# ---------------------------------------------------------------------------
# SignatureScore construction
# ---------------------------------------------------------------------------

class TestSignatureScoreConstruction:

    def test_minimal_construction(self):
        sig = SignatureScore(
            mode="diagnostic",
            weights={"G1": 0.5, "G2": -0.3},
            intercept=0.0,
            normalization={"G1": (5.0, 1.0), "G2": (4.0, 0.8)},
            cutoffs={"threshold": 0.2},
            risk_labels=["low", "high"],
        )
        assert sig.mode == "diagnostic"
        assert sig.panel_genes == ["G1", "G2"]
        assert sig.performance == {}

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="Invalid mode"):
            SignatureScore(
                mode="nonsense",
                weights={"G1": 0.5},
                intercept=0.0,
                normalization={"G1": (5.0, 1.0)},
                cutoffs={"threshold": 0.0},
                risk_labels=["low", "high"],
            )

    def test_mismatched_genes_raises(self):
        """weights and normalization must have the same keys."""
        with pytest.raises(ValueError, match="same gene keys"):
            SignatureScore(
                mode="diagnostic",
                weights={"G1": 0.5, "G2": 0.3},
                intercept=0.0,
                normalization={"G1": (5.0, 1.0)},  # missing G2
                cutoffs={"threshold": 0.0},
                risk_labels=["low", "high"],
            )

    def test_too_few_risk_labels_raises(self):
        with pytest.raises(ValueError, match="at least 2"):
            SignatureScore(
                mode="diagnostic",
                weights={"G1": 0.5},
                intercept=0.0,
                normalization={"G1": (5.0, 1.0)},
                cutoffs={"threshold": 0.0},
                risk_labels=["only_one"],
            )


# ---------------------------------------------------------------------------
# build_signature_score — diagnostic mode
# ---------------------------------------------------------------------------

class TestBuildDiagnostic:

    def test_two_group_diagnostic_runs(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        assert sig.mode == "diagnostic"
        assert sig.panel_genes == panel
        assert set(sig.weights.keys()) == set(panel)
        assert sig.risk_labels == ["low", "high"]
        assert "threshold" in sig.cutoffs
        assert "auc" in sig.performance
        # Given the synthetic signal, AUC should be high
        assert sig.performance["auc"] > 0.90

    def test_three_group_diagnostic_runs(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(
            X, y, panel_genes=panel, mode="diagnostic", risk_groups=3,
        )
        assert sig.risk_labels == ["low", "medium", "high"]
        assert {"low", "high"}.issubset(sig.cutoffs)
        assert sig.cutoffs["low"] < sig.cutoffs["high"]

    def test_weights_have_expected_sign(self, synthetic_diagnostic_data):
        """Panel genes were shifted UP in cases; their weights should be positive."""
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        for g in panel:
            assert sig.weights[g] > 0, f"Expected positive weight for {g}"

    def test_missing_panel_gene_raises(self, synthetic_diagnostic_data):
        X, y, _ = synthetic_diagnostic_data
        with pytest.raises(ValueError, match="missing required panel genes"):
            build_signature_score(
                X, y, panel_genes=["GENE_A", "NOT_IN_X"], mode="diagnostic",
            )

    def test_invalid_mode_raises(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        with pytest.raises(ValueError, match="Invalid mode"):
            build_signature_score(X, y, panel_genes=panel, mode="nonsense")

    def test_non_binary_labels_raise(self, synthetic_diagnostic_data):
        X, _, panel = synthetic_diagnostic_data
        y_bad = np.array([0, 1, 2] * 13 + [0])  # 3-class, length 40
        with pytest.raises(ValueError, match="binary"):
            build_signature_score(X, y_bad, panel_genes=panel, mode="diagnostic")

    def test_invalid_risk_groups_raises(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        with pytest.raises(ValueError, match="risk_groups must be 2 or 3"):
            build_signature_score(
                X, y, panel_genes=panel, mode="diagnostic", risk_groups=4,
            )


# ---------------------------------------------------------------------------
# build_signature_score — prognostic mode
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _LIFELINES_AVAILABLE, reason="lifelines not installed")
class TestBuildPrognostic:

    def test_two_group_prognostic_runs(self, synthetic_prognostic_data):
        X, time, event, panel = synthetic_prognostic_data
        sig = build_signature_score(
            X, y=None, panel_genes=panel, mode="prognostic",
            time=time, event=event,
        )
        assert sig.mode == "prognostic"
        assert "c_index" in sig.performance
        assert 0.5 < sig.performance["c_index"] <= 1.0

    def test_missing_time_raises(self, synthetic_prognostic_data):
        X, _, event, panel = synthetic_prognostic_data
        with pytest.raises(ValueError, match="time.*event"):
            build_signature_score(
                X, y=None, panel_genes=panel, mode="prognostic", event=event,
            )


# Guard: prognostic mode without lifelines should raise ImportError cleanly
class TestPrognosticWithoutLifelines:

    def test_raises_import_error_if_lifelines_missing(
        self, synthetic_prognostic_data, monkeypatch,
    ):
        from raptor.biomarker_discovery import signature_score as sig_mod

        monkeypatch.setattr(sig_mod, "_LIFELINES_AVAILABLE", False)

        X, time, event, panel = synthetic_prognostic_data
        with pytest.raises(ImportError, match="lifelines"):
            build_signature_score(
                X, y=None, panel_genes=panel, mode="prognostic",
                time=time, event=event,
            )


# ---------------------------------------------------------------------------
# .score() — applying to new data
# ---------------------------------------------------------------------------

class TestScoreApplication:

    def test_score_returns_series_with_correct_index(
        self, synthetic_diagnostic_data,
    ):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        scores = sig.score(X)
        assert isinstance(scores, pd.Series)
        assert list(scores.index) == list(X.index)
        assert len(scores) == len(X)

    def test_score_separates_classes(self, synthetic_diagnostic_data):
        """Cases (y=1) should have higher average scores than controls (y=0)."""
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        scores = sig.score(X)
        mean_case = scores[y == 1].mean()
        mean_ctrl = scores[y == 0].mean()
        assert mean_case > mean_ctrl

    def test_score_ignores_extra_genes(self, synthetic_diagnostic_data):
        """Extra genes in new data should not affect score output."""
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        X_extra = X.copy()
        X_extra["UNRELATED_GENE"] = 100.0
        scores_a = sig.score(X)
        scores_b = sig.score(X_extra)
        pd.testing.assert_series_equal(scores_a, scores_b)

    def test_score_missing_panel_gene_raises(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        X_incomplete = X.drop(columns=[panel[0]])
        with pytest.raises(ValueError, match="missing required panel genes"):
            sig.score(X_incomplete)


# ---------------------------------------------------------------------------
# .stratify() — risk group assignment
# ---------------------------------------------------------------------------

class TestStratify:

    def test_two_group_stratify(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        scores = sig.score(X)
        groups = sig.stratify(scores)
        assert isinstance(groups, pd.Series)
        assert set(groups.unique()).issubset({"low", "high"})
        # Controls should mostly be 'low', cases mostly 'high'
        ctrl_low_frac = (groups[y == 0] == "low").mean()
        case_high_frac = (groups[y == 1] == "high").mean()
        assert ctrl_low_frac > 0.7
        assert case_high_frac > 0.7

    def test_three_group_stratify(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(
            X, y, panel_genes=panel, mode="diagnostic", risk_groups=3,
        )
        scores = sig.score(X)
        groups = sig.stratify(scores)
        assert set(groups.unique()).issubset({"low", "medium", "high"})

    def test_stratify_with_bad_cutoffs_raises(self):
        """Three-group labels but two-group cutoffs -> should raise."""
        sig = SignatureScore(
            mode="diagnostic",
            weights={"G1": 0.5},
            intercept=0.0,
            normalization={"G1": (5.0, 1.0)},
            cutoffs={"threshold": 0.0},  # two-group cutoffs
            risk_labels=["low", "medium", "high"],  # three-group labels
        )
        scores = pd.Series([0.1, 0.5, 0.9])
        with pytest.raises(ValueError, match="Three-group"):
            sig.stratify(scores)

    def test_stratify_inverted_cutoffs_raises(self):
        """cutoffs['low'] >= cutoffs['high'] should raise."""
        sig = SignatureScore(
            mode="diagnostic",
            weights={"G1": 0.5},
            intercept=0.0,
            normalization={"G1": (5.0, 1.0)},
            cutoffs={"low": 0.8, "high": 0.2},  # inverted
            risk_labels=["low", "medium", "high"],
        )
        scores = pd.Series([0.1, 0.5, 0.9])
        with pytest.raises(ValueError, match="must be <"):
            sig.stratify(scores)


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

class TestSerialization:

    def test_to_dict_contains_all_fields(self, synthetic_diagnostic_data):
        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        d = sig.to_dict()
        expected = {
            "mode", "weights", "intercept", "normalization", "cutoffs",
            "risk_labels", "performance", "panel_genes",
        }
        assert set(d.keys()) == expected

    def test_to_dict_values_are_json_safe_types(
        self, synthetic_diagnostic_data,
    ):
        """to_dict output should not contain numpy scalars or arrays."""
        import json

        X, y, panel = synthetic_diagnostic_data
        sig = build_signature_score(X, y, panel_genes=panel, mode="diagnostic")
        d = sig.to_dict()
        # If this serializes without error, all values are JSON-safe
        s = json.dumps(d)
        assert isinstance(s, str)
        assert len(s) > 0


# ---------------------------------------------------------------------------
# Module-level sanity
# ---------------------------------------------------------------------------

def test_valid_modes_constant():
    assert VALID_MODES == ("diagnostic", "prognostic")