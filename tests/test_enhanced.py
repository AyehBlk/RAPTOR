"""Tests for the enhanced biomarker analysis integration."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
from unittest.mock import MagicMock

from raptor.biomarker_discovery.enhanced import (
    EnhancedBiomarkerResult,
    enhance_biomarker_result,
)
from raptor.biomarker_discovery.intent import BiomarkerIntent
from raptor.biomarker_discovery.signature_score import SignatureScore
from raptor.biomarker_discovery.direction_patterns import DirectionPattern
from raptor.biomarker_discovery.ratio_biomarkers import RatioSearchResult


# ---------------------------------------------------------------------------
# Mock BiomarkerResult (avoids running the full discovery pipeline)
# ---------------------------------------------------------------------------

@dataclass
class _MockClassificationResult:
    model_name: str = "logistic_regression"
    accuracy: float = 0.90
    auc: float = 0.92
    sensitivity: float = 0.88
    specificity: float = 0.92
    f1: float = 0.89
    metrics_per_fold: List[Dict] = field(default_factory=list)
    confusion_matrix: Optional[Any] = None
    roc_data: Optional[Dict] = None
    feature_importance: Optional[pd.DataFrame] = None
    trained_model: Optional[Any] = None


@dataclass
class _MockPanelOptimization:
    optimal_panel: List[str] = field(default_factory=list)
    optimal_size: int = 5
    optimal_auc: float = 0.90
    panel_curve: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame({"panel_size": [3, 5], "auc_mean": [0.85, 0.90]})
    )


@dataclass
class _MockBiomarkerResult:
    ranked_genes: pd.DataFrame = field(default_factory=lambda: pd.DataFrame())
    panel: List[str] = field(default_factory=list)
    panel_size: int = 0
    selection_results: Dict = field(default_factory=dict)
    classification_results: Dict = field(default_factory=dict)
    best_classifier: str = ""
    panel_optimization: Optional[Any] = None
    annotation_result: Optional[Any] = None
    study_design: str = "binary"
    validation_strategy: str = "nested_cv"
    n_samples: int = 0
    n_initial_candidates: int = 0
    parameters: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)
    timestamp: str = "2026-04-17T00:00:00"

    def summary(self) -> str:
        return f"MockResult: {self.panel_size} genes, AUC={self.panel_optimization.optimal_auc if self.panel_optimization else 'N/A'}"

    def save(self, output_dir):
        pass


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_setup() -> dict:
    """
    Build a synthetic expression dataset + mock BiomarkerResult with a
    trained logistic regression model. Returns all components needed
    for enhance_biomarker_result().
    """
    from sklearn.linear_model import LogisticRegression

    rng = np.random.default_rng(42)
    n_per = 30
    n_genes = 20

    # Build expression: first 5 genes are UP in disease, rest are noise
    base = rng.normal(5.0, 1.0, size=(n_per, n_genes))
    dis = rng.normal(5.0, 1.0, size=(n_per, n_genes))
    dis[:, 0:5] += 3.0

    X_arr = np.vstack([base, dis])
    gene_names = [f"g{i:02d}" for i in range(n_genes)]
    expression = pd.DataFrame(X_arr, columns=gene_names)
    labels = np.array([0] * n_per + [1] * n_per)

    # The "discovered" panel is the 5 signal genes
    panel_genes = [f"g{i:02d}" for i in range(5)]

    # Train a real logistic regression on the panel genes
    X_panel = expression[panel_genes].to_numpy()
    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X_panel, labels)

    clf_result = _MockClassificationResult(
        model_name="logistic_regression",
        auc=0.95,
        sensitivity=0.90,
        specificity=0.93,
        trained_model=lr,
    )

    mock_result = _MockBiomarkerResult(
        panel=panel_genes,
        panel_size=len(panel_genes),
        classification_results={"logistic_regression": clf_result},
        best_classifier="logistic_regression",
        panel_optimization=_MockPanelOptimization(
            optimal_panel=panel_genes, optimal_size=5
        ),
        n_samples=n_per * 2,
        n_initial_candidates=n_genes,
    )

    return {
        "base_result": mock_result,
        "expression": expression,
        "labels": labels,
        "group_names": ("healthy", "disease"),  # 0=healthy, 1=disease
        "panel_genes": panel_genes,
    }


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestEnhanceBasic:

    def test_returns_enhanced_result(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert isinstance(result, EnhancedBiomarkerResult)

    def test_intent_stored(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.intent.intent == "diagnostic"

    def test_accepts_intent_object(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        intent_obj = BiomarkerIntent("diagnostic")
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent=intent_obj,
        )
        assert result.intent is intent_obj

    def test_base_result_accessible(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.panel == s["panel_genes"]
        assert result.panel_size == 5
        assert result.best_classifier == "logistic_regression"


# ---------------------------------------------------------------------------
# Signature Score
# ---------------------------------------------------------------------------

class TestEnhancedSignature:

    def test_diagnostic_builds_signature(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.signature is not None
        assert isinstance(result.signature, SignatureScore)
        assert result.signature.mode == "diagnostic"

    def test_exploratory_skips_signature(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="exploratory",
        )
        assert result.signature is None

    def test_translational_skips_signature(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        # Need species args for translational
        intent = BiomarkerIntent("translational",
                                 source_species="mouse",
                                 target_species="human")
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent=intent,
        )
        assert result.signature is None


# ---------------------------------------------------------------------------
# Direction Pattern
# ---------------------------------------------------------------------------

class TestEnhancedDirection:

    def test_diagnostic_builds_direction(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.direction_pattern is not None
        assert isinstance(result.direction_pattern, DirectionPattern)

    def test_direction_has_panel_genes(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        dp = result.direction_pattern
        # At least some panel genes should be in the pattern
        assert dp.n_genes > 0

    def test_monitoring_skips_direction(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="monitoring",
        )
        assert result.direction_pattern is None


# ---------------------------------------------------------------------------
# Clinical Metrics
# ---------------------------------------------------------------------------

class TestEnhancedClinical:

    def test_diagnostic_builds_clinical(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.clinical_metrics is not None
        assert "youdens" in result.clinical_metrics
        assert "bootstrap_ci" in result.clinical_metrics
        assert "ppv_npv" in result.clinical_metrics
        assert "decision_curve" in result.clinical_metrics

    def test_youdens_j_positive(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.clinical_metrics["youdens"]["youdens_j"] > 0.5

    def test_bootstrap_ci_contains_auc(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        ci = result.clinical_metrics["bootstrap_ci"]
        assert ci["ci_lower"] <= ci["point_estimate"] <= ci["ci_upper"]

    def test_ppv_at_custom_prevalence(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic", prevalence=0.10,
        )
        assert result.clinical_metrics["ppv_npv"]["prevalence"] == 0.10

    def test_exploratory_skips_clinical(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="exploratory",
        )
        assert result.clinical_metrics is None

    def test_no_model_skips_clinical(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        # Remove trained model
        clf = s["base_result"].classification_results["logistic_regression"]
        clf.trained_model = None
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        # Clinical metrics should still exist but be empty or None
        # because no y_score could be computed
        assert result.clinical_metrics is None or "youdens" not in result.clinical_metrics


# ---------------------------------------------------------------------------
# Ratio Biomarkers
# ---------------------------------------------------------------------------

class TestEnhancedRatios:

    def test_diagnostic_builds_ratios(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        assert result.ratio_result is not None
        assert isinstance(result.ratio_result, RatioSearchResult)

    def test_ratio_uses_panel_genes(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        # Pairs should only involve panel genes
        for pair in result.ratio_result.pairs:
            assert pair.gene_a in s["panel_genes"]
            assert pair.gene_b in s["panel_genes"]

    def test_monitoring_skips_ratios(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="monitoring",
        )
        assert result.ratio_result is None


# ---------------------------------------------------------------------------
# Summary & save
# ---------------------------------------------------------------------------

class TestEnhancedOutput:

    def test_summary_includes_enhanced(self, synthetic_setup: dict) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        summary = result.summary()
        assert "ENHANCED ANALYSIS" in summary
        assert "Diagnostic" in summary

    def test_save_creates_enhanced_dir(self, synthetic_setup: dict, tmp_path) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        result.save(tmp_path)
        enhanced_dir = tmp_path / "enhanced"
        assert enhanced_dir.exists()

    def test_save_creates_direction_csv(self, synthetic_setup: dict, tmp_path) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        result.save(tmp_path)
        assert (tmp_path / "enhanced" / "direction_pattern.csv").exists()

    def test_save_creates_ratio_csv(self, synthetic_setup: dict, tmp_path) -> None:
        s = synthetic_setup
        result = enhance_biomarker_result(
            s["base_result"], s["expression"], s["labels"],
            s["group_names"], intent="diagnostic",
        )
        result.save(tmp_path)
        assert (tmp_path / "enhanced" / "ratio_biomarkers.csv").exists()