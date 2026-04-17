"""Tests for direction_patterns module."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery.direction_patterns import (
    DIRECTION_DOWN,
    DIRECTION_UP,
    VALID_DIRECTIONS,
    DirectionPattern,
    build_direction_pattern,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_pattern() -> DirectionPattern:
    """A small hand-built pattern: geneA UP, geneB DOWN, geneC UP."""
    return DirectionPattern(
        gene_directions={"geneA": "UP", "geneB": "DOWN", "geneC": "UP"},
        fold_changes={"geneA": 2.0, "geneB": -1.5, "geneC": 1.2},
        confidence={"geneA": 5.0, "geneB": 3.0, "geneC": 4.0},
        reference_group="disease",
        baseline_group="healthy",
        baseline_means={"geneA": 5.0, "geneB": 10.0, "geneC": 3.0},
        baseline_stds={"geneA": 1.0, "geneB": 1.0, "geneC": 1.0},
        n_samples=40,
    )


@pytest.fixture
def synthetic_cohort() -> tuple[pd.DataFrame, pd.Series]:
    """
    60 samples x 20 genes. Genes 0-4 are UP in disease, 5-9 are DOWN,
    10-19 are noise. Returned as (expression, labels).
    """
    rng = np.random.default_rng(42)
    n_per_group = 30
    n_genes = 20

    # Baseline expression ~ N(5, 1)
    base = rng.normal(5.0, 1.0, size=(n_per_group, n_genes))
    dis = rng.normal(5.0, 1.0, size=(n_per_group, n_genes))
    # Strong UP in disease for genes 0-4
    dis[:, 0:5] += 3.0
    # Strong DOWN in disease for genes 5-9
    dis[:, 5:10] -= 3.0
    # Genes 10-19 stay as noise

    X = np.vstack([base, dis])
    gene_names = [f"g{i:02d}" for i in range(n_genes)]
    expr = pd.DataFrame(X, columns=gene_names)
    labels = pd.Series(["healthy"] * n_per_group + ["disease"] * n_per_group)
    return expr, labels


# ---------------------------------------------------------------------------
# Construction / validation
# ---------------------------------------------------------------------------

def test_valid_construction(simple_pattern: DirectionPattern) -> None:
    assert simple_pattern.n_genes == 3
    assert simple_pattern.reference_group == "disease"
    assert simple_pattern.test_used == "t-test"


def test_invalid_direction_raises() -> None:
    with pytest.raises(ValueError, match="Invalid direction"):
        DirectionPattern(
            gene_directions={"g1": "SIDEWAYS"},
            fold_changes={"g1": 1.0},
            confidence={"g1": 2.0},
            reference_group="A", baseline_group="B",
            baseline_means={"g1": 0.0}, baseline_stds={"g1": 1.0},
            n_samples=10,
        )


def test_gene_set_mismatch_raises() -> None:
    with pytest.raises(ValueError, match="Gene set mismatch"):
        DirectionPattern(
            gene_directions={"g1": "UP", "g2": "DOWN"},
            fold_changes={"g1": 1.0},  # missing g2
            confidence={"g1": 2.0, "g2": 3.0},
            reference_group="A", baseline_group="B",
            baseline_means={"g1": 0.0, "g2": 0.0},
            baseline_stds={"g1": 1.0, "g2": 1.0},
            n_samples=10,
        )


def test_valid_directions_constant() -> None:
    assert DIRECTION_UP in VALID_DIRECTIONS
    assert DIRECTION_DOWN in VALID_DIRECTIONS
    assert len(VALID_DIRECTIONS) == 2


# ---------------------------------------------------------------------------
# Properties
# ---------------------------------------------------------------------------

def test_up_down_counts(simple_pattern: DirectionPattern) -> None:
    assert simple_pattern.n_up == 2
    assert simple_pattern.n_down == 1
    assert simple_pattern.n_up + simple_pattern.n_down == simple_pattern.n_genes


def test_genes_property(simple_pattern: DirectionPattern) -> None:
    assert set(simple_pattern.genes) == {"geneA", "geneB", "geneC"}


# ---------------------------------------------------------------------------
# concordance
# ---------------------------------------------------------------------------

def test_concordance_perfect_match(simple_pattern: DirectionPattern) -> None:
    # UP genes (A, C) elevated; DOWN gene (B) suppressed
    patient = pd.Series({"geneA": 10.0, "geneB": 5.0, "geneC": 8.0})
    score = simple_pattern.concordance(patient)
    assert score > 0.9


def test_concordance_perfect_inversion(simple_pattern: DirectionPattern) -> None:
    # UP genes low, DOWN gene high -> opposite of pattern
    patient = pd.Series({"geneA": 0.0, "geneB": 15.0, "geneC": -2.0})
    score = simple_pattern.concordance(patient)
    assert score < -0.9


def test_concordance_baseline_like(simple_pattern: DirectionPattern) -> None:
    # Patient sits exactly at baseline means -> score near 0
    patient = pd.Series({"geneA": 5.0, "geneB": 10.0, "geneC": 3.0})
    score = simple_pattern.concordance(patient)
    assert abs(score) < 1e-6


def test_concordance_accepts_dict(simple_pattern: DirectionPattern) -> None:
    score = simple_pattern.concordance({"geneA": 10.0, "geneB": 5.0, "geneC": 8.0})
    assert isinstance(score, float)
    assert score > 0.9


def test_concordance_dataframe_returns_series(simple_pattern: DirectionPattern) -> None:
    patients = pd.DataFrame({
        "geneA": [10.0, 0.0, 5.0],
        "geneB": [5.0, 15.0, 10.0],
        "geneC": [8.0, -2.0, 3.0],
    })
    scores = simple_pattern.concordance(patients)
    assert isinstance(scores, pd.Series)
    assert len(scores) == 3
    assert scores.iloc[0] > 0.9    # match
    assert scores.iloc[1] < -0.9   # inverted
    assert abs(scores.iloc[2]) < 1e-6  # baseline


def test_concordance_missing_genes_warns(simple_pattern: DirectionPattern) -> None:
    patient = pd.Series({"geneA": 10.0, "geneC": 8.0})  # geneB missing
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        score = simple_pattern.concordance(patient)
        assert any("missing" in str(w.message).lower() for w in caught)
    assert score > 0.9   # still computes on available genes


def test_concordance_no_overlap_raises(simple_pattern: DirectionPattern) -> None:
    patient = pd.Series({"otherGene1": 1.0, "otherGene2": 2.0})
    with pytest.raises(ValueError, match="No genes.*match this pattern"):
        simple_pattern.concordance(patient)


def test_concordance_unweighted_option(simple_pattern: DirectionPattern) -> None:
    # A patient that matches UP genes strongly but slightly opposes DOWN.
    # With weighting (geneA highest confidence=5), weighted > unweighted.
    patient = pd.Series({"geneA": 100.0, "geneB": 11.0, "geneC": 4.0})
    w = simple_pattern.concordance(patient, weighted=True)
    u = simple_pattern.concordance(patient, weighted=False)
    assert w != u
    # weighted should lean more positive because geneA (weight 5) dominates
    assert w > u


def test_concordance_zero_std_does_not_crash() -> None:
    p = DirectionPattern(
        gene_directions={"g1": "UP"},
        fold_changes={"g1": 2.0},
        confidence={"g1": 3.0},
        reference_group="R", baseline_group="B",
        baseline_means={"g1": 5.0},
        baseline_stds={"g1": 0.0},  # pathological
        n_samples=10,
    )
    score = p.concordance(pd.Series({"g1": 10.0}))
    assert np.isfinite(score)
    assert score > 0   # still in the right direction


# ---------------------------------------------------------------------------
# cross_cohort_check
# ---------------------------------------------------------------------------

def test_cross_cohort_perfect_agreement(simple_pattern: DirectionPattern) -> None:
    twin = DirectionPattern(
        gene_directions=dict(simple_pattern.gene_directions),
        fold_changes=dict(simple_pattern.fold_changes),
        confidence=dict(simple_pattern.confidence),
        reference_group="disease", baseline_group="healthy",
        baseline_means=dict(simple_pattern.baseline_means),
        baseline_stds=dict(simple_pattern.baseline_stds),
        n_samples=40,
    )
    report = simple_pattern.cross_cohort_check(twin)
    assert report["n_common"] == 3
    assert report["n_agree"] == 3
    assert report["n_disagree"] == 0
    assert report["agreement_fraction"] == 1.0
    assert report["disagreements"] == []


def test_cross_cohort_perfect_disagreement(simple_pattern: DirectionPattern) -> None:
    flipped = DirectionPattern(
        gene_directions={"geneA": "DOWN", "geneB": "UP", "geneC": "DOWN"},
        fold_changes={"geneA": -2.0, "geneB": 1.5, "geneC": -1.2},
        confidence={"geneA": 5.0, "geneB": 3.0, "geneC": 4.0},
        reference_group="disease", baseline_group="healthy",
        baseline_means={"geneA": 5.0, "geneB": 10.0, "geneC": 3.0},
        baseline_stds={"geneA": 1.0, "geneB": 1.0, "geneC": 1.0},
        n_samples=40,
    )
    report = simple_pattern.cross_cohort_check(flipped)
    assert report["n_agree"] == 0
    assert report["n_disagree"] == 3
    assert report["agreement_fraction"] == 0.0
    # Highly unlikely under chance (p=0.5), even with n=3
    assert report["binomial_p"] < 0.3


def test_cross_cohort_partial_agreement(simple_pattern: DirectionPattern) -> None:
    partial = DirectionPattern(
        gene_directions={"geneA": "UP", "geneB": "UP", "geneC": "UP"},  # B flipped
        fold_changes={"geneA": 2.0, "geneB": 1.0, "geneC": 1.2},
        confidence={"geneA": 5.0, "geneB": 3.0, "geneC": 4.0},
        reference_group="disease", baseline_group="healthy",
        baseline_means={"geneA": 5.0, "geneB": 10.0, "geneC": 3.0},
        baseline_stds={"geneA": 1.0, "geneB": 1.0, "geneC": 1.0},
        n_samples=40,
    )
    report = simple_pattern.cross_cohort_check(partial)
    assert report["n_agree"] == 2
    assert report["n_disagree"] == 1
    assert report["agreements"] == ["geneA", "geneC"]
    assert report["disagreements"] == [("geneB", "DOWN", "UP")]


def test_cross_cohort_disjoint_gene_sets(simple_pattern: DirectionPattern) -> None:
    disjoint = DirectionPattern(
        gene_directions={"geneX": "UP"},
        fold_changes={"geneX": 2.0},
        confidence={"geneX": 3.0},
        reference_group="disease", baseline_group="healthy",
        baseline_means={"geneX": 1.0}, baseline_stds={"geneX": 1.0},
        n_samples=10,
    )
    report = simple_pattern.cross_cohort_check(disjoint)
    assert report["n_common"] == 0
    assert np.isnan(report["agreement_fraction"])
    assert np.isnan(report["binomial_p"])
    assert set(report["only_in_this"]) == {"geneA", "geneB", "geneC"}
    assert report["only_in_other"] == ["geneX"]


def test_cross_cohort_rejects_non_pattern(simple_pattern: DirectionPattern) -> None:
    with pytest.raises(TypeError):
        simple_pattern.cross_cohort_check({"not": "a pattern"})


# ---------------------------------------------------------------------------
# build_direction_pattern
# ---------------------------------------------------------------------------

def test_build_detects_up_down_signal(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    pattern = build_direction_pattern(
        expr, labels,
        reference_group="disease", baseline_group="healthy",
    )
    # Genes 0-4 were implanted UP, 5-9 DOWN, 10-19 noise
    for i in range(5):
        assert pattern.gene_directions[f"g{i:02d}"] == "UP"
    for i in range(5, 10):
        assert pattern.gene_directions[f"g{i:02d}"] == "DOWN"
    # Noise genes mostly shouldn't appear
    noise_found = [g for g in pattern.genes if g[1:].isdigit() and int(g[1:]) >= 10]
    assert len(noise_found) < 5


def test_build_respects_p_threshold(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    loose = build_direction_pattern(expr, labels, "disease", "healthy",
                                    p_threshold=0.5, fc_threshold=0.0)
    strict = build_direction_pattern(expr, labels, "disease", "healthy",
                                     p_threshold=1e-10, fc_threshold=0.0)
    assert strict.n_genes < loose.n_genes


def test_build_respects_fc_threshold(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    loose = build_direction_pattern(expr, labels, "disease", "healthy",
                                    p_threshold=0.05, fc_threshold=0.1)
    strict = build_direction_pattern(expr, labels, "disease", "healthy",
                                     p_threshold=0.05, fc_threshold=3.3)
    assert strict.n_genes < loose.n_genes


def test_build_invalid_group_raises(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    with pytest.raises(ValueError, match="not found in labels"):
        build_direction_pattern(expr, labels, "martian", "healthy")


def test_build_label_length_mismatch_raises(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    with pytest.raises(ValueError, match="does not match"):
        build_direction_pattern(expr, labels.iloc[:10],
                                "disease", "healthy")


def test_build_too_few_samples_raises() -> None:
    expr = pd.DataFrame({"g1": [1.0, 2.0, 3.0], "g2": [4.0, 5.0, 6.0]})
    labels = pd.Series(["A", "A", "B"])  # only 1 sample in group B
    with pytest.raises(ValueError, match="at least 2 samples"):
        build_direction_pattern(expr, labels, "A", "B")


def test_build_invalid_test_raises(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    with pytest.raises(ValueError, match="'t-test' or 'mann-whitney'"):
        build_direction_pattern(expr, labels, "disease", "healthy",
                                test="chi-square")


def test_build_mann_whitney_works(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    pattern = build_direction_pattern(
        expr, labels, "disease", "healthy", test="mann-whitney",
    )
    assert pattern.test_used == "mann-whitney"
    assert pattern.n_genes >= 8   # should still find the implanted signals


def test_build_no_genes_passing_raises(synthetic_cohort: tuple) -> None:
    expr, labels = synthetic_cohort
    with pytest.raises(ValueError, match="No genes passed thresholds"):
        build_direction_pattern(expr, labels, "disease", "healthy",
                                p_threshold=1e-300, fc_threshold=100.0)


# ---------------------------------------------------------------------------
# to_dataframe
# ---------------------------------------------------------------------------

def test_to_dataframe_shape(simple_pattern: DirectionPattern) -> None:
    df = simple_pattern.to_dataframe()
    assert set(df.columns) == {
        "direction", "log2FC", "neg_log10_p", "baseline_mean", "baseline_std"
    }
    assert set(df.index) == {"geneA", "geneB", "geneC"}
    assert df.loc["geneA", "direction"] == "UP"