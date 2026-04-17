"""Tests for ratio_biomarkers module."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery.ratio_biomarkers import (
    RatioBiomarkerSearcher,
    RatioPair,
    RatioSearchResult,
    apply_ratios,
    build_ratio_features,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_data() -> tuple:
    """
    60 samples × 10 genes.
    Genes g00/g01 form a discriminating ratio: g00 is UP in disease,
    g01 is DOWN, so g00/g01 is strongly elevated in disease.
    The rest are noise.
    """
    rng = np.random.default_rng(42)
    n_per = 30
    n_genes = 10

    base = rng.normal(10.0, 1.0, size=(n_per, n_genes))
    dis = rng.normal(10.0, 1.0, size=(n_per, n_genes))

    # Implant signal: g00 UP in disease, g01 DOWN in disease
    dis[:, 0] += 5.0   # g00 elevated
    dis[:, 1] -= 5.0   # g01 suppressed
    # Clip to prevent negatives (expression can't be negative)
    dis[:, 1] = np.clip(dis[:, 1], 0.1, None)

    X = np.vstack([base, dis])
    gene_names = [f"g{i:02d}" for i in range(n_genes)]
    expr = pd.DataFrame(X, columns=gene_names)
    labels = pd.Series(["healthy"] * n_per + ["disease"] * n_per)
    return expr, labels


@pytest.fixture
def simple_pair() -> RatioPair:
    return RatioPair(
        gene_a="gA", gene_b="gB", auc=0.92,
        direction="higher_in_reference",
        mean_ratio_reference=2.5,
        mean_ratio_baseline=0.8,
    )


# ---------------------------------------------------------------------------
# RatioPair
# ---------------------------------------------------------------------------

class TestRatioPair:

    def test_name_property(self, simple_pair: RatioPair) -> None:
        assert simple_pair.name == "gA/gB"

    def test_to_dict(self, simple_pair: RatioPair) -> None:
        d = simple_pair.to_dict()
        assert d["gene_a"] == "gA"
        assert d["gene_b"] == "gB"
        assert d["auc"] == 0.92
        assert d["name"] == "gA/gB"


# ---------------------------------------------------------------------------
# RatioBiomarkerSearcher — construction
# ---------------------------------------------------------------------------

class TestSearcherConstruction:

    def test_defaults(self) -> None:
        s = RatioBiomarkerSearcher()
        assert s.top_k == 20
        assert s.min_auc == 0.70
        assert s.pseudocount == 1.0

    def test_invalid_top_k_raises(self) -> None:
        with pytest.raises(ValueError, match="top_k"):
            RatioBiomarkerSearcher(top_k=0)

    def test_invalid_min_auc_raises(self) -> None:
        with pytest.raises(ValueError, match="min_auc"):
            RatioBiomarkerSearcher(min_auc=0.3)

    def test_invalid_pseudocount_raises(self) -> None:
        with pytest.raises(ValueError, match="pseudocount"):
            RatioBiomarkerSearcher(pseudocount=-1.0)


# ---------------------------------------------------------------------------
# RatioBiomarkerSearcher.search
# ---------------------------------------------------------------------------

class TestSearch:

    def test_finds_implanted_ratio(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=5, min_auc=0.80)
        result = searcher.search(expr, labels, "disease", "healthy")
        # g00/g01 should be the top pair (or close to it)
        top_genes = {(p.gene_a, p.gene_b) for p in result.pairs[:3]}
        assert ("g00", "g01") in top_genes or ("g01", "g00") in top_genes

    def test_best_pair_has_high_auc(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=5, min_auc=0.70)
        result = searcher.search(expr, labels, "disease", "healthy")
        assert result.best_pair is not None
        assert result.best_pair.auc > 0.90

    def test_result_metadata(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=5, min_auc=0.70)
        result = searcher.search(expr, labels, "disease", "healthy")
        assert result.reference_group == "disease"
        assert result.baseline_group == "healthy"
        assert result.n_genes_input == 10
        assert result.n_pairs_tested == 45  # C(10,2)
        assert result.top_k == 5

    def test_min_auc_filters_pairs(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        loose = RatioBiomarkerSearcher(top_k=100, min_auc=0.50)
        strict = RatioBiomarkerSearcher(top_k=100, min_auc=0.95)
        r_loose = loose.search(expr, labels, "disease", "healthy")
        r_strict = strict.search(expr, labels, "disease", "healthy")
        assert r_strict.n_found <= r_loose.n_found

    def test_top_k_limits_output(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=3, min_auc=0.50)
        result = searcher.search(expr, labels, "disease", "healthy")
        assert result.n_found <= 3

    def test_pairs_sorted_by_auc(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=10, min_auc=0.50)
        result = searcher.search(expr, labels, "disease", "healthy")
        aucs = [p.auc for p in result.pairs]
        assert aucs == sorted(aucs, reverse=True)

    def test_direction_is_valid(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=10, min_auc=0.50)
        result = searcher.search(expr, labels, "disease", "healthy")
        valid = {"higher_in_reference", "higher_in_baseline"}
        for p in result.pairs:
            assert p.direction in valid

    def test_all_aucs_above_threshold(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(top_k=100, min_auc=0.75)
        result = searcher.search(expr, labels, "disease", "healthy")
        for p in result.pairs:
            assert p.auc >= 0.75

    def test_invalid_group_raises(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher()
        with pytest.raises(ValueError, match="not in labels"):
            searcher.search(expr, labels, "martian", "healthy")

    def test_label_length_mismatch_raises(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher()
        with pytest.raises(ValueError, match="labels length"):
            searcher.search(expr, labels.iloc[:5], "disease", "healthy")

    def test_candidate_genes_restricts_search(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(
            top_k=10, min_auc=0.50,
            candidate_genes=["g00", "g01", "g02"],
        )
        result = searcher.search(expr, labels, "disease", "healthy")
        assert result.n_genes_input == 3
        assert result.n_pairs_tested == 3  # C(3,2)

    def test_candidate_genes_missing_warns(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(
            candidate_genes=["g00", "g01", "FAKE_GENE"],
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            searcher.search(expr, labels, "disease", "healthy")
            assert any("not found" in str(w.message).lower() for w in caught)

    def test_too_few_candidates_raises(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        searcher = RatioBiomarkerSearcher(
            candidate_genes=["g00", "FAKE1", "FAKE2"],
        )
        with pytest.raises(ValueError, match="at least 2"):
            searcher.search(expr, labels, "disease", "healthy")


# ---------------------------------------------------------------------------
# RatioSearchResult
# ---------------------------------------------------------------------------

class TestSearchResult:

    def test_to_dataframe(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result = RatioBiomarkerSearcher(top_k=5, min_auc=0.50).search(
            expr, labels, "disease", "healthy"
        )
        df = result.to_dataframe()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == result.n_found
        assert "auc" in df.columns
        assert "gene_a" in df.columns

    def test_to_dataframe_empty(self) -> None:
        result = RatioSearchResult(
            pairs=[], reference_group="A", baseline_group="B",
            n_pairs_tested=0, n_genes_input=0, top_k=10, min_auc=0.99,
        )
        df = result.to_dataframe()
        assert len(df) == 0
        assert "auc" in df.columns

    def test_summary_string(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result = RatioBiomarkerSearcher(top_k=5, min_auc=0.50).search(
            expr, labels, "disease", "healthy"
        )
        s = result.summary()
        assert "Ratio Biomarker" in s
        assert "disease" in s

    def test_best_pair_none_when_empty(self) -> None:
        result = RatioSearchResult(
            pairs=[], reference_group="A", baseline_group="B",
            n_pairs_tested=0, n_genes_input=0, top_k=10, min_auc=0.99,
        )
        assert result.best_pair is None


# ---------------------------------------------------------------------------
# apply_ratios
# ---------------------------------------------------------------------------

class TestApplyRatios:

    def test_basic_application(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result = RatioBiomarkerSearcher(top_k=3, min_auc=0.70).search(
            expr, labels, "disease", "healthy"
        )
        ratio_df = apply_ratios(expr, result.pairs)
        assert isinstance(ratio_df, pd.DataFrame)
        assert len(ratio_df) == len(expr)
        assert ratio_df.shape[1] == result.n_found

    def test_column_names_match_pairs(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result = RatioBiomarkerSearcher(top_k=3, min_auc=0.70).search(
            expr, labels, "disease", "healthy"
        )
        ratio_df = apply_ratios(expr, result.pairs)
        expected_names = [p.name for p in result.pairs]
        assert list(ratio_df.columns) == expected_names

    def test_ratios_are_positive(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result = RatioBiomarkerSearcher(top_k=3, min_auc=0.70).search(
            expr, labels, "disease", "healthy"
        )
        ratio_df = apply_ratios(expr, result.pairs)
        assert (ratio_df > 0).all().all()

    def test_missing_gene_raises(self, simple_pair: RatioPair) -> None:
        expr = pd.DataFrame({"gA": [1.0, 2.0], "gC": [3.0, 4.0]})
        with pytest.raises(ValueError, match="missing genes"):
            apply_ratios(expr, [simple_pair])

    def test_empty_pairs_raises(self) -> None:
        expr = pd.DataFrame({"g1": [1.0]})
        with pytest.raises(ValueError, match="No ratio pairs"):
            apply_ratios(expr, [])

    def test_pseudocount_zero(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result = RatioBiomarkerSearcher(
            top_k=3, min_auc=0.70, pseudocount=0.0
        ).search(expr, labels, "disease", "healthy")
        # Should still work since our synthetic data has no zeros
        ratio_df = apply_ratios(expr, result.pairs, pseudocount=0.0)
        assert ratio_df.shape[1] == result.n_found


# ---------------------------------------------------------------------------
# build_ratio_features
# ---------------------------------------------------------------------------

class TestBuildRatioFeatures:

    def test_returns_result_and_dataframe(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result, ratio_df = build_ratio_features(
            expr, labels, "disease", "healthy", top_k=5, min_auc=0.70,
        )
        assert isinstance(result, RatioSearchResult)
        assert isinstance(ratio_df, pd.DataFrame)
        assert len(ratio_df) == len(expr)

    def test_no_pairs_found_raises(self) -> None:
        rng = np.random.default_rng(42)
        # Pure noise — no pair can discriminate
        expr = pd.DataFrame(rng.normal(10, 1, size=(40, 5)),
                            columns=[f"g{i}" for i in range(5)])
        labels = pd.Series(["ctrl"] * 20 + ["case"] * 20)
        with pytest.raises(ValueError, match="No ratio pairs found"):
            build_ratio_features(
                expr, labels, "case", "ctrl",
                min_auc=0.99,
            )

    def test_candidate_genes_passthrough(self, synthetic_data: tuple) -> None:
        expr, labels = synthetic_data
        result, ratio_df = build_ratio_features(
            expr, labels, "disease", "healthy",
            top_k=5, min_auc=0.50,
            candidate_genes=["g00", "g01", "g02", "g03"],
        )
        assert result.n_genes_input == 4


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:

    def test_two_genes_only(self) -> None:
        rng = np.random.default_rng(99)
        expr = pd.DataFrame({
            "gA": np.concatenate([rng.normal(5, 1, 20), rng.normal(10, 1, 20)]),
            "gB": np.concatenate([rng.normal(10, 1, 20), rng.normal(5, 1, 20)]),
        })
        labels = pd.Series(["ctrl"] * 20 + ["case"] * 20)
        result = RatioBiomarkerSearcher(top_k=5, min_auc=0.50).search(
            expr, labels, "case", "ctrl"
        )
        assert result.n_pairs_tested == 1
        assert result.n_found >= 1

    def test_constant_gene_filtered(self) -> None:
        rng = np.random.default_rng(99)
        expr = pd.DataFrame({
            "gA": np.concatenate([rng.normal(5, 1, 20), rng.normal(10, 1, 20)]),
            "gB": np.concatenate([rng.normal(10, 1, 20), rng.normal(5, 1, 20)]),
            "gConst": np.full(40, 7.0),  # zero variance
        })
        labels = pd.Series(["ctrl"] * 20 + ["case"] * 20)
        result = RatioBiomarkerSearcher(top_k=5, min_auc=0.50).search(
            expr, labels, "case", "ctrl"
        )
        # gConst should be filtered out, leaving C(2,2)=1 pair
        assert result.n_genes_input == 2

    def test_pseudocount_prevents_division_by_zero(self) -> None:
        expr = pd.DataFrame({
            "gA": [0.0, 0.0, 5.0, 5.0, 10.0, 10.0],
            "gB": [5.0, 5.0, 0.0, 0.0, 0.0, 0.0],
        })
        labels = pd.Series(["ctrl", "ctrl", "ctrl", "case", "case", "case"])
        searcher = RatioBiomarkerSearcher(top_k=5, min_auc=0.50, pseudocount=1.0)
        # Should not crash
        result = searcher.search(expr, labels, "case", "ctrl")
        assert result.n_pairs_tested == 1