"""
Tests for M2: bootstrap CI semantics for clinical metrics.

M2 is a reporting-layer fix with three touch-points:

    1. enhanced.py populates clin['optimism_gap'] with the correct
       pair of quantities (pooled OOF AUC vs true full-data training
       AUC). M2 did NOT change this; M1 introduced it. We assert on
       this structure so a future regression cannot silently remove
       it.

    2. The dashboard's _resolve_optimism_tiles helper routes reads
       from the optimism_gap struct when available, and falls back to
       the bootstrap-point-estimate approximation only when the
       optimism_gap is absent (pre-M1 pickled results). We test the
       resolution logic as a pure function, no Streamlit runtime
       needed.

    3. EnhancedBiomarkerResult.summary() labels the bootstrap CI line
       honestly: "OOF AUC bootstrap CI" when oof_used, else
       "Training AUC bootstrap CI". Symmetric, reader-friendly.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from types import ModuleType, SimpleNamespace
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import pytest


# =============================================================================
# Load dashboard helpers without triggering Streamlit runtime.
# Same approach as test_dashboard_panel_stability.py.
# =============================================================================


def _load_dashboard_helpers() -> ModuleType:
    here = Path(__file__).resolve().parent.parent
    dashboard_path = (
        here / "raptor" / "dashboard" / "pages"
        / "08_\U0001f9ec_Biomarker_Discovery.py"
    )
    if not dashboard_path.exists():
        pages_dir = here / "raptor" / "dashboard" / "pages"
        candidates = list(pages_dir.glob("08_*Biomarker_Discovery*.py"))
        if candidates:
            dashboard_path = candidates[0]
        else:
            raise FileNotFoundError(
                f"Could not find 08_*Biomarker_Discovery*.py under "
                f"{pages_dir}."
            )
    source = dashboard_path.read_text(encoding="utf-8")
    marker = "# PAGE CONFIG"
    idx = source.find(marker)
    helper_marker = "# Dashboard polish helpers"
    helper_idx = source.find(helper_marker)
    helpers_only = source[helper_idx:idx]

    mod = ModuleType("dashboard_helpers_under_test_m2")

    class _NoOp:
        def __getattr__(self, _):
            return _NoOp()
        def __call__(self, *args, **kwargs):
            return _NoOp()

    mod.__dict__["st"] = _NoOp()
    mod.__dict__["px"] = _NoOp()
    mod.__dict__["go"] = _NoOp()
    mod.__dict__["np"] = np
    mod.__dict__["pd"] = pd
    from pathlib import Path as _Path
    from typing import Dict, List, Optional, Tuple
    mod.__dict__.update({
        "Path": _Path, "Dict": Dict, "List": List,
        "Optional": Optional, "Tuple": Tuple,
        "io": __import__("io"), "sys": __import__("sys"),
    })
    exec(compile(helpers_only, str(dashboard_path), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def dash():
    return _load_dashboard_helpers()


# =============================================================================
# Shared small-fixture for end-to-end M1 contract tests
# =============================================================================


@pytest.fixture
def signal_discovery_fixture():
    """40 samples x 200 genes with strong 2.5 sigma signal in 10 genes.
    Small enough to run in ~30s under nested_cv."""
    rng = np.random.default_rng(42)
    n_samples, n_genes, n_signal = 40, 200, 10
    y_label = ['disease'] * 20 + ['healthy'] * 20
    X = rng.standard_normal((n_samples, n_genes)) * 0.8
    X[:20, :n_signal] += 2.5
    sample_ids = [f'S{i:03d}' for i in range(n_samples)]
    gene_ids = [f'GENE_{i:04d}' for i in range(n_genes)]
    counts = pd.DataFrame(
        np.clip(X.T * 100 + 500, 0, None).astype(int),
        index=gene_ids, columns=sample_ids,
    )
    metadata = pd.DataFrame(
        {'sample_id': sample_ids, 'condition': y_label}
    ).set_index('sample_id')
    return counts, metadata


# =============================================================================
# Group 8.1 — optimism_gap populated correctly under M1 (3 tests)
# =============================================================================


class TestClinOptimismGapPopulatedUnderM1:
    """enhanced.py's optimism_gap struct must carry the correct pair.

    These tests defend the M1+ contract: when nested_cv runs,
    clinical_metrics.optimism_gap has cv_auc_oof, training_auc, gap;
    training_auc comes from a full-data fit (hence on strong-signal
    fixtures it's >= the pooled OOF AUC).
    """

    def test_optimism_gap_present_with_oof(self, signal_discovery_fixture):
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata = signal_discovery_fixture
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                intent='diagnostic',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m2_test_gap_present',
                random_state=42, verbose=False,
                calibrate_consensus=True, alpha=0.05,
            )
        assert r.clinical_metrics is not None, \
            "enhanced.py must populate clinical_metrics under diagnostic intent"
        og = r.clinical_metrics.get('optimism_gap')
        assert og is not None, \
            "optimism_gap key must be present under nested_cv + OOF path"
        assert 'cv_auc_oof' in og
        assert 'training_auc' in og
        assert 'gap' in og

    def test_optimism_gap_training_auc_is_from_full_data_fit(
        self, signal_discovery_fixture
    ):
        """On strong signal, training AUC >= pooled OOF AUC.

        If training_auc were actually a mislabeled OOF estimate (the
        pre-M2 dashboard bug), this invariant would not hold
        systematically. By requiring training_auc >= cv_auc_oof on a
        fixture we know separates cleanly, we confirm training_auc is
        from a fit that saw every sample.
        """
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata = signal_discovery_fixture
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                intent='diagnostic',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m2_test_train_gte_oof',
                random_state=42, verbose=False,
            )
        og = r.clinical_metrics['optimism_gap']
        # Allow a tiny epsilon because on perfectly-separable signal
        # fixtures both training and OOF AUCs can saturate at 1.0.
        assert og['training_auc'] >= og['cv_auc_oof'] - 1e-9, (
            f"training_auc ({og['training_auc']:.4f}) < "
            f"cv_auc_oof ({og['cv_auc_oof']:.4f}) on a strong-signal "
            f"fixture. Either training_auc is mislabeled (came from "
            f"OOF, not full-data fit) or the test fixture has drifted."
        )
        # Gap matches difference
        assert abs(og['gap'] - (og['training_auc'] - og['cv_auc_oof'])) < 1e-9

    def test_optimism_gap_skipped_under_loocv(self, signal_discovery_fixture):
        """LOOCV path does not populate optimism_gap.

        enhanced.py gates the optimism_gap block on both oof_used and
        trained_model being present. LOOCV produces OOF preds but
        does NOT preserve trained_model. So optimism_gap is absent,
        and M2's dashboard helper falls through to the 'fallback' or
        'unavailable' branch.

        Note: as of M6, validation='loocv' emits a DeprecationWarning.
        We catch warnings so the test doesn't fail on that.
        """
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata = signal_discovery_fixture
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                intent='diagnostic',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='loocv',
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m2_test_loocv_no_optgap',
                random_state=42, verbose=False,
            )
        # Clinical metrics might still exist (bootstrap_ci, youdens...)
        # but optimism_gap should be absent because trained_model isn't
        # preserved under LOOCV.
        if r.clinical_metrics is not None:
            assert 'optimism_gap' not in r.clinical_metrics, (
                "optimism_gap should not be populated under LOOCV "
                "(trained_model not preserved). If this fails, check "
                "enhanced.py's gating at the optimism_gap block."
            )


# =============================================================================
# Group 8.2 — dashboard _resolve_optimism_tiles helper (3 tests)
# =============================================================================


class TestResolveOptimismTiles:
    """_resolve_optimism_tiles is the single place the dashboard reads
    optimism-diagnostic values. Test all three branches."""

    def _make_base(self, best_clf_name, best_auc):
        """Synthesize a minimal BiomarkerResult-like object for testing."""
        clf_result = SimpleNamespace(auc=best_auc)
        return SimpleNamespace(
            best_classifier=best_clf_name,
            classification_results={best_clf_name: clf_result},
        )

    def test_reads_optimism_gap_when_available(self, dash):
        """M1+ path: clin has optimism_gap -> use it verbatim."""
        clin = {
            "optimism_gap": {
                "cv_auc_oof": 0.85,
                "training_auc": 0.92,
                "gap": 0.07,
            },
            # Also has bootstrap_ci, but optimism_gap takes priority
            "bootstrap_ci": {"point_estimate": 0.88},
        }
        base = self._make_base("logistic_regression", 0.84)
        cv, tr, gap, source = dash._resolve_optimism_tiles(clin, base)
        assert source == "optimism_gap"
        assert cv == 0.85
        assert tr == 0.92
        assert gap == 0.07
        # NOT the classification_results auc (0.84) or bootstrap point (0.88)
        assert cv != 0.84
        assert tr != 0.88

    def test_falls_back_to_bootstrap_when_optimism_gap_absent(self, dash):
        """Pre-M1 pickled result: no optimism_gap key. Use
        (classification_results[best].auc, bootstrap_ci.point_estimate)
        and compute gap as their difference."""
        clin = {
            "bootstrap_ci": {"point_estimate": 0.88},
            # explicitly NO optimism_gap key
        }
        base = self._make_base("logistic_regression", 0.84)
        cv, tr, gap, source = dash._resolve_optimism_tiles(clin, base)
        assert source == "fallback"
        assert cv == 0.84
        assert tr == 0.88
        assert abs(gap - 0.04) < 1e-9

    def test_returns_unavailable_when_nothing_extractable(self, dash):
        """Neither optimism_gap nor (best_classifier + bootstrap_ci)
        present -> return all-None with source='unavailable'. Dashboard
        uses this to skip the optimism-diagnostic block entirely.
        """
        # Empty clin
        cv, tr, gap, source = dash._resolve_optimism_tiles({}, self._make_base(None, 0.0))
        assert source == "unavailable"
        assert cv is None and tr is None and gap is None

        # clin has bootstrap_ci but no best_classifier set on base
        clin = {"bootstrap_ci": {"point_estimate": 0.88}}
        base_no_best = SimpleNamespace(
            best_classifier=None, classification_results={},
        )
        cv2, tr2, gap2, source2 = dash._resolve_optimism_tiles(clin, base_no_best)
        assert source2 == "unavailable"
        assert cv2 is None


# =============================================================================
# Group 8.3 — EnhancedBiomarkerResult.summary() label (2 tests)
# =============================================================================


class TestSummaryLabelOofAuc:
    """summary() bootstrap-CI line labels the AUC honestly."""

    def test_summary_says_oof_auc_under_m1(self, signal_discovery_fixture):
        """When oof_used=True, summary contains 'OOF AUC bootstrap CI'."""
        from raptor.biomarker_discovery import discover_biomarkers
        counts, metadata = signal_discovery_fixture
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                intent='diagnostic',
                methods=['univariate_filter'],
                min_panel=5, max_panel=5, target_panel_size=5,
                validation='nested_cv', n_folds=5, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/m2_test_summary_oof',
                random_state=42, verbose=False,
            )
        assert r.clinical_metrics.get('oof_used') is True
        text = r.summary()
        assert "OOF AUC bootstrap CI" in text, (
            "Under M1 OOF path, summary() must label the bootstrap CI "
            "as 'OOF AUC bootstrap CI' (not plain 'AUC bootstrap CI' "
            "which is ambiguous)."
        )

    def test_summary_label_mechanics(self, dash):
        """Construct an EnhancedBiomarkerResult directly and verify the
        conditional label logic. Faster than the end-to-end version."""
        from raptor.biomarker_discovery import enhance_biomarker_result  # noqa: F401
        # Build a minimal EnhancedBiomarkerResult. We don't need the full
        # discovery pipeline for this test -- just the summary()
        # formatting branch. The enhanced module's label logic reads
        # cm.get("oof_used", False) and switches the label. We verify
        # both branches by constructing clin dicts with and without
        # oof_used=True and checking the output text.
        from raptor.biomarker_discovery.enhanced import EnhancedBiomarkerResult

        # Minimal base-result-like object that EnhancedBiomarkerResult
        # can wrap. Only the summary() method is exercised, so only
        # summary()-relevant fields need to be mock-friendly.
        class _MockBase:
            def summary(self):
                return "MOCK BASE SUMMARY"

        class _MockIntent:
            output_label = "Diagnostic"
            intent = "diagnostic"

        oof_clin = {
            "oof_used": True,
            "bootstrap_ci": {
                "point_estimate": 0.92,
                "ci_lower": 0.85, "ci_upper": 0.98,
                "ci_level": 0.95, "n_bootstrap": 2000,
            },
        }
        train_clin = {
            "oof_used": False,
            "bootstrap_ci": {
                "point_estimate": 0.95,
                "ci_lower": 0.88, "ci_upper": 0.99,
                "ci_level": 0.95, "n_bootstrap": 2000,
            },
        }

        r_oof = EnhancedBiomarkerResult(
            base_result=_MockBase(),
            intent=_MockIntent(),
            signature=None, direction_pattern=None,
            clinical_metrics=oof_clin, ratio_result=None,
        )
        r_train = EnhancedBiomarkerResult(
            base_result=_MockBase(),
            intent=_MockIntent(),
            signature=None, direction_pattern=None,
            clinical_metrics=train_clin, ratio_result=None,
        )

        assert "OOF AUC bootstrap CI" in r_oof.summary()
        assert "Training AUC bootstrap CI" in r_train.summary()
        # And the labels don't bleed across
        assert "Training AUC bootstrap CI" not in r_oof.summary()
        assert "OOF AUC bootstrap CI" not in r_train.summary()


# =============================================================================
# Group 8.4 — regression: stale LOOCV auto-override banner must stay removed
# =============================================================================


class TestStaleLOOCVBannerRemoved:
    """Guardrail against accidental re-introduction of a dead banner.

    Context: pre-M6 core.py had an n<20 LOOCV auto-override. The
    dashboard had a banner that fired on every small-cohort run to
    explain that auto-switch. M6 removed the auto-override but the
    banner remained, firing on every small-cohort nested_cv run and
    contradicting the populated clinical_metrics shown immediately
    below it (e.g. 10a at n=10 showed the banner + full clinical
    metrics + optimism_gap simultaneously, which would have been
    confusing in paper figures).

    The banner was removed in the M2 cleanup. This test pins that
    removal in place: the specific banner strings must not re-appear
    in the dashboard source. If a future small-n UX feature needs a
    banner, it should use a different trigger (e.g. classifier
    collapse detection) and different wording — not the stale
    n_samples < 20 check whose original trigger no longer exists.
    """

    def _dashboard_source(self) -> str:
        """Read the dashboard page source as text."""
        here = Path(__file__).resolve().parent.parent
        dashboard_path = (
            here / "raptor" / "dashboard" / "pages"
            / "08_\U0001f9ec_Biomarker_Discovery.py"
        )
        if not dashboard_path.exists():
            pages_dir = here / "raptor" / "dashboard" / "pages"
            candidates = list(pages_dir.glob("08_*Biomarker_Discovery*.py"))
            assert candidates, (
                f"Could not find 08_*Biomarker_Discovery*.py under {pages_dir}"
            )
            dashboard_path = candidates[0]
        return dashboard_path.read_text(encoding="utf-8")

    def test_no_auto_override_banner_string(self):
        """The "(auto-override)" parenthetical must not reappear."""
        source = self._dashboard_source()
        assert "LOOCV used (auto-override)" not in source, (
            "The stale LOOCV auto-override banner string is present in "
            "the dashboard source. This banner was removed in April "
            "2026 because M6 eliminated the n<20 auto-override it was "
            "trying to explain. If you genuinely need a small-cohort "
            "advisory, use a different trigger and wording."
        )

    def test_no_loocv_used_branch_string(self):
        """The else-branch "**LOOCV used**" wording must not reappear."""
        source = self._dashboard_source()
        # Look for the exact markdown-bold form from the removed else
        # branch. A generic "LOOCV used" could plausibly show up in
        # legitimate future documentation (e.g. in the Reports tab
        # describing what the user selected), so we pin the stronger
        # markdown-bold + trailing "for n=" signature.
        assert "**LOOCV used** for n=" not in source, (
            "The stale 'LOOCV used for n=...' banner string is present "
            "in the dashboard source. Removed April 2026; see "
            "TestStaleLOOCVBannerRemoved docstring for rationale."
        )

    def test_no_n_samples_lt_20_conditional_banner(self):
        """The trigger condition itself (n_samples < 20 gating an
        st.info) must not reappear near the overview tiles. We accept
        base.n_samples < 20 in other contexts (e.g. warnings, feature
        gating, low-sample-count paths in legitimate code), so we
        specifically look for the banner block's identifying comment
        plus its conditional on the same few lines."""
        source = self._dashboard_source()
        # The removed banner was preceded by a header comment; that
        # comment is now replaced with "(removed April 2026)". If
        # someone re-adds the banner they'd likely restore a header
        # matching the old one. Assert the old header is gone.
        assert "# ---- LOOCV auto-override banner ----" not in source, (
            "The old LOOCV auto-override banner header comment is "
            "present in dashboard source. The block below it was "
            "removed in April 2026 and should remain removed."
        )


# =============================================================================
# Group 8.5 — dashboard bootstrap CI label matches summary() convention
# =============================================================================


class TestDashboardBootstrapLabel:
    """Pin dashboard/summary label symmetry on the bootstrap CI.

    Context: M2 added a conditional label to
    EnhancedBiomarkerResult.summary() so the bootstrap CI line reads
    'OOF AUC bootstrap CI' under M1's OOF path and
    'Training AUC bootstrap CI' otherwise. That closed the labeling
    ambiguity in the text output. But the Streamlit dashboard's
    Clinical Metrics tab had its own hardcoded label 'AUC bootstrap
    95% CI' that did not mirror the conditional. Spotted during the
    10a paper-figure review (screenshot shows 'AUC bootstrap 95% CI'
    while summary() on the same result says 'OOF AUC bootstrap CI').

    The fix: the dashboard block now picks its label from
    clin["oof_used"] using the same convention. This test guards both
    directions: the right labels are present, the wrong (pre-fix)
    label is not.

    Like TestStaleLOOCVBannerRemoved, this is a pure text scan of the
    dashboard source. No Streamlit runtime required.
    """

    def _dashboard_source(self) -> str:
        here = Path(__file__).resolve().parent.parent
        dashboard_path = (
            here / "raptor" / "dashboard" / "pages"
            / "08_\U0001f9ec_Biomarker_Discovery.py"
        )
        if not dashboard_path.exists():
            pages_dir = here / "raptor" / "dashboard" / "pages"
            candidates = list(pages_dir.glob("08_*Biomarker_Discovery*.py"))
            assert candidates
            dashboard_path = candidates[0]
        return dashboard_path.read_text(encoding="utf-8")

    def test_oof_label_string_present(self):
        """When oof_used is True, the dashboard must render 'OOF AUC
        bootstrap 95% CI'. The exact string must be in the source."""
        source = self._dashboard_source()
        assert "**OOF AUC bootstrap 95% CI**" in source, (
            "Dashboard source missing the M1-path bootstrap label "
            "'OOF AUC bootstrap 95% CI'. This label must match "
            "EnhancedBiomarkerResult.summary() convention for "
            "consistency between the text output and the dashboard."
        )

    def test_training_label_string_present(self):
        """When oof_used is False, the dashboard must render 'Training
        AUC bootstrap 95% CI'. The exact string must be in the source."""
        source = self._dashboard_source()
        assert "**Training AUC bootstrap 95% CI**" in source, (
            "Dashboard source missing the fallback-path bootstrap "
            "label 'Training AUC bootstrap 95% CI'. This label must "
            "be present to cover the pre-M1 pickled-result path "
            "where oof_used is False."
        )

    def test_old_ambiguous_label_is_absent(self):
        """The pre-M2 unconditional label must not reappear.

        The markdown-bold form with no OOF/Training prefix is the
        pre-fix label. If someone reverts the conditional back to an
        unconditional st.markdown, this test catches it.
        """
        source = self._dashboard_source()
        # The exact pre-fix form, as it appeared in the dashboard.
        pre_fix = '"**AUC bootstrap 95% CI**"'
        assert pre_fix not in source, (
            "Dashboard source contains the pre-M2 unconditional label "
            "'**AUC bootstrap 95% CI**'. This label is ambiguous about "
            "whether the bootstrap ran on OOF or training predictions. "
            "Use the conditional on clin['oof_used'] that picks "
            "between 'OOF AUC bootstrap 95% CI' and "
            "'Training AUC bootstrap 95% CI' instead."
        )

    def test_labels_match_summary_convention(self):
        """Cross-check: enhanced.py summary() must use the same
        {OOF AUC, Training AUC} vocabulary that the dashboard uses.

        enhanced.py builds its label with an f-string
        (f"{_auc_label} bootstrap CI: ...") so the literal
        "OOF AUC bootstrap" / "Training AUC bootstrap" strings don't
        appear in the source directly. We instead confirm the two
        vocabulary elements -- "OOF AUC" and "Training AUC" -- are
        both present in enhanced.py's summary() region AND the
        dashboard's Clinical Metrics tab. If either side drops one
        of the vocab words, this test fails.
        """
        here = Path(__file__).resolve().parent.parent
        enhanced_path = (
            here / "raptor" / "biomarker_discovery" / "enhanced.py"
        )
        enhanced_src = enhanced_path.read_text(encoding="utf-8")

        # Both labels are string literals in the f-string conditional
        # _auc_label = "OOF AUC" if ... else "Training AUC"
        # and both must remain present.
        assert '"OOF AUC"' in enhanced_src or "'OOF AUC'" in enhanced_src, (
            "enhanced.py does not contain the literal 'OOF AUC' label "
            "token. The summary() conditional must expose this token."
        )
        assert (
            '"Training AUC"' in enhanced_src
            or "'Training AUC'" in enhanced_src
        ), (
            "enhanced.py does not contain the literal 'Training AUC' "
            "label token. The summary() conditional must expose this "
            "token."
        )

        # Dashboard uses the same vocabulary but with "bootstrap 95% CI"
        # suffix instead of "bootstrap CI (95%)".
        dashboard_src = self._dashboard_source()
        assert "OOF AUC bootstrap 95% CI" in dashboard_src, (
            "Dashboard source missing 'OOF AUC bootstrap 95% CI'"
        )
        assert "Training AUC bootstrap 95% CI" in dashboard_src, (
            "Dashboard source missing 'Training AUC bootstrap 95% CI'"
        )