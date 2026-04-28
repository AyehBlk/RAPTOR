"""
Tests for dashboard polish batch helpers (M6 + M4 diagnostic rendering).

Streamlit rendering is hard to unit-test directly — the st.* calls
depend on a running server. But the branching logic that picks colors,
interpretation text, banner levels, and chance expectations is pure
Python and CAN be tested without Streamlit runtime. These tests cover
the four extracted helpers at their decision boundaries so the logic
can be refactored with confidence.

The helpers under test live at the top of
raptor/dashboard/pages/08___Biomarker_Discovery.py. Because that file
is a Streamlit page (not a module in the normal sense), we import it
by loading its source and executing just the helper definitions up to
the st.set_page_config call — avoiding any actual Streamlit runtime.
"""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
from types import ModuleType

import pytest


# =============================================================================
# Load dashboard helpers without triggering Streamlit runtime
# =============================================================================


def _load_dashboard_helpers() -> ModuleType:
    """Load only the helper functions and constants from the dashboard page.

    The dashboard file calls st.set_page_config() at module level, which
    would attempt to initialize Streamlit and hit a 'no scriptrun'
    error in a pytest context. We avoid that by loading the source,
    truncating at the '# PAGE CONFIG' boundary (which is how the helpers
    are delimited from the page runtime), and exec'ing the truncated
    source into a fresh module.
    """
    # Dashboard pages in this repo follow a "NN_<emoji>_<Name>.py" naming
    # convention (see all siblings in raptor/dashboard/pages/). The
    # emoji is literally in the filename.
    here = Path(__file__).resolve().parent.parent
    dashboard_path = (
        here / "raptor" / "dashboard" / "pages"
        / "08_\U0001f9ec_Biomarker_Discovery.py"
    )
    if not dashboard_path.exists():
        # Fallback search: any 08_*Biomarker_Discovery*.py in pages/
        pages_dir = here / "raptor" / "dashboard" / "pages"
        candidates = list(pages_dir.glob("08_*Biomarker_Discovery*.py"))
        if candidates:
            dashboard_path = candidates[0]
        else:
            raise FileNotFoundError(
                f"Could not find 08_*_Biomarker_Discovery.py under "
                f"{pages_dir}. Checked: {dashboard_path}"
            )
    source = dashboard_path.read_text(encoding="utf-8")

    # Truncate at '# PAGE CONFIG' - everything up to that line is
    # imports + helpers + constants. Past that line is Streamlit page
    # execution code (st.set_page_config and below).
    marker = "# PAGE CONFIG"
    idx = source.find(marker)
    if idx < 0:
        raise RuntimeError(
            f"Could not find '{marker}' section divider in dashboard source. "
            f"The dashboard page layout may have changed."
        )
    helper_source = source[:idx]

    # Exec into a synthetic module; suppress the streamlit/plotly imports
    # by substituting mock objects first.
    mod = ModuleType("dashboard_helpers_under_test")

    # Provide minimal mocks for third-party imports that the header block
    # pulls in; the helpers themselves only use numpy + pandas + pure Python.
    class _NoOp:
        def __getattr__(self, _):
            return _NoOp()
        def __call__(self, *args, **kwargs):
            return _NoOp()

    mod.__dict__["st"] = _NoOp()
    mod.__dict__["px"] = _NoOp()
    mod.__dict__["go"] = _NoOp()
    # numpy and pandas are real
    import numpy as np  # noqa: F401
    import pandas as pd  # noqa: F401
    mod.__dict__["np"] = np
    mod.__dict__["pd"] = pd
    # Path and typing used by helpers
    from pathlib import Path as _Path
    from typing import Dict, List, Optional, Tuple
    mod.__dict__.update({
        "Path": _Path, "Dict": Dict, "List": List,
        "Optional": Optional, "Tuple": Tuple,
        "io": __import__("io"),
        "sys": __import__("sys"),
    })

    # Strip the top-of-file imports - they reference raptor modules that
    # may not be importable cleanly in the test harness
    # Execute line by line, skipping import statements and the try-block
    # imports that the dashboard uses for graceful fallback.
    # Simplest approach: just exec what's after the last import block,
    # guarded by try/except for module loading.
    # We'll exec only the lines that define constants or functions we
    # care about, found by searching for the helper section marker.
    helper_marker = "# Dashboard polish helpers"
    helper_idx = source.find(helper_marker)
    if helper_idx < 0:
        raise RuntimeError(
            f"Could not find '{helper_marker}' marker. Dashboard layout may "
            f"have changed."
        )
    # The helper section starts at helper_idx and extends up to `# PAGE CONFIG`.
    helpers_only = source[helper_idx:idx]

    # Verify it parses
    ast.parse(helpers_only)

    exec(compile(helpers_only, str(dashboard_path), "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def dash():
    """Loaded dashboard helpers module."""
    return _load_dashboard_helpers()


# =============================================================================
# _stability_color — traffic-light boundaries
# =============================================================================


class TestStabilityColor:
    """Verify the traffic-light mapping at the boundaries."""

    def test_phi_at_excellent_boundary(self, dash):
        # 0.75 is inclusive lower bound for green.
        assert dash._stability_color(0.75) == "green"
        assert dash._stability_color(0.76) == "green"

    def test_phi_just_below_excellent(self, dash):
        assert dash._stability_color(0.74) == "amber"
        assert dash._stability_color(0.7499) == "amber"

    def test_phi_at_intermediate_boundary(self, dash):
        # 0.50 is inclusive lower bound for amber.
        assert dash._stability_color(0.50) == "amber"
        assert dash._stability_color(0.51) == "amber"

    def test_phi_just_below_intermediate(self, dash):
        assert dash._stability_color(0.49) == "red"
        assert dash._stability_color(0.4999) == "red"

    def test_phi_at_extremes(self, dash):
        assert dash._stability_color(1.00) == "green"
        assert dash._stability_color(0.00) == "red"

    def test_paper_scenarios(self, dash):
        """S2 (0.836), 10a (0.526), 10c (0.260) — the three figures."""
        assert dash._stability_color(0.836) == "green"
        assert dash._stability_color(0.526) == "amber"
        assert dash._stability_color(0.260) == "red"


# =============================================================================
# _stability_interpretation — text per category
# =============================================================================


class TestStabilityInterpretation:
    """Each category gets its specified one-sentence interpretation."""

    def test_excellent_mentions_robust(self, dash):
        text = dash._stability_interpretation("excellent")
        assert "agree strongly" in text.lower()
        assert "robust" in text.lower()

    def test_intermediate_mentions_partial_agreement(self, dash):
        text = dash._stability_interpretation("intermediate")
        assert "partial agreement" in text.lower()
        assert "validation" in text.lower()

    def test_poor_mentions_weak_signal(self, dash):
        text = dash._stability_interpretation("poor")
        assert "disagree" in text.lower()
        assert "weak" in text.lower() or "absent" in text.lower()

    def test_unknown_label_falls_back_to_poor(self, dash):
        # Defensive: unexpected labels should land on the 'poor' text
        # rather than raising, so the dashboard never dies on a string
        # it doesn't recognize.
        text = dash._stability_interpretation("garbage")
        assert "disagree" in text.lower()


# =============================================================================
# _optimism_banner_level — six-branch decision tree
# =============================================================================


class TestOptimismBannerLevel:
    """The decision tree. Six branches tested at representative points."""

    def test_strong_signal_tight_gap_is_green(self, dash):
        """CV AUC >= 0.75 AND |gap| < 0.05 -> green."""
        level, text = dash._optimism_banner_level(cv_auc=0.90, gap=0.02)
        assert level == "green"
        assert "genuine signal" in text.lower()

    def test_chance_level_tight_gap_is_amber(self, dash):
        """CV AUC < 0.65 AND |gap| < 0.05 -> amber (was the reporting bug).

        This is 10c's scenario: CV AUC 0.583, gap -0.013. Pre-fix, this
        landed on green. Post-fix, it must land on amber with a
        chance-level message.
        """
        level, text = dash._optimism_banner_level(cv_auc=0.583, gap=-0.013)
        assert level == "amber"
        assert "chance level" in text.lower()

    def test_borderline_auc_tight_gap_is_amber(self, dash):
        """0.65 <= CV AUC < 0.75 AND |gap| < 0.05 -> amber."""
        level, text = dash._optimism_banner_level(cv_auc=0.70, gap=0.02)
        assert level == "amber"
        assert "borderline" in text.lower()

    def test_moderate_gap_is_amber(self, dash):
        """|gap| in [0.05, 0.10) -> amber regardless of CV AUC."""
        level, text = dash._optimism_banner_level(cv_auc=0.90, gap=0.07)
        assert level == "amber"
        assert "moderate" in text.lower()
        # Also on chance-level AUC
        level2, text2 = dash._optimism_banner_level(cv_auc=0.55, gap=0.07)
        assert level2 == "amber"

    def test_large_positive_gap_is_red(self, dash):
        """gap >= 0.10 -> red (overfitting)."""
        level, text = dash._optimism_banner_level(cv_auc=0.90, gap=0.15)
        assert level == "red"
        assert "large" in text.lower() or "substantially exceeds" in text.lower()

    def test_large_negative_gap_is_amber(self, dash):
        """gap <= -0.10 -> amber (atypical negative gap)."""
        level, text = dash._optimism_banner_level(cv_auc=0.70, gap=-0.15)
        assert level == "amber"
        assert "atypical" in text.lower() or "negative" in text.lower()

    def test_boundary_auc_0p75_tight_gap(self, dash):
        """Exactly 0.75 -> green (inclusive lower bound)."""
        level, _ = dash._optimism_banner_level(cv_auc=0.75, gap=0.02)
        assert level == "green"

    def test_boundary_auc_0p65_tight_gap(self, dash):
        """Exactly 0.65 -> amber borderline (inclusive lower bound, not chance)."""
        level, text = dash._optimism_banner_level(cv_auc=0.65, gap=0.02)
        assert level == "amber"
        assert "borderline" in text.lower()

    def test_boundary_gap_0p05_is_moderate(self, dash):
        """Exactly |gap|=0.05 -> moderate amber."""
        level, text = dash._optimism_banner_level(cv_auc=0.90, gap=0.05)
        assert level == "amber"
        assert "moderate" in text.lower()

    def test_boundary_gap_0p10_is_red(self, dash):
        """Exactly gap=0.10 -> red."""
        level, _ = dash._optimism_banner_level(cv_auc=0.90, gap=0.10)
        assert level == "red"


# =============================================================================
# _m4_chance_expectation — alpha * n_genes, rounded
# =============================================================================


class TestM4ChanceExpectation:
    """Trivial arithmetic helper; verify the rounding."""

    def test_typical_values(self, dash):
        assert dash._m4_chance_expectation(0.05, 200) == 10
        assert dash._m4_chance_expectation(0.05, 20000) == 1000
        assert dash._m4_chance_expectation(0.01, 200) == 2

    def test_rounds_to_nearest(self, dash):
        # 0.05 * 17 = 0.85 -> rounds to 1
        assert dash._m4_chance_expectation(0.05, 17) == 1
        # 0.05 * 30 = 1.5 -> Python banker's round to 2
        assert dash._m4_chance_expectation(0.05, 30) == 2

    def test_zero_genes(self, dash):
        assert dash._m4_chance_expectation(0.05, 0) == 0


# =============================================================================
# Module-level constants sanity
# =============================================================================


class TestConstants:
    """The thresholds must match what the scoping doc specified."""

    def test_strong_signal_threshold(self, dash):
        assert dash.OPTIMISM_AUC_STRONG_SIGNAL == 0.75

    def test_chance_level_threshold(self, dash):
        assert dash.OPTIMISM_AUC_CHANCE_LEVEL == 0.65

    def test_excellent_threshold_matches_nogueira(self, dash):
        assert dash.STABILITY_EXCELLENT_THRESHOLD == 0.75

    def test_intermediate_threshold(self, dash):
        assert dash.STABILITY_INTERMEDIATE_THRESHOLD == 0.50


# =============================================================================
# Small-n advisory (Enhancement A): tiered banner by cohort size
# =============================================================================


class TestSmallNAdvisory:
    """Boundary tests for the cohort-size-dependent advisory.

    Tiers (from the dashboard helper):
        n < 20            -> ('strong', ...)   strong advisory
        20 <= n < 50      -> ('moderate', ...) moderate advisory
        n >= 50           -> None              no advisory

    These boundaries are calibrated against the seed-stability
    investigation findings: at n=10 CV AUC range was 0.42 across 5
    seeds, and at n=60 CV AUC range was 0.000. The tiers must match
    those findings.
    """

    def test_n_below_strong_threshold(self, dash):
        """n=19: strong advisory."""
        result = dash._small_n_advisory(19)
        assert result is not None
        level, message = result
        assert level == "strong"
        assert "n=19" in message
        # The strong advisory must specifically tell users that
        # CV AUC point estimates are unreliable.
        assert "unreliable" in message.lower() or "variance" in message.lower()
        # And specifically tell users to trust panel composition + Phi.
        assert "panel composition" in message.lower()

    def test_n_at_strong_boundary(self, dash):
        """n=20: moderate advisory (20 is inclusive lower bound for moderate)."""
        result = dash._small_n_advisory(20)
        assert result is not None
        level, _ = result
        assert level == "moderate"

    def test_n_below_moderate_boundary(self, dash):
        """n=49: moderate advisory."""
        result = dash._small_n_advisory(49)
        assert result is not None
        level, _ = result
        assert level == "moderate"

    def test_n_at_moderate_boundary(self, dash):
        """n=50: no advisory (50 is inclusive lower bound for 'safe')."""
        result = dash._small_n_advisory(50)
        assert result is None

    def test_n_well_above_moderate(self, dash):
        """n=200: no advisory."""
        assert dash._small_n_advisory(200) is None

    def test_n_paper_scenarios(self, dash):
        """Paper-scenario calibration:
        S2  (n=60): no advisory (above the 50 cutoff)
        10a (n=10): strong advisory
        10c (n=60): no advisory
        """
        # S2 / 10c
        assert dash._small_n_advisory(60) is None
        # 10a
        result = dash._small_n_advisory(10)
        assert result is not None
        assert result[0] == "strong"

    def test_n_zero_or_negative_handled(self, dash):
        """Defensive: n=0 must not crash. Treat as strong advisory."""
        result = dash._small_n_advisory(0)
        assert result is not None
        assert result[0] == "strong"

    def test_constants_exposed(self, dash):
        """The threshold constants must be module-level for tunability."""
        assert dash.SMALL_N_THRESHOLD_STRONG == 20
        assert dash.SMALL_N_THRESHOLD_MODERATE == 50

    def test_default_verification_seeds_distinct(self, dash):
        """The DEFAULT_VERIFICATION_SEEDS must be three distinct values."""
        seeds = dash.DEFAULT_VERIFICATION_SEEDS
        assert len(seeds) == 3
        assert len(set(seeds)) == 3


# =============================================================================
# Multi-seed verification (Enhancement B): _summarize_seed_runs
# =============================================================================


class TestSummarizeSeedRuns:
    """Verify the invariant/variant split that the dashboard surfaces.

    These tests check the helper that powers the "trust the
    invariants, doubt the variants" interpretation summary at the
    bottom of the multi-seed comparison table. The helper must
    correctly identify what is stable across seeds vs what is not.
    """

    def _make_row(
        self, seed: int, cv_auc: float, gap: float, banner: str,
        phi: float, phi_label: str, panel: str,
    ) -> dict:
        """Construct a row matching run_one_seed's return shape."""
        return {
            "seed": seed, "cv_auc": cv_auc, "gap": gap,
            "banner": banner, "phi": phi, "phi_label": phi_label,
            "panel": panel,
            "train_auc": 1.0, "boot_lo": 0.5, "boot_hi": 1.0,
            "best_clf": "random_forest",
        }

    def test_empty_input_returns_safe_defaults(self, dash):
        s = dash._summarize_seed_runs([])
        assert s["n_runs"] == 0
        assert s["invariant_genes"] == []
        assert s["variant_genes"] == []
        assert s["invariant_phi_label"] is None
        assert s["invariant_banner"] is None

    def test_three_runs_with_invariant_phi_and_genes(self, dash):
        """10a-style: Phi 'poor' across all runs, 3 cancer genes
        invariant, companion genes variable, banner varies."""
        rows = [
            self._make_row(0, 0.90, +0.10, "amber (moderate)",
                           0.32, "poor", "TP53,BRCA1,MYC,GENE_032,GENE_046"),
            self._make_row(1, 1.00, +0.00, "green (genuine signal)",
                           0.34, "poor", "TP53,BRCA1,MYC,GENE_046,GENE_078"),
            self._make_row(456, 0.58, +0.42, "red (overfit)",
                           0.12, "poor", "TP53,BRCA1,MYC,GENE_032,GENE_091"),
        ]
        s = dash._summarize_seed_runs(rows)

        assert s["n_runs"] == 3
        assert s["invariant_phi_label"] == "poor"
        assert s["invariant_banner"] is None  # banner varies

        # The 3 cancer genes are in every panel
        assert set(s["invariant_genes"]) == {"TP53", "BRCA1", "MYC"}
        # Variants: at least the union of companion genes minus invariants
        assert "GENE_032" in s["variant_genes"]
        assert "GENE_046" in s["variant_genes"]
        # CV AUC range captures the spread
        assert s["cv_auc_range"] == (0.58, 1.00)
        # Gap range captures the spread
        assert s["gap_range"] == (0.0, 0.42)
        # Phi range
        assert s["phi_range"][0] == 0.12
        assert s["phi_range"][1] == 0.34

    def test_s2_style_full_invariance(self, dash):
        """S2-style: everything invariant across seeds."""
        rows = [
            self._make_row(0, 1.00, 0.00, "green (genuine signal)",
                           1.00, "excellent", "TP53,BRCA1,MYC,RB1,STAT3"),
            self._make_row(42, 1.00, 0.00, "green (genuine signal)",
                           1.00, "excellent", "TP53,BRCA1,MYC,RB1,STAT3"),
            self._make_row(123, 1.00, 0.00, "green (genuine signal)",
                           0.92, "excellent", "TP53,BRCA1,MYC,RB1,STAT3"),
        ]
        s = dash._summarize_seed_runs(rows)
        assert s["invariant_phi_label"] == "excellent"
        assert s["invariant_banner"] == "green (genuine signal)"
        assert set(s["invariant_genes"]) == {
            "TP53", "BRCA1", "MYC", "RB1", "STAT3",
        }
        assert s["variant_genes"] == []
        assert s["cv_auc_range"] == (1.0, 1.0)

    def test_no_overlap_between_panels_yields_no_invariants(self, dash):
        """If every run's panel is disjoint, no invariant genes."""
        rows = [
            self._make_row(0, 0.6, 0.1, "amber", 0.3, "poor", "A,B,C"),
            self._make_row(1, 0.6, 0.1, "amber", 0.3, "poor", "D,E,F"),
        ]
        s = dash._summarize_seed_runs(rows)
        assert s["invariant_genes"] == []
        assert sorted(s["variant_genes"]) == ["A", "B", "C", "D", "E", "F"]

    def test_phi_label_varies_returns_none(self, dash):
        """When phi_label changes, invariant_phi_label is None."""
        rows = [
            self._make_row(0, 0.8, 0.0, "green", 0.3, "poor", "TP53"),
            self._make_row(1, 0.8, 0.0, "green", 0.6, "intermediate", "TP53"),
        ]
        s = dash._summarize_seed_runs(rows)
        assert s["invariant_phi_label"] is None

    def test_single_run_treats_everything_as_invariant(self, dash):
        """One run: everything in it is trivially invariant."""
        rows = [
            self._make_row(42, 1.0, 0.0, "green", 1.0,
                           "excellent", "TP53,BRCA1,MYC"),
        ]
        s = dash._summarize_seed_runs(rows)
        assert s["n_runs"] == 1
        assert s["invariant_phi_label"] == "excellent"
        assert s["invariant_banner"] == "green"
        assert set(s["invariant_genes"]) == {"TP53", "BRCA1", "MYC"}
        assert s["variant_genes"] == []
        assert s["cv_auc_range"] == (1.0, 1.0)

    def test_10c_style_same_banner_different_phi(self, dash):
        """10c-style: all 'red (overfit)' but Phi values vary."""
        rows = [
            self._make_row(0, 0.63, 0.26, "red (overfit)",
                           0.36, "poor", "G_137,G_075,G_074,G_090,G_139"),
            self._make_row(42, 0.63, 0.37, "red (overfit)",
                           0.36, "poor", "G_137,G_075,G_074,G_090,G_038"),
            self._make_row(456, 0.44, 0.47, "red (overfit)",
                           0.20, "poor", "G_137,G_075,G_074,G_090,G_139"),
        ]
        s = dash._summarize_seed_runs(rows)
        assert s["invariant_banner"] == "red (overfit)"
        assert s["invariant_phi_label"] == "poor"
        # 4 noise genes appear in all 3 runs (this is what 10c showed
        # in the actual investigation -- panel composition is mostly
        # stable on noise too at adequate n, just nothing biological)
        assert "G_137" in s["invariant_genes"]
        assert "G_075" in s["invariant_genes"]


# =============================================================================
# Smoke tests: verify the multi-seed render helper exists and is wired
# =============================================================================


class TestMultiSeedHelpersWired:
    """Confirm the multi-seed helpers are present and importable.

    Full integration testing of the button + progressive rendering
    requires Streamlit's runtime which we don't have in unit tests.
    These tests confirm the helpers exist with the right signatures
    and that the button code path is reachable.
    """

    def test_render_helper_exists(self, dash):
        # The progressive-render helper is defined in the dashboard
        # source but lives outside the polish-helpers section, so it
        # is not in the loaded `dash` module from
        # _load_dashboard_helpers (which truncates at PAGE CONFIG).
        # Instead we scan the source for its definition.
        from pathlib import Path
        here = Path(__file__).resolve().parent.parent
        dashboard_path = (
            here / "raptor" / "dashboard" / "pages"
            / "08_\U0001f9ec_Biomarker_Discovery.py"
        )
        if not dashboard_path.exists():
            pages_dir = here / "raptor" / "dashboard" / "pages"
            candidates = list(
                pages_dir.glob("08_*Biomarker_Discovery*.py")
            )
            assert candidates
            dashboard_path = candidates[0]
        source = dashboard_path.read_text(encoding="utf-8")
        assert "def _render_multi_seed_verification" in source
        assert "def _draw_multiseed_table" in source
        assert "def _draw_invariants_summary" in source

    def test_button_key_present(self, dash):
        """The verify button has a stable key for testability."""
        from pathlib import Path
        here = Path(__file__).resolve().parent.parent
        dashboard_path = (
            here / "raptor" / "dashboard" / "pages"
            / "08_\U0001f9ec_Biomarker_Discovery.py"
        )
        if not dashboard_path.exists():
            pages_dir = here / "raptor" / "dashboard" / "pages"
            candidates = list(
                pages_dir.glob("08_*Biomarker_Discovery*.py")
            )
            assert candidates
            dashboard_path = candidates[0]
        source = dashboard_path.read_text(encoding="utf-8")
        assert 'key="btn_multi_seed_verify"' in source

    def test_session_state_save_for_multiseed(self, dash):
        """The discovery success path saves the inputs needed by the
        multi-seed verification button into st.session_state.

        The save block must avoid widget-owned keys (m10_min_panel,
        m10_max_panel, m10_seed) because Streamlit raises
        'cannot be modified after the widget ... is instantiated' on
        those. Non-widget keys use the _resolved suffix to make the
        non-collision intent explicit.
        """
        from pathlib import Path
        here = Path(__file__).resolve().parent.parent
        dashboard_path = (
            here / "raptor" / "dashboard" / "pages"
            / "08_\U0001f9ec_Biomarker_Discovery.py"
        )
        if not dashboard_path.exists():
            pages_dir = here / "raptor" / "dashboard" / "pages"
            candidates = list(
                pages_dir.glob("08_*Biomarker_Discovery*.py")
            )
            assert candidates
            dashboard_path = candidates[0]
        source = dashboard_path.read_text(encoding="utf-8")
        # Required saves in the success path
        assert "st.session_state['biomarker_counts']" in source
        assert "st.session_state['biomarker_metadata']" in source
        assert "m10_group_column_resolved" in source
        assert "m10_reference_group_resolved" in source
        assert "m10_baseline_group_resolved" in source
        assert "multi_seed_results" in source  # cache key reference

        # Forbidden: widget-owned-key writes in the success path
        # (these triggered the 'cannot be modified after the widget'
        # crash on the first deployment).
        forbidden_writes = [
            "st.session_state['m10_min_panel'] =",
            "st.session_state['m10_max_panel'] =",
            "st.session_state['m10_seed'] =",
        ]
        for forbidden in forbidden_writes:
            assert forbidden not in source, (
                f"Forbidden write to widget-owned key: {forbidden!r}. "
                f"This causes a Streamlit 'cannot be modified after "
                f"the widget ... is instantiated' error."
            )

    def test_render_call_does_not_pass_clin_kwarg(self):
        """Regression: the call site for _render_multi_seed_verification
        is on the Panel tab, where the local `clin` variable is NOT in
        scope (clin is defined inside the Clinical Metrics tab block).
        The first deployment of B passed `clin=clin` and crashed at
        button-click time with NameError. The fix derives clin from
        `result` inside the function. This test guards against the bug
        coming back."""
        from pathlib import Path
        here = Path(__file__).resolve().parent.parent
        dashboard_path = (
            here / "raptor" / "dashboard" / "pages"
            / "08_\U0001f9ec_Biomarker_Discovery.py"
        )
        if not dashboard_path.exists():
            pages_dir = here / "raptor" / "dashboard" / "pages"
            candidates = list(
                pages_dir.glob("08_*Biomarker_Discovery*.py")
            )
            assert candidates
            dashboard_path = candidates[0]
        source = dashboard_path.read_text(encoding="utf-8")

        # The call site must NOT pass clin as a kwarg.
        assert "clin=clin" not in source, (
            "The Panel-tab call to _render_multi_seed_verification "
            "must not pass `clin=clin`. `clin` is defined in the "
            "Clinical Metrics tab scope and is not visible on the "
            "Panel tab; passing it triggers NameError at button-"
            "click time. Derive clin from `result` inside the "
            "function instead."
        )

        # The function definition must derive clin from result.
        assert (
            'getattr(result, "clinical_metrics", None) or {}' in source
            or "getattr(result, 'clinical_metrics', None) or {}" in source
        ), (
            "_render_multi_seed_verification must derive `clin` from "
            "`result` (e.g. via getattr(result, 'clinical_metrics', "
            "None) or {{}}) so it can be invoked from any tab without "
            "depending on a tab-local `clin` variable."
        )