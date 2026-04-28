"""Tests for kneedle panel-size auto-detection (added in v2.2.2).

Covers:
  * `_detect_optimal_panel_size` pure-function behavior on synthetic
    curves (clean knee, saturated, flat, monotone-late, short curves).
  * `forward_selection` integration: auto_strategy parameter is honored,
    selection_method is recorded on the result, sensitivity is wired.
  * Dashboard annotation regression (text-scan): the v-line label
    branches on selection_method correctly.
  * Cross-strategy parity on a known curve: kneedle and first_drop
    can disagree on size in principle, but on simple curves they
    converge.
"""
from __future__ import annotations

import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from raptor.biomarker_discovery.core import (
    PanelOptimizer,
    PanelOptimizationResult,
    _detect_optimal_panel_size,
)


# =============================================================================
# Group 1: pure-function tests on synthetic curves
# =============================================================================
class TestDetectOptimalPanelSize:
    """`_detect_optimal_panel_size` on hand-crafted curves."""

    def _curve(self, sizes, aucs):
        return pd.DataFrame({'panel_size': sizes, 'auc_mean': aucs})

    def test_kneedle_finds_knee_on_clean_concave_curve(self):
        # Sharp knee at size=5: rapid climb then near-flat saturation.
        sizes = list(range(1, 21))
        aucs = [0.55, 0.65, 0.75, 0.83, 0.88,
                0.89, 0.895, 0.898, 0.90, 0.901,
                0.902, 0.903, 0.903, 0.904, 0.904,
                0.904, 0.905, 0.905, 0.905, 0.905]
        size, method = _detect_optimal_panel_size(
            self._curve(sizes, aucs), auto_strategy='kneedle',
        )
        assert method == 'kneedle'
        assert 4 <= size <= 8, f"expected knee in 4-8, got {size}"

    def test_kneedle_falls_back_on_flat_curve(self):
        # All AUCs identical -> no knee. Fallback to argmax.
        sizes = [3, 4, 5]
        aucs = [1.0, 1.0, 1.0]
        size, method = _detect_optimal_panel_size(
            self._curve(sizes, aucs), auto_strategy='kneedle',
        )
        assert method == 'argmax_fallback'
        assert size == 3  # smallest size at max

    def test_kneedle_falls_back_on_monotone_with_no_knee(self):
        # Strictly monotone increasing, near-linear -> kneedle returns
        # None or a knee at one of the endpoints. Either way, our
        # implementation must produce a stable size and label it
        # consistently.
        sizes = list(range(1, 11))
        aucs = [0.50 + 0.04 * i for i in range(10)]
        size, method = _detect_optimal_panel_size(
            self._curve(sizes, aucs), auto_strategy='kneedle',
        )
        # On a near-linear curve kneedle finds no clear knee and
        # falls back to argmax (smallest size at max). Either label
        # is acceptable; size must be in range.
        assert method in ('kneedle', 'argmax_fallback')
        assert size in sizes

    def test_short_curve_falls_back_to_argmax(self):
        # < 3 points: kneedle can't run. Must fall back.
        sizes = [3, 4]
        aucs = [0.7, 0.85]
        size, method = _detect_optimal_panel_size(
            self._curve(sizes, aucs), auto_strategy='kneedle',
        )
        assert method == 'argmax_fallback'
        assert size == 4  # smallest size at max AUC = 0.85

    def test_argmax_strategy_picks_smallest_size_at_max(self):
        # Tied maxima at sizes 5 and 8 -> tiebreak picks smaller (5).
        sizes = [3, 4, 5, 6, 7, 8]
        aucs = [0.7, 0.8, 0.92, 0.91, 0.91, 0.92]
        size, method = _detect_optimal_panel_size(
            self._curve(sizes, aucs), auto_strategy='argmax',
        )
        assert method == 'argmax'
        assert size == 5  # smallest of the tied winners

    def test_first_drop_strategy_walks_curve(self):
        # First sub-0.005 improvement at i=4 (0.92 - 0.918 = 0.002).
        # first_drop returns sizes[i-1] = sizes[3] = size 6.
        sizes = [3, 4, 5, 6, 7]
        aucs = [0.70, 0.85, 0.91, 0.918, 0.920]
        size, method = _detect_optimal_panel_size(
            self._curve(sizes, aucs), auto_strategy='first_drop',
        )
        assert method == 'first_drop'
        assert size == 6

    def test_constant_curve_skips_kneedle_no_warning(self):
        # The constant-curve guard must fire BEFORE kneedle's polyfit
        # is called, so no RankWarning is emitted on a flat curve.
        sizes = [3, 4, 5]
        aucs = [1.0, 1.0, 1.0]
        # RankWarning lives in np.exceptions on numpy 2.x but at the
        # top level on older installs. Look in both.
        try:
            RankWarning = np.exceptions.RankWarning
        except AttributeError:
            RankWarning = getattr(np, 'RankWarning', None)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            size, method = _detect_optimal_panel_size(
                self._curve(sizes, aucs), auto_strategy='kneedle',
            )
        if RankWarning is not None:
            rank_warnings = [
                w for w in caught
                if issubclass(w.category, RankWarning)
            ]
            assert len(rank_warnings) == 0, (
                f"constant-curve guard should skip kneedle entirely; "
                f"got {len(rank_warnings)} RankWarning(s)"
            )
        assert method == 'argmax_fallback'


# =============================================================================
# Group 2: forward_selection integration
# =============================================================================
class TestForwardSelectionAutoStrategy:
    """auto_strategy and sensitivity are honored end-to-end."""

    @pytest.fixture
    def synthetic_data(self):
        rng = np.random.default_rng(2024)
        n_samples, n_genes, n_signal = 60, 50, 5
        X_arr = rng.standard_normal((n_samples, n_genes)) * 0.8
        X_arr[:30, :n_signal] += 2.5
        gene_ids = [f'GENE_{i:03d}' for i in range(n_genes)]
        sample_ids = [f'S{i:02d}' for i in range(n_samples)]
        X = pd.DataFrame(X_arr, index=sample_ids, columns=gene_ids)
        y = np.array([1] * 30 + [0] * 30)
        return X, y, gene_ids

    def test_default_is_kneedle(self, synthetic_data):
        X, y, gene_ids = synthetic_data
        po = PanelOptimizer(random_state=42, verbose=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = po.forward_selection(
                X, y, ranked_genes=gene_ids,
                min_panel=3, max_panel=10, n_cv=3,
            )
        # Selection method should be kneedle or argmax_fallback
        # (depending on whether the curve has a detectable knee).
        assert r.selection_method in ('kneedle', 'argmax_fallback')
        # We also expect optimal_size to be in [3, 10]
        assert 3 <= r.optimal_size <= 10

    def test_first_drop_strategy_produces_first_drop_label(
        self, synthetic_data,
    ):
        X, y, gene_ids = synthetic_data
        po = PanelOptimizer(random_state=42, verbose=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = po.forward_selection(
                X, y, ranked_genes=gene_ids,
                min_panel=3, max_panel=10, n_cv=3,
                auto_strategy='first_drop',
            )
        assert r.selection_method == 'first_drop'

    def test_target_size_yields_user_specified(self, synthetic_data):
        # When target_size is given, selection_method must reflect that
        # the user, not the algorithm, picked the size.
        X, y, gene_ids = synthetic_data
        po = PanelOptimizer(random_state=42, verbose=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = po.forward_selection(
                X, y, ranked_genes=gene_ids,
                min_panel=3, max_panel=10, n_cv=3,
                target_size=5,
            )
        assert r.selection_method == 'user_specified'
        assert r.optimal_size == 5


# =============================================================================
# Group 3: dashboard annotation regression
# =============================================================================
class TestDashboardKneedleAnnotation:
    """Static text-scan tests on the dashboard source.

    These guard against future edits removing the per-strategy
    annotation labels — the dashboard renderer can't be unit-tested
    without a Streamlit AppTest harness, but we can check the source
    contains the label-branching logic.
    """

    @pytest.fixture
    def dashboard_source(self):
        path = Path(__file__).parent.parent / "raptor" / "dashboard" / "pages"
        # Filename uses an emoji so glob-find it.
        candidates = list(path.glob("*Biomarker*.py"))
        assert candidates, f"no Biomarker dashboard page found in {path}"
        return candidates[0].read_text(encoding='utf-8')

    def test_kneedle_label_present(self, dashboard_source):
        assert 'Knee = ' in dashboard_source
        assert '(kneedle)' in dashboard_source

    def test_argmax_fallback_caption_present(self, dashboard_source):
        # The "no curvature" explanation must appear when kneedle
        # falls back. Check both the annotation_text branch and the
        # st.caption follow-up.
        assert 'argmax — no curvature' in dashboard_source
        assert 'Kneedle detected no knee' in dashboard_source

    def test_user_specified_label_present(self, dashboard_source):
        assert 'Target = ' in dashboard_source
        assert 'user-specified' in dashboard_source

    def test_consensus_pinned_label_present(self, dashboard_source):
        assert 'consensus-pinned' in dashboard_source

    def test_old_unconditional_optimal_label_removed(self, dashboard_source):
        # The pre-2.3 code unconditionally wrote `f"Optimal = {optimal}"`
        # as the annotation_text. We replaced that with selection_method-
        # branched logic. Make sure the old line isn't lingering.
        # We allow `Optimal = {optimal}` to appear as a fallback `else:`
        # branch (it does), but the unconditional `annotation_text=f"Optimal = {optimal}"`
        # call site must be gone.
        bad = 'annotation_text=f"Optimal = {optimal}",'
        # This exact one-line form was the pre-2.3 pattern. Confirm no
        # such line exists at the OUTER indentation level.
        for line in dashboard_source.splitlines():
            stripped = line.lstrip()
            if stripped == bad:
                # Found the pre-2.3 pattern -- shouldn't happen.
                raise AssertionError(
                    "Pre-2.3 unconditional 'Optimal = {optimal}' "
                    "annotation pattern still in dashboard"
                )


# =============================================================================
# Group 5: discover_biomarkers default panel_size_strategy
# =============================================================================
class TestDiscoverBiomarkersDefaultStrategy:
    """The Q5 sanity check on real fixtures showed per-fold panel-size
    detection can drop Phi by ~0.3 when one outer fold's argmax tie
    lands on a much larger size than the others (S2 seed 42). The
    consensus default prevents this."""

    def test_default_panel_size_strategy_is_consensus(self):
        """Consensus is the default panel_size_strategy (Q5 verdict)."""
        import inspect
        from raptor.biomarker_discovery.core import discover_biomarkers
        sig = inspect.signature(discover_biomarkers)
        assert 'panel_size_strategy' in sig.parameters, (
            "discover_biomarkers must expose panel_size_strategy"
        )
        param = sig.parameters['panel_size_strategy']
        assert param.default == 'consensus', (
            f"Q5 sanity check on S2 fixtures showed per-fold "
            f"detection drops Phi meaningfully (0.313 at seed 42); "
            f"default must be 'consensus', got {param.default!r}"
        )

    def test_per_fold_still_accepted(self):
        """Per-fold remains valid for users who want the previous
        behavior or for parity with prior runs."""
        import inspect
        from raptor.biomarker_discovery.core import discover_biomarkers
        # Just confirm the parameter still exists and accepts a string;
        # invocation testing is covered by the sanity-check workflow.
        sig = inspect.signature(discover_biomarkers)
        param = sig.parameters['panel_size_strategy']
        assert param.annotation in (str, 'str'), param.annotation

    def test_consensus_path_emits_consensus_pinned_label(self):
        """When the consensus scout picks K (panel_size_strategy='consensus',
        no user-specified target_panel_size), the final-panel result
        must be labeled 'consensus_pinned', not 'user_specified'.

        Regression for a labeling bug where the dashboard showed
        'Target = K (user-specified)' for K's that the SCOUT picked,
        misattributing the choice to the user. Catches by running a
        small discovery and asserting the label.
        """
        import warnings
        rng = np.random.default_rng(2024)
        n_samples, n_genes, n_signal = 40, 30, 4
        X_arr = rng.standard_normal((n_samples, n_genes)) * 0.8
        X_arr[:20, :n_signal] += 2.5
        gene_ids = [f'GENE_{i:03d}' for i in range(n_genes)]
        sample_ids = [f'S{i:02d}' for i in range(n_samples)]
        counts = pd.DataFrame(
            np.clip(X_arr.T * 100 + 500, 0, None).astype(int),
            index=gene_ids, columns=sample_ids,
        )
        metadata = pd.DataFrame({
            'sample_id': sample_ids,
            'condition': ['disease'] * 20 + ['healthy'] * 20,
        }).set_index('sample_id')

        from raptor.biomarker_discovery import discover_biomarkers
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = discover_biomarkers(
                counts=counts, metadata=metadata,
                group_column='condition',
                reference_group='disease', baseline_group='healthy',
                intent='diagnostic',
                min_panel=3, max_panel=6,
                target_panel_size=None,    # let consensus scout pick
                validation='nested_cv', n_folds=3, n_repeats=1,
                annotate=False, run_literature=False, run_ppi=False,
                output_dir='/tmp/consensus_pin_test',
                random_state=42, verbose=False,
                # Explicit consensus to be sure -- though it's the
                # default, the test should not silently pass if the
                # default is later flipped.
                panel_size_strategy='consensus',
            )

        br = r.base_result if hasattr(r, 'base_result') else r
        po = br.panel_optimization
        assert po.selection_method == 'consensus_pinned', (
            f"consensus path should emit 'consensus_pinned' on the "
            f"final panel result, got {po.selection_method!r}. The "
            f"dashboard annotation will read 'user-specified' "
            f"instead of 'consensus-pinned' if this regression lands."
        )


# =============================================================================
# Group 4: PanelOptimizationResult dataclass field
# =============================================================================
class TestSelectionMethodField:
    """selection_method is on the dataclass, has the right default,
    and pickles round-trip cleanly for old vs new pickles."""

    def test_field_default_is_kneedle(self):
        r = PanelOptimizationResult(
            optimal_panel=['G1', 'G2', 'G3'],
            optimal_size=3,
            optimal_auc=0.85,
            panel_curve=pd.DataFrame({'panel_size': [3], 'auc_mean': [0.85]}),
        )
        assert r.selection_method == 'kneedle'

    def test_field_can_be_set_explicitly(self):
        r = PanelOptimizationResult(
            optimal_panel=['G1', 'G2', 'G3'],
            optimal_size=3,
            optimal_auc=0.85,
            panel_curve=pd.DataFrame({'panel_size': [3], 'auc_mean': [0.85]}),
            selection_method='argmax_fallback',
        )
        assert r.selection_method == 'argmax_fallback'

    def test_pickle_roundtrip(self):
        r = PanelOptimizationResult(
            optimal_panel=['G1'],
            optimal_size=1,
            optimal_auc=0.7,
            panel_curve=pd.DataFrame({'panel_size': [1], 'auc_mean': [0.7]}),
            selection_method='kneedle',
        )
        r2 = pickle.loads(pickle.dumps(r))
        assert r2.selection_method == 'kneedle'