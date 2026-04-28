"""
RAPTOR Module 10: CLI tests for kneedle panel-size detection.

Covers the three new flags on the ``raptor biomarker`` command:
    --auto-panel-strategy [kneedle|argmax|first_drop]
    --panel-sensitivity FLOAT
    --panel-size-strategy [per_fold|consensus]

Verifies:
    * Click parsing accepts valid values and rejects invalid ones
    * Each strategy produces the expected ``selection_method`` in the
      pickled result
    * --help advertises all three flags
    * Default behavior matches the dashboard (kneedle / consensus)

Run with:
    pytest tests/test_cli_kneedle.py -v --tb=short
    pytest tests/test_cli_kneedle.py -v -k "selection_method"

Author: Ayeh Bolouki
"""
import pickle
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from raptor.cli import main


# ============================================================================
# SHARED FIXTURE — small, fast synthetic with strong signal for stable runs
# ============================================================================
@pytest.fixture(scope="module")
def synth_data_dir():
    """Synthetic counts + metadata with strong signal so all strategies
    converge to a panel containing the planted DE genes.

    Strong signal keeps tests deterministic and fast: ~30s per CLI
    invocation. The fixture is module-scoped because every test rebuilds
    the data, but the data itself is deterministic (np.random.seed)."""
    tmpdir = Path(tempfile.mkdtemp(prefix="raptor_cli_kneedle_"))

    rng = np.random.default_rng(2024)
    n_samples, n_genes, n_signal = 40, 60, 5
    X = rng.standard_normal((n_samples, n_genes)) * 0.8
    # Half samples get a strong signal in the first n_signal genes
    X[:n_samples // 2, :n_signal] += 2.5

    gene_ids = [f"GENE_{i:03d}" for i in range(n_genes)]
    sample_ids = [f"S{i:02d}" for i in range(n_samples)]
    counts = pd.DataFrame(
        np.clip(X.T * 100 + 500, 0, None).astype(int),
        index=gene_ids, columns=sample_ids,
    )
    metadata = pd.DataFrame({
        'sample_id': sample_ids,
        'condition': (
            ['disease'] * (n_samples // 2)
            + ['healthy'] * (n_samples - n_samples // 2)
        ),
    })

    counts.to_csv(tmpdir / "counts.csv")
    metadata.to_csv(tmpdir / "metadata.csv", index=False)

    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


def _load_panel_optimization(output_dir: Path):
    """Load the panel_optimization object from a CLI run's pickle."""
    pkl = output_dir / "biomarker_result.pkl"
    assert pkl.exists(), f"biomarker_result.pkl not produced at {pkl}"
    with open(pkl, "rb") as f:
        result = pickle.load(f)
    base = result.base_result if hasattr(result, "base_result") else result
    return base.panel_optimization


def _run_cli(synth_data_dir, output_dir, *extra_args):
    """Invoke ``raptor biomarker`` with the given extra args, asserting
    a clean exit. Returns the CliRunner result for further inspection."""
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "biomarker",
            "--counts", str(synth_data_dir / "counts.csv"),
            "--metadata", str(synth_data_dir / "metadata.csv"),
            "--group-column", "condition",
            "--output", str(output_dir),
            "--no-annotate", "--no-literature", "--no-ppi",
            *extra_args,
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, (
        f"CLI exited with code {result.exit_code}.\n"
        f"stdout:\n{result.output}"
    )
    return result


# ============================================================================
# Group 1 — --help advertises new flags
# ============================================================================
class TestBiomarkerHelpAdvertisesKneedleFlags:
    """biomarker --help must list all three new flags so users can
    discover them."""

    def test_help_lists_auto_panel_strategy(self):
        runner = CliRunner()
        result = runner.invoke(main, ["biomarker", "--help"])
        assert result.exit_code == 0
        assert "--auto-panel-strategy" in result.output

    def test_help_lists_panel_sensitivity(self):
        runner = CliRunner()
        result = runner.invoke(main, ["biomarker", "--help"])
        assert result.exit_code == 0
        assert "--panel-sensitivity" in result.output

    def test_help_lists_panel_size_strategy(self):
        runner = CliRunner()
        result = runner.invoke(main, ["biomarker", "--help"])
        assert result.exit_code == 0
        assert "--panel-size-strategy" in result.output

    def test_help_documents_kneedle_default(self):
        """Help text should make clear kneedle is the default."""
        runner = CliRunner()
        result = runner.invoke(main, ["biomarker", "--help"])
        assert result.exit_code == 0
        # Help mentions kneedle by name in the auto-panel-strategy block
        assert "kneedle" in result.output.lower()


# ============================================================================
# Group 2 — Click rejects invalid values
# ============================================================================
class TestBiomarkerRejectsInvalidValues:
    """Click's type=Choice should reject typos cleanly."""

    def test_rejects_invalid_auto_panel_strategy(self, synth_data_dir, tmp_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "biomarker",
                "--counts", str(synth_data_dir / "counts.csv"),
                "--metadata", str(synth_data_dir / "metadata.csv"),
                "--group-column", "condition",
                "--output", str(tmp_path / "out"),
                "--auto-panel-strategy", "totally_bogus_strategy",
            ],
        )
        # Click exit code 2 = usage error
        assert result.exit_code == 2
        assert "auto-panel-strategy" in result.output.lower() \
            or "invalid" in result.output.lower()

    def test_rejects_invalid_panel_size_strategy(self, synth_data_dir, tmp_path):
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "biomarker",
                "--counts", str(synth_data_dir / "counts.csv"),
                "--metadata", str(synth_data_dir / "metadata.csv"),
                "--group-column", "condition",
                "--output", str(tmp_path / "out"),
                "--panel-size-strategy", "wrong_strategy_name",
            ],
        )
        assert result.exit_code == 2


# ============================================================================
# Group 3 — Strategies produce the expected selection_method
# ============================================================================
class TestSelectionMethodEndToEnd:
    """Each --auto-panel-strategy value should produce the matching
    selection_method in the pickled result.

    These tests run a real (small) discovery, so they're slower (~30s each).
    They're the most valuable coverage we have because they test the full
    CLI -> discover_biomarkers() -> PanelOptimizationResult chain."""

    def test_first_drop_produces_first_drop_label(
        self, synth_data_dir, tmp_path,
    ):
        out = tmp_path / "first_drop_out"
        _run_cli(
            synth_data_dir, out,
            "--auto-panel-strategy", "first_drop",
            # per_fold so first_drop fires (consensus would relabel)
            "--panel-size-strategy", "per_fold",
        )
        po = _load_panel_optimization(out)
        assert po.selection_method == "first_drop", (
            f"Expected 'first_drop' selection_method, got {po.selection_method!r}"
        )

    def test_argmax_produces_argmax_label(self, synth_data_dir, tmp_path):
        out = tmp_path / "argmax_out"
        _run_cli(
            synth_data_dir, out,
            "--auto-panel-strategy", "argmax",
            "--panel-size-strategy", "per_fold",
        )
        po = _load_panel_optimization(out)
        assert po.selection_method == "argmax", (
            f"Expected 'argmax' selection_method, got {po.selection_method!r}"
        )

    def test_consensus_default_produces_consensus_pinned_label(
        self, synth_data_dir, tmp_path,
    ):
        """When --panel-size-strategy is unspecified, the default is
        consensus, and the final-panel result should be labeled
        'consensus_pinned' (not 'user_specified'). Regression for the
        labeling bug discovered in the dashboard screenshot session."""
        out = tmp_path / "consensus_default_out"
        _run_cli(synth_data_dir, out)  # no flags = use defaults
        po = _load_panel_optimization(out)
        assert po.selection_method == "consensus_pinned", (
            f"Default consensus path should label 'consensus_pinned', "
            f"got {po.selection_method!r}. The dashboard annotation "
            f"will read 'user-specified' instead of 'consensus-pinned' "
            f"if this regresses."
        )

    def test_kneedle_per_fold_produces_kneedle_or_fallback(
        self, synth_data_dir, tmp_path,
    ):
        """Per-fold kneedle. On a saturated synthetic curve the result
        will be 'argmax_fallback' (no curvature). On a curve with a
        real knee it will be 'kneedle'. Either is acceptable here —
        we're testing that the path runs and produces ONE of the
        kneedle-branch labels, not first_drop or argmax."""
        out = tmp_path / "kneedle_per_fold_out"
        _run_cli(
            synth_data_dir, out,
            "--auto-panel-strategy", "kneedle",
            "--panel-size-strategy", "per_fold",
        )
        po = _load_panel_optimization(out)
        assert po.selection_method in ("kneedle", "argmax_fallback"), (
            f"Per-fold kneedle path should produce 'kneedle' or "
            f"'argmax_fallback', got {po.selection_method!r}"
        )


# ============================================================================
# Group 4 — Explicit --panel-size overrides auto-detection
# ============================================================================
class TestExplicitPanelSizeOverridesAutoDetection:
    """When the user sets --panel-size, kneedle flags become inert and
    the result is labeled 'user_specified'."""

    def test_panel_size_overrides_kneedle(self, synth_data_dir, tmp_path):
        out = tmp_path / "explicit_size_out"
        _run_cli(
            synth_data_dir, out,
            "--panel-size", "5",
            "--auto-panel-strategy", "kneedle",  # ignored
            "--panel-size-strategy", "consensus",  # ignored
        )
        po = _load_panel_optimization(out)
        assert po.optimal_size == 5, (
            f"Expected explicit size=5, got {po.optimal_size}"
        )
        assert po.selection_method == "user_specified", (
            f"Explicit --panel-size should produce 'user_specified', "
            f"got {po.selection_method!r}"
        )


# ============================================================================
# Group 5 — Sensitivity parameter is accepted
# ============================================================================
class TestPanelSensitivity:
    """--panel-sensitivity accepts floats and threads through to kneedle."""

    @pytest.mark.parametrize("sensitivity", [0.5, 1.0, 2.0])
    def test_sensitivity_values_accepted(
        self, synth_data_dir, tmp_path, sensitivity,
    ):
        out = tmp_path / f"sens_{sensitivity}_out"
        _run_cli(
            synth_data_dir, out,
            "--panel-sensitivity", str(sensitivity),
            "--auto-panel-strategy", "kneedle",
        )
        po = _load_panel_optimization(out)
        # We don't assert on optimal_size — sensitivity affects which knee
        # kneedle picks but the saturated curve will likely route to
        # argmax_fallback regardless. Just confirm the call succeeded
        # and produced a valid selection_method.
        assert po.selection_method in (
            "kneedle", "argmax_fallback", "consensus_pinned",
        ), f"unexpected selection_method {po.selection_method!r}"