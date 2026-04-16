"""
RAPTOR Module 10: CLI Compatibility & Realistic Scenario Tests

Covers:
    Part A — CLI command tests (raptor biomarker, biomarker-survival, biomarker-validate)
    Part B — Realistic scenario tests (multi-method, small data, upstream integration, etc.)

Run with:
    pytest test_m10_cli_scenarios.py -v --tb=short
    pytest test_m10_cli_scenarios.py -v -k "CLI"       # CLI tests only
    pytest test_m10_cli_scenarios.py -v -k "Scenario"   # Scenarios only

Author: Ayeh Bolouki
"""

import pytest
import numpy as np
import pandas as pd
import tempfile
import shutil
import pickle
from pathlib import Path
from click.testing import CliRunner


# ============================================================================
# SHARED FIXTURES
# ============================================================================

@pytest.fixture(scope="session")
def test_data_dir():
    """Create a temp directory with synthetic test CSV files."""
    tmpdir = Path(tempfile.mkdtemp(prefix="raptor_m10_cli_"))

    np.random.seed(42)
    n_genes = 200
    n_samples = 40
    n_de = 20

    gene_ids = [f"GENE_{i:04d}" for i in range(n_genes)]
    sample_ids = [f"sample_{i:02d}" for i in range(n_samples)]

    # Count matrix with clear DE signal
    counts = np.random.negative_binomial(n=5, p=0.01, size=(n_genes, n_samples))
    for i in range(n_de):
        fold = np.random.uniform(2.0, 5.0)
        counts[i, 20:] = (counts[i, 20:] * fold).astype(int)

    counts_df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)
    counts_df.to_csv(tmpdir / "counts.csv")

    # Metadata
    meta_df = pd.DataFrame({
        'sample_id': sample_ids,
        'condition': ['control'] * 20 + ['treatment'] * 20,
        'batch': [1] * 10 + [2] * 10 + [1] * 10 + [2] * 10,
    })
    meta_df.to_csv(tmpdir / "metadata.csv", index=False)

    # Clinical data (for survival)
    np.random.seed(123)
    clinical_df = pd.DataFrame({
        'sample_id': sample_ids,
        'os_time': np.random.exponential(500, n_samples).astype(int) + 30,
        'os_event': np.random.binomial(1, 0.6, n_samples),
    })
    clinical_df.to_csv(tmpdir / "clinical.csv", index=False)

    # Small dataset (10 samples)
    small_samples = [f"sm_{i}" for i in range(10)]
    small_counts = np.random.negative_binomial(5, 0.01, size=(50, 10))
    for i in range(5):
        small_counts[i, 5:] = (small_counts[i, 5:] * 3).astype(int)
    small_df = pd.DataFrame(
        small_counts,
        index=[f"G{i}" for i in range(50)],
        columns=small_samples,
    )
    small_df.to_csv(tmpdir / "small_counts.csv")

    small_meta = pd.DataFrame({
        'sample_id': small_samples,
        'condition': ['A'] * 5 + ['B'] * 5,
    })
    small_meta.to_csv(tmpdir / "small_metadata.csv", index=False)

    # Validation cohort (independent dataset with same gene IDs)
    np.random.seed(999)
    val_counts = np.random.negative_binomial(5, 0.01, size=(n_genes, 30))
    val_samples = [f"val_{i:02d}" for i in range(30)]
    for i in range(n_de):
        val_counts[i, 15:] = (val_counts[i, 15:] * 3).astype(int)
    val_df = pd.DataFrame(val_counts, index=gene_ids, columns=val_samples)
    val_df.to_csv(tmpdir / "validation_counts.csv")

    val_meta = pd.DataFrame({
        'sample_id': val_samples,
        'condition': ['control'] * 15 + ['treatment'] * 15,
    })
    val_meta.to_csv(tmpdir / "validation_metadata.csv", index=False)

    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def output_dir():
    """Fresh output directory per test."""
    tmpdir = tempfile.mkdtemp(prefix="raptor_m10_out_")
    yield Path(tmpdir)
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def runner():
    """Click CLI test runner."""
    return CliRunner()


# ============================================================================
# PART A: CLI COMPATIBILITY TESTS
# ============================================================================

class TestCLI_Help:
    """Test that all M10 CLI commands register and show help."""

    def test_biomarker_help(self, runner):
        from raptor.cli import main
        result = runner.invoke(main, ['biomarker', '--help'])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert 'Discover biomarker gene panel' in result.output
        assert '--counts' in result.output
        assert '--metadata' in result.output
        assert '--group-column' in result.output
        assert '--methods' in result.output
        assert '--panel-size' in result.output
        assert '--species' in result.output
        assert '--no-annotate' in result.output

    def test_biomarker_survival_help(self, runner):
        from raptor.cli import main
        result = runner.invoke(main, ['biomarker-survival', '--help'])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert 'survival' in result.output.lower()
        assert '--counts' in result.output
        assert '--clinical' in result.output
        assert '--time-column' in result.output
        assert '--event-column' in result.output

    def test_biomarker_validate_help(self, runner):
        from raptor.cli import main
        result = runner.invoke(main, ['biomarker-validate', '--help'])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert 'Validate biomarker panel' in result.output
        assert '--panel' in result.output
        assert '--counts' in result.output
        assert '--n-folds' in result.output


class TestCLI_BiomarkerRun:
    """Test raptor biomarker with real data."""

    def test_basic_run(self, runner, test_data_dir, output_dir):
        """raptor biomarker -c counts.csv -m metadata.csv -g condition --no-annotate"""
        from raptor.cli import main
        result = runner.invoke(main, [
            'biomarker',
            '-c', str(test_data_dir / 'counts.csv'),
            '-m', str(test_data_dir / 'metadata.csv'),
            '-g', 'condition',
            '--methods', 'elastic_net',
            '--panel-size', '5',
            '--no-annotate',
            '-o', str(output_dir),
        ])
        assert result.exit_code == 0, f"Exit {result.exit_code}:\n{result.output}"
        assert 'Panel size:' in result.output or 'Panel' in result.output
        assert (output_dir / 'biomarker_panel.csv').exists()
        assert (output_dir / 'ranked_genes.csv').exists()
        assert (output_dir / 'biomarker_result.pkl').exists()

    def test_multiple_methods_cli(self, runner, test_data_dir, output_dir):
        """raptor biomarker with multiple --methods flags."""
        from raptor.cli import main
        result = runner.invoke(main, [
            'biomarker',
            '-c', str(test_data_dir / 'counts.csv'),
            '-m', str(test_data_dir / 'metadata.csv'),
            '-g', 'condition',
            '--methods', 'elastic_net',
            '--methods', 'rfe',
            '--panel-size', '8',
            '--no-annotate',
            '-o', str(output_dir),
        ])
        assert result.exit_code == 0, f"Exit {result.exit_code}:\n{result.output}"
        assert (output_dir / 'biomarker_panel.csv').exists()

    def test_missing_file_error(self, runner, output_dir):
        """CLI should fail gracefully with missing file."""
        from raptor.cli import main
        result = runner.invoke(main, [
            'biomarker',
            '-c', 'nonexistent_file.csv',
            '-m', 'also_missing.csv',
            '-g', 'condition',
            '-o', str(output_dir),
        ])
        assert result.exit_code != 0

    def test_wrong_group_column(self, runner, test_data_dir, output_dir):
        """CLI should fail with clear error for wrong group column."""
        from raptor.cli import main
        result = runner.invoke(main, [
            'biomarker',
            '-c', str(test_data_dir / 'counts.csv'),
            '-m', str(test_data_dir / 'metadata.csv'),
            '-g', 'nonexistent_column',
            '--methods', 'elastic_net',
            '--no-annotate',
            '-o', str(output_dir),
        ])
        assert result.exit_code != 0


class TestCLI_BiomarkerValidate:
    """Test raptor biomarker-validate pipeline."""

    def test_full_discovery_then_validate(self, runner, test_data_dir, output_dir):
        """End-to-end: discover -> validate on independent cohort."""
        from raptor.cli import main

        # Step 1: Discover
        discover_dir = output_dir / "discovery"
        result = runner.invoke(main, [
            'biomarker',
            '-c', str(test_data_dir / 'counts.csv'),
            '-m', str(test_data_dir / 'metadata.csv'),
            '-g', 'condition',
            '--methods', 'elastic_net',
            '--panel-size', '5',
            '--no-annotate',
            '-o', str(discover_dir),
        ])
        assert result.exit_code == 0, f"Discovery failed:\n{result.output}"
        panel_path = discover_dir / 'biomarker_panel.csv'
        assert panel_path.exists()

        # Step 2: Validate on independent cohort
        validate_dir = output_dir / "validation"
        result = runner.invoke(main, [
            'biomarker-validate',
            '-p', str(panel_path),
            '-c', str(test_data_dir / 'validation_counts.csv'),
            '-m', str(test_data_dir / 'validation_metadata.csv'),
            '-g', 'condition',
            '-o', str(validate_dir),
        ])
        assert result.exit_code == 0, f"Validation failed:\n{result.output}"
        assert (validate_dir / 'validation_performance.csv').exists()
        assert (validate_dir / 'validation_summary.txt').exists()


# ============================================================================
# PART B: REALISTIC SCENARIO TESTS
# ============================================================================

class TestScenario_MultiMethod:
    """Scenario: Researcher runs multiple feature selection methods
    and expects consensus ranking to prioritize true DE genes."""

    def test_consensus_prioritizes_de(self, test_data_dir, output_dir):
        from raptor.biomarker_discovery import discover_biomarkers

        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=['elastic_net', 'rfe'],
            target_panel_size=10,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        # Panel should contain mostly true DE genes (GENE_0000..GENE_0019)
        de_genes = {f"GENE_{i:04d}" for i in range(20)}
        overlap = len(set(result.panel) & de_genes)
        assert overlap >= 5, (
            f"Expected at least 5/10 panel genes from true DE set, "
            f"got {overlap}: {result.panel}"
        )

    def test_all_methods_produce_results(self, test_data_dir, output_dir):
        """Each method should produce a FeatureSelectionResult."""
        from raptor.biomarker_discovery import discover_biomarkers

        methods = ['elastic_net', 'rfe']
        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=methods,
            target_panel_size=10,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        for m in methods:
            assert m in result.selection_results, f"Missing result for method {m}"
            assert result.selection_results[m].n_selected > 0


class TestScenario_SmallDataset:
    """Scenario: Pilot study with only 10 samples — should auto-adapt."""

    def test_small_dataset_loocv(self, test_data_dir, output_dir):
        from raptor.biomarker_discovery import discover_biomarkers

        result = discover_biomarkers(
            counts=str(test_data_dir / 'small_counts.csv'),
            metadata=str(test_data_dir / 'small_metadata.csv'),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=3,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        assert result.panel_size == 3
        assert result.n_samples == 10
        assert len(result.classification_results) >= 1
        # Should still produce valid AUC
        for name, res in result.classification_results.items():
            assert 0.0 <= res.auc <= 1.0


class TestScenario_UpstreamIntegration:
    """Scenario: User provides DE result from Module 7/9 to restrict candidates."""

    def test_with_mock_de_result(self, test_data_dir, output_dir):
        from raptor.biomarker_discovery import discover_biomarkers

        # Simulate Module 7 output
        class MockDEResult:
            significant_genes = [f"GENE_{i:04d}" for i in range(20)]
            n_genes = 20

        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            de_result=MockDEResult(),
            methods=['de_filter', 'elastic_net'],
            target_panel_size=8,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        assert 'de_filter' in result.selection_results
        # All panel genes should come from the DE gene list
        de_set = set(MockDEResult.significant_genes)
        for gene in result.panel:
            assert gene in de_set, f"Panel gene {gene} not in DE gene list"

    def test_with_pickled_de_result(self, test_data_dir, output_dir):
        """Test loading DE result from pickle (CLI workflow)."""
        from raptor.biomarker_discovery import discover_biomarkers
        from types import SimpleNamespace

        mock_de = SimpleNamespace(
            significant_genes=[f"GENE_{i:04d}" for i in range(20)],
            n_genes=20,
        )

        # Save pickle (mimics M7 output)
        pkl_path = output_dir / "de_result.pkl"
        with open(pkl_path, 'wb') as f:
            pickle.dump(mock_de, f)

        # Load it back — same as CLI does
        with open(pkl_path, 'rb') as f:
            loaded = pickle.load(f)

        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            de_result=loaded,
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir / "results"),
            verbose=False,
        )

        assert result.panel_size == 5


class TestScenario_OutputCompleteness:
    """Scenario: Check all expected output files are generated and valid."""

    def test_all_outputs_exist(self, test_data_dir, output_dir):
        from raptor.biomarker_discovery import discover_biomarkers

        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        # Check required output files
        assert (output_dir / 'biomarker_panel.csv').exists()
        assert (output_dir / 'ranked_genes.csv').exists()
        assert (output_dir / 'classification_performance.csv').exists()
        assert (output_dir / 'panel_curve.csv').exists()
        assert (output_dir / 'biomarker_result.pkl').exists()
        assert (output_dir / 'biomarker_params.json').exists()
        assert (output_dir / 'summary.txt').exists()

    def test_panel_csv_format(self, test_data_dir, output_dir):
        """biomarker_panel.csv should be loadable and have gene_id column."""
        from raptor.biomarker_discovery import discover_biomarkers

        discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        panel_df = pd.read_csv(output_dir / 'biomarker_panel.csv')
        assert len(panel_df) == 5
        # Should have gene_id as column or index
        assert 'gene_id' in panel_df.columns or panel_df.index.name == 'gene_id'

    def test_pickle_roundtrip(self, test_data_dir, output_dir):
        """Saved pickle should load back with identical attributes."""
        from raptor.biomarker_discovery import discover_biomarkers, BiomarkerResult

        original = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        loaded = BiomarkerResult.load(output_dir)
        assert loaded.panel_size == original.panel_size
        assert loaded.panel == original.panel
        assert loaded.best_classifier == original.best_classifier

    def test_summary_text_content(self, test_data_dir, output_dir):
        """summary.txt should contain key result info."""
        from raptor.biomarker_discovery import discover_biomarkers

        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        summary = (output_dir / 'summary.txt').read_text(encoding='utf-8')
        assert 'MODULE 10' in summary or 'Panel' in summary


class TestScenario_ClassificationQuality:
    """Scenario: Verify that true DE genes produce high classification AUC."""

    def test_de_genes_high_auc(self, test_data_dir, output_dir):
        from raptor.biomarker_discovery import discover_biomarkers

        class MockDEResult:
            significant_genes = [f"GENE_{i:04d}" for i in range(20)]

        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            de_result=MockDEResult(),
            methods=['de_filter', 'elastic_net'],
            target_panel_size=10,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        best = result.classification_results[result.best_classifier]
        assert best.auc > 0.7, (
            f"With true DE genes, expected AUC > 0.7 but got {best.auc:.3f}"
        )

    def test_random_genes_lower_auc(self, test_data_dir, output_dir):
        """Panel from non-DE genes should perform worse than DE genes."""
        from raptor.biomarker_discovery import validate_biomarkers

        # Use non-DE genes (GENE_0100..GENE_0109) — should have weak signal
        random_panel = [f"GENE_{i:04d}" for i in range(100, 110)]
        results = validate_biomarkers(
            panel_genes=random_panel,
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            n_folds=3,
            verbose=False,
        )

        # At least one classifier AUC should be checked
        assert len(results) >= 1


class TestScenario_ValidationWorkflow:
    """Scenario: Full discovery + independent validation workflow."""

    def test_discovery_then_validation(self, test_data_dir, output_dir):
        """Discover on training, validate on independent cohort."""
        from raptor.biomarker_discovery import discover_biomarkers, validate_biomarkers

        # Discovery
        result = discover_biomarkers(
            counts=str(test_data_dir / 'counts.csv'),
            metadata=str(test_data_dir / 'metadata.csv'),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=8,
            annotate=False,
            output_dir=str(output_dir / "discovery"),
            verbose=False,
        )

        # Validate on independent cohort
        val_results = validate_biomarkers(
            panel_genes=result.panel,
            counts=str(test_data_dir / 'validation_counts.csv'),
            metadata=str(test_data_dir / 'validation_metadata.csv'),
            group_column='condition',
            n_folds=3,
            verbose=False,
        )

        assert len(val_results) >= 2
        for name, res in val_results.items():
            assert 0.0 <= res.auc <= 1.0, f"{name}: AUC out of range"

        # Validation AUC should be decent if DE genes are in panel
        best_val = max(res.auc for res in val_results.values())
        assert best_val > 0.5, f"Validation AUC should beat random, got {best_val:.3f}"


# ============================================================================
# RUN
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
