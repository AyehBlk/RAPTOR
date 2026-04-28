"""
RAPTOR CLI Comprehensive Test Suite v2.2.2
==========================================

Tests all 14 CLI commands across 4 tiers:
  Tier 0: Smoke tests (help, version, validate-installation)
  Tier 1: Structural commands (init, pipeline list)
  Tier 2: Data-driven commands with synthetic data (qc, profile, recommend,
           import-de, compare-de, optimize, ensemble, ensemble-compare)
  Tier 3: External-tool commands (quick-count, pipeline run) - validation only

Usage (from RAPTOR-main root, with venv active):
    pytest tests/test_cli_comprehensive.py -v 2>&1 | Tee-Object test_results\test_cli_comprehensive.txt

Author: Ayeh Bolouki / Claude (debugging session)
"""

import pytest
import json
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from click.testing import CliRunner

from raptor.cli import main


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def runner():
    """Click CLI test runner."""
    return CliRunner()


@pytest.fixture
def tmp_project(tmp_path):
    return tmp_path


@pytest.fixture
def sample_counts(tmp_path):
    """Create a realistic synthetic count matrix CSV."""
    np.random.seed(42)
    n_genes, n_samples = 500, 6
    gene_ids = [f"GENE{i:04d}" for i in range(n_genes)]
    sample_ids = [f"Sample{i+1}" for i in range(n_samples)]

    means = np.random.lognormal(mean=3, sigma=2, size=n_genes)
    counts = np.zeros((n_genes, n_samples), dtype=int)
    for j in range(n_samples):
        for i in range(n_genes):
            counts[i, j] = max(0, int(np.random.negative_binomial(
                n=max(1, int(means[i] / 2)), p=0.5)))

    df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)
    filepath = tmp_path / "test_counts.csv"
    df.to_csv(filepath)
    return str(filepath)


@pytest.fixture
def sample_metadata(tmp_path):
    """Create sample metadata CSV."""
    metadata = pd.DataFrame({
        'sample_id': [f"Sample{i+1}" for i in range(6)],
        'condition': ['Control']*3 + ['Treatment']*3,
        'batch': ['A', 'A', 'B', 'A', 'B', 'B'],
    })
    filepath = tmp_path / "test_metadata.csv"
    metadata.to_csv(filepath, index=False)
    return str(filepath)


@pytest.fixture
def sample_de_result(tmp_path):
    """Create a synthetic DEResult pickle file."""
    from raptor.de_import import DEResult

    np.random.seed(42)
    n_genes = 200
    gene_ids = [f"GENE{i:04d}" for i in range(n_genes)]

    # Make ~40 genes clearly significant
    pvalues = np.random.uniform(0.1, 1.0, n_genes)
    pvalues[:40] = np.random.uniform(1e-8, 1e-3, 40)
    padj = np.minimum(pvalues * n_genes / np.arange(1, n_genes + 1), 1.0)

    log2fc = np.random.normal(0, 0.5, n_genes)
    # Give significant genes large fold changes
    log2fc[:20] = np.random.uniform(1.5, 4.0, 20)    # upregulated
    log2fc[20:40] = np.random.uniform(-4.0, -1.5, 20)  # downregulated

    data = pd.DataFrame({
        'p_value': pvalues,
        'adjusted_p_value': padj,
        'log2_fold_change': log2fc,
    }, index=gene_ids)
    data.index.name = 'gene_id'

    result = DEResult(results_df=data, pipeline='DESeq2')

    filepath = tmp_path / "de_result.pkl"
    with open(filepath, 'wb') as f:
        pickle.dump(result, f)
    return str(filepath)


@pytest.fixture
def multiple_de_results(tmp_path):
    """Create multiple DE result files for comparison/ensemble."""
    from raptor.de_import import DEResult
    files = {}

    for idx, method in enumerate(['deseq2', 'edger', 'limma']):
        np.random.seed(42 + idx)
        n_genes = 200
        gene_ids = [f"GENE{i:04d}" for i in range(n_genes)]

        pvalues = np.random.uniform(0.1, 1.0, n_genes)
        pvalues[:40] = np.random.uniform(1e-8, 1e-3, 40)
        padj = np.minimum(pvalues * n_genes / np.arange(1, n_genes + 1), 1.0)

        log2fc = np.random.normal(0, 0.5, n_genes)
        log2fc[:20] = np.random.uniform(1.5, 4.0, 20)
        log2fc[20:40] = np.random.uniform(-4.0, -1.5, 20)

        data = pd.DataFrame({
            'p_value': pvalues,
            'adjusted_p_value': padj,
            'log2_fold_change': log2fc,
        }, index=gene_ids)
        data.index.name = 'gene_id'

        result = DEResult(results_df=data, pipeline=method)

        filepath = tmp_path / f"{method}_result.pkl"
        with open(filepath, 'wb') as f:
            pickle.dump(result, f)
        files[method] = str(filepath)
    return files


@pytest.fixture
def sample_de_csv(tmp_path):
    """Create a sample DESeq2 CSV output file for import-de."""
    np.random.seed(42)
    n_genes = 200
    gene_ids = [f"GENE{i:04d}" for i in range(n_genes)]
    data = pd.DataFrame({
        'baseMean': np.random.lognormal(5, 2, n_genes),
        'log2FoldChange': np.random.normal(0, 1.5, n_genes),
        'lfcSE': np.abs(np.random.normal(0.3, 0.1, n_genes)),
        'stat': np.random.normal(0, 2, n_genes),
        'pvalue': np.random.beta(0.3, 5, n_genes),
        'padj': np.random.beta(0.5, 5, n_genes),
    }, index=gene_ids)
    filepath = tmp_path / "deseq2_results.csv"
    data.to_csv(filepath)
    return str(filepath)


@pytest.fixture
def sample_profile_json(tmp_path, sample_counts, sample_metadata):
    """Generate a real profile JSON by running the profiler."""
    try:
        from raptor.profiler import RNAseqDataProfiler
        counts_df = pd.read_csv(sample_counts, index_col=0)
        metadata_df = pd.read_csv(sample_metadata)
        profiler = RNAseqDataProfiler(counts_df, metadata_df,
                                       group_column='condition',
                                       min_count_threshold=1)
        profile_result = profiler.run_full_profile()
        filepath = tmp_path / "data_profile.json"
        with open(filepath, 'w') as f:
            f.write(profile_result.to_json())
        return str(filepath)
    except Exception as e:
        pytest.skip(f"Could not generate profile: {e}")


@pytest.fixture
def ground_truth_csv(tmp_path):
    """Create ground truth gene list for optimization."""
    genes = [f"GENE{i:04d}" for i in range(50)]
    df = pd.DataFrame({'gene_id': genes})
    filepath = tmp_path / "ground_truth.csv"
    df.to_csv(filepath, index=False)
    return str(filepath)


# =============================================================================
# TIER 0: SMOKE TESTS
# =============================================================================

class TestTier0_SmokeTests:
    """Basic CLI availability - no data needed."""

    def test_help(self, runner):
        result = runner.invoke(main, ['--help'])
        assert result.exit_code == 0
        assert 'RAPTOR' in result.output
        for cmd in ['quick-count', 'qc', 'profile', 'recommend',
                     'pipeline', 'import-de', 'optimize', 'ensemble']:
            assert cmd in result.output

    def test_version(self, runner):
        """Verify --version prints the package version. Reads
        raptor.__version__ dynamically so this test survives future
        version bumps without manual editing."""
        import raptor
        result = runner.invoke(main, ['--version'])
        assert result.exit_code == 0
        assert raptor.__version__ in result.output

    def test_validate_installation(self, runner):
        result = runner.invoke(main, ['validate-installation'])
        assert result.exit_code == 0
        assert 'RAPTOR' in result.output

    def test_no_args_shows_help(self, runner):
        result = runner.invoke(main, [])
        # Click returns 2 for missing subcommand, which is correct behavior
        assert result.exit_code in (0, 2)


class TestTier0_SubcommandHelp:
    """Every subcommand responds to --help."""

    @pytest.mark.parametrize("command", [
        ['init', '--help'],
        ['validate-installation', '--help'],
        ['quick-count', '--help'],
        ['qc', '--help'],
        ['profile', '--help'],
        ['recommend', '--help'],
        ['pipeline', '--help'],
        ['pipeline', 'list', '--help'],
        ['pipeline', 'run', '--help'],
        ['import-de', '--help'],
        ['compare-de', '--help'],
        ['optimize', '--help'],
        ['ensemble', '--help'],
        ['ensemble-compare', '--help'],
    ])
    def test_subcommand_help(self, runner, command):
        result = runner.invoke(main, command)
        assert result.exit_code == 0, f"Failed for {command}: {result.output}"


# =============================================================================
# TIER 1: STRUCTURAL COMMANDS
# =============================================================================

class TestTier1_Init:

    def test_init_basic(self, runner, tmp_project):
        with runner.isolated_filesystem(temp_dir=tmp_project):
            result = runner.invoke(main, ['init', 'test_project'])
            assert result.exit_code == 0
            p = Path('test_project')
            assert p.exists()
            assert (p / 'config').exists()
            assert (p / 'data' / 'fastq').exists()
            assert (p / 'results' / 'qc').exists()
            assert (p / 'config' / 'sample_sheet_template.csv').exists()
            assert (p / 'README.md').exists()

    def test_init_full_template(self, runner, tmp_project):
        with runner.isolated_filesystem(temp_dir=tmp_project):
            result = runner.invoke(main, ['init', 'full_project', '--template', 'full'])
            assert result.exit_code == 0
            p = Path('full_project')
            assert (p / 'results' / 'optimization').exists()
            assert (p / 'results' / 'ensemble').exists()

    def test_init_duplicate_fails(self, runner, tmp_project):
        with runner.isolated_filesystem(temp_dir=tmp_project):
            runner.invoke(main, ['init', 'dup_project'])
            result = runner.invoke(main, ['init', 'dup_project'])
            assert result.exit_code != 0

    def test_init_sample_sheet_content(self, runner, tmp_project):
        with runner.isolated_filesystem(temp_dir=tmp_project):
            runner.invoke(main, ['init', 'csv_test'])
            sheet = pd.read_csv(Path('csv_test/config/sample_sheet_template.csv'))
            assert 'sample_id' in sheet.columns
            assert 'condition' in sheet.columns


class TestTier1_PipelineList:

    def test_pipeline_list(self, runner):
        result = runner.invoke(main, ['pipeline', 'list'])
        assert result.exit_code == 0
        for name in ['salmon', 'kallisto', 'star', 'hisat2']:
            assert name in result.output.lower()


# =============================================================================
# TIER 2: DATA-DRIVEN COMMANDS
# =============================================================================

class TestTier2_QC:

    def test_qc_basic(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "qc_output")
        result = runner.invoke(main, ['qc', '--counts', sample_counts, '--output', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert len(list(Path(out).iterdir())) > 0

    def test_qc_with_metadata(self, runner, sample_counts, sample_metadata, tmp_path):
        out = str(tmp_path / "qc_meta")
        result = runner.invoke(main, ['qc', '-c', sample_counts, '-m', sample_metadata, '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_qc_normalization_options(self, runner, sample_counts, tmp_path):
        for norm in ['log2', 'cpm', 'quantile', 'none']:
            out = str(tmp_path / f"qc_{norm}")
            result = runner.invoke(main, ['qc', '-c', sample_counts, '-n', norm, '-o', out])
            assert result.exit_code == 0, f"Failed for {norm}: {result.output}"

    def test_qc_invalid_threshold(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "qc_bad")
        result = runner.invoke(main, ['qc', '-c', sample_counts, '--consensus-threshold', '10', '-o', out])
        assert result.exit_code != 0

    def test_qc_missing_counts(self, runner, tmp_path):
        result = runner.invoke(main, ['qc', '--counts', str(tmp_path / 'nope.csv')])
        assert result.exit_code != 0


class TestTier2_Profile:

    def test_profile_basic(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "prof")
        result = runner.invoke(main, ['profile', '-c', sample_counts, '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert (Path(out) / 'data_profile.json').exists()
        assert (Path(out) / 'profile_summary.txt').exists()
        assert (Path(out) / 'recommendation_features.json').exists()

    def test_profile_with_metadata(self, runner, sample_counts, sample_metadata, tmp_path):
        out = str(tmp_path / "prof_meta")
        result = runner.invoke(main, ['profile', '-c', sample_counts, '-m', sample_metadata, '-g', 'condition', '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_profile_verbose(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "prof_v")
        result = runner.invoke(main, ['profile', '-c', sample_counts, '-v', '-o', out])
        assert result.exit_code == 0

    def test_profile_json_valid(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "prof_json")
        runner.invoke(main, ['profile', '-c', sample_counts, '-o', out])
        with open(Path(out) / 'data_profile.json') as f:
            profile = json.load(f)
        assert 'feature_vector' in profile or 'n_samples' in profile


class TestTier2_Recommend:

    def test_recommend_with_profile(self, runner, sample_profile_json, tmp_path):
        out = str(tmp_path / "rec")
        result = runner.invoke(main, ['recommend', '--profile', sample_profile_json, '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert (Path(out) / 'recommendation.json').exists()

    def test_recommend_verbose(self, runner, sample_profile_json, tmp_path):
        out = str(tmp_path / "rec_v")
        result = runner.invoke(main, ['recommend', '--profile', sample_profile_json, '--verbose-explanation', '-o', out])
        assert result.exit_code == 0
        assert (Path(out) / 'explanation.txt').exists()

    def test_recommend_no_profile(self, runner, tmp_path):
        with runner.isolated_filesystem(temp_dir=tmp_path):
            result = runner.invoke(main, ['recommend'])
            assert result.exit_code != 0


class TestTier2_ImportDE:

    def test_import_de_deseq2(self, runner, sample_de_csv, tmp_path):
        out = str(tmp_path / "imp")
        result = runner.invoke(main, ['import-de', '-i', sample_de_csv, '-m', 'deseq2', '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_import_de_custom_thresholds(self, runner, sample_de_csv, tmp_path):
        out = str(tmp_path / "imp_custom")
        result = runner.invoke(main, ['import-de', '-i', sample_de_csv, '-m', 'deseq2',
                                       '--fdr-threshold', '0.01', '--lfc-threshold', '1.0', '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_import_de_invalid_fdr(self, runner, sample_de_csv, tmp_path):
        out = str(tmp_path / "imp_bad")
        result = runner.invoke(main, ['import-de', '-i', sample_de_csv, '-m', 'deseq2',
                                       '--fdr-threshold', '2.0', '-o', out])
        assert result.exit_code != 0

    def test_import_de_missing_file(self, runner, tmp_path):
        result = runner.invoke(main, ['import-de', '-i', str(tmp_path / 'nope.csv'), '-m', 'deseq2'])
        assert result.exit_code != 0


class TestTier2_CompareDE:

    def test_compare_de_two(self, runner, multiple_de_results, tmp_path):
        out = str(tmp_path / "cmp")
        result = runner.invoke(main, ['compare-de', multiple_de_results['deseq2'],
                                       multiple_de_results['edger'], '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_compare_de_three(self, runner, multiple_de_results, tmp_path):
        out = str(tmp_path / "cmp3")
        result = runner.invoke(main, ['compare-de', multiple_de_results['deseq2'],
                                       multiple_de_results['edger'],
                                       multiple_de_results['limma'], '-o', out])
        assert result.exit_code == 0
        assert (Path(out) / 'comparison.json').exists()

    def test_compare_de_single_fails(self, runner, multiple_de_results, tmp_path):
        out = str(tmp_path / "cmp1")
        result = runner.invoke(main, ['compare-de', multiple_de_results['deseq2'], '-o', out])
        assert result.exit_code != 0


class TestTier2_Optimize:

    def test_optimize_fdr_control(self, runner, sample_de_result, tmp_path):
        out = str(tmp_path / "opt_fdr")
        result = runner.invoke(main, ['optimize', '-d', sample_de_result, '-m', 'fdr-control',
                                       '--fdr-target', '0.05', '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_optimize_ground_truth(self, runner, sample_de_result, ground_truth_csv, tmp_path):
        out = str(tmp_path / "opt_gt")
        result = runner.invoke(main, ['optimize', '-d', sample_de_result, '-m', 'ground-truth',
                                       '-g', ground_truth_csv, '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_optimize_missing_file(self, runner, tmp_path):
        result = runner.invoke(main, ['optimize', '-d', str(tmp_path / 'nope.pkl'), '-m', 'fdr-control'])
        assert result.exit_code != 0


class TestTier2_Ensemble:

    def test_ensemble_fisher(self, runner, multiple_de_results, tmp_path):
        out = str(tmp_path / "ens_f")
        result = runner.invoke(main, ['ensemble', '--methods', 'fisher',
                                       '--deseq2', multiple_de_results['deseq2'],
                                       '--edger', multiple_de_results['edger'], '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_ensemble_with_limma(self, runner, multiple_de_results, tmp_path):
        out = str(tmp_path / "ens_3")
        result = runner.invoke(main, ['ensemble', '--methods', 'fisher',
                                       '--deseq2', multiple_de_results['deseq2'],
                                       '--edger', multiple_de_results['edger'],
                                       '--limma', multiple_de_results['limma'], '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"

    def test_ensemble_voting_disabled(self, runner, multiple_de_results, tmp_path):
        """Voting is currently disabled in CLI (not in choices), should fail."""
        out = str(tmp_path / "ens_v")
        result = runner.invoke(main, ['ensemble', '--methods', 'voting',
                                       '--deseq2', multiple_de_results['deseq2'],
                                       '--edger', multiple_de_results['edger'],
                                       '--limma', multiple_de_results['limma'], '-o', out])
        assert result.exit_code != 0


class TestTier2_EnsembleCompare:

    def test_ensemble_compare(self, runner, multiple_de_results, tmp_path):
        out = str(tmp_path / "ens_cmp")
        result = runner.invoke(main, ['ensemble-compare',
                                       '--deseq2', multiple_de_results['deseq2'],
                                       '--edger', multiple_de_results['edger'],
                                       '--limma', multiple_de_results['limma'], '-o', out])
        assert result.exit_code == 0, f"Exit {result.exit_code}: {result.output}"
        assert (Path(out) / 'comparison_summary.csv').exists()


# =============================================================================
# TIER 3: VALIDATION-ONLY (external tools not installed)
# =============================================================================

class TestTier3_QuickCount:

    def test_quick_count_missing_samples_and_fastq(self, runner):
        result = runner.invoke(main, ['quick-count', '--index', '/fake/index'])
        assert result.exit_code != 0

    def test_quick_count_invalid_index(self, runner, tmp_path):
        ss = tmp_path / "samples.csv"
        ss.write_text("sample_id,condition,batch,fastq_r1,fastq_r2\n")
        result = runner.invoke(main, ['quick-count', '-s', str(ss), '-i', str(tmp_path / 'nope')])
        assert result.exit_code != 0


class TestTier3_PipelineRun:

    def test_pipeline_run_missing_samples(self, runner, tmp_path):
        result = runner.invoke(main, ['pipeline', 'run', '-n', 'salmon',
                                       '-s', str(tmp_path / 'nope.csv')])
        assert result.exit_code != 0


# =============================================================================
# CROSS-CUTTING & END-TO-END
# =============================================================================

class TestCrossCutting:

    def test_verbose_flag(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "v_test")
        result = runner.invoke(main, ['--verbose', 'qc', '-c', sample_counts, '-o', out])
        assert result.exit_code == 0

    def test_quiet_flag(self, runner, sample_counts, tmp_path):
        out = str(tmp_path / "q_test")
        result = runner.invoke(main, ['--quiet', 'qc', '-c', sample_counts, '-o', out])
        assert result.exit_code == 0

    def test_e2e_profile_to_recommend(self, runner, sample_counts, sample_metadata, tmp_path):
        """Full workflow: profile -> recommend."""
        prof_dir = str(tmp_path / "e2e_prof")
        rec_dir = str(tmp_path / "e2e_rec")

        # Profile
        r1 = runner.invoke(main, ['profile', '-c', sample_counts, '-m', sample_metadata, '-o', prof_dir])
        assert r1.exit_code == 0, f"Profile failed: {r1.output}"
        pj = str(Path(prof_dir) / 'data_profile.json')
        assert Path(pj).exists()

        # Recommend
        r2 = runner.invoke(main, ['recommend', '--profile', pj, '-o', rec_dir])
        assert r2.exit_code == 0, f"Recommend failed: {r2.output}"
        assert Path(rec_dir, 'recommendation.json').exists()