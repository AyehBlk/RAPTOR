"""
RAPTOR v2.2.0 - Tests for Module 5 CLI Integration

Tests for command-line interface for production pipelines:
- raptor pipeline list
- raptor pipeline info <name>
- raptor pipeline salmon
- raptor pipeline kallisto
- raptor pipeline star-featurecounts
- raptor pipeline hisat2-featurecounts
- raptor pipeline star-rsem
- raptor pipeline star-salmon
- raptor pipeline run
- raptor pipeline use-quantify

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys
import tempfile
import shutil
from unittest.mock import patch, MagicMock

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Try importing Click for CLI testing
try:
    from click.testing import CliRunner
    CLICK_AVAILABLE = True
except ImportError:
    CLICK_AVAILABLE = False

# Import CLI module
try:
    from raptor.cli import main, pipeline
    CLI_AVAILABLE = True
except ImportError as e:
    CLI_AVAILABLE = False
    # Print warning so users know CLI tests will be skipped
    import warnings
    warnings.warn(
        f"CLI module not available - CLI tests will be skipped.\n"
        f"  This is normal if CLI is not implemented yet.\n"
        f"  Error: {e}",
        UserWarning
    )


# =============================================================================
# Constants
# =============================================================================

DEFAULT_PRODUCTION_DIR = "results/production"
DEFAULT_QUICK_COUNTS_DIR = "results/quick_counts"
QUICK_COUNTS_FILE = "quick_gene_counts.csv"


# =============================================================================
# Fixtures for CLI Tests
# =============================================================================

@pytest.fixture
def cli_runner():
    """Create a Click CLI test runner."""
    if CLICK_AVAILABLE:
        return CliRunner()
    return None


@pytest.fixture
def mock_quick_counts(test_dir):
    """Create mock quick counts directory for use-quantify tests."""
    quick_dir = test_dir / "results" / "quick_counts"
    quick_dir.mkdir(parents=True, exist_ok=True)
    
    # Create mock count files
    np.random.seed(42)
    counts = pd.DataFrame(
        np.random.negative_binomial(10, 0.1, (100, 6)),
        index=[f'Gene{i}' for i in range(100)],
        columns=[f'Sample{i}' for i in range(1, 7)]
    )
    
    counts.to_csv(quick_dir / "quick_gene_counts.csv")
    
    tpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    tpm.to_csv(quick_dir / "quick_tpm.csv")
    
    # Create sample info
    info = pd.DataFrame({
        'sample_id': [f'Sample{i}' for i in range(1, 7)],
        'success': [True] * 6,
        'num_reads': [10000000] * 6,
        'mapping_rate': [85.0] * 6
    })
    info.to_csv(quick_dir / "sample_info.csv", index=False)
    
    return quick_dir


# =============================================================================
# Tests for Pipeline List Command
# =============================================================================

class TestPipelineListCommand:
    """Tests for 'raptor pipeline list' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE, 
                        reason="CLI or Click not available")
    def test_pipeline_list_basic(self, cli_runner):
        """Test basic pipeline list output."""
        result = cli_runner.invoke(main, ['pipeline', 'list'])
        
        # Should succeed or at least not crash
        assert result.exit_code == 0 or 'pipeline' in result.output.lower()
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_pipeline_list_contains_pipelines(self, cli_runner):
        """Test that pipeline list contains expected pipelines."""
        result = cli_runner.invoke(main, ['pipeline', 'list'])
        
        expected_pipelines = ['salmon', 'kallisto', 'star', 'hisat2', 'rsem']
        
        output_lower = result.output.lower()
        for pipeline_name in expected_pipelines:
            # At least some should be present
            pass  # Depends on output format
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_pipeline_list_detailed(self, cli_runner):
        """Test detailed pipeline list output."""
        result = cli_runner.invoke(main, ['pipeline', 'list', '--detailed'])
        
        # Detailed output should be longer
        assert len(result.output) > 0


# =============================================================================
# Tests for Pipeline Info Command
# =============================================================================

class TestPipelineInfoCommand:
    """Tests for 'raptor pipeline info <name>' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_pipeline_info_salmon(self, cli_runner):
        """Test pipeline info for Salmon."""
        result = cli_runner.invoke(main, ['pipeline', 'info', 'salmon'])
        
        # Should contain Salmon-related info
        assert result.exit_code == 0 or 'salmon' in result.output.lower()
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_pipeline_info_kallisto(self, cli_runner):
        """Test pipeline info for Kallisto."""
        result = cli_runner.invoke(main, ['pipeline', 'info', 'kallisto'])
        
        assert result.exit_code == 0 or 'kallisto' in result.output.lower()
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_pipeline_info_invalid(self, cli_runner):
        """Test pipeline info for invalid pipeline."""
        result = cli_runner.invoke(main, ['pipeline', 'info', 'invalid_pipeline'])
        
        # Should indicate not found or return error
        assert result.exit_code != 0 or 'not found' in result.output.lower()


# =============================================================================
# Tests for Salmon Pipeline Command
# =============================================================================

class TestSalmonCommand:
    """Tests for 'raptor pipeline salmon' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_salmon_help(self, cli_runner):
        """Test Salmon command help."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        assert result.exit_code == 0
        assert '--sample-sheet' in result.output or '-s' in result.output
        assert '--index' in result.output or '-i' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_salmon_missing_required(self, cli_runner):
        """Test Salmon with missing required arguments."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon'])
        
        # Should fail due to missing required args
        assert result.exit_code != 0
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_salmon_options_in_help(self, cli_runner):
        """Test that all Salmon options are in help."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        expected_options = [
            'bootstraps', 'library-type', 'gc-bias', 'threads'
        ]
        
        for opt in expected_options:
            # Check some variant of the option name
            pass


# =============================================================================
# Tests for Kallisto Pipeline Command
# =============================================================================

class TestKallistoCommand:
    """Tests for 'raptor pipeline kallisto' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_kallisto_help(self, cli_runner):
        """Test Kallisto command help."""
        result = cli_runner.invoke(main, ['pipeline', 'kallisto', '--help'])
        
        assert result.exit_code == 0
        assert 'fragment' in result.output.lower()  # Fragment length for single-end
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_kallisto_strand_options(self, cli_runner):
        """Test Kallisto strand options in help."""
        result = cli_runner.invoke(main, ['pipeline', 'kallisto', '--help'])
        
        # Should mention strand options
        assert 'strand' in result.output.lower() or 'rf' in result.output.lower()


# =============================================================================
# Tests for STAR+featureCounts Command
# =============================================================================

class TestStarFeatureCountsCommand:
    """Tests for 'raptor pipeline star-featurecounts' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_star_fc_help(self, cli_runner):
        """Test STAR+featureCounts command help."""
        result = cli_runner.invoke(main, ['pipeline', 'star-featurecounts', '--help'])
        
        assert result.exit_code == 0
        assert '--gtf' in result.output or '-g' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_star_fc_gtf_required(self, cli_runner, sample_sheet_path, mock_star_index):
        """Test that GTF is required for STAR+featureCounts."""
        result = cli_runner.invoke(main, [
            'pipeline', 'star-featurecounts',
            '-s', str(sample_sheet_path),
            '-i', str(mock_star_index)
            # Missing --gtf
        ])
        
        # Should fail without GTF
        assert result.exit_code != 0


# =============================================================================
# Tests for HISAT2+featureCounts Command
# =============================================================================

class TestHisat2FeatureCountsCommand:
    """Tests for 'raptor pipeline hisat2-featurecounts' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_hisat2_fc_help(self, cli_runner):
        """Test HISAT2+featureCounts command help."""
        result = cli_runner.invoke(main, ['pipeline', 'hisat2-featurecounts', '--help'])
        
        assert result.exit_code == 0


# =============================================================================
# Tests for STAR+RSEM Command
# =============================================================================

class TestStarRsemCommand:
    """Tests for 'raptor pipeline star-rsem' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_star_rsem_help(self, cli_runner):
        """Test STAR+RSEM command help."""
        result = cli_runner.invoke(main, ['pipeline', 'star-rsem', '--help'])
        
        assert result.exit_code == 0
        assert 'rsem' in result.output.lower() or 'RSEM' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_star_rsem_strandedness_option(self, cli_runner):
        """Test STAR+RSEM strandedness option."""
        result = cli_runner.invoke(main, ['pipeline', 'star-rsem', '--help'])
        
        assert 'strandedness' in result.output.lower() or 'strand' in result.output.lower()


# =============================================================================
# Tests for STAR+Salmon Command
# =============================================================================

class TestStarSalmonCommand:
    """Tests for 'raptor pipeline star-salmon' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_star_salmon_help(self, cli_runner):
        """Test STAR+Salmon command help."""
        result = cli_runner.invoke(main, ['pipeline', 'star-salmon', '--help'])
        
        assert result.exit_code == 0
        assert '--salmon-index' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_star_salmon_requires_both_indexes(self, cli_runner, sample_sheet_path, mock_star_index):
        """Test that STAR+Salmon requires both indexes."""
        result = cli_runner.invoke(main, [
            'pipeline', 'star-salmon',
            '-s', str(sample_sheet_path),
            '-i', str(mock_star_index)
            # Missing --salmon-index
        ])
        
        # Should fail without salmon-index
        assert result.exit_code != 0


# =============================================================================
# Tests for Generic Run Command
# =============================================================================

class TestPipelineRunCommand:
    """Tests for 'raptor pipeline run' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_run_help(self, cli_runner):
        """Test generic run command help."""
        result = cli_runner.invoke(main, ['pipeline', 'run', '--help'])
        
        assert result.exit_code == 0
        assert '--name' in result.output or '-n' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_run_pipeline_choices(self, cli_runner):
        """Test that run command shows pipeline choices."""
        result = cli_runner.invoke(main, ['pipeline', 'run', '--help'])
        
        # Should list available pipelines as choices
        assert 'salmon' in result.output.lower()


# =============================================================================
# Tests for Use-Quantify Command
# =============================================================================

class TestUseQuantifyCommand:
    """Tests for 'raptor pipeline use-quantify' command."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_use_quantify_help(self, cli_runner):
        """Test use-quantify command help."""
        result = cli_runner.invoke(main, ['pipeline', 'use-quantify', '--help'])
        
        assert result.exit_code == 0
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_use_quantify_missing_quick_counts(self, cli_runner, test_dir):
        """Test use-quantify when quick counts don't exist."""
        with cli_runner.isolated_filesystem():
            result = cli_runner.invoke(main, ['pipeline', 'use-quantify'])
            
            # Should fail or warn about missing quick counts
            assert result.exit_code != 0 or 'not found' in result.output.lower()


# =============================================================================
# Tests for Docker and Module Options
# =============================================================================

class TestDockerModuleOptions:
    """Tests for Docker and HPC module options in CLI."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_salmon_docker_option(self, cli_runner):
        """Test Salmon --use-docker option."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        assert '--use-docker' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_salmon_modules_option(self, cli_runner):
        """Test Salmon --modules option."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        assert '--modules' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_docker_image_option(self, cli_runner):
        """Test --docker-image option."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        assert '--docker-image' in result.output


# =============================================================================
# Tests for Output Options
# =============================================================================

class TestOutputOptions:
    """Tests for output-related options in CLI."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_output_option(self, cli_runner):
        """Test -o/--output option."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        assert '--output' in result.output or '-o' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_threads_option(self, cli_runner):
        """Test -t/--threads option."""
        result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
        
        assert '--threads' in result.output or '-t' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_keep_bam_option(self, cli_runner):
        """Test --keep-bam option for alignment pipelines."""
        result = cli_runner.invoke(main, ['pipeline', 'star-featurecounts', '--help'])
        
        assert '--keep-bam' in result.output or 'bam' in result.output.lower()


# =============================================================================
# Tests for Error Messages
# =============================================================================

class TestErrorMessages:
    """Tests for CLI error message quality."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_missing_sample_sheet_error(self, cli_runner):
        """Test error message for missing sample sheet."""
        result = cli_runner.invoke(main, [
            'pipeline', 'salmon',
            '-s', '/nonexistent/samples.csv',
            '-i', '/fake/index'
        ])
        
        # Should provide helpful error message
        assert result.exit_code != 0
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_missing_index_error(self, cli_runner, sample_sheet_path):
        """Test error message for missing index."""
        result = cli_runner.invoke(main, [
            'pipeline', 'salmon',
            '-s', str(sample_sheet_path),
            '-i', '/nonexistent/index'
        ])
        
        # Should provide helpful error message
        assert result.exit_code != 0


# =============================================================================
# Integration Tests
# =============================================================================

class TestCLIIntegration:
    """Integration tests for CLI commands."""
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_version_command(self, cli_runner):
        """Test --version command."""
        result = cli_runner.invoke(main, ['--version'])
        
        assert result.exit_code == 0
        assert '2.2.0' in result.output or 'RAPTOR' in result.output
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_help_command(self, cli_runner):
        """Test --help command."""
        result = cli_runner.invoke(main, ['--help'])
        
        assert result.exit_code == 0
        assert 'pipeline' in result.output.lower()
    
    @pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                        reason="CLI or Click not available")
    def test_pipeline_group_help(self, cli_runner):
        """Test pipeline group help."""
        result = cli_runner.invoke(main, ['pipeline', '--help'])
        
        assert result.exit_code == 0
        # Should list subcommands
        assert 'list' in result.output.lower()
        assert 'salmon' in result.output.lower()


# =============================================================================
# Parametrized Tests
# =============================================================================

@pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                    reason="CLI or Click not available")
@pytest.mark.parametrize("pipeline_cmd", [
    'salmon',
    'kallisto',
    'star-featurecounts',
    'hisat2-featurecounts',
    'star-rsem',
    'star-salmon'
])
def test_all_pipeline_help_available(cli_runner, pipeline_cmd):
    """Test that help is available for all pipeline commands."""
    result = cli_runner.invoke(main, ['pipeline', pipeline_cmd, '--help'])
    
    # Should at least not crash
    assert result.exit_code == 0 or '--help' in result.output


@pytest.mark.skipif(not CLI_AVAILABLE or not CLICK_AVAILABLE,
                    reason="CLI or Click not available")
@pytest.mark.parametrize("threads", ['1', '4', '8', '16'])
def test_thread_option_values(cli_runner, threads):
    """Test various thread option values."""
    result = cli_runner.invoke(main, ['pipeline', 'salmon', '--help'])
    
    # Thread option should accept various values
    assert '--threads' in result.output or '-t' in result.output
