#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Module 6: Differential Expression Analysis (Python Wrapper)

This module provides Python wrappers for R-based differential expression
analysis scripts (DESeq2, edgeR, limma-voom).

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
"""

import subprocess
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class DEResult:
    """Differential expression analysis result."""
    method: str
    n_genes_tested: int
    n_significant: int
    n_upregulated: int
    n_downregulated: int
    pct_significant: float
    output_dir: Path
    results_file: Path
    summary_file: Path
    

class DEAnalysisWrapper:
    """
    Python wrapper for R-based differential expression analysis.
    
    Provides a unified interface to run DESeq2, edgeR, or limma-voom
    analysis from Python/RAPTOR.
    
    Parameters
    ----------
    method : str
        Analysis method: 'deseq2', 'edger', or 'limma'
    r_script_dir : str or Path, optional
        Directory containing R scripts (default: same as this file)
    
    Examples
    --------
    >>> from raptor.external_modules.module6_de_analysis import DEAnalysisWrapper
    >>> 
    >>> de = DEAnalysisWrapper(method='deseq2')
    >>> result = de.run(
    ...     counts='results/gene_counts.csv',
    ...     metadata='data/metadata.csv',
    ...     output_dir='results/de_analysis',
    ...     condition='treatment',
    ...     fdr=0.05
    ... )
    >>> print(f"Significant genes: {result.n_significant}")
    """
    
    VALID_METHODS = ['deseq2', 'edger', 'limma']
    
    def __init__(
        self,
        method: str = 'deseq2',
        r_script_dir: Optional[Union[str, Path]] = None
    ):
        self.method = method.lower()
        
        if self.method not in self.VALID_METHODS:
            raise ValueError(
                f"Invalid method: {method}. "
                f"Choose from: {', '.join(self.VALID_METHODS)}"
            )
        
        # Locate R scripts
        if r_script_dir is None:
            # R scripts are one level up from python_wrappers/ directory
            # Structure: module6_de_analysis/python_wrappers/run_de_analysis.py
            #           module6_de_analysis/r_scripts/run_deseq2.R
            self.r_script_dir = Path(__file__).parent.parent / 'r_scripts'
        else:
            self.r_script_dir = Path(r_script_dir)
        
        # Check R script exists
        self.r_script = self._get_r_script()
        if not self.r_script.exists():
            raise FileNotFoundError(
                f"R script not found: {self.r_script}\n"
                f"Expected location: {self.r_script_dir}"
            )
        
        logger.info(f"Initialized {self.method.upper()} wrapper")
    
    def _get_r_script(self) -> Path:
        """Get path to R script for chosen method."""
        script_map = {
            'deseq2': 'run_deseq2.R',
            'edger': 'run_edger.R',
            'limma': 'run_limma.R'
        }
        return self.r_script_dir / script_map[self.method]
    
    def run(
        self,
        counts: Union[str, Path],
        metadata: Union[str, Path],
        output_dir: Union[str, Path] = 'results/de_analysis',
        condition: str = 'condition',
        reference: Optional[str] = None,
        batch: Optional[str] = None,
        fdr: float = 0.05,
        lfc: float = 0.0,
        config: Optional[Union[str, Path]] = None,
        plots: bool = False,
        **kwargs
    ) -> DEResult:
        """
        Run differential expression analysis.
        
        Parameters
        ----------
        counts : str or Path
            Path to count matrix CSV/TSV file
        metadata : str or Path
            Path to sample metadata CSV file
        output_dir : str or Path
            Output directory for results
        condition : str
            Column name in metadata for condition
        reference : str, optional
            Reference level for comparison
        batch : str, optional
            Column name in metadata for batch correction
        fdr : float
            FDR threshold (default: 0.05)
        lfc : float
            Log2 fold change threshold (default: 0.0)
        config : str or Path, optional
            RAPTOR recommendation.yaml file
        plots : bool
            Generate QC plots (default: False)
        **kwargs
            Additional method-specific arguments
        
        Returns
        -------
        DEResult
            Analysis results with summary information
        
        Raises
        ------
        FileNotFoundError
            If input files don't exist
        RuntimeError
            If R script fails
        """
        # Validate inputs
        counts = Path(counts)
        metadata = Path(metadata)
        output_dir = Path(output_dir)
        
        if not counts.exists():
            raise FileNotFoundError(f"Count file not found: {counts}")
        if not metadata.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Running {self.method.upper()} analysis...")
        logger.info(f"  Counts: {counts}")
        logger.info(f"  Metadata: {metadata}")
        logger.info(f"  Output: {output_dir}")
        
        # Build R command
        cmd = self._build_command(
            counts=counts,
            metadata=metadata,
            output_dir=output_dir,
            condition=condition,
            reference=reference,
            batch=batch,
            fdr=fdr,
            lfc=lfc,
            config=config,
            plots=plots,
            **kwargs
        )
        
        # Run R script
        self._run_r_script(cmd)
        
        # Parse results
        result = self._parse_results(output_dir)
        
        logger.info(f"✅ Analysis complete!")
        logger.info(f"  Significant genes: {result.n_significant}")
        logger.info(f"  Upregulated: {result.n_upregulated}")
        logger.info(f"  Downregulated: {result.n_downregulated}")
        
        return result
    
    def _build_command(
        self,
        counts: Path,
        metadata: Path,
        output_dir: Path,
        condition: str,
        reference: Optional[str],
        batch: Optional[str],
        fdr: float,
        lfc: float,
        config: Optional[Path],
        plots: bool,
        **kwargs
    ) -> List[str]:
        """Build R script command with arguments."""
        cmd = [
            'Rscript',
            str(self.r_script),
            '--counts', str(counts),
            '--metadata', str(metadata),
            '--output', str(output_dir),
            '--condition', condition,
            '--fdr', str(fdr),
            '--lfc', str(lfc)
        ]
        
        # Optional arguments
        if reference:
            cmd.extend(['--reference', reference])
        if batch:
            cmd.extend(['--batch', batch])
        if config:
            cmd.extend(['--config', str(config)])
        if plots:
            cmd.append('--plots')
        
        # Add method-specific arguments
        for key, value in kwargs.items():
            # Convert Python naming to R naming
            r_key = key.replace('_', '-')
            
            # Handle boolean flags
            if isinstance(value, bool):
                if value:
                    cmd.append(f'--{r_key}')
            else:
                cmd.extend([f'--{r_key}', str(value)])
        
        return cmd
    
    def _run_r_script(self, cmd: List[str]) -> None:
        """Execute R script and handle errors."""
        logger.debug(f"Executing: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Log R output
            if result.stdout:
                logger.debug("R stdout:")
                for line in result.stdout.split('\n'):
                    if line.strip():
                        logger.debug(f"  {line}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"R script failed with exit code {e.returncode}")
            logger.error("R stderr:")
            for line in e.stderr.split('\n'):
                if line.strip():
                    logger.error(f"  {line}")
            
            raise RuntimeError(
                f"Differential expression analysis failed.\n"
                f"R script: {self.r_script}\n"
                f"Error: {e.stderr}"
            )
        except FileNotFoundError:
            raise RuntimeError(
                "Rscript not found. Please ensure R is installed and in PATH.\n"
                "Install R from: https://www.r-project.org/\n"
                "Then install required packages: Rscript install_packages.R"
            )
    
    def _parse_results(self, output_dir: Path) -> DEResult:
        """Parse analysis results from output files."""
        summary_file = output_dir / 'de_summary.json'
        results_file = output_dir / 'de_results.csv'
        
        if not summary_file.exists():
            raise RuntimeError(
                f"Summary file not found: {summary_file}\n"
                f"Analysis may have failed."
            )
        
        # Load summary
        with open(summary_file) as f:
            summary = json.load(f)
        
        # Create result object
        result = DEResult(
            method=self.method,
            n_genes_tested=summary['data']['n_genes_tested'],
            n_significant=summary['results']['n_significant'],
            n_upregulated=summary['results']['n_upregulated'],
            n_downregulated=summary['results']['n_downregulated'],
            pct_significant=summary['results']['pct_significant'],
            output_dir=output_dir,
            results_file=results_file,
            summary_file=summary_file
        )
        
        return result
    
    @staticmethod
    def check_r_packages() -> Dict[str, bool]:
        """
        Check if required R packages are installed.
        
        Returns
        -------
        dict
            Package name -> installed status
        """
        packages = [
            'DESeq2', 'edgeR', 'limma', 'apeglm', 'ashr',
            'optparse', 'jsonlite', 'yaml', 'ggplot2'
        ]
        
        r_code = f"""
        packages <- c({', '.join(f'"{p}"' for p in packages)})
        installed <- sapply(packages, requireNamespace, quietly = TRUE)
        cat(paste(installed, collapse = ','))
        """
        
        try:
            result = subprocess.run(
                ['Rscript', '-e', r_code],
                capture_output=True,
                text=True,
                check=True
            )
            
            installed = [s == 'TRUE' for s in result.stdout.strip().split(',')]
            return dict(zip(packages, installed))
            
        except (subprocess.CalledProcessError, FileNotFoundError):
            return {pkg: False for pkg in packages}


# =============================================================================
# CLI Interface
# =============================================================================

def main():
    """Command-line interface for DE analysis."""
    import argparse
    import click
    
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Module 6 - Differential Expression Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--method',
        choices=['deseq2', 'edger', 'limma'],
        default='deseq2',
        help='Analysis method (default: deseq2)'
    )
    parser.add_argument(
        '--counts',
        required=True,
        help='Count matrix CSV/TSV file'
    )
    parser.add_argument(
        '--metadata',
        required=True,
        help='Sample metadata CSV file'
    )
    parser.add_argument(
        '--output',
        default='results/de_analysis',
        help='Output directory (default: results/de_analysis)'
    )
    parser.add_argument(
        '--condition',
        default='condition',
        help='Metadata column for condition (default: condition)'
    )
    parser.add_argument(
        '--reference',
        help='Reference level for comparison'
    )
    parser.add_argument(
        '--batch',
        help='Metadata column for batch correction'
    )
    parser.add_argument(
        '--fdr',
        type=float,
        default=0.05,
        help='FDR threshold (default: 0.05)'
    )
    parser.add_argument(
        '--lfc',
        type=float,
        default=0.0,
        help='Log2 fold change threshold (default: 0.0)'
    )
    parser.add_argument(
        '--config',
        help='RAPTOR recommendation.yaml file'
    )
    parser.add_argument(
        '--plots',
        action='store_true',
        help='Generate QC plots'
    )
    parser.add_argument(
        '--check-packages',
        action='store_true',
        help='Check R package installation'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(levelname)s: %(message)s'
    )
    
    # Check packages
    if args.check_packages:
        click.echo("Checking R packages...")
        packages = DEAnalysisWrapper.check_r_packages()
        
        all_installed = all(packages.values())
        
        for pkg, installed in packages.items():
            status = "✓" if installed else "✗"
            click.echo(f"  {status} {pkg}")
        
        if not all_installed:
            click.echo("Some packages are missing.")
            click.echo("Run: Rscript install_packages.R")
            return 1
        else:
            click.echo("All packages installed!")
            return 0
    
    # Run analysis
    try:
        de = DEAnalysisWrapper(method=args.method)
        
        result = de.run(
            counts=args.counts,
            metadata=args.metadata,
            output_dir=args.output,
            condition=args.condition,
            reference=args.reference,
            batch=args.batch,
            fdr=args.fdr,
            lfc=args.lfc,
            config=args.config,
            plots=args.plots
        )
        
        click.echo("Analysis complete!")
        click.echo(f"  Method: {result.method.upper()}")
        click.echo(f"  Significant genes: {result.n_significant}")
        click.echo(f"  Upregulated: {result.n_upregulated}")
        click.echo(f"  Downregulated: {result.n_downregulated}")
        click.echo(f"  Results: {result.output_dir}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
