#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Module 7: Import DE Results

Imports and standardizes differential expression results from external R analysis.
Compatible with DESeq2, edgeR, limma-voom, and Wilcoxon test outputs.

Workflow Position:
    M6: External R Analysis → de_results.csv
    M7: Import DE Results (THIS MODULE) → DEResult object
    M8: Parameter Optimization
    M9: Ensemble Analysis
    M10: Biomarker Discovery

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
License: MIT
"""

import logging
import json
import pickle
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple, Union
from dataclasses import dataclass, field
import warnings

import numpy as np
import pandas as pd

# RAPTOR imports
try:
    from raptor.utils.validation import (
        validate_file_exists,
        validate_file_extension,
        validate_dataframe_columns,
        validate_value_range
    )
    from raptor.utils.errors import ValidationError, PipelineError, handle_errors
except ImportError:
    warnings.warn("RAPTOR utils not found. Using fallback validation.")
    
    # Fallback implementations
    class ValidationError(Exception):
        """Validation error."""
        def __init__(self, parameter: str, message: str = "", hint: str = ""):
            self.parameter = parameter
            self.message = message
            self.hint = hint
            super().__init__(f"{parameter}: {message}. {hint}")
    
    class PipelineError(Exception):
        """Pipeline error."""
        pass
    
    def handle_errors(func):
        """Simple error handling decorator."""
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logging.error(f"Error in {func.__name__}: {e}")
                raise
        return wrapper
    
    def validate_file_exists(filepath: Path, hint: str = "") -> None:
        if not filepath.exists():
            raise ValidationError('file', f"File not found: {filepath}", hint)
    
    def validate_file_extension(filepath: Path, extensions: List[str], hint: str = "") -> None:
        if filepath.suffix.lower() not in extensions:
            raise ValidationError('file', f"Invalid extension: {filepath.suffix}", hint)
    
    def validate_dataframe_columns(df: pd.DataFrame, required: List[str], hint: str = "") -> None:
        missing = set(required) - set(df.columns)
        if missing:
            raise ValidationError('columns', f"Missing columns: {missing}", hint)
    
    def validate_value_range(value: float, min_val: float, max_val: float, name: str, hint: str = "") -> None:
        if not (min_val <= value <= max_val):
            raise ValidationError(name, f"Value {value} outside range [{min_val}, {max_val}]", hint)


logger = logging.getLogger(__name__)


# =============================================================================
# Constants
# =============================================================================

# Column mappings for each R package
COLUMN_MAPPINGS = {
    'deseq2': {
        'log2_fold_change': ['log2FoldChange', 'log2FC', 'lfc', 'logFC'],
        'p_value': ['pvalue', 'pval', 'PValue', 'p.value'],
        'adjusted_p_value': ['padj', 'FDR', 'adj.P.Val', 'qvalue', 'q.value'],
        'base_mean': ['baseMean', 'AveExpr', 'logCPM', 'meanExpr'],
        'statistic': ['stat', 't', 'LR', 'F', 'z'],
        'se': ['lfcSE', 'SE', 'se']
    },
    'edger': {
        'log2_fold_change': ['logFC', 'log2FoldChange', 'lfc'],
        'p_value': ['PValue', 'pvalue', 'pval', 'p.value'],
        'adjusted_p_value': ['FDR', 'padj', 'adj.P.Val', 'qvalue'],
        'base_mean': ['logCPM', 'AveExpr', 'baseMean'],
        'statistic': ['LR', 'F', 'stat', 't']
    },
    'limma': {
        'log2_fold_change': ['logFC', 'log2FoldChange', 'lfc'],
        'p_value': ['P.Value', 'pvalue', 'PValue', 'pval'],
        'adjusted_p_value': ['adj.P.Val', 'padj', 'FDR', 'qvalue'],
        'base_mean': ['AveExpr', 'logCPM', 'baseMean'],
        'statistic': ['t', 'B', 'stat', 'F'],
        'b_statistic': ['B']
    },
    'wilcoxon': {
        'log2_fold_change': ['logFC', 'log2FoldChange', 'lfc'],
        'p_value': ['pvalue', 'PValue', 'pval', 'p.value'],
        'adjusted_p_value': ['padj', 'FDR', 'adj.P.Val', 'qvalue'],
        'base_mean': ['AveExpr', 'logCPM', 'baseMean'],
        'statistic': ['W', 'stat']
    }
}

# Required standardized columns
REQUIRED_COLUMNS = [
    'gene_id',
    'log2_fold_change',
    'p_value',
    'adjusted_p_value'
]

# CLI Parameters for RAPTOR integration
IMPORT_DE_CLI_PARAMS = {
    'de_file': {
        'type': str,
        'required': True,
        'help': 'Path to DE results CSV file from R analysis'
    },
    'output_dir': {
        'type': str,
        'default': 'results/de_imported',
        'help': 'Output directory for imported results'
    },
    'pipeline': {
        'type': str,
        'choices': ['auto', 'deseq2', 'edger', 'limma', 'wilcoxon'],
        'default': 'auto',
        'help': 'R package used (auto-detect if not specified)'
    },
    'fdr_threshold': {
        'type': float,
        'default': 0.05,
        'help': 'FDR threshold for significance (0-1)'
    },
    'lfc_threshold': {
        'type': float,
        'default': 0.0,
        'help': 'Log2 fold change threshold for significance'
    },
    'gene_id_column': {
        'type': str,
        'default': None,
        'help': 'Column name containing gene IDs (auto-detect if not specified)'
    }
}


# =============================================================================
# DEResult Data Class
# =============================================================================

@dataclass
class DEResult:
    """
    Standardized differential expression result object.
    
    Core data structure for RAPTOR Modules 7-10 workflow.
    Provides unified interface for parameter optimization, ensemble analysis,
    and biomarker discovery.
    
    Attributes
    ----------
    results_df : pd.DataFrame
        Standardized DE results with index as gene_id and columns:
        - log2_fold_change : float
        - p_value : float
        - adjusted_p_value : float
        - base_mean : float (optional)
        - statistic : float (optional)
        - is_significant : bool
        - direction : str ('up', 'down', 'unchanged')
    
    pipeline : str
        Name of DE pipeline used ('DESeq2', 'edgeR', 'limma-voom', 'Wilcoxon')
    
    parameters : dict
        Analysis parameters including:
        - fdr_threshold : float
        - lfc_threshold : float
        - pipeline_specific_params : dict
    
    metadata : dict
        Additional metadata including:
        - source_file : str
        - timestamp : str
        - raptor_version : str
        - n_samples : int (if available)
    
    Examples
    --------
    >>> from raptor.de_import import import_de_results
    >>> de_result = import_de_results('results/de_deseq2/de_results.csv')
    >>> print(f"Significant genes: {de_result.n_significant}")
    >>> top_genes = de_result.get_top_genes(n=50)
    >>> filtered = de_result.filter_by_threshold(fdr=0.01, lfc=1.0)
    """
    results_df: pd.DataFrame
    pipeline: str
    parameters: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate DEResult after initialization."""
        # Ensure index is gene_id
        if self.results_df.index.name != 'gene_id':
            if 'gene_id' in self.results_df.columns:
                self.results_df = self.results_df.set_index('gene_id')
            else:
                self.results_df.index.name = 'gene_id'
        
        # Ensure required columns exist
        required = ['log2_fold_change', 'p_value', 'adjusted_p_value']
        missing = set(required) - set(self.results_df.columns)
        if missing:
            raise ValidationError(
                parameter='columns',
                message=f"Missing required columns: {missing}",
                hint="DEResult requires log2_fold_change, p_value, and adjusted_p_value"
            )
    
    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------
    
    @property
    def n_genes(self) -> int:
        """Total number of genes tested."""
        return len(self.results_df)
    
    @property
    def n_significant(self) -> int:
        """Number of significant DE genes."""
        if 'is_significant' in self.results_df.columns:
            return int(self.results_df['is_significant'].sum())
        return 0
    
    @property
    def n_up(self) -> int:
        """Number of upregulated genes."""
        if 'direction' in self.results_df.columns and 'is_significant' in self.results_df.columns:
            return int(((self.results_df['direction'] == 'up') & 
                       self.results_df['is_significant']).sum())
        return 0
    
    @property
    def n_down(self) -> int:
        """Number of downregulated genes."""
        if 'direction' in self.results_df.columns and 'is_significant' in self.results_df.columns:
            return int(((self.results_df['direction'] == 'down') & 
                       self.results_df['is_significant']).sum())
        return 0
    
    @property
    def significant_genes(self) -> List[str]:
        """List of significant gene IDs."""
        if 'is_significant' in self.results_df.columns:
            return self.results_df[self.results_df['is_significant']].index.tolist()
        return []
    
    @property
    def upregulated_genes(self) -> List[str]:
        """List of upregulated gene IDs."""
        if 'direction' in self.results_df.columns and 'is_significant' in self.results_df.columns:
            mask = ((self.results_df['direction'] == 'up') & 
                   self.results_df['is_significant'])
            return self.results_df[mask].index.tolist()
        return []
    
    @property
    def downregulated_genes(self) -> List[str]:
        """List of downregulated gene IDs."""
        if 'direction' in self.results_df.columns and 'is_significant' in self.results_df.columns:
            mask = ((self.results_df['direction'] == 'down') & 
                   self.results_df['is_significant'])
            return self.results_df[mask].index.tolist()
        return []
    
    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------
    
    def get_top_genes(
        self,
        n: int = 50,
        by: str = 'adjusted_p_value',
        significant_only: bool = False
    ) -> pd.DataFrame:
        """
        Get top N genes by specified metric.
        
        Parameters
        ----------
        n : int
            Number of genes to return
        by : str
            Metric to sort by: 'adjusted_p_value', 'p_value', 'log2_fold_change', 'base_mean'
        significant_only : bool
            Only return significant genes
        
        Returns
        -------
        pd.DataFrame
            Top N genes
        """
        df = self.results_df
        
        if significant_only and 'is_significant' in df.columns:
            df = df[df['is_significant']]
        
        if by == 'adjusted_p_value':
            return df.nsmallest(n, 'adjusted_p_value')
        elif by == 'p_value':
            return df.nsmallest(n, 'p_value')
        elif by == 'log2_fold_change':
            return df.reindex(df['log2_fold_change'].abs().nlargest(n).index)
        elif by == 'base_mean' and 'base_mean' in df.columns:
            return df.nlargest(n, 'base_mean')
        else:
            return df.head(n)
    
    def filter_by_threshold(
        self,
        fdr_threshold: float = 0.05,
        lfc_threshold: float = 0.0
    ) -> 'DEResult':
        """
        Create new DEResult with different significance thresholds.
        
        Parameters
        ----------
        fdr_threshold : float
            FDR threshold (0-1)
        lfc_threshold : float
            Log2 fold change threshold (absolute value)
        
        Returns
        -------
        DEResult
            New DEResult object with updated significance calls
        """
        # Validate thresholds
        validate_value_range(fdr_threshold, 0, 1, 'fdr_threshold',
                           "FDR threshold must be between 0 and 1")
        validate_value_range(lfc_threshold, 0, 10, 'lfc_threshold',
                           "LFC threshold typically between 0 and 2")
        
        # Calculate new significance
        mask = (
            (self.results_df['adjusted_p_value'] < fdr_threshold) &
            (self.results_df['log2_fold_change'].abs() > lfc_threshold)
        )
        
        new_df = self.results_df.copy()
        new_df['is_significant'] = mask
        new_df['direction'] = np.where(
            new_df['log2_fold_change'] > 0, 'up',
            np.where(new_df['log2_fold_change'] < 0, 'down', 'unchanged')
        )
        new_df.loc[~mask, 'direction'] = 'unchanged'
        
        # Update parameters
        new_params = self.parameters.copy()
        new_params['fdr_threshold'] = fdr_threshold
        new_params['lfc_threshold'] = lfc_threshold
        
        return DEResult(
            results_df=new_df,
            pipeline=self.pipeline,
            parameters=new_params,
            metadata=self.metadata
        )
    
    def calculate_metrics(
        self,
        ground_truth: Optional[pd.DataFrame] = None
    ) -> Dict[str, float]:
        """
        Calculate performance metrics.
        
        Parameters
        ----------
        ground_truth : pd.DataFrame, optional
            Ground truth DE status (from simulation).
            Must have index as gene_id and column 'is_de' or 'DE' as boolean.
        
        Returns
        -------
        dict
            Performance metrics including:
            - Basic: n_genes, n_significant, n_up, n_down
            - With ground truth: sensitivity, specificity, precision, f1_score, actual_fdr
        """
        metrics = {
            'n_genes': self.n_genes,
            'n_significant': self.n_significant,
            'n_up': self.n_up,
            'n_down': self.n_down,
            'fdr_threshold': self.parameters.get('fdr_threshold', 0.05),
            'lfc_threshold': self.parameters.get('lfc_threshold', 0.0)
        }
        
        if ground_truth is not None:
            # Align gene IDs
            common_genes = self.results_df.index.intersection(ground_truth.index)
            
            if len(common_genes) == 0:
                logger.warning("No common genes between results and ground truth")
                return metrics
            
            predicted_de = set(self.significant_genes) & set(common_genes)
            
            # Get true DE genes
            if 'is_de' in ground_truth.columns:
                true_de = set(ground_truth[ground_truth['is_de']].index) & set(common_genes)
            elif 'DE' in ground_truth.columns:
                true_de = set(ground_truth[ground_truth['DE']].index) & set(common_genes)
            else:
                logger.warning("Ground truth must have 'is_de' or 'DE' column")
                return metrics
            
            true_non_de = set(common_genes) - true_de
            
            # Confusion matrix
            tp = len(predicted_de & true_de)
            fp = len(predicted_de & true_non_de)
            fn = len(true_de - predicted_de)
            tn = len(true_non_de - predicted_de)
            
            metrics.update({
                'true_positives': tp,
                'false_positives': fp,
                'false_negatives': fn,
                'true_negatives': tn,
                'sensitivity': tp / (tp + fn) if (tp + fn) > 0 else 0.0,
                'specificity': tn / (tn + fp) if (tn + fp) > 0 else 0.0,
                'precision': tp / (tp + fp) if (tp + fp) > 0 else 0.0,
                'actual_fdr': fp / (fp + tp) if (fp + tp) > 0 else 0.0,
            })
            
            # F1 score
            if metrics['precision'] + metrics['sensitivity'] > 0:
                metrics['f1_score'] = (
                    2 * metrics['precision'] * metrics['sensitivity'] /
                    (metrics['precision'] + metrics['sensitivity'])
                )
            else:
                metrics['f1_score'] = 0.0
            
            # Accuracy
            total = tp + tn + fp + fn
            metrics['accuracy'] = (tp + tn) / total if total > 0 else 0.0
        
        return metrics
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "",
            "╔" + "═" * 60 + "╗",
            "║" + "  🦖 RAPTOR DE Results Summary".center(60) + "║",
            "╠" + "═" * 60 + "╣",
            "",
            f"  Pipeline: {self.pipeline}",
            f"  Total genes tested: {self.n_genes:,}",
            "",
            f"  Significant genes: {self.n_significant:,} ({100*self.n_significant/self.n_genes:.1f}%)",
            f"    ↑ Upregulated:   {self.n_up:,}",
            f"    ↓ Downregulated: {self.n_down:,}",
            "",
            "  Thresholds:",
            f"    FDR: {self.parameters.get('fdr_threshold', 0.05)}",
            f"    LFC: {self.parameters.get('lfc_threshold', 0.0)}",
            "",
            "╚" + "═" * 60 + "╝",
            ""
        ]
        return "\n".join(lines)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'results_df': self.results_df,
            'pipeline': self.pipeline,
            'parameters': self.parameters,
            'metadata': self.metadata
        }
    
    def save(self, path: Union[str, Path]) -> None:
        """
        Save DEResult object.
        
        Parameters
        ----------
        path : str or Path
            Output path (.pkl, .pickle, or .json)
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        if path.suffix in ['.pkl', '.pickle']:
            with open(path, 'wb') as f:
                pickle.dump(self, f)
        elif path.suffix == '.json':
            # Save as JSON (DataFrame as records)
            data = {
                'results_df': self.results_df.reset_index().to_dict('records'),
                'pipeline': self.pipeline,
                'parameters': self.parameters,
                'metadata': self.metadata
            }
            with open(path, 'w') as f:
                json.dump(data, f, indent=2, default=str)
        else:
            raise ValidationError(
                parameter='path',
                message=f"Unsupported file extension: {path.suffix}",
                hint="Use .pkl, .pickle, or .json"
            )
        
        logger.info(f"DEResult saved to {path}")
    
    @classmethod
    def load(cls, path: Union[str, Path]) -> 'DEResult':
        """
        Load DEResult object.
        
        Parameters
        ----------
        path : str or Path
            Path to DEResult file (.pkl or .json)
        
        Returns
        -------
        DEResult
            Loaded DEResult object
        """
        path = Path(path)
        validate_file_exists(path, f"DEResult file not found: {path}")
        
        if path.suffix in ['.pkl', '.pickle']:
            with open(path, 'rb') as f:
                return pickle.load(f)
        elif path.suffix == '.json':
            with open(path, 'r') as f:
                data = json.load(f)
            
            # Reconstruct DataFrame
            df = pd.DataFrame(data['results_df']).set_index('gene_id')
            
            return cls(
                results_df=df,
                pipeline=data['pipeline'],
                parameters=data['parameters'],
                metadata=data['metadata']
            )
        else:
            raise ValidationError(
                parameter='path',
                message=f"Unsupported file extension: {path.suffix}",
                hint="Use .pkl, .pickle, or .json"
            )


# =============================================================================
# DE Import Functions
# =============================================================================

def detect_pipeline(df: pd.DataFrame) -> str:
    """
    Auto-detect which R package generated the results.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results DataFrame
    
    Returns
    -------
    str
        Detected pipeline: 'deseq2', 'edger', 'limma', or 'wilcoxon'
    """
    columns = set(df.columns)
    
    # DESeq2 signature
    if 'baseMean' in columns and 'lfcSE' in columns:
        return 'deseq2'
    
    # edgeR signature
    if 'logCPM' in columns and ('LR' in columns or 'F' in columns):
        return 'edger'
    
    # limma signature
    if 'AveExpr' in columns and 'B' in columns:
        return 'limma'
    
    # Wilcoxon or generic
    if 'logFC' in columns and 'pvalue' in columns:
        return 'wilcoxon'
    
    # Default
    logger.warning("Could not auto-detect pipeline. Assuming generic format.")
    return 'deseq2'  # Use DESeq2 mappings as default


def standardize_columns(
    df: pd.DataFrame,
    pipeline: str,
    gene_id_column: Optional[str] = None
) -> pd.DataFrame:
    """
    Standardize column names based on pipeline.
    
    Parameters
    ----------
    df : pd.DataFrame
        Raw DE results
    pipeline : str
        Pipeline name
    gene_id_column : str, optional
        Column containing gene IDs
    
    Returns
    -------
    pd.DataFrame
        Standardized DataFrame
    """
    df = df.copy()
    
    # Find gene ID column
    if gene_id_column is None:
        possible_gene_cols = ['gene_id', 'GeneID', 'gene', 'Gene', 'ID', 'id']
        for col in possible_gene_cols:
            if col in df.columns:
                gene_id_column = col
                break
        
        if gene_id_column is None and df.index.name not in [None, '']:
            # Use index as gene_id
            df['gene_id'] = df.index
            gene_id_column = 'gene_id'
    
    if gene_id_column is None:
        logger.warning("Could not find gene ID column. Using row numbers.")
        df['gene_id'] = [f"Gene_{i}" for i in range(len(df))]
        gene_id_column = 'gene_id'
    
    # Set gene_id as index
    if gene_id_column != 'gene_id':
        df['gene_id'] = df[gene_id_column]
    df = df.set_index('gene_id', drop=False)
    
    # Get column mappings
    mappings = COLUMN_MAPPINGS.get(pipeline, COLUMN_MAPPINGS['deseq2'])
    
    # Standardize columns
    standardized = {'gene_id': df['gene_id']}
    
    for standard_name, possible_names in mappings.items():
        for col_name in possible_names:
            if col_name in df.columns:
                standardized[standard_name] = df[col_name]
                break
    
    result_df = pd.DataFrame(standardized, index=df.index)
    
    # Verify required columns
    required = ['log2_fold_change', 'p_value', 'adjusted_p_value']
    missing = set(required) - set(result_df.columns)
    
    if missing:
        raise ValidationError(
            parameter='columns',
            message=f"Could not find required columns: {missing}",
            hint=f"Available columns: {list(df.columns)}"
        )
    
    return result_df


def calculate_significance(
    df: pd.DataFrame,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.0
) -> pd.DataFrame:
    """
    Calculate significance and direction.
    
    Parameters
    ----------
    df : pd.DataFrame
        Standardized results
    fdr_threshold : float
        FDR threshold
    lfc_threshold : float
        Log2FC threshold
    
    Returns
    -------
    pd.DataFrame
        Results with is_significant and direction columns
    """
    df = df.copy()
    
    # Calculate significance
    df['is_significant'] = (
        (df['adjusted_p_value'] < fdr_threshold) &
        (df['log2_fold_change'].abs() > lfc_threshold)
    )
    
    # Calculate direction
    df['direction'] = np.where(
        df['log2_fold_change'] > 0, 'up',
        np.where(df['log2_fold_change'] < 0, 'down', 'unchanged')
    )
    df.loc[~df['is_significant'], 'direction'] = 'unchanged'
    
    return df


@handle_errors
def import_de_results(
    de_file: Union[str, Path],
    output_dir: Union[str, Path] = 'results/de_imported',
    pipeline: str = 'auto',
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 0.0,
    gene_id_column: Optional[str] = None
) -> DEResult:
    """
    Import and standardize DE results from R analysis.
    
    Parameters
    ----------
    de_file : str or Path
        Path to DE results CSV file
    output_dir : str or Path
        Output directory
    pipeline : str
        Pipeline name or 'auto' for auto-detection
    fdr_threshold : float
        FDR threshold for significance
    lfc_threshold : float
        Log2FC threshold for significance
    gene_id_column : str, optional
        Column containing gene IDs
    
    Returns
    -------
    DEResult
        Standardized DE result object
    """
    de_file = Path(de_file)
    output_dir = Path(output_dir)
    
    logger.info("╔══════════════════════════════════════════════════════════════╗")
    logger.info("║       🦖 RAPTOR v2.2.0 - Import DE Results (Module 7)        ║")
    logger.info("╚══════════════════════════════════════════════════════════════╝")
    logger.info("")
    
    # Validation
    validate_file_exists(de_file, "DE results file not found")
    validate_file_extension(de_file, ['.csv', '.tsv', '.txt'],
                          "DE results must be CSV or TSV format")
    validate_value_range(fdr_threshold, 0, 1, 'fdr_threshold',
                       "FDR threshold must be between 0 and 1")
    
    # Load data
    logger.info(f"📂 Loading DE results from: {de_file}")
    if de_file.suffix == '.tsv':
        df = pd.read_csv(de_file, sep='\t')
    else:
        df = pd.read_csv(de_file)
    
    logger.info(f"   Loaded {len(df):,} genes")
    
    # Detect pipeline
    if pipeline == 'auto':
        pipeline = detect_pipeline(df)
        logger.info(f"   Auto-detected pipeline: {pipeline}")
    else:
        logger.info(f"   Using specified pipeline: {pipeline}")
    
    # Standardize columns
    logger.info("🔄 Standardizing column names...")
    df_standard = standardize_columns(df, pipeline, gene_id_column)
    
    # Calculate significance
    logger.info(f"📊 Calculating significance (FDR={fdr_threshold}, LFC={lfc_threshold})...")
    df_standard = calculate_significance(df_standard, fdr_threshold, lfc_threshold)
    
    # Create DEResult object
    parameters = {
        'fdr_threshold': fdr_threshold,
        'lfc_threshold': lfc_threshold,
        'source_pipeline': pipeline
    }
    
    metadata = {
        'source_file': str(de_file),
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.2.0',
        'module': 'M7'
    }
    
    de_result = DEResult(
        results_df=df_standard,
        pipeline=pipeline.upper(),
        parameters=parameters,
        metadata=metadata
    )
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save outputs
    logger.info(f"\n📁 Saving results to: {output_dir}/")
    
    # 1. Standardized results
    output_file = output_dir / 'de_standardized.csv'
    df_standard.to_csv(output_file)
    logger.info(f"   ✓ de_standardized.csv ({len(df_standard):,} genes)")
    
    # 2. Significant genes only
    sig_df = df_standard[df_standard['is_significant']]
    sig_file = output_dir / 'de_significant.csv'
    sig_df.to_csv(sig_file)
    logger.info(f"   ✓ de_significant.csv ({len(sig_df):,} genes)")
    
    # 3. Summary JSON
    metrics = de_result.calculate_metrics()
    summary_file = output_dir / 'de_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(metrics, f, indent=2, default=str)
    logger.info(f"   ✓ de_summary.json")
    
    # 4. DEResult object (for M8-M10)
    result_file = output_dir / 'de_result.pkl'
    de_result.save(result_file)
    logger.info(f"   ✓ de_result.pkl (for M8-M10)")
    
    # Print summary
    logger.info(de_result.summary())
    
    logger.info("✅ Import complete!")
    logger.info(f"\n  Next steps:")
    logger.info(f"    M8: raptor optimize --de-result {result_file}")
    logger.info(f"    M9: raptor ensemble --de-results [multiple_files]")
    logger.info(f"    M10: raptor biomarkers --de-result {result_file}")
    logger.info("")
    
    return de_result


# =============================================================================
# CLI Integration
# =============================================================================

def main():
    """CLI entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Module 7 - Import DE Results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    for param, config in IMPORT_DE_CLI_PARAMS.items():
        args = {'dest': param}
        if 'type' in config:
            args['type'] = config['type']
        if 'default' in config:
            args['default'] = config['default']
        if 'required' in config:
            args['required'] = config['required']
        if 'choices' in config:
            args['choices'] = config['choices']
        if 'help' in config:
            args['help'] = config['help']
        
        param_name = '--' + param.replace('_', '-')
        parser.add_argument(param_name, **args)
    
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(message)s'
    )
    
    # Import results
    try:
        de_result = import_de_results(
            de_file=args.de_file,
            output_dir=args.output_dir,
            pipeline=args.pipeline,
            fdr_threshold=args.fdr_threshold,
            lfc_threshold=args.lfc_threshold,
            gene_id_column=args.gene_id_column
        )
        return 0
    except Exception as e:
        logger.error(f"❌ Error: {e}")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
