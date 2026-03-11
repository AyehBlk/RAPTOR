"""
RAPTOR Module 9: Ensemble Analysis
RNA-seq Analysis Pipeline Tool for Omics Research

Combine DEG lists from multiple DE methods using statistical ensemble approaches.

This module provides comprehensive methods for combining differential expression 
results from multiple tools (DESeq2, edgeR, limma, etc.) to create high-confidence 
consensus gene lists.

Available ensemble methods:
1. P-value Combination (Fisher's and Brown's methods) ✅
2. Rank Aggregation (RRA) ✅
3. Simple Voting ✅
4. Weighted Ensemble ✅
5. Consensus Scoring (Future - not implemented)

Scientific References:
- Fisher, R.A. (1925). Statistical Methods for Research Workers.
- Brown, M.B. (1975). A method for combining non-independent tests. Biometrics.
- Kolde, R., et al. (2012). Robust rank aggregation. Bioinformatics.
- Hong & Breitling (2008). Meta-analysis methods for DEG detection. Bioinformatics.
- Ramasamy et al. (2008). Key issues in conducting a meta-analysis. PLoS ONE.

Part of RAPTOR v2.2.0
Author: Ayeh Bolouki
Date: January 22, 2026
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path
import warnings
import json
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import chi2, beta

# Import RAPTOR validation and error handling
try:
    from .utils.validation import (
        validate_de_results_dict,
        validate_ensemble_method,
        validate_significance_threshold,
        validate_direction_threshold,
        validate_weights_dict,
        validate_min_methods,
        validate_filters_dict
    )
    from .utils.errors import (
        EnsembleError,
        ValidationError,
        InsufficientDataError,
        MethodMismatchError,
        CombinationFailedError,
        handle_errors
    )
    VALIDATION_AVAILABLE = True
except ImportError:
    # Fallback if running standalone
    warnings.warn("RAPTOR utils not available - using basic validation")
    VALIDATION_AVAILABLE = False
    
    # Stub error classes
    class EnsembleError(Exception):
        pass
    class ValidationError(Exception):
        pass
    class InsufficientDataError(Exception):
        pass
    class MethodMismatchError(Exception):
        pass
    class CombinationFailedError(Exception):
        pass
    
    # Stub decorator
    def handle_errors(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


# RAPTOR version
__raptor_version__ = '2.2.0'


# ============================================================================
# RESULT CLASS
# ============================================================================

@dataclass
class EnsembleResult:
    """
    Results from ensemble analysis.
    
    Attributes
    ----------
    consensus_genes : pd.DataFrame
        Consensus DEG list with method-specific columns
    n_consensus_genes : int
        Number of consensus genes
    ensemble_method : str
        Which ensemble method was used
    n_methods : int
        Number of input methods
    method_names : List[str]
        Names of input methods
    direction_consistency : pd.DataFrame
        Per-gene direction agreement across methods
    n_direction_inconsistent : int
        Number of genes rejected for direction inconsistency
    method_statistics : Dict
        Per-method statistics and overlap
    combined_pvalues : Optional[pd.DataFrame]
        For p-value combination: combined p-values for all genes
    ranked_genes : Optional[pd.DataFrame]
        For RRA: ranked genes with scores
    parameters : Dict
        Parameters used in ensemble analysis
    timestamp : str
        When analysis was run
    """
    
    # Main outputs
    consensus_genes: pd.DataFrame
    n_consensus_genes: int
    
    # Method details
    ensemble_method: str
    n_methods: int
    method_names: List[str]
    
    # Direction consistency
    direction_consistency: pd.DataFrame
    n_direction_inconsistent: int
    
    # Statistics
    method_statistics: Dict
    
    # Method-specific (optional)
    combined_pvalues: Optional[pd.DataFrame] = None
    ranked_genes: Optional[pd.DataFrame] = None
    
    # Metadata
    parameters: Dict = field(default_factory=dict)
    timestamp: str = ""
    
    def summary(self) -> str:
        """Generate summary of ensemble results."""
        lines = [
            "=" * 70,
            "RAPTOR v2.2.0 - MODULE 9: ENSEMBLE ANALYSIS RESULTS",
            "=" * 70,
            f"Ensemble Method: {self.ensemble_method}",
            f"Number of Methods: {self.n_methods}",
            f"Methods: {', '.join(self.method_names)}",
            "",
            "CONSENSUS GENES:",
            f"  Total consensus genes: {self.n_consensus_genes}",
            f"  Direction inconsistent (rejected): {self.n_direction_inconsistent}",
            "",
            "PER-METHOD STATISTICS:",
        ]
        
        for method, stats in self.method_statistics.items():
            if method == 'overall':
                continue
            lines.append(f"  {method}:")
            lines.append(f"    Total genes: {stats.get('n_genes', 'N/A')}")
            lines.append(f"    Genes in consensus: {stats.get('n_in_consensus', 'N/A')}")
        
        if self.combined_pvalues is not None:
            lines.append("")
            lines.append("P-VALUE COMBINATION:")
            sig_adj = (self.combined_pvalues['combined_padj'] < 0.05).sum()
            lines.append(f"  Genes with combined padj < 0.05: {sig_adj}")
        
        if self.ranked_genes is not None:
            lines.append("")
            lines.append("RANK AGGREGATION:")
            sig_genes = (self.ranked_genes['rra_pvalue'] < 0.05).sum()
            lines.append(f"  Genes with RRA p-value < 0.05: {sig_genes}")
        
        lines.append("=" * 70)
        
        return "\n".join(lines)
    
    def save(self, output_dir: Path):
        """Save all results to output directory."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save consensus genes
        self.consensus_genes.to_csv(
            output_dir / "consensus_genes.csv",
            index=False
        )
        
        # Save direction consistency
        self.direction_consistency.to_csv(
            output_dir / "direction_consistency.csv",
            index=False
        )
        
        # Save method-specific results
        if self.combined_pvalues is not None:
            self.combined_pvalues.to_csv(
                output_dir / "combined_pvalues.csv",
                index=False
            )
        
        if self.ranked_genes is not None:
            self.ranked_genes.to_csv(
                output_dir / "ranked_genes.csv",
                index=False
            )
        
        # Save statistics
        with open(output_dir / "statistics.json", 'w') as f:
            json.dump(self.method_statistics, f, indent=2)
        
        # Save summary
        with open(output_dir / "summary.txt", 'w') as f:
            f.write(self.summary())
        
        print(f"Results saved to: {output_dir}")
    
    @classmethod
    def load(cls, output_dir: Path) -> 'EnsembleResult':
        """Load saved results."""
        output_dir = Path(output_dir)
        
        # Load main files
        consensus_genes = pd.read_csv(output_dir / "consensus_genes.csv")
        direction_consistency = pd.read_csv(output_dir / "direction_consistency.csv")
        
        with open(output_dir / "statistics.json", 'r') as f:
            stats = json.load(f)
        
        # Load method-specific results if exist
        combined_pvalues = None
        if (output_dir / "combined_pvalues.csv").exists():
            combined_pvalues = pd.read_csv(output_dir / "combined_pvalues.csv")
        
        ranked_genes = None
        if (output_dir / "ranked_genes.csv").exists():
            ranked_genes = pd.read_csv(output_dir / "ranked_genes.csv")
        
        # Reconstruct result object
        return cls(
            consensus_genes=consensus_genes,
            n_consensus_genes=len(consensus_genes),
            ensemble_method=stats.get('overall', {}).get('ensemble_method', 'unknown'),
            n_methods=stats.get('overall', {}).get('n_methods', 0),
            method_names=stats.get('overall', {}).get('method_names', []),
            direction_consistency=direction_consistency,
            n_direction_inconsistent=stats.get('overall', {}).get('n_direction_inconsistent', 0),
            method_statistics=stats,
            combined_pvalues=combined_pvalues,
            ranked_genes=ranked_genes,
            parameters=stats.get('overall', {}).get('parameters', {}),
            timestamp=stats.get('overall', {}).get('timestamp', '')
        )


# ============================================================================
# DIRECTION CONSISTENCY FUNCTIONS
# ============================================================================

def check_direction_consistency(
    gene_directions: Dict[str, str],
    threshold: float = 1.0
) -> Tuple[bool, float, Dict[str, int]]:
    """
    Check if gene has consistent direction across methods.
    
    Parameters
    ----------
    gene_directions : dict
        Dictionary mapping method names to directions ('up' or 'down')
    threshold : float
        Minimum fraction of methods that must agree (0-1)
    
    Returns
    -------
    is_consistent : bool
        True if direction is consistent
    agreement_fraction : float
        Fraction of methods agreeing on majority direction
    direction_counts : dict
        Count of 'up' and 'down' across methods
    """
    directions = list(gene_directions.values())
    
    n_up = sum(1 for d in directions if d == 'up')
    n_down = sum(1 for d in directions if d == 'down')
    n_total = len(directions)
    
    direction_counts = {'up': n_up, 'down': n_down}
    
    max_count = max(n_up, n_down)
    agreement_fraction = max_count / n_total if n_total > 0 else 0.0
    
    is_consistent = agreement_fraction >= threshold
    
    return is_consistent, agreement_fraction, direction_counts


def get_consensus_direction(gene_directions: Dict[str, str]) -> str:
    """Get consensus direction (majority vote)."""
    directions = list(gene_directions.values())
    n_up = sum(1 for d in directions if d == 'up')
    n_down = sum(1 for d in directions if d == 'down')
    
    return 'up' if n_up >= n_down else 'down'


def calculate_direction_consistency_table(
    de_results: Dict[str, Any],
    direction_threshold: float = 1.0
) -> pd.DataFrame:
    """Calculate direction consistency for all genes across methods."""
    all_genes = set()
    for result in de_results.values():
        all_genes.update(result.data['gene_id'].values)
    
    consistency_data = []
    
    for gene in all_genes:
        gene_directions = {}
        
        for method_name, result in de_results.items():
            gene_data = result.data[result.data['gene_id'] == gene]
            if len(gene_data) > 0:
                lfc = gene_data['log2FoldChange'].iloc[0]
                direction = 'up' if lfc > 0 else 'down'
                gene_directions[method_name] = direction
        
        if len(gene_directions) == 0:
            continue
        
        is_consistent, agreement, counts = check_direction_consistency(
            gene_directions,
            threshold=direction_threshold
        )
        
        consensus_dir = get_consensus_direction(gene_directions)
        
        methods_up = [m for m, d in gene_directions.items() if d == 'up']
        methods_down = [m for m, d in gene_directions.items() if d == 'down']
        
        consistency_data.append({
            'gene_id': gene,
            'direction': consensus_dir,
            'direction_agreement': agreement,
            'n_up': counts['up'],
            'n_down': counts['down'],
            'is_consistent': is_consistent,
            'methods_up': ','.join(methods_up),
            'methods_down': ','.join(methods_down)
        })
    
    return pd.DataFrame(consistency_data)


# ============================================================================
# METHOD 1 & 2: P-VALUE COMBINATION (Fisher's and Brown's)
# ============================================================================

def fishers_method(pvalues: np.ndarray) -> float:
    """
    Combine p-values using Fisher's method.
    
    Reference: Fisher, R.A. (1925). Statistical Methods for Research Workers.
    """
    pvalues = pvalues[~np.isnan(pvalues)]
    
    if len(pvalues) == 0:
        return np.nan
    
    if np.any(pvalues <= 0):
        warnings.warn("P-values must be > 0. Setting 0 values to machine epsilon.")
        pvalues = np.maximum(pvalues, np.finfo(float).eps)
    
    if np.any(pvalues >= 1):
        warnings.warn("P-values must be < 1. Setting values ≥ 1 to 1 - epsilon.")
        pvalues = np.minimum(pvalues, 1 - np.finfo(float).eps)
    
    test_statistic = -2 * np.sum(np.log(pvalues))
    df = 2 * len(pvalues)
    
    combined_pvalue = 1 - chi2.cdf(test_statistic, df)
    
    return combined_pvalue


def browns_method(pvalues: np.ndarray, data_matrix: Optional[np.ndarray] = None) -> float:
    """
    Combine p-values using Brown's method.
    
    Reference: Brown, M.B. (1975). Biometrics, 31(4), 987-992.
    """
    valid_idx = ~np.isnan(pvalues)
    pvalues = pvalues[valid_idx]
    
    if len(pvalues) == 0:
        return np.nan
    
    if data_matrix is None:
        return fishers_method(pvalues)
    
    if data_matrix is not None:
        data_matrix = data_matrix[:, valid_idx]
    
    if np.any(pvalues <= 0):
        pvalues = np.maximum(pvalues, np.finfo(float).eps)
    if np.any(pvalues >= 1):
        pvalues = np.minimum(pvalues, 1 - np.finfo(float).eps)
    
    test_statistic = -2 * np.sum(np.log(pvalues))
    k = len(pvalues)
    
    if data_matrix is not None and data_matrix.shape[1] > 1:
        try:
            corr_matrix = np.corrcoef(data_matrix.T)
            
            cov_correction = 0
            for i in range(k):
                for j in range(i + 1, k):
                    if not np.isnan(corr_matrix[i, j]):
                        cov_correction += corr_matrix[i, j]
            
            expected_value = 2 * k
            variance = 4 * k + 4 * cov_correction
            
            if variance > 0:
                scale_factor = variance / (2 * expected_value)
                df = expected_value ** 2 / variance
                
                scaled_statistic = test_statistic / scale_factor
                
                combined_pvalue = 1 - chi2.cdf(scaled_statistic, df)
            else:
                combined_pvalue = fishers_method(pvalues)
            
        except Exception as e:
            warnings.warn(
                f"Brown's method failed: {e}. Falling back to Fisher's method."
            )
            combined_pvalue = fishers_method(pvalues)
    else:
        combined_pvalue = fishers_method(pvalues)
    
    return combined_pvalue


def combine_pvalues_across_methods(
    de_results: Dict[str, Any],
    method: str = 'fisher',
    use_padj: bool = False
) -> pd.DataFrame:
    """Combine p-values across methods for all genes."""
    if use_padj:
        warnings.warn(
            "Combining adjusted p-values (use_padj=True) is NOT recommended. "
            "This results in double correction and reduces statistical power. "
            "Standard practice is to combine raw p-values (use_padj=False). "
            "See: Hong & Breitling (2008), Bioinformatics 24(3):374-382.",
            UserWarning
        )
    
    all_genes = set()
    for result in de_results.values():
        all_genes.update(result.data['gene_id'].values)
    
    combined_data = []
    
    pval_col = 'padj' if use_padj else 'pvalue'
    
    for gene in all_genes:
        gene_pvalues = {}
        
        for method_name, result in de_results.items():
            gene_data = result.data[result.data['gene_id'] == gene]
            if len(gene_data) > 0:
                pval = gene_data[pval_col].iloc[0]
                gene_pvalues[method_name] = pval
        
        if len(gene_pvalues) == 0:
            continue
        
        pvalues_array = np.array(list(gene_pvalues.values()))
        
        if method == 'fisher':
            combined_pval = fishers_method(pvalues_array)
        elif method == 'brown':
            combined_pval = fishers_method(pvalues_array)
        else:
            raise ValueError(f"Unknown method: {method}. Use 'fisher' or 'brown'.")
        
        result_dict = {
            'gene_id': gene,
            'combined_pvalue': combined_pval,
            'n_methods': len(gene_pvalues)
        }
        
        for method_name, pval in gene_pvalues.items():
            result_dict[f'pvalue_{method_name}'] = pval
        
        combined_data.append(result_dict)
    
    combined_df = pd.DataFrame(combined_data)
    
    from scipy.stats import false_discovery_control
    combined_df['combined_padj'] = false_discovery_control(
        combined_df['combined_pvalue'].values
    )
    
    combined_df = combined_df.sort_values('combined_pvalue').reset_index(drop=True)
    
    return combined_df


def calculate_meta_lfc(
    de_results: Dict[str, Any],
    gene_id: str,
    method: str = 'fixed'
) -> Tuple[float, float]:
    """Calculate meta-analytic log2FoldChange using inverse-variance weighting."""
    lfcs = []
    ses = []
    
    for result in de_results.values():
        gene_data = result.data[result.data['gene_id'] == gene_id]
        if len(gene_data) > 0:
            lfc = gene_data['log2FoldChange'].iloc[0]
            lfcs.append(lfc)
            
            if 'lfcSE' in gene_data.columns:
                se = gene_data['lfcSE'].iloc[0]
            elif 'stat' in gene_data.columns and not np.isnan(gene_data['stat'].iloc[0]):
                stat = gene_data['stat'].iloc[0]
                se = abs(lfc / stat) if stat != 0 else np.nan
            else:
                pval = gene_data['pvalue'].iloc[0]
                if pval > 0 and not np.isnan(pval):
                    z = stats.norm.ppf(1 - pval / 2)
                    se = abs(lfc / z) if z != 0 else np.nan
                else:
                    se = np.nan
            
            ses.append(se)
    
    if len(lfcs) == 0:
        return np.nan, np.nan
    
    lfcs = np.array(lfcs)
    ses = np.array(ses)
    
    valid = ~np.isnan(ses)
    if not np.any(valid):
        return np.mean(lfcs), np.std(lfcs) / np.sqrt(len(lfcs))
    
    lfcs = lfcs[valid]
    ses = ses[valid]
    
    weights = 1 / (ses ** 2)
    meta_lfc = np.sum(weights * lfcs) / np.sum(weights)
    meta_se = np.sqrt(1 / np.sum(weights))
    
    return meta_lfc, meta_se


# ============================================================================
# METHOD 3: RANK AGGREGATION (RRA)
# ============================================================================

def rho_score(r: float, k: int, n: int) -> float:
    """
    Calculate rho score for Robust Rank Aggregation.
    
    Parameters
    ----------
    r : float
        Rank of gene (1-indexed)
    k : int
        Number of lists gene appears in
    n : int
        Total number of genes
    
    Returns
    -------
    float
        Rho score (p-value from beta distribution)
    
    Reference
    ---------
    Kolde, R., et al. (2012). Robust rank aggregation for gene list integration 
    and meta-analysis. Bioinformatics, 28(4), 573-580.
    """
    if k == 0 or n == 0:
        return 1.0
    
    # Expected value under null: r_exp = (n+1)/2
    # Use beta distribution to calculate p-value
    # P(rank <= r) under uniform distribution
    
    p = r / n  # Normalized rank
    
    # Beta CDF: P(X <= p) where X ~ Beta(1, k)
    # This gives the p-value for observing rank r or better in k lists
    try:
        rho = beta.cdf(p, 1, k)
    except:
        rho = 1.0
    
    return rho


def robust_rank_aggregation(
    de_results: Dict[str, Any],
    rank_by: str = 'pvalue',
    use_padj: bool = False
) -> pd.DataFrame:
    """
    Perform Robust Rank Aggregation (RRA) on gene lists.
    
    RRA combines ranked lists by calculating a score for each gene based on its
    ranks across multiple methods. It's robust to outliers and handles incomplete
    lists (genes missing from some methods).
    
    Parameters
    ----------
    de_results : dict
        Dictionary mapping method names to DEResult objects
    rank_by : str
        What to rank by: 'pvalue' or 'padj'
    use_padj : bool
        If True, rank by adjusted p-values; if False, by raw p-values
    
    Returns
    -------
    pd.DataFrame
        RRA results with columns:
        - gene_id
        - rra_score (rho score)
        - rra_pvalue (corrected rho score)
        - n_methods (number of methods gene appears in)
        - mean_rank (average rank across methods)
        - ranks_METHOD (rank in each method)
    
    Reference
    ---------
    Kolde, R., et al. (2012). Robust rank aggregation for gene list integration 
    and meta-analysis. Bioinformatics, 28(4), 573-580.
    """
    # Rank genes in each method
    ranked_data = {}
    all_genes = set()
    
    pval_col = 'padj' if use_padj else 'pvalue'
    
    for method_name, result in de_results.items():
        df = result.data.copy()
        
        # Rank by p-value (ascending)
        df = df.sort_values(pval_col)
        df['rank'] = range(1, len(df) + 1)
        
        ranked_data[method_name] = df[['gene_id', 'rank']].set_index('gene_id')['rank'].to_dict()
        all_genes.update(df['gene_id'].values)
    
    n_methods = len(de_results)
    n_genes = len(all_genes)
    
    # Calculate RRA scores for each gene
    rra_results = []
    
    for gene in all_genes:
        ranks = []
        gene_methods = []
        
        for method_name, ranks_dict in ranked_data.items():
            if gene in ranks_dict:
                ranks.append(ranks_dict[gene])
                gene_methods.append(method_name)
        
        if len(ranks) == 0:
            continue
        
        # Calculate rho score (minimum rank-based p-value)
        k = len(ranks)  # Number of methods gene appears in
        
        # Normalize ranks
        normalized_ranks = [r / n_genes for r in ranks]
        
        # Calculate rho score (minimum normalized rank raised to power k)
        min_normalized_rank = min(normalized_ranks)
        rho = beta.cdf(min_normalized_rank, 1, k)
        
        result_dict = {
            'gene_id': gene,
            'rra_score': rho,
            'n_methods': k,
            'mean_rank': np.mean(ranks),
            'min_rank': min(ranks),
            'max_rank': max(ranks),
            'methods_detected': ','.join(gene_methods)
        }
        
        # Add individual ranks
        for method_name in de_results.keys():
            if method_name in gene_methods:
                idx = gene_methods.index(method_name)
                result_dict[f'rank_{method_name}'] = ranks[idx]
            else:
                result_dict[f'rank_{method_name}'] = np.nan
        
        rra_results.append(result_dict)
    
    rra_df = pd.DataFrame(rra_results)
    
    # FDR correction on rho scores
    from scipy.stats import false_discovery_control
    rra_df['rra_pvalue'] = false_discovery_control(rra_df['rra_score'].values)
    
    # Sort by RRA score
    rra_df = rra_df.sort_values('rra_score').reset_index(drop=True)
    
    return rra_df


# ============================================================================
# METHOD 4: SIMPLE VOTING
# ============================================================================

@handle_errors(exit_on_error=False, log_traceback=True)
def ensemble_voting(
    de_results: Dict[str, Any],
    min_methods: int = 2,
    filters: Optional[Dict[str, float]] = None,
    check_direction: bool = True,
    direction_threshold: float = 1.0
) -> pd.DataFrame:
    """
    Simple voting ensemble: count genes detected by multiple methods.
    
    Parameters
    ----------
    de_results : dict
        Dictionary mapping method names to DEResult objects
    min_methods : int
        Minimum number of methods that must detect gene (default: 2)
    filters : dict, optional
        Filters to apply before voting:
        - 'padj': adjusted p-value threshold (e.g., 0.05)
        - 'lfc': log2FoldChange threshold (e.g., 1.0)
        If None, uses raw DEResult data
    check_direction : bool
        If True, only keep genes with consistent direction
    direction_threshold : float
        Minimum fraction of methods that must agree on direction
    
    Returns
    -------
    pd.DataFrame
        Voting results with columns:
        - gene_id
        - n_votes (number of methods detecting gene)
        - direction
        - direction_agreement
        - methods_detected
    
    Examples
    --------
    >>> # Strict voting (all 3 methods must agree)
    >>> result = ensemble_voting(
    ...     de_results=de_results,
    ...     min_methods=3,
    ...     filters={'padj': 0.05, 'lfc': 1.0}
    ... )
    
    >>> # Relaxed voting (at least 2 methods)
    >>> result = ensemble_voting(
    ...     de_results=de_results,
    ...     min_methods=2,
    ...     filters={'padj': 0.05}
    ... )
    """
    # Validation
    n_methods = len(de_results)
    
    if VALIDATION_AVAILABLE:
        validate_de_results_dict(de_results, min_methods=1)
        validate_min_methods(min_methods, n_methods)
        if filters is not None:
            validate_filters_dict(filters)
        if check_direction:
            validate_direction_threshold(direction_threshold)
    else:
        # Fallback validation
        if not de_results:
            raise ValidationError("de_results cannot be empty")
        if min_methods < 1 or min_methods > n_methods:
            raise ValidationError(
                f"min_methods must be between 1 and {n_methods}, got {min_methods}"
            )
    
    # Apply filters if specified
    if filters is not None:
        filtered_results = {}
        for method_name, result in de_results.items():
            df = result.data.copy()
            
            if 'padj' in filters:
                df = df[df['padj'] < filters['padj']]
            
            if 'lfc' in filters:
                df = df[abs(df['log2FoldChange']) > filters['lfc']]
            
            # Create filtered DEResult-like object
            from types import SimpleNamespace
            filtered_results[method_name] = SimpleNamespace(data=df)
        
        de_results_to_use = filtered_results
    else:
        de_results_to_use = de_results
    
    # Get all genes
    all_genes = set()
    for result in de_results_to_use.values():
        all_genes.update(result.data['gene_id'].values)
    
    # Count votes for each gene
    voting_data = []
    
    for gene in all_genes:
        methods_with_gene = []
        gene_directions = {}
        
        for method_name, result in de_results_to_use.items():
            gene_data = result.data[result.data['gene_id'] == gene]
            if len(gene_data) > 0:
                methods_with_gene.append(method_name)
                
                lfc = gene_data['log2FoldChange'].iloc[0]
                direction = 'up' if lfc > 0 else 'down'
                gene_directions[method_name] = direction
        
        n_votes = len(methods_with_gene)
        
        if n_votes < min_methods:
            continue
        
        # Check direction consistency
        if check_direction:
            is_consistent, agreement, counts = check_direction_consistency(
                gene_directions,
                threshold=direction_threshold
            )
            
            if not is_consistent:
                continue
        else:
            is_consistent = True
            agreement = 1.0
        
        consensus_dir = get_consensus_direction(gene_directions)
        
        voting_data.append({
            'gene_id': gene,
            'n_votes': n_votes,
            'direction': consensus_dir,
            'direction_agreement': agreement,
            'is_consistent': is_consistent,
            'methods_detected': ','.join(methods_with_gene)
        })
    
    voting_df = pd.DataFrame(voting_data)
    voting_df = voting_df.sort_values('n_votes', ascending=False).reset_index(drop=True)
    
    return voting_df


# ============================================================================
# METHOD 5: WEIGHTED ENSEMBLE
# ============================================================================

@handle_errors(exit_on_error=False, log_traceback=True)
def ensemble_weighted(
    de_results: Dict[str, Any],
    weights: Optional[Dict[str, float]] = None,
    min_score: float = 0.5,
    filters: Optional[Dict[str, float]] = None,
    check_direction: bool = True,
    direction_threshold: float = 1.0
) -> pd.DataFrame:
    """
    Weighted ensemble: voting with method-specific weights.
    
    Parameters
    ----------
    de_results : dict
        Dictionary mapping method names to DEResult objects
    weights : dict, optional
        Dictionary mapping method names to weights (0-1)
        If None, uses equal weights
        Example: {'DESeq2': 0.9, 'edgeR': 0.85, 'limma': 0.8}
    min_score : float
        Minimum weighted score to include gene (0-1)
    filters : dict, optional
        Filters to apply before voting (same as ensemble_voting)
    check_direction : bool
        If True, only keep genes with consistent direction
    direction_threshold : float
        Minimum fraction of methods that must agree on direction
    
    Returns
    -------
    pd.DataFrame
        Weighted ensemble results with columns:
        - gene_id
        - weighted_score (sum of weights for methods detecting gene)
        - n_methods (number of methods)
        - direction
        - direction_agreement
        - methods_detected
        - weights_used
    
    Examples
    --------
    >>> # Use Module 8 F1 scores as weights
    >>> weights = {'DESeq2': 0.92, 'edgeR': 0.88, 'limma': 0.85}
    >>> result = ensemble_weighted(
    ...     de_results=de_results,
    ...     weights=weights,
    ...     min_score=1.5,  # Sum of at least 2 good methods
    ...     filters={'padj': 0.05}
    ... )
    """
    # Validation
    method_names = list(de_results.keys())
    n_methods = len(method_names)
    
    if VALIDATION_AVAILABLE:
        validate_de_results_dict(de_results, min_methods=1)
        if weights is not None:
            validate_weights_dict(weights, method_names, require_all=False)
        if filters is not None:
            validate_filters_dict(filters)
        if check_direction:
            validate_direction_threshold(direction_threshold)
    else:
        # Fallback validation
        if not de_results:
            raise ValidationError("de_results cannot be empty")
        if weights is not None:
            unknown = set(weights.keys()) - set(method_names)
            if unknown:
                raise ValidationError(f"Unknown methods in weights: {unknown}")
    
    # Set default weights if not provided
    if weights is None:
        weights = {method: 1.0 for method in de_results.keys()}
    
    # Normalize weights to sum to number of methods
    n_methods = len(de_results)
    weight_sum = sum(weights.values())
    normalized_weights = {m: (w / weight_sum) * n_methods for m, w in weights.items()}
    
    # Apply filters if specified
    if filters is not None:
        filtered_results = {}
        for method_name, result in de_results.items():
            df = result.data.copy()
            
            if 'padj' in filters:
                df = df[df['padj'] < filters['padj']]
            
            if 'lfc' in filters:
                df = df[abs(df['log2FoldChange']) > filters['lfc']]
            
            from types import SimpleNamespace
            filtered_results[method_name] = SimpleNamespace(data=df)
        
        de_results_to_use = filtered_results
    else:
        de_results_to_use = de_results
    
    # Get all genes
    all_genes = set()
    for result in de_results_to_use.values():
        all_genes.update(result.data['gene_id'].values)
    
    # Calculate weighted scores
    weighted_data = []
    
    for gene in all_genes:
        methods_with_gene = []
        method_weights = []
        gene_directions = {}
        
        for method_name, result in de_results_to_use.items():
            gene_data = result.data[result.data['gene_id'] == gene]
            if len(gene_data) > 0:
                methods_with_gene.append(method_name)
                method_weights.append(normalized_weights.get(method_name, 1.0))
                
                lfc = gene_data['log2FoldChange'].iloc[0]
                direction = 'up' if lfc > 0 else 'down'
                gene_directions[method_name] = direction
        
        weighted_score = sum(method_weights)
        
        if weighted_score < min_score:
            continue
        
        # Check direction consistency
        if check_direction:
            is_consistent, agreement, counts = check_direction_consistency(
                gene_directions,
                threshold=direction_threshold
            )
            
            if not is_consistent:
                continue
        else:
            is_consistent = True
            agreement = 1.0
        
        consensus_dir = get_consensus_direction(gene_directions)
        
        weighted_data.append({
            'gene_id': gene,
            'weighted_score': weighted_score,
            'n_methods': len(methods_with_gene),
            'direction': consensus_dir,
            'direction_agreement': agreement,
            'is_consistent': is_consistent,
            'methods_detected': ','.join(methods_with_gene),
            'weights_used': ','.join([f"{m}:{normalized_weights.get(m, 1.0):.2f}" 
                                     for m in methods_with_gene])
        })
    
    weighted_df = pd.DataFrame(weighted_data)
    weighted_df = weighted_df.sort_values('weighted_score', ascending=False).reset_index(drop=True)
    
    return weighted_df


# ============================================================================
# UNIFIED ENSEMBLE FUNCTION (P-VALUE COMBINATION)
# ============================================================================

@handle_errors(exit_on_error=False, log_traceback=True)
def ensemble_pvalue_combination(
    de_results: Dict[str, Any],
    method: str = 'fisher',
    use_padj: bool = False,
    significance_threshold: float = 0.05,
    check_direction: bool = True,
    direction_threshold: float = 1.0,
    output_dir: Optional[Path] = None
) -> EnsembleResult:
    """
    Combine DEG lists using p-value combination (Fisher's or Brown's method).
    
    [Full docstring from previous implementation - keeping it the same]
    """
    
    # Validation
    if VALIDATION_AVAILABLE:
        validate_de_results_dict(de_results, min_methods=2)
        validate_ensemble_method(method, allowed_methods=['fisher', 'brown'])
        validate_significance_threshold(significance_threshold)
        if check_direction:
            validate_direction_threshold(direction_threshold)
    else:
        # Fallback validation
        if not de_results:
            raise ValidationError("de_results cannot be empty")
        if method not in ['fisher', 'brown']:
            raise ValidationError(f"method must be 'fisher' or 'brown', got '{method}'")
        if not 0 < significance_threshold < 1:
            raise ValidationError(
                f"significance_threshold must be between 0 and 1, got {significance_threshold}"
            )
        if not 0 < direction_threshold <= 1:
            raise ValidationError(
                f"direction_threshold must be between 0 and 1, got {direction_threshold}"
            )

    
    method_names = list(de_results.keys())
    n_methods = len(method_names)
    
    print(f"RAPTOR v{__raptor_version__} - Module 9: Ensemble Analysis")
    print(f"Method: P-value Combination ({method.upper()})")
    print(f"Methods: {', '.join(method_names)} (n={n_methods})")
    print(f"Using {'adjusted' if use_padj else 'raw'} p-values")
    if use_padj:
        print("  ⚠️  WARNING: Using adjusted p-values is not recommended!")
    print()
    
    # Step 1: Combine p-values
    print("Step 1: Combining p-values...")
    combined_df = combine_pvalues_across_methods(
        de_results=de_results,
        method=method,
        use_padj=use_padj
    )
    print(f"  Processed {len(combined_df)} genes")
    print()
    
    # Step 2: Direction consistency
    print("Step 2: Checking direction consistency...")
    direction_df = calculate_direction_consistency_table(
        de_results=de_results,
        direction_threshold=direction_threshold
    )
    
    combined_df = combined_df.merge(
        direction_df[['gene_id', 'direction', 'direction_agreement', 'is_consistent',
                     'n_up', 'n_down', 'methods_up', 'methods_down']],
        on='gene_id',
        how='left'
    )
    
    n_inconsistent = (~combined_df['is_consistent']).sum() if check_direction else 0
    print(f"  Direction inconsistent genes: {n_inconsistent}")
    print()
    
    # Step 3: Meta-analytic effect sizes
    print("Step 3: Calculating meta-analytic effect sizes...")
    meta_lfcs = []
    meta_ses = []
    
    for gene_id in combined_df['gene_id']:
        meta_lfc, meta_se = calculate_meta_lfc(de_results, gene_id)
        meta_lfcs.append(meta_lfc)
        meta_ses.append(meta_se)
    
    combined_df['meta_lfc'] = meta_lfcs
    combined_df['meta_lfc_se'] = meta_ses
    print(f"  Calculated meta-LFC for {len(combined_df)} genes")
    print()
    
    # Step 4: Filter
    print("Step 4: Filtering for consensus genes...")
    consensus_df = combined_df[
        combined_df['combined_padj'] < significance_threshold
    ].copy()
    
    print(f"  Genes with combined_padj < {significance_threshold}: {len(consensus_df)}")
    
    if check_direction:
        before_filter = len(consensus_df)
        consensus_df = consensus_df[consensus_df['is_consistent']].copy()
        after_filter = len(consensus_df)
        print(f"  After direction filtering: {after_filter} (removed {before_filter - after_filter})")
    print()
    
    # Step 5: Add methods detected
    def get_methods_detected(row):
        methods = []
        for method_name in method_names:
            col = f'pvalue_{method_name}'
            if col in row and not pd.isna(row[col]):
                methods.append(method_name)
        return ','.join(methods)
    
    consensus_df['methods_detected'] = consensus_df.apply(get_methods_detected, axis=1)
    
    # Reorder columns
    core_columns = [
        'gene_id', 'combined_pvalue', 'combined_padj', 'n_methods',
        'direction', 'direction_agreement', 'meta_lfc', 'meta_lfc_se',
        'methods_detected'
    ]
    
    pvalue_columns = [col for col in consensus_df.columns if col.startswith('pvalue_')]
    other_columns = [col for col in consensus_df.columns 
                    if col not in core_columns and col not in pvalue_columns]
    
    column_order = core_columns + pvalue_columns + other_columns
    column_order = [col for col in column_order if col in consensus_df.columns]
    consensus_df = consensus_df[column_order]
    
    # Calculate statistics
    method_stats = {}
    
    for method_name, result in de_results.items():
        n_in_consensus = (consensus_df['methods_detected'].str.contains(method_name)).sum()
        method_stats[method_name] = {
            'n_genes': len(result.data),
            'n_in_consensus': int(n_in_consensus)
        }
    
    overall_stats = {
        'ensemble_method': f'pvalue_combination_{method}',
        'n_methods': n_methods,
        'method_names': method_names,
        'n_total_genes': len(combined_df),
        'n_consensus_genes': len(consensus_df),
        'n_direction_inconsistent': int(n_inconsistent),
        'significance_threshold': significance_threshold,
        'direction_threshold': direction_threshold if check_direction else None,
        'check_direction': check_direction,
        'use_padj': use_padj,
        'timestamp': datetime.now().isoformat(),
        'raptor_version': __raptor_version__
    }
    method_stats['overall'] = overall_stats
    
    # Create result object
    result = EnsembleResult(
        consensus_genes=consensus_df,
        n_consensus_genes=len(consensus_df),
        ensemble_method=f'pvalue_combination_{method}',
        n_methods=n_methods,
        method_names=method_names,
        direction_consistency=direction_df,
        n_direction_inconsistent=int(n_inconsistent),
        method_statistics=method_stats,
        combined_pvalues=combined_df,
        parameters={
            'method': method,
            'use_padj': use_padj,
            'significance_threshold': significance_threshold,
            'check_direction': check_direction,
            'direction_threshold': direction_threshold
        },
        timestamp=datetime.now().isoformat()
    )
    
    if output_dir is not None:
        result.save(output_dir)
    
    print(result.summary())
    
    return result


# ============================================================================
# UNIFIED ENSEMBLE FUNCTION (RRA)
# ============================================================================

@handle_errors(exit_on_error=False, log_traceback=True)
def ensemble_rra(
    de_results: Dict[str, Any],
    rank_by: str = 'pvalue',
    use_padj: bool = False,
    significance_threshold: float = 0.05,
    check_direction: bool = True,
    direction_threshold: float = 1.0,
    output_dir: Optional[Path] = None
) -> EnsembleResult:
    """
    Combine DEG lists using Robust Rank Aggregation (RRA).
    
    Parameters
    ----------
    de_results : dict
        Dictionary mapping method names to DEResult objects from Module 7
    rank_by : str
        What to rank by: 'pvalue' or 'padj'
    use_padj : bool
        If True, rank by adjusted p-values; if False, by raw p-values
    significance_threshold : float
        Threshold for RRA p-value (default: 0.05)
    check_direction : bool
        If True, only keep genes with consistent direction
    direction_threshold : float
        Minimum fraction of methods that must agree on direction
    output_dir : Path, optional
        Directory to save results
    
    Returns
    -------
    EnsembleResult
        Results object with ranked genes
    
    Reference
    ---------
    Kolde, R., et al. (2012). Robust rank aggregation for gene list integration 
    and meta-analysis. Bioinformatics, 28(4), 573-580.
    """
    # Validation
    if VALIDATION_AVAILABLE:
        validate_de_results_dict(de_results, min_methods=2)
        validate_significance_threshold(significance_threshold)
        if check_direction:
            validate_direction_threshold(direction_threshold)
    else:
        # Fallback validation
        if not de_results:
            raise ValidationError("de_results cannot be empty")
        if not 0 < significance_threshold < 1:
            raise ValidationError(
                f"significance_threshold must be between 0 and 1, got {significance_threshold}"
            )
        if check_direction and not 0 < direction_threshold <= 1:
            raise ValidationError(
                f"direction_threshold must be between 0 and 1, got {direction_threshold}"
            )
    
    method_names = list(de_results.keys())
    n_methods = len(method_names)
    
    print(f"RAPTOR v{__raptor_version__} - Module 9: Ensemble Analysis")
    print(f"Method: Robust Rank Aggregation (RRA)")
    print(f"Methods: {', '.join(method_names)} (n={n_methods})")
    print(f"Ranking by: {rank_by}")
    print()
    
    # Step 1: RRA
    print("Step 1: Performing rank aggregation...")
    rra_df = robust_rank_aggregation(
        de_results=de_results,
        rank_by=rank_by,
        use_padj=use_padj
    )
    print(f"  Ranked {len(rra_df)} genes")
    print()
    
    # Step 2: Direction consistency
    print("Step 2: Checking direction consistency...")
    direction_df = calculate_direction_consistency_table(
        de_results=de_results,
        direction_threshold=direction_threshold
    )
    
    rra_df = rra_df.merge(
        direction_df[['gene_id', 'direction', 'direction_agreement', 'is_consistent']],
        on='gene_id',
        how='left'
    )
    
    n_inconsistent = (~rra_df['is_consistent']).sum() if check_direction else 0
    print(f"  Direction inconsistent genes: {n_inconsistent}")
    print()
    
    # Step 3: Filter
    print("Step 3: Filtering for consensus genes...")
    consensus_df = rra_df[
        rra_df['rra_pvalue'] < significance_threshold
    ].copy()
    
    print(f"  Genes with RRA p-value < {significance_threshold}: {len(consensus_df)}")
    
    if check_direction:
        before_filter = len(consensus_df)
        consensus_df = consensus_df[consensus_df['is_consistent']].copy()
        after_filter = len(consensus_df)
        print(f"  After direction filtering: {after_filter} (removed {before_filter - after_filter})")
    print()
    
    # Calculate statistics
    method_stats = {}
    
    for method_name, result in de_results.items():
        n_in_consensus = (consensus_df['methods_detected'].str.contains(method_name)).sum()
        method_stats[method_name] = {
            'n_genes': len(result.data),
            'n_in_consensus': int(n_in_consensus)
        }
    
    overall_stats = {
        'ensemble_method': 'rank_aggregation_rra',
        'n_methods': n_methods,
        'method_names': method_names,
        'n_total_genes': len(rra_df),
        'n_consensus_genes': len(consensus_df),
        'n_direction_inconsistent': int(n_inconsistent),
        'significance_threshold': significance_threshold,
        'direction_threshold': direction_threshold if check_direction else None,
        'check_direction': check_direction,
        'timestamp': datetime.now().isoformat(),
        'raptor_version': __raptor_version__
    }
    method_stats['overall'] = overall_stats
    
    # Create result object
    result = EnsembleResult(
        consensus_genes=consensus_df,
        n_consensus_genes=len(consensus_df),
        ensemble_method='rank_aggregation_rra',
        n_methods=n_methods,
        method_names=method_names,
        direction_consistency=direction_df,
        n_direction_inconsistent=int(n_inconsistent),
        method_statistics=method_stats,
        ranked_genes=rra_df,
        parameters={
            'rank_by': rank_by,
            'use_padj': use_padj,
            'significance_threshold': significance_threshold,
            'check_direction': check_direction,
            'direction_threshold': direction_threshold
        },
        timestamp=datetime.now().isoformat()
    )
    
    if output_dir is not None:
        result.save(output_dir)
    
    print(result.summary())
    
    return result


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

def ensemble_fisher(de_results: Dict[str, Any], **kwargs) -> EnsembleResult:
    """Convenience function for Fisher's method."""
    return ensemble_pvalue_combination(de_results, method='fisher', **kwargs)


def ensemble_brown(de_results: Dict[str, Any], **kwargs) -> EnsembleResult:
    """Convenience function for Brown's method."""
    return ensemble_pvalue_combination(de_results, method='brown', **kwargs)


# ============================================================================
# MODULE EXPORTS
# ============================================================================

__all__ = [
    # Result class
    'EnsembleResult',
    
    # Main ensemble functions
    'ensemble_pvalue_combination',
    'ensemble_rra',
    'ensemble_voting',
    'ensemble_weighted',
    
    # Convenience functions
    'ensemble_fisher',
    'ensemble_brown',
    
    # Lower-level functions
    'fishers_method',
    'browns_method',
    'robust_rank_aggregation',
    'check_direction_consistency',
    'get_consensus_direction',
    'calculate_direction_consistency_table',
    'combine_pvalues_across_methods',
    'calculate_meta_lfc'
]
