"""
Adaptive Threshold Optimizer (ATO) for RNA-seq Differential Expression Analysis

This module provides data-driven methods for optimizing significance thresholds
in differential expression analysis, replacing arbitrary cutoffs with 
scientifically justified values.

Author: Ayeh Bolouki
License: MIT
Version: 2.1.1

References:
- Benjamini Y, Hochberg Y (1995) JRSS-B 57:289-300
- Benjamini Y, Yekutieli D (2001) Ann Stat 29:1165-1188
- Storey JD (2002) JRSS-B 64:479-498
- McCarthy DJ, Smyth GK (2009) Bioinformatics 25:765-771
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from typing import Dict, List, Tuple, Optional, Union, Any
from dataclasses import dataclass, field
import warnings


@dataclass
class ThresholdResult:
    """Container for threshold optimization results."""
    
    # Recommended thresholds
    logfc_cutoff: float
    padj_cutoff: float
    padj_method: str
    
    # Method details
    logfc_method: str
    logfc_reasoning: str
    padj_reasoning: str
    
    # Statistics
    n_significant_optimized: int
    n_significant_traditional: int
    pi0_estimate: float
    
    # Additional metrics
    metrics: Dict[str, Any] = field(default_factory=dict)
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 60,
            "ADAPTIVE THRESHOLD OPTIMIZER - RESULTS",
            "=" * 60,
            "",
            "RECOMMENDED THRESHOLDS:",
            f"  • Log2 Fold Change cutoff: |logFC| > {self.logfc_cutoff:.3f}",
            f"  • Adjusted p-value cutoff: padj < {self.padj_cutoff:.4f}",
            f"  • P-value adjustment method: {self.padj_method}",
            "",
            "REASONING:",
            f"  • LogFC method: {self.logfc_method}",
            f"    {self.logfc_reasoning}",
            "",
            f"  • Padj method: {self.padj_method}",
            f"    {self.padj_reasoning}",
            "",
            "IMPACT:",
            f"  • DE genes (optimized): {self.n_significant_optimized}",
            f"  • DE genes (traditional |logFC|>1, padj<0.05): {self.n_significant_traditional}",
            f"  • Estimated π₀ (proportion of true nulls): {self.pi0_estimate:.3f}",
            "=" * 60
        ]
        return "\n".join(lines)


class AdaptiveThresholdOptimizer:
    """
    Adaptive Threshold Optimizer for Differential Expression Analysis.
    
    Provides data-driven methods for selecting optimal significance thresholds
    based on the statistical properties of your specific dataset.
    
    Parameters
    ----------
    df : pd.DataFrame
        Differential expression results with columns:
        - 'log2FoldChange' or 'logFC': log2 fold changes
        - 'pvalue' or 'PValue': raw p-values  
        - 'padj' or 'adj.P.Val' (optional): adjusted p-values
        - 'baseMean' or 'AveExpr' (optional): mean expression
        - 'lfcSE' (optional): logFC standard error
        
    goal : str, default='discovery'
        Analysis goal: 'discovery', 'validation', or 'balanced'
        - 'discovery': Maximize true positives (FDR control)
        - 'validation': Minimize false positives (FWER control)
        - 'balanced': Balance sensitivity and specificity
        
    Examples
    --------
    >>> import pandas as pd
    >>> from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
    >>> 
    >>> # Load DESeq2 results
    >>> df = pd.read_csv('deseq2_results.csv')
    >>> 
    >>> # Optimize thresholds
    >>> ato = AdaptiveThresholdOptimizer(df, goal='discovery')
    >>> result = ato.optimize()
    >>> print(result.summary())
    >>> 
    >>> # Apply optimized thresholds
    >>> significant = ato.get_significant_genes()
    """
    
    # Column name mappings for different DE tools
    COLUMN_MAPPINGS = {
        'logfc': ['log2FoldChange', 'logFC', 'log2FC', 'lfc', 'LogFC'],
        'pvalue': ['pvalue', 'PValue', 'P.Value', 'pval', 'Pvalue'],
        'padj': ['padj', 'adj.P.Val', 'FDR', 'qvalue', 'q.value', 'BH', 'fdr'],
        'basemean': ['baseMean', 'AveExpr', 'logCPM', 'meanExpr', 'aveexpr'],
        'lfcse': ['lfcSE', 'SE', 'se', 'logFC.SE']
    }
    
    def __init__(
        self,
        df: pd.DataFrame,
        goal: str = 'discovery',
        verbose: bool = True
    ):
        self.original_df = df.copy()
        self.goal = goal.lower()
        self.verbose = verbose
        
        # Validate goal
        if self.goal not in ['discovery', 'validation', 'balanced']:
            raise ValueError("goal must be 'discovery', 'validation', or 'balanced'")
        
        # Standardize column names
        self.df = self._standardize_columns(df)
        
        # Remove rows with NaN in critical columns
        self._clean_data()
        
        # Initialize results storage
        self._results: Optional[ThresholdResult] = None
        self._metrics: Dict[str, Any] = {}
        
        if self.verbose:
            print(f"Loaded {len(self.df)} genes for threshold optimization")
            print(f"Analysis goal: {self.goal}")
    
    def _standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize column names from various DE tools."""
        df = df.copy()
        
        for standard_name, variants in self.COLUMN_MAPPINGS.items():
            for variant in variants:
                if variant in df.columns:
                    df = df.rename(columns={variant: standard_name})
                    break
        
        # Check required columns
        required = ['logfc', 'pvalue']
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}. "
                           f"Available columns: {list(df.columns)}")
        
        return df
    
    def _clean_data(self):
        """Remove rows with NaN values in critical columns."""
        initial_count = len(self.df)
        
        # Remove NaN in logfc and pvalue
        self.df = self.df.dropna(subset=['logfc', 'pvalue'])
        
        # Remove pvalue = 0 (causes issues with -log10)
        self.df = self.df[self.df['pvalue'] > 0]
        
        removed = initial_count - len(self.df)
        if removed > 0 and self.verbose:
            print(f"Removed {removed} genes with missing/invalid values")
    
    # =========================================================================
    # P-VALUE ADJUSTMENT METHODS
    # =========================================================================
    
    def benjamini_hochberg(self, pvalues: np.ndarray) -> np.ndarray:
        """
        Benjamini-Hochberg FDR correction.
        
        Controls FDR under independence or positive regression dependence (PRDS).
        
        Reference: Benjamini & Hochberg (1995)
        """
        n = len(pvalues)
        sorted_idx = np.argsort(pvalues)
        sorted_pvals = pvalues[sorted_idx]
        
        # BH adjustment: p_adj[i] = p[i] * n / rank[i]
        ranks = np.arange(1, n + 1)
        adjusted = sorted_pvals * n / ranks
        
        # Ensure monotonicity (cumulative minimum from right)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0, 1)
        
        # Restore original order
        result = np.empty(n)
        result[sorted_idx] = adjusted
        return result
    
    def benjamini_yekutieli(self, pvalues: np.ndarray) -> np.ndarray:
        """
        Benjamini-Yekutieli FDR correction.
        
        Controls FDR under ANY dependence structure.
        More conservative than BH.
        
        Reference: Benjamini & Yekutieli (2001)
        """
        n = len(pvalues)
        
        # Correction factor: c(n) = sum(1/i) for i=1 to n
        c_n = np.sum(1 / np.arange(1, n + 1))
        
        # BY is BH multiplied by c(n)
        bh_adjusted = self.benjamini_hochberg(pvalues)
        by_adjusted = np.clip(bh_adjusted * c_n, 0, 1)
        
        return by_adjusted
    
    def holm(self, pvalues: np.ndarray) -> np.ndarray:
        """
        Holm step-down procedure for FWER control.
        
        More powerful than Bonferroni, valid under any dependence.
        
        Reference: Holm (1979)
        """
        n = len(pvalues)
        sorted_idx = np.argsort(pvalues)
        sorted_pvals = pvalues[sorted_idx]
        
        # Holm: p_adj[i] = p[i] * (n - rank + 1)
        multipliers = n - np.arange(n)
        adjusted = sorted_pvals * multipliers
        
        # Ensure monotonicity (cumulative maximum from left)
        adjusted = np.maximum.accumulate(adjusted)
        adjusted = np.clip(adjusted, 0, 1)
        
        result = np.empty(n)
        result[sorted_idx] = adjusted
        return result
    
    def hochberg(self, pvalues: np.ndarray) -> np.ndarray:
        """
        Hochberg step-up procedure for FWER control.
        
        More powerful than Holm, requires PRDS.
        
        Reference: Hochberg (1988)
        """
        n = len(pvalues)
        sorted_idx = np.argsort(pvalues)[::-1]  # Descending
        sorted_pvals = pvalues[sorted_idx]
        
        # Hochberg: p_adj[i] = p[i] * (n - rank + 1), step-up
        multipliers = np.arange(1, n + 1)
        adjusted = sorted_pvals * multipliers
        
        # Ensure monotonicity (cumulative minimum from left)
        adjusted = np.minimum.accumulate(adjusted)
        adjusted = np.clip(adjusted, 0, 1)
        
        result = np.empty(n)
        result[sorted_idx] = adjusted
        return result
    
    def bonferroni(self, pvalues: np.ndarray) -> np.ndarray:
        """
        Bonferroni correction for FWER control.
        
        Most conservative, valid under any dependence.
        """
        return np.clip(pvalues * len(pvalues), 0, 1)
    
    def storey_qvalue(self, pvalues: np.ndarray, lambda_range: np.ndarray = None) -> Tuple[np.ndarray, float]:
        """
        Storey's q-value with π₀ estimation.
        
        More powerful than BH when many genes are truly DE.
        
        Parameters
        ----------
        pvalues : array
            Raw p-values
        lambda_range : array, optional
            Range of λ values for π₀ estimation
            
        Returns
        -------
        qvalues : array
            Estimated q-values (FDR)
        pi0 : float
            Estimated proportion of true nulls
            
        Reference: Storey (2002), Storey & Tibshirani (2003)
        """
        pvalues = np.asarray(pvalues)
        n = len(pvalues)
        
        if lambda_range is None:
            lambda_range = np.arange(0.05, 0.96, 0.05)
        
        # Estimate π₀ for each λ
        pi0_estimates = []
        for lam in lambda_range:
            # π₀(λ) = #{p > λ} / (n * (1 - λ))
            pi0_lam = np.sum(pvalues > lam) / (n * (1 - lam))
            pi0_estimates.append(min(pi0_lam, 1.0))
        
        pi0_estimates = np.array(pi0_estimates)
        
        # Smooth π₀ estimates using natural cubic spline
        # Simplified: use weighted average biased toward higher λ
        weights = lambda_range ** 2
        pi0 = np.average(pi0_estimates, weights=weights)
        pi0 = min(max(pi0, 0.01), 1.0)  # Bound between 0.01 and 1
        
        # Calculate q-values
        sorted_idx = np.argsort(pvalues)
        sorted_pvals = pvalues[sorted_idx]
        
        # q[i] = π₀ * n * p[i] / rank[i]
        ranks = np.arange(1, n + 1)
        qvalues = pi0 * n * sorted_pvals / ranks
        
        # Ensure monotonicity
        qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
        qvalues = np.clip(qvalues, 0, 1)
        
        result = np.empty(n)
        result[sorted_idx] = qvalues
        
        return result, pi0
    
    # =========================================================================
    # π₀ ESTIMATION
    # =========================================================================
    
    def estimate_pi0(self, method: str = 'storey') -> float:
        """
        Estimate π₀ (proportion of true null hypotheses).
        
        Parameters
        ----------
        method : str
            'storey': Storey's spline method
            'pounds': Pounds & Cheng method
            'histogram': Histogram-based estimation
            
        Returns
        -------
        pi0 : float
            Estimated proportion of true nulls (0-1)
        """
        pvalues = self.df['pvalue'].values
        
        if method == 'storey':
            _, pi0 = self.storey_qvalue(pvalues)
            
        elif method == 'pounds':
            # Pounds & Cheng (2006) method
            # Simple: 2 * mean(p) is estimator of π₀
            pi0 = min(2 * np.mean(pvalues), 1.0)
            
        elif method == 'histogram':
            # Histogram method: proportion in highest bin
            hist, edges = np.histogram(pvalues, bins=20)
            # Average of right half of histogram (p > 0.5)
            right_half = hist[10:]
            pi0 = 2 * np.mean(right_half) / len(pvalues) * 10
            pi0 = min(max(pi0, 0.01), 1.0)
            
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return pi0
    
    # =========================================================================
    # LOGFC CUTOFF OPTIMIZATION
    # =========================================================================
    
    def optimize_logfc_mad(self, k: float = 2.5) -> Tuple[float, str]:
        """
        MAD-based (Median Absolute Deviation) logFC cutoff.
        
        Uses non-significant genes to estimate null distribution,
        then sets cutoff at k * robust_sigma.
        
        Parameters
        ----------
        k : float
            Number of standard deviations (default 2.5)
            
        Returns
        -------
        cutoff : float
            Recommended |logFC| cutoff
        reasoning : str
            Explanation of the cutoff
        """
        logfc = self.df['logfc'].values
        
        # Use genes with high p-values as null set
        if 'padj' in self.df.columns:
            null_mask = self.df['padj'] > 0.5
        else:
            null_mask = self.df['pvalue'] > 0.5
        
        if null_mask.sum() < 100:
            # Not enough null genes, use all
            null_logfc = logfc
            warning_msg = " (using all genes due to limited null set)"
        else:
            null_logfc = logfc[null_mask]
            warning_msg = ""
        
        # Calculate MAD
        median_logfc = np.median(null_logfc)
        mad = np.median(np.abs(null_logfc - median_logfc))
        
        # Robust sigma (MAD * 1.4826 ≈ σ for normal distribution)
        robust_sigma = 1.4826 * mad
        
        # Cutoff
        cutoff = k * robust_sigma
        
        reasoning = (f"Based on null gene distribution: MAD={mad:.3f}, "
                    f"robust σ={robust_sigma:.3f}. Cutoff = {k}σ{warning_msg}")
        
        self._metrics['mad'] = mad
        self._metrics['robust_sigma'] = robust_sigma
        
        return cutoff, reasoning
    
    def optimize_logfc_mixture(self) -> Tuple[float, str]:
        """
        Mixture model-based logFC cutoff using EM algorithm.
        
        Fits a 2-component Gaussian mixture model:
        - Component 1: Null genes (centered near 0)
        - Component 2: DE genes (spread away from 0)
        
        Returns
        -------
        cutoff : float
            Recommended |logFC| cutoff
        reasoning : str
            Explanation of the cutoff
        """
        logfc = self.df['logfc'].values
        abs_logfc = np.abs(logfc)
        
        # Simple 2-component mixture model using EM
        # Initialize
        mu1, mu2 = 0.1, 0.5
        sigma1, sigma2 = 0.2, 0.5
        pi1 = 0.7  # Prior for null component
        
        # EM iterations
        for _ in range(50):
            # E-step: calculate responsibilities
            pdf1 = stats.norm.pdf(abs_logfc, mu1, sigma1)
            pdf2 = stats.norm.pdf(abs_logfc, mu2, sigma2)
            
            denom = pi1 * pdf1 + (1 - pi1) * pdf2 + 1e-10
            gamma1 = pi1 * pdf1 / denom
            gamma2 = (1 - pi1) * pdf2 / denom
            
            # M-step: update parameters
            n1 = gamma1.sum()
            n2 = gamma2.sum()
            
            pi1 = n1 / len(abs_logfc)
            
            mu1 = np.sum(gamma1 * abs_logfc) / (n1 + 1e-10)
            mu2 = np.sum(gamma2 * abs_logfc) / (n2 + 1e-10)
            
            sigma1 = np.sqrt(np.sum(gamma1 * (abs_logfc - mu1)**2) / (n1 + 1e-10))
            sigma2 = np.sqrt(np.sum(gamma2 * (abs_logfc - mu2)**2) / (n2 + 1e-10))
            
            # Ensure sigma > 0
            sigma1 = max(sigma1, 0.05)
            sigma2 = max(sigma2, 0.05)
        
        # Find crossover point where P(null) = 0.05
        def null_prob(x):
            pdf1 = stats.norm.pdf(x, mu1, sigma1)
            pdf2 = stats.norm.pdf(x, mu2, sigma2)
            return pi1 * pdf1 / (pi1 * pdf1 + (1 - pi1) * pdf2 + 1e-10)
        
        # Search for cutoff where null probability < 0.05
        test_points = np.linspace(0, 3, 300)
        null_probs = [null_prob(x) for x in test_points]
        
        try:
            cutoff_idx = np.where(np.array(null_probs) < 0.05)[0][0]
            cutoff = test_points[cutoff_idx]
        except IndexError:
            # Fallback to mu1 + 2*sigma1
            cutoff = mu1 + 2 * sigma1
        
        reasoning = (f"Mixture model: Null component (μ={mu1:.3f}, σ={sigma1:.3f}), "
                    f"DE component (μ={mu2:.3f}, σ={sigma2:.3f}). "
                    f"Cutoff where P(null)<0.05")
        
        self._metrics['mixture'] = {
            'mu1': mu1, 'sigma1': sigma1,
            'mu2': mu2, 'sigma2': sigma2,
            'pi_null': pi1
        }
        
        return cutoff, reasoning
    
    def optimize_logfc_power(
        self,
        n1: int = 3,
        n2: int = 3,
        power: float = 0.8,
        alpha: float = 0.05
    ) -> Tuple[float, str]:
        """
        Power-based minimum detectable effect size.
        
        Calculates the minimum logFC that can be reliably detected
        given sample size and desired power.
        
        Parameters
        ----------
        n1, n2 : int
            Sample sizes per group
        power : float
            Desired statistical power (default 0.8)
        alpha : float
            Significance level (default 0.05)
            
        Returns
        -------
        cutoff : float
            Minimum detectable |logFC|
        reasoning : str
            Explanation
        """
        # Estimate biological coefficient of variation from data
        if 'lfcse' in self.df.columns:
            # Use provided standard errors
            median_se = np.median(self.df['lfcse'])
        else:
            # Estimate from logFC distribution
            # SE ≈ |logFC| / |z-score| for significant genes
            if 'stat' in self.df.columns:
                se_estimates = np.abs(self.df['logfc'] / (self.df['stat'] + 1e-10))
                median_se = np.median(se_estimates[np.isfinite(se_estimates)])
            else:
                # Fallback: use MAD of logFC
                mad = np.median(np.abs(self.df['logfc'] - np.median(self.df['logfc'])))
                median_se = 1.4826 * mad / np.sqrt(n1 + n2)
        
        # z-scores for power calculation
        z_alpha = stats.norm.ppf(1 - alpha/2)  # Two-sided
        z_beta = stats.norm.ppf(power)
        
        # Minimum detectable effect
        # min_logfc = (z_alpha + z_beta) * SE
        min_logfc = (z_alpha + z_beta) * median_se
        
        reasoning = (f"Power analysis: n={n1}+{n2}, power={power}, α={alpha}. "
                    f"Median SE={median_se:.3f}. "
                    f"Minimum reliably detectable |logFC|")
        
        self._metrics['power'] = {
            'median_se': median_se,
            'z_alpha': z_alpha,
            'z_beta': z_beta,
            'n1': n1,
            'n2': n2
        }
        
        return min_logfc, reasoning
    
    def optimize_logfc_percentile(self, percentile: float = 95) -> Tuple[float, str]:
        """
        Percentile-based cutoff using null distribution.
        
        Parameters
        ----------
        percentile : float
            Percentile of null distribution (default 95)
            
        Returns
        -------
        cutoff : float
            Recommended |logFC| cutoff
        reasoning : str
            Explanation
        """
        logfc = self.df['logfc'].values
        
        # Use genes with high p-values as null
        if 'padj' in self.df.columns:
            null_mask = self.df['padj'] > 0.5
        else:
            null_mask = self.df['pvalue'] > 0.5
        
        if null_mask.sum() < 100:
            null_logfc = logfc
        else:
            null_logfc = logfc[null_mask]
        
        cutoff = np.percentile(np.abs(null_logfc), percentile)
        
        reasoning = f"{percentile}th percentile of null gene |logFC| distribution"
        
        return cutoff, reasoning
    
    # =========================================================================
    # PADJ METHOD SELECTION
    # =========================================================================
    
    def recommend_padj_method(self) -> Tuple[str, str]:
        """
        Recommend optimal p-value adjustment method based on data characteristics.
        
        Returns
        -------
        method : str
            Recommended method name
        reasoning : str
            Explanation
        """
        pvalues = self.df['pvalue'].values
        n = len(pvalues)
        
        # Estimate π₀
        pi0 = self.estimate_pi0(method='storey')
        self._metrics['pi0'] = pi0
        
        # For validation goal, use FWER control
        if self.goal == 'validation':
            return 'holm', f"FWER control for validation (Holm step-down, π₀={pi0:.3f})"
        
        # Check for potential correlation issues
        # Heuristic: if π₀ is very low, there's strong signal
        # If π₀ is near 1, most genes are null
        
        many_de_genes = pi0 < 0.8
        
        if self.goal == 'balanced':
            if many_de_genes:
                return 'qvalue', f"q-value recommended: strong signal detected (π₀={pi0:.3f})"
            else:
                return 'BH', f"BH recommended: moderate signal (π₀={pi0:.3f})"
        
        # Discovery goal
        if many_de_genes:
            method = 'qvalue'
            reasoning = (f"q-value recommended: π₀={pi0:.3f} indicates many DE genes. "
                        f"q-value provides ~{(1-pi0)*100:.0f}% more power than BH")
        else:
            method = 'BH'
            reasoning = f"BH recommended: π₀={pi0:.3f}, standard FDR control sufficient"
        
        return method, reasoning
    
    def recommend_padj_cutoff(self, method: str = 'BH') -> Tuple[float, str]:
        """
        Recommend adjusted p-value cutoff.
        
        For most analyses, 0.05 is appropriate. For more stringent
        control or very large datasets, may recommend 0.01 or 0.1.
        
        Returns
        -------
        cutoff : float
            Recommended padj cutoff
        reasoning : str
            Explanation
        """
        n = len(self.df)
        pi0 = self._metrics.get('pi0', self.estimate_pi0())
        
        if self.goal == 'validation':
            cutoff = 0.01
            reasoning = "Stringent cutoff for validation (0.01)"
        elif self.goal == 'discovery':
            if pi0 < 0.5:
                cutoff = 0.1
                reasoning = f"Relaxed cutoff (0.1) appropriate: strong signal (π₀={pi0:.3f})"
            else:
                cutoff = 0.05
                reasoning = "Standard FDR cutoff (0.05)"
        else:  # balanced
            cutoff = 0.05
            reasoning = "Standard FDR cutoff (0.05) for balanced analysis"
        
        return cutoff, reasoning
    
    # =========================================================================
    # MAIN OPTIMIZATION METHOD
    # =========================================================================
    
    def optimize(
        self,
        logfc_method: str = 'auto',
        n1: int = 3,
        n2: int = 3
    ) -> ThresholdResult:
        """
        Run full threshold optimization.
        
        Parameters
        ----------
        logfc_method : str
            Method for logFC optimization:
            - 'auto': Automatically select best method
            - 'mad': MAD-based robust estimation
            - 'mixture': Gaussian mixture model
            - 'power': Power-based minimum effect
            - 'percentile': Percentile-based
            
        n1, n2 : int
            Sample sizes (used for power analysis)
            
        Returns
        -------
        ThresholdResult
            Object containing optimized thresholds and reasoning
        """
        if self.verbose:
            print("\n" + "="*50)
            print("Running Adaptive Threshold Optimization")
            print("="*50)
        
        # 1. Optimize p-value adjustment method
        padj_method, padj_method_reasoning = self.recommend_padj_method()
        padj_cutoff, padj_cutoff_reasoning = self.recommend_padj_cutoff(padj_method)
        
        if self.verbose:
            print(f"\n[1/3] P-value adjustment: {padj_method}")
            print(f"      {padj_method_reasoning}")
        
        # 2. Calculate adjusted p-values with recommended method
        pvalues = self.df['pvalue'].values
        
        if padj_method == 'BH':
            self.df['padj_optimized'] = self.benjamini_hochberg(pvalues)
        elif padj_method == 'BY':
            self.df['padj_optimized'] = self.benjamini_yekutieli(pvalues)
        elif padj_method == 'holm':
            self.df['padj_optimized'] = self.holm(pvalues)
        elif padj_method == 'hochberg':
            self.df['padj_optimized'] = self.hochberg(pvalues)
        elif padj_method == 'bonferroni':
            self.df['padj_optimized'] = self.bonferroni(pvalues)
        elif padj_method == 'qvalue':
            self.df['padj_optimized'], _ = self.storey_qvalue(pvalues)
        else:
            # Default to BH
            self.df['padj_optimized'] = self.benjamini_hochberg(pvalues)
        
        # 3. Optimize logFC cutoff
        if logfc_method == 'auto':
            # Run multiple methods and pick consensus
            cutoffs = {}
            cutoffs['mad'], reason_mad = self.optimize_logfc_mad()
            cutoffs['mixture'], reason_mix = self.optimize_logfc_mixture()
            cutoffs['power'], reason_pow = self.optimize_logfc_power(n1, n2)
            cutoffs['percentile'], reason_pct = self.optimize_logfc_percentile()
            
            # Use median of methods as consensus
            logfc_cutoff = np.median(list(cutoffs.values()))
            logfc_method_used = 'consensus'
            logfc_reasoning = (f"Consensus of methods: MAD={cutoffs['mad']:.3f}, "
                             f"Mixture={cutoffs['mixture']:.3f}, "
                             f"Power={cutoffs['power']:.3f}, "
                             f"Percentile={cutoffs['percentile']:.3f}")
            
            self._metrics['logfc_methods'] = cutoffs
            
        elif logfc_method == 'mad':
            logfc_cutoff, logfc_reasoning = self.optimize_logfc_mad()
            logfc_method_used = 'MAD-based'
        elif logfc_method == 'mixture':
            logfc_cutoff, logfc_reasoning = self.optimize_logfc_mixture()
            logfc_method_used = 'Mixture model'
        elif logfc_method == 'power':
            logfc_cutoff, logfc_reasoning = self.optimize_logfc_power(n1, n2)
            logfc_method_used = 'Power-based'
        elif logfc_method == 'percentile':
            logfc_cutoff, logfc_reasoning = self.optimize_logfc_percentile()
            logfc_method_used = 'Percentile-based'
        else:
            raise ValueError(f"Unknown logfc_method: {logfc_method}")
        
        if self.verbose:
            print(f"\n[2/3] LogFC cutoff: |logFC| > {logfc_cutoff:.3f}")
            print(f"      Method: {logfc_method_used}")
            print(f"      {logfc_reasoning}")
        
        # 4. Count significant genes
        # Optimized thresholds
        sig_optimized = (
            (np.abs(self.df['logfc']) > logfc_cutoff) & 
            (self.df['padj_optimized'] < padj_cutoff)
        ).sum()
        
        # Traditional thresholds (|logFC| > 1, padj < 0.05)
        if 'padj' in self.df.columns:
            sig_traditional = (
                (np.abs(self.df['logfc']) > 1.0) & 
                (self.df['padj'] < 0.05)
            ).sum()
        else:
            sig_traditional = (
                (np.abs(self.df['logfc']) > 1.0) & 
                (self.df['padj_optimized'] < 0.05)
            ).sum()
        
        if self.verbose:
            print(f"\n[3/3] Results comparison:")
            print(f"      Optimized: {sig_optimized} DE genes")
            print(f"      Traditional: {sig_traditional} DE genes")
        
        # Create result object
        self._results = ThresholdResult(
            logfc_cutoff=logfc_cutoff,
            padj_cutoff=padj_cutoff,
            padj_method=padj_method,
            logfc_method=logfc_method_used,
            logfc_reasoning=logfc_reasoning,
            padj_reasoning=f"{padj_method_reasoning}. {padj_cutoff_reasoning}",
            n_significant_optimized=sig_optimized,
            n_significant_traditional=sig_traditional,
            pi0_estimate=self._metrics.get('pi0', 0),
            metrics=self._metrics.copy()
        )
        
        return self._results
    
    # =========================================================================
    # HELPER METHODS
    # =========================================================================
    
    def get_significant_genes(
        self,
        logfc_cutoff: float = None,
        padj_cutoff: float = None,
        use_optimized: bool = True
    ) -> pd.DataFrame:
        """
        Get significant genes based on thresholds.
        
        Parameters
        ----------
        logfc_cutoff : float, optional
            LogFC cutoff (uses optimized if None)
        padj_cutoff : float, optional
            Padj cutoff (uses optimized if None)
        use_optimized : bool
            Use optimized padj values (True) or original (False)
            
        Returns
        -------
        pd.DataFrame
            Subset of significant genes
        """
        if logfc_cutoff is None:
            if self._results is None:
                raise ValueError("Run optimize() first or provide cutoffs")
            logfc_cutoff = self._results.logfc_cutoff
        
        if padj_cutoff is None:
            if self._results is None:
                raise ValueError("Run optimize() first or provide cutoffs")
            padj_cutoff = self._results.padj_cutoff
        
        padj_col = 'padj_optimized' if use_optimized and 'padj_optimized' in self.df.columns else 'padj'
        
        if padj_col not in self.df.columns:
            padj_col = 'padj_optimized'
        
        mask = (np.abs(self.df['logfc']) > logfc_cutoff) & (self.df[padj_col] < padj_cutoff)
        
        return self.df[mask].copy()
    
    def compare_thresholds(
        self,
        logfc_values: List[float] = [0.5, 1.0, 1.5, 2.0],
        padj_values: List[float] = [0.01, 0.05, 0.1]
    ) -> pd.DataFrame:
        """
        Compare number of DE genes across different threshold combinations.
        
        Returns
        -------
        pd.DataFrame
            Comparison table with gene counts
        """
        padj_col = 'padj_optimized' if 'padj_optimized' in self.df.columns else 'padj'
        
        results = []
        for logfc in logfc_values:
            for padj in padj_values:
                n_sig = (
                    (np.abs(self.df['logfc']) > logfc) & 
                    (self.df[padj_col] < padj)
                ).sum()
                results.append({
                    'logFC_cutoff': logfc,
                    'padj_cutoff': padj,
                    'n_significant': n_sig
                })
        
        return pd.DataFrame(results)
    
    def get_adjustment_comparison(self) -> pd.DataFrame:
        """
        Compare different p-value adjustment methods.
        
        Returns
        -------
        pd.DataFrame
            Comparison of adjustment methods
        """
        pvalues = self.df['pvalue'].values
        
        results = {
            'gene': self.df.index,
            'pvalue': pvalues,
            'BH': self.benjamini_hochberg(pvalues),
            'BY': self.benjamini_yekutieli(pvalues),
            'Holm': self.holm(pvalues),
            'Hochberg': self.hochberg(pvalues),
            'Bonferroni': self.bonferroni(pvalues),
        }
        
        qvals, pi0 = self.storey_qvalue(pvalues)
        results['qvalue'] = qvals
        
        df_comp = pd.DataFrame(results)
        
        # Add summary row with counts at padj < 0.05
        summary = {
            'Method': ['BH', 'BY', 'Holm', 'Hochberg', 'Bonferroni', 'qvalue'],
            'Significant (padj<0.05)': [
                (df_comp['BH'] < 0.05).sum(),
                (df_comp['BY'] < 0.05).sum(),
                (df_comp['Holm'] < 0.05).sum(),
                (df_comp['Hochberg'] < 0.05).sum(),
                (df_comp['Bonferroni'] < 0.05).sum(),
                (df_comp['qvalue'] < 0.05).sum(),
            ]
        }
        
        if self.verbose:
            print("\nP-value adjustment method comparison (padj < 0.05):")
            for m, c in zip(summary['Method'], summary['Significant (padj<0.05)']):
                print(f"  {m}: {c} genes")
        
        return df_comp


# Convenience function
def optimize_thresholds(
    df: pd.DataFrame,
    goal: str = 'discovery',
    logfc_method: str = 'auto',
    n1: int = 3,
    n2: int = 3,
    verbose: bool = True
) -> ThresholdResult:
    """
    Convenience function for quick threshold optimization.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results with logFC and pvalue columns
    goal : str
        'discovery', 'validation', or 'balanced'
    logfc_method : str
        'auto', 'mad', 'mixture', 'power', or 'percentile'
    n1, n2 : int
        Sample sizes per group
    verbose : bool
        Print progress
        
    Returns
    -------
    ThresholdResult
        Optimization results
        
    Example
    -------
    >>> result = optimize_thresholds(df, goal='discovery')
    >>> print(result.summary())
    """
    ato = AdaptiveThresholdOptimizer(df, goal=goal, verbose=verbose)
    return ato.optimize(logfc_method=logfc_method, n1=n1, n2=n2)
