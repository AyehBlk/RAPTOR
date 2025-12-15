"""
Visualization module for Adaptive Threshold Optimizer.

Provides publication-ready plots for threshold optimization results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from typing import Optional, Tuple, List
import warnings


def plot_logfc_distribution(
    df: pd.DataFrame,
    logfc_col: str = 'logfc',
    padj_col: str = 'padj',
    optimized_cutoff: float = None,
    traditional_cutoff: float = 1.0,
    figsize: Tuple[int, int] = (12, 5),
    save_path: str = None
) -> plt.Figure:
    """
    Plot logFC distribution with cutoff comparison.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results
    logfc_col : str
        Column name for log2 fold change
    padj_col : str
        Column name for adjusted p-values
    optimized_cutoff : float
        Optimized |logFC| cutoff
    traditional_cutoff : float
        Traditional cutoff (default 1.0)
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.Figure
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    logfc = df[logfc_col].values
    
    # Left plot: LogFC histogram
    ax1 = axes[0]
    ax1.hist(logfc, bins=100, color='steelblue', alpha=0.7, edgecolor='white')
    ax1.axvline(0, color='black', linestyle='-', linewidth=1)
    
    if optimized_cutoff:
        ax1.axvline(optimized_cutoff, color='green', linestyle='--', linewidth=2, 
                   label=f'Optimized: ±{optimized_cutoff:.2f}')
        ax1.axvline(-optimized_cutoff, color='green', linestyle='--', linewidth=2)
    
    ax1.axvline(traditional_cutoff, color='red', linestyle=':', linewidth=2,
               label=f'Traditional: ±{traditional_cutoff:.2f}')
    ax1.axvline(-traditional_cutoff, color='red', linestyle=':', linewidth=2)
    
    ax1.set_xlabel('log2 Fold Change', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('LogFC Distribution with Cutoffs', fontsize=14)
    ax1.legend(loc='upper right')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Right plot: Cumulative distribution
    ax2 = axes[1]
    abs_logfc = np.abs(logfc)
    sorted_logfc = np.sort(abs_logfc)
    cumulative = np.arange(1, len(sorted_logfc) + 1) / len(sorted_logfc)
    
    ax2.plot(sorted_logfc, 1 - cumulative, color='steelblue', linewidth=2)
    ax2.set_xlabel('|log2 Fold Change| cutoff', fontsize=12)
    ax2.set_ylabel('Proportion of genes above cutoff', fontsize=12)
    ax2.set_title('Cumulative Distribution', fontsize=14)
    
    if optimized_cutoff:
        prop_opt = (abs_logfc > optimized_cutoff).mean()
        ax2.axvline(optimized_cutoff, color='green', linestyle='--', linewidth=2)
        ax2.axhline(prop_opt, color='green', linestyle='--', alpha=0.5)
        ax2.annotate(f'{prop_opt:.1%}', xy=(optimized_cutoff, prop_opt),
                    xytext=(optimized_cutoff + 0.3, prop_opt + 0.05),
                    fontsize=10, color='green')
    
    prop_trad = (abs_logfc > traditional_cutoff).mean()
    ax2.axvline(traditional_cutoff, color='red', linestyle=':', linewidth=2)
    ax2.axhline(prop_trad, color='red', linestyle=':', alpha=0.5)
    ax2.annotate(f'{prop_trad:.1%}', xy=(traditional_cutoff, prop_trad),
                xytext=(traditional_cutoff + 0.3, prop_trad + 0.05),
                fontsize=10, color='red')
    
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xlim(0, min(3, sorted_logfc.max()))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_pvalue_distribution(
    df: pd.DataFrame,
    pvalue_col: str = 'pvalue',
    pi0: float = None,
    figsize: Tuple[int, int] = (10, 5),
    save_path: str = None
) -> plt.Figure:
    """
    Plot p-value distribution with π₀ estimation.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results
    pvalue_col : str
        Column name for p-values
    pi0 : float, optional
        Estimated π₀
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.Figure
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    pvals = df[pvalue_col].values
    pvals = pvals[~np.isnan(pvals)]
    
    # Left: Histogram
    ax1 = axes[0]
    counts, bins, _ = ax1.hist(pvals, bins=50, color='steelblue', alpha=0.7, 
                               edgecolor='white', density=True)
    
    if pi0:
        ax1.axhline(pi0, color='red', linestyle='--', linewidth=2,
                   label=f'π₀ = {pi0:.3f}')
        ax1.legend()
    
    ax1.set_xlabel('P-value', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('P-value Distribution', fontsize=14)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Right: QQ plot
    ax2 = axes[1]
    n = len(pvals)
    expected = np.arange(1, n + 1) / (n + 1)
    observed = np.sort(pvals)
    
    ax2.scatter(-np.log10(expected), -np.log10(observed), 
               alpha=0.5, s=10, color='steelblue')
    
    # Diagonal line
    max_val = max(-np.log10(expected).max(), -np.log10(observed[observed > 0]).max())
    ax2.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Expected')
    
    ax2.set_xlabel('Expected -log10(p)', fontsize=12)
    ax2.set_ylabel('Observed -log10(p)', fontsize=12)
    ax2.set_title('QQ Plot', fontsize=14)
    ax2.legend()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_volcano(
    df: pd.DataFrame,
    logfc_col: str = 'logfc',
    padj_col: str = 'padj',
    logfc_cutoff: float = 1.0,
    padj_cutoff: float = 0.05,
    figsize: Tuple[int, int] = (10, 8),
    highlight_genes: List[str] = None,
    save_path: str = None
) -> plt.Figure:
    """
    Create volcano plot with significance thresholds.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results (must have index as gene names)
    logfc_col : str
        Column name for log2 fold change
    padj_col : str
        Column name for adjusted p-values
    logfc_cutoff : float
        |logFC| significance cutoff
    padj_cutoff : float
        Adjusted p-value cutoff
    figsize : tuple
        Figure size
    highlight_genes : list, optional
        Gene names to highlight
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    logfc = df[logfc_col].values
    padj = df[padj_col].values
    
    # Replace 0 padj values with minimum non-zero
    min_padj = padj[padj > 0].min()
    padj = np.where(padj == 0, min_padj / 10, padj)
    
    neg_log_padj = -np.log10(padj)
    
    # Classify points
    sig_up = (logfc > logfc_cutoff) & (padj < padj_cutoff)
    sig_down = (logfc < -logfc_cutoff) & (padj < padj_cutoff)
    not_sig = ~(sig_up | sig_down)
    
    # Plot points
    ax.scatter(logfc[not_sig], neg_log_padj[not_sig], 
              c='gray', alpha=0.5, s=20, label='Not significant')
    ax.scatter(logfc[sig_up], neg_log_padj[sig_up], 
              c='firebrick', alpha=0.7, s=30, label=f'Up ({sig_up.sum()})')
    ax.scatter(logfc[sig_down], neg_log_padj[sig_down], 
              c='steelblue', alpha=0.7, s=30, label=f'Down ({sig_down.sum()})')
    
    # Add threshold lines
    ax.axhline(-np.log10(padj_cutoff), color='black', linestyle='--', 
              linewidth=1, alpha=0.7)
    ax.axvline(logfc_cutoff, color='black', linestyle='--', 
              linewidth=1, alpha=0.7)
    ax.axvline(-logfc_cutoff, color='black', linestyle='--', 
              linewidth=1, alpha=0.7)
    
    # Highlight specific genes
    if highlight_genes:
        for gene in highlight_genes:
            if gene in df.index:
                x = df.loc[gene, logfc_col]
                y = -np.log10(df.loc[gene, padj_col])
                ax.annotate(gene, xy=(x, y), xytext=(5, 5),
                           textcoords='offset points', fontsize=8,
                           arrowprops=dict(arrowstyle='->', color='black'))
    
    ax.set_xlabel('log2 Fold Change', fontsize=12)
    ax.set_ylabel('-log10(adjusted p-value)', fontsize=12)
    ax.set_title(f'Volcano Plot (|logFC| > {logfc_cutoff}, padj < {padj_cutoff})', 
                fontsize=14)
    ax.legend(loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_threshold_comparison(
    df: pd.DataFrame,
    logfc_col: str = 'logfc',
    padj_col: str = 'padj',
    logfc_values: List[float] = [0.5, 1.0, 1.5, 2.0],
    padj_values: List[float] = [0.01, 0.05, 0.1],
    optimized_logfc: float = None,
    optimized_padj: float = None,
    figsize: Tuple[int, int] = (10, 6),
    save_path: str = None
) -> plt.Figure:
    """
    Heatmap comparing gene counts across threshold combinations.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results
    logfc_col : str
        Column name for log2 fold change
    padj_col : str
        Column name for adjusted p-values
    logfc_values : list
        LogFC cutoffs to compare
    padj_values : list
        Padj cutoffs to compare
    optimized_logfc : float, optional
        Optimized logFC cutoff to highlight
    optimized_padj : float, optional
        Optimized padj cutoff to highlight
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.Figure
    """
    # Calculate gene counts for each combination
    counts = np.zeros((len(padj_values), len(logfc_values)))
    
    for i, padj_cut in enumerate(padj_values):
        for j, logfc_cut in enumerate(logfc_values):
            n_sig = (
                (np.abs(df[logfc_col]) > logfc_cut) & 
                (df[padj_col] < padj_cut)
            ).sum()
            counts[i, j] = n_sig
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    im = ax.imshow(counts, cmap='YlOrRd', aspect='auto')
    
    # Add text annotations
    for i in range(len(padj_values)):
        for j in range(len(logfc_values)):
            text = ax.text(j, i, int(counts[i, j]),
                          ha='center', va='center', fontsize=12,
                          color='white' if counts[i, j] > counts.max()/2 else 'black')
    
    # Mark optimized thresholds
    if optimized_logfc and optimized_padj:
        # Find closest values
        logfc_idx = np.argmin(np.abs(np.array(logfc_values) - optimized_logfc))
        padj_idx = np.argmin(np.abs(np.array(padj_values) - optimized_padj))
        
        rect = mpatches.Rectangle((logfc_idx - 0.5, padj_idx - 0.5), 1, 1,
                                   fill=False, edgecolor='green', linewidth=3)
        ax.add_patch(rect)
    
    # Labels
    ax.set_xticks(range(len(logfc_values)))
    ax.set_xticklabels([f'{v}' for v in logfc_values])
    ax.set_yticks(range(len(padj_values)))
    ax.set_yticklabels([f'{v}' for v in padj_values])
    
    ax.set_xlabel('|log2 Fold Change| cutoff', fontsize=12)
    ax.set_ylabel('Adjusted p-value cutoff', fontsize=12)
    ax.set_title('Number of DE Genes by Threshold Combination', fontsize=14)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Number of DE genes', fontsize=11)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_adjustment_comparison(
    df: pd.DataFrame,
    pvalue_col: str = 'pvalue',
    figsize: Tuple[int, int] = (12, 5),
    save_path: str = None
) -> plt.Figure:
    """
    Compare different p-value adjustment methods.
    
    Parameters
    ----------
    df : pd.DataFrame
        DE results
    pvalue_col : str
        Column name for p-values
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.Figure
    """
    from .ato import AdaptiveThresholdOptimizer
    
    ato = AdaptiveThresholdOptimizer(df, verbose=False)
    pvals = df[pvalue_col].values
    
    # Calculate all adjustments
    methods = {
        'BH': ato.benjamini_hochberg(pvals),
        'BY': ato.benjamini_yekutieli(pvals),
        'Holm': ato.holm(pvals),
        'Hochberg': ato.hochberg(pvals),
        'Bonferroni': ato.bonferroni(pvals),
    }
    qvals, _ = ato.storey_qvalue(pvals)
    methods['q-value'] = qvals
    
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Left: Bar chart of significant genes
    ax1 = axes[0]
    sig_counts = {m: (v < 0.05).sum() for m, v in methods.items()}
    
    colors = ['steelblue', 'darkorange', 'green', 'red', 'purple', 'brown']
    bars = ax1.bar(sig_counts.keys(), sig_counts.values(), color=colors)
    
    ax1.set_ylabel('Number of significant genes (padj < 0.05)', fontsize=11)
    ax1.set_title('Method Comparison', fontsize=14)
    ax1.tick_params(axis='x', rotation=45)
    
    # Add value labels
    for bar, count in zip(bars, sig_counts.values()):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(count), ha='center', fontsize=10)
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Right: Scatter comparison of BH vs q-value
    ax2 = axes[1]
    ax2.scatter(methods['BH'], methods['q-value'], alpha=0.5, s=10)
    ax2.plot([0, 1], [0, 1], 'r--', linewidth=1)
    ax2.set_xlabel('BH adjusted p-value', fontsize=11)
    ax2.set_ylabel('q-value', fontsize=11)
    ax2.set_title('BH vs q-value Comparison', fontsize=14)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_optimization_summary(
    result,  # ThresholdResult
    df: pd.DataFrame,
    logfc_col: str = 'logfc',
    padj_col: str = 'padj',
    figsize: Tuple[int, int] = (14, 10),
    save_path: str = None
) -> plt.Figure:
    """
    Create comprehensive summary plot of optimization results.
    
    Parameters
    ----------
    result : ThresholdResult
        Output from AdaptiveThresholdOptimizer.optimize()
    df : pd.DataFrame
        DE results
    logfc_col : str
        Column name for log2 fold change
    padj_col : str
        Column name for adjusted p-values
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    matplotlib.Figure
    """
    fig = plt.figure(figsize=figsize)
    
    # Create grid
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    logfc = df[logfc_col].values
    padj = df[padj_col].values
    
    # 1. LogFC distribution (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(logfc, bins=80, color='steelblue', alpha=0.7, edgecolor='white')
    ax1.axvline(result.logfc_cutoff, color='green', linestyle='--', linewidth=2,
               label=f'Optimized: ±{result.logfc_cutoff:.2f}')
    ax1.axvline(-result.logfc_cutoff, color='green', linestyle='--', linewidth=2)
    ax1.axvline(1.0, color='red', linestyle=':', linewidth=2,
               label='Traditional: ±1.0')
    ax1.axvline(-1.0, color='red', linestyle=':', linewidth=2)
    ax1.set_xlabel('log2 Fold Change')
    ax1.set_ylabel('Count')
    ax1.set_title('LogFC Distribution')
    ax1.legend()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 2. P-value distribution (top right)
    ax2 = fig.add_subplot(gs[0, 1])
    pvals_clean = df['pvalue'].values[~np.isnan(df['pvalue'].values)]
    ax2.hist(pvals_clean, bins=50, color='steelblue', alpha=0.7, edgecolor='white')
    ax2.axhline(result.pi0_estimate * len(pvals_clean) / 50, color='red', 
               linestyle='--', linewidth=2, label=f'π₀ = {result.pi0_estimate:.3f}')
    ax2.set_xlabel('P-value')
    ax2.set_ylabel('Count')
    ax2.set_title(f'P-value Distribution (π₀ = {result.pi0_estimate:.3f})')
    ax2.legend()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # 3. Volcano plot (bottom left)
    ax3 = fig.add_subplot(gs[1, 0])
    min_padj = padj[padj > 0].min() if (padj > 0).any() else 1e-300
    padj_safe = np.where(padj == 0, min_padj / 10, padj)
    neg_log_padj = -np.log10(padj_safe)
    
    sig_up = (logfc > result.logfc_cutoff) & (padj < result.padj_cutoff)
    sig_down = (logfc < -result.logfc_cutoff) & (padj < result.padj_cutoff)
    not_sig = ~(sig_up | sig_down)
    
    ax3.scatter(logfc[not_sig], neg_log_padj[not_sig], c='gray', alpha=0.3, s=15)
    ax3.scatter(logfc[sig_up], neg_log_padj[sig_up], c='firebrick', alpha=0.7, s=25,
               label=f'Up ({sig_up.sum()})')
    ax3.scatter(logfc[sig_down], neg_log_padj[sig_down], c='steelblue', alpha=0.7, s=25,
               label=f'Down ({sig_down.sum()})')
    ax3.axhline(-np.log10(result.padj_cutoff), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(result.logfc_cutoff, color='black', linestyle='--', alpha=0.5)
    ax3.axvline(-result.logfc_cutoff, color='black', linestyle='--', alpha=0.5)
    ax3.set_xlabel('log2 Fold Change')
    ax3.set_ylabel('-log10(adjusted p-value)')
    ax3.set_title('Volcano Plot (Optimized Thresholds)')
    ax3.legend()
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    # 4. Summary text (bottom right)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')
    
    summary_text = f"""
    ADAPTIVE THRESHOLD OPTIMIZER RESULTS
    ════════════════════════════════════
    
    RECOMMENDED THRESHOLDS:
    • |log2FC| cutoff: {result.logfc_cutoff:.3f}
    • Adjusted p-value cutoff: {result.padj_cutoff}
    • P-value method: {result.padj_method}
    
    RESULTS:
    • DE genes (optimized): {result.n_significant_optimized}
    • DE genes (traditional): {result.n_significant_traditional}
    • Change: {result.n_significant_optimized - result.n_significant_traditional:+d} genes
    
    STATISTICS:
    • π₀ estimate: {result.pi0_estimate:.3f}
    • True DE genes: ~{(1-result.pi0_estimate)*100:.1f}%
    
    METHOD USED:
    • LogFC: {result.logfc_method}
    """
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Threshold Optimization Summary', fontsize=16, fontweight='bold', y=1.02)
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig
