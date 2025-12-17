"""
RAPTOR Threshold Optimizer Module

Adaptive Threshold Optimizer (ATO) for RNA-seq Differential Expression Analysis.
Provides data-driven methods for optimizing significance thresholds.

Example
-------
>>> from raptor.threshold_optimizer import AdaptiveThresholdOptimizer, optimize_thresholds
>>> 
>>> # Quick usage
>>> result = optimize_thresholds(df, goal='discovery')
>>> print(result.summary())
>>> 
>>> # Full control
>>> ato = AdaptiveThresholdOptimizer(df, goal='discovery')
>>> result = ato.optimize()
>>> significant_genes = ato.get_significant_genes()
"""

from .ato import (
    AdaptiveThresholdOptimizer,
    ThresholdResult,
    optimize_thresholds
)

from .visualization import (
    plot_logfc_distribution,
    plot_pvalue_distribution,
    plot_volcano,
    plot_threshold_comparison,
    plot_adjustment_comparison,
    plot_optimization_summary
)

__all__ = [
    'AdaptiveThresholdOptimizer',
    'ThresholdResult',
    'optimize_thresholds',
    'plot_logfc_distribution',
    'plot_pvalue_distribution',
    'plot_volcano',
    'plot_threshold_comparison',
    'plot_adjustment_comparison',
    'plot_optimization_summary'
]

__version__ = '0.1.0'
