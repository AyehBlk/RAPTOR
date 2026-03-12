#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Pipeline Recommender Module

ML-based pipeline recommendation system for RNA-seq differential expression.

This module recommends the optimal DE pipeline based on data characteristics
extracted by the profiler. Recommendations are based on extensive benchmarking
literature and empirical evidence.

============================================================================
DECISION FRAMEWORK (Literature-Based)
============================================================================

The recommendation system uses a hierarchical decision tree based on:

1. SAMPLE SIZE (Most Critical)
   - n < 3 per group: Limited options, DESeq2 preferred (conservative)
   - n = 3-7 per group: DESeq2 or edgeR (both appropriate)
   - n ≥ 8 per group: All methods work; Wilcoxon may have better FDR control
   - n > 20 per group: limma-voom recommended for efficiency

2. DISPERSION / BCV
   - Low BCV (<0.2): limma-voom works well
   - Moderate BCV (0.2-0.4): DESeq2 or edgeR
   - High BCV (>0.4): edgeR (handles overdispersion better)

3. OUTLIERS
   - Outliers + small samples: edgeR_robust
   - Outliers + large samples: limma-voom (robust)
   - No outliers: Standard methods

4. LOW COUNT GENES
   - Many low counts (>30%): edgeR preferred
   - Few low counts: All methods appropriate

5. DESIGN COMPLEXITY
   - Simple (2 groups): Any method
   - Complex (multiple factors): limma-voom preferred

============================================================================
PIPELINE CHARACTERISTICS
============================================================================

DESeq2:
- Best for: Small samples, batch effects, general use
- Statistical model: Negative binomial with shrinkage estimation
- Normalization: Median-of-ratios (RLE)
- Strengths: Conservative, handles batch effects, good documentation
- Weaknesses: Slower for very large datasets

edgeR:
- Best for: Low counts, overdispersed data, small samples
- Statistical model: Negative binomial with empirical Bayes
- Normalization: TMM (Trimmed Mean of M-values)
- Strengths: Fast, handles low counts well, robust mode available
- Weaknesses: Can be anti-conservative for some datasets

limma-voom:
- Best for: Large samples, complex designs, speed
- Statistical model: Linear model with precision weights
- Normalization: TMM + log-CPM transformation
- Strengths: Very fast, robust, handles complex designs
- Weaknesses: Needs sufficient replication (≥3 per group)

Wilcoxon (TMM + rank test):
- Best for: Very large samples (n ≥ 8 per group), FDR control
- Statistical model: Non-parametric
- Normalization: TMM
- Strengths: Best FDR control for large samples (Li et al. 2022)
- Weaknesses: Lower power for small samples, no covariate adjustment

============================================================================

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
import json

# RAPTOR v2.2.0 imports
from raptor.profiler import DataProfile, RNAseqDataProfiler

import logging

from raptor.utils.validation import (
    validate_count_matrix,
    validate_metadata,
    validate_group_column,
    validate_file_path,
    validate_directory_path,
    validate_numeric_range,
    validate_positive_integer,
    validate_probability
)

from raptor.utils.errors import (
    handle_errors,
    ValidationError,
    DependencyError,
    check_file_exists,
    validate_output_writable
)

logger = logging.getLogger(__name__)


# =============================================================================
# CLI PARAMETERS & MODULE METADATA - v2.2.0
# =============================================================================

RECOMMENDER_CLI_PARAMS = {
    'counts': {
        'flag': '--counts',
        'short': '-c',
        'required': True,
        'type': 'file',
        'help': 'Count matrix CSV file (genes x samples)',
        'validation': 'Must be CSV with gene IDs as rows, sample IDs as columns'
    },
    'metadata': {
        'flag': '--metadata',
        'short': '-m',
        'type': 'file',
        'help': 'Sample metadata CSV (optional, enables group-aware recommendations)',
        'validation': 'Must contain group/condition column'
    },
    'group_column': {
        'flag': '--group-column',
        'short': '-g',
        'default': 'condition',
        'type': 'str',
        'help': 'Column name for grouping samples (default: condition)'
    },
    'profile_file': {
        'flag': '--profile',
        'short': '-p',
        'type': 'file',
        'help': 'Pre-computed data profile JSON (alternative to counts/metadata)',
        'validation': 'Must be valid profile JSON from profiler module'
    },
    'output': {
        'flag': '--output',
        'short': '-o',
        'default': 'recommendation.json',
        'type': 'file',
        'help': 'Output recommendation JSON file'
    },
    'verbose': {
        'flag': '--verbose',
        'short': '-v',
        'action': 'store_true',
        'help': 'Print detailed recommendation reasoning'
    }
}

# Module metadata
MODULE_NAME = "Pipeline Recommender"
MODULE_VERSION = "2.2.0"
MODULE_ID = "M4"
MODULE_STAGE = 4

__all__ = [
    # Main class
    'PipelineRecommender',
    
    # Data class
    'Recommendation',
    
    # Convenience function
    'recommend_pipeline',
    
    # CLI integration
    'RECOMMENDER_CLI_PARAMS',
    'MODULE_NAME',
    'MODULE_VERSION',
    'MODULE_ID',
    'MODULE_STAGE'
]


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class PipelineInfo:
    """Information about a DE pipeline."""
    name: str
    full_name: str
    description: str
    model: str
    normalization: str
    strengths: List[str]
    weaknesses: List[str]
    min_replicates: int
    best_for: List[str]
    r_packages: List[str]
    reference: str


@dataclass
class Recommendation:
    """Pipeline recommendation with explanation."""
    
    # Primary recommendation
    primary_pipeline: str
    primary_score: float  # 0-100 confidence
    primary_reason: str
    
    # Alternative recommendation
    alternative_pipeline: str
    alternative_score: float
    alternative_reason: str
    
    # Third option if applicable
    third_option: Optional[str] = None
    third_score: Optional[float] = None
    
    # Decision factors
    decision_factors: Dict[str, Any] = field(default_factory=dict)
    
    # Warnings
    warnings: List[str] = field(default_factory=list)
    
    # R code snippets
    r_code_primary: str = ""
    r_code_alternative: str = ""
    
    # Full ranking of all pipelines
    all_scores: Dict[str, float] = field(default_factory=dict)
    
    def summary(self) -> str:
        """Generate human-readable recommendation summary."""
        lines = [
            "",
            "╔" + "═" * 68 + "╗",
            "║" + "  🦖 RAPTOR PIPELINE RECOMMENDATION".center(68) + "║",
            "╠" + "═" * 68 + "╣",
            "",
            f"  🥇 PRIMARY RECOMMENDATION: {self.primary_pipeline}",
            f"     Confidence: {self.primary_score:.0f}%",
            f"     Reason: {self.primary_reason}",
            "",
            f"  🥈 ALTERNATIVE: {self.alternative_pipeline}",
            f"     Confidence: {self.alternative_score:.0f}%",
            f"     Reason: {self.alternative_reason}",
        ]
        
        if self.third_option:
            lines.extend([
                "",
                f"  🥉 THIRD OPTION: {self.third_option}",
                f"     Confidence: {self.third_score:.0f}%"
            ])
        
        if self.warnings:
            lines.extend([
                "",
                "  ⚠️  WARNINGS:",
            ])
            for w in self.warnings:
                lines.append(f"     • {w}")
        
        lines.extend([
            "",
            "  📊 DECISION FACTORS:",
        ])
        for factor, value in self.decision_factors.items():
            if isinstance(value, float):
                lines.append(f"     {factor}: {value:.3f}")
            else:
                lines.append(f"     {factor}: {value}")
        
        lines.extend([
            "",
            "╚" + "═" * 68 + "╝",
            ""
        ])
        
        return "\n".join(lines)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'primary_pipeline': self.primary_pipeline,
            'primary_score': self.primary_score,
            'primary_reason': self.primary_reason,
            'alternative_pipeline': self.alternative_pipeline,
            'alternative_score': self.alternative_score,
            'alternative_reason': self.alternative_reason,
            'third_option': self.third_option,
            'third_score': self.third_score,
            'decision_factors': self.decision_factors,
            'warnings': self.warnings,
            'all_scores': self.all_scores
        }


# =============================================================================
# PIPELINE DEFINITIONS
# =============================================================================

PIPELINES = {
    'DESeq2': PipelineInfo(
        name='DESeq2',
        full_name='DESeq2 (Differential Expression analysis for Sequence count data 2)',
        description='Negative binomial GLM with shrinkage estimation for dispersion and fold changes',
        model='Negative Binomial with shrinkage',
        normalization='Median-of-ratios (RLE)',
        strengths=[
            'Conservative - good FDR control',
            'Handles batch effects well',
            'Automatic outlier detection (Cook\'s distance)',
            'Shrinkage of fold changes reduces noise',
            'Excellent documentation and support'
        ],
        weaknesses=[
            'Can be slow for very large datasets',
            'May be too conservative for some applications',
            'Requires ≥2 replicates per condition'
        ],
        min_replicates=2,
        best_for=[
            'Small sample sizes (3-10 per group)',
            'Studies with batch effects',
            'General-purpose RNA-seq analysis',
            'When reliability is more important than power'
        ],
        r_packages=['DESeq2'],
        reference='Love et al. (2014) Genome Biology 15:550'
    ),
    
    'edgeR': PipelineInfo(
        name='edgeR',
        full_name='edgeR (empirical analysis of DGE in R)',
        description='Negative binomial GLM with empirical Bayes moderation of dispersion',
        model='Negative Binomial with tagwise dispersion',
        normalization='TMM (Trimmed Mean of M-values)',
        strengths=[
            'Fast computation',
            'Handles low counts well',
            'Good for overdispersed data',
            'Robust mode available (edgeR_robust)',
            'Flexible GLM framework'
        ],
        weaknesses=[
            'Can be anti-conservative in some cases',
            'May have inflated FDR for large samples'
        ],
        min_replicates=2,
        best_for=[
            'Datasets with many low-count genes',
            'Highly overdispersed data',
            'When outliers are a concern (use robust)',
            'Fast exploratory analysis'
        ],
        r_packages=['edgeR'],
        reference='Robinson et al. (2010) Bioinformatics 26:139-140'
    ),
    
    'limma-voom': PipelineInfo(
        name='limma-voom',
        full_name='limma with voom transformation',
        description='Linear model with precision weights estimated from mean-variance trend',
        model='Linear model with empirical Bayes',
        normalization='TMM + log-CPM transformation',
        strengths=[
            'Very fast',
            'Excellent for large datasets',
            'Handles complex designs well',
            'Robust to outliers in large samples',
            'Well-established methodology'
        ],
        weaknesses=[
            'Needs sufficient replication (≥3 per group)',
            'Sensitive to outliers in small samples',
            'Assumes approximate normality after transformation'
        ],
        min_replicates=3,
        best_for=[
            'Large sample sizes (>20 per group)',
            'Complex experimental designs',
            'Multi-factor experiments',
            'When speed is important'
        ],
        r_packages=['limma', 'edgeR'],
        reference='Law et al. (2014) Genome Biology 15:R29'
    ),
    
    'Wilcoxon': PipelineInfo(
        name='Wilcoxon',
        full_name='Wilcoxon rank-sum test with TMM normalization',
        description='Non-parametric rank-based test after TMM normalization',
        model='Non-parametric (rank-based)',
        normalization='TMM',
        strengths=[
            'Best FDR control for large samples',
            'No distributional assumptions',
            'Robust to outliers',
            'Simple and interpretable'
        ],
        weaknesses=[
            'Lower power for small samples (<8 per group)',
            'Cannot adjust for covariates',
            'Requires large sample sizes'
        ],
        min_replicates=8,
        best_for=[
            'Large sample sizes (≥8 per group)',
            'Population-level studies',
            'When FDR control is critical',
            'Validation of parametric results'
        ],
        r_packages=['edgeR', 'stats'],
        reference='Li et al. (2022) Genome Biology 23:79'
    ),
    
    'edgeR_robust': PipelineInfo(
        name='edgeR_robust',
        full_name='edgeR with robust dispersion estimation',
        description='edgeR with observation weights to downweight outliers',
        model='Negative Binomial with robust estimation',
        normalization='TMM',
        strengths=[
            'Handles outliers well',
            'Maintains power while controlling FDR',
            'Better for heterogeneous datasets'
        ],
        weaknesses=[
            'Slightly more conservative than standard edgeR',
            'Requires edgeR version ≥ 3.8'
        ],
        min_replicates=2,
        best_for=[
            'Datasets with known outliers',
            'Small samples with suspected outliers',
            'Heterogeneous experimental conditions'
        ],
        r_packages=['edgeR'],
        reference='Zhou et al. (2014) Nucleic Acids Research'
    )
}


# =============================================================================
# MAIN RECOMMENDER CLASS
# =============================================================================

class PipelineRecommender:
    """
    ML-based pipeline recommendation system.
    
    Recommends the optimal DE pipeline based on data profile features.
    Uses a combination of rule-based decisions and scoring functions
    derived from benchmarking literature.
    
    Parameters
    ----------
    profile : DataProfile
        Data profile from RNAseqDataProfiler.
    
    Examples
    --------
    >>> profiler = RNAseqDataProfiler(counts, metadata)
    >>> profile = profiler.run_full_profile()
    >>> recommender = PipelineRecommender(profile)
    >>> recommendation = recommender.get_recommendation()
    >>> print(recommendation.summary())
    """
    
    @handle_errors(exit_on_error=False)
    def __init__(self, profile: DataProfile):
        """Initialize recommender with data profile."""
        # Validate profile
        if not profile:
            raise ValidationError("Profile cannot be empty")
        
        # Validate profile has required fields
        required_fields = ['n_samples', 'n_groups', 'bcv', 'min_group_size']
        if isinstance(profile, dict):
            missing = set(required_fields) - set(profile.keys())
            if missing:
                raise ValidationError(f"Profile missing required fields: {missing}")
        elif hasattr(profile, '__dict__'):
            # DataProfile object
            for field in required_fields:
                if not hasattr(profile, field) or getattr(profile, field) is None:
                    raise ValidationError(f"Profile missing required field: {field}")
        
        self.profile = profile
        self.scores = {}
        self.warnings = []
    
    def get_recommendation(self) -> Recommendation:
        """
        Generate pipeline recommendation.
        
        Returns
        -------
        Recommendation
            Complete recommendation with explanations.
        """
        logger.info("Generating pipeline recommendation...")
        
        # Calculate scores for each pipeline
        self._calculate_scores()
        
        # Generate warnings
        self._generate_warnings()
        
        # Sort pipelines by score
        sorted_pipelines = sorted(
            self.scores.items(), 
            key=lambda x: x[1], 
            reverse=True
        )
        
        # Create recommendation
        primary = sorted_pipelines[0]
        alternative = sorted_pipelines[1]
        third = sorted_pipelines[2] if len(sorted_pipelines) > 2 else (None, None)
        
        recommendation = Recommendation(
            primary_pipeline=primary[0],
            primary_score=primary[1],
            primary_reason=self._get_reason(primary[0]),
            alternative_pipeline=alternative[0],
            alternative_score=alternative[1],
            alternative_reason=self._get_reason(alternative[0]),
            third_option=third[0],
            third_score=third[1],
            decision_factors=self._get_decision_factors(),
            warnings=self.warnings,
            all_scores=self.scores
        )
        
        # Add R code
        recommendation.r_code_primary = self._get_r_code(primary[0])
        recommendation.r_code_alternative = self._get_r_code(alternative[0])
        
        return recommendation
    
    def _calculate_scores(self):
        """Calculate recommendation scores for each pipeline."""
        p = self.profile
        
        # Initialize base scores
        self.scores = {
            'DESeq2': 70.0,
            'edgeR': 70.0,
            'limma-voom': 70.0,
            'Wilcoxon': 50.0,  # Lower base - only good for specific cases
            'edgeR_robust': 60.0  # Lower base - specialized
        }
        
        # =====================================================================
        # SAMPLE SIZE SCORING (Most Critical Factor)
        # =====================================================================
        min_n = p.min_group_size
        
        if min_n < 3:
            # Very small samples - DESeq2 most conservative
            self.scores['DESeq2'] += 20
            self.scores['edgeR'] += 10
            self.scores['limma-voom'] -= 20  # Not recommended
            self.scores['Wilcoxon'] -= 30  # Definitely not
        elif min_n >= 3 and min_n < 8:
            # Standard small samples - all parametric methods OK
            self.scores['DESeq2'] += 15
            self.scores['edgeR'] += 15
            self.scores['limma-voom'] += 10
            self.scores['Wilcoxon'] -= 20  # Still too small
        elif min_n >= 8 and min_n < 20:
            # Medium samples - all methods appropriate
            self.scores['DESeq2'] += 10
            self.scores['edgeR'] += 10
            self.scores['limma-voom'] += 15
            self.scores['Wilcoxon'] += 15  # Now appropriate
        else:
            # Large samples - limma-voom and Wilcoxon shine
            self.scores['limma-voom'] += 20
            self.scores['Wilcoxon'] += 25  # Best FDR control
            self.scores['DESeq2'] += 5  # Still OK but slower
            self.scores['edgeR'] += 5  # May have inflated FDR
        
        # =====================================================================
        # DISPERSION / BCV SCORING
        # =====================================================================
        bcv = p.bcv
        
        if bcv < 0.2:
            # Low dispersion - limma-voom appropriate
            self.scores['limma-voom'] += 10
            self.scores['DESeq2'] += 5
        elif bcv < 0.4:
            # Moderate dispersion - all NB methods good
            self.scores['DESeq2'] += 10
            self.scores['edgeR'] += 10
            self.scores['limma-voom'] += 5
        else:
            # High dispersion - edgeR handles best
            self.scores['edgeR'] += 15
            self.scores['edgeR_robust'] += 10
            self.scores['DESeq2'] += 5
            self.scores['limma-voom'] -= 5
        
        # =====================================================================
        # OUTLIER SCORING
        # =====================================================================
        if p.has_outliers:
            severity = p.outlier_severity
            
            if severity == 'mild':
                self.scores['edgeR_robust'] += 10
                self.scores['DESeq2'] += 5  # Has Cook's distance
            elif severity == 'moderate':
                self.scores['edgeR_robust'] += 20
                self.scores['DESeq2'] += 5
                self.scores['edgeR'] -= 5
                self.scores['limma-voom'] -= 5 if min_n < 10 else 5  # Robust for large
            else:  # severe
                self.scores['edgeR_robust'] += 25
                self.scores['edgeR'] -= 10
                self.scores['limma-voom'] -= 10 if min_n < 10 else 0
        
        # =====================================================================
        # LOW COUNT SCORING
        # =====================================================================
        low_counts = p.low_count_proportion
        
        if low_counts > 0.3:
            # Many low counts - edgeR handles well
            self.scores['edgeR'] += 10
            self.scores['edgeR_robust'] += 5
            self.scores['DESeq2'] += 5
        elif low_counts > 0.5:
            # Very many low counts
            self.scores['edgeR'] += 15
            self.scores['edgeR_robust'] += 10
            self.scores['limma-voom'] -= 5
        
        # =====================================================================
        # SPARSITY SCORING
        # =====================================================================
        if p.sparsity > 0.6:
            # High sparsity - NB models handle better
            self.scores['DESeq2'] += 5
            self.scores['edgeR'] += 5
        
        # =====================================================================
        # LIBRARY SIZE VARIATION SCORING
        # =====================================================================
        lib_cv = p.library_size_cv
        
        if lib_cv > 0.3:
            # High variation - normalization critical
            self.scores['DESeq2'] += 5  # RLE robust
            self.scores['edgeR'] += 5  # TMM robust
        
        if p.library_size_range > 5:
            # Very different library sizes
            self.scores['DESeq2'] += 5
        
        # =====================================================================
        # BATCH EFFECT SCORING
        # =====================================================================
        if p.has_batch_effect:
            # Batch effects - DESeq2 handles well
            self.scores['DESeq2'] += 10
            self.scores['limma-voom'] += 5  # Can include in model
            self.scores['Wilcoxon'] -= 10  # Can't adjust for covariates
            
            if p.batch_confounded:
                # Severe warning - no method can really fix this
                for key in self.scores:
                    self.scores[key] -= 10
        
        # =====================================================================
        # DESIGN COMPLEXITY SCORING
        # =====================================================================
        if p.design_complexity == 'complex':
            self.scores['limma-voom'] += 15
            self.scores['DESeq2'] += 5
            self.scores['Wilcoxon'] -= 15  # Can't handle complex designs
        
        # =====================================================================
        # NORMALIZE SCORES
        # =====================================================================
        # Cap at 0-100 range
        for key in self.scores:
            self.scores[key] = max(0, min(100, self.scores[key]))
        
        # Ensure minimum replicates are met
        for pipeline, info in PIPELINES.items():
            if min_n < info.min_replicates:
                self.scores[pipeline] = max(0, self.scores[pipeline] - 30)
    
    def _generate_warnings(self):
        """Generate warnings based on data characteristics."""
        p = self.profile
        
        if p.min_group_size < 3:
            self.warnings.append(
                f"Low replication (n={p.min_group_size}): Results may be unreliable. "
                f"Consider increasing biological replicates."
            )
        
        if p.bcv > 0.5:
            self.warnings.append(
                f"High biological variation (BCV={p.bcv:.2f}): "
                f"Consider using edgeR or DESeq2 for robust dispersion estimation."
            )
        
        if p.has_outliers:
            self.warnings.append(
                f"Outliers detected ({p.outlier_severity}): "
                f"Consider removing outliers or using robust methods."
            )
        
        if p.has_batch_effect:
            self.warnings.append(
                "Batch effect detected: Include batch as covariate in DE model."
            )
            if p.batch_confounded:
                self.warnings.append(
                    "⚠️ CRITICAL: Batch is confounded with condition! "
                    "Results may not separate biological from technical effects."
                )
        
        if p.sparsity > 0.7:
            self.warnings.append(
                f"High sparsity ({p.sparsity:.1%} zeros): "
                f"Consider this is expected for your data type (e.g., low-input RNA-seq)."
            )
        
        if p.library_size_range > 10:
            self.warnings.append(
                f"Large library size variation ({p.library_size_range:.1f}x): "
                f"Ensure normalization is working correctly."
            )
    
    def _get_decision_factors(self) -> Dict[str, Any]:
        """Get key decision factors."""
        p = self.profile
        return {
            'sample_size_per_group': p.min_group_size,
            'n_groups': p.n_groups,
            'bcv': p.bcv,
            'bcv_category': p.bcv_category,
            'library_size_cv': p.library_size_cv,
            'low_count_proportion': p.low_count_proportion,
            'sparsity': p.sparsity,
            'has_outliers': p.has_outliers,
            'has_batch_effect': p.has_batch_effect,
            'design_complexity': p.design_complexity
        }
    
    def _get_reason(self, pipeline: str) -> str:
        """Get reason for recommending a pipeline."""
        p = self.profile
        
        reasons = {
            'DESeq2': self._get_deseq2_reason(),
            'edgeR': self._get_edger_reason(),
            'limma-voom': self._get_limma_reason(),
            'Wilcoxon': self._get_wilcoxon_reason(),
            'edgeR_robust': self._get_edger_robust_reason()
        }
        
        return reasons.get(pipeline, "General-purpose recommendation")
    
    def _get_deseq2_reason(self) -> str:
        p = self.profile
        reasons = []
        
        if p.min_group_size < 8:
            reasons.append("appropriate for small samples")
        if p.has_batch_effect:
            reasons.append("handles batch effects well")
        if p.bcv >= 0.2 and p.bcv <= 0.4:
            reasons.append("suitable for moderate dispersion")
        
        if not reasons:
            reasons.append("reliable general-purpose choice")
        
        return "DESeq2: " + ", ".join(reasons)
    
    def _get_edger_reason(self) -> str:
        p = self.profile
        reasons = []
        
        if p.low_count_proportion > 0.3:
            reasons.append("handles low counts well")
        if p.bcv > 0.4:
            reasons.append("good for high dispersion")
        
        if not reasons:
            reasons.append("fast and flexible NB model")
        
        return "edgeR: " + ", ".join(reasons)
    
    def _get_limma_reason(self) -> str:
        p = self.profile
        reasons = []
        
        if p.min_group_size >= 20:
            reasons.append("efficient for large samples")
        if p.design_complexity == 'complex':
            reasons.append("handles complex designs")
        if p.bcv < 0.2:
            reasons.append("appropriate for low dispersion")
        
        if not reasons:
            reasons.append("fast linear modeling approach")
        
        return "limma-voom: " + ", ".join(reasons)
    
    def _get_wilcoxon_reason(self) -> str:
        p = self.profile
        reasons = []
        
        if p.min_group_size >= 8:
            reasons.append("excellent FDR control for large samples")
        reasons.append("non-parametric, robust")
        
        return "Wilcoxon: " + ", ".join(reasons)
    
    def _get_edger_robust_reason(self) -> str:
        p = self.profile
        reasons = []
        
        if p.has_outliers:
            reasons.append("handles outliers well")
        reasons.append("robust dispersion estimation")
        
        return "edgeR_robust: " + ", ".join(reasons)
    
    def _get_r_code(self, pipeline: str) -> str:
        """Generate R code snippet for the pipeline."""
        
        code_templates = {
            'DESeq2': '''
# DESeq2 Analysis
library(DESeq2)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition  # Add batch if needed: ~ batch + condition
)

# Pre-filter low counts
keep <- rowSums(counts(dds) >= 10) >= min(table(metadata$condition))
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results
results <- results(dds, alpha = 0.05)
results <- lfcShrink(dds, coef = 2, type = "apeglm")  # Recommended

# Significant genes
sig_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)
''',
            
            'edgeR': '''
# edgeR Analysis
library(edgeR)

# Create DGEList
dge <- DGEList(counts = counts, group = metadata$condition)

# Filter low counts
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge, method = "TMM")

# Estimate dispersion
design <- model.matrix(~ condition, data = metadata)
dge <- estimateDisp(dge, design)

# Test for DE
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)

# Get results
results <- topTags(qlf, n = Inf)$table
sig_genes <- subset(results, FDR < 0.05 & abs(logFC) > 1)
''',
            
            'limma-voom': '''
# limma-voom Analysis
library(limma)
library(edgeR)

# Create DGEList
dge <- DGEList(counts = counts)

# Filter low counts
keep <- filterByExpr(dge, group = metadata$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix
design <- model.matrix(~ condition, data = metadata)

# voom transformation
v <- voom(dge, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef = 2, n = Inf)
sig_genes <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
''',
            
            'Wilcoxon': '''
# Wilcoxon rank-sum test (TMM normalized)
library(edgeR)

# Create DGEList and normalize
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")

# Get normalized counts
norm_counts <- cpm(dge, log = TRUE, prior.count = 1)

# Wilcoxon test for each gene
group <- metadata$condition
group1 <- which(group == levels(factor(group))[1])
group2 <- which(group == levels(factor(group))[2])

pvalues <- apply(norm_counts, 1, function(x) {
    wilcox.test(x[group1], x[group2])$p.value
})

# Adjust p-values
padj <- p.adjust(pvalues, method = "BH")

# Log fold change
logFC <- rowMeans(norm_counts[, group2]) - rowMeans(norm_counts[, group1])

# Results
results <- data.frame(logFC = logFC, pvalue = pvalues, padj = padj)
sig_genes <- subset(results, padj < 0.05 & abs(logFC) > 1)
''',
            
            'edgeR_robust': '''
# edgeR robust Analysis
library(edgeR)

# Create DGEList
dge <- DGEList(counts = counts, group = metadata$condition)

# Filter low counts
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge, method = "TMM")

# Design matrix
design <- model.matrix(~ condition, data = metadata)

# Estimate dispersion with robust method
dge <- estimateGLMRobustDisp(dge, design)

# Fit with robust method
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)

# Get results
results <- topTags(lrt, n = Inf)$table
sig_genes <- subset(results, FDR < 0.05 & abs(logFC) > 1)
'''
        }
        
        return code_templates.get(pipeline, "# Code not available")


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def recommend_pipeline(profile: DataProfile) -> Recommendation:
    """
    Get pipeline recommendation from data profile.
    
    Parameters
    ----------
    profile : DataProfile
        Data profile from profiler.
    
    Returns
    -------
    Recommendation
        Pipeline recommendation with explanations.
    
    Examples
    --------
    >>> from raptor import profile_data_quick, recommend_pipeline
    >>> profile = profile_data_quick(counts, metadata)
    >>> recommendation = recommend_pipeline(profile)
    >>> print(recommendation.summary())
    """
    recommender = PipelineRecommender(profile)
    return recommender.get_recommendation()


def get_pipeline_info(pipeline_name: str) -> Optional[PipelineInfo]:
    """
    Get detailed information about a pipeline.
    
    Parameters
    ----------
    pipeline_name : str
        Name of pipeline (DESeq2, edgeR, limma-voom, Wilcoxon, edgeR_robust)
    
    Returns
    -------
    PipelineInfo or None
        Pipeline information or None if not found.
    """
    return PIPELINES.get(pipeline_name)


def list_pipelines() -> List[str]:
    """List all available pipelines."""
    return list(PIPELINES.keys())


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("""
╔═══════════════════════════════════════════════════════════════════════════╗
║            🦖 RAPTOR Pipeline Recommender Module v2.2.0                   ║
╠═══════════════════════════════════════════════════════════════════════════╣
║                                                                           ║
║  ML-based pipeline recommendation for RNA-seq differential expression.   ║
║                                                                           ║
║  SUPPORTED PIPELINES:                                                     ║
║  ────────────────────                                                     ║
║  • DESeq2      - General purpose, conservative, handles batch effects    ║
║  • edgeR       - Fast, handles low counts and high dispersion            ║
║  • limma-voom  - Very fast, best for large samples, complex designs      ║
║  • Wilcoxon    - Non-parametric, best FDR control for n≥8               ║
║  • edgeR_robust - Handles outliers well                                  ║
║                                                                           ║
║  DECISION FACTORS:                                                        ║
║  ─────────────────                                                        ║
║  1. Sample size per group (most critical)                                ║
║  2. Biological coefficient of variation (BCV)                            ║
║  3. Presence of outliers                                                 ║
║  4. Low count gene proportion                                            ║
║  5. Design complexity                                                    ║
║  6. Batch effects                                                        ║
║                                                                           ║
║  USAGE:                                                                   ║
║  ──────                                                                   ║
║  from raptor import RNAseqDataProfiler, PipelineRecommender              ║
║                                                                           ║
║  profiler = RNAseqDataProfiler(counts, metadata)                         ║
║  profile = profiler.run_full_profile()                                   ║
║  recommender = PipelineRecommender(profile)                              ║
║  recommendation = recommender.get_recommendation()                       ║
║  print(recommendation.summary())                                         ║
║                                                                           ║
║  # R code for recommended pipeline                                       ║
║  print(recommendation.r_code_primary)                                    ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
    """)
