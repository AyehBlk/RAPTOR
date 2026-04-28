#!/usr/bin/env python3
"""
RAPTOR v2.2.2 - Module 10: Biomarker Discovery

Comprehensive biomarker discovery from RNA-seq differential expression data.
Provides multi-method feature selection, classification, gene panel optimization,
survival analysis, and biological annotation.

Workflow Position:
    M7: DE Import → DEResult
    M8: Parameter Optimization → OptimizationResult
    M9: Ensemble Analysis → EnsembleResult
    M10: Biomarker Discovery (THIS MODULE) → BiomarkerResult

Submodules:
==========
10A - Feature Selection & Ranking
    - DE-based filtering (bridge from M7/M8/M9)
    - LASSO / Elastic Net (embedded, penalized regression)
    - Boruta (wrapper, RF-based shadow features)
    - mRMR (filter, minimum redundancy maximum relevance)
    - Recursive Feature Elimination (wrapper, iterative)
    - SHAP-based ranking (model-agnostic interpretability)

10B - Classification & Prediction Models
    - Logistic Regression (regularized baseline)
    - Random Forest (non-linear, OOB estimates)
    - Support Vector Machine (linear/RBF kernels)
    - XGBoost / Gradient Boosting (high performance)
    - Nested cross-validation, LOOCV, bootstrap CI

10C - Gene Panel Optimization
    - Forward selection (greedy, incremental)
    - Backward elimination (greedy, decremental)
    - Stability selection (subsampled LASSO)
    - Panel size vs. performance curve

10D - Survival Analysis (optional: requires lifelines)
    - Cox proportional hazards with LASSO (CoxNet)
    - Kaplan-Meier validation of biomarker panels
    - Concordance index (C-index) evaluation

10E - Biological Annotation & Reporting
    - Gene annotation via MyGene.info
    - Pathway enrichment (Fisher's exact test)
    - Structured output and summary

Study Designs Supported:
    - Binary case-control (disease vs healthy)
    - Multi-class (tumor subtypes A/B/C)
    - Paired/longitudinal (before vs after)
    - Survival/prognostic (time-to-event)
    - Cross-cohort validation (discovery + replication)

Scientific References:
    [1]  Tibshirani, R. (1996). LASSO. JRSS-B, 58(1), 267-288.
    [2]  Zou & Hastie (2005). Elastic Net. JRSS-B, 67(2), 301-320.
    [3]  Kursa & Rudnicki (2010). Boruta. JOSS, 36(11), 1-13.
    [4]  Ding & Peng (2005). mRMR. IEEE TPAMI, 27(8), 1226-1238.
    [5]  Guyon et al. (2002). RFE with SVM. Machine Learning, 46, 389-422.
    [6]  Lundberg & Lee (2017). SHAP. NeurIPS 2017.
    [7]  Meinshausen & Buhlmann (2010). Stability selection. JRSS-B.
    [8]  Phan et al. (2025). BoMGene: Boruta-mRMR for gene expression. arXiv.
    [9]  Simon et al. (2011). CoxNet regularized Cox regression. JSS, 39(5).
    [10] Cox, D.R. (1972). Regression models and life-tables. JRSS-B.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.2
License: MIT
"""

import logging
import json
import pickle
import warnings
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple, Union
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

# =============================================================================
# RAPTOR IMPORTS (with fallback)
# =============================================================================

try:
    from raptor.utils.validation import (
        validate_count_matrix,
        validate_metadata,
        validate_group_column,
        validate_file_path,
        validate_directory_path,
    )
    from raptor.utils.errors import (
        RAPTORError,
        ValidationError,
        DependencyError,
        handle_errors,
    )
    _RAPTOR_UTILS_AVAILABLE = True
except ImportError:
    warnings.warn("RAPTOR utils not found. Using fallback validation.")
    _RAPTOR_UTILS_AVAILABLE = False

    class RAPTORError(Exception):
        pass

    class ValidationError(Exception):
        def __init__(self, parameter: str, message: str = "", hint: str = ""):
            self.parameter = parameter
            self.message = message
            self.hint = hint
            super().__init__(f"{parameter}: {message}. {hint}")

    class DependencyError(ImportError):
        pass

    def handle_errors(func):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logging.error(f"Error in {func.__name__}: {e}")
                raise
        return wrapper

    def validate_count_matrix(df, **kw):
        if not isinstance(df, pd.DataFrame):
            raise ValidationError('counts', 'Must be a DataFrame')
        return df

    def validate_metadata(df, **kw):
        if not isinstance(df, pd.DataFrame):
            raise ValidationError('metadata', 'Must be a DataFrame')
        return df

    def validate_group_column(df, col, **kw):
        if col not in df.columns:
            raise ValidationError('group_column', f"'{col}' not found in metadata")
        return col

    def validate_file_path(p, **kw):
        return Path(p)

    def validate_directory_path(p, **kw):
        Path(p).mkdir(parents=True, exist_ok=True)
        return Path(p)


logger = logging.getLogger(__name__)

# RAPTOR version
__raptor_version__ = '2.2.2'


# =============================================================================
# OPTIONAL DEPENDENCY CHECKS
# =============================================================================

def _check_sklearn():
    try:
        import sklearn  # noqa: F401
        return True
    except ImportError:
        return False


def _check_xgboost():
    try:
        import xgboost  # noqa: F401
        return True
    except ImportError:
        return False


def _check_shap():
    try:
        import shap  # noqa: F401
        return True
    except ImportError:
        return False


def _check_boruta():
    try:
        from boruta import BorutaPy  # noqa: F401
        return True
    except ImportError:
        return False


def _check_mrmr():
    try:
        import mrmr  # noqa: F401
        return True
    except ImportError:
        return False


def _check_lifelines():
    try:
        import lifelines  # noqa: F401
        return True
    except ImportError:
        return False


def _check_pywgcna():
    try:
        import PyWGCNA  # noqa: F401
        return True
    except ImportError:
        return False


_SKLEARN_AVAILABLE = _check_sklearn()
_XGBOOST_AVAILABLE = _check_xgboost()
_SHAP_AVAILABLE = _check_shap()
_BORUTA_AVAILABLE = _check_boruta()
_MRMR_AVAILABLE = _check_mrmr()
_LIFELINES_AVAILABLE = _check_lifelines()
_PYWGCNA_AVAILABLE = _check_pywgcna()


# =============================================================================
# CONSTANTS
# =============================================================================

# Feature selection methods
FEATURE_SELECTION_METHODS = [
    'de_filter',       # From M7/M8/M9 DE results
    'lasso',           # L1-penalized logistic regression
    'elastic_net',     # L1+L2 penalized logistic regression
    'boruta',          # Random Forest shadow feature selection
    'mrmr',            # Minimum Redundancy Maximum Relevance
    'rfe',             # Recursive Feature Elimination
    'shap',            # SHAP-based importance ranking
    'wgcna',           # WGCNA hub genes from co-expression modules
]

# Classification models
CLASSIFIER_NAMES = [
    'logistic_regression',
    'random_forest',
    'svm',
    'xgboost',
]

# Study design types
STUDY_DESIGNS = [
    'binary',        # Two-group comparison
    'multiclass',    # Multiple groups
    'paired',        # Paired/longitudinal
    'survival',      # Time-to-event
]

# Validation strategies
VALIDATION_STRATEGIES = [
    'nested_cv',     # Nested cross-validation (gold standard)
    'loocv',         # Leave-one-out CV (small samples)
    'split',         # Discovery/validation split
    'bootstrap',     # Bootstrap confidence intervals
]

# Default parameters
DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_LFC_THRESHOLD = 0.0
DEFAULT_N_FOLDS_OUTER = 5
DEFAULT_N_FOLDS_INNER = 3
DEFAULT_N_BOOTSTRAP = 100
DEFAULT_PANEL_MIN = 3
DEFAULT_PANEL_MAX = 50
DEFAULT_RANDOM_STATE = 42


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class FeatureSelectionResult:
    """
    Result from a single feature selection method.

    Attributes
    ----------
    method : str
        Name of the feature selection method.
    selected_genes : List[str]
        Genes selected by this method.
    gene_scores : pd.DataFrame
        DataFrame with gene_id index and 'score' column (importance/rank).
    n_selected : int
        Number of genes selected.
    parameters : Dict
        Parameters used for this method.
    """
    method: str
    selected_genes: List[str]
    gene_scores: pd.DataFrame
    n_selected: int
    parameters: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ClassificationResult:
    """
    Result from classifier evaluation.

    Attributes
    ----------
    model_name : str
        Name of the classifier.
    accuracy : float
        Mean accuracy across CV folds.
    auc : float
        Mean AUC-ROC across CV folds.
    sensitivity : float
        Mean sensitivity (recall for positive class).
    specificity : float
        Mean specificity (recall for negative class).
    f1 : float
        Mean F1-score.
    metrics_per_fold : List[Dict]
        Per-fold performance metrics.
    confusion_matrix : Optional[np.ndarray]
        Aggregated confusion matrix.
    roc_data : Optional[Dict]
        ROC curve data (fpr, tpr, thresholds).
    feature_importance : Optional[pd.DataFrame]
        Feature importance from the trained model.
    trained_model : Optional[Any]
        The fitted model object (from final full-data fit).
    oof_true : Optional[np.ndarray]
        Out-of-fold ground-truth labels (M1). One entry per sample,
        populated by ``evaluate_nested_cv`` (one OOF prediction per
        outer fold) and ``evaluate_loocv`` (one OOF prediction per
        left-out sample). These are the honest, held-out labels that
        correspond to ``oof_prob``. Downstream clinical metrics should
        use ``(oof_true, oof_prob)`` rather than predictions from the
        full-data ``trained_model`` to avoid reporting apparent
        (training-data) performance as if it were generalisation
        performance.
    oof_prob : Optional[np.ndarray]
        Out-of-fold predicted probabilities for the positive class
        (M1). Same length and order as ``oof_true``. Each entry is a
        prediction from a model that did NOT see the corresponding
        sample during training.
    """
    model_name: str
    accuracy: float = 0.0
    auc: float = 0.0
    sensitivity: float = 0.0
    specificity: float = 0.0
    f1: float = 0.0
    metrics_per_fold: List[Dict] = field(default_factory=list)
    confusion_matrix: Optional[np.ndarray] = None
    roc_data: Optional[Dict] = None
    feature_importance: Optional[pd.DataFrame] = None
    trained_model: Optional[Any] = None
    oof_true: Optional[np.ndarray] = None
    oof_prob: Optional[np.ndarray] = None


# -----------------------------------------------------------------------
# Deterministic "best classifier" selection with auditable tiebreak (M5).
# -----------------------------------------------------------------------

# Preference order when classifiers tie on both AUC and F1.
# Most interpretable first: LR > SVM > RF > XGBoost. The goal is to
# prefer simpler / more transparent models when performance is
# indistinguishable within numerical tolerance.
_CLASSIFIER_PREFERENCE_ORDER: Tuple[str, ...] = (
    'logistic_regression',
    'svm',
    'random_forest',
    'xgboost',
)


def _select_best_classifier(
    clf_results: Dict[str, ClassificationResult],
    auc_epsilon: float = 1e-6,
    f1_epsilon: float = 1e-6,
) -> str:
    """
    Pick the best classifier with a deterministic, auditable tiebreak.

    Selection rules, applied in order:

    1. Highest AUC.
    2. If two or more classifiers tie on AUC within ``auc_epsilon``,
       prefer the one with the highest F1.
    3. If still tied on F1 within ``f1_epsilon``, prefer the classifier
       highest in ``_CLASSIFIER_PREFERENCE_ORDER`` (LR > SVM > RF >
       XGBoost). Classifier names not in the preference tuple are
       ranked after all listed ones; further ties are broken
       alphabetically so the result is deterministic.

    When a tiebreak at step 2 or 3 actually fires, the decision is
    logged at INFO level so the choice is reproducible from the run
    log. Step 1 (single clear winner) is silent.

    Parameters
    ----------
    clf_results : dict of str -> ClassificationResult
        Mapping from classifier name to its evaluation result.
    auc_epsilon : float, default 1e-6
        AUC values within this distance of the max count as tied.
    f1_epsilon : float, default 1e-6
        F1 values within this distance of the max count as tied.

    Returns
    -------
    str
        The selected classifier key from ``clf_results``.

    Raises
    ------
    ValueError
        If ``clf_results`` is empty.
    """
    if not clf_results:
        raise ValueError(
            "clf_results is empty; cannot select a best classifier."
        )

    names = list(clf_results.keys())

    # Step 1: highest AUC (with epsilon tolerance for float ties).
    max_auc = max(clf_results[n].auc for n in names)
    auc_tied = [
        n for n in names
        if max_auc - clf_results[n].auc <= auc_epsilon
    ]

    if len(auc_tied) == 1:
        return auc_tied[0]

    # Step 2: tiebreak on F1.
    max_f1 = max(clf_results[n].f1 for n in auc_tied)
    f1_tied = [
        n for n in auc_tied
        if max_f1 - clf_results[n].f1 <= f1_epsilon
    ]

    if len(f1_tied) == 1:
        winner = f1_tied[0]
        logger.info(
            f"   Classifier tie on AUC={max_auc:.6f} among {auc_tied}; "
            f"resolved by F1 (winner: {winner}, F1={max_f1:.6f})."
        )
        return winner

    # Step 3: tiebreak on interpretability preference order.
    def _pref_key(name: str) -> Tuple[int, str]:
        if name in _CLASSIFIER_PREFERENCE_ORDER:
            return (_CLASSIFIER_PREFERENCE_ORDER.index(name), name)
        # Unknown classifier: push after all listed ones, alphabetical.
        return (len(_CLASSIFIER_PREFERENCE_ORDER), name)

    winner = min(f1_tied, key=_pref_key)
    logger.info(
        f"   Classifier tie on AUC={max_auc:.6f} and F1={max_f1:.6f} "
        f"among {f1_tied}; resolved by interpretability preference "
        f"(winner: {winner})."
    )
    return winner


# -----------------------------------------------------------------------
# M6: Nogueira et al. (2018 JMLR) feature-selection stability measure.
# -----------------------------------------------------------------------
#
# References
# ----------
# Nogueira, Sechidis & Brown (2018). "On the Stability of Feature
#   Selection Algorithms." Journal of Machine Learning Research 18:174,
#   pp. 1-54. https://jmlr.org/papers/v18/17-514.html
#
# The scalar stability measure is:
#
#     Phi_hat = 1 - ( (1/p) * sum_j [ M/(M-1) * (h_j/M) * (1 - h_j/M) ] )
#                   / ( (q/(Mp)) * (1 - q/(Mp)) )
#
# where:
#     p   = total number of features (candidate genes)
#     M   = number of feature subsets (one per fold; possibly >1 per
#           fold if repeated CV)
#     h_j = number of subsets containing gene j
#     q   = sum_j h_j = total "slots" used across all subsets
#
# Interpretation:
#     numerator   = average per-gene variance of selection indicators
#     denominator = the variance that would be observed under random
#                   uniform selection of q/M features from p total
#     Phi_hat     = 1 - (observed variance) / (chance variance)
#                 = 1 when selections are deterministic (identical folds)
#                 = 0 on average under random selection
#                 < 0 when selections are more dissimilar than chance
#
# Benchmark scale (Nogueira 2018, reiterated by stabm R package docs):
#     Phi_hat >= 0.75  excellent agreement beyond chance
#     0.4 <= Phi_hat < 0.75  intermediate to good
#     Phi_hat < 0.4    poor agreement

_NOGUEIRA_EXCELLENT_THRESHOLD = 0.75
_NOGUEIRA_POOR_THRESHOLD = 0.4


def _nogueira_stability_from_matrix(Z: np.ndarray) -> float:
    """
    Compute Nogueira et al. (2018) stability from an M x p binary
    selection matrix.

    Parameters
    ----------
    Z : np.ndarray of shape (M, p)
        Z[i, j] = 1 iff feature j was in subset i.

    Returns
    -------
    float
        Phi in [-1, 1]. Returns 0.0 in the degenerate cases where no
        gene is selected in any subset or where every gene is selected
        in every subset (both of which make the chance-variance term 0).
    """
    Z = np.asarray(Z, dtype=int)
    if Z.ndim != 2:
        raise ValueError(
            f"Expected 2-D selection matrix (M x p), got shape {Z.shape}."
        )
    M, p = Z.shape
    if M < 2:
        # Stability is not defined for a single subset; convention: 1.0
        # if the single subset is non-empty, else 0.0.
        return 1.0 if Z.sum() > 0 else 0.0
    if p == 0:
        return 0.0

    h = Z.sum(axis=0)                       # shape (p,), selection counts per gene
    q = int(h.sum())                        # total selections across all subsets

    # Chance-variance term: variance under random uniform selection of
    # k_bar = q/M features from p total (Bernoulli with p_sel = q/(Mp)).
    p_sel = q / (M * p)
    chance_var = p_sel * (1.0 - p_sel)
    if chance_var <= 0.0:
        # Everything selected, or nothing selected: the formula is
        # undefined. Convention: return 0.0 (no information).
        return 0.0

    # Observed per-gene variance, unbiased estimator with Bessel's correction.
    p_j = h / M                             # selection probability estimate per gene
    per_gene_var = (M / (M - 1.0)) * p_j * (1.0 - p_j)
    observed_var = float(per_gene_var.mean())

    return float(1.0 - observed_var / chance_var)


def _panels_to_selection_matrix(
    per_fold_panels: List[List[str]],
    gene_universe: List[str],
) -> Tuple[np.ndarray, pd.Series]:
    """
    Convert a list of per-fold gene panels to the M x p binary selection
    matrix used by the Nogueira stability formula, plus the per-gene
    selection frequency series.

    Parameters
    ----------
    per_fold_panels : list of lists of gene names
        One entry per fold (or fold-repeat).
    gene_universe : list of gene names
        The full set of candidate genes. Any gene in a fold panel that
        is NOT in gene_universe is silently ignored (this shouldn't
        happen for well-formed inputs, but the Nogueira formula needs a
        fixed column space).

    Returns
    -------
    Z : np.ndarray, shape (len(per_fold_panels), len(gene_universe))
    freq : pd.Series indexed by gene_universe, values h_j / M
    """
    if not per_fold_panels:
        return np.zeros((0, len(gene_universe)), dtype=int), pd.Series(
            0.0, index=gene_universe, name='selection_frequency'
        )

    gene_to_col = {g: i for i, g in enumerate(gene_universe)}
    M = len(per_fold_panels)
    p = len(gene_universe)
    Z = np.zeros((M, p), dtype=int)
    for i, panel in enumerate(per_fold_panels):
        for g in panel:
            col = gene_to_col.get(g)
            if col is not None:
                Z[i, col] = 1

    freq = pd.Series(
        Z.mean(axis=0), index=gene_universe, name='selection_frequency'
    )
    return Z, freq


def _bootstrap_nogueira_ci(
    Z: np.ndarray,
    n_bootstrap: int = 1000,
    ci: float = 0.95,
    seed: int = 42,
) -> Tuple[float, float]:
    """
    Bootstrap CI for the Nogueira stability over the fold axis.

    We resample M folds with replacement from the rows of Z and recompute
    stability on each resample. The (alpha/2, 1 - alpha/2) percentiles
    give the CI.

    Caveat: with small M (e.g. 5 folds), there are only 252 distinct
    bootstrap resamples of 5 items each, so the CI is intentionally
    wide. This honestly reflects the small-M uncertainty — use
    ``n_repeats > 1`` to tighten it.

    Parameters
    ----------
    Z : np.ndarray, shape (M, p)
        Binary selection matrix.
    n_bootstrap : int, default 1000
        Number of bootstrap resamples.
    ci : float, default 0.95
        Confidence level.
    seed : int, default 42
        RNG seed.

    Returns
    -------
    (lower, upper) : Tuple[float, float]
    """
    M = Z.shape[0]
    if M < 2:
        return (0.0, 0.0)

    rng = np.random.default_rng(seed)
    estimates: List[float] = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, M, size=M)
        Z_boot = Z[idx, :]
        estimates.append(_nogueira_stability_from_matrix(Z_boot))

    alpha = 1.0 - ci
    lo = float(np.percentile(estimates, 100.0 * alpha / 2.0))
    hi = float(np.percentile(estimates, 100.0 * (1.0 - alpha / 2.0)))
    return (lo, hi)


def _nogueira_benchmark_label(stability: float) -> str:
    """Benchmark-scale label per Nogueira et al. (2018)."""
    if stability >= _NOGUEIRA_EXCELLENT_THRESHOLD:
        return 'excellent'
    if stability < _NOGUEIRA_POOR_THRESHOLD:
        return 'poor'
    return 'intermediate'


def _compute_panel_stability(
    per_fold_panels: List[List[str]],
    per_fold_ranked_genes: List[pd.DataFrame],
    gene_universe: List[str],
    final_panel: List[str],
    n_folds: int,
    n_repeats: int,
    random_state: int = DEFAULT_RANDOM_STATE,
) -> 'PanelStabilityResult':
    """
    Build a ``PanelStabilityResult`` from per-fold panels and the final
    all-data panel.

    This is the single function that populates every field of the
    stability result. It is called once at the end of the pipeline-CV
    run in ``discover_biomarkers``.
    """
    Z, freq = _panels_to_selection_matrix(per_fold_panels, gene_universe)

    if Z.shape[0] < 2:
        # Degenerate: only one fold, cannot estimate stability.
        return PanelStabilityResult(
            per_fold_panels=per_fold_panels,
            per_fold_ranked_genes=per_fold_ranked_genes,
            gene_selection_frequency=freq,
            nogueira_stability=0.0,
            nogueira_stability_ci=(0.0, 0.0),
            n_folds=n_folds,
            n_repeats=n_repeats,
            benchmark_label='poor',
            final_panel_overlap=pd.Series(0.0, index=final_panel or []),
        )

    phi = _nogueira_stability_from_matrix(Z)
    ci = _bootstrap_nogueira_ci(Z, n_bootstrap=1000, ci=0.95, seed=random_state)
    label = _nogueira_benchmark_label(phi)

    # Per-gene overlap for the final panel specifically: how often did
    # each final-panel gene also appear in a fold panel?
    final_overlap = pd.Series(
        {g: float(freq.get(g, 0.0)) for g in (final_panel or [])},
        name='final_panel_overlap',
    )

    return PanelStabilityResult(
        per_fold_panels=per_fold_panels,
        per_fold_ranked_genes=per_fold_ranked_genes,
        gene_selection_frequency=freq,
        nogueira_stability=phi,
        nogueira_stability_ci=ci,
        n_folds=n_folds,
        n_repeats=n_repeats,
        benchmark_label=label,
        final_panel_overlap=final_overlap,
    )


# -----------------------------------------------------------------------
# M4: Significance-calibrated consensus score.
# -----------------------------------------------------------------------
#
# References
# ----------
# Kolde et al. (2012). "Robust rank aggregation for gene list integration
#   and meta-analysis." Bioinformatics 28(4):573-580. The RRA algorithm
#   computes a p-value for each gene's mean rank against a uniform-null.
#   M4 is a computationally trivial approximation: use the per-gene DE
#   p-value as a proxy for "is this gene distinguishable from chance at
#   all," and apply a two-tier multiplicative shrinkage.
#
# Efron (2004). "Large-scale simultaneous hypothesis testing: the choice
#   of a null hypothesis." JASA 99(465):96-104. Empirical-Bayes local
#   FDR shrinks non-significant effects by a data-driven factor. M4's
#   0.5 is a fixed simplification -- one knob, no mixture model.
#
# Conceptual shape
# ----------------
# For each gene:
#     p_value = Welch's t-test (or Mann-Whitney) between the two groups
#     weight  = weight_significant    if p_value <  alpha
#             = weight_nonsignificant otherwise
#     consensus_score_calibrated = consensus_score * weight
# The final consensus_rank is recomputed on consensus_score_calibrated
# so the dashboard-facing ordering reflects both method agreement AND
# a significance floor. Uncalibrated consensus_score is preserved for
# auditability.
#
# Scope: this runs ONCE on the final all-data consensus ranking in
# discover_biomarkers. It does NOT run inside M6's per-fold pipeline-CV
# loop -- the per-fold rankings drive panel selection, not reporting,
# so calibrating them would cost N_folds * DE computations for no
# statistical benefit.

def apply_significance_calibration(
    ranked_genes: pd.DataFrame,
    X: pd.DataFrame,
    y: np.ndarray,
    reference_group: Optional[str] = None,
    baseline_group: Optional[str] = None,
    alpha: float = 0.05,
    weight_significant: float = 1.0,
    weight_nonsignificant: float = 0.5,
    test: str = "welch",
) -> pd.DataFrame:
    """
    Add significance-calibration columns to a consensus-ranking DataFrame.

    For every gene in ``ranked_genes``, computes a two-sample DE test on
    the provided expression matrix, assigns a weight based on whether
    the p-value falls below ``alpha``, and rescales ``consensus_score``
    multiplicatively. The final ``consensus_rank`` column is recomputed
    on the calibrated score so downstream code that sorts by rank
    reflects the significance floor automatically.

    This is a post-processing step on a ranked-genes table; it does not
    touch the underlying FeatureSelector results.

    Parameters
    ----------
    ranked_genes : pd.DataFrame
        Output of ``FeatureSelector.consensus_ranking``. Must have at
        least the columns ``consensus_score`` and ``consensus_rank``;
        all other columns are preserved unchanged.
        Index is gene_id.
    X : pd.DataFrame (samples x genes)
        Expression matrix. Rows align with ``y``. Columns are gene_ids
        matching ``ranked_genes.index``. Any gene in ``ranked_genes``
        that is not in ``X.columns`` gets p_value=1.0, weight=
        weight_nonsignificant (honest: we cannot test it).
    y : np.ndarray of shape (n_samples,)
        Binary labels (0/1). Must contain exactly 2 unique classes;
        otherwise this function is a no-op and returns ``ranked_genes``
        unchanged with a logged warning (diagnostic/translational
        intents only, per scoping doc §6.5).
    reference_group : str, optional
        Label string for the 'disease' / 'reference' group. If None,
        inferred as the encoding of y==1.
    baseline_group : str, optional
        Label string for the control group. If None, inferred as y==0.
    alpha : float, default 0.05
        P-value threshold defining the significance tier.
    weight_significant : float, default 1.0
        Multiplier applied to ``consensus_score`` when p_value < alpha.
    weight_nonsignificant : float, default 0.5
        Multiplier applied to ``consensus_score`` when p_value >= alpha.
        Set to 1.0 to disable calibration (every gene keeps its raw
        score) while still reporting p-values.
    test : {'welch', 'mann_whitney'}, default 'welch'
        Two-sample test. Welch is vectorized and fast (~10ms for 20k
        genes); Mann-Whitney is looped (~2s) but more robust for
        small-n or heavy-outlier data.

    Returns
    -------
    pd.DataFrame
        Copy of ``ranked_genes`` with three new columns:
            p_value                    : float in [0, 1]
            weight                     : float (weight_significant or weight_nonsignificant)
            consensus_score_calibrated : float (= consensus_score * weight)
        And ``consensus_rank`` recomputed on ``consensus_score_calibrated``.
        Original ``consensus_score`` column is preserved unchanged for
        auditability.
    """
    from raptor.biomarker_discovery.univariate_de import compute_per_gene_de

    # --- Pre-validate inputs that, if wrong, mean we can't calibrate ---
    if 'consensus_score' not in ranked_genes.columns:
        raise ValueError(
            "ranked_genes must have a 'consensus_score' column; "
            "is it the output of FeatureSelector.consensus_ranking?"
        )

    y = np.asarray(y)
    unique_y = np.unique(y)
    if len(unique_y) != 2:
        # Non-binary labels: return unchanged with a logged note. Matches
        # DirectionPattern's convention; supports the prognostic/
        # monitoring/exploratory intents without forcing them through
        # a binary DE test that isn't well-defined.
        logger.warning(
            f"   apply_significance_calibration: y has {len(unique_y)} "
            f"unique classes; binary DE not defined. ranked_genes "
            f"returned unchanged."
        )
        return ranked_genes.copy()

    # --- Resolve reference/baseline group names and build label array ---
    # Two cases:
    #   (a) y is already string labels (e.g. from direction_patterns or
    #       an upstream wrapper): use y directly; reference/baseline
    #       group names must match y's unique values.
    #   (b) y is 0/1 encoding (the _prepare_expression_data convention):
    #       synthesize string labels from the provided group names.
    y_dtype_kind = np.asarray(y).dtype.kind
    if y_dtype_kind in ('U', 'S', 'O'):
        # String / object y: use as-is
        labels = np.asarray(y)
        unique_str = set(labels.tolist())
        if reference_group is None or baseline_group is None:
            # Infer: sorted gives a deterministic pick.
            srt = sorted(unique_str)
            baseline_group = baseline_group or srt[0]
            reference_group = reference_group or srt[1]
        if reference_group not in unique_str or baseline_group not in unique_str:
            logger.warning(
                f"   apply_significance_calibration: reference_group "
                f"{reference_group!r} or baseline_group "
                f"{baseline_group!r} not in y's unique values "
                f"{sorted(unique_str)}. ranked_genes returned unchanged."
            )
            return ranked_genes.copy()
    else:
        # Numeric y: synthesize string labels from group names (or
        # fallbacks) and encode 1 -> reference, 0 -> baseline.
        if reference_group is None:
            reference_group = 'class_1'
        if baseline_group is None:
            baseline_group = 'class_0'
        labels = np.where(np.asarray(y) == 1, reference_group, baseline_group)

    # --- Align X columns to ranked_genes index ---
    # We tolerate a mismatch (X may have been subset to the final
    # panel or filtered differently from the ranking gene list). Genes
    # in ranked_genes.index but missing from X.columns get p=1.0.
    common_genes = [g for g in ranked_genes.index if g in X.columns]
    missing = [g for g in ranked_genes.index if g not in X.columns]
    if missing and len(missing) > 10:
        logger.info(
            f"   apply_significance_calibration: {len(missing)} of "
            f"{len(ranked_genes)} ranked genes not in expression matrix; "
            f"they get p=1.0 (weight={weight_nonsignificant})."
        )
    elif missing:
        logger.debug(
            f"   apply_significance_calibration: genes missing from X: "
            f"{missing}"
        )

    # --- Compute per-gene DE on the tested subset ---
    if common_genes:
        X_common = X[common_genes]
        try:
            de_table = compute_per_gene_de(
                expression=X_common,
                labels=labels,
                reference_group=reference_group,
                baseline_group=baseline_group,
                test=test,
            )
        except Exception as e:
            logger.warning(
                f"   apply_significance_calibration: DE computation "
                f"failed ({e}); ranked_genes returned unchanged."
            )
            return ranked_genes.copy()
    else:
        # No genes to test: everything gets weight_nonsignificant.
        de_table = pd.DataFrame(
            columns=['p_value', 'neg_log10p', 'log2fc',
                     'mean_ref', 'mean_base', 'n_ref', 'n_base', 'skipped'],
        )
        de_table.index.name = 'gene_id'

    # --- Build the output DataFrame ---
    out = ranked_genes.copy()

    # p_value: default to 1.0 for genes not in de_table (missing from X)
    p_values = pd.Series(1.0, index=out.index, name='p_value', dtype=float)
    if not de_table.empty:
        common_idx = de_table.index.intersection(out.index)
        p_values.loc[common_idx] = de_table.loc[common_idx, 'p_value'].values
    out['p_value'] = p_values

    # weight: two-tier multiplicative
    out['weight'] = np.where(
        out['p_value'] < alpha,
        weight_significant,
        weight_nonsignificant,
    )

    # consensus_score_calibrated = consensus_score * weight
    out['consensus_score_calibrated'] = (
        out['consensus_score'] * out['weight']
    )

    # Recompute consensus_rank based on the calibrated score. Use
    # method='first' for a deterministic tie-break matching the
    # original consensus_ranking convention.
    out['consensus_rank'] = out['consensus_score_calibrated'].rank(
        ascending=False, method='first',
    ).astype(int)

    # Re-sort by new rank so downstream head(N) pulls the top by
    # calibrated score.
    out = out.sort_values('consensus_rank')

    n_sig = int((out['p_value'] < alpha).sum())
    n_nonsig = len(out) - n_sig
    logger.info(
        f"   M4 significance calibration: {n_sig} genes passed "
        f"alpha={alpha}; {n_nonsig} genes down-weighted by "
        f"{weight_nonsignificant}x (test={test})."
    )

    return out


@dataclass
class PanelOptimizationResult:
    """
    Result from gene panel size optimization.

    Attributes
    ----------
    optimal_panel : List[str]
        Genes in the optimal panel.
    optimal_size : int
        Optimal number of genes.
    optimal_auc : float
        AUC at optimal panel size.
    panel_curve : pd.DataFrame
        DataFrame with columns: panel_size, auc_mean, auc_std.
    all_panels : Dict[int, List[str]]
        Gene lists at each tested panel size.
    method : str
        Panel optimization method used.
    """
    optimal_panel: List[str]
    optimal_size: int
    optimal_auc: float
    panel_curve: pd.DataFrame
    all_panels: Dict[int, List[str]] = field(default_factory=dict)
    method: str = 'forward_selection'
    # selection_method records HOW optimal_size was chosen, separate
    # from `method` (which records WHICH panel-construction algorithm
    # ran, e.g. 'forward_selection' vs 'stability_selection'). Values:
    # 'kneedle', 'argmax_fallback', 'argmax', 'first_drop',
    # 'user_specified'. Surfaced on the dashboard so users can tell
    # whether the picked size came from a true knee, a saturated-curve
    # fallback, or an explicit target_panel_size override.
    selection_method: str = 'kneedle'


# -----------------------------------------------------------------------
# M6: Panel stability across pipeline-CV folds.
# -----------------------------------------------------------------------
#
# M6 pushes feature selection and panel optimization inside the outer CV
# loop (Ambroise & McLachlan 2002 PNAS 99:6562; Varma & Simon 2006 BMC
# Bioinf 7:91). A natural byproduct is one panel per fold. The degree to
# which these per-fold panels agree is itself a signal: stable panels
# indicate a robust biomarker set, unstable panels indicate the chosen
# genes are sensitive to small data perturbations.
#
# We quantify stability with the measure of Nogueira, Sechidis & Brown
# (JMLR 2018 18:174), which is the only similarity-based measure that
# satisfies all five of their axiomatic properties: fully defined, strict
# monotonicity, bounded, maximum stability iff deterministic selection,
# and correction for chance. Jaccard, Dice, and POG all fail
# correction-for-chance and systematically reward larger feature sets.
#
# The benchmark scale from the same paper (also used by stabm R package):
#     >= 0.75  excellent agreement beyond chance
#     0.4-0.75 intermediate to good
#     <  0.4   poor agreement

@dataclass
class PanelStabilityResult:
    """
    Cross-fold panel stability diagnostic from pipeline-CV (M6).

    Populated when ``discover_biomarkers`` runs with feature selection
    inside the outer CV loop. Reports how consistent per-fold biomarker
    panels are across folds, using the Nogueira et al. (2018 JMLR)
    stability measure — the only similarity-based measure that satisfies
    all five desirable axiomatic properties including correction for
    chance (which Jaccard, Dice, and POG lack).

    Attributes
    ----------
    per_fold_panels : List[List[str]]
        One gene list per outer CV fold (and per repeat, if n_repeats > 1).
        Length = n_outer * n_repeats. Order is (repeat_0_fold_0,
        repeat_0_fold_1, ..., repeat_r_fold_f).
    per_fold_ranked_genes : List[pd.DataFrame]
        Consensus ranking produced in each fold (same order as
        per_fold_panels). Useful for dashboard drill-down to inspect why
        a fold picked different genes.
    gene_selection_frequency : pd.Series
        For each gene in the universe, the fraction of folds in which it
        appeared in the fold panel. Values in [0, 1]. Genes not in any
        fold panel have frequency 0. A gene selected in every fold has
        frequency 1.
    nogueira_stability : float
        Scalar stability measure in [-1, 1]. Expected 0 under random
        selection, 1 for deterministic selection. Formula:
            Phi = 1 - (avg per-gene selection variance) / (chance variance)
    nogueira_stability_ci : Tuple[float, float]
        95% confidence interval for the stability measure, from a
        bootstrap over folds (1000 resamples). Note: for small M (e.g.
        5 folds) the CI is intentionally wide — this honestly reflects
        the small-sample uncertainty of the estimate.
    n_folds : int
        Number of outer folds per repeat.
    n_repeats : int
        Number of outer-CV repeats.
    benchmark_label : str
        One of {'excellent', 'intermediate', 'poor'} per the Nogueira 2018
        benchmark scale. 'excellent' = stability >= 0.75, 'poor' < 0.4,
        'intermediate' otherwise.
    final_panel_overlap : pd.Series
        For each gene in the final (all-data) panel, the fraction of folds
        in which that gene also appeared. 1.0 means the gene was picked
        in every fold; 0.0 means no fold picked it (which can happen when
        the final panel uses information from the full dataset that no
        single fold could see).
    """
    per_fold_panels: List[List[str]] = field(default_factory=list)
    per_fold_ranked_genes: List[pd.DataFrame] = field(default_factory=list)
    gene_selection_frequency: pd.Series = field(
        default_factory=lambda: pd.Series(dtype=float)
    )
    nogueira_stability: float = 0.0
    nogueira_stability_ci: Tuple[float, float] = (0.0, 0.0)
    n_folds: int = 0
    n_repeats: int = 1
    benchmark_label: str = 'poor'
    final_panel_overlap: pd.Series = field(
        default_factory=lambda: pd.Series(dtype=float)
    )

    def summary(self) -> str:
        """Human-readable one-paragraph summary."""
        lo, hi = self.nogueira_stability_ci
        return (
            f"Panel stability (Nogueira 2018): "
            f"Phi = {self.nogueira_stability:.3f} "
            f"[95% CI: {lo:.3f}, {hi:.3f}] "
            f"({self.benchmark_label}) "
            f"across {self.n_folds} folds x {self.n_repeats} repeats."
        )


@dataclass
class SurvivalResult:
    """
    Result from survival analysis.

    Attributes
    ----------
    significant_genes : List[str]
        Genes significantly associated with survival.
    cox_results : pd.DataFrame
        Cox regression results (HR, p-value, CI).
    c_index : float
        Concordance index for the panel.
    panel_genes : List[str]
        Genes in the survival panel (CoxNet selected).
    km_groups : Optional[Dict]
        Kaplan-Meier group assignments.
    """
    significant_genes: List[str] = field(default_factory=list)
    cox_results: Optional[pd.DataFrame] = None
    c_index: float = 0.0
    panel_genes: List[str] = field(default_factory=list)
    km_groups: Optional[Dict] = None


@dataclass
class AnnotationResult:
    """
    Biological annotation and enrichment results.

    Contains gene annotations, pathway enrichment, literature associations,
    and protein-protein interaction context for the biomarker panel.

    Attributes
    ----------
    gene_annotations : pd.DataFrame
        Per-gene annotations: symbol, name, GO terms, pathways.
    pathway_enrichment : pd.DataFrame
        Enriched pathways with p-values: pathway, source, p_value, padj,
        overlap_genes, overlap_count, gene_set_size, odds_ratio.
    literature_hits : pd.DataFrame
        PubMed/Europe PMC hits: gene_id, n_publications, top_terms,
        disease_associations, recent_pmids.
    ppi_network : Optional[Dict]
        STRING PPI network data: nodes, edges, enrichment_pvalue.
    background_size : int
        Number of background genes used for enrichment.
    species : str
        Species used for annotation.
    """
    gene_annotations: pd.DataFrame = field(default_factory=pd.DataFrame)
    pathway_enrichment: pd.DataFrame = field(default_factory=pd.DataFrame)
    literature_hits: pd.DataFrame = field(default_factory=pd.DataFrame)
    ppi_network: Optional[Dict] = None
    background_size: int = 0
    species: str = 'human'

    @property
    def n_annotated(self) -> int:
        return len(self.gene_annotations)

    @property
    def n_enriched_pathways(self) -> int:
        if self.pathway_enrichment.empty:
            return 0
        if 'adjusted_p_value' in self.pathway_enrichment.columns:
            return int((self.pathway_enrichment['adjusted_p_value'] < 0.05).sum())
        return len(self.pathway_enrichment)

    def save(self, output_dir: Path):
        """Save all annotation outputs."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not self.gene_annotations.empty:
            self.gene_annotations.to_csv(
                output_dir / "gene_annotations.csv", encoding='utf-8'
            )
        if not self.pathway_enrichment.empty:
            self.pathway_enrichment.to_csv(
                output_dir / "pathway_enrichment.csv", index=False, encoding='utf-8'
            )
        if not self.literature_hits.empty:
            self.literature_hits.to_csv(
                output_dir / "literature_hits.csv", index=False, encoding='utf-8'
            )
        if self.ppi_network is not None:
            with open(output_dir / "ppi_network.json", 'w', encoding='utf-8') as f:
                json.dump(self.ppi_network, f, indent=2, default=str)


@dataclass
class BiomarkerResult:
    """
    Comprehensive biomarker discovery result.

    Core data structure for RAPTOR Module 10. Contains all outputs from
    feature selection, classification, panel optimization, and annotation.

    Attributes
    ----------
    ranked_genes : pd.DataFrame
        All candidate genes with multi-method consensus ranking.
        Columns: gene_id (index), consensus_rank, consensus_score,
        n_methods_selected, per-method scores.
    panel : List[str]
        Final recommended gene panel.
    panel_size : int
        Number of genes in the panel.

    selection_results : Dict[str, FeatureSelectionResult]
        Per-method feature selection results.
    classification_results : Dict[str, ClassificationResult]
        Per-classifier evaluation results.
    best_classifier : str
        Name of the best performing classifier.

    panel_optimization : Optional[PanelOptimizationResult]
        Panel size optimization results.
    survival_result : Optional[SurvivalResult]
        Survival analysis results (if applicable).

    annotations : Optional[AnnotationResult]
        Biological annotation results (gene info, pathways, literature, PPI).

    study_design : str
        Study design type used.
    validation_strategy : str
        Validation strategy used.
    n_samples : int
        Number of samples in the dataset.
    n_initial_candidates : int
        Number of initial candidate genes.

    parameters : Dict
        All parameters used in the analysis.
    metadata : Dict
        Additional metadata.
    timestamp : str
        ISO format timestamp.
    """
    # Core outputs
    ranked_genes: pd.DataFrame
    panel: List[str]
    panel_size: int

    # Detailed results
    selection_results: Dict[str, FeatureSelectionResult] = field(default_factory=dict)
    classification_results: Dict[str, ClassificationResult] = field(default_factory=dict)
    best_classifier: str = ''

    # Optional results
    panel_optimization: Optional[PanelOptimizationResult] = None
    survival_result: Optional[SurvivalResult] = None

    # Annotations
    annotation_result: Optional[AnnotationResult] = None

    # M6: Panel stability across pipeline-CV folds (Nogueira 2018 JMLR).
    # Populated when discover_biomarkers runs with validation='nested_cv'.
    # None when LOOCV is used (no fold-to-fold comparison is defined).
    panel_stability: Optional[PanelStabilityResult] = None

    # Context
    study_design: str = 'binary'
    validation_strategy: str = 'nested_cv'
    n_samples: int = 0
    n_initial_candidates: int = 0

    # Metadata
    parameters: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 70,
            "RAPTOR v2.2.2 - MODULE 10: BIOMARKER DISCOVERY RESULTS",
            "=" * 70,
            f"Study Design: {self.study_design}",
            f"Validation: {self.validation_strategy}",
            f"Samples: {self.n_samples}",
            f"Initial candidates: {self.n_initial_candidates}",
            "",
            "FEATURE SELECTION:",
        ]

        for method, res in self.selection_results.items():
            lines.append(f"  {method}: {res.n_selected} genes selected")

        lines.append(f"\nCONSENSUS RANKING: {len(self.ranked_genes)} genes ranked")

        # M4: Significance-calibration headline, when the calibrated column is present
        if (
            self.ranked_genes is not None
            and 'weight' in self.ranked_genes.columns
            and 'p_value' in self.ranked_genes.columns
        ):
            alpha = self.parameters.get('alpha', 0.05)
            n_sig = int((self.ranked_genes['p_value'] < alpha).sum())
            n_nonsig = len(self.ranked_genes) - n_sig
            weight_nonsig = float(
                self.ranked_genes.loc[
                    self.ranked_genes['p_value'] >= alpha, 'weight'
                ].iloc[0]
            ) if n_nonsig > 0 else 0.5
            # Observed-vs-expected comparator: under the uniform-p null,
            # we expect round(alpha * n_genes) genes to pass by chance.
            # Showing this makes the ranked-genes table self-diagnostic:
            # "20 vs 10" reads as 2x excess = signal; "7 vs 10" reads as
            # at-or-below chance = likely noise. One integer, no
            # math to do mentally.
            n_expected_chance = int(round(alpha * len(self.ranked_genes)))
            lines.append(
                f"  Significance calibration: {n_sig} genes pass alpha={alpha} "
                f"[chance expectation: ~{n_expected_chance}]; "
                f"{n_nonsig} genes down-weighted by {weight_nonsig}x."
            )

        lines.append("\nCLASSIFICATION PERFORMANCE:")
        for name, res in self.classification_results.items():
            marker = " *" if name == self.best_classifier else ""
            lines.append(f"  {name}: AUC={res.auc:.3f}, F1={res.f1:.3f}{marker}")

        lines.append(f"\nRECOMMENDED PANEL: {self.panel_size} genes")
        if self.panel_optimization:
            lines.append(
                f"  Panel AUC: {self.panel_optimization.optimal_auc:.3f}"
            )

        if self.panel_stability is not None:
            ps = self.panel_stability
            lo, hi = ps.nogueira_stability_ci
            lines.append(f"\nPANEL STABILITY (Nogueira 2018, cross-fold):")
            lines.append(
                f"  Phi = {ps.nogueira_stability:.3f} "
                f"[95% CI: {lo:.3f}, {hi:.3f}] ({ps.benchmark_label})"
            )
            lines.append(
                f"  Folds observed: {len(ps.per_fold_panels)} "
                f"({ps.n_folds} folds x {ps.n_repeats} repeats)"
            )

        if self.survival_result and self.survival_result.c_index > 0:
            lines.append(f"\nSURVIVAL ANALYSIS:")
            lines.append(f"  C-index: {self.survival_result.c_index:.3f}")
            lines.append(
                f"  Prognostic genes: {len(self.survival_result.significant_genes)}"
            )

        if self.annotation_result is not None:
            lines.append(f"\nBIOLOGICAL ANNOTATION:")
            lines.append(f"  Annotated genes: {self.annotation_result.n_annotated}")
            lines.append(
                f"  Enriched pathways (FDR<0.05): "
                f"{self.annotation_result.n_enriched_pathways}"
            )
            if not self.annotation_result.literature_hits.empty:
                lines.append(
                    f"  Genes with literature hits: "
                    f"{len(self.annotation_result.literature_hits)}"
                )
            if self.annotation_result.ppi_network is not None:
                n_edges = len(
                    self.annotation_result.ppi_network.get('edges', [])
                )
                lines.append(f"  PPI network edges: {n_edges}")

        lines.append("")
        lines.append("=" * 70)
        return "\n".join(lines)

    def save(self, output_dir: Union[str, Path]):
        """Save all results to output directory."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Ranked genes
        self.ranked_genes.to_csv(output_dir / "ranked_genes.csv")

        # 2. Panel
        panel_df = pd.DataFrame({'gene_id': self.panel, 'rank': range(1, len(self.panel) + 1)})
        panel_df.to_csv(output_dir / "biomarker_panel.csv", index=False)

        # 3. Classification performance
        clf_rows = []
        for name, res in self.classification_results.items():
            clf_rows.append({
                'classifier': name,
                'accuracy': res.accuracy,
                'auc': res.auc,
                'sensitivity': res.sensitivity,
                'specificity': res.specificity,
                'f1': res.f1,
            })
        if clf_rows:
            pd.DataFrame(clf_rows).to_csv(
                output_dir / "classification_performance.csv", index=False
            )

        # 4. Panel optimization curve
        if self.panel_optimization:
            self.panel_optimization.panel_curve.to_csv(
                output_dir / "panel_curve.csv", index=False
            )

        # 5. Survival results
        if self.survival_result and self.survival_result.cox_results is not None:
            self.survival_result.cox_results.to_csv(
                output_dir / "cox_regression.csv"
            )

        # 6. Annotation results
        if self.annotation_result is not None:
            ann_dir = output_dir / "annotations"
            self.annotation_result.save(ann_dir)

        # 7. M6: Panel stability diagnostics
        if self.panel_stability is not None:
            ps = self.panel_stability
            # Per-fold panels: long-form CSV with fold_idx, gene_id
            rows = []
            for fold_i, panel in enumerate(ps.per_fold_panels):
                for g in panel:
                    rows.append({'fold_index': fold_i, 'gene_id': g})
            if rows:
                pd.DataFrame(rows).to_csv(
                    output_dir / "panel_stability_per_fold.csv", index=False
                )
            # Gene selection frequency across folds
            if not ps.gene_selection_frequency.empty:
                freq_df = ps.gene_selection_frequency.to_frame()
                freq_df.index.name = 'gene_id'
                freq_df = freq_df.sort_values(
                    'selection_frequency', ascending=False
                )
                freq_df.to_csv(output_dir / "panel_stability_frequency.csv")
            # Scalar diagnostics
            ci_lo, ci_hi = ps.nogueira_stability_ci
            stab_summary = {
                'nogueira_stability': ps.nogueira_stability,
                'nogueira_stability_ci_lower': ci_lo,
                'nogueira_stability_ci_upper': ci_hi,
                'benchmark_label': ps.benchmark_label,
                'n_folds': ps.n_folds,
                'n_repeats': ps.n_repeats,
                'n_panels_observed': len(ps.per_fold_panels),
                'reference': 'Nogueira, Sechidis & Brown (2018) JMLR 18:174',
            }
            with open(
                output_dir / "panel_stability_summary.json", 'w',
                encoding='utf-8',
            ) as f:
                json.dump(stab_summary, f, indent=2, default=str)

        # 8. Summary text
        with open(output_dir / "summary.txt", 'w', encoding='utf-8') as f:
            f.write(self.summary())

        # 9. Parameters JSON
        params_out = {
            'study_design': self.study_design,
            'validation_strategy': self.validation_strategy,
            'n_samples': self.n_samples,
            'n_initial_candidates': self.n_initial_candidates,
            'panel_size': self.panel_size,
            'panel_genes': self.panel,
            'best_classifier': self.best_classifier,
            'parameters': self.parameters,
            'timestamp': self.timestamp,
        }
        with open(output_dir / "biomarker_params.json", 'w', encoding='utf-8') as f:
            json.dump(params_out, f, indent=2, default=str)

        # 10. Pickle for downstream use
        with open(output_dir / "biomarker_result.pkl", 'wb') as f:
            pickle.dump(self, f)

        logger.info(f"Results saved to: {output_dir}")

    @classmethod
    def load(cls, output_dir: Union[str, Path]) -> 'BiomarkerResult':
        """Load saved BiomarkerResult from pickle."""
        pkl_path = Path(output_dir) / "biomarker_result.pkl"
        with open(pkl_path, 'rb') as f:
            return pickle.load(f)


# =============================================================================
# 10A: FEATURE SELECTION ENGINE
# =============================================================================

class FeatureSelector:
    """
    Multi-method feature selection for biomarker discovery.

    Applies multiple feature selection approaches and aggregates
    results into a consensus ranking using rank aggregation.

    Parameters
    ----------
    random_state : int
        Random seed for reproducibility.
    n_jobs : int
        Number of parallel jobs (-1 for all CPUs).
    verbose : bool
        Print progress messages.
    """

    def __init__(self, random_state: int = DEFAULT_RANDOM_STATE,
                 n_jobs: int = -1, verbose: bool = True):
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.results: Dict[str, FeatureSelectionResult] = {}

    def _log(self, msg):
        if self.verbose:
            logger.info(msg)

    # ---- DE-based filter (from M7/M8/M9) ----

    def select_de_filter(
        self,
        de_genes: List[str],
        all_genes: List[str],
        label: str = 'de_filter'
    ) -> FeatureSelectionResult:
        """
        Use DE-significant genes as initial filter.

        Parameters
        ----------
        de_genes : List[str]
            Gene IDs that are significant from M7/M8/M9.
        all_genes : List[str]
            All gene IDs in the expression matrix.

        Returns
        -------
        FeatureSelectionResult

        Leakage caveat
        --------------
        When ``de_genes`` was computed via DE analysis on the same samples
        that will later be cross-validated in this module, the DE step
        has already seen every label and those genes carry upstream
        selection bias that M10's pipeline-CV cannot fix. If possible,
        supply ``de_genes`` from an independent cohort, or prefer
        ``select_univariate_filter`` which recomputes per-fold on
        training data only (Lewis et al. 2023 nestedcv; Haury et al. 2011
        PLOS ONE).
        """
        self._log(f"   DE filter: {len(de_genes)} significant genes as candidates")
        scores = pd.DataFrame(
            {'score': [1.0 if g in de_genes else 0.0 for g in all_genes]},
            index=all_genes,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        result = FeatureSelectionResult(
            method=label,
            selected_genes=list(de_genes),
            gene_scores=scores,
            n_selected=len(de_genes),
            parameters={'source': 'de_results'},
        )
        self.results[label] = result
        return result

    # ---- Univariate filter (fold-safe replacement for de_filter) ----

    def select_univariate_filter(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_features: int = 50,
        test: str = 'welch',
        label: str = 'univariate_filter',
    ) -> FeatureSelectionResult:
        """
        Per-fold-safe univariate filter (Welch's t-test or Mann-Whitney).

        This method is the "fold-safe" replacement for ``select_de_filter``
        inside the pipeline-CV loop (M6). It runs a gene-by-gene two-sample
        test on the supplied expression matrix and labels, ranks genes by
        -log10(p-value), and selects the top ``n_features``. Because the
        computation uses only the (X, y) it was passed, it is safe to
        invoke inside a training fold with no upstream leakage.

        In cross-genome benchmarks Haury, Gestraud & Vert (2011, PLOS
        ONE) found that simple univariate filters performed as well as or
        better than complex wrapper / embedded methods on both predictive
        accuracy and stability, which is why this is the default filter
        in the ``nestedcv`` R package (Lewis et al. 2023, Bioinformatics
        Advances).

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes). Should already be on a
            sensible scale (raw counts are acceptable; the tests use
            rank-robust statistics so normalization is not strictly
            required, but log-transformed or variance-stabilized data is
            preferred).
        y : np.ndarray
            Binary labels (0/1).
        n_features : int, default 50
            Number of top genes to select by p-value.
        test : {'welch', 'mann_whitney'}, default 'welch'
            'welch': Welch's t-test (two-sided), robust to unequal variances.
            'mann_whitney': Non-parametric Mann-Whitney U test. Slower but
                    robust to outliers and non-normal distributions.
        label : str, default 'univariate_filter'
            Name to store this result under in ``self.results``.

        Returns
        -------
        FeatureSelectionResult
            Scores are -log10(p-value) so higher score == more discriminating.
            NaN p-values (from genes with zero variance in one or both
            classes) are converted to score 0.
        """
        from scipy import stats as sp_stats

        X_arr = np.asarray(X.values, dtype=float)
        y_arr = np.asarray(y).ravel()

        unique_y = np.unique(y_arr)
        if len(unique_y) != 2:
            raise ValueError(
                f"univariate_filter requires exactly 2 classes, got "
                f"{len(unique_y)} ({list(unique_y)})."
            )

        mask0 = y_arr == unique_y[0]
        mask1 = y_arr == unique_y[1]
        X0 = X_arr[mask0, :]
        X1 = X_arr[mask1, :]

        n_genes = X_arr.shape[1]
        p_values = np.full(n_genes, 1.0, dtype=float)

        if test == 'welch':
            # Vectorized Welch's t-test across all genes at once.
            # scipy.stats.ttest_ind handles the whole matrix with axis=0.
            with np.errstate(all='ignore'):
                t_stat, p_vals = sp_stats.ttest_ind(
                    X1, X0, axis=0, equal_var=False, nan_policy='omit'
                )
            p_values = np.asarray(p_vals, dtype=float)
        elif test == 'mann_whitney':
            # Mann-Whitney isn't vectorized over columns in scipy; loop.
            for j in range(n_genes):
                try:
                    with np.errstate(all='ignore'):
                        _, pv = sp_stats.mannwhitneyu(
                            X1[:, j], X0[:, j], alternative='two-sided'
                        )
                    p_values[j] = float(pv)
                except ValueError:
                    p_values[j] = 1.0
        else:
            raise ValueError(
                f"test must be 'welch' or 'mann_whitney', got {test!r}."
            )

        # Replace NaN p-values (zero-variance genes) with 1.0 (no signal)
        p_values = np.where(np.isnan(p_values), 1.0, p_values)
        # Floor to avoid log(0)
        p_floor = np.maximum(p_values, 1e-300)
        neg_log_p = -np.log10(p_floor)

        scores = pd.DataFrame(
            {'score': neg_log_p, 'p_value': p_values},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        n_select = min(int(n_features), n_genes)
        selected = scores.head(n_select).index.tolist()

        self._log(
            f"   Univariate filter ({test}): top {n_select} of {n_genes} genes "
            f"(min p = {p_values.min():.2e})"
        )

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores[['score']],
            n_selected=len(selected),
            parameters={'test': test, 'n_features': n_select},
        )
        self.results[label] = result
        return result

    # ---- LASSO / Elastic Net ----

    def select_lasso(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        alpha: float = 0.5,
        l1_ratio: float = 0.5,
        max_iter: int = 5000,
        label: str = 'elastic_net'
    ) -> FeatureSelectionResult:
        """
        Elastic Net feature selection via penalized logistic regression.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes).
        y : np.ndarray
            Binary labels (0/1).
        alpha : float
            Regularization strength (C = 1/alpha in sklearn).
        l1_ratio : float
            L1 vs L2 mixing: 1.0 = pure LASSO, 0.0 = pure Ridge.
        """
        if not _SKLEARN_AVAILABLE:
            raise DependencyError("scikit-learn is required for LASSO/Elastic Net")

        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        self._log(f"   Elastic Net (l1_ratio={l1_ratio}): fitting...")

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # sklearn uses C = 1/alpha; penalty='elasticnet' needs saga solver
        model = LogisticRegression(
            C=1.0 / max(alpha, 1e-6),
            l1_ratio=l1_ratio,
            solver='saga',
            max_iter=max_iter,
            random_state=self.random_state,
        )
        model.fit(X_scaled, y)

        # Gene importance = absolute coefficient
        coefs = np.abs(model.coef_).ravel()
        scores = pd.DataFrame(
            {'score': coefs},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        selected = scores[scores['score'] > 0].index.tolist()
        self._log(f"   Elastic Net: {len(selected)} non-zero coefficients")

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores,
            n_selected=len(selected),
            parameters={'alpha': alpha, 'l1_ratio': l1_ratio, 'max_iter': max_iter},
        )
        self.results[label] = result
        return result

    # ---- Boruta ----

    def select_boruta(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        max_iter: int = 100,
        label: str = 'boruta'
    ) -> FeatureSelectionResult:
        """
        Boruta feature selection using Random Forest shadow features.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes).
        y : np.ndarray
            Binary labels.
        max_iter : int
            Maximum Boruta iterations.
        """
        if not _BORUTA_AVAILABLE:
            raise DependencyError(
                "boruta_py is required: pip install boruta"
            )
        if not _SKLEARN_AVAILABLE:
            raise DependencyError("scikit-learn is required for Boruta")

        from boruta import BorutaPy
        from sklearn.ensemble import RandomForestClassifier

        self._log(f"   Boruta (max_iter={max_iter}): fitting...")

        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=5,
            random_state=self.random_state,
            n_jobs=self.n_jobs,
        )

        boruta = BorutaPy(
            rf,
            n_estimators='auto',
            max_iter=max_iter,
            random_state=self.random_state,
            verbose=0,
        )
        boruta.fit(X.values, y)

        # Ranking: 1 = confirmed, 2 = tentative, 3+ = rejected
        # Convert to scores: lower rank = higher score
        max_rank = boruta.ranking_.max()
        importance = (max_rank - boruta.ranking_ + 1) / max_rank

        scores = pd.DataFrame(
            {'score': importance},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        # Selected = confirmed + tentative
        selected_mask = boruta.support_ | boruta.support_weak_
        selected = X.columns[selected_mask].tolist()
        self._log(f"   Boruta: {len(selected)} features selected "
                  f"({int(boruta.support_.sum())} confirmed, "
                  f"{int(boruta.support_weak_.sum())} tentative)")

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores,
            n_selected=len(selected),
            parameters={'max_iter': max_iter},
        )
        self.results[label] = result
        return result

    # ---- mRMR ----

    def select_mrmr(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_features: int = 50,
        label: str = 'mrmr'
    ) -> FeatureSelectionResult:
        """
        Minimum Redundancy Maximum Relevance feature selection.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes).
        y : np.ndarray
            Binary labels.
        n_features : int
            Number of features to select.
        """
        if not _MRMR_AVAILABLE:
            raise DependencyError(
                "mrmr-selection is required: pip install mrmr-selection"
            )

        import mrmr

        self._log(f"   mRMR (K={n_features}): selecting...")

        # mrmr expects y as a Series
        y_series = pd.Series(y, index=X.index, name='target')
        selected = mrmr.mrmr_classif(X, y_series, K=n_features)

        # Build scores from selection order (first selected = highest score)
        n_total = len(X.columns)
        score_map = {gene: (n_features - i) / n_features
                     for i, gene in enumerate(selected)}

        scores = pd.DataFrame(
            {'score': [score_map.get(g, 0.0) for g in X.columns]},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        self._log(f"   mRMR: {len(selected)} features selected")

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores,
            n_selected=len(selected),
            parameters={'n_features': n_features},
        )
        self.results[label] = result
        return result

    # ---- Recursive Feature Elimination ----

    def select_rfe(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_features: int = 50,
        step: float = 0.1,
        label: str = 'rfe'
    ) -> FeatureSelectionResult:
        """
        Recursive Feature Elimination with Random Forest.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes).
        y : np.ndarray
            Binary labels.
        n_features : int
            Number of features to select.
        step : float
            Fraction of features to remove at each step.
        """
        if not _SKLEARN_AVAILABLE:
            raise DependencyError("scikit-learn is required for RFE")

        from sklearn.feature_selection import RFE
        from sklearn.ensemble import RandomForestClassifier

        self._log(f"   RFE (target={n_features}): running...")

        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=7,
            random_state=self.random_state,
            n_jobs=self.n_jobs,
        )

        rfe = RFE(
            estimator=rf,
            n_features_to_select=n_features,
            step=step,
        )
        rfe.fit(X.values, y)

        # ------------------------------------------------------------------
        # M3 fix: produce unique per-feature scores instead of mass ties.
        # ------------------------------------------------------------------
        # sklearn's ``RFE.ranking_`` assigns 1 to every one of the
        # ``n_features_to_select`` chosen features. Using ``1 / ranking_``
        # (or any monotonic transform of ranking_) therefore produces
        # ``n_features_to_select`` tied scores at the top. Downstream
        # consensus ranking then collapses all of those ties to a single
        # mid-rank (e.g. rank 25.5 when 50 features are selected), which
        # means RFE contributes effectively zero discrimination to the
        # consensus.
        #
        # We replace that with a two-tier score:
        #
        #   * Selected features (``support_ == True``): use the final
        #     estimator's per-feature importance signal. For tree-based
        #     estimators we read ``feature_importances_``; for linear
        #     estimators we read ``|coef_|``. Those values are rescaled
        #     into ``[1.0, 2.0]`` so selected features always outrank
        #     non-selected ones.
        #   * Non-selected features: use ``1 / ranking_``. sklearn's
        #     ``ranking_`` is unique among eliminated features (they get
        #     values 2, 3, 4, ... reflecting elimination order), so the
        #     reciprocal produces unique scores in ``(0, 0.5]`` that are
        #     strictly below any selected feature's score.
        #
        # Invariants:
        #   - Every selected feature's score > every non-selected score,
        #     so consensus ranks always place selected genes first.
        #   - Selected scores are unique (except in the degenerate case
        #     where the estimator assigns identical importance to every
        #     selected feature, in which case we fall back to a flat 1.5
        #     that still preserves the selected-vs-non-selected ordering).
        #   - Non-selected scores are always unique by construction.
        #
        # Note: a similar ties pattern exists in ``select_boruta`` (three
        # buckets: confirmed / tentative / rejected). That's a separate
        # issue and is currently less harmful because Boruta isn't in the
        # default method list; leaving it for a future pass.
        n_total = X.shape[1]
        selected_mask = rfe.support_

        estimator_importance = None
        if hasattr(rfe.estimator_, 'coef_'):
            # Linear estimator path: absolute coefficient magnitude.
            estimator_importance = np.abs(rfe.estimator_.coef_).ravel()
        elif hasattr(rfe.estimator_, 'feature_importances_'):
            # Tree estimator path (the default RFE estimator here).
            estimator_importance = rfe.estimator_.feature_importances_

        importance = np.empty(n_total, dtype=float)

        n_selected = int(selected_mask.sum())
        if (
            estimator_importance is not None
            and len(estimator_importance) == n_selected
        ):
            imp_min = float(estimator_importance.min())
            imp_max = float(estimator_importance.max())
            if imp_max > imp_min:
                # Linearly rescale into [1.0, 2.0].
                selected_scores = (
                    1.0
                    + (estimator_importance - imp_min) / (imp_max - imp_min)
                )
            else:
                # All selected features have identical importance — rare
                # but can happen on very small data. Flat fallback keeps
                # the selected-vs-non-selected ordering intact.
                selected_scores = np.full(n_selected, 1.5, dtype=float)
            importance[selected_mask] = selected_scores
        else:
            # Estimator exposes no importance signal (or shape mismatch
            # e.g. multiclass coef_). Flat fallback: all selected score
            # 1.5, all non-selected score < 1.0. Preserves the primary
            # invariant (selected outrank non-selected) and no worse than
            # the pre-M3 behavior within the selected block.
            importance[selected_mask] = 1.5

        # Non-selected features get 1 / elimination_rank, giving unique
        # scores in (0, 0.5] ordered by elimination order.
        importance[~selected_mask] = 1.0 / rfe.ranking_[~selected_mask]

        scores = pd.DataFrame(
            {'score': importance},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        selected = X.columns[rfe.support_].tolist()
        self._log(f"   RFE: {len(selected)} features selected")

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores,
            n_selected=len(selected),
            parameters={'n_features': n_features, 'step': step},
        )
        self.results[label] = result
        return result

    # ---- SHAP-based ranking ----

    def select_shap(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_features: int = 50,
        label: str = 'shap'
    ) -> FeatureSelectionResult:
        """
        SHAP-based feature selection using XGBoost or Random Forest.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes).
        y : np.ndarray
            Binary labels.
        n_features : int
            Number of top features to select.
        """
        if not _SHAP_AVAILABLE:
            raise DependencyError("shap is required: pip install shap")

        import shap

        self._log(f"   SHAP (top={n_features}): computing importance...")

        # Use XGBoost if available, else Random Forest
        if _XGBOOST_AVAILABLE:
            from xgboost import XGBClassifier
            model = XGBClassifier(
                n_estimators=100,
                max_depth=5,
                random_state=self.random_state,
                n_jobs=self.n_jobs,
                use_label_encoder=False,
                eval_metric='logloss',
                verbosity=0,
            )
        else:
            from sklearn.ensemble import RandomForestClassifier
            model = RandomForestClassifier(
                n_estimators=100,
                max_depth=5,
                random_state=self.random_state,
                n_jobs=self.n_jobs,
            )

        model.fit(X.values, y)

        # SHAP values
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X.values)

        # For binary: shap_values may be a list [class0, class1] or 2D array
        if isinstance(shap_values, list):
            shap_abs = np.abs(shap_values[1])  # Positive class
        elif shap_values.ndim == 3:
            shap_abs = np.abs(shap_values[:, :, 1])
        else:
            shap_abs = np.abs(shap_values)

        mean_shap = shap_abs.mean(axis=0)

        scores = pd.DataFrame(
            {'score': mean_shap},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        selected = scores.head(n_features).index.tolist()
        self._log(f"   SHAP: top {len(selected)} features selected")

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores,
            n_selected=len(selected),
            parameters={'n_features': n_features, 'base_model': type(model).__name__},
        )
        self.results[label] = result
        return result

    # ---- WGCNA hub genes ----

    def select_wgcna(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_hub_genes: int = 10,
        trait_corr_threshold: float = 0.3,
        min_module_size: int = 30,
        species: str = 'mus musculus',
        label: str = 'wgcna'
    ) -> FeatureSelectionResult:
        """
        WGCNA co-expression network analysis for hub gene selection.

        Builds a weighted co-expression network, identifies gene modules,
        correlates modules with the clinical trait, and extracts hub genes
        from significantly correlated modules.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes), log2-CPM normalized.
        y : np.ndarray
            Binary trait labels (0/1).
        n_hub_genes : int
            Number of top hub genes to extract per significant module.
        trait_corr_threshold : float
            Minimum absolute Pearson correlation between module eigengene
            and trait for a module to be considered significant.
        min_module_size : int
            Minimum number of genes per module.
        species : str
            Species name for PyWGCNA ('mus musculus', 'homo sapiens', etc.).

        Returns
        -------
        FeatureSelectionResult

        Scientific Basis
        ----------------
        Langfelder & Horvath (2008). WGCNA: an R package for weighted
        correlation network analysis. BMC Bioinformatics, 9, 559.

        Rezaie, Reese & Mortazavi (2023). PyWGCNA: a Python package for
        weighted gene co-expression network analysis. Bioinformatics, 39(7).
        """
        if not _PYWGCNA_AVAILABLE:
            raise DependencyError(
                "PyWGCNA is required: pip install PyWGCNA"
            )

        import PyWGCNA
        import anndata as ad

        self._log(f"   WGCNA: building co-expression network...")

        # PyWGCNA requires ≥15 samples
        n_samples = X.shape[0]
        if n_samples < 15:
            raise ValueError(
                f"WGCNA requires at least 15 samples, got {n_samples}. "
                f"Consider using other feature selection methods."
            )

        # Prepare AnnData input: PyWGCNA expects samples as obs, genes as var
        # X is already (samples x genes) from _prepare_expression_data
        adata = ad.AnnData(X)
        adata.obs['trait'] = y.astype(float)

        # Create PyWGCNA object
        # Note: PyWGCNA uses the WGCNA class, not pyWGCNA function
        wgcna_obj = PyWGCNA.WGCNA(
            name='raptor_m10',
            species=species,
            geneExp=adata,
            outputPath='',
            save=False,
        )

        # Preprocess: filter low-variance genes
        # (our data is already filtered, but PyWGCNA expects this step)
        try:
            wgcna_obj.preprocess(
                TPMcutoff=0,  # Already filtered upstream
                cut=0.0,      # No additional sample filtering
            )
        except Exception:
            # Some versions have different preprocessing API
            pass

        # Find modules (network construction + module detection)
        self._log(f"   WGCNA: detecting modules (min_size={min_module_size})...")
        try:
            wgcna_obj.findModules(
                minModuleSize=min_module_size,
                networkType='signed',
            )
        except Exception as e:
            self._log(f"   WGCNA findModules failed: {e}")
            self._log(f"   Trying with default parameters...")
            wgcna_obj.findModules()

        # Get module assignments from the processed AnnData
        datExpr = wgcna_obj.datExpr
        gene_info = datExpr.var

        if 'moduleColors' not in gene_info.columns and 'moduleLabels' not in gene_info.columns:
            self._log("   WGCNA: no modules detected, skipping")
            empty_scores = pd.DataFrame(
                {'score': 0.0}, index=X.columns
            )
            empty_scores.index.name = 'gene_id'
            return FeatureSelectionResult(
                method=label, selected_genes=[], gene_scores=empty_scores,
                n_selected=0, parameters={'error': 'no modules detected'},
            )

        # Module color assignments
        module_col = 'moduleColors' if 'moduleColors' in gene_info.columns else 'moduleLabels'
        module_assignments = gene_info[module_col]
        unique_modules = [m for m in module_assignments.unique() if m != 'grey']

        self._log(f"   WGCNA: found {len(unique_modules)} modules (excl. grey)")

        # Calculate module-trait correlation
        # Compute module eigengenes and correlate with trait
        trait = datExpr.obs['trait'].values if 'trait' in datExpr.obs.columns else y

        significant_modules = []
        module_correlations = {}

        for module in unique_modules:
            module_genes = module_assignments[module_assignments == module].index
            if len(module_genes) < 3:
                continue

            # Module eigengene = first PC of module gene expression
            module_expr = datExpr[:, module_genes].X
            if hasattr(module_expr, 'toarray'):
                module_expr = module_expr.toarray()
            if isinstance(module_expr, pd.DataFrame):
                module_expr = module_expr.values

            try:
                from sklearn.decomposition import PCA
                pca = PCA(n_components=1)
                eigengene = pca.fit_transform(module_expr).ravel()

                # Pearson correlation with trait
                corr, pval = stats.pearsonr(eigengene, trait)
                module_correlations[module] = {
                    'correlation': corr,
                    'p_value': pval,
                    'n_genes': len(module_genes),
                }

                if abs(corr) >= trait_corr_threshold and pval < 0.05:
                    significant_modules.append(module)

            except Exception:
                continue

        self._log(
            f"   WGCNA: {len(significant_modules)} modules significantly "
            f"correlated with trait (|r| >= {trait_corr_threshold}, p < 0.05)"
        )

        # Extract hub genes from significant modules
        # Hub gene = highest module membership (kME = correlation with eigengene)
        all_hub_genes = []
        gene_scores_dict = {}

        for module in significant_modules:
            module_genes = module_assignments[module_assignments == module].index
            module_expr = datExpr[:, module_genes].X
            if hasattr(module_expr, 'toarray'):
                module_expr = module_expr.toarray()

            try:
                # Compute eigengene
                pca = PCA(n_components=1)
                eigengene = pca.fit_transform(module_expr).ravel()

                # kME = correlation of each gene with module eigengene
                kme_scores = {}
                for i, gene in enumerate(module_genes):
                    gene_expr = module_expr[:, i]
                    corr, _ = stats.pearsonr(gene_expr, eigengene)
                    kme_scores[gene] = abs(corr)

                # Sort by kME, take top hub genes
                sorted_genes = sorted(
                    kme_scores.items(), key=lambda x: x[1], reverse=True
                )

                module_corr = abs(module_correlations[module]['correlation'])
                for gene, kme in sorted_genes[:n_hub_genes]:
                    # Score combines kME (intra-module connectivity)
                    # with module-trait correlation
                    gene_scores_dict[gene] = kme * module_corr
                    all_hub_genes.append(gene)

            except Exception:
                continue

        # Build scores DataFrame for all genes
        scores = pd.DataFrame(
            {'score': [gene_scores_dict.get(g, 0.0) for g in X.columns]},
            index=X.columns,
        )
        scores.index.name = 'gene_id'
        scores = scores.sort_values('score', ascending=False)

        # Deduplicate hub genes (gene may appear in multiple modules)
        selected = list(dict.fromkeys(all_hub_genes))

        self._log(
            f"   WGCNA: {len(selected)} hub genes from "
            f"{len(significant_modules)} significant modules"
        )

        result = FeatureSelectionResult(
            method=label,
            selected_genes=selected,
            gene_scores=scores,
            n_selected=len(selected),
            parameters={
                'n_hub_genes_per_module': n_hub_genes,
                'trait_corr_threshold': trait_corr_threshold,
                'min_module_size': min_module_size,
                'n_modules_total': len(unique_modules),
                'n_modules_significant': len(significant_modules),
                'module_correlations': module_correlations,
            },
        )
        self.results[label] = result
        return result

    # ---- Consensus ranking ----

    def consensus_ranking(
        self,
        all_genes: List[str],
    ) -> pd.DataFrame:
        """
        Aggregate feature selection results into consensus ranking.

        Uses rank aggregation (average rank across methods) to produce
        a unified gene ranking. Genes selected by more methods rank higher.

        Parameters
        ----------
        all_genes : List[str]
            Full list of genes to rank.

        Returns
        -------
        pd.DataFrame
            Columns: consensus_rank, consensus_score, n_methods_selected,
            plus per-method rank columns.
        """
        if not self.results:
            raise ValueError("No feature selection results available. Run methods first.")

        n_genes = len(all_genes)
        method_names = list(self.results.keys())
        rank_matrix = pd.DataFrame(index=all_genes)

        for method, res in self.results.items():
            # Convert scores to ranks (1 = best)
            method_scores = res.gene_scores.reindex(all_genes)
            method_scores = method_scores.fillna(0.0)
            # Rank: higher score = lower rank number (i.e. better)
            rank_matrix[f'rank_{method}'] = method_scores['score'].rank(
                ascending=False, method='average'
            )

        # Average rank across methods
        rank_cols = [c for c in rank_matrix.columns if c.startswith('rank_')]
        rank_matrix['mean_rank'] = rank_matrix[rank_cols].mean(axis=1)

        # Number of methods that selected each gene
        selection_count = pd.Series(0, index=all_genes, dtype=int)
        for method, res in self.results.items():
            for gene in res.selected_genes:
                if gene in selection_count.index:
                    selection_count[gene] += 1

        rank_matrix['n_methods_selected'] = selection_count

        # Consensus score: combine normalized rank with selection count
        # Higher score = better candidate
        max_rank = rank_matrix['mean_rank'].max()
        rank_matrix['consensus_score'] = (
            (max_rank - rank_matrix['mean_rank']) / max_rank * 0.7 +
            rank_matrix['n_methods_selected'] / len(method_names) * 0.3
        )

        # Final consensus rank
        rank_matrix['consensus_rank'] = rank_matrix['consensus_score'].rank(
            ascending=False, method='first'
        ).astype(int)

        rank_matrix = rank_matrix.sort_values('consensus_rank')
        rank_matrix.index.name = 'gene_id'

        return rank_matrix


# =============================================================================
# 10B: CLASSIFICATION & VALIDATION ENGINE
# =============================================================================

class ClassifierEvaluator:
    """
    Evaluate biomarker panels using classification models with
    rigorous cross-validation.

    Parameters
    ----------
    random_state : int
        Random seed.
    n_jobs : int
        Number of parallel jobs.
    verbose : bool
        Print progress.
    """

    def __init__(self, random_state: int = DEFAULT_RANDOM_STATE,
                 n_jobs: int = -1, verbose: bool = True):
        if not _SKLEARN_AVAILABLE:
            raise DependencyError("scikit-learn is required for classification")

        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

    def _log(self, msg):
        if self.verbose:
            logger.info(msg)

    def _get_classifiers(self) -> Dict[str, Any]:
        """Build dictionary of classifiers."""
        from sklearn.linear_model import LogisticRegression
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.svm import SVC

        clfs = {
            'logistic_regression': LogisticRegression(
                C=1.0, l1_ratio=0, max_iter=5000,
                random_state=self.random_state,
            ),
            'random_forest': RandomForestClassifier(
                n_estimators=200, max_depth=None,
                random_state=self.random_state, n_jobs=self.n_jobs,
            ),
            'svm': SVC(
                kernel='linear', C=1.0, probability=True,
                random_state=self.random_state,
            ),
        }

        if _XGBOOST_AVAILABLE:
            from xgboost import XGBClassifier
            clfs['xgboost'] = XGBClassifier(
                n_estimators=100, max_depth=5,
                random_state=self.random_state, n_jobs=self.n_jobs,
                use_label_encoder=False, eval_metric='logloss',
                verbosity=0,
            )

        return clfs

    def evaluate_nested_cv(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_outer: int = DEFAULT_N_FOLDS_OUTER,
        n_inner: int = DEFAULT_N_FOLDS_INNER,
        classifiers: Optional[List[str]] = None,
    ) -> Dict[str, ClassificationResult]:
        """
        Evaluate classifiers using nested cross-validation.

        Outer loop: performance estimation.
        Inner loop: hyperparameter tuning (not implemented in v1 — uses defaults).

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes), already restricted to panel.
        y : np.ndarray
            Binary labels.
        n_outer : int
            Number of outer CV folds.
        classifiers : List[str], optional
            Which classifiers to run. Default: all available.

        Returns
        -------
        Dict[str, ClassificationResult]
        """
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import (
            accuracy_score, roc_auc_score, f1_score,
            confusion_matrix, roc_curve,
        )
        from sklearn.preprocessing import StandardScaler

        all_clfs = self._get_classifiers()
        if classifiers:
            all_clfs = {k: v for k, v in all_clfs.items() if k in classifiers}

        results = {}
        outer_cv = StratifiedKFold(
            n_splits=n_outer, shuffle=True, random_state=self.random_state
        )

        for clf_name, clf_template in all_clfs.items():
            self._log(f"   Evaluating {clf_name}...")

            fold_metrics = []
            all_y_true = []
            all_y_prob = []

            for fold_i, (train_idx, test_idx) in enumerate(
                outer_cv.split(X, y)
            ):
                X_train = X.iloc[train_idx]
                X_test = X.iloc[test_idx]
                y_train = y[train_idx]
                y_test = y[test_idx]

                # Scale
                scaler = StandardScaler()
                X_train_s = scaler.fit_transform(X_train)
                X_test_s = scaler.transform(X_test)

                # Clone and fit
                from sklearn.base import clone
                clf = clone(clf_template)
                clf.fit(X_train_s, y_train)

                # Predict
                y_pred = clf.predict(X_test_s)
                y_prob = clf.predict_proba(X_test_s)[:, 1]

                # Metrics
                acc = accuracy_score(y_test, y_pred)
                try:
                    auc = roc_auc_score(y_test, y_prob)
                except ValueError:
                    auc = 0.5

                f1 = f1_score(y_test, y_pred, zero_division=0)

                cm = confusion_matrix(y_test, y_pred)
                tn, fp, fn, tp = cm.ravel() if cm.size == 4 else (0, 0, 0, 0)
                sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0

                fold_metrics.append({
                    'fold': fold_i,
                    'accuracy': acc,
                    'auc': auc,
                    'f1': f1,
                    'sensitivity': sens,
                    'specificity': spec,
                })

                all_y_true.extend(y_test)
                all_y_prob.extend(y_prob)

            # Aggregate
            mean_acc = np.mean([m['accuracy'] for m in fold_metrics])
            mean_auc = np.mean([m['auc'] for m in fold_metrics])
            mean_f1 = np.mean([m['f1'] for m in fold_metrics])
            mean_sens = np.mean([m['sensitivity'] for m in fold_metrics])
            mean_spec = np.mean([m['specificity'] for m in fold_metrics])

            # ROC on aggregated predictions
            roc_data = None
            try:
                fpr, tpr, thresholds = roc_curve(all_y_true, all_y_prob)
                roc_data = {
                    'fpr': fpr.tolist(),
                    'tpr': tpr.tolist(),
                }
            except Exception:
                pass

            # Feature importance from full-data fit. Store a Pipeline
            # (scaler + classifier) as trained_model so downstream calls
            # to predict_proba / decision_function automatically apply
            # the same StandardScaler the model was fit with. Storing a
            # bare estimator trained on scaled data causes probability
            # predictions on raw expression values to collapse to the
            # extreme tails of the sigmoid — a silent bug that only
            # surfaces in downstream clinical metrics (DCA, Youden's
            # threshold applied to new data, etc.).
            from sklearn.pipeline import Pipeline
            scaler_full = StandardScaler()
            clf_final = clone(clf_template)
            full_pipeline = Pipeline([
                ("scaler", scaler_full),
                ("clf", clf_final),
            ])
            full_pipeline.fit(X, y)

            feat_imp = None
            if hasattr(clf_final, 'feature_importances_'):
                feat_imp = pd.DataFrame({
                    'importance': clf_final.feature_importances_
                }, index=X.columns)
            elif hasattr(clf_final, 'coef_'):
                feat_imp = pd.DataFrame({
                    'importance': np.abs(clf_final.coef_).ravel()
                }, index=X.columns)
            if feat_imp is not None:
                feat_imp.index.name = 'gene_id'
                feat_imp = feat_imp.sort_values('importance', ascending=False)

            self._log(
                f"   {clf_name}: AUC={mean_auc:.3f}, "
                f"F1={mean_f1:.3f}, Sens={mean_sens:.3f}, Spec={mean_spec:.3f}"
            )

            results[clf_name] = ClassificationResult(
                model_name=clf_name,
                accuracy=mean_acc,
                auc=mean_auc,
                sensitivity=mean_sens,
                specificity=mean_spec,
                f1=mean_f1,
                metrics_per_fold=fold_metrics,
                roc_data=roc_data,
                feature_importance=feat_imp,
                trained_model=full_pipeline,
                # M1: persist honest out-of-fold predictions. These
                # come from outer-fold held-out evaluations (no sample
                # was used to train the model that predicted it). The
                # clinical metrics layer uses these rather than
                # full-data predictions to avoid reporting
                # training-data performance as generalisation.
                oof_true=np.asarray(all_y_true, dtype=int),
                oof_prob=np.asarray(all_y_prob, dtype=float),
            )

        return results

    def evaluate_loocv(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        classifiers: Optional[List[str]] = None,
    ) -> Dict[str, ClassificationResult]:
        """Leave-one-out cross-validation for small sample sizes.

        Populates accuracy, AUC, F1, sensitivity, and specificity on the
        returned ClassificationResult objects. Sensitivity and specificity
        are computed from the full LOOCV confusion matrix (labels=[0, 1])
        so results are correct even when a classifier collapses to a
        single predicted class.
        """
        from sklearn.model_selection import LeaveOneOut
        from sklearn.metrics import (
            accuracy_score, roc_auc_score, f1_score, confusion_matrix,
        )
        from sklearn.preprocessing import StandardScaler
        from sklearn.base import clone

        all_clfs = self._get_classifiers()
        if classifiers:
            all_clfs = {k: v for k, v in all_clfs.items() if k in classifiers}

        results = {}
        loo = LeaveOneOut()

        for clf_name, clf_template in all_clfs.items():
            self._log(f"   LOOCV: {clf_name}...")

            all_y_true = []
            all_y_pred = []
            all_y_prob = []

            for train_idx, test_idx in loo.split(X):
                X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                scaler = StandardScaler()
                X_train_s = scaler.fit_transform(X_train)
                X_test_s = scaler.transform(X_test)

                clf = clone(clf_template)
                clf.fit(X_train_s, y_train)

                all_y_true.append(y_test[0])
                all_y_pred.append(clf.predict(X_test_s)[0])
                all_y_prob.append(clf.predict_proba(X_test_s)[0, 1])

            all_y_true = np.array(all_y_true)
            all_y_pred = np.array(all_y_pred)
            all_y_prob = np.array(all_y_prob)

            acc = accuracy_score(all_y_true, all_y_pred)
            try:
                auc = roc_auc_score(all_y_true, all_y_prob)
            except ValueError:
                auc = 0.5
            f1 = f1_score(all_y_true, all_y_pred, zero_division=0)

            # Sensitivity and specificity from the full LOOCV confusion
            # matrix. Forcing labels=[0, 1] guarantees a 2x2 matrix even
            # when a classifier collapses to a single predicted class.
            cm = confusion_matrix(all_y_true, all_y_pred, labels=[0, 1])
            tn, fp, fn, tp = cm.ravel()
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

            self._log(
                f"   {clf_name}: AUC={auc:.3f}, Acc={acc:.3f}, "
                f"Sens={sensitivity:.3f}, Spec={specificity:.3f}"
            )

            results[clf_name] = ClassificationResult(
                model_name=clf_name,
                accuracy=acc,
                auc=auc,
                f1=f1,
                sensitivity=sensitivity,
                specificity=specificity,
                # M1: LOOCV produces one out-of-fold prediction per
                # sample (every sample is held out exactly once). These
                # predictions are even more rigorous than k-fold OOF
                # predictions because there is no randomness in fold
                # assignment and no sample overlap across training
                # sets. Clinical metrics downstream use these.
                oof_true=all_y_true.astype(int),
                oof_prob=all_y_prob.astype(float),
            )

        return results


# =============================================================================
# 10C: PANEL OPTIMIZATION
# =============================================================================
#
# Auto-detection of optimal panel size from a panel-size-vs-AUC curve.
#
# Pre-kneedle, RAPTOR used a "first-drop" heuristic that walked the curve
# forward and stopped at the first index where adding a gene gave less
# than 0.5% AUC improvement. This had three known failure modes documented
# in the kneedle scoping doc (April 2026):
#
#   1. Non-monotone curves: a single noisy down-tick early in the curve
#      terminated the loop prematurely, returning a panel size smaller
#      than the true elbow. Common at small n where each panel-size's
#      CV AUC has fold-split noise.
#   2. Saturated curves: when AUC saturates near 1.000 from the very
#      first panel size, the first improvement is < 0.005 and the loop
#      returned min_panel. User asked for up to N genes, got back the
#      floor.
#   3. Late-knee curves: real knees past size 30 were missed because a
#      single sub-threshold step elsewhere terminated the walk.
#
# Kneedle (Satopaa et al. 2011 ICDCSW, "Finding a Kneedle in a Haystack")
# detects the knee by maximum-curvature on a normalized, optionally-
# smoothed curve. Robust to local noise, saturation, and late knees.
# When kneedle returns no knee (e.g., on truly saturated curves with no
# curvature anywhere), the implementation falls back to argmax over the
# curve, with M5-style smallest-size tiebreak.

def _detect_optimal_panel_size(
    curve_df: pd.DataFrame,
    auto_strategy: str = 'kneedle',
    sensitivity: float = 1.0,
) -> Tuple[int, str]:
    """Detect optimal panel size from a panel-size-vs-AUC curve.

    Parameters
    ----------
    curve_df : pd.DataFrame
        Must have columns 'panel_size' and 'auc_mean'.
    auto_strategy : {'kneedle', 'argmax', 'first_drop'}
        Detection method.
            'kneedle'   : kneedle algorithm with polynomial smoothing.
                          Falls back to argmax if no knee is detected
                          or if kneed is unavailable or if there are
                          fewer than 3 points.
            'argmax'    : pick the size with the highest auc_mean. Ties
                          go to the smallest panel size.
            'first_drop': legacy behavior; walk forward, stop at first
                          index where improvement < 0.005. Retained for
                          backward-compat and debugging.
    sensitivity : float
        Kneedle's S parameter (default 1.0).

    Returns
    -------
    (optimal_size, method_used) where method_used is one of:
        'kneedle'         — kneedle found a knee
        'argmax_fallback' — kneedle returned None / unavailable; argmax
        'argmax'          — auto_strategy='argmax' explicitly
        'first_drop'      — auto_strategy='first_drop' explicitly
    """
    sizes = curve_df['panel_size'].values.astype(int)
    aucs = curve_df['auc_mean'].values.astype(float)

    if auto_strategy == 'argmax':
        max_auc = aucs.max()
        winners = sizes[aucs == max_auc]
        return int(winners.min()), 'argmax'

    if auto_strategy == 'first_drop':
        # Legacy behavior, kept verbatim from the pre-kneedle code:
        # optimal_idx tracks the most recent index where improvement was
        # at or above the threshold; stop at the first sub-threshold
        # gap and return optimal_idx (NOT i — that's already the bad
        # one).
        optimal_idx = 0
        for i in range(1, len(aucs)):
            if aucs[i] - aucs[i-1] < 0.005:
                break
            optimal_idx = i
        return int(sizes[optimal_idx]), 'first_drop'

    # auto_strategy == 'kneedle'
    if len(sizes) < 3:
        # Need at least 3 points for curvature analysis
        max_auc = aucs.max()
        winners = sizes[aucs == max_auc]
        return int(winners.min()), 'argmax_fallback'

    # Constant-curve guard: when AUC is essentially flat across all
    # tested sizes, kneedle's polynomial smoothing is poorly conditioned
    # and the answer is degenerate anyway (any size is "optimal"). Skip
    # the kneedle call to avoid a noisy RankWarning, and pick the
    # smallest size at the max AUC.
    if (aucs.max() - aucs.min()) < 1e-6:
        winners = sizes[aucs == aucs.max()]
        return int(winners.min()), 'argmax_fallback'

    try:
        from kneed import KneeLocator
    except ImportError:
        # Defensive: kneed should be in install_requires but if it isn't
        # available we fall back gracefully rather than crash.
        max_auc = aucs.max()
        winners = sizes[aucs == max_auc]
        return int(winners.min()), 'argmax_fallback'

    # Suppress kneedle's internal np.RankWarning. On near-flat or
    # short curves the polynomial fit may be poorly conditioned;
    # the warning is benign because the constant-curve guard above
    # has already handled the worst case, and a poorly-conditioned
    # fit on a noisy small-n curve still yields a usable knee.
    # RankWarning moved to np.exceptions in numpy 2.x; fall back to
    # the legacy attribute name on numpy <2.x installs.
    import warnings as _warnings
    try:
        _RankWarning = np.exceptions.RankWarning  # numpy 2.x+
    except AttributeError:
        _RankWarning = getattr(np, 'RankWarning', Warning)  # numpy <2
    try:
        with _warnings.catch_warnings():
            _warnings.filterwarnings('ignore', category=_RankWarning)
            kl = KneeLocator(
                sizes, aucs,
                curve='concave', direction='increasing',
                S=sensitivity,
                interp_method='polynomial',
            )
    except Exception:
        # Polynomial fitting can fail on degenerate curves; fall back.
        max_auc = aucs.max()
        winners = sizes[aucs == max_auc]
        return int(winners.min()), 'argmax_fallback'

    if kl.knee is None:
        max_auc = aucs.max()
        winners = sizes[aucs == max_auc]
        return int(winners.min()), 'argmax_fallback'

    return int(kl.knee), 'kneedle'


class PanelOptimizer:
    """
    Optimize gene panel size by evaluating classification performance
    at different panel sizes.

    Parameters
    ----------
    random_state : int
        Random seed.
    n_jobs : int
        Number of parallel jobs.
    verbose : bool
        Print progress.
    """

    def __init__(self, random_state: int = DEFAULT_RANDOM_STATE,
                 n_jobs: int = -1, verbose: bool = True):
        if not _SKLEARN_AVAILABLE:
            raise DependencyError("scikit-learn is required for panel optimization")
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

    def _log(self, msg):
        if self.verbose:
            logger.info(msg)

    def forward_selection(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        ranked_genes: List[str],
        min_panel: int = DEFAULT_PANEL_MIN,
        max_panel: int = DEFAULT_PANEL_MAX,
        step: int = 1,
        n_cv: int = 5,
        target_size: Optional[int] = None,
        auto_strategy: str = 'kneedle',
        sensitivity: float = 1.0,
    ) -> PanelOptimizationResult:
        """
        Greedy forward selection: add genes by consensus rank,
        evaluate AUC at each panel size.

        Parameters
        ----------
        X : pd.DataFrame
            Full expression matrix (samples x genes).
        y : np.ndarray
            Labels.
        ranked_genes : List[str]
            Genes ordered by consensus rank (best first).
        min_panel : int
            Minimum panel size to test.
        max_panel : int
            Maximum panel size to test.
        step : int
            Step size for panel sizes to evaluate.
        n_cv : int
            Number of CV folds for AUC estimation.
        target_size : int, optional
            If specified, the panel will be this exact size, bypassing
            auto-detection.
        auto_strategy : {'kneedle', 'argmax', 'first_drop'}
            Auto-detection method when target_size is not specified.
            Default 'kneedle' (Satopaa et al. 2011). See
            _detect_optimal_panel_size for full semantics.
        sensitivity : float
            Kneedle's S parameter (default 1.0); ignored for argmax /
            first_drop.
        """
        from sklearn.model_selection import StratifiedKFold, cross_val_score
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline

        self._log("Panel optimization (forward selection)...")

        # Ensure genes exist in X
        ranked_genes = [g for g in ranked_genes if g in X.columns]
        max_panel = min(max_panel, len(ranked_genes))

        sizes = list(range(min_panel, max_panel + 1, step))
        if sizes and sizes[-1] != max_panel:
            sizes.append(max_panel)

        curve_data = []
        all_panels = {}

        cv = StratifiedKFold(
            n_splits=n_cv, shuffle=True, random_state=self.random_state
        )

        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', RandomForestClassifier(
                n_estimators=100, random_state=self.random_state,
                n_jobs=self.n_jobs,
            ))
        ])

        for size in sizes:
            panel = ranked_genes[:size]
            X_panel = X[panel]

            scores = cross_val_score(
                pipeline, X_panel, y,
                cv=cv, scoring='roc_auc', n_jobs=self.n_jobs,
            )

            curve_data.append({
                'panel_size': size,
                'auc_mean': scores.mean(),
                'auc_std': scores.std(),
            })
            all_panels[size] = panel

        curve_df = pd.DataFrame(curve_data)

        if target_size is not None and target_size in all_panels:
            # User specified exact size
            optimal_size = target_size
            optimal_auc = curve_df[curve_df['panel_size'] == target_size]['auc_mean'].iloc[0]
            selection_method = 'user_specified'
        else:
            optimal_size, selection_method = _detect_optimal_panel_size(
                curve_df,
                auto_strategy=auto_strategy,
                sensitivity=sensitivity,
            )
            # optimal_size from the helper is guaranteed to be in
            # all_panels (drawn directly from curve_df['panel_size']).
            optimal_auc = float(
                curve_df[curve_df['panel_size'] == optimal_size]['auc_mean'].iloc[0]
            )

        optimal_panel = all_panels[optimal_size]

        self._log(
            f"   Optimal panel: {optimal_size} genes, AUC={optimal_auc:.3f} "
            f"(method={selection_method})"
        )

        return PanelOptimizationResult(
            optimal_panel=optimal_panel,
            optimal_size=optimal_size,
            optimal_auc=optimal_auc,
            panel_curve=curve_df,
            all_panels=all_panels,
            method='forward_selection',
            selection_method=selection_method,
        )

    def stability_selection(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        n_bootstrap: int = DEFAULT_N_BOOTSTRAP,
        threshold: float = 0.6,
        alpha_range: Optional[List[float]] = None,
    ) -> pd.DataFrame:
        """
        Stability selection: run LASSO on subsampled data repeatedly,
        report selection frequency per gene.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix.
        y : np.ndarray
            Labels.
        n_bootstrap : int
            Number of bootstrap iterations.
        threshold : float
            Minimum selection frequency to consider a gene stable.

        Returns
        -------
        pd.DataFrame
            Gene stability scores (selection_frequency, is_stable).
        """
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        self._log(f"   Stability selection ({n_bootstrap} iterations)...")

        if alpha_range is None:
            alpha_range = [0.01, 0.05, 0.1, 0.5, 1.0]

        selection_count = np.zeros(X.shape[1])
        total_runs = 0

        rng = np.random.RandomState(self.random_state)

        for i in range(n_bootstrap):
            # Subsample 50% of data
            n = len(X)
            idx = rng.choice(n, size=n // 2, replace=False)
            X_sub = X.iloc[idx]
            y_sub = y[idx]

            scaler = StandardScaler()
            X_sub_s = scaler.fit_transform(X_sub)

            for alpha in alpha_range:
                model = LogisticRegression(
                    C=1.0 / alpha, l1_ratio=1, solver='saga',
                    max_iter=3000, random_state=self.random_state,
                )
                try:
                    model.fit(X_sub_s, y_sub)
                    selected = np.abs(model.coef_).ravel() > 0
                    selection_count += selected.astype(float)
                    total_runs += 1
                except Exception:
                    continue

        freq = selection_count / max(total_runs, 1)

        stability_df = pd.DataFrame({
            'selection_frequency': freq,
            'is_stable': freq >= threshold,
        }, index=X.columns)
        stability_df.index.name = 'gene_id'
        stability_df = stability_df.sort_values(
            'selection_frequency', ascending=False
        )

        n_stable = stability_df['is_stable'].sum()
        self._log(f"   Stability selection: {n_stable} stable genes (freq >= {threshold})")

        return stability_df


# =============================================================================
# 10D: SURVIVAL ANALYSIS (optional)
# =============================================================================

class SurvivalAnalyzer:
    """
    Survival analysis for prognostic biomarker discovery.
    Requires lifelines package.

    Parameters
    ----------
    random_state : int
        Random seed.
    verbose : bool
        Print progress.
    """

    def __init__(self, random_state: int = DEFAULT_RANDOM_STATE,
                 verbose: bool = True):
        if not _LIFELINES_AVAILABLE:
            raise DependencyError(
                "lifelines is required for survival analysis: "
                "pip install lifelines"
            )
        self.random_state = random_state
        self.verbose = verbose

    def _log(self, msg):
        if self.verbose:
            logger.info(msg)

    def cox_univariate_screen(
        self,
        X: pd.DataFrame,
        time_col: np.ndarray,
        event_col: np.ndarray,
        fdr_threshold: float = 0.05,
    ) -> pd.DataFrame:
        """
        Screen genes individually via Cox proportional hazards.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix (samples x genes).
        time_col : np.ndarray
            Survival time.
        event_col : np.ndarray
            Event indicator (1=event, 0=censored).
        fdr_threshold : float
            FDR threshold for significance.

        Returns
        -------
        pd.DataFrame
            Gene-level Cox results: HR, p-value, adjusted_p, CI_lower, CI_upper.
        """
        from lifelines import CoxPHFitter
        from statsmodels.stats.multitest import multipletests

        self._log("   Cox univariate screen...")

        results_list = []

        for gene in X.columns:
            df = pd.DataFrame({
                'T': time_col,
                'E': event_col,
                gene: X[gene].values,
            })

            try:
                cph = CoxPHFitter()
                cph.fit(df, duration_col='T', event_col='E')
                summary = cph.summary

                results_list.append({
                    'gene_id': gene,
                    'hazard_ratio': np.exp(summary['coef'].iloc[0]),
                    'coef': summary['coef'].iloc[0],
                    'p_value': summary['p'].iloc[0],
                    'ci_lower': np.exp(summary['coef lower 95%'].iloc[0]),
                    'ci_upper': np.exp(summary['coef upper 95%'].iloc[0]),
                })
            except Exception:
                results_list.append({
                    'gene_id': gene,
                    'hazard_ratio': np.nan,
                    'coef': np.nan,
                    'p_value': np.nan,
                    'ci_lower': np.nan,
                    'ci_upper': np.nan,
                })

        results_df = pd.DataFrame(results_list)
        results_df = results_df.dropna(subset=['p_value'])

        # FDR correction
        if len(results_df) > 0:
            _, padj, _, _ = multipletests(
                results_df['p_value'], method='fdr_bh'
            )
            results_df['adjusted_p_value'] = padj
        else:
            results_df['adjusted_p_value'] = np.nan

        results_df['is_significant'] = results_df['adjusted_p_value'] < fdr_threshold
        results_df = results_df.sort_values('p_value')

        n_sig = results_df['is_significant'].sum()
        self._log(f"   Cox screen: {n_sig} significant genes (FDR < {fdr_threshold})")

        return results_df

    def cox_lasso_panel(
        self,
        X: pd.DataFrame,
        time_col: np.ndarray,
        event_col: np.ndarray,
        alpha: float = 0.1,
    ) -> SurvivalResult:
        """
        CoxNet (L1-penalized Cox regression) for panel selection.

        Parameters
        ----------
        X : pd.DataFrame
            Expression matrix restricted to candidate genes.
        time_col, event_col : np.ndarray
            Survival data.
        alpha : float
            Regularization strength.

        Returns
        -------
        SurvivalResult
        """
        from lifelines import CoxPHFitter

        self._log("   CoxNet panel selection...")

        df = X.copy()
        df['T'] = time_col
        df['E'] = event_col

        cph = CoxPHFitter(penalizer=alpha, l1_ratio=1.0)
        cph.fit(df, duration_col='T', event_col='E')

        # Non-zero coefficients are the panel
        cox_summary = cph.summary
        selected = cox_summary[cox_summary['coef'].abs() > 1e-6].index.tolist()

        c_index = cph.concordance_index_

        self._log(
            f"   CoxNet: {len(selected)} genes selected, C-index={c_index:.3f}"
        )

        return SurvivalResult(
            significant_genes=selected,
            cox_results=cox_summary,
            c_index=c_index,
            panel_genes=selected,
        )


# =============================================================================
# BIOLOGICAL ANNOTATION & REPORTING ENGINE
# =============================================================================

# Species mapping for STRING API
_STRING_SPECIES_MAP = {
    'human': 9606, 'mouse': 10090, 'rat': 10116,
    'zebrafish': 7955, 'drosophila': 7227, 'c_elegans': 6239,
}


class BiologicalAnnotator:
    """
    Comprehensive biological annotation for biomarker panels.

    Provides gene annotation, pathway enrichment, literature mining,
    and protein-protein interaction context.

    Parameters
    ----------
    species : str
        Species name (human, mouse, rat, etc.).
    verbose : bool
        Print progress messages.
    """

    def __init__(self, species: str = 'human', verbose: bool = True):
        self.species = species
        self.verbose = verbose

    def _log(self, msg):
        if self.verbose:
            logger.info(msg)

    # ---- 1. Gene Annotation (MyGene.info) ----

    def annotate_genes(
        self,
        gene_ids: List[str],
    ) -> pd.DataFrame:
        """
        Query MyGene.info for gene names, descriptions, GO terms,
        and pathway memberships.

        Parameters
        ----------
        gene_ids : List[str]
            Gene identifiers (Ensembl, Symbol, or Entrez).

        Returns
        -------
        pd.DataFrame
            Columns: symbol, name, entrezgene, type_of_gene, summary,
            go_bp (biological process), go_mf (molecular function),
            go_cc (cellular component), pathway_kegg, pathway_reactome.
        """
        try:
            import mygene
        except ImportError:
            self._log(
                "   mygene not installed. Install with: pip install mygene"
            )
            return pd.DataFrame(index=gene_ids)

        self._log(f"   Querying MyGene.info for {len(gene_ids)} genes...")

        mg = mygene.MyGeneInfo()

        try:
            results = mg.querymany(
                gene_ids,
                scopes='ensembl.gene,symbol,entrezgene',
                fields='symbol,name,entrezgene,type_of_gene,summary,'
                       'go.BP,go.MF,go.CC,'
                       'pathway.kegg,pathway.reactome',
                species=self.species,
                returnall=True,
            )
        except Exception as e:
            self._log(f"   MyGene.info query failed: {e}")
            return pd.DataFrame(index=gene_ids)

        rows = []
        for hit in results.get('out', []):
            row = {
                'query': hit.get('query', ''),
                'symbol': hit.get('symbol', ''),
                'name': hit.get('name', ''),
                'entrezgene': str(hit.get('entrezgene', '')),
                'type_of_gene': hit.get('type_of_gene', ''),
                'summary': hit.get('summary', ''),
            }

            # GO terms — extract term names
            for go_cat, go_key in [('go_bp', 'BP'), ('go_mf', 'MF'), ('go_cc', 'CC')]:
                go_data = hit.get('go', {})
                if isinstance(go_data, dict):
                    terms = go_data.get(go_key, [])
                    if isinstance(terms, list):
                        row[go_cat] = '; '.join(
                            t.get('term', '') for t in terms[:10]
                            if isinstance(t, dict)
                        )
                    elif isinstance(terms, dict):
                        row[go_cat] = terms.get('term', '')
                    else:
                        row[go_cat] = ''
                else:
                    row[go_cat] = ''

            # Pathways
            pw = hit.get('pathway', {})
            if isinstance(pw, dict):
                kegg = pw.get('kegg', [])
                if isinstance(kegg, list):
                    row['pathway_kegg'] = '; '.join(
                        p.get('name', '') for p in kegg if isinstance(p, dict)
                    )
                elif isinstance(kegg, dict):
                    row['pathway_kegg'] = kegg.get('name', '')
                else:
                    row['pathway_kegg'] = ''

                reactome = pw.get('reactome', [])
                if isinstance(reactome, list):
                    row['pathway_reactome'] = '; '.join(
                        p.get('name', '') for p in reactome if isinstance(p, dict)
                    )
                elif isinstance(reactome, dict):
                    row['pathway_reactome'] = reactome.get('name', '')
                else:
                    row['pathway_reactome'] = ''
            else:
                row['pathway_kegg'] = ''
                row['pathway_reactome'] = ''

            rows.append(row)

        ann_df = pd.DataFrame(rows)
        if not ann_df.empty:
            ann_df = ann_df.drop_duplicates(subset='query', keep='first')
            ann_df = ann_df.set_index('query')
            ann_df.index.name = 'gene_id'

        n_found = (ann_df['symbol'] != '').sum() if not ann_df.empty else 0
        self._log(f"   Annotated {n_found}/{len(gene_ids)} genes")

        return ann_df

    # ---- 2. Pathway Enrichment (Fisher's exact test ORA) ----

    def pathway_enrichment(
        self,
        gene_ids: List[str],
        background_genes: Optional[List[str]] = None,
        gene_sets: Optional[Dict[str, List[str]]] = None,
        databases: Optional[List[str]] = None,
        fdr_threshold: float = 0.05,
    ) -> pd.DataFrame:
        """
        Over-representation analysis using Fisher's exact test.

        If gseapy is installed, uses its enrich() function with built-in
        gene set libraries (KEGG, Reactome, GO). Otherwise, falls back
        to a manual Fisher's exact test if custom gene_sets are provided.

        Parameters
        ----------
        gene_ids : List[str]
            Biomarker gene set (symbols preferred).
        background_genes : List[str], optional
            Background gene universe. If None, uses the gene set library
            universe (gseapy) or requires gene_sets with explicit background.
        gene_sets : Dict[str, List[str]], optional
            Custom gene sets {pathway_name: [gene1, gene2, ...]}.
            Used when gseapy is not available.
        databases : List[str], optional
            Gene set databases to query. Default:
            ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023'].
        fdr_threshold : float
            FDR threshold for reporting.

        Returns
        -------
        pd.DataFrame
            Columns: pathway, source, p_value, adjusted_p_value,
            overlap_genes, overlap_count, gene_set_size, odds_ratio.
        """
        self._log("   Running pathway enrichment (ORA)...")

        if databases is None:
            databases = [
                'KEGG_2021_Human',
                'Reactome_2022',
                'GO_Biological_Process_2023',
            ]

        # Try gseapy first
        try:
            import gseapy as gp

            all_results = []

            for db in databases:
                try:
                    enr = gp.enrich(
                        gene_list=gene_ids,
                        gene_sets=db,
                        background=background_genes,
                        outdir=None,
                        no_plot=True,
                        verbose=False,
                    )

                    if enr.results is not None and not enr.results.empty:
                        df = enr.results.copy()
                        df['source'] = db
                        all_results.append(df)

                except Exception as e:
                    self._log(f"   Warning: {db} enrichment failed: {e}")
                    continue

            if all_results:
                combined = pd.concat(all_results, ignore_index=True)

                # Standardize column names
                col_map = {
                    'Term': 'pathway',
                    'P-value': 'p_value',
                    'Adjusted P-value': 'adjusted_p_value',
                    'Overlap': 'overlap_info',
                    'Genes': 'overlap_genes',
                    'Odds Ratio': 'odds_ratio',
                }
                combined = combined.rename(
                    columns={k: v for k, v in col_map.items()
                             if k in combined.columns}
                )

                # Parse overlap count
                if 'overlap_info' in combined.columns:
                    combined['overlap_count'] = combined['overlap_info'].apply(
                        lambda x: int(str(x).split('/')[0])
                        if '/' in str(x) else 0
                    )
                    combined['gene_set_size'] = combined['overlap_info'].apply(
                        lambda x: int(str(x).split('/')[1])
                        if '/' in str(x) else 0
                    )

                # Sort by p-value
                combined = combined.sort_values('p_value')

                n_sig = (combined['adjusted_p_value'] < fdr_threshold).sum() \
                    if 'adjusted_p_value' in combined.columns else 0
                self._log(
                    f"   Found {n_sig} enriched pathways (FDR < {fdr_threshold})"
                )

                return combined

        except ImportError:
            self._log("   gseapy not installed, trying manual ORA...")

        # Fallback: manual Fisher's exact test with custom gene_sets
        if gene_sets is None:
            self._log(
                "   No gene sets provided and gseapy not available. "
                "Install gseapy: pip install gseapy"
            )
            return pd.DataFrame()

        return self._manual_ora(
            gene_ids, gene_sets, background_genes, fdr_threshold,
        )

    def _manual_ora(
        self,
        gene_ids: List[str],
        gene_sets: Dict[str, List[str]],
        background_genes: Optional[List[str]],
        fdr_threshold: float,
    ) -> pd.DataFrame:
        """Manual over-representation analysis using Fisher's exact test."""
        from scipy.stats import fisher_exact
        from statsmodels.stats.multitest import multipletests

        query_set = set(gene_ids)

        if background_genes:
            universe = set(background_genes) | query_set
        else:
            all_gs_genes = set()
            for genes in gene_sets.values():
                all_gs_genes.update(genes)
            universe = all_gs_genes | query_set

        n_universe = len(universe)
        n_query = len(query_set & universe)

        results_list = []

        for pathway_name, pathway_genes in gene_sets.items():
            pathway_set = set(pathway_genes) & universe
            n_pathway = len(pathway_set)

            if n_pathway < 2:
                continue

            overlap = query_set & pathway_set
            n_overlap = len(overlap)

            if n_overlap == 0:
                continue

            # 2x2 contingency table
            a = n_overlap
            b = n_query - n_overlap
            c = n_pathway - n_overlap
            d = n_universe - n_query - n_pathway + n_overlap

            _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

            odds_ratio = (a * d) / max(b * c, 1)

            results_list.append({
                'pathway': pathway_name,
                'source': 'custom',
                'p_value': p_value,
                'overlap_genes': '; '.join(sorted(overlap)),
                'overlap_count': n_overlap,
                'gene_set_size': n_pathway,
                'odds_ratio': odds_ratio,
            })

        if not results_list:
            return pd.DataFrame()

        results_df = pd.DataFrame(results_list)

        # FDR correction
        _, padj, _, _ = multipletests(
            results_df['p_value'], method='fdr_bh',
        )
        results_df['adjusted_p_value'] = padj
        results_df = results_df.sort_values('p_value')

        n_sig = (results_df['adjusted_p_value'] < fdr_threshold).sum()
        self._log(f"   Manual ORA: {n_sig} enriched gene sets (FDR < {fdr_threshold})")

        return results_df

    # ---- 3. Literature Mining (NCBI Entrez / Europe PMC) ----

    def literature_search(
        self,
        gene_ids: List[str],
        disease_term: Optional[str] = None,
        max_per_gene: int = 5,
    ) -> pd.DataFrame:
        """
        Search PubMed/Europe PMC for publications associated with
        each biomarker gene.

        Uses the Europe PMC REST API (no authentication needed).

        Parameters
        ----------
        gene_ids : List[str]
            Gene symbols or IDs.
        disease_term : str, optional
            Disease context to add to queries (e.g., 'breast cancer').
        max_per_gene : int
            Maximum publications to retrieve per gene.

        Returns
        -------
        pd.DataFrame
            Columns: gene_id, n_publications, top_titles, pmids.
        """
        import urllib.request
        import urllib.parse

        self._log(f"   Literature search for {len(gene_ids)} genes...")

        base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

        results_list = []

        for gene in gene_ids:
            query = f'"{gene}" AND (RNA-seq OR transcriptomics OR biomarker)'
            if disease_term:
                query += f' AND "{disease_term}"'

            params = urllib.parse.urlencode({
                'query': query,
                'format': 'json',
                'pageSize': max_per_gene,
                'resultType': 'lite',
            })
            url = f"{base_url}?{params}"

            try:
                req = urllib.request.Request(
                    url,
                    headers={'User-Agent': 'RAPTOR/2.2.2 (biomarker discovery)'},
                )
                with urllib.request.urlopen(req, timeout=10) as response:
                    data = json.loads(response.read().decode('utf-8'))

                result_list = data.get('resultList', {}).get('result', [])
                n_total = data.get('hitCount', 0)

                titles = []
                pmids = []
                for r in result_list[:max_per_gene]:
                    titles.append(r.get('title', ''))
                    pmid = r.get('pmid', '')
                    if pmid:
                        pmids.append(str(pmid))

                results_list.append({
                    'gene_id': gene,
                    'n_publications': n_total,
                    'top_titles': ' | '.join(titles[:3]),
                    'pmids': ', '.join(pmids),
                })

            except Exception:
                results_list.append({
                    'gene_id': gene,
                    'n_publications': 0,
                    'top_titles': '',
                    'pmids': '',
                })

        results_df = pd.DataFrame(results_list)
        n_with_hits = (results_df['n_publications'] > 0).sum()
        self._log(f"   Literature: {n_with_hits}/{len(gene_ids)} genes have publications")

        return results_df

    # ---- 4. PPI Network (STRING API) ----

    def ppi_network(
        self,
        gene_ids: List[str],
        score_threshold: int = 400,
    ) -> Dict:
        """
        Query STRING database for protein-protein interactions.

        Uses the STRING REST API (no authentication for small queries).

        Parameters
        ----------
        gene_ids : List[str]
            Gene symbols.
        score_threshold : int
            Minimum interaction confidence score (0-1000).
            400 = medium confidence, 700 = high, 900 = highest.

        Returns
        -------
        Dict
            Keys: nodes, edges, n_nodes, n_edges, enrichment_pvalue,
            network_url (link to STRING visualization).
        """
        import urllib.request
        import urllib.parse

        self._log(f"   Querying STRING PPI for {len(gene_ids)} genes...")

        species_id = _STRING_SPECIES_MAP.get(self.species, 9606)

        # Get interactions
        base_url = "https://string-db.org/api/json"

        identifiers = '%0d'.join(gene_ids[:200])  # STRING limits to ~200

        params = urllib.parse.urlencode({
            'identifiers': identifiers,
            'species': species_id,
            'required_score': score_threshold,
            'caller_identity': 'RAPTOR_v2.2.2',
        })

        network_result = {
            'nodes': gene_ids,
            'edges': [],
            'n_nodes': len(gene_ids),
            'n_edges': 0,
            'enrichment_pvalue': None,
            'network_url': None,
        }

        # Fetch interactions
        try:
            url = f"{base_url}/network?{params}"
            req = urllib.request.Request(
                url,
                headers={'User-Agent': 'RAPTOR/2.2.2'},
            )
            with urllib.request.urlopen(req, timeout=15) as response:
                edges_data = json.loads(response.read().decode('utf-8'))

            edges = []
            for edge in edges_data:
                edges.append({
                    'protein_a': edge.get('preferredName_A', ''),
                    'protein_b': edge.get('preferredName_B', ''),
                    'score': edge.get('score', 0),
                })

            network_result['edges'] = edges
            network_result['n_edges'] = len(edges)

        except Exception as e:
            self._log(f"   STRING network query failed: {e}")

        # Fetch PPI enrichment p-value
        try:
            url = f"{base_url}/ppi_enrichment?{params}"
            req = urllib.request.Request(
                url,
                headers={'User-Agent': 'RAPTOR/2.2.2'},
            )
            with urllib.request.urlopen(req, timeout=10) as response:
                enr_data = json.loads(response.read().decode('utf-8'))

            if enr_data:
                network_result['enrichment_pvalue'] = enr_data[0].get(
                    'p_value', None
                )

        except Exception:
            pass

        # Build visualization URL
        try:
            gene_str = '%0d'.join(gene_ids[:50])
            network_result['network_url'] = (
                f"https://string-db.org/cgi/network?"
                f"identifiers={gene_str}&species={species_id}"
            )
        except Exception:
            pass

        n_edges = network_result['n_edges']
        ppi_p = network_result.get('enrichment_pvalue')
        ppi_msg = f"   STRING: {n_edges} interactions"
        if ppi_p is not None:
            ppi_msg += f", PPI enrichment p={ppi_p:.2e}"
        self._log(ppi_msg)

        return network_result

    # ---- 5. Report Generation ----

    def generate_report(
        self,
        panel_genes: List[str],
        annotation_result: 'AnnotationResult',
        classification_results: Optional[Dict[str, ClassificationResult]] = None,
        panel_optimization: Optional[PanelOptimizationResult] = None,
        output_path: Optional[Path] = None,
    ) -> str:
        """
        Generate a structured publication-ready biomarker report
        in Markdown format.

        Parameters
        ----------
        panel_genes : List[str]
            Final biomarker panel.
        annotation_result : AnnotationResult
            Biological annotation results.
        classification_results : Dict, optional
            Classifier performance metrics.
        panel_optimization : PanelOptimizationResult, optional
            Panel optimization results.
        output_path : Path, optional
            If provided, save report to this file.

        Returns
        -------
        str
            Markdown-formatted report.
        """
        self._log("   Generating biomarker report...")

        sections = []

        # Header
        sections.append("# RAPTOR Biomarker Discovery Report")
        sections.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
        sections.append(f"RAPTOR version: {__raptor_version__}")
        sections.append(f"Species: {self.species}")

        # Panel summary
        sections.append("\n## Biomarker Panel\n")
        sections.append(f"**Panel size:** {len(panel_genes)} genes\n")

        if not annotation_result.gene_annotations.empty:
            ann = annotation_result.gene_annotations
            sections.append("| Gene ID | Symbol | Name | Type |")
            sections.append("|---------|--------|------|------|")
            for gene in panel_genes:
                if gene in ann.index:
                    row = ann.loc[gene]
                    sections.append(
                        f"| {gene} | {row.get('symbol', '')} | "
                        f"{row.get('name', '')[:60]} | "
                        f"{row.get('type_of_gene', '')} |"
                    )
                else:
                    sections.append(f"| {gene} | - | - | - |")

        # Classification performance
        if classification_results:
            sections.append("\n## Classification Performance\n")
            sections.append("| Classifier | AUC | F1 | Sensitivity | Specificity |")
            sections.append("|------------|-----|-----|-------------|-------------|")
            for name, res in classification_results.items():
                sections.append(
                    f"| {name} | {res.auc:.3f} | {res.f1:.3f} | "
                    f"{res.sensitivity:.3f} | {res.specificity:.3f} |"
                )

        # Panel optimization
        if panel_optimization:
            sections.append("\n## Panel Optimization\n")
            sections.append(
                f"Optimal panel size: **{panel_optimization.optimal_size}** genes "
                f"(AUC = {panel_optimization.optimal_auc:.3f})"
            )

        # Pathway enrichment
        if not annotation_result.pathway_enrichment.empty:
            enr = annotation_result.pathway_enrichment
            sig_enr = enr[
                enr.get('adjusted_p_value', enr.get('p_value', pd.Series())) < 0.05
            ] if 'adjusted_p_value' in enr.columns else enr.head(10)

            if not sig_enr.empty:
                sections.append("\n## Enriched Pathways (FDR < 0.05)\n")
                sections.append("| Pathway | Source | P-value | Genes |")
                sections.append("|---------|--------|---------|-------|")
                for _, row in sig_enr.head(20).iterrows():
                    pval = row.get('adjusted_p_value', row.get('p_value', ''))
                    pval_str = f"{pval:.2e}" if isinstance(pval, float) else str(pval)
                    sections.append(
                        f"| {row.get('pathway', '')[:50]} | "
                        f"{row.get('source', '')} | "
                        f"{pval_str} | "
                        f"{row.get('overlap_count', '')} |"
                    )

        # Literature context
        if not annotation_result.literature_hits.empty:
            lit = annotation_result.literature_hits
            lit_with_hits = lit[lit['n_publications'] > 0]

            if not lit_with_hits.empty:
                sections.append("\n## Literature Context\n")
                sections.append(
                    f"{len(lit_with_hits)}/{len(lit)} panel genes have "
                    f"existing publications in PubMed."
                )
                sections.append("\n| Gene | Publications | Top Reference |")
                sections.append("|------|-------------|---------------|")
                for _, row in lit_with_hits.head(20).iterrows():
                    top_title = row['top_titles'].split(' | ')[0][:60] \
                        if row['top_titles'] else '-'
                    sections.append(
                        f"| {row['gene_id']} | {row['n_publications']} | "
                        f"{top_title} |"
                    )

        # PPI network
        if annotation_result.ppi_network is not None:
            ppi = annotation_result.ppi_network
            sections.append("\n## Protein-Protein Interactions\n")
            sections.append(
                f"STRING network: **{ppi['n_edges']}** interactions "
                f"among {ppi['n_nodes']} proteins"
            )
            if ppi.get('enrichment_pvalue') is not None:
                sections.append(
                    f"\nPPI enrichment p-value: {ppi['enrichment_pvalue']:.2e}"
                )
                if ppi['enrichment_pvalue'] < 0.05:
                    sections.append(
                        "The panel genes have significantly more interactions "
                        "than expected, suggesting functional coherence."
                    )
            if ppi.get('network_url'):
                sections.append(
                    f"\n[View interactive network on STRING]({ppi['network_url']})"
                )

        report_text = '\n'.join(sections)

        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(report_text)
            self._log(f"   Report saved: {output_path}")

        return report_text

    # ---- Full annotation pipeline ----

    def annotate_panel(
        self,
        panel_genes: List[str],
        background_genes: Optional[List[str]] = None,
        disease_term: Optional[str] = None,
        run_literature: bool = True,
        run_ppi: bool = True,
    ) -> AnnotationResult:
        """
        Run the complete annotation pipeline on a biomarker panel.

        Parameters
        ----------
        panel_genes : List[str]
            Gene panel to annotate.
        background_genes : List[str], optional
            Background gene universe for enrichment.
        disease_term : str, optional
            Disease context for literature search.
        run_literature : bool
            Whether to query PubMed/Europe PMC.
        run_ppi : bool
            Whether to query STRING for PPI.

        Returns
        -------
        AnnotationResult
        """
        self._log("Running biological annotation pipeline...")

        # 1. Gene annotation
        gene_ann = self.annotate_genes(panel_genes)

        # Use symbols for enrichment/literature if available
        symbols = []
        if not gene_ann.empty and 'symbol' in gene_ann.columns:
            for gene in panel_genes:
                if gene in gene_ann.index and gene_ann.loc[gene, 'symbol']:
                    symbols.append(gene_ann.loc[gene, 'symbol'])
                else:
                    symbols.append(gene)
        else:
            symbols = panel_genes

        # 2. Pathway enrichment
        enrichment_df = self.pathway_enrichment(
            symbols,
            background_genes=background_genes,
        )

        # 3. Literature mining
        lit_df = pd.DataFrame()
        if run_literature:
            lit_df = self.literature_search(
                symbols, disease_term=disease_term,
            )

        # 4. PPI network
        ppi_data = None
        if run_ppi:
            ppi_data = self.ppi_network(symbols)

        result = AnnotationResult(
            gene_annotations=gene_ann,
            pathway_enrichment=enrichment_df,
            literature_hits=lit_df,
            ppi_network=ppi_data,
            background_size=len(background_genes) if background_genes else 0,
            species=self.species,
        )

        self._log("   Annotation pipeline complete")
        return result


# =============================================================================
# MAIN CONVENIENCE FUNCTIONS (Public API)
# =============================================================================

def _prepare_expression_data(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    de_genes: Optional[List[str]] = None,
    min_count: int = 10,
    baseline_group: Optional[str] = None,
    reference_group: Optional[str] = None,
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Prepare expression matrix and labels for ML.

    - Filters low-count genes
    - Log2-CPM normalizes
    - Encodes labels as 0/1
    - Optionally restricts to DE genes

    Parameters
    ----------
    baseline_group : str, optional
        The group to encode as label 0 (control / untreated / healthy).
        If None, falls back to sorted-alphabetical order (legacy behaviour).
    reference_group : str, optional
        The group to encode as label 1 (case / treated / disease).
        If None, falls back to sorted-alphabetical order (legacy behaviour).

    Returns
    -------
    X : pd.DataFrame (samples x genes)
    y : np.ndarray (binary labels)
    gene_list : List[str] (gene order)
    """
    counts = counts.copy()

    # Filter low-count genes
    gene_means = counts.mean(axis=1)
    keep_genes = gene_means[gene_means >= min_count].index
    counts = counts.loc[keep_genes]

    # If DE genes provided, restrict to those that pass filter
    if de_genes is not None:
        de_in_data = [g for g in de_genes if g in counts.index]
        if len(de_in_data) < 10:
            logger.warning(
                f"Only {len(de_in_data)} DE genes in expression data. "
                f"Using all {len(counts)} filtered genes."
            )
        else:
            counts = counts.loc[de_in_data]

    # Log2-CPM normalization
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # Align samples between counts and metadata
    common_samples = list(
        set(log_cpm.columns) & set(metadata.iloc[:, 0].values
                                    if metadata.columns[0] != group_column
                                    else metadata.index.astype(str))
    )

    # Try to match by column name in metadata
    if 'sample_id' in metadata.columns:
        sample_col = 'sample_id'
    elif metadata.columns[0] != group_column:
        sample_col = metadata.columns[0]
    else:
        sample_col = None

    if sample_col:
        meta_samples = metadata[sample_col].astype(str).values
        count_samples = [str(s) for s in log_cpm.columns]
        common_samples = [s for s in count_samples if s in meta_samples]

        if not common_samples:
            # Try matching index
            common_samples = [
                str(s) for s in log_cpm.columns
                if str(s) in metadata.index.astype(str)
            ]
            if common_samples:
                sample_col = None  # Use index matching

    if not common_samples:
        raise ValidationError(
            'samples',
            "No matching samples between counts and metadata",
            "Ensure sample IDs match between files"
        )

    # Build aligned X and y
    if sample_col:
        meta_aligned = metadata.set_index(sample_col).loc[common_samples]
    else:
        meta_aligned = metadata.loc[common_samples]

    X = log_cpm[common_samples].T  # samples x genes
    X.index = common_samples

    # Encode labels
    groups = meta_aligned[group_column].values
    unique_groups = sorted(set(groups))

    if len(unique_groups) != 2:
        raise ValidationError(
            'group_column',
            f"Expected 2 groups for binary classification, got {len(unique_groups)}: {unique_groups}",
            "Multi-class support will be added in a future version"
        )

    label_map: Dict[str, int]
    if baseline_group is not None and reference_group is not None:
        # Respect the caller's explicit group assignment
        if baseline_group not in unique_groups:
            raise ValidationError(
                'baseline_group',
                f"baseline_group={baseline_group!r} not found in "
                f"metadata[{group_column!r}] values {unique_groups}",
                "Check the baseline_group argument matches a real group name",
            )
        if reference_group not in unique_groups:
            raise ValidationError(
                'reference_group',
                f"reference_group={reference_group!r} not found in "
                f"metadata[{group_column!r}] values {unique_groups}",
                "Check the reference_group argument matches a real group name",
            )
        if baseline_group == reference_group:
            raise ValidationError(
                'reference_group',
                f"baseline_group and reference_group must differ "
                f"(both were {baseline_group!r})",
                "Pick two distinct groups to compare",
            )
        label_map = {baseline_group: 0, reference_group: 1}
    elif baseline_group is not None or reference_group is not None:
        # Partial specification — treat as a programming error. Silently
        # falling back to alphabetical here would re-introduce the exact
        # polarity-flip bug the explicit kwargs are meant to prevent.
        raise ValidationError(
            'baseline_group/reference_group',
            "Either both or neither of baseline_group and reference_group "
            f"must be supplied (got baseline={baseline_group!r}, "
            f"reference={reference_group!r})",
            "Pass both group names together, or omit both to use the "
            "alphabetical-sort legacy default",
        )
    else:
        # Legacy fallback: alphabetical sort. This preserves backward
        # compatibility for callers who don't pass explicit groups.
        label_map = {unique_groups[0]: 0, unique_groups[1]: 1}

    y = np.array([label_map[g] for g in groups])

    # Log the mapping so it's in the pipeline output — makes label polarity
    # auditable from logs alone if a downstream result looks inverted.
    inv = {v: k for k, v in label_map.items()}
    logger.info(
        f"   Label encoding: 0={inv[0]!r} (baseline), 1={inv[1]!r} (reference)"
    )
    logger.info(
        f"   Prepared: {X.shape[0]} samples x {X.shape[1]} genes, "
        f"groups: {inv[0]}={int((y==0).sum())}, "
        f"{inv[1]}={int((y==1).sum())}"
    )

    return X, y, list(X.columns)


# =============================================================================
# M6: PIPELINE-CV ORCHESTRATOR (feature selection inside outer CV loop)
# =============================================================================
#
# The central theorem motivating M6 is Ambroise & McLachlan (2002, PNAS
# 99:6562): when feature selection uses outcome labels, cross-validation
# must be external to the selection process to produce an unbiased
# estimate of out-of-sample error. Varma & Simon (2006, BMC Bioinf 7:91)
# generalized this to any model selection step. Hastie, Tibshirani &
# Friedman (ESL 2nd ed., section 7.10.2) present this as the canonical
# "wrong way vs right way" pair.
#
# The fix is to run every step that reads y inside the outer training
# fold:
#     Stage 1 (feature selection) -- inside the fold
#     Stage 2 (panel optimization) -- inside the fold
#     Stage 3 (classifier fit + predict on held-out fold) -- inside the fold
#
# A single "final reported panel" is produced by re-running Stages 1+2
# on all data. This matches the nestedcv R package (Lewis et al. 2023)
# convention: the final deployable model uses all available samples,
# while the OOF AUC honestly estimates how well the procedure
# generalizes.
#
# Panel stability across folds (Nogueira et al. 2018 JMLR) is reported
# as an additional diagnostic. See module-level helpers at the top of
# this file for the stability math.

# The per-fold-safe default methods. Notably this REPLACES 'de_filter'
# with 'univariate_filter' when no de_genes list was supplied — because
# univariate_filter is the per-fold-safe equivalent (Lewis et al. 2023
# nestedcv; Haury et al. 2011 PLOS ONE).
_FOLD_SAFE_DEFAULT_METHODS = (
    'univariate_filter',
    'elastic_net',
    'rfe',
)


def _run_feature_selection(
    selector: 'FeatureSelector',
    methods: List[str],
    X: pd.DataFrame,
    y: np.ndarray,
    gene_list: List[str],
    de_genes: Optional[List[str]],
    species: str,
    n_select: int,
) -> None:
    """
    Populate ``selector.results`` by dispatching to each method in
    ``methods``. Errors in any one method are logged and skipped.

    Shared by the per-fold loop in ``_run_pipeline_cv`` and the final
    all-data pass in ``discover_biomarkers`` — both paths need to run
    the same method suite and produce the same ``FeatureSelectionResult``
    shape for downstream consensus ranking.
    """
    _wgcna_species_map = {
        'human': 'homo sapiens',
        'mouse': 'mus musculus',
        'rat': 'rattus norvegicus',
    }
    n_samples = X.shape[0]

    for method in methods:
        try:
            if method == 'de_filter' and de_genes:
                selector.select_de_filter(de_genes, gene_list)
            elif method == 'univariate_filter':
                selector.select_univariate_filter(
                    X, y, n_features=n_select, test='welch',
                )
            elif method in ('lasso', 'elastic_net'):
                l1 = 1.0 if method == 'lasso' else 0.5
                selector.select_lasso(X, y, l1_ratio=l1, label=method)
            elif method == 'boruta':
                selector.select_boruta(X, y)
            elif method == 'mrmr':
                selector.select_mrmr(X, y, n_features=n_select)
            elif method == 'rfe':
                selector.select_rfe(X, y, n_features=n_select)
            elif method == 'shap':
                selector.select_shap(X, y, n_features=n_select)
            elif method == 'wgcna':
                # WGCNA needs a minimum sample count. Skip rather than
                # crash when a training fold falls below it.
                if n_samples < 15:
                    logger.warning(
                        f"   Skipping wgcna: training fold has "
                        f"{n_samples} samples (<15 required). "
                        f"This is expected inside small-cohort CV folds."
                    )
                    continue
                wgcna_sp = _wgcna_species_map.get(species, species)
                selector.select_wgcna(X, y, species=wgcna_sp)
        except DependencyError as e:
            logger.warning(f"   Skipping {method}: {e}")
        except Exception as e:
            logger.warning(f"   Error in {method}: {e}")


def _resolve_pipeline_cv_methods(
    methods: Optional[List[str]],
    de_genes: Optional[List[str]],
) -> List[str]:
    """
    Return the method list to use per-fold, substituting the fold-safe
    ``univariate_filter`` for ``de_filter`` when no independent
    ``de_genes`` list was supplied.

    When ``de_genes`` IS supplied, we keep ``de_filter`` but emit a
    prominent warning about upstream leakage: if that gene list came
    from DE analysis on the same cohort, M10's pipeline-CV cannot
    compensate for the bias already baked into the gene list.
    """
    if methods is None:
        methods = list(_FOLD_SAFE_DEFAULT_METHODS)
        if _BORUTA_AVAILABLE:
            methods.append('boruta')
        if _MRMR_AVAILABLE:
            methods.append('mrmr')
        if _SHAP_AVAILABLE:
            methods.append('shap')

    resolved = list(methods)

    has_de_filter = 'de_filter' in resolved
    has_univariate = 'univariate_filter' in resolved

    if de_genes is not None and has_de_filter:
        # User explicitly wants to filter on precomputed de_genes. Honor
        # the request but flag the leakage risk.
        warnings.warn(
            "M10 pipeline-CV: 'de_filter' is being used with a precomputed "
            "de_genes list. If that list came from DE analysis on the same "
            "samples being evaluated here, the feature-selection step has "
            "already seen every label and the resulting CV estimate will "
            "still be optimistically biased -- M10 cannot fix leakage that "
            "happens upstream of M10. For fully fold-safe analysis, drop "
            "'de_filter' and use 'univariate_filter' instead, which "
            "recomputes per-fold on training data only "
            "(see Lewis et al. 2023 nestedcv; Haury et al. 2011 PLOS ONE).",
            UserWarning,
            stacklevel=3,
        )
    elif de_genes is None and has_de_filter and not has_univariate:
        # User requested de_filter but didn't provide de_genes; swap in
        # the fold-safe replacement.
        logger.info(
            "   M10 pipeline-CV: 'de_filter' requested without de_genes; "
            "using 'univariate_filter' (per-fold Welch's t-test) as the "
            "fold-safe replacement."
        )
        resolved = [m if m != 'de_filter' else 'univariate_filter' for m in resolved]

    return resolved


def _run_pipeline_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    gene_list: List[str],
    methods: List[str],
    de_genes: Optional[List[str]],
    species: str,
    n_outer: int = 5,
    n_repeats: int = 1,
    n_select: int = 50,
    min_panel: int = 5,
    max_panel: int = 50,
    target_panel_size: Optional[int] = None,
    random_state: int = 42,
    verbose: bool = True,
    classifiers: Optional[List[str]] = None,
    auto_panel_strategy: str = 'kneedle',
    panel_sensitivity: float = 1.0,
) -> Tuple[
    Dict[str, 'ClassificationResult'],
    List[List[str]],
    List[pd.DataFrame],
]:
    """
    Pipeline-CV orchestrator: runs feature selection + panel optimization
    INSIDE the outer CV loop, and evaluates classifiers on the held-out
    test fold.

    This implements the fix recommended by Ambroise & McLachlan (2002):
    every step that uses y is restricted to the training fold, so the
    OOF predictions stored on the returned ``ClassificationResult``
    objects are honest estimates of out-of-sample performance.

    Parameters
    ----------
    X : pd.DataFrame (samples x genes)
        Full expression matrix.
    y : np.ndarray
        Binary labels.
    gene_list : list of str
        All candidate gene names (X.columns).
    methods : list of str
        Feature-selection methods to run per fold. Should already have
        been passed through ``_resolve_pipeline_cv_methods`` so that
        'de_filter' vs 'univariate_filter' substitution has occurred.
    de_genes : list of str, optional
        Precomputed DE-significant genes (for 'de_filter' method).
    species : str
        Species name for WGCNA.
    n_outer : int, default 5
        Number of outer folds. Ambroise & McLachlan (2002) explicitly
        recommend 10-fold over LOOCV; 5 matches RAPTOR's current
        DEFAULT_N_FOLDS_OUTER and is used as the default to minimize
        behavioral change from pre-M6.
    n_repeats : int, default 1
        Number of outer-CV repeats (repeated k-fold). Per Bengio &
        Grandvalet (2004), repeated k-fold reduces the variance of the
        OOF AUC estimate. Default 1 preserves runtime; paper-grade
        results use 3-5.
    n_select : int, default 50
        Target number of features per selection method.
    min_panel, max_panel : int
        Bounds for the panel-size search.
    target_panel_size : int, optional
        If set, forces the exact panel size.
    random_state : int
        Seed. Each fold uses random_state + fold_idx + repeat_idx * n_outer
        so method-internal randomness (Boruta, RF, etc.) differs across
        folds but the whole run is reproducible.
    verbose : bool
        If True, log fold-by-fold progress.
    classifiers : list of str, optional
        Restrict to these classifier names. None = all available.

    Returns
    -------
    clf_results : dict mapping classifier name to ClassificationResult
        Each result has:
        * oof_true, oof_prob: honest held-out predictions (one per
          sample per repeat; len = n_samples * n_repeats)
        * accuracy, auc, f1, sensitivity, specificity: mean across folds
        * metrics_per_fold: per-fold dict with fold/repeat indices
        * roc_data: pooled ROC from all OOF predictions
        * trained_model: None (populated later from final all-data fit)
        * feature_importance: None (populated later)
    per_fold_panels : list of lists
        One gene list per (fold, repeat). Length = n_outer * n_repeats.
    per_fold_ranked_genes : list of pd.DataFrame
        One consensus ranking per fold.
    """
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import (
        accuracy_score, roc_auc_score, f1_score,
        confusion_matrix, roc_curve,
    )
    from sklearn.preprocessing import StandardScaler
    from sklearn.base import clone

    # Build the classifier pool once so the dispatch is identical to the
    # pre-M6 ClassifierEvaluator path.
    evaluator = ClassifierEvaluator(
        random_state=random_state, verbose=False,
    )
    all_clfs = evaluator._get_classifiers()
    if classifiers:
        all_clfs = {k: v for k, v in all_clfs.items() if k in classifiers}

    # Accumulators
    oof_true_all: List[int] = []
    oof_prob_by_clf: Dict[str, List[float]] = {n: [] for n in all_clfs}
    per_fold_metrics: Dict[str, List[Dict]] = {n: [] for n in all_clfs}
    per_fold_panels: List[List[str]] = []
    per_fold_ranked_genes: List[pd.DataFrame] = []

    for repeat_idx in range(n_repeats):
        repeat_seed = random_state + repeat_idx * 1000
        outer_cv = StratifiedKFold(
            n_splits=n_outer, shuffle=True, random_state=repeat_seed,
        )

        for fold_i, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
            fold_seed = repeat_seed + fold_i

            X_tr = X.iloc[train_idx]
            X_te = X.iloc[test_idx]
            y_tr = y[train_idx]
            y_te = y[test_idx]

            if verbose:
                logger.info(
                    f"   [pipeline-CV] repeat {repeat_idx+1}/{n_repeats}, "
                    f"fold {fold_i+1}/{n_outer} "
                    f"(n_train={len(train_idx)}, n_test={len(test_idx)})"
                )

            # --- Stage 1 on training fold only ---
            # n_jobs=1 inside folds: outer-fold parallelism is the right
            # level; spawning inner worker processes on tiny training
            # folds (<100 samples) is pure overhead (joblib process-fork
            # dominates compute). Users who want outer-fold parallelism
            # can enable it at a future _run_pipeline_cv parameter.
            fold_selector = FeatureSelector(
                random_state=fold_seed, n_jobs=1, verbose=False,
            )
            _run_feature_selection(
                fold_selector,
                methods=methods,
                X=X_tr, y=y_tr,
                gene_list=gene_list,
                de_genes=de_genes,
                species=species,
                n_select=min(n_select, X_tr.shape[1] // 2),
            )

            if not fold_selector.results:
                logger.warning(
                    f"   [pipeline-CV] fold {fold_i+1}: no feature "
                    f"selection methods succeeded; skipping fold."
                )
                # We still need to produce predictions for these samples
                # or we break the OOF invariant (every sample predicted
                # exactly once). Fall back to predicting the prior.
                prior = float(y_tr.mean())
                for clf_name in all_clfs:
                    oof_prob_by_clf[clf_name].extend([prior] * len(test_idx))
                oof_true_all.extend(y_te.tolist())
                per_fold_panels.append([])
                per_fold_ranked_genes.append(pd.DataFrame())
                continue

            fold_ranked = fold_selector.consensus_ranking(gene_list)
            fold_top = fold_ranked.head(max_panel).index.tolist()
            per_fold_ranked_genes.append(fold_ranked)

            # --- Stage 2 on training fold only ---
            # n_jobs=1 for the same reason as Stage 1: inner parallelism
            # on small folds loses to process-spawn overhead.
            fold_panel_optimizer = PanelOptimizer(
                random_state=fold_seed, n_jobs=1, verbose=False,
            )
            try:
                fold_panel_result = fold_panel_optimizer.forward_selection(
                    X_tr, y_tr,
                    ranked_genes=fold_top,
                    min_panel=min_panel,
                    max_panel=min(max_panel, len(fold_top)),
                    target_size=target_panel_size,
                    n_cv=min(5, len(train_idx) // 2),
                    auto_strategy=auto_panel_strategy,
                    sensitivity=panel_sensitivity,
                )
                fold_panel = fold_panel_result.optimal_panel
            except Exception as e:
                logger.warning(
                    f"   [pipeline-CV] fold {fold_i+1}: panel "
                    f"optimization failed ({e}); using top-{min_panel} "
                    f"from consensus ranking."
                )
                fold_panel = fold_top[:min_panel]

            per_fold_panels.append(fold_panel)

            # --- Stage 3: fit + predict per classifier ---
            if len(fold_panel) == 0:
                prior = float(y_tr.mean())
                for clf_name in all_clfs:
                    oof_prob_by_clf[clf_name].extend([prior] * len(test_idx))
                oof_true_all.extend(y_te.tolist())
                continue

            X_tr_panel = X_tr[fold_panel]
            X_te_panel = X_te[fold_panel]
            scaler = StandardScaler()
            X_tr_s = scaler.fit_transform(X_tr_panel)
            X_te_s = scaler.transform(X_te_panel)

            for clf_name, clf_template in all_clfs.items():
                clf = clone(clf_template)
                clf.fit(X_tr_s, y_tr)
                y_pred = clf.predict(X_te_s)
                y_prob = clf.predict_proba(X_te_s)[:, 1]

                oof_prob_by_clf[clf_name].extend(y_prob.tolist())

                # Per-fold metrics (matches pre-M6 ClassifierEvaluator
                # fold-metrics semantics for backward compat of
                # metrics_per_fold).
                acc = accuracy_score(y_te, y_pred)
                try:
                    auc_fold = roc_auc_score(y_te, y_prob)
                except ValueError:
                    auc_fold = 0.5
                f1_fold = f1_score(y_te, y_pred, zero_division=0)
                cm = confusion_matrix(y_te, y_pred, labels=[0, 1])
                tn, fp, fn, tp = cm.ravel()
                sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
                spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0

                per_fold_metrics[clf_name].append({
                    'repeat': repeat_idx,
                    'fold': fold_i,
                    'panel_size': len(fold_panel),
                    'accuracy': acc,
                    'auc': auc_fold,
                    'f1': f1_fold,
                    'sensitivity': sens,
                    'specificity': spec,
                })

            oof_true_all.extend(y_te.tolist())

    # --- Aggregate into ClassificationResult objects ---
    oof_true_arr = np.asarray(oof_true_all, dtype=int)
    clf_results: Dict[str, ClassificationResult] = {}

    for clf_name in all_clfs:
        oof_prob_arr = np.asarray(oof_prob_by_clf[clf_name], dtype=float)
        fold_rows = per_fold_metrics[clf_name]

        # Mean-across-folds scalar metrics (backward-compatible semantics
        # matching ClassifierEvaluator.evaluate_nested_cv's aggregate).
        if fold_rows:
            mean_acc = float(np.mean([m['accuracy'] for m in fold_rows]))
            mean_auc = float(np.mean([m['auc'] for m in fold_rows]))
            mean_f1 = float(np.mean([m['f1'] for m in fold_rows]))
            mean_sens = float(np.mean([m['sensitivity'] for m in fold_rows]))
            mean_spec = float(np.mean([m['specificity'] for m in fold_rows]))
        else:
            mean_acc = mean_auc = mean_f1 = mean_sens = mean_spec = 0.0

        # Pooled ROC from all OOF predictions
        roc_data = None
        try:
            fpr, tpr, _ = roc_curve(oof_true_arr, oof_prob_arr)
            roc_data = {'fpr': fpr.tolist(), 'tpr': tpr.tolist()}
        except Exception:
            pass

        clf_results[clf_name] = ClassificationResult(
            model_name=clf_name,
            accuracy=mean_acc,
            auc=mean_auc,
            sensitivity=mean_sens,
            specificity=mean_spec,
            f1=mean_f1,
            metrics_per_fold=fold_rows,
            roc_data=roc_data,
            feature_importance=None,   # populated later from final fit
            trained_model=None,        # populated later from final fit
            oof_true=oof_true_arr,
            oof_prob=oof_prob_arr,
        )

    return clf_results, per_fold_panels, per_fold_ranked_genes


@handle_errors(exit_on_error=False)
def discover_biomarkers(
    counts: Union[str, pd.DataFrame],
    metadata: Union[str, pd.DataFrame],
    group_column: str = 'condition',
    baseline_group: Optional[str] = None,
    reference_group: Optional[str] = None,
    de_result: Optional[Any] = None,
    ensemble_result: Optional[Any] = None,
    methods: Optional[List[str]] = None,
    target_panel_size: Optional[int] = None,
    min_panel: int = DEFAULT_PANEL_MIN,
    max_panel: int = DEFAULT_PANEL_MAX,
    validation: str = 'nested_cv',
    n_folds: int = DEFAULT_N_FOLDS_OUTER,
    n_repeats: int = 1,
    species: str = 'human',
    disease_term: Optional[str] = None,
    annotate: bool = True,
    run_literature: bool = True,
    run_ppi: bool = True,
    output_dir: Union[str, Path] = 'results/biomarkers',
    random_state: int = DEFAULT_RANDOM_STATE,
     verbose: bool = True,
    intent: Optional[Union[str, Any]] = None,
    prevalence: float = 0.05,
    calibrate_consensus: bool = True,
    alpha: float = 0.05,
    auto_panel_strategy: str = 'kneedle',
    panel_sensitivity: float = 1.0,
    panel_size_strategy: str = 'consensus',
) -> Union[BiomarkerResult, Any]:
    """
    Complete biomarker discovery pipeline.

    This is the main entry point for Module 10. It runs:
    1. Feature selection (multiple methods) — 10A
    2. Panel size optimization — 10C
    3. Classification evaluation — 10B
    4. Biological annotation & reporting — BiologicalAnnotator
    5. Assembly and output

    Parameters
    ----------
    counts : str or pd.DataFrame
        Count matrix (genes x samples). Path or DataFrame.
    metadata : str or pd.DataFrame
        Sample metadata with group_column. Path or DataFrame.
    group_column : str
        Column in metadata defining groups.
    de_result : DEResult, optional
        Module 7 DE result for initial filtering.
    ensemble_result : EnsembleResult, optional
        Module 9 ensemble result for consensus genes.
    methods : List[str], optional
        Feature selection methods to run. Default: auto-selected
        fold-safe methods (``univariate_filter``, ``elastic_net``,
        ``rfe``, plus optional ``boruta``/``mrmr``/``shap`` if available).
        Note: if ``de_result`` or ``ensemble_result`` is supplied,
        ``de_filter`` can be used with those upstream DE genes, but a
        ``UserWarning`` will flag the upstream leakage risk -- M10's
        pipeline-CV cannot correct for bias in a gene list that was
        itself computed from the same samples.
    target_panel_size : int, optional
        Target panel size. If None (default), auto-determined via
        ``auto_panel_strategy``.
    min_panel, max_panel : int
        Range for panel optimization.
    auto_panel_strategy : {'kneedle', 'argmax', 'first_drop'}
        Algorithm used to pick optimal_size when target_panel_size is
        not specified. Default 'kneedle' (Satopaa et al. 2011); falls
        back to argmax (smallest-size tiebreak) when no knee is
        detected. 'first_drop' is the legacy pre-kneedle behavior,
        retained for backward compatibility.
    panel_sensitivity : float
        Kneedle's S parameter (default 1.0); ignored for the other
        strategies.
    panel_size_strategy : {'per_fold', 'consensus'}
        Where panel-size detection runs in the nested_cv flow.
        'consensus' (default): run panel-size detection once at
        discovery level on the full data, then have all folds use
        that size as ``target_size`` (a 2-pass discovery). Per-fold
        size variability (especially when one fold's argmax tiebreak
        lands on a much larger size than the others) can drop Phi
        substantially even when the discovery-level K is stable;
        consensus mode prevents this. 'per_fold': each outer fold
        runs panel-size detection independently. Per-fold detection
        was an earlier candidate default but a Q5 sanity check on real
        fixtures (S2 seed 42: per-fold sizes (3,3,3,12,3) drove a
        0.313 Phi drop vs consensus) showed meaningful Phi
        instability, so consensus was made the default in v2.2.2.
    validation : str
        Validation strategy. 'nested_cv' (default) runs pipeline-CV
        with feature selection inside the outer fold, per Ambroise &
        McLachlan (2002 PNAS). 'loocv' is DEPRECATED as of v2.2.2:
        Ambroise & McLachlan explicitly recommend k-fold over LOOCV,
        because LOOCV has higher variance with no bias advantage once
        the external-CV leakage fix is applied; using it triggers a
        ``DeprecationWarning``.
    n_folds : int
        Number of outer CV folds. Default 5. For paper-grade results
        Ambroise & McLachlan (2002) recommend 10.
    n_repeats : int, default 1
        Number of outer-CV repeats. Per Bengio & Grandvalet (2004),
        repeated k-fold reduces the variance of the OOF AUC estimate;
        nestedcv (Lewis et al. 2023) uses 10x10-fold by default. Default
        1 to preserve runtime; set to 3 for dashboard-quality results or
        5 for paper-grade results. Runtime scales linearly with this
        parameter.
    species : str
        Species for annotation ('human', 'mouse', 'rat').
    disease_term : str, optional
        Disease context for literature search (e.g., 'breast cancer').
    annotate : bool
        Whether to run biological annotation pipeline.
    run_literature : bool
        Whether to query PubMed for gene-disease associations.
    run_ppi : bool
        Whether to query STRING for protein interactions.
    output_dir : str or Path
        Output directory.
    random_state : int
        Random seed.
    verbose : bool
        Print progress.

    Returns
    -------
    BiomarkerResult
        Complete biomarker discovery results. The new M6 field
        ``panel_stability`` (see ``PanelStabilityResult``) contains the
        Nogueira et al. (2018 JMLR) stability measure across per-fold
        panels and a benchmark label ('excellent', 'intermediate',
        'poor'). Honest OOF predictions are stored on each
        ``ClassificationResult.oof_true`` / ``.oof_prob`` (M1).

    Architecture (M6, post-2026-04)
    --------------------------------
    This function runs three logical stages:
        1. Pipeline-CV: in each outer fold, feature selection + panel
           optimization + classifier fit all happen on training data
           only; predictions for the held-out fold go into OOF arrays.
        2. Final-panel pass: Stages 1+2 re-run on ALL data to produce
           the single reported panel that a user would deploy.
        3. Annotation + reporting.
    OOF metrics (AUC, clinical metrics) are honest generalization
    estimates; the final panel is what you deploy. References:
    Ambroise & McLachlan 2002 PNAS 99:6562; Varma & Simon 2006 BMC
    Bioinf 7:91; Lewis et al. 2023 Bioinformatics Advances (nestedcv).

    Examples
    --------
    >>> from raptor.biomarker_discovery import discover_biomarkers
    >>> result = discover_biomarkers(
    ...     counts='counts.csv',
    ...     metadata='metadata.csv',
    ...     group_column='condition',
    ...     output_dir='results/biomarkers/',
    ... )
    >>> print(f"Panel: {result.panel}")
    >>> print(f"OOF AUC: {result.classification_results[result.best_classifier].auc:.3f}")
    >>> if result.panel_stability is not None:
    ...     print(result.panel_stability.summary())
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- M6: LOOCV deprecation warning ---
    # Ambroise & McLachlan (2002 PNAS) explicitly recommend k-fold over
    # LOOCV for gene-expression biomarker CV: LOOCV has higher variance
    # with no bias advantage once the external-CV leakage fix is
    # applied. We honor the user's request but warn once.
    if validation == 'loocv':
        warnings.warn(
            "validation='loocv' is DEPRECATED as of RAPTOR v2.2.2. "
            "Ambroise & McLachlan (2002 PNAS 99:6562) explicitly "
            "recommend k-fold cross-validation over leave-one-out: "
            "LOOCV has high variance and offers no bias advantage once "
            "the external-CV leakage fix is in place. Use "
            "validation='nested_cv' with n_folds=5 or n_folds=10 "
            "(Ambroise & McLachlan's recommendation) instead. "
            "LOOCV is still supported for now but may be removed in a "
            "future release.",
            DeprecationWarning,
            stacklevel=2,
        )

    logger.info("=" * 70)
    logger.info("RAPTOR v2.2.2 - MODULE 10: BIOMARKER DISCOVERY")
    logger.info("=" * 70)
    logger.info("")

    # --- Load data ---
    if isinstance(counts, (str, Path)):
        logger.info(f"Loading counts from: {counts}")
        counts = pd.read_csv(counts, index_col=0)
    if isinstance(metadata, (str, Path)):
        logger.info(f"Loading metadata from: {metadata}")
        metadata = pd.read_csv(metadata)

    # --- Get DE gene list if available ---
    de_genes = None
    if ensemble_result is not None and hasattr(ensemble_result, 'consensus_genes'):
        cg = ensemble_result.consensus_genes
        if 'gene_id' in cg.columns:
            de_genes = cg['gene_id'].tolist()
        else:
            de_genes = cg.index.tolist()
        logger.info(f"Using {len(de_genes)} consensus genes from M9 ensemble")
    elif de_result is not None and hasattr(de_result, 'significant_genes'):
        de_genes = de_result.significant_genes
        logger.info(f"Using {len(de_genes)} significant genes from M7 DE result")

    # --- Prepare data ---
    logger.info("Preparing expression data...")
    X, y, gene_list = _prepare_expression_data(
        counts, metadata, group_column, de_genes=de_genes,
        baseline_group=baseline_group,
        reference_group=reference_group,
    )
    n_samples = X.shape[0]
    n_genes = X.shape[1]

    # --- Determine available methods ---
    if methods is None:
        # Default method set: fold-safe choices first, with optional
        # heavyweights when their libraries are installed. Note that
        # 'de_filter' is NOT in the default list -- use 'univariate_filter'
        # instead (fold-safe per-fold Welch's t-test, matches nestedcv
        # convention per Lewis et al. 2023).
        methods = ['univariate_filter', 'elastic_net', 'rfe']
        if de_genes:
            # User supplied precomputed DE genes -- honor the request but
            # _resolve_pipeline_cv_methods below will emit the leakage
            # warning.
            methods.insert(0, 'de_filter')
        if _BORUTA_AVAILABLE:
            methods.append('boruta')
        if _MRMR_AVAILABLE:
            methods.append('mrmr')
        if _SHAP_AVAILABLE:
            methods.append('shap')
        if _PYWGCNA_AVAILABLE and n_samples >= 15:
            methods.append('wgcna')

    # Resolve methods for pipeline-CV: substitute univariate_filter for
    # de_filter when no independent de_genes is supplied; warn when
    # precomputed de_genes is supplied (upstream leakage risk).
    fold_methods = _resolve_pipeline_cv_methods(methods, de_genes)

    n_select = min(50, n_genes // 2)  # Target features per method

    # =====================================================================
    # STAGE 1+2+3 (M6 PIPELINE-CV): feature selection + panel optimization
    # + classifier evaluation, with every y-using step confined to the
    # training fold. See _run_pipeline_cv docstring and Ambroise &
    # McLachlan 2002 PNAS for the theoretical justification.
    # =====================================================================
    logger.info("")
    logger.info(
        f"STAGE 1-3 (pipeline-CV): feature selection + panel + "
        f"classifier evaluation inside {n_folds}-fold outer CV "
        f"x {n_repeats} repeat(s)"
    )
    logger.info("-" * 70)
    logger.info(
        "   Ambroise & McLachlan (2002 PNAS 99:6562): feature selection "
        "must be external to CV for unbiased error estimates."
    )

    # For LOOCV users (deprecated), route to the classical path unchanged
    # -- pipeline-LOOCV is explicitly out of scope per scoping doc Q1.
    if validation == 'loocv':
        logger.info("   Using classical LOOCV path (deprecated, see warning above).")
        # Classical Stage 1+2 on all data, Stage 3 LOOCV.
        selector = FeatureSelector(
            random_state=random_state, verbose=verbose,
        )
        _run_feature_selection(
            selector, methods=fold_methods, X=X, y=y,
            gene_list=gene_list, de_genes=de_genes, species=species,
            n_select=n_select,
        )
        logger.info("")
        logger.info("Computing consensus ranking...")
        ranked_genes = selector.consensus_ranking(gene_list)

        # --- M4: Significance-calibrated consensus ranking (LOOCV path) ---
        if calibrate_consensus:
            try:
                ranked_genes = apply_significance_calibration(
                    ranked_genes=ranked_genes,
                    X=X, y=y,
                    reference_group=reference_group,
                    baseline_group=baseline_group,
                    alpha=alpha,
                )
            except Exception as e:
                logger.warning(
                    f"   M4 calibration failed ({e}); ranked_genes "
                    f"retains uncalibrated consensus_score."
                )

        top_candidates = ranked_genes.head(max_panel).index.tolist()

        panel_optimizer = PanelOptimizer(
            random_state=random_state, verbose=verbose,
        )
        panel_result = panel_optimizer.forward_selection(
            X, y,
            ranked_genes=top_candidates,
            min_panel=min_panel,
            max_panel=min(max_panel, len(top_candidates)),
            target_size=target_panel_size,
            n_cv=min(n_folds, n_samples // 2),
            auto_strategy=auto_panel_strategy,
            sensitivity=panel_sensitivity,
        )
        panel_genes = panel_result.optimal_panel

        evaluator = ClassifierEvaluator(
            random_state=random_state, verbose=verbose,
        )
        X_panel = X[panel_genes]
        clf_results = evaluator.evaluate_loocv(X_panel, y)
        panel_stability_result = None  # no stability under LOOCV path

    else:
        # --- Consensus scout (panel_size_strategy='consensus', default) ---
        # Each outer fold uses the SAME panel size, determined by a
        # quick discovery-level pass on the full data BEFORE the
        # pipeline-CV. Per-fold detection is the alternative
        # (panel_size_strategy='per_fold') and skips this block, but
        # it can drop Phi substantially when a single fold's
        # auto-detector lands on a much larger size than the others
        # (e.g., on saturated argmax_fallback curves where ties at
        # different sizes resolve differently per fold). Q5 sanity
        # check on S2 (seed 42) showed a 0.313 Phi drop under
        # per_fold vs consensus driven by exactly this pattern.
        effective_target_panel_size = target_panel_size
        # Track whether effective_target_panel_size came from the
        # consensus scout (for selection_method labeling on the
        # final-panel pass below).
        target_from_consensus_scout = False
        if (
            panel_size_strategy == 'consensus'
            and target_panel_size is None
        ):
            logger.info("")
            logger.info(
                "Consensus scout: detecting panel size on full data "
                "before pipeline-CV"
            )
            logger.info("-" * 50)
            scout_selector = FeatureSelector(
                random_state=random_state, verbose=False,
            )
            _run_feature_selection(
                scout_selector, methods=fold_methods, X=X, y=y,
                gene_list=gene_list, de_genes=de_genes, species=species,
                n_select=n_select,
            )
            scout_ranked = scout_selector.consensus_ranking(gene_list)
            scout_top = scout_ranked.head(max_panel).index.tolist()
            scout_optimizer = PanelOptimizer(
                random_state=random_state, n_jobs=-1, verbose=False,
            )
            scout_result = scout_optimizer.forward_selection(
                X, y,
                ranked_genes=scout_top,
                min_panel=min_panel,
                max_panel=min(max_panel, len(scout_top)),
                target_size=None,
                n_cv=min(n_folds, len(y) // 2),
                auto_strategy=auto_panel_strategy,
                sensitivity=panel_sensitivity,
            )
            effective_target_panel_size = scout_result.optimal_size
            target_from_consensus_scout = True
            logger.info(
                f"   Consensus K = {effective_target_panel_size} "
                f"(method={scout_result.selection_method}); "
                f"all folds will use this size."
            )

        # --- Pipeline-CV: Stages 1+2+3 inside the outer fold ---
        clf_results, per_fold_panels, per_fold_ranked = _run_pipeline_cv(
            X=X, y=y,
            gene_list=gene_list,
            methods=fold_methods,
            de_genes=de_genes,
            species=species,
            n_outer=n_folds,
            n_repeats=n_repeats,
            n_select=n_select,
            min_panel=min_panel,
            max_panel=max_panel,
            target_panel_size=effective_target_panel_size,
            random_state=random_state,
            verbose=verbose,
            auto_panel_strategy=auto_panel_strategy,
            panel_sensitivity=panel_sensitivity,
        )

        # --- Final-panel pass: Stages 1+2 on ALL data ---
        # This is the panel users see as "the panel" -- the most data
        # you can use to pick a deployable model. The OOF AUC from the
        # loop above estimates generalization of THE PROCEDURE; the
        # final panel is what you would actually deploy. This matches
        # the nestedcv R package convention (Lewis et al. 2023).
        logger.info("")
        logger.info("Final-panel pass: Stages 1+2 on all data")
        logger.info("-" * 40)

        selector = FeatureSelector(
            random_state=random_state, verbose=verbose,
        )
        _run_feature_selection(
            selector, methods=fold_methods, X=X, y=y,
            gene_list=gene_list, de_genes=de_genes, species=species,
            n_select=n_select,
        )
        logger.info("")
        logger.info("Computing consensus ranking on all data...")
        ranked_genes = selector.consensus_ranking(gene_list)

        # --- M4: Significance-calibrated consensus ranking ---
        # Runs on the final all-data ranking (not inside the pipeline-CV
        # loop, per scoping doc §4). Multiplies consensus_score by a
        # two-tier weight based on per-gene DE p-value; re-sorts by the
        # calibrated score. Forward selection below then draws its top
        # candidates from the calibrated ordering, so noise genes that
        # happen to score high on the raw formula but fail p < alpha
        # get demoted.
        if calibrate_consensus:
            try:
                ranked_genes = apply_significance_calibration(
                    ranked_genes=ranked_genes,
                    X=X, y=y,
                    reference_group=reference_group,
                    baseline_group=baseline_group,
                    alpha=alpha,
                )
            except Exception as e:
                logger.warning(
                    f"   M4 calibration failed ({e}); ranked_genes "
                    f"retains uncalibrated consensus_score. Forward "
                    f"selection will use the raw ranking."
                )

        top_candidates = ranked_genes.head(max_panel).index.tolist()
        logger.info(f"   Top candidate genes: {len(top_candidates)}")

        panel_optimizer = PanelOptimizer(
            random_state=random_state, verbose=verbose,
        )
        # Use the same effective_target_panel_size resolved above:
        # under panel_size_strategy='consensus' this is the scout's K
        # so the deployable panel matches what the per-fold estimator
        # was tuned for. Under 'per_fold' (default) it's None, and
        # the final pass auto-detects on the full data.
        panel_result = panel_optimizer.forward_selection(
            X, y,
            ranked_genes=top_candidates,
            min_panel=min_panel,
            max_panel=min(max_panel, len(top_candidates)),
            target_size=effective_target_panel_size,
            n_cv=n_folds,
            auto_strategy=auto_panel_strategy,
            sensitivity=panel_sensitivity,
        )
        # When effective_target_panel_size came from the consensus
        # scout (not from the user), forward_selection labels it
        # 'user_specified' because it can't tell the difference at
        # its level. Override here so the dashboard annotation shows
        # 'Target = K (consensus-pinned)' rather than misattributing
        # the choice to the user.
        if target_from_consensus_scout:
            panel_result.selection_method = 'consensus_pinned'
        panel_genes = panel_result.optimal_panel

        # --- Populate trained_model + feature_importance per classifier ---
        # The OOF predictions on clf_results came from the pipeline-CV
        # loop above. For downstream code that uses `trained_model`
        # (enhanced.py clinical metrics, ratio biomarker searches,
        # applying the panel to new data), we now fit each classifier
        # on ALL data with the final panel. This is the M1-introduced
        # pattern -- preserved here so M1 + M6 compose correctly.
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline
        from sklearn.base import clone

        final_evaluator = ClassifierEvaluator(
            random_state=random_state, verbose=False,
        )
        all_clfs_final = final_evaluator._get_classifiers()
        X_final_panel = X[panel_genes] if panel_genes else X

        for clf_name, clf_template in all_clfs_final.items():
            if clf_name not in clf_results:
                # Skip classifiers that weren't run in the CV loop
                continue
            clf_result = clf_results[clf_name]
            try:
                scaler_final = StandardScaler()
                clf_final = clone(clf_template)
                full_pipeline = Pipeline([
                    ("scaler", scaler_final),
                    ("clf", clf_final),
                ])
                full_pipeline.fit(X_final_panel, y)

                feat_imp = None
                if hasattr(clf_final, 'feature_importances_'):
                    feat_imp = pd.DataFrame(
                        {'importance': clf_final.feature_importances_},
                        index=X_final_panel.columns,
                    )
                elif hasattr(clf_final, 'coef_'):
                    feat_imp = pd.DataFrame(
                        {'importance': np.abs(clf_final.coef_).ravel()},
                        index=X_final_panel.columns,
                    )
                if feat_imp is not None:
                    feat_imp.index.name = 'gene_id'
                    feat_imp = feat_imp.sort_values(
                        'importance', ascending=False
                    )

                clf_result.trained_model = full_pipeline
                clf_result.feature_importance = feat_imp
            except Exception as e:
                logger.warning(
                    f"   Final-panel fit failed for {clf_name}: {e}. "
                    f"trained_model will be None; OOF predictions still valid."
                )

        # --- Panel stability diagnostic (Nogueira 2018) ---
        panel_stability_result = _compute_panel_stability(
            per_fold_panels=per_fold_panels,
            per_fold_ranked_genes=per_fold_ranked,
            gene_universe=gene_list,
            final_panel=panel_genes,
            n_folds=n_folds,
            n_repeats=n_repeats,
            random_state=random_state,
        )
        logger.info("")
        logger.info(panel_stability_result.summary())

    # Best classifier (deterministic tiebreak: AUC -> F1 -> interpretability preference).
    # See _select_best_classifier for the full rule set.
    best_clf = _select_best_classifier(clf_results)

    # --- Biological Annotation ---
    annotation_res = None

    if annotate:
        logger.info("")
        logger.info("STAGE 4: Biological Annotation & Reporting")
        logger.info("-" * 40)

        try:
            annotator = BiologicalAnnotator(species=species, verbose=verbose)
            annotation_res = annotator.annotate_panel(
                panel_genes=panel_genes,
                background_genes=gene_list,
                disease_term=disease_term,
                run_literature=run_literature,
                run_ppi=run_ppi,
            )

            # Generate report
            report_path = output_dir / "biomarker_report.md"
            annotator.generate_report(
                panel_genes=panel_genes,
                annotation_result=annotation_res,
                classification_results=clf_results,
                panel_optimization=panel_result,
                output_path=report_path,
            )

        except Exception as e:
            logger.warning(f"   Annotation stage encountered an error: {e}")
            logger.warning("   Continuing without annotations.")

    # --- Assemble result ---
    logger.info("")
    logger.info("STAGE 5: Assembling results")
    logger.info("-" * 40)

    result = BiomarkerResult(
        ranked_genes=ranked_genes,
        panel=panel_genes,
        panel_size=len(panel_genes),
        selection_results=selector.results,
        classification_results=clf_results,
        best_classifier=best_clf,
        panel_optimization=panel_result,
        annotation_result=annotation_res,
        panel_stability=panel_stability_result,
        study_design='binary',
        validation_strategy=validation,
        n_samples=n_samples,
        n_initial_candidates=n_genes,
        parameters={
            'methods': methods,
            'target_panel_size': target_panel_size,
            'min_panel': min_panel,
            'max_panel': max_panel,
            'n_folds': n_folds,
            'n_repeats': n_repeats,
            'random_state': random_state,
            'group_column': group_column,
            'species': species,
            'disease_term': disease_term,
            'annotate': annotate,
            # M4 parameters
            'calibrate_consensus': calibrate_consensus,
            'alpha': alpha,
        },
        metadata={
            'raptor_version': __raptor_version__,
            'module': 'M10',
        },
    )

    # --- Save ---
    result.save(output_dir)

    logger.info("")
    logger.info(result.summary())
    logger.info("")
    logger.info(f"Results saved to: {output_dir}/")
    logger.info("")
    logger.info("Output files:")
    logger.info(f"   biomarker_panel.csv     - {panel_result.optimal_size}-gene panel")
    logger.info(f"   ranked_genes.csv        - Full consensus ranking")
    logger.info(f"   classification_performance.csv")
    logger.info(f"   panel_curve.csv         - Panel size vs AUC")
    if annotation_res is not None:
        logger.info(f"   biomarker_report.md     - Publication-ready report")
        logger.info(f"   annotations/            - Gene info, pathways, literature, PPI")
    logger.info(f"   biomarker_result.pkl    - Complete result (for downstream)")
    logger.info("")

    # --- Enhanced analysis (if intent is set) ---
    if intent is not None:
        try:
            from raptor.biomarker_discovery.enhanced import enhance_biomarker_result
            # Use caller-specified groups when available so baseline/reference
            # interpretation matches the y encoding from _prepare_expression_data.
            # Previously this always used sorted(unique_groups), which silently
            # disagreed with the caller's baseline/reference choice whenever
            # the reference group sorted before the baseline alphabetically
            # (e.g. "disease" < "healthy"), flipping every direction/ratio/
            # signature-coefficient interpretation downstream.
            if baseline_group is not None and reference_group is not None:
                group_names = (baseline_group, reference_group)
            else:
                unique_groups = sorted(metadata[group_column].dropna().unique())
                group_names = (unique_groups[0], unique_groups[1])
            result = enhance_biomarker_result(
                base_result=result,
                expression=X,
                labels=y,
                group_names=group_names,
                intent=intent,
                prevalence=prevalence,
                random_state=random_state,
                verbose=verbose,
            )
        except Exception as e:
            logger.warning(f"Enhanced analysis failed: {e}")
            logger.warning("Returning base BiomarkerResult.")

    return result


def discover_survival_biomarkers(
    counts: Union[str, pd.DataFrame],
    clinical: Union[str, pd.DataFrame],
    time_column: str = 'os_time',
    event_column: str = 'os_event',
    de_genes: Optional[List[str]] = None,
    fdr_threshold: float = 0.05,
    alpha: float = 0.1,
    output_dir: Union[str, Path] = 'results/survival_biomarkers',
    verbose: bool = True,
) -> SurvivalResult:
    """
    Survival-based biomarker discovery using Cox regression.

    Parameters
    ----------
    counts : str or pd.DataFrame
        Expression matrix (genes x samples).
    clinical : str or pd.DataFrame
        Clinical data with time and event columns.
    time_column : str
        Column name for survival time.
    event_column : str
        Column name for event indicator (1=event, 0=censored).
    de_genes : List[str], optional
        Restrict analysis to these candidate genes.
    fdr_threshold : float
        FDR threshold for univariate Cox screen.
    alpha : float
        Regularization for CoxNet panel selection.
    output_dir : str or Path
        Output directory.
    verbose : bool
        Print progress.

    Returns
    -------
    SurvivalResult
    """
    if not _LIFELINES_AVAILABLE:
        raise DependencyError(
            "lifelines is required for survival analysis: pip install lifelines"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("RAPTOR v2.2.2 - MODULE 10: SURVIVAL BIOMARKER DISCOVERY")
    logger.info("=" * 70)

    # Load data
    if isinstance(counts, (str, Path)):
        counts = pd.read_csv(counts, index_col=0)
    if isinstance(clinical, (str, Path)):
        clinical = pd.read_csv(clinical)

    # Normalize
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # Align samples
    if 'sample_id' in clinical.columns:
        clinical = clinical.set_index('sample_id')

    common = list(set(log_cpm.columns.astype(str)) & set(clinical.index.astype(str)))
    if len(common) < 10:
        raise ValidationError(
            'samples',
            f"Only {len(common)} matching samples",
            "Need at least 10 samples for survival analysis"
        )

    X = log_cpm[common].T
    time_arr = clinical.loc[common, time_column].values.astype(float)
    event_arr = clinical.loc[common, event_column].values.astype(int)

    # Restrict to DE genes if provided
    if de_genes:
        valid = [g for g in de_genes if g in X.columns]
        if len(valid) >= 10:
            X = X[valid]
            logger.info(f"Restricted to {len(valid)} DE candidate genes")

    # Univariate screen
    analyzer = SurvivalAnalyzer(verbose=verbose)
    cox_screen = analyzer.cox_univariate_screen(
        X, time_arr, event_arr, fdr_threshold=fdr_threshold,
    )

    sig_genes = cox_screen[cox_screen['is_significant']]['gene_id'].tolist()

    # CoxNet panel
    survival_result = SurvivalResult(significant_genes=sig_genes)

    if len(sig_genes) >= 3:
        X_sig = X[sig_genes]
        survival_result = analyzer.cox_lasso_panel(
            X_sig, time_arr, event_arr, alpha=alpha,
        )
        survival_result.significant_genes = sig_genes
        survival_result.cox_results = cox_screen.set_index('gene_id')

    # Save
    cox_screen.to_csv(output_dir / "cox_univariate_screen.csv", index=False)

    if survival_result.cox_results is not None:
        panel_df = pd.DataFrame({'gene_id': survival_result.panel_genes})
        panel_df.to_csv(output_dir / "survival_panel.csv", index=False)

    with open(output_dir / "survival_summary.json", 'w', encoding='utf-8') as f:
        json.dump({
            'n_genes_screened': len(X.columns),
            'n_significant': len(sig_genes),
            'n_panel': len(survival_result.panel_genes),
            'c_index': survival_result.c_index,
            'fdr_threshold': fdr_threshold,
            'alpha': alpha,
        }, f, indent=2)

    logger.info(f"\nResults saved to: {output_dir}/")
    return survival_result


def validate_biomarkers(
    panel_genes: List[str],
    counts: Union[str, pd.DataFrame],
    metadata: Union[str, pd.DataFrame],
    group_column: str = 'condition',
    baseline_group: Optional[str] = None,
    reference_group: Optional[str] = None,
    n_folds: int = 5,
    random_state: int = DEFAULT_RANDOM_STATE,
    verbose: bool = True,
) -> Dict[str, ClassificationResult]:
    """
    Validate a biomarker panel on an independent dataset.

    Parameters
    ----------
    panel_genes : List[str]
        Gene panel to validate.
    counts : str or pd.DataFrame
        Validation cohort expression matrix.
    metadata : str or pd.DataFrame
        Validation cohort metadata.
    group_column : str
        Group column in metadata.
    baseline_group, reference_group : str, optional
        Explicit baseline (label 0) and reference (label 1) group names.
        If both omitted, falls back to alphabetical sort. Should match the
        polarity used during discovery so panel directionality is
        interpreted consistently across cohorts.
    n_folds : int
        Number of CV folds.

    Returns
    -------
    Dict[str, ClassificationResult]
        Classification performance on validation data.
    """
    if isinstance(counts, (str, Path)):
        counts = pd.read_csv(counts, index_col=0)
    if isinstance(metadata, (str, Path)):
        metadata = pd.read_csv(metadata)

    X, y, _ = _prepare_expression_data(
        counts, metadata, group_column, de_genes=panel_genes,
        baseline_group=baseline_group,
        reference_group=reference_group,
    )

    # Restrict to panel genes present in validation data
    available = [g for g in panel_genes if g in X.columns]
    if len(available) < len(panel_genes):
        logger.warning(
            f"Only {len(available)}/{len(panel_genes)} panel genes "
            f"found in validation data"
        )

    X_panel = X[available]

    evaluator = ClassifierEvaluator(
        random_state=random_state, verbose=verbose,
    )

    return evaluator.evaluate_nested_cv(X_panel, y, n_outer=n_folds)


# =============================================================================
# MODULE INFORMATION
# =============================================================================

def get_dependencies_status() -> Dict[str, bool]:
    """Check availability of all optional dependencies."""
    return {
        'scikit-learn': _SKLEARN_AVAILABLE,
        'xgboost': _XGBOOST_AVAILABLE,
        'shap': _SHAP_AVAILABLE,
        'boruta': _BORUTA_AVAILABLE,
        'mrmr': _MRMR_AVAILABLE,
        'lifelines': _LIFELINES_AVAILABLE,
        'PyWGCNA': _PYWGCNA_AVAILABLE,
    }


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

def main():
    """CLI entry point for standalone execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description='RAPTOR Module 10 - Biomarker Discovery',
    )

    parser.add_argument('--counts', '-c', required=True,
                        help='Count matrix CSV (genes x samples)')
    parser.add_argument('--metadata', '-m', required=True,
                        help='Sample metadata CSV')
    parser.add_argument('--group-column', '-g', default='condition',
                        help='Group column in metadata')
    parser.add_argument('--de-result', '-d',
                        help='DE result pickle from M7')
    parser.add_argument('--ensemble-result', '-e',
                        help='Ensemble result pickle from M9')
    parser.add_argument('--methods', nargs='+',
                        help='Feature selection methods')
    parser.add_argument('--panel-size', type=int,
                        help='Target panel size')
    parser.add_argument('--species', default='human',
                        help='Species for annotation (human, mouse, rat)')
    parser.add_argument('--disease-term',
                        help='Disease context for literature search')
    parser.add_argument('--no-annotate', action='store_true',
                        help='Skip biological annotation')
    parser.add_argument('--no-literature', action='store_true',
                        help='Skip literature mining')
    parser.add_argument('--no-ppi', action='store_true',
                        help='Skip STRING PPI query')
    parser.add_argument('--output', '-o', default='results/biomarkers',
                        help='Output directory')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--intent', choices=['diagnostic', 'exploratory'],
                        default=None, help='Biomarker intent for enhanced analyses')
    parser.add_argument('--prevalence', type=float, default=0.05,
                        help='Disease prevalence for PPV/NPV (default: 0.05)')

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(message)s',
    )

    # Load optional inputs
    de_res = None
    if args.de_result:
        with open(args.de_result, 'rb') as f:
            de_res = pickle.load(f)

    ens_res = None
    if args.ensemble_result:
        with open(args.ensemble_result, 'rb') as f:
            ens_res = pickle.load(f)

    result = discover_biomarkers(
        counts=args.counts,
        metadata=args.metadata,
        group_column=args.group_column,
        de_result=de_res,
        ensemble_result=ens_res,
        methods=args.methods,
        target_panel_size=args.panel_size,
        species=args.species,
        disease_term=args.disease_term,
        annotate=not args.no_annotate,
        run_literature=not args.no_literature,
        run_ppi=not args.no_ppi,
        output_dir=args.output,
        verbose=args.verbose,
        intent=args.intent,
        prevalence=args.prevalence,
    )

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())