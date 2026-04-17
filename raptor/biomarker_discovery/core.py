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

        lines.append("\nCLASSIFICATION PERFORMANCE:")
        for name, res in self.classification_results.items():
            marker = " *" if name == self.best_classifier else ""
            lines.append(f"  {name}: AUC={res.auc:.3f}, F1={res.f1:.3f}{marker}")

        lines.append(f"\nRECOMMENDED PANEL: {self.panel_size} genes")
        if self.panel_optimization:
            lines.append(
                f"  Panel AUC: {self.panel_optimization.optimal_auc:.3f}"
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

        # Ranking: 1 = selected
        max_rank = rfe.ranking_.max()
        importance = (max_rank - rfe.ranking_ + 1) / max_rank

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

            # Feature importance from full-data fit
            scaler_full = StandardScaler()
            X_full_s = scaler_full.fit_transform(X)
            clf_full = clone(clf_template)
            clf_full.fit(X_full_s, y)

            feat_imp = None
            if hasattr(clf_full, 'feature_importances_'):
                feat_imp = pd.DataFrame({
                    'importance': clf_full.feature_importances_
                }, index=X.columns)
            elif hasattr(clf_full, 'coef_'):
                feat_imp = pd.DataFrame({
                    'importance': np.abs(clf_full.coef_).ravel()
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
                trained_model=clf_full,
            )

        return results

    def evaluate_loocv(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        classifiers: Optional[List[str]] = None,
    ) -> Dict[str, ClassificationResult]:
        """Leave-one-out cross-validation for small sample sizes."""
        from sklearn.model_selection import LeaveOneOut
        from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
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

            self._log(f"   {clf_name}: AUC={auc:.3f}, Acc={acc:.3f}")

            results[clf_name] = ClassificationResult(
                model_name=clf_name,
                accuracy=acc,
                auc=auc,
                f1=f1,
            )

        return results


# =============================================================================
# 10C: PANEL OPTIMIZATION
# =============================================================================

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
            If specified, the panel will be this exact size.
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
        else:
            # Auto-detect optimal via elbow: find where adding genes
            # gives < 0.5% AUC improvement
            optimal_idx = 0
            for i in range(1, len(curve_data)):
                improvement = curve_data[i]['auc_mean'] - curve_data[i-1]['auc_mean']
                if improvement < 0.005:
                    break
                optimal_idx = i
            optimal_size = curve_data[optimal_idx]['panel_size']
            optimal_auc = curve_data[optimal_idx]['auc_mean']

        optimal_panel = all_panels[optimal_size]

        self._log(
            f"   Optimal panel: {optimal_size} genes, AUC={optimal_auc:.3f}"
        )

        return PanelOptimizationResult(
            optimal_panel=optimal_panel,
            optimal_size=optimal_size,
            optimal_auc=optimal_auc,
            panel_curve=curve_df,
            all_panels=all_panels,
            method='forward_selection',
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
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """
    Prepare expression matrix and labels for ML.

    - Filters low-count genes
    - Log2-CPM normalizes
    - Encodes labels as 0/1
    - Optionally restricts to DE genes

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

    label_map = {unique_groups[0]: 0, unique_groups[1]: 1}
    y = np.array([label_map[g] for g in groups])

    logger.info(
        f"   Prepared: {X.shape[0]} samples x {X.shape[1]} genes, "
        f"groups: {unique_groups[0]}={int((y==0).sum())}, "
        f"{unique_groups[1]}={int((y==1).sum())}"
    )

    return X, y, list(X.columns)


@handle_errors(exit_on_error=False)
def discover_biomarkers(
    counts: Union[str, pd.DataFrame],
    metadata: Union[str, pd.DataFrame],
    group_column: str = 'condition',
    de_result: Optional[Any] = None,
    ensemble_result: Optional[Any] = None,
    methods: Optional[List[str]] = None,
    target_panel_size: Optional[int] = None,
    min_panel: int = DEFAULT_PANEL_MIN,
    max_panel: int = DEFAULT_PANEL_MAX,
    validation: str = 'nested_cv',
    n_folds: int = DEFAULT_N_FOLDS_OUTER,
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
        Feature selection methods to run. Default: all available.
    target_panel_size : int, optional
        Target panel size. If None, auto-determined.
    min_panel, max_panel : int
        Range for panel optimization.
    validation : str
        Validation strategy: 'nested_cv', 'loocv'.
    n_folds : int
        Number of CV folds.
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
        Complete biomarker discovery results.

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
    >>> print(f"AUC: {result.panel_optimization.optimal_auc:.3f}")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

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
    )
    n_samples = X.shape[0]
    n_genes = X.shape[1]

    # --- Determine available methods ---
    if methods is None:
        methods = ['elastic_net', 'rfe']
        if de_genes:
            methods.insert(0, 'de_filter')
        if _BORUTA_AVAILABLE:
            methods.append('boruta')
        if _MRMR_AVAILABLE:
            methods.append('mrmr')
        if _SHAP_AVAILABLE:
            methods.append('shap')
        if _PYWGCNA_AVAILABLE and n_samples >= 15:
            methods.append('wgcna')

    # --- 10A: Feature Selection ---
    logger.info("")
    logger.info("STAGE 1: Feature Selection")
    logger.info("-" * 40)

    selector = FeatureSelector(
        random_state=random_state, verbose=verbose,
    )

    n_select = min(50, n_genes // 2)  # Target features per method

    # Species mapping for WGCNA
    _wgcna_species_map = {
        'human': 'homo sapiens', 'mouse': 'mus musculus',
        'rat': 'rattus norvegicus',
    }

    for method in methods:
        try:
            if method == 'de_filter' and de_genes:
                selector.select_de_filter(de_genes, gene_list)
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
                wgcna_sp = _wgcna_species_map.get(species, species)
                selector.select_wgcna(X, y, species=wgcna_sp)
        except DependencyError as e:
            logger.warning(f"   Skipping {method}: {e}")
        except Exception as e:
            logger.warning(f"   Error in {method}: {e}")

    # Consensus ranking
    logger.info("")
    logger.info("Computing consensus ranking...")
    ranked_genes = selector.consensus_ranking(gene_list)
    top_candidates = ranked_genes.head(max_panel).index.tolist()
    logger.info(f"   Top candidate genes: {len(top_candidates)}")

    # --- 10C: Panel Optimization ---
    logger.info("")
    logger.info("STAGE 2: Panel Optimization")
    logger.info("-" * 40)

    panel_optimizer = PanelOptimizer(
        random_state=random_state, verbose=verbose,
    )

    panel_result = panel_optimizer.forward_selection(
        X, y,
        ranked_genes=top_candidates,
        min_panel=min_panel,
        max_panel=min(max_panel, len(top_candidates)),
        target_size=target_panel_size,
        n_cv=n_folds,
    )

    panel_genes = panel_result.optimal_panel

    # --- 10B: Classification Evaluation ---
    logger.info("")
    logger.info("STAGE 3: Classification Evaluation")
    logger.info("-" * 40)

    evaluator = ClassifierEvaluator(
        random_state=random_state, verbose=verbose,
    )

    X_panel = X[panel_genes]

    if validation == 'loocv' or n_samples < 20:
        logger.info("   Using Leave-One-Out CV (small sample size)")
        clf_results = evaluator.evaluate_loocv(X_panel, y)
    else:
        clf_results = evaluator.evaluate_nested_cv(
            X_panel, y, n_outer=n_folds,
        )

    # Best classifier
    best_clf = max(clf_results.keys(), key=lambda k: clf_results[k].auc)

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
            'random_state': random_state,
            'group_column': group_column,
            'species': species,
            'disease_term': disease_term,
            'annotate': annotate,
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
