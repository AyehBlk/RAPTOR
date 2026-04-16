"""
Signature Score Engine
======================

Converts a biomarker gene panel into a weighted risk score per patient.

A signature score is a single number that summarizes a patient's risk or
disease state using a panel of biomarker genes. For each patient, expression
values are z-score normalized, multiplied by learned weights, and summed.
The result is a scalar score that can be compared against learned cutoffs
to stratify patients into risk groups (low / medium / high).

This module provides:
    - build_signature_score(): factory that fits a score from training data
    - SignatureScore: the fitted-score object with .score() and .stratify()

Supported modes
---------------
diagnostic : weights from logistic regression, cutoff from Youden's J
prognostic : weights from Cox regression, cutoff from log-rank split
             (requires optional `lifelines` dependency)

Usage
-----
>>> from raptor.biomarker_discovery import build_signature_score
>>> # X_train: samples x genes DataFrame; y_train: 0/1 labels
>>> sig = build_signature_score(
...     X_train, y_train, panel_genes=["TAU1", "APP", "PSEN1"],
...     mode="diagnostic",
... )
>>> # Apply to a new cohort
>>> scores = sig.score(X_new)
>>> groups = sig.stratify(scores)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Optional dependency check
# ---------------------------------------------------------------------------

try:
    from lifelines import CoxPHFitter  # noqa: F401
    _LIFELINES_AVAILABLE = True
except ImportError:
    _LIFELINES_AVAILABLE = False


# ---------------------------------------------------------------------------
# Supported modes
# ---------------------------------------------------------------------------

VALID_MODES = ("diagnostic", "prognostic")


# ---------------------------------------------------------------------------
# Main dataclass
# ---------------------------------------------------------------------------

@dataclass
class SignatureScore:
    """
    A fitted signature score. Apply to new expression data to get per-patient
    risk scores, or stratify into risk groups using learned cutoffs.

    Fields
    ------
    mode : str
        Which mode this score was built for: 'diagnostic' or 'prognostic'.
    weights : Dict[str, float]
        Gene symbol -> learned coefficient (sign matters: positive = risk-up).
    intercept : float
        Model intercept (added to weighted sum before returning score).
    normalization : Dict[str, Tuple[float, float]]
        Gene symbol -> (mean, std) used for z-score normalization in .score().
    cutoffs : Dict[str, float]
        Named cutoffs used by .stratify(). For two-group: {"threshold": v}.
        For three-group: {"low": v1, "high": v2} with v1 < v2.
    risk_labels : List[str]
        Labels in ascending-risk order, e.g. ["low", "high"] or
        ["low", "medium", "high"].
    performance : Dict[str, float]
        Summary metrics captured at fit time, e.g. {"auc": 0.87, "youden_j": 0.6}.
    panel_genes : List[str]
        Ordered list of gene symbols in the panel (matches weights.keys()).
    """

    mode: str
    weights: Dict[str, float]
    intercept: float
    normalization: Dict[str, Tuple[float, float]]
    cutoffs: Dict[str, float]
    risk_labels: List[str]
    performance: Dict[str, float] = field(default_factory=dict)
    panel_genes: List[str] = field(init=False)

    def __post_init__(self) -> None:
        if self.mode not in VALID_MODES:
            raise ValueError(
                f"Invalid mode {self.mode!r}. Must be one of: {', '.join(VALID_MODES)}."
            )
        # panel_genes derived from weights (preserves insertion order)
        self.panel_genes = list(self.weights.keys())

        # Sanity checks
        if set(self.panel_genes) != set(self.normalization.keys()):
            raise ValueError(
                "weights and normalization must have the same gene keys."
            )
        if len(self.risk_labels) < 2:
            raise ValueError("risk_labels must have at least 2 entries.")

    # -----------------------------------------------------------------------
    # Apply to new data
    # -----------------------------------------------------------------------

    def score(self, expression: pd.DataFrame) -> pd.Series:
        """
        Compute per-sample signature scores from new expression data.

        Parameters
        ----------
        expression : pd.DataFrame
            Samples x genes. Must contain every gene in self.panel_genes.
            Extra genes are ignored.

        Returns
        -------
        pd.Series
            Index = samples (from expression.index), values = scores.
        """
        missing = [g for g in self.panel_genes if g not in expression.columns]
        if missing:
            raise ValueError(
                f"Expression data is missing required panel genes: {missing}"
            )

        # Restrict to panel, in the stored order
        X = expression[self.panel_genes].astype(float)

        # Z-score normalize using stored (mean, std)
        X_norm = X.copy()
        for gene, (mu, sd) in self.normalization.items():
            # Guard against zero variance in training (extremely rare, but safe)
            denom = sd if sd > 1e-12 else 1.0
            X_norm[gene] = (X[gene] - mu) / denom

        # Weighted sum + intercept
        w = np.array([self.weights[g] for g in self.panel_genes], dtype=float)
        raw_score = X_norm.values @ w + self.intercept

        return pd.Series(raw_score, index=expression.index, name="signature_score")

    # -----------------------------------------------------------------------
    # Stratify into risk groups
    # -----------------------------------------------------------------------

    def stratify(self, scores: pd.Series) -> pd.Series:
        """
        Assign each score to a risk group based on learned cutoffs.

        Two-group (default): cutoffs = {"threshold": v}
            -> "low" if score < v, else "high"
        Three-group: cutoffs = {"low": v1, "high": v2} with v1 < v2
            -> "low" if score < v1; "medium" if v1 <= score < v2; "high" if >= v2

        Parameters
        ----------
        scores : pd.Series
            Output of self.score().

        Returns
        -------
        pd.Series
            Index = samples, values = risk group labels (categorical).
        """
        if len(self.risk_labels) == 2:
            if "threshold" not in self.cutoffs:
                raise ValueError(
                    "Two-group stratification requires cutoffs={'threshold': v}."
                )
            t = self.cutoffs["threshold"]
            labels = np.where(scores.values < t, self.risk_labels[0], self.risk_labels[1])
        elif len(self.risk_labels) == 3:
            if not {"low", "high"}.issubset(self.cutoffs):
                raise ValueError(
                    "Three-group stratification requires cutoffs={'low': v1, 'high': v2}."
                )
            lo, hi = self.cutoffs["low"], self.cutoffs["high"]
            if lo >= hi:
                raise ValueError(
                    f"cutoffs['low'] ({lo}) must be < cutoffs['high'] ({hi})."
                )
            s = scores.values
            labels = np.where(
                s < lo, self.risk_labels[0],
                np.where(s < hi, self.risk_labels[1], self.risk_labels[2]),
            )
        else:
            raise ValueError(
                f"Only 2 or 3 risk groups supported, got {len(self.risk_labels)}."
            )

        return pd.Series(
            pd.Categorical(labels, categories=self.risk_labels, ordered=True),
            index=scores.index,
            name="risk_group",
        )

    # -----------------------------------------------------------------------
    # Serialize
    # -----------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Serialize to a plain dict (useful for JSON logs or saving)."""
        return {
            "mode": self.mode,
            "weights": dict(self.weights),
            "intercept": self.intercept,
            "normalization": {g: list(v) for g, v in self.normalization.items()},
            "cutoffs": dict(self.cutoffs),
            "risk_labels": list(self.risk_labels),
            "performance": dict(self.performance),
            "panel_genes": list(self.panel_genes),
        }


# ---------------------------------------------------------------------------
# Factory: build a signature score from training data
# ---------------------------------------------------------------------------

def build_signature_score(
    X: pd.DataFrame,
    y: Union[pd.Series, np.ndarray, Sequence],
    panel_genes: List[str],
    mode: str = "diagnostic",
    *,
    time: Optional[Union[pd.Series, np.ndarray]] = None,
    event: Optional[Union[pd.Series, np.ndarray]] = None,
    risk_groups: int = 2,
    random_state: int = 42,
) -> SignatureScore:
    """
    Fit a signature score from training data.

    Parameters
    ----------
    X : pd.DataFrame
        Samples x genes. Must contain all panel_genes.
    y : array-like, required for mode='diagnostic'
        Binary labels (0/1). Ignored for prognostic mode.
    panel_genes : List[str]
        Genes to include in the signature. Order is preserved.
    mode : str
        'diagnostic' or 'prognostic'.
    time : array-like, required for mode='prognostic'
        Survival time for each sample.
    event : array-like, required for mode='prognostic'
        Event indicator (1 = event occurred, 0 = censored).
    risk_groups : int
        2 for low/high, 3 for low/medium/high. Default 2.
    random_state : int
        Seed for reproducibility.

    Returns
    -------
    SignatureScore
        Fitted score object.

    Raises
    ------
    ValueError
        For invalid mode or missing genes.
    ImportError
        If mode='prognostic' and lifelines is not installed.
    """
    if mode not in VALID_MODES:
        raise ValueError(
            f"Invalid mode {mode!r}. Must be one of: {', '.join(VALID_MODES)}."
        )

    if risk_groups not in (2, 3):
        raise ValueError(f"risk_groups must be 2 or 3, got {risk_groups}.")

    missing = [g for g in panel_genes if g not in X.columns]
    if missing:
        raise ValueError(f"X is missing required panel genes: {missing}")

    # Restrict and compute z-score parameters on training data
    X_panel = X[panel_genes].astype(float)
    means = X_panel.mean(axis=0)
    stds = X_panel.std(axis=0, ddof=0)
    # Replace zero std with 1.0 to avoid division issues (constant gene)
    stds_safe = stds.where(stds > 1e-12, other=1.0)
    X_norm = (X_panel - means) / stds_safe

    normalization = {g: (float(means[g]), float(stds[g])) for g in panel_genes}

    if mode == "diagnostic":
        return _fit_diagnostic(
            X_norm, y, panel_genes, normalization, risk_groups, random_state
        )
    else:  # prognostic
        if not _LIFELINES_AVAILABLE:
            raise ImportError(
                "Prognostic mode requires the 'lifelines' package. "
                "Install with: pip install lifelines"
            )
        if time is None or event is None:
            raise ValueError(
                "Prognostic mode requires both 'time' and 'event' arguments."
            )
        return _fit_prognostic(
            X_norm, time, event, panel_genes, normalization, risk_groups, random_state
        )


# ---------------------------------------------------------------------------
# Internal fitters
# ---------------------------------------------------------------------------

def _fit_diagnostic(
    X_norm: pd.DataFrame,
    y: Union[pd.Series, np.ndarray, Sequence],
    panel_genes: List[str],
    normalization: Dict[str, Tuple[float, float]],
    risk_groups: int,
    random_state: int,
) -> SignatureScore:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, roc_curve

    y_arr = np.asarray(y).astype(int)
    if set(np.unique(y_arr)) - {0, 1}:
        raise ValueError("Diagnostic mode requires binary (0/1) labels.")

# Fit logistic regression on normalized expression.
    # sklearn >=1.8 deprecated the `penalty` argument in favor of `l1_ratio`:
    # l1_ratio=0 means pure L2, l1_ratio=1 means pure L1. We combine l1_ratio=0
    # with a very large C to get effectively unregularized coefficients, which
    # is what we want for a signature built on a small curated panel.
    model = LogisticRegression(
        l1_ratio=0,
        C=1e10,
        solver="saga",
        max_iter=5000,
        random_state=random_state,
    )
    model.fit(X_norm.values, y_arr)

    weights = {g: float(c) for g, c in zip(panel_genes, model.coef_[0])}
    intercept = float(model.intercept_[0])

    # Compute training scores for threshold selection
    train_scores = X_norm.values @ model.coef_[0] + intercept

    # Performance metrics
    try:
        auc = float(roc_auc_score(y_arr, train_scores))
    except ValueError:
        auc = float("nan")  # e.g. only one class present

    # Youden's J optimal threshold
    fpr, tpr, thresholds = roc_curve(y_arr, train_scores)
    j_scores = tpr - fpr
    best_idx = int(np.argmax(j_scores))
    best_threshold = float(thresholds[best_idx])
    best_j = float(j_scores[best_idx])

    performance = {"auc": auc, "youden_j": best_j}

    if risk_groups == 2:
        cutoffs = {"threshold": best_threshold}
        risk_labels = ["low", "high"]
    else:
        # Three-group: use tertiles of training scores
        q33, q66 = np.quantile(train_scores, [1 / 3, 2 / 3])
        cutoffs = {"low": float(q33), "high": float(q66)}
        risk_labels = ["low", "medium", "high"]

    return SignatureScore(
        mode="diagnostic",
        weights=weights,
        intercept=intercept,
        normalization=normalization,
        cutoffs=cutoffs,
        risk_labels=risk_labels,
        performance=performance,
    )


def _fit_prognostic(
    X_norm: pd.DataFrame,
    time: Union[pd.Series, np.ndarray],
    event: Union[pd.Series, np.ndarray],
    panel_genes: List[str],
    normalization: Dict[str, Tuple[float, float]],
    risk_groups: int,
    random_state: int,
) -> SignatureScore:
    from lifelines import CoxPHFitter
    from lifelines.statistics import logrank_test

    # Build the dataframe lifelines wants
    df = X_norm.copy()
    df["_time"] = np.asarray(time, dtype=float)
    df["_event"] = np.asarray(event, dtype=int)

    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(df, duration_col="_time", event_col="_event", show_progress=False)

    coefs = cph.params_
    weights = {g: float(coefs[g]) for g in panel_genes}
    intercept = 0.0  # Cox has no intercept; linear predictor is sum of log-HRs

    # Training linear predictor
    train_scores = X_norm.values @ np.array(
        [weights[g] for g in panel_genes], dtype=float
    )

    performance = {"c_index": float(cph.concordance_index_)}

    if risk_groups == 2:
        # Pick the split that maximizes the log-rank chi-square
        sorted_scores = np.sort(np.unique(train_scores))
        best_t, best_chi = float(np.median(train_scores)), -np.inf
        event_arr = np.asarray(event, dtype=int)
        time_arr = np.asarray(time, dtype=float)
        # Search candidate thresholds (skip extremes with tiny groups)
        for t in sorted_scores[1:-1]:
            high_mask = train_scores >= t
            if high_mask.sum() < 5 or (~high_mask).sum() < 5:
                continue
            try:
                res = logrank_test(
                    time_arr[high_mask], time_arr[~high_mask],
                    event_observed_A=event_arr[high_mask],
                    event_observed_B=event_arr[~high_mask],
                )
                if res.test_statistic > best_chi:
                    best_chi = float(res.test_statistic)
                    best_t = float(t)
            except Exception:
                continue
        cutoffs = {"threshold": best_t}
        risk_labels = ["low", "high"]
        performance["logrank_chi2"] = float(best_chi) if best_chi > -np.inf else float("nan")
    else:
        q33, q66 = np.quantile(train_scores, [1 / 3, 2 / 3])
        cutoffs = {"low": float(q33), "high": float(q66)}
        risk_labels = ["low", "medium", "high"]

    return SignatureScore(
        mode="prognostic",
        weights=weights,
        intercept=intercept,
        normalization=normalization,
        cutoffs=cutoffs,
        risk_labels=risk_labels,
        performance=performance,
    )