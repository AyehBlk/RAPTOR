"""
Clinical Metrics Module
=======================

Functions that answer "how useful is this biomarker in the clinic?"

Five analyses are provided:

1. **ppv_npv_at_prevalence** — Positive/Negative predictive values adjusted
   for a realistic disease prevalence (Bayes' theorem).
2. **bootstrap_ci** — Non-parametric bootstrap confidence intervals for any
   scalar metric (AUC, sensitivity, etc.).
3. **youdens_optimal_threshold** — The biomarker cutoff that maximises
   sensitivity + specificity − 1 (Youden, 1950).
4. **decision_curve_analysis** — Net benefit across a range of threshold
   probabilities (Vickers & Elkin, 2006).
5. **net_reclassification_improvement** — How many patients move to a more
   correct risk category when a new biomarker is added (Pencina et al., 2008).

All functions are self-contained: they accept numpy arrays or pandas objects,
return plain dicts or DataFrames, and have no side effects.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from sklearn.metrics import roc_curve, roc_auc_score


# ═══════════════════════════════════════════════════════════════════════════
# 1. PPV / NPV at prevalence
# ═══════════════════════════════════════════════════════════════════════════

def ppv_npv_at_prevalence(
    sensitivity: float,
    specificity: float,
    prevalence: float,
) -> Dict[str, float]:
    """
    Compute Positive and Negative Predictive Values at a given prevalence.

    In a case-control study, the raw PPV/NPV reflect the study's 50/50
    mix, not the real world.  This function applies Bayes' theorem to
    translate sensitivity and specificity into clinically meaningful
    predictive values at a realistic disease prevalence.

    Parameters
    ----------
    sensitivity : float in (0, 1]
        True positive rate (TP / (TP + FN)).
    specificity : float in (0, 1]
        True negative rate (TN / (TN + FP)).
    prevalence : float in (0, 1)
        Disease prevalence in the target population.

    Returns
    -------
    dict with keys: ppv, npv, prevalence, sensitivity, specificity
    """
    for name, val in [("sensitivity", sensitivity), ("specificity", specificity)]:
        if not 0 < val <= 1:
            raise ValueError(f"{name} must be in (0, 1], got {val}.")
    if not 0 < prevalence < 1:
        raise ValueError(f"prevalence must be in (0, 1), got {prevalence}.")

    # Bayes' theorem
    ppv = (sensitivity * prevalence) / (
        sensitivity * prevalence + (1 - specificity) * (1 - prevalence)
    )
    npv = (specificity * (1 - prevalence)) / (
        specificity * (1 - prevalence) + (1 - sensitivity) * prevalence
    )

    return {
        "ppv": ppv,
        "npv": npv,
        "prevalence": prevalence,
        "sensitivity": sensitivity,
        "specificity": specificity,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 2. Bootstrap confidence intervals
# ═══════════════════════════════════════════════════════════════════════════

def bootstrap_ci(
    y_true: np.ndarray,
    y_score: np.ndarray,
    metric_fn: Optional[Callable] = None,
    n_bootstrap: int = 2000,
    ci: float = 0.95,
    seed: int = 42,
) -> Dict[str, float]:
    """
    Non-parametric bootstrap confidence interval for a metric.

    Parameters
    ----------
    y_true : array-like, shape (n,)
        Binary ground-truth labels (0/1).
    y_score : array-like, shape (n,)
        Continuous scores (e.g. predicted probabilities).
    metric_fn : callable, optional
        Function(y_true, y_score) -> float.  Defaults to
        ``sklearn.metrics.roc_auc_score``.
    n_bootstrap : int, default 2000
        Number of bootstrap resamples.
    ci : float, default 0.95
        Confidence level (e.g. 0.95 for 95 % CI).
    seed : int, default 42
        Random seed for reproducibility.

    Returns
    -------
    dict with keys: point_estimate, ci_lower, ci_upper, ci_level, n_bootstrap
    """
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)

    if len(y_true) != len(y_score):
        raise ValueError(
            f"y_true length ({len(y_true)}) != y_score length ({len(y_score)})."
        )
    if not 0.50 <= ci < 1.0:
        raise ValueError(f"ci must be in [0.50, 1.0), got {ci}.")

    if metric_fn is None:
        metric_fn = roc_auc_score

    point = float(metric_fn(y_true, y_score))

    rng = np.random.default_rng(seed)
    n = len(y_true)
    estimates: List[float] = []

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        bt_true = y_true[idx]
        bt_score = y_score[idx]
        # Skip degenerate resamples (only one class present)
        if len(np.unique(bt_true)) < 2:
            continue
        try:
            estimates.append(float(metric_fn(bt_true, bt_score)))
        except Exception:
            continue

    if len(estimates) < 100:
        raise ValueError(
            f"Only {len(estimates)} valid bootstrap resamples out of "
            f"{n_bootstrap}. Data may be too small or too imbalanced."
        )

    alpha = 1.0 - ci
    lower = float(np.percentile(estimates, 100 * alpha / 2))
    upper = float(np.percentile(estimates, 100 * (1 - alpha / 2)))

    return {
        "point_estimate": point,
        "ci_lower": lower,
        "ci_upper": upper,
        "ci_level": ci,
        "n_bootstrap": len(estimates),
    }


# ═══════════════════════════════════════════════════════════════════════════
# 3. Youden's optimal threshold
# ═══════════════════════════════════════════════════════════════════════════

def youdens_optimal_threshold(
    y_true: np.ndarray,
    y_score: np.ndarray,
) -> Dict[str, float]:
    """
    Find the threshold maximising Youden's J = sensitivity + specificity − 1.

    Parameters
    ----------
    y_true : array-like, shape (n,)
        Binary labels (0/1).
    y_score : array-like, shape (n,)
        Continuous scores.

    Returns
    -------
    dict with keys: threshold, sensitivity, specificity, youdens_j, auc
    """
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)

    if len(np.unique(y_true)) < 2:
        raise ValueError("y_true must contain both 0 and 1.")

    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    j_scores = tpr + (1 - fpr) - 1          # = tpr - fpr = Youden's J
    best_idx = int(np.argmax(j_scores))

    auc = float(roc_auc_score(y_true, y_score))

    return {
        "threshold": float(thresholds[best_idx]),
        "sensitivity": float(tpr[best_idx]),
        "specificity": float(1 - fpr[best_idx]),
        "youdens_j": float(j_scores[best_idx]),
        "auc": auc,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 4. Decision Curve Analysis
# ═══════════════════════════════════════════════════════════════════════════

def decision_curve_analysis(
    y_true: np.ndarray,
    y_score: np.ndarray,
    thresholds: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Decision curve analysis (Vickers & Elkin, 2006).

    Computes net benefit for the biomarker model, a "treat all" strategy,
    and a "treat none" strategy across a range of threshold probabilities.

    Parameters
    ----------
    y_true : array-like, shape (n,)
        Binary outcomes (0/1).
    y_score : array-like, shape (n,)
        Predicted probabilities in [0, 1].
    thresholds : array-like, optional
        Threshold probabilities to evaluate. Defaults to 0.01–0.99 in
        steps of 0.01.

    Returns
    -------
    pd.DataFrame with columns:
        threshold, net_benefit_model, net_benefit_treat_all, net_benefit_treat_none
    """
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)

    if len(np.unique(y_true)) < 2:
        raise ValueError("y_true must contain both 0 and 1.")

    if thresholds is None:
        thresholds = np.arange(0.01, 1.0, 0.01)
    else:
        thresholds = np.asarray(thresholds, dtype=float)

    n = len(y_true)
    event_rate = y_true.mean()

    rows: List[Dict[str, float]] = []
    for pt in thresholds:
        if pt <= 0 or pt >= 1:
            continue

        # Model net benefit
        predicted_positive = (y_score >= pt)
        tp = int(np.sum(predicted_positive & (y_true == 1)))
        fp = int(np.sum(predicted_positive & (y_true == 0)))
        weight = pt / (1 - pt)
        nb_model = (tp / n) - (fp / n) * weight

        # Treat all: everyone is "positive"
        nb_all = event_rate - (1 - event_rate) * weight

        rows.append({
            "threshold": float(pt),
            "net_benefit_model": nb_model,
            "net_benefit_treat_all": nb_all,
            "net_benefit_treat_none": 0.0,
        })

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════════════
# 5. Net Reclassification Improvement
# ═══════════════════════════════════════════════════════════════════════════

def net_reclassification_improvement(
    y_true: np.ndarray,
    risk_old: np.ndarray,
    risk_new: np.ndarray,
    cutoffs: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """
    Category-based Net Reclassification Improvement (Pencina et al., 2008).

    Measures how many patients are reclassified into more correct risk
    categories when switching from an old model/biomarker to a new one.

    Parameters
    ----------
    y_true : array-like, shape (n,)
        Binary outcomes (0 = non-event, 1 = event).
    risk_old : array-like, shape (n,)
        Predicted probabilities from the old/baseline model.
    risk_new : array-like, shape (n,)
        Predicted probabilities from the new model (with biomarker).
    cutoffs : list of float, optional
        Risk category boundaries. Default: [0.2, 0.5] producing three
        categories (low / medium / high risk).

    Returns
    -------
    dict with keys:
        nri, nri_events, nri_nonevents, z_score, p_value,
        n_events, n_nonevents, events_up, events_down,
        nonevents_up, nonevents_down
    """
    y_true = np.asarray(y_true, dtype=int)
    risk_old = np.asarray(risk_old, dtype=float)
    risk_new = np.asarray(risk_new, dtype=float)

    if not (len(y_true) == len(risk_old) == len(risk_new)):
        raise ValueError("All input arrays must have the same length.")
    if len(np.unique(y_true)) < 2:
        raise ValueError("y_true must contain both 0 and 1.")

    if cutoffs is None:
        cutoffs = [0.2, 0.5]
    cutoffs = sorted(cutoffs)

    def _categorize(probs: np.ndarray) -> np.ndarray:
        cats = np.zeros(len(probs), dtype=int)
        for c in cutoffs:
            cats += (probs >= c).astype(int)
        return cats

    cat_old = _categorize(risk_old)
    cat_new = _categorize(risk_new)

    events = y_true == 1
    nonevents = y_true == 0

    n_events = int(events.sum())
    n_nonevents = int(nonevents.sum())

    if n_events == 0 or n_nonevents == 0:
        raise ValueError("Need at least one event and one non-event.")

    # Events: moving UP in category is good
    events_up = int(np.sum(cat_new[events] > cat_old[events]))
    events_down = int(np.sum(cat_new[events] < cat_old[events]))
    nri_events = (events_up - events_down) / n_events

    # Non-events: moving DOWN in category is good
    nonevents_up = int(np.sum(cat_new[nonevents] > cat_old[nonevents]))
    nonevents_down = int(np.sum(cat_new[nonevents] < cat_old[nonevents]))
    nri_nonevents = (nonevents_down - nonevents_up) / n_nonevents

    nri = nri_events + nri_nonevents

   # Asymptotic z-test (Pencina 2008, exact multinomial variance)
    p_up_e = events_up / n_events
    p_dn_e = events_down / n_events
    var_events = (p_up_e + p_dn_e - (p_up_e - p_dn_e) ** 2) / n_events

    p_up_ne = nonevents_up / n_nonevents
    p_dn_ne = nonevents_down / n_nonevents
    var_nonevents = (p_up_ne + p_dn_ne - (p_up_ne - p_dn_ne) ** 2) / n_nonevents

    se = np.sqrt(var_events + var_nonevents) if (var_events + var_nonevents) > 0 else 1e-10
    z = nri / se
    p_value = float(2 * sp_stats.norm.sf(abs(z)))

    return {
        "nri": nri,
        "nri_events": nri_events,
        "nri_nonevents": nri_nonevents,
        "z_score": z,
        "p_value": p_value,
        "n_events": n_events,
        "n_nonevents": n_nonevents,
        "events_up": events_up,
        "events_down": events_down,
        "nonevents_up": nonevents_up,
        "nonevents_down": nonevents_down,
    }