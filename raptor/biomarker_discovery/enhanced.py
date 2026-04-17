"""
Enhanced Biomarker Analysis
===========================

Wires the five M10 enhancement modules (intent, signature score, direction
patterns, clinical metrics, ratio biomarkers) into a unified post-discovery
analysis triggered by setting ``intent`` in ``discover_biomarkers()``.

The function ``enhance_biomarker_result()`` takes a completed
``BiomarkerResult`` plus the expression data and labels, and runs every
applicable enhanced analysis, returning an ``EnhancedBiomarkerResult``.

Usage
-----
>>> from raptor.biomarker_discovery import discover_biomarkers
>>> result = discover_biomarkers(counts, metadata, intent="diagnostic")
>>> result.signature          # SignatureScore object
>>> result.direction_pattern  # DirectionPattern object
>>> result.clinical_metrics   # dict with Youden, bootstrap CI, DCA, PPV/NPV
>>> result.ratio_result       # RatioSearchResult object
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from raptor.biomarker_discovery.intent import BiomarkerIntent
from raptor.biomarker_discovery.signature_score import (
    SignatureScore,
    build_signature_score,
)
from raptor.biomarker_discovery.direction_patterns import (
    DirectionPattern,
    build_direction_pattern,
)
from raptor.biomarker_discovery.clinical_metrics import (
    ppv_npv_at_prevalence,
    bootstrap_ci,
    youdens_optimal_threshold,
    decision_curve_analysis,
    net_reclassification_improvement,
)
from raptor.biomarker_discovery.ratio_biomarkers import (
    RatioBiomarkerSearcher,
    RatioSearchResult,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Enhanced result container
# ---------------------------------------------------------------------------

@dataclass
class EnhancedBiomarkerResult:
    """
    Extended biomarker result that wraps the base ``BiomarkerResult`` and
    adds outputs from the five M10 enhancement modules.

    Access everything from the original result via ``base_result``, or
    use the convenience properties (``panel``, ``ranked_genes``, etc.)
    which delegate to the base.

    Attributes
    ----------
    base_result : BiomarkerResult
        The original discovery result (panel, rankings, CV metrics, etc.).
    intent : BiomarkerIntent
        The intent object that configured the enhanced analyses.
    signature : SignatureScore, optional
        Learned weighted risk score with per-patient scoring and
        stratification cutoffs.
    direction_pattern : DirectionPattern, optional
        Per-gene UP/DOWN pattern with fold changes and confidence.
    clinical_metrics : dict, optional
        Dict with keys: 'youdens', 'bootstrap_ci', 'ppv_npv',
        'decision_curve'. Each holds the output of the corresponding
        clinical_metrics function.
    ratio_result : RatioSearchResult, optional
        Top discriminating gene-ratio pairs.
    """

    base_result: Any   # BiomarkerResult (Any to avoid circular import)
    intent: BiomarkerIntent
    signature: Optional[SignatureScore] = None
    direction_pattern: Optional[DirectionPattern] = None
    clinical_metrics: Optional[Dict[str, Any]] = None
    ratio_result: Optional[RatioSearchResult] = None

    # -- convenience delegates to base_result --------------------------------

    @property
    def panel(self) -> List[str]:
        return self.base_result.panel

    @property
    def panel_size(self) -> int:
        return self.base_result.panel_size

    @property
    def ranked_genes(self) -> pd.DataFrame:
        return self.base_result.ranked_genes

    @property
    def classification_results(self) -> Dict:
        return self.base_result.classification_results

    @property
    def best_classifier(self) -> str:
        return self.base_result.best_classifier

    @property
    def panel_optimization(self):
        return self.base_result.panel_optimization

    # -- summary -------------------------------------------------------------

    def summary(self) -> str:
        """Full summary: base result + enhanced analyses."""
        lines = [self.base_result.summary()]
        lines.append("")
        lines.append(f"ENHANCED ANALYSIS ({self.intent.output_label})")
        lines.append("-" * 50)

        if self.signature is not None:
            lines.append(f"  Signature score: {self.signature.mode} mode, "
                         f"{len(self.signature.panel_genes)} genes")
            lines.append(f"    Risk groups: {self.signature.risk_labels}")

        if self.direction_pattern is not None:
            dp = self.direction_pattern
            lines.append(f"  Direction pattern: {dp.n_up} UP, {dp.n_down} DOWN "
                         f"({dp.n_genes} total)")

        if self.clinical_metrics is not None:
            cm = self.clinical_metrics
            if "youdens" in cm:
                y = cm["youdens"]
                lines.append(f"  Youden's optimal: threshold={y['threshold']:.3f}, "
                             f"J={y['youdens_j']:.3f}")
            if "bootstrap_ci" in cm:
                b = cm["bootstrap_ci"]
                lines.append(f"  AUC bootstrap CI: {b['ci_lower']:.3f} - "
                             f"{b['ci_upper']:.3f} (95%)")
            if "ppv_npv" in cm:
                p = cm["ppv_npv"]
                lines.append(f"  PPV at prevalence {p['prevalence']:.1%}: "
                             f"{p['ppv']:.3f}")
                lines.append(f"  NPV at prevalence {p['prevalence']:.1%}: "
                             f"{p['npv']:.3f}")

        if self.ratio_result is not None:
            rr = self.ratio_result
            lines.append(f"  Ratio biomarkers: {rr.n_found} pairs found "
                         f"(top AUC={rr.best_pair.auc:.3f})"
                         if rr.best_pair else
                         f"  Ratio biomarkers: no pairs passed threshold")

        lines.append("")
        return "\n".join(lines)

    # -- save delegates to base + saves enhanced outputs ---------------------

    def save(self, output_dir):
        """Save base result plus enhanced analysis outputs."""
        from pathlib import Path
        output_dir = Path(output_dir)

        # Save base result
        self.base_result.save(output_dir)

        # Save enhanced outputs
        enhanced_dir = output_dir / "enhanced"
        enhanced_dir.mkdir(parents=True, exist_ok=True)

        if self.direction_pattern is not None:
            self.direction_pattern.to_dataframe().to_csv(
                enhanced_dir / "direction_pattern.csv"
            )

        if self.clinical_metrics is not None:
            if "decision_curve" in self.clinical_metrics:
                self.clinical_metrics["decision_curve"].to_csv(
                    enhanced_dir / "decision_curve.csv", index=False
                )

        if self.ratio_result is not None:
            self.ratio_result.to_dataframe().to_csv(
                enhanced_dir / "ratio_biomarkers.csv", index=False
            )

        if self.signature is not None:
            import json
            with open(enhanced_dir / "signature_score.json", "w",
                       encoding="utf-8") as f:
                json.dump(self.signature.to_dict(), f, indent=2, default=str)

        logger.info(f"Enhanced results saved to: {enhanced_dir}/")


# ---------------------------------------------------------------------------
# Main integration function
# ---------------------------------------------------------------------------

def enhance_biomarker_result(
    base_result: Any,
    expression: pd.DataFrame,
    labels: np.ndarray,
    group_names: Tuple[str, str],
    intent: Union[str, "BiomarkerIntent"],
    prevalence: float = 0.05,
    ratio_top_k: int = 10,
    ratio_min_auc: float = 0.70,
    random_state: int = 42,
    verbose: bool = True,
) -> EnhancedBiomarkerResult:
    """
    Run enhanced analyses on a completed BiomarkerResult.

    This is called automatically by ``discover_biomarkers()`` when
    ``intent`` is set, but can also be called standalone on any
    existing BiomarkerResult.

    Parameters
    ----------
    base_result : BiomarkerResult
        Completed discovery result with panel, classifiers, etc.
    expression : pd.DataFrame (samples × genes)
        The same expression matrix used for discovery (log2-CPM).
    labels : np.ndarray
        Binary labels (0/1), same order as expression rows.
    group_names : tuple of (str, str)
        (baseline_group_name, reference_group_name) where baseline
        maps to label 0 and reference maps to label 1.
    intent : str or BiomarkerIntent
        Intent mode (e.g. 'diagnostic') or a pre-built BiomarkerIntent.
    prevalence : float, default 0.05
        Disease prevalence for PPV/NPV calculation.
    ratio_top_k : int, default 10
        Max ratio pairs to return.
    ratio_min_auc : float, default 0.70
        Minimum AUC for ratio pairs.
    random_state : int, default 42
    verbose : bool, default True

    Returns
    -------
    EnhancedBiomarkerResult
    """
    if isinstance(intent, str):
        intent_obj = BiomarkerIntent(intent)
    else:
        intent_obj = intent

    baseline_group, reference_group = group_names
    panel_genes = base_result.panel

    if verbose:
        logger.info("")
        logger.info(f"ENHANCED ANALYSIS: {intent_obj.output_label}")
        logger.info("-" * 50)

    # -- 1. Signature Score --------------------------------------------------
    sig = None
    if intent_obj.intent in ("diagnostic", "prognostic"):
        try:
            mode = intent_obj.intent  # 'diagnostic' or 'prognostic'
            # Prognostic mode requires survival data which we don't have
            # in the standard discover_biomarkers flow. Fall back to
            # diagnostic mode with a note.
            if mode == "prognostic":
                if verbose:
                    logger.info("  Signature: falling back to diagnostic mode "
                                "(survival data not available in this flow)")
                mode = "diagnostic"

            sig = build_signature_score(
                X=expression,
                y=labels,
                panel_genes=panel_genes,
                mode=mode,
            )
            if verbose:
                logger.info(f"  Signature score: {mode} mode, "
                            f"{len(sig.genes)} genes, "
                            f"cutoffs={sig.cutoffs}")
        except Exception as e:
            if verbose:
                logger.warning(f"  Signature score failed: {e}")

    # -- 2. Direction Pattern ------------------------------------------------
    dp = None
    if intent_obj.intent in ("diagnostic", "prognostic", "predictive",
                              "exploratory", "translational"):
        try:
            # Build labels as a string Series for direction pattern
            str_labels = pd.Series(
                [reference_group if lbl == 1 else baseline_group
                 for lbl in labels],
                index=expression.index,
            )
            dp = build_direction_pattern(
                expression=expression[panel_genes],
                labels=str_labels,
                reference_group=reference_group,
                baseline_group=baseline_group,
                p_threshold=0.05,
                fc_threshold=0.0,  # panel genes already passed selection
            )
            if verbose:
                logger.info(f"  Direction pattern: {dp.n_up} UP, "
                            f"{dp.n_down} DOWN ({dp.n_genes} genes)")
        except Exception as e:
            if verbose:
                logger.warning(f"  Direction pattern failed: {e}")

    # -- 3. Clinical Metrics -------------------------------------------------
    clin = None
    if intent_obj.intent in ("diagnostic", "prognostic", "predictive"):
        try:
            clin = {}

            # Get predicted probabilities from the best classifier
            best_clf_name = base_result.best_classifier
            best_clf_result = base_result.classification_results.get(best_clf_name)
            y_score = None

            if best_clf_result and best_clf_result.trained_model is not None:
                model = best_clf_result.trained_model
                X_panel = expression[panel_genes]
                if hasattr(model, "predict_proba"):
                    y_score = model.predict_proba(X_panel)[:, 1]
                elif hasattr(model, "decision_function"):
                    raw = model.decision_function(X_panel)
                    # Sigmoid transform for SVM
                    y_score = 1.0 / (1.0 + np.exp(-raw))

            if y_score is not None:
                y_true = labels.astype(int)

                # Youden's optimal threshold
                clin["youdens"] = youdens_optimal_threshold(y_true, y_score)
                if verbose:
                    logger.info(f"  Youden's J: {clin['youdens']['youdens_j']:.3f} "
                                f"(threshold={clin['youdens']['threshold']:.3f})")

                # Bootstrap CI for AUC
                clin["bootstrap_ci"] = bootstrap_ci(
                    y_true, y_score, seed=random_state,
                )
                if verbose:
                    b = clin["bootstrap_ci"]
                    logger.info(f"  AUC 95% CI: [{b['ci_lower']:.3f}, "
                                f"{b['ci_upper']:.3f}]")

                # PPV/NPV at prevalence
                sens = clin["youdens"]["sensitivity"]
                spec = clin["youdens"]["specificity"]
                if sens > 0 and spec > 0:
                    clin["ppv_npv"] = ppv_npv_at_prevalence(
                        sensitivity=sens,
                        specificity=spec,
                        prevalence=prevalence,
                    )
                    if verbose:
                        p = clin["ppv_npv"]
                        logger.info(f"  PPV at {prevalence:.1%} prevalence: "
                                    f"{p['ppv']:.3f}")

                # Decision curve analysis
                clin["decision_curve"] = decision_curve_analysis(
                    y_true, y_score,
                )
                if verbose:
                    logger.info("  Decision curve analysis: computed")

            else:
                if verbose:
                    logger.warning("  Clinical metrics: no trained model with "
                                   "predict_proba available")

        except Exception as e:
            if verbose:
                logger.warning(f"  Clinical metrics failed: {e}")
            clin = None

    # -- 4. Ratio Biomarkers -------------------------------------------------
    rr = None
    if intent_obj.intent in ("diagnostic", "prognostic", "predictive",
                              "exploratory", "translational"):
        try:
            # Use panel genes as candidates (already validated)
            # But also search among all genes for ratio pairs
            str_labels = pd.Series(
                [reference_group if lbl == 1 else baseline_group
                 for lbl in labels],
                index=expression.index,
            )
            searcher = RatioBiomarkerSearcher(
                top_k=ratio_top_k,
                min_auc=ratio_min_auc,
                candidate_genes=panel_genes if len(panel_genes) >= 2 else None,
            )
            rr = searcher.search(
                expression, str_labels,
                reference_group=reference_group,
                baseline_group=baseline_group,
            )
            if verbose:
                if rr.best_pair:
                    logger.info(f"  Ratio biomarkers: {rr.n_found} pairs "
                                f"(best {rr.best_pair.name} AUC={rr.best_pair.auc:.3f})")
                else:
                    logger.info("  Ratio biomarkers: no pairs above threshold")
        except Exception as e:
            if verbose:
                logger.warning(f"  Ratio biomarkers failed: {e}")

    # -- Assemble enhanced result --------------------------------------------

    return EnhancedBiomarkerResult(
        base_result=base_result,
        intent=intent_obj,
        signature=sig,
        direction_pattern=dp,
        clinical_metrics=clin,
        ratio_result=rr,
    )