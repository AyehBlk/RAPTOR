"""
RAPTOR Dashboard - Biomarker Discovery (Module 10)

End-to-end biomarker discovery from count matrix + metadata, with optional
DE / ensemble inputs and the full M10 enhancement suite (intent, signature
score, direction patterns, clinical metrics, ratio biomarkers).

Author: Ayeh Bolouki
Version: 2.2.2
"""

import io
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# -----------------------------------------------------------------------------
# Imports (with graceful fallbacks)
# -----------------------------------------------------------------------------

try:
    from raptor.biomarker_discovery import (
        discover_biomarkers,
        validate_biomarkers,
        BiomarkerResult,
        BiomarkerIntent,
        VALID_INTENTS,
        apply_ratios,
        build_direction_pattern,
        _BORUTA_AVAILABLE,
        _MRMR_AVAILABLE,
        _SHAP_AVAILABLE,
    )
    BIOMARKER_AVAILABLE = True
    biomarker_import_error = None
except ImportError as e:
    BIOMARKER_AVAILABLE = False
    biomarker_import_error = str(e)
    VALID_INTENTS = (
        "diagnostic", "prognostic", "predictive",
        "monitoring", "exploratory", "translational",
    )
    _BORUTA_AVAILABLE = False
    _MRMR_AVAILABLE = False
    _SHAP_AVAILABLE = False

# Optional ensemble result type (for type-checking downstream consensus genes)
try:
    from raptor.ensemble import EnsembleResult  # noqa: F401
    ENSEMBLE_AVAILABLE = True
except ImportError:
    ENSEMBLE_AVAILABLE = False


# -----------------------------------------------------------------------------
# Helpers: extract best classifier info from either BiomarkerResult or
# EnhancedBiomarkerResult (the latter wraps the former in .base_result)
# -----------------------------------------------------------------------------

def _base(result):
    """Return the underlying BiomarkerResult regardless of whether result is
    plain or enhanced."""
    return result.base_result if hasattr(result, "base_result") else result


def _best_clf_name(result) -> str:
    b = _base(result)
    return b.best_classifier or "(none)"


def _best_clf_auc(result) -> float:
    b = _base(result)
    name = b.best_classifier
    if name and name in b.classification_results:
        return float(b.classification_results[name].auc)
    return 0.0


# -----------------------------------------------------------------------------
# Dashboard polish helpers (M6 + M4 diagnostic rendering)
# -----------------------------------------------------------------------------
# Thresholds defining the optimism-banner decision tree. Kept as
# module-level constants so future calibration studies can tune them
# without grep-hunting through render logic. These values are the
# outcome of the post-M6 / post-M4 dashboard-scoping discussion:
# "chance level" = AUC < 0.65, "strong signal" = AUC >= 0.75.
OPTIMISM_AUC_STRONG_SIGNAL = 0.75
OPTIMISM_AUC_CHANCE_LEVEL = 0.65

# Phi thresholds per Nogueira et al. 2018 JMLR 18:174 benchmark scale.
# Mirrored in core.py (_NOGUEIRA_EXCELLENT_THRESHOLD / _NOGUEIRA_POOR_THRESHOLD);
# duplicated here for the color mapper to avoid a reach-into-private-names
# import. If those ever change, update both.
STABILITY_EXCELLENT_THRESHOLD = 0.75
STABILITY_INTERMEDIATE_THRESHOLD = 0.50


def _stability_color(phi: float) -> str:
    """Traffic-light color for a Nogueira stability value.

    Returns 'green' (phi >= 0.75), 'amber' (0.5 <= phi < 0.75), or
    'red' (phi < 0.5). This is narrower than the Nogueira benchmark
    (which draws its 'poor' line at 0.4) because the dashboard is a
    user-facing signal and a wider red zone is the safer reporting
    choice for small cohorts.
    """
    if phi >= STABILITY_EXCELLENT_THRESHOLD:
        return "green"
    if phi >= STABILITY_INTERMEDIATE_THRESHOLD:
        return "amber"
    return "red"


def _stability_interpretation(benchmark_label: str) -> str:
    """One-sentence interpretation text matched to a benchmark label.

    benchmark_label comes from PanelStabilityResult.benchmark_label
    ('excellent' | 'intermediate' | 'poor').
    """
    if benchmark_label == "excellent":
        return (
            "Per-fold panels agree strongly. Panel composition is robust "
            "to sampling variation."
        )
    if benchmark_label == "intermediate":
        return (
            "Per-fold panels show partial agreement. Real signal likely "
            "exists, but specific panel membership has legitimate "
            "uncertainty. Typical at n < 30. Consider independent "
            "cohort validation."
        )
    return (
        "Per-fold panels disagree substantially. Different data subsets "
        "produce different panels. This indicates weak or absent signal, "
        "or severe small-cohort variability."
    )


def _optimism_banner_level(cv_auc: float, gap: float) -> Tuple[str, str]:
    """Decision tree for the Clinical Metrics optimism-gap banner.

    Returns (level, text) where level is 'green' / 'amber' / 'red'.
    Encodes the six-branch decision tree from the M6+M4 dashboard
    polish scoping doc, keyed on both |gap| AND CV AUC to prevent the
    pre-fix reporting bug where a small gap on chance-level AUC
    silently read as 'genuine signal'.
    """
    abs_gap = abs(gap)
    # Large negative gap — atypical but possible with small n + tree-based clfs
    if gap <= -0.10:
        return (
            "amber",
            f"Atypical negative optimism gap ({gap:+.3f}). Training AUC "
            f"is substantially below CV AUC. This can occur with small "
            f"cohorts and tree-based classifiers; review per-classifier "
            f"AUC values.",
        )
    # Large positive gap — overfitting
    if gap >= 0.10:
        return (
            "red",
            f"Large optimism gap ({gap:+.3f}). Training AUC "
            f"substantially exceeds CV AUC. Model fits training data "
            f"better than it generalizes. Consider larger n or stronger "
            f"regularization.",
        )
    # Moderate gap — regardless of AUC level
    if abs_gap >= 0.05:
        return (
            "amber",
            f"Moderate optimism gap ({gap:+.3f}). Training AUC exceeds "
            f"CV AUC by a visible margin. Use CV AUC as the honest "
            f"generalization estimate.",
        )
    # Gap is small (|gap| < 0.05). Now the AUC level decides the message.
    if cv_auc >= OPTIMISM_AUC_STRONG_SIGNAL:
        return (
            "green",
            f"CV and training AUC agree closely ({gap:+.3f}), consistent "
            f"with a panel capturing genuine signal.",
        )
    if cv_auc < OPTIMISM_AUC_CHANCE_LEVEL:
        return (
            "amber",
            f"CV and training AUC agree at chance level (CV AUC="
            f"{cv_auc:.3f}, gap {gap:+.3f}), consistent with absent or "
            f"very weak signal. Review Panel stability and clinical "
            f"metrics to confirm.",
        )
    # Borderline zone: 0.65 <= AUC < 0.75 with small gap
    return (
        "amber",
        f"CV and training AUC agree at a borderline level (CV AUC="
        f"{cv_auc:.3f}, gap {gap:+.3f}). Panel captures some signal but "
        f"clinical utility is uncertain. Review PPV at your target "
        f"prevalence.",
    )


def _m4_chance_expectation(alpha: float, n_genes: int) -> int:
    """Expected count of genes passing alpha under the uniform-p null.

    Under H0 the p-values are uniform on [0, 1], so the expected count
    with p < alpha is alpha * n_genes. Rounded to the nearest integer
    for display.
    """
    return int(round(alpha * n_genes))


def _render_stability_badge_markdown(
    phi: float, ci: Tuple[float, float], label: str,
) -> str:
    """Compact colored-text markdown for the Panel-tab header metric.

    Returns a single markdown line with the Phi value, its CI, and the
    benchmark label, colored by the traffic-light mapping. Used in the
    fifth overview-tile slot so the header reads at a glance.
    """
    color_map = {"green": "#2E7D32", "amber": "#F57C00", "red": "#C62828"}
    color = color_map[_stability_color(phi)]
    lo, hi = ci
    return (
        f"<div style='color:{color}; font-size:1.4em; font-weight:600; "
        f"line-height:1.2'>"
        f"Φ = {phi:.3f}"
        f"</div>"
        f"<div style='color:{color}; font-size:0.85em; opacity:0.85'>"
        f"[{lo:.3f}, {hi:.3f}] &middot; {label}"
        f"</div>"
    )


# -----------------------------------------------------------------------------
# M2 helper — resolve optimism-diagnostic tile values from the result.
# -----------------------------------------------------------------------------
# Why this exists: pre-M2, the dashboard read clin["bootstrap_ci"]
# ["point_estimate"] and labeled it "Training-data AUC" on the
# optimism-diagnostic tile. Under M1's OOF path, bootstrap_ci is
# computed on (oof_true, oof_prob), so that point estimate is actually
# the pooled OOF AUC — a different CV estimate, not a training-data
# estimate. The resulting "gap" was ~0 on honest signals (two OOF
# estimates of the same quantity) and flipped sign on small-n
# tree-based classifiers via the M5 tiebreak (producing the confusing
# "atypical negative gap" banner on 10a before M2).
#
# The fix: enhanced.py already computes a proper training AUC at its
# optimism_gap block (line ~484), comparing the OOF pooled AUC against
# a full-data re-fit on the final panel. M2 routes the dashboard to
# that purpose-built struct instead.


def _resolve_optimism_tiles(clin, base) -> Tuple[
    Optional[float], Optional[float], Optional[float], str,
]:
    """Resolve the three optimism-diagnostic tile values plus a source tag.

    Returns (cv_auc, train_auc, gap, source) where source is:

        'optimism_gap'  -- M1+ path. Reads clin['optimism_gap'] which
                           enhanced.py populates with a correctly-paired
                           comparison: cv_auc = pooled OOF AUC;
                           train_auc = AUC from a full-data fit on the
                           final panel; gap = their difference. This
                           is the real Harrell-style optimism gap.
        'fallback'      -- Pre-M1 path. No optimism_gap key on clin
                           (older pickled result generated before the
                           optimism-gap instrumentation was added).
                           Approximates using
                           classification_results[best].auc as cv_auc
                           and bootstrap_ci['point_estimate'] as
                           train_auc. Under M1 both of these are OOF
                           quantities so the gap is noise rather than
                           real optimism; the dashboard shows a
                           caption flagging this.
        'unavailable'   -- Neither source available; all three
                           returned values are None. The optimism
                           diagnostic block renders nothing.
    """
    # --- M1+ path: purpose-built optimism_gap struct ---
    og = clin.get("optimism_gap") if clin else None
    if og is not None:
        cv_auc = og.get("cv_auc_oof")
        train_auc = og.get("training_auc")
        gap = og.get("gap")
        if cv_auc is not None and train_auc is not None and gap is not None:
            return (cv_auc, train_auc, gap, "optimism_gap")

    # --- Fallback: older pickled result ---
    if clin and "bootstrap_ci" in clin and getattr(base, 'best_classifier', None):
        best_res = base.classification_results.get(base.best_classifier)
        if best_res is not None:
            cv_auc = getattr(best_res, "auc", None)
            train_auc = clin["bootstrap_ci"].get("point_estimate")
            if cv_auc is not None and train_auc is not None:
                return (cv_auc, train_auc, train_auc - cv_auc, "fallback")

    return (None, None, None, "unavailable")


# -----------------------------------------------------------------------------
# Small-n advisory — surface the run-to-run variance scaling with cohort size.
# -----------------------------------------------------------------------------
# Why this exists: at small n, single-seed CV AUC point estimates have
# substantial run-to-run variance. Empirically (see seed-stability
# investigation, April 2026):
#   * n=10 cohort: CV AUC range across 5 seeds = 0.42, banner spans
#     green/amber/red on the same data. Phi stays "poor" across runs;
#     the 3 ground-truth signal genes appear in 5/5 panels.
#   * n=60 cohort with real signal: CV AUC range 0.000, all diagnostics
#     converge.
# Run-to-run variance scales roughly as 1/sqrt(n) for fold splits, so
# the cutoffs below are calibrated to the regimes where the user
# should respectively trust, partially trust, or not trust a single
# CV AUC reading.
#
# Tiers, with rationale:
#   strong  (n < 20):  CV AUC point estimates are unreliable; report
#                      panel composition and Phi as primary diagnostics.
#                      Multi-seed verification strongly recommended.
#   moderate (20 <= n < 50):  Meaningful variance possible; multi-seed
#                             verification recommended.
#   none    (n >= 50): No advisory; single-seed reading is reasonably
#                      trustworthy.
SMALL_N_THRESHOLD_STRONG = 20
SMALL_N_THRESHOLD_MODERATE = 50


def _small_n_advisory(n_samples: int) -> Optional[Tuple[str, str]]:
    """Return (level, message) for a tiered advisory or None.

    level is 'strong' or 'moderate' so the dashboard can route to
    st.warning vs st.info appropriately.
    """
    if n_samples < SMALL_N_THRESHOLD_STRONG:
        return (
            "strong",
            f"**Small-cohort advisory (n={n_samples}).** At this "
            f"cohort size, single-seed cross-validation AUC point "
            f"estimates have substantial run-to-run variance — "
            f"empirically up to a 0.4 spread on the same data across "
            f"5 random seeds. Report **panel composition** (which "
            f"genes appear) and **panel stability Φ** as primary "
            f"diagnostics; treat single-seed CV AUC as one draw from "
            f"a wide distribution. Use the **\"Verify across 3 "
            f"seeds\"** button on the Panel tab before reporting "
            f"results.",
        )
    if n_samples < SMALL_N_THRESHOLD_MODERATE:
        return (
            "moderate",
            f"**Cohort-size advisory (n={n_samples}).** CV AUC point "
            f"estimates may show meaningful run-to-run variance at "
            f"this cohort size. Multi-seed verification (Panel tab "
            f"\"Verify across 3 seeds\" button) is recommended before "
            f"reporting clinical performance.",
        )
    return None


# -----------------------------------------------------------------------------
# Multi-seed verification — invariants vs variants summary.
# -----------------------------------------------------------------------------
# UX philosophy: surface the interpretation, not just the data. Users
# at small n need to be told which diagnostics are stable across seeds
# (trust them) and which are not (treat with caution). This helper
# computes the invariant/variant split so the dashboard can render it
# directly under the per-seed table without requiring the user to scan
# columns and infer.

# Default seed list for the "Verify across 3 seeds" button. Three
# additional seeds is enough to expose typical variance at small n
# (the seed-stability investigation used 5 for paper figures, but 3 is
# sufficient for a verification UX). The seed values themselves are
# arbitrary distinct integers; the user's current discovery seed is
# always also included if it isn't one of these three.
DEFAULT_VERIFICATION_SEEDS = [0, 1, 456]


def _summarize_seed_runs(rows: list) -> dict:
    """Compute invariants and variants across a list of seed-run dicts.

    Each row is expected to have keys matching run_one_seed's return
    shape: 'cv_auc', 'gap', 'banner', 'phi', 'phi_label', 'panel'
    (comma-separated string).

    Returns a dict the dashboard can consume:
        cv_auc_range:          (min, max)
        gap_range:             (min, max)
        phi_range:             (min, max)
        invariant_phi_label:   None if phi_label varies, else the
                               common label (e.g. 'poor', 'excellent')
        invariant_genes:       sorted list of genes appearing in EVERY
                               run's panel
        variant_genes:         sorted list of genes appearing in some
                               but not all runs
        invariant_banner:      None if banner varies, else the common
                               banner string
        n_runs:                len(rows)
    """
    n = len(rows)
    if n == 0:
        return {
            "cv_auc_range": (float("nan"), float("nan")),
            "gap_range": (float("nan"), float("nan")),
            "phi_range": (float("nan"), float("nan")),
            "invariant_phi_label": None,
            "invariant_genes": [],
            "variant_genes": [],
            "invariant_banner": None,
            "n_runs": 0,
        }

    cv_aucs = [r["cv_auc"] for r in rows if r.get("cv_auc") is not None]
    gaps = [r["gap"] for r in rows if r.get("gap") is not None]
    phis = [r["phi"] for r in rows if r.get("phi") is not None]
    phi_labels = [r["phi_label"] for r in rows if r.get("phi_label")]
    banners = [r["banner"] for r in rows if r.get("banner")]

    # Invariant phi label: same value across all runs?
    invariant_phi_label = (
        phi_labels[0] if len(set(phi_labels)) == 1 and phi_labels else None
    )
    # Invariant banner: same value across all runs?
    invariant_banner = (
        banners[0] if len(set(banners)) == 1 and banners else None
    )

    # Panel intersection: genes in EVERY panel = invariant
    panel_lists = [
        set(r["panel"].split(",")) for r in rows if r.get("panel")
    ]
    if panel_lists:
        all_genes = set().union(*panel_lists)
        invariant_genes = sorted(set.intersection(*panel_lists))
        variant_genes = sorted(all_genes - set(invariant_genes))
    else:
        invariant_genes = []
        variant_genes = []

    return {
        "cv_auc_range": (min(cv_aucs), max(cv_aucs)) if cv_aucs else (float("nan"), float("nan")),
        "gap_range": (min(gaps), max(gaps)) if gaps else (float("nan"), float("nan")),
        "phi_range": (min(phis), max(phis)) if phis else (float("nan"), float("nan")),
        "invariant_phi_label": invariant_phi_label,
        "invariant_genes": invariant_genes,
        "variant_genes": variant_genes,
        "invariant_banner": invariant_banner,
        "n_runs": n,
    }


# =============================================================================
# PAGE CONFIG
# =============================================================================

st.set_page_config(
    page_title="RAPTOR - Biomarker Discovery",
    page_icon="🧬",
    layout="wide",
)

init_session_state()
render_sidebar()

st.title("Biomarker Discovery — Module 10")
st.caption(
    "Multi-method feature selection, classification, panel optimization, and "
    "clinical-grade enhancement (signature score, ratios, direction patterns, "
    "clinical metrics)."
)

# -----------------------------------------------------------------------------
# Help section
# -----------------------------------------------------------------------------

with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Turn a differential expression signal into a validated,
    clinically actionable biomarker panel.

    ## What this module does

    Module 10 takes your count matrix plus group labels and runs a complete
    biomarker discovery pipeline:

    1. **Feature selection (10A)** — Multi-method ranking via DE filtering,
       LASSO / Elastic Net, RFE, Boruta, mRMR, SHAP
    2. **Classification (10B)** — Logistic regression, Random Forest, SVM,
       XGBoost with nested cross-validation
    3. **Panel optimization (10C)** — Finds the smallest panel with
       maximal AUC
    4. **Biological annotation (10E)** — MyGene.info, pathway enrichment,
       optional literature (PubMed) and PPI (STRING)

    ## Enhanced analysis (M10 suite)

    When you pick an **intent**, the page runs the full enhancement
    suite after discovery:

    - **Intent classification** — diagnostic, prognostic, predictive,
      monitoring, exploratory, translational
    - **Signature score** — weighted risk score with learned stratification
      cutoffs (applicable to new patients)
    - **Direction pattern** — per-gene UP/DOWN direction with fold change
      and confidence, for cross-cohort replication
    - **Clinical metrics** — Youden's optimal threshold, bootstrap 95% CI
      for AUC, PPV/NPV at realistic prevalence, decision curve analysis
    - **Ratio biomarkers** — Top scoring pairs (TSP) search for
      self-normalizing gene ratios

    ## Inputs

    - **Count matrix + metadata** — from Module 2 session state, or upload
    - **DE result** *(optional)* — from Module 7, pre-filters candidates
    - **Ensemble consensus genes** *(optional)* — from Module 9, stronger
      pre-filter

    ## Typical workflow

    1. Upload or load counts/metadata
    2. Pick group column and the two groups to compare
    3. Choose **intent** (this controls the enhancement suite)
    4. Tick feature-selection methods, pick panel size range
    5. Click **Run** — results populate below
    6. Export the panel CSV, signature JSON, and report
    """)

st.markdown("---")

# -----------------------------------------------------------------------------
# Availability check
# -----------------------------------------------------------------------------

if not BIOMARKER_AVAILABLE:
    st.error(f"""
    **Module 10 is not available.**

    Import error: `{biomarker_import_error}`

    Install missing dependencies with:
    ```
    pip install scikit-learn xgboost
    ```

    Optional (for extra feature-selection methods):
    ```
    pip install boruta mrmr-selection shap lifelines
    ```
    """)
    st.stop()


# =============================================================================
# SECTION 1: DATA SOURCE
# =============================================================================

st.markdown("## 1. Data Source")

source_col1, source_col2 = st.columns([1, 1])

with source_col1:
    count_source = st.radio(
        "**Count matrix**",
        ["Use Module 2 session data", "Upload CSV/TSV"],
        key="m10_count_source",
    )

with source_col2:
    meta_source = st.radio(
        "**Metadata**",
        ["Use Module 2 session data", "Upload CSV/TSV"],
        key="m10_meta_source",
    )

counts_df: Optional[pd.DataFrame] = None
meta_df: Optional[pd.DataFrame] = None

# --- counts ---
if count_source == "Use Module 2 session data":
    counts_df = st.session_state.get('m2_count_data')
    if counts_df is None:
        st.warning(
            "No count matrix in Module 2 session. "
            "Either run Module 2 (Quality) first or switch to upload."
        )
    else:
        st.success(
            f"Loaded counts from Module 2: **{counts_df.shape[0]} genes × "
            f"{counts_df.shape[1]} samples**"
        )
else:
    count_file = st.file_uploader(
        "Count matrix (genes × samples)",
        type=["csv", "tsv", "txt"],
        key="m10_count_upload",
    )
    if count_file is not None:
        try:
            sep = '\t' if count_file.name.endswith(('.tsv', '.txt')) else ','
            counts_df = pd.read_csv(count_file, sep=sep, index_col=0)
            st.success(
                f"Loaded: **{counts_df.shape[0]} genes × "
                f"{counts_df.shape[1]} samples**"
            )
        except Exception as e:
            st.error(f"Failed to read count matrix: {e}")

# --- metadata ---
if meta_source == "Use Module 2 session data":
    meta_df = st.session_state.get('m2_metadata')
    if meta_df is None:
        st.warning(
            "No metadata in Module 2 session. "
            "Either run Module 2 first or switch to upload."
        )
    else:
        st.success(
            f"Loaded metadata from Module 2: **{meta_df.shape[0]} samples × "
            f"{meta_df.shape[1]} columns**"
        )
else:
    meta_file = st.file_uploader(
        "Metadata (samples × columns, sample IDs in first column)",
        type=["csv", "tsv", "txt"],
        key="m10_meta_upload",
    )
    if meta_file is not None:
        try:
            sep = '\t' if meta_file.name.endswith(('.tsv', '.txt')) else ','
            meta_df = pd.read_csv(meta_file, sep=sep, index_col=0)
            st.success(
                f"Loaded: **{meta_df.shape[0]} samples × "
                f"{meta_df.shape[1]} columns**"
            )
        except Exception as e:
            st.error(f"Failed to read metadata: {e}")

# --- optional upstream inputs ---
st.markdown("### Optional upstream inputs")

has_de = st.session_state.get('m7_complete', False)
has_ensemble = st.session_state.get('m9_complete', False)

up_col1, up_col2 = st.columns(2)

with up_col1:
    use_de = st.checkbox(
        "Use DE result from Module 7",
        value=has_de,
        disabled=not has_de,
        help=(
            "Pre-filters biomarker candidates to DE-significant genes. "
            "Faster and usually stronger than all-genes search."
            if has_de else
            "Run Module 7 (Import DE) to enable."
        ),
        key="m10_use_de",
    )
    if has_de:
        de_dict = st.session_state.get('m7_de_results', {})
        de_method_choice = st.selectbox(
            "Which DE method to use?",
            list(de_dict.keys()) if de_dict else ["(none available)"],
            key="m10_de_method",
        )
    else:
        de_method_choice = None

with up_col2:
    use_ensemble = st.checkbox(
        "Use ensemble consensus from Module 9",
        value=has_ensemble,
        disabled=not has_ensemble,
        help=(
            "Pre-filters to multi-method consensus genes. Strongest prior — "
            "recommended when you have several DE runs."
            if has_ensemble else
            "Run Module 9 (Ensemble Analysis) to enable."
        ),
        key="m10_use_ensemble",
    )

# --- resolve upstream objects for later ---
de_result_obj = None
ensemble_result_obj = None
if use_de and has_de and de_method_choice and de_method_choice != "(none available)":
    de_result_obj = st.session_state['m7_de_results'].get(de_method_choice)
if use_ensemble and has_ensemble:
    m9_pkg = st.session_state.get('m9_result', {})
    ensemble_result_obj = m9_pkg.get('ensemble_result')


# =============================================================================
# SECTION 2: GROUP COLUMN + COMPARISON
# =============================================================================

st.markdown("---")
st.markdown("## 2. Group Setup")

group_column: Optional[str] = None
baseline_group: Optional[str] = None
reference_group: Optional[str] = None

if meta_df is not None:
    candidate_cols = [c for c in meta_df.columns if meta_df[c].nunique() >= 2]
    if not candidate_cols:
        st.error(
            "No metadata column has at least 2 distinct values. "
            "Cannot define groups."
        )
    else:
        # Prefer common names
        default_idx = 0
        for preferred in ("condition", "group", "treatment", "Condition", "Group"):
            if preferred in candidate_cols:
                default_idx = candidate_cols.index(preferred)
                break

        group_column = st.selectbox(
            "**Group column** in metadata",
            candidate_cols,
            index=default_idx,
            key="m10_group_col",
        )

        unique_groups = sorted(meta_df[group_column].dropna().unique().tolist())
        st.caption(
            f"Column `{group_column}` has {len(unique_groups)} distinct value(s): "
            + ", ".join(f"`{g}`" for g in unique_groups[:10])
            + ("…" if len(unique_groups) > 10 else "")
        )

        if len(unique_groups) < 2:
            st.error("Need at least 2 groups to compare.")
        else:
            g_col1, g_col2 = st.columns(2)
            with g_col1:
                baseline_group = st.selectbox(
                    "**Baseline group** (label = 0)",
                    unique_groups,
                    index=0,
                    help="Usually the control / healthy / untreated group.",
                    key="m10_baseline",
                )
            with g_col2:
                ref_options = [g for g in unique_groups if g != baseline_group]
                reference_group = st.selectbox(
                    "**Reference group** (label = 1)",
                    ref_options,
                    index=0,
                    help="Usually the disease / treated / condition of interest.",
                    key="m10_reference",
                )

            # Sample counts
            n_base = int((meta_df[group_column] == baseline_group).sum())
            n_ref = int((meta_df[group_column] == reference_group).sum())
            st.info(
                f"Samples in comparison: **{n_base}** {baseline_group} vs "
                f"**{n_ref}** {reference_group}"
            )
else:
    st.info("Load metadata above to configure groups.")


# =============================================================================
# SECTION 3: INTENT
# =============================================================================

st.markdown("---")
st.markdown("## 3. Biomarker Intent")
st.caption(
    "The *goal* of your analysis. This controls which enhancement modules "
    "run and what validation strategy is recommended."
)

intent_labels = {
    "none":          "— None (base pipeline only, no enhancement modules) —",
    "diagnostic":    "Diagnostic — Does the patient have the disease?",
    "prognostic":    "Prognostic — How will the patient fare over time?",
    "predictive":    "Predictive — How will the patient respond to treatment?",
    "monitoring":    "Monitoring — Is the disease progressing/regressing?",
    "exploratory":   "Exploratory — What patterns exist? (unsupervised)",
    "translational": "Translational — Will mouse findings replicate in humans?",
}

intent_options = ["none"] + list(VALID_INTENTS)

intent_key = st.radio(
    "Choose intent",
    intent_options,
    format_func=lambda k: intent_labels.get(k, k),
    horizontal=False,
    key="m10_intent",
)

source_species: Optional[str] = None
target_species: Optional[str] = None
custom_label: Optional[str] = None

if intent_key == "translational":
    ts_col1, ts_col2 = st.columns(2)
    with ts_col1:
        source_species = st.selectbox(
            "Source species",
            ["mouse", "rat", "human", "zebrafish", "other"],
            index=0,
            key="m10_source_species",
        )
    with ts_col2:
        target_species = st.selectbox(
            "Target species",
            ["human", "mouse", "rat", "other"],
            index=0,
            key="m10_target_species",
        )

with st.expander("Advanced intent options"):
    custom_label = st.text_input(
        "Custom output label (optional)",
        value="",
        help="Overrides the default label shown in reports and dashboards.",
        key="m10_custom_label",
    )
    if not custom_label:
        custom_label = None

# Preview
if intent_key == "none":
    st.info(
        "**Base pipeline only.** No enhancement modules will run. "
        "The classic biomarker discovery pipeline (feature selection, "
        "classification, panel optimization, annotation) will execute, "
        "producing a plain `BiomarkerResult`. Signature score, direction "
        "pattern, clinical metrics, and ratio biomarker tabs will remain "
        "empty. Choose a specific intent above to enable them."
    )
else:
    try:
        if intent_key == "translational":
            _preview = BiomarkerIntent(
                intent=intent_key,
                source_species=source_species,
                target_species=target_species,
                custom_label=custom_label,
            )
        else:
            _preview = BiomarkerIntent(
                intent=intent_key,
                custom_label=custom_label,
            )
        with st.expander("Intent configuration preview", expanded=False):
            st.code(_preview.describe(), language="text")
    except Exception as e:
        st.warning(f"Intent preview failed: {e}")


# =============================================================================
# SECTION 4: FEATURE SELECTION + PANEL CONFIG
# =============================================================================

st.markdown("---")
st.markdown("## 4. Feature Selection & Panel")

# Build available method list based on installed packages + upstream inputs
available_methods = ["elastic_net", "lasso", "rfe"]
if de_result_obj is not None or ensemble_result_obj is not None:
    available_methods.insert(0, "de_filter")
if _BORUTA_AVAILABLE:
    available_methods.append("boruta")
if _MRMR_AVAILABLE:
    available_methods.append("mrmr")
if _SHAP_AVAILABLE:
    available_methods.append("shap")

# Default: sensible small set
default_methods = [m for m in available_methods if m in ("de_filter", "elastic_net", "rfe")]

fs_col1, fs_col2 = st.columns([1.3, 1])
with fs_col1:
    methods_selected = st.multiselect(
        "**Feature selection methods** (consensus across methods)",
        options=available_methods,
        default=default_methods,
        help=(
            "Each method votes for candidate genes; final ranking is a "
            "multi-method consensus. More methods = slower but more robust."
        ),
        key="m10_methods",
    )
    # Show what's unavailable
    unavailable = []
    if not _BORUTA_AVAILABLE:
        unavailable.append("boruta (install: `pip install boruta`)")
    if not _MRMR_AVAILABLE:
        unavailable.append("mrmr (install: `pip install mrmr-selection`)")
    if not _SHAP_AVAILABLE:
        unavailable.append("shap (install: `pip install shap`)")
    if unavailable:
        st.caption("Unavailable: " + "; ".join(unavailable))

with fs_col2:
    min_panel = st.number_input(
        "Min panel size", min_value=1, max_value=100, value=5,
        help=(
            "Lowest panel size tested by the optimizer. Set to 1 to permit "
            "single-gene panels (useful for clinical-metrics demonstrations "
            "where imperfect sens/spec is required)."
        ),
        key="m10_min_panel",
    )
    max_panel = st.number_input(
        "Max panel size", min_value=min_panel, max_value=200, value=30,
        key="m10_max_panel",
    )
    target_panel_size = st.number_input(
        "Target panel size (0 = auto)", min_value=0, max_value=200, value=0,
        key="m10_target_panel",
    )

val_col1, val_col2 = st.columns(2)
with val_col1:
    validation = st.selectbox(
        "Validation strategy",
        ["nested_cv", "loocv"],
        index=0,
        help=(
            "nested_cv: nested k-fold cross-validation (recommended). "
            "loocv: leave-one-out (use for very small datasets, n < 30)."
        ),
        key="m10_validation",
    )
with val_col2:
    n_folds = st.number_input(
        "CV folds", min_value=3, max_value=10, value=5,
        disabled=(validation == "loocv"),
        key="m10_n_folds",
    )


# =============================================================================
# SECTION 5: ANNOTATION
# =============================================================================

st.markdown("---")
st.markdown("## 5. Biological Annotation")

ann_col1, ann_col2, ann_col3 = st.columns(3)
with ann_col1:
    annotate = st.checkbox(
        "Run annotation pipeline", value=True,
        help="MyGene.info + pathway enrichment",
        key="m10_annotate",
    )
with ann_col2:
    run_literature = st.checkbox(
        "PubMed literature search", value=False,
        disabled=not annotate,
        help="Slower — queries PubMed for each panel gene.",
        key="m10_run_lit",
    )
with ann_col3:
    run_ppi = st.checkbox(
        "STRING PPI network", value=False,
        disabled=not annotate,
        help="Fetches protein-protein interactions from STRING.",
        key="m10_run_ppi",
    )

sp_col1, sp_col2 = st.columns(2)
with sp_col1:
    species = st.selectbox(
        "Species (for annotation)",
        ["human", "mouse", "rat"],
        index=0,
        disabled=not annotate,
        key="m10_species",
    )
with sp_col2:
    disease_term = st.text_input(
        "Disease term (for literature search)",
        value="",
        placeholder="e.g. 'breast cancer', 'Alzheimer disease'",
        disabled=not (annotate and run_literature),
        key="m10_disease",
    )


# =============================================================================
# SECTION 6: ADVANCED
# =============================================================================

with st.expander("Advanced options", expanded=False):
    adv_col1, adv_col2 = st.columns(2)
    with adv_col1:
        prevalence = st.number_input(
            "Disease prevalence (for PPV/NPV)",
            min_value=0.001, max_value=0.999, value=0.05, step=0.01,
            format="%.3f",
            help=(
                "Realistic population prevalence — used to translate "
                "case-control sensitivity/specificity into clinically "
                "meaningful PPV and NPV."
            ),
            key="m10_prevalence",
        )
    with adv_col2:
        random_state = st.number_input(
            "Random seed", min_value=0, max_value=99999, value=42,
            key="m10_seed",
        )
    output_dir = st.text_input(
        "Output directory",
        value="results/biomarkers",
        key="m10_output_dir",
    )


# =============================================================================
# SECTION 7: RUN
# =============================================================================

st.markdown("---")
st.markdown("## 6. Run Discovery")

# Validate prerequisites
ready = (
    counts_df is not None
    and meta_df is not None
    and group_column is not None
    and baseline_group is not None
    and reference_group is not None
    and baseline_group != reference_group
    and len(methods_selected) >= 1
)

if not ready:
    missing = []
    if counts_df is None:
        missing.append("count matrix")
    if meta_df is None:
        missing.append("metadata")
    if group_column is None:
        missing.append("group column")
    if baseline_group == reference_group:
        missing.append("distinct baseline and reference groups")
    if len(methods_selected) < 1:
        missing.append("at least one feature-selection method")
    st.warning(f"**Not ready yet.** Missing: {', '.join(missing)}.")

run_clicked = st.button(
    "🚀 Run Biomarker Discovery",
    type="primary",
    disabled=not ready,
    use_container_width=True,
)

if run_clicked and ready:
    # Subset metadata to only the two selected groups — M10 assumes binary
    mask = meta_df[group_column].isin([baseline_group, reference_group])
    meta_subset = meta_df.loc[mask].copy()
    common_samples = [s for s in meta_subset.index if s in counts_df.columns]
    if len(common_samples) < len(meta_subset):
        st.warning(
            f"Metadata has {len(meta_subset)} samples in the selected groups, "
            f"but only {len(common_samples)} are in the count matrix columns. "
            f"Using the {len(common_samples)} intersecting samples."
        )
    meta_subset = meta_subset.loc[common_samples]
    counts_subset = counts_df[common_samples]

    # Build intent object
    try:
        if intent_key == "none":
            intent_obj = None
        elif intent_key == "translational":
            intent_obj = BiomarkerIntent(
                intent=intent_key,
                source_species=source_species,
                target_species=target_species,
                custom_label=custom_label,
            )
        else:
            intent_obj = BiomarkerIntent(
                intent=intent_key,
                custom_label=custom_label,
            )
    except Exception as e:
        st.error(f"Invalid intent configuration: {e}")
        st.stop()

    target = int(target_panel_size) if target_panel_size > 0 else None

    progress = st.progress(0, text="Preparing…")
    status = st.empty()

    try:
        status.info("Running Module 10 pipeline — this may take a few minutes.")
        progress.progress(10, text="Loading data and running feature selection…")

        result = discover_biomarkers(
            counts=counts_subset,
            metadata=meta_subset,
            group_column=group_column,
            # Pass the UI selection through so _prepare_expression_data and
            # the enhancement pipeline encode labels in the user's intended
            # polarity (baseline=0, reference=1). Without these, core.py
            # would fall back to alphabetical sort — flipping direction /
            # signature / ratio interpretation whenever the reference group
            # happens to sort before the baseline (e.g. "disease" < "healthy").
            baseline_group=baseline_group,
            reference_group=reference_group,
            de_result=de_result_obj,
            ensemble_result=ensemble_result_obj,
            methods=list(methods_selected),
            target_panel_size=target,
            min_panel=int(min_panel),
            max_panel=int(max_panel),
            validation=validation,
            n_folds=int(n_folds),
            species=species,
            disease_term=(disease_term or None),
            annotate=annotate,
            run_literature=run_literature,
            run_ppi=run_ppi,
            output_dir=output_dir,
            random_state=int(random_state),
            verbose=False,
            intent=intent_obj,
            prevalence=float(prevalence),
        )

        progress.progress(100, text="Done.")
        status.empty()

        # Store in session state (using m10_ convention)
        st.session_state['m10_result'] = result
        st.session_state['m10_complete'] = True
        st.session_state['m10_running'] = False
        st.session_state['m10_error'] = False

        # Save the derived inputs that the multi-seed verification
        # button needs to re-run discovery at additional seeds.
        # Within-session only; session_state is ephemeral by design.
        #
        # IMPORTANT: do NOT re-write widget-owned keys here. Keys like
        # m10_min_panel and m10_max_panel are bound to st.number_input
        # widgets earlier in the form; Streamlit raises
        # "cannot be modified after the widget ... is instantiated"
        # if we assign to them. Multi-seed verification reads those
        # widget-owned keys directly via session_state.get(...) at
        # button-click time, so they're already accessible.
        #
        # The values we DO save below are either (a) derived computed
        # state (counts/metadata subsets that aren't directly widget-
        # bound) or (b) plain values stored under non-widget keys
        # (suffix _resolved) so they don't collide with widget keys.
        st.session_state['biomarker_counts'] = counts_subset
        st.session_state['biomarker_metadata'] = meta_subset
        st.session_state['m10_group_column_resolved'] = group_column
        st.session_state['m10_reference_group_resolved'] = reference_group
        st.session_state['m10_baseline_group_resolved'] = baseline_group
        st.session_state['m10_target_panel_resolved'] = int(target) if target else int(max_panel)
        # Also clear any prior multi-seed cache so a fresh discovery
        # doesn't reuse stale results.
        st.session_state.pop('multi_seed_results', None)
        st.session_state.pop('multi_seed_run_requested', None)
        # NOTE: key is 'm10_intent_config', not 'm10_intent' — the latter is
        # already bound as the intent radio widget's key, so Streamlit raises
        # "cannot be modified after the widget with key ... is instantiated"
        # when we try to overwrite it here.
        st.session_state['m10_intent_config'] = (
            intent_obj.to_dict() if intent_obj is not None else None
        )
        st.session_state['m10_groups'] = {
            "baseline": baseline_group,
            "reference": reference_group,
            "group_column": group_column,
        }

        base_panel = (
            result.panel if hasattr(result, "panel") else result.base_result.panel
        )
        st.success(
            f"**Complete!** Final panel: **{len(base_panel)} genes** | "
            f"Best classifier: **{_best_clf_name(result)}** | "
            f"AUC: **{_best_clf_auc(result):.3f}**"
        )
        st.balloons()

    except Exception as e:
        progress.empty()
        status.empty()
        st.session_state['m10_error'] = True
        st.session_state['m10_running'] = False
        st.error(f"""
**Biomarker discovery failed.**

Error: `{e}`

**Troubleshooting:**
- Ensure the selected groups have at least 3 samples each
- Check that count values are non-negative integers (raw counts)
- Try fewer feature-selection methods (start with elastic_net + rfe)
- Reduce max_panel if the dataset is small
""")


# =============================================================================
# Multi-seed verification: progressive rendering helper (Option B2)
# =============================================================================
# Imported lazily inside the function so the dashboard doesn't pay the
# cost (and the import-error risk) unless the user actually clicks the
# verify button.

def _render_multi_seed_verification(
    base, result, verify_seeds: List[int],
) -> None:
    """Render the multi-seed verification table progressively.

    Each row appears as the corresponding seed completes. The first
    row reuses the existing discovery result (Cache option 1); the
    others call run_one_seed from seed_stability_investigation.py.
    Caches per-seed rows in st.session_state so navigating away and
    back preserves the results within the session.

    Note on `clin`: derived from `result` here rather than passed in,
    because this function is invoked from the Panel tab where `clin`
    isn't in scope (it's defined inside the Clinical Metrics tab
    block). EnhancedBiomarkerResult exposes clinical_metrics directly;
    bare BiomarkerResult does not.
    """
    import sys as _sys
    from pathlib import Path as _Path
    import warnings as _warnings

    # Derive clin from the result. EnhancedBiomarkerResult has it as
    # a dataclass field; a bare BiomarkerResult doesn't. Fall back to
    # an empty dict so downstream lookups don't crash.
    clin = getattr(result, "clinical_metrics", None) or {}

    # Resolve and import the script. The script lives at the repo
    # root; the dashboard is launched from there. Robust-fallback to
    # several relative locations.
    _script_paths = [
        _Path.cwd() / "seed_stability_investigation.py",
        _Path(__file__).resolve().parent.parent.parent.parent
            / "seed_stability_investigation.py",
    ]
    _script_path = None
    for _p in _script_paths:
        if _p.exists():
            _script_path = _p
            break

    if _script_path is None:
        st.error(
            "Multi-seed verification requires "
            "`seed_stability_investigation.py` at the RAPTOR repo "
            "root. Could not find it. Looked in: "
            + ", ".join(str(p) for p in _script_paths)
        )
        return

    # Import the script as a module
    if str(_script_path.parent) not in _sys.path:
        _sys.path.insert(0, str(_script_path.parent))
    try:
        import seed_stability_investigation as _ssi
    except Exception as e:
        st.error(f"Failed to import seed_stability_investigation: {e}")
        return

    # Resolve the input data and discovery parameters from
    # session_state. The user's existing discovery already ran with
    # specific values for these; we need to reuse them so the new
    # seed runs are comparable.
    #
    # Reads keys carefully:
    #   - biomarker_counts/metadata: set by the discovery success path
    #     (computed slices not directly bound to a widget)
    #   - m10_*_resolved: also set by the success path under
    #     non-widget names so they don't collide with widget keys
    #   - m10_min_panel / m10_max_panel: widget-owned keys, read
    #     directly via session_state.get; safe because we're reading,
    #     not writing
    counts_df = st.session_state.get("biomarker_counts")
    metadata_df = st.session_state.get("biomarker_metadata")
    group_col = st.session_state.get("m10_group_column_resolved", "condition")
    ref = st.session_state.get("m10_reference_group_resolved", "disease")
    baseline = st.session_state.get("m10_baseline_group_resolved", "healthy")
    min_panel = int(st.session_state.get("m10_min_panel", 3))
    max_panel = int(st.session_state.get("m10_max_panel", 5))
    target_panel = int(st.session_state.get(
        "m10_target_panel_resolved", max_panel,
    ))

    if counts_df is None or metadata_df is None:
        st.error(
            "Multi-seed verification needs the original counts and "
            "metadata in session state. Re-run discovery and try "
            "again."
        )
        return

    # Reuse cache from prior verify clicks within this session.
    cache_key = "multi_seed_results"
    if cache_key not in st.session_state:
        st.session_state[cache_key] = {}
    cached: dict = st.session_state[cache_key]

    # Set up the row builder for the current-seed reuse path.
    def _row_from_existing_result() -> dict:
        og = clin.get("optimism_gap", {}) if clin else {}
        bc = clin.get("bootstrap_ci", {}) if clin else {}
        stab = getattr(base, "panel_stability", None)
        cv_auc = og.get("cv_auc_oof", float("nan"))
        gap = og.get("gap", float("nan"))
        return {
            "seed": verify_seeds[0],
            "best_clf": getattr(base, "best_classifier", "n/a"),
            "cv_auc": cv_auc,
            "train_auc": og.get("training_auc", float("nan")),
            "gap": gap,
            "banner": _ssi.banner_category(cv_auc, gap),
            "phi": (
                stab.nogueira_stability if stab else float("nan")
            ),
            "phi_label": (
                stab.benchmark_label if stab else "n/a"
            ),
            "boot_lo": bc.get("ci_lower", float("nan")),
            "boot_hi": bc.get("ci_upper", float("nan")),
            "panel": ",".join(getattr(base, "panel", [])),
        }

    # Progressive rendering: a single placeholder we re-render after
    # each seed completes. Streamlit's `st.empty()` is the standard
    # idiom for this in non-fragment code.
    progress_area = st.empty()
    rows: list = []

    # Process seeds in order: current first (reuse), then new ones.
    outdir_root = _Path("/tmp/raptor_dashboard_multiseed")
    outdir_root.mkdir(parents=True, exist_ok=True)

    for i, seed in enumerate(verify_seeds):
        # Cache hit?
        if seed in cached:
            rows.append(cached[seed])
            with progress_area.container():
                _draw_multiseed_table(rows, total=len(verify_seeds))
            continue

        # First seed = reuse existing discovery
        if i == 0:
            row = _row_from_existing_result()
            cached[seed] = row
            rows.append(row)
            with progress_area.container():
                _draw_multiseed_table(rows, total=len(verify_seeds))
            continue

        # New seed: run discovery
        with progress_area.container():
            _draw_multiseed_table(
                rows, total=len(verify_seeds),
                pending_seed=seed,
            )

        try:
            with _warnings.catch_warnings():
                _warnings.simplefilter("ignore")
                row = _ssi.run_one_seed(
                    counts=counts_df, metadata=metadata_df,
                    seed=seed, group_col=group_col,
                    ref=ref, base=baseline,
                    outdir_root=outdir_root,
                    min_panel=min_panel,
                    max_panel=max_panel,
                    target=target_panel,
                )
        except Exception as e:
            st.error(f"Discovery failed at seed={seed}: {e}")
            continue

        cached[seed] = row
        rows.append(row)
        with progress_area.container():
            _draw_multiseed_table(rows, total=len(verify_seeds))

    # All seeds done -- render the invariants/variants summary
    if len(rows) >= 2:
        with progress_area.container():
            _draw_multiseed_table(rows, total=len(verify_seeds))
            _draw_invariants_summary(rows)


def _draw_multiseed_table(
    rows: list, total: int, pending_seed: Optional[int] = None,
) -> None:
    """Render the per-seed comparison table. Helper for progressive UI."""
    if not rows:
        st.info(f"Running 0/{total} seeds ...")
        return

    df = pd.DataFrame(rows)
    display_df = df[[
        "seed", "cv_auc", "train_auc", "gap", "banner",
        "phi", "phi_label", "boot_lo", "boot_hi", "best_clf",
    ]].copy()
    display_df.columns = [
        "Seed", "CV AUC", "Train AUC", "Gap", "Banner",
        "Φ", "Φ label", "Boot CI lo", "Boot CI hi", "Best clf",
    ]

    # Float formatting
    for col in ["CV AUC", "Train AUC", "Gap", "Φ", "Boot CI lo", "Boot CI hi"]:
        display_df[col] = display_df[col].map(
            lambda v: f"{v:.3f}" if pd.notna(v) else "n/a"
        )

    if pending_seed is not None:
        st.info(
            f"Running {len(rows)}/{total} seeds. "
            f"In progress: seed={pending_seed} ..."
        )
    else:
        st.success(f"Completed {len(rows)}/{total} seeds.")

    st.dataframe(display_df, use_container_width=True, hide_index=True)


def _draw_invariants_summary(rows: list) -> None:
    """Render the invariant/variant interpretation summary."""
    summary = _summarize_seed_runs(rows)

    st.markdown("#### What's stable across seeds (trust) vs what varies (caution)")

    # Two-column layout: invariants left, variants right
    inv_col, var_col = st.columns(2)

    with inv_col:
        st.markdown("**Invariants — trust these**")
        if summary["invariant_phi_label"]:
            st.markdown(
                f"- Φ category: **{summary['invariant_phi_label']}** "
                f"(across all {summary['n_runs']} runs)"
            )
        if summary["invariant_banner"]:
            st.markdown(
                f"- Banner: **{summary['invariant_banner']}** "
                f"(across all {summary['n_runs']} runs)"
            )
        if summary["invariant_genes"]:
            st.markdown(
                f"- Genes selected in **all {summary['n_runs']} runs**: "
                + ", ".join(f"`{g}`" for g in summary["invariant_genes"])
            )
        else:
            st.markdown(
                "- *No genes appeared in every run's panel — "
                "panel composition is unstable.*"
            )

    with var_col:
        st.markdown("**Variants — treat with caution**")
        cv_lo, cv_hi = summary["cv_auc_range"]
        gap_lo, gap_hi = summary["gap_range"]
        phi_lo, phi_hi = summary["phi_range"]
        if pd.notna(cv_lo):
            st.markdown(
                f"- CV AUC range: **{cv_lo:.3f}–{cv_hi:.3f}** "
                f"(spread {cv_hi - cv_lo:.3f})"
            )
        if pd.notna(gap_lo):
            st.markdown(
                f"- Gap range: **{gap_lo:+.3f} to {gap_hi:+.3f}**"
            )
        if pd.notna(phi_lo):
            st.markdown(
                f"- Φ range: **{phi_lo:.3f}–{phi_hi:.3f}**"
            )
        if not summary["invariant_banner"]:
            st.markdown(
                "- Banner category **varies across seeds** "
                "(see table)"
            )
        if summary["variant_genes"]:
            st.markdown(
                f"- Genes appearing in some but not all runs "
                f"({len(summary['variant_genes'])}): "
                + ", ".join(
                    f"`{g}`" for g in summary["variant_genes"][:8]
                )
                + (" ..." if len(summary["variant_genes"]) > 8 else "")
            )

    # Pedagogical caption -- the takeaway in one sentence.
    if summary["invariant_genes"] and not summary["invariant_banner"]:
        st.info(
            f"💡 **Trust the invariants, doubt the variants.** "
            f"The {len(summary['invariant_genes'])} genes selected "
            f"in every run reflect the real biology recovered by "
            f"this pipeline. The CV AUC variance reflects fold-split "
            f"sensitivity at this cohort size — single-seed point "
            f"estimates do not reliably capture clinical "
            f"performance. Report panel composition as the primary "
            f"discovery output."
        )


# =============================================================================
# SECTION 8: RESULTS
# =============================================================================

if st.session_state.get('m10_complete', False):
    st.markdown("---")
    st.markdown("## 7. Results")

    result = st.session_state['m10_result']
    base = _base(result)
    is_enhanced = hasattr(result, "base_result")

    # ---- Overview tiles ----
    _panel_stab = getattr(base, "panel_stability", None)
    if _panel_stab is not None and len(_panel_stab.per_fold_panels) > 0:
        # Five-column layout: add Phi as the final header metric with
        # a traffic-light color pulled from _stability_color. Using
        # raw markdown instead of st.metric so we can color the value
        # directly — st.metric's color knob only accepts 'normal' /
        # 'inverse' / 'off' which is not enough for our 3-state scheme.
        o1, o2, o3, o4, o5 = st.columns(5)
        o1.metric("Panel size", f"{base.panel_size} genes")
        o2.metric("Best classifier", _best_clf_name(result))
        o3.metric("Best AUC", f"{_best_clf_auc(result):.3f}")
        o4.metric("Samples", f"{base.n_samples}")
        with o5:
            st.caption(
                "Panel stability (Φ)",
                help=(
                    "Panel stability (Nogueira 2018) measures cross-fold "
                    "agreement of per-fold feature-selection panels. "
                    "Values range 0–1. Interpret alongside CV AUC: "
                    "high AUC with low Φ indicates small-cohort "
                    "variability; low AUC with low Φ indicates absent "
                    "signal."
                ),
            )
            st.markdown(
                _render_stability_badge_markdown(
                    _panel_stab.nogueira_stability,
                    _panel_stab.nogueira_stability_ci,
                    _panel_stab.benchmark_label,
                ),
                unsafe_allow_html=True,
            )
    else:
        # Pre-M6 pickled result or LOOCV path without stability data:
        # fall back to the original four-column layout untouched.
        o1, o2, o3, o4 = st.columns(4)
        o1.metric("Panel size", f"{base.panel_size} genes")
        o2.metric("Best classifier", _best_clf_name(result))
        o3.metric("Best AUC", f"{_best_clf_auc(result):.3f}")
        o4.metric("Samples", f"{base.n_samples}")

    # ---- LOOCV auto-override banner (removed April 2026) ----
    # Historical note: pre-M6 core.py had an n < 20 auto-override that
    # silently applied LOOCV regardless of the UI's validation
    # selection. This banner existed to explain that auto-switch and
    # flag that trained_model wouldn't be preserved. M6 removed the
    # auto-override; nested_cv now runs for any n, including very
    # small cohorts. The banner's trigger condition (base.n_samples
    # < 20) therefore has no legitimate path in the current pipeline.
    # In practice it was firing on every small-cohort nested_cv run
    # (e.g. the 10a paper-figure scenario, n=10) and contradicting
    # the populated clinical metrics, panel stability, and optimism
    # gap displayed immediately below it. Removed entirely rather
    # than repurposed because validation='loocv' is deprecated as of
    # M6 (DeprecationWarning in core.py); users who still choose it
    # explicitly don't need a nag banner. A regression test in
    # test_biomarker_optimism_diagnostic_m2.py asserts the banner
    # strings do not re-appear in this source file.

    # ---- Text summary ----
    with st.expander("Full text summary", expanded=False):
        st.code(
            result.summary() if is_enhanced else base.summary(),
            language="text",
        )

    # ---- Tabs ----
    tab_labels = [
        "Panel",
        "Ranked Genes",
        "Classification",
        "Panel Curve",
        "Signature Score",
        "Direction Pattern",
        "Clinical Metrics",
        "Ratio Biomarkers",
        "Annotations",
        "Downloads",
    ]
    tabs = st.tabs(tab_labels)

    # ---------- Tab 1: Panel ----------
    with tabs[0]:
        st.markdown(f"### Recommended panel ({base.panel_size} genes)")
        panel_df = pd.DataFrame({
            "rank": range(1, len(base.panel) + 1),
            "gene": base.panel,
        })
        # Enrich with consensus score if available
        if "consensus_score" in base.ranked_genes.columns:
            panel_df["consensus_score"] = panel_df["gene"].map(
                base.ranked_genes["consensus_score"]
            )
        if "n_methods_selected" in base.ranked_genes.columns:
            panel_df["n_methods"] = panel_df["gene"].map(
                base.ranked_genes["n_methods_selected"]
            )
        st.dataframe(panel_df, use_container_width=True, hide_index=True)

        # ---- Panel Stability subsection (M6 + polish batch) ----
        _stab = getattr(base, "panel_stability", None)
        if _stab is not None and len(_stab.per_fold_panels) > 0:
            st.markdown("---")
            st.markdown("### Panel Stability")
            st.caption(
                "How consistent is this panel across cross-validation folds? "
                "(Nogueira et al. 2018, JMLR 18:174)"
            )

            color_map = {
                "green": "#2E7D32",
                "amber": "#F57C00",
                "red": "#C62828",
            }
            _color = color_map[_stability_color(_stab.nogueira_stability)]
            _lo, _hi = _stab.nogueira_stability_ci

            # Large Phi display + CI + category label on the left, a block
            # of interpretation prose on the right.
            ps_left, ps_right = st.columns([1, 2])
            with ps_left:
                st.markdown(
                    f"<div style='color:{_color}; font-size:3em; "
                    f"font-weight:700; line-height:1'>Φ = "
                    f"{_stab.nogueira_stability:.3f}</div>"
                    f"<div style='color:{_color}; font-size:1em; "
                    f"margin-top:0.3em'>"
                    f"95% CI [{_lo:.3f}, {_hi:.3f}]</div>"
                    f"<div style='color:{_color}; font-size:1.2em; "
                    f"font-weight:600; margin-top:0.5em'>"
                    f"{_stab.benchmark_label.capitalize()}</div>",
                    unsafe_allow_html=True,
                )
            with ps_right:
                st.markdown(
                    "**Interpretation.** "
                    + _stability_interpretation(_stab.benchmark_label)
                )
                st.caption(
                    f"Computed across {_stab.n_folds} folds × "
                    f"{_stab.n_repeats} repeats = "
                    f"{len(_stab.per_fold_panels)} per-fold panels. "
                    f"The Nogueira measure is the only similarity-based "
                    f"stability measure that satisfies correction-for-"
                    f"chance (Jaccard, Dice, POG all inflate with panel "
                    f"size)."
                )

            # Per-gene cross-fold selection frequency
            freq_series = _stab.gene_selection_frequency
            if freq_series is not None and len(freq_series) > 0:
                st.markdown("#### Per-gene cross-fold selection frequency")
                st.caption(
                    "Fraction of folds in which each gene appeared in the "
                    "per-fold panel. Genes in the final recommended panel "
                    "(above) are highlighted."
                )
                freq_df = (
                    freq_series.sort_values(ascending=False)
                    .rename("selection_frequency")
                    .to_frame()
                )
                freq_df.index.name = "gene"
                freq_df["in_final_panel"] = freq_df.index.isin(base.panel)

                # Cap the displayed rows at 15, but offer a show-all toggle
                n_total = len(freq_df)
                show_all = st.toggle(
                    f"Show all {n_total} genes "
                    f"(default: top 15)",
                    value=False,
                    key="show_all_stability_genes",
                )
                _display_df = freq_df if show_all else freq_df.head(15)

                # Highlight final-panel genes with a soft green background
                def _highlight_panel(row):
                    if row["in_final_panel"]:
                        return ["background-color: #E8F5E9"] * len(row)
                    return [""] * len(row)

                st.dataframe(
                    _display_df.style.apply(_highlight_panel, axis=1).format(
                        {"selection_frequency": "{:.2f}"}
                    ),
                    use_container_width=True,
                )

                # CSV download for the full frequency table
                st.download_button(
                    "📥 panel_stability_frequency.csv",
                    data=freq_df.reset_index().to_csv(index=False).encode("utf-8"),
                    file_name="panel_stability_frequency.csv",
                    mime="text/csv",
                    key="dl_stability_freq",
                )

        # ---- Multi-seed verification ----
        # The "Verify across N seeds" feature: runs discover_biomarkers
        # at additional seeds and shows a comparison table that
        # surfaces invariants (trust these) vs variants (treat with
        # caution). Implementation reuses run_one_seed from the
        # standalone seed-stability investigation script
        # (seed_stability_investigation.py at the repo root) so the
        # dashboard and the paper-figure tooling stay in sync.
        st.markdown("---")
        st.markdown("### Multi-seed verification")
        st.caption(
            "Re-run discovery at additional random seeds to "
            "characterize run-to-run variance. The comparison table "
            "below the button shows which diagnostics are stable "
            "across seeds (trust these) vs which vary (treat with "
            "caution). Especially recommended at small n; see the "
            "Clinical Metrics tab for cohort-size advisory."
        )

        # Resolve the existing discovery seed so we can include it as
        # one of the verification runs without re-computing (Cache
        # option 1). Best-effort: try multiple session-state keys.
        _existing_seed = None
        for _key in ("m10_seed", "biomarker_seed", "random_state"):
            _candidate = st.session_state.get(_key)
            if _candidate is not None:
                try:
                    _existing_seed = int(_candidate)
                    break
                except (TypeError, ValueError):
                    pass
        if _existing_seed is None:
            _existing_seed = 42  # match the dashboard widget default

        # Seeds to run: existing + DEFAULT_VERIFICATION_SEEDS, dedupe
        # while preserving order so the existing seed appears first.
        _verify_seeds: List[int] = [_existing_seed]
        for _s in DEFAULT_VERIFICATION_SEEDS:
            if _s != _existing_seed:
                _verify_seeds.append(_s)
        _verify_seeds = _verify_seeds[:4]  # cap at 4 (existing + 3 new)

        st.caption(
            f"Will run at seeds: {_verify_seeds[0]} (current), "
            + ", ".join(str(s) for s in _verify_seeds[1:])
            + ". The current seed reuses your existing result; the "
            "others trigger fresh discovery runs (~10-60s each "
            "depending on cohort size and methods)."
        )

        if st.button(
            "🔄 Verify across multiple seeds",
            key="btn_multi_seed_verify",
            help="Re-runs feature selection, panel optimization, and "
                 "classification at additional seeds. Within-session "
                 "results are cached.",
        ):
            st.session_state["multi_seed_run_requested"] = True

        if st.session_state.get("multi_seed_run_requested"):
            # Render the comparison table progressively (Option B2).
            # Each seed's row appears as it completes. Existing seed
            # reuses the current displayed result; new seeds trigger
            # discovery via run_one_seed from the script.
            #
            # Note: clin is intentionally NOT passed here. clin is
            # only defined inside the Clinical Metrics tab scope; the
            # render function derives it from `result` itself, which
            # lets this button work from the Panel tab without a
            # cross-tab variable reference.
            _render_multi_seed_verification(
                base=base, result=result,
                verify_seeds=_verify_seeds,
            )

    # ---------- Tab 2: Ranked Genes ----------
    with tabs[1]:
        st.markdown(f"### All ranked candidates ({len(base.ranked_genes)} genes)")

        # ---- M4 chance-expectation header ----
        # When the ranked_genes table has the M4 calibration columns,
        # surface the observed-vs-expected comparator above the table.
        # Under H0 (uniform p-values), about alpha * n_genes genes will
        # pass by chance; showing observed-vs-expected makes the table
        # self-diagnostic without the user having to do any math.
        if (
            "p_value" in base.ranked_genes.columns
            and "weight" in base.ranked_genes.columns
        ):
            _alpha = base.parameters.get("alpha", 0.05)
            _n_genes = len(base.ranked_genes)
            _n_sig = int((base.ranked_genes["p_value"] < _alpha).sum())
            _n_nonsig = _n_genes - _n_sig
            _n_expected = _m4_chance_expectation(_alpha, _n_genes)
            st.markdown(
                f"**Significance calibration.** "
                f"**{_n_sig}** genes pass α={_alpha} "
                f"[chance expectation: ~{_n_expected} under the null "
                f"hypothesis]. Scores for the remaining "
                f"**{_n_nonsig}** genes are down-weighted by 0.5×."
            )

        st.caption(
            "Multi-method consensus ranking. `n_methods_selected` counts how "
            "many feature-selection methods chose this gene."
        )
        st.dataframe(
            base.ranked_genes.head(200),
            use_container_width=True,
            height=500,
        )
        if len(base.ranked_genes) > 200:
            st.caption(f"Showing top 200 of {len(base.ranked_genes)}. Download for full table.")

    # ---------- Tab 3: Classification ----------
    with tabs[2]:
        st.markdown("### Classification performance")

        # Cross-reference the Panel Stability diagnostic when present.
        # An inline one-liner is much more reliable than a Streamlit
        # anchor link; users who want more detail click the Panel tab.
        _stab_clf = getattr(base, "panel_stability", None)
        if _stab_clf is not None and len(_stab_clf.per_fold_panels) > 0:
            _phi_clf = _stab_clf.nogueira_stability
            _lbl_clf = _stab_clf.benchmark_label
            st.caption(
                f"Related diagnostic — Panel stability: "
                f"Φ={_phi_clf:.3f} ({_lbl_clf}). See Panel tab for full "
                f"interpretation and per-gene frequencies."
            )

        if not base.classification_results:
            st.info("No classification results available.")
        else:
            perf_rows = []
            for name, res in base.classification_results.items():
                # NOTE: ClassificationResult has sensitivity/specificity, not
                # precision/recall — the latter don't exist on the dataclass.
                # Sensitivity/specificity are also more standard in clinical
                # biomarker reporting.
                perf_rows.append({
                    "classifier": name,
                    "AUC": getattr(res, "auc", float("nan")),
                    "F1": getattr(res, "f1", float("nan")),
                    "accuracy": getattr(res, "accuracy", float("nan")),
                    "sensitivity": getattr(res, "sensitivity", float("nan")),
                    "specificity": getattr(res, "specificity", float("nan")),
                    "best": (name == base.best_classifier),
                })
            perf_df = pd.DataFrame(perf_rows).sort_values("AUC", ascending=False)
            st.dataframe(
                perf_df.style.format({
                    "AUC": "{:.3f}", "F1": "{:.3f}", "accuracy": "{:.3f}",
                    "sensitivity": "{:.3f}", "specificity": "{:.3f}",
                }),
                use_container_width=True,
                hide_index=True,
            )

            # Detect classifier collapse pattern. Common in LOOCV on small
            # balanced cohorts: default-hyperparameter tree ensembles
            # collapse to majority-class prediction, which produces AUC~0
            # under LOOCV (every held-out sample is predicted as the
            # opposite class since the training set's majority is always
            # the complement). Signal: at least one classifier has AUC<0.3
            # while at least one other reaches >0.7.
            _collapsed = perf_df[perf_df["AUC"] < 0.3]
            _succeeded = perf_df[perf_df["AUC"] > 0.7]
            if not _collapsed.empty and not _succeeded.empty:
                _names = ", ".join(_collapsed["classifier"].tolist())
                st.warning(
                    f"**Classifier collapse detected:** `{_names}` scored "
                    f"AUC<0.3 while others reached AUC>0.7. This is the "
                    f"classic LOOCV + balanced-classes + default-"
                    f"hyperparameter pathology — the classifier collapses "
                    f"to majority-class prediction, which is always wrong "
                    f"for the held-out sample under balanced LOOCV. This "
                    f"is a default-hyperparameter artifact on tiny data, "
                    f"not evidence that the gene panel is flawed. The "
                    f"successful classifiers above are the reliable "
                    f"evidence for this panel."
                )

            # Downloadable classification metrics CSV
            perf_csv = perf_df.drop(columns=["best"]).to_csv(
                index=False
            ).encode("utf-8")
            st.download_button(
                "📥 classification_performance.csv",
                data=perf_csv,
                file_name="classification_performance.csv",
                mime="text/csv",
                key="m10_dl_classification_perf_inline",
            )

            # Bar chart of AUCs
            fig = px.bar(
                perf_df,
                x="classifier",
                y="AUC",
                color="best",
                color_discrete_map={True: "#2E86DE", False: "#A4B0BE"},
                title="Classifier AUC comparison",
            )
            fig.update_layout(showlegend=False, yaxis_range=[0, 1])
            st.plotly_chart(fig, use_container_width=True)

    # ---------- Tab 4: Panel Curve ----------
    with tabs[3]:
        st.markdown("### Panel size vs. performance")
        if base.panel_optimization is None:
            st.info("Panel optimization was not run.")
        else:
            po = base.panel_optimization
            # Try common field names
            curve = None
            for attr in ("curve", "panel_curve", "results_df"):
                if hasattr(po, attr):
                    candidate = getattr(po, attr)
                    if isinstance(candidate, pd.DataFrame) and not candidate.empty:
                        curve = candidate
                        break
            if curve is not None:
                x_col = "panel_size" if "panel_size" in curve.columns else curve.columns[0]
                y_col = "auc" if "auc" in curve.columns else curve.columns[1]
                fig = px.line(
                    curve, x=x_col, y=y_col, markers=True,
                    title="Panel size vs. AUC",
                )
                # Clamp y-axis to [0, 1] — plotly's autoscaler nudges above 1.0
                # when all values are tied at the maximum, which looks wrong.
                fig.update_yaxes(range=[0, 1.05])
                optimal = getattr(po, "optimal_size", base.panel_size)

                # Branch the annotation label on selection_method (v2.3.0+).
                # The same vertical line conveys different meaning depending
                # on how the size was chosen: a knee, a fallback to argmax,
                # a user-specified target, or the legacy first-drop heuristic.
                # Old PanelOptimizationResult pickles without selection_method
                # default to a generic "Optimal" label.
                selection_method = getattr(po, "selection_method", None)
                if selection_method == "kneedle":
                    annotation_text = f"Knee = {optimal} (kneedle)"
                elif selection_method == "argmax_fallback":
                    annotation_text = f"Knee = {optimal} (argmax — no curvature)"
                elif selection_method == "argmax":
                    annotation_text = f"Optimal = {optimal} (argmax)"
                elif selection_method == "user_specified":
                    annotation_text = f"Target = {optimal} (user-specified)"
                elif selection_method == "consensus_pinned":
                    annotation_text = f"Target = {optimal} (consensus-pinned)"
                elif selection_method == "first_drop":
                    annotation_text = f"Optimal = {optimal} (first-drop, legacy)"
                else:
                    annotation_text = f"Optimal = {optimal}"

                fig.add_vline(
                    x=optimal, line_dash="dash", line_color="red",
                    annotation_text=annotation_text,
                )
                st.plotly_chart(fig, use_container_width=True)

                # When kneedle returned no knee and we fell back to argmax,
                # surface a brief explanation so users understand the line
                # marks the maximum AUC, not a curvature point.
                if selection_method == "argmax_fallback":
                    st.caption(
                        "Kneedle detected no knee in this curve "
                        "(saturated or near-flat). Panel size was chosen "
                        "as the smallest size at the maximum CV AUC."
                    )
            else:
                opt_size = getattr(po, "optimal_size", base.panel_size)
                opt_auc = getattr(po, "optimal_auc", _best_clf_auc(result))
                st.metric("Optimal size", opt_size)
                st.metric("Optimal AUC", f"{opt_auc:.3f}")

    # ---------- Tab 5: Signature Score ----------
    with tabs[4]:
        st.markdown("### Signature score")
        sig = getattr(result, "signature", None) if is_enhanced else None
        if sig is None:
            st.info(
                "No signature score available. "
                "Re-run with intent='diagnostic' or 'prognostic' to enable."
            )
        else:
            sc1, sc2, sc3 = st.columns(3)
            sc1.metric("Mode", sig.mode)
            sc2.metric("Genes in signature", len(sig.panel_genes))
            sc3.metric(
                "Risk groups",
                " / ".join(sig.risk_labels),
            )

            if sig.performance:
                st.markdown("**Performance at fit**")
                perf_df = pd.DataFrame([sig.performance]).T.reset_index()
                perf_df.columns = ["metric", "value"]
                st.dataframe(perf_df, use_container_width=True, hide_index=True)

            # Weights
            st.markdown("**Gene weights**")
            weights_df = pd.DataFrame({
                "gene": list(sig.weights.keys()),
                "weight": list(sig.weights.values()),
            }).sort_values("weight", key=abs, ascending=False)
            fig = px.bar(
                weights_df, x="weight", y="gene", orientation="h",
                color="weight", color_continuous_scale="RdBu_r",
                title="Signature coefficients (z-scored features)",
                color_continuous_midpoint=0,
            )
            fig.update_layout(height=max(300, 22 * len(weights_df)))
            st.plotly_chart(fig, use_container_width=True)

            # Cutoffs
            st.markdown("**Stratification cutoffs**")
            st.json(sig.cutoffs)

    # ---------- Tab 6: Direction Pattern ----------
    with tabs[5]:
        st.markdown("### Direction pattern (UP / DOWN)")
        dp = getattr(result, "direction_pattern", None) if is_enhanced else None
        if dp is None:
            st.info("No direction pattern available.")
        else:
            d1, d2, d3 = st.columns(3)
            d1.metric("UP in reference", dp.n_up)
            d2.metric("DOWN in reference", dp.n_down)
            d3.metric("Total genes", dp.n_genes)

            # Build a dataframe for display
            try:
                dir_df = dp.to_dataframe()
                # to_dataframe() uses 'log2FC' and 'neg_log10_p' as column
                # names and puts gene symbols in the index. Normalize to a
                # predictable schema so the volcano-plot check below matches.
                if dir_df.index.name is None or dir_df.index.name != "gene":
                    dir_df = dir_df.reset_index().rename(
                        columns={dir_df.index.name or "index": "gene"}
                    )
                else:
                    dir_df = dir_df.reset_index()
                dir_df = dir_df.rename(columns={
                    "log2FC": "log2fc",
                    "neg_log10_p": "neg_log10p",
                })
            except Exception:
                # Fallback: reconstruct from dicts with matching column names
                dir_df = pd.DataFrame({
                    "gene": list(dp.gene_directions.keys()),
                    "direction": list(dp.gene_directions.values()),
                    "log2fc": [dp.fold_changes.get(g, np.nan) for g in dp.gene_directions],
                    "neg_log10p": [dp.confidence.get(g, np.nan) for g in dp.gene_directions],
                })

            st.dataframe(dir_df, use_container_width=True, height=400)

            # Volcano-style visualization
            if {"log2fc", "neg_log10p", "direction"}.issubset(dir_df.columns):
                fig = px.scatter(
                    dir_df, x="log2fc", y="neg_log10p",
                    color="direction",
                    color_discrete_map={"UP": "#E74C3C", "DOWN": "#3498DB"},
                    hover_data=["gene"] if "gene" in dir_df.columns else None,
                    title=f"Direction pattern: {dp.reference_group} vs {dp.baseline_group}",
                )
                fig.add_vline(x=0, line_dash="dash", line_color="gray")
                fig.update_xaxes(title="log2 fold change (reference vs baseline)")
                fig.update_yaxes(title="−log10(p)")
                st.plotly_chart(fig, use_container_width=True)

    # ---------- Tab 7: Clinical Metrics ----------
    with tabs[6]:
        st.markdown("### Clinical metrics")

        # Small-n advisory: tiered by cohort size, surfaces run-to-run
        # variance scaling before the user reads the optimism banner
        # and clinical metrics. Helps users at small n avoid
        # over-interpreting a single-seed reading. See _small_n_advisory
        # docstring and the seed-stability investigation (April 2026).
        _adv = _small_n_advisory(base.n_samples)
        if _adv is not None:
            _level, _msg = _adv
            if _level == "strong":
                st.warning(_msg)
            else:
                st.info(_msg)
        clin = getattr(result, "clinical_metrics", None) if is_enhanced else None
        if not clin:
            st.info(
                "No clinical metrics available. "
                "Requires a diagnostic / prognostic / predictive intent "
                "and a classifier with predict_proba."
            )
        else:
            # ---- Optimism gap diagnostic ----
            # Triangulate the honest cross-validated AUC against the
            # training-data AUC. A large gap (>0.05) indicates that the
            # clinical metrics in the rest of this tab (Youden's J, PPV/NPV,
            # DCA, bootstrap CI) may be optimistic because they are all
            # computed from the final model's training-data predictions.
            # This is the user-facing version of Harrell's optimism
            # (apparent − cross-validated performance). 10c revealed that a
            # gap >0.10 on noise-only data silently endorses spurious
            # biomarkers as clinically valuable; exposing the gap at the
            # top of the tab lets users see the risk before reading metrics.
            #
            # M2 (April 2026): the extraction is now a single call to
            # _resolve_optimism_tiles. Under the M1+ path this reads
            # from enhanced.py's purpose-built optimism_gap struct (OOF
            # CV AUC vs true full-data training AUC). Older pickled
            # results fall through to a bootstrap-point-estimate
            # approximation with an explicit user caption.
            _cv_auc, _train_auc, _gap, _source = _resolve_optimism_tiles(
                clin, base,
            )

            if _cv_auc is not None and _train_auc is not None and _gap is not None:
                st.markdown("**Optimism diagnostic** — CV vs training AUC")
                g1, g2, g3 = st.columns(3)
                g1.metric("Cross-validated AUC", f"{_cv_auc:.3f}")
                g2.metric("Training-data AUC", f"{_train_auc:.3f}")
                g3.metric("Optimism gap", f"{_gap:+.3f}")

                # M2: fallback-path caption. When the result predates the
                # optimism-gap instrumentation, the tiles above are an
                # approximation: both values are OOF-derived under M1 so
                # the "gap" is statistical noise between two averaging
                # strategies, not real optimism. Transparency is the
                # point of the whole diagnostic; silent fallback would
                # invite cross-result confusion.
                if _source == "fallback":
                    st.caption(
                        "This result was generated before the optimism-gap "
                        "instrumentation was added. The Training-data AUC "
                        "tile shows the bootstrap CI point estimate as a "
                        "proxy for training AUC. For the full diagnostic "
                        "with a true full-data fit, re-run discovery."
                    )

                # Six-branch decision tree on (cv_auc, gap). Unlike the
                # pre-fix banner (which only checked |gap| and could
                # silently endorse chance-level AUC as "genuine
                # signal"), this version reads BOTH the gap magnitude
                # AND the CV AUC level so chance-level results never
                # land in the green banner.
                _level, _text = _optimism_banner_level(_cv_auc, _gap)
                if _level == "green":
                    st.success(f"**{_text}**")
                elif _level == "red":
                    st.error(
                        f"**High overfitting risk detected.** {_text} "
                        f"On genuine biomarker signal this gap is "
                        f"typically <0.03. A gap this size commonly "
                        f"occurs when feature selection memorizes "
                        f"noise-label correlations that nested CV "
                        f"correctly discounts. **Do not interpret the "
                        f"numbers below as validated clinical "
                        f"performance without independent-cohort "
                        f"verification (Section 8 validation workflow).**"
                    )
                else:  # amber
                    st.warning(f"**{_text}**")

                # Cross-reference line to the Panel Stability subsection
                # on the Panel tab. Streamlit anchor linking is fragile,
                # so we inline the Phi value here instead.
                _stab_cm = getattr(base, "panel_stability", None)
                if _stab_cm is not None and len(_stab_cm.per_fold_panels) > 0:
                    _phi_cm = _stab_cm.nogueira_stability
                    _lbl_cm = _stab_cm.benchmark_label
                    st.caption(
                        f"Related diagnostic — Panel stability: "
                        f"Φ={_phi_cm:.3f} ({_lbl_cm}). See Panel tab for "
                        f"full interpretation and per-gene frequencies."
                    )
                st.markdown("---")

            # Youden's + bootstrap CI
            yc_col, bc_col = st.columns(2)

            with yc_col:
                if "youdens" in clin:
                    y = clin["youdens"]
                    st.markdown("**Youden's optimal threshold**")
                    st.metric("Threshold", f"{y['threshold']:.3f}")
                    st.metric("Sensitivity", f"{y['sensitivity']:.3f}")
                    st.metric("Specificity", f"{y['specificity']:.3f}")
                    st.metric("Youden's J", f"{y['youdens_j']:.3f}")

            with bc_col:
                if "bootstrap_ci" in clin:
                    b = clin["bootstrap_ci"]
                    # M2 label consistency: match EnhancedBiomarkerResult.
                    # summary()'s convention exactly. When oof_used is
                    # True the bootstrap ran on (oof_true, oof_prob) so
                    # this CI is an OOF AUC CI. When oof_used is False
                    # (pre-M1 results or rare fallback), the bootstrap
                    # ran on training predictions so it's a training
                    # AUC CI. Symmetric explicit labels close the last
                    # sibling inconsistency after M2: summary() already
                    # labels it conditionally; the dashboard now does
                    # too. A test in
                    # test_biomarker_optimism_diagnostic_m2.py pins
                    # this consistency.
                    _oof_used = clin.get("oof_used", False)
                    _label = (
                        "**OOF AUC bootstrap 95% CI**" if _oof_used
                        else "**Training AUC bootstrap 95% CI**"
                    )
                    st.markdown(_label)
                    st.metric(
                        "Point estimate",
                        f"{b['point_estimate']:.3f}",
                    )
                    st.metric(
                        "95% CI",
                        f"[{b['ci_lower']:.3f}, {b['ci_upper']:.3f}]",
                    )
                    st.caption(f"n_bootstrap = {b['n_bootstrap']}")

            # PPV / NPV
            if "ppv_npv" in clin:
                st.markdown("**PPV / NPV at realistic prevalence**")
                p = clin["ppv_npv"]
                pp_col1, pp_col2, pp_col3 = st.columns(3)
                pp_col1.metric("Prevalence", f"{p['prevalence']:.1%}")
                pp_col2.metric("PPV", f"{p['ppv']:.3f}")
                pp_col3.metric("NPV", f"{p['npv']:.3f}")
                st.caption(
                    f"Derived from sensitivity={p['sensitivity']:.3f}, "
                    f"specificity={p['specificity']:.3f} via Bayes' theorem."
                )

            # Decision curve
            if "decision_curve" in clin:
                st.markdown("**Decision curve analysis**")
                dca = clin["decision_curve"]
                if isinstance(dca, pd.DataFrame) and not dca.empty:
                    # The default DCA sweeps threshold 0.01–0.99. At extreme
                    # thresholds the "treat all" term (event_rate − (1−p)·pt/(1−pt))
                    # goes sharply negative, which dominates the y-axis and
                    # makes the plot unreadable. Truncate to a clinically
                    # meaningful threshold window. Users who want the full
                    # curve can grab decision_curve.csv from the Downloads tab.
                    # Range slider avoids the low>high edge case that two
                    # independent number inputs would allow.
                    pt_low, pt_high = st.slider(
                        "Threshold probability range (plot window)",
                        min_value=0.01, max_value=0.99,
                        value=(0.05, 0.50), step=0.01,
                        help=(
                            "Net-benefit plot is only meaningful over a "
                            "clinically plausible range of threshold "
                            "probabilities. Default 0.05–0.50 covers most "
                            "diagnostic use cases."
                        ),
                        key="m10_dca_range",
                    )
                    dca_plot = dca[
                        (dca["threshold"] >= pt_low) &
                        (dca["threshold"] <= pt_high)
                    ].copy()

                    fig = go.Figure()
                    fig.add_trace(go.Scatter(
                        x=dca_plot["threshold"], y=dca_plot["net_benefit_model"],
                        mode="lines", name="Biomarker model",
                        line=dict(color="#2E86DE", width=3),
                    ))
                    fig.add_trace(go.Scatter(
                        x=dca_plot["threshold"], y=dca_plot["net_benefit_treat_all"],
                        mode="lines", name="Treat all",
                        line=dict(color="#E67E22", dash="dash"),
                    ))
                    fig.add_trace(go.Scatter(
                        x=dca_plot["threshold"], y=dca_plot["net_benefit_treat_none"],
                        mode="lines", name="Treat none",
                        line=dict(color="#7F8C8D", dash="dot"),
                    ))
                    fig.update_layout(
                        xaxis_title="Threshold probability",
                        yaxis_title="Net benefit",
                        title=(
                            f"Decision curve analysis "
                            f"(threshold {pt_low:.2f}–{pt_high:.2f})"
                        ),
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    st.caption(
                        "Net benefit at each threshold probability. "
                        "Biomarker model above 'Treat all' in the clinically "
                        "relevant range = net clinical value. "
                        "Full-range curve is in `decision_curve.csv`."
                    )

    # ---------- Tab 8: Ratio Biomarkers ----------
    with tabs[7]:
        st.markdown("### Ratio biomarkers (Top Scoring Pairs)")
        rr = getattr(result, "ratio_result", None) if is_enhanced else None
        if rr is None or rr.n_found == 0:
            st.info(
                "No ratio biomarkers found. "
                "Either ratios are not applicable to this intent or no pairs "
                "passed the minimum AUC threshold."
            )
        else:
            r1, r2, r3 = st.columns(3)
            r1.metric("Pairs found", rr.n_found)
            r2.metric("Pairs tested", f"{rr.n_pairs_tested:,}")
            if rr.best_pair:
                r3.metric("Best pair AUC", f"{rr.best_pair.auc:.3f}")

            ratio_df = rr.to_dataframe()
            st.dataframe(
                ratio_df.style.format({
                    "auc": "{:.3f}",
                    "mean_ratio_reference": "{:.3f}",
                    "mean_ratio_baseline": "{:.3f}",
                }),
                use_container_width=True,
                hide_index=True,
            )

            # Visualize top ratios
            top_n = min(10, len(ratio_df))
            top_df = ratio_df.head(top_n)
            fig = px.bar(
                top_df, x="auc", y="name", orientation="h",
                color="direction",
                color_discrete_map={
                    "higher_in_reference": "#E74C3C",
                    "higher_in_baseline":  "#3498DB",
                },
                title=f"Top {top_n} ratio biomarkers by AUC",
            )
            fig.update_layout(yaxis={'categoryorder': 'total ascending'})
            st.plotly_chart(fig, use_container_width=True)

    # ---------- Tab 9: Annotations ----------
    with tabs[8]:
        st.markdown("### Biological annotation")
        ann = base.annotation_result
        if ann is None:
            st.info("Annotation pipeline was not run.")
        else:
            a1, a2, a3 = st.columns(3)
            a1.metric("Genes annotated", ann.n_annotated)
            a2.metric("Enriched pathways (FDR<0.05)", ann.n_enriched_pathways)
            a3.metric(
                "Literature hits",
                len(ann.literature_hits) if hasattr(ann, "literature_hits")
                and ann.literature_hits is not None else 0,
            )

            # Gene info
            if hasattr(ann, "gene_info") and ann.gene_info is not None and not ann.gene_info.empty:
                st.markdown("**Gene annotations**")
                st.dataframe(ann.gene_info, use_container_width=True, height=300)

            # Pathway enrichment
            if hasattr(ann, "pathway_enrichment") and ann.pathway_enrichment is not None:
                pw = ann.pathway_enrichment
                if isinstance(pw, pd.DataFrame) and not pw.empty:
                    st.markdown("**Pathway enrichment**")
                    st.dataframe(pw.head(50), use_container_width=True, height=400)

            # PPI
            if hasattr(ann, "ppi_network") and ann.ppi_network is not None:
                n_edges = len(ann.ppi_network.get("edges", []))
                n_nodes = len(ann.ppi_network.get("nodes", []))
                st.caption(f"PPI network: {n_nodes} nodes, {n_edges} edges")

    # ---------- Tab 10: Downloads ----------
    with tabs[9]:
        st.markdown("### Download results")

        # Panel CSV
        panel_csv = pd.DataFrame({
            "rank": range(1, len(base.panel) + 1),
            "gene_id": base.panel,
        }).to_csv(index=False).encode("utf-8")
        st.download_button(
            "📥 biomarker_panel.csv",
            data=panel_csv,
            file_name="biomarker_panel.csv",
            mime="text/csv",
        )

        # Ranked genes CSV
        ranked_csv = base.ranked_genes.to_csv().encode("utf-8")
        st.download_button(
            "📥 ranked_genes.csv",
            data=ranked_csv,
            file_name="ranked_genes.csv",
            mime="text/csv",
        )

        # Classification performance CSV
        if base.classification_results:
            perf_rows_dl = []
            for name, res in base.classification_results.items():
                perf_rows_dl.append({
                    "classifier": name,
                    "AUC": getattr(res, "auc", float("nan")),
                    "F1": getattr(res, "f1", float("nan")),
                    "accuracy": getattr(res, "accuracy", float("nan")),
                    "sensitivity": getattr(res, "sensitivity", float("nan")),
                    "specificity": getattr(res, "specificity", float("nan")),
                    "is_best": (name == base.best_classifier),
                })
            perf_csv_dl = pd.DataFrame(perf_rows_dl).to_csv(index=False).encode("utf-8")
            st.download_button(
                "📥 classification_performance.csv",
                data=perf_csv_dl,
                file_name="classification_performance.csv",
                mime="text/csv",
                key="m10_dl_classification_perf_downloads_tab",
            )

        # Signature JSON
        if is_enhanced and getattr(result, "signature", None) is not None:
            import json
            sig_json = json.dumps(
                result.signature.to_dict(), indent=2, default=str
            ).encode("utf-8")
            st.download_button(
                "📥 signature_score.json",
                data=sig_json,
                file_name="signature_score.json",
                mime="application/json",
            )

        # Direction pattern CSV
        if is_enhanced and getattr(result, "direction_pattern", None) is not None:
            try:
                dp_csv = result.direction_pattern.to_dataframe().to_csv().encode("utf-8")
                st.download_button(
                    "📥 direction_pattern.csv",
                    data=dp_csv,
                    file_name="direction_pattern.csv",
                    mime="text/csv",
                )
            except Exception:
                pass

        # Ratio biomarkers CSV
        if is_enhanced and getattr(result, "ratio_result", None) is not None:
            rr = result.ratio_result
            if rr and rr.n_found > 0:
                ratios_csv = rr.to_dataframe().to_csv(index=False).encode("utf-8")
                st.download_button(
                    "📥 ratio_biomarkers.csv",
                    data=ratios_csv,
                    file_name="ratio_biomarkers.csv",
                    mime="text/csv",
                )

        # Decision curve CSV
        if is_enhanced:
            clin = getattr(result, "clinical_metrics", None)
            if clin and "decision_curve" in clin:
                dca = clin["decision_curve"]
                if isinstance(dca, pd.DataFrame) and not dca.empty:
                    dca_csv = dca.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        "📥 decision_curve.csv",
                        data=dca_csv,
                        file_name="decision_curve.csv",
                        mime="text/csv",
                    )

        st.caption(
            f"Results also saved to disk at: `{output_dir}/`. "
            "Use Module 8 (Reports) for a full biomarker discovery report."
        )

else:
    if ready:
        st.markdown("---")
        st.info(
            "Configure your analysis above and click **Run Biomarker Discovery** "
            "to populate results."
        )


# =============================================================================
# SECTION 8: VALIDATION (Scenario 7 — independent-cohort validation)
# =============================================================================
# Only visible once a discovery result exists. Lets the user upload a second
# cohort and apply the already-discovered panel (plus, for enhanced runs, the
# signature score and ratio biomarkers) to that cohort — producing side-by-side
# discovery-vs-validation comparison tables.
#
# This is the "does my biomarker replicate?" workflow. Without it, users would
# have to drop to the CLI; with it, validation is one upload away.
# =============================================================================

if st.session_state.get('m10_complete', False):
    st.markdown("---")
    st.markdown("## 8. Validation on Independent Cohort")
    st.caption(
        "Apply the discovered panel — and, when available, the signature score "
        "and ratio biomarkers — to an independent validation cohort. "
        "Validation AUC below discovery AUC is expected (mild overfitting); "
        "large drops suggest the panel may not generalize."
    )

    with st.expander("ℹ️ How validation works"):
        st.markdown("""
        **What this does:**
        1. Reads your validation counts + metadata
        2. Maps your validation groups to the same (baseline, reference) polarity
           as discovery
        3. Re-fits and cross-validates each classifier *from scratch* on the
           validation cohort using only the discovered panel genes
        4. If the discovery run had a diagnostic/prognostic intent:
           - Applies the stored signature score (weights + normalization) to
             validation samples and computes validation AUC of the score
           - Applies the discovered ratio-biomarker pairs and computes their
             validation AUC individually
           - Computes a direction pattern on validation data and reports
             concordance with discovery directions

        **What validation CAN'T tell you here:**
        - Whether the validation cohort has the same technical characteristics
          as discovery (batch, platform, library prep). Strong differences
          require upstream correction (M2 QC + ComBat).
        - Whether the effect is truly biological vs. a discovery-cohort artifact.
          A cleanly-replicating panel across independent cohorts is the best
          evidence, but not proof.
        """)

    val_discovery_result = st.session_state['m10_result']
    val_base = _base(val_discovery_result)
    val_is_enhanced = hasattr(val_discovery_result, "base_result")
    discovery_groups = st.session_state.get('m10_groups', {})
    discovery_baseline = discovery_groups.get('baseline')
    discovery_reference = discovery_groups.get('reference')
    discovery_panel = val_base.panel

    # --- Validation data source ---
    st.markdown("### 8.1 Validation Data")

    v_src_col1, v_src_col2 = st.columns(2)
    with v_src_col1:
        val_count_source = st.radio(
            "**Validation count matrix**",
            ["Upload CSV/TSV"],
            key="m10_val_count_source",
        )
    with v_src_col2:
        val_meta_source = st.radio(
            "**Validation metadata**",
            ["Upload CSV/TSV"],
            key="m10_val_meta_source",
        )

    val_counts_df = None
    val_meta_df = None

    val_count_file = st.file_uploader(
        "Validation count matrix (genes × samples)",
        type=["csv", "tsv", "txt"],
        key="m10_val_count_upload",
    )
    if val_count_file is not None:
        try:
            sep = '\t' if val_count_file.name.endswith(('.tsv', '.txt')) else ','
            val_counts_df = pd.read_csv(val_count_file, sep=sep, index_col=0)
            st.success(
                f"Loaded: **{val_counts_df.shape[0]} genes × "
                f"{val_counts_df.shape[1]} samples**"
            )
        except Exception as e:
            st.error(f"Failed to read validation count matrix: {e}")

    val_meta_file = st.file_uploader(
        "Validation metadata (samples × columns, sample IDs in first column)",
        type=["csv", "tsv", "txt"],
        key="m10_val_meta_upload",
    )
    if val_meta_file is not None:
        try:
            sep = '\t' if val_meta_file.name.endswith(('.tsv', '.txt')) else ','
            val_meta_df = pd.read_csv(val_meta_file, sep=sep, index_col=0)
            st.success(
                f"Loaded: **{val_meta_df.shape[0]} samples × "
                f"{val_meta_df.shape[1]} columns**"
            )
        except Exception as e:
            st.error(f"Failed to read validation metadata: {e}")

    # --- Group setup for validation cohort ---
    val_group_column = None
    val_baseline_group = None
    val_reference_group = None

    if val_meta_df is not None:
        st.markdown("### 8.2 Validation Group Setup")
        val_candidate_cols = [
            c for c in val_meta_df.columns if val_meta_df[c].nunique() >= 2
        ]
        if not val_candidate_cols:
            st.error(
                "No validation metadata column has at least 2 distinct values. "
                "Cannot define groups for validation."
            )
        else:
            # Prefer the same column name as discovery for consistency
            default_group_idx = 0
            disc_group_col = discovery_groups.get('group_column')
            if disc_group_col and disc_group_col in val_candidate_cols:
                default_group_idx = val_candidate_cols.index(disc_group_col)

            val_group_column = st.selectbox(
                "**Validation group column**",
                val_candidate_cols,
                index=default_group_idx,
                help=(
                    "Column in validation metadata defining the two groups. "
                    "Can differ from discovery column name, but the two groups "
                    "should correspond to the same biological contrast."
                ),
                key="m10_val_group_col",
            )

            val_unique_groups = sorted(
                val_meta_df[val_group_column].dropna().unique().tolist()
            )
            st.caption(
                f"Column `{val_group_column}` has {len(val_unique_groups)} "
                f"distinct value(s): "
                + ", ".join(f"`{g}`" for g in val_unique_groups[:10])
            )

            if len(val_unique_groups) < 2:
                st.error("Need at least 2 groups in validation data.")
            else:
                # Try to match discovery polarity automatically
                if discovery_baseline in val_unique_groups:
                    default_baseline_idx = val_unique_groups.index(discovery_baseline)
                else:
                    default_baseline_idx = 0

                vg_col1, vg_col2 = st.columns(2)
                with vg_col1:
                    val_baseline_group = st.selectbox(
                        "**Baseline group** (label = 0)",
                        val_unique_groups,
                        index=default_baseline_idx,
                        help=(
                            f"Match discovery baseline (`{discovery_baseline}`) "
                            "for consistent direction interpretation."
                        ),
                        key="m10_val_baseline",
                    )
                with vg_col2:
                    val_ref_options = [
                        g for g in val_unique_groups if g != val_baseline_group
                    ]
                    default_ref_idx = 0
                    if discovery_reference in val_ref_options:
                        default_ref_idx = val_ref_options.index(discovery_reference)
                    val_reference_group = st.selectbox(
                        "**Reference group** (label = 1)",
                        val_ref_options,
                        index=default_ref_idx,
                        help=(
                            f"Match discovery reference (`{discovery_reference}`) "
                            "for consistent direction interpretation."
                        ),
                        key="m10_val_reference",
                    )

                # Polarity alignment warning
                polarity_match = (
                    val_baseline_group == discovery_baseline
                    and val_reference_group == discovery_reference
                )
                if not polarity_match:
                    st.warning(
                        f"⚠️ Validation polarity (baseline=`{val_baseline_group}`, "
                        f"reference=`{val_reference_group}`) differs from "
                        f"discovery (baseline=`{discovery_baseline}`, "
                        f"reference=`{discovery_reference}`). "
                        "Direction-pattern concordance will reflect this "
                        "mismatch — ensure it is intentional."
                    )

                n_val_baseline = int(
                    (val_meta_df[val_group_column] == val_baseline_group).sum()
                )
                n_val_ref = int(
                    (val_meta_df[val_group_column] == val_reference_group).sum()
                )
                st.info(
                    f"Validation samples: **{n_val_baseline}** "
                    f"{val_baseline_group} vs **{n_val_ref}** {val_reference_group}"
                )

    # --- Panel coverage check ---
    if val_counts_df is not None:
        st.markdown("### 8.3 Panel Coverage in Validation Data")
        panel_in_val = [g for g in discovery_panel if g in val_counts_df.index]
        panel_missing = [g for g in discovery_panel if g not in val_counts_df.index]

        if panel_missing:
            st.warning(
                f"**{len(panel_missing)} of {len(discovery_panel)} panel genes "
                f"are missing from validation data:** "
                + ", ".join(f"`{g}`" for g in panel_missing)
                + ". Validation will proceed on the "
                f"**{len(panel_in_val)} available panel gene(s)**. "
                "This weakens the comparison — consider aligning gene IDs "
                "(e.g. via MyGene.info cross-referencing) before validating."
            )
        else:
            st.success(
                f"✅ All **{len(discovery_panel)}** panel genes present in "
                f"validation data: "
                + ", ".join(f"`{g}`" for g in discovery_panel)
            )

    # --- Run validation ---
    st.markdown("### 8.4 Run Validation")

    validation_ready = (
        val_counts_df is not None
        and val_meta_df is not None
        and val_group_column is not None
        and val_baseline_group is not None
        and val_reference_group is not None
        and val_baseline_group != val_reference_group
    )

    if not validation_ready:
        missing_val = []
        if val_counts_df is None:
            missing_val.append("validation count matrix")
        if val_meta_df is None:
            missing_val.append("validation metadata")
        if val_baseline_group == val_reference_group:
            missing_val.append("distinct baseline and reference groups")
        st.warning(
            f"**Not ready to validate.** Missing: {', '.join(missing_val)}."
        )

    validation_clicked = st.button(
        "🔬 Run Validation",
        type="primary",
        disabled=not validation_ready,
        use_container_width=True,
        key="m10_validation_button",
    )

    if validation_clicked and validation_ready:
        val_progress = st.progress(0, text="Preparing validation…")
        val_status = st.empty()

        try:
            # --- (a) Base-pipeline validation: cross-validated classification
            #         on the panel genes within the validation cohort ---
            val_status.info("Running classifiers on validation cohort…")
            val_progress.progress(20, text="Running validation classifiers…")

            val_clf_results = validate_biomarkers(
                panel_genes=discovery_panel,
                counts=val_counts_df,
                metadata=val_meta_df,
                group_column=val_group_column,
                baseline_group=val_baseline_group,
                reference_group=val_reference_group,
                n_folds=min(5, int(min(n_val_baseline, n_val_ref) / 2)),
                random_state=int(random_state),
                verbose=False,
            )

            val_progress.progress(60, text="Applying enhancement modules…")

            # --- (b) Enhancement modules (if the discovery run was enhanced) ---
            val_signature_auc = None
            val_ratio_aucs = None
            val_direction_concordance = None

            if val_is_enhanced:
                # Build log2-CPM on validation data, aligned to discovery sample IDs
                val_lib = val_counts_df.sum(axis=0)
                val_cpm = val_counts_df.div(val_lib, axis=1) * 1e6
                val_log_cpm = np.log2(val_cpm + 1).T  # samples x genes

                # Align to metadata samples
                val_meta_sample_ids = val_meta_df.index.astype(str).tolist()
                common = [s for s in val_log_cpm.index if s in val_meta_sample_ids]
                val_log_cpm_aligned = val_log_cpm.loc[common]

                # Binary labels for validation
                val_labels_aligned = (
                    val_meta_df.loc[common, val_group_column] == val_reference_group
                ).astype(int).to_numpy()

                # Signature score validation
                sig = getattr(val_discovery_result, "signature", None)
                if sig is not None:
                    try:
                        # Subset to panel genes actually present
                        sig_genes_in_val = [
                            g for g in sig.panel_genes if g in val_log_cpm_aligned.columns
                        ]
                        if len(sig_genes_in_val) == len(sig.panel_genes):
                            from sklearn.metrics import roc_auc_score
                            val_scores = sig.score(val_log_cpm_aligned)
                            val_signature_auc = float(roc_auc_score(
                                val_labels_aligned, val_scores.values
                            ))
                    except Exception as e:
                        st.warning(f"Signature-score validation failed: {e}")

                # Ratio biomarker validation
                rr = getattr(val_discovery_result, "ratio_result", None)
                if rr is not None and rr.n_found > 0:
                    try:
                        from sklearn.metrics import roc_auc_score
                        # Need the RAW (not log) expression for ratios, since
                        # ratios of log values are less interpretable.
                        # Use CPM (linear scale) with pseudocount handled inside apply_ratios.
                        val_ratio_features = apply_ratios(
                            val_cpm.T.loc[common],  # samples x genes, linear CPM
                            rr.pairs,
                            pseudocount=1.0,
                        )
                        val_ratio_aucs = []
                        for pair in rr.pairs:
                            ratio_vals = val_ratio_features[pair.name].values
                            auc = float(roc_auc_score(val_labels_aligned, ratio_vals))
                            if auc < 0.5:
                                auc = 1 - auc
                            val_ratio_aucs.append({
                                "pair": pair.name,
                                "discovery_auc": pair.auc,
                                "validation_auc": auc,
                            })
                    except Exception as e:
                        st.warning(f"Ratio-biomarker validation failed: {e}")

                # Direction pattern concordance
                dp = getattr(val_discovery_result, "direction_pattern", None)
                if dp is not None:
                    try:
                        str_labels_val = pd.Series(
                            [val_reference_group if lbl == 1 else val_baseline_group
                             for lbl in val_labels_aligned],
                            index=common,
                        )
                        # DirectionPattern exposes .genes (property), not
                        # .panel_genes — that's a SignatureScore attribute.
                        # Conflating the two was the Scenario 7 bug.
                        dp_genes = dp.genes
                        val_dp = build_direction_pattern(
                            expression=val_log_cpm_aligned[dp_genes],
                            labels=str_labels_val,
                            reference_group=val_reference_group,
                            baseline_group=val_baseline_group,
                            p_threshold=0.05,
                            fc_threshold=0.0,
                        )
                        # Compare directions gene-by-gene
                        agree = 0
                        total = 0
                        disagreements = []
                        filtered_out = []
                        for g in dp_genes:
                            if g in val_dp.gene_directions:
                                total += 1
                                if dp.gene_directions[g] == val_dp.gene_directions[g]:
                                    agree += 1
                                else:
                                    disagreements.append(
                                        (g, dp.gene_directions[g],
                                         val_dp.gene_directions[g])
                                    )
                            else:
                                # Panel gene didn't reach p<0.05 in validation —
                                # excluded from concordance statistic. Capture
                                # the validation-cohort log2FC for user context.
                                val_fc = val_log_cpm_aligned.loc[
                                    val_labels_aligned == 1, g
                                ].mean() - val_log_cpm_aligned.loc[
                                    val_labels_aligned == 0, g
                                ].mean()
                                filtered_out.append({
                                    "gene": g,
                                    "discovery_direction": dp.gene_directions[g],
                                    "validation_log2FC": float(val_fc),
                                })
                        val_direction_concordance = {
                            "agree": agree,
                            "total": total,
                            "fraction": agree / total if total > 0 else 0,
                            "disagreements": disagreements,
                            "filtered_out": filtered_out,
                            "n_panel": len(dp_genes),
                        }
                    except Exception as e:
                        st.warning(f"Direction-pattern validation failed: {e}")

            val_progress.progress(100, text="Done.")
            val_status.empty()

            # Persist results for re-display
            st.session_state['m10_validation_clf'] = val_clf_results
            st.session_state['m10_validation_signature_auc'] = val_signature_auc
            st.session_state['m10_validation_ratio_aucs'] = val_ratio_aucs
            st.session_state['m10_validation_direction'] = val_direction_concordance
            st.session_state['m10_validation_complete'] = True
            st.session_state['m10_val_panel_coverage'] = (
                len(panel_in_val), len(discovery_panel)
            )

            st.success("**Validation complete.** Results below.")

        except Exception as e:
            val_progress.empty()
            val_status.empty()
            st.error(
                f"**Validation failed.**\n\nError: `{e}`\n\n"
                "**Troubleshooting:**\n"
                "- Ensure validation metadata sample IDs match counts columns\n"
                "- Check that panel genes are present in validation counts\n"
                "- Verify baseline/reference selection matches discovery polarity"
            )

    # --- Display validation results ---
    if st.session_state.get('m10_validation_complete', False):
        st.markdown("### 8.5 Validation Results")

        val_clf = st.session_state['m10_validation_clf']
        val_sig_auc = st.session_state.get('m10_validation_signature_auc')
        val_ratio_aucs = st.session_state.get('m10_validation_ratio_aucs')
        val_dir = st.session_state.get('m10_validation_direction')
        coverage_n, coverage_total = st.session_state.get(
            'm10_val_panel_coverage', (0, 0)
        )

        val_tabs = st.tabs([
            "Classification", "Signature Score", "Ratio Biomarkers",
            "Direction Concordance", "Summary",
        ])

        # --- Classification tab: Discovery vs Validation side-by-side ---
        with val_tabs[0]:
            st.markdown("#### Classification performance: Discovery vs Validation")
            rows = []
            for clf_name in val_clf.keys():
                disc_res = val_base.classification_results.get(clf_name)
                val_res = val_clf[clf_name]
                rows.append({
                    "classifier": clf_name,
                    "discovery_AUC": (
                        getattr(disc_res, "auc", float("nan"))
                        if disc_res is not None else float("nan")
                    ),
                    "validation_AUC": getattr(val_res, "auc", float("nan")),
                    "discovery_sens": (
                        getattr(disc_res, "sensitivity", float("nan"))
                        if disc_res is not None else float("nan")
                    ),
                    "validation_sens": getattr(val_res, "sensitivity", float("nan")),
                    "discovery_spec": (
                        getattr(disc_res, "specificity", float("nan"))
                        if disc_res is not None else float("nan")
                    ),
                    "validation_spec": getattr(val_res, "specificity", float("nan")),
                })
            comp_df = pd.DataFrame(rows)
            comp_df["AUC_drop"] = (
                comp_df["discovery_AUC"] - comp_df["validation_AUC"]
            )
            st.dataframe(
                comp_df.style.format({
                    "discovery_AUC": "{:.3f}", "validation_AUC": "{:.3f}",
                    "discovery_sens": "{:.3f}", "validation_sens": "{:.3f}",
                    "discovery_spec": "{:.3f}", "validation_spec": "{:.3f}",
                    "AUC_drop": "{:+.3f}",
                }),
                use_container_width=True,
                hide_index=True,
            )
            st.caption(
                "AUC_drop = discovery − validation. Small positive drops (< 0.10) "
                "are typical mild overfitting. Negative drops (validation > "
                "discovery) can happen when the validation cohort is better-matched "
                "to the panel's signal than the discovery set."
            )

            # Discovery vs Validation AUC bar chart
            plot_df = comp_df.melt(
                id_vars="classifier",
                value_vars=["discovery_AUC", "validation_AUC"],
                var_name="cohort",
                value_name="AUC",
            )
            plot_df["cohort"] = plot_df["cohort"].str.replace("_AUC", "")
            fig = px.bar(
                plot_df, x="classifier", y="AUC", color="cohort",
                barmode="group",
                title="AUC: Discovery vs Validation",
                color_discrete_map={"discovery": "#2E86DE", "validation": "#E67E22"},
            )
            fig.update_yaxes(range=[0, 1.05])
            st.plotly_chart(fig, use_container_width=True)

            # Download button
            dl_csv = comp_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "📥 validation_classification_comparison.csv",
                data=dl_csv,
                file_name="validation_classification_comparison.csv",
                mime="text/csv",
                key="m10_dl_val_comparison",
            )

        # --- Signature Score tab ---
        with val_tabs[1]:
            st.markdown("#### Signature score on validation cohort")
            if val_sig_auc is None:
                st.info(
                    "Signature score not available. "
                    "Discovery run must use diagnostic/prognostic intent "
                    "for this to be computed."
                )
            else:
                sig = getattr(val_discovery_result, "signature", None)
                disc_sig_auc = sig.performance.get("auc", float("nan")) if sig else float("nan")
                sc1, sc2, sc3 = st.columns(3)
                sc1.metric("Discovery signature AUC", f"{disc_sig_auc:.3f}")
                sc2.metric("Validation signature AUC", f"{val_sig_auc:.3f}")
                sc3.metric("Drop", f"{disc_sig_auc - val_sig_auc:+.3f}")
                st.caption(
                    "Signature scores applied to validation cohort using the "
                    "stored per-gene z-score normalization and panel weights. "
                    "Validation AUC reflects how well the pre-fit scoring "
                    "function generalizes — no re-fitting on validation data."
                )

        # --- Ratio Biomarkers tab ---
        with val_tabs[2]:
            st.markdown("#### Ratio biomarkers on validation cohort")
            if val_ratio_aucs is None or len(val_ratio_aucs) == 0:
                st.info(
                    "No ratio biomarkers to validate. "
                    "Discovery run must find ratio pairs for this to display."
                )
            else:
                r_df = pd.DataFrame(val_ratio_aucs)
                r_df["AUC_drop"] = r_df["discovery_auc"] - r_df["validation_auc"]
                st.dataframe(
                    r_df.style.format({
                        "discovery_auc": "{:.3f}",
                        "validation_auc": "{:.3f}",
                        "AUC_drop": "{:+.3f}",
                    }),
                    use_container_width=True,
                    hide_index=True,
                )
                st.caption(
                    "Each ratio pair evaluated on validation data as a standalone "
                    "discriminator. AUC_drop shows how much each ratio "
                    "generalizes or fails to generalize."
                )

        # --- Direction Concordance tab ---
        with val_tabs[3]:
            st.markdown("#### Direction pattern concordance")
            if val_dir is None:
                st.info(
                    "Direction pattern not available. "
                    "Requires a discovery run with diagnostic/exploratory/"
                    "predictive/translational intent."
                )
            else:
                n_panel = val_dir.get('n_panel', val_dir['total'])
                filtered_out = val_dir.get('filtered_out', [])

                dc1, dc2, dc3 = st.columns(3)
                dc1.metric(
                    "Genes in agreement",
                    f"{val_dir['agree']}/{val_dir['total']}",
                    help=(
                        "Numerator: discovery direction matches validation "
                        "direction. Denominator: panel genes that reached "
                        "p<0.05 in validation cohort."
                    ),
                )
                dc2.metric("Concordance", f"{val_dir['fraction']:.1%}")
                from scipy.stats import binomtest
                try:
                    p_value = binomtest(
                        val_dir['agree'], val_dir['total'], p=0.5,
                        alternative='greater',
                    ).pvalue
                    dc3.metric("p-value (binomial)", f"{p_value:.3g}")
                except Exception:
                    dc3.metric("p-value (binomial)", "n/a")

                # Context: were any panel genes filtered out by p-threshold?
                if filtered_out:
                    st.warning(
                        f"**{len(filtered_out)} of {n_panel} panel gene(s) did "
                        f"not reach p<0.05 in the validation cohort** and are "
                        "excluded from the concordance statistic above. "
                        "This is common when validation cohorts are smaller "
                        "or individual gene effects are subtle. Table below "
                        "shows the observed validation-cohort fold change for "
                        "each filtered gene — directionally consistent "
                        "effects that just miss significance are different "
                        "from genuinely lost signal."
                    )
                    filt_df = pd.DataFrame(filtered_out)
                    st.dataframe(
                        filt_df.style.format({"validation_log2FC": "{:+.2f}"}),
                        use_container_width=True,
                        hide_index=True,
                    )

                if val_dir['disagreements']:
                    st.markdown("**Genes with discordant direction** (significant in both cohorts, direction flipped):")
                    dis_df = pd.DataFrame(
                        val_dir['disagreements'],
                        columns=["gene", "discovery_direction", "validation_direction"],
                    )
                    st.dataframe(dis_df, hide_index=True, use_container_width=True)
                else:
                    st.success(
                        f"✅ All {val_dir['total']} panel gene(s) with significant "
                        "validation signal show the same direction as in discovery."
                    )
                st.caption(
                    "Binomial p-value tests H0: genes agree by chance (p=0.5). "
                    "Significant p indicates concordance exceeds chance. "
                    "Small denominators limit power — 4/4 gives p=0.0625, 5/5 gives p=0.031."
                )

        # --- Summary tab ---
        with val_tabs[4]:
            st.markdown("#### Validation summary")
            lines = [
                f"**Discovery panel:** {len(discovery_panel)} genes "
                f"({coverage_n}/{coverage_total} available in validation)",
                f"**Discovery groups:** baseline=`{discovery_baseline}`, "
                f"reference=`{discovery_reference}`",
                f"**Validation groups:** baseline=`{val_baseline_group}`, "
                f"reference=`{val_reference_group}`",
            ]

            # Best classifier summary
            best_val_clf = max(val_clf.keys(), key=lambda k: val_clf[k].auc)
            best_val_auc = val_clf[best_val_clf].auc
            best_disc_clf = val_base.best_classifier
            best_disc_auc = val_base.classification_results[best_disc_clf].auc
            lines.append("")
            lines.append(
                f"**Best discovery classifier:** `{best_disc_clf}` "
                f"(AUC = {best_disc_auc:.3f})"
            )
            lines.append(
                f"**Best validation classifier:** `{best_val_clf}` "
                f"(AUC = {best_val_auc:.3f})"
            )
            drop = best_disc_auc - best_val_auc
            lines.append(f"**AUC drop:** {drop:+.3f}")
            if drop < 0.05:
                lines.append("*Interpretation:* Panel replicates strongly "
                             "(AUC drop < 0.05).")
            elif drop < 0.15:
                lines.append("*Interpretation:* Panel replicates with mild "
                             "overfitting (AUC drop 0.05–0.15).")
            else:
                lines.append("*Interpretation:* Panel shows substantial drop "
                             "(AUC drop > 0.15). Consider whether the "
                             "discovery cohort had systematic differences.")

            if val_sig_auc is not None:
                sig = getattr(val_discovery_result, "signature", None)
                disc_sig_auc = sig.performance.get("auc", float("nan")) if sig else float("nan")
                lines.append("")
                lines.append(
                    f"**Signature score:** discovery AUC = {disc_sig_auc:.3f}, "
                    f"validation AUC = {val_sig_auc:.3f}"
                )

            if val_dir is not None:
                lines.append("")
                n_panel_total = val_dir.get('n_panel', val_dir['total'])
                if val_dir['total'] < n_panel_total:
                    lines.append(
                        f"**Direction concordance:** {val_dir['agree']}/"
                        f"{val_dir['total']} genes agree "
                        f"({val_dir['fraction']:.1%}); "
                        f"{n_panel_total - val_dir['total']} of {n_panel_total} "
                        f"panel genes did not reach p<0.05 in validation."
                    )
                else:
                    lines.append(
                        f"**Direction concordance:** {val_dir['agree']}/"
                        f"{val_dir['total']} genes agree "
                        f"({val_dir['fraction']:.1%})"
                    )

            st.markdown("\n".join(lines))