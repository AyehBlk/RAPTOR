"""
Biomarker Discovery subpackage (Module 10).

Layout:
    core.py              -- original biomarker_discovery module (all existing logic)
    intent.py            -- intent classification
    signature_score.py   -- weighted risk score engine (coming next)

Backward compatibility: every public name from the old
`raptor/biomarker_discovery.py` file is re-exported here, so every existing
import like `from raptor.biomarker_discovery import discover_biomarkers`
continues to work unchanged.
"""

# Re-export everything public from core.py
from raptor.biomarker_discovery.core import *  # noqa: F401, F403

# Explicit re-export of M6-new public names so IDE autocomplete and
# `from raptor.biomarker_discovery import PanelStabilityResult` both work.
from raptor.biomarker_discovery.core import (
    PanelStabilityResult,
)  # noqa: F401

# Explicit re-export of M4-new public name.
from raptor.biomarker_discovery.core import (
    apply_significance_calibration,
)  # noqa: F401

# M4: shared univariate DE utility. Public so PTERO and external scripts
# can reuse it.
from raptor.biomarker_discovery.univariate_de import (
    compute_per_gene_de,
    VALID_TESTS as UNIVARIATE_DE_VALID_TESTS,
)  # noqa: F401

# Re-export private helpers and availability flags that tests and internal code need.
# These are underscore-prefixed, so `import *` above skips them — we must list them explicitly.
from raptor.biomarker_discovery.core import (
    _prepare_expression_data,
    _BORUTA_AVAILABLE,
    _MRMR_AVAILABLE,
    _SHAP_AVAILABLE,
    # M6: pipeline-CV orchestrator and stability helpers
    _run_pipeline_cv,
    _run_feature_selection,
    _resolve_pipeline_cv_methods,
    _nogueira_stability_from_matrix,
    _panels_to_selection_matrix,
    _bootstrap_nogueira_ci,
    _nogueira_benchmark_label,
    _compute_panel_stability,
    _NOGUEIRA_EXCELLENT_THRESHOLD,
    _NOGUEIRA_POOR_THRESHOLD,
    _FOLD_SAFE_DEFAULT_METHODS,
)  # noqa: F401

# M10 enhancement modules
from raptor.biomarker_discovery.intent import (
    BiomarkerIntent,
    VALID_INTENTS,
)  # noqa: F401

from raptor.biomarker_discovery.signature_score import (
    SignatureScore,
    build_signature_score,
)  # noqa: F401

from raptor.biomarker_discovery.direction_patterns import (
    DirectionPattern,
    build_direction_pattern,
    VALID_DIRECTIONS,
    DIRECTION_UP,
    DIRECTION_DOWN,
)  # noqa: F401

from raptor.biomarker_discovery.clinical_metrics import (
    ppv_npv_at_prevalence,
    bootstrap_ci,
    youdens_optimal_threshold,
    decision_curve_analysis,
    net_reclassification_improvement,
)  # noqa: F401

from raptor.biomarker_discovery.ratio_biomarkers import (
    RatioBiomarkerSearcher,
    RatioPair,
    RatioSearchResult,
    apply_ratios,
    build_ratio_features,
)  # noqa: F401

from raptor.biomarker_discovery.enhanced import (
    EnhancedBiomarkerResult,
    enhance_biomarker_result,
)  # noqa: F401