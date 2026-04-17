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

# Re-export private helpers and availability flags that tests and internal code need.
# These are underscore-prefixed, so `import *` above skips them — we must list them explicitly.
from raptor.biomarker_discovery.core import (
    _prepare_expression_data,
    _BORUTA_AVAILABLE,
    _MRMR_AVAILABLE,
    _SHAP_AVAILABLE,
)  # noqa: F401

# New M10 enhancement modules
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