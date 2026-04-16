"""
Biomarker Intent Classification
================================

Defines the *goal* of a biomarker discovery analysis. Different goals require
different statistical defaults, validation strategies, and output labels.

Six intent modes are supported:
    - diagnostic:    Does the patient have the disease? (classification)
    - prognostic:    How will the patient fare over time? (survival)
    - predictive:    How will the patient respond to treatment? (interaction)
    - monitoring:    Is the disease progressing or regressing? (longitudinal)
    - exploratory:   What patterns exist? (unsupervised, no preset outcome)
    - translational: Will mouse findings replicate in humans? (cross-species)

Usage
-----
>>> intent = BiomarkerIntent("diagnostic")
>>> intent.required_metrics
['AUC', 'sensitivity', 'specificity', 'PPV', 'NPV']
>>> intent.minimum_samples
30

>>> intent = BiomarkerIntent("translational", source_species="mouse",
...                          target_species="human")
>>> intent.required_data_columns
['expression', 'group', 'ortholog_map']
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


# ---------------------------------------------------------------------------
# Supported intent modes (single source of truth)
# ---------------------------------------------------------------------------

VALID_INTENTS = (
    "diagnostic",
    "prognostic",
    "predictive",
    "monitoring",
    "exploratory",
    "translational",
)


# ---------------------------------------------------------------------------
# Per-intent default configurations
# ---------------------------------------------------------------------------
# Each entry describes: required data columns, minimum sample count,
# which metrics should appear in the output report, and a human-readable
# label used in reports and dashboards.
#
# These are defaults, not hard rules -- users can override any field.

_INTENT_DEFAULTS: Dict[str, Dict[str, Any]] = {
    "diagnostic": {
        "study_design": "case_control",
        "validation_strategy": "nested_cv",
        "required_data_columns": ["expression", "group"],
        "minimum_samples": 30,
        "required_metrics": ["AUC", "sensitivity", "specificity", "PPV", "NPV"],
        "output_label": "Diagnostic Biomarker Candidates",
    },
    "prognostic": {
        "study_design": "cohort",
        "validation_strategy": "cox_cv",
        "required_data_columns": ["expression", "time", "event"],
        "minimum_samples": 50,
        "required_metrics": ["c_index", "HR", "HR_95CI", "logrank_p"],
        "output_label": "Prognostic Biomarker Candidates",
    },
    "predictive": {
        "study_design": "interaction",
        "validation_strategy": "treatment_interaction_cv",
        "required_data_columns": ["expression", "treatment", "response"],
        "minimum_samples": 40,
        "required_metrics": ["interaction_p", "subgroup_AUC", "HR_by_arm"],
        "output_label": "Predictive (Treatment-Response) Biomarker Candidates",
    },
    "monitoring": {
        "study_design": "longitudinal",
        "validation_strategy": "mixed_effects_cv",
        "required_data_columns": ["expression", "subject_id", "timepoint"],
        "minimum_samples": 20,  # per subject min 2 timepoints assumed
        "required_metrics": ["within_subject_change", "icc", "auc_trajectory"],
        "output_label": "Disease Monitoring Biomarker Candidates",
    },
    "exploratory": {
        "study_design": "unsupervised",
        "validation_strategy": "bootstrap_stability",
        "required_data_columns": ["expression"],
        "minimum_samples": 20,
        "required_metrics": ["cluster_stability", "explained_variance"],
        "output_label": "Exploratory Biomarker Candidates",
    },
    "translational": {
        "study_design": "cross_species",
        "validation_strategy": "ortholog_concordance",
        "required_data_columns": ["expression", "group", "ortholog_map"],
        "minimum_samples": 20,
        "required_metrics": ["direction_concordance", "ortholog_coverage", "cross_species_AUC"],
        "output_label": "Translational Biomarker Candidates",
    },
}


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

@dataclass
class BiomarkerIntent:
    """
    Describes the goal of a biomarker discovery analysis.

    Parameters
    ----------
    intent : str
        One of: diagnostic, prognostic, predictive, monitoring,
        exploratory, translational.
    source_species : str, optional
        For 'translational' intent, the species data was collected in
        (e.g. 'mouse'). Ignored for other intents.
    target_species : str, optional
        For 'translational' intent, the species to translate findings to
        (e.g. 'human'). Ignored for other intents.
    custom_label : str, optional
        Override the default output label shown in reports and dashboards.

    Attributes
    ----------
    study_design : str
        Statistical study design implied by the intent.
    validation_strategy : str
        Recommended cross-validation approach.
    required_data_columns : List[str]
        Column names that must exist in the input metadata.
    minimum_samples : int
        Minimum total sample count for reliable results.
    required_metrics : List[str]
        Metrics that should appear in the final report.
    output_label : str
        Human-readable label for reports and dashboards.
    """

    intent: str
    source_species: Optional[str] = None
    target_species: Optional[str] = None
    custom_label: Optional[str] = None

    # Fields filled in by __post_init__ from the defaults table
    study_design: str = field(init=False)
    validation_strategy: str = field(init=False)
    required_data_columns: List[str] = field(init=False)
    minimum_samples: int = field(init=False)
    required_metrics: List[str] = field(init=False)
    output_label: str = field(init=False)

    def __post_init__(self) -> None:
        # --- Validate intent value ---
        if self.intent not in VALID_INTENTS:
            raise ValueError(
                f"Invalid intent {self.intent!r}. "
                f"Must be one of: {', '.join(VALID_INTENTS)}."
            )

        # --- Validate translational-specific arguments ---
        if self.intent == "translational":
            if not self.source_species or not self.target_species:
                raise ValueError(
                    "Translational intent requires both 'source_species' "
                    "and 'target_species' (e.g. source_species='mouse', "
                    "target_species='human')."
                )
        else:
            # Silently ignore species args for non-translational intents,
            # but warn if they were set unnecessarily.
            if self.source_species or self.target_species:
                import warnings
                warnings.warn(
                    f"source_species/target_species are only used with "
                    f"intent='translational' (got intent={self.intent!r}). "
                    f"These arguments will be ignored.",
                    UserWarning,
                    stacklevel=2,
                )

        # --- Load defaults from the table ---
        defaults = _INTENT_DEFAULTS[self.intent]
        self.study_design = defaults["study_design"]
        self.validation_strategy = defaults["validation_strategy"]
        self.required_data_columns = list(defaults["required_data_columns"])  # copy
        self.minimum_samples = defaults["minimum_samples"]
        self.required_metrics = list(defaults["required_metrics"])  # copy
        self.output_label = self.custom_label or defaults["output_label"]

    # -----------------------------------------------------------------------
    # Helper methods
    # -----------------------------------------------------------------------

    def describe(self) -> str:
        """Return a human-readable multi-line summary of this intent."""
        lines = [
            f"Biomarker Intent: {self.intent}",
            f"  Output label:         {self.output_label}",
            f"  Study design:         {self.study_design}",
            f"  Validation strategy:  {self.validation_strategy}",
            f"  Required columns:     {', '.join(self.required_data_columns)}",
            f"  Minimum samples:      {self.minimum_samples}",
            f"  Required metrics:     {', '.join(self.required_metrics)}",
        ]
        if self.intent == "translational":
            lines.insert(1, f"  Source species:       {self.source_species}")
            lines.insert(2, f"  Target species:       {self.target_species}")
        return "\n".join(lines)

    def validate_data(self, metadata_columns: List[str]) -> None:
        """
        Check that `metadata_columns` contains all required columns for this
        intent. Raises ValueError listing any missing columns.
        """
        missing = [c for c in self.required_data_columns if c not in metadata_columns]
        if missing:
            raise ValueError(
                f"Intent {self.intent!r} requires columns {self.required_data_columns} "
                f"but the following are missing from metadata: {missing}"
            )

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dict (useful for JSON logs or dashboard state)."""
        return {
            "intent": self.intent,
            "source_species": self.source_species,
            "target_species": self.target_species,
            "custom_label": self.custom_label,
            "study_design": self.study_design,
            "validation_strategy": self.validation_strategy,
            "required_data_columns": self.required_data_columns,
            "minimum_samples": self.minimum_samples,
            "required_metrics": self.required_metrics,
            "output_label": self.output_label,
        }