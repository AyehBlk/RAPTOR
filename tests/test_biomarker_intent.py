"""
Unit tests for raptor.biomarker_discovery.intent

Covers:
    - All 6 intent modes construct correctly
    - Invalid intent strings raise ValueError
    - Translational intent requires species arguments
    - Non-translational intents warn when species arguments are given
    - Data validation catches missing required columns
    - describe() and to_dict() produce expected output
    - custom_label overrides the default output_label
"""

import pytest
import warnings

from raptor.biomarker_discovery import BiomarkerIntent, VALID_INTENTS


# ---------------------------------------------------------------------------
# Basic construction
# ---------------------------------------------------------------------------

class TestBasicConstruction:
    """Each valid intent should construct without error and populate defaults."""

    @pytest.mark.parametrize("intent_name", VALID_INTENTS)
    def test_all_intents_construct(self, intent_name):
        """Every valid intent name should produce a working instance."""
        # Translational needs species args; give defaults so this param test works.
        if intent_name == "translational":
            intent = BiomarkerIntent(
                intent_name,
                source_species="mouse",
                target_species="human",
            )
        else:
            intent = BiomarkerIntent(intent_name)

        assert intent.intent == intent_name
        assert intent.study_design  # non-empty
        assert intent.validation_strategy  # non-empty
        assert len(intent.required_data_columns) >= 1
        assert intent.minimum_samples > 0
        assert len(intent.required_metrics) >= 1
        assert intent.output_label  # non-empty

    def test_diagnostic_defaults(self):
        intent = BiomarkerIntent("diagnostic")
        assert intent.study_design == "case_control"
        assert intent.validation_strategy == "nested_cv"
        assert "AUC" in intent.required_metrics
        assert intent.minimum_samples == 30

    def test_prognostic_defaults(self):
        intent = BiomarkerIntent("prognostic")
        assert intent.validation_strategy == "cox_cv"
        assert "time" in intent.required_data_columns
        assert "event" in intent.required_data_columns
        assert "c_index" in intent.required_metrics

    def test_translational_defaults(self):
        intent = BiomarkerIntent(
            "translational",
            source_species="mouse",
            target_species="human",
        )
        assert intent.source_species == "mouse"
        assert intent.target_species == "human"
        assert "ortholog_map" in intent.required_data_columns
        assert "direction_concordance" in intent.required_metrics


# ---------------------------------------------------------------------------
# Validation of input
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_invalid_intent_raises(self):
        with pytest.raises(ValueError, match="Invalid intent"):
            BiomarkerIntent("nonsense")

    def test_translational_requires_source_species(self):
        with pytest.raises(ValueError, match="source_species"):
            BiomarkerIntent("translational", target_species="human")

    def test_translational_requires_target_species(self):
        with pytest.raises(ValueError, match="target_species|source_species"):
            BiomarkerIntent("translational", source_species="mouse")

    def test_species_args_warn_for_nontranslational(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            BiomarkerIntent("diagnostic", source_species="mouse")
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)
            assert "translational" in str(w[0].message)


# ---------------------------------------------------------------------------
# Data validation helper
# ---------------------------------------------------------------------------

class TestValidateData:

    def test_valid_columns_pass(self):
        intent = BiomarkerIntent("diagnostic")
        # Should not raise
        intent.validate_data(["expression", "group", "extra_col"])

    def test_missing_column_raises(self):
        intent = BiomarkerIntent("prognostic")
        with pytest.raises(ValueError, match="missing"):
            intent.validate_data(["expression", "time"])  # missing 'event'

    def test_empty_columns_raises(self):
        intent = BiomarkerIntent("diagnostic")
        with pytest.raises(ValueError, match="missing"):
            intent.validate_data([])


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

class TestOutputHelpers:

    def test_describe_returns_multiline_string(self):
        intent = BiomarkerIntent("diagnostic")
        text = intent.describe()
        assert "Biomarker Intent: diagnostic" in text
        assert "Output label" in text
        assert "\n" in text  # multi-line

    def test_describe_translational_includes_species(self):
        intent = BiomarkerIntent(
            "translational",
            source_species="mouse",
            target_species="human",
        )
        text = intent.describe()
        assert "mouse" in text
        assert "human" in text

    def test_to_dict_contains_all_fields(self):
        intent = BiomarkerIntent("diagnostic")
        d = intent.to_dict()
        expected_keys = {
            "intent", "source_species", "target_species", "custom_label",
            "study_design", "validation_strategy", "required_data_columns",
            "minimum_samples", "required_metrics", "output_label",
        }
        assert set(d.keys()) == expected_keys

    def test_custom_label_overrides_default(self):
        intent = BiomarkerIntent(
            "diagnostic",
            custom_label="My Custom Panel",
        )
        assert intent.output_label == "My Custom Panel"

    def test_custom_label_none_uses_default(self):
        intent = BiomarkerIntent("diagnostic")
        assert intent.output_label == "Diagnostic Biomarker Candidates"