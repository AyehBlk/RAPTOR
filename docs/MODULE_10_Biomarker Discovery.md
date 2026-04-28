# RAPTOR Module 10: Biomarker Discovery - Comprehensive User Guide

**Version:** 2.2.2  
**Author:** Ayeh Bolouki  
**Last Updated:** April 2026  

---

## Table of Contents

1. [Overview](#overview)
2. [Scientific Background](#scientific-background)
3. [Architecture](#architecture)
4. [Installation & Setup](#installation--setup)
5. [CLI Usage](#cli-usage)
6. [Python API Usage](#python-api-usage)
7. [Input Requirements](#input-requirements)
8. [Output Files](#output-files)
9. [The Intent System](#the-intent-system)
10. [Enhanced Analyses](#enhanced-analyses)
11. [Complete Examples](#complete-examples)
12. [Validation Workflow](#validation-workflow)
13. [Method Selection Guide](#method-selection-guide)
14. [Parameter Tuning](#parameter-tuning)
15. [Troubleshooting](#troubleshooting)
16. [Integration with Other Modules](#integration-with-other-modules)
17. [Scientific References](#scientific-references)

---

## Overview

Module 10 provides end-to-end biomarker discovery from RNA-seq data. Starting
from a count matrix and sample metadata, it discovers a gene panel that
discriminates between conditions, evaluates its clinical utility, and produces
a publication-ready report.

### Why Use Module 10?

**Problem:** Biomarker discovery requires expertise across multiple fields:
- Feature selection (machine learning)
- Classification and cross-validation (statistics)
- Panel size optimization (model selection)
- Clinical metrics (Bayes' theorem, decision theory)
- Biological annotation (gene ontology, pathways)

Researchers typically cobble together 5-10 separate tools, each with different
input formats, programming languages, and statistical assumptions.

**Solution:** Module 10 provides a unified pipeline where:
- One function call runs the complete workflow
- Eight feature selection methods with consensus ranking
- Four classifiers with nested cross-validation and feature selection inside the outer fold (Ambroise & McLachlan 2002 leakage fix)
- Automatic panel size detection (kneedle algorithm) with consensus pinning across CV folds
- Clinical utility metrics (PPV/NPV at prevalence, Youden's J, decision curves, bootstrap CI)
- Panel stability index (Nogueira 2018 Phi) reported with bootstrap CI and benchmark labels
- Self-normalizing ratio biomarkers
- Direction pattern validation for translational work
- Tiered small-cohort advisory and multi-seed verification
- Publication-ready Markdown report with biological annotation

### Key Features

- **Multi-method feature selection**: Univariate filter, Elastic Net, LASSO, Boruta, mRMR, RFE, SHAP, WGCNA
- **Rigorous validation**: Nested cross-validation with feature selection inside the outer fold (Ambroise & McLachlan 2002 PNAS), bootstrap confidence intervals, optional LOOCV
- **Panel-size auto-detection**: Kneedle algorithm (Satopaa et al. 2011 ICDCSW) with consensus pinning across CV folds
- **Panel stability**: Nogueira (2018 JMLR) stability index Phi computed across per-fold panels, with bootstrap CI and benchmark labels
- **Optimism diagnostic**: Cross-validated AUC vs training-set AUC with the gap surfaced as an overfitting signal
- **Clinical utility metrics**: Out-of-fold PPV/NPV at user-supplied prevalence, Youden's optimal threshold, sensitivity/specificity, decision curve analysis (Vickers & Elkin 2006), net reclassification improvement
- **Direction patterns**: UP/DOWN gene patterns with cross-cohort concordance checking
- **Ratio biomarkers**: Self-normalizing gene-pair ratios robust to batch effects (Geman et al. 2004 TSP)
- **Biological annotation**: MyGene.info gene info, KEGG/Reactome/GO pathway enrichment, Europe PMC literature mining, STRING protein-protein interactions
- **Small-cohort tooling**: Tiered advisory at n<20/<50, multi-seed verification button to surface single-seed instability
- **Intent system**: Six study-design modes (diagnostic, prognostic, predictive, monitoring, exploratory, translational) auto-configure metrics, validation strategy, and required metadata columns

### Pipeline Stages

The pipeline is organized so that every step which uses the outcome
labels is restricted to the training data of each outer cross-validation
fold. This is the **Ambroise & McLachlan (2002) fix**: when feature
selection sees the labels of samples that will later be evaluated, the
cross-validated AUC is optimistically biased. The same correction is
applied in the `nestedcv` R package (Lewis et al. 2023, Bioinformatics
Advances), Varma & Simon (2006 BMC Bioinformatics 7:91), and Hastie
Tibshirani & Friedman (Elements of Statistical Learning, section
7.10.2).

```
Input: Count matrix + Metadata
    |
    v
Outer cross-validation loop (5 folds x 1 repeat by default)
    For each outer fold:
        |
        +-- Stage 1: Multi-method feature selection on training data only
        |       (Welch's t-test filter, Elastic Net, LASSO, Boruta,
        |        mRMR, RFE, SHAP, WGCNA -- whichever are requested
        |        and have their dependencies installed)
        |
        +-- Stage 2: Consensus gene ranking on training data only
        |       Average rank across methods + selection frequency
        |
        +-- Stage 3: Panel-size auto-detection (kneedle algorithm)
        |       Either per-fold detection, or consensus pinning across
        |       all folds (default: consensus). See "Panel-Size
        |       Auto-Detection" section below for details.
        |
        +-- Stage 4: Classifier fit on training data, predict held-out fold
                Four classifiers: Logistic Regression, Random Forest,
                SVM, XGBoost (with deterministic tiebreak when AUCs
                coincide). Held-out predictions accumulate into honest
                out-of-fold (OOF) arrays.
    |
    v
Final-panel pass: Stages 1-3 re-run on ALL data
    Produces the deployable single panel that users see as "the panel".
    The OOF AUC from the loop above estimates how well the procedure
    generalizes; the final panel is what you would actually use.
    |
    v
Stage 5: Panel stability diagnostic (Nogueira 2018 Phi)
    Quantifies how much the per-fold panels agree, with bootstrap
    95% CI and benchmark labels (excellent/intermediate/poor).
    |
    v
Stage 6: Biological annotation
    MyGene.info gene info, KEGG/Reactome/GO pathway enrichment,
    Europe PMC literature mining, STRING PPI network query.
    |
    v
Stage 7 (optional, when intent is set): Enhanced analyses
    +-- Signature score      (weighted risk score per patient)
    +-- Direction pattern    (UP/DOWN gene direction validation)
    +-- Clinical metrics     (Youden, OOF PPV/NPV, DCA, bootstrap CI)
    +-- Ratio biomarkers     (self-normalizing gene-pair search)
    |
    v
Output: BiomarkerResult / EnhancedBiomarkerResult
```

---

## Scientific Background

### Feature Selection Methods

Module 10 uses multiple feature selection approaches because no single method
is universally best. Each method has different strengths:

| Method | Type | Approach | Strength | Reference |
|--------|------|----------|----------|-----------|
| **Univariate filter** | Filter | Per-fold Welch's t-test or Mann-Whitney U | Fold-safe; default replacement for DE filter | Haury et al. 2011 PLOS ONE |
| **DE filter** | Filter | Significant genes from M7/M8/M9 | Bridges upstream DE analysis (warns about upstream leakage) | — |
| **Elastic Net** | Embedded | L1+L2 penalized logistic regression | Handles correlated genes | Zou & Hastie, 2005 |
| **LASSO** | Embedded | L1 penalized logistic regression | Sparse solutions | Tibshirani, 1996 |
| **Boruta** | Wrapper | Random Forest shadow features | Identifies all relevant features | Kursa & Rudnicki, 2010 |
| **mRMR** | Filter | Minimum Redundancy Maximum Relevance | Selects non-redundant features | Ding & Peng, 2005 |
| **RFE** | Wrapper | Recursive Feature Elimination | Iterative refinement | Guyon et al., 2002 |
| **SHAP** | Model-agnostic | Shapley value importance | Interpretable rankings | Lundberg & Lee, 2017 |
| **WGCNA** | Network | Co-expression hub genes | Captures gene modules | Langfelder & Horvath, 2008 |

The default method set is `univariate_filter`, `elastic_net`, `rfe`, with `boruta`,
`mrmr`, and `shap` added when their optional dependencies are installed.

### Leakage Fix (Ambroise & McLachlan 2002 PNAS)

In a naive nested cross-validation pipeline, feature selection runs
once on the full dataset and the selected genes are then evaluated by
CV. Ambroise & McLachlan (2002 PNAS 99:6562) showed this gives
optimistically biased CV AUC estimates: by the time the validation
fold is evaluated, the selected genes have already "seen" the
validation labels indirectly through the selection step.

The fix is to push every step that uses the outcome labels inside the
training-side of the outer fold:

- Feature selection: training fold only
- Consensus ranking: training fold only
- Panel-size detection: training fold only (or pinned via consensus, see below)
- Classifier fit: training fold only
- Prediction: held-out fold

The held-out predictions across all folds form an honest out-of-fold
(OOF) array. AUC, Phi, and clinical metrics computed from OOF
predictions reflect generalization performance, not in-sample fit.

A separate "final-panel pass" re-runs feature selection and panel
optimization on **all** data to produce the single panel that gets
deployed. This is the convention used by the `nestedcv` R package
(Lewis et al. 2023): the OOF AUC tells you how well the procedure
generalizes, the final panel tells you what to deploy.

The `de_filter` method is kept as an option for users who supply
externally-computed DE genes, but it triggers a warning: any DE list
computed on the same cohort that will later be CV-evaluated has
already seen the labels and carries upstream selection bias that
M10's pipeline-CV cannot correct. The fold-safe replacement is
`univariate_filter`, which recomputes per-fold on training data only.

### Consensus Ranking

After running multiple methods, Module 10 aggregates results using rank
aggregation. Each gene receives a consensus score combining:
- **Normalized mean rank** across methods (70% weight)
- **Selection frequency** -- how many methods selected it (30% weight)

This produces a single ranked gene list where the top genes are those
consistently identified across multiple statistical frameworks.

### Significance Calibration of the Consensus Score

The raw consensus score rewards genes that score well across methods,
but it does not enforce a significance floor. A gene that scores 4th
overall but has a univariate p-value of 0.7 contributes noise to the
panel. To address this, the consensus ranking is post-processed by a
two-tier multiplicative weighting:

```
weight = 1.0  if p_value < alpha
weight = 0.5  otherwise
consensus_score_calibrated = consensus_score * weight
```

The default alpha is 0.05. Genes that fail the significance floor are
not removed, but they are demoted relative to genes that pass; forward
selection then draws its top candidates from the calibrated ordering.
The original `consensus_score` column is preserved for auditability.
This design borrows the spirit of the Robust Rank Aggregation framework
(Kolde et al. 2012 Bioinformatics) and Empirical-Bayes shrinkage
(Efron 2004 JASA), with one fixed shrinkage factor instead of a fitted
mixture model.

### Panel-Size Auto-Detection (Kneedle)

After consensus ranking, the panel is built by walking down the ranked
gene list and evaluating CV AUC at each candidate panel size. The
panel-size-vs-AUC curve usually has a knee: small panels miss signal,
large panels add noise, and somewhere in between is a natural
inflection point.

Module 10 detects this knee using the **kneedle algorithm** (Satopaa
et al. 2011 ICDCSW, "Finding a Kneedle in a Haystack"). Kneedle
detects the maximum-curvature point on a normalized, polynomially
smoothed version of the curve. Three failure modes of the legacy
"first-drop" heuristic are handled correctly by kneedle:

- **Non-monotone curves**: a single noisy down-tick early in the
  curve no longer terminates the walk prematurely.
- **Saturated curves**: when AUC saturates near 1.0 from the very
  first panel size, kneedle recognizes there is no inflection and
  falls back to argmax (smallest panel size at maximum AUC).
- **Late knees**: real knees past panel size 30 are no longer
  missed by an early sub-threshold step.

Three strategies are exposed via `auto_panel_strategy`:

- `kneedle` (default): Satopaa et al. 2011 with polynomial smoothing
  and argmax fallback for saturated curves.
- `argmax`: smallest panel size at maximum CV AUC.
- `first_drop`: legacy pre-v2.2.2 heuristic, retained for backward
  compatibility.

### Consensus Pinning Across CV Folds

When panel-size detection runs independently inside each outer fold
(`panel_size_strategy='per_fold'`), small fluctuations in the
training-fold curve can cause some folds to land on a much larger or
smaller panel size than others. This destabilizes the panel stability
estimate downstream: a single fold whose argmax tiebreak picks a large
size can drop Phi by 0.3 or more even when the discovery-level K is
stable.

The `consensus` strategy (default) addresses this by running a
discovery-level scout pass on the full data before the CV loop. The
scout's auto-detected K is then pinned across all outer folds.
Per-fold detection is still available as an opt-in for users who
specifically want it. The change to consensus-as-default was driven
by sanity-check evidence on real fixtures: on a 60-sample cohort
under a single seed, per-fold detection produced fold sizes
(3, 3, 3, 12, 3) and dropped Phi by 0.313 vs consensus mode without
changing the discovery-level K.

### Panel Stability (Nogueira 2018 Phi)

Pipeline-CV produces one gene panel per outer fold, by construction.
The degree to which these per-fold panels agree is itself a signal:
stable panels indicate a robust biomarker set; unstable panels
indicate the chosen genes are sensitive to small data perturbations.

Module 10 quantifies stability using the measure of Nogueira, Sechidis
& Brown (2018, JMLR 18:174). Nogueira Phi is the only similarity-based
stability measure that satisfies all five of their axiomatic
properties: fully defined, strictly monotonic, bounded, maximum
stability iff deterministic selection, and **corrected for chance**.
Jaccard, Dice, and POG all fail correction-for-chance and reward
larger feature sets artificially.

The benchmark scale (used by the `stabm` R package and Nogueira's
original paper):

- `Phi >= 0.75` -- excellent agreement beyond chance
- `0.4 <= Phi < 0.75` -- intermediate to good
- `Phi < 0.4` -- poor agreement

Phi is reported with a bootstrap 95% confidence interval over the fold
axis. With small M (e.g. 5 folds), the CI is intentionally wide and
honestly reflects the small-sample uncertainty of the estimate.

### Optimism Diagnostic

A trained classifier almost always achieves higher AUC on its training
data than on held-out data. The size of the gap tells you whether the
model is overfitting:

- **Small gap (e.g. < 0.05)**: model generalizes well; CV AUC is a
  reasonable estimate of deployment performance.
- **Large gap (e.g. > 0.15)**: model memorized the training data;
  deployment AUC will be closer to the CV estimate than the training
  estimate.

Module 10 reports both:

- **Training AUC**: the AUC the final-panel classifier achieves on the
  data it was fit on. Always optimistic.
- **CV AUC** (with bootstrap 95% CI): the AUC averaged across outer
  folds. Honest estimate of deployment performance.
- **Gap = Training AUC - CV AUC**: if this is large, the dashboard
  surfaces an amber/red banner.

### Clinical Utility Metrics

Beyond classification accuracy, clinical biomarkers require:

| Metric | Question It Answers | Reference |
|--------|---------------------|-----------|
| **PPV/NPV at prevalence** | If a patient tests positive, what's the probability they actually have the disease? | Altman & Bland, BMJ 1994 |
| **Youden's J** | At what threshold should we split positive vs negative? | Youden, 1950 |
| **Bootstrap CI** | How uncertain is our AUC estimate? | Efron, 1979 |
| **Decision Curve Analysis** | Does using this biomarker improve clinical decisions? | Vickers & Elkin, 2006 |
| **Net Reclassification** | Does the new biomarker move patients to correct risk categories? | Pencina et al., 2008 |

### Ratio Biomarkers

Gene ratios (geneA/geneB) are powerful biomarkers because they are
**self-normalizing**: batch effects, library size differences, and platform
variation cancel when dividing one gene by another measured in the same sample.
Based on the Top Scoring Pair (TSP) method (Geman et al., 2004).

### Direction Patterns

A biomarker signature isn't just "these genes matter" — it's "these genes go
UP and these go DOWN" in disease. Direction patterns capture per-gene
direction, fold change, and confidence, enabling cross-cohort validation
and translational checks.

### Small-Cohort Tooling

At small sample sizes (n < 50, and especially n < 20), single-seed CV AUC
estimates can be substantially unstable: re-running discovery with a
different random seed can produce a different gene panel, a different
classifier, and a noticeably different CV AUC. This is not a bug -- it's
the honest variance of nested cross-validation at small n. Reporting a
single-seed result without flagging this gives a misleading picture of
how reliable the panel is.

Module 10 surfaces this in two places:

- **Tiered advisory at the top of the Clinical Metrics tab**: a
  strong advisory below n=20, a moderate advisory below n=50, and no
  advisory at n>=50. The advisory text describes which downstream
  metrics carry the most run-to-run variance at the current n.

- **Multi-seed verification button on the Panel tab**: re-runs
  discovery at additional random seeds and produces a side-by-side
  truth table showing which diagnostics are stable across seeds (the
  panel banner color, the Phi label, the top-N genes) and which vary
  (the exact CV AUC, the best classifier, occasional gene swaps in
  positions 4-5 of a 5-gene panel). Users can trust the stable
  diagnostics and treat the variable ones with caution.

The multi-seed view is especially relevant for paper figures: a
single-seed CV AUC reported as a point estimate without a confidence
interval misleads readers about the precision of the result. The
verification button makes the run-to-run variance visible, and the
bootstrap 95% CI on CV AUC (always reported) makes the within-seed
sampling variance visible.

---

## Architecture

### Subpackage Structure

```
raptor/biomarker_discovery/
├── __init__.py            — Public API re-exports
├── core.py                — Main pipeline: discover_biomarkers()
│   ├── FeatureSelector        — 8 feature selection methods
│   ├── ClassifierEvaluator    — 4 classifiers with nested CV
│   ├── PanelOptimizer         — Forward selection + stability
│   ├── SurvivalAnalyzer       — Cox regression (optional)
│   ├── BiologicalAnnotator    — MyGene, pathways, PPI, literature
│   └── BiomarkerResult        — Main result container
├── intent.py              — BiomarkerIntent (6 study design modes)
├── signature_score.py     — SignatureScore + build_signature_score()
├── direction_patterns.py  — DirectionPattern + build_direction_pattern()
├── clinical_metrics.py    — 5 clinical utility functions
├── ratio_biomarkers.py    — RatioBiomarkerSearcher + apply_ratios()
└── enhanced.py            — EnhancedBiomarkerResult + integration
```

### Data Flow

```
discover_biomarkers(counts, metadata, intent="diagnostic")
    │
    ├─→ _prepare_expression_data()     → X (samples×genes), y (0/1 labels)
    ├─→ FeatureSelector                → ranked gene list
    ├─→ PanelOptimizer                 → optimal gene panel
    ├─→ ClassifierEvaluator            → AUC, F1, sensitivity, specificity
    ├─→ BiologicalAnnotator            → gene info, pathways, literature
    │
    └─→ enhance_biomarker_result()     [if intent is set]
        ├─→ build_signature_score()    → per-patient risk score
        ├─→ build_direction_pattern()  → UP/DOWN pattern
        ├─→ clinical metrics           → Youden, CI, PPV/NPV, DCA
        └─→ RatioBiomarkerSearcher     → top ratio pairs
```

---

## Installation & Setup

### Basic Installation

```bash
pip install raptor-rnaseq
```

### Verify Module 10

```python
from raptor.biomarker_discovery import discover_biomarkers, get_dependencies_status

# Check optional dependencies
deps = get_dependencies_status()
for name, available in deps.items():
    status = "available" if available else "missing"
    print(f"  {status} {name}")
```

### Optional Dependencies

| Package | Required For | Install |
|---------|-------------|---------|
| scikit-learn | All of M10 (core requirement) | `pip install scikit-learn` |
| boruta | Boruta feature selection | `pip install boruta` |
| mrmr-selection | mRMR feature selection | `pip install mrmr-selection` |
| shap | SHAP-based feature ranking | `pip install shap` |
| xgboost | XGBoost classifier + SHAP | `pip install xgboost` |
| lifelines | Survival analysis (prognostic) | `pip install lifelines` |
| PyWGCNA | WGCNA hub gene selection | `pip install PyWGCNA` |
| gseapy | Pathway enrichment (KEGG, GO) | `pip install gseapy` |
| mygene | Gene annotation | `pip install mygene` |

### Install All Optional Dependencies

```bash
pip install raptor-rnaseq[biomarker]
# or individually:
pip install boruta mrmr-selection shap xgboost lifelines gseapy mygene
```

---

## CLI Usage

### Basic Biomarker Discovery

```bash
# Minimal (uses defaults: elastic_net + RFE, auto panel size)
raptor biomarker -c counts.csv -m metadata.csv -g condition

# With upstream DE results
raptor biomarker -c counts.csv -m metadata.csv -g condition \
    -d de_result.pkl

# With M9 ensemble result
raptor biomarker -c counts.csv -m metadata.csv -g condition \
    -e ensemble_result.pkl --panel-size 10

# With intent (enables enhanced analyses)
raptor biomarker -c counts.csv -m metadata.csv -g condition \
    --intent diagnostic --prevalence 0.05

# Fast run (skip annotation and network queries)
raptor biomarker -c counts.csv -m metadata.csv -g condition \
    --no-annotate
```

### CLI Options

| Option | Flag | Default | Description |
|--------|------|---------|-------------|
| **Counts** | `-c` | Required | Count matrix CSV (genes x samples) |
| **Metadata** | `-m` | Required | Sample metadata CSV |
| **Group column** | `-g` | `condition` | Column defining groups |
| **DE result** | `-d` | None | Module 7 DE result pickle |
| **Ensemble result** | `-e` | None | Module 9 ensemble pickle |
| **Methods** | `--methods` | Auto | Feature selection methods |
| **Panel size** | `--panel-size` | Auto | Target panel size; overrides auto-detection |
| **Auto panel strategy** | `--auto-panel-strategy` | `kneedle` | `kneedle`, `argmax`, or `first_drop` (legacy). Used when `--panel-size` is not specified. |
| **Panel sensitivity** | `--panel-sensitivity` | `1.0` | Kneedle's S parameter; lower picks earlier knees, higher picks later. Ignored unless `--auto-panel-strategy=kneedle`. |
| **Panel size strategy** | `--panel-size-strategy` | `consensus` | `consensus` (detect once on full data, pin all CV folds to that K) or `per_fold` (each fold detects independently). |
| **Species** | `--species` | `human` | Species for annotation |
| **Disease term** | `--disease-term` | None | PubMed search context |
| **Intent** | `--intent` | None | `diagnostic` or `exploratory` |
| **Prevalence** | `--prevalence` | 0.05 | Disease prevalence for PPV/NPV |
| **No annotate** | `--no-annotate` | False | Skip annotation |
| **No literature** | `--no-literature` | False | Skip PubMed search |
| **No PPI** | `--no-ppi` | False | Skip STRING query |
| **Output** | `-o` | `results/biomarkers` | Output directory |

### Survival Biomarkers

```bash
raptor biomarker-survival \
    -c counts.csv \
    --clinical clinical.csv \
    --time-column os_time \
    --event-column os_event \
    -o results/survival_biomarkers
```

### Validate on Independent Cohort

```bash
raptor biomarker-validate \
    --panel results/biomarkers/biomarker_panel.csv \
    --counts validation_counts.csv \
    --metadata validation_metadata.csv \
    -o results/validation
```

---

## Python API Usage

### Quick Start

```python
from raptor.biomarker_discovery import discover_biomarkers

# Basic discovery
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    group_column="condition",
)

print(f"Panel: {result.panel}")
print(f"AUC: {result.classification_results[result.best_classifier].auc:.3f}")
```

### With Intent (Enhanced Analyses)

```python
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    group_column="condition",
    intent="diagnostic",
    prevalence=0.05,
    output_dir="results/diagnostic",
)

# Base results (always available)
print(f"Panel: {result.panel}")
print(f"Best classifier: {result.best_classifier}")

# Enhanced results (available because intent was set)
if result.signature:
    print(f"Signature: {result.signature.mode}, "
          f"cutoffs={result.signature.cutoffs}")

if result.direction_pattern:
    dp = result.direction_pattern
    print(f"Pattern: {dp.n_up} UP, {dp.n_down} DOWN")

if result.clinical_metrics:
    cm = result.clinical_metrics
    print(f"Youden's J: {cm['youdens']['youdens_j']:.3f}")
    print(f"AUC 95% CI: [{cm['bootstrap_ci']['ci_lower']:.3f}, "
          f"{cm['bootstrap_ci']['ci_upper']:.3f}]")
    print(f"PPV at 5%: {cm['ppv_npv']['ppv']:.3f}")

if result.ratio_result and result.ratio_result.best_pair:
    bp = result.ratio_result.best_pair
    print(f"Best ratio: {bp.name} (AUC={bp.auc:.3f})")
```

### Using Individual Modules Standalone

Each enhanced module can be used independently without running the full
discovery pipeline.

#### Signature Score

```python
from raptor.biomarker_discovery import build_signature_score
import pandas as pd
import numpy as np

# Your expression data (samples × genes)
X = pd.read_csv("expression.csv", index_col=0)
y = np.array([0]*30 + [1]*30)  # binary labels
panel = ["BRCA1", "TP53", "MYC", "EGFR", "HER2"]

# Build signature
sig = build_signature_score(X=X, y=y, panel_genes=panel, mode="diagnostic")

# Score new patients
new_patient_scores = sig.score(new_expression_data)
print(f"Patient scores: {new_patient_scores}")

# Stratify into risk groups
risk_groups = sig.stratify(new_patient_scores)
print(f"Risk groups: {risk_groups.value_counts()}")
```

#### Direction Patterns

```python
from raptor.biomarker_discovery import build_direction_pattern

# Learn pattern from discovery cohort
pattern = build_direction_pattern(
    expression=expression_df,
    labels=pd.Series(["healthy"]*30 + ["disease"]*30),
    reference_group="disease",
    baseline_group="healthy",
    p_threshold=0.05,
    fc_threshold=1.0,
)

print(f"Pattern: {pattern.n_up} UP, {pattern.n_down} DOWN")
print(pattern.to_dataframe())

# Score a new patient
score = pattern.concordance(new_patient_expression)
print(f"Concordance: {score:.3f}")
# +1 = fully matches disease, -1 = opposite, 0 = baseline

# Compare with another cohort's pattern
report = pattern.cross_cohort_check(other_pattern)
print(f"Agreement: {report['agreement_fraction']:.1%}")
print(f"Binomial p: {report['binomial_p']:.4f}")
```

#### Clinical Metrics

```python
from raptor.biomarker_discovery import (
    ppv_npv_at_prevalence,
    bootstrap_ci,
    youdens_optimal_threshold,
    decision_curve_analysis,
)

# Your classifier predictions
y_true = np.array([0]*50 + [1]*50)
y_score = model.predict_proba(X_test)[:, 1]

# 1. Optimal threshold
youden = youdens_optimal_threshold(y_true, y_score)
print(f"Optimal threshold: {youden['threshold']:.3f}")
print(f"Sensitivity: {youden['sensitivity']:.3f}")
print(f"Specificity: {youden['specificity']:.3f}")
print(f"Youden's J: {youden['youdens_j']:.3f}")

# 2. Bootstrap confidence interval
ci = bootstrap_ci(y_true, y_score, n_bootstrap=2000)
print(f"AUC: {ci['point_estimate']:.3f} "
      f"[{ci['ci_lower']:.3f}, {ci['ci_upper']:.3f}]")

# 3. PPV/NPV at different prevalences
for prev in [0.01, 0.05, 0.10]:
    result = ppv_npv_at_prevalence(
        sensitivity=youden['sensitivity'],
        specificity=youden['specificity'],
        prevalence=prev,
    )
    print(f"Prevalence {prev:.0%}: "
          f"PPV={result['ppv']:.3f}, NPV={result['npv']:.3f}")

# 4. Decision curve analysis
dca = decision_curve_analysis(y_true, y_score)
# dca is a DataFrame with columns:
# threshold, net_benefit_model, net_benefit_treat_all, net_benefit_treat_none
```

#### Ratio Biomarkers

```python
from raptor.biomarker_discovery import (
    RatioBiomarkerSearcher,
    apply_ratios,
    build_ratio_features,
)

# Search for discriminating gene ratios
searcher = RatioBiomarkerSearcher(top_k=10, min_auc=0.80)
result = searcher.search(
    expression=expression_df,
    labels=pd.Series(["healthy"]*30 + ["disease"]*30),
    reference_group="disease",
    baseline_group="healthy",
)

print(f"Found {result.n_found} ratio pairs")
for pair in result.pairs[:5]:
    print(f"  {pair.name}: AUC={pair.auc:.3f} ({pair.direction})")

# Apply to new data
ratio_features = apply_ratios(new_expression, result.pairs)
# ratio_features is a DataFrame with one column per ratio pair

# Or use the convenience wrapper
result, ratio_df = build_ratio_features(
    expression_df, labels, "disease", "healthy",
    top_k=10, min_auc=0.80,
)
```

---

## Input Requirements

### Count Matrix (`counts.csv`)

```csv
gene_id,sample_1,sample_2,sample_3,sample_4,sample_5,sample_6
ENSG00000141510,1523,1872,2103,456,523,678
ENSG00000146648,456,523,678,1523,1872,2103
ENSG00000012048,789,834,901,234,267,312
```

- **Rows** = genes (gene IDs as row names/first column)
- **Columns** = samples
- **Values** = raw integer counts (not normalized, not log-transformed)
- Module 10 internally applies log2-CPM normalization

### Metadata (`metadata.csv`)

```csv
sample_id,condition,batch
sample_1,healthy,1
sample_2,healthy,1
sample_3,healthy,2
sample_4,disease,1
sample_5,disease,2
sample_6,disease,2
```

- Must have a column matching sample IDs in count matrix
- Must have a group column (specified via `--group-column` or `group_column=`)
- Group column must contain exactly 2 unique values for binary classification
- Additional columns (batch, age, sex) are preserved but not used by M10

### Optional: DE Result Pickle

```python
# From Module 7
from raptor import import_deseq2
de_result = import_deseq2("deseq2_results.csv")
de_result.save("de_result.pkl")

# Pass to Module 10
result = discover_biomarkers(counts, metadata, de_result=de_result)
```

### Optional: Ensemble Result Pickle

```python
# From Module 9
from raptor import ensemble_brown
ens_result = ensemble_brown({"deseq2": d, "edger": e, "limma": l})
ens_result.save("ensemble_result.pkl")

# Pass to Module 10
result = discover_biomarkers(counts, metadata, ensemble_result=ens_result)
```

---

## Output Files

### Standard Output

```
results/biomarkers/
├── biomarker_panel.csv              # Final gene panel with ranks
├── ranked_genes.csv                 # Full consensus ranking (all genes)
├── classification_performance.csv   # AUC, F1, sensitivity, specificity per classifier
├── panel_curve.csv                  # Panel size vs AUC curve
├── biomarker_report.md              # Publication-ready Markdown report
├── biomarker_result.pkl             # Complete result object (for downstream)
├── biomarker_params.json            # Parameters used
├── summary.txt                      # Human-readable summary
└── annotations/                     # Biological annotation outputs
    ├── gene_annotations.csv         # Gene symbols, names, GO terms
    ├── pathway_enrichment.csv       # Enriched pathways (KEGG, Reactome, GO)
    ├── literature_hits.csv          # PubMed publications per gene
    └── ppi_network.json             # STRING protein interactions
```

### Enhanced Output (when `intent` is set)

```
results/biomarkers/enhanced/
├── direction_pattern.csv            # Per-gene direction, fold change, confidence
├── ratio_biomarkers.csv             # Top ratio pairs with AUC
├── signature_score.json             # Learned weights, cutoffs, performance
└── decision_curve.csv               # Net benefit across thresholds
```

### File Descriptions

| File | Format | Contents |
|------|--------|----------|
| `biomarker_panel.csv` | CSV | gene_id, rank — the final recommended panel |
| `ranked_genes.csv` | CSV | All genes with consensus_rank, consensus_score, n_methods_selected, per-method ranks |
| `classification_performance.csv` | CSV | classifier, accuracy, auc, sensitivity, specificity, f1 |
| `panel_curve.csv` | CSV | panel_size, auc_mean, auc_std — for plotting the optimization curve |
| `biomarker_report.md` | Markdown | Complete publication-ready report with tables |
| `direction_pattern.csv` | CSV | direction (UP/DOWN), log2FC, neg_log10_p, baseline_mean, baseline_std |
| `ratio_biomarkers.csv` | CSV | gene_a, gene_b, name, auc, direction, mean_ratio per group |
| `signature_score.json` | JSON | weights, intercept, normalization, cutoffs, risk_labels, performance |

---

## The Intent System

### What Is an Intent?

An intent describes the **goal** of your biomarker discovery analysis. Different
goals require different statistical approaches, validation strategies, and
output metrics.

Setting `intent` in `discover_biomarkers()` auto-configures the enhanced
analyses and labels your output appropriately.

### Available Intents

| Intent | Study Design | Question | CLI Available |
|--------|-------------|----------|---------------|
| `diagnostic` | Case-control | Does the patient have the disease? | Available |
| `exploratory` | Unsupervised | What patterns exist in the data? | Available |
| `prognostic` | Cohort/survival | How will the patient fare over time? | Planned for v2.3.0 |
| `predictive` | Treatment interaction | Will the patient respond to treatment? | Planned |
| `monitoring` | Longitudinal | Is the disease progressing? | Planned |
| `translational` | Cross-species | Do mouse findings replicate in humans? | Planned |

### Intent Configuration Details

```python
from raptor.biomarker_discovery import BiomarkerIntent

# Each intent auto-configures these attributes:
intent = BiomarkerIntent("diagnostic")
print(intent.study_design)          # "case_control"
print(intent.validation_strategy)   # "nested_cv"
print(intent.required_data_columns) # ["expression", "group"]
print(intent.minimum_samples)       # 30
print(intent.required_metrics)      # ["AUC", "sensitivity", "specificity", "PPV", "NPV"]
print(intent.output_label)          # "Diagnostic Biomarker Candidates"
```

### Which Enhanced Analyses Run Per Intent

| Intent | Signature Score | Direction Pattern | Clinical Metrics | Ratio Biomarkers |
|--------|:-:|:-:|:-:|:-:|
| diagnostic | Yes | Yes | Yes | Yes |
| exploratory | — | Yes | — | Yes |
| prognostic | Yes (falls back to diagnostic) | Yes | Yes | Yes |
| predictive | — | Yes | Yes | Yes |
| monitoring | — | — | Yes | — |
| translational | — | Yes | — | Yes |

---

## Enhanced Analyses

### 1. Signature Score

Turns a gene panel into a single per-patient risk score.

**How it works:**
1. Fits logistic regression on the panel genes (diagnostic mode)
2. Learns per-gene weights (positive = risk-increasing, negative = risk-decreasing)
3. Finds optimal cutoff via Youden's J statistic
4. New patients: z-score normalize → weighted sum → compare to cutoff

**Output:** `SignatureScore` object with:
- `weights`: per-gene coefficients
- `cutoffs`: risk group boundaries
- `risk_labels`: ["low", "high"] or ["low", "medium", "high"]
- `score()`: apply to new patient data
- `stratify()`: assign risk groups

### 2. Direction Pattern

Captures which genes go UP and which go DOWN in disease.

**How it works:**
1. For each gene: compute mean(disease) - mean(healthy) = log2FC
2. Statistical test (t-test or Mann-Whitney) for significance
3. Filter by p-value and fold-change thresholds
4. Store direction, magnitude, and confidence per gene

**Key methods:**
- `concordance(patient)`: score in [-1, +1] how well a patient matches the pattern
- `cross_cohort_check(other_pattern)`: compare direction agreement between cohorts
- `to_dataframe()`: inspect the full pattern

### 3. Clinical Metrics

Five functions answering "how useful is this biomarker in the clinic?"

| Function | What It Computes |
|----------|-----------------|
| `ppv_npv_at_prevalence()` | PPV and NPV adjusted for real-world disease prevalence |
| `bootstrap_ci()` | Non-parametric confidence intervals for AUC or any metric |
| `youdens_optimal_threshold()` | The cutoff maximizing sensitivity + specificity - 1 |
| `decision_curve_analysis()` | Net benefit compared to treat-all and treat-none |
| `net_reclassification_improvement()` | How many patients move to correct risk category |

### 4. Ratio Biomarkers

Gene ratios (geneA/geneB) that are robust to batch effects.

**How it works:**
1. Test all pairwise gene ratios among panel genes
2. Compute AUC for each ratio as a single-feature classifier
3. Rank by AUC, return top-k pairs above min_auc threshold
4. A pseudocount (default 1.0) prevents division by zero

**Why ratios matter:** If Hospital A's expression values are systematically
2× higher than Hospital B's, individual gene biomarkers fail across sites.
But the ratio geneA/geneB cancels the systematic difference because both
genes are equally affected.

---

## Complete Examples

### Example 1: Diagnostic Biomarker Discovery (CLI)

```bash
# Step 1: Run QC
raptor qc -c counts.csv -m metadata.csv -o qc/

# Step 2: Profile data
raptor profile -c counts.csv -m metadata.csv -g condition

# Step 3: Discover biomarkers with diagnostic intent
raptor biomarker -c counts.csv -m metadata.csv -g condition \
    --intent diagnostic \
    --prevalence 0.05 \
    --species human \
    --disease-term "breast cancer" \
    -o results/diagnostic/

# Step 4: Validate on independent cohort
raptor biomarker-validate \
    --panel results/diagnostic/biomarker_panel.csv \
    --counts validation_counts.csv \
    --metadata validation_metadata.csv \
    -o results/validation/
```

### Example 2: Full Python Workflow

```python
from raptor.biomarker_discovery import (
    discover_biomarkers,
    validate_biomarkers,
    build_direction_pattern,
)
import pandas as pd

# 1. Discover biomarkers
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    group_column="condition",
    intent="diagnostic",
    prevalence=0.05,
    species="human",
    disease_term="Alzheimer's disease",
    output_dir="results/alzheimer_biomarkers",
)

# 2. Print summary
print(result.summary())

# 3. Access the gene panel
print(f"\nRecommended panel ({result.panel_size} genes):")
for i, gene in enumerate(result.panel, 1):
    print(f"  {i}. {gene}")

# 4. Classification performance
for name, clf in result.classification_results.items():
    marker = " ← best" if name == result.best_classifier else ""
    print(f"  {name}: AUC={clf.auc:.3f}, F1={clf.f1:.3f}{marker}")

# 5. Enhanced analyses
if result.clinical_metrics:
    cm = result.clinical_metrics
    
    # Youden's optimal threshold
    y = cm['youdens']
    print(f"\nYouden's optimal threshold: {y['threshold']:.3f}")
    print(f"  Sensitivity: {y['sensitivity']:.3f}")
    print(f"  Specificity: {y['specificity']:.3f}")
    print(f"  J statistic: {y['youdens_j']:.3f}")
    
    # Bootstrap CI
    b = cm['bootstrap_ci']
    print(f"\nAUC: {b['point_estimate']:.3f} "
          f"[{b['ci_lower']:.3f}, {b['ci_upper']:.3f}] 95% CI")
    
    # PPV/NPV at prevalence
    p = cm['ppv_npv']
    print(f"\nAt {p['prevalence']:.1%} prevalence:")
    print(f"  PPV: {p['ppv']:.3f}")
    print(f"  NPV: {p['npv']:.3f}")

# 6. Direction pattern
if result.direction_pattern:
    dp = result.direction_pattern
    print(f"\nDirection pattern: {dp.n_up} UP, {dp.n_down} DOWN")
    print(dp.to_dataframe())

# 7. Ratio biomarkers
if result.ratio_result:
    print(f"\nTop ratio biomarkers:")
    for pair in result.ratio_result.pairs[:5]:
        print(f"  {pair.name}: AUC={pair.auc:.3f}")

# 8. Validate on independent data
val_results = validate_biomarkers(
    panel_genes=result.panel,
    counts="validation_counts.csv",
    metadata="validation_metadata.csv",
    group_column="condition",
)
for name, clf in val_results.items():
    print(f"  Validation {name}: AUC={clf.auc:.3f}")
```

### Example 3: Exploratory Analysis

```python
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    group_column="condition",
    intent="exploratory",
    annotate=False,    # Skip annotation for speed
    output_dir="results/exploratory",
)

# Exploratory runs direction patterns + ratio biomarkers
# but skips signature score and clinical metrics
if result.direction_pattern:
    print(f"Direction pattern found: {result.direction_pattern.n_genes} genes")

if result.ratio_result:
    print(f"Ratio pairs: {result.ratio_result.n_found}")
```

### Example 4: Using Upstream Module Results

```python
import pickle
from raptor.biomarker_discovery import discover_biomarkers

# Load M9 ensemble result
with open("results/ensemble/ensemble_result.pkl", "rb") as f:
    ensemble_result = pickle.load(f)

# Discover biomarkers using consensus genes from M9
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    group_column="condition",
    ensemble_result=ensemble_result,  # Uses consensus genes as candidates
    intent="diagnostic",
    target_panel_size=10,             # Aim for exactly 10 genes
    output_dir="results/biomarkers_from_ensemble",
)

print(f"Panel from ensemble: {result.panel}")
```

---

## Validation Workflow

### Discovery → Validation Pipeline

```
Discovery Cohort                   Validation Cohort
    │                                    │
    ├─→ discover_biomarkers()            │
    │       ↓                            │
    │   BiomarkerResult                  │
    │       ↓                            │
    │   result.panel                     │
    │       ↓                            ↓
    └───────────────→ validate_biomarkers(panel, val_counts, val_metadata)
                            ↓
                    Validation AUC, F1, etc.
```

### Cross-Cohort Direction Validation

```python
from raptor.biomarker_discovery import build_direction_pattern

# Learn pattern from cohort 1
pattern1 = build_direction_pattern(
    cohort1_expression, cohort1_labels,
    reference_group="disease", baseline_group="healthy",
)

# Learn pattern from cohort 2
pattern2 = build_direction_pattern(
    cohort2_expression, cohort2_labels,
    reference_group="disease", baseline_group="healthy",
)

# Check agreement
report = pattern1.cross_cohort_check(pattern2)
print(f"Direction agreement: {report['agreement_fraction']:.1%}")
print(f"p-value: {report['binomial_p']:.4f}")

if report['agreement_fraction'] > 0.8 and report['binomial_p'] < 0.05:
    print("Strong direction consistency across cohorts")
else:
    print("Weak direction consistency -- investigate disagreements")
    for gene, dir1, dir2 in report['disagreements']:
        print(f"  {gene}: cohort1={dir1}, cohort2={dir2}")
```

---

## Method Selection Guide

### Feature Selection Methods

| Your Situation | Recommended Methods |
|----------------|-------------------|
| **Small samples (n < 20)** | `elastic_net`, `rfe` |
| **Large samples (n > 50)** | All methods including `boruta`, `mrmr`, `shap` |
| **Have DE results from M7/M9** | Include `de_filter` |
| **Want module-level features** | Include `wgcna` (needs n ≥ 15) |
| **Want interpretability** | Include `shap` |
| **Maximum rigor** | Use all available methods |

### Classifier Selection

Module 10 automatically runs all available classifiers. The best one is
selected by AUC:

| Classifier | Strengths | When It Excels |
|-----------|-----------|----------------|
| Logistic Regression | Interpretable, fast | Linear separability |
| Random Forest | Handles interactions | Complex patterns |
| SVM | Good with small samples | High-dimensional data |
| XGBoost | Highest accuracy | Large datasets |

### Intent Selection

| Your Study | Use This Intent |
|------------|----------------|
| Case-control (disease vs healthy) | `diagnostic` |
| Pattern mining, no specific hypothesis | `exploratory` |
| Survival outcome | `prognostic` (planned) |
| Treatment response prediction | `predictive` (planned) |
| Longitudinal monitoring | `monitoring` (planned) |
| Mouse model → human translation | `translational` (planned) |

---

## Parameter Tuning

### Key Parameters

| Parameter | Default | Guidance |
|-----------|---------|----------|
| `target_panel_size` | Auto (kneedle) | Set explicitly to 5-20 for clinical translation; leave as None to use kneedle auto-detection |
| `auto_panel_strategy` | `kneedle` | Use `argmax` to force smallest-size-at-max-AUC; use `first_drop` only for backward compatibility |
| `panel_sensitivity` | 1.0 | Lower (e.g. 0.5) picks earlier knees, higher (e.g. 2.0) picks later. Most users should leave at default. |
| `panel_size_strategy` | `consensus` | Pin all CV folds to one K (default, recommended). `per_fold` lets each fold detect independently. |
| `min_panel` | 3 | Minimum genes to evaluate |
| `max_panel` | 50 | Maximum genes to evaluate |
| `n_folds` | 5 | Use 10 for paper-grade results (Ambroise & McLachlan recommendation) |
| `n_repeats` | 1 | Set to 3-5 for paper-grade results to reduce CV-variance |
| `validation` | `nested_cv` | LOOCV is deprecated as of v2.2.2 (Ambroise & McLachlan recommend k-fold) |
| `calibrate_consensus` | True | Multiplies consensus_score by significance weight (1.0 if p<alpha, 0.5 otherwise) |
| `alpha` | 0.05 | Significance threshold for the calibration weighting |
| `prevalence` | 0.05 | Set to actual disease prevalence in target population |
| `random_state` | 42 | Change for sensitivity analysis; use the multi-seed verification button on the dashboard for systematic checks |

### When to Override Kneedle Defaults

The `kneedle` + `consensus` defaults are tuned for the typical case
(small to medium cohort, moderate signal strength). Override when:

- You have an **a priori panel size** dictated by a clinical
  constraint (e.g. "panels of 5-10 genes are translatable to
  qPCR, panels >20 are not"). Set `target_panel_size` explicitly.
- The kneedle algorithm picks an unreasonably large panel because
  the curve has a late shallow knee. Try
  `auto_panel_strategy='argmax'` to force the smallest size at
  maximum AUC.
- You specifically want to study fold-to-fold panel-size variance
  itself. Set `panel_size_strategy='per_fold'` and inspect the
  `per_fold_panels` attribute on the result.
- You're comparing v2.2.2 results to a pre-kneedle pipeline. Set
  `auto_panel_strategy='first_drop'` to recover the legacy
  behavior.

### Prevalence Guidelines

| Disease Type | Typical Prevalence |
|-------------|-------------------|
| Rare disease | 0.001 - 0.01 (0.1% - 1%) |
| Common cancer | 0.01 - 0.05 (1% - 5%) |
| Common disease | 0.05 - 0.15 (5% - 15%) |
| Screening population | 0.01 - 0.10 |
| Enriched cohort | 0.20 - 0.50 |

**Impact on PPV:** A biomarker with 90% sensitivity and 90% specificity:
- At 1% prevalence: PPV = 8.3% (most positives are false!)
- At 5% prevalence: PPV = 32.1%
- At 50% prevalence: PPV = 90.0%

---

## Troubleshooting

### Common Issues

#### "Too few genes after filtering"

```
ValueError: Only 5 DE genes in expression data. Using all filtered genes.
```

**Solution:** Your DE gene list has very few genes overlapping with the count
matrix. Check that gene IDs match (Ensembl vs Symbol vs Entrez). Use Module
6b's Gene ID Conversion if needed.

#### "Expected 2 groups for binary classification"

```
ValidationError: Expected 2 groups, got 3: ['control', 'early', 'late']
```

**Solution:** Module 10 currently supports binary classification only.
Combine groups (e.g., merge 'early' and 'late' into 'disease') in your
metadata before running.

#### "No feature selection results available"

**Solution:** All feature selection methods failed. Check that:
- scikit-learn is installed
- Count matrix has enough genes (> 50 recommended)
- Both groups have at least 2 samples

#### "Enhanced analysis failed"

This is non-fatal — the base BiomarkerResult is still returned. Check the
warning message for specifics (usually a missing dependency or data issue).

#### "Annotation stage encountered an error"

Annotation requires internet access (MyGene.info, STRING, PubMed). Use
`--no-annotate` for offline runs.

### Performance Tips

| Tip | Impact |
|-----|--------|
| Use `--no-annotate --no-literature --no-ppi` | 10× faster (skips web queries) |
| Set `target_panel_size` explicitly | Faster panel optimization |
| Use fewer methods (e.g., `--methods elastic_net rfe`) | Faster feature selection |
| Reduce `n_folds` to 3 | Faster CV (less rigorous) |

---

## Integration with Other Modules

### Complete RAPTOR Workflow

```
Module 2: QC → Remove outliers, detect batch effects
    ↓
Module 3: Profile → Extract 32 features, compute BCV
    ↓
Module 4: Recommend → Select DE pipeline (DESeq2/edgeR/limma)
    ↓
[Run DE analysis in R]
    ↓
Module 7: Import → Standardize DE results
    ↓
Module 8: Optimize → Find optimal FDR/LFC thresholds
    ↓
Module 9: Ensemble → Combine multiple DE methods
    ↓
Module 10: Biomarker → Discover and validate gene panel
```

### Passing Results Between Modules

```python
# Module 7 → Module 10
from raptor import import_deseq2
de_result = import_deseq2("deseq2_results.csv")

# Module 9 → Module 10
from raptor import ensemble_brown
ens = ensemble_brown({"deseq2": d, "edger": e})

# Module 10 uses upstream results
result = discover_biomarkers(
    counts="counts.csv",
    metadata="metadata.csv",
    de_result=de_result,           # M7 genes as candidates
    ensemble_result=ens,           # M9 consensus genes (takes priority)
    intent="diagnostic",
)
```

---

## Scientific References

### Feature Selection

1. Tibshirani, R. (1996). Regression shrinkage and selection via the LASSO. *JRSS-B*, 58(1), 267-288.
2. Zou, H. & Hastie, T. (2005). Regularization and variable selection via the elastic net. *JRSS-B*, 67(2), 301-320.
3. Kursa, M.B. & Rudnicki, W.R. (2010). Feature selection with the Boruta package. *Journal of Statistical Software*, 36(11), 1-13.
4. Ding, C. & Peng, H. (2005). Minimum redundancy feature selection from microarray gene expression data. *IEEE TPAMI*, 27(8), 1226-1238.
5. Guyon, I. et al. (2002). Gene selection for cancer classification using support vector machines. *Machine Learning*, 46, 389-422.
6. Lundberg, S.M. & Lee, S.I. (2017). A unified approach to interpreting model predictions. *NeurIPS 2017*.
7. Langfelder, P. & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559.

### Clinical Metrics

8. Youden, W.J. (1950). Index for rating diagnostic tests. *Cancer*, 3(1), 32-35.
9. Efron, B. (1979). Bootstrap methods: another look at the jackknife. *Annals of Statistics*, 7(1), 1-26.
10. Vickers, A.J. & Elkin, E.B. (2006). Decision curve analysis: a novel method for evaluating prediction models. *Medical Decision Making*, 26(6), 565-574.
11. Pencina, M.J. et al. (2008). Evaluating the added predictive ability of a new marker. *Statistics in Medicine*, 27(2), 157-172.
12. Altman, D.G. & Bland, J.M. (1994). Diagnostic tests 2: predictive values. *BMJ*, 309(6947), 102.

### Ratio Biomarkers

13. Geman, D. et al. (2004). Classifying gene expression profiles from pairwise mRNA comparisons. *Statistical Applications in Genetics and Molecular Biology*, 3(1).

### Survival Analysis

14. Cox, D.R. (1972). Regression models and life-tables. *JRSS-B*, 34(2), 187-220.
15. Simon, N. et al. (2011). Regularization paths for Cox's proportional hazards model via coordinate descent. *Journal of Statistical Software*, 39(5).

---

**Version:** 2.2.2  
**Module:** M10 (Biomarker Discovery)  
**Status:** Production Ready (diagnostic + exploratory intents)  
**Tests:** 232 tests (224 pass, 8 skip for optional dependencies)