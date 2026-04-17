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
- ✅ One function call runs the complete workflow
- ✅ Six feature selection methods with consensus ranking
- ✅ Four classifiers with nested cross-validation
- ✅ Automatic panel size optimization
- ✅ Clinical utility metrics (PPV/NPV, decision curves, bootstrap CI)
- ✅ Self-normalizing ratio biomarkers
- ✅ Direction pattern validation for translational work
- ✅ Publication-ready Markdown report with annotation

### Key Features

- 🔬 **Multi-Method Feature Selection**: DE filter, Elastic Net, LASSO, Boruta, mRMR, RFE, SHAP, WGCNA
- 📊 **Rigorous Validation**: Nested cross-validation, LOOCV, bootstrap confidence intervals
- 🎯 **Panel Optimization**: Forward selection with elbow detection
- 🏥 **Clinical Metrics**: PPV/NPV at prevalence, Youden's J, decision curve analysis, NRI
- 🧬 **Direction Patterns**: UP/DOWN gene patterns with cross-cohort concordance checking
- 📈 **Ratio Biomarkers**: Self-normalizing gene ratios robust to batch effects
- 📋 **Biological Annotation**: MyGene.info, pathway enrichment, literature mining, STRING PPI
- 🎯 **Intent System**: Auto-configures analyses based on study design

### Pipeline Stages

```
Input: Count matrix + Metadata
    ↓
Stage 1: Multi-method feature selection (6-8 methods)
    ↓
Stage 2: Consensus gene ranking
    ↓
Stage 3: Panel size optimization (forward selection)
    ↓
Stage 4: Classification evaluation (4 classifiers, nested CV)
    ↓
Stage 5: Biological annotation & reporting
    ↓
Stage 6: Enhanced analyses (if intent is set)
    ├── Signature score (weighted risk scoring)
    ├── Direction pattern (UP/DOWN validation)
    ├── Clinical metrics (PPV, Youden, DCA, bootstrap CI)
    └── Ratio biomarkers (self-normalizing pairs)
    ↓
Output: BiomarkerResult / EnhancedBiomarkerResult
```

---

## Scientific Background

### Feature Selection Methods

Module 10 uses multiple feature selection approaches because no single method
is universally best. Each method has different strengths:

| Method | Type | Approach | Strength | Reference |
|--------|------|----------|----------|-----------|
| **DE Filter** | Filter | Uses significant genes from M7/M8/M9 | Bridges upstream DE analysis | — |
| **Elastic Net** | Embedded | L1+L2 penalized logistic regression | Handles correlated genes | Zou & Hastie, 2005 |
| **LASSO** | Embedded | L1 penalized logistic regression | Sparse solutions | Tibshirani, 1996 |
| **Boruta** | Wrapper | Random Forest shadow features | Identifies all relevant features | Kursa & Rudnicki, 2010 |
| **mRMR** | Filter | Minimum Redundancy Maximum Relevance | Selects non-redundant features | Ding & Peng, 2005 |
| **RFE** | Wrapper | Recursive Feature Elimination | Iterative refinement | Guyon et al., 2002 |
| **SHAP** | Model-agnostic | Shapley value importance | Interpretable rankings | Lundberg & Lee, 2017 |
| **WGCNA** | Network | Co-expression hub genes | Captures gene modules | Langfelder & Horvath, 2008 |

### Consensus Ranking

After running multiple methods, Module 10 aggregates results using rank
aggregation. Each gene receives a consensus score combining:
- **Normalized mean rank** across methods (70% weight)
- **Selection frequency** — how many methods selected it (30% weight)

This produces a single ranked gene list where the top genes are those
consistently identified across multiple statistical frameworks.

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
    status = "✅" if available else "❌"
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
| **Counts** | `-c` | Required | Count matrix CSV (genes × samples) |
| **Metadata** | `-m` | Required | Sample metadata CSV |
| **Group column** | `-g` | `condition` | Column defining groups |
| **DE result** | `-d` | None | Module 7 DE result pickle |
| **Ensemble result** | `-e` | None | Module 9 ensemble pickle |
| **Methods** | `--methods` | Auto | Feature selection methods |
| **Panel size** | `--panel-size` | Auto | Target panel size |
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
| `diagnostic` | Case-control | Does the patient have the disease? | ✅ Yes |
| `exploratory` | Unsupervised | What patterns exist in the data? | ✅ Yes |
| `prognostic` | Cohort/survival | How will the patient fare over time? | 🔜 Planned |
| `predictive` | Treatment interaction | Will the patient respond to treatment? | 🔜 Planned |
| `monitoring` | Longitudinal | Is the disease progressing? | 🔜 Planned |
| `translational` | Cross-species | Do mouse findings replicate in humans? | 🔜 Planned |

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
| diagnostic | ✅ | ✅ | ✅ | ✅ |
| exploratory | — | ✅ | — | ✅ |
| prognostic | ✅ (falls back to diagnostic) | ✅ | ✅ | ✅ |
| predictive | — | ✅ | ✅ | ✅ |
| monitoring | — | — | ✅ | — |
| translational | — | ✅ | — | ✅ |

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
    print("✅ Strong direction consistency across cohorts")
else:
    print("⚠️ Weak direction consistency — investigate disagreements")
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
| `target_panel_size` | Auto (elbow) | Set to 5-20 for clinical translation |
| `min_panel` | 3 | Minimum genes to evaluate |
| `max_panel` | 50 | Maximum genes to evaluate |
| `n_folds` | 5 | Increase to 10 for small datasets |
| `validation` | `nested_cv` | Use `loocv` for n < 20 |
| `prevalence` | 0.05 | Set to actual disease prevalence |
| `random_state` | 42 | Change for sensitivity analysis |

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