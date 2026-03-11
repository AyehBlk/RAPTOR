# RAPTOR API Documentation

**Python API Reference for RAPTOR v2.2.0**

---

## 📚 **Table of Contents**

1. [Installation](#installation)
2. [Core Modules](#core-modules)
3. [Module 2: Quality Assessment API](#module-2-quality-assessment-api)
4. [Module 3: Data Profiler API](#module-3-data-profiler-api)
5. [Module 4: Recommender API](#module-4-recommender-api)
6. [Module 7: DE Import API](#module-7-de-import-api)
7. [Module 8: Optimization API](#module-8-optimization-api)
8. [Module 9: Ensemble API](#module-9-ensemble-api)
9. [Utilities API](#utilities-api)
10. [Error Handling](#error-handling)

---

## Installation

```python
# Install
pip install raptor-rnaseq

# Import
import raptor
print(raptor.__version__)  # '2.2.0'
```

---

## Core Modules

### **Package Structure**

```python
raptor/
├── quality_assessment.py    # Module 2
├── profiler.py             # Module 3
├── recommender.py          # Module 4 (rule-based)
├── ml_recommender.py       # Module 4 (ML-based)
├── de_import.py            # Module 7
├── parameter_optimization.py  # Module 8
├── ensemble.py             # Module 9
├── pipelines/              # Module 5
└── utils/                  # Utilities
```

---

## Module 2: Quality Assessment API

### **Classes**

#### `DataQualityAssessor`

**Constructor:**
```python
from raptor import DataQualityAssessor

assessor = DataQualityAssessor(
    normalization='tpm',        # 'tpm', 'fpkm', 'cpm', 'counts'
    consensus_threshold=3,      # Min methods to flag outlier
    plot_output=None,           # Path for plots
    verbose=True                # Print progress
)
```

**Methods:**
```python
# Assess quality
report = assessor.assess_quality(
    counts: pd.DataFrame,       # Count matrix (genes × samples)
    metadata: pd.DataFrame = None  # Sample metadata
) -> QualityReport

# Individual methods
outliers_mad = assessor.detect_outliers_mad(counts)
outliers_if = assessor.detect_outliers_isolation_forest(counts)
outliers_lof = assessor.detect_outliers_lof(counts)
outliers_pca = assessor.detect_outliers_pca(counts)
outliers_clustering = assessor.detect_outliers_clustering(counts)
outliers_statistical = assessor.detect_outliers_statistical(counts)

# Batch effect detection
batch_effects = assessor.assess_batch_effects(counts, metadata, batch_column='batch')

# Generate plots
assessor.generate_qc_plots(counts, metadata, output_path='qc_plots.pdf')
```

#### `QualityReport`

**Attributes:**
```python
report.outliers: List[str]              # List of outlier sample IDs
report.quality_issues: Dict             # Detected quality issues
report.recommendations: List[str]       # Recommended actions
report.metrics: pd.DataFrame            # QC metrics per sample
report.consensus_methods: Dict          # Methods agreeing per sample
```

**Methods:**
```python
# Export
report.to_json('qc_report.json')
report.to_html('qc_report.html')
report.summary()  # Print summary
```

### **Convenience Functions**

#### `quick_quality_check()`

```python
from raptor import quick_quality_check

report = quick_quality_check(
    counts: Union[str, pd.DataFrame],
    metadata: Union[str, pd.DataFrame] = None,
    output_dir: str = 'qc_results',
    normalization: str = 'tpm',
    consensus_threshold: int = 3,
    generate_plots: bool = True
) -> QualityReport
```

**Example:**
```python
report = quick_quality_check(
    counts='counts.csv',
    metadata='metadata.csv',
    output_dir='qc_output/',
    normalization='tpm',
    generate_plots=True
)

print(f"Outliers: {report.outliers}")
print(f"Issues: {report.quality_issues}")
```

---

## Module 3: Data Profiler API

### **Classes**

#### `RNAseqDataProfiler`

**Constructor:**
```python
from raptor import RNAseqDataProfiler

profiler = RNAseqDataProfiler(
    min_count=1,               # Minimum count threshold
    verbose=True               # Print progress
)
```

**Methods:**
```python
# Profile data
profile = profiler.profile(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str = 'condition'
) -> DataProfile

# Calculate individual features
bcv = profiler.calculate_bcv(counts, metadata, group_column)
sparsity = profiler.calculate_sparsity(counts)
library_sizes = profiler.calculate_library_sizes(counts)
```

#### `DataProfile`

**Attributes:**
```python
profile.features: Dict[str, float]      # All 32 features
profile.bcv: float                      # BCV (most important!)
profile.bcv_category: str               # 'low', 'moderate', 'high'
profile.n_samples: int                  # Number of samples
profile.n_genes: int                    # Number of genes
profile.min_group_size: int            # Minimum group size
profile.sparsity: float                # Sparsity (0-1)
profile.recommendations: List[str]      # Analysis recommendations
```

**Methods:**
```python
# Export
profile.to_json('profile.json')
profile.to_dict()
profile.summary()  # Print summary

# Get recommendation features
rec_features = profile.get_recommendation_features()
```

### **Convenience Functions**

#### `profile_data_quick()`

```python
from raptor import profile_data_quick

profile = profile_data_quick(
    counts: Union[str, pd.DataFrame],
    metadata: Union[str, pd.DataFrame],
    group_column: str = 'condition',
    output_dir: str = 'profile_results',
    min_count: int = 1
) -> DataProfile
```

**Example:**
```python
profile = profile_data_quick(
    counts='counts.csv',
    metadata='metadata.csv',
    group_column='condition',
    output_dir='profile/'
)

print(f"BCV: {profile.bcv:.3f}")
print(f"Category: {profile.bcv_category}")
print(f"Sample size: {profile.n_samples}")
print(f"Sparsity: {profile.sparsity:.1%}")

# Access all 32 features
for feature, value in profile.features.items():
    print(f"{feature}: {value}")
```

---

## Module 4: Recommender API

### **Rule-Based Recommender**

#### `PipelineRecommender`

**Constructor:**
```python
from raptor import PipelineRecommender

recommender = PipelineRecommender()
```

**Methods:**
```python
# Get recommendation
recommendation = recommender.recommend(
    profile: DataProfile,
    constraints: Dict = None  # e.g., {'memory_gb': 16, 'time_hours': 2}
) -> Recommendation
```

#### `Recommendation`

**Attributes:**
```python
recommendation.pipeline_name: str       # Recommended pipeline
recommendation.confidence: float        # Confidence score (0-1)
recommendation.reasoning: str           # Why this pipeline
recommendation.alternatives: List       # Alternative options
```

### **ML-Based Recommender**

#### `MLRecommender`

**Constructor:**
```python
from raptor import MLRecommender

ml_recommender = MLRecommender(
    model_path: str = None  # Path to trained model (optional)
)
```

**Methods:**
```python
# Get recommendation
recommendation = ml_recommender.recommend(
    profile: DataProfile
) -> MLRecommendation

# Train new model
ml_recommender.train(
    training_data: pd.DataFrame,
    labels: pd.Series
)

# Save/load model
ml_recommender.save_model('model.pkl')
ml_recommender.load_model('model.pkl')
```

#### `MLRecommendation`

**Attributes:**
```python
recommendation.pipeline_name: str
recommendation.confidence: float
recommendation.feature_importance: Dict    # Feature contributions
recommendation.probabilities: Dict         # All pipeline probabilities
```

### **Convenience Functions**

#### `recommend_pipeline()`

```python
from raptor import recommend_pipeline

recommendation = recommend_pipeline(
    profile_file: str,           # Path to profile JSON
    method: str = 'ml',          # 'ml', 'rule-based', or 'both'
    model_path: str = None
) -> Union[Recommendation, MLRecommendation]
```

**Example:**
```python
# ML-based (recommended)
rec = recommend_pipeline(
    profile_file='profile.json',
    method='ml'
)

print(f"Pipeline: {rec.pipeline_name}")
print(f"Confidence: {rec.confidence:.2f}")
print(f"Reasoning: {rec.feature_importance}")

# Rule-based
rec = recommend_pipeline(
    profile_file='profile.json',
    method='rule-based'
)
```

---

## Module 7: DE Import API

### **Import Functions**

#### Import Specific Formats

```python
from raptor import (
    import_deseq2,
    import_edger,
    import_limma,
    import_wilcoxon
)

# DESeq2
deseq2_result = import_deseq2(
    filepath: str,
    gene_column: str = 'gene_id',
    pvalue_column: str = 'pvalue',
    padj_column: str = 'padj',
    lfc_column: str = 'log2FoldChange',
    basemean_column: str = 'baseMean'
) -> DEResult

# edgeR
edger_result = import_edger(
    filepath: str,
    gene_column: str = 'gene_id',
    # Auto-detects other columns
) -> DEResult

# limma
limma_result = import_limma(
    filepath: str,
    gene_column: str = 'gene_id',
    # Auto-detects other columns
) -> DEResult

# Wilcoxon
wilcox_result = import_wilcoxon(
    filepath: str,
    gene_column: str = 'gene_id',
    pvalue_column: str = 'pvalue'
) -> DEResult
```

#### Generic Import

```python
from raptor import import_de_result

result = import_de_result(
    filepath: str,
    method: str = 'auto',  # 'auto', 'deseq2', 'edger', 'limma', 'custom'
    gene_column: str = 'gene_id',
    pvalue_column: str = None,  # Auto-detect if None
    padj_column: str = None,
    lfc_column: str = None,
    basemean_column: str = None
) -> DEResult
```

### **Classes**

#### `DEResult`

**Attributes:**
```python
result.data: pd.DataFrame           # Standardized DE results
result.method: str                  # Method name
result.n_genes: int                 # Number of genes
result.n_significant: int           # Significant at default threshold
```

**Methods:**
```python
# Filter
sig_genes = result.filter_significant(
    padj_threshold=0.05,
    lfc_threshold=0.0
)

# Export
result.to_csv('standardized_results.csv')
result.to_json('standardized_results.json')

# Summary
result.summary()
```

### **Comparison Functions**

#### `compare_de_results()`

```python
from raptor import compare_de_results

comparison = compare_de_results(
    **de_results: DEResult,  # Named DE results
    threshold: float = 0.05
) -> DEComparison

# Example
comparison = compare_de_results(
    deseq2=deseq2_result,
    edger=edger_result,
    limma=limma_result,
    threshold=0.05
)
```

#### `DEComparison`

**Attributes:**
```python
comparison.overlap_matrix: pd.DataFrame    # Gene overlap matrix
comparison.agreement_stats: Dict           # Agreement statistics
comparison.unique_genes: Dict              # Method-specific genes
comparison.consensus_genes: List           # Genes found by all
```

**Methods:**
```python
# Generate Venn diagram
comparison.plot_venn(output='venn.png')

# Get intersection
shared = comparison.get_intersection(min_methods=2)
```

#### `merge_de_results()`

```python
from raptor import merge_de_results

merged = merge_de_results(
    de_results: List[DEResult],
    how: str = 'outer'  # 'inner' or 'outer'
) -> pd.DataFrame
```

---

## Module 8: Optimization API

### **Optimization Functions**

#### Method 1: Ground Truth

```python
from raptor import optimize_with_ground_truth

result = optimize_with_ground_truth(
    de_result: DEResult,
    ground_truth: pd.DataFrame,  # DataFrame with true positives
    output_dir: str = 'optimization/ground_truth',
    metric: str = 'f1_score'  # 'f1_score', 'precision', 'recall'
) -> OptimizationResult
```

#### Method 2: FDR Control

```python
from raptor import optimize_with_fdr_control

result = optimize_with_fdr_control(
    de_result: DEResult,
    fdr_target: float = 0.05,
    output_dir: str = 'optimization/fdr'
) -> OptimizationResult
```

#### Method 3: Stability

```python
from raptor import optimize_with_stability

result = optimize_with_stability(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: str = 'optimization/stability',
    n_folds: int = 5
) -> OptimizationResult
```

#### Method 4: Reproducibility

```python
from raptor import optimize_with_reproducibility

result = optimize_with_reproducibility(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    cohort2: pd.DataFrame,  # Independent cohort counts
    output_dir: str = 'optimization/reproducibility'
) -> OptimizationResult
```

### **Classes**

#### `OptimizationResult`

**Attributes:**
```python
result.optimal_threshold: Dict          # {'padj': 0.05, 'lfc': 0.58}
result.metrics: Dict                    # Performance at optimal
result.search_history: pd.DataFrame     # Optimization trajectory
result.recommendations: List[str]       # Suggested thresholds
```

**Methods:**
```python
# Export
result.to_json('optimization_result.json')
result.plot_trajectory(output='optimization_plot.png')
result.summary()
```

---

## Module 9: Ensemble API

### **Ensemble Methods**

#### Fisher's Method

```python
from raptor import ensemble_fisher

result = ensemble_fisher(
    de_results_dict: Dict[str, DEResult],
    significance_threshold: float = 0.05,
    direction_threshold: float = 0.0,
    output_dir: str = 'ensemble/fisher'
) -> EnsembleResult
```

#### Brown's Method

```python
from raptor import ensemble_brown

result = ensemble_brown(
    de_results_dict: Dict[str, DEResult],
    significance_threshold: float = 0.05,
    output_dir: str = 'ensemble/brown'
) -> EnsembleResult
```

#### Robust Rank Aggregation

```python
from raptor import ensemble_rra

result = ensemble_rra(
    de_results_dict: Dict[str, DEResult],
    significance_threshold: float = 0.05,
    output_dir: str = 'ensemble/rra'
) -> EnsembleResult
```

#### Voting Consensus

```python
from raptor import ensemble_voting

result = ensemble_voting(
    de_results_dict: Dict[str, DEResult],
    min_methods: int = 2,  # Gene must be in ≥2 methods
    significance_threshold: float = 0.05,
    output_dir: str = 'ensemble/voting'
) -> EnsembleResult
```

#### Weighted Ensemble

```python
from raptor import ensemble_weighted

result = ensemble_weighted(
    de_results_dict: Dict[str, DEResult],
    weights: Dict[str, float],  # e.g., {'deseq2': 0.4, 'edger': 0.3, 'limma': 0.3}
    significance_threshold: float = 0.05,
    output_dir: str = 'ensemble/weighted'
) -> EnsembleResult
```

### **Classes**

#### `EnsembleResult`

**Attributes:**
```python
result.consensus_genes: List[str]       # Consensus DE genes
result.combined_pvalues: pd.DataFrame   # Combined p-values
result.meta_lfc: pd.DataFrame          # Meta-analysis log2FC
result.direction_consistency: pd.DataFrame  # Direction agreement
result.method_agreement: pd.DataFrame   # Per-gene agreement
```

**Methods:**
```python
# Export
result.to_csv('consensus_genes.csv')
result.to_json('ensemble_result.json')

# Plots
result.plot_venn(output='venn.png')
result.plot_upset(output='upset.png')

# Summary
result.summary()
```

### **Lower-Level Functions**

```python
from raptor import (
    fishers_method,
    browns_method,
    robust_rank_aggregation,
    check_direction_consistency,
    get_consensus_direction,
    calculate_meta_lfc
)

# Fisher's combination
combined_p = fishers_method([0.01, 0.03, 0.05])

# Brown's combination (with correlation)
combined_p = browns_method(
    pvalues=[0.01, 0.03, 0.05],
    correlation_matrix=corr_matrix
)

# Direction consistency
consistent, direction = check_direction_consistency({
    'deseq2': 1.5,  # log2FC
    'edger': 1.2,
    'limma': 1.4
})

# Meta log2FC
meta_lfc = calculate_meta_lfc(
    lfc_dict={'deseq2': 1.5, 'edger': 1.2},
    weights={'deseq2': 0.6, 'edger': 0.4}
)
```

---

## Utilities API

### **Validation**

```python
from raptor.utils.validation import (
    validate_count_matrix,
    validate_metadata,
    validate_group_column,
    validate_file_path,
    validate_positive_integer,
    validate_probability
)

# Validate count matrix
validate_count_matrix(
    counts: pd.DataFrame,
    min_genes: int = 100,
    min_samples: int = 2
)

# Validate metadata
validate_metadata(
    metadata: pd.DataFrame,
    counts: pd.DataFrame = None
)

# Validate group column
validate_group_column(
    metadata: pd.DataFrame,
    column: str,
    min_groups: int = 2,
    min_samples_per_group: int = 2
)
```

### **Error Handling**

```python
from raptor.utils.errors import (
    RAPTORError,
    ValidationError,
    PipelineError,
    OptimizationError,
    EnsembleError
)

try:
    result = ensemble_fisher(de_results)
except EnsembleError as e:
    print(f"Ensemble failed: {e}")
except ValidationError as e:
    print(f"Input validation failed: {e}")
```

---

## Error Handling

### **Exception Hierarchy**

```
RAPTORError (base)
├── ValidationError
├── DependencyError
├── PipelineError
├── OptimizationError
│   ├── GroundTruthError
│   └── InsufficientDataError
└── EnsembleError
    ├── MethodMismatchError
    ├── DirectionInconsistencyError
    └── CombinationFailedError
```

### **Error Handling Examples**

```python
from raptor import ensemble_fisher, ValidationError, EnsembleError

try:
    result = ensemble_fisher({
        'deseq2': deseq2_result,
        'edger': edger_result
    })
except ValidationError as e:
    print(f"Input validation failed: {e}")
    # Check your input data
except EnsembleError as e:
    print(f"Ensemble analysis failed: {e}")
    # Try different ensemble method
except Exception as e:
    print(f"Unexpected error: {e}")
    # Debug
```

---

## Complete Example Workflows

### **Workflow 1: QC → Profile → Recommend**

```python
from raptor import quick_quality_check, profile_data_quick, recommend_pipeline

# 1. QC
qc_report = quick_quality_check('counts.csv', 'metadata.csv')
if len(qc_report.outliers) > 0:
    print(f"Warning: {len(qc_report.outliers)} outliers detected")

# 2. Profile
profile = profile_data_quick('counts.csv', 'metadata.csv', group_column='condition')
print(f"BCV: {profile.bcv:.3f} ({profile.bcv_category})")

# 3. Recommend
rec = recommend_pipeline(profile_file='profile_results/data_profile.json', method='ml')
print(f"Recommended: {rec.pipeline_name} (confidence: {rec.confidence:.2f})")
```

### **Workflow 2: Import → Optimize → Ensemble**

```python
from raptor import import_deseq2, import_edger, import_limma
from raptor import optimize_with_fdr_control, ensemble_brown

# 1. Import DE results
deseq2 = import_deseq2('deseq2_results.csv')
edger = import_edger('edger_results.csv')
limma = import_limma('limma_results.csv')

# 2. Optimize thresholds
opt_result = optimize_with_fdr_control(deseq2, fdr_target=0.05)
print(f"Optimal FDR: {opt_result.optimal_threshold['padj']}")

# 3. Ensemble analysis
ensemble_result = ensemble_brown({
    'deseq2': deseq2,
    'edger': edger,
    'limma': limma
})
print(f"Consensus genes: {len(ensemble_result.consensus_genes)}")

# 4. Export
ensemble_result.to_csv('consensus_genes.csv')
```

---

## API Reference Summary

**Full API documentation for each module available in:**
- [Module 2 API](docs/MODULE_2_Quality Assessment & Outlier Detection.md#api-reference)
- [Module 3 API](docs/MODULE_3_Data_Profiling.md#api-reference)
- [Module 4 API](docs/MODULE_4_Pipeline_Recommender.md#api-reference)
- [Module 7 API](docs/MODULE_7_DE_Import.md#api-reference)
- [Module 8 API](docs/MODULE_8_Parameter_Optimization.md#api-reference)
- [Module 9 API](docs/MODULE_9_Ensemble_Analysis.md#api-reference)

**Version:** 2.2.0  
**Last Updated:** March 2026
