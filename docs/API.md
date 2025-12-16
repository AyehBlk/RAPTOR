# RAPTOR v2.1.1 Python API Reference

Complete API documentation for programmatic use with ML-powered recommendations and adaptive threshold optimization.

## Installation
```python
pip install raptor-rnaseq>=2.1.1
import raptor
```

## What's New in v2.1.1

-  **Adaptive Threshold Optimizer** - Data-driven threshold selection for DE analysis
-  **Publication Methods Text** - Auto-generated methods paragraphs
-  **Threshold Visualization** - Volcano plots with optimized cutoffs

---

## 🎯 Threshold Optimizer Module (NEW)

### AdaptiveThresholdOptimizer

Data-driven threshold selection for differential expression analysis.

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
import pandas as pd

# Load DE results from any tool (DESeq2, edgeR, limma, etc.)
df = pd.read_csv('deseq2_results.csv')

# Create optimizer
ato = AdaptiveThresholdOptimizer(
    data=df,
    logfc_col='log2FoldChange',  # Column name for logFC
    pvalue_col='pvalue'           # Column name for p-values
)

# Optimize thresholds
result = ato.optimize(
    goal='discovery',           # 'discovery', 'balanced', or 'validation'
    padj_method='BH',           # P-value adjustment method
    logfc_method='auto',        # LogFC optimization method
    fdr_threshold=0.05          # Target FDR level
)

# Access results
print(f"Optimal logFC threshold: {result.logfc_threshold:.3f}")
print(f"Optimal p-value threshold: {result.pvalue_threshold}")
print(f"Significant genes: {result.n_significant}")
print(f"π₀ estimate: {result.pi0:.3f}")
print(f"Up-regulated: {result.n_up}")
print(f"Down-regulated: {result.n_down}")

# Get publication-ready methods text
print(result.methods_text)

# Get results dataframe with significance flags
results_df = result.results_df
results_df.to_csv('optimized_results.csv')
```

**Constructor Parameters:**
- `data` (pd.DataFrame): DE results with logFC and p-value columns
- `logfc_col` (str): Column name for log fold change
- `pvalue_col` (str): Column name for p-values

**optimize() Parameters:**
- `goal` (str): Analysis goal - 'discovery', 'balanced', or 'validation'
- `padj_method` (str): P-value adjustment - 'BH', 'BY', 'storey', 'holm', 'hochberg', 'bonferroni'
- `logfc_method` (str): LogFC method - 'auto', 'mad', 'mixture', 'power', 'percentile'
- `fdr_threshold` (float): Target FDR level (default: 0.05)

**ThresholdResult Fields:**
```python
@dataclass
class ThresholdResult:
    logfc_threshold: float      # Optimal |logFC| cutoff
    pvalue_threshold: float     # P-value cutoff used
    padj_method: str            # Adjustment method used
    logfc_method: str           # LogFC method used
    pi0: float                  # Estimated proportion of true nulls
    n_significant: int          # Number of significant genes
    n_up: int                   # Up-regulated count
    n_down: int                 # Down-regulated count
    goal: str                   # Analysis goal used
    methods_text: str           # Publication-ready paragraph
    results_df: pd.DataFrame    # Full results with flags
```

---

### optimize_thresholds() Convenience Function

Quick threshold optimization without class instantiation.

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load data
df = pd.read_csv('edger_results.csv')

# Quick optimization
result = optimize_thresholds(
    data=df,
    logfc_col='logFC',      # edgeR column name
    pvalue_col='PValue',    # edgeR column name
    goal='balanced'
)

print(f"Use |logFC| > {result.logfc_threshold:.2f}, padj < {result.pvalue_threshold}")
```

---

### Column Name Reference

Different tools use different column names:

```python
# DESeq2
result = optimize_thresholds(df, logfc_col='log2FoldChange', pvalue_col='pvalue')

# edgeR
result = optimize_thresholds(df, logfc_col='logFC', pvalue_col='PValue')

# limma
result = optimize_thresholds(df, logfc_col='logFC', pvalue_col='P.Value')

# NOISeq (uses probability, not p-value)
# Note: NOISeq requires special handling
result = optimize_thresholds(df, logfc_col='M', pvalue_col='prob')
```

---

### Analysis Goals Explained

```python
# DISCOVERY - Maximize sensitivity (exploratory)
# - More permissive thresholds
# - Catches more true positives
# - Higher false positive rate
result = ato.optimize(goal='discovery')

# BALANCED - Standard analysis (publication)
# - Balanced sensitivity/specificity
# - Standard FDR control
# - Recommended for most analyses
result = ato.optimize(goal='balanced')

# VALIDATION - Maximize specificity (confirmation)
# - Stringent thresholds
# - Fewer false positives
# - May miss some true positives
result = ato.optimize(goal='validation')
```

---

### P-value Adjustment Methods

```python
# Benjamini-Hochberg (default, recommended)
result = ato.optimize(padj_method='BH')

# Benjamini-Yekutieli (for correlated tests)
result = ato.optimize(padj_method='BY')

# Storey q-value (more power with π₀ estimation)
result = ato.optimize(padj_method='storey')

# Holm (FWER control, step-down)
result = ato.optimize(padj_method='holm')

# Hochberg (FWER control, step-up)
result = ato.optimize(padj_method='hochberg')

# Bonferroni (most conservative)
result = ato.optimize(padj_method='bonferroni')
```

---

### LogFC Methods

```python
# Auto/Consensus (default, recommended)
# Uses consensus of all methods
result = ato.optimize(logfc_method='auto')

# MAD-based (robust to outliers)
# threshold = median + k * MAD
result = ato.optimize(logfc_method='mad')

# Mixture model (Gaussian mixture)
# Separates DE and non-DE distributions
result = ato.optimize(logfc_method='mixture')

# Power-based (statistical power)
# Minimum effect size for desired power
result = ato.optimize(logfc_method='power')

# Percentile (distribution-based)
# 95th percentile of null distribution
result = ato.optimize(logfc_method='percentile')
```

---

### Visualization

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
from raptor.threshold_optimizer.visualization import (
    plot_volcano,
    plot_pvalue_distribution,
    plot_logfc_distribution,
    plot_threshold_comparison
)

# Create optimizer and get result
ato = AdaptiveThresholdOptimizer(df, 'log2FoldChange', 'pvalue')
result = ato.optimize(goal='discovery')

# Volcano plot with optimized thresholds
fig = plot_volcano(
    result.results_df,
    logfc_threshold=result.logfc_threshold,
    pvalue_threshold=result.pvalue_threshold,
    title="Volcano Plot with Optimized Thresholds"
)
fig.savefig('volcano.png', dpi=300)

# P-value distribution with π₀
fig = plot_pvalue_distribution(df['pvalue'], result.pi0)
fig.savefig('pvalue_dist.png', dpi=300)

# LogFC distribution with threshold
fig = plot_logfc_distribution(df['log2FoldChange'], result.logfc_threshold)
fig.savefig('logfc_dist.png', dpi=300)

# Threshold comparison heatmap
fig = plot_threshold_comparison(df, 'log2FoldChange', 'pvalue')
fig.savefig('threshold_comparison.png', dpi=300)
```

---

## Core Classes

### RNAseqDataProfiler

Profile RNA-seq count data for ML-powered pipeline recommendation.

```python
from raptor import RNAseqDataProfiler
import pandas as pd

# Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Create profiler with ML features
profiler = RNAseqDataProfiler(
    counts, 
    metadata,
    use_ml=True,
    monitor_resources=True
)

# Run profiling with quality assessment
profile = profiler.profile(
    quality_check=True,
    detect_outliers=True
)

# Access results
print(f"BCV: {profile['bcv']:.3f}")
print(f"Quality Score: {profile['quality_score']}/100")
print(f"ML Confidence: {profile['ml_confidence']:.2f}")
```

---

### MLPipelineRecommender

ML-powered pipeline recommendations with confidence scores.

```python
from raptor import MLPipelineRecommender

# Initialize
recommender = MLPipelineRecommender(
    model_path='default',
    use_ensemble=True
)

# Get recommendations
recommendations = recommender.recommend(
    profile,
    n=3,
    explain=True
)

# Access recommendation
top = recommendations[0]
print(f"Pipeline: {top['pipeline_name']}")
print(f"Confidence: {top['confidence']}")
print(f"Success Rate: {top['historical_success']:.1%}")

# After pipeline runs, use ATO for thresholds (NEW)
from raptor.threshold_optimizer import optimize_thresholds
de_results = pd.read_csv('pipeline_results.csv')
thresholds = optimize_thresholds(de_results, goal='balanced')
```

---

### QualityAssessor

Comprehensive data quality assessment.

```python
from raptor import QualityAssessor

assessor = QualityAssessor()

qc_results = assessor.assess(
    counts,
    metadata,
    check_contamination=True,
    check_batch_effects=True,
    check_outliers=True
)

print(f"Overall Score: {qc_results['overall_score']}/100")
```

---

### ResourceMonitor

Real-time resource usage tracking.

```python
from raptor import ResourceMonitor

monitor = ResourceMonitor(interval=5, log_file='resources.log')
monitor.start()

# Run analysis
# ...

stats = monitor.stop()
print(f"Peak Memory: {stats['peak_memory_gb']:.1f} GB")
```

---

### EnsembleAnalyzer

Combine results from multiple pipelines.

```python
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

# Combine pipeline results
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(
    results,
    method='vote',
    weights={'salmon_edgeR': 0.4, 'star_deseq2': 0.4, 'kallisto_sleuth': 0.2}
)

# Use ATO for uniform thresholds across ensemble (NEW)
thresholds = optimize_thresholds(ensemble['combined_results'], goal='balanced')
```

---

## Complete Workflow Examples

### Example 1: Full Analysis with ATO (NEW)

```python
from raptor import RNAseqDataProfiler, MLPipelineRecommender
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# 1. Load and profile data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile(quality_check=True)
print(f"Quality Score: {profile['quality_score']}/100")

# 2. Get ML recommendation
recommender = MLPipelineRecommender()
recommendations = recommender.recommend(profile, n=3)
print(f"Recommended: {recommendations[0]['pipeline_name']}")

# 3. [Run recommended pipeline - produces DE results]
# raptor run --pipeline 3 ...

# 4. Optimize thresholds (NEW in v2.1.1)
de_results = pd.read_csv('deseq2_results.csv')
result = optimize_thresholds(
    de_results,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)

print(f"\n🎯 Optimized Thresholds:")
print(f"   LogFC: |{result.logfc_threshold:.3f}|")
print(f"   P-value: {result.pvalue_threshold}")
print(f"   Significant: {result.n_significant} genes")

# 5. Save results with methods text
result.results_df.to_csv('final_results.csv')
with open('methods.txt', 'w') as f:
    f.write(result.methods_text)
```

### Example 2: Publication-Ready Analysis

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
from raptor.threshold_optimizer.visualization import plot_volcano
import pandas as pd

# Load DE results
df = pd.read_csv('deseq2_results.csv')

# Optimize for publication
ato = AdaptiveThresholdOptimizer(df, 'log2FoldChange', 'pvalue')
result = ato.optimize(goal='balanced')

# Generate all outputs
# 1. Methods text for paper
print("=== METHODS SECTION ===")
print(result.methods_text)

# 2. Results summary
print(f"\n=== RESULTS ===")
print(f"Significant genes: {result.n_significant}")
print(f"Up-regulated: {result.n_up}")
print(f"Down-regulated: {result.n_down}")
print(f"π₀ estimate: {result.pi0:.3f}")

# 3. Volcano plot
fig = plot_volcano(
    result.results_df,
    result.logfc_threshold,
    result.pvalue_threshold
)
fig.savefig('Figure_1_volcano.png', dpi=300, bbox_inches='tight')

# 4. Supplementary table
result.results_df.to_csv('Supplementary_Table_DE_genes.csv')
```

### Example 3: Comparing Different Goals

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
import pandas as pd

df = pd.read_csv('de_results.csv')
ato = AdaptiveThresholdOptimizer(df, 'log2FoldChange', 'pvalue')

# Compare all goals
for goal in ['discovery', 'balanced', 'validation']:
    result = ato.optimize(goal=goal)
    print(f"\n{goal.upper()}:")
    print(f"  LogFC threshold: {result.logfc_threshold:.3f}")
    print(f"  Significant genes: {result.n_significant}")
```

---

## Configuration

### v2.1.1 Configuration Options

```python
from raptor.utils import load_config, save_config

config = load_config('config.yaml')

# NEW: Threshold Optimizer settings
config['threshold_optimizer'] = {
    'enabled': True,
    'goal': 'balanced',
    'default_padj_method': 'BH',
    'default_logfc_method': 'auto',
    'fdr_threshold': 0.05
}

# NEW: Use adaptive thresholds in statistics
config['statistics'] = {
    'use_adaptive_thresholds': True,  # NEW
    'fdr_threshold': 0.05,
    'log2fc_threshold': 1.0  # Fallback if ATO disabled
}

save_config(config, 'config_v2.1.1.yaml')
```

---

## Error Handling

```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer, optimize_thresholds

try:
    result = optimize_thresholds(df, 'log2FoldChange', 'pvalue', goal='balanced')
    
except KeyError as e:
    print(f"Column not found: {e}")
    print(f"Available columns: {df.columns.tolist()}")
    
except ValueError as e:
    print(f"Invalid parameter: {e}")
    
except Exception as e:
    print(f"Optimization failed: {e}")
```

---

## Version Info

```python
import raptor

print(raptor.__version__)  # "2.1.1"

# Check ATO availability
try:
    from raptor.threshold_optimizer import optimize_thresholds
    print("✅ Threshold Optimizer available")
except ImportError:
    print("❌ Threshold Optimizer not available - update RAPTOR")
```

---

## Migration from v2.1.0

### Key Changes:

```python
# v2.1.0 - Manual thresholds
fdr_threshold = 0.05
log2fc_threshold = 1.0  # Arbitrary!

# v2.1.1 - Data-driven thresholds (NEW!)
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(de_results, goal='balanced')
fdr_threshold = result.pvalue_threshold
log2fc_threshold = result.logfc_threshold  # Optimized!
```

### Recommended Upgrades:

1. **Use ATO for thresholds**: Replace arbitrary cutoffs
2. **Include methods text**: Use `result.methods_text` in papers
3. **Enable in config**: Set `use_adaptive_thresholds: true`
4. **Try dashboard**: Use 🎯 Threshold Optimizer page

---

## See Also

- [THRESHOLD_OPTIMIZER.md](THRESHOLD_OPTIMIZER.md) - Comprehensive ATO guide
- [DASHBOARD.md](DASHBOARD.md) - Interactive dashboard with ATO page
- [ML_TRAINING.md](ML_TRAINING.md) - Train custom ML models
- [ENSEMBLE.md](ENSEMBLE.md) - Ensemble analysis with uniform thresholds

---

**RAPTOR v2.1.1**  
**Author**: Ayeh Bolouki  
**License**: MIT

**Upgrade today for data-driven thresholds!** 🎯🦖
