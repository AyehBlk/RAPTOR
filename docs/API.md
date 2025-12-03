# RAPTOR v2.1.0 Python API Reference

Complete API documentation for programmatic use with ML-powered recommendations.

## Installation
```python
pip install raptor-rnaseq>=2.1.0
import raptor
```

## What's New in v2.1.0

-  **ML-Powered Recommendations** - Machine learning models for smarter pipeline selection
-  **Interactive Dashboard** - Web-based visualization and monitoring  
-  **Resource Monitoring** - Real-time CPU/memory tracking
-  **Parameter Optimization** - Automated parameter tuning
-  **Quality Assessment** - Comprehensive data QC
-  **Ensemble Analysis** - Combine results from multiple pipelines

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
    use_ml=True,  # NEW in v2.1.0: Enable ML recommendations
    monitor_resources=True  # NEW: Track resource usage
)

# Run profiling with quality assessment
profile = profiler.profile(
    quality_check=True,  # NEW: Run QC
    detect_outliers=True  # NEW: Identify problematic samples
)

# Access results
print(f"BCV: {profile['bcv']:.3f}")
print(f"Mean depth: {profile['mean_depth']:.0f}")
print(f"Samples: {profile['n_samples']}")
print(f"Quality Score: {profile['quality_score']}/100")  # NEW
print(f"ML Confidence: {profile['ml_confidence']:.2f}")  # NEW
```

**New Methods in v2.1.0:**
- `profile(quality_check=False, detect_outliers=False)`: Enhanced profiling
- `assess_quality()`: Run comprehensive quality assessment
- `detect_contamination()`: Check for contamination
- `check_batch_effects()`: Analyze batch effects
- `monitor_resources()`: Track resource usage
- `export_for_ml()`: Prepare data for ML model training

**Profile Output (Enhanced):**
```python
{
    # Original metrics
    'n_samples': 12,
    'n_genes': 15234,
    'bcv': 0.42,
    'bcv_category': 'medium',
    'mean_depth': 25000000,
    'depth_category': 'high',
    'zero_inflation': 0.15,
    'library_size_cv': 0.18,
    
    # NEW in v2.1.0
    'quality_score': 87,  # Overall quality (0-100)
    'quality_flags': ['minor_batch_effect'],  # Warnings
    'ml_confidence': 0.89,  # ML recommendation confidence
    'contamination_detected': False,
    'outliers': ['Sample_7'],  # Flagged samples
    'batch_effect_severity': 'minor',  # None/minor/moderate/severe
    'resource_usage': {
        'cpu_percent': 45,
        'memory_gb': 2.3,
        'runtime_seconds': 8.5
    }
}
```

---

### MLPipelineRecommender (NEW)

ML-powered pipeline recommendations with confidence scores.

```python
from raptor import MLPipelineRecommender

# Initialize with pre-trained model
recommender = MLPipelineRecommender(
    model_path='default',  # or path to custom model
    use_ensemble=True  # NEW: Use ensemble of models
)

# Get ML-powered recommendations
recommendations = recommender.recommend(
    profile,
    n=3,  # Top 3
    explain=True  # NEW: Get feature importance
)

# Access recommendation
top = recommendations[0]
print(f"Pipeline: {top['pipeline_name']}")
print(f"ML Score: {top['ml_score']:.3f}")
print(f"Confidence: {top['confidence']}")  # high/medium/low
print(f"Similar Projects: {top['similar_count']}")  # NEW
print(f"Success Rate: {top['historical_success']:.1%}")  # NEW
```

**Methods:**
- `recommend(profile, n=3, explain=False)`: Get ML recommendations
- `predict_proba(profile)`: Get probability for each pipeline
- `explain_recommendation(profile, pipeline_id)`: Feature importance
- `get_similar_projects(profile, n=10)`: Find similar past analyses
- `train_custom_model(X, y, save_path)`: Train on your data

**Recommendation Output (Enhanced):**
```python
[
    {
        'pipeline_id': 3,
        'pipeline_name': 'Salmon-edgeR',
        
        # ML scores (NEW)
        'ml_score': 0.89,  # ML model prediction
        'confidence': 'high',  # Based on training data similarity
        'similar_count': 1247,  # Similar projects in training
        'historical_success': 0.87,  # Success rate on similar data
        
        # Original scores
        'score': 0.88,  # Rule-based score
        'reasoning': 'Excellent balance...',
        'expected_runtime': 1.2,
        'expected_memory': 10,
        'expected_accuracy': 0.90,
        
        # Feature importance (NEW)
        'top_features': [
            ('bcv', 0.35),
            ('n_samples', 0.22),
            ('depth', 0.18)
        ]
    },
    ...
]
```

---

### QualityAssessor (NEW)

Comprehensive data quality assessment.

```python
from raptor import QualityAssessor

assessor = QualityAssessor()

# Run quality assessment
qc_results = assessor.assess(
    counts,
    metadata,
    check_contamination=True,
    check_batch_effects=True,
    check_outliers=True
)

# Access results
print(f"Overall Score: {qc_results['overall_score']}/100")
print(f"Issues Found: {len(qc_results['issues'])}")

# Check specific metrics
if qc_results['library_size_cv'] > 0.3:
    print("WARNING: High library size variation")

if qc_results['contamination_detected']:
    print(f"ALERT: Contamination detected - {qc_results['contamination_pct']:.1%}")
```

**Methods:**
- `assess(counts, metadata, **checks)`: Run full QC
- `check_library_sizes()`: Analyze library size distribution
- `check_contamination()`: Detect contamination
- `check_batch_effects()`: Identify batch effects
- `detect_outliers()`: Find outlier samples
- `generate_qc_report()`: Create HTML report

---

### ResourceMonitor (NEW)

Real-time resource usage tracking.

```python
from raptor import ResourceMonitor

# Start monitoring
monitor = ResourceMonitor(
    interval=5,  # Check every 5 seconds
    log_file='resources.log'
)

monitor.start()

# Run analysis
# ... your code here ...

# Stop monitoring
stats = monitor.stop()

print(f"Peak CPU: {stats['peak_cpu_percent']}%")
print(f"Peak Memory: {stats['peak_memory_gb']:.1f} GB")
print(f"Avg CPU: {stats['avg_cpu_percent']}%")
print(f"Total Runtime: {stats['total_runtime']:.1f}s")

# Generate resource report
monitor.generate_report('resource_report.html')
```

**Methods:**
- `start()`: Begin monitoring
- `stop()`: Stop and return statistics
- `get_current_usage()`: Get instant snapshot
- `generate_report(output)`: Create visualization
- `export_timeseries()`: Export data for plotting
- `predict_requirements(profile)`: Estimate future needs

---

### ParameterOptimizer (NEW)

Automated parameter tuning for optimal results.

```python
from raptor import ParameterOptimizer

optimizer = ParameterOptimizer(
    method='bayesian',  # or 'grid', 'random'
    metric='f1_score'  # What to optimize
)

# Optimize parameters
best_params = optimizer.optimize(
    counts,
    metadata,
    ground_truth=validation_data,  # Optional
    pipeline_id=3,
    n_iterations=30
)

print(f"Best FDR: {best_params['fdr_threshold']}")
print(f"Best Log2FC: {best_params['log2fc_threshold']}")
print(f"Expected F1: {best_params['expected_f1']:.3f}")

# Run sensitivity analysis
sensitivity = optimizer.analyze_sensitivity(
    counts,
    metadata,
    parameters=['fdr_threshold', 'log2fc_threshold']
)
```

**Methods:**
- `optimize(counts, metadata, **kwargs)`: Find best parameters
- `grid_search(param_grid)`: Exhaustive search
- `random_search(param_distributions, n_iter)`: Random sampling
- `bayesian_optimize(bounds, n_iterations)`: Smart search
- `analyze_sensitivity(parameters)`: Parameter impact
- `cross_validate(params, folds=5)`: Validate parameters

---

### EnsembleAnalyzer (NEW)

Combine results from multiple pipelines.

```python
from raptor import EnsembleAnalyzer

# Run multiple pipelines
results = {
    'salmon_edgeR': pipeline3_results,
    'star_deseq2': pipeline1_results,
    'kallisto_sleuth': pipeline4_results
}

# Combine results
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(
    results,
    method='vote',  # or 'rank', 'score', 'weighted'
    weights={'salmon_edgeR': 0.4, 'star_deseq2': 0.4, 'kallisto_sleuth': 0.2}
)

# Get consensus genes
consensus_genes = ensemble.get_consensus(
    min_agreement=2  # In at least 2/3 pipelines
)

print(f"Consensus DE genes: {len(consensus_genes)}")
print(f"High confidence: {ensemble['high_confidence_count']}")
```

**Methods:**
- `combine(results, method, weights)`: Merge multiple results
- `get_consensus(min_agreement)`: Find agreed genes
- `calculate_confidence()`: Per-gene confidence scores
- `identify_discordant()`: Genes with disagreement
- `generate_venn()`: Create Venn diagrams
- `generate_report()`: Comprehensive comparison

---

### DashboardServer (NEW)

Interactive web-based dashboard.

```python
from raptor import DashboardServer

# Start dashboard
dashboard = DashboardServer(
    port=8501,
    results_dir='results/',
    auto_refresh=True
)

dashboard.start()
# Navigate to http://localhost:8501

# Dashboard features:
# - Real-time resource monitoring
# - Interactive quality plots
# - Parameter exploration
# - Result visualization
# - Live progress tracking
```

**Methods:**
- `start()`: Launch dashboard
- `stop()`: Shut down server
- `add_plot(name, fig)`: Add custom visualization
- `update_data(results)`: Refresh displayed data
- `export_session()`: Save dashboard state

---

### PipelineBenchmark

Enhanced benchmarking with resource tracking (Updated for v2.1.0).

```python
from raptor import PipelineBenchmark

benchmark = PipelineBenchmark(
    data_dir='fastq/',
    output_dir='results/',
    threads=8,
    memory='32G',
    monitor_resources=True,  # NEW: Track resources
    use_ml_ranking=True  # NEW: ML-based ranking
)

# Run pipelines with monitoring
results = benchmark.run_pipelines(
    [1, 3, 4],
    quality_check=True,  # NEW: QC before/after
    track_intermediate=True  # NEW: Save checkpoints
)

# Get ML-ranked results (NEW)
ranked = benchmark.rank_results(
    results,
    use_ml=True,
    criteria=['accuracy', 'speed', 'resources']
)

# Save comprehensive results
benchmark.save_results(
    results,
    include_resources=True,  # NEW
    include_qc=True  # NEW
)
```

**New Methods in v2.1.0:**
- `run_with_monitoring()`: Track resources during run
- `rank_results(use_ml=False)`: ML-powered ranking
- `generate_comparison_dashboard()`: Interactive comparison
- `estimate_resources(profile)`: Predict requirements
- `save_for_ml_training()`: Export for model training

---

## Utility Functions

### Enhanced for v2.1.0

```python
from raptor.utils import *

# Original utilities
ensure_dir('output/')
check_file_exists('data.csv')
check_command_exists('STAR')
check_required_tools(['STAR', 'salmon', 'kallisto'])

# NEW in v2.1.0: Resource utilities
available_memory = get_available_memory()  # Returns GB
available_cpus = get_cpu_count()  
available_disk = get_available_disk_space()

# NEW: ML utilities
features = extract_ml_features(counts, metadata)
scaled_features = scale_features(features)
predictions = load_ml_model('default').predict(scaled_features)

# NEW: Quality check utilities
qc_passed = quick_quality_check(counts)
outliers = detect_outliers_zscore(counts)
contamination_pct = estimate_contamination(counts, organism='human')

# NEW: Parameter utilities
optimal_threads = suggest_threads(n_samples=12)
optimal_memory = suggest_memory(n_genes=20000, n_samples=12)

# Data validation
validate_count_matrix(counts)
validate_metadata(metadata, counts)

# Configuration
config = load_config('config.yaml')
save_config(config, 'my_config.yaml')

# Formatting
format_time(3665)  # "1h 1m 5s"
format_bytes(1536000)  # "1.5 MB"
format_percent(0.8542)  # "85.4%"
```

---

## Complete Workflow Examples

### Example 1: ML-Powered Analysis (NEW)

```python
import pandas as pd
from raptor import (
    RNAseqDataProfiler,
    MLPipelineRecommender,
    QualityAssessor,
    ResourceMonitor,
    ReportGenerator
)

# 1. Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# 2. Quality assessment (NEW)
assessor = QualityAssessor()
qc_results = assessor.assess(counts, metadata)

if qc_results['overall_score'] < 70:
    print("WARNING: Low quality data!")
    print(qc_results['issues'])

# 3. Profile with ML (NEW)
monitor = ResourceMonitor()
monitor.start()

profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile(quality_check=True)

resource_stats = monitor.stop()

# 4. Get ML recommendations (NEW)
ml_recommender = MLPipelineRecommender()
recommendations = ml_recommender.recommend(
    profile, 
    n=3, 
    explain=True
)

# 5. Show confidence and reasoning
top = recommendations[0]
print(f"🥇 {top['pipeline_name']}")
print(f"   ML Score: {top['ml_score']:.3f}")
print(f"   Confidence: {top['confidence']}")
print(f"   Based on {top['similar_count']} similar projects")
print(f"   Historical success: {top['historical_success']:.1%}")

# 6. Generate comprehensive report (NEW features)
report_gen = ReportGenerator()
report_gen.generate_ml_report(
    profile,
    recommendations,
    qc_results,
    resource_stats,
    output='ml_analysis_report.html'
)
```

### Example 2: Ensemble Analysis (NEW)

```python
from raptor import (
    PipelineBenchmark,
    EnsembleAnalyzer,
    ReportGenerator
)

# 1. Run top 3 pipelines
benchmark = PipelineBenchmark(
    data_dir='fastq/',
    output_dir='ensemble_results/',
    monitor_resources=True
)

results = benchmark.run_pipelines([1, 3, 4])

# 2. Combine with ensemble (NEW)
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(
    results,
    method='weighted',
    weights={1: 0.4, 3: 0.4, 4: 0.2}  # Based on ML recommendations
)

# 3. Get high-confidence genes
consensus = ensemble.get_consensus(min_agreement=2)
high_conf = ensemble.get_high_confidence(threshold=0.8)

print(f"Consensus genes: {len(consensus)}")
print(f"High confidence: {len(high_conf)}")

# 4. Generate ensemble report (NEW)
report_gen = ReportGenerator()
report_gen.generate_ensemble_report(
    ensemble,
    output='ensemble_report.html'
)
```

### Example 3: Parameter Optimization (NEW)

```python
from raptor import ParameterOptimizer, PipelineBenchmark

# 1. Optimize parameters
optimizer = ParameterOptimizer(method='bayesian')

best_params = optimizer.optimize(
    counts,
    metadata,
    ground_truth=validated_genes,
    pipeline_id=3,
    n_iterations=30
)

print(f"Optimized FDR: {best_params['fdr_threshold']}")
print(f"Optimized Log2FC: {best_params['log2fc_threshold']}")
print(f"Expected F1-Score: {best_params['expected_f1']:.3f}")

# 2. Run with optimized parameters
benchmark = PipelineBenchmark(data_dir='fastq/')
results = benchmark.run_single_pipeline(
    pipeline_id=3,
    custom_params=best_params
)

# 3. Validate improvement
print(f"Baseline F1: 0.85")
print(f"Optimized F1: {results['f1_score']:.3f}")
print(f"Improvement: {((results['f1_score'] - 0.85) / 0.85 * 100):.1f}%")
```

### Example 4: Interactive Dashboard (NEW)

```python
from raptor import DashboardServer, ResourceMonitor

# Start real-time monitoring dashboard
dashboard = DashboardServer(
    port=8501,
    results_dir='results/',
    auto_refresh=True
)

dashboard.start()
# Open browser to http://localhost:8501

# Dashboard shows:
# - Real-time resource usage graphs
# - Live progress updates
# - Interactive quality plots
# - Parameter exploration tools
# - Result comparisons
```

---

## Configuration

### v2.1.0 Configuration Options (NEW)

```python
from raptor.utils import load_config, save_config

# Load config
config = load_config('config/config.yaml')

# NEW ML settings
config['ml'] = {
    'use_ml_recommendations': True,
    'model_path': 'default',  # or custom path
    'min_confidence': 0.7,
    'ensemble_models': True,
    'explain_predictions': True
}

# NEW monitoring settings
config['monitoring'] = {
    'track_resources': True,
    'interval_seconds': 5,
    'log_to_file': True,
    'generate_plots': True
}

# NEW quality settings
config['quality'] = {
    'run_qc_checks': True,
    'min_quality_score': 70,
    'detect_contamination': True,
    'check_batch_effects': True,
    'outlier_threshold': 3.0
}

# NEW optimization settings
config['optimization'] = {
    'enable_param_optimization': False,
    'optimization_method': 'bayesian',
    'n_iterations': 30,
    'cross_validate': True
}

# Modify existing settings
config['resources']['default_threads'] = 16
config['statistics']['fdr_threshold'] = 0.01

# Save
save_config(config, 'my_v2.1_config.yaml')

# Use in analysis
profiler = RNAseqDataProfiler(
    counts, 
    metadata,
    config_file='my_v2.1_config.yaml'
)
```

---

## Error Handling

```python
from raptor import RNAseqDataProfiler, RaptorError, QualityError
from raptor.utils import validate_count_matrix

try:
    # Validate first
    validate_count_matrix(counts)
    
    # Quality check (NEW)
    qc_results = QualityAssessor().assess(counts, metadata)
    if qc_results['overall_score'] < 60:
        raise QualityError(f"Low quality: {qc_results['issues']}")
    
    # Profile with ML
    profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
    profile = profiler.profile()
    
except FileNotFoundError as e:
    print(f"File not found: {e}")
except ValueError as e:
    print(f"Invalid data: {e}")
except QualityError as e:  # NEW
    print(f"Quality check failed: {e}")
except RaptorError as e:  # NEW: Base exception
    print(f"RAPTOR error: {e}")
except Exception as e:
    print(f"Unexpected error: {e}")
```

---

## Logging

```python
from raptor.utils import setup_logging
import logging

# Setup logging with more levels (NEW)
setup_logging(
    level='DEBUG',
    log_file='raptor.log',
    log_resource_usage=True,  # NEW
    log_ml_predictions=True  # NEW
)

# Use logger
logger = logging.getLogger('raptor')
logger.info("Starting ML-powered analysis...")
logger.debug(f"ML model loaded: {model_path}")
logger.warning(f"Low confidence: {confidence:.2f}")
```

---

## Version Info

```python
import raptor

print(raptor.__version__)  # "2.1.0"
print(raptor.__author__)
print(raptor.__email__)

# Get environment info (enhanced)
from raptor.utils import get_environment_info
info = get_environment_info(include_ml=True)  # NEW flag
print(info)

# Check ML model info (NEW)
from raptor.ml import get_model_info
model_info = get_model_info()
print(f"Model version: {model_info['version']}")
print(f"Training size: {model_info['training_samples']}")
print(f"Accuracy: {model_info['validation_accuracy']:.3f}")
```

---

## Type Hints

All classes use type hints (enhanced in v2.1.0):

```python
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import pandas as pd

def my_ml_analysis(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: Path,
    use_ml: bool = True,  # NEW
    optimize_params: bool = False  # NEW
) -> Dict[str, float]:
    """Run ML-powered analysis.
    
    Args:
        counts: Gene expression matrix
        metadata: Sample information
        output_dir: Output directory
        use_ml: Whether to use ML recommendations (NEW)
        optimize_params: Whether to optimize parameters (NEW)
        
    Returns:
        Dictionary with performance metrics
    """
    profiler = RNAseqDataProfiler(counts, metadata, use_ml=use_ml)
    return profiler.profile()
```

---

## Migration from v2.0.0

### Key Changes:

```python
# v2.0.0 (OLD)
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()
recommender = PipelineRecommender()
recommendations = recommender.recommend(profile)

# v2.1.0 (NEW - backward compatible!)
profiler = RNAseqDataProfiler(
    counts, 
    metadata,
    use_ml=True  # NEW: Enable ML
)
profile = profiler.profile(
    quality_check=True  # NEW: Run QC
)
recommender = MLPipelineRecommender()  # NEW class
recommendations = recommender.recommend(
    profile,
    explain=True  # NEW: Get explanations
)

# All v2.0.0 code still works!
# New features are opt-in
```

### Recommended Upgrades:

1. **Enable ML recommendations**: Set `use_ml=True`
2. **Add quality checks**: Use `quality_check=True`
3. **Monitor resources**: Use `ResourceMonitor()`
4. **Try ensemble**: Use `EnsembleAnalyzer` for critical analyses
5. **Optimize parameters**: Use `ParameterOptimizer` for best results

---

## Performance Tips (NEW)

```python
# Tip 1: Use ML for faster recommendations
profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
# ML model is 10x faster than full benchmarking

# Tip 2: Monitor to identify bottlenecks
monitor = ResourceMonitor()
monitor.start()
# ... analysis ...
stats = monitor.stop()
print(f"Bottleneck: {stats['bottleneck']}")

# Tip 3: Optimize parameters once, reuse
if not os.path.exists('optimized_params.yaml'):
    optimizer = ParameterOptimizer()
    best_params = optimizer.optimize(counts, metadata)
    save_config(best_params, 'optimized_params.yaml')
else:
    best_params = load_config('optimized_params.yaml')

# Tip 4: Use dashboard for interactive exploration
# Much faster than re-running analyses
dashboard = DashboardServer(results_dir='results/')
dashboard.start()
```

---

## See Also

- [UNDERSTANDING_ML.md](UNDERSTANDING_ML.md) - ML concepts explained
- [ML_TRAINING.md](ML_TRAINING.md) - Train custom ML models
- [DASHBOARD.md](DASHBOARD.md) - Interactive dashboard guide
- [RESOURCE_MONITORING.md](RESOURCE_MONITORING.md) - Resource tracking
- [PARAMETER_OPTIMIZATION.md](PARAMETER_OPTIMIZATION.md) - Parameter tuning
- [QUALITY_ASSESSMENT.md](QUALITY_ASSESSMENT.md) - Data QC guide
- [ENSEMBLE.md](ENSEMBLE.md) - Ensemble analysis

---

**RAPTOR v2.1.0**  
**Author**: Ayeh Bolouki   
**License**: MIT

**Upgrade today for ML-powered RNA-seq analysis!** 🤖🦖
