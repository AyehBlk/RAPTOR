# RAPTOR Python API Reference

Complete API documentation for programmatic use.

## Installation
```python
pip install raptor-rnaseq
import raptor
```

## Core Classes

### RNAseqDataProfiler

Profile RNA-seq count data for pipeline recommendation.

```python
from raptor import RNAseqDataProfiler
import pandas as pd

# Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Create profiler
profiler = RNAseqDataProfiler(counts, metadata)

# Run profiling
profile = profiler.profile()

# Access results
print(f"BCV: {profile['bcv']:.3f}")
print(f"Mean depth: {profile['mean_depth']:.0f}")
print(f"Samples: {profile['n_samples']}")
```

**Methods:**
- `profile()`: Run complete profiling, returns dict
- `calculate_bcv()`: Calculate biological variation
- `assess_depth()`: Analyze sequencing depth
- `detect_outliers()`: Identify problematic samples
- `generate_plots()`: Create diagnostic plots

**Profile Output:**
```python
{
    'n_samples': 12,
    'n_genes': 15234,
    'bcv': 0.42,
    'bcv_category': 'medium',
    'mean_depth': 25000000,
    'depth_category': 'high',
    'zero_inflation': 0.15,
    'library_size_cv': 0.18,
    'outliers': [],
    'quality_flags': []
}
```

### PipelineRecommender

Recommend optimal pipelines based on data profile.

```python
from raptor import PipelineRecommender

recommender = PipelineRecommender()

# Get recommendations
recommendations = recommender.recommend(profile)

# Top recommendation
top = recommendations[0]
print(f"Pipeline: {top['pipeline_name']}")
print(f"Score: {top['score']:.2f}")
print(f"Reasoning: {top['reasoning']}")
```

**Methods:**
- `recommend(profile, n=3)`: Get top N recommendations
- `score_pipeline(pipeline_id, profile)`: Score specific pipeline
- `compare_pipelines(profile)`: Compare all pipelines

**Recommendation Output:**
```python
[
    {
        'pipeline_id': 3,
        'pipeline_name': 'Salmon-edgeR',
        'score': 0.88,
        'confidence': 'high',
        'reasoning': 'Excellent balance...',
        'expected_runtime': 1.2,
        'expected_memory': 10,
        'expected_accuracy': 0.90
    },
    ...
]
```

### PipelineBenchmark

Run and benchmark multiple pipelines.

```python
from raptor import PipelineBenchmark

benchmark = PipelineBenchmark(
    data_dir='fastq/',
    output_dir='results/',
    threads=8,
    memory='32G'
)

# Run pipelines
results = benchmark.run_pipelines([1, 3, 4])

# Save results
benchmark.save_results(results)
```

**Methods:**
- `run_pipelines(pipeline_ids)`: Run specified pipelines
- `run_single_pipeline(pipeline_id)`: Run one pipeline
- `save_results(results)`: Save to JSON
- `load_results()`: Load from JSON

### DataSimulator

Generate simulated RNA-seq data.

```python
from raptor import DataSimulator

simulator = DataSimulator(
    n_genes=2000,
    n_samples=6,
    n_de=400,
    fold_changes=[0.5, 2.0],
    seed=42
)

summary = simulator.generate_data('sim_data/')
```

**Methods:**
- `generate_data(output_dir)`: Create simulated data
- Quick presets: `quick_simulate(output_dir, size='small')`

### ReportGenerator

Generate HTML/PDF reports.

```python
from raptor import ReportGenerator

generator = ReportGenerator()

# Profile report
generator.generate_profile_report(
    profile,
    recommendations,
    output='report.html'
)

# Benchmark report  
generator.generate_benchmark_report(
    results,
    output='comparison.html'
)
```

**Methods:**
- `generate_profile_report()`: Profiling + recommendations
- `generate_benchmark_report()`: Pipeline comparisons
- `generate_full_report()`: Complete analysis

## Utility Functions

```python
from raptor.utils import *

# File operations
ensure_dir('output/')
check_file_exists('data.csv')

# System checks
check_command_exists('STAR')
check_required_tools(['STAR', 'salmon', 'kallisto'])
get_available_memory()  # Returns GB
get_cpu_count()  # Returns cores

# Data validation
validate_count_matrix(counts)
validate_metadata(metadata, counts)

# Configuration
config = load_config('config.yaml')
save_config(config, 'my_config.yaml')

# Formatting
format_time(3665)  # "1h 1m 5s"
format_bytes(1536000)  # "1.5 MB"
```

## Complete Workflow Example

```python
import pandas as pd
from raptor import (
    RNAseqDataProfiler,
    PipelineRecommender,
    PipelineBenchmark,
    ReportGenerator
)

# 1. Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# 2. Profile data
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# 3. Get recommendations
recommender = PipelineRecommender()
recommendations = recommender.recommend(profile, n=3)

# 4. Generate report
report_gen = ReportGenerator()
report_gen.generate_profile_report(
    profile,
    recommendations,
    output='my_recommendation.html'
)

# 5. Optional: Run benchmark
top_3 = [r['pipeline_id'] for r in recommendations]
benchmark = PipelineBenchmark(
    data_dir='fastq/',
    output_dir='benchmark_results/'
)
results = benchmark.run_pipelines(top_3)

# 6. Generate comparison report
report_gen.generate_benchmark_report(
    results,
    output='my_comparison.html'
)
```

## Configuration

Load and modify configuration:

```python
from raptor.utils import load_config, save_config

# Load config
config = load_config('config/config.yaml')

# Modify settings
config['resources']['default_threads'] = 16
config['statistics']['fdr_threshold'] = 0.01

# Save
save_config(config, 'my_custom_config.yaml')

# Use in profiler
profiler = RNAseqDataProfiler(
    counts, 
    metadata,
    config_file='my_custom_config.yaml'
)
```

## Error Handling

```python
from raptor import RNAseqDataProfiler
from raptor.utils import validate_count_matrix

try:
    # Validate first
    validate_count_matrix(counts)
    
    # Then profile
    profiler = RNAseqDataProfiler(counts, metadata)
    profile = profiler.profile()
    
except FileNotFoundError as e:
    print(f"File not found: {e}")
except ValueError as e:
    print(f"Invalid data: {e}")
except Exception as e:
    print(f"Unexpected error: {e}")
```

## Logging

```python
from raptor.utils import setup_logging
import logging

# Setup logging
setup_logging(level='DEBUG', log_file='raptor.log')

# Use logger
logger = logging.getLogger('raptor')
logger.info("Starting analysis...")
```

## Version Info

```python
import raptor

print(raptor.__version__)
print(raptor.__author__)
print(raptor.__email__)

# Get environment info
from raptor.utils import get_environment_info
info = get_environment_info()
print(info)
```

## Type Hints

All classes use type hints:

```python
from typing import Dict, List
from pathlib import Path

def my_analysis(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    output_dir: Path
) -> Dict[str, float]:
    profiler = RNAseqDataProfiler(counts, metadata)
    return profiler.profile()
```

See source code for complete type annotations.

