# RAPTOR Profile & Recommend Guide

Get intelligent pipeline recommendations in 30 seconds!

## Quick Start
```bash
# Basic profiling
raptor profile --counts your_data.csv

# With metadata  
raptor profile --counts data.csv --metadata samples.csv

# Save HTML report
raptor profile --counts data.csv --output report.html
```

## How It Works

1. **Analyzes your data**: BCV, depth, zeros, library sizes
2. **Scores all 8 pipelines**: Accuracy, speed, memory, compatibility
3. **Recommends top 3**: With detailed reasoning
4. **Generates report**: Interactive HTML with visualizations

## Understanding Output

```
#1 RECOMMENDED: Pipeline 3 (Salmon-edgeR)
   Score: 0.88 | Confidence: High
   
   Why?
   ✓ Excellent balance for your data
   ✓ Handles medium BCV well (yours: 0.42)
   ✓ Fast runtime (~1hr for 12 samples)
   ✓ Low memory (8GB)
```

## Key Metrics Explained

**BCV (Biological Coefficient of Variation)**
- Low (<0.2): Cell lines, controlled
- Medium (0.2-0.6): Typical studies  
- High (>0.6): Clinical, complex

**Sequencing Depth**
- Low (<10M): May miss genes
- Medium (10-25M): Adequate
- High (>25M): Excellent

## Advanced Usage

```bash
# Prioritize accuracy
raptor profile --counts data.csv \
  --weight-accuracy 0.7 --weight-speed 0.1

# Resource constraints
raptor profile --counts data.csv \
  --max-memory 16G --max-runtime 4h

# Exclude pipelines
raptor profile --counts data.csv \
  --exclude-pipelines 7,8
```

## Python API
```python
from raptor import RNAseqDataProfiler, PipelineRecommender

# Profile
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# Recommend
recommender = PipelineRecommender()
recommendations = recommender.recommend(profile)
```

## Examples

**Standard DE study** (12 samples, 20M reads):  
→ Recommend: Pipeline 3 (Salmon-edgeR)

**Large cohort** (100 samples):  
→ Recommend: Pipeline 4 (Kallisto-Sleuth)

**Small pilot** (4 samples, need accuracy):  
→ Recommend: Pipeline 1 (STAR-RSEM-DESeq2)

See full examples and interpretation guide in complete documentation.

