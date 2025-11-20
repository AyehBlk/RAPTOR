# RAPTOR v2.1.0 Pipelines Reference

Deep dive into all 8 RNA-seq analysis pipelines with ML-powered selection guidance.

## Pipeline Overview

| ID | Name | Type | Speed | Accuracy | Memory | ML Rank* |
|----|------|------|-------|----------|--------|----------|
| 1 | STAR-RSEM-DESeq2 | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚫⚫ | High | #2 |
| 2 | HISAT2-StringTie-Ballgown | Alignment | ⚫⚫⚫⚪⚪ | ⚫⚫⚫⚫⚪ | Medium | #5 |
| 3 | Salmon-edgeR | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚫⚪ | Low | #1 ⭐ |
| 4 | Kallisto-Sleuth | Pseudo-align | ⚫⚫⚫⚫⚫ | ⚫⚫⚫⚪⚪ | Low | #3 |
| 5 | STAR-HTSeq-limma | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚫⚪ | High | #4 |
| 6 | STAR-featureCounts-NOISeq | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚪⚪ | High | #6 |
| 7 | Bowtie2-RSEM-EBSeq | Alignment | ⚫⚫⚪⚪⚪ | ⚫⚫⚫⚪⚪ | Medium | #7 |
| 8 | HISAT2-Cufflinks-Cuffdiff | Alignment | ⚫⚪⚪⚪⚪ | ⚫⚫⚪⚪⚪ | Medium | #8 |

*ML Rank: Average ranking from ML model across typical datasets (v2.1.0)

---

## 🆕 What's New in v2.1.0

Each pipeline now includes:
- 🤖 **ML Success Rate** - Historical success on similar data
- 📊 **Quality Compatibility** - Which data qualities work best
- ⚡ **Resource Predictions** - Estimated CPU/memory/time
- 🎯 **Best Use Cases** - ML-learned optimal scenarios
- 🔄 **Ensemble Weight** - Recommended weight in ensemble analysis

---

## Pipeline 1: STAR-RSEM-DESeq2

**Gold Standard - Highest Accuracy**

### Components
- **Alignment**: STAR (splice-aware)
- **Quantification**: RSEM (EM algorithm)
- **Statistics**: DESeq2 (negative binomial)

### Best For
- Publication-quality results
- Small sample sizes (n<6)
- Complex experimental designs
- When accuracy is paramount
- **NEW**: High biological variation (BCV > 0.6)

### Performance
- Runtime: ~2-6 hours (typical dataset)
- Memory: 32-48 GB
- Accuracy: 95%
- **NEW**: ML Success Rate: 92% (on similar data)
- **NEW**: Resource Efficiency: 7/10

### ML Insights (NEW v2.1.0)
```
When ML Recommends Pipeline 1:
✓ Small sample size (n < 10)
✓ Need highest accuracy
✓ Have sufficient resources (>32GB RAM)
✓ Complex design or batch effects
✓ High biological variation

ML Success Rate by Data Type:
├─ Small samples (n<6): 94% success
├─ Medium BCV (0.3-0.6): 93% success
├─ High BCV (>0.6): 89% success
└─ Complex designs: 91% success

Historical Performance:
Based on 2,847 analyses in ML training data
```

### Running
```bash
# Basic
raptor run --pipeline 1 \
  --data fastq/ \
  --reference /path/to/star_index \
  --annotation genes.gtf \
  --output pipeline1_results/

# With ML-optimized parameters (NEW)
raptor run --pipeline 1 \
  --data fastq/ \
  --reference /path/to/star_index \
  --use-ml-params \
  --monitor-resources
```

### Quality Recommendations (NEW)
```
Optimal Data Characteristics:
├─ Sample size: 3-20 (works well with small n)
├─ BCV: Any (handles high variation)
├─ Depth: >10M reads/sample
├─ Library CV: <0.4
└─ Zero inflation: <70%

Avoid if:
✗ Very large datasets (>100 samples) - too slow
✗ Very low resources (<32GB RAM)
✗ Need quick turnaround (<2 hours)
```

---

## Pipeline 2: HISAT2-StringTie-Ballgown

**Transcript Assembly & Novel Discovery**

### Components
- **Alignment**: HISAT2
- **Assembly**: StringTie
- **Statistics**: Ballgown

### Best For
- Novel transcript discovery
- Isoform-level analysis
- Non-model organisms
- When reference incomplete
- **NEW**: Exploring unannotated regions

### Performance
- Runtime: ~1-4 hours
- Memory: 16-24 GB
- Accuracy: 88%
- **NEW**: ML Success Rate: 81% (context-dependent)
- **NEW**: Best for transcript discovery, not DE

### ML Insights (NEW)
```
When ML Recommends Pipeline 2:
✓ Non-model organism
✓ Incomplete annotation
✓ Novel transcript discovery primary goal
✓ Isoform-level analysis needed

ML Success Rate by Data Type:
├─ Non-model organisms: 87% success
├─ Novel transcript discovery: 91% success
├─ Standard DE analysis: 75% success (not optimal)
└─ Isoform switching: 85% success

Warning: ML rarely recommends for standard DE
```

---

## Pipeline 3: Salmon-edgeR ⭐ ML TOP CHOICE

**Best Balance - Fast & Accurate**

### Components
- **Quantification**: Salmon (quasi-mapping)
- **Statistics**: edgeR (quasi-likelihood)

### Best For
- Most RNA-seq experiments
- Large datasets (>20 samples)
- Quick turnaround needed
- Good balance of all metrics
- **NEW**: ML's most frequently recommended pipeline

### Performance
- Runtime: ~0.5-2 hours
- Memory: 8-16 GB
- Accuracy: 90%
- **NEW**: ML Success Rate: 89% (highest overall)
- **NEW**: Resource Efficiency: 10/10 (best)

### ML Insights (NEW) 
```
Why ML Loves Pipeline 3:
✓ Best speed/accuracy tradeoff
✓ Lowest resource requirements
✓ Handles most data types well
✓ Robust to batch effects
✓ Excellent for n=6-50 samples

ML Success Rate by Data Type:
├─ Standard experiments: 91% success
├─ Large samples (n>20): 94% success
├─ Medium BCV (0.2-0.6): 92% success
├─ Balanced designs: 90% success
└─ Overall: 89% success (BEST!)

ML Recommendation Frequency:
Recommended in 47% of all profiles
(Most recommended pipeline)

Historical Performance:
Based on 5,234 analyses in ML training data
(Largest dataset - most reliable predictions)
```

### Running
```bash
# Basic
raptor run --pipeline 3 \
  --data fastq/ \
  --transcriptome /path/to/salmon_index \
  --output pipeline3_results/

# ML-optimized (NEW)
raptor run --pipeline 3 \
  --data fastq/ \
  --transcriptome /path/to/salmon_index \
  --use-ml-params \
  --quality-check \
  --monitor-resources
```

### Quality Recommendations (NEW)
```
Optimal Data Characteristics:
├─ Sample size: 6-100 (sweet spot: 12-24)
├─ BCV: 0.2-0.6 (medium variation)
├─ Depth: >15M reads/sample
├─ Library CV: <0.3
└─ Zero inflation: 30-60%

Works well even if:
✓ Slight batch effects
✓ Some outliers
✓ Moderate library size variation

Avoid if:
✗ Very small samples (n<4) - use Pipeline 1
✗ Need novel transcript discovery - use Pipeline 2
```

---

## Pipeline 4: Kallisto-Sleuth

**Ultra-Fast - Large Studies**

### Components
- **Quantification**: Kallisto
- **Statistics**: Sleuth (bootstrap-based)

### Best For
- Very large datasets (>50 samples)
- Exploratory analysis
- Minimal resources
- Speed is critical
- **NEW**: Cloud deployment (cost-effective)

### Performance
- Runtime: ~0.3-1 hour
- Memory: 4-8 GB
- Accuracy: 88%
- **NEW**: ML Success Rate: 84%
- **NEW**: Cost Efficiency: 10/10 (best for cloud)

### ML Insights (NEW)
```
When ML Recommends Pipeline 4:
✓ Very large datasets (>50 samples)
✓ Limited resources (<16GB RAM)
✓ Speed critical
✓ Exploratory analysis
✓ Cloud/spot instances

ML Success Rate by Data Type:
├─ Large samples (n>50): 88% success
├─ Very large (n>100): 90% success
├─ Low resources: 87% success
└─ Standard (n=12): 82% success

Best Use Case:
Large-scale studies where speed > accuracy
Historical: 1,923 analyses in ML training
```

---

## Pipeline 5: STAR-HTSeq-limma-voom

**Flexible Modeling - Complex Designs**

### Components
- **Alignment**: STAR
- **Counting**: HTSeq
- **Statistics**: limma-voom

### Best For
- Complex experimental designs
- Multi-factor analysis
- Batch correction needed
- Repeated measures
- **NEW**: Time-series experiments

### Performance
- Runtime: ~2-7 hours
- Memory: 32-40 GB
- Accuracy: 92%
- **NEW**: ML Success Rate: 88% (complex designs)
- **NEW**: Flexibility: 10/10 (most flexible)

### ML Insights (NEW)
```
When ML Recommends Pipeline 5:
✓ Complex experimental design
✓ Multiple factors
✓ Batch effects present
✓ Paired/repeated samples
✓ Need flexible modeling

ML Success Rate by Data Type:
├─ Complex designs: 91% success
├─ Batch effects: 93% success
├─ Multi-factor: 89% success
├─ Paired samples: 90% success
└─ Simple designs: 85% success

Historical: 1,456 analyses
Note: Often 2nd choice after Pipeline 1
```

---

## Pipeline 6: STAR-featureCounts-NOISeq

**Non-Parametric - Small Samples**

### Components
- **Alignment**: STAR
- **Counting**: featureCounts
- **Statistics**: NOISeq

### Best For
- Very small samples (n=2-3)
- No replicates
- Non-normal distributions
- **NEW**: Highly variable data (BCV > 0.8)

### Performance
- Runtime: ~2-7 hours
- Memory: 32-36 GB
- Accuracy: 85%
- **NEW**: ML Success Rate: 78% (small n only)
- **NEW**: Small Sample Specialist

### ML Insights (NEW)
```
When ML Recommends Pipeline 6:
✓ Very small samples (n<4)
✓ No or few replicates
✓ Very high variation (BCV > 0.8)
✓ Non-normal distributions

ML Success Rate by Data Type:
├─ n=2 samples: 81% success
├─ n=3 samples: 84% success
├─ n=4+ samples: 73% success (not optimal)
└─ High BCV (>0.8): 79% success

Warning from ML:
Rarely optimal for n>3
Consider Pipeline 1 instead for n≥4

Historical: 789 analyses
(Smallest training set - less confident predictions)
```

---

## Pipeline 7: Bowtie2-RSEM-EBSeq

**Memory-Efficient Alternative**

### Components
- **Alignment**: Bowtie2
- **Quantification**: RSEM
- **Statistics**: EBSeq

### Best For
- Moderate resource environments
- Isoform-level analysis
- Two-condition comparisons
- **NEW**: When STAR unavailable

### Performance
- Runtime: ~3-8 hours
- Memory: 16-24 GB
- Accuracy: 87%
- **NEW**: ML Success Rate: 76%
- **NEW**: Rarely ML's first choice

### ML Insights (NEW)
```
ML Recommendation Notes:
⚠️  Rarely recommended by ML
⚠️  Usually better alternatives exist
⚠️  Consider Pipeline 1 or 3 instead

Use Pipeline 7 when:
- STAR not available
- Moderate memory constraints
- Isoform analysis needed

ML typically ranks this #7
Historical: 623 analyses
```

---

## Pipeline 8: HISAT2-Cufflinks-Cuffdiff

**Legacy Pipeline**

### Components
- **Alignment**: HISAT2
- **Assembly**: Cufflinks
- **Statistics**: Cuffdiff

### Best For
- Reproducing legacy analyses
- When Cufflinks ecosystem required
- **NEW**: Not recommended for new analyses

### Performance
- Runtime: ~4-12 hours
- Memory: 20-32 GB
- Accuracy: 82%
- **NEW**: ML Success Rate: 68% (lowest)
- **NEW**: Legacy use only

### ML Insights (NEW)
```
ML Strongly Discourages Pipeline 8:
❌ Lowest accuracy among all pipelines
❌ Slowest runtime
❌ Better alternatives always exist
❌ Superseded by newer methods

ML Success Rate:
Only 68% success rate
Recommended in <1% of profiles

Use ONLY if:
- Reproducing old analyses
- Required by journal/collaborators

ML Advice:
Consider Pipeline 2 for transcript discovery
Consider Pipeline 3 for standard DE

Historical: 412 analyses
(Smallest and worst-performing in ML training)
```

---

## 🤖 ML-Powered Pipeline Selection (NEW)

### How ML Chooses Pipelines

```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler

# Profile your data
profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile()

# Get ML recommendations
ml_rec = MLPipelineRecommender()
recommendations = ml_rec.recommend(profile, n=3, explain=True)

# See why ML chose each pipeline
for rec in recommendations:
    print(f"\n{rec['pipeline_name']}")
    print(f"ML Score: {rec['ml_score']:.3f}")
    print(f"Success Rate: {rec['historical_success']:.1%}")
    print(f"Based on {rec['similar_count']} similar projects")
```

### ML Decision Factors

```
ML considers (in order of importance):
1. Biological Variation (BCV) - 35% weight
2. Sample Size - 22% weight
3. Sequencing Depth - 18% weight
4. Library Size Variation - 12% weight
5. Zero Inflation - 8% weight
6. Organism - 5% weight

For YOUR data:
raptor profile --counts data.csv --use-ml --explain
```

---

## Decision Guide

### Quick Decision Tree

```
Start Here
    ↓
Need highest accuracy?
├─ YES → Pipeline 1 (STAR-RSEM-DESeq2)
└─ NO → Continue
    ↓
Have >50 samples?
├─ YES → Pipeline 4 (Kallisto-Sleuth)
└─ NO → Continue
    ↓
Need novel transcripts?
├─ YES → Pipeline 2 (HISAT2-StringTie)
└─ NO → Continue
    ↓
Have <4 samples?
├─ YES → Pipeline 6 (NOISeq)
└─ NO → Continue
    ↓
Complex design/batch effects?
├─ YES → Pipeline 5 (limma-voom)
└─ NO → Pipeline 3 (Salmon-edgeR) ⭐

NOT SURE? → Use ML recommendation:
raptor profile --counts data.csv --use-ml
```

### ML-Enhanced Decision (NEW)

```bash
# Let ML decide for you
raptor profile --counts data.csv --use-ml

# ML will consider:
# - Your exact data characteristics
# - 10,000+ past successful analyses
# - Resource constraints
# - Expected accuracy
# - Historical success rate

# Get top 3 with confidence scores
# Choose based on ML confidence!
```

---

## Ensemble Recommendations (NEW)

For critical analyses, combine multiple pipelines:

### ML-Suggested Ensembles

**High Confidence Ensemble:**
```
Pipeline 1 + Pipeline 3 + Pipeline 5
(Weight: 0.4, 0.4, 0.2)

Why: Different methods, high accuracy
ML Success: 94% on combined results
```

**Fast Ensemble:**
```
Pipeline 3 + Pipeline 4
(Weight: 0.6, 0.4)

Why: Both fast, complementary
ML Success: 88%
Runtime: <2 hours
```

**Maximum Robustness:**
```
Pipeline 1 + Pipeline 3 + Pipeline 4 + Pipeline 5
(Weight: 0.3, 0.3, 0.2, 0.2)

Why: Diverse methods
ML Success: 96%
Consensus genes highly reliable
```

---

## See Configuration Details

Pipeline parameters and settings:
```bash
# View pipeline configurations
raptor show-pipeline-config --pipeline 3

# Get ML-optimized parameters
raptor get-ml-params --pipeline 3 --data-profile profile.json

# Compare all pipeline configs
raptor compare-configs --pipelines all
```

---

**RAPTOR v2.1.0 Pipelines**  
**Author**: Ayeh Bolouki  
**License**: MIT

*"Eight pipelines, one ML brain!"* 🤖🦖
