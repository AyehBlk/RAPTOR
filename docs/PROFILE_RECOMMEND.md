# RAPTOR v2.1.0 Profile & Recommend Guide

Get intelligent ML-powered pipeline recommendations in seconds!

## Quick Start
```bash
# Basic profiling with ML (NEW!)
raptor profile --counts your_data.csv --use-ml

# With metadata and quality check (NEW!)
raptor profile \
  --counts data.csv \
  --metadata samples.csv \
  --use-ml \
  --quality-check

# Save interactive HTML report
raptor profile \
  --counts data.csv \
  --use-ml \
  --output report.html
```

## What's New in v2.1.0

⭐ **ML-Powered Recommendations** - Learn from 10,000+ past analyses  
⭐ **Confidence Scores** - Know how reliable recommendations are  
⭐ **Similar Projects** - See what worked for similar data  
⭐ **Quality Assessment** - Automatic data QC  
⭐ **Feature Importance** - Understand WHY each pipeline was recommended  
⭐ **Resource Prediction** - Estimate CPU/memory/time needs  

---

## How It Works

### v2.0.0 (Rule-Based)
```
1. Analyzes your data → BCV, depth, zeros, library sizes
2. Applies fixed rules → If BCV > 0.4 then...
3. Scores pipelines → Based on predetermined weights
4. Recommends top 3 → With reasoning
```

### v2.1.0 (ML-Enhanced)
```
1. Analyzes your data → BCV, depth, zeros, library sizes
2. Runs quality checks → Contamination, outliers, batch effects (NEW!)
3. ML model prediction → Based on 10,000+ similar analyses (NEW!)
4. Confidence scoring → How sure is the model? (NEW!)
5. Find similar projects → What worked before? (NEW!)
6. Recommends top 3 → With ML scores, confidence, and reasoning
```

---

## Understanding Output

### Traditional Output (v2.0.0)
```
#1 RECOMMENDED: Pipeline 3 (Salmon-edgeR)
   Score: 0.88 | Reasoning: Excellent balance for your data
   
   Why?
   ✓ Excellent balance for your data
   ✓ Handles medium BCV well (yours: 0.42)
   ✓ Fast runtime (~1hr for 12 samples)
   ✓ Low memory (8GB)
```

### NEW ML-Enhanced Output (v2.1.0)
```
#1 RECOMMENDED: Pipeline 3 (Salmon-edgeR)
   Rule-based Score: 0.88 | ML Score: 0.92 ⭐
   Confidence: HIGH (89%)
   
   Why this pipeline?
   ✓ Excellent balance for your data
   ✓ Handles medium BCV well (yours: 0.42)
   ✓ Fast runtime (~22 min estimated)
   ✓ Low memory (12GB peak predicted)
   
   ML Insights (NEW):
   🤖 Based on 1,247 similar analyses
   📊 87.3% historical success rate on similar data
   🎯 Top features: BCV (35%), Sample size (22%), Depth (18%)
   ✅ High confidence - This is a very reliable recommendation!
   
   Similar Successful Projects:
   • Project #8472: 12 samples, BCV=0.41, Human → Success ✓
   • Project #2391: 10 samples, BCV=0.38, Human → Success ✓
   • Project #5674: 14 samples, BCV=0.45, Mouse → Success ✓
   
   Resource Predictions:
   • Runtime: 22 minutes (±5 min)
   • Peak Memory: 12 GB (±2 GB)
   • CPU Usage: 85% average
   • Disk Space: 25 GB needed
```

---

## Key Metrics Explained

### Data Characteristics

**BCV (Biological Coefficient of Variation)**
- Low (<0.2): Cell lines, controlled conditions
- Medium (0.2-0.6): Typical experiments  
- High (>0.6): Clinical samples, complex biology
- **NEW**: ML learns which pipelines handle each BCV best

**Sequencing Depth**
- Low (<10M): May miss genes
- Medium (10-25M): Adequate
- High (>25M): Excellent
- **NEW**: ML predicts if depth is sufficient for each pipeline

**Quality Score (NEW v2.1.0)**
- 90-100: Excellent data quality
- 80-89: Good quality
- 70-79: Acceptable, minor issues
- <70: Poor quality, address issues first

---

## ML Confidence Levels (NEW)

### High Confidence (>80%)
```
🟢 HIGH CONFIDENCE

What it means:
✓ Your data is very similar to many past successful analyses
✓ ML model has strong pattern matches
✓ Historical success rate is high (>85%)
✓ Multiple similar projects agree on this pipeline

What you should do:
→ Trust this recommendation confidently
→ This pipeline is very likely to work well
→ No need to benchmark - go ahead and use it!

Example:
ML Score: 0.92
Confidence: 89%
Similar Projects: 1,247
Historical Success: 87.3%
```

### Medium Confidence (60-80%)
```
🟡 MEDIUM CONFIDENCE

What it means:
⚠️  Your data is somewhat similar to past analyses
⚠️  ML model found reasonable patterns
⚠️  Historical success rate is good (70-85%)
⚠️  Some similar projects, but not many

What you should do:
→ This is a good recommendation
→ Consider the top 2-3 recommendations
→ Review alternatives
→ Optional: Quick benchmark top 3 pipelines

Example:
ML Score: 0.78
Confidence: 72%
Similar Projects: 342
Historical Success: 78.2%
```

### Low Confidence (<60%)
```
🔴 LOW CONFIDENCE

What it means:
❌ Your data is unusual or novel
❌ ML model hasn't seen many similar examples
❌ Historical data is limited
❌ Predictions are uncertain

What you should do:
→ Don't rely solely on ML
→ Consider running benchmark
→ Review rule-based recommendation too
→ Your data might be genuinely unique

Example:
ML Score: 0.58
Confidence: 54%
Similar Projects: 23
Historical Success: 64.1%

Recommendation: Run benchmark on subset to validate
```

---

## Advanced Usage

### Prioritize Different Factors

```bash
# Prioritize accuracy over speed
raptor profile \
  --counts data.csv \
  --use-ml \
  --weight-accuracy 0.7 \
  --weight-speed 0.1

# Prioritize speed (quick results)
raptor profile \
  --counts data.csv \
  --use-ml \
  --weight-speed 0.6 \
  --weight-accuracy 0.2

# Balance all factors
raptor profile \
  --counts data.csv \
  --use-ml \
  --balanced-weights
```

### Resource Constraints

```bash
# Limited resources
raptor profile \
  --counts data.csv \
  --use-ml \
  --max-memory 16G \
  --max-runtime 2h

# Only get recommendations that fit these constraints
# ML will filter out resource-intensive pipelines
```

### Exclude Pipelines

```bash
# Don't recommend certain pipelines
raptor profile \
  --counts data.csv \
  --use-ml \
  --exclude-pipelines 7,8  # Exclude legacy pipelines

# Only consider specific pipelines
raptor profile \
  --counts data.csv \
  --use-ml \
  --only-pipelines 1,3,4  # Only top performers
```

### Explain ML Decision (NEW)

```bash
# Get detailed ML explanation
raptor profile \
  --counts data.csv \
  --use-ml \
  --explain \
  --verbose

# Output includes:
# - Feature importance scores
# - Similar project details
# - Model confidence breakdown
# - Alternative recommendations
```

---

## Quality Assessment Integration (NEW)

```bash
# Run comprehensive quality check
raptor profile \
  --counts data.csv \
  --metadata metadata.csv \
  --quality-check

# QC checks include:
# - Library size distribution
# - Gene detection rate
# - Contamination check
# - Batch effect detection
# - Outlier identification
```

**QC Output:**
```
Quality Assessment Summary
═══════════════════════════

Overall Score: 87/100 🟢 GOOD

Detailed Metrics:
├─ Library Sizes:     92/100 ✅ (CV = 12%)
├─ Gene Detection:    90/100 ✅ (18,234 genes)
├─ Zero Inflation:    88/100 ✅ (42%)
├─ Batch Effects:     75/100 ⚠️  (Minor detected)
├─ Outliers:          80/100 ⚠️  (1 flagged: Sample_7)
└─ Contamination:    100/100 ✅ (None detected)

Issues Found:
⚠️  Minor batch effect detected (PVR = 8%)
   → Recommendation: Include batch in model

⚠️  Sample_7 is potential outlier
   → Recommendation: Review or exclude

✅ Data quality is good for analysis
```

---

## Python API

### Basic Usage
```python
from raptor import RNAseqDataProfiler, MLPipelineRecommender
import pandas as pd

# Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Profile with ML
profiler = RNAseqDataProfiler(
    counts, 
    metadata,
    use_ml=True  # Enable ML (NEW)
)
profile = profiler.profile(quality_check=True)  # Run QC (NEW)

# Get ML recommendations
ml_recommender = MLPipelineRecommender()
recommendations = ml_recommender.recommend(
    profile, 
    n=3,
    explain=True  # Get feature importance (NEW)
)

# Access results
top = recommendations[0]
print(f"Pipeline: {top['pipeline_name']}")
print(f"ML Score: {top['ml_score']:.3f}")
print(f"Confidence: {top['confidence']}")
print(f"Success Rate: {top['historical_success']:.1%}")
```

### Get Similar Projects (NEW)
```python
# Find similar past analyses
similar = ml_recommender.get_similar_projects(profile, n=10)

for project in similar:
    print(f"Project {project['id']}")
    print(f"  Similarity: {project['similarity']:.3f}")
    print(f"  Pipeline used: {project['pipeline']}")
    print(f"  Success: {project['success']}")
```

### Explain Recommendation (NEW)
```python
# Get feature importance
explanation = ml_recommender.explain_recommendation(
    profile, 
    pipeline_id=3
)

print("Feature Importance:")
for feature, importance in explanation['features']:
    print(f"  {feature}: {importance:.3f}")

print("\nWhy this works for your data:")
for reason in explanation['reasons']:
    print(f"  • {reason}")
```

---

## Understanding ML Feature Importance (NEW)

When ML recommends a pipeline, it shows which features mattered most:

```
Feature Importance for Pipeline 3 (Salmon-edgeR)
═════════════════════════════════════════════════

1. BCV (0.348) ████████████████████████████████████
   → Your BCV (0.42) is perfect for this pipeline
   
2. Sample Size (0.221) ██████████████████████
   → Your n=12 is in the sweet spot
   
3. Sequencing Depth (0.182) ██████████████████
   → Your 25M reads/sample is excellent
   
4. Zero Inflation (0.127) ████████████
   → Your 42% is normal and well-handled
   
5. Library Size CV (0.086) ████████
   → Your 12% variation is manageable
   
6. Organism (0.036) ███
   → Human data, commonly analyzed

Total explained: 100%

Key Insight:
Your BCV and sample size are the main reasons this 
pipeline was recommended. They match the optimal 
range for Salmon-edgeR based on 1,247 similar projects.
```

---

## Examples

### Standard DE Study (12 samples, 20M reads)

**Input:**
```bash
raptor profile --counts data.csv --use-ml
```

**Output:**
```
🥇 Pipeline 3 (Salmon-edgeR)
   ML Score: 0.92 | Confidence: HIGH (89%)
   
   Why: Perfect for your data
   - Medium BCV (0.42) → Optimal range
   - 12 samples → Ideal for this pipeline
   - 25M depth → Excellent
   
   Success Rate: 87.3% (on 1,247 similar analyses)
   Runtime: ~22 min | Memory: ~12 GB

→ Recommendation: Use Pipeline 3 with confidence!
```

---

### Large Cohort (100 samples)

**Input:**
```bash
raptor profile --counts large_data.csv --use-ml
```

**Output:**
```
🥇 Pipeline 4 (Kallisto-Sleuth)
   ML Score: 0.90 | Confidence: HIGH (88%)
   
   Why: Optimized for large datasets
   - 100 samples → Kallisto scales well
   - Fast runtime critical for this size
   - Low memory usage (6GB vs 45GB for STAR)
   
   Success Rate: 90.1% (on 234 large studies)
   Runtime: ~1 hour | Memory: ~6 GB
   
   Alternative: Pipeline 3 (Salmon-edgeR)
   Runtime: ~4 hours | Memory: ~15 GB
   More accurate but slower

→ Recommendation: Use Pipeline 4 for speed, or 
   Pipeline 3 if you can wait and want slightly 
   better accuracy.
```

---

### Small Pilot (4 samples, need accuracy)

**Input:**
```bash
raptor profile --counts pilot.csv --use-ml --weight-accuracy 0.8
```

**Output:**
```
🥇 Pipeline 1 (STAR-RSEM-DESeq2)
   ML Score: 0.88 | Confidence: HIGH (85%)
   
   Why: Best for small samples
   - Only 4 samples → Need robust method
   - DESeq2 handles small n well
   - Accuracy prioritized (your request)
   
   Success Rate: 89.7% (on small sample studies)
   Runtime: ~3 hours | Memory: ~45 GB
   
   ⚠️  Warning: Low Confidence for Pipeline 3
   Pipeline 3 (Salmon-edgeR) scored 0.76 but with
   LOW confidence (62%) for n=4. Not recommended.

→ Recommendation: Use Pipeline 1. Worth the extra time!
```

---

### Novel Organism (Rare fish species)

**Input:**
```bash
raptor profile --counts rare_fish.csv --use-ml
```

**Output:**
```
🔴 LOW CONFIDENCE WARNING

Top Recommendation:
🥇 Pipeline 3 (Salmon-edgeR)
   ML Score: 0.64 | Confidence: LOW (51%) ⚠️
   
   Why confidence is low:
   ❌ Only 8 similar non-model organisms in training
   ❌ No fish species in training data
   ❌ ML is essentially guessing
   
   Rule-based Score: 0.82 (More reliable here)
   
→ Recommendation: 
   Your data is genuinely novel. ML can't help much.
   
   Suggested approach:
   1. Run quick benchmark (Pipelines 1,3,4)
   2. Compare results
   3. Choose based on concordance
   
   Or trust rule-based recommendation (Pipeline 3)
   with the understanding that it's based on logic,
   not learned patterns.
```

---

## Troubleshooting

### Problem: All recommendations have low confidence

**This means:** Your data is unusual/novel

**Solutions:**
1. Check if data characteristics are extreme
2. Run quality check to identify issues
3. Consider benchmarking instead of single pipeline
4. Trust rule-based recommendations
5. Train custom ML model on your lab's data

---

### Problem: ML disagrees with rule-based recommendation

**Example:**
```
Rule-based: Pipeline 1 (Score: 0.82)
ML: Pipeline 3 (ML Score: 0.88, Confidence: HIGH)
```

**This means:** ML found patterns rules didn't catch

**What to do:**
- If ML confidence is HIGH → Trust ML
- If ML confidence is MEDIUM → Consider both
- If ML confidence is LOW → Trust rules
- When in doubt → Quick benchmark both

---

### Problem: Quality score is low

**Solutions:**
```bash
# Get detailed quality report
raptor profile \
  --counts data.csv \
  --quality-check \
  --verbose

# Address identified issues:
# - Remove outlier samples
# - Filter low-quality genes
# - Account for batch effects
# - Check for contamination

# Then re-profile
raptor profile --counts cleaned_data.csv --use-ml
```

---

## Best Practices

### Do's:
✅ Always use `--use-ml` for v2.1.0  
✅ Check ML confidence scores  
✅ Review quality assessment  
✅ Read why ML recommended each pipeline  
✅ Consider top 2-3 recommendations  
✅ Look at similar successful projects  

### Don'ts:
❌ Blindly follow low-confidence recommendations  
❌ Ignore quality warnings  
❌ Skip quality checks for critical analyses  
❌ Forget that ML is a tool, not a dictator  
❌ Ignore your domain expertise  

---

## See Also

- [UNDERSTANDING_ML.md](UNDERSTANDING_ML.md) - ML concepts explained
- [ML_TRAINING.md](ML_TRAINING.md) - Advanced ML usage
- [QUALITY_ASSESSMENT.md](QUALITY_ASSESSMENT.md) - Data QC guide
- [tutorial_04_ml_recommendations.md](tutorials/tutorial_04_ml_recommendations.md) - Hands-on tutorial
- [PIPELINES.md](PIPELINES.md) - Detailed pipeline information
- [DASHBOARD.md](DASHBOARD.md) - Interactive visualization

---

**RAPTOR v2.1.0 Profile & Recommend**  
**Author**: Ayeh Bolouki  
**Affiliation**: University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**License**: MIT

*"Profile once, recommend smartly!"* 🎯🤖🦖
