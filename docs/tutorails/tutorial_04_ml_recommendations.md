# Tutorial 4: ML-Powered Pipeline Recommendations

**Level**: Intermediate  
**Time**: 45 minutes  
**Goal**: Master machine learning-powered pipeline selection

---

## What You'll Learn

- How RAPTOR's ML recommendation system works
- How to interpret ML confidence scores
- When to trust ML vs rule-based recommendations
- How to customize ML recommendations
- Training custom models on your data
- **NEW in v2.1.1**: Combining ML recommendations with threshold optimization

---

## 🆕 What's New in v2.1.1

**ML + ATO Integration**: After getting ML pipeline recommendations, optimize your thresholds:

```bash
# 1. Get ML recommendation
raptor profile --counts data.csv --recommend

# 2. Run recommended pipeline
raptor run --pipeline RECOMMENDED --counts data.csv --output results/

# 3. Optimize thresholds with ATO
raptor optimize-thresholds --results results/degs.csv --goal balanced
```

This combines the best of both systems: ML picks the optimal pipeline, ATO determines the optimal thresholds.

- Completed [Tutorial 1: Getting Started](tutorial_01_getting_started.md)
- RAPTOR v2.1.1+ installed
- Python programming basics
- Understanding of RNA-seq analysis

---

## Why ML Recommendations?

### Traditional Approach (v2.0.0)
```
Your Data Characteristics
         ↓
    Rule-Based Logic
    (if-then rules)
         ↓
    Recommendation
```

**Problems:**
- Fixed rules can't capture complex patterns
- Doesn't learn from past successes
- Can't adapt to new scenarios

### ML Approach (v2.1.1)
```
Your Data Characteristics
         ↓
 ML Model (trained on 10,000+ analyses)
         ↓
  Prediction + Confidence + Reasoning
```

**Benefits:**
- ✅ Learns from real-world data
- ✅ Captures complex patterns
- ✅ Provides confidence scores
- ✅ Shows similar past projects
- ✅ Improves over time

---

## Step 1: Your First ML Recommendation

Let's compare rule-based vs ML recommendations:

```python
import pandas as pd
from raptor import RNAseqDataProfiler, PipelineRecommender, MLPipelineRecommender

# Load data
counts = pd.read_csv('data/counts.csv', index_col=0)
metadata = pd.read_csv('data/metadata.csv')

# Profile data
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# 1. Traditional rule-based recommendation
traditional_rec = PipelineRecommender()
trad_results = traditional_rec.recommend(profile, n=3)

print("=== TRADITIONAL (Rule-Based) ===")
for i, rec in enumerate(trad_results, 1):
    print(f"\n#{i} {rec['pipeline_name']}")
    print(f"   Score: {rec['score']:.3f}")
    print(f"   Reasoning: {rec['reasoning']}")

# 2. ML-powered recommendation (NEW!)
ml_rec = MLPipelineRecommender()
ml_results = ml_rec.recommend(profile, n=3, explain=True)

print("\n\n=== ML-POWERED ===")
for i, rec in enumerate(ml_results, 1):
    print(f"\n#{i} {rec['pipeline_name']}")
    print(f"   ML Score: {rec['ml_score']:.3f}")
    print(f"   Confidence: {rec['confidence']}")  # high/medium/low
    print(f"   Based on: {rec['similar_count']} similar projects")
    print(f"   Success rate: {rec['historical_success']:.1%}")
    print(f"   Top features: {rec['top_features'][:3]}")
```

**Output Example:**
```
=== TRADITIONAL (Rule-Based) ===

#1 Salmon-edgeR
   Score: 0.880
   Reasoning: Excellent balance for your data...

=== ML-POWERED ===

#1 Salmon-edgeR
   ML Score: 0.892
   Confidence: high
   Based on: 1,247 similar projects
   Success rate: 87.3%
   Top features: [('bcv', 0.35), ('n_samples', 0.22), ('depth', 0.18)]
```

---

## Step 2: Understanding ML Confidence

Confidence tells you how sure the ML model is:

```python
from raptor import MLPipelineRecommender

ml_rec = MLPipelineRecommender()
results = ml_rec.recommend(profile, n=3, explain=True)

for rec in results:
    confidence = rec['confidence']
    similar_count = rec['similar_count']
    
    print(f"\nPipeline: {rec['pipeline_name']}")
    print(f"Confidence: {confidence}")
    
    if confidence == 'high':
        print("✅ RELIABLE: Model has seen many similar datasets")
        print(f"   Based on {similar_count} similar projects")
        print("   → Safe to follow this recommendation")
        
    elif confidence == 'medium':
        print("⚠️  MODERATE: Some similar datasets in training")
        print(f"   Based on {similar_count} similar projects")
        print("   → Good recommendation, consider alternatives")
        
    else:  # low
        print("🔴 UNCERTAIN: Few similar datasets in training")
        print(f"   Based on only {similar_count} similar projects")
        print("   → Consider running benchmark to validate")
```

### What Affects Confidence?

```
HIGH Confidence (>0.8):
├─ Your data is similar to many training examples
├─ Model has strong pattern matches
└─ Recommendations are very reliable

MEDIUM Confidence (0.6-0.8):
├─ Some similar examples in training data
├─ Model found reasonable patterns
└─ Recommendations are good but verify

LOW Confidence (<0.6):
├─ Your data is unusual or novel
├─ Few similar examples in training
└─ Consider benchmarking instead
```

---

## Step 3: Finding Similar Projects

See what projects your data is similar to:

```python
ml_rec = MLPipelineRecommender()

# Find similar past analyses
similar_projects = ml_rec.get_similar_projects(profile, n=10)

print(f"Found {len(similar_projects)} similar projects:")
for i, proj in enumerate(similar_projects[:5], 1):
    print(f"\n{i}. Project: {proj['project_id']}")
    print(f"   Similarity: {proj['similarity']:.3f}")
    print(f"   Samples: {proj['n_samples']}")
    print(f"   BCV: {proj['bcv']:.3f}")
    print(f"   Organism: {proj['organism']}")
    print(f"   Pipeline used: {proj['pipeline_used']}")
    print(f"   Success: {'✓' if proj['success'] else '✗'}")
```

**Example Output:**
```
Found 10 similar projects:

1. Project: proj_8472
   Similarity: 0.94
   Samples: 12
   BCV: 0.41
   Organism: Homo sapiens
   Pipeline used: Salmon-edgeR
   Success: ✓

2. Project: proj_2391
   Similarity: 0.91
   Samples: 10
   BCV: 0.38
   Organism: Homo sapiens
   Pipeline used: Salmon-edgeR
   Success: ✓
```

**This tells you:**
- Your data is very similar to successful past projects
- They all used Salmon-edgeR
- High chance it will work for you too!

---

## Step 4: Feature Importance

Understand WHY the ML model made its recommendation:

```python
# Get detailed explanation
ml_rec = MLPipelineRecommender()
results = ml_rec.recommend(profile, n=1, explain=True)

top_rec = results[0]

print(f"Recommended: {top_rec['pipeline_name']}")
print(f"\nFeature Importance (why this pipeline?):\n")

for feature, importance in top_rec['top_features']:
    print(f"{feature:20s}: {importance:.3f}")
    
    # Explain what this means
    if feature == 'bcv':
        print(f"  → Your BCV ({profile['bcv']:.3f}) is typical for this pipeline")
    elif feature == 'n_samples':
        print(f"  → Your sample size ({profile['n_samples']}) is ideal for this pipeline")
    elif feature == 'depth':
        print(f"  → Your sequencing depth is well-matched")
```

**Output:**
```
Recommended: Salmon-edgeR

Feature Importance (why this pipeline?):

bcv                 : 0.348
  → Your BCV (0.42) is typical for this pipeline
n_samples           : 0.221
  → Your sample size (12) is ideal for this pipeline
mean_depth          : 0.182
  → Your sequencing depth is well-matched
zero_inflation      : 0.127
library_size_cv     : 0.086
organism_human      : 0.036
```

**Interpretation:**
- BCV is the most important factor (34.8%)
- Your BCV of 0.42 perfectly suits Salmon-edgeR
- Sample size (12) is also optimal
- These two features alone explain 57% of the decision

---

## Step 5: When to Trust ML vs Rules

### Decision Framework

```python
def should_i_trust_ml(ml_results):
    """Decision helper for ML recommendations."""
    
    top = ml_results[0]
    
    # Check 1: Confidence
    if top['confidence'] == 'high':
        print("✅ HIGH CONFIDENCE - Trust the ML recommendation")
        return True
    
    # Check 2: Agreement with rules
    # (Compare with traditional recommendation)
    traditional_rec = PipelineRecommender()
    trad_results = traditional_rec.recommend(profile, n=1)
    
    if top['pipeline_id'] == trad_results[0]['pipeline_id']:
        print("✅ ML AND RULES AGREE - Very reliable")
        return True
    
    # Check 3: Success rate
    if top['historical_success'] > 0.85:
        print("✅ HIGH SUCCESS RATE - ML has proven track record")
        return True
    
    # Check 4: Similar projects
    if top['similar_count'] > 500:
        print("✅ MANY SIMILAR PROJECTS - ML has lots of evidence")
        return True
    
    # If none of above
    if top['confidence'] == 'medium':
        print("⚠️  MEDIUM CONFIDENCE - ML is probably right, but verify")
        return "maybe"
    else:
        print("🔴 LOW CONFIDENCE - Run benchmark to be sure")
        return False

# Use it
ml_results = ml_rec.recommend(profile, n=3, explain=True)
trust_level = should_i_trust_ml(ml_results)

if trust_level == True:
    print("\n→ Go with ML recommendation!")
elif trust_level == "maybe":
    print("\n→ ML recommendation is good, but consider alternatives")
else:
    print("\n→ Run full benchmark to validate")
```

---

## Step 6: Customizing ML Recommendations

### Adjust Confidence Threshold

```python
# Require higher confidence
ml_rec = MLPipelineRecommender(min_confidence=0.8)

# This will only return high-confidence recommendations
results = ml_rec.recommend(profile, n=3)

# If no high-confidence recommendations exist, you'll get fewer results
if len(results) == 0:
    print("No high-confidence recommendations available")
    print("Your data is novel - consider benchmarking")
```

### Weight Different Criteria

```python
# Prioritize historical success over ML score
ml_rec = MLPipelineRecommender()

results = ml_rec.recommend(
    profile,
    n=3,
    weight_success=0.7,  # 70% weight on historical success
    weight_ml_score=0.3  # 30% weight on ML prediction
)
```

### Filter by Resource Constraints

```python
# Only show recommendations that fit your resources
results = ml_rec.recommend(
    profile,
    n=3,
    max_memory_gb=16,  # Maximum 16 GB
    max_runtime_hours=2  # Maximum 2 hours
)

for rec in results:
    print(f"{rec['pipeline_name']}")
    print(f"  Memory: {rec['expected_memory']} GB ✓")
    print(f"  Runtime: {rec['expected_runtime']:.1f} hours ✓")
```

---

## Step 7: Training Custom ML Models

If you have your own benchmark data, train a custom model:

```python
from raptor.ml import train_custom_model
import pandas as pd

# 1. Prepare your training data
# You need: profiles of past analyses + outcomes
training_data = pd.read_csv('my_past_analyses.csv')

# training_data should have:
# - Features: bcv, n_samples, depth, etc.
# - Target: best_pipeline_id (what worked best)
# - Outcome: success (True/False)

X = training_data[['bcv', 'n_samples', 'mean_depth', ...]]
y = training_data['best_pipeline_id']

# 2. Train model
model = train_custom_model(
    X, y,
    model_type='random_forest',  # or 'xgboost', 'neural_network'
    n_estimators=500,
    save_path='my_custom_model.pkl'
)

print(f"Training accuracy: {model.score(X, y):.3f}")

# 3. Use your custom model
ml_rec = MLPipelineRecommender(model_path='my_custom_model.pkl')
results = ml_rec.recommend(profile, n=3)

print("Using YOUR custom model trained on YOUR data!")
```

### What Makes Good Training Data?

```
Minimum Requirements:
├─ 50+ past analyses (100+ recommended)
├─ Diverse data types (different organisms, sample sizes, etc.)
├─ Known outcomes (which pipeline worked best)
└─ Consistent feature extraction

Ideal Training Data:
├─ 500+ analyses
├─ Covers all 8 pipelines
├─ Multiple organisms
├─ Range of sample sizes (3-100)
├─ Range of BCV (0.1-1.0)
└─ Validated results (qPCR, etc.)
```

---

## Step 8: Combining ML with Benchmarking

Best of both worlds: Use ML to narrow down, then benchmark:

```python
from raptor import MLPipelineRecommender, PipelineBenchmark

# 1. Get ML recommendations
ml_rec = MLPipelineRecommender()
ml_results = ml_rec.recommend(profile, n=3)

# 2. Extract top 3 pipeline IDs
top_3_pipelines = [r['pipeline_id'] for r in ml_results]

print(f"ML recommends pipelines: {top_3_pipelines}")
print("Now benchmarking only these 3 instead of all 8...")

# 3. Benchmark only the top 3 (much faster!)
benchmark = PipelineBenchmark(data_dir='fastq/')
results = benchmark.run_pipelines(top_3_pipelines)

print("\nBenchmark validates ML recommendation!")
```

**Time savings:**
- Full benchmark (8 pipelines): 12 hours
- ML-guided benchmark (3 pipelines): 4.5 hours
- **Saved: 7.5 hours (62%)**

---

## Step 9: ML Recommendation Report

Generate a comprehensive ML report:

```python
from raptor import MLPipelineRecommender, ReportGenerator

# Get ML recommendations
ml_rec = MLPipelineRecommender()
ml_results = ml_rec.recommend(profile, n=3, explain=True)

# Generate ML report
report_gen = ReportGenerator()
report_gen.generate_ml_report(
    profile,
    ml_results,
    output='ml_recommendation_report.html',
    include_similar_projects=True,
    include_feature_importance=True,
    include_confidence_analysis=True
)

print("ML report generated!")
print("Open: ml_recommendation_report.html")
```

**Report includes:**
- ML scores and confidence levels
- Feature importance plots
- Similar projects table
- Historical success rates
- Comparison with rule-based recommendations
- Decision confidence visualization

---

## Real-World Examples

### Example 1: High Confidence, Clear Winner

**Your data:**
- 12 samples, BCV=0.38, depth=25M

**ML Output:**
```
#1 Salmon-edgeR
   ML Score: 0.92
   Confidence: HIGH
   Similar projects: 1,847
   Success rate: 89.2%

→ Decision: Trust ML, use Salmon-edgeR
```

### Example 2: Medium Confidence, Two Close Options

**Your data:**
- 6 samples, BCV=0.52, depth=15M

**ML Output:**
```
#1 Salmon-edgeR
   ML Score: 0.79
   Confidence: MEDIUM
   
#2 STAR-DESeq2
   ML Score: 0.77
   Confidence: MEDIUM

→ Decision: Both are good. Try Salmon first (faster).
             If issues, switch to STAR-DESeq2.
```

### Example 3: Low Confidence, Novel Data

**Your data:**
- 3 samples, BCV=0.91, rare organism

**ML Output:**
```
#1 NOISeq
   ML Score: 0.58
   Confidence: LOW
   Similar projects: 23 (very few!)

→ Decision: ML is guessing. Run benchmark on subset,
             or use conservative pipeline (STAR-DESeq2).
```

---

## Troubleshooting

### Problem: "ML model not found"

```bash
# Download default ML model
raptor download-ml-model

# Or specify path
export RAPTOR_ML_MODEL=/path/to/model.pkl
```

### Problem: All recommendations have low confidence

**This means:** Your data is unlike anything in the training set.

**Solutions:**
1. Check if data is actually unusual (rare organism, extreme BCV, etc.)
2. Run benchmarking instead
3. Train custom model on your lab's data
4. Use rule-based recommendations

### Problem: ML disagrees with my expertise

**Remember:** ML learns from data, but you know your experiment!

**Consider:**
- Is your experiment truly novel?
- Does ML lack context (e.g., clinical importance)?
- Are there constraints ML doesn't know (e.g., tool availability)?

**Action:** Use ML as a guide, not a dictator. Your expertise + ML = best results!

---

## Best Practices

### Do's:
✅ Check confidence scores before trusting recommendations  
✅ Review similar projects to understand ML reasoning  
✅ Use ML to narrow down options for benchmarking  
✅ Train custom models with your lab's data  
✅ Combine ML with your domain expertise  

### Don'ts:
❌ Blindly follow low-confidence recommendations  
❌ Ignore rule-based recommendations entirely  
❌ Skip benchmarking for critical analyses  
❌ Train models on tiny datasets (<50 analyses)  
❌ Forget that ML is a tool, not a replacement for thinking  

---

## Summary

You've learned to:
- ✅ Use ML-powered pipeline recommendations
- ✅ Interpret confidence scores and feature importance
- ✅ Find and leverage similar past projects
- ✅ Decide when to trust ML vs rule-based recommendations
- ✅ Customize ML recommendations for your needs
- ✅ Train custom ML models on your data
- ✅ Combine ML with benchmarking for best results

---

## Next Steps

1. **Try Tutorial 5**: [Interactive Dashboard](tutorial_05_dashboard.md)
2. **Try Tutorial 6**: [Ensemble Analysis](tutorial_06_ensemble.md)
3. **Try Tutorial 7**: [Threshold Optimization](tutorial_07_threshold_optimizer.md) 🎯 (NEW!)
4. **Read**: [ML_TRAINING.md](../ML_TRAINING.md) - Deep dive into ML
5. **Read**: [UNDERSTANDING_ML.md](../UNDERSTANDING_ML.md) - ML concepts
6. **Read**: [THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md) - ATO guide

---

**Tutorial by Ayeh Bolouki**  
For RAPTOR v2.1.1

*"Let machines learn, so you can focus on biology!"* 🤖🦖
