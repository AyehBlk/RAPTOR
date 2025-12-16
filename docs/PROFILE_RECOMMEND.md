# RAPTOR v2.1.1 Profile & Recommend Guide

Get intelligent ML-powered pipeline recommendations and **data-driven thresholds** in seconds!

## Quick Start

```bash
# Basic profiling with ML
raptor profile --counts your_data.csv --use-ml

# With metadata and quality check
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

## 🆕 What's New in v2.1.1

### 🎯 Adaptive Threshold Optimizer (ATO)

After getting your pipeline recommendation, use ATO for data-driven thresholds:

```python
from raptor.threshold_optimizer import optimize_thresholds

# After running your recommended pipeline
result = optimize_thresholds(de_results, goal='balanced')
print(f"Optimal |logFC| > {result.logfc_threshold:.2f}")
print(result.methods_text)  # Copy to paper!
```

---

## All Features

⭐ **ML-Powered Recommendations** - Learn from 10,000+ past analyses  
⭐ **Confidence Scores** - Know how reliable recommendations are  
⭐ **Similar Projects** - See what worked for similar data  
⭐ **Quality Assessment** - Automatic data QC  
⭐ **Feature Importance** - Understand WHY each pipeline was recommended  
⭐ **Resource Prediction** - Estimate CPU/memory/time needs  
⭐ **🎯 Threshold Optimization** - Data-driven thresholds (NEW v2.1.1!)

---

## How It Works

### Complete v2.1.1 Workflow

```
1. Analyzes your data → BCV, depth, zeros, library sizes
2. Runs quality checks → Contamination, outliers, batch effects
3. ML model prediction → Based on 10,000+ similar analyses
4. Confidence scoring → How sure is the model?
5. Find similar projects → What worked before?
6. Recommends top 3 → With ML scores, confidence, and reasoning
7. 🎯 Run pipeline → Get DE results
8. 🎯 Optimize thresholds → Data-driven cutoffs (NEW!)
9. 🎯 Get methods text → Publication-ready (NEW!)
```

---

## Understanding Output

### ML-Enhanced Output

```
#1 RECOMMENDED: Pipeline 3 (Salmon-edgeR)
   Rule-based Score: 0.88 | ML Score: 0.92 ⭐
   Confidence: HIGH (89%)
   
   Why this pipeline?
   ✓ Excellent balance for your data
   ✓ Handles medium BCV well (yours: 0.42)
   ✓ Fast runtime (~22 min estimated)
   ✓ Low memory (12GB peak predicted)
   
   ML Insights:
   → Based on 1,247 similar analyses
   → 87.3% historical success rate on similar data
   → Top features: BCV (35%), Sample size (22%), Depth (18%)
   
   Similar Successful Projects:
   • Project #8472: 12 samples, BCV=0.41, Human → Success ✓
   • Project #2391: 10 samples, BCV=0.38, Human → Success ✓
   
   Resource Predictions:
   • Runtime: 22 minutes (±5 min)
   • Peak Memory: 12 GB (±2 GB)

   🎯 Next Step (NEW v2.1.1):
   After running this pipeline, use ATO for optimal thresholds:
   → raptor optimize-thresholds --input results.csv --goal balanced
```

---

## 🎯 Complete Workflow with ATO (NEW)

### Step 1: Profile & Get Recommendation

```bash
raptor profile --counts data.csv --use-ml
# Output: Recommended Pipeline 3 (Salmon-edgeR)
```

### Step 2: Run Recommended Pipeline

```bash
raptor run --pipeline 3 --data fastq/ --output results/
```

### Step 3: Optimize Thresholds (NEW!)

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# Load DE results
df = pd.read_csv('results/edger_results.csv')

# Optimize thresholds
result = optimize_thresholds(
    df,
    logfc_col='logFC',      # edgeR column name
    pvalue_col='PValue',    # edgeR column name
    goal='balanced'
)

# View results
print(f"Optimal |logFC| threshold: {result.logfc_threshold:.3f}")
print(f"Optimal p-value threshold: {result.pvalue_threshold}")
print(f"Significant genes: {result.n_significant}")

# Get publication methods text
print("\n📝 Methods text for your paper:")
print(result.methods_text)

# Save optimized results
result.results_df.to_csv('results/optimized_results.csv')
```

### Step 4: Generate Report

```bash
raptor report --results results/ --output final_report.html
```

---

## Key Metrics Explained

### Data Characteristics

**BCV (Biological Coefficient of Variation)**
- Low (<0.2): Cell lines, controlled conditions
- Medium (0.2-0.6): Typical experiments  
- High (>0.6): Clinical samples, complex biology
- **ML learns which pipelines handle each BCV best**

**Sequencing Depth**
- Low (<10M): May miss genes
- Medium (10-25M): Adequate
- High (>25M): Excellent

**Quality Score**
- 90-100: Excellent data quality
- 80-89: Good quality
- 70-79: Acceptable, minor issues
- <70: Poor quality, address issues first

---

## ML Confidence Levels

### High Confidence (>80%)
```
🟢 HIGH CONFIDENCE

What it means:
✓ Your data is very similar to many past successful analyses
✓ ML model has strong pattern matches
✓ Historical success rate is high (>85%)

What you should do:
→ Trust this recommendation confidently
→ Run pipeline, then use ATO for thresholds (NEW!)
```

### Medium Confidence (60-80%)
```
🟡 MEDIUM CONFIDENCE

What it means:
⚠️  Your data is somewhat similar to past analyses
⚠️  Historical success rate is good (70-85%)

What you should do:
→ Consider the top 2-3 recommendations
→ Use ATO with 'balanced' goal for safe thresholds (NEW!)
```

### Low Confidence (<60%)
```
🔴 LOW CONFIDENCE

What it means:
❌ Your data is unusual or novel
❌ ML model hasn't seen many similar examples

What you should do:
→ Run benchmark on multiple pipelines
→ Use ATO with 'discovery' goal to cast wider net (NEW!)
```

---

## 🎯 ATO Goals for Different Scenarios (NEW)

After running your pipeline, choose the right ATO goal:

### Discovery Goal
```python
result = optimize_thresholds(de_results, goal='discovery')
```
**Use when:**
- Exploratory analysis
- Want more candidate genes
- Will validate later
- Pilot study

**Effect:** More permissive thresholds

### Balanced Goal (Recommended)
```python
result = optimize_thresholds(de_results, goal='balanced')
```
**Use when:**
- Standard publication
- Most typical analyses
- Good balance needed

**Effect:** Standard FDR control

### Validation Goal
```python
result = optimize_thresholds(de_results, goal='validation')
```
**Use when:**
- Clinical applications
- Confirmation study
- Need high confidence

**Effect:** Stringent thresholds

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
```

### Resource Constraints

```bash
# Limited resources
raptor profile \
  --counts data.csv \
  --use-ml \
  --max-memory 16G \
  --max-runtime 2h
```

### Exclude Pipelines

```bash
# Don't recommend certain pipelines
raptor profile \
  --counts data.csv \
  --use-ml \
  --exclude-pipelines 7,8
```

---

## Quality Assessment Integration

```bash
# Run comprehensive quality check
raptor profile \
  --counts data.csv \
  --metadata metadata.csv \
  --quality-check
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
⚠️  Minor batch effect detected
   → Recommendation: Include batch in model

✅ Data quality is good for analysis
   → Use ATO with 'balanced' goal after pipeline runs
```

---

## Python API

### Complete Workflow

```python
from raptor import RNAseqDataProfiler, MLPipelineRecommender
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

# 1. Load data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# 2. Profile with ML
profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile(quality_check=True)

# 3. Get ML recommendations
ml_recommender = MLPipelineRecommender()
recommendations = ml_recommender.recommend(profile, n=3, explain=True)

# 4. View top recommendation
top = recommendations[0]
print(f"Recommended: {top['pipeline_name']}")
print(f"Confidence: {top['confidence']}")

# 5. [Run the recommended pipeline]
# raptor run --pipeline 3 ...

# 6. Load DE results and optimize thresholds (NEW!)
de_results = pd.read_csv('results/edger_results.csv')
threshold_result = optimize_thresholds(
    de_results,
    logfc_col='logFC',
    pvalue_col='PValue',
    goal='balanced'
)

print(f"\n🎯 Optimized Thresholds:")
print(f"   |logFC| > {threshold_result.logfc_threshold:.3f}")
print(f"   Significant: {threshold_result.n_significant} genes")

# 7. Get methods text for publication
print(f"\n📝 Methods:\n{threshold_result.methods_text}")

# 8. Save optimized results
threshold_result.results_df.to_csv('final_results.csv')
```

### Get Similar Projects

```python
# Find similar past analyses
similar = ml_recommender.get_similar_projects(profile, n=10)

for project in similar:
    print(f"Project {project['id']}")
    print(f"  Similarity: {project['similarity']:.3f}")
    print(f"  Pipeline used: {project['pipeline']}")
    print(f"  Success: {project['success']}")
```

### Explain Recommendation

```python
# Get feature importance
explanation = ml_recommender.explain_recommendation(profile, pipeline_id=3)

print("Feature Importance:")
for feature, importance in explanation['features']:
    print(f"  {feature}: {importance:.3f}")
```

---

## Examples

### Standard DE Study (12 samples)

**Step 1: Profile**
```bash
raptor profile --counts data.csv --use-ml
```

**Output:**
```
✅ Pipeline 3 (Salmon-edgeR)
   ML Score: 0.92 | Confidence: HIGH (89%)
   
   Why: Perfect for your data
   - Medium BCV (0.42) → Optimal range
   - 12 samples → Ideal for this pipeline
   
   🎯 Next: Run pipeline, then optimize thresholds with ATO
```

**Step 2: Run & Optimize**
```python
# After running pipeline
result = optimize_thresholds(de_results, goal='balanced')
print(f"Use: |logFC| > {result.logfc_threshold:.2f}")
# Output: Use: |logFC| > 0.73
```

---

### Large Cohort (100 samples)

**Profile:**
```bash
raptor profile --counts large_data.csv --use-ml
```

**Output:**
```
✅ Pipeline 4 (Kallisto-DESeq2)
   ML Score: 0.90 | Confidence: HIGH (88%)
   
   Why: Optimized for large datasets
   - 100 samples → Kallisto scales well
   - Fast runtime critical for this size
   
   🎯 Next: Use ATO with 'discovery' goal for exploratory analysis
```

**Optimize:**
```python
result = optimize_thresholds(de_results, goal='discovery')
# More permissive for large exploratory study
```

---

### Small Pilot (4 samples)

**Profile:**
```bash
raptor profile --counts pilot.csv --use-ml --weight-accuracy 0.8
```

**Output:**
```
✅ Pipeline 1 (STAR-RSEM-DESeq2)
   ML Score: 0.88 | Confidence: HIGH (85%)
   
   Why: Best for small samples
   - Only 4 samples → Need robust method
   - DESeq2 handles small n well
   
   🎯 Next: Use ATO with 'validation' goal for stringent thresholds
```

**Optimize:**
```python
result = optimize_thresholds(de_results, goal='validation')
# Stringent thresholds for small, important study
```

---

## Troubleshooting

### Problem: All recommendations have low confidence

**Solutions:**
1. Check if data characteristics are extreme
2. Run quality check to identify issues
3. Consider benchmarking instead of single pipeline
4. **Use ATO with 'discovery' goal to be more permissive**

### Problem: ML disagrees with rule-based

**What to do:**
- If ML confidence is HIGH → Trust ML
- If ML confidence is MEDIUM → Consider both
- If ML confidence is LOW → Trust rules
- **Always use ATO afterward for optimal thresholds**

### Problem: Quality score is low

```bash
# Get detailed quality report
raptor profile --counts data.csv --quality-check --verbose

# Address issues, then re-profile
raptor profile --counts cleaned_data.csv --use-ml

# After pipeline, ATO will still optimize for your data quality
```

---

## Best Practices

### Do's:
✅ Always use `--use-ml` for recommendations  
✅ Check ML confidence scores  
✅ Review quality assessment  
✅ **Use ATO after running pipeline (NEW!)**  
✅ **Include ATO methods text in publications (NEW!)**  
✅ Consider top 2-3 recommendations  

### Don'ts:
❌ Blindly follow low-confidence recommendations  
❌ Ignore quality warnings  
❌ **Use arbitrary thresholds - use ATO instead! (NEW!)**  
❌ Skip quality checks for critical analyses  
❌ Forget that ML is a tool, not a dictator  

---

## See Also

- [THRESHOLD_OPTIMIZER.md](THRESHOLD_OPTIMIZER.md) - 🎯 Complete ATO guide (NEW!)
- [UNDERSTANDING_ML.md](UNDERSTANDING_ML.md) - ML concepts explained
- [QUALITY_ASSESSMENT.md](QUALITY_ASSESSMENT.md) - Data QC guide
- [PIPELINES.md](PIPELINES.md) - Detailed pipeline information
- [DASHBOARD.md](DASHBOARD.md) - Interactive visualization + ATO page

---

**RAPTOR v2.1.1 Profile & Recommend**  
**Author**: Ayeh Bolouki  
**License**: MIT

*"Profile smartly, recommend confidently, threshold optimally!"* 🎯🦖
