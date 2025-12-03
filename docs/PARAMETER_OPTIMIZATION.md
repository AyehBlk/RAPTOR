#  RAPTOR v2.1.0 Parameter Optimization Guide

**Automated Parameter Tuning for Optimal Results**

Learn how to find the best parameters for your RNA-seq analysis using automated optimization techniques.

---

##  Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [What Can Be Optimized](#what-can-be-optimized)
4. [Optimization Methods](#optimization-methods)
5. [Grid Search](#grid-search)
6. [Bayesian Optimization](#bayesian-optimization)
7. [Sensitivity Analysis](#sensitivity-analysis)
8. [Multi-Objective Optimization](#multi-objective-optimization)
9. [Interpreting Results](#interpreting-results)
10. [Best Practices](#best-practices)

---

##  Overview

### What is Parameter Optimization?

Parameter optimization finds the best settings for your analysis to maximize:
- **Accuracy** - Find more true positives
- **Sensitivity** - Detect subtle differences
- **Specificity** - Reduce false positives
- **Reproducibility** - Consistent results

### Why Optimize Parameters?

**Default parameters:**
```
FDR threshold: 0.05
Log2FC threshold: 1.0
Min counts: 10
→ Results: 1,247 DE genes
```

**Optimized parameters:**
```
FDR threshold: 0.03
Log2FC threshold: 0.75
Min counts: 5
→ Results: 1,689 DE genes (with better accuracy!)
```

**Benefit:** 35% more genes detected, 12% better validation rate

---

##  Quick Start

### Automatic Optimization

```bash
# Optimize all parameters automatically
raptor optimize \
  --counts data/counts.csv \
  --metadata data/metadata.csv \
  --ground-truth validated_genes.csv \
  --output optimized_params/
```

**Output:**
```
Parameter Optimization Complete!
═══════════════════════════════

Optimized Parameters:
├─ FDR threshold: 0.032 (was 0.05)
├─ Log2FC threshold: 0.78 (was 1.0)
├─ Min count filter: 7 (was 10)
├─ Normalization: TMM (was default)
└─ Prior.df: 3.5 (was 4.0)

Performance Improvement:
├─ Sensitivity: 0.82 → 0.89 (+8.5%)
├─ Precision: 0.85 → 0.88 (+3.5%)
├─ F1-Score: 0.83 → 0.88 (+6.0%)

Results saved to: optimized_params/best_parameters.yaml

To use optimized parameters:
raptor analyze --config optimized_params/best_parameters.yaml
```

### Quick Optimization (Subset of Parameters)

```bash
# Optimize just thresholds (fast)
raptor optimize-quick \
  --counts data/counts.csv \
  --parameters fdr,log2fc
```

---

##  What Can Be Optimized

### Statistical Parameters

```yaml
statistical:
  # Significance thresholds
  fdr_threshold: [0.01, 0.05, 0.10]
  pvalue_threshold: [0.001, 0.01, 0.05]
  
  # Effect size thresholds
  log2fc_threshold: [0.5, 0.75, 1.0, 1.5]
  min_fold_change: [1.5, 2.0, 2.5]
  
  # Multiple testing correction
  adjustment_method: [BH, bonferroni, qvalue]
```

### Filtering Parameters

```yaml
filtering:
  # Count filtering
  min_counts: [1, 5, 10, 20]
  min_total_counts: [10, 50, 100]
  min_samples_expressed: [2, 3, 4]
  
  # Gene filtering
  min_cpm: [0.5, 1.0, 2.0]
  keep_top_n_genes: [10000, 15000, 20000]
```

### Normalization Parameters

```yaml
normalization:
  # Method selection
  method: [TMM, RLE, upperquartile]
  
  # Method-specific parameters
  tmm_log_ratio_trim: [0.2, 0.3, 0.4]
  tmm_sum_trim: [0.05, 0.1, 0.15]
```

### Pipeline-Specific Parameters

**edgeR:**
```yaml
edger:
  prior_df: [2.0, 3.0, 4.0, 5.0]
  robust: [true, false]
  trend_method: [none, locfit, loess]
```

**DESeq2:**
```yaml
deseq2:
  fit_type: [parametric, local, mean]
  beta_prior: [true, false]
  alpha: [0.05, 0.1]
  cook_cutoff: [true, false]
```

**limma-voom:**
```yaml
limma:
  robust: [true, false]
  trend: [true, false]
  span: [0.3, 0.5, 0.7]
```

---

##  Optimization Methods

### 1. Grid Search (Thorough)

**How it works:** Try every combination of parameters

```
Parameter Space:
FDR: [0.01, 0.05, 0.10]
Log2FC: [0.5, 1.0, 1.5]

Grid Search tries:
├─ FDR=0.01, Log2FC=0.5 → F1=0.84
├─ FDR=0.01, Log2FC=1.0 → F1=0.82
├─ FDR=0.01, Log2FC=1.5 → F1=0.78
├─ FDR=0.05, Log2FC=0.5 → F1=0.87 ✓ Best!
├─ FDR=0.05, Log2FC=1.0 → F1=0.85
├─ FDR=0.05, Log2FC=1.5 → F1=0.81
├─ FDR=0.10, Log2FC=0.5 → F1=0.86
├─ FDR=0.10, Log2FC=1.0 → F1=0.84
└─ FDR=0.10, Log2FC=1.5 → F1=0.80

Best: FDR=0.05, Log2FC=0.5, F1=0.87
```

**Pros:** Guaranteed to find best combination  
**Cons:** Slow for many parameters

---

### 2. Random Search (Faster)

**How it works:** Try random combinations

```
Instead of all 9 combinations,
try 5 random ones:

├─ FDR=0.08, Log2FC=0.7 → F1=0.85
├─ FDR=0.03, Log2FC=1.2 → F1=0.81
├─ FDR=0.06, Log2FC=0.6 → F1=0.86
├─ FDR=0.04, Log2FC=0.9 → F1=0.84
└─ FDR=0.05, Log2FC=0.5 → F1=0.87 ✓

Best: FDR=0.05, Log2FC=0.5, F1=0.87
```

**Pros:** Much faster than grid search  
**Cons:** Might miss optimal combination

---

### 3. Bayesian Optimization (Smart)

**How it works:** Learn from previous tries to pick next best

```
Iteration 1: Try FDR=0.05, Log2FC=1.0 → F1=0.85
Model predicts: "Try lower Log2FC"

Iteration 2: Try FDR=0.05, Log2FC=0.5 → F1=0.87 ✓
Model predicts: "Try nearby values"

Iteration 3: Try FDR=0.04, Log2FC=0.6 → F1=0.86
Model predicts: "Iteration 2 is best"

Stop: Found optimum in 3 tries!
(Grid search would take 9 tries)
```

**Pros:** Very efficient, finds optimum quickly  
**Cons:** More complex to set up

---

##  Grid Search

### Basic Grid Search

```bash
raptor optimize grid \
  --counts data/counts.csv \
  --ground-truth validated.csv \
  --parameters grid_config.yaml \
  --output grid_results/
```

**grid_config.yaml:**
```yaml
optimization:
  method: grid_search
  
  parameters:
    fdr_threshold:
      values: [0.01, 0.05, 0.10]
    
    log2fc_threshold:
      values: [0.5, 0.75, 1.0, 1.5]
    
    min_counts:
      values: [5, 10, 20]
  
  metric: f1_score  # What to optimize
  cv_folds: 5       # Cross-validation
```

### Results

```
Grid Search Results
═══════════════════

Tested 36 combinations (3 × 4 × 3)

Best Parameters:
├─ FDR threshold: 0.05
├─ Log2FC threshold: 0.75
└─ Min counts: 10

Performance:
├─ F1-Score: 0.88 (±0.03)
├─ Sensitivity: 0.89
├─ Precision: 0.87

Top 5 Combinations:
┌────────┬──────────┬────────────┬──────────┐
│ FDR    │ Log2FC   │ Min Counts │ F1-Score │
├────────┼──────────┼────────────┼──────────┤
│ 0.05   │ 0.75     │ 10         │ 0.88     │
│ 0.05   │ 0.50     │ 10         │ 0.87     │
│ 0.05   │ 0.75     │ 5          │ 0.87     │
│ 0.10   │ 0.75     │ 10         │ 0.86     │
│ 0.05   │ 1.00     │ 10         │ 0.85     │
└────────┴──────────┴────────────┴──────────┘

Total time: 45 minutes
```

### Visualization

```bash
raptor optimize plot \
  --results grid_results/ \
  --output heatmap.png
```

**Output: Parameter Heatmap**
```
F1-Score by Parameter Combination

Log2FC
  ↑
1.5 │ 0.78  0.80  0.81
1.0 │ 0.82  0.85  0.84
0.75│ 0.86  0.88  0.87  ← Best!
0.5 │ 0.85  0.87  0.86
    └────────────────→ FDR
    0.01  0.05  0.10
```

---

##  Bayesian Optimization

### Setup

```bash
raptor optimize bayesian \
  --counts data/counts.csv \
  --ground-truth validated.csv \
  --parameters bayes_config.yaml \
  --iterations 30 \
  --output bayes_results/
```

**bayes_config.yaml:**
```yaml
optimization:
  method: bayesian
  
  parameters:
    fdr_threshold:
      type: continuous
      range: [0.001, 0.20]
      log_scale: true
    
    log2fc_threshold:
      type: continuous
      range: [0.1, 2.0]
    
    min_counts:
      type: integer
      range: [1, 50]
    
    prior_df:
      type: continuous
      range: [1.0, 10.0]
  
  acquisition_function: expected_improvement
  n_initial_points: 5
  random_state: 42
```

### Optimization Process

```
Bayesian Optimization Progress
════════════════════════════════

Iteration 1: [Random initialization]
  FDR=0.05, Log2FC=1.0, MinCount=10, Prior=4.0
  F1-Score: 0.85

Iteration 2: [Random initialization]
  FDR=0.10, Log2FC=0.5, MinCount=20, Prior=3.0
  F1-Score: 0.83

...

Iteration 6: [Model-guided]
  Model predicts: Try FDR=0.032
  FDR=0.032, Log2FC=0.78, MinCount=7, Prior=3.5
  F1-Score: 0.89 ✓ New best!

Iteration 7: [Model-guided]
  Model predicts: Try nearby
  FDR=0.028, Log2FC=0.82, MinCount=6, Prior=3.8
  F1-Score: 0.88 (worse)

...

Iteration 30: [Convergence]
  No improvement in last 10 iterations
  Best found: F1=0.89
  Stopping early.

Best Parameters:
├─ FDR threshold: 0.032
├─ Log2FC threshold: 0.78
├─ Min counts: 7
└─ Prior.df: 3.5

Total evaluations: 30 (vs 1000+ for full grid)
Time: 15 minutes (vs 5+ hours for grid)
```

---

##  Sensitivity Analysis

### What is Sensitivity Analysis?

Test how sensitive results are to parameter changes.

```bash
raptor optimize sensitivity \
  --counts data/counts.csv \
  --output sensitivity_results/
```

### Results

```
Parameter Sensitivity Analysis
═══════════════════════════════

1. FDR Threshold
   Range tested: 0.01 to 0.10
   
   Impact on Number of DE Genes:
   0.01: 543 genes   ▂
   0.02: 892 genes   ▄
   0.03: 1245 genes  ▆
   0.04: 1534 genes  ███
   0.05: 1789 genes  ████  ← Standard
   0.06: 2012 genes  █████
   0.08: 2456 genes  ██████
   0.10: 2789 genes  ███████
   
   Sensitivity: HIGH
   ├─ Small changes = big impact
   └─ Be careful with this parameter!

2. Log2FC Threshold
   Range tested: 0.5 to 2.0
   
   Impact on Number of DE Genes:
   0.5:  2134 genes  ██████
   0.75: 1789 genes  ████  
   1.0:  1456 genes  ███   ← Standard
   1.5:  892 genes   ▂
   2.0:  423 genes   ▁
   
   Sensitivity: HIGH

3. Min Count Filter
   Range tested: 1 to 50
   
   Impact on Number of DE Genes:
   1:    1812 genes  ████
   5:    1798 genes  ████
   10:   1789 genes  ████  ← Standard
   20:   1756 genes  ███
   50:   1634 genes  ███
   
   Sensitivity: LOW
   ├─ Robust to changes
   └─ Safe to adjust

4. Prior.df (edgeR)
   Range tested: 1.0 to 10.0
   
   Impact on results:
   1.0:  More liberal
   4.0:  Standard
   10.0: More conservative
   
   Sensitivity: MEDIUM

Recommendations:
├─ Use conservative FDR (0.01-0.03) for discovery
├─ Use standard FDR (0.05) for publication
├─ Log2FC threshold matters more than min counts
└─ Prior.df: keep near default (3-5)
```

---

##  Multi-Objective Optimization

### Optimize Multiple Goals

```bash
raptor optimize multi \
  --counts data/counts.csv \
  --objectives sensitivity,precision,speed \
  --output multi_results/
```

**Trade-offs:**
```
Multi-Objective Optimization
═══════════════════════════

Objective 1: Maximize Sensitivity (find all true positives)
Objective 2: Maximize Precision (minimize false positives)
Objective 3: Minimize Runtime (fast analysis)

Pareto Front: (no single "best", depends on priorities)

Solution A: High Sensitivity
├─ Sensitivity: 0.95 ✓✓✓
├─ Precision: 0.72 ✓
├─ Runtime: 45 min ✓
└─ Use when: Don't want to miss anything

Solution B: Balanced
├─ Sensitivity: 0.89 ✓✓
├─ Precision: 0.88 ✓✓
├─ Runtime: 25 min ✓✓
└─ Use when: Good all-around (RECOMMENDED)

Solution C: High Precision
├─ Sensitivity: 0.78 ✓
├─ Precision: 0.95 ✓✓✓
├─ Runtime: 15 min ✓✓✓
└─ Use when: Need to be very sure (follow-up expensive)

Choose based on your goals!
```

---

##  Interpreting Results

### Understanding Metrics

**Sensitivity (Recall):**
```
Sensitivity = True Positives / (True Positives + False Negatives)

Example:
200 real DE genes exist
180 found by our method
Sensitivity = 180/200 = 0.90 (90%)

High sensitivity = Find most true DE genes
Low sensitivity = Miss many true DE genes
```

**Precision:**
```
Precision = True Positives / (True Positives + False Positives)

Example:
200 genes called DE
180 are truly DE
20 are false positives
Precision = 180/200 = 0.90 (90%)

High precision = Few false positives
Low precision = Many false positives
```

**F1-Score:**
```
F1 = 2 × (Precision × Sensitivity) / (Precision + Sensitivity)

Balanced measure of both
F1 = 1.0 is perfect
F1 = 0.5 is poor
```

### Trade-offs

```
Parameter      Effect on         Effect on
               Sensitivity       Precision
─────────────────────────────────────────────
↑ FDR          ↑ (more genes)    ↓ (more FP)
↓ FDR          ↓ (fewer genes)   ↑ (fewer FP)

↑ Log2FC       ↓ (fewer genes)   ↑ (stronger)
↓ Log2FC       ↑ (more genes)    ↓ (weaker)

↑ Min counts   ↓ (fewer genes)   ↑ (reliable)
↓ Min counts   ↑ (more genes)    ↓ (noisy)
```

**Example:**
```
Conservative (High Precision):
├─ FDR = 0.01
├─ Log2FC = 1.5
└─ Result: 500 genes, 95% precision, 75% sensitivity
   → Use when: Validation is expensive

Liberal (High Sensitivity):
├─ FDR = 0.10
├─ Log2FC = 0.5
└─ Result: 2000 genes, 75% precision, 95% sensitivity
   → Use when: Don't want to miss anything

Balanced:
├─ FDR = 0.05
├─ Log2FC = 0.75
└─ Result: 1200 genes, 85% precision, 87% sensitivity
   → Use for: Most situations
```

---

##  Best Practices

### 1. Always Use Validation Data

```bash
# WRONG: Optimize on same data you'll analyze
raptor optimize --counts data.csv --ground-truth data.csv

# RIGHT: Use separate validation set
raptor optimize \
  --counts training_data.csv \
  --ground-truth validation_data.csv
```

### 2. Use Cross-Validation

```bash
raptor optimize \
  --cv-folds 5 \  # 5-fold cross-validation
  --stratified    # Maintain class balance
```

### 3. Start Simple

```bash
# Step 1: Optimize just thresholds (fast)
raptor optimize-quick --parameters fdr,log2fc

# Step 2: If needed, optimize more
raptor optimize-full --parameters all
```

### 4. Check Sensitivity

```bash
# Always check how sensitive results are
raptor optimize sensitivity
```

### 5. Document Everything

```yaml
# Save optimization settings
optimization_run:
  date: 2025-11-19
  method: bayesian
  iterations: 30
  best_f1: 0.89
  best_params:
    fdr: 0.032
    log2fc: 0.78
  validation_data: validated_genes.csv
```

### 6. Consider Your Goals

**Exploratory research:**
- Higher sensitivity
- Accept more FP
- Use liberal thresholds

**Validation studies:**
- Higher precision
- Minimize FP
- Use conservative thresholds

**Publications:**
- Balance both
- Use standard thresholds (FDR=0.05)
- Show sensitivity analysis

---

##  Example Workflows

### Workflow 1: Quick Threshold Optimization

```bash
# 5-minute optimization
raptor optimize-quick \
  --counts data/counts.csv \
  --parameters fdr,log2fc \
  --n-trials 20

# Use optimized parameters
raptor analyze \
  --counts data/counts.csv \
  --fdr 0.032 \
  --log2fc 0.78
```

### Workflow 2: Comprehensive Optimization

```bash
# 1. Full parameter optimization (1-2 hours)
raptor optimize bayesian \
  --counts data/counts.csv \
  --ground-truth validated.csv \
  --parameters all \
  --iterations 50 \
  --output optimization/

# 2. Sensitivity analysis
raptor optimize sensitivity \
  --results optimization/ \
  --output sensitivity/

# 3. Use best parameters
raptor analyze --config optimization/best_parameters.yaml
```

### Workflow 3: Multi-Objective for Different Use Cases

```bash
# Optimize for different scenarios
raptor optimize multi \
  --counts data/counts.csv \
  --objectives sensitivity,precision \
  --output multi_objective/

# Get parameters for different goals:
# - discovery_params.yaml (high sensitivity)
# - validation_params.yaml (high precision)
# - balanced_params.yaml (balanced)

# Use appropriate one
raptor analyze --config multi_objective/balanced_params.yaml
```

---

##  Summary

Parameter optimization provides:
- ✅ **Better accuracy** - Find more true positives
- ✅ **Fewer false positives** - More reliable results
- ✅ **Automated tuning** - No guessing
- ✅ **Validated parameters** - Proven on your data
- ✅ **Understanding** - Know parameter effects
- ✅ **Customization** - Tailored to your needs

**Optimize once, benefit forever!** 🎯🚀

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.0  
**License:** MIT

---

