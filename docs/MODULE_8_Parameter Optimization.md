# RAPTOR Module 8: Parameter Optimization - User Guide

**A Comprehensive Guide to Optimizing Your Differential Expression Analysis**

**Version**: 2.2.0  
**Last Updated**: March 9, 2026  
**Difficulty**: Intermediate  
**Estimated Time**: 15-30 minutes per optimization

---

## 📚 Table of Contents

1. [Introduction](#introduction)
2. [Why Optimize Parameters?](#why-optimize-parameters)
3. [The Four Optimization Methods](#the-four-optimization-methods)
4. [Getting Started](#getting-started)
5. [Method 1: Ground Truth Optimization](#method-1-ground-truth-optimization)
6. [Method 2: FDR Control Optimization](#method-2-fdr-control-optimization)
7. [Method 3: Stability Optimization](#method-3-stability-optimization)
8. [Method 4: Reproducibility Optimization](#method-4-reproducibility-optimization)
9. [Understanding Your Results](#understanding-your-results)
10. [Choosing the Right Method](#choosing-the-right-method)
11. [Advanced Options](#advanced-options)
12. [Integration with Your Workflow](#integration-with-your-workflow)
13. [Best Practices](#best-practices)
14. [Troubleshooting](#troubleshooting)
15. [Frequently Asked Questions](#frequently-asked-questions)
16. [Real-World Examples](#real-world-examples)
17. [Scientific Background](#scientific-background)

---

## Introduction

### What is Module 8?

Module 8 is RAPTOR's **Parameter Optimization** system that helps you find the best FDR (False Discovery Rate) threshold and log-fold-change cutoff for your differential expression analysis.

**The Problem**: Most researchers use default thresholds (FDR < 0.05, |log2FC| > 0) without knowing if these are optimal for their specific dataset.

**The Solution**: Module 8 uses scientifically validated methods to find the thresholds that work best for YOUR data.

### What Gets Optimized?

Module 8 optimizes **two critical parameters**:

1. **Alpha (FDR threshold)**: The significance cutoff (e.g., 0.05, 0.01, 0.001)
   - Controls how many genes are called "differentially expressed"
   - Too lenient → many false positives
   - Too strict → miss real changes

2. **LFC threshold**: The minimum fold-change (e.g., 0, 0.5, 1.0)
   - Filters genes by biological effect size
   - 0 = no filter
   - 0.5 = 1.4-fold change minimum
   - 1.0 = 2-fold change minimum

### Where Does Module 8 Fit?

```
Your RNA-seq Analysis Pipeline:

Step 1: Quality Control (Module 2)
         ↓
Step 2: Data Profiling (Module 3)
         ↓
Step 3: Pipeline Selection (Module 4)
         ↓
Step 4: Run DE Analysis (External - DESeq2/edgeR/limma)
         ↓
Step 5: Import Results (Module 6)
         ↓
Step 6: OPTIMIZE PARAMETERS (Module 8) ⭐ YOU ARE HERE
         ↓
Step 7: Apply Optimal Thresholds
         ↓
Step 8: Ensemble Analysis (Module 9)
         ↓
Step 9: Biomarker Discovery (Module 10)
```

---

## Why Optimize Parameters?

### The Problem with Default Thresholds

**Most RNA-seq papers use**: FDR < 0.05, |log2FC| > 0

**But these defaults**:
- ✗ May not be appropriate for your sample size
- ✗ Don't account for your biological question
- ✗ Can lead to too many or too few genes
- ✗ Aren't validated for your specific data

### Real-World Impact

**Example 1: Cancer Study**
- Default thresholds: 2,500 DEGs
- Optimized thresholds: 850 DEGs
- **Result**: More focused gene set, better pathway enrichment

**Example 2: Drug Response**
- Default thresholds: 120 DEGs
- Optimized thresholds: 450 DEGs
- **Result**: Recovered known drug targets that were missed

**Example 3: Developmental Study**
- Default thresholds: FDR < 0.05
- Optimized thresholds: FDR < 0.0089
- **Result**: Better replication in validation cohort (85% → 94%)

### Benefits of Optimization

1. ✅ **More reliable results** - Validated thresholds for your data
2. ✅ **Better reproducibility** - Optimized for stability
3. ✅ **Justify your choices** - Scientific rationale for thresholds
4. ✅ **Improved downstream analysis** - Better gene sets for enrichment
5. ✅ **Publication quality** - Reviewers appreciate rigorous methods

---

## The Four Optimization Methods

Module 8 offers **4 different approaches** to find optimal parameters. Choose based on what data you have available:

### Quick Comparison Table

| Method | What You Need | When to Use | Computational Time |
|--------|--------------|-------------|-------------------|
| **Ground Truth** | Validated gene list (10-50+ genes) | You have known DE genes from literature | ~5 minutes |
| **FDR Control** | Nothing extra! | General purpose, no validation data | ~3 minutes |
| **Stability** | Count matrix + metadata | Want robust, reproducible results | ~15-30 minutes |
| **Reproducibility** | Two independent cohorts | Have validation dataset | ~5 minutes |

### Decision Tree

```
START: Which method should I use?

Do you have a list of validated DE genes from experiments or literature?
├─ YES → Use Ground Truth (gold standard)
│         Need: CSV file with gene IDs
│
└─ NO → Continue...

Do you have two independent cohorts (discovery + validation)?
├─ YES → Use Reproducibility
│         Need: DE results from both cohorts
│
└─ NO → Continue...

Do you have the original count matrix and metadata?
├─ YES → Use Stability
│         Need: counts.csv + metadata.csv
│
└─ NO → Use FDR Control
          Need: Only the DE results (you already have this!)
```

---

## Getting Started

### Prerequisites

Before running Module 8, you need:

1. ✅ **DE results imported** (from Module 6)
   ```bash
   raptor import-de -i deseq2_results.csv -m deseq2 -o results/de_imported
   ```

2. ✅ **DE result pickle file**: `de_result.pkl`
   - Created automatically by Module 6
   - Location: `results/de_imported/de_result.pkl`

### Quick Start

**Simplest command** (uses FDR control - no extra data needed):
```bash
raptor optimize -d results/de_imported/de_result.pkl -m fdr-control
```

That's it! Module 8 will:
- Analyze your DE results
- Search for optimal parameters
- Generate optimized thresholds
- Save results to `results/optimization/`

### Expected Output

After running, you'll see:
```
🦖 RAPTOR v2.2.0 - Module 8: Parameter Optimization

📂 Loading DE result...
   ✓ Loaded 15000 genes

🔬 Running fdr-control optimization...
   ✓ Optimization complete

📊 Optimal Parameters:
============================================================
Alpha (FDR) threshold: 0.0123
LFC threshold: 0.50
Best Score (fdr_distance): 0.0011
DEG genes at optimal: 1250

✅ Results saved:
   results/optimization/optimized_params.yaml
   results/optimization/deg_genes.csv (1250 genes)
   results/optimization/optimization_result.pkl
   results/optimization/convergence_history.json
```

---

## Method 1: Ground Truth Optimization

**Use when**: You have validated genes from literature, experiments, or databases

### What is Ground Truth Optimization?

This method finds parameters that maximize overlap with genes you **know** are differentially expressed. It's the **gold standard** when you have validated data.

### When to Use This

✅ **Perfect for**:
- You have confirmed DE genes from qRT-PCR validation
- Literature reports specific genes for your condition
- You have genes from public databases (e.g., CancerGeneNet, DisGeNET)
- Previous studies on same biological system

❌ **Don't use when**:
- You don't have any validated genes
- Your gene list is biased (only looked at certain genes)
- Validated genes are from a different biological context

### Required Data

**Ground truth file format** (`validated_genes.csv`):

```csv
gene_id
TP53
BRCA1
EGFR
MYC
KRAS
...
```

**Requirements**:
- **Minimum**: 10 genes (absolute minimum)
- **Recommended**: 20-30 genes (better optimization)
- **Ideal**: 50+ genes (robust results)

### Step-by-Step Tutorial

#### Step 1: Prepare Your Validated Gene List

Create a CSV file with one column: `gene_id`

**Example sources**:
- Literature review of your condition
- qRT-PCR validated genes
- Public databases (MSigDB, KEGG, Reactome)
- Previous RNA-seq studies

```bash
# Example: validated_genes.csv
echo "gene_id" > validated_genes.csv
echo "TP53" >> validated_genes.csv
echo "BRCA1" >> validated_genes.csv
echo "EGFR" >> validated_genes.csv
# ... add your genes
```

#### Step 2: Run Ground Truth Optimization

```bash
raptor optimize \
    -d results/de_imported/de_result.pkl \
    -m ground-truth \
    -g validated_genes.csv \
    -o results/opt_ground_truth
```

**What happens**:
1. Loads your DE results (15,000 genes)
2. Loads your validated genes (e.g., 25 genes)
3. Tries different FDR/LFC thresholds
4. Finds combination that best captures validated genes
5. Balances precision (avoiding false positives) and recall (finding true positives)

#### Step 3: Review Results

Check `results/opt_ground_truth/optimized_params.yaml`:

```yaml
best_parameters:
  alpha: 0.0234
  lfc_threshold: 0.75
best_score: 0.8456
optimization_method: ground_truth
metric: f1_score
n_deg_genes: 1250
performance_metrics:
  precision: 0.85
  recall: 0.84
  f1_score: 0.85
  true_positives: 21
  false_positives: 4
  false_negatives: 4
```

**Interpretation**:
- **Best alpha**: 0.0234 (stricter than default 0.05)
- **Best LFC**: 0.75 (require ~1.7-fold change)
- **F1-score**: 0.85 (excellent balance)
- **Found**: 21 of your 25 validated genes (84% recall)
- **Precision**: 85% of called genes are likely true

### Understanding the Metric: F1-Score

Ground truth optimization uses **F1-score** to balance precision and recall:

```
Precision = True Positives / (True Positives + False Positives)
           = How many of the genes we called are actually DE?

Recall = True Positives / (True Positives + False Negatives)
       = How many of the known DE genes did we find?

F1-score = 2 × (Precision × Recall) / (Precision + Recall)
         = Balance between precision and recall
```

**Good F1-scores**:
- 0.7-0.8: Good
- 0.8-0.9: Excellent
- 0.9+: Outstanding

### Example Output

```bash
$ raptor optimize -d de_result.pkl -m ground-truth -g validated_genes.csv

🦖 RAPTOR v2.2.0 - Module 8: Parameter Optimization

📋 Configuration:
   DE result: results/de_imported/de_result.pkl
   Method: ground-truth
   Ground truth: validated_genes.csv
   Output: results/optimization

📂 Loading DE result...
   ✓ Loaded 15000 genes

📂 Loading ground truth...
   ✓ Loaded 25 validated genes

🔬 Running ground-truth optimization...
   Strategy: grid search (5×5 = 25 evaluations)
   Metric: F1-score
   
   Progress: [========================================] 25/25
   
   Best F1-score: 0.8456
   Found at: alpha=0.0234, lfc=0.75
   ✓ Optimization complete

📊 Optimal Parameters:
============================================================
Alpha (FDR) threshold: 0.0234
LFC threshold: 0.75
Best Score (f1_score): 0.8456
DEG genes at optimal: 1250

Performance Metrics:
  Precision: 0.85 (85% of called genes are true)
  Recall: 0.84 (found 21 of 25 validated genes)
  F1-score: 0.85 (excellent balance)

✅ Results saved to: results/optimization/
```

### Tips for Ground Truth Optimization

1. **Gene ID matching**: Make sure gene IDs match your DE results
   - If DE results use Ensembl IDs (ENSG...), use Ensembl in ground truth
   - If using symbols (TP53, BRCA1), use symbols in ground truth

2. **More genes = better optimization**: 
   - 10 genes: Minimum, may be unstable
   - 20-30 genes: Good, recommended minimum
   - 50+ genes: Excellent, robust optimization

3. **Avoid bias**: Don't only include well-studied genes
   - Bad: Only highly-cited cancer genes
   - Good: Mix of well-known and less-studied genes

4. **Source matters**: Where did your validated genes come from?
   - Best: Your own qRT-PCR/Western blot validation
   - Good: Multiple independent studies
   - Okay: Single literature report
   - Caution: Computational predictions only

---

## Method 2: FDR Control Optimization

**Use when**: You don't have validated genes (most common scenario)

### What is FDR Control Optimization?

This method finds parameters that achieve a specific False Discovery Rate (FDR) target using statistical theory. **No validation data needed!**

### When to Use This

✅ **Perfect for**:
- No validated gene list available (most common)
- Want to control FDR at specific level (e.g., 1%, 5%)
- Exploratory analysis
- General-purpose optimization

✅ **Advantages**:
- No extra data required
- Fast (~3 minutes)
- Statistically principled
- Works for any dataset

### How It Works

Uses the **Storey & Tibshirani (2003) method** to estimate the true FDR:

1. Estimates π₀ (proportion of truly non-DE genes)
2. Calculates FDR for different threshold combinations
3. Finds parameters that achieve your target FDR

**Key insight**: Not all p-values below 0.05 have the same FDR. Some datasets need stricter thresholds to achieve true 5% FDR.

### Step-by-Step Tutorial

#### Step 1: Decide Your Target FDR

Common choices:
- **0.05 (5%)**: Standard for most studies
- **0.01 (1%)**: Strict, for high-confidence genes
- **0.10 (10%)**: Lenient, for exploratory analysis

#### Step 2: Run FDR Control Optimization

**Standard (5% FDR)**:
```bash
raptor optimize \
    -d results/de_imported/de_result.pkl \
    -m fdr-control \
    --fdr-target 0.05
```

**Strict (1% FDR)**:
```bash
raptor optimize \
    -d results/de_imported/de_result.pkl \
    -m fdr-control \
    --fdr-target 0.01 \
    -o results/opt_strict
```

**What happens**:
1. Loads DE results
2. Estimates π₀ (true null proportion)
3. Calculates FDR for different thresholds
4. Finds combination closest to target FDR
5. Returns optimal parameters

#### Step 3: Review Results

Check `results/optimization/optimized_params.yaml`:

```yaml
best_parameters:
  alpha: 0.0089
  lfc_threshold: 0.25
best_score: 0.0011
optimization_method: fdr_control
metric: fdr_distance
target_fdr: 0.01
n_deg_genes: 2150
performance_metrics:
  estimated_fdr: 0.0111
  pi0_estimate: 0.87
  fdr_distance: 0.0011
```

**Interpretation**:
- **Target FDR**: 0.01 (1%)
- **Achieved FDR**: 0.0111 (1.11%)
- **Distance**: 0.0011 (very close to target!)
- **Alpha**: 0.0089 (stricter than default)
- **π₀**: 0.87 (87% of genes are truly non-DE)

### Understanding the Results

**π₀ (Pi-zero) estimate**:
- Proportion of genes that are truly NOT differentially expressed
- Range: 0.5 - 1.0
- Higher = more genes are true nulls

**Typical values**:
- 0.7-0.8: Many DE genes (strong effect)
- 0.8-0.9: Moderate number of DE genes (typical)
- 0.9-1.0: Few DE genes (weak effect or small sample)

**FDR distance**:
- How close achieved FDR is to target
- <0.005: Excellent
- 0.005-0.01: Good
- >0.01: May need more genes for robust estimation

### Example Output

```bash
$ raptor optimize -d de_result.pkl -m fdr-control --fdr-target 0.01

🦖 RAPTOR v2.2.0 - Module 8: Parameter Optimization

📋 Configuration:
   DE result: results/de_imported/de_result.pkl
   Method: fdr-control
   FDR target: 0.01
   Output: results/optimization

📂 Loading DE result...
   ✓ Loaded 15000 genes

🔬 Running fdr-control optimization...
   Estimating π₀ (null proportion)...
   ✓ π₀ estimate: 0.87
   
   Searching for optimal parameters...
   Strategy: grid search (5×5 = 25 evaluations)
   
   Progress: [========================================] 25/25
   
   Best FDR: 0.0111 (target: 0.01)
   Found at: alpha=0.0089, lfc=0.25
   ✓ Optimization complete

📊 Optimal Parameters:
============================================================
Alpha (FDR) threshold: 0.0089
LFC threshold: 0.25
Best Score (fdr_distance): 0.0011
DEG genes at optimal: 2150

FDR Metrics:
  Target FDR: 0.01 (1%)
  Achieved FDR: 0.0111 (1.11%)
  Distance from target: 0.0011 ✓ Excellent
  π₀ estimate: 0.87

✅ Results saved to: results/optimization/
```

### When FDR Control Works Best

✅ **Works well when**:
- Large number of genes (>5,000)
- Many genes tested (~10,000-20,000)
- Clear differential expression signal

⚠️ **May struggle when**:
- Very few genes (<1,000)
- Extremely weak signal (π₀ > 0.95)
- Heavy batch effects

### Tips for FDR Control

1. **Choose target FDR carefully**:
   - 0.01 (1%): Very strict, high confidence
   - 0.05 (5%): Standard, balanced
   - 0.10 (10%): Lenient, exploratory

2. **Check π₀ estimate**:
   - If π₀ > 0.95: Very few DE genes, consider biology
   - If π₀ < 0.5: Unusual, check data quality

3. **Minimum genes**: Need >1,000 genes for robust FDR estimation

4. **Compare with default**: How different is optimized from FDR < 0.05?

---

## Method 3: Stability Optimization

**Use when**: You want robust, reproducible results

### What is Stability Optimization?

This method uses **bootstrap resampling** to find parameters that give the most stable gene selection. Genes that are consistently selected across bootstrap samples are more reliable.

### When to Use This

✅ **Perfect for**:
- Want robust, reproducible results
- Concerned about sampling variability
- Have original count matrix and metadata
- Planning follow-up validation experiments

✅ **Advantages**:
- No validation genes needed
- Measures reproducibility
- Accounts for sampling variability
- Great for small sample sizes

❌ **Don't use when**:
- Don't have count matrix (only DE results)
- Very limited computational time
- Fewer than 6 samples total

### How It Works

**Bootstrap stability selection** (Meinshausen & Bühlmann, 2010):

1. **Bootstrap**: Randomly resample your samples (with replacement)
2. **Re-analyze**: Rerun DE analysis on bootstrap sample
3. **Select genes**: Apply candidate thresholds
4. **Repeat**: Do this 100 times
5. **Measure stability**: How often is each gene selected?
6. **Optimize**: Find thresholds that maximize average stability

**Key insight**: Genes selected >90% of the time are much more reliable than genes selected 50% of the time.

### Required Data

You need **3 files**:

1. **DE results** (from Module 6): `de_result.pkl`
2. **Count matrix**: `counts.csv`
3. **Sample metadata**: `metadata.csv`

**Count matrix format** (`counts.csv`):
```csv
gene_id,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6
GENE1,1250,2340,1890,980,2100,1560
GENE2,450,390,420,510,380,440
GENE3,8900,9200,8700,9500,8800,9100
...
```

**Metadata format** (`metadata.csv`):
```csv
sample_id,condition,batch
Sample1,Control,Batch1
Sample2,Control,Batch1
Sample3,Control,Batch2
Sample4,Treatment,Batch1
Sample5,Treatment,Batch1
Sample6,Treatment,Batch2
```

### Step-by-Step Tutorial

#### Step 1: Prepare Your Data

Make sure you have:
- ✅ Count matrix (raw or normalized counts)
- ✅ Metadata with 'condition' or 'group' column
- ✅ DE results (from Module 6)

#### Step 2: Run Stability Optimization

**Standard (100 bootstrap iterations)**:
```bash
raptor optimize \
    -d results/de_imported/de_result.pkl \
    -m stability \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    -o results/opt_stability
```

**Fast (50 iterations for testing)**:
```bash
raptor optimize \
    -d de_result.pkl \
    -m stability \
    --counts counts.csv \
    --metadata metadata.csv \
    --n-bootstrap 50
```

**Thorough (200 iterations for publication)**:
```bash
raptor optimize \
    -d de_result.pkl \
    -m stability \
    --counts counts.csv \
    --metadata metadata.csv \
    --n-bootstrap 200
```

**What happens**:
1. Loads count matrix and metadata
2. Creates 100 bootstrap samples (resample with replacement)
3. Runs DE analysis on each bootstrap sample
4. Tests different FDR/LFC thresholds
5. Measures how stable gene selection is
6. Returns parameters with highest stability

#### Step 3: Review Results

Check `results/opt_stability/optimized_params.yaml`:

```yaml
best_parameters:
  alpha: 0.0234
  lfc_threshold: 0.75
best_score: 0.9123
optimization_method: stability
metric: stability_score
n_bootstrap: 100
n_deg_genes: 850
performance_metrics:
  stability_score: 0.9123
  mean_selection_probability: 0.91
  high_stability_genes: 775  # >90% selection rate
  medium_stability_genes: 75  # 70-90% selection rate
```

**Interpretation**:
- **Stability score**: 0.91 (excellent!)
- **Average gene**: Selected in 91% of bootstrap samples
- **High stability**: 775 genes selected >90% of time
- **Medium stability**: 75 genes selected 70-90% of time
- **Parameters**: alpha=0.0234, LFC=0.75

### Understanding Stability Score

**Stability score** = Average selection probability across genes

**Interpretation**:
- **0.9-1.0**: Excellent - Very stable gene selection
- **0.8-0.9**: Good - Reasonably stable
- **0.7-0.8**: Moderate - Some variability
- **<0.7**: Low - High variability, consider more stringent thresholds

**Why does this matter?**
- High stability (>0.9) → Genes will likely replicate in validation
- Low stability (<0.7) → Results may not replicate well

### Example Output

```bash
$ raptor optimize -d de_result.pkl -m stability \
    --counts counts.csv --metadata metadata.csv

🦖 RAPTOR v2.2.0 - Module 8: Parameter Optimization

📋 Configuration:
   DE result: de_result.pkl
   Method: stability
   Counts: counts.csv
   Metadata: metadata.csv
   Bootstrap iterations: 100
   Output: results/optimization

📂 Loading data...
   ✓ DE result: 15000 genes
   ✓ Count matrix: 15000 genes × 12 samples
   ✓ Metadata: 12 samples (6 Control, 6 Treatment)

🔬 Running stability optimization...
   Bootstrap progress: [====                ] 20/100 (20%)
   Bootstrap progress: [========            ] 40/100 (40%)
   Bootstrap progress: [============        ] 60/100 (60%)
   Bootstrap progress: [================    ] 80/100 (80%)
   Bootstrap progress: [====================] 100/100 (100%)
   
   Testing parameter combinations...
   Strategy: grid search (5×5 = 25 evaluations)
   
   Progress: [========================================] 25/25
   
   Best stability: 0.9123
   Found at: alpha=0.0234, lfc=0.75
   ✓ Optimization complete

📊 Optimal Parameters:
============================================================
Alpha (FDR) threshold: 0.0234
LFC threshold: 0.75
Best Score (stability_score): 0.9123
DEG genes at optimal: 850

Stability Metrics:
  Average selection probability: 91.2%
  High stability genes (>90%): 775 genes
  Medium stability genes (70-90%): 75 genes
  Low stability genes (<70%): 0 genes

✅ Results saved to: results/optimization/

⏱️  Total time: 18 minutes
```

### Computational Considerations

**Time estimates**:
- 50 bootstrap iterations: ~10 minutes
- 100 iterations (recommended): ~18-25 minutes
- 200 iterations (thorough): ~35-50 minutes

**Speed depends on**:
- Number of samples (more samples = longer)
- Number of genes (more genes = longer)
- Number of bootstrap iterations

**Tip**: Start with 50 iterations to test, then rerun with 100-200 for final results.

### Tips for Stability Optimization

1. **Bootstrap iterations**:
   - 50: Quick test, less robust
   - 100: Standard, good balance
   - 200: Thorough, for publication

2. **Sample requirements**:
   - Minimum: 6 samples (3 per group)
   - Recommended: 10+ samples (5+ per group)
   - Ideal: 20+ samples

3. **Check stability distribution**:
   - Most genes should have >0.8 selection probability
   - If many genes <0.5, results may not be robust

4. **Compare with other methods**:
   - Stability often gives stricter thresholds than defaults
   - This is good! More robust genes.

---

## Method 4: Reproducibility Optimization

**Use when**: You have two independent cohorts

### What is Reproducibility Optimization?

This method finds parameters that maximize gene overlap between **discovery** and **validation** cohorts. Ensures your findings replicate in independent data.

### When to Use This

✅ **Perfect for**:
- You have two independent cohorts (e.g., discovery + validation)
- Multi-center studies
- Two time points from same population
- Independent technical replicates

✅ **Advantages**:
- Direct measure of replicability
- Real-world validation
- Reduces false positives
- Publication-ready

❌ **Don't use when**:
- Only have one cohort
- Cohorts are from very different populations
- Different biological conditions

### Required Data

You need **DE results from TWO cohorts**:

1. **Discovery cohort**: `discovery_de_result.pkl`
2. **Validation cohort**: `validation_de_result.pkl`

**Both must**:
- Come from Module 6 (import-de)
- Test same biological comparison
- Use same gene IDs

### How It Works

**Jaccard index optimization**:

1. Apply candidate thresholds to Discovery cohort → Gene set A
2. Apply same thresholds to Validation cohort → Gene set B
3. Calculate overlap: Jaccard = |A ∩ B| / |A ∪ B|
4. Find thresholds that maximize Jaccard index

**Key insight**: Parameters that work for one cohort should work for another. Maximizing overlap finds generalizable thresholds.

### Step-by-Step Tutorial

#### Step 1: Prepare Two DE Results

**Discovery cohort**:
```bash
raptor import-de \
    -i discovery_deseq2_results.csv \
    -m deseq2 \
    -o results/discovery
```

**Validation cohort**:
```bash
raptor import-de \
    -i validation_deseq2_results.csv \
    -m deseq2 \
    -o results/validation
```

#### Step 2: Run Reproducibility Optimization

```bash
raptor optimize \
    -d results/discovery/de_result.pkl \
    -m reproducibility \
    --cohort2 results/validation/de_result.pkl \
    -o results/opt_reproducibility
```

**What happens**:
1. Loads both DE results
2. Tests different FDR/LFC thresholds
3. Applies thresholds to both cohorts
4. Calculates gene overlap (Jaccard index)
5. Returns parameters with best overlap

#### Step 3: Review Results

Check `results/opt_reproducibility/optimized_params.yaml`:

```yaml
best_parameters:
  alpha: 0.0145
  lfc_threshold: 0.50
best_score: 0.7234
optimization_method: reproducibility
metric: jaccard_index
n_deg_genes: 1420
performance_metrics:
  jaccard_index: 0.7234
  overlap_count: 1028
  discovery_only: 392
  validation_only: 445
  cohort1_total: 1420
  cohort2_total: 1473
```

**Interpretation**:
- **Jaccard index**: 0.72 (good overlap!)
- **Shared genes**: 1,028 genes in both cohorts
- **Discovery only**: 392 genes (27%)
- **Validation only**: 445 genes (30%)
- **Parameters**: alpha=0.0145, LFC=0.50

### Understanding Jaccard Index

**Jaccard Index** = Overlap / Union = |A ∩ B| / |A ∪ B|

**Range**: 0 (no overlap) to 1 (perfect overlap)

**Interpretation**:
- **0.8-1.0**: Excellent reproducibility
- **0.6-0.8**: Good reproducibility
- **0.4-0.6**: Moderate reproducibility
- **<0.4**: Poor reproducibility (check if cohorts are comparable)

**Example**:
```
Discovery cohort:  1420 DEGs → Gene set A
Validation cohort: 1473 DEGs → Gene set B

Overlap (A ∩ B):   1028 genes
Union (A ∪ B):     1420 + 1473 - 1028 = 1865 genes

Jaccard = 1028 / 1865 = 0.55 (moderate)
```

### Example Output

```bash
$ raptor optimize -d discovery.pkl -m reproducibility --cohort2 validation.pkl

🦖 RAPTOR v2.2.0 - Module 8: Parameter Optimization

📋 Configuration:
   DE result (discovery): results/discovery/de_result.pkl
   Method: reproducibility
   Cohort 2 (validation): results/validation/de_result.pkl
   Output: results/optimization

📂 Loading DE results...
   ✓ Discovery cohort: 15000 genes, 1850 significant (default)
   ✓ Validation cohort: 15000 genes, 1920 significant (default)
   ✓ Gene overlap: 14500 genes in both datasets

🔬 Running reproducibility optimization...
   Testing parameter combinations...
   Strategy: grid search (5×5 = 25 evaluations)
   
   Progress: [========================================] 25/25
   
   Best Jaccard: 0.7234
   Found at: alpha=0.0145, lfc=0.50
   ✓ Optimization complete

📊 Optimal Parameters:
============================================================
Alpha (FDR) threshold: 0.0145
LFC threshold: 0.50
Best Score (jaccard_index): 0.7234
DEG genes at optimal: 1420 (discovery), 1473 (validation)

Reproducibility Metrics:
  Jaccard index: 0.7234 ✓ Good
  Shared genes: 1028 (72% of discovery genes)
  Discovery-only genes: 392 (28%)
  Validation-only genes: 445 (30%)
  
Venn Diagram:
  Discovery:  [=========1420=========]
  Validation:       [========1473========]
  Overlap:          [===1028===]

✅ Results saved to: results/optimization/
```

### Venn Diagram Interpretation

```
     Discovery (1420)      Validation (1473)
          ┌─────────────┬─────────────┐
          │             │             │
          │     392     │    1028     │    445
          │  Discovery  │   Shared    │ Validation
          │    only     │             │    only
          │             │             │
          └─────────────┴─────────────┘
```

**What this means**:
- **1,028 genes**: Replicate in both cohorts (high confidence!)
- **392 genes**: Only in discovery (may be false positives or cohort-specific)
- **445 genes**: Only in validation (may have been missed in discovery)

### Tips for Reproducibility Optimization

1. **Cohort requirements**:
   - Same biological comparison (e.g., both "Treatment vs Control")
   - Same gene annotation (Ensembl, Symbol, etc.)
   - Comparable sample sizes (within 2-fold)

2. **What's a good Jaccard index?**:
   - >0.7: Excellent, cohorts agree well
   - 0.5-0.7: Good, reasonable overlap
   - 0.3-0.5: Moderate, cohorts somewhat different
   - <0.3: Poor, check if cohorts are comparable

3. **If overlap is low (<0.4)**:
   - Check if biological conditions match
   - Verify gene ID matching
   - Consider batch effects
   - Check sample quality in both cohorts

4. **Use shared genes**:
   - The 1,028 shared genes are your high-confidence set
   - Focus downstream analysis on these
   - Cohort-specific genes need extra validation

---

## Understanding Your Results

### Output Files Explained

After optimization, you'll have **4 files** in your output directory:

#### 1. `optimized_params.yaml`

**What it is**: Human-readable summary of optimal parameters

**Example**:
```yaml
best_parameters:
  alpha: 0.0123
  lfc_threshold: 0.50
best_score: 0.8456
optimization_method: ground_truth
search_strategy: grid
metric: f1_score
n_iterations: 25
n_deg_genes: 1250
performance_metrics:
  precision: 0.85
  recall: 0.84
  f1_score: 0.85
timestamp: '2026-03-09T14:30:00'
```

**How to use it**:
```bash
# Extract optimal values for next steps
ALPHA=$(grep "alpha:" optimized_params.yaml | awk '{print $2}')
LFC=$(grep "lfc_threshold:" optimized_params.yaml | awk '{print $2}')

echo "Use alpha=$ALPHA and lfc=$LFC for filtering"
```

#### 2. `deg_genes.csv`

**What it is**: Filtered genes using optimal thresholds

**Example**:
```csv
gene_id,log2_fold_change,p_value,adjusted_p_value,is_significant,direction,base_mean
TP53,2.45,0.0001,0.0012,True,up,5234.5
BRCA1,-1.87,0.0002,0.0018,True,down,3456.2
EGFR,1.92,0.0003,0.0023,True,up,8901.3
...
```

**How to use it**:
- Import into Excel/R for further analysis
- Use for pathway enrichment (KEGG, GO, Reactome)
- Validate top genes with qRT-PCR
- Generate heatmaps

#### 3. `optimization_result.pkl`

**What it is**: Complete Python object with all results

**How to use it**:
```python
from raptor.parameter_optimization import OptimizationResult

# Load results
result = OptimizationResult.load('optimization_result.pkl')

# Access any field
print(f"Best alpha: {result.best_parameters['alpha']}")
print(f"Best LFC: {result.best_parameters['lfc_threshold']}")
print(f"Score: {result.best_score}")

# Get convergence history
history = result.convergence_history
scores = [h['score'] for h in history]

# Plot optimization
import matplotlib.pyplot as plt
plt.plot(scores)
plt.xlabel('Iteration')
plt.ylabel('Score')
plt.title('Optimization Convergence')
plt.show()
```

#### 4. `convergence_history.json`

**What it is**: Complete optimization trajectory

**Example**:
```json
[
  {
    "iteration": 1,
    "alpha": 0.001,
    "lfc_threshold": 0.0,
    "score": 0.45,
    "n_deg_genes": 4500
  },
  {
    "iteration": 2,
    "alpha": 0.001,
    "lfc_threshold": 0.5,
    "score": 0.52,
    "n_deg_genes": 2800
  },
  ...
  {
    "iteration": 25,
    "alpha": 0.0123,
    "lfc_threshold": 0.50,
    "score": 0.8456,
    "n_deg_genes": 1250
  }
]
```

**How to use it**:
- Verify optimization converged
- Check if local optimum vs global optimum
- Understand parameter space
- Create visualizations

### Interpreting Optimal Parameters

#### Alpha (FDR Threshold)

**What you might see**:
```yaml
alpha: 0.0123
```

**Interpretation**:
- **0.001-0.01**: Very strict (0.1-1% FDR)
  - Few false positives
  - High confidence genes
  - May miss some true positives

- **0.01-0.05**: Moderate (1-5% FDR)
  - Balanced approach
  - Standard for most studies
  - Good for exploration

- **0.05-0.20**: Lenient (5-20% FDR)
  - More genes called
  - Higher false positive rate
  - Good for hypothesis generation

**Comparison with default**:
```
Default:   FDR < 0.05  →  1,850 genes
Optimized: FDR < 0.0123 →  1,250 genes (more stringent!)
```

#### LFC Threshold

**What you might see**:
```yaml
lfc_threshold: 0.50
```

**Interpretation**:
- **0.0**: No fold-change filter
  - Statistical significance only
  - May include small biological changes

- **0.5**: ~1.4-fold change minimum
  - Balances stats and biology
  - Good for most applications

- **1.0**: 2-fold change minimum
  - Strong biological effect
  - Standard for many studies

- **1.5-2.0**: 3-4 fold change
  - Very large effects only
  - May miss important genes

**Biological meaning**:
```
LFC = 0.5  →  1.4-fold change (2^0.5 = 1.41)
LFC = 1.0  →  2.0-fold change (2^1.0 = 2.00)
LFC = 1.5  →  2.8-fold change (2^1.5 = 2.83)
LFC = 2.0  →  4.0-fold change (2^2.0 = 4.00)
```

### Comparing with Default Thresholds

**Always check**: How do optimized thresholds compare to defaults?

**Example comparison**:
```
                 Default    Optimized   Difference
Alpha (FDR):     0.05       0.0123      -75% (stricter)
LFC threshold:   0.0        0.50        +0.5 (more stringent)
DEG count:       1,850      1,250       -600 genes (-32%)
```

**What this tells you**:
- Optimized is **more stringent** than default
- Default thresholds would give **600 false positives**
- Your optimal gene set is more **reliable** for validation

### Checking Optimization Quality

#### 1. Did optimization converge?

Check `convergence_history.json`:

```python
import json
import matplotlib.pyplot as plt

with open('convergence_history.json') as f:
    history = json.load(f)

scores = [h['score'] for h in history]
plt.plot(scores)
plt.xlabel('Iteration')
plt.ylabel('Score')
plt.title('Optimization Convergence')
plt.show()
```

**Good convergence**:
```
Score ↑
  |     .....****
  |   ....*
  | ...*
  |..*
  +-------------→ Iteration
```

**Poor convergence**:
```
Score ↑
  | .*..*..*.*
  |.*..*..*..*
  |*..*..*.*..
  |..*..*.*...
  +-------------→ Iteration
```

If poorly converged:
- Try more iterations
- Try different search strategy
- Check if enough data for optimization

#### 2. Is the score reasonable?

**Method-specific expectations**:

| Method | Metric | Good Score |
|--------|--------|-----------|
| Ground Truth | F1-score | >0.7 |
| FDR Control | FDR distance | <0.01 |
| Stability | Stability score | >0.8 |
| Reproducibility | Jaccard index | >0.6 |

#### 3. How many genes are called?

**Typical ranges**:
- Too few (<100): May be too strict, missing signal
- Reasonable (100-2000): Good for most analyses
- Many (2000-5000): Broad discovery, need validation
- Too many (>5000): May be too lenient, high FDR

**Context matters**:
- Small samples (n=3/group): Expect 100-500 genes
- Moderate samples (n=10/group): Expect 500-2000 genes
- Large samples (n=50/group): Expect 1000-5000 genes

---

## Choosing the Right Method

### Decision Framework

Use this flowchart to choose the best method for your study:

```
START

Q1: Do you have validated genes from literature/experiments?
├─ YES (10+ genes) → GROUND TRUTH
│                     ↓
│                     Provides: Gold standard validation
│                     Time: ~5 minutes
│                     Confidence: Highest
│
└─ NO → Continue to Q2

Q2: Do you have two independent cohorts?
├─ YES → REPRODUCIBILITY
│         ↓
│         Provides: Real-world replication
│         Time: ~5 minutes
│         Confidence: High
│
└─ NO → Continue to Q3

Q3: Do you have count matrix + metadata?
├─ YES → STABILITY
│         ↓
│         Provides: Bootstrap robustness
│         Time: ~20 minutes
│         Confidence: High
│
└─ NO → FDR CONTROL
          ↓
          Provides: Statistical FDR control
          Time: ~3 minutes
          Confidence: Good

END
```

### Method Comparison Table

| Criterion | Ground Truth | FDR Control | Stability | Reproducibility |
|-----------|-------------|-------------|-----------|----------------|
| **Data Required** | Validated genes | None | Counts + metadata | 2 cohorts |
| **Time** | Fast (5 min) | Fast (3 min) | Slow (20 min) | Fast (5 min) |
| **Confidence** | ★★★★★ | ★★★☆☆ | ★★★★☆ | ★★★★☆ |
| **Robustness** | ★★★★☆ | ★★★☆☆ | ★★★★★ | ★★★★☆ |
| **Generalizability** | ★★★☆☆ | ★★★★☆ | ★★★★☆ | ★★★★★ |
| **Best For** | Hypothesis-driven | General purpose | Robustness | Replication |

### Real-World Scenarios

#### Scenario 1: Cancer Biomarker Study

**Situation**:
- You have RNA-seq from 30 tumor samples (15 responders, 15 non-responders)
- Found 50 known cancer genes in literature
- Planning qRT-PCR validation of top genes

**Recommended**:
```bash
raptor optimize -d de_result.pkl -m ground-truth -g cancer_genes.csv
```

**Why**: Gold standard with validated genes, best for identifying biomarkers

---

#### Scenario 2: Drug Response Exploratory Study

**Situation**:
- First study of new drug in this disease
- 20 samples (10 treated, 10 control)
- No prior validated genes
- Have count matrix

**Recommended**:
```bash
raptor optimize -d de_result.pkl -m stability \
    --counts counts.csv --metadata metadata.csv
```

**Why**: No validated genes, stability ensures robust discovery

---

#### Scenario 3: Small Pilot Study

**Situation**:
- Only 6 samples (3 per group)
- Limited budget
- Quick preliminary results needed

**Recommended**:
```bash
raptor optimize -d de_result.pkl -m fdr-control --fdr-target 0.01
```

**Why**: Fast, no extra data needed, strict FDR for small sample

---

#### Scenario 4: Multi-Center Clinical Trial

**Situation**:
- Discovery cohort: 50 samples (Site A)
- Validation cohort: 45 samples (Site B)
- Need results that replicate

**Recommended**:
```bash
raptor optimize -d discovery.pkl -m reproducibility --cohort2 validation.pkl
```

**Why**: Direct measure of cross-site replication

---

### Can I Use Multiple Methods?

**Yes!** It's actually recommended to try 2-3 methods and compare:

```bash
# Method 1: Ground truth (if you have validated genes)
raptor optimize -d de_result.pkl -m ground-truth -g validated.csv -o opt_gt

# Method 2: FDR control
raptor optimize -d de_result.pkl -m fdr-control -o opt_fdr

# Method 3: Stability
raptor optimize -d de_result.pkl -m stability \
    --counts counts.csv --metadata metadata.csv -o opt_stability
```

**Then compare**:
```
Method        Alpha    LFC    DEGs   Score
Ground Truth  0.0123   0.50   1250   0.85 (F1)
FDR Control   0.0089   0.25   2150   0.001 (FDR dist)
Stability     0.0234   0.75   850    0.91 (Stability)
```

**Choose final parameters based on**:
- Which method fits your study design best
- Consistency across methods
- Biological knowledge
- Downstream application (validation vs discovery)

---

## Advanced Options

### Search Strategies

Module 8 supports **3 search strategies** to find optimal parameters:

#### 1. Grid Search (Default)

**What it does**: Exhaustively tests all combinations in a grid

**Advantages**:
- ✅ Guaranteed to find global optimum
- ✅ Systematic, reproducible
- ✅ Good for visualization

**Disadvantages**:
- ⚠️ Slower for fine grids
- ⚠️ Number of evaluations grows quickly

**When to use**: Default choice, works well for most cases

**Example**:
```bash
raptor optimize -d de_result.pkl -m fdr-control
# Uses default 5×5 grid = 25 evaluations
```

**Customize grid**:
```bash
raptor optimize -d de_result.pkl -m fdr-control --grid-points 10
# Uses 10×10 grid = 100 evaluations
```

**Grid layout**:
```
LFC ↑
    |
2.0 | • • • • •
1.5 | • • • • •
1.0 | • • • • •  ← Each • is one parameter combination
0.5 | • • • • •
0.0 | • • • • •
    +------------→ Alpha
    0.001 0.05 0.10 0.15 0.20
```

---

#### 2. Random Search

**What it does**: Randomly samples parameter combinations

**Advantages**:
- ✅ Faster than grid search
- ✅ Often finds good solutions quickly
- ✅ Better for high-dimensional spaces

**Disadvantages**:
- ⚠️ Not guaranteed to find exact optimum
- ⚠️ Results vary slightly between runs

**When to use**:
- Large parameter spaces
- Time constraints
- Exploratory analysis

**Example**:
```bash
raptor optimize -d de_result.pkl -m fdr-control \
    --strategy random --n-iterations 50
```

**Scientific basis**: Bergstra & Bengio (2012) showed random search often outperforms grid search in practice.

---

#### 3. Differential Evolution

**What it does**: Population-based global optimization algorithm

**Advantages**:
- ✅ Good at avoiding local optima
- ✅ Adaptive search
- ✅ Robust to complex landscapes

**Disadvantages**:
- ⚠️ More iterations needed
- ⚠️ Slightly slower per iteration

**When to use**:
- Complex optimization landscape
- Want global optimum guarantee
- Have computational time

**Example**:
```bash
raptor optimize -d de_result.pkl -m ground-truth \
    --strategy differential_evolution --n-iterations 100
```

**How it works**:
1. Initialize population of parameter sets
2. Create new candidates by combining existing ones
3. Test candidates and keep best
4. Repeat until convergence

---

### Strategy Comparison

| Strategy | Speed | Optimality | Reproducibility | Best Use |
|----------|-------|-----------|----------------|----------|
| **Grid** | Slow | Guaranteed | 100% | Default, visualization |
| **Random** | Fast | Good | Variable | Time limits, exploration |
| **Differential Evolution** | Medium | Best | High | Complex problems |

### Choosing a Strategy

**For most users**: Use default **grid search**
```bash
raptor optimize -d de_result.pkl -m fdr-control
```

**If optimization is slow**: Try **random search**
```bash
raptor optimize -d de_result.pkl -m stability \
    --strategy random --n-iterations 50
```

**For publication**: Use **differential evolution** for best results
```bash
raptor optimize -d de_result.pkl -m ground-truth \
    --strategy differential_evolution --n-iterations 200
```

---

### Custom Parameter Ranges

**Default ranges**:
```python
alpha: 0.001 - 0.20
lfc_threshold: 0.0 - 2.0
```

**To customize** (Python API):
```python
from raptor.parameter_optimization import ParameterSpace, GroundTruthOptimizer

# Define custom ranges
param_space = ParameterSpace(
    alpha_min=0.001,
    alpha_max=0.10,    # Only search up to 10% FDR
    lfc_min=0.5,       # Minimum 1.4-fold change
    lfc_max=1.5        # Maximum 2.8-fold change
)

# Create optimizer with custom space
optimizer = GroundTruthOptimizer(
    de_result=de_df,
    ground_truth=gt_df,
    parameter_space=param_space
)

# Optimize
result = optimizer.optimize(strategy='grid')
```

**When to customize**:
- You know biologically relevant fold-change ranges
- Want stricter FDR bounds (e.g., 0.001-0.05 instead of 0.001-0.20)
- Have specific requirements (e.g., must be >1.5-fold)

---

### Bootstrap Iterations (Stability Only)

**Default**: 100 iterations

**Customize**:
```bash
# Fast test (50 iterations)
raptor optimize -d de_result.pkl -m stability \
    --counts counts.csv --metadata metadata.csv \
    --n-bootstrap 50

# Standard (100 iterations)
raptor optimize -d de_result.pkl -m stability \
    --counts counts.csv --metadata metadata.csv \
    --n-bootstrap 100

# Thorough (200 iterations for publication)
raptor optimize -d de_result.pkl -m stability \
    --counts counts.csv --metadata metadata.csv \
    --n-bootstrap 200
```

**Recommendations**:
- **50**: Quick test, ~10 minutes
- **100**: Standard, good balance, ~20 minutes
- **200**: Thorough, for publication, ~40 minutes

**More is better**, but diminishing returns after 100-200.

---

## Integration with Your Workflow

### Complete RAPTOR Workflow

Here's how Module 8 fits into a complete RNA-seq analysis:

#### Step 1: Quality Control (Module 2)
```bash
raptor qc \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    --group condition \
    -o results/qc
```

#### Step 2: Data Profiling (Module 3)
```bash
raptor profile \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    -o results/profile
```

#### Step 3: Pipeline Recommendation (Module 4)
```bash
raptor recommend \
    --profile results/profile/data_profile.json \
    -o results/recommendation
```

#### Step 4: Run DE Analysis (External)
```r
# Run DESeq2 in R (example)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "deseq2_results.csv")
```

#### Step 5: Import DE Results (Module 6)
```bash
raptor import-de \
    -i deseq2_results.csv \
    -m deseq2 \
    -o results/de_imported
```

#### Step 6: Optimize Parameters (Module 8) ⭐ THIS MODULE
```bash
raptor optimize \
    -d results/de_imported/de_result.pkl \
    -m fdr-control \
    --fdr-target 0.01 \
    -o results/optimization
```

#### Step 7: Apply Optimal Thresholds
```bash
# Extract optimal parameters
ALPHA=$(grep "alpha:" results/optimization/optimized_params.yaml | awk '{print $2}')
LFC=$(grep "lfc_threshold:" results/optimization/optimized_params.yaml | awk '{print $2}')

# The optimized genes are already in deg_genes.csv!
# Or re-filter with optimal thresholds if needed
```

#### Step 8: Ensemble Analysis (Module 9)
```bash
raptor ensemble \
    --methods fisher brown rra \
    --deseq2 results/optimization/deg_genes.csv \
    --edger results/edger/deg_genes.csv \
    --limma results/limma/deg_genes.csv \
    -o results/ensemble
```

#### Step 9: Biomarker Discovery (Module 10)
```bash
raptor biomarkers \
    --de-result results/ensemble/ensemble_result.pkl \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    -o results/biomarkers
```

---

### Using Optimized Parameters

After optimization, you have **optimized thresholds**. Here's how to use them:

#### Option 1: Use the Pre-Filtered Genes

**Easiest**: Module 8 already filtered genes for you!

```bash
# The file deg_genes.csv contains genes at optimal thresholds
head results/optimization/deg_genes.csv

# Use directly for downstream analysis
```

**For pathway enrichment**:
```r
library(clusterProfiler)

# Load optimized genes
genes <- read.csv("results/optimization/deg_genes.csv")
gene_list <- genes$gene_id

# Run enrichment
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH"
)
```

#### Option 2: Re-Import with Optimal Thresholds

If you want to **re-run DE import** with optimal thresholds:

```bash
# Get optimal values
ALPHA=0.0123
LFC=0.50

# Re-import with optimal thresholds
raptor import-de \
    -i deseq2_results.csv \
    -m deseq2 \
    --fdr-threshold $ALPHA \
    --lfc-threshold $LFC \
    -o results/de_optimized
```

#### Option 3: Filter in R/Python

**In R**:
```r
# Load original DE results
res <- read.csv("deseq2_results.csv")

# Apply optimal thresholds
alpha <- 0.0123
lfc_threshold <- 0.50

filtered <- res[
    res$padj < alpha & abs(res$log2FoldChange) > lfc_threshold,
]

# Save filtered results
write.csv(filtered, "deseq2_optimized.csv")
```

**In Python**:
```python
import pandas as pd

# Load DE results
res = pd.read_csv("deseq2_results.csv")

# Apply optimal thresholds
alpha = 0.0123
lfc_threshold = 0.50

filtered = res[
    (res['padj'] < alpha) & 
    (abs(res['log2FoldChange']) > lfc_threshold)
]

# Save
filtered.to_csv("deseq2_optimized.csv")
```

---

### Integration with External Tools

#### Pathway Enrichment (clusterProfiler)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# Load optimized genes
genes <- read.csv("results/optimization/deg_genes.csv")
gene_list <- genes$gene_id

# GO enrichment
ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)

# KEGG enrichment
kegg <- enrichKEGG(
    gene = gene_list,
    organism = 'hsa',
    pAdjustMethod = "BH"
)

# Plot
dotplot(ego, showCategory=20)
dotplot(kegg, showCategory=20)
```

#### Gene Set Enrichment (GSEA)

```r
library(fgsea)

# Load optimized genes with LFC
genes <- read.csv("results/optimization/deg_genes.csv")

# Create ranked list
ranked_genes <- setNames(
    genes$log2_fold_change,
    genes$gene_id
)
ranked_genes <- sort(ranked_genes, decreasing=TRUE)

# Run GSEA
fgseaRes <- fgsea(
    pathways = hallmark_pathways,
    stats = ranked_genes,
    minSize = 15,
    maxSize = 500
)
```

#### Visualization (DESeq2)

```r
library(DESeq2)
library(ggplot2)

# Load DESeq2 object and optimized genes
dds <- readRDS("dds.rds")
genes <- read.csv("results/optimization/deg_genes.csv")

# Transform counts
vsd <- vst(dds, blind=FALSE)

# Heatmap of optimized genes
library(pheatmap)
pheatmap(
    assay(vsd)[genes$gene_id, ],
    scale = "row",
    show_rownames = FALSE
)

# MA plot with optimal thresholds
plotMA(dds, alpha=0.0123, ylim=c(-5,5))
abline(h=c(-0.5, 0.5), col="red", lty=2)  # LFC threshold
```

---

## Best Practices

### Before Running Optimization

1. ✅ **Check data quality first** (Module 2)
   - Remove outliers
   - Check batch effects
   - Ensure adequate sample size

2. ✅ **Profile your data** (Module 3)
   - Understand sample size
   - Check biological coefficient of variation
   - Verify group sizes

3. ✅ **Run preliminary DE analysis**
   - Use recommended pipeline (Module 4)
   - Check if you get reasonable number of genes
   - Verify DE analysis completed successfully

4. ✅ **Import DE results properly** (Module 6)
   - Ensure all required columns present
   - Check gene ID format matches
   - Verify significance columns correct

### During Optimization

1. ✅ **Choose appropriate method**
   - Use decision tree in this guide
   - Consider available data
   - Match method to study goals

2. ✅ **Start with defaults**
   - Use default search strategy (grid)
   - Use default iterations
   - Only customize if needed

3. ✅ **Monitor progress**
   - Check that optimization runs to completion
   - Watch for error messages
   - Verify reasonable runtime

4. ✅ **Save all outputs**
   - Keep all 4 output files
   - Document which method you used
   - Record optimization settings

### After Optimization

1. ✅ **Review results carefully**
   - Check optimal parameters are reasonable
   - Compare with defaults
   - Verify optimization converged

2. ✅ **Validate findings**
   - Check top genes make biological sense
   - Look for known genes in your system
   - Compare with literature

3. ✅ **Document your choice**
   - Record which method used
   - Note why you chose that method
   - Save optimized parameters for paper

4. ✅ **Use consistently**
   - Apply same thresholds to all DE results
   - Don't optimize separately for each comparison
   - Maintain consistency for reproducibility

### For Publication

1. ✅ **Use conservative method**
   - Ground truth if possible (gold standard)
   - Stability for robustness
   - Reproducibility if you have validation cohort

2. ✅ **Report thoroughly**
   ```
   "Parameter optimization was performed using RAPTOR v2.2.0 
   Module 8 with the ground truth method. A list of 25 
   validated differentially expressed genes from [citation] 
   was used as ground truth. Grid search optimization 
   identified optimal thresholds of FDR < 0.0123 and 
   |log2FC| > 0.50, achieving an F1-score of 0.85."
   ```

3. ✅ **Include supplementary info**
   - Upload `optimized_params.yaml` as supplementary file
   - Include convergence plot
   - Provide gene lists (deg_genes.csv)

4. ✅ **Compare with defaults**
   - Show both default and optimized results
   - Justify why optimized is better
   - Demonstrate robustness

### Common Mistakes to Avoid

❌ **Don't**:
1. Skip data QC before optimization
2. Optimize separately for each comparison
3. Cherry-pick method based on results
4. Ignore biological knowledge
5. Use too strict thresholds (0 genes)
6. Use too lenient thresholds (all genes)
7. Forget to document settings
8. Change thresholds after seeing results

✅ **Do**:
1. QC data first
2. Optimize once, apply consistently
3. Choose method before seeing results
4. Sanity-check with biology
5. Aim for reasonable gene counts (100-2000)
6. Balance stringency and sensitivity
7. Document everything
8. Pre-register analysis plan if possible

---

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: "ImportError: cannot import module"

**Error**:
```
ImportError: cannot import name 'optimize_with_ground_truth' from 'raptor.parameter_optimization'
```

**Solutions**:
1. Check RAPTOR is installed:
   ```bash
   pip list | grep raptor
   ```

2. Reinstall if needed:
   ```bash
   pip install --upgrade raptor
   ```

3. Check Python version (requires 3.8+):
   ```bash
   python --version
   ```

---

#### Issue 2: "Ground truth file has no genes in DE results"

**Error**:
```
ValidationError: No overlap between ground truth genes and DE results
```

**Causes**:
- Gene ID mismatch (symbols vs Ensembl IDs)
- Wrong organism annotation
- Ground truth genes not in count matrix

**Solutions**:
1. Check gene ID format:
   ```bash
   # In DE results
   head -1 deseq2_results.csv
   
   # In ground truth
   head validated_genes.csv
   ```

2. Convert gene IDs if needed:
   ```r
   library(biomaRt)
   
   # Convert symbols to Ensembl
   mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
   genes_converted <- getBM(
       attributes=c('hgnc_symbol', 'ensembl_gene_id'),
       filters='hgnc_symbol',
       values=ground_truth_genes,
       mart=mart
   )
   ```

---

#### Issue 3: "Optimization found 0 genes"

**Error**:
```
Warning: Optimal parameters result in 0 DEG genes
```

**Causes**:
- Very strict thresholds
- Weak differential expression signal
- Small sample size

**Solutions**:
1. Check original DE results:
   ```bash
   # How many genes at default thresholds?
   raptor import-de -i results.csv -m deseq2
   ```

2. Try more lenient target:
   ```bash
   # Instead of 1% FDR, try 5%
   raptor optimize -d de_result.pkl -m fdr-control --fdr-target 0.05
   ```

3. Consider biology:
   - Is weak signal expected?
   - Do you need more samples?
   - Are there strong batch effects?

---

#### Issue 4: "Optimization score is very low"

**Example**:
```
Best Score (f1_score): 0.23
```

**Interpretation by method**:

**Ground Truth**:
- F1 < 0.5: Poor overlap with validated genes
- Possible causes:
  - Validated genes not truly DE in your data
  - Sample size too small
  - Different biological context

**FDR Control**:
- FDR distance > 0.05: Can't achieve target FDR
- Possible causes:
  - Not enough genes for robust π₀ estimation
  - Extreme FDR target (too strict or lenient)

**Stability**:
- Stability < 0.6: High variability
- Possible causes:
  - Small sample size (n<8)
  - High biological variability
  - Batch effects

**Reproducibility**:
- Jaccard < 0.3: Poor overlap between cohorts
- Possible causes:
  - Different populations
  - Batch effects
  - Low power in one cohort

**Solutions**:
1. Review data quality
2. Check if expectations are realistic
3. Try different optimization method
4. Consider if more data is needed

---

#### Issue 5: "Stability optimization is too slow"

**Problem**: Takes >1 hour to complete

**Solutions**:
1. Reduce bootstrap iterations:
   ```bash
   raptor optimize -d de_result.pkl -m stability \
       --counts counts.csv --metadata metadata.csv \
       --n-bootstrap 50  # Instead of default 100
   ```

2. Use random search:
   ```bash
   raptor optimize -d de_result.pkl -m stability \
       --counts counts.csv --metadata metadata.csv \
       --strategy random --n-iterations 20
   ```

3. Use parallel processing (if available):
   ```bash
   # Set environment variable for parallel processing
   export RAPTOR_N_JOBS=4
   raptor optimize ...
   ```

---

#### Issue 6: "Optimal parameters same as defaults"

**Result**:
```
Optimal: alpha=0.05, lfc=0.0  (same as defaults!)
```

**Interpretation**:
- Defaults are actually optimal for your data (good!)
- OR optimization method didn't have enough signal

**What to do**:
1. This is fine! Defaults can be optimal
2. Document that you tested and confirmed
3. Consider trying different method to verify
4. Check optimization score - if high, defaults are truly optimal

---

#### Issue 7: "ValueError: Not enough genes"

**Error**:
```
ValueError: DE results must have at least 1000 genes for FDR estimation
```

**Cause**: Too few genes for robust optimization

**Solutions**:
1. Use less stringent pre-filtering
2. Increase gene count in DE analysis
3. Use different optimization method (ground truth if possible)

---

### Getting Help

If you encounter issues not covered here:

1. **Check logs**: Look for detailed error messages
   ```bash
   raptor optimize ... --verbose
   ```

2. **Minimal example**: Try with test data
   ```bash
   # RAPTOR includes test data
   raptor optimize -d test_data/de_result.pkl -m fdr-control
   ```

3. **GitHub Issues**: Report bugs at https://github.com/AyehBlk/RAPTOR/issues

4. **Include details**:
   - RAPTOR version
   - Full command used
   - Complete error message
   - Data summary (n genes, n samples)

---

## Frequently Asked Questions

### General Questions

**Q1: Do I have to optimize parameters?**

A: No, but it's highly recommended. Default thresholds (FDR < 0.05, |LFC| > 0) work for many studies, but optimization can:
- Reduce false positives
- Improve reproducibility
- Provide scientific justification
- Better match your specific data

**Q2: Which method should I use if I'm unsure?**

A: Start with **FDR Control**:
```bash
raptor optimize -d de_result.pkl -m fdr-control
```

It requires no extra data and works well for most cases. If results seem unreasonable, try **Stability** (if you have counts) or get validated genes for **Ground Truth**.

**Q3: How long does optimization take?**

A: Typical times:
- Ground Truth: 3-5 minutes
- FDR Control: 2-3 minutes  
- Stability: 15-30 minutes (depends on bootstrap iterations)
- Reproducibility: 3-5 minutes

**Q4: Can I use Module 8 with non-RNA-seq data?**

A: Yes! Module 8 works with any differential analysis results that have:
- Gene/feature IDs
- Log fold changes
- P-values
- Adjusted p-values (FDR)

Examples: proteomics, metabolomics, microarray, ChIP-seq peaks.

**Q5: What if my optimal thresholds are very different from defaults?**

A: This is actually common and good! It means:
- Defaults weren't appropriate for your data
- Optimization helped find better thresholds
- Your results will be more reliable

Document the difference and justify in your paper.

---

### Method-Specific Questions

**Q6: How many validated genes do I need for Ground Truth?**

A:
- **Minimum**: 10 genes (absolute minimum, may be unstable)
- **Recommended**: 20-30 genes (good optimization)
- **Ideal**: 50+ genes (robust, reliable)

More genes = better optimization. But even 10-15 well-validated genes are better than no validation.

**Q7: What's a good F1-score for Ground Truth?**

A:
- **0.9-1.0**: Excellent (rare, very high confidence)
- **0.8-0.9**: Very good (publication quality)
- **0.7-0.8**: Good (reasonable balance)
- **0.6-0.7**: Acceptable (may need more validation)
- **<0.6**: Poor (check validated genes and data)

**Q8: Why is my π₀ estimate so high (>0.9)?**

A: High π₀ (>0.9) means most genes are NOT differentially expressed. This could mean:
- Weak biological effect
- Small sample size (low power)
- Appropriate stringent thresholds needed

This is fine! Not all experiments show massive differential expression.

**Q9: How many bootstrap iterations should I use for Stability?**

A:
- **50**: Quick test, good for exploring (~10 min)
- **100**: Standard, recommended for most studies (~20 min)
- **200**: Thorough, for publication (~40 min)
- **500+**: Overkill, diminishing returns

100 iterations is the sweet spot for most applications.

**Q10: What's a good Jaccard index for Reproducibility?**

A:
- **0.8-1.0**: Excellent replication
- **0.6-0.8**: Good replication (typical for well-designed studies)
- **0.4-0.6**: Moderate replication (check for batch effects)
- **0.2-0.4**: Poor replication (cohorts may be too different)
- **<0.2**: Very poor (investigate discrepancies)

Even 0.5-0.6 is reasonable for cross-site or cross-platform studies.

---

### Technical Questions

**Q11: Can I optimize multiple comparisons together?**

A: No, optimize each comparison separately:

```bash
# Comparison 1: Treatment vs Control
raptor optimize -d treatment_vs_control.pkl -m fdr-control -o opt_trt_vs_ctrl

# Comparison 2: Time2 vs Time1
raptor optimize -d time2_vs_time1.pkl -m fdr-control -o opt_time2_vs_time1
```

Each comparison may have different optimal thresholds.

**Q12: Should I optimize before or after batch correction?**

A: Optimize **after** batch correction:

```
Data → Batch correction → DE analysis → Import → Optimize
```

Batch effects can affect optimal thresholds.

**Q13: Can I use the same optimal thresholds for different datasets?**

A: Generally no. Optimal thresholds are dataset-specific and depend on:
- Sample size
- Biological variability
- Effect size
- Technical factors

Always optimize for each new dataset.

**Q14: What if optimization gives stricter thresholds than I expected?**

A: Stricter thresholds are often better! They mean:
- Lower false positive rate
- More reliable genes
- Better replication
- Higher confidence results

Don't be alarmed if optimal FDR is 0.01 instead of 0.05.

**Q15: Can I use Module 8 with limma-voom or Wilcoxon results?**

A: Yes! Module 8 works with DE results from:
- DESeq2
- edgeR
- limma-voom
- Wilcoxon
- Any tool that provides log2FC, p-value, and FDR

Just import with Module 6 first, then optimize.

---

### Results Interpretation

**Q16: My optimized LFC threshold is 0. Is that okay?**

A: Yes, LFC = 0 means:
- Statistical significance is sufficient
- All fold changes included
- Optimization found no benefit to fold-change filtering

This is common when:
- Small fold changes are biologically relevant
- Sample size is large (good statistical power)
- You're in discovery mode

**Q17: Optimization found very few genes (N=50). Is this correct?**

A: Could be correct if:
- Weak biological signal
- Small sample size
- Very strict optimization for robustness

Check:
1. Do top genes make biological sense?
2. Is low gene count expected for your system?
3. Try different method to verify

If consistently low, may reflect true biology.

**Q18: Should I report both default and optimized results?**

A: For publication, report:
- **Main results**: Optimized parameters
- **Supplementary**: Comparison with defaults
- **Methods**: How optimization was performed

Example:
```
"We identified 1,250 differentially expressed genes using 
optimized thresholds (FDR < 0.0123, |log2FC| > 0.50) compared 
to 1,850 genes using default thresholds (FDR < 0.05)."
```

**Q19: What if Ground Truth and FDR Control give different results?**

A: Different methods may give different optima. This is normal:
- Ground Truth optimizes for known genes
- FDR Control optimizes for statistical FDR

**Choose based on**:
- Study goals (validation vs discovery)
- Data availability
- Which metric is more important

Or report both and focus on shared genes.

**Q20: How do I know if optimization actually improved my results?**

A: Compare:

1. **Overlap with validated genes** (if available):
   ```
   Default: 15/25 validated genes found (60%)
   Optimized: 21/25 validated genes found (84%)
   → 24% improvement!
   ```

2. **Pathway enrichment**:
   ```
   Default: 5 significant pathways (p < 0.05)
   Optimized: 12 significant pathways (p < 0.01)
   → Better biological signal
   ```

3. **Reproducibility** (if validation cohort):
   ```
   Default: Jaccard = 0.45
   Optimized: Jaccard = 0.68
   → 51% improvement in replication
   ```

---

## Real-World Examples

### Example 1: Breast Cancer Study

**Study Design**:
- 60 tumor samples (30 ER+, 30 ER-)
- Goal: Find biomarkers for ER status
- Have 35 known ER-associated genes from literature

**Approach**:
```bash
# Step 1: Import DE results
raptor import-de \
    -i deseq2_er_positive_vs_negative.csv \
    -m deseq2 \
    -o results/de_er

# Step 2: Optimize with ground truth
raptor optimize \
    -d results/de_er/de_result.pkl \
    -m ground-truth \
    -g known_er_genes.csv \
    -o results/opt_er
```

**Results**:
```yaml
best_parameters:
  alpha: 0.0145
  lfc_threshold: 0.75
best_score: 0.87  # F1-score
n_deg_genes: 980
performance_metrics:
  precision: 0.88
  recall: 0.86
  f1_score: 0.87
  true_positives: 30  # Found 30/35 known ER genes
```

**Interpretation**:
- Found 30 of 35 known ER genes (86% recall)
- 88% of 980 called genes likely true (precision)
- More stringent than default (0.0145 vs 0.05)
- Requires >1.7-fold change (LFC > 0.75)

**Next Steps**:
- Validate top 20 novel genes with qRT-PCR
- Pathway enrichment on 980 genes
- Use for classifier development

---

### Example 2: Drug Response Pilot

**Study Design**:
- 16 cell lines (8 sensitive, 8 resistant to drug X)
- Small pilot study
- No validated genes available
- Need robust results for follow-up

**Approach**:
```bash
# Step 1: Import DE results
raptor import-de \
    -i deseq2_sensitive_vs_resistant.csv \
    -m deseq2 \
    -o results/de_drug

# Step 2: Optimize for stability (no validated genes)
raptor optimize \
    -d results/de_drug/de_result.pkl \
    -m stability \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    --n-bootstrap 100 \
    -o results/opt_drug
```

**Results**:
```yaml
best_parameters:
  alpha: 0.0089
  lfc_threshold: 1.0
best_score: 0.91  # Stability score
n_deg_genes: 245
performance_metrics:
  stability_score: 0.91
  high_stability_genes: 223  # 91% selected >90% of time
```

**Interpretation**:
- Very strict thresholds (FDR < 0.009, LFC > 1.0)
- 245 genes with high stability (91%)
- 223 genes selected in >90% of bootstrap samples
- 2-fold change minimum ensures biological relevance

**Next Steps**:
- Focus on 223 high-stability genes
- Plan validation experiments
- Apply same thresholds to larger cohort

---

### Example 3: Multi-Center Clinical Trial

**Study Design**:
- Discovery: 80 patients (Site A)
- Validation: 75 patients (Site B)
- Same disease, same treatment
- Need parameters that replicate across sites

**Approach**:
```bash
# Import both cohorts
raptor import-de -i site_a_deseq2.csv -m deseq2 -o results/site_a
raptor import-de -i site_b_deseq2.csv -m deseq2 -o results/site_b

# Optimize for reproducibility
raptor optimize \
    -d results/site_a/de_result.pkl \
    -m reproducibility \
    --cohort2 results/site_b/de_result.pkl \
    -o results/opt_multicenter
```

**Results**:
```yaml
best_parameters:
  alpha: 0.0234
  lfc_threshold: 0.50
best_score: 0.68  # Jaccard index
n_deg_genes: 1450 (Site A), 1520 (Site B)
performance_metrics:
  jaccard_index: 0.68
  overlap_count: 987
  site_a_only: 463
  site_b_only: 533
```

**Interpretation**:
- Good cross-site replication (Jaccard = 0.68)
- 987 genes replicate across both sites
- Site-specific genes may reflect population or technical differences
- Moderate thresholds balance sensitivity and specificity

**Next Steps**:
- Focus on 987 shared genes for publication
- Investigate site-specific genes (batch effects?)
- Use shared genes for clinical predictor

---

### Example 4: Time Course Developmental Study

**Study Design**:
- 5 time points (0h, 6h, 12h, 24h, 48h)
- 4 replicates per time point
- Comparing each time point to baseline (0h)
- No validated genes

**Approach**:
```bash
# For each comparison, optimize separately
for TIME in 6h 12h 24h 48h; do
    raptor optimize \
        -d results/de_${TIME}_vs_0h/de_result.pkl \
        -m fdr-control \
        --fdr-target 0.01 \
        -o results/opt_${TIME}_vs_0h
done
```

**Results**:
| Time | Optimal Alpha | Optimal LFC | DEGs |
|------|--------------|-------------|------|
| 6h   | 0.0089 | 0.25 | 145 |
| 12h  | 0.0078 | 0.50 | 892 |
| 24h  | 0.0123 | 0.75 | 2150 |
| 48h  | 0.0156 | 1.00 | 3420 |

**Interpretation**:
- Progressive increase in DEGs over time
- Stricter FDR for later time points (more signal)
- Increasing LFC threshold (larger changes)
- Clear developmental progression

**Next Steps**:
- Cluster genes by expression pattern
- Identify early vs late response genes
- Pathway analysis for each time point

---

## Scientific Background

### Why Parameter Optimization Matters

**Traditional approach**: Use arbitrary thresholds
```
FDR < 0.05 + |log2FC| > 0
```

**Problems**:
1. Not validated for your specific data
2. Ignores sample size
3. Doesn't account for biological variability
4. May be too lenient or strict

**Optimization approach**: Find thresholds based on data
```
Search parameter space → Maximize objective → Validated thresholds
```

**Benefits**:
1. Data-driven, not arbitrary
2. Accounts for your specific study
3. Maximizes desired metric (replication, FDR control, etc.)
4. Scientifically defensible

---

### The Four Methods Explained

#### Method 1: Ground Truth Optimization

**Scientific Basis**: Supervised machine learning optimization

**Key Papers**:
- Binary classification evaluation (Powers, 2011)
- F1-score optimization (Lipton et al., 2014)

**How it works**:
```
1. True Positives (TP):  Known DE genes that pass thresholds
2. False Positives (FP): Other genes that pass thresholds
3. False Negatives (FN): Known DE genes that fail thresholds
4. True Negatives (TN):  Other genes that fail thresholds

Precision = TP / (TP + FP)  # How many called are true?
Recall = TP / (TP + FN)     # How many true are called?
F1 = 2 × (Precision × Recall) / (Precision + Recall)
```

**Optimize**: Find alpha and LFC that maximize F1-score

---

#### Method 2: FDR Control Optimization

**Scientific Basis**: False Discovery Rate control theory

**Key Papers**:
- Benjamini & Hochberg (1995): FDR control
- Storey & Tibshirani (2003): π₀ estimation
- Storey (2002): qvalue methodology

**How it works**:

1. **Estimate π₀** (proportion of true nulls):
   ```
   For lambda in [0.05, 0.1, ..., 0.95]:
       pi0(lambda) = #{p > lambda} / (n × (1-lambda))
   
   Smooth curve and extrapolate to lambda=1
   ```

2. **Estimate FDR** for threshold alpha:
   ```
   Estimated FDR = (pi0 × n × alpha) / #{p < alpha}
   ```

3. **Find alpha** closest to target FDR

**Advantages**:
- Statistically principled
- No validation data needed
- Well-established theory

---

#### Method 3: Stability Optimization

**Scientific Basis**: Bootstrap stability selection

**Key Papers**:
- Meinshausen & Bühlmann (2010): Stability selection
- Efron & Tibshirani (1994): Bootstrap methods

**How it works**:

```
For b = 1 to B (100 bootstrap samples):
    1. Resample samples with replacement
    2. Rerun DE analysis
    3. Apply candidate thresholds
    4. Record selected genes
    
For each gene:
    Selection_probability = #{times selected} / B
    
Stability_score = Average selection probability
```

**Optimize**: Find thresholds that maximize stability

**Why it works**: Genes selected consistently across samples are robust to sampling variability.

---

#### Method 4: Reproducibility Optimization

**Scientific Basis**: Replication science

**Key Papers**:
- Ioannidis (2005): Why most research findings are false
- Button et al. (2013): Power failure in neuroscience
- Collaboration (2015): Reproducibility project

**How it works**:

```
For candidate thresholds:
    1. Apply to Discovery cohort → Gene set A
    2. Apply to Validation cohort → Gene set B
    3. Calculate overlap:
       
       Jaccard = |A ∩ B| / |A ∪ B|
```

**Optimize**: Find thresholds that maximize Jaccard index

**Why it works**: Parameters that maximize overlap between independent cohorts are more likely to generalize.

---

### Search Strategies Explained

#### Grid Search

**Algorithm**:
```
For alpha in [0.001, 0.05, 0.10, 0.15, 0.20]:
    For lfc in [0.0, 0.5, 1.0, 1.5, 2.0]:
        Score = Evaluate(alpha, lfc)
        If Score > BestScore:
            BestScore = Score
            BestParams = (alpha, lfc)
```

**Complexity**: O(n × m) where n, m are grid points

**Guarantees**: Finds global optimum within grid

---

#### Random Search

**Algorithm** (Bergstra & Bengio, 2012):
```
For i = 1 to N iterations:
    alpha = Random(0.001, 0.20)
    lfc = Random(0.0, 2.0)
    Score = Evaluate(alpha, lfc)
    If Score > BestScore:
        BestScore = Score
        BestParams = (alpha, lfc)
```

**Complexity**: O(N)

**Why it works**: Often finds good solutions faster than grid search, especially in high dimensions

---

#### Differential Evolution

**Algorithm** (Storn & Price, 1997):
```
Initialize population P of size NP
For generation = 1 to MaxGen:
    For each x in P:
        # Mutation
        a, b, c = Random 3 members from P
        v = a + F × (b - c)
        
        # Crossover
        u = Crossover(x, v)
        
        # Selection
        If Score(u) > Score(x):
            P = P - {x} + {u}
```

**Complexity**: O(NP × MaxGen)

**Advantages**: Good at finding global optimum, adaptive

---

### Statistical Guarantees

#### FDR Control

**Benjamini-Hochberg procedure** guarantees:
```
FDR ≤ alpha
```

Under independence or positive dependence.

**Storey's q-value** refines this:
```
qvalue(p) = min FDR over all thresholds ≥ p
```

More powerful while maintaining FDR control.

---

#### Stability Selection

**Meinshausen & Bühlmann (2010)** show:
```
E[V] ≤ (1/(2π_thr - 1)) × q² × p / E[S]
```

Where:
- V = false discoveries
- π_thr = selection threshold
- q = threshold parameter
- p = total variables
- E[S] = expected selections

**Practical**: High stability (>0.8) → Low expected false positives

---

### Computational Complexity

| Method | Complexity | Typical Time |
|--------|-----------|--------------|
| Ground Truth | O(iter × n_genes) | ~5 min |
| FDR Control | O(iter × n_genes) | ~3 min |
| Stability | O(B × iter × DE_time) | ~20 min |
| Reproducibility | O(iter × n_genes) | ~5 min |

Where:
- iter = optimization iterations
- n_genes = number of genes
- B = bootstrap samples
- DE_time = time to run DE analysis

---

## Conclusion

Module 8 provides **scientifically validated parameter optimization** for differential expression analysis:

✅ **4 Methods**: Choose based on available data
✅ **3 Search Strategies**: Balance speed and optimality
✅ **Automated**: Easy-to-use CLI and Python API
✅ **Publication-Ready**: Scientifically justified thresholds
✅ **Reproducible**: Documented, versioned results

### Quick Takeaways

1. **Always optimize** - Don't rely on arbitrary defaults
2. **Choose method carefully** - Use decision tree in this guide
3. **Validate results** - Check that optimal parameters make sense
4. **Document thoroughly** - For reproducibility and publication
5. **Apply consistently** - Use same thresholds across comparisons

### Next Steps

After optimization:
1. ✅ Review optimal parameters
2. ✅ Use filtered gene list (deg_genes.csv)
3. ✅ Continue to Module 9 (Ensemble Analysis)
4. ✅ Or proceed to downstream analysis (pathways, validation)

---

**Thank you for using RAPTOR Module 8!**

For questions, issues, or feedback:
- **GitHub**: https://github.com/AyehBlk/RAPTOR
- **Email**: ayehbolouki1988@gmail.com
- **Documentation**: Full API docs at RAPTOR website

**Happy optimizing!** 🦖📊

---

**Version**: 2.2.0  
**Last Updated**: March 9, 2026  
**Author**: Ayeh Bolouki  
**License**: MIT
