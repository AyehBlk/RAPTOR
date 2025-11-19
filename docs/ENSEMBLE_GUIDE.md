# RAPTOR v2.1.0 Ensemble Analysis Guide

**Combine multiple pipelines for robust differential expression detection**

---

## 🎯 What is Ensemble Analysis?

Instead of trusting a single pipeline, **ensemble analysis combines results from multiple pipelines** to identify genes that are consistently called as differentially expressed across different methods.

### The Problem

Different RNA-seq pipelines can give different results:

```
Pipeline 1: Found 1,245 DE genes
Pipeline 3: Found 1,089 DE genes
Pipeline 4: Found 982 DE genes

Overlap: Only 734 genes called by all three (59%)
```

**Which genes should you trust?**

### The Solution

Ensemble analysis provides:
- ✅ **Higher confidence** - Genes must be called by multiple methods
- ✅ **Reduced false positives** - Method-specific artifacts eliminated
- ✅ **Robust biomarkers** - Genes that survive multiple analyses
- ✅ **Evidence-based decisions** - See agreement across methods

---

## 🚀 Quick Start

### Basic Ensemble Analysis

```bash
# Run ensemble with 3 recommended pipelines
raptor ensemble \
  --counts data.csv \
  --metadata metadata.csv \
  --pipelines 1,3,5 \
  --output ensemble_results/

# Or let RAPTOR choose pipelines automatically
raptor ensemble \
  --counts data.csv \
  --metadata metadata.csv \
  --auto-select \
  --output ensemble_results/
```

### Python API

```python
from raptor import EnsembleAnalyzer, PipelineRunner

# Run multiple pipelines
pipelines = [1, 3, 5]
results = {}
for p in pipelines:
    runner = PipelineRunner(pipeline_id=p)
    results[p] = runner.run(counts, metadata)

# Combine with ensemble
ensemble = EnsembleAnalyzer(method='meta_analysis')
consensus = ensemble.combine(results)

# View high-confidence genes
high_conf = consensus[consensus['confidence'] > 0.8]
print(f"High confidence genes: {len(high_conf)}")
```

---

## 🎓 When to Use Ensemble Analysis

### ✅ Use Ensemble When:

1. **High-stakes decisions**
   - Clinical biomarker discovery
   - Drug target validation
   - Regulatory submissions

2. **Conflicting results**
   - Different methods give different answers
   - Need to resolve discordance
   - Validating previous findings

3. **Publication requirements**
   - Reviewers request robustness checks
   - Multi-method validation needed
   - High-impact journals

4. **Maximum confidence needed**
   - Follow-up experiments are expensive
   - Limited validation resources
   - Irreversible clinical decisions

### ❌ Don't Use Ensemble When:

- Exploratory analysis (use single fast pipeline)
- Limited computational resources
- Time-critical decisions
- Data quality is poor (fix data first)

---

## 🔬 Ensemble Methods Available

RAPTOR v2.1.0 provides **6 ensemble methods**:

### 1. Rank Aggregation

**How it works:** Combine gene rankings from multiple pipelines

```python
ensemble = EnsembleAnalyzer(method='rank_aggregation')
consensus = ensemble.combine(results)
```

**Best for:**
- When absolute p-values differ but rankings agree
- Comparing different statistical methods
- Prioritizing genes for follow-up

**Methods:**
- **Borda count:** Simple ranking sum
- **Robust Rank Aggregation (RRA):** Statistical combination

**Example:**
```
Pipeline 1 ranking: Gene1(1st), Gene2(2nd), Gene3(5th)
Pipeline 3 ranking: Gene1(2nd), Gene2(1st), Gene3(4th)
Pipeline 5 ranking: Gene1(1st), Gene2(3rd), Gene3(6th)

Consensus ranking: Gene1(1st), Gene2(2nd), Gene3(5th)
```

### 2. Vote Counting

**How it works:** Simple majority vote - gene is DE if N pipelines call it

```python
ensemble = EnsembleAnalyzer(
    method='vote_counting',
    vote_threshold=0.6  # 60% must agree
)
consensus = ensemble.combine(results)
```

**Best for:**
- Simple, interpretable results
- When pipelines have similar performance
- Conservative gene lists

**Example:**
```
5 pipelines tested:
- Gene1: Called DE by 5/5 pipelines → Include (100% agreement)
- Gene2: Called DE by 4/5 pipelines → Include (80% agreement)
- Gene3: Called DE by 3/5 pipelines → Include (60% agreement)
- Gene4: Called DE by 2/5 pipelines → Exclude (40% agreement)
```

**Configuration:**
```yaml
ensemble:
  vote_counting:
    vote_threshold: 0.6  # Require 60% agreement
    weighted: true  # Weight by pipeline performance
```

### 3. Weighted Average

**How it works:** Combine p-values weighted by pipeline accuracy

```python
ensemble = EnsembleAnalyzer(
    method='weighted_average',
    weights={1: 1.2, 3: 1.0, 5: 0.8}  # Based on benchmark F1 scores
)
consensus = ensemble.combine(results)
```

**Best for:**
- When pipelines have different reliabilities
- You have benchmark data
- Incorporating prior knowledge

**Example:**
```
Pipeline weights (from benchmarks):
  P1: 1.2 (highest accuracy)
  P3: 1.0 (baseline)
  P5: 0.8 (lower for your data type)

Gene1 p-values:
  P1: 0.001 (weight 1.2) → effective: 0.0008
  P3: 0.005 (weight 1.0) → effective: 0.005
  P5: 0.01  (weight 0.8) → effective: 0.0125

Weighted combined p-value: 0.003
```

### 4. Meta-Analysis (Fisher's Method) ⭐ RECOMMENDED

**How it works:** Statistically combine p-values using Fisher's method

```python
ensemble = EnsembleAnalyzer(method='meta_analysis')
consensus = ensemble.combine(results)
```

**Best for:**
- Most rigorous statistical combination
- Publication-quality analysis
- When p-values are reliable

**How it works:**
```
Fisher's combined statistic: χ² = -2 Σ ln(p_i)

Pipeline p-values for Gene1:
  P1: 0.001
  P3: 0.005
  P5: 0.002

χ² = -2[ln(0.001) + ln(0.005) + ln(0.002)] = 29.4
Combined p-value = 4.2e-6 (highly significant)
```

**Advantages:**
- ✅ Well-established statistical method
- ✅ Handles different sample sizes
- ✅ Proper p-value combination
- ✅ Conservative (reduces false positives)

### 5. Intersection (Most Conservative)

**How it works:** Keep only genes called by ALL pipelines

```python
ensemble = EnsembleAnalyzer(method='intersection')
consensus = ensemble.combine(results)
```

**Best for:**
- Highest confidence genes only
- Expensive follow-up experiments
- Clinical applications

**Example:**
```
Pipeline 1: 1,245 genes
Pipeline 3: 1,089 genes  
Pipeline 5: 1,156 genes

Intersection: 734 genes (called by all three)
```

**Trade-off:**
- ✅ Very low false positive rate
- ❌ May miss true positives
- ❌ Small gene lists

### 6. Union (Most Liberal)

**How it works:** Include gene if ANY pipeline calls it

```python
ensemble = EnsembleAnalyzer(method='union')
consensus = ensemble.combine(results)
```

**Best for:**
- Exploratory screening
- Don't want to miss candidates
- Hypothesis generation

**Example:**
```
Pipeline 1: 1,245 genes
Pipeline 3: 1,089 genes
Pipeline 5: 1,156 genes

Union: 1,834 genes (called by at least one)
```

**Trade-off:**
- ✅ Catches all potential candidates
- ❌ Higher false positive rate
- ❌ Large gene lists need filtering

---

## 📊 Understanding Ensemble Results

### Confidence Scores

RAPTOR assigns confidence scores (0-1) to each gene:

```python
High Confidence (>0.8):    Called by 4-5 pipelines, low p-values
Medium Confidence (0.6-0.8): Called by 3 pipelines
Low Confidence (<0.6):      Called by 1-2 pipelines, borderline
```

**Interpreting confidence:**

| Score | Meaning | Recommendation |
|-------|---------|----------------|
| **0.9-1.0** | All methods agree strongly | Trust for clinical use |
| **0.8-0.9** | Strong consensus | Good biomarker candidates |
| **0.6-0.8** | Moderate agreement | Validate experimentally |
| **0.4-0.6** | Weak consensus | Further investigation |
| **<0.4** | High discordance | Be cautious |

### Example Output

```
Gene_ID          Pipelines  Confidence  Mean_LogFC  Combined_P  Status
─────────────────────────────────────────────────────────────────────────
ENSG00000000003  5/5        0.95        2.34        1.2e-12     ✓ High
ENSG00000000005  4/5        0.82        1.87        3.4e-08     ✓ High  
ENSG00000000419  3/5        0.67        1.45        2.1e-05     ⚠ Medium
ENSG00000000457  2/5        0.48        1.12        0.002       ⚠ Low
ENSG00000000460  5/5        0.93        -2.12       5.6e-11     ✓ High
```

### Visualizations Generated

1. **Venn Diagram** (2-3 pipelines)
```
     Pipeline 1
        ╱     ╲
     512   234   378
      │    ╱ ╲    │
      │  734  │   │
      │   │   │   │
    Pipeline 3  Pipeline 5
```

2. **UpSet Plot** (4+ pipelines)
```
Intersection Size
    ↑
500 │           █
400 │       █   █
300 │   █   █   █
200 │   █   █   █   █
100 │   █   █   █   █   █
    └───────────────────────→
        P1  P3  P5  P1  All
        only    only P3  five
```

3. **Concordance Heatmap**
```
Pipeline Agreement Matrix:
     P1   P3   P5   P7
P1  1.00 0.87 0.82 0.79
P3  0.87 1.00 0.91 0.84
P5  0.82 0.91 1.00 0.88
P7  0.79 0.84 0.88 1.00
```

4. **Confidence Distribution**
```
Number of Genes
    ↑
200 │     █
150 │     █
100 │     █  █
 50 │  █  █  █  █
  0 │  █  █  █  █  █
    └─────────────────→
      0.2 0.4 0.6 0.8 1.0
        Confidence Score
```

---

## 🔧 Advanced Configuration

### Custom Ensemble Configuration

Create `config/ensemble_config.yaml`:

```yaml
ensemble:
  enabled: true
  
  # Pipeline selection
  auto_select: true
  min_pipelines: 3
  max_pipelines: 5
  
  # Ensure diversity
  ensure_diversity: true
  diversity_criteria:
    - "alignment_method"     # STAR vs Salmon vs kallisto
    - "quantification_method"
    - "de_method"            # DESeq2 vs edgeR vs limma
  
  # Primary ensemble method
  methods:
    - "meta_analysis"        # Use Fisher's method
    - "vote_counting"        # Also compute votes
    - "rank_aggregation"     # And rankings
  
  primary_method: "meta_analysis"
  
  # Meta-analysis settings
  meta_analysis:
    method: "fishers"
    p_value_threshold: 0.05
    adjust_method: "BH"
    heterogeneity_test: true
  
  # Vote counting settings
  vote_counting:
    vote_threshold: 0.6
    weighted: true
    confidence_levels:
      high: 0.8
      medium: 0.6
      low: 0.4
  
  # Consensus requirements
  consensus:
    min_agreement: 0.6       # 60% of pipelines
    require_direction_consensus: true  # Same up/down regulation
    direction_agreement_threshold: 0.8
  
  # Output
  generate_ensemble_report: true
  save_individual_results: true
  plot_comparisons: true
```

### Running with Custom Config

```bash
raptor ensemble \
  --counts data.csv \
  --metadata metadata.csv \
  --config config/ensemble_config.yaml \
  --output ensemble_results/
```

---

## 💡 Best Practices

### Pipeline Selection

**Good combinations (diversity):**
```
✓ Pipelines 1, 3, 5  # STAR+DESeq2, Salmon+edgeR, STAR+limma
✓ Pipelines 1, 4, 7  # STAR+DESeq2, kallisto+sleuth, Bowtie2+EBSeq
```

**Poor combinations (too similar):**
```
✗ Pipelines 1, 5, 6  # All use STAR alignment
✗ Pipelines 3, 4     # Both pseudo-alignment only
```

**Recommended starter sets:**

| Data Type | Recommended Pipelines | Rationale |
|-----------|----------------------|-----------|
| Standard (n=6-12) | 1, 3, 5 | Balance of methods |
| Large (n>20) | 3, 4, 5 | Fast, robust methods |
| Small (n<6) | 1, 5, 6 | Methods for small n |
| Clinical | 1, 3, 5, 7 | Maximum validation |

### Confidence Thresholds

**For different applications:**

```yaml
# Discovery phase (exploratory)
confidence_threshold: 0.5

# Validation phase
confidence_threshold: 0.7

# Clinical/biomarker
confidence_threshold: 0.9
```

### Handling Discordance

When pipelines disagree:

1. **Check data quality**
   ```bash
   raptor qc --counts data.csv
   ```

2. **Examine discordant genes**
   ```python
   # Genes with high variance across pipelines
   discordant = consensus[consensus['agreement_cv'] > 0.5]
   ```

3. **Investigate reasons**
   - Different normalization methods
   - Outlier samples affecting some pipelines
   - Method-specific biases

4. **Use more pipelines**
   - If 3 pipelines disagree, try 5-6

---

## 📈 Real-World Examples

### Example 1: Cancer Biomarker Discovery

**Scenario:**
- 24 tumor samples vs 24 normal
- Need highly confident biomarkers
- Expensive validation (qPCR)

**Approach:**
```bash
# Run 5 diverse pipelines
raptor ensemble \
  --counts data.csv \
  --metadata metadata.csv \
  --pipelines 1,3,4,5,7 \
  --method intersection \
  --output biomarker_discovery/

# Filter for very high confidence
# Only genes called by all 5 pipelines
```

**Results:**
```
Total genes tested: 15,234
Pipeline 1: 1,456 DE genes
Pipeline 3: 1,289 DE genes
Pipeline 4: 1,123 DE genes
Pipeline 5: 1,367 DE genes
Pipeline 7: 1,198 DE genes

Intersection: 678 genes (called by all 5)
High confidence (>0.9): 456 genes

→ Select top 20 by effect size for validation
```

### Example 2: Resolving Conflicting Literature

**Scenario:**
- Previous study used Pipeline 1, found Gene X significant
- Your replication with Pipeline 3: Gene X not significant
- Need to resolve discrepancy

**Approach:**
```bash
# Run both pipelines plus neutral third party
raptor ensemble \
  --counts data.csv \
  --metadata metadata.csv \
  --pipelines 1,3,5 \
  --method meta_analysis \
  --output conflict_resolution/
```

**Results:**
```
Gene X results:
  Pipeline 1: p=0.02, logFC=1.5  ✓ Significant
  Pipeline 3: p=0.08, logFC=1.2  ✗ Not significant
  Pipeline 5: p=0.04, logFC=1.3  ✓ Significant

Meta-analysis: p=0.03, combined evidence suggests significance
Confidence: 0.67 (Medium)

Conclusion: Gene X is likely DE, but borderline
Recommendation: Validate experimentally
```

### Example 3: Large Clinical Study

**Scenario:**
- 96 patient samples
- Multi-center study (batch effects)
- Regulatory submission planned

**Approach:**
```bash
# Use methods good for large samples + batch correction
raptor ensemble \
  --counts data.csv \
  --metadata metadata.csv \
  --pipelines 1,3,5 \
  --method meta_analysis \
  --batch-column center \
  --output clinical_study/
```

**Results:**
```
Robust biomarkers identified:
- High confidence (>0.8): 234 genes
- Direction concordance: 98%
- Low heterogeneity (I²<25%): 89%

Ready for:
✓ Regulatory submission
✓ Independent validation cohort
✓ Clinical assay development
```

---

## 🐛 Troubleshooting

### Issue: Low agreement between pipelines

**Possible causes:**
- Poor data quality
- Batch effects
- Outlier samples
- Wrong pipeline choices

**Solutions:**
```bash
# 1. Check data quality
raptor qc --counts data.csv --detailed

# 2. Remove outliers
raptor ensemble \
  --counts data.csv \
  --exclude-samples outlier1,outlier2

# 3. Try more pipelines
raptor ensemble \
  --counts data.csv \
  --auto-select \
  --min-pipelines 5  # Use 5 instead of 3
```

### Issue: Ensemble takes too long

**Speed optimization:**
```bash
# Run pipelines in parallel
raptor ensemble \
  --counts data.csv \
  --pipelines 1,3,5 \
  --parallel \
  --max-jobs 3

# Use faster pipelines
raptor ensemble \
  --counts data.csv \
  --pipelines 3,4,5  # All fast methods

# Use pre-computed results
raptor ensemble \
  --results pipeline1_results/,pipeline3_results/,pipeline5_results/ \
  --combine-only
```

### Issue: Too few consensus genes

**If intersection is too small:**

```bash
# Lower consensus threshold
raptor ensemble \
  --counts data.csv \
  --pipelines 1,3,5 \
  --method vote_counting \
  --vote-threshold 0.5  # 50% instead of 60%

# Or use union then filter
raptor ensemble \
  --counts data.csv \
  --pipelines 1,3,5 \
  --method union \
  --post-filter confidence:0.7
```

---

## 📚 Statistical Background

### Why Combine P-values?

**Fisher's Method:**

$$\chi^2 = -2\sum_{i=1}^{k} \ln(p_i)$$

Follows chi-squared distribution with 2k degrees of freedom.

**Advantages:**
- Properly accounts for multiple testing
- More powerful than intersection
- Well-established theory

**Assumptions:**
- Tests are independent
- P-values are uniformly distributed under H₀
- Same null hypothesis tested

### Heterogeneity Assessment

**I² statistic:**
- I² = 0%: No heterogeneity (pipelines agree)
- I² = 25%: Low heterogeneity
- I² = 50%: Moderate heterogeneity
- I² = 75%: High heterogeneity (pipelines disagree)

**When I² is high:**
- Consider why pipelines disagree
- Use random effects model
- Report heterogeneity in results

---

## 📊 Reporting Ensemble Results

### For Publications

**Methods Section:**
```
Ensemble differential expression analysis was performed using 
RAPTOR v2.1.0 (Bolouki, 2025). Three independent pipelines 
were executed: (1) STAR-RSEM-DESeq2, (2) Salmon-edgeR, and 
(3) STAR-HTSeq-limma-voom. P-values were combined using 
Fisher's method with Benjamini-Hochberg adjustment (α=0.05). 
Genes were considered high-confidence if called by all three 
pipelines (intersection) with combined FDR<0.05 and |log2FC|>1.
```

**Results Section:**
```
Ensemble analysis identified 456 high-confidence differentially 
expressed genes (agreement score >0.8). Pipeline concordance was 
high (κ=0.82, I²=18%), indicating robust results. The intersection 
of all three pipelines yielded 234 genes, representing conservative 
biomarker candidates (Supplementary Table S1).
```

### Supplementary Materials

Include:
1. **Table S1:** All ensemble genes with confidence scores
2. **Table S2:** Individual pipeline results
3. **Figure S1:** Venn diagram or UpSet plot
4. **Figure S2:** Concordance heatmap
5. **Table S3:** Discordant genes and reasons

---

## 🎓 Summary

### Key Takeaways

✅ **When to use ensemble:**
- High-stakes decisions
- Need maximum confidence
- Publication requirements

✅ **Best methods:**
- Meta-analysis (Fisher's) for rigor
- Vote counting for simplicity
- Intersection for highest confidence

✅ **Pipeline selection:**
- 3-5 pipelines recommended
- Ensure diversity (alignment, stats)
- Include at least one slow/accurate method

✅ **Interpreting results:**
- High confidence (>0.8): Trust for action
- Medium confidence (0.6-0.8): Validate
- Low confidence (<0.6): Investigate

---

## 📧 Support

Questions about ensemble analysis?

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Discussions:** https://github.com/AyehBlk/RAPTOR/discussions
- **Email:** ayehbolouki1988@gmail.com

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**Version:** 2.1.0  
**License:** MIT

---

*"Not all pipelines are created equal, but their consensus might be."* 🎯
