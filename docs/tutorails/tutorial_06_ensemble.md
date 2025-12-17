# Tutorial 6: Ensemble Analysis - Combining Multiple Pipelines

**Level**: Advanced  
**Time**: 2-3 hours (mostly automated)  
**Goal**: Master combining results from multiple pipelines for maximum confidence

---

## What You'll Learn

- When and why to use ensemble analysis
- How to combine results from multiple pipelines
- How to identify high-confidence genes
- How to handle discordant results
- Best practices for robust analysis
- **NEW in v2.1.1**: Using ATO for uniform thresholds across pipelines

---

## 🆕 What's New in v2.1.1

**ATO + Ensemble Integration**: Apply the same data-driven threshold across all pipelines for fair comparison:

```bash
raptor ensemble --counts data.csv --metadata meta.txt \
                --pipelines STAR-DESeq2,Salmon-DESeq2,Kallisto-DESeq2 \
                --use-ato \
                --ato-goal balanced \
                --output ensemble_ato/
```

This ensures all pipelines use the **same statistically-justified threshold** rather than arbitrary |logFC| > 1. Results include auto-generated methods text for publications.

- Completed Tutorial 1 and Tutorial 2
- Understanding of differential expression analysis
- RAPTOR v2.1.1+ installed
- Sufficient computational resources to run multiple pipelines

---

## What is Ensemble Analysis?

### The Problem: Pipeline Uncertainty

```
Your Data → Pipeline 1 → 1,234 DE genes
Your Data → Pipeline 2 → 1,567 DE genes  
Your Data → Pipeline 3 → 1,089 DE genes

Which result do you trust? 
```

### The Solution: Ensemble Approach

```
Your Data → Pipeline 1 ┐
Your Data → Pipeline 2 ├→ COMBINE → High-confidence genes
Your Data → Pipeline 3 ┘
```

**Key insight:** Genes found by multiple pipelines are more likely to be true positives!

---

## Why Use Ensemble Analysis?

### Benefits

✅ **Higher confidence** - Multiple methods agree  
✅ **Reduced false positives** - Consensus filtering  
✅ **Method robustness** - Not dependent on single pipeline  
✅ **Better for publication** - Shows comprehensive validation  
✅ **Catches edge cases** - Different methods find different things  

### When to Use

**Use ensemble when:**
-  Critical research decisions
-  Publishing high-impact papers
-  Expensive follow-up experiments planned
-  Contradictory results from single pipeline
-  Need maximum confidence

**Skip ensemble when:**
-  Exploratory analysis only
-  Time pressure
-  Limited computational resources
-  Single pipeline gives clear results

---

## Step 1: Choose Pipelines for Ensemble

### Strategy 1: Top 3 from ML Recommendation

```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler
import pandas as pd

# Get ML recommendations
counts = pd.read_csv('data/counts.csv', index_col=0)
metadata = pd.read_csv('data/metadata.csv')

profiler = RNAseqDataProfiler(counts, metadata, use_ml=True)
profile = profiler.profile()

ml_rec = MLPipelineRecommender()
recommendations = ml_rec.recommend(profile, n=3)

# Use top 3 for ensemble
ensemble_pipelines = [r['pipeline_id'] for r in recommendations]
print(f"Ensemble will use pipelines: {ensemble_pipelines}")
```

**Example output:**
```
Ensemble will use pipelines: [3, 1, 4]
- Pipeline 3: Salmon-edgeR (fast, accurate)
- Pipeline 1: STAR-RSEM-DESeq2 (most accurate)
- Pipeline 4: Kallisto-Sleuth (very fast)
```

### Strategy 2: Methodological Diversity

Choose pipelines that use different approaches:

```python
# Diverse ensemble
ensemble_pipelines = [
    1,  # Alignment-based + DESeq2 (negative binomial)
    3,  # Pseudo-alignment + edgeR (quasi-likelihood)
    5   # Alignment-based + limma (voom transformation)
]
```

**Why this works:**
- Different alignment methods (STAR vs Salmon vs Kallisto)
- Different statistical approaches (DESeq2 vs edgeR vs limma)
- If all three agree → Very confident!

---

## Step 2: Run Multiple Pipelines

### Quick Ensemble (3 pipelines)

```bash
# Run top 3 recommended pipelines
raptor compare \
  --data fastq/ \
  --output ensemble_results/ \
  --pipelines 3 1 4 \
  --mode quick \
  --threads 16 \
  --memory 64G
```

**Time:** ~2-4 hours (vs 12+ hours for all 8)

### Full Ensemble (All 8 pipelines)

```bash
# For maximum confidence
raptor compare \
  --data fastq/ \
  --output full_ensemble/ \
  --pipelines all \
  --mode full \
  --threads 32 \
  --memory 128G
```

**Time:** 12-24 hours

### From Existing Results

Already have results from multiple pipelines?

```python
from raptor import EnsembleAnalyzer

# Load existing results
results = {
    'salmon_edgeR': 'results/pipeline3/diff_expr.csv',
    'star_deseq2': 'results/pipeline1/diff_expr.csv',
    'kallisto_sleuth': 'results/pipeline4/diff_expr.csv'
}

analyzer = EnsembleAnalyzer()
ensemble = analyzer.load_and_combine(results)
```

### Ensemble with ATO (NEW in v2.1.1)

Apply uniform data-driven thresholds for fair comparison:

```bash
# Command line - ensemble with ATO
raptor ensemble \
  --counts data.csv \
  --metadata meta.txt \
  --pipelines STAR-DESeq2,Salmon-DESeq2,Kallisto-DESeq2 \
  --use-ato \
  --ato-goal balanced \
  --output ensemble_ato/
```

```python
# Python API - ensemble with ATO
from raptor import EnsembleAnalyzer
from raptor.threshold_optimizer import optimize_thresholds

# First, combine all pipeline results
combined = pd.concat([
    pd.read_csv('results/star_deseq2.csv'),
    pd.read_csv('results/salmon_deseq2.csv'),
    pd.read_csv('results/kallisto_deseq2.csv')
])

# Get uniform threshold from ATO
ato_result = optimize_thresholds(combined, goal='balanced')
uniform_threshold = ato_result.logfc_threshold

print(f"Uniform threshold: {uniform_threshold:.3f}")
print(f"Methods text: {ato_result.methods_text}")

# Apply to ensemble analysis
analyzer = EnsembleAnalyzer()
ensemble = analyzer.combine(
    results_dir='results/',
    logfc_threshold=uniform_threshold,  # Use ATO threshold
    fdr_threshold=0.05
)
```

**Why ATO + Ensemble?**
- **Fair comparison**: All pipelines use the same statistically-justified threshold
- **No arbitrary cutoffs**: Data determines the optimal |logFC|
- **Publication-ready**: Auto-generated methods text explains your approach

### Basic Combination (Voting)

```python
from raptor import EnsembleAnalyzer

# Load results from all pipelines
analyzer = EnsembleAnalyzer()

# Combine using voting (default)
ensemble = analyzer.combine(
    results_dir='ensemble_results/',
    method='vote',  # Majority vote
    min_agreement=2  # At least 2 out of 3 must agree
)

print(f"Total unique genes across all pipelines: {ensemble['total_genes']}")
print(f"Consensus genes (2+ pipelines): {ensemble['consensus_count']}")
print(f"High confidence (all 3 pipelines): {ensemble['high_conf_count']}")
```

**Example output:**
```
Total unique genes across all pipelines: 2,345
Consensus genes (2+ pipelines): 1,234
High confidence (all 3 pipelines): 892
```

### Weighted Combination

```python
# Weight pipelines based on ML confidence
weights = {
    'salmon_edgeR': 0.4,    # ML score: 0.89
    'star_deseq2': 0.4,     # ML score: 0.85
    'kallisto_sleuth': 0.2  # ML score: 0.78
}

ensemble = analyzer.combine(
    results_dir='ensemble_results/',
    method='weighted',
    weights=weights
)
```

### Rank-Based Combination

```python
# Combine based on gene rankings
ensemble = analyzer.combine(
    results_dir='ensemble_results/',
    method='rank',  # Use rank aggregation
    top_n=2000  # Consider top 2000 from each pipeline
)
```

---

## Step 4: Analyze Consensus

### Consensus Genes (High Confidence)

```python
# Get genes found by multiple pipelines
consensus_genes = ensemble.get_consensus(min_agreement=2)

print(f"Found {len(consensus_genes)} consensus genes")
print("\nTop 10 consensus genes:")
print(consensus_genes.sort_values('confidence', ascending=False).head(10))
```

**Output:**
```
Found 1,234 consensus genes

Top 10 consensus genes:
gene_id        confidence  pipelines_found  avg_log2fc  avg_pvalue
ENSG00000001   1.00        3/3              2.45        1.2e-45
ENSG00000002   1.00        3/3             -1.89        3.4e-38
ENSG00000003   0.95        3/3              1.67        8.9e-32
ENSG00000004   0.92        3/3             -2.12        2.1e-28
...
```

### Confidence Scores

```python
# Calculate per-gene confidence
confidence_df = ensemble.calculate_confidence()

# Categorize by confidence
high_conf = confidence_df[confidence_df['confidence'] >= 0.9]
med_conf = confidence_df[(confidence_df['confidence'] >= 0.7) & 
                         (confidence_df['confidence'] < 0.9)]
low_conf = confidence_df[confidence_df['confidence'] < 0.7]

print(f"High confidence: {len(high_conf)} genes")
print(f"Medium confidence: {len(med_conf)} genes")
print(f"Low confidence: {len(low_conf)} genes")
```

---

## Step 5: Understand Discordance

### Identify Discordant Genes

```python
# Find genes where pipelines disagree
discordant = ensemble.identify_discordant()

print(f"Discordant genes: {len(discordant)}")
print("\nReasons for discordance:")
for reason, count in discordant.groupby('reason').size().items():
    print(f"  {reason}: {count} genes")
```

**Common reasons for discordance:**
```
Statistical method differences: 345 genes
Quantification differences: 234 genes
Fold change disagreement: 123 genes
Borderline significance: 89 genes
```

### Investigate Specific Discordant Gene

```python
# Check specific gene
gene_id = 'ENSG00001234'
gene_info = ensemble.get_gene_details(gene_id)

print(f"\nGene: {gene_id}")
print("\nResults across pipelines:")
for pipeline, result in gene_info['pipeline_results'].items():
    print(f"\n{pipeline}:")
    print(f"  Log2FC: {result['log2fc']:.2f}")
    print(f"  P-value: {result['pvalue']:.2e}")
    print(f"  FDR: {result['fdr']:.2e}")
    print(f"  Called DE: {result['is_de']}")

print(f"\nConsensus: {gene_info['consensus']}")
print(f"Confidence: {gene_info['confidence']:.2f}")
```

---

## Step 6: Visualization

### Venn Diagram

```python
# Create Venn diagram
ensemble.generate_venn_diagram(
    output='ensemble_venn.png',
    title='DE Genes Overlap'
)
```

**Shows:**
```
       Pipeline 1 (1,234)
            /  \
           /    \
    456  /  892  \  321
        /________\
   Pipeline 2    Pipeline 3
   (1,567)       (1,089)
   
Overlap (all 3): 892 genes
```

### UpSet Plot (Better for >3 pipelines)

```python
# For 4+ pipelines, use UpSet plot
ensemble.generate_upset_plot(
    output='ensemble_upset.png',
    title='Pipeline Intersections'
)
```

### Concordance Matrix

```python
# How similar are pipelines?
ensemble.plot_concordance_matrix(
    output='concordance.png'
)
```

**Output:**
```
           P1    P2    P3    P4
Pipeline1 1.00  0.87  0.82  0.79
Pipeline2 0.87  1.00  0.84  0.81
Pipeline3 0.82  0.84  1.00  0.76
Pipeline4 0.79  0.81  0.76  1.00
```

### Confidence Distribution

```python
import matplotlib.pyplot as plt

confidence_df = ensemble.calculate_confidence()

plt.figure(figsize=(10, 6))
plt.hist(confidence_df['confidence'], bins=50, edgecolor='black')
plt.xlabel('Confidence Score')
plt.ylabel('Number of Genes')
plt.title('Distribution of Gene Confidence Scores')
plt.axvline(0.9, color='red', linestyle='--', label='High confidence cutoff')
plt.legend()
plt.savefig('confidence_distribution.png', dpi=300)
```

---

## Step 7: Filter and Export

### Export High-Confidence Genes

```python
# Get high-confidence genes only
high_conf_genes = ensemble.get_consensus(
    min_agreement=3,  # All pipelines agree
    min_confidence=0.9
)

# Export for downstream analysis
high_conf_genes.to_csv('high_confidence_genes.csv', index=False)

print(f"Exported {len(high_conf_genes)} high-confidence genes")
```

### Create Multiple Gene Lists

```python
# Different confidence levels for different uses
lists = {
    'very_high': ensemble.get_consensus(min_agreement=3, min_confidence=0.95),
    'high': ensemble.get_consensus(min_agreement=3, min_confidence=0.85),
    'medium': ensemble.get_consensus(min_agreement=2, min_confidence=0.75),
}

for category, genes in lists.items():
    filename = f'genes_{category}_confidence.csv'
    genes.to_csv(filename, index=False)
    print(f"{category.upper()}: {len(genes)} genes → {filename}")
```

**Output:**
```
VERY_HIGH: 687 genes → genes_very_high_confidence.csv
HIGH: 892 genes → genes_high_confidence.csv
MEDIUM: 1,234 genes → genes_medium_confidence.csv
```

---

## Step 8: Generate Ensemble Report

```python
from raptor import ReportGenerator

report_gen = ReportGenerator()

# Comprehensive ensemble report
report_gen.generate_ensemble_report(
    ensemble,
    output='ensemble_report.html',
    include_venn=True,
    include_concordance=True,
    include_gene_table=True
)

print("Ensemble report generated: ensemble_report.html")
```

**Report includes:**
- Summary statistics
- Venn diagram
- Concordance matrix
- High-confidence gene table
- Discordant genes analysis
- Parameter comparison across pipelines
- Methods section for manuscript

---

## Step 9: Validate Ensemble Results

### Compare with Single Pipeline

```python
# How many extra genes did ensemble find?
single_pipeline = pd.read_csv('results/pipeline3/diff_expr.csv')
single_de = single_pipeline[single_pipeline['FDR'] < 0.05]

ensemble_consensus = ensemble.get_consensus(min_agreement=2)

# Genes found by ensemble but not single pipeline
ensemble_only = set(ensemble_consensus['gene_id']) - set(single_de['gene_id'])
both = set(ensemble_consensus['gene_id']) & set(single_de['gene_id'])

print(f"Single pipeline: {len(single_de)} genes")
print(f"Ensemble consensus: {len(ensemble_consensus)} genes")
print(f"Found by both: {len(both)} genes")
print(f"Ensemble added: {len(ensemble_only)} genes")
print(f"Ensemble removed: {len(set(single_de['gene_id']) - set(ensemble_consensus['gene_id']))} genes")
```

### Cross-Validation

```python
# Split data and test ensemble stability
from raptor import EnsembleAnalyzer

analyzer = EnsembleAnalyzer()
stability = analyzer.cross_validate_ensemble(
    results_dir='ensemble_results/',
    n_splits=5,
    method='vote'
)

print(f"Ensemble stability: {stability['mean_jaccard']:.3f}")
print(f"Stable genes: {len(stability['stable_genes'])}")
```

---

## Real-World Example

### Case Study: Critical Drug Target Discovery

**Scenario:** Identifying genes for expensive drug development

**Approach:**
```python
# 1. Run top 3 pipelines
pipelines = [1, 3, 5]  # STAR-DESeq2, Salmon-edgeR, STAR-limma

# 2. Require unanimous agreement
ensemble = analyzer.combine(
    results_dir='results/',
    method='vote',
    min_agreement=3  # ALL pipelines must agree
)

consensus = ensemble.get_consensus(
    min_agreement=3,
    min_confidence=0.95,  # Very high confidence
    min_abs_log2fc=1.5  # Strong effect size
)

# 3. Additional filtering
# Only genes with consistent direction
consistent_dir = consensus[consensus['direction_agreement'] == 1.0]

# Only genes significant in all individual analyses
high_sig = consistent_dir[consistent_dir['max_fdr'] < 0.01]

print(f"Final high-confidence targets: {len(high_sig)}")
high_sig.to_csv('drug_targets_high_confidence.csv', index=False)
```

**Result:**
```
Started with: 2,345 unique genes across 3 pipelines
After ensemble (3/3 agreement): 892 genes
After confidence filter (>0.95): 687 genes
After fold change filter (|FC|>1.5): 234 genes
After consistent direction: 234 genes (no change)
After stringent FDR (<0.01): 178 genes

Final high-confidence targets: 178 genes

→ These 178 genes are VERY reliable for expensive follow-up!
```

---

## Best Practices

### Do's:
✅ Use 3-5 pipelines for ensemble (sweet spot)  
✅ Choose diverse methodological approaches  
✅ Weight by ML confidence scores  
✅ Report both consensus and discordant genes  
✅ Validate key findings independently  
✅ Document all pipeline parameters  

### Don'ts:
❌ Use >8 pipelines (diminishing returns)  
❌ Use only similar pipelines (e.g., all pseudo-alignment)  
❌ Ignore discordant results completely  
❌ Apply ensemble to poor quality data  
❌ Trust ensemble blindly without validation  
❌ Mix incompatible pipeline versions  

---

## Troubleshooting

### Problem: Low Consensus

**Symptoms:** < 50% genes overlap between pipelines

**Possible causes:**
1. Poor data quality
2. Weak biological signal
3. Incompatible parameters across pipelines

**Solutions:**
```python
# Check quality first
from raptor import QualityAssessor

assessor = QualityAssessor()
qc = assessor.assess(counts, metadata)

if qc['overall_score'] < 70:
    print("WARNING: Low quality data")
    print("Consider data cleanup before ensemble")

# Check if parameters are consistent
ensemble.check_parameter_consistency()

# Use more lenient consensus
ensemble = analyzer.combine(min_agreement=2)  # Instead of 3
```

### Problem: Too Many Discordant Genes

**Investigation:**
```python
# Analyze discordance patterns
discordant = ensemble.identify_discordant()

# Group by reason
print("\nDiscordance breakdown:")
discordant.groupby('reason').size().sort_values(ascending=False)

# Check if specific pipeline is outlier
pipeline_agreement = ensemble.calculate_pipeline_agreement()
print(pipeline_agreement)
```

**Common solutions:**
- Remove outlier pipeline if one disagrees with all others
- Adjust parameters for consistency
- Use weighted ensemble giving less weight to problematic pipeline

---

## Advanced: Custom Ensemble Methods

### Custom Weighting Function

```python
def confidence_weighted_score(gene_results, pipeline_weights):
    """
    Custom scoring: weight by both pipeline quality and gene confidence.
    """
    total_score = 0
    total_weight = 0
    
    for pipeline, result in gene_results.items():
        weight = pipeline_weights[pipeline]
        # Boost weight if gene has low p-value in this pipeline
        if result['pvalue'] < 0.001:
            weight *= 2
        
        score = -np.log10(result['pvalue'])
        total_score += score * weight
        total_weight += weight
    
    return total_score / total_weight if total_weight > 0 else 0

# Use custom function
ensemble = analyzer.combine(
    results_dir='results/',
    method='custom',
    custom_function=confidence_weighted_score,
    weights={'pipeline1': 0.4, 'pipeline2': 0.3, 'pipeline3': 0.3}
)
```

---

## Summary

You've learned to:
- ✅ Choose appropriate pipelines for ensemble
- ✅ Combine results using multiple methods
- ✅ Calculate and interpret confidence scores
- ✅ Identify and understand discordant results
- ✅ Visualize ensemble overlap and concordance
- ✅ Export high-confidence gene lists
- ✅ Validate ensemble results
- ✅ Apply ensemble for critical research decisions

---

## Next Steps

1. **Apply to your data**: Run ensemble on your critical projects
2. **Read**: [ENSEMBLE.md](../ENSEMBLE.md) - Complete ensemble guide
3. **Try Tutorial 7**: [Resource Optimization](tutorial_07_resources.md)
4. **Advanced**: Create custom ensemble methods for your needs

---

## Citation

If you use ensemble analysis in your research:

```
We employed an ensemble approach combining results from [3] 
RNA-seq analysis pipelines (Salmon-edgeR, STAR-RSEM-DESeq2, 
and Kallisto-Sleuth) using RAPTOR v2.1.1 (Bolouki, 2025). 
Genes identified as differentially expressed by at least 2 
out of 3 pipelines (consensus threshold) were considered 
high-confidence candidates. This ensemble approach yielded 
[N] high-confidence genes (FDR < 0.05, |log2FC| > 1), 
providing robust results independent of methodological choice.
```

---

**Tutorial by Ayeh Bolouki**  
For RAPTOR v2.1.1

*"In statistical consensus, there is confidence!"* 
