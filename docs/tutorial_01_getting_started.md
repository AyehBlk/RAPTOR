# Tutorial 1: Getting Started with RAPTOR

**Level**: Beginner  
**Time**: 30 minutes  
**Goal**: Learn the basics of RAPTOR and run your first analysis

---

## What You'll Learn

- How to use RAPTOR's command-line interface
- How to profile your RNA-seq data
- How to get pipeline recommendations
- How to interpret the results

---

## Prerequisites

- RAPTOR installed (see [INSTALLATION.md](../INSTALLATION.md))
- Basic command-line knowledge
- RNA-seq count matrix (or use our simulated data)

---

## Step 1: Verify Installation

First, let's make sure RAPTOR is properly installed:

```bash
# Check RAPTOR version
raptor --version

# See available commands
raptor --help

# Check installed tools
raptor --check-tools
```

**Expected output:**
```
RAPTOR v2.0.0
RNA-seq Analysis Pipeline Testing and Optimization Resource
```

---

## Step 2: Generate Simulated Data

For this tutorial, we'll use simulated data with known ground truth:

```bash
# Create a test directory
mkdir raptor_tutorial
cd raptor_tutorial

# Generate simulated data (small dataset)
raptor simulate --output simulated_data/ --size small

# Check what was created
ls simulated_data/
```

**You should see:**
- `counts.csv` - Gene expression matrix
- `metadata.csv` - Sample information
- `ground_truth.csv` - True differentially expressed genes

### Understanding the Simulated Data

```bash
# View the count matrix (first few lines)
head simulated_data/counts.csv

# View metadata
cat simulated_data/metadata.csv
```

The data contains:
- **1000 genes** × **6 samples** (3 control, 3 treatment)
- **200 differentially expressed genes** (20%)
- **Fold changes**: 0.5× (down-regulated) and 2.0× (up-regulated)

---

## Step 3: Profile Your Data

Now let's use RAPTOR to analyze the data characteristics:

```bash
# Run data profiling
raptor profile \
  --counts simulated_data/counts.csv \
  --metadata simulated_data/metadata.csv \
  --output results/
```

### What's Happening?

RAPTOR analyzes:
1. **Library sizes** - Total counts per sample
2. **Biological variation** - Coefficient of variation
3. **Zero-inflation** - Percentage of zero counts
4. **Sample size** - Number of samples per condition
5. **Sequencing depth** - Average read depth
6. **Gene detection** - Number of expressed genes

---

## Step 4: Review the Results

RAPTOR generates an HTML report. Let's look at it:

```bash
# Open the report (Linux)
xdg-open results/raptor_profile_report.html

# macOS
open results/raptor_profile_report.html

# Windows (WSL)
explorer.exe results/raptor_profile_report.html
```

### The Report Contains:

1. **Executive Summary**
   - Quick overview of your data
   - Key characteristics
   - **Recommended pipeline**

2. **Data Overview**
   - Library size distribution
   - Sample clustering (PCA plot)
   - Correlation heatmap

3. **Profiling Results**
   - Detailed statistics
   - Quality metrics
   - Data characteristics table

4. **Pipeline Recommendations**
   - Top 3 recommended pipelines
   - Reasoning for each recommendation
   - Expected performance

---

## Step 5: Understanding the Recommendation

Let's examine the recommendation more carefully.

**Example output:**
```
═══════════════════════════════════════════════════
           PIPELINE RECOMMENDATION
═══════════════════════════════════════════════════

🥇 RECOMMENDED: Pipeline 3 (Salmon-edgeR)

   Score: 0.89 / 1.00

   Why this pipeline?
   ✓ Excellent balance of speed and accuracy
   ✓ Well-suited for your sample size (n=6)
   ✓ Handles your biological variation well
   ✓ Low resource requirements

   Expected Performance:
   • Runtime: ~15 minutes
   • Memory: ~8 GB
   • Accuracy: F1 = 0.88

───────────────────────────────────────────────────

🥈 Alternative: Pipeline 1 (STAR-RSEM-DESeq2)

   Score: 0.85 / 1.00

   Why consider this?
   ✓ Highest accuracy (F1 = 0.95)
   ✓ Best for difficult datasets
   ⚠ Slower (2-3 hours)
   ⚠ Higher memory (32 GB)

───────────────────────────────────────────────────

🥉 Alternative: Pipeline 4 (Kallisto-Sleuth)

   Score: 0.82 / 1.00

   Why consider this?
   ✓ Fastest option (~10 minutes)
   ✓ Very low memory (4 GB)
   ⚠ Lower accuracy for small samples
```

### Key Information:

- **Score**: Overall suitability (0-1 scale)
- **Why**: Reasoning based on your data characteristics
- **Expected Performance**: Estimates for your data size
- **Alternatives**: Other good options to consider

---

## Step 6: Run the Recommended Pipeline (Optional)

If you want to actually run the analysis:

```bash
# Run Pipeline 3 on the simulated data
raptor run \
  --pipeline 3 \
  --data simulated_data/ \
  --output pipeline3_results/ \
  --threads 4
```

**Note**: This requires FASTQ files. For count matrices, RAPTOR starts from the quantification step.

---

## Step 7: Compare with Ground Truth

Since we used simulated data, we can check accuracy:

```bash
# Compare results with known truth
raptor evaluate \
  --results pipeline3_results/diff_expr.csv \
  --truth simulated_data/ground_truth.csv \
  --output evaluation/
```

**Output:**
```
Performance Metrics:
  Sensitivity: 0.87
  Specificity: 0.92
  Precision: 0.85
  F1-Score: 0.86
  True Positives: 174 / 200
  False Positives: 32
```

---

## Step 8: Try with Your Own Data

Now that you understand the basics, try with your real data:

```bash
# Profile your data
raptor profile \
  --counts /path/to/your/counts.csv \
  --metadata /path/to/your/metadata.csv \
  --output my_analysis/
```

### Data Format Requirements:

**Count Matrix (`counts.csv`):**
```csv
gene_id,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6
Gene1,523,612,498,1250,1180,1340
Gene2,89,95,102,210,185,198
Gene3,2341,2567,2234,2401,2389,2456
...
```

**Metadata (`metadata.csv`):**
```csv
sample,condition,replicate
Sample1,Control,1
Sample2,Control,2
Sample3,Control,3
Sample4,Treatment,1
Sample5,Treatment,2
Sample6,Treatment,3
```

---

## Common Questions

### Q1: What if I only have FASTQ files?

You can still use RAPTOR to run complete pipelines from raw data.

### Q2: How long does profiling take?

Very fast! Usually < 1 minute for datasets with up to 50 samples and 25,000 genes.

### Q3: Can I customize the recommendation criteria?

Yes! Edit the configuration file:
```bash
# Copy default config
cp config/config.yaml my_config.yaml

# Edit weights
nano my_config.yaml

# Use custom config
raptor profile --counts data.csv --config my_config.yaml
```

### Q4: What if I disagree with the recommendation?

The recommendation is a suggestion based on data characteristics. You can:
1. Review alternative recommendations
2. Run full benchmarking to compare all pipelines
3. Choose any pipeline manually with `raptor run --pipeline X`

---

## Next Steps

Now that you've completed this tutorial:

1. **Try Tutorial 2**: [Running Full Benchmarking](tutorial_02_benchmarking.md)
2. **Explore different data types**: Try with different sample sizes and conditions
3. **Read the documentation**: [PROFILE_RECOMMEND.md](../PROFILE_RECOMMEND.md)
4. **Understand pipelines**: [PIPELINES.md](../PIPELINES.md)

---

## Troubleshooting

**Problem**: `raptor: command not found`
```bash
# Solution: Add to PATH
export PATH=$PATH:~/.local/bin
```

**Problem**: HTML report won't open
```bash
# Solution: Open manually in browser
firefox results/raptor_profile_report.html
```

**Problem**: Out of memory error
```bash
# Solution: Use less memory
raptor profile --counts data.csv --memory 8G
```

See [TROUBLESHOOTING.md](../TROUBLESHOOTING.md) for more help.

---

## Summary

You've learned to:
- ✅ Generate simulated RNA-seq data
- ✅ Profile data characteristics
- ✅ Get intelligent pipeline recommendations
- ✅ Interpret RAPTOR's HTML reports
- ✅ Evaluate results against ground truth

**Congratulations!** 🎉 You're ready to use RAPTOR for your RNA-seq analysis!

---

**Tutorial by Ayeh Bolouki**  
University of Namur, Belgium  
For RAPTOR v2.0.0
