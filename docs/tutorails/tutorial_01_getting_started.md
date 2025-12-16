# Tutorial 01: Getting Started with RAPTOR v2.1.1

**Level**: Beginner  
**Time**: 30-45 minutes  
**Goal**: Learn basic RAPTOR operations with ML recommendations and interactive dashboard

---

##  What You'll Learn

By the end of this tutorial, you'll be able to:
- ✅ Run RAPTOR's ML-powered pipeline recommendations
- ✅ Use the interactive web dashboard
- ✅ Profile your data with quality assessment
- ✅ Understand pipeline selection based on data characteristics
- ✅ Compare multiple pipelines automatically
- ✅ **NEW in v2.1.1**: Apply data-driven threshold optimization

---

## 🆕 What's New in v2.1.1

**Adaptive Threshold Optimizer (ATO)** - Replace arbitrary thresholds with data-driven cutoffs:

```bash
# Quick ATO example
raptor optimize-thresholds --results deseq2_results.csv --goal balanced
```

ATO determines optimal |logFC| thresholds based on YOUR data characteristics, not arbitrary values like "1.0". Learn more in [Tutorial 07](tutorial_07_threshold_optimizer.md).

---

##  Prerequisites

**Required**:
- RAPTOR v2.1.1 installed ([see INSTALLATION.md](../INSTALLATION.md))
- Basic command line knowledge
- Count matrix or FASTQ files

**Optional**:
- Python 3.9+ for advanced features
- Web browser for dashboard

**Verify ATO is available**:
```bash
raptor --version
# Should show: RAPTOR v2.1.1

python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

---

##  Test Dataset

We'll use a simulated dataset similar to real RNA-seq experiments:

```bash
# Download test data (or use your own)
mkdir -p ~/raptor_tutorial
cd ~/raptor_tutorial

# Create a simple count matrix
cat > counts_simple.csv << 'EOF'
gene_id,sample1,sample2,sample3,sample4,sample5,sample6
GENE001,1250,1180,1340,45,52,48
GENE002,890,920,875,2100,2050,2180
GENE003,2300,2400,2250,2280,2350,2420
GENE004,45,52,48,1250,1180,1340
GENE005,150,145,160,155,148,162
EOF
```

This represents:
- 5 genes
- 6 samples (3 control + 3 treated)
- Mix of expression patterns (up/down/stable)

---

##  Part 1: Quick Start with ML Recommendations

### Step 1: Get Smart Pipeline Recommendations

The simplest way to start - let RAPTOR's ML analyze your data and suggest the best pipeline:

```bash
raptor profile --counts counts_simple.csv \
               --groups control,control,control,treated,treated,treated \
               --recommend \
               --output profile_report.txt
```

**What just happened?**
- ✅ RAPTOR analyzed your data characteristics
- ✅ ML model predicted which pipeline performs best
- ✅ Generated a detailed recommendation report

### Step 2: Read Your Recommendation

```bash
cat profile_report.txt
```

You'll see something like:

```
=== RAPTOR v2.1.1 Profile Report ===
Dataset: counts_simple.csv

 DATA CHARACTERISTICS:
  - Genes: 5
  - Samples: 6
  - Library size: 960-1064 (low variance)
  - Expression range: High dynamic range
  - Zero inflation: 0% (excellent)
  
 ML RECOMMENDATION:
  Recommended Pipeline: STAR-DESeq2
  Confidence: 87%
  
  Why this pipeline?
  ✓ Optimal for your sample size (n=3 per group)
  ✓ Handles high dynamic range well
  ✓ Best performance for count-based analysis
  
 TOP 3 PIPELINES (ML-ranked):
  1. STAR-DESeq2 (score: 8.7/10)
  2. Salmon-DESeq2 (score: 8.5/10)
  3. Kallisto-DESeq2 (score: 8.2/10)
```

**Key insights from your profile:**
- **Library size**: How many reads per sample
- **Zero inflation**: Proportion of genes with zero counts
- **Dynamic range**: Difference between high/low expressed genes
- **ML confidence**: How certain the model is about its recommendation

---

##  Part 2: Interactive Dashboard

### Launch the Dashboard

The new v2.1.0 dashboard makes RAPTOR accessible without command line:

```bash
raptor dashboard
```

This opens in your web browser at `http://localhost:5000`

### Dashboard Features

**1. Data Upload Tab**:
- Drag-and-drop your count matrix
- Specify sample groups visually
- See instant data preview

**2. Profile & Recommend Tab**:
- Click "Analyze Data"
- View ML recommendations
- Interactive visualization of data characteristics
- Download recommendations as PDF

**3. Run Analysis Tab**:
- Select recommended pipeline(s)
- Set parameters with helpful tooltips
- Monitor progress in real-time
- See resource usage (CPU, memory)

**4. Results Tab**:
- Interactive MA plots
- Volcano plots
- Download results tables
- Export publication-ready figures

**5. 🎯 Threshold Optimizer Tab** (NEW in v2.1.1):
- Upload DE results
- Select analysis goal (discovery/balanced/validation)
- Get data-driven thresholds
- Download methods text for publications

### Try It Out

1. **Upload Data**:
   - Click "Choose File" → select `counts_simple.csv`
   - Enter groups: `control,control,control,treated,treated,treated`

2. **Get Recommendations**:
   - Click "Profile & Recommend"
   - View the ML-powered suggestions
   - See confidence scores

3. **Run Analysis** (optional):
   - Select "STAR-DESeq2" (recommended)
   - Click "Run Pipeline"
   - Watch real-time progress

**Pro tip**: The dashboard saves your session, so you can close it and come back later!

---

## 🎯 Part 3: Threshold Optimization (NEW in v2.1.1)

After running your pipeline, optimize your significance thresholds:

### Why Data-Driven Thresholds?

Traditional approach uses arbitrary cutoffs:
```python
# Why 1.0? Why not 0.5 or 1.5? 🤷
significant = results[(abs(results['logFC']) > 1.0) & (results['padj'] < 0.05)]
```

ATO determines optimal thresholds based on YOUR data:

### Command Line

```bash
raptor optimize-thresholds --results results/degs.csv \
                           --goal balanced \
                           --output threshold_results/
```

### Python

```python
from raptor.threshold_optimizer import optimize_thresholds
import pandas as pd

de_results = pd.read_csv('results/degs.csv')

result = optimize_thresholds(
    de_results,
    logfc_col='log2FoldChange',
    pvalue_col='pvalue',
    goal='balanced'
)

print(f"Optimal |logFC|: {result.logfc_threshold:.3f}")
print(f"Significant genes: {result.n_significant}")
print(f"\nMethods text:\n{result.methods_text}")
```

### Dashboard

In the 🎯 Threshold Optimizer tab:
1. Upload your DE results
2. Select the appropriate columns
3. Choose your goal
4. Click "Optimize"
5. Download results and methods text

**Learn more**: [Tutorial 07: Threshold Optimization](tutorial_07_threshold_optimizer.md)

---

##  Part 4: Command-Line Workflow (Traditional)

If you prefer command-line or need automation:

### Profile Your Data

```bash
raptor profile --counts counts_simple.csv \
               --groups control,control,control,treated,treated,treated \
               --output profile_detailed.json \
               --format json
```

### Get Detailed ML Insights

```bash
raptor ml-recommend --profile profile_detailed.json \
                    --top-n 3 \
                    --explain
```

Output:
```json
{
  "top_recommendations": [
    {
      "pipeline": "STAR-DESeq2",
      "score": 8.7,
      "confidence": 0.87,
      "reasons": [
        "Optimal for sample size n=3",
        "Best for count-based data",
        "Robust to outliers"
      ],
      "expected_sensitivity": 0.89,
      "expected_specificity": 0.91
    }
  ],
  "data_profile": {
    "difficulty": "easy",
    "noise_level": "low",
    "complexity": "simple"
  }
}
```

### Run the Recommended Pipeline

```bash
raptor run --pipeline STAR-DESeq2 \
           --counts counts_simple.csv \
           --groups control,control,control,treated,treated,treated \
           --output results/
```

---

## ✅ Part 5: Quality Assessment (New in v2.1.0)

Check if your data meets quality standards:

```bash
raptor assess-quality --counts counts_simple.csv \
                      --output quality_report.html
```

Open `quality_report.html` in your browser to see:
- ✅ Library size distribution
- ✅ Gene detection rates
- ✅ Sample correlation heatmap
- ✅ PCA plot
- ✅ Quality flags and warnings

**Common quality issues**:
- ❌ Low library size (< 1 million reads)
- ❌ High zero inflation (> 70%)
- ❌ Outlier samples
- ❌ Batch effects

---

##  Part 6: Compare Multiple Pipelines

Want to see how different pipelines perform?

### Quick Comparison (3 pipelines)

```bash
raptor compare --counts counts_simple.csv \
               --groups control,control,control,treated,treated,treated \
               --pipelines STAR-DESeq2,Salmon-DESeq2,Kallisto-DESeq2 \
               --output comparison/
```

### Full Benchmark (all 8 pipelines)

```bash
raptor benchmark --counts counts_simple.csv \
                 --groups control,control,control,treated,treated,treated \
                 --output benchmark_results/
```

**What you get**:
- Concordance matrix (how similar are results?)
- Performance metrics for each pipeline
- ML-based ranking
- Recommendation based on YOUR data

---

##  Understanding Your Results

### Key Output Files

```
results/
├── degs.csv                    # Differentially expressed genes
├── plots/
│   ├── ma_plot.png            # MA plot
│   ├── volcano_plot.png       # Volcano plot
│   └── heatmap.png            # Top genes heatmap
├── quality_report.html        # Data quality assessment
└── ml_recommendation.json     # ML analysis results
```

### Interpreting DEGs

Look at `degs.csv`:

```csv
gene_id,baseMean,log2FoldChange,pvalue,padj,significant
GENE001,701.7,-4.71,0.0001,0.0005,TRUE
GENE002,1519.2,4.68,0.0001,0.0005,TRUE
GENE003,2333.3,0.12,0.8234,0.9234,FALSE
```

- **log2FoldChange > 1**: Upregulated (treated > control)
- **log2FoldChange < -1**: Downregulated (treated < control)
- **padj < 0.05**: Statistically significant
- **significant = TRUE**: Meets both criteria

**Note**: With ATO (v2.1.1), your logFC threshold may differ from 1.0 based on your data!

---

##  Best Practices

### 1. Always Profile First
```bash
# Don't skip this!
raptor profile --counts yourdata.csv --recommend
```

### 2. Check Quality
```bash
# Catch problems early
raptor assess-quality --counts yourdata.csv
```

### 3. Use Data-Driven Thresholds (v2.1.1)
```bash
# Replace arbitrary |logFC| > 1
raptor optimize-thresholds --results degs.csv --goal balanced
```

### 4. Use the Dashboard for Exploration
```bash
# Great for interactive analysis
raptor dashboard
```

### 5. Use Command Line for Automation
```bash
# For pipelines and scripts
raptor run --pipeline $(raptor recommend --counts data.csv --best)
```

### 6. Compare When Unsure
```bash
# If ML confidence < 75%, compare top 3
raptor compare --pipelines top3 --counts data.csv
```

---

## ❓ Common Issues & Solutions

### "No recommendations available"
**Problem**: Not enough data characteristics  
**Solution**: 
```bash
raptor profile --counts yourdata.csv --verbose --detailed
```

### "Low quality score"
**Problem**: Data quality issues detected  
**Solution**: Check quality report
```bash
raptor assess-quality --counts yourdata.csv
# Look for warnings about library size, outliers, etc.
```

### Dashboard won't start
**Problem**: Port already in use  
**Solution**:
```bash
raptor dashboard --port 5001
# Or kill existing process:
pkill -f "raptor dashboard"
```

### ML recommendation seems wrong
**Problem**: Edge case not in training data  
**Solution**: Run ensemble analysis
```bash
raptor ensemble --counts yourdata.csv --pipelines all
# This combines multiple pipelines for robust results
```

### Threshold seems too permissive/stringent
**Problem**: Default ATO goal doesn't match your needs  
**Solution**: Try different goals
```bash
raptor optimize-thresholds --results degs.csv --goal discovery   # More permissive
raptor optimize-thresholds --results degs.csv --goal validation  # More stringent
```

---

##  What's Next?

Now that you know the basics, continue learning:

### Tutorial 02: Advanced Benchmarking
- Compare all 8 pipelines systematically
- Understand concordance analysis
- Use ensemble methods
- [→ tutorial_02_benchmarking.md](tutorial_02_benchmarking.md)

### Tutorial 03: Real Data Analysis
- Work with actual RNA-seq datasets
- Handle complex experimental designs
- Optimize for your specific needs
- [→ tutorial_03_real_data.md](tutorial_03_real_data.md)

### Tutorial 04: ML Deep Dive
- Understand how ML recommendations work
- Train custom models on your data
- Fine-tune for specific organism/tissue
- [→ tutorial_04_ml_recommendations.md](tutorial_04_ml_recommendations.md)

### Tutorial 05: Dashboard Mastery
- Build custom analysis workflows
- Create shareable reports
- Integrate with other tools
- [→ tutorial_05_dashboard.md](tutorial_05_dashboard.md)

### Tutorial 06: Ensemble Methods
- Combine multiple pipelines for robust results
- Use ATO with ensemble analysis
- [→ tutorial_06_ensemble.md](tutorial_06_ensemble.md)

### Tutorial 07: Threshold Optimization 🎯 (NEW)
- Master data-driven threshold selection
- Generate publication methods text
- [→ tutorial_07_threshold_optimizer.md](tutorial_07_threshold_optimizer.md)

---

##  Reference Documentation

- **[PROFILE_RECOMMEND.md](../PROFILE_RECOMMEND.md)**: Complete profiling guide
- **[ML_GUIDE.md](../ML_GUIDE.md)**: Machine learning details
- **[DASHBOARD_GUIDE.md](../DASHBOARD_GUIDE.md)**: Dashboard features
- **[QUALITY_GUIDE.md](../QUALITY_GUIDE.md)**: Quality assessment
- **[THRESHOLD_OPTIMIZER.md](../THRESHOLD_OPTIMIZER.md)**: ATO complete guide
- **[API.md](../API.md)**: Python API reference

---

##  Quick Reference Commands

```bash
# Get ML recommendation
raptor profile --counts data.csv --recommend

# Launch dashboard
raptor dashboard

# Check quality
raptor assess-quality --counts data.csv

# Run recommended pipeline
raptor run --pipeline STAR-DESeq2 --counts data.csv --groups ...

# Optimize thresholds (NEW in v2.1.1)
raptor optimize-thresholds --results degs.csv --goal balanced

# Compare top 3 pipelines
raptor compare --pipelines top3 --counts data.csv --groups ...

# Full benchmark
raptor benchmark --counts data.csv --groups ...

# Ensemble analysis
raptor ensemble --counts data.csv --groups ... --method weighted

# Ensemble with ATO (v2.1.1)
raptor ensemble --counts data.csv --groups ... --use-ato --ato-goal balanced
```

---

## 🆘 Need Help?

- **Quick questions**: See [FAQ.md](../FAQ.md)
- **Problems**: Check [TROUBLESHOOTING.md](../TROUBLESHOOTING.md)
- **GitHub Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

## ✅ Checklist: Tutorial Complete!

Did you:
- [ ] Install RAPTOR v2.1.1
- [ ] Verify ATO is available
- [ ] Get ML recommendations for test data
- [ ] Try the interactive dashboard
- [ ] Run a quality assessment
- [ ] Try threshold optimization (NEW!)
- [ ] Compare multiple pipelines
- [ ] Understand the output files
- [ ] Know where to get help

**Congratulations! You're ready for real data analysis! 🎉**

---

**Tutorial 01 - Getting Started**  
*RAPTOR v2.1.1*  
Created by Ayeh Bolouki  
Last updated: December 2025
