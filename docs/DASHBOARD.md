# 🦖 RAPTOR v2.1.1 Dashboard Guide

**Interactive Web Interface for RNA-seq Analysis**

The RAPTOR Dashboard provides a powerful, no-code interface for RNA-seq pipeline selection and analysis. Perfect for users who prefer visual, interactive workflows!

---

## 🆕 What's New in v2.1.1

###  Threshold Optimizer Page
A brand new dashboard page for **data-driven threshold selection**:
- Upload DE results from any pipeline
- Optimize significance thresholds automatically
- Interactive volcano plots with optimized cutoffs
- Download publication-ready methods text
- Compare multiple adjustment methods

---

##  Table of Contents

1. [Quick Start](#quick-start)
2. [Dashboard Features](#dashboard-features)
3. [Getting Started](#getting-started)
4. [Data Upload](#data-upload)
5. [Data Profiling](#data-profiling)
6. [ML Recommendations](#ml-recommendations)
7. [🎯 Threshold Optimizer](#threshold-optimizer) - **NEW!**
8. [Pipeline Comparison](#pipeline-comparison)
9. [Results Visualization](#results-visualization)
10. [Advanced Features](#advanced-features)
11. [Multi-User Setup](#multi-user-setup)
12. [Customization](#customization)
13. [Troubleshooting](#troubleshooting)

---

##  Quick Start

### Launch Dashboard (Local)

```bash
# Start dashboard
raptor dashboard

# Dashboard opens at http://localhost:8501
```

### Launch Dashboard (Custom Port)

```bash
# Use different port
raptor dashboard --port 8502

# Accessible from network
raptor dashboard --host 0.0.0.0 --port 8501
```

### Launch Dashboard (With Data)

```bash
# Pre-load results
raptor dashboard --results /path/to/raptor_results/

# Pre-load configuration
raptor dashboard --config my_config.yaml
```

**That's it!** Your browser will open automatically. 

---

##  Dashboard Features

### Core Capabilities

✅ **Data Upload** - Drag & drop count matrices  
✅ **Visual Profiling** - Interactive data exploration  
✅ **ML Recommendations** - Get pipeline suggestions  
✅ **🎯 Threshold Optimizer** - Data-driven significance thresholds (**NEW!**)  
✅ **Real-Time Analysis** - Monitor running pipelines  
✅ **Result Visualization** - Interactive plots and tables  
✅ **Report Generation** - Export publication-ready reports  
✅ **No Coding Required** - 100% point-and-click interface  
✅ **Mobile Responsive** - Works on tablets and phones  

### Advanced Features

✅ **Multi-Dataset Comparison** - Compare multiple experiments  
✅ **Custom Themes** - Light/dark mode  
✅ **Collaboration** - Share dashboards with team  
✅ **Version Control** - Track analysis history  
✅ **Cloud Integration** - Deploy analyses to AWS/GCP/Azure  
✅ **API Integration** - Programmatic access  

---

##  Getting Started

### 1. Launch Dashboard

```bash
raptor dashboard
```

**You'll see:**
```
  Welcome to RAPTOR Dashboard! 🦖

  You can now view your dashboard in your browser.

  Local URL: http://localhost:8501
  Network URL: http://192.168.1.100:8501

  ✅ Threshold Optimizer: Available (v2.1.1)

  For help, visit: https://github.com/AyehBlk/RAPTOR
```

### 2. Dashboard Home Screen

The dashboard opens with:

```
┌─────────────────────────────────────────────────┐
│  🦖 RAPTOR Dashboard v2.1.1                     │
│  RNA-seq Analysis Pipeline Selection            │
├─────────────────────────────────────────────────┤
│                                                  │
│  📁 Upload Data                                  │
│  📊 Profile & Recommend                          │
│  🎯 Threshold Optimizer  ← NEW!                  │
│  🔬 Run Analysis                                 │
│  📈 View Results                                 │
│  ⚙️  Settings                                    │
│                                                  │
└─────────────────────────────────────────────────┘
```

### 3. Quick Tutorial

Click **"Take Tour"** for an interactive walkthrough:
- Upload sample data
- Get recommendations
- Optimize thresholds (**NEW!**)
- View results
- Generate reports

**Time:** 5 minutes

---

##  Data Upload

### Upload Count Matrix

**Step 1:** Click **"Upload Data"** tab

**Step 2:** Drag & drop your files:
- `counts.csv` - Gene expression matrix
- `metadata.csv` - Sample information

**Or click "Browse Files"** to select manually.

### Supported Formats

**Count Matrix:**
```csv
gene_id,Sample1,Sample2,Sample3
ENSG00000000003,523,612,498
ENSG00000000005,89,95,102
```

**Metadata:**
```csv
sample,condition,replicate,batch
Sample1,Control,1,Batch1
Sample2,Control,2,Batch1
Sample3,Treatment,1,Batch2
```

**Also supports:**
- Tab-separated (TSV)
- Excel (XLSX) - automatically converted
- Compressed (GZ) - automatically decompressed
- RDS files (from R) - with conversion

---

## 📊 Data Profiling

Click **"Profile & Recommend"** to analyze your data.

### Interactive Profiling

**The dashboard shows:**
- Library Size Distribution
- Sample Clustering (PCA)
- Correlation Heatmap
- Quality Metrics

### Quality Metrics

```
Data Quality Assessment

┌────────────────────────────────────┐
│ Overall Quality Score: 87/100 🟢   │
├────────────────────────────────────┤
│ ✅ Library Sizes:        92/100    │
│ ✅ BCV (Variation):      85/100    │
│ ✅ Zero Inflation:       88/100    │
│ ✅ Gene Detection:       90/100    │
│ ⚠️  Outlier Samples:     75/100    │
│ ✅ Batch Effects:        95/100    │
└────────────────────────────────────┘
```

---

## 🤖 ML Recommendations

After profiling, the dashboard shows ML-powered recommendations:

```
┌─────────────────────────────────────────────────────────┐
│  ML-POWERED PIPELINE RECOMMENDATIONS                  │
├─────────────────────────────────────────────────────────┤
│                                                          │
│  #1 RECOMMENDED: Pipeline 3 (Salmon-edgeR)            │
│                                                          │
│    ML Confidence: 89% 🟢 HIGH                           │
│    Overall Score: 0.88/1.00                             │
│                                                          │
│    ✓ Excellent for your data characteristics           │
│    ✓ Handles medium BCV well (yours: 0.42)             │
│    ✓ Fast turnaround time                              │
│                                                          │
│    💡 Next step: Optimize thresholds with ATO          │
│                                                          │
└─────────────────────────────────────────────────────────┘
```

---

## 🎯 Threshold Optimizer

**NEW in v2.1.1!** The Threshold Optimizer page provides data-driven threshold selection for differential expression analysis.

### Accessing the Threshold Optimizer

Click **"🎯 Threshold Optimizer"** in the sidebar navigation.

### Page Overview

```
┌─────────────────────────────────────────────────────────┐
│  🎯 Adaptive Threshold Optimizer                        │
│  Data-Driven Threshold Selection                        │
├─────────────────────────────────────────────────────────┤
│                                                          │
│  📤 Upload DE Results                                    │
│  ┌─────────────────────────────────────────────────┐    │
│  │  Drag and drop your DE results file here         │    │
│  │  (CSV, TSV, or TXT with logFC, pvalue columns)   │    │
│  └─────────────────────────────────────────────────┘    │
│                                                          │
│  🎲 Or generate demo data: [Generate Demo Data]         │
│                                                          │
└─────────────────────────────────────────────────────────┘
```

### Step-by-Step Usage

#### Step 1: Upload or Generate Data

**Upload your DE results:**
- Supports DESeq2, edgeR, limma output formats
- Auto-detects column names (log2FoldChange, logFC, pvalue, etc.)
- Preview data before analysis

**Or use demo data:**
```
🎲 Generate Demo Data

Number of genes: [10000    ▼]

[Generate Demo Data]

✅ Generated 10,000 genes with:
   • 800 true DE genes (8%)
   • π₀ ≈ 0.92
   • Mixed effect sizes
```

#### Step 2: Configure Analysis

```
⚙️ Analysis Settings

Analysis Goal:
◉ Discovery (maximize sensitivity)
○ Balanced (standard analysis)  
○ Validation (maximize specificity)

LogFC Method:
[Auto (Consensus) ▼]
• Auto (Consensus) - Recommended
• MAD-based
• Mixture Model
• Power-based
• Percentile

Column Mapping:
LogFC Column:  [log2FoldChange ▼]
P-value Column: [pvalue ▼]
```

#### Step 3: Run Optimization

```
[🚀 Optimize Thresholds]

⏳ Running optimization...
   ├─ Estimating π₀... ✓
   ├─ Calculating logFC threshold... ✓
   ├─ Applying p-value adjustment... ✓
   └─ Generating visualizations... ✓

✅ Optimization Complete!
```

#### Step 4: View Results

```
📊 Optimization Results

┌─────────────────────────────────────────────────────────┐
│                    SUMMARY                               │
├─────────────────────────────────────────────────────────┤
│  📈 Optimal logFC Threshold:     0.73                   │
│  📉 P-value Threshold:           0.05 (BH adjusted)     │
│  🧬 Significant Genes:           1,247                  │
│  📊 π₀ Estimate:                 0.82                   │
│  🎯 Goal:                        Discovery              │
└─────────────────────────────────────────────────────────┘
```

### Interactive Visualizations

#### Volcano Plot with Optimized Thresholds

Interactive volcano plot showing:
- Significant genes (up/down) highlighted
- Optimized threshold lines
- Hover for gene details
- Zoom and pan controls

#### P-value Distribution

Histogram showing p-value distribution with π₀ estimation line.

#### LogFC Distribution

Distribution plot with optimized threshold marked.

### Download Options

```
📥 Downloads

┌─────────────────────────────────────────────────────────┐
│  [📊 Download Optimized Results (CSV)]                  │
│      • Full results with significance flags             │
│                                                          │
│  [📝 Download Methods Text]                              │
│      • Publication-ready paragraph                       │
│                                                          │
│  [📋 Download Full Report (HTML)]                        │
│      • Complete analysis report                          │
└─────────────────────────────────────────────────────────┘
```

### Methods Text Preview

```
📝 Publication Methods Text

"Differential expression significance thresholds were 
determined using the Adaptive Threshold Optimizer (ATO) 
with the 'discovery' analysis goal. The proportion of 
true null hypotheses (π₀) was estimated at 0.82 using 
Storey's spline method. An adjusted p-value threshold 
of 0.05 (Benjamini-Hochberg FDR correction) and log₂ 
fold change threshold of 0.73 (determined by MAD-based 
robust estimation from the data) were applied, 
identifying 1,247 differentially expressed genes 
(623 upregulated, 624 downregulated)."

[Copy to Clipboard] [Download as TXT]
```

---

## 📈 Advanced Features

### 1. Ensemble Analysis

Combine results from multiple pipelines for robust consensus.

### 2. Automated Parameter Tuning

```
 Automated Parameter Tuning

💡 Tip: Use Threshold Optimizer for data-driven thresholds!

Pipeline: [Salmon-edgeR ▼]
Optimization Method: ◉ Grid Search ○ Bayesian

[Start Optimization →]
```

### 3. Cloud Deployment

Deploy to AWS, GCP, or Azure for large-scale analyses.

---

## 👥 Multi-User Setup

### Team Deployment

```bash
# Deploy on server
raptor dashboard \
  --host 0.0.0.0 \
  --port 8501 \
  --auth-required \
  --multi-user
```

### User Roles

| Role | Upload Data | Run Analysis | Threshold Optimizer | Admin |
|------|-------------|--------------|---------------------|-------|
| Admin | ✅ | ✅ | ✅ | ✅ |
| Analyst | ✅ | ✅ | ✅ | ❌ |
| Viewer | ❌ | ❌ | ✅ (view only) | ❌ |

---

## 🎨 Customization

### Themes

Light/Dark mode and custom color schemes available in Settings.

### Custom Plots

Add custom visualizations with filters and settings.

---

## 🔧 Troubleshooting

### Common Issues

#### Issue: Dashboard won't start

```bash
# Use different port
raptor dashboard --port 8502

# Or kill existing process
lsof -ti:8501 | xargs kill -9
raptor dashboard
```

#### Issue: Threshold Optimizer not available

```bash
# Reinstall/update RAPTOR
pip install --upgrade raptor-rnaseq

# Verify installation
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

#### Issue: Upload fails

- Check file is CSV/TSV
- Ensure UTF-8 encoding
- No special characters in column names
- Save as CSV from Excel before uploading

#### Issue: Slow performance

```bash
# Limit data displayed
raptor dashboard --max-samples 50

# Use summary mode
raptor dashboard --lightweight
```

### Debug Mode

```bash
# Enable debug logging
raptor dashboard --debug

# Check logs
tail -f ~/.raptor/dashboard.log
```

---

## 🔌 API Integration

```python
from raptor.dashboard import DashboardAPI
from raptor.threshold_optimizer import optimize_thresholds

# Start dashboard
api = DashboardAPI()
api.start(port=8501)

# Run threshold optimization programmatically
import pandas as pd
de_results = pd.read_csv('de_results.csv')
result = optimize_thresholds(de_results, goal='discovery')
print(f"Optimal thresholds: |logFC| > {result.logfc_threshold:.2f}")

# Stop dashboard
api.stop()
```

---

## 📚 Best Practices

### For Lab Managers

1. Set up team dashboard on shared server
2. Configure user roles appropriately
3. Document standard workflows
4. Train users on Threshold Optimizer

### For Analysts

1. Always validate uploaded data
2. Review QC plots before proceeding
3. **Use Threshold Optimizer** for data-driven thresholds
4. Compare multiple pipelines when uncertain
5. Export and save analyses regularly

### For Publications

1. Use Threshold Optimizer for defensible thresholds
2. Include the auto-generated methods text
3. Export high-resolution figures
4. Document all analysis parameters

---

## 📖 Learning Resources

### Documentation

- [User Guide](USER_GUIDE.md)
- [Threshold Optimizer Guide](THRESHOLD_OPTIMIZER.md) - **NEW!**
- [API Reference](API.md)
- [FAQ](FAQ.md)

### Interactive Demo

```bash
# Launch demo with sample data
raptor dashboard --demo
```

---

## 🆘 Support

**Dashboard not working?**

1. Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
2. Read [FAQ](FAQ.md)
3. GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
4. Email: ayehbolouki1988@gmail.com

---

## 📋 Summary

The RAPTOR Dashboard provides:
- ✅ **No-code interface** for RNA-seq analysis
- ✅ **Interactive visualizations** for data exploration
- ✅ **ML-powered recommendations** with explanations
- ✅ **🎯 Threshold Optimizer** for data-driven thresholds (**NEW!**)
- ✅ **Real-time monitoring** of analyses
- ✅ **Team collaboration** features
- ✅ **Export capabilities** for publications
- ✅ **Cloud integration** for large projects

**Perfect for researchers who prefer visual, interactive workflows!** 🚀

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.1  
**License:** MIT

---

*"Making RNA-seq accessible to everyone, one click at a time!"* 🦖
