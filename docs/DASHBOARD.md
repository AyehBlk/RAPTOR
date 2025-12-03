#  RAPTOR v2.1.0 Dashboard Guide

**Interactive Web Interface for RNA-seq Analysis**

The RAPTOR Dashboard provides a powerful, no-code interface for RNA-seq pipeline selection and analysis. Perfect for users who prefer visual, interactive workflows!

---

##  Table of Contents

1. [Quick Start](#quick-start)
2. [Dashboard Features](#dashboard-features)
3. [Getting Started](#getting-started)
4. [Data Upload](#data-upload)
5. [Data Profiling](#data-profiling)
6. [ML Recommendations](#ml-recommendations)
7. [Pipeline Comparison](#pipeline-comparison)
8. [Results Visualization](#results-visualization)
9. [Advanced Features](#advanced-features)
10. [Multi-User Setup](#multi-user-setup)
11. [Customization](#customization)
12. [Troubleshooting](#troubleshooting)

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

**That's it!** Your browser will open automatically. 🎉

---

##  Dashboard Features

### Core Capabilities

✅ **Data Upload** - Drag & drop count matrices  
✅ **Visual Profiling** - Interactive data exploration  
✅ **ML Recommendations** - Get pipeline suggestions  
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

  For help, visit: https://github.com/AyehBlk/RAPTOR
```

### 2. Dashboard Home Screen

The dashboard opens with:

```
┌─────────────────────────────────────────────────┐
│  🦖 RAPTOR Dashboard v2.1.0                     │
│  RNA-seq Analysis Pipeline Selection            │
├─────────────────────────────────────────────────┤
│                                                  │
│  📁 Upload Data                                  │
│  📊 Profile & Recommend                          │
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

### File Format Validation

The dashboard automatically checks:
- ✅ Correct format (CSV/TSV)
- ✅ Genes as rows, samples as columns
- ✅ No missing values
- ✅ Metadata matches count matrix
- ✅ Valid column names

**If issues found:**
```
  Validation Issues Found:

1. Missing sample in metadata: Sample_7
2. Non-integer values detected in counts
3. 5 genes have all zeros

 Fix these issues or click "Ignore" to proceed
```

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

### Data Preview

After upload, you'll see:

```
✅ Data Loaded Successfully!

 Dataset Overview:
   • Genes: 20,134
   • Samples: 12
   • Conditions: 2 (Control, Treatment)
   • Total Counts: 245,678,901

 Preview (first 5 rows, 3 columns):
┌──────────────────┬──────────┬──────────┬──────────┐
│ gene_id          │ Sample1  │ Sample2  │ Sample3  │
├──────────────────┼──────────┼──────────┼──────────┤
│ ENSG00000000003  │ 523      │ 612      │ 498      │
│ ENSG00000000005  │ 89       │ 95       │ 102      │
│ ENSG00000000419  │ 2341     │ 2567     │ 2234     │
└──────────────────┴──────────┴──────────┴──────────┘

[Continue to Profiling →]
```

---

##  Data Profiling

Click **"Profile & Recommend"** to analyze your data.

### Interactive Profiling

**The dashboard shows:**

#### 1. Library Size Distribution

```
Library Sizes (Million Reads)
    ↑
30M │     ▂▄█▆▃
25M │   ▃█████▆▂
20M │  ▅████████▄
15M │ ▇███████████▇
10M │██████████████
    └───────────────→ Samples
    
 Statistics:
   Mean: 24.5M
   Median: 25.1M
   CV: 12.3% ✅ Good variation
```

#### 2. Sample Clustering (PCA)

```
Interactive PCA Plot
(Click to zoom, drag to rotate)

    PC2 (15%)
    ↑
    │    ●Control
    │  ●●●
    │●●●
    │      ▲▲▲Treatment
    │    ▲▲▲
    └─────────────→ PC1 (68%)

 Samples cluster by condition ✅
```

#### 3. Correlation Heatmap

```
Sample Correlation
(Hover for values)

         S1   S2   S3   S4   S5   S6
    S1 [1.0  0.95 0.93 0.45 0.42 0.48]
    S2 [0.95 1.0  0.94 0.43 0.44 0.46]
    S3 [0.93 0.94 1.0  0.41 0.43 0.45]
    S4 [0.45 0.43 0.41 1.0  0.96 0.94]
    S5 [0.42 0.44 0.43 0.96 1.0  0.95]
    S6 [0.48 0.46 0.45 0.94 0.95 1.0 ]

   Within-group: High ✅
   Between-group: Low ✅
```

#### 4. Quality Metrics

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

⚠️  Warning: Sample_7 is a potential outlier
   Consider reviewing or excluding
```

### Detailed Statistics

**Expandable sections:**

```
▼ Biological Coefficient of Variation (BCV)
  Value: 0.42
  Category: Medium ✅
  
  What this means:
  Your data shows moderate biological variation,
  typical for well-controlled experiments. Most
  pipelines will perform well.
  
  📖 Learn more about BCV

▼ Sequencing Depth
  Mean: 24.5M reads
  Category: High ✅
  
  What this means:
  Excellent sequencing depth. You have sufficient
  coverage for accurate quantification of most genes.
  
   View depth distribution

▼ Zero Inflation
  Percentage: 42%
  Category: Normal ✅
  
  What this means:
  Typical proportion of zero counts for RNA-seq data.
  No special handling required.
```

---

##  ML Recommendations

After profiling, the dashboard shows ML-powered recommendations:

### Recommendation Panel

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
│     Expected Runtime: ~22 minutes                     │
│     Expected Memory: ~12 GB                           │
│     Expected Accuracy: F1 = 0.88                      │
│                                                          │
│    ✓ Excellent for your data characteristics           │
│    ✓ Handles medium BCV well (yours: 0.42)             │
│    ✓ Fast turnaround time                              │
│    ✓ Low resource requirements                         │
│    ✓ 85-90% accuracy typical                           │
│                                                          │
│    [Why this pipeline? ▼] [Run Analysis →]            │
│                                                          │
├─────────────────────────────────────────────────────────┤
│                                                          │
│  #2 Alternative: Pipeline 1 (STAR-RSEM-DESeq2)        │
│                                                          │
│    ML Confidence: 82% 🟢 HIGH                           │
│    Overall Score: 0.85/1.00                             │
│                                                          │
│     Expected Runtime: ~3.5 hours                      │
│     Expected Memory: ~45 GB                           │
│     Expected Accuracy: F1 = 0.92                      │
│                                                          │
│    Consider if:                                         │
│    • You need highest accuracy                         │
│    • You have sufficient resources                     │
│    • Time is not critical                              │
│                                                          │
│    [Details ▼] [Run Analysis →]                        │
│                                                          │
├─────────────────────────────────────────────────────────┤
│                                                          │
│  #3 Alternative: Pipeline 4 (Kallisto-Sleuth)         │
│                                                          │
│    ML Confidence: 78% 🟡 MEDIUM                         │
│    Overall Score: 0.82/1.00                             │
│                                                          │
│    Fastest option if speed is priority                 │
│                                                          │
└─────────────────────────────────────────────────────────┘

[View All 8 Pipelines →] [Compare Recommendations →]
```

### Understanding Confidence Scores

Click **"Why this pipeline?"** to see:

```
 ML Model Explanation

Top Features Influencing Recommendation:

1. Sample Size (n=12)           Impact: +12%
   ├─ Moderate size favors balanced pipelines
   └─ Salmon-edgeR optimal for 6-20 samples

2. BCV (0.42 - Medium)         Impact: +18%
   ├─ Not too high, not too low
   └─ Salmon-edgeR handles this well

3. Sequencing Depth (24.5M)    Impact: +15%
   ├─ High depth = good accuracy potential
   └─ Salmon quasi-mapping efficient

4. Library Size CV (12.3%)     Impact: +8%
   ├─ Low variation = easier normalization
   └─ All pipelines suitable

5. Zero Inflation (42%)        Impact: +5%
   ├─ Normal level
   └─ No special considerations

 Model Confidence: 89%
   (Based on similarity to 1,247 training examples)

 Similar successful projects:
   • Project A: Cancer vs Normal (n=16, BCV=0.38)
   • Project B: Treatment response (n=10, BCV=0.45)
   • Project C: Time series (n=12, BCV=0.40)
```

### What-If Scenarios

Interactive sliders to explore:

```
 Explore Different Scenarios

Sample Size: [■■■■■□□□□□] 12 samples
├─ Current: Pipeline 3 recommended
├─ If 6 samples: Pipeline 6 recommended
└─ If 50 samples: Pipeline 4 recommended

BCV: [■■■■■□□□□□] 0.42 (Medium)
├─ Current: Pipeline 3 recommended
├─ If 0.2 (Low): Pipeline 4 recommended  
└─ If 0.7 (High): Pipeline 1 recommended

Resource Limit: [■■■■■■□□□□] 16 GB
├─ Current: Multiple options available
├─ If 8 GB: Pipeline 3 or 4 only
└─ If 64 GB: All pipelines available
```

---

##  Pipeline Comparison

Click **"Compare All Pipelines"** for detailed comparison:

### Comparison Table

```
Pipeline Comparison (Interactive - Click to sort)

Pipeline              Score  Confidence  Runtime  Memory  Accuracy
─────────────────────────────────────────────────────────────────
3. Salmon-edgeR       0.88   89% 🟢      22m     12 GB   0.88
1. STAR-RSEM-DESeq2   0.85   82% 🟢      3.5h    45 GB   0.92
5. STAR-HTSeq-limma   0.83   79% 🟡      3.2h    42 GB   0.89
4. Kallisto-Sleuth    0.82   78% 🟡      15m     6 GB    0.83
2. HISAT2-StringTie   0.75   71% 🟡      2.1h    28 GB   0.82
6. STAR-NOISeq        0.72   68% 🟡      3.0h    38 GB   0.80
7. Bowtie2-RSEM       0.70   65% 🟡      2.8h    32 GB   0.81
8. HISAT2-Cufflinks   0.62   58% 🟡      4.2h    35 GB   0.75

[Export Table] [Detailed View] [Custom Weights]
```

### Visual Comparisons

**Speed vs Accuracy:**
```
Interactive Scatter Plot
(Hover for details, click for info)

Accuracy
    ↑
1.0 │       ●1
0.9 │     ●5 ●3
0.8 │   ●2   ●4
0.7 │ ●8
    └───────────────→ Runtime (log)
      10m   1h   4h

Click any point for pipeline details
```

**Resource Efficiency:**
```
Memory Usage vs Accuracy

Memory (GB)
    ↑
50  │ ●1
40  │ ●5 ●6
30  │ ●2 ●7 ●8
20  │
10  │ ●3
5   │ ●4
    └───────────────→ Accuracy
    0.7   0.8   0.9   1.0
```

### Custom Weighting

```
 Adjust Your Priorities

Accuracy:    [■■■■■■□□□□] 60%
Speed:       [■■■□□□□□□□] 30%
Memory:      [■■□□□□□□□□] 10%

[Apply Weights] → Recommendations update in real-time!
```

---

##  Results Visualization

After running analysis (or loading previous results):

### Results Dashboard

```
┌─────────────────────────────────────────────────────────┐
│  ANALYSIS RESULTS                                     │
├─────────────────────────────────────────────────────────┤
│                                                          │
│ Analysis: Salmon-edgeR                                  │
│ Completed: 2025-11-19 14:30:15                         │
│ Runtime: 22 minutes 14 seconds                         │
│                                                          │
│  1,247 genes differentially expressed                │
│    • 623 up-regulated                                   │
│    • 624 down-regulated                                 │
│                                                          │
│ [Summary] [Plots] [Table] [Export] [Report]           │
│                                                          │
└─────────────────────────────────────────────────────────┘
```

### Interactive Plots

#### Volcano Plot
```
Interactive Volcano Plot
(Zoom, pan, hover for gene names)

-log10(FDR)
    ↑
8   │     ●
6   │   ●●●●●
4   │  ●●●●●●●
2   │●●●●●●●●●●
    └───────────────→ log2(Fold Change)
   -4  -2   0   2   4

🔴 Up-regulated (623)
🔵 Down-regulated (624)
⚫ Not significant

[Download Plot] [Customize] [Gene Labels]
```

#### MA Plot
```
MA Plot
(Click genes for details)

log2(FC)
    ↑
4   │  ●●●
2   │●●●●●●
0   │━━━━━━━
-2  │●●●●●●
-4  │  ●●●
    └───────────────→ log10(Mean Expression)
    0   2   4   6

[Download] [Export Gene List]
```

#### Heatmap
```
Top 50 Differentially Expressed Genes

              Control           Treatment
           S1  S2  S3  S4   S5  S6  S7  S8
Gene1    [██ ██ ██ ██] [░░ ░░ ░░ ░░] -3.2
Gene2    [██ ██ ██ ██] [░░ ░░ ░░ ░░] -2.8
...
Gene49   [░░ ░░ ░░ ░░] [██ ██ ██ ██] +2.9
Gene50   [░░ ░░ ░░ ░░] [██ ██ ██ ██] +3.1

[Customize] [Cluster] [Download]
```

### Interactive Data Table

```
Differential Expression Results
(Search, filter, sort - all interactive)

Search: [_______________] 

Gene ID          log2FC  FDR        Mean Expr  Status
──────────────────────────────────────────────────────
ENSG00000111640  3.45    1.2e-08    1250       ● Up
ENSG00000087086  3.21    2.3e-07    980        ● Up
ENSG00000148584  -3.12   1.5e-07    1450       ● Down
ENSG00000183878  -2.98   3.2e-06    820        ● Down
...

Showing 1-20 of 1,247 results

[Export CSV] [Export Excel] [Copy to Clipboard]
```

---

##  Advanced Features

### 1. Multi-Dataset Comparison

```
Compare Multiple Experiments

┌─────────────────────────────────────┐
│ Loaded Datasets:                    │
├─────────────────────────────────────┤
│ ✓ Experiment_A (12 samples)         │
│ ✓ Experiment_B (24 samples)         │
│ ✓ Experiment_C (18 samples)         │
├─────────────────────────────────────┤
│ [Add Dataset] [Compare] [Meta-Analy]│
└─────────────────────────────────────┘

Venn Diagram: Overlapping DE Genes
      
       ┌─────────┐
       │   A     │ 234
    ┌──┼─────┬───┼──┐
    │  │ 123 │89 │  │
    │B └─────┴───┘ C│
    │   167    145  │
    └───────────────┘

[Detailed Overlap] [Export Lists]
```

### 2. Real-Time Resource Monitoring

```
 Live Resource Monitor

CPU Usage:  [████████░░] 82%
Memory:     [██████░░░░] 58% (7.2 / 12 GB)
Disk I/O:   [███░░░░░░░] 34 MB/s

Time Elapsed: 00:15:23
Est. Remaining: 00:06:37

Current Step: Differential Expression Testing
Progress: [██████████░░░░] 73%

[Detailed View] [History] [Alerts]
```

### 3. Ensemble Analysis

```
 Ensemble Analysis Mode

Run multiple pipelines and combine results:

Selected Pipelines:
☑ Pipeline 1: STAR-RSEM-DESeq2
☑ Pipeline 3: Salmon-edgeR
☑ Pipeline 4: Kallisto-Sleuth

Ensemble Method:
◉ Weighted Average (by accuracy)
○ Majority Vote
○ Conservative (intersection)
○ Liberal (union)

[Run Ensemble Analysis →]

Results will show:
• Consensus DE genes
• Pipeline agreement scores
• Confidence per gene
• Robust findings
```

### 4. Parameter Optimization

```
 Automated Parameter Tuning

Optimize parameters for your data:

Pipeline: [Salmon-edgeR ▼]

Parameters to Optimize:
☑ FDR Threshold [0.01 - 0.10]
☑ Log2FC Threshold [0.5 - 2.0]
☑ Min Count Filter [1 - 10]

Optimization Method:
◉ Grid Search (thorough)
○ Bayesian Optimization (fast)

[Start Optimization →]

Est. Time: ~2 hours
```

### 5. Cloud Deployment

```
☁️ Deploy to Cloud

Platform: [AWS ▼]
Region: [us-east-1 ▼]

Instance Type: [r5.4xlarge ▼]
├─ CPU: 16 cores
├─ Memory: 128 GB
└─ Cost: ~$1.20/hour

☑ Use Spot Instances (70% savings)
☑ Auto-shutdown when complete
☑ Email notification

Estimated Cost: $8-12 for this analysis

[Configure] [Deploy →]
```

---

##  Multi-User Setup

### Team Deployment

**1. Server Setup:**
```bash
# Deploy on server
raptor dashboard \
  --host 0.0.0.0 \
  --port 8501 \
  --auth-required \
  --multi-user

# Configure authentication
raptor dashboard config \
  --add-user john@lab.edu --role admin \
  --add-user mary@lab.edu --role analyst \
  --add-user guest@lab.edu --role viewer
```

**2. User Roles:**

| Role | Upload Data | Run Analysis | View Results | Admin |
|------|-------------|--------------|--------------|-------|
| Admin | ✅ | ✅ | ✅ | ✅ |
| Analyst | ✅ | ✅ | ✅ | ❌ |
| Viewer | ❌ | ❌ | ✅ | ❌ |

**3. Collaboration Features:**

```
Share Analysis

Analysis: experiment_2025_11_19

Share with:
☑ john@lab.edu (Can edit)
☑ mary@lab.edu (Can view)

Generate shareable link:
https://raptor.lab.edu/share/abc123xyz

[Copy Link] [Send Email]
```

---

##  Customization

### Themes

```
 Dashboard Settings

Appearance:
◉ Light Mode 
○ Dark Mode
○ Auto (system preference)

Color Scheme:
[Scientific Blue ▼]
• Scientific Blue (default)
• Forest Green
• Sunset Orange
• Monochrome
• Custom...

[Apply] [Reset to Default]
```

### Custom Plots

```
➕ Add Custom Visualization

Plot Type: [Scatter Plot ▼]

X-axis: [Mean Expression ▼]
Y-axis: [Log2 Fold Change ▼]
Color by: [FDR Category ▼]
Size by: [Base Mean ▼]

Filters:
FDR < [0.05]
|log2FC| > [1.0]

[Preview] [Add to Dashboard]
```

### Export Templates

```
📄 Report Templates

Built-in Templates:
• Standard Report (PDF/HTML)
• Publication Supplement
• Lab Notebook Format
• Custom Template...

Custom Template Editor:
[Load Template] [Edit] [Preview] [Save]

Include:
☑ Methods section
☑ All plots
☑ Summary statistics
☑ Gene lists
☑ QC metrics
☐ Raw data tables

[Generate Report →]
```

---

##  Troubleshooting

### Common Issues

#### Issue: Dashboard won't start

**Error:**
```
OSError: [Errno 98] Address already in use
```

**Solution:**
```bash
# Use different port
raptor dashboard --port 8502

# Or kill existing process
lsof -ti:8501 | xargs kill -9
raptor dashboard
```

#### Issue: Upload fails

**Error:**
```
File validation failed: Invalid format
```

**Solution:**
- Check file is CSV/TSV
- Ensure UTF-8 encoding
- No special characters in column names
- Save as CSV from Excel before uploading

#### Issue: Plots not showing

**Error:**
```
Plot rendering failed
```

**Solution:**
```bash
# Clear browser cache
# Or try different browser

# Update dashboard
pip install --upgrade raptor-rnaseq

# Reinstall plotly
pip install --upgrade plotly
```

#### Issue: Slow performance

**Solution:**
```bash
# Limit data displayed
raptor dashboard --max-samples 50

# Use summary mode
raptor dashboard --lightweight

# Increase memory
raptor dashboard --server.maxMemory 4000
```

### Debug Mode

```bash
# Enable debug logging
raptor dashboard --debug

# Check logs
tail -f ~/.raptor/dashboard.log

# Get system info
raptor dashboard --sysinfo
```

---

##  Advanced Usage

### API Integration

Access dashboard programmatically:

```python
from raptor.dashboard import DashboardAPI

# Start dashboard
api = DashboardAPI()
api.start(port=8501)

# Upload data programmatically
api.upload_data(
    counts='data/counts.csv',
    metadata='data/metadata.csv'
)

# Get recommendations
recommendations = api.get_recommendations()

# Run analysis
results = api.run_pipeline(pipeline_id=3)

# Stop dashboard
api.stop()
```

### Custom Extensions

```python
# Add custom visualization
from raptor.dashboard import add_custom_plot

@add_custom_plot
def my_custom_plot(data):
    import plotly.graph_objects as go
    
    fig = go.Figure(data=go.Scatter(
        x=data['x'],
        y=data['y'],
        mode='markers'
    ))
    
    return fig

# Appears in dashboard automatically!
```

---

##  Best Practices

### For Lab Managers

1. **Set up team dashboard** on shared server
2. **Configure user roles** appropriately
3. **Regular backups** of analyses
4. **Document standard workflows**
5. **Train users** on dashboard features

### For Analysts

1. **Always validate uploaded data**
2. **Review QC plots** before proceeding
3. **Understand ML confidence** scores
4. **Compare multiple pipelines** when uncertain
5. **Export and save** analyses regularly

### For Collaborators

1. **Share analyses** via links, not files
2. **Use comments** to document decisions
3. **Export reports** for publications
4. **Keep analysis history**

---

##  Learning Resources

### Video Tutorials

Coming soon:
- Getting Started (5 min)
- Advanced Features (15 min)
- Team Setup (10 min)

### Interactive Demo

```bash
# Launch demo with sample data
raptor dashboard --demo
```

### Documentation

- [User Guide](USER_GUIDE.md)
- [API Reference](API.md)
- [FAQ](FAQ.md)

---

##  Support

**Dashboard not working?**

1. Check [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
2. Read [FAQ](FAQ.md)
3. GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
4. Email: ayehbolouki1988@gmail.com

---

##  Summary

The RAPTOR Dashboard provides:
- ✅ **No-code interface** for RNA-seq analysis
- ✅ **Interactive visualizations** for data exploration
- ✅ **ML-powered recommendations** with explanations
- ✅ **Real-time monitoring** of analyses
- ✅ **Team collaboration** features
- ✅ **Export capabilities** for publications
- ✅ **Cloud integration** for large projects

**Perfect for researchers who prefer visual, interactive workflows!** 🚀

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.0  
**License:** MIT

---

*"Making RNA-seq accessible to everyone, one click at a time!"* 
