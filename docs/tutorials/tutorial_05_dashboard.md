# Tutorial 5: Using the Interactive Dashboard

**Level**: Beginner  
**Time**: 30 minutes  
**Goal**: Master RAPTOR's web-based interactive dashboard

---

## What You'll Learn

- How to launch and navigate the dashboard
- Real-time resource monitoring
- Interactive data exploration
- Parameter tuning without re-running analysis
- Generating reports from the dashboard

---

## Prerequisites

- RAPTOR v2.1.0+ installed
- Basic web browser skills
- Completed Tutorial 1 (recommended)

---

## Why Use the Dashboard?

### Command-Line (Traditional)
```
Run analysis → Wait → See results → Want to try different parameters → Run again → Wait...
```
 **Time-consuming and inefficient**

### Dashboard (v2.1.0)
```
Run analysis ONCE → Explore interactively → Try different thresholds → See results instantly
```
 **Fast, interactive, visual**

---

## Step 1: Launch the Dashboard

### Quick Start

```bash
# Start dashboard in your results directory
cd my_raptor_project/
raptor dashboard --results results/

# Dashboard starts on http://localhost:8501
```

**You'll see:**
```
RAPTOR Interactive Dashboard v2.1.0
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

✓ Loading results from: results/
✓ Dashboard ready!

→ Open your browser to: http://localhost:8501

Press Ctrl+C to stop the dashboard
```

### Custom Port

```bash
# Use different port if 8501 is busy
raptor dashboard --results results/ --port 8080
```

### From Python

```python
from raptor import DashboardServer

dashboard = DashboardServer(
    results_dir='results/',
    port=8501,
    auto_refresh=True,  # Update every 5 seconds
    theme='dark'  # or 'light'
)

dashboard.start()
# Open browser to http://localhost:8501

# When done
dashboard.stop()
```

---

## Step 2: Dashboard Overview

When you open the dashboard, you'll see:

```
┌─────────────────────────────────────────────────────────────┐
│ RAPTOR Interactive Dashboard            [Results loaded ✓] │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│   Overview  │   Recommendations  │   Explore Data  │ │
│   Monitor   │   Parameters       │   Reports      │ │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Six Main Tabs:

1. ** Overview** - Summary of your analysis
2. ** Recommendations** - Pipeline suggestions
3. ** Explore Data** - Interactive plots
4. ** Monitor** - Resource usage (live!)
5. ** Parameters** - Tune thresholds
6. ** Reports** - Generate documents

---

## Step 3: Overview Tab

The Overview tab shows your analysis summary:

```
┌─────────────────────────────────────────────┐
│  Dataset Information                         │
├─────────────────────────────────────────────┤
│   Samples: 12 (6 control, 6 treatment)    │
│   Genes: 18,234 detected                  │
│   Quality Score: 87/100 🟢                │
│  ⚠️  Issues: 1 minor (batch effect)         │
└─────────────────────────────────────────────┘

┌─────────────────────────────────────────────┐
│  Key Metrics                                 │
├─────────────────────────────────────────────┤
│  BCV: 0.42 (Medium variation)               │
│  Depth: 25.3M reads/sample (High)           │
│  Library CV: 0.12 (Good)                    │
│  Zero inflation: 42% (Normal)               │
└─────────────────────────────────────────────┘

Quick Actions:
[View Recommendations] [Quality Report] [Resource Monitor]
```

**Interactive elements:**
- Click on metrics to see detailed plots
- Hover over quality score for breakdown
- Click "View Recommendations" to jump to that tab

---

## Step 4: Recommendations Tab

See and compare pipeline recommendations:

```
┌────────────────────────────────────────────────────────┐
│   #1: Salmon-edgeR                                   │
├────────────────────────────────────────────────────────┤
│  ML Score: 0.89 │ Confidence: HIGH │ Success: 87%     │
│                                                         │
│  Why this pipeline?                                    │
│  ✓ Perfect for your sample size (12)                  │
│  ✓ Handles medium BCV (0.42) excellently              │
│  ✓ Fast runtime (~22 min estimated)                   │
│                                                         │
│   Based on 1,247 similar projects                   │
│  [View Similar] [Feature Importance] [Run Pipeline]   │
└────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────┐
│   #2: STAR-RSEM-DESeq2                              │
├────────────────────────────────────────────────────────┤
│  ML Score: 0.85 │ Confidence: HIGH │ Success: 92%     │
│  [Details ▼]                                           │
└────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────┐
│  Compare All Recommendations                            │
├────────────────────────────────────────────────────────┤
│  [Interactive Comparison Chart]                        │
└────────────────────────────────────────────────────────┘
```

**Interactive features:**
- Click "[View Similar]" to see similar past projects
- Click "[Feature Importance]" for explanation
- Use sliders to adjust resource constraints
- Compare multiple recommendations side-by-side

---

## Step 5: Explore Data Tab

Interactive data visualization - the most powerful tab!

### 5.1: Library Size Distribution

```
┌──────────────────────────────────────────────┐
│  Library Sizes                               │
├──────────────────────────────────────────────┤
│  [Interactive histogram]                     │
│                                               │
│  Options:                                    │
│  □ Log scale                                 │
│  ☑ Show mean/median                         │
│  □ Color by condition                        │
│  [Download Plot] [Download Data]            │
└──────────────────────────────────────────────┘
```

**Try this:**
1. Check "Color by condition" - Do groups have similar library sizes?
2. Hover over bars to see exact values
3. Click "[Download Plot]" to save as PNG

### 5.2: PCA Plot (Interactive!)

```
┌──────────────────────────────────────────────┐
│  Principal Component Analysis                │
├──────────────────────────────────────────────┤
│  [Interactive 2D scatter plot]               │
│                                               │
│  Color by: [Condition ▼] [Batch] [Replicate]│
│  PC X-axis: [PC1 (68%) ▼]                   │
│  PC Y-axis: [PC2 (15%) ▼]                   │
│                                               │
│   Click points to identify samples        │
│   Zoom with mouse wheel                   │
│  [3D View] [Download]                        │
└──────────────────────────────────────────────┘
```

**Try this:**
1. Color by "Batch" - Any batch effects?
2. Switch to PC2 vs PC3 - Different patterns?
3. Click on outlier points - Identify problem samples
4. Try "[3D View]" to see PC1, PC2, PC3 together

### 5.3: Gene Expression Heatmap

```
┌──────────────────────────────────────────────┐
│  Expression Heatmap                          │
├──────────────────────────────────────────────┤
│  Show: [Top 50 ▼] variable genes            │
│                                               │
│  [Interactive heatmap with zoom]            │
│                                               │
│  Clustering:                                 │
│  ☑ Cluster samples                          │
│  ☑ Cluster genes                            │
│  □ Show dendrogram                          │
│                                               │
│   Color scale: [Red-Blue ▼]              │
│  [Download]                                  │
└──────────────────────────────────────────────┘
```

**Try this:**
1. Increase to "Top 100" genes
2. Check "Show dendrogram" to see relationships
3. Click on gene names to see expression across samples
4. Hover over cells for exact values

### 5.4: Correlation Matrix

```
┌──────────────────────────────────────────────┐
│  Sample Correlation                          │
├──────────────────────────────────────────────┤
│  [Interactive correlation heatmap]          │
│                                               │
│  Correlation: [Pearson ▼] [Spearman]       │
│                                               │
│  Options:                                    │
│  ☑ Show values                              │
│  ☑ Cluster samples                          │
│  □ Color by condition                        │
│  [Download]                                  │
└──────────────────────────────────────────────┘
```

**Look for:**
- High within-group correlation (>0.85) ✅ Good!
- Low between-group correlation (<0.70) ✅ Good!
- Outlier samples (low correlation with everyone)

---

## Step 6: Monitor Tab (Live!)

Real-time resource monitoring - watch your analysis run:

```
┌──────────────────────────────────────────────────┐
│  Resource Monitor                    [LIVE ●]   │
├──────────────────────────────────────────────────┤
│                                                   │
│  CPU Usage (updating every 5s)                  │
│  [Animated line graph]                          │
│  Current: 82% │ Peak: 95% │ Avg: 78%          │
│                                                   │
│  Memory Usage                                    │
│  [Animated line graph]                          │
│  Current: 18.6 GB │ Peak: 22.1 GB │ Max: 32 GB│
│  [████████████████░░░░░░] 69%                   │
│                                                   │
│  Disk I/O                                        │
│  Read:  156 MB/s │ Write: 34 MB/s              │
│  [Mini graph]                                    │
│                                                   │
│  Status: Salmon quantification (Sample 8/12)   │
│  Estimated time remaining: 23 minutes           │
│                                                   │
│  ⚠️  Alert: Memory >90%                         │
│  [Pause Analysis] [Adjust Resources]           │
└──────────────────────────────────────────────────┘
```

**Features:**
- **Live updates** - See resources in real-time
- **Alerts** - Get notified of issues
- **Estimates** - Know when analysis will finish
- **Historical view** - See resource usage over time

**Try during analysis:**
1. Watch CPU/memory graphs update
2. Check estimated completion time
3. Set alerts (e.g., "Alert if memory >90%")
4. Download resource usage data

---

## Step 7: Parameters Tab

Explore different parameter settings WITHOUT re-running!

```
┌──────────────────────────────────────────────┐
│  Parameter Explorer                           │
├──────────────────────────────────────────────┤
│                                               │
│  FDR Threshold: [━━━●━━━] 0.05               │
│                 0.01      0.10                │
│                                               │
│  Log2 Fold Change: [━━━━●━] 1.0             │
│                    0.5       2.0              │
│                                               │
│  Results with current settings:              │
│   DE Genes: 1,789                          │
│   Up-regulated: 924                        │
│   Down-regulated: 865                      │
│                                               │
│  [Reset to Defaults] [Apply & Save]         │
└──────────────────────────────────────────────┘

┌──────────────────────────────────────────────┐
│  Volcano Plot (Interactive)                  │
│  [Plot updates as you move sliders!]        │
│                                               │
│  🔴 Red dots: Significant genes              │
│  ⚪ Gray dots: Not significant               │
│  [Click dots to see gene names]             │
└──────────────────────────────────────────────┘
```

**Try this exercise:**

1. **Start conservative:**
   - FDR: 0.01
   - Log2FC: 1.5
   - **Result:** 543 genes

2. **Go more liberal:**
   - FDR: 0.10
   - Log2FC: 0.5
   - **Result:** 2,789 genes

3. **Find balance:**
   - FDR: 0.05
   - Log2FC: 1.0
   - **Result:** 1,789 genes ✓

**Watch the volcano plot update in real-time!**

### Parameter Sensitivity

```
┌──────────────────────────────────────────────┐
│  Parameter Sensitivity Analysis              │
├──────────────────────────────────────────────┤
│  How sensitive are results to changes?      │
│                                               │
│  FDR Threshold:         ██████ HIGH          │
│  Log2FC Threshold:      ████ MEDIUM         │
│  Min Count Filter:      ▌ LOW               │
│                                               │
│  Recommendation:                          │
│  Be careful with FDR - small changes cause  │
│  big differences in gene count.             │
└──────────────────────────────────────────────┘
```

---

## Step 8: Reports Tab

Generate publication-ready reports from the dashboard:

```
┌──────────────────────────────────────────────┐
│  Generate Reports                             │
├──────────────────────────────────────────────┤
│                                               │
│  Report Type:                                │
│  ○ Profile Report (Recommendations)          │
│  ● Analysis Report (Results + QC)           │
│  ○ Comparison Report (Multiple pipelines)   │
│  ○ Methods Section (For papers)             │
│                                               │
│  Include:                                    │
│  ☑ Executive summary                        │
│  ☑ Quality assessment                       │
│  ☑ Data visualizations                      │
│  ☑ Statistical details                      │
│  ☑ Parameter settings                       │
│  □ Raw data tables                          │
│                                               │
│  Format: [HTML ▼] [PDF] [Markdown]         │
│                                               │
│  [Generate Report]                           │
└──────────────────────────────────────────────┘
```

**Generated report includes:**
- All plots from Explore tab
- Current parameter settings
- Quality metrics
- Resource usage
- Recommendations
- Methods text (ready to paste into manuscript!)

---

## Step 9: Real-Time Collaboration

Share your dashboard with colleagues:

### Local Network

```bash
# Start dashboard with network access
raptor dashboard \
  --results results/ \
  --host 0.0.0.0 \  # Allow external access
  --port 8501

# Share this URL with colleagues:
# http://YOUR_IP:8501
```

**Use case:** Show results to PI during lab meeting

### Remote Server

```bash
# On server
raptor dashboard --results results/ --port 8501

# On your computer
ssh -L 8501:localhost:8501 user@server

# Open browser to: http://localhost:8501
```

**Use case:** Analyze data on HPC, view dashboard on laptop

---

## Step 10: Dashboard Keyboard Shortcuts

Make navigation faster:

```
Global Shortcuts:
├─ Tab        : Next tab
├─ Shift+Tab  : Previous tab
├─ Ctrl+R     : Refresh data
├─ Ctrl+D     : Download current view
├─ Ctrl+P     : Print/Save as PDF
└─ Ctrl+Q     : Quit dashboard

In Plots:
├─ Click+Drag : Pan
├─ Scroll     : Zoom
├─ Double-Click : Reset view
├─ Ctrl+Click : Select multiple
└─ Shift+Click : Box select
```

---

## Advanced Features

### 1. Compare Multiple Analyses

```python
from raptor import DashboardServer

# Load multiple result sets
dashboard = DashboardServer(
    results_dirs=[
        'experiment1/results/',
        'experiment2/results/',
        'experiment3/results/'
    ],
    comparison_mode=True
)

dashboard.start()
```

**Dashboard now shows:**
- Side-by-side comparisons
- Shared genes across experiments
- Batch effect analysis across experiments

### 2. Custom Plots

```python
import plotly.graph_objects as go
from raptor import DashboardServer

# Create custom plot
fig = go.Figure()
fig.add_trace(go.Scatter(x=x_data, y=y_data))

# Add to dashboard
dashboard = DashboardServer(results_dir='results/')
dashboard.add_custom_plot('My Custom Plot', fig)
dashboard.start()
```

### 3. Export Session

```bash
# Save current dashboard state
Ctrl+S in dashboard

# Or from command line
raptor dashboard \
  --results results/ \
  --export-session my_session.pkl

# Reload later
raptor dashboard \
  --load-session my_session.pkl
```

---

## Troubleshooting

### Problem: Dashboard won't start

```bash
# Check if port is in use
lsof -i :8501

# Kill existing process
kill -9 <PID>

# Or use different port
raptor dashboard --results results/ --port 8080
```

### Problem: Dashboard is slow

**Solutions:**
```bash
# 1. Reduce auto-refresh frequency
raptor dashboard \
  --results results/ \
  --refresh-interval 30  # 30 seconds instead of 5

# 2. Disable live monitoring
raptor dashboard \
  --results results/ \
  --no-live-monitor

# 3. Use lighter theme
raptor dashboard \
  --results results/ \
  --theme light \
  --minimal-plots
```

### Problem: Can't see plots

**Check browser:**
- Use Chrome, Firefox, or Safari (not IE)
- Enable JavaScript
- Clear browser cache
- Try incognito mode

---

## Best Practices

### Do's:
✅ Use dashboard for interactive exploration  
✅ Save interesting parameter combinations  
✅ Share dashboard with collaborators  
✅ Generate reports directly from dashboard  
✅ Monitor resources during long analyses  

### Don'ts:
❌ Leave dashboard running indefinitely (wastes resources)  
❌ Share dashboard on public network (security risk)  
❌ Rely solely on dashboard (keep command-line skills)  
❌ Modify files while dashboard is running  
❌ Try to run multiple analyses simultaneously via dashboard  

---

## Summary

You've learned to:
- ✅ Launch and navigate the interactive dashboard
- ✅ Explore data with interactive plots
- ✅ Monitor resources in real-time
- ✅ Tune parameters without re-running
- ✅ Generate reports from the dashboard
- ✅ Share results with collaborators
- ✅ Use advanced dashboard features

---

## Next Steps

1. **Try Tutorial 6**: [Ensemble Analysis](tutorial_06_ensemble.md)
2. **Try Tutorial 7**: [Resource Optimization](tutorial_07_resources.md)
3. **Read**: [DASHBOARD.md](../DASHBOARD.md) - Complete dashboard guide

---

**Tutorial by Ayeh Bolouki**  
For RAPTOR v2.1.0

*"See your data come alive!"* 
