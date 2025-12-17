# RAPTOR Interactive Dashboard 🦖

Web-based interface for RAPTOR v2.1.1 ML features and analysis tools.

---

## 🆕 What's New in v2.1.1

### Adaptive Threshold Optimizer (ATO)

The new **Threshold Optimizer** page provides data-driven threshold selection for differential expression analysis:

- **Data-driven logFC cutoffs** instead of arbitrary |logFC| > 1
- **Multiple p-value adjustment methods**: BH, BY, Holm, Hochberg, Bonferroni, q-value
- **π₀ estimation** for true null proportion
- **Interactive visualizations**: Volcano plots, distributions, heatmaps
- **Publication-ready methods text** generation
- **Export optimized results** as CSV

---

## ⚡ Quick Start

### Installation

```bash
# Install from PyPI (recommended)
pip install raptor-rnaseq[dashboard]

# Or install all features
pip install raptor-rnaseq[all]
```

### Launch Dashboard

```bash
# Start dashboard on default port (8501)
streamlit run dashboard.py

# Or specify custom port
streamlit run dashboard.py --server.port 8502

# Access on local network
streamlit run dashboard.py --server.address 0.0.0.0
```

The dashboard will open automatically in your browser at `http://localhost:8501`

---

## 📋 Features

### 1. 🏠 Home
- Welcome page with quick start guide
- System status overview (including ATO availability)
- Recent activity tracking
- Quick links to main features

### 2. 🤖 ML Recommender
**Get AI-powered pipeline recommendations**

- **Upload Data**: CSV count matrix or generate sample data
- **Profile**: Automatic data profiling with key metrics
- **Recommend**: ML model suggests optimal pipeline
- **Visualize**: Confidence scores, feature importance, alternatives

### 3. 🎯 Threshold Optimizer *(NEW in v2.1.1)*
**Data-driven threshold optimization for DE analysis**

- **Upload DE Results**: CSV/TSV from DESeq2, edgeR, or limma
- **Set Analysis Goal**: Discovery, Balanced, or Validation
- **Optimize**: Get recommended logFC and padj thresholds
- **Visualize**: Volcano plots, distributions, comparison heatmaps
- **Export**: Significant genes, summary report, methods text

**Supported input formats:**
- DESeq2 output (log2FoldChange, pvalue, padj)
- edgeR output (logFC, PValue, FDR)
- limma-voom output (logFC, P.Value, adj.P.Val)

**Example workflow:**
```
1. Upload your DESeq2 results (CSV)
2. Select goal: "discovery" 
3. Click "Optimize Thresholds"
4. View recommended: |logFC| > 0.45, padj < 0.05
5. Download significant genes and methods text
```

### 4. 📊 Resource Monitor
**Track system resources in real-time**

- **CPU Usage**: Real-time CPU percentage
- **Memory**: RAM usage and available memory
- **Disk I/O**: Read/write operations
- **Charts**: Interactive time-series plots
- **Export**: Download monitoring data as CSV

### 5. 🔬 Ensemble Analysis
**Combine results from multiple pipelines**

- **Select Pipelines**: Choose 3-6 pipelines to combine
- **Ensemble Methods**: Vote, rank product, p-value combination, weighted
- **Results**: Consensus scores, agreement heatmap, top genes
- **Export**: CSV/TXT files

### 6. 📈 Benchmarks
**Compare pipeline performance**

- Performance metrics (accuracy, precision, recall, F1)
- Runtime comparisons
- Trade-off analysis (accuracy vs speed)
- Detailed comparison tables

### 7. ⚙️ Settings
**Configure dashboard preferences**

- ML model selection and paths
- Data and output directories
- Performance settings (threads, memory)
- **Threshold Optimizer defaults** *(NEW)*
- Dashboard preferences (theme, auto-refresh)

---

## 📁 File Structure

```
RAPTOR/
├── launch_dashboard.py           # Dashboard launcher (root)
├── dashboard/                    # Dashboard folder
│   ├── dashboard.py              # Main dashboard application
│   ├── dashboard_README.md       # This file
│   └── docs/
│       └── THRESHOLD_OPTIMIZER.md
└── raptor/                       # Python package
    ├── __init__.py
    ├── ml_recommender.py
    ├── threshold_optimizer/      # ATO module (v2.1.1)
    │   ├── __init__.py
    │   ├── ato.py                # Main optimizer class
    │   └── visualization.py      # Plotting functions
    └── ... other modules

```

---

## 🎯 Threshold Optimizer Details

### Input Requirements

| Column | Aliases | Required |
|--------|---------|----------|
| log2FoldChange | logFC, log2FC, lfc | ✅ Yes |
| pvalue | PValue, P.Value, pval | ✅ Yes |
| padj | adj.P.Val, FDR, qvalue | Optional |
| baseMean | AveExpr, logCPM | Optional |
| lfcSE | SE | Optional |

### Analysis Goals

| Goal | Error Control | Use Case |
|------|--------------|----------|
| **Discovery** | FDR | Exploratory analysis, maximize true positives |
| **Balanced** | FDR | Publication, balance sensitivity/specificity |
| **Validation** | FWER | Confirming targets, minimize false positives |

### LogFC Methods

| Method | Description |
|--------|-------------|
| **auto** | Consensus of all methods (recommended) |
| **mad** | MAD-based robust estimation |
| **mixture** | Gaussian mixture model |
| **power** | Power-based minimum effect |
| **percentile** | 95th percentile of null |

### P-value Adjustment Methods

- **BH**: Benjamini-Hochberg (standard FDR)
- **BY**: Benjamini-Yekutieli (any dependence)
- **q-value**: Storey's q-value with π₀ estimation
- **Holm**: Step-down FWER control
- **Hochberg**: Step-up FWER control
- **Bonferroni**: Most conservative

---

## 💡 Usage Examples

### Example 1: Optimize DE Thresholds

```bash
1. Launch dashboard: streamlit run dashboard.py
2. Go to "🎯 Threshold Optimizer" page
3. Upload your DESeq2 results CSV
4. Select goal: "discovery"
5. Click "🚀 Optimize Thresholds"
6. View results:
   - Recommended |logFC| cutoff
   - Recommended padj method and cutoff
   - π₀ estimate
   - Number of DE genes
7. Download significant genes CSV
```

### Example 2: Get Pipeline Recommendation

```bash
1. Go to "🤖 ML Recommender" page
2. Click "Generate Sample Data" (or upload your CSV)
3. Click "Profile Data" 
4. Click "Get ML Recommendation"
5. View recommended pipeline with confidence score
```

### Example 3: Monitor Resources

```bash
1. Go to "📊 Resource Monitor" page
2. Click "▶️ Start Monitoring"
3. Run your analysis in another terminal
4. Watch real-time CPU/Memory usage
5. Click "💾 Export Data" to save metrics
```

---

## 🔧 Customization

### Change Port

```bash
streamlit run dashboard.py --server.port 8888
```

### Change Theme

Edit Streamlit config: `~/.streamlit/config.toml`

```toml
[theme]
primaryColor = "#2E7D32"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
textColor = "#262730"
font = "sans serif"
```

---

## 🐛 Troubleshooting

### Problem: "Threshold Optimizer not available"

**Solution:**
```bash
# Ensure RAPTOR v2.1.1+ is installed
pip install raptor-rnaseq>=2.1.1

# Or copy threshold_optimizer/ folder manually
```

### Problem: "Model not found"

**Solution:**
```bash
# Train a model first
python example_ml_workflow.py --n-datasets 200
```

### Problem: "psutil not installed"

**Solution:**
```bash
pip install raptor-rnaseq[dashboard]
```

### Problem: "Port 8501 already in use"

**Solution:**
```bash
# Use different port
streamlit run dashboard.py --server.port 8502

# Or kill existing process
lsof -ti:8501 | xargs kill -9
```

---

## 📚 Documentation

- **RAPTOR Main Docs**: See `docs/` folder
- **Threshold Optimizer**: `docs/THRESHOLD_OPTIMIZER.md`
- **ML Guide**: `docs/ML_GUIDE.md`
- **Ensemble Guide**: `docs/ENSEMBLE_GUIDE.md`

---

## 📖 Citation

If you use RAPTOR Dashboard in your research:

```bibtex
@software{bolouki2025raptor,
  author = {Bolouki, Ayeh},
  title = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year = {2025},
  version = {2.1.1},
  doi = {10.5281/zenodo.17607161},
  url = {https://github.com/AyehBlk/RAPTOR}
}
```

---

## 🆘 Support

- **PyPI**: https://pypi.org/project/raptor-rnaseq/
- **Email**: ayehbolouki1988@gmail.com
- **GitHub Issues**: https://github.com/AyehBlk/RAPTOR/issues

---

## 📜 License

MIT License - see LICENSE file

---

## 🙏 Acknowledgments

Built with:
- [Streamlit](https://streamlit.io/) - Web framework
- [Plotly](https://plotly.com/) - Interactive visualizations
- [Pandas](https://pandas.pydata.org/) - Data manipulation
- [SciPy](https://scipy.org/) - Statistical functions

---

**Making free science for everybody around the world** 🌍

---

**Author**: Ayeh Bolouki  
**Version**: 2.1.1  
**PyPI**: https://pypi.org/project/raptor-rnaseq/
