# RAPTOR Interactive Dashboard 🦖

Web-based interface for RAPTOR v2.1.0 ML features and analysis tools.

---

##  Quick Start

### Installation

```bash
# Install dashboard dependencies
pip install streamlit plotly psutil

# Or install all RAPTOR ML dependencies
pip install -r requirements_ml.txt
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

##  Features

### 1.  Home
- Welcome page with quick start guide
- System status overview
- Recent activity tracking
- Quick links to main features

### 2.  ML Recommender
**Get AI-powered pipeline recommendations**

- **Upload Data**: CSV count matrix or generate sample data
- **Profile**: Automatic data profiling with key metrics
- **Recommend**: ML model suggests optimal pipeline
- **Visualize**: Confidence scores, feature importance, alternatives

**Example workflow:**
```
1. Upload counts.csv (genes × samples)
2. Click "Profile Data" 
3. Click "Get ML Recommendation"
4. View recommended pipeline with confidence score
```

### 3.  Resource Monitor
**Track system resources in real-time**

- **CPU Usage**: Real-time CPU percentage
- **Memory**: RAM usage and available memory
- **Disk I/O**: Read/write operations
- **Charts**: Interactive time-series plots
- **Export**: Download monitoring data as CSV

**Controls:**
- ▶️ Start: Begin monitoring
- ⏸️ Pause: Pause data collection
- 🔄 Reset: Clear collected data
- 💾 Export: Download CSV

### 4.  Ensemble Analysis
**Combine results from multiple pipelines**

- **Select Pipelines**: Choose 3-6 pipelines to combine
- **Ensemble Methods**:
  - Vote counting (majority vote)
  - Rank product (combine rankings)
  - P-value combination (Fisher's method)
  - Weighted average (by accuracy)
  - Combined (all methods)
- **Results**:
  - Consensus score distribution
  - Pipeline agreement heatmap
  - Top consensus genes
  - Export to CSV/TXT

### 5.  Benchmarks
**Compare pipeline performance**

- Performance metrics (accuracy, precision, recall, F1)
- Runtime comparisons
- Trade-off analysis (accuracy vs speed)
- Detailed comparison tables

### 6.  Settings
**Configure dashboard preferences**

- ML model selection and paths
- Data and output directories
- Performance settings (threads, memory)
- Dashboard preferences (theme, auto-refresh)
- Save/load settings

---

##  File Structure

```
dashboard/
├── dashboard.py           # Main dashboard application
├── README.md             # This file
└── config/
    └── dashboard_settings.json  # Saved settings
```

---

##  Prerequisites

### Required Files

Before using the dashboard, ensure you have:

1. **ML Model** (for recommendations):
   ```bash
   # Train a model using example workflow
   python example_ml_workflow.py --n-datasets 200
   ```
   This creates: `models/raptor_rf_model.pkl`

2. **RAPTOR Modules** in Python path:
   - `ml_recommender.py`
   - `synthetic_benchmarks.py` (optional)

### Python Packages

- streamlit >= 1.28.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- plotly >= 5.17.0
- psutil >= 5.9.0 (for resource monitoring)

---

##  Usage Examples

### Example 1: Get Pipeline Recommendation

```bash
1. Launch dashboard: streamlit run dashboard.py
2. Go to "ML Recommender" page
3. Click "Generate Sample Data" (or upload your CSV)
4. Click "Profile Data"
5. Click "Get ML Recommendation"
6. View recommended pipeline and confidence score
```

### Example 2: Monitor Resource Usage

```bash
1. Go to "Resource Monitor" page
2. Click "▶️ Start Monitoring"
3. Run your analysis in another terminal
4. Watch real-time CPU/Memory usage
5. Click "💾 Export Data" to save metrics
```

### Example 3: Ensemble Analysis

```bash
1. Go to "Ensemble Analysis" page
2. Select multiple pipelines (e.g., 1, 3, 5)
3. Choose ensemble method (e.g., "p_value_combination")
4. Click "Run Ensemble Analysis"
5. View consensus genes and agreement heatmap
6. Download high-confidence genes
```

---

##  Customization

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

### Model Path

In Settings page, change "Model Directory" to your model location.

---

##  Troubleshooting

### Problem: "Model not found"

**Solution:**
```bash
# Train a model first
python example_ml_workflow.py --n-datasets 200

# Or specify correct model path in Settings
```

### Problem: "psutil not installed"

**Solution:**
```bash
pip install psutil
```

### Problem: "Port 8501 already in use"

**Solution:**
```bash
# Use different port
streamlit run dashboard.py --server.port 8502

# Or kill existing process
# On Linux/Mac:
lsof -ti:8501 | xargs kill -9

# On Windows:
netstat -ano | findstr :8501
taskkill /PID <PID> /F
```

### Problem: "Import error: ml_recommender"

**Solution:**
```bash
# Ensure files are in Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Or run from RAPTOR root directory
cd /path/to/RAPTOR
streamlit run dashboard/dashboard.py
```

### Problem: Dashboard is slow

**Solution:**
- Reduce monitoring frequency in Resource Monitor
- Clear browser cache
- Use smaller datasets
- Increase memory limit in Settings

---

##  Security Notes

### For Production Deployment

1. **Change default settings** in dashboard_settings.json
2. **Use authentication** if deploying publicly
3. **Restrict access** to specific IPs
4. **Enable HTTPS** for remote access
5. **Set memory limits** to prevent resource exhaustion

### Running on Server

```bash
# With authentication
streamlit run dashboard.py --server.headless true

# With password protection (requires config)
# Edit ~/.streamlit/credentials.toml
```

---

##  Documentation

- **RAPTOR Main Docs**: See `docs/` folder
- **ML Guide**: `docs/ML_GUIDE.md`
- **Ensemble Guide**: `docs/ENSEMBLE_GUIDE.md`
- **Troubleshooting**: `docs/TROUBLESHOOTING.md`

---

##  Contributing

Found a bug or have a feature request?

1. Check existing issues: https://github.com/AyehBlk/RAPTOR/issues
2. Open a new issue with:
   - Dashboard version
   - Browser and OS
   - Steps to reproduce
   - Screenshots (if applicable)

---

##  Performance Tips

### For Large Datasets

- Use data subsampling for profiling
- Enable caching in Settings
- Increase memory limit
- Use batch processing mode (if available)

### For Slow Connections

- Reduce auto-refresh frequency
- Disable real-time monitoring when not needed
- Use CSV export instead of interactive plots

---

##  Tutorial

### First-Time Users

**Step 1: Check System**
- Go to Home page
- Verify all status indicators are green ✅
- If model not found, follow training instructions

**Step 2: Try Sample Data**
- Go to ML Recommender
- Click "Use sample data"
- Click "Generate Sample Data"
- This creates 1000 genes × 6 samples

**Step 3: Get Recommendation**
- Click "Profile Data"
- Wait for profiling (5-10 seconds)
- Click "Get ML Recommendation"
- View confidence score and reasons

**Step 4: Explore Features**
- Try Resource Monitor (watch live metrics)
- Try Ensemble Analysis (combine pipelines)
- Check Benchmarks (compare pipelines)

---

##  Citation

If you use RAPTOR Dashboard in your research:

```bibtex
@software{raptor_dashboard_2025,
  author = {Bolouki, Ayeh},
  title = {RAPTOR Interactive Dashboard: Web Interface for 
           ML-Powered RNA-seq Analysis},
  year = {2025},
  version = {2.1.0},
  url = {https://github.com/AyehBlk/RAPTOR}
}
```

---

##  Support

- **Email**: ayehbolouki1988@gmail.com
- **GitHub Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Documentation**: See `docs/` folder

---

##  License

MIT License - see LICENSE file

---

##  Acknowledgments

Built with:
- [Streamlit](https://streamlit.io/) - Web framework
- [Plotly](https://plotly.com/) - Interactive visualizations
- [Pandas](https://pandas.pydata.org/) - Data manipulation

---

**Making free science for everybody around the world** 

---

**Author**: Ayeh Bolouki
**Version**: 2.1.0  
**Last Updated**: November 2025
