"""
RAPTOR Dashboard v2.2.2
Main Entry Point (Using Style Loader Helpers)

This version uses the style_loader.py helper functions for cleaner code.

Author: Ayeh Bolouki
Date: 14 April 2026
"""

import streamlit as st
from pathlib import Path
import sys

# Add current directory to path
sys.path.append(str(Path(__file__).parent))

# Import style helpers
try:
    from assets.styles.style_loader import (
        load_custom_css,
        load_logo,
        add_footer
    )
    STYLES_AVAILABLE = True
except ImportError:
    STYLES_AVAILABLE = False
    # Fallback functions if style_loader not available
    def load_custom_css():
        st.markdown("""<style>
            .main h1 { color: #2E7D32; }
        </style>""", unsafe_allow_html=True)
    
    def load_logo():
        return None
    
    def add_footer():
        st.caption("© 2026 Ayeh Bolouki")

# Page configuration
st.set_page_config(
    page_title="RAPTOR Dashboard",
    page_icon="🦖",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Load custom styles globally
load_custom_css()

# ============================================================================
# SIDEBAR WITH LOGO
# ============================================================================

with st.sidebar:
    # Add logo at top of sidebar
    logo_path = load_logo()
    if logo_path:
        st.image(logo_path, use_container_width=True)
        st.markdown("---")
    else:
        st.markdown("## 🦖 RAPTOR v2.2.2")
        st.markdown("---")
    
    # Navigation info
    st.markdown("### 📋 Navigation")
    st.info("""
    **Dashboard Features:**
    
    📡 **Data Acquisition (Module 6b):**
    • Search GEO, SRA, TCGA, ArrayExpress
    • Download & upload datasets
    • Pool datasets with batch correction
    • Gene ID conversion
    
    📊 **Data Assessment (Modules 2-4):**
    • Quality Assessment - QC & outlier detection
    • Data Profiler - Comprehensive profiling
    • ML Recommender - Pipeline suggestions
    
    🔬 **DE Analysis Workflow (Modules 7-9):**
    • Import DE - Load DESeq2/edgeR/limma results
    • Parameter Optimization - Find optimal thresholds
    • Ensemble Analysis - Combine multiple methods
    
    📈 **Visualization & Output:**
    • Interactive Visualizations
    • Publication-ready Reports
    • Dashboard Settings
    
    💡 **Note:** Run quantification (Modules 1 & 5) via CLI before 
    using dashboard. Start with count matrices from Salmon/Kallisto.
    """)

# ============================================================================
# MAIN CONTENT
# ============================================================================

# Header
st.title("🦖 RAPTOR Dashboard v2.2.2")
st.markdown("**RNA-seq Analysis Pipeline - Professional Interface**")
st.markdown("---")

# Welcome message
st.info("""
👈 **Select a page from the sidebar to begin!**

This dashboard provides a user-friendly interface for RAPTOR modules:
- 📡 Search & download from GEO, SRA, TCGA, ArrayExpress
- 📊 Quality assessment and data profiling
- 📥 Import differential expression results
- 🔬 Run ensemble analysis  
- 📈 Visualize results
- 📋 Export publication-ready reports
""")

# Quick stats
col1, col2, col3 = st.columns(3)

with col1:
    st.metric(
        label="🦖 Version",
        value="2.2.2",
        help="Current RAPTOR version"
    )

with col2:
    st.metric(
        label="📦 Dashboard Modules",
        value="8",
        help="Active modules in dashboard (M2-M4, M6b, M7-M9 + extras)"
    )

with col3:
    st.metric(
        label="✅ Status",
        value="Ready",
        help="Dashboard operational status"
    )

st.markdown("---")

# Getting started guide
with st.expander("📚 Getting Started Guide", expanded=False):
    st.markdown("""
    ### For New Users:
    
    **Starting Point:** You can either search & download datasets from
    public repositories (GEO, SRA, TCGA) using Data Acquisition, or
    start with count matrices from Salmon/Kallisto (run via CLI using
    Modules 1 & 5).
    
    **Dashboard Workflow:**
    
    0. **Data Acquisition** (Module 6b) [Optional]
       - Search GEO, SRA, TCGA, ArrayExpress
       - Download datasets or upload your own
       - Pool multiple datasets with batch correction
       - Gene ID conversion
    
    1. **Quality Assessment** (Module 2)
       - Upload count matrix
       - Detect outliers (6 methods)
       - Check batch effects
       
    2. **Data Profiling** (Module 3)
       - Characterize your data
       - Calculate BCV, sparsity
       - Get 32 statistical features
       
    3. **Get Recommendation** (Module 4)
       - ML-based pipeline suggestion
       - Optimal parameters
       
    4. **Import DE Results** (Module 7)
       - Upload DESeq2, edgeR, or limma results
       - Auto-detection and validation
       
    5. **Optimize Parameters** (Module 8) [Optional]
       - Find optimal FDR/logFC thresholds
       - 4 scientific methods
       
    6. **Run Ensemble Analysis** (Module 9)
       - Combine multiple DE methods
       - Get robust consensus genes
       
    7. **Visualize & Export**
       - Interactive volcano plots, heatmaps
       - PDF reports, CSV gene lists
    
    ### Quick Links:
    - 📖 [Documentation](https://github.com/AyehBlk/RAPTOR)
    - 💬 [Support](https://github.com/AyehBlk/RAPTOR/issues)
    - 📧 [Email](mailto:ayehbolouki1988@gmail.com)
    """)

# Feature cards
st.markdown("### 🎯 Key Features")

col1, col2, col3, col4 = st.columns(4)

with col1:
    with st.container():
        st.markdown("#### 📡 Data Acquisition")
        st.write("""
        Search GEO, SRA, TCGA, and ArrayExpress. 
        Download, pool, and prepare datasets for analysis.
        """)

with col2:
    with st.container():
        st.markdown("#### 📊 Quality Control")
        st.write("""
        Comprehensive QC with 6-method outlier detection, 
        batch effect analysis, and quality scoring.
        """)

with col3:
    with st.container():
        st.markdown("#### 🧬 Ensemble Analysis")
        st.write("""
        Combine DESeq2, edgeR, and limma results using 
        Fisher's method, Brown's method, or RRA.
        """)

with col4:
    with st.container():
        st.markdown("#### ⚙️ Optimization")
        st.write("""
        4 scientific methods to find optimal FDR and 
        log2FC thresholds for your data.
        """)

# Recent activity
st.markdown("---")
st.markdown("### 📋 Recent Activity")
st.info("No recent activity. Start with Quality Assessment or Import DE!")

# Workflow diagram
with st.expander("🔄 Complete Workflow"):
    st.markdown("""
    ### RAPTOR Analysis Workflow
    
    ```
    STAGE 0: Data Acquisition (Dashboard) [Optional]
    ──────────────────────────────────────────────────
    Step 0: Get Data (Module 6b)
       - Search GEO, SRA, TCGA, ArrayExpress
       - Download datasets or upload your own
       - Pool multiple studies with batch correction
       - Gene ID conversion
    
    STAGE 1: Quantification (CLI)
    ──────────────────────────────
    Step 1: Run Salmon/Kallisto
       - Fast quantification via CLI
       - Generates count matrix
    
    STAGE 2: Assessment (Dashboard)
    ──────────────────────────────
    Step 2: Quality Assessment (Module 2)
       - QC metrics
       - Outlier detection
       ↓
    Step 3: Data Profiler (Module 3)
       - Characterize your data
       - Calculate BCV, sparsity
       ↓
    Step 4: ML Recommender (Module 4)
       - Get pipeline recommendation
       - Optimal parameters
    
    STAGE 3: DE Analysis (External R + Dashboard)
    ──────────────────────────────
    Step 5: Run DE Analysis (R/Bioconductor)
       - DESeq2, edgeR, limma (external)
       ↓
    Step 6: Import DE Results (Module 7)
       - Load DESeq2/edgeR/limma
       - Standardize format
       ↓
    Step 7: Parameter Optimization (Module 8) [Optional]
       - Find optimal thresholds
       - 4 scientific methods
       ↓
    Step 8: Ensemble Analysis (Module 9)
       - Combine methods
       - Consensus genes
       ↓
    Step 9: Visualizations & Reports
       - Publication-ready outputs
       - Interactive plots
    ```
    """)

# Tips section
with st.expander("💡 Quick Tips"):
    st.markdown("""
    ### Getting the Best Results:
    
    **Quality Control:**
    - Always check for outliers before DE analysis
    - Review batch effects and apply correction if needed
    - Ensure library sizes are consistent
    
    **Ensemble Analysis:**
    - Use at least 2 methods for consensus
    - Higher consensus (all 3 methods) = higher confidence
    - Check direction agreement for each gene
    
    **Parameter Optimization:**
    - Use ground truth method if you have validated genes
    - FDR control method works without validation data
    - Stability method good for replicates
    
    **Visualization:**
    - Start with volcano plots for overview
    - Use MA plots to check for bias
    - PCA plots reveal sample structure
    """)

# Footer
add_footer()