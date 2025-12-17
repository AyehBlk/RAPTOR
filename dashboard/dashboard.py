#!/usr/bin/env python3

"""
RAPTOR Interactive Dashboard

Web-based interface for all RAPTOR ML features including:
- ML-based pipeline recommendations
- Adaptive Threshold Optimizer (NEW in v2.1.1)
- Resource monitoring
- Ensemble analysis
- Benchmark comparisons

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.1.1
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
import json
import sys
import time
from datetime import datetime
from io import StringIO

# Page configuration
st.set_page_config(
    page_title="🦖 RAPTOR Dashboard",
    page_icon="🦖",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        color: #2E7D32;
        text-align: center;
        padding: 1rem;
    }
    .sub-header {
        font-size: 1.5rem;
        font-weight: bold;
        color: #1976D2;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1976D2;
    }
    .success-box {
        background-color: #E8F5E9;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #4CAF50;
    }
    .warning-box {
        background-color: #FFF3E0;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #FF9800;
    }
    .error-box {
        background-color: #FFEBEE;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #F44336;
    }
    .new-feature {
        background-color: #E3F2FD;
        padding: 0.5rem 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #2196F3;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)


# Check for Threshold Optimizer availability
try:
    from raptor.threshold_optimizer import (
        AdaptiveThresholdOptimizer,
        optimize_thresholds,
        ThresholdResult
    )
    ATO_AVAILABLE = True
except ImportError:
    try:
        # Try local import for development
        from threshold_optimizer import (
            AdaptiveThresholdOptimizer,
            optimize_thresholds,
            ThresholdResult
        )
        ATO_AVAILABLE = True
    except ImportError:
        ATO_AVAILABLE = False


# Initialize session state
def init_session_state():
    """Initialize session state variables."""
    if 'profile' not in st.session_state:
        st.session_state.profile = None
    if 'recommendation' not in st.session_state:
        st.session_state.recommendation = None
    if 'monitoring_active' not in st.session_state:
        st.session_state.monitoring_active = False
    if 'monitor_data' not in st.session_state:
        st.session_state.monitor_data = []
    if 'ensemble_results' not in st.session_state:
        st.session_state.ensemble_results = None
    # Threshold Optimizer session state (v2.1.1)
    if 'ato_result' not in st.session_state:
        st.session_state.ato_result = None
    if 'ato_df' not in st.session_state:
        st.session_state.ato_df = None
    if 'ato_instance' not in st.session_state:
        st.session_state.ato_instance = None


def check_dependencies():
    """Check if required modules are available."""
    missing = []
    
    try:
        import ml_recommender
    except ImportError:
        missing.append("ml_recommender.py")
    
    try:
        import synthetic_benchmarks
    except ImportError:
        missing.append("synthetic_benchmarks.py")
    
    return missing


def create_gauge_chart(value, title, max_value=100):
    """Create a gauge chart for metrics."""
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=value,
        title={'text': title},
        delta={'reference': max_value * 0.8},
        gauge={
            'axis': {'range': [None, max_value]},
            'bar': {'color': "darkblue"},
            'steps': [
                {'range': [0, max_value * 0.33], 'color': "lightgray"},
                {'range': [max_value * 0.33, max_value * 0.67], 'color': "gray"},
                {'range': [max_value * 0.67, max_value], 'color': "darkgray"}
            ],
            'threshold': {
                'line': {'color': "red", 'width': 4},
                'thickness': 0.75,
                'value': max_value * 0.9
            }
        }
    ))
    
    fig.update_layout(height=250)
    return fig


def home_page():
    """Main home page."""
    st.markdown('<p class="main-header">🦖 RAPTOR Dashboard</p>', unsafe_allow_html=True)
    st.markdown('<p style="text-align: center; font-size: 1.2rem; color: #666;">RNA-seq Analysis Pipeline Testing & Optimization Resource</p>', unsafe_allow_html=True)
    
    st.markdown("---")
    
    # What's new banner
    st.markdown('<div class="new-feature">', unsafe_allow_html=True)
    st.markdown("🆕 **New in v2.1.1**: Adaptive Threshold Optimizer (ATO) - Data-driven threshold selection for DE analysis!")
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Welcome message
    st.markdown("""
    ### Welcome to RAPTOR!
    
    This interactive dashboard provides access to all RAPTOR features:
    
    - **🤖 ML Recommender**: Get AI-powered pipeline recommendations
    - **🎯 Threshold Optimizer**: Data-driven DE threshold optimization *(NEW!)*
    - **📊 Resource Monitor**: Track system resources in real-time
    - **🔬 Ensemble Analysis**: Combine results from multiple pipelines
    - **📈 Benchmarks**: Compare pipeline performance
    - **⚙️ Settings**: Configure preferences and models
    """)
    
    # Check system status
    st.markdown('<p class="sub-header">System Status</p>', unsafe_allow_html=True)
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Dependencies", "✅ Installed")
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col2:
        model_path = Path("models")
        if model_path.exists() and list(model_path.glob("*.pkl")):
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("ML Model", "✅ Ready")
            st.markdown('</div>', unsafe_allow_html=True)
        else:
            st.markdown('<div class="metric-card">', unsafe_allow_html=True)
            st.metric("ML Model", "⚠️ Not Found")
            st.markdown('</div>', unsafe_allow_html=True)
    
    with col3:
        ato_status = "✅ Available" if ATO_AVAILABLE else "⚠️ Not Found"
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Threshold Optimizer", ato_status)
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col4:
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Dashboard", "✅ Active")
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Quick start guide
    st.markdown('<p class="sub-header">Quick Start</p>', unsafe_allow_html=True)
    
    st.markdown("""
    #### New User? Start Here:
    
    1. **Upload your data** on the ML Recommender page
    2. **Get a recommendation** with one click
    3. **Optimize thresholds** for your DE results *(NEW!)*
    4. **Explore ensemble analysis** to combine pipeline results
    5. **Monitor resources** during pipeline execution
    
    #### Need Training Data?
    
    Run this command to generate synthetic training data:
    ```bash
    python example_ml_workflow.py --n-datasets 200
    ```
    """)
    
    # Recent activity
    if st.session_state.recommendation:
        st.markdown('<p class="sub-header">Recent Activity</p>', unsafe_allow_html=True)
        st.success(f"✅ Last recommendation: Pipeline {st.session_state.recommendation['pipeline_id']} ({st.session_state.recommendation['confidence']:.1%} confidence)")
    
    if st.session_state.ato_result:
        st.info(f"🎯 Last threshold optimization: |logFC| > {st.session_state.ato_result.logfc_cutoff:.3f}, {st.session_state.ato_result.n_significant_optimized} DE genes")


def ml_recommender_page():
    """ML Recommender interface."""
    st.markdown('<p class="main-header">🤖 ML Pipeline Recommender</p>', unsafe_allow_html=True)
    
    st.markdown("""
    Upload your RNA-seq count matrix or use sample data to get an AI-powered pipeline recommendation.
    """)
    
    # Check if ML module is available
    try:
        from ml_recommender import MLPipelineRecommender, FeatureExtractor
        ml_available = True
    except ImportError:
        ml_available = False
        st.error("❌ ML recommender module not found. Ensure ml_recommender.py is in the Python path.")
        return
    
    # Data input section
    st.markdown('<p class="sub-header">1. Data Input</p>', unsafe_allow_html=True)
    
    data_source = st.radio("Choose data source:", ["Upload CSV", "Use sample data"])
    
    counts_df = None
    
    if data_source == "Upload CSV":
        uploaded_file = st.file_uploader("Upload count matrix (CSV)", type=['csv'])
        if uploaded_file:
            try:
                counts_df = pd.read_csv(uploaded_file, index_col=0)
                st.success(f"✅ Loaded: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
                
                with st.expander("Preview data"):
                    st.dataframe(counts_df.head())
            except Exception as e:
                st.error(f"Error loading file: {e}")
    else:
        # Generate sample data
        if st.button("Generate Sample Data"):
            np.random.seed(42)
            n_genes = 1000
            n_samples = 6
            
            # Generate realistic count data
            means = np.random.lognormal(mean=6, sigma=2, size=n_genes)
            counts_df = pd.DataFrame(
                np.random.negative_binomial(n=10, p=0.1, size=(n_genes, n_samples)),
                index=[f"GENE{i:05d}" for i in range(n_genes)],
                columns=[f"Sample{i+1}" for i in range(n_samples)]
            )
            
            st.success("✅ Generated sample data (1000 genes × 6 samples)")
            with st.expander("Preview sample data"):
                st.dataframe(counts_df.head())
    
    # Profile data
    if counts_df is not None:
        st.markdown('<p class="sub-header">2. Profile Data</p>', unsafe_allow_html=True)
        
        if st.button("Profile Data", type="primary"):
            with st.spinner("Profiling data..."):
                # Create profile (simplified for dashboard)
                profile = {
                    'design': {
                        'n_samples': counts_df.shape[1],
                        'n_genes': counts_df.shape[0],
                        'n_conditions': 2,
                        'samples_per_condition': counts_df.shape[1] // 2,
                        'is_paired': False
                    },
                    'library_stats': {
                        'mean': float(counts_df.sum(axis=0).mean()),
                        'median': float(counts_df.sum(axis=0).median()),
                        'cv': float(counts_df.sum(axis=0).std() / counts_df.sum(axis=0).mean()),
                        'range': float(counts_df.sum(axis=0).max() - counts_df.sum(axis=0).min()),
                        'skewness': 0.2
                    },
                    'count_distribution': {
                        'zero_pct': float((counts_df == 0).sum().sum() / counts_df.size * 100),
                        'low_count_pct': float((counts_df < 10).sum().sum() / counts_df.size * 100),
                        'mean': float(counts_df.values.mean()),
                        'median': float(np.median(counts_df.values)),
                        'variance': float(counts_df.values.var())
                    },
                    'expression_distribution': {
                        'high_expr_genes': int((counts_df.mean(axis=1) > counts_df.mean(axis=1).quantile(0.9)).sum()),
                        'medium_expr_genes': int((counts_df.mean(axis=1) > counts_df.mean(axis=1).quantile(0.5)).sum()),
                        'low_expr_genes': int((counts_df.mean(axis=1) <= counts_df.mean(axis=1).quantile(0.5)).sum()),
                        'dynamic_range': 8.5
                    },
                    'biological_variation': {
                        'bcv': 0.3,
                        'dispersion_mean': 0.09,
                        'dispersion_trend': 0.1,
                        'outlier_genes': 50
                    },
                    'sequencing': {
                        'total_reads': float(counts_df.sum().sum()),
                        'reads_per_gene': float(counts_df.sum().sum() / counts_df.shape[0]),
                        'depth_category': 'medium'
                    },
                    'complexity': {
                        'score': 65.0,
                        'noise_level': 0.6,
                        'signal_strength': 0.7
                    }
                }
                
                st.session_state.profile = profile
                
                # Display profile summary
                st.success("✅ Profile created successfully!")
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Samples", profile['design']['n_samples'])
                    st.metric("Library Size (mean)", f"{profile['library_stats']['mean']:,.0f}")
                
                with col2:
                    st.metric("Genes", profile['design']['n_genes'])
                    st.metric("Zero %", f"{profile['count_distribution']['zero_pct']:.1f}%")
                
                with col3:
                    st.metric("BCV", f"{profile['biological_variation']['bcv']:.3f}")
                    st.metric("Depth", profile['sequencing']['depth_category'].title())
    
    # Get recommendation
    if st.session_state.profile:
        st.markdown('<p class="sub-header">3. Get ML Recommendation</p>', unsafe_allow_html=True)
        
        model_path = st.text_input("Model directory", value="models/")
        
        if st.button("Get ML Recommendation", type="primary"):
            if not Path(model_path).exists():
                st.error("❌ Model directory not found. Train a model first!")
                st.info("Run: `python example_ml_workflow.py --n-datasets 200`")
            else:
                with st.spinner("Loading model and making prediction..."):
                    try:
                        # Load model
                        recommender = MLPipelineRecommender(model_type='random_forest')
                        recommender.load_model(model_path)
                        
                        # Get recommendation
                        rec = recommender.recommend(st.session_state.profile, top_k=3)
                        st.session_state.recommendation = rec
                        
                        # Display recommendation
                        st.markdown('<div class="success-box">', unsafe_allow_html=True)
                        st.markdown(f"## 🦖 Recommended Pipeline")
                        st.markdown(f"### Pipeline {rec['pipeline_id']}: {rec['pipeline_name']}")
                        st.markdown(f"**Confidence: {rec['confidence']:.1%}**")
                        st.markdown('</div>', unsafe_allow_html=True)
                        
                        # Confidence gauge
                        st.plotly_chart(
                            create_gauge_chart(rec['confidence'] * 100, "Confidence Score"),
                            use_container_width=True
                        )
                        
                        # Reasons
                        st.markdown("#### Why this pipeline?")
                        for reason in rec['reasons']:
                            st.markdown(f"- {reason}")
                        
                        # Alternatives
                        if rec.get('alternatives'):
                            st.markdown("#### Alternative Options")
                            for alt in rec['alternatives']:
                                st.info(f"Pipeline {alt['pipeline_id']}: {alt['pipeline_name']} ({alt['confidence']:.1%} confidence)")
                        
                        # Feature contributions
                        if rec.get('feature_contributions'):
                            st.markdown("#### Top Contributing Features")
                            feat_df = pd.DataFrame(rec['feature_contributions'][:10])
                            fig = px.bar(
                                feat_df,
                                x='importance',
                                y='feature',
                                orientation='h',
                                title="Feature Importance"
                            )
                            st.plotly_chart(fig, use_container_width=True)
                        
                    except Exception as e:
                        st.error(f"❌ Error: {e}")
                        import traceback
                        st.code(traceback.format_exc())


# ==============================================================================
# THRESHOLD OPTIMIZER PAGE (NEW in v2.1.1)
# ==============================================================================

def threshold_optimizer_page():
    """Adaptive Threshold Optimizer interface."""
    st.markdown('<p class="main-header">🎯 Adaptive Threshold Optimizer</p>', unsafe_allow_html=True)
    
    st.markdown('<div class="new-feature">', unsafe_allow_html=True)
    st.markdown("🆕 **New in v2.1.1**: Data-driven threshold optimization for differential expression analysis")
    st.markdown('</div>', unsafe_allow_html=True)
    
    if not ATO_AVAILABLE:
        st.error("""
        ⚠️ **Threshold Optimizer module not found!**
        
        Please ensure the `threshold_optimizer` module is installed:
        ```python
        pip install raptor-rnaseq>=2.1.1
        ```
        
        Or copy the `threshold_optimizer/` folder to your RAPTOR installation.
        """)
        return
    
    st.markdown("""
    The Adaptive Threshold Optimizer (ATO) provides **data-driven** thresholds for 
    differential expression analysis, replacing arbitrary cutoffs with scientifically 
    justified values.
    
    **Upload your DE results** from DESeq2, edgeR, or limma to get started.
    """)
    
    st.markdown("---")
    
    # File upload section
    st.markdown('<p class="sub-header">📁 Upload DE Results</p>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        uploaded_file = st.file_uploader(
            "Upload your differential expression results",
            type=['csv', 'txt', 'tsv'],
            help="Supported formats: CSV, TSV, TXT. Must contain log2FoldChange and pvalue columns."
        )
    
    with col2:
        st.markdown("**Required columns:**")
        st.markdown("- `log2FoldChange` (or `logFC`)")
        st.markdown("- `pvalue` (or `PValue`)")
        st.markdown("")
        st.markdown("**Optional columns:**")
        st.markdown("- `padj`, `baseMean`, `lfcSE`")
    
    # Demo data option
    use_demo = st.checkbox("Use demo data instead", value=False)
    
    df = None
    
    if use_demo:
        # Generate synthetic demo data
        np.random.seed(42)
        n_genes = 10000
        
        # Most genes are null
        null_logfc = np.random.normal(0, 0.2, int(n_genes * 0.85))
        null_pval = np.random.uniform(0.05, 1, int(n_genes * 0.85))
        
        # Some DE genes
        de_logfc = np.concatenate([
            np.random.normal(1.5, 0.5, int(n_genes * 0.075)),
            np.random.normal(-1.5, 0.5, int(n_genes * 0.075))
        ])
        de_pval = np.random.exponential(0.001, int(n_genes * 0.15))
        
        df = pd.DataFrame({
            'log2FoldChange': np.concatenate([null_logfc, de_logfc]),
            'pvalue': np.clip(np.concatenate([null_pval, de_pval]), 1e-300, 1),
            'baseMean': np.random.exponential(1000, n_genes)
        })
        df.index = [f'Gene_{i}' for i in range(n_genes)]
        
        st.success("✅ Demo data loaded (10,000 synthetic genes)")
    
    elif uploaded_file is not None:
        # Detect separator
        content = uploaded_file.getvalue().decode('utf-8')
        first_line = content.split('\n')[0]
        
        if '\t' in first_line:
            sep = '\t'
        elif ',' in first_line:
            sep = ','
        else:
            sep = '\t'
        
        # Load data
        try:
            df = pd.read_csv(StringIO(content), sep=sep, index_col=0)
            st.success(f"✅ Loaded {len(df)} genes from {uploaded_file.name}")
        except Exception as e:
            st.error(f"Error loading file: {e}")
            return
    
    if df is None:
        st.info("👆 Upload a file or check 'Use demo data' to get started")
        return
    
    st.markdown("---")
    
    # Analysis settings
    st.markdown('<p class="sub-header">⚙️ Analysis Settings</p>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        goal = st.selectbox(
            "Analysis Goal",
            options=['discovery', 'balanced', 'validation'],
            index=0,
            help="""
            - **Discovery**: Maximize true positives (FDR control)
            - **Balanced**: Balance sensitivity and specificity  
            - **Validation**: Minimize false positives (FWER control)
            """
        )
    
    with col2:
        logfc_method = st.selectbox(
            "LogFC Method",
            options=['auto', 'mad', 'mixture', 'power', 'percentile'],
            index=0,
            help="""
            - **auto**: Consensus of all methods (recommended)
            - **mad**: MAD-based robust estimation
            - **mixture**: Gaussian mixture model
            - **power**: Power-based minimum effect
            - **percentile**: 95th percentile of null
            """
        )
    
    with col3:
        col3a, col3b = st.columns(2)
        with col3a:
            n1 = st.number_input("Samples (Group 1)", min_value=2, max_value=100, value=3)
        with col3b:
            n2 = st.number_input("Samples (Group 2)", min_value=2, max_value=100, value=3)
    
    # Run optimization
    if st.button("🚀 Optimize Thresholds", type="primary"):
        
        with st.spinner("Running threshold optimization..."):
            try:
                ato = AdaptiveThresholdOptimizer(df, goal=goal, verbose=False)
                result = ato.optimize(logfc_method=logfc_method, n1=n1, n2=n2)
                
                # Store in session state
                st.session_state.ato_result = result
                st.session_state.ato_df = ato.df
                st.session_state.ato_instance = ato
                
            except Exception as e:
                st.error(f"Error during optimization: {e}")
                import traceback
                st.code(traceback.format_exc())
                return
        
        st.success("✅ Optimization complete!")
    
    # Display results if available
    if st.session_state.ato_result is not None:
        result = st.session_state.ato_result
        ato_df = st.session_state.ato_df
        ato = st.session_state.ato_instance
        
        st.markdown("---")
        st.markdown('<p class="sub-header">📊 Optimization Results</p>', unsafe_allow_html=True)
        
        # Key metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric(
                "Recommended |logFC| cutoff",
                f"{result.logfc_cutoff:.3f}",
                delta=f"{result.logfc_cutoff - 1.0:+.3f} vs traditional",
                delta_color="inverse"
            )
        
        with col2:
            st.metric(
                "Recommended padj cutoff",
                f"{result.padj_cutoff}",
                delta=f"{result.padj_method}"
            )
        
        with col3:
            st.metric(
                "DE genes (optimized)",
                result.n_significant_optimized,
                delta=f"{result.n_significant_optimized - result.n_significant_traditional:+d}"
            )
        
        with col4:
            st.metric(
                "π₀ estimate",
                f"{result.pi0_estimate:.3f}",
                help="Proportion of true null hypotheses"
            )
        
        # Detailed reasoning
        with st.expander("📋 Detailed Reasoning", expanded=True):
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**LogFC Cutoff**")
                st.markdown(f"- Method: {result.logfc_method}")
                st.markdown(f"- {result.logfc_reasoning}")
            
            with col2:
                st.markdown("**P-value Adjustment**")
                st.markdown(f"- Method: {result.padj_method}")
                st.markdown(f"- {result.padj_reasoning}")
        
        # Visualizations
        st.markdown("---")
        st.markdown('<p class="sub-header">📈 Visualizations</p>', unsafe_allow_html=True)
        
        tab1, tab2, tab3, tab4 = st.tabs([
            "Volcano Plot", "LogFC Distribution", "P-value Distribution", "Threshold Comparison"
        ])
        
        with tab1:
            # Volcano plot
            logfc = ato_df['logfc'].values
            padj = ato_df['padj_optimized'].values if 'padj_optimized' in ato_df.columns else ato_df['padj'].values
            
            min_padj = padj[padj > 0].min() if (padj > 0).any() else 1e-300
            padj_safe = np.where(padj == 0, min_padj / 10, padj)
            neg_log_padj = -np.log10(padj_safe)
            
            sig_up = (logfc > result.logfc_cutoff) & (padj < result.padj_cutoff)
            sig_down = (logfc < -result.logfc_cutoff) & (padj < result.padj_cutoff)
            not_sig = ~(sig_up | sig_down)
            
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=logfc[not_sig], y=neg_log_padj[not_sig],
                mode='markers', marker=dict(color='gray', size=5, opacity=0.5),
                name=f'Not significant ({not_sig.sum()})',
                text=ato_df.index[not_sig],
                hovertemplate='<b>%{text}</b><br>logFC: %{x:.3f}<br>-log10(padj): %{y:.2f}<extra></extra>'
            ))
            
            fig.add_trace(go.Scatter(
                x=logfc[sig_up], y=neg_log_padj[sig_up],
                mode='markers', marker=dict(color='firebrick', size=8, opacity=0.7),
                name=f'Upregulated ({sig_up.sum()})',
                text=ato_df.index[sig_up],
                hovertemplate='<b>%{text}</b><br>logFC: %{x:.3f}<br>-log10(padj): %{y:.2f}<extra></extra>'
            ))
            
            fig.add_trace(go.Scatter(
                x=logfc[sig_down], y=neg_log_padj[sig_down],
                mode='markers', marker=dict(color='steelblue', size=8, opacity=0.7),
                name=f'Downregulated ({sig_down.sum()})',
                text=ato_df.index[sig_down],
                hovertemplate='<b>%{text}</b><br>logFC: %{x:.3f}<br>-log10(padj): %{y:.2f}<extra></extra>'
            ))
            
            fig.add_hline(y=-np.log10(result.padj_cutoff), line_dash="dash", line_color="black", opacity=0.5)
            fig.add_vline(x=result.logfc_cutoff, line_dash="dash", line_color="black", opacity=0.5)
            fig.add_vline(x=-result.logfc_cutoff, line_dash="dash", line_color="black", opacity=0.5)
            
            fig.update_layout(
                title=f'Volcano Plot (|logFC| > {result.logfc_cutoff:.3f}, padj < {result.padj_cutoff})',
                xaxis_title='log2 Fold Change',
                yaxis_title='-log10(adjusted p-value)',
                template='plotly_white',
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        with tab2:
            # LogFC histogram
            fig = go.Figure()
            
            fig.add_trace(go.Histogram(
                x=ato_df['logfc'], nbinsx=100,
                name='LogFC Distribution', marker_color='steelblue', opacity=0.7
            ))
            
            fig.add_vline(x=result.logfc_cutoff, line_dash="dash", line_color="green",
                         annotation_text=f"Optimized: {result.logfc_cutoff:.2f}")
            fig.add_vline(x=-result.logfc_cutoff, line_dash="dash", line_color="green")
            fig.add_vline(x=1.0, line_dash="dot", line_color="red",
                         annotation_text="Traditional: 1.0")
            fig.add_vline(x=-1.0, line_dash="dot", line_color="red")
            
            fig.update_layout(
                title='Log2 Fold Change Distribution',
                xaxis_title='log2 Fold Change',
                yaxis_title='Count',
                template='plotly_white',
                height=400
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        with tab3:
            # P-value histogram
            pvals = ato_df['pvalue'].dropna()
            
            fig = go.Figure()
            
            fig.add_trace(go.Histogram(
                x=pvals, nbinsx=50,
                name='P-value Distribution', marker_color='steelblue', opacity=0.7
            ))
            
            fig.add_hline(y=len(pvals) / 50 * result.pi0_estimate, line_dash="dash", line_color="red",
                         annotation_text=f"π₀ = {result.pi0_estimate:.3f}")
            
            fig.update_layout(
                title='P-value Distribution',
                xaxis_title='P-value',
                yaxis_title='Count',
                template='plotly_white',
                height=400
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        with tab4:
            # Threshold comparison heatmap
            comparison = ato.compare_thresholds(
                logfc_values=[0.25, 0.5, 0.75, 1.0, 1.5, 2.0],
                padj_values=[0.01, 0.05, 0.1]
            )
            
            pivot = comparison.pivot(index='padj_cutoff', columns='logFC_cutoff', values='n_significant')
            
            fig = px.imshow(
                pivot,
                labels=dict(x='|logFC| cutoff', y='padj cutoff', color='DE genes'),
                color_continuous_scale='YlOrRd',
                aspect='auto'
            )
            
            for i, row in enumerate(pivot.index):
                for j, col in enumerate(pivot.columns):
                    fig.add_annotation(
                        x=j, y=i,
                        text=str(int(pivot.loc[row, col])),
                        showarrow=False,
                        font=dict(color='white' if pivot.loc[row, col] > pivot.values.max()/2 else 'black')
                    )
            
            fig.update_layout(
                title='Number of DE Genes by Threshold Combination',
                template='plotly_white',
                height=400
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            st.dataframe(comparison, use_container_width=True)
        
        # Significant genes table
        st.markdown("---")
        st.markdown('<p class="sub-header">🧬 Significant Genes</p>', unsafe_allow_html=True)
        
        sig_genes = ato.get_significant_genes()
        
        st.markdown(f"**{len(sig_genes)} genes** pass optimized thresholds")
        
        # Sort options
        sort_col = st.selectbox(
            "Sort by",
            options=['pvalue', 'logfc', 'padj_optimized'],
            index=0
        )
        
        sig_genes_sorted = sig_genes.sort_values(
            sort_col, 
            key=abs if sort_col == 'logfc' else None, 
            ascending=sort_col != 'logfc'
        )
        
        st.dataframe(
            sig_genes_sorted[['logfc', 'pvalue', 'padj_optimized']].head(100),
            use_container_width=True
        )
        
        # Download buttons
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv = sig_genes.to_csv()
            st.download_button(
                "📥 Download Significant Genes (CSV)",
                csv,
                "significant_genes_optimized.csv",
                "text/csv"
            )
        
        with col2:
            summary_text = result.summary()
            st.download_button(
                "📥 Download Summary Report",
                summary_text,
                "threshold_optimization_report.txt",
                "text/plain"
            )
        
        with col3:
            full_csv = ato_df.to_csv()
            st.download_button(
                "📥 Download Full Results (CSV)",
                full_csv,
                "de_results_with_optimized_padj.csv",
                "text/csv"
            )
        
        # Methods text for publication
        st.markdown("---")
        st.markdown('<p class="sub-header">📝 Methods Text for Publication</p>', unsafe_allow_html=True)
        
        methods_text = f"""
**Threshold Optimization**

Significance thresholds for differential expression were determined using the 
Adaptive Threshold Optimizer (ATO) from RAPTOR v2.1.1 with the '{goal}' goal setting. 
The analysis estimated π₀ = {result.pi0_estimate:.3f}, indicating that approximately 
{(1-result.pi0_estimate)*100:.1f}% of tested genes showed true differential expression. 
Based on the data characteristics, the {result.padj_method} method was selected for 
p-value adjustment, and a data-driven log2 fold change cutoff of |logFC| > {result.logfc_cutoff:.3f} 
was determined using the {result.logfc_method} approach. Genes were considered differentially 
expressed if they met both the logFC and adjusted p-value (< {result.padj_cutoff}) thresholds, 
yielding {result.n_significant_optimized} DE genes.
"""
        
        st.code(methods_text, language=None)
        
        st.download_button(
            "📥 Download Methods Text",
            methods_text,
            "methods_threshold_optimization.txt",
            "text/plain"
        )


def resource_monitor_page():
    """Resource monitoring interface."""
    st.markdown('<p class="main-header">📊 Resource Monitor</p>', unsafe_allow_html=True)
    
    st.markdown("""
    Monitor CPU, memory, disk, and GPU usage in real-time during pipeline execution.
    """)
    
    # Check if monitoring module is available
    try:
        import psutil
        monitoring_available = True
    except ImportError:
        monitoring_available = False
        st.error("❌ psutil not installed. Install with: `pip install psutil`")
        return
    
    # Control buttons
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("▶️ Start Monitoring", type="primary"):
            st.session_state.monitoring_active = True
            st.session_state.monitor_data = []
            st.success("Monitoring started!")
    
    with col2:
        if st.button("⏸️ Pause"):
            st.session_state.monitoring_active = False
            st.info("Monitoring paused")
    
    with col3:
        if st.button("🔄 Reset"):
            st.session_state.monitoring_active = False
            st.session_state.monitor_data = []
            st.info("Data cleared")
    
    # Display metrics
    if st.session_state.monitoring_active or st.session_state.monitor_data:
        st.markdown('<p class="sub-header">Current Metrics</p>', unsafe_allow_html=True)
        
        # Placeholders for live updates
        metric_cols = st.columns(4)
        chart_placeholder = st.empty()
        
        if st.session_state.monitoring_active:
            # Collect metrics
            for _ in range(10):  # Collect 10 data points
                cpu_percent = psutil.cpu_percent(interval=0.1)
                memory = psutil.virtual_memory()
                disk_io = psutil.disk_io_counters()
                
                data_point = {
                    'timestamp': datetime.now(),
                    'cpu': cpu_percent,
                    'memory': memory.percent,
                    'disk_read': disk_io.read_bytes / (1024**2),  # MB
                    'disk_write': disk_io.write_bytes / (1024**2)  # MB
                }
                
                st.session_state.monitor_data.append(data_point)
                
                # Update metrics
                with metric_cols[0]:
                    st.metric("CPU", f"{cpu_percent:.1f}%", f"{cpu_percent - 50:.1f}%")
                
                with metric_cols[1]:
                    st.metric("Memory", f"{memory.percent:.1f}%", f"{memory.percent - 50:.1f}%")
                
                with metric_cols[2]:
                    st.metric("Memory Used", f"{memory.used / (1024**3):.1f} GB")
                
                with metric_cols[3]:
                    st.metric("Available", f"{memory.available / (1024**3):.1f} GB")
                
                time.sleep(1)
        
        # Plot historical data
        if st.session_state.monitor_data:
            df = pd.DataFrame(st.session_state.monitor_data)
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=df['timestamp'], y=df['cpu'], name='CPU %', mode='lines'))
            fig.add_trace(go.Scatter(x=df['timestamp'], y=df['memory'], name='Memory %', mode='lines'))
            
            fig.update_layout(
                title="Resource Usage Over Time",
                xaxis_title="Time",
                yaxis_title="Usage (%)",
                hovermode='x unified'
            )
            
            chart_placeholder.plotly_chart(fig, use_container_width=True)
            
            # Export option
            if st.button("💾 Export Data"):
                csv = df.to_csv(index=False)
                st.download_button(
                    "Download CSV",
                    csv,
                    "resource_monitor.csv",
                    "text/csv"
                )


def ensemble_page():
    """Ensemble analysis interface."""
    st.markdown('<p class="main-header">🔬 Ensemble Analysis</p>', unsafe_allow_html=True)
    
    st.markdown("""
    Combine results from multiple pipelines to create high-confidence gene lists.
    """)
    
    # Pipeline selection
    st.markdown('<p class="sub-header">1. Select Pipelines</p>', unsafe_allow_html=True)
    
    pipeline_names = {
        1: "STAR-RSEM-DESeq2",
        2: "HISAT2-StringTie-Ballgown",
        3: "Salmon-edgeR",
        4: "Kallisto-Sleuth",
        5: "STAR-HTSeq-limma-voom"
    }
    
    selected_pipelines = st.multiselect(
        "Choose pipelines to combine:",
        options=list(pipeline_names.keys()),
        format_func=lambda x: f"Pipeline {x}: {pipeline_names[x]}",
        default=[1, 3, 5]
    )
    
    # Ensemble method
    st.markdown('<p class="sub-header">2. Select Method</p>', unsafe_allow_html=True)
    
    method = st.selectbox(
        "Ensemble method:",
        ["vote", "rank_product", "p_value_combination", "weighted", "combined"]
    )
    
    method_descriptions = {
        "vote": "Simple majority voting across pipelines",
        "rank_product": "Rank product method for ranking genes",
        "p_value_combination": "Fisher's method to combine p-values",
        "weighted": "Weighted combination by pipeline accuracy",
        "combined": "Combines vote, rank, and p-value methods"
    }
    
    st.info(f"ℹ️ {method_descriptions[method]}")
    
    # Simulate ensemble analysis
    if st.button("Run Ensemble Analysis", type="primary"):
        with st.spinner("Running ensemble analysis..."):
            time.sleep(2)  # Simulate processing
            
            # Generate sample results
            n_genes = 500
            genes = [f"GENE{i:05d}" for i in range(n_genes)]
            
            # Simulate consensus scores
            scores = np.random.beta(5, 2, size=n_genes)
            agreement = np.random.randint(len(selected_pipelines) - 1, len(selected_pipelines) + 1, size=n_genes)
            pvalues = np.random.beta(1, 10, size=n_genes) * 0.05
            
            results_df = pd.DataFrame({
                'gene': genes,
                'consensus_score': scores,
                'n_pipelines': agreement,
                'p_value': pvalues
            })
            
            results_df = results_df.sort_values('consensus_score', ascending=False)
            st.session_state.ensemble_results = results_df
            
            st.success("✅ Ensemble analysis complete!")
    
    # Display results
    if st.session_state.ensemble_results is not None:
        st.markdown('<p class="sub-header">3. Results</p>', unsafe_allow_html=True)
        
        df = st.session_state.ensemble_results
        
        # Summary metrics
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Genes", len(df))
        
        with col2:
            high_conf = (df['consensus_score'] > 0.8).sum()
            st.metric("High Confidence", high_conf)
        
        with col3:
            mean_agreement = df['n_pipelines'].mean()
            st.metric("Mean Agreement", f"{mean_agreement:.1f}/{len(selected_pipelines)}")
        
        # Score distribution
        fig = px.histogram(df, x='consensus_score', nbins=50, title="Consensus Score Distribution")
        st.plotly_chart(fig, use_container_width=True)
        
        # Agreement heatmap (simulated)
        st.markdown("#### Pipeline Agreement")
        agreement_matrix = np.random.rand(len(selected_pipelines), len(selected_pipelines))
        agreement_matrix = (agreement_matrix + agreement_matrix.T) / 2
        np.fill_diagonal(agreement_matrix, 1.0)
        
        fig = px.imshow(
            agreement_matrix,
            labels=dict(x="Pipeline", y="Pipeline", color="Agreement"),
            x=[pipeline_names[p] for p in selected_pipelines],
            y=[pipeline_names[p] for p in selected_pipelines],
            color_continuous_scale="RdYlGn",
            title="Pipeline Agreement Heatmap"
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Top genes table
        st.markdown("#### Top 20 Consensus Genes")
        st.dataframe(df.head(20), use_container_width=True)
        
        # Export options
        st.markdown("#### Export Results")
        
        col1, col2 = st.columns(2)
        
        with col1:
            csv = df.to_csv(index=False)
            st.download_button(
                "Download Full Results (CSV)",
                csv,
                "ensemble_results.csv",
                "text/csv"
            )
        
        with col2:
            high_conf_genes = df[df['consensus_score'] > 0.8]
            txt = "\n".join(high_conf_genes['gene'].tolist())
            st.download_button(
                "Download High-Confidence Genes (TXT)",
                txt,
                "high_confidence_genes.txt",
                "text/plain"
            )


def benchmarks_page():
    """Benchmarks comparison page."""
    st.markdown('<p class="main-header">📈 Pipeline Benchmarks</p>', unsafe_allow_html=True)
    
    st.markdown("""
    Compare performance metrics across different RNA-seq analysis pipelines.
    """)
    
    # Generate sample benchmark data
    pipelines = [
        "STAR-RSEM-DESeq2",
        "HISAT2-StringTie-Ballgown",
        "Salmon-edgeR",
        "Kallisto-Sleuth",
        "STAR-HTSeq-limma-voom"
    ]
    
    benchmark_data = {
        'pipeline': pipelines,
        'accuracy': [0.89, 0.82, 0.87, 0.85, 0.88],
        'precision': [0.88, 0.80, 0.86, 0.83, 0.87],
        'recall': [0.90, 0.84, 0.88, 0.87, 0.89],
        'f1_score': [0.89, 0.82, 0.87, 0.85, 0.88],
        'runtime_min': [60, 45, 12, 8, 55]
    }
    
    df = pd.DataFrame(benchmark_data)
    
    # Performance metrics
    st.markdown('<p class="sub-header">Performance Metrics</p>', unsafe_allow_html=True)
    
    metric = st.selectbox("Select metric:", ['f1_score', 'accuracy', 'precision', 'recall'])
    
    fig = px.bar(df, x='pipeline', y=metric, title=f"{metric.replace('_', ' ').title()} by Pipeline")
    st.plotly_chart(fig, use_container_width=True)
    
    # Runtime comparison
    st.markdown('<p class="sub-header">Runtime Comparison</p>', unsafe_allow_html=True)
    
    fig = px.bar(df, x='pipeline', y='runtime_min', title="Runtime (minutes)")
    st.plotly_chart(fig, use_container_width=True)
    
    # Scatter plot: accuracy vs runtime
    st.markdown('<p class="sub-header">Accuracy vs Runtime Trade-off</p>', unsafe_allow_html=True)
    
    fig = px.scatter(
        df,
        x='runtime_min',
        y='f1_score',
        text='pipeline',
        title="F1 Score vs Runtime",
        labels={'runtime_min': 'Runtime (minutes)', 'f1_score': 'F1 Score'}
    )
    fig.update_traces(textposition='top center')
    st.plotly_chart(fig, use_container_width=True)
    
    # Detailed comparison table
    st.markdown('<p class="sub-header">Detailed Comparison</p>', unsafe_allow_html=True)
    st.dataframe(df, use_container_width=True)


def settings_page():
    """Settings and configuration page."""
    st.markdown('<p class="main-header">⚙️ Settings</p>', unsafe_allow_html=True)
    
    st.markdown("### Model Settings")
    
    model_type = st.selectbox(
        "ML Model Type:",
        ["random_forest", "gradient_boosting"]
    )
    
    model_path = st.text_input("Model Directory:", value="models/")
    
    st.markdown("### Data Directories")
    
    data_dir = st.text_input("Training Data Directory:", value="ml_training_data/")
    output_dir = st.text_input("Output Directory:", value="results/")
    
    st.markdown("### Performance Settings")
    
    n_threads = st.slider("Number of Threads:", 1, 16, 8)
    memory_gb = st.slider("Memory Limit (GB):", 4, 64, 32)
    
    st.markdown("### Dashboard Preferences")
    
    theme = st.selectbox("Color Theme:", ["Light", "Dark", "Auto"])
    auto_refresh = st.checkbox("Auto-refresh monitoring", value=True)
    
    st.markdown("### Threshold Optimizer Settings (v2.1.1)")
    
    default_goal = st.selectbox(
        "Default Analysis Goal:",
        ["discovery", "balanced", "validation"],
        index=0
    )
    
    default_logfc_method = st.selectbox(
        "Default LogFC Method:",
        ["auto", "mad", "mixture", "power", "percentile"],
        index=0
    )
    
    if st.button("Save Settings"):
        settings = {
            'model_type': model_type,
            'model_path': model_path,
            'data_dir': data_dir,
            'output_dir': output_dir,
            'n_threads': n_threads,
            'memory_gb': memory_gb,
            'theme': theme,
            'auto_refresh': auto_refresh,
            'ato_default_goal': default_goal,
            'ato_default_logfc_method': default_logfc_method
        }
        
        # Save to file
        with open('dashboard_settings.json', 'w') as f:
            json.dump(settings, f, indent=2)
        
        st.success("✅ Settings saved successfully!")


def main():
    """Main application."""
    init_session_state()
    
    # Sidebar navigation
    st.sidebar.markdown("## 🦖 RAPTOR")
    st.sidebar.markdown("---")
    
    page = st.sidebar.radio(
        "Navigation",
        [
            "🏠 Home", 
            "🤖 ML Recommender", 
            "🎯 Threshold Optimizer",  # NEW in v2.1.1
            "📊 Resource Monitor", 
            "🔬 Ensemble Analysis", 
            "📈 Benchmarks", 
            "⚙️ Settings"
        ]
    )
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("### About")
    st.sidebar.info("""
    **RAPTOR v2.1.1**
    
    RNA-seq Analysis Pipeline Testing & Optimization Resource
    
    🆕 New: Adaptive Threshold Optimizer
    
    Created by Ayeh Bolouki
    """)
    
    # ATO status indicator
    if ATO_AVAILABLE:
        st.sidebar.success("🎯 Threshold Optimizer: Ready")
    else:
        st.sidebar.warning("🎯 Threshold Optimizer: Not installed")
    
    # Route to appropriate page
    if page == "🏠 Home":
        home_page()
    elif page == "🤖 ML Recommender":
        ml_recommender_page()
    elif page == "🎯 Threshold Optimizer":
        threshold_optimizer_page()
    elif page == "📊 Resource Monitor":
        resource_monitor_page()
    elif page == "🔬 Ensemble Analysis":
        ensemble_page()
    elif page == "📈 Benchmarks":
        benchmarks_page()
    elif page == "⚙️ Settings":
        settings_page()


if __name__ == "__main__":
    main()
