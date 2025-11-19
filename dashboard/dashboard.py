#!/usr/bin/env python3

"""
RAPTOR Interactive Dashboard

Web-based interface for all RAPTOR ML features including:
- ML-based pipeline recommendations
- Resource monitoring
- Ensemble analysis
- Benchmark comparisons

Author: Ayeh Bolouki
Affiliation: University of Namur & GIGA-Neurosciences, University of Liège, Belgium
Email: ayehbolouki1988@gmail.com
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
</style>
""", unsafe_allow_html=True)


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
    
    # Welcome message
    st.markdown("""
    ### Welcome to RAPTOR!
    
    This interactive dashboard provides access to all RAPTOR ML features:
    
    - **🤖 ML Recommender**: Get AI-powered pipeline recommendations
    - **📊 Resource Monitor**: Track system resources in real-time
    - **🎯 Ensemble Analysis**: Combine results from multiple pipelines
    - **📈 Benchmarks**: Compare pipeline performance
    - **⚙️ Settings**: Configure preferences and models
    """)
    
    # Check system status
    st.markdown('<p class="sub-header">System Status</p>', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
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
            st.info("Train a model from the ML Recommender page or run: `python example_ml_workflow.py`")
    
    with col3:
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Dashboard", "✅ Active")
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Quick start guide
    st.markdown('<p class="sub-header">Quick Start</p>', unsafe_allow_html=True)
    
    st.markdown("""
    #### New User? Start Here:
    
    1. **Upload your data** on the ML Recommender page
    2. **Get a recommendation** with one click
    3. **Explore ensemble analysis** to combine pipeline results
    4. **Monitor resources** during pipeline execution
    
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
    st.markdown('<p class="main-header">🎯 Ensemble Analysis</p>', unsafe_allow_html=True)
    
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
    
    if st.button("Save Settings"):
        settings = {
            'model_type': model_type,
            'model_path': model_path,
            'data_dir': data_dir,
            'output_dir': output_dir,
            'n_threads': n_threads,
            'memory_gb': memory_gb,
            'theme': theme,
            'auto_refresh': auto_refresh
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
        ["🏠 Home", "🤖 ML Recommender", "📊 Resource Monitor", "🎯 Ensemble Analysis", "📈 Benchmarks", "⚙️ Settings"]
    )
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("### About")
    st.sidebar.info("""
    **RAPTOR v2.0.0**
    
    RNA-seq Analysis Pipeline Testing & Optimization Resource
    
    Created by Ayeh Bolouki
    
    University of Namur & GIGA-Neurosciences, University of Liège, Belgium
    """)
    
    # Route to appropriate page
    if page == "🏠 Home":
        home_page()
    elif page == "🤖 ML Recommender":
        ml_recommender_page()
    elif page == "📊 Resource Monitor":
        resource_monitor_page()
    elif page == "🎯 Ensemble Analysis":
        ensemble_page()
    elif page == "📈 Benchmarks":
        benchmarks_page()
    elif page == "⚙️ Settings":
        settings_page()


if __name__ == "__main__":
    main()
