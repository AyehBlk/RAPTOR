"""
RAPTOR Dashboard - Interactive Visualizations

Comprehensive visualization suite for differential expression analysis with
publication-ready interactive plots.

Author: Ayeh Bolouki
Version: 2.2.0
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
from pathlib import Path
from typing import Dict, List, Optional
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# Page config
st.set_page_config(
    page_title="RAPTOR - Visualizations",
    page_icon="📊",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("📊 RAPTOR Visualization Suite")
st.markdown("### Publication-ready interactive plots for differential expression analysis")

# Help section
with st.expander("ℹ️ How to use visualizations"):
    st.markdown("""
    **Purpose:** Create publication-quality visualizations from your RAPTOR analysis.
    
    ## 📊 **Available Visualizations**
    
    ### **With DE Results (Module 7):**
    
    **1. Volcano Plot** 🌋
    - Shows Log2FC vs -log10(p-value)
    - Highlights significant genes
    - Interactive gene labeling
    - **Use for:** Overview of DE results
    
    **2. MA Plot** 📊
    - Shows Log2FC vs mean expression
    - Identifies expression-dependent bias
    - Interactive exploration
    - **Use for:** Quality control, bias detection
    
    ### **With Count Matrix (Upload):**
    
    **3. PCA Plot** 🎨
    - Principal component analysis
    - Sample clustering visualization
    - Color by conditions
    - **Use for:** Sample quality, batch effects
    
    **4. Heatmap** 🔥
    - Hierarchical clustering
    - Top variable genes
    - Sample/gene clustering
    - **Use for:** Pattern identification
    
    **5. Gene Expression Profile** 📈
    - Individual gene expression
    - Across all samples
    - Grouped by condition
    - **Use for:** Specific gene validation
    
    **6. Sample Correlation** 🔗
    - Sample-to-sample correlation
    - Heatmap with clustering
    - Identifies outliers
    - **Use for:** QC, batch effects
    
    ## 💡 **Tips**
    
    - Adjust thresholds in real-time
    - Export in PNG, SVG, or PDF format
    - All plots are interactive (zoom, pan, select)
    - Use for publications and presentations
    """)

st.markdown("---")

# Check data availability
st.markdown("## 📊 Available Data")

has_de_results = st.session_state.get('m7_complete', False)
has_counts = False  # Count matrix usually not stored in session state

col1, col2 = st.columns(2)

with col1:
    if has_de_results:
        st.success("✅ DE Results available (Module 7)")
        de_results = st.session_state.get('m7_de_results', {})
        st.info(f"**Files:** {len(de_results)} DE result file(s)")
    else:
        st.warning("⬜ No DE results")
        st.caption("Run Module 7 (Import DE) first")

with col2:
    st.info("ℹ️ Count matrices via file upload")
    st.caption("Upload normalized count matrix for PCA/Heatmap")

st.markdown("---")

# Data selection
st.markdown("## 1️⃣ Select Data Source")

# For DE visualizations
de_source = None
de_data = None

if has_de_results:
    de_results_dict = st.session_state.get('m7_de_results', {})
    
    if len(de_results_dict) > 0:
        selected_file = st.selectbox(
            "Select DE result file:",
            list(de_results_dict.keys()),
            help="Choose which DE analysis to visualize"
        )
        
        de_source = "session"
        de_result = de_results_dict[selected_file]
        de_data = de_result.results_df.copy()
        
        st.success(f"✅ Using {selected_file} ({len(de_data):,} genes)")

# For count matrix visualizations (file upload)
st.markdown("### 📤 Upload Count Matrix (Optional)")
st.caption("For PCA, Heatmap, Expression Profile, and Correlation plots")

counts_file = st.file_uploader(
    "Upload normalized count matrix (CSV)",
    type=['csv'],
    help="Genes as rows, samples as columns"
)

metadata_file = st.file_uploader(
    "Upload sample metadata (CSV) - Optional",
    type=['csv'],
    help="Must have 'sample' column"
)

counts_data = None
metadata_data = None

if counts_file is not None:
    try:
        counts_data = pd.read_csv(counts_file, index_col=0)
        st.success(f"✅ Loaded count matrix: {counts_data.shape[0]} genes × {counts_data.shape[1]} samples")
    except Exception as e:
        st.error(f"❌ Error loading counts: {str(e)}")

if metadata_file is not None:
    try:
        metadata_data = pd.read_csv(metadata_file)
        if 'sample' not in metadata_data.columns:
            st.warning("⚠️ Metadata must have 'sample' column")
        else:
            st.success(f"✅ Loaded metadata for {len(metadata_data)} samples")
    except Exception as e:
        st.error(f"❌ Error loading metadata: {str(e)}")

st.markdown("---")

# Plot settings
st.markdown("## 2️⃣ Plot Settings")

col1, col2, col3 = st.columns(3)

with col1:
    fdr_threshold = st.slider(
        "FDR Threshold",
        0.001, 0.20, 0.05, 0.001,
        format="%.3f",
        help="Significance threshold"
    )

with col2:
    lfc_threshold = st.slider(
        "Log2FC Threshold",
        0.0, 3.0, 1.0, 0.1,
        format="%.1f",
        help="Fold change threshold"
    )
    fc = 2 ** lfc_threshold
    st.caption(f"Fold change: {fc:.2f}×")

with col3:
    point_size = st.slider(
        "Point Size",
        2, 10, 5, 1,
        help="Size of points in plots"
    )

st.markdown("---")

# Plot type selection
st.markdown("## 3️⃣ Select Visualization")

available_plots = []
if de_data is not None:
    available_plots.extend(["🌋 Volcano Plot", "📊 MA Plot"])
if counts_data is not None:
    available_plots.extend(["🎨 PCA Plot", "🔥 Heatmap", "📈 Gene Expression Profile", "🔗 Sample Correlation"])

if len(available_plots) == 0:
    st.warning("⚠️ No data available for visualization. Please:")
    st.markdown("- Run Module 7 (Import DE) for Volcano/MA plots")
    st.markdown("- Upload count matrix for PCA/Heatmap/Expression/Correlation plots")
    st.stop()

plot_type = st.selectbox(
    "Choose visualization:",
    available_plots
)

st.markdown("---")

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_volcano_plot(df, fdr_threshold, lfc_threshold, point_size):
    """Create interactive volcano plot."""
    
    # Calculate -log10(padj)
    df = df.copy()
    df['-log10(padj)'] = -np.log10(df['adjusted_p_value'].replace(0, 1e-300))
    
    # Assign colors
    colors = []
    for _, row in df.iterrows():
        if pd.isna(row['adjusted_p_value']) or pd.isna(row['log2_fold_change']):
            colors.append('Not significant')
        elif row['adjusted_p_value'] < fdr_threshold and abs(row['log2_fold_change']) >= lfc_threshold:
            if row['log2_fold_change'] > 0:
                colors.append('Upregulated')
            else:
                colors.append('Downregulated')
        else:
            colors.append('Not significant')
    
    df['Significance'] = colors
    
    # Create figure
    fig = go.Figure()
    
    # Add traces for each category
    color_map = {
        'Upregulated': '#d62728',
        'Downregulated': '#1f77b4',
        'Not significant': '#7f7f7f'
    }
    
    for category in ['Not significant', 'Downregulated', 'Upregulated']:
        df_cat = df[df['Significance'] == category]
        
        fig.add_trace(go.Scatter(
            x=df_cat['log2_fold_change'],
            y=df_cat['-log10(padj)'],
            mode='markers',
            name=f"{category} ({len(df_cat):,})",
            marker=dict(
                color=color_map[category],
                size=point_size,
                opacity=0.6 if category == 'Not significant' else 0.8
            ),
            text=df_cat.index if hasattr(df, 'index') else df_cat.get('gene_id', ''),
            hovertemplate='<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}<extra></extra>'
        ))
    
    # Add threshold lines
    fig.add_hline(
        y=-np.log10(fdr_threshold),
        line_dash="dash",
        line_color="gray",
        annotation_text=f"FDR={fdr_threshold}",
        annotation_position="right"
    )
    
    fig.add_vline(
        x=lfc_threshold,
        line_dash="dash",
        line_color="gray"
    )
    
    fig.add_vline(
        x=-lfc_threshold,
        line_dash="dash",
        line_color="gray"
    )
    
    fig.update_layout(
        title="Volcano Plot: Differential Expression",
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10(Adjusted P-value)",
        template="plotly_white",
        hovermode='closest',
        height=600
    )
    
    return fig

def create_ma_plot(df, fdr_threshold, lfc_threshold, point_size):
    """Create interactive MA plot."""
    
    df = df.copy()
    
    # Assign colors
    colors = []
    for _, row in df.iterrows():
        if pd.isna(row['adjusted_p_value']) or pd.isna(row['log2_fold_change']):
            colors.append('Not significant')
        elif row['adjusted_p_value'] < fdr_threshold and abs(row['log2_fold_change']) >= lfc_threshold:
            if row['log2_fold_change'] > 0:
                colors.append('Upregulated')
            else:
                colors.append('Downregulated')
        else:
            colors.append('Not significant')
    
    df['Significance'] = colors
    
    # Create figure
    fig = go.Figure()
    
    color_map = {
        'Upregulated': '#d62728',
        'Downregulated': '#1f77b4',
        'Not significant': '#7f7f7f'
    }
    
    for category in ['Not significant', 'Downregulated', 'Upregulated']:
        df_cat = df[df['Significance'] == category]
        
        fig.add_trace(go.Scatter(
            x=np.log10(df_cat['base_mean'] + 1),
            y=df_cat['log2_fold_change'],
            mode='markers',
            name=f"{category} ({len(df_cat):,})",
            marker=dict(
                color=color_map[category],
                size=point_size,
                opacity=0.6 if category == 'Not significant' else 0.8
            ),
            text=df_cat.index if hasattr(df, 'index') else df_cat.get('gene_id', ''),
            hovertemplate='<b>%{text}</b><br>Log10(Mean): %{x:.2f}<br>Log2FC: %{y:.2f}<extra></extra>'
        ))
    
    # Add threshold lines
    fig.add_hline(y=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_hline(y=-lfc_threshold, line_dash="dash", line_color="gray")
    
    fig.update_layout(
        title="MA Plot: Log2FC vs Mean Expression",
        xaxis_title="Log10(Mean Expression)",
        yaxis_title="Log2 Fold Change",
        template="plotly_white",
        hovermode='closest',
        height=600
    )
    
    return fig

def create_pca_plot(counts, metadata=None):
    """Create PCA plot."""
    
    # Standardize
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts.T)
    
    # PCA
    pca = PCA(n_components=min(10, counts.shape[1]))
    pca_coords = pca.fit_transform(counts_scaled)
    
    # Create dataframe
    pca_df = pd.DataFrame(
        pca_coords[:, :2],
        columns=['PC1', 'PC2'],
        index=counts.columns
    )
    
    # Add metadata if available
    if metadata is not None and 'sample' in metadata.columns:
        metadata = metadata.set_index('sample')
        for col in metadata.columns:
            if col in ['condition', 'group', 'treatment']:
                pca_df[col] = metadata.loc[pca_df.index, col]
                break
    
    # Create plot
    if 'condition' in pca_df.columns or 'group' in pca_df.columns or 'treatment' in pca_df.columns:
        color_col = next((c for c in ['condition', 'group', 'treatment'] if c in pca_df.columns), None)
        fig = px.scatter(
            pca_df,
            x='PC1',
            y='PC2',
            color=color_col,
            text=pca_df.index,
            title="PCA Plot: Sample Clustering"
        )
    else:
        fig = px.scatter(
            pca_df,
            x='PC1',
            y='PC2',
            text=pca_df.index,
            title="PCA Plot: Sample Clustering"
        )
    
    fig.update_traces(
        textposition='top center',
        marker=dict(size=12)
    )
    
    fig.update_layout(
        xaxis_title=f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)",
        yaxis_title=f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)",
        template="plotly_white",
        height=600
    )
    
    return fig

def create_heatmap(counts, n_genes=50):
    """Create heatmap of top variable genes."""
    
    # Select top variable genes
    gene_var = counts.var(axis=1).sort_values(ascending=False)
    top_genes = gene_var.head(n_genes).index
    
    subset = counts.loc[top_genes]
    
    # Z-score normalize
    subset_norm = (subset - subset.mean(axis=1, keepdims=True)) / (subset.std(axis=1, keepdims=True) + 1e-8)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=subset_norm.values,
        x=subset_norm.columns,
        y=subset_norm.index,
        colorscale='RdBu_r',
        zmid=0,
        hovertemplate='Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>'
    ))
    
    fig.update_layout(
        title=f"Heatmap: Top {n_genes} Most Variable Genes",
        xaxis_title="Samples",
        yaxis_title="Genes",
        template="plotly_white",
        height=800
    )
    
    return fig

def create_expression_profile(counts, gene_id):
    """Create gene expression profile."""
    
    if gene_id not in counts.index:
        st.error(f"Gene {gene_id} not found in count matrix")
        return None
    
    expression = counts.loc[gene_id]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=expression.index,
        y=expression.values,
        marker_color='#2E7D32',
        hovertemplate='Sample: %{x}<br>Expression: %{y:.1f}<extra></extra>'
    ))
    
    fig.update_layout(
        title=f"Expression Profile: {gene_id}",
        xaxis_title="Samples",
        yaxis_title="Expression (Normalized Counts)",
        template="plotly_white",
        height=500
    )
    
    return fig

def create_correlation_plot(counts):
    """Create sample correlation heatmap."""
    
    # Calculate correlation
    corr_matrix = counts.T.corr()
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=corr_matrix.values,
        x=corr_matrix.columns,
        y=corr_matrix.index,
        colorscale='RdBu_r',
        zmid=corr_matrix.values.mean(),
        zmin=corr_matrix.values.min(),
        zmax=1.0,
        hovertemplate='Sample 1: %{y}<br>Sample 2: %{x}<br>Correlation: %{z:.3f}<extra></extra>'
    ))
    
    fig.update_layout(
        title="Sample-to-Sample Correlation",
        xaxis_title="Samples",
        yaxis_title="Samples",
        template="plotly_white",
        height=600
    )
    
    return fig

# =============================================================================
# PLOT RENDERING
# =============================================================================

if "Volcano" in plot_type and de_data is not None:
    st.markdown("## 🌋 Volcano Plot")
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        fig = create_volcano_plot(de_data, fdr_threshold, lfc_threshold, point_size)
        st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        st.markdown("### 📊 Statistics")
        
        sig_mask = (de_data['adjusted_p_value'] < fdr_threshold) & (abs(de_data['log2_fold_change']) >= lfc_threshold)
        sig_genes = de_data[sig_mask]
        up_genes = sig_genes[sig_genes['log2_fold_change'] > 0]
        down_genes = sig_genes[sig_genes['log2_fold_change'] < 0]
        
        st.metric("Total Genes", f"{len(de_data):,}")
        st.metric("Significant", f"{len(sig_genes):,}")
        st.metric("Upregulated", f"{len(up_genes):,}")
        st.metric("Downregulated", f"{len(down_genes):,}")
        
        pct = (len(sig_genes) / len(de_data) * 100) if len(de_data) > 0 else 0
        st.metric("% Significant", f"{pct:.1f}%")

elif "MA Plot" in plot_type and de_data is not None:
    st.markdown("## 📊 MA Plot")
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        fig = create_ma_plot(de_data, fdr_threshold, lfc_threshold, point_size)
        st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        st.markdown("### 📊 Statistics")
        
        sig_mask = (de_data['adjusted_p_value'] < fdr_threshold) & (abs(de_data['log2_fold_change']) >= lfc_threshold)
        sig_genes = de_data[sig_mask]
        
        st.metric("Total Genes", f"{len(de_data):,}")
        st.metric("Significant", f"{len(sig_genes):,}")
        st.metric("Median Expr", f"{de_data['base_mean'].median():.1f}")
        st.metric("Median LFC", f"{de_data['log2_fold_change'].median():.2f}")

elif "PCA" in plot_type and counts_data is not None:
    st.markdown("## 🎨 PCA Plot")
    
    fig = create_pca_plot(counts_data, metadata_data)
    st.plotly_chart(fig, use_container_width=True)
    
    st.info("💡 **Tip:** Upload metadata with 'condition', 'group', or 'treatment' column to color samples")

elif "Heatmap" in plot_type and counts_data is not None:
    st.markdown("## 🔥 Heatmap")
    
    n_genes = st.slider("Number of genes to display", 20, 100, 50, 10)
    
    fig = create_heatmap(counts_data, n_genes)
    st.plotly_chart(fig, use_container_width=True)

elif "Expression Profile" in plot_type and counts_data is not None:
    st.markdown("## 📈 Gene Expression Profile")
    
    gene_id = st.text_input(
        "Enter gene ID:",
        help="Type the exact gene ID from your count matrix"
    )
    
    if gene_id:
        fig = create_expression_profile(counts_data, gene_id)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("👆 Enter a gene ID to view its expression profile")
        
        # Show available genes
        with st.expander("📋 Available genes"):
            st.dataframe(pd.DataFrame({'Gene ID': counts_data.index}), height=300)

elif "Correlation" in plot_type and counts_data is not None:
    st.markdown("## 🔗 Sample Correlation")
    
    fig = create_correlation_plot(counts_data)
    st.plotly_chart(fig, use_container_width=True)
    
    st.caption("**Interpretation:** High correlation (red) = similar samples. Low correlation (blue) = different samples or potential outliers.")

# Footer
st.markdown("---")
st.caption("💡 **Tip:** All plots are interactive - zoom, pan, and hover for details")
st.caption("📊 **Export:** Right-click plot → 'Save image as' for quick export")
