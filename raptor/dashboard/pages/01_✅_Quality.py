"""
RAPTOR Dashboard - Quality Assessment (Module 2)

Comprehensive quality control and assessment for RNA-seq count data.
Uses RAPTOR's DataQualityAssessor with 6-method outlier detection.

Author: Ayeh Bolouki
Version: 2.2.0
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# Import RAPTOR quality assessment module
try:
    from raptor.quality_assessment import (
        DataQualityAssessor,
        quick_quality_check,
        OutlierResult
    )
    QUALITY_MODULE_AVAILABLE = True
except ImportError as e:
    QUALITY_MODULE_AVAILABLE = False
    import_error = str(e)

# Page config
st.set_page_config(
    page_title="RAPTOR - Quality Assessment",
    page_icon="✅",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("Quality Assessment — Module 2")
st.caption("Comprehensive quality control and outlier detection for RNA-seq count data")

# Help section
with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Perform comprehensive quality assessment on count data to identify and filter problematic samples.
    
    **Quality Metrics Calculated:**
    1. **Library Quality** - Sequencing depth, library size distribution
    2. **Gene Detection** - Detection rates, expression levels
    3. **Outlier Detection** - 6 advanced methods (PCA, Isolation Forest, LOF, etc.)
    4. **Batch Effects** - Statistical detection with F-tests
    5. **Variance Structure** - BCV, coefficient of variation
    6. **Biological Signal** - Signal-to-noise ratio
    7. **Overall Quality Score** - Integrated 0-100 score
    
    **Input:**
    - Count matrix (genes × samples, raw counts)
    - Sample metadata (optional, for batch detection)
    
    **Output:**
    - Comprehensive QC report with scores
    - Outlier samples identified (consensus of 6 methods)
    - Batch effect detection results
    - Interactive visualizations
    - Filtered sample recommendations
    
    **6 Outlier Detection Methods:**
    1. **PCA + Mahalanobis Distance** - Multivariate outliers in PC space
    2. **Isolation Forest** - Tree-based anomaly detection
    3. **Local Outlier Factor (LOF)** - Density-based detection
    4. **Elliptic Envelope** - Robust covariance estimation
    5. **Correlation-based** - Low correlation with other samples
    6. **Library Size Z-score** - Extreme library sizes
    
    **Workflow:**
    ```
    Upload Data → Run QC → Review Metrics → Check Outliers → Filter Samples → Continue to M3
    ```
    """)

# Check if RAPTOR module is available
if not QUALITY_MODULE_AVAILABLE:
    st.error(f"""
    ❌ **RAPTOR quality_assessment module not available**
    
    Error: {import_error}
    
    Please ensure RAPTOR is properly installed:
    ```bash
    cd RAPTOR/
    pip install -e .
    ```
    """)
    st.stop()

st.markdown("---")

# Section 1: Load Data
st.markdown("## 1. Load Count Data")

# Initialize variables
count_data = None
metadata = None
data_source = None

# Check if data is from Module 1 (not available in dashboard)
if st.session_state.get('m1_complete', False):
    st.success("✅ Using count data from Module 1 (Quick Count)")
    count_data = st.session_state.get('m1_counts')
    metadata = st.session_state.get('m1_sample_sheet')
    data_source = "module1"
else:
    st.info("ℹ️ Module 1 (Quick Count) is CLI-only. Please upload count matrix below.")

# Upload count matrix
with st.expander("Upload count matrix", expanded=(count_data is None)):
    uploaded_counts = st.file_uploader(
        "Upload count matrix (CSV/TSV)",
        type=['csv', 'tsv', 'txt'],
        help="Genes as rows, samples as columns. First column should be gene IDs."
    )
    
    if uploaded_counts:
        try:
            # Read file
            if uploaded_counts.name.endswith('.csv'):
                count_data = pd.read_csv(uploaded_counts, index_col=0)
            else:
                count_data = pd.read_csv(uploaded_counts, sep='\t', index_col=0)
            
            st.success(f"✅ Loaded {count_data.shape[0]} genes × {count_data.shape[1]} samples")
            data_source = "upload"
            
            # Update session state
            st.session_state['n_samples'] = count_data.shape[1]
            st.session_state['n_genes'] = count_data.shape[0]
            
        except Exception as e:
            st.error(f"❌ Error loading count matrix: {str(e)}")
            count_data = None

# Upload metadata (optional)
with st.expander("Upload metadata (optional)", expanded=False):
    st.markdown("""
    **Metadata helps detect:**
    - Batch effects
    - Sample grouping
    - Experimental conditions
    
    **Required columns:**
    - `sample_id` - Must match count matrix column names
    - `condition` or `group` - Sample groups
    - `batch` - Sequencing batch (optional)
    """)
    
    uploaded_metadata = st.file_uploader(
        "Upload sample metadata (CSV)",
        type=['csv'],
        help="Must include 'sample_id' column matching count matrix columns"
    )
    
    if uploaded_metadata:
        try:
            metadata = pd.read_csv(uploaded_metadata)
            
            # Validate
            if 'sample_id' not in metadata.columns:
                st.error("❌ Metadata must have 'sample_id' column")
                metadata = None
            else:
                st.success(f"✅ Loaded metadata for {len(metadata)} samples")
                
        except Exception as e:
            st.error(f"❌ Error loading metadata: {str(e)}")
            metadata = None

# Show preview
if count_data is not None:
    with st.expander("Preview data"):
        st.markdown("**Count Matrix Preview:**")
        st.dataframe(count_data.head(10), use_container_width=True)
        
        if metadata is not None:
            st.markdown("**Metadata Preview:**")
            st.dataframe(metadata.head(10), use_container_width=True)

st.markdown("---")

# Section 2: Configuration and Run
if count_data is not None:
    st.markdown("## 2. Configuration")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Normalization Method:**")
        normalization = st.selectbox(
            "Select normalization",
            options=['log2', 'cpm', 'quantile', 'none'],
            index=0,
            help="""
            - log2: log2(counts + 1) - Default, most common
            - cpm: Counts per million + log2 - Use if library sizes vary >5x
            - quantile: Forces same distribution - Aggressive normalization
            - none: No transformation - Use if already normalized (TPM, FPKM)
            """
        )
    
    with col2:
        st.markdown("**Outlier Detection:**")
        outlier_consensus = st.slider(
            "Consensus threshold (# methods)",
            min_value=1,
            max_value=6,
            value=3,
            help="How many of the 6 methods must agree to call a sample an outlier"
        )
    
    # Run button
    if st.button("Run Quality Assessment", type="primary", use_container_width=True):
        
        with st.spinner("Running comprehensive quality assessment..."):
            try:
                # Initialize assessor
                assessor = DataQualityAssessor(
                    counts=count_data,
                    metadata=metadata,
                    normalization=normalization
                )
                
                # Run full quality assessment
                st.session_state['m2_running'] = True
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                status_text.text("Assessing library quality...")
                progress_bar.progress(0.15)
                
                # Run assessment (this calls all internal methods)
                qc_report = assessor.assess_quality()
                
                status_text.text("Detecting outliers with 6 methods...")
                progress_bar.progress(0.50)
                
                # Run advanced outlier detection
                outlier_result = assessor.detect_outliers_advanced(
                    consensus_threshold=outlier_consensus
                )
                
                status_text.text("Finalizing results...")
                progress_bar.progress(0.90)
                
                # Store results in session state
                st.session_state['m2_qc_report'] = qc_report
                st.session_state['m2_outlier_result'] = outlier_result
                st.session_state['m2_assessor'] = assessor
                st.session_state['m2_count_data'] = count_data
                st.session_state['m2_metadata'] = metadata
                st.session_state['m2_complete'] = True
                st.session_state['m2_running'] = False
                st.session_state['m2_error'] = False
                
                progress_bar.progress(1.0)
                status_text.text("✅ Quality assessment complete!")
                
                st.success("""
                ✅ **Quality assessment complete!**
                
                Scroll down to see:
                - Overall quality score
                - Component scores
                - Outlier analysis
                - Batch effect detection
                - Visualizations
                """)
                
            except Exception as e:
                st.error(f"""
                ❌ **Error during quality assessment**
                
                {str(e)}
                
                **Troubleshooting:**
                - Ensure count matrix has samples as columns, genes as rows
                - Check that counts are non-negative integers
                - Verify metadata sample_id matches count matrix columns
                """)
                st.session_state['m2_error'] = True
                st.session_state['m2_running'] = False

# Section 3: Results
if st.session_state.get('m2_complete', False) and not st.session_state.get('m2_running', False):
    st.markdown("---")
    st.markdown("## Quality Assessment Results")
    
    qc_report = st.session_state['m2_qc_report']
    outlier_result = st.session_state['m2_outlier_result']
    
    # Professional styling for result cards
    st.markdown("""<style>
    .qc-card { 
        background: #f8f9fa; border: 1px solid #e0e0e0; border-radius: 8px; 
        padding: 16px 20px; text-align: center; 
    }
    .qc-card-label { 
        font-size: 0.8rem; color: #666; text-transform: uppercase; 
        letter-spacing: 0.5px; margin-bottom: 6px; 
    }
    .qc-card-value { font-size: 1.6rem; font-weight: 600; color: #333; }
    .qc-card-value.excellent { color: #2e7d32; }
    .qc-card-value.good { color: #f57f17; }
    .qc-card-value.poor { color: #c62828; }
    .qc-status { 
        display: inline-block; width: 10px; height: 10px; border-radius: 50%; 
        margin-right: 6px; vertical-align: middle; 
    }
    .qc-status.pass { background: #2e7d32; }
    .qc-status.warn { background: #f57f17; }
    .qc-status.fail { background: #c62828; }
    .comp-row { 
        display: flex; gap: 12px; margin-top: 12px; 
    }
    .comp-card { 
        flex: 1; background: #fff; border: 1px solid #e8e8e8; border-radius: 6px; 
        padding: 12px 16px; border-left: 4px solid #ccc; 
    }
    .comp-card.pass { border-left-color: #2e7d32; }
    .comp-card.warn { border-left-color: #f57f17; }
    .comp-card.fail { border-left-color: #c62828; }
    .comp-name { font-size: 0.75rem; color: #888; text-transform: uppercase; letter-spacing: 0.3px; }
    .comp-score { font-size: 1.25rem; font-weight: 600; color: #333; margin-top: 2px; }
    </style>""", unsafe_allow_html=True)
    
    # ---- Overall Quality Score ----
    st.markdown("#### Overall Quality Score")
    
    overall_score = qc_report['overall']['score']
    overall_status = qc_report['overall']['status']
    
    quality_map = {
        'good': ('Excellent', 'excellent'),
        'acceptable': ('Good', 'good'),
        'poor': ('Poor', 'poor')
    }
    quality_label, quality_class = quality_map.get(overall_status, ('Unknown', 'good'))
    
    # Issues count
    components_data = qc_report.get('components', {})
    issues_count = sum(1 for comp in components_data.values() 
                      if comp.get('status') in ['warning', 'poor'])
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown(f"""<div class="qc-card">
            <div class="qc-card-label">Quality Score</div>
            <div class="qc-card-value">{overall_score:.1f} / 100</div>
        </div>""", unsafe_allow_html=True)
    with col2:
        st.markdown(f"""<div class="qc-card">
            <div class="qc-card-label">Rating</div>
            <div class="qc-card-value {quality_class}">{quality_label}</div>
        </div>""", unsafe_allow_html=True)
    with col3:
        issues_text = "None" if issues_count == 0 else str(issues_count)
        issues_cls = "excellent" if issues_count == 0 else "poor"
        st.markdown(f"""<div class="qc-card">
            <div class="qc-card-label">Issues Found</div>
            <div class="qc-card-value {issues_cls}">{issues_text}</div>
        </div>""", unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # ---- Component Scores ----
    st.markdown("#### Component Scores")
    
    comp_items = {
        'Library Quality': components_data.get('library_quality', {}).get('score', 0),
        'Gene Detection': components_data.get('gene_detection', {}).get('score', 0),
        'Biological Signal': components_data.get('biological_signal', {}).get('score', 0),
        'Variance Structure': components_data.get('variance_structure', {}).get('score', 0)
    }
    
    comp_html = '<div class="comp-row">'
    for comp_name, comp_score in comp_items.items():
        cls = "pass" if comp_score >= 70 else "warn" if comp_score >= 50 else "fail"
        comp_html += f"""<div class="comp-card {cls}">
            <div class="comp-name">{comp_name}</div>
            <div class="comp-score">{comp_score:.1f}</div>
        </div>"""
    comp_html += '</div>'
    st.markdown(comp_html, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # ---- Outlier Detection ----
    st.markdown("#### Outlier Detection")
    
    st.markdown(
        f"**Consensus threshold:** {outlier_result.consensus_threshold} methods &nbsp;|&nbsp; "
        f"**Outliers detected:** {outlier_result.n_outliers} samples"
    )
    
    if outlier_result.n_outliers > 0:
        st.warning(
            f"**{outlier_result.n_outliers} outlier sample(s) detected:** "
            f"{', '.join(outlier_result.outlier_samples)}\n\n"
            f"{outlier_result.summary()}\n\n"
            "Review these samples before proceeding. Consider removing or investigating technical causes."
        )
        
        with st.expander("Detailed outlier breakdown"):
            outlier_data = []
            for sample in outlier_result.outlier_samples:
                methods = [method for method, samples_list in outlier_result.method_results.items()
                           if sample in samples_list]
                outlier_data.append({
                    'Sample': sample,
                    'Methods Flagged': len(methods),
                    'Methods': ', '.join(methods)
                })
            st.dataframe(pd.DataFrame(outlier_data), use_container_width=True)
    else:
        st.success("All samples pass quality checks. No outliers detected.")
    
    # ---- Batch Effect Detection ----
    batch_info = qc_report.get('components', {}).get('batch_effects', {})
    
    if batch_info.get('batch_detected', False):
        st.markdown("#### Batch Effect Detection")
        
        st.warning(
            f"**Batch effects detected** — "
            f"Strength: {batch_info.get('strength', 'Unknown')}, "
            f"F-statistic: {batch_info.get('f_statistic', 0):.2f}, "
            f"P-value: {batch_info.get('p_value', 1):.2e}\n\n"
            f"{batch_info.get('recommendation', 'Consider batch correction in DE analysis.')}"
        )
        
        if batch_info.get('confounded', False):
            st.error(
                "**Confounded design detected.** "
                "Batch is perfectly correlated with biological condition, "
                "making it impossible to separate batch effects from biological effects."
            )
    
    # ---- Library Quality Details ----
    with st.expander("Library Quality Details"):
        lib_qual = qc_report.get('components', {}).get('library_quality', {})
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Mean Library Size", f"{lib_qual.get('mean_size', 0):,.0f}")
            st.metric("Median Library Size", f"{lib_qual.get('median_size', 0):,.0f}")
        with col2:
            st.metric("CV", f"{lib_qual.get('cv', 0):.3f}")
            st.metric("Category", lib_qual.get('category', 'Unknown'))
        with col3:
            st.metric("Min Library Size", f"{lib_qual.get('min_size', 0):,.0f}")
            st.metric("Max Library Size", f"{lib_qual.get('max_size', 0):,.0f}")
    
    # ---- Gene Detection Details ----
    with st.expander("Gene Detection Details"):
        gene_det = qc_report.get('components', {}).get('gene_detection', {})
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Mean Detection Rate", f"{gene_det.get('mean_detection_rate', 0)*100:.1f}%")
            st.metric("Genes Detected (avg)", f"{gene_det.get('mean_detected', 0):,.0f}")
        with col2:
            st.metric("Well-Detected Genes", f"{gene_det.get('well_detected', 0):,.0f}")
            st.caption("Genes with mean count ≥ 10")
        with col3:
            st.metric("Detection Category", gene_det.get('category', 'Unknown'))
    
    # Visualizations
    st.markdown("### Visualizations")
    
    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
        "Library Sizes", "PCA Plot", "Sample Correlation",
        "Expression Distribution", "RLE Plot", "Dendrogram", "Mean-Variance"
    ])
    
    with tab1:
        # Library size distribution
        count_data = st.session_state['m2_count_data']
        library_sizes = count_data.sum(axis=0)
        
        fig = go.Figure()
        
        # Color by outlier status
        colors = ['red' if sample in outlier_result.outlier_samples else 'green' 
                  for sample in library_sizes.index]
        
        fig.add_trace(go.Bar(
            x=library_sizes.index,
            y=library_sizes.values,
            marker_color=colors,
            name='Library Size'
        ))
        
        fig.update_layout(
            title="Library Sizes by Sample",
            xaxis_title="Sample",
            yaxis_title="Total Reads",
            height=500,
            showlegend=False
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        # Enhanced PCA plot with metadata coloring, PC selection, and 3D option
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            
            # Get count data and apply log2 normalization
            pca_counts = st.session_state['m2_count_data']
            log_data = np.log2(pca_counts + 1)
            
            # Transpose to samples × genes, then scale
            data = log_data.T
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
            
            # Compute up to 4 PCs (or fewer if samples < 4)
            n_components = min(4, data_scaled.shape[0], data_scaled.shape[1])
            pca = PCA(n_components=n_components)
            pca_result = pca.fit_transform(data_scaled)
            
            # Build DataFrame with all PCs
            pc_cols = [f'PC{i+1}' for i in range(n_components)]
            pca_df = pd.DataFrame(pca_result, columns=pc_cols, index=data.index)
            
            # Variance explained labels
            var_labels = {f'PC{i+1}': f'PC{i+1} ({pca.explained_variance_ratio_[i]*100:.1f}%)'
                          for i in range(n_components)}
            
            # Add metadata columns for coloring
            meta = st.session_state.get('m2_metadata')
            color_options = ['Outlier Status']
            
            pca_df['Outlier Status'] = ['⚠️ Outlier' if s in outlier_result.outlier_samples 
                                        else '✓ Normal' for s in pca_df.index]
            
            if meta is not None and 'sample_id' in meta.columns:
                meta_indexed = meta.set_index('sample_id')
                for col in meta_indexed.columns:
                    if meta_indexed[col].nunique() <= 20:  # categorical-friendly
                        pca_df[col] = meta_indexed[col].reindex(pca_df.index)
                        color_options.append(col)
            
            # PCA controls
            st.markdown("**PCA Settings**")
            ctrl_col1, ctrl_col2, ctrl_col3, ctrl_col4 = st.columns(4)
            
            with ctrl_col1:
                color_by = st.selectbox("Color by", color_options,
                                        index=len(color_options)-1 if len(color_options) > 1 else 0)
            with ctrl_col2:
                available_pcs = pc_cols
                x_pc = st.selectbox("X axis", available_pcs, index=0)
            with ctrl_col3:
                y_pc = st.selectbox("Y axis", available_pcs, index=min(1, n_components - 1))
            with ctrl_col4:
                plot_mode = st.selectbox("Dimension", ["2D", "3D"] if n_components >= 3 else ["2D"])
            
            if plot_mode == "3D" and n_components >= 3:
                z_pc = st.selectbox("Z axis", available_pcs, index=min(2, n_components - 1))
            
            # Define a clean color palette
            palette = px.colors.qualitative.Set2
            
            if plot_mode == "2D":
                fig = px.scatter(
                    pca_df,
                    x=x_pc, y=y_pc,
                    color=color_by,
                    color_discrete_sequence=palette,
                    text=pca_df.index,
                    hover_data=[c for c in color_options if c != color_by],
                    title=f"PCA — {var_labels.get(x_pc, x_pc)} vs {var_labels.get(y_pc, y_pc)}"
                )
                fig.update_traces(
                    marker=dict(size=12, line=dict(width=1.5, color='white')),
                    textposition='top center',
                    textfont=dict(size=10)
                )
                fig.update_layout(
                    height=650,
                    xaxis_title=var_labels.get(x_pc, x_pc),
                    yaxis_title=var_labels.get(y_pc, y_pc),
                    legend_title_text=color_by,
                    plot_bgcolor='#fafafa',
                    font=dict(size=12)
                )
                fig.update_xaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
                fig.update_yaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
            else:
                fig = px.scatter_3d(
                    pca_df,
                    x=x_pc, y=y_pc, z=z_pc,
                    color=color_by,
                    color_discrete_sequence=palette,
                    text=pca_df.index,
                    hover_data=[c for c in color_options if c != color_by],
                    title=f"3D PCA — {var_labels.get(x_pc, x_pc)} vs {var_labels.get(y_pc, y_pc)} vs {var_labels.get(z_pc, z_pc)}"
                )
                fig.update_traces(
                    marker=dict(size=8, line=dict(width=0.5, color='white')),
                    textfont=dict(size=9)
                )
                fig.update_layout(
                    height=750,
                    scene=dict(
                        xaxis_title=var_labels.get(x_pc, x_pc),
                        yaxis_title=var_labels.get(y_pc, y_pc),
                        zaxis_title=var_labels.get(z_pc, z_pc)
                    ),
                    legend_title_text=color_by,
                    font=dict(size=12)
                )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Scree plot
            with st.expander("📊 Variance Explained (Scree Plot)"):
                scree_df = pd.DataFrame({
                    'Component': pc_cols,
                    'Variance (%)': pca.explained_variance_ratio_ * 100,
                    'Cumulative (%)': np.cumsum(pca.explained_variance_ratio_) * 100
                })
                fig_scree = go.Figure()
                fig_scree.add_trace(go.Bar(
                    x=scree_df['Component'], y=scree_df['Variance (%)'],
                    name='Individual', marker_color='#66c2a5'
                ))
                fig_scree.add_trace(go.Scatter(
                    x=scree_df['Component'], y=scree_df['Cumulative (%)'],
                    name='Cumulative', mode='lines+markers',
                    marker_color='#fc8d62', line=dict(width=2)
                ))
                fig_scree.update_layout(
                    title='Variance Explained per Component',
                    yaxis_title='Variance (%)', height=350,
                    plot_bgcolor='#fafafa', font=dict(size=12)
                )
                st.plotly_chart(fig_scree, use_container_width=True)
            
        except Exception as e:
            st.error(f"Could not generate PCA plot: {str(e)}")
    
    with tab3:
        # Enhanced sample correlation heatmap
        count_data = st.session_state['m2_count_data']
        
        # Log2-transform before correlation
        log_counts = np.log2(count_data + 1)
        correlation = log_counts.corr()
        
        # Correlation settings
        corr_col1, corr_col2 = st.columns(2)
        with corr_col1:
            corr_method = st.selectbox("Correlation method", 
                                       ["pearson", "spearman", "kendall"], index=0)
        with corr_col2:
            corr_colorscale = st.selectbox("Color scheme",
                                           ["Viridis", "Teal-Orange", "Purple-Green", "Blue-Red"],
                                           index=0)
        
        # Recompute if method changed
        if corr_method != "pearson":
            correlation = log_counts.corr(method=corr_method)
        
        # Color scale mapping
        colorscale_map = {
            "Viridis": "Viridis",
            "Teal-Orange": [[0, '#009392'], [0.5, '#f6edbd'], [1, '#d0587e']],
            "Purple-Green": [[0, '#7b3294'], [0.5, '#f7f7f7'], [1, '#008837']],
            "Blue-Red": "RdBu_r"
        }
        chosen_scale = colorscale_map[corr_colorscale]
        
        # Dynamic zmin: focus on the actual range for better contrast
        min_corr = correlation.values[~np.eye(len(correlation), dtype=bool)].min()
        zmin_val = max(0, np.floor(min_corr * 20) / 20)  # round down to nearest 0.05
        
        fig = go.Figure(data=go.Heatmap(
            z=correlation.values,
            x=correlation.columns,
            y=correlation.index,
            colorscale=chosen_scale,
            zmin=zmin_val,
            zmax=1.0,
            text=correlation.values,
            texttemplate='%{text:.2f}',
            textfont={"size": 9, "color": "white"},
            colorbar=dict(
                title=dict(text="Correlation", font=dict(size=12)),
                tickfont=dict(size=10)
            ),
            hovertemplate='%{x} vs %{y}<br>Correlation: %{z:.4f}<extra></extra>'
        ))
        
        fig.update_layout(
            title=f"Sample-to-Sample Correlation ({corr_method.capitalize()})",
            height=650,
            xaxis_title="Sample",
            yaxis_title="Sample",
            font=dict(size=12),
            xaxis=dict(tickangle=45)
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Quick stats
        off_diag = correlation.values[~np.eye(len(correlation), dtype=bool)]
        st.caption(
            f"Range: {off_diag.min():.3f} – {off_diag.max():.3f}  |  "
            f"Mean: {off_diag.mean():.3f}  |  "
            f"Std: {off_diag.std():.3f}"
        )
    
    # ================================================================
    # TAB 4: Expression Box/Violin Plot
    # ================================================================
    with tab4:
        st.markdown("**Per-sample expression distribution** — detects shifted or skewed samples")
        
        count_data = st.session_state['m2_count_data']
        log_counts = np.log2(count_data + 1)
        
        plot_type = st.radio("Plot type", ["Box Plot", "Violin Plot"], horizontal=True,
                             key="expr_dist_type")
        
        # Melt to long format for plotly
        melted = log_counts.melt(var_name='Sample', value_name='Log2(Count+1)')
        
        # Add metadata coloring if available
        meta = st.session_state.get('m2_metadata')
        if meta is not None and 'sample_id' in meta.columns and 'condition' in meta.columns:
            cond_map = meta.set_index('sample_id')['condition'].to_dict()
            melted['Condition'] = melted['Sample'].map(cond_map)
            color_col = 'Condition'
        else:
            color_col = None
        
        if plot_type == "Box Plot":
            fig = px.box(
                melted, x='Sample', y='Log2(Count+1)',
                color=color_col,
                color_discrete_sequence=px.colors.qualitative.Set2,
                title="Per-Sample Expression Distribution (Box Plot)"
            )
        else:
            fig = px.violin(
                melted, x='Sample', y='Log2(Count+1)',
                color=color_col,
                color_discrete_sequence=px.colors.qualitative.Set2,
                box=True, points=False,
                title="Per-Sample Expression Distribution (Violin Plot)"
            )
        
        fig.update_layout(
            height=550, xaxis_tickangle=45,
            plot_bgcolor='#fafafa', font=dict(size=12),
            yaxis_title="Log2(Count + 1)"
        )
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=True, gridcolor='#eee')
        st.plotly_chart(fig, use_container_width=True)
        
        st.caption(
            "**Interpretation:** All samples should have similar distributions. "
            "A sample with a noticeably shifted or compressed box suggests a library preparation "
            "or sequencing problem."
        )
    
    # ================================================================
    # TAB 5: RLE (Relative Log Expression) Plot
    # ================================================================
    with tab5:
        st.markdown(
            "**Relative Log Expression** — gold standard for detecting unwanted variation. "
            "Each box shows a sample's deviation from the gene-wise median."
        )
        
        count_data = st.session_state['m2_count_data']
        log_counts = np.log2(count_data + 1)
        
        # RLE = log expression - gene-wise median across samples
        gene_medians = log_counts.median(axis=1)
        rle = log_counts.subtract(gene_medians, axis=0)
        
        # Melt
        rle_melted = rle.melt(var_name='Sample', value_name='RLE')
        
        # Add metadata coloring
        meta = st.session_state.get('m2_metadata')
        if meta is not None and 'sample_id' in meta.columns and 'condition' in meta.columns:
            cond_map = meta.set_index('sample_id')['condition'].to_dict()
            rle_melted['Condition'] = rle_melted['Sample'].map(cond_map)
            color_col = 'Condition'
        else:
            color_col = None
        
        fig = px.box(
            rle_melted, x='Sample', y='RLE',
            color=color_col,
            color_discrete_sequence=px.colors.qualitative.Set2,
            title="Relative Log Expression (RLE) Plot"
        )
        
        # Add reference line at 0
        fig.add_hline(y=0, line_dash="dash", line_color="red", line_width=1.5,
                      annotation_text="Median = 0", annotation_position="top left")
        
        fig.update_layout(
            height=550, xaxis_tickangle=45,
            plot_bgcolor='#fafafa', font=dict(size=12),
            yaxis_title="Relative Log Expression"
        )
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=True, gridcolor='#eee')
        st.plotly_chart(fig, use_container_width=True)
        
        # RLE summary stats
        rle_medians = rle.median(axis=0)
        rle_iqrs = rle.quantile(0.75, axis=0) - rle.quantile(0.25, axis=0)
        
        st.caption(
            "**Interpretation:** Well-normalized data has boxes centered at 0 with similar IQR. "
            "Boxes shifted away from 0 indicate systematic bias. Wide boxes indicate high variability."
        )
        
        # Flag problematic samples
        shifted = rle_medians[rle_medians.abs() > 0.5]
        wide = rle_iqrs[rle_iqrs > rle_iqrs.median() * 2]
        
        if len(shifted) > 0 or len(wide) > 0:
            with st.expander("⚠️ Flagged samples"):
                if len(shifted) > 0:
                    st.warning(f"**Shifted median (|median| > 0.5):** {', '.join(shifted.index)}")
                if len(wide) > 0:
                    st.warning(f"**High variability (IQR > 2× median IQR):** {', '.join(wide.index)}")
    
    # ================================================================
    # TAB 6: Sample Dendrogram (Hierarchical Clustering)
    # ================================================================
    with tab6:
        st.markdown(
            "**Hierarchical clustering** — shows which samples are most similar. "
            "Samples should cluster by biological condition, not by batch."
        )
        
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import pdist
        
        count_data = st.session_state['m2_count_data']
        log_counts = np.log2(count_data + 1).T  # samples × genes
        
        dendro_col1, dendro_col2 = st.columns(2)
        with dendro_col1:
            dist_metric = st.selectbox("Distance metric",
                                       ["euclidean", "correlation", "cosine"],
                                       index=1, key="dendro_dist")
        with dendro_col2:
            link_method = st.selectbox("Linkage method",
                                       ["ward", "complete", "average", "single"],
                                       index=0, key="dendro_link")
        
        # Compute distances
        if dist_metric == "correlation":
            # 1 - Pearson correlation as distance
            dists = pdist(log_counts.values, metric='correlation')
        else:
            dists = pdist(log_counts.values, metric=dist_metric)
        
        # Ward requires euclidean
        if link_method == "ward" and dist_metric != "euclidean":
            link_method_used = "complete"
            st.caption("ℹ️ Ward linkage requires Euclidean distance — using 'complete' instead.")
        else:
            link_method_used = link_method
        
        Z = linkage(dists, method=link_method_used)
        
        # Compute dendrogram data (scipy returns dict with plotting info)
        import matplotlib
        matplotlib.use('Agg')  # non-interactive backend
        import matplotlib.pyplot as plt
        
        fig_mpl, ax = plt.subplots(figsize=(12, 5))
        dendro_data = dendrogram(
            Z, labels=log_counts.index.tolist(), ax=ax,
            leaf_rotation=45, leaf_font_size=10
        )
        ax.set_ylabel("Distance", fontsize=12)
        ax.set_title("Sample Dendrogram", fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Color the x-tick labels by condition if metadata available
        meta = st.session_state.get('m2_metadata')
        if meta is not None and 'sample_id' in meta.columns and 'condition' in meta.columns:
            cond_map = meta.set_index('sample_id')['condition'].to_dict()
            palette_map = {}
            unique_conds = sorted(set(cond_map.values()))
            colors_list = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
            for i, cond in enumerate(unique_conds):
                palette_map[cond] = colors_list[i % len(colors_list)]
            
            for lbl in ax.get_xticklabels():
                sample = lbl.get_text()
                cond = cond_map.get(sample, None)
                if cond and cond in palette_map:
                    lbl.set_color(palette_map[cond])
                    lbl.set_fontweight('bold')
            
            # Add legend
            from matplotlib.patches import Patch
            legend_handles = [Patch(facecolor=palette_map[c], label=c) for c in unique_conds]
            ax.legend(handles=legend_handles, loc='upper right', fontsize=9)
        
        plt.tight_layout()
        st.pyplot(fig_mpl)
        plt.close(fig_mpl)
        
        st.caption(
            "**Interpretation:** Samples from the same condition should cluster together. "
            "If batches dominate the clustering instead of biology, batch correction is needed."
        )
    
    # ================================================================
    # TAB 7: Mean-Variance (BCV) Plot
    # ================================================================
    with tab7:
        st.markdown(
            "**Mean-Variance relationship** — shows how variance scales with expression level. "
            "The Biological Coefficient of Variation (BCV) is the square root of dispersion."
        )
        
        count_data = st.session_state['m2_count_data']
        
        # Calculate per-gene mean and variance on raw counts
        gene_means = count_data.mean(axis=1)
        gene_vars = count_data.var(axis=1)
        
        # Filter out zero-mean genes
        mask = gene_means > 0
        gene_means = gene_means[mask]
        gene_vars = gene_vars[mask]
        
        # BCV = sqrt(variance / mean^2) — biological coefficient of variation
        bcv = np.sqrt(gene_vars / (gene_means ** 2))
        
        # Log-transform for plotting
        log_mean = np.log10(gene_means)
        log_var = np.log10(gene_vars)
        
        mv_plot = st.radio("View", ["Mean-Variance", "BCV vs Mean"], horizontal=True,
                           key="mv_view")
        
        if mv_plot == "Mean-Variance":
            fig = go.Figure()
            
            fig.add_trace(go.Scattergl(
                x=log_mean.values, y=log_var.values,
                mode='markers',
                marker=dict(size=3, color='#3288bd', opacity=0.4),
                name='Genes',
                hovertemplate='Mean: 10^%{x:.1f}<br>Var: 10^%{y:.1f}<extra></extra>'
            ))
            
            # Add Poisson line (variance = mean) for reference
            x_range = np.linspace(log_mean.min(), log_mean.max(), 100)
            fig.add_trace(go.Scatter(
                x=x_range, y=x_range,
                mode='lines', line=dict(color='red', dash='dash', width=1.5),
                name='Poisson (var = mean)'
            ))
            
            fig.update_layout(
                title="Mean-Variance Relationship",
                xaxis_title="Log10(Mean Count)",
                yaxis_title="Log10(Variance)",
                height=600, plot_bgcolor='#fafafa', font=dict(size=12),
                legend=dict(x=0.02, y=0.98)
            )
            fig.update_xaxes(showgrid=True, gridcolor='#eee')
            fig.update_yaxes(showgrid=True, gridcolor='#eee')
            
        else:  # BCV vs Mean
            fig = go.Figure()
            
            fig.add_trace(go.Scattergl(
                x=log_mean.values, y=bcv.values,
                mode='markers',
                marker=dict(size=3, color='#d53e4f', opacity=0.4),
                name='Gene BCV',
                hovertemplate='Mean: 10^%{x:.1f}<br>BCV: %{y:.2f}<extra></extra>'
            ))
            
            # Add median BCV line
            median_bcv = bcv.median()
            fig.add_hline(y=median_bcv, line_dash="dash", line_color="#666", line_width=1,
                          annotation_text=f"Median BCV = {median_bcv:.2f}",
                          annotation_position="top right")
            
            fig.update_layout(
                title="Biological Coefficient of Variation (BCV) vs Mean Expression",
                xaxis_title="Log10(Mean Count)",
                yaxis_title="BCV (√dispersion)",
                height=600, plot_bgcolor='#fafafa', font=dict(size=12)
            )
            fig.update_xaxes(showgrid=True, gridcolor='#eee')
            fig.update_yaxes(showgrid=True, gridcolor='#eee')
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Summary stats
        st.caption(
            f"**Summary:** {len(gene_means):,} genes plotted  |  "
            f"Median BCV: {bcv.median():.3f}  |  "
            f"Mean BCV: {bcv.mean():.3f}  |  "
            f"BCV range: {bcv.min():.3f} – {bcv.quantile(0.99):.3f} (99th percentile)"
        )
        
        with st.expander("ℹ️ How to interpret"):
            st.markdown("""
            **Mean-Variance plot:**
            - Points above the red Poisson line show **overdispersion** (expected in RNA-seq)
            - More overdispersion at low counts is normal
            - If all points sit on the Poisson line, the data may be technical replicates
            
            **BCV plot:**
            - BCV < 0.1 → very low variation (technical replicates or very homogeneous)
            - BCV 0.1–0.4 → typical for well-controlled experiments
            - BCV 0.4–1.0 → high variation (outbred organisms, clinical samples)
            - BCV > 1.0 → very noisy data (check for technical issues)
            
            Low-expression genes (left side) typically show higher BCV due to count noise.
            """)
    
    # Export options
    st.markdown("### Export Results")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # Export filtered count matrix (outliers removed)
        if outlier_result.n_outliers > 0:
            filtered_counts = count_data.drop(columns=outlier_result.outlier_samples)
            
            csv = filtered_counts.to_csv()
            st.download_button(
                "📥 Download Filtered Counts (outliers removed)",
                csv,
                "filtered_counts.csv",
                "text/csv",
                help=f"Count matrix with {outlier_result.n_outliers} outliers removed"
            )
        else:
            st.info("No outliers to filter")
    
    with col2:
        # Export QC report
        import json
        
        report_json = json.dumps(qc_report, indent=2, default=str)
        st.download_button(
            "📥 Download QC Report (JSON)",
            report_json,
            "qc_report.json",
            "application/json"
        )
    
    # Next steps
    st.markdown("---")
    st.info("""
    ✅ **Quality assessment complete!**
    
    **Next steps:**
    1. Review outliers and decide whether to remove them
    2. Note any batch effects for downstream analysis
    3. Proceed to **Module 3: Data Profiler** to characterize your data
    4. Use cleaned data for differential expression analysis
    """)

else:
    if not st.session_state.get('m2_complete', False):
        st.info("Upload count data and run quality assessment to see results")
