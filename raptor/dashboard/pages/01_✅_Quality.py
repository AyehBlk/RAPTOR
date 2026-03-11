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
st.title("✅ Quality Assessment - Module 2")
st.markdown("### Comprehensive quality control for RNA-seq data")

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
st.markdown("## 1️⃣ Load Count Data")

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
with st.expander("📤 Upload count matrix", expanded=(count_data is None)):
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
with st.expander("📋 Upload metadata (optional)", expanded=False):
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
    with st.expander("👁️ Preview data"):
        st.markdown("**Count Matrix Preview:**")
        st.dataframe(count_data.head(10), use_container_width=True)
        
        if metadata is not None:
            st.markdown("**Metadata Preview:**")
            st.dataframe(metadata.head(10), use_container_width=True)

st.markdown("---")

# Section 2: Configuration and Run
if count_data is not None:
    st.markdown("## 2️⃣ Quality Assessment Configuration")
    
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
    if st.button("🚀 Run Quality Assessment", type="primary", use_container_width=True):
        
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
    st.markdown("## 3️⃣ Quality Assessment Results")
    
    qc_report = st.session_state['m2_qc_report']
    outlier_result = st.session_state['m2_outlier_result']
    
    # Overall Quality Score
    st.markdown("### 🎯 Overall Quality Score")
    
    overall_score = qc_report['overall']['score']
    overall_status = qc_report['overall']['status']
    
    # Map status to quality rating
    quality_map = {
        'good': ('EXCELLENT', '🟢'),
        'acceptable': ('GOOD', '🟡'),
        'poor': ('POOR', '🔴')
    }
    quality_rating, quality_icon = quality_map.get(overall_status, ('UNKNOWN', '⚪'))
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Quality Score", f"{overall_score:.1f}/100")
    
    with col2:
        st.metric("Quality Rating", f"{quality_icon} {quality_rating}")
    
    with col3:
        # Check components for issues
        components = qc_report.get('components', {})
        issues_count = sum(1 for comp in components.values() 
                          if comp.get('status') in ['warning', 'poor'])
        if issues_count == 0:
            st.metric("Issues Found", "✅ None")
        else:
            st.metric("Issues Found", f"⚠️ {issues_count}")
    
    # Component Scores
    st.markdown("### 📊 Component Scores")
    
    components_data = qc_report.get('components', {})
    
    components = {
        'Library Quality': components_data.get('library_quality', {}).get('score', 0),
        'Gene Detection': components_data.get('gene_detection', {}).get('score', 0),
        'Biological Signal': components_data.get('biological_signal', {}).get('score', 0),
        'Variance Structure': components_data.get('variance_structure', {}).get('score', 0)
    }
    
    col1, col2, col3, col4 = st.columns(4)
    
    for col, (comp_name, comp_score) in zip([col1, col2, col3, col4], components.items()):
        with col:
            status = "✅" if comp_score >= 70 else "⚠️" if comp_score >= 50 else "❌"
            col.metric(comp_name, f"{status} {comp_score:.1f}/100")
    
    # Outlier Detection Results
    st.markdown("### 🔍 Outlier Detection Results")
    
    st.markdown(f"""
    **Consensus Threshold:** {outlier_result.consensus_threshold} methods
    **Outliers Detected:** {outlier_result.n_outliers} samples
    """)
    
    if outlier_result.n_outliers > 0:
        st.warning(f"""
        ⚠️ **{outlier_result.n_outliers} outlier sample(s) detected:**
        
        {', '.join(outlier_result.outlier_samples)}
        
        **Detection Summary:**
        {outlier_result.summary()}
        
        **Recommendation:** Review these samples carefully before proceeding to analysis.
        Consider removing outliers or investigating technical issues.
        """)
        
        # Show which methods flagged each outlier
        with st.expander("📋 Detailed outlier detection breakdown"):
            outlier_data = []
            for sample in outlier_result.outlier_samples:
                methods = []
                for method, samples_list in outlier_result.method_results.items():
                    if sample in samples_list:
                        methods.append(method)
                outlier_data.append({
                    'Sample': sample,
                    'Methods Flagged': len(methods),
                    'Methods': ', '.join(methods)
                })
            
            outlier_df = pd.DataFrame(outlier_data)
            st.dataframe(outlier_df, use_container_width=True)
    
    else:
        st.success("✅ No outliers detected! All samples pass quality checks.")
    
    # Batch Effect Detection
    batch_info = qc_report.get('components', {}).get('batch_effects', {})
    
    if batch_info.get('batch_detected', False):
        st.markdown("### 🔬 Batch Effect Detection")
        
        st.warning(f"""
        ⚠️ **Batch effects detected**
        
        **Strength:** {batch_info.get('strength', 'Unknown')}
        **F-statistic:** {batch_info.get('f_statistic', 0):.2f}
        **P-value:** {batch_info.get('p_value', 1):.2e}
        
        **Recommendation:**
        {batch_info.get('recommendation', 'Consider batch correction in DE analysis')}
        """)
        
        if batch_info.get('confounded', False):
            st.error("""
            🚨 **Confounded Design Detected**
            
            Batch is perfectly correlated with biological condition.
            This makes it impossible to separate batch effects from biological effects.
            """)
    
    # Library Quality Details
    with st.expander("📚 Library Quality Details"):
        lib_qual = qc_report.get('components', {}).get('library_quality', {})
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Mean Library Size", f"{lib_qual.get('mean_size', 0):,.0f} reads")
            st.metric("Median Library Size", f"{lib_qual.get('median_size', 0):,.0f} reads")
        
        with col2:
            st.metric("CV (Coefficient of Variation)", f"{lib_qual.get('cv', 0):.3f}")
            st.metric("Category", lib_qual.get('category', 'Unknown'))
        
        with col3:
            st.metric("Min Library Size", f"{lib_qual.get('min_size', 0):,.0f} reads")
            st.metric("Max Library Size", f"{lib_qual.get('max_size', 0):,.0f} reads")
    
    # Gene Detection Details
    with st.expander("🧬 Gene Detection Details"):
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
    st.markdown("### 📈 Interactive Visualizations")
    
    tab1, tab2, tab3 = st.tabs(["📊 Library Sizes", "🎯 PCA Plot", "🔗 Sample Correlation"])
    
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
        # PCA plot from assessor
        assessor = st.session_state['m2_assessor']
        
        try:
            # Perform PCA on normalized data
            from sklearn.decomposition import PCA
            
            # Get normalized data
            data = assessor.data.T  # Transpose to samples x genes
            
            # Run PCA
            pca = PCA(n_components=2)
            pca_result = pca.fit_transform(data)
            
            # Create DataFrame
            pca_df = pd.DataFrame(
                pca_result,
                columns=['PC1', 'PC2'],
                index=data.index
            )
            
            # Add outlier status
            pca_df['Outlier'] = pca_df.index.isin(outlier_result.outlier_samples)
            
            # Plot
            fig = px.scatter(
                pca_df,
                x='PC1',
                y='PC2',
                color='Outlier',
                color_discrete_map={True: 'red', False: 'green'},
                text=pca_df.index,
                title=f"PCA Plot (PC1: {pca.explained_variance_ratio_[0]*100:.1f}%, PC2: {pca.explained_variance_ratio_[1]*100:.1f}%)"
            )
            
            fig.update_traces(textposition='top center')
            fig.update_layout(height=600)
            
            st.plotly_chart(fig, use_container_width=True)
            
        except Exception as e:
            st.error(f"Could not generate PCA plot: {str(e)}")
    
    with tab3:
        # Sample correlation heatmap
        count_data = st.session_state['m2_count_data']
        
        # Calculate correlation
        correlation = count_data.corr()
        
        fig = go.Figure(data=go.Heatmap(
            z=correlation.values,
            x=correlation.columns,
            y=correlation.index,
            colorscale='RdBu_r',
            zmid=0,
            text=correlation.values,
            texttemplate='%{text:.2f}',
            textfont={"size": 8},
            colorbar=dict(title="Correlation")
        ))
        
        fig.update_layout(
            title="Sample-to-Sample Correlation",
            height=600,
            xaxis_title="Sample",
            yaxis_title="Sample"
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    # Export options
    st.markdown("### 💾 Export Results")
    
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
        st.info("👆 Upload count data and run quality assessment to see results")
