"""
RAPTOR Dashboard - Data Profiler (Module 3)

Comprehensive data profiling and exploratory analysis for RNA-seq data.
Uses RAPTOR's RNAseqDataProfiler with 32-feature extraction.

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

# Import RAPTOR profiler module
try:
    from raptor.profiler import (
        RNAseqDataProfiler,
        DataProfile,
        profile_data_quick,
        get_key_characteristics
    )
    PROFILER_AVAILABLE = True
except ImportError as e:
    PROFILER_AVAILABLE = False
    import_error = str(e)

# Page config
st.set_page_config(
    page_title="RAPTOR - Data Profiler",
    page_icon="📊",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("📊 Data Profiler - Module 3")
st.markdown("### Comprehensive exploratory data analysis")

# Help section
with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Extract 32 statistical features from RNA-seq data to characterize your dataset and guide pipeline selection.
    
    **What this module analyzes:**
    
    **📊 Sample Characteristics (8 features):**
    - Sample sizes and group balance
    - Replicate sufficiency
    - Design complexity
    
    **📚 Library Metrics (8 features):**
    - Library size statistics (mean, median, CV)
    - Sequencing depth category
    - Library size variability
    
    **🧬 Gene Detection (3 features):**
    - Expressed genes count
    - Reliably expressed genes (mean count ≥ 10)
    - Detection rates
    
    **📈 Expression Distribution (7 features):**
    - Mean, median, variance on log scale
    - Skewness and kurtosis
    - Expression range and IQR
    
    **🎯 Dispersion Estimates (4 features):**
    - **BCV (Biological Coefficient of Variation)** - CRITICAL for pipeline choice
    - Common dispersion
    - Trended dispersion
    - Overdispersion ratio
    
    **🕳️ Sparsity Metrics (2 features):**
    - Zero proportion
    - Zero inflation index
    
    **Input:**
    - Count matrix (genes × samples, raw counts)
    - Sample metadata with 'condition' or 'group' column (optional)
    
    **Output:**
    - Complete DataProfile with 32 features
    - Interactive visualizations
    - Recommendations for Module 4 (pipeline selection)
    - JSON export for documentation
    
    **Workflow:**
    ```
    Quality-Filtered Data (M2) → Profile Extraction → Feature Analysis → M4 (Recommender)
    ```
    """)

# Check if RAPTOR module is available
if not PROFILER_AVAILABLE:
    st.error(f"""
    ❌ **RAPTOR profiler module not available**
    
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
st.markdown("## 1️⃣ Load Data")

# Initialize variables
count_data = None
metadata = None
data_source = None

# Check for data from previous modules
if st.session_state.get('m2_complete', False):
    st.success("✅ Using quality-filtered data from Module 2")
    count_data = st.session_state.get('m2_count_data')
    metadata = st.session_state.get('m2_metadata')
    data_source = "module2"
    
    # Update session state
    if count_data is not None:
        st.session_state['n_samples'] = count_data.shape[1]
        st.session_state['n_genes'] = count_data.shape[0]
    
elif st.session_state.get('m1_complete', False):
    st.info("ℹ️ Using data from Module 1 (Quality check skipped)")
    count_data = st.session_state.get('m1_counts')
    metadata = st.session_state.get('m1_sample_sheet')
    data_source = "module1"
else:
    st.info("ℹ️ Modules 1 & 2 are optional. Please upload count matrix below.")

# Upload count matrix
with st.expander("📤 Upload count matrix", expanded=(count_data is None)):
    uploaded_counts = st.file_uploader(
        "Upload count matrix (CSV/TSV)",
        type=['csv', 'tsv', 'txt'],
        help="Genes as rows, samples as columns. Raw counts (integers) required."
    )
    
    if uploaded_counts:
        try:
            if uploaded_counts.name.endswith('.csv'):
                count_data = pd.read_csv(uploaded_counts, index_col=0)
            else:
                count_data = pd.read_csv(uploaded_counts, sep='\t', index_col=0)
            
            st.success(f"✅ Loaded {count_data.shape[0]} genes × {count_data.shape[1]} samples")
            data_source = "upload"
            
            st.session_state['n_samples'] = count_data.shape[1]
            st.session_state['n_genes'] = count_data.shape[0]
            
        except Exception as e:
            st.error(f"❌ Error loading count matrix: {str(e)}")
            count_data = None

# Upload metadata (optional)
with st.expander("📋 Upload metadata (optional)", expanded=False):
    st.markdown("""
    **Metadata is used to:**
    - Calculate group sizes and balance
    - Assess replicate sufficiency
    - Determine design complexity
    
    **Required columns:**
    - `sample` or `sample_id` - Sample names (must match count matrix columns)
    - `condition` or `group` - Experimental groups
    - `batch` - Batch information (optional)
    """)
    
    uploaded_metadata = st.file_uploader(
        "Upload sample metadata (CSV)",
        type=['csv'],
        help="Must include 'sample' and 'condition' or 'group' columns"
    )
    
    if uploaded_metadata:
        try:
            metadata = pd.read_csv(uploaded_metadata)
            
            # Validate required columns
            if 'sample' not in metadata.columns and 'sample_id' not in metadata.columns:
                st.error("❌ Metadata must have 'sample' or 'sample_id' column")
                metadata = None
            elif 'condition' not in metadata.columns and 'group' not in metadata.columns:
                st.warning("⚠️ No 'condition' or 'group' column found. Group analysis will be limited.")
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
        
        # Quick stats
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Genes", f"{count_data.shape[0]:,}")
        
        with col2:
            st.metric("Samples", count_data.shape[1])
        
        with col3:
            total_counts = count_data.sum().sum()
            st.metric("Total Counts", f"{total_counts:,.0f}")
        
        with col4:
            sparsity = (count_data == 0).sum().sum() / count_data.size * 100
            st.metric("Sparsity", f"{sparsity:.1f}%")

st.markdown("---")

# Section 2: Configuration and Run
if count_data is not None:
    st.markdown("## 2️⃣ Profiling Configuration")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Group Column:**")
        
        # Determine available group columns
        group_options = ['condition']  # Default
        if metadata is not None:
            available_cols = [col for col in metadata.columns 
                            if col in ['condition', 'group', 'treatment', 'type']]
            if available_cols:
                group_options = available_cols
        
        group_column = st.selectbox(
            "Select group/condition column",
            options=group_options,
            help="Column in metadata that defines experimental groups"
        )
    
    with col2:
        st.markdown("**Detection Threshold:**")
        min_count_threshold = st.number_input(
            "Minimum count for 'expressed'",
            min_value=0,
            max_value=100,
            value=1,
            help="Minimum count to consider a gene expressed (default: 1)"
        )
    
    # Run button
    if st.button("🚀 Run Data Profiling", type="primary", use_container_width=True):
        
        with st.spinner("Running comprehensive data profiling..."):
            try:
                # Initialize profiler
                profiler = RNAseqDataProfiler(
                    counts=count_data,
                    metadata=metadata,
                    group_column=group_column,
                    min_count_threshold=min_count_threshold
                )
                
                st.session_state['m3_running'] = True
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                status_text.text("Profiling samples and groups...")
                progress_bar.progress(0.15)
                
                # Run full profiling (8 different analyses)
                profile = profiler.run_full_profile()
                
                status_text.text("Extracting statistical features...")
                progress_bar.progress(0.80)
                
                # Get recommendation features
                features = profile.get_recommendation_features()
                
                status_text.text("Finalizing profile...")
                progress_bar.progress(0.95)
                
                # Store results
                st.session_state['m3_profile'] = profile
                st.session_state['m3_features'] = features
                st.session_state['m3_profiler'] = profiler
                st.session_state['m3_count_data'] = count_data
                st.session_state['m3_metadata'] = metadata
                st.session_state['m3_complete'] = True
                st.session_state['m3_running'] = False
                st.session_state['m3_error'] = False
                
                progress_bar.progress(1.0)
                status_text.text("✅ Data profiling complete!")
                
                st.success("""
                ✅ **Data profiling complete!**
                
                Extracted 32 statistical features:
                - 8 sample characteristics
                - 8 library metrics
                - 3 gene detection features
                - 7 expression distribution features
                - 4 dispersion estimates
                - 2 sparsity metrics
                
                Scroll down to see results!
                """)
                
            except Exception as e:
                st.error(f"""
                ❌ **Error during profiling**
                
                {str(e)}
                
                **Troubleshooting:**
                - Ensure count matrix has samples as columns, genes as rows
                - Counts must be non-negative integers
                - If using metadata, ensure group column exists
                - Need at least 2 samples for profiling
                """)
                st.session_state['m3_error'] = True
                st.session_state['m3_running'] = False

# Section 3: Results
if st.session_state.get('m3_complete', False) and not st.session_state.get('m3_running', False):
    st.markdown("---")
    st.markdown("## 3️⃣ Data Profile Results")
    
    profile = st.session_state['m3_profile']
    features = st.session_state['m3_features']
    
    # Overall Summary
    st.markdown("### 📋 Profile Summary")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Samples", profile.n_samples)
        st.metric("Total Genes", f"{profile.n_genes:,}")
    
    with col2:
        st.metric("Expressed Genes", f"{profile.n_expressed_genes:,}")
        st.metric("Reliably Expressed", f"{profile.n_reliably_expressed:,}")
        st.caption("Mean count ≥ 10")
    
    with col3:
        st.metric("Library Size (Mean)", f"{profile.library_size_mean:,.0f}")
        st.metric("Sequencing Depth", profile.sequencing_depth_category.upper())
    
    with col4:
        st.metric("Sparsity", f"{profile.sparsity*100:.1f}%")
        st.metric("Sparsity Category", profile.sparsity_category.upper())
    
    # Key Features for Pipeline Selection
    st.markdown("### 🎯 Key Features for Pipeline Selection")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("#### Biological Coefficient of Variation (BCV)")
        bcv_status = "✅" if profile.bcv <= 0.4 else "⚠️" if profile.bcv <= 0.6 else "❌"
        st.metric("BCV", f"{bcv_status} {profile.bcv:.3f}")
        st.caption(f"Category: {profile.bcv_category.upper()}")
        
        st.markdown("""
        **BCV Categories:**
        - Low < 0.2: Clean data, simple pipelines work
        - Moderate 0.2-0.4: Standard pipelines
        - High > 0.4: Need robust methods
        """)
    
    with col2:
        st.markdown("#### Sample Design")
        st.metric("Number of Groups", profile.n_groups)
        st.metric("Min Group Size", profile.min_group_size)
        st.metric("Sample Balance", f"{profile.sample_balance:.2f}")
        
        replicates_status = "✅" if profile.has_sufficient_replicates else "⚠️"
        st.metric("Sufficient Replicates", f"{replicates_status} {profile.has_sufficient_replicates}")
        st.caption("≥3 per group recommended")
    
    with col3:
        st.markdown("#### Overdispersion")
        st.metric("Common Dispersion", f"{profile.common_dispersion:.3f}")
        st.metric("Overdispersion Ratio", f"{profile.overdispersion_ratio:.2f}")
        st.caption(f"Category: {profile.overdispersion_category.upper()}")
    
    # Sample Characteristics
    with st.expander("📊 Detailed Sample Characteristics"):
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Group Information:**")
            st.write(f"- Groups: {profile.n_groups}")
            st.write(f"- Group sizes: {profile.group_sizes}")
            st.write(f"- Min/Max group size: {profile.min_group_size}/{profile.max_group_size}")
            st.write(f"- Mean group size: {profile.mean_group_size:.1f}")
            st.write(f"- Balance ratio: {profile.sample_balance:.2f}")
        
        with col2:
            st.markdown("**Design Assessment:**")
            st.write(f"- Has replicates (≥2): {'✅ Yes' if profile.has_replicates else '❌ No'}")
            st.write(f"- Sufficient replicates (≥3): {'✅ Yes' if profile.has_sufficient_replicates else '⚠️ No'}")
            st.write(f"- Design complexity: {profile.design_complexity.upper()}")
    
    # Library Metrics
    with st.expander("📚 Detailed Library Metrics"):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Mean Library Size", f"{profile.library_size_mean:,.0f}")
            st.metric("Median Library Size", f"{profile.library_size_median:,.0f}")
            st.metric("SD Library Size", f"{profile.library_size_std:,.0f}")
        
        with col2:
            st.metric("CV (Coefficient of Variation)", f"{profile.library_size_cv:.3f}")
            st.metric("Library Size Range", f"{profile.library_size_range:.2f}x")
            st.caption("Max/Min ratio")
        
        with col3:
            st.metric("IQR", f"{profile.library_size_iqr:,.0f}")
            st.metric("Category", profile.library_size_category.upper())
            st.metric("Depth Category", profile.sequencing_depth_category.upper())
    
    # Expression Distribution
    with st.expander("📈 Expression Distribution Details"):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Mean Expression (log2)", f"{profile.expression_mean:.2f}")
            st.metric("Median Expression (log2)", f"{profile.expression_median:.2f}")
            st.metric("SD Expression", f"{profile.expression_std:.2f}")
        
        with col2:
            st.metric("Skewness", f"{profile.expression_skewness:.3f}")
            st.metric("Kurtosis", f"{profile.expression_kurtosis:.3f}")
            st.metric("IQR", f"{profile.expression_iqr:.2f}")
        
        with col3:
            st.metric("Expression Range", f"{profile.expression_range:.2f}")
            st.metric("25th Percentile", f"{profile.expression_percentiles.get(25, 0):.2f}")
            st.metric("75th Percentile", f"{profile.expression_percentiles.get(75, 0):.2f}")
    
    # Dispersion Details
    with st.expander("🎯 Detailed Dispersion Estimates"):
        st.markdown("""
        **Dispersion** measures biological variability after accounting for technical noise.
        Critical for choosing the right differential expression method.
        """)
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Mean Dispersion", f"{profile.dispersion_mean:.3f}")
            st.metric("Median Dispersion", f"{profile.dispersion_median:.3f}")
        
        with col2:
            st.metric("Common Dispersion", f"{profile.common_dispersion:.3f}")
            st.metric("BCV (sqrt of common)", f"{profile.bcv:.3f}")
        
        with col3:
            st.metric("Trended Slope", f"{profile.trended_dispersion_slope:.4f}")
            st.metric("Trended Intercept", f"{profile.trended_dispersion_intercept:.3f}")
    
    # Feature Vector for ML
    with st.expander("🤖 ML Feature Vector (32 features)"):
        st.markdown("**Complete feature vector for machine learning pipeline recommendation:**")
        
        features_df = pd.DataFrame([features]).T
        features_df.columns = ['Value']
        features_df.index.name = 'Feature'
        
        st.dataframe(features_df, use_container_width=True)
    
    # Visualizations
    st.markdown("### 📊 Interactive Visualizations")
    
    tab1, tab2, tab3, tab4 = st.tabs([
        "📚 Library Sizes", 
        "🧬 Gene Detection", 
        "📉 Mean-Variance",
        "🎨 Dispersion"
    ])
    
    with tab1:
        # Library size distribution
        count_data = st.session_state['m3_count_data']
        library_sizes = count_data.sum(axis=0)
        
        fig = go.Figure()
        
        fig.add_trace(go.Bar(
            x=library_sizes.index,
            y=library_sizes.values,
            marker_color='steelblue',
            name='Library Size'
        ))
        
        # Add mean line
        fig.add_hline(
            y=profile.library_size_mean,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Mean: {profile.library_size_mean:,.0f}",
            annotation_position="top right"
        )
        
        fig.update_layout(
            title="Library Sizes by Sample",
            xaxis_title="Sample",
            yaxis_title="Total Reads",
            height=500,
            showlegend=False
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        # Gene detection across samples
        detection_per_sample = (count_data > min_count_threshold).sum(axis=0)
        
        fig = go.Figure()
        
        fig.add_trace(go.Bar(
            x=detection_per_sample.index,
            y=detection_per_sample.values,
            marker_color='green',
            name='Detected Genes'
        ))
        
        fig.add_hline(
            y=profile.genes_per_sample_mean,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Mean: {profile.genes_per_sample_mean:,.0f}",
            annotation_position="top right"
        )
        
        fig.update_layout(
            title=f"Genes Detected per Sample (count > {min_count_threshold})",
            xaxis_title="Sample",
            yaxis_title="Number of Genes",
            height=500,
            showlegend=False
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        # Mean-variance relationship
        gene_means = count_data.mean(axis=1)
        gene_vars = count_data.var(axis=1)
        
        # Sample for plotting (too many points slow)
        if len(gene_means) > 5000:
            sample_idx = np.random.choice(len(gene_means), 5000, replace=False)
            plot_means = gene_means.iloc[sample_idx]
            plot_vars = gene_vars.iloc[sample_idx]
        else:
            plot_means = gene_means
            plot_vars = gene_vars
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=plot_means,
            y=plot_vars,
            mode='markers',
            marker=dict(
                size=3,
                color=plot_means,
                colorscale='Viridis',
                showscale=True,
                colorbar=dict(title="Mean")
            ),
            name='Genes'
        ))
        
        # Add y=x line (Poisson expectation)
        max_val = max(plot_means.max(), plot_vars.max())
        fig.add_trace(go.Scatter(
            x=[0, max_val],
            y=[0, max_val],
            mode='lines',
            line=dict(color='red', dash='dash'),
            name='Var = Mean (Poisson)'
        ))
        
        fig.update_layout(
            title="Mean-Variance Relationship",
            xaxis_title="Mean Expression",
            yaxis_title="Variance",
            xaxis_type="log",
            yaxis_type="log",
            height=600
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        st.caption("""
        **Interpretation:** 
        - Points on red line = Poisson distribution (var = mean)
        - Points above line = overdispersed (biological variability)
        - Points below line = underdispersed (technical artifacts)
        """)
    
    with tab4:
        # Dispersion estimates
        # Get per-gene dispersions from profiler if available
        st.markdown("**Biological Coefficient of Variation (BCV) Distribution**")
        
        # Create synthetic dispersion distribution for visualization
        # (In real implementation, would get from profiler.dispersions if available)
        fig = go.Figure()
        
        fig.add_vline(
            x=profile.bcv,
            line_dash="dash",
            line_color="red",
            annotation_text=f"BCV: {profile.bcv:.3f}",
            annotation_position="top"
        )
        
        # Add category zones
        fig.add_vrect(x0=0, x1=0.2, fillcolor="green", opacity=0.1, annotation_text="Low")
        fig.add_vrect(x0=0.2, x1=0.4, fillcolor="yellow", opacity=0.1, annotation_text="Moderate")
        fig.add_vrect(x0=0.4, x1=1.0, fillcolor="red", opacity=0.1, annotation_text="High")
        
        fig.update_layout(
            title="BCV Category Zones",
            xaxis_title="BCV",
            height=400
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        st.markdown(f"""
        **Your Data:**
        - BCV: **{profile.bcv:.3f}** ({profile.bcv_category})
        - Common Dispersion: **{profile.common_dispersion:.3f}**
        - Overdispersion Ratio: **{profile.overdispersion_ratio:.2f}** ({profile.overdispersion_category})
        """)
    
    # Export options
    st.markdown("### 💾 Export Profile")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        # Export full profile as JSON
        profile_json = profile.to_json(indent=2)
        st.download_button(
            "📥 Download Full Profile (JSON)",
            profile_json,
            "data_profile.json",
            "application/json",
            help="Complete profile with all 32 features"
        )
    
    with col2:
        # Export feature vector as CSV
        features_df = pd.DataFrame([features])
        features_csv = features_df.to_csv(index=False)
        st.download_button(
            "📥 Download Features (CSV)",
            features_csv,
            "profile_features.csv",
            "text/csv",
            help="32 features for ML pipeline recommendation"
        )
    
    with col3:
        # Export summary as text
        summary_text = profile.summary()
        st.download_button(
            "📥 Download Summary (TXT)",
            summary_text,
            "profile_summary.txt",
            "text/plain",
            help="Human-readable summary"
        )
    
    # Next steps
    st.markdown("---")
    st.info("""
    ✅ **Data profiling complete!**
    
    **Key Findings:**
    - BCV: {bcv:.3f} ({bcv_cat})
    - Sample design: {n_groups} groups, {min_size}-{max_size} samples per group
    - Expressed genes: {expressed:,} / {total:,} ({rate:.1f}%)
    - Sparsity: {sparsity:.1f}%
    
    **Next steps:**
    1. Review the 32 extracted features above
    2. Note BCV and sample characteristics
    3. Proceed to **Module 4: ML Recommender** for pipeline suggestion
    4. The feature vector will be used for ML-based recommendation
    """.format(
        bcv=profile.bcv,
        bcv_cat=profile.bcv_category,
        n_groups=profile.n_groups,
        min_size=profile.min_group_size,
        max_size=profile.max_group_size,
        expressed=profile.n_expressed_genes,
        total=profile.n_genes,
        rate=profile.detection_rate*100,
        sparsity=profile.sparsity*100
    ))

else:
    if not st.session_state.get('m3_complete', False):
        st.info("👆 Upload count data and run profiling to see results")
