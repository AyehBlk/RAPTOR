"""
RAPTOR Dashboard - Parameter Optimization (Module 8)

Interactive threshold optimization for differential expression results.
Integrates with DEResult objects from Module 7.

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

# Import RAPTOR modules
try:
    from raptor.de_import import DEResult
    DE_IMPORT_AVAILABLE = True
except ImportError as e:
    DE_IMPORT_AVAILABLE = False
    de_import_error = str(e)

try:
    from raptor.parameter_optimization import (
        ParameterOptimizer,
        ParameterSpace,
        OptimizationResult
    )
    OPTIMIZER_AVAILABLE = True
except ImportError as e:
    OPTIMIZER_AVAILABLE = False
    optimizer_error = str(e)

# Page config
st.set_page_config(
    page_title="RAPTOR - Parameter Optimization",
    page_icon="⚙️",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("⚙️ Parameter Optimization - Module 8")
st.markdown("### Interactive threshold tuning for differential expression")

# Help section
with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Optimize FDR and Log2FC thresholds to identify differentially expressed genes.
    
    ## 🎯 **Key Parameters**
    
    ### **Adjusted P-value (FDR)** - Statistical Significance
    - **Standard:** 0.05 (5% false discovery rate)
    - **Stringent:** 0.01 (1% FDR) - fewer false positives
    - **Exploratory:** 0.10 (10% FDR) - more discoveries
    
    ### **Log2 Fold Change** - Biological Significance
    - **Lenient:** ±0.58 (1.5-fold change)
    - **Standard:** ±1.0 (2-fold change)
    - **Stringent:** ±1.5 (3-fold change)
    - **Very Stringent:** ±2.0 (4-fold change)
    
    ## 📊 **Optimization Approaches**
    
    ### **1. Interactive Tuning** (This Module)
    - Use sliders to explore thresholds
    - See real-time impact on DEG counts
    - Visualize with volcano plot
    - Explore threshold combinations
    
    ### **2. Visual Guidance**
    - **Volcano Plot:** See DEG distribution
    - **Histograms:** Check p-value and fold change distributions
    - **Heatmap:** Explore threshold combinations
    - **Counts:** Track up/down regulated genes
    
    ### **3. Decision Criteria**
    
    **More Stringent (Fewer, Higher-Confidence Genes):**
    - Lower FDR (e.g., 0.01)
    - Higher Log2FC (e.g., 1.5 or 2.0)
    - Use when: Validation resources limited, need high confidence
    
    **Less Stringent (More Genes, Exploratory):**
    - Higher FDR (e.g., 0.10)
    - Lower Log2FC (e.g., 0.5 or 0.58)
    - Use when: Exploratory phase, hypothesis generation
    
    **Balanced (Standard):**
    - FDR: 0.05
    - Log2FC: 1.0
    - Use when: Typical differential expression study
    
    ## 🔄 **Workflow**
    ```
    Module 7: Import DE Results
         ↓
    Module 8: Optimize Thresholds ← YOU ARE HERE
         ↓
    Module 9: Ensemble Analysis
         ↓
    Robust Gene List
    ```
    
    ## 💡 **Tips**
    - Check p-value histogram for uniform distribution under null
    - Volcano plot should show symmetric wings if no batch effects
    - More stringent = fewer false positives, may miss true positives
    - Less stringent = more discoveries, more false positives
    - Use biological knowledge to set fold change thresholds
    - Standard thresholds (FDR<0.05, |Log2FC|>1.0) are reasonable defaults
    """)

# Check module availability
if not DE_IMPORT_AVAILABLE:
    st.error(f"""
    ❌ **RAPTOR de_import module not available**
    
    Error: {de_import_error}
    
    Please ensure RAPTOR is properly installed:
    ```bash
    cd RAPTOR/
    pip install -e .
    ```
    """)
    st.stop()

st.markdown("---")

# Check if data is available from Module 7
if not st.session_state.get('m7_complete', False):
    st.warning("⚠️ **No DE results found!**")
    st.info("Please import results in Module 7 first, then return here for optimization.")
    
    if st.button("📥 Go to Module 7 (Import DE)", type="primary"):
        try:
            st.switch_page("pages/04_📥_Import_DE.py")
        except:
            st.info("Navigate manually to Module 7 (Import DE)")
    
    st.stop()

# Get results from Module 7 (correct session state key)
de_results = st.session_state.get('m7_de_results', {})

if not de_results:
    st.error("❌ No DE results found in session state!")
    st.info("Please re-import your results in Module 7")
    st.stop()

# Section 1: Data Selection
st.markdown("## 1️⃣ Select Dataset for Optimization")

st.info(f"✅ Loaded **{len(de_results)}** DE result file(s) from Module 7")

# Show files info
files_info = []
for file_name, de_result in de_results.items():
    files_info.append({
        'File': file_name,
        'Pipeline': de_result.pipeline,
        'Genes': f"{de_result.n_genes:,}",
        'Current Sig.': f"{de_result.n_significant:,}"
    })

files_df = pd.DataFrame(files_info)
st.dataframe(files_df, use_container_width=True, hide_index=True)

# Select which result to optimize
if len(de_results) == 1:
    selected_file = list(de_results.keys())[0]
    st.success(f"📊 Using: **{selected_file}**")
else:
    selected_file = st.selectbox(
        "Select dataset to optimize:",
        list(de_results.keys()),
        help="Choose one result file for threshold optimization"
    )

# Get the DEResult object
de_result = de_results[selected_file]

# Access the standardized dataframe
df = de_result.results_df.copy()
pipeline = de_result.pipeline

# Get current parameters
current_fdr = de_result.parameters.get('fdr_threshold', 0.05)
current_lfc = de_result.parameters.get('lfc_threshold', 0.0)

st.success(f"""
✅ **Dataset:** {selected_file}  
📊 **Pipeline:** {pipeline}  
🧬 **Genes:** {de_result.n_genes:,}  
⚙️ **Current Thresholds:** FDR<{current_fdr}, |Log2FC|>{current_lfc}  
🎯 **Currently Significant:** {de_result.n_significant:,} genes
""")

st.markdown("---")

# Section 2: Current Statistics
st.markdown("## 2️⃣ Dataset Statistics")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Total Genes", f"{de_result.n_genes:,}")

with col2:
    st.metric("Current Sig.", f"{de_result.n_significant:,}")

with col3:
    st.metric("Up-regulated", f"{de_result.n_up:,}")

with col4:
    st.metric("Down-regulated", f"{de_result.n_down:,}")

# Distribution stats
with st.expander("📊 Distribution Statistics"):
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Adjusted P-value:**")
        st.write(f"• Min: {df['adjusted_p_value'].min():.2e}")
        st.write(f"• Median: {df['adjusted_p_value'].median():.2e}")
        st.write(f"• Max: {df['adjusted_p_value'].max():.2e}")
        st.write(f"• % < 0.05: {(df['adjusted_p_value'] < 0.05).sum() / len(df) * 100:.1f}%")
    
    with col2:
        st.markdown("**Log2 Fold Change:**")
        st.write(f"• Min: {df['log2_fold_change'].min():.2f}")
        st.write(f"• Median: {df['log2_fold_change'].median():.2f}")
        st.write(f"• Max: {df['log2_fold_change'].max():.2f}")
        st.write(f"• Mean |Log2FC|: {df['log2_fold_change'].abs().mean():.2f}")

st.markdown("---")

# Section 3: Interactive Threshold Tuning
st.markdown("## 3️⃣ Interactive Threshold Tuning")

st.info("""
💡 **How to use:**
1. Adjust sliders below to set FDR and Log2FC thresholds
2. Watch volcano plot and counts update in real-time
3. Explore different combinations
4. When satisfied, save parameters at the bottom
""")

col1, col2 = st.columns(2)

with col1:
    fdr_threshold = st.slider(
        "🎯 Adjusted P-value Threshold (FDR)",
        min_value=0.001,
        max_value=0.20,
        value=0.05,
        step=0.001,
        format="%.3f",
        help="Lower = more stringent (fewer genes)"
    )
    
    st.caption(f"📊 FDR: {fdr_threshold*100:.1f}% false discovery rate")
    
    # Category
    if fdr_threshold <= 0.01:
        st.caption("✅ Very Stringent")
    elif fdr_threshold <= 0.05:
        st.caption("✅ Standard")
    elif fdr_threshold <= 0.10:
        st.caption("⚠️ Lenient")
    else:
        st.caption("⚠️ Very Lenient")

with col2:
    lfc_threshold = st.slider(
        "🎯 Log2 Fold Change Threshold",
        min_value=0.0,
        max_value=3.0,
        value=1.0,
        step=0.1,
        format="%.1f",
        help="Higher = larger changes required"
    )
    
    fc = 2 ** lfc_threshold
    st.caption(f"📊 Fold change: {fc:.2f}×")
    
    # Category
    if lfc_threshold < 0.58:
        st.caption("⚠️ Very Lenient (< 1.5×)")
    elif lfc_threshold < 1.0:
        st.caption("⚠️ Lenient (< 2×)")
    elif lfc_threshold == 1.0:
        st.caption("✅ Standard (2×)")
    elif lfc_threshold <= 1.5:
        st.caption("✅ Stringent (≤ 3×)")
    else:
        st.caption("✅ Very Stringent (> 3×)")

# Calculate DEGs with current thresholds
mask_sig = df['adjusted_p_value'] < fdr_threshold
mask_up = (df['log2_fold_change'] > lfc_threshold) & mask_sig
mask_down = (df['log2_fold_change'] < -lfc_threshold) & mask_sig

n_up = mask_up.sum()
n_down = mask_down.sum()
n_total_deg = n_up + n_down

# Display results with current thresholds
st.markdown("### 🎯 Results with Current Thresholds")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        "⬆️ Upregulated",
        f"{n_up:,}",
        help=f"log2FC > {lfc_threshold} & FDR < {fdr_threshold}"
    )

with col2:
    st.metric(
        "⬇️ Downregulated",
        f"{n_down:,}",
        help=f"log2FC < -{lfc_threshold} & FDR < {fdr_threshold}"
    )

with col3:
    st.metric(
        "📊 Total DEGs",
        f"{n_total_deg:,}",
        help="Total differentially expressed genes"
    )

with col4:
    pct_deg = n_total_deg / de_result.n_genes * 100 if de_result.n_genes > 0 else 0
    st.metric(
        "% of Total",
        f"{pct_deg:.1f}%",
        help="Percentage of all genes"
    )

st.markdown("---")

# Section 4: Volcano Plot
st.markdown("## 4️⃣ Interactive Volcano Plot")

# Prepare data for plotting
plot_df = df.copy()
plot_df['-log10(FDR)'] = -np.log10(plot_df['adjusted_p_value'] + 1e-300)  # Avoid log(0)
plot_df['Category'] = 'Not Significant'
plot_df.loc[mask_up, 'Category'] = 'Upregulated'
plot_df.loc[mask_down, 'Category'] = 'Downregulated'

# Create volcano plot
fig = go.Figure()

# Not significant (gray)
not_sig = plot_df[plot_df['Category'] == 'Not Significant']
if len(not_sig) > 0:
    fig.add_trace(go.Scatter(
        x=not_sig['log2_fold_change'],
        y=not_sig['-log10(FDR)'],
        mode='markers',
        name='Not Significant',
        marker=dict(color='lightgray', size=4, opacity=0.4),
        hovertemplate='log2FC: %{x:.2f}<br>-log10(FDR): %{y:.2f}<extra></extra>'
    ))

# Upregulated (red)
up_df = plot_df[plot_df['Category'] == 'Upregulated']
if len(up_df) > 0:
    fig.add_trace(go.Scatter(
        x=up_df['log2_fold_change'],
        y=up_df['-log10(FDR)'],
        mode='markers',
        name=f'Upregulated ({n_up:,})',
        marker=dict(color='#d32f2f', size=6, opacity=0.7),
        hovertemplate='<b>Gene</b><br>log2FC: %{x:.2f}<br>-log10(FDR): %{y:.2f}<extra></extra>'
    ))

# Downregulated (blue)
down_df = plot_df[plot_df['Category'] == 'Downregulated']
if len(down_df) > 0:
    fig.add_trace(go.Scatter(
        x=down_df['log2_fold_change'],
        y=down_df['-log10(FDR)'],
        mode='markers',
        name=f'Downregulated ({n_down:,})',
        marker=dict(color='#1976d2', size=6, opacity=0.7),
        hovertemplate='<b>Gene</b><br>log2FC: %{x:.2f}<br>-log10(FDR): %{y:.2f}<extra></extra>'
    ))

# Add threshold lines
# Vertical lines for log2FC
fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="red", line_width=2,
              annotation_text=f"log2FC = {lfc_threshold}", annotation_position="top")
fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="blue", line_width=2,
              annotation_text=f"log2FC = -{lfc_threshold}", annotation_position="top")

# Horizontal line for FDR
fdr_line = -np.log10(fdr_threshold)
fig.add_hline(y=fdr_line, line_dash="dash", line_color="gray", line_width=2,
              annotation_text=f"FDR = {fdr_threshold}", annotation_position="right")

fig.update_layout(
    title=f"Volcano Plot - {selected_file}",
    xaxis_title="Log2 Fold Change",
    yaxis_title="-Log10(Adjusted P-value)",
    template="plotly_white",
    height=600,
    hovermode='closest',
    showlegend=True,
    legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01
    )
)

st.plotly_chart(fig, use_container_width=True)

st.caption("""
💡 **How to read:**
- **X-axis:** Magnitude of change (log2 scale)
- **Y-axis:** Statistical significance (-log10 scale, higher = more significant)
- **Red dots:** Significantly upregulated genes
- **Blue dots:** Significantly downregulated genes
- **Gray dots:** Not significant
- **Dashed lines:** Your current thresholds
""")

st.markdown("---")

# Section 5: Distributions
st.markdown("## 5️⃣ Parameter Distributions")

col1, col2 = st.columns(2)

with col1:
    # P-value distribution
    fig1 = go.Figure()
    
    fig1.add_trace(go.Histogram(
        x=df['adjusted_p_value'],
        nbinsx=50,
        name='All Genes',
        marker_color='#757575',
        opacity=0.7
    ))
    
    fig1.add_vline(x=fdr_threshold, line_dash="dash", line_color="red", line_width=2,
                   annotation_text=f"Threshold: {fdr_threshold}")
    
    fig1.update_layout(
        title="Adjusted P-value Distribution",
        xaxis_title="Adjusted P-value",
        yaxis_title="Count",
        template="plotly_white",
        height=400,
        showlegend=False
    )
    
    st.plotly_chart(fig1, use_container_width=True)
    
    st.caption(f"📊 {(df['adjusted_p_value'] < fdr_threshold).sum():,} genes below threshold")

with col2:
    # Log2FC distribution
    fig2 = go.Figure()
    
    fig2.add_trace(go.Histogram(
        x=df['log2_fold_change'],
        nbinsx=50,
        name='All Genes',
        marker_color='#757575',
        opacity=0.7
    ))
    
    fig2.add_vline(x=lfc_threshold, line_dash="dash", line_color="red", line_width=2,
                   annotation_text=f"+{lfc_threshold}")
    fig2.add_vline(x=-lfc_threshold, line_dash="dash", line_color="blue", line_width=2,
                   annotation_text=f"-{lfc_threshold}")
    
    fig2.update_layout(
        title="Log2 Fold Change Distribution",
        xaxis_title="Log2 Fold Change",
        yaxis_title="Count",
        template="plotly_white",
        height=400,
        showlegend=False
    )
    
    st.plotly_chart(fig2, use_container_width=True)
    
    st.caption(f"📊 {(df['log2_fold_change'].abs() > lfc_threshold).sum():,} genes above threshold")

st.markdown("---")

# Section 6: Threshold Effects Analysis
st.markdown("## 6️⃣ Threshold Combinations Explorer")

st.info("""
💡 **Explore different threshold combinations:**
This heatmap shows how many DEGs you'd get with different FDR and Log2FC thresholds.
Use it to understand the trade-off between stringency and discovery.
""")

# Calculate DEG counts for various thresholds
fdr_range = [0.001, 0.01, 0.05, 0.10, 0.20]
lfc_range = [0.0, 0.5, 1.0, 1.5, 2.0]

threshold_results = []

for fdr in fdr_range:
    for lfc in lfc_range:
        mask = (df['adjusted_p_value'] < fdr) & (df['log2_fold_change'].abs() > lfc)
        n_genes = mask.sum()
        threshold_results.append({
            'FDR': fdr,
            'Log2FC': lfc,
            'n_genes': n_genes
        })

threshold_df = pd.DataFrame(threshold_results)

# Create heatmap
pivot_df = threshold_df.pivot(index='Log2FC', columns='FDR', values='n_genes')

fig = go.Figure(data=go.Heatmap(
    z=pivot_df.values,
    x=[f"{x:.3f}" for x in pivot_df.columns],
    y=[f"{y:.1f}" for y in pivot_df.index],
    colorscale='Viridis',
    text=pivot_df.values,
    texttemplate='%{text}',
    textfont={"size": 12},
    colorbar=dict(title="# DEGs"),
    hovertemplate='FDR: %{x}<br>Log2FC: %{y}<br>DEGs: %{z}<extra></extra>'
))

# Mark current selection
# Find closest indices
fdr_idx = min(range(len(fdr_range)), key=lambda i: abs(fdr_range[i] - fdr_threshold))
lfc_idx = min(range(len(lfc_range)), key=lambda i: abs(lfc_range[i] - lfc_threshold))

fig.add_trace(go.Scatter(
    x=[f"{fdr_range[fdr_idx]:.3f}"],
    y=[f"{lfc_range[lfc_idx]:.1f}"],
    mode='markers',
    marker=dict(size=20, color='red', symbol='x', line=dict(width=2)),
    name='Current Selection',
    showlegend=True
))

fig.update_layout(
    title="DEG Count by Threshold Combinations",
    xaxis_title="FDR Threshold",
    yaxis_title="Log2FC Threshold",
    template="plotly_white",
    height=500
)

st.plotly_chart(fig, use_container_width=True)

st.caption("""
💡 **Reading the heatmap:**
- **Upper-left:** More lenient (many genes)
- **Lower-right:** More stringent (few genes)
- **Red X:** Your current selection
- **Numbers:** How many DEGs at each combination
""")

st.markdown("---")

# Section 7: Save & Export
st.markdown("## 7️⃣ Save Optimized Parameters & Export")

st.info("""
📌 **When you're happy with the thresholds:**
1. Click "Save Optimized Parameters" to store for Module 9
2. Download filtered gene lists as needed
3. Continue to Ensemble Analysis
""")

col1, col2 = st.columns([1, 3])

with col1:
    if st.button("💾 Save Optimized Parameters", type="primary", use_container_width=True):
        # Create new DEResult with optimized thresholds
        optimized_df = df.copy()
        
        # Recalculate significance with new thresholds
        optimized_df['is_significant'] = (
            (optimized_df['adjusted_p_value'] < fdr_threshold) &
            (optimized_df['log2_fold_change'].abs() > lfc_threshold)
        )
        
        # Update direction
        optimized_df['direction'] = 'unchanged'
        optimized_df.loc[
            (optimized_df['log2_fold_change'] > lfc_threshold) & optimized_df['is_significant'],
            'direction'
        ] = 'up'
        optimized_df.loc[
            (optimized_df['log2_fold_change'] < -lfc_threshold) & optimized_df['is_significant'],
            'direction'
        ] = 'down'
        
        # Create new DEResult
        from raptor.de_import import DEResult
        
        optimized_result = DEResult(
            results_df=optimized_df,
            pipeline=de_result.pipeline,
            parameters={
                'fdr_threshold': fdr_threshold,
                'lfc_threshold': lfc_threshold,
                'source_pipeline': de_result.pipeline,
                'optimized': True
            },
            metadata={
                **de_result.metadata,
                'optimized_in_module': 'M8',
                'source_file': selected_file
            }
        )
        
        # Store in session state
        st.session_state['m8_optimized_result'] = optimized_result
        st.session_state['m8_optimized_params'] = {
            'fdr_threshold': fdr_threshold,
            'lfc_threshold': lfc_threshold,
            'n_up': int(n_up),
            'n_down': int(n_down),
            'n_total': int(n_total_deg),
            'source_file': selected_file
        }
        st.session_state['m8_complete'] = True
        st.session_state['m8_running'] = False
        st.session_state['m8_error'] = False
        
        st.success("✅ Parameters saved successfully!")
        st.balloons()

with col2:
    if st.session_state.get('m8_complete', False):
        params = st.session_state['m8_optimized_params']
        st.success(f"""
        ✅ **Saved!**  
        DEGs: **{params['n_total']:,}** ({params['n_up']:,} up, {params['n_down']:,} down)  
        FDR < {params['fdr_threshold']}, |Log2FC| > {params['lfc_threshold']}
        """)

# Export options
st.markdown("### 📥 Export Gene Lists")

col1, col2, col3 = st.columns(3)

with col1:
    # Export all DEGs
    if n_total_deg > 0:
        export_df = df[mask_sig & (df['log2_fold_change'].abs() > lfc_threshold)].copy()
        csv_data = export_df.to_csv(index=True)
        
        st.download_button(
            "📥 All DEGs (CSV)",
            csv_data,
            f"degs_fdr{fdr_threshold}_lfc{lfc_threshold}.csv",
            "text/csv",
            use_container_width=True,
            help=f"{n_total_deg:,} genes"
        )
    else:
        st.button("📥 All DEGs (CSV)", disabled=True, use_container_width=True)
        st.caption("No DEGs with current thresholds")

with col2:
    # Export upregulated only
    if n_up > 0:
        up_df = df[mask_up].copy()
        csv_data = up_df.to_csv(index=True)
        
        st.download_button(
            "📥 Upregulated Only",
            csv_data,
            f"upregulated_fdr{fdr_threshold}_lfc{lfc_threshold}.csv",
            "text/csv",
            use_container_width=True,
            help=f"{n_up:,} genes"
        )
    else:
        st.button("📥 Upregulated Only", disabled=True, use_container_width=True)
        st.caption("No upregulated genes")

with col3:
    # Export downregulated only
    if n_down > 0:
        down_df = df[mask_down].copy()
        csv_data = down_df.to_csv(index=True)
        
        st.download_button(
            "📥 Downregulated Only",
            csv_data,
            f"downregulated_fdr{fdr_threshold}_lfc{lfc_threshold}.csv",
            "text/csv",
            use_container_width=True,
            help=f"{n_down:,} genes"
        )
    else:
        st.button("📥 Downregulated Only", disabled=True, use_container_width=True)
        st.caption("No downregulated genes")

# Continue to next module
if st.session_state.get('m8_complete', False):
    st.markdown("---")
    st.markdown("### ➡️ Next Step")
    
    if st.button("Continue to Ensemble Analysis (Module 9)", type="primary", use_container_width=True):
        try:
            st.switch_page("pages/05_🔬_Ensemble.py")
        except:
            st.info("Navigate manually to Module 9 (Ensemble Analysis)")

# Footer
st.markdown("---")
st.info(f"""
💡 **Current Selection:** FDR < {fdr_threshold}, |Log2FC| > {lfc_threshold} → **{n_total_deg:,} DEGs**

**Guidelines:**
- Standard thresholds: FDR < 0.05, |Log2FC| > 1.0
- More stringent = fewer false positives, may miss true positives
- Less stringent = more discoveries, more false positives
- Adjust based on your research goals and validation resources
""")

st.caption("📚 **Next:** Module 9 - Ensemble Analysis (if you have results from multiple methods)")
