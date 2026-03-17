"""
RAPTOR Dashboard - Ensemble Analysis (Module 9)

Combine multiple DE methods for robust consensus gene identification.
Uses Fisher's, Brown's, and RRA methods (functional implementations only).

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
from typing import Dict, List

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
    from raptor.ensemble import (
        EnsembleResult,
        ensemble_fisher,
        ensemble_brown,
        ensemble_rra
    )
    ENSEMBLE_AVAILABLE = True
except ImportError as e:
    ENSEMBLE_AVAILABLE = False
    ensemble_error = str(e)

# Page config
st.set_page_config(
    page_title="RAPTOR - Ensemble Analysis",
    page_icon="🔬",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("Ensemble Analysis — Module 9")
st.caption("Combine multiple DE methods for robust consensus gene identification")

# Help section
with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Combine DE results from multiple tools to identify robust consensus genes.
    
    ## **Available Ensemble Methods**
    
    ### **1. Fisher's Method**
    - Combines p-values using Fisher's combined probability test
    - Assumes methods are independent
    - Use when: Methods use different algorithms
    
    ### **2. Brown's Method** (Recommended)
    - Like Fisher's but accounts for correlation
    - Best for typical RNA-seq (DESeq2 + edgeR)
    - Use when: Methods use same data
    
    ### **3. RRA (Robust Rank Aggregation)**
    - Ranks genes in each method, aggregates robustly
    - Resistant to outliers
    - Use when: Methods produce different p-value distributions
    
    ## ℹ️ **Note**
    Additional methods (Voting, Weighted) are being upgraded and will be available soon.
    """)

# Check module availability
if not DE_IMPORT_AVAILABLE:
    st.error(f"RAPTOR de_import module not available: {de_import_error}")
    st.stop()

if not ENSEMBLE_AVAILABLE:
    st.error(f"RAPTOR ensemble module not available: {ensemble_error}")
    st.stop()

st.markdown("---")

# Check prerequisites
if not st.session_state.get('m7_complete', False):
    st.warning("**No DE results found!**")
    st.info("Please import results in Module 7 first.")
    
    if st.button("Go to Module 7 (Import DE)", type="primary"):
        try:
            st.switch_page("pages/04_📥_Import_DE.py")
        except:
            st.info("Navigate manually to Module 7")
    st.stop()

# Get results - CORRECT KEY!
de_results = st.session_state.get('m7_de_results', {})

if not de_results:
    st.error("No DE results in session state!")
    st.stop()

# Section 1: Data Selection
st.markdown("## 1. Select DE Results")

st.info(f"Loaded **{len(de_results)}** result(s) from Module 7")

# Show summary
result_summary = []
for file_name, de_result in de_results.items():
    result_summary.append({
        'File': file_name,
        'Pipeline': de_result.pipeline,
        'Total': f"{de_result.n_genes:,}",
        'Significant': f"{de_result.n_significant:,}",
        'Up': f"{de_result.n_up:,}",
        'Down': f"{de_result.n_down:,}"
    })

summary_df = pd.DataFrame(result_summary)
st.dataframe(summary_df, use_container_width=True, hide_index=True)

if len(de_results) < 2:
    st.warning(f"Only {len(de_results)} result. Ensemble works best with ≥2 methods.")

selected_files = st.multiselect(
    "Select results for ensemble:",
    list(de_results.keys()),
    default=list(de_results.keys())
)

if len(selected_files) == 0:
    st.stop()

selected_de_results = {name: de_results[name] for name in selected_files}

# Quick overlap preview
if len(selected_files) >= 2:
    with st.expander("Quick Overlap"):
        sig_genes = {n: set(de_results[n].significant_genes) for n in selected_files}
        files = list(sig_genes.keys())
        
        overlaps = []
        for i in range(len(files)):
            for j in range(i+1, len(files)):
                overlap = len(sig_genes[files[i]] & sig_genes[files[j]])
                union = len(sig_genes[files[i]] | sig_genes[files[j]])
                jaccard = overlap / union if union > 0 else 0
                overlaps.append({
                    'Method 1': files[i], 
                    'Method 2': files[j],
                    'Overlap': overlap, 
                    'Jaccard': f"{jaccard:.2f}"
                })
        
        if overlaps:
            st.dataframe(pd.DataFrame(overlaps), use_container_width=True, hide_index=True)
            avg_j = np.mean([float(r['Jaccard']) for r in overlaps])
            if avg_j >= 0.6:
                st.success(f"Good overlap (avg: {avg_j:.2f})")
            elif avg_j >= 0.4:
                st.info(f"ℹ️ Moderate overlap (avg: {avg_j:.2f})")
            else:
                st.warning(f"Low overlap (avg: {avg_j:.2f})")

st.markdown("---")

# Section 2: Method Selection
st.markdown("## 2. Select Ensemble Method")

st.info("""
ℹ️ **Available:** Fisher's, Brown's, RRA

**Coming Soon:** Voting and Weighted methods (being upgraded)
""")

ensemble_method = st.radio(
    "Choose ensemble method:",
    [
        "Fisher's Method",
        "Brown's Method (Recommended)",
        "RRA (Robust Rank Aggregation)"
    ],
    index=1,
    help="Brown's recommended for typical RNA-seq. Fisher's for independent methods. RRA for heterogeneous methods."
)

# Method explanation
if "Fisher" in ensemble_method:
    st.info("**Fisher's Method:** Combines p-values assuming independence. Use for truly independent methods.")
elif "Brown" in ensemble_method:
    st.success("**Brown's Method:** Accounts for correlation. Best for DESeq2+edgeR (recommended).")
elif "RRA" in ensemble_method:
    st.info("**RRA:** Rank-based aggregation. Robust to outliers. Good for heterogeneous methods.")

st.markdown("---")

# Section 3: Configuration
st.markdown("## 3. Configuration")

col1, col2 = st.columns(2)

with col1:
    require_same_direction = st.checkbox(
        "Require same direction",
        value=True,
        help="Gene must change in same direction across methods"
    )
    
    if require_same_direction:
        st.caption("Recommended - reduces artifacts")
    else:
        st.caption("Not recommended")

with col2:
    use_filters = st.checkbox(
        "Apply pre-filters (optional)",
        value=False,
        help="Filter genes before ensemble"
    )

if use_filters:
    st.markdown("**Filter Settings:**")
    col1, col2 = st.columns(2)
    
    with col1:
        filter_fdr = st.number_input(
            "Max FDR", 0.001, 1.0, 0.05, 0.01, format="%.3f"
        )
    
    with col2:
        filter_lfc = st.number_input(
            "Min |Log2FC|", 0.0, 5.0, 1.0, 0.5, format="%.1f"
        )
    
    filters = {'padj': filter_fdr, 'lfc': filter_lfc}
else:
    filters = None

st.markdown("---")

# Section 4: Run Analysis
st.markdown("## 4. Run Ensemble Analysis")

st.info(f"""
**Configuration:**
- Methods: {len(selected_files)}
- Ensemble: {ensemble_method.split('(')[0].strip()}
- Direction check: {'YES' if require_same_direction else 'NO'}
- Filters: {'FDR<' + str(filter_fdr) + ', |Log2FC|>' + str(filter_lfc) if use_filters else 'None'}
""")

if st.button("Run Ensemble Analysis", type="primary", use_container_width=True):
    
    if len(selected_files) < 2:
        st.warning("Only 1 method - ensemble works best with 2+")
    
    with st.spinner("Running ensemble analysis..."):
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        try:
            status_text.text("Preparing data...")
            progress_bar.progress(20)
            
            # Prepare input for ensemble functions
            from types import SimpleNamespace
            
            # Column name mapping: Module 7 standardized → ensemble expected
            col_rename = {
                'log2_fold_change': 'log2FoldChange',
                'p_value': 'pvalue',
                'adjusted_p_value': 'padj',
                'base_mean': 'baseMean',
            }
            
            ensemble_input = {}
            for name, de_result in selected_de_results.items():
                df = de_result.results_df.copy()
                # Ensure gene_id is a column, not just the index
                if 'gene_id' not in df.columns:
                    df = df.reset_index()
                # Rename standardized columns back to what ensemble expects
                df = df.rename(columns={k: v for k, v in col_rename.items() if k in df.columns})
                result_obj = SimpleNamespace(data=df)
                ensemble_input[name] = result_obj
            
            status_text.text(f"Running {ensemble_method.split('(')[0].strip()}...")
            progress_bar.progress(50)
            
            # Pre-filter data if filters are set
            if filters is not None:
                for name, obj in ensemble_input.items():
                    df = obj.data.copy()
                    if 'padj' in df.columns and 'padj' in filters:
                        df = df[df['padj'] <= filters['padj']]
                    if 'log2FoldChange' in df.columns and 'lfc' in filters:
                        df = df[df['log2FoldChange'].abs() >= filters['lfc']]
                    obj.data = df
            
            # Run appropriate method
            if "Fisher" in ensemble_method:
                ensemble_result = ensemble_fisher(
                    de_results=ensemble_input,
                    check_direction=require_same_direction
                )
                ensemble_type = "Fisher's Method"
            elif "Brown" in ensemble_method:
                ensemble_result = ensemble_brown(
                    de_results=ensemble_input,
                    check_direction=require_same_direction
                )
                ensemble_type = "Brown's Method"
            elif "RRA" in ensemble_method:
                ensemble_result = ensemble_rra(
                    de_results=ensemble_input,
                    check_direction=require_same_direction
                )
                ensemble_type = "RRA"
            
            status_text.text("Processing results...")
            progress_bar.progress(80)
            
            # Extract consensus genes
            consensus_df = ensemble_result.consensus_genes.copy()
            
            # Determine score column
            if 'combined_padj' in consensus_df.columns:
                score_col = 'combined_padj'
            elif 'rra_pvalue' in consensus_df.columns:
                score_col = 'rra_pvalue'
            else:
                score_col = None
            
            # Add tiers
            if score_col:
                consensus_df['tier'] = consensus_df[score_col].apply(
                    lambda x: 'Tier 1 (p<0.001)' if x < 0.001
                    else ('Tier 2 (p<0.01)' if x < 0.01 else 'Tier 3 (p<0.05)')
                )
                
                tier1_count = (consensus_df[score_col] < 0.001).sum()
                tier2_count = ((consensus_df[score_col] >= 0.001) & (consensus_df[score_col] < 0.01)).sum()
                tier3_count = ((consensus_df[score_col] >= 0.01) & (consensus_df[score_col] < 0.05)).sum()
            else:
                consensus_df['tier'] = 'Tier 2 (p<0.01)'
                tier1_count = 0
                tier2_count = len(consensus_df)
                tier3_count = 0
            
            status_text.text("Finalizing...")
            progress_bar.progress(100)
            
            # Store results
            result_package = {
                'ensemble_result': ensemble_result,
                'consensus_genes': consensus_df,
                'n_consensus': len(consensus_df),
                'tier1_count': int(tier1_count),
                'tier2_count': int(tier2_count),
                'tier3_count': int(tier3_count),
                'ensemble_method': ensemble_type,
                'n_methods': len(selected_files),
                'method_names': selected_files,
                'score_column': score_col,
                'parameters': {
                    'require_same_direction': require_same_direction,
                    'filters': filters
                }
            }
            
            st.session_state['m9_result'] = result_package
            st.session_state['m9_complete'] = True
            st.session_state['m9_running'] = False
            st.session_state['m9_error'] = False
            
            progress_bar.empty()
            status_text.empty()
            
            st.success(f"""
            **Complete!** Found **{len(consensus_df):,}** consensus genes
            
            - Tier 1: {tier1_count:,} genes (p<0.001)
            - Tier 2: {tier2_count:,} genes (p<0.01)
            - Tier 3: {tier3_count:,} genes (p<0.05)
            """)
            st.balloons()
            
        except Exception as e:
            progress_bar.empty()
            status_text.empty()
            
            st.error(f"""
            **Error:** {str(e)}
            
            **Troubleshooting:**
            - Check gene IDs are consistent
            - Try without filters
            - Verify same gene ID format
            """)
            st.session_state['m9_error'] = True

# Section 5: Results
if st.session_state.get('m9_complete', False):
    st.markdown("---")
    st.markdown("## 5. Ensemble Results")
    
    result = st.session_state['m9_result']
    consensus_df = result['consensus_genes']
    score_col = result['score_column']
    
    # Overall summary
    st.markdown("### Summary")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Consensus", f"{result['n_consensus']:,}")
    
    with col2:
        st.metric("Methods", result['n_methods'])
    
    with col3:
        st.metric("Method", result['ensemble_method'])
    
    with col4:
        if result['n_methods'] > 1:
            avg = sum(de_results[n].n_significant for n in result['method_names']) / result['n_methods']
            rate = (result['n_consensus'] / avg * 100) if avg > 0 else 0
            st.metric("Consensus Rate", f"{rate:.0f}%")
    
    # Tiers
    st.markdown("### Confidence Tiers")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("🥇 Tier 1", f"{result['tier1_count']:,}", help="p<0.001")
        st.caption("✅ qPCR validation")
    
    with col2:
        st.metric("🥈 Tier 2", f"{result['tier2_count']:,}", help="p<0.01")
        st.caption("✅ Pathway analysis")
    
    with col3:
        st.metric("🥉 Tier 3", f"{result['tier3_count']:,}", help="p<0.05")
        st.caption("ℹ️ Exploratory")
    
    # Visualizations
    st.markdown("---")
    st.markdown("### Visualizations")
    
    tab1, tab2, tab3 = st.tabs(["Tiers", "P-values", "🔝 Top Genes"])
    
    with tab1:
        # Pie chart
        tier_counts = consensus_df['tier'].value_counts()
        
        fig = go.Figure(data=[go.Pie(
            labels=tier_counts.index,
            values=tier_counts.values,
            marker=dict(colors=['#2E7D32', '#FFA726', '#5C6BC0']),
            hole=0.3,
            textinfo='label+value+percent'
        )])
        
        fig.update_layout(
            title=f"Consensus Genes by Tier ({len(consensus_df):,} total)",
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        # P-value distribution
        if score_col:
            plot_df = consensus_df.copy()
            plot_df['-log10(p)'] = -np.log10(plot_df[score_col] + 1e-300)
            
            fig = go.Figure()
            
            fig.add_trace(go.Histogram(
                x=plot_df['-log10(p)'],
                nbinsx=50,
                marker_color='#2E7D32',
                opacity=0.7
            ))
            
            # Threshold lines
            fig.add_vline(x=-np.log10(0.05), line_dash="dash", line_color="gray",
                         annotation_text="p=0.05")
            fig.add_vline(x=-np.log10(0.01), line_dash="dash", line_color="orange",
                         annotation_text="p=0.01")
            fig.add_vline(x=-np.log10(0.001), line_dash="dash", line_color="red",
                         annotation_text="p=0.001")
            
            fig.update_layout(
                title="Distribution of Combined P-values",
                xaxis_title="-Log10(P-value)",
                yaxis_title="Count",
                template="plotly_white",
                height=500,
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        # Top genes
        st.markdown("**Top 20 Genes:**")
        
        if score_col:
            top_genes = consensus_df.nsmallest(20, score_col)
        else:
            top_genes = consensus_df.head(20)
        
        # Select columns
        display_cols = ['gene_id']
        if score_col:
            display_cols.append(score_col)
        if 'direction' in top_genes.columns:
            display_cols.append('direction')
        if 'meta_lfc' in top_genes.columns:
            display_cols.append('meta_lfc')
        display_cols.append('tier')
        
        display_df = top_genes[display_cols].copy()
        
        # Format p-values
        if score_col:
            display_df[score_col] = display_df[score_col].apply(lambda x: f"{x:.2e}")
        
        st.dataframe(display_df, use_container_width=True, height=400)
    
    # Full table
    st.markdown("---")
    st.markdown("### All Consensus Genes")
    
    col1, col2 = st.columns([1, 3])
    
    with col1:
        tier_filter = st.selectbox(
            "Tier:",
            ["All"] + sorted(consensus_df['tier'].unique().tolist())
        )
    
    with col2:
        if 'direction' in consensus_df.columns:
            dir_filter = st.selectbox("Direction:", ["All", "up", "down"])
        else:
            dir_filter = "All"
    
    # Apply filters
    display_df = consensus_df.copy()
    
    if tier_filter != "All":
        display_df = display_df[display_df['tier'] == tier_filter]
    
    if dir_filter != "All" and 'direction' in display_df.columns:
        display_df = display_df[display_df['direction'] == dir_filter]
    
    st.dataframe(display_df, use_container_width=True, height=400)
    st.caption(f"Showing {len(display_df):,} of {len(consensus_df):,} genes")
    
    # Export
    st.markdown("---")
    st.markdown("## 6️⃣ Export")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        # Tier 1
        if result['tier1_count'] > 0:
            tier1_df = consensus_df[consensus_df['tier'] == 'Tier 1 (p<0.001)']
            st.download_button(
                "Tier 1",
                tier1_df.to_csv(index=False),
                "tier1_genes.csv",
                "text/csv",
                use_container_width=True
            )
        else:
            st.button("Tier 1", disabled=True, use_container_width=True)
    
    with col2:
        # All
        st.download_button(
            "All Genes",
            consensus_df.to_csv(index=False),
            "all_consensus.csv",
            "text/csv",
            use_container_width=True
        )
    
    with col3:
        # Upregulated
        if 'direction' in consensus_df.columns:
            up_df = consensus_df[consensus_df['direction'] == 'up']
            if len(up_df) > 0:
                st.download_button(
                    "Up",
                    up_df.to_csv(index=False),
                    "upregulated.csv",
                    "text/csv",
                    use_container_width=True
                )
            else:
                st.button("Up", disabled=True, use_container_width=True)
        else:
            st.button("Up", disabled=True, use_container_width=True)
    
    with col4:
        # Downregulated
        if 'direction' in consensus_df.columns:
            down_df = consensus_df[consensus_df['direction'] == 'down']
            if len(down_df) > 0:
                st.download_button(
                    "Down",
                    down_df.to_csv(index=False),
                    "downregulated.csv",
                    "text/csv",
                    use_container_width=True
                )
            else:
                st.button("Down", disabled=True, use_container_width=True)
        else:
            st.button("Down", disabled=True, use_container_width=True)
    
    # Summary
    st.markdown("---")
    st.info(f"""
    **Ensemble complete!**
    
    - **{result['n_consensus']:,}** consensus genes
    - **{result['tier1_count']:,}** highest confidence (Tier 1)
    - Method: {result['ensemble_method']}
    
    **Next:** Use Tier 1 for validation, Tier 1+2 for pathways
    """)

# Footer
st.markdown("---")
st.caption("Tier 1 genes (p<0.001) are most reliable for validation")
st.caption("Brown's Method recommended for typical RNA-seq (DESeq2+edgeR)")
