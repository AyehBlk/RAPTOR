"""
RAPTOR Dashboard - Import DE Results (Module 7)

Import and standardize differential expression results from
DESeq2, edgeR, limma-voom, and Wilcoxon.

Creates DEResult objects for downstream ensemble analysis.

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
from io import BytesIO
import json

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# Import RAPTOR de_import module
try:
    from raptor.de_import import (
        import_de_results,
        detect_pipeline,
        standardize_columns,
        calculate_significance,
        DEResult,
        COLUMN_MAPPINGS
    )
    DE_IMPORT_AVAILABLE = True
except ImportError as e:
    DE_IMPORT_AVAILABLE = False
    import_error = str(e)

# Page config
st.set_page_config(
    page_title="RAPTOR - Import DE",
    page_icon="📥",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# Main content
st.title("Import DE Results — Module 7")
st.caption("Import and standardize DE results from DESeq2, edgeR, limma, and Wilcoxon")

# Help section
with st.expander("ℹ️ How to use this module"):
    st.markdown("""
    **Purpose:** Import DE results from popular R packages and standardize them for ensemble analysis.
    
    **Supported Tools:**
    - **DESeq2** - Most widely used (~60% of publications)
    - **edgeR** - Excellent for low counts (~25% of publications)
    - **limma-voom** - Fast for large samples (~10% of publications)
    - **Wilcoxon** - Non-parametric test (~5% of publications)
    
    **What This Module Does:**
    1. **Auto-detects** which R package was used
    2. **Standardizes** column names across different formats
    3. **Calculates** significance based on your thresholds
    4. **Creates** DEResult objects for Module 9 (Ensemble Analysis)
    
    **Input File Format:**
    - CSV or TSV file
    - First column: Gene IDs (Ensembl IDs or gene symbols)
    - Required columns vary by tool (see below)
    
    **DESeq2 Expected Columns:**
    - `log2FoldChange` or `log2FC`
    - `pvalue` or `pval`
    - `padj` (adjusted p-value)
    - `baseMean` (optional)
    
    **edgeR Expected Columns:**
    - `logFC`
    - `PValue`
    - `FDR`
    - `logCPM` (optional)
    
    **limma-voom Expected Columns:**
    - `logFC`
    - `P.Value`
    - `adj.P.Val`
    - `AveExpr` (optional)
    
    **Wilcoxon Expected Columns:**
    - `logFC` or `log2FC`
    - `pvalue`
    - `padj` or `FDR`
    
    **Standardized Output:**
    All results are converted to:
    - `log2_fold_change`
    - `p_value`
    - `adjusted_p_value`
    - `base_mean` (if available)
    - `is_significant` (calculated)
    - `direction` ('up', 'down', 'unchanged')
    
    **Tips:**
    - Upload at least **2 files from different methods** for meaningful ensemble analysis
    - Gene IDs must be **consistent across files** (all Ensembl or all symbols)
    - Files should be from the **same comparison** (e.g., all treatment vs control)
    
    **Workflow:**
    ```
    R Analysis (External) → CSV files
         ↓
    Module 7 (This) → Import & Standardize
         ↓
    Module 9 → Ensemble Analysis
         ↓
    Robust Gene List
    ```
    """)

# Check module availability
if not DE_IMPORT_AVAILABLE:
    st.error(f"""
    **RAPTOR de_import module not available**
    
    Error: {import_error}
    
    Please ensure RAPTOR is properly installed:
    ```bash
    cd RAPTOR/
    pip install -e .
    ```
    """)
    st.stop()

st.markdown("---")

# Section 1: File Upload
st.markdown("## 1. Upload DE Result Files")

st.info("""
**What to upload:**
- Results from your R analysis (DESeq2, edgeR, limma-voom, or Wilcoxon)
- CSV or TSV format
- Gene IDs as first column or index
- For ensemble analysis: Upload results from **2+ different methods** on the **same comparison**
""")

uploaded_files = st.file_uploader(
    "Select DE result files (CSV or TSV)",
    type=['csv', 'tsv', 'txt'],
    accept_multiple_files=True,
    help="Upload one or more DE result files from R analysis"
)

if uploaded_files:
    st.success(f"Uploaded {len(uploaded_files)} file(s)")
    
    # Section 2: Format Detection & Validation
    st.markdown("---")
    st.markdown("## 2. Format Detection & Validation")
    
    detection_results = []
    file_data = {}
    
    with st.spinner("Detecting file formats..."):
        for uploaded_file in uploaded_files:
            try:
                # Read file
                if uploaded_file.name.endswith('.csv'):
                    df = pd.read_csv(uploaded_file)
                else:
                    df = pd.read_csv(uploaded_file, sep='\t')
                
                # Detect pipeline
                try:
                    pipeline = detect_pipeline(df)
                    detected_tool = pipeline.upper()
                    valid = True
                    error = None
                except Exception as e:
                    detected_tool = "Unknown"
                    valid = False
                    error = str(e)
                
                # Store results
                detection_results.append({
                    'name': uploaded_file.name,
                    'tool': detected_tool,
                    'valid': valid,
                    'error': error,
                    'n_genes': len(df),
                    'columns': list(df.columns)[:5]  # First 5 columns
                })
                
                file_data[uploaded_file.name] = df
                
            except Exception as e:
                detection_results.append({
                    'name': uploaded_file.name,
                    'tool': 'ERROR',
                    'valid': False,
                    'error': str(e),
                    'n_genes': 0,
                    'columns': []
                })
    
    # Display detection results
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Create detection summary table
        detection_df = pd.DataFrame([
            {
                'File': r['name'],
                'Detected Tool': r['tool'],
                'Genes': r['n_genes'],
                'Status': 'Valid' if r['valid'] else 'Invalid'
            }
            for r in detection_results
        ])
        
        st.dataframe(detection_df, use_container_width=True, hide_index=True)
    
    with col2:
        # Summary metrics
        valid_count = sum(r['valid'] for r in detection_results)
        st.metric("Valid Files", f"{valid_count}/{len(detection_results)}")
        
        if valid_count > 0:
            tools = [r['tool'] for r in detection_results if r['valid']]
            unique_tools = len(set(tools))
            st.metric("Unique Tools", unique_tools)
            
            if unique_tools >= 2:
                st.success("Good for ensemble!")
            else:
                st.warning("Consider adding another method")
    
    # Show column mapping reference
    with st.expander("Column Mapping Reference"):
        st.markdown("**How RAPTOR recognizes different formats:**")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**DESeq2:**")
            st.code("""
log2FoldChange → log2_fold_change
pvalue → p_value
padj → adjusted_p_value
baseMean → base_mean
            """, language='text')
            
            st.markdown("**edgeR:**")
            st.code("""
logFC → log2_fold_change
PValue → p_value
FDR → adjusted_p_value
logCPM → base_mean
            """, language='text')
        
        with col2:
            st.markdown("**limma-voom:**")
            st.code("""
logFC → log2_fold_change
P.Value → p_value
adj.P.Val → adjusted_p_value
AveExpr → base_mean
            """, language='text')
            
            st.markdown("**Wilcoxon:**")
            st.code("""
logFC → log2_fold_change
pvalue → p_value
padj → adjusted_p_value
            """, language='text')
    
    # Show detailed errors if any
    invalid_files = [r for r in detection_results if not r['valid']]
    if invalid_files:
        with st.expander("Invalid files - click for details"):
            for file_info in invalid_files:
                st.error(f"""
                **{file_info['name']}**
                
                Error: {file_info['error']}
                
                Detected columns: {', '.join(file_info['columns']) if file_info['columns'] else 'Could not read'}
                
                **Troubleshooting:**
                - Ensure file has required columns for the tool you used
                - Check that gene IDs are in first column or as index
                - Verify file format (CSV vs TSV)
                """)
    
    # Section 3: Significance Thresholds
    st.markdown("---")
    st.markdown("## 3. Significance Thresholds")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        fdr_threshold = st.number_input(
            "FDR Threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
            help="Adjusted p-value cutoff (default: 0.05)"
        )
    
    with col2:
        lfc_threshold = st.number_input(
            "Log2FC Threshold",
            min_value=0.0,
            max_value=10.0,
            value=0.0,
            step=0.5,
            help="Minimum |log2FC| for significance (default: 0 = no fold change filter)"
        )
    
    with col3:
        st.markdown("**Significance Criteria:**")
        st.caption(f"• FDR < {fdr_threshold}")
        st.caption(f"• |Log2FC| > {lfc_threshold}")
        st.caption("• Both must be TRUE")
    
    # Section 4: Data Preview
    st.markdown("---")
    st.markdown("## 4. Data Preview & Quality Check")
    
    if file_data:
        selected_file = st.selectbox(
            "Select file to preview:",
            list(file_data.keys()),
            help="Preview raw data before import"
        )
        
        if selected_file:
            preview_df = file_data[selected_file]
            file_info = next(r for r in detection_results if r['name'] == selected_file)
            
            # Show file info
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Genes", f"{len(preview_df):,}")
            
            with col2:
                st.metric("Detected Tool", file_info['tool'])
            
            with col3:
                # Count potential significant genes
                padj_cols = [c for c in preview_df.columns 
                           if any(x in c.lower() for x in ['padj', 'fdr', 'adj.p.val'])]
                if padj_cols:
                    n_sig = (preview_df[padj_cols[0]] < fdr_threshold).sum()
                    st.metric(f"Potential Sig.", f"{n_sig:,}")
                else:
                    st.metric("Potential Sig.", "N/A")
            
            with col4:
                # Calculate mean |log2FC|
                lfc_cols = [c for c in preview_df.columns 
                          if any(x in c.lower() for x in ['logfc', 'log2fc', 'log2foldchange'])]
                if lfc_cols:
                    mean_lfc = preview_df[lfc_cols[0]].abs().mean()
                    st.metric("Mean |Log2FC|", f"{mean_lfc:.2f}")
                else:
                    st.metric("Mean |Log2FC|", "N/A")
            
            # Data table
            st.markdown("**Raw Data Preview:**")
            st.dataframe(preview_df.head(20), use_container_width=True)
            
            # Quick visualization
            with st.expander("Quick Visualization"):
                if padj_cols and lfc_cols:
                    # Create volcano plot preview
                    plot_df = preview_df[[lfc_cols[0], padj_cols[0]]].copy()
                    plot_df.columns = ['logFC', 'padj']
                    plot_df = plot_df.dropna()
                    plot_df['-log10(padj)'] = -np.log10(plot_df['padj'] + 1e-300)
                    
                    # Determine significance
                    plot_df['Significant'] = (
                        (plot_df['padj'] < fdr_threshold) & 
                        (plot_df['logFC'].abs() > lfc_threshold)
                    )
                    
                    fig = px.scatter(
                        plot_df,
                        x='logFC',
                        y='-log10(padj)',
                        color='Significant',
                        color_discrete_map={True: 'red', False: 'gray'},
                        opacity=0.6,
                        title=f"Volcano Plot Preview - {selected_file}"
                    )
                    
                    # Add threshold lines
                    fig.add_hline(
                        y=-np.log10(fdr_threshold),
                        line_dash="dash",
                        line_color="blue",
                        annotation_text=f"FDR={fdr_threshold}"
                    )
                    
                    if lfc_threshold > 0:
                        fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="blue")
                        fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="blue")
                    
                    fig.update_layout(height=500)
                    st.plotly_chart(fig, use_container_width=True)
                    
                    st.caption(f"Red points: Significant (FDR<{fdr_threshold}, |Log2FC|>{lfc_threshold})")
                else:
                    st.info("Cannot create volcano plot - missing required columns")
    
    # Section 5: Import & Standardize
    st.markdown("---")
    st.markdown("## 5. Import & Standardize All Files")
    
    st.info("""
    **Import Process:**
    1. Detect pipeline for each file
    2. Standardize column names
    3. Calculate significance based on your thresholds
    4. Create DEResult objects
    5. Store for Module 9 (Ensemble Analysis)
    """)
    
    valid_files = [r['name'] for r in detection_results if r['valid']]
    
    if not valid_files:
        st.warning("No valid files to import. Please check file formats above.")
    else:
        st.success(f"Ready to import {len(valid_files)} valid file(s)")
        
        if st.button("Import & Standardize All Files", type="primary", use_container_width=True):
            
            with st.spinner(f"Importing and standardizing {len(valid_files)} file(s)..."):
                
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                try:
                    de_results_dict = {}
                    import_errors = []
                    
                    for idx, file_name in enumerate(valid_files):
                        try:
                            status_text.text(f"Processing {file_name}... ({idx+1}/{len(valid_files)})")
                            progress_bar.progress((idx + 0.5) / len(valid_files))
                            
                            df = file_data[file_name]
                            
                            # Detect pipeline
                            pipeline = detect_pipeline(df)
                            
                            # Standardize columns
                            df_standard = standardize_columns(df, pipeline, gene_id_column=None)
                            
                            # Calculate significance
                            df_standard = calculate_significance(
                                df_standard,
                                fdr_threshold=fdr_threshold,
                                lfc_threshold=lfc_threshold
                            )
                            
                            # Create DEResult object
                            de_result = DEResult(
                                results_df=df_standard,
                                pipeline=pipeline.upper(),
                                parameters={
                                    'fdr_threshold': fdr_threshold,
                                    'lfc_threshold': lfc_threshold,
                                    'source_pipeline': pipeline
                                },
                                metadata={
                                    'source_file': file_name,
                                    'n_genes': len(df_standard),
                                    'module': 'M7',
                                    'dashboard_import': True,
                                    'import_timestamp': pd.Timestamp.now().isoformat()
                                }
                            )
                            
                            de_results_dict[file_name] = de_result
                            progress_bar.progress((idx + 1) / len(valid_files))
                            
                        except Exception as e:
                            import_errors.append(f"**{file_name}**: {str(e)}")
                            continue
                    
                    # Store results
                    if de_results_dict:
                        st.session_state['m7_de_results'] = de_results_dict
                        st.session_state['m7_complete'] = True
                        st.session_state['m7_fdr_threshold'] = fdr_threshold
                        st.session_state['m7_lfc_threshold'] = lfc_threshold
                        st.session_state['m7_running'] = False
                        st.session_state['m7_error'] = False
                        
                        progress_bar.progress(1.0)
                        status_text.text("Import complete!")
                        
                        st.success(f"""
                        **Successfully imported {len(de_results_dict)} file(s)!**
                        
                        Files processed:
                        {chr(10).join(f'  • {name}' for name in de_results_dict.keys())}
                        
                        Scroll down to see results summary!
                        """)
                        st.balloons()
                        
                        if import_errors:
                            with st.expander("Some files had errors"):
                                for error in import_errors:
                                    st.warning(error)
                    else:
                        st.error("No files were successfully imported!")
                        if import_errors:
                            st.markdown("**Errors:**")
                            for error in import_errors:
                                st.error(error)
                        st.session_state['m7_error'] = True
                    
                except Exception as e:
                    st.error(f"""
                    **Error during import**
                    
                    {str(e)}
                    
                    **Troubleshooting:**
                    - Check file formats above
                    - Ensure gene IDs are consistent
                    - Verify required columns exist
                    """)
                    st.session_state['m7_error'] = True

else:
    st.info("Upload one or more DE result files to begin")

# Section 6: Results Summary & Analysis
if st.session_state.get('m7_complete', False):
    st.markdown("---")
    st.markdown("## 6️⃣ Results Summary & Analysis")
    
    de_results = st.session_state['m7_de_results']
    fdr_threshold = st.session_state.get('m7_fdr_threshold', 0.05)
    lfc_threshold = st.session_state.get('m7_lfc_threshold', 0.0)
    
    # Overall metrics
    st.markdown("### Overall Metrics")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Files Imported", len(de_results))
    
    with col2:
        unique_pipelines = len(set(r.pipeline for r in de_results.values()))
        st.metric("Unique Pipelines", unique_pipelines)
        if unique_pipelines >= 2:
            st.caption("✅ Good for ensemble")
        else:
            st.caption("⚠️ Consider more methods")
    
    with col3:
        total_sig = sum(r.n_significant for r in de_results.values())
        st.metric("Total Significant", f"{total_sig:,}")
    
    with col4:
        avg_sig_pct = np.mean([r.n_significant / max(r.n_genes, 1) * 100 
                               for r in de_results.values()])
        st.metric("Avg % Significant", f"{avg_sig_pct:.1f}%")
    
    # Per-file summary table
    st.markdown("### Per-File Summary")
    
    summary_data = []
    for file_name, de_result in de_results.items():
        n_up = int((de_result.results_df.get('direction', pd.Series()) == 'up').sum()) if 'direction' in de_result.results_df.columns else 0
        n_down = int((de_result.results_df.get('direction', pd.Series()) == 'down').sum()) if 'direction' in de_result.results_df.columns else 0
        pct_sig = de_result.n_significant / max(de_result.n_genes, 1) * 100
        summary_data.append({
            'File': file_name,
            'Pipeline': de_result.pipeline,
            'Total Genes': f"{de_result.n_genes:,}",
            'Significant': f"{de_result.n_significant:,}",
            'Up-regulated': f"{n_up:,}",
            'Down-regulated': f"{n_down:,}",
            '% Significant': f"{pct_sig:.1f}%"
        })
    
    summary_df = pd.DataFrame(summary_data)
    st.dataframe(summary_df, use_container_width=True, hide_index=True)
    
    # Visualizations
    st.markdown("### Comparative Visualizations")
    
    tab1, tab2, tab3 = st.tabs(["Significant Genes", "Overlap Analysis", "Distribution Comparison"])
    
    with tab1:
        # Bar chart of significant genes
        fig = go.Figure()
        
        file_names = list(de_results.keys())
        n_sig = [de_results[f].n_significant for f in file_names]
        n_up = [de_results[f].n_up for f in file_names]
        n_down = [de_results[f].n_down for f in file_names]
        
        fig.add_trace(go.Bar(
            name='Up-regulated',
            x=file_names,
            y=n_up,
            marker_color='red'
        ))
        
        fig.add_trace(go.Bar(
            name='Down-regulated',
            x=file_names,
            y=n_down,
            marker_color='blue'
        ))
        
        fig.update_layout(
            title="Significant DEGs per File",
            xaxis_title="File",
            yaxis_title="Number of Genes",
            barmode='stack',
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        if len(de_results) >= 2:
            # Gene overlap analysis
            st.markdown("**Significant Gene Overlap Between Files:**")
            
            # Get significant genes from each file
            sig_genes = {}
            for file_name, de_result in de_results.items():
                sig_genes[file_name] = set(de_result.significant_genes)
            
            # Calculate pairwise overlaps
            file_names = list(sig_genes.keys())
            n_files = len(file_names)
            
            # Create overlap matrix
            overlap_matrix = np.zeros((n_files, n_files))
            jaccard_matrix = np.zeros((n_files, n_files))
            
            for i in range(n_files):
                for j in range(n_files):
                    overlap = len(sig_genes[file_names[i]] & sig_genes[file_names[j]])
                    union = len(sig_genes[file_names[i]] | sig_genes[file_names[j]])
                    
                    overlap_matrix[i, j] = overlap
                    jaccard_matrix[i, j] = overlap / union if union > 0 else 0
            
            # Plot heatmap
            fig = go.Figure(data=go.Heatmap(
                z=jaccard_matrix,
                x=file_names,
                y=file_names,
                colorscale='RdYlGn',
                text=overlap_matrix.astype(int),
                texttemplate='%{text}',
                colorbar=dict(title="Jaccard<br>Index")
            ))
            
            fig.update_layout(
                title="Gene Overlap Analysis (Jaccard Index)",
                height=500,
                xaxis={'side': 'bottom'},
                yaxis={'side': 'left'}
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            st.caption("""
            **Jaccard Index** = Overlap / Union
            - 1.0 = Perfect agreement
            - 0.5 = Moderate overlap
            - 0.0 = No overlap
            
            Numbers in cells = Actual gene overlap count
            """)
            
            # Venn diagram info (for 2-3 files)
            if len(de_results) == 2:
                files = list(sig_genes.keys())
                set1, set2 = sig_genes[files[0]], sig_genes[files[1]]
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric(f"{files[0]} only", len(set1 - set2))
                with col2:
                    st.metric("Overlap", len(set1 & set2))
                with col3:
                    st.metric(f"{files[1]} only", len(set2 - set1))
            
        else:
            st.info("Upload at least 2 files to see overlap analysis")
    
    with tab3:
        # Distribution comparison
        st.markdown("**Log2FC Distribution Comparison:**")
        
        fig = go.Figure()
        
        for file_name, de_result in de_results.items():
            sig_genes = de_result.results_df[de_result.results_df['is_significant']]
            
            fig.add_trace(go.Violin(
                y=sig_genes['log2_fold_change'],
                name=file_name,
                box_visible=True,
                meanline_visible=True
            ))
        
        fig.update_layout(
            title="Log2FC Distribution (Significant Genes Only)",
            yaxis_title="Log2 Fold Change",
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    # Export options
    st.markdown("---")
    st.markdown("### Export Options")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**📥 Combined Results (CSV)**")
        # Combine all standardized results
        combined_dfs = []
        for file_name, de_result in de_results.items():
            df = de_result.results_df.copy()
            df['source_file'] = file_name
            df['pipeline'] = de_result.pipeline
            combined_dfs.append(df)
        
        combined_df = pd.concat(combined_dfs, axis=0)
        csv = combined_df.to_csv(index=True)
        
        st.download_button(
            label="Download Combined CSV",
            data=csv,
            file_name="raptor_m7_combined_results.csv",
            mime="text/csv",
            use_container_width=True,
            help="All results in standardized format"
        )
    
    with col2:
        st.markdown("**📊 Summary Report (JSON)**")
        summary_dict = {
            'metadata': {
                'n_files': len(de_results),
                'fdr_threshold': fdr_threshold,
                'lfc_threshold': lfc_threshold,
                'pipelines': list(set(r.pipeline for r in de_results.values()))
            },
            'files': {}
        }
        
        for file_name, de_result in de_results.items():
            summary_dict['files'][file_name] = de_result.calculate_metrics()
        
        json_str = json.dumps(summary_dict, indent=2, default=str)
        
        st.download_button(
            label="Download Summary JSON",
            data=json_str,
            file_name="raptor_m7_summary.json",
            mime="application/json",
            use_container_width=True,
            help="Metrics and statistics"
        )
    
    with col3:
        st.markdown("**🔬 Continue to Ensemble**")
        if len(de_results) >= 2:
            st.success("Ready for ensemble")
            if st.button("➡️ Go to Module 9", type="primary", use_container_width=True):
                # Try to navigate to ensemble page
                try:
                    st.switch_page("pages/05_🔬_Ensemble.py")
                except:
                    st.info("Navigate manually to Module 9 (Ensemble Analysis)")
        else:
            st.warning("2+ files recommended")
            st.caption("You can continue with 1 file, but ensemble works best with multiple methods")
            if st.button("➡️ Continue Anyway", use_container_width=True):
                try:
                    st.switch_page("pages/05_🔬_Ensemble.py")
                except:
                    st.info("Navigate manually to Module 9 (Ensemble Analysis)")
    
    # Next steps info
    st.markdown("---")
    st.info(f"""
    **Import complete! {len(de_results)} file(s) ready for ensemble analysis.**
    
    **What's next:**
    1. Review the summary and visualizations above
    2. Download exports if needed for documentation
    3. Proceed to **Module 9 (Ensemble Analysis)** to:
       - Combine evidence across methods
       - Calculate consensus scores
       - Identify robust DEGs
       - Generate final biomarker lists
    
    **Thresholds used:**
    - FDR < {fdr_threshold}
    - |Log2FC| > {lfc_threshold}
    """)

# Footer
st.markdown("---")
st.caption("**Tip:** For robust results, use 2+ methods (e.g., DESeq2 + edgeR)")
st.caption("**Next:** Module 9 - Ensemble Analysis combines evidence across methods")
