"""
RAPTOR Dashboard - Visualization Suite
Publication-ready interactive plots for RNA-seq differential expression analysis.

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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

st.set_page_config(page_title="RAPTOR - Visualizations", page_icon="V", layout="wide")
init_session_state()
render_sidebar()

st.title("Visualization Suite")
st.caption("Publication-ready interactive plots for differential expression analysis")

with st.expander("\u2139\ufe0f How to use"):
    st.markdown("""
**Available visualizations depend on your data:**

**From DE results (Module 7):** Volcano Plot, MA Plot

**From count matrix (upload or Module 2):** PCA (2D/3D), Heatmap with clustering,
Gene Expression Profiles (box/violin/strip/bar), Sample Correlation

Upload a count matrix and metadata to unlock all visualizations.
All plots are interactive \u2014 zoom, pan, hover for gene names, and export as PNG/SVG.
    """)

st.markdown("---")

# ============================================================================
# DATA LOADING
# ============================================================================
st.markdown("## 1. Data Sources")

has_de = st.session_state.get('m7_complete', False)
de_results_dict = st.session_state.get('m7_de_results', {}) if has_de else {}
counts_from_session = st.session_state.get('m2_count_data')
meta_from_session = st.session_state.get('m2_metadata')

col_de, col_counts = st.columns(2)
selected_de = None
de_data = None

with col_de:
    if has_de and de_results_dict:
        st.success(f"DE results available \u2014 {len(de_results_dict)} file(s)")
        selected_de = st.selectbox("Select DE result:", list(de_results_dict.keys()))
        de_result = de_results_dict[selected_de]
        de_data = de_result.results_df.copy()
        if 'gene_id' not in de_data.columns and de_data.index.name:
            de_data = de_data.reset_index()
    else:
        st.info("No DE results. Run Module 7 for Volcano/MA plots.")

with col_counts:
    if counts_from_session is not None:
        st.success(f"Count matrix from Module 2 \u2014 {counts_from_session.shape[0]:,} genes x {counts_from_session.shape[1]} samples")
    counts_upload = st.file_uploader("Or upload count matrix (CSV)", type=['csv'], key="viz_counts")
    meta_upload = st.file_uploader("Upload metadata (CSV, optional)", type=['csv'], key="viz_meta")

counts_data = None
metadata = None

if counts_upload is not None:
    try:
        counts_data = pd.read_csv(counts_upload, index_col=0)
        st.success(f"Loaded: {counts_data.shape[0]:,} genes x {counts_data.shape[1]} samples")
    except Exception as e:
        st.error(f"Error: {e}")
elif counts_from_session is not None:
    counts_data = counts_from_session

if meta_upload is not None:
    try:
        metadata = pd.read_csv(meta_upload)
    except:
        pass
elif meta_from_session is not None:
    metadata = meta_from_session

if metadata is not None:
    for id_col in ['sample_id', 'sample', 'Sample']:
        if id_col in metadata.columns:
            metadata = metadata.set_index(id_col)
            break

st.markdown("---")

# ============================================================================
# SETTINGS
# ============================================================================
st.markdown("## 2. Settings")
set_col1, set_col2, set_col3, set_col4 = st.columns(4)

with set_col1:
    fdr_thresh = st.slider("FDR threshold", 0.001, 0.20, 0.05, 0.001, format="%.3f")
with set_col2:
    lfc_thresh = st.slider("|Log2FC| threshold", 0.0, 3.0, 1.0, 0.1, format="%.1f")
with set_col3:
    point_size = st.slider("Point size", 2, 12, 5)
with set_col4:
    color_palette = st.selectbox("Color palette", ["Set2", "Dark2", "Pastel1", "Bold", "Viridis"])

palette_map = {
    "Set2": px.colors.qualitative.Set2,
    "Dark2": px.colors.qualitative.Dark2,
    "Pastel1": px.colors.qualitative.Pastel1,
    "Bold": px.colors.qualitative.Bold,
    "Viridis": px.colors.sequential.Viridis
}
colors = palette_map[color_palette]

st.markdown("---")

# ============================================================================
# PLOT TABS
# ============================================================================
st.markdown("## 3. Visualizations")

tab_names = []
if de_data is not None:
    tab_names.extend(["Volcano Plot", "MA Plot"])
if counts_data is not None:
    tab_names.extend(["PCA", "Heatmap", "Gene Expression", "Sample Correlation"])

if not tab_names:
    st.warning("No data available. Import DE results or upload a count matrix.")
    st.stop()

tabs = st.tabs(tab_names)
tab_idx = 0

# ============================================================================
# VOLCANO PLOT
# ============================================================================
if de_data is not None:
    with tabs[tab_idx]:
        df = de_data.copy()
        lfc_col = next((c for c in ['log2_fold_change', 'log2FoldChange', 'logFC'] if c in df.columns), None)
        pval_col = next((c for c in ['adjusted_p_value', 'padj', 'FDR', 'adj.P.Val'] if c in df.columns), None)
        gene_col = next((c for c in ['gene_id', 'Gene', 'gene'] if c in df.columns), None)

        if lfc_col and pval_col:
            df['_lfc'] = df[lfc_col]
            df['_padj'] = df[pval_col]
            df['_neglog10'] = -np.log10(df['_padj'].replace(0, 1e-300).clip(lower=1e-300))
            df['_gene'] = df[gene_col] if gene_col else df.index.astype(str)

            up_mask = (df['_padj'] < fdr_thresh) & (df['_lfc'] >= lfc_thresh)
            down_mask = (df['_padj'] < fdr_thresh) & (df['_lfc'] <= -lfc_thresh)
            df['Status'] = 'Not significant'
            df.loc[up_mask, 'Status'] = f'Up ({up_mask.sum()})'
            df.loc[down_mask, 'Status'] = f'Down ({down_mask.sum()})'

            color_map = {
                f'Up ({up_mask.sum()})': '#d62728',
                f'Down ({down_mask.sum()})': '#1f77b4',
                'Not significant': '#cccccc'
            }

            n_label = st.slider("Label top N genes", 0, 30, 10, key="vol_nlabel")
            sig_genes = df[up_mask | down_mask].nlargest(n_label, '_neglog10')

            fig = px.scatter(
                df, x='_lfc', y='_neglog10', color='Status',
                color_discrete_map=color_map,
                hover_data={'_gene': True, '_lfc': ':.2f', '_padj': ':.2e', '_neglog10': False, 'Status': False},
                labels={'_lfc': 'Log2 Fold Change', '_neglog10': '-Log10(Adjusted P-value)'},
            )
            fig.update_traces(marker=dict(size=point_size, opacity=0.7))
            fig.add_vline(x=lfc_thresh, line_dash="dash", line_color="#888", line_width=1)
            fig.add_vline(x=-lfc_thresh, line_dash="dash", line_color="#888", line_width=1)
            fig.add_hline(y=-np.log10(fdr_thresh), line_dash="dash", line_color="#888", line_width=1)

            if n_label > 0 and len(sig_genes) > 0:
                fig.add_trace(go.Scatter(
                    x=sig_genes['_lfc'], y=sig_genes['_neglog10'],
                    mode='text', text=sig_genes['_gene'],
                    textposition='top center', textfont=dict(size=9, color='black'),
                    showlegend=False
                ))

            fig.update_layout(
                title=f"Volcano Plot \u2014 {selected_de}",
                height=650, plot_bgcolor='#fafafa', font=dict(size=12),
                legend=dict(title='', orientation='h', yanchor='bottom', y=1.02)
            )
            fig.update_xaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
            fig.update_yaxes(showgrid=True, gridcolor='#eee')

            plot_col, stat_col = st.columns([4, 1])
            with plot_col:
                st.plotly_chart(fig, use_container_width=True)
            with stat_col:
                st.markdown("#### Summary")
                st.metric("Total", f"{len(df):,}")
                st.metric("Up", f"{up_mask.sum():,}")
                st.metric("Down", f"{down_mask.sum():,}")
                pct = (up_mask.sum() + down_mask.sum()) / max(len(df), 1) * 100
                st.metric("% DE", f"{pct:.1f}%")
        else:
            st.error("Missing required columns (log2FC and adjusted p-value)")
    tab_idx += 1

# ============================================================================
# MA PLOT
# ============================================================================
if de_data is not None:
    with tabs[tab_idx]:
        df = de_data.copy()
        lfc_col = next((c for c in ['log2_fold_change', 'log2FoldChange', 'logFC'] if c in df.columns), None)
        pval_col = next((c for c in ['adjusted_p_value', 'padj', 'FDR'] if c in df.columns), None)
        mean_col = next((c for c in ['base_mean', 'baseMean', 'logCPM', 'AveExpr'] if c in df.columns), None)

        if lfc_col and mean_col:
            df['_lfc'] = df[lfc_col]
            df['_mean'] = np.log10(df[mean_col] + 1) if df[mean_col].max() > 100 else df[mean_col]
            df['Status'] = 'Not significant'
            if pval_col:
                sig = (df[pval_col] < fdr_thresh) & (df['_lfc'].abs() >= lfc_thresh)
                df.loc[sig & (df['_lfc'] > 0), 'Status'] = 'Up'
                df.loc[sig & (df['_lfc'] < 0), 'Status'] = 'Down'

            fig = px.scatter(
                df, x='_mean', y='_lfc', color='Status',
                color_discrete_map={'Up': '#d62728', 'Down': '#1f77b4', 'Not significant': '#cccccc'},
                labels={'_mean': 'Log10(Mean Expression)', '_lfc': 'Log2 Fold Change'},
            )
            fig.update_traces(marker=dict(size=point_size, opacity=0.6))
            fig.add_hline(y=0, line_color='black', line_width=1)
            fig.add_hline(y=lfc_thresh, line_dash="dash", line_color="#888")
            fig.add_hline(y=-lfc_thresh, line_dash="dash", line_color="#888")

            fig.update_layout(
                title=f"MA Plot \u2014 {selected_de}",
                height=600, plot_bgcolor='#fafafa', font=dict(size=12),
                legend=dict(title='', orientation='h', yanchor='bottom', y=1.02)
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.error("Missing required columns (log2FC and mean expression)")
    tab_idx += 1

# ============================================================================
# PCA
# ============================================================================
if counts_data is not None:
    with tabs[tab_idx]:
        log_data = np.log2(counts_data + 1)
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(log_data.T)
        n_comp = min(6, data_scaled.shape[0], data_scaled.shape[1])
        pca_model = PCA(n_components=n_comp)
        pca_result = pca_model.fit_transform(data_scaled)
        pc_cols = [f'PC{i+1}' for i in range(n_comp)]
        var_pct = {f'PC{i+1}': pca_model.explained_variance_ratio_[i]*100 for i in range(n_comp)}
        var_labels = {pc: f'{pc} ({var_pct[pc]:.1f}%)' for pc in pc_cols}
        pca_df = pd.DataFrame(pca_result, columns=pc_cols, index=log_data.columns)

        color_options = ['None']
        if metadata is not None:
            for col in metadata.columns:
                if metadata[col].nunique() <= 20:
                    pca_df[col] = metadata[col].reindex(pca_df.index)
                    color_options.append(col)

        ctrl1, ctrl2, ctrl3, ctrl4 = st.columns(4)
        with ctrl1:
            pca_color = st.selectbox("Color by", color_options, index=min(1, len(color_options)-1), key="pca_color")
        with ctrl2:
            pca_x = st.selectbox("X axis", pc_cols, index=0, key="pca_x")
        with ctrl3:
            pca_y = st.selectbox("Y axis", pc_cols, index=min(1, n_comp-1), key="pca_y")
        with ctrl4:
            pca_dim = st.selectbox("Dimension", ["2D", "3D"] if n_comp >= 3 else ["2D"], key="pca_dim")

        if pca_dim == "3D":
            pca_z = st.selectbox("Z axis", pc_cols, index=min(2, n_comp-1), key="pca_z")

        color_arg = pca_color if pca_color != 'None' else None

        if pca_dim == "2D":
            fig = px.scatter(pca_df, x=pca_x, y=pca_y, color=color_arg,
                            color_discrete_sequence=colors, text=pca_df.index)
            fig.update_traces(marker=dict(size=12, line=dict(width=1.5, color='white')),
                            textposition='top center', textfont=dict(size=9))
            fig.update_layout(title=f"PCA \u2014 {var_labels[pca_x]} vs {var_labels[pca_y]}",
                            xaxis_title=var_labels[pca_x], yaxis_title=var_labels[pca_y],
                            height=650, plot_bgcolor='#fafafa', font=dict(size=12))
            fig.update_xaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
            fig.update_yaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
        else:
            fig = px.scatter_3d(pca_df, x=pca_x, y=pca_y, z=pca_z, color=color_arg,
                               color_discrete_sequence=colors, text=pca_df.index)
            fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color='white')))
            fig.update_layout(title=f"3D PCA \u2014 {var_labels[pca_x]} vs {var_labels[pca_y]} vs {var_labels[pca_z]}",
                            height=750, font=dict(size=12),
                            scene=dict(xaxis_title=var_labels[pca_x], yaxis_title=var_labels[pca_y], zaxis_title=var_labels[pca_z]))

        st.plotly_chart(fig, use_container_width=True)

        with st.expander("Variance explained (Scree plot)"):
            scree = pd.DataFrame({'Component': pc_cols, 'Variance (%)': [var_pct[pc] for pc in pc_cols],
                                  'Cumulative (%)': list(np.cumsum([var_pct[pc] for pc in pc_cols]))})
            fig_s = go.Figure()
            fig_s.add_bar(x=scree['Component'], y=scree['Variance (%)'], name='Individual', marker_color='#66c2a5')
            fig_s.add_scatter(x=scree['Component'], y=scree['Cumulative (%)'], name='Cumulative', mode='lines+markers', marker_color='#fc8d62')
            fig_s.update_layout(yaxis_title='Variance (%)', height=300, plot_bgcolor='#fafafa')
            st.plotly_chart(fig_s, use_container_width=True)
    tab_idx += 1

    # ========================================================================
    # HEATMAP
    # ========================================================================
    with tabs[tab_idx]:
        hm1, hm2, hm3 = st.columns(3)
        with hm1:
            n_genes_hm = st.slider("Top variable genes", 20, 200, 50, 10, key="hm_ngenes")
        with hm2:
            hm_cs = st.selectbox("Color scale", ["RdBu_r", "Viridis", "Plasma", "PiYG", "BrBG"], key="hm_cs")
        with hm3:
            hm_cluster = st.checkbox("Cluster rows and columns", value=True, key="hm_cluster")

        log_c = np.log2(counts_data + 1)
        gene_var = log_c.var(axis=1).sort_values(ascending=False)
        top_genes = gene_var.head(n_genes_hm).index
        subset = log_c.loc[top_genes]
        row_m = subset.mean(axis=1)
        row_s = subset.std(axis=1).replace(0, 1)
        subset_z = subset.subtract(row_m, axis=0).divide(row_s, axis=0)

        if hm_cluster and len(top_genes) > 2 and subset_z.shape[1] > 2:
            r_link = linkage(pdist(subset_z.values, 'correlation'), method='average')
            c_link = linkage(pdist(subset_z.T.values, 'correlation'), method='average')
            subset_z = subset_z.iloc[leaves_list(r_link), leaves_list(c_link)]

        fig = go.Figure(go.Heatmap(
            z=subset_z.values, x=subset_z.columns.tolist(), y=subset_z.index.tolist(),
            colorscale=hm_cs, zmid=0,
            hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>',
            colorbar=dict(title=dict(text="Z-score"), len=0.7)
        ))
        fig.update_layout(title=f"Heatmap \u2014 Top {n_genes_hm} Variable Genes (Z-score)",
                         height=max(500, n_genes_hm * 12), font=dict(size=11), plot_bgcolor='white',
                         yaxis=dict(tickfont=dict(size=max(6, min(10, 500 // n_genes_hm)))))
        st.plotly_chart(fig, use_container_width=True)

        if metadata is not None:
            with st.expander("Sample annotations"):
                anno_cols = [c for c in metadata.columns if metadata[c].nunique() <= 10]
                if anno_cols:
                    sel_anno = st.multiselect("Show:", anno_cols, default=anno_cols[:2])
                    for ac in sel_anno:
                        vals = metadata[ac].reindex(subset_z.columns)
                        uv = sorted(vals.dropna().unique())
                        nv = max(len(uv) - 1, 1)
                        v2n = {v: i for i, v in enumerate(uv)}
                        fig_a = go.Figure(go.Heatmap(
                            z=[[v2n.get(v, -1) for v in vals]], x=subset_z.columns.tolist(), y=[ac],
                            colorscale=[[i/nv, colors[i % len(colors)]] for i in range(len(uv))],
                            showscale=False, hovertemplate='%{x}: %{text}<extra></extra>',
                            text=[[str(v) for v in vals]]
                        ))
                        fig_a.update_layout(height=80, margin=dict(l=100, r=20, t=5, b=5),
                                           xaxis=dict(showticklabels=False))
                        st.plotly_chart(fig_a, use_container_width=True)
                        st.caption(f"{ac}: " + ", ".join(str(v) for v in uv))

        with st.expander("Highlight specific genes"):
            gs = st.text_input("Gene names (comma-separated):", key="hm_search")
            if gs:
                sg = [g.strip() for g in gs.split(',')]
                found = [g for g in sg if g in subset_z.index]
                nf = [g for g in sg if g not in subset_z.index]
                if found:
                    st.success(f"Found: {', '.join(found)}")
                    st.dataframe(subset_z.loc[found], use_container_width=True)
                if nf:
                    st.warning(f"Not in top {n_genes_hm}: {', '.join(nf)}")
    tab_idx += 1

    # ========================================================================
    # GENE EXPRESSION
    # ========================================================================
    with tabs[tab_idx]:
        st.markdown("Explore individual gene expression across samples and conditions.")
        gene_input = st.text_input("Gene ID(s) \u2014 comma-separated:", key="expr_genes")

        with st.expander("Browse genes"):
            gl = counts_data.index.tolist()
            sf = st.text_input("Filter:", key="gene_filter")
            if sf:
                gl = [g for g in gl if sf.lower() in g.lower()]
            st.dataframe(pd.DataFrame({'Gene': gl[:500]}), height=200)

        if gene_input:
            qg = [g.strip() for g in gene_input.split(',')]
            vg = [g for g in qg if g in counts_data.index]
            ig = [g for g in qg if g not in counts_data.index]
            if ig:
                st.warning(f"Not found: {', '.join(ig)}")
            if vg:
                le = np.log2(counts_data.loc[vg] + 1)
                ec1, ec2 = st.columns(2)
                with ec1:
                    pt = st.selectbox("Plot type", ["Box Plot", "Violin Plot", "Strip Plot", "Bar Plot", "Box + Strip"], key="expr_type")
                with ec2:
                    gc = None
                    if metadata is not None:
                        go2 = ['None'] + [c for c in metadata.columns if metadata[c].nunique() <= 10]
                        gs2 = st.selectbox("Group by", go2, index=min(1, len(go2)-1), key="expr_group")
                        if gs2 != 'None':
                            gc = gs2

                for gene in vg:
                    expr = le.loc[gene]
                    pdf = pd.DataFrame({'Sample': expr.index, 'Expression': expr.values})
                    if gc and metadata is not None:
                        pdf['Group'] = metadata[gc].reindex(pdf['Sample'].values).values
                    ca = 'Group' if 'Group' in pdf.columns else None

                    if pt == "Box Plot":
                        fig = px.box(pdf, x='Group' if ca else 'Sample', y='Expression', color=ca,
                                    color_discrete_sequence=colors, points='all')
                    elif pt == "Violin Plot":
                        fig = px.violin(pdf, x='Group' if ca else 'Sample', y='Expression', color=ca,
                                       color_discrete_sequence=colors, box=True, points='all')
                    elif pt == "Strip Plot":
                        fig = px.strip(pdf, x='Group' if ca else 'Sample', y='Expression', color=ca,
                                      color_discrete_sequence=colors)
                    elif pt == "Bar Plot":
                        fig = px.bar(pdf, x='Sample', y='Expression', color=ca, color_discrete_sequence=colors)
                    elif pt == "Box + Strip":
                        fig = go.Figure()
                        if ca:
                            for grp in pdf['Group'].unique():
                                gd = pdf[pdf['Group'] == grp]
                                fig.add_trace(go.Box(y=gd['Expression'], name=str(grp), boxpoints='all', jitter=0.3, pointpos=-1.5, marker=dict(size=6)))
                        else:
                            fig.add_trace(go.Box(y=pdf['Expression'], x=pdf['Sample'], boxpoints='all', jitter=0.3))

                    fig.update_layout(title=f"{gene} \u2014 Log2(Count + 1)", yaxis_title="Log2(Expression + 1)",
                                    height=450, plot_bgcolor='#fafafa', font=dict(size=12), showlegend=bool(ca))
                    fig.update_xaxes(tickangle=45 if not ca else 0)
                    st.plotly_chart(fig, use_container_width=True)

                if len(vg) > 1:
                    st.markdown("#### Multi-Gene Comparison")
                    melted = le.T.melt(var_name='Gene', value_name='Expression', ignore_index=False)
                    melted['Sample'] = melted.index
                    fig = px.box(melted, x='Gene', y='Expression', color='Gene', color_discrete_sequence=colors, points='all')
                    fig.update_layout(title="Expression Comparison", yaxis_title="Log2(Expression + 1)",
                                    height=450, plot_bgcolor='#fafafa', font=dict(size=12))
                    st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Enter a gene ID above to view expression profiles.")
    tab_idx += 1

    # ========================================================================
    # SAMPLE CORRELATION
    # ========================================================================
    with tabs[tab_idx]:
        cc1, cc2, cc3 = st.columns(3)
        with cc1:
            cm = st.selectbox("Method", ["pearson", "spearman", "kendall"], key="viz_corr_m")
        with cc2:
            ccs = st.selectbox("Color scale", ["Viridis", "RdBu_r", "Plasma", "Teal-Orange"], key="viz_corr_cs")
        with cc3:
            ccl = st.checkbox("Cluster samples", value=True, key="viz_corr_cl")

        lc = np.log2(counts_data + 1)
        corr = lc.corr(method=cm)

        if ccl and corr.shape[0] > 2:
            d = pdist(corr.values, 'correlation')
            Z = linkage(d, method='average')
            o = leaves_list(Z)
            corr = corr.iloc[o, o]

        csm = {"Viridis": "Viridis", "RdBu_r": "RdBu_r", "Plasma": "Plasma",
               "Teal-Orange": [[0, '#009392'], [0.5, '#f6edbd'], [1, '#d0587e']]}
        od = corr.values[~np.eye(len(corr), dtype=bool)]
        zmin = max(0, np.floor(od.min() * 20) / 20)

        fig = go.Figure(go.Heatmap(
            z=corr.values, x=corr.columns.tolist(), y=corr.index.tolist(),
            colorscale=csm[ccs], zmin=zmin, zmax=1.0,
            text=np.round(corr.values, 2), texttemplate='%{text}',
            textfont=dict(size=max(6, min(9, 120 // len(corr)))),
            hovertemplate='%{x} vs %{y}<br>r = %{z:.4f}<extra></extra>',
            colorbar=dict(title=dict(text=cm.capitalize()))
        ))
        fig.update_layout(title=f"Sample Correlation ({cm.capitalize()})", height=650, font=dict(size=12),
                         xaxis=dict(tickangle=45))
        st.plotly_chart(fig, use_container_width=True)
        st.caption(f"Range: {od.min():.3f} \u2013 {od.max():.3f}  |  Mean: {od.mean():.3f}  |  Std: {od.std():.3f}")

        if metadata is not None:
            with st.expander("Overlay sample metadata"):
                for col in [c for c in metadata.columns if metadata[c].nunique() <= 10][:3]:
                    vals = metadata[col].reindex(corr.columns)
                    st.markdown(f"**{col}:** " + " | ".join(f"{s}: {v}" for s, v in zip(corr.columns, vals)))

st.markdown("---")
st.caption("All plots are interactive \u2014 zoom, pan, hover for details. Right-click a plot to save as PNG/SVG.")