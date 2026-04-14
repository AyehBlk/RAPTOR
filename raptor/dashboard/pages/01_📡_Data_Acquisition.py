"""
RAPTOR Dashboard - Data Acquisition & Import

Search public repositories (GEO, TCGA, ArrayExpress, SRA),
upload your own data, manage a local data library, pool datasets
with batch correction, and run quality checks before analysis.

Author: Ayeh Bolouki
Version: 2.2.2
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
from typing import Optional
import json
import logging
from datetime import datetime

logger = logging.getLogger(__name__)

# Add project root to path (4 levels up from pages/)
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# Import RAPTOR acquisition module
try:
    from raptor.external_modules.acquisition import (
        AcquiredDataset,
        PooledDataset,
        CacheManager,
        DataCatalog,
        GEOConnector,
        TCGAConnector,
        ArrayExpConnector,
        SRAConnector,
        GeneIDMapper,
        PoolingEngine,
        get_available_components,
        DEFAULT_CACHE_DIR,
    )
    ACQUISITION_AVAILABLE = True
except ImportError as e:
    ACQUISITION_AVAILABLE = False
    import traceback
    import_error = str(e)
    import_traceback = traceback.format_exc()

# Import QC module for pooled data quality checks
try:
    from raptor.quality_assessment import DataQualityAssessor
    QC_AVAILABLE = True
except ImportError:
    QC_AVAILABLE = False

# Page config
st.set_page_config(
    page_title="RAPTOR - Data Acquisition",
    page_icon="📡",
    layout="wide"
)

# Initialize
init_session_state()
render_sidebar()

# =============================================================================
# PERFORMANCE HELPERS — Cached computation for large datasets
# =============================================================================

N_TOP_GENES_QC = 5000  # Use top variable genes for QC (not all 60K)
_PARQUET_CACHE_DIR = Path.home() / '.raptor' / 'parquet_cache'
_PARQUET_CACHE_DIR.mkdir(parents=True, exist_ok=True)


# --- Data-type-aware labels (for multi-omic TCGA support) ---

def _feature_label(ds, plural=True) -> str:
    """Return the correct feature label based on dataset type (e.g. 'genes', 'probes', 'miRNAs')."""
    gid = getattr(ds, 'gene_id_type', 'unknown')
    labels = {
        'mirna': ('miRNAs', 'miRNA'),
        'probe': ('probes', 'probe'),
        'segment': ('segments', 'segment'),
        'genomic_region': ('genomic regions', 'genomic region'),
        'protein': ('proteins', 'protein'),
    }
    pair = labels.get(gid, ('genes', 'gene'))
    return pair[0] if plural else pair[1]


def _matrix_label(ds) -> str:
    """Return the correct matrix name based on dataset type."""
    gid = getattr(ds, 'gene_id_type', 'unknown')
    labels = {
        'mirna': 'expression matrix',
        'probe': 'beta value matrix',
        'segment': 'segment matrix',
        'genomic_region': 'CNV matrix',
        'protein': 'protein expression matrix',
    }
    return labels.get(gid, 'count matrix')


# --- Parquet I/O (4-10x faster than CSV for large matrices) ---

def _save_parquet(df: pd.DataFrame, name: str) -> Path:
    """Save DataFrame as Parquet for fast reload. Returns path."""
    path = _PARQUET_CACHE_DIR / f"{name}.parquet"
    try:
        df.to_parquet(path, engine='pyarrow', compression='snappy')
    except ImportError:
        df.to_parquet(path, engine='fastparquet', compression='snappy')
    except Exception:
        # Fallback: save as feather (also columnar, fast)
        try:
            path = _PARQUET_CACHE_DIR / f"{name}.feather"
            df.to_feather(path)
        except Exception:
            path = _PARQUET_CACHE_DIR / f"{name}.csv"
            df.to_csv(path)
    return path


def _load_parquet(name: str) -> Optional[pd.DataFrame]:
    """Load DataFrame from Parquet cache. Returns None if not found."""
    for ext in ['.parquet', '.feather', '.csv']:
        path = _PARQUET_CACHE_DIR / f"{name}{ext}"
        if path.exists():
            try:
                if ext == '.parquet':
                    return pd.read_parquet(path)
                elif ext == '.feather':
                    return pd.read_feather(path)
                else:
                    return pd.read_csv(path, index_col=0)
            except Exception:
                continue
    return None


# --- Sparse matrix support (50-70% memory reduction for RNA-seq) ---

def _to_sparse_counts(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert count matrix to sparse-backed DataFrame.
    RNA-seq has ~17-40% zeros — sparse format saves significant memory.
    For 60K x 884 matrix: ~400MB dense → ~150MB sparse.
    """
    try:
        sparse_df = counts_df.astype(pd.SparseDtype("float64", fill_value=0.0))
        return sparse_df
    except Exception:
        return counts_df


def _sparse_sum(df, axis=0):
    """Sum that works on both dense and sparse DataFrames."""
    try:
        return df.sum(axis=axis)
    except Exception:
        return df.to_dense().sum(axis=axis) if hasattr(df, 'to_dense') else df.sum(axis=axis)


def _sparse_var(df, axis=0):
    """Variance that works on both dense and sparse DataFrames."""
    try:
        return df.var(axis=axis)
    except Exception:
        return df.to_dense().var(axis=axis) if hasattr(df, 'to_dense') else df.var(axis=axis)


# --- Lazy loading wrapper ---

class _LazyDataset:
    """
    Lightweight wrapper that keeps count matrix on disk (as Parquet)
    and only loads into memory when accessed.
    
    Usage:
        lazy = _LazyDataset(counts_df, "TCGA-LIHC")
        lazy.counts_df  # loads from Parquet on first access, cached afterward
    """
    def __init__(self, counts_df: pd.DataFrame, name: str):
        self._name = name
        self._counts = None
        self._path = _save_parquet(counts_df, f"lazy_{name}")
        self._shape = counts_df.shape
        self._columns = counts_df.columns.tolist()
        self._index_sample = counts_df.index[:5].tolist()
    
    @property
    def counts_df(self):
        if self._counts is None:
            self._counts = _load_parquet(f"lazy_{self._name}")
        return self._counts
    
    @property
    def shape(self):
        return self._shape
    
    def release(self):
        """Free memory — data will reload from disk on next access."""
        self._counts = None


# --- Cached computations (using top 5K genes for speed) ---

@st.cache_data(show_spinner=False, ttl=600)
def _cached_pca(data_hash, counts_values, n_samples, n_comp=6, gene_id_type='ensembl'):
    """Cached PCA computation on top variable features."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    df = pd.DataFrame(counts_values)

    # Data-type-aware preprocessing
    is_beta = gene_id_type == 'probe'  # methylation beta values (0–1)
    is_count = gene_id_type in ('ensembl', 'symbol', 'entrez', 'refseq', 'unknown')

    # Handle NaN: drop features with >50% missing, impute remainder
    if df.isna().any().any():
        nan_per_row = df.isna().sum(axis=1) / df.shape[1]
        df = df.loc[nan_per_row < 0.5]  # keep probes with <50% NaN
        if is_beta:
            # Impute with row median (biologically neutral for beta values)
            row_medians = df.median(axis=1)
            for col in df.columns:
                mask = df[col].isna()
                if mask.any():
                    df.loc[mask, col] = row_medians[mask]
        else:
            df = df.fillna(0)

    if is_count:
        log_df = np.log2(df + 1)
    else:
        # For beta values, RPM, protein expression: use values directly
        log_df = df

    # Use only top variable features for speed
    gene_var = log_df.var(axis=1)
    top_idx = gene_var.nlargest(min(N_TOP_GENES_QC, len(gene_var))).index
    log_top = log_df.loc[top_idx]

    var_mask = log_top.var(axis=1) > 0
    scaled = StandardScaler().fit_transform(log_top.loc[var_mask].T.values)

    n_comp = min(n_comp, scaled.shape[0], scaled.shape[1])
    pca = PCA(n_components=n_comp)
    result = pca.fit_transform(scaled)
    return result, pca.explained_variance_ratio_, n_comp


@st.cache_data(show_spinner=False, ttl=600)
def _cached_correlation(data_hash, counts_values, columns, method='pearson', gene_id_type='ensembl'):
    """Cached sample correlation on top variable features."""
    df = pd.DataFrame(counts_values, columns=columns)

    # Handle NaN: drop features with >50% missing, fill remainder
    if df.isna().any().any():
        nan_per_row = df.isna().sum(axis=1) / df.shape[1]
        df = df.loc[nan_per_row < 0.5]
        df = df.fillna(df.median(axis=1).reindex(df.index), axis=0)

    gene_var = df.var(axis=1)
    top_genes = gene_var.nlargest(min(N_TOP_GENES_QC, len(gene_var))).index

    is_count = gene_id_type in ('ensembl', 'symbol', 'entrez', 'refseq', 'unknown')
    if is_count:
        log_top = np.log2(df.loc[top_genes] + 1)
    else:
        log_top = df.loc[top_genes]

    return log_top.corr(method=method).values


@st.cache_data(show_spinner=False, ttl=600)
def _cached_lib_stats(data_hash, counts_values, gene_id_type='ensembl'):
    """Cached library size and feature detection stats."""
    df = pd.DataFrame(counts_values)
    is_beta = gene_id_type == 'probe'

    if is_beta:
        # For methylation: "library size" = mean beta, "detected" = non-NaN probes
        lib_sizes = df.mean(axis=0).values  # mean beta per sample
        gene_detected = df.notna().sum(axis=0).values  # probes with signal
        total_cells = df.shape[0] * df.shape[1]
        nan_pct = df.isna().sum().sum() / total_cells * 100 if total_cells > 0 else 0
        return lib_sizes, gene_detected, nan_pct
    else:
        lib_sizes = df.sum(axis=0).values
        gene_detected = (df > 0).sum(axis=0).values
        total_cells = df.shape[0] * df.shape[1]
        zero_pct = (df.values == 0).sum() / total_cells * 100 if total_cells > 0 else 0
        return lib_sizes, gene_detected, zero_pct


def _data_hash(counts_df):
    """Fast hash for cache key — uses shape + corner values, not full data."""
    shape_str = f"{counts_df.shape[0]}x{counts_df.shape[1]}"
    corner = ""
    if counts_df.shape[0] > 0 and counts_df.shape[1] > 0:
        try:
            v1 = counts_df.iloc[0, 0]
            v2 = counts_df.iloc[-1, -1]
            v3 = counts_df.iloc[0, -1]
            corner = f"{v1}_{v2}_{v3}"
        except Exception:
            corner = "na"
    col_hash = hash(tuple(counts_df.columns[:5].tolist() + counts_df.columns[-5:].tolist()))
    return f"{shape_str}_{corner}_{col_hash}"


# --- Publication-quality export configuration ---

PUB_FONT = dict(family="Arial, Helvetica, sans-serif", size=14, color="black")

def _pub_config(filename="RAPTOR_plot"):
    """Plotly chart config for publication-quality SVG/PNG export."""
    return {
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': filename,
            'height': 900,
            'width': 1200,
            'scale': 3,  # 3x resolution for crisp text
        },
        'displayModeBar': True,
        'modeBarButtonsToAdd': ['drawline', 'eraseshape'],
    }


def _pub_layout(fig, title="", width=1000, height=700):
    """Apply publication-quality styling to a Plotly figure."""
    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", font=dict(size=18, family="Arial")),
        font=PUB_FONT,
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=width,
        height=height,
        legend=dict(
            font=dict(size=13, family="Arial"),
            bordercolor="gray", borderwidth=1,
            bgcolor="rgba(255,255,255,0.9)",
        ),
        margin=dict(l=80, r=40, t=80, b=80),
    )
    fig.update_xaxes(
        showline=True, linewidth=1.5, linecolor='black',
        showgrid=True, gridcolor='#E8E8E8', gridwidth=0.5,
        tickfont=dict(size=12, family="Arial"),
        title_font=dict(size=14, family="Arial"),
        zeroline=True, zerolinecolor='#CCCCCC', zerolinewidth=1,
    )
    fig.update_yaxes(
        showline=True, linewidth=1.5, linecolor='black',
        showgrid=True, gridcolor='#E8E8E8', gridwidth=0.5,
        tickfont=dict(size=12, family="Arial"),
        title_font=dict(size=14, family="Arial"),
        zeroline=True, zerolinecolor='#CCCCCC', zerolinewidth=1,
    )
    return fig


def _show_plot(fig, filename="RAPTOR_plot", key_suffix=""):
    """Display a Plotly figure with publication-quality download buttons."""
    st.plotly_chart(fig, use_container_width=True, config=_pub_config(filename))
    dl1, dl2 = st.columns(2)
    with dl1:
        try:
            png_bytes = fig.to_image(format="png", width=1200, height=900, scale=3, engine="kaleido")
            st.download_button(
                "📥 Download PNG (300 DPI)", data=png_bytes,
                file_name=f"{filename}.png", mime="image/png",
                key=f"dl_png_{filename}_{key_suffix}",
            )
        except Exception:
            pass
    with dl2:
        try:
            svg_bytes = fig.to_image(format="svg", width=1200, height=900, scale=3, engine="kaleido")
            st.download_button(
                "📥 Download SVG (vector)", data=svg_bytes,
                file_name=f"{filename}.svg", mime="image/svg+xml",
                key=f"dl_svg_{filename}_{key_suffix}",
            )
        except Exception:
            pass

# Main content
st.title("Data Acquisition & Import")
st.caption("Search public repositories, upload your own data, or pool datasets for downstream analysis")

st.info("""
**No coding required!** This interactive dashboard lets you search thousands of public RNA-seq 
datasets, download them with one click, combine studies from different labs, and check data quality 
— all through this visual interface. Just follow the tabs from left to right.
""")

# Help section
with st.expander("ℹ️ How to use this page (click to expand)"):
    st.markdown("""
    ### Getting Started
    
    This page is your starting point for any RNA-seq analysis in RAPTOR. You can either:
    - **Search** public databases for existing datasets (no data of your own needed!)
    - **Upload** your own count matrix from your RNA-seq experiment
    - **Combine both** — mix your data with public datasets for meta-analysis
    
    ### Step-by-Step Guide
    
    **Step 1 — Search or Upload** (first tab)
    > Type keywords like "breast cancer RNA-seq" or "Alzheimer human" and click Search.
    > You'll see a table of matching datasets. Click **Download** to get one.
    > Or scroll down to upload your own CSV file — no command line needed.
    
    **Step 2 — Data Library** (second tab)
    > All your downloaded and uploaded datasets appear here.
    > You can preview the data, check how many genes and samples each has,
    > and manage what's stored on your computer.
    
    **Step 3 — Pool Datasets** (third tab)
    > Select 2 or more datasets and merge them into one combined dataset.
    > RAPTOR handles gene ID matching and batch correction automatically.
    > This is powerful for finding biomarkers that are consistent across studies.
    
    **Step 4 — Quality Check** (fourth tab)
    > Before continuing your analysis, check if the pooled data looks good.
    > You'll see library sizes, batch effects, PCA plots, and expression distributions.
    > RAPTOR flags any issues and suggests what to do.
    
    **Step 5 — Export** (fifth tab)
    > Download your data as CSV files, or pass it directly to the next
    > dashboard pages (Quality Assessment, Import DE, etc.).
    
    ### Why Pool Your Own Data with Public Studies?
    
    This is one of the most powerful workflows for **biomarker discovery**:
    
    > **The problem:** You run an RNA-seq experiment and find 200 differentially expressed
    > genes. But are these real biomarkers, or artifacts of your specific cohort, your
    > sequencing batch, or your sample demographics?
    >
    > **The solution:** Upload your own count matrix, then pool it with 2–3 public studies
    > of the same disease from GEO or TCGA. RAPTOR harmonizes gene IDs across studies,
    > applies batch correction (ComBat, quantile normalization) to remove technical
    > differences between labs, and merges everything into one large dataset.
    >
    > **The result:** A gene that is consistently differentially expressed across your data
    > AND independent public cohorts is a far more credible biomarker candidate. You get
    > the statistical power of hundreds of samples instead of your original 20–30, and any
    > biomarker that survives batch correction across 3+ labs has already passed a form of
    > external validation.
    >
    > **For publications:** Reviewers increasingly expect cross-cohort validation. Being
    > able to say "we validated our findings across N independent cohorts" strengthens any
    > manuscript — and RAPTOR makes this workflow accessible without custom integration code.
    
    ### Supported Repositories
    
    | Repository | What it has | Best for |
    |---|---|---|
    | **GEO** (NCBI) | 200,000+ datasets | General gene expression studies |
    | **TCGA** (NCI) | 33 cancer types | Cancer research |
    | **ArrayExpress** (EBI) | European studies | Alternative to GEO |
    | **SRA** (NCBI) | Raw sequencing data | When you need FASTQ files |
    
    ### Tips
    - Start with **GEO** — it has the most data and is easiest to use
    - For cancer research, try **TCGA** — every dataset downloads reliably
    - You don't need to install anything extra for basic searching and uploading
    - Pooling 3+ studies from different labs gives the most robust results
    - Upload your own data alongside public datasets for cross-cohort biomarker validation
    """)

# Check module availability
if not ACQUISITION_AVAILABLE:
    st.error(f"""
    **RAPTOR acquisition module not available**
    
    Error: {import_error}
    
    Install dependencies:
    ```bash
    pip install requests pyarrow GEOparse biopython mygene
    ```
    """)
    with st.expander("Full error traceback"):
        st.code(import_traceback)
        st.code(f"Python: {sys.executable}")
        st.code(f"sys.path: {sys.path[:5]}")
    st.stop()

# Show component status
components = get_available_components()
missing = [k for k, v in components.items() if not v]
if missing:
    st.warning(f"Some components unavailable: {', '.join(missing)}. Install optional deps for full functionality.")

st.markdown("---")

# =============================================================================
# Initialize shared objects in session state
# =============================================================================

if 'acq_cache' not in st.session_state:
    st.session_state['acq_cache'] = CacheManager()

if 'acq_catalog' not in st.session_state:
    st.session_state['acq_catalog'] = DataCatalog(DEFAULT_CACHE_DIR)

if 'acq_datasets' not in st.session_state:
    st.session_state['acq_datasets'] = {}

if 'acq_pooled' not in st.session_state:
    st.session_state['acq_pooled'] = None

cache = st.session_state['acq_cache']
catalog = st.session_state['acq_catalog']

# =============================================================================
# TABS
# =============================================================================

tab_search, tab_library, tab_pool, tab_qc, tab_export = st.tabs([
    "🔍 Search Repositories",
    "📚 Data Library",
    "🔗 Pool Datasets",
    "✅ Quality Check",
    "📥 Export",
])

# =============================================================================
# TAB 1: Search Repositories
# =============================================================================

with tab_search:
    st.markdown("## Get Data by Accession")
    st.caption("Enter a GEO, TCGA, ArrayExpress, or SRA accession number")
    
    col_acc, col_preview_btn, col_dl_btn = st.columns([3, 1, 1])
    
    with col_acc:
        direct_accession = st.text_input(
            "Accession number",
            placeholder="e.g., GSE306761, TCGA-LIHC, E-MTAB-1234, SRP123456",
            help="Enter a dataset accession to preview or download"
        )
    
    with col_preview_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        preview_clicked = st.button("Preview Info", use_container_width=True, key='preview_btn')
    
    with col_dl_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        direct_download = st.button("Download", type="primary", use_container_width=True, key='direct_dl')
    
    # Auto-detect repository from accession
    def _detect_repo(accession):
        accession = accession.strip()
        if accession.startswith('GSE') or accession.startswith('GSM') or accession.startswith('GPL'):
            return 'GEO'
        elif accession.startswith('TCGA') or accession.startswith('TARGET'):
            return 'TCGA'
        elif accession.startswith('E-'):
            return 'ArrayExpress'
        elif accession.startswith('SRP') or accession.startswith('ERP') or accession.startswith('DRP'):
            return 'SRA'
        return 'GEO'
    
    def _get_connector(repo):
        if repo == 'GEO':
            return GEOConnector(cache=True, cache_dir=cache.cache_dir)
        elif repo == 'TCGA':
            return TCGAConnector(cache=True, cache_dir=cache.cache_dir)
        elif repo == 'ArrayExpress':
            return ArrayExpConnector(cache=True, cache_dir=cache.cache_dir)
        else:
            return SRAConnector(cache=True, cache_dir=cache.cache_dir)
    
    def _repo_url(accession, repo=None):
        """Get the correct web URL for any repository/accession pair."""
        if repo is None:
            repo = _detect_repo(accession)
        if repo == 'GEO':
            return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}", "GEO"
        elif repo == 'TCGA':
            return f"https://portal.gdc.cancer.gov/projects/{accession}", "GDC Portal"
        elif repo == 'ArrayExpress':
            return f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}", "BioStudies"
        elif repo == 'SRA':
            return f"https://www.ncbi.nlm.nih.gov/sra/?term={accession}", "SRA"
        return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}", "Repository"
    
    def _repo_link_markdown(accession, repo=None):
        """Return a markdown link to the repository for any accession."""
        url, label = _repo_url(accession, repo)
        return f"[Open {accession} in {label}]({url})"
    
    # ---- PREVIEW INFO ----
    if preview_clicked and direct_accession:
        accession = direct_accession.strip()
        repo = _detect_repo(accession)
        
        with st.spinner(f"Fetching info for {accession} from {repo}..."):
            try:
                connector = _get_connector(repo)
                if hasattr(connector, 'get_study_info'):
                    info = connector.get_study_info(accession)
                else:
                    info = connector.info(accession)
                
                if info:
                    st.session_state['preview_info'] = info
                    st.session_state['preview_accession'] = accession
                    st.session_state['preview_repo'] = repo
                else:
                    st.warning(f"No info found for {accession}")
            except Exception as e:
                st.error(f"Failed to fetch info: {e}")
    
    # Display preview if available
    if 'preview_info' in st.session_state and st.session_state.get('preview_accession') == direct_accession.strip():
        info = st.session_state['preview_info']
        
        st.markdown("---")
        st.markdown(f"### {info.get('title', 'No title')}")
        st.markdown(_repo_link_markdown(direct_accession.strip()))
        
        # Key info row (compact, not oversized)
        n_cases = info.get('n_samples', info.get('n_cases', 'N/A'))
        n_samples_total = len(info.get('samples', [])) if isinstance(info.get('samples'), pd.DataFrame) else ''
        organism = info.get('organism', 'unknown')
        platform = info.get('platform_id', 'unknown')
        data_type = info.get('data_type', 'unknown')
        
        cases_str = f"**{n_cases}** cases"
        if n_samples_total and n_samples_total != n_cases:
            cases_str += f" ({n_samples_total} samples)"
        
        st.markdown(
            f"{cases_str} · **{organism}** · {platform} · {data_type}",
        )
        
        # Platform name
        platform_name = info.get('platform_name', '')
        if platform_name:
            st.caption(f"Platform: {platform_name}")
        
        # ---- TCGA-specific: sample type breakdown ----
        sample_type_counts = info.get('sample_type_counts', {})
        if sample_type_counts:
            with st.expander("Sample Type Breakdown", expanded=True):
                stc_df = pd.DataFrame([
                    {'Sample Type': k, 'Count': v,
                     'Category': classify_sample_type(k) if 'classify_sample_type' in dir() else ''}
                    for k, v in sorted(sample_type_counts.items(), key=lambda x: -x[1])
                ])
                st.dataframe(stc_df, use_container_width=True, hide_index=True)
                
                # Show category totals
                from raptor.external_modules.acquisition.tcga import classify_sample_type as _cst
                tumor_total = sum(v for k, v in sample_type_counts.items() if _cst(k) == 'Tumor')
                normal_total = sum(v for k, v in sample_type_counts.items() if _cst(k) == 'Normal')
                normal_types = [f"{k} ({v})" for k, v in sample_type_counts.items() if _cst(k) == 'Normal']
                
                st.caption(
                    f"**Total tumor samples: {tumor_total}** | "
                    f"**Total normal samples: {normal_total}** "
                    + (f"({' + '.join(normal_types)})" if len(normal_types) > 1 else "")
                )
        
        # ---- TCGA-specific: paired tumor/normal samples ----
        samples_df = info.get('samples')
        if samples_df is not None and isinstance(samples_df, pd.DataFrame) and not samples_df.empty:
            if 'case_id' in samples_df.columns and 'sample_type' in samples_df.columns:
                from raptor.external_modules.acquisition.tcga import classify_sample_type
                samples_df['_category'] = samples_df['sample_type'].apply(classify_sample_type)
                
                # Find cases with both tumor and normal
                cases_with_tumor = set(samples_df.loc[samples_df['_category'] == 'Tumor', 'case_id'])
                cases_with_normal = set(samples_df.loc[samples_df['_category'] == 'Normal', 'case_id'])
                paired_cases = cases_with_tumor & cases_with_normal
                tumor_only = cases_with_tumor - cases_with_normal
                normal_only = cases_with_normal - cases_with_tumor
                
                with st.expander(f"Paired Tumor/Normal ({len(paired_cases)} patients)"):
                    st.markdown(
                        f"**{len(paired_cases)}** patients have both tumor and normal samples "
                        f"(ideal for differential expression)  \n"
                        f"**{len(tumor_only)}** patients have tumor only  \n"
                        f"**{len(normal_only)}** patients have normal only"
                    )
                    
                    if paired_cases:
                        # Build pairs table
                        pair_rows = []
                        for case_id in sorted(paired_cases):
                            case_samples = samples_df[samples_df['case_id'] == case_id]
                            tumor_samples = case_samples[case_samples['_category'] == 'Tumor']
                            normal_samples = case_samples[case_samples['_category'] == 'Normal']
                            
                            t_ids = ', '.join(tumor_samples['sample_id'].tolist()) if 'sample_id' in tumor_samples.columns else ''
                            t_types = ', '.join(tumor_samples['sample_type'].unique())
                            n_ids = ', '.join(normal_samples['sample_id'].tolist()) if 'sample_id' in normal_samples.columns else ''
                            n_types = ', '.join(normal_samples['sample_type'].unique())
                            
                            pair_rows.append({
                                'Patient': case_id,
                                'Tumor Sample(s)': t_ids,
                                'Tumor Type': t_types,
                                'Normal Sample(s)': n_ids,
                                'Normal Type': n_types,
                            })
                        
                        pairs_df = pd.DataFrame(pair_rows)
                        st.dataframe(pairs_df, use_container_width=True, hide_index=True)
                        
                        st.caption(
                            "Paired samples from the same patient let you compare tumor vs. normal "
                            "tissue directly — somatic mutations and differentially expressed genes "
                            "found this way are specific to the cancer, not inherited variation."
                        )
                    else:
                        st.caption(
                            "No matched pairs found. You can still do tumor vs. normal analysis "
                            "using unpaired samples, but matched pairs give stronger results."
                        )
                
                # Clean up temp column
                samples_df = samples_df.drop(columns=['_category'], errors='ignore')
        
        # ---- TCGA-specific: clinical summary ----
        clinical_summary = info.get('clinical_summary', {})
        if clinical_summary:
            with st.expander("Clinical Summary"):
                cs_cols = st.columns(min(len(clinical_summary), 4))
                col_idx = 0
                for key, value in clinical_summary.items():
                    with cs_cols[col_idx % len(cs_cols)]:
                        if isinstance(value, dict) and 'mean' in value:
                            # Age statistics
                            st.markdown(f"**{key.replace('_', ' ').title()}**")
                            st.caption(
                                f"Mean: {value['mean']} | Median: {value['median']} | "
                                f"Range: {value['min']}–{value['max']}"
                            )
                        elif isinstance(value, dict):
                            st.markdown(f"**{key.replace('_', ' ').title()}**")
                            for k, v in sorted(value.items(), key=lambda x: -x[1]):
                                st.caption(f"{k}: {v}")
                    col_idx += 1
        
        # ---- TCGA-specific: available data categories ----
        data_categories = info.get('data_categories', {})
        if data_categories:
            with st.expander("Available Data Types"):
                dc_df = pd.DataFrame([
                    {'Data Category': k, 'Files': v}
                    for k, v in sorted(data_categories.items(), key=lambda x: -x[1])
                ])
                st.dataframe(dc_df, use_container_width=True, hide_index=True)
        
        # Summary (GEO)
        summary = info.get('summary', '')
        if summary:
            with st.expander("Study Summary", expanded=True):
                st.write(summary)
        
        # Overall design
        design = info.get('overall_design', '')
        if design:
            with st.expander("Experimental Design"):
                st.write(design)
        
        # Sample table (grouped by case for TCGA)
        samples_df = info.get('samples')
        if samples_df is not None and isinstance(samples_df, pd.DataFrame) and not samples_df.empty:
            n_cases = samples_df['case_id'].nunique() if 'case_id' in samples_df.columns else len(samples_df)
            n_samples = len(samples_df)
            
            # Check if this is TCGA-style data with case_id grouping
            is_tcga_style = 'case_id' in samples_df.columns and 'sample_type' in samples_df.columns
            
            if is_tcga_style:
                from raptor.external_modules.acquisition.tcga import classify_sample_type as _cls
                
                st.info(
                    f"**{n_cases} patients** with **{n_samples} total biospecimens** across all data types "
                    f"(WXS, RNA-seq, methylation, miRNA, etc.). "
                    f"Not all biospecimens have RNA-seq data — only samples with STAR count files will be downloaded."
                )
                
                # ---- Case Summary Table (one row per patient) ----
                with st.expander(f"Case Summary ({n_cases} patients)", expanded=True):
                    case_rows = []
                    for case_id, group in samples_df.groupby('case_id'):
                        tumor_samples = group[group['sample_type'].apply(lambda x: _cls(x) == 'Tumor')]
                        normal_samples = group[group['sample_type'].apply(lambda x: _cls(x) == 'Normal')]
                        
                        tumor_ids = ', '.join(tumor_samples['sample_id'].tolist()) if 'sample_id' in group.columns else ''
                        normal_ids = ', '.join(normal_samples['sample_id'].tolist()) if 'sample_id' in group.columns else ''
                        tumor_types = ', '.join(tumor_samples['sample_type'].unique()) if len(tumor_samples) > 0 else ''
                        normal_types = ', '.join(normal_samples['sample_type'].unique()) if len(normal_samples) > 0 else ''
                        
                        # Get patient-level clinical info (same across samples)
                        first = group.iloc[0]
                        has_pair = len(tumor_samples) > 0 and len(normal_samples) > 0
                        
                        case_rows.append({
                            'Case': case_id,
                            'Paired': 'Yes' if has_pair else '',
                            'Tumor': f"{len(tumor_samples)} ({tumor_types})" if len(tumor_samples) > 0 else '',
                            'Normal': f"{len(normal_samples)} ({normal_types})" if len(normal_samples) > 0 else '',
                            'Gender': first.get('gender', 'unknown'),
                            'Vital Status': first.get('vital_status', 'unknown'),
                            'Age': first.get('age_years', ''),
                            'Diagnosis': first.get('primary_diagnosis', 'unknown'),
                            'Stage': first.get('tumor_stage', ''),
                            'Days to Death': first.get('days_to_death', ''),
                        })
                    
                    case_summary = pd.DataFrame(case_rows)
                    
                    # Color paired cases
                    def _color_paired(row):
                        if row.get('Paired') == 'Yes':
                            return ['background-color: #E3F2FD'] * len(row)
                        return [''] * len(row)
                    
                    styled = case_summary.style.apply(_color_paired, axis=1)
                    st.dataframe(styled, use_container_width=True, hide_index=True)
                    
                    n_paired = sum(1 for r in case_rows if r['Paired'] == 'Yes')
                    st.caption(
                        f"Blue rows = paired cases (tumor + normal from same patient, {n_paired} patients). "
                        f"'Tumor 2' means 2 tumor samples (e.g., different vials or replicate extractions)."
                    )
                
                # ---- Data Availability per case ----
                with st.expander("Data Availability per Case (RNA-Seq, miRNA, Methylation, etc.)"):
                    st.caption(
                        "Shows which experimental strategies are available for each patient. "
                        "Only open-access files are counted."
                    )
                    
                    preview_acc = direct_accession.strip() if direct_accession else ''
                    if not preview_acc:
                        preview_acc = st.session_state.get('preview_accession', '')
                    avail_key = f'data_avail_{preview_acc}'
                    
                    if st.button("Load Data Availability", key='load_avail_btn'):
                        with st.spinner("Querying GDC for per-case data types..."):
                            try:
                                connector = _get_connector('TCGA')
                                avail_df = connector.get_sample_data_availability(preview_acc)
                                st.session_state[avail_key] = avail_df
                            except Exception as e:
                                st.error(f"Failed: {e}")
                    
                    if avail_key in st.session_state:
                        avail_df = st.session_state[avail_key]
                        if not avail_df.empty:
                            meta_cols = {'case_id', 'n_samples', 'tumor_samples', 'normal_samples'}
                            strat_cols = [c for c in avail_df.columns if c not in meta_cols]
                            
                            # Summary counts
                            strat_counts = {}
                            for col in strat_cols:
                                n = (avail_df[col] == 'Yes').sum()
                                if n > 0:
                                    strat_counts[col] = n
                            
                            if strat_counts:
                                tags = [f"**{k}**: {v}/{len(avail_df)} cases" for k, v in sorted(strat_counts.items(), key=lambda x: -x[1])]
                                st.markdown(" · ".join(tags))
                            
                            # Display with checkmarks
                            display_avail = avail_df.copy()
                            for col in strat_cols:
                                display_avail[col] = display_avail[col].replace({'Yes': '✓', '': '—'})
                            
                            st.dataframe(display_avail, use_container_width=True, hide_index=True)
                            
                            st.caption(
                                "✓ = open-access data exists. RAPTOR downloads RNA-Seq by default. "
                                "To download other types, use **TCGA Download Options** above "
                                "(select miRNA, Methylation, CNV, or Mutations) or "
                                "**TCGA Advanced Tools** below."
                            )
                
                # ---- Detailed per-case view ----
                with st.expander("Detailed Sample View (per case)"):
                    case_select = st.selectbox(
                        "Select a case to inspect",
                        sorted(samples_df['case_id'].unique()),
                        key='case_detail_select',
                    )
                    case_group = samples_df[samples_df['case_id'] == case_select]
                    
                    # Patient-level info
                    first = case_group.iloc[0]
                    st.markdown(
                        f"**{case_select}** — "
                        f"{first.get('gender', '?')}, "
                        f"age {first.get('age_years', '?')}, "
                        f"{first.get('vital_status', '?')}"
                        + (f", died day {int(first['days_to_death'])}" if pd.notna(first.get('days_to_death')) else "")
                    )
                    
                    # Samples for this case
                    show_cols = [c for c in ['sample_id', 'sample_type', 'tissue_type']
                                 if c in case_group.columns]
                    
                    def _color_st(row):
                        st_val = str(row.get('sample_type', '')).lower()
                        if 'tumor' in st_val or 'cancer' in st_val:
                            return ['background-color: #FFEBEE'] * len(row)
                        elif 'normal' in st_val:
                            return ['background-color: #E8F5E9'] * len(row)
                        return [''] * len(row)
                    
                    st.dataframe(
                        case_group[show_cols].style.apply(_color_st, axis=1),
                        use_container_width=True, hide_index=True,
                    )
                    
                    # Decode barcodes
                    from raptor.external_modules.acquisition.tcga import parse_tcga_barcode
                    for _, row in case_group.iterrows():
                        sid = row.get('sample_id', '')
                        if sid.startswith('TCGA-'):
                            bc = parse_tcga_barcode(sid)
                            st.caption(
                                f"`{sid}` → type **{bc.get('sample_type_code', '?')}** "
                                f"({bc.get('sample_type', '?')}), "
                                f"vial {bc.get('vial', '?')}"
                            )
            
            else:
                # Non-TCGA: plain table
                with st.expander(f"Sample Details ({n_samples} samples)", expanded=True):
                    st.dataframe(samples_df, use_container_width=True, hide_index=True)
            
            # ---- Download customization (TCGA) ----
            if is_tcga_style:
                samples_df['_cat'] = samples_df['sample_type'].apply(_cls)
                
                cases_tumor = set(samples_df.loc[samples_df['_cat'] == 'Tumor', 'case_id'])
                cases_normal = set(samples_df.loc[samples_df['_cat'] == 'Normal', 'case_id'])
                paired = cases_tumor & cases_normal
                
                with st.expander("Download Options"):
                    st.markdown("**Choose which samples to download:**")
                    st.caption(
                        "Note: only samples with RNA-seq STAR count files will be downloaded. "
                        "The actual number may be less than shown here."
                    )
                    dl_option = st.radio(
                        "Sample selection",
                        [
                            f"All RNA-seq samples (default)",
                            f"Tumor only ({len(samples_df[samples_df['_cat'] == 'Tumor'])} biospecimens)",
                            f"Normal only ({len(samples_df[samples_df['_cat'] == 'Normal'])} biospecimens)",
                            f"Paired cases only ({len(paired)} cases with both tumor + normal)",
                            "Custom selection",
                        ],
                        key='tcga_dl_option',
                    )
                    
                    if dl_option.startswith("Custom"):
                        avail_types = samples_df['sample_type'].unique().tolist()
                        custom_types = st.multiselect(
                            "Select sample types",
                            avail_types,
                            default=avail_types,
                            key='tcga_custom_types',
                        )
                    
                    st.markdown("---")
                    st.caption(
                        "**About raw data:** TCGA provides processed gene expression counts (open access) "
                        "and aligned BAM files (controlled access — requires dbGaP approval via "
                        "[GDC portal](https://portal.gdc.cancer.gov)). "
                        "Raw FASTQ files are not directly available. "
                        "Use **Generate Manifest** in TCGA Advanced Tools to get a manifest for "
                        "the GDC Data Transfer Tool (`gdc-client`) if you need BAM files."
                    )
                
                samples_df = samples_df.drop(columns=['_cat'], errors='ignore')
        
        # Supplementary files
        supp_files = info.get('supplementary_files', [])
        if supp_files:
            with st.expander(f"Supplementary Files ({len(supp_files)})"):
                for f in supp_files:
                    fname = str(f).split('/')[-1] if '/' in str(f) else str(f)
                    st.text(fname)
        
        # PubMed
        pubmed = info.get('pubmed_ids', [])
        if pubmed and pubmed != ['']:
            pids = [str(p) for p in pubmed if p]
            if pids:
                st.caption("PubMed: " + ", ".join(
                    f"[{pid}](https://pubmed.ncbi.nlm.nih.gov/{pid})" for pid in pids
                ))
        
        # Dates
        sub_date = info.get('submission_date', '')
        upd_date = info.get('last_update', '')
        if sub_date or upd_date:
            st.caption(f"Submitted: {sub_date} | Updated: {upd_date}")
        
        st.markdown("")
    
    # ---- DOWNLOAD ----
    if direct_download and direct_accession:
        accession = direct_accession.strip()
        repo = _detect_repo(accession)
        connector = _get_connector(repo)
        
        with st.spinner(f"Downloading {accession} from {repo}... This may take a few minutes."):
            try:
                ds = connector.download(accession)
                st.session_state['acq_datasets'][accession] = ds
                
                catalog.register_dataset(
                    repository=ds.repository,
                    accession=ds.accession,
                    organism=ds.organism,
                    n_genes=ds.n_genes,
                    n_samples=ds.n_samples,
                    data_type=ds.data_type,
                )
                
                st.success(
                    f"Downloaded {accession}: "
                    f"{ds.n_genes:,} {_feature_label(ds)}, {ds.n_samples} samples "
                    f"({ds.organism})"
                )
                
                with st.expander("Preview downloaded data"):
                    st.dataframe(ds.counts_df.head(10), use_container_width=True)
                    st.caption(f"Gene IDs: {ds.gene_id_type} | Data type: {ds.data_type}")
                    
            except Exception as e:
                st.error(f"Download failed: {e}")
    
    st.markdown("---")
    st.markdown("## Search by Keywords")
    st.caption("Search across repositories using keywords")
    
    col_repo, col_query = st.columns([1, 3])
    
    with col_repo:
        repository = st.selectbox(
            "Repository",
            ["GEO", "TCGA", "ArrayExpress", "SRA"],
            help="Select which repository to search"
        )
    
    with col_query:
        search_query = st.text_input(
            "Search query",
            placeholder="e.g., pancreatic cancer, PS19 tau, breast cancer BRCA1",
            help="Keywords to search for datasets"
        )
    
    col_org, col_type, col_max, col_search_btn = st.columns([2, 1.5, 1, 1])
    
    with col_org:
        organism = st.text_input(
            "Organism filter (optional)",
            placeholder="e.g., Homo sapiens",
            help="Filter results by organism"
        )
    
    with col_type:
        # Data type options depend on repository
        if repository == "GEO":
            type_options = ["RNA-seq", "Microarray", "All expression", "Any"]
        elif repository == "SRA":
            type_options = ["RNA-seq"]  # SRA is sequencing only
        elif repository == "TCGA":
            type_options = [
                "RNA-seq", "miRNA-Seq", "Methylation Array",
                "Copy Number Variation", "Somatic Mutations", "Protein Expression (RPPA)", "All"
            ]
        else:  # ArrayExpress
            type_options = ["RNA-seq", "Microarray", "Any"]
        
        dataset_type = st.selectbox(
            "Data type",
            type_options,
            help="Filter by experiment type" if len(type_options) > 1
                else "SRA contains sequencing data only"
        )
    
    with col_max:
        max_results = st.number_input("Max results", min_value=5, max_value=200, value=20)
    
    with col_search_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        search_clicked = st.button("Search", type="primary", use_container_width=True)
    
    # SRA-specific filter
    if repository == "SRA":
        col_layout, _ = st.columns([1, 3])
        with col_layout:
            sra_layout = st.selectbox("Library layout", ["Any", "PAIRED", "SINGLE"],
                help="Filter by paired-end or single-end reads")
    else:
        sra_layout = "Any"
    
    if search_clicked and search_query:
        with st.spinner(f"Searching {repository}..."):
            try:
                # Create connector
                if repository == "GEO":
                    connector = GEOConnector(cache=True, cache_dir=cache.cache_dir)
                elif repository == "TCGA":
                    connector = TCGAConnector(cache=True, cache_dir=cache.cache_dir)
                elif repository == "ArrayExpress":
                    connector = ArrayExpConnector(cache=True, cache_dir=cache.cache_dir)
                else:
                    connector = SRAConnector(cache=True, cache_dir=cache.cache_dir)
                
                results = connector.search(
                    search_query,
                    organism=organism if organism else None,
                    max_results=max_results,
                    dataset_type=dataset_type,
                    library_layout=sra_layout if repository == "SRA" else None,
                )
                
                st.session_state['search_results'] = results
                st.session_state['search_connector'] = connector
                st.session_state['search_repo'] = repository
                
            except Exception as e:
                st.error(f"Search failed: {e}")
                results = []
    
    # Display search results
    if 'search_results' in st.session_state and st.session_state['search_results']:
        results = st.session_state['search_results']
        search_repo = st.session_state.get('search_repo', 'GEO')
        
        st.success(f"Found {len(results)} datasets")
        
        # Build rich results table (columns depend on repository)
        table_rows = []
        for r in results:
            if search_repo == 'SRA':
                gse = r.extra.get('gse', '')
                geo_hint = r.extra.get('geo_hint', '')
                # Show GSE if found, or hint that GEO link exists
                geo_display = gse if gse else geo_hint
                table_rows.append({
                    'Accession': r.accession,
                    'GEO Link': geo_display,
                    'Title': r.title[:70] + ('...' if len(r.title) > 70 else ''),
                    'Organism': r.organism,
                    'Samples': r.n_samples,
                    'Runs': r.extra.get('n_runs', ''),
                    'Layout': r.extra.get('library_layout', ''),
                    'Avg Reads': r.extra.get('avg_reads', ''),
                    'Platform': r.platform,
                })
            else:
                pubmed = r.extra.get('pubmed_ids', [])
                gpl_id = r.extra.get('gpl_id', '')
                # TCGA shows cases (patients), GEO shows samples
                count_label = 'Cases' if search_repo == 'TCGA' else 'Samples'
                table_rows.append({
                    'Accession': r.accession,
                    'Title': r.title[:70] + ('...' if len(r.title) > 70 else ''),
                    'Type': r.extra.get('gds_type', ''),
                    'Organism': r.organism,
                    count_label: r.n_samples,
                    'Platform': r.platform,
                    'GPL': f"GPL{gpl_id}" if gpl_id and not gpl_id.startswith('GPL') else gpl_id,
                    'PubMed': pubmed[0] if pubmed else '',
                    'Cached': '✅' if cache.is_cached(search_repo, r.accession) else '',
                })
        
        results_df = pd.DataFrame(table_rows)
        # Drop columns that have no useful data (all empty, all zero, or all NaN)
        for col in list(results_df.columns):
            vals = results_df[col]
            if vals.dtype in ('int64', 'float64'):
                if (vals == 0).all() or vals.isna().all():
                    results_df = results_df.drop(columns=[col])
            else:
                if vals.astype(str).isin(['', 'nan', 'None']).all():
                    results_df = results_df.drop(columns=[col])

        st.dataframe(results_df, use_container_width=True, hide_index=True)
        
        # Dataset details on selection
        st.markdown("### Dataset Actions")
        
        accession_options = [r.accession for r in results]
        selected_accession = st.selectbox("Select dataset", accession_options)
        
        # Show description of selected dataset
        selected_result = next((r for r in results if r.accession == selected_accession), None)
        if selected_result:
            if search_repo == 'SRA':
                n_runs = selected_result.extra.get('n_runs', 0)
                layout = selected_result.extra.get('library_layout', '')
                avg_reads = selected_result.extra.get('avg_reads', '')
                sra_gse = selected_result.extra.get('gse', '')
                has_geo = selected_result.extra.get('has_geo_link', False)
                info_parts = [
                    f"**{n_runs} run(s)**",
                    f"Layout: {layout}",
                    f"Avg reads/run: {avg_reads}",
                    f"Platform: {selected_result.platform}",
                ]
                st.info(" | ".join(info_parts))
                if sra_gse:
                    st.caption(
                        f"Linked GEO study: "
                        f"[{sra_gse}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sra_gse}) "
                        f"— download processed counts from GEO instead of raw FASTQs"
                    )
                elif has_geo:
                    st.caption(
                        "This study has GEO sample accessions (GSM IDs) — "
                        "the linked GSE will be resolved automatically when you download the run table."
                    )
            else:
                design = selected_result.extra.get('design_hint', '')
                if design:
                    st.info(design)
                elif selected_result.description:
                    st.caption(selected_result.description[:300])
                
                pubmed = selected_result.extra.get('pubmed_ids', [])
                if pubmed:
                    links = ", ".join(
                        f"[{pid}](https://pubmed.ncbi.nlm.nih.gov/{pid})" for pid in pubmed
                    )
                    st.caption(f"PubMed: {links}")
        
        # SRA: Show Run Table / Others: Show Sample Types
        if search_repo == 'SRA':
            if st.button("Show Run Table", key='sra_run_table_btn'):
                connector = st.session_state.get('search_connector')
                if connector:
                    with st.spinner(f"Fetching run table for {selected_accession}..."):
                        try:
                            run_table = connector.get_run_table(selected_accession)
                            st.session_state['sra_run_table'] = run_table
                            st.session_state['sra_run_table_acc'] = selected_accession
                        except Exception as e:
                            st.error(f"Failed: {e}")
            if (st.session_state.get('sra_run_table_acc') == selected_accession
                    and 'sra_run_table' in st.session_state):
                rt = st.session_state['sra_run_table']
                st.markdown(f"**Run Table** ({len(rt)} runs)")
                show_cols = [c for c in ['run_accession', 'sample_title', 'sample_alias',
                    'instrument_model', 'library_layout', 'library_selection',
                    'read_count', 'base_count', 'scientific_name'] if c in rt.columns]
                st.dataframe(rt[show_cols] if show_cols else rt, use_container_width=True, hide_index=True)
                st.download_button("Download Run Table (CSV)", data=rt.to_csv(index=False),
                    file_name=f"{selected_accession}_run_table.csv", mime="text/csv")
        else:
            # GEO/TCGA/ArrayExpress: Show Sample Types
            if st.button("Show Sample Types", use_container_width=False, key='sample_types_btn'):
                connector = st.session_state.get('search_connector')
                if connector and hasattr(connector, 'get_sample_types'):
                    with st.spinner(f"Fetching sample info for {selected_accession}..."):
                        try:
                            stypes = connector.get_sample_types(selected_accession)
                            st.session_state['sample_types_result'] = stypes
                            st.session_state['sample_types_accession'] = selected_accession
                        except Exception as e:
                            st.error(f"Failed: {e}")
            if (st.session_state.get('sample_types_accession') == selected_accession
                    and 'sample_types_result' in st.session_state):
                stypes = st.session_state['sample_types_result']
                type_counts = stypes.get('type_counts', {})
                if type_counts:
                    tags = [f"**{name}** ({count})" for name, count in type_counts.items()]
                    st.markdown("**Sample types:** " + " · ".join(tags))
                    # TCGA-specific: tumor/normal counts
                    tumor_n = stypes.get('tumor_count', 0)
                    normal_n = stypes.get('normal_count', 0)
                    if tumor_n or normal_n:
                        st.caption(f"Tumor: {tumor_n} | Normal: {normal_n}")
        
        # ---- TCGA download options ----
        if search_repo == 'TCGA':
            with st.expander("TCGA Download Options"):
                tcga_opt_col1, tcga_opt_col2 = st.columns(2)
                with tcga_opt_col1:
                    # Sample type filter
                    stypes_result = st.session_state.get('sample_types_result', {})
                    avail_types = stypes_result.get('unique_types', [])
                    if not avail_types:
                        avail_types = ['Primary Tumor', 'Solid Tissue Normal', 'Blood Derived Normal', 'Metastatic']
                    tcga_sample_filter = st.multiselect(
                        "Filter by sample type (leave empty for all)",
                        avail_types,
                        key='tcga_sample_filter',
                        help="Select which sample types to download"
                    )
                with tcga_opt_col2:
                    from raptor.external_modules.acquisition.tcga import QUANT_COLUMNS
                    tcga_quant = st.selectbox(
                        "Quantification",
                        list(QUANT_COLUMNS.keys()),
                        format_func=lambda x: QUANT_COLUMNS[x],
                        key='tcga_quant_select',
                        help="Which column to extract from STAR-Counts files"
                    )
                
                tcga_data_type = st.selectbox(
                    "Data type to download",
                    ["Gene Expression", "miRNA Expression", "Copy Number (gene-level)",
                     "Copy Number (segments)", "DNA Methylation", "Somatic Mutations",
                     "Protein Expression (RPPA)"],
                    key='tcga_data_type_select',
                    help="TCGA has multiple data types beyond gene expression"
                )
        
        col_preview, col_dl, col_info = st.columns(3)
        
        with col_preview:
            if search_repo != 'SRA':
                if st.button("Preview Info", use_container_width=True, key='search_preview'):
                    connector = st.session_state.get('search_connector')
                    if connector and hasattr(connector, 'get_study_info'):
                        with st.spinner(f"Fetching info for {selected_accession}..."):
                            try:
                                info = connector.get_study_info(selected_accession)
                                st.session_state['preview_info'] = info
                                st.session_state['preview_accession'] = selected_accession
                            except Exception as e:
                                st.error(f"Failed: {e}")
        
        with col_dl:
            dl_label = "Get Run Table" if search_repo == 'SRA' else "Download Dataset"
            if st.button(dl_label, type="primary", use_container_width=True, key='search_dl'):
                connector = st.session_state.get('search_connector')
                if connector:
                    # Progress bar for TCGA (large downloads)
                    if search_repo == 'TCGA':
                        prog_col, cancel_col = st.columns([4, 1])
                        with prog_col:
                            progress_bar = st.progress(0, text=f"Preparing download for {selected_accession}...")
                        with cancel_col:
                            cancel_placeholder = st.empty()
                        
                        # Cancel file flag — progress callback checks this
                        cancel_flag_path = Path.home() / '.raptor' / '_cancel_download'
                        cancel_flag_path.unlink(missing_ok=True)  # clear any stale flag
                        
                        # Show cancel button (clicking it reruns page, stopping download)
                        if cancel_placeholder.button("Cancel", key='cancel_dl', type="secondary"):
                            cancel_flag_path.touch()
                            st.warning("Download cancelled. Cached data on disk is preserved.")
                            st.stop()

                        def tcga_progress(downloaded, total):
                            # Check cancel flag
                            if cancel_flag_path.exists():
                                cancel_flag_path.unlink(missing_ok=True)
                                raise Exception("Download cancelled by user")
                            pct = min(downloaded / max(total, 1), 1.0)
                            progress_bar.progress(pct, text=f"Downloading {selected_accession}: {downloaded}/{total} files ({pct*100:.0f}%)")

                        try:
                            sample_filter = st.session_state.get('tcga_sample_filter', [])
                            quant = st.session_state.get('tcga_quant_select', 'unstranded')
                            data_type_sel = st.session_state.get('tcga_data_type_select', 'Gene Expression')

                            has_filters = bool(sample_filter) or quant != 'unstranded'

                            if sample_filter:
                                filter_tag = '_'.join(s.replace(' ', '') for s in sorted(sample_filter))
                                ds_key = f"{selected_accession}_{filter_tag}"
                            else:
                                ds_key = selected_accession

                            if data_type_sel == 'Gene Expression':
                                ds = connector.download(
                                    selected_accession,
                                    sample_types=sample_filter if sample_filter else None,
                                    quant_column=quant,
                                    force=has_filters,
                                    progress_callback=tcga_progress,
                                )
                            elif data_type_sel == 'miRNA Expression':
                                ds = connector.download_mirna(
                                    selected_accession,
                                    sample_types=sample_filter if sample_filter else None,
                                    progress_callback=tcga_progress,
                                )
                                ds_key = f"{selected_accession}_miRNA"
                            elif data_type_sel.startswith('Copy Number (gene'):
                                ds = connector.download_cnv(
                                    selected_accession, cnv_type='gene_level',
                                    sample_types=sample_filter if sample_filter else None,
                                    progress_callback=tcga_progress,
                                )
                                ds_key = f"{selected_accession}_CNV"
                            elif data_type_sel.startswith('Copy Number (seg'):
                                ds = connector.download_cnv(
                                    selected_accession, cnv_type='masked',
                                    sample_types=sample_filter if sample_filter else None,
                                    progress_callback=tcga_progress,
                                )
                                ds_key = f"{selected_accession}_CNV_seg"
                            elif data_type_sel == 'DNA Methylation':
                                ds = connector.download_methylation(
                                    selected_accession,
                                    sample_types=sample_filter if sample_filter else None,
                                    progress_callback=tcga_progress,
                                )
                                ds_key = f"{selected_accession}_methylation"
                            elif data_type_sel == 'Somatic Mutations':
                                maf_df = connector.download_mutations(
                                    selected_accession,
                                    progress_callback=tcga_progress,
                                )
                                # Build metadata from MAF barcodes
                                maf_meta = pd.DataFrame()
                                if 'Tumor_Sample_Barcode' in maf_df.columns:
                                    barcodes = maf_df['Tumor_Sample_Barcode'].unique()
                                    maf_meta = pd.DataFrame({
                                        'sample_id': barcodes,
                                        'case_id': [b[:12] if len(b) >= 12 else b for b in barcodes],
                                        'mutations': [
                                            (maf_df['Tumor_Sample_Barcode'] == b).sum() for b in barcodes
                                        ],
                                    })
                                    maf_meta = maf_meta.set_index('sample_id')
                                
                                ds = AcquiredDataset(
                                    counts_df=maf_df,
                                    metadata=maf_meta,
                                    source_info={
                                        'repository': 'TCGA',
                                        'accession': selected_accession,
                                        'data_type': 'somatic_mutation',
                                        'organism': 'Homo sapiens',
                                    },
                                    gene_id_type='symbol',
                                )
                                ds_key = f"{selected_accession}_mutations"
                            elif data_type_sel == 'Protein Expression (RPPA)':
                                ds = connector.download_rppa(
                                    selected_accession,
                                    sample_types=sample_filter if sample_filter else None,
                                    progress_callback=tcga_progress,
                                )
                                ds_key = f"{selected_accession}_RPPA"
                            else:
                                ds = connector.download(selected_accession, progress_callback=tcga_progress)
                                ds_key = selected_accession

                            st.session_state['acq_datasets'][ds_key] = ds
                            progress_bar.progress(1.0, text="Download complete!")

                            catalog.register_dataset(
                                repository=ds.repository, accession=ds.accession,
                                organism=ds.organism, n_genes=ds.n_genes,
                                n_samples=ds.n_samples, data_type=ds.data_type,
                            )
                            filter_desc = f" (filter: {', '.join(sample_filter)})" if sample_filter else ""
                            st.success(f"Downloaded {ds_key}: {ds.n_genes:,} {_feature_label(ds)}, {ds.n_samples} cases{filter_desc}")

                        except Exception as e:
                            progress_bar.empty()
                            err_msg = str(e)
                            st.error(f"Download failed: Failed to download {selected_accession} from TCGA: {err_msg}")
                            logger.error(f"TCGA download error: {e}", exc_info=True)
                            if "end-of-stream" in err_msg.lower() or "compressed" in err_msg.lower():
                                st.info(
                                    "This usually means the download was cut off mid-transfer "
                                    "(the file was incomplete when decompression started). "
                                    "Try again — GDC servers can be intermittent for large projects. "
                                    "If it keeps failing, try filtering to fewer sample types to reduce download size."
                                )

                    else:
                        # Non-TCGA download (GEO, ArrayExpress, SRA)
                        with st.spinner(f"{'Fetching' if search_repo == 'SRA' else 'Downloading'} {selected_accession}..."):
                            try:
                                ds = connector.download(selected_accession)

                                # Carry over GSE link from search result into source_info
                                if search_repo == 'SRA' and selected_result:
                                    gse_from_search = selected_result.extra.get('gse', '')
                                    if gse_from_search and hasattr(ds, 'source_info'):
                                        ds.source_info['gse'] = gse_from_search

                                st.session_state['acq_datasets'][selected_accession] = ds

                                catalog.register_dataset(
                                    repository=ds.repository, accession=ds.accession,
                                    organism=ds.organism, n_genes=ds.n_genes,
                                    n_samples=ds.n_samples, data_type=ds.data_type,
                                )
                                if search_repo == 'SRA':
                                    st.success(f"Run table loaded: {selected_accession} ({ds.n_samples} runs)")
                                else:
                                    st.success(f"Downloaded {selected_accession}: {ds.n_genes:,} {_feature_label(ds)}, {ds.n_samples} samples")
                            except Exception as e:
                                st.error(f"Failed: {e}")
                                
                                # ArrayExpress: show helpful navigation links
                                if search_repo == 'ArrayExpress':
                                    links = []
                                    links.append(f"[Open {selected_accession} in BioStudies](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{selected_accession})")
                                    links.append(f"[Browse raw data on ENA](https://www.ebi.ac.uk/ena/browser/view/{selected_accession})")
                                    
                                    # E-GEOD accessions have a linked GEO study
                                    if selected_accession.startswith('E-GEOD-'):
                                        gse_id = 'GSE' + selected_accession.replace('E-GEOD-', '')
                                        links.append(f"[Download from GEO ({gse_id})](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}) — likely has processed counts")
                                    
                                    st.markdown(" · ".join(links))
        
        with col_info:
            if search_repo == 'SRA':
                sra_gse_link = selected_result.extra.get('gse', '') if selected_result else ''
                if sra_gse_link:
                    st.markdown(
                        f"[Open in SRA](https://www.ncbi.nlm.nih.gov/sra/?term={selected_accession}) "
                        f"| [Open in GEO ({sra_gse_link})](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sra_gse_link})"
                    )
                else:
                    st.markdown(f"[Open {selected_accession} in SRA](https://www.ncbi.nlm.nih.gov/sra/?term={selected_accession})")
            else:
                st.markdown(_repo_link_markdown(selected_accession, search_repo))
        
        # Display preview info from search results
        if (st.session_state.get('preview_accession') == selected_accession
                and 'preview_info' in st.session_state):
            info = st.session_state['preview_info']
            
            st.markdown("---")
            st.markdown(f"### {info.get('title', 'No title')}")
            st.markdown(_repo_link_markdown(selected_accession, search_repo))
            
            n_cases = info.get('n_samples', info.get('n_cases', 'N/A'))
            n_samp = len(info.get('samples', [])) if isinstance(info.get('samples'), pd.DataFrame) else ''
            if search_repo == 'TCGA':
                cases_str = f"**{n_cases}** patients"
                if n_samp and n_samp != n_cases:
                    cases_str += f" ({n_samp} biospecimens across all assays)"
            else:
                cases_str = f"**{n_cases}** cases"
                if n_samp and n_samp != n_cases:
                    cases_str += f" ({n_samp} samples)"
            st.markdown(
                f"{cases_str} · **{info.get('organism', 'unknown')}** · "
                f"{info.get('platform_id', 'unknown')} · {info.get('data_type', 'unknown')}"
            )
            
            pname = info.get('platform_name', '')
            if pname:
                st.caption(f"Platform: {pname}")
            
            summary = info.get('summary', '')
            if summary:
                with st.expander("Study Summary", expanded=True):
                    st.write(summary)
            
            samples_df = info.get('samples')
            if samples_df is not None and isinstance(samples_df, pd.DataFrame) and not samples_df.empty:
                is_tcga = 'case_id' in samples_df.columns and 'sample_type' in samples_df.columns
                n_cases_prev = samples_df['case_id'].nunique() if 'case_id' in samples_df.columns else len(samples_df)
                
                if is_tcga:
                    from raptor.external_modules.acquisition.tcga import classify_sample_type as _cls2
                    
                    # Clear info banner
                    st.info(
                        f"**TCGA project = multi-omic.** This project has **{n_cases_prev} patients** with "
                        f"**{len(samples_df)} total biospecimens** across all assays (RNA-seq, miRNA, methylation, WXS, etc.). "
                        f"When you click Download, only the data type you selected in **TCGA Download Options** "
                        f"(default: RNA-seq gene expression) will be downloaded."
                    )
                    
                    # Case summary table (most useful — expanded by default)
                    stc2 = info.get('sample_type_counts', {})
                    case_rows2 = []
                    for cid, grp in samples_df.groupby('case_id'):
                        t_samp = grp[grp['sample_type'].apply(lambda x: _cls2(x) == 'Tumor')]
                        n_samp2 = grp[grp['sample_type'].apply(lambda x: _cls2(x) == 'Normal')]
                        first2 = grp.iloc[0]
                        has_pair2 = len(t_samp) > 0 and len(n_samp2) > 0
                        case_rows2.append({
                            'Case': cid,
                            'Paired': 'Yes' if has_pair2 else '',
                            'Tumor': f"{len(t_samp)} ({', '.join(t_samp['sample_type'].unique())})" if len(t_samp) > 0 else '',
                            'Normal': f"{len(n_samp2)} ({', '.join(n_samp2['sample_type'].unique())})" if len(n_samp2) > 0 else '',
                            'Gender': first2.get('gender', ''),
                            'Age': first2.get('age_years', ''),
                            'Vital Status': first2.get('vital_status', ''),
                            'Days to Death': first2.get('days_to_death', ''),
                        })
                    
                    n_paired2 = sum(1 for r in case_rows2 if r['Paired'] == 'Yes')
                    with st.expander(f"Case Summary — {n_cases_prev} patients, {n_paired2} paired (tumor+normal)", expanded=True):
                        cs_df2 = pd.DataFrame(case_rows2)
                        def _cp2(row):
                            if row.get('Paired') == 'Yes':
                                return ['background-color: #E3F2FD'] * len(row)
                            return [''] * len(row)
                        st.dataframe(cs_df2.style.apply(_cp2, axis=1), use_container_width=True, hide_index=True)
                        st.caption("Blue rows = patients with both tumor and normal samples (ideal for DE analysis)")
                    
                    # Sample type breakdown
                    if stc2:
                        with st.expander("Biospecimen Breakdown (all assays, not just RNA-seq)"):
                            stc_df2 = pd.DataFrame([
                                {'Sample Type': k, 'Count': v, 'Category': _cls2(k)}
                                for k, v in sorted(stc2.items(), key=lambda x: -x[1])
                            ])
                            st.dataframe(stc_df2, use_container_width=True, hide_index=True)
                            tumor_t = sum(v for k, v in stc2.items() if _cls2(k) == 'Tumor')
                            normal_t = sum(v for k, v in stc2.items() if _cls2(k) == 'Normal')
                            normal_detail = [f"{k} ({v})" for k, v in stc2.items() if _cls2(k) == 'Normal']
                            st.caption(
                                f"**Tumor: {tumor_t}** | **Normal: {normal_t}** "
                                + (f"({' + '.join(normal_detail)})" if len(normal_detail) > 1 else "")
                            )
                    
                    # Data availability
                    with st.expander("Data Availability per Case (which assays exist for each patient)"):
                        avail_key2 = f'data_avail_search_{selected_accession}'
                        if st.button("Load Data Availability", key='load_avail_search_btn'):
                            with st.spinner("Querying GDC for per-case data types..."):
                                try:
                                    connector = st.session_state.get('search_connector')
                                    if connector and hasattr(connector, 'get_sample_data_availability'):
                                        avail_df2 = connector.get_sample_data_availability(selected_accession)
                                        st.session_state[avail_key2] = avail_df2
                                except Exception as e:
                                    st.error(f"Failed: {e}")
                        
                        if avail_key2 in st.session_state:
                            avail_df2 = st.session_state[avail_key2]
                            if not avail_df2.empty:
                                meta_cols2 = {'case_id', 'n_samples', 'tumor_samples', 'normal_samples'}
                                strat_cols2 = [c for c in avail_df2.columns if c not in meta_cols2]
                                counts2 = {col: (avail_df2[col] == 'Yes').sum() for col in strat_cols2}
                                counts2 = {k: v for k, v in counts2.items() if v > 0}
                                if counts2:
                                    tags2 = [f"**{k}**: {v}/{len(avail_df2)}" for k, v in sorted(counts2.items(), key=lambda x: -x[1])]
                                    st.markdown(" · ".join(tags2))
                                disp2 = avail_df2.copy()
                                for col in strat_cols2:
                                    disp2[col] = disp2[col].replace({'Yes': '✓', '': '—'})
                                st.dataframe(disp2, use_container_width=True, hide_index=True)
                                st.caption(
                                    "✓ = open-access data available. "
                                    "To download: expand **TCGA Download Options** above → "
                                    "select data type → click **Download Dataset**."
                                )
                    
                    # Clinical summary
                    cs_info = info.get('clinical_summary', {})
                    if cs_info:
                        with st.expander("Clinical Summary"):
                            cs_c = st.columns(min(len(cs_info), 4))
                            ci = 0
                            for key, val in cs_info.items():
                                with cs_c[ci % len(cs_c)]:
                                    if isinstance(val, dict) and 'mean' in val:
                                        st.markdown(f"**{key.replace('_', ' ').title()}**")
                                        st.caption(f"Mean: {val['mean']} | Median: {val['median']} | Range: {val['min']}–{val['max']}")
                                    elif isinstance(val, dict):
                                        st.markdown(f"**{key.replace('_', ' ').title()}**")
                                        for k2, v2 in sorted(val.items(), key=lambda x: -x[1]):
                                            st.caption(f"{k2}: {v2}")
                                ci += 1
                else:
                    with st.expander(f"Sample Details ({len(samples_df)} samples)"):
                        st.dataframe(samples_df, use_container_width=True, hide_index=True)
            
            supp_files = info.get('supplementary_files', [])
            if supp_files:
                with st.expander(f"Supplementary Files ({len(supp_files)})"):
                    for f in supp_files:
                        fname = str(f).split('/')[-1] if '/' in str(f) else str(f)
                        st.text(fname)
    
    # TCGA project browser
    if repository == "TCGA":
        with st.expander("Browse TCGA/TARGET Projects"):
            from raptor.external_modules.acquisition.tcga import TCGA_PROJECTS, TARGET_PROJECTS
            
            prog_filter = st.radio("Program", ["TCGA (33 cancers)", "TARGET (pediatric)", "All"],
                                   horizontal=True, key='tcga_prog_filter')
            if prog_filter.startswith("TCGA"):
                proj_dict = TCGA_PROJECTS
            elif prog_filter.startswith("TARGET"):
                proj_dict = TARGET_PROJECTS
            else:
                proj_dict = {**TCGA_PROJECTS, **TARGET_PROJECTS}
            
            proj_df = pd.DataFrame([
                {'Project': k, 'Description': v}
                for k, v in proj_dict.items()
            ])
            st.dataframe(proj_df, use_container_width=True, hide_index=True)
            
            st.caption(f"{len(proj_dict)} projects. Select one above and click Search, or enter the project ID directly.")
        
        # TCGA advanced tools
        with st.expander("TCGA Advanced Tools"):
            tcga_tool = st.selectbox(
                "Tool",
                ["Paired Tumor/Normal", "Survival Data", "Top Mutated Genes",
                 "Multi-Project Download", "Generate Manifest", "Compare Projects", "GDC Status"],
                key='tcga_adv_tool',
            )
            
            if tcga_tool == "GDC Status":
                if st.button("Check GDC API Status", key='gdc_status_btn'):
                    tcga = TCGAConnector(cache=False)
                    status = tcga.check_status()
                    st.json(status)
            
            elif tcga_tool == "Paired Tumor/Normal":
                pair_project = st.text_input("Project ID", value="TCGA-LIHC", key='pair_proj')
                if st.button("Find Paired Samples", key='pair_btn'):
                    with st.spinner("Querying GDC for tumor/normal pairs..."):
                        try:
                            tcga = TCGAConnector(cache=False)
                            pairs = tcga.get_paired_samples(pair_project)
                            if pairs.empty:
                                st.warning("No paired tumor/normal samples found.")
                            else:
                                st.success(
                                    f"Found {len(pairs)} pairs from "
                                    f"{pairs['patient_barcode'].nunique()} patients"
                                )
                                st.dataframe(pairs, use_container_width=True, hide_index=True)
                                csv_data = pairs.to_csv(index=False)
                                st.download_button(
                                    "Download Pairs (CSV)", data=csv_data,
                                    file_name=f"{pair_project}_paired_samples.csv",
                                    mime="text/csv",
                                )
                        except Exception as e:
                            st.error(f"Failed: {e}")
            
            elif tcga_tool == "Survival Data":
                surv_project = st.text_input("Project ID", value="TCGA-LIHC", key='surv_proj')
                if st.button("Get Survival Data", key='surv_btn'):
                    with st.spinner("Fetching survival data..."):
                        try:
                            tcga = TCGAConnector(cache=False)
                            surv = tcga.get_survival_data(surv_project)
                            n_events = (surv['overall_survival_event'] == 1).sum()
                            st.success(f"{len(surv)} patients, {n_events} events (deaths)")
                            st.dataframe(surv.head(50), use_container_width=True, hide_index=True)
                            csv_data = surv.to_csv(index=False)
                            st.download_button(
                                "Download Survival Data (CSV)", data=csv_data,
                                file_name=f"{surv_project}_survival.csv",
                                mime="text/csv",
                            )
                        except Exception as e:
                            st.error(f"Failed: {e}")
            
            elif tcga_tool == "Top Mutated Genes":
                mut_project = st.text_input("Project ID", value="TCGA-BRCA", key='mut_proj')
                n_top = st.slider("Number of genes", 10, 100, 30, key='mut_n')
                if st.button("Get Mutation Frequencies", key='mut_btn'):
                    with st.spinner("Querying GDC analysis endpoint..."):
                        try:
                            tcga = TCGAConnector(cache=False)
                            freq = tcga.get_gene_mutation_frequency(mut_project, n_genes=n_top)
                            if freq.empty:
                                st.warning("No mutation data available.")
                            else:
                                st.dataframe(freq, use_container_width=True, hide_index=True)
                                # Bar chart
                                if len(freq) > 0 and 'symbol' in freq.columns and 'num_cases' in freq.columns:
                                    fig = go.Figure(data=[
                                        go.Bar(x=freq['symbol'].head(20),
                                               y=freq['num_cases'].head(20),
                                               marker_color='#C62828')
                                    ])
                                    fig.update_layout(
                                        title=f"Top Mutated Genes — {mut_project}",
                                        xaxis_title="Gene", yaxis_title="Cases with mutation",
                                        height=400,
                                    )
                                    st.plotly_chart(fig, use_container_width=True)
                        except Exception as e:
                            st.error(f"Failed: {e}")
            
            elif tcga_tool == "Multi-Project Download":
                st.caption("Download and merge gene expression from multiple projects (e.g., all lung cancers)")
                multi_projects = st.text_input(
                    "Project IDs (comma-separated)",
                    value="TCGA-LUAD, TCGA-LUSC",
                    key='multi_proj_input',
                )
                multi_sample_filter = st.multiselect(
                    "Sample type filter (optional)",
                    ['Primary Tumor', 'Solid Tissue Normal', 'Blood Derived Normal', 'Metastatic'],
                    key='multi_sample_filter',
                )
                if st.button("Download & Merge", type="primary", key='multi_dl_btn'):
                    proj_list = [p.strip() for p in multi_projects.split(',') if p.strip()]
                    if len(proj_list) < 2:
                        st.warning("Enter at least 2 project IDs")
                    else:
                        mprog_col, mcancel_col = st.columns([4, 1])
                        with mprog_col:
                            multi_progress = st.progress(0, text=f"Starting download of {len(proj_list)} projects...")
                        with mcancel_col:
                            mcancel_placeholder = st.empty()
                        
                        cancel_flag_path = Path.home() / '.raptor' / '_cancel_download'
                        cancel_flag_path.unlink(missing_ok=True)
                        
                        if mcancel_placeholder.button("Cancel", key='cancel_multi_dl', type="secondary"):
                            cancel_flag_path.touch()
                            st.warning("Download cancelled. Cached data on disk is preserved.")
                            st.stop()
                        
                        try:
                            tcga = TCGAConnector(cache=True, cache_dir=cache.cache_dir)

                            n_projects = len(proj_list)
                            
                            def multi_callback(downloaded, total):
                                if cancel_flag_path.exists():
                                    cancel_flag_path.unlink(missing_ok=True)
                                    raise Exception("Download cancelled by user")
                                pct = min(downloaded / max(total, 1), 1.0)
                                if total <= n_projects:
                                    # Project-level progress (merge phase)
                                    multi_progress.progress(pct, text=f"Merging project {downloaded}/{total}...")
                                else:
                                    # File-level progress (download phase)
                                    multi_progress.progress(pct, text=f"Downloading: {downloaded}/{total} files ({pct*100:.0f}%)")

                            ds = tcga.download_multiple_projects(
                                proj_list,
                                sample_types=multi_sample_filter if multi_sample_filter else None,
                                progress_callback=multi_callback,
                            )
                            
                            # Get actual projects that succeeded (not requested)
                            actual_projects = []
                            study_labels_series = pd.Series('unknown', index=ds.metadata.index)
                            if 'project' in ds.metadata.columns:
                                study_labels_series = ds.metadata['project']
                                actual_projects = ds.metadata['project'].unique().tolist()
                            
                            pool_name = '+'.join(actual_projects) if actual_projects else '+'.join(proj_list)
                            
                            from raptor.external_modules.acquisition import PooledDataset
                            pooled_result = PooledDataset(
                                counts_df=ds.counts_df,
                                metadata=ds.metadata,
                                study_labels=study_labels_series,
                                pooling_info={'method': 'multi_project_download', 'projects': actual_projects},
                                component_datasets=actual_projects,
                            )
                            st.session_state['acq_pooled'] = pooled_result
                            multi_progress.progress(1.0, text="Download complete!")
                            
                            n_actual = len(actual_projects)
                            n_requested = len(proj_list)
                            if n_actual < n_requested:
                                failed = [p for p in proj_list if p not in actual_projects]
                                st.warning(
                                    f"Only {n_actual}/{n_requested} projects downloaded successfully. "
                                    f"Failed: {', '.join(failed)}. Try downloading them individually."
                                )
                            st.success(
                                f"Merged: {ds.n_genes:,} {_feature_label(ds)}, {ds.n_samples} samples "
                                f"from {n_actual} projects ({', '.join(actual_projects)})"
                            )
                        except Exception as e:
                            multi_progress.empty()
                            st.error(f"Failed: {e}")
            
            elif tcga_tool == "Generate Manifest":
                st.caption("Generate a manifest TSV for the GDC Data Transfer Tool (gdc-client)")
                man_project = st.text_input("Project ID", value="TCGA-LIHC", key='man_proj')
                if st.button("Generate Manifest", key='man_btn'):
                    with st.spinner("Generating manifest..."):
                        try:
                            tcga = TCGAConnector(cache=False)
                            manifest = tcga.generate_manifest(man_project)
                            st.success(f"Manifest: {len(manifest)} files")
                            st.dataframe(manifest.head(20), use_container_width=True, hide_index=True)
                            csv_data = manifest.to_csv(sep='\t', index=False)
                            st.download_button(
                                "Download Manifest (TSV)", data=csv_data,
                                file_name=f"{man_project}_manifest.tsv",
                                mime="text/tab-separated-values",
                            )
                            st.code(f"gdc-client download -m {man_project}_manifest.tsv", language="bash")
                        except Exception as e:
                            st.error(f"Failed: {e}")
            
            elif tcga_tool == "Compare Projects":
                comp_projects = st.text_input(
                    "Project IDs (comma-separated)",
                    value="TCGA-BRCA, TCGA-LUAD, TCGA-LIHC",
                    key='comp_proj_input',
                )
                if st.button("Compare", key='comp_btn'):
                    proj_list = [p.strip() for p in comp_projects.split(',') if p.strip()]
                    with st.spinner(f"Comparing {len(proj_list)} projects..."):
                        try:
                            tcga = TCGAConnector(cache=False)
                            comp = tcga.compare_projects(proj_list)
                            st.dataframe(comp, use_container_width=True, hide_index=True)
                        except Exception as e:
                            st.error(f"Failed: {e}")
    
    # Upload your own data
    st.markdown("---")
    st.markdown("### Or Upload Your Own Data")
    st.caption(
        "Have your own RNA-seq count matrix? Upload it here, then pool it with public studies "
        "in the Pool tab for cross-cohort biomarker validation — no command line needed."
    )

    col_up1, col_up2 = st.columns(2)
    
    with col_up1:
        uploaded_counts = st.file_uploader(
            "Count matrix (CSV/TSV)",
            type=['csv', 'tsv', 'txt'],
            key='user_counts_upload',
            help="Genes as rows, samples as columns"
        )
    
    with col_up2:
        uploaded_meta = st.file_uploader(
            "Sample metadata (CSV, optional)",
            type=['csv', 'tsv', 'txt'],
            key='user_meta_upload',
            help="Sample IDs as first column"
        )
    
    if uploaded_counts:
        if st.button("Import User Data"):
            try:
                # Save to temp and load
                import tempfile, os
                with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as f:
                    f.write(uploaded_counts.getvalue())
                    counts_path = f.name
                
                meta_path = None
                if uploaded_meta:
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as f:
                        f.write(uploaded_meta.getvalue())
                        meta_path = f.name
                
                ds = AcquiredDataset.from_csv(counts_path, meta_path)
                name = uploaded_counts.name.rsplit('.', 1)[0]
                st.session_state['acq_datasets'][name] = ds
                st.success(f"Imported: {ds.n_genes:,} {_feature_label(ds)}, {ds.n_samples} samples")
                
                # Cleanup
                os.unlink(counts_path)
                if meta_path:
                    os.unlink(meta_path)
                    
            except Exception as e:
                st.error(f"Import failed: {e}")

# =============================================================================
# TAB 2: Data Library
# =============================================================================

with tab_library:
    st.markdown("## Data Library")
    st.caption("Your personal collection of datasets — everything you've downloaded or uploaded")
    
    # Session datasets
    st.markdown("### Session Datasets")
    session_ds = st.session_state.get('acq_datasets', {})
    
    if session_ds:
        # Separate expression datasets from SRA run tables
        expr_rows = []
        sra_rows = []
        for name, ds in session_ds.items():
            is_sra = (
                ds.data_type == 'run_metadata'
                or (hasattr(ds, 'source_info') and ds.source_info.get('repository') == 'SRA')
            )
            if is_sra:
                src = ds.source_info if hasattr(ds, 'source_info') else {}
                gse = src.get('gse', '')
                gsm_ids = src.get('gsm_ids', [])
                if gse:
                    geo_col = gse
                elif gsm_ids:
                    geo_col = f'{len(gsm_ids)} GSMs found'
                else:
                    geo_col = ''
                sra_rows.append({
                    'Name': name,
                    'Repository': 'SRA',
                    'Runs': src.get('n_runs', ds.n_samples),
                    'Platform': src.get('platform', 'unknown'),
                    'Organism': ds.organism,
                    'GEO Link': geo_col,
                    'Data Type': 'Run metadata',
                })
            else:
                count_label = 'Cases' if (hasattr(ds, 'repository') and ds.repository == 'TCGA') else 'Samples'
                feat_label = _feature_label(ds).capitalize()
                expr_rows.append({
                    'Name': name,
                    'Repository': ds.repository,
                    feat_label: f"{ds.n_genes:,}",
                    count_label: ds.n_samples,
                    'Organism': ds.organism,
                    'Gene IDs': ds.gene_id_type,
                    'Data Type': ds.data_type,
                })
        
        if expr_rows:
            st.dataframe(pd.DataFrame(expr_rows), use_container_width=True, hide_index=True)
        
        if sra_rows:
            if expr_rows:
                st.markdown("**SRA Run Tables** (raw read metadata — not count matrices)")
            st.dataframe(pd.DataFrame(sra_rows), use_container_width=True, hide_index=True)
        
        # Dataset details
        selected_ds = st.selectbox(
            "Select dataset for details",
            list(session_ds.keys()),
            key='library_detail_select'
        )
        
        if selected_ds:
            ds = session_ds[selected_ds]
            
            is_sra_data = (
                ds.data_type == 'run_metadata'
                or (hasattr(ds, 'source_info') and ds.source_info.get('repository') == 'SRA')
            )
            
            if is_sra_data:
                # ----- SRA-specific detail view -----
                src = ds.source_info if hasattr(ds, 'source_info') else {}
                run_meta = ds.metadata
                
                # Compute total reads
                total_reads = 0
                if 'read_count' in run_meta.columns:
                    try:
                        total_reads = int(pd.to_numeric(run_meta['read_count'], errors='coerce').sum())
                    except Exception:
                        pass
                if total_reads > 1e9:
                    reads_str = f"{total_reads / 1e9:.1f}B"
                elif total_reads > 1e6:
                    reads_str = f"{total_reads / 1e6:.0f}M"
                else:
                    reads_str = f"{total_reads:,}"
                
                n_runs = src.get('n_runs', len(run_meta))
                platform = src.get('platform', 'unknown')
                organism = src.get('organism', ds.organism)
                
                # Professional info bar (consistent font size)
                st.markdown(
                    f"""<div style="
                        display: flex; gap: 2rem; padding: 0.8rem 0;
                        border-bottom: 1px solid #e0e0e0; margin-bottom: 0.8rem;
                    ">
                        <div><span style="color:#888; font-size:0.85rem;">Runs</span><br>
                            <span style="font-size:1.2rem; font-weight:600;">{n_runs}</span></div>
                        <div><span style="color:#888; font-size:0.85rem;">Platform</span><br>
                            <span style="font-size:1.2rem; font-weight:600;">{platform}</span></div>
                        <div><span style="color:#888; font-size:0.85rem;">Organism</span><br>
                            <span style="font-size:1.2rem; font-weight:600;">{organism}</span></div>
                        <div><span style="color:#888; font-size:0.85rem;">Total Reads</span><br>
                            <span style="font-size:1.2rem; font-weight:600;">{reads_str}</span></div>
                    </div>""",
                    unsafe_allow_html=True,
                )
                
                # Run table preview with relevant columns
                with st.expander("Run table", expanded=True):
                    run_display_cols = [
                        c for c in [
                            'run_accession', 'sample_title', 'sample_alias',
                            'instrument_model', 'library_layout', 'library_selection',
                            'read_count', 'base_count', 'scientific_name',
                        ] if c in run_meta.columns
                    ]
                    display_df = run_meta[run_display_cols] if run_display_cols else run_meta
                    st.dataframe(display_df.head(30), use_container_width=True, hide_index=True)
                    
                    if len(run_meta) > 30:
                        st.caption(f"Showing 30 of {len(run_meta)} runs")
                    
                    # Download run table as CSV
                    csv_data = run_meta.to_csv(index=True)
                    st.download_button(
                        "Download full run table (CSV)",
                        data=csv_data,
                        file_name=f"{selected_ds}_run_table.csv",
                        mime="text/csv",
                    )
                
                # SRA guidance with GSE link
                st.markdown("---")
                gse_linked = src.get('gse', '')
                if gse_linked:
                    st.info(
                        f"**This is an SRA run table** — it contains sequencing run metadata, "
                        f"not gene expression counts.\n\n"
                        f"**This study has a linked GEO accession: "
                        f"[{gse_linked}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_linked})**\n\n"
                        f"You can download processed counts directly from GEO using the accession "
                        f"`{gse_linked}` in the Search tab — no FASTQ download or quantification needed."
                    )
                else:
                    gsm_ids = src.get('gsm_ids', [])
                    if gsm_ids:
                        st.info(
                            f"**This is an SRA run table** — it contains sequencing run metadata, "
                            f"not gene expression counts.\n\n"
                            f"**{len(gsm_ids)} GEO sample(s) detected** (e.g., {gsm_ids[0]}). "
                            f"This study likely has processed counts on GEO."
                        )
                        if st.button(
                            "Look up linked GSE",
                            type="primary",
                            key=f'lookup_gse_{selected_ds}'
                        ):
                            with st.spinner("Querying NCBI for linked GEO Series..."):
                                try:
                                    connector = SRAConnector(cache=True, cache_dir=cache.cache_dir)
                                    gse_found = connector._gse_from_gsm(gsm_ids[0])
                                    if gse_found:
                                        ds.source_info['gse'] = gse_found
                                        st.success(
                                            f"Found linked GEO study: "
                                            f"[{gse_found}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_found}) "
                                            f"— use this accession in the Search tab to download processed counts."
                                        )
                                        st.rerun()
                                    else:
                                        st.warning(
                                            "Could not find a linked GSE. The GEO entry may not "
                                            "exist yet, or Biopython may not be installed."
                                        )
                                except Exception as e:
                                    st.error(f"GSE lookup failed: {e}")
                    else:
                        st.info(
                            "**This is an SRA run table** — it contains sequencing run metadata, "
                            "not gene expression counts. To get a count matrix:\n\n"
                            "1. Download FASTQ files using the SRA Toolkit (see below)\n"
                            "2. Quantify with RAPTOR Module 1 (Salmon/Kallisto) or Module 5 (STAR)\n"
                            "3. Upload the count matrix back here via **Search → Upload Your Own Data**\n"
                            "4. Pool with other studies, check batch effects, and continue analysis\n\n"
                            "Or search for the corresponding GSE on GEO for processed counts."
                        )
                
                # Download commands generation
                with st.expander("Generate FASTQ download commands"):
                    st.caption(
                        "Generate ready-to-run commands for downloading raw FASTQ files. "
                        "After downloading, quantify with RAPTOR Module 1 (Salmon/Kallisto) "
                        "or Module 5 (STAR + featureCounts) to produce a count matrix."
                    )
                    
                    col_os, col_tool, col_threads, col_outdir = st.columns([1, 1, 1, 2])
                    with col_os:
                        dl_os = st.selectbox(
                            "Platform",
                            ["Windows (PowerShell)", "Linux / macOS"],
                            key=f'sra_dl_os_{selected_ds}',
                            help="Choose your operating system"
                        )
                    with col_tool:
                        dl_tool = st.selectbox(
                            "Tool",
                            ["fasterq-dump", "fastq-dump"],
                            key=f'sra_dl_tool_{selected_ds}',
                            help="fasterq-dump is faster (multi-threaded) but outputs "
                                 "uncompressed files. fastq-dump can output gzipped files directly."
                        )
                    with col_threads:
                        dl_threads = st.number_input(
                            "Threads", min_value=1, max_value=32, value=4,
                            key=f'sra_dl_threads_{selected_ds}'
                        )
                    with col_outdir:
                        if 'Windows' in dl_os:
                            default_dir = ".\\fastq"
                        else:
                            default_dir = "./fastq"
                        dl_outdir = st.text_input(
                            "Output directory",
                            value=default_dir,
                            key=f'sra_dl_outdir_{selected_ds}'
                        )
                    
                    # Installation instructions
                    with st.expander("SRA Toolkit installation"):
                        if 'Windows' in dl_os:
                            st.markdown(
                                "**Windows setup:**\n\n"
                                "1. Download from [github.com/ncbi/sra-tools]"
                                "(https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) "
                                "(Windows 64-bit .zip)\n"
                                "2. Extract to a folder, e.g. `C:\\sra-toolkit`\n"
                                "3. Add to PATH in PowerShell:\n"
                                "   ```\n"
                                "   $env:PATH += \";C:\\sra-toolkit\\bin\"\n"
                                "   ```\n"
                                "4. Configure (first time only):\n"
                                "   ```\n"
                                "   vdb-config --interactive\n"
                                "   ```\n"
                                "   Enable 'Remote Access', then save and exit.\n"
                                "5. Test: `fasterq-dump --version`"
                            )
                        else:
                            st.markdown(
                                "**Linux / macOS setup:**\n\n"
                                "```bash\n"
                                "# Option A: Conda (recommended)\n"
                                "conda install -c bioconda sra-tools\n\n"
                                "# Option B: Manual\n"
                                "# Download from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit\n"
                                "# Extract and add bin/ to PATH\n"
                                "```\n"
                                "Then configure: `vdb-config --interactive` (enable Remote Access)"
                            )
                    
                    if st.button("Generate Commands", key=f'sra_gen_cmds_{selected_ds}'):
                        run_accs = []
                        if 'run_accession' in run_meta.columns:
                            run_accs = run_meta['run_accession'].tolist()
                        elif run_meta.index.name == 'sample_id':
                            run_accs = run_meta.index.tolist()
                        
                        if run_accs:
                            is_windows = 'Windows' in dl_os
                            cmds = []
                            for acc in run_accs:
                                if dl_tool == 'fasterq-dump':
                                    cmds.append(
                                        f"fasterq-dump --split-files --threads {dl_threads} "
                                        f"--outdir {dl_outdir} {acc}"
                                    )
                                else:
                                    cmds.append(
                                        f"fastq-dump --split-files --gzip "
                                        f"--outdir {dl_outdir} {acc}"
                                    )
                            
                            if is_windows:
                                script = f"# Download {len(cmds)} FASTQ files for {selected_ds}\n"
                                script += f"# Run in PowerShell after installing SRA Toolkit\n\n"
                                script += f"New-Item -ItemType Directory -Force -Path {dl_outdir}\n\n"
                                for i, cmd in enumerate(cmds, 1):
                                    script += f"Write-Host \"[{i}/{len(cmds)}] Downloading...\"\n"
                                    script += cmd + "\n\n"
                                script += f'Write-Host "Done. {len(cmds)} runs downloaded to {dl_outdir}"\n'
                                
                                st.code(script, language="powershell")
                                st.download_button(
                                    "Download as PowerShell script",
                                    data=script,
                                    file_name=f"download_{selected_ds}.ps1",
                                    mime="text/plain",
                                )
                            else:
                                script = "#!/bin/bash\n"
                                script += f"# Download {len(cmds)} FASTQ files for {selected_ds}\n"
                                script += f"set -e\n\n"
                                script += f"mkdir -p {dl_outdir}\n\n"
                                for i, cmd in enumerate(cmds, 1):
                                    script += f"echo \"[{i}/{len(cmds)}] Downloading...\"\n"
                                    script += cmd + "\n\n"
                                script += f'echo "Done. {len(cmds)} runs downloaded to {dl_outdir}"\n'
                                
                                st.code(script, language="bash")
                                st.download_button(
                                    "Download as shell script",
                                    data=script,
                                    file_name=f"download_{selected_ds}.sh",
                                    mime="text/x-shellscript",
                                )
                            
                            # Next steps guidance
                            st.markdown("---")
                            st.markdown("**Next steps after downloading FASTQs:**")
                            st.markdown(
                                f"1. Run the script above to download {len(run_accs)} FASTQ file(s)\n"
                                f"2. Quantify reads into a count matrix:\n"
                                f"   - **Module 1** — fast pseudo-alignment: "
                                f"`raptor quantify --tool salmon` or `raptor quantify --tool kallisto`\n"
                                f"   - **Module 5** — splice-aware alignment: "
                                f"STAR + featureCounts (for novel splice junctions or genome-level analysis)\n"
                                f"3. Upload the resulting count matrix back here in "
                                f"**Data Acquisition → Search → Upload Your Own Data**\n"
                                f"4. If combining with other studies, use the **Pool Datasets** tab to merge, "
                                f"check coverage, detect batch effects, and apply correction\n"
                                f"5. Continue with **Quality Assessment**, **Differential Expression**, etc."
                            )
                        else:
                            st.warning("No run accessions found in the run table")
            
            else:
                # ----- Dataset detail view (expression, miRNA, methylation, CNV, protein) -----
                is_tcga_ds = hasattr(ds, 'repository') and ds.repository == 'TCGA'
                is_large = ds.n_samples > 200  # Large dataset guard
                is_mutation = hasattr(ds, 'data_type') and 'mutation' in str(ds.data_type).lower()
                feat = _feature_label(ds)
                mat = _matrix_label(ds)
                
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric(feat.capitalize(), f"{ds.n_genes:,}")
                with col2:
                    st.metric("Cases" if is_tcga_ds else "Samples", ds.n_samples)
                with col3:
                    st.metric("Organism", ds.organism)
                with col4:
                    if is_large or is_mutation:
                        st.metric("Integrity", "—")
                        if is_mutation:
                            st.caption("N/A for mutation data")
                        else:
                            st.caption("Skipped for large dataset")
                    else:
                        integrity = ds.validate_integrity()
                        n_errs = len(integrity['errors'])
                        status = "OK" if integrity['valid'] else f"{n_errs} {'error' if n_errs == 1 else 'errors'}"
                        st.metric("Integrity", status)
                
                if not is_large and not is_mutation and 'integrity' in dir() and integrity.get('warnings'):
                    for w in integrity['warnings']:
                        st.warning(w)
                
                if is_mutation:
                    st.info(
                        f"**Somatic mutation data (MAF format).** "
                        f"{ds.n_genes:,} mutations across {ds.n_samples} samples. "
                        f"Use this for mutation frequency analysis, not for pooling or QC plots."
                    )
                elif ds.gene_id_type == 'probe':
                    st.info(
                        f"**DNA methylation data** — {ds.n_genes:,} CpG probes across {ds.n_samples} samples. "
                        f"Values are beta values (0–1), not read counts. "
                        f"NaN values indicate probes that failed detection in a sample — this is expected."
                    )
                elif ds.gene_id_type == 'mirna':
                    st.info(
                        f"**miRNA expression data** — {ds.n_genes:,} miRNAs across {ds.n_samples} samples. "
                        f"Values are reads per million (RPM)."
                    )
                elif ds.gene_id_type == 'protein':
                    st.info(
                        f"**Protein expression data (RPPA)** — {ds.n_genes:,} proteins across {ds.n_samples} samples."
                    )
                elif is_large:
                    st.info(
                        f"**Large dataset ({ds.n_samples} samples).** "
                        f"Use the **Quality Check** tab for PCA, heatmaps, and correlation analysis. "
                        f"Preview limited to avoid browser freeze."
                    )
                
                # Preview — limit for large datasets, collapsed by default
                if is_mutation:
                    maf = ds.counts_df
                    maf_cols = maf.columns.tolist()
                    
                    with st.expander("Preview mutation table (MAF)"):
                        # Show key columns only
                        key_maf_cols = [c for c in [
                            'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
                            'Variant_Classification', 'Variant_Type', 'Reference_Allele',
                            'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'HGVSp_Short',
                        ] if c in maf_cols]
                        if key_maf_cols:
                            st.dataframe(maf[key_maf_cols].head(20), use_container_width=True)
                            st.caption(f"Showing {len(key_maf_cols)} key columns of {len(maf_cols)} total")
                        else:
                            st.dataframe(maf.head(20), use_container_width=True)
                    
                    with st.expander("Mutation Summary"):
                        if 'Hugo_Symbol' in maf_cols:
                            top_genes = maf['Hugo_Symbol'].value_counts().head(20)
                            st.markdown("**Top 20 most mutated genes:**")
                            st.dataframe(top_genes.reset_index().rename(
                                columns={'index': 'Gene', 'Hugo_Symbol': 'Gene', 'count': 'Mutations'}
                            ), use_container_width=True, hide_index=True)
                        if 'Variant_Classification' in maf_cols:
                            st.markdown("**Variant types:**")
                            st.dataframe(maf['Variant_Classification'].value_counts().reset_index().rename(
                                columns={'index': 'Type', 'Variant_Classification': 'Type', 'count': 'Count'}
                            ), use_container_width=True, hide_index=True)
                        if 'Tumor_Sample_Barcode' in maf_cols:
                            n_patients = maf['Tumor_Sample_Barcode'].nunique()
                            mutations_per_patient = maf['Tumor_Sample_Barcode'].value_counts()
                            st.markdown(f"**Patients with mutations:** {n_patients}")
                            st.markdown(f"**Median mutations per patient:** {mutations_per_patient.median():.0f}")
                            st.markdown(f"**Range:** {mutations_per_patient.min()} – {mutations_per_patient.max()}")
                elif not is_large:
                    with st.expander(f"Preview {mat}"):
                        st.dataframe(ds.counts_df.head(20), use_container_width=True)
                    
                    with st.expander("Preview metadata"):
                        st.dataframe(ds.metadata.head(20), use_container_width=True)
                else:
                    with st.expander(f"Preview {mat} (20×20)"):
                        st.dataframe(ds.counts_df.iloc[:20, :20], use_container_width=True)
                    
                    with st.expander("Preview metadata (first 20 rows, key columns)"):
                        # Only show key columns for large datasets
                        key_cols = [c for c in ds.metadata.columns if c in (
                            'sample_id', 'sample_type', 'sample_category', 'case_id',
                            'gender', 'vital_status', 'ajcc_stage', 'project'
                        )]
                        if not key_cols:
                            key_cols = ds.metadata.columns[:8].tolist()
                        st.dataframe(ds.metadata[key_cols].head(20), use_container_width=True)
                        st.caption(f"Showing {len(key_cols)} of {len(ds.metadata.columns)} columns")
                
                # ---- Sample Filtering (TCGA) ----
                meta = ds.metadata
                is_tcga_data = hasattr(ds, 'repository') and ds.repository == 'TCGA'
                has_pairing_info = (
                    not meta.empty
                    and 'sample_category' in meta.columns
                    and 'case_id' in meta.columns
                )
                
                if has_pairing_info:
                    st.markdown("---")
                    st.markdown("### Sample Filtering")
                    
                    if is_large:
                        # For large datasets, compute only on button click
                        if st.button("Analyze sample pairing", key=f'analyze_pairing_{selected_ds}'):
                            with st.spinner("Analyzing paired samples..."):
                                tumor_cases = set(meta.loc[meta['sample_category'] == 'Tumor', 'case_id'])
                                normal_cases = set(meta.loc[meta['sample_category'] == 'Normal', 'case_id'])
                                paired_cases = tumor_cases & normal_cases
                                tumor_only_cases = tumor_cases - normal_cases
                                normal_only_cases = normal_cases - tumor_cases
                                st.session_state[f'pairing_{selected_ds}'] = {
                                    'tumor_cases': tumor_cases, 'normal_cases': normal_cases,
                                    'paired_cases': paired_cases, 'tumor_only_cases': tumor_only_cases,
                                    'normal_only_cases': normal_only_cases,
                                }
                        
                        pairing_data = st.session_state.get(f'pairing_{selected_ds}')
                        if pairing_data is None:
                            n_tumor = (meta['sample_category'] == 'Tumor').sum() if 'sample_category' in meta.columns else 0
                            n_normal = (meta['sample_category'] == 'Normal').sum() if 'sample_category' in meta.columns else 0
                            st.markdown(
                                f"**Current dataset:** {len(meta)} samples "
                                f"({n_tumor} tumor, {n_normal} normal). "
                                f"Click above to analyze paired samples."
                            )
                            has_pairing_info = False  # skip rest of section
                        else:
                            tumor_cases = pairing_data['tumor_cases']
                            normal_cases = pairing_data['normal_cases']
                            paired_cases = pairing_data['paired_cases']
                            tumor_only_cases = pairing_data['tumor_only_cases']
                            normal_only_cases = pairing_data['normal_only_cases']
                    else:
                        # Small dataset — compute immediately
                        tumor_cases = set(meta.loc[meta['sample_category'] == 'Tumor', 'case_id'])
                        normal_cases = set(meta.loc[meta['sample_category'] == 'Normal', 'case_id'])
                        paired_cases = tumor_cases & normal_cases
                        tumor_only_cases = tumor_cases - normal_cases
                        normal_only_cases = normal_cases - tumor_cases
                    
                    # Only show filtering UI if pairing data is computed
                    if 'paired_cases' in dir() and paired_cases is not None:
                        n_tumor = (meta['sample_category'] == 'Tumor').sum()
                        n_normal = (meta['sample_category'] == 'Normal').sum()
                    
                        paired_mask = meta['case_id'].isin(paired_cases)
                        paired_tumor = ((meta['sample_category'] == 'Tumor') & paired_mask).sum()
                        paired_normal = ((meta['sample_category'] == 'Normal') & paired_mask).sum()
                        
                        st.markdown(
                            f"**Current dataset:** {len(meta)} samples "
                            f"({n_tumor} tumor, {n_normal} normal) from "
                            f"{meta['case_id'].nunique()} patients"
                        )
                        
                        col_p1, col_p2, col_p3 = st.columns(3)
                        with col_p1:
                            st.metric("Paired Patients", len(paired_cases))
                            st.caption(f"{paired_tumor} tumor + {paired_normal} normal")
                        with col_p2:
                            st.metric("Tumor-only Patients", len(tumor_only_cases))
                        with col_p3:
                            st.metric("Normal-only Patients", len(normal_only_cases))
                        
                        filter_action = st.radio(
                            "Filter samples",
                            [
                                "Keep all (no filter)",
                                f"Paired only — {len(paired_cases)} patients ({paired_tumor + paired_normal} samples)",
                                f"Tumor only — {n_tumor} samples",
                                f"Normal only — {n_normal} samples",
                                "Custom — choose patients",
                            ],
                            key=f'sample_filter_action_{selected_ds}',
                        )
                        
                        if filter_action.startswith("Custom"):
                            all_patients = sorted(meta['case_id'].unique())
                            selected_patients = st.multiselect(
                                "Select patients to keep",
                                all_patients,
                                default=list(paired_cases) if paired_cases else all_patients,
                                format_func=lambda p: f"{p} [paired]" if p in paired_cases
                                    else f"{p} [{'tumor' if p in tumor_only_cases else 'normal'}]",
                                key=f'custom_patient_select_{selected_ds}',
                            )
                        
                        if not filter_action.startswith("Keep all"):
                            if filter_action.startswith("Paired"):
                                keep_mask = meta['case_id'].isin(paired_cases)
                                new_suffix = "paired"
                            elif filter_action.startswith("Tumor"):
                                keep_mask = meta['sample_category'] == 'Tumor'
                                new_suffix = "tumor"
                            elif filter_action.startswith("Normal"):
                                keep_mask = meta['sample_category'] == 'Normal'
                                new_suffix = "normal"
                            elif filter_action.startswith("Custom") and 'selected_patients' in dir():
                                keep_mask = meta['case_id'].isin(selected_patients)
                                new_suffix = "custom"
                            else:
                                keep_mask = pd.Series([True] * len(meta), index=meta.index)
                                new_suffix = "filtered"
                            
                            keep_samples = meta.index[keep_mask].tolist()
                            n_keep = len(keep_samples)
                            
                            st.info(f"Will keep **{n_keep}** of {len(meta)} samples")
                            
                            if n_keep > 0 and n_keep < len(meta):
                                with st.expander(f"Preview filtered samples ({n_keep})"):
                                    filtered_meta = meta.loc[keep_mask]
                                    st.dataframe(filtered_meta.head(50), use_container_width=True)
                            
                            new_name = f"{selected_ds}_{new_suffix}"
                            
                            if st.button(
                                f"Create filtered dataset: {new_name}",
                                type="primary",
                                key=f'apply_filter_{selected_ds}',
                            ):
                                valid_samples = [s for s in keep_samples if s in ds.counts_df.columns]
                                if valid_samples:
                                    new_counts = ds.counts_df[valid_samples]
                                    new_meta = meta.loc[meta.index.isin(valid_samples)]
                                    
                                    new_ds = AcquiredDataset(
                                        counts_df=new_counts,
                                        metadata=new_meta,
                                        source_info={
                                            **ds.source_info,
                                            'filtered_from': selected_ds,
                                            'filter': new_suffix,
                                            'n_original_samples': len(meta),
                                        },
                                        gene_id_type=ds.gene_id_type,
                                    )
                                    
                                    st.session_state['acq_datasets'][new_name] = new_ds
                                    st.success(
                                        f"Created **{new_name}**: "
                                        f"{new_ds.n_genes:,} {_feature_label(new_ds)}, {new_ds.n_samples} samples "
                                        f"({new_ds.metadata['sample_category'].value_counts().to_dict()})"
                                    )
                                    st.rerun()
                                else:
                                    st.error("No matching samples found in data matrix")
            
            if not is_sra_data and not is_mutation:
                st.markdown("---")
                st.markdown("### Gene ID Conversion")
            
                detected_type = ds.gene_id_type
                gene_sample = ds.counts_df.index[:5].tolist()
            
                col_det, col_target, col_convert = st.columns([2, 2, 1])
            
                with col_det:
                    st.markdown(f"**Current type:** `{detected_type}`")
                    st.caption(f"Examples: {', '.join(str(g) for g in gene_sample)}")
            
                with col_target:
                    # Only show conversion options that differ from current
                    target_options = [t for t in ['symbol', 'ensembl', 'entrez'] if t != detected_type]
                    target_type = st.selectbox(
                        "Convert to",
                        target_options,
                        key=f'geneid_target_{selected_ds}',
                        help="Choose the target gene ID type"
                    )
            
                with col_convert:
                    st.markdown("<br>", unsafe_allow_html=True)
                    convert_clicked = st.button(
                        "Convert",
                        type="primary",
                        use_container_width=True,
                        key=f'convert_btn_{selected_ds}'
                    )
            
                # Show previous conversion result if available
                conv_key = f'conversion_result_{selected_ds}'
            
                if convert_clicked:
                    with st.spinner(f"Converting {detected_type} → {target_type} via MyGene.info..."):
                        try:
                            mapper = GeneIDMapper(species=ds.organism)
                        
                            gene_ids = ds.counts_df.index.astype(str).tolist()
                            mapping = mapper.convert(gene_ids, detected_type, target_type)
                        
                            # Calculate stats
                            total = len(mapping)
                            mapped = sum(1 for v in mapping.values() if v is not None)
                            failed = total - mapped
                            pct = 100 * mapped / total if total > 0 else 0
                        
                            # Find duplicates (multiple source IDs → same target)
                            target_vals = [v for v in mapping.values() if v is not None]
                            n_unique = len(set(target_vals))
                            n_dups = mapped - n_unique
                        
                            # Failed genes
                            failed_genes = [k for k, v in mapping.items() if v is None]
                        
                            # Store result
                            st.session_state[conv_key] = {
                                'mapping': mapping,
                                'total': total,
                                'mapped': mapped,
                                'failed': failed,
                                'pct': pct,
                                'n_dups': n_dups,
                                'failed_genes': failed_genes,
                                'from_type': detected_type,
                                'to_type': target_type,
                            }
                        
                        except Exception as e:
                            st.error(f"Conversion failed: {e}")
            
                # Display conversion results
                if conv_key in st.session_state:
                    cr = st.session_state[conv_key]
                
                    # Stats row
                    col_s1, col_s2, col_s3, col_s4 = st.columns(4)
                    with col_s1:
                        st.metric("Total Genes", f"{cr['total']:,}")
                    with col_s2:
                        color = "normal" if cr['pct'] > 80 else "off"
                        st.metric("Mapped", f"{cr['mapped']:,}", f"{cr['pct']:.1f}%")
                    with col_s3:
                        st.metric("Failed", f"{cr['failed']:,}")
                    with col_s4:
                        st.metric("Duplicates", f"{cr['n_dups']:,}")
                
                    # Verdict
                    if cr['pct'] >= 90:
                        st.success(
                            f"Excellent conversion: {cr['pct']:.1f}% mapped "
                            f"({cr['from_type']} → {cr['to_type']})"
                        )
                    elif cr['pct'] >= 70:
                        st.warning(
                            f"Moderate conversion: {cr['pct']:.1f}% mapped. "
                            f"{cr['failed']:,} genes could not be converted."
                        )
                    else:
                        st.error(
                            f"Low conversion: {cr['pct']:.1f}% mapped. "
                            f"Check that the detected ID type ({cr['from_type']}) is correct "
                            f"and the organism matches."
                        )
                
                    # Duplicates info
                    if cr['n_dups'] > 0:
                        st.info(
                            f"{cr['n_dups']} genes map to the same {cr['to_type']} ID. "
                            f"These will be aggregated (mean) during conversion."
                        )
                
                    # Failed genes
                    if cr['failed_genes']:
                        with st.expander(f"Failed genes ({cr['failed']:,})"):
                            failed_df = pd.DataFrame({
                                'Gene ID': cr['failed_genes'][:500],
                            })
                            if len(cr['failed_genes']) > 500:
                                st.caption(f"Showing first 500 of {len(cr['failed_genes']):,}")
                            st.dataframe(failed_df, use_container_width=True, hide_index=True)
                        
                            # Download failed genes
                            csv_failed = '\n'.join(cr['failed_genes'])
                            st.download_button(
                                "Download failed genes list",
                                data=csv_failed,
                                file_name=f"{selected_ds}_unmapped_genes.txt",
                                mime="text/plain",
                            )
                
                    # Apply conversion button
                    st.markdown("")
                    col_apply, col_preview_conv = st.columns(2)
                
                    with col_apply:
                        if st.button(
                            f"Apply: Convert {selected_ds} to {cr['to_type']}",
                            type="primary",
                            use_container_width=True,
                            key=f'apply_conv_{selected_ds}'
                        ):
                            with st.spinner("Applying conversion..."):
                                try:
                                    mapper = GeneIDMapper(species=ds.organism)
                                    new_counts = mapper.convert_index(
                                        ds.counts_df,
                                        from_type=cr['from_type'],
                                        to_type=cr['to_type'],
                                        drop_unmapped=True,
                                        aggregate='mean',
                                    )
                                
                                    # Update the dataset in session
                                    from raptor.external_modules.acquisition import AcquiredDataset
                                    new_source = {
                                        **ds.source_info,
                                        'accession': selected_ds,
                                        'original_gene_id_type': cr['from_type'],
                                        'converted_to': cr['to_type'],
                                    }
                                    new_ds = AcquiredDataset(
                                        counts_df=new_counts,
                                        metadata=ds.metadata.copy(),
                                        source_info=new_source,
                                        gene_id_type=cr['to_type'],
                                    )
                                
                                    st.session_state['acq_datasets'][selected_ds] = new_ds
                                
                                    # Also save converted version to cache so it persists
                                    try:
                                        cache.save_dataset(
                                            counts_df=new_counts,
                                            metadata=ds.metadata.copy(),
                                            source_info=new_source,
                                            gene_id_type=cr['to_type'],
                                        )
                                    except Exception:
                                        pass  # cache save is best-effort
                                
                                    # Clear conversion result
                                    del st.session_state[conv_key]
                                
                                    st.success(
                                        f"Converted & saved {selected_ds}: "
                                        f"{new_counts.shape[0]:,} genes ({cr['to_type']}). "
                                        f"Lost {cr['failed']:,} unmapped + {cr['n_dups']} duplicates merged."
                                    )
                                    st.rerun()
                                
                                except Exception as e:
                                    st.error(f"Conversion failed: {e}")
                
                    with col_preview_conv:
                        if st.button(
                            "Preview converted genes",
                            use_container_width=True,
                            key=f'preview_conv_{selected_ds}'
                        ):
                            mapping = cr['mapping']
                            preview_rows = []
                            for orig, converted in list(mapping.items())[:20]:
                                preview_rows.append({
                                    f'Original ({cr["from_type"]})': orig,
                                    f'Converted ({cr["to_type"]})': converted if converted else '❌ unmapped',
                                })
                            st.dataframe(
                                pd.DataFrame(preview_rows),
                                use_container_width=True,
                                hide_index=True,
                            )
    else:
        st.info("No datasets loaded yet. Use the Search tab to download data or upload your own.")
    
    
    # =========================================================================
    # Metadata Editor
    # =========================================================================
    
    # Filter out mutation datasets from metadata editor
    expression_ds = {
        name: ds for name, ds in session_ds.items()
        if not (hasattr(ds, 'data_type') and 'mutation' in str(ds.data_type).lower())
    }
    
    if expression_ds:
        st.markdown("---")
        st.markdown("### Sample Metadata Editor")
        st.caption(
            "Create a unified metadata table across all datasets. "
            "Assign conditions, batches, and groups. Remove unwanted samples."
        )
        
        # Build combined metadata from all session datasets
        if 'meta_editor_df' not in st.session_state or st.button(
            "Refresh from datasets", key='meta_refresh'
        ):
            combined_rows = []
            for ds_name, ds in session_ds.items():
                meta = ds.metadata.copy()
                meta = meta.reset_index()
                
                # Ensure sample_id column exists
                if 'sample_id' not in meta.columns:
                    if meta.index.name == 'sample_id':
                        meta = meta.reset_index()
                    else:
                        meta.insert(0, 'sample_id', meta.index.astype(str))
                
                # Add dataset source column
                meta.insert(0, 'dataset', ds_name)
                
                combined_rows.append(meta)
            
            if combined_rows:
                combined_meta = pd.concat(combined_rows, ignore_index=True)
            else:
                combined_meta = pd.DataFrame(columns=['dataset', 'sample_id'])
            
            # Add default user columns if not present
            if 'batch' not in combined_meta.columns:
                # Auto-assign batch: each study = unique batch with clean label
                ds_names = combined_meta['dataset'].unique().tolist()
                batch_map = {ds: f"Batch_{i+1}" for i, ds in enumerate(ds_names)}
                combined_meta['batch'] = combined_meta['dataset'].map(batch_map)
            if 'condition' not in combined_meta.columns:
                combined_meta['condition'] = ''
            if 'group' not in combined_meta.columns:
                combined_meta['group'] = ''
            
            # Initialize excluded samples
            if 'excluded_samples' not in st.session_state:
                st.session_state['excluded_samples'] = set()
            
            st.session_state['meta_editor_df'] = combined_meta
        
        meta_df = st.session_state['meta_editor_df']
        
        # =================================================================
        # Add new column
        # =================================================================
        col_add1, col_add2, col_add3 = st.columns([2, 2, 1])
        
        with col_add1:
            new_col_name = st.text_input(
                "New column name",
                placeholder="e.g., treatment, timepoint, tissue",
                key='new_col_input'
            )
        
        with col_add2:
            default_value = st.text_input(
                "Default value (optional)",
                placeholder="e.g., Control",
                key='new_col_default'
            )
        
        with col_add3:
            st.markdown("<br>", unsafe_allow_html=True)
            if st.button("Add Column", use_container_width=True, key='add_col_btn'):
                if new_col_name and new_col_name not in meta_df.columns:
                    meta_df[new_col_name] = default_value if default_value else ''
                    st.session_state['meta_editor_df'] = meta_df
                    st.rerun()
                elif new_col_name in meta_df.columns:
                    st.warning(f"Column '{new_col_name}' already exists")
                else:
                    st.warning("Enter a column name")
        
        with st.expander("Manage columns"):
            user_columns = [c for c in meta_df.columns if c not in ['dataset', 'sample_id']]
            cols_to_remove = st.multiselect(
                "Select columns to remove",
                user_columns,
                key='cols_remove_select'
            )
            if cols_to_remove and st.button("Remove selected columns", key='remove_cols_btn'):
                meta_df = meta_df.drop(columns=cols_to_remove)
                st.session_state['meta_editor_df'] = meta_df
                st.rerun()
        
        # =================================================================
        # Color-coded display
        # =================================================================
        st.markdown("---")
        
        # Color palette for datasets
        DATASET_COLORS = [
            '#E3F2FD', '#FFF3E0', '#E8F5E9', '#FCE4EC', '#F3E5F5',
            '#E0F7FA', '#FFFDE7', '#EFEBE9', '#F1F8E9', '#EDE7F6',
        ]
        
        unique_datasets = meta_df['dataset'].unique().tolist()
        dataset_color_map = {
            ds: DATASET_COLORS[i % len(DATASET_COLORS)]
            for i, ds in enumerate(unique_datasets)
        }
        
        # Color legend
        if len(unique_datasets) > 1:
            legend_parts = []
            for ds in unique_datasets:
                color = dataset_color_map[ds]
                n_samples = len(meta_df[meta_df['dataset'] == ds])
                legend_parts.append(
                    f'<span style="background-color:{color}; padding:2px 8px; '
                    f'border-radius:3px; margin-right:8px;">'
                    f'<b>{ds}</b> ({n_samples})</span>'
                )
            st.markdown("**Datasets:** " + " ".join(legend_parts), unsafe_allow_html=True)
        
        # Mark excluded samples
        excluded = st.session_state.get('excluded_samples', set())
        display_df = meta_df.copy()
        display_df['included'] = display_df.apply(
            lambda r: '❌' if (r['dataset'], r['sample_id']) in excluded else '✅', axis=1
        )
        
        # Styled view
        def color_by_dataset(row):
            color = dataset_color_map.get(row['dataset'], '#FFFFFF')
            if (row['dataset'], row['sample_id']) in excluded:
                return ['background-color: #F5F5F5; color: #BBBBBB; text-decoration: line-through'] * len(row)
            return [f'background-color: {color}'] * len(row)
        
        styled_df = display_df.style.apply(color_by_dataset, axis=1)
        st.dataframe(styled_df, use_container_width=True, hide_index=True, height=400)
        
        # =================================================================
        # Bulk Edit (select rows → apply value)
        # =================================================================
        st.markdown("### Bulk Edit Samples")
        st.caption("Select multiple samples, then assign a value to any column at once")
        
        # Filter by dataset first for easier selection
        bulk_dataset_filter = st.selectbox(
            "Filter by dataset",
            ['All datasets'] + unique_datasets,
            key='bulk_ds_filter'
        )
        
        if bulk_dataset_filter == 'All datasets':
            available_samples = meta_df[['dataset', 'sample_id']].apply(
                lambda r: f"{r['dataset']} | {r['sample_id']}", axis=1
            ).tolist()
        else:
            mask = meta_df['dataset'] == bulk_dataset_filter
            available_samples = meta_df.loc[mask, ['dataset', 'sample_id']].apply(
                lambda r: f"{r['dataset']} | {r['sample_id']}", axis=1
            ).tolist()
        
        selected_samples = st.multiselect(
            "Select samples",
            available_samples,
            key='bulk_sample_select',
            help="Select the samples you want to edit"
        )
        
        if selected_samples:
            st.caption(f"{len(selected_samples)} sample(s) selected")
            
            editable_cols = [c for c in meta_df.columns if c not in ['dataset', 'sample_id']]
            
            col_col, col_val, col_apply = st.columns([2, 2, 1])
            
            with col_col:
                edit_column = st.selectbox(
                    "Column to edit",
                    editable_cols,
                    key='bulk_edit_col'
                )
            
            with col_val:
                edit_value = st.text_input(
                    "Value to assign",
                    placeholder="e.g., Disorder, Treatment, Brain",
                    key='bulk_edit_val'
                )
            
            with col_apply:
                st.markdown("<br>", unsafe_allow_html=True)
                if st.button("Apply to Selected", type="primary", use_container_width=True, key='bulk_apply_btn'):
                    if edit_value is not None and edit_column:
                        # Parse selected samples back to dataset + sample_id
                        for sample_str in selected_samples:
                            ds_name, sample_id = sample_str.split(' | ', 1)
                            mask = (meta_df['dataset'] == ds_name) & (meta_df['sample_id'] == sample_id)
                            meta_df.loc[mask, edit_column] = edit_value
                        
                        st.session_state['meta_editor_df'] = meta_df
                        st.success(f"Set '{edit_column}' = '{edit_value}' for {len(selected_samples)} samples")
                        st.rerun()
            
            # =============================================================
            # Exclude / Include samples
            # =============================================================
            col_exc, col_inc = st.columns(2)
            
            with col_exc:
                if st.button("Exclude Selected Samples", use_container_width=True, key='exclude_btn'):
                    for sample_str in selected_samples:
                        ds_name, sample_id = sample_str.split(' | ', 1)
                        excluded.add((ds_name, sample_id))
                    st.session_state['excluded_samples'] = excluded
                    st.success(f"Excluded {len(selected_samples)} samples")
                    st.rerun()
            
            with col_inc:
                if st.button("Include Selected Samples", use_container_width=True, key='include_btn'):
                    for sample_str in selected_samples:
                        ds_name, sample_id = sample_str.split(' | ', 1)
                        excluded.discard((ds_name, sample_id))
                    st.session_state['excluded_samples'] = excluded
                    st.success(f"Re-included {len(selected_samples)} samples")
                    st.rerun()
        
        # Show exclusion summary
        if excluded:
            st.warning(
                f"{len(excluded)} sample(s) excluded. "
                f"These will be removed from pooling and count matrices."
            )
            
            if st.button("Apply Exclusions to Datasets", type="primary", key='apply_exclusions_btn'):
                for ds_name, ds in list(session_ds.items()):
                    # Find excluded samples for this dataset
                    ds_excluded = [s for d, s in excluded if d == ds_name]
                    if ds_excluded:
                        # Remove from counts
                        cols_to_keep = [c for c in ds.counts_df.columns if c not in ds_excluded]
                        if cols_to_keep:
                            new_counts = ds.counts_df[cols_to_keep]
                            # Remove from metadata
                            new_meta = ds.metadata.copy()
                            if hasattr(new_meta, 'index'):
                                new_meta = new_meta.loc[new_meta.index.isin(cols_to_keep)]
                            
                            from raptor.external_modules.acquisition import AcquiredDataset
                            new_ds = AcquiredDataset(
                                counts_df=new_counts,
                                metadata=new_meta,
                                source_info=ds.source_info,
                                gene_id_type=ds.gene_id_type,
                            )
                            st.session_state['acq_datasets'][ds_name] = new_ds
                
                # Clear exclusions and refresh
                st.session_state['excluded_samples'] = set()
                if 'meta_editor_df' in st.session_state:
                    del st.session_state['meta_editor_df']
                st.success("Exclusions applied — samples removed from count matrices and metadata")
                st.rerun()
        
        # =================================================================
        # Manual edit (fallback for individual cells)
        # =================================================================
        with st.expander("Edit individual cells"):
            column_config = {
                'dataset': st.column_config.TextColumn('Dataset', disabled=True),
                'sample_id': st.column_config.TextColumn('Sample ID', disabled=True),
            }
            
            edited_df = st.data_editor(
                meta_df,
                use_container_width=True,
                hide_index=True,
                column_config=column_config,
                num_rows="fixed",
                key='meta_data_editor',
            )
            st.session_state['meta_editor_df'] = edited_df
        
        # =================================================================
        # Export
        # =================================================================
        export_df = st.session_state.get('meta_editor_df', meta_df)
        col_exp1, col_exp2, col_exp3 = st.columns(3)
        
        with col_exp1:
            csv_data = export_df.to_csv(index=False)
            st.download_button(
                "Download Metadata (CSV)",
                data=csv_data,
                file_name="sample_metadata.csv",
                mime="text/csv",
                use_container_width=True,
            )
        
        with col_exp2:
            if st.button("Store for Downstream", use_container_width=True, key='store_meta_btn'):
                st.session_state['m6b_combined_metadata'] = export_df
                st.success("Metadata stored for downstream pages")
        
        with col_exp3:
            n_datasets = export_df['dataset'].nunique()
            n_samples = len(export_df)
            n_excluded = len(excluded)
            st.caption(f"{n_datasets} dataset(s) | {n_samples} sample(s) | {n_excluded} excluded")
    
    # =========================================================================
    # Mutation Annotation — cross-reference MAF with expression metadata
    # =========================================================================
    
    # Detect expression + mutation pairs for the same TCGA project
    session_ds_all = st.session_state.get('acq_datasets', {})
    
    expr_datasets = {}
    mut_datasets = {}
    for name, ds in session_ds_all.items():
        is_mut = hasattr(ds, 'data_type') and 'mutation' in str(ds.data_type).lower()
        is_sra = (ds.data_type == 'run_metadata' 
                  or (hasattr(ds, 'source_info') and ds.source_info.get('repository') == 'SRA'))
        if is_mut:
            # Extract project ID (e.g., "TCGA-BRCA" from "TCGA-BRCA_mutations")
            proj = ds.source_info.get('accession', name.replace('_mutations', ''))
            mut_datasets[proj] = (name, ds)
        elif not is_sra and ds.source_info.get('repository') == 'TCGA':
            proj = ds.source_info.get('accession', name)
            expr_datasets[proj] = (name, ds)
    
    # Find matching pairs
    matched_projects = set(expr_datasets.keys()) & set(mut_datasets.keys())
    
    if matched_projects:
        st.markdown("---")
        st.markdown("### Mutation Annotation")
        st.caption(
            "Cross-reference somatic mutations with expression data to add mutation status "
            "columns (e.g., TP53_status = Mutant/Wild-type) to your metadata. "
            "These columns become available in QC → Color by for mutation-driven analysis."
        )
        
        for proj in sorted(matched_projects):
            expr_name, expr_ds = expr_datasets[proj]
            mut_name, mut_ds = mut_datasets[proj]
            maf_df = mut_ds.counts_df  # MAF DataFrame stored as counts_df
            
            # Validate MAF has required columns
            if 'Hugo_Symbol' not in maf_df.columns or 'Tumor_Sample_Barcode' not in maf_df.columns:
                st.warning(f"{proj}: MAF data missing Hugo_Symbol or Tumor_Sample_Barcode columns.")
                continue
            
            st.markdown(f"**{proj}** — {expr_ds.n_samples} expression samples + "
                        f"{maf_df['Tumor_Sample_Barcode'].nunique()} mutation cases")
            
            # Top mutated genes from MAF
            gene_mut_counts = maf_df['Hugo_Symbol'].value_counts()
            
            # Count unique patients per gene (not just total mutations)
            patient_per_gene = (
                maf_df.assign(
                    patient=maf_df['Tumor_Sample_Barcode'].str[:12]
                )
                .groupby('Hugo_Symbol')['patient']
                .nunique()
                .sort_values(ascending=False)
            )
            
            # Show top genes with patient counts
            top_genes = patient_per_gene.head(30)
            n_mut_patients = maf_df['Tumor_Sample_Barcode'].str[:12].nunique()
            
            # Gene selection
            gene_options = top_genes.index.tolist()
            default_genes = []
            for g in ['TP53', 'PIK3CA', 'CDH1', 'GATA3', 'KRAS', 'BRAF', 'EGFR', 'PTEN']:
                if g in gene_options:
                    default_genes.append(g)
                if len(default_genes) >= 3:
                    break
            
            selected_genes = st.multiselect(
                f"Select genes to annotate ({proj})",
                gene_options,
                default=default_genes,
                key=f'mut_annot_genes_{proj}',
                help="Each selected gene adds a column to the expression metadata "
                     "(e.g., TP53_status = Mutant or Wild-type)",
            )
            
            if selected_genes:
                # Show preview: how many patients are mutant for each gene
                preview_rows = []
                for gene in selected_genes:
                    n_cases = patient_per_gene.get(gene, 0)
                    pct = 100 * n_cases / n_mut_patients if n_mut_patients > 0 else 0
                    preview_rows.append({
                        'Gene': gene,
                        'Mutant Cases': n_cases,
                        'Total MAF Cases': n_mut_patients,
                        'Mutation %': f"{pct:.1f}%",
                        'Total Mutations': gene_mut_counts.get(gene, 0),
                    })
                st.dataframe(pd.DataFrame(preview_rows), use_container_width=True, hide_index=True)
                
                st.caption(
                    f"Note: MAF covers {n_mut_patients} cases (open-access). "
                    f"Expression has {expr_ds.n_samples} samples. "
                    f"Samples not in MAF will be labeled 'Unknown (not in MAF)'."
                )
            
            if selected_genes and st.button(
                f"Annotate {expr_name} with Mutation Status",
                type="primary", key=f'mut_annot_btn_{proj}',
                use_container_width=True,
            ):
                with st.spinner("Cross-referencing mutations with expression samples..."):
                    try:
                        meta = expr_ds.metadata.copy()
                        
                        # Extract patient barcodes from expression samples
                        # Try patient_barcode column first, then parse from sample_id
                        if 'patient_barcode' in meta.columns:
                            expr_patients = meta['patient_barcode'].astype(str)
                        else:
                            # Parse from index (sample_id like TCGA-BH-A18H-01A-...)
                            expr_patients = pd.Series(
                                [str(sid)[:12] if str(sid).startswith('TCGA-') else str(sid)
                                 for sid in meta.index],
                                index=meta.index,
                            )
                        
                        # Extract patient barcodes from MAF
                        maf_patients = maf_df['Tumor_Sample_Barcode'].str[:12]
                        maf_patient_set = set(maf_patients.unique())
                        
                        n_annotated = 0
                        for gene in selected_genes:
                            # Find patients with mutations in this gene
                            gene_maf = maf_df[maf_df['Hugo_Symbol'] == gene]
                            mutant_patients = set(gene_maf['Tumor_Sample_Barcode'].str[:12].unique())
                            
                            # Classify each expression sample
                            col_name = f"{gene}_status"
                            statuses = []
                            for sid in meta.index:
                                patient = expr_patients.get(sid, str(sid)[:12])
                                if patient in mutant_patients:
                                    statuses.append('Mutant')
                                elif patient in maf_patient_set:
                                    statuses.append('Wild-type')
                                else:
                                    statuses.append('Unknown (not in MAF)')
                            
                            meta[col_name] = statuses
                            
                            # Count results
                            vc = meta[col_name].value_counts()
                            n_mut = vc.get('Mutant', 0)
                            n_wt = vc.get('Wild-type', 0)
                            n_unk = vc.get('Unknown (not in MAF)', 0)
                            n_annotated += 1
                            
                            st.info(
                                f"**{gene}:** {n_mut} Mutant, {n_wt} Wild-type, "
                                f"{n_unk} Unknown (not in MAF)"
                            )
                        
                        # Update the dataset's metadata in session state
                        expr_ds.metadata = meta
                        st.session_state['acq_datasets'][expr_name] = expr_ds
                        
                        st.success(
                            f"Added {n_annotated} mutation status column(s) to {expr_name}. "
                            f"Go to **Quality Check** → Color by → select `{selected_genes[0]}_status` "
                            f"to see mutation-driven separation in PCA."
                        )
                        
                        # Show the new columns
                        with st.expander("Preview annotated metadata"):
                            show_cols = ['patient_barcode'] if 'patient_barcode' in meta.columns else []
                            show_cols += [f"{g}_status" for g in selected_genes]
                            show_cols += ['sample_type', 'gender', 'tumor_stage']
                            show_cols = [c for c in show_cols if c in meta.columns]
                            st.dataframe(meta[show_cols].head(30), use_container_width=True)
                        
                    except Exception as e:
                        st.error(f"Annotation failed: {e}")
                        import traceback
                        st.code(traceback.format_exc())
    
    # Cached datasets
    st.markdown("---")
    st.markdown("### Cached on Disk")
    
    cached_list = cache.list_cached_datasets()
    
    # Also scan TCGA parquet_cache for non-expression data (mutations, miRNA, CNV, etc.)
    parquet_cache_dir = Path.home() / '.raptor' / 'parquet_cache'
    _PARQUET_DATA_TYPE_MAP = {
        '_mutations': 'somatic_mutation',
        '_mirna': 'mirna_rpm',
        '_cnv_gene_level': 'cnv_gene_level',
        '_cnv_segment': 'cnv_segment',
        '_cnv_masked': 'cnv_masked',
        '_methylation': 'methylation_beta',
        '_rppa': 'protein_expression_rppa',
    }
    
    if parquet_cache_dir.exists():
        seen_accessions = {d['accession'] for d in cached_list}
        
        for fpath in sorted(parquet_cache_dir.iterdir()):
            fname = fpath.stem  # e.g., "TCGA-BRCA_mutations" or "TCGA-BRCA_unstranded"
            # Strip .csv from .csv.gz files
            if fname.endswith('.csv'):
                fname = fname[:-4]
            
            # Match known non-expression data types
            for suffix, dtype in _PARQUET_DATA_TYPE_MAP.items():
                if suffix in fname:
                    accession_part = fname.split(suffix)[0]  # e.g., "TCGA-BRCA"
                    cache_key = f"{accession_part}{suffix}"
                    
                    # Skip if already listed via CacheManager
                    if cache_key in seen_accessions:
                        break
                    
                    # Get file size for sample/gene count estimate
                    try:
                        file_size_mb = fpath.stat().st_size / (1024 * 1024)
                        from datetime import datetime as _dt
                        cached_date = _dt.fromtimestamp(fpath.stat().st_mtime).isoformat()
                    except Exception:
                        file_size_mb = 0
                        cached_date = 'unknown'
                    
                    # Try to read shape without loading full data
                    n_genes = 0
                    n_samples = 0
                    try:
                        if fpath.suffix == '.parquet':
                            _peek = pd.read_parquet(fpath, engine='pyarrow', columns=[])
                            n_genes = len(_peek)
                            # Can't easily get columns without loading, estimate from file size
                        elif str(fpath).endswith('.csv.gz'):
                            import gzip as _gz
                            with _gz.open(fpath, 'rt') as _f:
                                header = _f.readline()
                                n_samples = header.count(',')
                                # Count a few lines for row estimate
                                for i, _ in enumerate(_f):
                                    if i > 100:
                                        break
                                n_genes = i + 1
                    except Exception:
                        pass
                    
                    cached_list.append({
                        'repository': 'TCGA',
                        'accession': cache_key,
                        'organism': 'Homo sapiens',
                        'n_genes': n_genes,
                        'n_samples': n_samples,
                        'data_type': dtype,
                        'cached_date': cached_date,
                        '_parquet_cache_path': str(fpath),  # internal: for loading
                    })
                    seen_accessions.add(cache_key)
                    break
    
    # Initialize project tags in session state
    if 'project_tags' not in st.session_state:
        st.session_state['project_tags'] = {}
    
    if cached_list:
        # Add project tags to display
        cached_df = pd.DataFrame(cached_list)
        cached_df['project'] = cached_df['accession'].map(
            lambda a: st.session_state['project_tags'].get(a, '')
        )
        
        # Fix misleading columns for SRA run tables and TCGA non-expression data
        # SRA entries have data_type='run_metadata' and n_genes is actually run count
        _DATA_TYPE_LABELS = {
            'somatic_mutation': 'Somatic Mutations (MAF)',
            'mirna_rpm': 'miRNA Expression',
            'cnv_gene_level': 'CNV (gene-level)',
            'cnv_segment': 'CNV (segments)',
            'cnv_masked': 'CNV (masked segments)',
            'methylation_beta': 'DNA Methylation',
            'protein_expression_rppa': 'Protein Expression (RPPA)',
            'run_metadata': 'Run metadata',
            'raw_counts': 'Gene Expression',
        }
        
        display_rows = []
        for _, row in cached_df.iterrows():
            r = dict(row)
            dtype = r.get('data_type', '')
            display_type = _DATA_TYPE_LABELS.get(dtype, dtype)
            
            if r.get('data_type') == 'run_metadata' or r.get('repository') == 'SRA':
                display_rows.append({
                    'project': r.get('project', ''),
                    'repository': r.get('repository', 'SRA'),
                    'accession': r.get('accession', ''),
                    'organism': r.get('organism', ''),
                    'runs': r.get('n_genes', ''),
                    'data_type': 'Run metadata',
                })
            else:
                display_rows.append({
                    'project': r.get('project', ''),
                    'repository': r.get('repository', ''),
                    'accession': r.get('accession', ''),
                    'organism': r.get('organism', ''),
                    'genes': f"{r.get('n_genes', 0):,}" if r.get('n_genes') else '',
                    'samples': r.get('n_samples', '') if r.get('n_samples') else '',
                    'data_type': display_type,
                })
        
        display_df = pd.DataFrame(display_rows)
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        
        # Project tagging
        with st.expander("Tag datasets with project names"):
            st.caption("Assign project names to remember what each dataset is for")
            cached_accessions = [d['accession'] for d in cached_list]
            
            col_tag_acc, col_tag_name, col_tag_btn = st.columns([2, 2, 1])
            with col_tag_acc:
                tag_acc = st.selectbox(
                    "Dataset",
                    cached_accessions,
                    key='tag_acc_select'
                )
            with col_tag_name:
                current_tag = st.session_state['project_tags'].get(tag_acc, '')
                tag_name = st.text_input(
                    "Project name",
                    value=current_tag,
                    placeholder="e.g., PS19, liver_cancer, my_thesis",
                    key='tag_name_input'
                )
            with col_tag_btn:
                st.markdown("<br>", unsafe_allow_html=True)
                if st.button("Set Tag", use_container_width=True, key='set_tag_btn'):
                    st.session_state['project_tags'][tag_acc] = tag_name
                    st.rerun()
            
            # Quick batch tag — apply same tag to multiple
            batch_tag = st.text_input(
                "Or tag all selected datasets below with:",
                placeholder="e.g., PS19_study",
                key='batch_tag_input'
            )
        
        # Multi-select for loading
        st.markdown("**Load datasets into session:**")
        cached_accessions = [d['accession'] for d in cached_list]
        
        # Show accessions with project tags for easier identification
        display_labels = []
        for acc in cached_accessions:
            tag = st.session_state['project_tags'].get(acc, '')
            label = f"{acc} [{tag}]" if tag else acc
            display_labels.append(label)
        
        selected_to_load = st.multiselect(
            "Select datasets to load",
            cached_accessions,
            format_func=lambda a: f"{a} [{st.session_state['project_tags'].get(a, '')}]"
                if st.session_state['project_tags'].get(a) else a,
            key='cache_multi_load',
        )
        
        col_load_sel, col_load_all, col_del = st.columns(3)
        
        def _load_cached_entry(entry):
            """Load a dataset from either CacheManager or parquet_cache."""
            acc = entry['accession']
            
            # Check if this is a parquet_cache entry
            pq_path = entry.get('_parquet_cache_path')
            if pq_path:
                pq_path = Path(pq_path)
                if not pq_path.exists():
                    return None
                try:
                    dtype = entry.get('data_type', '')
                    if str(pq_path).endswith('.csv.gz'):
                        df = pd.read_csv(pq_path, compression='gzip', index_col=0)
                    elif str(pq_path).endswith('.parquet'):
                        df = pd.read_parquet(pq_path, engine='pyarrow')
                    else:
                        return None
                    
                    # Build metadata for the loaded data
                    # Extract project ID from accession (e.g., "TCGA-BRCA_mutations" → "TCGA-BRCA")
                    proj_id = acc
                    for suffix in _PARQUET_DATA_TYPE_MAP:
                        if suffix in proj_id:
                            proj_id = proj_id.split(suffix)[0]
                            break
                    
                    # Build appropriate metadata based on data type
                    meta = pd.DataFrame()
                    gene_id_type = 'symbol'
                    
                    if dtype == 'somatic_mutation' and 'Tumor_Sample_Barcode' in df.columns:
                        barcodes = df['Tumor_Sample_Barcode'].unique()
                        meta = pd.DataFrame({
                            'sample_id': barcodes,
                            'case_id': [b[:12] if len(b) >= 12 else b for b in barcodes],
                            'mutations': [(df['Tumor_Sample_Barcode'] == b).sum() for b in barcodes],
                        }).set_index('sample_id')
                    elif dtype == 'mirna_rpm':
                        gene_id_type = 'mirna'
                    elif dtype == 'protein_expression_rppa':
                        gene_id_type = 'protein'
                    elif dtype.startswith('cnv'):
                        gene_id_type = 'ensembl' if 'gene_level' in dtype else 'segment'
                    elif dtype == 'methylation_beta':
                        gene_id_type = 'probe'
                    
                    # Also try loading companion metadata file
                    meta_path = pq_path.parent / f"{pq_path.stem.replace('.csv', '')}_meta.csv.gz"
                    if not meta_path.exists():
                        meta_path = pq_path.parent / f"{pq_path.stem}_meta.csv.gz"
                    if meta_path.exists():
                        try:
                            meta = pd.read_csv(meta_path, compression='gzip', index_col=0)
                        except Exception:
                            pass
                    
                    return AcquiredDataset(
                        counts_df=df,
                        metadata=meta if not meta.empty else pd.DataFrame(index=df.columns[:0]),
                        source_info={
                            'repository': 'TCGA',
                            'accession': proj_id,
                            'data_type': dtype,
                            'organism': 'Homo sapiens',
                        },
                        gene_id_type=gene_id_type,
                    )
                except Exception as e:
                    logger.warning(f"Failed to load {acc} from parquet_cache: {e}")
                    return None
            
            # Standard CacheManager load
            repo = entry.get('repository', 'unknown')
            data = cache.load_dataset(repo, acc)
            if data:
                return AcquiredDataset(
                    counts_df=data['counts_df'],
                    metadata=data['metadata'],
                    source_info=data['source_info'],
                    gene_id_type=data['gene_id_type'],
                )
            return None
        
        with col_load_sel:
            if st.button("Load Selected", use_container_width=True, key='load_selected_btn'):
                if selected_to_load:
                    loaded = 0
                    for acc in selected_to_load:
                        entry = next((d for d in cached_list if d['accession'] == acc), None)
                        if entry is None:
                            continue
                        ds = _load_cached_entry(entry)
                        if ds:
                            st.session_state['acq_datasets'][acc] = ds
                            # Apply batch tag if set
                            if batch_tag:
                                st.session_state['project_tags'][acc] = batch_tag
                            loaded += 1
                    st.success(f"Loaded {loaded} dataset(s) into session")
                    st.rerun()
                else:
                    st.warning("Select at least one dataset")
        
        with col_load_all:
            if st.button("Load All", use_container_width=True, key='load_all_btn'):
                loaded = 0
                for entry in cached_list:
                    acc = entry['accession']
                    ds = _load_cached_entry(entry)
                    if ds:
                        st.session_state['acq_datasets'][acc] = ds
                        loaded += 1
                st.success(f"Loaded all {loaded} dataset(s) into session")
                st.rerun()
        
        with col_del:
            if selected_to_load and st.button("Delete Selected", use_container_width=True, key='delete_selected_btn'):
                deleted = 0
                for acc in selected_to_load:
                    entry = next((d for d in cached_list if d['accession'] == acc), None)
                    if entry is None:
                        continue
                    
                    # Check if this is a parquet_cache entry (miRNA, methylation, etc.)
                    pq_path = entry.get('_parquet_cache_path')
                    if pq_path:
                        try:
                            pq_file = Path(pq_path)
                            pq_dir = pq_file.parent
                            # Derive base name (handle .csv.gz double extension)
                            base = pq_file.stem
                            if base.endswith('.csv'):
                                base = base[:-4]
                            # Delete main file and any associated meta files
                            for f in pq_dir.glob(f"{base}*"):
                                f.unlink()
                            catalog.unregister_dataset(entry['repository'], acc)
                            if acc in st.session_state.get('project_tags', {}):
                                del st.session_state['project_tags'][acc]
                            deleted += 1
                        except Exception as e:
                            st.warning(f"Failed to delete {acc}: {e}")
                    else:
                        # Standard CacheManager entry
                        repo = entry['repository']
                        if cache.delete_dataset(repo, acc):
                            catalog.unregister_dataset(repo, acc)
                            if acc in st.session_state.get('project_tags', {}):
                                del st.session_state['project_tags'][acc]
                            deleted += 1
                st.success(f"Deleted {deleted} dataset(s) from cache")
                st.rerun()
    else:
        st.info("No datasets cached on disk yet.")
    
    # Cache info
    with st.expander("Cache Information"):
        st.code(cache.cache_summary())

# =============================================================================
# TAB 3: Pool Datasets
# =============================================================================

with tab_pool:
    st.markdown("## Pool Datasets")
    st.caption("Combine multiple studies into one powerful dataset — RAPTOR handles the technical details for you")
    
    with st.expander("Why pool datasets? (click to expand)"):
        st.markdown("""
        ### Cross-Cohort Validation for Biomarker Discovery
        
        Pooling is the key step that turns single-study findings into robust, publishable biomarker 
        candidates. Here's why:
        
        **The challenge with single studies:** A single RNA-seq experiment might identify hundreds 
        of differentially expressed genes, but many of these will be artifacts of your specific 
        cohort — driven by batch effects, sample demographics, sequencing platform differences, 
        or chance variation in a small sample size.
        
        **What pooling gives you:**
        
        - **Statistical power** — instead of 20–30 samples from one lab, you analyze hundreds of 
          samples from multiple independent cohorts
        - **Cross-cohort validation** — a gene that is consistently differentially expressed across 
          your data AND 2–3 independent public studies is a far more credible biomarker candidate
        - **Batch effect detection** — the Quality Check tab reveals whether technical differences 
          between labs dominate the signal (bad) or whether biological signal is consistent across 
          studies (good)
        - **Publication strength** — reviewers increasingly expect cross-cohort validation; 
          "we validated our findings across N independent cohorts" is a strong statement
        
        ### Recommended Workflow
        
        1. **Upload your own count matrix** from your RNA-seq experiment (Search tab → Upload)
        2. **Download 2–3 public studies** of the same disease/condition from GEO or TCGA
        3. **Pool all datasets here** — RAPTOR harmonizes gene IDs and applies batch correction
        4. **Check quality** — verify studies intermingle by condition (not by batch) in PCA
        5. **Run differential expression** on the pooled data for robust biomarker candidates
        6. **Biomarker Discovery** (Module 10) will rank candidates by cross-study consistency
        
        ### What Happens Under the Hood
        
        When you click "Pool Selected Datasets", RAPTOR:
        - Converts all gene IDs to a common type (e.g., gene symbols) via MyGene.info
        - Finds the intersection (or union) of genes across all studies
        - Merges count matrices and metadata into a single dataset
        - Applies batch correction (ComBat, quantile normalization, or median ratio) 
          to remove lab-specific technical variation while preserving biological signal
        """)
    
    session_ds = st.session_state.get('acq_datasets', {})
    
    # Filter out SRA run tables and mutation data — they can't be pooled
    poolable_ds = {}
    sra_skipped = []
    for name, ds in session_ds.items():
        is_sra = (
            ds.data_type == 'run_metadata'
            or (hasattr(ds, 'source_info') and ds.source_info.get('repository') == 'SRA')
        )
        is_mut = hasattr(ds, 'data_type') and 'mutation' in str(ds.data_type).lower()
        if is_sra:
            sra_skipped.append(name)
        elif is_mut:
            pass  # silently skip mutation data from pooling
        else:
            poolable_ds[name] = ds
    
    if sra_skipped:
        st.caption(
            f"SRA run tables excluded from pooling ({', '.join(sra_skipped)}) "
            f"— they contain run metadata, not gene expression counts. "
            f"Download the linked GSE from GEO for count data."
        )
    
    if len(poolable_ds) < 2:
        st.warning(
            f"You need at least 2 expression datasets to pool. "
            f"Currently have {len(poolable_ds)} "
            f"(plus {len(sra_skipped)} SRA run table(s) which cannot be pooled). "
            f"Go to the Search tab to download more."
        )
    else:
        # Dataset selection
        st.markdown("### 1. Select Datasets to Pool")
        
        available_names = list(poolable_ds.keys())
        selected_for_pool = st.multiselect(
            "Select datasets (minimum 2)",
            available_names,
            default=available_names[:2] if len(available_names) >= 2 else available_names,
            help="Select the datasets you want to merge"
        )
        
        if len(selected_for_pool) >= 2:
            # Show selection summary
            sel_rows = []
            for name in selected_for_pool:
                ds = poolable_ds[name]
                sel_rows.append({
                    'Dataset': name,
                    'Genes': ds.n_genes,
                    'Samples': ds.n_samples,
                    'Gene IDs': ds.gene_id_type,
                    'Organism': ds.organism,
                })
            st.dataframe(pd.DataFrame(sel_rows), use_container_width=True, hide_index=True)
            
            # Gene overlap preview
            gene_sets = [set(poolable_ds[n].gene_ids) for n in selected_for_pool]
            common = gene_sets[0]
            all_genes = gene_sets[0]
            for gs in gene_sets[1:]:
                common = common & gs
                all_genes = all_genes | gs
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Common Genes", f"{len(common):,}")
            with col2:
                st.metric("Total Unique Genes", f"{len(all_genes):,}")
            with col3:
                overlap_pct = 100 * len(common) / len(all_genes) if all_genes else 0
                st.metric("Overlap %", f"{overlap_pct:.1f}%")
            
            # Pooling options
            st.markdown("### 2. Pooling Options")
            
            col_opt1, col_opt2, col_opt3 = st.columns(3)
            
            with col_opt1:
                join_method = st.selectbox(
                    "Join method",
                    ["inner", "outer"],
                    help="**inner** = keep only genes found in ALL datasets (recommended, safer). "
                         "**outer** = keep all genes, fill missing values with 0."
                )
            
            with col_opt2:
                batch_correction = st.selectbox(
                    "Batch correction",
                    ["none", "combat", "quantile", "median_ratio"],
                    help="Different labs produce slightly different measurements. "
                         "Batch correction adjusts for these differences. "
                         "**combat** is the gold standard. **none** is fine if studies are similar."
                )
            
            with col_opt3:
                pool_name = st.text_input(
                    "Pool name",
                    value=f"pool_{'_'.join(selected_for_pool[:2])}",
                    help="Give your combined dataset a memorable name"
                )
            
            # Smart default: if all datasets have same gene ID type, default to that (skip conversion)
            all_gene_types = [poolable_ds[n].gene_id_type for n in selected_for_pool]
            most_common_type = max(set(all_gene_types), key=all_gene_types.count) if all_gene_types else 'symbol'
            gene_type_options = ["symbol", "ensembl", "entrez"]
            default_idx = gene_type_options.index(most_common_type) if most_common_type in gene_type_options else 0
            
            target_gene_id = st.selectbox(
                "Target gene ID type",
                gene_type_options,
                index=default_idx,
                help="Harmonize all datasets to this ID type. If all datasets already use the same type, keep it to skip conversion (faster)."
            )
            
            # Run pooling
            st.markdown("### 3. Run Pooling")
            
            if st.button("Pool Selected Datasets", type="primary", use_container_width=True):
                datasets_to_pool = [poolable_ds[n] for n in selected_for_pool]
                
                # Check if gene ID conversion can be skipped
                all_same_type = len(set(d.gene_id_type for d in datasets_to_pool)) == 1
                source_type = datasets_to_pool[0].gene_id_type
                skip_conversion = all_same_type and source_type == target_gene_id
                
                if skip_conversion:
                    spinner_msg = f"Pooling {len(datasets_to_pool)} datasets... (same gene IDs — skipping conversion)"
                else:
                    spinner_msg = f"Pooling {len(datasets_to_pool)} datasets... (converting gene IDs to {target_gene_id}, merging)"
                
                with st.spinner(spinner_msg):
                    try:
                        engine = PoolingEngine(
                            target_gene_id=target_gene_id,
                            species=datasets_to_pool[0].organism,
                        )
                        
                        pool = engine.merge(
                            datasets_to_pool,
                            pool_name=pool_name,
                            method=join_method,
                            batch_correction=batch_correction,
                        )
                        
                        st.session_state['acq_pooled'] = pool
                        
                        # Cache the pool
                        cache.save_pool(
                            pool_name=pool_name,
                            counts_df=pool.counts_df,
                            metadata=pool.metadata,
                            study_labels=pool.study_labels,
                            pooling_info=pool.pooling_info,
                            component_datasets=pool.component_datasets,
                        )
                        
                        st.success(
                            f"Pooling complete: {pool.n_genes:,} genes, "
                            f"{pool.n_samples} samples from {pool.n_studies} studies"
                        )
                        
                    except Exception as e:
                        st.error(f"Pooling failed: {e}")
    
    # Show pooled result if available
    pool = st.session_state.get('acq_pooled')
    
    if pool is not None:
        st.markdown("---")
        st.markdown("### Pooled Dataset Summary")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Studies", pool.n_studies)
        with col2:
            st.metric("Genes", f"{pool.n_genes:,}")
        with col3:
            st.metric("Samples", pool.n_samples)
        with col4:
            st.metric("Batch Correction", pool.pooling_info.get('batch_correction', 'none'))
        
        # Samples per study chart
        sps = pool.samples_per_study
        fig = go.Figure(data=[
            go.Bar(x=list(sps.keys()), y=list(sps.values()), marker_color='#2E7D32')
        ])
        fig.update_layout(
            title="Samples per Study",
            xaxis_title="Study",
            yaxis_title="Number of Samples",
            height=350,
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Preview
        with st.expander("Preview pooled count matrix"):
            st.dataframe(pool.counts_df.head(20), use_container_width=True)

# =============================================================================
# TAB 4: Quality Check
# =============================================================================

with tab_qc:
    st.markdown("## Quality Check")
    st.caption(
        "Explore any dataset before analysis — PCA, library sizes, expression distributions, "
        "and more. Works on GEO, TCGA, ArrayExpress, uploaded data, or pooled datasets."
    )
    
    session_ds = st.session_state.get('acq_datasets', {})
    pool = st.session_state.get('acq_pooled')
    
    # Filter out SRA run tables
    qc_datasets = {
        name: ds for name, ds in session_ds.items()
        if ds.data_type != 'run_metadata'
        and not (hasattr(ds, 'source_info') and ds.source_info.get('repository') == 'SRA')
        and not (hasattr(ds, 'data_type') and 'mutation' in str(ds.data_type).lower())
    }
    
    if not qc_datasets and pool is None:
        st.info("Download a dataset first (Search tab), then come here to check quality.")
    else:
        # Choose what to QC
        qc_options = []
        if qc_datasets:
            qc_options.extend([f"Dataset: {name}" for name in qc_datasets.keys()])
        if pool is not None:
            qc_options.append("Pooled Data")
        
        qc_selection = st.selectbox("Select data to check", qc_options, key='qc_select')
        
        # ---- SINGLE DATASET QC ----
        if qc_selection and qc_selection.startswith("Dataset:"):
            ds_name = qc_selection.replace("Dataset: ", "")
            ds = qc_datasets[ds_name]
            counts = ds.counts_df
            metadata = ds.metadata
            
            st.success(f"**{ds_name}** — {ds.n_genes:,} {_feature_label(ds)}, {ds.n_samples} samples, {ds.organism}")
            
            # Choose color-by column from metadata
            color_options = ['None']
            if not metadata.empty:
                # Find columns with 2-20 unique values (good for coloring)
                for col in metadata.columns:
                    n_unique = metadata[col].nunique()
                    if 2 <= n_unique <= 20 and col not in ('case_id', 'patient_barcode', 'morphology'):
                        color_options.append(col)
            
            color_by = st.selectbox(
                "Color samples by",
                color_options,
                index=color_options.index('sample_type') if 'sample_type' in color_options
                    else (color_options.index('sample_category') if 'sample_category' in color_options else 0),
                key='qc_color_by',
                help="Choose a metadata column to group and color samples"
            )
            
            # Basic QC metrics (cached)
            st.markdown("### Overview")
            dhash = _data_hash(counts)
            lib_vals, det_vals, zero_pct = _cached_lib_stats(dhash, counts.values, gene_id_type=ds.gene_id_type)
            lib_sizes = pd.Series(lib_vals, index=counts.columns)
            gene_detected = pd.Series(det_vals, index=counts.columns)
            
            is_beta_qc = ds.gene_id_type == 'probe'
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                if is_beta_qc:
                    st.metric("Median Beta", f"{lib_sizes.median():.3f}")
                else:
                    st.metric("Median Library Size", f"{lib_sizes.median():,.0f}")
            with col2:
                cv_val = lib_sizes.std() / lib_sizes.mean() if lib_sizes.mean() != 0 else 0
                st.metric("Beta CV" if is_beta_qc else "Library Size CV", f"{cv_val:.2f}")
            with col3:
                if is_beta_qc:
                    st.metric("Probes with Signal", f"{gene_detected.median():,.0f}")
                else:
                    st.metric(f"Median {_feature_label(ds).capitalize()} Detected", f"{gene_detected.median():,.0f}")
            with col4:
                st.metric("NaN %" if is_beta_qc else "Zero %", f"{zero_pct:.1f}%")
            
            # Build color labels
            if color_by != 'None' and not metadata.empty and color_by in metadata.columns:
                # Match sample order
                color_labels = metadata.reindex(counts.columns)[color_by].fillna('unknown')
            else:
                color_labels = pd.Series(['all'] * counts.shape[1], index=counts.columns)
            
            # Option to exclude "Unknown (not in MAF)" when viewing mutation status
            is_mutation_color = color_by.endswith('_status') and 'Unknown (not in MAF)' in color_labels.values
            if is_mutation_color:
                n_unknown = (color_labels == 'Unknown (not in MAF)').sum()
                exclude_unknown = st.checkbox(
                    f"Exclude {n_unknown} 'Unknown (not in MAF)' samples from QC",
                    value=False,
                    key='qc_exclude_unknown',
                    help="These samples have no open-access mutation data in GDC. "
                         "Excluding them shows a cleaner Mutant vs Wild-type comparison.",
                )
                if exclude_unknown:
                    keep_mask = color_labels != 'Unknown (not in MAF)'
                    keep_samples = color_labels[keep_mask].index.tolist()
                    counts = counts[keep_samples]
                    color_labels = color_labels[keep_mask]
                    metadata = metadata.loc[metadata.index.isin(keep_samples)]
                    st.caption(
                        f"Showing {len(keep_samples)} samples "
                        f"({(color_labels == 'Mutant').sum()} Mutant, "
                        f"{(color_labels == 'Wild-type').sum()} Wild-type)"
                    )
                    # Recompute metrics for filtered data
                    dhash = _data_hash(counts)
                    lib_vals, det_vals, zero_pct = _cached_lib_stats(dhash, counts.values, gene_id_type=ds.gene_id_type)
                    lib_sizes = pd.Series(lib_vals, index=counts.columns)
                    gene_detected = pd.Series(det_vals, index=counts.columns)
            
            # QC sub-tabs
            qc_t1, qc_t2, qc_t3, qc_t4, qc_t5, qc_t6 = st.tabs([
                "Library Sizes", "PCA", "Gene Expression",
                "Heatmap", "Sample Correlation", "QC Metrics"
            ])
            
            with qc_t1:
                lib_df = pd.DataFrame({
                    'Sample': lib_sizes.index,
                    'Library Size': lib_sizes.values,
                    'Group': color_labels.values,
                })
                
                fig = px.box(
                    lib_df, x='Group', y='Library Size',
                    color='Group', points='all',
                    title="Library Size Distribution",
                )
                fig.update_layout(height=450, showlegend=False)
                _show_plot(fig, "RAPTOR_qc_1", "1")
                
                # Genes detected per sample
                det_df = pd.DataFrame({
                    'Sample': gene_detected.index,
                    'Genes Detected': gene_detected.values,
                    'Group': color_labels.values,
                })
                fig2 = px.box(
                    det_df, x='Group', y='Genes Detected',
                    color='Group', points='all',
                    title="Genes Detected per Sample (count > 0)",
                )
                fig2.update_layout(height=400, showlegend=False)
                _show_plot(fig2, "RAPTOR_qc_2", "2")
                
                # Assess
                if color_labels.nunique() > 1:
                    group_medians = lib_df.groupby('Group')['Library Size'].median()
                    fold = group_medians.max() / group_medians.min() if group_medians.min() > 0 else float('inf')
                    if fold > 5:
                        st.warning(f"Library size varies {fold:.1f}-fold across groups. Consider normalization.")
                    elif fold > 2:
                        st.info(f"Library size varies {fold:.1f}-fold across groups. Moderate variation.")
                    else:
                        st.success(f"Library sizes are consistent ({fold:.1f}-fold range).")
            
            with qc_t2:
                st.markdown("**PCA of samples**")
                
                colors_palette = px.colors.qualitative.Set2
                
                # Cached PCA on top 5K variable features
                pca_result, var_ratios, n_comp = _cached_pca(
                    dhash, counts.values, counts.shape[1], n_comp=6,
                    gene_id_type=ds.gene_id_type,
                )
                
                pc_cols = [f'PC{i+1}' for i in range(n_comp)]
                var_pct = {f'PC{i+1}': var_ratios[i] * 100 for i in range(n_comp)}
                var_labels = {pc: f'{pc} ({var_pct[pc]:.1f}%)' for pc in pc_cols}
                
                pca_df = pd.DataFrame(pca_result, columns=pc_cols, index=counts.columns)
                pca_df['Sample'] = pca_df.index
                
                # Add all metadata columns for coloring
                pca_color_options = ['None']
                if not metadata.empty:
                    for col in metadata.columns:
                        if metadata[col].nunique() <= 20 and col not in ('case_id', 'patient_barcode', 'morphology'):
                            pca_df[col] = metadata[col].reindex(pca_df.index)
                            pca_color_options.append(col)
                
                # Controls
                ctrl1, ctrl2, ctrl3, ctrl4 = st.columns(4)
                with ctrl1:
                    pca_color = st.selectbox(
                        "Color by", pca_color_options,
                        index=pca_color_options.index('sample_type') if 'sample_type' in pca_color_options
                            else min(1, len(pca_color_options) - 1),
                        key="qc_pca_color"
                    )
                with ctrl2:
                    pca_x = st.selectbox("X axis", pc_cols, index=0, key="qc_pca_x")
                with ctrl3:
                    pca_y = st.selectbox("Y axis", pc_cols, index=min(1, n_comp - 1), key="qc_pca_y")
                with ctrl4:
                    pca_dim = st.selectbox("Dimension", ["2D", "3D"] if n_comp >= 3 else ["2D"], key="qc_pca_dim")
                
                if pca_dim == "3D":
                    pca_z = st.selectbox("Z axis", pc_cols, index=min(2, n_comp - 1), key="qc_pca_z")
                
                show_labels = st.checkbox("Show sample labels", value=len(counts.columns) <= 30, key="qc_pca_labels")
                
                color_arg = pca_color if pca_color != 'None' else None
                
                if pca_dim == "2D":
                    fig = px.scatter(
                        pca_df, x=pca_x, y=pca_y, color=color_arg,
                        color_discrete_sequence=colors_palette,
                        text=pca_df.index if show_labels else None,
                        hover_data=['Sample'],
                    )
                    fig.update_traces(
                        marker=dict(size=12, line=dict(width=1.5, color='white')),
                        textposition='top center', textfont=dict(size=9),
                    )
                    fig.update_layout(
                        title=f"PCA — {var_labels[pca_x]} vs {var_labels[pca_y]}",
                        xaxis_title=var_labels[pca_x],
                        yaxis_title=var_labels[pca_y],
                        height=650, plot_bgcolor='#fafafa', font=dict(size=12),
                    )
                    fig.update_xaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
                    fig.update_yaxes(showgrid=True, gridcolor='#eee', zeroline=True, zerolinecolor='#ccc')
                else:
                    fig = px.scatter_3d(
                        pca_df, x=pca_x, y=pca_y, z=pca_z, color=color_arg,
                        color_discrete_sequence=colors_palette,
                        text=pca_df.index if show_labels else None,
                        hover_data=['Sample'],
                    )
                    fig.update_traces(marker=dict(size=8, line=dict(width=0.5, color='white')))
                    fig.update_layout(
                        title=dict(text=f"<b>3D PCA — {var_labels[pca_x]} vs {var_labels[pca_y]} vs {var_labels[pca_z]}</b>",
                                   font=dict(size=16, family="Arial")),
                        height=750, font=PUB_FONT,
                        paper_bgcolor='white',
                        scene=dict(
                            xaxis_title=var_labels[pca_x],
                            yaxis_title=var_labels[pca_y],
                            zaxis_title=var_labels[pca_z],
                        ),
                        legend=dict(font=dict(size=13), bordercolor="gray", borderwidth=1),
                    )
                
                # Apply publication layout for 2D plots
                if pca_dim == "2D":
                    _pub_layout(fig, f"PCA — {var_labels[pca_x]} vs {var_labels[pca_y]}", width=1000, height=700)
                
                _show_plot(fig, "RAPTOR_PCA", "3")
                
                # Direct download buttons for publication
                dl1, dl2 = st.columns(2)
                with dl1:
                    try:
                        png_bytes = fig.to_image(format="png", width=1200, height=900, scale=3, engine="kaleido")
                        st.download_button(
                            "Download PNG (300 DPI)", data=png_bytes,
                            file_name="RAPTOR_PCA.png", mime="image/png",
                            key="dl_pca_png",
                        )
                    except Exception:
                        st.caption("Install kaleido for PNG export: `pip install kaleido`")
                with dl2:
                    try:
                        svg_bytes = fig.to_image(format="svg", width=1200, height=900, scale=3, engine="kaleido")
                        st.download_button(
                            "Download SVG (vector)", data=svg_bytes,
                            file_name="RAPTOR_PCA.svg", mime="image/svg+xml",
                            key="dl_pca_svg",
                        )
                    except Exception:
                        st.caption("Install kaleido for SVG export: `pip install kaleido`")
                
                if color_arg and color_labels.nunique() > 1:
                    st.caption(
                        "If groups form distinct clusters, the biological signal is strong. "
                        "For tumor vs. normal: clear separation = differential expression exists."
                    )
                
                # Scree plot
                with st.expander("Variance explained (Scree plot)"):
                    scree = pd.DataFrame({
                        'Component': pc_cols,
                        'Variance (%)': [var_pct[pc] for pc in pc_cols],
                        'Cumulative (%)': list(np.cumsum([var_pct[pc] for pc in pc_cols])),
                    })
                    fig_s = go.Figure()
                    fig_s.add_bar(
                        x=scree['Component'], y=scree['Variance (%)'],
                        name='Individual', marker_color='#66c2a5',
                    )
                    fig_s.add_scatter(
                        x=scree['Component'], y=scree['Cumulative (%)'],
                        name='Cumulative', mode='lines+markers', marker_color='#fc8d62',
                    )
                    fig_s.update_layout(
                        yaxis_title='Variance (%)', height=300,
                        plot_bgcolor='#fafafa', title="Scree Plot",
                    )
                    _show_plot(fig_s, "RAPTOR_qc_4", "4")
                    st.dataframe(scree, use_container_width=True, hide_index=True)
                
                # Export PCA coordinates
                with st.expander("Export PCA data"):
                    export_cols = pc_cols + (['Sample'] + [c for c in pca_df.columns if c not in pc_cols and c != 'Sample'])
                    export_pca = pca_df[[c for c in export_cols if c in pca_df.columns]]
                    st.dataframe(export_pca.head(10), use_container_width=True)
                    csv_pca = export_pca.to_csv()
                    st.download_button(
                        "Download PCA coordinates (CSV)",
                        data=csv_pca,
                        file_name="RAPTOR_PCA_coordinates.csv",
                        mime="text/csv",
                    )
                
                st.caption(
                    "All plots are interactive — zoom, pan, hover for details. "
                    "Click the camera icon (top right of plot) to save as SVG/PNG."
                )
            
            with qc_t3:
                st.markdown("**Expression distributions**")
                sample_genes = counts.index[:2000] if len(counts) > 2000 else counts.index
                
                fig = go.Figure()
                for group in color_labels.unique():
                    group_samples = [s for s, g in zip(counts.columns, color_labels) if g == group]
                    if group_samples:
                        mean_expr = np.log2(counts.loc[sample_genes, group_samples].mean(axis=1) + 1)
                        fig.add_trace(go.Violin(
                            y=mean_expr, name=str(group),
                            box_visible=True, meanline_visible=True,
                        ))
                
                fig.update_layout(
                    title="Log2 Mean Expression Distribution",
                    yaxis_title="Log2(mean count + 1)",
                    height=450,
                )
                _show_plot(fig, "RAPTOR_qc_5", "5")
                
                # Per-sample density
                st.markdown("**Per-sample expression density**")
                fig_density = go.Figure()
                # Show up to 20 samples
                sample_subset = list(counts.columns[:20])
                for s in sample_subset:
                    vals = np.log2(counts.loc[sample_genes, s] + 1)
                    vals = vals[vals > 0]
                    fig_density.add_trace(go.Violin(
                        y=vals, name=s[:15],
                        box_visible=False, meanline_visible=True,
                        line=dict(width=1),
                    ))
                fig_density.update_layout(
                    title=f"Per-Sample Expression Density (first {len(sample_subset)} samples)",
                    yaxis_title="Log2(count + 1)",
                    height=450, showlegend=False,
                )
                _show_plot(fig_density, "RAPTOR_qc_6", "6")
            
            with qc_t4:
                from scipy.cluster.hierarchy import linkage, leaves_list
                from scipy.spatial.distance import pdist
                
                feat = _feature_label(ds)
                st.markdown(f"**Heatmap — Top Variable {feat.capitalize()}**")
                hm1, hm2, hm3 = st.columns(3)
                with hm1:
                    n_genes_hm = st.slider(f"Top variable {feat}", 20, 200, 50, 10, key="qc_hm_ngenes")
                with hm2:
                    hm_cs = st.selectbox("Color scale", ["RdBu_r", "Viridis", "Plasma", "PiYG", "BrBG"], key="qc_hm_cs")
                with hm3:
                    hm_cluster = st.checkbox("Cluster rows and columns", value=True, key="qc_hm_cluster")
                
                is_count_hm = ds.gene_id_type in ('ensembl', 'symbol', 'entrez', 'refseq', 'unknown')
                hm_data = counts.copy()

                # For data with NaN (methylation): drop features with >50% missing before variance calc
                if hm_data.isna().any().any():
                    nan_frac = hm_data.isna().sum(axis=1) / hm_data.shape[1]
                    hm_data = hm_data.loc[nan_frac < 0.5]
                    # Impute remaining NaN with row median for variance/display
                    row_med = hm_data.median(axis=1)
                    hm_data = hm_data.apply(lambda col: col.fillna(row_med))

                if is_count_hm:
                    log_c_hm = np.log2(hm_data + 1)
                else:
                    log_c_hm = hm_data

                gene_var_hm = log_c_hm.var(axis=1).sort_values(ascending=False)
                top_genes_hm = gene_var_hm.head(n_genes_hm).index
                subset = log_c_hm.loc[top_genes_hm]
                
                # Z-score normalization (row-wise)
                row_m = subset.mean(axis=1)
                row_s = subset.std(axis=1).replace(0, 1)
                subset_z = subset.subtract(row_m, axis=0).divide(row_s, axis=0)
                
                # Hierarchical clustering
                if hm_cluster and len(top_genes_hm) > 2 and subset_z.shape[1] > 2:
                    try:
                        r_dist = pdist(subset_z.values, 'correlation')
                        c_dist = pdist(subset_z.T.values, 'correlation')
                        r_dist = np.nan_to_num(r_dist, nan=0.0)
                        c_dist = np.nan_to_num(c_dist, nan=0.0)
                        r_link = linkage(r_dist, method='average')
                        c_link = linkage(c_dist, method='average')
                        subset_z = subset_z.iloc[leaves_list(r_link), leaves_list(c_link)]
                    except Exception:
                        pass
                
                fig = go.Figure(go.Heatmap(
                    z=subset_z.values,
                    x=subset_z.columns.tolist(),
                    y=subset_z.index.tolist(),
                    colorscale=hm_cs, zmid=0,
                    hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>',
                    colorbar=dict(title=dict(text="Z-score"), len=0.7),
                ))
                fig.update_layout(
                    title=f"Top {n_genes_hm} Variable {feat.capitalize()} (Z-score normalized)",
                    height=max(500, n_genes_hm * 12),
                    font=dict(size=11), plot_bgcolor='white',
                    xaxis=dict(tickangle=45, tickfont=dict(size=7)),
                    yaxis=dict(tickfont=dict(size=max(6, min(10, 500 // n_genes_hm)))),
                )
                _show_plot(fig, "RAPTOR_qc_7", "7")
                
                # Annotation bars from metadata
                if not metadata.empty:
                    with st.expander("Sample annotations"):
                        anno_cols = [c for c in metadata.columns if metadata[c].nunique() <= 10
                                     and c not in ('case_id', 'patient_barcode', 'morphology')]
                        if anno_cols:
                            colors_palette = px.colors.qualitative.Set2
                            sel_anno = st.multiselect("Show:", anno_cols,
                                default=[c for c in ['sample_type', 'sample_category', 'gender'] if c in anno_cols][:2],
                                key='qc_hm_anno')
                            for ac in sel_anno:
                                vals = metadata[ac].reindex(subset_z.columns)
                                uv = sorted(vals.dropna().unique())
                                nv = max(len(uv) - 1, 1)
                                v2n = {v: i for i, v in enumerate(uv)}
                                fig_a = go.Figure(go.Heatmap(
                                    z=[[v2n.get(v, -1) for v in vals]],
                                    x=subset_z.columns.tolist(), y=[ac],
                                    colorscale=[[i/nv, colors_palette[i % len(colors_palette)]] for i in range(len(uv))],
                                    showscale=False,
                                    hovertemplate='%{x}: %{text}<extra></extra>',
                                    text=[[str(v) for v in vals]],
                                ))
                                fig_a.update_layout(height=80, margin=dict(l=100, r=20, t=5, b=5),
                                                   xaxis=dict(showticklabels=False))
                                _show_plot(fig_a, "RAPTOR_qc_8", "8")
                                st.caption(f"{ac}: " + ", ".join(str(v) for v in uv))
                
                # Feature search
                with st.expander(f"Search specific {feat}"):
                    gs = st.text_input(f"{_feature_label(ds, plural=False).capitalize()} IDs (comma-separated):", key="qc_hm_search")
                    if gs:
                        sg = [g.strip() for g in gs.split(',')]
                        found = [g for g in sg if g in subset_z.index]
                        if found:
                            st.success(f"Found: {', '.join(found)}")
                            st.dataframe(subset_z.loc[found], use_container_width=True)
                        nf = [g for g in sg if g not in subset_z.index]
                        if nf:
                            st.warning(f"Not in top {n_genes_hm}: {', '.join(nf)}")
                
                # Table
                with st.expander(f"Top variable {feat} table"):
                    feat_singular = _feature_label(ds, plural=False).capitalize()
                    val_label = "Mean Beta" if ds.gene_id_type == 'probe' else "Mean Expression"
                    var_df = pd.DataFrame({
                        feat_singular: gene_var_hm.head(n_genes_hm).index,
                        'Variance': gene_var_hm.head(n_genes_hm).values.round(2),
                        val_label: log_c_hm.loc[gene_var_hm.head(n_genes_hm).index].mean(axis=1).round(2).values,
                    })
                    st.dataframe(var_df, use_container_width=True, hide_index=True)
            
            with qc_t5:
                st.markdown("**Sample-to-sample correlation**")
                
                # Guard: subsample large datasets to prevent browser freeze
                MAX_SAMPLES_CORR = 150
                corr_counts = counts
                corr_color_labels = color_labels if color_by != 'None' else None
                
                if counts.shape[1] > MAX_SAMPLES_CORR:
                    st.warning(
                        f"Dataset has {counts.shape[1]} samples. Subsampling to {MAX_SAMPLES_CORR} "
                        f"for correlation heatmap (a {counts.shape[1]}x{counts.shape[1]} matrix "
                        f"would freeze the browser)."
                    )
                    # Stratified subsampling: proportional from each group
                    if color_by != 'None' and color_labels is not None and color_labels.nunique() > 1:
                        sub_idx = []
                        for grp in color_labels.unique():
                            grp_samples = [s for s in counts.columns if color_labels.get(s, '') == grp]
                            n_take = max(1, int(MAX_SAMPLES_CORR * len(grp_samples) / counts.shape[1]))
                            sub_idx.extend(grp_samples[:n_take])
                        sub_idx = sub_idx[:MAX_SAMPLES_CORR]
                    else:
                        # Uniform subsampling
                        step = max(1, counts.shape[1] // MAX_SAMPLES_CORR)
                        sub_idx = list(counts.columns[::step])[:MAX_SAMPLES_CORR]
                    corr_counts = counts[sub_idx]
                
                cc1, cc2, cc3 = st.columns(3)
                with cc1:
                    corr_method = st.selectbox("Method", ["pearson", "spearman"], key="qc_corr_method")
                with cc2:
                    corr_cs = st.selectbox("Color scale", ["RdYlGn", "Viridis", "RdBu_r", "Plasma"], key="qc_corr_cs")
                with cc3:
                    corr_cluster = st.checkbox("Cluster samples", value=True, key="qc_corr_cluster")
                
                corr_hash = _data_hash(corr_counts)
                corr_vals = _cached_correlation(corr_hash, corr_counts.values, corr_counts.columns.tolist(), corr_method, gene_id_type=ds.gene_id_type)
                corr = pd.DataFrame(corr_vals, index=corr_counts.columns, columns=corr_counts.columns)
                
                # Cluster
                if corr_cluster and corr.shape[0] > 2:
                    try:
                        d = pdist(corr.values, 'correlation')
                        d = np.nan_to_num(d, nan=0.0)
                        Z = linkage(d, method='average')
                        o = leaves_list(Z)
                        corr = corr.iloc[o, o]
                    except Exception:
                        pass
                
                # Label with group info
                if color_by != 'None':
                    corr_labels = [f"{s} ({color_labels.get(s, '')})" for s in corr.columns]
                else:
                    corr_labels = corr.columns.tolist()
                
                od = corr.values[~np.eye(len(corr), dtype=bool)]
                zmin = max(0, np.floor(od.min() * 20) / 20)
                
                # Only show text annotations for small matrices
                show_text = len(corr) <= 50
                
                heatmap_kwargs = dict(
                    z=corr.values, x=corr_labels, y=corr_labels,
                    colorscale=corr_cs, zmin=zmin, zmax=1.0,
                    hovertemplate='%{x} vs %{y}<br>r = %{z:.4f}<extra></extra>',
                    colorbar=dict(title=dict(text=corr_method.capitalize())),
                )
                if show_text:
                    heatmap_kwargs['text'] = np.round(corr.values, 2)
                    heatmap_kwargs['texttemplate'] = '%{text}'
                    heatmap_kwargs['textfont'] = dict(size=max(6, min(9, 120 // len(corr))))
                
                fig = go.Figure(go.Heatmap(**heatmap_kwargs))
                fig.update_layout(
                    title=f"Sample Correlation ({corr_method.capitalize()})",
                    height=max(500, min(len(corr) * 35, 1200)),
                    font=dict(size=11),
                    xaxis=dict(tickangle=45, tickfont=dict(size=max(5, min(7, 900 // len(corr))))),
                    yaxis=dict(tickfont=dict(size=max(5, min(7, 900 // len(corr))))),
                )
                _show_plot(fig, "RAPTOR_qc_9", "9")
                
                st.caption(f"Range: {od.min():.3f} – {od.max():.3f} | Mean: {od.mean():.3f} | Std: {od.std():.3f}")
                
                # Smart within/between group assessment
                if color_by != 'None' and color_labels is not None and color_labels.nunique() > 1:
                    # Use vectorized approach for speed
                    labels_ord = np.array([color_labels.get(s, '') for s in corr.columns])
                    triu_idx = np.triu_indices(len(corr), k=1)
                    triu_vals = corr.values[triu_idx]
                    same_group = labels_ord[triu_idx[0]] == labels_ord[triu_idx[1]]
                    
                    within_corrs = triu_vals[same_group]
                    between_corrs = triu_vals[~same_group]
                    
                    if len(within_corrs) > 0 and len(between_corrs) > 0:
                        mw = np.mean(within_corrs)
                        mb = np.mean(between_corrs)
                        col_w, col_b = st.columns(2)
                        with col_w:
                            st.metric("Within-group", f"{mw:.3f}")
                        with col_b:
                            st.metric("Between-group", f"{mb:.3f}")
                        
                        if mw > 0.8 and mb < mw - 0.1:
                            st.success(
                                f"Strong group structure: within-group ({mw:.3f}) >> between-group ({mb:.3f}). "
                                f"Clear biological differences between groups."
                            )
                        elif mw > mb:
                            st.info(f"Moderate group structure: within ({mw:.3f}) > between ({mb:.3f}).")
                        else:
                            st.warning(f"Weak structure: within ({mw:.3f}) ≈ between ({mb:.3f}).")
            
            with qc_t6:
                st.markdown("**QC Metrics**")
                st.caption("Gene-level quality indicators for RNA-seq data")
                
                gene_ids = counts.index.astype(str)
                
                # Mitochondrial genes (MT- prefix for symbols, or ENSG known MT genes)
                mt_mask = gene_ids.str.startswith('MT-') | gene_ids.str.startswith('mt-')
                # For Ensembl IDs: known MT genes
                mt_ensembl = {
                    'ENSG00000198888', 'ENSG00000198763', 'ENSG00000198804',
                    'ENSG00000198712', 'ENSG00000228253', 'ENSG00000198899',
                    'ENSG00000198938', 'ENSG00000198840', 'ENSG00000212907',
                    'ENSG00000198886', 'ENSG00000198786', 'ENSG00000198695',
                    'ENSG00000198727',
                }
                mt_mask = mt_mask | gene_ids.isin(mt_ensembl)
                
                # Ribosomal genes (RPL/RPS prefix)
                ribo_mask = gene_ids.str.startswith('RPL') | gene_ids.str.startswith('RPS')
                ribo_mask = ribo_mask | gene_ids.str.startswith('ENSG000001544') # common ribo prefix pattern
                
                if mt_mask.any():
                    mt_counts = counts.loc[mt_mask].sum(axis=0)
                    mt_pct = (mt_counts / lib_sizes * 100)
                    
                    mt_df = pd.DataFrame({
                        'Sample': mt_pct.index,
                        'MT %': mt_pct.values,
                        'Group': color_labels.values,
                    })
                    fig_mt = px.box(mt_df, x='Group', y='MT %', color='Group', points='all',
                        title="Mitochondrial Gene Expression (%)")
                    fig_mt.update_layout(height=400, showlegend=False)
                    _show_plot(fig_mt, "RAPTOR_qc_10", "10")
                    
                    median_mt = mt_pct.median()
                    if median_mt > 20:
                        st.warning(f"High mitochondrial expression (median {median_mt:.1f}%). May indicate cell stress or dying cells.")
                    elif median_mt > 10:
                        st.info(f"Moderate MT expression (median {median_mt:.1f}%).")
                    else:
                        st.success(f"MT expression OK (median {median_mt:.1f}%).")
                else:
                    st.caption("Mitochondrial genes not detected (may need gene symbol conversion).")
                
                # Genes detected vs library size
                st.markdown("**Library Complexity**")
                complex_df = pd.DataFrame({
                    'Library Size': lib_sizes.values,
                    'Genes Detected': gene_detected.values,
                    'Group': color_labels.values,
                    'Sample': lib_sizes.index,
                })
                fig_complex = px.scatter(
                    complex_df, x='Library Size', y='Genes Detected',
                    color='Group', hover_data=['Sample'],
                    title="Library Size vs. Genes Detected",
                )
                fig_complex.update_layout(height=400)
                fig_complex.update_traces(marker=dict(size=8))
                _show_plot(fig_complex, "RAPTOR_qc_11", "11")
                st.caption(
                    "Samples should follow a saturation curve. Outliers far below the trend "
                    "may have degraded RNA or library preparation issues."
                )
                
                # Expression density — zero-count genes per sample
                st.markdown("**Zero-count genes per sample**")
                zero_per_sample = (counts == 0).sum(axis=0)
                zero_pct_per_sample = zero_per_sample / counts.shape[0] * 100
                
                zero_df = pd.DataFrame({
                    'Sample': zero_pct_per_sample.index,
                    'Zero Genes (%)': zero_pct_per_sample.values,
                    'Group': color_labels.values,
                })
                fig_zero = px.box(zero_df, x='Group', y='Zero Genes (%)', color='Group', points='all',
                    title="Percentage of Zero-Count Genes per Sample")
                fig_zero.update_layout(height=400, showlegend=False)
                _show_plot(fig_zero, "RAPTOR_qc_12", "12")
            
            # Quality Verdict
            st.markdown("---")
            st.markdown("### Quality Verdict")
            issues = []
            if lib_sizes.std() / lib_sizes.mean() > 1.0:
                issues.append("High library size variability (CV > 1.0)")
            if zero_pct > 80:
                issues.append("Very high sparsity (>80% zeros)")
            if gene_detected.median() < 10000:
                issues.append("Low gene detection rate (median < 10,000)")
            
            if not issues:
                st.success("Your dataset looks good! Proceed to DE analysis or pooling.")
            else:
                st.warning(
                    f"Found {len(issues)} potential issue(s):\n" +
                    "\n".join(f"- {issue}" for issue in issues) +
                    "\n\nConsider filtering low-count genes or checking for problematic samples."
                )
        
        # ---- POOLED DATA QC (existing) ----
        elif qc_selection == "Pooled Data" and pool is not None:
            st.success(
                f"Analyzing pooled data: {pool.n_genes:,} genes, "
                f"{pool.n_samples} samples, {pool.n_studies} studies"
            )
            
            counts = pool.counts_df
            
            # Build combined metadata from pooled data
            pool_meta = pd.DataFrame({'study': pool.study_labels})
            # Add columns from pool's own metadata (multi-project download stores clinical data here)
            if hasattr(pool, 'metadata') and pool.metadata is not None and not pool.metadata.empty:
                for col in pool.metadata.columns:
                    if col not in pool_meta.columns and col != 'study':
                        vals = pool.metadata[col].reindex(pool_meta.index)
                        if vals.notna().sum() > 0:
                            pool_meta[col] = vals
            # Also try to get metadata from original datasets in session
            for ds_name, ds_obj in st.session_state.get('acq_datasets', {}).items():
                if hasattr(ds_obj, 'metadata') and not ds_obj.metadata.empty:
                    for col in ds_obj.metadata.columns:
                        if col not in pool_meta.columns:
                            matching = ds_obj.metadata[col].reindex(pool_meta.index)
                            if matching.notna().sum() > 0:
                                pool_meta[col] = matching
            
            # Color-by selector
            pool_color_options = ['study']
            for col in pool_meta.columns:
                if col != 'study' and pool_meta[col].nunique() >= 2 and pool_meta[col].nunique() <= 20:
                    pool_color_options.append(col)
            
            pool_color_by = st.selectbox(
                "Color samples by",
                pool_color_options,
                index=0,
                key='pool_qc_color_by',
                help="Choose metadata column to group samples. 'study' shows project origin."
            )
            
            pool_color_labels = pool_meta[pool_color_by].reindex(counts.columns).fillna('unknown')
            
            # Cached metrics
            pool_dh = _data_hash(counts)
            lib_vals_p, det_vals_p, zero_pct = _cached_lib_stats(pool_dh, counts.values)
            lib_sizes = pd.Series(lib_vals_p, index=counts.columns)
            gene_detected = pd.Series(det_vals_p, index=counts.columns)
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Median Library Size", f"{lib_sizes.median():,.0f}")
            with col2:
                st.metric("Library Size CV", f"{lib_sizes.std() / lib_sizes.mean():.2f}")
            with col3:
                st.metric("Median Genes Detected", f"{gene_detected.median():,.0f}")
            with col4:
                st.metric("Zero %", f"{zero_pct:.1f}%")
            
            qc_tab1, qc_tab2, qc_tab3, qc_tab4, qc_tab5, qc_tab6 = st.tabs([
                "Library Sizes", "Batch Assessment", "PCA", "Gene Expression",
                "Heatmap", "Sample Correlation"
            ])
            
            with qc_tab1:
                lib_df = pd.DataFrame({
                    'Sample': lib_sizes.index,
                    'Library Size': lib_sizes.values,
                    'Group': pool_color_labels.values,
                })
                fig = px.box(lib_df, x='Group', y='Library Size', color='Group',
                             title=f"Library Size Distribution by {pool_color_by}", points='all')
                fig.update_layout(height=450, showlegend=False)
                _show_plot(fig, "RAPTOR_qc_13", "13")
                
                group_medians = lib_df.groupby('Group')['Library Size'].median()
                fold_range = group_medians.max() / group_medians.min() if group_medians.min() > 0 else float('inf')
                if fold_range > 5:
                    st.warning(f"Library size varies {fold_range:.1f}-fold across groups.")
                elif fold_range > 2:
                    st.info(f"Library size varies {fold_range:.1f}-fold across groups.")
                else:
                    st.success(f"Library sizes are consistent ({fold_range:.1f}-fold range).")
            
            with qc_tab2:
                st.markdown("**Study-level batch assessment**")
                study_means = {}
                for study in pool.studies:
                    study_samples = pool.get_study_samples(study)
                    study_cols = [s for s in study_samples if s in counts.columns]
                    if study_cols:
                        study_means[study] = counts[study_cols].mean(axis=1)
                
                if len(study_means) >= 2:
                    means_df = pd.DataFrame(study_means)
                    corr_matrix = means_df.corr()
                    fig = go.Figure(data=go.Heatmap(
                        z=corr_matrix.values, x=corr_matrix.columns, y=corr_matrix.index,
                        colorscale='RdYlGn', text=np.round(corr_matrix.values, 3),
                        texttemplate='%{text}', colorbar=dict(title="Correlation"), zmin=0, zmax=1,
                    ))
                    fig.update_layout(title="Mean Expression Correlation Between Studies", height=400)
                    _show_plot(fig, "RAPTOR_qc_14", "14")
                    
                    min_corr = corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].min()
                    if min_corr < 0.8:
                        st.warning(f"Lowest inter-study correlation: {min_corr:.3f}. Batch correction recommended.")
                    elif min_corr < 0.95:
                        st.info(f"Lowest inter-study correlation: {min_corr:.3f}. Moderate batch effects.")
                    else:
                        st.success(f"Lowest inter-study correlation: {min_corr:.3f}. Well-correlated.")
                
                sps = pool.samples_per_study
                min_s, max_s = min(sps.values()), max(sps.values())
                imbalance = max_s / min_s if min_s > 0 else float('inf')
                if imbalance > 5:
                    st.warning(f"Sample imbalance: {imbalance:.1f}x ({min_s} to {max_s}).")
                else:
                    st.success(f"Sample balance OK ({min_s} to {max_s} per study).")
            
            with qc_tab3:
                colors_pal = px.colors.qualitative.Set2
                
                # Cached PCA on top 5K variable genes
                pool_dhash = _data_hash(counts)
                pcs_p, var_ratios_p, n_comp_p = _cached_pca(
                    pool_dhash, counts.values, counts.shape[1], n_comp=6
                )
                
                pc_cols_p = [f'PC{i+1}' for i in range(n_comp_p)]
                var_pct_p = {f'PC{i+1}': var_ratios_p[i]*100 for i in range(n_comp_p)}
                var_lab_p = {pc: f'{pc} ({var_pct_p[pc]:.1f}%)' for pc in pc_cols_p}
                
                pca_df_p = pd.DataFrame(pcs_p, columns=pc_cols_p, index=counts.columns)
                pca_df_p['Sample'] = pca_df_p.index
                pca_df_p['study'] = pool.study_labels.reindex(pca_df_p.index).values
                
                # Add metadata columns for coloring
                pca_pool_color_opts = ['study']
                for col in pool_meta.columns:
                    if col != 'study' and pool_meta[col].nunique() >= 2 and pool_meta[col].nunique() <= 20:
                        pca_df_p[col] = pool_meta[col].reindex(pca_df_p.index)
                        if col not in pca_pool_color_opts:
                            pca_pool_color_opts.append(col)
                
                ctrl1, ctrl2, ctrl3, ctrl4 = st.columns(4)
                with ctrl1:
                    pca_pool_color = st.selectbox("Color by", pca_pool_color_opts,
                        index=0, key="pool_pca_color")
                with ctrl2:
                    pca_px = st.selectbox("X axis", pc_cols_p, index=0, key="pool_pca_x")
                with ctrl3:
                    pca_py = st.selectbox("Y axis", pc_cols_p, index=min(1, n_comp_p-1), key="pool_pca_y")
                with ctrl4:
                    pca_pdim = st.selectbox("Dimension", ["2D","3D"] if n_comp_p >= 3 else ["2D"], key="pool_pca_dim")
                
                if pca_pdim == "3D":
                    pca_pz = st.selectbox("Z axis", pc_cols_p, index=min(2, n_comp_p-1), key="pool_pca_z")
                
                if pca_pdim == "2D":
                    fig = px.scatter(pca_df_p, x=pca_px, y=pca_py, color=pca_pool_color,
                        color_discrete_sequence=colors_pal, hover_data=['Sample'],
                        title=f"PCA — {var_lab_p[pca_px]} vs {var_lab_p[pca_py]}")
                    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='white')))
                    _pub_layout(fig, f"PCA — {var_lab_p[pca_px]} vs {var_lab_p[pca_py]}", width=1000, height=700)
                else:
                    fig = px.scatter_3d(pca_df_p, x=pca_px, y=pca_py, z=pca_pz, color=pca_pool_color,
                        color_discrete_sequence=colors_pal, hover_data=['Sample'],
                        title=f"3D PCA — {var_lab_p[pca_px]} vs {var_lab_p[pca_py]} vs {var_lab_p[pca_pz]}")
                    fig.update_traces(marker=dict(size=6, line=dict(width=0.5, color='white')))
                    fig.update_layout(
                        title=dict(text=f"<b>3D PCA — {var_lab_p[pca_px]} vs {var_lab_p[pca_py]} vs {var_lab_p[pca_pz]}</b>",
                                   font=dict(size=16, family="Arial")),
                        height=750, font=PUB_FONT, paper_bgcolor='white',
                        scene=dict(xaxis_title=var_lab_p[pca_px], yaxis_title=var_lab_p[pca_py], zaxis_title=var_lab_p[pca_pz]),
                        legend=dict(font=dict(size=13), bordercolor="gray", borderwidth=1),
                    )
                
                _show_plot(fig, "RAPTOR_pooled_PCA", "15")
                
                dl1p, dl2p = st.columns(2)
                with dl1p:
                    try:
                        png_bytes = fig.to_image(format="png", width=1200, height=900, scale=3, engine="kaleido")
                        st.download_button(
                            "Download PNG (300 DPI)", data=png_bytes,
                            file_name="RAPTOR_pooled_PCA.png", mime="image/png",
                            key="dl_pool_pca_png",
                        )
                    except Exception:
                        st.caption("Install kaleido for PNG export: `pip install kaleido`")
                with dl2p:
                    try:
                        svg_bytes = fig.to_image(format="svg", width=1200, height=900, scale=3, engine="kaleido")
                        st.download_button(
                            "Download SVG (vector)", data=svg_bytes,
                            file_name="RAPTOR_pooled_PCA.svg", mime="image/svg+xml",
                            key="dl_pool_pca_svg",
                        )
                    except Exception:
                        st.caption("Install kaleido for SVG export: `pip install kaleido`")
                
                with st.expander("Scree plot"):
                    scree = pd.DataFrame({'Component': pc_cols_p, 'Variance (%)': [var_pct_p[pc] for pc in pc_cols_p],
                        'Cumulative (%)': list(np.cumsum([var_pct_p[pc] for pc in pc_cols_p]))})
                    fig_s = go.Figure()
                    fig_s.add_bar(x=scree['Component'], y=scree['Variance (%)'], name='Individual', marker_color='#66c2a5')
                    fig_s.add_scatter(x=scree['Component'], y=scree['Cumulative (%)'], name='Cumulative', mode='lines+markers', marker_color='#fc8d62')
                    fig_s.update_layout(yaxis_title='Variance (%)', height=300, plot_bgcolor='#fafafa')
                    _show_plot(fig_s, "RAPTOR_qc_16", "16")
            
            with qc_tab4:
                st.markdown("**Expression distributions**")
                sample_genes = counts.index[:2000] if len(counts) > 2000 else counts.index
                fig = go.Figure()
                for grp in pool_color_labels.unique():
                    grp_cols = [s for s, g in zip(counts.columns, pool_color_labels) if g == grp]
                    if grp_cols:
                        mean_expr = np.log2(counts.loc[sample_genes, grp_cols].mean(axis=1) + 1)
                        fig.add_trace(go.Violin(y=mean_expr, name=str(grp), box_visible=True, meanline_visible=True))
                fig.update_layout(title=f"Log2 Mean Expression by {pool_color_by}", yaxis_title="Log2(mean count + 1)", height=450)
                _show_plot(fig, "RAPTOR_qc_17", "17")
            
            with qc_tab5:
                from scipy.cluster.hierarchy import linkage, leaves_list
                from scipy.spatial.distance import pdist
                
                st.markdown("**Heatmap — Top Variable Genes**")
                
                MAX_SAMPLES_HEATMAP = 200
                hm_counts = counts
                if counts.shape[1] > MAX_SAMPLES_HEATMAP:
                    st.warning(
                        f"Dataset has {counts.shape[1]} samples. Subsampling to {MAX_SAMPLES_HEATMAP} "
                        f"for heatmap rendering (stratified by group). Biological patterns are preserved."
                    )
                    # Stratified subsample
                    subsample_idx = []
                    for grp in pool_color_labels.unique():
                        grp_samples = [s for s, g in zip(counts.columns, pool_color_labels) if g == grp]
                        n_take = max(1, int(MAX_SAMPLES_HEATMAP * len(grp_samples) / counts.shape[1]))
                        subsample_idx.extend(grp_samples[:n_take])
                    hm_counts = counts[subsample_idx[:MAX_SAMPLES_HEATMAP]]
                
                hm1, hm2, hm3 = st.columns(3)
                with hm1:
                    n_genes_hm = st.slider("Top variable genes", 20, 200, 50, 10, key="pool_hm_ngenes")
                with hm2:
                    hm_cs = st.selectbox("Color scale", ["RdBu_r", "Viridis", "Plasma", "PiYG"], key="pool_hm_cs")
                with hm3:
                    hm_cluster = st.checkbox("Cluster rows & columns", value=True, key="pool_hm_cluster")
                
                log_hm = np.log2(hm_counts + 1)
                gene_var_hm = log_hm.var(axis=1).nlargest(n_genes_hm)
                subset = log_hm.loc[gene_var_hm.index]
                row_m = subset.mean(axis=1)
                row_s = subset.std(axis=1).replace(0, 1)
                subset_z = subset.subtract(row_m, axis=0).divide(row_s, axis=0)
                
                if hm_cluster and subset_z.shape[0] > 2 and subset_z.shape[1] > 2:
                    try:
                        r_d = pdist(subset_z.values, 'correlation')
                        c_d = pdist(subset_z.T.values, 'correlation')
                        r_d = np.nan_to_num(r_d, nan=0.0)
                        c_d = np.nan_to_num(c_d, nan=0.0)
                        subset_z = subset_z.iloc[leaves_list(linkage(r_d, 'average')),
                                                  leaves_list(linkage(c_d, 'average'))]
                    except Exception:
                        pass
                
                fig = go.Figure(go.Heatmap(
                    z=subset_z.values, x=subset_z.columns.tolist(), y=subset_z.index.tolist(),
                    colorscale=hm_cs, zmid=0,
                    hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>',
                    colorbar=dict(title=dict(text="Z-score"), len=0.7),
                ))
                fig.update_layout(
                    title=f"Top {n_genes_hm} Variable Genes (Z-score)",
                    height=max(500, n_genes_hm * 12),
                    xaxis=dict(tickangle=45, tickfont=dict(size=6)),
                    yaxis=dict(tickfont=dict(size=max(6, min(10, 500 // n_genes_hm)))),
                )
                _show_plot(fig, "RAPTOR_qc_18", "18")
            
            with qc_tab6:
                st.markdown("**Sample-to-sample correlation**")
                
                MAX_SAMPLES_CORR = 150
                corr_counts = counts
                corr_labels_sub = pool_color_labels
                
                if counts.shape[1] > MAX_SAMPLES_CORR:
                    st.warning(
                        f"Dataset has {counts.shape[1]} samples. Subsampling to {MAX_SAMPLES_CORR} "
                        f"for correlation heatmap (a {counts.shape[1]}×{counts.shape[1]} matrix would freeze the browser)."
                    )
                    sub_idx = []
                    for grp in pool_color_labels.unique():
                        grp_samples = [s for s, g in zip(counts.columns, pool_color_labels) if g == grp]
                        n_take = max(1, int(MAX_SAMPLES_CORR * len(grp_samples) / counts.shape[1]))
                        sub_idx.extend(grp_samples[:n_take])
                    sub_idx = sub_idx[:MAX_SAMPLES_CORR]
                    corr_counts = counts[sub_idx]
                    corr_labels_sub = pool_color_labels.loc[sub_idx]
                
                cc1, cc2 = st.columns(2)
                with cc1:
                    corr_m = st.selectbox("Method", ["pearson", "spearman"], key="pool_corr_m")
                with cc2:
                    corr_cl = st.checkbox("Cluster samples", value=True, key="pool_corr_cl")
                
                corr_vals = _cached_correlation(
                    _data_hash(corr_counts), corr_counts.values,
                    corr_counts.columns.tolist(), corr_m
                )
                corr_df = pd.DataFrame(corr_vals, index=corr_counts.columns, columns=corr_counts.columns)
                
                if corr_cl and corr_df.shape[0] > 2:
                    try:
                        d = pdist(corr_df.values, 'correlation')
                        d = np.nan_to_num(d, nan=0.0)
                        o = leaves_list(linkage(d, 'average'))
                        corr_df = corr_df.iloc[o, o]
                    except Exception:
                        pass
                
                od = corr_df.values[~np.eye(len(corr_df), dtype=bool)]
                zmin = max(0, np.floor(od.min() * 20) / 20)
                
                fig = go.Figure(go.Heatmap(
                    z=corr_df.values, x=corr_df.columns.tolist(), y=corr_df.index.tolist(),
                    colorscale='RdYlGn', zmin=zmin, zmax=1.0,
                    hovertemplate='%{x} vs %{y}<br>r = %{z:.3f}<extra></extra>',
                    colorbar=dict(title=dict(text=corr_m.capitalize())),
                ))
                fig.update_layout(
                    title=f"Sample Correlation ({corr_m.capitalize()})",
                    height=max(500, len(corr_df) * 4),
                    xaxis=dict(tickangle=45, tickfont=dict(size=6)),
                    yaxis=dict(tickfont=dict(size=6)),
                )
                _show_plot(fig, "RAPTOR_qc_19", "19")
                st.caption(f"Range: {od.min():.3f} – {od.max():.3f} | Mean: {od.mean():.3f}")
                
                # Within vs between group
                if pool_color_labels.nunique() > 1:
                    within, between = [], []
                    labs = [corr_labels_sub.get(s, '') for s in corr_df.columns]
                    for i in range(len(labs)):
                        for j in range(i+1, len(labs)):
                            if labs[i] == labs[j]:
                                within.append(corr_df.values[i, j])
                            else:
                                between.append(corr_df.values[i, j])
                    if within and between:
                        cw, cb = st.columns(2)
                        with cw:
                            st.metric("Within-group", f"{np.mean(within):.3f}")
                        with cb:
                            st.metric("Between-group", f"{np.mean(between):.3f}")

# =============================================================================
# TAB 5: Export
# =============================================================================

with tab_export:
    st.markdown("## Export Data")
    st.caption("Download your datasets as spreadsheet-compatible CSV files, or pass them to the next analysis step")
    
    session_ds = st.session_state.get('acq_datasets', {})
    pool = st.session_state.get('acq_pooled')
    
    # Export individual datasets
    st.markdown("### Individual Datasets")
    
    if session_ds:
        export_name = st.selectbox(
            "Select dataset to export",
            list(session_ds.keys()),
            key='export_ds_select'
        )
        
        ds = session_ds[export_name]
        
        # Detect data type for appropriate display
        is_mutation = hasattr(ds, 'data_type') and 'mutation' in str(ds.data_type).lower()
        dtype = ds.source_info.get('data_type', 'raw_counts')
        
        _EXPORT_TYPE_LABELS = {
            'somatic_mutation': 'MAF table',
            'mirna_rpm': 'miRNA expression matrix',
            'cnv_gene_level': 'CNV gene-level matrix',
            'cnv_segment': 'CNV segment table',
            'methylation_beta': 'Methylation beta-value matrix',
            'protein_expression_rppa': 'RPPA protein expression matrix',
        }
        type_label = _EXPORT_TYPE_LABELS.get(dtype, 'Count matrix')
        
        # Show dataset size info with appropriate label
        mem_mb = ds.counts_df.memory_usage(deep=True).sum() / 1024 / 1024
        is_large = mem_mb > 50  # Large dataset warning threshold
        
        if is_mutation:
            n_mutations = len(ds.counts_df)
            n_cases = ds.counts_df['Tumor_Sample_Barcode'].nunique() if 'Tumor_Sample_Barcode' in ds.counts_df.columns else ds.n_samples
            st.caption(f"MAF table: {n_mutations:,} mutations from {n_cases} cases | {len(ds.counts_df.columns)} columns | Memory: {mem_mb:.1f} MB")
        else:
            st.caption(f"{type_label}: {ds.n_genes:,} × {ds.n_samples} | Memory: {mem_mb:.1f} MB")
        
        if is_large:
            st.caption("Large dataset — exports are generated on-demand to avoid freezing.")
        
        if is_mutation:
            # ---- Mutation-specific export ----
            
            # Key columns (small, fast)
            key_maf_cols = [c for c in [
                'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position',
                'Variant_Classification', 'Variant_Type', 'Reference_Allele',
                'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'HGVSp_Short',
            ] if c in ds.counts_df.columns]
            
            col1, col2 = st.columns(2)
            
            with col1:
                if st.button(f"Prepare MAF — Key Columns ({len(key_maf_cols)} cols)", 
                             key='exp_maf_key', use_container_width=True):
                    with st.spinner("Generating CSV..."):
                        st.session_state['_export_maf_key'] = ds.counts_df[key_maf_cols].to_csv(index=False)
                
                if '_export_maf_key' in st.session_state:
                    st.download_button(
                        "Download Key Columns CSV",
                        data=st.session_state['_export_maf_key'],
                        file_name=f"{export_name}_key_columns.csv",
                        mime="text/csv", use_container_width=True,
                    )
            
            with col2:
                if st.button(f"Prepare Full MAF ({len(ds.counts_df.columns)} cols)",
                             key='exp_maf_full', use_container_width=True):
                    with st.spinner("Generating CSV (this may take a minute)..."):
                        st.session_state['_export_maf_full'] = ds.counts_df.to_csv(index=False)
                
                if '_export_maf_full' in st.session_state:
                    st.download_button(
                        "Download Full MAF CSV",
                        data=st.session_state['_export_maf_full'],
                        file_name=f"{export_name}_full.csv",
                        mime="text/csv", use_container_width=True,
                    )
            
            # Mutation summary (always fast — just gene counts)
            if 'Hugo_Symbol' in ds.counts_df.columns:
                with st.expander("Mutation summary by gene"):
                    gene_counts = ds.counts_df['Hugo_Symbol'].value_counts().head(30)
                    st.dataframe(
                        pd.DataFrame({'Gene': gene_counts.index, 'Mutations': gene_counts.values}),
                        use_container_width=True, hide_index=True,
                    )
                    summary_csv = gene_counts.reset_index()
                    summary_csv.columns = ['Gene', 'Mutations']
                    st.download_button(
                        "Download Mutation Summary (CSV)",
                        data=summary_csv.to_csv(index=False),
                        file_name=f"{export_name}_gene_summary.csv",
                        mime="text/csv",
                    )
        
        else:
            # ---- Standard expression/matrix export ----
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if is_large:
                    if st.button(f"Prepare {type_label} (CSV)", key='exp_csv', use_container_width=True):
                        with st.spinner(f"Generating CSV ({mem_mb:.0f} MB)..."):
                            st.session_state['_export_csv'] = ds.counts_df.to_csv()
                    if '_export_csv' in st.session_state:
                        st.download_button(
                            f"Download {type_label} (CSV)",
                            data=st.session_state['_export_csv'],
                            file_name=f"{export_name}_counts.csv",
                            mime="text/csv", use_container_width=True,
                        )
                else:
                    csv_data = ds.counts_df.to_csv()
                    st.download_button(
                        f"Download {type_label} (CSV)",
                        data=csv_data,
                        file_name=f"{export_name}_counts.csv",
                        mime="text/csv", use_container_width=True,
                    )
            
            with col2:
                try:
                    if is_large:
                        if st.button(f"Prepare {type_label} (Parquet)", key='exp_pq', use_container_width=True):
                            with st.spinner("Generating Parquet..."):
                                buf = BytesIO()
                                ds.counts_df.to_parquet(buf, engine='pyarrow', compression='snappy')
                                st.session_state['_export_pq'] = buf.getvalue()
                        if '_export_pq' in st.session_state:
                            st.download_button(
                                f"Download {type_label} (Parquet) ⚡",
                                data=st.session_state['_export_pq'],
                                file_name=f"{export_name}_counts.parquet",
                                mime="application/octet-stream", use_container_width=True,
                            )
                    else:
                        buf = BytesIO()
                        ds.counts_df.to_parquet(buf, engine='pyarrow', compression='snappy')
                        st.download_button(
                            f"Download {type_label} (Parquet) ⚡",
                            data=buf.getvalue(),
                            file_name=f"{export_name}_counts.parquet",
                            mime="application/octet-stream", use_container_width=True,
                            help="Parquet is 4-10x faster to load than CSV and 50% smaller."
                        )
                except Exception:
                    st.caption("Parquet export requires pyarrow: pip install pyarrow")
            
            with col3:
                meta_csv = ds.metadata.to_csv()
                st.download_button(
                    "Download Metadata (CSV)",
                    data=meta_csv,
                    file_name=f"{export_name}_metadata.csv",
                    mime="text/csv", use_container_width=True,
                )
    else:
        st.info("No datasets loaded. Use the Search tab first.")
    
    # Export pooled dataset
    st.markdown("---")
    st.markdown("### Pooled Dataset")
    
    if pool is not None:
        st.success(
            f"Pooled: {pool.n_genes:,} genes, {pool.n_samples} samples, {pool.n_studies} studies"
        )
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            pool_csv = pool.counts_df.to_csv()
            st.download_button(
                "Download Pooled Counts (CSV)",
                data=pool_csv,
                file_name="pooled_counts.csv",
                mime="text/csv",
                use_container_width=True,
            )
        
        with col2:
            pool_meta = pool.metadata.to_csv()
            st.download_button(
                "Download Pooled Metadata (CSV)",
                data=pool_meta,
                file_name="pooled_metadata.csv",
                mime="text/csv",
                use_container_width=True,
            )
        
        with col3:
            pool_info = json.dumps({
                'pooling_info': pool.pooling_info,
                'component_datasets': pool.component_datasets,
                'study_labels': pool.study_labels.to_dict(),
                'n_studies': pool.n_studies,
                'n_genes': pool.n_genes,
                'n_samples': pool.n_samples,
            }, indent=2, default=str)
            st.download_button(
                "Download Pool Info (JSON)",
                data=pool_info,
                file_name="pool_info.json",
                mime="application/json",
                use_container_width=True,
            )
        
        # Store in session state for Module 10
        st.markdown("---")
        st.markdown("### Pass to Downstream Pages")
        st.info(
            "The pooled dataset is stored in session state. "
            "Quality Assessment, Import DE, and other downstream pages "
            "will automatically detect it."
        )
        
        if st.button("Store for Downstream Analysis", type="primary"):
            st.session_state['m6b_pooled_dataset'] = pool
            st.session_state['m6b_pooled_counts'] = pool.counts_df
            st.session_state['m6b_pooled_metadata'] = pool.metadata
            st.session_state['m6b_study_labels'] = pool.study_labels
            st.success("Stored in session state. Navigate to Quality Assessment or Import DE.")
    else:
        st.info("No pooled dataset available. Use the Pool tab first.")

# =============================================================================
# Footer
# =============================================================================

st.markdown("---")
st.caption("**Tip:** Upload your own data + pool with 2–3 public GEO studies of the same condition for robust cross-cohort biomarker validation")
st.caption("**Next step:** Go to Quality Assessment to check your data, or Import DE to load differential expression results")
