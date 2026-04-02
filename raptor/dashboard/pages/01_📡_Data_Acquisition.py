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
import json
from datetime import datetime

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
    
    # ---- PREVIEW INFO ----
    if preview_clicked and direct_accession:
        accession = direct_accession.strip()
        repo = _detect_repo(accession)
        
        with st.spinner(f"Fetching info for {accession} from {repo}..."):
            try:
                if repo == 'GEO':
                    connector = GEOConnector(cache=True, cache_dir=cache.cache_dir)
                    info = connector.get_study_info(accession)
                else:
                    connector = _get_connector(repo)
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
        
        # Key metrics row
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Samples", info.get('n_samples', 'N/A'))
        with col2:
            st.metric("Organism", info.get('organism', 'unknown'))
        with col3:
            platform = info.get('platform_id', 'unknown')
            st.metric("Platform", platform)
        with col4:
            data_type = info.get('data_type', 'unknown')
            st.metric("Data Type", data_type)
        
        # Platform name
        platform_name = info.get('platform_name', '')
        if platform_name:
            st.caption(f"Platform: {platform_name}")
        
        # Summary
        summary = info.get('summary', '')
        if summary:
            with st.expander("Study Summary", expanded=True):
                st.write(summary)
        
        # Overall design
        design = info.get('overall_design', '')
        if design:
            with st.expander("Experimental Design"):
                st.write(design)
        
        # Sample table
        samples_df = info.get('samples')
        if samples_df is not None and isinstance(samples_df, pd.DataFrame) and not samples_df.empty:
            with st.expander(f"Sample Details ({len(samples_df)} samples)", expanded=True):
                st.dataframe(samples_df, use_container_width=True, hide_index=True)
        
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
                    f"{ds.n_genes:,} genes, {ds.n_samples} samples "
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
            type_options = ["RNA-seq"]  # TCGA genomic data
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
                table_rows.append({
                    'Accession': r.accession,
                    'Title': r.title[:70] + ('...' if len(r.title) > 70 else ''),
                    'Type': r.extra.get('gds_type', ''),
                    'Organism': r.organism,
                    'Samples': r.n_samples,
                    'Platform': r.platform,
                    'GPL': f"GPL{gpl_id}" if gpl_id and not gpl_id.startswith('GPL') else gpl_id,
                    'PubMed': pubmed[0] if pubmed else '',
                    'Cached': '✅' if cache.is_cached(search_repo, r.accession) else '',
                })
        
        results_df = pd.DataFrame(table_rows)
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
                                st.success(f"Downloaded {selected_accession}: {ds.n_genes:,} genes, {ds.n_samples} samples")
                        except Exception as e:
                            st.error(f"Failed: {e}")
        
        with col_info:
            if search_repo == 'SRA':
                sra_gse_link = selected_result.extra.get('gse', '') if selected_result else ''
                if sra_gse_link:
                    if st.button("Open in GEO", use_container_width=True, key='search_link'):
                        st.markdown(
                            f"[Open {sra_gse_link} in GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sra_gse_link}) "
                            f"| [Open {selected_accession} in SRA](https://www.ncbi.nlm.nih.gov/sra/?term={selected_accession})"
                        )
                else:
                    if st.button("Open in SRA", use_container_width=True, key='search_link'):
                        st.markdown(f"[Open {selected_accession} in SRA](https://www.ncbi.nlm.nih.gov/sra/?term={selected_accession})")
            else:
                if st.button("Open in GEO", use_container_width=True, key='search_link'):
                    st.markdown(f"[Open {selected_accession} in GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={selected_accession})")
        
        # Display preview info from search results
        if (st.session_state.get('preview_accession') == selected_accession
                and 'preview_info' in st.session_state):
            info = st.session_state['preview_info']
            
            st.markdown("---")
            st.markdown(f"### {info.get('title', 'No title')}")
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Samples", info.get('n_samples', 'N/A'))
            with col2:
                st.metric("Organism", info.get('organism', 'unknown'))
            with col3:
                st.metric("Platform", info.get('platform_id', 'unknown'))
            with col4:
                st.metric("Data Type", info.get('data_type', 'unknown'))
            
            pname = info.get('platform_name', '')
            if pname:
                st.caption(f"Platform: {pname}")
            
            summary = info.get('summary', '')
            if summary:
                with st.expander("Study Summary", expanded=True):
                    st.write(summary)
            
            samples_df = info.get('samples')
            if samples_df is not None and isinstance(samples_df, pd.DataFrame) and not samples_df.empty:
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
        with st.expander("Browse TCGA Projects"):
            tcga = TCGAConnector(cache=False)
            projects = tcga.list_projects()
            proj_df = pd.DataFrame([
                {'Project': k, 'Description': v}
                for k, v in projects.items()
            ])
            st.dataframe(proj_df, use_container_width=True, hide_index=True)
    
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
                st.success(f"Imported: {ds.n_genes:,} genes, {ds.n_samples} samples")
                
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
                expr_rows.append({
                    'Name': name,
                    'Repository': ds.repository,
                    'Genes': f"{ds.n_genes:,}",
                    'Samples': ds.n_samples,
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
                # ----- Expression dataset detail view -----
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Genes", f"{ds.n_genes:,}")
                with col2:
                    st.metric("Samples", ds.n_samples)
                with col3:
                    st.metric("Organism", ds.organism)
                with col4:
                    integrity = ds.validate_integrity()
                    status = "OK" if integrity['valid'] else f"{len(integrity['errors'])} errors"
                    st.metric("Integrity", status)
                
                if integrity['warnings']:
                    for w in integrity['warnings']:
                        st.warning(w)
                
                # Preview
                with st.expander("Preview count matrix"):
                    st.dataframe(ds.counts_df.head(20), use_container_width=True)
                
                with st.expander("Preview metadata"):
                    st.dataframe(ds.metadata.head(20), use_container_width=True)
            
            if not is_sra_data:
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
    
    if session_ds:
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
    
    # Cached datasets
    st.markdown("---")
    st.markdown("### Cached on Disk")
    
    cached_list = cache.list_cached_datasets()
    
    # Initialize project tags in session state
    if 'project_tags' not in st.session_state:
        st.session_state['project_tags'] = {}
    
    if cached_list:
        # Add project tags to display
        cached_df = pd.DataFrame(cached_list)
        cached_df['project'] = cached_df['accession'].map(
            lambda a: st.session_state['project_tags'].get(a, '')
        )
        
        # Fix misleading columns for SRA run tables
        # SRA entries have data_type='run_metadata' and n_genes is actually run count
        display_rows = []
        for _, row in cached_df.iterrows():
            r = dict(row)
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
                    'samples': r.get('n_samples', ''),
                    'data_type': r.get('data_type', ''),
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
        
        with col_load_sel:
            if st.button("Load Selected", use_container_width=True, key='load_selected_btn'):
                if selected_to_load:
                    loaded = 0
                    for acc in selected_to_load:
                        repo = next(d['repository'] for d in cached_list if d['accession'] == acc)
                        data = cache.load_dataset(repo, acc)
                        if data:
                            ds = AcquiredDataset(
                                counts_df=data['counts_df'],
                                metadata=data['metadata'],
                                source_info=data['source_info'],
                                gene_id_type=data['gene_id_type'],
                            )
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
                for d in cached_list:
                    acc = d['accession']
                    data = cache.load_dataset(d['repository'], acc)
                    if data:
                        ds = AcquiredDataset(
                            counts_df=data['counts_df'],
                            metadata=data['metadata'],
                            source_info=data['source_info'],
                            gene_id_type=data['gene_id_type'],
                        )
                        st.session_state['acq_datasets'][acc] = ds
                        loaded += 1
                st.success(f"Loaded all {loaded} dataset(s) into session")
                st.rerun()
        
        with col_del:
            if selected_to_load and st.button("Delete Selected", use_container_width=True, key='delete_selected_btn'):
                deleted = 0
                for acc in selected_to_load:
                    repo = next(d['repository'] for d in cached_list if d['accession'] == acc)
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
    
    # Filter out SRA run tables — they contain metadata, not expression counts
    poolable_ds = {}
    sra_skipped = []
    for name, ds in session_ds.items():
        is_sra = (
            ds.data_type == 'run_metadata'
            or (hasattr(ds, 'source_info') and ds.source_info.get('repository') == 'SRA')
        )
        if is_sra:
            sra_skipped.append(name)
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
            
            target_gene_id = st.selectbox(
                "Target gene ID type",
                ["symbol", "ensembl", "entrez"],
                help="Harmonize all datasets to this ID type (requires mygene)"
            )
            
            # Run pooling
            st.markdown("### 3. Run Pooling")
            
            if st.button("Pool Selected Datasets", type="primary", use_container_width=True):
                datasets_to_pool = [poolable_ds[n] for n in selected_for_pool]
                
                with st.spinner("Pooling datasets... (harmonizing gene IDs, merging, correcting batch effects)"):
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
    st.markdown("## Quality Check on Pooled Data")
    st.caption("Make sure your combined dataset is ready for analysis — RAPTOR checks for common issues automatically")
    
    pool = st.session_state.get('acq_pooled')
    
    if pool is None:
        st.info("Pool datasets first (Pool tab), then come here to check quality.")
    else:
        st.success(
            f"Analyzing pooled data: {pool.n_genes:,} genes, "
            f"{pool.n_samples} samples, {pool.n_studies} studies"
        )
        
        # Basic QC metrics
        st.markdown("### Basic QC Metrics")
        
        counts = pool.counts_df
        
        # Library sizes
        lib_sizes = counts.sum(axis=0)
        gene_detected = (counts > 0).sum(axis=0)
        zero_pct = (counts == 0).sum().sum() / (counts.shape[0] * counts.shape[1]) * 100
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Median Library Size", f"{lib_sizes.median():,.0f}")
        with col2:
            st.metric("Library Size CV", f"{lib_sizes.std() / lib_sizes.mean():.2f}")
        with col3:
            st.metric("Median Genes Detected", f"{gene_detected.median():,.0f}")
        with col4:
            st.metric("Zero %", f"{zero_pct:.1f}%")
        
        # Tabs for different QC plots
        qc_tab1, qc_tab2, qc_tab3, qc_tab4 = st.tabs([
            "Library Sizes", "Batch Assessment", "PCA", "Gene Expression"
        ])
        
        with qc_tab1:
            # Library size distribution by study
            lib_df = pd.DataFrame({
                'Sample': lib_sizes.index,
                'Library Size': lib_sizes.values,
                'Study': pool.study_labels.loc[lib_sizes.index].values,
            })
            
            fig = px.box(
                lib_df, x='Study', y='Library Size',
                color='Study',
                title="Library Size Distribution by Study",
                points='all',
            )
            fig.update_layout(height=450, showlegend=False)
            st.plotly_chart(fig, use_container_width=True)
            
            # Detect imbalance
            study_medians = lib_df.groupby('Study')['Library Size'].median()
            fold_range = study_medians.max() / study_medians.min() if study_medians.min() > 0 else float('inf')
            
            if fold_range > 5:
                st.warning(
                    f"Library size varies {fold_range:.1f}-fold across studies. "
                    f"Consider normalization or batch correction."
                )
            elif fold_range > 2:
                st.info(f"Library size varies {fold_range:.1f}-fold across studies. Moderate variation.")
            else:
                st.success(f"Library sizes are consistent ({fold_range:.1f}-fold range).")
        
        with qc_tab2:
            # Batch effect assessment
            st.markdown("**Study-level batch assessment**")
            
            # Per-gene mean expression by study
            study_means = {}
            for study in pool.studies:
                study_samples = pool.get_study_samples(study)
                study_cols = [s for s in study_samples if s in counts.columns]
                if study_cols:
                    study_means[study] = counts[study_cols].mean(axis=1)
            
            if len(study_means) >= 2:
                # Correlation between study means
                means_df = pd.DataFrame(study_means)
                corr_matrix = means_df.corr()
                
                fig = go.Figure(data=go.Heatmap(
                    z=corr_matrix.values,
                    x=corr_matrix.columns,
                    y=corr_matrix.index,
                    colorscale='RdYlGn',
                    text=np.round(corr_matrix.values, 3),
                    texttemplate='%{text}',
                    colorbar=dict(title="Correlation"),
                    zmin=0, zmax=1,
                ))
                fig.update_layout(
                    title="Mean Expression Correlation Between Studies",
                    height=400,
                )
                st.plotly_chart(fig, use_container_width=True)
                
                # Assess
                min_corr = corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].min()
                if min_corr < 0.8:
                    st.warning(
                        f"Lowest inter-study correlation: {min_corr:.3f}. "
                        f"Strong batch effects detected. Batch correction recommended."
                    )
                elif min_corr < 0.95:
                    st.info(f"Lowest inter-study correlation: {min_corr:.3f}. Moderate batch effects.")
                else:
                    st.success(f"Lowest inter-study correlation: {min_corr:.3f}. Studies are well-correlated.")
            
            # Sample count balance
            st.markdown("**Sample balance across studies**")
            sps = pool.samples_per_study
            min_s, max_s = min(sps.values()), max(sps.values())
            imbalance = max_s / min_s if min_s > 0 else float('inf')
            
            if imbalance > 5:
                st.warning(f"Sample imbalance: {imbalance:.1f}x ({min_s} to {max_s} samples). This may bias results.")
            else:
                st.success(f"Sample balance OK ({min_s} to {max_s} samples per study).")
        
        with qc_tab3:
            # PCA colored by study
            st.markdown("**PCA of pooled samples (colored by study)**")
            
            with st.spinner("Computing PCA..."):
                from sklearn.decomposition import PCA
                from sklearn.preprocessing import StandardScaler
                
                # Log-transform and scale
                log_counts = np.log2(counts.T + 1)
                
                # Handle constant features
                var_mask = log_counts.var() > 0
                log_filtered = log_counts.loc[:, var_mask]
                
                scaler = StandardScaler()
                scaled = scaler.fit_transform(log_filtered)
                
                pca = PCA(n_components=min(3, scaled.shape[1]))
                pcs = pca.fit_transform(scaled)
                
                pca_df = pd.DataFrame({
                    'PC1': pcs[:, 0],
                    'PC2': pcs[:, 1],
                    'Study': pool.study_labels.loc[counts.columns].values,
                    'Sample': counts.columns,
                })
                
                if pcs.shape[1] >= 3:
                    pca_df['PC3'] = pcs[:, 2]
            
            fig = px.scatter(
                pca_df, x='PC1', y='PC2',
                color='Study',
                hover_data=['Sample'],
                title=f"PCA — PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%) vs PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)",
            )
            fig.update_layout(height=500)
            fig.update_traces(marker=dict(size=8))
            st.plotly_chart(fig, use_container_width=True)
            
            # Check if studies cluster separately (batch effect indicator)
            st.caption(
                "If studies form distinct clusters, batch effects are likely present. "
                "Check if batch correction resolved this."
            )
            
            # Variance explained
            var_explained = pca.explained_variance_ratio_[:3] * 100
            st.markdown(f"**Variance explained:** PC1={var_explained[0]:.1f}%, PC2={var_explained[1]:.1f}%"
                        + (f", PC3={var_explained[2]:.1f}%" if len(var_explained) >= 3 else ""))
        
        with qc_tab4:
            # Gene expression distributions
            st.markdown("**Gene expression distributions by study**")
            
            # Sample a subset for speed
            sample_genes = counts.index[:2000] if len(counts) > 2000 else counts.index
            
            fig = go.Figure()
            for study in pool.studies:
                study_samples = pool.get_study_samples(study)
                study_cols = [s for s in study_samples if s in counts.columns]
                if study_cols:
                    mean_expr = np.log2(counts.loc[sample_genes, study_cols].mean(axis=1) + 1)
                    fig.add_trace(go.Violin(
                        y=mean_expr,
                        name=study,
                        box_visible=True,
                        meanline_visible=True,
                    ))
            
            fig.update_layout(
                title="Log2 Mean Expression Distribution by Study",
                yaxis_title="Log2(mean count + 1)",
                height=450,
            )
            st.plotly_chart(fig, use_container_width=True)
        
        # Overall quality verdict
        st.markdown("---")
        st.markdown("### Quality Verdict")
        
        issues = []
        if fold_range > 5:
            issues.append("Large library size differences between studies")
        if len(study_means) >= 2 and min_corr < 0.8:
            issues.append("Strong batch effects detected")
        if zero_pct > 80:
            issues.append("Very high sparsity (>80% zeros)")
        if imbalance > 5:
            issues.append("Severe sample imbalance across studies")
        
        if not issues:
            st.success("Your pooled dataset looks great! You can confidently move to the next analysis step.")
        else:
            st.warning(
                f"Found {len(issues)} potential issue(s) to review:\n" +
                "\n".join(f"- {issue}" for issue in issues) +
                "\n\nDon't worry — you can re-pool with different settings (try enabling batch correction) "
                "or remove problematic studies from the selection."
            )

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
        
        col1, col2 = st.columns(2)
        
        with col1:
            # CSV export
            csv_data = ds.counts_df.to_csv()
            st.download_button(
                "Download Counts (CSV)",
                data=csv_data,
                file_name=f"{export_name}_counts.csv",
                mime="text/csv",
                use_container_width=True,
            )
        
        with col2:
            # Metadata export
            meta_csv = ds.metadata.to_csv()
            st.download_button(
                "Download Metadata (CSV)",
                data=meta_csv,
                file_name=f"{export_name}_metadata.csv",
                mime="text/csv",
                use_container_width=True,
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
