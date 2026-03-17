"""
RAPTOR Dashboard - Settings & Configuration
============================================

Configure dashboard preferences and RAPTOR workflow settings.

Author: Ayeh Bolouki
Version: 2.2.0
"""

import streamlit as st
import json
from pathlib import Path
from typing import Dict, Any
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from raptor.dashboard.components.sidebar import render_sidebar, init_session_state

# =============================================================================
# PAGE CONFIGURATION
# =============================================================================

st.set_page_config(
    page_title="RAPTOR - Settings",
    page_icon="⚙️",
    layout="wide"
)

# =============================================================================
# DEFAULT SETTINGS
# =============================================================================

DEFAULT_SETTINGS = {
    'general': {
        'default_output_dir': 'output/',
        'auto_save': True,
        'show_advanced_options': False
    },
    'thresholds': {
        'fdr_threshold': 0.05,
        'lfc_threshold': 1.0,
        'min_count': 10,
        'min_samples': 3
    },
    'quality_control': {
        'outlier_consensus_threshold': 3,
        'min_library_size': 1000000,
        'max_library_size_cv': 0.5,
        'batch_correction': True
    },
    'ensemble': {
        'default_method': 'brown',  # FIXED: brown is recommended
        'min_methods_consensus': 2,
        'direction_threshold': 0.8
    },
    'visualization': {
        'default_point_size': 5,
        'default_color_scheme': 'RdBu_r',
        'figure_width': 1200,
        'figure_height': 800,
        'dpi': 300
    },
    'reporting': {
        'include_plots': True,
        'include_statistics': True,
        'format': 'markdown'
    }
}

# =============================================================================
# INITIALIZE SESSION STATE
# =============================================================================

# Initialize RAPTOR dashboard session state
init_session_state()

# Initialize settings-specific state
if 'raptor_settings' not in st.session_state:
    st.session_state.raptor_settings = DEFAULT_SETTINGS.copy()
if 'settings_modified' not in st.session_state:
    st.session_state.settings_modified = False

# =============================================================================
# SETTINGS MANAGEMENT FUNCTIONS
# =============================================================================

def save_settings_to_file(settings: Dict[str, Any], filepath: Path):
    """Save settings to JSON file."""
    try:
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, 'w') as f:
            json.dump(settings, f, indent=2)
        return True
    except Exception as e:
        st.error(f"Error saving settings: {str(e)}")
        return False

def load_settings_from_file(filepath: Path) -> Dict[str, Any]:
    """Load settings from JSON file."""
    try:
        if filepath.exists():
            with open(filepath, 'r') as f:
                return json.load(f)
        else:
            return DEFAULT_SETTINGS.copy()
    except Exception as e:
        st.warning(f"Could not load settings: {str(e)}. Using defaults.")
        return DEFAULT_SETTINGS.copy()

def reset_to_defaults():
    """Reset all settings to defaults."""
    st.session_state.raptor_settings = DEFAULT_SETTINGS.copy()
    st.session_state.settings_modified = True

# =============================================================================
# MAIN UI
# =============================================================================

# Render sidebar
render_sidebar()

# Header
st.title("RAPTOR Settings & Configuration")
st.markdown("### Configure your RAPTOR dashboard preferences and analysis defaults")

# Help section
with st.expander("ℹ️ About Settings"):
    st.markdown("""
    **Purpose:** Configure default parameters and preferences for RAPTOR analysis.
    
    ## **Settings Categories**
    
    ### **General**
    - Output directories
    - Auto-save preferences
    - Interface options
    
    ### **Thresholds**
    - Default FDR threshold
    - Default Log2FC threshold
    - Minimum count/sample filters
    
    ### **Quality Control**
    - Outlier detection settings
    - Library size requirements
    - Batch correction options
    
    ### **Ensemble Analysis**
    - Default ensemble method
    - Consensus requirements
    - Direction threshold
    
    ### **Visualization**
    - Plot appearance
    - Figure dimensions
    - Export resolution
    
    ### **Reporting**
    - Report content
    - Default export format
    
    ## **Saving Settings**
    
    - Click "Save Settings" to persist changes
    - Settings saved to: `~/.raptor/settings.json`
    - Export/import settings for sharing
    - Reset to defaults anytime
    """)

st.markdown("---")

# Settings tabs
tabs = st.tabs([
    "🏠 General",
    "Thresholds",
    "Quality Control",
    "Ensemble",
    "Visualization",
    "Reporting",
    "About"
])

# =============================================================================
# TAB 1: GENERAL SETTINGS
# =============================================================================

with tabs[0]:
    st.header("🏠 General Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📁 File Paths")
        
        default_output = st.text_input(
            "Default Output Directory",
            value=st.session_state.raptor_settings['general']['default_output_dir'],
            help="Default location for saving analysis results"
        )
        if default_output != st.session_state.raptor_settings['general']['default_output_dir']:
            st.session_state.raptor_settings['general']['default_output_dir'] = default_output
            st.session_state.settings_modified = True
        
        auto_save = st.checkbox(
            "Auto-save results",
            value=st.session_state.raptor_settings['general']['auto_save'],
            help="Automatically save results after each analysis step"
        )
        if auto_save != st.session_state.raptor_settings['general']['auto_save']:
            st.session_state.raptor_settings['general']['auto_save'] = auto_save
            st.session_state.settings_modified = True
    
    with col2:
        st.subheader("Interface")
        
        show_advanced = st.checkbox(
            "Show advanced options",
            value=st.session_state.raptor_settings['general']['show_advanced_options'],
            help="Show advanced configuration options in analysis modules"
        )
        if show_advanced != st.session_state.raptor_settings['general']['show_advanced_options']:
            st.session_state.raptor_settings['general']['show_advanced_options'] = show_advanced
            st.session_state.settings_modified = True

# =============================================================================
# TAB 2: THRESHOLDS
# =============================================================================

with tabs[1]:
    st.header("Default Thresholds")
    
    st.info("These are the default thresholds used when you start a new analysis. You can always adjust them in individual modules.")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Statistical Thresholds")
        
        fdr_threshold = st.slider(
            "FDR Threshold",
            min_value=0.001,
            max_value=0.20,
            value=st.session_state.raptor_settings['thresholds']['fdr_threshold'],
            step=0.001,
            format="%.3f",
            help="Default adjusted p-value threshold (0.05 is standard)"
        )
        if fdr_threshold != st.session_state.raptor_settings['thresholds']['fdr_threshold']:
            st.session_state.raptor_settings['thresholds']['fdr_threshold'] = fdr_threshold
            st.session_state.settings_modified = True
        
        lfc_threshold = st.slider(
            "Log2 Fold Change Threshold",
            min_value=0.0,
            max_value=3.0,
            value=st.session_state.raptor_settings['thresholds']['lfc_threshold'],
            step=0.1,
            format="%.1f",
            help="Default log2 fold change threshold (1.0 = 2-fold change)"
        )
        if lfc_threshold != st.session_state.raptor_settings['thresholds']['lfc_threshold']:
            st.session_state.raptor_settings['thresholds']['lfc_threshold'] = lfc_threshold
            st.session_state.settings_modified = True
        
        fc = 2 ** lfc_threshold
        st.caption(f"Fold change: {fc:.2f}×")
    
    with col2:
        st.subheader("🔢 Filtering Thresholds")
        
        min_count = st.number_input(
            "Minimum Read Count",
            min_value=1,
            max_value=100,
            value=st.session_state.raptor_settings['thresholds']['min_count'],
            help="Minimum number of reads per gene"
        )
        if min_count != st.session_state.raptor_settings['thresholds']['min_count']:
            st.session_state.raptor_settings['thresholds']['min_count'] = min_count
            st.session_state.settings_modified = True
        
        min_samples = st.number_input(
            "Minimum Sample Count",
            min_value=1,
            max_value=10,
            value=st.session_state.raptor_settings['thresholds']['min_samples'],
            help="Minimum number of samples with min count"
        )
        if min_samples != st.session_state.raptor_settings['thresholds']['min_samples']:
            st.session_state.raptor_settings['thresholds']['min_samples'] = min_samples
            st.session_state.settings_modified = True

# =============================================================================
# TAB 3: QUALITY CONTROL
# =============================================================================

with tabs[2]:
    st.header("Quality Control Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Outlier Detection")
        
        outlier_threshold = st.slider(
            "Consensus Threshold",
            min_value=1,
            max_value=6,
            value=st.session_state.raptor_settings['quality_control']['outlier_consensus_threshold'],
            help="Number of methods that must agree to flag outlier"
        )
        if outlier_threshold != st.session_state.raptor_settings['quality_control']['outlier_consensus_threshold']:
            st.session_state.raptor_settings['quality_control']['outlier_consensus_threshold'] = outlier_threshold
            st.session_state.settings_modified = True
        
        st.caption(f"Outlier if detected by ≥{outlier_threshold} methods")
    
    with col2:
        st.subheader("Library Quality")
        
        min_lib_size = st.number_input(
            "Minimum Library Size",
            min_value=100000,
            max_value=10000000,
            value=st.session_state.raptor_settings['quality_control']['min_library_size'],
            step=100000,
            help="Minimum total reads per sample"
        )
        if min_lib_size != st.session_state.raptor_settings['quality_control']['min_library_size']:
            st.session_state.raptor_settings['quality_control']['min_library_size'] = min_lib_size
            st.session_state.settings_modified = True
        
        max_cv = st.slider(
            "Max Library Size CV",
            min_value=0.1,
            max_value=1.0,
            value=st.session_state.raptor_settings['quality_control']['max_library_size_cv'],
            step=0.05,
            format="%.2f",
            help="Maximum coefficient of variation for library sizes"
        )
        if max_cv != st.session_state.raptor_settings['quality_control']['max_library_size_cv']:
            st.session_state.raptor_settings['quality_control']['max_library_size_cv'] = max_cv
            st.session_state.settings_modified = True
    
    batch_correction = st.checkbox(
        "Enable Batch Correction Detection",
        value=st.session_state.raptor_settings['quality_control']['batch_correction'],
        help="Check for batch effects in quality control"
    )
    if batch_correction != st.session_state.raptor_settings['quality_control']['batch_correction']:
        st.session_state.raptor_settings['quality_control']['batch_correction'] = batch_correction
        st.session_state.settings_modified = True

# =============================================================================
# TAB 4: ENSEMBLE ANALYSIS
# =============================================================================

with tabs[3]:
    st.header("Ensemble Analysis Settings")
    
    st.info("Configure default settings for ensemble analysis (Module 9)")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Method Selection")
        
        # FIXED: Only show functional methods!
        default_method = st.selectbox(
            "Default Ensemble Method",
            options=['brown', 'fisher', 'rra'],
            index=['brown', 'fisher', 'rra'].index(
                st.session_state.raptor_settings['ensemble']['default_method']
            ),
            help="Default method for ensemble analysis"
        )
        if default_method != st.session_state.raptor_settings['ensemble']['default_method']:
            st.session_state.raptor_settings['ensemble']['default_method'] = default_method
            st.session_state.settings_modified = True
        
        st.caption("""
        **Brown's** - Recommended (accounts for correlation)  
        **Fisher's** - Independent methods only  
        **RRA** - Rank-based, robust to outliers
        """)
    
    with col2:
        st.subheader("Consensus Settings")
        
        min_consensus = st.slider(
            "Minimum Methods for Consensus",
            min_value=1,
            max_value=5,
            value=st.session_state.raptor_settings['ensemble']['min_methods_consensus'],
            help="Minimum number of methods that must agree"
        )
        if min_consensus != st.session_state.raptor_settings['ensemble']['min_methods_consensus']:
            st.session_state.raptor_settings['ensemble']['min_methods_consensus'] = min_consensus
            st.session_state.settings_modified = True
        
        direction_threshold = st.slider(
            "Direction Agreement Threshold",
            min_value=0.5,
            max_value=1.0,
            value=st.session_state.raptor_settings['ensemble']['direction_threshold'],
            step=0.05,
            format="%.2f",
            help="Fraction of methods that must agree on direction"
        )
        if direction_threshold != st.session_state.raptor_settings['ensemble']['direction_threshold']:
            st.session_state.raptor_settings['ensemble']['direction_threshold'] = direction_threshold
            st.session_state.settings_modified = True
    
    # FIXED: Only describe functional methods
    with st.expander("Ensemble Method Descriptions"):
        st.markdown("""
        ### **Available Ensemble Methods**
        
        **Brown's Method** (Recommended)
        - Combines p-values accounting for correlation
        - Best for: Typical RNA-seq (DESeq2 + edgeR + limma-voom)
        - Use when: Methods use same normalized count data
        - **This is the recommended default!**
        
        **Fisher's Method**
        - Combines p-values assuming independence
        - Best for: Truly independent methods
        - Use when: Methods use different algorithms (e.g., DESeq2 + Wilcoxon)
        
        **RRA (Robust Rank Aggregation)**
        - Ranks genes across methods, aggregates robustly
        - Best for: Finding top genes regardless of p-values
        - Use when: Methods have very different p-value distributions
        
        ---
        
        **Note:** Voting and Weighted ensemble methods are being upgraded 
        and will be available in a future release.
        """)

# =============================================================================
# TAB 5: VISUALIZATION
# =============================================================================

with tabs[4]:
    st.header("Visualization Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Plot Appearance")
        
        point_size = st.slider(
            "Default Point Size",
            min_value=2,
            max_value=10,
            value=st.session_state.raptor_settings['visualization']['default_point_size'],
            help="Size of points in scatter plots"
        )
        if point_size != st.session_state.raptor_settings['visualization']['default_point_size']:
            st.session_state.raptor_settings['visualization']['default_point_size'] = point_size
            st.session_state.settings_modified = True
        
        color_scheme = st.selectbox(
            "Default Color Scheme",
            options=['RdBu_r', 'viridis', 'plasma', 'coolwarm', 'seismic'],
            index=['RdBu_r', 'viridis', 'plasma', 'coolwarm', 'seismic'].index(
                st.session_state.raptor_settings['visualization']['default_color_scheme']
            ),
            help="Color scheme for heatmaps and plots"
        )
        if color_scheme != st.session_state.raptor_settings['visualization']['default_color_scheme']:
            st.session_state.raptor_settings['visualization']['default_color_scheme'] = color_scheme
            st.session_state.settings_modified = True
    
    with col2:
        st.subheader("📐 Export Settings")
        
        fig_width = st.number_input(
            "Figure Width (px)",
            min_value=600,
            max_value=2400,
            value=st.session_state.raptor_settings['visualization']['figure_width'],
            step=100
        )
        if fig_width != st.session_state.raptor_settings['visualization']['figure_width']:
            st.session_state.raptor_settings['visualization']['figure_width'] = fig_width
            st.session_state.settings_modified = True
        
        fig_height = st.number_input(
            "Figure Height (px)",
            min_value=400,
            max_value=1600,
            value=st.session_state.raptor_settings['visualization']['figure_height'],
            step=100
        )
        if fig_height != st.session_state.raptor_settings['visualization']['figure_height']:
            st.session_state.raptor_settings['visualization']['figure_height'] = fig_height
            st.session_state.settings_modified = True
        
        dpi = st.number_input(
            "DPI (Resolution)",
            min_value=72,
            max_value=600,
            value=st.session_state.raptor_settings['visualization']['dpi'],
            step=50,
            help="Dots per inch for exported figures"
        )
        if dpi != st.session_state.raptor_settings['visualization']['dpi']:
            st.session_state.raptor_settings['visualization']['dpi'] = dpi
            st.session_state.settings_modified = True

# =============================================================================
# TAB 6: REPORTING
# =============================================================================

with tabs[5]:
    st.header("Reporting Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📄 Report Content")
        
        include_plots = st.checkbox(
            "Include Plots in Reports",
            value=st.session_state.raptor_settings['reporting']['include_plots'],
            help="Embed plots in generated reports"
        )
        if include_plots != st.session_state.raptor_settings['reporting']['include_plots']:
            st.session_state.raptor_settings['reporting']['include_plots'] = include_plots
            st.session_state.settings_modified = True
        
        include_stats = st.checkbox(
            "Include Detailed Statistics",
            value=st.session_state.raptor_settings['reporting']['include_statistics'],
            help="Include comprehensive statistical tables"
        )
        if include_stats != st.session_state.raptor_settings['reporting']['include_statistics']:
            st.session_state.raptor_settings['reporting']['include_statistics'] = include_stats
            st.session_state.settings_modified = True
    
    with col2:
        st.subheader("Export Format")
        
        report_format = st.selectbox(
            "Default Report Format",
            options=['markdown', 'html', 'text'],
            index=['markdown', 'html', 'text'].index(
                st.session_state.raptor_settings['reporting']['format']
            ),
            help="Default format for exported reports"
        )
        if report_format != st.session_state.raptor_settings['reporting']['format']:
            st.session_state.raptor_settings['reporting']['format'] = report_format
            st.session_state.settings_modified = True

# =============================================================================
# TAB 7: ABOUT
# =============================================================================

with tabs[6]:
    st.header("About RAPTOR")
    
    st.markdown("""
    ## RAPTOR v2.2.0
    
    **RNA-seq Analysis Pipeline Testing and Optimization Resource**
    
    ### Author
    **Ayeh Bolouki**  
    Email: ayehbolouki1988@gmail.com
    
    ### Features
    
    RAPTOR provides a complete workflow for RNA-seq differential expression analysis:
    
    - **Quality Control** - 6-method outlier detection, batch effect analysis
    - **Data Profiling** - 32 features including BCV for pipeline selection
    - **Hybrid Recommendation** - Rule-based + ML pipeline selection
    - **DE Import** - Standardize results from DESeq2, edgeR, limma-voom, Wilcoxon
    - **Parameter Optimization** - Interactive threshold tuning
    - **Ensemble Analysis** - Combine multiple methods (Fisher's, Brown's, RRA)
    - **Reports** - Generate publication-ready documentation
    - **Settings** - Configurable defaults and preferences
    
    ### Documentation
    
    - **GitHub:** [https://github.com/AyehBlk/RAPTOR](https://github.com/AyehBlk/RAPTOR)
    - **Modules:** See module documentation files (MODULE*_*.md)
    - **CLI:** Run `raptor --help` for command-line usage
    
    ### Version History
    
    **v2.2.0** (Current - February 2026)
    - Complete dashboard interface
    - Hybrid ML+Rule recommendation system
    - Enhanced ensemble analysis (Fisher's, Brown's, RRA)
    - Interactive parameter optimization
    - Improved batch effect detection
    - 32-feature data profiler
    
    **v2.1.0**
    - ML-based pipeline recommendation
    - Parameter optimization module
    - Production pipeline runner
    
    **v2.0.0**
    - Complete refactor
    - Modular architecture
    - Dashboard interface
    
    ### License
    
    MIT License - See repository for details
    
    ### Citation
    
    If you use RAPTOR in your research, please cite:
    
    ```
    Bolouki, A. (2026). RAPTOR: RNA-seq Analysis Pipeline 
    Testing and Optimization Resource. v2.2.0.
    GitHub: https://github.com/AyehBlk/RAPTOR
    ```
    
    ### Support
    
    For issues, questions, or contributions:
    - **GitHub Issues:** [Report bugs](https://github.com/AyehBlk/RAPTOR/issues)
    - **Email:** ayehbolouki1988@gmail.com
    
    ---
    
    © 2026 Ayeh Bolouki. All rights reserved.
    """)

# =============================================================================
# SETTINGS MANAGEMENT (Bottom of page)
# =============================================================================

st.markdown("---")
st.markdown("## Settings Management")

col1, col2, col3, col4 = st.columns(4)

settings_file = Path.home() / '.raptor' / 'settings.json'

with col1:
    if st.button("Save Settings", use_container_width=True, type="primary"):
        if save_settings_to_file(st.session_state.raptor_settings, settings_file):
            st.success("Settings saved!")
            st.session_state.settings_modified = False
            st.rerun()

with col2:
    if st.button("📂 Load Settings", use_container_width=True):
        st.session_state.raptor_settings = load_settings_from_file(settings_file)
        st.session_state.settings_modified = False
        st.success("Settings loaded!")
        st.rerun()

with col3:
    if st.button("Reset to Defaults", use_container_width=True):
        reset_to_defaults()
        st.success("Reset to defaults!")
        st.rerun()

with col4:
    # Export settings
    settings_json = json.dumps(st.session_state.raptor_settings, indent=2)
    st.download_button(
        label="Export",
        data=settings_json,
        file_name="raptor_settings.json",
        mime="application/json",
        use_container_width=True
    )

# Import settings
st.markdown("### Import Settings")
uploaded_file = st.file_uploader(
    "Upload settings file (JSON)",
    type=['json'],
    help="Import previously exported settings"
)

if uploaded_file:
    try:
        imported_settings = json.load(uploaded_file)
        st.session_state.raptor_settings = imported_settings
        st.session_state.settings_modified = True
        st.success("Settings imported successfully!")
        st.rerun()
    except json.JSONDecodeError:
        st.error("Invalid settings file format")

# Status
st.markdown("---")
st.info(f"""
**Settings Location:** `{settings_file}`

Settings are automatically saved when you click "Save Settings" and persist across sessions.
""")

if st.session_state.settings_modified:
    st.warning("**You have unsaved changes.** Click 'Save Settings' above to persist your changes.")
else:
    st.success("All changes saved")

# Footer
st.markdown("---")
st.caption("**RAPTOR v2.2.0** | Settings | Configure your analysis preferences")
