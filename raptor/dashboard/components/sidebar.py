"""
Unified Sidebar Component

Reusable sidebar that appears on all dashboard pages.
Shows workflow progress, module status, and quick stats.

Author: Ayeh Bolouki
"""

import streamlit as st
from pathlib import Path

def render_sidebar():
    """
    Render the unified sidebar with workflow tracking and navigation.
    
    This should be called at the top of every dashboard page.
    """
    
    # Logo (if exists)
    logo_path = Path(__file__).parent.parent / "assets" / "images" / "logo.png"
    if logo_path.exists():
        st.sidebar.image(str(logo_path), width=200)
    else:
        st.sidebar.markdown("# 🦖 RAPTOR")
    
    st.sidebar.markdown("### v2.2.0")
    st.sidebar.markdown("---")
    
    # Workflow Progress
    st.sidebar.markdown("### 📊 Workflow Progress")
    
    # Calculate progress from session state
    progress = calculate_workflow_progress()
    st.sidebar.progress(progress)
    st.sidebar.caption(f"{int(progress*100)}% Complete")
    
    st.sidebar.markdown("---")
    
    # Module Status (Dashboard modules only - M2, M3, M4, M7, M8, M9)
    st.sidebar.markdown("### 🔄 Module Status")
    
    modules = {
        "M2: Quality": get_module_status(2),
        "M3: Profiler": get_module_status(3),
        "M4: Recommender": get_module_status(4),
        "M7: Import DE": get_module_status(7),
        "M8: Optimization": get_module_status(8),
        "M9: Ensemble": get_module_status(9),
    }
    
    for module, status in modules.items():
        if status == "complete":
            icon = "✅"
            color = "green"
        elif status == "running":
            icon = "⏳"
            color = "orange"
        elif status == "error":
            icon = "❌"
            color = "red"
        else:  # not started
            icon = "⭕"
            color = "gray"
        
        st.sidebar.markdown(f"{icon} {module}")
    
    st.sidebar.markdown("---")
    
    # Quick Stats
    st.sidebar.markdown("### 📈 Quick Stats")
    
    n_samples = st.session_state.get('n_samples', 0)
    n_genes = st.session_state.get('n_genes', 0)
    n_de_genes = st.session_state.get('n_de_genes', 0)
    
    st.sidebar.metric("Samples", n_samples)
    st.sidebar.metric("Total Genes", n_genes)
    st.sidebar.metric("DE Genes", n_de_genes)
    
    st.sidebar.markdown("---")
    
    # Help Section
    st.sidebar.markdown("### 💡 Need Help?")
    st.sidebar.info("Click ℹ️ icons for inline help on each page")
    
    if st.sidebar.button("📖 Documentation"):
        st.sidebar.markdown("[View Docs](https://github.com/AyehBlk/RAPTOR)")
    
    if st.sidebar.button("💬 Report Issue"):
        st.sidebar.markdown("[GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)")


def calculate_workflow_progress():
    """
    Calculate overall workflow progress based on completed modules.
    
    Dashboard tracks 6 modules: M2, M3, M4, M7, M8, M9
    (M1 and M5 are CLI-only, not in dashboard)
    
    Returns:
        float: Progress value between 0.0 and 1.0
    """
    completed_modules = []
    
    # Check which dashboard modules are complete
    if st.session_state.get('m2_complete', False):
        completed_modules.append(2)
    if st.session_state.get('m3_complete', False):
        completed_modules.append(3)
    if st.session_state.get('m4_complete', False):
        completed_modules.append(4)
    if st.session_state.get('m7_complete', False):
        completed_modules.append(7)
    if st.session_state.get('m8_complete', False):
        completed_modules.append(8)
    if st.session_state.get('m9_complete', False):
        completed_modules.append(9)
    
    # Calculate progress (6 dashboard modules total)
    total_modules = 6
    progress = len(completed_modules) / total_modules
    
    return progress


def get_module_status(module_num):
    """
    Get the status of a specific module.
    
    Args:
        module_num (int): Module number (1-9)
    
    Returns:
        str: Status ('not_started', 'running', 'complete', 'error')
    """
    # Check session state for module status
    if st.session_state.get(f'm{module_num}_error', False):
        return 'error'
    elif st.session_state.get(f'm{module_num}_running', False):
        return 'running'
    elif st.session_state.get(f'm{module_num}_complete', False):
        return 'complete'
    else:
        return 'not_started'


def init_session_state():
    """
    Initialize session state variables for the dashboard.
    
    Call this at the start of each page to ensure all variables exist.
    
    Note: Initializes all module flags (1-9) for compatibility,
    but dashboard only displays M2, M3, M4, M7, M8, M9.
    """
    # Module completion flags (initialize all for compatibility)
    for i in [1, 2, 3, 4, 5, 7, 8, 9]:
        if f'm{i}_complete' not in st.session_state:
            st.session_state[f'm{i}_complete'] = False
        if f'm{i}_running' not in st.session_state:
            st.session_state[f'm{i}_running'] = False
        if f'm{i}_error' not in st.session_state:
            st.session_state[f'm{i}_error'] = False
    
    # Data storage
    if 'n_samples' not in st.session_state:
        st.session_state['n_samples'] = 0
    if 'n_genes' not in st.session_state:
        st.session_state['n_genes'] = 0
    if 'n_de_genes' not in st.session_state:
        st.session_state['n_de_genes'] = 0
    
    # Module-specific data
    if 'm7_results' not in st.session_state:
        st.session_state['m7_results'] = None
    if 'm9_result' not in st.session_state:
        st.session_state['m9_result'] = None
