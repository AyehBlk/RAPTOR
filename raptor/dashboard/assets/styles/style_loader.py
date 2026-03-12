"""
RAPTOR Dashboard - Style Loader
================================

Helper functions to load custom CSS and styling for the dashboard.

Author: Ayeh Bolouki
Version: 2.2.0
Date: January 2026

Usage:
------
In your dashboard page:

```python
from assets.styles.style_loader import load_custom_css

# Load custom styles
load_custom_css()
```
"""

import streamlit as st
from pathlib import Path


def load_custom_css(css_file: str = "custom_styles.css") -> None:
    """
    Load custom CSS file into Streamlit app.
    
    Parameters
    ----------
    css_file : str
        Name of CSS file to load (default: custom_styles.css)
    
    Examples
    --------
    >>> load_custom_css()  # Load default styles
    >>> load_custom_css("dark_theme.css")  # Load dark theme
    """
    # Get the path to the CSS file
    css_path = Path(__file__).parent / css_file
    
    if css_path.exists():
        with open(css_path, 'r') as f:
            css_content = f.read()
        
        # Inject CSS into Streamlit
        st.markdown(f'<style>{css_content}</style>', unsafe_allow_html=True)
    else:
        st.warning(f"⚠️ Custom CSS file not found: {css_file}")


def load_logo(logo_file: str = "raptor_logo.png") -> str:
    """
    Load and return path to logo file.
    
    Parameters
    ----------
    logo_file : str
        Name of logo file (default: raptor_logo.png)
    
    Returns
    -------
    str
        Path to logo file, or None if not found
    
    Examples
    --------
    >>> logo_path = load_logo()
    >>> if logo_path:
    >>>     st.image(logo_path, width=200)
    """
    logo_path = Path(__file__).parent / logo_file
    
    if logo_path.exists():
        return str(logo_path)
    else:
        return None


def apply_raptor_branding() -> None:
    """
    Apply RAPTOR branding to the current page.
    
    This adds:
    - Custom title with 🦖 emoji
    - Subtitle with version info
    - Custom CSS styling
    
    Examples
    --------
    >>> apply_raptor_branding()
    >>> st.title("My Custom Page")  # Will be styled with RAPTOR theme
    """
    load_custom_css()
    
    # Add custom title styling
    st.markdown("""
        <style>
        .raptor-header {
            text-align: center;
            padding: 1rem 0;
            margin-bottom: 2rem;
        }
        .raptor-logo {
            font-size: 4rem;
        }
        .raptor-title {
            font-size: 2.5rem;
            font-weight: 700;
            color: #2E7D32;
            margin: 0.5rem 0;
        }
        .raptor-version {
            font-size: 1rem;
            color: #666;
        }
        </style>
    """, unsafe_allow_html=True)


def create_metric_card(label: str, value: str, delta: str = None, 
                       delta_color: str = "normal") -> None:
    """
    Create a styled metric card.
    
    Parameters
    ----------
    label : str
        Metric label
    value : str
        Metric value
    delta : str, optional
        Change value (e.g., "+10%")
    delta_color : str
        Color for delta ("normal", "inverse", "off")
    
    Examples
    --------
    >>> create_metric_card("Total Genes", "20,543", "+1,234", "normal")
    >>> create_metric_card("Quality Score", "87/100")
    """
    st.markdown(f"""
        <div class="module-card">
            <div class="stat-item">
                <div class="label">{label}</div>
                <div class="value">{value}</div>
                {f'<div style="color: {"green" if delta_color == "normal" else "red"}; font-size: 0.875rem;">{delta}</div>' if delta else ''}
            </div>
        </div>
    """, unsafe_allow_html=True)


def create_status_badge(text: str, status: str = "info") -> None:
    """
    Create a status badge.
    
    Parameters
    ----------
    text : str
        Badge text
    status : str
        Badge type: "success", "warning", "error", "info"
    
    Examples
    --------
    >>> create_status_badge("Analysis Complete", "success")
    >>> create_status_badge("Low Quality", "warning")
    """
    st.markdown(f"""
        <span class="status-badge {status}">{text}</span>
    """, unsafe_allow_html=True)


def create_info_box(title: str, content: str, box_type: str = "info") -> None:
    """
    Create a styled info box.
    
    Parameters
    ----------
    title : str
        Box title
    content : str
        Box content (can include HTML)
    box_type : str
        Type: "info", "success", "warning", "error"
    
    Examples
    --------
    >>> create_info_box(
    >>>     "Important Note",
    >>>     "Make sure to check quality metrics before proceeding.",
    >>>     "warning"
    >>> )
    """
    icon_map = {
        "info": "ℹ️",
        "success": "✅",
        "warning": "⚠️",
        "error": "❌"
    }
    
    icon = icon_map.get(box_type, "ℹ️")
    
    st.markdown(f"""
        <div class="module-card" style="border-left: 4px solid 
            {'#2196F3' if box_type == 'info' else 
             '#4CAF50' if box_type == 'success' else 
             '#FF9800' if box_type == 'warning' else 
             '#F44336'};">
            <h4>{icon} {title}</h4>
            <p>{content}</p>
        </div>
    """, unsafe_allow_html=True)


def add_footer() -> None:
    """
    Add RAPTOR footer to the page.
    
    Examples
    --------
    >>> add_footer()  # Add at the end of your page
    """
    st.markdown("""
        <div class="raptor-footer">
            <p><strong>🦖 RAPTOR v2.2.0</strong></p>
            <p>RNA-seq Analysis Pipeline Testing and Optimization Resource</p>
            <p>Author: Ayeh Bolouki | Email: ayehbolouki1988@gmail.com</p>
            <p style="font-size: 0.75rem; color: #999; margin-top: 0.5rem;">
                © 2026 Ayeh Bolouki. All rights reserved.
            </p>
        </div>
    """, unsafe_allow_html=True)


def get_color_palette() -> dict:
    """
    Get RAPTOR color palette.
    
    Returns
    -------
    dict
        Dictionary of color codes for RAPTOR branding
    
    Examples
    --------
    >>> colors = get_color_palette()
    >>> primary_color = colors['primary']
    """
    return {
        'primary': '#2E7D32',      # Dark green
        'primary_light': '#81C784', # Light green
        'primary_dark': '#1B5E20',  # Very dark green
        'accent': '#4CAF50',        # Accent green
        'success': '#4CAF50',       # Success green
        'info': '#2196F3',          # Info blue
        'warning': '#FF9800',       # Warning orange
        'error': '#F44336',         # Error red
        'text': '#262730',          # Text dark gray
        'text_light': '#666666',    # Light gray
        'background': '#FFFFFF',    # White
        'background_secondary': '#F0F2F6'  # Light gray
    }


# Example usage
if __name__ == "__main__":
    # This would be in your dashboard page
    import streamlit as st
    
    # Apply RAPTOR branding
    apply_raptor_branding()
    
    # Your page content
    st.title("Example Dashboard Page")
    
    # Create metric cards
    col1, col2, col3 = st.columns(3)
    with col1:
        create_metric_card("Total Genes", "20,543", "+1,234")
    with col2:
        create_metric_card("Quality Score", "87/100")
    with col3:
        create_metric_card("Outliers", "2", "-3")
    
    # Create info box
    create_info_box(
        "Getting Started",
        "Upload your count matrix to begin analysis.",
        "info"
    )
    
    # Add footer
    add_footer()
