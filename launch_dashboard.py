#!/usr/bin/env python3
"""
RAPTOR Dashboard Launcher

Quick launcher for the interactive web dashboard.

Author: Ayeh Bolouki
"""

import subprocess
import sys
from pathlib import Path


def _get_raptor_version() -> str:
    """Read the RAPTOR version from the installed package, falling
    back gracefully if it can't be imported. Single source of truth
    is raptor/__init__.py's __version__."""
    try:
        from raptor import __version__
        return __version__
    except ImportError:
        return "unknown"


def check_dependencies():
    """Check if required packages are installed."""
    required = {
        'streamlit': 'streamlit',
        'plotly': 'plotly',
        'pandas': 'pandas',
        'numpy': 'numpy',
    }
    
    missing = []
    for package, import_name in required.items():
        try:
            __import__(import_name)
        except ImportError:
            missing.append(package)
    
    return missing

def check_raptor_modules():
    """Check which RAPTOR modules are available."""
    modules_to_check = [
        ('quality_assessment', 'Module 2: Quality Assessment'),
        ('profiler', 'Module 3: Data Profiler'),
        ('recommender', 'Module 4: ML Recommender'),
        ('external_modules.acquisition', 'Module 6b: Data Acquisition'),
        ('de_import', 'Module 7: DE Import'),
        ('parameter_optimization', 'Module 8: Parameter Optimization'),
        ('ensemble', 'Module 9: Ensemble Analysis'),
    ]
    
    available = []
    missing = []
    
    for module_name, description in modules_to_check:
        try:
            __import__(f'raptor.{module_name}')
            available.append(description)
        except ImportError:
            missing.append(description)
    
    return available, missing

def main():
    """Launch the RAPTOR dashboard."""
    
    # Check dependencies first
    missing = check_dependencies()
    if missing:
        print(f"""
    ❌ Missing required packages: {', '.join(missing)}
    
    Install with:
      pip install {' '.join(missing)}
    
    Or install all dashboard dependencies:
      pip install streamlit plotly pandas numpy
        """)
        sys.exit(1)
    
    # Look for dashboard in raptor/dashboard/
    base_path = Path(__file__).parent
    
    possible_paths = [
        base_path / "raptor" / "dashboard" / "app.py",     # Primary location
        base_path / "dashboard" / "app.py",                # Alternative
        base_path / "raptor" / "dashboard" / "dashboard.py", # Alternative name
    ]
    
    dashboard_path = None
    for path in possible_paths:
        if path.exists():
            dashboard_path = path
            break
    
    if dashboard_path is None:
        print("❌ Error: Dashboard not found!")
        print("\nSearched in:")
        for path in possible_paths:
            print(f"  • {path}")
        print("\nMake sure the dashboard folder exists at raptor/dashboard/")
        print("\nExpected structure:")
        print("  RAPTOR/")
        print("  ├── launch_dashboard.py  (this file)")
        print("  └── raptor/")
        print("      └── dashboard/")
        print("          ├── app.py")
        print("          └── pages/")
        sys.exit(1)
    
    # Check if core modules are available
    available_modules, missing_modules = check_raptor_modules()
    
    if available_modules:
        if missing_modules:
            core_status = f"✅ {len(available_modules)} available"
            if len(missing_modules) <= 3:
                core_status += f" | ⚠️ Missing: {', '.join([m.split(':')[0] for m in missing_modules])}"
        else:
            core_status = f"✅ All {len(available_modules)} modules available"
    else:
        core_status = "⚠️ No RAPTOR modules found (dashboard will work but with limited functionality)"
    
    # Display launch information
    raptor_version = _get_raptor_version()
    print(f"""
    ╔═══════════════════════════════════════════════════════════════╗
    ║   🦖 Launching RAPTOR v{raptor_version} — Dashboard
    ║      Professional RNA-seq Analysis Interface
    ╚═══════════════════════════════════════════════════════════════╝
    
    Dashboard Location: {dashboard_path.parent}
    Core Modules: {core_status}
    
    NEW in v{raptor_version}:
    ────────────────────────────────────────────────────────────────
    - Data Acquisition: Search & download from GEO, SRA, TCGA
    - Pool datasets with batch correction (ComBat, quantile)
    - Gene ID conversion (Ensembl/Symbol/Entrez via MyGene.info)
    - Quality checks on pooled data (PCA, correlations, RLE)
    - SRA run tables with FASTQ download scripts
    - TCGA multi-omic support (miRNA, methylation, CNV, RPPA)
    - Biomarker discovery (Module 10) with kneedle panel-size
      selection and consensus pinning across CV folds
    
    📋 Available Features:
    ────────────────────────────────────────────────────────────────
    • 🏠 Home - Overview & quick start
    • 📡 Data Acquisition - GEO/SRA/TCGA search & download (Module 6b)
    • ✅ Quality Assessment - QC metrics (Module 2)
    • 📊 Data Profiler - Comprehensive profiling (Module 3)
    • 🤖 ML Recommender - Pipeline recommendations (Module 4)
    • 📥 Import DE - Import DE results (Module 7)
    • ⚙️ Parameter Optimization - Threshold optimization (Module 8)
    • 🧬 Ensemble Analysis - Consensus genes (Module 9)
    • 🧬 Biomarker Discovery - Panel + clinical metrics (Module 10)
    • 📊 Visualizations - Interactive plots
    • 📋 Reports - Publication-ready exports
    • ⚙️ Settings - Dashboard configuration
    
    💡 Note: Modules 1 & 5 (quantification) available via CLI
          Run Salmon/Kallisto before using dashboard
    
    🌐 The dashboard will open in your default web browser.
    📍 URL: http://localhost:8501
    
    ⌨️  Press Ctrl+C to stop the server.
    ═══════════════════════════════════════════════════════════════
    """)
    
    try:
        # Launch Streamlit
        subprocess.run([
            sys.executable, "-m", "streamlit", "run",
            str(dashboard_path),
            "--server.headless", "true",
            "--browser.gatherUsageStats", "false",
            "--theme.primaryColor", "#2E7D32",
            "--theme.backgroundColor", "#FFFFFF",
            "--theme.secondaryBackgroundColor", "#F0F2F6"
        ])
    except KeyboardInterrupt:
        print("\n\n✅ Dashboard stopped. Thank you for using RAPTOR!")
    except FileNotFoundError:
        print("\n❌ Streamlit not installed!")
        print("\nInstall with:")
        print("  pip install streamlit plotly pandas numpy")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error launching dashboard: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()