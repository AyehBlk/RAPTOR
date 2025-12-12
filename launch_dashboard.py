#!/usr/bin/env python3
"""
RAPTOR Dashboard Launcher

Quick launcher for the interactive web dashboard.

Author: Ayeh Bolouki
Version: 2.1.0
"""

import subprocess
import sys
from pathlib import Path

def main():
    """Launch the RAPTOR dashboard."""
    # Look for dashboard in multiple possible locations
    base_path = Path(__file__).parent
    
    possible_paths = [
        base_path / "dashboard" / "app.py",           # Primary location
        base_path / "dashboard" / "dashboard.py",     # Alternative name
        base_path / "raptor" / "dashboard" / "app.py", # Inside raptor package
        base_path / "dashboard.py",                    # Same directory (fallback)
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
        print("\nMake sure the dashboard folder exists with app.py or dashboard.py")
        sys.exit(1)
    
    print(f"""
    ╔══════════════════════════════════════════════════════════════╗
    ║              🦖 Launching RAPTOR v2.1.0 Dashboard            ║
    ╚══════════════════════════════════════════════════════════════╝
    
    Dashboard: {dashboard_path.name}
    Location:  {dashboard_path.parent}
    
    The dashboard will open in your default web browser.
    
    Features available:
    • 🤖 ML-based pipeline recommendations
    • 📊 Real-time resource monitoring  
    • 🎯 Ensemble analysis
    • 📈 Benchmark comparisons
    • 📋 Quality assessment
    • 📄 Automated reporting
    
    Press Ctrl+C to stop the server.
    """)
    
    try:
        subprocess.run([
            sys.executable, "-m", "streamlit", "run",
            str(dashboard_path),
            "--server.headless", "true",
            "--browser.gatherUsageStats", "false"
        ])
    except KeyboardInterrupt:
        print("\n\n✅ Dashboard stopped. Thank you for using RAPTOR!")
    except FileNotFoundError:
        print("\n❌ Streamlit not installed!")
        print("\nInstall with:")
        print("  pip install streamlit plotly")
        print("\nOr install all RAPTOR dependencies:")
        print("  pip install -r requirements.txt")
        sys.exit(1)


if __name__ == "__main__":
    main()
