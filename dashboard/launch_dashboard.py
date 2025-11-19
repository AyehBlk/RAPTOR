#!/usr/bin/env python3

"""
RAPTOR Dashboard Launcher

Quick launcher for the interactive web dashboard.

Author: Ayeh Bolouki
"""

import subprocess
import sys
from pathlib import Path

def main():
    """Launch the dashboard."""
    dashboard_path = Path(__file__).parent / "dashboard.py"
    
    if not dashboard_path.exists():
        print("❌ Error: dashboard.py not found")
        sys.exit(1)
    
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║              🦖 Launching RAPTOR Dashboard                   ║
    ╚══════════════════════════════════════════════════════════════╝
    
    The dashboard will open in your default web browser.
    
    Features available:
    • 🤖 ML-based pipeline recommendations
    • 📊 Real-time resource monitoring
    • 🎯 Ensemble analysis
    • 📈 Benchmark comparisons
    
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
        print("\n\n✅ Dashboard stopped.")
    except FileNotFoundError:
        print("\n❌ Streamlit not installed!")
        print("Install with: pip install streamlit plotly")
        sys.exit(1)


if __name__ == "__main__":
    main()
