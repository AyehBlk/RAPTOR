#!/usr/bin/env python3
"""
RAPTOR Dashboard Launcher

Quick launcher for the interactive web dashboard.

Author: Ayeh Bolouki
Version: 2.1.1
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
    
    # Check for threshold optimizer
    ato_status = "✅ Available"
    try:
        from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
    except ImportError:
        try:
            from threshold_optimizer import AdaptiveThresholdOptimizer
        except ImportError:
            ato_status = "⚠️ Not installed"
    
    print(f"""
    ╔═══════════════════════════════════════════════════════════════╗
    ║              🦖 Launching RAPTOR v2.1.1 Dashboard             ║
    ╠═══════════════════════════════════════════════════════════════╣
    ║                                                               ║
    ║  🆕 NEW: Adaptive Threshold Optimizer (ATO)                   ║
    ║     Data-driven threshold selection for DE analysis           ║
    ║                                                               ║
    ╚═══════════════════════════════════════════════════════════════╝
    
    Dashboard: {dashboard_path.name}
    Location:  {dashboard_path.parent}
    
    Threshold Optimizer: {ato_status}
    
    The dashboard will open in your default web browser.
    
    Features available:
    • 🤖 ML-based pipeline recommendations
    • 🎯 Adaptive Threshold Optimizer (NEW!)
    • 📊 Real-time resource monitoring  
    • 🔬 Ensemble analysis
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
        print("  pip install raptor-rnaseq[dashboard]")
        sys.exit(1)


if __name__ == "__main__":
    main()
