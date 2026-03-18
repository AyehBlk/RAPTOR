"""
RAPTOR Dashboard Launcher

Usage (for pip users):
    python -m raptor.launch_dashboard

This module allows users who installed via pip to launch
the dashboard without needing the source repository.
"""
import os
import subprocess
import sys


def main():
    """Launch the RAPTOR Streamlit dashboard."""
    app_path = os.path.join(os.path.dirname(__file__), 'dashboard', 'app.py')

    if not os.path.exists(app_path):
        print("Error: Dashboard app.py not found at:", app_path)
        print("Make sure you installed with: pip install raptor-rnaseq[dashboard]")
        sys.exit(1)

    try:
        import streamlit  # noqa: F401
    except ImportError:
        print("Streamlit is not installed.")
        print("Install with: pip install raptor-rnaseq[dashboard]")
        sys.exit(1)

    print("Launching RAPTOR Dashboard...")
    print("URL: http://localhost:8501")
    print("Press Ctrl+C to stop.\n")

    subprocess.run([
        sys.executable, "-m", "streamlit", "run",
        app_path,
        "--server.headless", "true",
        "--browser.gatherUsageStats", "false",
    ])


if __name__ == "__main__":
    main()
