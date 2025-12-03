#!/usr/bin/env python3

"""
RAPTOR v2.1.0 - Master Installation & Verification Script

This script will:
1. Check Python version
2. Install all dependencies
3. Run comprehensive tests
4. Optionally generate training data
5. Launch the dashboard

Author: Ayeh Bolouki
Version: 2.1.0
"""

import sys
import subprocess
import os
from pathlib import Path

def print_header(text):
    """Print formatted header."""
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70 + "\n")

def print_success(text):
    """Print success message."""
    print(f"✅ {text}")

def print_error(text):
    """Print error message."""
    print(f"❌ {text}")

def print_info(text):
    """Print info message."""
    print(f"ℹ️  {text}")

def check_python_version():
    """Check if Python version is adequate."""
    print_header("Checking Python Version")
    
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print_error("Python 3.8 or higher is required!")
        print_info("Please upgrade Python: https://www.python.org/downloads/")
        return False
    
    print_success(f"Python {version.major}.{version.minor} is compatible")
    return True

def install_dependencies():
    """Install required packages."""
    print_header("Installing Dependencies")
    
    # Try different requirements files
    requirements_files = ["requirements.txt", "requirements_ml.txt"]
    requirements_file = None
    
    for rf in requirements_files:
        if Path(rf).exists():
            requirements_file = Path(rf)
            break
    
    if requirements_file is None:
        print_error("No requirements file found!")
        print_info("Expected: requirements.txt or requirements_ml.txt")
        return False
    
    print_info(f"Installing packages from {requirements_file}...")
    print_info("This may take a few minutes...\n")
    
    try:
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-r", str(requirements_file), "--quiet"],
            check=True
        )
        print_success("All dependencies installed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to install dependencies: {e}")
        print_info(f"Try manually: pip install -r {requirements_file}")
        return False

def verify_installation():
    """Verify all modules can be imported."""
    print_header("Verifying Installation")
    
    modules = [
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("sklearn", "scikit-learn"),
        ("scipy", "scipy"),
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("joblib", "joblib"),
        ("streamlit", "streamlit"),
        ("plotly", "plotly"),
        ("psutil", "psutil"),
        ("statsmodels", "statsmodels"),
        ("tqdm", "tqdm"),
        ("yaml", "pyyaml"),
    ]
    
    all_ok = True
    
    for module_name, display_name in modules:
        try:
            __import__(module_name)
            print_success(f"{display_name} installed")
        except ImportError:
            print_error(f"{display_name} not found")
            all_ok = False
    
    return all_ok

def run_tests():
    """Run the test suite."""
    print_header("Running Test Suite")
    
    test_files = ["test_ml_system.py", "tests/test_ml_system.py"]
    test_file = None
    
    for tf in test_files:
        if Path(tf).exists():
            test_file = Path(tf)
            break
    
    if test_file is None:
        print_info("Test file not found, skipping tests")
        return True
    
    print_info("Running comprehensive tests...")
    print_info("This will test all components of the system.\n")
    
    try:
        result = subprocess.run(
            [sys.executable, str(test_file)],
            capture_output=False,
            text=True
        )
        
        if result.returncode == 0:
            print_success("All tests passed!")
            return True
        else:
            print_error("Some tests failed. Check output above.")
            return False
    except Exception as e:
        print_error(f"Error running tests: {e}")
        return False

def generate_training_data():
    """Generate synthetic training data."""
    print_header("Training Data Generation")
    
    choice = input("Generate synthetic training data? (recommended for first use) [Y/n]: ").strip().lower()
    
    if choice in ['', 'y', 'yes']:
        workflow_files = ["example_ml_workflow.py", "examples/example_ml_workflow.py"]
        workflow_file = None
        
        for wf in workflow_files:
            if Path(wf).exists():
                workflow_file = Path(wf)
                break
        
        if workflow_file is None:
            print_info("Workflow file not found, skipping data generation")
            return True
        
        print_info("Generating 200 synthetic datasets...")
        print_info("This will take 2-5 minutes...\n")
        
        try:
            subprocess.run(
                [sys.executable, str(workflow_file), "--n-datasets", "200"],
                check=True
            )
            print_success("Training data generated successfully!")
            print_success("Model trained and saved to models/ directory")
            return True
        except subprocess.CalledProcessError as e:
            print_error(f"Failed to generate data: {e}")
            return False
    else:
        print_info("Skipping training data generation")
        print_info("You can generate it later with: python example_ml_workflow.py")
        return True

def launch_dashboard():
    """Launch the interactive dashboard."""
    print_header("Launch Dashboard")
    
    choice = input("Launch the interactive dashboard now? [Y/n]: ").strip().lower()
    
    if choice in ['', 'y', 'yes']:
        launcher_files = ["launch_dashboard.py", "raptor/launch_dashboard.py"]
        dashboard_files = ["dashboard.py", "raptor/dashboard.py"]
        
        launcher_file = None
        for lf in launcher_files:
            if Path(lf).exists():
                launcher_file = Path(lf)
                break
        
        dashboard_file = None
        for df in dashboard_files:
            if Path(df).exists():
                dashboard_file = Path(df)
                break
        
        if launcher_file:
            print_info("Launching dashboard...")
            try:
                subprocess.run([sys.executable, str(launcher_file)], check=False)
            except KeyboardInterrupt:
                print("\n")
                print_success("Dashboard stopped")
        elif dashboard_file:
            print_info("Launching dashboard...")
            print_info("The dashboard will open in your browser at http://localhost:8501")
            print_info("Press Ctrl+C to stop the dashboard\n")
            
            try:
                subprocess.run(
                    [sys.executable, "-m", "streamlit", "run", str(dashboard_file)],
                    check=False
                )
            except KeyboardInterrupt:
                print("\n")
                print_success("Dashboard stopped")
        else:
            print_error("Dashboard files not found!")
            return
    else:
        print_info("You can launch it later with: python launch_dashboard.py")

def print_summary():
    """Print installation summary."""
    print_header("🎉 Installation Complete!")
    
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║         RAPTOR v2.1.0 - Ready to Use! 🦖                             ║
╚══════════════════════════════════════════════════════════════════════╝

What's available:

🎨 Interactive Dashboard:
   python launch_dashboard.py
   → Web-based interface at http://localhost:8501

🤖 ML Recommendations:
   python raptor_ml_cli.py profile --counts data.csv --use-ml

📊 Quality Assessment:
   from raptor.data_quality_assessment import DataQualityAssessor

🎯 Ensemble Analysis:
   from raptor.ensemble_analysis import EnsembleAnalyzer

⚡ Resource Monitoring:
   from raptor.resource_monitoring import ResourceMonitor

📚 Documentation:
   - README.md              - Quick start guide
   - COMPLETE_README.md     - Master guide
   - QUICK_START.md         - Get running in 5 minutes
   - DASHBOARD_GUIDE.md     - Web interface guide

🧪 Test System:
   python test_ml_system.py

🔧 Support:
   Email: ayehbolouki1988@gmail.com
   GitHub: https://github.com/AyehBlk/RAPTOR

🚀 Quick Start:
   1. Launch dashboard: python launch_dashboard.py
   2. Upload data or use sample
   3. Get AI-powered recommendation
   4. Enjoy! 🎉

""")

def main():
    """Main installation workflow."""
    print("""
    ╔══════════════════════════════════════════════════════════════════╗
    ║              🦖 RAPTOR v2.1.0 Installer                          ║
    ║                                                                  ║
    ║     Machine Learning-Powered RNA-seq Analysis System             ║
    ║              with Interactive Dashboard                          ║
    ║                                                                  ║
    ║                   Version 2.1.0                                  ║
    ╚══════════════════════════════════════════════════════════════════╝
    """)
    
    # Step 1: Check Python
    if not check_python_version():
        sys.exit(1)
    
    # Step 2: Install dependencies
    if not install_dependencies():
        print_error("Installation failed at dependency stage")
        sys.exit(1)
    
    # Step 3: Verify installation
    if not verify_installation():
        print_error("Some dependencies are missing")
        print_info("Try: pip install -r requirements.txt")
        sys.exit(1)
    
    # Step 4: Run tests
    if not run_tests():
        print_error("Tests failed")
        print_info("You can continue, but some features may not work properly")
        choice = input("Continue anyway? [y/N]: ").strip().lower()
        if choice not in ['y', 'yes']:
            sys.exit(1)
    
    # Step 5: Generate training data (optional)
    generate_training_data()
    
    # Step 6: Print summary
    print_summary()
    
    # Step 7: Launch dashboard (optional)
    launch_dashboard()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️  Installation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
