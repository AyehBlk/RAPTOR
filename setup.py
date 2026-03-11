"""
RAPTOR v2.2.0 - RNA-seq Analysis Pipeline Testing and Optimization Resource
Setup Configuration

This setup.py enables RAPTOR to be installed as a Python package with full
CLI support and automatic dependency management.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
Repository: https://github.com/AyehBlk/RAPTOR
"""

from setuptools import setup, find_packages
from pathlib import Path

# =============================================================================
# Read README for long description
# =============================================================================

this_directory = Path(__file__).parent
long_description = ""

readme_file = this_directory / "README.md"
if readme_file.exists():
    long_description = readme_file.read_text(encoding='utf-8')

# =============================================================================
# Version Information
# =============================================================================

VERSION = "2.2.1"
PYTHON_REQUIRES = ">=3.8"

# =============================================================================
# Core Dependencies
# =============================================================================

INSTALL_REQUIRES = [
    # Scientific Computing Core
    "numpy>=1.19.0",
    "pandas>=1.1.0",
    "scipy>=1.5.0",
    
    # Visualization
    "plotly>=5.0.0",
    "matplotlib>=3.3.0",
    "seaborn>=0.11.0",
    
    # Statistics and Machine Learning
    "scikit-learn>=0.24.0",
    "statsmodels>=0.12.0",
    
    # CLI Framework
    "click>=7.0",
    
    # Configuration Management
    "pyyaml>=5.3.0",
    "toml>=0.10.0",
    
    # Utilities
    "tqdm>=4.50.0",        # Progress bars
    "colorama>=0.4.0",     # Colored terminal output
    "joblib>=1.0.0",       # Parallel processing
]

# =============================================================================
# Optional Dependencies
# =============================================================================

DASHBOARD_REQUIRES = [
    "streamlit>=1.20.0",
    "altair>=4.0.0",
]

DEV_REQUIRES = [
    "pytest>=6.0.0",
    "pytest-cov>=2.10.0",
    "pytest-mock>=3.6.0",
    "black>=21.0",
    "flake8>=3.8.0",
    "isort>=5.0.0",
    "mypy>=0.900",
]

DOCS_REQUIRES = [
    "sphinx>=3.5.0",
    "sphinx-rtd-theme>=0.5.0",
    "myst-parser>=0.15.0",
    "sphinx-click>=3.0.0",
]

# Combined extras
EXTRAS_REQUIRE = {
    "dashboard": DASHBOARD_REQUIRES,
    "dev": DEV_REQUIRES,
    "docs": DOCS_REQUIRES,
    "all": DASHBOARD_REQUIRES + DEV_REQUIRES + DOCS_REQUIRES,
}

# =============================================================================
# Package Data (Non-Python Files to Include)
# =============================================================================

PACKAGE_DATA = {
    "raptor": [
        # R scripts for differential expression (Module 6)
        "external_modules/module6_de_analysis/r_scripts/*.R",
        
        # Pipeline configuration files
        "pipelines/*/config.yaml",
        
        # Shell scripts for pipelines
        "pipelines/*/*.sh",
        
        # Documentation files in submodules
        "external_modules/module6_de_analysis/documentation/*.md",
        "external_modules/module6_de_analysis/*.md",
        "pipelines/star_salmon/*.md",
        
        # Dashboard assets
        "dashboard/assets/*.html",
        "dashboard/assets/styles/*.css",
        "dashboard/assets/styles/*.png",
        "dashboard/assets/styles/*.PNG",
        
        # Any other data files
        "*.yaml",
        "*.toml",
    ],
}

# =============================================================================
# Entry Points (CLI Commands)
# =============================================================================

ENTRY_POINTS = {
    "console_scripts": [
        # Main RAPTOR CLI
        "raptor=raptor.cli:main",
    ],
}

# =============================================================================
# PyPI Classifiers
# =============================================================================

CLASSIFIERS = [
    # Development Status
    "Development Status :: 4 - Beta",
    
    # Intended Audience
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    "Intended Audience :: Healthcare Industry",
    
    # Topic
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Information Analysis",
    "Topic :: Scientific/Engineering :: Medical Science Apps.",
    
    # License
    "License :: OSI Approved :: MIT License",
    
    # Python Versions
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    
    # Other Languages
    "Programming Language :: R",
    
    # Operating Systems
    "Operating System :: OS Independent",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows",
    
    # Additional
    "Natural Language :: English",
    "Environment :: Console",
    "Environment :: Web Environment",
]

# =============================================================================
# Keywords for PyPI Search
# =============================================================================

KEYWORDS = [
    "bioinformatics",
    "RNA-seq",
    "transcriptomics",
    "differential expression",
    "genomics",
    "data analysis",
    "machine learning",
    "pipeline",
    "ensemble methods",
    "DESeq2",
    "edgeR",
    "limma",
    "quality control",
    "parameter optimization",
    "biomarker discovery",
    "kallisto",
    "salmon",
    "STAR",
]

# =============================================================================
# Main Setup Configuration
# =============================================================================

setup(
    # -------------------------------------------------------------------------
    # Basic Information
    # -------------------------------------------------------------------------
    name="raptor-rnaseq",
    version=VERSION,
    description="RNA-seq Analysis Pipeline Testing and Optimization Resource",
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    # -------------------------------------------------------------------------
    # Author Information
    # -------------------------------------------------------------------------
    author="Ayeh Bolouki",
    author_email="ayehbolouki1988@gmail.com",
    
    # -------------------------------------------------------------------------
    # URLs
    # -------------------------------------------------------------------------
    url="https://github.com/AyehBlk/RAPTOR",
    project_urls={
        "Documentation": "https://github.com/AyehBlk/RAPTOR/tree/main/docs",
        "Source": "https://github.com/AyehBlk/RAPTOR",
        "Tracker": "https://github.com/AyehBlk/RAPTOR/issues",
    },
    
    # -------------------------------------------------------------------------
    # License
    # -------------------------------------------------------------------------
    license="MIT",
    
    # -------------------------------------------------------------------------
    # Package Discovery
    # -------------------------------------------------------------------------
    packages=find_packages(exclude=["tests", "tests.*", "docs", "examples"]),
    
    # -------------------------------------------------------------------------
    # Package Data
    # -------------------------------------------------------------------------
    package_data=PACKAGE_DATA,
    include_package_data=True,
    
    # -------------------------------------------------------------------------
    # Dependencies
    # -------------------------------------------------------------------------
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    
    # -------------------------------------------------------------------------
    # Entry Points (CLI)
    # -------------------------------------------------------------------------
    entry_points=ENTRY_POINTS,
    
    # -------------------------------------------------------------------------
    # Classification
    # -------------------------------------------------------------------------
    classifiers=CLASSIFIERS,
    keywords=", ".join(KEYWORDS),
    
    # -------------------------------------------------------------------------
    # Additional Options
    # -------------------------------------------------------------------------
    zip_safe=False,
    platforms=["any"],
)

# =============================================================================
# Post-Install Message
# =============================================================================

print("""
╔════════════════════════════════════════════════════════════════════════╗
║                                                                        ║
║  🦖 RAPTOR v2.2.0 Installation Complete!                             ║
║                                                                        ║
║  RNA-seq Analysis Pipeline Testing and Optimization Resource          ║
║                                                                        ║
╚════════════════════════════════════════════════════════════════════════╝

✅ Core package installed successfully!

📋 Quick Start:
   • CLI: raptor --help
   • Dashboard: streamlit run raptor/dashboard/app.py
   • Python API: from raptor import ensemble_fisher

📚 Documentation: https://github.com/AyehBlk/RAPTOR/tree/main/docs

🔬 Available Modules:
   ✓ Module 2: Quality Assessment
   ✓ Module 3: Data Profiler
   ✓ Module 4: Pipeline Recommender (Rule-based + ML)
   ✓ Module 5: Production Pipelines (Kallisto, Salmon, STAR, etc.)
   ✓ Module 6: Differential Expression (DESeq2, edgeR, limma)
   ✓ Module 7: DE Import & Comparison
   ✓ Module 8: Parameter Optimization
   ✓ Module 9: Ensemble Analysis

⚙️  Pipeline Tools Available:
   • Kallisto (quick & full)
   • Salmon (quick & full)
   • STAR + FeatureCounts
   • STAR + RSEM
   • STAR + Salmon
   • HISAT2 + FeatureCounts

⚠️  External Dependencies Required:
   
   R Packages (for Module 6 - DE Analysis):
   • DESeq2, edgeR, limma
   • Install with: Rscript raptor/external_modules/module6_de_analysis/r_scripts/install_packages.R
   
   Alignment Tools (for Module 5 - Pipelines):
   • STAR, HISAT2 (alignment)
   • Kallisto, Salmon (pseudoalignment)
   • Samtools, FeatureCounts, RSEM (counting)

💡 Test Installation:
   python -c "import raptor; print(raptor.get_info())"
   raptor --version

📦 Optional Features:
   • Dashboard: pip install .[dashboard]
   • Development: pip install .[dev]
   • Documentation: pip install .[docs]
   • Everything: pip install .[all]

🆘 Support:
   • Email: ayehbolouki1988@gmail.com
   • GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues

Happy analyzing! 🧬
""")
