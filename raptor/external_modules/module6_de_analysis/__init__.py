"""
RAPTOR Module 6 - Differential Expression Analysis (External R Module)

This module provides R-based differential expression analysis using
DESeq2, edgeR, or limma-voom.

Available Methods:
    - DESeq2: Recommended for moderate sample sizes (n > 6 per group)
    - edgeR: Good for small sample sizes (n < 6 per group)
    - limma-voom: Good for larger datasets with many samples

Example:
    >>> from raptor.external_modules.module6_de_analysis import DEAnalysisWrapper
    >>> de = DEAnalysisWrapper(method='deseq2')
    >>> result = de.run(counts='gene_counts.csv', metadata='metadata.csv')
"""

__version__ = "2.2.0"
__module__ = "M6"

from pathlib import Path
import warnings

# Module paths
MODULE_DIR = Path(__file__).parent
R_SCRIPTS_DIR = MODULE_DIR / "r_scripts"
PYTHON_WRAPPERS_DIR = MODULE_DIR / "python_wrappers"

# Available R scripts
DESEQ2_SCRIPT = R_SCRIPTS_DIR / "run_deseq2.R"
EDGER_SCRIPT = R_SCRIPTS_DIR / "run_edger.R"
LIMMA_SCRIPT = R_SCRIPTS_DIR / "run_limma.R"
INSTALL_SCRIPT = R_SCRIPTS_DIR / "install_packages.R"

# Try to import Python wrappers (may not be available in minimal installations)
try:
    from .python_wrappers import DEAnalysisWrapper
    WRAPPERS_AVAILABLE = True
except ImportError as e:
    WRAPPERS_AVAILABLE = False
    warnings.warn(
        f"Python wrappers not available: {e}\n"
        f"You can still use R scripts directly, but the Python API won't work.\n"
        f"To enable wrappers, ensure python_wrappers/ directory exists and contains proper code.",
        ImportWarning
    )
    DEAnalysisWrapper = None

# Try to import DEResult if it exists
try:
    from .python_wrappers import DEResult
except ImportError:
    DEResult = None


def _validate_r_scripts():
    """
    Validate that R scripts exist.
    
    Raises warnings if scripts are missing but doesn't fail,
    as this allows documentation building and development without R installed.
    """
    required_scripts = {
        'DESeq2': DESEQ2_SCRIPT,
        'edgeR': EDGER_SCRIPT,
        'limma': LIMMA_SCRIPT,
        'Package installer': INSTALL_SCRIPT
    }
    
    missing = {name: path for name, path in required_scripts.items() if not path.exists()}
    
    if missing:
        missing_str = '\n'.join([f"  - {name}: {path}" for name, path in missing.items()])
        warnings.warn(
            f"Missing R scripts:\n{missing_str}\n"
            f"Module 6 functionality will be limited.\n"
            f"Run installation script or check your RAPTOR installation.",
            UserWarning
        )


def get_available_methods():
    """
    Get list of available DE analysis methods.
    
    Returns:
        list: Available method names ['deseq2', 'edger', 'limma']
        
    Example:
        >>> from raptor.external_modules.module6_de_analysis import get_available_methods
        >>> methods = get_available_methods()
        >>> print(methods)
        ['deseq2', 'edger', 'limma']
    """
    methods = []
    if DESEQ2_SCRIPT.exists():
        methods.append('deseq2')
    if EDGER_SCRIPT.exists():
        methods.append('edger')
    if LIMMA_SCRIPT.exists():
        methods.append('limma')
    return methods


def check_r_installation():
    """
    Check if R is installed and accessible.
    
    Returns:
        bool: True if R is available, False otherwise
        
    Example:
        >>> from raptor.external_modules.module6_de_analysis import check_r_installation
        >>> if check_r_installation():
        ...     print("R is installed!")
    """
    import shutil
    return shutil.which('Rscript') is not None


def list_r_scripts():
    """
    List all available R scripts with their status.
    
    Returns:
        dict: Script names mapped to (path, exists) tuples
        
    Example:
        >>> from raptor.external_modules.module6_de_analysis import list_r_scripts
        >>> scripts = list_r_scripts()
        >>> for name, (path, exists) in scripts.items():
        ...     status = "✓" if exists else "✗"
        ...     print(f"{status} {name}: {path}")
    """
    scripts = {
        'deseq2': (DESEQ2_SCRIPT, DESEQ2_SCRIPT.exists()),
        'edger': (EDGER_SCRIPT, EDGER_SCRIPT.exists()),
        'limma': (LIMMA_SCRIPT, LIMMA_SCRIPT.exists()),
        'install': (INSTALL_SCRIPT, INSTALL_SCRIPT.exists())
    }
    return scripts


def get_script_path(method):
    """
    Get the path to a specific R script.
    
    Parameters:
        method (str): Method name ('deseq2', 'edger', 'limma', or 'install')
        
    Returns:
        Path: Path to the R script
        
    Raises:
        ValueError: If method is not recognized
        FileNotFoundError: If script doesn't exist
        
    Example:
        >>> from raptor.external_modules.module6_de_analysis import get_script_path
        >>> script = get_script_path('deseq2')
        >>> print(script)
        PosixPath('.../r_scripts/run_deseq2.R')
    """
    method = method.lower()
    
    script_map = {
        'deseq2': DESEQ2_SCRIPT,
        'edger': EDGER_SCRIPT,
        'limma': LIMMA_SCRIPT,
        'install': INSTALL_SCRIPT
    }
    
    if method not in script_map:
        raise ValueError(
            f"Unknown method: {method}. "
            f"Available methods: {', '.join(script_map.keys())}"
        )
    
    script_path = script_map[method]
    
    if not script_path.exists():
        raise FileNotFoundError(
            f"R script not found: {script_path}\n"
            f"Run the installation script or check your RAPTOR installation."
        )
    
    return script_path


# Validate R scripts on import (with warnings, not errors)
_validate_r_scripts()


# Export everything
__all__ = [
    # Version info
    '__version__',
    '__module__',
    
    # Paths
    'MODULE_DIR',
    'R_SCRIPTS_DIR',
    'PYTHON_WRAPPERS_DIR',
    
    # R script paths
    'DESEQ2_SCRIPT',
    'EDGER_SCRIPT',
    'LIMMA_SCRIPT',
    'INSTALL_SCRIPT',
    
    # Python wrappers (if available)
    'DEAnalysisWrapper',
    'DEResult',
    'WRAPPERS_AVAILABLE',
    
    # Helper functions
    'get_available_methods',
    'check_r_installation',
    'list_r_scripts',
    'get_script_path',
]
