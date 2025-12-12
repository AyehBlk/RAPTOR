"""
Utility Functions

Helper functions and utilities used across RAPTOR modules.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import os
import sys
import shutil
import logging
from pathlib import Path
from typing import Optional, List, Dict
import subprocess

logger = logging.getLogger(__name__)


# =============================================================================
# File and Directory Utilities
# =============================================================================

def ensure_dir(directory: str) -> Path:
    """
    Create directory if it doesn't exist.
    
    Parameters
    ----------
    directory : str
        Directory path
    
    Returns
    -------
    Path
        Path object for the directory
    
    Examples
    --------
    >>> output_dir = ensure_dir('results/analysis')
    >>> print(output_dir.exists())
    True
    """
    dir_path = Path(directory)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def check_file_exists(filepath: str, file_type: str = 'file') -> bool:
    """
    Check if file exists and is readable.
    
    Parameters
    ----------
    filepath : str
        Path to file
    file_type : str
        Type description for error message
    
    Returns
    -------
    bool
        True if file exists and is readable
    
    Raises
    ------
    FileNotFoundError
        If file doesn't exist
    PermissionError
        If file is not readable
    """
    path = Path(filepath)
    
    if not path.exists():
        raise FileNotFoundError(f"{file_type} not found: {filepath}")
    
    if not os.access(path, os.R_OK):
        raise PermissionError(f"{file_type} is not readable: {filepath}")
    
    return True


def get_file_size(filepath: str) -> float:
    """
    Get file size in MB.
    
    Parameters
    ----------
    filepath : str
        Path to file
    
    Returns
    -------
    float
        File size in MB
    """
    size_bytes = Path(filepath).stat().st_size
    size_mb = size_bytes / (1024 * 1024)
    return size_mb


def clean_directory(directory: str, pattern: str = '*'):
    """
    Remove files matching pattern from directory.
    
    Parameters
    ----------
    directory : str
        Directory path
    pattern : str
        Glob pattern for files to remove
    """
    dir_path = Path(directory)
    
    if not dir_path.exists():
        logger.warning(f"Directory does not exist: {directory}")
        return
    
    for item in dir_path.glob(pattern):
        if item.is_file():
            item.unlink()
            logger.debug(f"Removed file: {item}")
        elif item.is_dir():
            shutil.rmtree(item)
            logger.debug(f"Removed directory: {item}")


# =============================================================================
# System and Resource Utilities
# =============================================================================

def check_command_exists(command: str) -> bool:
    """
    Check if command exists in PATH.
    
    Parameters
    ----------
    command : str
        Command name
    
    Returns
    -------
    bool
        True if command exists
    
    Examples
    --------
    >>> has_star = check_command_exists('STAR')
    >>> print(has_star)
    True
    """
    return shutil.which(command) is not None


def check_required_tools(tools: List[str]) -> Dict[str, bool]:
    """
    Check if required tools are installed.
    
    Parameters
    ----------
    tools : list of str
        List of required tool names
    
    Returns
    -------
    dict
        Dictionary mapping tool names to availability
    
    Examples
    --------
    >>> tools = ['STAR', 'salmon', 'kallisto']
    >>> availability = check_required_tools(tools)
    >>> print(availability)
    {'STAR': True, 'salmon': True, 'kallisto': False}
    """
    availability = {}
    
    for tool in tools:
        availability[tool] = check_command_exists(tool)
        if not availability[tool]:
            logger.warning(f"Tool not found: {tool}")
    
    return availability


def run_command(cmd: List[str], cwd: Optional[str] = None, 
                capture_output: bool = True) -> subprocess.CompletedProcess:
    """
    Run shell command with error handling.
    
    Parameters
    ----------
    cmd : list of str
        Command and arguments
    cwd : str, optional
        Working directory
    capture_output : bool
        Whether to capture stdout/stderr
    
    Returns
    -------
    subprocess.CompletedProcess
        Completed process object
    
    Raises
    ------
    subprocess.CalledProcessError
        If command fails
    
    Examples
    --------
    >>> result = run_command(['echo', 'Hello World'])
    >>> print(result.stdout.decode().strip())
    Hello World
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=capture_output,
            check=True,
            text=True
        )
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Error: {e.stderr}")
        raise


def get_available_memory() -> float:
    """
    Get available system memory in GB.
    
    Returns
    -------
    float
        Available memory in GB
    
    Examples
    --------
    >>> mem = get_available_memory()
    >>> print(f"{mem:.1f} GB available")
    15.8 GB available
    """
    try:
        import psutil
        mem = psutil.virtual_memory()
        return mem.available / (1024 ** 3)
    except ImportError:
        logger.warning("psutil not available, cannot determine memory")
        return 0.0


def get_cpu_count() -> int:
    """
    Get number of available CPU cores.
    
    Returns
    -------
    int
        Number of CPU cores
    
    Examples
    --------
    >>> cores = get_cpu_count()
    >>> print(f"{cores} cores available")
    8 cores available
    """
    return os.cpu_count() or 1


# =============================================================================
# Data Validation Utilities
# =============================================================================

def validate_count_matrix(counts, min_genes: int = 100, min_samples: int = 2):
    """
    Validate count matrix format and content.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix
    min_genes : int
        Minimum number of genes
    min_samples : int
        Minimum number of samples
    
    Raises
    ------
    ValueError
        If validation fails
    """
    import pandas as pd
    import numpy as np
    
    if not isinstance(counts, pd.DataFrame):
        raise ValueError("Counts must be a pandas DataFrame")
    
    n_genes, n_samples = counts.shape
    
    if n_genes < min_genes:
        raise ValueError(f"Too few genes: {n_genes} < {min_genes}")
    
    if n_samples < min_samples:
        raise ValueError(f"Too few samples: {n_samples} < {min_samples}")
    
    # Check for negative values
    if (counts < 0).any().any():
        raise ValueError("Count matrix contains negative values")
    
    # Check for non-numeric values
    if not np.issubdtype(counts.values.dtype, np.number):
        raise ValueError("Count matrix must contain only numeric values")
    
    # Warn about all-zero genes
    zero_genes = (counts.sum(axis=1) == 0).sum()
    if zero_genes > 0:
        logger.warning(f"{zero_genes} genes have zero counts across all samples")
    
    logger.info(f"Count matrix validated: {n_genes} genes × {n_samples} samples")


def validate_metadata(metadata, counts, required_columns: Optional[List[str]] = None):
    """
    Validate sample metadata.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Sample metadata
    counts : pd.DataFrame
        Count matrix (for sample name checking)
    required_columns : list of str, optional
        Required column names
    
    Raises
    ------
    ValueError
        If validation fails
    """
    import pandas as pd
    
    if not isinstance(metadata, pd.DataFrame):
        raise ValueError("Metadata must be a pandas DataFrame")
    
    # Check sample correspondence
    if 'sample' in metadata.columns:
        meta_samples = set(metadata['sample'])
        count_samples = set(counts.columns)
        
        if meta_samples != count_samples:
            missing = count_samples - meta_samples
            extra = meta_samples - count_samples
            
            if missing:
                logger.warning(f"Samples in counts but not metadata: {missing}")
            if extra:
                logger.warning(f"Samples in metadata but not counts: {extra}")
    
    # Check required columns
    if required_columns:
        missing_cols = set(required_columns) - set(metadata.columns)
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
    
    logger.info(f"Metadata validated: {len(metadata)} samples")


# =============================================================================
# Progress and Logging Utilities
# =============================================================================

def setup_logging(level: str = 'INFO', log_file: Optional[str] = None):
    """
    Setup logging configuration.
    
    Parameters
    ----------
    level : str
        Logging level: 'DEBUG', 'INFO', 'WARNING', 'ERROR'
    log_file : str, optional
        Path to log file
    """
    log_level = getattr(logging, level.upper())
    
    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Add file handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        file_handler.setFormatter(formatter)
        logging.getLogger('').addHandler(file_handler)
        
        logger.info(f"Logging to file: {log_file}")


def print_banner(text: str, char: str = '='):
    """
    Print text banner.
    
    Parameters
    ----------
    text : str
        Banner text
    char : str
        Character for border
    
    Examples
    --------
    >>> print_banner("RAPTOR Analysis")
    ================
    RAPTOR Analysis
    ================
    """
    border = char * len(text)
    print(f"\n{border}\n{text}\n{border}\n")


def format_time(seconds: float) -> str:
    """
    Format time duration in human-readable format.
    
    Parameters
    ----------
    seconds : float
        Time in seconds
    
    Returns
    -------
    str
        Formatted time string
    
    Examples
    --------
    >>> print(format_time(3665))
    1h 1m 5s
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        mins = int(seconds / 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds / 3600)
        mins = int((seconds % 3600) / 60)
        secs = int(seconds % 60)
        return f"{hours}h {mins}m {secs}s"


def format_bytes(bytes: float) -> str:
    """
    Format bytes in human-readable format.
    
    Parameters
    ----------
    bytes : float
        Size in bytes
    
    Returns
    -------
    str
        Formatted size string
    
    Examples
    --------
    >>> print(format_bytes(1536000))
    1.5 MB
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes < 1024.0:
            return f"{bytes:.1f} {unit}"
        bytes /= 1024.0
    return f"{bytes:.1f} PB"


# =============================================================================
# Configuration Utilities
# =============================================================================

def load_config(config_file: str) -> Dict:
    """
    Load configuration from YAML file.
    
    Parameters
    ----------
    config_file : str
        Path to config file
    
    Returns
    -------
    dict
        Configuration dictionary
    """
    import yaml
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    logger.info(f"Loaded config: {config_file}")
    return config


def save_config(config: Dict, config_file: str):
    """
    Save configuration to YAML file.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary
    config_file : str
        Path to output file
    """
    import yaml
    
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    logger.info(f"Saved config: {config_file}")


# =============================================================================
# Version and Environment Info
# =============================================================================

def get_environment_info() -> Dict:
    """
    Get information about Python environment and dependencies.
    
    Returns
    -------
    dict
        Environment information
    """
    import platform
    
    info = {
        'python_version': sys.version,
        'platform': platform.platform(),
        'processor': platform.processor(),
    }
    
    # Try to get package versions
    try:
        import numpy
        info['numpy_version'] = numpy.__version__
    except ImportError:
        info['numpy_version'] = 'not installed'
    
    try:
        import pandas
        info['pandas_version'] = pandas.__version__
    except ImportError:
        info['pandas_version'] = 'not installed'
    
    return info


if __name__ == '__main__':
    print("RAPTOR Utilities")
    print("===============")
    print("\nUtility functions for RAPTOR package.")
    print("\nAvailable utilities:")
    print("  • File and directory management")
    print("  • System resource checking")
    print("  • Data validation")
    print("  • Logging and progress tracking")
    print("  • Configuration handling")
