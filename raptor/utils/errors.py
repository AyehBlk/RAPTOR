"""
RAPTOR v2.2.0 - Error Handling Utilities

Comprehensive error handling and user-friendly error messages.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import sys
import logging
import traceback
from pathlib import Path
from typing import Optional, Callable, Any
from functools import wraps

logger = logging.getLogger(__name__)


# =============================================================================
# Custom RAPTOR Exceptions
# =============================================================================

class RAPTORError(Exception):
    """Base exception for RAPTOR."""
    pass


class ValidationError(RAPTORError):
    """Raised when input validation fails."""
    pass


class PipelineError(RAPTORError):
    """Raised when pipeline execution fails."""
    pass


class FileFormatError(RAPTORError):
    """Raised when file format is invalid."""
    pass


class DependencyError(RAPTORError):
    """Raised when required dependency is missing."""
    pass


class ConfigurationError(RAPTORError):
    """Raised when configuration is invalid."""
    pass


# =============================================================================
# Module 8: Parameter Optimization Errors
# =============================================================================

class OptimizationError(RAPTORError):
    """Raised when parameter optimization fails."""
    pass


class GroundTruthError(RAPTORError):
    """
    Raised when ground truth validation fails.
    
    This is specifically for ground truth issues like:
    - Too few validated genes (< 10)
    - Missing required columns
    - Invalid gene identifiers
    """
    pass


class InsufficientDataError(RAPTORError):
    """
    Raised when insufficient data for optimization.
    
    Examples:
    - Too few genes for FDR estimation (< 1000)
    - Too few bootstrap iterations (< 50)
    - Empty DE results
    """
    pass


# =============================================================================
# Error Handler Decorator
# =============================================================================

def handle_errors(exit_on_error: bool = True, 
                 log_traceback: bool = True):
    """
    Decorator for comprehensive error handling.
    
    Parameters
    ----------
    exit_on_error : bool
        If True, exit program on error
    log_traceback : bool
        If True, log full traceback for unexpected errors
    
    Examples
    --------
    @handle_errors(exit_on_error=True)
    def my_function():
        # Function code
        pass
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            try:
                return func(*args, **kwargs)
            
            # Known RAPTOR errors - user-friendly messages
            except ValidationError as e:
                logger.error(f"❌ Validation Error: {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except GroundTruthError as e:
                logger.error(f"❌ Ground Truth Error: {e}")
                logger.info("    See Module 8 documentation for ground truth requirements")
                logger.info("    Minimum: 10 genes, Recommended: 20+, Ideal: 50+")
                logger.info("    Alternative: Use FDR control, stability, or reproducibility methods")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except InsufficientDataError as e:
                logger.error(f"❌ Insufficient Data: {e}")
                logger.info("    Check data requirements for the optimization method")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except OptimizationError as e:
                logger.error(f"❌ Optimization Error: {e}")
                logger.info("    Try a different optimization strategy or method")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except FileNotFoundError as e:
                logger.error(f"❌ File Not Found: {e}")
                logger.info("    Check file path and permissions")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except FileFormatError as e:
                logger.error(f"❌ File Format Error: {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except DependencyError as e:
                logger.error(f"❌ Missing Dependency: {e}")
                logger.info("    Install required packages with:")
                logger.info("    pip install -r requirements.txt")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except ConfigurationError as e:
                logger.error(f"❌ Configuration Error: {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except PipelineError as e:
                logger.error(f"❌ Pipeline Error: {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except ValueError as e:
                logger.error(f"❌ Invalid Value: {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except TypeError as e:
                logger.error(f"❌ Type Error: {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except KeyError as e:
                logger.error(f"❌ Key Error: Missing key {e}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except PermissionError as e:
                logger.error(f"❌ Permission Denied: {e}")
                logger.info("    Check file/directory permissions")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            except ImportError as e:
                logger.error(f"❌ Import Error: {e}")
                logger.info("    Install required package:")
                logger.info(f"    pip install {str(e).split()[-1]}")
                if exit_on_error:
                    sys.exit(1)
                raise
            
            # Unexpected errors - full traceback
            except Exception as e:
                logger.error(f"❌ Unexpected Error: {type(e).__name__}: {e}")
                if log_traceback:
                    logger.error("\n" + "="*60)
                    logger.error("FULL TRACEBACK:")
                    logger.error("="*60)
                    logger.error(traceback.format_exc())
                    logger.error("="*60)
                
                logger.info("\n💡 This is an unexpected error. Please:")
                logger.info("   1. Check your inputs are correct")
                logger.info("   2. Report this issue: https://github.com/AyehBlk/RAPTOR/issues")
                logger.info("   3. Include the error message and traceback above")
                
                if exit_on_error:
                    sys.exit(1)
                raise
        
        return wrapper
    return decorator


# =============================================================================
# Specific Error Helpers
# =============================================================================

def check_file_exists(filepath: Path, file_description: str = "File") -> Path:
    """
    Check if file exists with informative error.
    
    Parameters
    ----------
    filepath : Path
        File path to check
    file_description : str
        Description for error message
    
    Returns
    -------
    Path
        File path if exists
    
    Raises
    ------
    FileNotFoundError
        If file doesn't exist
    """
    if not filepath.exists():
        raise FileNotFoundError(
            f"{file_description} not found: {filepath}\n"
            f"    Current directory: {Path.cwd()}\n"
            f"    Please check the path and try again."
        )
    return filepath


def check_dependency(module_name: str, 
                     package_name: Optional[str] = None,
                     install_name: Optional[str] = None) -> None:
    """
    Check if Python module is available.
    
    Parameters
    ----------
    module_name : str
        Module name to import
    package_name : str, optional
        Package name (if different from module)
    install_name : str, optional
        Name to use for pip install
    
    Raises
    ------
    DependencyError
        If module not available
    """
    try:
        __import__(module_name)
    except ImportError:
        install_pkg = install_name or package_name or module_name
        raise DependencyError(
            f"Required module '{module_name}' not available.\n"
            f"    Install with: pip install {install_pkg}"
        )


def check_tool_available(tool_name: str, 
                        install_instructions: Optional[str] = None) -> None:
    """
    Check if external tool is available in PATH.
    
    Parameters
    ----------
    tool_name : str
        Tool name (e.g., 'salmon', 'kallisto')
    install_instructions : str, optional
        Installation instructions
    
    Raises
    ------
    DependencyError
        If tool not available
    """
    import shutil
    
    if not shutil.which(tool_name):
        msg = f"Required tool '{tool_name}' not found in PATH."
        if install_instructions:
            msg += f"\n    {install_instructions}"
        else:
            msg += f"\n    Please install {tool_name} and ensure it's in your PATH."
        raise DependencyError(msg)


def validate_output_writable(output_path: Path) -> Path:
    """
    Check if output path is writable.
    
    Parameters
    ----------
    output_path : Path
        Output file or directory path
    
    Returns
    -------
    Path
        Validated path
    
    Raises
    ------
    PermissionError
        If not writable
    """
    # Check parent directory is writable
    parent = output_path.parent if output_path.suffix else output_path
    
    if parent.exists():
        if not parent.is_dir():
            raise ValidationError(f"Parent path exists but is not a directory: {parent}")
        
        # Try to create a test file
        test_file = parent / ".raptor_write_test"
        try:
            test_file.touch()
            test_file.unlink()
        except PermissionError:
            raise PermissionError(
                f"Cannot write to directory: {parent}\n"
                f"    Check directory permissions."
            )
    else:
        # Try to create parent directory
        try:
            parent.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            raise PermissionError(
                f"Cannot create directory: {parent}\n"
                f"    Check parent directory permissions."
            )
    
    return output_path


# =============================================================================
# User-Friendly Error Messages
# =============================================================================

def format_user_error(error_type: str, 
                     message: str, 
                     suggestion: Optional[str] = None,
                     details: Optional[str] = None) -> str:
    """
    Format user-friendly error message.
    
    Parameters
    ----------
    error_type : str
        Error type (e.g., "File Not Found")
    message : str
        Error message
    suggestion : str, optional
        Suggestion for fixing
    details : str, optional
        Additional details
    
    Returns
    -------
    str
        Formatted error message
    """
    lines = [
        "",
        "=" * 60,
        f"❌ ERROR: {error_type}",
        "=" * 60,
        "",
        message,
    ]
    
    if details:
        lines.extend([
            "",
            "Details:",
            f"  {details}",
        ])
    
    if suggestion:
        lines.extend([
            "",
            "💡 Suggestion:",
            f"  {suggestion}",
        ])
    
    lines.extend([
        "",
        "=" * 60,
        ""
    ])
    
    return "\n".join(lines)


# =============================================================================
# Progress and Status Helpers
# =============================================================================

class ProgressIndicator:
    """Simple progress indicator for long operations."""
    
    def __init__(self, total: int, description: str = "Processing"):
        """
        Initialize progress indicator.
        
        Parameters
        ----------
        total : int
            Total number of items
        description : str
            Description of task
        """
        self.total = total
        self.description = description
        self.current = 0
        self.last_percent = -1
    
    def update(self, n: int = 1):
        """Update progress by n items."""
        self.current += n
        percent = int(100 * self.current / self.total)
        
        if percent != self.last_percent and percent % 10 == 0:
            logger.info(f"  {self.description}: {percent}% ({self.current}/{self.total})")
            self.last_percent = percent
    
    def finish(self):
        """Mark as finished."""
        logger.info(f"  ✓ {self.description}: Complete ({self.total}/{self.total})")


# =============================================================================
# Testing Functions
# =============================================================================

def test_error_handling():
    """Test error handling functions."""
    
    print("Testing RAPTOR Error Handling")
    print("=" * 60)
    
    # Test custom exceptions
    print("\n1. Testing custom exceptions:")
    try:
        raise ValidationError("Test validation error")
    except ValidationError as e:
        print(f"   ✓ ValidationError: {e}")
    
    # Test file checking
    print("\n2. Testing file checking:")
    try:
        check_file_exists(Path("/nonexistent/file.txt"), "Test file")
    except FileNotFoundError as e:
        print(f"   ✓ Caught: {str(e).split(chr(10))[0]}")
    
    # Test dependency checking
    print("\n3. Testing dependency checking:")
    try:
        check_dependency("nonexistent_module", install_name="fake-package")
    except DependencyError as e:
        print(f"   ✓ Caught: {str(e).split(chr(10))[0]}")
    
    print("\n" + "=" * 60)
    print("✓ All tests passed!")


if __name__ == '__main__':
    test_error_handling()
# =============================================================================
# Module 9: Ensemble Analysis Errors
# =============================================================================

class EnsembleError(RAPTORError):
    """
    Base exception for ensemble analysis errors.
    
    Raised when ensemble analysis fails for any reason.
    """
    pass


class MethodMismatchError(EnsembleError):
    """
    Raised when DE methods don't have compatible results.
    
    Examples:
    - Different gene sets with no overlap
    - Incompatible column names
    - Missing required columns in some methods
    """
    pass


class DirectionInconsistencyError(EnsembleError):
    """
    Raised when genes have inconsistent direction across methods.
    
    This is a warning error - genes with direction inconsistency
    are flagged but may still be included depending on settings.
    """
    pass


class CombinationFailedError(EnsembleError):
    """
    Raised when p-value or rank combination fails.
    
    Examples:
    - Invalid p-values (p=0, p=1, p=NaN)
    - Too few genes for RRA
    - Correlation matrix computation fails
    """
    pass
