# Contributing to RAPTOR

Thank you for considering contributing to RAPTOR! 🦖

---

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [How Can I Contribute?](#how-can-i-contribute)
3. [Development Setup](#development-setup)
4. [Coding Standards](#coding-standards)
5. [Submitting Changes](#submitting-changes)
6. [Testing](#testing)

---

## Code of Conduct

This project follows a code of conduct. By participating, you are expected to uphold this code. Please be respectful and constructive.

---

## How Can I Contribute?

### Reporting Bugs

**Before submitting a bug report:**
- Check existing issues to avoid duplicates
- Collect information: Python version, RAPTOR version, error messages

**Submit a bug report:**
- Use the bug report template
- Provide a clear title and description
- Include steps to reproduce
- Share error messages and logs
- Mention your environment (OS, Python version)

### Suggesting Features

**Feature requests are welcome!**
- Use the feature request template
- Explain the use case
- Describe the desired behavior
- Consider implementation complexity

### Code Contributions

**We welcome:**
- Bug fixes
- New pipeline implementations
- Performance improvements
- Documentation improvements
- Test coverage improvements

---

## Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- R 4.0+ (for Module 6)
- Virtual environment tool (venv or conda)

### Setup Steps

```bash
# 1. Fork and clone
git clone https://github.com/YOUR-USERNAME/RAPTOR.git
cd RAPTOR

# 2. Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
venv\Scripts\activate     # Windows

# 3. Install in editable mode with dev dependencies
pip install -e .[dev]

# 4. Install R packages (for Module 6)
Rscript raptor/external_modules/module6_de_analysis/r_scripts/install_packages.R

# 5. Run tests to verify
pytest tests/
```

---

## Coding Standards

### Python Style

**Follow PEP 8** with these specifics:

```python
# Line length: 100 characters (not 79)
# Imports: sorted with isort
# Formatting: use black

# Format code before committing
black raptor/
isort raptor/

# Check with flake8
flake8 raptor/
```

### Documentation

**All public functions must have docstrings:**

```python
def example_function(param1: str, param2: int) -> dict:
    """
    Brief description of function.
    
    Longer description if needed. Explain what the function does,
    any important details, edge cases, etc.
    
    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int
        Description of param2
    
    Returns
    -------
    dict
        Description of return value
    
    Examples
    --------
    >>> result = example_function('test', 42)
    >>> print(result)
    {'status': 'success'}
    
    Raises
    ------
    ValueError
        If param2 is negative
    """
    if param2 < 0:
        raise ValueError("param2 must be non-negative")
    
    return {'status': 'success', 'param1': param1, 'param2': param2}
```

### Type Hints

**Use type hints for all function signatures:**

```python
from typing import List, Dict, Optional, Union
from pathlib import Path

def process_files(
    file_paths: List[Path],
    config: Optional[Dict[str, str]] = None
) -> Dict[str, int]:
    """Process multiple files."""
    # Implementation
    pass
```

### Testing

**All new code must have tests:**

```python
# In tests/test_new_feature.py
import pytest
from raptor.new_module import new_function

def test_new_function_basic():
    """Test basic functionality."""
    result = new_function(input_data)
    assert result == expected_output

def test_new_function_edge_case():
    """Test edge case handling."""
    with pytest.raises(ValueError):
        new_function(invalid_input)
```

---

## Submitting Changes

### Git Workflow

```bash
# 1. Create a new branch
git checkout -b feature/my-new-feature
# or
git checkout -b fix/bug-description

# 2. Make changes and commit
git add .
git commit -m "Add feature: description"

# Use conventional commits:
# feat: New feature
# fix: Bug fix
# docs: Documentation changes
# test: Test additions
# refactor: Code refactoring
# style: Formatting changes

# 3. Push to your fork
git push origin feature/my-new-feature

# 4. Create Pull Request on GitHub
```

### Pull Request Process

**Before submitting:**
1. ✅ Run tests: `pytest tests/`
2. ✅ Format code: `black raptor/ && isort raptor/`
3. ✅ Check linting: `flake8 raptor/`
4. ✅ Update documentation if needed
5. ✅ Add tests for new features
6. ✅ Update CHANGELOG.md

**PR Description should include:**
- What changes were made
- Why the changes were needed
- How to test the changes
- Related issue numbers (if any)

**Example PR template:**
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement

## Testing
How to test these changes

## Checklist
- [ ] Tests pass
- [ ] Code formatted with black
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
```

---

## Testing

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=raptor --cov-report=html

# Run specific test file
pytest tests/test_ensemble.py

# Run specific test
pytest tests/test_ensemble.py::test_fisher_method

# Run with verbose output
pytest tests/ -v
```

### Writing Tests

**Test file structure:**

```
tests/
├── conftest.py           # Shared fixtures
├── test_module_name.py   # Tests for module_name
└── test_data/           # Test data files
```

**Example test:**

```python
import pytest
from raptor.ensemble import ensemble_fisher

@pytest.fixture
def sample_de_results():
    """Create sample DE results for testing."""
    # Setup test data
    return de_results

def test_ensemble_fisher_basic(sample_de_results):
    """Test basic Fisher's method functionality."""
    result = ensemble_fisher(sample_de_results)
    
    assert result is not None
    assert len(result.consensus_genes) > 0
    assert result.ensemble_method == 'fisher'

def test_ensemble_fisher_threshold():
    """Test Fisher's method with custom threshold."""
    result = ensemble_fisher(
        de_results,
        significance_threshold=0.01
    )
    assert all(result.consensus_genes['combined_pvalue'] < 0.01)
```

---

## Adding New Modules

### Module Structure

```
raptor/
└── new_module/
    ├── __init__.py       # Exports
    ├── core.py           # Main functionality
    ├── utils.py          # Helper functions
    └── tests/
        └── test_new_module.py
```

### Integration Checklist

When adding a new module:

- [ ] Add to `raptor/__init__.py` imports
- [ ] Update `setup.py` if new dependencies needed
- [ ] Create documentation in `docs/MODULE_X_Name.md`
- [ ] Add example script in `examples/0X_module_name.py`
- [ ] Write comprehensive tests
- [ ] Update README.md
- [ ] Update CHANGELOG.md

---

## Adding New Pipelines

### Pipeline Structure

```
raptor/pipelines/new_pipeline/
├── __init__.py
├── pipeline.py        # Main pipeline class
├── config.yaml        # Default configuration
└── run_pipeline.sh    # Shell script
```

### Pipeline Checklist

- [ ] Inherit from `BasePipeline`
- [ ] Implement all required methods
- [ ] Add to `raptor/pipelines/__init__.py` registry
- [ ] Update `PIPELINE_METADATA` dictionary
- [ ] Add config.yaml with parameters
- [ ] Write validation checks (10+ recommended)
- [ ] Add to CLI in `raptor/cli.py`
- [ ] Create example usage script
- [ ] Write tests (minimum 5 tests)
- [ ] Update pipeline documentation

---

## Documentation Guidelines

### Module Documentation

Each module should have:

1. **README section** describing purpose
2. **API documentation** for all public functions
3. **Usage examples** showing common workflows
4. **Parameter descriptions** explaining all options

### Code Comments

```python
# Good comment: Explains WHY
# Use Fisher's method here because it handles correlation better
result = fishers_method(pvalues)

# Bad comment: States the obvious
# Calculate mean
mean = sum(values) / len(values)
```

---

## Release Process

(For maintainers only)

### Version Bump

```bash
# 1. Update version in setup.py
VERSION = "2.3.0"

# 2. Update CHANGELOG.md
## [2.3.0] - 2026-XX-XX
### Added
- New feature description

# 3. Commit version bump
git commit -am "Bump version to 2.3.0"

# 4. Tag release
git tag -a v2.3.0 -m "Release version 2.3.0"
git push origin v2.3.0

# 5. Build and publish to PyPI
python -m build
twine upload dist/*
```

---

## Questions?

- 📧 Email: ayehbolouki1988@gmail.com
- 🐛 Issues: https://github.com/AyehBlk/RAPTOR/issues
- 💬 Discussions: https://github.com/AyehBlk/RAPTOR/discussions

---

**Thank you for contributing to RAPTOR!** 🙏
