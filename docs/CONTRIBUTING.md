# 🤝 Contributing to RAPTOR

Thank you for your interest in contributing to RAPTOR! This guide will help you get started.

---

## 📋 Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Setup](#development-setup)
4. [How to Contribute](#how-to-contribute)
5. [Coding Standards](#coding-standards)
6. [Testing](#testing)
7. [Documentation](#documentation)
8. [Pull Request Process](#pull-request-process)
9. [Community](#community)

---

## 🤗 Code of Conduct

### Our Pledge

We are committed to making participation in RAPTOR a harassment-free experience for everyone, regardless of:
- Age, body size, disability, ethnicity
- Gender identity and expression
- Level of experience
- Nationality, personal appearance
- Race, religion
- Sexual identity and orientation

### Our Standards

**Positive behaviors:**
- ✅ Using welcoming and inclusive language
- ✅ Being respectful of differing viewpoints
- ✅ Gracefully accepting constructive criticism
- ✅ Focusing on what's best for the community
- ✅ Showing empathy toward others

**Unacceptable behaviors:**
- ❌ Trolling, insulting/derogatory comments
- ❌ Public or private harassment
- ❌ Publishing others' private information
- ❌ Other unethical or unprofessional conduct

### Enforcement

Report unacceptable behavior to: ayehbolouki1988@gmail.com

---

## 🚀 Getting Started

### Prerequisites

- Python 3.8 or higher
- Git
- GitHub account
- Basic knowledge of RNA-seq analysis (helpful but not required)

### First Steps

1. **Star the repository** ⭐
2. **Fork the repository** 🍴
3. **Clone your fork** 📂
4. **Set up development environment** 💻
5. **Find an issue to work on** 🎯

---

## 💻 Development Setup

### 1. Fork and Clone

```bash
# Fork on GitHub first, then:
git clone https://github.com/YOUR_USERNAME/RAPTOR.git
cd RAPTOR
```

### 2. Create Virtual Environment

```bash
# Using venv
python -m venv raptor_dev
source raptor_dev/bin/activate  # On Windows: raptor_dev\Scripts\activate

# OR using conda
conda create -n raptor_dev python=3.9
conda activate raptor_dev
```

### 3. Install Development Dependencies

```bash
# Install in development mode
pip install -e .

# Install development dependencies
pip install -r requirements-dev.txt
```

**requirements-dev.txt includes:**
- pytest (testing)
- pytest-cov (coverage)
- black (formatting)
- flake8 (linting)
- mypy (type checking)
- sphinx (documentation)
- pre-commit (git hooks)

### 4. Set Up Pre-commit Hooks

```bash
pre-commit install
```

This automatically checks:
- Code formatting
- Linting
- Type hints
- Trailing whitespace
- YAML validation

### 5. Configure Git

```bash
git config user.name "Your Name"
git config user.email "your.email@example.com"

# Set upstream remote
git remote add upstream https://github.com/AyehBlk/RAPTOR.git
```

### 6. Verify Installation

```bash
# Run tests
pytest

# Check formatting
black --check raptor/

# Run linting
flake8 raptor/

# Type checking
mypy raptor/
```

---

## 🎯 How to Contribute

### Types of Contributions

#### 1. 🐛 Bug Reports

**Before submitting:**
- Check [existing issues](https://github.com/AyehBlk/RAPTOR/issues)
- Try latest version
- Gather necessary info

**Include in report:**
```markdown
**Bug Description**
Clear description of the bug

**To Reproduce**
1. Step 1
2. Step 2
3. See error

**Expected Behavior**
What should happen

**System Info**
- OS: [e.g., Ubuntu 22.04]
- Python: [e.g., 3.9.7]
- RAPTOR: [e.g., 2.1.0]

**Logs**
```
Paste relevant logs here
```

**Config**
```yaml
# Paste config.yaml
```
```

#### 2. 💡 Feature Requests

**Before suggesting:**
- Check [existing requests](https://github.com/AyehBlk/RAPTOR/discussions)
- Consider if it fits RAPTOR's scope
- Think about implementation

**Include in request:**
- Clear use case
- Expected behavior
- Alternative approaches considered
- Willingness to implement

#### 3. 📝 Documentation

**What to improve:**
- Fix typos
- Clarify confusing sections
- Add examples
- Translate to other languages
- Improve API docs

**How to contribute:**
```bash
# Edit .md files directly
git checkout -b docs/improve-installation
# Make changes
git commit -m "docs: clarify installation steps"
git push origin docs/improve-installation
# Create pull request
```

#### 4. 🔧 Code Contributions

**Good first issues:**
- Look for [`good first issue`](https://github.com/AyehBlk/RAPTOR/labels/good%20first%20issue) label
- Check [`help wanted`](https://github.com/AyehBlk/RAPTOR/labels/help%20wanted) label

**Types of contributions:**
- Bug fixes
- New features
- Performance improvements
- Refactoring
- Test coverage

---

## 📋 Coding Standards

### Python Style Guide

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) with some modifications:

#### Formatting

**Use Black formatter:**
```bash
black raptor/
```

**Line length:** 100 characters (not 79)

**Imports:**
```python
# Standard library
import os
import sys

# Third-party
import numpy as np
import pandas as pd

# Local
from raptor.core import RAPTORAnalysis
from raptor.utils import validate_config
```

#### Naming Conventions

```python
# Classes: PascalCase
class PipelineAnalyzer:
    pass

# Functions: snake_case
def run_analysis():
    pass

# Variables: snake_case
sample_count = 10
fastq_files = []

# Constants: UPPER_CASE
MAX_THREADS = 32
DEFAULT_CONFIG = "config.yaml"

# Private: leading underscore
def _internal_function():
    pass

_internal_var = 42
```

#### Type Hints

**Use type hints for all functions:**

```python
from typing import List, Dict, Optional, Union
from pathlib import Path

def analyze_samples(
    fastq_dir: Path,
    output_dir: Path,
    threads: int = 8,
    pipelines: Optional[List[str]] = None
) -> Dict[str, pd.DataFrame]:
    """
    Analyze RNA-seq samples.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        output_dir: Output directory for results
        threads: Number of threads to use
        pipelines: List of pipelines to run (None = all)
        
    Returns:
        Dictionary mapping pipeline names to result DataFrames
        
    Raises:
        FileNotFoundError: If fastq_dir doesn't exist
        ValueError: If invalid pipeline specified
    """
    pass
```

#### Docstrings

**Use Google-style docstrings:**

```python
def calculate_tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """Calculate TPM (Transcripts Per Million) values.
    
    This function normalizes raw count data to TPM, accounting for
    both sequencing depth and gene length.
    
    Args:
        counts: Raw count matrix (genes × samples)
        lengths: Gene lengths in base pairs
        
    Returns:
        TPM-normalized expression matrix
        
    Example:
        >>> counts = pd.DataFrame({'sample1': [100, 200]})
        >>> lengths = pd.Series([1000, 2000])
        >>> tpm = calculate_tpm(counts, lengths)
        
    Note:
        Genes with zero length are excluded from calculations.
    """
    # Implementation
    pass
```

### Code Quality

#### Write Clean Code

```python
# ❌ BAD
def f(x,y,z):
    return x+y+z if z>0 else x-y

# ✅ GOOD
def calculate_sum(
    value_a: float,
    value_b: float,
    condition: float
) -> float:
    """Calculate sum or difference based on condition."""
    if condition > 0:
        return value_a + value_b + condition
    else:
        return value_a - value_b
```

#### Use Meaningful Names

```python
# ❌ BAD
def proc(d):
    for i in d:
        x = i * 2
        y.append(x)
    return y

# ✅ GOOD
def normalize_counts(raw_counts: List[int]) -> List[float]:
    """Normalize count data by doubling each value."""
    normalized = []
    for count in raw_counts:
        normalized_count = count * 2
        normalized.append(normalized_count)
    return normalized
```

#### Avoid Magic Numbers

```python
# ❌ BAD
if quality_score < 30:
    filter_sample()

# ✅ GOOD
MINIMUM_QUALITY_SCORE = 30

if quality_score < MINIMUM_QUALITY_SCORE:
    filter_sample()
```

#### Handle Errors Gracefully

```python
# ❌ BAD
def read_file(path):
    return open(path).read()

# ✅ GOOD
def read_file(path: Path) -> str:
    """Read file contents safely."""
    try:
        with open(path, 'r') as f:
            return f.read()
    except FileNotFoundError:
        logger.error(f"File not found: {path}")
        raise
    except PermissionError:
        logger.error(f"Permission denied: {path}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error reading {path}: {e}")
        raise
```

---

## 🧪 Testing

### Writing Tests

**Use pytest framework:**

```python
# tests/test_analysis.py
import pytest
import pandas as pd
from raptor.analysis import RAPTORAnalysis

def test_basic_analysis():
    """Test basic analysis workflow."""
    analysis = RAPTORAnalysis()
    result = analysis.run(
        fastq_dir="tests/data/fastq",
        output_dir="tests/output"
    )
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0

def test_invalid_config():
    """Test that invalid config raises error."""
    with pytest.raises(ValueError):
        RAPTORAnalysis(config="invalid.yaml")

@pytest.fixture
def sample_data():
    """Fixture providing sample test data."""
    return pd.DataFrame({
        'gene': ['GENE1', 'GENE2'],
        'count': [100, 200]
    })

def test_normalization(sample_data):
    """Test TPM normalization."""
    from raptor.normalization import calculate_tpm
    lengths = pd.Series([1000, 2000], index=['GENE1', 'GENE2'])
    tpm = calculate_tpm(sample_data, lengths)
    assert tpm['count'].sum() == pytest.approx(1e6)
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_analysis.py

# Run specific test
pytest tests/test_analysis.py::test_basic_analysis

# Run with coverage
pytest --cov=raptor --cov-report=html

# Run fast tests only
pytest -m "not slow"

# Verbose output
pytest -v

# Stop on first failure
pytest -x
```

### Test Coverage

**Aim for ≥80% coverage:**

```bash
# Generate coverage report
pytest --cov=raptor --cov-report=term-missing

# View HTML report
open htmlcov/index.html
```

### Test Categories

```python
import pytest

@pytest.mark.unit
def test_unit():
    """Unit test - fast, isolated."""
    pass

@pytest.mark.integration
def test_integration():
    """Integration test - multiple components."""
    pass

@pytest.mark.slow
def test_slow():
    """Slow test - long running."""
    pass

@pytest.mark.cloud
def test_cloud():
    """Cloud test - requires cloud credentials."""
    pass
```

---

## 📚 Documentation

### Types of Documentation

#### 1. Code Documentation

**Inline comments for complex logic:**
```python
# Calculate normalized expression accounting for batch effects
# using ComBat-Seq algorithm (Zhang et al., 2020)
normalized = combat_seq(raw_counts, batch_info)
```

**Docstrings for all public functions:**
```python
def analyze_differential_expression(
    counts: pd.DataFrame,
    groups: List[str]
) -> pd.DataFrame:
    """Identify differentially expressed genes.
    
    Uses DESeq2-like normalization and statistical testing.
    
    Args:
        counts: Raw count matrix
        groups: Sample group labels
        
    Returns:
        DataFrame with statistics (log2FC, p-value, FDR)
    """
    pass
```

#### 2. User Documentation

**Markdown files in docs/:**
- User guides
- Tutorials
- API reference
- Configuration guides
- Troubleshooting

**Build documentation:**
```bash
cd docs/
make html
open _build/html/index.html
```

#### 3. Examples

**Add working examples:**
```python
# examples/basic_analysis.py
"""
Example: Basic RNA-seq Analysis with RAPTOR

This example demonstrates a simple RNA-seq analysis workflow.
"""

from raptor import RAPTORAnalysis

# Initialize analysis
raptor = RAPTORAnalysis(config="config.yaml")

# Run analysis
results = raptor.run(
    fastq_dir="data/fastq",
    output_dir="results"
)

# Generate report
raptor.generate_report(output="report.html")
```

---

## 🔄 Pull Request Process

### 1. Create a Branch

```bash
# Update your fork
git checkout main
git pull upstream main

# Create feature branch
git checkout -b feature/add-new-pipeline

# OR for bug fix
git checkout -b fix/memory-leak-issue-42
```

**Branch naming:**
- `feature/description` - New features
- `fix/description` - Bug fixes
- `docs/description` - Documentation
- `refactor/description` - Code refactoring
- `test/description` - Test improvements

### 2. Make Changes

```bash
# Make your changes
vim raptor/pipelines/new_pipeline.py

# Add tests
vim tests/test_new_pipeline.py

# Update documentation
vim docs/pipelines.md
```

### 3. Test Your Changes

```bash
# Run tests
pytest

# Check formatting
black --check raptor/

# Run linting
flake8 raptor/

# Type checking
mypy raptor/

# Build docs
cd docs && make html
```

### 4. Commit Changes

**Follow conventional commits:**

```bash
# Format: <type>(<scope>): <description>

git commit -m "feat(pipelines): add support for CellRanger"
git commit -m "fix(dashboard): resolve plot rendering issue #42"
git commit -m "docs(api): improve calculate_tpm documentation"
git commit -m "test(ml): add tests for model training"
git commit -m "refactor(core): simplify config validation"
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `test`: Tests
- `refactor`: Code refactoring
- `perf`: Performance improvement
- `style`: Formatting
- `chore`: Maintenance

### 5. Push to GitHub

```bash
git push origin feature/add-new-pipeline
```

### 6. Create Pull Request

**On GitHub:**
1. Click "New Pull Request"
2. Select your branch
3. Fill out PR template:

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation
- [ ] Refactoring

## Testing
- [ ] Tests pass locally
- [ ] Added new tests
- [ ] Updated documentation

## Checklist
- [ ] Code follows style guide
- [ ] Self-reviewed code
- [ ] Commented complex sections
- [ ] Updated documentation
- [ ] No breaking changes (or documented)

## Related Issues
Closes #42
```

### 7. Code Review Process

**What to expect:**
1. Automated checks run (CI/CD)
2. Maintainer reviews code
3. Feedback and requests for changes
4. Approval and merge

**Respond to feedback:**
```bash
# Make requested changes
git add .
git commit -m "refactor: address review feedback"
git push origin feature/add-new-pipeline
```

### 8. After Merge

```bash
# Update your fork
git checkout main
git pull upstream main

# Delete feature branch
git branch -d feature/add-new-pipeline
git push origin --delete feature/add-new-pipeline
```

---

## 👥 Community

### Communication Channels

- **GitHub Issues**: Bug reports, feature requests
- **GitHub Discussions**: Questions, ideas, general discussion
- **Email**: ayehbolouki1988@gmail.com

### Getting Help

**Before asking:**
1. Check [FAQ](FAQ.md)
2. Search [existing issues](https://github.com/AyehBlk/RAPTOR/issues)
3. Read [documentation](docs/)

**When asking:**
- Be specific
- Provide context
- Include minimal reproducible example
- Share relevant logs/errors

### Recognition

**Contributors will be:**
- Listed in [CONTRIBUTORS.md](CONTRIBUTORS.md)
- Mentioned in release notes
- Added to Zenodo citation (for major contributions)

### Maintainers

Current maintainers:
- **Ayeh Bolouki** (@AyehBlk) - Lead Developer

---

## 📋 Checklist for Contributions

Before submitting, ensure:

- [ ] Code follows style guide
- [ ] Added/updated tests
- [ ] All tests pass
- [ ] Updated documentation
- [ ] Added docstrings
- [ ] Included type hints
- [ ] Committed with conventional commits
- [ ] Referenced related issues
- [ ] Requested review from maintainer

---

## 🎓 Learning Resources

### Python
- [Python Documentation](https://docs.python.org/3/)
- [Real Python](https://realpython.com/)
- [Python Type Hints](https://docs.python.org/3/library/typing.html)

### Testing
- [pytest Documentation](https://docs.pytest.org/)
- [Test-Driven Development](https://www.obeythetestinggoat.com/)

### Git
- [Pro Git Book](https://git-scm.com/book/en/v2)
- [GitHub Guides](https://guides.github.com/)

### RNA-seq
- [RNA-seq Tutorial](https://github.com/griffithlab/rnaseq_tutorial)
- [RNA-seqlopedia](https://rnaseq.uoregon.edu/)

---

## 📞 Contact

**Questions about contributing?**
- Email: ayehbolouki1988@gmail.com
- GitHub: @AyehBlk

---

## 🙏 Thank You!

Your contributions make RAPTOR better for everyone in the RNA-seq community!

Every contribution, no matter how small, is valued and appreciated. 💙

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**License:** MIT

---

*"Alone we can do so little; together we can do so much."* - Helen Keller 🤝
