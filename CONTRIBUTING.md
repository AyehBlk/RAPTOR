# Contributing to RAPTOR 🦖

Thank you for your interest in contributing to RAPTOR! This project thrives on community input and we welcome contributions from researchers, bioinformaticians, and developers worldwide.

##  Ways to Contribute

### 1. Report Bugs 🐛

Found a bug? Please help us fix it!

**Before reporting:**
- Check if the issue already exists in [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
- Make sure you're using the latest version

**What to include:**
- Clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- Your environment (OS, Python version, RAPTOR version)
- Error messages or logs
- Example data if possible (or synthetic example)

**Template:**
```markdown
**Bug Description:**
Brief description

**Steps to Reproduce:**
1. Step one
2. Step two
3. Step three

**Expected Behavior:**
What should happen

**Actual Behavior:**
What actually happens

**Environment:**
- OS: Ubuntu 22.04
- Python: 3.10
- RAPTOR: 2.0.0

**Error Message:**
```
Paste error here
```
```

### 2. Suggest Features 

Have an idea to improve RAPTOR?

**Good feature requests include:**
- Clear description of the feature
- Why it would be useful
- How it relates to RNA-seq analysis
- Example use case
- Any references or papers supporting the idea

Open a GitHub Discussion or Issue with label `enhancement`.

### 3. Improve Documentation 

Documentation improvements are always welcome!

**Areas that need help:**
- Clarifying existing documentation
- Adding examples
- Fixing typos
- Adding tutorials
- Translating documentation
- Creating video walkthroughs

**How to contribute:**
1. Fork the repository
2. Edit files in `docs/` folder
3. Submit a Pull Request

### 4. Add New Pipelines 

Want to add a new RNA-seq pipeline to RAPTOR?

**Requirements:**
- Complete workflow (alignment/quantification + statistics)
- Widely used or novel method with publication
- Reproducible implementation
- Tests demonstrating it works

**Process:**
1. Open an Issue to discuss the pipeline
2. Create a new folder in `pipelines/`
3. Follow the pipeline template structure
4. Add documentation
5. Include test data
6. Submit Pull Request

See [Pipeline Development Guide](docs/PIPELINE_DEVELOPMENT.md) for details.

### 5. Share Benchmark Results 

Ran RAPTOR on your data? Share your results!

**What to share:**
- Dataset characteristics (size, organism, design)
- Pipelines compared
- Performance results
- Any insights or surprises
- Publication reference if applicable

This helps improve recommendations for the community!

### 6. Fix Issues 

Want to fix an existing issue?

**Good first issues:**
- Look for `good first issue` label
- Issues labeled `help wanted`
- Documentation improvements
- Test coverage

**Before starting:**
- Comment on the issue to claim it
- Ask questions if anything is unclear
- Discuss approach if it's a big change

---

##  Development Process

### Setting Up Development Environment

```bash
# 1. Fork and clone
git clone https://github.com/YOUR-USERNAME/RAPTOR.git
cd RAPTOR

# 2. Create development environment
conda env create -f environment_dev.yml
conda activate raptor-dev

# 3. Install in development mode
pip install -e .

# 4. Run tests to verify setup
pytest tests/
```

### Making Changes

1. **Create a branch:**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/bug-description
   ```

2. **Make your changes:**
   - Write clean, readable code
   - Follow existing code style
   - Add comments where helpful
   - Update documentation if needed

3. **Test your changes:**
   ```bash
   # Run all tests
   pytest tests/
   
   # Run specific test
   pytest tests/test_profiler.py
   
   # Check code style
   flake8 raptor/
   
   # Check type hints
   mypy raptor/
   ```

4. **Commit your changes:**
   ```bash
   git add .
   git commit -m "Add feature: clear description"
   ```
   
   **Good commit messages:**
   - Clear and descriptive
   - Present tense ("Add feature" not "Added feature")
   - Reference issue numbers when applicable
   
   Examples:
   - ✅ `Add zero-inflation detection to profiler (#42)`
   - ✅ `Fix memory leak in benchmark module`
   - ✅ `Update documentation for profile command`
   - ❌ `fixed stuff`
   - ❌ `updates`

5. **Push to your fork:**
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Submit Pull Request:**
   - Go to GitHub and create Pull Request
   - Fill in the PR template
   - Link related issues
   - Describe what changed and why

---

##  Pull Request Guidelines

### Before Submitting

- ✅ Code follows project style
- ✅ Tests pass (`pytest tests/`)
- ✅ Documentation updated if needed
- ✅ No unnecessary files included
- ✅ Commits are clean and logical
- ✅ Branch is up to date with main

### PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Related Issues
Fixes #(issue number)

## Changes Made
- Change 1
- Change 2
- Change 3

## Testing
How did you test this?

## Screenshots (if applicable)
Before/after screenshots

## Checklist
- [ ] Tests pass
- [ ] Documentation updated
- [ ] Code follows style guidelines
- [ ] Commits are clean
```

### Review Process

1. Maintainers will review your PR
2. You may be asked to make changes
3. Once approved, PR will be merged
4. You'll be added to contributors list! 🎉

---

##  Code Style Guidelines

### Python Code

**Follow PEP 8:**
- Use 4 spaces for indentation
- Max line length: 88 characters (Black formatter default)
- Use meaningful variable names
- Add docstrings to functions

**Example:**
```python
def calculate_library_size_cv(counts: pd.DataFrame) -> float:
    """
    Calculate coefficient of variation for library sizes.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes × samples)
    
    Returns
    -------
    float
        Coefficient of variation (std/mean)
    
    Examples
    --------
    >>> counts = pd.DataFrame({'S1': [100, 200], 'S2': [150, 250]})
    >>> cv = calculate_library_size_cv(counts)
    >>> print(f"{cv:.2f}")
    0.12
    """
    library_sizes = counts.sum(axis=0)
    return library_sizes.std() / library_sizes.mean()
```

### R Code

**Follow Bioconductor style:**
- Use `<-` for assignment
- CamelCase for function names
- Meaningful variable names
- Roxygen2 documentation

**Example:**
```r
#' Calculate DEGs using DESeq2
#'
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @return DESeq2 results object
#' @export
runDESeq2Analysis <- function(counts, metadata) {
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata,
        design = ~ condition
    )
    dds <- DESeq(dds)
    return(results(dds))
}
```

### Shell Scripts

**Best practices:**
- Use `#!/bin/bash` shebang
- Quote variables: `"$variable"`
- Check exit codes
- Add comments

---

##  Testing Guidelines

### Writing Tests

**Good tests are:**
- Independent (can run in any order)
- Repeatable (same result every time)
- Fast (avoid slow operations)
- Clear (easy to understand what's tested)

**Example:**
```python
import pytest
import pandas as pd
from raptor.profiler import RNAseqDataProfiler

def test_library_size_calculation():
    """Test that library sizes are calculated correctly."""
    # Create test data
    counts = pd.DataFrame({
        'S1': [100, 200, 300],
        'S2': [150, 250, 350]
    })
    
    # Expected library sizes
    expected = pd.Series([600, 750], index=['S1', 'S2'])
    
    # Calculate
    profiler = RNAseqDataProfiler(counts)
    result = profiler.calculate_library_sizes()
    
    # Assert
    pd.testing.assert_series_equal(result, expected)

def test_handles_zero_inflation():
    """Test profiler handles highly zero-inflated data."""
    # Create zero-inflated data
    counts = pd.DataFrame({
        'S1': [0, 0, 0, 100],
        'S2': [0, 0, 0, 150]
    })
    
    profiler = RNAseqDataProfiler(counts)
    zero_pct = profiler.calculate_zero_percentage()
    
    assert zero_pct == 75.0  # 6 zeros out of 8 values
```

### Running Tests

```bash
# All tests
pytest tests/

# With coverage
pytest --cov=raptor tests/

# Verbose output
pytest -v tests/

# Specific test
pytest tests/test_profiler.py::test_library_size_calculation

# Stop on first failure
pytest -x tests/
```

---

##  Documentation Guidelines

### Docstring Format

Use NumPy style docstrings:

```python
def recommend_pipeline(profile: dict, priority: str = 'balanced') -> dict:
    """
    Recommend optimal pipeline based on data profile.
    
    This function analyzes data characteristics and matches them to
    pipeline strengths using a scoring system.
    
    Parameters
    ----------
    profile : dict
        Data profile containing statistical characteristics
    priority : str, optional
        Optimization priority: 'accuracy', 'speed', 'memory', or 'balanced'
        Default is 'balanced'
    
    Returns
    -------
    dict
        Recommendation with structure:
        {
            'pipeline_id': int,
            'pipeline_name': str,
            'score': float,
            'reasoning': list of str
        }
    
    Raises
    ------
    ValueError
        If priority is not one of the valid options
    
    Examples
    --------
    >>> profile = {'library_size_cv': 0.3, 'n_samples': 6}
    >>> rec = recommend_pipeline(profile, priority='accuracy')
    >>> print(rec['pipeline_name'])
    'STAR-RSEM-DESeq2'
    
    Notes
    -----
    The scoring system weighs different factors based on priority:
    - accuracy: Emphasizes sensitivity and precision
    - speed: Prioritizes fast methods
    - memory: Selects low-memory pipelines
    - balanced: Data-driven weighting
    
    See Also
    --------
    RNAseqDataProfiler : For generating profiles
    PipelineBenchmark : For validating recommendations
    """
    # Implementation here
    pass
```

### README Updates

When adding features:
1. Update main README.md
2. Add to appropriate section
3. Update table of contents if needed
4. Add example usage
5. Update Quick Start if it changes workflow

---

##  Recognition

### Contributors

All contributors will be:
- Listed in AUTHORS.md
- Mentioned in release notes
- Credited in documentation
- Added to Zenodo author list (for DOI)

### Significant Contributions

Major contributions may result in:
- Co-authorship on future papers
- Maintainer status
- Your name in the tool itself

---

##  Getting Help

### Questions?

- 💬 **GitHub Discussions**: For general questions
- 🐛 **GitHub Issues**: For bugs and feature requests
- 📧 **Email**: ayehbolouki1988@gmail.com for private matters

### Communication Guidelines

- Be respectful and professional
- Stay on topic
- Search before asking (question may be answered)
- Provide context and details
- Be patient - maintainers are volunteers

---

##  Code of Conduct

### Our Pledge

RAPTOR is committed to providing a welcoming, inclusive environment for all contributors regardless of:
- Background or identity
- Experience level
- Geographic location
- Institutional affiliation

### Expected Behavior

- Use welcoming and inclusive language
- Respect differing viewpoints
- Accept constructive criticism gracefully
- Focus on what's best for the community
- Show empathy toward others

### Unacceptable Behavior

- Harassment or discrimination
- Trolling or inflammatory comments
- Personal or political attacks
- Publishing others' private information
- Unprofessional conduct

### Enforcement

Violations may result in:
1. Warning
2. Temporary ban
3. Permanent ban

Report issues to: ayehbolouki1988@gmail.com

---

##  Contribution Priorities

### High Priority

- Adding new pipelines
- Improving recommendation algorithm
- Performance optimizations
- Bug fixes
- Documentation improvements

### Medium Priority

- Additional visualizations
- New metrics
- Extended format support
- Web interface

### Future Plans

- Single-cell RNA-seq support
- Machine learning enhancements
- Cloud deployment options
- Interactive dashboard

---

##  Resources

### Helpful Links

- [GitHub Flow Guide](https://guides.github.com/introduction/flow/)
- [Writing Good Commit Messages](https://chris.beams.io/posts/git-commit/)
- [PEP 8 Style Guide](https://www.python.org/dev/peps/pep-0008/)
- [Bioconductor Coding Style](https://bioconductor.org/developers/how-to/coding-style/)

### Learning Resources

- [Python Testing with pytest](https://docs.pytest.org/)
- [Git Basics](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics)
- [RNA-seq Analysis Methods](https://www.nature.com/articles/nprot.2016.095)

---

## 🙏 Thank You!

Every contribution, no matter how small, helps make RAPTOR better for the entire research community. Thank you for being part of this open science initiative!

**Let's make free science for everybody around the world!** 🦖

---

##  License

By contributing to RAPTOR, you agree that your contributions will be licensed under the MIT License.
