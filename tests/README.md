# RAPTOR Tests

Comprehensive test suite for RAPTOR RNA-seq analysis tool.

## Test Structure

```
tests/
├── README.md                 # This file
├── test_profiler.py         # Unit tests for RNAseqDataProfiler
├── test_recommender.py      # Unit tests for PipelineRecommender
└── test_pipelines.sh        # Integration tests for pipelines
```

## Running Tests

### Quick Test (Recommended First)
```bash
# Run all Python unit tests
cd tests/
pytest -v

# Run specific test file
pytest test_profiler.py -v
pytest test_recommender.py -v
```

### Full Test Suite
```bash
# Run all tests including integration
chmod +x test_pipelines.sh
./test_pipelines.sh
```

### With Coverage
```bash
# Generate coverage report
pytest --cov=raptor --cov-report=html

# View coverage
xdg-open htmlcov/index.html
```

## Test Files

### test_profiler.py
Tests for the RNAseqDataProfiler class.

**Coverage:**
- Data loading and validation
- BCV calculation
- Sequencing depth analysis
- Zero inflation metrics
- Library size CV
- Outlier detection
- Quality flag generation
- Edge cases (negative values, missing data, etc.)

**Run:**
```bash
pytest test_profiler.py -v

# Run specific test class
pytest test_profiler.py::TestRNAseqDataProfiler -v

# Run specific test
pytest test_profiler.py::TestRNAseqDataProfiler::test_bcv_calculation -v
```

### test_recommender.py
Tests for the PipelineRecommender class.

**Coverage:**
- Recommendation generation
- Scoring algorithms
- Pipeline ranking
- Custom weights
- Resource constraints
- Edge cases

**Run:**
```bash
pytest test_recommender.py -v

# With verbose output
pytest test_recommender.py -vv

# Stop on first failure
pytest test_recommender.py -x
```

### test_pipelines.sh
Integration tests for pipeline installation and execution.

**Coverage:**
- RAPTOR installation
- Bioinformatics tools availability
- R/Bioconductor packages
- Configuration files
- Data simulation
- End-to-end profiling workflow
- Output file validation

**Run:**
```bash
chmod +x test_pipelines.sh
./test_pipelines.sh

# View test log
cat pipeline_tests_*/test.log
```

## Prerequisites

### For Python Tests:
```bash
# Install pytest and dependencies
pip install pytest pytest-cov

# Install RAPTOR
pip install raptor-rnaseq

# Or install in development mode
pip install -e .
```

### For Integration Tests:
```bash
# Basic tools (check with test script)
# - RAPTOR
# - Python 3.8+
# - Bioinformatics tools (STAR, Salmon, etc.)
```

## Test Markers

Tests use pytest markers for categorization:

```bash
# Run only fast tests
pytest -m "not slow"

# Run only slow tests
pytest -m slow

# Run specific category
pytest -m "profiler"
```

## Continuous Integration

These tests are run automatically on:
- Every push to main branch
- All pull requests
- Daily scheduled runs

See `.github/workflows/` for CI configuration.

## Writing New Tests

### Adding a Test to test_profiler.py:
```python
class TestRNAseqDataProfiler:
    def test_my_new_feature(self, simple_counts, simple_metadata):
        """Test description"""
        profiler = RNAseqDataProfiler(simple_counts, simple_metadata)
        result = profiler.my_new_method()
        
        assert result is not None
        assert result > 0
```

### Adding a Test to test_recommender.py:
```python
class TestPipelineRecommender:
    def test_my_recommendation_feature(self, sample_profile_small):
        """Test description"""
        recommender = PipelineRecommender()
        recommendations = recommender.recommend(sample_profile_small)
        
        assert len(recommendations) > 0
        assert recommendations[0]['score'] > 0.5
```

### Adding an Integration Test:
```bash
# In test_pipelines.sh, add new test section:

print_header "Test X: My New Feature"

print_info "Testing new feature..."

if my_test_command; then
    test_result 0 "Feature works correctly"
else
    test_result 1 "Feature failed"
fi
```

## Test Data

Tests use:
- **Fixtures:** Defined in test files (@pytest.fixture)
- **Simulated data:** Generated programmatically
- **Mock objects:** For external dependencies

## Troubleshooting

### Problem: "RAPTOR not installed"
```bash
pip install raptor-rnaseq
# Or for development:
pip install -e .
```

### Problem: Import errors
```bash
# Install test dependencies
pip install -r requirements-dev.txt

# Or install individually
pip install pytest pytest-cov numpy pandas scipy
```

### Problem: Tests fail unexpectedly
```bash
# Run with verbose output
pytest -vv

# Show print statements
pytest -s

# Show full traceback
pytest --tb=long
```

### Problem: Slow tests timeout
```bash
# Skip slow tests
pytest -m "not slow"

# Increase timeout
pytest --timeout=300
```

## Test Coverage Goals

Target coverage levels:
- **Overall:** >80%
- **Core modules:** >90%
- **Critical functions:** 100%

Current coverage:
```bash
pytest --cov=raptor --cov-report=term-missing
```

## Best Practices

1. **Write tests first** (TDD approach)
2. **Keep tests independent** (no test depends on another)
3. **Use descriptive names** (test_should_do_something)
4. **Test edge cases** (empty data, extreme values)
5. **Mock external dependencies** (file I/O, network calls)
6. **Use fixtures** for common setup
7. **Keep tests fast** (mark slow tests)
8. **Document complex tests**

## Example Test Session

```bash
# Navigate to tests directory
cd tests/

# Run all tests with coverage
pytest --cov=raptor --cov-report=html --cov-report=term

# View results
cat coverage_report.txt

# Run integration tests
./test_pipelines.sh

# Check all passed
echo $?  # Should be 0 if all passed
```

## Contributing Tests

When contributing to RAPTOR:

1. **Write tests** for new features
2. **Ensure all tests pass** before submitting PR
3. **Aim for >80% coverage** for new code
4. **Document test purpose** in docstrings
5. **Use existing patterns** from current tests

## Resources

- **pytest documentation:** https://docs.pytest.org/
- **RAPTOR documentation:** [../docs/](../docs/)
- **Contributing guide:** [../CONTRIBUTING.md](../CONTRIBUTING.md)

## Support

Questions about tests?
- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Discussions:** https://github.com/AyehBlk/RAPTOR/discussions
- **Email:** ayehbolouki1988@gmail.com

---

**Testing ensures quality! 🧪🦖**
