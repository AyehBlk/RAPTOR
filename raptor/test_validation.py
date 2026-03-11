#!/usr/bin/env python3
"""
Test Script for RAPTOR Validation

Tests all validation functions to ensure they work correctly.

Author: Ayeh Bolouki
Version: 2.2.0
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from validation import (
    validate_count_matrix,
    validate_metadata,
    validate_group_column,
    validate_file_path,
    validate_directory_path,
    validate_numeric_range,
    validate_probability,
    validate_positive_integer
)
from errors import (
    ValidationError,
    check_file_exists,
    validate_output_writable
)


def test_count_matrix_validation():
    """Test count matrix validation."""
    print("\n1. Testing Count Matrix Validation")
    print("-" * 50)
    
    # Test 1: Valid count matrix
    counts = pd.DataFrame(
        np.random.poisson(100, (1000, 6)),
        index=[f"gene_{i}" for i in range(1000)],
        columns=[f"sample_{i}" for i in range(6)]
    )
    try:
        validate_count_matrix(counts)
        print("   ✓ Valid count matrix accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 2: Negative counts (should fail)
    counts_neg = counts.copy()
    counts_neg.iloc[0, 0] = -1
    try:
        validate_count_matrix(counts_neg)
        print("   ✗ FAILED: Negative counts NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Negative counts rejected correctly")
    
    # Test 3: Too few genes (should fail)
    counts_small = counts.iloc[:10, :]
    try:
        validate_count_matrix(counts_small, min_genes=100)
        print("   ✗ FAILED: Too few genes NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Too few genes rejected correctly")
    
    # Test 4: Too few samples (should fail)
    counts_few_samples = counts.iloc[:, :1]
    try:
        validate_count_matrix(counts_few_samples, min_samples=2)
        print("   ✗ FAILED: Too few samples NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Too few samples rejected correctly")
    
    # Test 5: Duplicate gene names (should fail)
    counts_dup = counts.copy()
    counts_dup.index = [f"gene_{i % 500}" for i in range(1000)]
    try:
        validate_count_matrix(counts_dup)
        print("   ✗ FAILED: Duplicate gene names NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Duplicate gene names rejected correctly")
    
    # Test 6: NaN values (should fail)
    counts_nan = counts.copy()
    counts_nan.iloc[0, 0] = np.nan
    try:
        validate_count_matrix(counts_nan)
        print("   ✗ FAILED: NaN values NOT rejected!")
        return False
    except ValueError:
        print("   ✓ NaN values rejected correctly")
    
    return True


def test_metadata_validation():
    """Test metadata validation."""
    print("\n2. Testing Metadata Validation")
    print("-" * 50)
    
    # Create valid count matrix
    counts = pd.DataFrame(
        np.random.poisson(100, (1000, 6)),
        index=[f"gene_{i}" for i in range(1000)],
        columns=[f"sample_{i}" for i in range(6)]
    )
    
    # Test 1: Valid metadata
    metadata = pd.DataFrame({
        'sample_id': [f'sample_{i}' for i in range(6)],
        'condition': ['A', 'A', 'A', 'B', 'B', 'B'],
        'batch': [1, 1, 2, 1, 2, 2]
    })
    try:
        validate_metadata(metadata, counts)
        print("   ✓ Valid metadata accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 2: Missing sample_id column (should fail)
    metadata_no_id = metadata.drop('sample_id', axis=1)
    try:
        validate_metadata(metadata_no_id)
        print("   ✗ FAILED: Missing sample_id NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Missing sample_id rejected correctly")
    
    # Test 3: Duplicate sample IDs (should fail)
    metadata_dup = metadata.copy()
    metadata_dup.iloc[1, 0] = metadata_dup.iloc[0, 0]
    try:
        validate_metadata(metadata_dup)
        print("   ✗ FAILED: Duplicate sample IDs NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Duplicate sample IDs rejected correctly")
    
    # Test 4: Samples in counts missing from metadata (should fail)
    metadata_missing = metadata.iloc[:4, :]
    try:
        validate_metadata(metadata_missing, counts)
        print("   ✗ FAILED: Missing samples NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Missing samples rejected correctly")
    
    return True


def test_group_column_validation():
    """Test group column validation."""
    print("\n3. Testing Group Column Validation")
    print("-" * 50)
    
    # Test 1: Valid group column
    metadata = pd.DataFrame({
        'sample_id': [f'sample_{i}' for i in range(6)],
        'condition': ['A', 'A', 'A', 'B', 'B', 'B']
    })
    try:
        validate_group_column(metadata, 'condition', min_groups=2, min_samples_per_group=3)
        print("   ✓ Valid group column accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 2: Non-existent column (should fail)
    try:
        validate_group_column(metadata, 'nonexistent')
        print("   ✗ FAILED: Non-existent column NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Non-existent column rejected correctly")
    
    # Test 3: Too few groups (should fail)
    metadata_one_group = pd.DataFrame({
        'sample_id': [f'sample_{i}' for i in range(6)],
        'condition': ['A'] * 6
    })
    try:
        validate_group_column(metadata_one_group, 'condition', min_groups=2)
        print("   ✗ FAILED: Too few groups NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Too few groups rejected correctly")
    
    # Test 4: Too few samples per group (should fail)
    metadata_small_groups = pd.DataFrame({
        'sample_id': [f'sample_{i}' for i in range(6)],
        'condition': ['A', 'A', 'B', 'B', 'C', 'C']
    })
    try:
        validate_group_column(metadata_small_groups, 'condition', min_samples_per_group=3)
        print("   ✗ FAILED: Too few samples per group NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Too few samples per group rejected correctly")
    
    return True


def test_numerical_validation():
    """Test numerical validation."""
    print("\n4. Testing Numerical Validation")
    print("-" * 50)
    
    # Test 1: Valid positive integer
    try:
        validate_positive_integer(10, 'n_samples')
        print("   ✓ Valid positive integer accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 2: Negative integer (should fail)
    try:
        validate_positive_integer(-5, 'n_samples')
        print("   ✗ FAILED: Negative integer NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Negative integer rejected correctly")
    
    # Test 3: Zero (should fail)
    try:
        validate_positive_integer(0, 'n_samples')
        print("   ✗ FAILED: Zero NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Zero rejected correctly")
    
    # Test 4: Valid probability
    try:
        validate_probability(0.05, 'alpha')
        print("   ✓ Valid probability accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 5: Invalid probability > 1 (should fail)
    try:
        validate_probability(1.5, 'alpha')
        print("   ✗ FAILED: Probability > 1 NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Probability > 1 rejected correctly")
    
    # Test 6: Invalid probability < 0 (should fail)
    try:
        validate_probability(-0.1, 'alpha')
        print("   ✗ FAILED: Probability < 0 NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Probability < 0 rejected correctly")
    
    # Test 7: Numeric range validation
    try:
        validate_numeric_range(5, 'value', min_val=0, max_val=10)
        print("   ✓ Valid numeric range accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 8: Out of range (should fail)
    try:
        validate_numeric_range(15, 'value', min_val=0, max_val=10)
        print("   ✗ FAILED: Out of range value NOT rejected!")
        return False
    except ValueError:
        print("   ✓ Out of range value rejected correctly")
    
    return True


def test_file_path_validation():
    """Test file path validation."""
    print("\n5. Testing File Path Validation")
    print("-" * 50)
    
    # Test 1: Non-existent file (should fail)
    try:
        validate_file_path("/nonexistent/file.csv", must_exist=True)
        print("   ✗ FAILED: Missing file NOT rejected!")
        return False
    except FileNotFoundError:
        print("   ✓ Missing file rejected correctly")
    
    # Test 2: Path without existence check (should pass)
    try:
        path = validate_file_path("output.csv", must_exist=False)
        print("   ✓ Path without existence check accepted")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return False
    
    # Test 3: Create a temp file and validate it
    temp_file = Path("/tmp/test_raptor_validation.txt")
    temp_file.write_text("test")
    try:
        path = validate_file_path(temp_file, must_exist=True)
        print("   ✓ Existing file accepted")
        temp_file.unlink()  # Clean up
    except Exception as e:
        print(f"   ✗ Error: {e}")
        if temp_file.exists():
            temp_file.unlink()
        return False
    
    # Test 4: Directory path creation
    temp_dir = Path("/tmp/test_raptor_dir")
    try:
        path = validate_directory_path(temp_dir, create_if_missing=True)
        print("   ✓ Directory created successfully")
        temp_dir.rmdir()  # Clean up
    except Exception as e:
        print(f"   ✗ Error: {e}")
        if temp_dir.exists():
            temp_dir.rmdir()
        return False
    
    return True


def run_all_tests():
    """Run all validation tests."""
    print("\n" + "=" * 60)
    print("RAPTOR v2.2.0 - Validation Test Suite")
    print("=" * 60)
    
    results = []
    
    # Run each test
    results.append(("Count Matrix", test_count_matrix_validation()))
    results.append(("Metadata", test_metadata_validation()))
    results.append(("Group Column", test_group_column_validation()))
    results.append(("Numerical", test_numerical_validation()))
    results.append(("File Path", test_file_path_validation()))
    
    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    for test_name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{test_name:20s}: {status}")
    
    all_passed = all(passed for _, passed in results)
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✓ ALL TESTS PASSED!")
        print("=" * 60)
        return 0
    else:
        print("✗ SOME TESTS FAILED!")
        print("=" * 60)
        return 1


if __name__ == '__main__':
    exit_code = run_all_tests()
    sys.exit(exit_code)
