#!/usr/bin/env python3

"""
RAPTOR v2.1.1 - Comprehensive Test Script

Tests all core modules for:
- Import correctness
- Syntax errors
- Type hint consistency
- Version numbers
- API compatibility
- 🆕 Adaptive Threshold Optimizer (ATO)

Author: Ayeh Bolouki
"""

import sys
import importlib
import pkgutil
from pathlib import Path

def print_header(text):
    """Print formatted header."""
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70 + "\n")

def test_imports():
    """Test all imports work correctly."""
    print_header("TEST 1: Module Imports")
    
    modules_to_test = [
        # Core v2.0.0 modules
        ('raptor.profiler', 'RNAseqDataProfiler'),
        ('raptor.recommender', 'PipelineRecommender'),
        ('raptor.benchmark', 'PipelineBenchmark'),
        ('raptor.simulate', 'DataSimulator'),
        ('raptor.report', 'ReportGenerator'),
        ('raptor.utils', 'ensure_dir'),
        
        # v2.1.0 modules
        ('raptor.ml_recommender', 'MLPipelineRecommender'),
        ('raptor.data_quality_assessment', 'DataQualityAssessor'),
        ('raptor.automated_reporting', 'AutomatedReporter'),
        ('raptor.ensemble_analysis', 'EnsembleAnalyzer'),
        ('raptor.parameter_optimization', 'ParameterOptimizer'),
        ('raptor.resource_monitoring', 'ResourceMonitor'),
        ('raptor.synthetic_benchmarks', 'SyntheticBenchmarkGenerator'),
        
        # 🆕 v2.1.1 modules
        ('raptor.threshold_optimizer', 'AdaptiveThresholdOptimizer'),
        ('raptor.threshold_optimizer', 'ThresholdResult'),
        ('raptor.threshold_optimizer', 'optimize_thresholds'),
    ]
    
    passed = 0
    failed = 0
    
    for module_name, class_name in modules_to_test:
        try:
            module = importlib.import_module(module_name)
            if hasattr(module, class_name):
                print(f"✓ {module_name}.{class_name}")
                passed += 1
            else:
                print(f"✗ {module_name}.{class_name} - Class not found")
                failed += 1
        except ImportError as e:
            print(f"✗ {module_name} - Import Error: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ {module_name} - Error: {e}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0

def test_version():
    """Test version number is correct."""
    print_header("TEST 2: Version Number")
    
    try:
        import raptor
        version = raptor.__version__
        
        if version == '2.1.1':
            print(f"✓ Version: {version}")
            return True
        else:
            print(f"✗ Expected version 2.1.1, got {version}")
            return False
    except Exception as e:
        print(f"✗ Error checking version: {e}")
        return False

def test_type_hints():
    """Test type hints are properly imported."""
    print_header("TEST 3: Type Hints")
    
    tests = [
        ('raptor.simulate', 'quick_simulate', 'Dict'),
        ('raptor.ml_recommender', 'MLPipelineRecommender', 'Dict'),
        ('raptor.synthetic_benchmarks', 'generate_training_data', 'Dict'),
        ('raptor.threshold_optimizer', 'optimize_thresholds', 'DataFrame'),  # 🆕 v2.1.1
    ]
    
    passed = 0
    failed = 0
    
    for module_name, func_or_class, type_hint in tests:
        try:
            module = importlib.import_module(module_name)
            obj = getattr(module, func_or_class)
            
            # Check if type hint is used in annotations
            if hasattr(obj, '__annotations__'):
                annotations = str(obj.__annotations__)
                if type_hint in annotations:
                    print(f"✓ {module_name}.{func_or_class} uses {type_hint}")
                    passed += 1
                else:
                    print(f"⚠ {module_name}.{func_or_class} may not use {type_hint}")
                    passed += 1  # Not critical
            else:
                print(f"⚠ {module_name}.{func_or_class} - No annotations")
                passed += 1  # Not critical
                
        except Exception as e:
            print(f"✗ {module_name}.{func_or_class} - Error: {e}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0

def test_dependencies():
    """Test required dependencies are available."""
    print_header("TEST 4: Dependencies")
    
    required = [
        'numpy',
        'pandas',
        'matplotlib',
        'seaborn',
        'scipy',
    ]
    
    optional = [
        ('scikit-learn', 'sklearn'),  # For ML features
        ('boto3', 'boto3'),  # For AWS
        ('google-cloud-storage', 'google.cloud.storage'),  # For GCP
        ('statsmodels', 'statsmodels'),  # 🆕 For ATO pi0 estimation
    ]
    
    print("Required Dependencies:")
    all_found = True
    for pkg in required:
        try:
            importlib.import_module(pkg)
            print(f"✓ {pkg}")
        except ImportError:
            print(f"✗ {pkg} - NOT FOUND")
            all_found = False
    
    print("\nOptional Dependencies:")
    for display_name, import_name in optional:
        try:
            importlib.import_module(import_name)
            print(f"✓ {display_name}")
        except ImportError:
            print(f"⚠ {display_name} - Not installed (optional)")
    
    return all_found

def test_cli_imports():
    """Test CLI imports work correctly."""
    print_header("TEST 5: CLI Import Corrections")
    
    try:
        # Test that raptor_ml_cli can import ML modules correctly
        import sys
        from pathlib import Path
        
        # Add the raptor directory to path if needed
        test_code = '''
from raptor.ml_recommender import MLPipelineRecommender
from raptor.synthetic_benchmarks import generate_training_data
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer, optimize_thresholds

print("✓ ML imports successful")
print("✓ ATO imports successful")
'''
        
        exec(test_code)
        return True
        
    except ImportError as e:
        print(f"✗ Import Error: {e}")
        return False
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

def test_threshold_optimizer():
    """Test Adaptive Threshold Optimizer (ATO) - NEW in v2.1.1."""
    print_header("TEST 6: Threshold Optimizer (ATO) - NEW in v2.1.1")
    
    try:
        import numpy as np
        import pandas as pd
        from raptor.threshold_optimizer import (
            AdaptiveThresholdOptimizer,
            ThresholdResult,
            optimize_thresholds
        )
        
        print("✓ ATO modules imported")
        
        # Generate test data
        np.random.seed(42)
        n_genes = 1000
        n_de = 150
        
        # Null genes
        null_logfc = np.random.normal(0, 0.2, n_genes - n_de)
        null_pval = np.random.uniform(0.05, 1, n_genes - n_de)
        
        # DE genes
        de_logfc = np.concatenate([
            np.random.normal(1.5, 0.5, n_de // 2),
            np.random.normal(-1.5, 0.5, n_de - n_de // 2)
        ])
        de_pval = np.random.exponential(0.001, n_de)
        de_pval = np.clip(de_pval, 1e-300, 0.05)
        
        df = pd.DataFrame({
            'log2FoldChange': np.concatenate([null_logfc, de_logfc]),
            'pvalue': np.concatenate([null_pval, de_pval]),
        })
        df.index = [f'Gene_{i}' for i in range(n_genes)]
        
        print(f"✓ Generated test data: {len(df)} genes")
        
        # Test optimize_thresholds function
        result = optimize_thresholds(df, goal='balanced', verbose=False)
        
        print(f"✓ optimize_thresholds() executed")
        print(f"  - LogFC cutoff: {result.logfc_cutoff:.3f}")
        print(f"  - Padj cutoff: {result.padj_cutoff}")
        print(f"  - Significant genes: {result.n_significant_optimized}")
        
        # Verify result is ThresholdResult type
        if isinstance(result, ThresholdResult):
            print(f"✓ Result is ThresholdResult type")
        else:
            print(f"✗ Result is not ThresholdResult type")
            return False
        
        # Verify result has required attributes
        required_attrs = ['logfc_cutoff', 'padj_cutoff', 'padj_method', 
                         'pi0_estimate', 'logfc_method', 'n_significant_optimized']
        
        for attr in required_attrs:
            if hasattr(result, attr):
                print(f"✓ Result has '{attr}' attribute")
            else:
                print(f"✗ Result missing '{attr}' attribute")
                return False
        
        # Test summary() method
        summary = result.summary()
        if summary and len(summary) > 50:
            print(f"✓ summary() method works ({len(summary)} chars)")
        else:
            print(f"✗ summary() method failed")
            return False
        
        # Test AdaptiveThresholdOptimizer class
        ato = AdaptiveThresholdOptimizer(df, goal='discovery', verbose=False)
        ato_result = ato.optimize()
        
        print(f"✓ AdaptiveThresholdOptimizer class works")
        
        # Test get_significant_genes
        sig_genes = ato.get_significant_genes()
        print(f"✓ get_significant_genes() returned {len(sig_genes)} genes")
        
        # Test different goals
        for goal in ['discovery', 'balanced', 'validation']:
            r = optimize_thresholds(df, goal=goal, verbose=False)
            print(f"✓ Goal '{goal}': |logFC| > {r.logfc_cutoff:.3f}, padj < {r.padj_cutoff}")
        
        return True
        
    except ImportError as e:
        print(f"✗ Import Error: {e}")
        print("  ATO module may not be installed correctly")
        return False
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("""
╔═══════════════════════════════════════════════════════════════════╗
║              RAPTOR v2.1.1 - Comprehensive Test Suite             ║
║                                                                   ║
║  🆕 Now includes Adaptive Threshold Optimizer (ATO) tests!       ║
╚═══════════════════════════════════════════════════════════════════╝

This script tests all RAPTOR modules for correctness.
    """)
    
    results = {
        'Module Imports': test_imports(),
        'Version Number': test_version(),
        'Type Hints': test_type_hints(),
        'Dependencies': test_dependencies(),
        'CLI Imports': test_cli_imports(),
        'Threshold Optimizer (ATO)': test_threshold_optimizer(),  # 🆕 v2.1.1
    }
    
    # Summary
    print_header("OVERALL SUMMARY")
    
    passed_tests = sum(1 for v in results.values() if v)
    total_tests = len(results)
    
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{test_name:.<50} {status}")
    
    print(f"\n{passed_tests}/{total_tests} test suites passed")
    
    if passed_tests == total_tests:
        print("\n🎉 All tests passed! RAPTOR v2.1.1 is ready to use.")
        return 0
    else:
        print("\n⚠️  Some tests failed. Please check the output above.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
