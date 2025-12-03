#!/usr/bin/env python3

"""
Quick Test Script for RAPTOR ML Recommender

Tests all components to ensure proper installation and functionality.

Author: Ayeh Bolouki
"""

import sys
from pathlib import Path

def test_imports():
    """Test if all required modules can be imported."""
    print("Testing imports...")
    
    try:
        import numpy as np
        print("  ✓ numpy")
    except ImportError:
        print("  ✗ numpy - install with: pip install numpy")
        return False
    
    try:
        import pandas as pd
        print("  ✓ pandas")
    except ImportError:
        print("  ✗ pandas - install with: pip install pandas")
        return False
    
    try:
        from sklearn.ensemble import RandomForestClassifier
        print("  ✓ scikit-learn")
    except ImportError:
        print("  ✗ scikit-learn - install with: pip install scikit-learn")
        return False
    
    try:
        import matplotlib.pyplot as plt
        print("  ✓ matplotlib")
    except ImportError:
        print("  ✗ matplotlib - install with: pip install matplotlib")
        return False
    
    try:
        import seaborn as sns
        print("  ✓ seaborn")
    except ImportError:
        print("  ✗ seaborn - install with: pip install seaborn")
        return False
    
    try:
        import joblib
        print("  ✓ joblib")
    except ImportError:
        print("  ✗ joblib - install with: pip install joblib")
        return False
    
    return True


def test_ml_modules():
    """Test custom ML modules."""
    print("\nTesting custom modules...")
    
    try:
        from raptor.ml_recommender import MLPipelineRecommender, FeatureExtractor
        print("  ✓ ml_recommender")
    except ImportError as e:
        print(f"  ✗ ml_recommender - {e}")
        return False
    
    try:
        from raptor.synthetic_benchmarks import SyntheticBenchmarkGenerator
        print("  ✓ synthetic_benchmarks")
    except ImportError as e:
        print(f"  ✗ synthetic_benchmarks - {e}")
        return False
    
    return True


def test_feature_extraction():
    """Test feature extraction."""
    print("\nTesting feature extraction...")
    
    from raptor.ml_recommender import FeatureExtractor
    
    # Create dummy profile
    profile = {
        'design': {
            'n_samples': 6,
            'n_genes': 20000,
            'n_conditions': 2,
            'samples_per_condition': 3,
            'is_paired': False
        },
        'library_stats': {
            'mean': 20000000,
            'median': 19500000,
            'cv': 0.15,
            'range': 5000000,
            'skewness': 0.2
        },
        'count_distribution': {
            'zero_pct': 40.0,
            'low_count_pct': 55.0,
            'mean': 1000,
            'median': 700,
            'variance': 50000
        },
        'expression_distribution': {
            'high_expr_genes': 2000,
            'medium_expr_genes': 6000,
            'low_expr_genes': 12000,
            'dynamic_range': 8.5
        },
        'biological_variation': {
            'bcv': 0.3,
            'dispersion_mean': 0.09,
            'dispersion_trend': 0.1,
            'outlier_genes': 500
        },
        'sequencing': {
            'total_reads': 120000000,
            'reads_per_gene': 6000,
            'depth_category': 'high'
        },
        'complexity': {
            'score': 65.0,
            'noise_level': 0.6,
            'signal_strength': 0.7
        }
    }
    
    try:
        features = FeatureExtractor.extract_features(profile)
        print(f"  ✓ Extracted {len(features.columns)} features")
        return True
    except Exception as e:
        print(f"  ✗ Feature extraction failed: {e}")
        return False


def test_data_generation():
    """Test synthetic data generation."""
    print("\nTesting data generation...")
    
    from raptor.synthetic_benchmarks import SyntheticBenchmarkGenerator
    import tempfile
    import shutil
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp())
    
    try:
        generator = SyntheticBenchmarkGenerator(n_datasets=5, seed=42)
        summary = generator.generate_benchmarks(str(temp_dir))
        
        print(f"  ✓ Generated {summary['n_datasets']} datasets")
        
        # Check files exist
        dataset_dirs = list(temp_dir.glob('dataset_*'))
        if len(dataset_dirs) == 5:
            print(f"  ✓ All dataset directories created")
        else:
            print(f"  ✗ Expected 5 directories, found {len(dataset_dirs)}")
            return False
        
        # Check file contents
        first_dataset = dataset_dirs[0]
        if (first_dataset / 'data_profile.json').exists():
            print(f"  ✓ Profile files created")
        else:
            print(f"  ✗ Profile files missing")
            return False
        
        if (first_dataset / 'benchmark_results.json').exists():
            print(f"  ✓ Benchmark files created")
        else:
            print(f"  ✗ Benchmark files missing")
            return False
        
        return True
        
    except Exception as e:
        print(f"  ✗ Data generation failed: {e}")
        return False
    
    finally:
        # Cleanup
        shutil.rmtree(temp_dir)


def test_model_training():
    """Test model training."""
    print("\nTesting model training...")
    
    from raptor.ml_recommender import MLPipelineRecommender
    from raptor.synthetic_benchmarks import generate_training_data
    import tempfile
    import shutil
    
    # Create temporary directories
    data_dir = Path(tempfile.mkdtemp())
    model_dir = Path(tempfile.mkdtemp())
    
    try:
        # Generate small dataset
        print("  Generating 20 training samples...")
        generate_training_data(n_datasets=20, output_dir=str(data_dir), seed=42)
        
        # Train model
        print("  Training Random Forest...")
        recommender = MLPipelineRecommender(model_type='random_forest')
        results = recommender.train_from_benchmarks(str(data_dir))
        
        print(f"  ✓ Training completed")
        print(f"    Test accuracy: {results['test_score']:.3f}")
        
        # Save model
        recommender.save_model(str(model_dir))
        print(f"  ✓ Model saved")
        
        # Load model
        recommender2 = MLPipelineRecommender(model_type='random_forest')
        recommender2.load_model(str(model_dir))
        print(f"  ✓ Model loaded")
        
        return True
        
    except Exception as e:
        print(f"  ✗ Model training failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Cleanup
        shutil.rmtree(data_dir)
        shutil.rmtree(model_dir)


def test_prediction():
    """Test prediction."""
    print("\nTesting prediction...")
    
    from raptor.ml_recommender import MLPipelineRecommender, FeatureExtractor
    from raptor.synthetic_benchmarks import SyntheticBenchmarkGenerator, generate_training_data
    import tempfile
    import shutil
    
    # Create temporary directories
    data_dir = Path(tempfile.mkdtemp())
    model_dir = Path(tempfile.mkdtemp())
    
    try:
        # Generate and train
        print("  Training model...")
        generate_training_data(n_datasets=20, output_dir=str(data_dir), seed=42)
        
        recommender = MLPipelineRecommender(model_type='random_forest')
        recommender.train_from_benchmarks(str(data_dir))
        
        # Generate test profile
        print("  Generating test profile...")
        generator = SyntheticBenchmarkGenerator(n_datasets=1, seed=999)
        profile = generator._generate_profile()
        
        # Make prediction
        print("  Making prediction...")
        recommendation = recommender.recommend(profile, top_k=3)
        
        print(f"  ✓ Prediction successful")
        print(f"    Recommended: Pipeline {recommendation['pipeline_id']}")
        print(f"    Confidence: {recommendation['confidence']:.1%}")
        print(f"    Alternatives: {len(recommendation['alternatives'])}")
        
        # Verify structure
        required_keys = ['pipeline_id', 'pipeline_name', 'confidence', 'reasons']
        if all(key in recommendation for key in required_keys):
            print(f"  ✓ Recommendation structure valid")
        else:
            print(f"  ✗ Recommendation structure invalid")
            return False
        
        return True
        
    except Exception as e:
        print(f"  ✗ Prediction failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Cleanup
        shutil.rmtree(data_dir)
        shutil.rmtree(model_dir)


def main():
    """Run all tests."""
    print("=" * 70)
    print("RAPTOR ML Recommender - System Test")
    print("=" * 70)
    
    tests = [
        ("Import dependencies", test_imports),
        ("Custom modules", test_ml_modules),
        ("Feature extraction", test_feature_extraction),
        ("Data generation", test_data_generation),
        ("Model training", test_model_training),
        ("Prediction", test_prediction)
    ]
    
    results = []
    
    for test_name, test_func in tests:
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"\n✗ {test_name} crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    for test_name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{status:8s} {test_name}")
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    print(f"\nResults: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! System is ready to use.")
        return 0
    else:
        print("\n⚠️  Some tests failed. Check errors above.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
