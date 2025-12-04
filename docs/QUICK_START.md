# 🦖 RAPTOR ML Upgrade - Quick Start Guide

## What's Included

This package contains a complete **Machine Learning-based Pipeline Recommendation System** for RAPTOR:

### Core Files (8 total)

1. **ml_recommender.py** - Main ML recommender implementation
2. **synthetic_benchmarks.py** - Training data generator
3. **example_ml_workflow.py** - Complete demo workflow
4. **raptor_ml_cli.py** - Enhanced command-line interface
5. **test_ml_system.py** - Comprehensive test suite
6. **requirements_ml.txt** - Python dependencies
7. **ML_RECOMMENDER_README.md** - Full documentation
8. **IMPLEMENTATION_SUMMARY.md** - Technical overview (this file)

## Installation (5 Minutes)

### Step 1: Install Dependencies

```bash
pip install numpy pandas scikit-learn scipy matplotlib seaborn joblib
```

Or use the requirements file:
```bash
pip install -r requirements_ml.txt
```

### Step 2: Verify Installation

```bash
python test_ml_system.py
```

Expected output:
```
✓ PASS Import dependencies
✓ PASS Custom modules  
✓ PASS Feature extraction
✓ PASS Data generation
✓ PASS Model training
✓ PASS Prediction

🎉 All tests passed! System is ready to use.
```

## First Run (2 Minutes)

### See It In Action

```bash
python example_ml_workflow.py --n-datasets 50
```

This will:
1. Generate 50 synthetic benchmark datasets
2. Train a Random Forest model
3. Evaluate performance
4. Create visualizations
5. Test predictions

Output includes:
- Training accuracy: ~88%
- Test accuracy: ~85%
- Figures in `figures/` directory
- Trained model in `models/` directory

## Basic Usage

### 1. Train a Model

```python
from ml_recommender import MLPipelineRecommender

# Create and train
recommender = MLPipelineRecommender(model_type='random_forest')
results = recommender.train_from_benchmarks('ml_training_data/')

print(f"Accuracy: {results['test_score']:.3f}")

# Save for later
recommender.save_model('my_models/')
```

### 2. Make Predictions

```python
from ml_recommender import MLPipelineRecommender

# Load trained model
recommender = MLPipelineRecommender()
recommender.load_model('my_models/')

# Get recommendation (profile is from RNAseqDataProfiler)
recommendation = recommender.recommend(profile, top_k=3)

print(f"Best: Pipeline {recommendation['pipeline_id']}")
print(f"Confidence: {recommendation['confidence']:.1%}")
```

### 3. Use with RAPTOR

```bash
# Profile your data with ML recommendation
python raptor_ml_cli.py profile \
    --counts your_data.csv \
    --use-ml \
    --ml-model models/
```

## Integration with Your RAPTOR Installation

### Option A: Standalone Usage

Use the files as-is alongside your existing RAPTOR installation:

```python
# In your scripts
from ml_recommender import MLPipelineRecommender
from raptor.profiler import RNAseqDataProfiler

# Profile data
counts = pd.read_csv('data.csv', index_col=0)
profiler = RNAseqDataProfiler(counts)
profile = profiler.run_full_profile()

# Get ML recommendation
recommender = MLPipelineRecommender()
recommender.load_model('models/')
rec = recommender.recommend(profile)
```

### Option B: Full Integration

1. Copy `ml_recommender.py` and `synthetic_benchmarks.py` to your RAPTOR package:
```bash
cp ml_recommender.py /path/to/raptor/raptor/
cp synthetic_benchmarks.py /path/to/raptor/raptor/
```

2. Update `raptor/__init__.py`:
```python
from raptor.ml_recommender import MLPipelineRecommender
```

3. Enhance existing `recommender.py`:
```python
# In raptor/recommender.py
try:
    from raptor.ml_recommender import MLPipelineRecommender
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False

class PipelineRecommender:
    def recommend(self, use_ml=False, ml_model_path=None):
        if use_ml and ML_AVAILABLE and ml_model_path:
            ml_rec = MLPipelineRecommender()
            ml_rec.load_model(ml_model_path)
            return ml_rec.recommend(self.profile)
        else:
            # Fall back to rule-based
            return self._rule_based_recommend()
```

## Directory Structure

After running the workflow, you'll have:

```
raptor_ml_upgrade/
├── ml_recommender.py              # Core ML module
├── synthetic_benchmarks.py        # Data generator
├── example_ml_workflow.py         # Demo workflow
├── raptor_ml_cli.py              # CLI interface
├── test_ml_system.py             # Tests
├── requirements_ml.txt           # Dependencies
├── ML_RECOMMENDER_README.md      # Full docs
└── IMPLEMENTATION_SUMMARY.md     # This file

Generated during use:
├── ml_training_data/             # Training datasets
│   ├── dataset_0000/
│   │   ├── data_profile.json
│   │   └── benchmark_results.json
│   ├── dataset_0001/
│   └── ...
├── models/                       # Trained models
│   ├── ml_recommender_random_forest.pkl
│   ├── scaler_random_forest.pkl
│   └── metadata_random_forest.json
└── figures/                      # Visualizations
    ├── confusion_matrix.png
    ├── feature_importance.png
    └── pipeline_metrics.png
```

## Common Use Cases

### Use Case 1: Quick Recommendation

```bash
# One-liner to get recommendation
python raptor_ml_cli.py profile --counts data.csv --use-ml
```

### Use Case 2: Train on Your Benchmarks

```bash
# You have real benchmark results
python raptor_ml_cli.py train-ml \
    --benchmark-dir my_benchmarks/ \
    --output production_models/
```

### Use Case 3: Compare Models

```bash
# Compare Random Forest vs Gradient Boosting
python example_ml_workflow.py --compare-models
```

### Use Case 4: Batch Processing

```python
from ml_recommender import MLPipelineRecommender
from raptor.profiler import RNAseqDataProfiler
import pandas as pd
from pathlib import Path

# Load model once
recommender = MLPipelineRecommender()
recommender.load_model('models/')

# Process multiple datasets
for data_file in Path('datasets/').glob('*.csv'):
    counts = pd.read_csv(data_file, index_col=0)
    profiler = RNAseqDataProfiler(counts)
    profile = profiler.run_full_profile()
    
    rec = recommender.recommend(profile)
    print(f"{data_file.name}: Pipeline {rec['pipeline_id']}")
```

## Performance Tips

### For Speed
```python
# Use multiprocessing
recommender.model.n_jobs = -1

# Reduce feature set (if needed)
# Focus on top 15 most important features
```

### For Accuracy
```python
# Train on more data
generate_training_data(n_datasets=500)

# Use ensemble
rf_rec = rf_model.recommend(profile)
gb_rec = gb_model.recommend(profile)
# Combine predictions
```

### For Production
```python
# Load model once at startup
_MODEL = MLPipelineRecommender()
_MODEL.load_model('models/')

def get_recommendation(profile):
    return _MODEL.recommend(profile)
```

## Troubleshooting

### "Model not found"
```bash
# Train a model first
python example_ml_workflow.py --n-datasets 100
```

### "Import error"
```bash
# Install dependencies
pip install -r requirements_ml.txt
```

### "Low confidence predictions"
- Train on more diverse data
- Use ensemble of models
- Check data quality
- Consider rule-based fallback

## What's Next?

1. **Collect Real Benchmarks**: Run RAPTOR comparisons and save results
2. **Retrain Model**: Incorporate your data for better accuracy
3. **Share Model**: Export for collaborators or publication
4. **Monitor Performance**: Track recommendation accuracy over time

## Need Help?

- **Full Documentation**: See ML_RECOMMENDER_README.md
- **Technical Details**: See IMPLEMENTATION_SUMMARY.md
- **Run Tests**: `python test_ml_system.py`
- **Email**: ayehbolouki1988@gmail.com

## Quick Command Reference

```bash
# Test system
python test_ml_system.py

# Run complete demo
python example_ml_workflow.py

# Train model
python raptor_ml_cli.py train-ml --benchmark-dir data/

# Get recommendation
python raptor_ml_cli.py profile --counts data.csv --use-ml

# Generate training data
python raptor_ml_cli.py generate-data --n-datasets 200
```

## Success Criteria

✓ All tests pass
✓ Model trains in <5 seconds
✓ Test accuracy >80%
✓ Predictions in <0.1 seconds
✓ Confidence scores provided
✓ Visualizations generated

---

**Ready to upgrade RAPTOR? Start with:**
```bash
python example_ml_workflow.py
```

🦖 Happy analyzing!
