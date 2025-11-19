# RAPTOR v2.1.0 Machine Learning Guide

**Complete guide to ML-based pipeline recommendations**

---

## 🎯 What's New in v2.1.0

RAPTOR v2.1.0 introduces **AI-powered pipeline recommendations** using machine learning models trained on extensive benchmark data. Instead of rule-based heuristics, ML models learn from hundreds of real analyses to predict which pipeline will perform best for your specific data.

---

## 🤖 How It Works

### Traditional Approach (v2.0.0)
```
Your Data → Rule-based Scoring → Recommendation
```
- Simple rules (if BCV > 0.4, prefer edgeR)
- Limited by human-defined thresholds
- Doesn't learn from new data

### ML Approach (v2.1.0)
```
Your Data → Feature Extraction → ML Model → Confidence Score → Recommendation
```
- Learns patterns from 500+ benchmark datasets
- Adapts to complex interactions
- Provides confidence scores
- Explains predictions

---

## 🚀 Quick Start

### Basic ML Recommendation

```bash
# Enable ML recommendations (default in v2.1.0)
raptor profile --counts data.csv --use-ml

# View ML model information
raptor ml-info

# Train custom model (advanced)
raptor ml-train --benchmarks ./training_data/
```

### Python API

```python
from raptor import MLPipelineRecommender, RNAseqDataProfiler

# Profile your data
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# Get ML recommendation
ml_recommender = MLPipelineRecommender(model_type='random_forest')
recommendation = ml_recommender.recommend(profile)

print(f"Recommended: Pipeline {recommendation['pipeline_id']}")
print(f"Confidence: {recommendation['confidence']:.1%}")
print(f"Reasoning: {recommendation['reasoning']}")
```

---

## 📊 Understanding ML Recommendations

### Recommendation Output

```python
{
    'pipeline_id': 3,
    'pipeline_name': 'Salmon-edgeR',
    'confidence': 0.87,  # 87% confidence
    'confidence_level': 'high',  # high, medium, low
    
    # Why this pipeline?
    'reasoning': [
        'Excellent balance for medium BCV (0.42)',
        'Optimal for your sample size (n=12)',
        'Fast runtime expected (~1 hour)',
        'Low memory requirements (8-12 GB)'
    ],
    
    # Feature importance
    'key_features': {
        'bcv': 0.35,
        'sample_size': 0.25,
        'sequencing_depth': 0.20,
        'zero_inflation': 0.12,
        'library_size_cv': 0.08
    },
    
    # Alternative suggestions
    'alternatives': [
        {
            'pipeline_id': 1,
            'pipeline_name': 'STAR-RSEM-DESeq2',
            'confidence': 0.78,
            'reason': 'Higher accuracy, but slower'
        },
        {
            'pipeline_id': 4,
            'pipeline_name': 'Kallisto-Sleuth',
            'confidence': 0.71,
            'reason': 'Faster, but slightly lower accuracy'
        }
    ],
    
    # Expected performance
    'expected_performance': {
        'accuracy': 0.88,
        'runtime_hours': 1.2,
        'memory_gb': 10,
        'sensitivity': 0.89,
        'precision': 0.87
    }
}
```

### Confidence Levels

| Confidence | Meaning | Action |
|------------|---------|--------|
| **High (>80%)** | Model is very confident | Safe to use this pipeline |
| **Medium (60-80%)** | Reasonable confidence | Consider alternatives |
| **Low (<60%)** | Model uncertain | Run benchmarking or try multiple |

---

## 🔍 Features Used by ML Models

### Core Features (Always Used)

1. **Biological Coefficient of Variation (BCV)**
   - Most important feature (~35% importance)
   - Measures biological variability
   - Range: 0.1 (low) to 1.0+ (high)

2. **Sample Size**
   - Number of samples per condition
   - Critical for statistical power
   - Small (n<6), Medium (6-12), Large (>12)

3. **Sequencing Depth**
   - Total reads per sample
   - Affects gene detection
   - Low (<10M), Medium (10-30M), High (>30M)

4. **Zero Inflation**
   - Percentage of zero counts
   - Indicates data sparsity
   - Typical: 30-60%

5. **Library Size Coefficient of Variation**
   - Variability in total counts
   - Good: <0.2, Moderate: 0.2-0.4, High: >0.4

### Advanced Features (Optional)

6. **Dispersion Pattern**
   - How variance relates to mean
   - Affects method choice

7. **Batch Effect Score**
   - Strength of batch effects
   - Important for multi-batch studies

8. **Gene Detection Rate**
   - Proportion of expressed genes
   - Quality indicator

9. **Dynamic Range**
   - Span of expression levels
   - Affects normalization

10. **Outlier Count**
    - Number of problematic samples
    - Influences robust method selection

---

## 🎓 ML Models Available

### 1. Random Forest (Default)

**Best for:** General use, robust predictions

```python
recommender = MLPipelineRecommender(model_type='random_forest')
recommender.load_model('models/')
recommendation = recommender.recommend(profile)
```

**Characteristics:**
- ✅ Fast predictions (<1 second)
- ✅ Handles non-linear relationships
- ✅ Provides feature importance
- ✅ Robust to outliers
- ⚠️ Can overfit small datasets

**Parameters:**
```yaml
ml_recommendation:
  random_forest:
    n_estimators: 200
    max_depth: 15
    min_samples_split: 5
    min_samples_leaf: 2
```

### 2. Gradient Boosting

**Best for:** Highest accuracy when well-tuned

```python
recommender = MLPipelineRecommender(model_type='gradient_boosting')
```

**Characteristics:**
- ✅ Often highest accuracy
- ✅ Captures complex patterns
- ✅ Feature importance
- ⚠️ Slower training
- ⚠️ Requires careful tuning

**Parameters:**
```yaml
ml_recommendation:
  gradient_boosting:
    n_estimators: 150
    learning_rate: 0.1
    max_depth: 10
```

### 3. Neural Network

**Best for:** Very large training datasets, complex patterns

```python
recommender = MLPipelineRecommender(model_type='neural_network')
```

**Characteristics:**
- ✅ Captures very complex relationships
- ✅ Scales to large datasets
- ⚠️ Requires more training data (500+)
- ⚠️ Less interpretable

**Parameters:**
```yaml
ml_recommendation:
  neural_network:
    hidden_layers: [128, 64, 32]
    activation: 'relu'
    dropout_rate: 0.3
    epochs: 100
```

### 4. Ensemble (All Models Combined)

**Best for:** Maximum reliability, production use

```python
recommender = MLPipelineRecommender(model_type='ensemble')
```

**Characteristics:**
- ✅ Combines predictions from multiple models
- ✅ Most robust
- ✅ Best confidence estimates
- ⚠️ Slower predictions
- ⚠️ Requires all models trained

---

## 🧪 Training Custom ML Models

### When to Train Custom Models

✅ **Train custom models if:**
- You have unique data characteristics
- Working with non-standard organisms
- Have extensive benchmark data (100+)
- Need specialized optimization

❌ **Use default models if:**
- Standard RNA-seq experiments
- Limited training data
- First time using RAPTOR

### Training Workflow

#### Step 1: Generate Training Data

```bash
# Option A: Use existing benchmarks
raptor ml-train \
  --benchmarks ./my_benchmarks/ \
  --output models/custom/

# Option B: Generate synthetic data
raptor simulate-training-data \
  --n-datasets 500 \
  --output training_data/ \
  --diversity high

raptor ml-train \
  --training-data training_data/ \
  --output models/custom/
```

#### Step 2: Train Models

```python
from raptor import MLPipelineRecommender, BenchmarkDataLoader

# Load training data
loader = BenchmarkDataLoader('./training_data/')
X_train, y_train = loader.load_features_and_labels()

# Train model
recommender = MLPipelineRecommender(model_type='random_forest')
results = recommender.train(X_train, y_train, cv_folds=5)

# View training results
print(f"Cross-validation accuracy: {results['cv_accuracy']:.2%}")
print(f"Best parameters: {results['best_params']}")

# Save model
recommender.save_model('models/my_custom_model/')
```

#### Step 3: Evaluate Model

```python
# Load test data
X_test, y_test = loader.load_test_data()

# Evaluate
metrics = recommender.evaluate(X_test, y_test)

print(f"Test Accuracy: {metrics['accuracy']:.2%}")
print(f"F1 Score: {metrics['f1_score']:.2%}")
print(f"Confusion Matrix:\n{metrics['confusion_matrix']}")
```

#### Step 4: Use Custom Model

```bash
# Use your custom model
raptor profile \
  --counts data.csv \
  --use-ml \
  --model-path models/my_custom_model/
```

---

## 📈 Model Performance Metrics

### Understanding Model Quality

```python
from raptor import MLPipelineRecommender

recommender = MLPipelineRecommender()
recommender.load_model('models/')

# Get model performance
performance = recommender.get_performance_metrics()

print(performance)
```

**Output:**
```python
{
    'training_samples': 523,
    'test_samples': 131,
    'accuracy': 0.87,
    'balanced_accuracy': 0.85,
    'f1_score': 0.86,
    'precision': 0.88,
    'recall': 0.84,
    
    # Per-pipeline performance
    'per_pipeline': {
        1: {'precision': 0.92, 'recall': 0.88},
        2: {'precision': 0.78, 'recall': 0.75},
        3: {'precision': 0.91, 'recall': 0.89},
        4: {'precision': 0.85, 'recall': 0.82},
        5: {'precision': 0.88, 'recall': 0.86},
        6: {'precision': 0.80, 'recall': 0.77},
        7: {'precision': 0.82, 'recall': 0.79},
        8: {'precision': 0.75, 'recall': 0.71}
    },
    
    # Training details
    'training_time': 45.3,  # seconds
    'model_size_mb': 12.4,
    'features_used': 15,
    'cv_folds': 5
}
```

### Benchmarking Against Ground Truth

```bash
# Compare ML recommendations with known best pipeline
raptor ml-benchmark \
  --test-data benchmarks/test_set/ \
  --output ml_evaluation/
```

**Results:**
```
ML Model Evaluation Results
════════════════════════════════════

Overall Performance:
  Accuracy: 87.2%
  Top-3 Accuracy: 96.5%  # Correct pipeline in top 3
  
Pipeline-Specific:
  Pipeline 1: 92% accuracy (46/50 correct)
  Pipeline 3: 91% accuracy (89/98 correct)
  Pipeline 4: 85% accuracy (42/49 correct)
  ...
  
Common Mistakes:
  • Confused P1 ↔ P3 in 8 cases (similar performance)
  • Underestimated P5 in 3 cases (complex designs)
  
Confidence Calibration:
  High confidence: 94% correct
  Medium confidence: 79% correct
  Low confidence: 61% correct
```

---

## 🔧 Advanced Configuration

### Fine-tuning ML Behavior

Create `config/ml_config.yaml`:

```yaml
ml_recommendation:
  # Model selection
  enabled: true
  model_type: "random_forest"
  model_path: "./models/raptor_rf_model.pkl"
  
  # Feature engineering
  features:
    - "bcv"
    - "sequencing_depth"
    - "num_samples"
    - "num_replicates"
    - "dispersion"
    - "library_size_cv"
    - "zero_inflation"
    - "batch_effect_score"
  
  # Feature transformations
  feature_engineering:
    polynomial_features: true
    polynomial_degree: 2
    log_transform: ["sequencing_depth", "num_samples"]
    feature_selection:
      enabled: true
      method: "mutual_info"
      n_features: 15
  
  # Prediction settings
  prediction:
    confidence_threshold: 0.7
    top_n_recommendations: 3
    explain_predictions: true
    
    # Uncertainty estimation
    uncertainty_estimation:
      enabled: true
      method: "dropout_monte_carlo"
      n_iterations: 100
  
  # Fallback behavior
  fallback:
    use_rule_based: true  # If ML confidence < threshold
    min_confidence: 0.5
```

### Hyperparameter Tuning

```python
from raptor import MLPipelineRecommender, HyperparameterOptimizer

# Define search space
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 15, 20],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

# Optimize
optimizer = HyperparameterOptimizer(
    model_type='random_forest',
    param_grid=param_grid,
    cv_folds=5
)

best_params = optimizer.optimize(X_train, y_train)
print(f"Best parameters: {best_params}")

# Train with best parameters
recommender = MLPipelineRecommender(
    model_type='random_forest',
    **best_params
)
recommender.train(X_train, y_train)
```

---

## 🎯 Use Cases & Examples

### Example 1: Standard RNA-seq Study

**Your data:**
- 12 samples (6 vs 6)
- BCV: 0.38 (moderate)
- Depth: 25M reads/sample
- No batch effects

**ML Recommendation:**
```
Pipeline 3 (Salmon-edgeR)
Confidence: 89% (High)

Reasoning:
✓ Optimal for moderate BCV
✓ Good sample size (n=12)
✓ Balanced speed/accuracy
✓ Low resource requirements

Expected: F1=0.88, Runtime=1.2h, Memory=10GB
```

### Example 2: Small Pilot Study

**Your data:**
- 6 samples (3 vs 3)
- BCV: 0.52 (high)
- Depth: 15M reads/sample

**ML Recommendation:**
```
Pipeline 1 (STAR-RSEM-DESeq2)
Confidence: 72% (Medium)

Reasoning:
✓ Best for small samples
✓ Handles high variation
⚠ Will take longer
⚠ Requires 32GB RAM

Alternative: Pipeline 3 (75% as good, much faster)
```

### Example 3: Large Clinical Cohort

**Your data:**
- 48 samples (24 vs 24)
- BCV: 0.68 (very high)
- Batch effects present

**ML Recommendation:**
```
Pipeline 5 (STAR-HTSeq-limma-voom)
Confidence: 81% (High)

Reasoning:
✓ Excellent for large samples
✓ Built-in batch correction
✓ Handles high variation
✓ Flexible modeling

Expected: F1=0.86, Runtime=4h, Memory=40GB
```

### Example 4: Uncertain Case

**Your data:**
- Unusual characteristics
- Mixed quality samples

**ML Recommendation:**
```
Pipeline 3 (Salmon-edgeR)
Confidence: 58% (Low)

⚠️ Low confidence recommendation

Suggested actions:
1. Run ensemble analysis (combine multiple pipelines)
2. Perform full benchmarking
3. Review data quality carefully

Alternatives:
- Pipeline 1: 56% confidence
- Pipeline 4: 54% confidence
```

---

## 🐛 Troubleshooting ML Features

### Issue: Model not found

**Error:**
```
FileNotFoundError: Model file not found at models/raptor_rf_model.pkl
```

**Solution:**
```bash
# Option 1: Download pre-trained model
raptor ml-download-model

# Option 2: Train your own
raptor ml-train --training-data ./data/ --quick

# Option 3: Use rule-based (fallback)
raptor profile --counts data.csv --no-ml
```

### Issue: Low confidence predictions

**Possible causes:**
- Data outside training distribution
- Unusual characteristics
- Model needs retraining

**Solutions:**
```bash
# Check data characteristics
raptor profile --counts data.csv --detailed

# Use ensemble for higher confidence
raptor profile --counts data.csv --use-ml --ensemble

# Train model with similar data
raptor ml-train --training-data ./similar_benchmarks/
```

### Issue: Wrong recommendations

**If ML consistently gives poor recommendations:**

```python
# Validate model performance
from raptor import MLPipelineRecommender

recommender = MLPipelineRecommender()
recommender.load_model('models/')

# Check feature importance
importance = recommender.get_feature_importance()
print(importance)

# If features seem wrong, retrain
recommender.train(X_new, y_new)
recommender.save_model('models/updated/')
```

---

## 📚 Further Reading

### Papers & Methods

- **Random Forests:** Breiman, L. (2001). Machine Learning, 45(1), 5-32.
- **Feature Importance:** Lundberg & Lee (2017). SHAP values for model interpretation
- **RNA-seq Benchmarking:** Soneson & Delorenzi (2013). Genome Biology, 14:R95

### RAPTOR Resources

- **Benchmarking Guide:** [BENCHMARKING.md](BENCHMARKING.md)
- **API Reference:** [API.md](API.md)
- **Configuration:** [config/ml_config.yaml](../config/examples/config_ml_advanced.yaml)

---

## 🎓 Best Practices

### ✅ Do:
- Trust high-confidence recommendations (>80%)
- Review reasoning provided
- Consider alternatives for medium confidence
- Use ensemble for critical analyses
- Document ML model version used

### ❌ Don't:
- Ignore low confidence warnings
- Use outdated models for new data types
- Skip data quality checks
- Assume ML is always right
- Forget to cite RAPTOR and ML methods

---

## 📊 ML Model Changelog

### v2.1.0 (Current)
- Initial ML recommendation system
- Random Forest (default)
- Gradient Boosting option
- Neural Network support
- Ensemble mode
- Trained on 523 benchmark datasets

### Future (v2.2.0)
- Transfer learning for new organisms
- Online learning (model updates)
- Multi-task learning
- Confidence intervals
- Explainable AI enhancements

---

## 📧 Support

Questions about ML features?

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Discussions:** https://github.com/AyehBlk/RAPTOR/discussions
- **Email:** ayehbolouki1988@gmail.com

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**Version:** 2.1.0  
**License:** MIT  

---

*"From rule-based heuristics to machine learning - making pipeline selection smarter."* 🤖
