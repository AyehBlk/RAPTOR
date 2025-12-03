#  RAPTOR v2.1.0 ML Training Guide

**Train Custom Machine Learning Models for Your Lab**

Learn how to train custom ML models tailored to your specific RNA-seq data characteristics and analysis requirements.

---

##  Table of Contents

1. [Overview](#overview)
2. [When to Train Custom Models](#when-to-train-custom-models)
3. [Prerequisites](#prerequisites)
4. [Quick Start](#quick-start)
5. [Data Collection](#data-collection)
6. [Data Preparation](#data-preparation)
7. [Model Training](#model-training)
8. [Model Evaluation](#model-evaluation)
9. [Model Deployment](#model-deployment)
10. [Advanced Topics](#advanced-topics)
11. [Troubleshooting](#troubleshooting)
12. [Best Practices](#best-practices)

---

##  Overview

RAPTOR v2.1.0 includes a pre-trained ML model trained on **10,000+ real-world RNA-seq analyses**. However, you can train custom models for:

- **Lab-specific workflows** - Optimize for your specific data types
- **Novel organisms** - Better predictions for non-model species
- **Specialized applications** - Single-cell, spatial, long-read RNA-seq
- **Resource constraints** - Optimize for available hardware
- **Publication requirements** - Document custom validation

### How It Works

```
Your Past Analyses
        ↓
   [Feature Extraction]
        ↓
   [Model Training]
        ↓
  Custom ML Model
        ↓
   [Validation]
        ↓
  Deploy in RAPTOR
```

---

##  When to Train Custom Models

### Train Custom Model When:

✅ **Large analysis history** (20+ completed analyses)  
✅ **Consistent data types** (same organism, tissue types)  
✅ **Specialized workflow** (unusual experimental designs)  
✅ **Better than default** (can measure improvement)  
✅ **Publication requirement** (need documented validation)  

### Use Default Model When:

❌ **Few analyses** (<20 projects)  
❌ **Diverse data types** (mixed organisms/tissues)  
❌ **Standard workflows** (default works well)  
❌ **Limited time/resources** (training takes effort)  
❌ **First-time user** (learn with defaults first)  

---

##  Prerequisites

### Required Knowledge

- Basic RNA-seq analysis experience
- Understanding of ML concepts (helpful but not required)
- Python basics (for advanced customization)

### Software Requirements

```bash
# Install RAPTOR with ML tools
pip install raptor-rnaseq[ml]

# Or install dependencies separately
pip install scikit-learn>=1.3.0
pip install xgboost>=2.0.0
pip install shap>=0.42.0
```

### Data Requirements

**Minimum:**
- 20 completed RNA-seq analyses
- Known pipeline outcomes
- Documented data characteristics

**Recommended:**
- 50+ analyses
- Multiple pipeline comparisons
- Ground truth (if available)

---

##  Quick Start

### 1. Export Your Analysis History

```bash
# Export all past RAPTOR analyses
raptor ml export \
  --results-dir /path/to/all_raptor_results/ \
  --output training_data.csv
```

### 2. Train Model

```bash
# Train custom model
raptor ml train \
  --data training_data.csv \
  --output my_custom_model.pkl \
  --validate
```

### 3. Use Custom Model

```bash
# Get recommendations with custom model
raptor ml recommend \
  --counts data/counts.csv \
  --model my_custom_model.pkl
```

**That's it!** 🎉 Your custom model is ready to use.

---

##  Data Collection

### Automatic Collection

RAPTOR automatically tracks analysis metadata:

```bash
# Check what data is available
raptor ml check-data

# Expected output:
Found 47 completed analyses
├─ With ground truth: 23
├─ With benchmarks: 15
├─ With user feedback: 31
└─ Complete metadata: 47

Sufficient data for training ✅
```

### Manual Collection

If you have analyses from other tools:

```bash
# Create training data template
raptor ml create-template --output training_template.csv

# Opens template:
# analysis_id,organism,tissue,n_samples,bcv,depth,...
```

**Fill in template with your data:**

```csv
analysis_id,organism,tissue,n_samples,bcv,sequencing_depth,library_prep,read_length,best_pipeline,accuracy_f1
project_001,Homo sapiens,brain,12,0.42,25000000,polyA,150,3,0.88
project_002,Mus musculus,liver,8,0.38,30000000,totalRNA,100,1,0.92
project_003,Homo sapiens,blood,24,0.65,20000000,polyA,150,5,0.85
...
```

### Collecting Ground Truth

```python
# Extract ground truth from your analyses
from raptor.ml import extract_ground_truth

# From DESeq2 results
ground_truth = extract_ground_truth(
    results_file='deseq2_results.csv',
    fdr_threshold=0.05,
    log2fc_threshold=1.0
)

# From validated genes
ground_truth = extract_ground_truth(
    validated_genes='qpcr_validated.csv',
    match_to='raptor_results.csv'
)
```

---

##  Data Preparation

### 1. Load and Validate Data

```python
from raptor.ml import TrainingDataValidator

# Load data
validator = TrainingDataValidator('training_data.csv')

# Check data quality
report = validator.validate()

print(report)
```

**Example output:**
```
Training Data Validation Report
════════════════════════════════

✅ Data Quality: PASS

Records: 47
Features: 25
Target variable: best_pipeline

Missing Values:
├─ organism: 0 (0%)
├─ bcv: 2 (4%) - Will impute
├─ sequencing_depth: 1 (2%) - Will impute
└─ All other features: Complete

Class Balance:
├─ Pipeline 1: 8 (17%)
├─ Pipeline 3: 18 (38%)
├─ Pipeline 4: 12 (26%)
├─ Pipeline 5: 6 (13%)
└─ Others: 3 (6%)

⚠️  Warning: Limited data for pipelines 2, 6, 7, 8
   Recommendation: Collect more diverse examples

Data Quality Score: 85/100 ✅
Ready for training!
```

### 2. Feature Engineering

```python
from raptor.ml import FeatureEngineer

# Automatic feature engineering
engineer = FeatureEngineer()

# Add derived features
engineer.add_feature('depth_per_sample', 
                     lambda df: df['sequencing_depth'] / df['n_samples'])
engineer.add_feature('cv_adjusted_bcv',
                     lambda df: df['bcv'] * (1 + df['library_size_cv']))

# One-hot encode categoricals
engineer.encode_categorical(['organism', 'tissue', 'library_prep'])

# Transform data
X_train, y_train = engineer.fit_transform(training_data)
```

### 3. Handle Imbalanced Classes

```python
from raptor.ml import balance_classes

# Balance training data
X_balanced, y_balanced = balance_classes(
    X_train, 
    y_train,
    method='smote',  # or 'oversample', 'undersample', 'class_weight'
    min_samples_per_class=10
)
```

---

##  Model Training

### Basic Training

```bash
# Train with defaults
raptor ml train \
  --data training_data.csv \
  --output models/my_model.pkl \
  --algorithm random_forest
```

### Advanced Training

```python
from raptor.ml import ModelTrainer

# Initialize trainer
trainer = ModelTrainer(
    algorithm='random_forest',
    n_estimators=200,
    max_depth=15,
    min_samples_split=5,
    random_state=42
)

# Train model
trainer.fit(X_train, y_train)

# Save model
trainer.save('models/my_custom_model.pkl')
```

### Hyperparameter Optimization

```python
from raptor.ml import hyperparameter_search

# Define search space
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 15, 20, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

# Grid search
best_model, best_params = hyperparameter_search(
    X_train, y_train,
    param_grid=param_grid,
    cv=5,
    scoring='accuracy',
    n_jobs=-1
)

print(f"Best parameters: {best_params}")
print(f"Best CV score: {best_model.best_score_:.3f}")
```

### Bayesian Optimization

```python
from raptor.ml import bayesian_optimize

# Bayesian optimization (faster than grid search)
best_model = bayesian_optimize(
    X_train, y_train,
    n_iterations=50,
    random_state=42
)

# Much faster for large search spaces!
```

### Training with Cross-Validation

```python
from sklearn.model_selection import cross_val_score

# 5-fold cross-validation
cv_scores = cross_val_score(
    trainer.model,
    X_train, y_train,
    cv=5,
    scoring='accuracy'
)

print(f"CV Accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
```

---

##  Model Evaluation

### Comprehensive Evaluation

```bash
# Evaluate on test set
raptor ml evaluate \
  --model models/my_model.pkl \
  --test-data test_data.csv \
  --output evaluation_report.html
```

### Evaluation Metrics

```python
from raptor.ml import evaluate_model

# Load model
model = load_model('models/my_model.pkl')

# Evaluate
metrics = evaluate_model(
    model,
    X_test, y_test,
    output_file='evaluation_report.html'
)

print(metrics)
```

**Output:**
```
Model Evaluation Report
═══════════════════════

Overall Accuracy: 87.2%

Per-Pipeline Performance:
┌──────────┬───────────┬───────────┬──────────┐
│ Pipeline │ Precision │ Recall    │ F1-Score │
├──────────┼───────────┼───────────┼──────────┤
│ 1        │ 0.90      │ 0.85      │ 0.87     │
│ 3        │ 0.92      │ 0.90      │ 0.91     │
│ 4        │ 0.85      │ 0.88      │ 0.86     │
│ 5        │ 0.80      │ 0.75      │ 0.77     │
│ Others   │ 0.70      │ 0.65      │ 0.67     │
└──────────┴───────────┴───────────┴──────────┘

Confusion Matrix:
         Predicted
       1   3   4   5  Other
A  1 [17   2   1   0   0 ]
c  3 [ 1  27   2   0   0 ]
t  4 [ 0   1  21   2   0 ]
u  5 [ 1   0   2  15   2 ]
a Other[ 0   0   1   1   3 ]

✅ Model performs well on common pipelines
⚠️  Limited performance on rare pipelines

Recommendation: Deploy model ✅
```

### Feature Importance

```python
from raptor.ml import explain_model
import matplotlib.pyplot as plt

# Get feature importance
importance = explain_model(model, feature_names=X.columns)

# Plot
importance.plot(kind='barh', figsize=(10, 8))
plt.title('Feature Importance')
plt.xlabel('Importance Score')
plt.tight_layout()
plt.savefig('feature_importance.png')
```

**Example output:**
```
Feature Importance (Top 10)

bcv                    ████████████████ 0.18
sequencing_depth       █████████████    0.15
n_samples              ██████████       0.12
library_size_cv        ████████         0.10
zero_inflation         ███████          0.08
organism_human         ██████           0.07
read_length            █████            0.06
tissue_brain           ████             0.05
depth_per_sample       ████             0.05
library_prep_polyA     ███              0.04
```

### Model Explainability (SHAP)

```python
from raptor.ml import explain_prediction
import shap

# SHAP explainer
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# Visualize for specific prediction
explain_prediction(
    model,
    X_test.iloc[0],
    feature_names=X.columns,
    output_file='explanation.png'
)
```

---

##  Model Deployment

### Deploy Custom Model

```bash
# Set as default model
raptor ml deploy \
  --model models/my_model.pkl \
  --set-default

# Verify deployment
raptor ml list-models
```

**Output:**
```
Available Models:
├─ default (RAPTOR v2.1.0)
│  └─ Trained on 10,000+ analyses
│
└─ my_custom_model ⭐ DEFAULT
   ├─ Trained: 2025-11-19
   ├─ Accuracy: 87.2%
   ├─ Training samples: 47
   └─ Last used: 5 minutes ago
```

### Use Custom Model

```bash
# Get recommendations
raptor profile \
  --counts data/counts.csv \
  --model my_custom_model

# Or specify model path
raptor profile \
  --counts data/counts.csv \
  --model /path/to/my_model.pkl
```

### Comparison with Default

```bash
# Compare custom vs default model
raptor ml compare \
  --model1 default \
  --model2 my_custom_model \
  --test-data test_data.csv
```

**Output:**
```
Model Comparison
════════════════

Test Set: 15 samples

Metric           Default    Custom     Winner
─────────────────────────────────────────────
Accuracy         0.85       0.87       Custom
Precision        0.84       0.86       Custom
Recall           0.86       0.88       Custom
F1-Score         0.85       0.87       Custom

Top-1 Match      80%        87%        Custom
Top-3 Match      93%        93%        Tie

Recommendation: Use custom model ✅
```

---

##  Advanced Topics

### Ensemble Models

```python
from raptor.ml import EnsembleTrainer

# Train ensemble of models
ensemble = EnsembleTrainer(
    models=[
        ('rf', RandomForestClassifier()),
        ('xgb', XGBClassifier()),
        ('lgbm', LGBMClassifier())
    ],
    voting='soft'  # or 'hard'
)

# Train ensemble
ensemble.fit(X_train, y_train)

# Ensemble often performs better!
```

### Multi-Output Models

```python
# Predict multiple outcomes simultaneously
from raptor.ml import MultiOutputTrainer

# Train to predict: best_pipeline, runtime, memory
trainer = MultiOutputTrainer(
    targets=['best_pipeline', 'runtime_minutes', 'memory_gb']
)

trainer.fit(X_train, y_train_multi)

# Get comprehensive predictions
predictions = trainer.predict(X_new)
```

### Online Learning

```python
# Update model with new data
from raptor.ml import update_model

# Load existing model
model = load_model('models/my_model.pkl')

# Update with new analyses
model_updated = update_model(
    model,
    X_new, y_new,
    learning_rate=0.1
)

# Save updated model
model_updated.save('models/my_model_v2.pkl')
```

### Transfer Learning

```python
# Start from default model, fine-tune on your data
from raptor.ml import transfer_learn

# Load pre-trained default model
base_model = load_model('default')

# Fine-tune on your data
custom_model = transfer_learn(
    base_model,
    X_train, y_train,
    freeze_layers=5,  # Keep first 5 layers frozen
    epochs=50
)

# Often works better with limited data!
```

### Active Learning

```python
# Intelligently select which analyses to add
from raptor.ml import active_learning_suggest

# Get suggestions for which analyses would improve model most
suggestions = active_learning_suggest(
    model,
    candidate_analyses=unlabeled_data,
    n_suggestions=10,
    strategy='uncertainty'  # or 'diversity', 'expected_improvement'
)

print("Run these analyses next for maximum model improvement:")
for analysis in suggestions:
    print(f"  - {analysis}")
```

---

##  Troubleshooting

### Issue: Poor Model Performance

**Symptoms:**
```
Model Accuracy: 62%
Warning: Model performs poorly
```

**Possible Causes & Solutions:**

```python
# 1. Insufficient data
print(f"Training samples: {len(X_train)}")
# Solution: Collect more analyses (need 50+ ideally)

# 2. Imbalanced classes
print(y_train.value_counts())
# Solution: Balance classes or adjust class weights

# 3. Feature quality
# Check for high correlation, missing values
X_train.corr().abs().max()
# Solution: Feature engineering, remove correlated features

# 4. Overfitting
# Compare train vs test accuracy
print(f"Train: {train_acc:.3f}, Test: {test_acc:.3f}")
# Solution: Regularization, simpler model, more data
```

### Issue: Model Overconfident

**Symptoms:**
```
Confidence: 98% 
Actual: Wrong prediction
```

**Solution - Calibration:**
```python
from sklearn.calibration import CalibratedClassifierCV

# Calibrate model
calibrated_model = CalibratedClassifierCV(
    model,
    method='isotonic',  # or 'sigmoid'
    cv=5
)

calibrated_model.fit(X_train, y_train)

# Now probabilities are more reliable!
```

### Issue: Slow Training

**Solution:**
```python
# 1. Reduce data size (if very large)
from sklearn.model_selection import train_test_split
X_sample, _, y_sample, _ = train_test_split(
    X_train, y_train, 
    train_size=0.5,
    stratify=y_train
)

# 2. Use simpler model
model = RandomForestClassifier(
    n_estimators=100,  # Instead of 300
    max_depth=10       # Instead of None
)

# 3. Parallelize
model = RandomForestClassifier(
    n_jobs=-1  # Use all CPU cores
)

# 4. Use GPU (if available)
import xgboost as xgb
model = xgb.XGBClassifier(
    tree_method='gpu_hist'  # GPU acceleration
)
```

### Issue: Model Won't Load

**Error:**
```
VersionError: Model trained with scikit-learn 1.2, 
current version is 1.3
```

**Solution:**
```bash
# Option 1: Use same version
pip install scikit-learn==1.2.0

# Option 2: Retrain model
raptor ml train --data training_data.csv --output new_model.pkl

# Option 3: Convert model
raptor ml convert --model old_model.pkl --output new_model.pkl
```

---

##  Best Practices

### Data Collection

✅ **Collect diverse examples** - Various data types  
✅ **Document everything** - All metadata  
✅ **Include failures** - Unsuccessful pipelines too  
✅ **Regular updates** - Add new analyses  
✅ **Quality over quantity** - Good labels matter  

### Training

✅ **Start simple** - Basic model first  
✅ **Cross-validate** - Always!  
✅ **Track experiments** - MLflow, Weights & Biases  
✅ **Version models** - my_model_v1, v2, etc.  
✅ **Document hyperparameters** - For reproducibility  

### Evaluation

✅ **Hold-out test set** - Never touch during training  
✅ **Multiple metrics** - Not just accuracy  
✅ **Per-class performance** - Some pipelines harder  
✅ **Explain predictions** - Use SHAP, feature importance  
✅ **Compare to baseline** - Better than default?  

### Deployment

✅ **A/B testing** - Compare with default model  
✅ **Monitor performance** - Track predictions  
✅ **Collect feedback** - Were predictions good?  
✅ **Regular retraining** - With new data  
✅ **Version control** - Git for models too  

---

##  Example Workflow

### Complete Training Pipeline

```bash
#!/bin/bash
# complete_training.sh

# 1. Export historical data
echo "Exporting training data..."
raptor ml export \
  --results-dir /data/all_raptor_results/ \
  --output training_data.csv

# 2. Validate data
echo "Validating data..."
raptor ml validate --data training_data.csv

# 3. Train model
echo "Training model..."
raptor ml train \
  --data training_data.csv \
  --algorithm random_forest \
  --optimize-hyperparameters \
  --cv 5 \
  --output models/custom_model.pkl

# 4. Evaluate
echo "Evaluating model..."
raptor ml evaluate \
  --model models/custom_model.pkl \
  --test-data test_data.csv \
  --output evaluation_report.html

# 5. Compare with default
echo "Comparing with default..."
raptor ml compare \
  --model1 default \
  --model2 models/custom_model.pkl \
  --test-data test_data.csv

# 6. Deploy if better
echo "Deploying model..."
raptor ml deploy \
  --model models/custom_model.pkl \
  --set-default

echo "✅ Training complete!"
```

---

##  Learning Resources

### Tutorials

- [Tutorial: Basic Model Training](tutorials/ml_training_basic.md)
- [Tutorial: Advanced Techniques](tutorials/ml_training_advanced.md)
- [Tutorial: Production Deployment](tutorials/ml_production.md)

### Examples

```bash
# Example 1: Train for specific organism
examples/ml_training/organism_specific.py

# Example 2: Train for single-cell data
examples/ml_training/single_cell.py

# Example 3: Multi-output prediction
examples/ml_training/multi_output.py
```

### External Resources

- [Scikit-learn Documentation](https://scikit-learn.org/)
- [XGBoost Guide](https://xgboost.readthedocs.io/)
- [SHAP Documentation](https://shap.readthedocs.io/)

---

##  Support

**Need help with ML training?**

1. Check [FAQ.md](FAQ.md) - ML section
2. Read [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
3. GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
4. Email: ayehbolouki1988@gmail.com

---

##  Summary

Custom ML model training allows you to:
- ✅ **Optimize for your lab** - Better predictions
- ✅ **Handle specialized data** - Novel organisms, methods
- ✅ **Improve over time** - Learn from your analyses
- ✅ **Document rigorously** - For publications
- ✅ **Share with community** - Help others

**With 50+ quality training examples, custom models typically achieve 85-90% accuracy!** 🚀

---

**Author:** Ayeh Bolouki  
**Version:** 2.1.0  
**License:** MIT

---

*"Train once, predict forever - make RAPTOR truly yours!"*
