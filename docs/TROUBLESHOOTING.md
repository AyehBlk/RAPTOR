# 🔧 RAPTOR v2.1.1 Troubleshooting Guide

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

A comprehensive guide to diagnosing and resolving common issues with RAPTOR v2.1.1.

---

## 📋 Table of Contents

1. [Installation Issues](#installation-issues)
2. [Configuration Problems](#configuration-problems)
3. [ML Recommendation Issues](#ml-recommendation-issues)
4. [🎯 Threshold Optimizer Issues](#threshold-optimizer-issues) - **NEW!**
5. [Dashboard Problems](#dashboard-problems)
6. [Pipeline Execution Errors](#pipeline-execution-errors)
7. [Resource Monitoring Issues](#resource-monitoring-issues)
8. [Ensemble Analysis Problems](#ensemble-analysis-problems)
9. [Data Quality Issues](#data-quality-issues)
10. [Cloud Deployment Issues](#cloud-deployment-issues)
11. [Performance Problems](#performance-problems)
12. [Getting Help](#getting-help)

---

## 🆕 What's New in v2.1.1 Troubleshooting

- Added [Threshold Optimizer Issues](#threshold-optimizer-issues) section
- Updated installation verification for ATO
- Added ATO-specific dashboard troubleshooting

---

## 📦 Installation Issues

### Issue: pip install fails

**Symptoms:**
```bash
ERROR: Could not find a version that satisfies the requirement raptor
```

**Solutions:**

1. **Update pip:**
```bash
pip install --upgrade pip
```

2. **Check Python version:**
```bash
python --version  # Should be 3.8+
```

3. **Install from GitHub:**
```bash
pip install git+https://github.com/AyehBlk/RAPTOR.git@v2.1.1
```

4. **Use virtual environment:**
```bash
python -m venv raptor_env
source raptor_env/bin/activate  # On Windows: raptor_env\Scripts\activate
pip install raptor-rnaseq
```

---

### Issue: Dependency conflicts

**Symptoms:**
```bash
ERROR: pip's dependency resolver does not currently take into account all packages...
```

**Solutions:**

1. **Create fresh environment:**
```bash
conda create -n raptor python=3.9
conda activate raptor
pip install raptor-rnaseq
```

2. **Install specific versions:**
```bash
pip install scikit-learn==1.3.0 pandas==2.0.0
pip install raptor-rnaseq
```

---

### Issue: Threshold Optimizer not available (NEW)

**Symptoms:**
```bash
ModuleNotFoundError: No module named 'raptor.threshold_optimizer'
```

**Solutions:**

1. **Update RAPTOR:**
```bash
pip install --upgrade raptor-rnaseq
```

2. **Verify installation:**
```bash
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

3. **Check version:**
```bash
python -c "import raptor; print(raptor.__version__)"
# Should show 2.1.1
```

---

## ⚙️ Configuration Problems

### Issue: Config file not found

**Symptoms:**
```bash
FileNotFoundError: config.yaml not found
```

**Solutions:**

1. **Create default config:**
```bash
raptor config --create
```

2. **Specify config path:**
```bash
raptor analyze --config /path/to/config.yaml
```

---

### Issue: Invalid threshold_optimizer config (NEW)

**Symptoms:**
```bash
ValueError: Invalid threshold_optimizer configuration
```

**Solutions:**

1. **Check config structure:**
```yaml
threshold_optimizer:
  enabled: true
  goal: "discovery"  # Must be: discovery, balanced, or validation
  default_logfc_method: "auto"
```

2. **Validate goal value:**
```python
# Valid goals
valid_goals = ["discovery", "balanced", "validation"]
```

3. **Use config template:**
```bash
raptor config --template > config.yaml
```

---

## 🤖 ML Recommendation Issues

### Issue: Model not found

**Symptoms:**
```bash
FileNotFoundError: ML model file not found
```

**Solutions:**

1. **Train model first:**
```bash
raptor ml train --data training_data.csv
```

2. **Download pre-trained model:**
```bash
raptor ml download --model default
```

---

### Issue: Low confidence predictions

**Symptoms:**
```
Warning: Low confidence (45%) - manual review recommended
```

**Solutions:**

1. **Provide more metadata**
2. **Use ensemble mode**
3. **Consider using Threshold Optimizer for data-driven thresholds** (NEW in v2.1.1)

---

## 🎯 Threshold Optimizer Issues (NEW)

### Issue: ATO module not found

**Symptoms:**
```bash
ImportError: cannot import name 'AdaptiveThresholdOptimizer'
```

**Solutions:**

1. **Update to v2.1.1:**
```bash
pip install --upgrade raptor-rnaseq
```

2. **Verify version:**
```bash
python -c "import raptor; print(raptor.__version__)"
# Must be 2.1.1 or higher
```

3. **Check import:**
```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer, optimize_thresholds
```

---

### Issue: Column not found in data

**Symptoms:**
```bash
KeyError: 'log2FoldChange' not found in columns
```

**Solutions:**

1. **Check column names:**
```python
import pandas as pd
df = pd.read_csv('results.csv')
print(df.columns.tolist())
```

2. **Specify correct column names:**
```python
ato = AdaptiveThresholdOptimizer(
    df,
    logfc_col='logFC',      # edgeR uses 'logFC'
    pvalue_col='PValue'     # edgeR uses 'PValue'
)
```

3. **Common column names by tool:**
```python
# DESeq2
logfc_col='log2FoldChange', pvalue_col='pvalue'

# edgeR
logfc_col='logFC', pvalue_col='PValue'

# limma
logfc_col='logFC', pvalue_col='P.Value'
```

---

### Issue: π₀ estimation fails

**Symptoms:**
```bash
Warning: π₀ estimation failed, using default value 0.9
```

**Solutions:**

1. **Check p-value distribution:**
```python
import matplotlib.pyplot as plt
plt.hist(df['pvalue'].dropna(), bins=50)
plt.show()
# Should show uniform distribution with spike at 0
```

2. **Ensure sufficient data:**
```python
# Need at least 100 genes with valid p-values
valid_pvalues = df['pvalue'].dropna()
print(f"Valid p-values: {len(valid_pvalues)}")  # Should be >100
```

3. **Check for NA values:**
```python
print(f"NA p-values: {df['pvalue'].isna().sum()}")
# Remove or impute NA values
```

4. **Try different π₀ method:**
```python
result = ato.optimize(goal='discovery', pi0_method='histogram')
```

---

### Issue: LogFC threshold is too extreme

**Symptoms:**
```
Optimal logFC threshold: 5.2  # Seems too high
Significant genes: 3  # Very few genes
```

**Solutions:**

1. **Check data distribution:**
```python
print(f"LogFC range: {df['log2FoldChange'].min():.2f} to {df['log2FoldChange'].max():.2f}")
print(f"LogFC std: {df['log2FoldChange'].std():.2f}")
```

2. **Try different method:**
```python
# Use percentile method for extreme distributions
result = ato.optimize(goal='discovery', logfc_method='percentile')
```

3. **Use manual threshold:**
```python
# If auto methods fail, set reasonable default
result = ato.optimize(goal='discovery')
if result.logfc_threshold > 3:
    print("Warning: Using fallback threshold")
    # Consider using result.logfc_threshold = 1.0
```

---

### Issue: No significant genes found

**Symptoms:**
```
Significant genes: 0
```

**Solutions:**

1. **Check thresholds:**
```python
print(f"LogFC threshold: {result.logfc_threshold}")
print(f"P-value threshold: {result.pvalue_threshold}")
```

2. **Try more permissive goal:**
```python
# Change from validation to discovery
result = ato.optimize(goal='discovery')  # More permissive
```

3. **Check data quality:**
```python
# Ensure p-values are not all 1.0 or NA
print(df['pvalue'].describe())
```

4. **Verify data is differential expression results:**
```python
# Should have logFC and pvalue columns with proper values
print(df[['log2FoldChange', 'pvalue']].head())
```

---

### Issue: Methods text not generated

**Symptoms:**
```python
print(result.methods_text)
# Output: None or empty string
```

**Solutions:**

1. **Check result object:**
```python
print(type(result))  # Should be ThresholdResult
print(result)  # View all fields
```

2. **Access directly:**
```python
from raptor.threshold_optimizer import AdaptiveThresholdOptimizer
ato = AdaptiveThresholdOptimizer(df, 'log2FoldChange', 'pvalue')
result = ato.optimize(goal='balanced')
print(result.methods_text)  # Should have text
```

3. **Generate manually:**
```python
methods = f"""Thresholds determined using ATO v2.1.1 with '{result.goal}' goal.
LogFC threshold: {result.logfc_threshold:.3f}
P-value threshold: {result.pvalue_threshold}
Significant genes: {result.n_significant}"""
```

---

### Issue: Dashboard ATO page not working

**Symptoms:**
```
Dashboard loads but Threshold Optimizer page shows error
```

**Solutions:**

1. **Check ATO availability:**
```bash
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('OK')"
```

2. **Update RAPTOR:**
```bash
pip install --upgrade raptor-rnaseq
```

3. **Clear browser cache:**
- Chrome: Ctrl+Shift+Del
- Firefox: Ctrl+Shift+Del

4. **Restart dashboard:**
```bash
# Kill existing process
lsof -ti:8501 | xargs kill -9

# Restart
raptor dashboard
```

---

## 📊 Dashboard Problems

### Issue: Dashboard won't start

**Symptoms:**
```bash
OSError: [Errno 98] Address already in use
```

**Solutions:**

1. **Change port:**
```bash
raptor dashboard --port 8502
```

2. **Kill existing process:**
```bash
lsof -ti:8501 | xargs kill -9
```

---

### Issue: Threshold Optimizer page missing (NEW)

**Symptoms:**
```
No "🎯 Threshold Optimizer" option in sidebar
```

**Solutions:**

1. **Verify v2.1.1:**
```bash
python -c "import raptor; print(raptor.__version__)"
```

2. **Update RAPTOR:**
```bash
pip install --upgrade raptor-rnaseq
```

3. **Check dashboard version:**
Look for "v2.1.1" in dashboard header

---

### Issue: Dashboard shows no data

**Symptoms:**
```
No results found. Run analysis first.
```

**Solutions:**

1. **Run analysis first:**
```bash
raptor analyze --config config.yaml
```

2. **Specify results directory:**
```bash
raptor dashboard --results /path/to/results
```

---

## 🔬 Pipeline Execution Errors

### Issue: Pipeline fails immediately

**Symptoms:**
```bash
Error: Pipeline 'salmon' failed at initialization
```

**Solutions:**

1. **Check tool installation:**
```bash
salmon --version
star --version
```

2. **Add tools to PATH:**
```bash
export PATH=$PATH:/path/to/salmon/bin
```

---

### Issue: Out of memory

**Symptoms:**
```bash
java.lang.OutOfMemoryError: Java heap space
```

**Solutions:**

1. **Increase memory limits:**
```yaml
resources:
  max_memory: "32GB"
```

2. **Reduce parallelism:**
```yaml
resources:
  threads: 8  # Reduce from 16
```

---

## 📈 Resource Monitoring Issues

### Issue: Monitoring not starting

**Symptoms:**
```bash
Warning: Resource monitoring disabled
```

**Solutions:**

1. **Install psutil:**
```bash
pip install psutil>=5.9.0
```

2. **Enable in config:**
```yaml
resource_monitoring:
  enabled: true
```

---

## 🎭 Ensemble Analysis Problems

### Issue: Ensemble fails to combine

**Symptoms:**
```bash
Error: Cannot combine results from different pipelines
```

**Solutions:**

1. **Ensure consistent settings:**
```yaml
ensemble:
  normalize_method: "TMM"
```

2. **Use ATO for uniform thresholds (NEW):**
```yaml
threshold_optimizer:
  enabled: true
  ensemble_mode:
    uniform_thresholds: true
```

---

## 📊 Data Quality Issues

### Issue: Low quality scores

**Symptoms:**
```
Warning: Quality score 45/100 - data may be problematic
```

**Solutions:**

1. **Run detailed QC:**
```bash
raptor qc --fastq data/ --detailed
```

2. **Consider ATO for appropriate thresholds:**
```python
# ATO can help select thresholds appropriate for your data quality
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(de_results, goal='discovery')
```

---

## ☁️ Cloud Deployment Issues

### Issue: Authentication failed

**Symptoms:**
```bash
Error: Unable to authenticate with AWS
```

**Solutions:**

**AWS:**
```bash
aws configure
```

**GCP:**
```bash
gcloud auth login
```

---

## 🚀 Performance Problems

### Issue: Analysis very slow

**Symptoms:**
```
Expected: 2 hours, Actual: 12 hours
```

**Solutions:**

1. **Increase parallelism:**
```yaml
resources:
  threads: 32
```

2. **Use faster pipelines:**
```yaml
pipelines:
  - "salmon"  # Faster than STAR
```

3. **ATO is fast (NEW):**
```python
# ATO optimization takes <1 second
result = optimize_thresholds(df, goal='balanced')
```

---

## 🆘 Getting Help

### Before Asking for Help

1. **Check logs:**
```bash
cat raptor.log
cat raptor_error.log
```

2. **Enable debug mode:**
```yaml
logging:
  level: "DEBUG"
```

3. **Collect system info:**
```bash
raptor --version
python --version
pip list | grep raptor
```

4. **Verify ATO (v2.1.1):**
```bash
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('ATO: OK')"
```

### Reporting Issues

**Include in your report:**

1. **RAPTOR version:**
```bash
raptor --version
```

2. **Configuration:**
```bash
cat config.yaml
```

3. **Error message:**
```bash
tail -n 50 raptor_error.log
```

4. **For ATO issues, include:**
```python
print(df.columns.tolist())
print(df.shape)
print(df[['logfc_col', 'pvalue_col']].describe())
```

### Support Channels

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Email:** ayehbolouki1988@gmail.com
- **Documentation:** https://github.com/AyehBlk/RAPTOR/wiki

---

## 📚 Additional Resources

- [Installation Guide](INSTALLATION.md)
- [Configuration Guide](CONFIGURATION.md)
- [Threshold Optimizer Guide](THRESHOLD_OPTIMIZER.md) - **NEW!**
- [API Documentation](API.md)
- [FAQ](FAQ.md)

---

**Author:** Ayeh Bolouki   
**Version:** 2.1.1  
**License:** MIT

---

*"Every bug is a feature waiting to be understood."* 🐛
