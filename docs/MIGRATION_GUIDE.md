# RAPTOR Migration Guide: v2.0.0 → v2.1.0

**Complete guide for upgrading to RAPTOR v2.1.0**

---

## 🎯 Overview

RAPTOR v2.1.0 is a **major feature release** that adds 8 new capabilities while maintaining **100% backward compatibility** with v2.0.0.

### What's New in v2.1.0

1. 🤖 **ML-Based Pipeline Recommendations** - AI-powered selection
2. 📊 **Advanced Quality Assessment** - Comprehensive QC metrics
3. ⚡ **Real-Time Resource Monitoring** - Track CPU, memory, GPU
4. 🎯 **Ensemble Analysis** - Combine multiple pipelines
5. 🌐 **Interactive Dashboard** - Web-based interface
6. 🔧 **Parameter Optimization** - Automated tuning
7. 📄 **Automated Reporting** - Enhanced reports with interpretation
8. ☁️ **Cloud Integration** - AWS/GCP/Azure support

---

## ✅ Is Migration Right for You?

### Migrate to v2.1.0 if:
- ✅ You want ML-powered recommendations
- ✅ You need ensemble analysis for robustness
- ✅ You want the interactive dashboard
- ✅ You're starting new analyses
- ✅ You have time to explore new features

### Stay on v2.0.0 if:
- ⚠️ You're in the middle of a publication (wait until submitted)
- ⚠️ Your current workflow is working perfectly
- ⚠️ You have strict computational constraints
- ⚠️ You need to reproduce exact v2.0.0 results

**Good news:** v2.1.0 is **100% backward compatible**, so all v2.0.0 commands still work!

---

## 📦 Installation

### Option 1: Upgrade Existing Installation

```bash
# If installed via pip
pip install --upgrade raptor-rnaseq

# If installed via conda
conda update -c bioconda raptor

# Verify version
raptor --version  # Should show v2.1.0
```

### Option 2: Fresh Installation (Recommended)

```bash
# Create new conda environment
conda create -n raptor-v2.1 python=3.9
conda activate raptor-v2.1

# Install RAPTOR v2.1.0 with all features
conda install -c bioconda raptor

# Or from GitHub
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
git checkout v2.1.0
pip install -e .

# Verify
raptor --version
```

### Option 3: Keep Both Versions (Safest)

```bash
# Keep v2.0.0 environment
conda activate raptor-v2.0  # Your existing environment

# Create separate v2.1.0 environment
conda create -n raptor-v2.1 python=3.9
conda activate raptor-v2.1
pip install raptor-rnaseq==2.1.0

# Switch between versions as needed
conda activate raptor-v2.0  # For old analyses
conda activate raptor-v2.1  # For new analyses
```

---

## 🔧 New Dependencies

### Required (Core Features)

```bash
# These are automatically installed
pip install scikit-learn>=1.3.0  # For ML features
pip install psutil>=5.9.0        # For resource monitoring
```

### Optional (Advanced Features)

```bash
# For interactive dashboard
pip install streamlit>=1.28.0
pip install plotly>=5.17.0

# For cloud integration
pip install boto3>=1.28.0      # AWS
pip install google-cloud>=2.0  # GCP
pip install azure-storage>=12.0 # Azure

# For parameter optimization
pip install optuna>=3.3.0      # Bayesian optimization

# Install all optional dependencies
pip install raptor-rnaseq[all]
```

---

## 📝 Configuration Files

### v2.0.0 Config Files Still Work!

Your existing `config.yaml` files are **100% compatible** with v2.1.0.

```bash
# Your old config works as-is
raptor profile \
  --counts data.csv \
  --config my_v2.0_config.yaml  # ✓ Still works!
```

### New v2.1.0 Config Options

To use new features, add these sections to your config:

```yaml
# Add to your existing config.yaml

# ML Recommendations (NEW)
ml_recommendation:
  enabled: true
  model_type: "random_forest"
  confidence_threshold: 0.7

# Quality Assessment (ENHANCED)
quality_assessment:
  enabled: true
  metrics:
    - "library_size"
    - "detected_genes"
    - "mitochondrial_content"
    - "mapping_rate"

# Resource Monitoring (NEW)
resource_monitoring:
  enabled: true
  sampling_interval: 1.0

# Ensemble Analysis (NEW)
ensemble:
  enabled: false  # Enable when needed
  methods:
    - "meta_analysis"
    - "vote_counting"

# Dashboard (NEW)
dashboard:
  enabled: false  # Launch separately
  port: 8501

# Automated Reporting (ENHANCED)
automated_reporting:
  enabled: true
  interpretation:
    enabled: true
    databases: ["GO", "KEGG", "Reactome"]
```

**Or use example configs:**

```bash
# Copy example config for your use case
cp config/examples/config_ml_advanced.yaml my_config.yaml
cp config/examples/config_ensemble.yaml my_ensemble_config.yaml
```

---

## 🔄 Command Changes

### Good News: No Breaking Changes!

All v2.0.0 commands work identically in v2.1.0:

```bash
# These all work exactly the same
raptor profile --counts data.csv
raptor compare --data fastq/ --output results/
raptor run --pipeline 3 --data fastq/
```

### New Commands in v2.1.0

```bash
# NEW: ML-powered recommendations
raptor profile --counts data.csv --use-ml

# NEW: Ensemble analysis
raptor ensemble --counts data.csv --pipelines 1,3,5

# NEW: Launch dashboard
raptor dashboard

# NEW: Train custom ML model
raptor ml-train --benchmarks ./data/

# NEW: Monitor resources
raptor monitor --pid 12345

# NEW: Parameter optimization
raptor optimize --counts data.csv --pipeline 3
```

---

## 🎯 Migration Scenarios

### Scenario 1: Basic User (Just Profiling)

**v2.0.0 workflow:**
```bash
raptor profile --counts data.csv --output report.html
```

**v2.1.0 equivalent (exactly the same):**
```bash
raptor profile --counts data.csv --output report.html
```

**v2.1.0 enhanced (optional):**
```bash
# Enable ML for smarter recommendations
raptor profile --counts data.csv --use-ml --output report.html
```

**Migration effort:** None (or 5 minutes to try ML)

### Scenario 2: Power User (Running Pipelines)

**v2.0.0 workflow:**
```bash
# Run profiling
raptor profile --counts data.csv

# Run recommended pipeline
raptor run --pipeline 3 --data fastq/ --output results/

# Generate report
raptor report --results results/ --output final_report.html
```

**v2.1.0 workflow (same commands work):**
```bash
# Everything works identically
raptor profile --counts data.csv
raptor run --pipeline 3 --data fastq/ --output results/
raptor report --results results/ --output final_report.html
```

**v2.1.0 enhanced (with new features):**
```bash
# Use ML recommendations
raptor profile --counts data.csv --use-ml

# Run with resource monitoring
raptor run --pipeline 3 --data fastq/ --monitor

# Generate enhanced report with interpretation
raptor report --results results/ --interpret --output report.html
```

**Migration effort:** 30 minutes to explore new features

### Scenario 3: Researcher (Benchmarking)

**v2.0.0 workflow:**
```bash
# Run benchmark
raptor compare --data fastq/ --pipelines 1,3,5 --output benchmark/

# Analyze results
raptor report --results benchmark/ --output comparison.html
```

**v2.1.0 workflow (works identically):**
```bash
raptor compare --data fastq/ --pipelines 1,3,5 --output benchmark/
raptor report --results benchmark/ --output comparison.html
```

**v2.1.0 enhanced (with ensemble):**
```bash
# Run benchmark (same as before)
raptor compare --data fastq/ --pipelines 1,3,5 --output benchmark/

# NEW: Combine results with ensemble
raptor ensemble --results benchmark/ --method meta_analysis

# Generate comprehensive report
raptor report --results benchmark/ --ensemble --output report.html
```

**Migration effort:** 1 hour to learn ensemble analysis

### Scenario 4: Developer (Python API)

**v2.0.0 code:**
```python
from raptor import RNAseqDataProfiler, PipelineRecommender

# Profile data
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# Get recommendation
recommender = PipelineRecommender()
recommendation = recommender.recommend(profile)
```

**v2.1.0 code (backward compatible):**
```python
# Your v2.0.0 code works unchanged!
from raptor import RNAseqDataProfiler, PipelineRecommender

profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

recommender = PipelineRecommender()
recommendation = recommender.recommend(profile)
```

**v2.1.0 enhanced (with ML):**
```python
# Use new ML recommender
from raptor import RNAseqDataProfiler, MLPipelineRecommender

profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.profile()

# NEW: ML-based recommendation
ml_recommender = MLPipelineRecommender(model_type='random_forest')
recommendation = ml_recommender.recommend(profile)
print(f"Confidence: {recommendation['confidence']:.1%}")
```

**Migration effort:** 2 hours to update scripts with new features

---

## 📊 Results Compatibility

### Will My Results Change?

**Short answer:** No, if you use the same commands.

**Detailed answer:**

| Feature | v2.0.0 Result | v2.1.0 Result | Notes |
|---------|---------------|---------------|-------|
| **Pipeline execution** | Identical | Identical | Same pipelines, same parameters |
| **Rule-based recommendations** | Identical | Identical | Same algorithm |
| **Benchmark comparisons** | Identical | Identical | Same metrics |
| **ML recommendations** | N/A | New feature | Different from rule-based |

### Reproducing v2.0.0 Results

If you need exact v2.0.0 results:

```bash
# Disable new features
raptor profile \
  --counts data.csv \
  --no-ml \
  --legacy-mode \
  --output results/

# Or use v2.0.0 config
raptor profile \
  --counts data.csv \
  --config my_v2.0_config.yaml
```

---

## 🚨 Potential Issues & Solutions

### Issue 1: Import Errors

**Error:**
```python
ImportError: cannot import name 'MLPipelineRecommender'
```

**Cause:** Old RAPTOR version cached

**Solution:**
```bash
# Clear Python cache
pip cache purge
pip uninstall raptor-rnaseq
pip install raptor-rnaseq==2.1.0

# Or force reinstall
pip install --force-reinstall raptor-rnaseq==2.1.0
```

### Issue 2: Config Validation Warnings

**Warning:**
```
Warning: Unknown config parameter 'ml_recommendation'
```

**Cause:** Using v2.1.0 config with v2.0.0 installation

**Solution:**
```bash
# Verify version
raptor --version

# If showing v2.0.0, upgrade
pip install --upgrade raptor-rnaseq

# Or ignore warning (v2.0.0 will skip unknown parameters)
```

### Issue 3: Dashboard Won't Start

**Error:**
```bash
raptor dashboard
# ModuleNotFoundError: No module named 'streamlit'
```

**Solution:**
```bash
# Install dashboard dependencies
pip install streamlit plotly

# Or install all optional dependencies
pip install raptor-rnaseq[all]
```

### Issue 4: ML Model Not Found

**Error:**
```
FileNotFoundError: Model file not found at models/raptor_rf_model.pkl
```

**Solution:**
```bash
# Download pre-trained model
raptor ml-download-model

# Or disable ML temporarily
raptor profile --counts data.csv --no-ml

# Or train your own
raptor ml-train --benchmarks ./data/ --quick
```

### Issue 5: Memory Issues with New Features

**Problem:** Resource monitoring uses extra memory

**Solution:**
```yaml
# Adjust in config
resource_monitoring:
  enabled: true
  sampling_interval: 5.0  # Less frequent sampling
  track_metrics:
    - "cpu_percent"
    - "memory_percent"
    # Remove less critical metrics
```

---

## 📚 Updated Documentation

### New Documentation Files

Read these to learn new features:

1. **[ML_GUIDE.md](ML_GUIDE.md)** - ML recommendations
2. **[ENSEMBLE_GUIDE.md](ENSEMBLE_GUIDE.md)** - Ensemble analysis
3. **[DASHBOARD_GUIDE.md](dashboard/README.md)** - Interactive dashboard
4. **[QUALITY_GUIDE.md](QUALITY_GUIDE.md)** - Quality assessment
5. **[RESOURCE_MONITOR_GUIDE.md](RESOURCE_MONITOR_GUIDE.md)** - Monitoring
6. **[PARAMETER_OPTIMIZATION.md](PARAMETER_OPTIMIZATION.md)** - Tuning
7. **[CLOUD_DEPLOYMENT.md](CLOUD_DEPLOYMENT.md)** - Cloud integration

### Updated Documentation

These have been updated for v2.1.0:

- **INSTALLATION.md** - New dependencies
- **API.md** - New classes and methods
- **FAQ.md** - v2.1.0 questions
- **TROUBLESHOOTING.md** - New feature issues

---

## 🎓 Learning Path

### Week 1: Basic Migration
- Day 1: Install v2.1.0
- Day 2: Test existing workflows
- Day 3: Try ML recommendations
- Day 4: Explore dashboard
- Day 5: Read new documentation

### Week 2: Advanced Features
- Day 1: Learn ensemble analysis
- Day 2: Try resource monitoring
- Day 3: Explore quality assessment
- Day 4: Test parameter optimization
- Day 5: Practice with real data

### Month 1: Full Adoption
- Week 1-2: Basic migration (above)
- Week 3: Integrate into workflows
- Week 4: Train team members

---

## ✅ Migration Checklist

Before migrating to v2.1.0:

### Preparation
- [ ] Read this migration guide
- [ ] Back up existing results
- [ ] Document current workflow
- [ ] Test on non-critical data first

### Installation
- [ ] Install v2.1.0
- [ ] Verify version (`raptor --version`)
- [ ] Test basic commands
- [ ] Check dependencies

### Validation
- [ ] Run simple test analysis
- [ ] Compare results with v2.0.0 (optional)
- [ ] Test new features you'll use
- [ ] Update scripts if using API

### Adoption
- [ ] Update documentation
- [ ] Train team members
- [ ] Update workflows
- [ ] Monitor for issues

---

## 🆘 Getting Help

### Resources

1. **Documentation:** Read [ML_GUIDE.md](ML_GUIDE.md) and [ENSEMBLE_GUIDE.md](ENSEMBLE_GUIDE.md)
2. **FAQ:** Check [FAQ.md](FAQ.md) for common questions
3. **Troubleshooting:** See [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
4. **Examples:** Look at [config/examples/](../config/examples/)

### Community Support

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Discussions:** https://github.com/AyehBlk/RAPTOR/discussions
- **Email:** ayehbolouki1988@gmail.com

### Reporting Migration Issues

If you encounter problems:

```markdown
**Issue Title:** v2.1.0 Migration: [Brief Description]

**Information:**
- Previous version: v2.0.0
- New version: v2.1.0
- Operating system: [e.g., Ubuntu 22.04]
- Installation method: [pip/conda/source]

**Problem:**
[Describe what's not working]

**Expected:**
[What you expected to happen]

**Steps to reproduce:**
1. [First step]
2. [Second step]
3. [...]

**Error messages:**
```
[Paste error messages here]
```

**Already tried:**
- [ ] Reinstalled RAPTOR
- [ ] Checked documentation
- [ ] Searched existing issues
```

---

## 📊 Feature Comparison

### v2.0.0 vs v2.1.0

| Feature | v2.0.0 | v2.1.0 | Notes |
|---------|--------|--------|-------|
| **8 Pipeline implementations** | ✅ | ✅ | Identical |
| **Data profiling** | ✅ | ✅ Enhanced | More metrics |
| **Rule-based recommendations** | ✅ | ✅ | Same algorithm |
| **ML recommendations** | ❌ | ✅ | **NEW** |
| **Benchmarking** | ✅ | ✅ | Same |
| **Ensemble analysis** | ❌ | ✅ | **NEW** |
| **Interactive dashboard** | ❌ | ✅ | **NEW** |
| **Resource monitoring** | Basic | ✅ Advanced | **ENHANCED** |
| **Quality assessment** | Basic | ✅ Comprehensive | **ENHANCED** |
| **Parameter optimization** | Manual | ✅ Automated | **NEW** |
| **Automated reporting** | ✅ | ✅ Enhanced | Interpretation added |
| **Cloud integration** | ❌ | ✅ | **NEW** |


---

## 📅 Version Timeline

- **v2.0.0** (Released: Jan 2025) - Initial release
- **v2.1.0** (Released: Nov 2025) - Major feature update
- **v2.1.1** (Planned: Dec 2025) - Bug fixes
- **v2.2.0** (Planned: Q1 2026) - Additional features

---

## 🎯 Conclusion

### Key Takeaways

✅ **v2.1.0 is backward compatible** - Your v2.0.0 workflows work unchanged

✅ **Migration is optional** - But highly recommended for new projects

✅ **New features are opt-in** - Use what you need, ignore the rest

✅ **Comprehensive documentation** - Guides for every new feature

✅ **Community support** - Help available via GitHub and email

### Recommended Migration Timeline

- **Immediate:** If starting new project
- **Within 1 month:** For active researchers
- **Within 3 months:** For established workflows
- **When convenient:** For published work

---

**Welcome to RAPTOR v2.1.0!** 🦖✨

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & University of Liège, Belgium  
**Version:** 2.1.0  
**Date:** November 2025  
**License:** MIT

---

*"Evolution, not revolution - upgrading RAPTOR while preserving what works."* 🔄
