# RAPTOR v2.1.0 Configuration Examples

This directory contains example configuration files for different RAPTOR usage scenarios. Each configuration is optimized for specific use cases and workflows.

---

## 📁Available Configurations

### 1. **config.yaml** - Complete Reference Configuration
**Purpose**: Full configuration file with all possible settings documented

**When to use**:
- Learning about all RAPTOR features
- Reference for available options
- Customizing your own configuration

**Features**: All v2.1.0 features enabled with detailed documentation

**Computational Requirements**:
- CPU: 8 cores
- RAM: 32 GB
- Time: 2-4 hours (typical dataset)

**Command**:
```bash
raptor profile --counts data.csv --config config.yaml
```

---

### 2. **config_minimal.yaml** - Quick Start Configuration
**Purpose**: Minimal configuration for basic RNA-seq analysis

**When to use**:
- First time using RAPTOR
- Quick exploratory analysis
- Limited computational resources
- Standard RNA-seq workflow

**Features**:
- ✅ Basic profiling
- ✅ ML recommendations
- ✅ Quality assessment
- ✅ Resource monitoring
- ❌ Ensemble analysis
- ❌ Parameter optimization
- ❌ Cloud integration

**Computational Requirements**:
- CPU: 4 cores
- RAM: 16 GB
- Time: 1-2 hours

**Perfect for**:
- Students and beginners
- Laptop analysis
- Quick pipeline testing
- Proof-of-concept studies

**Command**:
```bash
raptor profile --counts data.csv --config config_minimal.yaml
```

---

### 3. **config_ml_advanced.yaml** - Machine Learning Research
**Purpose**: Advanced ML configuration for model training and hyperparameter tuning

**When to use**:
- Training custom ML models
- Hyperparameter optimization experiments
- Model comparison studies
- ML research and development
- Building production models

**Features**:
- ✅ Multiple ML models (RF, XGBoost, Neural Networks)
- ✅ Hyperparameter tuning (Grid, Random, Bayesian)
- ✅ Ensemble ML methods
- ✅ SHAP explanations
- ✅ Comprehensive model evaluation
- ✅ Experiment tracking (MLflow ready)

**Computational Requirements**:
- CPU: 12 cores (16+ recommended)
- RAM: 32 GB (64 GB for large datasets)
- GPU: Optional but recommended
- Time: Variable (4-24 hours for training)

**Perfect for**:
- Data scientists
- ML researchers
- Bioinformaticians developing new methods
- Labs with computational clusters

**Command**:
```bash
# Train models
raptor ml-train --benchmarks ./benchmark_data/ --config config_ml_advanced.yaml

# Use trained model
raptor profile --counts data.csv --config config_ml_advanced.yaml --use-ml
```

**Notes**:
- Requires diverse training data (500+ datasets recommended)
- GPU acceleration available for Neural Networks and XGBoost
- Supports MLflow for experiment tracking

---

### 4. **config_cloud_aws.yaml** - AWS Cloud Deployment
**Purpose**: Configuration for running RAPTOR on AWS infrastructure

**When to use**:
- Large-scale batch processing
- Institutional cloud deployment
- Cost-effective spot instances
- Auto-scaling workloads
- Remote data storage (S3)

**Features**:
- ✅ AWS Batch integration
- ✅ S3 data storage
- ✅ EC2 spot instances (60-90% cost savings)
- ✅ CloudWatch monitoring
- ✅ Auto-scaling
- ✅ Container deployment

**Computational Requirements**:
- Instance: c5.4xlarge (16 vCPU, 32 GB RAM)
- Storage: S3 bucket
- Cost: ~$0.20-0.40/hour (spot instances)

**Perfect for**:
- Large cohort studies (100+ samples)
- Multi-project analyses
- Labs with AWS accounts
- Cost-sensitive large-scale analyses

**Prerequisites**:
```bash
# AWS CLI configured
aws configure

# S3 bucket created
aws s3 mb s3://your-raptor-bucket

# Docker image available
docker pull ayehblk/raptor:2.1.0
```

**Command**:
```bash
raptor profile \
  --counts s3://your-raptor-bucket/data.csv \
  --config config_cloud_aws.yaml \
  --cloud
```

**Cost Estimation**:
- Typical RNA-seq analysis: $2-10 per sample
- 100 samples: ~$200-1000
- Storage: ~$0.023/GB/month

**Notes**:
- Remember to update bucket names in config
- Configure IAM roles appropriately
- Monitor costs with CloudWatch alarms

---

### 5. **config_publication_ready.yaml** - Publication Standards
**Purpose**: Configuration for manuscript preparation and publication

**When to use**:
- Preparing manuscripts
- Generating supplementary materials
- Peer review responses
- Grant applications
- Conference presentations

**Features**:
- ✅ Comprehensive quality control
- ✅ Biological interpretation (GO, KEGG, Reactome, etc.)
- ✅ Publication-quality figures (300 DPI)
- ✅ Detailed methods section generation
- ✅ Reproducibility documentation
- ✅ Statistical rigor
- ✅ Supplementary materials

**Computational Requirements**:
- CPU: 8 cores
- RAM: 32 GB
- Time: 3-6 hours

**Perfect for**:
- Academic researchers
- Graduate students
- Postdocs preparing papers
- PI labs

**Output Includes**:
- High-resolution figures (PNG, PDF, SVG)
- Comprehensive HTML/PDF reports
- Supplementary tables (all DE genes, enrichment results)
- Methods section text
- Software version documentation
- Analysis provenance tracking

**Command**:
```bash
raptor profile --counts data.csv --config config_publication_ready.yaml
```

**Publication Checklist**:
- ☑️ Random seed set (reproducibility)
- ☑️ QC metrics documented
- ☑️ Pipeline selection justified
- ☑️ Statistical methods described
- ☑️ Multiple testing correction applied
- ☑️ Effect sizes reported
- ☑️ 300 DPI figures
- ☑️ Software versions documented
- ☑️ Code available
- ☑️ Data archived (GEO/SRA)

**Citation**:
Includes automated citation generation for RAPTOR and dependencies.

---

### 6. **config_ensemble.yaml** - Ensemble Analysis
**Purpose**: Combine multiple pipelines for robust differential expression

**When to use**:
- High-stakes clinical decisions
- Biomarker discovery
- Discordant results between pipelines
- Maximum confidence required
- Meta-analysis across studies

**Features**:
- ✅ Multiple pipeline execution (3-6 pipelines)
- ✅ Statistical meta-analysis
- ✅ Rank aggregation
- ✅ Weighted voting
- ✅ Concordance assessment
- ✅ Confidence scoring
- ✅ Discordance analysis

**Computational Requirements**:
- CPU: 12-16 cores
- RAM: 48-64 GB
- Time: 4-8 hours (runs multiple pipelines)
- Disk: 50-100 GB (intermediate files)

**Perfect for**:
- Clinical/diagnostic applications
- Drug target validation
- Conflicting literature findings
- High-impact publications
- Regulatory submissions

**Ensemble Methods**:
1. **Rank Aggregation**: Combine gene rankings
2. **Vote Counting**: Majority vote across pipelines
3. **Weighted Average**: Weight by pipeline performance
4. **Meta-Analysis**: Statistical combination (Fisher's, Stouffer's)
5. **Intersection**: Most conservative (all agree)
6. **Union**: Most liberal (any calls)

**Confidence Levels**:
- **Very High (>90%)**: Clinical decisions
- **High (75-90%)**: Robust biomarkers
- **Medium (60-75%)**: Validation recommended
- **Low (<60%)**: Investigate discordance

**Command**:
```bash
raptor ensemble --counts data.csv --config config_ensemble.yaml
```

**Output Includes**:
- Venn diagrams (pipeline overlap)
- UpSet plots (multi-way intersections)
- Concordance heatmaps
- Confidence scores per gene
- Discordance analysis report

---

##  Quick Selection Guide

**Choose based on your priority:**

| Priority | Configuration | Use Case |
|----------|--------------|----------|
| **Speed** | `config_minimal.yaml` | Quick results, learning |
| **Robustness** | `config_ensemble.yaml` | Maximum confidence |
| **Publication** | `config_publication_ready.yaml` | Manuscript preparation |
| **ML Research** | `config_ml_advanced.yaml` | Model development |
| **Scale** | `config_cloud_aws.yaml` | Large datasets, cloud |
| **Learning** | `config.yaml` | See all options |

---

##  Customization Tips

### Starting from Minimal
```yaml
# Start with minimal.yaml
# Then enable features as needed:

ml_recommendation:
  enabled: true  # Add ML recommendations

ensemble:
  enabled: true  # Add ensemble analysis
  min_pipelines: 3
```

### Combining Configs
```bash
# Use cloud infrastructure with publication quality
raptor profile \
  --counts s3://bucket/data.csv \
  --config config_cloud_aws.yaml \
  --output-format publication
```

### Override Specific Settings
```bash
# Use minimal config but change thresholds
raptor profile \
  --counts data.csv \
  --config config_minimal.yaml \
  --fdr-threshold 0.01 \
  --lfc-threshold 1.5
```

---

##  Best Practices

### 1. **Always Start Simple**
Begin with `config_minimal.yaml` to understand your data, then move to advanced configs.

### 2. **Document Your Choices**
```yaml
# In your custom config, add comments:
fdr_threshold: 0.01  # More stringent due to large cohort
```

### 3. **Version Control Your Configs**
```bash
git add config_custom.yaml
git commit -m "Analysis config for Project X"
```

### 4. **Test Before Production**
```bash
# Test with small dataset first
raptor profile --counts test_data.csv --config config_ensemble.yaml
```

### 5. **Save Configs with Results**
```bash
# RAPTOR automatically saves the config used
# Find it in: raptor_output/analysis_parameters.yaml
```

---

##  Comparison Table

| Feature | Minimal | ML Advanced | Cloud AWS | Publication | Ensemble |
|---------|---------|-------------|-----------|-------------|----------|
| **ML Recommendation** | ✅ Basic | ✅ Advanced | ✅ Basic | ✅ Basic | ✅ Advanced |
| **Quality Assessment** | ✅ Basic | ✅ Advanced | ✅ Basic | ✅ Comprehensive | ✅ Basic |
| **Resource Monitoring** | ✅ | ✅ | ✅ CloudWatch | ✅ | ✅ |
| **Ensemble Analysis** | ❌ | ❌ | ❌ | Optional | ✅ |
| **Parameter Optimization** | ❌ | ✅ | ❌ | ❌ | ❌ |
| **Cloud Integration** | ❌ | ❌ | ✅ AWS | ❌ | ❌ |
| **Biological Interpretation** | ❌ | ❌ | ✅ Basic | ✅ Comprehensive | ✅ Basic |
| **Publication Figures** | Basic | Advanced | Basic | 300 DPI | Advanced |
| **Time (typical)** | 1-2h | 4-24h | 2-4h | 3-6h | 4-8h |
| **RAM Required** | 16 GB | 32-64 GB | 32 GB | 32 GB | 48-64 GB |
| **Complexity** | ⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ |

---

##  Getting Help

### Configuration Issues
```bash
# Validate your config
raptor validate-config --config your_config.yaml

# Get default config
raptor generate-config --output default_config.yaml
```

### Common Questions

**Q: Can I combine multiple configs?**
A: Yes! Use `--config` multiple times. Later configs override earlier ones:
```bash
raptor profile --config base.yaml --config custom.yaml --counts data.csv
```

**Q: How do I know which pipelines are selected?**
A: Check the log file or report. RAPTOR documents all decisions.

**Q: What if I run out of memory?**
A: Use `config_minimal.yaml` or reduce `default_threads` in your config.

**Q: Can I use ensemble + cloud?**
A: Absolutely! Combine settings from both configs.

---

##  Additional Resources

- **RAPTOR Documentation**: See `docs/` folder
- **Tutorials**: `docs/tutorials/`
- **GitHub Issues**: https://github.com/AyehBlk/RAPTOR/issues
- **Email**: ayehbolouki1988@gmail.com

---

##  Version History

- **v2.1.0** (2025-11): Initial example configs with 8 major new features
- **v2.0.0** (2023): Base configuration system

---

**Author**: Ayeh Bolouki  
**License**: MIT  

*Making free science for everybody around the world* 🌍
