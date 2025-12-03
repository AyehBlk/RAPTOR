#  RAPTOR v2.1.0 Frequently Asked Questions (FAQ)

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

Your questions answered about RAPTOR v2.1.0!

---

##  Table of Contents

1. [General Questions](#general-questions)
2. [Installation & Setup](#installation--setup)
3. [ML Recommendations](#ml-recommendations)
4. [Dashboard](#dashboard)
5. [Pipeline Selection](#pipeline-selection)
6. [Resource Requirements](#resource-requirements)
7. [Cloud Computing](#cloud-computing)
8. [Results & Output](#results--output)
9. [Troubleshooting](#troubleshooting)
10. [Advanced Usage](#advanced-usage)

---

##  General Questions

### What is RAPTOR?

RAPTOR (RNA-seq Analysis Pipeline Testing and Optimization Resource) is a comprehensive framework for testing, comparing, and optimizing RNA-seq analysis pipelines. Version 2.1.0 adds ML-based recommendations, interactive dashboards, and cloud integration.

**Key Features:**
-  ML-based pipeline recommendations (85-90% accuracy)
-  Interactive web dashboard (no coding required)
-  Real-time resource monitoring
-  Ensemble analysis across pipelines
-  Cloud computing support (AWS/GCP/Azure)

---

### Who should use RAPTOR?

**Perfect for:**
- ✅ Researchers new to RNA-seq
- ✅ Core facilities processing diverse samples
- ✅ Labs wanting to optimize pipelines
- ✅ Anyone comparing multiple analysis methods
- ✅ Users wanting evidence-based tool selection

**Not ideal for:**
- ❌ Those needing only single-cell RNA-seq (scRNA-seq)*
- ❌ Users with only web-based tools (no command line)
- ❌ Projects requiring specialized, custom pipelines

*Single-cell support planned for v2.2.0

---

### What's new in v2.1.0?

**8 Major New Features:**

1. **ML Recommendations** - Smart pipeline selection
2. **Quality Assessment** - Comprehensive data QC
3. **Resource Monitoring** - Real-time CPU/memory tracking
4. **Ensemble Analysis** - Combine multiple pipelines
5. **Interactive Dashboard** - Visual, no-code interface
6. **Parameter Optimization** - Automated tuning
7. **Automated Reporting** - Publication-ready outputs
8. **Cloud Integration** - AWS/GCP/Azure support

[See complete changelog](CHANGELOG.md)

---

### Is RAPTOR free?

**Yes!** RAPTOR is completely free and open-source under the MIT License.

- ✅ Free to use
- ✅ Free to modify
- ✅ Free for commercial use
- ✅ No registration required
- ✅ No data sharing required

However, cloud computing (AWS/GCP/Azure) has costs if you choose to use it.

---

### How do I cite RAPTOR?

```
Bolouki, A. (2025). RAPTOR v2.1.0: RNA-seq Analysis Pipeline 
Testing and Optimization Resource with Machine Learning 
Recommendations. GitHub. https://github.com/AyehBlk/RAPTOR
```

Full publication coming soon!

---

##  Installation & Setup

### What are the system requirements?

**Minimum:**
- Python 3.8+
- 8 GB RAM
- 50 GB free disk space
- Linux/macOS/Windows WSL

**Recommended:**
- Python 3.9+
- 32 GB RAM
- 500 GB free disk space
- 16+ CPU cores

[Detailed requirements](INSTALLATION.md#system-requirements)

---

### How long does installation take?

**Quick install (pip):** 5-10 minutes
```bash
pip install raptor-rnaseq
```

**Full install with tools:** 30-60 minutes
```bash
conda install -c bioconda raptor salmon star kallisto
```

[Installation guide](INSTALLATION.md)

---

### Can I install without conda?

Yes! Multiple installation methods:

```bash
# Option 1: pip only
pip install raptor-rnaseq

# Option 2: Docker
docker pull ayehblk/raptor:2.1.0

# Option 3: From source
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -e .
```

---

### Do I need to install all pipelines?

**No!** Install only what you need:

```bash
# Just Salmon (fast, accurate)
conda install -c bioconda salmon

# Salmon + Kallisto (ensemble mode)
conda install -c bioconda salmon kallisto

# All tools (full testing)
conda install -c bioconda salmon star kallisto rsem hisat2
```

RAPTOR will detect available tools automatically.

---

## 🤖 ML Recommendations

### How accurate are ML recommendations?

**Very accurate!** Based on extensive testing:

- ✅ 85-90% accuracy overall
- ✅ 95%+ for common organisms
- ✅ 80-85% for novel conditions

**Confidence scores:**
- High (>70%): Trust the recommendation
- Medium (50-70%): Review suggestion
- Low (<50%): Manual selection advised

---

### What data does the ML model need?

**Required:**
- Organism (e.g., Homo sapiens)
- Tissue/cell type
- Read length

**Recommended:**
- Sequencing depth
- Library preparation method
- Number of samples
- Biological question

**Example:**
```yaml
metadata:
  organism: "Homo sapiens"
  tissue_type: "brain"
  read_length: 150
  sequencing_depth: "50M"
  library_prep: "polyA"
  samples: 20
  analysis_type: "differential_expression"
```

---

### Can I trust ML over my expertise?

**Use ML as a guide, not gospel!**

**Advantages:**
- ✅ Learned from 10,000+ analyses
- ✅ Considers many factors simultaneously
- ✅ Reduces bias from limited experience

**Limitations:**
- ❌ May not know YOUR specific needs
- ❌ Can't consider unpublished methods
- ❌ Black box (can't always explain why)

**Best practice:** Use ML suggestion as starting point, then apply your expertise!

---

### What if ML confidence is low?

**Don't panic!** Low confidence means:
- Novel experimental design
- Unusual organism/tissue
- Limited training data

**What to do:**
1. **Review metadata** - Add more details
2. **Check alternatives** - ML shows top 3 options
3. **Run ensemble** - Test multiple pipelines
4. **Manual selection** - Use your expertise

```bash
raptor ml recommend --explain
# Shows why confidence is low
```

---

### Can I retrain the model?

**Yes!** Customize for your lab:

```bash
# Create training data
raptor ml export --results raptor_results/

# Add your data
# Edit training_data.csv

# Retrain
raptor ml train --data training_data.csv --output my_model.pkl

# Use custom model
raptor analyze --ml-model my_model.pkl
```

[ML training guide](docs/ML_TRAINING.md)

---

## 📊 Dashboard

### Do I need coding to use the dashboard?

**No coding required!** The dashboard is fully interactive:

- ✅ Point-and-click interface
- ✅ Drag-and-drop for file upload
- ✅ Automatic plot generation
- ✅ Export results with button click

**Start dashboard:**
```bash
raptor dashboard
# Open browser to http://localhost:8501
```

---

### Can multiple people access the dashboard?

**Yes!** Deploy for team access:

**Local network:**
```bash
raptor dashboard --host 0.0.0.0 --port 8501
# Access from other computers: http://YOUR_IP:8501
```

**Public server (with authentication):**
```bash
raptor dashboard --host 0.0.0.0 --auth-required
# Set up password protection
```

[Dashboard deployment guide](docs/DASHBOARD.md#team-deployment)

---

### What if dashboard is slow?

**Common solutions:**

1. **Reduce data displayed:**
```yaml
dashboard:
  max_samples: 50
  lazy_loading: true
```

2. **Use summary view:**
- Click "Summary" tab instead of "All Data"

3. **Upgrade server:**
```bash
# Allocate more memory
raptor dashboard --server.maxMemory 4000
```

---

### Can I embed dashboard in reports?

**Yes!** Export interactive HTML:

```bash
raptor dashboard --export report.html
```

Or embed specific plots:
```python
from raptor.dashboard import export_plot
export_plot(results, plot_type="volcano", output="volcano.html")
```

---

##  Pipeline Selection

### Which pipeline should I use?

**Use ML recommendations!**

```bash
raptor ml recommend --config config.yaml
```

**Or follow these guidelines:**

**Best all-around:** Salmon
- Fast, accurate, easy to use
- Good for most projects

**Most accurate quantification:** Salmon or Kallisto
- Quasi-mapping algorithms
- Similar results, both excellent

**Traditional mapping:** STAR + RSEM
- Required by some journals
- Slower but thorough

**Novel transcripts:** StringTie + Salmon
- De novo transcript discovery

[Detailed comparison](docs/PIPELINE_COMPARISON.md)

---

### Can I use multiple pipelines?

**Absolutely! Recommended!**

**Ensemble mode:**
```yaml
ensemble:
  enabled: true
  pipelines:
    - salmon
    - kallisto
    - star_rsem
  method: "weighted_average"
```

**Benefits:**
- ✅ More robust results
- ✅ Catch pipeline-specific artifacts
- ✅ Increase confidence
- ✅ Satisfy reviewers

---

### What about legacy pipelines?

**RAPTOR supports:**
- ✅ Salmon (2019+)
- ✅ Kallisto (2016+)
- ✅ STAR (2013+)
- ✅ RSEM (2011+)
- ✅ HTSeq (2010+)
- ✅ HISAT2 (2015+)

Older tools (TopHat, Bowtie1) not recommended - outdated algorithms.

---

##  Resource Requirements

### How much memory do I need?

**Depends on genome and pipeline:**

**Human genome:**
- Salmon: 8-16 GB
- Kallisto: 8-12 GB
- STAR: 32-48 GB

**Mouse genome:**
- Salmon: 6-12 GB
- Kallisto: 6-10 GB
- STAR: 32-40 GB

**Reduce memory:**
```yaml
resources:
  star_index_sparse: true  # Use sparse index for STAR
```

---

### How long does analysis take?

**Per sample (human, 50M reads):**

- Salmon: 5-10 minutes
- Kallisto: 3-7 minutes
- STAR: 15-30 minutes
- Full pipeline: 20-45 minutes

**For 100 samples:**
- Serial: 33-75 hours
- Parallel (16 cores): 2-5 hours
- Cloud: 1-2 hours

---

### Can I run on a laptop?

**Yes, but...**

**✅ Works fine:**
- Small projects (<20 samples)
- Fast pipelines (Salmon/Kallisto)
- Lower model organisms

**❌ Challenging:**
- Large projects (100+ samples)
- STAR alignment (needs 32+ GB RAM)
- Multiple pipelines simultaneously

**Recommendation:** Use cloud for large projects!

---

##  Cloud Computing

### Should I use cloud computing?

**Consider cloud if:**
- ✅ Large projects (50+ samples)
- ✅ No powerful local computer
- ✅ Need results quickly
- ✅ Want to test many pipelines
- ✅ Running infrequently

**Stick with local if:**
- ❌ Small projects (<20 samples)
- ❌ Have powerful workstation/cluster
- ❌ Privacy/security concerns
- ❌ Limited budget

---

### How much does cloud cost?

**Approximate costs:**

**AWS (50 samples, ensemble mode):**
- On-demand: $30-50
- Spot instances: $12-20

**GCP (50 samples):**
- Standard: $25-45
- Preemptible: $10-18

**Azure (50 samples):**
- Pay-as-you-go: $32-55
- Spot: $13-22

[Detailed pricing](docs/CLOUD_DEPLOYMENT.md#cost-estimation)

---

### How do I control cloud costs?

**8 Cost-Saving Strategies:**

1. **Use spot/preemptible instances** (70% savings)
2. **Set budget alerts**
3. **Auto-shutdown when done**
4. **Choose cheaper regions**
5. **Delete results after download**
6. **Compress data**
7. **Use cheaper storage tiers**
8. **Run during off-peak hours**

[Cost optimization guide](docs/CLOUD_DEPLOYMENT.md#cost-optimization)

---

### Is my data safe in the cloud?

**Security measures:**

1. **Encryption at rest and in transit**
2. **Private buckets/storage**
3. **IAM access controls**
4. **VPC isolation**
5. **Automatic deletion**

**Best practices:**
- De-identify sensitive data
- Use temporary credentials
- Enable logging
- Delete data when done
- Check compliance requirements

---

##  Results & Output

### What outputs does RAPTOR generate?

**Standard outputs:**
- Raw counts matrix
- Normalized counts (TPM/FPKM)
- Quality control metrics
- Pipeline comparison report
- Interactive plots

**v2.1.0 additions:**
- ML recommendation report
- Resource usage statistics
- Ensemble consensus results
- Automated publication-ready figures
- Interactive HTML dashboard

[Output formats guide](docs/OUTPUT_FORMATS.md)

---

### Can I use RAPTOR results in other tools?

**Yes! Compatible with:**
- ✅ DESeq2 (R)
- ✅ edgeR (R)
- ✅ limma (R)
- ✅ Seurat (R - via conversion)
- ✅ Scanpy (Python - via conversion)
- ✅ Pandas (Python)

**Export formats:**
- CSV
- TSV
- HDF5
- RDS (R)
- AnnData (Python)

```bash
raptor export --format deseq2 --output results.rds
```

---

### How do I share results with collaborators?

**Multiple options:**

1. **Dashboard link:**
```bash
raptor dashboard --host 0.0.0.0 --share
# Generates shareable link
```

2. **Static HTML report:**
```bash
raptor report --output report.html
# Email or host on website
```

3. **Cloud storage:**
```bash
raptor upload --destination s3://mybucket/results/
# Share bucket access
```

4. **GitHub:**
```bash
# Upload to repository
git add raptor_results/
git commit -m "Add RNA-seq results"
git push
```

---

### How long should I keep results?

**Recommended retention:**

- **Raw FASTQ:** 1-2 years (or until publication)
- **Aligned BAM:** 6 months - 1 year
- **Count matrices:** Keep forever (small files)
- **QC reports:** Keep forever
- **Intermediate files:** Delete after validation

**Storage tips:**
```yaml
output:
  keep_intermediate: false  # Auto-delete intermediate
  compress_outputs: true    # Save space
```

---

##  Troubleshooting

### My analysis failed - what do I do?

**Step-by-step troubleshooting:**

1. **Check the log:**
```bash
tail -n 50 raptor_error.log
```

2. **Run validation:**
```bash
raptor validate --config config.yaml
```

3. **Test with small dataset:**
```bash
raptor analyze --test-mode
```

4. **Enable debug mode:**
```yaml
logging:
  level: "DEBUG"
```

5. **Search GitHub issues:**
https://github.com/AyehBlk/RAPTOR/issues

[Full troubleshooting guide](TROUBLESHOOTING.md)

---

### Why is my analysis slow?

**Common causes:**

1. **Insufficient parallelization**
```yaml
resources:
  threads: 32  # Increase
  parallel_samples: 8
```

2. **Slow I/O**
- Use SSD instead of HDD
- Store data locally (not network drive)

3. **Memory swapping**
```bash
free -h  # Check if swapping
```

4. **Wrong pipeline choice**
- Use Salmon/Kallisto (faster)
- Avoid STAR for routine work

---

### Results don't match published data

**Possible reasons:**

1. **Different preprocessing**
- Trimming parameters
- Quality filtering

2. **Different normalization**
- TPM vs FPKM vs CPM
- Library size calculation

3. **Different genome version**
- GRCh37 vs GRCh38
- Gencode vs Ensembl

4. **Different filtering**
- Low-count gene filtering
- Biotype filtering

**Solution:** Match their methods exactly!

---

##  Advanced Usage

### Can I customize pipelines?

**Yes! Full customization:**

```yaml
custom_pipeline:
  name: "my_salmon"
  steps:
    - trim_galore:
        quality: 25
        length: 36
    - salmon:
        libType: "A"
        validateMappings: true
        gcBias: true
    - tximport:
        countsFromAbundance: "lengthScaledTPM"
```

[Pipeline customization guide](docs/CUSTOM_PIPELINES.md)

---

### Can I integrate with my pipeline?

**Yes! Multiple integration points:**

**Python API:**
```python
from raptor import RAPTORAnalysis
raptor = RAPTORAnalysis(config="config.yaml")
results = raptor.run()
```

**Command-line:**
```bash
# Use RAPTOR as a module
my_pipeline.sh | raptor process --stdin
```

**R integration:**
```r
# Import RAPTOR results
library(reticulate)
raptor <- import("raptor")
results <- raptor$load_results("raptor_results/")
```

---

### Can I add new features?

**Absolutely! RAPTOR is open-source:**

1. **Fork repository:**
```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
git checkout -b my-feature
```

2. **Add feature:**
```python
# Edit source code
```

3. **Submit pull request:**
```bash
git push origin my-feature
# Create PR on GitHub
```

[Contributing guide](CONTRIBUTING.md)

---

### How do I benchmark new pipelines?

**Use RAPTOR's benchmarking mode:**

```bash
raptor benchmark \
  --pipelines salmon,kallisto,my_new_tool \
  --metrics accuracy,speed,memory \
  --test-data test_dataset/
```

**Generates:**
- Performance comparison
- Resource usage
- Accuracy metrics
- Publication-ready plots

---

##  Additional Resources

### Where can I learn more?

**Documentation:**
- [User Guide](USER_GUIDE.md) - Comprehensive tutorial
- [API Reference](API.md) - Python API docs
- [Configuration Guide](CONFIGURATION.md) - All settings explained

**Examples:**
- [Tutorial 1: Basic Analysis](examples/01_basic/)
- [Tutorial 2: Ensemble Mode](examples/02_ensemble/)
- [Tutorial 3: Cloud Deployment](examples/03_cloud/)

**Video Tutorials:**
- Coming soon!

---

### How do I get help?

**Quick questions:**
- Check this FAQ
- Read [Troubleshooting Guide](TROUBLESHOOTING.md)

**Bug reports:**
- Open [GitHub issue](https://github.com/AyehBlk/RAPTOR/issues)
- Include log files and config

**Feature requests:**
- Open [GitHub discussion](https://github.com/AyehBlk/RAPTOR/discussions)

**Direct contact:**
- Email: ayehbolouki1988@gmail.com

**Response time:** Usually 1-3 days

---

### Can I contribute?

**Yes please!** Contributions welcome:

-  Documentation improvements
-  Bug reports
-  Feature suggestions
-  Code contributions
-  Benchmark datasets
-  Tutorials

[Contributing guide](CONTRIBUTING.md)

---

### What's coming in future versions?

**Planned for v2.2.0:**
- Single-cell RNA-seq support
- Spatial transcriptomics
- Long-read RNA-seq (ONT/PacBio)
- Enhanced visualization
- Multi-language support

**Planned for v2.3.0:**
- Real-time analysis
- GPU acceleration
- Advanced QC features
- Integration with more tools

[Roadmap](ROADMAP.md)

---

### How can I stay updated?

**Stay in the loop:**

-  Star on [GitHub](https://github.com/AyehBlk/RAPTOR)
-  Watch releases
-  Subscribe to mailing list
-  Follow @AyehBolou (coming soon!)

---

##  Quick Reference

### Common Commands

```bash
# Install
pip install raptor-rnaseq

# Get ML recommendation
raptor ml recommend --config config.yaml

# Run analysis
raptor analyze --config config.yaml

# Start dashboard
raptor dashboard

# Generate report
raptor report --output results.html

# Export for DESeq2
raptor export --format deseq2

# Run in cloud
raptor cloud deploy --platform aws

# Get help
raptor --help
```

---

### Quick Troubleshooting

```bash
# Validate config
raptor validate

# Check logs
tail -f raptor.log

# Test installation
raptor test

# System info
raptor sysinfo

# Version
raptor --version
```

---

##  Contact

**Author:** Ayeh Bolouki  
**Email:** ayehbolouki1988@gmail.com  
**GitHub:** https://github.com/AyehBlk/RAPTOR  
**Version:** 2.1.0  
**License:** MIT

---

**Still have questions?**  
Don't hesitate to reach out! We're here to help. 💙

*"The only stupid question is the one not asked!"* 
