# 🔧 RAPTOR v2.1.0 Troubleshooting Guide

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

A comprehensive guide to diagnosing and resolving common issues with RAPTOR v2.1.0.

---

## 📋 Table of Contents

1. [Installation Issues](#installation-issues)
2. [Configuration Problems](#configuration-problems)
3. [ML Recommendation Issues](#ml-recommendation-issues)
4. [Dashboard Problems](#dashboard-problems)
5. [Pipeline Execution Errors](#pipeline-execution-errors)
6. [Resource Monitoring Issues](#resource-monitoring-issues)
7. [Ensemble Analysis Problems](#ensemble-analysis-problems)
8. [Data Quality Issues](#data-quality-issues)
9. [Cloud Deployment Issues](#cloud-deployment-issues)
10. [Performance Problems](#performance-problems)

---

## 🔧 Installation Issues

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
pip install git+https://github.com/AyehBlk/RAPTOR.git@v2.1.0
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

3. **Use requirements file:**
```bash
pip install -r requirements.txt
```

---

### Issue: Missing system dependencies

**Symptoms:**
```bash
ImportError: libgomp.so.1: cannot open shared object file
```

**Solutions:**

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential libgomp1
```

**CentOS/RHEL:**
```bash
sudo yum install gcc gcc-c++ libgomp
```

**macOS:**
```bash
brew install gcc libomp
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

3. **Check working directory:**
```bash
pwd  # Should be in project root
ls -la  # Should see config.yaml
```

---

### Issue: Invalid configuration

**Symptoms:**
```bash
ValueError: Invalid pipeline name: slamonm
```

**Solutions:**

1. **Validate config:**
```bash
raptor config --validate
```

2. **Check pipeline names:**
```yaml
pipelines:
  - "salmon"  # Correct
  - "star"    # Correct
  # NOT "slamonm" or "str"
```

3. **Use config template:**
```bash
raptor config --template > config.yaml
```

---

### Issue: Path errors in config

**Symptoms:**
```bash
FileNotFoundError: /data/fastq/sample.fq.gz not found
```

**Solutions:**

1. **Use absolute paths:**
```yaml
data:
  fastq_dir: "/home/user/project/data/fastq"  # Not "~/data/fastq"
```

2. **Verify paths exist:**
```bash
ls -l /data/fastq/  # Check before running
```

3. **Use wildcards correctly:**
```yaml
data:
  fastq_pattern: "*_R{1,2}.fastq.gz"  # Correct
  # NOT: "*.fastq.gz"  # Too broad
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

3. **Specify model path:**
```yaml
ml_recommendation:
  model_path: "/path/to/model.pkl"
```

---

### Issue: Low confidence predictions

**Symptoms:**
```
Warning: Low confidence (45%) - manual review recommended
```

**Solutions:**

1. **Provide more metadata:**
```yaml
metadata:
  organism: "Homo sapiens"
  tissue_type: "brain"
  read_length: 150
  sequencing_depth: "50M"
  library_prep: "polyA"
```

2. **Retrain with similar data:**
```bash
raptor ml train --data my_similar_projects.csv
```

3. **Use ensemble mode:**
```yaml
ml_recommendation:
  use_ensemble: true
  min_confidence: 0.6
```

---

### Issue: Prediction takes too long

**Symptoms:**
```
ML prediction running for > 5 minutes...
```

**Solutions:**

1. **Reduce feature extraction:**
```yaml
ml_recommendation:
  quick_mode: true
  skip_advanced_features: true
```

2. **Use cached predictions:**
```yaml
ml_recommendation:
  cache_predictions: true
```

3. **Increase timeout:**
```yaml
ml_recommendation:
  timeout: 600  # 10 minutes
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

3. **Use different host:**
```bash
raptor dashboard --host 0.0.0.0 --port 8501
```

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

3. **Check output path:**
```yaml
output:
  base_dir: "./raptor_results"  # Ensure this exists
```

---

### Issue: Plots not rendering

**Symptoms:**
```
Error loading plot: JavaScript heap out of memory
```

**Solutions:**

1. **Reduce plot complexity:**
```yaml
dashboard:
  max_points_per_plot: 1000
  simplify_plots: true
```

2. **Clear browser cache:**
- Chrome: Ctrl+Shift+Del
- Firefox: Ctrl+Shift+Del
- Safari: Cmd+Option+E

3. **Update browsers:**
```bash
# Latest Chrome, Firefox, or Safari recommended
```

---

### Issue: Dashboard slow/unresponsive

**Symptoms:**
```
Dashboard loads but interactions are slow
```

**Solutions:**

1. **Limit data displayed:**
```yaml
dashboard:
  max_samples_displayed: 50
  lazy_loading: true
```

2. **Increase memory:**
```bash
raptor dashboard --server.maxUploadSize 1000
```

3. **Use lighter theme:**
```yaml
dashboard:
  theme: "minimal"
  disable_animations: true
```

---

## 🚀 Pipeline Execution Errors

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

3. **Install missing tools:**
```bash
conda install -c bioconda salmon star kallisto
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
  java_opts: "-Xmx24g"
```

2. **Reduce parallelism:**
```yaml
resources:
  threads: 8  # Reduce from 16
```

3. **Process samples in batches:**
```yaml
execution:
  batch_size: 10
  sequential_mode: true
```

---

### Issue: Pipeline hangs/freezes

**Symptoms:**
```
Pipeline running for hours with no progress...
```

**Solutions:**

1. **Enable verbose logging:**
```yaml
logging:
  level: "DEBUG"
  log_file: "raptor_debug.log"
```

2. **Check system resources:**
```bash
top  # Check CPU/memory
df -h  # Check disk space
```

3. **Set timeout:**
```yaml
execution:
  max_runtime: 3600  # 1 hour per sample
  timeout_action: "skip"
```

---

### Issue: Index not found

**Symptoms:**
```bash
Error: Genome index not found at /data/genome/index
```

**Solutions:**

1. **Build index:**
```bash
raptor index --genome genome.fa --gtf genes.gtf
```

2. **Download pre-built:**
```bash
raptor index --download hg38
```

3. **Specify correct path:**
```yaml
reference:
  index: "/data/genomes/hg38/salmon_index"
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

2. **Enable monitoring:**
```yaml
resource_monitoring:
  enabled: true
  interval: 10
```

3. **Check permissions:**
```bash
# On Linux, may need to run with elevated permissions
sudo raptor analyze --config config.yaml
```

---

### Issue: Inaccurate metrics

**Symptoms:**
```
CPU: 0%, Memory: 0% (clearly wrong)
```

**Solutions:**

1. **Restart monitoring:**
```bash
raptor monitor --restart
```

2. **Update psutil:**
```bash
pip install --upgrade psutil
```

3. **Platform-specific fix:**

**macOS:**
```bash
sudo purge  # Clear cached memory stats
```

**Linux:**
```bash
sync; echo 3 > /proc/sys/vm/drop_caches
```

---

### Issue: Monitoring overhead

**Symptoms:**
```
Analysis slower with monitoring enabled
```

**Solutions:**

1. **Reduce monitoring frequency:**
```yaml
resource_monitoring:
  interval: 60  # Check every 60 seconds instead of 10
```

2. **Disable detailed tracking:**
```yaml
resource_monitoring:
  track_per_process: false
  detailed_metrics: false
```

3. **Use lightweight mode:**
```yaml
resource_monitoring:
  lightweight_mode: true
```

---

## 🔬 Ensemble Analysis Problems

### Issue: Ensemble fails to combine

**Symptoms:**
```bash
Error: Cannot combine results from different pipelines
```

**Solutions:**

1. **Ensure consistent settings:**
```yaml
ensemble:
  normalize_method: "TMM"  # Same for all pipelines
  filter_low_counts: true
```

2. **Check output formats:**
```bash
# All pipelines must produce compatible output
ls raptor_results/salmon/
ls raptor_results/kallisto/
```

3. **Validate individually first:**
```bash
raptor analyze --pipeline salmon --validate-only
raptor analyze --pipeline kallisto --validate-only
```

---

### Issue: Ensemble results inconsistent

**Symptoms:**
```
Warning: High variance between pipeline results
```

**Solutions:**

1. **Check input data quality:**
```bash
raptor qc --fastq data/fastq/
```

2. **Review pipeline parameters:**
```yaml
ensemble:
  outlier_detection: true
  consensus_threshold: 0.7
```

3. **Investigate discrepancies:**
```python
from raptor.ensemble import compare_pipelines
compare_pipelines(['salmon', 'kallisto', 'star'])
```

---

### Issue: Missing pipeline results

**Symptoms:**
```bash
Error: Ensemble requires ≥2 pipelines, found 1
```

**Solutions:**

1. **Run missing pipelines:**
```bash
raptor analyze --pipeline kallisto --ensemble-mode
```

2. **Check for failed pipelines:**
```bash
ls raptor_results/*/
grep "ERROR" raptor.log
```

3. **Lower minimum requirement:**
```yaml
ensemble:
  min_pipelines: 2  # Ensure this matches your runs
```

---

## 🔍 Data Quality Issues

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

2. **Check FastQC reports:**
```bash
fastqc data/fastq/*.fastq.gz
```

3. **Trim low-quality reads:**
```yaml
preprocessing:
  trim_adapters: true
  quality_threshold: 20
  min_length: 36
```

---

### Issue: Contamination detected

**Symptoms:**
```
Warning: Possible contamination - 15% reads unmapped
```

**Solutions:**

1. **Run contamination screen:**
```bash
raptor qc --screen-contamination
```

2. **Check species:**
```yaml
quality_assessment:
  expected_species: "Homo sapiens"
  contamination_threshold: 0.05
```

3. **Filter contaminated reads:**
```bash
# Use tools like BBDuk or Kraken2
```

---

### Issue: Batch effects detected

**Symptoms:**
```
Warning: Significant batch effects detected
```

**Solutions:**

1. **Include batch info:**
```yaml
metadata:
  batch_column: "sequencing_run"
  correct_batch: true
```

2. **Use ComBat:**
```python
from raptor.quality import correct_batch_effects
corrected_data = correct_batch_effects(data, batch_info)
```

3. **Include in design:**
```yaml
analysis:
  design: "~ batch + condition"
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
# Enter Access Key ID and Secret Access Key
```

**GCP:**
```bash
gcloud auth login
gcloud config set project PROJECT_ID
```

**Azure:**
```bash
az login
az account set --subscription SUBSCRIPTION_ID
```

---

### Issue: Insufficient permissions

**Symptoms:**
```bash
Error: Access denied to S3 bucket
```

**Solutions:**

1. **Check IAM permissions:**
```bash
aws iam get-user
aws iam list-attached-user-policies --user-name USERNAME
```

2. **Add required permissions:**
```json
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Action": [
      "s3:GetObject",
      "s3:PutObject",
      "batch:SubmitJob"
    ],
    "Resource": "*"
  }]
}
```

---

### Issue: Instance launch failed

**Symptoms:**
```bash
Error: Could not launch instance - quota exceeded
```

**Solutions:**

1. **Request quota increase:**
```bash
# AWS
aws service-quotas request-service-quota-increase \
  --service-code batch --quota-code L-34F8E9F6 \
  --desired-value 100
```

2. **Use different region:**
```yaml
cloud:
  aws:
    region: "us-west-2"  # Try different region
```

3. **Use smaller instances:**
```yaml
cloud:
  aws:
    instance_type: "m5.large"  # Instead of m5.4xlarge
```

---

### Issue: High cloud costs

**Symptoms:**
```
Unexpected AWS bill of $500 for 100 samples
```

**Solutions:**

1. **Use spot instances:**
```yaml
cloud:
  aws:
    use_spot: true
    max_spot_price: 0.50
```

2. **Set budget alerts:**
```bash
aws budgets create-budget --budget file://budget.json
```

3. **Enable auto-shutdown:**
```yaml
cloud:
  auto_shutdown: true
  max_idle_time: 1800  # 30 minutes
```

---

### Issue: Data transfer slow

**Symptoms:**
```
Data upload: 0.5 MB/s (expected 10+ MB/s)
```

**Solutions:**

1. **Use transfer acceleration:**
```yaml
cloud:
  aws:
    s3_transfer_acceleration: true
```

2. **Compress data:**
```bash
tar -czf data.tar.gz data/
```

3. **Parallel uploads:**
```bash
aws s3 sync data/ s3://bucket/ --parallel
```

---

## ⚡ Performance Problems

### Issue: Analysis very slow

**Symptoms:**
```
Expected: 2 hours, Actual: 12 hours for 50 samples
```

**Solutions:**

1. **Increase parallelism:**
```yaml
resources:
  threads: 32  # Use more cores
  parallel_samples: 8
```

2. **Use faster pipelines:**
```yaml
pipelines:
  - "salmon"  # Faster than STAR
  - "kallisto"  # Faster than RSEM
```

3. **Reduce replicates:**
```yaml
parameter_optimization:
  bootstrap_iterations: 50  # Instead of 100
```

---

### Issue: High memory usage

**Symptoms:**
```
System using 95% RAM, causing swapping
```

**Solutions:**

1. **Limit memory per process:**
```yaml
resources:
  max_memory_per_job: "8GB"
```

2. **Process in batches:**
```yaml
execution:
  batch_size: 5
  sequential_mode: true
```

3. **Use memory-efficient tools:**
```yaml
pipelines:
  - "salmon"  # More memory-efficient than STAR
```

---

### Issue: Disk space running out

**Symptoms:**
```bash
Error: No space left on device
```

**Solutions:**

1. **Clean up intermediate files:**
```yaml
output:
  keep_intermediate: false
  auto_cleanup: true
```

2. **Compress outputs:**
```yaml
output:
  compress_outputs: true
  compression_level: 9
```

3. **Monitor disk usage:**
```bash
du -sh raptor_results/
df -h
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
  verbose: true
```

3. **Collect system info:**
```bash
raptor --version
python --version
pip list | grep raptor
uname -a
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

4. **System info:**
```bash
raptor sysinfo
```

5. **Steps to reproduce:**
```
1. Set up config.yaml with...
2. Run: raptor analyze...
3. Error occurs at...
```

### Support Channels

- **GitHub Issues:** https://github.com/AyehBlk/RAPTOR/issues
- **Email:** ayehbolouki1988@gmail.com
- **Documentation:** https://github.com/AyehBlk/RAPTOR/wiki

---

## 📚 Additional Resources

- [Installation Guide](INSTALLATION.md)
- [Configuration Guide](CONFIGURATION.md)
- [API Documentation](API.md)
- [FAQ](FAQ.md)

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**Version:** 2.1.0  
**License:** MIT

---

*"Every bug is a feature waiting to be understood."* 🐛✨
