# ⚡ RAPTOR v2.1.0 Resource Monitoring Guide

**Real-Time Tracking of Computational Resources**

Monitor CPU, memory, disk, and network usage during RNA-seq analysis to optimize performance and predict resource needs.

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Real-Time Monitoring](#real-time-monitoring)
4. [Resource Metrics](#resource-metrics)
5. [Performance Analysis](#performance-analysis)
6. [Resource Prediction](#resource-prediction)
7. [Optimization Tips](#optimization-tips)
8. [Troubleshooting](#troubleshooting)
9. [Best Practices](#best-practices)

---

## 🎯 Overview

### What is Resource Monitoring?

Resource monitoring tracks how much CPU, memory, disk, and network your analysis uses in real-time. This helps you:

✅ **Optimize performance** - Identify bottlenecks  
✅ **Predict requirements** - Estimate needs for future runs  
✅ **Prevent failures** - Avoid out-of-memory errors  
✅ **Choose hardware** - Right-size your compute  
✅ **Track costs** - Cloud cost estimation  
✅ **Debug issues** - Find performance problems  

### What Gets Monitored

```
┌─────────────────────────────────────────────┐
│                                              │
│  CPU Usage          [████████░░] 82%        │
│  Memory             [██████░░░░] 58%        │
│  Disk I/O           [███░░░░░░░] 34 MB/s    │
│  Network            [█░░░░░░░░░] 12 MB/s    │
│                                              │
│  Per-Pipeline Breakdown:                    │
│  • Salmon:          CPU 45%, RAM 12GB       │
│  • edgeR:           CPU 15%, RAM 4GB        │
│                                              │
│  Estimated Remaining: 23 minutes            │
│                                              │
└─────────────────────────────────────────────┘
```

---

## ⚡ Quick Start

### Enable Monitoring (Automatic)

```bash
# Resource monitoring is enabled by default in v2.1.0
raptor analyze --config config.yaml

# Monitoring runs in background automatically
```

### View Live Monitoring

```bash
# Option 1: Terminal dashboard
raptor monitor --live

# Option 2: Web dashboard
raptor dashboard
# Open browser, click "Resource Monitor"

# Option 3: Real-time plot
raptor monitor --plot-live
```

### View Historical Data

```bash
# View past analysis
raptor monitor --results /path/to/results/

# Generate report
raptor monitor --report --output monitoring_report.html
```

---

## 📊 Real-Time Monitoring

### Terminal Dashboard

```bash
raptor monitor --live
```

**Output:**
```
╔════════════════════════════════════════════════════════════╗
║  RAPTOR Resource Monitor - Live View                       ║
╠════════════════════════════════════════════════════════════╣
║                                                             ║
║  Analysis: Salmon-edgeR                                    ║
║  Started: 2025-11-19 14:30:15                              ║
║  Runtime: 00:15:23                                         ║
║  Status: Running - Quantification (Sample 8/12)            ║
║                                                             ║
║  ┌─────────────────────────────────────────────────────┐  ║
║  │ CPU Usage (16 cores)                                │  ║
║  │ [████████████████░░░░] 82%                          │  ║
║  │                                                      │  ║
║  │ Per-Core: ▇▇▇▇▇▇▇▇▇▇▇▇▇▆▆▅                           │  ║
║  └─────────────────────────────────────────────────────┘  ║
║                                                             ║
║  ┌─────────────────────────────────────────────────────┐  ║
║  │ Memory Usage (32 GB available)                      │  ║
║  │ [██████████████░░░░░] 58% (18.6 GB used)            │  ║
║  │                                                      │  ║
║  │ Breakdown:                                          │  ║
║  │  • Salmon:    12.4 GB                               │  ║
║  │  • edgeR:      4.2 GB                               │  ║
║  │  • System:     2.0 GB                               │  ║
║  └─────────────────────────────────────────────────────┘  ║
║                                                             ║
║  ┌─────────────────────────────────────────────────────┐  ║
║  │ Disk I/O                                            │  ║
║  │ Read:  [████░░░░░░] 156 MB/s                        │  ║
║  │ Write: [██░░░░░░░░]  34 MB/s                        │  ║
║  │                                                      │  ║
║  │ Total Read:  24.3 GB                                │  ║
║  │ Total Write:  8.7 GB                                │  ║
║  └─────────────────────────────────────────────────────┘  ║
║                                                             ║
║  ┌─────────────────────────────────────────────────────┐  ║
║  │ Estimated Completion: 23 minutes                    │  ║
║  │ Peak Memory Expected: 24 GB                         │  ║
║  │ Disk Space Needed: ~45 GB                           │  ║
║  └─────────────────────────────────────────────────────┘  ║
║                                                             ║
║  Press 'q' to quit, 'r' to refresh, 'p' to pause         ║
╚════════════════════════════════════════════════════════════╝
```

### Web Dashboard

```bash
raptor dashboard
```

**Features:**
- 📊 Real-time graphs (updating every 5 seconds)
- 📈 Historical trends
- 🎯 Per-pipeline breakdown
- 🔍 Drill-down into specific samples
- 💾 Export data
- ⚠️ Alert configuration

**Dashboard View:**
```
┌──────────────────────────────────────────────────────┐
│ Resource Monitor                          [Live] ●   │
├──────────────────────────────────────────────────────┤
│                                                       │
│  CPU Usage (Last 10 minutes)                         │
│  100% ┤                                              │
│   80% ┤     ╭─╮   ╭──╮                              │
│   60% ┤  ╭──╯ ╰───╯  ╰─╮                            │
│   40% ┤──╯             ╰───                         │
│    0% └───────────────────────────────────          │
│                                                       │
│  Memory Usage (Last 10 minutes)                      │
│  32GB ┤                                              │
│  24GB ┤           ╭────────                         │
│  16GB ┤      ╭────╯                                 │
│   8GB ┤──────╯                                       │
│   0GB └───────────────────────────────────          │
│                                                       │
│  Current Status:                                     │
│  ├─ CPU: 82% (13/16 cores active)                   │
│  ├─ Memory: 18.6 GB / 32 GB (58%)                   │
│  ├─ Disk I/O: Read 156 MB/s, Write 34 MB/s          │
│  └─ Bottleneck: None detected ✓                     │
│                                                       │
└──────────────────────────────────────────────────────┘
```

---

## 📈 Resource Metrics

### CPU Metrics

**What's Tracked:**
```yaml
cpu:
  total_usage: 82%           # Overall CPU usage
  per_core_usage:            # Individual core usage
    - core_0: 95%
    - core_1: 87%
    - core_2: 92%
    # ... for all cores
  process_breakdown:
    salmon: 45%
    star: 0%
    edgeR: 15%
    system: 22%
  efficiency: 0.82           # How well parallelized (0-1)
  context_switches: 1234     # Process switching rate
```

**Key Insights:**
- **High usage (>90%)**: Good! Using available resources
- **Low usage (<50%)**: Inefficient, increase parallelism
- **Uneven cores**: Some cores idle, others maxed
- **High context switching**: Too many parallel processes

### Memory Metrics

**What's Tracked:**
```yaml
memory:
  total: 32 GB
  used: 18.6 GB
  available: 13.4 GB
  percent: 58%
  
  breakdown:
    salmon_quantification: 12.4 GB
    edgeR_analysis: 4.2 GB
    system_overhead: 2.0 GB
  
  peak_usage: 22.1 GB          # Maximum seen
  swap_used: 0 GB              # Good! No swapping
  
  predictions:
    estimated_peak: 24 GB
    safety_margin: 8 GB        # 32 - 24 = 8 GB buffer
```

**Key Insights:**
- **Swap usage > 0**: Bad! Need more RAM
- **Usage > 90%**: Risk of out-of-memory
- **Large safety margin**: Can reduce RAM or run more samples
- **Gradual increase**: Normal, stable
- **Sudden spikes**: Investigate which step

### Disk I/O Metrics

**What's Tracked:**
```yaml
disk:
  read_speed: 156 MB/s         # Current read speed
  write_speed: 34 MB/s         # Current write speed
  
  total_read: 24.3 GB
  total_written: 8.7 GB
  
  iops:                        # Operations per second
    read: 3456
    write: 892
  
  latency:
    read: 2.3 ms               # Average latency
    write: 5.1 ms
  
  by_location:
    /data: 18.2 GB read
    /tmp: 6.1 GB read/write
    /results: 8.7 GB written
```

**Key Insights:**
- **Low speeds (<100 MB/s)**: Disk bottleneck!
- **High latency (>10 ms)**: Slow disk or network storage
- **Many small operations**: Consider buffering
- **Temp directory heavy usage**: Ensure fast disk

### Network Metrics (Cloud/Cluster)

**What's Tracked:**
```yaml
network:
  download_speed: 12 MB/s
  upload_speed: 3 MB/s
  
  total_downloaded: 2.1 GB
  total_uploaded: 0.5 GB
  
  connections:
    s3_bucket: active
    reference_db: completed
  
  latency: 45 ms
```

---

## 🔍 Performance Analysis

### Bottleneck Detection

```bash
raptor monitor --analyze-bottlenecks
```

**Output:**
```
Performance Bottleneck Analysis
═══════════════════════════════

Current Analysis: Salmon-edgeR (Sample 8/12)

✅ CPU: No bottleneck
   ├─ Usage: 82% (optimal range)
   ├─ All cores active
   └─ Efficiency: 0.82 (good)

✅ Memory: No bottleneck
   ├─ Usage: 58% (safe)
   ├─ No swapping
   └─ Peak expected: 24GB (within limits)

⚠️  DISK I/O: Minor bottleneck detected
   ├─ Read speed: 156 MB/s
   ├─ Expected: 300+ MB/s
   └─ Recommendation: Use SSD or local storage
      (Currently using network storage)

✅ Network: No bottleneck

Overall Performance: 8.5/10 ⭐
Primary Limitation: Disk I/O
Estimated Speedup with SSD: 25-30%
```

### Per-Pipeline Comparison

```bash
raptor monitor --compare-pipelines
```

**Output:**
```
Pipeline Resource Comparison
═══════════════════════════

Pipeline          Peak CPU  Peak Memory  Runtime   Efficiency
──────────────────────────────────────────────────────────────
Salmon-edgeR      95%       14 GB        22 min    ⭐⭐⭐⭐⭐
Kallisto-Sleuth   88%       8 GB         15 min    ⭐⭐⭐⭐⭐
STAR-RSEM         78%       45 GB        3.5 hr    ⭐⭐⭐☆☆
STAR-HTSeq        76%       42 GB        3.2 hr    ⭐⭐⭐☆☆

Recommendations:
• Salmon-edgeR: Excellent resource efficiency
• STAR-based: CPU underutilized, increase threads
• For 32GB RAM: Avoid STAR-RSEM, use Salmon
```

### Timeline Analysis

```bash
raptor monitor --timeline --output timeline.html
```

**Generates interactive timeline showing:**
```
Time →
0min     10min    20min    30min
├────────┼────────┼────────┼────────┤

CPU:     [██████████████████████░░]
Memory:  [░░░░████████████████░░░░]
Disk:    [████████░░░░░░░░░░░░████]

Phases:
├─ 0-8min:   Index loading    [High disk, low CPU]
├─ 8-25min:  Quantification   [High CPU, medium mem]
└─ 25-32min: DE analysis      [Medium CPU, high mem]

Bottlenecks:
⚠️  0-5min: Disk read limited (150 MB/s max)
✅  8-25min: Optimal performance
⚠️  25-30min: Memory approaching limit (28/32 GB)
```

---

## 🔮 Resource Prediction

### Predict Requirements for New Analysis

```bash
raptor monitor predict \
  --samples 50 \
  --pipeline 3 \
  --based-on previous_run_results/
```

**Output:**
```
Resource Prediction for 50 Samples (Salmon-edgeR)
═══════════════════════════════════════════════════

Based on analysis of: 12 samples (reference run)

Predicted Requirements:
┌──────────────────────────────────────────────┐
│ CPU Cores:        16 recommended             │
│                   (32 optimal)               │
│                                               │
│ Memory:           24-28 GB peak              │
│                   32 GB recommended          │
│                   (with 4GB safety margin)   │
│                                               │
│ Disk Space:       180-200 GB                 │
│                   ├─ Input: 85 GB            │
│                   ├─ Temp: 45 GB             │
│                   └─ Output: 50 GB           │
│                                               │
│ Runtime:          45-60 minutes              │
│                   (with 16 cores)            │
│                                               │
│ Network:          2.5 GB download            │
│ (if cloud)        0.8 GB upload              │
└──────────────────────────────────────────────┘

Confidence: High (based on similar analysis)

Suitable Systems:
✅ Your current system (32GB RAM, 16 cores)
✅ AWS r5.2xlarge
✅ GCP n2-highmem-8
⚠️  NOT recommended: Systems with <24GB RAM

Estimated Costs:
• Local: Free (you have resources)
• AWS (spot): $8-12
• GCP (preemptible): $6-10
```

### Scaling Predictions

```bash
# How would resources scale with more samples?
raptor monitor predict-scaling \
  --samples 10,20,50,100,200 \
  --pipeline 3
```

**Output:**
```
Resource Scaling Analysis (Salmon-edgeR)
═══════════════════════════════════════

Samples  CPU    Memory    Disk      Runtime
──────────────────────────────────────────────
10       45%    12 GB     80 GB     18 min
20       65%    16 GB     140 GB    35 min
50       85%    24 GB     320 GB    85 min
100      90%    32 GB     620 GB    165 min
200      95%    48 GB     1.2 TB    320 min

Key Observations:
• Memory: Linear scaling (~0.5 GB per sample)
• Disk: Linear scaling (~6 GB per sample)
• Runtime: Linear with samples (1.6 min/sample)
• CPU: Plateaus at ~16 cores

Recommendations:
• 10-50 samples: Current system perfect
• 100 samples: Upgrade to 32GB minimum
• 200 samples: Consider 64GB or cloud
```

---

## 💡 Optimization Tips

### CPU Optimization

**Problem: Low CPU usage (<50%)**
```bash
# Solution 1: Increase threads
raptor analyze --threads 32  # Instead of 16

# Solution 2: Parallel samples
raptor analyze --parallel-samples 4

# Solution 3: Use more efficient pipeline
raptor analyze --pipeline 3  # Salmon uses CPU better than STAR
```

**Problem: Uneven CPU usage**
```bash
# Check if hyperthreading is on
lscpu | grep Thread

# If hyperthreading: Use physical cores only
raptor analyze --threads 16 --no-hyperthreading
```

### Memory Optimization

**Problem: Using too much memory**
```bash
# Solution 1: Process samples sequentially
raptor analyze --sequential

# Solution 2: Reduce memory per job
raptor analyze --max-memory-per-sample 8G

# Solution 3: Use memory-efficient pipeline
raptor analyze --pipeline 3  # Salmon uses less RAM than STAR
```

**Problem: Running out of memory**
```bash
# Enable swap (emergency only)
sudo fallocate -l 8G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# Better: Upgrade RAM or use cloud
raptor cloud deploy --instance-type r5.4xlarge
```

### Disk I/O Optimization

**Problem: Slow disk access**
```bash
# Solution 1: Use SSD for temp files
raptor analyze --temp-dir /mnt/ssd/temp/

# Solution 2: Increase buffer size
raptor analyze --io-buffer-size 128M

# Solution 3: Compress on the fly
raptor analyze --compress-intermediate

# Solution 4: Use ramdisk for temp (if enough RAM)
sudo mkdir /mnt/ramdisk
sudo mount -t tmpfs -o size=32G tmpfs /mnt/ramdisk
raptor analyze --temp-dir /mnt/ramdisk/
```

### Pipeline-Specific Optimizations

**Salmon:**
```bash
raptor analyze --pipeline 3 \
  --salmon-threads 16 \
  --salmon-no-version-check \
  --salmon-no-frag-length-dist
```

**STAR:**
```bash
raptor analyze --pipeline 1 \
  --star-threads 16 \
  --star-limit-bam-sort-ram 32000000000 \
  --star-genomeSAsparseD 2  # Reduce memory
```

---

## 🔧 Troubleshooting

### Issue: "Out of Memory" Error

**Symptoms:**
```
Error: std::bad_alloc
Killed
Process terminated: Out of memory
```

**Solutions:**
```bash
# 1. Check current memory usage
raptor monitor --check-memory

# 2. Reduce parallel jobs
raptor analyze --parallel-samples 1

# 3. Use memory-efficient pipeline
raptor analyze --pipeline 3  # or 4

# 4. Enable sequential processing
raptor analyze --sequential

# 5. Increase swap (temporary)
sudo swapon --show
sudo fallocate -l 16G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# 6. Use cloud with more RAM
raptor cloud deploy --instance-type r5.4xlarge
```

### Issue: Very Slow Performance

**Diagnosis:**
```bash
raptor monitor --diagnose
```

**Possible causes and fixes:**
```
1. Disk Bottleneck (I/O wait >30%)
   → Use SSD
   → Move to local disk (not NFS)
   
2. Memory Swapping
   → Check: free -h
   → Fix: Add RAM or reduce parallel jobs
   
3. Network Storage
   → Check: df -h (is it NFS?)
   → Fix: Copy data locally first
   
4. CPU Underutilized
   → Check: top
   → Fix: Increase --threads

5. Low Parallelization
   → Fix: Use --parallel-samples
```

### Issue: Monitoring Not Working

**Solutions:**
```bash
# 1. Check if psutil is installed
pip list | grep psutil

# 2. Install/upgrade psutil
pip install --upgrade psutil

# 3. Check permissions (on Linux)
# Monitoring needs read access to /proc/
ls -l /proc/

# 4. Enable verbose logging
raptor analyze --monitor-verbose

# 5. Manual monitoring
# Use system tools instead
htop
iotop -o
```

---

## 📚 Best Practices

### Before Analysis

✅ **Predict requirements**
```bash
raptor monitor predict --samples 50 --pipeline 3
```

✅ **Check available resources**
```bash
# Memory
free -h

# Disk space
df -h

# CPU
lscpu
```

✅ **Use appropriate hardware**
- SSD for temp files
- Local storage (not network)
- Sufficient RAM (with 25% buffer)

### During Analysis

✅ **Monitor in real-time**
```bash
raptor monitor --live
```

✅ **Watch for warnings**
- Memory >90%: Risk of crash
- Swap usage >0: Performance hit
- CPU <50%: Inefficient

✅ **Enable alerts**
```bash
raptor monitor --alert-memory 90 \
  --alert-disk 95 \
  --notify your.email@example.com
```

### After Analysis

✅ **Review resource usage**
```bash
raptor monitor --results results/ --report
```

✅ **Learn for next time**
- Peak memory used
- Actual runtime
- Bottlenecks identified

✅ **Archive monitoring data**
```bash
raptor monitor --export results/monitoring_data.json
```

---

## 📊 Example Monitoring Report

```bash
raptor monitor --report --output report.html
```

**Generated Report Contents:**

```markdown
# Resource Monitoring Report

## Analysis Summary
- Pipeline: Salmon-edgeR
- Samples: 12
- Runtime: 32 minutes
- Date: 2025-11-19

## Resource Usage

### CPU
- Average: 78%
- Peak: 95%
- Efficiency: 0.82 (good)
- Bottleneck: None

### Memory
- Average: 14.2 GB
- Peak: 18.6 GB
- Max Available: 32 GB
- Safety Margin: 13.4 GB ✅

### Disk I/O
- Total Read: 24.3 GB
- Total Written: 8.7 GB
- Average Read Speed: 165 MB/s
- Average Write Speed: 38 MB/s
- Bottleneck: Minor (network storage)

## Timeline
[Interactive chart showing resource usage over time]

## Recommendations
1. Consider local SSD for 25-30% speedup
2. Current hardware well-suited for this analysis
3. Could handle up to 50 samples with same resources

## Predictions for Larger Runs
- 50 samples: Need 24-28 GB RAM, ~85 minutes
- 100 samples: Need 32-36 GB RAM, ~165 minutes
```

---

## 🎉 Summary

Resource monitoring provides:
- ✅ **Real-time visibility** - Know what's happening
- ✅ **Bottleneck identification** - Find slowdowns
- ✅ **Resource prediction** - Plan future analyses
- ✅ **Cost estimation** - Budget for cloud
- ✅ **Performance optimization** - Run faster
- ✅ **Failure prevention** - Avoid crashes

**Monitor your analysis, optimize your workflow!** 📊⚡

---

**Author:** Ayeh Bolouki  
**Affiliation:** University of Namur & GIGA-Neurosciences, University of Liège, Belgium  
**Version:** 2.1.0  
**License:** MIT

---

*"You can't optimize what you don't measure!"* ⚡📊
