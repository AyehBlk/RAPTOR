# RAPTOR 🦖

**RNA-seq Analysis Pipeline Testing and Optimization Resource**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![GitHub Stars](https://img.shields.io/github/stars/AyehBlk/RAPTOR?style=social)](https://github.com/AyehBlk/RAPTOR)

> **Making RNA-seq pipeline selection evidence-based, not guesswork.**

Created by **Ayeh Bolouki**  
📧 ayehbolouki1988@gmail.com | ayehgeek@gmail.com  
🏛️ Postdoctoral Researcher, GIGA Center, University of Liège, Belgium

---

## 📖 Table of Contents

- [The Challenge](#-the-challenge)
- [The RAPTOR Solution](#-the-raptor-solution)
- [Quick Start](#-quick-start)
- [Pipeline Implementations](#-pipeline-implementations)
- [Benchmarking Framework](#-benchmarking-framework)
- [Intelligent Recommendations (NEW!)](#-intelligent-recommendations-new)
- [Installation](#-installation)
- [Usage Examples](#-usage-examples)
- [Documentation](#-documentation)
- [Citation](#-citation)
- [Contributing](#-contributing)
- [License](#-license)

---

## 🔬 The Challenge

Over **50+ bioinformatics tools** exist for RNA-seq differential expression analysis, and choosing the right combination is overwhelming:

- **Alignment**: STAR? HISAT2? Salmon? Kallisto?
- **Quantification**: RSEM? StringTie? HTSeq? featureCounts?
- **Statistics**: DESeq2? edgeR? limma? Sleuth? Ballgown?

Different pipelines have different strengths and trade-offs:
- Some prioritize **accuracy** but are slow
- Others are **ultra-fast** but may miss subtle differences
- Some work better with **small sample sizes**
- Others excel with **large datasets**

**The real question:** *"Which pipeline should I use for MY data?"*

Traditional approaches:
- ❌ **Guesswork** based on lab traditions
- ❌ **Following tutorials** without understanding trade-offs
- ❌ **Spending weeks** benchmarking manually
- ❌ **Hoping** you made the right choice

---

## 🦖 The RAPTOR Solution

RAPTOR provides **two complementary approaches** to help you make informed, evidence-based decisions:

### 1️⃣ Comprehensive Pipeline Benchmarking

**What it does:**
- Implements **8 complete RNA-seq pipelines** spanning different methodological approaches
- Tests each pipeline on **standardized datasets** with known ground truth
- Measures **accuracy** (sensitivity, precision, F1), **efficiency** (speed, memory), and **reproducibility**
- Generates **publication-quality comparison reports**

**Why it matters:**
- Understand the **objective performance** of each pipeline
- See **real trade-offs** between accuracy, speed, and resources
- Make decisions based on **evidence**, not assumptions
- Validate pipeline choices for your publications

### 2️⃣ Intelligent Data-Driven Recommendations ✨ NEW!

**What it does:**
- Analyzes **20+ statistical characteristics** of YOUR data
- Profiles library sizes, zero-inflation, biological variation, sample size, depth
- Matches data characteristics to pipeline strengths
- Provides **personalized recommendation** with clear reasoning

**Why it matters:**
- Get the **optimal pipeline for YOUR specific data** in 30 seconds
- Understand **WHY** it's recommended
- No need to benchmark all pipelines manually
- Confidence in your analytical choices

---

##  Quick Start

Choose your path based on what you need:

###  Path A: Fast Track - Get Recommendation (Recommended)

**Best for:** Most users who want quick, data-driven guidance

```bash
# Install
pip install raptor-rnaseq

# Analyze your data and get recommendation
raptor profile --counts your_counts.csv --metadata metadata.csv --output report.html

# Opens HTML report showing:
#  RECOMMENDED: Pipeline 3 (Salmon-edgeR)
#     Excellent for your large sample size
#     Fast methods sufficient at high depth
#     Good balance of speed and accuracy

# Run the recommended pipeline
raptor run --pipeline 3 --data your_fastq/ --output results/
```

**Time:** 30 seconds for recommendation + pipeline runtime

###  Path B: Deep Dive - Full Benchmarking

**Best for:** Method developers, publications requiring extensive validation

```bash
# Install with all dependencies
conda install -c bioconda raptor

# Run all 8 pipelines on your data
raptor compare --data your_fastq/ --output comparison/

# Generate comprehensive comparison report
raptor report --results comparison/ --output full_report.html
```

**Time:** Several hours (depending on data size)

###  Path C: Quick Demo

**Best for:** First-time users exploring RAPTOR

```bash
raptor demo
```

**Time:** 5 minutes

---

##  Pipeline Implementations

RAPTOR implements **8 complete RNA-seq analysis workflows**, each representing different methodological approaches used in the research community:

### Pipeline 1: STAR-RSEM-DESeq2 (Gold Standard)

**Approach:** Alignment-based with robust statistics

- **Alignment**: STAR (2-pass mode)
- **Quantification**: RSEM (EM algorithm)
- **Normalization**: DESeq2 median-of-ratios
- **Statistics**: DESeq2 (negative binomial with shrinkage)

**Strengths:**
- ✅ Highest accuracy for differential expression
- ✅ Excellent with low replication (n<3)
- ✅ Robust normalization for high variation
- ✅ Handles zero-inflation well

**Trade-offs:**
-  Slowest runtime (20-30 min for typical dataset)
-  Highest memory usage (~32GB)

**Best for:** Publication-quality results, difficult datasets, when accuracy is paramount

---

### Pipeline 2: HISAT2-StringTie-Ballgown

**Approach:** Alignment-based with transcript assembly

- **Alignment**: HISAT2
- **Assembly & Quantification**: StringTie
- **Normalization**: TPM
- **Statistics**: Ballgown

**Strengths:**
- ✅ Novel transcript discovery
- ✅ Isoform-level analysis
- ✅ Good for non-model organisms

**Trade-offs:**
-  Medium speed
-  Lower accuracy for gene-level DE

**Best for:** Transcriptome assembly, non-model organisms, isoform analysis

---

### Pipeline 3: Salmon-edgeR (Recommended for Most Users)

**Approach:** Alignment-free pseudo-alignment with robust statistics

- **Pseudo-alignment**: Salmon (quasi-mapping)
- **Import**: tximport (gene-level aggregation)
- **Normalization**: TMM (edgeR)
- **Statistics**: edgeR (quasi-likelihood F-test)

**Strengths:**
- ✅ Excellent accuracy (F1: 0.88)
- ✅ 3-5× faster than alignment-based
- ✅ Low memory usage (~8GB)
- ✅ Good balance of all metrics

**Trade-offs:**
-  Slightly lower accuracy than STAR-RSEM-DESeq2 for difficult data

**Best for:** Most RNA-seq experiments, large datasets, when speed and accuracy both matter

---

### Pipeline 4: Kallisto-Sleuth (Ultra-Fast)

**Approach:** Alignment-free with integrated statistics

- **Pseudo-alignment**: Kallisto
- **Normalization & Statistics**: Sleuth (integrated)

**Strengths:**
- ✅ Fastest option (5-10 minutes)
- ✅ Very low memory (~4GB)
- ✅ Good accuracy at high depth

**Trade-offs:**
-  Lower accuracy for low-depth or difficult data
-  May miss subtle differences

**Best for:** Exploratory analysis, large cohort studies, when speed is critical

---

### Pipeline 5: STAR-HTSeq-limma-voom

**Approach:** Alignment-based with flexible modeling

- **Alignment**: STAR
- **Quantification**: HTSeq
- **Normalization**: TMM + voom transformation
- **Statistics**: limma (empirical Bayes)

**Strengths:**
- ✅ Flexible statistical modeling
- ✅ Excellent for complex designs
- ✅ Handles batch effects well

**Trade-offs:**
-  Medium-slow speed
-  High memory usage

**Best for:** Complex experimental designs, batch correction, multi-factor analysis

---

### Pipeline 6: STAR-featureCounts-NOISeq

**Approach:** Non-parametric statistics

- **Alignment**: STAR
- **Quantification**: featureCounts
- **Statistics**: NOISeq (non-parametric)

**Strengths:**
- ✅ No distribution assumptions
- ✅ Robust to outliers

**Trade-offs:**
-  Lower power than parametric methods

**Best for:** Data not fitting standard distributions, exploratory analysis

---

### Pipeline 7: Bowtie2-RSEM-EBSeq

**Approach:** Bayesian statistical framework

- **Alignment**: Bowtie2
- **Quantification**: RSEM
- **Statistics**: EBSeq (empirical Bayes)

**Strengths:**
- ✅ Bayesian posterior probabilities
- ✅ Handles isoform uncertainty

**Trade-offs:**
-  Very slow
-  Conservative for gene-level analysis

**Best for:** Isoform switching analysis, when Bayesian framework preferred

---

### Pipeline 8: HISAT2-Cufflinks-Cuffdiff

**Approach:** Legacy pipeline for comparison

- **Alignment**: HISAT2
- **Assembly & Quantification**: Cufflinks
- **Statistics**: Cuffdiff

**Strengths:**
-  Historical reference
-  Widely published (older studies)

**Trade-offs:**
-  Outdated methodology
-  Lower accuracy than modern methods

**Best for:** Comparison with legacy studies, methodological benchmarks

---

##  Benchmarking Framework

RAPTOR evaluates pipelines across **three critical dimensions** to provide comprehensive performance assessment:

###  Accuracy Metrics

Evaluated on datasets with **known ground truth** (simulated data + SEQC/MAQC validated datasets):

- **Sensitivity (Recall)** - Ability to detect true positives
- **Specificity** - Ability to avoid false positives
- **Precision (PPV)** - Proportion of true positives among calls
- **F1 Score** - Harmonic mean of precision and recall
- **ROC AUC** - Overall discriminative ability

###  Efficiency Metrics

Measured on standardized hardware:

- **Wall-clock Runtime** - Total time to completion
- **Peak Memory Usage** - Maximum RAM required
- **CPU Utilization** - Multi-threading efficiency
- **Disk Space** - Intermediate and final file sizes

###  Reproducibility Metrics

Tested across independent runs:

- **Between-run Concordance** - DEG list overlap
- **Fold Change Correlation** - Consistency of effect sizes
- **P-value Correlation** - Statistical reproducibility

###  Benchmark Datasets

**1. Simulated Data (Polyester)**
- Controlled ground truth
- Known DE genes
- Variable parameters: fold changes, sample size, depth

**2. SEQC/MAQC Benchmark**
- Real RNA-seq data
- qPCR validation (gold standard)
- Published ground truth

**3. User Data (Optional)**
- Test on your own datasets
- Compare to your existing results

---

##  Intelligent Recommendations (NEW!)

The **Profile & Recommend** feature analyzes your data and suggests the optimal pipeline based on its characteristics.

### What Gets Analyzed

#### 1. Library Statistics
- Mean library size
- Library size variation (coefficient of variation)
- Dynamic range (max/min ratio)

#### 2. Count Distribution
- **Zero-inflation**: Percentage of zero counts
- **Low-count genes**: Genes with minimal expression
- **Expression range**: Log2 fold difference between highest and lowest

#### 3. Biological Variation
- **BCV (Biological Coefficient of Variation)**: Expected variation between replicates
- **Overdispersion**: Variance exceeding Poisson expectation
- **Mean-variance relationship**: How variability scales with expression

#### 4. Experimental Design
- Number of replicates per condition
- Total sample size
- Design balance
- Complexity (paired, multi-factor, etc.)

#### 5. Sequencing Characteristics
- Sequencing depth category (low/medium/high/very high)
- Coverage uniformity
- Overall quality indicators

### How Recommendations Work

RAPTOR uses a **scoring system** that matches data characteristics to pipeline strengths:

```
Score Components:
├── Data difficulty      (40%) - Variation, zero-inflation, sample size
├── Sequencing quality   (30%) - Depth, uniformity
├── Priority weight      (20%) - User-specified (accuracy/speed/memory)
└── Design complexity    (10%) - Experimental design considerations
```

Each pipeline receives a score (0-200) for your specific data, with **higher = better match**.

### Example Recommendations

#### Scenario 1: Challenging Dataset

**Data Characteristics:**
- High library size variation (CV > 0.5)
- Low replication (n=2 per group)
- High zero-inflation (>60% zeros)
- Low sequencing depth (<10M reads)

**→ RAPTOR Recommends: Pipeline 1 (STAR-RSEM-DESeq2)**

**Score: 165/200**

**Reasoning:**
- ✅ DESeq2's robust normalization handles high variation excellently
- ✅ Shrinkage estimators work well with low replication
- ✅ Alignment-based methods more accurate at low depth
- ✅ Handles zero-inflation better than count-based methods

---

#### Scenario 2: Large, High-Quality Dataset

**Data Characteristics:**
- Low library size variation (CV < 0.2)
- Good replication (n=10 per group)
- Low zero-inflation (<30% zeros)
- High sequencing depth (>40M reads)

**→ RAPTOR Recommends: Pipeline 3 (Salmon-edgeR)**

**Score: 145/200**

**Reasoning:**
- ✅ Pseudo-alignment sufficient at high depth and quality
- ✅ 3-5× faster for large sample sizes
- ✅ Low memory requirements enable parallel processing
- ✅ edgeR's TMM normalization excellent for low variation

---

#### Scenario 3: Speed-Critical Analysis

**Data Characteristics:**
- Moderate quality metrics
- Large cohort (n=50 samples)
- Time-sensitive project

**→ RAPTOR Recommends: Pipeline 4 (Kallisto-Sleuth)**

**Score: 130/200**

**Reasoning:**
- ✅ Fastest option (5-10 minutes vs hours)
- ✅ Minimal memory enables processing many samples
- ✅ Sufficient accuracy for well-powered studies
- ✅ Enables rapid iteration during exploratory phase

---

### Priority Modes

You can specify what matters most for your analysis:

```bash
# Maximize accuracy (default for difficult data)
raptor profile --counts data.csv --priority accuracy

# Maximize speed (large datasets, time constraints)
raptor profile --counts data.csv --priority speed

# Minimize memory (limited computational resources)
raptor profile --counts data.csv --priority memory

# Balanced approach (data-driven only, default)
raptor profile --counts data.csv --priority balanced
```

### Output Formats

1. **Console Summary** - Quick text-based recommendation
2. **HTML Report** - Interactive dashboard with:
   - Data profiling visualizations
   - Recommendation with reasoning
   - Alternative pipeline options
   - Quick-start commands
3. **JSON Export** - For programmatic access and automation

---

##  Installation

### System Requirements

- **OS**: Linux, macOS, WSL2 (Windows)
- **Python**: 3.8 or higher
- **Memory**: 8GB minimum, 32GB recommended for full benchmarking
- **Storage**: 10GB for core tools, 50GB for complete installation

### Quick Installation Options

#### Option 1: Profile & Recommend Only (Lightweight)

For users who only need data analysis and recommendations:

```bash
pip install raptor-rnaseq
```

**Size:** ~200MB  
**Includes:** Data profiling and recommendation system only  
**No pipeline tools** (STAR, Salmon, etc.) included

#### Option 2: Complete Installation (Conda)

For users who want full benchmarking capabilities:

```bash
# Create environment with all dependencies
conda env create -f environment.yml
conda activate raptor

# Verify installation
raptor --version
raptor test
```

**Size:** ~5GB  
**Includes:** All 8 pipelines with dependencies

#### Option 3: Docker (Portable)

For reproducible environments:

```bash
docker pull ayehblk/raptor:latest
docker run -v $(pwd):/data ayehblk/raptor raptor profile --counts /data/counts.csv
```

#### Option 4: From Source (Developers)

```bash
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR
pip install -e .
```

### Dependencies

**Core Python packages** (auto-installed):
```
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
scikit-learn>=0.24.0
```

**Pipeline tools** (for benchmarking, Conda installation):
```
STAR>=2.7.10
Salmon>=1.5.0
Kallisto>=0.48.0
HISAT2>=2.2.1
Bowtie2>=2.4.0
RSEM>=1.3.3
StringTie>=2.2.0
HTSeq>=2.0.0
featureCounts (subread)>=2.0.0
```

**R packages** (for statistics):
```
DESeq2>=1.34.0
edgeR>=3.36.0
limma>=3.50.0
Sleuth>=0.30.0
Ballgown>=2.26.0
NOISeq>=2.38.0
EBSeq>=1.34.0
```

---

##  Usage Examples

### Example 1: First-Time User

```bash
# Install lightweight version
pip install raptor-rnaseq

# Get recommendation for your data
raptor profile --counts my_counts.csv --metadata my_metadata.csv --output recommendation.html

# View the HTML report in your browser
# It shows: " Recommended: Pipeline 3 (Salmon-edgeR)"

# Run the recommended pipeline (requires full installation)
raptor run --pipeline 3 --data my_fastq/ --output results/
```

---

### Example 2: Compare Top Recommendations

```bash
# Get recommendation
raptor profile --counts counts.csv --output report.html

# Report shows top 3 options:
# 1. Pipeline 3 (Score: 145) - Recommended
# 2. Pipeline 1 (Score: 138) - High accuracy alternative
# 3. Pipeline 4 (Score: 125) - Fast alternative

# Benchmark only these three to validate
raptor compare --pipelines 1,3,4 --data fastq/ --output validation/

# Generate comparison report
raptor report --results validation/ --output validation_report.html
```

---

### Example 3: Different Priority Modes

```bash
# For publication (maximize accuracy)
raptor profile --counts data.csv --priority accuracy
# → Likely recommends: STAR-RSEM-DESeq2

# For large cohort study (maximize speed)
raptor profile --counts data.csv --priority speed
# → Likely recommends: Kallisto-Sleuth or Salmon-edgeR

# For cloud computing (minimize cost/memory)
raptor profile --counts data.csv --priority memory
# → Likely recommends: Kallisto-Sleuth
```

---

### Example 4: Automated Workflow (Python API)

```python
from raptor import RNAseqDataProfiler, PipelineRecommender
import pandas as pd

# Load your data
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Profile the data
profiler = RNAseqDataProfiler(counts, metadata)
profile = profiler.run_full_profile()

# Get recommendation
recommender = PipelineRecommender(profile)
recommendation = recommender.get_recommendation(priority='balanced')

# Print results
print(f"Recommended Pipeline: {recommendation['primary']['pipeline_name']}")
print(f"Score: {recommendation['primary']['score']}/200")
print(f"\nReasoning:")
for reason in recommendation['primary']['reasoning']:
    print(f"   {reason}")

# Export for documentation
recommendation.to_json('pipeline_recommendation.json')
```

---

### Example 5: Full Benchmarking Workflow

```bash
# Step 1: Simulate test data with known ground truth
raptor simulate --n-genes 2000 --n-samples 6 --n-de 400 --output simulated/

# Step 2: Run all 8 pipelines
raptor compare --data simulated/fastq/ --output benchmark_results/

# Step 3: Generate comprehensive report
raptor report --results benchmark_results/ --output full_benchmark.html

# Step 4: Extract specific metrics
raptor metrics --results benchmark_results/ --format csv --output metrics.csv
```

---

### Example 6: Custom Configuration

```bash
# Create custom config
cat > my_config.yaml << EOF
output_dir: "my_results/"
threads: 16
memory: "64GB"

# Analysis parameters
fdr_threshold: 0.01  # More stringent
min_counts: 20       # Higher filtering

# Profiling weights
weights:
  accuracy: 0.5
  speed: 0.3
  memory: 0.2
EOF

# Run with custom config
raptor profile --counts data.csv --config my_config.yaml --output report.html
```

---

## 📖 Documentation

### Complete Guides

- 📘 **[Profile & Recommend Guide](docs/PROFILE_RECOMMEND.md)** - Detailed documentation of recommendation system
- 📗 **[Benchmarking Guide](docs/BENCHMARKING.md)** - How to compare pipelines systematically
- 📙 **[Installation Guide](docs/INSTALLATION.md)** - Detailed setup for different systems
- 📕 **[Pipeline Details](docs/PIPELINES.md)** - In-depth description of each pipeline
- 📔 **[API Reference](docs/API.md)** - Python API documentation

### Quick References

- **[FAQ](docs/FAQ.md)** - Common questions answered
- **[Troubleshooting](docs/TROUBLESHOOTING.md)** - Solving common problems
- **[Best Practices](docs/BEST_PRACTICES.md)** - Tips for optimal results
- **[Metrics Explained](docs/METRICS.md)** - Understanding performance metrics

### Tutorials

- **[Quick Start Tutorial](docs/tutorials/QUICKSTART.md)** - 5-minute introduction
- **[First Analysis](docs/tutorials/FIRST_ANALYSIS.md)** - Step-by-step walkthrough
- **[Advanced Usage](docs/tutorials/ADVANCED.md)** - Complex scenarios
- **[Integration Guide](docs/tutorials/INTEGRATION.md)** - Incorporating into existing workflows

---

##  Citation

If you use RAPTOR in your research, please cite:

### Software Citation

```bibtex
@software{raptor2025,
  author = {Ayeh Bolouki},
  title = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/AyehBlk/RAPTOR},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

### For Profile & Recommend Feature

```bibtex
@software{raptor_profile2025,
  author = {Ayeh Bolouki},
  title = {RAPTOR Profile \& Recommend: Intelligent Pipeline Selection for RNA-seq Analysis},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/AyehBlk/RAPTOR}
}
```

### For Published Paper (Coming Soon)

```bibtex
@article{raptor2025paper,
  author = {Ayeh Bolouki},
  title = {RAPTOR: A Comprehensive Framework for RNA-seq Pipeline Benchmarking and Intelligent Selection},
  journal = {Bioinformatics},
  year = {2025},
  note = {In preparation}
}
```

---

## 🤝 Contributing

We welcome contributions from the bioinformatics community! RAPTOR thrives on collaborative input.

### Ways to Contribute

1. ** Add New Pipelines**
   - Implement additional analysis workflows
   - Contribute cutting-edge methods
   - Share your lab's optimized pipeline

2. ** Improve Documentation**
   - Write tutorials
   - Translate documentation
   - Share use cases and examples

3. ** Report Bugs**
   - Found an issue? Let us know!
   - Include reproducible examples
   - Help us improve stability

4. ** Suggest Features**
   - Ideas for improvements
   - New metrics or visualizations
   - Enhanced recommendation logic

5. ** Share Benchmark Results**
   - Contribute results from your data
   - Help improve recommendations
   - Build community knowledge

### How to Contribute

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

See **[CONTRIBUTING.md](CONTRIBUTING.md)** for detailed guidelines, code standards, and testing requirements.

---

##  Acknowledgments

RAPTOR builds upon the excellent work of the bioinformatics community:

### Alignment & Quantification Tools
- **STAR** - Dobin et al., Bioinformatics 2013
- **HISAT2** - Kim et al., Nature Methods 2015
- **Salmon** - Patro et al., Nature Methods 2017
- **Kallisto** - Bray et al., Nature Biotechnology 2016
- **RSEM** - Li & Dewey, BMC Bioinformatics 2011
- **StringTie** - Pertea et al., Nature Biotechnology 2015
- **HTSeq** - Anders et al., Bioinformatics 2015
- **featureCounts** - Liao et al., Bioinformatics 2014

### Statistical Analysis Tools
- **DESeq2** - Love et al., Genome Biology 2014
- **edgeR** - Robinson et al., Bioinformatics 2010
- **limma** - Ritchie et al., Nucleic Acids Research 2015
- **Sleuth** - Pimentel et al., Nature Methods 2017
- **Ballgown** - Frazee et al., Nature Biotechnology 2015
- **NOISeq** - Tarazona et al., Nucleic Acids Research 2015
- **EBSeq** - Leng et al., Bioinformatics 2013

### Infrastructure
- **Bioconductor** - Huber et al., Nature Methods 2015
- **NumPy & SciPy** - Harris et al., Nature 2020
- **pandas** - McKinney, Python for Data Analysis 2012
- **scikit-learn** - Pedregosa et al., JMLR 2011

**Special thanks** to the open-source bioinformatics community for making computational research accessible to everyone, everywhere.

---

##  Support & Contact

### Get Help

-  **Documentation**: Browse [docs/](docs/) folder
-  **Discussions**: [GitHub Discussions](https://github.com/AyehBlk/RAPTOR/discussions)
-  **Bug Reports**: [GitHub Issues](https://github.com/AyehBlk/RAPTOR/issues)
-  **Email**: ayehbolouki1988@gmail.com | ayehgeek@gmail.com

### Contact Author

**Ayeh Bolouki**
-  Postdoctoral Researcher
-  GIGA Center, University of Liège, Belgium
-  GitHub: [@AyehBlk](https://github.com/AyehBlk)
-  LinkedIn: [Ayeh Bolouki](https://linkedin.com/in/ayeh-bolouki)

---

##  License

RAPTOR is released under the **MIT License**.

```
MIT License

Copyright (c) 2025 Ayeh Bolouki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

See **[LICENSE](LICENSE)** file for complete text.

---

## 📈 Project Status

- ✅ **v2.0**: Current release with intelligent recommendation system
- 🚧 **v2.1**: Machine learning-enhanced recommendations (In progress)
- 🔜 **v3.0**: Single-cell RNA-seq support (Planned)

---

## 🌟 Star Us on GitHub!

If RAPTOR helps your research, please ⭐ **star this repository**!

It helps other researchers discover the tool and supports continued development.

[![GitHub stars](https://img.shields.io/github/stars/AyehBlk/RAPTOR?style=social)](https://github.com/AyehBlk/RAPTOR/stargazers)

---

## 🌍 Mission: Free Science for Everybody

RAPTOR is committed to **open science** and making computational tools accessible to all researchers, regardless of:

-  Geographic location
-  Institutional affiliation
-  Funding availability
-  Technical background

**"Let's make free science for everybody around the world."**

All tools are free and open-source.

---

## 🦖 Project Structure

```
RAPTOR/
├── README.md                          # This file
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation metadata
├── CONTRIBUTING.md                    # Contribution guidelines
├── environment.yml                    # Conda environment
├── setup.py                           # Python package setup
│
├── raptor/                            # Main Python package
│   ├── __init__.py
│   ├── cli.py                        # Command-line interface
│   ├── profiler.py                   # Data profiling
│   ├── recommender.py                # Pipeline recommendation
│   ├── benchmark.py                  # Benchmarking framework
│   ├── simulate.py                   # Data simulation
│   ├── report.py                     # Report generation
│   └── utils.py                      # Utility functions
│
├── pipelines/                         # Pipeline implementations
│   ├── pipeline1_star_rsem_deseq2/
│   ├── pipeline2_hisat2_stringtie_ballgown/
│   ├── pipeline3_salmon_edger/
│   ├── pipeline4_kallisto_sleuth/
│   ├── pipeline5_star_htseq_limma/
│   ├── pipeline6_star_featurecounts_noiseq/
│   ├── pipeline7_bowtie2_rsem_ebseq/
│   └── pipeline8_hisat2_cufflinks_cuffdiff/
│
├── scripts/                           # Analysis and utility scripts
│   ├── 00_simulate_data.R
│   ├── 01_run_all_pipelines.sh
│   ├── 02_profile_data.py
│   ├── 03_compare_results.R
│   └── 04_visualize_comparison.R
│
├── config/                            # Configuration files
│   ├── config.yaml                   # Default configuration
│   └── pipelines.yaml                # Pipeline parameters
│
├── docs/                              # Documentation
│   ├── INSTALLATION.md
│   ├── PROFILE_RECOMMEND.md
│   ├── BENCHMARKING.md
│   ├── PIPELINES.md
│   ├── API.md
│   ├── FAQ.md
│   ├── TROUBLESHOOTING.md
│   └── tutorials/
│
├── examples/                          # Example workflows
│   ├── demo.sh
│   ├── quick_profile.sh
│   └── full_benchmark.sh
│
└── tests/                             # Unit tests
    ├── test_profiler.py
    ├── test_recommender.py
    └── test_pipelines.sh
```

---

<div align="center">

## 🦖 Make confident, evidence-based decisions for your RNA-seq analysis

**⭐ Star • 🍴 Fork • 📢 Share**

Made with ❤️ for the open science community

[Get Started](#-quick-start) • [Documentation](docs/) • [Report Issue](https://github.com/AyehBlk/RAPTOR/issues)

</div>
