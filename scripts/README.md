# RAPTOR Analysis Scripts

This directory contains utility scripts for RNA-seq data simulation, pipeline execution, and results comparison.

## 📋 Scripts Overview

| Script | Purpose | Language |
|--------|---------|----------|
| `00_simulate_data.R` | Simulate RNA-seq count data for testing | R |
| `01_run_all_pipelines.sh` | Run all 8 pipelines on the same dataset | Bash |
| `02_profile_data.py` | Profile data characteristics for pipeline selection | Python |
| `03_compare_results.R` | Compare differential expression results | R |
| `04_visualize_comparison.R` | Create comparison visualizations | R |

## 🚀 Quick Start

### 1. Simulate Test Data

```bash
# Simulate RNA-seq data with 6 samples and 1000 DE genes
Rscript scripts/00_simulate_data.R \
  --output simulated_data \
  --nsamples 6 \
  --ndiff 1000 \
  --foldchange 3
```

**Output:**
- `simulated_data/reads/*.fastq.gz` - Simulated FASTQ files
- `simulated_data/transcriptome.fa` - Reference transcriptome
- `simulated_data/annotation.gtf` - Gene annotation
- `simulated_data/sample_metadata.csv` - Sample information
- `simulated_data/truth_set.csv` - Ground truth DE genes

### 2. Profile Your Data

```bash
# Profile count matrix to get data characteristics
python scripts/02_profile_data.py \
  data/count_matrix.csv \
  --output data_profile.json
```

**Output:**
- Data statistics (library size, zero percentage, etc.)
- Pipeline recommendations based on data characteristics
- JSON profile for programmatic use

### 3. Run All Pipelines

```bash
# Run all 8 pipelines on your data
bash scripts/01_run_all_pipelines.sh \
  simulated_data/reads \
  results/benchmark \
  simulated_data
```

**Arguments:**
- `input_dir` - Directory with FASTQ files
- `output_dir` - Base output directory
- `reference_dir` - Directory with genome/transcriptome

**Output:**
- `results/benchmark/pipeline1_star_rsem_deseq2/` - Pipeline 1 results
- `results/benchmark/pipeline2_hisat2_stringtie_ballgown/` - Pipeline 2 results
- ... (all 8 pipelines)
- `results/benchmark/comparison_summary.txt` - Summary report

### 4. Compare Results

```bash
# Compare results from all pipelines
Rscript scripts/03_compare_results.R \
  results/benchmark \
  --truth simulated_data/truth_set.csv
```

**Output:**
- `comparison_summary.csv` - Number of DEGs per pipeline
- `overlap_matrix.csv` - Pairwise overlaps between pipelines
- `performance_metrics.csv` - Precision, recall, F1 scores (if truth provided)

### 5. Visualize Comparison

```bash
# Create comparison plots
Rscript scripts/04_visualize_comparison.R \
  results/benchmark
```

**Output:**
- `plots/01_significant_genes_barplot.pdf` - DEG counts per pipeline
- `plots/02_overlap_heatmap.pdf` - Pipeline similarity heatmap
- `plots/03_venn_diagram.pdf` - Overlap visualization
- `plots/04a_precision_recall.pdf` - Performance comparison
- `plots/04b_f1_scores.pdf` - F1 score comparison

## 📖 Detailed Documentation

### Script 1: 00_simulate_data.R

Simulate RNA-seq count data using the polyester package.

**Usage:**
```bash
Rscript scripts/00_simulate_data.R [options]
```

**Options:**
- `-o, --output DIR` - Output directory (default: simulated_data)
- `-n, --nsamples N` - Number of samples, must be even (default: 6)
- `-g, --ngenes N` - Number of genes (default: 10000)
- `-r, --reads N` - Reads per sample (default: 1000000)
- `-d, --ndiff N` - Number of DE genes (default: 1000)
- `-f, --foldchange N` - Fold change for DE genes (default: 3)
- `-s, --seed N` - Random seed (default: 42)

**Example:**
```bash
# Small test dataset
Rscript scripts/00_simulate_data.R \
  --output test_data \
  --nsamples 4 \
  --ngenes 5000 \
  --ndiff 500

# Large realistic dataset
Rscript scripts/00_simulate_data.R \
  --output large_data \
  --nsamples 12 \
  --ngenes 20000 \
  --ndiff 2000 \
  --reads 5000000
```

**Output Structure:**
```
output_dir/
├── reads/
│   ├── control_1_R1.fastq.gz
│   ├── control_1_R2.fastq.gz
│   └── ...
├── transcriptome.fa
├── annotation.gtf
├── sample_metadata.csv
├── truth_set.csv
└── simulation_summary.txt
```

### Script 2: 02_profile_data.py

Analyze RNA-seq count data characteristics and recommend suitable pipelines.

**Usage:**
```bash
python scripts/02_profile_data.py <count_matrix.csv> [options]
```

**Options:**
- `-o, --output FILE` - Output JSON file (default: data_profile.json)
- `-v, --verbose` - Verbose output

**Input Format:**
CSV file with genes as rows and samples as columns:
```csv
gene_id,sample1,sample2,sample3
GENE001,100,120,95
GENE002,50,55,48
...
```

**Profile Metrics:**
- **Sample size**: Number of samples
- **Library size**: Mean, median, CV
- **Gene detection**: Genes detected per sample
- **Zero inflation**: Percentage of zero counts
- **Expression characteristics**: Mean, variance, dynamic range
- **Dispersion**: Biological coefficient of variation
- **Sample correlation**: Inter-sample similarity
- **Batch effects**: PCA-based outlier detection

**Example:**
```bash
python scripts/02_profile_data.py \
  data/counts.csv \
  --output profile.json \
  --verbose
```

**Output:**
```json
{
  "n_samples": 6,
  "n_genes": 10000,
  "mean_library_size": 1234567,
  "library_size_cv": 0.15,
  "zero_percentage": 45.2,
  "mean_sample_correlation": 0.92,
  "recommendations": [
    "Pipeline 1 (STAR-RSEM-DESeq2) - gold standard for adequate samples",
    "Pipeline 5 (limma-voom) - excellent for complex designs"
  ]
}
```

### Script 3: 01_run_all_pipelines.sh

Execute all 8 RAPTOR pipelines on the same dataset for benchmarking.

**Usage:**
```bash
bash scripts/01_run_all_pipelines.sh <input_dir> <output_dir> <reference_dir>
```

**Arguments:**
1. `input_dir` - Directory containing FASTQ files and sample sheet
2. `output_dir` - Base directory for all pipeline outputs
3. `reference_dir` - Directory with reference genome/transcriptome

**Reference Directory Structure:**
```
reference_dir/
├── genome.fa           # For alignment-based pipelines
├── annotation.gtf      # Gene annotation
└── transcriptome.fa    # For alignment-free pipelines
```

**Features:**
- Automatically updates pipeline configurations with correct paths
- Runs pipelines in sequence
- Tracks successes and failures
- Generates comparison summary
- Creates detailed log file

**Example:**
```bash
bash scripts/01_run_all_pipelines.sh \
  data/fastq \
  results/comparison \
  references/human_grch38
```

**Output:**
```
results/comparison/
├── pipeline1_star_rsem_deseq2/
│   ├── qc/
│   ├── results/
│   └── logs/
├── pipeline2_hisat2_stringtie_ballgown/
├── ... (all 8 pipelines)
├── comparison_summary.txt
└── pipeline_comparison_YYYYMMDD_HHMMSS.log
```

### Script 4: 03_compare_results.R

Compare differential expression results across pipelines.

**Usage:**
```bash
Rscript scripts/03_compare_results.R <results_dir> [--truth truth_set.csv]
```

**Arguments:**
- `results_dir` - Directory containing pipeline output folders
- `--truth FILE` - (Optional) Ground truth file for accuracy assessment

**Features:**
- Finds and loads DEG results from all pipelines
- Standardizes column names across different formats
- Calculates pairwise overlaps (Jaccard index)
- Identifies consensus genes
- Computes performance metrics if truth set provided

**Output:**
- `comparison_summary.csv` - DEG counts per pipeline
- `overlap_matrix.csv` - Pairwise overlap matrix
- `performance_metrics.csv` - Precision, recall, F1 scores

**Performance Metrics:**
- **TP (True Positives)**: Correctly identified DE genes
- **FP (False Positives)**: Incorrectly identified as DE
- **FN (False Negatives)**: Missed DE genes
- **Precision**: TP / (TP + FP)
- **Recall**: TP / (TP + FN)
- **F1 Score**: Harmonic mean of precision and recall

**Example:**
```bash
# Without ground truth
Rscript scripts/03_compare_results.R results/benchmark

# With ground truth
Rscript scripts/03_compare_results.R \
  results/benchmark \
  --truth simulated_data/truth_set.csv
```

### Script 5: 04_visualize_comparison.R

Create comprehensive visualizations of pipeline comparison results.

**Usage:**
```bash
Rscript scripts/04_visualize_comparison.R <results_dir>
```

**Prerequisites:**
Must run `03_compare_results.R` first to generate comparison files.

**Generated Plots:**

1. **Bar Plot** (`01_significant_genes_barplot.pdf`)
   - Number of significant genes per pipeline
   - Shows pipeline performance at a glance

2. **Overlap Heatmap** (`02_overlap_heatmap.pdf`)
   - Jaccard similarity between all pipeline pairs
   - Identifies which pipelines produce similar results

3. **Venn Diagram** (`03_venn_diagram.pdf`)
   - Overlap visualization (for 2-5 pipelines)
   - Shows unique and shared DE genes

4. **Precision-Recall Plot** (`04a_precision_recall.pdf`)
   - Performance comparison against ground truth
   - Higher and to the right is better

5. **F1 Score Comparison** (`04b_f1_scores.pdf`)
   - Overall performance ranking
   - Balances precision and recall

6. **Confusion Components** (`04c_confusion_components.pdf`)
   - TP, FP, FN breakdown per pipeline
   - Shows error types

7. **Summary Comparison** (`05_summary_comparison.pdf`)
   - Multi-metric radar plot
   - Overall performance profile

**Example:**
```bash
Rscript scripts/04_visualize_comparison.R results/benchmark
```

**Output:**
```
results/benchmark/plots/
├── 01_significant_genes_barplot.pdf
├── 02_overlap_heatmap.pdf
├── 03_venn_diagram.pdf
├── 04a_precision_recall.pdf
├── 04b_f1_scores.pdf
├── 04c_confusion_components.pdf
└── 05_summary_comparison.pdf
```

## 🔄 Complete Workflow

Here's a complete workflow from data simulation to visualization:

```bash
# Step 1: Simulate test data
Rscript scripts/00_simulate_data.R \
  --output sim_data \
  --nsamples 6 \
  --ndiff 1000

# Step 2: Profile the simulated data (optional)
python scripts/02_profile_data.py \
  sim_data/counts.csv \
  --output sim_data/profile.json

# Step 3: Run all pipelines
bash scripts/01_run_all_pipelines.sh \
  sim_data/reads \
  results/test_run \
  sim_data

# Step 4: Compare results
Rscript scripts/03_compare_results.R \
  results/test_run \
  --truth sim_data/truth_set.csv

# Step 5: Visualize comparison
Rscript scripts/04_visualize_comparison.R \
  results/test_run
```

## 📊 Understanding Results

### Interpreting Performance Metrics

**Precision (Positive Predictive Value)**:
- What percentage of identified DEGs are true DE genes?
- High precision = Few false positives
- Important when you want confidence in results

**Recall (Sensitivity)**:
- What percentage of true DE genes are identified?
- High recall = Few false negatives
- Important when you don't want to miss anything

**F1 Score**:
- Harmonic mean of precision and recall
- Balances both metrics
- Use for overall performance ranking

**Jaccard Index**:
- Measures similarity between two sets
- Range: 0 (no overlap) to 1 (identical)
- Values > 0.5 indicate good agreement

### Choosing Best Pipeline

Consider multiple factors:

1. **Performance Metrics** (if truth available):
   - Highest F1 score
   - Good balance of precision and recall

2. **Consensus**:
   - Pipelines with high Jaccard similarity
   - Genes found by multiple pipelines (high confidence)

3. **Your Priorities**:
   - Need speed? → Salmon or Kallisto
   - Need accuracy? → STAR-RSEM-DESeq2
   - Complex design? → limma-voom
   - Small sample? → EBSeq or NOISeq

## 🔧 Dependencies

### R Packages
```r
# Data manipulation
install.packages("dplyr")
install.packages("tidyr")
install.packages("purrr")

# Visualization
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("ggVennDiagram")

# Bioconductor
BiocManager::install("polyester")
BiocManager::install("Biostrings")
BiocManager::install("rtracklayer")
```

### Python Packages
```bash
pip install pandas numpy scipy scikit-learn
```

## 💡 Tips and Best Practices

### For Simulation

- Start with smaller datasets for testing
- Use `--seed` for reproducibility
- Match read depth to your real data
- Adjust `--foldchange` based on expected effect size

### For Pipeline Comparison

- Use same parameters across pipelines for fair comparison
- Ensure adequate computational resources
- Monitor log files for errors
- Keep intermediate files for troubleshooting

### For Analysis

- Always check QC metrics first
- Look for consensus among pipelines
- Don't rely on single pipeline
- Consider biological context, not just statistics

## 📧 Support

- GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
- Email: ayehbolouki1988@gmail.com

## 📜 License

MIT License - See LICENSE file

---

**Created for RAPTOR Project**
**Author**: Ayeh Bolouki
**Organization**: GIGA Center, University of Liège
