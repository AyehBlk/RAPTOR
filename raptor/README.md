# RAPTOR Examples

Comprehensive examples demonstrating RAPTOR v2.1.0 workflows and new features.

## 🆕 What's New in v2.1.0

- **ML-Based Pipeline Recommendations**: Machine learning models trained on synthetic benchmarks
- **Interactive Dashboard**: Real-time visualization and monitoring
- **Advanced Quality Assessment**: Batch effect detection and comprehensive quality scoring
- **Real-Time Resource Monitoring**: Track CPU, memory, and disk usage during analysis
- **Ensemble Analysis Methods**: Multiple statistical approaches for robust results

## Quick Start

```bash
# Run comprehensive demo (includes v2.1.0 features)
chmod +x demo.sh && ./demo.sh

# Profile your data with ML recommendations
chmod +x quick_profile.sh && ./quick_profile.sh my_counts.csv

# Full benchmark comparison
chmod +x full_benchmark.sh && ./full_benchmark.sh -d fastq_data/

# ML workflow example
python example_ml_workflow.py --n-datasets 200

# Quality assessment example
python example_quality_assessment.py
```

## Available Examples

### Shell Scripts

#### 1. `demo.sh` - Comprehensive Demo
**Duration:** ~30 minutes  
**Description:** Complete demonstration of RAPTOR's main features

**Features demonstrated:**
- Data simulation
- Data profiling with ML recommendations
- Interactive dashboard usage
- Python API examples
- Quality assessment

**Usage:**
```bash
./demo.sh
```

**Output:**
- Simulated RNA-seq data
- Profiling reports with ML recommendations
- HTML interactive reports
- Python API examples

---

#### 2. `quick_profile.sh` - Fast Profiling Workflow
**Duration:** ~5 minutes  
**Description:** Rapid data profiling and pipeline recommendations

**Features:**
- Data validation
- Quick quality check
- ML-based recommendations
- Resource estimation

**Usage:**
```bash
# Basic profiling
./quick_profile.sh my_counts.csv

# With metadata
./quick_profile.sh my_counts.csv my_metadata.csv

# Custom output
./quick_profile.sh counts.csv -o my_results/ -t 8
```

**Options:**
- `-o, --output DIR`: Output directory (default: `quick_profile_results/`)
- `-t, --threads N`: Number of threads (default: 4)
- `-h, --help`: Show help message

---

#### 3. `full_benchmark.sh` - Comprehensive Pipeline Comparison
**Duration:** 2-24 hours (depends on mode)  
**Description:** Rigorous benchmarking of multiple RNA-seq pipelines

**Features:**
- Compare up to 8 different pipelines
- Runtime and memory tracking
- Accuracy evaluation with ground truth
- Parallel pipeline execution
- Comprehensive HTML reports

**Usage:**
```bash
# Quick benchmark (2-3 hours)
./full_benchmark.sh -d data/ --mode quick

# Standard benchmark (4-6 hours)
./full_benchmark.sh -d data/ --mode standard -t 16

# Full benchmark with all pipelines (12-24 hours)
./full_benchmark.sh -d data/ --mode full -t 32 -m 64G

# With ground truth for accuracy
./full_benchmark.sh -d data/ --ground-truth truth.csv

# Parallel execution
./full_benchmark.sh -d data/ -p 1,3,4,5 --parallel 2
```

**Benchmark Modes:**
- `quick`: 2-3 fastest pipelines (~1-2 hours)
- `standard`: 4 recommended pipelines (~4-6 hours)
- `full`: All 8 pipelines (~12-24 hours)

**Pipeline Options:**
1. STAR-RSEM-DESeq2 (highest accuracy)
2. HISAT2-StringTie-Ballgown (novel transcripts)
3. Salmon-edgeR (recommended, balanced)
4. Kallisto-Sleuth (fastest)
5. STAR-HTSeq-limma (complex designs)
6. STAR-featureCounts-NOISeq (small samples)
7. Bowtie2-RSEM-EBSeq (memory efficient)
8. HISAT2-Cufflinks-Cuffdiff (legacy)

**Options:**
- `-d, --data DIR`: FASTQ data directory (required)
- `-o, --output DIR`: Output directory (default: `benchmark_results/`)
- `-p, --pipelines LIST`: Pipeline IDs (default: 1,3,4,5)
- `-t, --threads N`: Number of threads (default: 8)
- `-m, --memory SIZE`: Memory limit (default: 32G)
- `--mode MODE`: Benchmark mode (quick|standard|full)
- `--parallel N`: Run N pipelines in parallel (default: 1)
- `--ground-truth FILE`: Ground truth for accuracy evaluation

---

### Python Examples

#### 4. `example_ml_workflow.py` - ML Recommender Workflow
**Duration:** ~15-30 minutes  
**Description:** Complete machine learning pipeline recommendation workflow

**Features demonstrated:**
- Synthetic training data generation
- ML model training (Random Forest, Gradient Boosting)
- Model evaluation and visualization
- Feature importance analysis
- Predictions on new data
- Model comparison

**Usage:**
```bash
# Basic workflow with 200 datasets
python example_ml_workflow.py

# With more training data
python example_ml_workflow.py --n-datasets 500

# Using existing data
python example_ml_workflow.py --skip-generation --data-dir my_data/

# Compare models
python example_ml_workflow.py --compare-models

# Use Gradient Boosting
python example_ml_workflow.py --model-type gradient_boosting
```

**Options:**
- `--n-datasets N`: Number of synthetic datasets (default: 200)
- `--model-type`: Model type (random_forest|gradient_boosting)
- `--skip-generation`: Skip data generation, use existing
- `--compare-models`: Compare Random Forest vs Gradient Boosting
- `--data-dir DIR`: Training data directory
- `--model-dir DIR`: Model output directory

**Output:**
- Trained ML models
- Performance metrics
- Confusion matrices
- Feature importance plots
- Model comparison charts

---

#### 5. `example_quality_assessment.py` - Advanced Quality Assessment
**Duration:** ~5 minutes  
**Description:** Comprehensive data quality assessment with visualizations

**Features demonstrated:**
- Basic quality assessment
- Batch effect detection
- Outlier identification
- Poor quality data detection
- Multi-dataset comparison

**Usage:**
```bash
# Run all examples
python example_quality_assessment.py
```

**Examples included:**
1. **Basic Usage**: Quick quality check with visualization
2. **Batch Detection**: Identify and quantify batch effects
3. **Poor Quality**: Detect problematic datasets
4. **Dataset Comparison**: Compare quality across multiple datasets

**Output:**
- Quality assessment reports
- Batch effect visualizations
- Outlier detection plots
- Comparison charts

---

## Requirements

### Shell Scripts
- RAPTOR installed: `pip install raptor-rnaseq`
- For `full_benchmark.sh`:
  - All bioinformatics tools (STAR, Salmon, Kallisto, etc.)
  - Reference genome indices
  - FASTQ data files

### Python Scripts
- Python 3.8+
- RAPTOR with dependencies:
  ```bash
  pip install raptor-rnaseq[ml]
  ```
- Required packages: pandas, numpy, matplotlib, seaborn, scikit-learn

---

## File Organization

```
examples/
├── README.md                          # This file
├── demo.sh                           # Comprehensive demo
├── quick_profile.sh                  # Fast profiling
├── full_benchmark.sh                 # Pipeline comparison
├── example_ml_workflow.py            # ML recommender workflow
├── example_quality_assessment.py     # Quality assessment
└── outputs/                          # Generated results (gitignored)
```

---

## Expected Outputs

### From Shell Scripts
- **demo.sh**:
  - `raptor_demo_YYYYMMDD_HHMMSS/`
    - Simulated data
    - Profile results
    - Python API examples
    
- **quick_profile.sh**:
  - `quick_profile_results/`
    - `raptor_profile_report.html`
    - `recommendations.json`
    - `QUICK_SUMMARY.txt`
    
- **full_benchmark.sh**:
  - `benchmark_results/`
    - `benchmark_comparison.html`
    - `benchmark_results.json`
    - `BENCHMARK_SUMMARY.txt`
    - Pipeline-specific outputs

### From Python Scripts
- **example_ml_workflow.py**:
  - `ml_training_data/`: Synthetic datasets
  - `models/`: Trained ML models
  - `figures/`: Performance visualizations
  
- **example_quality_assessment.py**:
  - `example1_quality.png`
  - `example2_quality_with_batch.png`
  - `example3_poor_quality.png`
  - `example4_comparison.png`

---

## Troubleshooting

### Common Issues

**1. RAPTOR not found**
```bash
pip install raptor-rnaseq
# or with ML features
pip install raptor-rnaseq[ml]
```

**2. Missing bioinformatics tools (for benchmarking)**
- Check tool installation: `raptor check --tools`
- Install missing tools or use Docker/Singularity containers

**3. Memory errors**
- Reduce number of threads: `-t 4`
- Reduce memory limit appropriately
- Use smaller datasets for testing

**4. Python import errors**
```bash
# Ensure RAPTOR is properly installed
pip install --upgrade raptor-rnaseq[ml]

# Check installation
python -c "import raptor; print(raptor.__version__)"
```

**5. Permission denied**
```bash
chmod +x *.sh
```

---

## Best Practices

### For Shell Scripts
1. **Start small**: Use `demo.sh` before `full_benchmark.sh`
2. **Test first**: Use `--mode quick` for benchmarking
3. **Check resources**: Ensure sufficient CPU/memory
4. **Use metadata**: Provide sample information when available
5. **Parallel execution**: Use `--parallel` for faster benchmarks

### For Python Scripts
1. **Virtual environment**: Use isolated Python environment
2. **Dependencies**: Install all requirements first
3. **Save outputs**: Preserve generated models and figures
4. **Experiment**: Modify parameters to explore features
5. **Documentation**: Read inline comments for details

---

## Integration with RAPTOR CLI

All example workflows can be integrated into standard RAPTOR usage:

```bash
# 1. Profile your data (uses ML recommendations)
raptor profile --counts counts.csv --metadata metadata.csv

# 2. Get recommendations
raptor recommend --profile profile.json

# 3. Run selected pipeline
raptor run --pipeline 3 --data fastq_data/

# 4. Compare pipelines
raptor compare --data fastq_data/ --pipelines 1,3,4,5

# 5. Generate report
raptor report --results results/ --output report.html
```

---

## Additional Resources

### Documentation
- [RAPTOR Documentation](https://github.com/AyehBlk/RAPTOR/tree/main/docs)
- [Tutorials](https://github.com/AyehBlk/RAPTOR/tree/main/docs/tutorials)
- [API Reference](https://github.com/AyehBlk/RAPTOR/tree/main/docs/api)

### Example Datasets
- Small test dataset: Use `demo.sh` to generate
- Real datasets: Check `docs/datasets.md`

### Videos & Tutorials
- Coming soon: Video walkthroughs
- Interactive Jupyter notebooks

---

## Support

**Questions or Issues?**
- GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
- Email: ayehbolouki1988@gmail.com
- Documentation: https://github.com/AyehBlk/RAPTOR/tree/main/docs

---

## Citation

If you use RAPTOR in your research, please cite:

```bibtex
@software{raptor2024,
  author = {Bolouki, Ayeh},
  title = {RAPTOR: RNA-seq Analysis Pipeline Testing and Optimization Resource},
  year = {2024},
  version = {2.1.0},
  url = {https://github.com/AyehBlk/RAPTOR}
}
```

---

## License

MIT License - See LICENSE file for details

---

**Last Updated:** December 2024  
**RAPTOR Version:** 2.1.0  
**Author:** Ayeh Bolouki  
