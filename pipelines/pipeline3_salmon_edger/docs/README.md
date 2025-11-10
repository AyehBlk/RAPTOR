# Salmon-edgeR

## Overview

**Fast alignment-free with robust dispersion estimation**

### Pipeline Components

- **Alignment**: None (alignment-free)
- **Quantification**: Salmon
- **Statistics**: edgeR

## Quick Start

### 1. Configure

Edit `config/pipeline_config.yaml`:

```yaml
io:
  input_dir: "/path/to/fastq"
  sample_sheet: "/path/to/samples.csv"
  output_dir: "/path/to/output"
  
resources:
  threads: 8
```

### 2. Run

```bash
bash scripts/run_pipeline.sh config/pipeline_config.yaml
```

### 3. Results

```
output/
├── qc/                    # Quality control reports
├── results/               # Differential expression results
└── logs/                  # Pipeline logs
```

## Documentation

See main [RAPTOR documentation](../../README.md) for detailed usage.

## Support

- GitHub Issues: https://github.com/AyehBlk/RAPTOR/issues
- Email: ayehbolouki1988@gmail.com
