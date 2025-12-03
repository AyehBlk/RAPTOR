# 🦖 RAPTOR Tests

Test scripts for RAPTOR v2.1.0.

## Test Scripts

| Script | Type | Description |
|--------|------|-------------|
| `test_ml_system.py` | Unit | Tests ML recommender components (imports, feature extraction, training, prediction) |
| `test_recommender.py` | Unit | Pytest unit tests for PipelineRecommender class |
| `test_raptor_v2_1_0.py` | Integration | Comprehensive v2.1.0 tests (modules, version, dependencies, CLI) |
| `test_pipelines.sh` | Integration | Pipeline installation and tool availability tests |

## Running Tests

```bash
# Run all Python tests
python test_ml_system.py
python test_raptor_v2_1_0.py

# Run pytest tests
pytest test_recommender.py -v

# Run shell tests
chmod +x test_pipelines.sh
./test_pipelines.sh
```

## Test Coverage

- **Module imports:** All v2.0.0 and v2.1.0 modules
- **ML system:** Feature extraction, data generation, model training, prediction
- **Recommender:** Initialization, recommendation generation, structure validation
- **Dependencies:** Required and optional packages
- **Bioinformatics tools:** STAR, Salmon, Kallisto, HISAT2, etc.

## Requirements

```bash
pip install raptor-rnaseq[ml]
pip install pytest
```

---
*RAPTOR v2.1.0 - Making free science for everybody around the world 🌍*
