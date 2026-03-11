# STAR + Salmon Pipeline v2.2.0 - Complete Upgrade Summary

## 🎯 Overview

The **STAR + Salmon** pipeline has been fully upgraded to v2.2.0 with comprehensive validation and error handling. **This is the 6th pipeline complete (86% progress)!**

This is a **HYBRID** pipeline combining:
- **STAR alignment** → Genome BAM + Transcriptome BAM
- **Salmon quantification** → Gene counts + TPM + Bootstraps

🔥 **Best of both worlds:** BAM files for visualization + Bootstrap uncertainty for DE!

---

## 📊 File Summary

| File | Size | Lines | Purpose |
|------|------|-------|---------|
| **pipeline.py** | 26KB | 713 | Hybrid STAR + Salmon pipeline |
| **__init__.py** | 3.6KB | - | CLI exports + metadata |
| **config.yaml** | 8.3KB | 262 | Configuration |

**Total:** 3 files ready to use

---

## ⚠️  CRITICAL: Two-Index Requirement!

### This Pipeline is DIFFERENT!

**STAR + Salmon requires TWO indexes:**

```bash
# Index 1: STAR genome index
--index /path/to/star_index/

# Index 2: Salmon transcriptome index (REQUIRED!)
--salmon-index /path/to/salmon_index/
```

**Why two indexes?**
1. **STAR** aligns to the **genome** → produces genome BAM + transcriptome BAM
2. **Salmon** quantifies from **transcriptome** BAM → produces counts + TPM

**If you forget salmon_index:**
```python
raise ValidationError(
    f"Salmon index not found: {salmon_index}\n"
    f"  This hybrid pipeline requires BOTH:\n"
    f"    1. STAR genome index (provided with --index)\n"
    f"    2. Salmon transcriptome index (provide with --salmon-index)\n"
    f"  \n"
    f"  Build Salmon index with:\n"
    f"    salmon index -t transcripts.fa -i salmon_index"
)
```

---

## 🔄 What Changed

### 1. **pipeline.py** - Hybrid Logic ✅

#### Key Additions:
```python
# NEW: Comprehensive validation imports
from raptor.utils.validation import (
    validate_positive_integer,
    validate_directory_path,
)
from raptor.utils.errors import (
    handle_errors,
    ValidationError,
    PipelineError,
)

# NEW: Error handling decorators
@handle_errors(exit_on_error=False)
def __init__(self, ...):
def run_sample(self, ...):
def combine_counts(self, ...):
```

#### Hybrid-Specific Validations (15+ parameters):
- ✅ **salmon_index** - CRITICAL! Must exist
- ✅ **STAR two_pass_mode** - Must be 'None' or 'Basic'
- ✅ **Salmon library_type** - Must be one of 13 valid types
- ✅ **bootstraps** - Must be >= 0
- ✅ **STAR numeric parameters** - Positive integers
- ✅ **Two-stage validation** - Both STAR and Salmon must succeed

---

## 🎯 Two-Stage Process

### Stage 1: STAR Alignment
```python
# Run STAR alignment
cmd_star = [
    'STAR',
    '--genomeDir', str(index_path),
    '--quantMode', 'TranscriptomeSAM',  # CRITICAL: For Salmon!
    ...
]

# Outputs:
#   - Aligned.sortedByCoord.out.bam (genome BAM)
#   - Aligned.toTranscriptome.out.bam (for Salmon!)
```

### Stage 2: Salmon Quantification
```python
# Run Salmon on transcriptome BAM
cmd_salmon = [
    'salmon', 'quant',
    '-i', str(self.salmon_index),
    '-a', str(tx_bam),  # Read from transcriptome BAM!
    '--numBootstraps', str(self.bootstraps),
    ...
]

# Outputs:
#   - quant.sf (gene counts + TPM)
#   - Bootstraps for uncertainty
```

---

## 🔍 Critical Validations Implemented

### 1. Salmon Index Validation (MOST CRITICAL!)
```python
try:
    self.salmon_index = validate_directory_path(
        salmon_index,
        must_exist=True,
        description="Salmon transcriptome index"
    )
except Exception as e:
    raise ValidationError(
        f"Salmon index not found: {salmon_index}\n"
        f"  This hybrid pipeline requires BOTH:\n"
        f"    1. STAR genome index (provided with --index)\n"
        f"    2. Salmon transcriptome index (provide with --salmon-index)\n"
        f"  \n"
        f"  Build Salmon index with:\n"
        f"    salmon index -t transcripts.fa -i salmon_index"
    )
```

### 2. STAR Two-Pass Mode
```python
valid_two_pass = ['None', 'Basic']
if two_pass_mode not in valid_two_pass:
    raise ValidationError(
        f"Invalid two_pass_mode: '{two_pass_mode}'\n"
        f"  Must be one of: {', '.join(valid_two_pass)}\n"
        f"  Recommended: 'Basic' for better sensitivity"
    )
```

### 3. Salmon Library Type (Reused from Salmon pipeline!)
```python
valid_library_types = [
    'A',      # Auto-detect (recommended)
    'IU',     # Inward, unstranded
    'ISF',    # Inward, stranded forward
    'ISR',    # Inward, stranded reverse (dUTP - most common)
    # ... 13 total types
]

if library_type not in valid_library_types:
    raise ValidationError(
        f"Invalid library_type: '{library_type}'\n"
        f"  Must be one of: {', '.join(valid_library_types)}\n"
        f"  Common options:\n"
        f"    'A'   = auto-detect (recommended)\n"
        f"    'ISR' = inward, stranded reverse (dUTP - most common)"
    )
```

### 4. Bootstraps Validation
```python
if bootstraps < 0:
    raise ValidationError(
        f"bootstraps must be >= 0, got {bootstraps}\n"
        f"  Use 0 for no bootstraps (fastest)\n"
        f"  Use 30+ for uncertainty estimation\n"
        f"  Use 100+ for differential expression (recommended)"
    )
```

### 5. Transcriptome BAM Validation
```python
tx_bam = sample_star / "Aligned.toTranscriptome.out.bam"
if not tx_bam.exists():
    raise PipelineError(
        f"Transcriptome BAM not found for {sample.sample_id}\n"
        f"  Expected: {tx_bam}\n"
        f"  STAR must output transcriptome BAM for Salmon.\n"
        f"  Check: --quantMode TranscriptomeSAM was used."
    )
```

---

## 📋 Progress: 6/7 Pipelines Complete!

### ✅ **Completed** (86%):
1. **HISAT2 + featureCounts** - Low memory alternative
2. **Kallisto** - Fastest
3. **Salmon** - Best balance
4. **STAR + featureCounts** - Gold standard for genes
5. **STAR + RSEM** - Gold standard for isoforms
6. **STAR + Salmon** - Hybrid (BAM + bootstraps) ⭐ NEW!

### ⏳ **Remaining** (14%):
7. Quick pipelines (Module 1) - Simplified versions

---

## 🌟 Six-Way Comparison

| Feature | HISAT2 | Kallisto | Salmon | STAR-FC | STAR-RSEM | STAR-Salmon |
|---------|--------|----------|--------|---------|-----------|-------------|
| **Speed** | Moderate | Fastest | Fast | Slow | Slowest | Slow |
| **Memory** | 16GB | 4GB | 8GB | 32GB | 32GB | 32GB |
| **BAM Files** | Yes | No | No | Yes | Yes | ✅ Yes |
| **Bootstraps** | No | Yes | Yes | No | No | ✅ Yes |
| **Isoforms** | No | No | No | No | Yes | No |
| **Use Case** | Need BAM | Quick | Std DE | Gene DE | Isoform | BAM+Boot |

**STAR + Salmon Unique:** ONLY pipeline with BOTH BAM files AND bootstraps!

---

## ✨ What Makes STAR + Salmon Special

### Best of Both Worlds:
- ✅ **BAM files** (from STAR) - For IGV visualization
- ✅ **Bootstrap uncertainty** (from Salmon) - For sleuth DE
- ✅ **Bias correction** - GC + sequence bias
- ✅ **Two-pass mode** - Novel junction discovery
- ✅ **Widely compatible** - Works with standard tools

### When to Use:
- Need BAM files AND uncertainty estimates
- sleuth differential expression
- IGV visualization + uncertainty quantification
- Publication with both genome browser figures and bootstrap DE
- Want flexibility (both alignment and quantification)

### When NOT to Use:
- Only need gene counts (use STAR + featureCounts)
- Only need uncertainty (use Salmon alone)
- Don't need BAM files (use Salmon)
- Limited memory (<32GB)
- Want simplest workflow (use Salmon)

---

## 🔧 Enhanced Error Messages

### STAR Failure:
```python
raise PipelineError(
    f"STAR alignment failed for {sample.sample_id}\n"
    f"  Check log: {star_log}\n"
    f"  Common issues:\n"
    f"    - Index version mismatch\n"
    f"    - Corrupted FASTQ files\n"
    f"    - Insufficient memory ({self.memory_gb}GB allocated, need 32GB+)\n"
    f"    - Disk space full"
)
```

### Salmon Failure:
```python
raise PipelineError(
    f"Salmon quantification failed for {sample.sample_id}\n"
    f"  Check log: {salmon_log}\n"
    f"  Common issues:\n"
    f"    - Salmon index doesn't match transcripts in STAR index\n"
    f"    - Wrong library type (try -l A for auto-detect)\n"
    f"    - Corrupted transcriptome BAM\n"
    f"    - Salmon version mismatch with index"
)
```

---

## 📚 CLI Compatibility

### CLI Parameters Defined:
```python
STAR_SALMON_CLI_PARAMS = {
    'salmon-index': {
        'flag': '--salmon-index',
        'type': str,
        'required': True,  # ⚠️  REQUIRED!
        'help': 'Salmon transcriptome index (REQUIRED!)',
        'param_name': 'salmon_index'
    },
    'two-pass-mode': {
        'flag': '--two-pass-mode',
        'type': str,
        'default': 'Basic',
        'choices': ['None', 'Basic'],
        'help': 'STAR two-pass mode',
        'param_name': 'two_pass_mode'
    },
    # ... 5 total CLI parameters
}
```

### MODULE_INFO:
```python
MODULE_INFO = {
    'id': 'M7',
    'name': 'star_salmon',
    'stage': 2,
    'description': 'STAR + Salmon (BAM + bootstraps)',
    'version': '2.2.0',
    'requires': ['STAR>=2.7.0', 'Salmon>=1.9.0', 'samtools>=1.10'],
    'special_requirements': 'TWO indexes (STAR genome + Salmon transcriptome)',
}
```

---

## 🎯 Config File Highlights

### Two-Index Configuration:
```yaml
# ⚠️  REQUIRES TWO INDEXES:
#     1. STAR genome index (--index)
#     2. Salmon transcriptome index (--salmon-index)

star:
  two_pass_mode: Basic
  out_filter_multimap_nmax: 20

salmon:
  library_type: "A"    # Auto-detect
  bootstraps: 100      # For uncertainty
  gc_bias: true        # Bias correction
```

### Index Building Guide:
```yaml
# Building the indexes:
# STAR genome index:
#   STAR --runMode genomeGenerate \
#        --genomeDir star_index \
#        --genomeFastaFiles genome.fa \
#        --sjdbGTFfile annotation.gtf
#
# Salmon transcriptome index:
#   salmon index \
#          -t transcripts.fa \
#          -i salmon_index
```

---

## ✅ Verification Checklist

### All Parameters Validated:
- [x] salmon_index - Must exist (CRITICAL!)
- [x] STAR two_pass_mode - 'None' or 'Basic'
- [x] STAR numeric parameters - Positive integers
- [x] Salmon library_type - One of 13 valid types
- [x] bootstraps - >= 0
- [x] num_gibbs_samples - >= 0

### Two-Stage Process:
- [x] STAR outputs transcriptome BAM
- [x] Salmon reads from transcriptome BAM
- [x] Both stages have error handling
- [x] Both stages validated independently

### Error Handling:
- [x] @handle_errors on __init__
- [x] @handle_errors on run_sample
- [x] @handle_errors on combine_counts
- [x] Clear error messages for both stages

---

## 💡 Migration Steps

```bash
# Navigate to pipeline directory
cd pipelines/star_salmon/

# Backup
cp pipeline.py pipeline.py.bak
cp __init__.py __init__.py.bak
cp config.yaml config.yaml.bak

# Install v2.2.0
cp star_salmon_pipeline_v2.2.0_UPDATED.py pipeline.py
cp star_salmon___init___v2.2.0_UPDATED.py __init__.py
cp star_salmon_config_v2.2.0_UPDATED.yaml config.yaml

# Test
python -c "from raptor.pipelines.star_salmon import StarSalmonPipeline; print('✅')"
```

---

## 🔥 Use Case Examples

### Example 1: sleuth Differential Expression
```python
# STAR + Salmon is perfect for sleuth!
pipeline = StarSalmonPipeline(
    output_dir='results/',
    salmon_index='salmon_index/',  # Required!
    bootstraps=100,  # For sleuth
    library_type='A'
)
result = pipeline.run('samples.csv', 'star_index/')

# Use with sleuth:
# - Salmon output has bootstraps
# - sleuth uses bootstraps for uncertainty
# - You also get BAM files for IGV
```

### Example 2: Publication with Both Visualization and DE
```python
# Get BAM files for IGV figures
# Get bootstraps for uncertainty in DE
pipeline = StarSalmonPipeline(
    output_dir='publication/',
    salmon_index='salmon_index/',
    bootstraps=200,  # Publication quality
    gc_bias=True,
    two_pass_mode='Basic'
)
```

---

## 📊 vs Other Hybrid Approaches

### STAR + Salmon vs STAR + RSEM:
| Feature | STAR+Salmon | STAR+RSEM |
|---------|-------------|-----------|
| **Isoform-level** | No | Yes |
| **Bootstraps** | ✅ Yes | CI (optional) |
| **Speed** | Faster | Slower |
| **Use with sleuth** | ✅ Yes | No |
| **Complexity** | Medium | High |

**When to choose STAR+Salmon:** Need BAM + bootstraps, use sleuth
**When to choose STAR+RSEM:** Need isoform quantification

---

## 🎉 Achievement: 6/7 Complete!

**86% DONE - ONLY ONE MORE TO GO!**

All 6 completed pipelines have:
- ✅ Full input validation (10-17 checks each)
- ✅ Clear, actionable error messages
- ✅ @handle_errors decorator
- ✅ CLI-compatible metadata
- ✅ Comprehensive configs
- ✅ QC thresholds
- ✅ Professional quality

---

## 📞 Final Sprint!

**Only ONE pipeline left:**
- Quick pipelines (Module 1)
  - Simplified versions of existing pipelines
  - Less parameters to validate
  - Faster to complete
  - Estimated: 1-2 hours

**We're in the home stretch! 86% complete!** 🚀

---

**Version:** 2.2.0  
**Date:** January 2026  
**Author:** Ayeh Bolouki  
**Email:** ayehbolouki1988@gmail.com

**🎉 SIX PIPELINES COMPLETE! ONE MORE TO GO! 86% DONE! 🚀**
