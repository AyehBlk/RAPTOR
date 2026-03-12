# 🔬 RAPTOR Module 6 - R Scripts Analysis & Recommendations

**Date:** January 19, 2026  
**Module:** M6 - Differential Expression Analysis (External R Module)  
**Status:** ✅ Scripts are comprehensive but need additional utilities

---

## 📋 TABLE OF CONTENTS

1. [Current Scripts Analysis](#current-scripts-analysis)
2. [What's Good](#whats-good)
3. [What's Missing](#whats-missing)
4. [Additional Scripts Needed](#additional-scripts-needed)
5. [Recommendations for Improvements](#recommendations-for-improvements)
6. [Complete File Structure](#complete-file-structure)

---

## 📊 CURRENT SCRIPTS ANALYSIS

### ✅ Scripts You Have (3 files)

| Script | Lines | Status | Quality |
|--------|-------|--------|---------|
| run_deseq2.R | 472 | ✅ Comprehensive | Excellent |
| run_edger.R | 391 | ✅ Comprehensive | Excellent |
| run_limma.R | 402 | ✅ Comprehensive | Excellent |

---

## 🎯 WHAT'S GOOD

### 1. **Excellent Structure** ✅
All three scripts follow consistent patterns:
- Clear header with documentation
- Command-line argument parsing with optparse
- RAPTOR configuration integration (recommendation.yaml)
- Comprehensive error checking
- Standardized output format (de_results.csv, de_significant.csv, de_summary.json)
- Optional plot generation

### 2. **RAPTOR Integration** ✅
- Load RAPTOR recommendation.yaml
- Override defaults with RAPTOR parameters
- Standardized output for Module 7 import
- Consistent JSON metadata format

### 3. **Comprehensive Options** ✅

**DESeq2:**
- ✅ Batch correction support
- ✅ LFC shrinkage (apeglm, ashr, normal)
- ✅ Fit type selection (parametric, local, mean)
- ✅ Min count filtering
- ✅ Parallel threads support
- ✅ QC plots (MA, dispersion, PCA, volcano)

**edgeR:**
- ✅ Multiple normalization methods (TMM, RLE, upperquartile)
- ✅ Multiple test methods (QLF, LRT, exact)
- ✅ Robust dispersion estimation
- ✅ Batch correction support
- ✅ QC plots (BCV, MD, MDS)

**limma-voom:**
- ✅ voomWithQualityWeights option
- ✅ Robust empirical Bayes
- ✅ Mean-variance trend
- ✅ Batch correction support
- ✅ QC plots (voom, MD, SA, volcano)

### 4. **Output Standardization** ✅
All scripts produce identical output structure:
```
results/de_analysis/
├── de_results.csv        # Full results with standardized columns
├── de_significant.csv    # Significant genes only
├── de_summary.json       # Metadata for RAPTOR M7
└── de_plots.pdf          # Optional QC plots
```

**Standardized columns:**
- gene_id
- log2FoldChange
- pvalue
- padj
- baseMean
- stat
- direction (up/down/unchanged)
- is_significant

### 5. **Robust Error Handling** ✅
- File existence checks
- Package availability checks (apeglm, ashr, yaml)
- Sample alignment validation
- Condition column validation
- Reference level handling

---

## ⚠️ WHAT'S MISSING

### 1. **Helper/Utility Scripts** ❌

Your main analysis scripts are excellent, but you need supporting utilities:

#### Missing: **install_packages.R**
Purpose: Install all required Bioconductor packages
- DESeq2, edgeR, limma, apeglm, ashr
- Check versions, dependencies
- Test installations

#### Missing: **validate_inputs.R**
Purpose: Pre-validate inputs before running analysis
- Check count matrix format
- Check metadata format
- Verify sample alignment
- Check for outliers/low counts
- Suggest appropriate method

#### Missing: **compare_methods.R**
Purpose: Run all 3 methods and compare results
- Venn diagrams of significant genes
- Correlation of log fold changes
- Method concordance analysis
- Help users choose best method

#### Missing: **generate_report.Rmd**
Purpose: Generate comprehensive HTML report
- QC plots from all methods
- DE results summary
- Top genes tables
- Gene set enrichment preview
- Export-ready figures

### 2. **Python Integration Scripts** ❌

Since RAPTOR is primarily Python, you need Python wrappers:

#### Missing: **run_de_analysis.py**
Purpose: Python wrapper to run R scripts
- Call R scripts from Python
- Pass RAPTOR parameters
- Handle errors
- Parse JSON output
- Import results to RAPTOR

#### Missing: **de_pipeline.py**
Purpose: High-level DE pipeline manager
- Select appropriate method based on data
- Run analysis
- Generate reports
- Import to RAPTOR M7

### 3. **Configuration Files** ❌

#### Missing: **de_config.yaml**
Purpose: Central configuration for DE analysis
```yaml
deseq2:
  default_shrinkage: apeglm
  default_fit_type: parametric
  default_min_count: 10

edger:
  default_normalization: TMM
  default_method: QLF
  robust: true

limma:
  default_voom_method: voom
  robust_ebayes: false
  trend: false

thresholds:
  fdr: 0.05
  lfc: 0.0
  min_count: 10
```

### 4. **Testing Scripts** ❌

#### Missing: **test_de_scripts.R**
Purpose: Unit tests for R scripts
- Test with synthetic data
- Test error handling
- Test output formats
- Verify RAPTOR compatibility

### 5. **Documentation** ❌

#### Missing: **README.md**
Purpose: Module 6 documentation
- Overview of each method
- When to use which method
- Installation instructions
- Usage examples
- Troubleshooting

#### Missing: **METHODS_COMPARISON.md**
Purpose: Help users choose methods
- DESeq2 vs edgeR vs limma-voom comparison
- Pros/cons of each
- Performance benchmarks
- Recommendation flowchart

---

## 📦 ADDITIONAL SCRIPTS NEEDED

### Priority 1 (Essential): 🔴

#### 1. **install_packages.R** (Essential)
```r
#!/usr/bin/env Rscript
# Install all required packages for RAPTOR Module 6

# Check for BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Required packages
packages <- c(
    "DESeq2",
    "edgeR", 
    "limma",
    "apeglm",     # For DESeq2 shrinkage
    "ashr",       # For DESeq2 shrinkage
    "optparse",   # CLI parsing
    "jsonlite",   # JSON output
    "yaml",       # Config files
    "ggplot2"     # Plotting
)

# Install missing packages
for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing:", pkg, "\n")
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
        cat("✓", pkg, "already installed\n")
    }
}

# Print versions
cat("\n📊 Package Versions:\n")
for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  ", pkg, ":", 
            as.character(packageVersion(pkg)), "\n")
    }
}
```

#### 2. **run_de_analysis.py** (Essential for RAPTOR integration)
```python
#!/usr/bin/env python3
"""
RAPTOR Module 6 - R Script Wrapper

Calls R scripts for differential expression analysis.
"""

import subprocess
import json
from pathlib import Path

class DEAnalysis:
    """Wrapper for R differential expression scripts."""
    
    def __init__(self, method='deseq2'):
        self.method = method.lower()
        self.script_dir = Path(__file__).parent
        
    def run(self, counts, metadata, output_dir, **kwargs):
        """Run DE analysis with specified method."""
        
        # Map to R script
        scripts = {
            'deseq2': 'run_deseq2.R',
            'edger': 'run_edger.R',
            'limma': 'run_limma.R'
        }
        
        script = self.script_dir / scripts[self.method]
        
        # Build command
        cmd = [
            'Rscript', str(script),
            '--counts', counts,
            '--metadata', metadata,
            '--output', output_dir
        ]
        
        # Add optional arguments
        for key, value in kwargs.items():
            cmd.extend([f'--{key.replace("_", "-")}', str(value)])
        
        # Run R script
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"DE analysis failed: {result.stderr}")
        
        # Parse results
        summary_file = Path(output_dir) / 'de_summary.json'
        with open(summary_file) as f:
            summary = json.load(f)
        
        return summary

# CLI interface
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='RAPTOR M6 - Differential Expression Analysis'
    )
    parser.add_argument('--method', choices=['deseq2', 'edger', 'limma'],
                       default='deseq2')
    parser.add_argument('--counts', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--output', default='results/de_analysis')
    parser.add_argument('--condition', default='condition')
    parser.add_argument('--fdr', type=float, default=0.05)
    
    args = parser.parse_args()
    
    analysis = DEAnalysis(method=args.method)
    summary = analysis.run(
        counts=args.counts,
        metadata=args.metadata,
        output_dir=args.output,
        condition=args.condition,
        fdr=args.fdr
    )
    
    print(f"✅ Analysis complete!")
    print(f"   Significant genes: {summary['results']['n_significant']}")
```

#### 3. **compare_methods.R** (Very useful)
```r
#!/usr/bin/env Rscript
# Compare results from DESeq2, edgeR, and limma-voom

library(ggplot2)

# Load all three result files
deseq2 <- read.csv("results/deseq2/de_results.csv")
edger <- read.csv("results/edger/de_results.csv")
limma <- read.csv("results/limma/de_results.csv")

# Find significant genes in each
sig_deseq2 <- deseq2$gene_id[deseq2$is_significant]
sig_edger <- edger$gene_id[edger$is_significant]
sig_limma <- limma$gene_id[limma$is_significant]

# Venn diagram
library(VennDiagram)
venn.diagram(
    x = list(
        DESeq2 = sig_deseq2,
        edgeR = sig_edger,
        limma = sig_limma
    ),
    filename = "method_comparison_venn.png",
    main = "Significant Genes: Method Comparison"
)

# Correlation of log fold changes
common_genes <- Reduce(intersect, list(
    deseq2$gene_id, edger$gene_id, limma$gene_id
))

lfc_comparison <- data.frame(
    gene_id = common_genes,
    deseq2_lfc = deseq2$log2FoldChange[match(common_genes, deseq2$gene_id)],
    edger_lfc = edger$log2FoldChange[match(common_genes, edger$gene_id)],
    limma_lfc = limma$log2FoldChange[match(common_genes, limma$gene_id)]
)

# Correlation matrix
cor_matrix <- cor(lfc_comparison[, -1], use = "complete.obs")
print("LFC Correlation Matrix:")
print(round(cor_matrix, 3))

# Generate comparison report
cat("\n═══════════════════════════════════════════════\n")
cat("  DE Methods Comparison Summary\n")
cat("═══════════════════════════════════════════════\n\n")
cat("Significant genes:\n")
cat("  DESeq2:", length(sig_deseq2), "\n")
cat("  edgeR:", length(sig_edger), "\n")
cat("  limma-voom:", length(sig_limma), "\n")
cat("\nAgreement:\n")
cat("  All 3 methods:", length(Reduce(intersect, list(sig_deseq2, sig_edger, sig_limma))), "\n")
cat("  DESeq2 & edgeR:", length(intersect(sig_deseq2, sig_edger)), "\n")
cat("  DESeq2 & limma:", length(intersect(sig_deseq2, sig_limma)), "\n")
cat("  edgeR & limma:", length(intersect(sig_edger, sig_limma)), "\n")
```

### Priority 2 (Important): 🟡

#### 4. **validate_inputs.R**
```r
#!/usr/bin/env Rscript
# Pre-validate inputs before DE analysis

validate_count_matrix <- function(counts_file) {
    # Check format, dimensions, missing values
    # Check for low counts
    # Recommend filtering
}

validate_metadata <- function(metadata_file, counts_file) {
    # Check sample alignment
    # Check condition column
    # Check for batch effects
    # Suggest covariates
}

recommend_method <- function(counts, metadata) {
    # Based on sample size, counts, design
    # Return recommended method with rationale
}
```

#### 5. **generate_report.Rmd**
```r
---
title: "RAPTOR Module 6 - Differential Expression Report"
output: 
  html_document:
    toc: true
    theme: united
---

## Overview
Analysis: `r params$method`
Date: `r Sys.Date()`

## QC Plots
[Include all QC plots]

## Results Summary
[Summary tables]

## Top DE Genes
[Top 50 up/down regulated]
```

### Priority 3 (Nice to have): 🟢

#### 6. **batch_correction.R**
For explicit batch correction with ComBat-seq or similar

#### 7. **gsea_enrichment.R**
Quick gene set enrichment analysis on DE results

#### 8. **export_for_tools.R**
Export DE results for external tools (IPA, GSEA, etc.)

---

## 💡 RECOMMENDATIONS FOR IMPROVEMENTS

### Minor Issues in Current Scripts:

#### 1. **DESeq2 Script - Line 207** ⚠️
```r
# Missing code - truncated view
# Need to see filtering implementation
```
**Action:** Verify pre-filtering code is complete

#### 2. **All Scripts - Package Dependencies** ⚠️
```r
# Currently checks at runtime
# Should check all packages at start
```
**Add to all scripts:**
```r
# Check required packages
required_pkgs <- c("optparse", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
    stop(paste("Missing required packages:", paste(missing_pkgs, collapse = ", "),
               "\nRun: Rscript install_packages.R"))
}
```

#### 3. **Multi-Factor Designs** 🆕
Current scripts support basic batch correction, but consider adding:
- Multiple factor designs (~treatment + time + batch)
- Interaction terms
- Paired/blocked designs
- Time series support

**Add to all scripts:**
```r
make_option(c("--formula"), type = "character", default = NULL,
            help = "Custom model formula (advanced) [optional]")

# Then in design:
if (!is.null(args$formula)) {
    design <- model.matrix(as.formula(args$formula), data = metadata)
} else {
    # Current simple design
}
```

#### 4. **Sample Outlier Detection** 🆕
Add automatic outlier detection:
```r
# Add to all scripts after normalization
detect_outliers <- function(dds) {
    # Use PCA distance or Cook's distance
    # Flag potential outliers
    # Allow --remove-outliers flag
}
```

#### 5. **Gene Filtering Reporting** 🆕
Add more detail about filtering:
```r
# After filtering, report:
cat("   Genes filtered:\n")
cat("     • Zero counts:", n_zero, "\n")
cat("     • Low counts:", n_low, "\n")
cat("     • Passed:", n_kept, "\n")
cat("     • Filtering stringency:", stringency, "\n")
```

---

## 📂 COMPLETE FILE STRUCTURE

### Recommended Module 6 Structure:

```
RAPTOR/raptor/external_modules/
└── module6_de_analysis/
    ├── README.md                          # Module documentation
    ├── METHODS_COMPARISON.md              # Help users choose methods
    ├── de_config.yaml                     # Central configuration
    │
    ├── r_scripts/                         # R analysis scripts
    │   ├── run_deseq2.R                   # ✅ You have
    │   ├── run_edger.R                    # ✅ You have
    │   ├── run_limma.R                    # ✅ You have
    │   ├── install_packages.R             # 🔴 NEED
    │   ├── validate_inputs.R              # 🟡 NEED
    │   ├── compare_methods.R              # 🔴 NEED
    │   ├── batch_correction.R             # 🟢 Optional
    │   └── gsea_enrichment.R              # 🟢 Optional
    │
    ├── python_wrappers/                   # Python integration
    │   ├── __init__.py
    │   ├── run_de_analysis.py             # 🔴 NEED
    │   ├── de_pipeline.py                 # 🟡 NEED
    │   └── de_importer.py                 # 🟡 NEED (for M7)
    │
    ├── reports/                           # Report generation
    │   ├── generate_report.Rmd            # 🟡 NEED
    │   └── report_template.html
    │
    ├── tests/                             # Testing
    │   ├── test_de_scripts.R              # 🟡 NEED
    │   ├── test_python_wrappers.py        # 🟡 NEED
    │   └── synthetic_data/
    │       ├── test_counts.csv
    │       └── test_metadata.csv
    │
    └── examples/                          # Usage examples
        ├── example_deseq2.sh
        ├── example_edger.sh
        ├── example_limma.sh
        └── example_compare.sh
```

---

## ✅ SUMMARY & NEXT STEPS

### What You Have: ✅
- ✅ Excellent main analysis scripts (DESeq2, edgeR, limma)
- ✅ Proper RAPTOR integration
- ✅ Standardized outputs
- ✅ Comprehensive options
- ✅ Good error handling

### What You Need: 🎯

**Priority 1 (Essential):**
1. 🔴 **install_packages.R** - Install dependencies
2. 🔴 **run_de_analysis.py** - Python wrapper for RAPTOR
3. 🔴 **compare_methods.R** - Method comparison
4. 🔴 **README.md** - Documentation

**Priority 2 (Important):**
5. 🟡 **validate_inputs.R** - Pre-validation
6. 🟡 **de_pipeline.py** - High-level pipeline
7. 🟡 **generate_report.Rmd** - HTML reports
8. 🟡 **de_config.yaml** - Central config

**Priority 3 (Nice to have):**
9. 🟢 **batch_correction.R** - Advanced batch correction
10. 🟢 **gsea_enrichment.R** - Quick enrichment
11. 🟢 Test scripts and examples

### Minor Improvements to Current Scripts:
- Add package checking at start
- Add multi-factor design support
- Add outlier detection
- Enhanced filtering reports
- More detailed logging

---

## 🎯 RECOMMENDATIONS

### Immediate Actions:

1. **Create install_packages.R** - Users need this first
2. **Create run_de_analysis.py** - Essential for RAPTOR integration
3. **Create README.md** - Document the module
4. **Add package checks** to existing R scripts
5. **Test scripts** with sample data

### Future Enhancements:

1. **Multi-factor designs** - More complex experiments
2. **Time series support** - For temporal data
3. **Single-cell support** - For scRNA-seq
4. **Web interface** - Interactive DE analysis
5. **Pre-built reports** - Publication-ready

---

## 📊 OVERALL ASSESSMENT

**Rating: 8/10** ⭐⭐⭐⭐⭐⭐⭐⭐☆☆

**Strengths:**
- ✅ Core analysis scripts are excellent
- ✅ Proper RAPTOR integration
- ✅ Comprehensive options
- ✅ Standardized outputs

**Areas for Improvement:**
- ⚠️ Missing supporting utilities
- ⚠️ No Python integration yet
- ⚠️ No documentation
- ⚠️ No testing framework

**Recommendation:**
Your R scripts are **production-ready**, but you need to build the **ecosystem** around them:
- Installation scripts
- Python wrappers
- Documentation
- Testing

**Timeline:**
- Priority 1 items: 1-2 days
- Priority 2 items: 2-3 days  
- Priority 3 items: 3-5 days
- **Total: ~1 week for complete M6**

---

**Status:** ✅ Core scripts are excellent - now build the supporting infrastructure!

**Next:** I can help you create any of the missing scripts. Which would you like first?

