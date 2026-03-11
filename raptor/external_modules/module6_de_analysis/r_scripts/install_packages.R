#!/usr/bin/env Rscript
# =============================================================================
# RAPTOR v2.2.0 - Module 6: Package Installation
# =============================================================================
#
# This script installs all required R/Bioconductor packages for Module 6
# differential expression analysis.
#
# Author: Ayeh Bolouki
# Email: ayehbolouki1988@gmail.com
# Version: 2.2.0
# License: MIT
#
# Usage:
#   Rscript install_packages.R [--test]
#
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║    🦖 RAPTOR v2.2.0 - Installing Module 6 Packages          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# =============================================================================
# Check and Install BiocManager
# =============================================================================

cat("📦 Checking BiocManager...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("   Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
} else {
    cat("   ✓ BiocManager already installed\n")
}

# Set Bioconductor version
BiocManager::install(version = "3.18", update = FALSE, ask = FALSE)

# =============================================================================
# Required Packages
# =============================================================================

# Core DE analysis packages
core_packages <- c(
    "DESeq2",      # DESeq2 analysis
    "edgeR",       # edgeR analysis
    "limma"        # limma-voom analysis
)

# Shrinkage methods for DESeq2
shrinkage_packages <- c(
    "apeglm",      # Adaptive t prior shrinkage
    "ashr"         # Adaptive shrinkage
)

# Utility packages
utility_packages <- c(
    "optparse",    # Command-line parsing
    "jsonlite",    # JSON output
    "yaml"         # YAML config files
)

# Plotting packages
plotting_packages <- c(
    "ggplot2",     # Advanced plotting
    "pheatmap",    # Heatmaps
    "RColorBrewer" # Color palettes
)

# Optional packages for enhanced functionality
optional_packages <- c(
    "VennDiagram",       # Venn diagrams for comparison
    "ggrepel",           # Better plot labels
    "ComplexHeatmap",    # Advanced heatmaps
    "clusterProfiler",   # Gene set enrichment
    "org.Hs.eg.db",      # Human gene annotations
    "org.Mm.eg.db"       # Mouse gene annotations
)

# All required packages
all_required <- c(
    core_packages,
    shrinkage_packages,
    utility_packages,
    plotting_packages
)

# =============================================================================
# Installation Function
# =============================================================================

install_if_missing <- function(packages, optional = FALSE) {
    for (pkg in packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            version <- as.character(packageVersion(pkg))
            cat(sprintf("   ✓ %-20s %s\n", pkg, version))
        } else {
            if (optional) {
                cat(sprintf("   ○ %-20s (optional, skipping)\n", pkg))
            } else {
                cat(sprintf("   ⏳ Installing %s...\n", pkg))
                
                tryCatch({
                    BiocManager::install(pkg, 
                                       update = FALSE, 
                                       ask = FALSE,
                                       quiet = TRUE)
                    version <- as.character(packageVersion(pkg))
                    cat(sprintf("   ✓ %-20s %s (installed)\n", pkg, version))
                }, error = function(e) {
                    cat(sprintf("   ✗ %-20s FAILED: %s\n", pkg, e$message))
                })
            }
        }
    }
}

# =============================================================================
# Install Packages
# =============================================================================

cat("\n📦 Installing required packages...\n")

cat("\n  Core Analysis Packages:\n")
install_if_missing(core_packages)

cat("\n  Shrinkage Methods:\n")
install_if_missing(shrinkage_packages)

cat("\n  Utility Packages:\n")
install_if_missing(utility_packages)

cat("\n  Plotting Packages:\n")
install_if_missing(plotting_packages)

if (!test_mode) {
    cat("\n  Optional Packages (for enhanced functionality):\n")
    install_if_missing(optional_packages, optional = TRUE)
}

# =============================================================================
# Test Installations
# =============================================================================

if (test_mode || "--test" %in% args) {
    cat("\n🧪 Testing installations...\n")
    
    test_results <- list()
    
    for (pkg in all_required) {
        test_results[[pkg]] <- requireNamespace(pkg, quietly = TRUE)
    }
    
    # Summary
    n_success <- sum(unlist(test_results))
    n_total <- length(test_results)
    
    if (n_success == n_total) {
        cat(sprintf("   ✅ All %d required packages loaded successfully!\n", n_total))
    } else {
        cat(sprintf("   ⚠️ %d/%d packages loaded successfully\n", n_success, n_total))
        
        failed <- names(test_results)[!unlist(test_results)]
        cat("   Failed packages:\n")
        for (pkg in failed) {
            cat(sprintf("     ✗ %s\n", pkg))
        }
    }
    
    # Test core functionality
    cat("\n🔬 Testing core functionality...\n")
    
    # Test DESeq2
    tryCatch({
        library(DESeq2, quietly = TRUE)
        cat("   ✓ DESeq2 loads correctly\n")
    }, error = function(e) {
        cat("   ✗ DESeq2 failed:", e$message, "\n")
    })
    
    # Test edgeR
    tryCatch({
        library(edgeR, quietly = TRUE)
        cat("   ✓ edgeR loads correctly\n")
    }, error = function(e) {
        cat("   ✗ edgeR failed:", e$message, "\n")
    })
    
    # Test limma
    tryCatch({
        library(limma, quietly = TRUE)
        cat("   ✓ limma loads correctly\n")
    }, error = function(e) {
        cat("   ✗ limma failed:", e$message, "\n")
    })
}

# =============================================================================
# Version Report
# =============================================================================

cat("\n📊 Installed Package Versions:\n")
cat("════════════════════════════════════════════════════════════\n")

for (pkg in all_required) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        version <- as.character(packageVersion(pkg))
        cat(sprintf("  %-25s %s\n", pkg, version))
    }
}

# =============================================================================
# System Information
# =============================================================================

cat("\n💻 System Information:\n")
cat("════════════════════════════════════════════════════════════\n")
cat("  R version:", R.version.string, "\n")
cat("  Bioconductor version:", as.character(BiocManager::version()), "\n")
cat("  Platform:", R.version$platform, "\n")

# =============================================================================
# Final Summary
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  ✅ PACKAGE INSTALLATION COMPLETE!\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")
cat("  🔬 You can now run RAPTOR Module 6 DE analysis:\n")
cat("     • Rscript run_deseq2.R --help\n")
cat("     • Rscript run_edger.R --help\n")
cat("     • Rscript run_limma.R --help\n")
cat("\n")
cat("  📚 Optional packages not installed?\n")
cat("     Rerun without --test flag:\n")
cat("     Rscript install_packages.R\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Making free science for everybody around the world 🌍\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n")
