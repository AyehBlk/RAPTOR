#!/usr/bin/env Rscript
# =============================================================================
# RAPTOR v2.2.0 - Module 6: DESeq2 Differential Expression Analysis
# =============================================================================
#
# This script performs differential expression analysis using DESeq2.
# It is designed to work with RAPTOR's recommendation output.
#
# INPUT:
#   - Count matrix (CSV/TSV): genes × samples
#   - Metadata (CSV): sample_id, condition [, batch, ...]
#   - Optional: recommendation.yaml from RAPTOR M4
#
# OUTPUT:
#   - de_results.csv: Full DE results for RAPTOR M7 import
#   - de_significant.csv: Significant genes only
#   - de_summary.json: Analysis summary
#   - de_plots.pdf: QC plots (optional)
#
# Author: Ayeh Bolouki
# Email: ayehbolouki1988@gmail.com
# Version: 2.2.0
# License: MIT
#
# Usage:
#   Rscript run_deseq2.R \
#       --counts results/gene_counts.csv \
#       --metadata data/metadata.csv \
#       --output results/de_analysis \
#       --condition condition \
#       --reference Control \
#       --fdr 0.05 \
#       --lfc 0
#
# =============================================================================

# Suppress package startup messages
suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(jsonlite)
})

# =============================================================================
# Command Line Arguments
# =============================================================================

option_list <- list(
    make_option(c("-c", "--counts"), type = "character", default = NULL,
                help = "Count matrix CSV/TSV file [required]"),
    make_option(c("-m", "--metadata"), type = "character", default = NULL,
                help = "Sample metadata CSV file [required]"),
    make_option(c("-o", "--output"), type = "character", default = "results/de_analysis",
                help = "Output directory [default: results/de_analysis]"),
    make_option(c("--config"), type = "character", default = NULL,
                help = "RAPTOR recommendation.yaml file [optional]"),
    make_option(c("--condition"), type = "character", default = "condition",
                help = "Metadata column for condition [default: condition]"),
    make_option(c("--reference"), type = "character", default = NULL,
                help = "Reference level for comparison [default: alphabetically first]"),
    make_option(c("--batch"), type = "character", default = NULL,
                help = "Metadata column for batch correction [optional]"),
    make_option(c("--fdr"), type = "double", default = 0.05,
                help = "FDR threshold [default: 0.05]"),
    make_option(c("--lfc"), type = "double", default = 0,
                help = "Log2 fold change threshold [default: 0]"),
    make_option(c("--shrinkage"), type = "character", default = "apeglm",
                help = "LFC shrinkage method: apeglm, ashr, normal, none [default: apeglm]"),
    make_option(c("--fit-type"), type = "character", default = "parametric",
                help = "Dispersion fit type: parametric, local, mean [default: parametric]"),
    make_option(c("--min-count"), type = "integer", default = 10,
                help = "Minimum count for filtering [default: 10]"),
    make_option(c("--plots"), action = "store_true", default = FALSE,
                help = "Generate QC plots [default: FALSE]"),
    make_option(c("--threads"), type = "integer", default = 1,
                help = "Number of parallel threads [default: 1]")
)

parser <- OptionParser(
    usage = "%prog [options]",
    option_list = option_list,
    description = paste(
        "\n🦖 RAPTOR v2.2.0 - DESeq2 Differential Expression Analysis",
        "\n\nModule 6 of RAPTOR workflow (Stage 3: DE Analysis)",
        "\nPerforms differential expression analysis using DESeq2.",
        "\n\nOutput will be standardized for RAPTOR Module 7 import.",
        sep = ""
    )
)

args <- parse_args(parser)

# =============================================================================
# Validation
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║       🦖 RAPTOR v2.2.0 - DESeq2 Analysis (Module 6)         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Check required arguments
if (is.null(args$counts)) {
    stop("ERROR: --counts is required. Use -h for help.")
}
if (is.null(args$metadata)) {
    stop("ERROR: --metadata is required. Use -h for help.")
}

# Check files exist
if (!file.exists(args$counts)) {
    stop(paste("ERROR: Count file not found:", args$counts))
}
if (!file.exists(args$metadata)) {
    stop(paste("ERROR: Metadata file not found:", args$metadata))
}

# =============================================================================
# Load RAPTOR Configuration (if provided)
# =============================================================================

raptor_config <- NULL
if (!is.null(args$config) && file.exists(args$config)) {
    cat("📂 Loading RAPTOR configuration:", args$config, "\n")
    
    if (requireNamespace("yaml", quietly = TRUE)) {
        raptor_config <- yaml::read_yaml(args$config)
        
        # Override defaults with RAPTOR recommendations
        if (!is.null(raptor_config$fdr_threshold)) {
            args$fdr <- raptor_config$fdr_threshold
            cat("   Using RAPTOR FDR:", args$fdr, "\n")
        }
        if (!is.null(raptor_config$lfc_threshold)) {
            args$lfc <- raptor_config$lfc_threshold
            cat("   Using RAPTOR LFC:", args$lfc, "\n")
        }
        if (!is.null(raptor_config$min_count)) {
            args$min_count <- raptor_config$min_count
            cat("   Using RAPTOR min_count:", args$min_count, "\n")
        }
    } else {
        cat("   ⚠️ yaml package not installed, skipping config\n")
    }
}

# =============================================================================
# Load Data
# =============================================================================

cat("\n📂 Loading data...\n")

# Determine delimiter
if (grepl("\\.tsv$", args$counts)) {
    counts <- read.table(args$counts, header = TRUE, sep = "\t", 
                        row.names = 1, check.names = FALSE)
} else {
    counts <- read.csv(args$counts, row.names = 1, check.names = FALSE)
}

metadata <- read.csv(args$metadata, row.names = 1, stringsAsFactors = FALSE)

cat("   Counts:", nrow(counts), "genes ×", ncol(counts), "samples\n")
cat("   Metadata:", nrow(metadata), "samples\n")

# Verify sample alignment
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0) {
    # Try matching on sample_id column
    if ("sample_id" %in% colnames(metadata)) {
        rownames(metadata) <- metadata$sample_id
        common_samples <- intersect(colnames(counts), rownames(metadata))
    }
}

if (length(common_samples) < ncol(counts)) {
    cat("   ⚠️ Subsetting to", length(common_samples), "matching samples\n")
}

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Ensure condition column exists
if (!args$condition %in% colnames(metadata)) {
    stop(paste("ERROR: Condition column not found:", args$condition))
}

# Convert to factor
metadata[[args$condition]] <- factor(metadata[[args$condition]])

# Set reference level if specified
if (!is.null(args$reference)) {
    if (args$reference %in% levels(metadata[[args$condition]])) {
        metadata[[args$condition]] <- relevel(metadata[[args$condition]], 
                                              ref = args$reference)
        cat("   Reference level:", args$reference, "\n")
    } else {
        cat("   ⚠️ Reference level not found, using:", 
            levels(metadata[[args$condition]])[1], "\n")
    }
}

cat("   Conditions:", paste(levels(metadata[[args$condition]]), collapse = " vs "), "\n")

# =============================================================================
# Pre-filtering
# =============================================================================

cat("\n🔍 Filtering low-count genes...\n")

# Keep genes with sufficient counts
keep <- rowSums(counts >= args$min_count) >= 2  # At least 2 samples
counts_filtered <- counts[keep, ]

cat("   Before:", nrow(counts), "genes\n")
cat("   After:", nrow(counts_filtered), "genes\n")
cat("   Removed:", nrow(counts) - nrow(counts_filtered), "low-count genes\n")

# =============================================================================
# DESeq2 Analysis
# =============================================================================

cat("\n⚙️ Running DESeq2 analysis...\n")

# Create design formula
if (!is.null(args$batch) && args$batch %in% colnames(metadata)) {
    design_formula <- as.formula(paste0("~ ", args$batch, " + ", args$condition))
    cat("   Design: ~ ", args$batch, " + ", args$condition, " (batch-corrected)\n")
} else {
    design_formula <- as.formula(paste0("~ ", args$condition))
    cat("   Design: ~ ", args$condition, "\n")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts_filtered),
    colData = metadata,
    design = design_formula
)

# Set parallel processing
if (args$threads > 1) {
    cat("   Using", args$threads, "threads\n")
    BPPARAM <- BiocParallel::MulticoreParam(args$threads)
} else {
    BPPARAM <- BiocParallel::SerialParam()
}

# Run DESeq2
cat("   Running DESeq2 pipeline...\n")
dds <- DESeq(dds, fitType = args$`fit-type`, parallel = TRUE, BPPARAM = BPPARAM)

# =============================================================================
# Extract Results
# =============================================================================

cat("\n📊 Extracting results...\n")

# Get results with LFC threshold if specified
if (args$lfc > 0) {
    res <- results(dds, alpha = args$fdr, lfcThreshold = args$lfc)
    cat("   Using lfcThreshold =", args$lfc, "\n")
} else {
    res <- results(dds, alpha = args$fdr)
}

# Apply shrinkage if requested
if (args$shrinkage != "none") {
    cat("   Applying LFC shrinkage:", args$shrinkage, "\n")
    
    # Get coefficient name
    coef_name <- resultsNames(dds)[2]  # First contrast after intercept
    
    if (args$shrinkage == "apeglm") {
        if (!requireNamespace("apeglm", quietly = TRUE)) {
            cat("   ⚠️ apeglm not installed, using normal shrinkage\n")
            res_shrunk <- lfcShrink(dds, coef = coef_name, type = "normal")
        } else {
            res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
        }
    } else if (args$shrinkage == "ashr") {
        if (!requireNamespace("ashr", quietly = TRUE)) {
            cat("   ⚠️ ashr not installed, using normal shrinkage\n")
            res_shrunk <- lfcShrink(dds, coef = coef_name, type = "normal")
        } else {
            res_shrunk <- lfcShrink(dds, coef = coef_name, type = "ashr")
        }
    } else {
        res_shrunk <- lfcShrink(dds, coef = coef_name, type = "normal")
    }
    
    # Keep original p-values but use shrunk LFCs
    res$log2FoldChange <- res_shrunk$log2FoldChange
    if (!is.null(res_shrunk$lfcSE)) {
        res$lfcSE <- res_shrunk$lfcSE
    }
}

# =============================================================================
# Prepare Output
# =============================================================================

cat("\n📁 Preparing output...\n")

# Create output directory
output_dir <- args$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Reorder columns for RAPTOR compatibility
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", 
                     "stat", "pvalue", "padj")]

# Add direction column
res_df$direction <- ifelse(res_df$log2FoldChange > 0, "up",
                          ifelse(res_df$log2FoldChange < 0, "down", "unchanged"))

# Mark significant genes
res_df$is_significant <- !is.na(res_df$padj) & 
                        res_df$padj < args$fdr &
                        abs(res_df$log2FoldChange) > args$lfc

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# Reset row names
rownames(res_df) <- res_df$gene_id

# =============================================================================
# Save Results
# =============================================================================

# Full results
output_file <- file.path(output_dir, "de_results.csv")
write.csv(res_df, output_file, row.names = FALSE)
cat("   ✓ de_results.csv\n")

# Significant genes only
sig_df <- res_df[res_df$is_significant, ]
sig_file <- file.path(output_dir, "de_significant.csv")
write.csv(sig_df, sig_file, row.names = FALSE)
cat("   ✓ de_significant.csv (", nrow(sig_df), " genes)\n", sep = "")

# Summary JSON
n_up <- sum(sig_df$direction == "up", na.rm = TRUE)
n_down <- sum(sig_df$direction == "down", na.rm = TRUE)

summary_list <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    raptor_version = "2.2.0",
    module = "M6",
    pipeline = "DESeq2",
    deseq2_version = as.character(packageVersion("DESeq2")),
    parameters = list(
        fdr_threshold = args$fdr,
        lfc_threshold = args$lfc,
        shrinkage = args$shrinkage,
        fit_type = args$`fit-type`,
        min_count = args$min_count,
        batch_correction = !is.null(args$batch)
    ),
    data = list(
        n_samples = ncol(counts_filtered),
        n_genes_input = nrow(counts),
        n_genes_tested = nrow(counts_filtered),
        condition_column = args$condition,
        reference_level = levels(metadata[[args$condition]])[1],
        comparison = paste(rev(levels(metadata[[args$condition]])), collapse = "_vs_")
    ),
    results = list(
        n_significant = nrow(sig_df),
        n_upregulated = n_up,
        n_downregulated = n_down,
        pct_significant = round(100 * nrow(sig_df) / nrow(res_df), 2)
    )
)

summary_file <- file.path(output_dir, "de_summary.json")
write_json(summary_list, summary_file, pretty = TRUE, auto_unbox = TRUE)
cat("   ✓ de_summary.json\n")

# =============================================================================
# Generate Plots (Optional)
# =============================================================================

if (args$plots) {
    cat("\n📊 Generating QC plots...\n")
    
    plot_file <- file.path(output_dir, "de_plots.pdf")
    pdf(plot_file, width = 10, height = 8)
    
    # MA plot
    DESeq2::plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))
    abline(h = c(-args$lfc, args$lfc), col = "blue", lty = 2)
    
    # Dispersion plot
    DESeq2::plotDispEsts(dds, main = "Dispersion Estimates")
    
    # PCA plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        vsd <- vst(dds, blind = FALSE)
        pca_data <- DESeq2::plotPCA(vsd, intgroup = args$condition, returnData = TRUE)
        percentVar <- round(100 * attr(pca_data, "percentVar"))
        
        p <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, color = .data[[args$condition]])) +
            ggplot2::geom_point(size = 3) +
            ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
            ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
            ggplot2::ggtitle("PCA of VST-transformed counts") +
            ggplot2::theme_minimal()
        print(p)
    }
    
    # Volcano plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        volcano_df <- res_df
        volcano_df$neg_log10_padj <- -log10(volcano_df$padj)
        volcano_df$significance <- ifelse(volcano_df$is_significant,
                                         ifelse(volcano_df$direction == "up", "Up", "Down"),
                                         "NS")
        
        p <- ggplot2::ggplot(volcano_df, 
                            ggplot2::aes(log2FoldChange, neg_log10_padj, color = significance)) +
            ggplot2::geom_point(alpha = 0.5, size = 1) +
            ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
            ggplot2::geom_vline(xintercept = c(-args$lfc, args$lfc), linetype = "dashed") +
            ggplot2::geom_hline(yintercept = -log10(args$fdr), linetype = "dashed") +
            ggplot2::xlab("log2 Fold Change") +
            ggplot2::ylab("-log10(adjusted p-value)") +
            ggplot2::ggtitle("Volcano Plot") +
            ggplot2::theme_minimal()
        print(p)
    }
    
    dev.off()
    cat("   ✓ de_plots.pdf\n")
}

# =============================================================================
# Final Summary
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  ✅ DESeq2 ANALYSIS COMPLETE!\n")
cat("═══════════════════════════════════════════════════════════════\n")

cat("\n  📊 Summary:\n")
cat("     • Total genes tested:", nrow(res_df), "\n")
cat("     • Significant genes:", nrow(sig_df), 
    "(", round(100 * nrow(sig_df) / nrow(res_df), 1), "%)\n")
cat("     • Upregulated:", n_up, "\n")
cat("     • Downregulated:", n_down, "\n")

cat("\n  📁 Output Directory:", output_dir, "\n")
cat("     • de_results.csv (full results for RAPTOR M7)\n")
cat("     • de_significant.csv (significant genes only)\n")
cat("     • de_summary.json (analysis metadata)\n")
if (args$plots) cat("     • de_plots.pdf (QC plots)\n")

cat("\n  🔜 Next Steps:\n")
cat("     Import results into RAPTOR:\n")
cat("     raptor import-de --de-file", output_file, "\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Making free science for everybody around the world 🌍\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
