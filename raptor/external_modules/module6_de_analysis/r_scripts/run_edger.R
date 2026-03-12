#!/usr/bin/env Rscript
# =============================================================================
# RAPTOR v2.2.0 - Module 6: edgeR Differential Expression Analysis
# =============================================================================
#
# This script performs differential expression analysis using edgeR.
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
#
# Author: Ayeh Bolouki
# Email: ayehbolouki1988@gmail.com
# Version: 2.2.0
# License: MIT
#
# Usage:
#   Rscript run_edger.R \
#       --counts results/gene_counts.csv \
#       --metadata data/metadata.csv \
#       --output results/de_analysis \
#       --condition condition \
#       --reference Control \
#       --fdr 0.05
#
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(edgeR)
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
    make_option(c("--normalization"), type = "character", default = "TMM",
                help = "Normalization method: TMM, RLE, upperquartile, none [default: TMM]"),
    make_option(c("--min-count"), type = "integer", default = 10,
                help = "Minimum count for filterByExpr [default: 10]"),
    make_option(c("--min-total-count"), type = "integer", default = 15,
                help = "Minimum total count [default: 15]"),
    make_option(c("--robust"), action = "store_true", default = FALSE,
                help = "Use robust dispersion estimation [default: FALSE]"),
    make_option(c("--method"), type = "character", default = "QLF",
                help = "Test method: QLF, LRT, exact [default: QLF]"),
    make_option(c("--plots"), action = "store_true", default = FALSE,
                help = "Generate QC plots [default: FALSE]")
)

parser <- OptionParser(
    usage = "%prog [options]",
    option_list = option_list,
    description = paste(
        "\n🦖 RAPTOR v2.2.0 - edgeR Differential Expression Analysis",
        "\n\nModule 6 of RAPTOR workflow (Stage 3: DE Analysis)",
        "\nPerforms differential expression analysis using edgeR.",
        "\n\nOutput will be standardized for RAPTOR Module 7 import.",
        sep = ""
    )
)

args <- parse_args(parser)

# =============================================================================
# Banner
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║       🦖 RAPTOR v2.2.0 - edgeR Analysis (Module 6)          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Check required arguments
if (is.null(args$counts)) stop("ERROR: --counts is required")
if (is.null(args$metadata)) stop("ERROR: --metadata is required")

if (!file.exists(args$counts)) stop(paste("ERROR: Count file not found:", args$counts))
if (!file.exists(args$metadata)) stop(paste("ERROR: Metadata file not found:", args$metadata))

# =============================================================================
# Load RAPTOR Configuration
# =============================================================================

if (!is.null(args$config) && file.exists(args$config)) {
    cat("📂 Loading RAPTOR configuration:", args$config, "\n")
    
    if (requireNamespace("yaml", quietly = TRUE)) {
        raptor_config <- yaml::read_yaml(args$config)
        
        if (!is.null(raptor_config$fdr_threshold)) args$fdr <- raptor_config$fdr_threshold
        if (!is.null(raptor_config$lfc_threshold)) args$lfc <- raptor_config$lfc_threshold
        if (!is.null(raptor_config$min_count)) args$min_count <- raptor_config$min_count
        if (!is.null(raptor_config$normalization)) args$normalization <- raptor_config$normalization
    }
}

# =============================================================================
# Load Data
# =============================================================================

cat("\n📂 Loading data...\n")

if (grepl("\\.tsv$", args$counts)) {
    counts <- read.table(args$counts, header = TRUE, sep = "\t", 
                        row.names = 1, check.names = FALSE)
} else {
    counts <- read.csv(args$counts, row.names = 1, check.names = FALSE)
}

metadata <- read.csv(args$metadata, row.names = 1, stringsAsFactors = FALSE)

cat("   Counts:", nrow(counts), "genes ×", ncol(counts), "samples\n")
cat("   Metadata:", nrow(metadata), "samples\n")

# Align samples
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0) {
    if ("sample_id" %in% colnames(metadata)) {
        rownames(metadata) <- metadata$sample_id
        common_samples <- intersect(colnames(counts), rownames(metadata))
    }
}

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Ensure condition column exists
if (!args$condition %in% colnames(metadata)) {
    stop(paste("ERROR: Condition column not found:", args$condition))
}

metadata[[args$condition]] <- factor(metadata[[args$condition]])

# Set reference level
if (!is.null(args$reference)) {
    if (args$reference %in% levels(metadata[[args$condition]])) {
        metadata[[args$condition]] <- relevel(metadata[[args$condition]], 
                                              ref = args$reference)
    }
}

cat("   Conditions:", paste(levels(metadata[[args$condition]]), collapse = " vs "), "\n")

# =============================================================================
# Create DGEList and Filter
# =============================================================================

cat("\n⚙️ Creating DGEList...\n")

y <- DGEList(counts = as.matrix(counts), 
             group = metadata[[args$condition]])

# Filter low-expression genes
cat("🔍 Filtering low-count genes...\n")
keep <- filterByExpr(y, min.count = args$`min-count`, 
                    min.total.count = args$`min-total-count`)
y <- y[keep, , keep.lib.sizes = FALSE]

cat("   Before:", nrow(counts), "genes\n")
cat("   After:", sum(keep), "genes\n")

# =============================================================================
# Normalization
# =============================================================================

cat("\n📊 Normalizing (", args$normalization, ")...\n", sep = "")

if (args$normalization != "none") {
    y <- calcNormFactors(y, method = args$normalization)
}

# =============================================================================
# Design Matrix
# =============================================================================

if (!is.null(args$batch) && args$batch %in% colnames(metadata)) {
    design <- model.matrix(~ metadata[[args$batch]] + metadata[[args$condition]])
    cat("   Design with batch correction\n")
} else {
    design <- model.matrix(~ metadata[[args$condition]])
}

# =============================================================================
# Dispersion Estimation
# =============================================================================

cat("\n⚙️ Estimating dispersions...\n")

y <- estimateDisp(y, design, robust = args$robust)

cat("   Common dispersion:", round(y$common.dispersion, 4), "\n")
cat("   BCV:", round(sqrt(y$common.dispersion), 4), "\n")

# =============================================================================
# Differential Expression Testing
# =============================================================================

cat("\n⚙️ Testing for DE (", args$method, ")...\n", sep = "")

if (args$method == "QLF") {
    fit <- glmQLFit(y, design, robust = args$robust)
    results <- glmQLFTest(fit, coef = ncol(design))
} else if (args$method == "LRT") {
    fit <- glmFit(y, design)
    results <- glmLRT(fit, coef = ncol(design))
} else {
    # Exact test (only for simple designs)
    results <- exactTest(y)
}

# =============================================================================
# Extract and Format Results
# =============================================================================

cat("\n📊 Extracting results...\n")

res_table <- topTags(results, n = Inf, adjust.method = "BH", sort.by = "PValue")$table

res_df <- data.frame(
    gene_id = rownames(res_table),
    logFC = res_table$logFC,
    logCPM = res_table$logCPM,
    stringsAsFactors = FALSE
)

# Add test statistic based on method
if (args$method == "QLF") {
    res_df$F <- res_table$F
    res_df$stat <- res_table$F
} else if (args$method == "LRT") {
    res_df$LR <- res_table$LR
    res_df$stat <- res_table$LR
} else {
    res_df$stat <- NA
}

res_df$PValue <- res_table$PValue
res_df$FDR <- res_table$FDR

# Add standardized column names for RAPTOR compatibility
res_df$log2FoldChange <- res_df$logFC
res_df$pvalue <- res_df$PValue
res_df$padj <- res_df$FDR
res_df$baseMean <- 2^res_df$logCPM  # Convert logCPM to linear scale

# Add direction and significance
res_df$direction <- ifelse(res_df$logFC > 0, "up",
                          ifelse(res_df$logFC < 0, "down", "unchanged"))

res_df$is_significant <- res_df$FDR < args$fdr & abs(res_df$logFC) > args$lfc

# Sort by FDR
res_df <- res_df[order(res_df$FDR), ]
rownames(res_df) <- res_df$gene_id

# =============================================================================
# Save Results
# =============================================================================

cat("\n📁 Saving results...\n")

output_dir <- args$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Full results
output_file <- file.path(output_dir, "de_results.csv")
write.csv(res_df, output_file, row.names = FALSE)
cat("   ✓ de_results.csv\n")

# Significant genes
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
    pipeline = "edgeR",
    edger_version = as.character(packageVersion("edgeR")),
    parameters = list(
        fdr_threshold = args$fdr,
        lfc_threshold = args$lfc,
        normalization = args$normalization,
        test_method = args$method,
        robust = args$robust,
        min_count = args$`min-count`
    ),
    data = list(
        n_samples = ncol(y$counts),
        n_genes_input = nrow(counts),
        n_genes_tested = nrow(y$counts),
        common_dispersion = round(y$common.dispersion, 4),
        bcv = round(sqrt(y$common.dispersion), 4)
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
    
    # BCV plot
    plotBCV(y, main = "Biological Coefficient of Variation")
    
    # MD plot
    plotMD(results, main = "edgeR MD Plot")
    abline(h = c(-args$lfc, args$lfc), col = "blue", lty = 2)
    
    # MDS plot
    plotMDS(y, col = as.numeric(y$samples$group), 
           main = "MDS Plot")
    
    dev.off()
    cat("   ✓ de_plots.pdf\n")
}

# =============================================================================
# Final Summary
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  ✅ edgeR ANALYSIS COMPLETE!\n")
cat("═══════════════════════════════════════════════════════════════\n")

cat("\n  📊 Summary:\n")
cat("     • Total genes tested:", nrow(res_df), "\n")
cat("     • Significant genes:", nrow(sig_df), 
    "(", round(100 * nrow(sig_df) / nrow(res_df), 1), "%)\n")
cat("     • Upregulated:", n_up, "\n")
cat("     • Downregulated:", n_down, "\n")

cat("\n  📁 Output Directory:", output_dir, "\n")

cat("\n  🔜 Next Steps:\n")
cat("     raptor import-de --de-file", output_file, "\n")
cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Making free science for everybody around the world 🌍\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
