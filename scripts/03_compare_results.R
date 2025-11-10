#!/usr/bin/env Rscript

###############################################################################
# 03_compare_results.R
# Compare differential expression results from multiple pipelines
#
# Usage: Rscript 03_compare_results.R <results_dir> [--truth truth_set.csv]
#
# Author: Ayeh Bolouki
# Organization: RAPTOR Project
# License: MIT
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript 03_compare_results.R <results_dir> [--truth truth_set.csv]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  results_dir  : Directory containing pipeline output folders\n")
  cat("  --truth      : Optional ground truth file for accuracy assessment\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript 03_compare_results.R results/benchmark\n")
  cat("  Rscript 03_compare_results.R results/benchmark --truth simulated_data/truth_set.csv\n")
  quit(status = 1)
}

results_dir <- args[1]
truth_file <- NULL

# Check for truth file
if (length(args) >= 3 && args[2] == "--truth") {
  truth_file <- args[3]
}

cat("=== RAPTOR Pipeline Comparison ===\n")
cat(sprintf("Results directory: %s\n", results_dir))
if (!is.null(truth_file)) {
  cat(sprintf("Truth set: %s\n", truth_file))
}
cat("\n")

# Find all pipeline result directories
pipeline_dirs <- list.dirs(results_dir, recursive = FALSE)
pipeline_names <- basename(pipeline_dirs)

cat(sprintf("Found %d pipeline result directories\n", length(pipeline_dirs)))
cat("\n")

# Function to find and load DEG results from a pipeline
load_deg_results <- function(pipeline_dir) {
  # Common file patterns for DEG results
  patterns <- c(
    "**/DEG_results*.csv",
    "**/results*.csv",
    "**/differential_expression/*significant*.csv",
    "**/de_genes*.csv"
  )
  
  for (pattern in patterns) {
    files <- list.files(pipeline_dir, pattern = pattern, 
                       recursive = TRUE, full.names = TRUE)
    
    # Prefer files with "significant" in name
    sig_files <- grep("significant", files, value = TRUE)
    if (length(sig_files) > 0) {
      files <- sig_files
    }
    
    if (length(files) > 0) {
      # Take first file found
      tryCatch({
        df <- read.csv(files[1])
        return(df)
      }, error = function(e) {
        return(NULL)
      })
    }
  }
  
  return(NULL)
}

# Load results from all pipelines
cat("Loading pipeline results...\n")
all_results <- list()

for (i in seq_along(pipeline_dirs)) {
  pipeline_name <- pipeline_names[i]
  cat(sprintf("  %d. %s... ", i, pipeline_name))
  
  results <- load_deg_results(pipeline_dirs[i])
  
  if (!is.null(results)) {
    all_results[[pipeline_name]] <- results
    cat(sprintf("✓ %d genes\n", nrow(results)))
  } else {
    cat("✗ No results found\n")
  }
}

cat(sprintf("\nSuccessfully loaded results from %d pipelines\n\n", 
            length(all_results)))

if (length(all_results) == 0) {
  cat("Error: No pipeline results could be loaded\n")
  quit(status = 1)
}

# Standardize column names across pipelines
standardize_columns <- function(df) {
  df_cols <- tolower(names(df))
  
  # Gene ID
  gene_col <- which(df_cols %in% c("gene_id", "geneid", "gene", "id", "x"))
  if (length(gene_col) > 0) {
    names(df)[gene_col[1]] <- "gene_id"
  } else if ("X" %in% names(df)) {
    names(df)[names(df) == "X"] <- "gene_id"
  } else {
    df$gene_id <- rownames(df)
  }
  
  # Log fold change
  lfc_col <- which(df_cols %in% c("log2foldchange", "logfc", "log2fc", "lfc"))
  if (length(lfc_col) > 0) {
    names(df)[lfc_col[1]] <- "log2fc"
  }
  
  # P-value
  pval_col <- which(df_cols %in% c("pvalue", "pval", "p.value"))
  if (length(pval_col) > 0) {
    names(df)[pval_col[1]] <- "pvalue"
  }
  
  # Adjusted p-value
  padj_col <- which(df_cols %in% c("padj", "fdr", "adj.p.value", "p.adj", "qvalue"))
  if (length(padj_col) > 0) {
    names(df)[padj_col[1]] <- "padj"
  }
  
  return(df)
}

# Standardize all results
all_results <- map(all_results, standardize_columns)

# Extract significant genes from each pipeline
get_sig_genes <- function(df, padj_threshold = 0.05) {
  if ("padj" %in% names(df)) {
    sig <- df %>%
      filter(!is.na(padj) & padj < padj_threshold) %>%
      pull(gene_id)
    return(sig)
  }
  return(character(0))
}

sig_genes <- map(all_results, get_sig_genes)

# Summary statistics
cat("=== Significant Genes per Pipeline ===\n")
for (name in names(sig_genes)) {
  cat(sprintf("  %-35s: %5d genes\n", name, length(sig_genes[[name]])))
}
cat("\n")

# Pairwise comparisons
cat("=== Pairwise Overlaps ===\n")
pipeline_names_vec <- names(sig_genes)
n_pipelines <- length(pipeline_names_vec)

overlap_matrix <- matrix(0, nrow = n_pipelines, ncol = n_pipelines,
                        dimnames = list(pipeline_names_vec, pipeline_names_vec))

for (i in 1:n_pipelines) {
  for (j in 1:n_pipelines) {
    if (i <= j) {
      overlap <- length(intersect(sig_genes[[i]], sig_genes[[j]]))
      overlap_matrix[i, j] <- overlap
      overlap_matrix[j, i] <- overlap
      
      if (i < j) {
        jaccard <- overlap / length(union(sig_genes[[i]], sig_genes[[j]]))
        cat(sprintf("  %s vs %s: %d genes (Jaccard: %.3f)\n",
                   pipeline_names_vec[i], pipeline_names_vec[j],
                   overlap, jaccard))
      }
    }
  }
}
cat("\n")

# Core consensus genes (found in majority of pipelines)
all_genes_list <- sig_genes[sapply(sig_genes, length) > 0]
if (length(all_genes_list) > 0) {
  gene_counts <- table(unlist(all_genes_list))
  threshold <- ceiling(length(all_genes_list) / 2)
  
  consensus_genes <- names(gene_counts[gene_counts >= threshold])
  cat(sprintf("Consensus genes (in >= %d pipelines): %d\n", 
              threshold, length(consensus_genes)))
  cat("\n")
}

# Compare to ground truth if provided
if (!is.null(truth_file) && file.exists(truth_file)) {
  cat("=== Comparison to Ground Truth ===\n")
  
  truth <- read.csv(truth_file)
  true_de <- truth %>%
    filter(is_de == TRUE) %>%
    pull(gene_id)
  
  cat(sprintf("True DE genes: %d\n\n", length(true_de)))
  
  # Calculate performance metrics for each pipeline
  performance <- data.frame(
    pipeline = character(),
    TP = numeric(),
    FP = numeric(),
    FN = numeric(),
    precision = numeric(),
    recall = numeric(),
    f1_score = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(sig_genes)) {
    predicted <- sig_genes[[name]]
    
    tp <- length(intersect(predicted, true_de))
    fp <- length(setdiff(predicted, true_de))
    fn <- length(setdiff(true_de, predicted))
    
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1 <- ifelse(precision + recall > 0, 
                2 * (precision * recall) / (precision + recall), 0)
    
    performance <- rbind(performance, data.frame(
      pipeline = name,
      TP = tp,
      FP = fp,
      FN = fn,
      precision = precision,
      recall = recall,
      f1_score = f1
    ))
  }
  
  # Sort by F1 score
  performance <- performance %>%
    arrange(desc(f1_score))
  
  cat("Performance Metrics:\n")
  print(performance, row.names = FALSE)
  cat("\n")
  
  # Save performance metrics
  perf_file <- file.path(results_dir, "performance_metrics.csv")
  write.csv(performance, perf_file, row.names = FALSE)
  cat(sprintf("✓ Performance metrics saved: %s\n\n", perf_file))
}

# Create comparison summary
comparison_summary <- data.frame(
  pipeline = names(sig_genes),
  n_significant = sapply(sig_genes, length),
  stringsAsFactors = FALSE
)

# Save comparison results
summary_file <- file.path(results_dir, "comparison_summary.csv")
write.csv(comparison_summary, summary_file, row.names = FALSE)

overlap_file <- file.path(results_dir, "overlap_matrix.csv")
write.csv(overlap_matrix, overlap_file)

cat("=== Output Files ===\n")
cat(sprintf("  Comparison summary: %s\n", summary_file))
cat(sprintf("  Overlap matrix: %s\n", overlap_file))
if (!is.null(truth_file) && file.exists(truth_file)) {
  cat(sprintf("  Performance metrics: %s\n", 
              file.path(results_dir, "performance_metrics.csv")))
}

cat("\n✓ Pipeline comparison complete!\n")
cat("\nNext step: Rscript scripts/04_visualize_comparison.R\n")
