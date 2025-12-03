#!/usr/bin/env Rscript

###############################################################################
# 04_visualize_comparison.R
# Create visualizations comparing pipeline results
#
# Usage: Rscript 04_visualize_comparison.R <results_dir>
#
# Author: Ayeh Bolouki
# Organization: RAPTOR Project
# License: MIT
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggVennDiagram)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript 04_visualize_comparison.R <results_dir>\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  results_dir : Directory containing comparison results\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript 04_visualize_comparison.R results/benchmark\n")
  quit(status = 1)
}

results_dir <- args[1]

cat("=== RAPTOR Pipeline Visualization ===\n")
cat(sprintf("Results directory: %s\n\n", results_dir))

# Create plots directory
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Check for required files
summary_file <- file.path(results_dir, "comparison_summary.csv")
overlap_file <- file.path(results_dir, "overlap_matrix.csv")

if (!file.exists(summary_file)) {
  cat("Error: comparison_summary.csv not found\n")
  cat("Please run 03_compare_results.R first\n")
  quit(status = 1)
}

# Load data
cat("Loading comparison data...\n")
comparison <- read.csv(summary_file)
overlap_matrix <- as.matrix(read.csv(overlap_file, row.names = 1))

# Check for performance metrics
perf_file <- file.path(results_dir, "performance_metrics.csv")
has_performance <- file.exists(perf_file)
if (has_performance) {
  performance <- read.csv(perf_file)
  cat("  Performance metrics found\n")
}

cat(sprintf("  %d pipelines found\n\n", nrow(comparison)))

# 1. Bar plot: Number of significant genes per pipeline
cat("Creating visualizations...\n")
cat("  1. Bar plot - Significant genes per pipeline\n")

p1 <- ggplot(comparison, aes(x = reorder(pipeline, n_significant), 
                             y = n_significant)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n_significant), hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "Number of Significant Genes per Pipeline",
       x = "Pipeline",
       y = "Number of Significant Genes (padj < 0.05)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10))

ggsave(file.path(plots_dir, "01_significant_genes_barplot.pdf"), 
       p1, width = 10, height = 6)
ggsave(file.path(plots_dir, "01_significant_genes_barplot.png"), 
       p1, width = 10, height = 6, dpi = 300)

# 2. Overlap heatmap
cat("  2. Heatmap - Pairwise overlaps\n")

# Calculate Jaccard index for heatmap
n_pipelines <- nrow(overlap_matrix)
jaccard_matrix <- matrix(0, nrow = n_pipelines, ncol = n_pipelines)
rownames(jaccard_matrix) <- rownames(overlap_matrix)
colnames(jaccard_matrix) <- colnames(overlap_matrix)

for (i in 1:n_pipelines) {
  for (j in 1:n_pipelines) {
    if (i == j) {
      jaccard_matrix[i, j] <- 1
    } else {
      # Union size = size_i + size_j - intersection
      size_i <- comparison$n_significant[comparison$pipeline == rownames(overlap_matrix)[i]]
      size_j <- comparison$n_significant[comparison$pipeline == colnames(overlap_matrix)[j]]
      intersection <- overlap_matrix[i, j]
      union_size <- size_i + size_j - intersection
      jaccard_matrix[i, j] <- ifelse(union_size > 0, intersection / union_size, 0)
    }
  }
}

pdf(file.path(plots_dir, "02_overlap_heatmap.pdf"), width = 10, height = 9)
pheatmap(jaccard_matrix,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
         main = "Pipeline Similarity (Jaccard Index)",
         fontsize = 10,
         fontsize_number = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         angle_col = 45)
dev.off()

png(file.path(plots_dir, "02_overlap_heatmap.png"), 
    width = 10, height = 9, units = "in", res = 300)
pheatmap(jaccard_matrix,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
         main = "Pipeline Similarity (Jaccard Index)",
         fontsize = 10,
         fontsize_number = 8,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         angle_col = 45)
dev.off()

# 3. Venn diagram (for top 5 pipelines by number of DEGs)
if (nrow(comparison) >= 2 && nrow(comparison) <= 5) {
  cat("  3. Venn diagram - Overlap between pipelines\n")
  
  # Load actual gene lists
  gene_lists <- list()
  
  for (pipeline in comparison$pipeline) {
    # Find result files
    pipeline_dir <- file.path(results_dir, pipeline)
    result_files <- list.files(pipeline_dir, 
                               pattern = "*significant*.csv|*results*.csv",
                               recursive = TRUE, full.names = TRUE)
    
    if (length(result_files) > 0) {
      tryCatch({
        df <- read.csv(result_files[1])
        # Get gene IDs
        if ("gene_id" %in% names(df)) {
          genes <- df$gene_id
        } else if ("X" %in% names(df)) {
          genes <- df$X
        } else {
          genes <- rownames(df)
        }
        gene_lists[[pipeline]] <- genes
      }, error = function(e) {
        # Skip if error
      })
    }
  }
  
  if (length(gene_lists) >= 2) {
    # Take top 5 for visualization
    if (length(gene_lists) > 5) {
      top_pipelines <- comparison %>%
        arrange(desc(n_significant)) %>%
        slice(1:5) %>%
        pull(pipeline)
      gene_lists <- gene_lists[top_pipelines]
    }
    
    tryCatch({
      p3 <- ggVennDiagram(gene_lists, 
                         label = "count",
                         label_alpha = 0) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = "Overlap of Significant Genes Between Pipelines") +
        theme(plot.title = element_text(face = "bold", size = 14))
      
      ggsave(file.path(plots_dir, "03_venn_diagram.pdf"), 
             p3, width = 10, height = 8)
      ggsave(file.path(plots_dir, "03_venn_diagram.png"), 
             p3, width = 10, height = 8, dpi = 300)
    }, error = function(e) {
      cat("    Warning: Could not create Venn diagram\n")
    })
  }
}

# 4. Performance comparison (if truth set was provided)
if (has_performance) {
  cat("  4. Performance comparison plots\n")
  
  # Precision-Recall plot
  p4a <- ggplot(performance, 
                aes(x = recall, y = precision, label = pipeline)) +
    geom_point(size = 4, color = "steelblue") +
    geom_text(hjust = -0.1, vjust = 0.5, size = 3, angle = 0) +
    xlim(0, 1) + ylim(0, 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "gray") +
    labs(title = "Pipeline Performance: Precision vs Recall",
         x = "Recall (Sensitivity)",
         y = "Precision") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave(file.path(plots_dir, "04a_precision_recall.pdf"), 
         p4a, width = 10, height = 8)
  ggsave(file.path(plots_dir, "04a_precision_recall.png"), 
         p4a, width = 10, height = 8, dpi = 300)
  
  # F1 score comparison
  p4b <- ggplot(performance, 
                aes(x = reorder(pipeline, f1_score), y = f1_score)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3f", f1_score)), 
              hjust = -0.2, size = 3.5) +
    coord_flip() +
    ylim(0, 1) +
    labs(title = "Pipeline Performance: F1 Score",
         x = "Pipeline",
         y = "F1 Score") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 10))
  
  ggsave(file.path(plots_dir, "04b_f1_scores.pdf"), 
         p4b, width = 10, height = 6)
  ggsave(file.path(plots_dir, "04b_f1_scores.png"), 
         p4b, width = 10, height = 6, dpi = 300)
  
  # Confusion matrix visualization
  performance_long <- performance %>%
    select(pipeline, TP, FP, FN) %>%
    pivot_longer(cols = c(TP, FP, FN), 
                names_to = "metric", 
                values_to = "count")
  
  p4c <- ggplot(performance_long, 
                aes(x = pipeline, y = count, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("TP" = "darkgreen", 
                                 "FP" = "orange", 
                                 "FN" = "red"),
                     labels = c("True Positives", 
                               "False Positives", 
                               "False Negatives")) +
    labs(title = "Pipeline Performance: Confusion Matrix Components",
         x = "Pipeline",
         y = "Count",
         fill = "Metric") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  
  ggsave(file.path(plots_dir, "04c_confusion_components.pdf"), 
         p4c, width = 12, height = 8)
  ggsave(file.path(plots_dir, "04c_confusion_components.png"), 
         p4c, width = 12, height = 8, dpi = 300)
}

# 5. Summary comparison plot
cat("  5. Summary comparison plot\n")

# Normalize n_significant to 0-1 scale for visualization
comparison_viz <- comparison %>%
  mutate(norm_significant = n_significant / max(n_significant))

if (has_performance) {
  comparison_viz <- comparison_viz %>%
    left_join(performance %>% select(pipeline, precision, recall, f1_score),
              by = "pipeline")
  
  comparison_long <- comparison_viz %>%
    select(pipeline, norm_significant, precision, recall, f1_score) %>%
    pivot_longer(cols = -pipeline, 
                names_to = "metric", 
                values_to = "value")
  
  metric_labels <- c(
    "norm_significant" = "Normalized\nDEG Count",
    "precision" = "Precision",
    "recall" = "Recall",
    "f1_score" = "F1 Score"
  )
  
  p5 <- ggplot(comparison_long, 
               aes(x = metric, y = value, group = pipeline, color = pipeline)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_discrete(labels = metric_labels) +
    ylim(0, 1) +
    labs(title = "Pipeline Performance Summary",
         x = "",
         y = "Value (normalized)",
         color = "Pipeline") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 10),
          legend.position = "right")
  
  ggsave(file.path(plots_dir, "05_summary_comparison.pdf"), 
         p5, width = 12, height = 8)
  ggsave(file.path(plots_dir, "05_summary_comparison.png"), 
         p5, width = 12, height = 8, dpi = 300)
}

cat("\n=== Visualization Complete! ===\n")
cat(sprintf("All plots saved to: %s\n\n", plots_dir))

cat("Generated plots:\n")
cat("  1. 01_significant_genes_barplot - Number of DEGs per pipeline\n")
cat("  2. 02_overlap_heatmap - Jaccard similarity between pipelines\n")
if (file.exists(file.path(plots_dir, "03_venn_diagram.pdf"))) {
  cat("  3. 03_venn_diagram - Overlap visualization\n")
}
if (has_performance) {
  cat("  4a. 04a_precision_recall - Precision vs Recall\n")
  cat("  4b. 04b_f1_scores - F1 score comparison\n")
  cat("  4c. 04c_confusion_components - TP/FP/FN breakdown\n")
  cat("  5. 05_summary_comparison - Overall performance radar\n")
}

cat("\n✓ All visualizations complete!\n")
