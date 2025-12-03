#!/usr/bin/env Rscript

###############################################################################
# 00_simulate_data.R
# Simulate RNA-seq count data for pipeline testing and benchmarking
#
# Usage: Rscript 00_simulate_data.R [options]
#
# Author: Ayeh Bolouki
# Organization: RAPTOR Project
# License: MIT
###############################################################################

suppressPackageStartupMessages({
  library(polyester)
  library(Biostrings)
  library(rtracklayer)
  library(optparse)
})

# Command line options
option_list <- list(
  make_option(c("-o", "--output"), type="character", default="simulated_data",
              help="Output directory [default: %default]"),
  make_option(c("-n", "--nsamples"), type="integer", default=6,
              help="Number of samples (must be even) [default: %default]"),
  make_option(c("-g", "--ngenes"), type="integer", default=10000,
              help="Number of genes to simulate [default: %default]"),
  make_option(c("-r", "--reads"), type="integer", default=1000000,
              help="Number of reads per sample [default: %default]"),
  make_option(c("-d", "--ndiff"), type="integer", default=1000,
              help="Number of differentially expressed genes [default: %default]"),
  make_option(c("-f", "--foldchange"), type="numeric", default=3,
              help="Fold change for DE genes [default: %default]"),
  make_option(c("-s", "--seed"), type="integer", default=42,
              help="Random seed for reproducibility [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Simulate RNA-seq data for RAPTOR pipeline testing")
opt <- parse_args(opt_parser)

# Validate inputs
if (opt$nsamples %% 2 != 0) {
  stop("Number of samples must be even (equal groups for case/control)")
}

# Set seed for reproducibility
set.seed(opt$seed)

cat("=== RAPTOR RNA-seq Data Simulation ===\n")
cat(sprintf("Output directory: %s\n", opt$output))
cat(sprintf("Number of samples: %d\n", opt$nsamples))
cat(sprintf("Number of genes: %d\n", opt$ngenes))
cat(sprintf("Reads per sample: %d\n", opt$reads))
cat(sprintf("DE genes: %d\n", opt$ndiff))
cat(sprintf("Fold change: %.1f\n", opt$foldchange))
cat(sprintf("Random seed: %d\n\n", opt$seed))

# Create output directory
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

# Generate transcriptome
cat("Step 1: Generating synthetic transcriptome...\n")
transcript_length <- 1000
transcripts <- DNAStringSet(replicate(opt$ngenes, {
  paste(sample(c("A", "C", "G", "T"), transcript_length, replace = TRUE), 
        collapse = "")
}))
names(transcripts) <- paste0("GENE", sprintf("%05d", 1:opt$ngenes))

# Save transcriptome
fasta_file <- file.path(opt$output, "transcriptome.fa")
writeXStringSet(transcripts, fasta_file)
cat(sprintf("  Saved transcriptome: %s\n", fasta_file))

# Create GTF annotation
cat("Step 2: Creating annotation file...\n")
gtf_data <- data.frame(
  seqname = names(transcripts),
  source = "simulation",
  feature = "exon",
  start = 1,
  end = transcript_length,
  score = ".",
  strand = "+",
  frame = ".",
  attribute = paste0('gene_id "', names(transcripts), '"; transcript_id "', 
                    names(transcripts), '_1";')
)

gtf_file <- file.path(opt$output, "annotation.gtf")
write.table(gtf_data, gtf_file, quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)
cat(sprintf("  Saved annotation: %s\n", gtf_file))

# Setup simulation parameters
cat("Step 3: Setting up expression parameters...\n")
n_per_group <- opt$nsamples / 2

# Baseline expression levels (negative binomial parameters)
baseline_mean <- rnorm(opt$ngenes, mean = 100, sd = 50)
baseline_mean[baseline_mean < 10] <- 10  # Minimum expression

# Identify DE genes
de_genes <- sample(1:opt$ngenes, opt$ndiff)
n_up <- floor(opt$ndiff / 2)
n_down <- opt$ndiff - n_up
up_genes <- de_genes[1:n_up]
down_genes <- de_genes[(n_up + 1):opt$ndiff]

# Create fold change matrix
fold_changes <- matrix(1, nrow = opt$ngenes, ncol = opt$nsamples)
# Upregulated in treatment group
fold_changes[up_genes, (n_per_group + 1):opt$nsamples] <- opt$foldchange
# Downregulated in treatment group  
fold_changes[down_genes, (n_per_group + 1):opt$nsamples] <- 1/opt$foldchange

cat(sprintf("  DE genes: %d upregulated, %d downregulated\n", n_up, n_down))

# Simulate reads
cat("Step 4: Simulating RNA-seq reads...\n")
cat("  This may take a few minutes...\n")

simulate_experiment_countmat(
  fasta = fasta_file,
  gtf = gtf_file,
  seqpath = opt$output,
  outdir = file.path(opt$output, "reads"),
  num_reps = c(n_per_group, n_per_group),
  reads_per_transcript = baseline_mean,
  fold_changes = fold_changes,
  readlen = 100,
  paired = TRUE,
  seed = opt$seed
)

# Rename files to meaningful names
cat("Step 5: Organizing output files...\n")
reads_dir <- file.path(opt$output, "reads")

for (i in 1:n_per_group) {
  file.rename(
    file.path(reads_dir, paste0("sample_", sprintf("%02d", i), "_1.fasta")),
    file.path(reads_dir, paste0("control_", i, "_R1.fastq"))
  )
  file.rename(
    file.path(reads_dir, paste0("sample_", sprintf("%02d", i), "_2.fasta")),
    file.path(reads_dir, paste0("control_", i, "_R2.fastq"))
  )
}

for (i in 1:n_per_group) {
  j <- i + n_per_group
  file.rename(
    file.path(reads_dir, paste0("sample_", sprintf("%02d", j), "_1.fasta")),
    file.path(reads_dir, paste0("treatment_", i, "_R1.fastq"))
  )
  file.rename(
    file.path(reads_dir, paste0("sample_", sprintf("%02d", j), "_2.fasta")),
    file.path(reads_dir, paste0("treatment_", i, "_R2.fastq"))
  )
}

# Compress FASTQ files
cat("Step 6: Compressing FASTQ files...\n")
fastq_files <- list.files(reads_dir, pattern = "\\.fastq$", full.names = TRUE)
for (f in fastq_files) {
  system(paste("gzip", f))
}

# Create sample metadata
cat("Step 7: Creating sample metadata...\n")
metadata <- data.frame(
  sample_id = c(paste0("control_", 1:n_per_group), 
                paste0("treatment_", 1:n_per_group)),
  condition = rep(c("control", "treatment"), each = n_per_group),
  batch = rep("batch1", opt$nsamples),
  replicate = rep(1:n_per_group, 2),
  fastq_1 = c(paste0("control_", 1:n_per_group, "_R1.fastq.gz"),
              paste0("treatment_", 1:n_per_group, "_R1.fastq.gz")),
  fastq_2 = c(paste0("control_", 1:n_per_group, "_R2.fastq.gz"),
              paste0("treatment_", 1:n_per_group, "_R2.fastq.gz"))
)

metadata_file <- file.path(opt$output, "sample_metadata.csv")
write.csv(metadata, metadata_file, row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved metadata: %s\n", metadata_file))

# Save truth set (DE genes)
cat("Step 8: Saving ground truth information...\n")
truth_data <- data.frame(
  gene_id = names(transcripts),
  is_de = FALSE,
  direction = "none",
  true_fc = 1
)

truth_data$is_de[up_genes] <- TRUE
truth_data$direction[up_genes] <- "up"
truth_data$true_fc[up_genes] <- opt$foldchange

truth_data$is_de[down_genes] <- TRUE
truth_data$direction[down_genes] <- "down"
truth_data$true_fc[down_genes] <- 1/opt$foldchange

truth_file <- file.path(opt$output, "truth_set.csv")
write.csv(truth_data, truth_file, row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved ground truth: %s\n", truth_file))

# Create summary report
cat("Step 9: Creating summary report...\n")
summary_file <- file.path(opt$output, "simulation_summary.txt")
sink(summary_file)
cat("=== RAPTOR RNA-seq Simulation Summary ===\n\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Random seed: %d\n\n", opt$seed))
cat("Parameters:\n")
cat(sprintf("  Samples: %d (%d control, %d treatment)\n", 
            opt$nsamples, n_per_group, n_per_group))
cat(sprintf("  Genes: %d\n", opt$ngenes))
cat(sprintf("  DE genes: %d (%d up, %d down)\n", opt$ndiff, n_up, n_down))
cat(sprintf("  Fold change: %.1f\n", opt$foldchange))
cat(sprintf("  Reads per sample: %d\n", opt$reads))
cat(sprintf("  Read length: 100 bp paired-end\n\n"))
cat("Output files:\n")
cat(sprintf("  Transcriptome: %s\n", fasta_file))
cat(sprintf("  Annotation: %s\n", gtf_file))
cat(sprintf("  Reads: %s/*.fastq.gz\n", reads_dir))
cat(sprintf("  Metadata: %s\n", metadata_file))
cat(sprintf("  Ground truth: %s\n", truth_file))
sink()

cat(sprintf("\n=== Simulation Complete! ===\n"))
cat(sprintf("Output directory: %s\n", opt$output))
cat(sprintf("Summary report: %s\n\n", summary_file))

cat("Next steps:\n")
cat("  1. Use simulated FASTQ files to test pipelines\n")
cat("  2. Compare results to ground truth (truth_set.csv)\n")
cat("  3. Evaluate pipeline performance\n\n")
