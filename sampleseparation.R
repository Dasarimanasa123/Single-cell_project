library(Matrix)
library(data.table)
library(dplyr)

# Paths to your input files
mtx_file <- "GSE178318_matrix.mtx"
barcodes_file <- "GSE178318_barcodes.tsv"
features_file <- "GSE178318_genes.tsv"

# Load the Matrix Market file (.mtx)
combined_matrix <- readMM("GSE178318_matrix.mtx.gz")

# Load barcodes and features
barcodes <- fread("GSE178318_barcodes.tsv.gz", header = FALSE)
features <- fread("GSE178318_genes.tsv.gz", header = FALSE)

# Add a sample ID column to barcodes (modify based on your barcode structure)
# Example: If sample IDs are suffixes after a hyphen
barcodes[, sample_id := tstrsplit(V1, "-", keep = 2)]
barcodes[, sample_id := tstrsplit(V1, "_", keep = 2)]

# Get unique sample IDs
unique_samples <- unique(barcodes$sample_id)

# Create a directory to save output
output_dir <- "separated_samples"
dir.create(output_dir, showWarnings = FALSE)

# Split and save the matrix for each sample
for (sample in unique_samples) {
  cat("Processing sample:", sample, "\n")
  
  # Get indices of barcodes belonging to this sample
  sample_indices <- which(barcodes$sample_id == sample)
  
  # Subset the matrix
  sample_matrix <- combined_matrix[, sample_indices]
  
  # Save the sample-specific matrix, barcodes, and features
  sample_dir <- file.path(output_dir, sample)
  dir.create(sample_dir, showWarnings = FALSE)
  
  # Write the sample-specific matrix
  writeMM(sample_matrix, file.path(sample_dir, "matrix.mtx"))
  
  # Write the barcodes and features
  fwrite(barcodes[sample_indices, .(V1)], file.path(sample_dir, "barcodes.tsv"), col.names = FALSE, sep = "\t")
  fwrite(features, file.path(sample_dir, "features.tsv"), col.names = FALSE, sep = "\t")
  
  cat("Saved data for sample:", sample, "in", sample_dir, "\n")
}
