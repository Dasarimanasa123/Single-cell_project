library(Seurat)
library(tidyverse)

data <- readRDS("GSE158520_M1_Tumor.rds")
data_2 <- readRDS("GSE158520_M4_Tumor.rds")

# Load Matrix package
library(Matrix)

# Extract the raw count matrix
count_matrix <- data_2@assays$RNA@counts

# Save the count matrix as a sparse matrix (if it's not already)
sparse_matrix <- as(count_matrix, "CsparseMatrix")

# Write the sparse matrix to a .mtx file
writeMM(sparse_matrix, "count_matrix.mtx")

# Save the barcodes (cells) to a file
write.csv(colnames(count_matrix), "barcodes.tsv", row.names = FALSE)

# Save the feature names (genes) to a file
write.csv(rownames(count_matrix), "features.tsv", row.names = FALSE)

COL12 <- ReadMtx(
  mtx = "GSM4802190/count_matrix.mtx", 
  features = "GSM4802190/features.tsv", 
  cells = "GSM4802190/barcodes.tsv", 
  feature.column = 1  # Specify the correct column for features
)

