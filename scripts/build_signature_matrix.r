# scripts/build_signature_matrix.r
# Description: Builds a DWLS signature matrix from a Scanpy .h5ad snapshot.
# Dependencies: Seurat, DWLS, zellkonverter, SingleCellExperiment, dotenv

library("DWLS")
library("zellkonverter")
library("SingleCellExperiment")
library("Seurat")
library("dotenv")

# 1. Setup Environment
# Load environment variables (assumes .env file exists in project root)
load_dot_env()
data_folder <- Sys.getenv("DATA_FOLDER")

# Define paths
input_file <- file.path(data_folder, "sc_deconv_snapshot.h5ad")
output_path <- file.path(data_folder, "dwls_results")

# Ensure output directory exists
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# 2. Load and Convert Data
# Reading h5ad directly into SingleCellExperiment
message("Reading h5ad file: ", input_file)
sce <- readH5AD(input_file, use_hdf5 = TRUE)

# Convert to Seurat object
# Assuming 'X' in h5ad contains raw counts (based on preproc_sc.qmd)
message("Converting to Seurat object...")
seu <- as.Seurat(sce, counts = "X", data = NULL)

# 3. Preprocessing
# DWLS requires normalized data for DE analysis.
# If 'X' was raw counts, we must normalize.
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Set the identity class to the cell type column
# Ensure 'celltype' exists in metadata (imported from .obs)
if (!"celltype" %in% colnames(seu@meta.data)) {
  stop("Column 'celltype' not found in metadata.")
}
Idents(seu) <- "celltype"

# 4. Build Signature Matrix
message("Building DWLS signature matrix...")

# Arguments:
# scdata: The Seurat object
# id: Vector of cell type labels corresponding to cells (required by DWLS)
# path: Output directory
# diff.cutoff: LogFC cutoff (default 0.5)
# pval.cutoff: P-value cutoff (default 0.01)
# f: Max number of genes per group (default 200)

buildSignatureMatrixUsingSeurat(
  scdata = seu,
  id = as.character(Idents(seu)),
  path = output_path,
  diff.cutoff = 0.5,
  pval.cutoff = 0.01
)

message("Signature matrix generation complete. Results saved to: ", output_path)