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

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
output_dir <- normalizePath(output_path) # Convert to full absolute path

# 2. Load and Convert Data
message("Reading h5ad file: ", input_file)
sce <- readH5AD(input_file, use_hdf5 = TRUE)

message("Converting to Seurat object...")
seu <- as.Seurat(sce, counts = "X", data = NULL)

# 3. Preprocessing
# DWLS with MAST works best on Log-Normalized data
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Ensure celltype annotation exists
if (!"celltype" %in% colnames(seu@meta.data)) {
  stop("Column 'celltype' not found in metadata.")
}
Idents(seu) <- "celltype"

# 4. Prepare Data for DWLS Core Function
# We extract the normalized matrix and ID vector manually to bypass the version conflict
message("Extracting data matrix and identities...")

# Extract sparse matrix (using 'data' slot for normalized counts)
# Note: For very large datasets, this conversion to a regular matrix might be memory intensive.
# DWLS/MAST often requires a dense matrix or specific sparse handling. 
# We convert to a standard matrix to ensure compatibility with older DWLS code.
data_matrix <- as.matrix(GetAssayData(seu, layer = "data"))

# Extract identities
cell_ids <- as.character(Idents(seu))

# 5. Build Signature Matrix
message("Building DWLS signature matrix using MAST...")
message("Intermediate results will be saved to: ", output_dir)

# Capture the returned matrix object
sig_mat <- buildSignatureMatrixMAST(
  scdata = data_matrix,
  id = cell_ids,
  path = output_dir,   # DWLS uses this for DE stats files
  diff.cutoff = 0.5,
  pval.cutoff = 0.01,
  f = 250              # make it comparable in gene size to BulkTrajBlend
)

# 6. Save Final Output
# Explicitly write the result to a file
final_output_file <- file.path(output_dir, "signature_matrix.txt")
message("Writing final signature matrix to: ", final_output_file)

write.table(
  sig_mat, 
  file = final_output_file, 
  sep = "\t", 
  quote = FALSE, 
  col.names = NA
)

message("Done.")
