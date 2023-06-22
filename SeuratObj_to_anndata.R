setwd("<path/to/dir>")

# Read the Seurat object generated from coreSC
library(Seurat)
library(scrattch.io)

# TODO address scrattch.io dep

seurat_obj <- readRDS("<path/to/.rds/file>", refhook = NULL)

# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[, 1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[, 2]

# Take a look at the metadata file
metadata <- data.frame(seurat_obj@meta.data)
metadata$object <- gsub(".*X", "", metadata$object)
metadata$cluster <- metadata$seurat_clusters

write.csv(
  metadata,
  file = "metadata.csv",
  quote = FALSE,
  row.names = FALSE
)

# write expression counts matrix
counts_matrix <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  slot = "counts"
)

# write counts matrix to file in memory considerate maner
write_dgCMatrix_csv(
  counts_matrix,
  "counts.mtx",
  col1_name = "",
  chunk_size = 1000
)

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(
  seurat_obj@reductions$pca@cell.embeddings,
  file = "pca.csv",
  quote = FALSE,
  row.names = FALSE
)

# write gene names
write.table(
  data.frame("gene" = rownames(counts_matrix)),
  file = "gene_names.csv",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
