


# Load required libraries
library(Seurat)
library(loupeR)


# Load merged object (after SCTransform integration and clustering)
merged_obj <- readRDS("C:/scRNAseq/Shiny_App_Project/outputs/merged_scRNA-Seq_Final.rds")


# OPTIONAL: Restore RNA assay counts if needed
merged_obj[["RNA"]] <- CreateAssayObject(
  counts = GetAssayData(merged_obj, assay = "RNA", slot = "counts")
)


# Run KMeans clustering (based on PCA embeddings)
for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  kmeans_result <- kmeans(Embeddings(merged_obj, "pca"), centers = k)$cluster
  merged_obj[[cname]] <- as.factor(kmeans_result)
}


# Assign seurat_clusters (if not already present)
merged_obj$seurat_clusters <- Idents(merged_obj)


# Build the misc$clusters list
cluster_list <- list(
  active_cluster  = "seurat_clusters",
  seurat_clusters = merged_obj$seurat_clusters,
  orig.ident      = merged_obj$orig.ident,
  graph_based     = merged_obj$seurat_clusters
)


# Add KMeans clusters to misc$clusters
for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  vec <- merged_obj@meta.data[[cname]]
  vec <- factor(vec)
  names(vec) <- colnames(merged_obj)
  cluster_list[[cname]] <- vec
}


# Assign the full cluster list to misc
merged_obj@misc$clusters <- cluster_list


# Export to Loupe .cloupe file
create_loupe_from_seurat(
  obj = merged_obj,
  output_dir = "C:/scRNAseq/Shiny_App_Project/outputs",
  output_name = "merged_obj_KMeans_cloupe",
  dedup_clusters = FALSE,
  force = TRUE
)


