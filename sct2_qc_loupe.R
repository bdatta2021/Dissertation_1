


# For exporting clope file of SCT2


# Part_1

library(Seurat)



sct2_path  <- "C:/scRNAseq/Shiny_App_Project/scRNASeq/SCT2"

sct2_data <- Read10X(data.dir = sct2_path)

sct2 <- CreateSeuratObject(counts = sct2_data, project = "SCT2")






# Part_2

sct2[["percent.mt"]] <- PercentageFeatureSet(sct2, pattern = "^MT-")

VlnPlot(sct2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(sct2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(sct2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sct2_qc <- subset(sct2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)


# part_3


sct2_qc <- NormalizeData(sct2_qc)
sct2_qc <- FindVariableFeatures(sct2_qc)
sct2_qc <- ScaleData(sct2_qc)
sct2_qc <- RunPCA(sct2_qc, npcs = 30)

library(DoubletFinder)

nExp <- round(0.075 * ncol(sct2_qc))
sct2_qc <- doubletFinder(sct2_qc, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp)

sct2_qc <- subset(sct2_qc, subset = DF.classifications_0.25_0.09_931 == 'Singlet')



# part_4




sct2_qc <- ScaleData(sct2_qc)
sct2_qc <- RunPCA(sct2_qc, npcs = 30)
sct2_qc <- RunTSNE(sct2_qc, dims = 1:30)
sct2_qc <- FindNeighbors(sct2_qc, dims = 1:30)
sct2_qc <- FindClusters(sct2_qc, resolution = 0.57)
sct2_qc  <- RunUMAP(sct2_qc, dims = 1:30)

sct2_qc$seurat_clusters <- Idents(sct2_qc)  # ✅ This is critical!



# part_5


for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  kmeans_result <- kmeans(Embeddings(sct2_qc, "pca"), centers = k)$cluster
  sct2_qc[[cname]] <- as.factor(kmeans_result)
}



# part_6

#cluster_list <- list()


cluster_list <- list(
  active_cluster  = "seurat_clusters",
  seurat_clusters = sct1_qc$seurat_clusters,
  orig.ident      = sct1_qc$orig.ident,
  graph_based     = sct1_qc$seurat_clusters
)

# Add KMeans clusters
for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  vec <- sct2_qc@meta.data[[cname]]  # ✅ Fix: extract as vector, not data.frame
  vec <- factor(vec)
  names(vec) <- colnames(sct2_qc)
  cluster_list[[cname]] <- vec
}

# ✅ Final assignment
sct2_qc@misc$clusters <- cluster_list




# part_7


library(loupeR)

create_loupe_from_seurat(
  obj = sct2_qc,
  output_dir = "C:/scRNAseq/Shiny_App_Project/outputs",
  output_name = "SCT2_KMeans_18_cloupe",
  dedup_clusters = FALSE,
  force = TRUE
)



