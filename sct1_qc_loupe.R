

# Part_1

library(Seurat)

sct1_path <- "C:/scRNAseq/Shiny_App_Project/scRNASeq/SCT1"
sct1_data <- Read10X(data.dir = sct1_path)
sct1 <- CreateSeuratObject(counts = sct1_data, project = "SCT1")


# Part_2

sct1[["percent.mt"]] <- PercentageFeatureSet(sct1, pattern = "^MT-")

VlnPlot(sct1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(sct1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(sct1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sct1_qc <- subset(sct1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)


# part_3


sct1_qc <- NormalizeData(sct1_qc)
sct1_qc <- FindVariableFeatures(sct1_qc)
sct1_qc <- ScaleData(sct1_qc)
sct1_qc <- RunPCA(sct1_qc, npcs = 30)

library(DoubletFinder)
nExp <- round(0.075 * ncol(sct1_qc))
sct1_qc <- doubletFinder(sct1_qc, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp)

sct1_qc <- subset(sct1_qc, subset = DF.classifications_0.25_0.09_843 == 'Singlet')



# part_4




sct1_qc <- ScaleData(sct1_qc)
sct1_qc <- RunPCA(sct1_qc, npcs = 30)
sct1_qc <- RunTSNE(sct1_qc, dims = 1:30)
sct1_qc <- FindNeighbors(sct1_qc, dims = 1:30)
sct1_qc <- FindClusters(sct1_qc, resolution = 0.6)
sct1_qc  <- RunUMAP(sct1_qc, dims = 1:30)

sct1_qc$seurat_clusters <- Idents(sct1_qc)  # ✅ This is critical!



# part_5


for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  kmeans_result <- kmeans(Embeddings(sct1_qc, "pca"), centers = k)$cluster
  sct1_qc[[cname]] <- as.factor(kmeans_result)
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
  vec <- sct1_qc@meta.data[[cname]]  # ✅ Fix: extract as vector, not data.frame
  vec <- factor(vec)
  names(vec) <- colnames(sct1_qc)
  cluster_list[[cname]] <- vec
}

# ✅ Final assignment
sct1_qc@misc$clusters <- cluster_list




# part_7


library(loupeR)

create_loupe_from_seurat(
  obj = sct1_qc,
  output_dir = "C:/scRNAseq/Shiny_App_Project/outputs",
  output_name = "SCT1_KMeans_17_cloupe",
  dedup_clusters = FALSE,
  force = TRUE
)



