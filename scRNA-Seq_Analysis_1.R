


install.packages("Seurat")

library(Seurat)

# SCT1

sct1_path  <- "C:/scRNAseq/Shiny_App_Project/scRNASeq/SCT1"

sct1_data <- Read10X(data.dir = sct1_path)

sct1 <- CreateSeuratObject(counts = sct1_data, project = "SCT1")

sct1

head(sct1@meta.data)



sct1_counts <- Read10X(data.dir = "C:/scRNAseq/Shiny_App_Project/scRNASeq/SCT1")





sct1_assay <- CreateAssayObject(counts = sct1_counts)






# Step 3: Create Seurat objects with assays

sct1_obj <- CreateSeuratObject(counts = sct1_assay, assay = "RNA", project = "SCT1")





# Add mitochondrial percentage

sct1_obj[["percent.mt"]] <- PercentageFeatureSet(sct1_obj, pattern = "^MT-")


# Visualize the QC metrices before fltering


VlnPlot(sct1_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Step_5: Scattar plots for more insight

FeatureScatter(sct1_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(sct1_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# Step_6: Apply standard filtering (no saving yet )

# adjusting thresholds as needed

sct1_qc <- subset(sct1_obj, subset = nFeature_RNA> 200 & nFeature_RNA < 6000)



#  if i want to filter percent.mt too:

# sct1_qc <- subset(sct1_qc, subset = percent.mt < 10)



# Step_7: check how many cells were kept

dim(sct1_obj)

dim(sct1_qc)


#  Viewing metadata for QC cells

head(sct1_qc@meta.data)


# 2. Normalize & PCA for SCT1

sct1_qc <- NormalizeData(sct1_qc)
sct1_qc <- FindVariableFeatures(sct1_qc)
sct1_qc <- ScaleData(sct1_qc)
sct1_qc <- RunPCA(sct1_qc, npcs = 30)


# 3. Run Doublet Finder

# pick pN, pK,nExp based on guidence


library(DoubletFinder)

nExp <- round(0.075 * ncol(sct1_qc))


sct1_qc <- doubletFinder(sct1_qc, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp)

table(sct1_qc@meta.data$DF.classifications_0.25_0.09_843)

colnames(sct1_qc@meta.data)



# 4. Remove doublets
sct1_qc <- subset(sct1_qc, subset = DF.classifications_0.25_0.09_843 == 'Singlet')



# 5. Final normalization and clustering 

sct1_qc <- ScaleData(sct1_qc)
sct1_qc <- RunPCA(sct1_qc, npcs = 30)

sct1_qc  <- RunTSNE(sct1_qc, dims = 1:30)
sct1_qc  <- FindNeighbors(sct1_qc, dims = 1:30)
sct1_qc  <- FindClusters(sct1_qc, resolution = 0.6)  # or whichever value you choose
sct1_qc  <- RunUMAP(sct1_qc, dims = 1:30)



DimPlot(sct1_qc, reduction = "umap", label = TRUE)


# identify marker genes for each cluster


# markers_sct1 <- FindAllMarkers(sct1_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



# DimPlot(sct2_obj, reduction = "umap", label = TRUE)










# SCT2

sct2_path  <- "C:/scRNAseq/Shiny_App_Project/scRNASeq/SCT2"

sct2_data <- Read10X(data.dir = sct2_path)

sct2 <- CreateSeuratObject(counts = sct2_data, project = "SCT2")

sct2

head(sct2@meta.data)




sct2_counts <- Read10X(data.dir = "C:/scRNAseq/Shiny_App_Project/scRNASeq/SCT2")


sct2_assay <- CreateAssayObject(counts = sct2_counts)


sct2_obj <- CreateSeuratObject(counts = sct2_assay, assay = "RNA", project = "SCT2")


sct2_obj[["percent.mt"]] <- PercentageFeatureSet(sct2_obj, pattern = "^MT-")



# Visualize the QC metrices before fltering


VlnPlot(sct2_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Step_5: Scattar plots for more insight

FeatureScatter(sct2_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(sct2_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")




# Step_6: Apply standard filtering (no saving yet )

# adjusting thresholds as needed

sct2_qc <- subset(sct2_obj, subset = nFeature_RNA> 200 & nFeature_RNA < 6000)






# if i want to filter percent.mt too:

# sct2_qc <- subset(sct2_qc, subset = percent.mt < 10)



# Step_7: check how many cells were kept

dim(sct2_obj)

dim(sct2_qc)


#  Viewing metadata for QC cells

head(sct2_qc@meta.data)




# 2. Normalize & PCA for SCT2

sct2_qc <- NormalizeData(sct2_qc)
sct2_qc <- FindVariableFeatures(sct2_qc)
sct2_qc <- ScaleData(sct2_qc)
sct2_qc <- RunPCA(sct2_qc, npcs = 30)


# 3. Run Doublet Finder

# pick pN, pK,nExp based on guidence

library(DoubletFinder)

nExp <- round(0.075 * ncol(sct2_qc))


sct2_qc <- doubletFinder(sct2_qc, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp)

colnames(sct2_qc@meta.data)


table(sct2_qc@meta.data$DF.classifications_0.25_0.09_931)




# 4. Remove doublets
sct2_qc <- subset(sct2_qc, subset = DF.classifications_0.25_0.09_931 == 'Singlet')



# 5. Final normalization and clustering 

sct2_qc <- ScaleData(sct2_qc)
sct2_qc <- RunPCA(sct2_qc, npcs = 30)

sct2_qc <- RunTSNE(sct2_qc, dims = 1:30)
sct2_qc <- FindNeighbors(sct2_qc, dims = 1:30)
sct2_qc <- FindClusters(sct2_qc, resolution = 0.57)  # or whichever value you choose
sct2_qc <- RunUMAP(sct2_qc, dims = 1:30)



DimPlot(sct2_qc, reduction = "umap", label = TRUE)


# identify marker genes for each cluster


# markers_sct2 <- FindAllMarkers(sct2_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)






merged_obj <- merge(sct1_qc, y = sct2_qc, add.cell.ids = c("SCT1", "SCT2"), project = "merged_scRNASeq")




# Split the merged object by sample
split_obj <- SplitObject(merged_obj, split.by = "orig.ident")

# Normalize each using SCTransform (recommended for batch correction)
split_obj <- lapply(split_obj, SCTransform, verbose = FALSE)

# Integration steps
features <- SelectIntegrationFeatures(object.list = split_obj, nfeatures = 3000)

split_obj <- PrepSCTIntegration(object.list = split_obj, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = split_obj, normalization.method = "SCT", anchor.features = features)

merged_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")




DefaultAssay(merged_obj) <- "integrated"
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 30)
merged_obj <- RunUMAP(merged_obj, dims = 1:30)
merged_obj <- RunTSNE(merged_obj, dims = 1:30)
merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.6)  # or adjust as needed


cluster_markers <- FindAllMarkers(merged_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)





#Cell Type Annotation

# Step 1: If not already installed
install.packages("BiocManager")

# Step 2: Use BiocManager to install SingleR and celldex
BiocManager::install("SingleR")
BiocManager::install("celldex")

library(SingleR)
library(celldex)

ref <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = GetAssayData(merged_obj, assay = "integrated", slot = "data"),
                ref = ref,
                labels = ref$label.main)
merged_obj$SingleR_label <- pred$labels

DimPlot(merged_obj, group.by = "SingleR_label", label = TRUE)




#FeaturePlot / Heatmap of top markers (Visualize cluster-defining genes)

library(dplyr)


top10 <- cluster_markers %>% group_by(cluster) %>% top_n(1, avg_log2FC)
FeaturePlot(merged_obj, features = top10$gene)
DoHeatmap(merged_obj, features = top10$gene) + NoLegend()



# Save the object

#saveRDS(sct1_qc, file = "C:/scRNAseq/Shiny_App_Project/outputs/sct1_obj_before_merge.rds")

#saveRDS(sct2_qc, file = "C:/scRNAseq/Shiny_App_Project/outputs/sct2_obj_before_merge.rds")


saveRDS(merged_obj, "outputs/merged_scRNA-Seq_Final.rds")




library(loupeR)



library(Seurat)


# Now this will work
cluster_ids <- Idents(merged_obj)

#Convert cluster identities to a factor with cell names
cluster_ids <- Idents(merged_obj)
cluster_vector <- setNames(as.character(cluster_ids), colnames(merged_obj))


# Register clustering identities
merged_obj@misc$clusters <- list("RNA_snn_res.0.6" = setNames(as.character(Idents(merged_obj)), colnames(merged_obj)))

export_loupe(obj = merged_obj, file = "merged_obj_output.cloupe")

# Add RNA assay using SCT counts
merged_obj[["RNA"]] <- CreateAssayObject(
  counts = GetAssayData(merged_obj, assay = "RNA", slot = "counts")
)







# Step 1: Run KMeans Clustering on PCA embeddings
for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  sct1_qc[[cname]] <- kmeans(Embeddings(sct1_qc, "pca"), centers = k)$cluster
}

# Step 2: Construct misc$clusters
cluster_list <- list(
  active_cluster  = "seurat_clusters",
  seurat_clusters = sct1_qc$seurat_clusters,
  orig.ident      = sct1_qc$orig.ident,
  graph_based     =sct1_qc$seurat_clusters
)

# Add K Means 2-10 to cluster_list

for (k in 2:10) {
  cname <- paste0("kmeans_", k)
  vec <- sct1_qc[[cname]] [, 1]
  vec <- factor(vec)
  names(vec) <- colnames(sct1_qc)
  cluster_list[[cname]] <- vec
}


# Assign to misc

sct1_qc@misc$clusters <- cluster_list

# Step 3: Export to .cloupe
create_loupe_from_seurat(
  obj = sct1_qc,
  output_dir = "C:/scRNAseq/Shiny_App_Project/outputs",
  output_name = "SCT1_KMeans_cloupe",
  dedup_clusters = FALSE,
  force = TRUE
)







# Step : Export .cloupe files for Loupe Browser

#library(LoupeR)

#export_loupe(obj = sct1_qc, file = "output_SCT1.cloupe")

#export_loupe(obj = sct2_qc, file = "output_SCT2.cloupe")

#export_loupe(obj = merged_obj, file = "output_merged_obj.cloupe")







