


install.packages("Seurat")

# scripts/build_slim_sct2.R
library(Seurat)
library(DoubletFinder)

raw_path <- "C:/outputs/SCT2"
sct2     <- CreateSeuratObject(Read10X(data.dir = raw_path),
                               project = "SCT2")

sct2[["percent.mt"]] <- PercentageFeatureSet(sct2, pattern = "^MT-")

VlnPlot(sct2, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
FeatureScatter(sct2, "nCount_RNA", "percent.mt")
FeatureScatter(sct2, "nCount_RNA", "nFeature_RNA")

sct2 <- subset(sct2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

sct2 <- NormalizeData(sct2)
sct2 <- FindVariableFeatures(sct2)
sct2 <- ScaleData(sct2)
sct2 <- RunPCA(sct2, npcs = 30)

nExp <- round(0.075 * ncol(sct2))
sct2 <- doubletFinder(sct2, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp)
sct2 <- subset(sct2, subset = DF.classifications_0.25_0.09_931 == "Singlet")

sct2 <- ScaleData(sct2)
sct2 <- RunPCA(sct2, npcs = 30)
sct2 <- RunTSNE(sct2, dims = 1:30)
sct2 <- FindNeighbors(sct2, dims = 1:30)
sct2 <- FindClusters(sct2, resolution = 0.57)
sct2 <- RunUMAP(sct2, dims = 1:30)

DimPlot(sct2, reduction = "umap", label = TRUE)

# Diet
DefaultAssay(sct2) <- "RNA"
slim_sct2 <- DietSeurat(
  object    = sct2,
  assays    = "RNA",
  features  = rownames(sct2),
  dimreducs = names(sct2@reductions),
  graphs    = NULL,
  misc      = sct2@misc
)

dir.create("data", showWarnings = FALSE)


saveRDS(slim_sct1, file = "data/SCT2_slim.rds")

