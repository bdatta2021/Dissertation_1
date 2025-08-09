



# scripts/build_slim_sct1.R
#
# 1) Library setup & read raw 10X
library(Seurat)
library(DoubletFinder)

raw_path <- "C:/outputs/SCT1"
sct1     <- CreateSeuratObject(Read10X(data.dir = raw_path),
                               project = "SCT1")

# 2) Compute QC metrics
sct1[["percent.mt"]] <- PercentageFeatureSet(sct1, pattern = "^MT-")

# 3) Visual QC (optional when developing)
VlnPlot(sct1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(sct1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(sct1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 4) Subset / filter
sct1 <- subset(sct1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# 5) Normalize → find HVGs → scale → PCA
sct1 <- NormalizeData(sct1)
sct1 <- FindVariableFeatures(sct1)
sct1 <- ScaleData(sct1)
sct1 <- RunPCA(sct1, npcs = 30)

# 6) Doublet detection
nExp <- round(0.075 * ncol(sct1))
sct1 <- doubletFinder(sct1, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp)
sct1 <- subset(sct1, subset = DF.classifications_0.25_0.09_843 == "Singlet")

# 7) Clustering & embeddings
sct1 <- ScaleData(sct1)
sct1 <- RunPCA(sct1, npcs = 30)
sct1 <- RunTSNE(sct1, dims = 1:30)
sct1 <- FindNeighbors(sct1, dims = 1:30)
sct1 <- FindClusters(sct1, resolution = 0.6)
sct1 <- RunUMAP(sct1, dims = 1:30)

# 8) Optionally inspect final plot
DimPlot(sct1, reduction = "umap", label = TRUE)

# ────────────────────────────────────────────────────────────────────────────
# 9) Now *diet* down to only what the Shiny app will visualize:
DefaultAssay(sct1) <- "RNA"     # or "integrated" if you were using that assay
slim_sct1 <- DietSeurat(
  object    = sct1,
  assays    = "RNA",
  features  = rownames(sct1),              # or a shorter vector of marker genes
  dimreducs = names(sct1@reductions),      # e.g. c("pca","umap","tsne")
  graphs    = NULL,                        # drop neighbor graphs
  misc      = sct1@misc                    # keep your cluster lists & metadata
)

# 10) Save the *tiny* object (<50 MB) into your app folder
dir.create("data", showWarnings = FALSE)

saveRDS(slim_sct1, file = "data/SCT1_slim.rds")


