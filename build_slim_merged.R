


# scripts/build_slim_merged.R
#
# 1) Load the two per-sample slim objects you already built
sct1 <- readRDS("data/SCT1_slim.rds")
sct2 <- readRDS("data/SCT2_slim.rds")

# ────────────────────────────────────────────────────────────────────────────
# 2) Merge the two samples
merged <- merge(
  x = sct1,
  y = sct2,
  add.cell.ids = c("SCT1", "SCT2"),
  project = "merged_scRNASeq"
)

# 3) Split and run SCTransform on each (batch correction prep)
split_list <- SplitObject(merged, split.by = "orig.ident")
split_list <- lapply(split_list, SCTransform, verbose = FALSE)

# 4) Find integration features & anchors
features <- SelectIntegrationFeatures(object.list = split_list, nfeatures = 3000)
split_list <- PrepSCTIntegration(object.list = split_list,
                                 anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = split_list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

# 5) Integrate
merged_int <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# 6) Scale, PCA, UMAP, t-SNE, neighbors, clustering
DefaultAssay(merged_int) <- "integrated"
merged_int <- ScaleData(merged_int, verbose = FALSE)
merged_int <- RunPCA(merged_int, npcs = 30, verbose = FALSE)
merged_int <- RunUMAP(merged_int, dims = 1:30, verbose = FALSE)
merged_int <- RunTSNE(merged_int, dims = 1:30, verbose = FALSE)
merged_int <- FindNeighbors(merged_int, dims = 1:30, verbose = FALSE)
merged_int <- FindClusters(merged_int, resolution = 0.6, verbose = FALSE)

# 7) (Optional) Cell-type annotation / marker finding here if you like—
#    but remember: heavyweight steps belong upstream of Shiny deployment.

# ────────────────────────────────────────────────────────────────────────────
# 8) Diet down to only what the app displays
#    keep: integrated assay, chosen embeddings, and meta data
slim_merged <- DietSeurat(
  object    = merged_int,
  assays    = "integrated",
  features  = rownames(merged_int),           # or a shorter vector of marker genes
  dimreducs = names(merged_int@reductions),   # e.g. c("pca","umap","tsne")
  graphs    = NULL,                           # drop any neighbor‐graph slots
  misc      = merged_int@misc                 # cluster lists, SingleR labels, etc.
)


# 9) Save the tiny merged object (<50 MB)
dir.create("data", showWarnings = FALSE, recursive = TRUE)
saveRDS(slim_merged, file = "data/merged_slim.rds")




