


library(Seurat)

# 1) point to your parent directory of those 10 folders
vis_dir   <- "C:/scRNASeq/Shiny_App_Project/Visium"
out_dir   <- file.path(vis_dir, "processed_rds")
dir.create(out_dir, showWarnings = FALSE)

# 2) get each sample folder (no recursion)
sample_dirs <- list.dirs(vis_dir, recursive = FALSE, full.names = TRUE)

for (d in sample_dirs) {
  message("Loading ", basename(d), " …")
  obj <- Load10X_Spatial(data.dir = d)
  save_file <- file.path(out_dir, paste0(basename(d), ".rds"))
  message("Saving to ", save_file)
  saveRDS(obj, save_file)
}

