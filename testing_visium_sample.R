



# Replace with wherever your ZIP lives:
zip_path <- "C:/Testing file/11425933A3.zip"

# 1) unpack to a temp directory
tmp <- file.path(tempdir(), "visium_quick")
unlink(tmp, recursive=TRUE); dir.create(tmp, recursive=TRUE)
unzip(zip_path, exdir=tmp)

# 2) find the actual data folder (descend if there’s a single top‐level dir)
roots    <- list.dirs(tmp, full.names=TRUE, recursive=FALSE)
data_dir <- if(length(roots)==1) roots else tmp

# 3) load with Seurat and count genes & spots
library(Seurat); library(SeuratData)
obj <- Load10X_Spatial(data.dir = data_dir)
cts <- GetAssayData(obj, assay="Spatial", slot="counts")
message("Genes: ", nrow(cts),  " —  Cells: ", ncol(cts))

