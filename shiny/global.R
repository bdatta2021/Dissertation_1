
 
# library(Seurat)
# library(loupeR)


# --- global.R (runs before ui/server) ---

# Allow very large uploads
options(shiny.maxRequestSize = 5 * 1024^3)   # 5 GB


# Make sure Bioconductor repos are available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
options(repos = BiocManager::repositories())



# Auto-install any missing packages (CRAN + Bioc)


need <- function(pkgs, bioc = FALSE) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) {
    if (bioc) BiocManager::install(missing, ask = FALSE, update = FALSE)
    else      install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}



# ---- LIST the REAL DEPENDENCIES HERE ----

cran_pkgs <- c(
  "shiny","ggplot2","dplyr","tidyr","Matrix","DT","plotly",
  "stringr","readr","magrittr","shinythemes","shinycssloaders"
)
bioc_pkgs <- c(
  "SeuratObject","SingleCellExperiment","BiocParallel"  # add/remove to match your code
)


need(cran_pkgs, bioc = FALSE)
need(bioc_pkgs,  bioc = TRUE)


