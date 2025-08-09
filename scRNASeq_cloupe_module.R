


# scRNASeq_cloupe_module.R

library(shiny)
library(loupeR)

# — UI Module ——————————————————————————————————————————————————————
scRNASeqCloupeUI <- function(id) {
  ns <- NS(id)
  tagList(
    # **no fileInput here**: we use the one in app.R
    textInput(
      ns("output_path"), 
      "Output Path (.cloupe):",
      placeholder = "C:/scRNAseq/Shiny_App_Project/output.cloupe",
      value       = "output.cloupe"
    ),
    actionButton(ns("generate_cloupe"), "Generate .cloupe File"),
    verbatimTextOutput(ns("cloupe_status"))
  )
}

# — Server Module ———————————————————————————————————————————————————
# now takes a `files` argument
scRNASeqCloupeServer <- function(id, files) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$generate_cloupe, {
      
      tbl <- files()
      if (is.null(tbl) || nrow(tbl) == 0) {
        output$cloupe_status <- renderText("Please upload at least one .rds in the Export tab.")
        return()
      }
      # take the first uploaded .rds
      rds_path <- tbl$datapath[[1]]
      
      # decide output filename
      out   <- input$output_path
      if (!grepl("\\.cloupe$", out, ignore.case=TRUE)) {
        out <- paste0(out, ".cloupe")
      }
      
      # load
      seurat_obj <- tryCatch(readRDS(rds_path), error=function(e){
        output$cloupe_status <- renderText(paste("Failed to read .rds:", e$message))
        NULL
      })
      if (is.null(seurat_obj)) return()
      
      # build cluster list
      for (k in 2:10) {
        cname <- paste0("kmeans_", k)
        km    <- kmeans(Embeddings(seurat_obj, "pca"), centers = k)$cluster
        seurat_obj[[cname]] <- factor(km)
      }
      cl_list <- list(
        active_cluster  = "seurat_clusters",
        seurat_clusters = seurat_obj$seurat_clusters,
        orig.ident      = seurat_obj$orig.ident
      )
      for (k in 2:10) {
        cname <- paste0("kmeans_", k)
        vec   <- seurat_obj@meta.data[[cname]]
        names(vec) <- colnames(seurat_obj)
        cl_list[[cname]] <- vec
      }
      seurat_obj@misc$clusters <- cl_list
      
      # export
      tryCatch({
        create_loupe_from_seurat(
          obj          = seurat_obj,
          output_name  = tools::file_path_sans_ext(basename(out)),
          output_dir   = dirname(out),
          dedup_clusters = FALSE,
          force          = TRUE
        )
        output$cloupe_status <- renderText(paste("Written .cloupe to", out))
      }, error=function(e){
        output$cloupe_status <- renderText(paste("Export failed:", e$message))
      })
    })
  })
}



