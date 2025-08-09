# 
# # visium_cloupe_module.R
# 
# visiumCloupeUI <- function(id) {
#   ns <- NS(id)
#   tagList(
#     fileInput(ns("rds_file"), "Upload Merged Seurat Object (.rds):"),
#     textInput(ns("output_path"), "Output Path (.cloupe):", value = "C:/scRNAseq/Shiny_App_Project/output_visium"),
#     actionButton(ns("generate_btn"), "Generate .cloupe File"),
#     verbatimTextOutput(ns("status"))
#   )
# }
# 
# visiumCloupeServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#     output$status <- renderText({ "" })
# 
#     observeEvent(input$generate_btn, {
#       req(input$rds_file)
#       req(input$output_path)
# 
#       tryCatch({
#         library(Seurat)
#         library(loupeR)
# 
#         # Load uploaded Seurat object
#         combined_visium <- readRDS(input$rds_file$datapath)
# 
#         # Run KMeans clustering (2 to 10)
#         for (k in 2:10) {
#           cname <- paste0("kmeans_", k)
#           kmeans_result <- kmeans(Embeddings(combined_visium, "pca"), centers = k)$cluster
#           combined_visium[[cname]] <- as.factor(kmeans_result)
#         }
# 
#         # Prepare LoupeR cluster groups
#         cluster_list <- list(
#           active_cluster  = "seurat_clusters",
#           seurat_clusters = combined_visium$seurat_clusters,
#           orig.ident      = combined_visium$orig.ident,
#           graph_based     = combined_visium$seurat_clusters
#         )
# 
#         for (k in 2:10) {
#           cname <- paste0("kmeans_", k)
#           vec <- combined_visium@meta.data[[cname]]
#           vec <- factor(vec)
#           names(vec) <- colnames(combined_visium)
#           cluster_list[[cname]] <- vec
#         }
# 
#         # Assign to misc slot
#         combined_visium@misc$clusters <- cluster_list
# 
#         # Export using LoupeR
#         create_loupe_from_seurat(
#           obj = combined_visium,
#           output_dir = dirname(input$output_path),
#           output_name = tools::file_path_sans_ext(basename(input$output_path)),
#           dedup_clusters = FALSE,
#           force = TRUE
#         )
# 
#       #  output$status <- renderText("Loupe .cloupe file generated with KMeans_2 to KMeans_10.")
# 
#        output$status <- renderText("Successfully exported .cloupe file to C:/scRNAseq/Shiny_App_Project/output_visium  ")
#       }, error = function(e) {
#         output$status <- renderText(paste("Failed to export .cloupe:", e$message))
#       })
#     })
#   })
# }







# its for the professor's wanted to see in this way and starting from july 28th at evening. The just above one used till july 28th till morning.







# visium_cloupe_module.R

library(shiny)
library(Seurat)
library(loupeR)

# — UI Module ——————————————————————————————————————————————————————
visiumCloupeUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(
      ns("rds_file"),
      "Upload Merged Visium Seurat Object (.rds):",
      accept = ".rds"
    ),
    textInput(
      ns("output_path"),
      "Output Path (no extension):",
      value = "output_visium"
    ),
    actionButton(ns("generate_btn"), "Generate .cloupe File"),
    verbatimTextOutput(ns("status"))
  )
}

# — Server Module ———————————————————————————————————————————————————
# Now takes a second argument `files` (reactive from the app)
visiumCloupeServer <- function(id, files) {
  moduleServer(id, function(input, output, session) {

    output$status <- renderText({ "" })

    observeEvent(input$generate_btn, {
      tbl <- files()
      req(tbl)
      # We only take the first uploaded .rds
      path_in  <- tbl$datapath[1]
      name_in  <- tbl$name[1]
      out_base <- input$output_path
      if (!nzchar(out_base)) {
        output$status("Please provide an output filename.")
        return()
      }

      tryCatch({
        # Load Seurat object
        vis_obj <- readRDS(path_in)

        # Compute K-means clusters 2–10
        for (k in 2:10) {
          cname <- paste0("kmeans_", k)
          km   <- kmeans(Embeddings(vis_obj, "pca"), centers = k)$cluster
          vis_obj[[cname]] <- factor(km)
        }

        # Build cluster list for LoupeR
        clust_list <- list(
          active_cluster  = "seurat_clusters",
          seurat_clusters = vis_obj$seurat_clusters,
          orig.ident      = vis_obj$orig.ident,
          graph_based     = vis_obj$seurat_clusters
        )
        for (k in 2:10) {
          cname <- paste0("kmeans_", k)
          vec   <- vis_obj@meta.data[[cname]]
          names(vec) <- colnames(vis_obj)
          clust_list[[cname]] <- vec
        }
        vis_obj@misc$clusters <- clust_list

        # Write out .cloupe
        create_loupe_from_seurat(
          obj         = vis_obj,
          output_dir  = dirname(out_base),
          output_name = tools::file_path_sans_ext(basename(out_base)),
          dedup_clusters = FALSE,
          force = TRUE
        )

        output$status <- renderText({
          paste0(".cloupe file written to ", out_base, ".cloupe")
        })
      }, error = function(e) {
        output$status <- renderText({
          paste("Failed to export .cloupe:", e$message)
        })
      })
    })
  })
}




