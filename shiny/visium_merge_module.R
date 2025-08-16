





# Final_4



library(shiny)
library(Seurat)
library(DT)

# ── UI Module ──────────────────────────────────────────────────────────────────
visiumMergeUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      # ── Sidebar (inputs + buttons) ────────────────────────────────────────────
      column(
        width = 4,
        fileInput(ns("sample1"),  "Upload Sample 1 (.rds):",   accept = ".rds"),
        fileInput(ns("sample2"),  "Upload Sample 2 (.rds):",   accept = ".rds"),
        fileInput(ns("sample3"),  "Upload Sample 3 (.rds):",   accept = ".rds"),
        fileInput(ns("sample4"),  "Upload Sample 4 (.rds):",   accept = ".rds"),
        fileInput(ns("sample5"),  "Upload Sample 5 (.rds):",   accept = ".rds"),
        fileInput(ns("sample6"),  "Upload Sample 6 (.rds):",   accept = ".rds"),
        fileInput(ns("sample7"),  "Upload Sample 7 (.rds):",   accept = ".rds"),
        fileInput(ns("sample8"),  "Upload Sample 8 (.rds):",   accept = ".rds"),
        fileInput(ns("sample9"),  "Upload Sample 9 (.rds):",   accept = ".rds"),
       # fileInput(ns("sample10"), "Upload Sample 10 (.rds):",  accept = ".rds"),
        actionButton(ns("merge_btn"),  "Merge & Normalize All Visium Samples"),
        downloadButton(ns("download_rds"), "Download Merged .rds")
      ),
      
      # ── Main (plot + table) ───────────────────────────────────────────────────
      column(
        width = 8,
        plotOutput(ns("spatial_plot")),
        DTOutput(ns("summary_table"))
      )
    )
  )
}

# ── Server Module ─────────────────────────────────────────────────────────────
visiumMergeServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # holds the merged object
    merged_visium <- reactiveVal(NULL)
    
    observeEvent(input$merge_btn, {
      # require all ten uploads
      req(input$sample1, input$sample2, input$sample3,
          input$sample4, input$sample5, input$sample6,
          input$sample7, input$sample8, input$sample9)
          
      # read each .rds into a list
      paths <- c(
        input$sample1$datapath, input$sample2$datapath,
        input$sample3$datapath, input$sample4$datapath,
        input$sample5$datapath, input$sample6$datapath,
        input$sample7$datapath, input$sample8$datapath,
        input$sample9$datapath 
      )
      sample_list <- lapply(paths, readRDS)
      names(sample_list) <- paste0("Sample", seq_along(sample_list))
      
      # merge them all
      merged <- Reduce(function(x, y) merge(x, y = y), sample_list)
      
      # ensure valid spatial spots (if present)
      if ("nCount_Spatial" %in% colnames(merged@meta.data)) {
        keep <- !is.na(merged$nCount_Spatial) & is.finite(merged$nCount_Spatial)
        merged <- subset(merged, cells = colnames(merged)[keep])
      }
      
      # SCTransform (fallback to Normalize+Scale)
      tryCatch({
        merged <- SCTransform(merged, assay = "Spatial", verbose = FALSE)
      }, error = function(e) {
        showNotification(
          "SCTransform failed — using NormalizeData + ScaleData",
          type     = "warning",
          duration = 5
        )
        DefaultAssay(merged) <- "Spatial"
        merged <- NormalizeData(merged)
        merged <- ScaleData(merged)
      })
      
      # dimensionality reduction
      merged <- RunPCA(merged)
      merged <- RunUMAP(merged, dims = 1:10)
      merged <- RunTSNE(merged, dims = 1:10)
      
      # store for outputs
      merged_visium(merged)
    })
    
    # ── Plot ─────────────────────────────────────────────────────────────────
    output$spatial_plot <- renderPlot({
      req(merged_visium())
      obj <- merged_visium()
      if ("spatial" %in% names(obj@images)) {
        SpatialDimPlot(obj, label = TRUE)
      } else {
        DimPlot(obj, reduction = "umap", label = TRUE)
      }
    })
    
    # ── Summary table ─────────────────────────────────────────────────────────
    output$summary_table <- renderDT({
      req(merged_visium())
      obj <- merged_visium()
      
      # grab counts from the first valid assay
      mat <- NULL
      for (a in Assays(obj)) {
        mat <- tryCatch(GetAssayData(obj[[a]], slot = "counts"),
                        error = function(e) NULL)
        if (!is.null(mat) && nrow(mat) && ncol(mat)) break
      }
      req(mat)
      
      df <- data.frame(
        Cells           = ncol(mat),
        Features        = nrow(mat),
        SparsityPercent = round(100 * sum(mat == 0) / length(mat), 2),
        NonZeroEntries  = sum(mat != 0)
      )
      datatable(df, options = list(dom = 't', paging = FALSE))
    })
    
    # ── Download Handler ─────────────────────────────────────────────────────
    output$download_rds <- downloadHandler(
      filename = function() {
        paste0("merged_visium_", Sys.Date(), ".rds")
      },
      content = function(file) {
        req(merged_visium())
        saveRDS(merged_visium(), file)
      }
    )
  })
}













