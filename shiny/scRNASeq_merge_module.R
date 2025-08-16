









# scRNASeq_merge_module.R



library(shiny)
library(Seurat)
library(SeuratObject)
library(DT)

# ── UI ─────────────────────────────────────────────────────────────────────────
scRNASeqMergeUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 4,
        # two separate uploads
        fileInput(ns("sample1"), "Upload Sample 1 (.rds):", accept = ".rds"),
        fileInput(ns("sample2"), "Upload Sample 2 (.rds):", accept = ".rds"),
        actionButton(ns("merge_btn"), "Merge & Normalize"),
        downloadButton(ns("download_rds"), "Download Merged .rds")
      ),
      column(
        width = 8,
        plotOutput(ns("umap_plot")),
        DTOutput(ns("summary_table"))
      )
    )
  )
}

# ── SERVER ─────────────────────────────────────────────────────────────────────
scRNASeqMergeServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    merged_obj <- reactiveVal(NULL)

    observeEvent(input$merge_btn, {
      req(input$sample1, input$sample2)
      obj1 <- readRDS(input$sample1$datapath)
      obj2 <- readRDS(input$sample2$datapath)

      DefaultAssay(obj1) <- "RNA"
      DefaultAssay(obj2) <- "RNA"

      m <- merge(obj1, y = obj2, add.cell.ids = c("S1","S2"))
      DefaultAssay(m) <- "RNA"
      m <- NormalizeData(m)
      m <- FindVariableFeatures(m)
      m <- ScaleData(m)
      m <- RunPCA(m)
      m <- RunUMAP(m, dims = 1:10)
      m <- RunTSNE(m, dims = 1:10)

      merged_obj(m)
    })

    output$umap_plot <- renderPlot({
      req(merged_obj())
      DimPlot(merged_obj(), reduction = "umap", label = TRUE)
    })

    output$summary_table <- renderDT({
      req(merged_obj())
      obj <- merged_obj()
      mat <- GetAssayData(obj, assay="RNA", slot="counts")
      df <- data.frame(
        Cells    = ncol(mat),
        Features = nrow(mat),
        Sparsity = round(100 * sum(mat==0)/(length(mat)),2)
      )
      datatable(df, options=list(dom="t", paging=FALSE))
    })

    output$download_rds <- downloadHandler(
      filename = function() {
        paste0("merged_scRNASeq_", Sys.Date(), ".rds")
      },
      content = function(file) {
        req(merged_obj())
        saveRDS(merged_obj(), file)
      }
    )
  })
}




# Edition_1

# library(shiny)
# library(Seurat)
# library(SeuratObject)
# library(DT)
#
# # ── UI ─────────────────────────────────────────────────────────────────────────
# scRNASeqMergeUI <- function(id) {
#   ns <- NS(id)
#   tagList(
#     fluidRow(
#
#       # ── Sidebar column ─────────────────────────────────────────────────────────
#       column(
#         width = 4,
#         # two separate uploads
#         fileInput(  ns("sample1"), "Upload Sample 1 (.rds):", accept = ".rds"),
#         fileInput(  ns("sample2"), "Upload Sample 2 (.rds):", accept = ".rds"),
#         actionButton(ns("merge_btn"),  "Merge & Normalize"),
#
#         ## ─── YOUR THREE NEW SPINNERS ────────────────────────────────
#         numericInput(
#           inputId = ns("minCellsPerGene"),
#           label   = "Minimum Cells per Gene:",
#           value   = 3,
#           min     = 0,
#           step    = 1
#         ),  ## ← NEW
#
#         numericInput(
#           inputId = ns("minFeaturesPerCell"),
#           label   = "Minimum Features per Cell:",
#           value   = 200,
#           min     = 0,
#           step    = 1
#         ),  ## ← NEW
#
#         numericInput(
#           inputId = ns("nRows"),
#           label   = "Rows to show in table:",
#           value   = 5,
#           min     = 1,
#           max     = 100,
#           step    = 1
#         ),  ## ← NEW
#         ## ─────────────────────────────────────────────────────────────
#
#         downloadButton(ns("download_rds"), "Download Merged .rds")
#       ),
#
#       # ── Main column ────────────────────────────────────────────────────────────
#       column(
#         width = 8,
#         plotOutput(ns("umap_plot")),
#         DTOutput(ns("summary_table"))
#       )
#     )
#   )
# }
#
# # ── SERVER ─────────────────────────────────────────────────────────────────────
# scRNASeqMergeServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#
#     merged_obj <- reactiveVal(NULL)
#
#     observeEvent(input$merge_btn, {
#       req(input$sample1, input$sample2)
#
#       # 1) read the two Seurat objects
#       obj1 <- readRDS(input$sample1$datapath)
#       obj2 <- readRDS(input$sample2$datapath)
#
#       DefaultAssay(obj1) <- "RNA"
#       DefaultAssay(obj2) <- "RNA"
#
#       # 2) merge them
#       m <- merge(obj1, y = obj2, add.cell.ids = c("S1","S2"))
#       DefaultAssay(m) <- "RNA"
#
#       ## ─── APPLY QC FILTERS FROM YOUR NEW SPINNERS ────────────────────────────
#       # a) drop genes seen in fewer than minCellsPerGene cells
#       counts <- GetAssayData(m, assay = "RNA", slot = "counts")
#       keep_genes <- Matrix::rowSums(counts > 0) >= input$minCellsPerGene
#       m <- subset(m, features = rownames(counts)[keep_genes])
#
#       # b) drop cells with fewer than minFeaturesPerCell features
#       m <- subset(m, subset = nFeature_RNA >= input$minFeaturesPerCell)
#       ## ────────────────────────────────────────────────────────────────────────
#
#       # 3) standard preprocessing
#       m <- NormalizeData(m)
#       m <- FindVariableFeatures(m)
#       m <- ScaleData(m)
#       m <- RunPCA(m)
#       m <- RunUMAP(m, dims = 1:10)
#       m <- RunTSNE(m, dims = 1:10)
#
#       # store for plotting & table
#       merged_obj(m)
#     })
#
#     # ── UMAP PLOT ────────────────────────────────────────────────────────────────
#     output$umap_plot <- renderPlot({
#       req(merged_obj())
#       DimPlot(merged_obj(), reduction = "umap", label = TRUE)
#     })
#
#     # ── SUMMARY TABLE ────────────────────────────────────────────────────────────
 #   output$summary_table <- renderDT({
#       req(merged_obj())
#       obj <- merged_obj()
#       mat <- GetAssayData(obj, assay = "RNA", slot = "counts")
#       df <- data.frame(
#         Cells    = ncol(mat),
#         Features = nrow(mat),
#         Sparsity = round(100 * sum(mat == 0) / length(mat), 2)
#       )
#       datatable(
#         df,
#         options = list(
#           dom          = "t",
#           pageLength   = input$nRows,   ## ← NEW: pagesize follows spinner
#           lengthChange = FALSE
#         )
#       )
#     })
#
#     # ── DOWNLOAD HANDLER ────────────────────────────────────────────────────────
#     output$download_rds <- downloadHandler(
#       filename = function() {
#         paste0("merged_scRNASeq_", Sys.Date(), ".rds")
#       },
#       content = function(file) {
#         req(merged_obj())
#         saveRDS(merged_obj(), file)
#       }
#     )
#
#   })
# }
#
