


# Testing_2




# visium_module.R

library(shiny)
library(Seurat)
library(SeuratData)   # for Load10X_Spatial()
library(DT)
library(patchwork)    # for wrap_plots()

# --------------------------
# UI Module
# --------------------------
visiumUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(  ns("zip"),               "Upload visium ZIP",      accept = ".zip"    ),
      actionButton(ns("go"),                "Generate"                                ),
      numericInput(ns("minFeaturesPerSpot"),"Min genes",             200, min = 0, step = 1),
      numericInput(ns("minCountsPerSpot"),  "Min spot UMIs",         500, min = 0, step = 1),
      numericInput(ns("maxPercentMito"),    "Max % mito",            20,  min = 0, max = 100),
      numericInput(ns("nRows"),             "Rows to show in table:", 5,  min = 1, max = 100),
      textOutput(ns("error_msg"))
    ),
    mainPanel(
      navlistPanel(
        id      = ns("spatial_tabs"),
        widths  = c(2, 10),
        tabPanel("Data Summary",   DTOutput( ns("summary_table") )),
        tabPanel("QC Distribution",plotOutput(ns("qc_violin"))),
        tabPanel("QC Scatter",     plotOutput(ns("qc_scatter"))),
        tabPanel("UMAP",           plotOutput(ns("umap_plot")) )
      )
    )
  )
}

# --------------------------
# Server Module
# --------------------------
visiumServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns        <- session$ns
    error_msg <- reactiveVal("Waiting for upload…")
    
    visium_obj <- eventReactive(input$go, {
      zipf <- input$zip
      
      # ── 1) bail if no file uploaded ─────────────────────
      if (is.null(zipf) || nrow(zipf) == 0) {
        error_msg("Please upload a zipped Visium directory first.")
        return(NULL)
      }
      
      # ── 2) unpack into a fresh tempdir ──────────────────
      tmpdir <- file.path(tempdir(), id)
      unlink(tmpdir, recursive = TRUE)
      dir.create(tmpdir, recursive = TRUE)
      unzip(zipf$datapath, exdir = tmpdir)
      
      # ── 3) if there’s one top‐level folder, descend into it ──
      roots    <- list.dirs(tmpdir, full.names = TRUE, recursive = FALSE)
      data_dir <- if (length(roots) == 1) roots else tmpdir
      
      # ── 4) load ────────────────────────────────────────────
      error_msg("Loading Visium data…")
      obj <- tryCatch(
        Load10X_Spatial(data.dir = data_dir),
        error = function(e) {
          error_msg(paste("Error loading data:", e$message))
          NULL
        }
      )
      if (is.null(obj)) return(NULL)
      error_msg("")  # clear
      
      # ── 5) compute percent.mt ────────────────────────────
      obj[["percent.mt"]] <- PercentageFeatureSet(obj,
                                                  pattern = "^MT-",
                                                  assay   = "Spatial")
      
      # ── 6) user‐defined QC filters ───────────────────────
      obj <- subset(obj, subset = nFeature_Spatial >= input$minFeaturesPerSpot)
      counts     <- GetAssayData(obj, assay = "Spatial", slot = "counts")
      keep_spots <- colSums(counts) >= input$minCountsPerSpot
      obj        <- subset(obj, cells = colnames(counts)[keep_spots])
      obj        <- subset(obj, subset = percent.mt <= input$maxPercentMito)
      
      obj
    })
    
    # ── status text ───────────────────────────────────────
    output$error_msg <- renderText(error_msg())
    
    # ── 2) Data summary table ─────────────────────────────
    output$summary_table <- renderDT({
      obj  <- req(visium_obj())
      cnts <- GetAssayData(obj, assay = "Spatial", slot = "counts")
      md   <- obj@meta.data
      df <- data.frame(
        Sample           = unique(md$orig.ident),
        Spots            = ncol(cnts),
        Features         = nrow(cnts),
        MedianFeatures   = round(median(md$nFeature_Spatial), 1),
        MedianCounts     = round(median(md$nCount_Spatial), 1),
        MeanPercent.MT   = round(mean(md$percent.mt), 2),
        `Sparsity (% nz)`= round(100 * sum(cnts != 0) / length(cnts), 2)
      )
      datatable(df, rownames = FALSE,
                options = list(dom = "t",
                               pageLength   = input$nRows,
                               lengthChange = FALSE))
    })
    
    # ── 3) QC violin ───────────────────────────────────────
    output$qc_violin <- renderPlot({
      obj <- req(visium_obj())
      VlnPlot(obj,
              features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
              assay    = "Spatial",
              ncol     = 3)
    })
    
    # ── 4) QC scatter ─────────────────────────────────────
    output$qc_scatter <- renderPlot({
      obj <- req(visium_obj())
      DefaultAssay(obj) <- "Spatial"        # ← make sure Spatial is active
      p1 <- FeatureScatter(obj,
                           feature1 = "nCount_Spatial",
                           feature2 = "percent.mt")
      p2 <- FeatureScatter(obj,
                           feature1 = "nCount_Spatial",
                           feature2 = "nFeature_Spatial")
      wrap_plots(p1, p2)
    })
    
    # ── 5) UMAP ────────────────────────────────────────────
    output$umap_plot <- renderPlot({
      obj <- req(visium_obj())
      obj <- obj |>
        NormalizeData()       |>
        FindVariableFeatures()|>
        ScaleData()           |>
        RunPCA(npcs = 20)     |>
        RunUMAP(dims = 1:20)
      DimPlot(obj, reduction = "umap")
    })
    
  })
}




