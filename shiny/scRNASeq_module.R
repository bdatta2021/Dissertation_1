









 # Edition_2:
  
  
  # Trying to add more subtabs into Data tab for "scRNASeq Analysis " tab
  
  
  library(shiny)
  library(Seurat)
  library(DT)
  
  # --------------------------
  # UI Module
  # --------------------------
  scRNASeqUI <- function(id) {
    ns <- NS(id)
    tagList(
      sidebarLayout(
        sidebarPanel(
          actionButton(ns("go"), "Generate"),      # ← only the Generate button
          
          ## ─── YOUR THREE SPINNERS ───
          numericInput(
            inputId = ns("minCellsPerGene"),
            label   = "Minimum Cells per Gene:",
            value   = 3,
            min     = 0,
            step    = 1
          ),                                        # ← NEW
          numericInput(
            inputId = ns("minFeaturesPerCell"),
            label   = "Minimum Features per Cell:",
            value   = 200,
            min     = 0,
            step    = 1
          ),                                        # ← NEW
          numericInput(
            inputId = ns("nRows"),
            label   = "Rows to show in table:",
            value   = 5,
            min     = 1,
            max     = 100,
            step    = 1
          ),                                        # ← NEW
          ## ────────────────────────────
          
          textOutput(ns("error_message"))
        ),
       
        #  mainPanel(
        #  plotOutput(ns("umap_plot")),
        #  DTOutput(ns("summary_table"))
        #  )
        
        
        mainPanel(
          
          tabsetPanel(
            
            id = ns("qc_tabs"),
            
            
            
            tabPanel(
              
              title = "Data Summary",
              
              DTOutput(ns("summary_table"))
            ),
            
            
            
            tabPanel(
              
              title = "QC Distribution",
              
              plotOutput(ns("qc_violin"))
              
            ),
            
            
            tabPanel(
              
              title = "QC Scatter",
              
              plotOutput(ns("qc_scatter"))
              
            ),
            
            
            tabPanel(
              
              title = "UMAP",
              
              plotOutput(ns("umap_plot"))
              
            )
            
                      )
        
          )
        
        
      )
    )
  }
  
  
  
  
  # --------------------------
  # Server Module
  # --------------------------
  scRNASeqServer <- function(id, files) {
    
    moduleServer(id, function(input, output, session) {
      
      ns <- session$ns
      
      
      
      # A little message for the user
      
      error_msg <- reactiveVal("Waiting to generate…")
      
      
      
      # Load and build the suerat object only when they click “Generate”
      
      
      seurat_obj <- eventReactive(input$go, {
        
        tbl <- files()
        
        if (is.null(tbl) || nrow(tbl) == 0) {
          
          error_msg("Please upload your three 10X files first.")
          
          return(NULL)
          
        }
        
        
        
        # checking the presence of the three core filenames
        
        
        required <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
        
        missing  <- setdiff(required, tbl$name)
        
        if (length(missing)) {
          
          error_msg(paste("Missing file(s):", paste(missing, collapse = ", ")))
          
          return(NULL)
          
        }
        
        
        
        # copy into a clean temp folder for Read10X()
        
        tmpdir <- file.path(tempdir(), id)
        
        unlink(tmpdir, recursive = TRUE)
        
        dir.create(tmpdir, recursive = TRUE)
        
        file.copy(tbl$datapath, file.path(tmpdir, tbl$name))
        
        
        
        error_msg("Loading data…")
        
        tryCatch({
          
          counts <- Read10X(data.dir = tmpdir)
          
          # ─── USE YOUR SPINNER VALUES HERE ───
          
          obj <- CreateSeuratObject(
            
            counts       = counts,
            
            min.cells    = input$minCellsPerGene,     # Spinner #1
            
            min.features = input$minFeaturesPerCell   # Spinner #2
            
          )
          
          
          obj[["percent.mt"]] <- PercentageFeatureSet(
            
            obj,
            
            pattern = "^MT-"
            
            
          )
        
            error_msg("")  # clear any error
          
          obj
          
        }, error = function(e) {
          
          error_msg(paste("Error loading data:", conditionMessage(e)))
          
          NULL
          
        })
        
      })
      
      
      # render status
      output$error_message <- renderText(error_msg())
      
      
      
      # UMAP plot (this drives the "umap" tab)
      
      output$umap_plot <- renderPlot({
        
        obj <- seurat_obj()
        
        req(obj)      # Nothing renders until we have an object
        
        obj <- obj |>
          
          NormalizeData() |>
          
          FindVariableFeatures() |>
          
          ScaleData() |>
          
          RunPCA() |>
          
          RunUMAP(dims = 1:10)
        
        DimPlot(obj, reduction = "umap", label = TRUE)
        
      })
      
      
      
      
      # Data summary table - (this drives the "Data summary tab")
      
      
      output$summary_table <- renderDT({
        
        obj <- seurat_obj()
        
        req(obj)
        
        
        
        # Pull out counts + meta data
        
        cnts <- GetAssayData(obj, slot = "counts")
        
        md <- obj@meta.data
        
        
        
        # Build richer a one-row summary
        
        df <- data.frame(
          
          `Sample`  = unique(md$orig.ident),
          
          `Cells`   = ncol(cnts),
          
          `Features`  = nrow(cnts),
          
          `MedianFeatures` = round(median(md$nFeature_RNA),1),
          
          `MedianCounts`   = round(median(md$nCount_RNA),1),
          
          `MeanPercent.MT` = round(mean(md$percent.mt),2),
          
          
          `Sparsity (% non-zero)` = round(100 * sum(cnts != 0) / length(cnts), 2)
          
          
        )
        
        
        datatable(df,
                  
          options = list(
            
            pageLength  = input$nRows,   # "Rows to show " spinner # 3
            
            lengthChange = FALSE
            
          )
          
        )
        
      })
      
    
    
  
  
  
  
  
  # QC distribution violins (for "QC Distributions" tab)
  
  
  output$qc_violin <- renderPlot({
    
    obj <- seurat_obj()
    
    req(obj)
    
    VlnPlot(
      
      obj,
      
      features = c("nFeature_RNA","nCount_RNA","percent.mt"),
      
      ncol = 3
      
    )
    
  })
  
  
  # QC sacatter plots (for "QC scatter " tab)
  
  
  
  output$qc_scatter <- renderPlot({
    
 obj <- seurat_obj()
 
 req(obj)
 
 p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
 
 p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
 
 patchwork::wrap_plots(p1, p2)  # puts them side by side
 
 
 })   

  
                                  
})    # closing the moduleServer
  

}
  
  
  
