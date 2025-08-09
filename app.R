
# allow uploads up to 5 GB (5 × 1024^3 bytes)
options(shiny.maxRequestSize = 5 * 1024^3)




#Testing_1

# app.R
library(shiny)

# — Make sure Bioc repos are set —
if (requireNamespace("BiocManager", quietly=TRUE)) {
  options(repos = BiocManager::repositories())
}


# — Allow large uploads —
#options(shiny.maxRequestSize = 3000 * 1024^2)

# — Source modules —
source("scRNASeq_module.R")
source("visium_module.R")
source("scRNASeq_merge_module.R")
source("visium_merge_module.R")
source("scRNASeq_visual_module.R")
source("visium_visual_module.R")
source("scRNASeq_cloupe_module.R")
source("visium_cloupe_module.R")






ui <- fluidPage(
  
  tags$h2("Single‐Cell & Spatial Transcriptomics Viewer",
          style="font‐size:26px; font‐weight:bold; text‐align:center; margin‐top:10px;"),
  tabsetPanel(
    tabPanel("Data",
             tabsetPanel(
               tabPanel("scRNA‐seq Analysis",
                        fileInput("scrna_files",
                                  "Upload 10X scRNA‐seq files (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)",
                                  multiple = TRUE,
                                  accept   = c(".mtx.gz",".tsv.gz",".rds")),
                        scRNASeqUI("scRNA")
               ),
               tabPanel("Visium Analysis",
                        # fileInput(
                        #   "visium_files",
                        #   "Upload your entire Visium directory as a .zip",
                        #   multiple = FALSE,
                        #   accept   = ".zip"
                        # ),
                        visiumUI("visium")    
               )
               
             )
    ),
    
    
    # ── MERGE TAB 
    tabPanel("Merge",
             tabsetPanel(
               tabPanel("Merge scRNA-seq", scRNASeqMergeUI("merge")),
               tabPanel("Merge Visium",    visiumMergeUI("visium_merge"))
             )
    ),
    
    # ── VISUALIZATION TAB 
    tabPanel("Visualization",
             tabsetPanel(
               tabPanel("scRNA-seq",
                        fluidRow(
                          column(3,
                                 fileInput("vis_scrna_files",
                                           "Upload processed scRNA-seq object (.rds)",
                                           accept=".rds")
                          ),
                          column(9, scRNASeqVisualUI("scrna_vis"))
                        )
               ),
               tabPanel("Visium",
                        fluidRow(
                          column(3,
                                 fileInput("vis_visium_files",
                                           "Upload processed Visium object (.rds)",
                                           accept=".rds")
                          ),
                          column(9, visiumVisualUI("visium_vis"))
                        )
               )
             )
    ),
    
    # ── EXPORT TAB 
    tabPanel("Export",
             tabsetPanel(
               tabPanel("scRNA-seq",
                        fileInput("export_scrna_files",
                                  "Upload scRNA-seq export files (.cloupe, .rds)",
                                  multiple=TRUE,
                                  accept=c(".cloupe",".rds")),
                        scRNASeqCloupeUI("cloupe")
               ),
               tabPanel("Visium",
                        fileInput("export_visium_files",
                                  "Upload Visium export files (.cloupe, .rds)",
                                  multiple=TRUE,
                                  accept=c(".cloupe",".rds")),
                        visiumCloupeUI("visium_cloupe_ui")
               )
             )
             
             
    )
    
    
  )
)


# server part

server <- function(input, output, session) {
  # scRNA
  scRNASeqServer("scRNA", files = reactive(input$scrna_files))
  
  # Visium 
  # visiumServer("visium", files = reactive(input$visium_files))
  visiumServer("visium")
  
  
  
  # — MERGE MODULES —
  scRNASeqMergeServer("merge")
  visiumMergeServer("visium_merge")
  
  
  # — VISUALIZATION MODULES —
  scRNASeqVisualServer("scrna_vis", file = reactive(input$vis_scrna_files))
  visiumVisualServer  ("visium_vis", file = reactive(input$vis_visium_files))
  
  
  # — CLOUPE / EXPORT MODULES —
  scRNASeqCloupeServer("cloupe", files = reactive(input$export_scrna_files))
  visiumCloupeServer  (
    id    = "visium_cloupe_ui",
    files = reactive(input$export_visium_files)
  )
  
  
}




# finally launch
shinyApp(ui, server)



