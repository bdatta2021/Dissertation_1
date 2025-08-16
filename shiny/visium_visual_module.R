

# visium_visual_module.R

library(shiny)
library(Seurat)
library(DT)

# — UI Module ——————————————————————————————————————————————————————
visiumVisualUI <- function(id) {
  ns <- NS(id)
  tagList(
    # no fileInput here anymore – it's in app.R
    selectInput(ns("reduction"), "Select Reduction", choices = NULL),
    selectInput(ns("group_by"),  "Color By",       choices = NULL),
    actionButton(ns("plot_btn"), "Generate Plot"),
    plotOutput(ns("dim_plot")),
    plotOutput(ns("spatial_plot"))
  )
}

# — Server Module ———————————————————————————————————————————————————
# note the added `file` parameter
visiumVisualServer <- function(id, file) {
  moduleServer(id, function(input, output, session) {
    #  wait for the uploaded .rds
    obj <- eventReactive(file(), {
      req(file())
      readRDS(file()$datapath)
    })

    #  once loaded, populate the two dropdowns
    observeEvent(obj(), {
      md  <- obj()@meta.data
      reds <- names(obj()@reductions)
      updateSelectInput(session, "group_by", choices   = colnames(md))
      updateSelectInput(session, "reduction", choices  = reds)
    })

    #  render the UMAP/PCA/tSNE plot
    output$dim_plot <- renderPlot({
      req(input$plot_btn, obj(), input$reduction, input$group_by)
      DimPlot(
        obj(),
        reduction = input$reduction,
        group.by  = input$group_by,
        label     = TRUE
      )
    })

    #  render the spatial‐coordinates plot
    output$spatial_plot <- renderPlot({
      req(input$plot_btn, obj(), input$group_by)
      SpatialDimPlot(
        obj(),
        group.by = input$group_by,
        label    = TRUE
      )
    })
  })
}


