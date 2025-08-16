







# -------------------
# UI Module (no longer needs fileInput)
# -------------------
scRNASeqVisualUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(ns("reduction"), "Select Reduction", choices = NULL),
    selectInput(ns("color_by"),  "Color By",       choices = NULL),
    actionButton(ns("plot_btn"), "Generate Plot"),
    plotOutput(ns("plot")),
    DTOutput(ns("meta_table"))
  )
}

# -------------------
# Server Module
# -------------------
scRNASeqVisualServer <- function(id, file) {
  moduleServer(id, function(input, output, session) {
    # when the user uploads:
    obj <- eventReactive(file(), {
      req(file())
      readRDS(file()$datapath)
    })

    # once we have an object, populate the two dropdowns
    observeEvent(obj(), {
      md <- obj()@meta.data
      updateSelectInput(session, "color_by", choices     = colnames(md))
      updateSelectInput(session, "reduction", choices   = names(obj()@reductions))
    })

    # render the plot
    output$plot <- renderPlot({
      req(obj(), input$color_by, input$reduction, input$plot_btn)
      DimPlot(
        obj(),
        reduction = input$reduction,
        group.by  = input$color_by,
        label     = TRUE
      )
    })

    # render the metadata table
    output$meta_table <- renderDT({
      req(obj())
      datatable(
        obj()@meta.data,
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })
  })
}




