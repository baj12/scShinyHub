# source("moduleServer.R", local = TRUE)
# source("reactives.R", local = TRUE)

# TODO: verify that this anything and then integrate in DUMMY
myZippedReportFiles <- c("gqcProjections.csv")



gQC_X1 <<- "tsne1"
gQC_X2 <<- "tsne2"
gQC_X3 <<- "tsne3"
gQC_col <<- "sampleNames"
observe({
  gQC_X1 <<- input$gQC_dim3D_x
})
observe({
  gQC_X2 <<- input$gQC_dim3D_y
})
observe({
  gQC_X3 <<- input$gQC_dim3D_z
})
observe({
  gQC_col <<- input$gQC_col3D
})

# gQC_update3DInput ----
#' gQC_update3DInput
#' update axes for tsne display
gQC_update3DInput <- reactive({
  projections <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  # choices = colnames(projections)[unlist(lapply(colnames(projections), function(x) !is.factor(projections[,x])))]
  choices = colnames(projections)
  # Can also set the label and select items
  updateSelectInput(session, "gQC_dim3D_x",
                    choices = choices,
                    selected = gQC_X1
  )
  
  updateSelectInput(session, "gQC_dim3D_y",
                    choices = choices,
                    selected = gQC_X2
  )
  updateSelectInput(session, "gQC_dim3D_z",
                    choices = choices,
                    selected = gQC_X3
  )
  updateSelectInput(session, "gQC_col3D",
                    choices = colnames(projections),
                    selected = gQC_col
  )
})

# gQC_tsne_main ----
output$gQC_tsne_main <- renderPlotly({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "gQC_tsne_main")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gQC_tsne_main", id = "gQC_tsne_main", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$gQC_tsne_main\n")
  
  upI <- gQC_update3DInput()
  projections <- projections()
  dimX <- input$gQC_dim3D_x
  dimY <- input$gQC_dim3D_y
  dimZ <- input$gQC_dim3D_z
  dimCol <- input$gQC_col3D
  scols <- sampleCols$colPal
  
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$gQC_tsne_main:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gQC_tsne_main.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gQC_tsne_main.RData")
  
  retVal <- tsnePlot(projections, dimX, dimY, dimZ, dimCol, scols)
  
  printTimeEnd(start.time, "tsnePlot")
  exportTestValues(tsnePlot = {str(retVal)})  
  return(layout(retVal))
})

#' tsnePlot
#' function that plots in 3D the tsne projection
tsnePlot <- function(projections, dimX, dimY, dimZ, dimCol, scols) {
  projections <- as.data.frame(projections)
  projections$dbCluster <- as.factor(projections$dbCluster)
  
  if (dimCol == "sampleNames") {
    myColors <- scols
  } else {
    myColors <- NULL
  }
  
  p <-
    plot_ly(
      projections,
      x = formula(paste("~ ", dimX)),
      y = formula(paste("~ ", dimY)),
      z = formula(paste("~ ", dimZ)),
      type = "scatter3d",
      color = formula(paste("~ ", dimCol)),
      colors = myColors,
      hoverinfo = "text",
      text = paste("Cluster:", as.numeric(as.character(projections$dbCluster))),
      mode = "markers",
      marker =
        list(
          line = list(width = 0),
          size = rep(10, nrow(projections)),
          sizeref = 3
        )
    )
  return(p)
}

# gQC_umap_main 2D plot ----
callModule(
  clusterServer,
  "gQC_umap_main",
  projections
)

# gQC_projectionTableMod ----
callModule(
  tableSelectionServer, 
  "gQC_projectionTableMod", 
  projectionTable)

# gQC_plotUmiHist ----
output$gQC_plotUmiHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output_gQC_plotUmiHist\n")
  scEx <- scEx()
  scols <- sampleCols$colPal
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gQC_plotUmiHist.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/gQC_plotUmiHist.RData")

  dat <- data.frame(counts = Matrix::colSums(assays(scEx)[["counts"]]))
  dat$sample <- colData(scEx)$sampleNames
  ggplot(data = dat, aes(counts, fill = sample)) +
    geom_histogram(bins = 50) +
    labs(title = "Histogram for raw counts", x = "count", y = "Frequency") +
    scale_fill_manual(values = scols, aesthetics = "fill")
  
})

output$gQC_plotSampleHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output_sampleHist\n")
  sampleInf <- sampleInfo()
  scols <- sampleCols$colPal
  
  if (is.null(sampleInf)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/sampleHist.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/sampleHist.RData")
  gQC_sampleHistFunc(sampleInf, scols)
})

output$gQC_variancePCA <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output$gQC_variancePCA\n")
  h2("Variances of PCs")
  pca <- pca()
  if (is.null(pca)) {
    return(NULL)
  }
  barplot(pca$var_pcs, main = "Variance captured by first PCs")
})
