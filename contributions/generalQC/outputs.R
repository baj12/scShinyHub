# source("moduleServer.R", local = TRUE)
# source("reactives.R", local = TRUE)

# TODO: verify that this anything and then integrate in DUMMY
myZippedReportFiles <- c("gqcProjections.csv")

# update3DInput ----
# TODO see module on how to remember selections
#' update3DInput
#' update axes for tsne display
update3DInput <- reactive({
  projections <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "dim3D_x",
                    choices = colnames(projections),
                    selected = colnames(projections)[1]
  )
  
  updateSelectInput(session, "dim3D_y",
                    choices = colnames(projections),
                    selected = colnames(projections)[2]
  )
  updateSelectInput(session, "dim3D_z",
                    choices = colnames(projections),
                    selected = colnames(projections)[3]
  )
  updateSelectInput(session, "col3D",
                    choices = colnames(projections),
                    selected = colnames(projections)[3]
  )
})

# tsne_main ----
output$tsne_main <- renderPlotly({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "tsne_main")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne_main", id = "tsne_main", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$tsne_main\n")
  
  upI <- update3DInput()
  projections <- projections()
  dimX <- input$dim3D_x
  dimY <- input$dim3D_y
  dimZ <- input$dim3D_z
  dimCol <- input$col3D
  scols <- sampleCols$colPal
  
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$tsne_main:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne_main.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/tsne_main.RData")
  
  retVal <- tsnePlot(projections, dimX, dimY, dimZ, dimCol, scols)
  
  printTimeEnd(tsnePlot, "tsnePlot")
  exportTestValues(tsnePlot = {str(retVal)})  
  return(layout(retVal))
})

#' tsnePlot
#' function that plots in 3D the tsne projection
tsnePlot <- function() {
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

# umap_main 2D plot ----
callModule(
  clusterServer,
  "umap_main",
  projections
)

# projectionTableMod ----
callModule(
  tableSelectionServer, 
  "projectionTableMod", 
  projectionTable)

# plotUmiHist ----
output$plotUmiHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output_plotUmiHist\n")
  scEx <- scEx()
  scols <- sampleCols$colPal
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/plotUmiHist.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/plotUmiHist.RData")

  dat <- data.frame(counts = Matrix::colSums(assays(scEx)[["counts"]]))
  dat$sample <- colData(scEx)$sampleNames
  ggplot(data = dat, aes(counts, fill = sample)) +
    geom_histogram(bins = 50) +
    labs(title = "Histogram for raw counts", x = "count", y = "Frequency") +
    scale_fill_manual(values = scols, aesthetics = "fill")
  
})

output$plotSampleHist <- renderPlot({
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
  sampleHistFunc(sampleInf, scols)
})

output$variancePCA <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output$variancePCA\n")
  h2("Variances of PCs")
  pca <- pca()
  if (is.null(pca)) {
    return(NULL)
  }
  barplot(pca$var_pcs, main = "Variance captured by first PCs")
})
