source("moduleServer.R", local = TRUE)
source("reactives.R", local = TRUE)


myZippedReportFiles <- c("gqcProjections.csv")


update3DInput <- reactive({
  tsneData <- projections()

  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }

  # Can also set the label and select items
  updateSelectInput(session, "dim3D_x",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[1]
  )

  updateSelectInput(session, "dim3D_y",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[2]
  )
  updateSelectInput(session, "dim3D_z",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[3]
  )
  updateSelectInput(session, "col3D",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[3]
  )
})


output$tsne_main <- renderPlotly({
  upI <- update3DInput()
  if (DEBUG) cat(file = stderr(), "output$tsne_main\n")
  projections <- projections()
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$tsne_main:NULL\n")
    return(NULL)
  }
  dimX <- input$dim3D_x
  dimY <- input$dim3D_y
  dimZ <- input$dim3D_z
  dimCol <- input$col3D


  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne_main.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/tsne_main.RData")

  projections <- as.data.frame(projections)
  # cat(stderr(),colnames(projections)[1:5])
  projections$dbCluster <- as.factor(projections$dbCluster)

  p <-
    plot_ly(
      projections,
      x = formula(paste("~ ", dimX)),
      y = formula(paste("~ ", dimY)),
      z = formula(paste("~ ", dimZ)),
      type = "scatter3d",
      color = formula(paste("~ ", dimCol)),
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
  # layout(p)
  if (DEBUG) cat(file = stderr(), "output$tsne_main: done\n")
  return(layout(p))
})

# output$umap_main <- renderPlotly({
#   embedding <- umapReact()
#   gbmlog <- gbm_log()
#   pointSize <- 1
#   if (is.null(embedding)) {
#     if (DEBUG) cat(file = stderr(), "output$umap_main:NULL\n")
#     return(NULL)
#   }
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/umap_main.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   
#   # outTab$UMAP1 = embedding$UMAP1
#   # outTab$UMAP2 = embedding$UMAP2
#   embedding %>%
#      ggplot(aes_string(UMAP1, UMAP2)) + geom_point(size = pointSize)
# })

  callModule(
    clusterServer,
    "umap_main",
    projections
    # ,
    # defaultValues = c("UMAP1", "UMAP2")
  )



r <- callModule(tableSelectionServer, "cellSelectionTSNEMod", inputTSNESample)


output$plotUmiHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output_plotUmiHist\n")
  gbm <- gbm()
  if (is.null(gbm)) {
    return(NULL)
  }
  hist(Matrix::colSums(exprs(gbm)), breaks = 50, main = "histogram of number of UMIs per cell")
})

output$plotSampleHist <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output_sampleHist\n")
  sampleInf <- sampleInfo()
  if (is.null(sampleInf)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/sampleHist.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/sampleHist.RData")
  sampleHistFunc(sampleInf)
})

output$variancePCA <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output$variancePCA\n")
  h2("hello")
  pca <- pca()
  if (is.null(pca)) {
    return(NULL)
  }
  barplot(pca$var_pcs, main = "Variance captured by first PCs")
})
