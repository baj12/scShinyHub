
myZippedReportFiles <- c("output_topExpGenes.csv")


# All clusters heat map ------
callModule(pHeatMapModule, "coExpHeatmapModule", heatmapReactive)

# 2D plot with selection of cells ------
# assigning it to a variable allows us to interact with the plot and collect selection events
selctedCluster <-
  callModule(
    clusterServer,
    "selected",
    projections,
    reactive(input$gene_id_sch)
  )

# selected clusters heatmap module -----
callModule(
  pHeatMapModule,
  "heatmapSelectedModule",
  heatmapSelectedReactive
)

# max expressed genes ----
callModule(
  tableSelectionServer,
  "topExpGenes",
  topExpGenesTable
)

# SOM heatmap module -----
callModule(
  pHeatMapModule,
  "heatmapSOM",
  heatmapSOMReactive
)

# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
output$geneGrp_vio_plot <- renderPlot({
  start.time <- base::Sys.time()
  if (DEBUG) cat(file = stderr(), "output$geneGrp_vio_plot\n")
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "geneGrp_vio_plot")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("geneGrp_vio_plot", id = "geneGrp_vio_plot", duration = NULL)
  }
  
  projections <- projections()
  scEx <- scEx()
  geneListStr <- input$geneGrpVioIds
  projectionVar <- input$dimension_xVioiGrp
  minExpr <- input$coEminExpr
  showPermutations <- input$showPermutations
  colPal = geneGrp_vioFunc
  sampCol = sampleCols$colPal
  
  upI <- updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$geneGrp_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/geneGrp_vio_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/geneGrp_vio_plot.RData")
  
  featureData <- rowData(scEx)
  retVal <- geneGrp_vioFunc(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx,
    featureData = featureData,
    minExpr = minExpr,
    dbCluster = projectionVar,
    showPermutations = showPermutations,
    sampCol = sampCol
  )
  
  printTimeEnd(start.time, "geneGrp_vio_plot")
  exportTestValues(geneGrp_vio_plot = {retVal})  
  return(retVal)
})
