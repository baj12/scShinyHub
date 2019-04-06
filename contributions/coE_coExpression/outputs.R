
myZippedReportFiles <- c("output_coE_topExpGenes.csv")


# All clusters heat map ------
callModule(
  pHeatMapModule, 
  "coExpHeatmapModule", 
  coE_heatmapReactive)

# 2D plot with selection of cells ------
# assigning it to a variable allows us to interact with the plot and collect selection events
coE_selctedCluster <-
  callModule(
    clusterServer,
    "coE_selected",
    projections,
    reactive(input$coE_gene_id_sch)
  )

# selected clusters heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSelectedModule",
  coE_heatmapSelectedReactive
)

# max expressed genes ----
callModule(
  tableSelectionServer,
  "coE_topExpGenes",
  coE_topExpGenesTable
)

# SOM heatmap module -----
callModule(
  pHeatMapModule,
  "coE_heatmapSOM",
  coE_coE_heatmapSOMReactive
)

# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
output$coE_geneGrp_vio_plot <- renderPlot({
  start.time <- base::Sys.time()
  if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot\n")
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "coE_geneGrp_vio_plot")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("coE_geneGrp_vio_plot", id = "coE_geneGrp_vio_plot", duration = NULL)
  }
  
  projections <- projections()
  scEx <- scEx()
  geneListStr <- input$coE_geneGrpVioIds
  projectionVar <- input$coE_dimension_xVioiGrp
  minExpr <- input$coEminExpr
  coE_showPermutations <- input$coE_showPermutations
  colPal = coE_geneGrp_vioFunc
  sampCol = sampleCols$colPal
  
  upI <- coE_updateInputXviolinPlot() # no need to check because this is done in projections
  if (is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$coE_geneGrp_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/coE_geneGrp_vio_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/coE_geneGrp_vio_plot.RData")
  
  featureData <- rowData(scEx)
  retVal <- coE_geneGrp_vioFunc(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx,
    featureData = featureData,
    minExpr = minExpr,
    dbCluster = projectionVar,
    coE_showPermutations = coE_showPermutations,
    sampCol = sampCol
  )
  
  printTimeEnd(start.time, "coE_geneGrp_vio_plot")
  exportTestValues(coE_geneGrp_vio_plot = {retVal})  
  return(retVal)
})
