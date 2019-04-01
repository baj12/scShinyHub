
myZippedReportFiles <- c("output_topExpGenes.csv")


# updateInputXviolinPlot ---------
# Update x axis selection possibilities for violin plot
updateInputXviolinPlot <- reactive({
  tsneData <- projections()

  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }

  # Can also set the label and select items
  updateSelectInput(
    session,
    "dimension_x3",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[1]
  )
  updateSelectInput(
    session,
    "dimension_y3",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[2]
  )

  coln <- colnames(tsneData)
  choices <- c()
  for (cn in coln) {
    if (length(levels(as.factor(tsneData[, cn]))) < 20) {
      choices <- c(choices, cn)
    }
  }
  if (length(choices) == 0) {
    choices <- c("no valid columns")
  }
  updateSelectInput(
    session,
    "dimension_xVioiGrp",
    choices = choices,
    selected = choices[1]
  )
})


# heatmapReactive -------
# reactive for module pHeatMapModule
# for all clusters menu item
heatmapReactive <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "heatmap")
  )
  if (DEBUG) cat(file = stderr(), "output$heatmap\n")
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$heatmap_geneids
  sampCol = sampleCols$colPal
  
  if (is.null(scEx_log) | is.null(projections)) {
    return(list(
      src = "empty.png",
      contentType = "image/png",
      width = 96,
      height = 96,
      alt = "heatmap should be here"
    ))
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("heatmap", id = "heatmap", duration = NULL)
  }

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/heatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/heatmap.RData")
  featureData <- rowData(scEx_log)
  scEx_matrix <- as.matrix(assays(scEx_log)[["logcounts"]])
  retval <- coE_heatmapFunc(
    featureData = featureData, scEx_matrix = scEx_matrix,
    projections = projections, genesin = genesin, cells = colnames(scEx_matrix),
    sampCol = sampCol
  )

  return(retval)
})



# All clusters heat map ------
callModule(pHeatMapModule, "coExpHeatmapModule", heatmapReactive)




# 2D plot with selection of cells ------
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
  if (DEBUG) {
    cat(file = stderr(), "output$geneGrp_vio_plot\n")
  }
  # if (v$doPlot == FALSE)
  #   return()
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
    if (DEBUG) {
      cat(file = stderr(), "output$geneGrp_vio_plot:NULL\n")
    }
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
  if (DEBUG) {
    cat(file = stderr(), "output$plotCoExpression:done\n")
  }
  retVal
})
