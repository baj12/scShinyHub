
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
  featureData <- featureDataReact()
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$heatmap_geneids
  sampCol = sampleCols$colPal
  
  if (is.null(featureData) | is.null(scEx_log) | is.null(projections)) {
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
  scEx_matrix <- assays(scEx_log)[[1]]
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



# plotCoExpression ----
# binarized 2D plot
# TODO module?
# output$plotCoExpression <- renderPlot({
#   if (DEBUG) {
#     cat(file = stderr(), "output$plotCoExpression\n")
#   }
#   # if (vvvvvvv$doPlot == FALSE)
#   #   return()
#   featureData <- featureDataReact()
#   scEx_log <- log2cpm()
#   upI <- updateInputXviolinPlot() # no need to check because this is done in projections
#   projections <- projections()
#   if (is.null(featureData) |
#     is.null(projections) |
#     is.null(log2cpm) | is.null(input$clusters3)) {
#     return(NULL)
#   }
# 
#   genesin <- input$mclustids
#   cl3 <- input$clusters3
#   dimx3 <- input$dimension_x3
#   dimy3 <- input$dimension_y3
#   # posCells <- positiveCells$positiveCells # we use this variable to be able to save the global variable in this context
#   # posCellsAll <- positiveCells$positiveCellsAll
# 
# 
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/plotCoExpression.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/scShinyHubDebug/plotCoExpression.RData")
#   p1 <- plotCoExpressionFunc(
#     featureData,
#     scEx_log,
#     upI,
#     projections,
#     genesin,
#     cl3,
#     dimx3,
#     dimy3
#   )
#   return(p1)
# })
# 

# downloadExpressionOnOff -----
# binarized
# TODO module download
# output$downloadExpressionOnOff <- downloadHandler(
#   filename = function() {
#     paste(input$clusters3, "PositiveCells.csv", sep = "_")
#   },
#   content = function(file) {
#     featureData <- featureDataReact()
#     log2cpm <- log2cpm()
#
#     if (is.null(featureData) | is.null(log2cpm) | is.null(positiveCells$positiveCells)) {
#       return(NULL)
#     }
#
#     cells <- positiveCells$positiveCells
#     # if(DEBUG)cat(file=stderr(),cells[1:5])
#
#     if (length(cells) == 1) {
#       subsetExpression <- log2cpm[, cells]
#       subsetExpression <-
#         as.data.frame(subsetExpression, row.names = rownames(log2cpm))
#       colnames(subsetExpression) <- cells
#       subsetExpression$Associated.Gene.Name <-
#         featureData[rownames(subsetExpression), "Associated.Gene.Name"]
#       write.csv(subsetExpression, file)
#     }
#     else {
#       subsetExpression <- log2cpm[, cells]
#       # cat(stderr(),colnames(subsetExpression)[1:5])
#       subsetExpression$Associated.Gene.Name <-
#         featureData[rownames(subsetExpression), "Associated.Gene.Name"]
#       # cat(stderr(),colnames(subsetExpression))
#       write.csv(subsetExpression, file)
#     }
#   }
# )

# binarized on / off table ------
# TODO do we need it?
# output$onOffTable <- DT::renderDataTable({
#   if (DEBUG) {
#     cat(file = stderr(), "output$onOffTable\n")
#   }
#   projections <- projections()
#   posCellsAll <- positiveCells$positiveCellsAll # we use this variable to be able to save the global variable in this context
#
#   if (is.null(projections)) {
#     return(NULL)
#   }
#
#
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/onOffTable.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/scShinyHubDebug/onOffTable.RData")
#
#   merge <- projections
#   if (DEBUG) {
#     cat(
#       file = stderr(),
#       paste("positiveCells$positiveCellsAll:---", posCellsAll, "---\n")
#     )
#   }
#
#   merge$CoExpression <- posCellsAll
#   df <-
#     as.data.frame(table(merge[, c("dbCluster", "CoExpression")]))
#   dfOut <- cast(df, dbCluster ~ CoExpression)
#   colnames(dfOut) <- c("Cluster", "OFF", "ON")
#   rownames(dfOut) <- dfOut$Cluster
#   dfOut["Sum", ] <- c("", sum(dfOut$OFF), sum(dfOut$ON))
#   if (DEBUG) {
#     cat(file = stderr(), "output$onOffTable:done\n")
#   }
#   DT::datatable(dfOut)
# })

# TODO as module
# coexpression binarized
# output$clusters3 <- renderUI({
#   if (DEBUG) {
#     cat(file = stderr(), "output$clusters3\n")
#   }
#   projections <- projections()
#   if (is.null(projections)) {
#     HTML("Please load data first")
#     return(NULL)
#   }
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/clusters3.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file="~/scShinyHubDebug/clusters3.RData")
# 
#   noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
#   selectizeInput(
#     "clusters3",
#     label = "Cluster",
#     choices = c(1:noOfClusters),
#     selected = 1,
#     multiple = TRUE
#   )
# })

# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# TODO module for violin plot  ??
output$geneGrp_vio_plot <- renderPlot({
  if (DEBUG) {
    cat(file = stderr(), "output$geneGrp_vio_plot\n")
  }
  # if (v$doPlot == FALSE)
  #   return()
  featureData <- featureDataReact()
  projections <- projections()
  scEx <- scEx()
  geneListStr <- input$geneGrpVioIds
  projectionVar <- input$dimension_xVioiGrp
  minExpr <- input$coEminExpr
  showPermutations <- input$showPermutations
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

  retVal <- geneGrp_vioFunc(
    genesin = geneListStr,
    projections = projections,
    scEx = scEx,
    featureData = featureData,
    minExpr = minExpr,
    dbCluster = projectionVar,
    showPermutations = showPermutations
  )
  if (DEBUG) {
    cat(file = stderr(), "output$plotCoExpression:done\n")
  }
  retVal
})
