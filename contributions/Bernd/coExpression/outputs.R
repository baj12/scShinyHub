updateInputx3 = reactive({
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
  
  coln = colnames(tsneData)
  choices = c()
  for(cn in coln){
    if(length(levels(as.factor(tsneData[,cn]))) < 20)
      choices = c(choices, cn)
  }
  if(length(choices)==0){
    choices = c("no valid columns")
  }
  updateSelectInput(
    session,
    "dimension_xVioiGrp",
    choices = choices,
    selected = choices[1]
  )
})

heatmapFunc <- function(featureData, gbm_matrix, projections, genesin, cells){
  genesin = geneName2Index(genesin, featureData)
  # genesin <- toupper(genesin)
  # genesin <- gsub(" ", "", genesin, fixed = TRUE)
  # genesin <- strsplit(genesin, ',')
  # 
  # map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
  # cat(file = stderr(), length(map))
  # if(!is.null(getDefaultReactiveDomain())){
  #   removeNotification( id="heatmapWarning")
  #   removeNotification( id="heatmapNotFound")
  # }
  # if(length(genesin) == 0){
  #   if(!is.null(getDefaultReactiveDomain())){
  #     showNotification("no genes found", id = "heatmapWarning", type = "warning", duration = 20)
  #   }
  #   return(list(src = "empty.png",
  #        contentType = 'image/png',
  #        width = 96,
  #        height = 96,
  #        alt = "heatmap should be here")
  #        )
  # }
  expression <- gbm_matrix[genesin, cells]
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or genes not expressed'
  ))
  
  # # display genes not found
  # notFound = genesin[[1]][which(!genesin[[1]] %in% featureData$Associated.Gene.Name)]
  # if(length(notFound)>0){
  #   if(!is.null(getDefaultReactiveDomain())){
  #     showNotification(paste("following genes were not found", notFound,collapse = " "), id="heatmapNotFound", type = "warning", duration = 20)
  #   }
  # }

  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]
  
  # expression <- expression[, rownames(projections)]
  expression <- expression[complete.cases(expression), ]
  
  if(!("sample" %in% colnames(projections))){
    projections$sample=1
  }
  annotation <- data.frame(projections[cells, c("dbCluster", "sample")])
  rownames(annotation) <- colnames(expression)
  colnames(annotation) <- c('Cluster', 'sample')
  
  # For high-res displays, this will be greater than 1
  pixelratio <- session$clientData$pixelratio
  if(is.null(pixelratio)) pixelratio = 1
  width  <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  if (is.null(width)) {width = 96*7} # 7x7 inch output
  if (is.null(height)) {height = 96*7}
  outfile <- paste0(tempdir(), '/heatmap', base::sample(1:10000, 1), '.png')
  cat(file = stderr(), paste("saving to: ", outfile, '\n'))
  # this can fail with na/inf in hclust error message if there is a row with all the same values
  # med = median(as.vector(as.matrix(expression)))
  # stDev = sd(as.vector(as.matrix(expression)))
  # minBreak = max(0, med - 3* stDev)
  # maxBreak = med + 3* stDev
  # stepBreak = (maxBreak - minBreak) / 6
  pheatmap(
      as.matrix(expression)[,order(annotation[,1], annotation[,2])],
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      scale = 'row',
      fontsize_row = 10,
      labels_col = colnames(expression),
      labels_row = featureData[rownames(expression), 'Associated.Gene.Name'],
      show_rownames = TRUE,
      annotation_col = annotation,
      show_colnames = FALSE,
      annotation_legend = TRUE,
      # breaks = seq(minBreak, maxBreak, by = stepBreak),
      # filename = 'test.png',
      filename = normalizePath(outfile),
      colorRampPalette(rev(brewer.pal(
        n = 6, name =
          "RdBu"
      )))(6)
      
    )
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification( id = "heatmap")
  }
  return(list(src = normalizePath(outfile),
              contentType = 'image/png',
              width = width,
              height = height,
              alt = "heatmap should be here"))
  
}

# TODO mnodule for heatmap?  
output$heatmap <- renderImage({
  if(DEBUG)cat(file=stderr(), "output$heatmap\n")
  featureData = featureDataReact()
  gbm_log = gbm_log()
  projections = projections()
  genesin <- input$heatmap_geneids
  if(is.null(featureData) | is.null(gbm_log) | is.null(projections)){
    return(list(src = "empty.png",
                contentType = 'image/png',
                width = 96,
                height = 96,
                alt = "heatmap should be here")
    )
  }
  if(!is.null(getDefaultReactiveDomain())){
    showNotification("heatmap", id="heatmap", duration = NULL)
  }
  
  if(DEBUGSAVE)
    save(file = "~/scShinyHubDebug/heatmap.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file = "~/scShinyHubDebug/heatmap.RData")
  gbm_matrix = as.matrix(exprs(gbm_log))
   retval = heatmapFunc(featureData = featureData, gbm_matrix = gbm_matrix, 
                        projections = projections, genesin = genesin, cells = colnames(gbm_matrix))
    
   if(!is.null(getDefaultReactiveDomain())){
      removeNotification( id="heatmap")
   }
   return(retval)
})


# TODO module for cluster plot?

selctedCluster <-
  callModule(clusterServer,
             "selected",
             projections,
             reactive(input$gene_id_sch))



# TODO module for heatmap?
output$selectedHeatmap <- renderImage({
  if (DEBUG)
    cat(file = stderr(), "output$selectedHeatmap\n")
  featureData = featureDataReact()
  gbm_matrix = gbm_matrix()
  projections = projections()
  genesin <- input$heatmap_geneids2
  sc = selctedCluster()
  scCL = sc$cluster
  # scBP = sc$brushedPs()
  scCells = sc$selectedCells()
 
  if(DEBUGSAVE)
    save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
  if (is.null(featureData) |
      is.null(gbm_matrix) |
      is.null(projections) | is.null(scCells) | length(scCells) == 0) {
    return(
      list(
        src = "empty.png",
        contentType = 'image/png',
        width = 96,
        height = 96,
        alt = "heatmap should be here"
      )
    )
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("selectedheatmap", id = "selectedHeatmap", duration = NULL)
  }
  
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
  
  # subsetData <-
  #   subset(projections, as.numeric(as.character(projections$dbCluster)) %in% scCL)
  # cells.1 <- rownames(brushedPoints(subsetData, scBP))
  cells.1 <- scCells
  retval = heatmapFunc(featureData, gbm_matrix, projections, genesin, cells = cells.1)
  
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "selectedHeatmap")
  }
  return(retval)
  
})

plotCoExpressionFunc <-
  function(featureData,
           gbm_log,
           upI,
           projections,
           genesin,
           cl3,
           dimx3,
           dimy3) {
    
    genesin = geneName2Index(genesin, featureData)
    # genesin <- toupper(genesin)
    # genesin <- gsub(" ", "", genesin, fixed = TRUE)
    # genesin <- strsplit(genesin, ',')
    
    subsetData <- subset(projections, dbCluster %in% cl3)
    # cells.1 <- rownames(subsetData)
    
    
    # map <-
    #   rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
    #if(DEBUG)cat(file=stderr(),map[1])
    
    expression <- gbm_log[genesin, ]
    #if(DEBUG)cat(file=stderr(),rownames(expression))
    
    #expression<-expression[complete.cases(expression),]
    #if(DEBUG)cat(file=stderr(),rownames(expression))
    
    # display genes not found
    # notFound = genesin[[1]][which(!genesin[[1]] %in% featureData$Associated.Gene.Name)]
    # if (length(notFound) > 0) {
    #   if (!is.null(getDefaultReactiveDomain())) {
    #     showNotification(
    #       paste("following genes were not found", notFound, collapse = " "),
    #       id = "plotCoExpressionNotFound",
    #       type = "warning",
    #       duration = 20
    #     )
    #   }
    # }
    
    validate(need(
      is.na(sum(expression)) != TRUE,
      'Gene symbol incorrect or genes not expressed'
    ))
    
    bin <- expression
    bin[] <- 0
    
    for (i in 1:nrow(expression)) {
      x <- Mclust(expression[i, ], G = 2)
      bin[i, ] <- x$classification
    }
    bin <- bin - 1
    allexprs <- apply(bin, 2, sum)
    plotexprs <- allexprs
    plotexprs[] <- 0
    plotexprs[allexprs >= length(rownames(bin))] <- 1
    # TODO positiveCells is changing too often. Maybe this can be controlled a bit more? check if changed? Use global variable not a reactive value?
    positiveCells$positiveCells <- allexprs >= length(rownames(bin))
    positiveCells$positiveCellsAll <- plotexprs
    
    mergeExprs <- plotexprs[rownames(subsetData)]
    #if(DEBUG)cat(file=stderr(),length(mergeExprs))
    
    subsetData$CoExpression <- factor(mergeExprs)
    subsetData$dbCluster <- as.factor(subsetData$dbCluster)
    p1 <-
      ggplot(subsetData,
             aes_string(x = dimx3, y = dimy3)) +
      geom_point(aes_string(
        shape = "sample",
        # alpha = 'CoExpression',
        color = "dbCluster"
      ),
      size = 4) +
      theme_bw()
    
    if (DEBUG)
      cat(file = stderr(), "output$plotCoExpression:done\n")
    return(p1)
    
  }
# TODO module?
output$plotCoExpression <- renderPlot({
  if (DEBUG)
    cat(file = stderr(), "output$plotCoExpression\n")
  # if (vvvvvvv$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  gbm_log = log2cpm()
  upI = updateInputx3() # no need to check because this is done in projections
  projections = projections()
  if (is.null(featureData) |
      is.null(projections) |
      is.null(log2cpm) | is.null(input$clusters3)) {
    return(NULL)
  }
  
  genesin <- input$mclustids
  cl3 = input$clusters3
  dimx3 = input$dimension_x3
  dimy3 = input$dimension_y3
  posCells = positiveCells$positiveCells # we use this variable to be able to save the global variable in this context
  posCellsAll = positiveCells$positiveCellsAll
  
  
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/plotCoExpression.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/plotCoExpression.RData")
  p1 = plotCoExpressionFunc(featureData,
                            gbm_log,
                            upI,
                            projections,
                            genesin,
                            cl3,
                            dimx3,
                            dimy3)
  return(p1)
})



# TODO module download
output$downloadExpressionOnOff <- downloadHandler(
  filename = function() {
    paste(input$clusters3, "PositiveCells.csv", sep = '_')
  },
  content = function(file) {
    featureData = featureDataReact()
    log2cpm = log2cpm()
    
    if(is.null(featureData) | is.null(log2cpm) | is.null(positiveCells$positiveCells)){
      return(NULL)
    }
    
    cells <- positiveCells$positiveCells
    #if(DEBUG)cat(file=stderr(),cells[1:5])
    
    if (length(cells) == 1) {
      subsetExpression <- log2cpm[, cells]
      subsetExpression <-
        as.data.frame(subsetExpression, row.names = rownames(log2cpm))
      colnames(subsetExpression) <- cells
      subsetExpression$Associated.Gene.Name <-
        featureData[rownames(subsetExpression), 'Associated.Gene.Name']
      write.csv(subsetExpression, file)
    }
    else{
      subsetExpression <- log2cpm[, cells]
      #cat(stderr(),colnames(subsetExpression)[1:5])
      subsetExpression$Associated.Gene.Name <-
        featureData[rownames(subsetExpression), 'Associated.Gene.Name']
      #cat(stderr(),colnames(subsetExpression))
      write.csv(subsetExpression, file)
    }
  }
)

# TODO do we need it?
output$onOffTable <- DT::renderDataTable({
  if (DEBUG)
    cat(file = stderr(), "output$onOffTable\n")
  projections = projections()
  posCellsAll = positiveCells$positiveCellsAll # we use this variable to be able to save the global variable in this context
  
  if (is.null(projections)) {
    return(NULL)
  }
  
  
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/onOffTable.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/onOffTable.RData")
  
  merge <- projections
  if (DEBUG)
    cat(file = stderr(),
        paste("positiveCells$positiveCellsAll:---", posCellsAll, "---\n"))

  merge$CoExpression <- posCellsAll
  df <-
    as.data.frame(table(merge[, c('dbCluster', 'CoExpression')]))
  dfOut <- cast(df, dbCluster ~ CoExpression)
  colnames(dfOut) <- c("Cluster", 'OFF', 'ON')
  rownames(dfOut) <- dfOut$Cluster
  dfOut['Sum', ] <- c('', sum(dfOut$OFF), sum(dfOut$ON))
  if (DEBUG)
    cat(file = stderr(), "output$onOffTable:done\n")
  DT::datatable(dfOut)
  
})

# TODO as module
# coexpression binarized
output$clusters3 <- renderUI({
  if (DEBUG)
    cat(file = stderr(), "output$clusters3\n")
  projections = projections()
  if (is.null(projections)) {
    HTML("Please load data first")
    return(NULL)
  }
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/clusters3.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/clusters3.RData")
  
  noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
  selectizeInput(
    "clusters3",
    label = "Cluster",
    choices = c(0:noOfClusters),
    selected = 0,
    multiple = TRUE
  )
  
})

geneGrp_vioFunc <- function(genesin, projections, gbm, featureData, minExpr=1, dbCluster) {
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ',')
  
  map <-
    rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
  if (DEBUG)
    cat(file = stderr(), length(map))
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmapWarning")
    removeNotification(id = "heatmapNotFound")
  }
  if (length(map) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no genes found",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(
      NULL
    )
  }
  
  expression <- colSums(as.matrix(exprs(gbm[map, ])) >= minExpr)
  
  
 
  projections <- cbind(projections, coExpVal = expression)
  # if(class(projections[,dbCluster])=="factor"){
  p1 <-
    ggplot(projections, aes_string(factor(projections[,dbCluster]), "coExpVal", fill = factor(projections[,dbCluster]))) +
    geom_violin(scale = "width") +
    stat_summary(
      fun.y = median,
      geom = "point",
      size = 5,
      color = 'black'
    ) +
    stat_summary(fun.data = n_fun, geom = "text") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 60,
        size = 12,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      legend.position = "none"
    ) +
    xlab(dbCluster) +
    ylab('number genes from list') 
  # }else{
  #   return(NULL)
  # }
  if (DEBUG)
    cat(file = stderr(), "output$gene_vio_plot:done\n")
  return(p1)
  
}

# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# TODO module for violin plot  ??
output$geneGrp_vio_plot <- renderPlot({
  if (DEBUG)
    cat(file = stderr(), "output$geneGrp_vio_plot\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  projections = projections()
  gbm = gbm()
  geneListStr = input$geneGrpVioIds
  projectionVar = input$dimension_xVioiGrp
  minExpr = input$coEminExpr
  upI = updateInputx3() # no need to check because this is done in projections
  if (is.null(projections)) {
    if (DEBUG)
      cat(file = stderr(), "output$geneGrp_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/geneGrp_vio_plot.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/geneGrp_vio_plot.RData")
  
  retVal = geneGrp_vioFunc(genesin = geneListStr,
                           projections = projections,
                           gbm = gbm,
                           featureData = featureData,
                           minExpr = minExpr,
                           dbCluster = projectionVar)
  if (DEBUG)
    cat(file = stderr(), "output$plotCoExpression:done\n")
  return(retVal)
  
})
