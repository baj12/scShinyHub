source("reactives.R")

# myHeavyCalculations = list(c("scaterPNG", "scaterPNG"))

# Expression ------------------------------------------------------------------
expCluster <- callModule(clusterServer, "expclusters", projections, reactive(input$gene_id))

# these observes should be independant of projections since this will be then executed by default for any changes
updateInputx4 = reactive({
  tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "dimension_x4",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[1]
  )
  
  # Can also set the label and select items
  updateSelectInput(session, "dimension_y4",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[2]
  )
  return(TRUE)
})

output$NumberOfGenesInclude = renderText({
  idx = scGeneIdxInclude()
  paste("Number of genes to be included: ", length(idx))
})

output$NumberOfGenesExclude = renderText({
  idx = scGeneIdxExclude()
  paste("Number of genes to be included: ", length(idx))
})


# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# TODO module for violin plot  ??
output$gene_vio_plot <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$gene_vio_plot\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  # log2cpm = log2cpm()
  gbm_log = gbm_log()
  projections = projections()
  g_id = input$gene_id
  
  if ( is.null(featureData) | is.null(gbm_log) | is.null(projections)) {
    if ( DEBUG ) cat(file = stderr(), "output$gene_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/gene_vio_plot.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/gene_vio_plot.RData")
  
  geneid = geneName2Index(g_id, featureData)  
  
  # geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
  #                                        toupper(input$gene_id)), ])[1]
  
  # expression <- exprs(gbm_log)[geneid, ]
  if (length(geneid) == 1) {
    expression = exprs(gbm_log)[geneid, ]
  }else{
    expression = Matrix::colSums(exprs(gbm_log)[geneid, ])
  }
  
  validate(need(is.na(sum(expression)) != TRUE, ''))
  
  projections <- cbind(projections, expression)
  names(projections)[length(projections)] <- 'values'
  
  p1 <-
    ggplot(projections, aes(factor(dbCluster), values, fill = factor(dbCluster))) +
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
        angle = 90,
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
    xlab('Cluster') +
    ylab('Expression') +
    ggtitle(paste(toupper(featureData[geneid,"Associated.Gene.Name"]),collapse = ", "))
  if(DEBUG)cat(file=stderr(), "output$gene_vio_plot:done\n")
  return(p1)
  # })
})

# EXPLORE TABL DOWNLOAD SELECTED WITH BRUSH ------------------------------------------------------------------
# TODO move to were it belongs  
# TODO module for download?  
# TODO either integrate in module or use input$bi as selected cell names as a parameter/reactive/global variable.
#      this function will not work as expected as of now
output$downloadExpression <- downloadHandler(
  filename = function() {
    paste(input$cluster, "Selected_Expression_table.csv", sep = '_')
  },
  content = function(file) {
    featureData = featureDataReact()
    # log2cpm = log2cpm()
    gbm_log = gbm_log()
    projections = projections()
    if(is.null(featureData) | is.null(gbm_log) | is.null(projections)){
      return(NULL)
    }
    if (DEBUGSAVE) 
      save(file = "~/scShinyHubDebug/downloadExpression.RData", list = c(ls(),ls(envir = globalenv())))
    # load(file="~/scShinyHubDebug/downloadExpression.RData")
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                           toupper(input$gene_id)), ])[1]
    
    expression <- exprs(gbm_log)[geneid, ]
    #cat(stderr(),colnames(expression)[1:5])
    projections <- cbind(projections, t(expression))
    #if(DEBUG)cat(file=stderr(),grep('^T_',rownames(projections)))
    
    names(projections)[names(projections) == geneid] <- 'values'
    
    #if(DEBUG)cat(file=stderr(),grep('^T_',rownames(projections)))
    
    subsetData <- subset(projections, dbCluster == input$cluster)
    #if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
    cells.names <- brushedPoints(subsetData, input$b1, allRows = T)
    #if(DEBUG)cat(file=stderr(),colnames(cells.names))
    cells <-
      rownames(subsetData[which(cells.names$selected_ == TRUE), ])
    #if(DEBUG)cat(file=stderr(),cells[1:5])
    
    if (length(cells) == 1) {
      subsetExpression <- exprs(gbm_log)[, cells]
      subsetExpression <-
        as.data.frame(subsetExpression, row.names = rownames(gbm_log))
      colnames(subsetExpression) <- cells
      subsetExpression$Associated.Gene.Name <-
        featureData[rownames(subsetExpression), 'Associated.Gene.Name']
      write.csv(subsetExpression, file)
    }
    else{
      subsetExpression <- exprs(gbm_log)[, cells]
      #cat(stderr(),colnames(subsetExpression)[1:5])
      
      subsetExpression$Associated.Gene.Name <-
        featureData[rownames(subsetExpression), 'Associated.Gene.Name']
      #cat(stderr(),colnames(subsetExpression))
      write.csv(subsetExpression, file)
    }
  }
)

##############################
### Panel Plot
# TODO as module
# data expression panel plot 
output$clusters4 <- renderUI({
  if(DEBUG)cat(file=stderr(), "output$clusters4\n")
  projections = projections()
  upI = updateInputx4()
  if(is.null(projections)){
    HTML("Please load data firts")
  }else{
    noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectInput(
      "clusters4",
      label = "Cluster",
      choices = c(c('All'),c(0:noOfClusters)),
      selected = 0
    )
  }
})  

output$cvHist <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$cvHist\n")
  
  gbm_log = gbm_log()
  if(is.null(gbm_log)){
    return(NULL)
  }
  
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/cvHist.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/cvHist.RData")
  
  retVal = apply(exprs(gbm_log),1,FUN = function(x){mean(x)/sd(x)})
  return(hist(retVal, breaks = 100))
  
})


cvHeatmapFunc <- function(featureData, gbm_matrix, projections, genesin, cells){
  genesin = geneName2Index(genesin, featureData)
  
  expression <- gbm_matrix[genesin, cells]
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or genes not expressed'
  ))
  
  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]
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
  nonZeroRows = which(rowSums(expression)>0)
  pheatmap(
    as.matrix(expression)[nonZeroRows,order(annotation[,1], annotation[,2])],
    cluster_rows = FALSE,
    cluster_cols = TRUE,
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
output$cvHeatMap <- renderImage({
  if(DEBUG)cat(file=stderr(), "output$heatmap\n")
  featureData = featureDataReact()
  gbm_log = gbm_log()
  projections = projections()
  genesin <- input$cvHeatmap_geneids
  
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
    save(file = "~/scShinyHubDebug/cvHeatMap.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file = "~/scShinyHubDebug/cvHeatMap.RData")
  gbm_matrix = as.matrix(exprs(gbm_log))
  
  if (nchar(genesin)>0) {
    ensNames = geneName2Index(genesin, featureData)
    gbm_matrix = gbm_matrix[ensNames,]
  }
  retVal = apply(gbm_matrix,1,FUN = function(x){mean(x)/sd(x)})
  retValOrder = order(retVal, decreasing = FALSE)
  if( (nchar(genesin)<1) & (length(retValOrder)>100)){
    retValOrder = retValOrder[1:100]
    genesin = 
  }
  retval = cvHeatmapFunc(featureData = featureData[retValOrder, ], 
                         gbm_matrix = gbm_matrix[retValOrder,], 
                         projections = projections[retValOrder,], 
                         genesin = genesin, 
                         cells = colnames(gbm_matrix))
  
  if(!is.null(getDefaultReactiveDomain())){
    removeNotification( id="heatmap")
  }
  return(retval)
})



output$panelPlot <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$panelPlot\n")
  
  featureData = featureDataReact()
  # log2cpm = log2cpm()
  gbm_log = gbm_log()
  projections = projections()
  if(is.null(featureData) | is.null(gbm_log) | is.null(projections)){
    return(NULL)
  }
  
  genesin <- input$panelplotids
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ',')
  genesin<-genesin[[1]]
  cl4 = input$clusters4
  dimx4 = input$dimension_x4
  dimy4 = input$dimension_y4
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/panelPlot.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/panelPlot.RData")
  
  if(DEBUG)cat(file=stderr(),length(genesin))
  par(mfrow=c(ceiling(length(genesin)/4),4), mai = c(0, 0., 0., 0.))
  rbPal <- colorRampPalette(c('#f0f0f0','red'))
  if(DEBUG)cat(file=stderr(),cl4)
  
  if (cl4 == 'All') 
  {
    for (i in 1:length(genesin)){
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              exprs(gbm_log)[
                rownames(featureData[which(featureData$Associated.Gene.Name==genesin[i]),])
                ,]
            ),breaks = 10))]
      plot(projections[, dimx4],projections[, dimy4],col=Col,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
      title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
      if(DEBUG)cat(file=stderr(),genesin[i])
    }
  }
  else{
    for (i in 1:length(genesin)){
      
      subsetTSNE <- subset(projections, dbCluster == cl4)
      
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              exprs(gbm_log)[
                rownames(featureData[which(featureData$Associated.Gene.Name==genesin[i]),])
                ,]
            ),breaks = 10))]
      
      names(Col)<-rownames(projections)
      plotCol<-Col[rownames(subsetTSNE)]
      plot(subsetTSNE[, dimx4],subsetTSNE[, dimy4],col=plotCol,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
      title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
      if(DEBUG)cat(file=stderr(), cl4)
    }
  }
})


##############################
### Scater QC


output$scaterQC <- renderImage({
  if(DEBUG)cat(file=stderr(), "output$scaterQC\n")
  scaterReads = scaterReads()
  if(is.null(scaterReads)){
    return(NULL)
  }
  
  scaterPNG()
}
)

output$tsne_plt <- renderPlotly({
  if(DEBUG)cat(file=stderr(), "output$tsne_plt\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  # log2cpm = log2cpm()
  gbm_log = gbm_log()
  g_id = input$gene_id
  projections = projections()
  
  if (is.null(featureData) | is.null(gbm_log) | is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/tsne_plt.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/tsne_plt.RData")
  
  
  geneid = geneName2Index(g_id, featureData)  
  
  if (length(geneid) == 1) {
    expression = exprs(gbm_log)[geneid, ]
  }else{
    expression = Matrix::colSums(exprs(gbm_log)[geneid, ])
  }
  
  # expression <- log2cpm[geneid, ]
  # cat(file = stderr(), rownames(expression))
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or gene not expressed'
  ))
  
  projections <- cbind(projections, expression)
  names(projections)[ncol(projections)] <- 'values'
  
  p <-
    plot_ly(
      projections,
      x = ~ tsne1,
      y = ~ tsne2,
      z = ~ tsne3,
      type = "scatter3d",
      hoverinfo = "text",
      text = paste('Cluster:', as.numeric(as.character(projections$dbCluster))),
      mode = 'markers',
      marker = list(
        size = 2,
        line = list(width = 0),
        color =  ~ values,
        colors = 'Greens'
      )
    )
  layout(p, title = paste(toupper(featureData[geneid,"Associated.Gene.Name"]),collapse = ", "))
  # })
})








