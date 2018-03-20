# Expression ------------------------------------------------------------------
expCluster <- callModule(clusterServer, "expclusters", tsne.data, reactive(input$gene_id))




# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# TODO module for violin plot  ??
output$gene_vio_plot <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$gene_vio_plot\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
    if(DEBUG)cat(file=stderr(), "output$gene_vio_plot:NULL\n")
    return(NULL)
  }
  
  geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                         toupper(input$gene_id)), ])[1]
  
  expression <- log2cpm[geneid, ]
  
  validate(need(is.na(sum(expression)) != TRUE, ''))
  
  tsne.data <- cbind(tsne.data, t(expression))
  names(tsne.data)[names(tsne.data) == geneid] <- 'values'
  
  p1 <-
    ggplot(tsne.data, aes(factor(dbCluster), values, fill = factor(dbCluster))) +
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
    ggtitle(toupper(input$gene_id))
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
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                           toupper(input$gene_id)), ])[1]
    
    expression <- log2cpm[geneid, ]
    #cat(stderr(),colnames(expression)[1:5])
    tsne.data <- cbind(tsne.data, t(expression))
    #if(DEBUG)cat(file=stderr(),grep('^T_',rownames(tsne.data)))
    
    names(tsne.data)[names(tsne.data) == geneid] <- 'values'
    
    #if(DEBUG)cat(file=stderr(),grep('^T_',rownames(tsne.data)))
    
    subsetData <- subset(tsne.data, dbCluster == input$cluster)
    #if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
    cells.names <- brushedPoints(subsetData, input$b1, allRows = T)
    #if(DEBUG)cat(file=stderr(),colnames(cells.names))
    cells <-
      rownames(subsetData[which(cells.names$selected_ == TRUE), ])
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

##############################
### Panel Plot

output$panelPlot <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$panelPlot\n")
  
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
    return(NULL)
  }
  
  genesin <- input$panelplotids
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ',')
  genesin<-genesin[[1]]
  
  if(DEBUG)cat(file=stderr(),length(genesin))
  par(mfrow=c(ceiling(length(genesin)/4),4), mai = c(0, 0., 0., 0.))
  rbPal <- colorRampPalette(c('#f0f0f0','red'))
  if(DEBUG)cat(file=stderr(),input$clusters4)
  
  if (input$clusters4 == 'All') 
  {
    for (i in 1:length(genesin)){
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              log2cpm[
                rownames(featureData[which(featureData$Associated.Gene.Name==genesin[i]),])
                ,]
            ),breaks = 10))]
      plot(tsne.data[,input$dimension_x4],tsne.data[,input$dimension_y4],col=Col,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
      title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
      if(DEBUG)cat(file=stderr(),genesin[i])
    }
  }
  else{
    for (i in 1:length(genesin)){
      
      subsetTSNE <- subset(tsne.data, dbCluster == input$clusters4)
      
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              log2cpm[
                rownames(featureData[which(featureData$Associated.Gene.Name==genesin[i]),])
                ,]
            ),breaks = 10))]
      
      names(Col)<-rownames(tsne.data)
      plotCol<-Col[rownames(subsetTSNE)]
      plot(subsetTSNE[,input$dimension_x4],subsetTSNE[,input$dimension_y4],col=plotCol,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
      title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
      if(DEBUG)cat(file=stderr(),input$clusters4)
    }
  }
})


##############################
### Scater QC

output$scaterQC <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$scaterQC\n")
  scaterReads = scaterReads()
  if(is.null(scaterReads)){
    return(NULL)
  }
  withProgress(message = 'calculating scater plot', value = 0, {
    p1 <- scater::plotQC(scaterReads, type = "highest-expression", col_by_variable="fixed")
  })
  p1
  #if(DEBUG)cat(file=stderr(), "DiffExpTest\n")
}
)

output$tsne_plt <- renderPlotly({
  if(DEBUG)cat(file=stderr(), "output$tsne_plt\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
    return(NULL)
  }
  geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                         toupper(input$gene_id)), ])[1]
  
  expression <- log2cpm[geneid, ]
  cat(file = stderr(), rownames(expression))
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or gene not expressed'
  ))
  
  tsne.data <- cbind(tsne.data, t(expression))
  names(tsne.data)[names(tsne.data) == geneid] <- 'values'
  
  p <-
    plot_ly(
      tsne.data,
      x = ~ V1,
      y = ~ V2,
      z = ~ V3,
      type = "scatter3d",
      hoverinfo = "text",
      text = paste('Cluster:', tsne.data$dbCluster),
      mode = 'markers',
      marker = list(
        size = 2,
        line = list(width = 0),
        color =  ~ values,
        colors = 'Greens'
      )
    )
  layout(p, title = toupper(input$gene_id))
  # })
})








