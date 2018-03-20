# TODO mnodule for heatmap?  
output$heatmap <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$heatmap\n")
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
    return(NULL)
  }
  
  genesin <- input$heatmap_geneids
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ',')
  
  map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
  cat(file = stderr(), length(map))
  
  expression <- log2cpm[map, ]
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or genes not expressed'
  ))
  
  tsne.data <- tsne.data[order(tsne.data$dbCluster), ]
  
  expression <- expression[, rownames(tsne.data)]
  expression <- expression[complete.cases(expression), ]
  
  annotation <- data.frame(factor(tsne.data$dbCluster))
  rownames(annotation) <- colnames(expression)
  colnames(annotation) <- c('Cluster')
  
  h <-
    pheatmap(
      as.matrix(expression),
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
      breaks = seq(-6, 6, by = .12),
      colorRampPalette(rev(brewer.pal(
        n = 6, name =
          "RdBu"
      )))(100)
      
    )
  h
  
  # })
})


# TODO module for cluster plot?  

selctedCluster <- callModule(clusterServer, "selected", tsne.data, reactive(input$gene_id_sch))

# output$clusterPlot2 <- renderPlot({
#   if(DEBUG)cat(file=stderr(), "output$clusterPlot2\n")
#   featureData = featureDataReact()
#   log2cpm = log2cpm()
#   tsne.data = tsne.data()
#   if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
#     return(NULL)
#   }
#   
#   # isolate({
#   geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
#                                          toupper(input$gene_id_sch)), ])[1]
#   
#   expression <- log2cpm[geneid, ]
#   
#   validate(need(is.na(sum(expression)) != TRUE, ''))
#   
#   tsne.data <- cbind(tsne.data, t(expression))
#   names(tsne.data)[names(tsne.data) == geneid] <- 'values'
#   
#   if(DEBUG)cat(file=stderr(), paste("output$dge_plot1:---",input$clusters2,"---\n"))
#   subsetData <- subset(tsne.data, dbCluster %in% input$clusters2)
#   p1 <-
#     ggplot(subsetData,
#            aes_string(x = input$dimension_x2, y = input$dimension_y2)) +
#     geom_point(aes_string(size = 2, color = 'values')) +
#     geom_point(shape = 1,
#                size = 4,
#                aes(colour = dbCluster)) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(
#         angle = 90,
#         size = 12,
#         vjust = 0.5
#       ),
#       axis.text.y = element_text(size = 10),
#       strip.text.x = element_text(size = 16),
#       strip.text.y = element_text(size = 14),
#       axis.title.x = element_text(face = "bold", size = 16),
#       axis.title.y = element_text(face = "bold", size = 16),
#       legend.position = "none"
#     ) +
#     ggtitle(paste(toupper(input$gene_id_sch), input$clusters2, sep =
#                     '-Cluster')) +
#     scale_colour_gradient2(low = 'grey50', high = "red")
#   p1
#   # })
# })


# TODO module for heatmap?  
output$selectedHeatmap <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$selectedHeatmap\n")
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
    return(NULL)
  }
  
  genesin <- input$heatmap_geneids2
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ',')
  
  sc =selctedCluster()
  
  subsetData <-
    subset(tsne.data, tsne.data$dbCluster %in% sc$cluster)
  cells.1 <- rownames(brushedPoints(subsetData, sc$brushedPs()))

  map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
  #if(DEBUG)cat(file=stderr(),map[1])
  
  expression <- log2cpm[map, cells.1]
  if(DEBUG)cat(file = stderr(), rownames(expression))
  
  expression <- expression[complete.cases(expression), ]
  if(DEBUG)cat(file = stderr(), rownames(expression))
  mColor <- max(expression)
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or genes not expressed'
  ))
  

  h <-
    pheatmap(
      as.matrix(expression),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = 'row',
      fontsize_row = 10,
      labels_col = colnames(expression),
      labels_row = featureData[rownames(expression), 'Associated.Gene.Name'],
      show_rownames = TRUE,
      show_colnames = FALSE,
      breaks = seq(-6, 6, by = .12),
      colorRampPalette(rev(brewer.pal(
        n = 6, name =
          "RdBu"
      )))(100)
      
    )
  h
})



# TODO module?  
output$plotCoExpression <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$plotCoExpression\n")
  # if (vvvvvvv$doPlot == FALSE)
  #   return()
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(tsne.data) | is.null(log2cpm) | is.null(input$clusters3)){
    return(NULL)
  }
  
  # isolate({
  genesin <- input$mclustids
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ',')
  
  subsetData <-
    subset(tsne.data, dbCluster %in% input$clusters3)
  cells.1 <- rownames(subsetData)
  
  
  map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
  #if(DEBUG)cat(file=stderr(),map[1])
  
  expression <- log2cpm[map, ]
  #if(DEBUG)cat(file=stderr(),rownames(expression))
  
  #expression<-expression[complete.cases(expression),]
  #if(DEBUG)cat(file=stderr(),rownames(expression))
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    'Gene symbol incorrect or genes not expressed'
  ))
  
  bin <- expression
  bin[] <- 0
  
  for (i in 1:nrow(expression))
  {
    x <- Mclust(expression[i, ], G = 2)
    bin[i, ] <- x$classification
  }
  bin <- bin - 1
  allexprs <- apply(bin, 2, sum)
  plotexprs <- allexprs
  plotexprs[] <- 0
  plotexprs[allexprs >= length(rownames(bin))] <- 1
  positiveCells$positiveCells <- allexprs >= length(rownames(bin))
  positiveCells$positiveCellsAll <- plotexprs
  #save(subsetData,bin,allexprs,file='~/Desktop/test.Rds')
  #if(DEBUG)cat(file=stderr(),names(allexprs))
  
  mergeExprs <- plotexprs[rownames(subsetData)]
  #if(DEBUG)cat(file=stderr(),length(mergeExprs))
  
  subsetData$CoExpression <- mergeExprs
  #if(DEBUG)cat(file=stderr(),colnames(subsetData))
  
  p1 <-
    ggplot(subsetData,
           aes_string(x = input$dimension_x3, y = input$dimension_y3)) +
    geom_point(aes_string(size = 2, color = 'CoExpression')) +
    geom_point(shape = 1,
               size = 4,
               aes(colour = dbCluster)) +
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
    #ggtitle(paste(toupper(input$gene_id),input$cluster,sep='-Cluster'))+
    scale_colour_gradient2(low = 'grey50', high = "red")
  p1
  # })
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
  if(DEBUG)cat(file=stderr(), "output$onOffTable\n")
  tsne.data = tsne.data()
  
  if( is.null(tsne.data | is.null(positiveCells$positiveCellsAll)) ){
    return(NULL)
  }
  
  merge <- tsne.data
  if(DEBUG)cat(file=stderr(), paste("positiveCells$positiveCellsAll:---",positiveCells$positiveCellsAll,"---\n"))
  
  merge$CoExpression <- positiveCells$positiveCellsAll
  df <-
    as.data.frame(table(merge[, c('dbCluster', 'CoExpression')]))
  dfOut <- cast(df, dbCluster ~ CoExpression)
  colnames(dfOut) <- c("Cluster", 'OFF', 'ON')
  rownames(dfOut) <- dfOut$Cluster
  dfOut['Sum', ] <- c('', sum(dfOut$OFF), sum(dfOut$ON))
  DT::datatable(dfOut)
  
})

