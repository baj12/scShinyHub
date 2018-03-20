dge <- reactive({
  if(DEBUG)cat(file=stderr(), "dge\n")
  featureData = featureDataReact()
  log2cpm = log2cpm()
  tsne.data = tsne.data()
  if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
    return(NULL)
  }
  
  subsetData <- subset(tsne.data, dbCluster %in% input$clusters1)
  cells.1 <- rownames(brushedPoints(subsetData, input$db1))
  cells.2 <- rownames(brushedPoints(subsetData, input$db2))
  
  subsetExpression <- log2cpm[, union(cells.1, cells.2)]
  
  genes.use <- rownames(subsetExpression)
  data.1 = apply(subsetExpression[genes.use, cells.1], 1, expMean)
  data.2 = apply(subsetExpression[genes.use, cells.2], 1, expMean)
  total.diff = (data.1 - data.2)
  
  genes.diff = names(which(abs(total.diff) > .2))
  genes.use = ainb(genes.use, genes.diff)
  
  toReturn <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  toReturn[, "avg_diff"] = total.diff[rownames(toReturn)]
  toReturn$Associated.Gene.Name <-
    featureData[rownames(toReturn), 'Associated.Gene.Name']
  selectedDge <- toReturn
  cat(stderr(), rownames(toReturn)[1:5])
  return(toReturn)
  
  # })
})
