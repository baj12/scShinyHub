output$tsne_main <- renderPlotly({
  if(DEBUG)cat(file=stderr(), "output$tsne_main\n")
  tsne.data = tsne.data()
  if(is.null(tsne.data)){
    if(DEBUG)cat(file=stderr(), "output$tsne_main:NULL\n")
    return(NULL)
  }
  tsne.data <- as.data.frame(tsne.data)
  #cat(stderr(),colnames(tsne.data)[1:5])
  tsne.data$dbCluster <- as.factor(tsne.data$dbCluster)
  
  p <-
    plot_ly(
      tsne.data,
      x = ~ V1,
      y = ~ V2,
      z = ~ V3,
      type = "scatter3d",
      color =  ~ dbCluster,
      hoverinfo = "text",
      text = paste('Cluster:', tsne.data$dbCluster),
      mode = 'markers',
      marker =
        list(
          line = list(width = 0),
          size = rep(10, nrow(tsne.data)),
          sizeref = 3
        )
    )
  if(DEBUG)cat(file=stderr(), "output$tsne_main: done\n")
  return(layout(p))
  
  
})


output$plotUmiHist <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output_plotUmiHist\n")
  gbm = gbm()
  if(is.null(gbm))
    return(NULL)
  hist(colSums(as.matrix(exprs(gbm))), breaks = 50, main="histogram of number of UMI per cell")
})

# TODO move to generalQC
output$variancePCA <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$variancePCA\n")
  h2("hello")
  pca = pca()
  if(is.null(pca))
    return(NULL)
  barplot(pca$var_pcs, main="Variance captured by first PCs")
})


