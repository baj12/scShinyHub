output$tsne_main <- renderPlotly({
  if(DEBUG)cat(file=stderr(), "output$tsne_main\n")
  projections = projections()
  if(is.null(projections)){
    if(DEBUG)cat(file=stderr(), "output$tsne_main:NULL\n")
    return(NULL)
  }
  projections <- as.data.frame(projections)
  #cat(stderr(),colnames(projections)[1:5])
  projections$dbCluster <- as.factor(projections$dbCluster)
  
  p <-
    plot_ly(
      projections,
      x = ~ tsne1,
      y = ~ tsne2,
      z = ~ tsne3,
      type = "scatter3d",
      color =  ~ dbCluster,
      hoverinfo = "text",
      text = paste('Cluster:', as.numeric(as.character(projections$dbCluster))),
      mode = 'markers',
      marker =
        list(
          line = list(width = 0),
          size = rep(10, nrow(projections)),
          sizeref = 3
        )
    )
  if(DEBUG)cat(file=stderr(), "output$tsne_main: done\n")
  return(layout(p))
  
  
})
source("moduleServer.R", local=TRUE)
source("reactives.R", local=TRUE)


r<-callModule(tableSelectionServer, "cellSelectionTSNEMod", inputTSNESample)


output$plotUmiHist <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output_plotUmiHist\n")
  gbm = gbm()
  if(is.null(gbm))
    return(NULL)
  hist(colSums(as.matrix(exprs(gbm))), breaks = 50, main = "histogram of number of UMIs per cell")
})

output$plotSampleHist <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output_sampleHist\n")
  sampleInf = sampleInfo()
  if(is.null(sampleInf))
    return(NULL)
  if(DEBUGSAVE) save(file = "~/scShinyHubDebug/sampleHist.RData", list = ls())
  # load(file = "~/scShinyHubDebug/sampleHist.RData")
  sampleHistFunc(sampleInf)
})

output$variancePCA <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$variancePCA\n")
  h2("hello")
  pca = pca()
  if(is.null(pca))
    return(NULL)
  barplot(pca$var_pcs, main="Variance captured by first PCs")
})


