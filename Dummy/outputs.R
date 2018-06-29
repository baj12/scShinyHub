# The output type has to be in line with the tablist item. I.e. plotOutput in this case
output$Dummy_plot <- renderPlot({
  tsne.data = tsne.data()
  dummyNRow = DummyReactive()
  if( is.null(tsne.data) | is.null(dummyNRow) ){
    return(NULL)
  }
  
  if(DEBUG)cat(file=stderr(), paste("Dummy_plot:\n"))
  if(DEBUGSAVE) save(file="~/scShinyHubDebug/Dummy_plot.RData", list=ls())
  # load(file="~/scShinyHubDebug/dge_plot1.RData")
  
  plot(dummyNRow)
  
})


output$DummySavedPlot <- renderImage({
  if(DEBUG)cat(file=stderr(), paste("DummySavedPlot:\n"))
  gbm = gbm()
  if (is.null(gbm)){
    if(DEBUG)cat(file=stderr(), paste("DummySavedPlot:NULL\n"))
    return(NULL)
  }
  retVal = imageDummyPrecompute()
  if(DEBUG)cat(file=stderr(), paste("DummySavedPlot:done\n"))
  return(retVal)
})