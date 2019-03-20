# The output type has to be in line with the tablist item. I.e. plotOutput in this case
output$Dummy_plot <- renderPlot({
  projections <- projections()
  dummyNRow <- DummyReactive()
  if (is.null(projections) | is.null(dummyNRow)) {
    return(NULL)
  }
  
  if (DEBUG) cat(file = stderr(), paste("Dummy_plot:\n"))
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/Dummy_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/Dummy_plot.RData")
  
  plot(dummyNRow)
})


output$DummySavedPlot <- renderImage({
  if (DEBUG) cat(file = stderr(), paste("DummySavedPlot:\n"))
  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) cat(file = stderr(), paste("DummySavedPlot:NULL\n"))
    return(NULL)
  }
  retVal <- imageDummyPrecompute()
  if (DEBUG) cat(file = stderr(), paste("DummySavedPlot:done\n"))
  return(retVal)
})
