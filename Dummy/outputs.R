# list of files to be included in the report.
myZippedReportFiles <- c("dummyTableOutput.csv")


# The output type has to be in line with the tablist item. I.e. plotOutput in this case
#' here, we prepare the data in the DummyReactive and just provide the plot 
output$Dummy_plot <- renderPlot({
  # load reactive data
  dummyNRow <- DummyReactive()
  # return if nothing to be computed
  if (is.null(projections) | is.null(dummyNRow)) {
    return(NULL)
  }
  
  # some debugging messages
  if (DEBUG) cat(file = stderr(), paste("Dummy_plot:\n"))
  # for development and debugging purposes
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/Dummy_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/Dummy_plot.RData")
  
  # create and return the plot
  plot(dummyNRow)

})

observe({
  if (DEBUG) cat(file = stderr(), paste0("observe: clusters\n"))
  mod_cl1 <<- input$clusters
})

output$DummySavedPlot <- renderImage({
  # track how much time is spent here
  start.time <- base::Sys.time()
  
  # some debugging messages
  if (DEBUG) cat(file = stderr(), paste("DummySavedPlot:\n"))

  retVal <- imageDummyPrecompute()
  # return if nothing to be computed
  if (is.null(retVal)) {
    if (DEBUG) cat(file = stderr(), paste("DummySavedPlot:NULL\n"))
    return(NULL)
  }
  # print debugging information on the console
  printTimeEnd(start.time, "inputData")
  # for automated shiny testing using shinytest
  exportTestValues(DummySavedPlot = {retVal })  
  # I prefer calling return to return a value, this way I know what is happeing... (Am I too old?)
  return(retVal)
})
