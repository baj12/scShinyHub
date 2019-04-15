#' DummyFunc ----
#' function doing the actual work
#' being used below for the shiny widget and in the report.
#' in the report all reactives that have been used before will be accessible by name
DummyFunc <- function(scEx_log) {
  nrow(scEx_log)
}

# DummyReactive ----
#' DummyReactive
#' controls the calculation for the shiny widget
DummyReactive <- reactive({
  # some debugging messages saying that we start the process
  if (DEBUG) cat(file = stderr(), "DummyReactive started.\n")
  # track how much time is spent here
  start.time <- base::Sys.time()
  
  # remove any notification on exit that we don't want
  on.exit({
    printTimeEnd(start.time, "DummyReactive")
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DummyFunc")
  })
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("loading", id = "DummyFunc", duration = NULL)
  }
  # remove any permanant notification if we rerun reactive
  if (!is.null(getDefaultReactiveDomain()))
    removeNotification(id = "DummyFuncPerm")
  
  
  # call dependancies (reactives)
  scEx <- scEx()                 # raw data, filtered by genes/cells
  scEx_log <- scEx_log()         # normalized data
  prj <- projections()           # projections, includes manually set groups
  inputData <- inputData()       # raw, unfiltered data
  pca <- pca()                   # pca projections, loadings...
  tsne <- tsne()                 # tSNE projections
  dbCluster <- dbCluster()       # clustering results, also available through projections
  scaterReads <- scaterReads()   #singleCellExperiment object/reactive with QC metrix
  # check if they are available
  if (is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "pca:NULL\n")
    return(NULL)
  }
  
  # for development and debugging purposes
  # this is run after loading all reactive values
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DummyReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/DummyReactive.RData')
  
  # actual calculation
  # we separate the actual calculations from anything related to shiny to be able to use it in the reports
  retVal <- DummyFunc(scEx_log)
  
  if (retVal == 0 & !is.null(getDefaultReactiveDomain())) {
    showNotification("Dummy is 0", type = "warning", duration = NULL, id = "DummyFuncPerm") # has to be removed by use, no removeNotification is following.
    return(NULL)
  }
  
  # print debugging information on the console
  printTimeEnd(start.time, "DummyReactive")
  # for automated shiny testing using shinytest
  exportTestValues(DummyReactive = {retVal})  
  # I prefer calling return to return a value, this way I know what is happeing... (Am I too old?)
  return(retVal)
})

# imageDummyPrecompute ----
#' imageDummyPrecompute
#' this is an example of calculating an image to a file and returning the local reference
imageDummyPrecompute <- reactive({
  # track how much time is spent here
  start.time <- base::Sys.time()
  
  # remove any notification on exit that we don't want
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "imageDummyPrecompute")
  )
  # remove any permanant notification if we rerun reactive
  if (!is.null(getDefaultReactiveDomain()))
    removeNotification(id = "imageDummyPrecomputePerm")
  
  # some debugging messages
  if (DEBUG) cat(file = stderr(), "imageDummyPrecompute started.\n")
  
  # call dependancies (reactives)
  # pick the ones that are needed and remove others
  scEx <- scEx()                 # raw data, filtered by genes/cells
  scEx_log <- scEx_log()         # normalized data
  prj <- projections()           # projections, includes manually set groups
  inputData <- inputData()       # raw, unfiltered data
  pca <- pca()                   # pca projections, loadings...
  dbCluster <- dbCluster()       # clustering results, also available through projections
  
  # check if they are available
  if (is.null(scEx)) {
    if (DEBUG) cat(file = stderr(), "pca:NULL\n")
    return(NULL)
  }
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("preparing", id = "imageDummyPrecompute", duration = NULL)
  }
  
  # for development and debugging purposes
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DummyReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/DummyReactive.RData')
  
  ### actual work is done starting here
  width  <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  # For high-res displays, this will be greater than 1
  pixelratio <- session$clientData$pixelratio
  if (is.null(pixelratio)) pixelratio <- 1
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  myPNGwidth <- width / 96
  myPNGheight <- height / 96
  
  outfile <- paste0(tempdir(), "/dummy.png")
  if (DEBUG) cat(file = stderr(), paste("output file: ", outfile, "\n"))
  if (DEBUG) cat(file = stderr(), paste("output file normalized: ", normalizePath(outfile, mustWork = FALSE), "\n"))
  m <- data.frame("V1" = Matrix::colSums(assays(scEx)[["counts"]]))
  p <- ggplot(m, aes(V1)) + geom_bar()
  ggsave(file = normalizePath(outfile, mustWork = FALSE), plot = p, width = myPNGwidth, height = myPNGheight, units = "in")
  
  retVal <-  list(
    src = normalizePath(outfile, mustWork = FALSE),
    contentType = "image/png",
    width = width,
    height = height,
    alt = "Dummy should be here"
  )
  
  ### actual work is done
  
  printTimeEnd(start.time, "inputData")
  exportTestValues(imageDummyPrecompute = {retVal})  
  return(retVal)
})
