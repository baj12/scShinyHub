
DummyFunc <- function(scEx_log) {
  # here we perform the calculations and provide the resulting data.
  # run_pca comes from the cellranger package.
  nrow(scEx_log)
}

# here we define reactive values/variables
# e.g.
DummyReactive <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DummyFunc")
  )
  
  # some debugging messages
  if (DEBUG) cat(file = stderr(), "pca\n")
  # call dependancies (reactives)
  scEx_log <- scEx_log()
  prj <- projections()
  scEx <- scEx()
  
  # check if they are available
  if (is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "pca:NULL\n")
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("loading", id = "DummyFunc", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DummyReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/DummyReactive.RData')
  
  # actual calculation
  retVal <- DummyFunc(scEx_log)
  
  if (retVal == 0 & !is.null(getDefaultReactiveDomain())) {
    showNotification("Dummy is 0", type = "warning", duration = NULL) # has to be removed by use, no removeNotification is following.
    return(NULL)
  }
  if (DEBUG) cat(file = stderr(), "inputData: done\n")
  return(retVal)
})

# declare heavy calculations
myHeavyCalculations <- list(c("running DummyReactive", "DummyReactive"))


# this will never be executed as it won't be called...

imageDummyPrecompute <- reactive({
  if (DEBUG) cat(file = stderr(), "imageDummyPrecompute\n")
  
  scEx <- scEx()
  # check if they are available
  if (is.null(scEx)) {
    if (DEBUG) cat(file = stderr(), "imageDummyPrecompute:NULL\n")
    return(NULL)
  }
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  
  
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("loading", id = "DummyFunc", duration = NULL)
  }
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
  m <- data.frame("V1" = Matrix::colSums(assays(scEx)[[1]]))
  p <- ggplot(m, aes(V1)) + geom_bar()
  ggsave(file = normalizePath(outfile, mustWork = FALSE), plot = p, width = myPNGwidth, height = myPNGheight, units = "in")
  
  if (DEBUG) cat(file = stderr(), "imageDummyPrecompute:done\n")
  
  return(list(
    src = normalizePath(outfile, mustWork = FALSE),
    contentType = "image/png",
    width = width,
    height = height,
    alt = "Dummy should be here"
  ))
})
