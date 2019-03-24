require(ggplot2)

scaterPNG <- reactive({
  if (DEBUG) cat(file = stderr(), "scaterPNG\n")
  start.time <- base::Sys.time()
  scaterReads <- scaterReads()
  if (is.null(scaterReads)) {
    return(NULL)
  }


  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scater.Rmd", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/scater.Rmd')


  if (is.null(width)) {
    width <- 96 * 7
  }
  if (is.null(height)) {
    height <- 96 * 7
  }

  myPNGwidth <- width / 96
  myPNGheight <- height / 96

  outfile <- paste0(tempdir(), "/scaterPlot.png")
  # outfile <- paste0("~/scShinyHubDebug",'/scaterPlot.png')
  if (DEBUG) cat(file = stderr(), paste("output file: ", outfile, "\n"))
  if (DEBUG) cat(file = stderr(), paste("output file normalized: ", normalizePath(outfile, mustWork = FALSE), "\n"))
  n <- min(nrow(scaterReads), 50)
  # use plotHighestExprs instead of plotQC
  # p1 <- scater::plotQC(scaterReads, type = "highest-expression", colour_cells_by = "fixed", n = n)
  p1 <- scater::plotHighestExprs(scaterReads, colour_cells_by = "log10_total_counts", n=n)
  tryCatch(
    ggsave(file = normalizePath(outfile, mustWork = FALSE), plot = p1, width = myPNGwidth, height = myPNGheight, units = "in"),
    error = function(e) {
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("Problem saving ggplot", type = "warning", duration = NULL)
      }
      return(NULL)
    }
  )

  if (DEBUG) cat(file = stderr(), "done:scaterPNG\n")
  end.time <- base::Sys.time()
  cat(file = stderr(), paste("this took: ", difftime(end.time, start.time, units = "min"), " min\n"))
  
  return(list(
    src = normalizePath(outfile, mustWork = FALSE),
    contentType = "image/png",
    width = width,
    height = height,
    alt = "Scater plot should be here"
  ))
})
