# LIBRARIES -----------------------------------------------------------------
library(shiny)
# library(reactlog)
library(shinyTree)
library(tibble)
library(shinyBS)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)
library(pheatmap)
library(threejs)
library(sm)
library(RColorBrewer)
library(mclust)
library(reshape)
library(ggplot2)
library(knitr)
library(kableExtra)
library(shinyWidgets)
library(scater)
library(shinyMCE)
library(kohonen)
library(Rsomoclu)
library(gtools)
library(SingleCellExperiment)
library(Matrix)
library(colourpicker)
library(shinytest)
library(scran)
library(callr)
library(debugme)


if (file.exists("defaultValues.R")) {
  base::source(file = "defaultValues.R")
} else {
  base::warning("no defaultsValues.R file")
  base::stop("stop")
}

# list available colors for samples and clusters, other colors are defined independantly.
if (!exists("allowedColors")) {
  allowedColors = unique(c("#8c510a","#d8b365","#f6e8c3","#c7eae5","#5ab4ac","#01665e","#c51b7d","#e9a3c9",
                           "#fde0ef","#e6f5d0","#a1d76a","#4d9221","#762a83","#af8dc3","#e7d4e8","#d9f0d3",
                           "#7fbf7b","#1b7837","#b35806","#f1a340","#fee0b6","#d8daeb","#998ec3","#542788",
                           "#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac","#b2182b","#ef8a62",
                           "#fddbc7","#e0e0e0","#999999","#4d4d4d"))
}
Sys.setenv(DEBUGME = ".")
base::source("serverFunctions.R")

# TODO needs to be an option
seed <- 2

# enableBookmarking(store = "server")

# global variable with directory where to store files to be included in reports
reportTempDir <<- base::tempdir()

shinyServer(function(input, output, session) {
  base::set.seed(seed)
  # check that directory is availabl, otherwise create it
  if (DEBUG) {
    if (!dir.exists("~/scShinyHubDebug")) {
      base::dir.create("~/scShinyHubDebug")
    }
    # TODO ??? clean directory??
  }
  
  # files to be included in report
  # developers can add in outputs.R a variable called "myZippedReportFiles"
  zippedReportFiles <- c("report.html", "sessionData.RData", 
                         "normalizedCounts.csv", "variables.used.txt", 
                         "inputUsed.RData")

  base::options(shiny.maxRequestSize = 2000 * 1024^2)
  
  # TODO check if file exists
  # TODO have this as an option to load other files
  base::load(file = "geneLists.RData")
  
  if (DEBUG) base::cat(file = stderr(), "ShinyServer running\n")
  # base calculations that are quite expensive to calculate
  # display name, reactive name to be executed
  heavyCalculations <- list(
    c("pca", "pca"),
    c("scran_Cluster", "scran_Cluster"),
    c("projections", "projections")
  )
  
  # base projections
  # display name, reactive to calculate projections
  projectionFunctions <<- list(
    c("sampleNames", "sample"),
    c("Gene count", "geneCount"),
    c("UMI count", "umiCount"),
    c("before filter", "beforeFilterPrj")
  )
  
  # differential expression functions
  # used in subcluster analysis
  diffExpFunctions <<- list()
  
  # load global reactives, modules, etc ----
  base::source("reactives.R", local = TRUE)
  base::source("outputs.R", local = TRUE)
  base::source("modulesUI.R", local = TRUE)
  base::source("moduleServer.R", local = TRUE)
  
  # bookmarking ----
  # couldn't get bookmarking to work, esp. with the input file
  # setBookmarkExclude(c("bookmark1"))
  # observeEvent(input$bookmark1, {
  #   if (DEBUG) cat(file = stderr(), paste("bookmarking: \n"))
  #   if (DEBUG) cat(file = stderr(), paste(names(input), collapse = "\n"))
  #
  #   session$doBookmark()
  #   if (DEBUG) cat(file = stderr(), paste("bookmarking: DONE\n"))
  # })
  # Need to exclude the buttons from themselves being bookmarked
  
  # load contribution reactives ----
  # parse all reactives.R files under contributions to include in application
  uiFiles <- base::dir(
    path = "contributions", pattern = "reactives.R",
    full.names = TRUE, recursive = TRUE
  )
  for (fp in uiFiles) {
    if (DEBUG) base::cat(file = stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations <- NULL
    myProjections <- NULL
    myDiffExpFunctions <- NULL
    base::source(fp, local = TRUE)
    
    heavyCalculations <- append2list(myHeavyCalculations, heavyCalculations)
    projectionFunctions <- append2list(myProjections, projectionFunctions)
    diffExpFunctions <- append2list(myDiffExpFunctions, diffExpFunctions)
  }
  
  # update diffExpression radiobutton
  dgeChoices = c()
  if (length(diffExpFunctions) > 0) {
    for (li in 1:length(diffExpFunctions)) {
      liVal <- diffExpFunctions[[li]]
      if (length(liVal) == 2) {
        dgeChoices = c(dgeChoices, liVal[1])
      } else {
        # shouldn't happen
        error("number of values for normalization function is not 2\n")
      }
    }
  }
  updateRadioButtons(session = session, inputId = "sCA_dgeRadioButton",
                     choices = dgeChoices)
  # make variable global
  diffExpFunctions <<- diffExpFunctions
  
  # load contribution outputs ----
  # parse all outputs.R files under contributions to include in application
  uiFiles <- base::dir(path = "contributions", pattern = "outputs.R", 
                       full.names = TRUE, recursive = TRUE)
  for (fp in uiFiles) {
    if (DEBUG) cat(file = stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations <- NULL
    myProjections <- NULL
    myZippedReportFiles <- c()
    base::source(fp, local = TRUE)
    heavyCalculations <- append2list(myHeavyCalculations, heavyCalculations)
    projectionFunctions <<- append2list(myProjections, projectionFunctions)
    zippedReportFiles <- c(zippedReportFiles, myZippedReportFiles)
  }

  
}) # END SERVER

# shiny::showReactLog()

# enableBookmarking(store = "server")
