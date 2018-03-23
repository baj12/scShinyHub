

# LIBRARY -----------------------------------------------------------------

library(shiny)
library(shinyTree)
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
library(cellrangerRkit)
library(SCORPIUS)
library(ggplot2)
library(knitr)
library(kableExtra)
library(shinyWidgets)
library(scater)

source("serverFunctions.R")
source("privatePlotFunctions.R")
DEBUG=TRUE

#needs to be an option
seed=1

shinyServer(function(input, output, session) {
  # TODO create a UI element for seed
  set.seed(seed)
  
  
  options(shiny.maxRequestSize = 2000 * 1024 ^ 2)
  
  # TODO check if file exists
  # TODO have this as an option to load other files
  load(file = "geneLists.RData")
  
  if(DEBUG) cat(file=stderr(), "ShinyServer running\n")
  
  # base calculations that are quite expensive to calculate
  heavyCalculations = list(c("pca", "pca"),
                           c("kmClustering", "kmClustering"),
                           c("tsne", "tsne")
  )
  
  # ------------------------------------------------------------------------------------------------------------
  # load global reactives, modules, etc
  source("reactives.R", local = TRUE)
  source("outputs.R", local = TRUE)
  source("modulesUI.R", local = TRUE)  
  source("moduleServer.R", local = TRUE)
  
  
  # ------------------------------------------------------------------------------------------------------------
  # load contribution reactives
  # parse all reactives.R files under contributions to include in application
  uiFiles = dir(path = "contributions", pattern = "reactives.R", full.names = TRUE, recursive = TRUE)
  for(fp in uiFiles){
    if(DEBUG)cat(file=stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations = NULL
    source(fp, local = TRUE)
    heavyCalculations = appendHeavyCalculations(myHeavyCalculations, heavyCalculations)
  }
  # load contribution outputs
  # parse all outputs.R files under contributions to include in application
  uiFiles = dir(path = "contributions", pattern = "outputs.R", full.names = TRUE, recursive = TRUE)
  for(fp in uiFiles){
    if(DEBUG)cat(file=stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations = NULL
    source(fp, local = TRUE)
    heavyCalculations = appendHeavyCalculations(myHeavyCalculations, heavyCalculations)
  }
  
  positiveCells <- reactiveValues(positiveCells = NULL,
                                  positiveCellsAll = NULL)
  selectedDge <- reactiveValues()
  
  # ------------------------------------------------------------------------------------------------------------
  # handling expensive calcualtions
  forceCalc <-observe({
    input$goCalc
    isolate({
      if(DEBUG)cat(file=stderr(), "forceCalc\n")
      # list of output variable and function name
      
      withProgress(message = 'Performing heavy calculations', value = 0, {
        n=length(heavyCalculations)
        for(calc in heavyCalculations){
          incProgress(1/n, detail = paste("Creating ", calc[1]))
          if(DEBUG)cat(file=stderr(), paste("forceCalc ", calc[1],"\n"))
          assign(calc[1], eval(parse(text = paste0( calc[2], "()"))))
        }
      })
    })
  })
  
  # ------------------------------------------------------------------------------------------------------------
  
  
  # TODO as module
  # sub cluster analysis ( used for 2 panels )
  output$clusters1 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters1\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectizeInput(
        "clusters1",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0, 
        multiple = TRUE
      )
    }
  })
  
  # # TODO as module
  # output$clusters2 <- renderUI({
  #   if(DEBUG)cat(file=stderr(), "output$clusters2\n")
  #   tsne.data = tsne.data()
  #   if(is.null(tsne.data)){
  #     HTML("Please load data firts")
  #   }else{
  #     noOfClusters <- max(tsne.data$dbCluster)
  #     selectizeInput(
  #       "clusters2",
  #       label = "Cluster",
  #       choices = c(0:noOfClusters),
  #       selected = 0, 
  #       multiple = TRUE
  #     )
  #   }
  # })
  
  # TODO as module
  # coexpression binarized
  output$clusters3 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters3\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectizeInput(
        "clusters3",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0, 
        multiple = TRUE
      )
    }
  })
  
  # TODO as module
  # data expression panel plot 
  output$clusters4 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters4\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectInput(
        "clusters4",
        label = "Cluster",
        choices = c(c('All'),c(0:noOfClusters)),
        selected = 0
      )
    }
  })  
  
  
  
  
  # Report creation ------------------------------------------------------------------
  output$report <- downloadHandler(
    filename = "report.html",
    
    content = function(file) {
      
      ip = inputData()
      if( is.null(ip) ){
        if(DEBUG)cat(file=stderr(), "output$report:NULL\n")
        return(NULL)
      }
      tDir = tempdir()

      # ------------------------------------------------------------------------------------------------------------
      # the reactive.R cam hold functions that can be used in the report to reduce the possibility of code replication
      # we copy them to the temp directory and load them in the markdown
      uiFiles = dir(path = "contributions", pattern = "reactives.R", full.names = TRUE, recursive = TRUE)
      reactiveFiles = ""
      for(fp in c("reactives.R", uiFiles)){
        if(DEBUG)cat(file=stderr(), paste("loading: ", fp, "\n"))
        tmpFile = tempfile(pattern = "file", tmpdir = tDir, fileext = ".R")
        file.copy(fp, tmpFile, overwrite = TRUE)
        reactiveFiles = paste0(reactiveFiles, "source(\"", tmpFile,"\")\n", collapse = "\n")
      }
      reactiveFiles = paste0("\n\n```{r load-reactives}\n", reactiveFiles, "\n```\n\n")
      
      
      
            # ------------------------------------------------------------------------------------------------------------
      # handle plugin reports
      # load contribution reports
      # parse all report.Rmd files under contributions to include in application
      uiFiles = dir(path = "contributions", pattern = "report.Rmd", full.names = TRUE, recursive = TRUE)
      pluginReportsString = ""
      fpRidx = 1
      for(fp in uiFiles){
        if(DEBUG)cat(file=stderr(), paste("loading: ", fp, "\n"))
        tmpFile = tempfile(pattern = "file", tmpdir = tDir, fileext = ".Rmd")
        file.copy(fp, tmpFile, overwrite = TRUE)
        pluginReportsString = paste0(pluginReportsString, 
                                     "\n\n```{r child-report-",fpRidx,", child = '" ,tmpFile,"'}\n```\n\n")
        fpRidx = fpRidx + 1
      }
      
      
      
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tDir, "report.Rmd")
      
      tempServerFunctions <- file.path(tDir, "serverFunctions.R")
      file.copy("serverFunctions.R", tempServerFunctions, overwrite = TRUE)
      tempprivatePlotFunctions <- file.path(tDir, "privatePlotFunctions.R")
      file.copy("privatePlotFunctions.R", tempprivatePlotFunctions, overwrite = TRUE)
      
      # create a new list of all parameters that can be passed to the markdown doc.
      inputNames = names(input)
      params <- list(
        tempServerFunctions = tempServerFunctions,
        tempprivatePlotFunctions = tempprivatePlotFunctions,
        calledFromShiny = TRUE # this is to notify the markdown that we are running the script from shiny. used for debugging/development
      )
      for (idx in 1:length(names(input))){
        params[[inputNames[idx]]] = input[[inputNames[idx]]]
      }
      
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # read the template and replace parameters placeholder with list 
      # of paramters
      x <- readLines(tempReport)
      # x <- readLines("report.Rmd")
      paramString = paste0("  ", names(params),": NA", collapse = "\n")
      y <- gsub( "#__PARAMPLACEHOLDER__", paramString, x )
      y <- gsub( "__CHILDREPORTS__", pluginReportsString, y )
      y <- gsub( "__LOAD_REACTIVES__", reactiveFiles, y )
      # cat(y, file="tempReport.Rmd", sep="\n")
      cat(y, file=tempReport, sep="\n")
      
      if(DEBUG)cat(file=stderr(), "output$report:gbm:\n")
      if(DEBUG)cat(file=stderr(), paste("\n", tempReport,"\n"))
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
})# END SERVER





