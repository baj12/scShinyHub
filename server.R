

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

if (file.exists("defaultValues.R")) {
  source(file = "defaultValues.R")
} else {
  warning("no defaultsValues.R file")
  stop("stop")
} 

source("serverFunctions.R")


# create large example files from split
# this is needed to overcome the size limit in GitHub
if (!file.exists("Examples/PBMC-Apheresis.new.Rds")) {
  xaaName <- "Examples/PBMC.xaa"
  xabName <- "Examples/PBMC.xab"
  contents <- readBin(xaaName, "raw", file.info(xaaName)$size)
  contents2 <- readBin(xabName, "raw", file.info(xabName)$size)
  outFile <- file("Examples/PBMC-Apheresis.new.Rds", 'ab')
  writeBin(contents, outFile)
  writeBin(contents2, outFile )
  close(outFile)
}

#needs to be an option
seed=2

enableBookmarking(store = "server")

shinyServer(function(input, output, session) {
  
  # TODO-BJ create a UI element for seed
  set.seed(seed)
  # check that directory is availabl, otherwise create it
  if(DEBUG){
    if(!dir.exists("~/scShinyHubDebug")){
      dir.create("~/scShinyHubDebug")
    }
    # TODO ??? clean directory??
  }
  
  options(shiny.maxRequestSize = 2000 * 1024 ^ 2)
  
  # TODO check if file exists
  # TODO have this as an option to load other files
  load(file = "geneLists.RData")
  
  if(DEBUG) cat(file=stderr(), "ShinyServer running\n")
  
  # base calculations that are quite expensive to calculate
  # display name, reactive name to be executed
  heavyCalculations = list(c("pca", "pca"),
                           c("kmClustering", "kmClustering"),
                           c("projections", "projections")
                           
  )
  
  # base projections 
  # display name, reactive to calculate projections
  projectionFunctions <<- list(
                             c("sampleNames", "sample"),
                             c("Gene count", "geneCount"),
                             c("UMI count", "umiCount")
  )
  
  
  
  # ------------------------------------------------------------------------------------------------------------
  # load global reactives, modules, etc
  # why not import them  earlier? I rember that there was an issue. could be documented
  source("reactives.R", local = TRUE)
  source("outputs.R", local = TRUE)
  source("modulesUI.R", local = TRUE)  
  source("moduleServer.R", local = TRUE)
  
  # ------------------------------------------------------------------------------------------------------------
  # bookmarking
  setBookmarkExclude(c("bookmark1"))
  observeEvent(input$bookmark1, {
    if(DEBUG)cat(file=stderr(), paste("bookmarking: \n"))
    if(DEBUG)cat(file=stderr(), paste(names(input), collapse = "\n"))
    
    session$doBookmark()
    if(DEBUG)cat(file=stderr(), paste("bookmarking: DONE\n"))
  })
  # Need to exclude the buttons from themselves being bookmarked
  
  # ------------------------------------------------------------------------------------------------------------
  # load contribution reactives
  # parse all reactives.R files under contributions to include in application
  uiFiles = dir(path = "contributions", pattern = "reactives.R", 
                full.names = TRUE, recursive = TRUE)
  for(fp in uiFiles){
    if(DEBUG)cat(file=stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations = NULL
    myProjections = NULL
    source(fp, local = TRUE)
    heavyCalculations = appendHeavyCalculations(myHeavyCalculations, heavyCalculations)
    projectionFunctions = appendHeavyCalculations(myProjections, projectionFunctions)
  }
  # load contribution outputs
  # parse all outputs.R files under contributions to include in application
  uiFiles = dir(path = "contributions", pattern = "outputs.R", full.names = TRUE, recursive = TRUE)
  for(fp in uiFiles){
    if(DEBUG)cat(file=stderr(), paste("loading: ", fp, "\n"))
    myHeavyCalculations = NULL
    myProjections = NULL
    source(fp, local = TRUE)
    heavyCalculations = appendHeavyCalculations(myHeavyCalculations, heavyCalculations)
    projectionFunctions <<- appendHeavyCalculations(myProjections, projectionFunctions)
  }
  
  # TODO move coexpression for binarized needed
  # in reactives., report, server, coexpression/output
  positiveCells <- reactiveValues(positiveCells = NULL,
                                  positiveCellsAll = NULL)
  
  # ------------------------------------------------------------------------------------------------------------
  # handling expensive calcualtions
  forceCalc <- observe({
    input$goCalc
    start.time <- Sys.time()
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
    end.time <- Sys.time()
    # tfmt <- "%Hh %Mm %Ss"
    # t1 <- strptime(end.time - start.time, format=tfmt)
    cat(file=stderr(), paste("this took: ", difftime(end.time, start.time, units="min") ," min\n"))
    updateMemUse$update = isolate(updateMemUse$update) + 1
  })
  
  
  output$countscsv <- downloadHandler(
    filename = paste0("counts.",Sys.Date(),".csv"),
    content = function(file){
      if(DEBUG)cat(file=stderr(), paste("countcsv: \n"))
      gbmlog = gbm_log()
      if(is.null(gbmlog)){
        return(NULL)
      }
      write.csv(as.matrix(exprs(gbmlog)), file)
    }
  )
  
  output$RDSsave<- downloadHandler(
    filename = paste0("project.",Sys.Date(),".Rds"),
    content = function(file){
      if(DEBUG)cat(file=stderr(), paste("RDSsave: \n"))
      gbm = gbm()
      featuredata <- featureDataReact()
      
      if(is.null(gbm) | is.null(featuredata)){
        return(NULL)
      }
      if (DEBUGSAVE)
        save(file = "~/scShinyHubDebug/RDSsave.RData", list = c(ls(),ls(envir = globalenv())))
      # load(file='~/scShinyHubDebug/RDSsave.RData')
      
      save(file = file, list = c("featuredata",  "gbm"))
      # write.csv(as.matrix(exprs(gbm)), file)
    }
  )
  
  # Report creation ------------------------------------------------------------------
  output$report <- downloadHandler(
    filename = "report.html",
    
    content = function(file) {
      
      ip = inputData()
      if( is.null(ip) ){
        if(DEBUG) cat(file=stderr(), "output$report:NULL\n")
        return(NULL)
      }
      tDir = tempdir()
      reactiveFiles = ""
      
      #-----------
      # fixed files
      tmpFile = tempfile(pattern = "file", tmpdir = tDir, fileext = ".RData")
      file.copy("geneLists.RData", tmpFile, overwrite = TRUE)
      reactiveFiles = paste0(reactiveFiles, "load(file=\"", tmpFile,"\")\n", collapse = "\n")
      
      #----- Projections
      # projections can contain mannually annotated groups of cells and different normalizations.
      # to reduce complexity we are going to save those in a separate RData file
      tmpPrjFile = tempfile(pattern = 'file', tmpdir = tDir, fileext = ".RData")
      projections = projections()
      gbm_log = gbm_log()
      gNames = groupNames$namesDF
      base::save(file = tmpPrjFile, list = c("projections", "gbm_log", "gNames"))
      
      # ------------------------------------------------------------------------------------------------------------
      # the reactive.R can hold functions that can be used in the report to reduce the possibility of code replication
      # we copy them to the temp directory and load them in the markdown
      uiFiles = dir(path = "contributions", pattern = "reactives.R", full.names = TRUE, recursive = TRUE)
      for (fp in c("reactives.R", uiFiles)) {
        if (DEBUG) cat(file = stderr(), paste("loading: ", fp, "\n"))
        tmpFile = tempfile(pattern = "file", tmpdir = tDir, fileext = ".R")
        file.copy(fp, tmpFile, overwrite = TRUE)
        reactiveFiles = paste0(reactiveFiles, "source(\"", tmpFile,"\")\n", collapse = "\n")
      }
      # otherwise reactie might overwrite projections...
      reactiveFiles = paste0(reactiveFiles, "load(file=\"", tmpPrjFile,"\")\n", collapse = "\n")
      # encapsulte the load files in an R block
      reactiveFiles = paste0("\n\n```{r load-reactives, include=FALSE}\n", reactiveFiles, "\n```\n\n")
      
      
      
      
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

      # create a new list of all parameters that can be passed to the markdown doc.
      inputNames = names(input)
      params <- list(
        tempServerFunctions = tempServerFunctions,
        # tempprivatePlotFunctions = tempprivatePlotFunctions,
        calledFromShiny = TRUE # this is to notify the markdown that we are running the script from shiny. used for debugging/development
      )
      for (idx in 1:length(names(input))) {
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
      cat(y, file = tempReport, sep = "\n")
      
      if (DEBUG) cat(file = stderr(), "output$report:gbm:\n")
      if (DEBUG) cat(file = stderr(), paste("\n", tempReport,"\n"))
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app)
      renderEnv = new.env(parent = globalenv())
      if (DEBUG) file.copy(tempReport, '~/scShinyHubDebug/tempReport.Rmd')
      myparams = params
      if (DEBUG) save(file = '~/scShinyHubDebug/tempReport.RData', list = c("myparams","renderEnv"))
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = renderEnv
      )
    }
  )
  
})# END SERVER


# enableBookmarking(store = "server")
