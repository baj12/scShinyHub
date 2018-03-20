

# LIBRARY -----------------------------------------------------------------

library(shiny)
library(shinyTree)
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
  
  
  # TODO integrate with plugin structure
  # create template for input parameters
  # source other reports
  output$report <- downloadHandler(
    filename = "report.html",
    
    content = function(file) {
      
      ip = inputData()
      if( is.null(ip) ){
        if(DEBUG)cat(file=stderr(), "output$report:NULL\n")
        return(NULL)
      }
      
      
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      tempServerFunctions <- file.path(tempdir(), "serverFunctions.R")
      file.copy("serverFunctions.R", tempServerFunctions, overwrite = TRUE)
      tempprivatePlotFunctions <- file.path(tempdir(), "privatePlotFunctions.R")
      file.copy("privatePlotFunctions.R", tempprivatePlotFunctions, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        tempServerFunctions = tempServerFunctions,
        tempprivatePlotFunctions = tempprivatePlotFunctions,
        b1 = input$b1,
        cluster = input$cluster,
        cluster5 = input$cluster5,
        clusters = input$clusters,
        clusters1 = input$clusters1,
        clusters2 = input$clusters2,
        clusters3 = input$clusters3,
        clusters4 = input$clusters4,
        db1 = input$db1,
        db2 = input$db2,
        dimension_x = input$dimension_x,
        dimension_x1 = input$dimension_x1,
        dimension_x2 = input$dimension_x2,
        dimension_x3 = input$dimension_x3,
        dimension_x4 = input$dimension_x4,
        dimension_y = input$dimension_y,
        dimension_y1 = input$dimension_y1,
        dimension_y2 = input$dimension_y2,
        dimension_y3 = input$dimension_y3,
        dimension_y4 = input$dimension_y4,
        file1 = input$file1,
        gene_id = input$gene_id,
        gene_id_sch = input$gene_id_sch,
        geneListSelection = input$geneListSelection,
        heatmap_geneids = input$heatmap_geneids,
        heatmap_geneids2 = input$heatmap_geneids2,
        maxGenes = input$maxGenes,
        mclustids = input$mclustids,
        minExpGenes = input$minExpGenes,
        minGenes = input$minGenes,
        minGenesGS = input$minGenesGS,
        panelplotids = input$panelplotids,
        # positiveCells = positiveCells$positiveCells,
        # positiveCellsAll = positiveCells$positiveCellsAll,
        scb1 = input$scb1,
        selectIds = input$selectIds
      )
      if(DEBUG)cat(file=stderr(), "output$report:gbm:\n")
      if(DEBUG)cat(file=stderr(), str(gbm))
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





