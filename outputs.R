# SUMMARY STATS ----------------------------------------------------------------
source("moduleServer.R", local=TRUE)
source("reactives.R", local=TRUE)


#################################
# Parameters / normalization 
output$normalizationRadioButtonValue <- renderPrint({ input$normalizationRadioButton })

normaliztionParameters = list(raw = "no Parameters needed")
parFiles = dir(path = "contributions", pattern = "parameters.R", full.names = TRUE, recursive = TRUE)
for(fp in parFiles){
  myNormalizationParameters = list()
  source(fp, local = TRUE)
  if(DEBUGSAVE)
    save(file = "~/scShinyHubDebug/normalizationsParameters.RData", 
         list = c("normaliztionParameters",ls(),ls(envir = globalenv())))
  # load(file = "~/scShinyHubDebug/normalizationsParameters.RData")
  
  for (li in 1:length(myNormalizationParameters)){
    lVal = myNormalizationParameters[[li]]
    if(length(lVal)>0){
      if(DEBUG){
        cat(file=stderr(), paste("normalization Choice: ", 
                                        names(myNormalizationParameters)[li], " ",
                   lVal, "\n"))
        cat(file=stderr(), paste("class: ", 
                                 class(myNormalizationParameters[[li]]), " ",
                                 lVal, "\n"))
      }
      oldNames = names(normaliztionParameters)
      normaliztionParameters[[length(normaliztionParameters) + 1]] = lVal
      names(normaliztionParameters ) = c(oldNames, names(myNormalizationParameters)[li])
    }
  }
}

output$normalizationsParametersDynamic <- renderUI({
  if(is.null(input$normalizationRadioButton))
    return(NULL)
  selectedChoice = input$normalizationRadioButton
  
  if (DEBUG) {
    cat(file=stderr(), paste(class(normaliztionParameters)),"\n")
    cat(file=stderr(), paste(length(normaliztionParameters)),"\n")
    for (li in normaliztionParameters){
      if(length(li)>0){
        if(DEBUG)cat(file=stderr(), paste("normaliztionParameters: ", li, "\n"))
        if(DEBUG)cat(file=stderr(), paste("normaliztionParameters: ", names(li), "\n"))
      }
    }
  }
   if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/normalizationsParametersDynamic.RData", 
         list = c("normaliztionParameters",ls(),ls(envir = globalenv())))
  # load(file = "~/scShinyHubDebug/normalizationsParametersDynamic.RData")
  do.call("switch", args = c(selectedChoice,
         normaliztionParameters,
          h3('no parameters provided')
  ))
})

# End of Parameters / normalization
#################################

output$summaryStatsSideBar<-renderUI({
  if(DEBUG)cat(file=stderr(), "output$summaryStatsSideBar\n")
  gbm = gbm_matrix()
  if(is.null(gbm) ){
    if(DEBUG)cat(file=stderr(), "output$summaryStatsSideBar:NULL\n")
    return(NULL)
  }
  if(input$noStats){
    if(DEBUG)cat(file=stderr(), "output$summaryStatsSideBar:off\n")
    return(NULL)
  }
  
  line1<-paste('No. of cells: ', dim(gbm)[2],sep='\t')
  line2<-paste('No. of genes: ' ,  dim(gbm)[1],sep='\t')
  line3<-paste('Median UMIs per cell: ', medianUMI(),sep='\t')
  line4<-paste('Median Genes with min 1 UMI: ', medianENSG(),sep='\t')
  line5<-paste('Total number of reads: ' , sum(gbm))
  line6<-paste("Memory used:", getMemoryUsed())
  HTML(
    paste0("Summary statistics of this dataset:", '<br/>','<br/>',
           line1, '<br/>',
           line2, '<br/>',
           line3, '<br/>',
           line4, '<br/>',
           line5, '<br/>',
           line6
    )
  )
})


# Select Genes ----------------------------------------------------------------
# this is part of the basic functionality from this tools and thus, can stay in this file.
output$geneListSelection <- renderTree({ 
  geneLists
})


# ONOFF TAB RENDER TABLE ALL CELLS ------------------------------------------------------------------
# TODO module for DT
# this is part of the basic functionality from this tools and thus, can stay in this file.
output$selectedGenesTable <- DT::renderDataTable({
  if(DEBUG)cat(file=stderr(), "output$selectedGenesTable\n")
  dataTables = inputData()
  useGenes = useGenes()
  useCells = useCells()
  if(is.null(dataTables) | is.null(useGenes) | is.null(useCells))
    return(NULL)
  
  gbm = as.matrix(exprs(dataTables$gbm))
  fd = dataTables$featuredata
  dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
  dt$rowSums = rowSums(gbm[useGenes,useCells])
  dt$rowSamples = rowSums(gbm[useGenes,useCells]>0)
  DT::datatable(dt)
})

# TODO module for DT
# TODO move to were it belongs  
output$removedGenesTable <- DT::renderDataTable({
  if(DEBUG)cat(file=stderr(), "output$removedGenesTable\n")
  dataTables = inputData()
  useGenes = !useGenes()
  useCells = useCells()
  if(is.null(dataTables) | is.null(useGenes)  | is.null(useCells) )
    return(NULL)
  
  gbm = as.matrix(exprs(dataTables$gbm))
  fd = dataTables$featuredata
  dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
  dt$rowSums = rowSums(gbm[useGenes,useCells])
  dt$rowSamples = rowSums(gbm[useGenes,useCells]>0)
  DT::datatable(dt)
})


# TODO module of DT with selected names above
# Print names of selected genes for gene selection above table
output$gsSelectedGenes <- renderText({
  if(DEBUG)cat(file=stderr(), "gsSelectedGenes\n")
  dataTables = inputData()
  useGenes = useGenes()
  useCells = useCells()
  if(is.null(dataTables) | is.null(useGenes) | is.null(useCells))
    return(NULL)
  
  gbm = as.matrix(exprs(dataTables$gbm))
  fd = dataTables$featuredata
  dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
  
  paste0(dt$Associated.Gene.Name[input$selectedGenesTable_rows_selected],",")
})

# Print names of removed genes for gene selection  
output$gsrmGenes <- renderText({
  if(DEBUG)cat(file=stderr(), "gsrmGenes\n")
  dataTables = inputData()
  useGenes = useGenes()
  useCells = useCells()
  if(is.null(dataTables) | is.null(useGenes) | is.null(useCells))
    return(NULL)
  
  gbm = as.matrix(exprs(dataTables$gbm))
  fd = dataTables$featuredata
  dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
  
  paste0(dt$Associated.Gene.Name[input$removedGenesTable_rows_selected],",")
})

output$DEBUGSAVEstring <-  renderText({ DEBUGSAVE <<- input$DEBUGSAVE })

r<-callModule(tableSelectionServer, "cellSelectionMod", inputSample)

r<-callModule(tableSelectionServer, "normalizationResult", gbmLogMatrixDisplay)


