# SUMMARY STATS ----------------------------------------------------------------
source("moduleServer.R", local=TRUE)
source("reactives.R", local=TRUE)

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
  HTML(
    paste0("Summary statistics of this dataset:", '<br/>','<br/>',
           line1, '<br/>',
           line2, '<br/>',
           line3, '<br/>',
           line4, '<br/>',
           line5
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
  dataTables = inputData$inputData
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
  dataTables = inputData$inputData
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
  dataTables = inputData$inputData
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
  dataTables = inputData$inputData
  useGenes = !useGenes()
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



