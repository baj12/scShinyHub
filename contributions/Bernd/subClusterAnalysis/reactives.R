selectedDge <- reactiveValues()

dge_func <- function(projections, log2cpm, featureData, dbCluster, cl1, db1, db2){
  if(DEBUG)cat(file=stderr(), "dge1\n")
  subsetData <- subset(projections, dbCluster %in% cl1)
  if(DEBUG)cat(file=stderr(), "dge2\n")
  cells.1 <- rownames(brushedPoints(subsetData, db1))
  if(DEBUG)cat(file=stderr(), "dge3\n")
  cells.2 <- rownames(brushedPoints(subsetData, db2))
  
  if(DEBUG)cat(file=stderr(), "dge4\n")
  subsetExpression <- log2cpm[, union(cells.1, cells.2)]
  
  if(DEBUG)cat(file=stderr(), "dge5\n")
  genes.use <- rownames(subsetExpression)
  if(DEBUG)cat(file=stderr(), "dge6\n")
  data.1 = apply(subsetExpression[genes.use, cells.1], 1, expMean)
  if(DEBUG)cat(file=stderr(), "dge7\n")
  data.2 = apply(subsetExpression[genes.use, cells.2], 1, expMean)
  if(DEBUG)cat(file=stderr(), "dge8\n")
  total.diff = (data.1 - data.2)
  
  if(DEBUG)cat(file=stderr(), "dge9\n")
  genes.diff = names(which(abs(total.diff) > .2))
  if(DEBUG)cat(file=stderr(), "dge10\n")
  genes.use = ainb(genes.use, genes.diff)
  
  if(DEBUG)cat(file=stderr(), "dge11\n")
  toReturn <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  if(DEBUG)cat(file=stderr(), "dge12\n")
  toReturn[, "avg_diff"] = total.diff[rownames(toReturn)]
  if(DEBUG)cat(file=stderr(), "dge13\n")
  toReturn$Associated.Gene.Name <-
    featureData[rownames(toReturn), 'Associated.Gene.Name']
  return(toReturn)
}


dge <- reactive({
  if(DEBUG)cat(file=stderr(), "dge\n")
  featureData = featureDataReact()
  log2cpm = log2cpm()
  projections = projections()
  cl1 = input$clusters1
  db1 = input$db1
  db2 = input$db2
  # dbcl = dbCluster
  if(is.null(featureData) | is.null(log2cpm) | is.null(projections)){
    return(NULL)
  }
  if(!is.null(getDefaultReactiveDomain())){
    showNotification("running dge", id="dge", duration = NULL)
  }
  if(DEBUGSAVE) 
    save(file = "~/scShinyHubDebug/dge.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file='~/scShinyHubDebug/dge.RData')

  toReturn = dge_func(projections, log2cpm, featureData, dbcl, cl1, db1, db2)
  if(DEBUG)cat(file=stderr(), "dge14\n")
  if(nrow(toReturn)==0){
    if(DEBUG)cat(file=stderr(), "dge: nothing found\n")
    if(!is.null(getDefaultReactiveDomain())){
      showNotification("dge: nothing found", id="dgewarning", duration = 10, type = 'warning')
    }
  }
  if(DEBUG)cat(file=stderr(), "dge15\n")
  selectedDge <- toReturn
  cat(stderr(), rownames(toReturn)[1:5])
  if(DEBUG)cat(file=stderr(), "dge: done\n")
  if(!is.null(getDefaultReactiveDomain())){
    removeNotification( id="dge")
  }
  
  return(toReturn)
  
})
