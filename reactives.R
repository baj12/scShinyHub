# general reactive 
# (same as global variables, but they are reactive to manipulation and are lazy, i.e. they only get executed when needed.)

# reactive values  ------------------------------------------------------------------
# should only hold original data
inputDataFunc <- function(inFile){
  if(DEBUG)cat(file=stderr(), "DEBUG:inputData\n")
  
  load(inFile$datapath)
  
  cat(stderr(), 'Loaded')
  dataTables=list()
  dataTables$log2cpmOrg <- log2cpm
  dataTables$tsne.data <- tsne.data
  dataTables$tsne.dataOrg <- tsne.data
  featuredata <-
    featuredata[which(featuredata$Chromosome.Name %in% c(unlist(lapply(
      seq(1, 22, 1), toString
    )), c("X", "Y", "MT", "N"))), ]
  featuredata$Associated.Gene.Name <-
    toupper(featuredata$Associated.Gene.Name)
  featuredata<-featuredata[rownames(log2cpm),]
  dataTables$featuredataOrg <- featuredata
  dataTables$positiveCells <- NULL
  dataTables$positiveCellsAll <- NULL
  
  # take only genes that are in all tables
  rnames = rownames(featuredata)
  rnames = rnames[rnames %in% rownames(log2cpm)]
  rnames = rnames[rnames %in% rownames(gbm)]
  
  dataTables$log2cpm <- log2cpm[rnames,]
  dataTables$gbm = gbm[rnames,]
  dataTables$gbm_log = gbm_log[rnames,]
  dataTables$featuredata <- featuredata[rnames,]
  
  # some checks
  if(sum(is.infinite(as.matrix(exprs(gbm))))>0){
    showNotification("gbm contains infinite values")
    return(NULL)
  }
  if(DEBUG)cat(file=stderr(), "inputData: done\n")
  return(dataTables)
  
}

inputData = reactive({
  inFile <- input$file1
  if (is.null(inFile)){
    if(DEBUG)cat(file=stderr(), "inputData: NULL\n")
    return(NULL)
  }
  return(inputDataFunc(inFile))
})

medianENSGfunc <- function(log2cpm){
  geneC = colSums(log2cpm>0, na.rm = TRUE)
  if(DEBUG)cat(file=stderr(), "medianENSG:done\n")
  return(median(t(geneC)))
}

medianENSG <- reactive({
  if(DEBUG)cat(file=stderr(), "medianENSG\n")
  log2cpm = log2cpm()
  if(is.null(log2cpm) ){
    if(DEBUG)cat(file=stderr(), "medianENSG:NULL\n")
    return(0)
  }
  if(ncol(log2cpm)<=1 | nrow(log2cpm) < 1){
    return(0)
  }
  return(medianENSGfunc(log2cpm))
})

medianUMIfunc <- function(log2cpm){
  umiC = colSums(log2cpm, na.rm = TRUE)
  if(DEBUG)cat(file=stderr(), "medianUMI:done\n")
  return(median(t(umiC)))
}

medianUMI <- reactive({
  if(DEBUG)cat(file=stderr(), "medianUMI\n")
  log2cpm = log2cpm()
  if(is.null(log2cpm) ){
    if(DEBUG)cat(file=stderr(), "medianUMI:NULL\n")
    return(0)
  }
  return(medianUMIfunc(log2cpm))
})

useCellsFunc <- function(dataTables){
  if(DEBUG)cat(file=stderr(), "useCells2\n")
  minG = input$minGenes
  gbm = as.matrix(exprs(dataTables$gbm))
  goodCols = colSums(gbm) > minG
  
  maxG = input$maxGenes
  goodCols = goodCols & ( colSums(gbm,na.rm = TRUE) <= maxG)
  
  geneNames = input$minExpGenes
  ids = which(grepl(geneNames, dataTables$featuredata$Associated.Gene.Name))
  if(length(ids)==1){
    goodCols = goodCols &gbm[ids, ] > 0 
  }else if(length(ids) == 0){
    showNotification("not enough cells, check gene names for min coverage")
    return(NULL)
  }else{
    goodCols = goodCols & (colSums(gbm[ids, ]) > 0)
  }
  if(sum(goodCols)==0) {
    showNotification("not enough cells left")
    return(NULL)
  }
  if(DEBUG)cat(file=stderr(), "useCells:done\n")
  return(goodCols)
}

# collects information from all places where cells are being removed or specified
useCells <- reactive({
  if(DEBUG)cat(file=stderr(), "useCells\n")
  dataTables = inputData()
  if(!exists("dataTables") || is.null(dataTables)){
    if(DEBUG)cat(file=stderr(), "useCells:NULL\n")
    return(NULL)
  }
  return(useCellsFunc(dataTables))
})


featureDataReact = reactive({
  if(DEBUG)cat(file=stderr(), "featureData\n")
  dataTables = inputData()
  useGenes = useGenes()
  if(!exists("dataTables") | is.null(dataTables) | is.null(useGenes)){
    if(DEBUG)cat(file=stderr(), "featureData:NULL\n")
    return(NULL)
  }
  if(DEBUG)cat(file=stderr(), "featureData:done\n")
  return(dataTables$featuredata[useGenes, ])
})

useGenesFunc <- function(dataTables, useCells, ipIDs, geneListSelection, minGene){
  if(nchar(ipIDs)>0) 
  {
    rmIds = !grepl(ipIDs, dataTables$featuredata$Associated.Gene.Name)
  }
  else 
  {
    rmIds = rep(TRUE,nrow(dataTables$log2cpm))
  }
  
  # gene groups to be included
  if(!is.null(geneListSelection)){
    selectedgeneList = get_selected(geneListSelection)
    if(length(selectedgeneList)>0){
      selGenes = c()
      for(sIdx in 1:length(selectedgeneList)){
        print(sIdx)
        att = attr(selectedgeneList[[sIdx]], "ancestry")
        if(length(att)>0){
          selGenes = c(selGenes , geneLists[[att]][[selectedgeneList[[sIdx]]]])
        }
      }
      selGenes = unique(selGenes)
      rmIds = rownames(dataTables$log2cpm) %in% selGenes & rmIds
    }
  }
  
  # overall gene expression Min
  if(!is.null(minGene)){
    selGenes = rowSums(as.matrix(exprs(dataTables$gbm[,useCells]))) >=minGene
    rmIds = rmIds & selGenes
  }
  
  if(DEBUG)cat(file=stderr(), "useGenes: done\n")
  return(rmIds)
  
}

# collects information from all places where genes being removed or specified
useGenes <- reactive({
  if(DEBUG)cat(file=stderr(), "useGenes\n")
  dataTables = inputData()
  useCells = useCells()
  if(!exists("dataTables") | is.null(dataTables) | is.null(useCells) | length(dataTables$featuredata$Associated.Gene.Name)==0){
    if(DEBUG)cat(file=stderr(), "useGenes: NULL\n")
    return(NULL)
  }
  #genes to be explicitly removed
  ipIDs = input$selectIds
  geneListSelection = input$geneListSelection
  minGene <- input$minGenesGS
  return(useGenesFunc(dataTables, useCells, ipIDs, geneListSelection, minGene))
})


# individual values
gbm_log <- reactive({
  if(DEBUG)cat(file=stderr(), "gbm_log\n")
  dataTables = inputData()
  useCells = useCells()
  useGenes = useGenes()
  if(!exists("dataTables") | is.null(dataTables) | is.null(useGenes)){
    if(DEBUG)cat(file=stderr(), "gbm_log:NULL\n")
    return(NULL)
  }
  # gbm rownames are ENSG numbers
  # dataTables$gbm_log[useGenes, useCells]
  if(DEBUG)cat(file=stderr(), "gbm_log:Done\n")
  return(dataTables$gbm_log[useGenes, useCells])
  
})

gbm <- reactive({
  if(DEBUG)cat(file=stderr(), "gbm\n")
  dataTables = inputData()
  useCells = useCells()
  useGenes = useGenes()
  if(!exists("dataTables") | is.null(dataTables) | is.null(useGenes)| is.null(useCells)){
    if(DEBUG)cat(file=stderr(), "gbm: NULL\n")
    return(NULL)
  }
  # gbm rownames are ENSG numbers
  # dataTables$gbm_log[useGenes, useCells]
  if(DEBUG)cat(file=stderr(), "gbm:DONE\n")
  return(dataTables$gbm[useGenes, useCells])
})


log2cpm <- reactive({
  if(DEBUG)cat(file=stderr(), "log2cpm\n")
  dataTables = inputData()
  useCells = useCells()
  useGenes = useGenes()
  if(!exists("dataTables") | is.null(dataTables) | is.null(useGenes)){
    if(DEBUG)cat(file=stderr(), "log2cpm: NULL\n")
    return(NULL)}
  # TODO useGenes/useCells
  if(DEBUG)cat(file=stderr(), "log2cpm:Done\n")
  return(dataTables$log2cpm[useGenes, useCells])
})


pcaFunc <- function(gbm_log){
  if(DEBUG)cat(file=stderr(), "pca:done\n")
  pca = tryCatch({
    run_pca(gbm_log)},
    error = function(e){
      if(!is.null(getDefaultReactiveDomain())){
        showNotification("Problem with PCA, probably not enough cells?")
      }
      return(NULL)
    }
  )
  return(pca)
  
}

pca = reactive({
  if(DEBUG)cat(file=stderr(), "pca\n")
  gbm_log = gbm_log()
  if(is.null(gbm_log)){
    if(DEBUG)cat(file=stderr(), "pca:NULL\n")
    return(NULL)
  }
  return(pcaFunc(gbm_log))
})

kmClustering = reactive({
  if(DEBUG)cat(file=stderr(), "kmClustering\n")
  pca = pca()
  if(is.null(pca)){
    if(DEBUG)cat(file=stderr(), "kmClustering:NULL\n")
    return(NULL)
  }
  clustering=list()
  
  kNr = 10
  for(kNr in 2:10) {
    set.seed(seed = seed)
    km = run_kmeans_clustering(pca, k=kNr)
    clustering[[paste0("kmeans_",kNr,"_clusters")]] = data.frame("Barcode" = rownames(data.frame(km$cluster)), "Cluster" = km$cluster)
  }
  if(DEBUG)cat(file=stderr(), "kmClustering:done\n")
  return(clustering)
})


tsne = reactive({
  if(DEBUG)cat(file=stderr(), "tsne\n")
  pca = pca()
  if(is.null(pca)){
    if(DEBUG)cat(file=stderr(), "tsne: NULL\n")
    return(NULL)
  }
  if(DEBUG)cat(file=stderr(), "tsne: done\n")
  set.seed(seed = seed)
  return(run_tsne(pca, dims = 3, perplexity = 30, theta = 0.5))
})




tsne.data = reactive({
  if(DEBUG)cat(file=stderr(), "tsne.data\n")
  tsne = tsne()
  clustering = kmClustering()
  if(is.null(tsne) | is.null(clustering)){
    if(DEBUG)cat(file=stderr(), "tsne.data: NULL\n")
    return(NULL)
  }
  
  tsne.data = data.frame(tsne$Y)
  colnames(tsne.data) = c("V1", "V2", "V3")
  tsne.data$dbCluster = clustering$kmeans_10_clusters$Cluster-1
  rownames(tsne.data) = clustering$kmeans_10_clusters$Barcode
  if(DEBUG)cat(file=stderr(), "tsne.data: done\n")
  return(tsne.data)
})

