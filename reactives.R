# general reactive
# (same as global variables, but they are reactive to manipulation and are lazy, i.e. they only get executed when needed.)

# # default values
# # not working like this. need to update the values in each UI via an active update
#
# defaultValueSingleGeneReact <- reactive({
#   cat(stderr(), 'defaultValueSingleGene\n')
#   defaultValueSingleGene <<- input$defaultValueSingleGene
# })
#
# defaultValueMultiGenesReact <- reactive({
#   cat(stderr(), 'defaultValueMultiGenes\n')
#   defaultValueMultiGenes <<- input$defaultValueMultiGenes
# })



# reactive values  ------------------------------------------------------------------
# should only hold original data
# internal, should not be used by plug-ins
inputDataFunc <- function(inFile) {
  if (DEBUG) {
    cat(file = stderr(), "DEBUG:inputData\n")
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("loading", id = "inputDataFunc", duration = NULL)
  }
  start.time <- Sys.time()
  load(inFile$datapath)

  cat(stderr(), "Loaded")
  dataTables <- list()
  featuredata$Associated.Gene.Name <-
    toupper(featuredata$Associated.Gene.Name)
  featuredata <- featuredata[rownames(gbm), ]
  dataTables$featuredataOrg <- featuredata
  # dataTables$positiveCells <- NULL
  # dataTables$positiveCellsAll <- NULL

  # take only genes that are in all tables
  rnames <- rownames(featuredata)
  # rnames = rnames[rnames %in% rownames(log2cpm)]
  rnames <- rnames[rnames %in% rownames(gbm)]
  # rnames = rnames[rnames %in% rownames(gbm_log)]

  # cnames = colnames(log2cpm)
  # cnames =colnames(gbm)
  # cnames = cnames[cnames %in% colnames(gbm_log)]

  # dataTables$log2cpm <- log2cpm[rnames, cnames]
  dataTables$gbm <- gbm[rnames, ]
  # dataTables$gbm_log = gbm_log[rnames, cnames]
  dataTables$featuredata <- featuredata[rnames, ]

  if (is.null(gbm$barcode)) {
    showNotification("gbm doesn't contain barcode column")
    return(NULL)
  }
  # some checks
  if (sum(is.infinite(as.matrix(exprs(gbm)))) > 0) {
    showNotification("gbm contains infinite values")
    return(NULL)
  }
  if (sum(c("id", "symbol") %in% colnames(fData(gbm))) < 2) {
    showNotification("gbm - fData doesn't contain id and symbol columns", duration = NULL)
  }
  if (DEBUG) {
    cat(file = stderr(), "inputData: done\n")
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "inputDataFunc")
  }
  end.time <- Sys.time()
  cat(file = stderr(), paste("load data took: ", difftime(end.time, start.time, units = "min"), " min\n"))

  return(dataTables)
}

# internal, should not be used by plug-ins
inputData <- reactive({
  inFile <- input$file1
  if (is.null(inFile)) {
    if (DEBUG) {
      cat(file = stderr(), "inputData: NULL\n")
    }
    return(NULL)
  }
  return(inputDataFunc(inFile))
})

medianENSGfunc <- function(gbm) {
  geneC <- colSums(gbm > 0, na.rm = TRUE)
  return(median(t(geneC)))
}

medianENSG <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "medianENSG\n")
  }
  gbm <- gbm_matrix()
  if (is.null(gbm)) {
    if (DEBUG) {
      cat(file = stderr(), "medianENSG:NULL\n")
    }
    return(0)
  }
  if (ncol(gbm) <= 1 | nrow(gbm) < 1) {
    return(0)
  }
  retVal <- medianENSGfunc(gbm)
  if (DEBUG) {
    cat(file = stderr(), "medianENSG:done\n")
  }
  return(retVal)
})

medianUMIfunc <- function(gbm) {
  umiC <- colSums(gbm, na.rm = TRUE)
  return(median(t(umiC)))
}

medianUMI <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "medianUMI\n")
  }
  gbm <- gbm_matrix()
  if (is.null(gbm)) {
    if (DEBUG) {
      cat(file = stderr(), "medianUMI:NULL\n")
    }
    return(0)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/medianUMI.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/medianUMI.RData')
  retVal <- medianUMIfunc(gbm)
  if (DEBUG) {
    cat(file = stderr(), "medianUMI:done\n")
  }
  return(retVal)
})

# for now we don't have a way to specifically select cells
# we could cluster numbers or the like
# internal, should not be used by plug-ins
useCellsFunc <-
  function(dataTables,
             geneNames,
             rmCells,
             rmPattern,
             keepCells,
             cellKeepOnly) {
    if (DEBUG) {
      cat(file = stderr(), "useCells2\n")
    }
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/useCellsFunc.RData", list = c(ls()))
    }
    # load(file='~/scShinyHubDebug/useCellsFunc.Rdata')
    goodCols <- rep(TRUE, ncol(dataTables$gbm))
    gbm <- as.matrix(exprs(dataTables$gbm))
    #### start: cells with genes expressed
    # take only cells where these genes are expressed with at least one read
    genesin <- toupper(geneNames)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ",")
    genesin <- genesin[[1]]

    cellKeep <- toupper(keepCells)
    cellKeep <- gsub(" ", "", cellKeep, fixed = TRUE)
    cellKeep <- strsplit(cellKeep, ",")
    cellKeep <- cellKeep[[1]]

    cellKeepOnly <- toupper(cellKeepOnly)
    cellKeepOnly <- gsub(" ", "", cellKeepOnly, fixed = TRUE)
    cellKeepOnly <- strsplit(cellKeepOnly, ",")
    cellKeepOnly <- cellKeepOnly[[1]]

    # specifically remove cells
    if (nchar(rmCells) > 0) {
      cellsRM <- toupper(rmCells)
      cellsRM <- gsub(" ", "", cellsRM, fixed = TRUE)
      cellsRM <- strsplit(cellsRM, ",")
      cellsRM <- cellsRM[[1]]
      goodCols[which(toupper(colnames(dataTables$gbm)) %in% cellsRM)] <- FALSE
    }

    # remove cells by pattern
    if (nchar(rmPattern) > 0) {
      goodCols[grepl(rmPattern, colnames(dataTables$gbm))] <- FALSE
    }

    if (!length(cellKeep) == 0) {
      ids <- which(toupper(colnames(dataTables$gbm)) %in% cellKeep)
      goodCols[ids] <- TRUE
    }

    # genes that have to be expressed at least in one of them.
    selCols <- rep(FALSE, length(goodCols))
    if (!length(genesin) == 0) {
      ids <- which(toupper(dataTables$featuredata$Associated.Gene.Name) %in% genesin)
      if (length(ids) == 1) {
        selCols <- gbm[ids, ] > 0
      } else if (length(ids) == 0) {
        showNotification(
          "not enough cells, check gene names for min coverage",
          type = "warning",
          duration = NULL
        )
        return(NULL)
      } else {
        selCols <- colSums(gbm[ids, ]) > 0
      }
      goodCols <- goodCols & selCols
    }

    if (!length(cellKeepOnly) == 0) {
      goodCols[c(1:length(goodCols))] <- FALSE
      ids <- which(toupper(colnames(dataTables$gbm)) %in% cellKeepOnly)
      goodCols[ids] <- TRUE
    }

    #### end: cells with genes expressed

    if (DEBUG) {
      cat(file = stderr(), "useCells:done\n")
    }
    return(goodCols)
  }

# works on cells only
# internal, should not be used by plug-ins
useCells <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "useCells\n")
  }
  dataTables <- inputData()
  geneNames <- input$minExpGenes
  rmCells <- input$cellsFiltersOut
  rmPattern <- input$cellPatternRM
  keepCells <- input$cellKeep
  cellKeepOnly <- input$cellKeepOnly
  # useGenes = isolate(useGenes())
  if (!exists("dataTables") || is.null(dataTables)) {
    if (DEBUG) {
      cat(file = stderr(), "useCells:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("which cells to use", id = "useCells", duration = NULL)
  }
  retVal <- useCellsFunc(
    dataTables,
    geneNames,
    rmCells,
    rmPattern,
    keepCells,
    cellKeepOnly
  )
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "useCells")
  }
  return(retVal)
})

# TODO: check that it is ok that we use dataTables directly and not useGenes()
featureDataReact <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "featureData\n")
  }
  dataTables <- inputData()
  gbm <- gbm()
  if (!exists("dataTables") | is.null(dataTables) | is.null(gbm)) {
    if (DEBUG) {
      cat(file = stderr(), "featureData:NULL\n")
    }
    return(NULL)
  }
  useGenes <- rownames(dataTables$featuredata) %in% rownames(gbm)
  if (DEBUG) {
    cat(file = stderr(), "featureData:done\n")
  }
  return(dataTables$featuredata[useGenes, ])
})

useGenesFunc <-
  function(dataTables,
             ipIDs, # regular expression of genes to be removed
             geneListSelection,
             genesKeep,
             geneLists) {
    gList <- geneLists # global variable, assigning it locally ensures that it will be saved
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/useGenesFunc.Rdata", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file='~/scShinyHubDebug/useGenesFunc.Rdata')
    # regular expression with gene names to be removed
    if (nchar(ipIDs) > 0) {
      keepIDs <- !grepl(ipIDs, dataTables$featuredata$Associated.Gene.Name)
    } else {
      keepIDs <- rep(TRUE, nrow(dataTables$gbm))
    }
    genesKeep <- toupper(genesKeep)
    genesKeep <- gsub(" ", "", genesKeep, fixed = TRUE)
    genesKeep <- strsplit(genesKeep, ",")
    genesKeep <- genesKeep[[1]]
    keepGeneIds <- which(dataTables$featuredata$Associated.Gene.Name %in% genesKeep)

    # dataTables$featuredata$Associated.Gene.Name[keepIDs]
    # gene groups to be included
    if (!is.null(geneListSelection)) {
      selectedgeneList <- get_selected(geneListSelection)
      if (length(selectedgeneList) > 0) {
        selGenes <- c()
        for (sIdx in 1:length(selectedgeneList)) {
          print(sIdx)
          att <- attr(selectedgeneList[[sIdx]], "ancestry")
          if (length(att) > 0) {
            selGenes <- c(selGenes, gList[[att]][[selectedgeneList[[sIdx]]]])
          }
        }
        selGenes <- unique(selGenes)
        keepIDs <- (rownames(dataTables$gbm) %in% selGenes) & keepIDs
      }
    }

    # # overall gene expression Min
    # if(!is.null(minGene)){
    #   selGenes = rowSums(as.matrix(exprs(dataTables$gbm[,useCells]))) >=minGene
    #   keepIDs = keepIDs & selGenes
    # }
    keepIDs[keepGeneIds] <- TRUE
    return(keepIDs)
  }

# collects information from all places where genes being removed or specified
useGenes <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "useGenes\n")
  }
  dataTables <- inputData()
  # useCells = useCells()
  # minGene <- input$minGenesGS
  ipIDs <- input$selectIds # regular expression of genes to be removed
  genesKeep <- input$genesKeep
  geneListSelection <- input$geneListSelection

  if (!exists("dataTables") |
    is.null(dataTables) |
    length(dataTables$featuredata$Associated.Gene.Name) == 0) {
    if (DEBUG) {
      cat(file = stderr(), "useGenes: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("which genes to use", id = "useGenes", duration = NULL)
  }
  retVal <- useGenesFunc(dataTables, ipIDs, geneListSelection, genesKeep, geneLists)
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "useGenes")
  }
  if (DEBUG) {
    cat(file = stderr(), "useGenes: done\n")
  }
  return(retVal)
})


# recursive removal
# we need to have a first test to say if we can go

# will be called recursively to ensure that nothing changes when cells/genes are changing.
gbmFunc <-
  function(gbmOrg,
             useCells,
             useGenes,
             minGene,
             minG,
             maxG) {
    # save(file="~/scShinyHubDebug/gbmFunc.RData", list=ls())
    # load(file="~/scShinyHubDebug/gbmFunc.RData")
    if (DEBUG) {
      cat(file = stderr(), "gbmFunc\n")
    }
    # if(DEBUG)cat(file=stderr(), paste("col:",ncol(gbmOrg),"\n"))
    # if(DEBUG)cat(file=stderr(), paste("l useCells:",length(useCells),"\n"))
    # if(DEBUG)cat(file=stderr(), paste("row:",nrow(gbmOrg),"\n"))
    # if(DEBUG)cat(file=stderr(), paste("l useGenes:",length(useGenes),"\n"))
    # gbmOrg, useCells, useGenes cannot be NULL

    # change names to be hopefully a bit more clear
    changed <- FALSE # trace if something changed
    keepGenes <- useGenes
    keepCells <- useCells
    gbm <- as.matrix(exprs(gbmOrg))

    # overall gene expression Min
    if (!is.null(minGene)) {
      selGenes <- rowSums(gbm[, keepCells]) >= minGene
      selGenes <- keepGenes & selGenes
      if (!all(selGenes == keepGenes)) {
        keepGenes <- selGenes
        changed <- TRUE
      }
    }

    # min reads per cell
    if (!is.null(minG)) {
      selCols <- colSums(gbm[keepGenes, ], na.rm = FALSE) > minG
      selCols[is.na(selCols)] <- FALSE
      selCols <- keepCells & selCols
      if (!all(selCols == keepCells)) {
        keepCells <- selCols
        changed <- TRUE
      }
    }

    # max reads per cell
    if (!is.null(maxG)) {
      selCols <- colSums(gbm[keepGenes, ], na.rm = FALSE) <= maxG
      selCols[is.na(selCols)] <- FALSE
      selCols <- selCols & keepCells
      if (!all(selCols == keepCells)) {
        changed <- TRUE
        keepCells <- selCols
      }
    }

    if (sum(keepCells) == 0) {
      showNotification("not enough cells left",
        type = "warning",
        duration = NULL
      )
      return(NULL)
    }
    if (sum(keepGenes) == 0) {
      showNotification("not enough genes left",
        type = "warning",
        duration = NULL
      )
      return(NULL)
    }

    # if something changed, check that it doesn't change again
    gbmNew <- gbmOrg[keepGenes, keepCells]
    if (changed) {
      gbmNew <- gbmFunc(gbmOrg[keepGenes, keepCells], useCells[keepCells], useGenes[keepGenes], minGene, minG, maxG)
      if (is.null(gbmNew)) {
        return(NULL)
      }
    }

    pD <- pData(gbmNew)
    for (colN in colnames(pD)) {
      if (colN == "barcode") next()
      pD[, colN] <- factor(as.character(pD[, colN]))
    }
    pData(gbmNew) <- pD

    return(gbmNew)
  }


# apply filters that depend on genes & cells
# it is here that useCells and useGenes are combined and applied to select for
gbm <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gbm\n")
  }
  dataTables <- inputData()
  useCells <- useCells()
  useGenes <- useGenes()
  minGene <- input$minGenesGS # min number of reads per gene
  minG <- input$minGenes # min number of reads per cell
  maxG <- input$maxGenes # max number of reads per cell
  if (!exists("dataTables") |
    is.null(dataTables) | is.null(useGenes) | is.null(useCells)) {
    if (DEBUG) {
      cat(file = stderr(), "gbm: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gbm", id = "gbm", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gbm.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gbm.RData")

  retVal <- gbmFunc(
    gbmOrg = dataTables$gbm,
    useCells = useCells,
    useGenes = useGenes,
    minGene = minGene,
    minG = minG,
    maxG = maxG
  )
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "gbm")
  }
  if (DEBUG) {
    cat(file = stderr(), "gbm:DONE\n")
  }
  return(retVal)
})

# takes a lot of memory and should be avoided.
gbm_matrix <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gbm_matrix\n")
  }
  gbm <- gbm()
  if (is.null(gbm)) {
    if (DEBUG) {
      cat(file = stderr(), "gbm_matrix:NULL\n")
    }
    return(NULL)
  }
  if (ncol(gbm) <= 1 | nrow(gbm) < 1) {
    return(NULL)
  }
  retVal <- as.matrix(exprs(gbm))
  if (DEBUG) {
    cat(file = stderr(), "gbm_matrix:done\n")
  }
  return(retVal)
})


rawNormalization <- reactive({
  gbm <- gbm()
  if (DEBUG) {
    cat(file = stderr(), "rawNormalization\n")
  }
  return(gbm)
})

# individual values
gbm_log <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gbm_log\n")
  }
  # dataTables = inputData()
  # useCells = useCells()
  # useGenes = useGenes()
  gbm <- gbm()
  normMethod <- input$normalizationRadioButton

  if (is.null(gbm)) {
    if (DEBUG) {
      cat(file = stderr(), "gbm_log:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Normalizing data", id = "gbm_log", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gbm_log.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gbm_log.RData")

  gbm_log <- do.call(normMethod, args = list())

  # gbm rownames are ENSG numbers
  # dataTables$gbm_log[useGenes, useCells]
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "gbm_log")
  }
  if (DEBUG) {
    cat(file = stderr(), "gbm_log:Done\n")
  }
  return(gbm_log)
})

gbmLogMatrix <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gbmLogMatrix\n")
  }
  # dataTables = inputData()
  # useCells = useCells()
  # useGenes = useGenes()
  gbmLog <- gbm_log()
  if (is.null(gbmLog)) {
    if (DEBUG) {
      cat(file = stderr(), "gbmLogMatrix:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Calculating gbmLogmatrix",
      id = "gbmLogMatrix",
      duration = NULL
    )
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gbmLogMatrix.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gbmLogMatrix.RData")

  retVal <- as.data.frame(as.matrix(exprs(gbmLog)))
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "gbmLogMatrix")
  }
  if (DEBUG) {
    cat(file = stderr(), "gbmLogMatrix:Done\n")
  }
  return(retVal)
})

# gbmLog matrix with symbol as first column
# TODO
# we should probably just rename the rows and then have an option to tableSelectionServer that shows (or not) rownames
gbmLogMatrixDisplay <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gbmLogMatrixDisplay\n")
  }
  # dataTables = inputData()
  # useCells = useCells()
  # useGenes = useGenes()
  gbmLog <- gbm_log()
  if (is.null(gbmLog)) {
    if (DEBUG) {
      cat(file = stderr(), "gbmLogMatrixDisplay:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Calculating gbmLogmatrix",
      id = "gbmLogMatrixDisplay",
      duration = NULL
    )
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gbmLogMatrixDisplay.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gbmLogMatrixDisplay.RData")

  retVal <- as.data.frame(as.matrix(exprs(gbmLog)))
  rownames(retVal) <- make.names(fData(gbmLog)$symbol, unique = TRUE)

  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "gbmLogMatrixDisplay")
  }
  if (DEBUG) {
    cat(file = stderr(), "gbmLogMatrixDisplay:Done\n")
  }
  return(retVal)
})

pcaFunc <- function(gbm_log) {
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/pcaFunc.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/pcaFunc.RData")
  pca <- tryCatch({
    run_pca(gbm_log)
  },
  error = function(e) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "Problem with PCA, probably not enough cells?",
        type = "warning",
        duration = NULL
      )
    }
    return(NULL)
  }
  )
  return(pca)
}

pca <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "pca\n")
  }
  gbm_log <- gbm_log()
  if (is.null(gbm_log)) {
    if (DEBUG) {
      cat(file = stderr(), "pca:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("pca", id = "pca", duration = NULL)
  }
  retVal <- pcaFunc(gbm_log)
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "pca")
  }
  if (DEBUG) {
    cat(file = stderr(), "pca:done\n")
  }

  return(retVal)
})


kmClusteringFunc <- function(pca, seed, kNr = 10) {
  clustering <- list()

  # kNr = 10
  # for (kNr in 2:kNr) {
  set.seed(seed = seed)
  km <- run_kmeans_clustering(pca, k = kNr)
  clustering[[paste0("kmeans_", kNr, "_clusters")]] <- data.frame(
    "Barcode" = rownames(data.frame(km$cluster)),
    "Cluster" = km$cluster
  )
  # }
  return(clustering)
}

kmClustering <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "kmClustering\n")
  }
  pca <- pca()
  seed <- input$seed
  kNr <- input$kNr
  # kNr = 10
  if (is.null(pca)) {
    if (DEBUG) {
      cat(file = stderr(), "kmClustering:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/kmClustering.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/kmClustering.RData")
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("kmClustering", id = "kmClustering", duration = NULL)
  }

  if (is.null(seed)) {
    seed <- 1
  }
  retVal <- kmClusteringFunc(pca, seed = seed, kNr = kNr)
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "kmClustering")
  }
  if (DEBUG) {
    cat(file = stderr(), "kmClustering:done\n")
  }

  return(retVal)
})


# TODO separate  function from reactive
# projections
# each column is of length of number of cells
# if factor than it is categorical and can be cluster number of sample etc
# if numeric can be projection

# projections is a reactive and cannot be used in reports. Reports have to organize
# themselves as it is done here with tsne.data.

projections <- reactive({
  # gbm is the fundamental variable with the raw data, which is available after loading
  # data. Here we ensure that everything is loaded and all varialbles are set by waiting
  # input data being loaded
  gbm <- gbm()
  pca <- pca()
  if (!exists("gbm") | is.null(gbm) | !exists("pca") | is.null(pca)) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  projections <- data.frame(pca$x[, c(1, 2, 3)])
  if (DEBUG) {
    cat(file = stderr(), "projections\n")
  }

  # phenotypic data/ annotations of cells can already be included in the gbm object. We collect this information, but only for variable that hold information
  # i.e. length(levels) > 1 & < number of rows
  pd <- pData(gbm)
  if (ncol(pd) < 2) {
    cat(file = stderr(), "phenoData for gbm has less than 2 columns\n")
    return(NULL)
  }


  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/projections.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/projections.RData")
  withProgress(message = "Performing projections", value = 0, {
    n <- length(projectionFunctions)
    iter <- 1
    for (proj in projectionFunctions) {
      incProgress(1 / n, detail = paste("Creating ", proj[1]))
      if (DEBUG) cat(file = stderr(), paste("calculation projection:  ", proj[1], "\n"))
      assign("tmp", eval(parse(text = paste0(proj[2], "()"))))
      if (DEBUGSAVE) {
        save(file = paste0("~/scShinyHubDebug/projections.", iter, ".RData"), list = c("tmp"))
        iter <- iter + 1
      }
      # load(file="~/scShinyHubDebug/projections.1.RData")
      if (class(tmp) == "data.frame") {
        cn <- make.names(c(colnames(projections), colnames(tmp)))
      } else {
        cn <- make.names(c(colnames(projections), make.names(proj[1])))
      }
      if (length(tmp) == 0) {
        #   showNotification(
        #     paste("warning: ", proj[1], "didn't produce a result"),
        #     type = "warning",
        #     duration = NULL
        #   )
        next()
      }
      if (ncol(projections) == 0) {
        # never happening because we set pca first
        projections <- data.frame(tmp = tmp)
      } else {
        if (nrow(projections) == length(tmp)) {
          projections <- cbind(projections, tmp)
        } else if (nrow(projections) == nrow(tmp)) {
          projections <- cbind(projections, tmp)
        } else {
          stop("error: ", proj[1], "didn't produce a result")
        }
      }

      colnames(projections) <- cn
      observe(proj[2], quoted = TRUE)
    }
  })
  # add a column for gene specific information that will be filled/updated on demand
  projections$UmiCountPerGenes <- 0
  return(projections)
})

groupNames <- reactiveValues(
  namesDF = data.frame()
  # {
  # if(DEBUG)
  #   cat(file = stderr(), "groupNames\n")
  # projections = projections()
  # if(is.null(projections)){
  #   return(NULL)
  # }
  # if (DEBUGSAVE)
  #   save(file = "~/scShinyHubDebug/groupNames.RData", list = c(ls(),ls(envir = globalenv())))
  # # load(file="~/scShinyHubDebug/groupNames.RData")
  # retVal = c()
  # for(cn in colnames(projections)){
  #   if(class(projections[,cn])=="logical"){
  #     retVal = c(retVal, cn)
  #   }
  # }
  # return(retVal)
  # }
)

initializeGroupNames <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "initializeGroupNames\n")
  }
  gbm <- gbm()
  if (is.null(gbm)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/initializeGroupNames.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/initializeGroupNames.RData")
  isolate({
    df <- data.frame(all = rep(TRUE, ncol(gbm)), none = rep(FALSE, ncol(gbm)))
    rownames(df) <- colnames(gbm)
    cat(file = stderr(), "initializeGroupNames2\n")
    groupNames[["namesDF"]] <- df
    cat(file = stderr(), "initializeGroupNames3\n")
  })
})


dbCluster <- reactive({
  kNr <- input$kNr
  # kNr = 10
  if (DEBUG) {
    cat(file = stderr(), "dbCluster\n")
  }
  clustering <- kmClustering()
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/dbCluster.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/dbCluster.RData")

  if (is.null(clustering)) {
    if (DEBUG) {
      cat(file = stderr(), "dbCluster: NULL\n")
    }
    return(NULL)
  }

  dbCluster <- factor(clustering[[paste0("kmeans_", kNr, "_clusters")]]$Cluster - 1)

  return(dbCluster)
})

# sample --------
sample <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "sample\n")
  }
  gbm <- gbm()
  if (is.null(gbm)) {
    if (DEBUG) {
      cat(file = stderr(), "sample: NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/sample.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/sample.RData")
  # samp <- gsub(".*-(.*)", "\\1", colnames(gbm))
  # if (length(levels(as.factor(samp))) > 1) {
  #   sample <- samp
  # } else {
  #   sample <- rep("1", ncol(gbm))
  # }
  pd <- pData(gbm)
  retVal <- NULL
  for (pdColName in colnames(pd)) {
    if (length(levels(factor(pd[, pdColName]))) < 100) {
      if (is.null(retVal)) {
        retVal <- data.frame(pd[, pdColName])
        colnames(retVal) <- pdColName
      }
      retVal[, pdColName] <- factor(as.character(pd[, pdColName]))
    }
  }
  retVal
  # return(sample)
})

# geneCount --------
geneCount <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "geneCount\n")
  }
  gbm <- gbm()
  if (is.null(gbm)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/geneCount.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/geneCount.RData")
  retVal <- colSums(as.matrix(exprs(gbm)) > 0)
  return(retVal)
})

umiCount <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "umiCount\n")
  }
  gbm <- gbm()
  if (is.null(gbm)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/umiCount.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/umiCount.RData")
  retVal <- colSums(as.matrix(exprs(gbm)))
  return(retVal)
})


sampleInfoFunc <- function(gbm) {
  # gsub(".*-(.*)", "\\1", gbm$barcode)
  pData(gbm)$sampleNames
}


# sampleInfo -------
# sample information
sampleInfo <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "sampleInfo\n")
  }
  gbm <- gbm()
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/sampleInfo.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/sampleInfo.RData")
  if (!exists("gbm")) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }

  ret <- sampleInfoFunc(gbm)
  if (DEBUG) {
    cat(file = stderr(), "sampleInfo: done\n")
  }
  return(ret)
})


# table of input cells with sample information
inputSample <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "inputSample\n")
  }
  dataTables <- inputData()
  # sampInf = sampleInfo() # not working as this is relative to the activated cells
  if (is.null(dataTables)) {
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputSample", id = "inputSample", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/inputSample.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/inputSample.RData')
  sampInf <- gsub(".*-(.*)", "\\1", dataTables$gbm$barcode)
  cellIds <- data.frame(
    cellName = colnames(dataTables$gbm),
    sample = sampInf,
    ngenes = colSums(as.matrix(exprs(dataTables$gbm)))
  )

  if (DEBUG) {
    cat(file = stderr(), "inputSample: done\n")
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "inputSample")
  }
  if (dim(cellIds)[1] > 1) {
    return(cellIds)
  } else {
    return(NULL)
  }
})

updateMemUse <- reactiveValues(
  update = 1
)


getMemoryUsed <- reactive({
  require(pryr)
  if (DEBUG) {
    cat(file = stderr(), "getMemoryUsed\n")
  }
  umu <- updateMemUse$update
  paste(utils:::format.object_size(mem_used(), "auto"), umu)
})

# used in coExpression, subclusterAnalysis, moduleServer, generalQC, DataExploration
# TODO change to gbm_log everywhere and remove
log2cpm <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "log2cpm\n")
  }
  gbmLog <- gbm_log()
  if (is.null(gbmLog)) {
    if (DEBUG) {
      cat(file = stderr(), "log2cpm: NULL\n")
    }
    return(NULL)
  }
  log2cpm <- as.data.frame(as.matrix(exprs(gbmLog)))

  return(log2cpm)
})

#### plot2Dprojection ----------------
# used in moduleServer and reports
plot2Dprojection <- function(gbm_log, gbm, projections, g_id, featureData, geneNames, dimX, dimY, clId, grpN, legend.position) {
  geneid <- geneName2Index(g_id, featureData)
  grpNs <- groupNames$namesDF

  # if (length(geneid) == 1) {
  #   expression <- exprs(gbm_log)[geneid, ,drop=FALSE]
  # } else {
    expression <- Matrix::colSums(exprs(gbm_log)[geneid, ,drop=FALSE])
  # }
  validate(need(is.na(sum(expression)) != TRUE, ""))
  
  projections <- updateProjectionsWithUmiCount(dimY, dimX, geneNames, featureData, gbm, projections)


  projections <- cbind(projections, t(expression))
  names(projections)[ncol(projections)] <- "exprs"

  if (DEBUG) {
    cat(file = stderr(), paste("output$dge_plot1:---", clId, "---\n"))
  }
  subsetData <- subset(projections, dbCluster %in% clId)
  # subsetData$dbCluster = factor(subsetData$dbCluster)
  # if there are more than 18 samples ggplot cannot handle different shapes and we ignore the
  # sample information
  if (length(as.numeric(as.factor(subsetData$sample))) > 18) {
    subsetData$shape <- 1
  } else {
    subsetData$shape <- as.numeric(as.factor(subsetData$sample))
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/clusterPlot.RData", list = c(ls(), "legend.position", ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/clusterPlot.RData")
  if (nrow(subsetData) == 0) return(NULL)
  p1 <-
    ggplot(
      subsetData,
      aes_string(x = dimX, y = dimY)
    ) +
    geom_point(aes_string(shape = "shape", size = 2, color = "exprs"), show.legend = TRUE) +
    scale_shape_identity() +
    geom_point(
      shape = 1,
      size = 4,
      aes(colour = as.numeric(dbCluster))
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        size = 12,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      legend.position = legend.position
    ) +
    ggtitle(paste(toupper(g_id), clId, sep = "-Cluster", collapse = " ")) +
    scale_fill_continuous()
  selectedCells <- NULL
  if (length(grpN) > 0) {
    if (length(grpNs[rownames(subsetData), grpN]) > 0 & sum(grpNs[rownames(subsetData), grpN]) > 0) {
      grpNSub <- grpNs[rownames(subsetData), ]
      selectedCells <- rownames(grpNSub[grpNSub[, grpN], ])
    }
  }
  if (!is.null(selectedCells)) {
    shape <- rep("a", nrow(subsetData))
    selRows <- which(rownames(subsetData) %in% selectedCells)
    shape[selRows] <- "b"
    p1 <- p1 + geom_point(data = subsetData[selRows, ], mapping = aes(shape = shape, size = 4), colour = "red")
  }
  p1
}
