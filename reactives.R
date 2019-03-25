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

# library(cellrangerRkit)
inputFileStats <- reactiveValues(
  stats = NULL
)


sampleCols <- reactiveValues(
  colPal = c("1" = colorRampPalette(brewer.pal(
    n = 6, name =
      "RdYlBu"
  ))(1))
)
# reactive values  ------------------------------------------------------------------
# should only hold original data
# internal, should not be used by plug-ins
inputDataFunc <- function(inFile) {
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "inputDataFunc")
    }
  )
  start.time <- Sys.time()
  if (DEBUG) {
    cat(file = stderr(), "DEBUG:inputData\n")
    save(file = "inputDataFunc.RData", list = c("inFile"))
  }
  # load("inputDataFunc.RData")
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("loading", id = "inputDataFunc", duration = NULL)
  }
  start.time <- Sys.time()

  # load("test.RData")

  stats <- tibble(.rows = length(inFile$datapath))
  stats$names <- inFile$name
  stats$nFeatures <- 0
  stats$nCells <- 0

  #
  cat(file = stderr(), paste("reading", inFile$name[1], "\n"))
  fp <- inFile$datapath[1]
  # fp ="scEx.Rds"
  fpLs <- load(fp)
  scExFound <- FALSE
  for (varName in fpLs) {
    if ("SingleCellExperiment" %in% class(get(varName))) {
      scEx <- get(varName)
      scExFound <- TRUE
    }
  }
  if (!scExFound) {
    return(NULL)
  }
  # fd <- featuredata
  fdAll <- rowData(scEx)
  pdAll <- colData(scEx)
  exAll <- assays(scEx)[["counts"]]
  stats[1, "nFeatures"] <- nrow(fdAll)
  stats[1, "nCells"] <- nrow(pdAll)

  if (length(inFile$datapath) > 1) {
    for (fpIdx in 2:length(inFile$datapath)) {
      cat(file = stderr(), paste("reading", inFile$name[fpIdx], "\n"))
      fp <- inFile$datapath[fpIdx]
      fpLs <- load(fp)
      if (!"scEx" %in% fpLs) {
        next()
      }
      fdIdx <- intersect(rownames(fdAll), rownames(rowData(scEx)))
      # if (length(fdIdx) != nrow(fd)) {
      #   cat(file = stderr(), "Houston, there is a problem with the features\n")
      # }
      # fd <- featuredata[fdIdx, ]
      fdAll <- fdAll[fdIdx, ]
      pd1 <- colData(scEx)
      ex1 <- assays(scEx)[["counts"]][fdIdx, ]
      if (sum(rownames(pdAll) %in% rownames(pd1)) > 0) {
        cat(file = stderr(), "Houston, there are cells with the same name\n")
        rownames(pd1) <- paste0(rownames(pd1), "_", fpIdx)
        pdAll <- rbind(pdAll, pd1)
        colnames(ex1) <- rownames(pd1)
      } else {
        pdAll <- rbind(pdAll, pd1)
      }
      stats[fpIdx, "nFeatures"] <- nrow(fd)
      stats[fpIdx, "nCells"] <- nrow(pd1)

      exAll <- Matrix::cbind2(exAll[fdIdx, ], ex1)
    }
  }
  exAll <- as(exAll, "dgTMatrix")
  scEx <- SingleCellExperiment(
    assay = list(counts = exAll),
    colData = pdAll,
    rowData = fdAll
  )

  cat(stderr(), "Loaded")
  dataTables <- list()
  featuredata <- rowData(scEx)
  dataTables$featuredataOrg <- rowData(scEx)
  dataTables$scEx <- scEx
  dataTables$featuredata <- featuredata

  if (is.null(scEx$barcode)) {
    showNotification("scEx doesn't contain barcode column", type = "error")
    return(NULL)
  }
  # some checks

  if (sum(is.infinite(assays(scEx)[["counts"]])) > 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("scEx contains infinite values",
        type = "error"
      )
    }
    return(NULL)
  }

  if ("sampleNames" %in% names(colData(scEx))) {
    sampNames <- levels(colData(scEx)$sampleNames)
    isolate({
      sampleCols$colPal <- colorRampPalette(brewer.pal(
        n = 6, name =
          "PRGn"
      ))(length(sampNames))
      names(sampleCols$colPal) <- sampNames
    })
  } else {
    showNotification("scEx - colData doesn't contain sampleNames",
      duration = NULL, type = "error"
    )
  }

  if (sum(c("id", "symbol") %in% colnames(rowData(scEx))) < 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("scEx - rowData doesn't contain id and/or symbol columns",
        duration = NULL, type = "error"
      )
    }
  }

  if (!sum(c("Associated.Gene.Name", "Gene.Biotype", "Description") %in% colnames(featuredata)) == 3) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("featuredata - one of is missing: Associated.Gene.Name, Gene.Biotype, Description)",
        duration = NULL, type = "error"
      )
    }
    if (!"Gene.Biotype" %in% colnames(featuredata)) {
      featuredata$"Gene.Biotype" <- "not given"
    }
    if (!"Description" %in% colnames(featuredata)) {
      featuredata$"Description" <- "not given"
    }
    dataTables$featuredata <- featuredata
  }
  # if (is.null(rowData(dataTables$scEx)$symbol)){
  #
  # }
  if (DEBUG) {
    cat(file = stderr(), "inputData: done\n")
  }
  end.time <- Sys.time()
  cat(file = stderr(), paste("===load data:done", difftime(end.time, start.time, units = "min"), " min\n"))

  inputFileStats$stats <- stats
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

medianENSGfunc <- function(scEx) {
  geneC <- Matrix::colSums(scEx > 0, na.rm = TRUE)
  return(median(t(geneC)))
}

medianENSG <- reactive({
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "medianENSG\n")
  }
  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "medianENSG:NULL\n")
    }
    return(0)
  }
  scEx <- assays(scEx)[["counts"]]
  if (ncol(scEx) <= 1 | nrow(scEx) < 1) {
    return(0)
  }
  retVal <- medianENSGfunc(scEx)
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), paste("===medianENSG:done", difftime(end.time, start.time, units = "min"), "\n"))
  }
  return(retVal)
})

medianUMIfunc <- function(scEx) {
  umiC <- Matrix::colSums(scEx, na.rm = TRUE)
  return(median(t(umiC)))
}

medianUMI <- reactive({
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "medianUMI\n")
  }
  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "medianUMI:NULL\n")
    }
    return(0)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/medianUMI.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/medianUMI.RData')
  scEx <- assays(scEx)[["counts"]]
  retVal <- medianUMIfunc(scEx)
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===medianUMI:done\n", difftime(end.time, start.time, units = "min"), "\n")
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
    goodCols <- rep(TRUE, ncol(dataTables$scEx))
    scEx <- assays(dataTables$scEx)[[1]]
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
      goodCols[which(toupper(colnames(dataTables$scEx)) %in% cellsRM)] <- FALSE
    }

    # remove cells by pattern
    if (nchar(rmPattern) > 0) {
      goodCols[grepl(rmPattern, colnames(dataTables$scEx))] <- FALSE
    }

    if (!length(cellKeep) == 0) {
      ids <- which(toupper(colnames(dataTables$scEx)) %in% cellKeep)
      goodCols[ids] <- TRUE
    }

    # genes that have to be expressed at least in one of them.
    selCols <- rep(FALSE, length(goodCols))
    if (!length(genesin) == 0) {
      ids <- which(toupper(dataTables$featuredata$Associated.Gene.Name) %in% genesin)
      if (length(ids) == 1) {
        selCols <- scEx[ids, ] > 0
      } else if (length(ids) == 0) {
        showNotification(
          "not enough cells, check gene names for min coverage",
          type = "warning",
          duration = NULL
        )
        return(NULL)
      } else {
        selCols <- Matrix::colSums(scEx[ids, ]) > 0
      }
      goodCols <- goodCols & selCols
    }

    if (!length(cellKeepOnly) == 0) {
      goodCols[c(1:length(goodCols))] <- FALSE
      ids <- which(toupper(colnames(dataTables$scEx)) %in% cellKeepOnly)
      goodCols[ids] <- TRUE
    }

    #### end: cells with genes expressed

    return(goodCols)
  }

# works on cells only
# internal, should not be used by plug-ins
useCells <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "useCells")
    }
  )
  start.time <- Sys.time()
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
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===useCells:done", difftime(end.time, start.time, units = "min"), "\n")
  }

  return(retVal)
})

# TODO: check that it is ok that we use dataTables directly and not useGenes()
featureDataReact <- reactive({
  start.time <- Sys.time()
  if (DEBUG) {
    cat(file = stderr(), "featureData\n")
  }
  dataTables <- inputData()
  scEx <- scEx()
  if (!exists("dataTables") | is.null(dataTables) | is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "featureData:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/featureDataReact.Rdata", list = c(ls(), ls(envir = globalenv())))
  }

  useGenes <- rownames(dataTables$featuredata) %in% rownames(scEx)
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===featureData:done", difftime(end.time, start.time, units = "min"), "\n")
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
      keepIDs <- rep(TRUE, nrow(dataTables$scEx))
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
        keepIDs <- (rownames(dataTables$scEx) %in% selGenes) & keepIDs
      }
    }

    keepIDs[keepGeneIds] <- TRUE
    return(keepIDs)
  }

# before gene filtering
beforeFilterCounts <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "beforeFilterCounts")
    }
  )
  dataTables <- inputData()
  ipIDs <- input$selectIds # regular expression of genes to be removed
  if (!exists("dataTables") |
    is.null(dataTables) |
    length(dataTables$featuredata$Associated.Gene.Name) == 0) {
    if (DEBUG) {
      cat(file = stderr(), "beforeFilterCounts: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("beforeFilterCounts", id = "beforeFilterCounts", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/beforeFilterCounts.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/beforeFilterCounts.RData")

  geneIDs <- NULL
  if (nchar(ipIDs) > 0) {
    geneIDs <- grepl(ipIDs, dataTables$featuredata$Associated.Gene.Name)
  }
  if (is.null(geneIDs)) {
    return(rep(0, nrow(dataTables$featuredata)))
  }
  return(Matrix::colSums(assays(dataTables$scEx)[["counts"]][geneIDs, ]))
})

# collects information from all places where genes being removed or specified
useGenes <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "useGenes")
    }
  )
  start.time <- Sys.time()
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
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===useGenes: done", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(retVal)
})


# recursive removal
# we need to have a first test to say if we can go

# will be called recursively to ensure that nothing changes when cells/genes are changing.
scExFunc <-
  function(scExOrg,
             useCells,
             useGenes,
             minGene,
             minG,
             maxG) {
    save(file = "~/scShinyHubDebug/scExFunc.RData", list = ls())
    # load(file="~/scShinyHubDebug/scExFunc.RData")
    if (DEBUG) {
      cat(file = stderr(), "scExFunc\n")
    }
    # if(DEBUG)cat(file=stderr(), paste("col:",ncol(scExOrg),"\n"))
    # if(DEBUG)cat(file=stderr(), paste("l useCells:",length(useCells),"\n"))
    # if(DEBUG)cat(file=stderr(), paste("row:",nrow(scExOrg),"\n"))
    # if(DEBUG)cat(file=stderr(), paste("l useGenes:",length(useGenes),"\n"))
    # scExOrg, useCells, useGenes cannot be NULL

    # change names to be hopefully a bit more clear
    changed <- FALSE # trace if something changed
    keepGenes <- useGenes
    keepCells <- useCells
    scEx <- assays(scExOrg)[[1]]

    # overall gene expression Min
    if (!is.null(minGene)) {
      selGenes <- Matrix::rowSums(scEx[, keepCells]) >= minGene
      selGenes <- keepGenes & selGenes
      if (!all(selGenes == keepGenes)) {
        keepGenes <- selGenes
        changed <- TRUE
      }
    }

    # min reads per cell
    if (!is.null(minG)) {
      selCols <- Matrix::colSums(scEx[keepGenes, ], na.rm = FALSE) > minG
      selCols[is.na(selCols)] <- FALSE
      selCols <- keepCells & selCols
      if (!all(selCols == keepCells)) {
        keepCells <- selCols
        changed <- TRUE
      }
    }

    # max reads per cell
    if (!is.null(maxG)) {
      selCols <- Matrix::colSums(scEx[keepGenes, ], na.rm = FALSE) <= maxG
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
    scExNew <- scExOrg[keepGenes, keepCells]
    if (changed) {
      scExNew <- scExFunc(scExOrg[keepGenes, keepCells], useCells[keepCells], useGenes[keepGenes], minGene, minG, maxG)
      if (is.null(scExNew)) {
        return(NULL)
      }
    }

    pD <- colData(scExNew)
    for (colN in colnames(pD)) {
      if (colN == "barcode") next()
      if (class(pD[, colN]) %in% c("character")) {
        pD[, colN] <- factor(as.character(pD[, colN]))
      }
    }
    colData(scExNew) <- pD

    return(scExNew)
  }


# apply filters that depend on genes & cells
# it is here that useCells and useGenes are combined and applied to select for
scEx <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scEx")
    }
  )
  start.time <- Sys.time()
  if (DEBUG) {
    cat(file = stderr(), "scEx\n")
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
      cat(file = stderr(), "scEx: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scEx", id = "scEx", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scEx.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/scEx.RData")

  retVal <- scExFunc(
    scExOrg = dataTables$scEx,
    useCells = useCells,
    useGenes = useGenes,
    minGene = minGene,
    minG = minG,
    maxG = maxG
  )
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===scEx:DONE", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(retVal)
})

rawNormalization <- reactive({
  scEx <- scEx()
  if (DEBUG) {
    cat(file = stderr(), "rawNormalization\n")
  }
  return(scEx)
})

# individual values
scEx_log <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scEx_log")
    }
  )
  start.time <- Sys.time()
  if (DEBUG) {
    cat(file = stderr(), "scEx_log\n")
  }
  # dataTables = inputData()
  # useCells = useCells()
  # useGenes = useGenes()
  scEx <- scEx()
  normMethod <- input$normalizationRadioButton

  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "scEx_log:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Normalizing data", id = "scEx_log", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scEx_log.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/scEx_log.RData")

  scEx_log <- do.call(normMethod, args = list())

  # scEx rownames are ENSG numbers
  # dataTables$scEx_log[useGenes, useCells]
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===scEx_log:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(scEx_log)
})


# scExLog matrix with symbol as first column
# TODO
# we should probably just rename the rows and then have an option to tableSelectionServer that shows (or not) rownames
scExLogMatrixDisplay <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "scExLogMatrixDisplay")
    }
  )
  start.time <- Sys.time()
  if (DEBUG) {
    cat(file = stderr(), "scExLogMatrixDisplay\n")
  }
  # dataTables = inputData()
  # useCells = useCells()
  # useGenes = useGenes()
  scEx_log <- scEx_log()
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "scExLogMatrixDisplay:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("Calculating scExLogmatrix",
      id = "scExLogMatrixDisplay",
      duration = NULL
    )
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scExLogMatrixDisplay.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/scExLogMatrixDisplay.RData")

  # TODO
  if (ncol(scEx_log) > 20000) {

  }
  retVal <- as.data.frame(as.matrix(assays(scEx_log)[[1]]))
  rownames(retVal) <- make.names(rowData(scEx_log)$symbol, unique = TRUE)

  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===scExLogMatrixDisplay:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(retVal)
})

pcaFunc <- function(scEx_log) {
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/pcaFunc.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/pcaFunc.RData")
  scaterPCA <- tryCatch({
    # not sure, but this works on another with dgTMatrix
    if (class(assays(scEx_log)[["logcounts"]]) == "dgTMatrix") {
      assays(scEx_log)[["logcounts"]] <- as(assays(scEx_log)[["logcounts"]], "dgCMatrix")
    }
    scater::runPCA(scEx_log,
      ncomponents = 10, method = "irlba",
      ntop = 500, exprs_values = "logcounts"
    )
  },
  error = function(e) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "Problem with PCA, probably not enough cells?",
        type = "warning",
        duration = NULL
      )
    }
    cat(file = stderr(), "PCA FAILED!!!\n")
    return(NULL)
  }
  )
  if (is.null(scaterPCA)) return(NULL)
  # pca = reducedDim(scaterPCA, "PCA")
  # attr(pca,"percentVar")
  #
  return(list(
    x = SingleCellExperiment::reducedDim(scaterPCA, "PCA"),
    var_pcs = attr(SingleCellExperiment::reducedDim(scaterPCA, "PCA"), "percentVar")
  ))
}

pca <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "pca")
    }
  )
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "pca\n")
  }
  scEx_log <- scEx_log()
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "pca:NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("pca", id = "pca", duration = NULL)
  }
  retVal <- pcaFunc(scEx_log)
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===pca:donedone", difftime(end.time, start.time, units = "min"), "\n")
  }

  return(retVal)
})


kmClusteringFunc <- function(pca, seed, kNr = 10) {
  clustering <- list()

  # kNr = 10
  # for (kNr in 2:kNr) {
  set.seed(seed = seed)
  km <- cellrangerRkit::run_kmeans_clustering(pca, k = kNr)
  clustering[[paste0("kmeans_", kNr, "_clusters")]] <- data.frame(
    "Barcode" = rownames(data.frame(km$cluster)),
    "Cluster" = km$cluster
  )
  # }
  return(clustering)
}

kmClustering <- reactive({
  on.exit(
    removeNotification(id = "kmClustering")
  )
  start.time <- Sys.time()
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
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===kmClustering:done", difftime(end.time, start.time, units = "min"), "\n")
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


# Here, we store projections that are created during the session. These can be selections of cells or other values that
# are not possible to precalculate.
sessionProjections <- reactiveValues(
  prjs = data.frame()
)

projections <- reactive({
  start.time <- Sys.time()
  # scEx is the fundamental variable with the raw data, which is available after loading
  # data. Here we ensure that everything is loaded and all varialbles are set by waiting
  # input data being loaded
  scEx <- scEx()
  pca <- pca()
  prjs <- (sessionProjections$prjs)
  if (!exists("scEx") | is.null(scEx) | !exists("pca") | is.null(pca)) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  projections <- data.frame(pca$x[, c(1, 2, 3)])
  if (DEBUG) {
    cat(file = stderr(), "projections\n")
  }

  # phenotypic data/ annotations of cells can already be included in the scEx object. We collect this information, but only for variable that hold information
  # i.e. length(levels) > 1 & < number of rows
  pd <- colData(scEx)
  if (ncol(pd) < 2) {
    cat(file = stderr(), "phenoData for scEx has less than 2 columns\n")
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
      start.time1 <- Sys.time()
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
        if (DEBUG) {
          end.time <- Sys.time()
          cat(file = stderr(), "===", proj[1], ":done", difftime(end.time, start.time1, units = "min"), "\n")
        }
      }

      colnames(projections) <- cn
      observe(proj[2], quoted = TRUE)
    }
  })
  # add a column for gene specific information that will be filled/updated on demand
  projections$UmiCountPerGenes <- 0
  projections$UmiCountPerGenes2 <- 0
  for (pdIdx in colnames(pd)) {
    if (!pdIdx %in% colnames(projections)) {
      projections[, pdIdx] <- pd[, pdIdx]
    }
  }
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===projections:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  if (ncol(prjs) > 0 & nrow(prjs) == nrow(projections)) {
    projections <- cbind(projections, prjs)
  }
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
  scEx <- scEx()
  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/initializeGroupNames.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/initializeGroupNames.RData")
  isolate({
    df <- data.frame(all = rep(TRUE, ncol(scEx)), none = rep(FALSE, ncol(scEx)))
    rownames(df) <- colnames(scEx)
    cat(file = stderr(), "initializeGroupNames2\n")
    groupNames[["namesDF"]] <- df
    cat(file = stderr(), "initializeGroupNames3\n")
  })
})


dbCluster <- reactive({
  start.time <- Sys.time()

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

  # dbCluster <- factor(clustering[[paste0("kmeans_", kNr, "_clusters")]]$Cluster - 1)
  dbCluster <- factor(clustering[[paste0("kmeans_", kNr, "_clusters")]]$Cluster)

  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===dbCluster:done", difftime(end.time, start.time, units = "min"), "\n")
  }

  return(dbCluster)
})

# sample --------
sample <- reactive({
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "sample\n")
  }
  scEx <- scEx()
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "sample: NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/sample.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/sample.RData")
  # samp <- gsub(".*-(.*)", "\\1", colnames(scEx))
  # if (length(levels(as.factor(samp))) > 1) {
  #   sample <- samp
  # } else {
  #   sample <- rep("1", ncol(scEx))
  # }
  pd <- colData(scEx)
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
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===sample:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  retVal
  # return(sample)
})

# geneCount --------
geneCount <- reactive({
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "geneCount\n")
  }
  scEx <- scEx()
  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/geneCount.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/geneCount.RData")
  retVal <- Matrix::colSums(assays(scEx)[["counts"]] > 0)
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===geneCount:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(retVal)
})

beforeFilterPrj <- reactive({
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "umiCount\n")
  }
  scEx <- scEx()
  bfc <- beforeFilterCounts()
  if (is.null(scEx) | is.null(bfc)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/beforeFilterPrj.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/beforeFilterPrj.RData")
  cn <- colnames(scEx)
  retVal <- bfc[cn]
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===beforeFilterPrj:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(retVal)
})

umiCount <- reactive({
  start.time <- Sys.time()

  if (DEBUG) {
    cat(file = stderr(), "umiCount\n")
  }
  scEx <- scEx()
  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/umiCount.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/umiCount.RData")
  retVal <- Matrix::colSums(assays(scEx)[["counts"]])
  if (DEBUG) {
    end.time <- Sys.time()
    cat(file = stderr(), "===umiCount:done", difftime(end.time, start.time, units = "min"), "\n")
  }
  return(retVal)
})


sampleInfoFunc <- function(scEx) {
  # gsub(".*-(.*)", "\\1", scEx$barcode)
  colData(scEx)$sampleNames
}


# sampleInfo -------
# sample information
sampleInfo <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "sampleInfo\n")
  }
  scEx <- scEx()
  if (!exists("scEx")) {
    if (DEBUG) {
      cat(file = stderr(), "sampleInfo: NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/sampleInfo.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/sampleInfo.RData")

  ret <- sampleInfoFunc(scEx)
  if (DEBUG) {
    cat(file = stderr(), "sampleInfo: done\n")
  }
  return(ret)
})


# table of input cells with sample information
inputSample <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "inputSample")
    }
  )
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
  # TODO should come from sampleInfo
  sampInf <- gsub(".*-(.*)", "\\1", dataTables$scEx$barcode)
  cellIds <- data.frame(
    cellName = colnames(dataTables$scEx),
    sample = sampInf,
    ngenes = Matrix::colSums(assays(dataTables$scEx)[[1]])
  )

  if (DEBUG) {
    cat(file = stderr(), "inputSample: done\n")
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
# TODO change to scEx_log everywhere and remove
log2cpm <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "log2cpm\n")
  }
  scEx_log <- scEx_log()
  if (is.null(scEx_log)) {
    if (DEBUG) {
      cat(file = stderr(), "log2cpm: NULL\n")
    }
    return(NULL)
  }
  log2cpm <- as.data.frame(as.matrix(assays(scEx_log)[[1]]))

  return(log2cpm)
})


# dummy function to return NULL
returnNull <- function() {
  return(NULL)
}

#### plot2Dprojection ----------------
# used in moduleServer and reports
plot2Dprojection <- function(scEx_log, scEx, projections, g_id, featureData,
                             geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
                             logx = FALSE, logy = FALSE, divXBy = "None", divYBy = "None") {
  geneid <- geneName2Index(g_id, featureData)

  if (length(geneid) == 0) {
    return(NULL)
  }
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ,drop=FALSE]
  # } else {
  expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, , drop = FALSE])
  # }
  validate(need(is.na(sum(expression)) != TRUE, ""))
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ]
  # } else {
  #   expression <- Matrix::colSums(exprs(scEx_log)[geneid, ])
  # }
  # validate(need(is.na(sum(expression)) != TRUE, ""))

  # geneid <- geneName2Index(geneNames, featureData)
  projections <- updateProjectionsWithUmiCount(
    dimX = dimX, dimY = dimY,
    geneNames = geneNames,
    geneNames2 = geneNames2,
    featureData = featureData,
    scEx = scEx, projections = projections
  )


  projections <- cbind(projections, expression)
  names(projections)[ncol(projections)] <- "exprs"

  if (DEBUG) {
    cat(file = stderr(), paste("output$dge_plot1:---", clId[1], "---\n"))
  }
  subsetData <- subset(projections, dbCluster %in% clId)
  # subsetData$dbCluster = factor(subsetData$dbCluster)
  # if there are more than 18 samples ggplot cannot handle different shapes and we ignore the
  # sample information
  if (length(as.numeric(as.factor(subsetData$sample))) > 18) {
    subsetData$shape <- as.factor(1)
  } else {
    subsetData$shape <- as.numeric(as.factor(subsetData$sample))
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/clusterPlot.RData", list = c(ls(), "legend.position", ls(envir = globalenv())))
    cat(file = stderr(), paste("plot2Dprojection saving done.\n"))
  }
  # load(file="~/scShinyHubDebug/clusterPlot.RData")
  if (nrow(subsetData) == 0) return(NULL)
  # subsetData$shape = as.factor(1)
  gtitle <- paste(toupper(g_id), clId, sep = "-Cluster", collapse = " ")
  if (nchar(gtitle) > 50) {
    gtitle <- paste(substr(gtitle, 1, 50), "...")
  }

  require(plotly)
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  if (divXBy != "None") {
    subsetData[, dimX] <- subsetData[, dimX] / subsetData[, divXBy]
  }
  if (divYBy != "None") {
    subsetData[, dimY] <- subsetData[, dimY] / subsetData[, divYBy]
  }

  typeX <- typeY <- "linear"
  if (logx) {
    typeX <- "log"
  }
  if (logy) {
    typeY <- "log"
  }
  if (is.factor(subsetData[, dimX])) {
    typeX <- NULL
  }
  if (is.factor(subsetData[, dimY])) {
    typeY <- NULL
  }
  xAxis <- list(
    title = dimX,
    titlefont = f,
    type = typeX
  )
  yAxis <- list(
    title = dimY,
    titlefont = f,
    type = typeY
  )
  p1 <- plot_ly(data = subsetData, source = "subset") %>%
    add_trace(
      x = ~ get(dimX),
      y = ~ get(dimY),
      type = "scatter",
      mode = "markers",
      text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
      marker = list(
        color = subsetData[, "exprs"],
        size = 4
      )
    ) %>%
    # add_trace() %>%
    layout(
      xaxis = xAxis,
      yaxis = yAxis,
      title = gtitle,
      dragmode = "select"
    )

  selectedCells <- NULL
  if (length(grpN) > 0) {
    if (length(grpNs[rownames(subsetData), grpN]) > 0 & sum(grpNs[rownames(subsetData), grpN], na.rm = TRUE) > 0) {
      grpNSub <- grpNs[rownames(subsetData), ]
      selectedCells <- rownames(grpNSub[grpNSub[, grpN], ])
    }
  }
  if (!is.null(selectedCells)) {
    shape <- rep("a", nrow(subsetData))
    selRows <- which(rownames(subsetData) %in% selectedCells)
    shape[selRows] <- "b"
    x1 <- subsetData[selectedCells, dimX, drop = FALSE]
    y1 <- subsetData[selectedCells, dimY, drop = FALSE]
    p1 <- p1 %>%
      add_trace(
        x = x1[, 1], y = y1[, 1],
        marker = list(
          color = rep("green", nrow(x1)),
          size = 5
        ),
        text = ~ paste(
          rownames(subsetData[selectedCells, ]),
          "<br />", subsetData$exprs[selectedCells]
        ),
        type = "scatter",
        mode = "markers",
        name = "selected"
      )
  }
  p1
}

