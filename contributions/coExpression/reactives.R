# heatmapFunc ---------------------------------
# used by both selection and all to create input for heatmap module
coE_heatmapFunc <- function(featureData, scEx_matrix, projections, genesin, cells, sampCol) {
  start.time <- base::Sys.time()
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/coE_heatmapFunc.RData", list = c(ls(envir = globalenv()), ls()))
  }
  # load(file = "~/scShinyHubDebug/coE_heatmapFunc.RData")
  
  #  create parameters used for pheatmap module
  genesin <- geneName2Index(genesin, featureData)
  if (is.null(genesin) | is.null(cells)) {
    return(NULL)
  }
  expression <- scEx_matrix[genesin, cells]
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or genes not expressed"
  ))
  
  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]
  
  if (!("sampleNames" %in% colnames(projections))) {
    projections$sample <- 1
  }
  annotation <- data.frame(projections[cells, c("dbCluster", "sampleNames")])
  rownames(annotation) <- colnames(expression)
  colnames(annotation) <- c("Cluster", "sampleNames")
  
  # For high-res displays, this will be greater than 1
  # pixelratio <- session$clientData$pixelratio
  pixelratio <- 1
  if (is.null(pixelratio)) pixelratio <- 1
  # width <- session$clientData$output_plot_width
  # height <- session$clientData$output_plot_height
  width <- NULL
  height <- NULL
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  
  nonZeroRows <- which(Matrix::rowSums(expression) > 0)
  
  annCols <- list("sampleNames" = sampCol)
  
  retVal <- list(
    mat = expression[nonZeroRows, order(annotation[, 1], annotation[, 2])],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    fontsize_row = 14,
    labels_col = colnames(expression),
    labels_row = featureData[rownames(expression), "symbol"],
    show_rownames = TRUE,
    annotation_col = annotation,
    show_colnames = FALSE,
    annotation_legend = TRUE,
    # breaks = seq(minBreak, maxBreak, by = stepBreak),
    # filename = 'test.png',
    # filename = normalizePath(outfile),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(
      n = 6, name =
        "RdBu"
    )))(6),
    annotation_colors = annCols
  )
  
  # print debugging information on the console
  printTimeEnd(start.time, "inputData")
  # for automated shiny testing using shinytest
  exportTestValues(coE_heatmapFunc = {
    retVal
  })
  
  # this is what is run in the module
  # do.call(TRONCO::pheatmap, retVal)
  return(retVal)
}


# heatmapSelectedReactive ----
# reactive function for selected heatmap
heatmapSelectedReactive <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "selectedHeatmap")
    }
  )
  if (DEBUG) {
    cat(file = stderr(), "output$heatmapSelectedReactive\n")
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("selectedheatmap", id = "selectedHeatmap", duration = NULL)
  }
  
  scEx <- scEx()
  projections <- projections()
  genesin <- input$heatmap_geneids2
  sc <- selctedCluster()
  scCL <- sc$cluster
  scCells <- sc$selectedCells()
  sampCol <- sampleCols$colPal
  
  if (is.null(scEx) ||
      is.null(projections) || is.null(scCells) || length(scCells) == 0) {
    return(
      list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "heatmap should be here"
      )
    )
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
  
  
  scEx_matrix <- assays(scEx)[["counts"]]
  featureData <- rowData(scEx)
  
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
  
  retval <- coE_heatmapFunc(featureData, scEx_matrix, projections, genesin,
                            cells = scCells, sampCol = sampCol
  )
  # print debugging information on the console
  printTimeEnd(start.time, "inputData")
  # for automated shiny testing using shinytest
  exportTestValues(heatmapSelectedReactive = {retVal})  
  return(retval)
})

# topExpGenesTable ----
#' topExpGenesTable
#' in coexpressionSelected tab, showing the table of top expressed genes for a given 
#' selection
#' coEtgPerc = genes shown have to be expressed in at least X % of cells
#' coEtgMinExpr = genes shown have at least to X UMIs expressed 
topExpGenesTable <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "topExpGenesTable")
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("topExpGenesTable", id = "topExpGenesTable", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "topExpGenesTable\n")
  
  scEx_log <- scEx_log()
  coEtgPerc <- input$coEtgPerc
  coEtgminExpr <- input$coEtgMinExpr
  sc <- selctedCluster()
  scCL <- sc$cluster
  scCells <- sc$selectedCells()
  
  if (is.null(scEx_log) || is.null(scCells)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/output_topExpGenes.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/output_topExpGenes.RData")
  
  featureData <- rowData(scEx_log)
  # we only work on cells that have been selected
  mat <- assays(scEx_log)[[1]][, scCells]
  # only genes that express at least coEtgminExpr UMIs
  mat[mat < coEtgminExpr] <- 0
  # only genes that are expressed in coEtgPerc or more cells
  allexpressed <- Matrix::rowSums(mat > 0) / length(scCells) * 100 >= coEtgPerc
  mat <- mat[allexpressed, ]
  
  cv <- function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  }
  matCV <- apply(mat, 1, cv)
  # top.genes <- as.data.frame(exprs(scEx_log))
  maxRows <- min(nrow(mat), 200)
  top.genesOrder <- order(matCV, decreasing = TRUE)[1:maxRows]
  retVal = NULL
  if (dim(mat)[1] > 0) {
    mat <- mat[top.genesOrder, ]
    fd <- featureData[rownames(mat), c("symbol", "Description")]
    matCV <- matCV[rownames(mat)]
    fd <- cbind(fd, matCV)
    colnames(fd) <- c("gene", "description", "CV")
    # since we are returning a table to be plotted, we convert to regular table (non-sparse)
    outMat <- cbind2(fd, as.matrix(mat))
    rownames(outMat) <- make.unique(as.character(outMat$gene), sep = "___")
    retVal = as.data.frame(outMat)
  } 
  
  printTimeEnd(start.time, "topExpGenesTable")
  exportTestValues(topExpGenesTable = {retVal})  
  return(retVal)
})



#' geneGrp_vioFunc
#' generates a ggplot object with a violin plot
#' optionally creates all permutations.
geneGrp_vioFunc <- function(genesin, projections, scEx, featureData, minExpr = 1,
                            dbCluster, showPermutations = FALSE, sampCol) {
  if (DEBUG) cat(file = stderr(), "geneGrp_vioFunc\n")
  require(gtools)
  require(stringr)
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")[[1]]
  
  map <- rownames(featureData[which(featureData$symbol %in% genesin), ])
  
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/geneGrp_vioFunc.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/geneGrp_vioFunc.RData")
  
  if (length(map) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no genes found",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }
  
  expression <- Matrix::colSums(assays(scEx)[["counts"]][map, ] >= minExpr)
  ylabText <- "number genes from list"
  
  if (showPermutations) {
    perms <- rep("", length(expression))
    ylabText <- "Permutations"
    xPerm <- length(genesin)
    if (xPerm > 10) {
      xPerm <- 10
      warning("reducing number of permutations to 10")
    }
    # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
    for (r in 1:xPerm) {
      comb <- combinations(xPerm, r, genesin)
      for (cIdx in 1:nrow(comb)) {
        map <-
          rownames(featureData[which(featureData$symbol %in% comb[cIdx, ]), ])
        # permIdx <- Matrix::colSums(exprs(gbm[map, ]) >= minExpr) == length(comb[cIdx, ])
        
        permIdx <- Matrix::colSums(assays(scEx)[["counts"]][map, , drop = FALSE] >= minExpr) == length(comb[cIdx, ])
        perms[permIdx] <- paste0(comb[cIdx, ], collapse = "+")
      }
    }
    # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
    perms <- factor(perms)
    permsNames <- levels(perms)
    permsNum <- unlist(lapply(strsplit(permsNames, "\\+"), length))
    perms <- factor(as.character(perms), levels = permsNames[order(permsNum)])
    permsNames <- str_wrap(levels(perms))
    perms <- as.integer(perms)
    projections <- cbind(projections, coExpVal = perms)
  } else {
    projections <- cbind(projections, coExpVal = expression)
    permsNames <- as.character(1:max(expression))
  }
  
  prj <- factor(projections[, dbCluster])
  mycolPal <- colorRampPalette(brewer.pal(
    n = 6, name =
      "RdYlBu"
  ))(length(levels(prj)))
  
  if (dbCluster == "sampleNames") {
    mycolPal <- sampCol
  }
  
  p1 <-
    ggplot(projections, aes_string(prj, "coExpVal",
                                   fill = factor(projections[, dbCluster])
    )) +
    geom_violin(scale = "count") +
    scale_fill_manual(values = mycolPal, aesthetics = "fill") +
    stat_summary(
      fun.y = median,
      geom = "point",
      size = 5,
      color = "black"
    ) +
    stat_summary(fun.data = n_fun, geom = "text") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 60,
        size = 12,
        vjust = 0.5
      ),
      axis.text.y = element_text(size = 10),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 12),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      legend.position = "none"
    ) +
    xlab(dbCluster) +
    
    scale_y_continuous(breaks = 1:length(permsNames), labels = str_wrap(permsNames)) +
    ylab(ylabText)
  
  if (DEBUG) cat(file = stderr(), "output$gene_vio_plot:done\n")
  return(p1)
}

#' somFunction
#' iData = expression matrix, rows = genes
#' cluster genes in SOM 
#' returns genes associated together within a som-node
somFunction <- function(iData, nSom, geneName) {
  require(kohonen)
  require(Rsomoclu)
  if (sum(geneName %in% rownames(iData)) == 0) return(NULL)
  res2 <- Rsomoclu.train(
    input_data = iData,
    nSomX = nSom, nSomY = nSom,
    nEpoch = 10,
    radius0 = 0,
    radiusN = 0,
    radiusCooling = "linear",
    mapType = "planar",
    gridType = "rectangular",
    scale0 = 1,
    scaleN = 0.01,
    scaleCooling = "linear"
  )
  
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/somFunction.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/somFunction.RData")
  rownames(res2$globalBmus) <- make.unique(as.character(rownames(iData)), sep = "___")
  simGenes <- rownames(res2$globalBmus)[which(res2$globalBmus[, 1] == res2$globalBmus[geneName, 1] &
                                                res2$globalBmus[, 2] == res2$globalBmus[geneName, 2])]
  return(simGenes)
}

# heatmapSOMReactive ----
#' heatmapSOMReactive
#' calculates a self organizing map (SOM) on the expression data and identifies genes
#' that cluster together with a gene of interest
# TODO: check that we are using only raw counts and not normalized
heatmapSOMReactive <- reactive({
  start.time <- Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "heatmapSOMReactive")
    }
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("somheatmap", id = "heatmapSOMReactive", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$somReactive\n")
  
  scEx <- scEx()
  projections <- projections()
  genesin <- input$geneSOM
  nSOM <- input$dimSOM
  sampCol <- sampleCols$colPal
  
  if (is.null(scEx)) {
    return(
      list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "heatmap should be here"
      )
    )
  }

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/heatmapSOMReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/heatmapSOMReactive.RData")

  scEx_matrix <- as.matrix(assays(scEx)[["counts"]])
  featureData <- rowData(scEx)
  # go from readable gene name to ENSG number
  genesin <- geneName2Index(genesin, featureData)

  geneNames <- somFunction(iData = scEx_matrix, nSom = nSOM, geneName = genesin)

  # plot the genes found
  output$somGenes <- renderText({
    paste(featureData[geneNames, "symbol"], collapse = ", ", sep = ",")
  })
  if (length(geneNames) < 2) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no additional gene found. reduce size of SOM",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(NULL)
  }

  # create variables for heatmap module
  annotation <- data.frame(projections[, c("dbCluster", "sampleNames")])
  rownames(annotation) <- rownames(projections)
  colnames(annotation) <- c("Cluster", "sampleNames")
  annCols <- list("sampleNames" = sampCol)
  
  retVal <- list(
    mat = scEx_matrix[geneNames, ],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    fontsize_row = 14,
    labels_row = featureData[geneNames, "symbol"],
    show_rownames = TRUE,
    annotation_col = annotation,
    show_colnames = FALSE,
    annotation_legend = TRUE,
    # breaks = seq(minBreak, maxBreak, by = stepBreak),
    # filename = 'test.png',
    # filename = normalizePath(outfile),
    color = colorRampPalette(rev(brewer.pal(
      n = 6, name =
        "RdBu"
    )))(6),
    annotation_colors = annCols
  )
  # system.time(do.call(TRONCO::pheatmap, retVal))
  
  printTimeEnd(start.time, "inputData")
  exportTestValues(heatmapSOMReactive = {retVal })  
  return(retVal)
})

# updateInputXviolinPlot----
#' updateInputXviolinPlot
#' Update x/y axis selection possibilities for violin plot
#' could probably be an observer, but it works like this as well...
updateInputXviolinPlot <- reactive({
  tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(
    session,
    "dimension_x3",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[1]
  )
  updateSelectInput(
    session,
    "dimension_y3",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[2]
  )
  
  coln <- colnames(tsneData)
  choices <- c()
  for (cn in coln) {
    if (length(levels(as.factor(tsneData[, cn]))) < 20) {
      choices <- c(choices, cn)
    }
  }
  if (length(choices) == 0) {
    choices <- c("no valid columns")
  }
  updateSelectInput(
    session,
    "dimension_xVioiGrp",
    choices = choices,
    selected = choices[1]
  )
})


# heatmapReactive -------
# reactive for module pHeatMapModule
# for all clusters menu item
heatmapReactive <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "heatmap")
  )
  if (!is.null(getDefaultReactiveDomain()))
    removeNotification(id = "heatmapWarning")
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("heatmap", id = "heatmap", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$heatmap\n")

  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$heatmap_geneids
  sampCol = sampleCols$colPal
  
  if (is.null(scEx_log) | is.null(projections)) {
    return(list(
      src = "empty.png",
      contentType = "image/png",
      width = 96,
      height = 96,
      alt = "heatmap should be here"
    ))
  }
  
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/heatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/heatmap.RData")

  featureData <- rowData(scEx_log)
  scEx_matrix <- as.matrix(assays(scEx_log)[["logcounts"]])
  retVal <- coE_heatmapFunc(
    featureData = featureData, scEx_matrix = scEx_matrix,
    projections = projections, genesin = genesin, cells = colnames(scEx_matrix),
    sampCol = sampCol
  )
  
  printTimeEnd(start.time, "DummyReactive")
  exportTestValues(heatmapReactive = {retVal})  
  return(retVal)
})
