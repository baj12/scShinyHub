# heatmapFunc ---------------------------------
# used by bot selection and all
coE_heatmapFunc <- function(featureData, gbm_matrix, projections, genesin, cells) {
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "heatmap")
  )
  #  create parameters used for pheatmap module
  genesin <- geneName2Index(genesin, featureData)
  expression <- gbm_matrix[genesin, cells]

  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or genes not expressed"
  ))

  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]

  # expression <- expression[, rownames(projections)]
  # expression <- expression[GeneBCMatrix::complete.cases(expression), ]

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
  # outfile <- paste0(tempdir(), "/heatmap", base::sample(1:10000, 1), ".png")
  # cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  # this can fail with na/inf in hclust error message if there is a row with all the same values
  # med = median(as.vector(as.matrix(expression)))
  # stDev = sd(as.vector(as.matrix(expression)))
  # minBreak = max(0, med - 3* stDev)
  # maxBreak = med + 3* stDev
  # stepBreak = (maxBreak - minBreak) / 6
  nonZeroRows <- which(Matrix::rowSums(expression) > 0)
  retVal <- list(
    mat = expression[nonZeroRows, order(annotation[, 1], annotation[, 2])],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    fontsize_row = 14,
    labels_col = colnames(expression),
    labels_row = featureData[rownames(expression), "Associated.Gene.Name"],
    show_rownames = TRUE,
    annotation_col = annotation,
    show_colnames = FALSE,
    annotation_legend = TRUE,
    # breaks = seq(minBreak, maxBreak, by = stepBreak),
    # filename = 'test.png',
    # filename = normalizePath(outfile),
    colorRampPalette(rev(brewer.pal(
      n = 6, name =
        "RdBu"
    )))(6)
  )
  retVal
}


# heatmapSelectedReactive ----
# reactive function for selected heatmap
heatmapSelectedReactive <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "selectedHeatmap")
  )
  if (DEBUG) {
    cat(file = stderr(), "output$heatmapSelectedReactive\n")
  }
  featureData <- featureDataReact()
  gbm_matrix <- gbm()
  projections <- projections()
  genesin <- input$heatmap_geneids2
  sc <- selctedCluster()
  scCL <- sc$cluster
  # scBP = sc$brushedPs()
  scCells <- sc$selectedCells()

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
  if (is.null(featureData) ||
    is.null(gbm_matrix) ||
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
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("selectedheatmap", id = "selectedHeatmap", duration = NULL)
  }

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")

  # subsetData <-
  #   subset(projections, as.numeric(as.character(projections$dbCluster)) %in% scCL)
  # cells.1 <- rownames(brushedPoints(subsetData, scBP))
  cells.1 <- scCells
  retval <- coE_heatmapFunc(featureData, gbm_matrix, projections, genesin, cells = cells.1)

  return(retval)
})

topExpGenesTable <- reactive({
  if (DEBUG) cat(file = stderr(), "output$topExpGenes\n")
  featureData <- featureDataReact()
  gbm_log <- gbm_log()
  coEtgPerc <- input$coEtgPerc
  coEtgminExpr <- input$coEtgMinExpr
  sc <- selctedCluster()
  scCL <- sc$cluster
  # scBP = sc$brushedPs()
  scCells <- sc$selectedCells()

  if (is.null(featureData) || is.null(gbm_log) || is.null(scCells)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/output_topExpGenes.RData", list = c("scCells", "scCL", "sc", ls()))
  }
  # load(file="~/scShinyHubDebug/output_topExpGenes.RData")
  # we only work on cells that have been selected
  # mat <- as.matrix(exprs(gbm_log))[, scCells]
  mat <- exprs(gbm_log)[, scCells]
  # only genes that express at least coEtgminExpr UMIs
  mat[mat < coEtgminExpr] <- 0
  # only genes that are expressed in coEtgPerc or more cells
  allexpressed <- Matrix::rowSums(mat > 0) / length(scCells) * 100 >= coEtgPerc
  mat <- mat[allexpressed, ]

  cv <- function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  }
  matCV <- apply(mat, 1, cv)
  # top.genes <- as.data.frame(exprs(gbm_log))
  maxRows <- min(nrow(mat), 200)
  top.genesOrder <- order(matCV, decreasing = TRUE)[1:maxRows]
  if (dim(mat)[1] > 0) {
    mat <- mat[top.genesOrder, ]
    fd <- featureData[rownames(mat), c("Associated.Gene.Name", "Description")]
    matCV <- matCV[rownames(mat)]
    fd <- cbind(fd, matCV)
    colnames(fd) <- c("gene", "description", "CV")
    # since we are returning a table to be plotted, we convert to regular table (non-sparse)
    outMat <- cbind2(fd, as.matrix(mat))
    rownames(outMat) <- make.unique(as.character(outMat$gene), sep = "___")
    return(outMat)
  } else {
    return(NULL)
  }
})


# TODO module for heatmap?
# output$selectedHeatmap <- renderImage({
#   if (DEBUG) {
#     cat(file = stderr(), "output$selectedHeatmap\n")
#   }
#   featureData <- featureDataReact()
#   gbm_matrix <- gbm_matrix()
#   projections <- projections()
#   genesin <- input$heatmap_geneids2
#   sc <- selctedCluster()
#   scCL <- sc$cluster
#   # scBP = sc$brushedPs()
#   scCells <- sc$selectedCells()
#
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
#   if (is.null(featureData) |
#     is.null(gbm_matrix) |
#     is.null(projections) | is.null(scCells) | length(scCells) == 0) {
#     return(
#       list(
#         src = "empty.png",
#         contentType = "image/png",
#         width = 96,
#         height = 96,
#         alt = "heatmap should be here"
#       )
#     )
#   }
#   if (!is.null(getDefaultReactiveDomain())) {
#     showNotification("selectedheatmap", id = "selectedHeatmap", duration = NULL)
#   }
#
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/selectedHeatmap.RData", list = c(ls(), ls(envir = globalenv())))
#   }
#   # load(file = "~/scShinyHubDebug/selectedHeatmap.RData")
#
#   # subsetData <-
#   #   subset(projections, as.numeric(as.character(projections$dbCluster)) %in% scCL)
#   # cells.1 <- rownames(brushedPoints(subsetData, scBP))
#   cells.1 <- scCells
#   retval <- heatmapFunc(featureData, gbm_matrix, projections, genesin, cells = cells.1)
#
#   if (!is.null(getDefaultReactiveDomain())) {
#     removeNotification(id = "selectedHeatmap")
#   }
#   return(retval)
# })

# plotCoExpressionFunc ------
# used in binarized 2D plot
plotCoExpressionFunc <-
  function(featureData,
             gbm_log,
             upI,
             projections,
             genesin,
             cl3,
             dimx3,
             dimy3) {
    genesin <- geneName2Index(genesin, featureData)
    # genesin <- toupper(genesin)
    # genesin <- gsub(" ", "", genesin, fixed = TRUE)
    # genesin <- strsplit(genesin, ',')

    subsetData <- subset(projections, dbCluster %in% cl3)
    # cells.1 <- rownames(subsetData)


    # map <-
    #   rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
    # if(DEBUG)cat(file=stderr(),map[1])

    expression <- gbm_log[genesin, ]
    # if(DEBUG)cat(file=stderr(),rownames(expression))

    # expression<-expression[complete.cases(expression),]
    # if(DEBUG)cat(file=stderr(),rownames(expression))

    # display genes not found
    # notFound = genesin[[1]][which(!genesin[[1]] %in% featureData$Associated.Gene.Name)]
    # if (length(notFound) > 0) {
    #   if (!is.null(getDefaultReactiveDomain())) {
    #     showNotification(
    #       paste("following genes were not found", notFound, collapse = " "),
    #       id = "plotCoExpressionNotFound",
    #       type = "warning",
    #       duration = 20
    #     )
    #   }
    # }

    validate(need(
      is.na(sum(expression)) != TRUE,
      "Gene symbol incorrect or genes not expressed"
    ))

    bin <- expression
    bin[] <- 0

    for (i in 1:nrow(expression)) {
      x <- Mclust(expression[i, ], G = 2)
      bin[i, ] <- x$classification
    }
    bin <- bin - 1
    allexprs <- apply(bin, 2, sum)
    plotexprs <- allexprs
    plotexprs[] <- 0
    plotexprs[allexprs >= length(rownames(bin))] <- 1
    # TODO positiveCells is changing too often. Maybe this can be controlled a bit more? check if changed? Use global variable not a reactive value?
    # positiveCells$positiveCells <- allexprs >= length(rownames(bin))
    # positiveCells$positiveCellsAll <- plotexprs

    mergeExprs <- plotexprs[rownames(subsetData)]
    # if(DEBUG)cat(file=stderr(),length(mergeExprs))

    subsetData$CoExpression <- factor(mergeExprs)
    subsetData$dbCluster <- as.factor(subsetData$dbCluster)
    p1 <-
      ggplot(
        subsetData,
        aes_string(x = dimx3, y = dimy3)
      ) +
      geom_point(aes_string(
        shape = "sampleNames",
        # alpha = 'CoExpression',
        color = "dbCluster"
      ),
      size = 4
      ) +
      theme_bw()

    if (DEBUG) {
      cat(file = stderr(), "output$plotCoExpression:done\n")
    }
    return(p1)
  }

# geneGrp_vioFunc ------
geneGrp_vioFunc <- function(genesin, projections, gbm, featureData, minExpr = 1,
                            dbCluster, showPermutations = FALSE) {
  require(gtools)
  require(stringr)
  start.time <- Sys.time()
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")[[1]]
  # vIdx=1
  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
  
  map <-
    rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin), ])
  if (DEBUG) {
    cat(file = stderr(), length(map))
  }
  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
  if (length(map) == 0) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(
        "no genes found",
        id = "heatmapWarning",
        type = "warning",
        duration = 20
      )
    }
    return(
      NULL
    )
  }

  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
  expression <- Matrix::colSums(exprs(gbm[map, ]) >= minExpr)
  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
  ylabText <- "number genes from list"
  # projections = projections[,1:12]
  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
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
          rownames(featureData[which(featureData$Associated.Gene.Name %in% comb[cIdx, ]), ])

        permIdx <- Matrix::colSums(exprs(gbm[map, ]) >= minExpr) == length(comb[cIdx, ])
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
  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
  


  # if(class(projections[,dbCluster])=="factor"){

  p1 <-
    ggplot(projections, aes_string(factor(projections[, dbCluster]), "coExpVal",
      fill = factor(projections[, dbCluster])
    )) +
    geom_violin(scale = "count") +
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
  # }else{
  #   return(NULL)
  # }
  if (DEBUG) {
    cat(file = stderr(), "output$gene_vio_plot:done\n")
  }
  # cat(file = stderr(), paste("===violin-", vIdx,"-", difftime(Sys.time(), start.time, units = "min"), " min\n")); vIdx = vIdx+1;start.time <- Sys.time()
  return(p1)
}

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

  rownames(res2$globalBmus) <- make.unique(as.character(rownames(iData)), sep = "___")
  simGenes <- rownames(res2$globalBmus)[which(res2$globalBmus[, 1] == res2$globalBmus[geneName, 1] &
    res2$globalBmus[, 2] == res2$globalBmus[geneName, 2])]
  return(simGenes)
}

heatmapSOMReactive <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "heatmapSOMReactive")
  )

  start.time <- Sys.time()
  if (DEBUG) {
    cat(file = stderr(), "output$somReactive\n")
  }
  gbm_matrix <- as.matrix(exprs(gbm()))
  projections <- projections()
  genesin <- input$geneSOM
  nSOM <- input$dimSOM
  featureData <- featureDataReact()

  if (is.null(gbm_matrix)) {
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
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("somheatmap", id = "heatmapSOMReactive", duration = NULL)
  }
  # go from readable gene name to ENSG number
  genesin <- geneName2Index(genesin, featureData)

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/heatmapSOMReactive.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/heatmapSOMReactive.RData")

  # subsetData <-
  #   subset(projections, as.numeric(as.character(projections$dbCluster)) %in% scCL)
  # cells.1 <- rownames(brushedPoints(subsetData, scBP))

  # gbmOrg = gbm_matrix
  # gbm_matrix = as(gbm_matrix, "dgTMatrix")
  # rownames(gbm_matrix) = rownames(featureData)
  # dgTMatrix makes som crash R, I guess it is because it is calling a C function that is able to handle it...
  geneNames <- somFunction(iData = gbm_matrix, nSom = nSOM, geneName = genesin)
  output$somGenes <- renderText({
    paste(featureData[geneNames, "Associated.Gene.Name"], collapse = ", ", sep = ",")
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
  annotation <- data.frame(projections[, c("dbCluster", "sampleNames")])
  rownames(annotation) <- rownames(projections)
  colnames(annotation) <- c("Cluster", "sampleNames")

  retVal <- list(
    mat = gbm_matrix[geneNames, ],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    fontsize_row = 14,
    labels_row = featureData[geneNames, "Associated.Gene.Name"],
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
    )))(6)
  )
  # pheatmap(gbm_matrix[rownames(gbm_matrix)[1:10], ],
  #          cluster_rows = TRUE,
  #          cluster_cols = TRUE,
  #          scale = "row")
  # do.call(pheatmap, retVal)


  end.time <- Sys.time()
  cat(file = stderr(), paste("===heatmapSOMReactive:done", difftime(end.time, start.time, units = "min"), " min\n"))

  return(retVal)
})
