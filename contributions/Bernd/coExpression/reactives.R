# heatmapFunc ---------------------------------
# used by bot selection and all
coE_heatmapFunc <- function(featureData, gbm_matrix, projections, genesin, cells) {
  #  create parameters used for pheatmap module
  if(length(genesin)==0 | length(cells)==0) return(NULL)
  genesin <- geneName2Index(genesin, featureData)
  expression <- gbm_matrix[genesin, cells]

  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or genes not expressed"
  ))

  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]

  # expression <- expression[, rownames(projections)]
  expression <- expression[complete.cases(expression), ]

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
  nonZeroRows <- which(rowSums(expression) > 0)
  retVal <- list(
    mat = as.matrix(expression)[nonZeroRows, order(annotation[, 1], annotation[, 2])],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    fontsize_row = 10,
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
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmap")
  }
  retVal
}


# heatmapSelectedReactive ----
# reactive function for selected heatmap
heatmapSelectedReactive <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "output$selectedHeatmap\n")
  }
  featureData <- featureDataReact()
  gbm_matrix <- gbm_matrix()
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
  if (is.null(featureData) |
    is.null(gbm_matrix) |
    is.null(projections) | is.null(scCells) | length(scCells) == 0) {
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

  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "selectedHeatmap")
  }
  return(retval)
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
    positiveCells$positiveCells <- allexprs >= length(rownames(bin))
    positiveCells$positiveCellsAll <- plotexprs

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
geneGrp_vioFunc <- function(genesin, projections, gbm, featureData, minExpr = 1, dbCluster) {
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")

  map <-
    rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
  if (DEBUG) {
    cat(file = stderr(), length(map))
  }
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmapWarning")
    removeNotification(id = "heatmapNotFound")
  }
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

  expression <- colSums(as.matrix(exprs(gbm[map, ])) >= minExpr)



  projections <- cbind(projections, coExpVal = expression)
  # if(class(projections[,dbCluster])=="factor"){
  p1 <-
    ggplot(projections, aes_string(factor(projections[, dbCluster]), "coExpVal", fill = factor(projections[, dbCluster]))) +
    geom_violin(scale = "width") +
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
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 16),
      legend.position = "none"
    ) +
    xlab(dbCluster) +
    ylab("number genes from list")
  # }else{
  #   return(NULL)
  # }
  if (DEBUG) {
    cat(file = stderr(), "output$gene_vio_plot:done\n")
  }
  return(p1)
}
