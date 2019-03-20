source("reactives.R")

# myHeavyCalculations = list(c("scaterPNG", "scaterPNG"))

# Expression ------------------------------------------------------------------
expCluster <- callModule(clusterServer, "expclusters", projections, reactive(input$gene_id))

# these observes should be independant of projections since this will be then executed by default for any changes
updateInputx4 <- reactive({
  tsneData <- projections()

  # Can use character(0) to remove all choices
  if (is.null(tsneData)) {
    return(NULL)
  }

  # Can also set the label and select items
  updateSelectInput(session, "dimension_x4",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[1]
  )

  # Can also set the label and select items
  updateSelectInput(session, "dimension_y4",
    choices = colnames(tsneData),
    selected = colnames(tsneData)[2]
  )
  return(TRUE)
})

output$NumberOfGenesInclude <- renderText({
  idx <- scGeneIdxInclude()
  paste("Number of genes to be included: ", length(idx))
})

output$NumberOfGenesExclude <- renderText({
  idx <- scGeneIdxExclude()
  paste("Number of genes to be included: ", length(idx))
})


# EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
# TODO module for violin plot  ??
output$gene_vio_plot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output$gene_vio_plot\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData <- featureDataReact()
  # log2cpm = log2cpm()
  scEx_log <- scEx_log()
  projections <- projections()
  g_id <- input$gene_id

  if (is.null(featureData) | is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$gene_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gene_vio_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gene_vio_plot.RData")

  geneid <- geneName2Index(g_id, featureData)

  # geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
  #                                        toupper(input$gene_id)), ])[1]

  # expression <- exprs(scEx_log)[geneid, ]
  if (length(geneid) == 1) {
    expression <- assays(scEx_log)[[1]][geneid, ]
  } else {
    expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, ])
  }

  validate(need(is.na(sum(expression)) != TRUE, ""))

  projections <- cbind(projections, expression)
  names(projections)[length(projections)] <- "values"

  p1 <-
    ggplot(projections, aes(factor(dbCluster), values, fill = factor(dbCluster))) +
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
        angle = 90,
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
    xlab("Cluster") +
    ylab("Expression") +
    ggtitle(paste(toupper(featureData[geneid, "Associated.Gene.Name"]), collapse = ", "))
  if (DEBUG) cat(file = stderr(), "output$gene_vio_plot:done\n")
  return(p1)
  # })
})

# EXPLORE TABL DOWNLOAD SELECTED WITH BRUSH ------------------------------------------------------------------
# TODO move to were it belongs
# TODO module for download?
# TODO either integrate in module or use input$bi as selected cell names as a parameter/reactive/global variable.
#      this function will not work as expected as of now
output$downloadExpression <- downloadHandler(
  filename = function() {
    paste(input$cluster, "Selected_Expression_table.csv", sep = "_")
  },
  content = function(file) {
    featureData <- featureDataReact()
    # log2cpm = log2cpm()
    scEx_log <- scEx_log()
    projections <- projections()
    if (is.null(featureData) | is.null(scEx_log) | is.null(projections)) {
      return(NULL)
    }
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/downloadExpression.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/downloadExpression.RData")
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
      toupper(input$gene_id)), ])[1]

    expression <- assays(scEx_log)[[1]][geneid, ]
    # cat(stderr(),colnames(expression)[1:5])
    projections <- cbind(projections, t(expression))
    # if(DEBUG)cat(file=stderr(),grep('^T_',rownames(projections)))

    names(projections)[names(projections) == geneid] <- "values"

    # if(DEBUG)cat(file=stderr(),grep('^T_',rownames(projections)))

    subsetData <- subset(projections, dbCluster == input$cluster)
    # if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
    cells.names <- brushedPoints(subsetData, input$b1, allRows = T)
    # if(DEBUG)cat(file=stderr(),colnames(cells.names))
    cells <-
      rownames(subsetData[which(cells.names$selected_ == TRUE), ])
    # if(DEBUG)cat(file=stderr(),cells[1:5])

    if (length(cells) == 1) {
      subsetExpression <- assays(scEx_log)[[1]][, cells]
      subsetExpression <-
        as.data.frame(subsetExpression, row.names = rownames(scEx_log))
      colnames(subsetExpression) <- cells
      subsetExpression$Associated.Gene.Name <-
        featureData[rownames(subsetExpression), "Associated.Gene.Name"]
      write.csv(subsetExpression, file)
    }
    else {
      subsetExpression <- assays(scEx_log)[[1]][, cells]
      # cat(stderr(),colnames(subsetExpression)[1:5])

      subsetExpression$Associated.Gene.Name <-
        featureData[rownames(subsetExpression), "Associated.Gene.Name"]
      # cat(stderr(),colnames(subsetExpression))
      write.csv(subsetExpression, file)
    }
  }
)

##############################
### Panel Plot
# TODO as module
# data expression panel plot
output$clusters4 <- renderUI({
  if (DEBUG) cat(file = stderr(), "output$clusters4\n")
  projections <- projections()
  upI <- updateInputx4()
  if (is.null(projections)) {
    HTML("Please load data firts")
  } else {
    noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectInput(
      "clusters4",
      label = "Cluster",
      choices = c(c("All"), c(1:noOfClusters)),
      selected = "All"
    )
  }
})

output$panelPlot <- renderPlot({
  if (DEBUG) cat(file = stderr(), "output$panelPlot\n")

  featureData <- featureDataReact()
  # log2cpm = log2cpm()
  scEx_log <- scEx_log()
  scEx <- scEx()
  projections <- projections()
  if (is.null(featureData) | is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }

  genesin <- input$panelplotids
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")
  genesin <- genesin[[1]]
  genesin <- genesin[which(genesin %in% featureData$Associated.Gene.Name)]
  cl4 <- input$clusters4
  dimx4 <- input$dimension_x4
  dimy4 <- input$dimension_y4
  
  if (is.null(cl4)) return(NULL)
  
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/panelPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/panelPlot.RData")

  if (DEBUG) cat(file = stderr(), length(genesin))
  par(mfrow = c(ceiling(length(genesin) / 4), 4), mai = c(0., .3, .3, .3))
  rbPal <- colorRampPalette(c("#f0f0f0", "red"))
  if (DEBUG) cat(file = stderr(), cl4)
  ylim <- c(min(projections[, dimy4]), max(projections[, dimy4]))
  if (class(projections[, dimx4]) == "factor" & dimy4 == "UMI.count") {
    ymax <- 0
    for (i in 1:length(genesin)) {
      geneIdx <- which(featureData$Associated.Gene.Name == genesin[i])
      ymax <- max(ymax, max(Matrix::colSums(assays(scEx)[["counts"]][geneIdx, , drop = FALSE])))
    }
    ylim <- c(0, ymax)
  }
  if (cl4 == "All") {
    for (i in 1:length(genesin)) {
      geneIdx <- which(featureData$Associated.Gene.Name == genesin[i])
      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              assays(scEx_log)[[1]][
                rownames(featureData[geneIdx, ]),
              ]
            ),
            breaks = 10
          )
        )
      ]
      if (class(projections[, dimx4]) == "factor" & dimy4 == "UMI.count") {
        projections[, dimy4] <- Matrix::colSums(assays(scEx)[["counts"]][geneIdx, , drop = FALSE])
      }

      plot(projections[, dimx4], projections[, dimy4],
        col = Col, pch = 16, frame.plot = TRUE, ann = FALSE, ylim = ylim
      )
      title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
      if (DEBUG) cat(file = stderr(), genesin[i])
    }
  } else {
    for (i in 1:length(genesin)) {
      geneIdx <- which(featureData$Associated.Gene.Name == genesin[i])
      subsetTSNE <- subset(projections, dbCluster == cl4)

      Col <- rbPal(10)[
        as.numeric(
          cut(
            as.numeric(
              assays(scEx_log)[[1]][
                rownames(featureData[geneIdx, ]),
              ]
            ),
            breaks = 10
          )
        )
      ]

      names(Col) <- rownames(projections)
      plotCol <- Col[rownames(subsetTSNE)]
      if (class(projections[, dimx4]) == "factor" & dimy4 == "UMI.count") {
        projections[, dimy4] <- Matrix::colSums(assays(scEx)[["counts"]][geneIdx, , drop = FALSE])
        subsetTSNE <- subset(projections, dbCluster == cl4)
      }

      plot(subsetTSNE[, dimx4], subsetTSNE[, dimy4],
        col = plotCol, pch = 16, frame.plot = TRUE,
        ann = FALSE, ylim = ylim
      )
      title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
      if (DEBUG) cat(file = stderr(), cl4)
    }
  }
})


##############################
### Scater QC


output$scaterQC <- renderImage({
  if (DEBUG) cat(file = stderr(), "output$scaterQC\n")
  scaterReads <- scaterReads()
  if (is.null(scaterReads)) {
    return(NULL)
  }

  scaterPNG()
})

output$tsne_plt <- renderPlotly({
  if (DEBUG) cat(file = stderr(), "output$tsne_plt\n")
  # if (v$doPlot == FALSE)
  #   return()
  featureData <- featureDataReact()
  # log2cpm = log2cpm()
  scEx_log <- scEx_log()
  g_id <- input$gene_id
  projections <- projections()

  if (is.null(featureData) | is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne_plt.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/tsne_plt.RData")


  geneid <- geneName2Index(g_id, featureData)

  if (length(geneid) == 1) {
    expression <- assays(scEx_log)[[1]][geneid, ]
  } else {
    expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, ])
  }

  # expression <- log2cpm[geneid, ]
  # cat(file = stderr(), rownames(expression))

  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or gene not expressed"
  ))

  projections <- cbind(projections, expression)
  names(projections)[ncol(projections)] <- "values"

  p <-
    plot_ly(
      projections,
      x = ~tsne1,
      y = ~tsne2,
      z = ~tsne3,
      type = "scatter3d",
      hoverinfo = "text",
      text = paste("Cluster:", as.numeric(as.character(projections$dbCluster))),
      mode = "markers",
      marker = list(
        size = 2,
        line = list(width = 0),
        color = ~values,
        colors = "Greens"
      )
    )
  layout(p, title = paste(toupper(featureData[geneid, "Associated.Gene.Name"]), collapse = ", "))
  # })
})
