source("reactives.R")

# since DE_scaterPNG is not used frequently it is not included in the heavyCalculations
# list
# myHeavyCalculations = list(c("DE_scaterPNG", "DE_scaterPNG"))

# Expression ------------------------------------------------------------------
callModule(
  clusterServer,
  "DE_expclusters",
  projections,
  reactive(input$DE_gene_id)
)

# DE_updateInputExpPanel ----
#' DE_updateInputExpPanel
#' update x/y coordinates that can be chosen based on available
#' projections

DE_X1 <<- "tsne1"
observe({
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_x\n")
  DE_X1 <<- input$DE_dim_x
})
DE_Y1 <<- "tsne1"
observe({
  if (DEBUG) cat(file = stderr(), "observe: DE_dim_y\n")
  DE_Y1 <<- input$DE_dim_y
})

DE_updateInputExpPanel <- reactive({
  projections <- projections()

  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }

  # Can also set the label and select items
  updateSelectInput(session, "DE_dim_x",
    choices = colnames(projections),
    selected = DE_X1
  )

  # Can also set the label and select items
  updateSelectInput(session, "DE_dim_y",
    choices = colnames(projections),
    selected = DE_Y1
  )
  return(TRUE)
})




# EXPLORE TAB VIOLIN PLOT ----
# TODO module for violin plot  ??
output$DE_gene_vio_plot <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "DE_gene_vio_plot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_gene_vio_plot", id = "DE_gene_vio_plot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_gene_vio_plot\n")

  scEx_log <- scEx_log()
  projections <- projections()
  g_id <- input$DE_gene_id
  ccols <- clusterCols$colPal

  if (is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$DE_gene_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DE_gene_vio_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/DE_gene_vio_plot.RData")


  featureData <- rowData(scEx_log)
  geneid <- geneName2Index(g_id, featureData)

  if (length(geneid) == 0) {
    return(NULL)
  }

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
    scale_color_manual(values = ccols) +
    scale_fill_manual(values = ccols, aesthetics = "fill") +
    
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
    ggtitle(paste(toupper(featureData[geneid, "symbol"]), collapse = ", "))

  printTimeEnd(start.time, "DE_gene_vio_plot")
  exportTestValues(DE_gene_vio_plot = {
    p1
  })
  return(p1)
})


### Panel Plot ----
#' DE_clusterSelectionPanelPlot
#' update selection options for clusters
#' since we allow also "all" we have to have a different strategy for the input
#' it is debateable whether this is usefull to have a different strategy, but for now
#' we leave it as it.
DE_cl1 <<- "All"
observe({
  if (DEBUG) cat(file = stderr(), "observe: DE_clusterSelectionPanelPlot\n")
  DE_cl1 <<- input$DE_clusterSelectionPanelPlot
})
output$DE_clusterSelectionPanelPlot <- renderUI({
  if (DEBUG) cat(file = stderr(), "output$DE_clusterSelectionPanelPlot\n")
  projections <- projections()
  upI <- DE_updateInputExpPanel()
  if (is.null(projections)) {
    HTML("Please load data")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    selectInput(
      "DE_clusterSelectionPanelPlot",
      label = "Cluster",
      choices = c(c("All"), noOfClusters),
      selected = DE_cl1
    )
  }
})

# DE_panelPlot ----
#' DE_panelPlot
#' plot multiple panels for a given list of genes
#' If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to 
#' the count for that gene. Otherwise, all genes are used.
#' normalized counts are used for plotting
output$DE_panelPlot <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_panelPlot")
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_panelPlot", id = "DE_panelPlot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_panelPlot\n")
  
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$DE_panelplotids
  cl4 <- input$DE_clusterSelectionPanelPlot
  dimx4 <- input$DE_dim_x
  dimy4 <- input$DE_dim_y

  if (is.null(scEx_log) | is.null(projections) | is.null(cl4)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DE_panelPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/DE_panelPlot.RData")
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")
  genesin <- genesin[[1]]
 
  featureData <- rowData(scEx_log)
  featureData$symbol = toupper(featureData$symbol)
  genesin <- genesin[which(genesin %in% featureData$symbol)]

  par(mfrow = c(ceiling(length(genesin) / 4), 4), mai = c(0., .3, .3, .3))
  rbPal <- colorRampPalette(c("#f0f0f0", "red"))
  ylim <- c(min(projections[, dimy4]), max(projections[, dimy4]))
  if (class(projections[, dimx4]) == "factor" & dimy4 == "UMI.count") {
    ymax <- 0
    for (i in 1:length(genesin)) {
      geneIdx <- which(featureData$symbol == genesin[i])
      ymax <- max(ymax, max(Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])))
    }
    ylim <- c(0, ymax)
  }
  if (cl4 == "All") {
    for (i in 1:length(genesin)) {
      geneIdx <- which(featureData$symbol == genesin[i])
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
        projections[, dimy4] <- Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])
      }

      plot(projections[, dimx4], projections[, dimy4],
        col = Col, pch = 16, frame.plot = TRUE, ann = FALSE, ylim = ylim
      )
      title(genesin[i], line = -1.2, adj = 0.05, cex.main = 2)
      if (DEBUG) cat(file = stderr(), genesin[i])
    }
  } else {
    for (i in 1:length(genesin)) {
      geneIdx <- which(featureData$symbol == genesin[i])
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
        projections[, dimy4] <- Matrix::colSums(assays(scEx_log)[["logcounts"]][geneIdx, , drop = FALSE])
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
  
  printTimeEnd(start.time, "DE_panelPlot")
  exportTestValues(DE_panelPlot = {ls()})  
  
})


# Scater QC ----
output$DE_scaterQC <- renderImage({
  if (DEBUG) cat(file = stderr(), "output$DE_scaterQC\n")
  scaterReads <- scaterReads()
  if (is.null(scaterReads)) {
    return(NULL)
  }

  DE_scaterPNG()
})

# DE_tsne_plt ----
# tSNE plot within Data exploration - Expressoin
output$DE_tsne_plt <- renderPlotly({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_tsne_plt")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_tsne_plt", id = "DE_tsne_plt", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$DE_tsne_plt\n")
  
  scEx_log <- scEx_log()
  g_id <- input$DE_gene_id
  projections <- projections()

  if (is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DE_tsne_plt.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/DE_tsne_plt.RData")

  retVal <- DE_dataExpltSNEPlot(scEx_log, g_id, projections)

  printTimeEnd(start.time, "DE_dataExpltSNEPlot")
  exportTestValues(DE_dataExpltSNEPlot = {str(retVal)})  
  return(retVal)
})


