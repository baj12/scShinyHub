source("reactives.R")

# since scaterPNG is not used frequently it is not included in the heavyCalculations
# list
# myHeavyCalculations = list(c("scaterPNG", "scaterPNG"))

# Expression ------------------------------------------------------------------
expCluster <- callModule(
  clusterServer,
  "expclusters",
  projections,
  reactive(input$gene_id)
)

# updateInputExpPanel ----
#' updateInputExpPanel
#' update x/y coordinates that can be chosen based on available
#' projections
updateInputExpPanel <- reactive({
  projections <- projections()

  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }

  # Can also set the label and select items
  updateSelectInput(session, "de_dim_x",
    choices = colnames(projections),
    selected = colnames(projections)[1]
  )

  # Can also set the label and select items
  updateSelectInput(session, "dimension_y4",
    choices = colnames(projections),
    selected = colnames(projections)[2]
  )
  return(TRUE)
})

de_X1 <<- "tsne1"
observe({
  de_X1 <<- input$de_dim_x
})





# EXPLORE TAB VIOLIN PLOT ----
# TODO module for violin plot  ??
output$gene_vio_plot <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "gene_vio_plot")
    }
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("gene_vio_plot", id = "gene_vio_plot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$gene_vio_plot\n")

  scEx_log <- scEx_log()
  projections <- projections()
  g_id <- input$gene_id

  if (is.null(scEx_log) | is.null(projections)) {
    if (DEBUG) cat(file = stderr(), "output$gene_vio_plot:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gene_vio_plot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gene_vio_plot.RData")


  featureData <- rowData(scEx_log)
  geneid <- geneName2Index(g_id, featureData)

  if (length(geneid) == 0) {
    return(NULL)
  }
  # geneid <- rownames(featureData[which(featureData$symbol ==
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
    ggtitle(paste(toupper(featureData[geneid, "symbol"]), collapse = ", "))

  printTimeEnd(start.time, "gene_vio_plot")
  exportTestValues(gene_vio_plot = {
    p1
  })
  return(p1)
})


### Panel Plot ----
#' clusterSelectionPanelPlot
#' update selection options for clusters
#' since we allow also "all" we have to have a different strategy for the input
#' it is debateable whether this is usefull to have a different strategy, but for now
#' we leave it as it.
output$clusterSelectionPanelPlot <- renderUI({
  if (DEBUG) cat(file = stderr(), "output$clusterSelectionPanelPlot\n")
  projections <- projections()
  upI <- updateInputExpPanel()
  if (is.null(projections)) {
    HTML("Please load data")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    selectInput(
      "clusterSelectionPanelPlot",
      label = "Cluster",
      choices = c(c("All"), noOfClusters),
      selected = "All"
    )
  }
})

# panelPlot ----
#' panelPlot
#' plot multiple panels for a given list of genes
#' If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to 
#' the count for that gene. Otherwise, all genes are used.
#' normalized counts are used for plotting
output$panelPlot <- renderPlot({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "panelPlot")
  )
  # show in the app that this is running
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("panelPlot", id = "panelPlot", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$panelPlot\n")
  
  scEx_log <- scEx_log()
  projections <- projections()
  genesin <- input$panelplotids
  cl4 <- input$clusterSelectionPanelPlot
  dimx4 <- input$de_dim_x
  dimy4 <- input$dimension_y4
  
  if (is.null(scEx_log) | is.null(projections) | is.null(cl4)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/panelPlot.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/panelPlot.RData")
  
  genesin <- toupper(genesin)
  genesin <- gsub(" ", "", genesin, fixed = TRUE)
  genesin <- strsplit(genesin, ",")
  genesin <- genesin[[1]]
 
  featureData <- rowData(scEx_log)
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
  
  printTimeEnd(start.time, "panelPlot")
  exportTestValues(panelPlot = {ls()})  
  
})


# Scater QC ----
output$scaterQC <- renderImage({
  if (DEBUG) cat(file = stderr(), "output$scaterQC\n")
  scaterReads <- scaterReads()
  if (is.null(scaterReads)) {
    return(NULL)
  }

  scaterPNG()
})

# tsne_plt ----
# tSNE plot within Data exploration - Expressoin
output$tsne_plt <- renderPlotly({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "tsne_plt")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne_plt", id = "tsne_plt", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "output$tsne_plt\n")
  
  scEx_log <- scEx_log()
  g_id <- input$gene_id
  projections <- projections()

  if (is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne_plt.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/tsne_plt.RData")

  retVal <- dataExpltSNEPlot(scEx_log, g_id, projections)

  printTimeEnd(start.time, "dataExpltSNEPlot")
  exportTestValues(dataExpltSNEPlot = {str(retVal)})  
  return(retVal)
})


