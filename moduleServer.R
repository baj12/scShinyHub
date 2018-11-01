source("reactives.R")

#' clusterServer
#'
#' server side shiny module function for printing a 2D represenation of cells
#' it uses the global projections object for plotting

#' @param gene_id name of the gene to be plotted (comma separated list, will be set to upper case)
#' @param tData projections, dataframe with cluster numbers
#' @param DEBUG whether or not to plot debugging messages on stderr
#' @param selectedCells cells that should be marked as a triangle
#' @param legend.position "none", ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#'
#' uses global reactives featueDataRact
#'                       log2cpm # for defining the color of the cells





#
# TODO parameter gene_id should be able to handle multiple gene_ids
# TODO coloring based on number of genes selected (if NULL=color by cluster NR)
#      handle coloring for binarized plot
# TODO potentially integrate gene_id selection into module?
clusterServer <- function(input, output, session,
                          tData,
                          gene_id,
                          # selectedCells  = NULL,
                          legend.position = "none") {
  ns <- session$ns
  subsetData <- NULL
  selectedGroupName <- ""
  groupName <- ""
  
  updateInput <- reactive({
    tsneData <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(tsneData)) {
      return(NULL)
    }
    
    # Can also set the label and select items
    updateSelectInput(session, "dimension_x",
                      choices = colnames(tsneData),
                      selected = colnames(tsneData)[1]
    )
    
    updateSelectInput(session, "dimension_y",
                      choices = colnames(tsneData),
                      selected = colnames(tsneData)[2]
    )
  })
  
  
  
  returnValues <- reactiveValues(
    cluster = reactive(input$clusters),
    # cellNames = ifelse(is.null(subsetData),
    #                    NULL,
    #                    rownames(brushedPoints(subsetData, reactive(input$b1)))
    #                    ),
    selectedCells = reactive({
      if (DEBUG) {
        cat(file = stderr(), paste("clusterServers selectedCells\n"))
      }
      retVal <- rownames(selectedCellNames())
      if (length(retVal) == 0) {
        cat(file = stderr(), paste("selectedCellNames is null\n"))
        retVal <- NULL
      }
      grpN <- make.names(input$groupName)
      grpNs <- groupNames$namesDF
      if (length(grpN) == 0 | length(grpNs) == 0) {
        return(retVal)
      }
      inpClusters <- input$clusters
      projections <- projections()
      dimY <- input$dimension_y
      dimX <- input$dimension_x
      geneNames <- input$geneIds
      featureData <- featureDataReact()
      gbm <- gbm()
      
      if (DEBUGSAVE) {
        cat(file = stderr(), paste("selectedCell: saving\n"))
        save(file = "~/scShinyHubDebug/clusterServerreturnValues.RData", list = c(ls(), ls(envir = globalenv())))
      }
      # load(file="~/scShinyHubDebug/clusterServerreturnValues.RData")
      if (!is.null(projections)) {
        projections <- updateProjectionsWithUmiCount(dimY, dimX, geneNames, featureData, gbm, projections)
        
        subsetData <- subset(projections, dbCluster %in% inpClusters)
        grpSubset <- grpNs[rownames(subsetData), ]
        grpVal <- rownames(grpSubset[grpSubset[, grpN], ])
        if (length(grpVal) > 0) {
          return(grpVal)
        }
      }
      # subsetData <-
      #   subset(projections, as.numeric(as.character(projections$dbCluster)) %in% scCL)
      # cells.1 <- rownames(brushedPoints(subsetData, scBP))
      
      return(retVal)
    })
  )
  
  if (DEBUG) {
    cat(file = stderr(), paste("clusterServers", session$ns("clusters"), "\n"))
  }
  
  output$clusters <- renderUI({
    si <- NULL
    projections <- tData()
    upI <- updateInput() # needed to update input of this module
    ns <- session$ns
    if (is.null(projections)) {
      HTML("Please load data first")
    } else {
      noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
      si <- selectizeInput(
        ns("clusters"),
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0,
        multiple = TRUE
      )
    }
    si
  })
  
  output$clusterPlot <- renderPlot({
    if (DEBUG) {
      cat(file = stderr(), paste("Module: output$clusterPlot", session$ns(input$clusters), "\n"))
    }
    featureData <- featureDataReact()
    log2cpm <- log2cpm()
    gbm <- gbm()
    gbm_log <- gbm_log()
    projections <- tData()
    grpNs <- groupNames$namesDF
    grpN <- make.names(input$groupName)
    
    returnValues$cluster <- input$clusters
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    clId <- input$clusters
    g_id <- gene_id()
    geneNames <- input$geneIds
    
    if (is.null(featureData) | is.null(log2cpm) | is.null(gbm) | is.null(gbm_log) | is.null(projections) | is.null(g_id) | nchar(g_id) == 0) {
      if (DEBUG) cat(file = stderr(), paste("output$clusterPlot:NULL\n"))
      return(NULL)
    }
    
    if (DEBUGSAVE) {
      save(
        file = paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."),
        list = c(ls(envir = globalenv()), ls(), "legend.position")
      )
    }
    # load(file=paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."))
    
    geneid <- geneName2Index(g_id, featureData)
    
    if (length(geneid) == 1) {
      expression <- exprs(gbm_log)[geneid, ]
    } else {
      expression <- Matrix::colSums(exprs(gbm_log)[geneid, ])
    }
    validate(need(is.na(sum(expression)) != TRUE, ""))
    
    geneid <- geneName2Index(geneNames, featureData)
    projections <- updateProjectionsWithUmiCount(dimY, dimX, geneNames, featureData, gbm, projections)
    
    
    projections <- cbind(projections, t(expression))
    names(projections)[ncol(projections)] <- "exprs"
    
    if (DEBUG) {
      cat(file = stderr(), paste("output$dge_plot1:---", ns(clId), "---\n"))
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
  })
  
  # observe({
  #   updateTextInput(session = session, inputId = "groupName",
  #                   value = input$groupNames)
  # })
  
  
  observe({
    input$groupNames # dropdown list with names of cell groups
    isolate({
      updateTextInput(
        session = session, inputId = "groupName",
        value = input$groupNames
      )
    })
  })
  
  visibleCellNames <- reactive({
    if (DEBUG) {
      cat(file = stderr(), "cluster: selectedCellNames\n")
    }
    projections <- projections()
    if (is.null(projections)) {
      return(NULL)
    }
    inpClusters <- input$clusters
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    # if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
    return(subsetData)
  })
  
  
  selectedCellNames <- reactive({
    if (DEBUG) {
      cat(file = stderr(), "cluster: selectedCellNames\n")
    }
    brushedPs <- input$b1
    projections <- projections()
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    geneNames <- input$geneIds
    featureData <- featureDataReact()
    gbm <- gbm()
    if (is.null(projections)) {
      return(NULL)
    }
    inpClusters <- input$clusters
    
    if (DEBUGSAVE) {
      cat(file = stderr(), "cluster: selectedCellNames: saving\n")
      save(file = "~/scShinyHubDebug/selectedCellNames.RData", list = c(ls(), "legend.position", ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/selectedCellNames.RData")
    
    geneid <- geneName2Index(geneNames, featureData)
    projections <- updateProjectionsWithUmiCount(dimY, dimX, geneNames, featureData, gbm, projections)
    
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    # if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
    cells.names <- brushedPoints(subsetData, brushedPs)
    return(cells.names)
  })
  
  addToGroupValue <- FALSE
  
  observe({
    if (DEBUG) {
      cat(file = stderr(), "input$addToGroup changed\n")
    }
    if (!is.null(input$addToGroup)) {
      if (input$addToGroup) {
        addToGroupValue <<- TRUE
      } else {
        addToGroupValue <<- FALSE
      }
    }
  })
  
  observe({
    ns <- session$ns
    input$changeGroups # action button
    addToSelection <- addToGroupValue
    if (DEBUG) {
      cat(file = stderr(), "cluster: changeGroups\n")
    }
    
    isolate({
      brushedPs <- input$b1
      gbm <- gbm()
      inpClusters <- input$clusters
      grpN <- make.names(input$groupName)
      grpNs <- groupNames$namesDF
      cells.names <- selectedCellNames()
      visibleCells <- visibleCellNames()
      if (ncol(grpNs) == 0) {
        initializeGroupNames()
        grpNs <- groupNames$namesDF
      }
      if (is.null(gbm)) {
        return(NULL)
      }
    })
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/changeGroups.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/changeGroups.RData")
    if (!grpN %in% colnames(grpNs)) {
      grpNs[, grpN] <- FALSE
    }
    if (!addToSelection) {
      grpNs[rownames(visibleCells), grpN] <- FALSE
    }
    grpNs[rownames(cells.names), grpN] <- TRUE
    groupNames$namesDF <- grpNs
    updateSelectInput(session, ns("groupNames"),
                      choices = colnames(grpNs),
                      selected = grpN
    )
    updateTextInput(
      session = session, inputId = "groupName",
      value = grpN
    )
    
    selectedGroupName <<- grpN
  })
  
  # display the number of cells that belong to the group, but only from the visible ones
  output$nCellsVisibleSelected <- renderText({
    grpN <- make.names(input$groupName)
    grpNs <- groupNames$namesDF
    inpClusters <- input$clusters
    projections <- projections()
    if (is.null(projections)) {
      return(NULL)
    }
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/nCellsVisibleSelected.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/nCellsVisibleSelected.RData")
    
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    retVal <- paste("Number of visible cells in section", sum(grpNs[rownames(subsetData), grpN]))
    return(retVal)
  })
  
  # display the number of cells that belong to the group, including the cells from non visible clusters
  output$nCellsAllSelected <- renderText({
    grpNs <- groupNames$namesDF
    grpN <- make.names(input$groupName)
    retVal <- paste("number of cells in group over all cells", sum(grpNs[, grpN]))
    return(retVal)
  })
  
  
  output$additionalOptions <- renderUI({
    if (DEBUG) cat(file = stderr(), "cluster: additionalOptions\n")
    ns <- session$ns
    moreOptions <- input$moreOptions
    groupNs <- groupNames$namesDF
    if (!moreOptions) {
      return("")
    }
    
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/additionalOptions.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/additionalOptions.RData")
    
    tagList(
      checkboxInput(ns("addToGroup"), "Add to group/otherwise overwrite", addToGroupValue),
      textInput(ns(id = "groupName"), label = "name group", value = groupName),
      selectInput(
        ns("groupNames"),
        label = "group names",
        choices = colnames(groupNs),
        selected = selectedGroupName
      ),
      verbatimTextOutput(ns("nCellsVisibleSelected")),
      verbatimTextOutput(ns("nCellsAllSelected")),
      actionButton(ns("changeGroups"), "change current selection"),
      checkboxInput(ns("showCells"), "show cell names", FALSE),
      verbatimTextOutput(ns("cellSelection"))
    )
  })
  
  output$cellSelection <- renderText({
    if (DEBUG) cat(file = stderr(), "cluster: cellSelection\n")
    ns <- session$ns
    brushedPs <- (input$b1)
    projections <- projections()
    inpClusters <- (input$clusters)
    myshowCells <- (input$showCells)
    geneNames <- input$geneIds
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    featureData <- featureDataReact()
    gbm <- gbm()
    
    if (!myshowCells) {
      return("")
    }
    if (is.null(projections)) {
      return("")
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cluser cell Selection", id = "clustercellSelection", duration = NULL)
    }
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/clustercellSelection", list = c(ls(envir = globalenv()), ls()))
    }
    # load(file=paste0("~/scShinyHubDebug/clustercellSelection", "ns", ".RData", collapse = "."))
    # load(file=paste0("~/scShinyHubDebug/clustercellSelection"))
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    geneid <- geneName2Index(geneNames, featureData)
    subsetData <- updateProjectionsWithUmiCount(dimY, dimX, geneNames = geneNames, featureData, gbm, subsetData)
    
    cat(file = stderr(), paste(brushedPs$xmin, brushedPs$xmax, "\n"))
    for (axis in c("x", "y")) {
      if (class(subsetData[, brushedPs$mapping[axis][[1]]]) == "factor") {
        subsetData[, brushedPs$mapping[axis][[1]]] <- as.numeric(droplevels(subsetData[, brushedPs$mapping[axis][[1]]]))
      }
    }
    cat(file = stderr(), "cluster: cellSelection\n")
    cells.names <- brushedPoints(subsetData, brushedPs)
    retVal <- paste(rownames(cells.names), collapse = ", ")
    if (DEBUG) cat(file = stderr(), "cluster: cellSelection: done\n")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "clustercellSelection")
    }
    return(retVal)
  })
  
  return(reactive({
    returnValues
  }))
}


tableSelectionServer <- function(input, output, session,
                                 dataTab) {
  if (DEBUG) cat(file = stderr(), paste("tableSelectionServer", session$ns("test"), "\n"))
  ns <- session$ns
  
  output$cellSelection <- renderText({
    if (DEBUG) cat(file = stderr(), "cellSelection\n")
    ns <- session$ns
    
    dataTables <- dataTab()
    selectedRows <- input$cellNameTable_rows_selected
    if (is.null(dataTables)) {
      return(NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellSelection", id = "cellSelection", duration = NULL)
    }
    if (DEBUGSAVE) {
      save(
        file = paste0("~/scShinyHubDebug/cellSelection", "ns", ".RData", collapse = "."),
        list = c(ls(), ls(envir = globalenv()))
      )
    }
    # load(file=paste0("~/scShinyHubDebug/cellSelection", "ns", ".RData", collapse = "."))
    
    # in case there is a table with multiple same row ids (see crPrioGenesTable) the gene names has "___" appended plus a number
    # remove this here
    if (length(selectedRows) > 0) {
      retVal <- rownames(dataTables[selectedRows, ])
      retVal <- sub("(.*)___.*", "\\1", retVal)
      retVal <- unique(retVal)
      retVal <- paste0(retVal, collapse = ", ")
    } else {
      retVal <- NULL
    }
    if (DEBUG) cat(file = stderr(), "cellSelection: done\n")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellSelection")
    }
    return(retVal)
  })
  
  proxy <- dataTableProxy("cellNameTable")
  
  observeEvent(input$selectAll, {
    if (DEBUG) cat(file = stderr(), "input$selectAll\n")
    ipSelect <- input$selectAll
    prox <- proxy
    allrows <- input$cellNameTable_rows_all
    proxy %>% selectRows(NULL)
    if (DEBUGSAVE) {
      save(
        file = paste0("~/scShinyHubDebug/inputselectAll.RData", collapse = "."),
        list = c(ls(), ls(envir = globalenv()))
      )
    }
    # load(file=paste0("~/scShinyHubDebug/inputselectAll.RData", collapse = "."))
    if (ipSelect) {
      proxy %>% selectRows(allrows)
    }
  })
  
  output$columnSelection <- renderUI({
    
  })
  
  output$cellNameTable <- renderDT({
    if (DEBUG) cat(file = stderr(), "output$cellNameTable\n")
    dataTables <- dataTab()
    ns <- session$ns
    
    if (is.null(dataTables)) {
      return(NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellNameTable", id = "cellNameTable", duration = NULL)
    }
    if (DEBUGSAVE) {
      save(
        file = paste0("~/scShinyHubDebug/cellNameTable", "ns", ".RData", collapse = "."),
        list = c(ls(), ls(envir = globalenv()))
      )
    }
    # load(file=paste0("~/scShinyHubDebug/cellNameTable", "", ".RData", collapse = "."))
    
    if (DEBUG) cat(file = stderr(), "cellNameTable: done\n")
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "cellNameTable")
    }
    maxCol <- min(20, ncol(dataTables))
    if (dim(dataTables)[1] > 1) {
      return(DT::datatable(dataTables[, 1:maxCol],
                           rownames = F, filter = "top",
                           options = list(
                             orderClasses = TRUE,
                             autoWidth = TRUE
                           )
      ))
    } else {
      return(NULL)
    }
  })
}

heatmapModuleFunc <- function(
  featureData, 
  gbm_matrix, 
  projections, 
  # genesin, 
  cells
) {
  # genesin <- geneName2Index(genesin, featureData)
  # expression <- gbm_matrix[genesin, cells]
  expression <- gbm_matrix[, cells]
  
  validate(need(
    is.na(sum(expression)) != TRUE,
    "Gene symbol incorrect or genes not expressed"
  ))
  
  projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]
  
  expression <- expression[complete.cases(expression), ]
  
  if (!("sampleNames" %in% colnames(projections))) {
    projections$sample <- 1
  }
  annotation <- data.frame(projections[cells, c("dbCluster", "sampleNames")])
  rownames(annotation) <- colnames(expression)
  colnames(annotation) <- c("Cluster", "sampleNames")
  
  # For high-res displays, this will be greater than 1
  pixelratio <- session$clientData$pixelratio
  if (is.null(pixelratio)) pixelratio <- 1
  width <- session$clientData$output_plot_width
  height <- session$clientData$output_plot_height
  if (is.null(width)) {
    width <- 96 * 7
  } # 7x7 inch output
  if (is.null(height)) {
    height <- 96 * 7
  }
  outfile <- paste0(tempdir(), "/heatmap", base::sample(1:10000, 1), ".png")
  cat(file = stderr(), paste("saving to: ", outfile, "\n"))
  # this can fail with na/inf in hclust error message if there is a row with all the same values
  # med = median(as.vector(as.matrix(expression)))
  # stDev = sd(as.vector(as.matrix(expression)))
  # minBreak = max(0, med - 3* stDev)
  # maxBreak = med + 3* stDev
  # stepBreak = (maxBreak - minBreak) / 6
  nonZeroRows <- which(rowSums(expression) > 0)
  pheatmap(
    as.matrix(expression)[nonZeroRows, order(annotation[, 1], annotation[, 2])],
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
    filename = normalizePath(outfile),
    colorRampPalette(rev(brewer.pal(
      n = 6, name =
        "RdBu"
    )))(6)
  )
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "heatmap")
  }
  return(list(
    src = normalizePath(outfile),
    contentType = "image/png",
    width = width,
    height = height,
    alt = "heatmap should be here"
  ))
}



pHeatMapModule <- function(input, output, session,
                           pheatmapList # list with arguments for pheatmap
) {
  if (DEBUG) cat(file = stderr(), paste("pHeatMapModule", session$ns("test"), "\n"))
  ns <- session$ns
  
  updateInput <- reactive({
    proje <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(proje)) {
      return(NULL)
    }
    
    # Can also set the label and select items
    updateSelectInput(session, "ColNames",
                      choices = colnames(proje),
                      selected = "sampleNames"
    )
    
    # updateSelectInput(session, "dimension_y",
    #                   choices = colnames(proje),
    #                   selected = colnames(proje)[2]
    # )
  })
  
  output$pHeatMapPlot <- renderImage({
    ns <- session$ns
    heatmapData = pheatmapList()
    if (DEBUG) cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot\n")
    # genesin <- ns(input$heatmap_geneids)
    if (is.null(heatmapData)) {
      return(list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "pHeatMapPlot should be here"
      ))
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("pHeatMapPlotModule", id = "pHeatMapPlotModule", duration = NULL)
    }
    
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/pHeatMapPlotModule.RData", list = c(ls(), ls(envir = globalenv()), "heatmapData","input", "output", "session", "pheatmapList", "ns"))
    }
    # load(file = "~/scShinyHubDebug/pHeatMapPlotModule.RData")
    outfile <- paste0(tempdir(), "/heatmap", ns("debug"),base::sample(1:10000, 1), ".png")
    filename = normalizePath(outfile)
    heatmapData$filename = outfile
    do.call(pheatmap, heatmapData)
    
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "pHeatMapPlotModule")
    }
    pixelratio <- session$clientData$pixelratio
    if (is.null(pixelratio)) pixelratio <- 1
    width <- session$clientData$output_plot_width
    height <- session$clientData$output_plot_height
    if (is.null(width)) {
      width <- 96 * 7
    } # 7x7 inch output
    if (is.null(height)) {
      height <- 96 * 7
    }
    
    return(list(src = normalizePath(outfile),
                contentType = 'image/png',
                width = width,
                height = height,
                alt = "heatmap should be here"))
    
  })
  
  output$additionalOptions <- renderUI({
    if (DEBUG) cat(file = stderr(), "heatmapModule: additionalOptions\n")
    ns <- session$ns
    moreOptions <- (input$moreOptions)
    groupNs <- groupNames$namesDF
    if (!moreOptions) {
      return("")
    }
    
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/heatMapadditionalOptions.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/heatMapadditionalOptions.RData")
    
    tagList(
      selectInput(
        ns("ColNames"),
        label = "group names",
        choices = "sampleNames",
        selected = "sampleNames"
      )
    )
  })
  
}