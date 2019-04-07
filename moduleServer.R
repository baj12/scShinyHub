source("reactives.R")
library(psych)

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
                          gene_id = returnNull, # reactive
                          # selectedCells  = NULL,
                          legend.position = "none"
                          # ,
                          # defaultValues = c("tsne1", "tsne2")
) {
  # TODO
  ns <- session$ns
  subsetData <- NULL
  selectedGroupName <- ""
  groupName <- ""
  
  # dim1 <- defaultValues[1]
  # dim2 <- defaultValues[2]
  dim1 <- "PC1"
  dim2 <- "PC2"
  dimCol <- "Gene.count"
  divXBy <- "None"
  divYBy <- "None"
  mod_cl1 <- ""
  observe({
    mod_cl1 <<- input$clusters
  })
  
  observe({
    dim1 <<- input$dimension_x
  })
  observe({
    dim2 <<- input$dimension_y
  })
  observe({
    dimCol <<- input$dimension_col
  })
  observe({
    divXBy <<- input$devideXBy
  })
  observe({
    divYBy <<- input$devideYBy
  })
  
  updateInput <- reactive({
    tsneData <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(tsneData)) {
      return(NULL)
    }
    
    # Can also set the label and select items
    if (is.null(mod_cl1) || mod_cl1 == "") mod_cl1 = levels(tsneData$dbCluster)
    updateSelectInput(session, "clusters",
                      choices = levels(tsneData$dbCluster),
                      selected = mod_cl1
    )
    updateSelectInput(session, "dimension_x",
                      choices = colnames(tsneData),
                      selected = dim1
    )
    updateSelectInput(session, "dimension_y",
                      choices = colnames(tsneData),
                      selected = dim2
    )
    updateSelectInput(session, "dimension_col",
                      choices = colnames(tsneData),
                      selected = dimCol
    )
    
    updateSelectInput(session, "devideXBy",
                      choices = c("None",colnames(tsneData)),
                      selected = divXBy
    )
    updateSelectInput(session, "devideYBy",
                      choices = c("None",colnames(tsneData)),
                      selected = divYBy
    )
    
  })
  
  selectedCellNames <- reactive({
    start.time <- Sys.time()
    brushedPs <- event_data("plotly_selected", source = "subset")
    # brushedPs <- input$b1
    projections <- projections()
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    # scEx <- scEx()
    scEx_log <- scEx_log()
    if (DEBUG) {
      cat(file = stderr(), "+++cluster: selectedCellNames\n")
    }
    if (is.null(projections) | is.null(brushedPs)) {
      return(NULL)
    }
    inpClusters <- input$clusters
    
    if (DEBUGSAVE) {
      cat(file = stderr(), "cluster: selectedCellNames: saving\n")
      save(file = "~/scShinyHubDebug/selectedCellNames.RData", list = c(ls(), "legend.position", ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/selectedCellNames.RData")
    
    
    featureData <- rowData(scEx_log)
    geneid <- geneName2Index(geneNames, featureData)
    projections <- updateProjectionsWithUmiCount(
      dimX = dimX, dimY = dimY,
      geneNames = geneNames,
      geneNames2 = geneNames2,
      scEx = scEx_log, projections = projections
    )
    
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    cells.names <- rownames(projections)[subset(brushedPs, curveNumber == 0)$pointNumber + 1]
    
    printTimeEnd(start.time, "selectedCellNames")
    exportTestValues(selectedCellNames = { cells.names })
    return(cells.names)
  })
  
  returnValues <- reactiveValues(
    cluster = reactive(input$clusters),
    selectedCells = reactive({
      start.time <- Sys.time()
      
      if (DEBUG) {
        cat(file = stderr(), paste("clusterServers selectedCells\n"))
      }
      retVal <- selectedCellNames()
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
      geneNames2 <- input$geneIds2
      scEx_log <- scEx_log()
      
      if (DEBUGSAVE) {
        cat(file = stderr(), paste("selectedCell: saving\n"))
        base::save(file = "~/scShinyHubDebug/clusterServerreturnValues.RData", list = c(ls(), ls(envir = globalenv())))
      }
      # load(file="~/scShinyHubDebug/clusterServerreturnValues.RData")
      featureData <- rowData(scEx_log)
      if (!is.null(projections)) {
        projections <- updateProjectionsWithUmiCount(
          dimX = dimX, dimY = dimY,
          geneNames = geneNames,
          geneNames2 = geneNames2,
          scEx = scEx_log, projections = projections
        )
        
        subsetData <- subset(projections, dbCluster %in% inpClusters)
        grpSubset <- grpNs[rownames(subsetData), ]
        grpVal <- rownames(grpSubset[grpSubset[, grpN], ])
        if (length(grpVal) > 0) {
          return(grpVal)
        }
      }
      
      printTimeEnd(start.time, "clusterServerReturnVal")
      exportTestValues(clusterServerReturnVal = { retVal })
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
    cat(file = stderr(), paste("2observe: ns(input$clusters)", session$ns(input$clusters), "\n"))
    cat(file = stderr(), paste("2observe: ns(mod_cl1)", ns(mod_cl1), "\n"))
    if (is.null(projections)) {
      HTML("Please load data first")
    } else {
      noOfClusters <- levels(as.factor(projections$dbCluster))
      # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
      si <- selectizeInput(
        ns("clusters"),
        label = "Cluster",
        choices = noOfClusters,
        # selected = input$clusters, # not working because of stack, too slow and possible to create infinite loop
        selected = mod_cl1,
        multiple = TRUE
      )
    }
    si
  })
  
  output$clusterPlot <- renderPlotly({
    if (DEBUG) {
      cat(file = stderr(), paste("Module: output$clusterPlot", session$ns(input$clusters), "\n"))
    }
    scEx_log <- scEx_log()
    scEx_log <- scEx_log()
    projections <- tData()
    grpNs <- groupNames$namesDF
    grpN <- make.names(input$groupName)
    
    returnValues$cluster <- input$clusters
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    dimCol <- input$dimension_col
    clId <- input$clusters
    g_id <- gene_id()
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    logx <- input$logX
    logy <- input$logY
    divXBy <- input$devideXBy
    divYBy <- input$devideYBy
    scols <- sampleCols$colPal
    ccols <- clusterCols$colPal
    
    
    if (is.null(scEx_log) | is.null(scEx_log) | is.null(projections)) {
      if (DEBUG) cat(file = stderr(), paste("output$clusterPlot:NULL\n"))
      return(NULL)
    }
    
    featureData <- rowData(scEx_log)
    if (DEBUGSAVE) {
      cat(file = stderr(), paste("cluster plot saving\n"))
      save(
        file = paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."),
        list = c(ls(envir = globalenv()), ls(), "legend.position")
      )
      cat(file = stderr(), paste("cluster plot saving done\n"))
    }
    
    # load(file=paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."))
    if (is.null(g_id) || nchar(g_id) == 0) {
      g_id <- featureData$symbol
    }
    if (is.null(logx)) logx <- FALSE
    if (is.null(logy)) logy <- FALSE
    if (is.null(divXBy)) divXBy <- "None"
    if (is.null(divYBy)) divYBy <- "None"
    
    subsetData <- updateProjectionsWithUmiCount(
      dimX = dimX, dimY = dimY,
      geneNames = geneNames,
      geneNames2 = geneNames2,
      scEx = scEx_log, projections = projections
    )
    if (dimCol == "sampleNames") {
      myColors <- scols
    } else {
      myColors <- NULL
    }
    if (dimCol == "dbCluster") {
      myColors <- ccols
    }
    
    p1 <- plot2Dprojection(scEx_log, projections, g_id, featureData, geneNames,
                           geneNames2, dimX, dimY, clId, grpN, legend.position,
                           grpNs = grpNs, logx, logy, divXBy, divYBy, dimCol, colors = myColors
    )
    return(p1)
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
    
    # we isolate here because we only want to change if the button is clicked.
    # TODO what happens if new file is loaded??? => problem!
    isolate({
      # we should react on a changed filename, but that would imply calculating pca's etc directly after loading
      scEx <- scEx()
      prjs <- sessionProjections$prjs
      # brushedPs <- event_data("plotly_selected", source = "subset")
      # scEx <- scEx()
      inpClusters <- input$clusters
      grpN <- make.names(input$groupName)
      grpNs <- groupNames$namesDF
      cells.names <- selectedCellNames()
      visibleCells <- visibleCellNames()
      if (ncol(grpNs) == 0) {
        initializeGroupNames()
        grpNs <- groupNames$namesDF
      }
      if (is.null(scEx)) {
        return(NULL)
      }
    })
    if (length(grpN) == 0) {
      return(NULL)
    }
    if (DEBUGSAVE) {
      cat(file = stderr(), "save: changeGroups\n")
      save(file = "~/scShinyHubDebug/changeGroups.RData", list = c(ls(), ls(envir = globalenv())))
      cat(file = stderr(), "done save: changeGroups\n")
    }
    # load(file="~/scShinyHubDebug/changeGroups.RData")
    if (!grpN %in% colnames(grpNs)) {
      grpNs[, grpN] <- FALSE
    }
    if (!addToSelection) {
      grpNs[rownames(visibleCells), grpN] <- FALSE
    }
    grpNs[cells.names, grpN] <- TRUE
    groupNames$namesDF <- grpNs
    updateSelectInput(session, ns("groupNames"),
                      choices = colnames(grpNs),
                      selected = grpN
    )
    updateTextInput(
      session = session, inputId = "groupName",
      value = grpN
    )
    
    if (ncol(prjs) > 0) {
      # didn't find a way to easily overwrite columns
      for (cn in colnames(grpNs)) {
        if (cn %in% colnames(prjs)) {
          prjs[, cn] <- grpNs[, cn]
        } else {
          prjs <- base::cbind(prjs, grpNs[, cn], deparse.level = 0)
          colnames(prjs)[ncol(prjs)] <- cn
        }
      }
      
      sessionProjections$prjs <- prjs
    } else {
      sessionProjections$prjs <- grpNs
    }
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
    val = 0
    if (grpN %in% colnames(grpNs))
      val = sum(grpNs[, grpN])
    retVal <- paste("number of cells in group over all cells", val)
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
      fluidRow(
        column(
          3,
          checkboxInput(ns("logX"), "log transform X", value = FALSE)
        ),
        column(
          3,
          checkboxInput(ns("logY"), "log transform Y", value = FALSE)
        ),
        column(
          3,
          selectInput(
            ns("devideXBy"),
            label = "Devide X by",
            choices = c("None", "Gene.Count", "UMI.Count"),
            selected = "None"
          )
        ),
        column(
          3,
          selectInput(
            ns("devideYBy"),
            label = "Devide Y by",
            choices = c("None", "Gene.Count", "UMI.Count"),
            selected = "None"
          )
        )
        
        
      ),
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
    start.time <- Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "clustercellSelection")
      }
    )
    if (DEBUG) cat(file = stderr(), "cluster: cellSelection\n")
    ns <- session$ns
    brushedPs <- event_data("plotly_selected", source = "subset")
    projections <- projections()
    inpClusters <- (input$clusters)
    myshowCells <- (input$showCells)
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    scEx_log <- scEx_log()
    
    if (!myshowCells) {
      return("")
    }
    if (is.null(projections)) {
      return("")
    }
    if (is.null(brushedPs)) {
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
    featureData <- rowData(scEx_log)
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    geneid <- geneName2Index(geneNames, featureData)
    subsetData <- updateProjectionsWithUmiCount(
      dimX = dimX, dimY = dimY,
      geneNames = geneNames,
      geneNames2 = geneNames2,
      scEx = scEx_log, projections = projections
    )
    
    # cat(file = stderr(), paste(brushedPs$xmin, brushedPs$xmax, "\n"))
    # for (axis in c("x", "y")) {
    #   if (class(subsetData[, brushedPs$mapping[axis][[1]]]) == "factor") {
    #     subsetData[, brushedPs$mapping[axis][[1]]] <- as.numeric(droplevels(subsetData[, brushedPs$mapping[axis][[1]]]))
    #   }
    # }
    # cells.names <- brushedPoints(subsetData, brushedPs)
    cells.names <- rownames(projections)[subset(brushedPs, curveNumber == 0)$pointNumber + 1]
    
    retVal <- paste(cells.names, collapse = ", ")
    
    printTimeEnd(start.time, "cellSelection")
    exportTestValues(ClusterCellSelection = { retVal })
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
    start.time <- Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "cellSelection")
      }
    )
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
    
    printTimeEnd(start.time, "tableSelectionServer-cellSelection")
    exportTestValues(tableSelectionServercellSelection = { retVal })
    return(retVal)
  })
  
  proxy <- dataTableProxy("cellNameTable")
  
  observeEvent(input$selectAll, {
    if (DEBUG) cat(file = stderr(), "observe input$selectAll\n")
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
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "cellNameTable")
      }
    )
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
    # load(file=paste0("~/scShinyHubDebug/cellNameTable", "ns", ".RData", collapse = "."))
    
    if (DEBUG) cat(file = stderr(), "cellNameTable: done\n")
    maxCol <- min(20, ncol(dataTables))
    dataTables = as.data.frame(dataTables)
    if (dim(dataTables)[1] > 1) {
      return(DT::datatable(dataTables[, 1:maxCol],
                           rownames = F, filter = "top",
                           options = list(
                             orderClasses = TRUE,
                             autoWidth = TRUE
                           )
      ))
    } else {
      return(warning("test"))
    }
  })
  
  output$download_cellNameTable <- downloadHandler(
    filename = function() {
      paste("cellNameTable", "table.csv", sep = "_")
    },
    content = function(file) {
      dataTables <- dataTab()
      write.csv(dataTables, file)
    }
  )
}

# heatmapModuleFunc <- function(
#                               featureData,
#                               scEx_matrix,
#                               projections,
#                               # genesin,
#                               cells) {
#   # genesin <- geneName2Index(genesin, featureData)
#   # expression <- scEx_matrix[genesin, cells]
#   expression <- scEx_matrix[, cells]
# 
#   validate(need(
#     is.na(sum(expression)) != TRUE,
#     "Gene symbol incorrect or genes not expressed"
#   ))
# 
#   projections <- projections[order(as.numeric(as.character(projections$dbCluster))), ]
# 
#   expression <- expression[complete.cases(expression), ]
# 
#   if (!("sampleNames" %in% colnames(projections))) {
#     projections$sample <- 1
#   }
#   annotation <- data.frame(projections[cells, c("dbCluster", "sampleNames")])
#   rownames(annotation) <- colnames(expression)
#   colnames(annotation) <- c("Cluster", "sampleNames")
# 
#   # For high-res displays, this will be greater than 1
#   pixelratio <- session$clientData$pixelratio
#   if (is.null(pixelratio)) pixelratio <- 1
#   width <- session$clientData$output_plot_width
#   height <- session$clientData$output_plot_height
#   if (is.null(width)) {
#     width <- 96 * 7
#   } # 7x7 inch output
#   if (is.null(height)) {
#     height <- 96 * 7
#   }
#   outfile <- paste0(tempdir(), "/heatmap", base::sample(1:10000, 1), ".png")
#   cat(file = stderr(), paste("saving to: ", outfile, "\n"))
#   # this can fail with na/inf in hclust error message if there is a row with all the same values
#   # med = median(as.vector(as.matrix(expression)))
#   # stDev = sd(as.vector(as.matrix(expression)))
#   # minBreak = max(0, med - 3* stDev)
#   # maxBreak = med + 3* stDev
#   # stepBreak = (maxBreak - minBreak) / 6
#   nonZeroRows <- which(rowSums(expression) > 0)
#   mycolPal <- colorRampPalette(brewer.pal(
#     n = 6, name =
#       "RdYlBu"
#   ))(6)
#   if (dbCluster == "sampleNames"){
#     mycolPal = sampleCols$colPal
#   }
#   
#   TRONCO::pheatmap(
#     expression[nonZeroRows, order(annotation[, 1], annotation[, 2])],
#     cluster_rows = TRUE,
#     cluster_cols = TRUE,
#     scale = "row",
#     fontsize_row = 14,
#     labels_col = colnames(expression),
#     labels_row = featureData[rownames(expression), "symbol"],
#     show_rownames = TRUE,
#     annotation_col = annotation,
#     show_colnames = FALSE,
#     annotation_legend = TRUE,
#     # breaks = seq(minBreak, maxBreak, by = stepBreak),
#     # filename = 'test.png',
#     filename = normalizePath(outfile, mustWork = FALSE),
#     color = mycolPal
#   )
#   return(list(
#     src = normalizePath(outfile, mustWork = FALSE),
#     contentType = "image/png",
#     width = width,
#     height = height,
#     alt = "heatmap should be here"
#   ))
# }
# 
# 

pHeatMapModule <- function(input, output, session,
                           pheatmapList # list with arguments for pheatmap
) {
  if (DEBUG) cat(file = stderr(), paste("pHeatMapModule", session$ns("test"), "\n"))
  ns <- session$ns
  
  outfilePH <- NULL
  
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
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        removeNotification(id = "pHeatMapPlotModule")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("pHeatMapPlotModule", id = "pHeatMapPlotModule", duration = NULL)
    }
    ns <- session$ns
    heatmapData <- pheatmapList()
    addColNames <- input$ColNames
    orderColNames <- input$orderNames
    moreOptions <- input$moreOptions
    colTree <- input$showColTree
    
    proje <- projections()
    if (DEBUG) cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot\n")
    if (DEBUGSAVE) {
      cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot saving\n")
      save(file = "~/scShinyHubDebug/pHeatMapPlotModule.RData", list = c(ls(), ls(envir = globalenv()), "heatmapData", "input", "output", "session", "pheatmapList", "ns"))
      cat(file = stderr(), "output$pHeatMapModule:pHeatMapPlot saving done\n")
    }
    # load(file = "~/scShinyHubDebug/pHeatMapPlotModule.RData")
    
    if (is.null(heatmapData) | is.null(proje) | is.null(heatmapData$mat)) {
      return(list(
        src = "empty.png",
        contentType = "image/png",
        width = 96,
        height = 96,
        alt = "pHeatMapPlot should be here"
      ))
    }
    
    outfile <- paste0(tempdir(), "/heatmap", ns("debug"), base::sample(1:10000, 1), ".png")
    outfile <- normalizePath(outfile, mustWork = FALSE)
    heatmapData$filename <- outfile
    
    if (length(addColNames) > 0 & moreOptions) {
      heatmapData$annotation_col <- proje[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
    }
    if (sum(orderColNames %in% colnames(proje)) > 0 & moreOptions) {
      heatmapData$cluster_cols <- FALSE
      colN <- rownames(psych::dfOrder(proje, orderColNames))
      colN <- colN[colN %in% colnames(heatmapData$mat)]
      heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
      # return()
    }
    if (moreOptions) {
      heatmapData$cluster_cols = colTree
    }
    # orgMat = heatmapData$mat
    
    # heatmapData$mat = orgMat
    # system.time(do.call(pheatmap, heatmapData))
    # heatmapData$mat = as(orgMat, "dgTMatrix")
    heatmapData$fontsize <- 14
    # heatmapData$fontsize_row = 18
    # heatmapData$filename=NULL
    system.time(do.call(TRONCO::pheatmap, heatmapData))
    
    pixelratio <- session$clientData$pixelratio
    if (is.null(pixelratio)) pixelratio <- 1
    # width <- session$clientData$output_plot_width
    # height <- session$clientData$output_plot_height
    # if (is.null(width)) {
    #   width <- 96 * 7
    # } # 7x7 inch output
    # if (is.null(height)) {
    #   height <- 96 * 7
    # }
    outfilePH <<- outfile
    return(list(
      src = outfilePH,
      contentType = "image/png",
      width = "100%",
      height = "100%",
      alt = "heatmap should be here"
    ))
  })
  
  output$additionalOptions <- renderUI({
    if (DEBUG) cat(file = stderr(), "heatmapModule: additionalOptions\n")
    ns <- session$ns
    moreOptions <- (input$moreOptions)
    groupNs <- groupNames$namesDF
    proje <- projections()
    if (!moreOptions | is.null(proje)) {
      return("")
    }
    
    
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/heatMapadditionalOptions.RData", list = c(ls(), ls(envir = globalenv())))
    }
    # load(file="~/scShinyHubDebug/heatMapadditionalOptions.RData")
    
    tagList(
      checkboxInput(ns("showColTree"), label = "Show tree for cells", value = FALSE),
      selectInput(
        ns("ColNames"),
        label = "group names",
        choices = colnames(proje),
        selected = "sampleNames",
        multiple = TRUE
      ),
      
      selectInput(
        ns("orderNames"),
        label = "order of columns",
        choices = colnames(proje),
        selected = "",
        multiple = TRUE
      )
    )
  })
  output$download_pHeatMapUI <- downloadHandler(
    filename = function() {
      paste("pHeatMap", "data.zip", sep = "_")
    },
    content = function(file) {
      heatmapData <- pheatmapList()
      addColNames <- input$ColNames
      orderColNames <- input$orderNames
      moreOptions <- input$moreOptions
      groupNs <- groupNames$namesDF
      proje <- projections()
      save(file = "~/scShinyHubDebug/download_pHeatMapUI.RData", list = c("outfilePH", ls(), ls(envir = globalenv())))
      dfilename <- paste0(reportTempDir, "/sessionData.RData")
      base::save(
        file = dfilename, list =
          c("heatmapData", "addColNames", "orderColNames", "moreOptions", "proje", "groupNs")
      )
      
      
      # 2do: not sure why I cannot find the original file...
      # maybe there is an intermediate session created?
      outfile <- paste0(tempdir(), "/heatmap", ns("debug"), base::sample(1:10000, 1), ".png")
      outfile <- normalizePath(outfile, mustWork = FALSE)
      heatmapData$filename <- outfile
      
      if (length(addColNames) > 0 & moreOptions) {
        heatmapData$annotation_col <- proje[rownames(heatmapData$annotation_col), addColNames, drop = FALSE]
      }
      if (sum(orderColNames %in% colnames(proje)) > 0 & moreOptions) {
        heatmapData$cluster_cols <- FALSE
        colN <- rownames(psych::dfOrder(proje, orderColNames))
        colN <- colN[colN %in% colnames(heatmapData$mat)]
        heatmapData$mat <- heatmapData$mat[, colN, drop = FALSE]
      }
      
      do.call(pheatmap, heatmapData)
      
      
      zippedReportFiles <- c(
        dfilename,
        outfile
      )
      zip(file, zippedReportFiles, flags = "-9Xj")
    }
  )
}
