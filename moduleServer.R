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
  ns <- session$ns
  subsetData <- NULL
  selectedGroupName <- ""
  groupName <- ""
  addToGroupValue <- FALSE
  
  
  # dim1 <- defaultValues[1]
  # dim2 <- defaultValues[2]
  dim1 <- "PC1"
  dim2 <- "PC2"
  dimCol <- "Gene.count"
  divXBy <- "None"
  divYBy <- "None"
  mod_cl1 <- ""
  observe({
    if (DEBUG) cat(file = stderr(), paste0("observe: clusters\n"))
    mod_cl1 <<- input$clusters
  })
  
  observe({
    if (DEBUG) cat(file = stderr(), paste0("observe: dimension_x\n"))
    dim1 <<- input$dimension_x
  })
  observe({
    if (DEBUG) cat(file = stderr(), paste0("observe: dimension_y\n"))
    dim2 <<- input$dimension_y
  })
  observe({
    if (DEBUG) cat(file = stderr(), paste0("observe: dimension_col\n"))
    dimCol <<- input$dimension_col
  })
  observe({
    if (DEBUG) cat(file = stderr(), paste0("observe: devideXBy\n"))
    divXBy <<- input$devideXBy
  })
  observe({
    if (DEBUG) cat(file = stderr(), paste0("observe: devideYBy\n"))
    divYBy <<- input$devideYBy
  })
  # clusterServer - observe input$groupNames ----
  observe({
    if (DEBUG) cat(file = stderr(), "observe input$groupNames \n")
    input$groupNames # dropdown list with names of cell groups
    isolate({
      updateTextInput(
        session = session, inputId = "groupName",
        value = input$groupNames
      )
    })
  })
  
  
  # clusterServer - updateInput ----
  updateInput <- reactive({
    if (DEBUG) cat(file = stderr(), paste0("updateInput\n"))
    tsneData <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(tsneData)) {
      return(NULL)
    }
    
    # Can also set the label and select items
    if (is.null(mod_cl1) || mod_cl1 == "") mod_cl1 = levels(tsneData$dbCluster)
    updateSelectInput(session, "clusters",
                      choices = levels(tsneData$dbCluster)
                      ,
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
  
  # clusterServer - selectedCellNames ----
  selectedCellNames <- reactive({
    start.time <- base::Sys.time()
    if (DEBUG) cat(file = stderr(), "+++cluster: selectedCellNames\n")
    on.exit(
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "selectedCellNames")
    )
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("selectedCellNames", id = "selectedCellNames", duration = NULL)
    }
    
    brushedPs <- event_data("plotly_selected", source = "subset")
    projections <- projections()
    dimY <- input$dimension_y
    dimX <- input$dimension_x
    geneNames <- input$geneIds
    geneNames2 <- input$geneIds2
    scEx_log <- scEx_log()
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
  
  # clusterServer - returnValues ----
  returnValues <- reactiveValues(
    cluster = reactive(input$clusters),
    selectedCells = reactive({
      if (DEBUG) cat(file = stderr(), "reactiveValues: selectedCells.\n")
      start.time <- Sys.time()
      
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
  
  # if (DEBUG) {
  #   cat(file = stderr(), paste("clusterServers", session$ns("clusters"), "\n"))
  # }
  
  # clusterServer - output$clusters ----
  output$clusters <- renderUI({
    if (DEBUG) cat(file = stderr(), paste("2observe: ns(input$clusters)", session$ns(input$clusters), "\n"))
    retVal <- NULL
    projections <- tData()
    upI <- updateInput() # needed to update input of this module
    ns <- session$ns
     if (DEBUG) cat(file = stderr(), paste("2observe: ns(mod_cl1)", ns(mod_cl1), "\n"))
    if (is.null(projections)) {
      HTML("Please load data first")
    } else {
      noOfClusters <- levels(as.factor(projections$dbCluster))
      # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
      retVal <- selectizeInput(
        ns("clusters"),
        label = "Cluster",
        choices = noOfClusters,
        # selected = input$clusters, # not working because of stack, too slow and possible to create infinite loop
        selected = mod_cl1,
        multiple = TRUE
      )
    }
    retVal
  })
  
  # clusterServer - output$clusterPlot ----
  output$clusterPlot <- renderPlotly({
    if (DEBUG) cat(file = stderr(), paste("Module: output$clusterPlot\n"))
    start.time <- base::Sys.time()
    
    # remove any notification on exit that we don't want
    on.exit(
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "clusterPlot")
    )
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("clusterPlot", id = "clusterPlot", duration = NULL)
    }
    
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
    
    # load(file=paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."));DEBUGSAVE=FALSE
    if (is.null(g_id) || nchar(g_id) == 0) {
      g_id <- featureData$symbol
    }
    if (is.null(logx)) logx <- FALSE
    if (is.null(logy)) logy <- FALSE
    if (is.null(divXBy)) divXBy <- "None"
    if (is.null(divYBy)) divYBy <- "None"
    
    updateProjectionsWithUmiCount(
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
    printTimeEnd(start.time, "clusterPlot")
    exportTestValues(clusterPlot = {p1})  
    return(p1)
  })
  
  # observe({
  #   updateTextInput(session = session, inputId = "groupName",
  #                   value = input$groupNames)
  # })
  
  
  # clusterServer - visibleCellNames ----
  visibleCellNames <- reactive({
    if (DEBUG) cat(file = stderr(), "cluster: selectedCellNames\n")
    start.time <- base::Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "visibleCellNames")
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("visibleCellNames", id = "visibleCellNames", duration = NULL)
    }
    
    projections <- projections()
    if (is.null(projections)) {
      return(NULL)
    }
    inpClusters <- input$clusters
    subsetData <- subset(projections, dbCluster %in% inpClusters)

    printTimeEnd(start.time, "visibleCellNames")
    exportTestValues(visibleCellNames = {subsetData})  
    return(subsetData)
  })

  # clusterServer - observe input$addToGroup ----
  observe({
    if (DEBUG) cat(file = stderr(), "observe input$addToGroup \n")
    if (!is.null(input$addToGroup)) {
      if (input$addToGroup) {
        addToGroupValue <<- TRUE
      } else {
        addToGroupValue <<- FALSE
      }
    }
  })
  
  # clusterServer - observe input$changeGroups ----
  observe({
    if (DEBUG) cat(file = stderr(), "observe input$changeGroups \n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "observe input$changeGroups")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "input-changeGroups")
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("observe: changeGroups", id = "input-changeGroups", duration = NULL)
    }
    
    ns <- session$ns
    input$changeGroups # action button
    addToSelection <- addToGroupValue

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
  
  # clusterServer - output$nCellsVisibleSelected ----
  # display the number of cells that belong to the group, but only from the visible ones
  output$nCellsVisibleSelected <- renderText({
    if (DEBUG) cat(file = stderr(), "nCellsVisibleSelected.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "nCellsVisibleSelected")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "nCellsVisibleSelected")
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("nCellsVisibleSelected", id = "nCellsVisibleSelected", duration = NULL)
    }
    
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
    
    exportTestValues(DummyReactive = {retVal})  
    return(retVal)
  })
  
  # clusterServer - output$nCellsAllSelected ----
  # display the number of cells that belong to the group, including the cells from non visible clusters
  output$nCellsAllSelected <- renderText({
    if (DEBUG) cat(file = stderr(), "nCellsAllSelected started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "nCellsAllSelected")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "nCellsAllSelected")
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("nCellsAllSelected", id = "nCellsAllSelected", duration = NULL)
    }

        grpNs <- groupNames$namesDF
    grpN <- make.names(input$groupName)
    val = 0
    if (grpN %in% colnames(grpNs))
      val = sum(grpNs[, grpN])
    retVal <- paste("number of cells in group over all cells", val)
    
    exportTestValues(DummyReactive = {retVal})  
    return(retVal)
  })
  
  
  # clusterServer - output$additionalOptions ----
  output$additionalOptions <- renderUI({
    if (DEBUG) cat(file = stderr(), "additionalOptions started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "additionalOptions")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "additionalOptions")
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("additionalOptions", id = "additionalOptions", duration = NULL)
    }

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
  
  # clusterServer - output$cellSelection ----
  output$cellSelection <- renderText({
    if (DEBUG) cat(file = stderr(), "cellSelection started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "cellSelection")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "cellSelection")
    })
    # show in the app that this is running
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellSelection", id = "cellSelection", duration = NULL)
    }

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
    
    cells.names <- rownames(projections)[subset(brushedPs, curveNumber == 0)$pointNumber + 1]
    retVal <- paste(cells.names, collapse = ", ")
    
    exportTestValues(ClusterCellSelection = { retVal })
    return(retVal)
  })
  
  # clusterServer - return ----
  return(reactive({
    returnValues
  }))
}


tableSelectionServer <- function(input, output, session,
                                 dataTab) {
  if (DEBUG) cat(file = stderr(), paste("tableSelectionServer", session$ns("test"), "\n"))
  ns <- session$ns
  modSelectedRows <- c()
  
  output$cellSelection <- renderText({
    if (DEBUG) cat(file = stderr(), "cellSelection\n")
    start.time <- Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        printTimeEnd(start.time, "cellSelection")
        removeNotification(id = "cellSelection")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellSelection", id = "cellSelection", duration = NULL)
    }
    
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
    # prox <- proxy
    allrows <- input$cellNameTable_rows_all
    
    if (DEBUGSAVE) {
      save(
        file = paste0("~/scShinyHubDebug/inputselectAll.RData", collapse = "."),
        list = c(ls(), ls(envir = globalenv()))
      )
    }
    # load(file=paste0("~/scShinyHubDebug/inputselectAll.RData", collapse = "."))
    if (ipSelect) {
      proxy %>% selectRows(allrows)
    } else {
      proxy %>% selectRows(NULL)
    }
  })
  
  # output$columnSelection <- renderUI({
  #   
  # })
  observe({
    if (DEBUG) cat(file = stderr(), "observe input$cellNameTable_rows_selected\n")
    modSelectedRows <<- input$cellNameTable_rows_selected
  })
  
  output$cellNameTable <- renderDT({
    if (DEBUG) cat(file = stderr(), "output$cellNameTable\n")
    start.time <- base::Sys.time()
    on.exit(
      if (!is.null(getDefaultReactiveDomain())) {
        printTimeEnd(start.time, "cellNameTable")
        removeNotification(id = "cellNameTable")
      }
    )
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("cellNameTable", id = "cellNameTable", duration = NULL)
    }
    
    dataTables <- dataTab()
    ns <- session$ns
    reorderCells <- input$reorderCells
    selectedRows <- input$cellNameTable_rows_selected
    
    if (is.null(dataTables)) {
      return(NULL)
    }
    if (DEBUGSAVE) {
      save(
        file = paste0("~/scShinyHubDebug/cellNameTable", "ns", ".RData", collapse = "."),
        list = c(ls(), ls(envir = globalenv()))
      )
    }
    # load(file=paste0("~/scShinyHubDebug/cellNameTable", "ns", ".RData", collapse = "."))
    
    
    maxCol <- min(20, ncol(dataTables))
    if (dim(dataTables)[1] > 1) {
      numericCols <- colnames(dataTables %>% select_if(is.numeric))
      nonNumericCols  <- which(!colnames(dataTables) %in% numericCols) # to keep non numeric columns...
      numericCols  <- which(colnames(dataTables) %in% numericCols)
      if (reorderCells &&  length(selectedRows) > 0) {
        csums <- colSums(dataTables[selectedRows, numericCols])
        cols2disp = numericCols[order(csums, decreasing = TRUE)]
      } else {
        cols2disp = numericCols
      }
      cols2disp = c(nonNumericCols, cols2disp)[1:maxCol]
      dataTables = as.data.frame(dataTables[, cols2disp])
      if (DEBUG) cat(file = stderr(), "cellNameTable: done\n")
      return(DT::datatable(dataTables,
                           rownames = F, 
                           filter = "top",
                           selection = list(mode = 'multiple', selected = modSelectedRows),
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
      if (DEBUG) cat(file = stderr(), "download_cellNameTable started.\n")
      start.time <- base::Sys.time()
      on.exit({
        printTimeEnd(start.time, "download_cellNameTable")
        if (!is.null(getDefaultReactiveDomain()))
          removeNotification(id = "download_cellNameTable")
      })
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("download_cellNameTable", id = "download_cellNameTable", duration = NULL)
      }
      
      dataTables <- dataTab()
      write.csv(dataTables, file)
    }
  )
}


pHeatMapModule <- function(input, output, session,
                           pheatmapList # list with arguments for pheatmap
) {
  if (DEBUG) cat(file = stderr(), paste("pHeatMapModule", session$ns("test"), "\n"))
  ns <- session$ns
  
  outfilePH <- NULL
  
  # pHeatMapModule - updateInput ----
  updateInput <- reactive({
    if (DEBUG) cat(file = stderr(), "updateInput started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "updateInput")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "updateInput")
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("updateInput", id = "updateInput", duration = NULL)
    }
    
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
  
  # pHeatMapModule - pHeatMapPlot ----
  output$pHeatMapPlot <- renderImage({
    if (DEBUG) cat(file = stderr(), "pHeatMapPlot started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "pHeatMapPlot")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "pHeatMapPlot")
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("pHeatMapPlot", id = "pHeatMapPlot", duration = NULL)
    }
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification( id = "pHeatMapPlotWARNING")
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
    if ( nrow(heatmapData$mat) > 100 ) {
      showNotification(
        "more than 100 row in heatmap. This can be very slow to display. Only showing first 100 rows",
        id = "pHeatMapPlotWARNING",
        type = "warning",
        duration = 20
      )
      heatmapData$mat = heatmapData$mat[1:100,]
    }
    
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
  
  # pHeatMapModule - additionalOptions ----
  output$additionalOptions <- renderUI({
    if (DEBUG) cat(file = stderr(), "additionalOptions started.\n")
    start.time <- base::Sys.time()
    on.exit({
      printTimeEnd(start.time, "additionalOptions")
      if (!is.null(getDefaultReactiveDomain()))
        removeNotification(id = "additionalOptions")
    })
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("additionalOptions", id = "additionalOptions", duration = NULL)
    }

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

  # pHeatMapModule - download_pHeatMapUI ----
  output$download_pHeatMapUI <- downloadHandler(
    filename = function() {
      paste("pHeatMap", "data.zip", sep = "_")
    },
    content = function(file) {
      if (DEBUG) cat(file = stderr(), "download_pHeatMapUI started.\n")
      start.time <- base::Sys.time()
      on.exit({
        printTimeEnd(start.time, "download_pHeatMapUI")
        if (!is.null(getDefaultReactiveDomain()))
          removeNotification(id = "download_pHeatMapUI")
      })
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("download_pHeatMapUI", id = "download_pHeatMapUI", duration = NULL)
      }
      
      heatmapData <- pheatmapList()
      addColNames <- input$ColNames
      orderColNames <- input$orderNames
      moreOptions <- input$moreOptions
      groupNs <- groupNames$namesDF
      proje <- projections()
      if (DEBUGSAVE)
        save(file = "~/scShinyHubDebug/download_pHeatMapUI.RData", list = c("outfilePH", ls(), ls(envir = globalenv())))
      # load("~/scShinyHubDebug/download_pHeatMapUI.RData")
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
