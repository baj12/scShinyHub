#' selectedDge
#' stores table with differentially expressed genes
#' is used to export file and write csv file
selectedDge <- reactiveValues(
  dgeTable = data.frame()
)

<<<<<<< HEAD
#' getCells
#' get list of cells from selctions 
getCells <- function(projections, cl1, db1, db2) {
  dbCluster = projections$dbCluster
=======

# 
dge_func <- function(projections, scEx_log, dbCluster, cl1, db1, db2) {
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
  subsetData <- subset(projections, dbCluster %in% cl1)
  cells.1 <- rownames(shiny::brushedPoints(subsetData, db1))
  cells.2 <- rownames(shiny::brushedPoints(subsetData, db2))
  retVal <- list(c1 = cells.1, c2 = cells.2)
  return(retVal)
}

#' define different methods for calculating diff. expressed genes
#' first entry is displayed in Radio box, second is function to be called.
myDiffExpFunctions = list(
  c("dge as per CellView", "dge_CellViewfunc"),
  c("t-test", "dge_ttest")
)

#' dge_CellViewfunc
#' calculate differentically expressed genes given 2 sets of cells
dge_CellViewfunc <- function(scEx_log, cells.1, cells.2) {
  
  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]),]
  genes.use <- rownames(subsetExpression)
  # expMean exponential mean
  data.1 <- apply(subsetExpression[genes.use, cells.1], 1, expMean)
  data.2 <- apply(subsetExpression[genes.use, cells.2], 1, expMean)
  total.diff <- (data.1 - data.2)
  
  genes.diff <- names(which(abs(total.diff) > .2))
  genes.use <- ainb(genes.use, genes.diff)
  
  retVal <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  retVal[, "avg_diff"] <- total.diff[rownames(retVal)]
  retVal$symbol <-
    featureData[rownames(retVal), "symbol"]
  return(retVal)
}

<<<<<<< HEAD
#' dge_ttest
#' t-test on selected cells
dge_ttest <- function(scEx_log, cells.1, cells.2) {
  featureData <- rowData(scEx_log)
  scEx_log <- as.matrix(assays(scEx_log)[[1]])
  subsetExpression <- scEx_log[complete.cases(scEx_log[, union(cells.1, cells.2)]),]
  genes.use <- rownames(subsetExpression)
  
  p_val <- apply(subsetExpression, 1, function(x) t.test(x[cells.1], x[cells.2])$p.value)
  p_val <- p_val[!is.na(p_val)]
  
  retVal = data.frame(p_val = p_val, symbol = featureData[names(p_val), "symbol"])
  return(retVal)
}

#' dge
#' manage calculation for differential expression analysis
=======
#' dge
#' calculate differential expression analysis
#' based on 
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
dge <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain())) {
      removeNotification(id = "dge")
    }
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("dge", id = "dge", duration = NULL)
  }
  if (!is.null(getDefaultReactiveDomain()))
    removeNotification(id = "dgewarning")
  if (DEBUG) cat(file = stderr(), "dge\n")
<<<<<<< HEAD
=======
  
  scEx_log <- scEx_log()
  prj <- projections()
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
  
  scEx_log <- scEx_log()
  projections <- projections()
  cl1 <- input$dgeClustersSelection
  db1 <- input$db1
  db2 <- input$db2
  method <- input$dgeRadioButton
  
<<<<<<< HEAD
  if (is.null(scEx_log) | is.null(projections)) {
    return(NULL)
  }
=======
  if (is.null(scEx_log) | is.null(prj)) {
    return(NULL)
  }
  
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/dge.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/dge.RData')
  
<<<<<<< HEAD
  methodIdx <- ceiling(which(unlist(diffExpFunctions)== method)/2)
  dgeFunc <- diffExpFunctions[[methodIdx]][2]
  gCells <- getCells(projections, cl1, db1, db2)
  retVal <- do.call(dgeFunc, args = list(scEx_log = scEx_log,
                                         cells.1 = gCells$c1, cells.2 = gCells$c2))
=======
  retVal <- dge_func(
    projections = prj, scEx_log = scEx_log,
    dbCluster = prj$dbCluster, cl1 = cl1, db1 = db1, db2 = db2
  )
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
  
  if (nrow(retVal) == 0) {
    if (DEBUG) cat(file = stderr(), "dge: nothing found\n")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: nothing found", id = "dgewarning", duration = 10, type = "warning")
    }
  }
<<<<<<< HEAD
  # update reactiveValue
  selectedDge$dgeTable <- retVal
  
  printTimeEnd(start.time, "dge")
  exportTestValues(dge = {retVal})  
  return(retVal)
})

# using these global variables allows us to store a set values even when the projections are changing
subClusterDim1 <- "PC1"
subClusterDim2 <- "PC2"
subClusterClusters <<- NULL
observe({
  subClusterDim1 <<- input$subscluster_x1
})
observe({
  subClusterDim2 <<- input$subscluster_y1
})
observe({
  projections <- projections()
  if (!is.null(projections)) {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    if (is.null(subClusterClusters)){
       subClusterClusters <<- noOfClusters
    }
  }
})
observe({
  subClusterClusters <<- input$dgeClustersSelection
})


# subcluster axes ----
# update axes in subcluster analysis
updateInputSubclusterAxes <- reactive({
  projections <- projections()
  # we combine the group names with the projections to add ability to select groups
  # gn <- groupNames$namesDF
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/updateInputSubclusterAxes.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/updateInputSubclusterAxes.RData")
  # if (length(gn) > 0) {
  #   projections <- cbind(projections, gn[rownames(projections), ] * 1)
  # }
  # Can also set the label and select items
  updateSelectInput(session, "subscluster_x1",
                    choices = colnames(projections),
                    selected = subClusterDim1
  )
  
  updateSelectInput(session, "subscluster_y1",
                    choices = colnames(projections),
                    selected = subClusterDim2
  )
})


# sub cluster analysis ( used for 2 panels )
output$dgeClustersSelection <- renderUI({
  projections <- projections()
  up1 <- updateInputSubclusterAxes()
  
  if (DEBUG) cat(file = stderr(), "output$dgeClustersSelection\n")
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/dgeClustersSelection.RData", list = c(ls(envir = globalenv(), ls(), "subClusterClusters")))
  }
  # load(file="~/scShinyHubDebug/dgeClustersSelection.RData")
  
  
  if (is.null(projections)) {
    tags$span(style="color:red", "Please load data first")
  } else {
    noOfClusters <- levels(as.factor(projections$dbCluster))
    # noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectizeInput(
      "dgeClustersSelection",
      label = "Cluster",
      choices = noOfClusters,
      selected = subClusterClusters,
      multiple = TRUE
    )
  }
})
=======
  selectedDge$dgeTable <- retVal
  
  printTimeEnd(start.time, "dge")
  exportTestValues(dge = {retVal})  
  return(retVal)
})

# using these global variables allows us to store a set values even when the projections are changing
subClusterDim1 <- "PC1"
subClusterDim2 <- "PC2"
subClusterClusters <- "1"
observe({
  subClusterDim1 <<- input$subscluster_x1
})
observe({
  subClusterDim2 <<- input$subscluster_y1
})
observe({
  subClusterClusters <<- input$clusters1
})


# subcluster axes ----
# update axes in subcluster analysis
updateInputSubclusterAxes <- reactive({
  projections <- projections()
  # we combine the group names with the projections to add ability to select groups
  # gn <- groupNames$namesDF
  # Can use character(0) to remove all choices
  if (is.null(projections)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/updateInputSubclusterAxes.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/updateInputSubclusterAxes.RData")
  # if (length(gn) > 0) {
  #   projections <- cbind(projections, gn[rownames(projections), ] * 1)
  # }
  # Can also set the label and select items
  updateSelectInput(session, "subscluster_x1",
                    choices = colnames(projections),
                    selected = subClusterDim1
  )
  
  updateSelectInput(session, "subscluster_y1",
                    choices = colnames(projections),
                    selected = subClusterDim2
  )
})


# sub cluster analysis ( used for 2 panels )
output$clusters1 <- renderUI({
  
  projections <- projections()
  up1 <- updateInputSubclusterAxes()
  
  if (DEBUG) cat(file = stderr(), "output$clusters1\n")
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/clusters1_dge.RData", list = c(ls(envir = globalenv(), ls(), "subClusterClusters")))
  }
  # load(file="~/scShinyHubDebug/clusters1_dge.RData")
  
  
  if (is.null(projections)) {
    HTML("Please load data first")
  } else {
    noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectizeInput(
      "clusters1",
      label = "Cluster",
      choices = c(1:noOfClusters),
      selected = subClusterClusters,
      multiple = TRUE
    )
  }
})
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20


#' subCluster2Dplot
#' plots cells in 2D for the subcluster anlaysis. The handling of the selection is done
#' outside this function
subCluster2Dplot <- function() {
  renderPlot({
    if (DEBUG) cat(file = stderr(), "output$dge_plot2\n")
    
    projections <- projections()
    x1 <- input$subscluster_x1
    y1 <- input$subscluster_y1
<<<<<<< HEAD
    c1 <- input$dgeClustersSelection
=======
    c1 <- input$clusters1
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
    
    if (is.null(projections)) {
      return(NULL)
    }
    if (DEBUGSAVE) {
      save(file = "~/scShinyHubDebug/dge_plot2.RData", list = c(ls(envir = globalenv(), ls())))
    }
    # load(file="~/scShinyHubDebug/dge_plot2.RData")
    
    
    subsetData <- subset(projections, dbCluster %in% c1)
    p1 <-
      ggplot(subsetData,
             aes_string(x = x1, y = y1),
             color = "dbCluster"
      ) +
      geom_point(aes(colour = dbCluster)) +
      geom_point(
        shape = 1,
        size = 4,
        aes(colour = dbCluster)
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
        legend.position = "none"
      ) +
      ggtitle(c1)
    p1
  })
}
