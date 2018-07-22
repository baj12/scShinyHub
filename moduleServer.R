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
                          selectedCells  = NULL,
                          legend.position = "none"){
  ns = session$ns 
  subsetData = NULL
  selectedGroupName = ""
  
  updateInput <- reactive({
    tsneData <- projections()
    
    # Can use character(0) to remove all choices
    if (is.null(tsneData)){
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
    brushedPs = reactive(input$b1)
  )
  if(DEBUG)cat(file=stderr(), paste("clusterServers", session$ns("clusters"), "\n"))
  
  output$clusters <- renderUI({
    si=NULL
    projections = tData()
    upI <- updateInput()
    ns = session$ns
    if(is.null(projections)){
      HTML("Please load data first")
    }else{
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
    if(DEBUG)cat(file=stderr(), paste("Module: output$clusterPlot",session$ns(input$clusters), "\n"))
    featureData = featureDataReact()
    log2cpm = log2cpm()
    projections = tData()
    
    returnValues$cluster = input$clusters
    dimY = input$dimension_y
    dimX = input$dimension_x
    clId = input$clusters
    g_id=gene_id()
    
    
    if(is.null(featureData) | is.null(log2cpm) | is.null(projections) | is.null(g_id) | nchar(g_id) == 0){
      if(DEBUG)cat(file=stderr(), paste("output$clusterPlot:NULL\n"))
      return(NULL)
    }
    
    if(DEBUGSAVE) 
      save(file = paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."),
           list = c(ls(),ls(envir = globalenv())))
    # load(file=paste0("~/scShinyHubDebug/clusterPlot", "ns", ".RData", collapse = "."))
    
    g_id <- toupper(g_id)
    g_id <- gsub(" ", "", g_id, fixed = TRUE)
    g_id <- strsplit(g_id, ',')
    g_id<-g_id[[1]]
    
    notFound = featureData[!toupper(g_id) %in% featureData$Associated.Gene.Name, "Associated.Gene.Name"]
    if (length(featureData$Associated.Gene.Name) == length(notFound)) { # in case there is only one gene that is not available. 
      notFound = g_id
    }
    if(length(notFound)>0){
      if(DEBUG)cat(file=stderr(), paste("gene names not found: ",notFound, "\n"))
      if(!is.null(getDefaultReactiveDomain())){
        showNotification(paste("following genes were not found", notFound,collapse = " "), id="moduleNotFound", type = "warning", duration = 20)
      }
      
    }
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name %in%
                                           toupper(g_id)), ])[1]
    
    expression <- log2cpm[geneid, ]
    
    validate(need(is.na(sum(expression)) != TRUE, ''))
    
    projections <- cbind(projections, t(expression))
    names(projections)[names(projections) == geneid] <- 'exprs'
    
    if(DEBUG)cat(file=stderr(), paste("output$dge_plot1:---",ns(clId),"---\n"))
    subsetData <- subset(projections, dbCluster %in% clId)
    # subsetData$dbCluster = factor(subsetData$dbCluster)
    subsetData$shape = as.numeric(as.factor(subsetData$sample))
    if(DEBUGSAVE) 
      save(file = "~/scShinyHubDebug/clusterPlot.RData", list = c(ls(),ls(envir = globalenv())))
    # load(file="~/scShinyHubDebug/clusterPlot.RData")
    p1 <-
      ggplot(subsetData,
             aes_string(x = dimX, y = dimY)) +
      geom_point(aes_string(shape="shape", size = 2, color = 'exprs'), show.legend = TRUE) +
      scale_shape_identity() +
      geom_point(shape = 1,
                 size = 4,
                 aes(colour = as.numeric(dbCluster))) +
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
      ggtitle(paste(toupper(g_id), clId, sep = '-Cluster', collapse = " ")) +
      scale_fill_continuous()
    if( ! is.null(selectedCells)){
      shape = rep('a', nrow(subsetData))
      selRows = which(rownames(subsetData) %in% selectedCells)
      shape[selRows] = 'b'
      p1 = p1 + geom_point(mapping=aes(shape=shape, coloir='red', size=4, colour = dbCluster))
    }
    p1
  })
  
  observe({
    ns = session$ns
    input$changeGroups
    
    if(DEBUG)cat(file=stderr(), "cluster: changeGroups\n")
    
    isolate({
      brushedPs = (input$b1)
      projections = projections()
      inpClusters = (input$clusters)
      grpN = make.names(input$groupName)
      grpNs = groupNames$namesDF
      if(ncol(grpNs) == 0){
        initializeGroupNames()
        grpNs = groupNames$namesDF
      }
      if(is.null(projections) ){
        return(NULL)
      }     
      subsetData <- subset(projections, dbCluster %in% inpClusters)
      #if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
      cells.names <- brushedPoints(subsetData, brushedPs)
      
      if(DEBUGSAVE) 
        save(file = "~/scShinyHubDebug/changeGroups.RData", list = c(ls(),ls(envir = globalenv())))
      # load(file="~/scShinyHubDebug/changeGroups.RData")
      cat(file=stderr(), "cluster: changeGroups2\n")
      grpNs[, grpN] = FALSE
      grpNs[rownames(cells.names), grpN] = TRUE
      cat(file=stderr(), "cluster: changeGroups3\n")
      # projections[,grpN] = FALSE
      # projections[brushedPs,grpN] = TRUE
    })
  })
  
  output$additionalOptions <- renderUI({
    if(DEBUG)cat(file=stderr(), "cluster: additionalOptions\n")
    ns = session$ns
    moreOptions = (input$moreOptions)
    groupNs = groupNames$namesDF
    if(! moreOptions){return("")}
    
    if(DEBUGSAVE) 
      save(file = "~/scShinyHubDebug/additionalOptions.RData", list = c(ls(),ls(envir = globalenv())))
    # load(file="~/scShinyHubDebug/additionalOptions.RData")

        tagList(
      textInput(ns("groupName"), "name group"),
      selectInput(
        ns('groupNames'),
        label = 'group names',
        choice = colnames(groupNs),
        selected = selectedGroupName
      ),
      actionButton(ns('changeGroups'), 'change current selection'),
      checkboxInput(ns("showCells"), "show cell names", FALSE),
      verbatimTextOutput(ns('cellSelection'))
    )
  })
  
  output$cellSelection <- renderText({
    if(DEBUG)cat(file=stderr(), "cluster: cellSelection\n")
    ns = session$ns
    brushedPs = (input$b1)
    projections = projections()
    inpClusters = (input$clusters)
    myshowCells = (input$showCells)
    if(! myshowCells){return("")}
    if(is.null(projections) ){
      return("")
    }
    if(!is.null(getDefaultReactiveDomain())){
      showNotification("cluser cell Selection", id="clustercellSelection", duration = NULL)
    }
    if(DEBUGSAVE) 
      save(file = "~/scShinyHubDebug/clustercellSelection", list = c(ls(),ls(envir = globalenv())))
    # load(file=paste0("~/scShinyHubDebug/clustercellSelection", "ns", ".RData", collapse = "."))
    subsetData <- subset(projections, dbCluster %in% inpClusters)
    #if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
    cells.names <- brushedPoints(subsetData, brushedPs)
    retVal = paste(rownames(cells.names), collapse = ", ")
    if(DEBUG)cat(file=stderr(), "cluster: cellSelection: done\n")
    if(!is.null(getDefaultReactiveDomain())){
      removeNotification( id="clustercellSelection")
    }
    return(retVal)
  })
  
  
  
  
  return(reactive({
    returnValues
  }))
}


tableSelectionServer <- function(input, output, session, 
                                 dataTab){
  
  if(DEBUG)cat(file=stderr(), paste("tableSelectionServer", session$ns("test"), "\n"))
  ns = session$ns 
  
  output$cellSelection <- renderText({
    if(DEBUG)cat(file=stderr(), "cellSelection\n")
    ns = session$ns
    
    dataTables = dataTab()
    selectedRows = input$cellNameTable_rows_selected
    if(is.null(dataTables)){
      return(NULL)
    }
    if(!is.null(getDefaultReactiveDomain())){
      showNotification("cellSelection", id="cellSelection", duration = NULL)
    }
    if(DEBUGSAVE) 
      save(file = paste0("~/scShinyHubDebug/cellSelection", "ns", ".RData", collapse = "."),
           list = c(ls(),ls(envir = globalenv())))
    # load(file=paste0("~/scShinyHubDebug/cellSelection", "ns", ".RData", collapse = "."))
    
    if(length(selectedRows)>0){
      retVal = paste0(rownames(dataTables[selectedRows,]), collapse = ", ")
    }else{
      retVal = NULL
    }
    if(DEBUG)cat(file=stderr(), "cellSelection: done\n")
    if(!is.null(getDefaultReactiveDomain())){
      removeNotification( id="cellSelection")
    }
    return(retVal)
  })
  
  proxy = dataTableProxy('cellNameTable')
  
  observeEvent(input$selectAll, {
    if(DEBUG)cat(file=stderr(), "input$selectAll\n")
    ipSelect = input$selectAll
    prox = proxy
    allrows = input$cellNameTable_rows_all
    proxy %>% selectRows( NULL )
    if(DEBUGSAVE) save(file=paste0("~/scShinyHubDebug/inputselectAll.RData", collapse = "."), 
                       list= c(ls(),ls(envir = globalenv())))
    # load(file=paste0("~/scShinyHubDebug/inputselectAll.RData", collapse = "."))
    if(ipSelect){
      proxy %>% selectRows( allrows )
    }
  })
  
  output$cellNameTable <- renderDT({
    if(DEBUG)cat(file=stderr(), "output$cellNameTable\n")
    dataTables = dataTab()
    ns = session$ns
    
    if(is.null(dataTables)){
      return(NULL)
    }
    if(!is.null(getDefaultReactiveDomain())){
      showNotification("cellNameTable", id="cellNameTable", duration = NULL)
    }
    if(DEBUGSAVE) 
      save(file=paste0("~/scShinyHubDebug/cellNameTable", "ns", ".RData", collapse = "."), 
           list= c(ls(),ls(envir = globalenv())))
    # load(file=paste0("~/scShinyHubDebug/cellNameTable", "", ".RData", collapse = "."))
    
    if(DEBUG)cat(file=stderr(), "cellNameTable: done\n")
    if(!is.null(getDefaultReactiveDomain())){
      removeNotification( id="cellNameTable")
    }
    if (dim(dataTables)[1] > 1) {
      return(DT::datatable(dataTables, filter = 'top',
                           options = list(
                             orderClasses = TRUE,
                             autoWidth=TRUE
                           )))
    }else{
      return(NULL)
    }
  })
}




