#' clusterServer
#' 
#' server side shiny module function for printing a 2D represenation of cells
#' it uses the global tsne.data object for plotting

#' @param gene_id name of the gene to be plotted (comma separated list, will be set to upper case)
#' @param tData tsne.data, dataframe with cluster numbers
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
                          legend.position = "none", 
                          DEBUG=TRUE){
  ns = session$ns 
  subsetData = NULL
  
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
    tsne.data = tData()
    ns = session$ns
    if(is.null(tsne.data)){
      HTML("Please load data first")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
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
    tsne.data = tData()
    returnValues$cluster = input$clusters
    
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      if(DEBUG)cat(file=stderr(), paste("output$clusterPlot:NULL\n"))
      return(NULL)
    }
    g_id=gene_id()
    g_id <- toupper(g_id)
    g_id <- gsub(" ", "", g_id, fixed = TRUE)
    g_id <- strsplit(g_id, ',')
    g_id<-g_id[[1]]
    
    notFound = featureData[!toupper(g_id) %in% featureData$Associated.Gene.Name, "Associated.Gene.Name"]
    if(length(notFound)>0){
      if(DEBUG)cat(file=stderr(), paste("gene names not found: ",notFound, "\n"))
    }
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name %in%
                                           toupper(g_id)), ])[1]
    
    expression <- log2cpm[geneid, ]
    
    validate(need(is.na(sum(expression)) != TRUE, ''))
    
    tsne.data <- cbind(tsne.data, t(expression))
    names(tsne.data)[names(tsne.data) == geneid] <- 'exprs'
    
    if(DEBUG)cat(file=stderr(), paste("output$dge_plot1:---",ns(input$clusters),"---\n"))
    subsetData <- subset(tsne.data, dbCluster %in% input$clusters)
    # subsetData$dbCluster = factor(subsetData$dbCluster)
    dimY = input$dimension_y
    dimX = input$dimension_x
    clId = input$clusters
    p1 <-
      ggplot(subsetData,
             aes_string(x = dimX, y = dimY)) +
      geom_point(aes_string(shape=18, size = 2, color = 'exprs'), show.legend = TRUE) +
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
    # save(file = "~/Desktop/clusterPlot.RData", list=c(ls(),"input"))
    # rm(list=ls())
    # load(file = "~/Desktop/clusterPlot.RData")
    if( ! is.null(selectedCells)){
      shape = rep('a', nrow(subsetData))
      selRows = which(rownames(subsetData) %in% selectedCells)
      shape[selRows] = 'b'
      p1 = p1 + geom_point(mapping=aes(shape=shape, size=4, colour = dbCluster))
    }
    p1
  })
  
  # TODO return selected cells
  
  return(reactive({
    returnValues
  }))
}

