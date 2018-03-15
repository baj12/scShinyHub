

# LIBRARY -----------------------------------------------------------------

library(shiny)
library(shinyTree)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)
library(pheatmap)
library(threejs)
library(sm)
library(RColorBrewer)
library(mclust)
library(reshape)
library(cellrangerRkit)
library(SCORPIUS)
library(ggplot2)
library(knitr)
library(kableExtra)
library(shinyWidgets)
library(scater)

source("serverFunctions.R")
source("privatePlotFunctions.R")
DEBUG=TRUE

#needs to be an option
seed=1

shinyServer(function(input, output, session) {
  set.seed(seed)
  
  # FUNCTIONS ------------------------------------------------------------------
  
  options(shiny.maxRequestSize = 2000 * 1024 ^ 2)
  
  # TODO check if file exists
  load(file = "geneLists.RData")
  
  cat(file=stderr(), "ShinyServer running\n")
  
  # load global reactives
  source("reactives.R", local = TRUE)
  
  # load contribution reactives
  # parse all reactives.R files under contributions to include in application
  uiFiles = dir(path = "contributions", pattern = "reactives.R", full.names = TRUE, recursive = TRUE)
  for(fp in uiFiles){
    source(fp, local = TRUE)
  }
  
  
  # load contribution outputs
  # parse all outputs.R files under contributions to include in application
  uiFiles = dir(path = "contributions", pattern = "outputs.R", full.names = TRUE, recursive = TRUE)
  for(fp in uiFiles){
    source(fp, local = TRUE)
  }
  
  
  
  
  
  
  
  positiveCells <- reactiveValues(positiveCells = NULL,
                                  positiveCellsAll = NULL)
  selectedDge <- reactiveValues()
  
  forceCalc <-observe({
    input$goCalc
    isolate({
      if(DEBUG)cat(file=stderr(), "forceCalc\n")
      withProgress(message = 'Performing heavy calculations', value = 0, {
        n=5
        incProgress(1/n, detail = paste("Creating pca"))
        if(DEBUG)cat(file=stderr(), "forceCalc PCA\n")
        pca = pca()
        incProgress(1/n, detail = paste("Creating kmClustering"))
        if(DEBUG)cat(file=stderr(), "forceCalc kmClustering\n")
        kmClustering = kmClustering()
        incProgress(1/n, detail = paste("Creating tsne"))
        if(DEBUG)cat(file=stderr(), "forceCalc tsne\n")
        tsne = tsne()
        incProgress(1/n, detail = paste("Creating prioritized_genes"))
        if(DEBUG)cat(file=stderr(), "forceCalc prioritized_genes\n")
        prioritized_genes()
        incProgress(1/n, detail = paste("Creating scaterReads"))
        if(DEBUG)cat(file=stderr(), "forceCalc scaterReads\n")
        scaterReads = scaterReads()
      })
    })
  })
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # RENDER UI  ------------------------------------------------------------------
  # TODO as module
  output$clusters <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data first")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectizeInput(
        "cluster",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0, 
        multiple = TRUE
      )
    }
  })
  
  # TODO as module
  output$clusters1 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters1\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectizeInput(
        "clusters1",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0, 
        multiple = TRUE
      )
    }
  })
  
  # TODO as module
  output$clusters2 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters2\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectizeInput(
        "clusters2",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0, 
        multiple = TRUE
      )
    }
  })
  
  # TODO as module
  output$clusters3 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters3\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectizeInput(
        "clusters3",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0, 
        multiple = TRUE
      )
    }
  })
  
  # TODO as module
  output$clusters4 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters4\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data firts")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectInput(
        "clusters4",
        label = "Cluster",
        choices = c(c('All'),c(0:noOfClusters)),
        selected = 0
      )
    }
  })  
  
  # TODO as module
  output$clusters5 <- renderUI({
    if(DEBUG)cat(file=stderr(), "output$clusters\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      HTML("Please load data first")
    }else{
      noOfClusters <- max(tsne.data$dbCluster)
      selectInput(
        "cluster5",
        label = "Cluster",
        choices = c(0:noOfClusters),
        selected = 0
      )
    }
  })
  
  # MAIN 3D PLOT ------------------------------------------------------------------
  # TODO move to generalQC
    output$tsne_main <- renderPlotly({
    if(DEBUG)cat(file=stderr(), "output$tsne_main\n")
    tsne.data = tsne.data()
    if(is.null(tsne.data)){
      if(DEBUG)cat(file=stderr(), "output$tsne_main:NULL\n")
      return(NULL)
    }
    tsne.data <- as.data.frame(tsne.data)
    #cat(stderr(),colnames(tsne.data)[1:5])
    tsne.data$dbCluster <- as.factor(tsne.data$dbCluster)
    
    p <-
      plot_ly(
        tsne.data,
        x = ~ V1,
        y = ~ V2,
        z = ~ V3,
        type = "scatter3d",
        color =  ~ dbCluster,
        hoverinfo = "text",
        text = paste('Cluster:', tsne.data$dbCluster),
        mode = 'markers',
        marker =
          list(
            line = list(width = 0),
            size = rep(10, nrow(tsne.data)),
            sizeref = 3
          )
      )
    if(DEBUG)cat(file=stderr(), "output$tsne_main: done\n")
    return(layout(p))
    
    
  })
  
  # SUMMARY STATS ----------------------------------------------------------------
  

  output$summaryStatsSideBar<-renderUI({
    if(DEBUG)cat(file=stderr(), "output$summaryStatsSideBar\n")
    log2cpm = log2cpm()
    if(is.null(log2cpm) ){
      if(DEBUG)cat(file=stderr(), "output$summaryStatsSideBar:NULL\n")
      return(NULL)
    }
    
    line1<-paste('No. of cells:', dim(log2cpm)[2],sep='\t')
    line2<-paste('No. of genes:',  dim(log2cpm)[1],sep='\t')
    line3<-paste('Median UMIs per cell:', medianUMI(),sep='\t')
    line4<-paste('Median Genes with min 1 UMI:', medianENSG(),sep='\t')
    HTML(
      paste0("Summary statistics of this dataset:", '<br/>','<br/>',
             line1, '<br/>',
             line2, '<br/>',
             line3, '<br/>',
             line4
      )
    )
  })
  
  # General QC ----------------------------------------------------------------
  # TODO move to generalQC
  output$plotUmiHist <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output_plotUmiHist\n")
    gbm = gbm()
    if(is.null(gbm))
      return(NULL)
    hist(colSums(as.matrix(exprs(gbm))), breaks = 50, main="histogram of number of UMI per cell")
  })
  
  # TODO move to generalQC
  output$variancePCA <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$variancePCA\n")
    h2("hello")
    pca = pca()
    if(is.null(pca))
      return(NULL)
    barplot(pca$var_pcs, main="Variance captured by first PCs")
  })
  
  
  # TODO move to generalQC
  output$plotHistogramsAll <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$plotHistogramsAll\n")
    h2("hello")
    pca = pca()
    if(is.null(pca))
      return(NULL)
    barplot(pca$var_pcs, main="Variance captured by first PCs")
    
  })
  # Select Genes ----------------------------------------------------------------
  
  output$geneListSelection <- renderTree({ 
    geneLists
  })
  

  # ONOFF TAB RENDER TABLE ALL CELLS ------------------------------------------------------------------
# TODO module for DT
# TODO move to were it belongs  
    output$selectedGenesTable <- DT::renderDataTable({
    if(DEBUG)cat(file=stderr(), "output$selectedGenesTable\n")
    dataTables = inputData()
    useGenes = useGenes()
    useCells = useCells()
    if(is.null(dataTables) | is.null(useGenes) | is.null(useCells))
      return(NULL)
    
    gbm = as.matrix(exprs(dataTables$gbm))
    fd = dataTables$featuredata
    dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
    dt$rowSums = rowSums(gbm[useGenes,useCells])
    dt$rowSamples = rowSums(gbm[useGenes,useCells]>0)
    DT::datatable(dt)
  })
  
  # TODO module for DT
  # TODO move to were it belongs  
  output$removedGenesTable <- DT::renderDataTable({
    if(DEBUG)cat(file=stderr(), "output$removedGenesTable\n")
    dataTables = inputData()
    useGenes = !useGenes()
    useCells = useCells()
    if(is.null(dataTables) | is.null(useGenes)  | is.null(useCells) )
      return(NULL)
    
    gbm = as.matrix(exprs(dataTables$gbm))
    fd = dataTables$featuredata
    dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
    dt$rowSums = rowSums(gbm[useGenes,useCells])
    dt$rowSamples = rowSums(gbm[useGenes,useCells]>0)
    DT::datatable(dt)
  })
  
  
  ### Cell selection
  
  
  # EXPLORE TAB 3D PLOT ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  output$tsne_plt <- renderPlotly({
    if(DEBUG)cat(file=stderr(), "output$tsne_plt\n")
    # if (v$doPlot == FALSE)
    #   return()
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                           toupper(input$gene_id)), ])[1]
    
    expression <- log2cpm[geneid, ]
    cat(file = stderr(), rownames(expression))
    
    validate(need(
      is.na(sum(expression)) != TRUE,
      'Gene symbol incorrect or gene not expressed'
    ))
    
    tsne.data <- cbind(tsne.data, t(expression))
    names(tsne.data)[names(tsne.data) == geneid] <- 'values'
    
    p <-
      plot_ly(
        tsne.data,
        x = ~ V1,
        y = ~ V2,
        z = ~ V3,
        type = "scatter3d",
        hoverinfo = "text",
        text = paste('Cluster:', tsne.data$dbCluster),
        mode = 'markers',
        marker = list(
          size = 2,
          line = list(width = 0),
          color =  ~ values,
          colors = 'Greens'
        )
      )
    layout(p, title = toupper(input$gene_id))
    # })
  })
  
  # EXPLORE TAB CLUSTER PLOT ------------------------------------------------------------------
  # TODO move to were it belongs  
  output$clusterPlot <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$clusterPlot\n")
    # if (v$doPlot == FALSE)
    #   return()
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                           toupper(input$gene_id)), ])[1]
    
    expression <- log2cpm[geneid, ]
    
    validate(need(is.na(sum(expression)) != TRUE, ''))
    
    tsne.data <- cbind(tsne.data, t(expression))
    names(tsne.data)[names(tsne.data) == geneid] <- 'values'
    
    if(DEBUG)cat(file=stderr(), paste("output$dge_plot1:---",input$cluster,"---\n"))
    subsetData <- subset(tsne.data, dbCluster %in% input$cluster)
    p1 <-
      ggplot(subsetData,
             aes_string(x = input$dimension_x, y = input$dimension_y)) +
      geom_point(aes_string(size = 2, color = 'values')) +
      geom_point(shape = 1,
                 size = 4,
                 aes(colour = dbCluster)) +
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
      ggtitle(paste(toupper(input$gene_id), input$cluster, sep = '-Cluster')) 
    # +
    #   scale_colour_gradientn(colours = terrain.colors(10))
    p1
  })
  
  # EXPLORE TAB VIOLIN PLOT ------------------------------------------------------------------
  # TODO move to were it belongs  
  # TODO module for violin plot  
  output$gene_vio_plot <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$gene_vio_plot\n")
    # if (v$doPlot == FALSE)
    #   return()
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      if(DEBUG)cat(file=stderr(), "output$gene_vio_plot:NULL\n")
      return(NULL)
    }
    
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                           toupper(input$gene_id)), ])[1]
    
    expression <- log2cpm[geneid, ]
    
    validate(need(is.na(sum(expression)) != TRUE, ''))
    
    tsne.data <- cbind(tsne.data, t(expression))
    names(tsne.data)[names(tsne.data) == geneid] <- 'values'
 
    p1 <-
      ggplot(tsne.data, aes(factor(dbCluster), values, fill = factor(dbCluster))) +
      geom_violin(scale = "width") +
      stat_summary(
        fun.y = median,
        geom = "point",
        size = 5,
        color = 'black'
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
      xlab('Cluster') +
      ylab('Expression') +
      ggtitle(toupper(input$gene_id))
    if(DEBUG)cat(file=stderr(), "output$gene_vio_plot:done\n")
    return(p1)
    # })
  })
  
  # EXPLORE TABL DOWNLOAD SELECTED WITH BRUSH ------------------------------------------------------------------
  # TODO move to were it belongs  
  # TODO module for download?  
  output$downloadExpression <- downloadHandler(
    filename = function() {
      paste(input$cluster, "Selected_Expression_table.csv", sep = '_')
    },
    content = function(file) {
      featureData = featureDataReact()
      log2cpm = log2cpm()
      tsne.data = tsne.data()
      if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
        return(NULL)
      }
      geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                             toupper(input$gene_id)), ])[1]
      
      expression <- log2cpm[geneid, ]
      #cat(stderr(),colnames(expression)[1:5])
      tsne.data <- cbind(tsne.data, t(expression))
      #if(DEBUG)cat(file=stderr(),grep('^T_',rownames(tsne.data)))
      
      names(tsne.data)[names(tsne.data) == geneid] <- 'values'
      
      #if(DEBUG)cat(file=stderr(),grep('^T_',rownames(tsne.data)))
      
      subsetData <- subset(tsne.data, dbCluster == input$cluster)
      #if(DEBUG)cat(file=stderr(),rownames(subsetData)[1:5])
      cells.names <- brushedPoints(subsetData, input$b1, allRows = T)
      #if(DEBUG)cat(file=stderr(),colnames(cells.names))
      cells <-
        rownames(subsetData[which(cells.names$selected_ == TRUE), ])
      #if(DEBUG)cat(file=stderr(),cells[1:5])
      
      if (length(cells) == 1) {
        subsetExpression <- log2cpm[, cells]
        subsetExpression <-
          as.data.frame(subsetExpression, row.names = rownames(log2cpm))
        colnames(subsetExpression) <- cells
        subsetExpression$Associated.Gene.Name <-
          featureData[rownames(subsetExpression), 'Associated.Gene.Name']
        write.csv(subsetExpression, file)
      }
      else{
        subsetExpression <- log2cpm[, cells]
        #cat(stderr(),colnames(subsetExpression)[1:5])
        
        subsetExpression$Associated.Gene.Name <-
          featureData[rownames(subsetExpression), 'Associated.Gene.Name']
        #cat(stderr(),colnames(subsetExpression))
        write.csv(subsetExpression, file)
      }
    }
  )
  # EXPLORE TAB PANEL PLOT------------------------------------------------------------------
  
  # TODO move to were it belongs  
  output$panelPlot <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$panelPlot\n")
    
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    
    genesin <- input$panelplotids
    genesin <- toupper(genesin)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ',')
    genesin<-genesin[[1]]
    
    if(DEBUG)cat(file=stderr(),length(genesin))
    par(mfrow=c(ceiling(length(genesin)/4),4), mai = c(0, 0., 0., 0.))
    rbPal <- colorRampPalette(c('#f0f0f0','red'))
    if(DEBUG)cat(file=stderr(),input$clusters4)
    
    if (input$clusters4 == 'All') 
    {
      for (i in 1:length(genesin)){
        Col <- rbPal(10)[
          as.numeric(
            cut(
              as.numeric(
                log2cpm[
                  rownames(featureData[which(featureData$Associated.Gene.Name==genesin[i]),])
                  ,]
              ),breaks = 10))]
        plot(tsne.data[,input$dimension_x4],tsne.data[,input$dimension_y4],col=Col,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
        title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
        if(DEBUG)cat(file=stderr(),genesin[i])
      }
    }
    else{
      for (i in 1:length(genesin)){
        
        subsetTSNE <- subset(tsne.data, dbCluster == input$clusters4)
        
        Col <- rbPal(10)[
          as.numeric(
            cut(
              as.numeric(
                log2cpm[
                  rownames(featureData[which(featureData$Associated.Gene.Name==genesin[i]),])
                  ,]
              ),breaks = 10))]
        
        names(Col)<-rownames(tsne.data)
        plotCol<-Col[rownames(subsetTSNE)]
        plot(subsetTSNE[,input$dimension_x4],subsetTSNE[,input$dimension_y4],col=plotCol,pch=16,axes = FALSE,frame.plot = TRUE, ann=FALSE)
        title(genesin[i],line=-1.2,adj = 0.05,cex.main=2)
        if(DEBUG)cat(file=stderr(),input$clusters4)
      }
    }
  })
  
  # Explore tab scater QC plot
  
  # TODO move to were it belongs  
  output$scaterQC <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$scaterQC\n")
    scaterReads = scaterReads()
    if(is.null(scaterReads)){
      return(NULL)
    }
    scater::plotQC(scaterReads, type = "highest-expression", col_by_variable="fixed")
    #if(DEBUG)cat(file=stderr(), "DiffExpTest\n")
  }
  )
  
  # CO-EXPRESSION HEATMAP ALL CLUSTERS ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO mnodule for heatmap?  
  output$heatmap <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$heatmap\n")
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    
    genesin <- input$heatmap_geneids
    genesin <- toupper(genesin)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ',')
    
    map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
    cat(file = stderr(), length(map))
    
    expression <- log2cpm[map, ]
    
    validate(need(
      is.na(sum(expression)) != TRUE,
      'Gene symbol incorrect or genes not expressed'
    ))
    
    tsne.data <- tsne.data[order(tsne.data$dbCluster), ]
    
    expression <- expression[, rownames(tsne.data)]
    expression <- expression[complete.cases(expression), ]
    
    annotation <- data.frame(factor(tsne.data$dbCluster))
    rownames(annotation) <- colnames(expression)
    colnames(annotation) <- c('Cluster')
    
    h <-
      pheatmap(
        as.matrix(expression),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        scale = 'row',
        fontsize_row = 10,
        labels_col = colnames(expression),
        labels_row = featureData[rownames(expression), 'Associated.Gene.Name'],
        show_rownames = TRUE,
        annotation_col = annotation,
        show_colnames = FALSE,
        annotation_legend = TRUE,
        breaks = seq(-6, 6, by = .12),
        colorRampPalette(rev(brewer.pal(
          n = 6, name =
            "RdBu"
        )))(100)
        
      )
    h
    
    # })
  })
  
  # CO EXPRESSION TAB CLUSTER PLOT ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO module for cluster plot?  
  output$clusterPlot2 <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$clusterPlot2\n")
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    
    # isolate({
    geneid <- rownames(featureData[which(featureData$Associated.Gene.Name ==
                                           toupper(input$gene_id_sch)), ])[1]
    
    expression <- log2cpm[geneid, ]
    
    validate(need(is.na(sum(expression)) != TRUE, ''))
    
    tsne.data <- cbind(tsne.data, t(expression))
    names(tsne.data)[names(tsne.data) == geneid] <- 'values'
    
    if(DEBUG)cat(file=stderr(), paste("output$dge_plot1:---",input$clusters2,"---\n"))
    subsetData <- subset(tsne.data, dbCluster %in% input$clusters2)
    p1 <-
      ggplot(subsetData,
             aes_string(x = input$dimension_x2, y = input$dimension_y2)) +
      geom_point(aes_string(size = 2, color = 'values')) +
      geom_point(shape = 1,
                 size = 4,
                 aes(colour = dbCluster)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          size = 12,
          vjust = 0.5
        ),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        legend.position = "none"
      ) +
      ggtitle(paste(toupper(input$gene_id_sch), input$clusters2, sep =
                      '-Cluster')) +
      scale_colour_gradient2(low = 'grey50', high = "red")
    p1
    # })
  })
  # CO EXPRESSION TAB SELECTED HEATMAP ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO module for heatmap?  
  output$selectedHeatmap <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$selectedHeatmap\n")
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    
    # isolate({
    genesin <- input$heatmap_geneids2
    genesin <- toupper(genesin)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ',')
    
    subsetData <-
      subset(tsne.data, dbCluster %in% input$clusters2)
    cells.1 <- rownames(brushedPoints(subsetData, input$scb1))
    
    
    map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
    #if(DEBUG)cat(file=stderr(),map[1])
    
    expression <- log2cpm[map, cells.1]
    if(DEBUG)cat(file = stderr(), rownames(expression))
    
    expression <- expression[complete.cases(expression), ]
    if(DEBUG)cat(file = stderr(), rownames(expression))
    mColor <- max(expression)
    
    validate(need(
      is.na(sum(expression)) != TRUE,
      'Gene symbol incorrect or genes not expressed'
    ))
    
    #annotation<-data.frame(factor(tsne.data$dbCluster))
    #rownames(annotation)<-colnames(expression)
    #colnames(annotation)<-c('Cluster')
    
    h <-
      pheatmap(
        as.matrix(expression),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = 'row',
        fontsize_row = 10,
        labels_col = colnames(expression),
        labels_row = featureData[rownames(expression), 'Associated.Gene.Name'],
        show_rownames = TRUE,
        show_colnames = FALSE,
        breaks = seq(-6, 6, by = .12),
        colorRampPalette(rev(brewer.pal(
          n = 6, name =
            "RdBu"
        )))(100)
        
      )
    h
    # })
  })
  
  # CO EXPRESSION TAB ON/OFF PLOT ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO module?  
  output$plotCoExpression <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$plotCoExpression\n")
    # if (vvvvvvv$doPlot == FALSE)
    #   return()
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(tsne.data) | is.null(log2cpm) | is.null(input$clusters3)){
      return(NULL)
    }
    
    # isolate({
    genesin <- input$mclustids
    genesin <- toupper(genesin)
    genesin <- gsub(" ", "", genesin, fixed = TRUE)
    genesin <- strsplit(genesin, ',')
    
    subsetData <-
      subset(tsne.data, dbCluster %in% input$clusters3)
    cells.1 <- rownames(subsetData)
    
    
    map <- rownames(featureData[which(featureData$Associated.Gene.Name %in% genesin[[1]]), ])
    #if(DEBUG)cat(file=stderr(),map[1])
    
    expression <- log2cpm[map, ]
    #if(DEBUG)cat(file=stderr(),rownames(expression))
    
    #expression<-expression[complete.cases(expression),]
    #if(DEBUG)cat(file=stderr(),rownames(expression))
    
    validate(need(
      is.na(sum(expression)) != TRUE,
      'Gene symbol incorrect or genes not expressed'
    ))
    
    bin <- expression
    bin[] <- 0
    
    for (i in 1:nrow(expression))
    {
      x <- Mclust(expression[i, ], G = 2)
      bin[i, ] <- x$classification
    }
    bin <- bin - 1
    allexprs <- apply(bin, 2, sum)
    plotexprs <- allexprs
    plotexprs[] <- 0
    plotexprs[allexprs >= length(rownames(bin))] <- 1
    positiveCells$positiveCells <- allexprs >= length(rownames(bin))
    positiveCells$positiveCellsAll <- plotexprs
    #save(subsetData,bin,allexprs,file='~/Desktop/test.Rds')
    #if(DEBUG)cat(file=stderr(),names(allexprs))
    
    mergeExprs <- plotexprs[rownames(subsetData)]
    #if(DEBUG)cat(file=stderr(),length(mergeExprs))
    
    subsetData$CoExpression <- mergeExprs
    #if(DEBUG)cat(file=stderr(),colnames(subsetData))
    
    p1 <-
      ggplot(subsetData,
             aes_string(x = input$dimension_x3, y = input$dimension_y3)) +
      geom_point(aes_string(size = 2, color = 'CoExpression')) +
      geom_point(shape = 1,
                 size = 4,
                 aes(colour = dbCluster)) +
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
      #ggtitle(paste(toupper(input$gene_id),input$cluster,sep='-Cluster'))+
      scale_colour_gradient2(low = 'grey50', high = "red")
    p1
    # })
  })
  
  
  # ONOFF TAB DOWNLOAD POSITIVECELLS ------------------------------------------------------------------
  # TODO move to were it belongs  
  # TODO module download  
  output$downloadExpressionOnOff <- downloadHandler(
    filename = function() {
      paste(input$clusters3, "PositiveCells.csv", sep = '_')
    },
    content = function(file) {
      featureData = featureDataReact()
      log2cpm = log2cpm()
      
      if(is.null(featureData) | is.null(log2cpm) | is.null(positiveCells$positiveCells)){
        return(NULL)
      }
      
      cells <- positiveCells$positiveCells
      #if(DEBUG)cat(file=stderr(),cells[1:5])
      
      if (length(cells) == 1) {
        subsetExpression <- log2cpm[, cells]
        subsetExpression <-
          as.data.frame(subsetExpression, row.names = rownames(log2cpm))
        colnames(subsetExpression) <- cells
        subsetExpression$Associated.Gene.Name <-
          featureData[rownames(subsetExpression), 'Associated.Gene.Name']
        write.csv(subsetExpression, file)
      }
      else{
        subsetExpression <- log2cpm[, cells]
        #cat(stderr(),colnames(subsetExpression)[1:5])
        subsetExpression$Associated.Gene.Name <-
          featureData[rownames(subsetExpression), 'Associated.Gene.Name']
        #cat(stderr(),colnames(subsetExpression))
        write.csv(subsetExpression, file)
      }
    }
  )
  
  # ONOFF TAB RENDER TABLE ALL CELLS ------------------------------------------------------------------
  # TODO move to were it belongs  
  # TODO do we need it?  
  output$onOffTable <- DT::renderDataTable({
    if(DEBUG)cat(file=stderr(), "output$onOffTable\n")
    tsne.data = tsne.data()
    
    if( is.null(tsne.data | is.null(positiveCells$positiveCellsAll)) ){
      return(NULL)
    }
    
    merge <- tsne.data
    if(DEBUG)cat(file=stderr(), paste("positiveCells$positiveCellsAll:---",positiveCells$positiveCellsAll,"---\n"))
    
    merge$CoExpression <- positiveCells$positiveCellsAll
    df <-
      as.data.frame(table(merge[, c('dbCluster', 'CoExpression')]))
    dfOut <- cast(df, dbCluster ~ CoExpression)
    colnames(dfOut) <- c("Cluster", 'OFF', 'ON')
    rownames(dfOut) <- dfOut$Cluster
    dfOut['Sum', ] <- c('', sum(dfOut$OFF), sum(dfOut$ON))
    DT::datatable(dfOut)
    
  })
  
  
  # SUBCLUSTER DGE PLOT1 ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO module?  
  output$dge_plot1 <- renderPlot({
    tsne.data = tsne.data()
    if( is.null(tsne.data) ){
      return(NULL)
    }
    
    x1 = input$dimension_x1
    y1 = input$dimension_y1
    c1 = input$clusters1
    subsetData <- subset(tsne.data, dbCluster %in% c1)
    p1 <-
      ggplot(subsetData,
             aes_string(x = x1, y = y1),
             colour = "dbCluster") +
      geom_point(aes(colour = dbCluster)) +
      geom_point(shape = 1,
                 size = 4,
                 aes(colour = dbCluster)) +
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
    #save(file="~/Desktop/test.RData", list=ls())
    
    p1
    
  })
  # SUBCLUSTER DGE PLOT2 ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO module?  
  output$dge_plot2 <- renderPlot({
    if(DEBUG)cat(file=stderr(), "output$dge_plot2\n")
    tsne.data = tsne.data()
    if( is.null(tsne.data) ){
      return(NULL)
    }
    
    subsetData <- subset(tsne.data, dbCluster %in% input$clusters1)
    p1 <-
      ggplot(subsetData,
             aes_string(x = input$dimension_x1, y = input$dimension_y1),
             color = "dbCluster") +
      geom_point(aes(colour = dbCluster)) +
      geom_point(shape = 1,
                 size = 4,
                 aes(colour = dbCluster)) +
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
      ggtitle(input$clusters1)
    p1
    # })
  })
  # SUBCLUSTER DGE ANALYSIS ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  dge <- reactive({
    if(DEBUG)cat(file=stderr(), "dge\n")
    featureData = featureDataReact()
    log2cpm = log2cpm()
    tsne.data = tsne.data()
    if(is.null(featureData) | is.null(log2cpm) | is.null(tsne.data)){
      return(NULL)
    }
    
    subsetData <- subset(tsne.data, dbCluster %in% input$clusters1)
    cells.1 <- rownames(brushedPoints(subsetData, input$db1))
    cells.2 <- rownames(brushedPoints(subsetData, input$db2))
    
    subsetExpression <- log2cpm[, union(cells.1, cells.2)]
    
    genes.use <- rownames(subsetExpression)
    data.1 = apply(subsetExpression[genes.use, cells.1], 1, expMean)
    data.2 = apply(subsetExpression[genes.use, cells.2], 1, expMean)
    total.diff = (data.1 - data.2)
    
    genes.diff = names(which(abs(total.diff) > .2))
    genes.use = ainb(genes.use, genes.diff)
    
    toReturn <-
      DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
    toReturn[, "avg_diff"] = total.diff[rownames(toReturn)]
    toReturn$Associated.Gene.Name <-
      featureData[rownames(toReturn), 'Associated.Gene.Name']
    selectedDge <- toReturn
    cat(stderr(), rownames(toReturn)[1:5])
    return(toReturn)
    
    # })
  })
  # SUBCLUSTER DGE OUTPUT TABLE ------------------------------------------------------------------
  
  # TODO move to were it belongs  
  # TODO module DT?  
  output$dge <- DT::renderDataTable({
    if(DEBUG)cat(file=stderr(), "output$dge\n")
    featureData = featureDataReact()
    if(is.null(featureData)){
      return(NULL)
    }
    
    # isolate({
    top.genes <- dge()
    top.genes$Associated.Gene.Name <-
      featureData[rownames(top.genes), 'Associated.Gene.Name']
    if (dim(top.genes)[1] > 1) {
      DT::datatable(top.genes,
                    options = list(
                      orderClasses = TRUE,
                      lengthMenu = c(10, 30, 50),
                      pageLength = 10
                    ))
    }
    # })
  })
  
  
  # TODO move to were it belongs  
  output$crSelectedGenes <- renderText({
    if(DEBUG)cat(file=stderr(), "crSelectedGenes\n")
    featureData = featureDataReact()
    if(is.null(featureData)){
      return(NULL)
    }
    top.genes <- dge()
    top.genes$Associated.Gene.Name <-
      featureData[rownames(top.genes), 'Associated.Gene.Name']
    
    paste0(top.genes$Associated.Gene.Name[input$dge_rows_selected],",")
  })
  
  # TODO move to were it belongs  
  output$gsSelectedGenes <- renderText({
    if(DEBUG)cat(file=stderr(), "gsSelectedGenes\n")
    dataTables = inputData()
    useGenes = useGenes()
    useCells = useCells()
    if(is.null(dataTables) | is.null(useGenes) | is.null(useCells))
      return(NULL)
    
    gbm = as.matrix(exprs(dataTables$gbm))
    fd = dataTables$featuredata
    dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
    
    paste0(dt$Associated.Gene.Name[input$selectedGenesTable_rows_selected],",")
  })
  
  # TODO move to were it belongs  
  output$gsrmGenes <- renderText({
    if(DEBUG)cat(file=stderr(), "gsrmGenes\n")
    dataTables = inputData()
    useGenes = !useGenes()
    useCells = useCells()
    if(is.null(dataTables) | is.null(useGenes) | is.null(useCells))
      return(NULL)
    
    gbm = as.matrix(exprs(dataTables$gbm))
    fd = dataTables$featuredata
    dt = fd[useGenes,c("Associated.Gene.Name", "Gene.Biotype", "Description")]
    
    paste0(dt$Associated.Gene.Name[input$removedGenesTable_rows_selected],",")
  })
  
  
  # SUBCLUSTER DGE DOWNLOADS ------------------------------------------------------------------
  # TODO move to were it belongs  
  # TODO module download?  
  output$download_dge_table <- downloadHandler(
    filename = function() {
      paste("SubCluster", "DGE_table.csv", sep = '_')
    },
    content = function(file) {
      write.csv(selectedDge, file)
    }
  )
  
  
  
  # TODO modularize:
  # create template for input parameters
  # source other reports
  output$report <- downloadHandler(
    filename = "report.html",
    
    content = function(file) {
      
      ip = inputData()
      if( is.null(ip) ){
        if(DEBUG)cat(file=stderr(), "output$report:NULL\n")
        return(NULL)
      }
      
      
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      tempServerFunctions <- file.path(tempdir(), "serverFunctions.R")
      file.copy("serverFunctions.R", tempServerFunctions, overwrite = TRUE)
      tempprivatePlotFunctions <- file.path(tempdir(), "privatePlotFunctions.R")
      file.copy("privatePlotFunctions.R", tempprivatePlotFunctions, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        tempServerFunctions = tempServerFunctions,
        tempprivatePlotFunctions = tempprivatePlotFunctions,
        b1 = input$b1,
        cluster = input$cluster,
        cluster5 = input$cluster5,
        clusters = input$clusters,
        clusters1 = input$clusters1,
        clusters2 = input$clusters2,
        clusters3 = input$clusters3,
        clusters4 = input$clusters4,
        db1 = input$db1,
        db2 = input$db2,
        dimension_x = input$dimension_x,
        dimension_x1 = input$dimension_x1,
        dimension_x2 = input$dimension_x2,
        dimension_x3 = input$dimension_x3,
        dimension_x4 = input$dimension_x4,
        dimension_y = input$dimension_y,
        dimension_y1 = input$dimension_y1,
        dimension_y2 = input$dimension_y2,
        dimension_y3 = input$dimension_y3,
        dimension_y4 = input$dimension_y4,
        file1 = input$file1,
        gene_id = input$gene_id,
        gene_id_sch = input$gene_id_sch,
        geneListSelection = input$geneListSelection,
        heatmap_geneids = input$heatmap_geneids,
        heatmap_geneids2 = input$heatmap_geneids2,
        maxGenes = input$maxGenes,
        mclustids = input$mclustids,
        minExpGenes = input$minExpGenes,
        minGenes = input$minGenes,
        minGenesGS = input$minGenesGS,
        panelplotids = input$panelplotids,
        # positiveCells = positiveCells$positiveCells,
        # positiveCellsAll = positiveCells$positiveCellsAll,
        scb1 = input$scb1,
        selectIds = input$selectIds
      )
      if(DEBUG)cat(file=stderr(), "output$report:gbm:\n")
      if(DEBUG)cat(file=stderr(), str(gbm))
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
})# END SERVER





