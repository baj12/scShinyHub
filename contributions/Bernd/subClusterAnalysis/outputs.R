updateInputx1 <- reactive({
  tsneData <- projections()
  
  # Can use character(0) to remove all choices
  if (is.null(tsneData)){
    return(NULL)
  }
  
  # Can also set the label and select items
  updateSelectInput(session, "dimension_x1",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[1]
  )
  
  updateSelectInput(session, "dimension_y1",
                    choices = colnames(tsneData),
                    selected = colnames(tsneData)[2]
  )
})



# TODO module?  
output$dge_plot1 <- renderPlot({
  projections = projections()
  up1 <- updateInputx1()
  x1 = input$dimension_x1
  y1 = input$dimension_y1
  c1 = input$clusters1
  if( is.null(projections) | length(c1)==0 | length(x1)==0 | length(y1)==0 ){
    return(NULL)
  }
  
  if(DEBUG)cat(file=stderr(), paste("dge_plot1: x1: ", x1,"\n"))
  if(DEBUG)cat(file=stderr(), paste("dge_plot1: y1: ", y1,"\n"))
  if(DEBUG)cat(file=stderr(), paste("dge_plot1: c1: ", c1,"\n"))
  if(DEBUGSAVE) save(file="~/scShinyHubDebug/dge_plot1.RData", list=ls())
  # load(file="~/scShinyHubDebug/dge_plot1.RData")
  
  subsetData <- subset(projections, dbCluster %in% c1)
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
  p1
  
})
# SUBCLUSTER DGE PLOT2 ------------------------------------------------------------------

# TODO move to were it belongs  
# TODO module?  
output$dge_plot2 <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$dge_plot2\n")
  projections = projections()
  
  if( is.null(projections) ){
    return(NULL)
  }
  
  subsetData <- subset(projections, dbCluster %in% input$clusters1)
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
  if (dim(top.genes)[1] > 0) {
    return(DT::datatable(top.genes,
                  options = list(
                    orderClasses = TRUE,
                    lengthMenu = c(10, 30, 50),
                    pageLength = 10
                  ))
           )
  }else{
    return(NULL)
  }
  # })
})



# TODO module download?  
output$download_dge_table <- downloadHandler(
  filename = function() {
    paste("SubCluster", "DGE_table.csv", sep = '_')
  },
  content = function(file) {
    write.csv(selectedDge, file)
  }
)


# TODO as module
# sub cluster analysis ( used for 2 panels )
output$clusters1 <- renderUI({
  if(DEBUG)cat(file=stderr(), "output$clusters1\n")
  projections = projections()
  up1 <- updateInputx1()
  if(is.null(projections)){
    HTML("Please load data firts")
  }else{
    noOfClusters <- max(as.numeric(as.character(projections$dbCluster)))
    selectizeInput(
      "clusters1",
      label = "Cluster",
      choices = c(0:noOfClusters),
      selected = 0, 
      multiple = TRUE
    )
  }
})

