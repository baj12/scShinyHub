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



# TODO module download?  
output$download_dge_table <- downloadHandler(
  filename = function() {
    paste("SubCluster", "DGE_table.csv", sep = '_')
  },
  content = function(file) {
    write.csv(selectedDge, file)
  }
)



