output$crHeat_plot1 <- renderPlot({
  if(DEBUG)cat(file=stderr(), "output$crHeat_plot1\n")
  gbm = gbm()
  tsne.data = tsne.data()
  prioritized_genes = prioritized_genes()
  if( is.null(gbm) | is.null(gbm) | is.null(prioritized_genes)){
    if(DEBUG)cat(file=stderr(), "output$crHeat_plot1:NULL\n")
    return(NULL)
  }
  
  example_K <- 10 
  example_Cols <- rev(brewer.pal(10,"Set3")) # customize plotting colors
  
  cells_to_plot <- order_cell_by_clusters(gbm, tsne.data$dbCluster)
  
  example_col = example_Cols[1:example_K]
  p = gbm_pheatmap(gbm=log_gene_bc_matrix(gbm), 
                   genes_to_plot=prioritized_genes, 
                   cells_to_plot=cells_to_plot,
                   n_genes=10, 
                   colour=example_col,
                   limits = c(-3, 3))
  
  if(DEBUG)cat(file=stderr(), "output$crHeat_plot1:done\n")
  p
})

output$crPrioGenes <- DT::renderDataTable({
  if(DEBUG)cat(file=stderr(), "output$crPrioGenes\n")
  prioritized_genes = prioritized_genes()
  if(is.null(prioritized_genes) )
    return(NULL)
  
  dt = prioritized_genes[[input$cluster5]]
  DT::datatable(dt)
})


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


# TODO as module ? 
# cell ranger output table
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

