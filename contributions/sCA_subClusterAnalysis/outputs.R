# TODO: make sure this is working
myZippedReportFiles <- c("DGE.csv")




# sCA_dge_plot1 ----
#' sCA_dge_plot1
#' left plot for selection of cells
output$sCA_dge_plot1 <- subCluster2Dplot()
# SUBCLUSTER DGE PLOT2 -----
#' sCA_dge_plot2
#' right plot
output$sCA_dge_plot2 <- subCluster2Dplot()

# sCA_dgeTable ----
#' sCA_dgeTable
#' Table with differential expressed genes

sCA_dgeTableReac <- reactive({
  if (DEBUG) cat(file = stderr(), "output$dge\n")
  scEx <- scEx()
  top.genes <- sCA_dge()
  
  if (is.null(scEx)) {
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/output_dge.RData", list = c(ls(envir = globalenv(), ls())))
  }
  # load(file="~/scShinyHubDebug/output_dge.RData")
  featureData <- rowData(scEx)
  
  top.genes$symbol <-
    featureData[rownames(top.genes), "symbol"]
  rownames(top.genes) <- top.genes$symbol
  if ("Description" %in% colnames(featureData)) {
    top.genes$Description <- featureData[rownames(top.genes), "Description"]
  }
  if (dim(top.genes)[1] > 0) {
    return(top.genes)
  } else {
    return(NULL)
  }
  
})

# dge table ----
callModule(
  tableSelectionServer, 
  "sCA_dgeTable", 
  sCA_dgeTableReac)





