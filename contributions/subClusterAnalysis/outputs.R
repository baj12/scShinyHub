# TODO: make sure this is working
myZippedReportFiles <- c("DGE.csv")




# dge_plot1 ----
#' dge_plot1
#' left plot for selection of cells
output$dge_plot1 <- subCluster2Dplot()
<<<<<<< HEAD

=======
  
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
# SUBCLUSTER DGE PLOT2 -----
#' dge_plot2
#' right plot
output$dge_plot2 <- subCluster2Dplot()

# dgeTable ----
#' dgeTable
#' Table with differential expressed genes
<<<<<<< HEAD
# output$dgeTable <- DT::renderDataTable({
#   if (DEBUG) cat(file = stderr(), "output$dge\n")
#   scEx <- scEx()
#   top.genes <- dge()
# 
#     if (is.null(scEx)) {
#     return(NULL)
#   }
#   if (DEBUGSAVE) {
#     save(file = "~/scShinyHubDebug/output_dge.RData", list = c(ls(envir = globalenv(), ls())))
#   }
#   # load(file="~/scShinyHubDebug/output_dge.RData")
# 
#   featureData <- rowData(scEx)
#   
#   top.genes$symbol <-
#     featureData[rownames(top.genes), "symbol"]
#   if ("Description" %in% colnames(featureData)) {
#     top.genes$Description <- featureData[rownames(top.genes), "Description"]
#   }
#   if (dim(top.genes)[1] > 0) {
#     return(DT::datatable(top.genes,
#       options = list(
#         orderClasses = TRUE,
#         lengthMenu = c(10, 30, 50),
#         pageLength = 10
#       )
#     ))
#   } else {
#     return(NULL)
#   }
# })

dgeTableReac <- reactive({
  if (DEBUG) cat(file = stderr(), "output$dge\n")
  scEx <- scEx()
  top.genes <- dge()
  
  if (is.null(scEx)) {
=======
output$dgeTable <- DT::renderDataTable({
  if (DEBUG) cat(file = stderr(), "output$dge\n")
  scEx <- scEx()
  top.genes <- dge()

    if (is.null(scEx)) {
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/output_dge.RData", list = c(ls(envir = globalenv(), ls())))
  }
  # load(file="~/scShinyHubDebug/output_dge.RData")
<<<<<<< HEAD
  
=======

>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
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

# projectionTableMod ----
callModule(
  tableSelectionServer, 
  "dgeTable", 
  dgeTableReac)


<<<<<<< HEAD

=======
>>>>>>> 012c5c303a29ae54541730abda33b3d4e1088f20
# download differentially expressed genes
output$download_dge_table <- downloadHandler(
  filename = function() {
    paste("SubCluster", "DGE_table.csv", sep = "_")
  },
  content = function(file) {
    write.csv(selectedDge$dgeTable, file)
  }
)



