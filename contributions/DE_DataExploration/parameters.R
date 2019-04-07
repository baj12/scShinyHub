
# normalization parameters

# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "DE_logNormalization"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  scEx_log = h5("no Parameters implemented")
)

# DE_logNormalization ----
#' DE_logNormalization
#' reactive for normalizing data according to seurat
DE_logNormalization <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "DE_logNormalization")
  )
  if (DEBUG) cat(file = stderr(), "DE_logNormalization\n")
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("DE_logNormalization", id = "DE_logNormalization", duration = NULL)
  }
  
  scEx <- scEx()
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "DE_logNormalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/DE_logNormalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/DE_logNormalization.RData")
  
  # TODO ?? define scaling factor somewhere else???
  sfactor = max(max(assays(scEx)[["counts"]]),1000)
  retVal <- DE_logNormalizationfunc(scEx, scalingFactor = sfactor)
  
  printTimeEnd(start.time, "DE_logNormalization")
  exportTestValues(DE_logNormalization = {assays(retVal)[["logcounts"]]})  
  return(retVal)
})

#' DE_logNormalizationfunc
#' actual computation of the normalization as it is done in seurat
DE_logNormalizationfunc <- function(scEx, scalingFactor = 10000) {
  use_genes <- sort(unique(1 + slot(as(assays(scEx)[[1]], "dgTMatrix"), 
                                    "i")))
  
  bc_sums <- Matrix::colSums(assays(scEx)[[1]])
  median_sum <- median(bc_sums)
  A <- as(assays(scEx)[[1]], "dgCMatrix")
  A@x <- A@x / Matrix::colSums(A)[assays(scEx)[[1]]@j + 1L]
  scEx_bcnorm <- SingleCellExperiment(assay = list(logcounts = as(A,"dgTMatrix")),
                                      colData = colData(scEx),
                                      rowData = rowData(scEx))
  
  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}
