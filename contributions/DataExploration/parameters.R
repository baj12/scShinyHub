
# normalization parameters

# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "scEx_logNormalization"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  scEx_log = h5("no Parameters implemented")
)

scEx_logNormalization <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "Dummy_Normalization")
  )
  if (DEBUG) cat(file = stderr(), "scEx_logNormalization\n")
  
  scEx <- scEx()
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "scEx_logNormalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scEx_logNormalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/scEx_logNormalization.RData")
  
  # TODO ?? define scaling factor somewhere else???
  
  retVal <- scEx_logNormalizationfunc(scEx)
  
  if (DEBUG) {
    cat(file = stderr(), "scEx_logNormalization:Done\n")
  }
  return(retVal)
})

#' scEx_logNormalizationfunc
#' actual computation of the normalization as it is done in seurat
scEx_logNormalizationfunc <- function(scEx, scalingFactor = 10000) {
  use_genes <- sort(unique(1 + slot(as(assays(scEx)[["counts"]], "dgTMatrix"), 
                                    "i")))
  
  bc_sums <- Matrix::colSums(assays(scEx)[["counts"]])
  median_sum <- median(bc_sums)
  A <- as(assays(scEx)[["counts"]], "dgCMatrix")
  A@x <- A@x / Matrix::colSums(A)[assays(scEx)[["counts"]]@j + 1L]
  # new_matrix <- sweep(exprs(scEx), 2, median_sum/bc_sums, "*")
  # scEx_bcnorm <- (newGeneBCMatrix(A, fData(scEx), pData(scEx),
  #   template = scEx
  # ))
  scEx_bcnorm <- SingleCellExperiment(assay = list(logcounts = as(A,"dgTMatrix")),
                                      colData = colData(scEx),
                                      rowData = rowData(scEx))
  
  # gbm_bcnorm <- normalize_barcode_sums_to_median(gbm)
  # gbm_log <- log_gene_bc_matrix(gbm_bcnorm, base = 10)
  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  return(scEx_bcnorm)
}
