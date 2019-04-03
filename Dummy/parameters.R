
#####
# normalization parameters
#####


# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "Dummy_Normalization"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  scEx_log = h4("no Parameters implemented")
)

Dummy_Normalization <- reactive({
  # just devide by number of cells and scale
  scalingFactor = 10000
  if (DEBUG) {
    cat(file = stderr(), "Dummy_Normalization\n")
  }
  scEx <- scEx()
  
  if (is.null(scEx)) {
    if (DEBUG) {
      cat(file = stderr(), "Dummy_Normalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/Dummy_Normalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/Dummy_Normalization.RData")
  
  # use_genes <- get_nonzero_genes(scEx)
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
  
  if (DEBUG) {
    cat(file = stderr(), "Dummy_Normalization:Done\n")
  }
  return(scEx_bcnorm)
})
