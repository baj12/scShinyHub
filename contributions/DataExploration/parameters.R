
# sub menue items for the parameter tab
# myPparameters = list(
#   menuSubItem("Normalization2", tabName = "mynormalizations")
# )

# tab content

# tabList = list(
#   tabItem("mynormalizations",list(
#     tags$h3("Test"),
#     fluidRow(column(10,
#                     tags$h4("h4 text as test")
#                     # 10, offset = 1,
#                     # plotOutput('plotUmiHist') %>% withSpinner()
#     ))
#   )
#   )
# )

# normalization parameters

# choice for the radio buttion
myNormalizationChoices <- list(
  scEx_log = "scEx_logNormalization"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  scEx_log = h4("no Parameters implemented")
)

scEx_logNormalization <- reactive({
  # TODO ?? define scaling factor somewhere else???
  scalingFactor = 10000
  if (DEBUG) {
    cat(file = stderr(), "scEx_logNormalization\n")
  }
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
  
  # use_genes <- get_nonzero_genes(scEx)
<<<<<<< HEAD
  use_genes <- sort(unique(1 + slot(as(assays(scEx)[["counts"]], "dgTMatrix"), 
=======
  use_genes <- sort(unique(1 + slot(as(assays(scEx)[[1]], "dgTMatrix"), 
>>>>>>> 15ba2245b451381ee5096a6fd814dfdefef83320
                                    "i")))
  
  bc_sums <- Matrix::colSums(assays(scEx)[["counts"]])
  median_sum <- median(bc_sums)
  A <- as(assays(scEx)[["counts"]], "dgCMatrix")
  A@x <- A@x / Matrix::colSums(A)[assays(scEx)[["counts"]]@j + 1L]
  # new_matrix <- sweep(exprs(scEx), 2, median_sum/bc_sums, "*")
  # scEx_bcnorm <- (newGeneBCMatrix(A, fData(scEx), pData(scEx),
  #   template = scEx
  # ))
<<<<<<< HEAD
  scEx_bcnorm <- SingleCellExperiment(assay = list(logcounts = as(A,"dgTMatrix")),
=======
  scEx_bcnorm <- SingleCellExperiment(assay = as(A,"dgTMatrix"),
>>>>>>> 15ba2245b451381ee5096a6fd814dfdefef83320
                                      colData = colData(scEx),
                                      rowData = rowData(scEx))
  
  # gbm_bcnorm <- normalize_barcode_sums_to_median(gbm)
  # gbm_log <- log_gene_bc_matrix(gbm_bcnorm, base = 10)
  x <- uniqTsparse(assays(scEx_bcnorm)[[1]])
  slot(x, "x") <- log(1 + slot(x, "x"), base = 2) * scalingFactor
  assays(scEx_bcnorm)[[1]] <- x
  
  if (DEBUG) {
    cat(file = stderr(), "scEx_logNormalization:Done\n")
  }
  return(scEx_bcnorm)
})
