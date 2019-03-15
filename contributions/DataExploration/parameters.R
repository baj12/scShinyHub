
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
  gbm_log = "gbm_logNormalization"
)

# value should be of class shiny.tag
# will be displayed via renderUI
myNormalizationParameters <- list(
  gbm_log = h4("no Parameters implemented")
)

gbm_logNormalization <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "gbm_logNormalization\n")
  }
  countES <- countES()

  if (is.null(countES)) {
    if (DEBUG) {
      cat(file = stderr(), "gbm_logNormalization:NULL\n")
    }
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/gbm_logNormalization.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/gbm_logNormalization.RData")

  use_genes <- get_nonzero_genes(countES)
  bc_sums <- Matrix::colSums(countES)
  median_sum <- median(bc_sums)
  A <- as(exprs(countES), "dgCMatrix")
  A@x <- A@x / Matrix::colSums(A)[exprs(countES)@j + 1L]
  # new_matrix <- sweep(exprs(gbm), 2, median_sum/bc_sums, "*")
  gbm_bcnorm <- (newGeneBCMatrix(A, fData(countES), pData(countES),
    template = countES
  ))
  # gbm_bcnorm <- normalize_barcode_sums_to_median(gbm)
  gbm_log <- log_gene_bc_matrix(gbm_bcnorm, base = 10)

  if (DEBUG) {
    cat(file = stderr(), "gbm_logNormalization:Done\n")
  }
  return(gbm_log)
})
