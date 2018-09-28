
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
myNormalizationChoices = list(
                              scater_norm = "scater_norm")

# value should be of class shiny.tag
myNormalizationParameters = list(
  scater_norm = h4("no Parameters implemented")
)


scater_norm <- reactive({
  if (DEBUG)
    cat(file = stderr(), "scaterProtocol\n")
  gbm = gbm()
  sampleInfo = inputSample()
  scaterReads = scaterReads()
  minExpression = ncol(scaterReads)/10
  
  require(scater)
  require(scran)
  
  if (is.null(gbm) | is.null(sampleInfo)) {
    if (DEBUG)
      cat(file = stderr(), "scater_norm:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE){
    cat(file = stderr(), "scater_norm:saving\n")
    save(file = "~/scShinyHubDebug/scater_norm.RData", list = c(ls(),ls(envir = globalenv())))
    cat(file = stderr(), "scater_norm:saving done\n")
    # exit()
  }
  # load(file="~/scShinyHubDebug/scater_norm.RData")
  t(as.matrix(table(sampleInfo$sample)))[1,]
  cdata = colData(scaterReads)
  cdata$Treatment = sampleInfo$sample
  colData(scaterReads) = cdata
  filter_by_expr_features <- (scaterReads$total_features_by_counts > minExpression)
  scaterReads$use <- (
    filter_by_expr_features 
    )
  
  scaterReads <- scran::computeSumFactors(scaterReads, sizes = seq(21, min(table(sampleInfo$sample)), 5), 
                                          clusters = sampleInfo$sample, subset.row =scaterReads$use)
  # summary(sizeFactors(scaterReads))
  plot(sizeFactors(scaterReads), scaterReads$total_counts/1e6, log="xy",
       ylab="Library size (millions)", xlab="Size factor")
  scaterReads <- normalize(scaterReads)
  plotExplanatoryVariables(scaterReads, variables=c("total_features_by_counts",
                                            "log10_total_features_by_counts", "pct_counts_in_top_100_features"))
  # # if (length(levels(sampleInfo$sample))>1){
  #   gbm_list = list()
  #   for (smpLvl in levels(sampleInfo$sample)){
  #     gbm_list[[length(gbm_list)+1]] = gbm[sampleInfo$sample == smpLvl]
  #   }
  #   gbm_list <- lapply(gbm_list,load_molecule_info)
  # }
  # set.seed(0)
  # gbm_list <- list(gbm1, gbm2)
  # gbm_list <- lapply(gbm_list,load_molecule_info) # load sample molecule information
  # gbm_list_equalized <- equalize_gbms(gbm_list) # equalize the gene-barcode matrices
  # merged_gbm <- concatenate_gene_bc_matrices(gbm_list_equalized)
  
  if (DEBUG)
    cat(file = stderr(), "gbm_logNormalization:Done\n")
  retVal = SummarizedExperiment::assays(scaterReads)$logcounts
  rownames(retVal) = rownames(gbm)
  retVal = newGeneBCMatrix(retVal, pd=pData(gbm), fd=fData(gbm))
  return(retVal)
  
})







