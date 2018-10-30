
# # sub menue items for the parameter tab
# myPparameters = list(
#   menuSubItem("Normalization2", tabName = "mynormalizations")
# )
# 
# # tab content
# 
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
myNormalizationChoices = list(gbm_log = "gbm_logNormalization" #,
                              #equalize_gbms = "equalize_gbms"
                              )

# value should be of class shiny.tag
myNormalizationParameters = list(gbm_log = h4("no Parameters implemented")
                                 #,
                                 #equalize_gbms = h4("no Parameters implemented")
                                 
)
                                 
# Seems to be depricated
#
# equalize_gbms <- reactive({
#   if (DEBUG)
#     cat(file = stderr(), "equalize_gbms\n")
#   gbm = gbm()
#   sampleInfo = inputSample()
#   
#   if (is.null(gbm) | is.null(sampleInfo)) {
#     if (DEBUG)
#       cat(file = stderr(), "equalize_gbms:NULL\n")
#     return(NULL)
#   }
#   if (DEBUGSAVE){
#     cat(file = stderr(), "equalize_gbms:saving\n")
#     save(file = "~/scShinyHubDebug/equalize_gbms", list = c(ls(),ls(envir = globalenv())))
#     cat(file = stderr(), "equalize_gbms:saving done\n")
#     exit()
#   }
#     # load(file="~/scShinyHubDebug/equalize_gbms")
#   
#   
#   if (length(levels(sampleInfo$sample))>1){
#     gbm_list = list()
#     for (smpLvl in levels(sampleInfo$sample)){
#       gbm_list[[length(gbm_list)+1]] = gbm[sampleInfo$sample == smpLvl]
#     }
#     gbm_list <- lapply(gbm_list,load_molecule_info)
#   }
#   # set.seed(0)
#   # gbm_list <- list(gbm1, gbm2)
#   # gbm_list <- lapply(gbm_list,load_molecule_info) # load sample molecule information
#   # gbm_list_equalized <- equalize_gbms(gbm_list) # equalize the gene-barcode matrices
#   # merged_gbm <- concatenate_gene_bc_matrices(gbm_list_equalized)
#   
#   if (DEBUG)
#     cat(file = stderr(), "gbm_logNormalization:Done\n")
#   return(gbm)
#   
# })

gbm_logNormalization <- reactive({
  if (DEBUG)
    cat(file = stderr(), "gbm_logNormalization\n")
  gbm = gbm()

  if (is.null(gbm)) {
    if (DEBUG)
      cat(file = stderr(), "gbm_logNormalization:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/gbm_logNormalization.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/gbm_logNormalization.RData")
  
  use_genes <- get_nonzero_genes(gbm)
  gbm_bcnorm <- normalize_barcode_sums_to_median(gbm)
  gbm_log <- log_gene_bc_matrix(gbm_bcnorm, base = 10)

  if (DEBUG)
    cat(file = stderr(), "gbm_logNormalization:Done\n")
  return(gbm_log)
  
  
})








