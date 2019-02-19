# normalization parameters

# choice for the radio buttion
# cat(file = stderr(), "parameters1\n")
#     myNormalizationChoices <- list(dca_log = "dca_impute" #,
#                                   #equalize_gbms = "equalize_gbms"
#     )
tryCatch(
  {
    # save(file = "test.RData", list = c(ls(), ls(envir = .GlobalEnv)))
#    load(file = "test.RData")
    system("dca -h > /dev/null")
    myNormalizationChoices <- list(dca_log = "dca_impute" #,
                                  #equalize_gbms = "equalize_gbms"
    )

    # value should be of class shiny.tag
    myNormalizationParameters <- list(dca_log = h4("no Parameters implemented")
                                     #,
                                     #equalize_gbms = h4("no Parameters implemented")

    )
  },
  warning = function(e){
    print("warning")
    myNormalizationChoices <<- c()
    myNormalizationParameters <<- list()
  },
  error = function(e){
    print("error")
    myNormalizationChoices <<- c()
    myNormalizationParameters <<- list()
  }
)

# cat(file = stderr(), "parameters2\n")


dca_impute <- reactive({
  if (DEBUG)
    cat(file = stderr(), "dca_impute\n")
  gbm_matrix <- gbm_matrix()
  gbm <- gbm()
  # load("/Users/bernd/scShinyHubDebug/heatmapSOMReactive.RData")
  
  if (is.null(gbm_matrix)) {
    if (DEBUG)
      cat(file = stderr(), "dca_impute:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE)
    save(file = "~/scShinyHubDebug/dca_impute.RData", list = c(ls(),ls(envir = globalenv())))
  # load(file="~/scShinyHubDebug/dca_impute.RData")
  
  tfile = tempfile(pattern = "dcaInput", tmpdir = tempdir(), fileext = ".csv")
  tdir = paste0(tempdir(),"/dcaresults")
  resFile = paste0(tdir,"/mean_norm.tsv")
  write.csv(x = gbm_matrix, file = tfile)
  system(paste("dca", tfile, tdir))
  if (!file.exists(resFile)) {
      cat(file = stderr(), "dca_impute:failed\n")
    return(NULL)
  }
  newmat = read.csv(resFile, sep = "\t", row.names = 1)

  colnames(newmat) = colnames(gbm_matrix)
  colnames(gbm_matrix) %in% colnames(newmat)
  gbm_matrix[rownames(newmat), colnames(newmat)] = as.matrix(newmat[rownames(newmat), colnames(newmat)])
  gbm_matrix[!rownames(newmat) %in% rownames(gbm_matrix), ] = 0
  if (DEBUG)
    cat(file = stderr(), "dca_impute:Done\n")
  
  newgbmLog = newGeneBCMatrix(mat = as(gbm_matrix, "dgTMatrix"), fd = fData(gbm), pd = pData(gbm))
  return(newgbmLog)
  
  
})



cat(file = stderr(), "parameters (Theislab) end\n")

