require(uwot)
require(tidyverse)
require(SingleCellExperiment)

# here we define reactive values/variables

# scaterReadsFunc ----
#' scaterReadsFunc
#' calculate the QC metrix and return updated singleCellExperiment object
scaterReadsFunc <- function(scEx) {

  if (class(assays(scEx)[["counts"]]) == "dgTMatrix") {
    assays(scEx)[["counts"]] = as(assays(scEx)[["counts"]], "dgCMatrix")
  }
  
  ercc <- rownames(scEx)[grepl("ERCC-", rownames(scEx))]
  
  mt <- rownames(scEx)[grepl("^MT", rownames(scEx))]
  
  scEx <- scater::calculateQCMetrics(
    scEx
  )
  filter_by_expr_features <- (scEx$total_features_by_counts > 200)
  scEx$use <- (
    # sufficient features (genes)
    filter_by_expr_features
    # sufficient molecules counted
    # filter_by_total_counts &
    # sufficient endogenous RNA
    # filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    # filter_by_MT
  )
  
  return(scEx)
}

# scaterReads ----
#' scaterReads
#' singleCellExperiment object/reactive with QC metrix
scaterReads <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "scaterReads")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("scaterReads", id = "scaterReads", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "scaterReads\n")
  
  scEx <- scEx()
  # scEx_log = scEx_log()
  if (is.null(scEx)) {
    return(NULL)
  }
  retVal <- scaterReadsFunc(scEx)

  printTimeEnd(start.time, "scaterReads")
  exportTestValues(scaterReads = {str(retVal)})  
  return(retVal)
})


# sampleHistFunc ----
#' sampleHistFunc
#' create a histogram from samples
sampleHistFunc <- function(samples, scols) {
  counts <- table(samples)
  barplot(counts,
          main = "histogram of number of cell per sample",
          xlab = "Samples",
          col=scols
  )
}


# projectionTable -----
#' projectionTable
#' input for tableSelectionServer
#' presents all projections
projectionTable <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "projectionTable")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("projectionTable", id = "projectionTable", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "projectionTable\n")
  
  projections <- projections()
  
  if (is.null(projections)) {
    return(NULL)
  }

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/projectionTable.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/projectionTable.RData")
  
  printTimeEnd(start.time, "projectionTable")
  exportTestValues(projectionTable = {projections})  
  return(projections)
})

# tsne ----
#' tsne
#' reactive calculating the tSNE projections
tsne <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "tsne")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne", id = "tsne", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "tsne\n")
  
  pca <- pca()
  tsneDim <- input$tsneDim
  tsnePerplexity <- input$tsnePerplexity
  tsneTheta <- input$tsneTheta
  tsneSeed <- input$tsneSeed
  
  if (is.null(pca)) {
    if (DEBUG) cat(file = stderr(), "tsne: NULL\n")
    return(NULL)
  }
  
  retVal <- tsneFunc(pca, tsneDim, tsnePerplexity, tsneTheta, tsneSeed)

  printTimeEnd(start.time, "tsne")
  exportTestValues(tsne = {retVal})  
  return(retVal)
})

tsneFunc <- function(pca, tsneDim, tsnePerplexity, tsneTheta, tsneSeed) {
  set.seed(seed = tsneSeed)
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/tsne.RData')
  require(parallel)
  require(Rtsne)
  np = dim(pca$x)[2]
  tsne <- tryCatch({
    Rtsne::Rtsne(
      pca$x[,1:np], pca = FALSE, dims = tsneDim,
      perplexity = tsnePerplexity,
      theta = tsneTheta,
      check_duplicates = FALSE, num_threads = detectCores()
    )
  },
  error = function(e) {
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(paste("Problem with tsne:", e),
                       type = "error",
                       duration = NULL
      )
    }
    return(NULL)
  }
  )
  retVal <- data.frame(tsne$Y)
  colnames(retVal) <- paste0("tsne", c(1:ncol(retVal)))
  return(retVal)
}


# umapReact ----
#' umapReact
#' reactive for calculating UMAP projection
umapReact <- reactive({
  start.time <- base::Sys.time()
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "umapReact")
  )
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("umapReact", id = "umapReact", duration = NULL)
  }
  if (DEBUG) cat(file = stderr(), "umapReact started.\n")
  
  scEx_log <- scEx_log()
  myseed <- input$um_randSeed
  xaxis <- input$um_xaxis
  yaxis <- input$um_yaxis
  cellT <- input$um_ct
  inputCT <- input$um_inputCT
  sampleRatio <- as.numeric(input$um_sampleRatio)
  sampleIds <- input$um_sampleIds
  InlevelOrd <- input$um_levelOrd
  UMAP1 <- input$um_umap1
  UMAP2 <- input$um_umap2
  runUMAP <- input$activateUMAP
  
  n_neighbors <- as.numeric(input$um_n_neighbors)
  n_components <- as.numeric(input$um_n_components)
  n_epochs <- as.numeric(input$um_n_epochs)
  alpha <- as.numeric(input$um_alpha)
  init <- input$um_init
  min_dist <- as.numeric(input$um_min_dist)
  set_op_mix_ratio <- as.numeric(input$um_set_op_mix_ratio)
  local_connectivity <- as.numeric(input$um_local_connectivity)
  bandwidth <- as.numeric(input$um_bandwidth)
  gamma <- as.numeric(input$um_gamma)
  negative_sample_rate <- as.numeric(input$um_negative_sample_rate)
  metric <- input$um_metric
  spread <- as.numeric(input$um_spread)
  
  if (is.null(scEx_log)) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/umap_react.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load("~/scShinyHubDebug/umap_react.RData")
  if (!runUMAP) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  umapData <- as.matrix(assays(scEx_log)[[1]])
  compCases = complete.cases(umapData)
  
  # TODO it might be possible to reuse nearest neighbor information to speeed up recomputations
  # with eg. new seed
  
  set.seed(myseed)
  embedding <- uwot::umap(t(as.matrix(assays(scEx_log)[[1]])),
                          n_neighbors = n_neighbors,
                          n_components = n_components, n_epochs = n_epochs,
                          # alpha = alpha,
                          init = init,
                          spread = spread,
                          min_dist = min_dist,
                          set_op_mix_ratio = set_op_mix_ratio,
                          local_connectivity = local_connectivity,
                          bandwidth = bandwidth,
                          # gamma = gamma,
                          negative_sample_rate = negative_sample_rate,
                          metric = metric,
                          n_threads = detectCores()
  )
  embedding = as.data.frame(embedding)
  colnames(embedding) = paste0("UMAP", 1:n_components)
  rownames(embedding) = colnames(scEx_log)
  
  printTimeEnd(start.time, "umapReact")
  exportTestValues(umapReact = {embedding})  
  return(embedding)
})


# myProjections ----
myProjections <- list(
  c("tsne", "tsne"),
  c("dbCluster", "dbCluster"),
  c("umap", "umapReact")
)

# myHeavyCalculations ----
# declare function as heavy
myHeavyCalculations <- list(
  c("scaterReads", "scaterReads"),
  c("tsne", "tsne")
)
