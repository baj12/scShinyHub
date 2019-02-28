require(uwot)
require(tidyverse)
require(SingleCellExperiment)

# here we define reactive values/variables

# scaterReadsFunc <- function(gbm, gbm_log, fd){
scaterReadsFunc <- function(gbm, fd) {
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/scaterReadsFunc.Rmd", list = c(ls()))
  }
  # load(file='~/scShinyHubDebug/scaterReadsFunc.Rmd')

  counts <- as.matrix(exprs(gbm))

  anno <- pData(gbm)
  anno$sample_id <- anno$sampleNames
  anno$fixed <- "red"
  # anno$individual= "NA1"
  # anno$replicate = "r1"
  # anno$well = "A01"
  # anno$batch = "b1"
  pheno_data <- new("AnnotatedDataFrame", anno)
  # rownames(pheno_data) <- pheno_dat

  reads <- as.matrix(counts)
  rownames(reads) <- make.unique(fd[rownames(reads), "Associated.Gene.Name"], sep = "___")
  rownames(reads)[is.na(rownames(reads)) ] <- "na"
  reads <- SingleCellExperiment(
    assays = list(counts = reads),
    colData = anno
  )
  ercc <- rownames(reads)[grepl("ERCC-", rownames(reads))]

  mt <- rownames(reads)[grepl("^MT", rownames(reads))]

  reads <- scater::calculateQCMetrics(
    reads
  )
  filter_by_expr_features <- (reads$total_features_by_counts > 200)
  reads$use <- (
    # sufficient features (genes)
    filter_by_expr_features
    # sufficient molecules counted
    # filter_by_total_counts &
    # sufficient endogenous RNA
    # filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    # filter_by_MT
  )
  return(reads)
}

scaterReads <- reactive({
  if (DEBUG) cat(file = stderr(), "scaterReads\n")
  gbm <- gbm()
  # gbm_log = gbm_log()
  fd <- featureDataReact()
  if (is.null(gbm) | is.null(gbm_log)) {
    return(NULL)
  }
  # return(scaterReadsFunc(gbm, gbm_log, fd))
  return(scaterReadsFunc(gbm, fd))
})


sampleHistFunc <- function(samples) {
  counts <- table(samples)
  barplot(counts,
    main = "histogram of number of cell per sample",
    xlab = "Samples"
  )
  # x <- hist(as.integer(as.factor(samples)),
  #   main = "histogram of number of cell per sample",
  #   labels = levels(as.factor(samples)),
  #   breaks = 0:length(levels(as.factor(samples))),
  #   xlab = "Samples"
  # )
}

inputTSNESample <- reactive({
  on.exit(
    removeNotification(id = "inputTSNESample")
  )
  if (DEBUG) cat(file = stderr(), "inputTSNESample\n")
  projections <- projections()
  if (is.null(projections)) {
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("inputTSNESample", id = "inputTSNESample", duration = NULL)
  }

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/inputTSNESample.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file = "~/scShinyHubDebug/inputTSNESample.RData")

  return(projections)
})

# TODO separate  function from reactive : done? run_tsne is already the function.
# Maybe we need a normalized name like tsneFunc?
tsne <- reactive({
  on.exit(
    removeNotification(id = "tsne")
  )
  if (DEBUG) {
    cat(file = stderr(), "tsne\n")
  }
  pca <- pca()
  tsneDim <- input$tsneDim
  tsnePerplexity <- input$tsnePerplexity
  tsneTheta <- input$tsneTheta
  tsneSeed <- input$tsneSeed
  if (is.null(pca)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne", id = "tsne", duration = NULL)
  }
  set.seed(seed = tsneSeed)
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/tsne.RData')
  require(parallel)

  retval <- tryCatch({
    run_tsne(
      pca,
      dims = tsneDim,
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
  if (DEBUG) {
    cat(file = stderr(), "tsne: done\n")
  }
  return(retval)
})
tsne1 <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "tsne1\n")
  }
  tsne.data <- tsne.data()
  if (is.null(tsne.data)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne1: NULL\n")
    }
    return(NULL)
  }
  return(tsne.data$tsne1)
})
tsne2 <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "tsne1\n")
  }
  tsne.data <- tsne.data()
  if (is.null(tsne.data)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne2: NULL\n")
    }
    return(NULL)
  }
  return(tsne.data$tsne2)
})
tsne3 <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "tsne1\n")
  }
  tsne.data <- tsne.data()
  if (is.null(tsne.data)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne3: NULL\n")
    }
    return(NULL)
  }
  return(tsne.data$tsne3)
})
tsne4 <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "tsne1\n")
  }
  tsne.data <- tsne.data()
  if (is.null(tsne.data)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne4: NULL\n")
    }
    return(NULL)
  }
  return(tsne.data$tsne4)
})
tsne5 <- reactive({
  if (DEBUG) {
    cat(file = stderr(), "tsne1\n")
  }
  tsne.data <- tsne.data()
  if (is.null(tsne.data)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne5: NULL\n")
    }
    return(NULL)
  }
  return(tsne.data$tsne5)
})


tsne.data <- reactive({
  on.exit(
    removeNotification(id = "tsne.data")
  )
  if (DEBUG) {
    cat(file = stderr(), "tsne.data\n")
  }
  tsne <- tsne()
  if (is.null(tsne)) {
    if (DEBUG) {
      cat(file = stderr(), "tsne.data: NULL\n")
    }
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("tsne.data", id = "tsne.data", duration = NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/tsne.data.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file="~/scShinyHubDebug/tsne.data.RData")
  tsne.data <- data.frame(tsne$Y)
  colnames(tsne.data) <- paste0("tsne", c(1:ncol(tsne.data)))

  if (DEBUG) {
    cat(file = stderr(), "tsne.data: done\n")
  }
  return(tsne.data)
})

umapReact <- reactive({
  gbmlog <- gbm_log()
  
  start.time <- base::Sys.time()
  set.seed(input$um_randSeed)
  xaxis <- input$um_xaxis
  yaxis <- input$um_yaxis
  cellT <- input$um_ct
  inputCT <- input$um_inputCT
  sampleRatio <- as.numeric(input$um_sampleRatio)
  sampleIds <- input$um_sampleIds
  InlevelOrd <- input$um_levelOrd
  UMAP1 <- input$um_umap1
  UMAP2 <- input$um_umap2
  
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
  
  if (is.null(gbmlog)) {
    if (DEBUG) cat(file = stderr(), "output$umap_react:NULL\n")
    return(NULL)
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/umap_react.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load("~/scShinyHubDebug/umap_react.RData")
  
  umapData <- as.matrix(exprs(gbmlog))
  compCases = complete.cases(umapData)
  
  # TODO it might be possible to reuse nearest neighbor information to speeed up recomputations
  # with eg. new seed
  
  embedding <- uwot::umap(t(as.matrix(exprs(gbmlog))),
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
  rownames(embedding) = colnames(gbmlog)
  
  end.time <- Sys.time()
  if (DEBUG)
    cat(file = stderr(), paste("umap took: ", difftime(end.time, start.time, units = "min"), " min\n"))
  return(embedding)
})



myProjections <- list(
  c("tsne1", "tsne1"),
  c("tsne2", "tsne2"),
  c("tsne3", "tsne3"),
  c("tsne4", "tsne4"),
  c("tsne5", "tsne5"),
  c("dbCluster", "dbCluster"),
  c("umap", "umapReact")
)


# declare function as heavy
myHeavyCalculations <- list(
  c("scaterReads", "scaterReads"),
  c("tsne", "tsne")
)
