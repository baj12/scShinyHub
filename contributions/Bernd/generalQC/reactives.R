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
  rownames(reads) <- make.unique(fd[rownames(reads), "Associated.Gene.Name"])
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
  barplot(counts, main="histogram of number of cell per sample", 
          xlab="Samples")
  # x <- hist(as.integer(as.factor(samples)),
  #   main = "histogram of number of cell per sample",
  #   labels = levels(as.factor(samples)),
  #   breaks = 0:length(levels(as.factor(samples))),
  #   xlab = "Samples"
  # )
}

inputTSNESample <- reactive({
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

  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "inputTSNESample")
  }
  return(projections)
})

# TODO separate  function from reactive : done? run_tsne is already the function.
# Maybe we need a normalized name like tsneFunc?
tsne <- reactive({
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
  retval <- tryCatch({
    run_tsne(
      pca,
      dims = tsneDim,
      perplexity = tsnePerplexity,
      theta = tsneTheta,
      check_duplicates = FALSE
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
  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "tsne")
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

  if (!is.null(getDefaultReactiveDomain())) {
    removeNotification(id = "tsne.data")
  }
  if (DEBUG) {
    cat(file = stderr(), "tsne.data: done\n")
  }
  return(tsne.data)
})


myProjections <- list(
  c("tsne1", "tsne1"),
  c("tsne2", "tsne2"),
  c("tsne3", "tsne3"),
  c("tsne4", "tsne4"),
  c("tsne5", "tsne5"),
  c("dbCluster", "dbCluster")
)


# declare function as heavy
myHeavyCalculations <- list(
  c("scaterReads", "scaterReads"),
  c("tsne", "tsne")
)