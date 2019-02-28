selectedDge <- reactiveValues(
  dgeTable = data.frame()
)


# TODO remove log2cpm change to gbm_log
dge_func <- function(projections, log2cpm, featureData, dbCluster, cl1, db1, db2) {
  subsetData <- subset(projections, dbCluster %in% cl1)
  cells.1 <- rownames(shiny::brushedPoints(subsetData, db1))
  cells.2 <- rownames(shiny::brushedPoints(subsetData, db2))

  subsetExpression <- log2cpm[, union(cells.1, cells.2)]

  genes.use <- rownames(subsetExpression)
  data.1 <- apply(subsetExpression[genes.use, cells.1], 1, expMean)
  data.2 <- apply(subsetExpression[genes.use, cells.2], 1, expMean)
  total.diff <- (data.1 - data.2)

  genes.diff <- names(which(abs(total.diff) > .2))
  genes.use <- ainb(genes.use, genes.diff)

  toReturn <-
    DiffExpTest(subsetExpression, cells.1, cells.2, genes.use = genes.use)
  toReturn[, "avg_diff"] <- total.diff[rownames(toReturn)]
  toReturn$Associated.Gene.Name <-
    featureData[rownames(toReturn), "Associated.Gene.Name"]
  return(toReturn)
}


dge <- reactive({
  on.exit(
    if (!is.null(getDefaultReactiveDomain()))
      removeNotification(id = "dge")
  )
  if (DEBUG) cat(file = stderr(), "dge\n")
  featureData <- featureDataReact()
  log2cpm <- log2cpm()
  prj <- projections()
  gn <- groupNames$namesDF

  cl1 <- input$clusters1
  db1 <- input$db1
  db2 <- input$db2

  # dbcl = dbCluster
  if (is.null(featureData) | is.null(log2cpm) | is.null(prj)) {
    return(NULL)
  }
  if (!is.null(getDefaultReactiveDomain())) {
    showNotification("running dge", id = "dge", duration = NULL)
  }
  if (length(gn) > 0) {
    prj <- cbind(prj, gn[rownames(prj), ] * 1)
  }

  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/dge.RData", list = c(ls(), ls(envir = globalenv())))
  }
  # load(file='~/scShinyHubDebug/dge.RData')

  toReturn <- dge_func(
    projections = prj, log2cpm = log2cpm,
    featureData = featureData,
    dbCluster = prj$dbCluster, cl1 = cl1, db1 = db1, db2 = db2
  )
  if (DEBUG) cat(file = stderr(), "dge14\n")
  if (nrow(toReturn) == 0) {
    if (DEBUG) cat(file = stderr(), "dge: nothing found\n")
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification("dge: nothing found", id = "dgewarning", duration = 10, type = "warning")
    }
  }
  if (DEBUG) cat(file = stderr(), "dge15\n")
  selectedDge$dgeTable <- toReturn
  cat(stderr(), rownames(toReturn)[1:5])
  if (DEBUG) cat(file = stderr(), "dge: done\n")

  return(toReturn)
})
