geneName2Index <- function(g_id, featureData) {
  if (DEBUG) {
    cat(file = stderr(), paste("geneName2Index\n"))
  }
  if(is.null(g_id)){
    return(NULL)
  }
  
  g_id <- toupper(g_id)
  g_id <- gsub(" ", "", g_id, fixed = TRUE)
  g_id <- strsplit(g_id, ",")
  g_id <- g_id[[1]]

  notFound <- g_id[!g_id %in% toupper(featureData$symbol)]
  if (length(featureData$symbol) == length(notFound)) {
    # in case there is only one gene that is not available.
    notFound <- g_id
  }
  if (length(notFound) > 0) {
    if (DEBUG) {
      cat(file = stderr(), paste("gene names not found: ", notFound, "\n"))
    }
    if (!is.null(getDefaultReactiveDomain())) {
      showNotification(paste("following genes were not found", notFound, collapse = " "),
        id = "moduleNotFound", type = "warning",
        duration = 20
      )
    }
  }

  geneid <- rownames(featureData[which(toupper(featureData$symbol) %in% toupper(g_id)), ])
  if (DEBUG) {
    cat(file = stderr(), paste("done: geneName2Index\n"))
  }
  return(geneid)
}


updateProjectionsWithUmiCount <- function(dimX, dimY, geneNames, geneNames2 = NULL, featureData, scEx, projections) {
  if ((dimY == "UmiCountPerGenes") | (dimX == "UmiCountPerGenes")) {
    geneNames <- geneName2Index(geneNames, featureData)
    if ((length(geneNames) > 0) & (length(geneNames[[1]]) > 0)) {
      if (length(geneNames) == 1) {
        projections$UmiCountPerGenes <- assays(scEx)[["counts"]][geneNames, ]
      } else {
        projections$UmiCountPerGenes <- Matrix::colSums(assays(scEx)[["counts"]][geneNames, ])
      }
    } else {
      projections$UmiCountPerGenes <- 0
    }
  }
  if ((dimY == "UmiCountPerGenes2") | (dimX == "UmiCountPerGenes2")) {
    geneNames <- geneName2Index(geneNames2, featureData)
    if ((length(geneNames) > 0) & (length(geneNames[[1]]) > 0)) {
      if (length(geneNames) == 1) {
        projections$UmiCountPerGenes2 <- assays(scEx)[["counts"]][geneNames, ]
      } else {
        projections$UmiCountPerGenes2 <- Matrix::colSums(assays(scEx)[["counts"]][geneNames, ])
      }
    } else {
      projections$UmiCountPerGenes2 <- 0
    }
  }
  return(projections)
}


# append to heavyCalculations
appendHeavyCalculations <- function(myHeavyCalculations, heavyCalculations) {
  for (hc in myHeavyCalculations) {
    if (length(hc) == 2 & is.character(hc[1]) & is.character(hc[2])) {
      heavyCalculations[[length(heavyCalculations) + 1]] <- hc
    } else {
      stop(paste("myHeavyCalculations is malformatted\nhc:", hc, "\nmyHeavyCalculations:", myHeavyCalculations, "\n"))
    }
  }
  return(heavyCalculations)
}

#### plot2Dprojection ----------------
# used in moduleServer and reports
plot2Dprojection <- function(scEx_log, scEx, projections, g_id, featureData,
                             geneNames, geneNames2, dimX, dimY, clId, grpN, legend.position, grpNs,
                             logx = FALSE, logy = FALSE, divXBy = "None", divYBy = "None") {
  geneid <- geneName2Index(g_id, featureData)
  
  if (length(geneid) == 0) {
    return(NULL)
  }
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ,drop=FALSE]
  # } else {
  expression <- Matrix::colSums(assays(scEx_log)[[1]][geneid, , drop = FALSE])
  # }
  validate(need(is.na(sum(expression)) != TRUE, ""))
  # if (length(geneid) == 1) {
  #   expression <- exprs(scEx_log)[geneid, ]
  # } else {
  #   expression <- Matrix::colSums(exprs(scEx_log)[geneid, ])
  # }
  # validate(need(is.na(sum(expression)) != TRUE, ""))
  
  # geneid <- geneName2Index(geneNames, featureData)
  projections <- updateProjectionsWithUmiCount(
    dimX = dimX, dimY = dimY,
    geneNames = geneNames,
    geneNames2 = geneNames2,
    featureData = featureData,
    scEx = scEx, projections = projections
  )
  
  
  projections <- cbind(projections, expression)
  names(projections)[ncol(projections)] <- "exprs"
  
  if (DEBUG) {
    cat(file = stderr(), paste("output$dge_plot1:---", clId[1], "---\n"))
  }
  subsetData <- subset(projections, dbCluster %in% clId)
  # subsetData$dbCluster = factor(subsetData$dbCluster)
  # if there are more than 18 samples ggplot cannot handle different shapes and we ignore the
  # sample information
  if (length(as.numeric(as.factor(subsetData$sample))) > 18) {
    subsetData$shape <- as.factor(1)
  } else {
    subsetData$shape <- as.numeric(as.factor(subsetData$sample))
  }
  if (DEBUGSAVE) {
    save(file = "~/scShinyHubDebug/clusterPlot.RData", list = c(ls(), "legend.position", ls(envir = globalenv())))
    cat(file = stderr(), paste("plot2Dprojection saving done.\n"))
  }
  # load(file="~/scShinyHubDebug/clusterPlot.RData")
  if (nrow(subsetData) == 0) return(NULL)
  # subsetData$shape = as.factor(1)
  gtitle <- paste(toupper(g_id), clId, sep = "-Cluster", collapse = " ")
  if (nchar(gtitle) > 50) {
    gtitle <- paste(substr(gtitle, 1, 50), "...")
  }
  
  require(plotly)
  f <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  if (divXBy != "None") {
    subsetData[, dimX] <- subsetData[, dimX] / subsetData[, divXBy]
  }
  if (divYBy != "None") {
    subsetData[, dimY] <- subsetData[, dimY] / subsetData[, divYBy]
  }
  
  typeX <- typeY <- "linear"
  if (logx) {
    typeX <- "log"
  }
  if (logy) {
    typeY <- "log"
  }
  if (is.factor(subsetData[, dimX])) {
    typeX <- NULL
  }
  if (is.factor(subsetData[, dimY])) {
    typeY <- NULL
  }
  xAxis <- list(
    title = dimX,
    titlefont = f,
    type = typeX
  )
  yAxis <- list(
    title = dimY,
    titlefont = f,
    type = typeY
  )
  p1 <- plot_ly(data = subsetData, source = "subset") %>%
    add_trace(
      x = ~ get(dimX),
      y = ~ get(dimY),
      type = "scatter",
      mode = "markers",
      text = ~ paste(1:nrow(subsetData), " ", rownames(subsetData), "<br />", subsetData$exprs),
      marker = list(
        color = subsetData[, "exprs"],
        size = 4
      )
    ) %>%
    # add_trace() %>%
    layout(
      xaxis = xAxis,
      yaxis = yAxis,
      title = gtitle,
      dragmode = "select"
    )
  
  selectedCells <- NULL
  if (length(grpN) > 0) {
    if (length(grpNs[rownames(subsetData), grpN]) > 0 & sum(grpNs[rownames(subsetData), grpN], na.rm = TRUE) > 0) {
      grpNSub <- grpNs[rownames(subsetData), ]
      selectedCells <- rownames(grpNSub[grpNSub[, grpN], ])
    }
  }
  if (!is.null(selectedCells)) {
    shape <- rep("a", nrow(subsetData))
    selRows <- which(rownames(subsetData) %in% selectedCells)
    shape[selRows] <- "b"
    x1 <- subsetData[selectedCells, dimX, drop = FALSE]
    y1 <- subsetData[selectedCells, dimY, drop = FALSE]
    p1 <- p1 %>%
      add_trace(
        x = x1[, 1], y = y1[, 1],
        marker = list(
          color = rep("green", nrow(x1)),
          size = 5
        ),
        text = ~ paste(
          rownames(subsetData[selectedCells, ]),
          "<br />", subsetData$exprs[selectedCells]
        ),
        type = "scatter",
        mode = "markers",
        name = "selected"
      )
  }
  p1
}


# functions should go in external file

n_fun <- function(x) {
  return(data.frame(y = -0.5, label = paste0(length(x), "\ncells")))
}

diffLRT <- function(x, y, xmin = 1) {
  lrtX <- bimodLikData(x)
  lrtY <- bimodLikData(y)
  lrtZ <- bimodLikData(c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(lrt_diff, 3, lower.tail = F))
}

bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- minmax(length(x2) / length(x), min = 1e-05, max = (1 - 1e-05))
  likA <- length(x1) * log(1 - xal)
  mysd <- sd(x2)
  if (length(x2) < 2) {
    mysd <- 1
  }
  likB <- length(x2) * log(xal) + sum(dnorm(x2, mean(x2), mysd, log = TRUE))
  return(likA + likB)
}

ainb <- function(a, b) {
  a2 <- a[a %in% b]
  return(a2)
}

minmax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}
set.ifnull <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  return(x)
}

expMean <- function(x) {
  return(log(mean(exp(x) - 1) + 1))
}


DiffExpTest <- function(expression, cells.1, cells.2, genes.use = NULL, print.bar = TRUE) {
  cat(file = stderr(), "DiffExpTest\n")
  genes.use <- set.ifnull(genes.use, rownames(expression))
  p_val <- unlist(lapply(genes.use, function(x) diffLRT(as.numeric(expression[x, cells.1]), as.numeric(expression[x, cells.2]))))
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}
