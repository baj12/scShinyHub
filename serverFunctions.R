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

  notFound <- g_id[!toupper(g_id) %in% featureData$Associated.Gene.Name]
  if (length(featureData$Associated.Gene.Name) == length(notFound)) {
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

  geneid <- rownames(featureData[which(featureData$Associated.Gene.Name %in% toupper(g_id)), ])
  return(geneid)
}


updateProjectionsWithUmiCount <- function(dimX, dimY, geneNames, geneNames2 = NULL, featureData, gbm, projections) {
  if ((dimY == "UmiCountPerGenes") | (dimX == "UmiCountPerGenes")) {
    geneNames <- geneName2Index(geneNames, featureData)
    if ((length(geneNames) > 0) & (length(geneNames[[1]]) > 0)) {
      if (length(geneNames) == 1) {
        projections$UmiCountPerGenes <- exprs(gbm)[geneNames, ]
      } else {
        projections$UmiCountPerGenes <- Matrix::colSums(exprs(gbm)[geneNames, ])
      }
    } else {
      projections$UmiCountPerGenes <- 0
    }
  }
  if ((dimY == "UmiCountPerGenes2") | (dimX == "UmiCountPerGenes2")) {
    geneNames <- geneName2Index(geneNames2, featureData)
    if ((length(geneNames) > 0) & (length(geneNames[[1]]) > 0)) {
      if (length(geneNames) == 1) {
        projections$UmiCountPerGenes2 <- exprs(gbm)[geneNames, ]
      } else {
        projections$UmiCountPerGenes2 <- Matrix::colSums(exprs(gbm)[geneNames, ])
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
