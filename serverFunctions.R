# functions should go in external file

n_fun <- function(x) {
  return(data.frame(y = -.5, label = paste0(length(x), "\ncells")))
}

diffLRT = function(x, y, xmin = 1) {
  lrtX = bimodLikData(x)
  lrtY = bimodLikData(y)
  lrtZ = bimodLikData(c(x, y))
  lrt_diff = 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(lrt_diff, 3, lower.tail = F))
}

bimodLikData = function(x, xmin = 0) {
  x1 = x[x <= xmin]
  x2 = x[x > xmin]
  xal = minmax(length(x2) / length(x),
               min = 1e-5,
               max = (1 - 1e-5))
  likA = length(x1) * log(1 - xal)
  mysd = sd(x2)
  if (length(x2) < 2) {
    mysd = 1
  }
  likB = length(x2) * log(xal) + sum(dnorm(x2, mean(x2), mysd, log = TRUE))
  return(likA + likB)
}

ainb = function(a, b) {
  a2 = a[a %in% b]
  return(a2)
}

minmax = function(data, min, max) {
  data2 = data
  data2[data2 > max] = max
  data2[data2 < min] = min
  return(data2)
}
set.ifnull = function(x, y) {
  if (is.null(x))
    return(y)
  return(x)
}

expMean = function(x) {
  return(log(mean(exp(x) - 1) + 1))
}


DiffExpTest = function(expression,
                       cells.1,
                       cells.2,
                       genes.use = NULL,
                       print.bar = TRUE) {
  cat(file=stderr(), "DiffExpTest\n")
  genes.use = set.ifnull(genes.use, rownames(expression))
  p_val = unlist(lapply(genes.use, function(x)
    diffLRT(
      as.numeric(expression[x, cells.1]), as.numeric(expression[x, cells.2])
    )))
  to.return = data.frame(p_val, row.names = genes.use)
  return(to.return)
}
