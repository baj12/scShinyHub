rm(list = ls())
library(cellrangerRkit)
packageVersion("cellrangerRkit")

# from CellView to create original data
log2cpm <- read.csv("Data/Expression.csv", row.names = 1, stringsAsFactors = F, as.is = T, check.names = F)
featuredata <- read.csv("Databases/HG19_v74_FeatureData.csv", row.names = 1, stringsAsFactors = F, as.is = T, sep = ",")
tsne.data <- read.csv("Data/TNSE_dbscan.csv", row.names = 1, stringsAsFactors = F, as.is = T)



################################# this loads the original CellView example data including featuredata, log2cpm, tsne.data
load("Examples/PBMC-Apheresis.Rds")


pd <- data.frame(row.names = colnames(log2cpm), sampleNames = rep("1", ncol(log2cpm)), barcode = colnames(log2cpm))
mat <- as.matrix(round(2^(log2cpm) - 1))
# here we handle infinite values that come from the ENSGUMI and ENSGGENES
infRows <- unique(which(is.infinite(mat), arr.ind = TRUE)[, 1])
if (length(infRows) > 0) {
  mat <- mat[-infRows, ]
}

mat <- mat[rownames(mat) %in% rownames(featuredata), ]
fd <- featuredata[rownames(mat), ]
gbm <- newGeneBCMatrix(mat = as(mat, "dgTMatrix"), fd = fd, pd = pd)
fData(gbm)$symbol <- fData(gbm)$symbol
fData(gbm)$id <- rownames(fData(gbm))
save(file = "Examples/PBMC-Apheresis.new.Rds", list = c("featuredata", "gbm"))

#################################
rm(list = ls())

# this loads the original CellView example data including featuredata, log2cpm, tsne.data
load("Examples/Normal_Pancreatic_Islets_10X.Rds")

pd <- data.frame(row.names = colnames(log2cpm), sample = rep("1", ncol(log2cpm)))
mat <- as.matrix(round(2^(log2cpm) - 1))
# here we handle infinite values that come from the ENSGUMI and ENSGGENES
infRows <- unique(which(is.infinite(mat), arr.ind = TRUE)[, 1])
if (length(infRows) > 0) {
  mat <- mat[-infRows, ]
}
fd <- featuredata[rownames(mat), ]
gbm <- newGeneBCMatrix(mat = as(mat, "dgTMatrix"), fd = fd, pd = pd)
fData(gbm)$symbol <- fData(gbm)$symbol
fData(gbm)$id <- rownames(fData(gbm))
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm)
gbm_log <- log_gene_bc_matrix(gbm_bcnorm, base = 10)
save(file = "Examples/PBMC-Normal_Pancreatic_Islets_10X.new.Rds", list = c("featuredata", "tsne.data", "log2cpm", "gbm", "gbm_log"))
